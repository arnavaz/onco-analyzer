#!/bin/bash
# ---------------------------------------------------------------------------
# sync_from_gcp.sh
#
# Detects new patients in the GCP bucket, groups them by cohort (derived from
# the donor ID prefix), downloads each cohort's patients selectively, then
# submits a separate oa_to_cbio conversion pipeline per cohort.
#
# Cohort prefix extraction examples:
#   MOHCCN-P7-00091 -> MOHCCN-P7
#   MOHCCN-P9-00567 -> MOHCCN-P9
#   PANCAN-G-0584   -> PANCAN-G
#
# Per-cohort output layout:
#   $DOWNLOAD_DIR/$COHORT/$PATIENT_ID/{purple,orange,linx,...}
#   $CBIO_OUT_DIR/$COHORT/import_folder/
#
# Pre-generated clinical files are looked up at:
#   $PRE_CLIN_DIR/$COHORT/patient_clinical.tsv  (patient-level)
#   $PRE_CLIN_DIR/$COHORT/sample_clinical.tsv   (sample-level)
# If not found, the conversion runs without pre-population.
#
# Usage:
#   bash sync_from_gcp.sh [OPTIONS]
#
#   --cohort <ID>     Process only this cohort prefix (e.g. MOHCCN-P7).
#                     Defaults to all cohorts found in the bucket.
#   -A|--account <ID> SLURM account to charge jobs to (default: pughlab).
#   -p|--patients|--patients-file <LIST_OR_FILE>
#                     Process only these patient IDs (space/comma string or file path).
#   --download-only   Download files only; skip SLURM job submission (Step 3).
#                     Use this on your Mac (with sshfs) to write files directly
#                     to the cluster filesystem, then run with --skip-download
#                     on H4H to submit the conversion jobs.
#   --merge-only      Skip per-sample jobs; run only the merge step per cohort.
#   --clinical-only   Process cohorts using only pre-generated clinical data from
#                     $PRE_CLIN_DIR (skips GCP download and per-sample genomic jobs;
#                     does not require $DOWNLOAD_DIR/$COHORT to exist).
#   --dry-run         Print what would happen without downloading or submitting.
#   --skip-download   Already have files locally; skip download, run conversion only.
#                     Use this on H4H after a --download-only run on your Mac.
#   -h|--help         Show this help.
#
# Two-machine workflow (Mac has gcloud + sshfs; H4H has SLURM):
#   Mac:  bash sync_from_gcp.sh --cohort PANCAN-G --download-only
#         (writes via sshfs directly to cluster filesystem)
#   H4H:  bash sync_from_gcp.sh --cohort PANCAN-G --skip-download
#         (submits oa_to_cbio SLURM jobs on already-downloaded data)
#
# Compatibility: bash 3.2+ (macOS default), bash 4+ (Linux/cluster)
#
# GCP Automation note:
#   This script can be invoked by a GCP Cloud Run Job triggered by Cloud
#   Scheduler (e.g., nightly) or a Pub/Sub notification on bucket object
#   finalization. When automating, upload pre-clinical files to a GCP bucket
#   and gsutil-copy them to $PRE_CLIN_DIR before calling this script.
# ---------------------------------------------------------------------------

set -euo pipefail
umask 0002

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OA_TO_CBIO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --- Configuration ----------------------------------------------------------
BUCKET="gs://hmf-mohccn-p1-reanalysis-output"
DOWNLOAD_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download"
CBIO_OUT_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/cbio_import"
PRE_CLIN_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/pre_gen_clin_data"
SLURM_ACCOUNT="pughlab"
MANIFEST="${DOWNLOAD_DIR}/.download_manifest.tsv"

COHORT_FILTER=""
PATIENT_FILTER_STR=""
MERGE_ONLY="false"
CLINICAL_ONLY="false"
DRY_RUN="false"
SKIP_DOWNLOAD="false"
DOWNLOAD_ONLY="false"
VERIFY="false"
CHECK_CHANGED="false"
ALL_FILES="false"
FORCE="false"

# --- Parse arguments --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort)         COHORT_FILTER="$2"; shift 2 ;;
        -A|--account)     SLURM_ACCOUNT="$2"; shift 2 ;;
        -p|--patients|--patients-file)
            if [[ -f "$2" ]]; then
                PATIENT_FILTER_STR=$(tr '\n\r,' '  ' < "$2" | xargs)
            else
                PATIENT_FILTER_STR=$(echo "$2" | tr ',' ' ' | xargs)
            fi
            shift 2 ;;
        --merge-only)     MERGE_ONLY="true";  shift ;;
        --clinical-only)  CLINICAL_ONLY="true"; SKIP_DOWNLOAD="true"; MERGE_ONLY="true"; shift ;;
        --dry-run)        DRY_RUN="true";     shift ;;
        --skip-download)  SKIP_DOWNLOAD="true"; shift ;;
        --download-only)  DOWNLOAD_ONLY="true"; shift ;;
        --verify)         VERIFY="true";      shift ;;
        --check-changed)  CHECK_CHANGED="true"; shift ;;
        --all-files)      ALL_FILES="true";   shift ;;
        --force)          FORCE="true";       shift ;;
        -h|--help)
            sed -n '2,55p' "$0" | grep '^#' | sed 's/^# //; s/^#$//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Derive cohort prefix from a patient ID.
# Strips the trailing -NNNNN numeric component.
# Works with sed (POSIX, no -P flag needed).
# ---------------------------------------------------------------------------
get_cohort() {
    echo "$1" | sed 's/-[0-9]*$//'
}

# ---------------------------------------------------------------------------
# Derive cBioPortal study metadata from a cohort prefix.
# Echoes three tab-separated fields: STUDY_ID  STUDY_NAME  TYPE_OF_CANCER
# Add new cohort mappings here as they arrive.
# ---------------------------------------------------------------------------
get_study_meta() {
    local cohort="$1"
    local study_id study_name type_of_cancer
    # Normalize to lowercase with underscores for study_id
    study_id=$(echo "$cohort" | tr '[:upper:]' '[:lower:]' | tr '-' '_')
    study_id="${study_id}_wgts"
    type_of_cancer="mixed"

    case "$cohort" in
        MOHCCN-P7) study_name="MOHCCN-P7 WGTS" ;;
        MOHCCN-P9) study_name="MOHCCN-P9 WGTS" ;;
        PANCAN-G)  study_name="PANCAN-G WGTS" ;;
        *)         study_name="${cohort} WGTS" ;;  # Generic fallback
    esac

    printf "%s\t%s\t%s\n" "$study_id" "$study_name" "$type_of_cancer"
}

# ---------------------------------------------------------------------------
# Check if a patient is already marked DONE in the manifest.
# Uses awk for portability (no grep -P).
# Manifest format: patient<TAB>cohort<TAB>STATUS<TAB>timestamp
# ---------------------------------------------------------------------------
is_done() {
    local patient="$1"
    local cohort
    cohort=$(get_cohort "$patient")

    # Check manifest status
    if ! awk -v p="$patient" '$1 == p && $3 == "DONE" { found=1 } END { exit !found }' "$MANIFEST" 2>/dev/null; then
        return 1
    fi

    # If --verify is set, check if all expected files are present locally
    if [[ "$VERIFY" == "true" ]]; then
        # Check purple files
        local purple_vcf cnv purity segment orange_json orange_pdf linx_drivers
        purple_vcf=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/purple/"*.purple.somatic.vcf.gz 2>/dev/null | head -n 1 || true)
        cnv=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/purple/"*.purple.cnv.gene.tsv 2>/dev/null | head -n 1 || true)
        purity=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/purple/"*.purple.purity.tsv 2>/dev/null | head -n 1 || true)
        segment=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/purple/"*.purple.segment.tsv 2>/dev/null | head -n 1 || true)
        orange_json=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/orange/"*.orange.json 2>/dev/null | head -n 1 || true)
        orange_pdf=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/orange/"*.orange.pdf 2>/dev/null | head -n 1 || true)
        linx_drivers=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/linx/somatic_annotations/"*.linx.drivers.tsv 2>/dev/null | head -n 1 || true)

        if [[ -z "$purple_vcf" ]] || [[ -z "$cnv" ]] || [[ -z "$purity" ]] || [[ -z "$segment" ]] || \
           [[ -z "$orange_json" ]] || [[ -z "$orange_pdf" ]] || [[ -z "$linx_drivers" ]]; then
            return 1 # Missing required files
        fi

        # If it is not a PANCAN-G patient, check isofox files as well
        if [[ "$cohort" != "PANCAN-G" ]]; then
            local isf_gene
            isf_gene=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/isofox/"*.isf.gene_data.tsv 2>/dev/null | head -n 1 || true)
            if [[ -z "$isf_gene" ]]; then
                return 1 # Missing expected RNA files
            fi
        fi

        # If --all-files is set, check key files from other subdirectories as well
        if [[ "$ALL_FILES" == "true" ]]; then
            local cobalt chord amber
            cobalt=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/cobalt/"*.cobalt.ratio.pcf 2>/dev/null | head -n 1 || true)
            chord=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/chord/"*.chord.prediction.tsv 2>/dev/null | head -n 1 || true)
            amber=$(ls "${DOWNLOAD_DIR}/${cohort}/${patient}/amber/"*.amber.contamination.tsv 2>/dev/null | head -n 1 || true)
            if [[ -z "$cobalt" ]] || [[ -z "$chord" ]] || [[ -z "$amber" ]]; then
                return 1 # Missing expected additional directories
            fi
        fi
    fi

    return 0
}

# ---------------------------------------------------------------------------
# Convert bytes to human-readable string (no numfmt dependency).
# ---------------------------------------------------------------------------
bytes_to_human() {
    local bytes=$1
    if [[ $bytes -ge 1073741824 ]]; then
        awk "BEGIN {printf \"%.1f GB\", ${bytes}/1073741824}"
    elif [[ $bytes -ge 1048576 ]]; then
        awk "BEGIN {printf \"%.1f MB\", ${bytes}/1048576}"
    elif [[ $bytes -ge 1024 ]]; then
        awk "BEGIN {printf \"%.1f KB\", ${bytes}/1024}"
    else
        echo "${bytes} B"
    fi
}

if [[ "$SKIP_DOWNLOAD" == "false" ]] && [[ "$CLINICAL_ONLY" == "false" ]]; then
    mkdir -p "$DOWNLOAD_DIR"
    chmod g+rwX "$DOWNLOAD_DIR" 2>/dev/null || true
    if [[ ! -e "$MANIFEST" ]]; then
        touch "$MANIFEST"
    fi
    chmod g+rw "$MANIFEST" 2>/dev/null || true
fi

echo "========================================================"
echo "GCP → cBioPortal multi-cohort sync"
echo "  Bucket:          ${BUCKET}"
echo "  Download root:   ${DOWNLOAD_DIR}"
echo "  cBio output:     ${CBIO_OUT_DIR}"
echo "  Pre-clin dir:    ${PRE_CLIN_DIR}"
echo "  SLURM Account:   ${SLURM_ACCOUNT}"
echo "  Cohort filter:   ${COHORT_FILTER:-<all>}"
echo "  Skip download:   ${SKIP_DOWNLOAD}"
echo "  Download only:   ${DOWNLOAD_ONLY}"
echo "  Merge only:      ${MERGE_ONLY}"
echo "  Clinical only:   ${CLINICAL_ONLY}"
echo "  Dry run:         ${DRY_RUN}"
echo "========================================================"

# ---------------------------------------------------------------------------
# Step 1: Discover patients and group by cohort.
#   --clinical-only: discovers cohorts from PRE_CLIN_DIR.
#   --skip-download (H4H, no gcloud): reads from local DOWNLOAD_DIR filesystem
#     (and falls back to PRE_CLIN_DIR for clinical-only cohorts).
#   default (Mac, has gcloud): queries the GCP bucket.
# ---------------------------------------------------------------------------
echo ""

# Use a temp directory to group patients by cohort (bash 3.2 compatible)
TMPDIR_COHORTS=$(mktemp -d)
trap 'rm -rf "${TMPDIR_COHORTS}"' EXIT

TMPDIR_NEW="${TMPDIR_COHORTS}/new"
TMPDIR_ALL="${TMPDIR_COHORTS}/all"
mkdir -p "$TMPDIR_NEW" "$TMPDIR_ALL"

if [[ "$CLINICAL_ONLY" == "true" ]]; then
    echo "Step 1: Discovering cohorts from pre-generated clinical directory (--clinical-only mode)..."
    echo "  Reading: ${PRE_CLIN_DIR}/"

    if [[ ! -d "$PRE_CLIN_DIR" ]]; then
        echo "  ERROR: Pre-clinical directory does not exist: ${PRE_CLIN_DIR}"
        exit 1
    fi

    for cohort_dir in "${PRE_CLIN_DIR}"/*/; do
        [[ -d "$cohort_dir" ]] || continue
        cohort=$(basename "$cohort_dir")
        case "$cohort" in
            .*|work|logs|bin|tmp) continue ;;
        esac
        if [[ -n "$COHORT_FILTER" ]] && [[ "$cohort" != "$COHORT_FILTER" ]]; then
            continue
        fi
        touch "${TMPDIR_ALL}/${cohort}"
    done

    if [[ -n "$COHORT_FILTER" ]] && [[ ! -f "${TMPDIR_ALL}/${COHORT_FILTER}" ]]; then
        echo "  ERROR: Cohort '${COHORT_FILTER}' not found in ${PRE_CLIN_DIR}/"
        exit 1
    fi

    CLIN_COHORT_COUNT=$(ls -1 "${TMPDIR_ALL}" 2>/dev/null | wc -l | tr -d ' ')
    echo "  Total clinical-only cohorts found: ${CLIN_COHORT_COUNT}"

elif [[ "$SKIP_DOWNLOAD" == "true" ]]; then
    echo "Step 1: Discovering patients from local filesystem (--skip-download mode)..."
    echo "  Reading: ${DOWNLOAD_DIR}/"

    if [[ ! -d "$DOWNLOAD_DIR" ]] && [[ ! -d "$PRE_CLIN_DIR" ]]; then
        echo "  ERROR: Neither download directory (${DOWNLOAD_DIR}) nor pre-clinical directory (${PRE_CLIN_DIR}) exists."
        echo "  Run on Mac first: sync_from_gcp.sh --download-only"
        exit 1
    fi

    # Walk local cohort subdirs: DOWNLOAD_DIR/COHORT/PATIENT/
    if [[ -d "$DOWNLOAD_DIR" ]]; then
        for cohort_dir in "${DOWNLOAD_DIR}"/*/; do
            [[ -d "$cohort_dir" ]] || continue
            cohort=$(basename "$cohort_dir")
            if [[ -n "$COHORT_FILTER" ]] && [[ "$cohort" != "$COHORT_FILTER" ]]; then
                continue
            fi
            for patient_dir in "${cohort_dir}"*/; do
                [[ -d "$patient_dir" ]] || continue
                patient=$(basename "$patient_dir")
                if [[ -n "$PATIENT_FILTER_STR" ]]; then
                    is_matched=false
                    for pid in $PATIENT_FILTER_STR; do
                        if [[ "$patient" == "$pid" ]]; then
                            is_matched=true
                            break
                        fi
                    done
                    if [[ "$is_matched" == "false" ]]; then
                        continue
                    fi
                fi
                echo "$patient" >> "${TMPDIR_ALL}/${cohort}"
            done
        done
    fi

    TOTAL=$(cat "${TMPDIR_ALL}"/* 2>/dev/null | wc -l | tr -d ' ')
    echo "  Total patients found locally: ${TOTAL}"

else
    echo "Step 1: Discovering patients in bucket..."

    # Extract patient IDs: strip trailing slash first, then strip GCS path prefix
    if ! ALL_PATIENTS=$(gcloud storage ls "${BUCKET}/" 2>&1); then
        echo "ERROR: Failed to list bucket contents from GCS."
        echo "This is usually because you are not logged in or your session has expired."
        echo "Please run:"
        echo "  gcloud auth login"
        echo ""
        echo "Error details:"
        echo "$ALL_PATIENTS"
        exit 1
    fi
    ALL_PATIENTS=$(echo "$ALL_PATIENTS" | sed 's|/$||' | sed 's|.*/||' | grep -v '^$' || true)

    TOTAL=$(echo "$ALL_PATIENTS" | grep -v '^$' | wc -l | tr -d ' ')
    echo "  Total patients in bucket: ${TOTAL}"

    for patient in $ALL_PATIENTS; do
        cohort=$(get_cohort "$patient")
        if [[ -n "$COHORT_FILTER" ]] && [[ "$cohort" != "$COHORT_FILTER" ]]; then
            continue
        fi
        if [[ -n "$PATIENT_FILTER_STR" ]]; then
            is_matched=false
            for pid in $PATIENT_FILTER_STR; do
                if [[ "$patient" == "$pid" ]]; then
                    is_matched=true
                    break
                fi
            done
            if [[ "$is_matched" == "false" ]]; then
                continue
            fi
        fi
        echo "$patient" >> "${TMPDIR_ALL}/${cohort}"
        if [[ "$FORCE" == "true" ]] || ! is_done "$patient"; then
            echo "$patient" >> "${TMPDIR_NEW}/${cohort}"
        fi
    done
fi

# Auto-detect clinical-only cohorts in PRE_CLIN_DIR when their GCP/download folder
# does not exist or is empty (0 genomic patients found).
CLIN_ONLY_COUNT=0
if [[ "$CLINICAL_ONLY" != "true" ]] && [[ -d "$PRE_CLIN_DIR" ]] && [[ -z "$PATIENT_FILTER_STR" ]]; then
    for cohort_dir in "${PRE_CLIN_DIR}"/*/; do
        [[ -d "$cohort_dir" ]] || continue
        cohort=$(basename "$cohort_dir")
        case "$cohort" in
            .*|work|logs|bin|tmp) continue ;;
        esac
        if [[ -n "$COHORT_FILTER" ]] && [[ "$cohort" != "$COHORT_FILTER" ]]; then
            continue
        fi
        existing_count=0
        [[ -f "${TMPDIR_ALL}/${cohort}" ]] && existing_count=$(wc -l < "${TMPDIR_ALL}/${cohort}" | tr -d ' ')
        if [[ $existing_count -eq 0 ]]; then
            has_clin=false
            for fname in data_clinical_patient.txt data_clinical_patients.txt patient_clinical.tsv \
                         data_clinical_sample.txt data_clinical_samples.txt sample_clinical.tsv; do
                if [[ -f "${cohort_dir}/${fname}" ]]; then
                    has_clin=true
                    break
                fi
            done
            if [[ "$has_clin" == "true" ]]; then
                echo "  [${cohort}] No GCP/download patients found -> enabling --clinical-only automatically."
                touch "${TMPDIR_ALL}/${cohort}"
                CLIN_ONLY_COUNT=$(( CLIN_ONLY_COUNT + 1 ))
            fi
        fi
    done
fi

echo ""
echo "  Cohorts found:"
for cohort_file in "${TMPDIR_ALL}"/*; do
    [[ -f "$cohort_file" ]] || continue
    cohort=$(basename "$cohort_file")
    total_c=$(wc -l < "$cohort_file" | tr -d ' ')
    new_c=0
    [[ -f "${TMPDIR_NEW}/${cohort}" ]] && new_c=$(wc -l < "${TMPDIR_NEW}/${cohort}" | tr -d ' ')
    if [[ "$CLINICAL_ONLY" == "true" ]] || [[ $total_c -eq 0 ]]; then
        echo "    ${cohort}: clinical-only (pre-generated clinical data)"
    else
        echo "    ${cohort}: ${total_c} total, ${new_c} new"
    fi
done

# Count total new across all cohorts
ALL_NEW_COUNT=0
for f in "${TMPDIR_NEW}"/*; do
    [[ -f "$f" ]] || continue
    n=$(wc -l < "$f" | tr -d ' ')
    ALL_NEW_COUNT=$(( ALL_NEW_COUNT + n ))
done

# ---------------------------------------------------------------------------
# Dry-run: estimate total download size by sampling one patient per cohort.
# Queries GCP for the subdirs we actually download; multiplies by patient count.
# Only runs when --dry-run is set and there are new patients to download.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" == "true" ]] && [[ "$SKIP_DOWNLOAD" == "false" ]] && [[ $ALL_NEW_COUNT -gt 0 ]]; then
    # These subdirs must match DOWNLOAD_SUBDIRS in download_from_gcp.sh
    DL_SUBDIRS="purple orange linx/somatic_annotations isofox"

    echo ""
    echo "  Estimating download size (sampling 1 patient per cohort from GCP)..."
    TOTAL_ESTIMATE=0

    for cohort_file in "${TMPDIR_NEW}"/*; do
        [[ -f "$cohort_file" ]] || continue
        cohort=$(basename "$cohort_file")
        n_new=$(wc -l < "$cohort_file" | tr -d ' ')
        [[ $n_new -eq 0 ]] && continue

        sample=$(head -1 "$cohort_file")
        patient_bytes=0

        for subdir in $DL_SUBDIRS; do
            # gcloud storage ls -l lists files with sizes; awk sums the byte column
            size=$( (gcloud storage ls -l \
                "${BUCKET}/${sample}/${subdir}/" 2>/dev/null || true) \
                | awk '/^[[:space:]]*[0-9]/ {sum += $1} END {print sum+0}')
            patient_bytes=$(( patient_bytes + size ))
        done

        cohort_bytes=$(( patient_bytes * n_new ))
        TOTAL_ESTIMATE=$(( TOTAL_ESTIMATE + cohort_bytes ))

        echo "    ${cohort}: ~$(bytes_to_human $patient_bytes)/patient × ${n_new} patients = ~$(bytes_to_human $cohort_bytes)"
    done

    echo ""
    echo "  ┌─────────────────────────────────────────────────────────"
    echo "  │  Estimated total download: ~$(bytes_to_human $TOTAL_ESTIMATE)"
    echo "  │  (based on sample patient per cohort; actual may vary)"
    echo "  └─────────────────────────────────────────────────────────"
fi

if [[ $ALL_NEW_COUNT -eq 0 ]] && [[ "$SKIP_DOWNLOAD" == "false" ]]; then
    if [[ $CLIN_ONLY_COUNT -gt 0 ]] && [[ "$DOWNLOAD_ONLY" == "false" ]]; then
        echo ""
        echo "  No new genomic patients to download; proceeding with ${CLIN_ONLY_COUNT} clinical-only cohort(s)."
    else
        echo ""
        echo "  No new patients found in any cohort. Nothing to download."
        echo "  To force a re-run, delete entries from: ${MANIFEST}"
        echo "  To run conversion on already-downloaded data, use --skip-download."
        exit 0
    fi
fi

# ---------------------------------------------------------------------------
# Step 2: Download new patients per cohort
# ---------------------------------------------------------------------------
if [[ "$SKIP_DOWNLOAD" == "false" ]] && [[ "$MERGE_ONLY" != "true" ]] && [[ $ALL_NEW_COUNT -gt 0 ]]; then
    echo ""
    echo "Step 2: Downloading new patients..."

    for cohort_file in "${TMPDIR_NEW}"/*; do
        [[ -f "$cohort_file" ]] || continue
        cohort=$(basename "$cohort_file")

        echo ""
        echo "  --- Cohort: ${cohort} ---"

        # Build comma-separated patient list for download_from_gcp.sh
        NEW_CSV=$(tr '\n' ',' < "$cohort_file" | sed 's/,$//')

        if [[ "$DRY_RUN" == "true" ]]; then
            echo "  [DRY RUN] Would download patients: ${NEW_CSV}"
        else
            # Forward force, verify, check-changed, and all-files flags if set
            force_arg=""
            [[ "$FORCE" == "true" ]] && force_arg="--force"
            verify_arg=""
            [[ "$VERIFY" == "true" ]] && verify_arg="--verify"
            check_changed_arg=""
            [[ "$CHECK_CHANGED" == "true" ]] && check_changed_arg="--check-changed"
            all_files_arg=""
            [[ "$ALL_FILES" == "true" ]] && all_files_arg="--all-files"

            bash "${SCRIPT_DIR}/download_from_gcp.sh" \
                --bucket   "$BUCKET" \
                --outdir   "$DOWNLOAD_DIR" \
                --cohort   "$cohort" \
                --patients "$NEW_CSV" \
                ${force_arg} \
                ${verify_arg} \
                ${check_changed_arg} \
                ${all_files_arg}
        fi
    done
else
    echo ""
    echo "Step 2: Skipped (--skip-download or --merge-only set)."
fi

# ---------------------------------------------------------------------------
# Step 3: Submit oa_to_cbio conversion per cohort
# ---------------------------------------------------------------------------
if [[ "$DOWNLOAD_ONLY" == "true" ]]; then
    echo ""
    echo "Step 3: Skipped (--download-only set)."
    echo ""
    echo "========================================================"
    echo "Download complete. Files written to:"
    for cohort_file in "${TMPDIR_ALL}"/*; do
        [[ -f "$cohort_file" ]] || continue
        cohort=$(basename "$cohort_file")
        echo "  ${DOWNLOAD_DIR}/${cohort}/"
    done
    echo ""
    echo "Now on H4H, run:"
    echo "  bash oa_to_cbio/scripts/sync_from_gcp.sh${COHORT_FILTER:+ --cohort ${COHORT_FILTER}} --skip-download"
    echo "========================================================"
    exit 0
fi

echo ""
echo "Step 3: Submitting oa_to_cbio conversion per cohort..."

for cohort_file in "${TMPDIR_ALL}"/*; do
    [[ -f "$cohort_file" ]] || continue
    cohort=$(basename "$cohort_file")

    COHORT_DOWNLOAD_DIR="${DOWNLOAD_DIR}/${cohort}"
    COHORT_CBIO_OUT="${CBIO_OUT_DIR}/${cohort}"
    COHORT_PATIENT_COUNT=$(wc -l < "$cohort_file" | tr -d ' ')

    # Pre-clinical files: look for cBioPortal-native names first, fall back to legacy tsv names.
    # Prefer cohort-specific dir, fall back to shared root.
    PRE_PATIENT=""
    PRE_SAMPLE=""
    for base in "${PRE_CLIN_DIR}/${cohort}" "${PRE_CLIN_DIR}"; do
        [[ -d "$base" ]] || continue
        # Patient clinical
        if [[ -z "$PRE_PATIENT" ]]; then
            for fname in data_clinical_patient.txt data_clinical_patients.txt patient_clinical.tsv; do
                [[ -f "${base}/${fname}" ]] && PRE_PATIENT="${base}/${fname}" && break
            done
        fi
        # Sample clinical
        if [[ -z "$PRE_SAMPLE" ]]; then
            for fname in data_clinical_sample.txt data_clinical_samples.txt sample_clinical.tsv; do
                [[ -f "${base}/${fname}" ]] && PRE_SAMPLE="${base}/${fname}" && break
            done
        fi
        # Stop as soon as both are found in the same directory
        [[ -n "$PRE_PATIENT" ]] && [[ -n "$PRE_SAMPLE" ]] && break
    done

    # Locate the pre-clinical dir containing the matched files (for timeline copy)
    PRE_CLIN_COHORT_DIR=""
    if [[ -d "${PRE_CLIN_DIR}/${cohort}" ]]; then
        PRE_CLIN_COHORT_DIR="${PRE_CLIN_DIR}/${cohort}"
    elif [[ -d "${PRE_CLIN_DIR}" ]]; then
        PRE_CLIN_COHORT_DIR="${PRE_CLIN_DIR}"
    fi

    # Determine whether this cohort has only clinical data (or GCP download folder is missing/empty)
    COHORT_CLINICAL_ONLY="false"
    LOCAL_PATIENT_DIRS=0
    if [[ -d "$COHORT_DOWNLOAD_DIR" ]]; then
        LOCAL_PATIENT_DIRS=$(find "$COHORT_DOWNLOAD_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [[ "$CLINICAL_ONLY" == "true" ]] || [[ ! -d "$COHORT_DOWNLOAD_DIR" ]] || [[ "$LOCAL_PATIENT_DIRS" -eq 0 ]] || [[ "$COHORT_PATIENT_COUNT" -eq 0 ]]; then
        if [[ "$CLINICAL_ONLY" != "true" ]]; then
            echo "  [${cohort}] GCP download folder (${COHORT_DOWNLOAD_DIR}) does not exist or is empty -> enabling --clinical-only automatically."
        fi
        COHORT_CLINICAL_ONLY="true"
    fi

    # Skip if neither genomic download data nor pre-clinical data exists
    if [[ "$COHORT_CLINICAL_ONLY" == "true" ]] && [[ -z "$PRE_PATIENT" ]] && [[ -z "$PRE_SAMPLE" ]]; then
        echo "  [${cohort}] No local genomic data at ${COHORT_DOWNLOAD_DIR} and no pre-clinical data found. Skipping conversion."
        continue
    fi

    # Derive study metadata
    META=$(get_study_meta "$cohort")
    STUDY_ID=$(echo "$META" | cut -f1)
    STUDY_NAME=$(echo "$META" | cut -f2)
    TYPE_OF_CANCER=$(echo "$META" | cut -f3)
    STUDY_DESC="${STUDY_NAME} data processed with Oncoanalyser v3 and received from GCP bucket."

    echo ""
    echo "  --- Cohort: ${cohort} ---"
    echo "  Study ID:    ${STUDY_ID}"
    echo "  Study Name:  ${STUDY_NAME}"
    echo "  Account:     ${SLURM_ACCOUNT}"
    if [[ "$COHORT_CLINICAL_ONLY" == "true" ]]; then
        echo "  Input:       <clinical-only> (${PRE_CLIN_COHORT_DIR:-none})"
    else
        echo "  Input:       ${COHORT_DOWNLOAD_DIR}"
    fi
    echo "  Output:      ${COHORT_CBIO_OUT}"

    mkdir -p "$COHORT_CBIO_OUT"

    OA_TO_CBIO_CMD=(
        bash "${OA_TO_CBIO_DIR}/oa_to_cbio.sh"
        -o "$COHORT_CBIO_OUT"
        -A "$SLURM_ACCOUNT"
        -s "$STUDY_ID"
        -n "$STUDY_NAME"
        -d "$STUDY_DESC"
        -t "$TYPE_OF_CANCER"
    )

    if [[ "$COHORT_CLINICAL_ONLY" == "false" ]] && [[ -d "$COHORT_DOWNLOAD_DIR" ]]; then
        OA_TO_CBIO_CMD+=(-i "$COHORT_DOWNLOAD_DIR")
    fi

    if [[ -n "$PRE_PATIENT" ]] && [[ -f "$PRE_PATIENT" ]]; then
        OA_TO_CBIO_CMD+=(-P "$PRE_PATIENT")
        echo "  Pre-patient clinical: ${PRE_PATIENT}"
    fi
    if [[ -n "$PRE_SAMPLE" ]] && [[ -f "$PRE_SAMPLE" ]]; then
        OA_TO_CBIO_CMD+=(-S "$PRE_SAMPLE")
        echo "  Pre-sample clinical:  ${PRE_SAMPLE}"
    fi
    if [[ -n "$PRE_CLIN_COHORT_DIR" ]]; then
        OA_TO_CBIO_CMD+=(-C "$PRE_CLIN_COHORT_DIR")
        echo "  Pre-clin dir:         ${PRE_CLIN_COHORT_DIR}"
    fi

    if [[ "$MERGE_ONLY" == "true" ]] || [[ "$COHORT_CLINICAL_ONLY" == "true" ]]; then
        OA_TO_CBIO_CMD+=(--merge-only)
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] Would run: ${OA_TO_CBIO_CMD[*]}"
    else
        "${OA_TO_CBIO_CMD[@]}"
    fi
done

echo ""
echo "========================================================"
echo "Sync complete."
echo ""
echo "Per-cohort cBioPortal output dirs:"
for cohort_file in "${TMPDIR_ALL}"/*; do
    [[ -f "$cohort_file" ]] || continue
    cohort=$(basename "$cohort_file")
    echo "  ${cohort}: ${CBIO_OUT_DIR}/${cohort}/logs/"
done
echo "========================================================"
