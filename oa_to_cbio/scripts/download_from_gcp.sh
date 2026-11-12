#!/bin/bash
# ---------------------------------------------------------------------------
# download_from_gcp.sh
#
# Selectively downloads only the files needed for cBioPortal conversion from
# a GCP bucket containing Oncoanalyser WGTS outputs.
#
# Patients are routed into per-cohort subdirectories derived from their donor
# ID prefix (e.g., MOHCCN-P7-00091 -> $OUTDIR/MOHCCN-P7/MOHCCN-P7-00091/).
#
# Usage (manual, per batch):
#   bash download_from_gcp.sh [OPTIONS]
#
# Options:
#   -b|--bucket    GCP bucket URI  (default: gs://hmf-mohccn-p1-reanalysis-output)
#   -o|--outdir    Local root output directory
#                  (default: /cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download)
#   -p|--patients  Space- or comma-separated patient IDs to download.
#                  If not set, auto-discovers all patients in the bucket.
#   -c|--cohort    Limit to a specific cohort prefix (e.g. MOHCCN-P7).
#                  Only patients whose derived cohort matches will be downloaded.
#   -f|--force     Re-download even if already marked complete in the manifest.
#   -n|--dry-run   Print what would be downloaded without downloading.
#   -h|--help      Show this help.
#
# Output layout:
#   $OUTDIR/
#     $COHORT/               <- derived from patient ID prefix
#       $PATIENT_ID/
#         purple/
#         orange/
#         linx/somatic_annotations/
#         isofox/
#
# Idempotency:
#   A manifest file ($OUTDIR/.download_manifest.tsv) records completed patients.
#   Re-running skips patients already listed as DONE unless --force is set.
#
# GCP Automation note:
#   This script can be called from a Cloud Run Job triggered by a Pub/Sub
#   notification on bucket object creation, or by Cloud Scheduler.
#   Pass --patients to restrict to the newly arrived patient IDs.
# ---------------------------------------------------------------------------

set -euo pipefail
umask 0002

# --- Defaults ---------------------------------------------------------------
BUCKET="gs://hmf-mohccn-p1-reanalysis-output"
OUTDIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download"
PATIENT_IDS=()
COHORT_FILTER=""
FORCE="false"
DRY_RUN="false"
VERIFY="false"
CHECK_CHANGED="false"
ALL_FILES="false"

# --- Parse arguments --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bucket)   BUCKET="$2";          shift 2 ;;
        -o|--outdir)   OUTDIR="$2";          shift 2 ;;
        -p|--patients|--patients-file)
            if [[ -f "$2" ]]; then
                PATIENT_IDS=($(tr '\n\r,' '  ' < "$2"))
            else
                IFS=', ' read -r -a PATIENT_IDS <<< "$2"
            fi
            shift 2 ;;
        -c|--cohort)   COHORT_FILTER="$2";   shift 2 ;;
        -f|--force)    FORCE="true";         shift 1 ;;
        -n|--dry-run)  DRY_RUN="true";       shift 1 ;;
        -v|--verify)   VERIFY="true";        shift 1 ;;
        --check-changed) CHECK_CHANGED="true"; shift 1 ;;
        --all-files)   ALL_FILES="true";     shift 1 ;;
        -h|--help)
            sed -n '2,50p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

BUCKET="${BUCKET%/}"  # Strip trailing slash
MANIFEST="${OUTDIR}/.download_manifest.tsv"

# Verify GCS access / authentication
if ! gcs_check=$(gcloud storage ls "${BUCKET}/" 2>&1); then
    echo "ERROR: Cannot access GCS bucket ${BUCKET}."
    echo "Please make sure you are authenticated with Google Cloud by running:"
    echo "  gcloud auth login"
    echo ""
    echo "Error details:"
    echo "$gcs_check"
    exit 1
fi

mkdir -p "$OUTDIR"
chmod g+rwX "$OUTDIR" 2>/dev/null || true
if [[ ! -e "$MANIFEST" ]]; then
    touch "$MANIFEST"
fi
chmod g+rw "$MANIFEST" 2>/dev/null || true

# ---------------------------------------------------------------------------
# Derive cohort prefix from a patient ID by stripping the trailing numeric part.
# Examples:
#   MOHCCN-P7-00091 -> MOHCCN-P7
#   MOHCCN-P9-01008 -> MOHCCN-P9
#   PANCAN-G-0584   -> PANCAN-G
# ---------------------------------------------------------------------------
get_cohort() {
    echo "$1" | sed 's/-[0-9]*$//'
}

# ---------------------------------------------------------------------------
# Files to download per patient.
# Two parallel arrays: DOWNLOAD_SUBDIRS[i] and DOWNLOAD_PATTERNS_LIST[i]
# (Associative arrays require bash 4+; parallel arrays work on bash 3.2+)
# Alignments, teal BAMs, sigs, neo, qsee, logs, and circos/plots are excluded.
# ---------------------------------------------------------------------------
DOWNLOAD_SUBDIRS=(
    "purple"
    "orange"
    "linx/somatic_annotations"
    "isofox"
    "teal"
)
DOWNLOAD_PATTERNS_LIST=(
    "*.purple.somatic.vcf.gz *.purple.somatic.vcf.gz.tbi *.purple.germline.vcf.gz *.purple.germline.vcf.gz.tbi *.purple.cnv.gene.tsv *.purple.segment.tsv *.purple.purity.tsv *.purple.qc"
    "*.orange.json *.orange.pdf"
    "*.linx.vis_sv_data.tsv *.linx.svs.tsv *.linx.breakend.tsv *.linx.fusion.tsv *.linx.drivers.tsv *.linx.clusters.tsv *.linx.links.tsv"
    "*.isf.gene_data.tsv *.isf.fusions.tsv *.isf.summary.tsv"
    "*.teal.tellength.tsv"
)

# ---------------------------------------------------------------------------
# Check if a patient is already marked DONE in the manifest.
# Uses awk for portability (no grep -P needed).
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
        purple_vcf=$(ls "${OUTDIR}/${cohort}/${patient}/purple/"*.purple.somatic.vcf.gz 2>/dev/null | head -n 1 || true)
        cnv=$(ls "${OUTDIR}/${cohort}/${patient}/purple/"*.purple.cnv.gene.tsv 2>/dev/null | head -n 1 || true)
        purity=$(ls "${OUTDIR}/${cohort}/${patient}/purple/"*.purple.purity.tsv 2>/dev/null | head -n 1 || true)
        segment=$(ls "${OUTDIR}/${cohort}/${patient}/purple/"*.purple.segment.tsv 2>/dev/null | head -n 1 || true)
        orange_json=$(ls "${OUTDIR}/${cohort}/${patient}/orange/"*.orange.json 2>/dev/null | head -n 1 || true)
        orange_pdf=$(ls "${OUTDIR}/${cohort}/${patient}/orange/"*.orange.pdf 2>/dev/null | head -n 1 || true)
        linx_drivers=$(ls "${OUTDIR}/${cohort}/${patient}/linx/somatic_annotations/"*.linx.drivers.tsv 2>/dev/null | head -n 1 || true)

        if [[ -z "$purple_vcf" ]] || [[ -z "$cnv" ]] || [[ -z "$purity" ]] || [[ -z "$segment" ]] || \
           [[ -z "$orange_json" ]] || [[ -z "$orange_pdf" ]] || [[ -z "$linx_drivers" ]]; then
            return 1 # Missing required files
        fi

        # If it is not a PANCAN-G patient, check isofox files as well
        if [[ "$cohort" != "PANCAN-G" ]]; then
            local isf_gene
            isf_gene=$(ls "${OUTDIR}/${cohort}/${patient}/isofox/"*.isf.gene_data.tsv 2>/dev/null | head -n 1 || true)
            if [[ -z "$isf_gene" ]]; then
                return 1 # Missing expected RNA files
            fi
        fi

        # If --all-files is set, check key files from other subdirectories as well
        if [[ "$ALL_FILES" == "true" ]]; then
            local cobalt chord amber
            cobalt=$(ls "${OUTDIR}/${cohort}/${patient}/cobalt/"*.cobalt.ratio.pcf 2>/dev/null | head -n 1 || true)
            chord=$(ls "${OUTDIR}/${cohort}/${patient}/chord/"*.chord.prediction.tsv 2>/dev/null | head -n 1 || true)
            amber=$(ls "${OUTDIR}/${cohort}/${patient}/amber/"*.amber.contamination.tsv 2>/dev/null | head -n 1 || true)
            if [[ -z "$cobalt" ]] || [[ -z "$chord" ]] || [[ -z "$amber" ]]; then
                return 1 # Missing expected additional directories
            fi
        fi
    fi

    return 0
}

mark_done() {
    local patient="$1" cohort="$2"
    grep -v "^${patient}	" "$MANIFEST" > "${MANIFEST}.tmp" && mv "${MANIFEST}.tmp" "$MANIFEST" || true
    printf "%s\t%s\tDONE\t%s\n" "$patient" "$cohort" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$MANIFEST"
    chmod g+rw "$MANIFEST" 2>/dev/null || true
}

mark_failed() {
    local patient="$1" cohort="$2"
    grep -v "^${patient}	" "$MANIFEST" > "${MANIFEST}.tmp" && mv "${MANIFEST}.tmp" "$MANIFEST" || true
    printf "%s\t%s\tFAILED\t%s\n" "$patient" "$cohort" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$MANIFEST"
    chmod g+rw "$MANIFEST" 2>/dev/null || true
}

# ---------------------------------------------------------------------------
# Auto-discover patients if none specified
# ---------------------------------------------------------------------------
if [[ ${#PATIENT_IDS[@]} -eq 0 ]]; then
    echo "No --patients specified. Auto-discovering from ${BUCKET}/..."
    # Use sed instead of grep -P for macOS/bash 3.2 compatibility
    while IFS= read -r line; do
        p=$(echo "$line" | sed 's|/$||' | sed 's|.*/||')
        [[ -n "$p" ]] && PATIENT_IDS+=("$p")
    done < <(gcloud storage ls "${BUCKET}/" 2>/dev/null || true)
    echo "Found ${#PATIENT_IDS[@]} patient(s) in bucket."
fi

# Apply cohort filter if specified
if [[ -n "$COHORT_FILTER" ]]; then
    FILTERED=()
    for p in "${PATIENT_IDS[@]}"; do
        if [[ "$(get_cohort "$p")" == "$COHORT_FILTER" ]]; then
            FILTERED+=("$p")
        fi
    done
    echo "Cohort filter '${COHORT_FILTER}': ${#FILTERED[@]} of ${#PATIENT_IDS[@]} patient(s) match."
    PATIENT_IDS=("${FILTERED[@]}")
fi

# ---------------------------------------------------------------------------
# Main download loop
# ---------------------------------------------------------------------------
DOWNLOADED=0
SKIPPED=0
FAILED=0

for patient in "${PATIENT_IDS[@]}"; do
    COHORT=$(get_cohort "$patient")

    echo ""
    echo "============================================================"
    echo "Patient: ${patient}  |  Cohort: ${COHORT}"
    echo "============================================================"

    if [[ "$FORCE" != "true" ]] && is_done "$patient"; then
        echo "  Already marked DONE in manifest. Skipping (use --force to re-download)."
        (( SKIPPED++ )) || true
        continue
    fi

    PATIENT_OUTDIR="${OUTDIR}/${COHORT}/${patient}"

    # Determine GCS prefix and layout (handles standard, nested, and flat bucket layouts dynamically)
    gcs_prefix="${BUCKET}/${patient}"
    is_flat="false"
    if gcloud storage ls "${BUCKET}/${patient}/purple/" &>/dev/null; then
        # Standard layout: ${BUCKET}/${patient}/purple/
        gcs_prefix="${BUCKET}/${patient}"
    elif gcloud storage ls "${BUCKET}/${patient}/${patient}/purple/" &>/dev/null || gcloud storage ls "${BUCKET}/${patient}/${patient}/" &>/dev/null; then
        # Nested layout: ${BUCKET}/${patient}/${patient}/purple/
        gcs_prefix="${BUCKET}/${patient}/${patient}"
    elif gcloud storage ls "${BUCKET}/${patient}/*.purple.*" &>/dev/null; then
        # Flat layout: files stored directly in ${BUCKET}/${patient}/ without subdirs
        gcs_prefix="${BUCKET}/${patient}"
        is_flat="true"
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] Would create: ${PATIENT_OUTDIR}/"
        if [[ "$ALL_FILES" == "true" ]]; then
            echo "  [DRY RUN] Would sync recursively: ${gcs_prefix} -> ${PATIENT_OUTDIR}"
            echo "            (excluding work, alignments, and BAMs/CRAMs/FASTQs)"
        else
            for i in "${!DOWNLOAD_SUBDIRS[@]}"; do
                subdir="${DOWNLOAD_SUBDIRS[$i]}"
                for pattern in ${DOWNLOAD_PATTERNS_LIST[$i]}; do
                    if [[ "$is_flat" == "true" ]]; then
                        echo "  [DRY RUN] Would download: ${gcs_prefix}/${pattern} -> ${PATIENT_OUTDIR}/${subdir}/"
                    else
                        echo "  [DRY RUN] Would download: ${gcs_prefix}/${subdir}/${pattern}"
                    fi
                done
            done
        fi
        continue
    fi

    success="true"
    if [[ "$ALL_FILES" == "true" ]]; then
        echo "  Syncing all files recursively from ${gcs_prefix} ..."
        mkdir -p "${PATIENT_OUTDIR}"

        # Build exclude filter for gcloud storage rsync
        # - Exclude work/ and alignments/ directories
        # - Exclude *.bam, *.bai, *.cram, *.crai, and *.fastq.gz files
        exclude_regex="(^|.*/)(work|alignments)(/.*|$)|\.(bam|bai|cram|crai|fastq\.gz)$"

        # Build clobber/no-clobber options
        extra_opts=()
        if [[ "$CHECK_CHANGED" == "false" ]]; then
            extra_opts+=("--no-clobber")
        fi

        if ! rsync_out=$(gcloud storage rsync "${gcs_prefix}" "${PATIENT_OUTDIR}" \
            --recursive \
            --exclude="${exclude_regex}" \
            "${extra_opts[@]}" 2>&1); then
            echo "  ERROR: Sync failed for ${patient}."
            echo "$rsync_out"
            success="false"
        fi
    else
        for i in "${!DOWNLOAD_SUBDIRS[@]}"; do
            subdir="${DOWNLOAD_SUBDIRS[$i]}"
            LOCAL_SUBDIR="${PATIENT_OUTDIR}/${subdir}"
            mkdir -p "$LOCAL_SUBDIR"

            for pattern in ${DOWNLOAD_PATTERNS_LIST[$i]}; do
                # Skip if file already exists locally, unless --check-changed is set
                if [[ "$CHECK_CHANGED" == "false" ]]; then
                    file_exists="false"
                    for f in "${LOCAL_SUBDIR}"/${pattern}; do
                        if [[ -f "$f" ]]; then
                            file_exists="true"
                            break
                        fi
                    done
                    if [[ "$file_exists" == "true" ]]; then
                        echo "  Skipping ${subdir}/${pattern} (already exists locally)"
                        continue
                    fi
                fi

                echo "  Downloading ${subdir}/${pattern} ..."
                if [[ "$is_flat" == "true" ]]; then
                    src_uri="${gcs_prefix}/${pattern}"
                else
                    src_uri="${gcs_prefix}/${subdir}/${pattern}"
                fi

                # || true: a missing optional file (e.g. isofox absent for DNA-only) won't abort
                gcloud storage cp \
                    "${src_uri}" \
                    "${LOCAL_SUBDIR}/" \
                    --no-clobber 2>&1 | grep -v "^$" || true
            done
        done
    fi

    # Verify minimum required files (somatic VCF + orange JSON are mandatory)
    SOMATIC_VCF=$(ls "${PATIENT_OUTDIR}/purple/"*.purple.somatic.vcf.gz 2>/dev/null | head -n 1 || true)
    ORANGE_JSON=$(ls "${PATIENT_OUTDIR}/orange/"*.orange.json 2>/dev/null | head -n 1 || true)

    if [[ "$success" == "false" ]] || [[ -z "$SOMATIC_VCF" ]] || [[ -z "$ORANGE_JSON" ]]; then
        echo "  ERROR: Sync failed or missing mandatory files (somatic VCF or orange JSON) for ${patient}."
        mark_failed "$patient" "$COHORT"
        (( FAILED++ )) || true
        success="false"
    fi

    if [[ "$success" == "true" ]]; then
        mark_done "$patient" "$COHORT"
        (( DOWNLOADED++ )) || true
        echo "  Done: ${patient} -> ${COHORT}/"
    fi
done

echo ""
echo "============================================================"
echo "Download summary:"
echo "  Downloaded: ${DOWNLOADED}"
echo "  Skipped (already done): ${SKIPPED}"
echo "  Failed: ${FAILED}"
echo "  Manifest: ${MANIFEST}"
echo "============================================================"

if [[ $FAILED -gt 0 ]]; then
    echo "WARNING: ${FAILED} patient(s) failed. Check logs above."
    exit 1
fi
