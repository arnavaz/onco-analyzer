#!/bin/bash
# ---------------------------------------------------------------------------
# resummarize_clinical.sh
#
# Fast, lightweight clinical re-summarization tool.
# Re-runs process_clinical.py and merge_clinical.py for one or all cohorts
# WITHOUT submitting heavy SLURM jobs or re-running vcf2maf / VEP / CNV / SV.
#
# Usage:
#   bash resummarize_clinical.sh --cohort MOHCCN-P7
#   bash resummarize_clinical.sh --cohort MOHCCN-P7 --clinical-only
#   bash resummarize_clinical.sh --all
#   bash resummarize_clinical.sh --cohort PANCAN-G --dry-run
# ---------------------------------------------------------------------------

set -euo pipefail
umask 0002

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOAD_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download"
CBIO_OUT_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/cbio_import"
PRE_CLIN_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/pre_gen_clin_data"

COHORT_FILTER=""
PROCESS_ALL=false
CLINICAL_ONLY=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort)         COHORT_FILTER="$2"; shift 2 ;;
        --all)            PROCESS_ALL=true;   shift ;;
        --clinical-only)  CLINICAL_ONLY=true; shift ;;
        --download-dir)   DOWNLOAD_DIR="$2";  shift 2 ;;
        --out-dir)        CBIO_OUT_DIR="$2";  shift 2 ;;
        --pre-clin-dir)   PRE_CLIN_DIR="$2";  shift 2 ;;
        --dry-run)        DRY_RUN=true;       shift ;;
        -h|--help)
            echo "Usage: $0 [--cohort <COHORT> | --all] [--clinical-only] [--download-dir <DIR>] [--out-dir <DIR>] [--pre-clin-dir <DIR>] [--dry-run]"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$COHORT_FILTER" && "$PROCESS_ALL" == "false" ]]; then
    echo "ERROR: Please specify --cohort <COHORT> or --all"
    exit 1
fi

add_cohort_unique() {
    local candidate="$1"
    local existing
    for existing in "${COHORTS[@]:-}"; do
        [[ "$existing" == "$candidate" ]] && return 0
    done
    COHORTS+=("$candidate")
}

# Find cohorts to process
COHORTS=()
if [[ -n "$COHORT_FILTER" ]]; then
    COHORTS+=("$COHORT_FILTER")
else
    if [[ "$CLINICAL_ONLY" == "false" && -d "$DOWNLOAD_DIR" ]]; then
        for c_dir in "${DOWNLOAD_DIR}"/*; do
            [[ -d "$c_dir" ]] || continue
            c_name=$(basename "$c_dir")
            case "$c_name" in
                .*|work|logs|bin|tmp) continue ;;
            esac
            add_cohort_unique "$c_name"
        done
    fi
    if [[ -d "$PRE_CLIN_DIR" ]]; then
        for c_dir in "${PRE_CLIN_DIR}"/*; do
            [[ -d "$c_dir" ]] || continue
            c_name=$(basename "$c_dir")
            case "$c_name" in
                .*|work|logs|bin|tmp) continue ;;
            esac
            add_cohort_unique "$c_name"
        done
    fi
fi

echo "============================================================"
echo "Starting Clinical Re-Summarization"
echo "Cohorts:       ${COHORTS[*]:-<none>}"
echo "Clinical only: ${CLINICAL_ONLY}"
echo "============================================================"

for cohort in "${COHORTS[@]:-}"; do
    [[ -n "$cohort" ]] || continue
    COHORT_DOWNLOAD="${DOWNLOAD_DIR}/${cohort}"
    COHORT_PRE_CLIN="${PRE_CLIN_DIR}/${cohort}"
    COHORT_CBIO="${CBIO_OUT_DIR}/${cohort}"
    CLINICAL_TMP="${COHORT_CBIO}/clinical_tmp"
    IMPORT_DIR="${COHORT_CBIO}/import_folder"

    if [[ ! -d "$COHORT_DOWNLOAD" && ! -d "$COHORT_PRE_CLIN" ]]; then
        echo "Warning: Neither download dir (${COHORT_DOWNLOAD}) nor pre-clinical dir (${COHORT_PRE_CLIN}) found. Skipping."
        continue
    fi

    echo ""
    echo ">>> Processing Cohort: ${cohort}"
    [[ "$DRY_RUN" == "false" ]] && mkdir -p "$CLINICAL_TMP" "$IMPORT_DIR"

    # 1. Loop through samples and run process_clinical.py (unless clinical-only or download dir is missing/empty)
    sample_count=0
    COHORT_CLINICAL_ONLY="$CLINICAL_ONLY"
    LOCAL_PATIENT_DIRS=0
    if [[ -d "$COHORT_DOWNLOAD" ]]; then
        LOCAL_PATIENT_DIRS=$(find "$COHORT_DOWNLOAD" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [[ "$COHORT_CLINICAL_ONLY" != "true" ]] && [[ ! -d "$COHORT_DOWNLOAD" || "$LOCAL_PATIENT_DIRS" -eq 0 ]]; then
        echo "  GCP download folder (${COHORT_DOWNLOAD}) does not exist or is empty -> enabling --clinical-only automatically."
        COHORT_CLINICAL_ONLY=true
    fi

    if [[ "$COHORT_CLINICAL_ONLY" == "true" ]]; then
        echo "  [--clinical-only] Skipping per-sample genomic clinical extraction."
    else
        for PATIENT_DIR in "${COHORT_DOWNLOAD}"/*; do
            [[ -d "$PATIENT_DIR" ]] || continue
            PATIENT_ID=$(basename "$PATIENT_DIR")
            case "$PATIENT_ID" in
                .*|fingerprint_metrics|gene_expression_files|logs|orange_links|runtime_summary|pipeline_info|work)
                    continue
                    ;;
            esac

            SAMPLE_DIRS=()
            if [[ -d "${PATIENT_DIR}/purple" || -d "${PATIENT_DIR}/orange" ]]; then
                SAMPLE_DIRS+=("$PATIENT_DIR")
            else
                for D2 in "${PATIENT_DIR}"/*; do
                    [[ -d "$D2" ]] || continue
                    if [[ -d "${D2}/purple" || -d "${D2}/orange" ]]; then
                        SAMPLE_DIRS+=("$D2")
                    fi
                done
            fi

            for SAMPLE_DIR in "${SAMPLE_DIRS[@]}"; do
                TUMOR_ID=$(basename "$SAMPLE_DIR")
                NORMAL_ID=""

                SOMATIC_VCF=$(ls "${SAMPLE_DIR}/purple/"*.purple.somatic.vcf.gz 2>/dev/null | head -n 1 || true)
                if [[ -n "$SOMATIC_VCF" ]]; then
                    VCF_BASE=$(basename "$SOMATIC_VCF")
                    TUMOR_ID="${VCF_BASE%.purple.somatic.vcf.gz}"
                fi

                DETECTED_NORMAL=$(ls "${SAMPLE_DIR}/teal/"*.teal.tellength.tsv 2>/dev/null | grep -v "/${TUMOR_ID}\.teal\." | head -n 1 || true)
                if [[ -n "$DETECTED_NORMAL" ]]; then
                    NORMAL_ID=$(basename "$DETECTED_NORMAL" | sed 's/\.teal\.tellength\.tsv$//')
                fi

                ORANGE_JSON="${SAMPLE_DIR}/orange/${TUMOR_ID}.orange.json"
                [[ -f "$ORANGE_JSON" ]] || ORANGE_JSON=$(ls "${SAMPLE_DIR}/orange/"*.orange.json 2>/dev/null | head -n 1 || true)
                [[ -z "$ORANGE_JSON" ]] && ORANGE_JSON="NA"

                CLINICAL_OUT="${CLINICAL_TMP}/${TUMOR_ID}.clinical.tsv"

                cmd=(python3 "${SCRIPT_DIR}/process_clinical.py" "$TUMOR_ID" "$ORANGE_JSON" "$CLINICAL_OUT")
                [[ -n "$NORMAL_ID" ]] && cmd+=(--normal-id "$NORMAL_ID")

                if [[ "$DRY_RUN" == "true" ]]; then
                    echo "  [DRY RUN] SAMPLE_DIR=${SAMPLE_DIR} ${cmd[*]}"
                else
                    export SAMPLE_DIR
                    "${cmd[@]}"
                fi
                (( sample_count++ )) || true
            done
        done
        echo "  Processed clinical metrics for ${sample_count} sample(s)."
    fi

    # 2. Locate pre-generated clinical files
    PRE_PATIENT=""
    PRE_SAMPLE=""
    for base in "${PRE_CLIN_DIR}/${cohort}" "${PRE_CLIN_DIR}" "${IMPORT_DIR}"; do
        [[ -d "$base" ]] || continue
        if [[ -z "$PRE_PATIENT" ]]; then
            for fname in data_clinical_patient.txt data_clinical_patients.txt patient_clinical.tsv; do
                [[ -f "${base}/${fname}" ]] && PRE_PATIENT="${base}/${fname}" && break
            done
        fi
        if [[ -z "$PRE_SAMPLE" ]]; then
            for fname in data_clinical_sample.txt data_clinical_samples.txt sample_clinical.tsv; do
                [[ -f "${base}/${fname}" ]] && PRE_SAMPLE="${base}/${fname}" && break
            done
        fi
        [[ -n "$PRE_PATIENT" && -n "$PRE_SAMPLE" ]] && break
    done

    # 3. Merge clinical data
    SAMPLE_OUT="${IMPORT_DIR}/data_clinical_sample.txt"
    PATIENT_OUT="${IMPORT_DIR}/data_clinical_patient.txt"

    merge_cmd=(
        python3 "${SCRIPT_DIR}/merge_clinical.py"
        --tmp-dir "$CLINICAL_TMP"
        --out-sample "$SAMPLE_OUT"
        --out-patient "$PATIENT_OUT"
    )
    [[ -n "$PRE_SAMPLE" ]] && merge_cmd+=(--pre-sample "$PRE_SAMPLE")
    [[ -n "$PRE_PATIENT" ]] && merge_cmd+=(--pre-patient "$PRE_PATIENT")

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] ${merge_cmd[*]}"
    else
        echo "  Merging into ${SAMPLE_OUT}..."
        "${merge_cmd[@]}"
    fi

    # 4. Copy supplemental pre-clinical files (timelines, meta_study.txt, etc.) and ensure metadata exists
    STUDY_ID="$(echo "$cohort" | tr '[:upper:]' '[:lower:]' | tr '-' '_')_wgts"
    STUDY_NAME="${cohort} WGTS"
    STUDY_DESC="${STUDY_NAME} data processed with Oncoanalyser v3 and received from GCP bucket."

    if [[ "$DRY_RUN" == "false" ]]; then
        if [[ -d "$COHORT_PRE_CLIN" ]]; then
            SKIP_BASENAMES="data_clinical_patient.txt data_clinical_patients.txt data_clinical_sample.txt data_clinical_samples.txt patient_clinical.tsv sample_clinical.tsv"
            for f in "${COHORT_PRE_CLIN}"/data_*.txt "${COHORT_PRE_CLIN}"/meta_*.txt; do
                [[ -f "$f" ]] || continue
                bname=$(basename "$f")
                skip=false
                for skip_name in $SKIP_BASENAMES; do
                    [[ "$bname" == "$skip_name" ]] && skip=true && break
                done
                [[ "$skip" == "true" ]] && continue
                # Do not overwrite an existing meta_study.txt in IMPORT_DIR unless it is missing
                if [[ "$bname" == "meta_study.txt" && -f "${IMPORT_DIR}/meta_study.txt" ]]; then
                    continue
                fi
                cp "$f" "${IMPORT_DIR}/${bname}"
            done
        fi

        if [[ ! -f "${IMPORT_DIR}/meta_study.txt" ]]; then
            cat <<EOF > "${IMPORT_DIR}/meta_study.txt"
type_of_cancer: mixed
cancer_study_identifier: ${STUDY_ID}
name: ${STUDY_NAME}
description: ${STUDY_DESC}
reference_genome: hg38
add_global_case_list: true
EOF
        fi

        python3 "${SCRIPT_DIR}/sanitize_clinical_timeline.py" "$IMPORT_DIR" --apply
        python3 "${SCRIPT_DIR}/generate_missing_meta.py" "$IMPORT_DIR" "$STUDY_ID"

        if [[ ! -f "$PATIENT_OUT" && -f "${IMPORT_DIR}/meta_clinical_patient.txt" ]]; then
            rm -f "${IMPORT_DIR}/meta_clinical_patient.txt"
        fi
        chmod -R g+rwX "$COHORT_CBIO" 2>/dev/null || true
    fi

    echo "  Cohort ${cohort} complete!"
done

echo ""
echo "============================================================"
echo "All done! Clinical files successfully re-summarized."
echo "============================================================"
