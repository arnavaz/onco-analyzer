#!/bin/bash
# ---------------------------------------------------------------------------
# run_mohccn_p7.sh
#
# Converts Oncoanalyser WGTS outputs for the MOHCCN Phase 7 cohort
# (downloaded from gs://hmf-mohccn-p1-reanalysis-output/) into a
# cBioPortal-ready study bundle.
#
# Workflow:
#   1. Download new patients first (if not already done):
#         bash oa_to_cbio/scripts/download_from_gcp.sh --cohort MOHCCN-P7
#      OR to auto-detect all cohorts and download + convert in one step:
#         bash oa_to_cbio/scripts/sync_from_gcp.sh --cohort MOHCCN-P7
#
#   2. Run this script to convert already-downloaded data:
#         bash runs/run_mohccn_p7.sh [--merge-only]
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OA_TO_CBIO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --- Paths ------------------------------------------------------------------
COHORT="MOHCCN-P7"
# Downloaded patients are routed to per-cohort subdirs by download_from_gcp.sh
DOWNLOAD_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download/${COHORT}"
OUTPUT_DIR="/cluster/projects/mohccn/oncoanalyser_cbioportal/cbio_import/${COHORT}"
SLURM_ACCOUNT="deperrotgroup"

# Pre-generated clinical files: check cohort-specific dir first, then shared root
PRE_CLIN_ROOT="/cluster/projects/mohccn/oncoanalyser_cbioportal/pre_gen_clin_data"
PRE_PATIENT_CLINICAL="${PRE_CLIN_ROOT}/${COHORT}/patient_clinical.tsv"
PRE_SAMPLE_CLINICAL="${PRE_CLIN_ROOT}/${COHORT}/sample_clinical.tsv"
# Fallback to shared files if cohort-specific ones don't exist
[[ ! -f "$PRE_PATIENT_CLINICAL" ]] && PRE_PATIENT_CLINICAL="${PRE_CLIN_ROOT}/patient_clinical.tsv"
[[ ! -f "$PRE_SAMPLE_CLINICAL"  ]] && PRE_SAMPLE_CLINICAL="${PRE_CLIN_ROOT}/sample_clinical.tsv"

# --- Study Metadata ---------------------------------------------------------
STUDY_ID="mohccn_p7_wgts"
STUDY_NAME="MOHCCN Phase 7 WGTS"
STUDY_DESC="MOHCCN whole-genome and transcriptome data processed with Oncoanalyser v3. Received via GCP bucket and converted for cBioPortal."
TYPE_OF_CANCER="mixed"

MERGE_ONLY="false"
if [[ "${1:-}" == "--merge-only" ]]; then
    MERGE_ONLY="true"
fi

echo "========================================================"
echo "Starting MOHCCN P7 cBioPortal conversion"
echo "  Input:        ${DOWNLOAD_DIR}"
echo "  Output:       ${OUTPUT_DIR}"
echo "  SLURM Account:${SLURM_ACCOUNT}"
echo "  Study ID:     ${STUDY_ID}"
echo "  Merge only:   ${MERGE_ONLY}"
echo "========================================================"

mkdir -p "$OUTPUT_DIR"

# Build the command
CMD=(
    bash "${OA_TO_CBIO_DIR}/oa_to_cbio.sh"
    -i "$DOWNLOAD_DIR"
    -o "$OUTPUT_DIR"
    -A "$SLURM_ACCOUNT"
    -s "$STUDY_ID"
    -n "$STUDY_NAME"
    -d "$STUDY_DESC"
    -t "$TYPE_OF_CANCER"
)

# Add pre-generated clinical enrichment if files exist
if [[ -f "$PRE_PATIENT_CLINICAL" ]]; then
    CMD+=(-P "$PRE_PATIENT_CLINICAL")
    echo "  Using pre-patient clinical: ${PRE_PATIENT_CLINICAL}"
fi
if [[ -f "$PRE_SAMPLE_CLINICAL" ]]; then
    CMD+=(-S "$PRE_SAMPLE_CLINICAL")
    echo "  Using pre-sample clinical:  ${PRE_SAMPLE_CLINICAL}"
fi

if [[ "$MERGE_ONLY" == "true" ]]; then
    CMD+=(--merge-only)
fi

echo ""
"${CMD[@]}"

echo "----------------------------------------"
echo "Job submission complete. Logs: ${OUTPUT_DIR}/logs/"
echo "Once all jobs finish, validate with:"
echo "  cbio_validator.sh ${OUTPUT_DIR}/import_folder"
echo "========================================================"
