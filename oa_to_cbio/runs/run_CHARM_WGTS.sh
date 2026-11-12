#!/bin/bash
# ---------------------------------------------------------------------------
# run_CHARM_WGTS.sh
# Runs the oa_to_cbio pipeline to convert nf-core/oncoanalyser outputs
# of the CHARM WGTS cohort into a cBioPortal-ready study bundle.
# ---------------------------------------------------------------------------

set -euo pipefail

# Determine script directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OA_TO_CBIO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Input and Output Directories
INPUT_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser"
OUTPUT_DIR="${INPUT_DIR}/cbio_import"
SLURM_ACCOUNT="deperrotgroup"

# cBioPortal Study Metadata Configuration
STUDY_ID="charm_wgts"
STUDY_NAME="CHARM WGTS"
STUDY_DESC="CHARM WGTS Case processed with nf-core/oncoanalyser"
TYPE_OF_CANCER="mixed"

MERGE_ONLY="false"
if [[ "${1:-}" == "--merge-only" ]]; then
    MERGE_ONLY="true"
fi

echo "========================================"
echo "Starting Oncoanalyser to cBioPortal pipeline for CHARM WGTS"
echo "========================================"
echo "Input Dir:      ${INPUT_DIR}"
echo "Output Dir:     ${OUTPUT_DIR}"
echo "SLURM Account:  ${SLURM_ACCOUNT}"
echo "Study ID:       ${STUDY_ID}"
echo "Study Name:     ${STUDY_NAME}"
echo "Merge only:     ${MERGE_ONLY}"
echo "========================================"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Build the command
CMD=(
    bash "${OA_TO_CBIO_DIR}/oa_to_cbio.sh"
    -i "${INPUT_DIR}"
    -o "${OUTPUT_DIR}"
    -A "${SLURM_ACCOUNT}"
    -s "${STUDY_ID}"
    -n "${STUDY_NAME}"
    -d "${STUDY_DESC}"
    -t "${TYPE_OF_CANCER}"
)

if [[ "$MERGE_ONLY" == "true" ]]; then
    CMD+=(--merge-only)
fi

# Run the master orchestrator script
"${CMD[@]}"

echo "----------------------------------------"
echo "cBioPortal import pipeline job submission complete."
echo "Check output logs in ${OUTPUT_DIR}/logs"
echo "========================================"
