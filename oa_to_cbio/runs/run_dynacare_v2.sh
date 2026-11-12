#!/bin/bash
# ---------------------------------------------------------------------------
# run_dynacare_v2.sh
# Runs the oa_to_cbio pipeline to convert nf-core/oncoanalyser outputs
# of the Dynacare_v2 cohort into a cBioPortal-ready study bundle.
# ---------------------------------------------------------------------------

set -euo pipefail

# Determine script directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OA_TO_CBIO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Input and Output Directories
INPUT_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results"
OUTPUT_DIR="${INPUT_DIR}/cbio_import"
SLURM_ACCOUNT="deperrotgroup"

# cBioPortal Study Metadata Configuration
STUDY_ID="multimmr_hmf_dynacare_v2"
STUDY_NAME="MultiMMR HMF Dynacare v2"
STUDY_DESC="Dynacare cohort captured using v2 of the MultiMMR panel and processed with nf-core/oncoanalyser"
TYPE_OF_CANCER="mixed"

echo "========================================"
echo "Starting Oncoanalyser to cBioPortal pipeline"
echo "========================================"
echo "Input Dir:      ${INPUT_DIR}"
echo "Output Dir:     ${OUTPUT_DIR}"
echo "SLURM Account:  ${SLURM_ACCOUNT}"
echo "Study ID:       ${STUDY_ID}"
echo "Study Name:     ${STUDY_NAME}"
echo "========================================"

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Run the master orchestrator script
bash "${OA_TO_CBIO_DIR}/oa_to_cbio.sh" \
  -i "${INPUT_DIR}" \
  -o "${OUTPUT_DIR}" \
  -A "${SLURM_ACCOUNT}" \
  -s "${STUDY_ID}" \
  -n "${STUDY_NAME}" \
  -d "${STUDY_DESC}" \
  -t "${TYPE_OF_CANCER}" \
  -g "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/final_panel_files/LOH_plot_genes.tsv" \
  --unfiltered-germline-maf

echo "----------------------------------------"
echo "cBioPortal import pipeline job submission complete."
echo "Check output logs in ${OUTPUT_DIR}/logs"
echo "========================================"
