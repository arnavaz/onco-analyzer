#!/bin/bash
# ==============================================================================
# create_clean_repo.sh
# Creates a clean, streamlined repository containing only the necessary files
# to run the MultiMMR panel workflow, omitting all extraneous code/scripts.
# ==============================================================================
set -euo pipefail

# Print usage
usage() {
    echo "Usage: $0 <destination_directory>"
    echo "Example: $0 /path/to/oncoanalyser_multimmr"
    exit 1
}

if [[ $# -ne 1 ]]; then
    usage
fi

DEST_DIR="$1"
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Verify the script is run from inside the correct source repository
if [[ ! -f "${SOURCE_DIR}/run_oncoanalyser.sh" ]]; then
    echo "Error: Must run this script from the root of the oncoanalyser_pughlab repository." >&2
    exit 1
fi

echo ">>> Creating clean repository directory at: ${DEST_DIR}"
mkdir -p "${DEST_DIR}"

echo ">>> Initializing new Git repository..."
git -C "${DEST_DIR}" init

echo ">>> Copying Nextflow pipeline..."
mkdir -p "${DEST_DIR}/pipeline"
cp -R "${SOURCE_DIR}/pipeline/" "${DEST_DIR}/pipeline/"

echo ">>> Copying configuration files..."
mkdir -p "${DEST_DIR}/conf"
cp "${SOURCE_DIR}/conf/h4h.config" "${DEST_DIR}/conf/"
cp "${SOURCE_DIR}/conf/signatures_etiology.tsv" "${DEST_DIR}/conf/"

echo ">>> Copying samplesheet setup helpers..."
mkdir -p "${DEST_DIR}/Setup"
cp "${SOURCE_DIR}/Setup/make_oncoanalyser_samplesheet.py" "${DEST_DIR}/Setup/"
cp "${SOURCE_DIR}/Setup/samplesheet_template.csv" "${DEST_DIR}/Setup/"
cp "${SOURCE_DIR}/Setup/injection_sample_information_template.txt" "${DEST_DIR}/Setup/"

echo ">>> Copying packaged v3 panel resources..."
mkdir -p "${DEST_DIR}/resources/panels/multimmr_v3"
cp -R "${SOURCE_DIR}/resources/panels/multimmr_v3/" "${DEST_DIR}/resources/panels/multimmr_v3/"

echo ">>> Copying pipeline submission wrappers..."
cp "${SOURCE_DIR}/run_oncoanalyser.sh" "${DEST_DIR}/"
cp "${SOURCE_DIR}/submit_all_patients.sh" "${DEST_DIR}/"

echo ">>> Copying mismatch repair and resource-generation scripts..."
mkdir -p "${DEST_DIR}/scripts"
cp "${SOURCE_DIR}/scripts/compile_mmr_summary.py" "${DEST_DIR}/scripts/"
cp "${SOURCE_DIR}/scripts/summarize_findings.py" "${DEST_DIR}/scripts/"
cp "${SOURCE_DIR}/scripts/generate_panel_snps.py" "${DEST_DIR}/scripts/"
cp "${SOURCE_DIR}/scripts/generate_msi_panel.py" "${DEST_DIR}/scripts/"
cp "${SOURCE_DIR}/scripts/package_mmr_reports.py" "${DEST_DIR}/scripts/"
cp "${SOURCE_DIR}/scripts/plot_gene_baf.py" "${DEST_DIR}/scripts/"

echo ">>> Setting up README.md in the new repository..."
cp "${SOURCE_DIR}/README_MultiMMR.md" "${DEST_DIR}/README.md"

echo ">>> Generating clean .gitignore..."
cat << 'EOF' > "${DEST_DIR}/.gitignore"
.DS_Store
.nextflow*
work/
nextflow.log*
.nextflow.log*
*.sif
*.img
hmf_temp/
logs/nf_*
logs/extract_*
EOF

echo ">>> Creating initial Git commit..."
git -C "${DEST_DIR}" add .
git -C "${DEST_DIR}" commit -m "Initial commit of clean MultiMMR panel workflow repository"

echo "========================================================================="
echo " Clean repository successfully created at: ${DEST_DIR}"
echo "========================================================================="
echo "To publish this repository to GitHub:"
echo "  1. Change directory: cd ${DEST_DIR}"
echo "  2. Run GitHub CLI command:"
echo "     gh repo create oncoanalyser_multimmr --private --source=. --remote=origin --push"
echo "  3. Or push manually:"
echo "     git remote add origin git@github.com:<YOUR-ORG-OR-USER>/oncoanalyser_multimmr.git"
echo "     git branch -M main"
echo "     git push -u origin main"
echo "========================================================================="
