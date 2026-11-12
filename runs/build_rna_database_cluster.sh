#!/bin/bash
#SBATCH -J build_rna_database
#SBATCH -p pughlab
#SBATCH -c 32
#SBATCH --mem=250G
#SBATCH -t 10:00:00
#SBATCH -o /cluster/projects/mohccn/oncoanalyser_cbioportal/logs/build_rna_db_%j.out
#SBATCH -e /cluster/projects/mohccn/oncoanalyser_cbioportal/logs/build_rna_db_%j.err

# ---------------------------------------------------------------------------
# build_rna_database_cluster.sh
#
# SLURM batch runner to build or rebuild the full multi-cohort DuckDB
# analytical database from Oncoanalyser / Isofox RNA-Seq runs in gcp_download
# and cBioPortal clinical metadata in cbio_import.
#
# Submit to SLURM:
#   sbatch [--account pughlab] runs/build_rna_database_cluster.sh
#
# Or run interactively within salloc allocation:
#   bash runs/build_rna_database_cluster.sh
# ---------------------------------------------------------------------------

set -euo pipefail
umask 0002

# Environment
module load python3/3.10.9 2>/dev/null || module load python3 2>/dev/null || true

SOURCE_REPO="${SOURCE_REPO:-/cluster/projects/pughlab/src/oncoanalyser_pughlab}"
GCP_DOWNLOAD_DIR="${GCP_DOWNLOAD_DIR:-/cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download}"
CBIO_IMPORT_DIR="${CBIO_IMPORT_DIR:-/cluster/projects/mohccn/oncoanalyser_cbioportal/cbio_import}"
OUTPUT_DB="${OUTPUT_DB:-/cluster/projects/mohccn/oncoanalyser_cbioportal/hmf_rnaseq_cohort.duckdb}"
LOGS_DIR="${LOGS_DIR:-/cluster/projects/mohccn/oncoanalyser_cbioportal/logs}"

mkdir -p "${LOGS_DIR}"
cd "${SOURCE_REPO}"

echo "======================================================================"
echo "HMF RNA-Seq Cohort Analytical Database Rebuild"
echo "  Start Time:  $(date)"
echo "  Source Repo: ${SOURCE_REPO}"
echo "  Results:     ${GCP_DOWNLOAD_DIR}/*"
echo "  Clinical:    ${CBIO_IMPORT_DIR}/*"
echo "  Output DB:   ${OUTPUT_DB}"
echo "  Threads:     28"
echo "======================================================================"

python3 scripts/build_rna_database.py \
    --results-dir "${GCP_DOWNLOAD_DIR}"/* \
    --cbio-dir "${CBIO_IMPORT_DIR}"/* \
    --output-db "${OUTPUT_DB}" \
    --batch-column BATCH_ID \
    --threads 28

echo "======================================================================"
echo "Database build complete: $(date)"
ls -lh "${OUTPUT_DB}"
echo ""
echo "Next step: Copy updated database to local laptop for RNA Studio:"
echo "  scp jbruce@h4huhndata1:${OUTPUT_DB} ~/Downloads/"
echo "======================================================================"
