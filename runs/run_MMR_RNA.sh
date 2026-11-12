#!/bin/bash
# ---------------------------------------------------------------------------
# Oncoanalyser RNA-Seq Only Pipeline Runner (MMR Cohort)
#
# Generates per-sample Oncoanalyser samplesheets from MMR_RNAseq.tsv and
# MMR_RNA_fastq_list.txt, then submits SLURM batch jobs for STAR alignment
# and Isofox RNA quantification in WGTS mode.
# ---------------------------------------------------------------------------

set -euo pipefail

# Environment & Accounts
module load python3/3.10.9 2>/dev/null || module load python3 2>/dev/null || true

ACCOUNT="${ACCOUNT:-pughlab}"
PARTITION="${PARTITION:-all}"
SOURCE_REPO="${SOURCE_REPO:-/cluster/projects/pughlab/src/oncoanalyser_pughlab}"

# Input & Output Paths
FASTQ_DIR="${FASTQ_DIR:-/cluster/projects/pughlab/hereditary/external_data/upload}"
FASTQ_PATH_LIST="${FASTQ_PATH_LIST:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/MMR_RNA/sample_info/MMR_RNA_fastq_list.txt}"
OUTPUT_DIR="${OUTPUT_DIR:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/MMR_RNA/Oncoanalyser}"
UNFORMATED_SAMPLE_SHEET="${UNFORMATED_SAMPLE_SHEET:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/MMR_RNA/sample_info/MMR_RNAseq.tsv}"
SHEETS_DIR="${OUTPUT_DIR}/sample_sheets"

echo "======================================================================"
echo "MMR RNA-Seq Only Pipeline Setup"
echo "  Source Repo:     ${SOURCE_REPO}"
echo "  Account:         ${ACCOUNT}"
echo "  Partition:       ${PARTITION}"
echo "  FASTQ List:      ${FASTQ_PATH_LIST}"
echo "  Metadata TSV:    ${UNFORMATED_SAMPLE_SHEET}"
echo "  Output Dir:      ${OUTPUT_DIR}"
echo "  Sample Sheets:   ${SHEETS_DIR}"
echo "======================================================================"

# Step 1: Generate per-sample Oncoanalyser samplesheets
echo ""
echo "Step 1: Generating RNA samplesheets..."
python3 "${SOURCE_REPO}/Setup/make_rna_samplesheet.py" \
  --fastq-list "${FASTQ_PATH_LIST}" \
  --metadata-tsv "${UNFORMATED_SAMPLE_SHEET}" \
  --out "${SHEETS_DIR}"

# Step 2: Configure offline Nextflow & Singularity environment
export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

# Setup user-controlled temporary directory to prevent compute node /tmp exhaustion
export TMPDIR="${OUTPUT_DIR}/tmp"
export TEMP="${OUTPUT_DIR}/tmp"
export TMP="${OUTPUT_DIR}/tmp"
export SINGULARITY_TMPDIR="${OUTPUT_DIR}/tmp"
export APPTAINER_TMPDIR="${OUTPUT_DIR}/tmp"
export SINGULARITYENV_TMPDIR="${OUTPUT_DIR}/tmp"
export SINGULARITYENV_TEMP="${OUTPUT_DIR}/tmp"
export SINGULARITYENV_TMP="${OUTPUT_DIR}/tmp"
mkdir -p "${TMPDIR}"

# Step 3: Submit Oncoanalyser jobs per sample
echo ""
echo "Step 2: Submitting Oncoanalyser RNA runs..."

SUBMITTED_COUNT=0
SKIPPED_COUNT=0

for sample_sheet in "${SHEETS_DIR}"/*/samplesheet.csv; do
  [[ -f "$sample_sheet" ]] || continue
  sample_name=$(basename "$(dirname "$sample_sheet")")
  sample_outdir="${OUTPUT_DIR}/${sample_name}"
  sample_tmpdir="${sample_outdir}/tmp"

  mkdir -p "${sample_outdir}/logs" "${sample_tmpdir}"

  # 1. Skip if already completed (isofox gene data output exists)
  if compgen -G "${sample_outdir}/isofox/*.isf.gene_data.*" > /dev/null || \
     compgen -G "${sample_outdir}/*/isofox/*.isf.gene_data.*" > /dev/null; then
    echo "  [${sample_name}] already completed; skipping."
    ((SKIPPED_COUNT++)) || true
    continue
  fi

  # 2. Skip if already running or queued in SLURM
  running_job=$(squeue -u "$USER" --name="onco_rna_${sample_name}" -h -o "%i" | head -n 1)
  if [[ -n "$running_job" ]]; then
    echo "  [${sample_name}] already running in SLURM (Job ID: ${running_job}); skipping."
    ((SKIPPED_COUNT++)) || true
    continue
  fi

  echo "  Submitting RNA run for: ${sample_name}..."
  job_id=$(sbatch --parsable \
    --partition="$PARTITION" \
    --time=24:00:00 \
    --cpus-per-task=1 \
    --mem=6G \
    --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,TMPDIR=${sample_tmpdir},TEMP=${sample_tmpdir},TMP=${sample_tmpdir},SINGULARITY_TMPDIR=${sample_tmpdir},APPTAINER_TMPDIR=${sample_tmpdir},SINGULARITYENV_TMPDIR=${sample_tmpdir},SINGULARITYENV_TEMP=${sample_tmpdir},SINGULARITYENV_TMP=${sample_tmpdir} \
    -A "$ACCOUNT" \
    --job-name="onco_rna_${sample_name}" \
    --output="${sample_outdir}/logs/onco_rna_${sample_name}_%j.out" \
    --error="${sample_outdir}/logs/onco_rna_${sample_name}_%j.err" \
    "${SOURCE_REPO}/run_oncoanalyser.sh" \
    "$sample_sheet" \
    "$sample_outdir" \
    wgts \
    GRCh38_hmf \
    "$ACCOUNT" \
    --metadata-tsv "${UNFORMATED_SAMPLE_SHEET}" \
    --force-nextflow)
  job_id=$(echo "$job_id" | cut -d';' -f1)
  echo "    -> Submitted Job ID: ${job_id} (${sample_name})"
  ((SUBMITTED_COUNT++)) || true
done

echo ""
echo "======================================================================"
echo "Submission Summary:"
echo "  Submitted: ${SUBMITTED_COUNT}"
echo "  Skipped:   ${SKIPPED_COUNT}"
echo ""
echo "When runs finish, generate the interactive HTML dashboard with:"
echo "module load python3/3.10.9"
echo "  python3 ${SOURCE_REPO}/scripts/compile_mmr_rna_dashboard.py \\"
echo "    --results-dir ${OUTPUT_DIR} \\"
echo "    --metadata-tsv ${UNFORMATED_SAMPLE_SHEET}"
echo "======================================================================"