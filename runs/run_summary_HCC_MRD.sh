#!/bin/bash
# ==============================================================================
# Interactive Runner: HCC WISP WGTS Summary & cBioPortal Imports
#
# Note: Designed as interactive blocks to run directly on the cluster.
# ==============================================================================

SOURCE_REPO="/cluster/projects/pughlab/src/oncoanalyser_pughlab"
export MPLCONFIGDIR="/tmp/${USER}_matplotlib"
mkdir -p "$MPLCONFIGDIR"


## ==============================================================================
## BLOCK 0: Inspect Partition & Node Availability (Find cores & memory)
## ==============================================================================

# 1. Find "mixed" nodes with >= 10 idle cores and >= 40 GB free RAM right now:
# sinfo -N -h -p pughlab,himem,superhimem,long -o "%P %N %T %C %e" | awk -F'[/ ]+' '$5 >= 10 && $8 >= 40000 {printf "%-12s %-12s State: %-5s Idle Cores: %-3s Free Mem: %5d GB\n", $1, $2, $3, $5, int($8/1024)}'

## ==============================================================================
## BLOCK 1: Summarize Cohort WGTS/WISP MRD Results: BATCH 1 ONLY
## ==============================================================================
# Note: Always pass -A pughlab explicitly to avoid default account mismatch
salloc -A pughlab -p himem -c 10 --time=10:00:00 --mem=40G

module load python3/3.10.9

python3 /cluster/projects/pughlab/src/oncoanalyser_pughlab/scripts/compile_wgts_summary.py \
  --results-dir /cluster/projects/pughlab/hcc/external_data/batch_1 \
  --output-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/Oncoanalyser/summary_batch_1 \
  --timepoints-tsv /cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/HCC_metadata.tsv \
  --cohort-name "HCC_WISP_B1" \
  --genes-of-interest TP53,CTNNB1,AXIN1,ARID1A,ARID2,BAP1,TERT,RB1,CDKN2A,PIK3CA,NFE2L2,KEAP1,ALB,APOB,PTEN \
  --fingerprint-stride 20 \
  --threads 10


## ==============================================================================
## BLOCK 2: Stage Merged Batch Outputs (Symlink Batches 1 & 2)
## ==============================================================================
MERGED_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/merged_OA_outputs"
mkdir -p "$MERGED_DIR"

# Symlink sample folders while filtering out pipeline logs and internal Nextflow folders
for batch in /cluster/projects/pughlab/hcc/external_data/batch_1 /cluster/projects/pughlab/hcc/external_data/batch_2; do
  for d in "$batch"/*; do
    [ -d "$d" ] || continue
    base=$(basename "$d")
    case "$base" in
      pipeline_info|work|logs|runtime_summary|fingerprint_metrics) continue ;;
      *) ln -sfn "$d" "$MERGED_DIR/" ;;
    esac
  done
done


## ==============================================================================
## BLOCK 3: Summarize Merged WGTS/WISP Cohort (Batches 1 & 2 Combined)
## ==============================================================================
# Allocate interactive node (use pughlab or himem if pughlab is busy; 40G is sufficient)
salloc -A pughlab -p pughlab -c 10 --time=10:00:00 --mem=40G
salloc -A deperrotgroup -p himem -c 10 --time=10:00:00 --mem=60G

module load python3/3.10.9

python3 /cluster/projects/pughlab/src/oncoanalyser_pughlab/scripts/compile_wgts_summary.py \
  --results-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/merged_OA_outputs \
  --output-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/Oncoanalyser/summary \
  --timepoints-tsv /cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/HCC_metadata.tsv \
  --cohort-name "HCC_WISP" \
  --genes-of-interest TP53,CTNNB1,AXIN1,ARID1A,ARID2,BAP1,TERT,RB1,CDKN2A,PIK3CA,NFE2L2,KEAP1,ALB,APOB,PTEN \
  --fingerprint-stride 20 \
  --threads 10


## ==============================================================================
## BLOCK 4: Generate cBioPortal Import Bundle
## ==============================================================================
# Run from login node (or interactive node). oa_to_cbio dispatches SLURM jobs per sample.
OA_TO_CBIO_DIR="/cluster/projects/pughlab/src/oncoanalyser_pughlab/oa_to_cbio"
INPUT_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/merged_OA_outputs"
OUTPUT_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/HCC_WISP/Oncoanalyser/cbio_import"
SLURM_ACCOUNT="pughlab"
STUDY_ID="hcc_wisp"
STUDY_NAME="HCC WISP"
STUDY_DESC="Comprehensive molecular profiles for HCC plasma WGS cohort. Primary and relapse samples collected at baseline and serial timepoints."
TYPE_OF_CANCER="HCC"

mkdir -p "${OUTPUT_DIR}"

bash "${OA_TO_CBIO_DIR}/oa_to_cbio.sh" \
  -i "${INPUT_DIR}" \
  -o "${OUTPUT_DIR}" \
  -A "${SLURM_ACCOUNT}" \
  -s "${STUDY_ID}" \
  -n "${STUDY_NAME}" \
  -d "${STUDY_DESC}" \
  -t "${TYPE_OF_CANCER}"

