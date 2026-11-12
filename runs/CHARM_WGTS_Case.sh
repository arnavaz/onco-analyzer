#!/bin/bash

module load python3/3.12.11

OUTDIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser
ACCOUNT=deperrotgroup

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab


FASTQ_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/all_fastq_links
MAPPING_TSV=/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Unarchived_fastqs/CHARMWGTS.release.2026-06-11_10:37.fastqs.map.tsv
SHEETS_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/sample_sheets

echo "Generating sample sheets..."
python3 "${SOURCE_REPO}/Setup/make_oncoanalyser_samplesheet.py" \
  --fastq-dir "$FASTQ_DIR" \
  --mapping-tsv "$MAPPING_TSV" \
  --out "$SHEETS_DIR"

export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
  [[ -f "$sample_sheet" ]] || continue
  patient=$(basename "$(dirname "$sample_sheet")")
  patient_outdir="${OUTDIR}/${patient}"
  
  mkdir -p "${patient_outdir}/logs"
  
  echo "Submitting oncoanalyser job for: ${patient}"
  sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true -A "$ACCOUNT" --partition=all --time=04:00:00 \
    --job-name="onco_${patient}" \
    --output="${patient_outdir}/logs/onco_${patient}_%j.out" \
    --error="${patient_outdir}/logs/onco_${patient}_%j.err" \
    "${SOURCE_REPO}/run_oncoanalyser.sh" \
    "$sample_sheet" \
    "$patient_outdir" \
    wgts \
    GRCh38_hmf \
    "$ACCOUNT" \
    --ichor-pon /cluster/projects/pughlab/references/panel_of_normals/hg38/IchorCNA/pon_ichorCNA__CHARM__median.rds \
    --ichorWGS \
    --skip-nextflow
done


## Once done, move on to Ultima plasm step using --mode purity_estimate

## Run again in tumour only mode
T_ONLY_SHEET_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/sample_sheets/CHM2-01-0238_T_only"
mkdir -p "${T_ONLY_SHEET_DIR}"
cat <<EOF > "${T_ONLY_SHEET_DIR}/samplesheet.csv"
group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
CHM2-01-0238,CHM2-01-0238,S25-22052,tumor,dna,cram,,/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238/S25-22052/alignments/dna/S25-22052.redux.cram
EOF

T_ONLY_OUTDIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238_T_only"
patient_outdir="${T_ONLY_OUTDIR}/CHM2-01-0238"
mkdir -p "${patient_outdir}/logs"

echo "Submitting oncoanalyser tumor-only job for: CHM2-01-0238"
sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true -A "$ACCOUNT" \
  --job-name="onco_CHM2-01-0238_T_only" \
  --output="${patient_outdir}/logs/onco_CHM2-01-0238_%j.out" \
  --error="${patient_outdir}/logs/onco_CHM2-01-0238_%j.err" \
  "${SOURCE_REPO}/run_oncoanalyser.sh" \
  "${T_ONLY_SHEET_DIR}/samplesheet.csv" \
  "${patient_outdir}" \
  wgts \
  GRCh38_hmf \
  "$ACCOUNT" \
  --ichor-pon /cluster/projects/pughlab/references/panel_of_normals/hg38/IchorCNA/pon_ichorCNA__CHARM__median.rds \
  --ichorWGS


## Ultima data

cram_dir=/cluster/projects/pughlab/hereditary/external_data/CHARM2UGpos_plasma080626/merged_crams
