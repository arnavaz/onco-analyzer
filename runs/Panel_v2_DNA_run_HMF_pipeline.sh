## Generate SNP panel augmentation

conda activate nf-core
#run locally
python3 scripts/generate_panel_snps.py \
  --target-regions-bed /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/target_regions_definition.38.bed.gz \
  --default-loci /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/hmf_pipeline_resources.38_v2.3.0--2/dna/copy_number/AmberGermlineSites.38.tsv.gz \
  --output /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/AmberGermlineSites.expanded.38.tsv.gz \
  --genes all_bed \
  --genome GRCh38_hmf

## Generate panel resources

# 1. Export caching and offline environment variables
export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true
# 2. Submit the panel resource creation job
sbatch --export=ALL,SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab,NXF_OFFLINE=true,RESUME=true -A deperrotgroup \
  /cluster/projects/pughlab/src/oncoanalyser_pughlab/run_oncoanalyser.sh \
  /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/sample_sheet/samplesheet.panel_resource_creation.csv \
  /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/panel_resource \
  panel_resource_creation \
  GRCh38_hmf \
  deperrotgroup \
  --driver-gene-panel /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/MMR_drivers.tsv \
  --target-regions-bed /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/target_regions_definition.38.bed.gz


## Generate off-target ichorCNA Panel of Normals (PoN) and Exclusion BED
salloc -c 10 -p superhimem -t 10:0:0 --mem 120G
conda activate nf-core

module load hmmcopy_utils/170718
module load R/4.1.0

python3 /cluster/projects/pughlab/src/oncoanalyser_pughlab/scripts/generate_ichorcna_panel_resources.py \
  --bams /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/PON_links/*.bam \
  --chrom-sizes /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/hg38_genome.chrom.sizes \
  --output-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources \
  --mosdepth-bin /cluster/projects/pughlab/bin/mosdepth/mosdepth \
  --padding 250 \
  --threads 10

## Submit all patients with updated germline modifications and correct panel files
conda activate nf-core

PANEL_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/final_panel_files
PANEL_NAME=multimmr_hmf
PANEL_CONFIG=$PANEL_DIR/panel_data.config

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab
SHEETS_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results
RESULTS_BASE=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results
ACCOUNT=deperrotgroup
UHN_ICHOR_DIR=/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/sWGS/IchorCNA/
OTB_ICHOR_DIR=/cluster/projects/pughlab/pughlab/projects/MultiMMR/OTB/sWGS/IchorCNA/
SKIP_SAMPLE=NONE

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
  sample_name=$(basename "$(dirname "$sample_sheet")")

  if [[ "$sample_name" == "$SKIP_SAMPLE" ]]; then
    echo "Skipping already-submitted sample: $sample_name"
    continue
  fi

  # Strip prefix if present (e.g. DYNA- or DYNA_)
  clean_name="${sample_name#DYNA-}"
  clean_name="${clean_name#DYNA_}"
  
  # Extract patient ID (e.g. UHN-27682)
  patient_id=$(echo "${clean_name}" | grep -oE "[A-Z]+-[0-9]+" | head -n 1 || echo "${clean_name}")
  
  ICHOR_DIR=""
  # Search across all specified parent directories
  for parent_dir in "${UHN_ICHOR_DIR}" "${OTB_ICHOR_DIR}"; do
      if [[ -d "${parent_dir}" ]]; then
          # 1. Try matching the nested structure patient_id/patient_id-T
          candidate_dir="${parent_dir}/${patient_id}/${patient_id}-T"
          
          # 2. Fallback: try matching on patient_id folder
          if [[ ! -d "${candidate_dir}" ]]; then
              candidate_dir="${parent_dir}/${patient_id}"
          fi
          
          # 3. Fallback: try matching clean_name/clean_name-T
          if [[ ! -d "${candidate_dir}" ]]; then
              candidate_dir="${parent_dir}/${clean_name}/${clean_name}-T"
          fi
          
          # 4. Fallback: try matching clean_name directly
          if [[ ! -d "${candidate_dir}" ]]; then
              candidate_dir="${parent_dir}/${sample_name}"
          fi
          
          # 5. Fallback: 2-level deep fuzzy find for patient_id containing "-T"
          if [[ ! -d "${candidate_dir}" ]]; then
              candidate_dir=$(find "${parent_dir}" -maxdepth 2 -type d -name "*${patient_id}*" 2>/dev/null | grep -E -- "-T$" | head -n 1)
          fi
          
          # If a valid directory was found on disk, use it and break out of parent loop
          if [[ -n "${candidate_dir}" && -d "${candidate_dir}" ]]; then
              ICHOR_DIR="${candidate_dir}"
              break
          fi
      fi
  done

  # Check if a valid directory was found on disk
  if [[ -n "${ICHOR_DIR}" && -d "${ICHOR_DIR}" ]]; then
      echo "Found Ichor directory: $ICHOR_DIR"
  else
      echo "WARNING: No matching sWGS ichorCNA directory found for sample $sample_name. Running without sWGS fallback."
      ICHOR_DIR=""
  fi

  PATIENT_CODE=$(echo $sample_name | sed 's/.*-//')

  ENDOMETRIAL_LIST="163243 39096 76452 41463 714363 714364 714576 714823 714834 714960 715035 715042 715465 715620 715628 715779 715900 715916 716627 716666 716757"
  
  TISSUE="colorectal"
  for i in $ENDOMETRIAL_LIST; do
    if [[ "$PATIENT_CODE" == "$i" ]]; then
      TISSUE="endometrial"
      break
    fi
  done
  echo "TISSUE: $TISSUE"
  
  outdir="$RESULTS_BASE/$sample_name"
  mkdir -p "$outdir/logs"

  echo "Submitting: $sample_name"
  sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false -A "$ACCOUNT" \
    --job-name="onco_${sample_name}" \
    --output="${outdir}/logs/onco_${sample_name}_%j.out" \
    --error="${outdir}/logs/onco_${sample_name}_%j.err" \
    ${SOURCE_REPO}/run_oncoanalyser.sh \
    "$sample_sheet" \
    "$outdir" \
    targeted \
    GRCh38_hmf \
    "$ACCOUNT" \
    --ichor-dir $ICHOR_DIR \
    --tissue $TISSUE \
    --methylation-data /cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/EMSeq/MethylDackel/2025-10-10_MultiMMR_CpG_methylation_matrix.RData \
    --ihc-tsv /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v3/DNASeq/cohort_ihc.tsv \
    --panel-name "$PANEL_NAME" \
    --panel-config "$PANEL_CONFIG" \
    --amber-expanded-snps \
    --force-panel \
    --ichor-sub-clones FALSE \
    --ichor-exclusion-bed /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_exclusion.bed \
    --ichor-pon /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_panelOfNormals.rds

done

## Run once above is complete

python3 ${SOURCE_REPO}/scripts/compile_mmr_summary.py \
  --results-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results

## Package for sharing with package_mmr_reports.py
python3 /cluster/projects/pughlab/src/oncoanalyser_pughlab/scripts/package_mmr_reports.py \
  --results-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results \
  --output-zip /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/multimmr_cohort_reports.zip