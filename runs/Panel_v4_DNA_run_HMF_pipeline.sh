
## taken from /MultiMMR_HMF/panel/v4/create_panel_resources.sh
## create bed file


### WHEN RE-BUILDING RESOURCES, CHECK QC first (don't use N5 from the first batch, coverage is too low)

module load tabix/0.2.6

cd /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/

sort -k1,1 -k2,2n Probes_merged_ok_shareable_UHN_MultiMMRv4_TE-95402177_hg38_MergedDesign_260128024205.bed \
| bgzip -c > MultiMMRv4_target_regions_definition.38.bed.gz

tabix -p bed MultiMMRv4_target_regions_definition.38.bed.gz

## 1. Create target link directories and symlink NL82 Normal and Tumour samples
mkdir -p /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/NL82_Normal_cram_links
mkdir -p /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/NL82_Tumour_cram_links
SRC_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/NL82_v3/DNA/results"
NORMAL_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/NL82_Normal_cram_links"
TUMOUR_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/NL82_Tumour_cram_links"

find "$SRC_DIR"/*/*/alignments/dna -type f \( -name "*.redux.bam" -o -name "*.redux.bam.bai" -o -name "*.redux.cram" -o -name "*.redux.cram.crai" \) 2>/dev/null | \
grep -E "(B[0-9]{2}-W|-[Ww]-|-[Ww]$)" | while read -r file; do
    ln -sf "$file" "$NORMAL_DIR/"
done

find "$SRC_DIR"/*/*/alignments/dna -type f \( -name "*.redux.bam" -o -name "*.redux.bam.bai" -o -name "*.redux.cram" -o -name "*.redux.cram.crai" \) 2>/dev/null | \
grep -E "(T[0-9]{2}|-[Ff]-|-[Ff]$|-[Rr]-|-[Rr]$|UWID|GI019|EN0100)" | while read -r file; do
    ln -sf "$file" "$TUMOUR_DIR/"
done


## 2. Build combined samplesheet for panel_resource_creation (unaligned FASTQ normals + pre-aligned BAM/CRAM normals)
(
  echo "group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath"
  
  # 1. Unaligned Dynacare Blood Normal FASTQs
  for fq1 in /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/fastq_files/devRun1MMR/dnaSeq/*_R1*.fastq.gz; do
    [ -e "$fq1" ] || continue
    fq2="${fq1/_R1/_R2}"
    sample=$(basename "$fq1" | sed 's/_.*//')
    sample_type=$(echo $sample | egrep -o '^.' | sed 's/^T/tumor/' | sed 's/^N/normal/')
    
    if [ "$sample_type" = "normal" ]; then
      subject_id=$(echo $sample | sed 's/-TE_.*//' | sed 's/[T|N]/P/')
      lane=$(echo "$fq1" | grep -oE '_L[0-9]+_' | tr -d '_L')
      if [ -z "$lane" ]; then lane="001"; fi
      echo "$subject_id,$subject_id,${sample},normal,dna,fastq,library_id:${sample};lane:${lane},$fq1;$fq2"
    fi
  done

  # 2. Pre-aligned NL82 Blood Normal BAMs/CRAMs
  for f in /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/NL82_Normal_cram_links/*.bam; do
    [ -e "$f" ] || continue
    fname=$(basename "$f" | sed -E 's/^NL82-//; s/\.redux\..*//')
    echo "$fname,$fname,$fname,normal,dna,bam,,$f"
  done
) > /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/samplesheet.panel_resource_creation.csv


## Run resource creation workflow for panel

module load java/21
module load singularity/4.3.2
export PATH="$HOME/.local/bin:$PATH"
export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

cd /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/

sbatch -p pughlab -A pughlab /cluster/projects/pughlab/src/oncoanalyser_multimmr/run_oncoanalyser.sh \
/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/samplesheet.panel_resource_creation.csv \
/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/panel_resource \
panel_resource_creation \
GRCh38_hmf \
pughlab \
--driver-gene-panel /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/MMR_drivers.tsv \
--target-regions-bed /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/MultiMMRv4_target_regions_definition.38.bed.gz \
-- -resume --process.queue pughlab

## Generate SNP panel augmentation


#run locally (because internet)
cd '/Users/jbruce/bin/oncoanalyser_multimmr/'

## Instead of   --default-loci /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/hmf_pipeline_resources.38_v2.3.0--2/dna/copy_number/AmberGermlineSites.38.tsv.gz \
## using previous version to speed things along.  For new panels use the hmf one

# python3 scripts/generate_panel_snps.py \
#   --target-regions-bed /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/MultiMMRv4_target_regions_definition.38.bed.gz \
#   --default-loci /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v3/AmberGermlineSites.expanded.38.tsv.gz \
#   --output /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/AmberGermlineSites.expanded.38.tsv.gz \
#   --genes all_bed \
#   --genome GRCh38_hmf

# ## Regenerate with all reported SNPs
# mv /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/AmberGermlineSites.expanded.38.tsv.gz /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/MAF_01_AmberGermlineSites.expanded.38.tsv.gz

python3 scripts/generate_panel_snps.py \
  --target-regions-bed /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/MultiMMRv4_target_regions_definition.38.bed.gz \
  --default-loci /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/hmf_pipeline_resources.38_v2.3.0--2/dna/copy_number/AmberGermlineSites.38.tsv.gz \
  --output /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/AmberGermlineSites.expanded.38.tsv.gz \
  --genes all_bed \
  --maf 0.0 \
  --no-cache \
  --genome GRCh38_hmf


## check how many SNPs were added

zgrep -v "^Chromosome" /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/AmberGermlineSites.expanded.38.tsv.gz | wc -l
zgrep -v "^Chromosome" /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/hmf_pipeline_resources.38_v2.3.0--2/dna/copy_number/AmberGermlineSites.38.tsv.gz | wc -l


## ==============================================================================
## Generate off-target ichorCNA Panel of Normals (PoN) and Empirical Exclusion BED
## ==============================================================================
# Run on an interactive compute node allocation (or sbatch):
# salloc -c 10 -p superhimem -t 10:0:0 --mem 120G
# conda activate nf-core
# module load hmmcopy_utils/170718
# module load R/4.1.0

ICHOR_RES_DIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/ichor_resources"
mkdir -p "${ICHOR_RES_DIR}"

python3 /cluster/projects/pughlab/src/oncoanalyser_multimmr/scripts/generate_ichorcna_panel_resources.py \
  --bams /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/NL82_Normal_cram_links/*.bam \
  --chrom-sizes /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/hg38_genome.chrom.sizes \
  --output-dir "${ICHOR_RES_DIR}" \
  --mosdepth-bin /cluster/projects/pughlab/bin/mosdepth/mosdepth \
  --padding 250 \
  --threads 10


## match tumour normal and generate per-patient samplesheets from Dyna_v4 FASTQs

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_multimmr

PANEL_DIR=${SOURCE_REPO}/resources/panels/multimmr_v4
PANEL_NAME=multimmr_hmf
PANEL_CONFIG=$PANEL_DIR/panel_data.config
OUTDIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/DNASeq/results
DNA_FASTQ_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/fastq_files/devRun1MMR/dnaSeq
SWGS_FASTQ_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/fastq_files/devRun1MMR/sWGS
EMSEQ_FASTQ_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/fastq_files/devRun1MMR/emSeq
ACCOUNT=deperrotgroup

mkdir -p "$OUTDIR"

# Generate per-patient samplesheets (DYNA-P1 through DYNA-P5)
python3 ${SOURCE_REPO}/scripts/make_dyna_samplesheets.py \
   --fastq-dir "${DNA_FASTQ_DIR}" \
   --out "${OUTDIR}" \
   --patient-prefix "DYNA-P"

## Submit patients with updated germline modifications and correct panel files

SHEETS_DIR="${OUTDIR}"
RESULTS_BASE="${OUTDIR}"

SKIP_SAMPLE=NONE

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
  sample_name=$(basename "$(dirname "$sample_sheet")")

  if [[ "$sample_name" == "$SKIP_SAMPLE" ]]; then
    echo "Skipping already-submitted sample: $sample_name"
    continue
  fi

  ICHOR_OPTS=()
  if [[ -n "${SWGS_FASTQ_DIR}" && -d "${SWGS_FASTQ_DIR}" ]]; then
      echo "Using sWGS FASTQ directory for on-the-fly ichorCNA: ${SWGS_FASTQ_DIR}"
      ICHOR_OPTS=(--swgs-fastq-dir "${SWGS_FASTQ_DIR}")
  else
      echo "WARNING: No sWGS FASTQ directory found. Running without sWGS fallback."
  fi

  EMSEQ_OPTS=()
  if [[ -n "${EMSEQ_FASTQ_DIR}" && -d "${EMSEQ_FASTQ_DIR}" ]]; then
      echo "Using EM-Seq FASTQ directory for methylation analysis: ${EMSEQ_FASTQ_DIR}"
      EMSEQ_OPTS=(--emseq-fastq-dir "${EMSEQ_FASTQ_DIR}")
  fi

  TISSUE="colorectal"

  outdir="$RESULTS_BASE/$sample_name"
  mkdir -p "$outdir/logs"

  METH_DATA="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/EMSeq/MethylDackel/MultiMMR_CpG_methylation_matrix.RData"
  if [[ ! -f "${METH_DATA}" ]]; then
      METH_DATA="${outdir}/EMSeq/MethylDackel/MultiMMR_CpG_methylation_matrix.RData"
  fi
  METH_OPTS=()
  if [[ -f "${METH_DATA}" ]]; then
      METH_OPTS=(--methylation-data "${METH_DATA}")
  fi

  ORTHOGONAL_TSV="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/metadata/orthogonal_test_mutations.tsv"
  ORTH_OPTS=()
  if [[ -f "${ORTHOGONAL_TSV}" ]]; then
      ORTH_OPTS=(--orthogonal-tsv "${ORTHOGONAL_TSV}")
  fi

  PANEL_ICHOR_OPTS=()
  EXCLUSION_BED="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/ichor_resources/custom_exclusion.bed"
  PON_RDS="/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v4/ichor_resources/custom_panelOfNormals.rds"
  if [[ -f "${EXCLUSION_BED}" && -f "${PON_RDS}" ]]; then
      PANEL_ICHOR_OPTS=(
          --ichor-exclusion-bed "${EXCLUSION_BED}"
          --ichor-pon "${PON_RDS}"
          --ichor-sub-clones FALSE
      )
  fi

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
    "${ICHOR_OPTS[@]}" \
    "${PANEL_ICHOR_OPTS[@]}" \
    "${EMSEQ_OPTS[@]}" \
    --tissue $TISSUE \
    --sequencing-platform illumina \
    "${METH_OPTS[@]}" \
    "${ORTH_OPTS[@]}" \
    --panel-name "$PANEL_NAME" \
    --panel-config "$PANEL_CONFIG" \
    --amber-expanded-snps \
    --force-panel \
    --force-nextflow

done

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_multimmr

python3 ${SOURCE_REPO}/scripts/compile_mmr_summary.py \
  --results-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v4/DNASeq/results
