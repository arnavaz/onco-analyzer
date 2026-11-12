#!/bin/bash
# ==============================================================================
# Test ichorCNA with FULL Complement of Short Reads (90-150 bp) against PoN
# ------------------------------------------------------------------------------
# Context:
#   In BFIERCE_0003 (sample 536878), the short fragment pool contains ~234M reads
#   (~9.07X coverage). Previous runs downsampled this pool to 1.0X (~25M reads),
#   throwing away 89% of short reads and drowning low ctDNA signal in noise.
#
#   This script evaluates ichorCNA on 100% of the short fragments (no downsampling)
#   against a Panel of Normals (PoN), avoiding plasma vs buffy-coat artifacts.
#
# Usage:
#   # 1. Run with standard Broad Institute hg38 1Mb PoN:
#   sbatch -A pughlab -p pughlab runs/test_BFIERCE_0003_full_short_ichor.sh
#
#   # 2. Run with in-house healthy blood control PoN:
#   sbatch -A pughlab -p pughlab runs/test_BFIERCE_0003_full_short_ichor.sh /path/to/inhouse_blood_pon.rds
#
#   Note: The generated WIG file (100% short reads) is cached. Re-running against
#   a new PoN takes ~60 seconds because read extraction is not repeated.
# ==============================================================================
#SBATCH --job-name=ichor_full_short_0003
#SBATCH --partition=pughlab
#SBATCH --account=pughlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3/BFIERCE_0003/test_full_short_pon/logs/test_%j.out
#SBATCH --error=/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3/BFIERCE_0003/test_full_short_pon/logs/test_%j.err

set -euo pipefail

# ------------------------------------------------------------------------------
# Modules & Environment
# ------------------------------------------------------------------------------
module load samtools/1.20
module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true
module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R 2>/dev/null || true

PATIENT="BFIERCE_0003"
SAMPLE="536878"
BAM_PATH="/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser/BFIERCE_0003/BFIERCE-FU-004/alignments/dna/536878.redux.bam"
REF_FASTA="/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta"
SOURCE_REPO="${SOURCE_REPO:-/cluster/projects/pughlab/src/oncoanalyser_pughlab}"
OUT_BASE="/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3/${PATIENT}/test_full_short_pon"
TMPDIR="${OUT_BASE}/tmp"

mkdir -p "${OUT_BASE}/logs" "${TMPDIR}"
export TMPDIR="${TMPDIR}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

# ==============================================================================
# Resolve Panel of Normals (PoN)
# ------------------------------------------------------------------------------
CHARM_PON="/cluster/projects/pughlab/references/panel_of_normals/hg38/IchorCNA/pon_ichorCNA__CHARM__median.rds"
BROAD_PON=$(Rscript -e "cat(system.file('extdata', 'HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds', package='ichorCNA'))" 2>/dev/null || true)

# Default to in-house CHARM PoN unless overridden by $1
PON_FILE="${1:-${CHARM_PON}}"

if [[ ! -f "${PON_FILE}" ]]; then
    echo "ERROR: PoN file not found: ${PON_FILE}" >&2
    exit 1
fi

PON_LABEL=$(basename "${PON_FILE}" .rds)
echo "=================================================================="
echo " Starting Full-Complement Short Read ichorCNA Test"
echo " Patient:      ${PATIENT}"
echo " Sample:       ${SAMPLE}"
echo " Input BAM:    ${BAM_PATH}"
echo " Primary PoN:  ${PON_FILE} (${PON_LABEL})"
echo " Output Dir:   ${OUT_BASE}/run_${PON_LABEL}"
echo " Date:         $(date)"
echo " Host:         $(hostname)"
echo "=================================================================="

# ------------------------------------------------------------------------------
# Step 1: Generate Full-Complement Short Fragment WIG (Cached)
# ------------------------------------------------------------------------------
WIG_FILE="${OUT_BASE}/${SAMPLE}.short_full_1Mb.wig"
CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

if [[ -f "${WIG_FILE}" && -s "${WIG_FILE}" ]]; then
    echo "[$(date)] Found existing full-short WIG: ${WIG_FILE}. Skipping BAM filtering."
else
    TEMP_SHORT_BAM="${TMPDIR}/${SAMPLE}.short_full_$$.bam"
    echo "[$(date)] Slicing 100% of short fragments (90-150 bp) from ${BAM_PATH}..."
    echo "         (Zero downsampling; retaining all ~234M short reads)"
    
    samtools view -@ 8 -T "${REF_FASTA}" -b \
        -e "(tlen >= 90 && tlen <= 150) || (tlen <= -90 && tlen >= -150)" \
        "${BAM_PATH}" > "${TEMP_SHORT_BAM}"

    echo "[$(date)] Indexing short fragment BAM..."
    samtools index -@ 8 "${TEMP_SHORT_BAM}"

    SHORT_READ_COUNT=$(samtools idxstats "${TEMP_SHORT_BAM}" | awk '{s+=$3} END {print s}')
    echo "[$(date)] Total short reads retained: ${SHORT_READ_COUNT}"

    echo "[$(date)] Running readCounter (1 Mb bins)..."
    readCounter --window 1000000 --quality 20 --chromosome "${CHRS}" "${TEMP_SHORT_BAM}" > "${WIG_FILE}"

    echo "[$(date)] WIG generation complete. Removing temporary BAM..."
    rm -f "${TEMP_SHORT_BAM}" "${TEMP_SHORT_BAM}.bai"
fi

# ------------------------------------------------------------------------------
# Step 2: Execute ichorCNA against PoN
# ------------------------------------------------------------------------------
GC_WIG=$(Rscript -e "cat(system.file('extdata', 'gc_hg38_1000kb.wig', package='ichorCNA'))" 2>/dev/null || true)
MAP_WIG=$(Rscript -e "cat(system.file('extdata', 'map_hg38_1000kb.wig', package='ichorCNA'))" 2>/dev/null || true)
CENTROMERE=$(Rscript -e "cat(system.file('extdata', 'GRCh38.GCA_000001405.2_centromere_acen.txt', package='ichorCNA'))" 2>/dev/null || true)

run_ichor_pon() {
    local target_pon="$1"
    local label=$(basename "${target_pon}" .rds)
    local run_dir="${OUT_BASE}/run_${label}"
    mkdir -p "${run_dir}"

    echo ""
    echo "------------------------------------------------------------------"
    echo "[$(date)] Executing ichorCNA with PoN: ${label}"
    echo " Directory: ${run_dir}"
    echo "------------------------------------------------------------------"

    (
        cd "${run_dir}"
        Rscript "${SOURCE_REPO}/resources/runIchorCNA.R" \
            --WIG "${WIG_FILE}" \
            --id "${SAMPLE}_short_full" \
            --genomeStyle "UCSC" \
            --genomeBuild "hg38" \
            --gcWig "${GC_WIG}" \
            --mapWig "${MAP_WIG}" \
            --centromere "${CENTROMERE}" \
            --normalPanel "${target_pon}" \
            --chrNormalize "paste0('chr', 1:22)" \
            --chrTrain "paste0('chr', 1:22)" \
            --chrs "c(paste0('chr', 1:22),'chrX')" \
            --normal "c(0.5,0.6,0.7,0.8,0.9)" \
            --ploidy "c(2,3)" \
            --scStates "NULL"
    )

    # Flatten nested directory if created by plotSolutions
    if [[ -d "${run_dir}/${SAMPLE}_short_full" ]]; then
        mv -f "${run_dir}/${SAMPLE}_short_full"/* "${run_dir}/" 2>/dev/null || true
        rmdir "${run_dir}/${SAMPLE}_short_full" 2>/dev/null || true
    fi

    local params_file="${run_dir}/${SAMPLE}_short_full.params.txt"
    if [[ -f "${params_file}" ]]; then
        echo " [Done] Metrics for ${label}:"
        head -n 12 "${params_file}"
    else
        echo " WARNING: Params file not found at ${params_file}"
    fi
}

# Run primary target PoN (CHARM PoN by default)
run_ichor_pon "${PON_FILE}"

# If primary was CHARM and Broad PoN is available, run Broad as well for side-by-side comparison
if [[ "${PON_FILE}" == "${CHARM_PON}" && -n "${BROAD_PON}" && -f "${BROAD_PON}" ]]; then
    run_ichor_pon "${BROAD_PON}"
fi

echo ""
echo "=================================================================="
echo " [$(date)] All ichorCNA Comparisons Complete!"
echo " Results stored in: ${OUT_BASE}"
echo "=================================================================="


