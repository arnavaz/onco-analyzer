#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_swgs_ichorcna.sh
# Aligns shallow WGS (sWGS) FASTQs and runs ichorCNA copy number profiling.
#
# Usage:
#   run_swgs_ichorcna.sh --swgs-fastq-dir <dir> --sample-id <id> --outdir <dir> [options]
#
# Options:
#   --swgs-fastq-dir <dir>   Directory containing sWGS FASTQs (required)
#   --sample-id <id>         Sample identifier (e.g. T1-19-00888-TE_S1, required)
#   --patient-id <id>        Patient identifier (e.g. DYNA-P1, optional)
#   --outdir <dir>           Output directory for ichorCNA results (required)
#   --genome <genome>        Reference genome (GRCh38_hmf or GRCh37_hmf, default: GRCh38_hmf)
#   --threads <num>          Number of threads for alignment (default: 4)
#   --keep-bam               Keep intermediate sWGS BAM (default: false)
#   --source-repo <dir>      Path to oncoanalyser repository root (optional)
#   -h, --help               Show this help message
# ---------------------------------------------------------------------------

set -euo pipefail

SWGS_FASTQ_DIR=""
SAMPLE_ID=""
PATIENT_ID=""
OUTDIR=""
GENOME="GRCh38_hmf"
THREADS=4
KEEP_BAM="false"
SOURCE_REPO=""

usage() {
    cat >&2 <<EOF
Usage:
  $(basename "$0") --swgs-fastq-dir <dir> --sample-id <id> --outdir <dir> [options]

Options:
  --swgs-fastq-dir <dir>   Directory containing sWGS FASTQs
  --sample-id <id>         Sample identifier (e.g. T1-19-00888-TE_S1)
  --patient-id <id>        Patient identifier (e.g. DYNA-P1)
  --outdir <dir>           Output directory for ichorCNA results
  --genome <genome>        Reference genome (GRCh38_hmf or GRCh37_hmf, default: GRCh38_hmf)
  --threads <num>          Number of threads for alignment (default: 4)
  --keep-bam               Keep intermediate sWGS BAM (default: false)
  --source-repo <dir>      Path to oncoanalyser repository root
  -h, --help               Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --swgs-fastq-dir|--swgs_fastq_dir)
            SWGS_FASTQ_DIR="$2"
            shift 2
            ;;
        --sample-id|--sample_id)
            SAMPLE_ID="$2"
            shift 2
            ;;
        --patient-id|--patient_id)
            PATIENT_ID="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --genome)
            GENOME="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --keep-bam|--keep_bam)
            KEEP_BAM="true"
            shift
            ;;
        --source-repo|--source_repo)
            SOURCE_REPO="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${SWGS_FASTQ_DIR}" || -z "${SAMPLE_ID}" || -z "${OUTDIR}" ]]; then
    echo "ERROR: --swgs-fastq-dir, --sample-id, and --outdir are all required." >&2
    usage
    exit 1
fi

if [[ ! -d "${SWGS_FASTQ_DIR}" ]]; then
    echo "WARNING: sWGS FASTQ directory not found: ${SWGS_FASTQ_DIR}. Skipping sWGS ichorCNA." >&2
    exit 0
fi

# Resolve source repo if not provided
if [[ -z "${SOURCE_REPO}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    SOURCE_REPO="$(dirname "${SCRIPT_DIR}")"
fi

mkdir -p "${OUTDIR}"

# Check if final params.txt already exists
if [[ -f "${OUTDIR}/${SAMPLE_ID}.params.txt" && -s "${OUTDIR}/${SAMPLE_ID}.params.txt" ]]; then
    echo "sWGS ichorCNA results already exist for ${SAMPLE_ID} at ${OUTDIR}/${SAMPLE_ID}.params.txt. Skipping execution." >&2
    exit 0
fi

# Discover sWGS FASTQs
# Build candidate tokens to search for
CANDIDATE_TOKENS=()

# 1. Exact match on sample ID
CANDIDATE_TOKENS+=("${SAMPLE_ID}")

# 2. Extract standard cohort donor pattern (e.g. 19-00888 from T1-19-00888-TE_S1)
pat_token=$(echo "${SAMPLE_ID}" | grep -oE "[0-9]{2}-[0-9]{4,5}" | head -n 1 || true)
if [[ -n "${pat_token}" ]]; then
    CANDIDATE_TOKENS+=("${pat_token}")
fi

# 3. Patient ID tokens
if [[ -n "${PATIENT_ID}" ]]; then
    CANDIDATE_TOKENS+=("${PATIENT_ID}")
    p_num=$(echo "${PATIENT_ID}" | grep -oE "P[0-9]+" | head -n 1 || true)
    if [[ -n "${p_num}" ]]; then
        CANDIDATE_TOKENS+=("${p_num}")
    fi
fi

MATCHED_TOKEN=""
mapfile -t R1_FILES < <(
    for tok in "${CANDIDATE_TOKENS[@]}"; do
        files=($(find "${SWGS_FASTQ_DIR}" -maxdepth 2 -type f \( -name "*${tok}*R1*.fastq.gz" -o -name "*${tok}*R1*.fq.gz" \) 2>/dev/null | sort))
        if [[ ${#files[@]} -gt 0 ]]; then
            MATCHED_TOKEN="${tok}"
            printf "%s\n" "${files[@]}"
            break
        fi
    done
)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "WARNING: No matching sWGS FASTQs found for sample '${SAMPLE_ID}' (patient: '${PATIENT_ID}') in ${SWGS_FASTQ_DIR}. Skipping sWGS ichorCNA." >&2
    exit 0
fi

echo "Found ${#R1_FILES[@]} sWGS R1 FASTQ file(s) for ${SAMPLE_ID} matching token '${MATCHED_TOKEN}'." >&2

# Build paired R2 files
R2_FILES=()
for r1 in "${R1_FILES[@]}"; do
    r2="${r1/_R1/_R2}"
    if [[ ! -f "${r2}" ]]; then
        r2="${r1/_R1./_R2.}"
    fi
    if [[ ! -f "${r2}" ]]; then
        echo "ERROR: Missing paired R2 file for ${r1}" >&2
        exit 1
    fi
    R2_FILES+=("${r2}")
done

# Load environment modules on HPC
module load bwa-mem2/2.2.1 2>/dev/null || module load bwa-mem2 2>/dev/null || true
module load samtools/1.20 2>/dev/null || module load samtools/1.19 2>/dev/null || module load samtools 2>/dev/null || true
module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true
module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R/4.5.2 2>/dev/null || module load r/4.3.1 2>/dev/null || module load R 2>/dev/null || module load r 2>/dev/null || true

# Resolve BWA command
BWA_BIN="bwa-mem2"
if ! command -v bwa-mem2 >/dev/null 2>&1; then
    if command -v bwa >/dev/null 2>&1; then
        BWA_BIN="bwa"
    else
        # Try Singularity cache
        SING_CACHE="${SINGULARITY_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR:-/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache}}"
        bwa_img=$(find "${SING_CACHE}" -name "*bwa-mem2*.img" -o -name "*bwa-mem2*.sif" 2>/dev/null | head -n 1 || true)
        if [[ -n "${bwa_img}" && -f "${bwa_img}" ]]; then
            BWA_BIN="singularity exec ${bwa_img} bwa-mem2"
        else
            echo "ERROR: Neither bwa-mem2 nor bwa found in PATH or Singularity cache." >&2
            exit 1
        fi
    fi
fi

# Resolve reference data
REF_DATA_BASE="/cluster/projects/pughlab/references/oncoanalyser_resources"
if [[ "${GENOME}" == "GRCh37_hmf" ]]; then
    REF_FASTA="${REF_DATA_BASE}/oncoanalyser_refs/reference_data/2.3.0/GRCh37/Homo_sapiens.GRCh37.GATK.illumina.fasta"
    BWA_INDEX="${REF_FASTA}"
    GENOME_BUILD="hg19"
    CENT_FILE="GRCh37.p13_centromere_UCSC-gapTable.txt"
    DEFAULT_PON_FILE="HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
else
    REF_FASTA="${REF_DATA_BASE}/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta"
    BWA_INDEX="${REF_DATA_BASE}/oncoanalyser_refs/reference_data/2.3.0/bwa-mem2_index-2.2.1/GRCh38_masked_exclusions_alts_hlas.fasta"
    GENOME_BUILD="hg38"
    CENT_FILE="GRCh38.GCA_000001405.2_centromere_acen.txt"
    DEFAULT_PON_FILE="HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds"
fi

SWGS_BAM="${OUTDIR}/${SAMPLE_ID}.swgs.bam"

# Align if BAM does not exist
if [[ ! -f "${SWGS_BAM}" || ! -s "${SWGS_BAM}" ]]; then
    # Split available threads between bwa-mem2 and samtools sort, and bound sort memory
    SORT_THREADS=4
    if [[ "${THREADS}" -le 4 ]]; then
        SORT_THREADS=2
    fi
    ALIGN_THREADS=$(( THREADS - SORT_THREADS ))
    if [[ "${ALIGN_THREADS}" -lt 2 ]]; then
        ALIGN_THREADS=2
    fi

    echo "Aligning sWGS FASTQs for ${SAMPLE_ID} using ${BWA_BIN} (${ALIGN_THREADS} align threads, ${SORT_THREADS} sort threads)..." >&2
    ${BWA_BIN} mem -t "${ALIGN_THREADS}" "${BWA_INDEX}" \
        <(zcat "${R1_FILES[@]}") \
        <(zcat "${R2_FILES[@]}") \
        | samtools sort -@ "${SORT_THREADS}" -m 1G -o "${SWGS_BAM}" -
    samtools index "${SWGS_BAM}"
fi

# Determine chromosome naming style from BAM
first_chr=$(samtools idxstats "${SWGS_BAM}" 2>/dev/null | awk '{print $1}' | head -n 1 || true)
if [[ "${first_chr}" =~ ^chr ]]; then
    CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    GENOME_STYLE="UCSC"
    CHR_NORM="paste0('chr', 1:22)"
    CHR_TRAIN="paste0('chr', 1:22)"
    CHR_ANALYZE="c(paste0('chr', 1:22),'chrX')"
else
    CHRS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
    GENOME_STYLE="NCBI"
    CHR_NORM="c(1:22)"
    CHR_TRAIN="c(1:22)"
    CHR_ANALYZE="c(1:22,'X')"
fi

# Run readCounter
SWGS_WIG="${OUTDIR}/${SAMPLE_ID}.wig"
if [[ ! -f "${SWGS_WIG}" || ! -s "${SWGS_WIG}" ]]; then
    echo "Running readCounter for ${SAMPLE_ID}..." >&2
    readCounter --window 1000000 --quality 20 --chromosome "${CHRS}" "${SWGS_BAM}" > "${SWGS_WIG}"
fi

# Prepare ichorCNA reference files
DEFAULT_PON=$(Rscript -e "cat(system.file('extdata', '${DEFAULT_PON_FILE}', package='ichorCNA'))" 2>/dev/null || true)
GC_WIG=$(Rscript -e "cat(system.file('extdata', 'gc_${GENOME_BUILD}_1000kb.wig', package='ichorCNA'))")
MAP_WIG=$(Rscript -e "cat(system.file('extdata', 'map_${GENOME_BUILD}_1000kb.wig', package='ichorCNA'))")
CENTROMERE=$(Rscript -e "cat(system.file('extdata', '${CENT_FILE}', package='ichorCNA'))")

PON_ARG=""
if [[ -n "${DEFAULT_PON}" && -f "${DEFAULT_PON}" ]]; then
    PON_ARG="--normalPanel ${DEFAULT_PON}"
fi

echo "Running runIchorCNA.R for ${SAMPLE_ID}..." >&2
(
    cd "${OUTDIR}"
    Rscript "${SOURCE_REPO}/resources/runIchorCNA.R" \
        --WIG "${SWGS_WIG}" \
        --id "${SAMPLE_ID}" \
        --genomeStyle "${GENOME_STYLE}" \
        --genomeBuild "${GENOME_BUILD}" \
        --gcWig "${GC_WIG}" \
        --mapWig "${MAP_WIG}" \
        --centromere "${CENTROMERE}" \
        ${PON_ARG} \
        --chrNormalize "${CHR_NORM}" \
        --chrTrain "${CHR_TRAIN}" \
        --chrs "${CHR_ANALYZE}" \
        --normal "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.0)" \
        --ploidy "c(2,3)" \
        --estimateScPrevalence TRUE \
        --scStates "c(1,3)"
)

# Extract read count and write QC stats
total_reads=$(samtools idxstats "${SWGS_BAM}" 2>/dev/null | awk '{sum+=$3} END {print sum}' || echo "0")
cov=0
if [[ "${total_reads}" =~ ^[0-9]+$ ]] && [[ "${total_reads}" -gt 0 ]]; then
    cov=$(awk -v reads="${total_reads}" 'BEGIN {printf "%.4f", (reads * 150) / 3100000000}')
fi

cat <<EOF > "${OUTDIR}/qc_stats.json"
{
  "source": "swgs",
  "total_reads": ${total_reads},
  "coverage": ${cov},
  "genome_coverage": ${cov},
  "off_target_coverage": ${cov}
}
EOF

# Clean up WIG file
rm -f "${SWGS_WIG}"

# Remove BAM unless requested to keep
if [[ "${KEEP_BAM}" != "true" ]]; then
    rm -f "${SWGS_BAM}" "${SWGS_BAM}.bai"
fi

echo "Completed sWGS ichorCNA for ${SAMPLE_ID} under ${OUTDIR}" >&2
exit 0
