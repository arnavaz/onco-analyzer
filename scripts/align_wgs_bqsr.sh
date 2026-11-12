#!/bin/bash
# =============================================================================
# align_wgs_bqsr.sh
#
# Mimics the OICR GSI bamMergePreprocessing pipeline:
#   BWA-MEM (per lane, parallel SLURM jobs)
#   -> Picard MarkDuplicates (OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 for NovaSeq X)
#   -> GATK BaseRecalibrator + ApplyBQSR
#
# Read group tags (@RG) are auto-parsed from Illumina FASTQ read names.
# The PU field is reconstructed as: {date}_{instrument}_{run}_{flowcell}-{i7}-{i5}_{lane}
#
# Usage:
#   bash align_wgs_bqsr.sh \
#       -s SAMPLE_ID       (SM tag; e.g. BDWGTS_0388_Ov_P_A15)
#       -l LIBRARY_ID      (LB tag; e.g. BDWGTS_0388_02_LB01-01)
#       -f FASTQ_DIR       (directory containing per-lane R1/R2 FASTQ.gz pairs)
#       -o OUTPUT_DIR
#       [-a SLURM_ACCOUNT] (default: deperrotgroup)
#       [-t TMP_DIR]       (default: $SLURM_TMPDIR or /tmp)
#       [--no-bqsr]        (stop after MarkDuplicates; skip BQSR)
#       [--dry-run]        (print commands without submitting)
#
# FASTQ pairing: finds all *_R1*.fastq.gz files and pairs each with the
# corresponding *_R2*.fastq.gz (same filename with R1 -> R2 substitution).
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Reference paths (H4H defaults — adjust if needed)
# ---------------------------------------------------------------------------
REFERENCE="/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa"

# BQSR known sites — GATK hg38 bundle
BUNDLE="/cluster/tools/data/genomes/human/hg38/hg38bundle"
KNOWN_SNPS="${BUNDLE}/Homo_sapiens_assembly38.dbsnp.vcf.gz"
KNOWN_MILLS="${BUNDLE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_INDELS="${BUNDLE}/Homo_sapiens_assembly38.known_indels.vcf.gz"
KNOWN_1KG="${BUNDLE}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# ---------------------------------------------------------------------------
# Modules (H4H available versions)
# ---------------------------------------------------------------------------
MODULE_BWA="bwa/0.7.15"         # same generation as OICR pipeline
MODULE_SAMTOOLS="samtools/1.20"
MODULE_GATK="gatk/4.1.8.1"      # closest to original 4.1.6.0; avoid 4.6+ (BQSR differences)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SAMPLE_ID=""
LIBRARY_ID=""
FASTQ_DIR=""
OUTPUT_DIR=""
SLURM_ACCOUNT="deperrotgroup"
TMP_DIR=""
SLURM_PARTITION=""
NO_BQSR="false"
DRY_RUN="false"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--sample)   SAMPLE_ID="$2";    shift 2 ;;
        -l|--library)  LIBRARY_ID="$2";   shift 2 ;;
        -f|--fastqs)   FASTQ_DIR="$2";    shift 2 ;;
        -o|--outdir)   OUTPUT_DIR="$2";   shift 2 ;;
        -a|--account)  SLURM_ACCOUNT="$2"; shift 2 ;;
        -p|--partition) SLURM_PARTITION="$2"; shift 2 ;;
        -t|--tmpdir)   TMP_DIR="$2";      shift 2 ;;
        --no-bqsr)     NO_BQSR="true";    shift ;;
        --dry-run)     DRY_RUN="true";    shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

for required in SAMPLE_ID LIBRARY_ID FASTQ_DIR OUTPUT_DIR; do
    [[ -z "${!required}" ]] && { echo "ERROR: --${required,,} is required."; exit 1; }
done

[[ ! -d "$FASTQ_DIR" ]] && { echo "ERROR: FASTQ dir not found: $FASTQ_DIR"; exit 1; }

if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${SLURM_TMPDIR:-${OUTPUT_DIR}/tmp}"
fi

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------
ALIGN_DIR="${OUTPUT_DIR}/per_lane_bams"
LOGS_DIR="${OUTPUT_DIR}/logs"
FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.filter.deduped.recalibrated.bam"
DEDUPED_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.filter.deduped.bam"
METRICS_FILE="${OUTPUT_DIR}/${SAMPLE_ID}.markdup.metrics"
RECAL_TABLE="${OUTPUT_DIR}/${SAMPLE_ID}.bqsr.table"

mkdir -p "$ALIGN_DIR" "$LOGS_DIR" "$TMP_DIR"

# SLURM Partition lines
PARTITION_LINE_ALIGN=""
PARTITION_LINE_MD="#SBATCH --partition=himem"
PARTITION_LINE_BQSR=""

if [[ -n "$SLURM_PARTITION" ]]; then
    PARTITION_LINE_ALIGN="#SBATCH --partition=${SLURM_PARTITION}"
    PARTITION_LINE_MD="#SBATCH --partition=${SLURM_PARTITION}"
    PARTITION_LINE_BQSR="#SBATCH --partition=${SLURM_PARTITION}"
fi

# ---------------------------------------------------------------------------
# Pre-flight check: exit early if the final target BAM already exists and is valid
# ---------------------------------------------------------------------------
TARGET_BAM="$FINAL_BAM"
[[ "$NO_BQSR" == "true" ]] && TARGET_BAM="$DEDUPED_BAM"

if [[ -f "$TARGET_BAM" && -s "$TARGET_BAM" && -f "${TARGET_BAM%.bam}.bai" ]]; then
    echo "============================================================"
    echo "SUCCESS: Target BAM already exists: ${TARGET_BAM}"
    echo "Nothing to do. Graceful resume skipped submission."
    echo "============================================================"
    exit 0
fi

echo "============================================================"
echo "WGS Alignment + BQSR Pipeline"
echo "  Sample:      ${SAMPLE_ID}"
echo "  Library:     ${LIBRARY_ID}"
echo "  FASTQs:      ${FASTQ_DIR}"
echo "  Output:      ${OUTPUT_DIR}"
echo "  Reference:   ${REFERENCE}"
echo "  SLURM acct:  ${SLURM_ACCOUNT}"
echo "  Dry run:     ${DRY_RUN}"
echo "============================================================"

# ---------------------------------------------------------------------------
# Helper: parse @RG fields from an Illumina FASTQ read name.
#
# Illumina format: @instrument:run:flowcell:lane:tile:x:y read:filtered:ctrl:index
# Constructs PU/ID as: {date}_{instrument}_{run_padded}_{flowcell}-{i7}-{i5}_{lane}
# ---------------------------------------------------------------------------
parse_rg_from_fastq() {
    local r1="$1"
    local sample="$2"
    local library="$3"

    local header
    header=$(zcat "$r1" | head -1)

    # Field 1 (before space): @instrument:run:flowcell:lane:...
    local field1
    field1=$(echo "$header" | awk '{print $1}')
    local instrument run flowcell lane
    instrument=$(echo "$field1" | cut -d: -f1 | sed 's/^@//')
    run=$(echo "$field1" | cut -d: -f2)
    flowcell=$(echo "$field1" | cut -d: -f3)
    lane=$(echo "$field1" | cut -d: -f4)

    # Zero-pad run to 4 digits
    local run_padded
    run_padded=$(printf "%04d" "$run")

    # Field 2 (after space): read:filtered:ctrl:index (e.g. 1:N:0:TCCGTATG+TGGATGGT)
    local index i7 i5
    index=$(echo "$header" | awk '{print $2}' | cut -d: -f4)
    i7=$(echo "$index" | cut -d+ -f1)
    i5=$(echo "$index" | cut -d+ -f2)
    # Handle single-index (no +): i5 will be empty
    if [[ -z "$i5" ]]; then
        index_tag="${i7}"
    else
        index_tag="${i7}-${i5}"
    fi

    local filename
    filename=$(basename "$r1")
    local date_tag="00000000"
    if [[ "$filename" =~ _([0-9]{8})_ ]]; then
        date_tag="${BASH_REMATCH[1]}"
    fi

    local rg_id="${date_tag}_${instrument}_${run_padded}_${flowcell}-${index_tag}_${lane}"
    local rg_pu="${date_tag}_${instrument}_${run_padded}_${flowcell}-${index_tag}_${lane}"

    # Return as tab-separated string: ID PU SM LB PL lane
    printf "%s\t%s\t%s\t%s\tILLUMINA\t%s" \
        "$rg_id" "$rg_pu" "$sample" "$library" "$lane"
}

# ---------------------------------------------------------------------------
# Step 1: Submit per-lane BWA-MEM alignment jobs
# ---------------------------------------------------------------------------
echo ""
echo "Step 1: Discovering FASTQ pairs and submitting alignment jobs..."

ALIGN_JOB_IDS=()
LANE_BAMS=()

# Determine FASTQ pattern: filter by LIBRARY_ID if files containing it exist,
# otherwise fall back to all R1 files in the directory.
R1_FILES=()
for pattern in "${FASTQ_DIR}"/*"${LIBRARY_ID}"*_R1*.fastq.gz "${FASTQ_DIR}"/*"${LIBRARY_ID}"*_R1*.fq.gz; do
    [[ -f "$pattern" ]] && R1_FILES+=("$pattern")
done
if [ ${#R1_FILES[@]} -eq 0 ]; then
    for pattern in "${FASTQ_DIR}"/*_R1*.fastq.gz "${FASTQ_DIR}"/*_R1*.fq.gz; do
        [[ -f "$pattern" ]] && R1_FILES+=("$pattern")
    done
fi

for R1 in ${R1_FILES[@]+"${R1_FILES[@]}"}; do

    # Find matching R2
    R2="${R1/_R1/_R2}"
    if [[ ! -f "$R2" ]]; then
        echo "  WARNING: No R2 found for $R1 — skipping"
        continue
    fi

    # Parse RG fields from FASTQ header
    RG_INFO=$(parse_rg_from_fastq "$R1" "$SAMPLE_ID" "$LIBRARY_ID")
    RG_ID=$(echo "$RG_INFO" | cut -f1)
    RG_PU=$(echo "$RG_INFO" | cut -f2)
    LANE=$(echo "$RG_INFO" | cut -f6)

    LANE_BAM="${ALIGN_DIR}/${RG_ID}.bam"
    LANE_BAMS+=("$LANE_BAM")

    echo "  Lane: ${LANE} | RG_ID: ${RG_ID}"
    echo "    R1: $(basename "$R1")"
    echo "    R2: $(basename "$R2")"

    # Extract flowcell tag to build unique job name
    FLOWCELL_TAG=$(echo "$RG_ID" | cut -d_ -f4 | cut -d- -f1)
    JOB_NAME="align_${SAMPLE_ID}_${FLOWCELL_TAG}_L${LANE}"

    # Check if this lane alignment has already finished successfully
    if [[ -f "$LANE_BAM" && -s "$LANE_BAM" && -f "${LANE_BAM}.bai" ]]; then
        echo "    [SKIP] Lane ${LANE} alignment already complete: $(basename "$LANE_BAM")"
        continue
    fi

    # Check if there is already a running or queued job for this lane
    existing_job_id=""
    if command -v squeue >/dev/null 2>&1; then
        existing_job_id=$(squeue -u "$USER" -n "$JOB_NAME" -h -o "%A" | head -n 1 || true)
    fi

    if [[ -n "$existing_job_id" ]]; then
        echo "    [SKIP] Job ${JOB_NAME} is already queued or running in SLURM: ${existing_job_id}"
        ALIGN_JOB_IDS+=("$existing_job_id")
        continue
    fi

    # Build @RG string for bwa mem (using printf to generate literal tabs)
    RG_TAG=$(printf "@RG\tID:%s\tSM:%s\tLB:%s\tPL:ILLUMINA\tPU:%s" "$RG_ID" "$SAMPLE_ID" "$LIBRARY_ID" "$RG_PU")

    # Inline alignment job script
    ALIGN_SCRIPT=$(cat <<SCRIPT
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --account=${SLURM_ACCOUNT}
${PARTITION_LINE_ALIGN}
#SBATCH --output=${LOGS_DIR}/align_${FLOWCELL_TAG}_L${LANE}.%j.out
#SBATCH --error=${LOGS_DIR}/align_${FLOWCELL_TAG}_L${LANE}.%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G

module load ${MODULE_BWA}
module load ${MODULE_SAMTOOLS}

set -euo pipefail

# Use local node scratch space if running in SLURM, otherwise fall back to TMP_DIR
TMP_DIR="\${SLURM_TMPDIR:-${TMP_DIR}}"

echo "Aligning lane ${LANE}: \$(date)"

bwa mem -t 16 -Y -R "${RG_TAG}" "${REFERENCE}" "${R1}" "${R2}" | samtools sort -@ 4 -m 4G -T "\${TMP_DIR}/${RG_ID}" -o "${LANE_BAM}"

samtools index "${LANE_BAM}"
echo "Lane ${LANE} done: \$(date)"
SCRIPT
)

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] Would submit: ${JOB_NAME}"
        ALIGN_JOB_IDS+=("DRY_RUN_JOB")
    else
        JOB_ID=$(echo "$ALIGN_SCRIPT" | sbatch --parsable)
        ALIGN_JOB_IDS+=("$JOB_ID")
        echo "  Submitted alignment job: ${JOB_ID}"
    fi
done

if [[ ${#LANE_BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No R1 FASTQ files found in ${FASTQ_DIR}"
    exit 1
fi

echo ""
echo "  ${#LANE_BAMS[@]} lanes discovered."

# Build --INPUT flags for MarkDuplicates
MD_INPUTS=""
for bam in ${LANE_BAMS[@]+"${LANE_BAMS[@]}"}; do
    MD_INPUTS="${MD_INPUTS} --INPUT ${bam}"
done

# Build SLURM dependency string for active alignment jobs
# Filters out any blank entries (from skipped lanes)
ALIGN_DEPS=()
for job in ${ALIGN_JOB_IDS[@]+"${ALIGN_JOB_IDS[@]}"}; do
    [[ -n "$job" ]] && ALIGN_DEPS+=("$job")
done
DEPENDENCY_STR=$(IFS=:; echo "${ALIGN_DEPS[*]}")

# ---------------------------------------------------------------------------
# Step 2: MarkDuplicates
# ---------------------------------------------------------------------------
echo ""

# Check if MarkDuplicates output already exists and is valid
SKIP_MARKDUP="false"
if [[ -f "$DEDUPED_BAM" && -s "$DEDUPED_BAM" && -f "${DEDUPED_BAM%.bam}.bai" ]]; then
    echo "Step 2: [SKIP] MarkDuplicates output already exists: ${DEDUPED_BAM}"
    SKIP_MARKDUP="true"
fi

if [[ "$SKIP_MARKDUP" == "true" ]]; then
    MARKDUP_JOB_ID=""
else
    echo "Step 2: Submitting MarkDuplicates job..."

    MARKDUP_SCRIPT=$(cat <<SCRIPT
#!/bin/bash
#SBATCH --job-name=markdup_${SAMPLE_ID}
#SBATCH --account=${SLURM_ACCOUNT}
${PARTITION_LINE_MD}
#SBATCH --output=${LOGS_DIR}/markdup.%j.out
#SBATCH --error=${LOGS_DIR}/markdup.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G

module load ${MODULE_GATK}

set -euo pipefail

# Use local node scratch space if running in SLURM, otherwise fall back to TMP_DIR
TMP_DIR="\${SLURM_TMPDIR:-${TMP_DIR}}"

echo "MarkDuplicates: \$(date)"

gatk --java-options "-Xmx108G" MarkDuplicates ${MD_INPUTS} --OUTPUT "${DEDUPED_BAM}" --METRICS_FILE "${METRICS_FILE}" --REMOVE_DUPLICATES false --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --SORTING_COLLECTION_SIZE_RATIO 0.15 --MAX_RECORDS_IN_RAM 150000 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --VALIDATION_STRINGENCY SILENT --CREATE_INDEX true --COMPRESSION_LEVEL 2 --TMP_DIR "\${TMP_DIR}" --CLEAR_DT true --READ_NAME_REGEX null

echo "MarkDuplicates done: \$(date)"
SCRIPT
)

fi

if [[ "$SKIP_MARKDUP" != "true" ]]; then
    DEP_OPT=""
    if [[ -n "$DEPENDENCY_STR" ]]; then
        DEP_OPT="--dependency=afterok:${DEPENDENCY_STR}"
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY RUN] Would submit: markdup_${SAMPLE_ID} (after: ${DEPENDENCY_STR:-none})"
        MARKDUP_JOB_ID="DRY_RUN_MD"
    elif [[ "$DEPENDENCY_STR" == "DRY_RUN_JOB"* ]]; then
        MARKDUP_JOB_ID="DRY_RUN_MD"
    else
        MARKDUP_JOB_ID=$(echo "$MARKDUP_SCRIPT" | sbatch --parsable $DEP_OPT)
        echo "  Submitted MarkDuplicates job: ${MARKDUP_JOB_ID}"
    fi
fi

# ---------------------------------------------------------------------------
# Stop here if --no-bqsr was specified
# ---------------------------------------------------------------------------
if [[ "$NO_BQSR" == "true" ]]; then
    echo ""
    echo "============================================================"
    echo "--no-bqsr set. Pipeline will stop after MarkDuplicates."
    echo "Output: ${DEDUPED_BAM}"
    echo "============================================================"
    exit 0
fi

# ---------------------------------------------------------------------------
# Step 3: BaseRecalibrator
# ---------------------------------------------------------------------------
echo ""
echo "Step 3: Submitting BaseRecalibrator job (dependency on MarkDuplicates)..."

BQSR_SCRIPT=$(cat <<SCRIPT
#!/bin/bash
#SBATCH --job-name=bqsr_${SAMPLE_ID}
#SBATCH --account=${SLURM_ACCOUNT}
${PARTITION_LINE_BQSR}
#SBATCH --output=${LOGS_DIR}/bqsr.%j.out
#SBATCH --error=${LOGS_DIR}/bqsr.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load ${MODULE_GATK}

set -euo pipefail

echo "BaseRecalibrator: \$(date)"

gatk BaseRecalibrator -I "${DEDUPED_BAM}" -R "${REFERENCE}" --known-sites "${KNOWN_SNPS}" --known-sites "${KNOWN_MILLS}" --known-sites "${KNOWN_INDELS}" --known-sites "${KNOWN_1KG}" -O "${RECAL_TABLE}"

echo "BaseRecalibrator done: \$(date)"

echo "ApplyBQSR: \$(date)"

gatk ApplyBQSR -I "${DEDUPED_BAM}" -R "${REFERENCE}" --bqsr-recal-file "${RECAL_TABLE}" --output "${FINAL_BAM}" --preserve-qscores-less-than 6 --use-original-qualities false --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --create-output-bam-index true --read-validation-stringency SILENT

echo "ApplyBQSR done: \$(date)"
echo "Final BAM: ${FINAL_BAM}"
SCRIPT
)

DEP_OPT=""
if [[ -n "$MARKDUP_JOB_ID" ]]; then
    if [[ "$DRY_RUN" == "true" || "$MARKDUP_JOB_ID" == "DRY_RUN_MD" ]]; then
        DEP_OPT=""
    else
        DEP_OPT="--dependency=afterok:${MARKDUP_JOB_ID}"
    fi
fi

if [[ "$DRY_RUN" == "true" || "${MARKDUP_JOB_ID:-}" == "DRY_RUN_MD" ]]; then
    echo "  [DRY RUN] Would submit: bqsr_${SAMPLE_ID} (after: ${MARKDUP_JOB_ID:-none})"
else
    BQSR_JOB_ID=$(echo "$BQSR_SCRIPT" | sbatch --parsable $DEP_OPT)
    echo "  Submitted BQSR job: ${BQSR_JOB_ID}"
fi

# ---------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "Pipeline submitted. Final output when complete:"
echo "  ${FINAL_BAM}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "============================================================"
