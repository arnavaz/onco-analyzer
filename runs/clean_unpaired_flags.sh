#!/bin/bash
# ==============================================================================
# Sanitize SAM Flags & Mate Fields on Already Coordinate-Sorted, MC-Tagged BAMs
# ------------------------------------------------------------------------------
# Why this is needed:
#   When BAMs have had alt contigs and unmapped reads stripped prior to
#   samtools fixmate, singletons whose mates were deleted get bit 0x1 (PAIRED)
#   cleared by fixmate while retaining paired-end flags (0x2, 0x8, 0x20, 0x40,
#   0x80), RNEXT, PNEXT, and TLEN. Strict HTSJDK readers in HMF tools (CIDER,
#   BAMTOOLS queryUnmapped, ESVEE FragmentSizeDistribution) abort with:
#     "ERROR::INVALID_FLAG_SECOND_OF_PAIR: Second of pair flag should not be set
#      for unpaired read."
#   Lenient readers (TEAL_PREP) emit multi-GB validation warnings to .command.err.
#
# Why this is fast:
#   Because input .mc.bam files are ALREADY coordinate-sorted and MC-tagged,
#   this runs a single linear streaming pass (no collate, no sort) preserving
#   exact coordinate order and paired reads' MC tags.
#
# Usage (submit tumor and normal as SEPARATE jobs so one cannot affect the other):
#   sbatch -A pughlab -p pughlab -J clean_tumor  runs/clean_unpaired_flags.sh <in.mc.bam> <out.clean.mc.bam>
#   sbatch -A pughlab -p pughlab -J clean_normal runs/clean_unpaired_flags.sh <in.mc.bam> <out.clean.mc.bam>
#
# Smoke test first (~30s; slices CIDER IGH region chr14:105580000-106880000):
#   sbatch -A pughlab -p pughlab --time=00:30:00 --cpus-per-task=8 --mem=16G \
#     runs/clean_unpaired_flags.sh --test <in.mc.bam> <test_out.clean.mc.bam>
# ==============================================================================
#SBATCH --job-name=clean_bam_flags
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=clean_bam_%x_%j.out
#SBATCH --error=clean_bam_%x_%j.err

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  clean_unpaired_flags.sh [--test] [--verify-htsjdk | --no-verify-htsjdk] <input.mc.bam> <output.clean.mc.bam> [tmp_dir]

Options:
  --test               Slice chr14:105580000-106880000 (CIDER IGH locus) for fast smoke test
  --verify-htsjdk      Verify output with cached HMF Cider (HTSJDK STRICT) container (default)
  --no-verify-htsjdk   Skip Singularity/HMF HTSJDK STRICT smoke check
  -h, --help           Show this help message
EOF
    exit 1
}

TEST_MODE=0
VERIFY_HTSJDK=1
declare -a POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test)
            TEST_MODE=1
            shift
            ;;
        --verify-htsjdk)
            VERIFY_HTSJDK=1
            shift
            ;;
        --no-verify-htsjdk)
            VERIFY_HTSJDK=0
            shift
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            usage
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL[@]}"
[[ $# -lt 2 || $# -gt 3 ]] && usage

# Load required cluster modules BEFORE calling samtools / singularity
module load samtools/1.20 2>/dev/null || module load samtools

IN_BAM="$1"
OUT_BAM="$2"
OUT_DIR="$(dirname "${OUT_BAM}")"
SAMPLE="$(basename "${OUT_BAM}" .bam)"
TMP_BASE="${3:-${OUT_DIR}/tmp}/${SAMPLE}.clean.${SLURM_JOB_ID:-$$}"
TEST_REGION="chr14:105580000-106880000"

# Resolve AWK sanitizer safely even when SLURM executes a spool copy of this script
DEFAULT_REPO="/cluster/projects/pughlab/src/oncoanalyser_pughlab"
SCRIPT_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." 2>/dev/null && pwd || true)"
SOURCE_REPO="${SOURCE_REPO:-${SCRIPT_REPO}}"
if [[ ! -f "${SOURCE_REPO}/scripts/sanitize_sam_flags.awk" ]]; then
    SOURCE_REPO="${DEFAULT_REPO}"
fi
AWK_SCRIPT="${SOURCE_REPO}/scripts/sanitize_sam_flags.awk"

[[ -s "${IN_BAM}" ]] || { echo "ERROR: Input BAM not found or empty: ${IN_BAM}" >&2; exit 1; }
[[ ! -e "${OUT_BAM}" ]] || { echo "ERROR: Output already exists, refusing to overwrite: ${OUT_BAM}" >&2; exit 1; }
[[ -f "${AWK_SCRIPT}" ]] || { echo "ERROR: Required AWK sanitizer not found: ${AWK_SCRIPT}" >&2; exit 1; }

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found after module load" >&2; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "ERROR: awk not found in PATH" >&2; exit 1; }

THREADS="${SLURM_CPUS_PER_TASK:-16}"
READ_T=$(( THREADS / 4 ));  [[ ${READ_T} -lt 1 ]] && READ_T=1
[[ ${READ_T} -gt 4 ]] && READ_T=4
WRITE_T=$(( THREADS - READ_T - 2 )); [[ ${WRITE_T} -lt 1 ]] && WRITE_T=1

mkdir -p "${OUT_DIR}" "${TMP_BASE}"
export TMPDIR="${TMP_BASE}"
export TEMP="${TMP_BASE}"
export TMP="${TMP_BASE}"

PARTIAL="${OUT_BAM%.bam}.partial.bam"
AUDIT_LOG="${TMP_BASE}/sanitize_audit.log"

# Verify index exists for --test region slicing and idxstats checks
if [[ ! -s "${IN_BAM}.bai" && ! -s "${IN_BAM%.bam}.bai" ]]; then
    echo "[$(date)] Input BAM index not found; indexing ${IN_BAM} with ${THREADS} threads..."
    samtools index -@ "${THREADS}" "${IN_BAM}"
fi

echo "=================================================================="
echo " Clean Unpaired Flags & Mate Fields (Single-Pass Streaming)"
echo " Host:      $(hostname)   Date: $(date)"
echo " Input:     ${IN_BAM}"
echo " Output:    ${OUT_BAM}"
echo " Partial:   ${PARTIAL}"
echo " Temp:      ${TMP_BASE}"
echo " AWK script:${AWK_SCRIPT}"
echo " Threads:   read=${READ_T} write=${WRITE_T} index=${THREADS}"
echo " Mode:      $([[ ${TEST_MODE} -eq 1 ]] && echo "TEST (${TEST_REGION})" || echo "FULL")"
echo " samtools:  $(samtools --version | head -1)"
echo "=================================================================="

stream_input() {
    if [[ ${TEST_MODE} -eq 1 ]]; then
        samtools view -h -@ "${READ_T}" "${IN_BAM}" "${TEST_REGION}"
    else
        samtools view -h -@ "${READ_T}" "${IN_BAM}"
    fi
}

echo "[$(date)] Streaming SAM flag sanitization into ${PARTIAL}..."
stream_input \
    | LC_ALL=C awk -f "${AWK_SCRIPT}" 2> "${AUDIT_LOG}" \
    | samtools view -@ "${WRITE_T}" -b -o "${PARTIAL}" -

echo "[$(date)] Sanitizer audit summary:"
cat "${AUDIT_LOG}"

TOTAL_RECS=$(awk -F'total=' 'NF>1 {split($2,a,"\t"); print a[1]}' "${AUDIT_LOG}" | tail -n 1 || echo "0")
if [[ ! "${TOTAL_RECS}" =~ ^[0-9]+$ || "${TOTAL_RECS}" -eq 0 ]]; then
    echo "ERROR: Sanitizer processed 0 records; aborting and leaving ${PARTIAL} for inspection." >&2
    exit 1
fi

echo "[$(date)] Indexing ${PARTIAL} with ${THREADS} threads..."
samtools index -@ "${THREADS}" "${PARTIAL}"
samtools quickcheck -v "${PARTIAL}"

if [[ ${TEST_MODE} -eq 0 ]]; then
    echo "[$(date)] Verifying total read count conservation (idxstats)..."
    IN_N=$(samtools idxstats "${IN_BAM}"  | awk '{s+=$3+$4} END {print s+0}')
    OUT_N=$(samtools idxstats "${PARTIAL}" | awk '{s+=$3+$4} END {print s+0}')
    echo "  input=${IN_N} output=${OUT_N}"
    if [[ "${IN_N}" -ne "${OUT_N}" ]]; then
        echo "ERROR: Read count mismatch (in=${IN_N}, out=${OUT_N}). Leaving ${PARTIAL} for inspection." >&2
        exit 1
    fi
fi

echo "[$(date)] Verifying MC tag retention on properly paired mapped primaries (up to 1M reads)..."
read -r CHECKED MISSING < <(
    { samtools view -f 1 -F 0x90C "${PARTIAL}" 2>/dev/null || true; } | head -n 1000000 \
        | awk '{n++; if ($0 !~ /\tMC:Z:/) m++} END {print n+0, m+0}'
)
echo "  checked=${CHECKED} missing_MC=${MISSING}"
if [[ "${CHECKED}" -eq 0 || "${MISSING}" -gt $(( CHECKED / 1000 )) ]]; then
    echo "ERROR: MC validation failed (>0.1% of paired mapped reads lack MC). Leaving ${PARTIAL} for inspection." >&2
    exit 1
fi

echo "[$(date)] Re-auditing head, chr14 IGH locus, chrY tail, and unplaced (*) records for zero violations..."
{
    { samtools view -h "${PARTIAL}" 2>/dev/null || true; } | head -n 200500
    { samtools view "${PARTIAL}" "chr14:105580000-106880000" 2>/dev/null || true; }
    if [[ ${TEST_MODE} -eq 0 ]]; then
        { samtools view "${PARTIAL}" "chrY:55000000-57227415" 2>/dev/null || true; }
        { samtools view "${PARTIAL}" "*" 2>/dev/null || true; }
    fi
} | LC_ALL=C awk -v AUDIT_ONLY=1 -f "${AWK_SCRIPT}"

if [[ ${VERIFY_HTSJDK} -eq 1 ]]; then
    CIDER_IMG="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache/depot.galaxyproject.org-singularity-hmftools-cider-1.2--hdfd78af_0.img"
    REF_FASTA="/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta"

    module load singularity/4.3.2 2>/dev/null || module load singularity 2>/dev/null || true
    if command -v singularity >/dev/null 2>&1; then
        if [[ -s "${CIDER_IMG}" && -s "${REF_FASTA}" ]]; then
            echo "[$(date)] Running HTSJDK STRICT verification via HMF Cider (${CIDER_IMG})..."
            CIDER_VERIFY_DIR="${TMP_BASE}/cider_verify"
            mkdir -p "${CIDER_VERIFY_DIR}"
            declare -a SING_BINDS=(
                -B "/cluster/projects:/cluster/projects"
                -B "${TMP_BASE}:/tmp"
                -B "${TMP_BASE}:/var/tmp"
            )
            if [[ "${OUT_DIR}" != /cluster/projects/* && "${OUT_DIR}" != "/cluster/projects" ]]; then
                SING_BINDS+=(-B "${OUT_DIR}:${OUT_DIR}")
            fi
            if [[ "${TMP_BASE}" != /cluster/projects/* && "${TMP_BASE}" != "/cluster/projects" && "${TMP_BASE}" != "${OUT_DIR}"* ]]; then
                SING_BINDS+=(-B "${TMP_BASE}:${TMP_BASE}")
            fi
            singularity exec \
                "${SING_BINDS[@]}" \
                --env TMPDIR="${TMP_BASE}",TEMP="${TMP_BASE}",TMP="${TMP_BASE}" \
                "${CIDER_IMG}" \
                cider \
                    -Xmx12G \
                    com.hartwig.hmftools.cider.CiderApplication \
                    -sample "${SAMPLE}" \
                    -bam "${PARTIAL}" \
                    -ref_genome_version 38 \
                    -ref_genome "${REF_FASTA}" \
                    -threads 4 \
                    -output_dir "${CIDER_VERIFY_DIR}"
            echo "[$(date)] HTSJDK STRICT verification via Cider PASSED!"
            rm -f "${CIDER_VERIFY_DIR}/"* 2>/dev/null || true
            rmdir "${CIDER_VERIFY_DIR}" 2>/dev/null || true
        else
            echo "[$(date)] NOTE: Cider container or reference FASTA not found locally; skipping Cider HTSJDK check."
        fi
    else
        echo "[$(date)] NOTE: singularity not available; skipping Cider HTSJDK check."
    fi
fi

mv -f "${PARTIAL}" "${OUT_BAM}"
mv -f "${PARTIAL}.bai" "${OUT_BAM}.bai"
rm -f "${AUDIT_LOG}" 2>/dev/null || true
rmdir "${TMP_BASE}" 2>/dev/null || echo "NOTE: temp dir not empty, left in place: ${TMP_BASE}"

echo "[$(date)] Done: ${OUT_BAM}"
