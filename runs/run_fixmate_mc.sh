#!/bin/bash
# ==============================================================================
# Add Mate-CIGAR (MC) tags to a pre-aligned, coordinate-sorted BAM
# ------------------------------------------------------------------------------
# Why: REDUX (even in -bqr_jitter_msi_only / generate_redux_tsvs_only mode)
#      samples the first reads and exits with "required mate CIGAR not set" when
#      MC tags are absent (e.g. GATK/Picard BAMs from external providers).
#
# How: a single streamed samtools pipeline (host samtools, no container/Java):
#        collate (group mates) -> fixmate (adds MC) -> sort (coordinate) -> index
#      Intermediate data is passed uncompressed through pipes; only temp spill
#      files from collate/sort touch disk (~1x input BAM size each).
#
# Usage (submit tumor and normal as SEPARATE jobs so one can't kill the other):
#   sbatch -A pughlab -p pughlab -J mc_tumor  run_fixmate_mc.sh <in.bam> <out.mc.bam>
#   sbatch -A pughlab -p pughlab -J mc_normal run_fixmate_mc.sh <in.bam> <out.mc.bam>
#
# Smoke test first (~minutes; slices one region via the index, no full scan):
#   sbatch -A pughlab -p pughlab --time=01:00:00 --cpus-per-task=8 --mem=16G \
#     run_fixmate_mc.sh --test <in.bam> <test_out.mc.bam>
# ==============================================================================
#SBATCH --job-name=fixmate_mc
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=fixmate_mc_%x_%j.out
#SBATCH --error=fixmate_mc_%x_%j.err

set -euo pipefail

usage() {
    echo "Usage: $0 [--test] <input.bam> <output.mc.bam> [tmp_dir]" >&2
    exit 1
}

TEST_MODE=0
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=1
    shift
fi
[[ $# -lt 2 || $# -gt 3 ]] && usage

IN_BAM="$1"
OUT_BAM="$2"
OUT_DIR="$(dirname "${OUT_BAM}")"
SAMPLE="$(basename "${OUT_BAM}" .bam)"
TMP_BASE="${3:-${OUT_DIR}/tmp}/${SAMPLE}.${SLURM_JOB_ID:-$$}"
TEST_REGION="chr22:20000000-22000000"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Split threads across the three stages (sort does the heavy compression)
SORT_T=$(( THREADS / 2 ));    [[ ${SORT_T} -lt 1 ]] && SORT_T=1
COLLATE_T=$(( THREADS / 4 )); [[ ${COLLATE_T} -lt 1 ]] && COLLATE_T=1
FIXMATE_T=$(( THREADS / 8 )); [[ ${FIXMATE_T} -lt 1 ]] && FIXMATE_T=1
# Per-thread sort memory: ~50% of the allocation spread over sort threads.
# The rest is headroom for collate, which loads one whole temp bucket into RAM
# at a time (bucket ~= input_size / COLLATE_N compressed, ~4x that in memory).
MEM_MB="${SLURM_MEM_PER_NODE:-16384}"
SORT_MEM_MB=$(( MEM_MB * 5 / 10 / SORT_T )); [[ ${SORT_MEM_MB} -lt 768 ]] && SORT_MEM_MB=768

# Number of collate temp buckets. Default 64 gave 6.7 GB buckets (~25-30 GB RAM)
# on 400 GB WGS BAMs and OOM-killed the job. 512 -> ~0.8 GB buckets (~3 GB RAM).
# collate holds all buckets open at once, so stay under the open-file limit.
COLLATE_N="${COLLATE_N:-512}"
NOFILE_LIMIT="$(ulimit -n)"
if [[ "${NOFILE_LIMIT}" != "unlimited" && $(( COLLATE_N + 64 )) -gt ${NOFILE_LIMIT} ]]; then
    ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
    NOFILE_LIMIT="$(ulimit -n)"
    if [[ "${NOFILE_LIMIT}" != "unlimited" && $(( COLLATE_N + 64 )) -gt ${NOFILE_LIMIT} ]]; then
        COLLATE_N=$(( NOFILE_LIMIT - 64 ))
        echo "WARNING: open-file limit ${NOFILE_LIMIT}; reducing COLLATE_N to ${COLLATE_N}" >&2
    fi
fi

[[ -s "${IN_BAM}" ]] || { echo "ERROR: input BAM not found: ${IN_BAM}" >&2; exit 1; }
[[ -e "${OUT_BAM}" ]] && { echo "ERROR: output already exists, refusing to overwrite: ${OUT_BAM}" >&2; exit 1; }

DEFAULT_REPO="/cluster/projects/pughlab/src/oncoanalyser_pughlab"
SCRIPT_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." 2>/dev/null && pwd || true)"
SOURCE_REPO="${SOURCE_REPO:-${SCRIPT_REPO}}"
if [[ ! -f "${SOURCE_REPO}/scripts/sanitize_sam_flags.awk" ]]; then
    SOURCE_REPO="${DEFAULT_REPO}"
fi
AWK_SCRIPT="${SOURCE_REPO}/scripts/sanitize_sam_flags.awk"
[[ -f "${AWK_SCRIPT}" ]] || { echo "ERROR: SAM flag sanitizer not found: ${AWK_SCRIPT}" >&2; exit 1; }

module load samtools/1.20

mkdir -p "${TMP_BASE}"
export TMPDIR="${TMP_BASE}"
PARTIAL="${OUT_BAM%.bam}.partial.bam"
AUDIT_LOG="${TMP_BASE}/sanitize_audit.log"

echo "=================================================================="
echo " fixmate MC tagging + unpaired SAM flag sanitization"
echo " Host:      $(hostname)   Date: $(date)"
echo " Input:     ${IN_BAM}"
echo " Output:    ${OUT_BAM}"
echo " Temp:      ${TMP_BASE}"
echo " AWK script:${AWK_SCRIPT}"
echo " Threads:   collate=${COLLATE_T} fixmate=${FIXMATE_T} sort=${SORT_T} (-m ${SORT_MEM_MB}M)"
echo " Collate:   ${COLLATE_N} temp buckets (open-file limit: ${NOFILE_LIMIT})"
echo " Mode:      $([[ ${TEST_MODE} -eq 1 ]] && echo "TEST (${TEST_REGION})" || echo FULL)"
echo " samtools:  $(samtools --version | head -1)"
echo "=================================================================="

collate_input() {
    if [[ ${TEST_MODE} -eq 1 ]]; then
        samtools view -u -@ "${COLLATE_T}" "${IN_BAM}" "${TEST_REGION}" \
          | samtools collate -@ "${COLLATE_T}" -n "${COLLATE_N}" -O -u -T "${TMP_BASE}/collate" -
    else
        samtools collate -@ "${COLLATE_T}" -n "${COLLATE_N}" -O -u -T "${TMP_BASE}/collate" "${IN_BAM}"
    fi
}

echo "[$(date)] Running collate | fixmate | sanitize SAM flags | sort ..."
collate_input \
  | samtools fixmate -@ "${FIXMATE_T}" -u - - \
  | samtools view -@ 2 -h - \
  | LC_ALL=C awk -f "${AWK_SCRIPT}" 2> "${AUDIT_LOG}" \
  | samtools sort -@ "${SORT_T}" -m "${SORT_MEM_MB}M" -T "${TMP_BASE}/sort" -o "${PARTIAL}" -

echo "[$(date)] Sanitizer audit summary:"
cat "${AUDIT_LOG}"

echo "[$(date)] Indexing ..."
samtools index -@ "${THREADS}" "${PARTIAL}"
samtools quickcheck "${PARTIAL}"

# Validation: paired, mapped reads with mapped mates (excl. secondary/supp)
# should all carry MC. Orphans whose mate was stripped upstream have unpaired
# flags sanitized (0x2|0x8|0x20|0x40|0x80 cleared, RNEXT=*, PNEXT=0, TLEN=0)
# for HTSJDK ValidationStringency.STRICT compatibility.
echo "[$(date)] Validating MC tags on first 1M properly mapped paired reads ..."
read -r CHECKED MISSING < <(
    { samtools view -f 1 -F 0x90C "${PARTIAL}" 2>/dev/null || true; } | head -n 1000000 \
      | awk '{n++; if ($0 !~ /\tMC:Z:/) m++} END {print n+0, m+0}'
)
echo "  checked=${CHECKED} missing_MC=${MISSING}"
if [[ "${CHECKED}" -eq 0 || "${MISSING}" -gt $(( CHECKED / 1000 )) ]]; then
    echo "ERROR: MC validation failed (>0.1% of checked reads lack MC). Leaving ${PARTIAL} for inspection." >&2
    exit 1
fi

echo "[$(date)] Auditing SAM flags across head, chr14 IGH, chrY tail, and unplaced (*) ..."
{
    { samtools view -h "${PARTIAL}" 2>/dev/null || true; } | head -n 200500
    { samtools view "${PARTIAL}" "chr14:105580000-106880000" 2>/dev/null || true; }
    if [[ ${TEST_MODE} -eq 0 ]]; then
        { samtools view "${PARTIAL}" "chrY:55000000-57227415" 2>/dev/null || true; }
        { samtools view "${PARTIAL}" "*" 2>/dev/null || true; }
    fi
} | LC_ALL=C awk -v AUDIT_ONLY=1 -f "${AWK_SCRIPT}"

echo "[$(date)] Read counts (input vs output, flagstat totals) ..."
if [[ ${TEST_MODE} -eq 0 ]]; then
    IN_N=$(samtools idxstats "${IN_BAM}"  | awk '{s+=$3+$4} END {print s}')
    OUT_N=$(samtools idxstats "${PARTIAL}" | awk '{s+=$3+$4} END {print s}')
    echo "  input=${IN_N} output=${OUT_N}"
    if [[ "${IN_N}" -ne "${OUT_N}" ]]; then
        echo "ERROR: read count mismatch. Leaving ${PARTIAL} for inspection." >&2
        exit 1
    fi
fi

mv "${PARTIAL}" "${OUT_BAM}"
mv "${PARTIAL}.bai" "${OUT_BAM}.bai"
rm -f "${AUDIT_LOG}" 2>/dev/null || true
rmdir "${TMP_BASE}" 2>/dev/null || echo "NOTE: temp dir not empty, left in place: ${TMP_BASE}"

echo "[$(date)] Done: ${OUT_BAM}"
