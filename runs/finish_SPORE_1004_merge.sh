#!/bin/bash
# ==============================================================================
# Finish Sorting and Merging Mate-CIGAR Fixed Chunks for SPORE_1004
# ------------------------------------------------------------------------------
# Why this script is needed:
#   BamChecker successfully processed all reads and wrote all 16 mate-CIGAR
#   fixed chunks (unsorted_00.bam through unsorted_15.bam) for both Tumor and
#   Normal. It failed only at the very end when trying to invoke host samtools
#   inside Singularity because /cluster/tools was not mounted.
#
# This script finishes the job natively on the host:
#   1. Sorts all 16 chunks per sample in parallel with host samtools.
#   2. Merges the 16 sorted chunks into <sample>.mc.bam.
#   3. Indexes the final BAMs.
#   4. Validates EOF and contig headers with samtools quickcheck.
#
# Resources:
#   Allocates 16 CPUs, 48 GB RAM on himem or pughlab. Takes ~45-60 min total.
# ==============================================================================
#SBATCH --job-name=spore_merge
#SBATCH --partition=himem
#SBATCH --account=deperrotgroup
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams/logs/merge_%j.out
#SBATCH --error=/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams/logs/merge_%j.err

set -euo pipefail

BAM_DIR="/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams"
cd "${BAM_DIR}"

module load samtools/1.20

echo "=================================================================="
echo " Starting Post-BamChecker Sort & Merge for SPORE_1004"
echo " Host:   $(hostname)"
echo " Date:   $(date)"
echo " Directory: ${BAM_DIR}"
echo "=================================================================="

process_sample() {
    local prefix="$1"
    local final_bam="$2"
    local sample_label="$3"

    echo ""
    echo "------------------------------------------------------------------"
    echo " [$(date)] Processing ${sample_label}..."
    echo "------------------------------------------------------------------"

    # Verify all 16 unsorted chunks exist
    local missing=0
    for i in $(seq -w 00 15); do
        local chunk="${prefix}.unsorted_${i}.bam"
        if [[ ! -s "${chunk}" ]]; then
            echo "ERROR: Missing or empty chunk: ${chunk}" >&2
            missing=$((missing + 1))
        fi
    done

    if [[ ${missing} -gt 0 ]]; then
        echo "Aborting ${sample_label}: ${missing} chunks are missing or 0 bytes." >&2
        return 1
    fi

    # 1. Sort all 16 chunks in parallel (4 concurrent jobs, each with 4 threads = 16 CPUs)
    echo "[$(date)] Sorting 16 unsorted chunks..."
    local sorted_files=()
    for i in $(seq -w 00 15); do
        local in_chunk="${prefix}.unsorted_${i}.bam"
        local out_chunk="${prefix}.sorted_${i}.bam"
        sorted_files+=("${out_chunk}")

        if [[ -s "${out_chunk}" ]] && samtools quickcheck "${out_chunk}" 2>/dev/null; then
            echo "  Chunk ${i} already sorted and valid, skipping."
            continue
        fi

        echo "  Sorting chunk ${i} -> ${out_chunk}..."
        samtools sort -@ 4 -m 2G -O bam "${in_chunk}" -o "${out_chunk}" &

        # Keep concurrency capped at 4 parallel jobs
        if [[ $(jobs -r -p | wc -l) -ge 4 ]]; then
            wait -n
        fi
    done
    wait
    echo "[$(date)] All 16 chunks sorted successfully."

    # 2. Merge sorted chunks into final BAM
    echo "[$(date)] Merging 16 sorted chunks into ${final_bam}..."
    samtools merge -@ 16 -f "${final_bam}" "${sorted_files[@]}"
    echo "[$(date)] Merge complete."

    # 3. Index final BAM
    echo "[$(date)] Indexing ${final_bam}..."
    samtools index -@ 8 "${final_bam}"

    # 4. Quickcheck validation
    echo "[$(date)] Validating ${final_bam} with samtools quickcheck..."
    samtools quickcheck "${final_bam}"
    echo "  [OK] ${final_bam} is valid and fully intact!"

    echo "[$(date)] Finished ${sample_label} successfully."
}

# --- Process Tumor ---
TUMOR_PREFIX="${BAM_DIR}/SPORE_1004_Bm_P_WG_CevoMayo4-1_DNA.filter.deduped.recalibrated"
TUMOR_FINAL="${BAM_DIR}/SPORE_1004_Bm_P_WG.mc.bam"
process_sample "${TUMOR_PREFIX}" "${TUMOR_FINAL}" "Tumor (SPORE_1004_Bm_P_WG)"

# --- Process Normal ---
NORMAL_PREFIX="${BAM_DIR}/SPORE_1004_Ly_R_WG_CevoMayo4_Germline.filter.deduped.recalibrated"
NORMAL_FINAL="${BAM_DIR}/SPORE_1004_Ly_R_WG.mc.bam"
process_sample "${NORMAL_PREFIX}" "${NORMAL_FINAL}" "Normal (SPORE_1004_Ly_R_WG)"

echo ""
echo "=================================================================="
echo " [$(date)] Both Tumor and Normal BAMs sorted, merged, and indexed!"
echo " Final output BAMs:"
echo "   ${TUMOR_FINAL}"
echo "   ${NORMAL_FINAL}"
echo "=================================================================="
