#!/bin/bash
# ==============================================================================
# Fix Missing Mate-CIGAR (MC Tag) on Pre-aligned WGS BAMs for HMF Redux
# ------------------------------------------------------------------------------
# Why this is required:
#   Pre-aligned BAMs stripped of alt and unmapped contigs are missing the
#   Mate-CIGAR ('MC') tag on paired reads. HMF Redux aborts immediately with:
#   "FATAL: required mate CIGAR not set, exiting"
#
# Solution:
#   Runs BamChecker from HMF bam-tools container to calculate and add MC tags.
#   Uses host samtools (via module load and -B /cluster/tools) for sort & merge.
# ==============================================================================
#SBATCH --job-name=fix_mate_cigar
#SBATCH --partition=pughlab
#SBATCH --account=pughlab
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams/logs/fix_mate_%j.out
#SBATCH --error=/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams/logs/fix_mate_%j.err

set -euo pipefail

# ------------------------------------------------------------------------------
# Environment & Modules
# ------------------------------------------------------------------------------
module load samtools/1.20
module load singularity/4.3.2

BAM_DIR="/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams"
REF_FASTA="/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta"
CONTAINER="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache/hartwigmedicalfoundation-bam-tools-1.7-beta.2.img"
SAMTOOLS_BIN="$(which samtools)"

mkdir -p "${BAM_DIR}/logs" "${BAM_DIR}/tmp"
export TMPDIR="${BAM_DIR}/tmp"
export TEMP="${BAM_DIR}/tmp"
export TMP="${BAM_DIR}/tmp"
export SINGULARITY_TMPDIR="${BAM_DIR}/tmp"
export APPTAINER_TMPDIR="${BAM_DIR}/tmp"
export SINGULARITYENV_TMPDIR="${BAM_DIR}/tmp"
export SINGULARITYENV_TEMP="${BAM_DIR}/tmp"
export SINGULARITYENV_TMP="${BAM_DIR}/tmp"

# Provide ample JVM memory for BamChecker
export JAVA_TOOL_OPTIONS="-Xmx48G -Djava.io.tmpdir=${TMPDIR}"

# Ensure Sambamba sorts write intermediate chunks strictly to project tmpdir
# (BamChecker calls 'bamtool sort ...' without --tmpdir; this wrapper guarantees
# --tmpdir is explicitly passed on the CLI, and filename contains 'sambamba'
# so BamToolName.fromPath() recognizes it properly)
SAMBAMBA_WRAPPER="${BAM_DIR}/sambamba"
cat << 'EOF' > "${SAMBAMBA_WRAPPER}"
#!/bin/bash
if [[ "${1:-}" == "sort" ]]; then
    exec /tools/sambamba sort --tmpdir="${TMPDIR:-/tmp}" "${@:2}"
else
    exec /tools/sambamba "$@"
fi
EOF
chmod +x "${SAMBAMBA_WRAPPER}"

echo "=================================================================="
echo " Starting BamChecker Mate-CIGAR Fix"
echo " Host:         $(hostname)"
echo " Date:         $(date)"
echo " Samtools:     ${SAMTOOLS_BIN}"
echo " Container:    ${CONTAINER}"
echo " Ref Genome:   ${REF_FASTA}"
echo "=================================================================="

# Allocate half of total CPUs to each concurrent sample
TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-32}"
THREADS_PER_SAMPLE=$(( TOTAL_CPUS / 2 ))
[[ $THREADS_PER_SAMPLE -lt 1 ]] && THREADS_PER_SAMPLE=1

fix_bam() {
    local in_bam="$1"
    local sample_name="$2"
    local threads="$3"
    local out_bam="${BAM_DIR}/${sample_name}.mc.bam"

    echo "[$(date)] Starting BamChecker on ${sample_name} with ${threads} threads..."
    singularity exec \
        -B /cluster/projects:/cluster/projects \
        -B "${TMPDIR}:/tmp" \
        -B "${TMPDIR}:/var/tmp" \
        --env TMPDIR="${TMPDIR}",TEMP="${TMPDIR}",TMP="${TMPDIR}" \
        "${CONTAINER}" \
        bamtools com.hartwig.hmftools.bamtools.checker.BamChecker \
            -bam_file "${in_bam}" \
            -ref_genome "${REF_FASTA}" \
            -bamtool "${SAMBAMBA_WRAPPER}" \
            -output_dir "${BAM_DIR}" \
            -output_bam "${out_bam}" \
            -threads "${threads}"

    echo "[$(date)] Indexing ${out_bam} with ${threads} threads..."
    samtools index -@ "${threads}" "${out_bam}"

    echo "[$(date)] Validating ${out_bam} with samtools quickcheck..."
    samtools quickcheck "${out_bam}"
    echo "[$(date)] Successfully finished and validated ${sample_name}."
}

# Run tumor and normal concurrently (${THREADS_PER_SAMPLE} threads each)
fix_bam "${BAM_DIR}/SPORE_1004_Bm_P_WG_CevoMayo4-1_DNA.filter.deduped.recalibrated.bam" "SPORE_1004_Bm_P_WG" "${THREADS_PER_SAMPLE}" &
PID_TUMOR=$!

fix_bam "${BAM_DIR}/SPORE_1004_Ly_R_WG_CevoMayo4_Germline.filter.deduped.recalibrated.bam" "SPORE_1004_Ly_R_WG" "${THREADS_PER_SAMPLE}" &
PID_NORMAL=$!

echo "Dispatched Tumor job (PID: ${PID_TUMOR}) and Normal job (PID: ${PID_NORMAL}). Waiting for completion..."

wait "${PID_TUMOR}"
wait "${PID_NORMAL}"

echo ""
echo "=================================================================="
echo " [$(date)] All BAMs fixed, indexed, and validated successfully!"
echo " Outputs:"
echo "   ${BAM_DIR}/SPORE_1004_Bm_P_WG.mc.bam"
echo "   ${BAM_DIR}/SPORE_1004_Ly_R_WG.mc.bam"
echo "=================================================================="

