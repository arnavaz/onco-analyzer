#!/bin/bash
# ==============================================================================
# Pipeline Runner: SPORE_1004 (Whole Genome Tumor-Normal Analysis)
# Target: HMF Oncoanalyser v3.1 (WGTS Mode, GRCh38_hmf)
# Input: Pre-aligned, duplicate-marked, recalibrated WGS BAMs (alt & unmapped stripped)
# ==============================================================================
set -euo pipefail

DRY_RUN="${DRY_RUN:-false}"
if [[ "${1:-}" == "--preflight-only" ]]; then
    DRY_RUN="true"
elif [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    echo "Usage: bash $0 [--preflight-only]" >&2
    exit 0
fi

# Load required cluster modules BEFORE calling samtools
module load samtools/1.20 2>/dev/null || module load samtools

# ------------------------------------------------------------------------------
# 1. Configuration & Paths
# ------------------------------------------------------------------------------
PROJECT_ID="SPORE_1004"
SAMPLE_LABEL="${PROJECT_ID}_Bm_P_WG"
NORMAL_LABEL="${PROJECT_ID}_Ly_R_WG"

ACCOUNT="${ACCOUNT:-pughlab}"
PARTITION="${PARTITION:-pughlab}"
WGTS_TIME="${WGTS_TIME:-5-00:00:00}"
SOURCE_REPO="${SOURCE_REPO:-/cluster/projects/pughlab/src/oncoanalyser_pughlab}"

# Unarchived BAM storage on cluster filesystem
BAM_DIR="/cluster/projects/stewartlab/SPORE_1004_response/Unarchived_bams"
TUMOUR_BAM="${BAM_DIR}/SPORE_1004_Bm_P_WG_CevoMayo4-1_DNA.filter.deduped.recalibrated.bam"
NORMAL_BAM="${BAM_DIR}/SPORE_1004_Ly_R_WG_CevoMayo4_Germline.filter.deduped.recalibrated.bam"

# Prefer sanitized mate-CIGAR fixed BAMs (.clean.mc.bam), then .mc.bam
if [[ -s "${BAM_DIR}/${SAMPLE_LABEL}.clean.mc.bam" ]]; then
    TUMOUR_BAM="${BAM_DIR}/${SAMPLE_LABEL}.clean.mc.bam"
    echo "Using flag-sanitized mate-CIGAR fixed tumour BAM: ${TUMOUR_BAM}"
elif [[ -s "${BAM_DIR}/${SAMPLE_LABEL}.mc.bam" ]]; then
    TUMOUR_BAM="${BAM_DIR}/${SAMPLE_LABEL}.mc.bam"
    echo "Using mate-CIGAR fixed tumour BAM: ${TUMOUR_BAM}"
fi

if [[ -s "${BAM_DIR}/${NORMAL_LABEL}.clean.mc.bam" ]]; then
    NORMAL_BAM="${BAM_DIR}/${NORMAL_LABEL}.clean.mc.bam"
    echo "Using flag-sanitized mate-CIGAR fixed normal BAM: ${NORMAL_BAM}"
elif [[ -s "${BAM_DIR}/${NORMAL_LABEL}.mc.bam" ]]; then
    NORMAL_BAM="${BAM_DIR}/${NORMAL_LABEL}.mc.bam"
    echo "Using mate-CIGAR fixed normal BAM: ${NORMAL_BAM}"
fi

# Oncoanalyser output & sample sheet locations
OUTDIR="/cluster/projects/stewartlab/SPORE_1004_response/Oncoanalyser"
SHEETS_DIR="${OUTDIR}/sample_sheets"
CONFIG_TSV="${SHEETS_DIR}/config_${PROJECT_ID}.tsv"
PATIENT_OUTDIR="${OUTDIR}/${PROJECT_ID}"
LOGS_DIR="${PATIENT_OUTDIR}/logs"
TMP_DIR="${PATIENT_OUTDIR}/tmp"

mkdir -p "${BAM_DIR}" "${SHEETS_DIR}" "${LOGS_DIR}" "${TMP_DIR}"

# ------------------------------------------------------------------------------
# 2. Download Pre-aligned BAMs & Indexes from UHN On-Prem ECS
# ------------------------------------------------------------------------------
# Note: Dell EMC ECS throttles multi-threaded 'aws s3 cp' leading to timeouts.
# We use single-stream 'aws s3api get-object' with --cli-read-timeout 0.
ECS_BUCKET="myelomagroup-ecsresgen"
TUMOUR_KEY="all_DNA_WG_bams/SPORE_1004_Bm_P_WG_CevoMayo4-1_DNA.filter.deduped.recalibrated.bam"
NORMAL_KEY="all_DNA_WG_bams/SPORE_1004_Ly_R_WG_CevoMayo4_Germline.filter.deduped.recalibrated.bam"

echo "=========================================================="
echo " [Step 1/5] Checking / Downloading BAMs from ECS..."
echo "=========================================================="

# Helper function to download an ECS object if missing or empty
download_ecs_file() {
    local key="$1"
    local dest="$2"
    if [[ -s "$dest" ]]; then
        echo "  [OK] Already downloaded: $dest ($(du -h "$dest" | cut -f1))"
    else
        echo "  --> Downloading $key -> $dest..."
        aws --endpoint-url=https://ecs.uhn.ca/ --profile=myeloma --cli-read-timeout 0 \
            s3api get-object --bucket "$ECS_BUCKET" --key "$key" "$dest"
    fi
}

download_ecs_file "$TUMOUR_KEY" "$TUMOUR_BAM"
download_ecs_file "$NORMAL_KEY" "$NORMAL_BAM"

# Download or generate BAM indexes (.bai)
for bam in "$TUMOUR_BAM" "$NORMAL_BAM"; do
    bai="${bam}.bai"
    alt_bai="${bam%.bam}.bai"
    if [[ -s "$bai" || -s "$alt_bai" ]]; then
        echo "  [OK] BAM index exists for $(basename "$bam")"
    elif [[ "$bam" == *.clean.mc.bam || "$bam" == *.mc.bam ]]; then
        echo "  --> Local fixed BAM index missing. Building index with samtools..."
        samtools index -@ 8 "$bam"
    else
        # Determine remote S3 key for raw recalibrated BAM index
        if [[ "$bam" == "$TUMOUR_BAM" ]]; then
            remote_bai="${TUMOUR_KEY}.bai"
        else
            remote_bai="all_DNA_WG_bams/SPORE_1004_Ly_R_WG_CevoMayo4_Germline.filter.deduped.recalibrated.bai"
        fi

        echo "  --> Downloading index: $remote_bai -> $bai..."
        if ! aws --endpoint-url=https://ecs.uhn.ca/ --profile=myeloma --cli-read-timeout 0 \
            s3api get-object --bucket "$ECS_BUCKET" --key "$remote_bai" "$bai" 2>/dev/null; then
            echo "  --> Remote index not found. Building index with samtools locally..."
            samtools index -@ 4 "$bam"
        fi
    fi
done

# ------------------------------------------------------------------------------
# 3. BAM Sanity & Preflight Compatibility Checks
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " [Step 2/5] Running Preflight Compatibility Checks on BAMs..."
echo "=========================================================="

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found after 'module load samtools'." >&2
    exit 1
}

echo "  Checking BAM truncation / EOF headers with samtools quickcheck..."
samtools quickcheck -v "$TUMOUR_BAM" "$NORMAL_BAM" || {
    echo "ERROR: BAM quickcheck failed! One or more BAMs are corrupted or incomplete." >&2
    exit 1
}

# Verify chromosome naming format (GRCh38_hmf expects 'chr' prefix: chr1..chr22, chrX, chrY)
# Note: 2>/dev/null prevents samtools SIGPIPE warning when grep -m 1 exits early
chr_prefix_check=$((samtools view -H "$TUMOUR_BAM" 2>/dev/null || true) | grep -m 1 "^@SQ" | grep -o "SN:chr" || true)
if [[ "$chr_prefix_check" == "SN:chr" ]]; then
    echo "  [OK] Chromosome contigs use 'chr' prefix (matches HMF GRCh38_hmf reference)."
else
    echo "  WARNING: Chromosome contigs in BAM do NOT appear to use 'chr' prefix!"
    echo "           HMF GRCh38_hmf expects 'chr1', 'chr2', ..., 'chrX'."
    echo "           Please verify contig names: samtools view -H $TUMOUR_BAM | grep '^@SQ' | head -n 5"
fi

echo "  Checking Mate-CIGAR (MC) tags and HTSJDK STRICT SAM flag compatibility..."
SCRIPT_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." 2>/dev/null && pwd || true)"
AWK_SCRIPT="${SOURCE_REPO}/scripts/sanitize_sam_flags.awk"
if [[ ! -f "${AWK_SCRIPT}" && -n "${SCRIPT_REPO}" ]]; then
    AWK_SCRIPT="${SCRIPT_REPO}/scripts/sanitize_sam_flags.awk"
fi
[[ -f "${AWK_SCRIPT}" ]] || {
    echo "ERROR: SAM flag audit script not found: ${AWK_SCRIPT}" >&2
    exit 1
}

for label_bam in "Tumour:${TUMOUR_BAM}:${SAMPLE_LABEL}" "Normal:${NORMAL_BAM}:${NORMAL_LABEL}"; do
    IFS=':' read -r role bam_path sample_lbl <<< "${label_bam}"
    read -r mc_checked mc_missing < <(
        { samtools view -f 1 -F 0x90C "${bam_path}" 2>/dev/null || true; } \
            | head -n 100000 \
            | awk '{n++; if ($0 !~ /\tMC:Z:/) m++} END {print n+0, m+0}'
    )
    echo "  [${role}] MC tag sample: checked=${mc_checked}, missing_MC=${mc_missing} ($(basename "${bam_path}"))"
    if [[ "${mc_checked}" -eq 0 || "${mc_missing}" -gt $(( mc_checked / 100 )) ]]; then
        echo "ERROR: ${role} BAM lacks MC tags (>1% missing). REDUX will fail." >&2
        echo "       Generate MC-tagged + flag-sanitized BAM first:" >&2
        echo "         sbatch -A ${ACCOUNT} -p ${PARTITION} -J mc_${role} ${SOURCE_REPO}/runs/run_fixmate_mc.sh ${bam_path} ${BAM_DIR}/${sample_lbl}.clean.mc.bam" >&2
        exit 1
    fi

    # Audit head 200k records + CIDER IGH locus (chr14) + chrY tail + unplaced (*) for HTSJDK STRICT compliance
    set +e
    {
        { samtools view -h "${bam_path}" 2>/dev/null || true; } | head -n 200500
        { samtools view "${bam_path}" "chr14:105580000-106880000" 2>/dev/null || true; }
        { samtools view "${bam_path}" "chrY:55000000-57227415" 2>/dev/null || true; }
        { samtools view "${bam_path}" "*" 2>/dev/null || true; }
    } | LC_ALL=C awk -v AUDIT_ONLY=1 -f "${AWK_SCRIPT}"
    audit_rc=$?
    set -e
    if [[ ${audit_rc} -ne 0 ]]; then
        echo "ERROR: ${role} BAM failed HTSJDK STRICT flag/mate audit (exit=${audit_rc})!" >&2
        echo "       Broken unpaired flags/mate fields crash CIDER/BAMTOOLS and corrupt ESVEE/TEAL." >&2
        echo "       Fix with single-pass streaming cleanup BEFORE re-running Oncoanalyser:" >&2
        echo "         sbatch -A ${ACCOUNT} -p ${PARTITION} -J clean_${role} ${SOURCE_REPO}/runs/clean_unpaired_flags.sh ${bam_path} ${BAM_DIR}/${sample_lbl}.clean.mc.bam" >&2
        exit 1
    fi
    echo "  [OK] ${role} BAM passed MC + HTSJDK STRICT flag audit."
done

# ==============================================================================
# PIPELINE ARCHITECTURE NOTES ON "STRIPPED ALT & UNMAPPED READS":
# ------------------------------------------------------------------------------
# 1. ESVEE (Structural Variant Calling):
#    - Esvee assembles breakpoints using split reads (SA tags), discordant read pairs,
#      and one-end-anchored (OEA) pairs (mapped read with an unmapped mate).
#    - Because unmapped reads were completely removed from this BAM, unmapped mates (0x8)
#      are absent. This reduces sensitivity for novel non-reference insertions and mobile
#      element assembly.
#    - However, standard deletions, tandem duplications, inversions, and translocations
#      between canonical chromosomes will still be called via discordant pairs and split reads.
#
# 2. VIRUS INTERPRETER:
#    - Relies on viral mappings or unmapped reads. We explicitly exclude it via
#      '--processes_exclude neo,virusinterpreter' to prevent pipeline halts.
#
# 3. AMBER & COBALT:
#    - Query canonical autosomes and sex chromosomes (chr1-22, X, Y) for BAF and depth.
#    - Alt contigs and unmapped reads are natively ignored; BAMs run cleanly.
#
# 4. REDUX RE-ALIGNMENT:
#    - Do NOT run with '--realign'. Extracting FASTQs from BAMs with unmapped reads
#      stripped produces orphaned singletons and drops read pairs.
#    - We reuse existing alignments directly via '--generate-redux-tsvs-only', which
#      computes jitter and microsatellite parameters directly from the BAM.
# ==============================================================================

# ------------------------------------------------------------------------------
# 4. Generate Oncoanalyser v3.1 Samplesheet
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " [Step 3/5] Generating Samplesheet via make_oncoanalyser_samplesheet.py..."
echo "=========================================================="

# Create metadata config TSV mapping BAMs to donor and tissue types
cat <<EOF > "${CONFIG_TSV}"
fastq file name	sample id	tissue type	library type	donor_id	tissue origin
${TUMOUR_BAM}	${SAMPLE_LABEL}	primary tumor	wg	${PROJECT_ID}	bone_marrow
${NORMAL_BAM}	${NORMAL_LABEL}	reference	wg	${PROJECT_ID}	blood
EOF

echo "  Config TSV generated: ${CONFIG_TSV}"
cat "${CONFIG_TSV}"
echo ""

# Run samplesheet generator with --generate-redux-tsvs-only to preserve BAM alignment
python3 "${SOURCE_REPO}/Setup/make_oncoanalyser_samplesheet.py" \
    --config-tsv "${CONFIG_TSV}" \
    --out "${SHEETS_DIR}" \
    --oncoanalyser-outdir "${OUTDIR}" \
    --generate-redux-tsvs-only

SAMPLE_SHEET="${SHEETS_DIR}/${PROJECT_ID}/samplesheet.csv"
if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: Expected samplesheet was not created at: ${SAMPLE_SHEET}" >&2
    exit 1
fi

echo "  [OK] Generated samplesheet: ${SAMPLE_SHEET}"
cat "${SAMPLE_SHEET}"
echo ""

# ------------------------------------------------------------------------------
# 5. Environment & Singularity Setup
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " [Step 4/5] Setting up Execution Environment & Singularity Cache..."
echo "=========================================================="

export SINGULARITY_CACHEDIR="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache"
export NXF_SINGULARITY_CACHEDIR="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache"
export NXF_OFFLINE=true

# Isolate temp and scratch to avoid filling compute node /tmp
export TMPDIR="${TMP_DIR}"
export TEMP="${TMP_DIR}"
export TMP="${TMP_DIR}"
export SINGULARITY_TMPDIR="${TMP_DIR}"
export APPTAINER_TMPDIR="${TMP_DIR}"
export SINGULARITYENV_TMPDIR="${TMP_DIR}"
export SINGULARITYENV_TEMP="${TMP_DIR}"
export SINGULARITYENV_TMP="${TMP_DIR}"

# ------------------------------------------------------------------------------
# 6. SLURM Head Job Submission
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " [Step 5/5] Submitting Oncoanalyser WGTS Run for ${PROJECT_ID}..."
echo "=========================================================="

purple_dir="${PATIENT_OUTDIR}/${SAMPLE_LABEL}/purple"
amber_dir="${PATIENT_OUTDIR}/${SAMPLE_LABEL}/amber"

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "  [DRY-RUN] Preflight and samplesheet generation passed (--preflight-only). Skipping SLURM submission."
    exit 0
fi

# Check if job is currently in queue
running_job=$(squeue -u "$USER" --name="onco_wgts_${PROJECT_ID}" -h -o "%i" | head -n 1 || true)

if [[ -n "$running_job" ]]; then
    echo "  [INFO] Job onco_wgts_${PROJECT_ID} is already queued/running in SLURM (Job ID: ${running_job})."
elif [[ -d "$purple_dir" && -d "$amber_dir" ]]; then
    echo "  [INFO] Outputs already exist at ${purple_dir}."
    echo "         Skipping SLURM submission. To force re-run, rename or remove the patient outdir."
else
    echo "  Submitting SLURM job..."
    job_id=$(sbatch --parsable \
        --partition="${PARTITION}" \
        --account="${ACCOUNT}" \
        --time="${WGTS_TIME}" \
        --job-name="onco_wgts_${PROJECT_ID}" \
        --output="${LOGS_DIR}/onco_wgts_${PROJECT_ID}_%j.out" \
        --error="${LOGS_DIR}/onco_wgts_${PROJECT_ID}_%j.err" \
        --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${TMP_DIR},TEMP=${TMP_DIR},TMP=${TMP_DIR},SINGULARITY_TMPDIR=${TMP_DIR},APPTAINER_TMPDIR=${TMP_DIR},SINGULARITYENV_TMPDIR=${TMP_DIR},SINGULARITYENV_TEMP=${TMP_DIR},SINGULARITYENV_TMP=${TMP_DIR} \
        "${SOURCE_REPO}/run_oncoanalyser.sh" \
        "${SAMPLE_SHEET}" \
        "${PATIENT_OUTDIR}" \
        wgts \
        GRCh38_hmf \
        "${ACCOUNT}" \
        --processes_exclude neo,virusinterpreter)
    job_id=$(echo "$job_id" | cut -d';' -f1)
    echo "  [SUCCESS] Submitted WGTS Pipeline Job ID: ${job_id}"
    echo "  Monitor job status with:  squeue -j ${job_id}"
    echo "  Tail run log with:        tail -f ${LOGS_DIR}/onco_wgts_${PROJECT_ID}_${job_id}.out"
fi

# ------------------------------------------------------------------------------
# 7. Post-Run Summary Compilation Command (for reference)
# ------------------------------------------------------------------------------
# Once the SLURM job completes, compile the multi-omics dashboard:
#
#   python3 ${SOURCE_REPO}/scripts/compile_wgts_summary.py \
#       --results-dir "${OUTDIR}" \
#       --output-dir "${OUTDIR}/wgts_summary" \
#       --cohort-name "SPORE_1004"
# ==============================================================================
