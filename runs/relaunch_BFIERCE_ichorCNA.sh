#!/bin/bash
# ==============================================================================
# Dedicated IchorCNA Relauncher for BFIERCE Cohort
# ------------------------------------------------------------------------------
# Generates BOTH:
#   1. Standard 1.0X Downsampled ichorCNA
#   2. Short Fragment (90-150 bp) Size-Selected 1.0X ichorCNA
#
# Directory Output Placement:
#   Root patient directory:
#     <patient>/ichorCNA_plasma_only/<sample_id>/
#       ├── standard/                  (Full standard 1.0x solution)
#       ├── short/                     (Full short 90-150 bp solution)
#       ├── <sample_id>.params.txt     (Points to short solution - default)
#       ├── <sample_id>.short.params.txt
#       ├── <sample_id>.standard.params.txt
#       ├── <sample_id>_genomeWide.pdf (Points to short solution - default)
#       ├── <sample_id>_short_genomeWide.pdf
#       └── <sample_id>_standard_genomeWide.pdf
#     <patient>/ichorCNA/<sample_id>   (Symlink to ichorCNA_plasma_only)
#
# Usage:
#   # Submit all BFIERCE patients:
#   bash runs/relaunch_BFIERCE_ichorCNA.sh
#
#   # Submit a single patient (e.g. BFIERCE_0001):
#   bash runs/relaunch_BFIERCE_ichorCNA.sh BFIERCE_0001
# ==============================================================================
set -euo pipefail

SOURCE_REPO="${SOURCE_REPO:-/cluster/projects/pughlab/src/oncoanalyser_pughlab}"
OUTDIR="${OUTDIR:-/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3}"
SHEETS_DIR="${OUTDIR}/sample_sheets"
ACCOUNT="${ACCOUNT:-pughlab}"
PARTITION="${PARTITION:-pughlab}"
TARGET_PATIENT="${1:-}"

export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

echo "=================================================================="
echo " BFIERCE Dual ichorCNA Relauncher (Standard 1.0X + Short 90-150 bp)"
echo " Results Directory: ${OUTDIR}"
echo " Account / Queue:   ${ACCOUNT} / ${PARTITION}"
echo " Target Patient:    ${TARGET_PATIENT:-All Patients}"
echo " Mode:              Post-processing only (--skip-nextflow, no SAGE)"
echo "=================================================================="

for pdir in "${SHEETS_DIR}"/BFIERCE_*; do
    [[ -d "$pdir" ]] || continue
    patient=$(basename "$pdir")

    if [[ -n "${TARGET_PATIENT}" && "${patient}" != "${TARGET_PATIENT}" ]]; then
        continue
    fi

    patient_outdir="${OUTDIR}/${patient}"
    patient_tmpdir="${patient_outdir}/tmp"
    mkdir -p "${patient_outdir}/logs" "${patient_tmpdir}"

    purity_sheets=()
    for ps in "${pdir}"/samplesheet.purity.*.csv; do
        [[ -f "$ps" ]] && purity_sheets+=("$ps")
    done
    if [[ ${#purity_sheets[@]} -eq 0 && -f "${pdir}/samplesheet.purity.csv" ]]; then
        purity_sheets=("${pdir}/samplesheet.purity.csv")
    fi

    if [[ ${#purity_sheets[@]} -eq 0 ]]; then
        echo "  [${patient}] No purity sample sheet found; skipping."
        continue
    fi

    for purity_sheet in "${purity_sheets[@]}"; do
        [[ -f "$purity_sheet" ]] || continue
        p_name=$(basename "$purity_sheet" .csv)
        p_name="${p_name#samplesheet.purity.}"
        p_name="${p_name#samplesheet.plasma.}"
        p_name="${p_name#samplesheet.}"
        if [[ -z "$p_name" || "$p_name" == "purity" ]]; then
            p_name=$(awk -F',' '$7 ~ /longitudinal_sample/ || $4 ~ /^[pP]lasma/ {print $3}' "$purity_sheet" | head -n 1)
            p_name="${p_name:-plasma}"
        fi

        plasma_outdir="${patient_outdir}/plasma/${p_name}"
        mkdir -p "${plasma_outdir}" "${patient_tmpdir}"

        # Clean up any zero-byte broken files or stray links
        find "${patient_outdir}" -type f -name "*.downsampled_1x.bam" -size 0 -delete 2>/dev/null || true
        # Remove any stray symlink inside plasma/ so it never pollutes the plasma subfolder
        rm -f "${plasma_outdir}/ichorCNA" 2>/dev/null || true

        echo "  Submitting dual ichorCNA job for: ${patient} (${p_name})..."
        job_id=$(sbatch --parsable \
            --partition="${PARTITION}" \
            -A "${ACCOUNT}" \
            --time=12:00:00 \
            --cpus-per-task=8 \
            --mem=32G \
            --job-name="ichor_${patient}_${p_name}" \
            --output="${patient_outdir}/logs/ichor_${patient}_${p_name}_%j.out" \
            --error="${patient_outdir}/logs/ichor_${patient}_${p_name}_%j.err" \
            --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,TMPDIR=${patient_tmpdir},TEMP=${patient_tmpdir},TMP=${patient_tmpdir},SINGULARITY_TMPDIR=${patient_tmpdir},APPTAINER_TMPDIR=${patient_tmpdir},SINGULARITYENV_TMPDIR=${patient_tmpdir},SINGULARITYENV_TEMP=${patient_tmpdir},SINGULARITYENV_TMP=${patient_tmpdir} \
            "${SOURCE_REPO}/run_oncoanalyser.sh" \
            "${purity_sheet}" \
            "${plasma_outdir}" \
            purity_estimate \
            GRCh38_hmf \
            "${ACCOUNT}" \
            --purity_estimate_mode wgts \
            --skip-nextflow \
            --force-ichor)
        job_id=$(echo "$job_id" | cut -d';' -f1)
        echo "    Submitted Job ID: ${job_id}"
    done
done

echo ""
echo "=================================================================="
echo " All ichorCNA post-processing jobs submitted."
echo " Monitor queue with:  squeue -u \$USER --name='ichor_*'"
echo ""
echo " When completed, recompile the cohort summary:"
echo "   python3 ${SOURCE_REPO}/scripts/compile_wgts_summary.py \\"
echo "     --results-dir ${OUTDIR} \\"
echo "     --cohort-name BFIERCE \\"
echo "     --genes-of-interest PTEN,TP53,PIK3CA,MLH1,MSH2,MSH6,PMS2,POLE,KRAS,ARID1A,CTNNB1,PPP2R1A"
echo "=================================================================="

