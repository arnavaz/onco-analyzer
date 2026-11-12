#!/bin/bash
# ==============================================================================
# Relaunch Step 2 (MRD / WISP Plasma Purity Estimate) for all BFIERCE Cases
# Re-runs WISP with -resume to compute LOH (AMBER_LOH), CNV, and SNV purity
# ==============================================================================

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab
OUTDIR="${OUTDIR:-/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3}"
SHEETS_DIR="${OUTDIR}/sample_sheets"
ACCOUNT="deperrotgroup"

export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

echo "=================================================================="
echo " Relaunching Step 2 (Plasma MRD / WISP) for all BFIERCE Cases"
echo " Results Directory: ${OUTDIR}"
echo " Account:           ${ACCOUNT}"
echo "=================================================================="

for patient_dir in "${SHEETS_DIR}"/*; do
    [[ -d "$patient_dir" ]] || continue
    patient=$(basename "$patient_dir")
    patient_outdir="${OUTDIR}/${patient}"
    patient_tmpdir="${patient_outdir}/tmp"
    
    purity_sheets=()
    for ps in "${patient_dir}"/samplesheet.purity.*.csv; do
        [[ -f "$ps" ]] && purity_sheets+=("$ps")
    done
    if [[ ${#purity_sheets[@]} -eq 0 && -f "${patient_dir}/samplesheet.purity.csv" ]]; then
        purity_sheets=("${patient_dir}/samplesheet.purity.csv")
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
        mkdir -p "${plasma_outdir}" "${patient_outdir}/logs" "${patient_tmpdir}"
        
        echo "Submitting Step 2 (Plasma / WISP) for: ${patient} (${p_name})"
        sbatch \
          --partition=all \
          --time=24:00:00 \
          --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,FORCE_NEXTFLOW=true,TMPDIR=${patient_tmpdir},TEMP=${patient_tmpdir},TMP=${patient_tmpdir},SINGULARITY_TMPDIR=${patient_tmpdir},APPTAINER_TMPDIR=${patient_tmpdir},SINGULARITYENV_TMPDIR=${patient_tmpdir},SINGULARITYENV_TEMP=${patient_tmpdir},SINGULARITYENV_TMP=${patient_tmpdir} \
          -A "$ACCOUNT" \
          --job-name="onco_plasma_${patient}_${p_name}" \
          --output="${patient_outdir}/logs/onco_plasma_${patient}_${p_name}_%j.out" \
          --error="${patient_outdir}/logs/onco_plasma_${patient}_${p_name}_%j.err" \
          "${SOURCE_REPO}/run_oncoanalyser.sh" \
          "$purity_sheet" \
          "$plasma_outdir" \
          purity_estimate \
          GRCh38_hmf \
          "$ACCOUNT" \
          --purity_estimate_mode wgts \
          --force-nextflow
    done
done

echo "=================================================================="
echo " All Step 2 plasma jobs queued with -resume."
echo " Monitor with: squeue -u \$USER --name='onco_plasma*'"
echo "=================================================================="
