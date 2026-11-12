#!/bin/bash

# ==============================================================================
# TRACE_OC MRD / Purity Estimation Runner Script (Ultima Genomics BAMs)
# 
# Pipeline Flow:
#   Step 1: WGTS Mode (Primary Tumor + Normal REDUX BAMs) -> generates PURPLE & AMBER
#   Step 2: Purity Estimate Mode (Longitudinal Plasma REDUX BAMs + Normal REDUX BAM
#           + primary PURPLE/AMBER dirs) -> runs SAGE_APPEND, AMBER, COBALT, and WISP.
# ==============================================================================

module load python3/3.12.11

PREV_OUTDIR="/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser"
OUTDIR="${OUTDIR:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3}"
ACCOUNT="${ACCOUNT:-pughlab}"
PARTITION="${PARTITION:-pughlab}"
SOURCE_REPO="/cluster/projects/pughlab/src/oncoanalyser_pughlab"

# Sample mapping TSV containing REDUX BAM paths and annotations
CONFIG_TSV="${1:-"${SOURCE_REPO}/runs/metadata/redux_config_TRACEOC.tsv"}"
SHEETS_DIR="${OUTDIR}/sample_sheets"

mkdir -p "${SHEETS_DIR}"
mkdir -p "${OUTDIR}"

echo "=========================================================="
echo "Generating Oncoanalyser sample sheets for TRACE_OC..."
echo "Config TSV: ${CONFIG_TSV}"
echo "Sheets Dir: ${SHEETS_DIR}"
echo "Outdir:     ${OUTDIR}"
echo "=========================================================="

python3 "${SOURCE_REPO}/Setup/make_oncoanalyser_samplesheet.py" \
  --config-tsv "$CONFIG_TSV" \
  --out "$SHEETS_DIR" \
  --oncoanalyser-outdir "$OUTDIR"

export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

# Setup user-controlled temporary directory to prevent compute node /tmp exhaustion
export TMPDIR="${OUTDIR}/tmp"
export TEMP="${OUTDIR}/tmp"
export TMP="${OUTDIR}/tmp"
export SINGULARITY_TMPDIR="${OUTDIR}/tmp"
export APPTAINER_TMPDIR="${OUTDIR}/tmp"
export SINGULARITYENV_TMPDIR="${OUTDIR}/tmp"
export SINGULARITYENV_TEMP="${OUTDIR}/tmp"
export SINGULARITYENV_TMP="${OUTDIR}/tmp"
mkdir -p "${TMPDIR}"

echo ""
echo "=========================================================="
echo "Submitting TRACE_OC Oncoanalyser pipeline jobs..."
echo "=========================================================="

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
  [[ -f "$sample_sheet" ]] || continue
  patient=$(basename "$(dirname "$sample_sheet")")
  patient_outdir="${OUTDIR}/${patient}"
  patient_tmpdir="${patient_outdir}/tmp"

  mkdir -p "${patient_outdir}/logs" "${patient_tmpdir}"

  # Determine tumor anchor
  tumor_anchor=$(grep -m 1 ",tumor," "$sample_sheet" | cut -d',' -f1)
  purple_dir="${patient_outdir}/${tumor_anchor}/purple"
  amber_dir="${patient_outdir}/${tumor_anchor}/amber"

  # Check if Step 1 (WGTS) is currently running or completed
  running_wgts=$(squeue -u "$USER" --name="onco_wgts_${patient}" -h -o "%i" | head -n 1)

  wgts_job_id=""
  if [[ -n "$running_wgts" ]]; then
    echo "[${patient}] Step 1 (WGTS) is currently running in SLURM queue (Job ID: ${running_wgts})"
    wgts_job_id="$running_wgts"
  elif [[ -d "$purple_dir" && -d "$amber_dir" ]]; then
    echo "[${patient}] Step 1 (WGTS) outputs already exist (${purple_dir}); skipping Step 1 submission."
  else
    echo "[${patient}] Submitting Step 1 (WGTS mode with generate_redux_tsvs_only)..."
    wgts_job_id=$(sbatch --parsable \
      --partition="$PARTITION" \
      --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${patient_tmpdir},TEMP=${patient_tmpdir},TMP=${patient_tmpdir},SINGULARITY_TMPDIR=${patient_tmpdir},APPTAINER_TMPDIR=${patient_tmpdir},SINGULARITYENV_TMPDIR=${patient_tmpdir},SINGULARITYENV_TEMP=${patient_tmpdir},SINGULARITYENV_TMP=${patient_tmpdir} \
      -A "$ACCOUNT" \
      --job-name="onco_wgts_${patient}" \
      --output="${patient_outdir}/logs/onco_wgts_${patient}_%j.out" \
      --error="${patient_outdir}/logs/onco_wgts_${patient}_%j.err" \
      "${SOURCE_REPO}/run_oncoanalyser.sh" \
      "$sample_sheet" \
      "$patient_outdir" \
      wgts \
      GRCh38_hmf \
      "$ACCOUNT" \
      --processes_exclude neo,virusinterpreter)
    wgts_job_id=$(echo "$wgts_job_id" | cut -d';' -f1)
    echo "  Submitted WGTS Job ID: ${wgts_job_id}"
  fi

  # Step 2: Longitudinal Plasma Purity Estimate (MRD / WISP)
  # Look for timepoint-specific purity sheets (e.g. samplesheet.purity.OV-P-02-T0.csv)
  purity_sheets=()
  for ps in "${SHEETS_DIR}/${patient}"/samplesheet.purity.*.csv; do
    [[ -f "$ps" ]] && purity_sheets+=("$ps")
  done

  # Fallback to single samplesheet.purity.csv if no timepoint-specific sheets found
  if [[ ${#purity_sheets[@]} -eq 0 && -f "${SHEETS_DIR}/${patient}/samplesheet.purity.csv" ]]; then
    purity_sheets=("${SHEETS_DIR}/${patient}/samplesheet.purity.csv")
  fi

  if [[ ${#purity_sheets[@]} -gt 0 ]]; then
    for purity_sheet in "${purity_sheets[@]}"; do
      [[ -f "$purity_sheet" ]] || continue
      p_name=$(basename "$purity_sheet" .csv)
      p_name="${p_name#samplesheet.purity.}"
      p_name="${p_name#samplesheet.plasma.}"
      p_name="${p_name#samplesheet.}"
      if [[ -z "$p_name" || "$p_name" == "purity" ]]; then
        p_name=$(awk -F',' '$3=="plasma" || $3=="Plasma" {print $2}' "$purity_sheet" | head -n 1)
        p_name="${p_name:-plasma}"
      fi

      plasma_outdir="${patient_outdir}/plasma/${p_name}"
      timepoint_tmpdir="${plasma_outdir}/tmp"
      mkdir -p "${plasma_outdir}" "${timepoint_tmpdir}"

      # 1. Skip if both WISP and ichorCNA are already complete
      has_wisp=$(compgen -G "${plasma_outdir}/*/wisp/*.wisp.summary.tsv" "${plasma_outdir}/*/wisp/*.summary.tsv" "${plasma_outdir}/*/wisp/*.tsv" > /dev/null && echo "true" || echo "false")
      has_ichor=$(compgen -G "${patient_outdir}/ichorCNA_plasma_vs_n/${p_name}/*.params.txt" "${patient_outdir}/ichorCNA_plasma_only/${p_name}/*.params.txt" "${plasma_outdir}/ichorCNA*/*.params.txt" > /dev/null && echo "true" || echo "false")
      if [[ "${has_wisp}" == "true" && "${has_ichor}" == "true" && "${FORCE_ICHOR:-false}" != "true" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already completed (both WISP and ichorCNA exist); skipping."
        continue
      fi

      # 2. Check if currently running in SLURM queue
      running_plasma=$(squeue -u "$USER" --name="onco_plasma_${patient}_${p_name}" -h -o "%i" | head -n 1)
      if [[ -n "$running_plasma" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already running in SLURM (Job ID: ${running_plasma}); skipping."
        continue
      fi

      # 3. If WISP is already done, skip Nextflow and run ichorCNA only
      walltime="24:00:00"
      extra_flags=()
      dep_args=()
      if [[ "${has_wisp}" == "true" ]]; then
        echo "  [${patient}] WISP already complete for ${p_name} -> Running ichorCNA post-processing only (--skip-nextflow)."
        walltime="02:00:00"
        extra_flags=(--skip-nextflow)
      elif [[ -n "$wgts_job_id" ]]; then
        dep_args=(--dependency=afterok:"${wgts_job_id}")
        echo "  [${patient}] Queuing Step 2 (${p_name}) with dependency on WGTS Job ${wgts_job_id}..."
      else
        echo "  [${patient}] Launching Step 2 (${p_name}) in parallel (WGTS complete)..."
      fi

      purity_job_id=$(sbatch --parsable "${dep_args[@]}" \
        --partition="$PARTITION" \
        --time="${walltime}" \
        --cpus-per-task=8 \
        --mem=32G \
        --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${timepoint_tmpdir},TEMP=${timepoint_tmpdir},TMP=${timepoint_tmpdir},SINGULARITY_TMPDIR=${timepoint_tmpdir},APPTAINER_TMPDIR=${timepoint_tmpdir},SINGULARITYENV_TMPDIR=${timepoint_tmpdir},SINGULARITYENV_TEMP=${timepoint_tmpdir},SINGULARITYENV_TMP=${timepoint_tmpdir} \
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
        --sequencing_platform ultima \
        "${extra_flags[@]}")
      purity_job_id=$(echo "$purity_job_id" | cut -d';' -f1)
      echo "    Submitted Plasma Job ID: ${purity_job_id}"
    done
  else
    echo "  [${patient}] No purity/MRD sample sheets found; skipping Step 2."
  fi
done

echo ""
echo "=========================================================="
echo "All jobs queued. Monitor with: squeue -u $USER"
echo "=========================================================="

# ==============================================================================
# Standalone: Rerun Step 2 (Purity only) when Step 1 is already complete / active
# ==============================================================================
# cd "${SOURCE_REPO}"
# OUTDIR="${OUTDIR:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3}"
# ACCOUNT="${ACCOUNT:-pughlab}"
# PARTITION="${PARTITION:-pughlab}"
# SHEETS_DIR="${OUTDIR}/sample_sheets"
# for purity_sheet in "${SHEETS_DIR}"/*/samplesheet.purity.*.csv; do
#     [[ -f "$purity_sheet" ]] || continue
#     patient=$(basename "$(dirname "$purity_sheet")")
#     patient_outdir="${OUTDIR}/${patient}"
#     p_name=$(basename "$purity_sheet" .csv)
#     p_name="${p_name#samplesheet.purity.}"
#     p_name="${p_name#samplesheet.plasma.}"
#     p_name="${p_name#samplesheet.}"
#     if [[ -z "$p_name" || "$p_name" == "purity" ]]; then
#         p_name=$(awk -F',' '$3=="plasma" || $3=="Plasma" {print $2}' "$purity_sheet" | head -n 1)
#         p_name="${p_name:-plasma}"
#     fi
#     plasma_outdir="${patient_outdir}/plasma/${p_name}"
#     mkdir -p "${plasma_outdir}" "${patient_outdir}/logs"
#     
#     # Check if Step 1 (WGTS) is still active in the SLURM queue
#     wgts_job_id=$(squeue -u "$USER" --name="onco_wgts_${patient}" -h -o "%i" | head -n 1)
#     
#     if [[ -n "$wgts_job_id" ]]; then
#         dep_arg="--dependency=afterok:${wgts_job_id}"
#         echo "[${patient}] Step 1 is still running (Job ID: ${wgts_job_id}) -> Queuing Step 2 with dependency."
#     else
#         dep_arg=""
#         echo "[${patient}] Step 1 is already finished -> Launching Step 2 immediately."
#     fi
#     sbatch ${dep_arg} \
#       --partition="${PARTITION}" \
#       --time=08:00:00 \
#       --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true \
#       -A "$ACCOUNT" \
#       --job-name="onco_plasma_${patient}_${p_name}" \
#       --output="${patient_outdir}/logs/onco_plasma_${patient}_${p_name}_%j.out" \
#       --error="${patient_outdir}/logs/onco_plasma_${patient}_${p_name}_%j.err" \
#       "${SOURCE_REPO}/run_oncoanalyser.sh" \
#       "$purity_sheet" \
#       "$plasma_outdir" \
#       purity_estimate \
#       GRCh38_hmf \
#       "$ACCOUNT" \
#       --purity_estimate_mode wgts \
#       --sequencing_platform ultima
# done

## ==============================================================================
## Summarize Cohort WGTS/WISP MRD Results
## ==============================================================================
# module load python3/3.12.11
# python3 "${SOURCE_REPO}/scripts/compile_wgts_summary.py" \
#   --results-dir "${OUTDIR}" \
#   --cohort-name "TRACE_OC"


## using Tobi's results
module load python3/3.10.9

SOURCE_REPO="/cluster/projects/pughlab/src/oncoanalyser_pughlab"

python3 "${SOURCE_REPO}/scripts/compile_wgts_summary.py" \
  --results-dir "/cluster/projects/pughlab/myeloma/projects/MM_cell_drugs/WGS_Pipeline/HARTWIG/oncoanalyser_output/" \
  --output-dir "/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Tobi_summary" \
  --cohort-name "TRACEOC" \
   --timepoints-tsv "/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Tobi_summary/TRACEOC_metadata.tsv" \
  --genes-of-interest BRCA1,BRCA2,TP53,ATM,CHEK2,KMT2C,SLX4,KMT2A,CTCF,SMC3,TET2,ZFHX3,CHD4,ATRX,MED12,MRE11A,MUC1,SMARCA4,PIK3CA,AMER1,PTPRS,BRD4,IRF2,MET,MTOR,PIK3C2G,PIK3CB,PIK3CD,PIK3R1,TCF7L2,PTCH1,NUP93,CYP1B1,TAP2

# ==============================================================================
# Run TRACEOC using Tobi's CRAMs with Oncoanalyser v3.1
# Cram list: /cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3/sample_info/TRACEOC_crams_samples.tsv
# ==============================================================================

module load python3/3.12.11

OUTDIR="${OUTDIR:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3}"
ACCOUNT="${ACCOUNT:-pughlab}"
PARTITION="${PARTITION:-pughlab}"
SOURCE_REPO="${SOURCE_REPO:-/cluster/projects/pughlab/src/oncoanalyser_pughlab}"
CRAM_CONFIG_TSV="/cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3/sample_info/TRACEOC_crams_samples.tsv"
SHEETS_DIR="${OUTDIR}/sample_sheets"

mkdir -p "${SHEETS_DIR}" "${OUTDIR}"

echo "=========================================================="
echo "Generating Oncoanalyser v3.1 sample sheets from CRAMs..."
echo "Config TSV: ${CRAM_CONFIG_TSV}"
echo "Sheets Dir: ${SHEETS_DIR}"
echo "Outdir:     ${OUTDIR}"
echo "=========================================================="

python3 "${SOURCE_REPO}/Setup/make_oncoanalyser_samplesheet.py" \
  --config-tsv "${CRAM_CONFIG_TSV}" \
  --out "${SHEETS_DIR}" \
  --oncoanalyser-outdir "${OUTDIR}"

export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

# Setup user-controlled temporary directory to prevent compute node /tmp exhaustion
export TMPDIR="${OUTDIR}/tmp"
export TEMP="${OUTDIR}/tmp"
export TMP="${OUTDIR}/tmp"
export SINGULARITY_TMPDIR="${OUTDIR}/tmp"
export APPTAINER_TMPDIR="${OUTDIR}/tmp"
export SINGULARITYENV_TMPDIR="${OUTDIR}/tmp"
export SINGULARITYENV_TEMP="${OUTDIR}/tmp"
export SINGULARITYENV_TMP="${OUTDIR}/tmp"
mkdir -p "${TMPDIR}"

echo ""
echo "=========================================================="
echo "Submitting TRACEOC v3.1 CRAM pipeline jobs..."
echo "=========================================================="

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
  [[ -f "$sample_sheet" ]] || continue
  patient=$(basename "$(dirname "$sample_sheet")")
  patient_outdir="${OUTDIR}/${patient}"
  patient_tmpdir="${patient_outdir}/tmp"

  mkdir -p "${patient_outdir}/logs" "${patient_tmpdir}"

  # Determine tumor anchor
  tumor_anchor=$(grep -m 1 ",tumor," "$sample_sheet" | cut -d',' -f1)
  purple_dir="${patient_outdir}/${tumor_anchor}/purple"
  amber_dir="${patient_outdir}/${tumor_anchor}/amber"

  # Check if Step 1 (WGTS) is currently running or completed
  running_wgts=$(squeue -u "$USER" --name="onco_wgts_${patient}" -h -o "%i" | head -n 1)

  wgts_job_id=""
  if [[ -n "$running_wgts" ]]; then
    echo "[${patient}] Step 1 (WGTS) is currently running in SLURM queue (Job ID: ${running_wgts})"
    wgts_job_id="$running_wgts"
  elif [[ -d "$purple_dir" && -d "$amber_dir" ]]; then
    echo "[${patient}] Step 1 (WGTS) outputs already exist (${purple_dir}); skipping Step 1 submission."
  else
    echo "[${patient}] Submitting Step 1 (WGTS mode with Ultima platform)..."
    wgts_job_id=$(sbatch --parsable \
      --partition="$PARTITION" \
      --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${patient_tmpdir},TEMP=${patient_tmpdir},TMP=${patient_tmpdir},SINGULARITY_TMPDIR=${patient_tmpdir},APPTAINER_TMPDIR=${patient_tmpdir},SINGULARITYENV_TMPDIR=${patient_tmpdir},SINGULARITYENV_TEMP=${patient_tmpdir},SINGULARITYENV_TMP=${patient_tmpdir} \
      -A "$ACCOUNT" \
      --job-name="onco_wgts_${patient}" \
      --output="${patient_outdir}/logs/onco_wgts_${patient}_%j.out" \
      --error="${patient_outdir}/logs/onco_wgts_${patient}_%j.err" \
      "${SOURCE_REPO}/run_oncoanalyser.sh" \
      "$sample_sheet" \
      "$patient_outdir" \
      wgts \
      GRCh38_hmf \
      "$ACCOUNT" \
      --processes_exclude neo,virusinterpreter \
      --sequencing_platform ultima)
    wgts_job_id=$(echo "$wgts_job_id" | cut -d';' -f1)
    echo "  Submitted WGTS Job ID: ${wgts_job_id}"
  fi

  # Step 2: Longitudinal Plasma Purity Estimate (MRD / WISP)
  purity_sheets=()
  for ps in "${SHEETS_DIR}/${patient}"/samplesheet.purity.*.csv; do
    [[ -f "$ps" ]] && purity_sheets+=("$ps")
  done

  if [[ ${#purity_sheets[@]} -eq 0 && -f "${SHEETS_DIR}/${patient}/samplesheet.purity.csv" ]]; then
    purity_sheets=("${SHEETS_DIR}/${patient}/samplesheet.purity.csv")
  fi

  if [[ ${#purity_sheets[@]} -gt 0 ]]; then
    prev_dep="${wgts_job_id}"
    for purity_sheet in "${purity_sheets[@]}"; do
      [[ -f "$purity_sheet" ]] || continue
      p_name=$(basename "$purity_sheet" .csv)
      p_name="${p_name#samplesheet.purity.}"
      p_name="${p_name#samplesheet.plasma.}"
      p_name="${p_name#samplesheet.}"
      if [[ -z "$p_name" || "$p_name" == "purity" ]]; then
        p_name=$(awk -F',' '$3=="plasma" || $3=="Plasma" {print $2}' "$purity_sheet" | head -n 1)
        p_name="${p_name:-plasma}"
      fi

      plasma_outdir="${patient_outdir}/plasma/${p_name}"

      # 1. Skip if both WISP and ichorCNA are already complete
      has_wisp=$(compgen -G "${plasma_outdir}/*/wisp/*.wisp.summary.tsv" > /dev/null && echo "true" || echo "false")
      has_ichor=$(compgen -G "${patient_outdir}/ichorCNA_plasma_vs_n/${p_name}/*.params.txt" "${patient_outdir}/ichorCNA_plasma_only/${p_name}/*.params.txt" > /dev/null && echo "true" || echo "false")
      if [[ "${has_wisp}" == "true" && "${has_ichor}" == "true" && "${FORCE_ICHOR:-false}" != "true" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already completed (both WISP and ichorCNA exist); skipping."
        continue
      fi

      # 2. Check if currently running in SLURM queue
      running_plasma=$(squeue -u "$USER" --name="onco_plasma_${patient}_${p_name}" -h -o "%i" | head -n 1)
      if [[ -n "$running_plasma" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already running in SLURM (Job ID: ${running_plasma}); chaining."
        prev_dep="${running_plasma}"
        continue
      fi

      walltime="12:00:00"
      extra_flags=()
      if [[ "${has_wisp}" == "true" ]]; then
        echo "  [${patient}] WISP already complete for ${p_name} -> Running ichorCNA post-processing only (--skip-nextflow)."
        walltime="04:00:00"
        extra_flags=(--skip-nextflow)
      fi

      dep_args=()
      if [[ -n "$prev_dep" ]]; then
        dep_args=(--dependency=afterok:"${prev_dep}")
        echo "  [${patient}] Queuing Step 2 (${p_name}) with dependency on Job ${prev_dep}..."
      else
        echo "  [${patient}] Launching Step 2 (${p_name}) immediately (WGTS already complete)..."
      fi

      mkdir -p "${plasma_outdir}"

      purity_job_id=$(sbatch --parsable "${dep_args[@]}" \
        --partition="$PARTITION" \
        --time="${walltime}" \
        --cpus-per-task=8 \
        --mem=32G \
        --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${patient_tmpdir},TEMP=${patient_tmpdir},TMP=${patient_tmpdir},SINGULARITY_TMPDIR=${patient_tmpdir},APPTAINER_TMPDIR=${patient_tmpdir},SINGULARITYENV_TMPDIR=${patient_tmpdir},SINGULARITYENV_TEMP=${patient_tmpdir},SINGULARITYENV_TMP=${patient_tmpdir} \
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
        --sequencing_platform ultima \
        "${extra_flags[@]}")
      purity_job_id=$(echo "$purity_job_id" | cut -d';' -f1)
      echo "    Submitted Plasma Job ID: ${purity_job_id}"

      # Chain consecutive timepoint runs to avoid concurrent Nextflow directory lock collisions
      prev_dep="${purity_job_id}"
    done
  else
    echo "  [${patient}] No purity/MRD sample sheets found; skipping Step 2."
  fi
done

echo ""
echo "=========================================================="
echo "All TRACEOC v3.1 CRAM jobs queued. Monitor with: squeue -u $USER"
echo "=========================================================="


## ==============================================================================
## Summarize Cohort WGTS/WISP MRD Results (Run after jobs complete)
## ==============================================================================
salloc -c 10 -p pughlab -t 10:0:0 --mem 40G

module load python3/3.10.9
python3 /cluster/projects/pughlab/src/oncoanalyser_pughlab/scripts/compile_wgts_summary.py \
  --results-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3 \
  --output-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Oncoanalyser_v3/summary \
  --cohort-name "TRACEOC_v3.1" \
  --timepoints-tsv /cluster/projects/mohccn/Sandbox/JEFF_TMP/TRACEOC/Tobi_summary/TRACEOC_metadata.tsv \
  --genes-of-interest BRCA1,BRCA2,TP53,ATM,CHEK2,KMT2C,SLX4,KMT2A,CTCF,SMC3,TET2,ZFHX3,CHD4,ATRX,MED12,MRE11A,MUC1,SMARCA4,PIK3CA,AMER1,PTPRS,BRD4,IRF2,MET,MTOR,PIK3C2G,PIK3CB,PIK3CD,PIK3R1,TCF7L2,PTCH1,NUP93,CYP1B1,TAP2 \
  --fingerprint-stride 20 \
  --fingerprint-threads 10
## Head to head to compare v3.0 vs. v3.1