#!/bin/bash

# ==============================================================================
# NETMRD Workflow Runner: Step 1 (WGTS) + Step 2 (Longitudinal Plasma WISP/MRD)
# ==============================================================================

TARGET_LINK_DIR="/cluster/projects/rskimlab/TEMP/NETMRD/all_fastq_links"
mkdir -p "${TARGET_LINK_DIR}"

# 1. Symlink using the deduplicated list (excludes OICR_data copies that exist elsewhere)
SOURCE_REPO="/cluster/projects/pughlab/src/oncoanalyser_pughlab"
FASTQ_LIST="/cluster/projects/rskimlab/TEMP/NETMRD/NETMRD_fastq_list.dedup.txt"
[[ -f "${FASTQ_LIST}" ]] || FASTQ_LIST="${SOURCE_REPO}/runs/metadata/NETMRD_fastq_list.dedup.txt"

if [[ -f "${FASTQ_LIST}" ]]; then
    while IFS= read -r source_path; do
        [[ -z "$source_path" || "$source_path" =~ ^# ]] && continue
        filename=$(basename "$source_path")
        ln -sf "$source_path" "${TARGET_LINK_DIR}/${filename}"
    done < "${FASTQ_LIST}"
fi

# 2. Configuration & Sample Sheets
MAPPING_TSV="/cluster/projects/rskimlab/TEMP/NETMRD/NETMRD.fastqs.map.tsv"
[[ -f "${MAPPING_TSV}" ]] || MAPPING_TSV="${SOURCE_REPO}/runs/metadata/NETMRD.fastqs.map.tsv"

OUTDIR="${OUTDIR:-/cluster/projects/rskimlab/TEMP/NETMRD/Oncoanalyser}"
SHEETS_DIR="${OUTDIR}/sample_sheets"

# SLURM Account & Partition configuration:
# Test C verified that account 'pughlab' on partition 'pughlab' (node137) is available immediately.
# We cap head-job walltime to 5 days so it also complies with 'himem' (max 7d) and 'all' (max 5d).
ACCOUNT="${ACCOUNT:-pughlab}"
if [[ -z "${PARTITION:-}" ]]; then
    if [[ "${ACCOUNT}" == "pughlab" ]]; then
        PARTITION="pughlab"
    elif [[ "${ACCOUNT}" == "deperrotgroup" ]]; then
        PARTITION="all"
    else
        PARTITION="all"
    fi
fi

# Guard against known invalid account/partition combinations on H4H
if [[ "${ACCOUNT}" == "deperrotgroup" && "${PARTITION}" == "pughlab" ]]; then
    echo "WARNING: Partition 'pughlab' is private to account 'pughlab'. Switching partition to 'all'."
    PARTITION="all"
fi

if [[ -z "${WGTS_TIME:-}" ]]; then
    if [[ "${PARTITION}" == "pughlab" || "${PARTITION}" == "long" ]]; then
        WGTS_TIME="14-00:00:00"
    elif [[ "${PARTITION}" == "himem" ]]; then
        WGTS_TIME="7-00:00:00"
    else
        WGTS_TIME="5-00:00:00"
    fi
fi

module load python3/3.12.11

echo "Generating Oncoanalyser sample sheets for NETMRD..."
python3 "${SOURCE_REPO}/Setup/make_oncoanalyser_samplesheet.py" \
  --fastq-dir "${TARGET_LINK_DIR}" \
  --mapping-tsv "${MAPPING_TSV}" \
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
echo "Submitting NETMRD Oncoanalyser pipeline jobs..."
echo "Account:   ${ACCOUNT}"
echo "Partition: ${PARTITION}"
echo "Outdir:    ${OUTDIR}"
echo "Time Limit:${WGTS_TIME}"
echo "=========================================================="

for sample_sheet in "${SHEETS_DIR}"/*/samplesheet.csv; do
  [[ -f "$sample_sheet" ]] || continue
  patient=$(basename "$(dirname "$sample_sheet")")
  patient_outdir="${OUTDIR}/${patient}"
  patient_tmpdir="${patient_outdir}/tmp"
  
  mkdir -p "${patient_outdir}/logs" "${patient_tmpdir}"
  
  # Determine tumor anchor and check Step 1 (WGTS) status
  tumor_anchor=$(grep -m 1 ",tumor," "$sample_sheet" | cut -d',' -f1 || true)
  purple_dir="${patient_outdir}/${tumor_anchor}/purple"
  amber_dir="${patient_outdir}/${tumor_anchor}/amber"

  # Check if Step 1 (WGTS) is currently running or completed
  running_wgts=$(squeue -u "$USER" --name="onco_wgts_${patient}" -h -o "%i" | head -n 1)

  has_purple="false"
  if compgen -G "${patient_outdir}/*/purple/*.purple.purity.tsv" > /dev/null || [[ -d "$purple_dir" && -d "$amber_dir" ]]; then
    has_purple="true"
  fi

  wgts_job_id=""
  if [[ -n "$running_wgts" ]]; then
    echo "[${patient}] Step 1 (WGTS) is currently running in SLURM queue (Job ID: ${running_wgts})"
    wgts_job_id="$running_wgts"
  elif [[ "${has_purple}" == "true" ]]; then
    echo "[${patient}] Step 1 (WGTS) outputs already exist (${purple_dir}); skipping Step 1 submission."
  else
    echo "[${patient}] Submitting Step 1 (WGTS mode)..."
    sbatch_out=$(sbatch --parsable \
      --partition="${PARTITION}" \
      --time="${WGTS_TIME}" \
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
      --processes_exclude neo,virusinterpreter 2>&1) || true

    if [[ "$sbatch_out" =~ ^[0-9]+ ]]; then
      wgts_job_id=$(echo "$sbatch_out" | cut -d';' -f1)
      echo "  Submitted WGTS job ID: ${wgts_job_id}"
    else
      echo "  ERROR: Failed to submit Step 1 (WGTS) for ${patient}: ${sbatch_out}"
    fi
  fi

  # Step 2: Longitudinal Plasma Purity Estimate (MRD / WISP & ichorCNA)
  purity_sheets=()
  for ps in "${SHEETS_DIR}/${patient}"/samplesheet.purity.*.csv; do
    [[ -f "$ps" ]] && purity_sheets+=("$ps")
  done
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
        p_name=$(awk -F',' '$7 ~ /longitudinal_sample/ || $4 ~ /^[pP]lasma/ {print $3}' "$purity_sheet" | head -n 1)
        p_name="${p_name:-plasma}"
      fi

      plasma_outdir="${patient_outdir}/plasma/${p_name}"
      timepoint_tmpdir="${plasma_outdir}/tmp"
      mkdir -p "${plasma_outdir}" "${patient_outdir}/logs" "${timepoint_tmpdir}"

      # 1. Skip if both WISP and ichorCNA are already complete
      has_wisp=$(compgen -G "${plasma_outdir}/*/wisp/*.wisp.summary.tsv" "${plasma_outdir}/*/wisp/*.wisp.tsv" > /dev/null && echo "true" || echo "false")
      has_ichor=$(compgen -G "${patient_outdir}/ichorCNA_plasma_vs_n/${p_name}/*.params.txt" "${patient_outdir}/ichorCNA_plasma_only/${p_name}/*.params.txt" "${plasma_outdir}/ichorCNA*/*.params.txt" > /dev/null && echo "true" || echo "false")
      if [[ "${has_wisp}" == "true" && "${has_ichor}" == "true" && "${FORCE_ICHOR:-false}" != "true" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already completed (both WISP and ichorCNA exist); skipping."
        continue
      fi

      # 2. Check if currently running in SLURM queue
      running_plasma=$(squeue -u "$USER" --name="onco_plasma_${patient}_${p_name}" -h -o "%i" | head -n 1)
      if [[ -n "$running_plasma" ]]; then
        plasma_reason=$(squeue -j "$running_plasma" -h -o "%r" 2>/dev/null || true)
        if [[ "$plasma_reason" == "DependencyNeverSatisfied" ]]; then
          echo "  [${patient}] Plasma (${p_name}) Job ${running_plasma} is stuck in DependencyNeverSatisfied; cancelling..."
          scancel "$running_plasma" 2>/dev/null || true
          running_plasma=""
        else
          echo "  [${patient}] Plasma (${p_name}) already running in SLURM (Job ID: ${running_plasma}); skipping."
          continue
        fi
      fi

      # 3. Dependency gating: Only launch Step 2 if Step 1 is running or already complete
      walltime="24:00:00"
      extra_flags=()
      dep_args=()
      if [[ "${has_wisp}" == "true" ]]; then
        echo "  [${patient}] WISP already complete for ${p_name} -> Running ichorCNA post-processing only (--skip-nextflow)."
        walltime="12:00:00"
        extra_flags=(--skip-nextflow)
      elif [[ -n "$wgts_job_id" ]]; then
        dep_args=(--dependency=afterok:"${wgts_job_id}")
        echo "  [${patient}] Queuing Step 2 (${p_name}) with dependency on WGTS Job ${wgts_job_id}..."
      elif [[ "${has_purple}" == "true" ]]; then
        echo "  [${patient}] Launching Step 2 (${p_name}) in parallel (WGTS outputs exist)..."
      else
        echo "  [${patient}] SKIPPING Step 2 (${p_name}): Step 1 (WGTS) has not completed and is not running."
        continue
      fi

      sbatch "${dep_args[@]}" \
        --partition="${PARTITION}" \
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
        --processes_exclude neo,virusinterpreter \
        "${extra_flags[@]}"
    done
  else
    echo "  No purity/MRD sample sheet found for ${patient}; skipping Step 2."
  fi
done

## Summarize
salloc -c 10 -p superhimem -t 10:0:0 --mem 100G



module load python3/3.10.9

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab
OUTDIR=/cluster/projects/rskimlab/TEMP/NETMRD/Oncoanalyser
python3 "${SOURCE_REPO}/scripts/compile_wgts_summary.py" \
  --results-dir "${OUTDIR}" \
  --cohort-name "NETMRD" \
  --genes-of-interest MEN1,DAXX,ATRX,TSC2,PTEN,CDKN1B,TP53,RB1,ARID1A,EIF1AX,SDHB,SDHC,SDHD,SDHAF2,RET,VHL,NF1,MAX,TMEM127,FH,APC,KRAS \
  --fingerprint-stride 20 \
  --fingerprint-threads 10

