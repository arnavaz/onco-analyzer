#!/bin/bash

## create fastq links
# 1. Create target directory if it doesn't already exist
mkdir -p /cluster/projects/rskimlab/projects/BFIERCE/all_fastq_links
# 2. Loop through the list, resolve paths, and symlink
while IFS= read -r line; do
    # Skip empty lines or comments
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    
    # Replace leading ./ with /cluster/projects/rskimlab/
    source_path="/cluster/projects/rskimlab/${line#./}"
    filename=$(basename "$source_path")
    
    # Create the symlink (force overwrite if link already exists)
    ln -sf "$source_path" "/cluster/projects/rskimlab/projects/BFIERCE/all_fastq_links/$filename"
done < /cluster/projects/rskimlab/fastq_list.txt

module load python3/3.12.11

PREV_OUTDIR="/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser"
OUTDIR="${OUTDIR:-/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3}"
ACCOUNT="deperrotgroup"
PARTITION="${PARTITION:-himem}"

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab

FASTQ_DIR=/cluster/projects/rskimlab/projects/BFIERCE/all_fastq_links
MAPPING_TSV=/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser/metadata/BFIERCE_BFIERCE2_MANUAL.fastqs.map.tsv
SHEETS_DIR="${OUTDIR}/sample_sheets"

echo "Generating sample sheets using existing alignments with generate_redux_tsvs_only..."
python3 "${SOURCE_REPO}/Setup/make_oncoanalyser_samplesheet.py" \
  --fastq-dir "$FASTQ_DIR" \
  --mapping-tsv "$MAPPING_TSV" \
  --use-existing-redux-bams "$PREV_OUTDIR" \
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

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
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
    echo "[${patient}] Submitting Step 1 (WGTS mode with generate_redux_tsvs_only)..."
    wgts_job_id=$(sbatch --parsable \
      --partition="${PARTITION}" \
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
    echo "  Submitted WGTS job ID: ${wgts_job_id}"
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

      # 1. Skip if both WISP and dual (short + standard) ichorCNA are already complete
      has_wisp=$(compgen -G "${plasma_outdir}/*/wisp/*.wisp.summary.tsv" "${plasma_outdir}/*/wisp/*.wisp.tsv" > /dev/null && echo "true" || echo "false")
      has_ichor_short=$(compgen -G "${patient_outdir}/ichorCNA_plasma_only/${p_name}/short/*.params.txt" "${patient_outdir}/ichorCNA_plasma_only/${p_name}/*.short.params.txt" "${patient_outdir}/ichorCNA_plasma_vs_n/${p_name}/short/*.params.txt" > /dev/null && echo "true" || echo "false")
      has_ichor_std=$(compgen -G "${patient_outdir}/ichorCNA_plasma_only/${p_name}/standard/*.params.txt" "${patient_outdir}/ichorCNA_plasma_only/${p_name}/*.standard.params.txt" "${patient_outdir}/ichorCNA_plasma_vs_n/${p_name}/standard/*.params.txt" > /dev/null && echo "true" || echo "false")

      if [[ "${has_wisp}" == "true" && "${has_ichor_short}" == "true" && "${has_ichor_std}" == "true" && "${FORCE_ICHOR:-false}" != "true" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already completed (WISP + both standard/short ichorCNA exist); skipping."
        continue
      fi

      # 2. Check if currently running in SLURM queue
      running_plasma=$(squeue -u "$USER" --name="onco_plasma_${patient}_${p_name}" -h -o "%i" | head -n 1)
      if [[ -n "$running_plasma" ]]; then
        echo "  [${patient}] Plasma (${p_name}) already running in SLURM (Job ID: ${running_plasma}); skipping."
        continue
      fi

      # 3. If WISP is already done, skip Nextflow and run ichorCNA only
      walltime="12:00:00"
      extra_flags=()
      dep_args=()
      if [[ "${has_wisp}" == "true" ]]; then
        echo "  [${patient}] WISP already complete for ${p_name} -> Running ichorCNA post-processing only (--skip-nextflow)."
        extra_flags=(--skip-nextflow)
        if [[ "${FORCE_ICHOR:-false}" == "true" ]]; then
          extra_flags+=(--force-ichor)
        fi
      elif [[ -n "$wgts_job_id" ]]; then
        dep_args=(--dependency=afterok:"${wgts_job_id}")
        echo "  [${patient}] Queuing Step 2 (${p_name}) with dependency on WGTS Job ${wgts_job_id}..."
      else
        echo "  [${patient}] Launching Step 2 (${p_name}) in parallel (WGTS complete)..."
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
        "${extra_flags[@]}"
    done
  else
    echo "  No purity/MRD sample sheet found for ${patient}; skipping Step 2."
  fi
done

## ==============================================================================
## Summarize Cohort WGTS/WISP MRD Results
## ==============================================================================

salloc -c 10 -p pughlab -t 10:0:0 --mem 40G


module load python3/3.10.9

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab
OUTDIR=/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3
python3 "${SOURCE_REPO}/scripts/compile_wgts_summary.py" \
  --results-dir "${OUTDIR}" \
  --cohort-name "BFIERCE" \
  --genes-of-interest PTEN,TP53,PIK3CA,MLH1,MSH2,MSH6,PMS2,POLE,KRAS,ARID1A,CTNNB1,PPP2R1A \
  --fingerprint-stride 20 \
  --fingerprint-threads 10

