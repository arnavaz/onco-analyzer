#!/bin/bash
# ==============================================================================
# run_BFIERCE_delfi_tf.sh
# ------------------------------------------------------------------------------
# Dedicated DELFI-TF Feature Extractor & Dispatcher for BFIERCE Cohort.
#
# Supports older and newer BFIERCE directory layouts:
#   - V2: /cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser/<patient>/...
#   - V3: /cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3/<patient>/...
#   - Links: /cluster/projects/rskimlab/projects/BFIERCE/fragmentomics/bam_links/<patient>/...
#
# Features:
#   - Automatically reuses 1.0X downsampled BAMs from ichorCNA if available (<60s runtime)
#   - Uses chunked streaming with --max-reads (default: 50M reads, ~2X) on deep WGS BAMs (no OOM)
#   - Isolates all logs and results under <patient>/DELFI/<sample_id>/logs/
#   - Built-in --compile mode to aggregate features and matched-normal truth for model training
#
# Usage:
#   # Submit all BFIERCE patients:
#   bash runs/run_BFIERCE_delfi_tf.sh
#
#   # Submit a single patient (e.g. BFIERCE_0003):
#   bash runs/run_BFIERCE_delfi_tf.sh BFIERCE_0003
#
#   # Compile extracted features into a DELFI-TF training manifest:
#   bash runs/run_BFIERCE_delfi_tf.sh --compile
# ==============================================================================
set -euo pipefail

TARGET_PATIENT="${1:-}"
ACCOUNT="${2:-pughlab}"
PARTITION="${3:-pughlab}"
MAX_READS="${4:-50000000}"

DELFI_REPO="${DELFI_REPO:-/cluster/projects/pughlab/src/DELFI-TF}"
V2_BASE="${V2_BASE:-/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser}"
V3_BASE="${V3_BASE:-/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3}"
BAM_LINKS="${BAM_LINKS:-/cluster/projects/rskimlab/projects/BFIERCE/fragmentomics/bam_links}"
SHEETS_DIR="${V3_BASE}/sample_sheets"
MODEL_PATH="${DELFI_REPO}/models/delfi_tf_default_model.rds"

# --- Compilation Mode ---
if [[ "${TARGET_PATIENT}" == "--compile" || "${TARGET_PATIENT}" == "-c" ]]; then
    echo "=================================================================="
    echo " Compiling BFIERCE DELFI-TF Training Manifest"
    echo " Source Directory: ${V3_BASE}"
    echo "=================================================================="

    OUT_MANIFEST="${V3_BASE}/DELFI/BFIERCE_delfi_training_manifest.csv"
    mkdir -p "$(dirname "${OUT_MANIFEST}")"

    python3 - <<'EOF'
import os, glob, re, csv, sys

v3_base = "/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3"
out_manifest = os.path.join(v3_base, "DELFI", "BFIERCE_delfi_training_manifest.csv")

# Known high-confidence matched-normal ichorCNA and PURPLE truth benchmark
KNOWN_TRUTH = {
    "536878": 0.0855,  # BFIERCE_0003
    "536882": 0.1340,  # BFIERCE_0007
    "536881": 0.1036,  # BFIERCE_0006
    "536884": 0.1030,  # BFIERCE_0009
    "536880": 0.0350,  # BFIERCE_0005
    "536879": 0.0400,  # BFIERCE_0004
}

rows = []
feature_files = glob.glob(os.path.join(v3_base, "*", "DELFI", "*", "*_delfi_tf_features.tsv"))
print(f"Found {len(feature_files)} extracted DELFI-TF feature files.")

for ff in sorted(feature_files):
    sample_id = os.path.basename(ff).replace("_delfi_tf_features.tsv", "")
    sample_dir = os.path.dirname(ff)
    patient_dir = os.path.dirname(os.path.dirname(sample_dir))
    patient_id = os.path.basename(patient_dir)

    # 1. Resolve BAM path
    bam_path = ""
    # Check for downsampled BAM first
    for cand in [
        os.path.join(patient_dir, f"ichorCNA_plasma_only/{sample_id}/{sample_id}.downsampled_1x.bam"),
        os.path.join(patient_dir, f"ichorCNA_plasma_vs_n/{sample_id}/{sample_id}.downsampled_1x.bam"),
        os.path.join(f"/cluster/projects/rskimlab/projects/BFIERCE/fragmentomics/bam_links/{patient_id}/plasma/{sample_id}.bam"),
    ]:
        if os.path.isfile(cand):
            bam_path = cand
            break

    if not bam_path:
        # Glob for full redux BAM
        matches = glob.glob(os.path.join(patient_dir, f"plasma/*/*/alignments/{sample_id}/{sample_id}.redux.bam"))
        if not matches:
            matches = glob.glob(os.path.join(f"/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser/{patient_id}/*/alignments/dna/{sample_id}.redux.bam"))
        if matches:
            bam_path = matches[0]

    # 2. Resolve Ground Truth TF
    gt_tf = None
    if sample_id in KNOWN_TRUTH:
        gt_tf = KNOWN_TRUTH[sample_id]
    else:
        # Check matched normal ichorCNA params
        param_candidates = glob.glob(os.path.join(patient_dir, f"ichorCNA_plasma_vs_n/{sample_id}/*params.txt"))
        if not param_candidates:
            param_candidates = glob.glob(os.path.join(patient_dir, f"ichorCNA_plasma_only/{sample_id}/*params.txt"))
        for pf in param_candidates:
            try:
                with open(pf, "r") as f:
                    for line in f:
                        if "tumor fraction" in line.lower():
                            parts = line.split(":")
                            if len(parts) >= 2:
                                gt_tf = float(parts[1].strip())
                                break
            except Exception:
                pass
            if gt_tf is not None:
                break

    if gt_tf is None:
        # Check PURPLE or WISP
        purple_files = glob.glob(os.path.join(patient_dir, f"plasma/*/{sample_id}*.purple.purity.tsv"))
        for pf in purple_files:
            try:
                with open(pf, "r") as f:
                    lines = [l.strip() for l in f if l.strip()]
                    if len(lines) >= 2:
                        gt_tf = float(lines[1].split()[0])
            except Exception:
                pass

    if gt_tf is not None and bam_path:
        rows.append({
            "sample_id": sample_id,
            "bam_path": bam_path,
            "ground_truth_tf": f"{gt_tf:.4f}",
            "cohort": "BFIERCE"
        })

if rows:
    with open(out_manifest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "bam_path", "ground_truth_tf", "cohort"])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print(f"Successfully compiled {len(rows)} training samples to: {out_manifest}")
    print("\nSamples compiled:")
    for r in rows:
        print(f"  - {r['sample_id']} (Truth: {r['ground_truth_tf']}) -> {r['bam_path']}")
else:
    print("Warning: No samples with both features and ground truth found yet.")
EOF
    exit 0
fi

echo "=================================================================="
echo " BFIERCE DELFI-TF Batch Dispatcher"
echo " Account / Partition: ${ACCOUNT} / ${PARTITION}"
echo " DELFI Repo:          ${DELFI_REPO}"
echo " Target Patient:      ${TARGET_PATIENT:-All Patients}"
echo " Max Fragments:       ${MAX_READS} (0 = Full Deep WGS)"
echo "=================================================================="

# Function to locate plasma BAM across V2, V3, and bam_links
resolve_bfierce_plasma() {
    local patient="$1"
    local sample_id="$2"

    # Priority 1: Downsampled 1.0X BAM from previous ichorCNA run (instant runtime)
    for cand in \
        "${V3_BASE}/${patient}/ichorCNA_plasma_only/${sample_id}/${sample_id}.downsampled_1x.bam" \
        "${V3_BASE}/${patient}/ichorCNA_plasma_vs_n/${sample_id}/${sample_id}.downsampled_1x.bam" \
        "${V3_BASE}/${patient}/ichorCNA/${sample_id}/${sample_id}.downsampled_1x.bam"; do
        if [[ -f "$cand" && -s "$cand" ]]; then
            echo "$cand"
            return 0
        fi
    done

    # Priority 2: Direct symlink in fragmentomics bam_links
    for cand in "${BAM_LINKS}/${patient}/plasma/${sample_id}.bam" \
                "${BAM_LINKS}/${patient}/plasma/"*.bam; do
        if [[ -f "$cand" && -s "$cand" ]]; then
            echo "$cand"
            return 0
        fi
    done

    # Priority 3: V3 redux BAM
    for cand in "${V3_BASE}/${patient}"/plasma/*/*/alignments/"${sample_id}"/"${sample_id}".redux.bam \
                "${V3_BASE}/${patient}"/*/alignments/"${sample_id}"/"${sample_id}".redux.bam; do
        if [[ -f "$cand" && -s "$cand" ]]; then
            echo "$cand"
            return 0
        fi
    done

    # Priority 4: V2 redux BAM
    for cand in "${V2_BASE}/${patient}"/*/alignments/dna/"${sample_id}".redux.bam \
                "${V2_BASE}/${patient}"/purity/samplesheet.purity/*/alignments/dna/"${sample_id}".redux.bam; do
        if [[ -f "$cand" && -s "$cand" ]]; then
            echo "$cand"
            return 0
        fi
    done

    return 1
}

# Scan patients
patients=()
if [[ -n "${TARGET_PATIENT}" ]]; then
    patients=("${TARGET_PATIENT}")
else
    for pdir in "${SHEETS_DIR}"/BFIERCE_*; do
        [[ -d "$pdir" ]] || continue
        patients+=("$(basename "$pdir")")
    done
fi

submitted=0
for patient in "${patients[@]}"; do
    echo ""
    echo "--- Checking Patient: ${patient} ---"
    p_sheet_dir="${SHEETS_DIR}/${patient}"

    # Extract plasma samples from purity sample sheets
    plasma_samples=()
    for ps in "${p_sheet_dir}"/samplesheet.purity.*.csv "${p_sheet_dir}/samplesheet.purity.csv"; do
        [[ -f "$ps" ]] || continue
        while IFS=',' read -r gid sub sid stype seqtype ftype info fpath || [ -n "$sid" ]; do
            sid=$(echo "${sid:-}" | tr -d '\r' | xargs)
            stype=$(echo "${stype:-}" | tr -d '\r' | xargs)
            info=$(echo "${info:-}" | tr -d '\r' | xargs)
            [[ -z "$sid" || "$sid" == "sample_id" ]] && continue
            if [[ "$info" == *"longitudinal_sample"* || "$stype" =~ ^[pP]lasma || "$stype" == "tumor" ]]; then
                # Avoid picking up matched normal ID
                if [[ ! "$sid" =~ Normal && ! "$sid" =~ -N ]]; then
                    plasma_samples+=("$sid")
                fi
            fi
        done < "$ps"
    done

    # If no sheet found, fallback to directory scan
    if [[ ${#plasma_samples[@]} -eq 0 ]]; then
        for cand in "${V3_BASE}/${patient}/ichorCNA_plasma_only/"* \
                    "${V3_BASE}/${patient}/ichorCNA_plasma_vs_n/"* \
                    "${BAM_LINKS}/${patient}/plasma/"*.bam; do
            [[ -e "$cand" ]] || continue
            fname=$(basename "$cand")
            sid="${fname%%.*}"
            [[ -n "$sid" && "$sid" != "*" ]] && plasma_samples+=("$sid")
        done
    fi

    # Deduplicate plasma sample IDs
    IFS=$'\n' plasma_samples=($(sort -u <<<"${plasma_samples[*]}")); unset IFS

    if [[ ${#plasma_samples[@]} -eq 0 ]]; then
        echo "  [WARN] No plasma samples resolved for ${patient}; skipping."
        continue
    fi

    for sid in "${plasma_samples[@]}"; do
        target_bam=$(resolve_bfierce_plasma "${patient}" "${sid}" || true)
        if [[ -z "${target_bam}" || ! -f "${target_bam}" ]]; then
            echo "  [SKIP] Could not locate BAM for ${patient} sample ${sid}."
            continue
        fi

        sample_outdir="${V3_BASE}/${patient}/DELFI/${sid}"
        mkdir -p "${sample_outdir}/logs"

        # Check if already processed
        if [[ -f "${sample_outdir}/${sid}_delfi_tf_features.tsv" ]]; then
            echo "  [EXISTS] DELFI features already extracted for ${patient} (${sid})."
            continue
        fi

        echo "  Submitting DELFI-TF job for ${patient} (${sid})..."
        echo "    Input BAM:  ${target_bam}"
        echo "    Output Dir: ${sample_outdir}"

        model_arg=""
        if [[ -f "${MODEL_PATH}" ]]; then
            model_arg="${MODEL_PATH}"
        fi

        job_id=$(sbatch --parsable \
            -A "${ACCOUNT}" \
            -p "${PARTITION}" \
            --output="${sample_outdir}/logs/delfi_${sid}_%j.out" \
            --error="${sample_outdir}/logs/delfi_${sid}_%j.err" \
            "${DELFI_REPO}/slurm/submit_single.sh" \
            "${sid}" \
            "${target_bam}" \
            "${sample_outdir}" \
            "${model_arg}" \
            "${MAX_READS}")

        job_id=$(echo "${job_id}" | cut -d';' -f1)
        echo "    Job ID: ${job_id} | Logs -> ${sample_outdir}/logs/"
        submitted=$((submitted + 1))
    done
done

echo ""
echo "=================================================================="
echo " Total DELFI-TF Jobs Submitted: ${submitted}"
echo " When jobs complete, compile training set with:"
echo "   bash runs/run_BFIERCE_delfi_tf.sh --compile"
echo "=================================================================="
