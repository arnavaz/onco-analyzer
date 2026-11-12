#!/bin/bash
# ---------------------------------------------------------------------------
# submit_all_patients.sh — submit one oncoanalyser SLURM head job per patient
#
# Expects per-patient samplesheets created by make_oncoanalyser_samplesheet.py:
#   <samplesheet_dir>/<patient>/samplesheet.csv
#
# Usage:
#   bash Setup/submit_all_patients.sh <samplesheet_dir> <results_base> <mode> <genome> [account] [options]
#
# Example:
#   bash Setup/submit_all_patients.sh \
#       /cluster/projects/mohccn/STAGING/PANCAN_F/PC3-PMC-OSCC \
#       /cluster/projects/pughlab/results/oncoanalyser/PC3-PMC-OSCC \
#       wgts \
#       GRCh38_hmf \
#       pughlab
#
# Modes:   wgts | targeted | purity_estimate | panel_resource_creation | panel(alias->targeted)
# Genomes: GRCh37_hmf | GRCh38_hmf
#
# To resume all interrupted runs:
#   RESUME=true bash Setup/submit_all_patients.sh ...
#
# Options (forwarded to run_oncoanalyser.sh):
#   --panel-name <name>
#   --panel-config <path>          (repeatable)
#   --force-panel
#   --driver-gene-panel <path>
#   --target-regions-bed <path>
#   --isofox-gene-ids <path>
#   --ichor-exclusion-bed <path>
#   --ichor-pon <path>
#   --ichorWGS
#   --skip-nextflow
# ---------------------------------------------------------------------------

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  bash submit_all_patients.sh <samplesheet_dir> <results_base> <mode> <genome> [account] [options]

Modes:
  wgts | targeted | purity_estimate | panel_resource_creation | panel(alias->targeted)

Options:
  --panel-name <name>
  --panel-config <path>          (repeatable)
  --force-panel
  --processes-exclude <list>     Comma-separated list of processes to exclude (e.g. peach,lilac)
  --driver-gene-panel <path>
  --target-regions-bed <path>
  --isofox-gene-ids <path>
  --ichor-exclusion-bed <path>   Exclusion BED file for off-target ichorCNA run
  --ichor-pon <path>             Panel of Normals (PoN) RDS file for off-target ichorCNA run
  --ichorWGS                     Run ichorCNA for WGS samples (downsampled to 1X; auto-enabled for WGS modes)
  --no-ichorWGS                  Disable automatic ichorCNA for WGS samples
  --skip-nextflow                Bypass the Nextflow execution step and jump directly to post-processing
  --force-nextflow               Force Nextflow execution even if pipeline outputs already exist
  -h, --help
  -- [nextflow options...]       Any arguments after -- are forwarded directly to Nextflow
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

declare -a POSITIONAL=()
declare -a FORWARD_OPTS=()
SKIP_NEXTFLOW="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --panel-name|--panel-config|--driver-gene-panel|--target-regions-bed|--isofox-gene-ids|--processes-exclude|--processes_exclude|--ichor-exclusion-bed|--ichor_exclusion_bed|--ichor-pon|--ichor_pon)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            FORWARD_OPTS+=("$1" "$2")
            shift 2
            ;;
        --force-panel)
            FORWARD_OPTS+=("$1")
            shift
            ;;
        --ichorWGS|--ichorwgs)
            FORWARD_OPTS+=("--ichorWGS")
            shift
            ;;
        --no-ichorWGS|--no-ichorwgs|--no-ichor-wgs|--skip-ichorWGS|--skip-ichorwgs|--skip-ichor-wgs)
            FORWARD_OPTS+=("--no-ichorWGS")
            shift
            ;;
        --skip-nextflow|--skip_nextflow)
            FORWARD_OPTS+=("--skip-nextflow")
            SKIP_NEXTFLOW="true"
            shift
            ;;
        --force-nextflow|--force_nextflow)
            FORWARD_OPTS+=("--force-nextflow")
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                POSITIONAL+=("$1")
                shift
            done
            ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL[@]}"
[[ $# -ge 2 ]] || { echo "ERROR: Missing required arguments" >&2; usage; exit 1; }

SAMPLESHEET_DIR="$1"
RESULTS_BASE="$2"
MODE="${3:-wgts}"
GENOME="${4:-GRCh38_hmf}"
ACCOUNT="${5:-${SLURM_ACCOUNT:-}}"

if [[ "${MODE}" == "panel" ]]; then
    MODE="targeted"
fi

SOURCE_REPO="${SOURCE_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
WRAPPER="${SOURCE_REPO}/run_oncoanalyser.sh"

if [[ ! -f "${WRAPPER}" ]]; then
    echo "ERROR: run_oncoanalyser.sh not found at ${WRAPPER}" >&2
    exit 1
fi

WRAPPER_PARTITION="$(awk -F= '/^#SBATCH --partition=/{print $2; exit}' "${WRAPPER}" | xargs)"
if [[ -z "${ACCOUNT}" && "${WRAPPER_PARTITION}" == "pughlab" ]]; then
    ACCOUNT="pughlab"
fi

declare -a EXTRA_ARGS=()
if [[ $# -gt 5 ]]; then
    EXTRA_ARGS+=("--" "${@:6}")
fi

submitted=0
skipped=0

for patient_dir in "${SAMPLESHEET_DIR}"/*/; do
    [[ -d "${patient_dir}" ]] || continue
    patient=$(basename "${patient_dir}")
    samplesheet="${patient_dir}samplesheet.csv"

    if [[ ! -f "${samplesheet}" ]]; then
        echo "SKIP (no samplesheet): ${patient}"
        (( skipped++ )) || true
        continue
    fi

    outdir="${RESULTS_BASE}/${patient}"
    mkdir -p "${outdir}/logs"

    PARTITION_OPT=""
    if [[ "${SKIP_NEXTFLOW}" == "true" ]]; then
        PARTITION_OPT="--partition=all --time=04:00:00"
    fi

    echo "Submitting: ${patient}"
    if [[ -n "${ACCOUNT}" ]]; then
        echo "  account: ${ACCOUNT}"
        sbatch ${PARTITION_OPT} --export=ALL,SOURCE_REPO="${SOURCE_REPO}",NXF_OFFLINE=true -A "${ACCOUNT}" \
            --job-name="onco_${patient}" \
            --output="${outdir}/logs/onco_${patient}_%j.out" \
            --error="${outdir}/logs/onco_${patient}_%j.err" \
            "${WRAPPER}" \
            "${samplesheet}" \
            "${outdir}" \
            "${MODE}" \
            "${GENOME}" \
            "${ACCOUNT}" \
            "${FORWARD_OPTS[@]}" \
            "${EXTRA_ARGS[@]}"
    else
        sbatch ${PARTITION_OPT} --export=ALL,SOURCE_REPO="${SOURCE_REPO}",NXF_OFFLINE=true \
            --job-name="onco_${patient}" \
            --output="${outdir}/logs/onco_${patient}_%j.out" \
            --error="${outdir}/logs/onco_${patient}_%j.err" \
            "${WRAPPER}" \
            "${samplesheet}" \
            "${outdir}" \
            "${MODE}" \
            "${GENOME}" \
            "" \
            "${FORWARD_OPTS[@]}" \
            "${EXTRA_ARGS[@]}"
    fi
    (( submitted++ )) || true
done

echo ""
echo "Submitted: ${submitted} patient(s)"
if (( skipped > 0 )); then
    echo "Skipped:   ${skipped} (no samplesheet.csv found)"
fi
