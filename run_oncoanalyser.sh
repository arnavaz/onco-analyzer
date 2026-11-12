#!/bin/bash
#SBATCH --job-name=oncoanalyser
#SBATCH --partition=long
#SBATCH --account=pughlab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=10-00:00:00          # 10 days wall-time for the Nextflow head job
#SBATCH --output=logs/nf_%j.out
#SBATCH --error=logs/nf_%j.err

# ---------------------------------------------------------------------------
# nf-core/oncoanalyser — SLURM head-job submission wrapper
# Usage: sbatch run_oncoanalyser.sh <samplesheet.csv> <outdir> <mode> <genome> [account] [options]
#
# Modes:    wgts | targeted | purity_estimate | panel_resource_creation | panel(alias->targeted)
# Genomes:  GRCh37_hmf | GRCh38_hmf
# Workdir:  automatically set to <outdir>/work
#
# Options:
#   --panel-name <name>            Panel name passed to --panel (e.g. my_custom_panel)
#   --panel-config <path>          Extra Nextflow config for panel_data_paths (repeatable)
#   --force-panel                  Pass --force_panel for non-built-in panels
#   --driver-gene-panel <path>     Required for panel_resource_creation
#   --target-regions-bed <path>    Required for panel_resource_creation
#   --isofox-gene-ids <path>       Optional for panel_resource_creation (RNA-enabled panels)
#   --ichor-dir <path>             Pre-computed ichorCNA directory
#   --swgs-fastq-dir <path>        Optional directory containing sWGS FASTQs to run ichorCNA
#   --ichor-exclusion-bed <path>   Exclusion BED file for off-target ichorCNA run (triggers run)
#   --ichor-pon <path>             Panel of Normals (PoN) RDS file for off-target ichorCNA run
#   --ichorWGS                     Run ichorCNA for WGS samples (downsampled to 1X; auto-enabled for WGS modes)
#   --no-ichorWGS                  Disable automatic ichorCNA for WGS samples
#   --skip-nextflow                Bypass Nextflow run and run post-processing directly
#   --force-nextflow               Force Nextflow run even if pipeline outputs already exist
#   -h, --help                     Show help
# ---------------------------------------------------------------------------

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  sbatch run_oncoanalyser.sh <samplesheet.csv> <outdir> <mode> <genome> [account] [options]

Modes:
  wgts | targeted | purity_estimate | panel_resource_creation | panel(alias->targeted)

Options:
  --panel-name <name>            Panel name for targeted mode (maps to --panel)
  --panel-config <path>          Additional Nextflow config for panel data (repeatable)
  --force-panel                  Pass --force_panel (required for non-built-in panels)
  --profile <profile>            Nextflow profile to use (default: singularity, e.g. docker, conda)
  --processes-exclude <list>     Comma-separated list of processes to exclude (e.g. peach,lilac)
  --driver-gene-panel <path>     Required in panel_resource_creation mode
  --target-regions-bed <path>    Required in panel_resource_creation mode
  --isofox-gene-ids <path>       Optional in panel_resource_creation mode (RNA panels)
  --realign                      Automatically extract FASTQs from BAMs, regenerate samplesheet, and align using HMF genome build
  --methylation-data <path>      Optional path to external gene methylation data TSV
  --ichor-dir <path>             Pre-computed ichorCNA directory
  --swgs-fastq-dir <path>        Optional directory containing sWGS FASTQs to run ichorCNA
  --emseq-fastq-dir <path>       Optional directory containing EM-Seq FASTQs to run methylation pipeline
  --emseq-dir <path>             Pre-computed EM-Seq output directory
  --ichor-exclusion-bed <path>   Exclusion BED file for off-target ichorCNA run (triggers run)
  --ichor-pon <path>             Panel of Normals (PoN) RDS file for off-target ichorCNA run
  --ichor-sub-clones             Enable subclonal fraction estimation in off-target ichorCNA (default: false)
  --ichor-size-select [min,max]  Enrich plasma WGS for short fragments in single stream (default: 90,150)
  --no-ichor-size-select         Disable short-fragment size selection for plasma WGS
  --force-ichor                  Force ichorCNA execution even if outputs already exist
  --ichorWGS                     Run ichorCNA for WGS samples (tumour/normal/plasma downsampled to 1X; default for WGS modes)
  --no-ichorWGS                  Disable automatic ichorCNA for WGS samples
  --skip-nextflow                Bypass the Nextflow execution step and jump directly to post-processing
  --force-nextflow               Force Nextflow execution even if pipeline outputs already exist
  --purity-estimate-mode <mode>  Purity estimate run mode (default: wgts, choices: wgts|targeted)
  --purity-source <src>          Purity source for clinical reports: purple (default), auto, swgs, panel, vaf, pathology
  --sequencing-platform <plat>   Sequencing technology: illumina (default), ultima, sbx
  --metadata-tsv <path>          Optional clinical/cohort metadata TSV (with HLA alleles) for RNA post-processing
  --mixcr-path <path>            Optional MiXCR .clns file or exported clonotype TSV for cap-TCR analysis
  --delfi-tf                     Enable DELFI-TF tumor fraction estimation for plasma WGS (default: true)
  --no-delfi-tf                  Disable DELFI-TF estimation
  --delfi-repo <path>            Path to DELFI-TF repository (default: /cluster/projects/pughlab/src/DELFI-TF)
  --delfi-model <path>           Pretrained Random Forest model .rds for DELFI-TF inference
  --delfi-max-reads <N>          Maximum fragments for DELFI feature extraction (default: 50,000,000; 0 for all)
  --keep-downsampled-bams        Permanently preserve downsampled 1.0X and short-selected BAMs (default: true)
  --no-keep-downsampled-bams     Delete downsampled BAMs after readCounter finishes
  -h, --help                     Show this help message
  -- [nextflow options...]       Any arguments after -- are forwarded directly to Nextflow
EOF
}

# --- Required/optional arguments -------------------------------------------
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

PANEL_NAME=""
FORCE_PANEL="false"
DRIVER_GENE_PANEL=""
TARGET_REGIONS_BED=""
ISOFOX_GENE_IDS=""
PROCESSES_EXCLUDE=""
REALIGN="false"
METHYLATION_DATA=""
IHC_TSV=""
IHC_FILE=""
ORTHOGONAL_TSV=""
AMBER_EXPANDED_SNPS="false"
AGE=""
TISSUE=""
ICHOR_DIR=""
SWGS_FASTQ_DIR=""
EMSEQ_FASTQ_DIR=""
EMSEQ_DIR=""
ICHOR_EXCLUSION_BED=""
ICHOR_PON=""
ICHOR_SUB_CLONES="false"
ICHOR_WGS="false"
EXPLICIT_NO_ICHOR="false"
SKIP_NEXTFLOW="false"
FORCE_NEXTFLOW="${FORCE_NEXTFLOW:-false}"
PANEL_ICHOR_DIR=""
QUEUE_ARG=""
PROFILE="${PROFILE:-singularity}"
PURITY_ESTIMATE_MODE="${PURITY_ESTIMATE_MODE:-wgts}"
PURITY_SOURCE="${PURITY_SOURCE:-purple}"
SEQUENCING_PLATFORM="${SEQUENCING_PLATFORM:-illumina}"
METADATA_TSV=""
MIXCR_PATH=""
ICHOR_SIZE_SELECT="${ICHOR_SIZE_SELECT:-true}"
ICHOR_MIN_TLEN="${ICHOR_MIN_TLEN:-90}"
ICHOR_MAX_TLEN="${ICHOR_MAX_TLEN:-150}"
FORCE_ICHOR="${FORCE_ICHOR:-false}"
ENABLE_DELFI="${ENABLE_DELFI:-true}"
DELFI_REPO="${DELFI_REPO:-/cluster/projects/pughlab/src/DELFI-TF}"
DELFI_MODEL="${DELFI_MODEL:-${DELFI_REPO}/models/delfi_tf_default_model.rds}"
DELFI_MAX_READS="${DELFI_MAX_READS:-50000000}"
KEEP_DOWNSAMPLED_BAMS="${KEEP_DOWNSAMPLED_BAMS:-true}"
declare -a PANEL_CONFIGS=()
declare -a POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --purity-source|--purity_source)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            PURITY_SOURCE="$2"
            shift 2
            ;;
        --queue)
            [[ $# -ge 2 ]] || { echo "ERROR: --queue requires a value" >&2; usage; exit 1; }
            QUEUE_ARG="$2"
            shift 2
            ;;
        --profile)
            [[ $# -ge 2 ]] || { echo "ERROR: --profile requires a value" >&2; usage; exit 1; }
            PROFILE="$2"
            shift 2
            ;;
        --sequencing-platform|--sequencing_platform)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            SEQUENCING_PLATFORM="$2"
            shift 2
            ;;
        --panel-name)
            [[ $# -ge 2 ]] || { echo "ERROR: --panel-name requires a value" >&2; usage; exit 1; }
            PANEL_NAME="$2"
            shift 2
            ;;
        --panel-config)
            [[ $# -ge 2 ]] || { echo "ERROR: --panel-config requires a value" >&2; usage; exit 1; }
            PANEL_CONFIGS+=("$2")
            shift 2
            ;;
        --force-panel)
            FORCE_PANEL="true"
            shift
            ;;
        --processes-exclude|--processes_exclude)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            PROCESSES_EXCLUDE="$2"
            shift 2
            ;;
        --driver-gene-panel)
            [[ $# -ge 2 ]] || { echo "ERROR: --driver-gene-panel requires a value" >&2; usage; exit 1; }
            DRIVER_GENE_PANEL="$2"
            shift 2
            ;;
        --target-regions-bed)
            [[ $# -ge 2 ]] || { echo "ERROR: --target-regions-bed requires a value" >&2; usage; exit 1; }
            TARGET_REGIONS_BED="$2"
            shift 2
            ;;
        --isofox-gene-ids)
            [[ $# -ge 2 ]] || { echo "ERROR: --isofox-gene-ids requires a value" >&2; usage; exit 1; }
            ISOFOX_GENE_IDS="$2"
            shift 2
            ;;
        --realign)
            REALIGN="true"
            shift
            ;;
        --methylation-data|--methylation_data)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            METHYLATION_DATA="$2"
            shift 2
            ;;
        --ihc-tsv|--ihc_tsv)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            IHC_TSV="$2"
            shift 2
            ;;
        --ihc-file|--ihc_file)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            IHC_FILE="$2"
            shift 2
            ;;
        --orthogonal-tsv|--orthogonal_tsv|--known-mutations-tsv|--orthogonal-test-tsv)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            ORTHOGONAL_TSV="$2"
            shift 2
            ;;
        --amber-expanded-snps|--amber_expanded_snps)
            AMBER_EXPANDED_SNPS="true"
            shift
            ;;
        --age)
            [[ $# -ge 2 ]] || { echo "ERROR: --age requires a value" >&2; usage; exit 1; }
            AGE="$2"
            shift 2
            ;;
        --tissue)
            [[ $# -ge 2 ]] || { echo "ERROR: --tissue requires a value" >&2; usage; exit 1; }
            TISSUE="$2"
            shift 2
            ;;
        --ichor-dir|--ichor_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            ICHOR_DIR="$2"
            shift 2
            ;;
        --swgs-fastq-dir|--swgs_fastq_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            SWGS_FASTQ_DIR="$2"
            shift 2
            ;;
        --emseq-fastq-dir|--emseq_fastq_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            EMSEQ_FASTQ_DIR="$2"
            shift 2
            ;;
        --emseq-dir|--emseq_dir)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            EMSEQ_DIR="$2"
            shift 2
            ;;
        --ichor-exclusion-bed|--ichor_exclusion_bed)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            ICHOR_EXCLUSION_BED="$2"
            shift 2
            ;;
        --ichor-pon|--ichor_pon)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            ICHOR_PON="$2"
            shift 2
            ;;
        --ichor-sub-clones|--ichor_sub_clones)
            if [[ $# -ge 2 ]] && [[ "$2" =~ ^(true|false|TRUE|FALSE)$ ]]; then
                if [[ "$2" == "true" || "$2" == "TRUE" ]]; then
                    ICHOR_SUB_CLONES="true"
                else
                    ICHOR_SUB_CLONES="false"
                fi
                shift 2
            else
                ICHOR_SUB_CLONES="true"
                shift
            fi
            ;;
        --ichor-size-select|--ichor_size_select)
            if [[ $# -ge 2 ]] && [[ "$2" =~ ^[0-9]+,[0-9]+$ ]]; then
                ICHOR_SIZE_SELECT="true"
                ICHOR_MIN_TLEN=$(echo "$2" | cut -d',' -f1)
                ICHOR_MAX_TLEN=$(echo "$2" | cut -d',' -f2)
                shift 2
            else
                ICHOR_SIZE_SELECT="true"
                shift
            fi
            ;;
        --no-ichor-size-select|--no_ichor_size_select)
            ICHOR_SIZE_SELECT="false"
            shift
            ;;
        --force-ichor|--force_ichor)
            FORCE_ICHOR="true"
            shift
            ;;
        --ichorWGS|--ichorwgs)
            ICHOR_WGS="true"
            shift
            ;;
        --no-ichorWGS|--no-ichorwgs|--no-ichor-wgs|--skip-ichorWGS|--skip-ichorwgs|--skip-ichor-wgs)
            ICHOR_WGS="false"
            EXPLICIT_NO_ICHOR="true"
            shift
            ;;
        --skip-nextflow|--skip_nextflow)
            SKIP_NEXTFLOW="true"
            shift
            ;;
        --force-nextflow|--force_nextflow)
            FORCE_NEXTFLOW="true"
            shift
            ;;
        --purity-estimate-mode|--purity_estimate_mode)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            PURITY_ESTIMATE_MODE="$2"
            shift 2
            ;;
        --metadata-tsv|--metadata_tsv)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            METADATA_TSV="$2"
            shift 2
            ;;
        --mixcr-path|--mixcr_path)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            MIXCR_PATH="$2"
            shift 2
            ;;
        --delfi-tf|--enable-delfi|--enable_delfi)
            ENABLE_DELFI="true"
            shift
            ;;
        --no-delfi-tf|--disable-delfi|--disable_delfi)
            ENABLE_DELFI="false"
            shift
            ;;
        --delfi-repo|--delfi_repo)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            DELFI_REPO="$2"
            shift 2
            ;;
        --delfi-model|--delfi_model)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            DELFI_MODEL="$2"
            shift 2
            ;;
        --delfi-max-reads|--delfi_max_reads)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; usage; exit 1; }
            DELFI_MAX_READS="$2"
            shift 2
            ;;
        --keep-downsampled-bams|--keep_downsampled_bams)
            KEEP_DOWNSAMPLED_BAMS="true"
            shift
            ;;
        --no-keep-downsampled-bams)
            KEEP_DOWNSAMPLED_BAMS="false"
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
SAMPLESHEET="$1"
OUTDIR="$2"
MODE="${3:-wgts}"
GENOME="${4:-GRCh38_hmf}"
ACCOUNT_ARG="${5:-}"
WORK_DIR="${OUTDIR}/work"
TMP_DIR="${TMPDIR:-${OUTDIR}/tmp}"

# Setup user-controlled temporary directories to prevent exhaustion of compute node /tmp
mkdir -p "${OUTDIR}" "${WORK_DIR}" "${TMP_DIR}"
export TMPDIR="${TMP_DIR}"
export TEMP="${TMP_DIR}"
export TMP="${TMP_DIR}"
export SINGULARITY_TMPDIR="${TMP_DIR}"
export APPTAINER_TMPDIR="${TMP_DIR}"
export SINGULARITYENV_TMPDIR="${TMP_DIR}"
export SINGULARITYENV_TEMP="${TMP_DIR}"
export SINGULARITYENV_TMP="${TMP_DIR}"

case "${MODE}" in
    panel)
        echo "WARN: mode 'panel' is deprecated in wrapper usage; mapping to 'targeted'" >&2
        MODE="targeted"
        ;;
    wgts|targeted|purity_estimate|panel_resource_creation)
        ;;
    *)
        echo "ERROR: Unsupported mode: ${MODE}" >&2
        usage
        exit 1
        ;;
esac

# Auto-enable ichorCNA for WGS modes unless explicitly disabled
if [[ "${EXPLICIT_NO_ICHOR}" != "true" ]]; then
    if [[ "${MODE}" == "wgts" || ( "${MODE}" == "purity_estimate" && "${PURITY_ESTIMATE_MODE}" == "wgts" ) ]]; then
        ICHOR_WGS="true"
    fi
fi

# --- Paths & environment ----------------------------------------------------
SOURCE_REPO="${SOURCE_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SING_CACHE="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache"
DEFAULT_SOURCE_REPO="/cluster/projects/pughlab/src/oncoanalyser_multimmr"

# If script is executed from a SLURM spool copy, BASH_SOURCE may resolve to
# /var/spool/slurmd/... and not the actual repository path.
if [[ ! -f "${SOURCE_REPO}/run_oncoanalyser.sh" && -f "${DEFAULT_SOURCE_REPO}/run_oncoanalyser.sh" ]]; then
    SOURCE_REPO="${DEFAULT_SOURCE_REPO}"
fi

PIPELINE_DIR="${SOURCE_REPO}/pipeline"
if [[ ! -d "${PIPELINE_DIR}" ]]; then
    PIPELINE_DIR="/cluster/projects/pughlab/bin/oncoanalyser"
fi

CUSTOM_CONFIG="${SOURCE_REPO}/conf/h4h.config"

# Validate mode-specific options and files up front for fast failure.
if [[ "${MODE}" == "targeted" ]]; then
    if [[ -n "${PANEL_NAME}" && "${PANEL_NAME}" != "tso500" ]]; then
        [[ "${FORCE_PANEL}" == "true" ]] || {
            echo "ERROR: Non-built-in panel '${PANEL_NAME}' requires --force-panel" >&2
            exit 1
        }
        [[ ${#PANEL_CONFIGS[@]} -gt 0 ]] || {
            echo "ERROR: Non-built-in panel '${PANEL_NAME}' requires at least one --panel-config file" >&2
            exit 1
        }
    fi
    if [[ "${FORCE_PANEL}" == "true" && -z "${PANEL_NAME}" ]]; then
        echo "ERROR: --force-panel requires --panel-name in targeted mode" >&2
        exit 1
    fi
fi

if [[ "${MODE}" == "panel_resource_creation" ]]; then
    [[ -n "${DRIVER_GENE_PANEL}" ]] || {
        echo "ERROR: panel_resource_creation requires --driver-gene-panel" >&2
        exit 1
    }
    [[ -n "${TARGET_REGIONS_BED}" ]] || {
        echo "ERROR: panel_resource_creation requires --target-regions-bed" >&2
        exit 1
    }
fi

# Load required modules
module load java/21
module load singularity/4.3.2

# Ensure Nextflow is on PATH (installed to ~/.local/bin)
export PATH="${HOME}/.local/bin:${PATH}"

# Preflight checks to fail early with clear messages.
[[ -f "${PIPELINE_DIR}/main.nf" ]] || { echo "ERROR: Pipeline not found: ${PIPELINE_DIR}/main.nf" >&2; exit 1; }
[[ -f "${SOURCE_REPO}/run_oncoanalyser.sh" ]] || { echo "ERROR: Source repo not found: ${SOURCE_REPO}" >&2; exit 1; }
[[ -f "${CUSTOM_CONFIG}" ]] || { echo "ERROR: Config not found: ${CUSTOM_CONFIG}" >&2; exit 1; }
[[ -f "${SAMPLESHEET}" ]] || { echo "ERROR: Samplesheet not found: ${SAMPLESHEET}" >&2; exit 1; }
for panel_cfg in "${PANEL_CONFIGS[@]}"; do
    [[ -f "${panel_cfg}" ]] || { echo "ERROR: Panel config not found: ${panel_cfg}" >&2; exit 1; }
done
[[ -z "${DRIVER_GENE_PANEL}" || -f "${DRIVER_GENE_PANEL}" ]] || { echo "ERROR: Driver gene panel not found: ${DRIVER_GENE_PANEL}" >&2; exit 1; }
[[ -z "${TARGET_REGIONS_BED}" || -f "${TARGET_REGIONS_BED}" ]] || { echo "ERROR: Target regions BED not found: ${TARGET_REGIONS_BED}" >&2; exit 1; }
[[ -z "${ISOFOX_GENE_IDS}" || -f "${ISOFOX_GENE_IDS}" ]] || { echo "ERROR: Isofox gene IDs not found: ${ISOFOX_GENE_IDS}" >&2; exit 1; }

# ---------------------------------------------------------------------------
# RNA-only samplesheet auto-detection
# ---------------------------------------------------------------------------
HAS_DNA_SAMPLES="false"
HAS_RNA_SAMPLES="false"
if awk -F',' 'NR>1 && tolower($5)=="dna" {exit 0} END {exit 1}' "${SAMPLESHEET}" 2>/dev/null; then
    HAS_DNA_SAMPLES="true"
fi
if awk -F',' 'NR>1 && tolower($5)=="rna" {exit 0} END {exit 1}' "${SAMPLESHEET}" 2>/dev/null; then
    HAS_RNA_SAMPLES="true"
fi

IS_RNA_ONLY="false"
if [[ "${HAS_RNA_SAMPLES}" == "true" && "${HAS_DNA_SAMPLES}" == "false" ]]; then
    IS_RNA_ONLY="true"
    echo "==================================================================" >&2
    echo "Detected RNA-only samplesheet: auto-configuring RNA execution mode." >&2
    echo "==================================================================" >&2

    # Auto-exclude DNA-dependent processes if not already excluded
    for p in neo virusinterpreter; do
        if ! [[ ",${PROCESSES_EXCLUDE}," =~ ,${p}, ]]; then
            if [[ -z "${PROCESSES_EXCLUDE}" ]]; then
                PROCESSES_EXCLUDE="${p}"
            else
                PROCESSES_EXCLUDE="${PROCESSES_EXCLUDE},${p}"
            fi
        fi
    done

    # Disable WGS ichorCNA (DNA-only)
    ICHOR_WGS="false"
fi

# Singularity binary path (loaded via module; override if needed)
export SINGULARITY_CACHEDIR="${SING_CACHE}"
export NXF_SINGULARITY_CACHEDIR="${SING_CACHE}"

# Prevent Nextflow from pulling containers or configs from the internet
# Comment this out only if you need to update containers or the pipeline
export NXF_OFFLINE=true

# Propagate a selected SLURM account to child task submissions.
# Precedence: wrapper arg > exported ONCO_SLURM_ACCOUNT > SLURM job account/env.
ONCO_SLURM_ACCOUNT_EFFECTIVE="${ACCOUNT_ARG:-${ONCO_SLURM_ACCOUNT:-${SLURM_JOB_ACCOUNT:-${SLURM_ACCOUNT:-}}}}"
export ONCO_SLURM_ACCOUNT="${ONCO_SLURM_ACCOUNT_EFFECTIVE}"

# Propagate a selected SLURM queue/partition to child task submissions only if explicitly specified.
# Precedence: wrapper arg > exported ONCO_SLURM_QUEUE. Do NOT inherit the head job's partition.
ONCO_SLURM_QUEUE_EFFECTIVE="${QUEUE_ARG:-${ONCO_SLURM_QUEUE:-}}"
export ONCO_SLURM_QUEUE="${ONCO_SLURM_QUEUE_EFFECTIVE}"

# ---------------------------------------------------------------------------
# Create required directories
# ---------------------------------------------------------------------------
mkdir -p "${OUTDIR}" "${WORK_DIR}" "${SING_CACHE}"
mkdir -p logs 2>/dev/null || true

# ---------------------------------------------------------------------------
# Pipeline completion check (auto-skip Nextflow if outputs already exist)
# ---------------------------------------------------------------------------
PIPELINE_ALREADY_COMPLETED="false"
if [[ "${FORCE_NEXTFLOW}" != "true" && "${SKIP_NEXTFLOW}" != "true" ]]; then
    if [[ "${IS_RNA_ONLY}" == "true" ]]; then
        isofox_done=$(find "${OUTDIR}" -maxdepth 4 -type f \( -name "*.isf.gene_data.*" -o -name "*.isf.summary.*" \) \
            ! -path "${OUTDIR}/work/*" \
            ! -path "${OUTDIR}/.nextflow/*" \
            ! -path "${OUTDIR}/logs/*" 2>/dev/null | head -n 1)
        if [[ -n "${isofox_done}" && -s "${isofox_done}" ]]; then
            PIPELINE_ALREADY_COMPLETED="true"
        fi
    elif [[ "${MODE}" == "wgts" || "${MODE}" == "targeted" ]]; then
        purple_done=$(find "${OUTDIR}" -maxdepth 4 -type f -name "*.purple.purity.tsv" \
            ! -path "${OUTDIR}/work/*" \
            ! -path "${OUTDIR}/.nextflow/*" \
            ! -path "${OUTDIR}/logs/*" 2>/dev/null | head -n 1)

        # In wgts and targeted modes, ORANGE report is the terminal output unless excluded
        orange_done=""
        if [[ "${PROCESSES_EXCLUDE:-}" =~ (^|,)[[:space:]]*(orange|ORANGE)[[:space:]]*(,|$) ]]; then
            orange_done="excluded"
        else
            orange_done=$(find "${OUTDIR}" -maxdepth 4 -type f \( -name "*.orange.pdf" -o -name "*.orange.json.gz" \) \
                ! -path "${OUTDIR}/work/*" \
                ! -path "${OUTDIR}/.nextflow/*" \
                ! -path "${OUTDIR}/logs/*" 2>/dev/null | head -n 1)
        fi

        if [[ -n "${purple_done}" && -s "${purple_done}" && -n "${orange_done}" ]]; then
            PIPELINE_ALREADY_COMPLETED="true"
        fi
    elif [[ "${MODE}" == "purity_estimate" ]]; then
        wisp_done=$(find "${OUTDIR}" -maxdepth 4 -type f \( -name "*.wisp.summary.tsv" -o -name "*.wisp.tsv" \) \
            ! -path "${OUTDIR}/work/*" \
            ! -path "${OUTDIR}/.nextflow/*" \
            ! -path "${OUTDIR}/logs/*" 2>/dev/null | head -n 1)
        if [[ -n "${wisp_done}" && -s "${wisp_done}" ]]; then
            PIPELINE_ALREADY_COMPLETED="true"
        fi
    fi

    if [[ "${PIPELINE_ALREADY_COMPLETED}" == "true" ]]; then
        echo "==================================================================" >&2
        echo "Pipeline outputs already completed successfully in ${OUTDIR}." >&2
        echo "Skipping Nextflow execution and proceeding to post-processing / ichorCNA." >&2
        echo "==================================================================" >&2
        SKIP_NEXTFLOW="true"
    fi
fi

# ---------------------------------------------------------------------------
# Realignment extraction preprocessing
# ---------------------------------------------------------------------------
if [[ "${REALIGN}" == "true" && "${SKIP_NEXTFLOW}" != "true" ]]; then
    echo "Initializing BAM/CRAM back-extraction and realignment..."
    FASTQ_OUTDIR="${OUTDIR}/realigned_fastqs"
    mkdir -p "${FASTQ_OUTDIR}"

    # Parse BAM files from samplesheet using inline python to avoid bash CSV parsing issues
    mapfile -t BAM_ROWS < <(python3 -c "
import csv, sys
with open(sys.argv[1], 'r') as f:
    reader = csv.DictReader(row for row in f if row.strip() and not row.startswith('#'))
    for row in reader:
        if row.get('filetype') in ['bam', 'cram']:
            print(f\"{row['sample_id']}\t{row['filepath']}\")
" "${SAMPLESHEET}")

    if [[ ${#BAM_ROWS[@]} -eq 0 ]]; then
        echo "No BAM/CRAM files found in samplesheet; skipping extraction."
    else
        echo "Found ${#BAM_ROWS[@]} BAM/CRAM file(s) for extraction."
        declare -a JOB_IDS=()
        for row in "${BAM_ROWS[@]}"; do
            sample_id=$(echo "$row" | cut -f1)
            filepath=$(echo "$row" | cut -f2)

            r1="${FASTQ_OUTDIR}/${sample_id}_L001_R1.fastq.gz"
            r2="${FASTQ_OUTDIR}/${sample_id}_L001_R2.fastq.gz"

            skip_extraction="false"
            if [[ -f "$r1" && -f "$r2" && -s "$r1" && -s "$r2" ]]; then
                echo "Found existing FASTQ files for sample ${sample_id}. Verifying read counts..."
                if r1_lines=$(gzip -dc "$r1" 2>/dev/null | wc -l) && \
                   r2_lines=$(gzip -dc "$r2" 2>/dev/null | wc -l) && \
                   [[ "$r1_lines" -eq "$r2_lines" && "$r1_lines" -gt 0 ]]; then
                    reads=$((r1_lines / 4))
                    echo "  Verification successful: both FASTQ files exist, are valid, and contain ${reads} reads. Skipping extraction."
                    skip_extraction="true"
                else
                    echo "  Verification failed: FASTQ files are mismatched, empty, or corrupted. Re-extracting."
                fi
            fi

            if [[ "${skip_extraction}" == "true" ]]; then
                continue
            fi

            echo "Submitting extraction job for sample: ${sample_id} (${filepath})"
            
            # Formulate the wrap script
            # We want to make sure it loads samtools/1.20 and runs sort and fastq
            WRAP_COMMAND="module load samtools/1.20 && samtools sort -n -@ 8 \"${filepath}\" -o \"${FASTQ_OUTDIR}/${sample_id}.name_sorted.bam\" && samtools fastq -@ 8 -1 \"${FASTQ_OUTDIR}/${sample_id}_L001_R1.fastq.gz\" -2 \"${FASTQ_OUTDIR}/${sample_id}_L001_R2.fastq.gz\" -0 /dev/null -s /dev/null -n \"${FASTQ_OUTDIR}/${sample_id}.name_sorted.bam\" && rm \"${FASTQ_OUTDIR}/${sample_id}.name_sorted.bam\""

            # Submit via sbatch
            # Extract account options
            ACCOUNT_OPT=""
            if [[ -n "${ONCO_SLURM_ACCOUNT}" ]]; then
                ACCOUNT_OPT="-A ${ONCO_SLURM_ACCOUNT}"
            fi

            # Submit and capture job ID
            job_id=$(sbatch --parsable --cpus-per-task=8 --mem=16G --time=12:00:00 --partition=all ${ACCOUNT_OPT} --job-name="extract_${sample_id}" --output="logs/extract_${sample_id}_%j.out" --error="logs/extract_${sample_id}_%j.err" --wrap="${WRAP_COMMAND}")
            echo "  Submitted job ID: ${job_id}"
            JOB_IDS+=("${job_id}")
        done

        # Wait for all extraction jobs to complete
        if [[ ${#JOB_IDS[@]} -gt 0 ]]; then
            echo "Waiting for extraction jobs: ${JOB_IDS[*]}"
            JOB_IDS_COMMA=$(IFS=,; echo "${JOB_IDS[*]}")
            while true; do
                active_jobs=$(squeue -j "${JOB_IDS_COMMA}" -h -o "%A" 2>/dev/null || true)
                if [[ -z "${active_jobs}" ]]; then
                    break
                fi
                sleep 30
            done

            # Verify exit status of all extraction jobs
            failed_jobs=()
            for jid in "${JOB_IDS[@]}"; do
                state=$(sacct -j "$jid" -n -o State | head -n 1 | xargs || true)
                if [[ -n "$state" && "$state" != "COMPLETED" ]]; then
                    failed_jobs+=("${jid} (${state})")
                fi
            done

            if [[ ${#failed_jobs[@]} -gt 0 ]]; then
                echo "ERROR: The following extraction jobs did not complete successfully: ${failed_jobs[*]}" >&2
                exit 1
            fi
            echo "All extraction jobs completed successfully."
        else
            echo "No extraction jobs were submitted (all samples skipped or no BAMs)."
        fi
        
        # Verify output files exist and are non-empty
        for row in "${BAM_ROWS[@]}"; do
            sample_id=$(echo "$row" | cut -f1)
            r1="${FASTQ_OUTDIR}/${sample_id}_L001_R1.fastq.gz"
            r2="${FASTQ_OUTDIR}/${sample_id}_L001_R2.fastq.gz"
            if [[ ! -s "$r1" || ! -s "$r2" ]]; then
                echo "ERROR: Expected output FASTQ files not found or empty for sample: ${sample_id}" >&2
                echo "  Expected: $r1 and $r2" >&2
                exit 1
            fi
        done
        echo "All extraction jobs completed successfully."
    fi

    # Generate the new FASTQ-based samplesheet
    REIMAGE_SAMPLESHEET="${OUTDIR}/samplesheet.realigned.csv"
    echo "Generating realigned samplesheet: ${REIMAGE_SAMPLESHEET}"
    python3 -c "
import csv, sys
fastq_dir = sys.argv[2]
with open(sys.argv[1], 'r') as f_in, open(sys.argv[3], 'w', newline='') as f_out:
    reader = csv.DictReader(row for row in f_in if row.strip() and not row.startswith('#'))
    writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames)
    writer.writeheader()
    for row in reader:
        if row.get('filetype') in ['bam', 'cram']:
            sample_id = row['sample_id']
            row['filetype'] = 'fastq'
            row['info'] = f'library_id:{sample_id};lane:001'
            row['filepath'] = f'{fastq_dir}/{sample_id}_L001_R1.fastq.gz;{fastq_dir}/{sample_id}_L001_R2.fastq.gz'
        writer.writerow(row)
" "${SAMPLESHEET}" "${FASTQ_OUTDIR}" "${REIMAGE_SAMPLESHEET}"

    # Use the new samplesheet for nextflow
    SAMPLESHEET="${REIMAGE_SAMPLESHEET}"
fi

# ---------------------------------------------------------------------------
# Optional: resume a previous run by passing -resume to nextflow
# Set RESUME=true before calling sbatch to resume:
#   RESUME=true sbatch run_oncoanalyser.sh ...
# ---------------------------------------------------------------------------
INITIAL_RESUME="${RESUME:-false}"

# Optional: automatic retry with -resume if Nextflow exits non-zero.
# Useful for transient scheduler failures (e.g. submit rate/QOS spikes).
#
# Controls:
#   AUTO_RESUME_ON_FAILURE=true|false   (default: true)
#   MAX_RESUME_ATTEMPTS=<int>           (default: 3)
#   RESUME_BACKOFF_SECONDS=<int>        (default: 300)
AUTO_RESUME_ON_FAILURE="${AUTO_RESUME_ON_FAILURE:-true}"
MAX_RESUME_ATTEMPTS="${MAX_RESUME_ATTEMPTS:-3}"
RESUME_BACKOFF_SECONDS="${RESUME_BACKOFF_SECONDS:-300}"

# Optional: convert final BAM outputs to CRAM after a successful pipeline run.
#
# Controls:
#   CONVERT_FINAL_BAM_TO_CRAM=true|false  (default: true)
CONVERT_FINAL_BAM_TO_CRAM="${CONVERT_FINAL_BAM_TO_CRAM:-true}"

if ! [[ "${MAX_RESUME_ATTEMPTS}" =~ ^[0-9]+$ ]] || [[ "${MAX_RESUME_ATTEMPTS}" -lt 1 ]]; then
    echo "ERROR: MAX_RESUME_ATTEMPTS must be an integer >= 1 (got: ${MAX_RESUME_ATTEMPTS})" >&2
    exit 1
fi
if ! [[ "${RESUME_BACKOFF_SECONDS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: RESUME_BACKOFF_SECONDS must be an integer >= 0 (got: ${RESUME_BACKOFF_SECONDS})" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------------
echo "========================================"
echo "  nf-core/oncoanalyser"
echo "  Mode:        ${MODE}"
echo "  Genome:      ${GENOME}"
echo "  Samplesheet: ${SAMPLESHEET}"
echo "  Output:      ${OUTDIR}"
echo "  Work dir:    ${WORK_DIR}"
echo "  NXF_HOME:    ${NXF_HOME:-$HOME/.nextflow}"
echo "  NXF_OFFLINE: ${NXF_OFFLINE}"
echo "  Account:     ${ONCO_SLURM_ACCOUNT:-auto}"
echo "  Date:        $(date)"
echo "  Auto-resume: ${AUTO_RESUME_ON_FAILURE} (attempts=${MAX_RESUME_ATTEMPTS}, backoff=${RESUME_BACKOFF_SECONDS}s)"
echo "  BAM->CRAM:   ${CONVERT_FINAL_BAM_TO_CRAM}"
echo "  Profile:     ${PROFILE}"
echo "  Panel name:  ${PANEL_NAME:-none}"
echo "  Realign:     ${REALIGN}"
echo "  Force panel: ${FORCE_PANEL}"
echo "  Ichor WGS:   ${ICHOR_WGS}"
echo "  Skip NF:     ${SKIP_NEXTFLOW}"
echo "  Panel cfgs:  ${#PANEL_CONFIGS[@]}"
if [[ "${MODE}" == "panel_resource_creation" ]]; then
    echo "  Driver gene: ${DRIVER_GENE_PANEL}"
    echo "  Regions BED: ${TARGET_REGIONS_BED}"
    echo "  Isofox IDs:  ${ISOFOX_GENE_IDS:-none}"
fi
echo "========================================"

# Run from OUTDIR so session state and transient files stay local to patient output.
cd "${OUTDIR}"

# Compute base patient directory (e.g. if OUTDIR is <patient>/plasma/<plasma_id>, PATIENT_BASEDIR is <patient>)
PATIENT_BASEDIR="${OUTDIR}"
if [[ "$(basename "$(dirname "${OUTDIR}")")" == "plasma" ]]; then
    PATIENT_BASEDIR="$(dirname "$(dirname "${OUTDIR}")")"
elif [[ "$(basename "${OUTDIR}")" == "plasma" ]]; then
    PATIENT_BASEDIR="$(dirname "${OUTDIR}")"
fi

# ---------------------------------------------------------------------------
# Launch parallel support workflows (sWGS ichorCNA & EM-Seq methylation)
# These run concurrently with Nextflow; the downstream summary step will
# wait for them to complete before generating the final report.
# ---------------------------------------------------------------------------
declare -a PARALLEL_SUPPORT_JOB_IDS=()
declare -a PARALLEL_SUPPORT_PIDS=()

ACCOUNT_OPT=""
if [[ -n "${ONCO_SLURM_ACCOUNT}" ]]; then
    ACCOUNT_OPT="-A ${ONCO_SLURM_ACCOUNT}"
fi

# Determine partition for high-memory support jobs (>30 GB RAM)
# On H4H: partition 'all' max is 30G. High-memory jobs MUST use 'himem' (or 'pughlab' if on pughlab account).
PARTITION_SUPPORT="himem"
if [[ "${ONCO_SLURM_ACCOUNT}" == "pughlab" ]]; then
    PARTITION_SUPPORT="pughlab,himem"
fi
if [[ -n "${ONCO_SLURM_QUEUE}" ]]; then
    PARTITION_SUPPORT="${ONCO_SLURM_QUEUE}"
fi

# 1. Parallel sWGS ichorCNA
if [[ -n "${SWGS_FASTQ_DIR}" && -d "${SWGS_FASTQ_DIR}" ]]; then
    if [[ -z "${ICHOR_DIR}" || ! -d "${ICHOR_DIR}" ]] || ! find "${ICHOR_DIR}" -type f -name "*params.txt" 2>/dev/null | grep -q .; then
        mapfile -t SWGS_TUMOR_SAMPLES < <(awk -F',' '$4=="tumor" && $5=="dna" {print $3}' "${SAMPLESHEET}" | sort -u)
        SWGS_SUBJECT_ID=$(awk -F',' 'NR>1 && $2!="" {print $2; exit}' "${SAMPLESHEET}")

        for t_id in "${SWGS_TUMOR_SAMPLES[@]}"; do
            SAMPLE_ICHOR_DIR="${OUTDIR}/ichorCNA_t_vs_n/${t_id}"
            if [[ -f "${SAMPLE_ICHOR_DIR}/${t_id}.params.txt" ]] || [[ -f "${OUTDIR}/ichorCNA/${t_id}/${t_id}.params.txt" ]] || [[ -f "${OUTDIR}/ichorCNA/${t_id}.params.txt" ]] || [[ -f "${OUTDIR}/ichorCNA/params.txt" ]]; then
                echo "sWGS ichorCNA results already present for ${t_id} (skipping execution)" >&2
                ICHOR_DIR="${OUTDIR}/ichorCNA_t_vs_n"
                continue
            fi

            mkdir -p "${OUTDIR}/logs" "${SAMPLE_ICHOR_DIR}"
            echo "Launching sWGS ichorCNA for ${t_id} in parallel with Nextflow..." >&2
            SWGS_CMD="bash '${SOURCE_REPO}/scripts/run_swgs_ichorcna.sh' --swgs-fastq-dir '${SWGS_FASTQ_DIR}' --sample-id '${t_id}' --patient-id '${SWGS_SUBJECT_ID}' --outdir '${SAMPLE_ICHOR_DIR}' --genome '${GENOME}' --threads 8 --source-repo '${SOURCE_REPO}'"

            if command -v sbatch >/dev/null 2>&1; then
                job_id=$(sbatch --parsable --cpus-per-task=8 --mem=58G --time=12:00:00 --partition=${PARTITION_SUPPORT} ${ACCOUNT_OPT} \
                    --job-name="swgs_${t_id}" \
                    --output="${OUTDIR}/logs/swgs_${t_id}_%j.out" \
                    --error="${OUTDIR}/logs/swgs_${t_id}_%j.err" \
                    --wrap="${SWGS_CMD}")
                echo "  Submitted parallel sWGS ichorCNA SLURM job ID: ${job_id}" >&2
                PARALLEL_SUPPORT_JOB_IDS+=("${job_id}")
            else
                bash -c "${SWGS_CMD}" > "${OUTDIR}/logs/swgs_${t_id}.log" 2>&1 &
                pid=$!
                echo "  Started parallel sWGS ichorCNA background PID: ${pid}" >&2
                PARALLEL_SUPPORT_PIDS+=("${pid}")
            fi
        done
    fi
fi

# 2. Parallel EM-Seq Methylation
if [[ -n "${EMSEQ_FASTQ_DIR}" && -d "${EMSEQ_FASTQ_DIR}" ]] || [[ -n "${EMSEQ_DIR}" && -d "${EMSEQ_DIR}" ]]; then
    CAND_METH=""
    if [[ -n "${METHYLATION_DATA}" && -f "${METHYLATION_DATA}" ]]; then
        CAND_METH="${METHYLATION_DATA}"
    elif [[ -n "${EMSEQ_DIR}" && -d "${EMSEQ_DIR}" ]]; then
        CAND_METH=$(find "${EMSEQ_DIR}" -type f -name "*CpG_methylation_matrix.RData" 2>/dev/null | head -n 1 || true)
    elif [[ -d "${OUTDIR}/MethylDackel" ]]; then
        CAND_METH=$(find "${OUTDIR}/MethylDackel" -type f -name "*CpG_methylation_matrix.RData" 2>/dev/null | head -n 1 || true)
    elif [[ -d "${OUTDIR}/EMSeq/MethylDackel" ]]; then
        CAND_METH=$(find "${OUTDIR}/EMSeq/MethylDackel" -type f -name "*CpG_methylation_matrix.RData" 2>/dev/null | head -n 1 || true)
    fi

    if [[ -n "${CAND_METH}" && -f "${CAND_METH}" ]]; then
        echo "EM-Seq methylation matrix already present: ${CAND_METH} (skipping execution)" >&2
        METHYLATION_DATA="${CAND_METH}"
    elif [[ -n "${EMSEQ_FASTQ_DIR}" && -d "${EMSEQ_FASTQ_DIR}" ]]; then
        EM_SUBJECT_ID=$(awk -F',' 'NR>1 && $2!="" {print $2; exit}' "${SAMPLESHEET}")

        # Check if matching EM-Seq FASTQs exist for this patient before submitting
        P_NUM=$(echo "${EM_SUBJECT_ID}" | grep -oE "P[0-9]+" | sed 's/P//' || true)
        HAS_EM_FASTQS="false"
        if [[ -n "${P_NUM}" ]] && find "${EMSEQ_FASTQ_DIR}" -maxdepth 2 -type f \( -name "*R1*.fastq.gz" -o -name "*R1*.fq.gz" \) 2>/dev/null | grep -qE "([TN]${P_NUM}-|-P${P_NUM}-|DYNA-P${P_NUM})"; then
            HAS_EM_FASTQS="true"
        elif [[ -f "${SAMPLESHEET}" ]]; then
            while IFS=, read -r _grp _sub smp_id _rest; do
                token=$(echo "${smp_id}" | grep -oE "[0-9]{2}-[0-9]{4,5}" || true)
                if [[ -n "${token}" ]] && find "${EMSEQ_FASTQ_DIR}" -maxdepth 2 -type f \( -name "*R1*.fastq.gz" -o -name "*R1*.fq.gz" \) 2>/dev/null | grep -q "${token}"; then
                    HAS_EM_FASTQS="true"
                    break
                fi
            done < <(grep -v "^#" "${SAMPLESHEET}" 2>/dev/null | grep -v "^group_id" || true)
        fi

        if [[ "${HAS_EM_FASTQS}" != "true" ]]; then
            echo "No matching EM-Seq FASTQs found for ${EM_SUBJECT_ID} in ${EMSEQ_FASTQ_DIR} (skipping execution)" >&2
        else
            EM_TARGET_OUT="${EMSEQ_DIR:-${OUTDIR}/EMSeq}"
            mkdir -p "${OUTDIR}/logs" "${EM_TARGET_OUT}"

            echo "Launching EM-Seq methylation pipeline for ${EM_SUBJECT_ID} in parallel with Nextflow..." >&2
            EMSEQ_CMD="bash '${SOURCE_REPO}/scripts/run_emseq_pipeline.sh' --emseq-fastq-dir '${EMSEQ_FASTQ_DIR}' --outdir '${EM_TARGET_OUT}' --patient-id '${EM_SUBJECT_ID}' --samplesheet '${SAMPLESHEET}' --project 'MultiMMR' --threads 8 --source-repo '${SOURCE_REPO}'"

            if command -v sbatch >/dev/null 2>&1; then
                job_id=$(sbatch --parsable --cpus-per-task=8 --mem=58G --time=24:00:00 --partition=${PARTITION_SUPPORT} ${ACCOUNT_OPT} \
                    --job-name="emseq_${EM_SUBJECT_ID}" \
                    --output="${OUTDIR}/logs/emseq_${EM_SUBJECT_ID}_%j.out" \
                    --error="${OUTDIR}/logs/emseq_${EM_SUBJECT_ID}_%j.err" \
                    --wrap="${EMSEQ_CMD}")
                echo "  Submitted parallel EM-Seq pipeline SLURM job ID: ${job_id}" >&2
                PARALLEL_SUPPORT_JOB_IDS+=("${job_id}")
            else
                bash -c "${EMSEQ_CMD}" > "${OUTDIR}/logs/emseq_${EM_SUBJECT_ID}.log" 2>&1 &
                pid=$!
                echo "  Started parallel EM-Seq pipeline background PID: ${pid}" >&2
                PARALLEL_SUPPORT_PIDS+=("${pid}")
            fi
        fi
    fi
fi

if [[ "${SKIP_NEXTFLOW}" == "true" ]]; then
    echo "Skipping Nextflow run (--skip-nextflow active)." >&2
else
    attempt=1
    while :; do
        RESUME_FLAG=""
        if [[ "${INITIAL_RESUME}" == "true" || "${attempt}" -gt 1 ]]; then
            RESUME_FLAG="-resume"
        fi

        declare -a NF_ARGS=()
        NF_ARGS+=(run "${PIPELINE_DIR}/main.nf")
        NF_ARGS+=(-profile "${PROFILE}")
        NF_ARGS+=(-c "${CUSTOM_CONFIG}")
        for panel_cfg in "${PANEL_CONFIGS[@]}"; do
            NF_ARGS+=(-c "${panel_cfg}")
        done
        NF_ARGS+=(-w "${WORK_DIR}")
        NF_ARGS+=(--validate_params false)
        NF_ARGS+=(--mode "${MODE}")
        NF_ARGS+=(--genome "${GENOME}")
        NF_ARGS+=(--input "${SAMPLESHEET}")
        NF_ARGS+=(--outdir "${OUTDIR}")

        if [[ "${MODE}" == "targeted" && -n "${PANEL_NAME}" ]]; then
            NF_ARGS+=(--panel "${PANEL_NAME}")
        fi
        if [[ "${MODE}" == "targeted" && "${FORCE_PANEL}" == "true" ]]; then
            NF_ARGS+=(--force_panel)
        fi
        if [[ "${MODE}" == "purity_estimate" ]]; then
            NF_ARGS+=(--purity_estimate_mode "${PURITY_ESTIMATE_MODE:-wgts}")
        fi
        if [[ "${AMBER_EXPANDED_SNPS}" == "true" ]]; then
            NF_ARGS+=(--amber_expanded_snps)
        fi
        if [[ -n "${PROCESSES_EXCLUDE}" ]]; then
            NF_ARGS+=(--processes_exclude "${PROCESSES_EXCLUDE}")
        fi
        if [[ -n "${SEQUENCING_PLATFORM}" ]]; then
            NF_ARGS+=(--sequencing_platform "${SEQUENCING_PLATFORM}")
        fi

        if [[ "${MODE}" == "panel_resource_creation" ]]; then
            NF_ARGS+=(--driver_gene_panel "${DRIVER_GENE_PANEL}")
            NF_ARGS+=(--target_regions_bed "${TARGET_REGIONS_BED}")
            if [[ -n "${ISOFOX_GENE_IDS}" ]]; then
                NF_ARGS+=(--isofox_gene_ids "${ISOFOX_GENE_IDS}")
            fi
        fi

        NF_ARGS+=(-with-report "${OUTDIR}/pipeline_info/execution_report.html")
        NF_ARGS+=(-with-timeline "${OUTDIR}/pipeline_info/execution_timeline.html")
        NF_ARGS+=(-with-trace "${OUTDIR}/pipeline_info/execution_trace.txt")
        NF_ARGS+=(-with-dag "${OUTDIR}/pipeline_info/pipeline_dag.html")
        if [[ -n "${RESUME_FLAG}" ]]; then
            NF_ARGS+=("${RESUME_FLAG}")
        fi

        # Forward any extra positional arguments (e.g. from after '--') to Nextflow
        if [[ $# -gt 5 ]]; then
            NF_ARGS+=("${@:6}")
        fi

        echo "Nextflow attempt ${attempt}/${MAX_RESUME_ATTEMPTS} (resume_flag='${RESUME_FLAG:-none}')"
        if nextflow "${NF_ARGS[@]}" \
            2>&1 | tee -a "${OUTDIR}/nextflow.log"; then
            break
        fi

        if [[ "${AUTO_RESUME_ON_FAILURE}" != "true" ]]; then
            echo "Nextflow failed and AUTO_RESUME_ON_FAILURE is disabled; exiting." >&2
            exit 1
        fi
        if [[ "${attempt}" -ge "${MAX_RESUME_ATTEMPTS}" ]]; then
            echo "Nextflow failed after ${MAX_RESUME_ATTEMPTS} attempt(s); exiting." >&2
            exit 1
        fi
        echo "Nextflow failed on attempt ${attempt}; retrying with -resume in ${RESUME_BACKOFF_SECONDS}s..." >&2
        sleep "${RESUME_BACKOFF_SECONDS}"
        attempt=$((attempt + 1))
    done
fi

# ---------------------------------------------------------------------------
# Optional post-run off-target ichorCNA run or WGS ichorCNA run
# ---------------------------------------------------------------------------
if [[ "${ICHOR_WGS}" == "true" ]]; then
    (
        set +e
        set +o pipefail
        echo "Starting post-run WGS ichorCNA execution (Tumours, Normals, and Plasmas)..." >&2
    ORIGINAL_ICHOR_DIR="${ICHOR_DIR}"

    # Compute base patient directory (e.g. if OUTDIR is <patient>/plasma/<plasma_id>, PATIENT_BASEDIR is <patient>)
    PATIENT_BASEDIR="${OUTDIR}"
    if [[ "$(basename "$(dirname "${OUTDIR}")")" == "plasma" ]]; then
        PATIENT_BASEDIR="$(dirname "$(dirname "${OUTDIR}")")"
    elif [[ "$(basename "${OUTDIR}")" == "plasma" ]]; then
        PATIENT_BASEDIR="$(dirname "${OUTDIR}")"
    fi
    
    # Load required modules
    module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true
    module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R/4.5.2 2>/dev/null || module load r/4.3.1 2>/dev/null || module load R 2>/dev/null || module load r 2>/dev/null || true
    module load samtools/1.20 2>/dev/null || module load samtools/1.19 2>/dev/null || module load samtools 2>/dev/null || true

    # Resolve FASTA path for samtools view decoding if CRAM files are present
    if [[ "${GENOME}" == "GRCh37_hmf" ]]; then
        REF_FASTA="$(awk -F"'" '/ref_data_genome_fasta\s*=/{print $2; exit}' "${CUSTOM_CONFIG}")"
        GENOME_BUILD="hg19"
        CENT_FILE="GRCh37.p13_centromere_UCSC-gapTable.txt"
    else
        REF_FASTA="$(awk -F"'" '/ref_data_genome_fasta\s*=/{count++; if(count==2){print $2; exit}}' "${CUSTOM_CONFIG}")"
        GENOME_BUILD="hg38"
        CENT_FILE="GRCh38.GCA_000001405.2_centromere_acen.txt"
    fi

    # Discover and pair all DNA samples (Tumours, Normals, Plasmas)
    mapfile -t PAIRED_SAMPLES < <(python3 -c "
import csv, os, sys

samplesheet_path = sys.argv[1]
outdir = sys.argv[2]
mode = sys.argv[3] if len(sys.argv) > 3 else ''

samples = []
with open(samplesheet_path, 'r', encoding='utf-8-sig') as f:
    lines = [l for l in f if l.strip() and not l.strip().startswith('#')]
    reader = csv.DictReader(lines)
    for row in reader:
        samples.append(row)

outdir_bams = {}
sample_bams_outdir = {}

for root, dirs, files in os.walk(outdir):
    dirs[:] = [d for d in dirs if d not in ('work', 'logs', 'pipeline_info', '.nextflow', 'ichorCNA', 'ichorCNA_tumor_only', 'ichorCNA_normal_only', 'ichorCNA_t_vs_n', 'ichorCNA_plasma_vs_n', 'ichorCNA_plasma_only')]
    for f in files:
        if f.endswith(('.bam', '.cram')) and not f.endswith('.downsampled.bam') and not f.endswith('.downsampled_1x.bam'):
            full_path = os.path.join(root, f)
            outdir_bams[f] = full_path
            if 'alignments' in root and 'dna' in root:
                sid = f
                for ext in ('.redux.bam', '.redux.cram', '.bam', '.cram'):
                    if sid.endswith(ext):
                        sid = sid[:-len(ext)]
                        break
                sample_bams_outdir[sid] = full_path

def resolve_bam(sample_id, row):
    fp = (row.get('filepath') or '').strip()
    ft = (row.get('filetype') or '').strip().lower()
    if ft in ('bam', 'cram', 'bam_redux', 'cram_redux') and fp:
        for p in fp.split(';'):
            p = p.strip()
            if not p:
                continue
            if os.path.isfile(p):
                return p
            if p.endswith('.bam') and os.path.isfile(p[:-4] + '.cram'):
                return p[:-4] + '.cram'
            if p.endswith('.cram') and os.path.isfile(p[:-5] + '.bam'):
                return p[:-5] + '.bam'

    if sample_id in sample_bams_outdir and os.path.isfile(sample_bams_outdir[sample_id]):
        return sample_bams_outdir[sample_id]

    for candidate in (
        f'{sample_id}.redux.bam',
        f'{sample_id}.redux.cram',
        f'{sample_id}.bam',
        f'{sample_id}.cram',
    ):
        if candidate in outdir_bams and os.path.isfile(outdir_bams[candidate]):
            return outdir_bams[candidate]

    for fname, fpath in outdir_bams.items():
        if fname.startswith(sample_id) and fname.endswith(('.redux.bam', '.redux.cram', '.bam', '.cram')):
            if os.path.isfile(fpath):
                return fpath

    # Fallback search for older/alternative directory structures (e.g. BFIERCE V2, V3, bam_links)
    patient_hint = (row.get('subject_id') or row.get('group_id') or '').strip()
    if patient_hint:
        candidate_dirs = [
            f"/cluster/projects/rskimlab/projects/BFIERCE/fragmentomics/bam_links/{patient_hint}/plasma",
            f"/cluster/projects/rskimlab/projects/BFIERCE/fragmentomics/bam_links/{patient_hint}/tumour",
            f"/cluster/projects/rskimlab/projects/BFIERCE/fragmentomics/bam_links/{patient_hint}/normal",
        ]
        for cdir in candidate_dirs:
            for ext in ('.bam', '.cram', '.redux.bam', '.redux.cram'):
                cp = os.path.join(cdir, f"{sample_id}{ext}")
                if os.path.isfile(cp):
                    return cp

        if patient_hint.startswith('BFIERCE'):
            import glob
            v2_v3_patterns = [
                f"/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser/{patient_hint}/*/alignments/dna/{sample_id}.redux.bam",
                f"/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3/{patient_hint}/*/alignments/{sample_id}/{sample_id}.redux.bam",
                f"/cluster/projects/rskimlab/projects/BFIERCE/Oncoanalyser_v3/{patient_hint}/plasma/*/*/alignments/{sample_id}/{sample_id}.redux.bam",
            ]
            for pat in v2_v3_patterns:
                matches = glob.glob(pat)
                if matches and os.path.isfile(matches[0]):
                    return matches[0]

    return None

def get_sample_category(row):
    st = (row.get('sample_type') or '').strip().lower()
    info = (row.get('info') or '').strip().lower()
    if st == 'normal':
        return 'normal'
    if mode == 'purity_estimate' or st == 'plasma' or 'longitudinal_sample' in info or 'plasma' in info or 'is_plasma' in info:
        return 'plasma'
    return 'tumor'

dna_samples = []
for r in samples:
    st = (r.get('sequence_type') or '').strip().lower()
    ft = (r.get('filetype') or '').strip().lower()
    if ft in ('purple_dir', 'amber_dir'):
        continue
    if not st or st == 'dna':
        dna_samples.append(r)

seen_samples = set()
unique_rows = []
for r in dna_samples:
    sid = r.get('sample_id', '').strip()
    if sid and sid not in seen_samples:
        seen_samples.add(sid)
        unique_rows.append(r)

normals_by_subject = {}
for r in unique_rows:
    if get_sample_category(r) == 'normal':
        sid = r.get('sample_id', '').strip()
        sub = r.get('subject_id', '').strip()
        grp = r.get('group_id', '').strip()
        path = resolve_bam(sid, r)
        if path:
            if sub and sub not in normals_by_subject:
                normals_by_subject[sub] = (sid, path)
            if grp and grp not in normals_by_subject:
                normals_by_subject[grp] = (sid, path)

for r in unique_rows:
    sid = r.get('sample_id', '').strip()
    cat = get_sample_category(r)
    spath = resolve_bam(sid, r)
    if not spath:
        print(f'WARN: Could not locate BAM/CRAM for DNA sample {sid}', file=sys.stderr)
        continue
    sub = r.get('subject_id', '').strip()
    grp = r.get('group_id', '').strip()

    nid, npath = 'NONE', 'NONE'
    if cat in ('tumor', 'plasma'):
        if sub in normals_by_subject:
            nid, npath = normals_by_subject[sub]
        elif grp in normals_by_subject:
            nid, npath = normals_by_subject[grp]

    print(f'{sid}\t{cat}\t{spath}\t{nid}\t{npath}')
" "${SAMPLESHEET}" "${OUTDIR}" "${MODE}")

    if [[ ${#PAIRED_SAMPLES[@]} -eq 0 ]]; then
        echo "WARNING: No WGS DNA samples found with BAM/CRAM files in ${OUTDIR} or ${SAMPLESHEET}." >&2
    else
        echo "Found ${#PAIRED_SAMPLES[@]} sample(s) for WGS ichorCNA processing." >&2

        GC_WIG=$(Rscript -e "cat(system.file('extdata', 'gc_${GENOME_BUILD}_1000kb.wig', package='ichorCNA'))" 2>/dev/null || true)
        MAP_WIG=$(Rscript -e "cat(system.file('extdata', 'map_${GENOME_BUILD}_1000kb.wig', package='ichorCNA'))" 2>/dev/null || true)
        CENTROMERE=$(Rscript -e "cat(system.file('extdata', '${CENT_FILE}', package='ichorCNA'))" 2>/dev/null || true)

        # Resolve PON arguments
        PON_ARG=""
        if [[ -n "${ICHOR_PON}" ]]; then
            PON_ARG="--normalPanel ${ICHOR_PON}"
        else
            if [[ "${GENOME_BUILD}" == "hg38" ]]; then
                DEFAULT_PON_FILE="HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds"
            else
                DEFAULT_PON_FILE="HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
            fi
            DEFAULT_PON=$(Rscript -e "cat(system.file('extdata', '${DEFAULT_PON_FILE}', package='ichorCNA'))" 2>/dev/null || true)
            if [[ -n "${DEFAULT_PON}" && -f "${DEFAULT_PON}" ]]; then
                PON_ARG="--normalPanel ${DEFAULT_PON}"
            fi
        fi

        SC_STATES_VAL="NULL"
        if [[ "${ICHOR_SUB_CLONES}" == "true" ]]; then
            SC_STATES_VAL="c(1,3)"
        fi

        # Helper function to downsample BAM/CRAM to 1.0X and run readCounter
        generate_sample_wig() {
            local s_id="$1"
            local s_path="$2"
            local s_dir="$3"
            local t_wig="$4"
            local s_chrs="$5"
            local s_cat="${6:-tumor}"
            local wig_type="${7:-standard}" # "standard" or "short"

            if [[ -f "${t_wig}" && -s "${t_wig}" && "${FORCE_ICHOR}" != "true" ]]; then
                local non_zero_counts
                non_zero_counts=$(grep -v "^#" "${t_wig}" 2>/dev/null | grep -v "^fixedStep" | awk '$1 > 0 {found=1; exit} END {print (found==1) ? "true" : "false"}')
                if [[ "${non_zero_counts}" == "true" ]]; then
                    echo "    WIG file already exists for ${s_id} (${wig_type}): ${t_wig}." >&2
                    return 0
                else
                    echo "    Existing WIG file ${t_wig} has 0 counts (from failed previous attempt). Regenerating..." >&2
                fi
            fi

            mkdir -p "${s_dir}"
            local total_reads
            total_reads=$( (samtools idxstats "${s_path}" 2>/dev/null || true) | awk '{sum+=$3} END {print sum}')
            local cov=0
            if [[ "${total_reads}" =~ ^[0-9]+$ ]] && [[ "${total_reads}" -gt 0 ]]; then
                cov=$(awk -v reads="${total_reads}" 'BEGIN {printf "%.4f", (reads * 150) / 3100000000}')
            fi
            echo "    Sample ${s_id} estimated coverage: ${cov}x (reads: ${total_reads})" >&2

            local downsample_bam="${s_path}"
            local is_temp="false"
            local nthreads="${SLURM_CPUS_PER_TASK:-4}"
            if [[ "${nthreads}" -lt 4 ]]; then
                nthreads=4
            fi
            local temp_bam="${TMPDIR:-$s_dir}/${s_id}_${wig_type}_$$.downsampled_1x.bam"
            mkdir -p "$(dirname "${temp_bam}")"

            if [[ "${wig_type}" == "short" ]]; then
                local final_short_bam="${s_dir}/${s_id}.short_selected.bam"
                if [[ -f "${final_short_bam}" && -s "${final_short_bam}" && -f "${final_short_bam}.bai" ]]; then
                    echo "    Found existing short fragment size-selected BAM: ${final_short_bam}. Reusing." >&2
                    downsample_bam="${final_short_bam}"
                    is_temp="false"
                else
                    local min_l="${ICHOR_MIN_TLEN:-90}"
                    local max_l="${ICHOR_MAX_TLEN:-150}"
                    echo "    Applying in-silico fragment size selection (${min_l}-${max_l} bp) for ${s_id}..." >&2

                    # Leverage Oncoanalyser WISP output if available to instantly compute short fragment proportion
                    local short_pct=""
                    local wisp_frag_file
                    wisp_frag_file=$(find "${OUTDIR}" "${PATIENT_BASEDIR}" -maxdepth 5 -type f -name "*${s_id}*.wisp.fragment_length.tsv" 2>/dev/null | head -n 1)
                    if [[ -z "${wisp_frag_file}" ]]; then
                        wisp_frag_file=$(find "${OUTDIR}" "${PATIENT_BASEDIR}" -maxdepth 5 -type f -name "*.wisp.fragment_length.tsv" 2>/dev/null | head -n 1)
                    fi

                    if [[ -n "${wisp_frag_file}" && -f "${wisp_frag_file}" ]]; then
                        short_pct=$(awk -v min_l="${min_l}" -v max_l="${max_l}" 'NR>1 {tot+=$2+$3; if($1>=min_l && $1<=max_l) short+=$2+$3} END {if(tot>0) printf "%.6f", short/tot}' "${wisp_frag_file}")
                    fi

                    if [[ -z "${short_pct}" || "${short_pct}" == "0.000000" ]]; then
                        short_pct="0.200000"
                    fi

                    local short_reads
                    short_reads=$(awk -v t="${total_reads}" -v p="${short_pct}" 'BEGIN {printf "%.0f", t * p}')
                    local short_cov
                    short_cov=$(awk -v r="${short_reads}" 'BEGIN {printf "%.4f", (r * 120) / 3100000000}')
                    echo "    ${s_id} estimated short fragment pool: ${short_reads} reads (~${short_cov}x, $(awk -v p="${short_pct}" 'BEGIN {printf "%.1f%%", p*100}'))" >&2

                    local target_cov_short="${ICHOR_SHORT_TARGET_COV:-10.0}"
                    local needs_downsample_short
                    needs_downsample_short=$(awk -v c="${short_cov}" -v tc="${target_cov_short}" 'BEGIN {print (c > tc) ? "true" : "false"}')
                    if [[ "${needs_downsample_short}" == "true" ]]; then
                        local fraction
                        fraction=$(awk -v c="${short_cov}" -v tc="${target_cov_short}" 'BEGIN {printf "%.6f", tc / c}')
                        local s_param
                        s_param=$(awk -v f="${fraction}" 'BEGIN {printf "%.6f", 100 + f}')
                        echo "    Streaming fast pre-downsampling (-s ${s_param}) + size selection (${min_l}-${max_l} bp) to ~${target_cov_short}x..." >&2

                        # Pre-subsample with -s in first pipe so samtools only evaluates filter expression on sampled reads (10x faster!)
                        samtools view -@ "${nthreads}" -T "${REF_FASTA}" -u -s "${s_param}" "${s_path}" \
                            | samtools view -@ "${nthreads}" -b \
                                -e "(tlen >= ${min_l} && tlen <= ${max_l}) || (tlen <= -${min_l} && tlen >= -${max_l})" - \
                            > "${temp_bam}"
                    else
                        echo "    Short fragment coverage <= ${target_cov_short}x (~${short_cov}x); streaming size selection with 100% retention..." >&2
                        samtools view -@ "${nthreads}" -T "${REF_FASTA}" -b \
                            -e "(tlen >= ${min_l} && tlen <= ${max_l}) || (tlen <= -${min_l} && tlen >= -${max_l})" "${s_path}" \
                            > "${temp_bam}"
                    fi

                    if [[ ! -s "${temp_bam}" ]]; then
                        echo "    WARNING: Size selection produced empty BAM for ${s_id}; falling back to standard downsampling without size filter..." >&2
                        samtools view -@ "${nthreads}" -T "${REF_FASTA}" -b -s 100.05 "${s_path}" > "${temp_bam}"
                    fi

                    samtools index -@ "${nthreads}" "${temp_bam}"
                    if [[ "${KEEP_DOWNSAMPLED_BAMS:-true}" == "true" ]]; then
                        mv -f "${temp_bam}" "${final_short_bam}"
                        mv -f "${temp_bam}.bai" "${final_short_bam}.bai"
                        downsample_bam="${final_short_bam}"
                        is_temp="false"
                        echo "    Permanently retained short-fragment BAM: ${final_short_bam}" >&2
                    else
                        downsample_bam="${temp_bam}"
                        is_temp="true"
                    fi
                fi
            else
                local final_std_bam="${s_dir}/${s_id}.downsampled_1x.bam"
                if [[ -f "${final_std_bam}" && -s "${final_std_bam}" && -f "${final_std_bam}.bai" ]]; then
                    echo "    Found existing 1.0x downsampled BAM: ${final_std_bam}. Reusing." >&2
                    downsample_bam="${final_std_bam}"
                    is_temp="false"
                else
                    local needs_downsample
                    needs_downsample=$(awk -v c="${cov}" 'BEGIN {print (c > 1.0) ? "true" : "false"}')

                    if [[ "${needs_downsample}" == "true" ]]; then
                        local fraction
                        fraction=$(awk -v c="${cov}" 'BEGIN {printf "%.6f", 1.0 / c}')
                        local s_param
                        s_param=$(awk -v f="${fraction}" 'BEGIN {printf "%.6f", 100 + f}')
                        echo "    Downsampling ${s_id} reads to ~1.0x (fraction: ${fraction}, -s ${s_param}) using ${nthreads} threads..." >&2

                        samtools view -@ "${nthreads}" -T "${REF_FASTA}" -b -s "${s_param}" "${s_path}" > "${temp_bam}"
                        samtools index -@ "${nthreads}" "${temp_bam}"
                        if [[ "${KEEP_DOWNSAMPLED_BAMS:-true}" == "true" ]]; then
                            mv -f "${temp_bam}" "${final_std_bam}"
                            mv -f "${temp_bam}.bai" "${final_std_bam}.bai"
                            downsample_bam="${final_std_bam}"
                            is_temp="false"
                            echo "    Permanently retained downsampled 1.0x BAM: ${final_std_bam}" >&2
                        else
                            downsample_bam="${temp_bam}"
                            is_temp="true"
                        fi
                    fi
                fi
            fi

            echo "    Running readCounter for ${s_id} (${wig_type})..." >&2
            readCounter --window 1000000 --quality 20 --chromosome "${s_chrs}" "${downsample_bam}" > "${t_wig}"

            if [[ "${is_temp}" == "true" && -f "${temp_bam}" ]]; then
                rm -f "${temp_bam}" "${temp_bam}.bai"
            fi
        }

        # Process each sample
        for row in "${PAIRED_SAMPLES[@]}"; do
            sample_id=$(echo "$row" | cut -f1)
            sample_cat=$(echo "$row" | cut -f2)
            sample_path=$(echo "$row" | cut -f3)
            normal_id=$(echo "$row" | cut -f4)
            normal_path=$(echo "$row" | cut -f5)

            # In purity_estimate mode, skip profiling the matched normal again (already profiled during Step 1 WGTS)
            if [[ "${MODE}" == "purity_estimate" && "${sample_cat}" == "normal" ]]; then
                echo "  [SKIPPING NORMAL] Normal ${sample_id} in purity_estimate mode already processed during baseline WGTS." >&2
                continue
            fi

            echo "--------------------------------------------------------" >&2
            echo "  Processing [${sample_cat^^}]: ${sample_id} (Matched Normal: ${normal_id})" >&2

            # Determine target directories for paired and unpaired analysis under PATIENT_BASEDIR
            PAIRED_DIR=""
            UNPAIRED_DIR=""
            if [[ "${sample_cat}" == "normal" ]]; then
                UNPAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_normal_only/${sample_id}"
            elif [[ "${sample_cat}" == "tumor" ]]; then
                if [[ "${normal_path}" != "NONE" ]]; then
                    PAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_t_vs_n/${sample_id}"
                    UNPAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_tumor_only/${sample_id}"
                else
                    UNPAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_tumor_only/${sample_id}"
                fi
            elif [[ "${sample_cat}" == "plasma" ]]; then
                if [[ "${normal_path}" != "NONE" ]]; then
                    PAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n/${sample_id}"
                    UNPAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_plasma_only/${sample_id}"
                else
                    UNPAIRED_DIR="${PATIENT_BASEDIR}/ichorCNA_plasma_only/${sample_id}"
                fi
            fi

            # Check if all runs for this sample are already done
            already_done="false"
            if [[ "${FORCE_ICHOR}" == "true" ]]; then
                already_done="false"
            elif [[ "${sample_cat}" == "normal" ]]; then
                if [[ -s "${UNPAIRED_DIR}/${sample_id}.params.txt" ]]; then
                    already_done="true"
                fi
            elif [[ "${sample_cat}" == "plasma" ]]; then
                # For plasma, both standard and short-fragment solutions must exist
                local has_std_u="false"
                local has_sh_u="false"
                [[ -s "${UNPAIRED_DIR}/standard/${sample_id}.params.txt" || -s "${UNPAIRED_DIR}/${sample_id}.standard.params.txt" ]] && has_std_u="true"
                [[ -s "${UNPAIRED_DIR}/short/${sample_id}.params.txt" || -s "${UNPAIRED_DIR}/${sample_id}.short.params.txt" ]] && has_sh_u="true"
                if [[ "${has_std_u}" == "true" && "${has_sh_u}" == "true" ]]; then
                    already_done="true"
                fi
            else
                if [[ -n "${PAIRED_DIR}" && -n "${UNPAIRED_DIR}" ]]; then
                    if [[ -s "${PAIRED_DIR}/${sample_id}.params.txt" && -s "${UNPAIRED_DIR}/${sample_id}.params.txt" ]]; then
                        already_done="true"
                    fi
                elif [[ -n "${UNPAIRED_DIR}" && -s "${UNPAIRED_DIR}/${sample_id}.params.txt" ]]; then
                    already_done="true"
                elif [[ -n "${PAIRED_DIR}" && -s "${PAIRED_DIR}/${sample_id}.params.txt" ]]; then
                    already_done="true"
                fi
            fi

            if [[ "${already_done}" == "true" ]]; then
                echo "    ichorCNA results (standard & short) already exist for ${sample_cat} ${sample_id}. Skipping execution." >&2
                continue
            fi

            # Auto-detect chromosome naming style from BAM/CRAM idxstats
            first_chr=$( (samtools idxstats "${sample_path}" 2>/dev/null || true) | awk 'NR==1 {print $1; exit}')
            if [[ "${first_chr}" =~ ^chr ]]; then
                CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
                GENOME_STYLE="UCSC"
                CHR_NORM="paste0('chr', 1:22)"
                CHR_TRAIN="paste0('chr', 1:22)"
                CHR_ANALYZE="c(paste0('chr', 1:22),'chrX')"
            else
                CHRS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
                GENOME_STYLE="NCBI"
                CHR_NORM="c(1:22)"
                CHR_TRAIN="c(1:22)"
                CHR_ANALYZE="c(1:22,'X')"
            fi

            # Primary work directory for this sample's WIG
            PRIMARY_DIR="${PAIRED_DIR:-${UNPAIRED_DIR}}"
            mkdir -p "${PRIMARY_DIR}"

            # Helper function to execute Rscript runIchorCNA.R in a directory
            execute_ichor_run() {
                local run_dir="$1"
                local in_wig="$2"
                local n_wig="$3"
                local run_label="$4"

                mkdir -p "${run_dir}"
                if [[ -s "${run_dir}/${sample_id}.params.txt" && "${FORCE_ICHOR}" != "true" ]]; then
                    echo "    ichorCNA ${run_label} results already exist in ${run_dir}. Skipping." >&2
                    return 0
                fi

                echo "    Running ichorCNA (${run_label}) for ${sample_id} in $(basename "${run_dir%/*}")/$(basename "${run_dir}")..." >&2
                local norm_arg=""
                if [[ -n "${n_wig}" && -f "${n_wig}" ]]; then
                    norm_arg="--NORMWIG ${n_wig}"
                fi

                (
                    cd "${run_dir}"
                    Rscript "${SOURCE_REPO}/resources/runIchorCNA.R" \
                        --WIG "${in_wig}" \
                        ${norm_arg} \
                        --id "${sample_id}" \
                        --genomeStyle "${GENOME_STYLE}" \
                        --genomeBuild "${GENOME_BUILD}" \
                        --gcWig "${GC_WIG}" \
                        --mapWig "${MAP_WIG}" \
                        --centromere "${CENTROMERE}" \
                        ${PON_ARG} \
                        --chrNormalize "${CHR_NORM}" \
                        --chrTrain "${CHR_TRAIN}" \
                        --chrs "${CHR_ANALYZE}" \
                        --normal "c(0.5,0.6,0.7,0.8,0.9)" \
                        --ploidy "c(2,3)" \
                        --scStates "${SC_STATES_VAL}"
                )
                # Flatten any nested sample subfolder created by plotSolutions
                if [[ -d "${run_dir}/${sample_id}" ]]; then
                    mv -f "${run_dir}/${sample_id}"/* "${run_dir}/" 2>/dev/null || true
                    rmdir "${run_dir}/${sample_id}" 2>/dev/null || true
                fi
            }

            # Generate normal WIG if matched normal exists
            normal_wig=""
            if [[ "${normal_path}" != "NONE" && -f "${normal_path}" ]]; then
                candidate_wigs=(
                    "${PATIENT_BASEDIR}/ichorCNA_normal_only/${normal_id}/${normal_id}.wig"
                    "${PATIENT_BASEDIR}/ichorCNA/${normal_id}/${normal_id}.wig"
                    "${OUTDIR}/../../ichorCNA_normal_only/${normal_id}/${normal_id}.wig"
                    "${OUTDIR}/../../ichorCNA/${normal_id}/${normal_id}.wig"
                    "${OUTDIR}/../ichorCNA/${normal_id}/${normal_id}.wig"
                )
                for cw in "${candidate_wigs[@]}"; do
                    if [[ -f "${cw}" && -s "${cw}" ]]; then
                        normal_wig="${cw}"
                        echo "    Found existing normal WIG: ${normal_wig}" >&2
                        break
                    fi
                done
                if [[ -z "${normal_wig}" ]]; then
                    NORMAL_ICHOR_DIR="${PATIENT_BASEDIR}/ichorCNA_normal_only/${normal_id}"
                    mkdir -p "${NORMAL_ICHOR_DIR}"
                    normal_wig="${NORMAL_ICHOR_DIR}/${normal_id}.wig"
                    generate_sample_wig "${normal_id}" "${normal_path}" "${NORMAL_ICHOR_DIR}" "${normal_wig}" "${CHRS}" "normal" "standard"
                fi
            fi

            if [[ "${sample_cat}" == "normal" ]]; then
                # --- Normal Run with PoN ---
                mkdir -p "${UNPAIRED_DIR}"
                sample_wig="${PRIMARY_DIR}/${sample_id}.wig"
                generate_sample_wig "${sample_id}" "${sample_path}" "${PRIMARY_DIR}" "${sample_wig}" "${CHRS}" "normal" "standard"
                execute_ichor_run "${UNPAIRED_DIR}" "${sample_wig}" "" "normal-PoN"
                mkdir -p "${PATIENT_BASEDIR}/ichorCNA"
                ln -sfn "../ichorCNA_normal_only/${sample_id}" "${PATIENT_BASEDIR}/ichorCNA/${sample_id}" 2>/dev/null || true
            elif [[ "${sample_cat}" == "plasma" ]]; then
                # --- Plasma ctDNA Run: Dual Standard 1.0X AND Short-Fragment 90-150 bp 1.0X ---
                sample_wig_std="${PRIMARY_DIR}/${sample_id}.standard.wig"
                sample_wig_short="${PRIMARY_DIR}/${sample_id}.short.wig"

                generate_sample_wig "${sample_id}" "${sample_path}" "${PRIMARY_DIR}" "${sample_wig_std}" "${CHRS}" "plasma" "standard"
                generate_sample_wig "${sample_id}" "${sample_path}" "${PRIMARY_DIR}" "${sample_wig_short}" "${CHRS}" "plasma" "short"

                # Process target directories (PoN unpaired and matched-normal paired)
                for target_base_dir in "${UNPAIRED_DIR}" "${PAIRED_DIR}"; do
                    [[ -z "${target_base_dir}" ]] && continue
                    is_paired_run="false"
                    active_norm_wig=""
                    if [[ "${target_base_dir}" == "${PAIRED_DIR}" && -n "${normal_wig}" ]]; then
                        is_paired_run="true"
                        active_norm_wig="${normal_wig}"
                    fi

                    mkdir -p "${target_base_dir}/standard" "${target_base_dir}/short"

                    # 1. Standard downsampled run
                    execute_ichor_run "${target_base_dir}/standard" "${sample_wig_std}" "${active_norm_wig}" "plasma-standard"

                    # 2. Short-fragment size-selected run
                    execute_ichor_run "${target_base_dir}/short" "${sample_wig_short}" "${active_norm_wig}" "plasma-short"

                    # 3. Create suffixed convenience links in directory root
                    (
                        cd "${target_base_dir}"
                        ln -sfn "standard/${sample_id}.params.txt" "${sample_id}.standard.params.txt" 2>/dev/null || true
                        ln -sfn "standard/${sample_id}_genomeWide.pdf" "${sample_id}_standard_genomeWide.pdf" 2>/dev/null || true
                        ln -sfn "standard/${sample_id}.seg.txt" "${sample_id}.standard.seg.txt" 2>/dev/null || true

                        ln -sfn "short/${sample_id}.params.txt" "${sample_id}.short.params.txt" 2>/dev/null || true
                        ln -sfn "short/${sample_id}_genomeWide.pdf" "${sample_id}_short_genomeWide.pdf" 2>/dev/null || true
                        ln -sfn "short/${sample_id}.seg.txt" "${sample_id}.short.seg.txt" 2>/dev/null || true

                        # Primary default files for downstream tools/tables: SHORT READ IS DEFAULT FOR PLASMA!
                        ln -sfn "short/${sample_id}.params.txt" "${sample_id}.params.txt" 2>/dev/null || true
                        ln -sfn "short/${sample_id}_genomeWide.pdf" "${sample_id}_genomeWide.pdf" 2>/dev/null || true
                        ln -sfn "short/${sample_id}.seg.txt" "${sample_id}.seg.txt" 2>/dev/null || true
                        ln -sfn "short/${sample_id}.RData" "${sample_id}.RData" 2>/dev/null || true
                    )
                done

                # Ensure PoN results are linked in PATIENT_BASEDIR/ichorCNA
                mkdir -p "${PATIENT_BASEDIR}/ichorCNA"
                rel_dir="../$(basename "${UNPAIRED_DIR%/*}")/${sample_id}"
                ln -sfn "${rel_dir}" "${PATIENT_BASEDIR}/ichorCNA/${sample_id}" 2>/dev/null || true

                # Clean up sample WIG files
                rm -f "${sample_wig_std}" "${sample_wig_short}" 2>/dev/null || true

                # --- DELFI-TF Feature Extraction & Inference ---
                if [[ "${ENABLE_DELFI}" == "true" ]]; then
                    DELFI_BIN=""
                    if [[ -n "${DELFI_REPO}" && -x "${DELFI_REPO}/bin/delfi_tf" ]]; then
                        DELFI_BIN="${DELFI_REPO}/bin/delfi_tf"
                    elif [[ -x "/cluster/projects/pughlab/src/DELFI-TF/bin/delfi_tf" ]]; then
                        DELFI_BIN="/cluster/projects/pughlab/src/DELFI-TF/bin/delfi_tf"
                    elif command -v delfi_tf >/dev/null 2>&1; then
                        DELFI_BIN=$(command -v delfi_tf)
                    fi

                    if [[ -n "${DELFI_BIN}" ]]; then
                        DELFI_DIR="${PATIENT_BASEDIR}/DELFI/${sample_id}"
                        mkdir -p "${DELFI_DIR}/logs"

                        # Prefer downsampled 1.0x BAM to run ultra-fast with minimal RAM
                        delfi_bam=""
                        if [[ -f "${PRIMARY_DIR}/${sample_id}.downsampled_1x.bam" && -s "${PRIMARY_DIR}/${sample_id}.downsampled_1x.bam" ]]; then
                            delfi_bam="${PRIMARY_DIR}/${sample_id}.downsampled_1x.bam"
                            echo "    [DELFI-TF] Using downsampled 1.0x BAM: ${delfi_bam}" >&2
                        else
                            delfi_bam="${sample_path}"
                            echo "    [DELFI-TF] Using primary BAM with streaming: ${delfi_bam}" >&2
                        fi

                        echo "    [DELFI-TF] Running DELFI-TF for ${sample_id}..." >&2
                        delfi_model_arg=""
                        if [[ -n "${DELFI_MODEL}" && -f "${DELFI_MODEL}" ]]; then
                            delfi_model_arg="--model ${DELFI_MODEL}"
                        fi

                        if [[ -n "${delfi_model_arg}" ]]; then
                            bash "${DELFI_BIN}" predict \
                                --bam "${delfi_bam}" \
                                --id "${sample_id}" \
                                ${delfi_model_arg} \
                                --max-reads "${DELFI_MAX_READS:-50000000}" \
                                --outdir "${DELFI_DIR}" \
                                > "${DELFI_DIR}/logs/delfi_${sample_id}.log" 2>&1 || \
                                echo "    WARNING: DELFI-TF predict encountered an issue; check ${DELFI_DIR}/logs/delfi_${sample_id}.log" >&2
                        else
                            bash "${DELFI_BIN}" extract \
                                --bam "${delfi_bam}" \
                                --id "${sample_id}" \
                                --max-reads "${DELFI_MAX_READS:-50000000}" \
                                --outdir "${DELFI_DIR}" \
                                > "${DELFI_DIR}/logs/delfi_${sample_id}.log" 2>&1 || \
                                echo "    WARNING: DELFI-TF extract encountered an issue; check ${DELFI_DIR}/logs/delfi_${sample_id}.log" >&2
                        fi

                        if [[ -f "${DELFI_DIR}/${sample_id}_delfi_tf_features.tsv" ]]; then
                            echo "    [DELFI-TF] Features generated successfully: ${DELFI_DIR}/${sample_id}_delfi_tf_features.tsv" >&2
                        fi
                    else
                        echo "    [DELFI-TF] delfi_tf executable not found. Set --delfi-repo to run DELFI-TF." >&2
                    fi
                fi
            else
                # --- Primary Tumor Run (Standard 1.0X) ---
                sample_wig="${PRIMARY_DIR}/${sample_id}.wig"
                generate_sample_wig "${sample_id}" "${sample_path}" "${PRIMARY_DIR}" "${sample_wig}" "${CHRS}" "tumor" "standard"

                for target_base_dir in "${UNPAIRED_DIR}" "${PAIRED_DIR}"; do
                    [[ -z "${target_base_dir}" ]] && continue
                    is_paired_run="false"
                    active_norm_wig=""
                    if [[ "${target_base_dir}" == "${PAIRED_DIR}" && -n "${normal_wig}" ]]; then
                        is_paired_run="true"
                        active_norm_wig="${normal_wig}"
                    fi

                    execute_ichor_run "${target_base_dir}" "${sample_wig}" "${active_norm_wig}" "tumor"
                done

                mkdir -p "${PATIENT_BASEDIR}/ichorCNA"
                rel_dir="../$(basename "${UNPAIRED_DIR%/*}")/${sample_id}"
                ln -sfn "${rel_dir}" "${PATIENT_BASEDIR}/ichorCNA/${sample_id}" 2>/dev/null || true
                rm -f "${sample_wig}" 2>/dev/null || true
            fi
        done

        # Point ICHOR_DIR to primary output folder so downstream summarize_findings.py uses it (defaulting to PoN)
        if [[ -d "${PATIENT_BASEDIR}/ichorCNA_tumor_only" ]]; then
            ICHOR_DIR="${PATIENT_BASEDIR}/ichorCNA_tumor_only"
        elif [[ -d "${PATIENT_BASEDIR}/ichorCNA_plasma_only" ]]; then
            ICHOR_DIR="${PATIENT_BASEDIR}/ichorCNA_plasma_only"
        elif [[ -d "${PATIENT_BASEDIR}/ichorCNA_t_vs_n" ]]; then
            ICHOR_DIR="${PATIENT_BASEDIR}/ichorCNA_t_vs_n"
        elif [[ -d "${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n" ]]; then
            ICHOR_DIR="${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n"
        else
            ICHOR_DIR="${OUTDIR}/ichorCNA"
        fi

        echo "Completed WGS ichorCNA execution. Sourced to summarize_findings: ${ICHOR_DIR}" >&2

        # Generate comparison TSV (preferring PoN directory)
        comp_dir="${PATIENT_BASEDIR}/ichorCNA_tumor_only"
        if [[ ! -d "${comp_dir}" ]]; then
            if [[ -d "${PATIENT_BASEDIR}/ichorCNA_plasma_only" ]]; then
                comp_dir="${PATIENT_BASEDIR}/ichorCNA_plasma_only"
            elif [[ -d "${PATIENT_BASEDIR}/ichorCNA_t_vs_n" ]]; then
                comp_dir="${PATIENT_BASEDIR}/ichorCNA_t_vs_n"
            elif [[ -d "${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n" ]]; then
                comp_dir="${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n"
            else
                comp_dir="${OUTDIR}/ichorCNA"
            fi
        fi
        mkdir -p "${comp_dir}"
        comparison_tsv="${comp_dir}/tumor_fraction_comparison.tsv"
        echo -e "sample_id\tsample_type\tpanel_tumor_fraction\tprovided_tumor_fraction\thmf_tumor_fraction\tdelfi_tumor_fraction\tdelfi_pa_score" > "${comparison_tsv}"
        for row in "${PAIRED_SAMPLES[@]}"; do
            sample_id=$(echo "$row" | cut -f1)
            sample_cat=$(echo "$row" | cut -f2)
            patient_id=$(echo "${sample_id}" | grep -oE "[A-Za-z0-9_-]+" | head -n 1 || echo "${sample_id}")

            # Extract WGS ichorCNA TF from candidate folders
            panel_tf="NA"
            s_cand_dirs=(
                "${PATIENT_BASEDIR}/ichorCNA_tumor_only/${sample_id}"
                "${PATIENT_BASEDIR}/ichorCNA_plasma_only/${sample_id}"
                "${PATIENT_BASEDIR}/ichorCNA_t_vs_n/${sample_id}"
                "${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n/${sample_id}"
                "${PATIENT_BASEDIR}/ichorCNA_normal_only/${sample_id}"
                "${OUTDIR}/ichorCNA/${sample_id}"
            )
            for scd in "${s_cand_dirs[@]}"; do
                if [[ -d "${scd}" ]]; then
                    panel_params=$(find "${scd}" -maxdepth 2 -type f \( -name "*params.txt" -o -name "*metrics.txt" \) 2>/dev/null | head -n 1)
                    if [[ -n "${panel_params}" && -f "${panel_params}" ]]; then
                        panel_tf=$(grep -i "tumor fraction" "${panel_params}" 2>/dev/null | awk -F: '{print $2}' | xargs 2>/dev/null || echo "NA")
                        if [[ "${panel_tf}" != "NA" && -n "${panel_tf}" ]]; then
                            break
                        fi
                    fi
                fi
            done

            # Extract provided TF
            provided_tf="NA"
            if [[ -n "${ORIGINAL_ICHOR_DIR}" && -d "${ORIGINAL_ICHOR_DIR}" ]]; then
                provided_params=$(find "${ORIGINAL_ICHOR_DIR}" -maxdepth 3 -type f \( -name "*${sample_id}*params.txt" -o -name "*${patient_id}*params.txt" \) 2>/dev/null | head -n 1)
                if [[ -z "${provided_params}" ]]; then
                    provided_params=$(find "${ORIGINAL_ICHOR_DIR}" -maxdepth 3 -type f \( -name "*params.txt" -o -name "*metrics.txt" \) 2>/dev/null | head -n 1)
                fi
                if [[ -n "${provided_params}" && -f "${provided_params}" ]]; then
                    provided_tf=$(grep -i "tumor fraction" "${provided_params}" 2>/dev/null | awk -F: '{print $2}' | xargs 2>/dev/null || echo "NA")
                fi
            fi

            # Extract HMF pipeline tumor fraction (PURPLE purity or WISP purity)
            hmf_tf="NA"
            purple_purity_file=$(find "${OUTDIR}" -maxdepth 4 -type f -name "*${sample_id}*.purple.purity.tsv" ! -path "*/work/*" ! -path "*/.nextflow/*" 2>/dev/null | head -n 1)
            if [[ -z "${purple_purity_file}" ]]; then
                purple_purity_file=$(find "${OUTDIR}" -maxdepth 4 -type f -name "*.purple.purity.tsv" ! -path "*/work/*" ! -path "*/.nextflow/*" 2>/dev/null | head -n 1)
            fi
            if [[ -n "${purple_purity_file}" && -f "${purple_purity_file}" ]]; then
                hmf_tf=$(awk 'NR==2 {print $1}' "${purple_purity_file}" 2>/dev/null || echo "NA")
            fi
            if [[ "${hmf_tf}" == "NA" ]]; then
                wisp_file=$(find "${OUTDIR}" -maxdepth 4 -type f -name "*${sample_id}*.wisp.summary.tsv" ! -path "*/work/*" ! -path "*/.nextflow/*" 2>/dev/null | head -n 1)
                if [[ -z "${wisp_file}" ]]; then
                    wisp_file=$(find "${OUTDIR}" -maxdepth 4 -type f -name "*.wisp.summary.tsv" ! -path "*/work/*" ! -path "*/.nextflow/*" 2>/dev/null | head -n 1)
                fi
                if [[ -n "${wisp_file}" && -f "${wisp_file}" ]]; then
                    hmf_tf=$(awk -F'\t' 'NR==2 {print $1}' "${wisp_file}" 2>/dev/null || echo "NA")
                fi
            fi

            # Extract DELFI-TF estimated tumor fraction and PA-score
            delfi_tf="NA"
            delfi_pa="NA"
            delfi_pred="${PATIENT_BASEDIR}/DELFI/${sample_id}/${sample_id}_delfi_tf_prediction.tsv"
            delfi_feat="${PATIENT_BASEDIR}/DELFI/${sample_id}/${sample_id}_delfi_tf_features.tsv"
            if [[ -f "${delfi_pred}" ]]; then
                delfi_tf=$(awk 'NR==2 {print $2}' "${delfi_pred}" 2>/dev/null || echo "NA")
                delfi_pa=$(awk 'NR==2 {print $4}' "${delfi_pred}" 2>/dev/null || echo "NA")
            elif [[ -f "${delfi_feat}" ]]; then
                delfi_pa=$(awk -F'\t' 'NR==2 {print $2}' "${delfi_feat}" 2>/dev/null || echo "NA")
            fi

            echo -e "${sample_id}\t${sample_cat}\t${panel_tf}\t${provided_tf}\t${hmf_tf}\t${delfi_tf}\t${delfi_pa}" >> "${comparison_tsv}"
        done
        echo "Saved tumor fraction comparison TSV to: ${comparison_tsv}" >&2
        mkdir -p "${PATIENT_BASEDIR}/ichorCNA"
        cp -f "${comparison_tsv}" "${PATIENT_BASEDIR}/ichorCNA/tumor_fraction_comparison.tsv" 2>/dev/null || true
        if [[ "${OUTDIR}" != "${PATIENT_BASEDIR}" ]]; then
            cp -f "${comparison_tsv}" "${OUTDIR}/tumor_fraction_comparison.tsv" 2>/dev/null || true
        fi
    fi
    ) || echo "WARNING: WGS ichorCNA post-processing encountered an error, but Nextflow outputs completed successfully." >&2
elif [[ -n "${ICHOR_EXCLUSION_BED}" && -n "${ICHOR_PON}" ]]; then
    echo "Starting post-run off-target ichorCNA execution..." >&2
    ORIGINAL_ICHOR_DIR="${ICHOR_DIR}"
    
    # Load required modules
    module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true
    module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R 2>/dev/null || true
    module load samtools/1.19 2>/dev/null || module load samtools 2>/dev/null || true
    
    # Locate all *.redux.bam files under the output directory (excluding work, logs, etc)
    COHORT_BAMS_FILE="${OUTDIR}/ichor_cohort_bams.txt"
    find "${OUTDIR}" -type f -name "*.redux.bam" \
        ! -path "${OUTDIR}/work/*" \
        ! -path "${OUTDIR}/pipeline_info/*" \
        ! -path "${OUTDIR}/logs/*" | while read -r filepath; do
            filename=$(basename "${filepath}")
            sample_id="${filename%.redux.bam}"
            echo -e "${sample_id}\t${filepath}"
        done > "${COHORT_BAMS_FILE}"

    if [[ ! -s "${COHORT_BAMS_FILE}" ]]; then
        echo "WARNING: No *.redux.bam files found in ${OUTDIR} to run ichorCNA." >&2
        rm -f "${COHORT_BAMS_FILE}"
    else
        # Pair tumor and normal samples using Python script
        PY_SCRIPT="${OUTDIR}/pair_samples.py"
        cat > "${PY_SCRIPT}" <<'EOF'
import csv
import sys

samplesheet = sys.argv[1]
bam_map_file = sys.argv[2]

bam_map = {}
with open(bam_map_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            bam_map[parts[0]] = parts[1]

samples = []
with open(samplesheet, 'r') as f:
    lines = [line for line in f if not line.strip().startswith('#')]
    reader = csv.DictReader(lines)
    for row in reader:
        samples.append(row)

for row in samples:
    if row.get('sample_type', '').lower() == 'tumor' and row.get('sequence_type', '').lower() == 'dna':
        tumor_id = row['sample_id']
        tumor_path = bam_map.get(tumor_id)
        if not tumor_path:
            # Try fuzzy match to handle sample_id prefixing/suffixing
            matched_key = next((k for k in bam_map.keys() if tumor_id in k or k in tumor_id), None)
            if matched_key:
                tumor_path = bam_map[matched_key]
                tumor_id = matched_key
            else:
                continue
        subject_id = row.get('subject_id')
        group_id = row.get('group_id')
        normal_id = 'NONE'
        normal_path = 'NONE'
        for nrow in samples:
            if nrow.get('sample_type', '').lower() == 'normal' and nrow.get('sequence_type', '').lower() == 'dna':
                if nrow.get('subject_id') == subject_id or (group_id and nrow.get('group_id') == group_id):
                    nid = nrow['sample_id']
                    matched_norm_key = next((k for k in bam_map.keys() if nid in k or k in nid), None)
                    if matched_norm_key:
                        normal_id = matched_norm_key
                        normal_path = bam_map[matched_norm_key]
                        break
        print(f"{tumor_id}\t{tumor_path}\t{normal_id}\t{normal_path}")
EOF

        mapfile -t PAIRED_SAMPLES < <(python3 "${PY_SCRIPT}" "${SAMPLESHEET}" "${COHORT_BAMS_FILE}")
        rm -f "${PY_SCRIPT}" "${COHORT_BAMS_FILE}"

        echo "Running paired/tumor-only off-target ichorCNA for targeted mode..." >&2
        for row in "${PAIRED_SAMPLES[@]}"; do
            tumor_id=$(echo "$row" | cut -f1)
            tumor_path=$(echo "$row" | cut -f2)
            normal_id=$(echo "$row" | cut -f3)
            normal_path=$(echo "$row" | cut -f4)

            echo "  Processing Tumor: ${tumor_id} (Matched Normal: ${normal_id})" >&2
            SAMPLE_ICHOR_DIR="${OUTDIR}/ichorCNA_panel/${tumor_id}"
            mkdir -p "${SAMPLE_ICHOR_DIR}"

            # Determine chromosomes list and genome build
            if [[ "${GENOME}" == "GRCh37_hmf" ]]; then
                CHRS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
                GENOME_BUILD="hg19"
                CHR_NORM="c(1:22)"
                CHR_TRAIN="c(1:22)"
                CHR_ANALYZE="c(1:22,'X')"
                GENOME_STYLE="NCBI"
            else
                CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
                GENOME_BUILD="hg38"
                CHR_NORM="paste0('chr', 1:22)"
                CHR_TRAIN="paste0('chr', 1:22)"
                CHR_ANALYZE="c(paste0('chr', 1:22),'chrX')"
                GENOME_STYLE="UCSC"
            fi

            # Estimate off-target coverage for tumor
            total_reads=$(samtools idxstats "${tumor_path}" | awk '{sum+=$3} END {print sum}')
            on_target_reads=$(samtools view -c -L "${ICHOR_EXCLUSION_BED}" "${tumor_path}")
            off_target_reads=$((total_reads - on_target_reads))
            if [[ ${off_target_reads} -lt 0 ]]; then off_target_reads=0; fi
            off_target_cov=$(awk -v reads="${off_target_reads}" 'BEGIN {print (reads * 150) / 3100000000}')
            echo "    Tumor off-target coverage: ${off_target_cov}x" >&2

            downsample_bam="${tumor_path}"
            is_temp_bam=false
            needs_downsample=$(awk -v cov="${off_target_cov}" 'BEGIN {print (cov > 3.0) ? "true" : "false"}')
            if [[ "${needs_downsample}" == "true" ]]; then
                fraction=$(awk -v cov="${off_target_cov}" 'BEGIN {print 3.0 / cov}')
                s_param=$(awk -v frac="${fraction}" 'BEGIN {printf "%.6f", 100 + frac}')
                echo "    Downsampling tumor off-target reads to 3.0x..." >&2
                temp_bam="${SAMPLE_ICHOR_DIR}/${tumor_id}.downsampled.bam"
                samtools view -b -s "${s_param}" "${tumor_path}" > "${temp_bam}"
                samtools index "${temp_bam}"
                downsample_bam="${temp_bam}"
                is_temp_bam=true
            fi

            # Count tumor reads
            raw_wig="${SAMPLE_ICHOR_DIR}/${tumor_id}.raw.wig"
            readCounter --window 1000000 --quality 20 --chromosome "${CHRS}" "${downsample_bam}" > "${raw_wig}"
            if [[ "${is_temp_bam}" == "true" ]]; then
                rm -f "${downsample_bam}" "${downsample_bam}.bai"
            fi

            # Filter tumor wig
            tumor_wig="${SAMPLE_ICHOR_DIR}/${tumor_id}.wig"
            python3 "${SOURCE_REPO}/scripts/filter_wig.py" "${raw_wig}" "${ICHOR_EXCLUSION_BED}" "${tumor_wig}"
            rm -f "${raw_wig}"

            # Process normal WIG if matched normal exists
            normal_wig=""
            if [[ "${normal_path}" != "NONE" ]]; then
                # Estimate off-target coverage for normal
                n_total_reads=$(samtools idxstats "${normal_path}" | awk '{sum+=$3} END {print sum}')
                n_on_target_reads=$(samtools view -c -L "${ICHOR_EXCLUSION_BED}" "${normal_path}")
                n_off_target_reads=$((n_total_reads - n_on_target_reads))
                if [[ ${n_off_target_reads} -lt 0 ]]; then n_off_target_reads=0; fi
                n_off_target_cov=$(awk -v reads="${n_off_target_reads}" 'BEGIN {print (reads * 150) / 3100000000}')
                echo "    Normal off-target coverage: ${n_off_target_cov}x" >&2

                n_downsample_bam="${normal_path}"
                n_is_temp_bam=false
                n_needs_downsample=$(awk -v cov="${n_off_target_cov}" 'BEGIN {print (cov > 3.0) ? "true" : "false"}')
                if [[ "${n_needs_downsample}" == "true" ]]; then
                    n_fraction=$(awk -v cov="${n_off_target_cov}" 'BEGIN {print 3.0 / cov}')
                    n_s_param=$(awk -v frac="${n_fraction}" 'BEGIN {printf "%.6f", 100 + frac}')
                    echo "    Downsampling normal off-target reads to 3.0x..." >&2
                    n_temp_bam="${SAMPLE_ICHOR_DIR}/${normal_id}.downsampled.bam"
                    samtools view -b -s "${n_s_param}" "${normal_path}" > "${n_temp_bam}"
                    samtools index "${n_temp_bam}"
                    n_downsample_bam="${n_temp_bam}"
                    n_is_temp_bam=true
                fi

                n_raw_wig="${SAMPLE_ICHOR_DIR}/${normal_id}.raw.wig"
                readCounter --window 1000000 --quality 20 --chromosome "${CHRS}" "${n_downsample_bam}" > "${n_raw_wig}"
                if [[ "${n_is_temp_bam}" == "true" ]]; then
                    rm -f "${n_downsample_bam}" "${n_downsample_bam}.bai"
                fi

                normal_wig="${SAMPLE_ICHOR_DIR}/${normal_id}.wig"
                python3 "${SOURCE_REPO}/scripts/filter_wig.py" "${n_raw_wig}" "${ICHOR_EXCLUSION_BED}" "${normal_wig}"
                rm -f "${n_raw_wig}"
            fi

            # Resolve GC, Map, and Centromere file paths dynamically from R package
            if [[ "${GENOME_BUILD}" == "hg38" ]]; then
                CENT_FILE="GRCh38.GCA_000001405.2_centromere_acen.txt"
            else
                CENT_FILE="GRCh37.p13_centromere_UCSC-gapTable.txt"
            fi
            GC_WIG=$(Rscript -e "cat(system.file('extdata', 'gc_${GENOME_BUILD}_1000kb.wig', package='ichorCNA'))")
            MAP_WIG=$(Rscript -e "cat(system.file('extdata', 'map_${GENOME_BUILD}_1000kb.wig', package='ichorCNA'))")
            CENTROMERE=$(Rscript -e "cat(system.file('extdata', '${CENT_FILE}', package='ichorCNA'))")

            SC_STATES_VAL="NULL"
            if [[ "${ICHOR_SUB_CLONES}" == "true" ]]; then
                SC_STATES_VAL="c(1,3)"
            fi

            NORMWIG_ARG=""
            PON_ARG_STANDARD=""
            if [[ -n "${normal_wig}" && -f "${normal_wig}" ]]; then
                NORMWIG_ARG="--NORMWIG ${normal_wig}"
            else
                PON_ARG_STANDARD="--normalPanel ${ICHOR_PON}"
            fi

            # Step 3: Run Rscript runIchorCNA.R inside the output directory
            (
                cd "${SAMPLE_ICHOR_DIR}"
                Rscript "${SOURCE_REPO}/resources/runIchorCNA.R" \
                    --WIG "${tumor_wig}" \
                    ${NORMWIG_ARG} \
                    --id "${tumor_id}" \
                    --genomeStyle "${GENOME_STYLE}" \
                    --genomeBuild "${GENOME_BUILD}" \
                    --gcWig "${GC_WIG}" \
                    --mapWig "${MAP_WIG}" \
                    --centromere "${CENTROMERE}" \
                    ${PON_ARG_STANDARD} \
                    --chrNormalize "${CHR_NORM}" \
                    --chrTrain "${CHR_TRAIN}" \
                    --chrs "${CHR_ANALYZE}" \
                    --normal "c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)" \
                    --ploidy "c(2,3)" \
                    --scStates "${SC_STATES_VAL}"
            )

            # Flatten any nested sample subfolder created by plotSolutions
            if [[ -d "${SAMPLE_ICHOR_DIR}/${tumor_id}" ]]; then
                mv -f "${SAMPLE_ICHOR_DIR}/${tumor_id}"/* "${SAMPLE_ICHOR_DIR}/" 2>/dev/null || true
                rmdir "${SAMPLE_ICHOR_DIR}/${tumor_id}" 2>/dev/null || true
            fi
            
            # Calculate target regions length from exclusion bed
            if [[ -n "${ICHOR_EXCLUSION_BED}" && -f "${ICHOR_EXCLUSION_BED}" ]]; then
                if [[ "${ICHOR_EXCLUSION_BED}" == *.gz ]]; then
                    target_length=$(zcat "${ICHOR_EXCLUSION_BED}" | awk 'BEGIN {sum=0} {sum+=($3-$2)} END {print sum}')
                else
                    target_length=$(awk 'BEGIN {sum=0} {sum+=($3-$2)} END {print sum}' "${ICHOR_EXCLUSION_BED}")
                fi
            else
                target_length=0
            fi
            
            # Calculate coverage
            if [[ ${target_length} -gt 0 ]]; then
                on_target_cov=$(awk -v reads="${on_target_reads}" -v len="${target_length}" 'BEGIN {print (reads * 150) / len}')
                off_target_cov=$(awk -v reads="${off_target_reads}" -v len="${target_length}" 'BEGIN {print (reads * 150) / (3088286401 - len)}')
            else
                on_target_cov=0.0
                off_target_cov=$(awk -v reads="${off_target_reads}" 'BEGIN {print (reads * 150) / 3100000000}')
            fi
            
            # Write QC JSON
            mkdir -p "${SAMPLE_ICHOR_DIR}/${tumor_id}"
            cat <<EOF > "${SAMPLE_ICHOR_DIR}/qc_stats.json"
{
  "total_reads": ${total_reads},
  "on_target_reads": ${on_target_reads},
  "on_target_coverage": ${on_target_cov},
  "off_target_reads": ${off_target_reads},
  "off_target_coverage": ${off_target_cov}
}
EOF
            cat <<EOF > "${SAMPLE_ICHOR_DIR}/${tumor_id}/qc_stats.json"
{
  "total_reads": ${total_reads},
  "on_target_reads": ${on_target_reads},
  "on_target_coverage": ${on_target_cov},
  "off_target_reads": ${off_target_reads},
  "off_target_coverage": ${off_target_cov}
}
EOF
            
            rm -f "${tumor_wig}"
            if [[ -n "${normal_wig}" && -f "${normal_wig}" ]]; then
                rm -f "${normal_wig}"
            fi
        done
        
        # Point PANEL_ICHOR_DIR to this output folder so downstream summarize_findings.py automatically uses it
        PANEL_ICHOR_DIR="${OUTDIR}/ichorCNA_panel"
        echo "Completed off-target ichorCNA execution. Sourced to summarize_findings: ${PANEL_ICHOR_DIR}" >&2
        
        # Write comparison TSV if an original ichor-dir was provided
        if [[ -n "${ORIGINAL_ICHOR_DIR}" && -d "${ORIGINAL_ICHOR_DIR}" ]]; then
            comparison_tsv="${OUTDIR}/ichorCNA_panel/tumor_fraction_comparison.tsv"
            echo -e "sample_id\tpanel_tumor_fraction\tprovided_tumor_fraction\thmf_tumor_fraction" > "${comparison_tsv}"
            for bam in "${COHORT_BAMS[@]}"; do
                sample_id="$(basename "${bam}" .redux.bam)"
                if [[ "${sample_id}" == *-N || "${sample_id}" == *-N-* || "${sample_id}" == *normal* || "${sample_id}" == *BC* ]]; then
                    continue
                fi
                
                # Extract patient ID (e.g., UHN-90202) for matching provided sWGS metrics
                patient_id=$(echo "${sample_id}" | grep -oE "[A-Z]+-[0-9]+" | head -n 1 || echo "${sample_id}")
                
                # Extract panel TF
                panel_tf="NA"
                if [[ -d "${OUTDIR}/ichorCNA_panel/${sample_id}" ]]; then
                    panel_params=$(find "${OUTDIR}/ichorCNA_panel/${sample_id}" -type f -name "*params.txt" 2>/dev/null | head -n 1)
                    if [[ -z "${panel_params}" ]]; then
                        panel_params=$(find "${OUTDIR}/ichorCNA_panel/${sample_id}" -type f -name "*metrics.txt" 2>/dev/null | head -n 1)
                    fi
                    if [[ -n "${panel_params}" && -f "${panel_params}" ]]; then
                        panel_tf=$(grep -i "tumor fraction" "${panel_params}" | awk -F: '{print $2}' | xargs 2>/dev/null || echo "NA")
                    fi
                fi
                
                # Extract provided TF
                provided_tf="NA"
                if [[ -d "${ORIGINAL_ICHOR_DIR}" ]]; then
                    # 1. Try finding params matching the patient ID
                    provided_params=$(find "${ORIGINAL_ICHOR_DIR}" -type f -name "*${patient_id}*params.txt" 2>/dev/null | head -n 1)
                    if [[ -z "${provided_params}" ]]; then
                        provided_params=$(find "${ORIGINAL_ICHOR_DIR}" -type f -name "*${patient_id}*metrics.txt" 2>/dev/null | head -n 1)
                    fi
                    # 2. Fallback: find any params/metrics file in the folder
                    if [[ -z "${provided_params}" ]]; then
                        provided_params=$(find "${ORIGINAL_ICHOR_DIR}" -type f -name "*params.txt" 2>/dev/null | head -n 1)
                    fi
                    if [[ -z "${provided_params}" ]]; then
                        provided_params=$(find "${ORIGINAL_ICHOR_DIR}" -type f -name "*metrics.txt" 2>/dev/null | head -n 1)
                    fi
                    
                    if [[ -n "${provided_params}" && -f "${provided_params}" ]]; then
                        provided_tf=$(grep -i "tumor fraction" "${provided_params}" | awk -F: '{print $2}' | xargs 2>/dev/null || echo "NA")
                    fi
                fi
                
                # Extract HMF pipeline tumor fraction (PURPLE purity)
                hmf_tf="NA"
                purple_purity_file=$(find "${OUTDIR}/${sample_id}/purple" -type f -name "*.purple.purity.tsv" 2>/dev/null | head -n 1)
                if [[ -n "${purple_purity_file}" && -f "${purple_purity_file}" ]]; then
                    hmf_tf=$(awk 'NR==2 {print $1}' "${purple_purity_file}" 2>/dev/null || echo "NA")
                fi
                
                echo -e "${sample_id}\t${panel_tf}\t${provided_tf}\t${hmf_tf}" >> "${comparison_tsv}"
            done
            echo "Saved tumor fraction comparison TSV to: ${comparison_tsv}" >&2
        fi
    fi
fi

# ---------------------------------------------------------------------------
# Post-success conversion: BAM -> CRAM (then remove BAMs)
# ---------------------------------------------------------------------------
if [[ "${CONVERT_FINAL_BAM_TO_CRAM}" == "true" ]]; then
    # Read FASTA path directly from h4h.config for consistent CRAM reference.
    if [[ "${GENOME}" == "GRCh37_hmf" ]]; then
        REF_FASTA="$(awk -F"'" '/ref_data_genome_fasta\s*=/{print $2; exit}' "${CUSTOM_CONFIG}")"
    else
        REF_FASTA="$(awk -F"'" '/ref_data_genome_fasta\s*=/{count++; if(count==2){print $2; exit}}' "${CUSTOM_CONFIG}")"
    fi
    if [[ -z "${REF_FASTA}" || ! -f "${REF_FASTA}" ]]; then
        echo "ERROR: Could not resolve ref_data_genome_fasta from ${CUSTOM_CONFIG}" >&2
        exit 1
    fi

    module load samtools/1.20 >/dev/null 2>&1 || module load samtools >/dev/null 2>&1 || true
    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: samtools is required for BAM->CRAM conversion but was not found in PATH" >&2
        exit 1
    fi

    mapfile -t BAM_FILES < <(find "${OUTDIR}" -type f -name "*.redux.bam" \
        ! -path "${OUTDIR}/work/*" \
        ! -path "${OUTDIR}/pipeline_info/*" \
        ! -path "${OUTDIR}/logs/*")

    if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
        echo "No BAM files found for post-run conversion under ${OUTDIR}" >&2
    else
        echo "Converting ${#BAM_FILES[@]} BAM file(s) to CRAM..." >&2
        converted=0
        local_threads="${SLURM_CPUS_PER_TASK:-4}"
        for bam in "${BAM_FILES[@]}"; do
            cram="${bam%.bam}.cram"
            echo "  BAM -> CRAM: ${bam}" >&2
            if samtools view -@ "${local_threads}" -C -T "${REF_FASTA}" -o "${cram}" "${bam}" 2>/dev/null; then
                samtools index -@ "${local_threads}" "${cram}" 2>/dev/null || true
                rm -f "${bam}" "${bam}.bai" "${bam%.bam}.bai" 2>/dev/null || true
                ((converted++)) || true
            elif samtools view -@ "${local_threads}" -C --output-fmt cram,embed_ref=2 -T "${REF_FASTA}" -o "${cram}" "${bam}" 2>/dev/null; then
                samtools index -@ "${local_threads}" "${cram}" 2>/dev/null || true
                rm -f "${bam}" "${bam}.bai" "${bam%.bam}.bai" 2>/dev/null || true
                ((converted++)) || true
            else
                echo "  WARNING: Could not convert ${bam} to CRAM; leaving original BAM intact." >&2
            fi
        done
        echo "Completed BAM->CRAM conversion for ${converted} file(s)." >&2
    fi
fi

# ---------------------------------------------------------------------------
# Synchronize parallel support workflows (sWGS ichorCNA & EM-Seq methylation)
# ---------------------------------------------------------------------------
if [[ ${#PARALLEL_SUPPORT_JOB_IDS[@]} -gt 0 ]]; then
    echo "Waiting for parallel support SLURM jobs to complete: ${PARALLEL_SUPPORT_JOB_IDS[*]}..." >&2
    JOB_IDS_COMMA=$(IFS=,; echo "${PARALLEL_SUPPORT_JOB_IDS[*]}")
    while true; do
        active_jobs=$(squeue -j "${JOB_IDS_COMMA}" -h -o "%A" 2>/dev/null || true)
        if [[ -z "${active_jobs}" ]]; then
            break
        fi
        sleep 20
    done

    # Check exit statuses
    for jid in "${PARALLEL_SUPPORT_JOB_IDS[@]}"; do
        state=$(sacct -j "$jid" -n -o State 2>/dev/null | head -n 1 | xargs || true)
        if [[ -n "$state" && "$state" != "COMPLETED" ]]; then
            echo "WARNING: Parallel support job ${jid} finished with non-COMPLETED state: ${state}" >&2
        else
            echo "Parallel support job ${jid} completed successfully." >&2
        fi
    done
fi

if [[ ${#PARALLEL_SUPPORT_PIDS[@]} -gt 0 ]]; then
    echo "Waiting for parallel support background processes to complete: ${PARALLEL_SUPPORT_PIDS[*]}..." >&2
    for pid in "${PARALLEL_SUPPORT_PIDS[@]}"; do
        wait "$pid" 2>/dev/null || true
    done
fi

# Resolve completed sWGS/WGS ichorCNA results
if [[ -z "${ICHOR_DIR}" || ! -d "${ICHOR_DIR}" ]]; then
    for cand_ichor in \
        "${PATIENT_BASEDIR}/ichorCNA_tumor_only" \
        "${PATIENT_BASEDIR}/ichorCNA_plasma_only" \
        "${PATIENT_BASEDIR}/ichorCNA_t_vs_n" \
        "${PATIENT_BASEDIR}/ichorCNA_plasma_vs_n" \
        "${PATIENT_BASEDIR}/ichorCNA" \
        "${OUTDIR}/ichorCNA_tumor_only" \
        "${OUTDIR}/ichorCNA_plasma_only" \
        "${OUTDIR}/ichorCNA_t_vs_n" \
        "${OUTDIR}/ichorCNA_plasma_vs_n" \
        "${OUTDIR}/ichorCNA"; do
        if [[ -d "${cand_ichor}" ]] && find "${cand_ichor}" -type f -name "*params.txt" 2>/dev/null | grep -q .; then
            ICHOR_DIR="${cand_ichor}"
            echo "Resolved ichorCNA results: ${ICHOR_DIR}" >&2
            break
        fi
    done
fi

# Resolve completed EM-Seq methylation results
if [[ -z "${METHYLATION_DATA}" || ! -f "${METHYLATION_DATA}" ]]; then
    EM_TARGET_OUT="${EMSEQ_DIR:-${OUTDIR}/EMSeq}"
    NEW_METH=$(find "${EM_TARGET_OUT}" -type f -name "*CpG_methylation_matrix.RData" 2>/dev/null | head -n 1 || true)
    if [[ -z "${NEW_METH}" && -d "${OUTDIR}/MethylDackel" ]]; then
        NEW_METH=$(find "${OUTDIR}/MethylDackel" -type f -name "*CpG_methylation_matrix.RData" 2>/dev/null | head -n 1 || true)
    fi
    if [[ -n "${NEW_METH}" && -f "${NEW_METH}" ]]; then
        METHYLATION_DATA="${NEW_METH}"
        echo "Resolved EM-Seq methylation matrix: ${METHYLATION_DATA}" >&2
    fi
fi

# ---------------------------------------------------------------------------
# Generate final summary report
# ---------------------------------------------------------------------------
if [[ "${IS_RNA_ONLY}" == "true" ]]; then
    echo "Running RNA-only post-processing (fusion neoantigen HLA binding and TCR repertoire analysis)..."
    RNA_POSTPROC_SCRIPT="${SOURCE_REPO}/scripts/run_rna_postprocessing.py"
    if [[ -f "${RNA_POSTPROC_SCRIPT}" ]]; then
        PYTHON_EXEC="python3"
        if [[ -x "${HOME}/.conda/envs/nf-core/bin/python3" ]]; then
            PYTHON_EXEC="${HOME}/.conda/envs/nf-core/bin/python3"
        elif command -v conda &>/dev/null && conda info --envs 2>/dev/null | grep -q "nf-core"; then
            PYTHON_EXEC="conda run -n nf-core python3"
        fi
        declare -a RNA_POST_ARGS=()
        RNA_POST_ARGS+=(--output-dir "${OUTDIR}")
        RNA_POST_ARGS+=(--samplesheet "${SAMPLESHEET}")
        if [[ -n "${METADATA_TSV}" ]]; then
            RNA_POST_ARGS+=(--metadata-tsv "${METADATA_TSV}")
        fi
        if [[ -n "${MIXCR_PATH}" ]]; then
            RNA_POST_ARGS+=(--mixcr-path "${MIXCR_PATH}")
        fi
        ${PYTHON_EXEC} "${RNA_POSTPROC_SCRIPT}" "${RNA_POST_ARGS[@]}" || echo "WARN: RNA post-processing script exited with non-zero status, continuing." >&2
    else
        echo "WARN: RNA post-processing script not found at ${RNA_POSTPROC_SCRIPT}" >&2
    fi
else
    echo "Generating final summary report for panel/DNA run..."
    SUMMARIZE_SCRIPT="${SOURCE_REPO}/scripts/summarize_findings.py"
    if [[ -f "${SUMMARIZE_SCRIPT}" ]]; then
        declare -a SUM_ARGS=()
        SUM_ARGS+=(--output-dir "${OUTDIR}")
        SUM_ARGS+=(--samplesheet "${SAMPLESHEET}")
    
    if [[ -n "${METHYLATION_DATA}" ]]; then
        SUM_ARGS+=(--methylation-data "${METHYLATION_DATA}")
    fi
    
    if [[ -n "${IHC_TSV}" ]]; then
        SUM_ARGS+=(--ihc-tsv "${IHC_TSV}")
    fi
    
    if [[ -n "${IHC_FILE}" ]]; then
        SUM_ARGS+=(--ihc-file "${IHC_FILE}")
    fi
    
    if [[ -n "${ORTHOGONAL_TSV}" ]]; then
        SUM_ARGS+=(--orthogonal-tsv "${ORTHOGONAL_TSV}")
    fi
    
    if [[ -n "${AGE}" ]]; then
        SUM_ARGS+=(--age "${AGE}")
    fi
    
    if [[ -n "${TISSUE}" ]]; then
        SUM_ARGS+=(--tissue "${TISSUE}")
    fi
    
    if [[ -n "${ICHOR_DIR}" ]]; then
        SUM_ARGS+=(--ichor-dir "${ICHOR_DIR}")
    fi

    if [[ -n "${PURITY_SOURCE}" ]]; then
        SUM_ARGS+=(--purity-source "${PURITY_SOURCE}")
    fi

    if [[ -n "${PANEL_ICHOR_DIR:-}" && -d "${PANEL_ICHOR_DIR:-}" ]]; then
        SUM_ARGS+=(--panel-ichor-dir "${PANEL_ICHOR_DIR}")
    elif [[ -d "${OUTDIR}/ichorCNA_panel" ]]; then
        SUM_ARGS+=(--panel-ichor-dir "${OUTDIR}/ichorCNA_panel")
    fi
    
    if [[ -n "${PANEL_NAME}" ]]; then
        SUM_ARGS+=(--project "${PANEL_NAME}")
    else
        SUM_ARGS+=(--project "oncoanalyser")
    fi

    if [[ -n "${GENOME}" ]]; then
        SUM_ARGS+=(--genome "${GENOME}")
    fi

    # Extract FASTA from config for IGV report
    REF_FASTA=""
    if [[ -f "${CUSTOM_CONFIG}" ]]; then
        if [[ "${GENOME}" == "GRCh37_hmf" ]]; then
            REF_FASTA="$(awk -F"'" '/ref_data_genome_fasta\s*=/{print $2; exit}' "${CUSTOM_CONFIG}")"
        else
            REF_FASTA="$(awk -F"'" '/ref_data_genome_fasta\s*=/{count++; if(count==2){print $2; exit}}' "${CUSTOM_CONFIG}")"
        fi
    fi
    if [[ -n "${REF_FASTA}" && -f "${REF_FASTA}" ]]; then
        SUM_ARGS+=(--fasta "${REF_FASTA}")
    fi

    # Locate MANE GTF for IGV track
    MANE_GTF="/cluster/projects/pughlab/references/RefSeq/SYMBOL_MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"
    if [[ -f "${MANE_GTF}" ]]; then
        SUM_ARGS+=(--gtf "${MANE_GTF}")
    fi
    
    # Try to load R module for Rscript (needed for RData methylation extraction)
    if command -v module &>/dev/null; then
        module load R/4.5.2 2>/dev/null || module load r/4.3.1 2>/dev/null || module load R 2>/dev/null || module load r 2>/dev/null || true
    fi

    # Try to activate nf-core conda environment for pandas/matplotlib (needed for summary and plotting)
    if command -v conda &>/dev/null; then
        CONDA_BASE=$(conda info --base 2>/dev/null)
        if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
            source "${CONDA_BASE}/etc/profile.d/conda.sh"
            conda activate nf-core 2>/dev/null || true
        fi
    fi

    # Determine Python execution command (prioritize absolute environment python path to bypass compute node shell issues)
    PYTHON_EXEC="python3"
    if [[ -x "${HOME}/.conda/envs/nf-core/bin/python3" ]]; then
        PYTHON_EXEC="${HOME}/.conda/envs/nf-core/bin/python3"
    elif command -v conda &>/dev/null && conda info --envs 2>/dev/null | grep -q "nf-core"; then
        PYTHON_EXEC="conda run -n nf-core python3"
    fi

        # Run the summarizer
        ${PYTHON_EXEC} "${SUMMARIZE_SCRIPT}" "${SUM_ARGS[@]}" || echo "WARN: Summarizer script failed, continuing." >&2
    else
        echo "WARN: Summarizer script not found at ${SUMMARIZE_SCRIPT}" >&2
    fi
fi

echo "Pipeline complete: $(date)"
