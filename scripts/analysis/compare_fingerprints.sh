#!/bin/bash

################################################################################
# compare_fingerprints.sh
#
# PURPOSE:
#   Compare germline variants across a cohort (all-vs-all) to identify related
#   or identical individuals using GATK4 GenotypeConcordance.
#
# INPUT:
#   - Cohort results directory with per-patient subdirectories
#   - Each patient expected to have: [PATIENT]/[TUMOURSAMPLE]/sage/germline/
#   - SAGE germline VCF files: <TUMOUR_SAMPLE_ID>.sage.germline.vcf.gz
#   - VCF files must include both tumor and normal sample genotypes
#
# OUTPUT:
#   - fingerprint_metrics/ directory containing:
#     * sample_inventory.tsv — discovered samples and source VCFs
#     * pairwise_concordance.tsv — all sample-vs-sample comparisons
#     * pairs/<id>.genotype_concordance_*_metrics — per-pair GATK outputs
#
# DEPENDENCIES:
#   - module load gatk/4.6.0.0
#
# USAGE:
#   bash compare_fingerprints.sh <cohort_results_dir> [options]
#
# OPTIONS:
#   --output-dir <path>    Output directory (default: fingerprint_metrics/)
#   --dry-run              Print commands without executing
#   --help                 Show this message
#
################################################################################

set -euo pipefail


# Defaults
COHORT_DIR=""
OUTPUT_DIR="fingerprint_metrics"
DRY_RUN=false

# Source shared utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

`
DEVELOPMENT:
SCRIPT_DIR="/cluster/projects/pughlab/src/oncoanalyser_pughlab/scripts/"
COHORT_DIR="/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser"
OUTPUT_DIR="/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/fingerprint_metrics"
DRY_RUN=false
`

source "${SCRIPT_DIR}/shared/utils.sh"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help)
            grep '^#' "$0" | sed 's/^# //'
            exit 0
            ;;
        -*)
            die "Unknown option: $1"
            ;;
        *)
            COHORT_DIR="$1"
            shift
            ;;
    esac
done

# Validation
[[ -z "$COHORT_DIR" ]] && die "Usage: $0 <cohort_results_dir>"
require_dir "$COHORT_DIR"

# Load modules
log_info "Loading GATK 4.6.0.0..."
module load gatk/4.6.0.0 || die "Failed to load gatk module"

# Verify GATK is available
require_command "gatk"
require_command "find"
require_command "zgrep"

# Setup output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/pairs"
log_info "Output directory: $OUTPUT_DIR"

# Find all SAGE germline VCFs
log_info "Scanning for SAGE germline VCFs in $COHORT_DIR..."
declare -a VCF_FILES

while IFS= read -r vcf; do
    if [[ -f "${vcf}.tbi" ]]; then
        VCF_FILES+=("$vcf")
    fi
done < <(find "$COHORT_DIR" -path "*/sage/germline/*.sage.germline.vcf.gz" -type f)

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    die "No SAGE germline VCFs found in $COHORT_DIR"
fi

log_info "Found ${#VCF_FILES[@]} germline VCF(s)"

# Build sample inventory (tumor + normal samples from each VCF)
declare -a SAMPLE_IDS
declare -a SAMPLE_VCFS

{
    echo -e "sample_id\tvcf_path"

    for vcf in "${VCF_FILES[@]}"; do
        # VCF header sample columns begin at field 10
        header_line="$(zgrep -m1 '^#CHROM' "$vcf")"
        sample_list="$(echo "$header_line" | cut -f10-)"
        IFS=$'\t' read -r -a samples <<< "$sample_list"

        if [[ ${#samples[@]} -eq 0 ]]; then
            log_warn "No samples found in VCF header: $vcf"
            continue
        fi

        for sample in "${samples[@]}"; do
            SAMPLE_IDS+=("$sample")
            SAMPLE_VCFS+=("$vcf")
            echo -e "${sample}\t${vcf}"
        done
    done
} > "$OUTPUT_DIR/sample_inventory.tsv"

if [[ ${#SAMPLE_IDS[@]} -lt 2 ]]; then
    die "Need at least 2 total samples across VCFs for all-vs-all comparison"
fi

log_info "Discovered ${#SAMPLE_IDS[@]} sample entries for comparison"

log_info "Computing all-vs-all genotype concordance..."
echo -e "call_sample\tcall_vcf\ttruth_sample\ttruth_vcf\tmetrics_prefix\tstatus" > "$OUTPUT_DIR/pairwise_concordance.tsv"

pair_count=0
for ((i = 0; i < ${#SAMPLE_IDS[@]}; i++)); do
    for ((j = i + 1; j < ${#SAMPLE_IDS[@]}; j++)); do
        call_sample="${SAMPLE_IDS[$i]}"
        call_vcf="${SAMPLE_VCFS[$i]}"
        truth_sample="${SAMPLE_IDS[$j]}"
        truth_vcf="${SAMPLE_VCFS[$j]}"

        safe_call="${call_sample//[^A-Za-z0-9._-]/_}"
        safe_truth="${truth_sample//[^A-Za-z0-9._-]/_}"
        metrics_prefix="$OUTPUT_DIR/pairs/${safe_call}__vs__${safe_truth}"

        cmd="gatk GenotypeConcordance \
            --CALL_VCF '$call_vcf' \
            --CALL_SAMPLE '$call_sample' \
            --TRUTH_VCF '$truth_vcf' \
            --TRUTH_SAMPLE '$truth_sample' \
            --OUTPUT '$metrics_prefix' \
            --IGNORE_FILTER_STATUS true"

        if $DRY_RUN; then
            log_info "[DRY-RUN] $cmd"
            echo -e "${call_sample}\t${call_vcf}\t${truth_sample}\t${truth_vcf}\t${metrics_prefix}\tDRY_RUN" >> "$OUTPUT_DIR/pairwise_concordance.tsv"
        else
            if eval "$cmd"; then
                echo -e "${call_sample}\t${call_vcf}\t${truth_sample}\t${truth_vcf}\t${metrics_prefix}\tOK" >> "$OUTPUT_DIR/pairwise_concordance.tsv"
            else
                log_warn "GenotypeConcordance failed: ${call_sample} vs ${truth_sample}"
                echo -e "${call_sample}\t${call_vcf}\t${truth_sample}\t${truth_vcf}\t${metrics_prefix}\tFAILED" >> "$OUTPUT_DIR/pairwise_concordance.tsv"
            fi
        fi

        ((pair_count++)) || true
    done
done

log_info "Complete! Results written to $OUTPUT_DIR/"
log_info "  - sample_inventory.tsv: discovered samples and VCF paths"
log_info "  - pairwise_concordance.tsv: all-vs-all run status by pair"
log_info "  - pairs/: per-pair GATK metric files"
log_info "Total comparisons attempted: $pair_count"
