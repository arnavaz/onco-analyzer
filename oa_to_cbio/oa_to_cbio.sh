#!/bin/bash
# Master orchestrator script for Oncoanalyser to cBioPortal conversion
# Usage: ./oa_to_cbio.sh -i <oncoanalyser_out_dir> -o <cbio_study_dir> -A <slurm_account> [-s study_id] [-n study_name] [-d study_desc]

set -e
umask 0002

# Default values
STUDY_ID="oa_study"
STUDY_NAME="Oncoanalyser Study"
STUDY_DESC="Data converted from nf-core/oncoanalyser."
ACCOUNT=""
TYPE_OF_CANCER="mixed"
SCRIPT_DIR_ORCH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="${SCRIPT_DIR_ORCH}/cbio.conf"

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -A <slurm_account> [-s study_id] [-n study_name] [-d study_desc] [-t type_of_cancer] [-c conf_file] [-S pre_sample_clinical] [-P pre_patient_clinical] [-C pre_clin_dir] [-u | --unfiltered-maf] [--non-germline-filter] [--unfiltered-germline-maf] [-g | --driver-gene-panel <path>] [-m | --merge-only]"
    exit 1
}

UNFILTERED_MAF="false"
NON_GERMLINE_FILTER="false"
UNFILTERED_GERMLINE_MAF="false"
DRIVER_GENE_PANEL=""
MERGE_ONLY="false"

while [ $# -gt 0 ]; do
    case "$1" in
        -i) IN_DIR="$2"; shift 2 ;;
        -o) OUT_DIR="$2"; shift 2 ;;
        -A) ACCOUNT="$2"; shift 2 ;;
        -s) STUDY_ID="$2"; shift 2 ;;
        -n) STUDY_NAME="$2"; shift 2 ;;
        -d) STUDY_DESC="$2"; shift 2 ;;
        -t) TYPE_OF_CANCER="$2"; shift 2 ;;
        -c) CONF_FILE="$2"; shift 2 ;;
        -S) PRE_SAMPLE="$2"; shift 2 ;;
        -P) PRE_PATIENT="$2"; shift 2 ;;
        -C) PRE_CLIN_DIR="$2"; shift 2 ;;
        -u|--unfiltered-maf) UNFILTERED_MAF="true"; shift 1 ;;
        --non-germline-filter) NON_GERMLINE_FILTER="true"; shift 1 ;;
        --unfiltered-germline-maf) UNFILTERED_GERMLINE_MAF="true"; shift 1 ;;
        -g|--driver-gene-panel) DRIVER_GENE_PANEL="$2"; shift 2 ;;
        -m|--merge-only|--merge_only|--clinical-only) MERGE_ONLY="true"; shift 1 ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
done

# Auto-enable merge-only (clinical-only) if IN_DIR does not exist or is empty
if [ "$MERGE_ONLY" = "false" ]; then
    if [ -z "$IN_DIR" ] || [ ! -d "$IN_DIR" ] || [ -z "$(find -L "$IN_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)" ]; then
        if [ -n "$PRE_SAMPLE" ] || [ -n "$PRE_PATIENT" ] || [ -n "$PRE_CLIN_DIR" ]; then
            echo "Input directory (${IN_DIR:-<none>}) does not exist or is empty -> enabling clinical-only (--merge-only) mode automatically."
            MERGE_ONLY="true"
            IN_DIR=""
        else
            echo "Error: Missing or empty input directory (-i) and no pre-clinical data provided."
            usage
        fi
    fi
fi

if [ -z "$OUT_DIR" ] || [ -z "$ACCOUNT" ]; then
    echo "Error: Missing required arguments."
    usage
fi

if [ -n "$IN_DIR" ] && [ -d "$IN_DIR" ]; then
    IN_DIR=$(realpath "$IN_DIR")
fi
mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

if [ -f "$CONF_FILE" ]; then
    echo "Loading configuration from $CONF_FILE"
    source "$CONF_FILE"
    export AMP_RATIO GAIN_RATIO HET_DEL_RATIO HOM_DEL_ABS GTF_FILE
else
    echo "Warning: Config file $CONF_FILE not found, using default CNV thresholds."
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
export SCRIPT_DIR

echo "Initializing output directory: $OUT_DIR"
mkdir -p "${OUT_DIR}/logs"

IMPORT_DIR="${OUT_DIR}/import_folder"
mkdir -p "$IMPORT_DIR"

# Pre-copy meta_study.txt from PRE_CLIN_DIR if present to prevent generate_meta.sh from overwriting it
if [ -n "$PRE_CLIN_DIR" ] && [ -f "${PRE_CLIN_DIR}/meta_study.txt" ] && [ ! -f "${IMPORT_DIR}/meta_study.txt" ]; then
    echo "Found pre-generated meta_study.txt in ${PRE_CLIN_DIR}. Pre-copying to preserve manual modifications."
    cp "${PRE_CLIN_DIR}/meta_study.txt" "${IMPORT_DIR}/meta_study.txt"
fi

# Generate meta files inside import_folder
"${SCRIPT_DIR}/generate_meta.sh" "$IMPORT_DIR" "$STUDY_ID" "$STUDY_NAME" "$STUDY_DESC" "$TYPE_OF_CANCER"

# We no longer write dummy headers here since merge_clinical.py handles full creation or merging

JOB_IDS=()

if [ "$MERGE_ONLY" = "false" ]; then
    echo "Finding patients and samples in $IN_DIR..."
    # Loop through patient directories
    for PATIENT_DIR in "$IN_DIR"/*/; do
        if [ ! -d "$PATIENT_DIR" ]; then continue; fi
        PATIENT_ID=$(basename "$PATIENT_DIR")
        
        # Skip non-patient directories
        case "$PATIENT_ID" in
            fingerprint_metrics|gene_expression_files|logs|orange_links|runtime_summary|pipeline_info|work)
                continue
                ;;
        esac
        
        # Handle both Nested (Oncoanalyser) and Flat (HMFTools) structures
        SAMPLE_DIRS=()
        if [ -d "${PATIENT_DIR}/purple" ] || [ -d "${PATIENT_DIR}/orange" ]; then
            # Flat structure: PATIENT_DIR is actually the SAMPLE_DIR
            SAMPLE_DIRS+=("$PATIENT_DIR")
        else
            # Nested structure: Look for sample subdirectories
            for D2 in "$PATIENT_DIR"/*/; do
                if [ ! -d "$D2" ]; then continue; fi
                D2_NAME=$(basename "$D2")
                if [ "$D2_NAME" = "pipeline_info" ] || [ "$D2_NAME" = "work" ]; then continue; fi
                if [ -d "${D2}/purple" ] || [ -d "${D2}/orange" ]; then
                    SAMPLE_DIRS+=("$D2")
                fi
            done
        fi
        
        for SAMPLE_DIR in "${SAMPLE_DIRS[@]}"; do
            TUMOR_ID=$(basename "$SAMPLE_DIR")
            NORMAL_ID="UNKNOWN_NORMAL"
            
            # Dynamically determine TUMOR_ID and NORMAL_ID from the somatic VCF
            SOMATIC_VCF=$(ls "${SAMPLE_DIR}/purple/"*.purple.somatic.vcf.gz 2>/dev/null | head -n 1 || true)
            if [ -n "$SOMATIC_VCF" ]; then
                VCF_BASE=$(basename "$SOMATIC_VCF")
                TUMOR_ID="${VCF_BASE%.purple.somatic.vcf.gz}"
                
                # Extract sample columns from VCF header (Cols 10+)
                SAMPLES=$(zgrep -m 1 "^#CHROM" "$SOMATIC_VCF" 2>/dev/null | cut -f 10- || true)
                if [ -n "$SAMPLES" ]; then
                    S1=$(echo "$SAMPLES" | awk '{print $1}')
                    S2=$(echo "$SAMPLES" | awk '{print $2}')
                    
                    if [ "$S1" = "$TUMOR_ID" ] && [ -n "$S2" ]; then
                        NORMAL_ID="$S2"
                    elif [ "$S2" = "$TUMOR_ID" ] && [ -n "$S1" ]; then
                        NORMAL_ID="$S1"
                    elif [ -n "$S1" ] && [ -n "$S2" ]; then
                        # Fallback if neither matches exactly (rare)
                        NORMAL_ID="$S2"
                    fi
                fi
            fi
            
            if [ -z "$NORMAL_ID" ] || [ "$NORMAL_ID" = "$TUMOR_ID" ]; then
                NORMAL_ID="UNKNOWN_NORMAL"
            fi
            
            echo "Found sample: Tumor=$TUMOR_ID, Normal=$NORMAL_ID"
            
            # Submit process_sample job
            export TUMOR_ID NORMAL_ID SAMPLE_DIR OUT_DIR UNFILTERED_MAF NON_GERMLINE_FILTER UNFILTERED_GERMLINE_MAF DRIVER_GENE_PANEL GTF_FILE
            export AMP_RATIO GAIN_RATIO HET_DEL_RATIO HOM_DEL_ABS
            JOB_ID=$(sbatch --parsable -A "$ACCOUNT" --job-name="oa2cbio_${TUMOR_ID}" -o "${OUT_DIR}/logs/${TUMOR_ID}.%j.out" -e "${OUT_DIR}/logs/${TUMOR_ID}.%j.err" --export=ALL -t 06:00:00 -c 2 --mem=28G "${SCRIPT_DIR}/process_sample.sh")
            
            if [ -n "$JOB_ID" ]; then
                echo "  Submitted job $JOB_ID for $TUMOR_ID"
                JOB_IDS+=("$JOB_ID")
            else
                echo "  Failed to submit job for $TUMOR_ID"
            fi
        done
    done
else
    echo "Merge-only mode enabled. Skipping sample processing jobs."
fi

# Sort and uniq steps removed since python merger handles deduplication

# Submit merge job with dependencies
export OUT_DIR IMPORT_DIR IN_DIR PRE_SAMPLE PRE_PATIENT PRE_CLIN_DIR STUDY_ID SCRIPT_DIR
if [ ${#JOB_IDS[@]} -gt 0 ]; then
    DEPENDENCY_STR=$(IFS=:; echo "${JOB_IDS[*]}")
    echo "Submitting merge job dependent on: $DEPENDENCY_STR"
    MERGE_JOB_ID=$(sbatch --parsable -A "$ACCOUNT" --job-name="oa2cbio_merge_${STUDY_ID}" -o "${OUT_DIR}/logs/merge_${STUDY_ID}.%j.out" -e "${OUT_DIR}/logs/merge_${STUDY_ID}.%j.err" --dependency=afterok:$DEPENDENCY_STR --export=OUT_DIR,IMPORT_DIR,IN_DIR,PRE_SAMPLE,PRE_PATIENT,PRE_CLIN_DIR,STUDY_ID,SCRIPT_DIR -t 01:00:00 -c 1 --mem=4G "${SCRIPT_DIR}/merge_study.sh")
    echo "Merge job ID: $MERGE_JOB_ID"
else
    echo "Submitting standalone merge job..."
    MERGE_JOB_ID=$(sbatch --parsable -A "$ACCOUNT" --job-name="oa2cbio_merge_${STUDY_ID}" -o "${OUT_DIR}/logs/merge_${STUDY_ID}.%j.out" -e "${OUT_DIR}/logs/merge_${STUDY_ID}.%j.err" --export=OUT_DIR,IMPORT_DIR,IN_DIR,PRE_SAMPLE,PRE_PATIENT,PRE_CLIN_DIR,STUDY_ID,SCRIPT_DIR -t 01:00:00 -c 1 --mem=4G "${SCRIPT_DIR}/merge_study.sh")
    echo "Merge job ID: $MERGE_JOB_ID"
fi

echo "Pipeline launched."
