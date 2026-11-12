#!/bin/bash
# SLURM script to merge all sample files
# Usage: sbatch -A <account> --export=... merge_study.sh

# Environment variables expected:
# OUT_DIR

umask 0002
echo "Merging study files in ${OUT_DIR}"

# 1. Merge MAFs
MAF_OUT="${IMPORT_DIR}/data_mutations.txt"
first=true
for maf in "${OUT_DIR}"/maf_tmp/*.maf; do
    if [ -f "$maf" ]; then
        case "$(basename "$maf")" in
            *.somatic.maf|*.germline.maf|*.unfiltered_sage.maf)
                continue
                ;;
        esac
        if [ "$first" = true ]; then
            grep -v "^#" "$maf" | head -n 1 > "$MAF_OUT"
            first=false
        fi
        grep -v "^#" "$maf" | tail -n +2 >> "$MAF_OUT"
    fi
done

# 2. Merge Segments
SEG_OUT="${IMPORT_DIR}/data_segment.txt"
first=true
for seg in "${OUT_DIR}"/seg_tmp/*.seg.tsv; do
    if [ -f "$seg" ]; then
        if [ "$first" = true ]; then
            head -n 1 "$seg" > "$SEG_OUT"
            first=false
        fi
        tail -n +2 "$seg" >> "$SEG_OUT"
    fi
done

# Output copy with .seg extension for easy loading in IGV
if [ -f "$SEG_OUT" ]; then
    cp "$SEG_OUT" "${IMPORT_DIR}/data_segment.seg"
fi


# 3. Merge Discrete CNA
CNA_OUT="${IMPORT_DIR}/data_CNA.txt"
if [ -d "${OUT_DIR}/cna_tmp" ]; then
    python3 "${SCRIPT_DIR}/merge_cna.py" "${OUT_DIR}/cna_tmp" "$CNA_OUT"
fi

# 4. Merge SVs
SV_OUT="${IMPORT_DIR}/data_structural_variants.txt"
first=true
for sv in "${OUT_DIR}"/sv_tmp/*.sv.tsv; do
    if [ -f "$sv" ]; then
        if [ "$first" = true ]; then
            head -n 1 "$sv" > "$SV_OUT"
            first=false
        fi
        tail -n +2 "$sv" >> "$SV_OUT"
    fi
done

# 5. Merge Clinical Data
CLINICAL_SAMPLE_OUT="${IMPORT_DIR}/data_clinical_sample.txt"
CLINICAL_PATIENT_OUT="${IMPORT_DIR}/data_clinical_patient.txt"

# Build arguments for merge_clinical.py
MERGE_CMD=("python3" "${SCRIPT_DIR}/merge_clinical.py" 
           "--tmp-dir" "${OUT_DIR}/clinical_tmp" 
           "--out-sample" "$CLINICAL_SAMPLE_OUT" 
           "--out-patient" "$CLINICAL_PATIENT_OUT")

if [ -n "$PRE_SAMPLE" ] && [ -f "$PRE_SAMPLE" ]; then
    MERGE_CMD+=("--pre-sample" "$PRE_SAMPLE")
fi

if [ -n "$PRE_PATIENT" ] && [ -f "$PRE_PATIENT" ]; then
    MERGE_CMD+=("--pre-patient" "$PRE_PATIENT")
fi

# Run python merger to handle dynamic columns and key-value formats
"${MERGE_CMD[@]}"

# 6. Merge Expression Data
EXPR_OUT="${IMPORT_DIR}/data_mrna_seq_v2_rsem.txt"
if [ -d "${OUT_DIR}/expr_tmp" ]; then
    python3 "${SCRIPT_DIR}/merge_expression.py" "${OUT_DIR}/expr_tmp" "$EXPR_OUT"
fi

# 6.5. Copy pre-generated clinical files (timelines, etc.) to import folder
# Patient and sample clinical data are already merged above via merge_clinical.py.
# Supplemental files are copied here; timeline data files are sanitized below.
# The filenames are expected to already follow cBioPortal conventions.
if [ -n "$PRE_CLIN_DIR" ] && [ -d "$PRE_CLIN_DIR" ]; then
    echo "Copying pre-generated clinical files from ${PRE_CLIN_DIR}..."
    SKIP_BASENAMES="data_clinical_patient.txt data_clinical_patients.txt data_clinical_sample.txt data_clinical_samples.txt patient_clinical.tsv sample_clinical.tsv"
    copied=0
    for f in "${PRE_CLIN_DIR}"/data_*.txt "${PRE_CLIN_DIR}"/meta_*.txt; do
        [ -f "$f" ] || continue
        bname=$(basename "$f")
        # Skip the patient/sample files already handled by merge_clinical.py
        skip=false
        for skip_name in $SKIP_BASENAMES; do
            [ "$bname" = "$skip_name" ] && skip=true && break
        done
        $skip && continue
        cp "$f" "${IMPORT_DIR}/${bname}"
        echo "  Copied: ${bname}"
        copied=$(( copied + 1 ))
    done
    echo "  Total pre-clinical files copied: ${copied}"
fi

if ! python3 "${SCRIPT_DIR}/sanitize_clinical_timeline.py" "${IMPORT_DIR}" --apply; then
    echo "Error: Failed to sanitize clinical timeline files." >&2
    exit 1
fi

# Cleanup metadata files for optional data types if no actual data was generated
# (This prevents cBioPortal validator from failing due to missing data files)
# Remove legacy meta_expression.txt to prevent duplicate stable ID issues
rm -f "${IMPORT_DIR}/meta_expression.txt"

# Patient Clinical
if [ ! -f "${IMPORT_DIR}/data_clinical_patient.txt" ]; then
    echo "No patient clinical data found. Removing patient clinical metadata..."
    rm -f "${IMPORT_DIR}/meta_clinical_patient.txt"
fi

# Structural Variants
if [ ! -f "${IMPORT_DIR}/data_structural_variants.txt" ] || [ $(grep -v "^#" "${IMPORT_DIR}/data_structural_variants.txt" | wc -l) -le 1 ]; then
    echo "No structural variant data found. Removing structural variant files..."
    rm -f "${IMPORT_DIR}/data_structural_variants.txt" "${IMPORT_DIR}/meta_structural_variants.txt"
fi

# mRNA Expression
if [ ! -f "${IMPORT_DIR}/data_mrna_seq_v2_rsem.txt" ] || [ $(wc -l < "${IMPORT_DIR}/data_mrna_seq_v2_rsem.txt") -le 1 ]; then
    echo "No expression data found. Removing expression files..."
    rm -f "${IMPORT_DIR}/data_mrna_seq_v2_rsem.txt" "${IMPORT_DIR}/meta_mrna_seq_v2_rsem.txt"
fi

# Discrete CNA
if [ ! -f "${IMPORT_DIR}/data_CNA.txt" ] || [ $(wc -l < "${IMPORT_DIR}/data_CNA.txt") -le 1 ]; then
    echo "No CNA data found. Removing CNA files..."
    rm -f "${IMPORT_DIR}/data_CNA.txt" "${IMPORT_DIR}/meta_CNA.txt"
fi

# Segments
if [ ! -f "${IMPORT_DIR}/data_segment.txt" ] || [ $(wc -l < "${IMPORT_DIR}/data_segment.txt") -le 1 ]; then
    echo "No segment data found. Removing segment files..."
    rm -f "${IMPORT_DIR}/data_segment.txt" "${IMPORT_DIR}/meta_segment.txt" "${IMPORT_DIR}/data_segment.seg"
fi

# Mutations
if [ ! -f "${IMPORT_DIR}/data_mutations.txt" ] || [ $(grep -v "^#" "${IMPORT_DIR}/data_mutations.txt" | wc -l) -le 1 ]; then
    echo "No mutation data found. Removing mutation files..."
    rm -f "${IMPORT_DIR}/data_mutations.txt" "${IMPORT_DIR}/meta_mutations.txt"
fi

# 7. Generate WISP ctDNA Liquid Biopsy Resource Panel (if plasma data present)
WISP_INPUT="${IN_DIR:-${OUT_DIR}}"
if [ -n "$WISP_INPUT" ] && [ -d "$WISP_INPUT" ]; then
    echo "Checking for plasma ctDNA and WISP results in ${WISP_INPUT}..."
    python3 "${SCRIPT_DIR}/generate_wisp_cbio_resource.py" \
        -i "$WISP_INPUT" \
        -o "$IMPORT_DIR" \
        -s "$STUDY_ID"
fi

echo "Generating Case Lists..."
python3 "${SCRIPT_DIR}/generate_case_lists.py" "$IMPORT_DIR" "$STUDY_ID"

# 8. Generate Missing Metadata Files
echo "Generating missing metadata files..."
python3 "${SCRIPT_DIR}/generate_missing_meta.py" "$IMPORT_DIR" "$STUDY_ID"

chmod -R g+rwX "$OUT_DIR" 2>/dev/null || true
echo "Study merge complete."

