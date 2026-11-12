#!/bin/bash
# SLURM script to process a single sample
# Usage: sbatch -A <account> --export=... process_sample.sh

# Environment variables expected:
# TUMOR_ID, NORMAL_ID, SAMPLE_DIR, OUT_DIR

umask 0002

module load vcf2maf/1.6.22
module load samtools/1.20  # Provides tabix
module load R/4.1.0         # Provides Rscript

# Ensure tmp dirs exist
mkdir -p "${OUT_DIR}/maf_tmp"
mkdir -p "${OUT_DIR}/cna_tmp"
mkdir -p "${OUT_DIR}/seg_tmp"
mkdir -p "${OUT_DIR}/sv_tmp"
mkdir -p "${OUT_DIR}/clinical_tmp"
mkdir -p "${OUT_DIR}/expr_tmp"

echo "Processing $TUMOR_ID"

SOMATIC_VCF_IN="${SAMPLE_DIR}/purple/${TUMOR_ID}.purple.somatic.vcf.gz"
GERMLINE_VCF_IN="${SAMPLE_DIR}/purple/${TUMOR_ID}.purple.germline.vcf.gz"
ORANGE_JSON="${SAMPLE_DIR}/orange/${TUMOR_ID}.orange.json"

MAF_OUT="${OUT_DIR}/maf_tmp/${TUMOR_ID}.maf"
SOMATIC_MAF="${OUT_DIR}/maf_tmp/${TUMOR_ID}.somatic.maf"
GERMLINE_MAF="${OUT_DIR}/maf_tmp/${TUMOR_ID}.germline.maf"
CNV_IN="${SAMPLE_DIR}/purple/${TUMOR_ID}.purple.cnv.gene.tsv"
SEG_IN="${SAMPLE_DIR}/purple/${TUMOR_ID}.purple.segment.tsv"

# Handle legacy resume
if [ -s "$MAF_OUT" ] && [ ! -s "$SOMATIC_MAF" ] && [ ! -s "${SOMATIC_MAF}.unannotated" ]; then
    echo "Legacy MAF found, migrating to somatic MAF..."
    mv "$MAF_OUT" "${SOMATIC_MAF}.unannotated"
fi

# Function to run vcf2maf and annotate
process_vcf() {
    local vcf_in="$1"
    local maf_final="$2"
    local status="$3"
    local vcf_type="$4"

    if [ -f "$vcf_in" ]; then
        if [ ! -s "${maf_final}.unannotated" ]; then
            if [ -s "$maf_final" ]; then
                echo "Found existing $vcf_type MAF. Skipping VEP and renaming for re-annotation..."
                mv "$maf_final" "${maf_final}.unannotated"
            else
                echo "Running vcf2maf for $vcf_type..."
                TMP_VCF="${OUT_DIR}/maf_tmp/${TUMOR_ID}.purple.${vcf_type}.vcf"
                zcat "$vcf_in" > "$TMP_VCF"
                
                rm -f "${TMP_VCF%.vcf}.vep.vcf"
                vcf2maf.pl \
                    --input-vcf "$TMP_VCF" \
                    --output-maf "${maf_final}.unannotated" \
                    --tumor-id "$TUMOR_ID" \
                    --normal-id "$NORMAL_ID" \
                    --vep-path /cluster/tools/software/centos7/vep/112 \
                    --vep-data /cluster/projects/pughlab/references/VEP_cache/112 \
                    --ref-fasta /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta \
                    --ncbi-build GRCh38
                    
                rm -f "$TMP_VCF"
            fi
        fi
        
        if [ -s "${maf_final}.unannotated" ]; then
            echo "Annotating $vcf_type MAF with status $status..."
            local annot_cmd=("python3" "${SCRIPT_DIR}/annotate_maf.py" "${maf_final}.unannotated" "$ORANGE_JSON" "$maf_final" "$status")
            if [ "$UNFILTERED_MAF" = "true" ]; then
                annot_cmd+=("--unfiltered")
            fi
            if [ "${NON_GERMLINE_FILTER:-false}" = "true" ]; then
                annot_cmd+=("--non-germline-filter")
            fi
            "${annot_cmd[@]}"
            rm -f "${maf_final}.unannotated"
        fi
    else
        echo "Warning: No $vcf_type VCF found at $vcf_in"
    fi
}

# 1. vcf2maf
process_vcf "$SOMATIC_VCF_IN" "$SOMATIC_MAF" "Somatic" "somatic"
process_vcf "$GERMLINE_VCF_IN" "$GERMLINE_MAF" "Germline" "germline"

# Merge MAFs
rm -f "$MAF_OUT"
if [ -s "$SOMATIC_MAF" ]; then
    cat "$SOMATIC_MAF" > "$MAF_OUT"
fi
if [ -s "$GERMLINE_MAF" ]; then
    if [ -s "$MAF_OUT" ]; then
        grep -v "^#" "$GERMLINE_MAF" | tail -n +2 >> "$MAF_OUT"
    else
        cat "$GERMLINE_MAF" > "$MAF_OUT"
    fi
fi

# 1.25 Optional: SAGE somatic VCF to unfiltered germline candidates MAF
SAGE_VCF_IN="${SAMPLE_DIR}/sage/somatic/${TUMOR_ID}.sage.somatic.vcf.gz"
UNFILTERED_SAGE_MAF="${OUT_DIR}/maf_tmp/${TUMOR_ID}.unfiltered_sage.maf"

if [ "${UNFILTERED_GERMLINE_MAF:-false}" = "true" ] && [ -f "$SAGE_VCF_IN" ]; then
    if [ ! -s "$UNFILTERED_SAGE_MAF" ]; then
        echo "Running vcf2maf for SAGE somatic VCF (unfiltered germline candidates)..."
        TMP_VCF="${OUT_DIR}/maf_tmp/${TUMOR_ID}.sage.somatic.vcf"
        zcat "$SAGE_VCF_IN" > "$TMP_VCF"
        
        rm -f "${TMP_VCF%.vcf}.vep.vcf"
        vcf2maf.pl \
            --input-vcf "$TMP_VCF" \
            --output-maf "$UNFILTERED_SAGE_MAF" \
            --tumor-id "$TUMOR_ID" \
            --normal-id "$NORMAL_ID" \
            --vep-path /cluster/tools/software/centos7/vep/112 \
            --vep-data /cluster/projects/pughlab/references/VEP_cache/112 \
            --ref-fasta /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta \
            --ncbi-build GRCh38
            
        rm -f "$TMP_VCF"
    else
        echo "Found existing unfiltered_sage MAF at $UNFILTERED_SAGE_MAF. Skipping..."
    fi
fi

# 1.5 Clinical Data Extraction
CLINICAL_OUT="${OUT_DIR}/clinical_tmp/${TUMOR_ID}.clinical.tsv"
ORANGE_JSON="${SAMPLE_DIR}/orange/${TUMOR_ID}.orange.json"

# Auto-detect matched normal ID from TEAL outputs if unset
if [ -z "$NORMAL_ID" ] || [ "$NORMAL_ID" = "UNKNOWN_NORMAL" ]; then
    DETECTED_NORMAL=$(ls "${SAMPLE_DIR}/teal/"*.teal.tellength.tsv 2>/dev/null | grep -v "/${TUMOR_ID}\.teal\." | head -n 1 || true)
    if [ -n "$DETECTED_NORMAL" ]; then
        NORMAL_ID=$(basename "$DETECTED_NORMAL" | sed 's/\.teal\.tellength\.tsv$//')
    fi
fi

python3 "${SCRIPT_DIR}/process_clinical.py" "$TUMOR_ID" "$ORANGE_JSON" "$CLINICAL_OUT" --normal-id "$NORMAL_ID"

# Extract Ploidy
PLOIDY=2.0
if [ -f "$CLINICAL_OUT" ]; then
    EXTRACTED_PLOIDY=$(awk -F'\t' '$1=="PLOIDY" {print $2}' "$CLINICAL_OUT")
    if [[ "$EXTRACTED_PLOIDY" != "NA" ]] && [[ "$EXTRACTED_PLOIDY" != "" ]]; then
        PLOIDY=$EXTRACTED_PLOIDY
    fi
fi

# 2. Discrete CNA
if [ -f "$CNV_IN" ]; then
    python3 "${SCRIPT_DIR}/parse_cna.py" "$CNV_IN" "$TUMOR_ID" "$PLOIDY" "${AMP_RATIO:-3.0}" "${GAIN_RATIO:-2.0}" "${OUT_DIR}/cna_tmp/${TUMOR_ID}.cna.tsv"
fi

# 3. Segments
if [ -f "$SEG_IN" ]; then
    awk -F'\t' -v sample="$TUMOR_ID" 'BEGIN {OFS="\t"} 
    NR==1 {
        cn_col = 0
        cn_col_is_ratio = 0
        status_col = 0
        for(i=1; i<=NF; i++) {
            h = $i
            if(h == "chromosome" || h == "chrom") chrom_col = i
            if(h == "start") start_col = i
            if(h == "end") end_col = i
            if(h == "bafCount" || h == "baf_count") baf_col = i
            if(h == "germlineStatus" || h == "germline_status") status_col = i
            if(h == "observedTumorRatio" || h == "observed_tumor_ratio") {
                cn_col = i
                cn_col_is_ratio = 1
            }
        }
        # Fallback to absolute copy number if raw ratio is not found
        if(cn_col == 0) {
            for(i=1; i<=NF; i++) {
                h = $i
                if(h == "copyNumber" || h == "copy_number" || h == "tumorCopyNumber") {
                    cn_col = i
                }
            }
        }
        print "ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"
        next
    } 
    {
        chrom = (chrom_col ? $chrom_col : $1)
        start = (start_col ? $start_col : $2)
        end = (end_col ? $end_col : $3)
        
        status = (status_col ? $status_col : "")
        if(status == "CENTROMETIC" || status == "EXCLUDED" || status == "NOISE") next
        
        marks = 100
        if(baf_col && $baf_col != "" && $baf_col != "NA") {
            val = $baf_col + 0
            if(val == $baf_col && val > 0) {
                marks = int(val)
            }
        }
        
        cn = ""
        if(cn_col) cn = $cn_col
        if(cn == "" || cn == "NA") next
        
        # Strip "chr" prefix
        sub(/^chr/, "", chrom)
        
        # Skip zero-width or empty/invalid segments
        if(start == end || start == "" || end == "") next
        
        if(cn_col_is_ratio) {
            if(cn <= 0.25) segmean = -2.0
            else segmean = log(cn)/log(2)
        } else {
            if(cn <= 0) segmean = -2.0
            else segmean = log(cn/2)/log(2)
        }
        
        print sample, chrom, start, end, marks, segmean
    }' "$SEG_IN" > "${OUT_DIR}/seg_tmp/${TUMOR_ID}.seg.tsv"
fi

# 4. SVs
SV_VIS_IN="${SAMPLE_DIR}/linx/somatic_annotations/${TUMOR_ID}.linx.vis_sv_data.tsv"
SV_SVS_IN="${SAMPLE_DIR}/linx/somatic_annotations/${TUMOR_ID}.linx.svs.tsv"
SV_BND_IN="${SAMPLE_DIR}/linx/somatic_annotations/${TUMOR_ID}.linx.breakend.tsv"
FUSION_IN=$(ls "${SAMPLE_DIR}/isofox/"*.isf.fusions.{csv,tsv} 2>/dev/null | head -n 1 || true)
SV_VCF_IN="${SAMPLE_DIR}/purple/${TUMOR_ID}.purple.sv.vcf.gz"

if [ -f "$SV_VIS_IN" ] && [ -f "$SV_SVS_IN" ] && [ -f "$SV_BND_IN" ]; then
    Rscript "${SCRIPT_DIR}/process_sv.R" "$SV_VIS_IN" "$SV_SVS_IN" "$SV_BND_IN" "$FUSION_IN" "$SV_VCF_IN" "${OUT_DIR}/sv_tmp/${TUMOR_ID}.sv.tsv" "$TUMOR_ID"
else
    echo "Warning: Missing LINX SV files for $TUMOR_ID"
fi

# 5. Expression Data
# Support both .csv (older local runs) and .tsv (newer GCP/OA3 outputs)
EXPR_IN=$(ls "${SAMPLE_DIR}/isofox/"*.isf.gene_data.{csv,tsv} 2>/dev/null | head -n 1 || true)

# If not found in isofox/, try the root SAMPLE_DIR (some HMFTools versions flatten this)
if [ -z "$EXPR_IN" ] || [ ! -f "$EXPR_IN" ]; then
    EXPR_IN=$(ls "${SAMPLE_DIR}/"*.isf.gene_data.{csv,tsv} 2>/dev/null | head -n 1 || true)
fi

if [ -n "$EXPR_IN" ] && [ -f "$EXPR_IN" ]; then
    # Extract GeneName, GeneId, and AdjTPM
    # Auto-detect separator: TSV files from OA3/GCP use tab; older CSV files use comma.
    # Column layout is identical: col1=GeneId, col2=GeneName, col10=AdjTPM
    awk 'BEGIN{FS=","} NR==1{if(index($0,"\t")>0) FS="\t"} NR>1{print $2 "\t" $1 "\t" $10}' "$EXPR_IN" > "${OUT_DIR}/expr_tmp/${TUMOR_ID}.expr.tsv"
else
    echo "Warning: Missing Isofox expression file for $TUMOR_ID"
fi

# 6. Generate Gene-focused BAF Plot
echo "Generating Gene-focused BAF Plot for $TUMOR_ID..."
mkdir -p "${OUT_DIR}/plots"
PLOT_MAF_IN=""
if [ -s "${OUT_DIR}/maf_tmp/${TUMOR_ID}.unfiltered_sage.maf" ]; then
    PLOT_MAF_IN="${OUT_DIR}/maf_tmp/${TUMOR_ID}.unfiltered_sage.maf"
elif [ -s "${OUT_DIR}/maf_tmp/${TUMOR_ID}.maf" ]; then
    PLOT_MAF_IN="${OUT_DIR}/maf_tmp/${TUMOR_ID}.maf"
fi

if [ -n "$PLOT_MAF_IN" ]; then
    DRIVERS_TSV="${DRIVER_GENE_PANEL:-/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/final_panel_files/MMR_drivers.tsv}"
    if [ -f "$DRIVERS_TSV" ]; then
        baf_cmd=("python3" "${SCRIPT_DIR}/plot_gene_baf.py" \
            --maf "$PLOT_MAF_IN" \
            --genes-file "$DRIVERS_TSV" \
            --output "${OUT_DIR}/plots/${TUMOR_ID}_gene_baf.png" \
            --title "Gene BAF Plot: ${TUMOR_ID}")
        
        if [ -f "$SEG_IN" ]; then
            baf_cmd+=(--purple-segments "$SEG_IN")
        fi
        if [ -n "${GTF_FILE:-}" ] && [ -f "$GTF_FILE" ]; then
            baf_cmd+=(--gtf "$GTF_FILE")
        fi
        
        "${baf_cmd[@]}"
    else
        echo "Warning: MMR_drivers.tsv not found at $DRIVERS_TSV. Skipping BAF plotting."
    fi
else
    echo "Warning: No MAF file found for BAF plotting."
fi

echo "Done processing $TUMOR_ID"
