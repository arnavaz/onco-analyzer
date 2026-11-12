#!/bin/bash
# Generate cBioPortal metadata files
# Usage: ./generate_meta.sh <out_dir> <study_id> <study_name> <study_desc> [type_of_cancer]

OUT_DIR=$1
STUDY_ID=$2
STUDY_NAME=$3
STUDY_DESC=$4
TYPE_OF_CANCER=${5:-mixed}

# meta_study.txt
if [ ! -f "${OUT_DIR}/meta_study.txt" ]; then
cat <<EOF > "${OUT_DIR}/meta_study.txt"
type_of_cancer: ${TYPE_OF_CANCER}
cancer_study_identifier: ${STUDY_ID}
name: ${STUDY_NAME}
description: ${STUDY_DESC}
reference_genome: hg38
add_global_case_list: true
EOF
else
    echo "Existing meta_study.txt found in ${OUT_DIR}. Skipping generation to preserve manual changes."
fi

# meta_clinical_patient.txt
cat <<EOF > "${OUT_DIR}/meta_clinical_patient.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: CLINICAL
datatype: PATIENT_ATTRIBUTES
data_filename: data_clinical_patient.txt
EOF

# meta_clinical_sample.txt
cat <<EOF > "${OUT_DIR}/meta_clinical_sample.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
EOF

# meta_mutations.txt
cat <<EOF > "${OUT_DIR}/meta_mutations.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: MUTATION_EXTENDED
datatype: MAF
stable_id: mutations
show_profile_in_analysis_tab: true
profile_name: Mutations
profile_description: Mutation data.
data_filename: data_mutations.txt
EOF

# meta_CNA.txt
cat <<EOF > "${OUT_DIR}/meta_CNA.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: DISCRETE
stable_id: cna
show_profile_in_analysis_tab: true
profile_name: Putative copy-number alterations
profile_description: Putative copy-number alterations. Values: -2 = homozygous deletion (minCopyNumber < 0.5); -1 = heterozygous deletion / LOH (minMinorAlleleCopyNumber < 0.5); 0 = neutral; 1 = gain (maxCopyNumber > 2.0 x ploidy); 2 = high-level amplification (maxCopyNumber > 3.0 x ploidy).
data_filename: data_CNA.txt
EOF

# meta_segment.txt
cat <<EOF > "${OUT_DIR}/meta_segment.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: SEG
reference_genome_id: hg38
description: Copy number segment data.
data_filename: data_segment.txt
EOF

# meta_structural_variants.txt
cat <<EOF > "${OUT_DIR}/meta_structural_variants.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: STRUCTURAL_VARIANT
datatype: SV
stable_id: structural_variants
show_profile_in_analysis_tab: true
profile_name: Structural variants
profile_description: Structural variant data.
data_filename: data_structural_variants.txt
EOF

# meta_mrna_seq_v2_rsem.txt
cat <<EOF > "${OUT_DIR}/meta_mrna_seq_v2_rsem.txt"
cancer_study_identifier: ${STUDY_ID}
genetic_alteration_type: MRNA_EXPRESSION
datatype: CONTINUOUS
stable_id: rna_seq_v2_mrna
show_profile_in_analysis_tab: true
profile_name: mRNA expression (TPM)
profile_description: Expression levels in Transcripts Per Million (TPM).
data_filename: data_mrna_seq_v2_rsem.txt
EOF

echo "Generated metadata files in ${OUT_DIR}"
