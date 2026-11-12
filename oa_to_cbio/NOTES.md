# cBioPortal Conversion Strategy (oncoanalyser/HMF)

## Purpose
This file captures the recommended strategy for converting nf-core/oncoanalyser outputs into cBioPortal-compatible study data.

## Core Conclusions
- No widely adopted public wrapper was identified that directly converts oncoanalyser/HMF outputs to cBioPortal import files end-to-end.
- The most robust path is a standards-based conversion using cBioPortal-native schemas (MAF, SV, CNA, clinical).
- Build and validate in phases: start with mutations, then add SV/CNA/expression and report-derived metadata.

## Preferred Data Sources by Data Type
- Somatic SNV/Indel (final callset): `purple/*.purple.somatic.vcf.gz`
- Germline SNV/Indel (if importing): `purple/*.purple.germline.vcf.gz`
- Raw somatic small variants (alternative/raw): `sage/somatic/*.sage.somatic.vcf.gz`
- Functional annotation of SAGE calls: `pave/somatic/*.pave.somatic.vcf.gz`
- Structural variants: `purple/*.purple.sv.vcf.gz` and/or `linx` outputs
- Copy number: PURPLE CNV/segment TSV outputs
- Expression (optional): ISOFOX outputs

## ORANGE JSON Guidance
- Treat ORANGE JSON as enrichment, not the primary variant source.
- Best use:
  - Populate sample/patient-level clinical attributes (purity, ploidy, TMB, MSI, HRD/WGD flags, summary interpretation fields).
  - Add report-level interpretation metadata.
- Do not rely on ORANGE JSON alone for variant-level import because cBioPortal expects structured mutation/SV/CNA files.

## Recommended Implementation Phases

### Phase 1: Minimal Viable Import (Highest ROI)
1. Convert somatic VCF to MAF (use `vcf2maf`).
2. Create required cBioPortal clinical files (patient/sample).
3. Validate dataset with cBioPortal validator/import checks.
4. Import mutation + clinical data first.

#### Phase 1 Deliverables (exact files)
- `meta_study.txt`
- `meta_clinical_patient.txt`
- `data_clinical_patient.txt`
- `meta_clinical_sample.txt`
- `data_clinical_sample.txt`
- `meta_mutations.txt`
- `data_mutations.txt` (MAF)

#### Phase 1 Required minimum fields
- `meta_study.txt`: `cancer_study_identifier`, `type_of_cancer`, `name`, `description`, `short_name`
- `data_clinical_patient.txt`: `PATIENT_ID` (+ optional demographics)
- `data_clinical_sample.txt`: `SAMPLE_ID`, `PATIENT_ID` (+ optional sample attributes)
- `data_mutations.txt` (MAF): `Tumor_Sample_Barcode`, `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`, `Variant_Type`

#### Phase 1 Mapping notes
- Small-variant source priority: `purple/*.purple.somatic.vcf.gz` -> MAF.
- Keep a stable `SAMPLE_ID` policy and apply it to all files before import.

### Phase 2: Enriched Genomics
1. Add SV file(s) mapped to cBioPortal SV schema.
2. Add CNA data (SEG and/or discrete CNA matrix).
3. Optionally add expression data.
4. Re-validate and import incrementally.

#### Phase 2 Deliverables (additive)
- `meta_structural_variants.txt`
- `data_structural_variants.txt`
- `meta_CNA.txt`
- `data_CNA.txt` (discrete CNA matrix) and/or
- `meta_segment.txt`
- `data_segment.txt` (SEG)
- Optional expression:
  - `meta_mrna_seq_v2_rsem.txt`
  - `data_mrna_seq_v2_rsem.txt`

#### Phase 2 Source mapping
- SV: `purple/*.purple.sv.vcf.gz` and/or LINX outputs mapped to cBioPortal SV table format.
- CNA: PURPLE CNV/segment TSV outputs mapped to discrete CNA and/or SEG.
- Expression: ISOFOX tables mapped to cBioPortal expression matrix.

#### Phase 2 Required consistency checks
- Every `SAMPLE_ID` in SV/CNA/expression exists in `data_clinical_sample.txt`.
- Gene symbols normalized consistently across MAF, CNA, expression.
- Genome build is consistent across all modalities.

### Phase 3: ORANGE-Derived Enrichment
1. Extract selected ORANGE JSON fields into clinical/sample attributes.
2. Keep original variant-level data in MAF/SV/CNA files.
3. Avoid duplicating variant semantics already represented in core genomic files.

#### Phase 3 Deliverables (attribute enrichment)
- Update existing clinical files:
  - `data_clinical_sample.txt`
  - `data_clinical_patient.txt`
- Optional custom profile for report-level values:
  - `meta_clinical_report.txt` (if using a separate profile)
  - `data_clinical_report.txt`

#### Phase 3 Recommended ORANGE-derived fields
- Sample-level: tumor purity, ploidy, TMB, MSI status, WGD flag.
- Patient/sample-level summaries: HRD status, key interpretation flags, report QC summaries.
- Keep actionability tiers and narrative interpretation as additive attributes, not replacements for genomic event tables.

## Import Bundle Layout (recommended)
- `study/`
  - `meta_study.txt`
  - `meta_clinical_patient.txt`
  - `data_clinical_patient.txt`
  - `meta_clinical_sample.txt`
  - `data_clinical_sample.txt`
  - `meta_mutations.txt`
  - `data_mutations.txt`
  - `meta_structural_variants.txt` (phase 2+)
  - `data_structural_variants.txt` (phase 2+)
  - `meta_CNA.txt` (phase 2+)
  - `data_CNA.txt` (phase 2+)
  - `meta_segment.txt` (phase 2+)
  - `data_segment.txt` (phase 2+)
  - `meta_mrna_seq_v2_rsem.txt` (optional)
  - `data_mrna_seq_v2_rsem.txt` (optional)

## Validation Gates by Phase
- Phase 1 gate: mutation + clinical import succeeds and portal pages render sample-level mutation views correctly.
- Phase 2 gate: SV/CNA/expression profiles load with no orphan sample IDs and expected case counts.
- Phase 3 gate: ORANGE-derived attributes are queryable/filterable and do not duplicate core event tables.

## Practical Notes
- Maintain strict sample ID consistency across all files (most common source of import failure).
- Preserve important HMF-specific annotations as extra columns where supported.
- Prefer iterative pilot imports on a small cohort before full-scale conversion.

## Tooling to Use
- `vcf2maf` for VCF -> MAF conversion.
- cBioPortal official file format docs and validator/import tools.
- Custom scripts for mapping LINX/PURPLE/ISOFOX outputs into cBioPortal tabular schemas.

## What Not to Do
- Do not attempt a one-shot direct import of raw oncoanalyser folders.
- Do not use ORANGE JSON as a replacement for mutation/SV/CNA files.
- Do not mix inconsistent sample identifiers across modalities.

# User Implementation Guidance and Comments
# ------------------------------------------
# - vcf2maf can be accessed on the cluster using: module load vcf2maf/1.6.22
# - Use SLURM job dependencies so that all vcf2maf jobs finish before combining MAFs
# - Add a parameter to set the SLURM account (-A [account]) and propagate to all jobs
# - All MAFs should be combined automatically once generated
# - Use bash scripts and SLURM for parallelization, not Python
# - Clinical file generation and other downstream steps should depend on successful MAF combination
