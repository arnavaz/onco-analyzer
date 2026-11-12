# Oncoanalyser to cBioPortal Pipeline

This repository contains a set of Bash and R scripts orchestrated via SLURM to convert standard nf-core/oncoanalyser outputs into a cBioPortal-ready study import bundle. 

## Features
- **Phase 1 (Core Imports):** Automatically generates clinical metadata files (`meta_study.txt`, `data_clinical_patient.txt`, `data_clinical_sample.txt`) and converts PURPLE somatic VCFs to MAF format using `vcf2maf`.
- **Phase 2 (Copy Number & Segments):** Parses PURPLE discrete CNA (`cnv.gene.tsv`) and segment (`segment.tsv`) files into cBioPortal-compatible matrices and lists.
- **Structural Variants:** Includes a placeholder R script (`process_sv.R`) ready to convert LINX SV outputs to the cBioPortal SV table format.

## Directory Structure

```text
/Users/jbruce/bin/oa_to_cbio/
├── oa_to_cbio.sh               # Master pipeline orchestrator
├── README.md                   # This documentation
├── NOTES.md                    # Original strategic notes
├── examples/                   # Directory containing example tree structures
└── scripts/
    ├── generate_meta.sh        # Generates cBioPortal meta_*.txt files
    ├── process_sample.sh       # SLURM script for converting a single sample (vcf2maf, CNA, SEG, SVs)
    ├── process_sv.R            # Placeholder R script to process SVs
    └── merge_study.sh          # SLURM script for joining all single-sample files together
```

## How It Works

1. **Initialization (`oa_to_cbio.sh`)**: 
   The orchestrator script creates the necessary output directories and delegates the generation of the `meta_*.txt` boilerplate files to `generate_meta.sh`.
   
2. **Clinical Data Harvesting**: 
   The orchestrator scans through your provided patient folders, determines the Tumor and Normal sample IDs based on the `.redux.bam` alignments, and constructs the required `data_clinical_patient.txt` and `data_clinical_sample.txt` tables.

3. **Parallel Sample Processing (`process_sample.sh`)**: 
   For every tumor sample, `oa_to_cbio.sh` dispatches a SLURM job. This parallel job:
   - Uses `vcf2maf.pl` (via `module load vcf2maf/1.6.22`) to build the MAF file from the somatic VCF.
   - Restructures the PURPLE CNV output into discrete CNA vectors (`-2`, `-1`, `0`, `1`, `2`) based on the thresholding of `minCopyNumber`.
   - Modifies the PURPLE Segment output into the specific coordinate schema cBioPortal expects (`log2(CN/2)`).
   - Executes the placeholder `process_sv.R` script for Structural Variants.

4. **Study Assembly (`merge_study.sh`)**: 
   A final merge SLURM job is submitted with an `afterok` dependency on all sample jobs. Once every sample job finishes, this script collects and consolidates the temporary MAFs, CNA vectors, and Segment lists into single combined `data_mutations.txt`, `data_CNA.txt`, and `data_segment.txt` files.

## Structural Variants (SV) Parsing

The structural variant processing is currently routed to `scripts/process_sv.R`. This file accepts `linx.svs.tsv` (or `purple.sv.vcf.gz`) as input. To implement SV conversions, drop your existing R logic into this file so that it maps the breakends to cBioPortal's required columns (e.g., `Site1_Hugo_Symbol`, `Site2_Hugo_Symbol`, `SV_Class`). 

## Usage Example

To launch the pipeline, run the master orchestrator `oa_to_cbio.sh` from the head node. You must provide the input directory, the desired output directory, and your SLURM account string.

```bash
cd /Users/jbruce/bin/oa_to_cbio/

./oa_to_cbio.sh \
  -i /cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/ \
  -o /path/to/desired/cbio_import_dir \
  -A your_slurm_account \
  -s "pancan_h" \
  -n "PANCAN-H Example Study" \
  -d "Processed from nf-core/oncoanalyser"
```

## Validation
Before fully importing into your cBioPortal instance, it is highly recommended to run the resulting output directory through the official [cBioPortal Data Validator](https://github.com/cBioPortal/cbioportal/tree/master/core/src/main/scripts/importer) (`validateData.py`) to ensure syntax, headers, and mapping are fully compliant.
