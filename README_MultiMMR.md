# MultiMMR Panel Workflow & Resource Generation Guide

This guide describes how to isolate this workflow into a clean, dedicated repository, run the MultiMMR panel workflow, compile cohort MMR reports, and generate custom panel resources from scratch when the panel capture design changes.

---

## 1. Up-Front Commands to Create a Clean Split Repository & Push to GitHub

To share this repository with a group outside of Pugh Lab without including extraneous files (like test runs, cBioPortal scripts, or WGTS wrappers), you can use the `create_clean_repo.sh` utility to generate a clean, streamlined repository. 

Follow the steps below:

### Step 1: Create the Clean Repository Local Copy
Run the `create_clean_repo.sh` script, providing a path where you want the new repository to be created (e.g., in a sibling folder `oncoanalyser_multimmr`):

```bash
# Run the script to extract only core pipeline files, Python scripts, wrappers, and MultiMMR v3 resources
bash create_clean_repo.sh /path/to/oncoanalyser_multimmr
```

### Step 2: Push the Clean Repository to GitHub
Once created, navigate into the new directory and push it to a new GitHub repository:

#### Option A: Using the GitHub CLI (`gh`) (Recommended)
```bash
cd /path/to/oncoanalyser_multimmr

# Authenticate with GitHub (if not already logged in)
gh auth login

# Create a private/public repo on GitHub and push the code directly
gh repo create oncoanalyser_multimmr --private --source=. --remote=origin --push
```

#### Option B: Using standard Git CLI
1. Go to [GitHub](https://github.com/) and create a new repository (e.g., `oncoanalyser_multimmr`). Do **not** initialize it with any files.
2. Run the following commands in the new directory:

```bash
cd /path/to/oncoanalyser_multimmr

# Add the new repository as origin remote
git remote add origin git@github.com:<YOUR-GITHUB-USERNAME-OR-ORG>/oncoanalyser_multimmr.git
# Or using HTTPS:
# git remote add origin https://github.com/<YOUR-GITHUB-USERNAME-OR-ORG>/oncoanalyser_multimmr.git

# Set main branch and push
git branch -M main
git push -u origin main
```

---

## 2. Packaged MultiMMR v3 Panel Resources

The MultiMMR v3 panel files are packaged directly within this repository to enable out-of-the-box execution:
- **Location:** `resources/panels/multimmr_v3/`
- **Configuration File:** [panel_data.config](file:///Users/jbruce/bin/oncoanalyser_pughlab/resources/panels/multimmr_v3/panel_data.config)

The configuration file uses relative paths utilizing the Nextflow `${projectDir}` parameter, which automatically resolves paths correctly relative to where Nextflow runs:

```nextflow
params {
  ref_data_panel_data_path = "${projectDir}/../resources/panels/multimmr_v3"
  ...
}
```

---

## 3. Running the MultiMMR Pipeline

### Step 1: Generate Patient Samplesheets
Use the helper script to create per-patient samplesheets from your input FASTQ files and a mapping sheet:

```bash
module load python3/3.12.11
python3 Setup/make_oncoanalyser_samplesheet.py \
  --fastq-dir /path/to/panel_fastqs \
  --mapping-tsv /path/to/injection_sheet.tsv \
  --out /path/to/sheets_per_patient
```

### Step 2: Submit a Patient Run
To run targeted analysis on a single patient using the packaged MultiMMR v3 panel, use the `run_oncoanalyser.sh` SLURM submission wrapper:

```bash
SOURCE_REPO=/path/to/your/cloned/repo
sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A <your-slurm-account> \
  ${SOURCE_REPO}/run_oncoanalyser.sh \
  /path/to/sheets_per_patient/PATIENT_X/samplesheet.csv \
  /path/to/results/PATIENT_X \
  targeted \
  GRCh38_hmf \
  <your-slurm-account> \
  --panel-name multimmr_hmf \
  --panel-config ${SOURCE_REPO}/resources/panels/multimmr_v3/panel_data.config \
  --amber-expanded-snps \
  --force-panel
```

---

## 4. Upstream MMR Summary & Cohort Compilation

### Patient-Level Summarization (Upstream)
When `run_oncoanalyser.sh` runs targeted mode, it automatically triggers [scripts/summarize_findings.py](file:///Users/jbruce/bin/oncoanalyser_pughlab/scripts/summarize_findings.py) on completion. This script reads the genomic outputs (copy number, structural variants, somatic mutations, germline modifications, and methylation matrix files) and writes a local MMR summary for the sample:
- `/path/to/results/PATIENT_X/final_summary.tsv`
- `/path/to/results/PATIENT_X/summary/PATIENT_X_clinical_mmr_report.html`

### Cohort-Level Summary (Downstream)
Once all patient runs in a cohort have completed, compile their findings into a single cohort HTML dashboard using [scripts/compile_mmr_summary.py](file:///Users/jbruce/bin/oncoanalyser_pughlab/scripts/compile_mmr_summary.py):

```bash
python3 scripts/compile_mmr_summary.py \
  --results-dir /path/to/results_base
```

This will output `cohort_mmr_dashboard.html` in the `/path/to/results_base` directory. The dashboard displays the clinical classification, MSI status, and individual MMR gene risk status for all patients, and includes clickable links opening each patient's individual `clinical_mmr_report.html`.

---

## 5. Generating Panel Resources from Scratch (For Custom Panels)

If the capture panel changes (e.g., v4), you must generate new panel resources. Follow these steps:

### Step 1: Prepare the Target Regions BED File
Sort and compress the design/probes BED file representing the capture target coordinates:
```bash
sort -k1,1 -k2,2n target_bases.bed | bgzip -c > my_panel_target_regions.38.bed.gz
tabix -p bed my_panel_target_regions.38.bed.gz
```

### Step 2: Generate SNP Panel Augmentation (Amber Loci)
To enable reliable copy number profiling from panel data, expand the standard Amber heterozygous sites file to include high-frequency SNPs within the panel's targeted regions. Run [scripts/generate_panel_snps.py](file:///Users/jbruce/bin/oncoanalyser_pughlab/scripts/generate_panel_snps.py):

```bash
python3 scripts/generate_panel_snps.py \
  --target-regions-bed my_panel_target_regions.38.bed.gz \
  --default-loci /path/to/reference_data/dna/copy_number/AmberGermlineSites.38.tsv.gz \
  --output AmberGermlineSites.expanded.38.tsv.gz \
  --genes all_bed \
  --genome GRCh38_hmf
```

### Step 3: Run Panel Resource Creation Workflow
Create a merged cohort-level samplesheet (at least 20 samples recommended) and run oncoanalyser in `panel_resource_creation` mode. This will generate panel-specific somatic artifact and normalisation files:

```bash
SOURCE_REPO=/path/to/your/cloned/repo
sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A <your-slurm-account> \
  ${SOURCE_REPO}/run_oncoanalyser.sh \
  /path/to/samplesheet.panel_resource_creation.csv \
  /path/to/panel_resource_output \
  panel_resource_creation \
  GRCh38_hmf \
  <your-slurm-account> \
  --driver-gene-panel /path/to/DriverGenePanel.38.tsv \
  --target-regions-bed my_panel_target_regions.38.bed.gz
```

Once complete, retrieve the generated files from the output directory:
- `cobalt.region_normalisation.38.tsv` (found in `panel_resource_output/panel_resources/`)
- `pave.somatic_artefacts.38.tsv` (found in `panel_resource_output/panel_resources/`)

### Step 4: Generate MSI Target Sites (Optional)
If you need to scan targeted BED intervals for homopolymer repeats to run MSI evaluations, run [scripts/generate_msi_panel.py](file:///Users/jbruce/bin/oncoanalyser_pughlab/scripts/generate_msi_panel.py):

```bash
python3 scripts/generate_msi_panel.py \
  --target-regions-bed my_panel_target_regions.38.bed.gz \
  --msi-database /path/to/msi_jitter_sites.38.tsv.gz \
  --min-mappability 0.9 \
  --output my_panel_msi_sites.tsv
```

### Step 5: Write the Panel Configuration
Put all generated resources in a single directory (e.g. `resources/panels/my_custom_panel/`) and write a `panel_data.config` file mapping them:

```nextflow
params {
  ref_data_panel_data_path = "${projectDir}/../resources/panels/my_custom_panel"

  panel_data_paths {
    my_custom_panel {
      "38" {
        driver_gene_panel           = "DriverGenePanel.38.tsv"
        pon_artefacts               = "pave.somatic_artefacts.38.tsv"
        target_region_bed           = "my_panel_target_regions.38.bed.gz"
        target_region_normalisation = "cobalt.region_normalisation.38.tsv"
        target_region_ratios        = "target_regions_ratios.38.tsv"
        target_region_msi_indels    = "my_panel_msi_sites.tsv"
        amber_heterozygous_sites    = "AmberGermlineSites.expanded.38.tsv.gz"
      }
    }
  }
}
```
Use this config file with `--panel-config` and `--panel-name my_custom_panel` during targeted analysis runs.
