# Running Oncoanalyser on HPC4Health (H4H)

Comprehensive guide, operational workflow, and cohort analysis documentation for Pugh Lab users.

Oncoanalyser is an `nf-core` Nextflow pipeline for cancer Whole Genome / Transcriptome Sequencing (WGS/WTS/WGTS) and targeted panel analysis based on the Hartwig Medical Foundation (HMF) software suite.

- **Upstream Repository**: [https://github.com/nf-core/oncoanalyser](https://github.com/nf-core/oncoanalyser)
- **H4H Cluster Source**: `/cluster/projects/pughlab/src/oncoanalyser_pughlab`
- **SLURM Partitions**: `pughlab`, `deperrotgroup`

---

## ⚡ Quickstart: Onboard a New Cohort in 1 Step

Whenever a new cohort arrives with FASTQ files and release map TSVs (even with multiple release batches, schema variations, or redundant archive folders), use **`Setup/prepare_cohort_run.py`** to ingest everything, deduplicate paths, verify metadata concordance, and generate a turnkey SLURM runner:

```bash
module load python3/3.12.11

python3 Setup/prepare_cohort_run.py \
  --cohort-name <COHORT> \
  --map-files /path/to/metadata/*.fastqs.map.tsv \
  --fastq-list /path/to/all_fastq_paths.txt \
  --exclude-duplicate-pattern "OICR_data/fastqs" \
  --account pughlab \
  --partition pughlab \
  --outdir /cluster/projects/<PROJECT>/TEMP/<COHORT>/Oncoanalyser \
  --target-link-dir /cluster/projects/<PROJECT>/TEMP/<COHORT>/all_fastq_links
```

### What this automatically handles:
- **Header Harmonization**: Reconciles changing column names across releases (`sample_id` vs `donor_id`, `case_id`, `group_id`).
- **FASTQ Deduplication**: Excludes duplicate files found across archive or transfer paths.
- **100% Concordance Assertion**: Verifies that every single FASTQ in the mapping file exists before submission.
- **Fail-Safe SLURM Runner (`runs/run_<COHORT>.sh`)**:
  - Builds symlink directory without collisions.
  - Automatically runs `Setup/make_oncoanalyser_samplesheet.py`.
  - Caps WGTS head jobs at safe walltimes (`--time=5-00:00:00`) to prevent SLURM rejection on `himem` / `all` partitions.
  - Enforces two-stage gating: plasma MRD/WISP jobs wait on primary tumor outputs (`--dependency=afterok:${wgts_job_id}`) or launch in parallel if tumor results already exist.
  - Integrates with Cluster Job Tracker (`track_batch_init` and `track_job`).

> [!TIP]
> **Working with AI Coding Assistants?**  
> This repository includes the **`cohort-run-generator`** agent skill (`.agents/skills/cohort-run-generator/SKILL.md`). You can simply provide your raw FASTQ list and map files to an assistant and instruct:  
> *"Use the cohort-run-generator skill to set up `<COHORT>`."*

---

## Table of Contents

- [⚡ Quickstart: Onboard a New Cohort in 1 Step](#-quickstart-onboard-a-new-cohort-in-1-step)
1. [Shared Resources & Infrastructure Profile](#shared-resources--infrastructure-profile)
2. [One-Time User Setup](#one-time-user-setup)
3. [Samplesheet Generation & Cohort Intake](#samplesheet-generation--cohort-intake)
   - [Automated Cohort Intake & Runner Generation (`prepare_cohort_run.py`)](#automated-cohort-intake--runner-generation-prepare_cohort_runpy)
   - [Minimum 4-Column Format](#minimum-required-columns-4-columns)
   - [Recommended Format with Longitudinal Plasma ctDNA (WISP)](#recommended-format-includes-patient-id--plasma-ctdna-wisp-tracking)
   - [RNA-Seq Only Samplesheets](#rna-seq-only-samplesheets)
4. [Production Pipeline Workflows](#production-pipeline-workflows)
   - [1. Whole Genome / Transcriptome Sequencing (WGTS)](#1-whole-genome--transcriptome-sequencing-wgts)
   - [2. Longitudinal Plasma ctDNA & Minimal Residual Disease (WISP MRD)](#2-longitudinal-plasma-ctdna--minimal-residual-disease-wisp-mrd)
   - [3. Targeted Custom Panels & Off-Target ichorCNA](#3-targeted-custom-panels--off-target-ichorcna)
   - [4. RNA-Seq Only & Multi-Omics Analysis](#4-rna-seq-only--multi-omics-analysis)
5. [Automated Resuming & Fault Tolerance](#automated-resuming--fault-tolerance)
6. [Cohort-Level Compilation & Reporting](#cohort-level-compilation--reporting)
   - [Master Cohort Summary (`compile_wgts_summary.py`)](#master-cohort-summary-compile_wgts_summarypy)
   - [Interactive Purity Curation Console (`curate_purity_fits.py`)](#interactive-purity-curation-console-curate_purity_fitspy)
   - [MMR RNA-Seq Dashboard (`compile_mmr_rna_dashboard.py`)](#mmr-rna-seq-dashboard-compile_mmr_rna_dashboardpy)
   - [High-Resolution BAF & Epigenomic Plots (`plot_gene_baf.py`)](#high-resolution-baf--epigenomic-plots-plot_gene_bafpy)
7. [cBioPortal Data Staging (`oa_to_cbio`)](#cbioportal-data-staging-oa_to_cbio)
8. [Git Workflow & Cluster Synchronization](#git-workflow--cluster-synchronization)
9. [Repository File Inventory](#repository-file-inventory)

---

## Shared Resources & Infrastructure Profile

On H4H, compute nodes and file systems mounting `/cluster/projects` have **no outbound internet access** (`NXF_OFFLINE=true`). All containers and reference assets must be pre-staged.

- **Singularity Image Cache**:  
  `/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache`  
  *(All Apptainer/Singularity container images, including HMF tools, Sage, Amber, Cobalt, Purple, Linx, Esvee, Pave, Cuppa, and Isofox).*
- **HMF Pipeline References (v2.3.0)**:  
  `/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/hmf_pipeline_resources.38_v2.3.0--2`
- **HMF Pipeline References (v3.1.0-beta)**:  
  `/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0/hmf_pipeline_resources.38_v3.1.0--1`
- **Shared Pipeline Root**:  
  `/cluster/projects/pughlab/src/oncoanalyser_pughlab`
- **Shared Cluster Config**:  
  `/cluster/projects/pughlab/src/oncoanalyser_pughlab/conf/h4h.config`

> [!NOTE]
> **Downloading v3.1 HMF Resources:**  
> When new references are released, download them externally and transfer to H4H:
> ```bash
> gcloud storage cp gs://hmf-public/HMFtools-Resources/pipeline/oncoanalyser/3.1-beta/pipeline-resources.hmftools.3_1_0-2.tar.gz /path/to/reference_data/
> ```

### SLURM Submission Policy: Always Use Inline CLI Flags

> [!IMPORTANT]
> **Always submit SLURM jobs using explicit inline CLI flags rather than relying on `#SBATCH` script headers alone:**
> ```bash
> sbatch -A pughlab -p pughlab [options] <script.sh>
> ```
> **Why this matters on H4H:**  
> When invoking `sbatch <script.sh>` without inline flags, SLURM validates partition access against your user ID's **default account** (which is `deperrotgroup` for many users) *before* or independently of reading `#SBATCH --account` directives inside the script. If the header specifies `#SBATCH --partition=pughlab`, SLURM evaluates `deperrotgroup` + `pughlab` (an unauthorized combination) and immediately fails with:
> ```
> sbatch: error: Batch job submission failed: Invalid account or account/partition combination specified
> ```
> Passing `-A <account> -p <partition>` directly on the command line establishes the account context upfront, avoiding this error and bypassing any cluster mount NFS attribute cache delays.

---

## One-Time User Setup

Create an `nf-core` environment on the cluster login node:

```bash
module load java/21
module load singularity/4.3.2

conda create -y -n nf-core python=3.11
conda activate nf-core
pip install nf-core pandas matplotlib
```

Ensure Nextflow is installed in your user path (`~/.local/bin`):

```bash
mkdir -p ~/.local/bin
curl -L https://github.com/nextflow-io/nextflow/releases/download/v25.10.5/nextflow -o ~/.local/bin/nextflow
chmod +x ~/.local/bin/nextflow
export PATH="$HOME/.local/bin:$PATH"
nextflow -version
```

Nextflow plugins must be installed while on the cluster login node (pinned to `nf-schema@2.5.1`):
```bash
nextflow plugin install nf-schema@2.5.1
```

---

## Samplesheet Generation & Cohort Intake

### Automated Cohort Intake & Runner Generation (`prepare_cohort_run.py`)

When starting a new cohort with multiple fastq releases, inconsistent mapping TSV headers, or redundant archive fastqs, use `Setup/prepare_cohort_run.py` to automate end-to-end intake and generate a turnkey SLURM runner:

```bash
module load python3/3.12.11

python3 Setup/prepare_cohort_run.py \
  --cohort-name NETMRD \
  --map-files metadata/NETMRD.release.*.fastqs.map.tsv \
  --fastq-list /cluster/projects/.../NETMRD_fastq_list.txt \
  --exclude-duplicate-pattern "OICR_data/fastqs" \
  --account pughlab \
  --partition pughlab \
  --outdir /cluster/projects/rskimlab/TEMP/NETMRD/Oncoanalyser \
  --target-link-dir /cluster/projects/rskimlab/TEMP/NETMRD/all_fastq_links
```

**What this produces automatically:**
1. **Unified Map TSV** (`runs/metadata/<COHORT>.fastqs.map.tsv`): Merges multi-batch mapping manifests, harmonizes field schemas, and resolves clinical donor vs sample IDs.
2. **Deduplicated FASTQ List** (`runs/metadata/<COHORT>_fastq_list.dedup.txt`): Filters out redundant fastqs matching `--exclude-duplicate-pattern` and verifies 100% concordance against the mapping file.
3. **Turnkey SLURM Runner** (`runs/run_<COHORT>.sh`): Generates a production script that:
   - Sets up symlinks without duplicates.
   - Generates per-patient samplesheets via `make_oncoanalyser_samplesheet.py`.
   - Submits Step 1 (WGTS) with safe `--time=5-00:00:00` walltime caps.
   - Enforces two-stage gating: queues Step 2 (MRD / WISP) with `--dependency=afterok:<job_id>` or launches immediately if PURPLE tumor outputs already exist.
   - Integrates with Cluster Job Tracker (`track_batch_init` and `track_job`).

---

### Manual / Direct Samplesheet Generation

Generate standardized per-patient samplesheets using `Setup/make_oncoanalyser_samplesheet.py`:

```bash
module load python3/3.12.11

python3 Setup/make_oncoanalyser_samplesheet.py \
  --fastq-dir /path/to/fastqs \
  --mapping-tsv COHORT_injection_sheet.tsv \
  --out /path/to/output/sample_sheets
```

`--mapping-tsv` accepts standard OICR/Pugh Lab injection manifests or minimal TSV mapping files.

### Minimum Required Columns (4 Columns)

| Column Header | Accepted Aliases | Purpose / Accepted Values |
| :--- | :--- | :--- |
| `fastq file name` | `files`, `filename` | FASTQ file name(s) (R1/R2 filenames, separated by `;` or listed across rows) |
| `sample id` | `sample_id`, `sample` | Unique sample identifier (e.g. `PAT01-T01`, `PAT01-N01`) |
| `tissue type` | `tissue_type`, `tissue type (eg. primary tumor, reference)` | `r` (Normal/Reference)<br>`primary tumor`, `tumor`, `metastatic` (Tumor) |
| `library type` | `library_type`, `library type (wg:whole genome, wt:whole transcriptome, pg:plasma, ex:exome)` | `wg` (WGS DNA), `wt` (RNA), `ex` (Exome), `pg` (Plasma DNA) |

### Recommended Format (Includes Patient ID & Plasma ctDNA WISP Tracking)

Providing `donor_id` and `tissue origin` enables automatic multi-sample grouping and longitudinal plasma sample tagging for WISP:

```tsv
fastq file name	sample id	tissue type	library type	donor_id	tissue origin
PAT01_T01_L001_R1_001.fastq.gz;PAT01_T01_L001_R2_001.fastq.gz	PAT01-T01	primary tumor	wg	PAT01	primary
PAT01_N01_L001_R1_001.fastq.gz;PAT01_N01_L001_R2_001.fastq.gz	PAT01-N01	r	wg	PAT01	blood
PAT01_PL01_L001_R1_001.fastq.gz;PAT01_PL01_L001_R2_001.fastq.gz	PAT01-PL01	plasma	wg	PAT01	Pl
PAT01_PL02_L001_R1_001.fastq.gz;PAT01_PL02_L001_R2_001.fastq.gz	PAT01-PL02	plasma	wg	PAT01	Pl
```

- **`donor_id`** (or `patient_id`, `subject_id`): Group normal, primary tumor, and longitudinal plasma samples under one patient.
- **`tissue origin`**: Setting this to **`Pl`** or **`plasma`** tags the sample as a longitudinal plasma sample and formats the `info` column as `library_id:<lib>;lane:<lane>;longitudinal_sample` for WISP ctDNA purity estimation.

### RNA-Seq Only Samplesheets

To generate sheets for RNA-only studies (e.g. MMR cohort):

```bash
python3 Setup/make_oncoanalyser_samplesheet.py \
  --fastq-dir /path/to/rna_fastqs \
  --mapping-tsv MMR_RNAseq_metadata.tsv \
  --out sample_sheets_rna \
  --default-sequence-type rna \
  --default-sample-type tumor
```

---

## Production Pipeline Workflows

All wrapper scripts read `conf/h4h.config`, enforce `NXF_OFFLINE=true`, and manage SLURM job submission.

### 1. Whole Genome / Transcriptome Sequencing (WGTS)

#### Submit a Whole Cohort (one head job per patient):
```bash
/cluster/projects/pughlab/src/oncoanalyser_pughlab/submit_all_patients.sh \
  /path/to/sample_sheets \
  /path/to/results_base \
  wgts \
  GRCh38_hmf \
  pughlab
```

#### Submit a Single Patient:
```bash
export NXF_OFFLINE=true
SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab

sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A pughlab \
  /cluster/projects/pughlab/src/oncoanalyser_pughlab/run_oncoanalyser.sh \
  /path/to/patient/samplesheet.csv \
  /path/to/patient/output \
  wgts \
  GRCh38_hmf \
  pughlab
```

---

### 2. Longitudinal Plasma ctDNA & Minimal Residual Disease (WISP MRD)

The MRD workflow operates in a two-tier execution pattern:
1. **Tier 1 (Tumor/Normal WGTS)**: Run the primary tumor and matched normal through `wgts` to generate baseline somatic SNVs (`sage/somatic/`), CNVs, and purity/ploidy models (`purple/`).
2. **Tier 2 (Plasma WISP)**: Process longitudinal plasma samples (`T0`, `T1`, `T2`, ...) in parallel using `--mode targeted` / WISP, which seeds somatic variants from Tier 1 into the plasma BAM to quantify ctDNA fraction across SNV, LOH, and copy number modalities.

Launch scripts are located in `runs/`:
- `runs/run_TRACEOC_MRD.sh` (Ovarian cancer MRD cohort)
- `runs/run_BFIERCE_MRD.sh` (Breast cancer MRD cohort)
- `runs/run_summary_HCC_MRD.sh` (Hepatocellular carcinoma cohort)

Key features:
- **Parallel Dispatching**: Evaluates all plasma timepoints concurrently.
- **Smart Skip**: Checks for completed `.wisp.summary.tsv` files and skips redundant executions.
- **AMBER BAF & LOH Integration**: Incorporates B-allele frequency tracking from Amber for sensitive copy-neutral loss of heterozygosity detection.

---

### 3. Targeted Custom Panels & Off-Target ichorCNA

The pipeline supports custom targeted panel workflows and off-target read profiling:

#### 1) Create Custom Panel Resources (`panel_resource_creation`):
```bash
sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A pughlab \
  /cluster/projects/pughlab/src/oncoanalyser_pughlab/run_oncoanalyser.sh \
  /path/to/samplesheet.panel_resource_creation.csv \
  /path/to/output/panel_resource_creation \
  panel_resource_creation \
  GRCh38_hmf \
  pughlab \
  --driver-gene-panel /path/to/DriverGenePanel.38.tsv \
  --target-regions-bed /path/to/target_regions_definition.38.bed.gz
```

#### 2) Run Targeted Analysis with Off-Target ichorCNA:
Targeted panel runs can estimate copy number alterations and tumor fractions from off-target reads:
```bash
sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A pughlab \
  /cluster/projects/pughlab/src/oncoanalyser_pughlab/run_oncoanalyser.sh \
  /path/to/patient/samplesheet.csv \
  /path/to/patient/output \
  targeted \
  GRCh38_hmf \
  pughlab \
  --panel-name my_custom_panel \
  --panel-config /path/to/panel_data.config \
  --force-panel \
  --ichor-exclusion-bed /path/to/custom_exclusion.bed \
  --ichor-pon /path/to/custom_panelOfNormals.rds
```

#### 3) sWGS ichorCNA Auto-Discovery & Caching:
`run_oncoanalyser.sh` automatically checks for pre-computed sWGS ichorCNA results (`*params.txt` under `ichorCNA/`) across 5 distinct modalities:
- `panel_tumor_normal`
- `panel_tumor_only`
- `panel_plasma_only`
- `panel_normal_only`
- `swgs_tumor_only`

If found, results are reused directly without re-running ichorCNA.

---

### 4. RNA-Seq Only & Multi-Omics Analysis

Run Oncoanalyser exclusively on RNA samples for gene expression quantification, fusion detection, and immune repertoire exploration:

```bash
bash runs/run_MMR_RNA.sh
```

Key subroutines:
- **Isofox**: Transcript/gene-level expression (TPM) and splicing analysis.
- **Isofox vs Hartwig Reference**: Percentile ranking against thousands of cancer-type matched tumors from the Hartwig Medical Foundation cohort (`gene_distribution.csv.gz`).
- **TRUST4 / Immune Repertoire**: TCR and BCR CDR3 clonotype extraction and diversity profiling.
- **Neoantigen Binding Model**: HLA class I binding prediction for fusion peptides (`scripts/hla_binding.py`, `scripts/build_tumor_antigen_reference.py`).

---

## Automated Resuming & Fault Tolerance

`run_oncoanalyser.sh` includes automatic failure detection and backoff retry with Nextflow `-resume`.

Default retry settings:
```bash
AUTO_RESUME_ON_FAILURE=true
MAX_RESUME_ATTEMPTS=3
RESUME_BACKOFF_SECONDS=300
```

To override for a run:
```bash
AUTO_RESUME_ON_FAILURE=true MAX_RESUME_ATTEMPTS=5 RESUME_BACKOFF_SECONDS=180 \
sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A pughlab ...
```

Upon successful pipeline completion:
- Output BAMs are automatically converted to CRAM using the reference FASTA to reduce disk space footprint.
- Intermediate BAMs are cleaned up (`CONVERT_FINAL_BAM_TO_CRAM=true`).

---

## Cohort-Level Compilation & Reporting

### Master Cohort Summary (`compile_wgts_summary.py`)

Compiles multi-modal genomic and longitudinal ctDNA MRD data across an entire cohort into publications figures and a standalone interactive HTML dashboard.

```bash
module load python3/3.10.9

python3 scripts/compile_wgts_summary.py \
  --results-dir /path/to/cohort/results \
  --output-dir /path/to/summary_output \
  --cohort-name "MY_COHORT" \
  --timepoints-tsv /path/to/metadata_or_timepoints.tsv \
  --genes-of-interest TP53,KRAS,PIK3CA,ARID1A,BAP1,NF2,SETD2,CDKN2A \
  --purity-overrides /path/to/purity_overrides.tsv
```

#### Key Capabilities:
- **Directory Agnostic**: Automatically discovers samples across nested patient hierarchies (`results/PATIENT/tumor`, `results/PATIENT/plasma/T0`) as well as completely flat per-sample run directories (`results/PATIENT`, `results/PATIENT-T0`, `results/PATIENT-T1`).
- **Flexible Timepoint / Metadata Ingestion**: Accepts 2-column or 3-column TSV/CSVs with or without headers (mapping `patient_id`, `sample_id`, `timepoint`, `date`, `days`).
- **Longitudinal ctDNA Trajectories**: Generates per-patient log-scale trajectory plots with shaded assay Limit of Detection (LOD) regions.
- **Multi-Modal WISP Decomposition**: Visualizes SNV, LOH, and Copy Number ctDNA purity estimates.
- **Sample Identity & Swap Detection**: Evaluates germline BAF concordance from Amber across all sample pairs (`CONFIRMED_SAME_DONOR`, `MISMATCH_SUSPECTED_SWAP`, `SINGLE_SAMPLE_NO_COMPARATOR`).
- **Embedded Orange PDF Reports**: Automatically discovers companion native HMF Orange clinical reports, stages them into `orange_reports/`, and links them directly inside the interactive dashboard.

#### Outputs Generated:
| Output File | Description |
| :--- | :--- |
| `cohort_wgts_summary.tsv` | Master patient table (clinical KPIs, purity, ploidy, TMB, MSI, drivers, ctDNA status) |
| `cohort_wgts_samples.tsv` | Tidy sample-level metadata (Tumor, Normal, Plasma) |
| `cohort_wgts_mrd_summary.tsv` | Comprehensive longitudinal plasma metrics (SNV/LOH/CNV/WISP purity, LODs, p-values) |
| `cohort_wgts_patient_matrix.tsv` | Patient $\times$ Timepoint longitudinal matrix |
| `cohort_wgts_oncoprint.tsv` | Gene $\times$ Sample alteration matrix for Oncoprints |
| `cohort_fingerprint_summary.tsv` | Intra-patient sample concordance and identity verification calls |
| `cohort_fingerprint_matrix.tsv` | All-vs-all pairwise Jaccard distance matrix |
| `cohort_wgts_summary.png` / `.pdf` | Multi-panel publication figure (purity/ploidy, TMB, depth, Oncoprint) |
| `cohort_wgts_mrd_multimodal.png` / `.pdf` | Longitudinal ctDNA trajectory overview across all patients |
| `cohort_wgts_dashboard.html` / `index.html` | Standalone interactive dashboard with dynamic search, filtering, and embedded plots |
| `patient_plots/` | Per-patient publication-quality log-scale MRD trajectory plots |
| `plasma_plots/` | Per-plasma multi-modal WISP decomposition plots |

---

### Interactive Purity Curation Console (`curate_purity_fits.py`)

PURPLE occasionally faces ambiguity between diploid and tetraploid whole-genome duplication solutions, especially in subclonal or low-purity tumors.

Launch the purity curation console:
```bash
python3 scripts/curate_purity_fits.py \
  --results-dir /path/to/cohort/results \
  --output-html purity_curation_console.html \
  --embed-images
```

Features:
- Visualizes BAF vs Copy Number fit contours and alternative local minima.
- Allows user curation of purity and ploidy with interactive notes.
- Exports `purity_overrides.tsv`, which is directly ingested by `compile_wgts_summary.py --purity-overrides`.

---

### MMR RNA-Seq Dashboard (`compile_mmr_rna_dashboard.py`)

Interactive clinical dashboard for cohort RNA-Seq data:
```bash
python3 scripts/compile_mmr_rna_dashboard.py \
  --results-dir /path/to/rna_results \
  --metadata-tsv MMR_RNAseq_metadata.tsv \
  --output-html cohort_mmr_rna_dashboard.html
```

Features:
- **Interactive Volcano Plot**: Differential expression with dynamic fold-change and FDR cutoffs.
- **Hartwig Reference Percentiles**: Cancer-type matched percentile rankings for every gene.
- **TCR / BCR Repertoire**: Top CDR3 clonotypes, V/J segment usage, and clonotype diversity indices.
- **Fusion-TCR & Antigen Explorer**: Predicted neoantigen peptides with predicted HLA class I binding affinities.

---

### High-Resolution BAF & Epigenomic Plots (`plot_gene_baf.py`)

Generates multi-track focal gene plots:
```bash
python3 scripts/plot_gene_baf.py \
  --amber-baf /path/to/sample.amber.baf.tsv.gz \
  --purple-segments /path/to/sample.purple.segment.tsv \
  --gene BRCA2 \
  --output BRCA2_baf.png \
  --output-html BRCA2_baf.html
```

Features:
- Track A: Somatic mutation ploidy and allele-specific copy number.
- Track B: Phased BAF with segmented minor allele fits and rolling MAD confidence bands.
- Track C: DNA methylation beta-values from EM-Seq (`MethylDackel`).
- Macro-Track: Whole-chromosome copy number context.

---

## cBioPortal Data Staging (`oa_to_cbio`)

The `oa_to_cbio/` sub-package converts Oncoanalyser outputs into standard cBioPortal staging formats for upload to cBioPortal instances.

```bash
cd oa_to_cbio
./oa_to_cbio.sh \
  -i /path/to/oncoanalyser_output \
  -o /path/to/cbio_study_dir \
  -A pughlab \
  -s my_study_id \
  -n "My Study Display Name"
```

Features:
- Dispatches parallel SLURM jobs to format mutations (MAF), CNVs (discrete and continuous), structural variants, and fusions.
- Flattens mutational signatures into clinical attributes.
- Annotates ecDNA amplicon-driven gene amplifications.
- Synchronizes cohorts from Google Cloud Storage via `oa_to_cbio/scripts/sync_from_gcp.sh`.

---

## Git Workflow & Cluster Synchronization

To maintain synchronization between GitHub and the offline-restricted H4H cluster:

- **Remotes Configured**:
  - `origin`: GitHub (`https://github.com/pughlab/oncoanalyser_pughlab.git`)
  - `h4h`: H4H cluster repository (`jbruce@h4huhndata1:/cluster/projects/pughlab/src/oncoanalyser_pughlab`)
- **Pushing Changes**:
  Use `git pushall` to push checked-out branches simultaneously to both GitHub and the cluster:
  ```bash
  git pushall
  ```
  *(The cluster repository has `receive.denyCurrentBranch = updateInstead` configured, allowing direct push updates to the working tree).*

> [!CAUTION]
> **DO NOT** manually copy modified code files into mounted cluster filesystems (`/Users/jbruce/Desktop/cluster/...`). Doing so creates uncommitted working-tree modifications on the cluster that block subsequent `git push` operations. Always deploy via `git pushall`.

---

## Repository File Inventory

```text
├── conf/
│   └── h4h.config                   # SLURM resource, container, and reference profile for H4H
├── Setup/
│   ├── make_oncoanalyser_samplesheet.py  # Automated samplesheet builder from FASTQs + metadata
│   ├── samplesheet_template.csv     # Example CSV samplesheet template
│   └── installation_commands.sh     # Reference environment bootstrap commands
├── submit_all_patients.sh           # Cohort SLURM head job submission orchestrator
├── run_oncoanalyser.sh              # Single-patient runner with auto-resume & CRAM conversion
├── scripts/
│   ├── compile_wgts_summary.py      # Master multi-omics & longitudinal MRD dashboard compiler
│   ├── curate_purity_fits.py        # Web-based PURPLE purity/ploidy curation console
│   ├── compile_mmr_rna_dashboard.py # Interactive RNA-Seq & TCR repertoire dashboard
│   ├── plot_gene_baf.py             # Multi-track focal gene BAF & methylation plotter
│   ├── run_swgs_ichorcna.sh         # sWGS copy number & tumor fraction pipeline
│   ├── benchmark_ichor_bin_sizes.py # Bin size optimization benchmarking for ichorCNA
│   ├── generate_ichorcna_panel_resources.py # Exclusion BED & PoN builder for panel ichorCNA
│   ├── build_tumor_antigen_reference.py     # Neoantigen database curator
│   ├── hla_binding.py               # HLA class I binding prediction engine
│   └── analysis/
│       ├── compare_fingerprints.sh  # SAGE germline concordance comparator
│       └── summarize_runtime.py     # SLURM wall time & task runtime profiler
├── runs/                            # Cohort run wrappers (TRACEOC, BFIERCE, HCC, MMR)
│   ├── run_TRACEOC_MRD.sh
│   ├── run_BFIERCE_MRD.sh
│   ├── run_summary_HCC_MRD.sh
│   └── run_MMR_RNA.sh
└── oa_to_cbio/                      # cBioPortal export pipeline & GCP sync tools
    ├── oa_to_cbio.sh
    └── scripts/sync_from_gcp.sh
```
