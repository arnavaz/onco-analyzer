
## ----- Install: Nextflow and OncoAnalyzer on the cluster -----------------------

# Install Nextflow using H4H build node
h4h # alias: ssh jbruce@h4huhnlogin1.uhnresearch.ca
h4hb # alias: salloc --partition=build -c 1 -t 1:0:0 --mem 2G

module load java/21
curl -L https://github.com/nextflow-io/nextflow/releases/download/v25.10.5/nextflow \
    -o ~/.local/bin/nextflow
chmod +x ~/.local/bin/nextflow

nextflow -version

# clone oncoanalyzer using git on local but to pughlab bin
cd /Users/jbruce/Desktop/cluster/projects/pughlab/bin
git clone https://github.com/nf-core/oncoanalyser.git && cd oncoanalyser && git checkout master && pwd

# On interactive node
h4h
h4hI

cd /cluster/projects/pughlab/bin/oncoanalyser
find conf -type f | head -20

## ----- One-time SIF image cache setup (cluster build node with internet) ------
# Run on the build node after requesting temporary outbound internet from your admin
h4hb

module load java/21

# Install nf-core into an isolated conda env (avoids Python 3.7 system conflicts)
conda create -y -n nf-core python=3.11
conda activate nf-core
pip install nf-core pandas matplotlib
nf-core --version   # should show 4.x

module load singularity/4.3.2
module load java/21

export NXF_SINGULARITY_CACHEDIR=/mnt/users/jeff_singularity/singularity_cache
export SINGULARITY_CACHEDIR=/mnt/users/jeff_singularity/singularity_cache

# Downloads all SIF images directly into the shared cache and produces
# a self-contained pipeline bundle alongside the already-cloned pipeline.
# --container-cache-utilisation amend: skip images already present in the cache
nf-core pipelines download nf-core/oncoanalyser \
    --revision 2.3.0 \
    --container-system singularity \
    --container-cache-utilisation amend \
    --download-configuration yes \
    --outdir /mnt/users/jeff_singularity/oncoanalyser_bundle

# Verify cache
du -sh /mnt/users/jeff_singularity/singularity_cache

## ----- One-time reference data download (same internet session, build node) --
# Downloads HMF tool data bundle + genome FASTA/indices for GRCh38
# ~300-500 GB — use a himem node if possible for the decompression steps
# NOTE: run this AFTER the SIF cache is populated (uses the same singularity containers)

mkdir -p /mnt/users/jeff_singularity/oncoanalyser_refs

nextflow run /mnt/users/jeff_singularity/oncoanalyser_bundle/2_3_0/main.nf \
    --mode prepare_reference \
    --genome GRCh38_hmf \
    --ref_data_types wgs,wts,bwamem2_index,dict,dna_alignment,fai,fasta,gridss_index,hmftools,img,panel,rna_alignment,star_index \
    --validate_params false \
    --outdir /mnt/users/jeff_singularity/oncoanalyser_refs \
    -w /mnt/users/jeff_singularity/oncoanalyser_refs/work

# Verify — key dirs that should now exist:
# /mnt/users/jeff_singularity/oncoanalyser_refs/reference_data/<version>/hmf_data/
# /mnt/users/jeff_singularity/oncoanalyser_refs/reference_data/<version>/GRCh38_hmf/
ls /mnt/users/jeff_singularity/oncoanalyser_refs/reference_data/
du -sh /mnt/users/jeff_singularity/oncoanalyser_refs/

# Clean the work dir after successful reference download (saves space)
nextflow clean -f /mnt/users/jeff_singularity/oncoanalyser_refs/work

## ----- Setup: directories and configs on the cluster -----------------------
# move shared directories to their permanent locations
mv /cluster/projects/pughlab/jeff_singularity/ /cluster/projects/pughlab/references/oncoanalyser_resources

## Fixed broken links

OLD=/mnt/users/jeff_singularity
NEW=/cluster/projects/pughlab/references/oncoanalyser_resources

# Preview links that still point to old location
find "$NEW" -type l -lname "${OLD}/*" | head
find "$NEW" -type l -lname "${OLD}/*" | wc -l

# Rewrite absolute symlink targets from OLD -> NEW
find "$NEW" -type l -lname "${OLD}/*" -print0 | while IFS= read -r -d '' L; do
  T="$(readlink "$L")"
  NT="${T/#$OLD/$NEW}"
  if [ -e "$NT" ] || [ -L "$NT" ]; then
    ln -sfn "$NT" "$L"
  else
    echo "Missing target after rewrite: $L -> $NT"
  fi
done

# Verify
find "$NEW" -type l -lname "${OLD}/*" | wc -l
find "$NEW" -xtype l | head


SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab

## ----- Test run (from H4H login node) --------------------------------------
# Uses built-in test profile with small public datasets — validates Singularity + SLURM

module load java/21
module load singularity/4.3.2
export NXF_OFFLINE=true
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache

cd /cluster/projects/pughlab/bin/oncoanalyser
mkdir -p test_run/output logs

# Preflight: confirm reference resources exist at configured location
ls /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/

conda activate nf-core

nextflow run main.nf \
    -profile singularity,test \
    -c conf/h4h.config \
    -w /cluster/projects/pughlab/scratch/oncoanalyser_work \
    --validate_params false \
    --mode wgts \
    --genome GRCh38_hmf \
    --outdir test_run/output

## ----- Production run (via sbatch) -----------------------------------------
# Edit samplesheet_template.csv with your real samples first
mkdir -p /cluster/projects/pughlab/results/oncoanalyser/PATIENT1/logs

sbatch run_oncoanalyser.sh \
    /cluster/projects/pughlab/data/PATIENT1/samplesheet.csv \
    /cluster/projects/pughlab/results/oncoanalyser/PATIENT1 \
    wgts \
    GRCh38_hmf
# work dir will be: /cluster/projects/pughlab/results/oncoanalyser/PATIENT1/work

## ----- Resume an interrupted run -------------------------------------------
RESUME=true sbatch run_oncoanalyser.sh \
    /cluster/projects/pughlab/data/PATIENT1/samplesheet.csv \
    /cluster/projects/pughlab/results/oncoanalyser/PATIENT1 \
    wgts \
    GRCh38_hmf

## ----- Monitoring ----------------------------------------------------------
squeue -u jbruce                           # view running SLURM jobs
tail -f logs/nf_<JOBID>.out               # follow Nextflow head job log
ls /cluster/projects/pughlab/results/oncoanalyser/PATIENT1/pipeline_info/  # HTML reports

## Additional dependancies for runtime plotting

# On build node
module load java/21
module load singularity/4.3.2
conda activate nf-core
python -m pip install --upgrade pip
python -m pip install pandas matplotlib

## For deveopment on laptop but syncing with cluster

# 1. On cluster, allow pushing into the checked-out branch:
cd /cluster/projects/pughlab/src/oncoanalyser_pughlab
git config receive.denyCurrentBranch updateInstead
git config receive.advertisePushOptions true

# 2. On local, add cluster as a remote:
cd ~/bin/oncoanalyser_pughlab
git remote add h4h jbruce@h4huhndata1:/cluster/projects/pughlab/src/oncoanalyser_pughlab

# 3. Create a convenient alias in your local repo:
git config --global alias.pushall '!git push h4h main && git push origin main'

# 4. Run alias to push to both
git pushall
Ido Nofech-Mozes et al PMID: 36959212


## Download GRCh37 references for TSO500
# 1. Define output reference directory
REF_DIR=/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh37
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# 2. Download the fasta and indices
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/Homo_sapiens.GRCh37.GATK.illumina.fasta
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/samtools_index/1.16/Homo_sapiens.GRCh37.GATK.illumina.fasta.dict
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh37_hmf/24.0/bwa_index_image/0.7.17-r1188/Homo_sapiens.GRCh37.GATK.illumina.fasta.img

# 3. Download the database tarballs
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/hmftools/hmf_pipeline_resources.37_v2.3.0--2.tar.gz
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.37_v2.3.0--2.tar.gz

# 4. Unpack tarballs
tar -xzvf hmf_pipeline_resources.37_v2.3.0--2.tar.gz
tar -xzvf hmf_panel_resources.tso500.37_v2.3.0--2.tar.gz

## NOT RUN YET

# 5. Clean up tarballs (optional, to save space)
rm *.tar.gz

## Add hg38

# 1. Navigate to your reference directory
cd /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/

# 2. Download the GRCh38 panel package
wget https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/hmf_reference_data/panels/hmf_panel_resources.tso500.38_v2.3.0--2.tar.gz

# 3. Unpack and clean up (on cluster for improved speed)
cd /cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/

tar -xzvf hmf_panel_resources.tso500.38_v2.3.0--2.tar.gz
rm hmf_panel_resources.tso500.38_v2.3.0--2.tar.gz



## updating singularity images on build node

## Zhibin mounted this dir for me to use temporarily
# 1. Set up directories
CACHE_DIR="/mnt/users/jeff_singularity/singularity_cache"
mkdir -p "${CACHE_DIR}"
cd "${CACHE_DIR}"

export SINGULARITY_CACHEDIR="${CACHE_DIR}/tmp_cache"
mkdir -p "${SINGULARITY_CACHEDIR}"

echo "=== 1. Pulling HMF Docker Hub Containers ==="
singularity pull --name hartwigmedicalfoundation-purple-4.5-beta.3.img docker://hartwigmedicalfoundation/purple:4.5-beta.3
singularity pull --name hartwigmedicalfoundation-sage-5.1-beta.2.img docker://hartwigmedicalfoundation/sage:5.1-beta.2
singularity pull --name hartwigmedicalfoundation-esvee-2.1-beta.2.img docker://hartwigmedicalfoundation/esvee:2.1-beta.2
singularity pull --name hartwigmedicalfoundation-bam-tools-1.7-beta.2.img docker://hartwigmedicalfoundation/bam-tools:1.7-beta.2
singularity pull --name hartwigmedicalfoundation-isofox-2.1-beta.2.img docker://hartwigmedicalfoundation/isofox:2.1-beta.2
singularity pull --name hartwigmedicalfoundation-linx-2.4-beta.2.img docker://hartwigmedicalfoundation/linx:2.4-beta.2
singularity pull --name hartwigmedicalfoundation-pave-2.0-beta.2.img docker://hartwigmedicalfoundation/pave:2.0-beta.2
singularity pull --name hartwigmedicalfoundation-wisp-1.4-beta.2.img docker://hartwigmedicalfoundation/wisp:1.4-beta.2
singularity pull --name hartwigmedicalfoundation-neo-1.4-beta.1.img docker://hartwigmedicalfoundation/neo:1.4-beta.1

## HERE HERE HERE


echo "=== 2. Pulling Upgraded Galaxy/BioContainers for 3.1.0 ==="
# Core HMF tools using Galaxy releases
singularity pull --name depot.galaxyproject.org-singularity-hmftools-amber-4.3--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-amber:4.3--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-cobalt-3.0--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-cobalt:3.0--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-orange-5.0.1--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-orange:5.0.1--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-redux-2.0.5--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-redux:2.0.5--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-teal-1.4--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-teal:1.4--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-cider-1.2--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-cider:1.2--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-cuppa-2.5.1--py311r42hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-cuppa:2.5.1--py311r42hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-lilac-2.0--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-lilac:2.0--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-qsee-1.0--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-qsee:1.0--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-hmftools-fastq-tools-1.0--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-fastq-tools:1.0--hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-r-linxreport-1.2.0--r45hdfd78af_0.img https://depot.galaxyproject.org/singularity/r-linxreport:1.2.0--r45hdfd78af_0
singularity pull --name depot.galaxyproject.org-singularity-ubuntu-24.04.img https://depot.galaxyproject.org/singularity/ubuntu:24.04

# SAGE Galaxy fallback image
singularity pull --name depot.galaxyproject.org-singularity-hmftools-sage-5.0.2--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-sage:5.0.2--hdfd78af_0

echo "=== 3. Standard Preprocessing Utilities ==="
singularity pull --name depot.galaxyproject.org-singularity-fastp-0.23.4--hadf994f_2.img https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2
singularity pull --name depot.galaxyproject.org-singularity-gatk4-4.6.1.0--py310hdfd78af_0.img https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0

# =============================================================================
# NEW (added after 3.1.0 container audit, 2026-10-10) -- run these on the build node
# =============================================================================
# Audit: every container referenced in pipeline/modules/** was checked against the cache.
# Only ONE image is required-but-missing: MULTIQC (3.1 uses multiqc 1.35 + the oncoanalyser
# multiqc plugin; the older 3.0 mulled-v2-4dde5019... image is NOT compatible).
# Without it every WGTS run fails at the very last step (MULTIQC_REPORTING:MULTIQC).
echo "=== NEW A. REQUIRED: MultiQC (3.1.0) ==="
singularity pull --name depot.galaxyproject.org-singularity-mulled-v2-d8bcaa0f18f0068460deb4841052ef5429108a27-f0d5b44eec0233d9bc8357f24755076bceef041a-0.img https://depot.galaxyproject.org/singularity/mulled-v2-d8bcaa0f18f0068460deb4841052ef5429108a27:f0d5b44eec0233d9bc8357f24755076bceef041a-0

## don't have the other dir mounted right now - using home dir for this single pull
SRC="/cluster/home/jbruce/singularity_images"
DEST="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache"

# Move new .img files into the permanent cache without overwriting unchanged ones
mv -vn "${SRC}"/*.img "${DEST}/"

# NEW B. OPTIONAL fallbacks -- NOT needed while pipeline/conf/hmf_containers.config is active
# (it overrides these processes with the beta images pulled in section 1). Only pull these if
# you ever disable hmf_containers.config; then these Galaxy defaults from the modules are used.
#echo "=== NEW B. OPTIONAL: module-default Galaxy images (fallback only) ==="
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-bam-tools-1.6.1--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-bam-tools:1.6.1--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-esvee-2.0.1--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-esvee:2.0.1--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-isofox-2.0.1--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-isofox:2.0.1--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-linx-2.3.1--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-linx:2.3.1--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-neo-1.3--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-neo:1.3--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-pave-1.9--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-pave:1.9--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-purple-4.4--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-purple:4.4--hdfd78af_0
#singularity pull --name depot.galaxyproject.org-singularity-hmftools-wisp-1.3.1--hdfd78af_0.img https://depot.galaxyproject.org/singularity/hmftools-wisp:1.3.1--hdfd78af_0
# =============================================================================
# END NEW
# =============================================================================

# Clean temporary layer cache
rm -rf "${SINGULARITY_CACHEDIR}"
echo "All Oncoanalyser 3.1.0 images successfully pulled to ${CACHE_DIR}!"

## Sync to cluster

SRC="/cluster/projects/pughlab/jeff_singularity/singularity_cache"
DEST="/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache"

# Move new .img files into the permanent cache without overwriting unchanged ones
mv -vn "${SRC}"/*.img "${DEST}/"

# Ensure Galaxy alias symlinks exist in the shared cache
cd "${DEST}"
for f in depot.galaxyproject.org-singularity-*.img; do
    # NEW (2026-10-10) safety guard: some old-version entries are depot.* SYMLINKS pointing at a real
    # file stored under the stripped name. Plain `ln -sf` would overwrite that real file with a
    # circular link. Only alias real files, and never replace an existing name.
    alias_name="${f#depot.galaxyproject.org-singularity-}"
    [[ -f "$f" && ! -L "$f" && ! -e "$alias_name" ]] && ln -s "$f" "$alias_name"
done

# =============================================================================
# NEW (2026-10-10): docker.io- aliases for the HMF beta images
# =============================================================================
# pipeline/conf/hmf_containers.config uses docker.io/hartwigmedicalfoundation/<tool>:<ver>, so
# Nextflow looks for docker.io-hartwigmedicalfoundation-<tool>-<ver>.img. The pulls in section 1
# are saved WITHOUT the docker.io- prefix, so Nextflow tried (and failed) to download them.
# Idempotent. Already done by hand on 2026-10-09; rerun after pulling any new beta image.
cd "${DEST}"
for f in hartwigmedicalfoundation-*.img; do
    [[ -f "$f" ]] && ln -sf "$f" "docker.io-$f"
done
# =============================================================================
# END NEW
# =============================================================================


## Download resources from Paul

TARGET_DIR="/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0"
mkdir -p "${TARGET_DIR}"

# 1. Download the new 3.1.0 HMF tool resource bundle (12.27 GB)
gcloud storage cp -r "gs://hmf-310-transfer/3.1.0-1/hmftools" "${TARGET_DIR}/"
mv "${TARGET_DIR}/hmftools" "${TARGET_DIR}/hmf_pipeline_resources.38_v3.1.0--1"

# 2. Copy the release metadata
gcloud storage cp "gs://hmf-310-transfer/3.1.0-1/release.json" "${TARGET_DIR}/"

# 3. Symlink the shared FASTA and STAR indexes from 2.3.0 (saves 116 GB of disk space)
REF_230="/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0"
cd "${TARGET_DIR}"
for f in "${REF_230}"/GRCh38_masked_exclusions_alts_hlas.fasta* \
         "${REF_230}"/bwa-mem2_index-2.2.1 \
         "${REF_230}"/gridss_index-2.13.2 \
         "${REF_230}"/star_index-gencode_38-2.7.3a \
         "${REF_230}"/hmf_panel_resources.tso500.38_v2.3.0--2; do
    ln -sf "$f" ./
done

## Done

