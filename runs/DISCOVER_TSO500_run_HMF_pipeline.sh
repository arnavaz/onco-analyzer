module load python3/3.12.11

PANEL_NAME=tso500
OUTDIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/DISCOVER_DISTINCTION_Shamini/HMF/results
ACCOUNT=deperrotgroup

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab


## Tumour-Only

python3 /cluster/projects/pughlab/src/oncoanalyser_pughlab/Setup/make_oncoanalyser_samplesheet.py \
   --bam-dir /cluster/projects/mohccn/Sandbox/JEFF_TMP/DISCOVER_DISTINCTION_Shamini/HMF/bam_links \
   --out $OUTDIR

  # Single-sample run:
PANEL_NAME=tso500

SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab
SAMPLE=/cluster/projects/mohccn/Sandbox/JEFF_TMP/DISCOVER_DISTINCTION_Shamini/HMF/results/TSOCF151_tumor/samplesheet.csv
OUTDIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/DISCOVER_DISTINCTION_Shamini/HMF/results/TSOCF151_tumor/
mkdir -p "$OUTDIR/logs"
sample_name="TSOCF151_tumor"

sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A deperrotgroup \
  --job-name="onco_${sample_name}" \
  --output="${OUTDIR}/logs/onco_${sample_name}_%j.out" \
  --error="${OUTDIR}/logs/onco_${sample_name}_%j.err" \
  ${SOURCE_REPO}/run_oncoanalyser.sh \
  "$SAMPLE" \
  "$OUTDIR" \
  targeted \
  GRCh38_hmf \
  deperrotgroup \
  --panel-name tso500 \
  --realign

## Loop through multiple remaining patients:


SOURCE_REPO=/cluster/projects/pughlab/src/oncoanalyser_pughlab
SHEETS_DIR=/cluster/projects/mohccn/Sandbox/JEFF_TMP/DISCOVER_DISTINCTION_Shamini/HMF/results
RESULTS_BASE=/cluster/projects/mohccn/Sandbox/JEFF_TMP/DISCOVER_DISTINCTION_Shamini/HMF/results
ACCOUNT=deperrotgroup
PANEL_NAME=tso500


SKIP_SAMPLE=TSOCF151_tumor

for sample_sheet in "$SHEETS_DIR"/*/samplesheet.csv; do
  sample_name=$(basename "$(dirname "$sample_sheet")")

  if [[ "$sample_name" == "$SKIP_SAMPLE" ]]; then
    echo "Skipping already-submitted sample: $sample_name"
    continue
  fi

  outdir="$RESULTS_BASE/$sample_name"
  mkdir -p "$outdir/logs"

  echo "Submitting: $sample_name"
  sbatch --export=ALL,SOURCE_REPO=${SOURCE_REPO},NXF_OFFLINE=true -A "$ACCOUNT" \
    --job-name="onco_${sample_name}" \
    --output="${outdir}/logs/onco_${sample_name}_%j.out" \
    --error="${outdir}/logs/onco_${sample_name}_%j.err" \
    ${SOURCE_REPO}/run_oncoanalyser.sh \
    "$sample_sheet" \
    "$outdir" \
    targeted \
    GRCh38_hmf \
    "$ACCOUNT" \
    --panel-name "$PANEL_NAME" \
    --realign

done