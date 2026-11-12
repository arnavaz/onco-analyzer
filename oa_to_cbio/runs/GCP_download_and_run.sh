
## On mac
cd /Users/jbruce/bin/oncoanalyser_pughlab
bash oa_to_cbio/scripts/sync_from_gcp.sh --cohort MOHCCN-P37 --download-only

## On cluster
cd /cluster/projects/pughlab/src/oncoanalyser_pughlab
bash oa_to_cbio/scripts/sync_from_gcp.sh --cohort PANCAN-G --skip-download


## On mac
cd /Users/jbruce/bin/oncoanalyser_pughlab

## On cluster
cd /cluster/projects/pughlab/src/oncoanalyser_pughlab
bash oa_to_cbio/scripts/sync_from_gcp.sh --cohort MOHCCN-P7 --skip-download

## Download the rest - on mac

cd /Users/jbruce/bin/oncoanalyser_pughlab
bash oa_to_cbio/scripts/sync_from_gcp.sh --download-only

## On cluster
cd /cluster/projects/pughlab/src/oncoanalyser_pughlab
bash oa_to_cbio/scripts/sync_from_gcp.sh --skip-download

cd /cluster/projects/pughlab/src/oncoanalyser_pughlab
bash oa_to_cbio/scripts/resummarize_clinical.sh --all

salloc -c 32 -p pughlab -t 10:0:0 --mem 250G

module load python3/3.10.9
 
cd /cluster/projects/pughlab/src/oncoanalyser_pughlab
## Note: Batches not annotated yet; applies cross-cohort TMM library scaling while skipping ComBat
python3 scripts/build_rna_database.py \
    --results-dir /cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download/* \
    --cbio-dir /cluster/projects/mohccn/oncoanalyser_cbioportal/cbio_import/* \
    --output-db /cluster/projects/mohccn/oncoanalyser_cbioportal/hmf_rnaseq_cohort.duckdb \
    --batch-column BATCH_ID \
    --threads 28

## Too slow on mount, copy to local

scp jbruce@h4huhndata1:/cluster/projects/mohccn/oncoanalyser_cbioportal/hmf_rnaseq_cohort.duckdb ~/Downloads

## Local visualization / web server
# cd ~/bin/oncoanalyser_pughlab
# python3 scripts/serve_rna_studio.py --db ~/Downloads/hmf_rnaseq_cohort.duckdb --port 8080



## outputs go to
/cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download/
/cluster/projects/mohccn/oncoanalyser_cbioportal/cbio_import/

mkdir -p ~/Downloads/CUP_orange_reports

for id in \
  MOHCCN-P16-00938 \
  MOHCCN-P16-01305 \
  MOHCCN-P16-04936 \
  MOHCCN-P2-02927 \
  PANCAN-D-0633 \
  MOHCCN-P2-02796 \
  MOHCCN-P2-03720 \
  MOHCCN-P2-05403 \
  MOHCCN-P2-05656 \
  MOHCCN-P17-01788 \
  MOHCCN-P2-02795 \
  MOHCCN-P2-03216 \
  MOHCCN-P2-03422 \
  MOHCCN-P36-02294 \
  MOHCCN-P2-03570 \
  MOHCCN-P1-00022 \
  MOHCCN-P16-02356 \
  MOHCCN-P2-00962 \
  MOHCCN-P2-01387 \
  MOHCCN-P2-02564 \
  MOHCCN-P2-03183 \
  MOHCCN-P8-00500 \
  PANCAN-B-0046 \
  PANCAN-B-0059 \
  MOHCCN-P2-03208 \
  MOHCCN-P2-03595 \
  MOHCCN-P2-03650 \
  MOHCCN-P2-00526 \
  MOHCCN-P2-01396 \
  MOHCCN-P2-01413 \
  MOHCCN-P2-01429 
do
    cp /cluster/projects/mohccn/oncoanalyser_cbioportal/gcp_download/*/"$id"/orange/"$id"-*.orange.pdf \
       ~/Downloads/CUP_orange_reports/
done

