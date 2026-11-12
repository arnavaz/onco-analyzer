import sys
import os
import glob
from collections import defaultdict

if len(sys.argv) < 3:
    print("Usage: python3 merge_expression.py <expr_tmp_dir> <output_file>")
    sys.exit(1)

tmp_dir = sys.argv[1]
out_file = sys.argv[2]

tsv_files = glob.glob(os.path.join(tmp_dir, "*.expr.tsv"))
if not tsv_files:
    print(f"No expression TSVs found in {tmp_dir}")
    sys.exit(0)

samples = []
# data[gene][sample] = value
data = defaultdict(dict)
# Map Hugo_Symbol -> Entrez_Gene_Id (using Ensembl GeneId as fallback since cBioPortal expects an ID or 0)
gene_ids = {}

for tsv in tsv_files:
    sample_id = os.path.basename(tsv).replace(".expr.tsv", "")
    samples.append(sample_id)
    with open(tsv, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                gene_name = parts[0]
                gene_id = parts[1]
                tpm = parts[2]
                
                # Keep the Ensembl ID or just use 0 if you want strict Entrez formatting
                gene_ids[gene_name] = gene_id
                data[gene_name][sample_id] = tpm

# Write output in cBioPortal format
with open(out_file, 'w') as fout:
    # Header
    header = ["Hugo_Symbol"] + samples
    fout.write("\t".join(header) + "\n")
    
    # Rows
    for gene in sorted(data.keys()):
        row = [gene]
        for s in samples:
            row.append(data[gene].get(s, "0.0"))
        fout.write("\t".join(row) + "\n")

print(f"Merged {len(samples)} expression profiles into {out_file}")
