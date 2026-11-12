#!/usr/bin/env python3
import sys
import os
import glob
from collections import defaultdict

if len(sys.argv) < 3:
    print("Usage: python3 merge_cna.py <cna_tmp_dir> <output_file>")
    sys.exit(1)

tmp_dir = sys.argv[1]
out_file = sys.argv[2]

tsv_files = glob.glob(os.path.join(tmp_dir, "*.cna.tsv"))
if not tsv_files:
    print(f"No CNA TSVs found in {tmp_dir}")
    sys.exit(0)

samples = []
# data[gene][sample] = value
data = defaultdict(dict)

for tsv in tsv_files:
    with open(tsv, 'r') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        sample_id = header[1] if len(header) >= 2 else os.path.basename(tsv).replace(".cna.tsv", "")
        samples.append(sample_id)
        
        for line in f:
            parts = line.rstrip('\r\n').split('\t')
            if len(parts) >= 2:
                gene_name = parts[0]
                val = parts[1]
                data[gene_name][sample_id] = val

# Write output in cBioPortal format
with open(out_file, 'w') as fout:
    # Header
    header = ["Hugo_Symbol"] + samples
    fout.write("\t".join(header) + "\n")
    
    # Rows
    for gene in sorted(data.keys()):
        row = [gene]
        for s in samples:
            row.append(data[gene].get(s, "0")) # Default to neutral (0)
        fout.write("\t".join(row) + "\n")

print(f"Merged {len(samples)} CNA profiles into {out_file}")
