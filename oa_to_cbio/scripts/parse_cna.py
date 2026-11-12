import sys
import csv
import os

if len(sys.argv) < 7:
    print("Usage: python3 parse_cna.py <cnv_in> <sample_id> <ploidy> <amp_ratio> <gain_ratio> <out_tsv>")
    sys.exit(1)

cnv_in = sys.argv[1]
sample_id = sys.argv[2]
try:
    ploidy = float(sys.argv[3])
except ValueError:
    ploidy = 2.0
try:
    amp_ratio = float(sys.argv[4])
except ValueError:
    amp_ratio = 3.0
try:
    gain_ratio = float(sys.argv[5])
except ValueError:
    gain_ratio = 2.0
out_tsv = sys.argv[6]

os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)

with open(cnv_in, 'r') as f_in, open(out_tsv, 'w') as f_out:
    reader = csv.DictReader(f_in, delimiter='\t')
    writer = csv.writer(f_out, delimiter='\t')
    
    writer.writerow(["Hugo_Symbol", sample_id])
    
    for row in reader:
        gene = row.get("gene", "").strip()
        if not gene:
            continue
        
        # Only process canonical transcript row if present (standard PURPLE format)
        is_canonical = row.get("isCanonical", "true").lower() == "true"
        if not is_canonical:
            continue
            
        try:
            min_cn = float(row.get("minCopyNumber", 2.0))
            max_cn = float(row.get("maxCopyNumber", 2.0))
            min_minor_cn = float(row.get("minMinorAlleleCopyNumber", 1.0))
            
            # Homozygous deletion: minCopyNumber < 0.5
            if min_cn < 0.5:
                val = -2
            # Heterozygous deletion / LOH: minMinorAlleleCopyNumber < 0.5
            elif min_minor_cn < 0.5:
                val = -1
            # High-level amplification: maxCopyNumber > amp_ratio * ploidy
            elif max_cn > (amp_ratio * ploidy):
                val = 2
            # Gain: maxCopyNumber > gain_ratio * ploidy
            elif max_cn > (gain_ratio * ploidy):
                val = 1
            else:
                val = 0
                
            writer.writerow([gene, val])
        except ValueError:
            continue
