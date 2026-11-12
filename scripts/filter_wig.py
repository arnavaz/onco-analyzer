#!/usr/bin/env python3
import sys
import os
from collections import defaultdict

def filter_wig_file(wig_path, exclusion_bed_path, out_wig_path):
    # Load exclusion BED
    exclusions = defaultdict(list)
    with open(exclusion_bed_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                exclusions[chrom].append((start, end))
                
    # Sort exclusions
    for chrom in exclusions:
        exclusions[chrom].sort(key=lambda x: x[0])
        
    def overlaps_exclusions(chrom, bin_start, bin_end):
        excl_list = exclusions.get(chrom, [])
        if not excl_list:
            return False
        import bisect
        idx = bisect.bisect_left(excl_list, (bin_start, 0))
        for i in [idx - 1, idx, idx + 1]:
            if 0 <= i < len(excl_list):
                estart, eend = excl_list[i]
                if max(bin_start, estart) < min(bin_end, eend):
                    return True
        return False

    with open(wig_path, 'r') as in_f, open(out_wig_path, 'w') as out_f:
        current_chrom = None
        current_start = None
        current_step = None
        current_span = None
        bin_idx = 0
        
        for line in in_f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith("track") or line_str.startswith("variableStep"):
                out_f.write(line)
            elif line_str.startswith("fixedStep"):
                out_f.write(line)
                parts = line_str.split()
                meta = {}
                for p in parts[1:]:
                    k, v = p.split("=")
                    meta[k] = v
                current_chrom = meta.get("chrom")
                current_start = int(meta.get("start", 1))
                current_step = int(meta.get("step", 1000000))
                current_span = int(meta.get("span", 1000000))
                bin_idx = 0
            else:
                bin_start = current_start + bin_idx * current_step
                bin_end = bin_start + current_span
                
                if overlaps_exclusions(current_chrom, bin_start, bin_end):
                    out_f.write("0\n")
                else:
                    out_f.write(line)
                bin_idx += 1

def main():
    if len(sys.argv) < 4:
        print("Usage: filter_wig.py <input.wig> <exclusion.bed> <output.wig>")
        sys.exit(1)
    filter_wig_file(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
