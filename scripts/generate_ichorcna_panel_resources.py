#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
import subprocess
import tempfile
import multiprocessing
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Generate empirical exclusion BED and custom Panel of Normals (PoN) for off-target ichorCNA.")
    parser.add_argument("--bams", nargs="+", required=True, help="Paths to normal BAM/CRAM files used to build resource")
    parser.add_argument("--chrom-sizes", required=True, help="Path to chromosome sizes file (chrom\\tsize)")
    parser.add_argument("--output-dir", required=True, help="Output directory to save generated files")
    parser.add_argument("--threshold", type=float, default=None, help="Coverage threshold to identify baited regions (default: auto-determined)")
    parser.add_argument("--target-fraction", type=float, default=0.02, help="Empirical genome coverage fraction to determine threshold (default: 0.02)")
    parser.add_argument("--padding", type=int, default=250, help="Base pair padding around identified baited regions (default: 250)")
    parser.add_argument("--window-size", type=int, default=1000000, help="Bin window size in bp for readCounter/ichorCNA (default: 1000000)")
    parser.add_argument("--quality", type=int, default=20, help="Minimum mapping quality filter for readCounter (default: 20)")
    parser.add_argument("--chrs", default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY", help="Comma-separated chromosomes to analyze")
    parser.add_argument("--mosdepth-bin", default="mosdepth", help="Path to mosdepth binary")
    parser.add_argument("--readcounter-bin", default="readCounter", help="Path to readCounter binary")
    default_rscript = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources", "runIchorCNA.R")
    parser.add_argument("--ichor-rscript", default=default_rscript, help="Path to runIchorCNA.R Rscript")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads/jobs (default: 4)")
    return parser.parse_args()

def load_chrom_sizes(path):
    sizes = {}
    with open(path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes

def run_mosdepth_sample(args_tuple):
    bam, threshold, mosdepth_bin, out_dir = args_tuple
    sample_name = os.path.splitext(os.path.basename(bam))[0]
    prefix = os.path.join(out_dir, f"tmp_{sample_name}_mosdepth")
    
    print(f"Running mosdepth for {sample_name}...")
    cmd = [
        mosdepth_bin,
        "--by", "100",
        prefix,
        bam
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running mosdepth for {sample_name}: {e.stderr.decode()}", file=sys.stderr)
        return []

    # Read and parse regions bed.gz
    regions_file = f"{prefix}.regions.bed.gz"
    intervals = defaultdict(list)
    if os.path.exists(regions_file):
        with gzip.open(regions_file, 'rt') as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    cov = float(parts[3])
                    if cov >= threshold:
                        intervals[chrom].append((start, end))
                        
    # Clean up temp mosdepth files
    for suffix in [".regions.bed.gz", ".regions.bed.gz.csi", ".mosdepth.global.dist.txt", ".mosdepth.summary.txt"]:
        fpath = prefix + suffix
        if os.path.exists(fpath):
            os.remove(fpath)
            
    return intervals

def merge_intervals(intervals_list):
    if not intervals_list:
        return []
    intervals_list.sort(key=lambda x: (x[0], x[1]))
    merged = [intervals_list[0]]
    for current in intervals_list[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return merged

def run_readcounter_sample(args_tuple):
    bam, readcounter_bin, window_size, quality, chrs, out_dir, exclusion_bed = args_tuple
    sample_name = os.path.splitext(os.path.basename(bam))[0]
    out_wig = os.path.join(out_dir, f"{sample_name}.wig")
    
    # 1. Estimate off-target coverage of the normal sample
    total_reads = 0
    try:
        res = subprocess.run(["samtools", "idxstats", bam], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for line in res.stdout.decode().splitlines():
            parts = line.split()
            if len(parts) >= 3:
                total_reads += int(parts[2])
    except Exception as e:
        print(f"Warning: could not estimate total reads for {sample_name}: {e}", file=sys.stderr)
        
    on_target_reads = 0
    try:
        res = subprocess.run(["samtools", "view", "-c", "-q", str(quality), "-L", exclusion_bed, bam], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        on_target_reads = int(res.stdout.decode().strip())
    except Exception as e:
        print(f"Warning: could not estimate on-target reads for {sample_name}: {e}", file=sys.stderr)
        
    off_target_reads = total_reads - on_target_reads
    if off_target_reads < 0:
        off_target_reads = 0
        
    # Estimate off-target coverage (assume 150bp read length and 3.1e9 genome size)
    off_target_cov = (off_target_reads * 150) / 3100000000
    print(f"Normal {sample_name} estimated off-target coverage: {off_target_cov:.3f}x")
    
    downsample_bam = bam
    is_temp_bam = False
    
    if off_target_cov > 1.0:
        fraction = 1.0 / off_target_cov
        s_param = f"100.{int(fraction * 1000000):06d}"
        print(f"Downsampling normal {sample_name} off-target reads to 1.0x (fraction: {fraction:.4f})...")
        
        temp_bam = os.path.join(out_dir, f"tmp_downsampled_{sample_name}.bam")
        try:
            with open(temp_bam, 'wb') as out_f:
                subprocess.run(["samtools", "view", "-b", "-s", s_param, "-q", str(quality), bam], check=True, stdout=out_f, stderr=subprocess.PIPE)
            subprocess.run(["samtools", "index", temp_bam], check=True)
            downsample_bam = temp_bam
            is_temp_bam = True
        except Exception as e:
            print(f"Warning: downsampling failed for normal {sample_name}: {e}", file=sys.stderr)
            if os.path.exists(temp_bam):
                os.remove(temp_bam)
            downsample_bam = bam
            is_temp_bam = False
            
    print(f"Running readCounter for {sample_name}...")
    cmd = [
        readcounter_bin,
        "--window", str(window_size),
        "--quality", str(quality),
        "--chromosome", chrs,
        downsample_bam
    ]
    try:
        with open(out_wig, 'w') as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error running readCounter for {sample_name}: {e.stderr.decode()}", file=sys.stderr)
        out_wig = None
    finally:
        if is_temp_bam:
            if os.path.exists(downsample_bam):
                os.remove(downsample_bam)
            if os.path.exists(downsample_bam + ".bai"):
                os.remove(downsample_bam + ".bai")
                
    return out_wig

def filter_wig_file(wig_path, exclusion_bed_path, out_wig_path):
    print(f"Filtering baited regions in {os.path.basename(wig_path)} using exclusion BED...")
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
        
    # Helper to check if a bin overlaps with any exclusion region
    # Since both bins and exclusions are sorted, we can search efficiently.
    def overlaps_exclusions(chrom, bin_start, bin_end):
        excl_list = exclusions.get(chrom, [])
        if not excl_list:
            return False
        # Binary search for the first exclusion starting after or near bin_start
        import bisect
        idx = bisect.bisect_left(excl_list, (bin_start, 0))
        # Check current index and previous index
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
                # Value row
                # Calculate coordinates for fixedStep bin
                bin_start = current_start + bin_idx * current_step
                bin_end = bin_start + current_span
                
                # Check overlap
                if overlaps_exclusions(current_chrom, bin_start, bin_end):
                    out_f.write("0\n")
                else:
                    out_f.write(line)
                bin_idx += 1

def determine_empirical_threshold(bam, mosdepth_bin, out_dir, target_fraction=0.02):
    sample_name = os.path.splitext(os.path.basename(bam))[0]
    prefix = os.path.join(out_dir, f"pilot_{sample_name}_mosdepth")
    cmd = [mosdepth_bin, "--by", "100", prefix, bam]
    print(f"Running pilot mosdepth on {sample_name} to determine empirical threshold...")
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Warning: Pilot mosdepth failed: {e.stderr.decode()}", file=sys.stderr)
        return 50.0  # Fallback threshold
        
    dist_file = prefix + ".mosdepth.global.dist.txt"
    if not os.path.exists(dist_file):
        print("Warning: Pilot mosdepth global distribution file not found.", file=sys.stderr)
        return 50.0
        
    empirical_threshold = 50.0
    with open(dist_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0] == "total":
                depth = float(parts[1])
                fraction = float(parts[2])
                if fraction <= target_fraction:
                    empirical_threshold = depth
                    break
                    
    # Clean up pilot files
    for suffix in [".regions.bed.gz", ".regions.bed.gz.csi", ".mosdepth.global.dist.txt", ".mosdepth.summary.txt"]:
        p = prefix + suffix
        if os.path.exists(p):
            os.remove(p)
            
    # Ensure we don't return too low a threshold (minimum 10.0 to prevent noise)
    empirical_threshold = max(10.0, empirical_threshold)
    print(f"Empirically determined threshold: {empirical_threshold} (covers <= {target_fraction*100}% of the genome)")
    return empirical_threshold

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. Load chrom sizes
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    
    exclusion_bed_path = os.path.join(args.output_dir, "custom_exclusion.bed")
    if os.path.exists(exclusion_bed_path):
        print(f"Step 1 & 2: Exclusion BED already exists at {exclusion_bed_path}. Skipping mosdepth run.")
    else:
        threshold = args.threshold
        if threshold is None:
            threshold = determine_empirical_threshold(args.bams[0], args.mosdepth_bin, args.output_dir, args.target_fraction)
            
        # 2. Parallel run mosdepth to discover baited regions
        mosdepth_args = [(bam, threshold, args.mosdepth_bin, args.output_dir) for bam in args.bams]
        print(f"Step 1: Identifying high-coverage baited regions in parallel using {args.threads} threads (threshold: {threshold})...")
        with multiprocessing.Pool(args.threads) as pool:
            results = pool.map(run_mosdepth_sample, mosdepth_args)
            
        # Merge regions across all samples
        all_regions = defaultdict(list)
        for res in results:
            for chrom, intervals in res.items():
                all_regions[chrom].extend(intervals)
                
        print("Step 2: Merging and padding identified baited regions...")
        with open(exclusion_bed_path, 'w') as out_f:
            for chrom, intervals in all_regions.items():
                # Merge overlapping intervals
                merged = merge_intervals(intervals)
                # Apply padding
                padded = []
                max_size = chrom_sizes.get(chrom, sys.maxsize)
                for start, end in merged:
                    padded_start = max(0, start - args.padding)
                    padded_end = min(max_size, end + args.padding)
                    padded.append((padded_start, padded_end))
                # Merge overlapping padded intervals
                final_merged = merge_intervals(padded)
                # Write to BED
                for start, end in final_merged:
                    out_f.write(f"{chrom}\t{start}\t{end}\n")
                    
        print(f"Custom exclusion BED saved to: {exclusion_bed_path}")
    
    # 3. Parallel run readCounter to generate WIG files for missing samples
    needed_bams = []
    filtered_wig_list = []
    for bam in args.bams:
        sample_name = os.path.splitext(os.path.basename(bam))[0]
        filtered_wig = os.path.join(args.output_dir, f"{sample_name}.filtered.wig")
        filtered_wig_list.append(filtered_wig)
        if not os.path.exists(filtered_wig):
            needed_bams.append(bam)
            
    if needed_bams:
        readcounter_args = [(bam, args.readcounter_bin, args.window_size, args.quality, args.chrs, args.output_dir, exclusion_bed_path) for bam in needed_bams]
        print(f"Step 3: Generating WIG files for {len(needed_bams)} missing samples in parallel using {args.threads} threads...")
        with multiprocessing.Pool(args.threads) as pool:
            wig_files = pool.map(run_readcounter_sample, readcounter_args)
            
        # 4. Filter WIG files using exclusion BED
        for wig_file in wig_files:
            if not wig_file:
                continue
            filtered_wig = wig_file.replace(".wig", ".filtered.wig")
            filter_wig_file(wig_file, exclusion_bed_path, filtered_wig)
            # Clean up unfiltered WIG
            if os.path.exists(wig_file):
                os.remove(wig_file)
    else:
        print("Step 3 & 4: All filtered WIG files already exist. Skipping readCounter and filtering.")
            
    # 5. Create Panel of Normals (PoN)
    print("Step 4: Creating Panel of Normals (PoN) via ichorCNA...")
    wig_list_path = os.path.join(args.output_dir, "normal_wig_files.txt")
    with open(wig_list_path, 'w') as f:
        for w in filtered_wig_list:
            f.write(w + "\n")
            
    pon_rds_path = os.path.join(args.output_dir, "custom_panelOfNormals.rds")
    
    # Determine window size string (e.g. 1000kb)
    window_size_kb = f"{args.window_size // 1000}kb"
    # Determine genome build (hg19 or hg38)
    genome_build = "hg38"
    if "37" in args.chrom_sizes or "19" in args.chrom_sizes:
        genome_build = "hg19"
        
    pon_cmd = [
        "Rscript", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources", "create_pon_custom.R"),
        "--filelist", wig_list_path,
        "--outfile", pon_rds_path,
        "--genomeBuild", genome_build,
        "--windowSize", window_size_kb,
        "--genomeStyle", "UCSC"
    ]
    print(f"Running command: {' '.join(pon_cmd)}")
    try:
        subprocess.run(pon_cmd, check=True)
        print(f"PoN RDS successfully created at: {pon_rds_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error creating Panel of Normals: {e}", file=sys.stderr)
        sys.exit(1)
        
    print("All tasks completed successfully!")

if __name__ == "__main__":
    main()
