#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
import bisect
import shutil
import subprocess
import tempfile
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Generate high-quality targeted panel MSI loci file.")
    parser.add_argument("--target-regions-bed", required=True, help="Path to panel targeted BED file (.bed or .bed.gz)")
    parser.add_argument("--genome-fasta", help="Path to reference genome FASTA file (for de novo scanning)")
    parser.add_argument("--msi-database", help="Path to genome-wide msi_jitter_sites.38.tsv.gz (for intersection)")
    parser.add_argument("--min-length", type=int, default=10, help="Minimum homopolymer length for de novo scan (default: 10)")
    parser.add_argument("--min-mappability", type=float, default=0.9, help="Minimum mappability for database sites (default: 0.9)")
    parser.add_argument("--output", required=True, help="Path to write the output TSV file")
    
    args = parser.parse_args()
    if not args.genome_fasta and not args.msi_database:
        parser.error("Must provide either --genome-fasta (for de novo scan) or --msi-database (for intersection)")
    return args

def normalize_chrom(chrom):
    """Normalize chromosome names to standard chr1-22, X, Y, M format."""
    c = chrom.strip()
    if not c.startswith("chr"):
        c = "chr" + c
    if c == "chrMT":
        c = "chrM"
    return c

def load_bed_regions(bed_path):
    """Parse BED file and return intervals grouped by chromosome as sorted lists of (start_1_indexed, end_1_indexed)."""
    intervals = defaultdict(list)
    
    open_func = gzip.open if bed_path.endswith(".gz") else open
    mode = "rt" if bed_path.endswith(".gz") else "r"
    
    print(f"Loading target regions from {bed_path}...")
    count = 0
    with open_func(bed_path, mode, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom = normalize_chrom(parts[0])
            try:
                # BED format: 0-indexed start, 1-indexed end
                # Represents genomic bases from start+1 to end (1-indexed)
                start = int(parts[1])
                end = int(parts[2])
                intervals[chrom].append((start + 1, end))
                count += 1
            except ValueError:
                continue
                
    # Sort intervals for binary search
    for chrom in intervals:
        intervals[chrom].sort()
        
    print(f"Loaded {count} intervals across {len(intervals)} chromosomes.")
    return intervals

def is_position_in_bed(intervals, chrom, pos):
    """Check if 1-indexed position falls within any BED interval using binary search."""
    chrom_intervals = intervals.get(chrom)
    if not chrom_intervals:
        return False
        
    # Search for the interval containing pos
    idx = bisect.bisect_right(chrom_intervals, (pos, float('inf'))) - 1
    if idx >= 0:
        start, end = chrom_intervals[idx]
        if start <= pos <= end:
            return True
    return False

def run_intersection(bed_intervals, db_path, min_mappability):
    """Intersect genome-wide database with BED intervals."""
    print(f"Intersecting with genome-wide MSI database {db_path}...")
    matched_sites = []
    
    open_func = gzip.open if db_path.endswith(".gz") else open
    mode = "rt" if db_path.endswith(".gz") else "r"
    
    with open_func(db_path, mode, encoding="utf-8") as f:
        # Read header
        header = f.readline()
        if not header:
            return matched_sites
            
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom = normalize_chrom(parts[0])
            try:
                start = int(parts[1])
                end = int(parts[2])
                mappability = float(parts[4])
                
                if mappability < min_mappability:
                    continue
                    
                # If the repeat overlaps with any BED interval
                # Let's check the start position of the repeat
                if is_position_in_bed(bed_intervals, chrom, start):
                    # In VCF format, the anchor base immediately precedes the repeat
                    anchor_pos = start - 1
                    matched_sites.append((chrom, anchor_pos))
            except ValueError:
                continue
                
    print(f"Found {len(matched_sites)} matching database loci within BED regions.")
    return matched_sites

def scan_fasta_de_novo(bed_intervals, fasta_path, min_length):
    """Scan FASTA file within targeted BED intervals to find homopolymers de novo."""
    matched_sites = []
    
    # Try using samtools faidx for high-performance extraction if available
    if shutil.which("samtools"):
        print("Found 'samtools' in PATH. Running high-performance faidx query...")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as temp_regions:
            for chrom, intervals in bed_intervals.items():
                for start, end in intervals:
                    temp_regions.write(f"{chrom}:{start}-{end}\n")
            temp_regions_path = temp_regions.name
            
        try:
            cmd = ["samtools", "faidx", "-r", temp_regions_path, fasta_path]
            print(f"Executing: {' '.join(cmd)}")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            current_header = None
            seq_lines = []
            
            def parse_header_range(header_str):
                # Header format is like "chr1:101-200"
                # Returns (chrom, start_pos_1_indexed)
                try:
                    parts = header_str.split(":")
                    chrom = normalize_chrom(parts[0])
                    start = int(parts[1].split("-")[0])
                    return chrom, start
                except Exception:
                    return None, None

            def process_region_sequence(header_str, seq_str):
                chrom, start = parse_header_range(header_str)
                if not chrom or start is None or not seq_str:
                    return
                # Scan for homopolymers
                i = 0
                n = len(seq_str)
                while i < n:
                    char = seq_str[i]
                    if char.upper() not in 'ACGT':
                        i += 1
                        continue
                    run_len = 1
                    while i + run_len < n and seq_str[i + run_len].upper() == char.upper():
                        run_len += 1
                    if run_len >= min_length:
                        # Anchor base coordinate is: start + i - 1
                        anchor_pos = start + i - 1
                        matched_sites.append((chrom, anchor_pos))
                    i += run_len

            for line in proc.stdout:
                line = line.strip()
                if line.startswith(">"):
                    if current_header:
                        process_region_sequence(current_header, "".join(seq_lines))
                        seq_lines = []
                    current_header = line[1:]
                else:
                    seq_lines.append(line)
                    
            if current_header:
                process_region_sequence(current_header, "".join(seq_lines))
                
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                print(f"Warning: samtools returned non-zero exit code: {proc.returncode}. Stderr: {stderr}")
        finally:
            if os.path.exists(temp_regions_path):
                os.remove(temp_regions_path)
                
        print(f"Found {len(matched_sites)} homopolymer repeat loci using samtools.")
        return matched_sites

    # Fallback to slow line-by-line parser if samtools is not available
    print("Warning: 'samtools' not found. Falling back to slow line-by-line genome-wide FASTA parsing...")
    current_chrom = None
    chrom_seq_lines = []
    
    def process_chromosome(chrom, seq_str):
        if not seq_str:
            return
        chrom_intervals = bed_intervals.get(chrom)
        if not chrom_intervals:
            return
            
        for start, end in chrom_intervals:
            start_idx = max(0, start - 1)
            end_idx = min(len(seq_str), end)
            
            subseq = seq_str[start_idx:end_idx]
            
            # Find homopolymers
            i = 0
            n = len(subseq)
            while i < n:
                char = subseq[i]
                if char.upper() not in 'ACGT':
                    i += 1
                    continue
                run_len = 1
                while i + run_len < n and subseq[i + run_len].upper() == char.upper():
                    run_len += 1
                if run_len >= min_length:
                    anchor_pos = start + i - 1
                    matched_sites.append((chrom, anchor_pos))
                i += run_len

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if current_chrom:
                    seq_str = "".join(chrom_seq_lines)
                    process_chromosome(current_chrom, seq_str)
                    chrom_seq_lines = []
                header = line.strip()[1:]
                current_chrom = normalize_chrom(header.split()[0])
                print(f"Reading chromosome {current_chrom}...")
            else:
                if current_chrom:
                    chrom_seq_lines.append(line.strip())
                    
        if current_chrom:
            seq_str = "".join(chrom_seq_lines)
            process_chromosome(current_chrom, seq_str)
            
    print(f"Found {len(matched_sites)} de novo homopolymer repeat loci.")
    return matched_sites

def main():
    args = parse_args()
    
    # Load BED target regions
    bed_intervals = load_bed_regions(args.target_regions_bed)
    
    # Generate list of sites
    sites = []
    if args.msi_database:
        sites = run_intersection(bed_intervals, args.msi_database, args.min_mappability)
    elif args.genome_fasta:
        sites = scan_fasta_de_novo(bed_intervals, args.genome_fasta, args.min_length)
        
    # Deduplicate and sort sites
    unique_sites = sorted(list(set(sites)))
    
    # Write output
    print(f"Writing {len(unique_sites)} unique loci to {args.output}...")
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("Chromosome\tPosition\n")
        for chrom, pos in unique_sites:
            f.write(f"{chrom}\t{pos}\n")
            
    print("Done!")

if __name__ == "__main__":
    main()
