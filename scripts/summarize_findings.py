#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import gzip
import json
import subprocess
import shlex
import time
import glob
import re

# Genes of interest in the exact order requested
GENES_OF_INTEREST = [
    'MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM', 'MLH3', 'POLD1', 'POLE', 'MUTYH', 'TP53', 'BRAF'
]

GENE_BOUNDS = {
    "GRCh38_hmf": {
        "MLH1":  {"chrom": "chr3",  "start": 36992834, "end": 37051505, "strand": 1},
        "MSH2":  {"chrom": "chr2",  "start": 47403067, "end": 47663146, "strand": 1},
        "MSH6":  {"chrom": "chr2",  "start": 47695530, "end": 47810109, "strand": 1},
        "PMS2":  {"chrom": "chr7",  "start": 5970925,  "end": 6009130,  "strand": -1},
        "EPCAM": {"chrom": "chr2",  "start": 47345158, "end": 47387601, "strand": 1},
        "MLH3":  {"chrom": "chr14", "start": 75013224, "end": 75051532, "strand": -1},
        "POLD1": {"chrom": "chr19", "start": 50384204, "end": 50418018, "strand": 1},
        "POLE":  {"chrom": "chr12", "start": 132623753, "end": 132687376, "strand": -1},
        "MUTYH": {"chrom": "chr1",  "start": 45329163, "end": 45340893, "strand": -1},
        "TP53":  {"chrom": "chr17", "start": 7661779,  "end": 7687546,  "strand": -1},
        "BRAF":  {"chrom": "chr7",  "start": 140719327, "end": 140925199, "strand": -1}
    },
    "GRCh37_hmf": {
        "MLH1":  {"chrom": "chr3",  "start": 37034823, "end": 37107380, "strand": 1},
        "MSH2":  {"chrom": "chr2",  "start": 47630108, "end": 47789450, "strand": 1},
        "MSH6":  {"chrom": "chr2",  "start": 47922669, "end": 48037240, "strand": 1},
        "PMS2":  {"chrom": "chr7",  "start": 6012870,  "end": 6048756,  "strand": -1},
        "EPCAM": {"chrom": "chr2",  "start": 47572297, "end": 47614740, "strand": 1},
        "MLH3":  {"chrom": "chr14", "start": 75480467, "end": 75518235, "strand": -1},
        "POLD1": {"chrom": "chr19", "start": 50887461, "end": 50921273, "strand": 1},
        "POLE":  {"chrom": "chr12", "start": 133200348, "end": 133263951, "strand": -1},
        "MUTYH": {"chrom": "chr1",  "start": 45794835, "end": 45806142, "strand": -1},
        "TP53":  {"chrom": "chr17", "start": 7565097,  "end": 7590856,  "strand": -1},
        "BRAF":  {"chrom": "chr7",  "start": 140419127, "end": 140624564, "strand": -1}
    }
}

def parse_samplesheet(samplesheet_path):
    """Parse Nextflow samplesheet CSV to extract tumor sample information.
    Returns a list of dicts with keys: group_id, subject_id, sample_id
    """
    tumor_samples = []
    if not os.path.exists(samplesheet_path):
        print(f"Error: Samplesheet not found at {samplesheet_path}", file=sys.stderr)
        return tumor_samples

    with open(samplesheet_path, mode='r', encoding='utf-8') as f:
        # Some samplesheets might have comments starting with #
        lines = [line for line in f if not line.strip().startswith('#')]
        reader = csv.DictReader(lines)
        for row in reader:
            if row.get('sample_type') == 'tumor' or row.get('sample_type', '').lower() == 'tumor':
                tumor_samples.append({
                    'group_id': row.get('group_id', '').strip(),
                    'subject_id': row.get('subject_id', '').strip(),
                    'sample_id': row.get('sample_id', '').strip()
                })
    
    # Deduplicate by group_id
    seen = set()
    dedup_samples = []
    for s in tumor_samples:
        if s['group_id'] and s['group_id'] not in seen:
            seen.add(s['group_id'])
            dedup_samples.append(s)
            
    return dedup_samples

def get_sample_filepaths(samplesheet_path, group_id, purple_dir):
    """Find all bam/cram track filepaths for the given group_id.
    Checks the samplesheet first, and if empty/fastq, falls back to scanning the alignments directory.
    """
    paths = []
    if os.path.exists(samplesheet_path):
        try:
            with open(samplesheet_path, mode='r', encoding='utf-8') as f:
                lines = [line for line in f if not line.strip().startswith('#')]
                reader = csv.DictReader(lines)
                for row in reader:
                    if row.get('group_id', '').strip() == group_id:
                        p_list = row.get('filepath', '').split(';')
                        for path in p_list:
                            path = path.strip()
                            if path.lower().endswith(('.bam', '.cram')) and os.path.exists(path):
                                paths.append(path)
        except Exception as e:
            print(f"Warning: Failed to parse samplesheet: {e}", file=sys.stderr)
            
    # Fallback to scanning the alignments directory in the results output folder
    if not paths and purple_dir:
        parent_dir = os.path.dirname(purple_dir)
        search_dirs = [
            os.path.join(parent_dir, 'alignments'),
            os.path.join(parent_dir, 'alignments', 'dna'),
            os.path.join(parent_dir, group_id, 'alignments'),
            os.path.join(parent_dir, group_id, 'alignments', 'dna'),
            os.path.join(os.path.dirname(parent_dir), 'alignments'),
            os.path.join(os.path.dirname(parent_dir), group_id, 'alignments'),
            os.path.join(os.path.dirname(parent_dir), group_id, 'alignments', 'dna'),
        ]
        for align_dir in search_dirs:
            if os.path.exists(align_dir):
                for root, dirs, files in os.walk(align_dir):
                    for f in files:
                        if f.lower().endswith(('.bam', '.cram')) and not f.startswith('.'):
                            paths.append(os.path.join(root, f))
                    
    # Filter to avoid duplicate bam/cram for the same sample
    grouped = {}
    for p in sorted(list(set(paths))):
        base, ext = os.path.splitext(p)
        ext = ext.lower()
        if ext in ['.bam', '.cram']:
            if base not in grouped:
                grouped[base] = []
            grouped[base].append(p)
            
    filtered_paths = []
    for base, plist in grouped.items():
        crams = [x for x in plist if x.lower().endswith('.cram')]
        if crams:
            filtered_paths.append(crams[0])
        else:
            filtered_paths.append(plist[0])
            
    return sorted(filtered_paths)

def generate_igv_bed_and_report(group_id, subject_id, purple_dir, pave_dir, samplesheet_path, genome, methylation_data, summary_dir, fasta_path, gtf_path, igv_env, tumor_sample_id=None, normal_sample_id=None):
    """Generate an IGV BED file and run create_report to create an HTML report."""
    # 1. Parse all variants
    all_variants = []
    
    # We must retrieve bounds for chrom
    bounds_subset = GENE_BOUNDS.get(genome, GENE_BOUNDS['GRCh38_hmf'])
    
    for gene in GENES_OF_INTEREST:
        variants = get_detailed_variants(purple_dir, pave_dir, group_id, gene, subject_id, samplesheet_path, genome, methylation_data, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id)
        if variants:
            chrom = bounds_subset.get(gene, {}).get('chrom', 'chr1')
            for var in variants:
                if 'methylation' in str(var.get('var_type', '')).lower():
                    continue
                all_variants.append({
                    'chrom': chrom,
                    'pos': var['pos'],
                    'gene': gene,
                    'label': var['label'],
                    'type': var['var_type']
                })
                
    if not all_variants:
        print(f"No variants found for sample {group_id}, skipping IGV report.", file=sys.stderr)
        return
        
    # Create summary directory if needed
    os.makedirs(os.path.join(summary_dir, 'igv'), exist_ok=True)
    bed_path = os.path.join(summary_dir, 'igv', f"{group_id}_igv_report.bed")
    
    # 2. Write BED file
    try:
        with open(bed_path, 'w', encoding='utf-8') as f:
            for var in all_variants:
                # 0-based start, 1-based end
                start = max(0, var['pos'] - 1)
                end = var['pos']
                name = "{}:{}_{}".format(var['gene'], var['label'].replace(" ", "_"), var['type'])
                f.write(f"{var['chrom']}\t{start}\t{end}\t{name}\n")
        print(f"IGV BED file successfully written to {bed_path}", file=sys.stderr)
    except Exception as e:
        print(f"Error: Failed to write IGV BED file {bed_path}: {e}", file=sys.stderr)
        return
        
    # 3. Locate tracks from samplesheet or alignments folder
    tracks = get_sample_filepaths(samplesheet_path, group_id, purple_dir)
    if not tracks:
        print(f"Warning: No BAM/CRAM track paths found in samplesheet for {group_id}, skipping IGV report.", file=sys.stderr)
        return
        
    # 4. Build command
    output_html = os.path.join(summary_dir, 'igv', f"{group_id}_igv_report.html")
    
    # Check for direct create_report binary to avoid conda 4.6.8 read-only root_prefix PermissionError
    create_report_bin = None
    candidate_bins = [
        os.path.expanduser(f"~/.conda/envs/{igv_env}/bin/create_report"),
        f"/cluster/home/jbruce/.conda/envs/{igv_env}/bin/create_report",
        f"/cluster/tools/software/centos7/python/3.7.2/envs/{igv_env}/bin/create_report",
        f"/cluster/tools/software/rocky9/envs/{igv_env}/bin/create_report"
    ]
    for cb in candidate_bins:
        if os.path.exists(cb):
            create_report_bin = cb
            break

    if create_report_bin:
        cmd = [
            create_report_bin, bed_path,
            "--flanking", "2000",
            "--output", output_html
        ]
    else:
        cmd = [
            "conda", "run", "-n", igv_env,
            "create_report", bed_path,
            "--flanking", "2000",
            "--output", output_html
        ]
    
    if fasta_path:
        cmd.extend(["--fasta", fasta_path])
        
    cmd.append("--tracks")
    cmd.extend(tracks)
    if gtf_path and os.path.exists(gtf_path):
        cmd.append(gtf_path)
    
    print(f"Running IGV report command: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            print(f"IGV HTML report successfully written to {output_html}", file=sys.stderr)
            inject_igv_navigation_script(output_html)
        else:
            print(f"Warning: create_report failed for {group_id}: {res.stderr}", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Failed to execute create_report for {group_id}: {e}", file=sys.stderr)

def get_tumor_sample_id(samplesheet_path, group_id):
    """Parse samplesheet to find the matching tumor sample_id for the given group_id."""
    if not os.path.exists(samplesheet_path):
        return None
    with open(samplesheet_path, mode='r', encoding='utf-8') as f:
        lines = [line for line in f if not line.strip().startswith('#')]
        reader = csv.DictReader(lines)
        for row in reader:
            r_group = row.get('group_id', '').strip()
            r_type = row.get('sample_type', '').strip().lower()
            if r_group == group_id and r_type == 'tumor':
                return row.get('sample_id', '').strip()
    return None

def get_normal_sample_id(samplesheet_path, group_id):
    """Parse samplesheet to find the matching normal sample_id for the given group_id."""
    if not os.path.exists(samplesheet_path):
        return None
    with open(samplesheet_path, mode='r', encoding='utf-8') as f:
        lines = [line for line in f if not line.strip().startswith('#')]
        reader = csv.DictReader(lines)
        for row in reader:
            r_group = row.get('group_id', '').strip()
            r_type = row.get('sample_type', '').strip().lower()
            if r_group == group_id and r_type == 'normal':
                return row.get('sample_id', '').strip()
    return None

def is_flexible_match(col_name, sample_id, sample_type):
    """Match a methylation matrix column header flexibly to a target samplesheet sample_id."""
    if not col_name or not sample_id:
        return False
    col = col_name.strip().lower()
    sid = sample_id.strip().lower()
    
    # 1. Exact Match
    if col == sid:
        return True
        
    # Remove non-alphanumeric chars for robust substring matching
    import re
    col_clean = re.sub(r'[^a-z0-9]', '', col)
    sid_clean = re.sub(r'[^a-z0-9]', '', sid)
    
    # Check sample type indicator matches
    col_parts = set(re.split(r'[-_]', col))
    is_tumor_col = (
        't' in col_parts or 'tumor' in col or 'tumour' in col or 
        col.endswith('t') or col.endswith('_t') or col.endswith('-t') or
        bool(re.match(r'^t\d*[-_]', col))
    )
    is_normal_col = (
        'n' in col_parts or 'normal' in col or 'control' in col or 
        col.endswith('n') or col.endswith('_n') or col.endswith('-n') or
        bool(re.match(r'^n\d*[-_]', col))
    )
    
    is_tumor_sid = sample_type.lower() == 'tumor'
    is_normal_sid = sample_type.lower() == 'normal'
    
    if is_tumor_sid and not is_tumor_col:
        return False
    if is_normal_sid and not is_normal_col:
        return False
    if is_tumor_sid and is_normal_col and not is_tumor_col:
        return False
    if is_normal_sid and is_tumor_col and not is_normal_col:
        return False
        
    # Check donor token match (e.g. 22-00256, 24-00003, etc.)
    col_donor = re.search(r'\d{2}-\d{4,5}', col)
    sid_donor = re.search(r'\d{2}-\d{4,5}', sid)
    if col_donor and sid_donor and col_donor.group(0) == sid_donor.group(0):
        return True

    # Patient number match with donor prefix (e.g. T3-.. vs T3-..)
    col_prefix = re.match(r'^([tn]\d+)', col)
    sid_prefix = re.match(r'^([tn]\d+)', sid)
    if col_prefix and sid_prefix and col_prefix.group(1) == sid_prefix.group(1):
        return True
        
    # Check if the clean base ID (excluding type suffixes) matches
    col_base_clean = col_clean.replace('tumor', '').replace('normal', '').replace('tumour', '').strip('tn')
    sid_base_clean = sid_clean.replace('tumor', '').replace('normal', '').replace('tumour', '').strip('tn')
    
    if col_base_clean and (col_base_clean in sid_clean or sid_base_clean in col_clean):
        return True
        
    # Check if numeric portions match (e.g. UHN-105729 vs 105729)
    col_nums = "".join(re.findall(r'\d+', col))
    sid_nums = "".join(re.findall(r'\d+', sid))
    if col_nums and sid_nums and (col_nums in sid_nums or sid_nums in col_nums):
        return True
        
    return False

def get_sample_ids_from_samplesheet(samplesheet_path, group_id):
    """Resolve tumor_sample_id and normal_sample_id for a given group_id from the samplesheet."""
    tumor_sample_id = None
    normal_sample_id = None
    if samplesheet_path and os.path.exists(samplesheet_path):
        try:
            with open(samplesheet_path, mode='r', encoding='utf-8-sig') as f:
                lines = [l for l in f if l.strip() and not l.strip().startswith('#')]
                reader = csv.DictReader(lines)
                for r in reader:
                    if r.get('group_id', '').strip() == group_id:
                        st = (r.get('sample_type') or '').strip().lower()
                        sid = (r.get('sample_id') or '').strip()
                        if st == 'tumor' and not tumor_sample_id:
                            tumor_sample_id = sid
                        elif st == 'normal' and not normal_sample_id:
                            normal_sample_id = sid
        except Exception:
            pass
    return tumor_sample_id, normal_sample_id

def find_purple_file(purple_dir, suffix, group_id, sample_id=None):
    """Locate a PURPLE file matching a given suffix (e.g. 'purity.tsv', 'driver.catalog.somatic.tsv').
    Searches for {sample_id}.purple.{suffix}, {group_id}.purple.{suffix}, or wildcards.
    """
    if not purple_dir or not os.path.isdir(purple_dir):
        return None
    candidates = []
    if sample_id:
        candidates.append(os.path.join(purple_dir, f"{sample_id}.purple.{suffix}"))
    if group_id and group_id != sample_id:
        candidates.append(os.path.join(purple_dir, f"{group_id}.purple.{suffix}"))
    for c in candidates:
        if os.path.isfile(c):
            return c
    # Glob fallback
    matches = glob.glob(os.path.join(purple_dir, f"*.purple.{suffix}"))
    if matches:
        return sorted(matches)[0]
    return None

def find_pave_vcf(pave_dir, vcf_type, group_id, sample_id=None):
    """Locate a PAVE VCF for vcf_type ('somatic' or 'germline').
    Supports nested subdirectories (pave/somatic/, pave/germline/) and flat pave/ dir.
    """
    if not pave_dir or not os.path.isdir(pave_dir):
        return None
    subdirs = [vcf_type, '']
    for sub in subdirs:
        d = os.path.join(pave_dir, sub) if sub else pave_dir
        if not os.path.isdir(d):
            continue
        candidates = []
        if sample_id:
            candidates.append(os.path.join(d, f"{sample_id}.pave.{vcf_type}.vcf.gz"))
        if group_id and group_id != sample_id:
            candidates.append(os.path.join(d, f"{group_id}.pave.{vcf_type}.vcf.gz"))
        for c in candidates:
            if os.path.isfile(c):
                return c
        matches = glob.glob(os.path.join(d, f"*.pave.{vcf_type}.vcf.gz"))
        if matches:
            return sorted(matches)[0]
    return None

def find_orange_file(orange_dir, suffix, group_id, sample_id=None):
    """Locate Orange report files (e.g. '.orange.pdf', '.orange.json')."""
    if not orange_dir or not os.path.isdir(orange_dir):
        return None
    candidates = []
    if sample_id:
        candidates.append(os.path.join(orange_dir, f"{sample_id}{suffix}"))
    if group_id and group_id != sample_id:
        candidates.append(os.path.join(orange_dir, f"{group_id}{suffix}"))
    for c in candidates:
        if os.path.isfile(c):
            return c
    matches = glob.glob(os.path.join(orange_dir, f"*{suffix}"))
    if matches:
        return sorted(matches)[0]
    return None

def get_purity(purple_dir, group_id, sample_id=None):
    """Parse purity from *.purple.purity.tsv"""
    purity_file = find_purple_file(purple_dir, "purity.tsv", group_id, sample_id=sample_id)
    if not purity_file or not os.path.exists(purity_file):
        return "NA"
    
    try:
        with open(purity_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                purity_val = row.get('purity')
                if purity_val is not None:
                    # Convert to float and format to 4 decimal places if possible
                    return str(round(float(purity_val), 5))
    except Exception as e:
        print(f"Warning: Failed to parse purity file {purity_file}: {e}", file=sys.stderr)
    
    return "NA"

def get_msi_status(purple_dir, group_id, sample_id=None):
    """Parse MSI status from *.purple.purity.tsv"""
    purity_file = find_purple_file(purple_dir, "purity.tsv", group_id, sample_id=sample_id)
    if not purity_file or not os.path.exists(purity_file):
        return "NA"
    
    try:
        with open(purity_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                msi_status = row.get('msStatus')
                if msi_status is not None:
                    return msi_status.strip()
    except Exception as e:
        print(f"Warning: Failed to parse MSI status from {purity_file}: {e}", file=sys.stderr)
    
    return "NA"

def get_msi_unstable_sites_proportion(purple_dir, group_id, sample_id=None):
    """Parse MSI site count and unstable count from purple.command.log.
    Returns: (msi_indels, total_msi_sites, proportion) or (None, None, None)
    """
    import re
    # Locate the purple command log
    parent_dir = os.path.dirname(purple_dir)
    logs_dir = os.path.join(parent_dir, 'logs')
    candidates = []
    if sample_id:
        candidates.extend([
            os.path.join(logs_dir, f"{sample_id}.purple.command.log"),
            os.path.join(purple_dir, f"{sample_id}.purple.command.log")
        ])
    candidates.extend([
        os.path.join(logs_dir, f"{group_id}.purple.command.log"),
        os.path.join(purple_dir, f"{group_id}.purple.command.log"),
        os.path.join(purple_dir, "purple.command.log")
    ])
    
    log_file = None
    for c in candidates:
        if os.path.exists(c):
            log_file = c
            break
    if not log_file:
        matches = glob.glob(os.path.join(purple_dir, "*.purple.command.log"))
        if not matches:
            matches = glob.glob(os.path.join(logs_dir, "*.purple.command.log"))
        if matches:
            log_file = matches[0]
            
    if not log_file or not os.path.exists(log_file):
        return None, None, None
        
    msi_indels = None
    total_msi_sites = None
    
    try:
        with open(log_file, 'r', encoding='utf-8') as f:
            for line in f:
                if "loaded" in line and "MSI indels from file" in line:
                    parts = line.strip().split()
                    for i, part in enumerate(parts):
                        if part == "loaded" and i + 1 < len(parts):
                            try:
                                total_msi_sites = int(parts[i+1])
                            except ValueError:
                                pass
                if "msiIndels(" in line:
                    # Look for load(...) msiIndels(Y perMb=...) line
                    match = re.search(r'msiIndels\((\d+(?:\.\d+)?)\s+perMb=', line)
                    if match:
                        try:
                            msi_indels = int(float(match.group(1)))
                        except ValueError:
                            pass
    except Exception as e:
        print(f"Warning: Failed to parse PURPLE log {log_file}: {e}", file=sys.stderr)
        
    if msi_indels is not None and total_msi_sites is not None and total_msi_sites > 0:
        proportion = msi_indels / total_msi_sites
        return msi_indels, total_msi_sites, proportion
        
    return None, None, None

def get_germline_msi_proportion(purple_dir, group_id, msi_sites_file=None, sample_id=None):
    """
    Parse the germline VCF file to count PASS indels at the MSI sites.
    Returns: (germline_indels, total_msi_sites, proportion) or (None, None, None)
    """
    if msi_sites_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        msi_sites_file = os.path.join(os.path.dirname(script_dir), "resources", "panels", "multimmr_v3", "MMR_msi_sites.v3.tsv")
        
    if not os.path.exists(msi_sites_file):
        # Fallback to v2 path if v3 is not found
        msi_sites_file = os.path.join(os.path.dirname(script_dir), "resources", "panels", "multimmr_v2", "MMR_msi_sites.v2.tsv")
        
    if not os.path.exists(msi_sites_file):
        # Additional fallback to locate anywhere in panels folder
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for root, dirs, files in os.walk(os.path.join(parent_dir, "resources", "panels")):
            for file in files:
                if "msi_sites" in file and file.endswith(".tsv"):
                    msi_sites_file = os.path.join(root, file)
                    break
        
    if not os.path.exists(msi_sites_file):
        return None, None, None
        
    msi_sites = set()
    try:
        with open(msi_sites_file, 'r', encoding='utf-8') as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    msi_sites.add((parts[0], int(parts[1])))
    except Exception as e:
        print(f"Warning: Failed to read MSI sites file {msi_sites_file}: {e}", file=sys.stderr)
        return None, None, None
        
    if not msi_sites:
        return None, None, None
        
    import gzip
    vcf_file = find_purple_file(purple_dir, "germline.vcf.gz", group_id, sample_id=sample_id)
    if not vcf_file or not os.path.exists(vcf_file):
        return None, len(msi_sites), 0.0
        
    germline_indels = 0
    try:
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                if len(ref) != len(alt):
                    if (chrom, pos) in msi_sites:
                        filt = parts[6]
                        if filt == 'PASS':
                            germline_indels += 1
    except Exception as e:
        print(f"Warning: Failed to parse germline VCF {vcf_file}: {e}", file=sys.stderr)
        return None, len(msi_sites), None
        
    proportion = germline_indels / len(msi_sites)
    return germline_indels, len(msi_sites), proportion

def parse_copynumber_and_loh(purple_dir, group_id, sample_id=None):
    """Parse cnv.gene.tsv for copy number changes and LOH.
    Returns a dict of gene -> {copynumber: status_str, LOH: status_str}
    """
    cnv_file = find_purple_file(purple_dir, "cnv.gene.tsv", group_id, sample_id=sample_id)
    results = {gene: {'copynumber': None, 'LOH': None, 'min_cn': 2.0, 'max_cn': 2.0} for gene in GENES_OF_INTEREST}
    
    if not cnv_file or not os.path.exists(cnv_file):
        return results
        
    try:
        with open(cnv_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row.get('gene', '').strip()
                # Check isCanonical if it's there
                is_canonical = row.get('isCanonical', 'true').lower() == 'true'
                
                if gene_name in results and is_canonical:
                    try:
                        min_cn = float(row.get('minCopyNumber', 2.0))
                        max_cn = float(row.get('maxCopyNumber', 2.0))
                        min_minor_cn = float(row.get('minMinorAlleleCopyNumber', 1.0))
                        
                        results[gene_name]['min_cn'] = min_cn
                        results[gene_name]['max_cn'] = max_cn
                        
                        # Copy number classification
                        if min_cn < 0.5:
                            results[gene_name]['copynumber'] = 'Homozygous Deletion'
                        elif min_cn < 1.7:
                            results[gene_name]['copynumber'] = 'Heterozygous Loss'
                        elif max_cn > 2.8:
                            results[gene_name]['copynumber'] = f'Gain (CN={max_cn:.2f})'
                        else:
                            results[gene_name]['copynumber'] = 'Normal'
                            
                        # LOH classification (minor allele copy number < 0.35 for panel sensitivity)
                        if min_minor_cn < 0.35:
                            results[gene_name]['LOH'] = 'LOH'
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Warning: Failed to parse CNV file {cnv_file}: {e}", file=sys.stderr)

    # Segment-level LOH fallback for panel data (where sparse SNPs or segment overshooting can smooth gene-level minor CN)
    segment_file = find_purple_file(purple_dir, "segment.tsv", group_id, sample_id=sample_id)
    if segment_file and os.path.exists(segment_file):
        try:
            with open(segment_file, mode='r', encoding='utf-8') as sf:
                reader = csv.DictReader(sf, delimiter='\t')
                for row in reader:
                    chrom = row.get('chromosome', '').strip()
                    chrom_norm = chrom[3:] if chrom.startswith('chr') else chrom
                    start = int(row.get('start', 0))
                    end = int(row.get('end', 0))
                    germ_status = row.get('germlineStatus', '').strip()
                    baf_count = int(row.get('bafCount', 0))
                    obs_baf = float(row.get('observedBAF', 0.5)) if row.get('observedBAF') else 0.5

                    # Check overlap with any target gene
                    for gene_name in results:
                        bounds = GENE_BOUNDS.get("GRCh38_hmf", {}).get(gene_name)
                        if bounds:
                            g_chrom = bounds['chrom']
                            g_chrom_norm = g_chrom[3:] if g_chrom.startswith('chr') else g_chrom
                            if chrom_norm == g_chrom_norm and start <= bounds['end'] and end >= bounds['start']:
                                if germ_status == 'HET_DELETION' or (baf_count > 0 and (obs_baf > 0.65 or obs_baf < 0.35)):
                                    results[gene_name]['LOH'] = 'LOH'
        except Exception as e:
            print(f"Warning: Failed to parse segment file for LOH fallback {segment_file}: {e}", file=sys.stderr)
        
    return results

def parse_germline_deletions(purple_dir, group_id, genome="GRCh38_hmf", sample_id=None):
    """Parse purple.germline.deletion.tsv or germline_amp_del.tsv to find genes of interest with germline copy number deletions.
    Only includes high-quality calls (filter == 'PASS') that are strictly sub-genic and have flanking diploid/non-deleted support.
    """
    if genome and (genome.startswith('T') or genome.startswith('N') or not any(g in str(genome).lower() for g in ['grch', 'hg', 'b37'])):
        sample_id, genome = genome, "GRCh38_hmf"
    genes_with_del = set()
    deletion_file = find_purple_file(purple_dir, "germline.deletion.tsv", group_id, sample_id=sample_id)
    if not deletion_file:
        deletion_file = find_purple_file(purple_dir, "germline_amp_del.tsv", group_id, sample_id=sample_id)
        
    if not deletion_file or not os.path.exists(deletion_file):
        return genes_with_del
        
    # Helpers
    def normalize_chrom(c):
        return c[3:] if c.startswith('chr') else c

    # Read deletion regions
    regions_by_gene = {}
    try:
        with open(deletion_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row.get('gene', '').strip()
                filt = row.get('filter', '').strip()
                event_type = row.get('type', row.get('driver', '')).upper()
                is_del = 'DEL' in event_type or 'DELETION' in event_type
                if gene in GENES_OF_INTEREST and (filt in ('PASS', '') or filt is None) and is_del:
                    chrom = row.get('chromosome', '').strip()
                    start = int(row.get('regionStart', row.get('start', 0)))
                    end = int(row.get('regionEnd', row.get('end', 0)))
                    if gene not in regions_by_gene:
                        regions_by_gene[gene] = []
                    regions_by_gene[gene].append((chrom, start, end))
    except Exception as e:
        print(f"Warning: Failed to parse germline deletions file {deletion_file}: {e}", file=sys.stderr)
        return genes_with_del

    # Read segment data once if we have candidates to check
    segments_by_chrom = {}
    if regions_by_gene:
        segment_file = find_purple_file(purple_dir, "segment.tsv", group_id, sample_id=sample_id)
        if segment_file and os.path.exists(segment_file):
            try:
                with open(segment_file, mode='r', encoding='utf-8') as sf:
                    reader = csv.DictReader(sf, delimiter='\t')
                    for row in reader:
                        chrom_norm = normalize_chrom(row.get('chromosome', '').strip())
                        if chrom_norm not in segments_by_chrom:
                            segments_by_chrom[chrom_norm] = []
                        segments_by_chrom[chrom_norm].append({
                            'start': int(row.get('start', 0)),
                            'end': int(row.get('end', 0)),
                            'germlineStatus': row.get('germlineStatus', '').strip()
                        })
                for c in segments_by_chrom:
                    segments_by_chrom[c].sort(key=lambda x: x['start'])
            except Exception as e:
                print(f"Warning: Failed to parse segments file {segment_file} for germline check: {e}", file=sys.stderr)

    # Evaluate each gene's deletion calls
    for gene, regions in regions_by_gene.items():
        if not regions:
            continue
            
        # 1. Group contiguous deletion regions on the same chromosome
        regions_by_chrom = {}
        for chrom, start, end in regions:
            chrom_norm = normalize_chrom(chrom)
            if chrom_norm not in regions_by_chrom:
                regions_by_chrom[chrom_norm] = []
            regions_by_chrom[chrom_norm].append((start, end))
            
        # Merge overlapping/adjacent regions on each chromosome
        merged_regions = []
        for chrom_norm, coords in regions_by_chrom.items():
            coords.sort(key=lambda x: x[0])
            cur_start, cur_end = coords[0]
            for next_start, next_end in coords[1:]:
                if next_start <= cur_end + 1:
                    cur_end = max(cur_end, next_end)
                else:
                    merged_regions.append((chrom_norm, cur_start, cur_end))
            merged_regions.append((chrom_norm, cur_start, cur_end))

        # Evaluate the merged regions
        gene_size = 100000  # Default 100 kb fallback
        bounds = GENE_BOUNDS.get(genome, {}).get(gene)
        if bounds:
            gene_size = bounds["end"] - bounds["start"]

        for chrom_norm, m_start, m_end in merged_regions:
            merged_span = m_end - m_start
            
            # Check 1: Size must be strictly sub-genic (< 50% of gene size and <= 20 kb)
            if merged_span >= 0.5 * gene_size or merged_span > 20000:
                continue
                
            # Check 2: Flanking support in segments (both sides must be DIPLOID, not UNKNOWN or deleted)
            has_flanking_support = True
            segs = segments_by_chrom.get(chrom_norm, [])
            if segs:
                upstream_seg = None
                downstream_seg = None
                for seg in segs:
                    if seg['end'] < m_start:
                        if upstream_seg is None or seg['end'] > upstream_seg['end']:
                            upstream_seg = seg
                    if seg['start'] > m_end:
                        if downstream_seg is None or seg['start'] < downstream_seg['start']:
                            downstream_seg = seg
                            
                if not upstream_seg or upstream_seg['germlineStatus'] != 'DIPLOID':
                    has_flanking_support = False
                if not downstream_seg or downstream_seg['germlineStatus'] != 'DIPLOID':
                    has_flanking_support = False
            else:
                has_flanking_support = False
                    
            if has_flanking_support:
                genes_with_del.add(gene)
                break  # If any merged segment passes, the gene is deleted

    return genes_with_del

def parse_germline_catalog(purple_dir, group_id, sample_id=None):
    """Parse purple.driver.catalog.germline.tsv for germline driver events.
    """
    genes_with_driver = set()
    catalog_file = find_purple_file(purple_dir, "driver.catalog.germline.tsv", group_id, sample_id=sample_id)
    if not catalog_file or not os.path.exists(catalog_file):
        return genes_with_driver
        
    try:
        with open(catalog_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row.get('gene', '').strip()
                if gene in GENES_OF_INTEREST:
                    genes_with_driver.add(gene)
    except Exception as e:
        print(f"Warning: Failed to parse germline driver catalog {catalog_file}: {e}", file=sys.stderr)
        
    return genes_with_driver

def load_exon_coordinates(gtf_path, genes_of_interest):
    """Parse a GTF file to get all exon coordinates for genes of interest.
    Returns a dict: {gene: [(start, end), ...]}
    """
    exons = {gene: [] for gene in genes_of_interest}
    if not gtf_path or not os.path.exists(gtf_path):
        return exons
        
    import gzip
    open_fn = gzip.open if gtf_path.endswith('.gz') else open
    try:
        with open_fn(gtf_path, mode='rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'exon':
                    continue
                
                info_str = parts[8]
                gene_name = None
                
                # Fast attribute parsing
                if 'gene_name "' in info_str:
                    gene_name = info_str.split('gene_name "')[1].split('"')[0]
                elif 'gene_symbol "' in info_str:
                    gene_name = info_str.split('gene_symbol "')[1].split('"')[0]
                elif 'gene_id "' in info_str:
                    gene_name = info_str.split('gene_id "')[1].split('"')[0]
                        
                if gene_name in genes_of_interest:
                    start = int(parts[3])
                    end = int(parts[4])
                    exons[gene_name].append((start, end))
                    
        # Sort exon boundaries
        for gene in exons:
            exons[gene] = sorted(list(set(exons[gene])))
    except Exception as e:
        print(f"Warning: Failed to parse GTF {gtf_path}: {e}", file=sys.stderr)
    return exons

def parse_pathogenic_svs(vcf_file, exons_dict, genome):
    """Parse somatic or germline SV VCF to find events disrupting coding exons of genes of interest.
    Returns a dict mapping gene to a list of SV coordinate/type strings.
    """
    disruptions = {gene: [] for gene in exons_dict}
    if not os.path.exists(vcf_file):
        return disruptions
        
    # Get standard chrom bounds for verification
    bounds_subset = GENE_BOUNDS.get(genome, GENE_BOUNDS['GRCh38_hmf'])
        
    try:
        import gzip
        open_fn = gzip.open if vcf_file.endswith('.gz') else open
        with open_fn(vcf_file, mode='rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                    
                filt = parts[6].strip()
                if filt != 'PASS':
                    continue
                    
                chrom = parts[0].strip()
                pos = int(parts[1])
                
                info_str = parts[7]
                info = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info[k] = v
                    else:
                        info[item] = True
                        
                sv_type = info.get('SVTYPE', 'BND')
                end_pos = pos
                if 'END' in info:
                    try:
                        end_pos = int(info['END'])
                    except ValueError:
                        pass
                        
                for gene, exons in exons_dict.items():
                    gene_chrom = bounds_subset.get(gene, {}).get('chrom')
                    if not gene_chrom:
                        continue
                        
                    # Standardize chrom string comparison
                    c_std = chrom.lower() if chrom.lower().startswith('chr') else 'chr' + chrom.lower()
                    gc_std = gene_chrom.lower() if gene_chrom.lower().startswith('chr') else 'chr' + gene_chrom.lower()
                    if c_std != gc_std:
                        continue
                        
                    overlaps = False
                    for (ex_start, ex_end) in exons:
                        if sv_type == 'BND':
                            if ex_start <= pos <= ex_end:
                                overlaps = True
                                break
                        else:
                            if max(pos, ex_start) <= min(end_pos, ex_end):
                                overlaps = True
                                break
                                
                    if overlaps:
                        desc = f"{sv_type} disrupting exon at {chrom}:{pos}"
                        if sv_type != 'BND' and end_pos != pos:
                            desc += f"-{end_pos}"
                        disruptions[gene].append(desc)
    except Exception as e:
        print(f"Warning: Failed to parse SV VCF {vcf_file}: {e}", file=sys.stderr)
    return disruptions

def parse_pave_vcf(pave_dir, group_id, vcf_type, sample_id=None):
    """Parse PAVE somatic or germline VCF to find genes of interest with pathogenic or likely pathogenic variants.
    """
    mutated_genes = set()
    vcf_file = find_pave_vcf(pave_dir, vcf_type, group_id, sample_id=sample_id)
    if not vcf_file or not os.path.exists(vcf_file):
        return mutated_genes
        
    high_impact_effects = {
        'frameshift_variant', 'stop_gained', 'splice_donor_variant', 
        'splice_acceptor_variant', 'start_lost', 'stop_lost', 'transcript_ablation'
    }
    
    try:
        with gzip.open(vcf_file, mode='rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                # Check filter is PASS
                filt = parts[6].strip()
                if filt != 'PASS':
                    continue
                    
                info_str = parts[7]
                info = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info[k] = v
                    else:
                        info[item] = True
                
                # 1. Determine Gene Name
                gene = None
                if 'IMPACT' in info:
                    imp_parts = info['IMPACT'].split(',')
                    if imp_parts:
                        gene = imp_parts[0].strip()
                if not gene and 'PAVE_TI' in info:
                    transcripts = info['PAVE_TI'].split(',')
                    for tx in transcripts:
                        tx_parts = tx.split('|')
                        if len(tx_parts) > 1:
                            gene = tx_parts[1].strip()
                            if gene:
                                break
                                
                if not gene or gene not in GENES_OF_INTEREST:
                    continue
                
                # 2. Check clinical significance (CLNSIG)
                is_pathogenic = False
                if 'CLNSIG' in info:
                    clnsig_val = info['CLNSIG'].lower()
                    if 'pathogenic' in clnsig_val and 'conflicting' not in clnsig_val:
                        is_pathogenic = True
                
                # 3. Check for high-impact mutation consequences (nonsense/frameshift/splice)
                # If the mutation is high impact, it is also considered pathogenic
                is_high_impact = False
                if not is_pathogenic:
                    effects = set()
                    # Parse effects from IMPACT
                    if 'IMPACT' in info:
                        imp_parts = info['IMPACT'].split(',')
                        if len(imp_parts) > 2:
                            effects.add(imp_parts[2].strip().lower())
                    # Parse effects from PAVE_TI
                    if 'PAVE_TI' in info:
                        transcripts = info['PAVE_TI'].split(',')
                        for tx in transcripts:
                            tx_parts = tx.split('|')
                            if len(tx_parts) > 3:
                                effects.add(tx_parts[3].strip().lower())
                    if effects & high_impact_effects:
                        is_high_impact = True
                        
                if is_pathogenic or is_high_impact:
                    mutated_genes.add(gene)
                    
    except Exception as e:
        print(f"Warning: Failed to parse PAVE VCF file {vcf_file}: {e}", file=sys.stderr)
        
    return mutated_genes

def parse_somatic_catalog(purple_dir, group_id, sample_id=None):
    """Parse driver.catalog.somatic.tsv for somatic mutation events.
    """
    catalog_file = find_purple_file(purple_dir, "driver.catalog.somatic.tsv", group_id, sample_id=sample_id)
    mutated_genes = set()
    
    if not catalog_file or not os.path.exists(catalog_file):
        return mutated_genes
        
    try:
        with open(catalog_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_name = row.get('gene', '').strip()
                if gene_name in GENES_OF_INTEREST:
                    # Check mutation counts
                    try:
                        missense = int(row.get('missense', 0))
                        nonsense = int(row.get('nonsense', 0))
                        splice = int(row.get('splice', 0))
                        inframe = int(row.get('inframe', 0))
                        frameshift = int(row.get('frameshift', 0))
                        
                        total_muts = missense + nonsense + splice + inframe + frameshift
                        driver_type = row.get('driver', '').strip()
                        
                        if total_muts > 0 or (driver_type and driver_type not in ['', '-', 'NONE']):
                            mutated_genes.add(gene_name)
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Warning: Failed to parse catalog file {catalog_file}: {e}", file=sys.stderr)
        
    return mutated_genes

def find_methylation_data(explicit_path, output_dir, group_id):
    """Find methylation matrix either from explicit argument or local directory structure."""
    if explicit_path and os.path.exists(explicit_path):
        return explicit_path
        
    search_dirs = [
        os.path.join(output_dir, group_id, 'EMSeq', 'MethylDackel'),
        os.path.join(output_dir, 'EMSeq', 'MethylDackel'),
        os.path.join(os.path.dirname(output_dir.rstrip('/')), 'EMSeq', 'MethylDackel'),
        os.path.join(output_dir, group_id, 'MethylDackel'),
        os.path.join(output_dir, 'MethylDackel')
    ]
    for d in search_dirs:
        if os.path.isdir(d):
            for cand_name in ['MultiMMR_CpG_methylation_matrix.RData', '2026-10-13_MultiMMR_CpG_methylation_matrix.RData']:
                cand = os.path.join(d, cand_name)
                if os.path.exists(cand):
                    return cand
            import glob
            rdata_matches = glob.glob(os.path.join(d, "*methylation_matrix.RData"))
            if rdata_matches:
                return rdata_matches[0]
            tsv_matches = glob.glob(os.path.join(d, "*promoter_cpgs.tsv"))
            if tsv_matches:
                return tsv_matches[0]
    return None

from functools import lru_cache

@lru_cache(maxsize=128)
def parse_methylation(methylation_path, group_id, subject_id, samplesheet_path, genome="GRCh38_hmf"):
    """Parse external methylation file to check if genes are methylated.
    Supports:
    1. Gene-level summary TSV/CSV (old format).
    2. CpG-level site TSV/CSV matrices.
    3. CpG-level RData/rda matrices (converts via temporary Rscript).
    Returns a dictionary of results: {gene: {'status': 'somatic'|'germline', 'tumor_beta': float, 'normal_beta': float}}
    """
    methylation_results = {}
    if not methylation_path or not os.path.exists(methylation_path):
        return methylation_results
        
    import tempfile
    import subprocess
    
    temp_tsv = None
    parsing_target = methylation_path
    
    if methylation_path.lower().endswith(('.rdata', '.rda')):
        cached_tsv = methylation_path + ".promoter_cpgs.tsv"
        if os.path.exists(cached_tsv) and os.path.getsize(cached_tsv) > 0:
            parsing_target = cached_tsv
        else:
            try:
                # Create temporary file for R output
                fd_out, temp_tsv = tempfile.mkstemp(suffix='.tsv')
                os.close(fd_out)
                
                # Create temporary file for promoter coords
                fd_prom, temp_prom = tempfile.mkstemp(suffix='.txt')
                os.close(fd_prom)
                
                # Write full gene region coordinates (promoter + gene body)
                with open(temp_prom, 'w', encoding='utf-8') as f_prom:
                    bounds_subset = GENE_BOUNDS.get(genome, GENE_BOUNDS['GRCh38_hmf'])
                    for gene in GENES_OF_INTEREST:
                        bounds = bounds_subset.get(gene)
                        if not bounds:
                            continue
                        g_start = bounds['start']
                        g_end = bounds['end']
                        g_chrom = bounds['chrom']
                        p_start = min(g_start, g_end) - 2000
                        p_end = max(g_start, g_end) + 2000
                        f_prom.write(f"{g_chrom}\t{p_start}\t{p_end}\t{gene}\n")
                
                # Write R script
                fd_r, r_script = tempfile.mkstemp(suffix='.R')
                os.close(fd_r)
                
                r_code = f"""
args <- commandArgs(trailingOnly = TRUE)
rdata_file <- args[1]
promoter_file <- args[2]
output_file <- args[3]

# Load promoter coordinates
promoters <- read.table(promoter_file, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "start", "end", "gene"))

# Load RData
load(rdata_file)

# Find matrix/data.frame
all_objs <- ls()
matrix_obj <- NULL
for (obj in all_objs) {{
  if (obj %in% c("rdata_file", "promoter_file", "output_file", "promoters", "all_objs", "obj")) next
  val <- get(obj)
  if (is.data.frame(val) || is.matrix(val)) {{
    matrix_obj <- val
    break
  }}
}}

if (is.null(matrix_obj)) {{
  stop("Error: No data.frame or matrix found in RData file")
}}

df <- as.data.frame(matrix_obj)
col_names <- colnames(df)
seqnames_col <- grep("seqnames|chrom|chr|chromosome", col_names, ignore.case = TRUE, value = TRUE)[1]
start_col <- grep("start|pos|position", col_names, ignore.case = TRUE, value = TRUE)[1]

if (is.na(seqnames_col) || is.na(start_col)) {{
  stop("Error: Coordinate columns (seqnames, start) not found in matrix")
}}

standardize_chr <- function(c) {{
  c <- as.character(c)
  c <- ifelse(grepl("^chr", c, ignore.case=TRUE), tolower(c), paste0("chr", tolower(c)))
  return(c)
}}

df_chr_std <- standardize_chr(df[[seqnames_col]])
promoters$chrom_std <- standardize_chr(promoters$chrom)

matching_rows <- c()
for (i in 1:nrow(promoters)) {{
  p_chrom <- promoters$chrom_std[i]
  p_start <- promoters$start[i]
  p_end <- promoters$end[i]
  
  rows_in_range <- which(df_chr_std == p_chrom & df[[start_col]] >= p_start & df[[start_col]] <= p_end)
  if (length(rows_in_range) > 0) {{
    df[rows_in_range, "gene_annotation"] <- promoters$gene[i]
    matching_rows <- c(matching_rows, rows_in_range)
  }}
}}

if (length(matching_rows) > 0) {{
  filtered_df <- df[unique(matching_rows), ]
  write.table(filtered_df, file = output_file, sep = "\\t", quote = FALSE, row.names = FALSE)
}} else {{
  write.table(df[0, ], file = output_file, sep = "\\t", quote = FALSE, row.names = FALSE)
}}
"""
                with open(r_script, 'w', encoding='utf-8') as f_r:
                    f_r.write(r_code)
                    
                # Run Rscript with module loading fallback for cluster environments
                full_r_cmd = (
                    'for init in "$MODULESHOME/init/bash" "/usr/share/lmod/lmod/init/bash" "/usr/share/Modules/init/bash" "/etc/profile.d/modules.sh" "/etc/profile.d/z00_lmod.sh" "/etc/profile.d/z01_lmod.sh"; do '
                    '  if [ -f "$init" ]; then source "$init"; break; fi; '
                    'done; '
                    'if ! command -v Rscript &>/dev/null; then module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R 2>/dev/null || true; fi; '
                    f'Rscript {shlex.quote(r_script)} {shlex.quote(methylation_path)} {shlex.quote(temp_prom)} {shlex.quote(temp_tsv)}'
                )
                res = subprocess.run(["bash", "-c", full_r_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if res.returncode != 0:
                    err_msg = res.stderr.strip() or res.stdout.strip() or "Unknown Rscript error"
                    raise Exception(err_msg)
                
                # Copy temp_tsv to cached_tsv for subsequent threads/runs
                try:
                    import shutil
                    shutil.copyfile(temp_tsv, cached_tsv)
                except Exception:
                    pass

                parsing_target = temp_tsv
                os.remove(temp_prom)
                os.remove(r_script)
                
            except Exception as e:
                print(f"Warning: Failed to extract RData using Rscript: {e}", file=sys.stderr)
                if temp_tsv and os.path.exists(temp_tsv):
                    os.remove(temp_tsv)
                return {}
            
    try:
        with open(parsing_target, mode='r', encoding='utf-8') as f:
            # Determine separator (TSV or CSV)
            sample_line = f.readline()
            delim = ',' if ',' in sample_line else '\t'
            f.seek(0)
            
            reader = csv.DictReader(f, delimiter=delim)
            
            # Detect coordinate columns
            chrom_col = None
            pos_col = None
            if reader.fieldnames:
                for col in reader.fieldnames:
                    col_lower = col.lower()
                    if any(x in col_lower for x in ['seqname', 'chrom', 'chr', 'chromosome']):
                        chrom_col = col
                        break
                for col in reader.fieldnames:
                    col_lower = col.lower()
                    if any(x in col_lower for x in ['start', 'pos', 'position']) and col != chrom_col:
                        pos_col = col
                        break
                        
            is_cpg = chrom_col is not None and pos_col is not None
            
            if is_cpg:
                # Retrieve actual sample IDs from samplesheet
                tumor_sample_id = get_tumor_sample_id(samplesheet_path, group_id)
                normal_sample_id = get_normal_sample_id(samplesheet_path, group_id)
                
                tumor_col = None
                normal_col = None
                
                print(f"DEBUG: Parsing CpG matrix methylation columns for sample {group_id}...", file=sys.stderr)
                if reader.fieldnames:
                    print(f"DEBUG: Found {len(reader.fieldnames)} columns in CpG matrix. First 15: {reader.fieldnames[:15]}", file=sys.stderr)
                
                # Match columns flexibly
                for col in reader.fieldnames:
                    if tumor_sample_id and is_flexible_match(col, tumor_sample_id, 'tumor'):
                        tumor_col = col
                    if normal_sample_id and is_flexible_match(col, normal_sample_id, 'normal'):
                        normal_col = col
                        
                # Fallback to direct group_id matching
                if not tumor_col and group_id:
                    for col in reader.fieldnames:
                        if group_id.lower() in col.lower() or col.lower() in group_id.lower():
                            tumor_col = col
                            break
                            
                # Fallback for normal (e.g. matching subject_id and -N/normal)
                if not normal_col and subject_id:
                    for col in reader.fieldnames:
                        col_l = col.lower()
                        if subject_id.lower() in col_l and (col_l.endswith('-n') or 'normal' in col_l or 'control' in col_l or col_l.endswith('n')):
                            normal_col = col
                            break
                            
                print(f"DEBUG: Matched columns -> tumor_col: {tumor_col}, normal_col: {normal_col}", file=sys.stderr)
                            
                # Prepare gene promoter ranges
                def std_chr(c):
                    c = c.strip().lower()
                    return c if c.startswith('chr') else 'chr' + c
                    
                promoters = {}
                bounds_subset = GENE_BOUNDS.get(genome, GENE_BOUNDS['GRCh38_hmf'])
                for gene in GENES_OF_INTEREST:
                    bounds = bounds_subset.get(gene)
                    if not bounds:
                        continue
                    strand = bounds['strand']
                    g_start = bounds['start']
                    g_end = bounds['end']
                    g_chrom = std_chr(bounds['chrom'])
                    if strand == 1:
                        p_start = g_start - 2000
                        p_end = g_start + 500
                    else:
                        p_start = g_end - 500
                        p_end = g_end + 2000
                    promoters[gene] = {
                        'chrom': g_chrom,
                        'p_start': p_start,
                        'p_end': p_end,
                        'full_start': min(g_start, g_end) - 2000,
                        'full_end': max(g_start, g_end) + 2000
                    }
                    
                gene_betas = {gene: {'tumor': [], 'normal': [], 'cpgs': []} for gene in GENES_OF_INTEREST}
                
                for row in reader:
                    chrom_raw = row.get(chrom_col)
                    pos_raw = row.get(pos_col)
                    if not chrom_raw or not pos_raw:
                        continue
                    try:
                        chrom = std_chr(chrom_raw)
                        pos = int(pos_raw)
                    except ValueError:
                        continue
                        
                    for gene, info in promoters.items():
                        if chrom == info['chrom'] and info['full_start'] <= pos <= info['full_end']:
                            is_promoter = info['p_start'] <= pos <= info['p_end']
                            region_name = 'promoter' if is_promoter else 'body'
                            
                            t_val = None
                            n_val = None
                            if tumor_col:
                                t_val_raw = row.get(tumor_col)
                                if t_val_raw is not None:
                                    try:
                                        t_val = float(t_val_raw)
                                        if is_promoter:
                                            gene_betas[gene]['tumor'].append(t_val)
                                    except ValueError:
                                        pass
                            if normal_col:
                                n_val_raw = row.get(normal_col)
                                if n_val_raw is not None:
                                    try:
                                        n_val = float(n_val_raw)
                                        if is_promoter:
                                            gene_betas[gene]['normal'].append(n_val)
                                    except ValueError:
                                        pass
                            if t_val is not None or n_val is not None:
                                gene_betas[gene]['cpgs'].append({
                                    'pos': pos,
                                    'tumor': t_val if t_val is not None else 0.0,
                                    'normal': n_val,
                                    'region': region_name
                                })
                                        
                # Compute averages and classify
                for gene in GENES_OF_INTEREST:
                    t_list = gene_betas[gene]['tumor']
                    n_list = gene_betas[gene]['normal']
                    if not t_list and not n_list:
                        continue
                        
                    t_avg = sum(t_list) / len(t_list) if t_list else 0.0
                    n_avg = sum(n_list) / len(n_list) if n_list else 0.0
                    
                    status = 'unmethylated'
                    if t_avg >= 0.20:
                        if n_list and n_avg >= 0.20:
                            status = 'germline'
                        elif n_list:
                            if n_avg < 0.1 and (t_avg - n_avg) >= 0.15:
                                status = 'somatic'
                            else:
                                status = 'unmethylated'
                        else:
                            status = 'somatic'
                    elif n_list and n_avg >= 0.20:
                        status = 'germline'
                        
                    # Always include parsed CpG details if data is present
                    methylation_results[gene] = {
                        'status': status,
                        'tumor_beta': t_avg,
                        'normal_beta': n_avg if n_list else None,
                        'cpgs': gene_betas[gene]['cpgs']
                    }
            else:
                # Fallback to old gene-level format
                f.seek(0)
                reader = csv.DictReader(f, delimiter=delim)
                for row in reader:
                    m_sample = row.get('Sample', '').strip()
                    m_patient = row.get('Patient', '').strip()
                    m_type = row.get('Type', '').strip().lower()
                    
                    match = False
                    if m_sample and group_id and m_sample.lower() == group_id.lower():
                        match = True
                    elif m_patient and subject_id and m_patient.lower() == subject_id.lower() and m_type in ['tumour', 'tumor']:
                        match = True
                        
                    if match:
                        for gene in GENES_OF_INTEREST:
                            val = row.get(gene)
                            if val is not None:
                                try:
                                    if float(val) > 0.1:
                                        methylation_results[gene] = {
                                            'status': 'somatic',
                                            'tumor_beta': float(val),
                                            'normal_beta': None
                                        }
                                except ValueError:
                                    pass
                        break
                        
    except Exception as e:
        print(f"Warning: Failed to parse methylation file {methylation_path}: {e}", file=sys.stderr)
        
    # Clean up temp output file
    if temp_tsv and os.path.exists(temp_tsv):
        os.remove(temp_tsv)
        
    return methylation_results

def parse_vcf_variants(vcf_file, target_gene, var_type, group_id, min_somatic_vaf=0.0, tumor_sample_id=None, normal_sample_id=None):
    variants = []
    if not vcf_file or not os.path.exists(vcf_file):
        return variants
        
    tumor_col_idx = 9  # default fallback
    normal_col_idx = None
    try:
        with gzip.open(vcf_file, mode='rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    cols = line.strip().split('\t')
                    if tumor_sample_id and tumor_sample_id in cols:
                        tumor_col_idx = cols.index(tumor_sample_id)
                    elif group_id in cols:
                        tumor_col_idx = cols.index(group_id)
                    else:
                        tumor_col_idx = 10 if len(cols) > 10 else 9
                        
                    if normal_sample_id and normal_sample_id in cols:
                        normal_col_idx = cols.index(normal_sample_id)
                    elif len(cols) > 9 and tumor_col_idx == 10:
                        normal_col_idx = 9
                    elif len(cols) > 10 and tumor_col_idx == 9:
                        normal_col_idx = 10
                    continue
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                filt = parts[6].strip()
                if filt != 'PASS':
                    continue
                    
                pos = int(parts[1])
                info_str = parts[7]
                info = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info[k] = v
                    else:
                        info[item] = True
                        
                # Determine gene
                gene = None
                if 'IMPACT' in info:
                    imp_parts = info['IMPACT'].split(',')
                    if imp_parts:
                        gene = imp_parts[0].strip()
                if not gene and 'PAVE_TI' in info:
                    transcripts = info['PAVE_TI'].split(',')
                    for tx in transcripts:
                        tx_parts = tx.split('|')
                        if len(tx_parts) > 1:
                            gene = tx_parts[1].strip()
                            if gene:
                                break
                                
                if gene != target_gene:
                    continue

                # Pathogenicity and Driver Filters
                if var_type == "somatic" and 'REPORTED' not in info:
                    continue
                    
                if var_type == "germline":
                    is_pathogenic = False
                    clnsig = info.get('CLNSIG', '')
                    clnsig_lower = clnsig.lower()
                    if 'pathogenic' in clnsig_lower or 'likely_pathogenic' in clnsig_lower:
                        if 'conflicting' not in clnsig_lower and 'benign' not in clnsig_lower:
                            is_pathogenic = True
                            
                    is_truncating = False
                    high_impact_consequences = [
                        'frameshift_variant', 'stop_gained', 
                        'splice_donor_variant', 'splice_acceptor_variant',
                        'start_lost', 'stop_lost', 'transcript_ablation'
                    ]
                    if 'IMPACT' in info:
                        imp_parts = info['IMPACT'].split(',')
                        if len(imp_parts) > 2:
                            consequence = imp_parts[2].strip()
                            if any(c in consequence for c in high_impact_consequences):
                                is_truncating = True
                    
                    pave_ti = info.get('PAVE_TI', '')
                    if not is_truncating and pave_ti:
                        # Fallback: check any transcript impact for this variant
                        for tx in pave_ti.split(','):
                            tx_parts = tx.split('|')
                            if len(tx_parts) > 3:
                                consequence = tx_parts[3].strip()
                                if any(c in consequence for c in high_impact_consequences):
                                    is_truncating = True
                                    break
                                        
                    if not is_pathogenic and not is_truncating and 'REPORTED' not in info:
                        continue
                    
                # Extract canonical transcript change (HGVS.p or HGVS.c)
                change = ""
                if 'PAVE_TI' in info:
                    transcripts = info['PAVE_TI'].split(',')
                    # Try canonical first
                    for tx in transcripts:
                        tx_parts = tx.split('|')
                        if len(tx_parts) > 6:
                            is_canonical = tx_parts[4].strip() in ['true', '1']
                            if is_canonical:
                                change = tx_parts[6].strip() or tx_parts[5].strip()
                                if change:
                                    break
                    # Try any if canonical has no change
                    if not change:
                        for tx in transcripts:
                            tx_parts = tx.split('|')
                            if len(tx_parts) > 6:
                                change = tx_parts[6].strip() or tx_parts[5].strip()
                                if change:
                                    break
                
                # Extract VAF from the tumor genotype field (between 0.0 and 1.0)
                vaf = None
                format_fields = parts[8].split(':')
                if 'AF' in format_fields and len(parts) > tumor_col_idx:
                    af_idx = format_fields.index('AF')
                    tumor_vals = parts[tumor_col_idx].split(':')
                    if len(tumor_vals) > af_idx:
                        af_val = tumor_vals[af_idx]
                        if ',' in af_val:
                            af_val = af_val.split(',')[0]
                        try:
                            vaf = float(af_val)
                        except ValueError:
                            pass
                            
                # Also extract normal VAF if normal column is known
                norm_vaf = None
                if normal_col_idx is not None and 'AF' in format_fields and len(parts) > normal_col_idx:
                    af_idx = format_fields.index('AF')
                    norm_vals = parts[normal_col_idx].split(':')
                    if len(norm_vals) > af_idx:
                        af_val = norm_vals[af_idx]
                        if ',' in af_val:
                            af_val = af_val.split(',')[0]
                        try:
                            norm_vaf = float(af_val)
                        except ValueError:
                            pass
                
                cn = 2.0
                if 'PURPLE_CN' in info:
                    try:
                        cn = float(info['PURPLE_CN'])
                    except ValueError:
                        pass
                        
                macn = 1.0
                if 'PURPLE_MACN' in info:
                    try:
                        macn = float(info['PURPLE_MACN'])
                    except ValueError:
                        pass
                        
                vcn = 0.0
                if 'PURPLE_VCN' in info:
                    try:
                        vcn = float(info['PURPLE_VCN'])
                    except ValueError:
                        pass
                        
                biallelic = 'BIALLELIC' in info
                
                lps = None
                if 'LPS' in info:
                    try:
                        lps = int(info['LPS'])
                    except ValueError:
                        pass
                        
                lps_rc = None
                if 'LPS_RC' in info:
                    try:
                        lps_rc = int(info['LPS_RC'])
                    except ValueError:
                        pass
                
                # Calculate Phasing Status
                # To prevent placeholders when there is no actual read-backed phasing,
                # somatic variants must have LPS or be biallelic to have a determined phasing status.
                if var_type == "somatic" and lps is None and not biallelic:
                    if macn < 0.2:
                        phasing = "Unphased (LOH region)"
                    elif cn < 1.7:
                        phasing = "Unphased (Loss region)"
                    elif cn > 2.8:
                        phasing = "Unphased (Gain region)"
                    else:
                        phasing = "Unphased"
                else:
                    phasing = "Diploid (no CNV/LOH)"
                    if macn < 0.2:  # LOH region
                        if biallelic or vcn >= 0.8:
                            phasing = "LOH (mutated allele retained, wild-type lost)"
                        else:
                            phasing = "LOH (mutated allele lost, wild-type retained)"
                    elif cn < 1.7:  # Loss
                        if biallelic or vcn >= 0.8:
                            phasing = "In Phase with Loss"
                        else:
                            phasing = "Out of Phase with Loss"
                    elif cn > 2.8:  # Gain
                        if vcn >= 1.5:
                            phasing = "In Phase with Gain (amplified)"
                        else:
                            phasing = "Out of Phase with Gain (non-amplified)"
                        
                if vaf is None:
                    vaf_pct = "NA"
                else:
                    vaf_val = vaf * 100
                    if vaf_val < 0.1:
                        vaf_pct = f"{vaf_val:.3f}%"
                    elif vaf_val < 1.0:
                        vaf_pct = f"{vaf_val:.2f}%"
                    elif vaf_val < 10.0:
                        vaf_pct = f"{vaf_val:.1f}%"
                    else:
                        vaf_pct = f"{int(round(vaf_val))}%"
                if var_type == "somatic" and vaf is not None and vaf < min_somatic_vaf:
                    continue
                    
                if var_type == "germline":
                    label_parts = [f"Tumor VAF: {vaf_pct}"]
                    if norm_vaf is not None:
                        n_vaf_val = norm_vaf * 100
                        if n_vaf_val < 0.1:
                            n_vaf_pct = f"{n_vaf_val:.3f}%"
                        elif n_vaf_val < 1.0:
                            n_vaf_pct = f"{n_vaf_val:.2f}%"
                        elif n_vaf_val < 10.0:
                            n_vaf_pct = f"{n_vaf_val:.1f}%"
                        else:
                            n_vaf_pct = f"{int(round(n_vaf_val))}%"
                        label_parts.append(f"Normal VAF: {n_vaf_pct}")
                    if cn is not None and abs(cn - 2.0) > 0.1:
                        label_parts.append(f"Tumor CN: {cn:.2f}")
                    label = f"Germline {target_gene} {change} ({', '.join(label_parts)})"
                else:
                    label_parts = [f"VAF: {vaf_pct}"]
                    if cn is not None and abs(cn - 2.0) > 0.1:
                        label_parts.append(f"CN: {cn:.2f}")
                    label = f"Somatic {target_gene} {change} ({', '.join(label_parts)})"
                
                variants.append({
                    'pos': pos,
                    'label': label,
                    'change': change,
                    'phasing': phasing,
                    'vaf': vaf,
                    'lps': lps,
                    'lps_rc': lps_rc,
                    'var_type': var_type
                })
    except Exception as e:
        print(f"Warning: Failed to parse variants from {vcf_file}: {e}", file=sys.stderr)
        
    return variants

def get_detailed_variants(purple_dir, pave_dir, group_id, gene_name, subject_id=None, samplesheet_path=None, genome="GRCh38_hmf", methylation_data=None, min_somatic_vaf=0.0, tumor_sample_id=None, normal_sample_id=None):
    """Parse somatic and germline VCFs from PURPLE/PAVE to extract variant details (pos, label, phasing, vaf) for a specific gene.
    """
    variants = []
    
    # 1. Parse somatic VCF from purple dir, falling back to pave dir
    somatic_vcf = find_purple_file(purple_dir, "somatic.vcf.gz", group_id, sample_id=tumor_sample_id)
    if not somatic_vcf:
        somatic_vcf = find_pave_vcf(pave_dir, "somatic", group_id, sample_id=tumor_sample_id)
        
    if somatic_vcf:
        variants.extend(parse_vcf_variants(somatic_vcf, gene_name, "somatic", group_id, min_somatic_vaf, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id))
    
    # 2. Parse germline VCF from purple dir, falling back to pave dir
    germline_vcf = find_purple_file(purple_dir, "germline.vcf.gz", group_id, sample_id=tumor_sample_id)
    if not germline_vcf:
        germline_vcf = find_pave_vcf(pave_dir, "germline", group_id, sample_id=tumor_sample_id)
    if not germline_vcf and normal_sample_id:
        germline_vcf = find_purple_file(purple_dir, "germline.vcf.gz", group_id, sample_id=normal_sample_id)
        if not germline_vcf:
            germline_vcf = find_pave_vcf(pave_dir, "germline", group_id, sample_id=normal_sample_id)
        
    if germline_vcf:
        variants.extend(parse_vcf_variants(germline_vcf, gene_name, "germline", group_id, min_somatic_vaf, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id))
    
    # 3. Add methylation if present
    if methylation_data:
        meth = parse_methylation(methylation_data, group_id, subject_id, samplesheet_path, genome)
        if gene_name in meth and meth[gene_name]['status'] != 'unmethylated':
            status = meth[gene_name]['status']
            t_beta = meth[gene_name]['tumor_beta']
            n_beta = meth[gene_name]['normal_beta']
            
            # Retrieve TSS/promoter coordinate to place the variant box at the start of the gene
            bounds = GENE_BOUNDS.get(genome, GENE_BOUNDS['GRCh38_hmf']).get(gene_name)
            if bounds:
                strand = bounds['strand']
                tss = bounds['start'] if strand == 1 else bounds['end']
                
                label_parts = [f"{status.capitalize()} Hypermethylation"]
                if n_beta is not None:
                    label_parts.append(f"(Beta: T={t_beta:.2f}, N={n_beta:.2f})")
                else:
                    label_parts.append(f"(Beta: {t_beta:.2f})")
                label = " ".join(label_parts)
                
                variants.append({
                    'pos': tss,
                    'label': label,
                    'phasing': 'N/A (epigenetic)',
                    'vaf': t_beta,
                    'lps': None,
                    'lps_rc': None,
                    'var_type': f"{status}_methylation"
                })
    
    # Deduplicate variants by position and label to prevent duplicates
    dedup_vars = []
    seen = set()
    for v in variants:
        key = (v['pos'], v['label'])
        if key not in seen:
            seen.add(key)
            dedup_vars.append(v)
            
    return dedup_vars

def inject_igv_navigation_script(igv_html_path):
    """Injects variant search query parameter scrolling JS into the generated IGV HTML report."""
    if not os.path.exists(igv_html_path):
        return
        
    js_code = """
<!-- INJECTED BY ONCOANALYSER POST-PROCESSING -->
<script>
window.addEventListener('load', function() {
    const urlParams = new URLSearchParams(window.location.search);
    const targetLocus = urlParams.get('locus') || window.location.hash.substring(1);
    if (targetLocus) {
        console.log("IGV Report Injected Script: Found target locus parameter: " + targetLocus);
        // Wait for table and IGV to finish rendering (1.5s delay is safe)
        setTimeout(function() {
            const rows = document.querySelectorAll('tr');
            let matched = false;
            for (let row of rows) {
                // Look for cells matching the coordinate
                if (row.textContent.includes(targetLocus)) {
                    console.log("IGV Report Injected Script: Found matching row for locus, clicking row.");
                    row.click();
                    row.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                // If coordinates didn't match directly, try clean digits matching
                const targetNums = targetLocus.replace(/\\D/g, "");
                if (targetNums.length > 3) {
                    for (let row of rows) {
                        const rowClean = row.textContent.replace(/\\D/g, "");
                        if (rowClean.includes(targetNums)) {
                            console.log("IGV Report Injected Script: Found partial coordinate match, clicking row.");
                            row.click();
                            row.scrollIntoView({ behavior: 'smooth', block: 'center' });
                            break;
                        }
                    }
                }
            }
        }, 1500);
    }
});
</script>
<!-- END INJECTED BY ONCOANALYSER POST-PROCESSING -->
"""
    try:
        with open(igv_html_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        if "<!-- INJECTED BY ONCOANALYSER POST-PROCESSING -->" in content:
            return  # Already injected
            
        # Insert before </body>
        if "</body>" in content:
            new_content = content.replace("</body>", f"{js_code}</body>")
        elif "</html>" in content:
            new_content = content.replace("</html>", f"{js_code}</html>")
        else:
            new_content = content + js_code
            
        with open(igv_html_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Successfully injected navigation script into {igv_html_path}", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Failed to inject navigation script into {igv_html_path}: {e}", file=sys.stderr)

def parse_clinical_metadata(ihc_tsv, ihc_file, group_id, subject_id, cmd_age=None, cmd_tissue=None):
    """
    Parse clinical metadata (IHC results, age, tissue) from various inputs.
    Returns a dict: {'ihc_results': {gene: status}, 'age': int, 'tissue': str}
    """
    ihc_results = {g: "No Data" for g in ['MLH1', 'MSH2', 'MSH6', 'PMS2']}
    age = 999
    tissue = "unknown"
    
    if cmd_age is not None:
        try:
            age = int(cmd_age)
        except (ValueError, TypeError):
            pass
            
    if cmd_tissue is not None:
        tissue = cmd_tissue.strip().lower()

    # Helper to normalize patient/sample IDs for soft matching
    def normalize_id(identifier):
        if not identifier:
            return ""
        import re
        val = identifier.upper().strip()
        val = re.sub(r'-(UHN|OTB|PMH|MOH|MOHCCN|MSS|MSI|MT|N|T)-', '-', val)
        return re.sub(r'[^A-Z0-9]', '', val).lower()

    # 1. Parse sample-level key-value text file if provided
    if ihc_file and os.path.exists(ihc_file):
        try:
            with open(ihc_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or ':' not in line:
                        continue
                    key, val = line.split(':', 1)
                    key = key.strip().upper()
                    val = val.strip()
                    for gene in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
                        if gene in key and "RESULT" in key:
                            ihc_results[gene] = val
                    if "AGE" in key and cmd_age is None:
                        try:
                            age = int(val)
                        except ValueError:
                            pass
                    if ("TISSUE" in key or "ORIGIN" in key or "CANCER" in key) and cmd_tissue is None:
                        tissue = val.lower()
        except Exception as e:
            print(f"Warning: Failed to parse clinical file {ihc_file}: {e}", file=sys.stderr)
            
    # 2. Parse cohort-wide TSV file if provided
    if ihc_tsv and os.path.exists(ihc_tsv):
        try:
            with open(ihc_tsv, 'r', encoding='utf-8') as f:
                first_line = f.readline()
                delim = ',' if ',' in first_line else '\t'
                f.seek(0)
                
                reader = csv.DictReader(f, delimiter=delim)
                sample_col = None
                patient_col = None
                gene_cols = {}
                age_col = None
                tissue_col = None
                
                if reader.fieldnames:
                    for col in reader.fieldnames:
                        col_l = col.lower()
                        if col_l in ['sample', 'sample_id', 'group_id', 'group']:
                            sample_col = col
                        elif col_l in ['patient', 'patient_id', 'subject_id', 'subject']:
                            patient_col = col
                        
                        # Find IHC columns
                        for gene in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
                            if gene in col.upper() and ("IHC" in col.upper() or "RESULT" in col.upper() or col.upper() == gene):
                                gene_cols[gene] = col
                                
                        # Find age/tissue columns
                        if 'age' in col_l:
                            age_col = col
                        if any(x in col_l for x in ['tissue', 'origin', 'cancer']):
                            tissue_col = col
                                
                norm_subject = normalize_id(subject_id)
                norm_group = normalize_id(group_id)
                
                for row in reader:
                    match = False
                    row_sample = row.get(sample_col, '').strip() if sample_col else ''
                    row_patient = row.get(patient_col, '').strip() if patient_col else ''
                    
                    norm_row_sample = normalize_id(row_sample)
                    norm_row_patient = normalize_id(row_patient)
                    
                    if row_sample and group_id and row_sample.lower() == group_id.lower():
                        match = True
                    elif row_patient and subject_id and row_patient.lower() == subject_id.lower():
                        match = True
                    elif row_patient and group_id and group_id.lower().startswith(row_patient.lower()):
                        match = True
                    elif norm_row_patient and norm_subject and norm_row_patient == norm_subject:
                        match = True
                    elif norm_row_patient and norm_group and norm_row_patient == norm_group:
                        match = True
                    elif norm_row_sample and norm_group and norm_row_sample == norm_group:
                        match = True
                        
                    if match:
                        for gene, col in gene_cols.items():
                            val = row.get(col, '').strip()
                            if val:
                                ihc_results[gene] = val
                        if age_col and cmd_age is None:
                            try:
                                age = int(row.get(age_col).strip())
                            except (ValueError, TypeError):
                                pass
                        if tissue_col and cmd_tissue is None:
                            val = row.get(tissue_col, '').strip()
                            if val:
                                tissue = val.lower()
                        break
        except Exception as e:
            print(f"Warning: Failed to parse clinical TSV {ihc_tsv}: {e}", file=sys.stderr)
            
    return {
        'ihc_results': ihc_results,
        'age': age,
        'tissue': tissue
    }

def parse_orthogonal_metadata(orthogonal_tsv, group_id, subject_id):
    """
    Parse orthogonal test / known mutations from a TSV or CSV file.
    Supports both standardized tidy format and raw multi-column spreadsheet formats.
    Returns a dict:
      {
        'pathology_purity': float or None,
        'test_result': str or None,
        'variants': list of dicts
      }
    """
    if not orthogonal_tsv or not os.path.exists(orthogonal_tsv):
        return None
        
    def normalize_id(identifier):
        if not identifier:
            return ""
        import re
        val = identifier.upper().strip()
        val = re.sub(r'-(UHN|OTB|PMH|MOH|MOHCCN|MSS|MSI|MT|N|T)-', '-', val)
        return re.sub(r'[^A-Z0-9]', '', val).lower()

    norm_group = normalize_id(group_id)
    norm_subject = normalize_id(subject_id)
    
    pathology_purity = None
    overall_test_result = None
    variants = []
    
    try:
        with open(orthogonal_tsv, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            delim = ',' if ',' in first_line else '\t'
            f.seek(0)
            
            reader = csv.DictReader(f, delimiter=delim)
            if not reader.fieldnames:
                return None
                
            field_map = {}
            for col in reader.fieldnames:
                col_clean = col.strip().lower()
                if col_clean in ['patient_id', 'patient', 'sample', 'sample_id', 'group_id', 'subject_id']:
                    field_map['patient'] = col
                elif 'sample_type' in col_clean or col_clean == 'type':
                    field_map['sample_type'] = col
                elif any(x in col_clean for x in ['pathology', 'cellularity', 'tumor su']):
                    field_map['pathology_purity'] = col
                elif col_clean in ['gene', 'target_gene']:
                    field_map['gene'] = col
                elif col_clean in ['comments', 'comment']:
                    field_map['comments'] = col
                elif col_clean in ['hgvsc', 'cdna']:
                    field_map['hgvsc'] = col
                elif col_clean in ['hgvsp', 'protein']:
                    field_map['hgvsp'] = col
                elif col_clean in ['variant', 'variants']:
                    field_map['variant'] = col
                elif col_clean in ['origin', 'germline_status', 'variant_type']:
                    field_map['origin'] = col
                elif col_clean in ['clinical_significance', 'significance', 'interpretation']:
                    field_map['clinical_significance'] = col
                elif any(x in col_clean for x in ['test result', 'test_result', 'result']):
                    field_map['test_result'] = col

            var_cols = [col for col in reader.fieldnames if re.match(r'^variant\s*\d+', col.strip(), re.I)]
            
            for row in reader:
                p_val = row.get(field_map.get('patient', ''), '').strip()
                norm_p = normalize_id(p_val)
                
                match = False
                if p_val and group_id and p_val.lower() == group_id.lower():
                    match = True
                elif p_val and subject_id and p_val.lower() == subject_id.lower():
                    match = True
                elif norm_p and norm_group and (norm_p == norm_group or norm_group.startswith(norm_p) or norm_p.startswith(norm_group)):
                    match = True
                elif norm_p and norm_subject and (norm_p == norm_subject or norm_subject.startswith(norm_p) or norm_p.startswith(norm_subject)):
                    match = True
                    
                if not match:
                    continue
                    
                s_type = row.get(field_map.get('sample_type', ''), '').strip().lower()
                h_c_raw = row.get(field_map.get('hgvsc', ''), '').strip()
                h_p_raw = row.get(field_map.get('hgvsp', ''), '').strip()
                v_str_raw = row.get(field_map.get('variant', ''), '').strip()
                
                h_c = h_c_raw if h_c_raw not in ['NA', '-', ''] else ''
                h_p = h_p_raw if h_p_raw not in ['NA', '-', ''] else ''
                v_str = v_str_raw if v_str_raw not in ['NA', '-', ''] else ''
                
                # Check for any valid variant in extra var_cols
                valid_var_cols = []
                for vc in var_cols:
                    val = row.get(vc, '').strip()
                    if val and val not in ['NA', '-']:
                        valid_var_cols.append((vc, val))

                if s_type == 'normal':
                    has_var = bool(valid_var_cols) or bool(h_c) or bool(h_p) or bool(v_str)
                    if not has_var:
                        continue
                        
                if pathology_purity is None and 'pathology_purity' in field_map:
                    raw_pur = row.get(field_map['pathology_purity'], '').strip()
                    if raw_pur and raw_pur not in ['NA', '-']:
                        try:
                            val_f = float(raw_pur)
                            if val_f > 1.0:
                                val_f = val_f / 100.0
                            pathology_purity = val_f
                        except ValueError:
                            pass
                            
                if 'test_result' in field_map:
                    tr = row.get(field_map['test_result'], '').strip()
                    if tr and tr not in ['NA', '-', 'Normal Control']:
                        overall_test_result = tr
                    elif not overall_test_result and tr and tr not in ['NA', '-']:
                        overall_test_result = tr

                gene_val = row.get(field_map.get('gene', ''), '').strip()
                if gene_val in ['NA', '-']:
                    gene_val = ''
                if not gene_val and 'comments' in field_map:
                    comm = row.get(field_map['comments'], '').strip()
                    for g in GENES_OF_INTEREST:
                        if g.lower() == comm.lower() or g in comm.upper():
                            gene_val = g
                            break
                    if not gene_val and 'no mutation' in comm.lower():
                        gene_val = 'None'
                        
                origin_val = row.get(field_map.get('origin', ''), '').strip()
                if origin_val in ['NA', '-']:
                    origin_val = ''
                clin_sig = row.get(field_map.get('clinical_significance', ''), '').strip()
                if clin_sig in ['NA', '-']:
                    clin_sig = ''
                if not clin_sig and overall_test_result:
                    clin_sig = overall_test_result
                if not origin_val and overall_test_result:
                    if 'germline' in overall_test_result.lower():
                        origin_val = 'germline'
                    elif 'pathogenic' in overall_test_result.lower() or 'vus' in overall_test_result.lower():
                        origin_val = 'somatic'

                if h_c or h_p or v_str:
                    if not h_c and v_str:
                        mc = re.search(r'(c\.[^\s,]+)', v_str)
                        if mc: h_c = mc.group(1)
                    if not h_p and v_str:
                        mp = re.search(r'(p\.[^\s,]+|p\.\([^\)]+\))', v_str)
                        if mp: h_p = mp.group(1)
                    full_var = v_str if v_str else f"{h_c} {h_p}".strip()
                    if gene_val != 'None':
                        variants.append({
                            'gene': gene_val or 'MMR',
                            'hgvsc': h_c,
                            'hgvsp': h_p,
                            'variant': full_var,
                            'origin': origin_val or 'unknown',
                            'clinical_significance': clin_sig or 'Unknown',
                            'test_result': overall_test_result or ''
                        })
                elif valid_var_cols:
                    for vc, raw_v in valid_var_cols:
                        mc = re.search(r'(c\.[^\s,]+)', raw_v)
                        mp = re.search(r'(p\.[^\s,]+|p\.\([^\)]+\))', raw_v)
                        if gene_val != 'None':
                            variants.append({
                                'gene': gene_val or 'MMR',
                                'hgvsc': mc.group(1) if mc else '',
                                'hgvsp': mp.group(1) if mp else '',
                                'variant': raw_v,
                                'origin': origin_val or 'unknown',
                                'clinical_significance': clin_sig or 'Unknown',
                                'test_result': overall_test_result or ''
                            })
                            
            if not variants and overall_test_result and 'no pathogenic' in overall_test_result.lower():
                variants.append({
                    'gene': 'None',
                    'hgvsc': '',
                    'hgvsp': '',
                    'variant': 'No mutation detected',
                    'origin': 'NA',
                    'clinical_significance': 'Negative',
                    'test_result': overall_test_result
                })
                
    except Exception as e:
        print(f"Warning: Failed to parse orthogonal test TSV {orthogonal_tsv}: {e}", file=sys.stderr)
        return None

    return {
        'pathology_purity': pathology_purity,
        'test_result': overall_test_result,
        'variants': variants
    }

def evaluate_orthogonal_concordance(gene, orth_metadata, gene_reports):
    """
    Evaluates concordance between orthogonal test findings and oncoanalyser findings for a gene.
    """
    if not orth_metadata:
        return {'has_data': False}
        
    orth_vars = orth_metadata.get('variants', [])
    overall_res = orth_metadata.get('test_result', '')
    
    gene_orth_vars = [v for v in orth_vars if v['gene'].upper() == gene.upper() or v['gene'] == 'MMR']
    
    rep = gene_reports.get(gene, {})
    s_mut = rep.get('s_mut', 'Wildtype')
    g_mut = rep.get('g_mut', 'Wildtype')
    onco_vars = []
    if s_mut and s_mut != 'Wildtype':
        onco_vars.append(f"Somatic: {s_mut}")
    if g_mut and g_mut != 'Wildtype':
        onco_vars.append(f"Germline: {g_mut}")
    onco_str = "; ".join(onco_vars) if onco_vars else "Wildtype"
    
    if gene_orth_vars:
        orth_str_list = [v['variant'] for v in gene_orth_vars]
        orth_str = "; ".join(orth_str_list)
        clin_sigs = list(dict.fromkeys(filter(None, [v.get('clinical_significance') for v in gene_orth_vars])))
        clin_sig_str = ", ".join(clin_sigs) if clin_sigs else "Unknown"
        
        matched = False
        import re
        for ov in gene_orth_vars:
            h_p = ov.get('hgvsp', '')
            h_c = ov.get('hgvsc', '')
            p_clean = re.sub(r'[^A-Za-z0-9]', '', h_p.replace('p.', '')) if h_p else ''
            
            if p_clean and p_clean.lower() in re.sub(r'[^A-Za-z0-9]', '', onco_str).lower():
                matched = True
                break
            m_codon = re.search(r'\d+', h_p) if h_p else None
            if m_codon and m_codon.group(0) in onco_str:
                matched = True
                break
            if h_c and h_c.lower() in onco_str.lower():
                matched = True
                break
                
        if matched:
            return {
                'has_data': True,
                'orthogonal_findings': orth_str,
                'clinical_significance': clin_sig_str,
                'onco_findings': onco_str,
                'concordance': "Concordant (Confirmed)",
                'badge_class': "badge-somatic",
                'is_concordant': True
            }
        else:
            if onco_str == 'Wildtype':
                return {
                    'has_data': True,
                    'orthogonal_findings': orth_str,
                    'clinical_significance': clin_sig_str,
                    'onco_findings': onco_str,
                    'concordance': "Discordant (Not Detected)",
                    'badge_class': "badge-discordant",
                    'is_concordant': False
                }
            else:
                return {
                    'has_data': True,
                    'orthogonal_findings': orth_str,
                    'clinical_significance': clin_sig_str,
                    'onco_findings': onco_str,
                    'concordance': "Partially Concordant",
                    'badge_class': "badge-loh",
                    'is_concordant': False
                }
    else:
        is_negative_patient = any(v['gene'] == 'None' for v in orth_vars) or (overall_res and 'no pathogenic' in overall_res.lower())
        if is_negative_patient:
            if onco_str == 'Wildtype':
                return {
                    'has_data': True,
                    'orthogonal_findings': "No Mutation Reported",
                    'clinical_significance': "Negative",
                    'onco_findings': "Wildtype",
                    'concordance': "Consistent (Wildtype)",
                    'badge_class': "badge-neutral",
                    'is_concordant': True
                }
            else:
                return {
                    'has_data': True,
                    'orthogonal_findings': "No Mutation Reported",
                    'clinical_significance': "Negative",
                    'onco_findings': onco_str,
                    'concordance': "Novel Finding (Not in Orthogonal Assay)",
                    'badge_class': "badge-loh",
                    'is_concordant': True
                }
                
    return {'has_data': False}

def resolve_swgs_ichor_dir(explicit_dir, output_dir, group_id, tumor_sample_id=None):
    """Auto-discover sWGS ichorCNA directory for a sample if not explicitly provided."""
    if explicit_dir and os.path.isdir(explicit_dir):
        if tumor_sample_id:
            for c in [
                os.path.join(explicit_dir, tumor_sample_id),
                os.path.join(explicit_dir, tumor_sample_id, tumor_sample_id),
                os.path.join(explicit_dir, "ichorCNA_t_vs_n", tumor_sample_id),
                os.path.join(explicit_dir, "ichorCNA_plasma_vs_n", tumor_sample_id),
                os.path.join(explicit_dir, "ichorCNA_tumor_only", tumor_sample_id),
                os.path.join(explicit_dir, "ichorCNA_plasma_only", tumor_sample_id),
                os.path.join(explicit_dir, "ichorCNA", tumor_sample_id),
            ]:
                if os.path.isdir(c):
                    return c
        if group_id:
            for c in [
                os.path.join(explicit_dir, group_id, "ichorCNA_t_vs_n"),
                os.path.join(explicit_dir, group_id, "ichorCNA_plasma_vs_n"),
                os.path.join(explicit_dir, group_id, "ichorCNA_tumor_only"),
                os.path.join(explicit_dir, group_id, "ichorCNA"),
                os.path.join(explicit_dir, group_id, "sWGS_ichorCNA"),
                os.path.join(explicit_dir, group_id),
            ]:
                if os.path.isdir(c):
                    return c
        return explicit_dir

    search_bases = [output_dir]
    if group_id:
        search_bases.append(os.path.join(output_dir, group_id))
    parent_dir = os.path.dirname(os.path.abspath(output_dir)) if output_dir else None
    if parent_dir and os.path.isdir(parent_dir):
        search_bases.append(parent_dir)
        if group_id:
            search_bases.append(os.path.join(parent_dir, group_id))

    for base in search_bases:
        if not base or not os.path.isdir(base):
            continue
        for sub in [
            "ichorCNA_t_vs_n",
            "ichorCNA_plasma_vs_n",
            "ichorCNA_tumor_only",
            "ichorCNA_plasma_only",
            "ichorCNA",
            "sWGS_ichorCNA",
            "sWGS/ichorCNA",
            "swgs/ichorCNA",
            "sWGS",
            "swgs",
        ]:
            cand = os.path.join(base, sub)
            if os.path.isdir(cand):
                if tumor_sample_id:
                    cand_sample = os.path.join(cand, tumor_sample_id)
                    cand_nested = os.path.join(cand, tumor_sample_id, tumor_sample_id)
                    if os.path.isdir(cand_sample):
                        return cand_sample
                    if os.path.isdir(cand_nested):
                        return cand_nested
                for f in os.listdir(cand):
                    sub_cand = os.path.join(cand, f)
                    if os.path.isdir(sub_cand) and not f.startswith('.'):
                        for s_file in os.listdir(sub_cand):
                            if s_file.endswith('.params.txt'):
                                return sub_cand
                return cand
    return None

def resolve_panel_ichor_dir(explicit_dir, output_dir, group_id, tumor_sample_id=None):
    """Auto-discover panel off-target ichorCNA directory for a sample."""
    if explicit_dir and os.path.isdir(explicit_dir):
        if tumor_sample_id:
            for c in [os.path.join(explicit_dir, tumor_sample_id), os.path.join(explicit_dir, tumor_sample_id, tumor_sample_id)]:
                if os.path.isdir(c):
                    return c
        if group_id:
            for c in [os.path.join(explicit_dir, group_id, "ichorCNA_panel"), os.path.join(explicit_dir, group_id, "panel_ichorCNA"), os.path.join(explicit_dir, group_id)]:
                if os.path.isdir(c):
                    return c
        return explicit_dir

    search_bases = [output_dir]
    if group_id:
        search_bases.append(os.path.join(output_dir, group_id))
    parent_dir = os.path.dirname(os.path.abspath(output_dir)) if output_dir else None
    if parent_dir and os.path.isdir(parent_dir):
        search_bases.append(parent_dir)
        if group_id:
            search_bases.append(os.path.join(parent_dir, group_id))

    for base in search_bases:
        if not base or not os.path.isdir(base):
            continue
        for sub in ["ichorCNA_panel", "panel_ichorCNA"]:
            cand = os.path.join(base, sub)
            if os.path.isdir(cand):
                if tumor_sample_id:
                    cand_sample = os.path.join(cand, tumor_sample_id)
                    cand_nested = os.path.join(cand, tumor_sample_id, tumor_sample_id)
                    if os.path.isdir(cand_sample):
                        return cand_sample
                    if os.path.isdir(cand_nested):
                        return cand_nested
                for f in os.listdir(cand):
                    sub_cand = os.path.join(cand, f)
                    if os.path.isdir(sub_cand) and not f.startswith('.'):
                        for s_file in os.listdir(sub_cand):
                            if s_file.endswith('.params.txt'):
                                return sub_cand
                return cand
    return None

def parse_ichorcna_metrics(ichor_dir, sample_id, patient_id=None):
    """
    Locates and parses the ichorCNA params/metrics file for a sample.
    Recursively scans the directory and parent directory if needed, prioritizing tumor (-T) subfolders and matching sample ID.
    Returns: (purity, ploidy, sex) as floats/strings, or (None, None, None).
    """
    if not ichor_dir or not os.path.exists(ichor_dir):
        return None, None, None
        
    dirs_to_check = [ichor_dir]
    p_dir = os.path.dirname(os.path.abspath(ichor_dir))
    if p_dir and os.path.isdir(p_dir) and p_dir != os.path.abspath(ichor_dir):
        dirs_to_check.append(p_dir)

    candidates = []
    seen_paths = set()
    for d in dirs_to_check:
        for root, dirs, files in os.walk(d):
            for f_name in files:
                if not f_name.endswith('.txt') or f_name.startswith('._'):
                    continue
                if 'params' in f_name.lower() or 'final_metrics' in f_name.lower():
                    f_path = os.path.join(root, f_name)
                    if f_path in seen_paths:
                        continue
                    seen_paths.add(f_path)
                    score = 0
                    if sample_id and sample_id in f_name:
                        score += 50
                    elif sample_id and sample_id in root:
                        score += 25
                    if patient_id and patient_id in f_name:
                        score += 20
                    elif patient_id and patient_id in root:
                        score += 10
                    if '-T' in root or '-T' in f_name or 'tumor' in root.lower() or 't_vs_n' in root.lower() or 'plasma' in root.lower():
                        score += 10
                    if 't_vs_n' in root.lower() or 'plasma_vs_n' in root.lower():
                        score += 5
                    if '-N' in root or '-N' in f_name or 'normal' in root.lower():
                        score -= 20
                    candidates.append((score, f_path, f_name))
                
    candidates.sort(key=lambda x: -x[0])
    
    purity = None
    ploidy_val = None
    gender = None
    
    for score, f_path, f_name in candidates:
        try:
            with open(f_path, 'r', encoding='utf-8') as f:
                if 'final_metrics' in f_name.lower():
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        tf = row.get('Tumour.Fraction') or row.get('Tumor.Fraction') or row.get('"Tumour.Fraction"') or row.get('"Tumor.Fraction"')
                        pl = row.get('Ploidy') or row.get('"Ploidy"')
                        if tf is not None and pl is not None:
                            tf = tf.replace('"', '').strip()
                            pl = pl.replace('"', '').strip()
                            purity = float(tf)
                            ploidy_val = float(pl)
                            break
                else:
                    lines = [line.strip() for line in f if line.strip()]
                    for line in lines:
                        if 'gender' in line.lower() and ':' in line:
                            gender = line.split(':', 1)[1].strip()
                        if 'tumor fraction' in line.lower() and ':' in line:
                            try:
                                purity = float(line.split(':', 1)[1].strip())
                            except ValueError:
                                pass
                        if 'ploidy' in line.lower() and ':' in line:
                            try:
                                ploidy_val = float(line.split(':', 1)[1].strip())
                            except ValueError:
                                pass
                                
                    for i, line in enumerate(lines):
                        if line.startswith('Sample\tTumor Fraction\tPloidy') or line.startswith('Sample\tTumour Fraction\tPloidy'):
                            if i + 1 < len(lines):
                                parts = lines[i+1].split('\t')
                                if len(parts) >= 3:
                                    try:
                                        purity = float(parts[1])
                                        ploidy_val = float(parts[2])
                                    except ValueError:
                                        pass
                                    
            if purity is not None or ploidy_val is not None or gender is not None:
                return purity, ploidy_val, gender
        except Exception as e:
            print(f"Warning: Failed to parse ichorCNA file {f_name}: {e}", file=sys.stderr)
            
    return purity, ploidy_val, gender

def compute_all_purity_sources(purple_dir, ichor_dir, panel_ichor_dir, group_id, subject_id=None, tumor_sample_id=None, orthogonal_data=None):
    """
    Computes and extracts purity and ploidy estimates across all available platforms:
    - PURPLE (*.purple.purity.tsv)
    - sWGS ichorCNA (ichorCNA/*.params.txt)
    - Panel off-target ichorCNA (ichorCNA_panel/*.params.txt)
    - Clonal Somatic VAF from high-depth panel sequencing (2 * median clonal somatic VAF)
    - Pathology Tumor Cellularity (from IHC or orthogonal data)
    """
    res = {
        'purple': None,
        'purple_ploidy': None,
        'purple_method': None,
        'swgs': None,
        'swgs_ploidy': None,
        'panel_ichor': None,
        'panel_ploidy': None,
        'vaf': None,
        'vaf_stats': None,
        'pathology': None,
        'gender': None
    }
    
    # 1. PURPLE
    if purple_dir and os.path.exists(purple_dir):
        purity_file = find_purple_file(purple_dir, "purity.tsv", group_id, sample_id=tumor_sample_id)
        if purity_file and os.path.exists(purity_file):
            try:
                with open(purity_file, mode='r', encoding='utf-8') as f:
                    for r in csv.DictReader(f, delimiter='\t'):
                        if r.get('purity') is not None and r.get('purity') != 'NA':
                            p_val = float(r['purity'])
                            if p_val > 0.0:
                                res['purple'] = p_val
                        if r.get('ploidy') is not None and r.get('ploidy') != 'NA':
                            res['purple_ploidy'] = float(r['ploidy'])
                        if r.get('gender') is not None:
                            res['gender'] = r['gender']
                        if r.get('method') is not None:
                            res['purple_method'] = r['method']
            except Exception as e:
                print(f"Warning: Failed to parse PURPLE purity file: {e}", file=sys.stderr)
                
    # 2. sWGS ichorCNA
    if ichor_dir and os.path.exists(ichor_dir):
        pur, pl, sx = parse_ichorcna_metrics(ichor_dir, tumor_sample_id or group_id, subject_id)
        if pur is not None:
            res['swgs'] = pur
        if pl is not None:
            res['swgs_ploidy'] = pl
        if sx is not None and not res['gender']:
            res['gender'] = sx
            
    # 3. Panel ichorCNA
    if panel_ichor_dir and os.path.exists(panel_ichor_dir):
        pur, pl, sx = parse_ichorcna_metrics(panel_ichor_dir, tumor_sample_id or group_id, subject_id)
        if pur is not None:
            res['panel_ichor'] = pur
        if pl is not None:
            res['panel_ploidy'] = pl
        if sx is not None and not res['gender']:
            res['gender'] = sx
            
    # 4. Clonal Somatic VAF from high-depth panel sequencing
    if purple_dir and os.path.exists(purple_dir):
        somatic_vcf = find_purple_file(purple_dir, "somatic.vcf.gz", group_id, sample_id=tumor_sample_id)
        if somatic_vcf and os.path.exists(somatic_vcf):
            try:
                import gzip
                pass_vafs = []
                with gzip.open(somatic_vcf, 'rt') as f:
                    tumor_col_idx = 10
                    for line in f:
                        if line.startswith('#'):
                            if line.startswith('#CHROM'):
                                headers = line.strip().split('\t')
                                for idx in range(9, len(headers)):
                                    h = headers[idx]
                                    if h.startswith('T') or '-T' in h or 'tumor' in h.lower():
                                        tumor_col_idx = idx
                            continue
                        parts = line.strip().split('\t')
                        if parts[6] != 'PASS':
                            continue
                        fmt = parts[8].split(':')
                        if len(parts) > tumor_col_idx:
                            tumor_sample = parts[tumor_col_idx].split(':')
                            s_dict = dict(zip(fmt, tumor_sample))
                            af_str = s_dict.get('AF')
                            dp_str = s_dict.get('DP')
                            if af_str and dp_str:
                                af = float(af_str)
                                dp = int(dp_str)
                                if dp >= 30 and af >= 0.05:
                                    pass_vafs.append(af)
                if pass_vafs:
                    pass_vafs.sort()
                    med = pass_vafs[len(pass_vafs)//2]
                    p75 = pass_vafs[int(len(pass_vafs)*0.75)]
                    p90 = pass_vafs[int(len(pass_vafs)*0.90)]
                    res['vaf'] = min(1.0, 2.0 * med)
                    res['vaf_stats'] = {
                        'n_pass': len(pass_vafs),
                        'median_vaf': med,
                        'p75_vaf': p75,
                        'p90_vaf': p90
                    }
            except Exception:
                pass

    # 5. Pathology cellularity
    if orthogonal_data and orthogonal_data.get('pathology_purity') is not None:
        res['pathology'] = orthogonal_data['pathology_purity']
        
    return res

def resolve_sample_purity(purity_dict, requested_source="purple"):
    """
    Selects the active purity and ploidy based on the requested source,
    evaluates consistency against clonal somatic VAFs, and constructs
    user-facing sanity checks and suggestions.
    
    Returns: (purity_float, ploidy_float, selected_source, purity_notice)
    """
    req = (requested_source or "purple").lower().strip()
    p_pur = purity_dict.get('purple')
    s_pur = purity_dict.get('swgs')
    pan_pur = purity_dict.get('panel_ichor')
    v_pur = purity_dict.get('vaf')
    v_stats = purity_dict.get('vaf_stats')
    path_pur = purity_dict.get('pathology')
    
    selected_source = req
    purity_float = None
    ploidy_float = None
    purity_notice = None
    
    # 1. Selection logic
    if req == 'auto':
        # Smart decision heuristic:
        if p_pur is not None and p_pur > 0.0:
            if v_pur is not None and v_stats and v_stats.get('n_pass', 0) >= 10:
                n_p = v_stats['n_pass']
                m_vaf = v_stats['median_vaf']
                if p_pur < 0.5 * v_pur:
                    selected_source = 'vaf'
                    purity_float = v_pur
                    ploidy_float = purity_dict.get('purple_ploidy') or 2.0
                    purity_notice = f"Auto-selected Clonal VAF Purity ({v_pur*100:.1f}%) because PURPLE ({p_pur*100:.1f}%) significantly under-called relative to {n_p} PASS clonal somatic variants (median VAF: {m_vaf*100:.1f}%)."
                elif p_pur > 2.5 * v_pur and v_pur > 0.1:
                    selected_source = 'vaf'
                    purity_float = v_pur
                    ploidy_float = purity_dict.get('purple_ploidy') or 2.0
                    purity_notice = f"Auto-selected Clonal VAF Purity ({v_pur*100:.1f}%) because PURPLE ({p_pur*100:.1f}%) was substantially higher than observed clonal VAF peak ({m_vaf*100:.1f}%)."
                else:
                    selected_source = 'purple'
                    purity_float = p_pur
                    ploidy_float = purity_dict.get('purple_ploidy')
            else:
                selected_source = 'purple'
                purity_float = p_pur
                ploidy_float = purity_dict.get('purple_ploidy')
        elif v_pur is not None and v_stats and v_stats.get('n_pass', 0) >= 5:
            selected_source = 'vaf'
            purity_float = v_pur
            ploidy_float = purity_dict.get('purple_ploidy') or 2.0
            purity_notice = f"Auto-selected Clonal VAF Purity ({v_pur*100:.1f}%) as primary fallback (PURPLE model unavailable)."
        elif s_pur is not None and s_pur >= 0.05:
            selected_source = 'swgs'
            purity_float = s_pur
            ploidy_float = purity_dict.get('swgs_ploidy')
        elif pan_pur is not None:
            selected_source = 'panel_ichor'
            purity_float = pan_pur
            ploidy_float = purity_dict.get('panel_ploidy')
        elif path_pur is not None:
            selected_source = 'pathology'
            purity_float = path_pur
            ploidy_float = 2.0
    elif req == 'swgs':
        if s_pur is not None:
            purity_float = s_pur
            ploidy_float = purity_dict.get('swgs_ploidy')
        else:
            purity_float = p_pur or pan_pur or v_pur
            selected_source = 'fallback'
    elif req in ['panel', 'panel_ichor']:
        if pan_pur is not None:
            purity_float = pan_pur
            ploidy_float = purity_dict.get('panel_ploidy')
        else:
            purity_float = p_pur or s_pur or v_pur
            selected_source = 'fallback'
    elif req == 'vaf':
        if v_pur is not None:
            purity_float = v_pur
            ploidy_float = purity_dict.get('purple_ploidy') or 2.0
        else:
            purity_float = p_pur or s_pur or pan_pur
            selected_source = 'fallback'
    elif req == 'pathology':
        if path_pur is not None:
            purity_float = path_pur
            ploidy_float = purity_dict.get('purple_ploidy') or 2.0
        else:
            purity_float = p_pur or s_pur or pan_pur
            selected_source = 'fallback'
    else:
        # Default: PURPLE
        selected_source = 'purple'
        if p_pur is not None and p_pur > 0.0:
            purity_float = p_pur
            ploidy_float = purity_dict.get('purple_ploidy')
        elif v_pur is not None:
            selected_source = 'vaf (fallback)'
            purity_float = v_pur
            ploidy_float = 2.0
            purity_notice = f"PURPLE purity was unavailable; fell back to Clonal Somatic VAF implied purity ({v_pur*100:.1f}%)."
        elif s_pur is not None:
            selected_source = 'swgs (fallback)'
            purity_float = s_pur
            ploidy_float = purity_dict.get('swgs_ploidy')
            purity_notice = f"PURPLE purity was unavailable; fell back to sWGS ichorCNA ({s_pur*100:.1f}%)."
        elif pan_pur is not None:
            selected_source = 'panel_ichor (fallback)'
            purity_float = pan_pur
            ploidy_float = purity_dict.get('panel_ploidy')
            purity_notice = f"PURPLE purity was unavailable; fell back to Panel ichorCNA ({pan_pur*100:.1f}%)."

    # 2. Sanity check for PURPLE (Hybrid check requested by user)
    if 'purple' in selected_source and p_pur is not None and v_pur is not None and v_stats:
        med_vaf = v_stats.get('median_vaf', 0)
        n_pass = v_stats.get('n_pass', 0)
        if n_pass >= 5:
            if p_pur < 0.6 * v_pur:
                swgs_info = f"sWGS ichorCNA: {s_pur*100:.1f}%, " if s_pur is not None else ""
                panel_info = f"Panel ichorCNA: {pan_pur*100:.1f}%, " if pan_pur is not None else ""
                purity_notice = (
                    f"Selected PURPLE purity is {p_pur*100:.1f}%, but observed clonal somatic VAF "
                    f"(median {med_vaf*100:.1f}%, N={n_pass}) implies ~{v_pur*100:.1f}% tumor cellularity "
                    f"({swgs_info}{panel_info}implied VAF: {v_pur*100:.1f}%). "
                    f"If copy number alterations appear compressed, consider re-running with '--purity-source vaf' or '--purity-source auto'."
                )
            elif p_pur > 2.2 * v_pur and v_pur > 0.15:
                purity_notice = (
                    f"Selected PURPLE purity is {p_pur*100:.1f}%, which is substantially higher than implied by "
                    f"median clonal somatic VAF ({med_vaf*100:.1f}%, implied {v_pur*100:.1f}%). "
                    f"If subclonal heterogeneity or copy number over-scaling is suspected, consider '--purity-source vaf' or '--purity-source auto'."
                )

    if ploidy_float is None:
        ploidy_float = purity_dict.get('purple_ploidy') or purity_dict.get('swgs_ploidy') or purity_dict.get('panel_ploidy') or 2.0

    return purity_float, ploidy_float, selected_source, purity_notice

def parse_ichorcna_mad(ichor_dir):
    """
    Parses the MAD score (copy number noise level) from the ichorCNA params.txt file.
    """
    if not ichor_dir or not os.path.exists(ichor_dir):
        return None
    dirs_to_check = [ichor_dir]
    p_dir = os.path.dirname(os.path.abspath(ichor_dir))
    if p_dir and os.path.isdir(p_dir) and p_dir != os.path.abspath(ichor_dir):
        dirs_to_check.append(p_dir)

    seen_paths = set()
    for d in dirs_to_check:
        for root, _, files in os.walk(d):
            for f in files:
                if f.endswith('.params.txt') and not f.startswith('._'):
                    f_path = os.path.join(root, f)
                    if f_path in seen_paths:
                        continue
                    seen_paths.add(f_path)
                    try:
                        with open(f_path, 'r') as fp:
                            for line in fp:
                                if ('mad' in line.lower() or 'gc-map correction mad' in line.lower()) and ':' in line:
                                    val_str = line.split(':', 1)[1].strip()
                                    if val_str and val_str.upper() != 'NA':
                                        return float(val_str)
                    except Exception:
                        pass
    return None

def parse_ichorcna_qc_stats(ichor_dir, allow_swgs=False):
    """
    Locates and parses the qc_stats.json file if available in the ichorCNA directory.
    Also attaches the parsed MAD score from params.txt if available.
    """
    if not ichor_dir or not os.path.exists(ichor_dir):
        return None
    res_stats = None
    direct_json = os.path.join(ichor_dir, "qc_stats.json")
    if os.path.exists(direct_json):
        try:
            with open(direct_json, 'r') as f:
                res_stats = json.load(f)
        except Exception as e:
            print(f"Warning: Failed to parse qc_stats.json at {direct_json}: {e}", file=sys.stderr)
    if not res_stats:
        for root, dirs, files in os.walk(ichor_dir):
            if "qc_stats.json" in files:
                try:
                    with open(os.path.join(root, "qc_stats.json"), 'r') as f:
                        res_stats = json.load(f)
                        break
                except Exception as e:
                    print(f"Warning: Failed to parse qc_stats.json: {e}", file=sys.stderr)
            
    if res_stats and not allow_swgs:
        if res_stats.get('source') == 'swgs':
            return None
        # Detect legacy sWGS qc_stats.json where off_target was set to 0.0 and on-target equaled total reads
        if res_stats.get('off_target_coverage') == 0.0 and res_stats.get('off_target_reads') == 0 and res_stats.get('total_reads') == res_stats.get('on_target_reads'):
            return None

    if res_stats is None:
        res_stats = {}
    mad = parse_ichorcna_mad(ichor_dir)
    if mad is not None:
        res_stats['mad'] = mad
    return res_stats if res_stats else None

def get_swgs_qc_stats(ichor_dir, group_id, subject_id):
    """
    Finds and calculates sWGS (shallow WGS) total reads and coverage if available.
    """
    if not ichor_dir or not os.path.exists(ichor_dir):
        return None
        
    stats = parse_ichorcna_qc_stats(ichor_dir, allow_swgs=True)
    if stats:
        # Determine coverage
        cov_val = stats.get('coverage') or stats.get('genome_coverage')
        # If not present or 0, check on_target_coverage (legacy artifact in swgs json)
        if (cov_val is None or float(cov_val) == 0.0) and stats.get('on_target_coverage'):
            try:
                if float(stats.get('on_target_coverage')) > 0.0:
                    cov_val = stats.get('on_target_coverage')
            except (ValueError, TypeError):
                pass
        if cov_val is not None:
            stats['coverage'] = float(cov_val)
            stats['off_target_coverage'] = float(cov_val)
            
        mad = parse_ichorcna_mad(ichor_dir)
        if mad is not None:
            stats['mad'] = mad
        return stats

    bam_file = None
    search_roots = [ichor_dir]
    for s_root in search_roots:
        if not s_root or not os.path.exists(s_root):
            continue
        for root, dirs, files in os.walk(s_root):
            for f in files:
                if f.endswith('.bam') and not f.startswith('.'):
                    if (group_id and group_id in f) or (subject_id and subject_id in f) or '-T' in f or 'tumor' in f.lower():
                        bam_file = os.path.join(root, f)
                        break
            if bam_file:
                break
        if bam_file:
            break

    res_dict = {}
    if bam_file:
        try:
            cmd = f"module load samtools 2>/dev/null || true; samtools idxstats {shlex.quote(bam_file)}"
            out = subprocess.check_output(["bash", "-c", cmd], stderr=subprocess.DEVNULL).decode()
            total_reads = sum(int(line.split('\t')[2]) for line in out.strip().split('\n') if line and len(line.split('\t')) >= 3)
            if total_reads > 0:
                cov = (total_reads * 150) / 3100000000.0
                res_dict = {
                    "total_reads": total_reads,
                    "coverage": cov,
                    "off_target_coverage": cov
                }
        except Exception:
            pass

    mad = parse_ichorcna_mad(ichor_dir)
    if mad is not None:
        res_dict['mad'] = mad
        
    return res_dict if res_dict else None

def get_panel_bam_qc_stats(purple_dir, group_id, samplesheet_path):
    """Calculates total reads for targeted panel from the tumor BAM/CRAM."""
    bams = get_sample_filepaths(samplesheet_path, group_id, purple_dir)
    tumor_bam = None
    for b in bams:
        base = os.path.basename(b).lower()
        if '-t' in base or 'tumor' in base or (group_id and group_id.lower() in base and not '-n' in base):
            tumor_bam = b
            break
    if not tumor_bam and bams:
        tumor_bam = bams[0]
        
    if not tumor_bam or not os.path.exists(tumor_bam):
        return None
        
    try:
        import subprocess, shlex
        cmd = f"module load samtools 2>/dev/null || true; samtools idxstats {shlex.quote(tumor_bam)}"
        out = subprocess.check_output(["bash", "-c", cmd], stderr=subprocess.DEVNULL).decode()
        total_reads = sum(int(line.split('\t')[2]) for line in out.strip().split('\n') if line and len(line.split('\t')) >= 3)
        if total_reads > 0:
            return {
                "total_reads": total_reads,
                "on_target_reads": None,
                "on_target_coverage": None,
                "off_target_reads": None,
                "off_target_coverage": None,
                "mad": None
            }
    except Exception:
        pass
    return None

def parse_ichorcna_segments(ichor_dir, sample_id, patient_id=None, genome="GRCh38_hmf"):
    """
    Parses the ichorCNA seg.txt file to determine copy number states for genes of interest.
    Recursively scans the directory, prioritizing tumor (-T) subfolders.
    Returns: dict of gene -> copy_number (int/float) or None.
    """
    results = {gene: None for gene in GENES_OF_INTEREST}
    if not ichor_dir or not os.path.exists(ichor_dir):
        return results
        
    candidates = []
    for root, dirs, files in os.walk(ichor_dir):
        for f_name in files:
            for ext in ['.seg.txt', '.seg', '.cna.seg']:
                if f_name.endswith(ext):
                    f_path = os.path.join(root, f_name)
                    score = 0
                    if '-T' in root or '-T' in f_name or 'tumor' in root.lower():
                        score += 10
                    if ext == '.seg.txt':
                        score += 2
                    candidates.append((score, f_path))
                    
    if not candidates:
        return results
        
    candidates.sort(key=lambda x: -x[0])
    seg_file = candidates[0][1]
        
    bounds_key = "GRCh38_hmf"
    if "37" in genome:
        bounds_key = "GRCh37_hmf"
    bounds = GENE_BOUNDS.get(bounds_key, {})
    
    try:
        with open(seg_file, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            delim = ',' if ',' in first_line else '\t'
            f.seek(0)
            
            reader = csv.DictReader(f, delimiter=delim)
            for row in reader:
                row_cleaned = {k.strip().lower().replace('.', '_').replace('"', ''): v.strip().replace('"', '') for k, v in row.items() if k is not None}
                
                chrom = None
                start_val = None
                end_val = None
                cn_val = None
                
                for k, v in row_cleaned.items():
                    if k == 'chrom' or k == 'chromosome' or k == 'chr':
                        chrom = v
                    elif k == 'start' or k == 'loc_start' or k.endswith('_start'):
                        start_val = v
                    elif k == 'end' or k == 'loc_end' or k.endswith('_end'):
                        end_val = v
                    elif k == 'copy_number' or k == 'cn' or k == 'corrected_copy_number' or k.endswith('_copy_number') or k.endswith('_corrected_copy_number') or k.endswith('_cn'):
                        cn_val = v
                
                if not chrom or not start_val or not end_val or cn_val is None:
                    continue
                    
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                if chrom == 'chr23':
                    chrom = 'chrX'
                elif chrom == 'chr24':
                    chrom = 'chrY'
                    
                try:
                    seg_start = int(start_val)
                    seg_end = int(end_val)
                    seg_cn = float(cn_val)
                except ValueError:
                    continue
                    
                for gene, g_info in bounds.items():
                    g_chrom = g_info['chrom']
                    g_start = g_info['start']
                    g_end = g_info['end']
                    
                    if chrom == g_chrom:
                        if max(seg_start, g_start) <= min(seg_end, g_end):
                            results[gene] = seg_cn
    except Exception as e:
        print(f"Warning: Failed to parse ichorCNA segment file {seg_file}: {e}", file=sys.stderr)
        
    return results

def evaluate_ihc_concordance(gene, ihc_results, gene_reports, purity_float=None):
    """Evaluates the clinico-genomic concordance for a MMR gene.
    Returns a dict with keys:
      - 'genomic_status': human readable genomic/epigenomic status
      - 'concordance': concordance label (e.g. 'Concordant (Loss)', 'Discordant', etc.)
      - 'class': badge class or status styling
      - 'explanation': brief scientific explanation of the concordance/discordance
      - 'is_discordant': boolean indicating if results are discordant
      - 'discordance_note': string describing the discordance
    """
    ihc_val = ihc_results.get(gene, "No Data") if ihc_results else "No Data"
    rep = gene_reports.get(gene, {})
    risk_val = rep.get('risk', '')
    meth_val = rep.get('meth', '')
    
    # Helper to check if a gene has genomic biallelic loss / promoter methylation
    def has_genomic_loss(g):
        r = gene_reports.get(g, {})
        if 'high risk' in str(r.get('risk', '')).lower():
            return True
            
        meth_str = str(r.get('meth', '')).lower()
        if 'promoter methylated' in meth_str:
            has_loh = r.get('loh') == 'LOH'
            has_mut = r.get('s_mut') != 'Wildtype' or r.get('g_mut') != 'Wildtype'
            if has_loh or has_mut:
                return True
                
            import re
            m = re.search(r't=(\d+\.\d+)', meth_str)
            if m and purity_float and purity_float > 0:
                t_beta = float(m.group(1))
                if t_beta >= (purity_float * 0.70):
                    return True
        return False

    # Helper to determine if a gene has IHC loss
    def has_ihc_loss(g):
        if not ihc_results:
            return False
        val = ihc_results.get(g, "No Data")
        return any(x in str(val).lower() for x in ['loss', 'deficient', 'absent'])

    # Check if EPCAM has a deletion/loss causing MSH2 silencing
    epcam_deleted = False
    if 'EPCAM' in gene_reports:
        epcam_rep = gene_reports['EPCAM']
        epcam_s_cn = str(epcam_rep.get('s_cn', '')).lower()
        epcam_g_cn = str(epcam_rep.get('g_cn', '')).lower()
        if 'loss' in epcam_s_cn or 'loss' in epcam_g_cn or 'deletion' in epcam_s_cn or 'deletion' in epcam_g_cn:
            msh2_rep = gene_reports.get('MSH2', {})
            msh2_s_cn = str(msh2_rep.get('s_cn', '')).lower()
            msh2_g_cn = str(msh2_rep.get('g_cn', '')).lower()
            msh2_deleted = 'loss' in msh2_s_cn or 'loss' in msh2_g_cn or 'deletion' in msh2_s_cn or 'deletion' in msh2_g_cn
            if not msh2_deleted:
                epcam_deleted = True

    has_mut = rep.get('s_mut') != 'Wildtype' or rep.get('g_mut') != 'Wildtype'
    has_loh = rep.get('loh') == 'LOH'
    is_methylated = 'promoter methylated' in str(meth_val).lower()

    if 'high risk' in str(risk_val).lower():
        genomic_status = "Biallelic Inactivation"
    elif is_methylated:
        if has_mut:
            genomic_status = "Biallelic Inactivation (Mutation & Methylation)"
        elif has_loh:
            genomic_status = "Biallelic Inactivation (Methylation & LOH)"
        else:
            genomic_status = "Promoter Hypermethylated"
    elif 'moderate risk' in str(risk_val).lower():
        genomic_status = "Monoallelic Alteration"
    else:
        genomic_status = "Wildtype"

    is_ihc_loss = any(x in str(ihc_val).lower() for x in ['loss', 'deficient', 'absent'])
    is_ihc_intact = any(x in str(ihc_val).lower() for x in ['intact', 'retained', 'normal', 'present'])
    
    gen_loss = has_genomic_loss(gene)
    
    concordance = "Not Evaluated"
    concordance_class = "badge-neutral"
    explanation = ""
    is_discordant = False
    discordance_note = ""

    if not is_ihc_loss and not is_ihc_intact:
        concordance = "Not Evaluated"
        concordance_class = "badge-neutral"
    elif is_ihc_loss:
        biologically_concordant = False
        if gen_loss:
            biologically_concordant = True
            concordance = "Concordant (Genomic Loss)"
            explanation = f"explained by primary genomic/epigenomic loss ({genomic_status})"
        elif gene == 'PMS2' and has_genomic_loss('MLH1') and has_ihc_loss('MLH1'):
            biologically_concordant = True
            concordance = "Concordant (MLH1 Loss - partner degradation)"
            explanation = "explained by secondary loss due to MLH1 inactivation"
        elif gene == 'MSH6' and has_genomic_loss('MSH2') and has_ihc_loss('MSH2'):
            biologically_concordant = True
            concordance = "Concordant (MSH2 Loss - partner degradation)"
            explanation = "explained by secondary loss due to MSH2 inactivation"
        elif gene == 'MSH2' and epcam_deleted:
            biologically_concordant = True
            concordance = "Concordant (EPCAM deletion - MSH2 silencing)"
            explanation = "explained by secondary silencing due to EPCAM deletion"
        elif gene == 'MSH6' and epcam_deleted and has_ihc_loss('MSH2'):
            biologically_concordant = True
            concordance = "Concordant (EPCAM deletion - partner degradation)"
            explanation = "explained by secondary loss due to MSH2 silencing via EPCAM deletion"

        if biologically_concordant:
            concordance_class = "badge-methylated"
        else:
            concordance = "Discordant"
            concordance_class = "badge-somatic"
            is_discordant = True
            if 'promoter methylated' in str(meth_val).lower():
                concordance = "Discordant (Single Hit)"
                discordance_note = "IHC shows loss of expression, but genomic profiling only identified a single somatic hit (Promoter Hypermethylation). The second inactivating hit is missing, which could be an undetected somatic mutation, structural variant, or microduplication of the unmethylated allele."
            elif has_mut:
                concordance = "Discordant (Single Hit)"
                discordance_note = f"IHC shows loss of expression, but genomic profiling only identified a single somatic hit (Monoallelic Mutation). The second inactivating hit is missing, which could be an undetected structural variant, deep intronic mutation, or promoter hypermethylation."
            elif has_loh:
                concordance = "Discordant (Single Hit)"
                discordance_note = "IHC shows loss of expression, but genomic profiling only identified a single somatic hit (LOH / Heterozygous Loss). The second inactivating hit is missing."
            else:
                discordance_note = f"IHC shows loss of expression, but genomic status is {genomic_status}. This suggests a potentially missed structural variant, non-coding mutation, or epigenetic event."
    elif is_ihc_intact:
        biologically_concordant = True
        if gen_loss:
            biologically_concordant = False
        elif gene == 'MSH2' and epcam_deleted:
            biologically_concordant = False
            
        if biologically_concordant:
            concordance = "Concordant (Intact)"
            concordance_class = "badge-methylated"
        else:
            concordance = "Discordant"
            concordance_class = "badge-somatic"
            is_discordant = True
            if gen_loss:
                discordance_note = f"Genomic profiling shows {genomic_status}, but IHC indicates intact nuclear expression. This can occur with functional loss variants (e.g., missense mutations) that retain epitope binding capacity."
            elif gene == 'MSH2' and epcam_deleted:
                discordance_note = f"EPCAM deletion is detected (which transcriptionally silences MSH2), but IHC indicates intact MSH2 expression."

    return {
        'genomic_status': genomic_status,
        'concordance': concordance,
        'class': concordance_class,
        'explanation': explanation,
        'is_discordant': is_discordant,
        'discordance_note': discordance_note
    }

def get_base64_img(img_path):
    import base64
    try:
        with open(img_path, "rb") as f:
            data = f.read()
            encoded = base64.b64encode(data).decode('utf-8')
            ext = os.path.splitext(img_path)[1].replace('.', '')
            return f"data:image/{ext};base64,{encoded}"
    except Exception as e:
        print(f"Warning: Failed to base64 encode {img_path}: {e}", file=sys.stderr)
def convert_pdf_to_png(pdf_path):
    import subprocess
    import shutil
    import tempfile
    
    gs_path = shutil.which("gs")
    if not gs_path and os.path.exists("/usr/local/bin/gs"):
        gs_path = "/usr/local/bin/gs"
        
    convert_path = shutil.which("convert")
    if not convert_path and os.path.exists("/opt/homebrew/bin/convert"):
        convert_path = "/opt/homebrew/bin/convert"
        
    temp_dir = tempfile.gettempdir()
    temp_png = os.path.join(temp_dir, f"temp_ichor_{os.path.basename(pdf_path)}.png")
    
    success = False
    if gs_path:
        try:
            cmd = [
                gs_path,
                "-dNOPAUSE",
                "-dBATCH",
                "-sDEVICE=png16m",
                "-r150",
                f"-sOutputFile={temp_png}",
                pdf_path
            ]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            success = os.path.exists(temp_png)
        except Exception as e:
            print(f"Warning: gs conversion failed: {e}", file=sys.stderr)
            
    if not success and convert_path:
        try:
            cmd = [
                convert_path,
                "-density", "150",
                pdf_path,
                temp_png
            ]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            success = os.path.exists(temp_png)
        except Exception as e:
            print(f"Warning: convert conversion failed: {e}", file=sys.stderr)
            
    if success:
        b64 = get_base64_img(temp_png)
        try:
            os.remove(temp_png)
        except Exception:
            pass
        return b64
        
    return None

def find_ichor_plot(dir_path, sample_id=None):
    if not dir_path or not os.path.exists(dir_path):
        return None
        
    ichor_b64 = None
    
    def score_candidate(path):
        score = 0
        base = os.path.basename(path)
        if base.startswith("._"):
            return -9999
        if sample_id and sample_id in base:
            score += 100
        if '-T' in base or '-T' in path or 'tumor' in path.lower():
            score += 50
        if '-N' in base or '-N' in path or 'normal' in path.lower():
            score -= 50
        if 'all_sols' in path:
            score -= 10
        return score

    # 1. Look for ichorCNA genomeWide PNG plot
    ichor_candidates = []
    for root, dirs, files in os.walk(dir_path):
        for name in files:
            if name.endswith(".png") and "genomeWide" in name and not name.startswith("._"):
                ichor_candidates.append(os.path.join(root, name))
                
    ichor_candidates = sorted(ichor_candidates, key=lambda x: -score_candidate(x))
    for c in ichor_candidates:
        if os.path.exists(c) and not os.path.basename(c).startswith("._"):
            ichor_b64 = get_base64_img(c)
            if ichor_b64:
                return ichor_b64
                
    # 2. Look for ichorCNA genomeWide PDF plot to convert
    pdf_candidates = []
    for root, dirs, files in os.walk(dir_path):
        for name in files:
            if (name.endswith("_genomeWide.pdf") or (name.endswith(".pdf") and "genomeWide" in name)) and not name.startswith("._"):
                pdf_candidates.append(os.path.join(root, name))
                    
    pdf_candidates = sorted(pdf_candidates, key=lambda x: -score_candidate(x))
    for c in pdf_candidates:
        if os.path.exists(c) and not os.path.basename(c).startswith("._"):
            ichor_b64 = convert_pdf_to_png(c)
            if ichor_b64:
                return ichor_b64
                
    return None

def find_genomic_plots(purple_dir, ichor_dir, group_id, panel_ichor_dir=None, sample_id=None):
    circos_b64 = None
    ichor_b64 = None
    panel_ichor_b64 = None
    
    # 1. Look for PURPLE Circos plot
    if purple_dir and os.path.exists(purple_dir):
        circos_candidates = []
        if sample_id:
            circos_candidates.extend([
                os.path.join(purple_dir, "plot", f"{sample_id}.circos.png"),
                os.path.join(purple_dir, f"{sample_id}.circos.png")
            ])
        circos_candidates.extend([
            os.path.join(purple_dir, "plot", f"{group_id}.circos.png"),
            os.path.join(purple_dir, f"{group_id}.circos.png")
        ])
        for root, dirs, files in os.walk(purple_dir):
            for name in files:
                if name.endswith("circos.png") and not name.startswith("._"):
                    circos_candidates.append(os.path.join(root, name))
                    
        for c in circos_candidates:
            if os.path.exists(c) and not os.path.basename(c).startswith("._"):
                circos_b64 = get_base64_img(c)
                if circos_b64:
                    break
                    
    # 2. Look for sWGS ichorCNA plot
    if ichor_dir:
        ichor_b64 = find_ichor_plot(ichor_dir, sample_id=sample_id or group_id)
        
    # 3. Look for panel ichorCNA plot
    if panel_ichor_dir:
        p_abs = os.path.abspath(panel_ichor_dir)
        i_abs = os.path.abspath(ichor_dir) if ichor_dir else ""
        i_abs_prefix = i_abs.rstrip(os.sep) + os.sep if i_abs else ""
        if not i_abs or (p_abs != i_abs and not p_abs.startswith(i_abs_prefix)):
            panel_ichor_b64 = find_ichor_plot(panel_ichor_dir, sample_id=sample_id or group_id)
            if panel_ichor_b64 == ichor_b64:
                panel_ichor_b64 = None
        
    return circos_b64, ichor_b64, panel_ichor_b64

def generate_clinical_html_report(group_id, subject_id, purity, ploidy, tmb, tmb_status, msi, msi_status, wgd, clinical_info, gene_reports, summary_dir, ihc_results=None, purity_float=None, purple_dir=None, ichor_dir=None, ichor_sex=None, panel_ichor_dir=None, qc_stats=None, swgs_qc_stats=None, orange_dir=None, tumor_sample_id=None, orthogonal_data=None, purity_dict=None, selected_source="purple", purity_notice=None):
    """Generate an interactive HTML clinical report with links to BAF plots and IGV reports."""
    html_file = os.path.join(summary_dir, f"{group_id}_clinical_mmr_report.html")
    
    is_deficient = "PROFICIENT" not in clinical_info['classification'].upper()
    status_class = "status-deficient" if is_deficient else "status-proficient"
    
    html = []
    html.append("<!DOCTYPE html>")
    html.append("<html lang='en'>")
    html.append("<head>")
    html.append("    <meta charset='UTF-8'>")
    html.append("    <meta name='viewport' content='width=device-width, initial-scale=1.0'>")
    html.append(f"    <title>MMR Clinical Genomics Report: {group_id}</title>")
    html.append("    <link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap' rel='stylesheet'>")
    html.append("    <style>")
    html.append("        body { font-family: 'Inter', sans-serif; background-color: #f8fafc; color: #1e293b; margin: 0; padding: 40px; }")
    html.append("        .container { max-width: 1200px; margin: 0 auto; background-color: #ffffff; padding: 40px; border-radius: 12px; box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.05), 0 2px 4px -2px rgb(0 0 0 / 0.05); }")
    html.append("        h1 { font-size: 28px; font-weight: 700; color: #0f172a; margin-top: 0; margin-bottom: 8px; border-bottom: 2px solid #e2e8f0; padding-bottom: 16px; }")
    html.append("        h2 { font-size: 20px; font-weight: 600; color: #1e293b; margin-top: 32px; margin-bottom: 16px; border-bottom: 1px solid #f1f5f9; padding-bottom: 8px; }")
    html.append("        h3 { font-size: 16px; font-weight: 600; color: #334155; margin-top: 24px; margin-bottom: 12px; }")
    html.append("        .meta { color: #64748b; font-size: 14px; margin-bottom: 24px; }")
    html.append("        table { width: 100%; border-collapse: collapse; margin-bottom: 24px; text-align: left; }")
    html.append("        th, td { padding: 12px 16px; border-bottom: 1px solid #e2e8f0; font-size: 14px; }")
    html.append("        th { background-color: #f8fafc; font-weight: 600; color: #475569; }")
    html.append("        tr:hover { background-color: #f8fafc; }")
    html.append("        .badge { display: inline-block; padding: 4px 8px; border-radius: 6px; font-size: 12px; font-weight: 500; }")
    html.append("        .badge-somatic { background-color: #fee2e2; color: #991b1b; }")
    html.append("        .badge-germline { background-color: #f3e8ff; color: #6b21a8; }")
    html.append("        .badge-methylated { background-color: #d1fae5; color: #065f46; }")
    html.append("        .badge-loh { background-color: #ffedd5; color: #c2410c; }")
    html.append("        .badge-neutral { background-color: #f1f5f9; color: #475569; }")
    html.append("        .badge-discordant { background-color: #fef3c7; color: #92400e; border: 1px solid #fde68a; }")
    html.append("        .executive-summary { padding: 24px; border-radius: 8px; margin-bottom: 32px; border-left: 6px solid; }")
    html.append("        .status-deficient { background-color: #fef2f2; border-color: #ef4444; color: #991b1b; }")
    html.append("        .status-proficient { background-color: #f0fdf4; border-color: #22c55e; color: #166534; }")
    html.append("        .summary-title { font-size: 18px; font-weight: 700; margin-top: 0; margin-bottom: 8px; }")
    html.append("        .summary-desc { font-size: 15px; margin: 0; line-height: 1.6; }")
    html.append("        a { color: #2563eb; text-decoration: none; font-weight: 500; }")
    html.append("        a:hover { text-decoration: underline; color: #1d4ed8; }")
    html.append("        .gene-link { font-weight: 600; color: #0f172a; border-bottom: 1px dotted #94a3b8; }")
    html.append("        .variant-link { font-family: monospace; font-size: 13px; background-color: #f1f5f9; padding: 2px 6px; border-radius: 4px; color: #2563eb; }")
    html.append("        ul { padding-left: 20px; margin: 0; }")
    html.append("        li { font-size: 14px; margin-bottom: 8px; line-height: 1.5; }")
    html.append("        .tab-container { margin-top: 16px; border: 1px solid #e2e8f0; border-radius: 8px; overflow: hidden; background-color: #ffffff; }")
    html.append("        .tabs { display: flex; background-color: #f8fafc; border-bottom: 1px solid #e2e8f0; }")
    html.append("        .tab-btn { background: none; border: none; padding: 12px 24px; font-size: 14px; font-weight: 600; color: #64748b; cursor: pointer; transition: all 0.2s ease; border-right: 1px solid #e2e8f0; }")
    html.append("        .tab-btn:hover { background-color: #f1f5f9; color: #0f172a; }")
    html.append("        .tab-btn.active { background-color: #ffffff; color: #2563eb; border-bottom: 2px solid #2563eb; padding-bottom: 10px; }")
    html.append("        .tab-content { display: none; }")
    html.append("        .tab-content.active-content { display: block; }")
    html.append("    </style>")
    html.append("</head>")
    html.append("<body>")
    html.append("    <div class='container'>")
    html.append(f"        <h1>CLINICAL GENOMICS REPORT: MISMATCH REPAIR (MMR) STATUS</h1>")
    html.append(f"        <div class='meta'>Report Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</div>")
    
    # 1. Patient & Sample Information
    html.append("        <h2>1. Patient & Sample Information</h2>")
    html.append("        <table>")
    html.append("            <tr><th>Attribute</th><th>Value</th></tr>")
    html.append(f"            <tr><td><strong>Patient ID</strong></td><td>{subject_id}</td></tr>")
    html.append(f"            <tr><td><strong>Sample ID</strong></td><td>{group_id}</td></tr>")
    if tumor_sample_id and tumor_sample_id != group_id:
        html.append(f"            <tr><td><strong>Tumor Biospecimen ID</strong></td><td>{tumor_sample_id}</td></tr>")
    if ichor_sex:
        html.append(f"            <tr><td><strong>Determined Sex</strong></td><td>{ichor_sex}</td></tr>")
    
    # Tumor Purity Row with Multi-Platform Chips
    purity_row_html = [f"<div style='display: flex; align-items: center; gap: 8px; margin-bottom: 6px;'><strong style='font-size: 15px; color: #0f172a;'>{purity}</strong><span class='badge' style='background-color: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; font-size: 11px;'>Active: {selected_source.upper()}</span></div>"]
    if purity_dict:
        chips = []
        if purity_dict.get('purple') is not None:
            p_val = purity_dict['purple']
            chips.append(f"<span style='background: #f8fafc; padding: 3px 8px; border-radius: 4px; border: 1px solid #e2e8f0; color: #334155;'><strong>PURPLE:</strong> {p_val*100:.1f}%</span>")
        if purity_dict.get('swgs') is not None:
            s_val = purity_dict['swgs']
            chips.append(f"<span style='background: #f8fafc; padding: 3px 8px; border-radius: 4px; border: 1px solid #e2e8f0; color: #334155;'><strong>sWGS ichorCNA:</strong> {s_val*100:.1f}%</span>")
        if purity_dict.get('panel_ichor') is not None:
            pan_val = purity_dict['panel_ichor']
            chips.append(f"<span style='background: #f8fafc; padding: 3px 8px; border-radius: 4px; border: 1px solid #e2e8f0; color: #334155;'><strong>Panel ichorCNA:</strong> {pan_val*100:.1f}%</span>")
        if purity_dict.get('vaf') is not None:
            v_val = purity_dict['vaf']
            v_stats = purity_dict.get('vaf_stats') or {}
            extra = f" (med: {v_stats['median_vaf']*100:.1f}%, N={v_stats['n_pass']})" if 'median_vaf' in v_stats else ""
            chips.append(f"<span style='background: #f8fafc; padding: 3px 8px; border-radius: 4px; border: 1px solid #e2e8f0; color: #334155;'><strong>Clonal VAF Implied:</strong> {v_val*100:.1f}%{extra}</span>")
        path_val = purity_dict.get('pathology')
        if path_val is None and orthogonal_data and orthogonal_data.get('pathology_purity') is not None:
            path_val = orthogonal_data['pathology_purity']
        if path_val is not None:
            chips.append(f"<span style='background: #f8fafc; padding: 3px 8px; border-radius: 4px; border: 1px solid #e2e8f0; color: #334155;'><strong>Pathology:</strong> {path_val*100:.0f}%</span>")
        if chips:
            purity_row_html.append("<div style='display: flex; flex-wrap: wrap; gap: 6px; font-size: 12px; margin-top: 4px;'>" + "".join(chips) + "</div>")
    html.append(f"            <tr><td><strong>Tumor Purity</strong></td><td>{''.join(purity_row_html)}</td></tr>")

    if orthogonal_data and orthogonal_data.get('pathology_purity') is not None:
        path_pur = orthogonal_data['pathology_purity']
        html.append(f"            <tr><td><strong>Pathology Tumor Cellularity</strong></td><td>{path_pur*100:.0f}%</td></tr>")
    html.append(f"            <tr><td><strong>Tumor Ploidy</strong></td><td>{ploidy}</td></tr>")
    html.append(f"            <tr><td><strong>Tumor Mutational Burden (TMB)</strong></td><td>{tmb} ({tmb_status})</td></tr>")
    html.append(f"            <tr><td><strong>Microsatellite Instability (MSI)</strong></td><td>{msi} ({msi_status})</td></tr>")
    html.append(f"            <tr><td><strong>Whole Genome Duplication (WGD)</strong></td><td>{wgd}</td></tr>")
    if orthogonal_data and orthogonal_data.get('test_result'):
        html.append(f"            <tr><td><strong>Orthogonal Assay Result</strong></td><td><strong style='color: #2563eb;'>{orthogonal_data['test_result']}</strong></td></tr>")
    
    # Check for Orange PDF report
    orange_pdf = None
    if orange_dir:
        orange_pdf = find_orange_file(orange_dir, ".orange.pdf", group_id, sample_id=tumor_sample_id)
    if not orange_pdf and purple_dir:
        p_orange = os.path.join(os.path.dirname(purple_dir), 'orange')
        if os.path.isdir(p_orange):
            orange_pdf = find_orange_file(p_orange, ".orange.pdf", group_id, sample_id=tumor_sample_id)
    if orange_pdf and os.path.exists(orange_pdf):
        rel_orange_pdf = os.path.relpath(orange_pdf, summary_dir)
        html.append(f"            <tr><td><strong>Orange Comprehensive Report</strong></td><td><a href='{rel_orange_pdf}' target='_blank' class='variant-link' style='font-weight: 600; padding: 4px 10px; font-size: 13px; background-color: #ffedd5; color: #c2410c; text-decoration: none; border-radius: 4px; display: inline-flex; align-items: center; gap: 6px;'>📄 Open Orange PDF Report</a></td></tr>")
    html.append("        </table>")
    
    # Tumor Purity Discrepancy Advisory Banner
    if purity_notice:
        html.append("        <div style='background-color: #fffbeb; border: 1px solid #fef3c7; border-left: 5px solid #f59e0b; padding: 14px 18px; border-radius: 8px; margin: 16px 0 24px 0; color: #92400e;'>")
        html.append("            <div style='font-weight: 700; font-size: 14px; margin-bottom: 4px; display: flex; align-items: center; gap: 6px;'>")
        html.append("                <span style='font-size: 16px;'>⚠️</span> <strong>Tumor Purity Advisory & Clonal VAF Consistency</strong>")
        html.append("            </div>")
        html.append(f"            <div style='font-size: 13.5px; line-height: 1.5;'>{purity_notice}</div>")
        html.append("        </div>")
    
    # 2. Executive Summary
    if purity_float is not None and purity_float < 0.3:
        html.append(f"        <div class='executive-summary' style='background-color: #fffbeb; border-left: 4px solid #d97706; padding: 12px; margin-bottom: 20px; color: #92400e; border-radius: 4px;'>")
        html.append(f"            <div style='font-weight: 700; margin-bottom: 4px;'>⚠️ Low Tumor Purity Warning ({purity_float*100:.1f}%)</div>")
        html.append(f"            <p style='margin: 0; font-size: 14px;'>Tumor fraction is below the 30% quality control threshold. Copy number profiling and somatic calls have low confidence and may be artifactual.</p>")
        html.append("        </div>")
        
    html.append(f"        <div class='executive-summary {status_class}'>")
    html.append(f"            <div class='summary-title'>Clinical Classification: {clinical_info['classification']}</div>")
    html.append(f"            <p class='summary-desc'><strong>Confidence Level:</strong> {clinical_info['confidence']}</p>")
    html.append(f"            <p class='summary-desc' style='margin-top: 8px;'><strong>Evidence & Notes:</strong> {clinical_info['notes']}</p>")
    html.append("        </div>")
    
    # helper to format mutation cell into HTML links
    def format_mut_cell_html(mut_str, vars_list, gene):
        if mut_str == 'Wildtype' or not mut_str:
            return "<span class='badge badge-neutral'>Wildtype</span>"
        
        import re
        parts = [p.strip() for p in re.split(r',\s*(?![^()]*\))', mut_str)]
        formatted = []
        for p in parts:
            var_pos = None
            m_change = re.match(r'(\S+)\s*\(', p)
            change_clean = m_change.group(1) if m_change else p
            
            for v in vars_list:
                if change_clean in v['label']:
                    var_pos = v['pos']
                    break
                    
            if var_pos:
                formatted.append("<a href='igv/{}_igv_report.html?locus={}' target='_blank' class='variant-link'>{}</a>".format(group_id, var_pos, p))
            else:
                formatted.append("<span class='badge badge-neutral'>{}</span>".format(p))
        return ", ".join(formatted)
        
    def format_cn_cell_html(cn_val):
        suffix = " (ichorCNA)" if "ichorcna" in str(cn_val).lower() else ""
        val_clean = str(cn_val).replace(" (ichorCNA)", "").strip()
        
        if val_clean in ['Normal', '2', '2.0'] or val_clean.startswith('Normal'):
            return "<span class='badge badge-neutral'>Normal (2){}</span>".format(suffix)
        elif 'homozygous deletion' in val_clean.lower() or 'hom.del' in val_clean.lower():
            return "<span class='badge badge-somatic' style='background-color: #ef4444; color: #ffffff; font-weight: 600;'>Homozygous Deletion{}</span>".format(suffix)
        elif 'heterozygous loss' in val_clean.lower() or 'loss' in val_clean.lower() or 'het.loss' in val_clean.lower():
            return "<span class='badge badge-somatic'>Heterozygous Loss{}</span>".format(suffix)
        elif 'gain' in val_clean.lower():
            return "<span class='badge badge-loh'>{}</span>".format(cn_val)
        return "<span class='badge badge-loh'>{}{}</span>".format(cn_val, suffix)
        
    def format_loh_cell_html(loh_val):
        if loh_val == 'LOH':
            return "<span class='badge badge-loh'>LOH</span>"
        return "<span class='badge badge-neutral'>Normal</span>"
        
    def format_meth_cell_html(meth_val):
        if not meth_val or meth_val == 'No Data':
            return "<span class='badge badge-neutral' style='background-color: #475569; color: #f8fafc;'>No Data</span>"
        if 'unmethylated' in str(meth_val).lower() or meth_val == 'Unmethylated':
            return "<span class='badge badge-neutral'>Unmethylated</span>"
        return "<span class='badge badge-methylated'>{}</span>".format(meth_val)
        
    def format_risk_cell_html(risk_val):
        if "HIGH RISK" in risk_val:
            return "<span class='badge badge-somatic' style='font-weight: 700;'>{}</span>".format(risk_val)
        elif "MODERATE" in risk_val:
            return "<span class='badge badge-loh' style='font-weight: 600;'>{}</span>".format(risk_val)
        elif "LOW" in risk_val:
            return "<span class='badge badge-neutral' style='font-weight: 600;'>{}</span>".format(risk_val)
        return "<span class='badge badge-neutral'>{}</span>".format(risk_val)

    # 3. IHC & Genomic Correlation Table (if IHC results are provided)
    has_ihc = False
    if ihc_results:
        has_ihc = any(v not in ["No Data", "Not Evaluated"] for v in ihc_results.values())
        
    sec_num = 3
    if has_ihc:
        html.append(f"        <h2>{sec_num}. Immunohistochemistry (IHC) & Genomic Correlation</h2>")
        html.append("        <table>")
        html.append("            <tr>")
        html.append("                <th>Gene</th><th>IHC Clinical Result</th><th>Genomic / Epigenomic Status</th><th>Concordance</th>")
        html.append("            </tr>")
        
        discordance_notes = []
        for g in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
            ihc_val = ihc_results.get(g, "No Data")
            res = evaluate_ihc_concordance(g, ihc_results, gene_reports, purity_float)
            
            # Check for expression loss or intact to render IHC badge
            is_ihc_loss = any(x in str(ihc_val).lower() for x in ['loss', 'deficient', 'absent'])
            is_ihc_intact = any(x in str(ihc_val).lower() for x in ['intact', 'retained', 'normal', 'present'])
            
            if is_ihc_loss:
                ihc_badge = f"<span class='badge badge-somatic' style='font-weight: 600;'>{ihc_val}</span>"
            elif is_ihc_intact:
                ihc_badge = f"<span class='badge badge-methylated'>{ihc_val}</span>"
            else:
                ihc_badge = f"<span class='badge badge-neutral'>{ihc_val}</span>"
                
            if res['is_discordant']:
                discordance_notes.append(f"<strong>{g}</strong>: {res['discordance_note']}")
                
            html.append("            <tr>")
            html.append(f"                <td><strong>{g}</strong></td>")
            html.append(f"                <td>{ihc_badge}</td>")
            html.append(f"                <td>{res['genomic_status']}</td>")
            html.append(f"                <td><span class='badge {res['class']}'>{res['concordance']}</span></td>")
            html.append("            </tr>")
            
        html.append("        </table>")
        
        if discordance_notes:
            html.append("        <div class='executive-summary status-deficient' style='border-color: #fb923c; background-color: #fffbeb; color: #c2410c;'>")
            html.append("            <div class='summary-title'>⚠️ Clinico-Genomic Discordance Alert</div>")
            html.append("            <ul>")
            for note in discordance_notes:
                html.append(f"                <li style='color: #9a3412;'>{note}</li>")
            html.append("            </ul>")
            html.append("        </div>")
            
        sec_num += 1

    # Orthogonal Assay & Known Mutation Concordance (if provided)
    has_orth = False
    if orthogonal_data and (orthogonal_data.get('variants') or orthogonal_data.get('test_result')):
        has_orth = True
        
    if has_orth:
        html.append(f"        <h2>{sec_num}. Orthogonal Assay & Known Mutation Concordance</h2>")
        html.append("        <table>")
        html.append("            <tr>")
        html.append("                <th>Gene</th><th>Orthogonal Assay Finding</th><th>Assay Interpretation</th><th>Oncoanalyser Genomic Finding</th><th>Concordance Status</th>")
        html.append("            </tr>")
        
        orth_genes = list(dict.fromkeys([v['gene'] for v in orthogonal_data.get('variants', []) if v['gene'] not in ['None', 'MMR']]))
        if not orth_genes:
            orth_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2']
        else:
            for cg in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
                if cg not in orth_genes:
                    orth_genes.append(cg)
                    
        for g in orth_genes:
            c_res = evaluate_orthogonal_concordance(g, orthogonal_data, gene_reports)
            if c_res.get('has_data'):
                html.append("            <tr>")
                html.append(f"                <td><strong>{g}</strong></td>")
                html.append(f"                <td>{c_res['orthogonal_findings']}</td>")
                html.append(f"                <td>{c_res['clinical_significance']}</td>")
                html.append(f"                <td>{c_res['onco_findings']}</td>")
                html.append(f"                <td><span class='badge {c_res['badge_class']}'>{c_res['concordance']}</span></td>")
                html.append("            </tr>")
        html.append("        </table>")
        sec_num += 1

    # Mismatch Repair Core Genes
    html.append(f"        <h2>{sec_num}. Mismatch Repair (MMR) Core Genes</h2>")
    html.append("        <table>")
    html.append("            <tr>")
    html.append("                <th>Gene</th><th>Somatic Mutation</th><th>Germline Mutation</th><th>Somatic CN</th><th>Germline CN</th><th>LOH</th><th>Promoter Methylation</th>")
    html.append("            </tr>")
    
    mmr_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM', 'MLH3']
    for g in mmr_genes:
        rep = gene_reports[g]
        g_link = "<a href='plots/{}_{}_baf.html' target='_blank' class='gene-link' title='Click to view interactive BAF & Copy Number Plot'>{}</a>".format(group_id, g, g)
        s_mut_html = format_mut_cell_html(rep['s_mut'], rep['vars_list'], g)
        g_mut_html = format_mut_cell_html(rep['g_mut'], rep['vars_list'], g)
        s_cn_html = format_cn_cell_html(rep['s_cn'])
        g_cn_html = format_cn_cell_html(rep['g_cn'])
        loh_html = format_loh_cell_html(rep['loh'])
        meth_html = format_meth_cell_html(rep['meth'])
        
        html.append(f"            <tr>")
        html.append(f"                <td>{g_link}</td><td>{s_mut_html}</td><td>{g_mut_html}</td><td>{s_cn_html}</td><td>{g_cn_html}</td><td>{loh_html}</td><td>{meth_html}</td>")
        html.append(f"            </tr>")
    html.append("        </table>")
    sec_num += 1
    
    # 5. Other Driver Genes
    html.append(f"        <h2>{sec_num}. Other Driver & Ancillary Genes</h2>")
    html.append("        <table>")
    html.append("            <tr>")
    html.append("                <th>Gene</th><th>Somatic Mutation</th><th>Germline Mutation</th><th>Somatic CN</th><th>Germline CN</th><th>LOH</th><th>Promoter Methylation</th>")
    html.append("            </tr>")
    
    other_genes = ['POLD1', 'POLE', 'MUTYH', 'TP53', 'BRAF']
    for g in other_genes:
        rep = gene_reports[g]
        g_link = "<a href='plots/{}_{}_baf.html' target='_blank' class='gene-link' title='Click to view interactive BAF & Copy Number Plot'>{}</a>".format(group_id, g, g)
        s_mut_html = format_mut_cell_html(rep['s_mut'], rep['vars_list'], g)
        g_mut_html = format_mut_cell_html(rep['g_mut'], rep['vars_list'], g)
        s_cn_html = format_cn_cell_html(rep['s_cn'])
        g_cn_html = format_cn_cell_html(rep['g_cn'])
        loh_html = format_loh_cell_html(rep['loh'])
        meth_html = format_meth_cell_html(rep['meth'])
        
        html.append(f"            <tr>")
        html.append(f"                <td>{g_link}</td><td>{s_mut_html}</td><td>{g_mut_html}</td><td>{s_cn_html}</td><td>{g_cn_html}</td><td>{loh_html}</td><td>{meth_html}</td>")
        html.append(f"            </tr>")
    html.append("        </table>")
    sec_num += 1
    
    # 6. Detailed Variant Phasing
    html.append(f"        <h2>{sec_num}. Detailed Variant Phasing & Annotations</h2>")
    has_vars = False
    
    for g in GENES_OF_INTEREST:
        rep = gene_reports[g]
        if rep['vars_list']:
            has_vars = True
            html.append(f"        <h3>{g}</h3>")
            html.append("        <ul>")
            for v in rep['vars_list']:
                v_link = "<a href='igv/{}_igv_report.html?locus={}' target='_blank' class='variant-link'>{}</a>".format(group_id, v['pos'], v['label'])
                html.append(f"            <li>")
                html.append(f"                <strong>Variant:</strong> {v_link}<br>")
                html.append(f"                Phasing status: {v['phasing']}<br>")
                if v.get('lps') is not None:
                    html.append(f"                Phase Set Support: ID={v['lps']}, read count={v['lps_rc']}")
                html.append(f"            </li>")
            html.append("        </ul>")
            
    if not has_vars:
        html.append("        <p>No pathogenic mutations or promoter methylation events detected.</p>")
        
    # 7. Quality Control & Genomic Visualization
    if not ichor_dir and summary_dir:
        parent_out = os.path.dirname(os.path.abspath(summary_dir))
        ichor_dir = resolve_swgs_ichor_dir(None, parent_out, group_id, tumor_sample_id=tumor_sample_id)
    if not panel_ichor_dir and summary_dir:
        parent_out = os.path.dirname(os.path.abspath(summary_dir))
        panel_ichor_dir = resolve_panel_ichor_dir(None, parent_out, group_id, tumor_sample_id=tumor_sample_id)
    if not swgs_qc_stats and ichor_dir:
        swgs_qc_stats = get_swgs_qc_stats(ichor_dir, group_id, subject_id)

    circos_b64, ichor_b64, panel_ichor_b64 = find_genomic_plots(purple_dir, ichor_dir, group_id, panel_ichor_dir, sample_id=tumor_sample_id)
    
    plots = []
    if qc_stats or swgs_qc_stats:
        plots.append(('seq-qc', '⚡ Sequencing & Alignment QC', 'qc_table'))
    # Always include sWGS ichorCNA tab; if data/plot is not available, render informative placeholder
    plots.append(('swgs-ichor', 'sWGS ichorCNA Genome Plot', ichor_b64 if ichor_b64 else 'swgs_na'))
    if panel_ichor_b64 and panel_ichor_b64 != ichor_b64:
        plots.append(('panel-ichor', 'Panel ichorCNA Genome Plot', panel_ichor_b64))
    if circos_b64:
        plots.append(('purple-circos', 'PURPLE Circos Plot', circos_b64))
        
    if plots:
        sec_num += 1
        html.append(f"        <h2>{sec_num}. Quality Control & Genomic Visualizations</h2>")
        
        if len(plots) >= 2:
            html.append("        <div class='tab-container'>")
            html.append("            <div class='tabs'>")
            for i, (tab_id, tab_title, _) in enumerate(plots):
                active_class = " active" if i == 0 else ""
                html.append(f"                <button class='tab-btn{active_class}' onclick=\"openTab(event, '{tab_id}')\">{tab_title}</button>")
            html.append("            </div>")
            
            for i, (tab_id, tab_title, item_data) in enumerate(plots):
                display_style = "display: block;" if i == 0 else "display: none;"
                active_content_class = " active-content" if i == 0 else ""
                if tab_id == 'seq-qc':
                    total_reads = f"{qc_stats.get('total_reads', 0):,}" if qc_stats and qc_stats.get('total_reads') and qc_stats.get('total_reads') != 'NA' else "N/A"
                    on_target_reads = f"{qc_stats.get('on_target_reads', 0):,}" if qc_stats and qc_stats.get('on_target_reads') and qc_stats.get('on_target_reads') != 'NA' else "N/A"
                    on_target_cov = f"{float(qc_stats['on_target_coverage']):.2f}x" if qc_stats and qc_stats.get('on_target_coverage') and qc_stats.get('on_target_coverage') != 'NA' else "N/A"
                    off_target_reads = f"{qc_stats.get('off_target_reads', 0):,}" if qc_stats and qc_stats.get('off_target_reads') and qc_stats.get('off_target_reads') != 'NA' else "N/A"
                    off_target_cov = f"{float(qc_stats['off_target_coverage']):.2f}x" if qc_stats and qc_stats.get('off_target_coverage') and qc_stats.get('off_target_coverage') != 'NA' else "N/A"
                    p_mad_val = qc_stats.get('mad') if qc_stats else None
                    panel_mad = f"{float(p_mad_val):.4f}" if p_mad_val is not None and p_mad_val != 'NA' else "N/A (Not Run)"
                    
                    swgs_reads = f"{swgs_qc_stats.get('total_reads', 0):,}" if swgs_qc_stats and swgs_qc_stats.get('total_reads') and swgs_qc_stats.get('total_reads') != 'NA' else "N/A"
                    swgs_cov = "N/A"
                    if swgs_qc_stats:
                        c_val = swgs_qc_stats.get('coverage') or swgs_qc_stats.get('genome_coverage') or swgs_qc_stats.get('off_target_coverage')
                        if c_val is not None and c_val != 'NA':
                            try:
                                if float(c_val) > 0.0:
                                    swgs_cov = f"{float(c_val):.2f}x"
                            except (ValueError, TypeError):
                                pass
                    s_mad_val = swgs_qc_stats.get('mad') if swgs_qc_stats else None
                    swgs_mad = f"{float(s_mad_val):.4f}" if s_mad_val is not None and s_mad_val != 'NA' else "N/A"

                    html.append(f"            <div id='seq-qc' class='tab-content{active_content_class}' style='padding: 24px; {display_style}'>")
                    html.append("                <h3 style='margin-top: 0; margin-bottom: 16px; font-size: 16px; color: #1e293b;'>Sequencing, Alignment & Copy Number Noise (ichorCNA MAD) QC</h3>")
                    html.append("                <table style='margin-bottom: 0;'>")
                    html.append("                    <tr><th>QC Metric</th><th>Targeted Panel</th><th>Shallow WGS (sWGS)</th></tr>")
                    html.append(f"                    <tr><td><strong>Total Reads</strong></td><td>{total_reads}</td><td>{swgs_reads}</td></tr>")
                    html.append(f"                    <tr><td><strong>On-Target Reads</strong></td><td>{on_target_reads}</td><td>N/A (Genome-Wide)</td></tr>")
                    html.append(f"                    <tr><td><strong>On-Target Coverage</strong></td><td>{on_target_cov}</td><td>N/A (Genome-Wide)</td></tr>")
                    html.append(f"                    <tr><td><strong>Off-Target / Genome Reads</strong></td><td>{off_target_reads}</td><td>{swgs_reads}</td></tr>")
                    html.append(f"                    <tr><td><strong>Off-Target / Genome Coverage</strong></td><td>{off_target_cov}</td><td>{swgs_cov}</td></tr>")
                    html.append(f"                    <tr style='background-color: #f8fafc;'><td><strong>ichorCNA MAD Score (Copy Number Noise)</strong></td><td><strong style='color: #2563eb;'>{panel_mad}</strong></td><td><strong style='color: #2563eb;'>{swgs_mad}</strong></td></tr>")
                    html.append("                </table>")
                    html.append("            </div>")
                elif tab_id == 'swgs-ichor' and (item_data == 'swgs_na' or not item_data):
                    html.append(f"            <div id='swgs-ichor' class='tab-content{active_content_class}' style='text-align: center; padding: 48px 24px; {display_style}'>")
                    html.append("                <div style='max-width: 520px; margin: 0 auto; background: #f8fafc; border: 1.5px dashed #cbd5e1; border-radius: 12px; padding: 36px 24px;'>")
                    html.append("                    <div style='font-size: 36px; margin-bottom: 12px;'>🧬</div>")
                    html.append("                    <div style='font-size: 16px; font-weight: 700; color: #334155; margin-bottom: 6px;'>sWGS data not available</div>")
                    html.append("                    <div style='font-size: 13.5px; color: #64748b; line-height: 1.5;'>Low-pass whole-genome sequencing (sWGS) was not performed or ichorCNA genome-wide copy number data is not available for this sample.</div>")
                    html.append("                </div>")
                    html.append("            </div>")
                else:
                    img_b64 = item_data
                    max_width = "90%" if tab_id == 'purple-circos' else "100%"
                    html.append(f"            <div id='{tab_id}' class='tab-content{active_content_class}' style='text-align: center; padding: 20px; {display_style}'>")
                    html.append(f"                <img src='{img_b64}' style='max-width: {max_width}; height: auto; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);' alt='{tab_title}'>")
                    html.append("            </div>")
            html.append("        </div>")
            
            html.append("        <script>")
            html.append("        function openTab(evt, tabName) {")
            html.append("            var i, tabcontent, tablinks;")
            html.append("            tabcontent = document.getElementsByClassName('tab-content');")
            html.append("            for (i = 0; i < tabcontent.length; i++) {")
            html.append("                tabcontent[i].style.display = 'none';")
            html.append("            }")
            html.append("            tablinks = document.getElementsByClassName('tab-btn');")
            html.append("            for (i = 0; i < tablinks.length; i++) {")
            html.append("                tablinks[i].className = tablinks[i].className.replace(' active', '');")
            html.append("            }")
            html.append("            document.getElementById(tabName).style.display = 'block';")
            html.append("            evt.currentTarget.className += ' active';")
            html.append("        }")
            html.append("        </script>")
        else:
            tab_id, tab_title, item_data = plots[0]
            if tab_id == 'swgs-ichor' and (item_data == 'swgs_na' or not item_data):
                html.append("        <div style='text-align: center; padding: 48px 24px;'>")
                html.append("            <div style='max-width: 520px; margin: 0 auto; background: #f8fafc; border: 1.5px dashed #cbd5e1; border-radius: 12px; padding: 36px 24px;'>")
                html.append("                <div style='font-size: 36px; margin-bottom: 12px;'>🧬</div>")
                html.append("                <div style='font-size: 16px; font-weight: 700; color: #334155; margin-bottom: 6px;'>sWGS data not available</div>")
                html.append("                <div style='font-size: 13.5px; color: #64748b; line-height: 1.5;'>Low-pass whole-genome sequencing (sWGS) was not performed or ichorCNA genome-wide copy number data is not available for this sample.</div>")
                html.append("            </div>")
                html.append("        </div>")
            else:
                max_width = "90%" if tab_id == 'purple-circos' else "100%"
                html.append("        <div style='text-align: center; padding: 20px;'>")
                html.append(f"            <h3 style='text-align: left;'>{tab_title}</h3>")
                html.append(f"            <img src='{item_data}' style='max-width: {max_width}; height: auto; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);' alt='{tab_title}'>")
                html.append("        </div>")
            
    html.append("    </div>")
    html.append("</body>")
    html.append("</html>")
    
    try:
        with open(html_file, 'w', encoding='utf-8') as f_out:
            f_out.write("\n".join(html))
        print(f"Clinical MMR HTML report successfully written to {html_file}")
    except Exception as e:
        print(f"Error: Failed to write clinical HTML report to {html_file}: {e}", file=sys.stderr)

def classify_lynch_vs_sporadic(patient_data: dict) -> dict:
    """Classifies a patient as having Lynch Syndrome, Sporadic MMRD, or Lynch-Like
    syndrome based on paired tumor-normal NGS data and optional IHC data.
    """
    msi_status = patient_data.get("msi_status", "MSS").upper()
    tmb = float(patient_data.get("tmb", 0.0))
    tissue = patient_data.get("tumor_tissue_origin", "unknown").lower()
    age = patient_data.get("age_at_diagnosis", 999)

    # Base check: Is the tumor actually mismatch repair deficient (dMMR)?
    is_dmmr_by_ngs = msi_status == "MSI-H" or tmb >= 10.0

    # Integrate IHC baseline if present
    ihc_data = patient_data.get("ihc_results", {})
    # Filter out empty or "No Data" fields
    valid_ihc = {k: v for k, v in ihc_data.items() if v and v != "No Data"}
    has_ihc = bool(valid_ihc)
    is_dmmr_by_ihc = any(
        status.lower() in ["absent", "loss", "loss of nuclear expression", "deficient"] 
        for status in valid_ihc.values()
    )

    # If both NGS and IHC confirm proficiency, or if NGS is proficient and no IHC exists
    if not is_dmmr_by_ngs and (not has_ihc or not is_dmmr_by_ihc):
        return {
            "classification": "MMR-Proficient (MSS / TMB-Low / Intact IHC)",
            "confidence": "High",
            "notes": "Tumor displays intact mismatch repair mechanisms.",
        }

    # 2. Rule in Hereditary Lynch Syndrome (Germline Variant Found)
    germline_variants = patient_data.get("germline_pathogenic_variants", [])
    if germline_variants:
        genes_affected = ", ".join(set(germline_variants))
        return {
            "classification": "CONFIRMED LYNCH SYNDROME (Hereditary)",
            "confidence": "High",
            "notes": f"Pathogenic germline variant detected in: {genes_affected}.",
        }

    # 2b. Rule in Hereditary Lynch Syndrome via Constitutional Epimutation
    germline_mlh1_methylated = patient_data.get("germline_mlh1_promoter_methylation", False)
    if germline_mlh1_methylated:
        return {
            "classification": "LYNCH SYNDROME (Constitutional MLH1 Epimutation)",
            "confidence": "High",
            "notes": "Constitutional (germline) MLH1 promoter hypermethylation detected in normal tissue baseline, confirming hereditary epimutation.",
        }

    # 3. Rule in Sporadic CRC via BRAF V600E (Only applicable to Colorectal)
    somatic_braf_v600e = patient_data.get("somatic_braf_v600e", False)
    if tissue == "colorectal" and somatic_braf_v600e:
        return {
            "classification": "SPORADIC MMRD (Somatic BRAF Driver)",
            "confidence": "High",
            "notes": "Colorectal tumor with BRAF V600E variant; strongly indicates sporadic MLH1 promoter methylation.",
        }

    # 4. Rule in Sporadic via Epigenetic Silencing (Direct Methylation Testing)
    mlh1_methylated = patient_data.get("somatic_mlh1_promoter_methylation", False)
    if mlh1_methylated:
        mlh1_mut = patient_data.get("somatic_mutation_labels", {}).get("MLH1")
        if mlh1_mut:
            notes_str = f"Biallelic inactivation of MLH1 in {tissue} tumor (somatic mutation {mlh1_mut} and promoter methylation)."
        else:
            notes_str = f"Confirmed epigenetic silencing of MLH1 in {tissue} tumor."
        
        # Clinical guidelines (NCCN / ACG / ACMG): MLH1 promoter methylation under age 50 is atypical
        # and requires clinical correlation or constitutional methylation assay.
        if age != 999 and age < 50:
            notes_str += f" Warning: Somatic MLH1 promoter hypermethylation in a patient under 50 ({age} years) is atypical; clinical correlation and germline/constitutional testing are recommended."
            return {
                "classification": "SPORADIC MMRD (Somatic MLH1 Promoter Hypermethylation)",
                "confidence": "Moderate",
                "notes": notes_str,
            }
        else:
            return {
                "classification": "SPORADIC MMRD (Somatic MLH1 Promoter Hypermethylation)",
                "confidence": "High",
                "notes": notes_str,
            }

    # 5. Evaluate Gene-Specific Somatic Double-Hits (Biallelic Inactivation)
    somatic_variants = patient_data.get("somatic_pathogenic_variants", [])
    loh_status = patient_data.get("loss_of_heterozygosity", {})

    mmr_genes = ["MLH1", "MSH2", "MSH6", "PMS2"]
    double_hit_gene = None

    for gene in mmr_genes:
        # Check if IHC data contradicts or guides us (if IHC is present)
        # E.g., if MSH2 protein is absent, a double hit must occur in MSH2 or its binding partner MSH6
        if has_ihc and valid_ihc.get(gene, "Intact").lower() not in ["absent", "loss", "loss of nuclear expression", "deficient"]:
            continue  # Skip checking double-hits on proteins verified as Intact by IHC

        # Calculate exact hits per individual gene
        mut_count = somatic_variants.count(gene)
        has_loh = loh_status.get(gene, False)

        # Condition A: 2 separate somatic mutations in the exact same gene
        # Condition B: 1 somatic mutation coupled with LOH of that exact same gene
        if mut_count >= 2 or (mut_count == 1 and has_loh):
            double_hit_gene = gene
            break  # Biallelic event found

    if double_hit_gene:
        gene_muts = patient_data.get("somatic_mutation_labels", {}).get(double_hit_gene)
        gene_loh = loh_status.get(double_hit_gene, False)
        if gene_muts and gene_loh:
            notes_str = f"Non-hereditary; biallelic somatic inactivation identified explicitly in {double_hit_gene} ({gene_muts} and LOH)."
        elif gene_muts:
            notes_str = f"Non-hereditary; biallelic somatic inactivation identified explicitly in {double_hit_gene} ({gene_muts})."
        else:
            notes_str = f"Non-hereditary; biallelic somatic inactivation identified explicitly in {double_hit_gene}."
        return {
            "classification": "SOMATIC DOUBLE-HIT (Lynch-Like Syndrome)",
            "confidence": "High",
            "notes": notes_str,
        }

    # 6. Unexplained dMMR: Fallback Logic based on Clinical Risk Corridors
    somatic_genes = list(set(somatic_variants))
    if somatic_genes:
        mut_desc_list = [f"{g} ({patient_data.get('somatic_mutation_labels', {}).get(g, '')})" for g in somatic_genes]
        mut_desc_str = ", ".join(mut_desc_list)
        return {
            "classification": "UNCLASSIFIED dMMR (Single Somatic Hit / Cryptic Second Hit)",
            "confidence": "Moderate",
            "notes": f"Non-hereditary; solitary monoallelic somatic mutation identified in {mut_desc_str} without evidence of a second hit (no LOH, copy loss, or promoter hypermethylation). Consider cryptic non-coding alterations, structural rearrangements, or secondary epigenetic events.",
        }

    if age != 999 and age < 50:
        return {
            "classification": "SUSPECT CONSTITUTIONAL EPIMUTATION / CRYPTIC LYNCH",
            "confidence": "Moderate",
            "notes": "No germline or somatic double-hits found, but young age requires normal tissue MLH1 methylation assay.",
        }

    age_str = f"aged {age}" if age != 999 else "with age unknown"
    return {
        "classification": "UNCLASSIFIED dMMR (Likely Cryptic Sporadic Event)",
        "confidence": "Low",
        "notes": f"Germline negative dMMR tumor in a patient {age_str} without direct evidence of a somatic double hit.",
    }

def check_germline_deletion_bias(purple_dir, group_id, purity_float, genome, cnv_loh, somatic_vafs_dict, high_quality_deletions=None, sample_id=None):
    """Checks purple.segment.tsv for any segment overlapping the genes of interest
    that has a germlineStatus of HET_DELETION or HOM_DELETION.
    If found, calculates the corrected gene-level copy number assuming normal baseline of 2.0.
    Returns a dict: {gene: {status: status, original_cn: cn, corrected_cn: corrected_cn}}
    """
    bias_info = {}
    if purity_float is None or purity_float <= 0.0 or purity_float >= 1.0:
        return bias_info

    segment_file = find_purple_file(purple_dir, "segment.tsv", group_id, sample_id=sample_id)
    if not segment_file or not os.path.exists(segment_file):
        return bias_info

    # Find the genome key
    genome_key = "GRCh38_hmf"
    if genome and "37" in str(genome):
        genome_key = "GRCh37_hmf"

    bounds = GENE_BOUNDS.get(genome_key, {})
    if not bounds:
        return bias_info

    try:
        with open(segment_file, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                chrom = row.get('chromosome')
                start_seg = int(row.get('start', 0))
                end_seg = int(row.get('end', 0))
                germline_status = row.get('germlineStatus', '').strip()

                if germline_status not in ['HET_DELETION', 'HOM_DELETION']:
                    continue

                for gene, info in bounds.items():
                    # Only apply correction if the germline deletion has been validated as high quality
                    if high_quality_deletions is not None and gene not in high_quality_deletions:
                        continue

                    gene_chrom = info['chrom']
                    # Handle chr naming differences
                    if not chrom.startswith('chr') and gene_chrom.startswith('chr'):
                        gene_chrom_clean = gene_chrom.replace('chr', '')
                    elif chrom.startswith('chr') and not gene_chrom.startswith('chr'):
                        gene_chrom_clean = 'chr' + gene_chrom
                    else:
                        gene_chrom_clean = gene_chrom

                    if chrom == gene_chrom_clean:
                        # Check overlap
                        if not (end_seg < info['start'] or start_seg > info['end']):
                            current_status = bias_info.get(gene, {}).get('status')
                            if current_status == 'HOM_DELETION':
                                continue # HOM_DELETION is more severe, keep it
                            bias_info[gene] = {
                                'status': germline_status
                            }

        # Calculate corrected copy number using the gene-level copy number
        for gene, info in list(bias_info.items()):
            gene_cn = cnv_loh.get(gene, {}).get('min_cn', 2.0)
            status = info['status']
            diff = 1.0 if status == 'HET_DELETION' else 2.0
            
            s_vafs = somatic_vafs_dict.get(gene, [])
            has_loh = cnv_loh.get(gene, {}).get('LOH') == 'LOH'
            
            if not has_loh:
                # If there is no LOH, the region is copy-number neutral (diploid normal)
                corrected_cn = 2.0
            else:
                # If there is LOH, it's a somatic copy number loss. We apply somatic VAF checks if variants exist.
                if s_vafs and purity_float is not None and purity_float > 0.0 and purity_float < 1.0:
                    # Adjust expected VAF baseline for the germline deletion in normal cells (2.0 - diff)
                    vaf_diploid = purity_float / (2.0 * purity_float + (2.0 - diff) * (1.0 - purity_float))
                    vaf_loss = purity_float / (1.0 * purity_float + (2.0 - diff) * (1.0 - purity_float))
                    
                    is_diploid_fit = False
                    for vaf in s_vafs:
                        dist_diploid = abs(vaf - vaf_diploid)
                        dist_loss = abs(vaf - vaf_loss)
                        if dist_diploid < dist_loss:
                            is_diploid_fit = True
                            break
                            
                    if is_diploid_fit:
                        corrected_cn = 2.0
                    else:
                        corrected_cn = gene_cn + ((1.0 - purity_float) / purity_float) * diff
                        corrected_cn = max(0.0, corrected_cn)
                else:
                    if purity_float is not None and purity_float > 0.0 and purity_float < 1.0:
                        corrected_cn = gene_cn + ((1.0 - purity_float) / purity_float) * diff
                    else:
                        corrected_cn = gene_cn
                    corrected_cn = max(0.0, corrected_cn)
            
            bias_info[gene] = {
                'status': status,
                'original_cn': max(0.0, gene_cn),
                'corrected_cn': max(0.0, corrected_cn)
            }
    except Exception as e:
        print(f"Warning: Failed to check germline deletion bias: {e}", file=sys.stderr)

    return bias_info

def generate_clinical_report(purple_dir, pave_dir, group_id, subject_id, samplesheet_path, genome, methylation_data, summary_dir, gtf_path=None, ihc_tsv=None, ihc_file=None, cmd_age=None, cmd_tissue=None, ichor_dir=None, min_somatic_vaf=0.0, panel_ichor_dir=None, tumor_sample_id=None, normal_sample_id=None, orange_dir=None, orthogonal_tsv=None, purity_source="purple"):
    """Generate a structured clinical genomics MMR report for the current sample.
    """
    if not ichor_dir and summary_dir:
        parent_out = os.path.dirname(os.path.abspath(summary_dir))
        ichor_dir = resolve_swgs_ichor_dir(None, parent_out, group_id, tumor_sample_id=tumor_sample_id)
    if not panel_ichor_dir and summary_dir:
        parent_out = os.path.dirname(os.path.abspath(summary_dir))
        panel_ichor_dir = resolve_panel_ichor_dir(None, parent_out, group_id, tumor_sample_id=tumor_sample_id)

    metadata = parse_clinical_metadata(ihc_tsv, ihc_file, group_id, subject_id, cmd_age, cmd_tissue)
    ihc_results = metadata['ihc_results']
    age = metadata['age']
    tissue = metadata['tissue']
    orthogonal_data = parse_orthogonal_metadata(orthogonal_tsv, group_id, subject_id) if orthogonal_tsv else None
    purity_file = find_purple_file(purple_dir, "purity.tsv", group_id, sample_id=tumor_sample_id)
    purity, ploidy, tmb, tmb_status, msi, msi_status, wgd = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
    purity_float = None
    ichor_sex = None
    
    if purity_file and os.path.exists(purity_file):
        try:
            with open(purity_file, mode='r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    pur_val = row.get('purity')
                    pl_val = row.get('ploidy')
                    tmb_val = row.get('tmbPerMb')
                    tmb_s = row.get('tmbStatus')
                    msi_val = row.get('msIndelsPerMb')
                    msi_s = row.get('msStatus')
                    wgd_val = row.get('wholeGenomeDuplication', '').lower()
                    
                    if pur_val is not None:
                        purity_float = float(pur_val)
                        purity = f"{purity_float*100:.1f}%"
                    if pl_val is not None:
                        ploidy = f"{float(pl_val):.2f}"
                    if tmb_val is not None:
                        tmb = f"{float(tmb_val):.2f}/Mb"
                    if tmb_s is not None:
                        tmb_status = tmb_s
                    if msi_val is not None:
                        msi = f"{float(msi_val):.2f}/Mb"
                    if msi_s is not None:
                        msi_status = msi_s
                    if wgd_val:
                        wgd = "Yes" if wgd_val in ['true', '1', 'yes'] else "No"
        except Exception as e:
            print(f"Warning: Failed to parse purity file for report: {e}", file=sys.stderr)
            
    # QC stats gathering
    qc_stats = None
    if panel_ichor_dir:
        qc_stats = parse_ichorcna_qc_stats(panel_ichor_dir)
    if not qc_stats:
        qc_stats = get_panel_bam_qc_stats(purple_dir, group_id, samplesheet_path)
        
    swgs_qc_stats = get_swgs_qc_stats(ichor_dir, group_id, subject_id) if ichor_dir else None
        
    # Multi-platform purity extraction and hybrid resolution
    purity_dict = compute_all_purity_sources(purple_dir, ichor_dir, panel_ichor_dir, group_id, subject_id=subject_id, tumor_sample_id=tumor_sample_id, orthogonal_data=orthogonal_data)
    purity_float, ploidy_float, selected_source, purity_notice = resolve_sample_purity(purity_dict, requested_source=purity_source)
    
    if purity_float is not None:
        purity = f"{purity_float*100:.1f}%"
    if ploidy_float is not None:
        ploidy = f"{ploidy_float:.2f}"
    if purity_dict.get('gender') and not ichor_sex:
        ichor_sex = purity_dict['gender']

    # Gather alterations
    cnv_loh = parse_copynumber_and_loh(purple_dir, group_id, sample_id=tumor_sample_id)

    # Recalculate absolute copy numbers ONLY if active purity source is NOT PURPLE
    # When PURPLE is used (default), preserve PURPLE's original calibrated copy numbers.
    if 'purple' not in selected_source and purity_float is not None:
        purple_pur = purity_dict.get('purple')
        if purple_pur is not None and purple_pur > 0.0 and abs(purple_pur - purity_float) > 0.0001:
            try:
                for gene in cnv_loh:
                    raw_cn = cnv_loh[gene]['min_cn']
                    max_cn = cnv_loh[gene]['max_cn']
                    if raw_cn is not None:
                        # Back-calculate raw coverage ratio R and re-scale using active purity
                        R = (raw_cn * purple_pur) / 2.0 + 1.0 - purple_pur
                        new_cn = 2.0 * (R - 1.0) / purity_float + 2.0
                        cnv_loh[gene]['min_cn'] = new_cn
                        
                        if max_cn is not None:
                            R_max = (max_cn * purple_pur) / 2.0 + 1.0 - purple_pur
                            new_max = 2.0 * (R_max - 1.0) / purity_float + 2.0
                            cnv_loh[gene]['max_cn'] = new_max
                        else:
                            new_max = new_cn
                            
                        # Re-classify and update copynumber status string based on un-capped values
                        if new_cn < 0.5:
                            cnv_loh[gene]['copynumber'] = 'Homozygous Deletion'
                        elif new_cn < 1.7:
                            cnv_loh[gene]['copynumber'] = 'Heterozygous Loss'
                        elif new_max > 2.8:
                            cnv_loh[gene]['copynumber'] = f'Gain (CN={max(0.0, new_max):.2f})'
                        else:
                            cnv_loh[gene]['copynumber'] = 'Normal'
            except Exception as e:
                print(f"Warning: Failed to re-scale copy numbers for {selected_source}: {e}", file=sys.stderr)
    
    # Pre-parse high-quality germline deletions to filter segment bias corrections
    germline_deletions = parse_germline_deletions(purple_dir, group_id, genome=genome, sample_id=tumor_sample_id)
    high_quality_deletions = germline_deletions
    
    # Pre-gather somatic variant VAFs for the bias check
    somatic_vafs_dict = {}
    for gene in GENES_OF_INTEREST:
        vars_list = get_detailed_variants(purple_dir, pave_dir, group_id, gene, subject_id, samplesheet_path, genome, methylation_data, min_somatic_vaf=min_somatic_vaf, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id)
        somatic_vafs = [v.get('vaf') for v in vars_list if v.get('var_type') == 'somatic' and v.get('vaf') is not None]
        somatic_vafs_dict[gene] = somatic_vafs
        
    bias_dict = check_germline_deletion_bias(purple_dir, group_id, purity_float, genome, cnv_loh, somatic_vafs_dict, high_quality_deletions=None, sample_id=tumor_sample_id)
    
    # Cap all cnv_loh copy numbers at 0.0 now that the bias check has completed
    for gene in cnv_loh:
        if cnv_loh[gene]['min_cn'] is not None:
            cnv_loh[gene]['min_cn'] = max(0.0, cnv_loh[gene]['min_cn'])
        if cnv_loh[gene]['max_cn'] is not None:
            cnv_loh[gene]['max_cn'] = max(0.0, cnv_loh[gene]['max_cn'])
    
    # Fallback to ichorCNA copy number segments to fill gaps (e.g. PMS2 mapping issues)
    ichor_segments = parse_ichorcna_segments(ichor_dir, group_id, subject_id, genome)
    for gene in GENES_OF_INTEREST:
        if cnv_loh[gene]['copynumber'] is None and gene in ichor_segments and ichor_segments[gene] is not None:
            ichor_cn = ichor_segments[gene]
            if ichor_cn < 0.5:
                cnv_loh[gene]['copynumber'] = 'Homozygous Deletion'
                cnv_loh[gene]['is_ichor_cn'] = True
            elif ichor_cn <= 1.7:
                cnv_loh[gene]['copynumber'] = 'Heterozygous Loss'
                cnv_loh[gene]['is_ichor_cn'] = True
            elif ichor_cn >= 2.8:
                cnv_loh[gene]['copynumber'] = f'Gain (CN={ichor_cn:.2f})'
                cnv_loh[gene]['is_ichor_cn'] = True
            else:
                cnv_loh[gene]['copynumber'] = 'Normal'
                cnv_loh[gene]['is_ichor_cn'] = True

    somatic_muts = parse_somatic_catalog(purple_dir, group_id, sample_id=tumor_sample_id)
    germline_muts = parse_pave_vcf(pave_dir, group_id, 'germline', sample_id=tumor_sample_id)
    methylation = parse_methylation(methylation_data, group_id, subject_id, samplesheet_path, genome)
    
    # Exon-disrupting SV parsing
    exons_dict = {}
    if gtf_path and os.path.exists(gtf_path):
        exons_dict = load_exon_coordinates(gtf_path, GENES_OF_INTEREST)
    somatic_sv_file = find_purple_file(purple_dir, "sv.vcf.gz", group_id, sample_id=tumor_sample_id)
    germline_sv_file = find_purple_file(purple_dir, "sv.germline.vcf.gz", group_id, sample_id=tumor_sample_id)
    somatic_svs = parse_pathogenic_svs(somatic_sv_file, exons_dict, genome) if (somatic_sv_file and exons_dict) else {g: [] for g in GENES_OF_INTEREST}
    germline_svs = parse_pathogenic_svs(germline_sv_file, exons_dict, genome) if (germline_sv_file and exons_dict) else {g: [] for g in GENES_OF_INTEREST}
    
    gene_reports = {}
    for gene in GENES_OF_INTEREST:
        vars_list = get_detailed_variants(purple_dir, pave_dir, group_id, gene, subject_id, samplesheet_path, genome, methylation_data, min_somatic_vaf=min_somatic_vaf, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id)
        
        s_muts = []
        g_muts = []
        
        is_ichor = cnv_loh[gene].get('is_ichor_cn', False)
        s_cn = cnv_loh[gene]['copynumber'] or 'Normal'
        
        bias = bias_dict.get(gene)
        if bias and not is_ichor:
            corrected_val = bias['corrected_cn']
            if corrected_val < 0.5:
                corrected_status = 'Homozygous Deletion'
            elif corrected_val < 1.7:
                corrected_status = 'Heterozygous Loss'
            elif corrected_val > 2.8:
                corrected_status = f'Gain (CN={corrected_val:.2f})'
            else:
                corrected_status = 'Normal'
            
            cnv_loh[gene]['copynumber'] = corrected_status
            s_cn = corrected_status
            s_cn_display = f"{corrected_status} (PURPLE CN: {bias['original_cn']:.2f}, Corrected CN: {bias['corrected_cn']:.2f} due to germline {bias['status']})"
        else:
            s_cn_display = f"{s_cn} (ichorCNA)" if (s_cn != 'Normal' and is_ichor) else s_cn
            
        g_cn = 'Deletion' if (gene in germline_deletions) else 'Normal'
        loh = 'LOH' if cnv_loh[gene]['LOH'] else 'Normal'
        meth_status = 'No Data'
        
        if gene in methylation:
            m_info = methylation[gene]
            status = m_info['status']
            t_beta = m_info['tumor_beta']
            n_beta = m_info['normal_beta']
            if status == 'unmethylated':
                meth_status = f"Promoter Unmethylated (T={t_beta:.2f}"
            else:
                meth_status = f"Promoter Methylated ({status.capitalize()}, T={t_beta:.2f}"
            if n_beta is not None:
                meth_status += f", N={n_beta:.2f}"
            meth_status += ")"
            
        # Add SVs to s_muts and g_muts list, and add them to vars_list for detailed annotations/BED mapping
        import re
        for sv_desc in somatic_svs.get(gene, []):
            s_muts.append(sv_desc)
            pos_val = 0
            m_pos = re.search(r':(\d+)', sv_desc)
            if m_pos:
                pos_val = int(m_pos.group(1))
            vars_list.append({
                'pos': pos_val,
                'label': f"Somatic SV: {sv_desc}",
                'phasing': 'N/A (structural variant)',
                'vaf': None,
                'lps': None,
                'lps_rc': None,
                'var_type': 'somatic_sv'
            })
            
        for sv_desc in germline_svs.get(gene, []):
            g_muts.append(sv_desc)
            pos_val = 0
            m_pos = re.search(r':(\d+)', sv_desc)
            if m_pos:
                pos_val = int(m_pos.group(1))
            vars_list.append({
                'pos': pos_val,
                'label': f"Germline SV: {sv_desc}",
                'phasing': 'N/A (structural variant)',
                'vaf': None,
                'lps': None,
                'lps_rc': None,
                'var_type': 'germline_sv'
            })

        for v in vars_list:
            v_type = v.get('var_type')
            label = v.get('label', '')
            if v_type in ['somatic_sv', 'germline_sv']:
                continue # Already handled above
            # Extract mutation change details
            # E.g. "Somatic MLH1 p.Arg123Ter (VAF: 45%)" -> "p.Arg123Ter (VAF: 45%)"
            if label.startswith("Somatic ") or label.startswith("Germline "):
                parts = label.split(' ', 2)
                change_desc = parts[2] if len(parts) > 2 else label
            else:
                change_desc = label
            
            if v_type == 'somatic':
                s_muts.append(change_desc)
            elif v_type == 'germline':
                g_muts.append(change_desc)
                
        s_mut_str = ", ".join(s_muts) if s_muts else 'Wildtype'
        g_mut_str = ", ".join(g_muts) if g_muts else 'Wildtype'
        
        # Risk assessment
        has_loh = loh == 'LOH'
        has_somatic_mut = len(s_muts) > 0
        has_germline_mut = len(g_muts) > 0
        has_methylation = (gene in methylation) and (methylation[gene]['status'] in ['somatic', 'germline'])
        has_germline_del = g_cn == 'Deletion'
        
        risk_status = "No Alterations"
        risk_score = 0
        
        # Check specific levels of loss / gain
        has_hom_del = 'homozygous deletion' in s_cn_display.lower()
        has_het_loss = 'heterozygous loss' in s_cn_display.lower()
        has_gain = 'gain' in s_cn_display.lower()
        
        if has_hom_del or (has_loh and (has_somatic_mut or has_germline_mut or has_methylation or has_germline_del)):
            risk_status = "HIGH RISK - Biallelic Inactivation"
            risk_score = 3
        elif len(s_muts) >= 2 or (has_somatic_mut and has_germline_mut):
            risk_status = "HIGH RISK - Compound Heterozygous"
            risk_score = 3
        elif has_somatic_mut or has_germline_mut or has_methylation or has_germline_del:
            risk_status = "MODERATE RISK - Monoallelic Alteration"
            risk_score = 2
        elif has_het_loss:
            risk_status = "MODERATE RISK - Somatic Copy Number Loss"
            risk_score = 2
        elif has_gain:
            risk_status = "LOW RISK - Somatic Copy Number Gain"
            risk_score = 1
            
        gene_reports[gene] = {
            's_mut': s_mut_str,
            'g_mut': g_mut_str,
            's_muts': s_muts,
            'g_muts': g_muts,
            's_cn': s_cn_display,
            'g_cn': g_cn,
            'loh': loh,
            'meth': meth_status,
            'risk': risk_status,
            'score': risk_score,
            'vars_list': vars_list
        }
        
    # Executive summary status recommendations
    mmr_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM', 'MLH3']
    dmmr_genes = []
    mono_mmr_genes = []
    
    for g in mmr_genes:
        rep = gene_reports[g]
        if rep['score'] == 3:
            dmmr_genes.append(g)
        elif rep['score'] == 2:
            mono_mmr_genes.append(g)
             # 4. Run classification engine
    # Prepare patient_data
    patient_data = {
        'msi_status': msi_status,
        'tmb': float(tmb.replace('/Mb', '')) if tmb != 'NA' else 0.0,
        'tumor_tissue_origin': tissue,
        'age_at_diagnosis': age,
        'germline_pathogenic_variants': [],
        'somatic_braf_v600e': False,
        'somatic_mlh1_promoter_methylation': False,
        'germline_mlh1_promoter_methylation': False,
        'somatic_pathogenic_variants': [],
        'loss_of_heterozygosity': {},
        'ihc_results': ihc_results,
        'somatic_mutation_labels': {}
    }
    
    for gene in ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'EPCAM']:
        rep = gene_reports.get(gene)
        if rep:
            if rep['g_mut'] != 'Wildtype':
                patient_data['germline_pathogenic_variants'].append(gene)
            elif rep['g_cn'] == 'Deletion':
                patient_data['germline_pathogenic_variants'].append(gene)
                
    braf_rep = gene_reports.get('BRAF')
    if braf_rep and braf_rep['s_mut'] != 'Wildtype':
        if 'V600E' in braf_rep['s_mut']:
            patient_data['somatic_braf_v600e'] = True
            
    mlh1_rep = gene_reports.get('MLH1')
    if mlh1_rep and mlh1_rep['meth'] != 'No Data':
        if 'promoter methylated' in mlh1_rep['meth'].lower():
            if 'germline' in mlh1_rep['meth'].lower():
                patient_data['germline_mlh1_promoter_methylation'] = True
            else:
                patient_data['somatic_mlh1_promoter_methylation'] = True
            
    for gene in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
        rep = gene_reports.get(gene)
        if rep:
            for _ in rep.get('s_muts', []):
                patient_data['somatic_pathogenic_variants'].append(gene)
            if rep.get('s_muts'):
                patient_data['somatic_mutation_labels'][gene] = rep['s_mut']
            patient_data['loss_of_heterozygosity'][gene] = (rep['loh'] == 'LOH')
            
    clinical_info = classify_lynch_vs_sporadic(patient_data)
    clinical_info['qc_stats'] = qc_stats
    clinical_info['swgs_qc_stats'] = swgs_qc_stats
    mean_diploid_vaf, mean_all_vaf = compute_driver_vaf_metrics(purple_dir, group_id, sample_id=tumor_sample_id)
    clinical_info['mean_driver_vaf_diploid'] = mean_diploid_vaf
    clinical_info['mean_driver_vaf_all'] = mean_all_vaf
    bias_notes = []
    for g, b in bias_dict.items():
        if b['corrected_cn'] < 1.7:
            bias_notes.append(f"Note: {g} has a germline {b['status']} call in the PURPLE baseline, distorting absolute tumor CN to {b['original_cn']:.2f}. Corrected CN (diploid normal baseline) is {b['corrected_cn']:.2f}.")
    if bias_notes:
        clinical_info['notes'] += " | " + " ".join(bias_notes)
        
    # Apply purity advisory and low tumor purity QC warnings
    if purity_notice:
        clinical_info['notes'] = f"Purity Advisory: {purity_notice} | " + clinical_info['notes']

    if purity_float is not None and purity_float < 0.3:
        clinical_info['confidence'] = "Low (Purity < 30%)"
        warning_msg = f"WARNING: Low tumor purity ({purity_float*100:.1f}%) is below the 30% quality control threshold. Copy number profiling and somatic calls have low confidence and may be artifactual."
        clinical_info['notes'] = warning_msg + " | " + clinical_info['notes']
        
    # Generate Markdown content
    md = []
    md.append(f"# CLINICAL GENOMICS REPORT: MISMATCH REPAIR (MMR) STATUS")
    md.append(f"**Generated on**: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    p_pur_str = f"{purity_dict['purple']*100:.1f}%" if purity_dict.get('purple') is not None else "NA"
    s_pur_str = f"{purity_dict['swgs']*100:.1f}%" if purity_dict.get('swgs') is not None else "NA"
    pan_pur_str = f"{purity_dict['panel_ichor']*100:.1f}%" if purity_dict.get('panel_ichor') is not None else "NA"
    vaf_pur_str = f"{purity_dict['vaf']*100:.1f}%" if purity_dict.get('vaf') is not None else "NA"
    path_val = purity_dict.get('pathology')
    if path_val is None and orthogonal_data and orthogonal_data.get('pathology_purity') is not None:
        path_val = orthogonal_data['pathology_purity']
    path_pur_str = f"{path_val*100:.0f}%" if path_val is not None else "NA"

    md.append(f"## 1. Patient & Sample Information")
    md.append(f"| Attribute | Value |")
    md.append(f"| :--- | :--- |")
    md.append(f"| **Patient ID** | {subject_id} |")
    md.append(f"| **Sample ID** | {group_id} |")
    if tumor_sample_id and tumor_sample_id != group_id:
        md.append(f"| **Tumor Biospecimen ID** | {tumor_sample_id} |")
    if ichor_sex:
        md.append(f"| **Determined Sex** | {ichor_sex} |")
    md.append(f"| **Tumor Purity (Active)** | **{purity}** (Source: {selected_source.upper()}) |")
    pur_platforms = [f"PURPLE: {p_pur_str}", f"sWGS: {s_pur_str}", f"Panel ichorCNA: {pan_pur_str}", f"Clonal VAF Implied: {vaf_pur_str}"]
    if path_pur_str != "NA":
        pur_platforms.append(f"Pathology: {path_pur_str}")
    pur_summary_str = " | ".join(pur_platforms)
    md.append(f"| **Purity Estimates (All Platforms)** | {pur_summary_str} |")
    if path_pur_str != "NA":
        md.append(f"| **Pathology Tumor Cellularity** | {path_pur_str} |")
    md.append(f"| **Tumor Ploidy** | {ploidy} |")
    md.append(f"| **Tumor Mutational Burden (TMB)** | {tmb} ({tmb_status}) |")
    md.append(f"| **Microsatellite Instability (MSI)** | {msi} ({msi_status}) |")
    md.append(f"| **Whole Genome Duplication (WGD)** | {wgd} |")
    md.append("")
    
    md.append(f"## 2. Executive Summary & Clinical Assessment")
    if purity_notice:
        md.append(f"> [!WARNING]")
        md.append(f"> **Tumor Purity Advisory**: {purity_notice}")
        md.append(f">")
    if purity_float is not None and purity_float < 0.3:
        md.append(f"> [!WARNING]")
        md.append(f"> **Low Tumor Purity Warning ({purity_float*100:.1f}%)**: The tumor fraction is below the 30% quality control threshold. Copy number profiling and somatic calls have low confidence and may be artifactual.")
        md.append(f">")
    md.append(f"> [!IMPORTANT]")
    md.append(f"> **Clinical Classification**: **{clinical_info['classification']}**")
    md.append(f">")
    md.append(f"> **Confidence Level**: {clinical_info['confidence']}")
    md.append(f">")
    md.append(f"> **Classification Evidence Details**:")
    for note in clinical_info['notes'].split(' | '):
        md.append(f"> * {note}")
    md.append("")
    
    md.append(f"## 3. Mismatch Repair (MMR) Core Genes")
    md.append(f"| Gene | Somatic Mutation | Germline Mutation | Somatic CN | Germline CN | LOH | Promoter Methylation |")
    md.append(f"| :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
    for g in mmr_genes:
        rep = gene_reports[g]
        md.append(f"| **{g}** | {rep['s_mut']} | {rep['g_mut']} | {rep['s_cn']} | {rep['g_cn']} | {rep['loh']} | {rep['meth']} |")
    md.append("")
    
    has_ihc = False
    if ihc_results:
        has_ihc = any(v not in ["No Data", "Not Evaluated"] for v in ihc_results.values())
        
    next_sec = 4
    if has_ihc:
        md.append(f"## {next_sec}. Immunohistochemistry (IHC) & Genomic Correlation")
        md.append(f"| Gene | IHC Clinical Result | Genomic / Epigenomic Status | Concordance |")
        md.append(f"| :--- | :--- | :--- | :--- |")
        for g in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
            res = evaluate_ihc_concordance(g, ihc_results, gene_reports, purity_float)
            md.append(f"| **{g}** | {ihc_results.get(g, 'No Data')} | {res['genomic_status']} | {res['concordance']} |")
        md.append("")
        next_sec += 1

    if orthogonal_data and (orthogonal_data.get('variants') or orthogonal_data.get('test_result')):
        md.append(f"## {next_sec}. Orthogonal Assay & Known Mutation Concordance")
        if orthogonal_data.get('test_result'):
            md.append(f"**Overall Orthogonal Result**: {orthogonal_data['test_result']}")
        if orthogonal_data.get('pathology_purity') is not None:
            md.append(f"**Pathology Tumor Cellularity**: {int(orthogonal_data['pathology_purity']*100)}%")
        md.append("")
        md.append(f"| Gene | Orthogonal Finding | Assay Interpretation | Oncoanalyser Genomic Finding | Concordance |")
        md.append(f"| :--- | :--- | :--- | :--- | :--- |")
        orth_genes = list(dict.fromkeys([v['gene'] for v in orthogonal_data.get('variants', []) if v['gene'] not in ['None', 'MMR']]))
        if not orth_genes:
            orth_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2']
        else:
            for cg in ['MLH1', 'MSH2', 'MSH6', 'PMS2']:
                if cg not in orth_genes:
                    orth_genes.append(cg)
        for g in orth_genes:
            c_res = evaluate_orthogonal_concordance(g, orthogonal_data, gene_reports)
            if c_res.get('has_data'):
                md.append(f"| **{g}** | {c_res['orthogonal_findings']} | {c_res.get('clinical_significance', 'Unknown')} | {c_res['onco_findings']} | **{c_res['concordance']}** |")
        md.append("")
        next_sec += 1
        
    md.append(f"## {next_sec}. Other Driver & Ancillary Genes")
    md.append(f"| Gene | Somatic Mutation | Germline Mutation | Somatic CN | Germline CN | LOH | Promoter Methylation |")
    md.append(f"| :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
    other_genes = ['POLD1', 'POLE', 'MUTYH', 'TP53', 'BRAF']
    for g in other_genes:
        rep = gene_reports[g]
        md.append(f"| **{g}** | {rep['s_mut']} | {rep['g_mut']} | {rep['s_cn']} | {rep['g_cn']} | {rep['loh']} | {rep['meth']} |")
    md.append("")
    next_sec += 1
    
    md.append(f"## {next_sec}. Detailed Variant Phasing & Annotations")
    has_vars = False
    for g in GENES_OF_INTEREST:
        rep = gene_reports[g]
        if rep['vars_list']:
            has_vars = True
            md.append(f"### {g}")
            for v in rep['vars_list']:
                md.append(f"- **Variant**: `{v['label']}`")
                md.append(f"  - Phasing status: {v['phasing']}")
                if v.get('lps') is not None:
                    md.append(f"  - Phase Set Support: ID={v['lps']}, read count={v['lps_rc']}")
    if not has_vars:
        md.append("No pathogenic mutations or promoter methylation events detected.")
        
    # Write report file
    report_file = os.path.join(summary_dir, f"{group_id}_clinical_mmr_report.md")
    try:
        with open(report_file, 'w', encoding='utf-8') as f_out:
            f_out.write("\n".join(md))
        print(f"Clinical MMR report successfully written to {report_file}")
    except Exception as e:
        print(f"Error: Failed to write clinical report to {report_file}: {e}", file=sys.stderr)
 
    # Write clinical HTML report
    try:
        generate_clinical_html_report(
            group_id=group_id,
            subject_id=subject_id,
            purity=purity,
            ploidy=ploidy,
            tmb=tmb,
            tmb_status=tmb_status,
            msi=msi,
            msi_status=msi_status,
            wgd=wgd,
            clinical_info=clinical_info,
            gene_reports=gene_reports,
            summary_dir=summary_dir,
            ihc_results=ihc_results,
            purity_float=purity_float,
            purple_dir=purple_dir,
            ichor_dir=ichor_dir,
            ichor_sex=ichor_sex,
            panel_ichor_dir=panel_ichor_dir,
            qc_stats=qc_stats,
            swgs_qc_stats=swgs_qc_stats,
            orange_dir=orange_dir,
            tumor_sample_id=tumor_sample_id,
            orthogonal_data=orthogonal_data,
            purity_dict=purity_dict,
            selected_source=selected_source,
            purity_notice=purity_notice
        )
    except Exception as e:
        print(f"Error: Failed to generate HTML report for {group_id}: {e}", file=sys.stderr)

    return clinical_info

def compute_driver_vaf_metrics(purple_dir, group_id, sample_id=None):
    """
    Parses PURPLE driver catalog (driver.catalog.somatic.tsv)
    and segments (segment.tsv) to compute:
      - mean driver VAF: diploid (in regions with CN=2)
      - mean driver VAF: all
    """
    if not purple_dir or not os.path.exists(purple_dir):
        return "NA", "NA"
        
    catalog_file = find_purple_file(purple_dir, "driver.catalog.somatic.tsv", group_id, sample_id=sample_id)
    segment_file = find_purple_file(purple_dir, "segment.tsv", group_id, sample_id=sample_id)
    
    # 1. Parse diploid segments
    diploid_segments = []
    if segment_file and os.path.exists(segment_file):
        try:
            with open(segment_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    chr_name = row.get('chromosome') or row.get('chr')
                    try:
                        start = int(row.get('start', 0))
                        end = int(row.get('end', 0))
                        cn = float(row.get('copyNumber', 0))
                        if 1.5 <= cn <= 2.5:
                            diploid_segments.append((chr_name, start, end))
                    except ValueError:
                        pass
        except Exception:
            pass

    def is_diploid(chr_n, pos):
        for s_chr, s_start, s_end in diploid_segments:
            if s_chr == chr_n and s_start <= pos <= s_end:
                return True
        return False

    vafs_all = []
    vafs_diploid = []
    
    # 2. Parse drivers from catalog
    if catalog_file and os.path.exists(catalog_file):
        try:
            with open(catalog_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    driver_type = row.get('driver', '').strip()
                    driver_lik = row.get('driverLikelihood', '0')
                    try:
                        lik_val = float(driver_lik)
                    except ValueError:
                        lik_val = 0.0
                    
                    if (driver_type and driver_type not in ['', '-', 'NONE']) or lik_val > 0:
                        vaf_val = row.get('vaf') or row.get('tumorVaf')
                        if vaf_val is not None:
                            try:
                                vaf_f = float(vaf_val)
                                chr_n = row.get('chromosome') or row.get('chr')
                                pos = int(row.get('start', 0))
                                
                                vafs_all.append(vaf_f)
                                if is_diploid(chr_n, pos):
                                    vafs_diploid.append(vaf_f)
                            except ValueError:
                                pass
        except Exception as e:
            print(f"Warning: Error parsing driver catalog for VAF: {e}", file=sys.stderr)
            
    # 3. Fallback to purple.somatic.vcf.gz if catalog is missing/empty
    if not vafs_all:
        somatic_vcf = find_purple_file(purple_dir, "somatic.vcf.gz", group_id, sample_id=sample_id)
        if os.path.exists(somatic_vcf):
            try:
                import gzip
                with gzip.open(somatic_vcf, 'rt') as f_vcf:
                    for line in f_vcf:
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:
                            info = parts[7]
                            if 'PURPLE_DRIVER' in info or 'REPORTED' in info:
                                chr_n = parts[0]
                                pos = int(parts[1])
                                vaf_f = None
                                for item in info.split(';'):
                                    if item.startswith('VAF='):
                                        try:
                                            vaf_f = float(item.split('=')[1])
                                        except ValueError:
                                            pass
                                if vaf_f is not None:
                                    vafs_all.append(vaf_f)
                                    if is_diploid(chr_n, pos):
                                        vafs_diploid.append(vaf_f)
            except Exception:
                pass
                
    mean_all_str = f"{sum(vafs_all)/len(vafs_all):.4f}" if vafs_all else "NA"
    mean_diploid_str = f"{sum(vafs_diploid)/len(vafs_diploid):.4f}" if vafs_diploid else "NA"
    
    return mean_diploid_str, mean_all_str

def main():
    parser = argparse.ArgumentParser(description='Summarize findings for MMR and cancer genes from oncoanalyser outputs.')
    parser.add_argument('--output-dir', required=True, help='Path to run output directory')
    parser.add_argument('--samplesheet', required=True, help='Path to Nextflow samplesheet CSV')
    parser.add_argument('--methylation-data', help='Path to external methylation gene data TSV')
    parser.add_argument('--project', default='oncoanalyser', help='Project name/prefix for the final summary file')
    parser.add_argument('--genome', default='GRCh38_hmf', help='Genome build (GRCh38_hmf or GRCh37_hmf)')
    parser.add_argument('--fasta', help='Path to reference genome FASTA file')
    parser.add_argument('--gtf', help='Path to transcript GTF/GFF file for gene models')
    parser.add_argument('--igv-reports-env', default='igvreports_env', help='Conda environment for igv-reports')
    parser.add_argument('--ihc-tsv', help='Path to cohort-wide IHC status TSV mapping file')
    parser.add_argument('--ihc-file', help='Path to single-sample key-value IHC text config file')
    parser.add_argument('--age', type=int, help='Patient age at diagnosis')
    parser.add_argument('--tissue', help='Tumor tissue origin (e.g. colorectal)')
    parser.add_argument('--ichor-dir', help='Path to ichorCNA output directory for ploidy/purity fallbacks and CNV gap filling')
    parser.add_argument('--panel-ichor-dir', help='Path to panel off-target ichorCNA directory')
    parser.add_argument('--min-somatic-vaf', type=float, default=0.05, help='Minimum somatic VAF filter threshold (default: 0.05)')
    parser.add_argument('--orthogonal-tsv', help='Path to orthogonal test / known mutations TSV mapping file')
    parser.add_argument('--purity-source', default='purple', choices=['purple', 'auto', 'swgs', 'panel', 'panel_ichor', 'vaf', 'pathology'], help='Purity estimation source (default: purple). Options: purple, auto, swgs, panel/panel_ichor, vaf, pathology')
    
    args = parser.parse_args()
    
    if not args.fasta:
        default_fasta = "/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/2.3.0/GRCh38_masked_exclusions_alts_hlas.fasta"
        if os.path.exists(default_fasta):
            args.fasta = default_fasta
    
    # 1. Parse samplesheet to find all tumor samples
    samples = parse_samplesheet(args.samplesheet)
    if not samples:
        print("No tumor samples found in the samplesheet.", file=sys.stderr)
        sys.exit(0)
        
    # 2. Build rows of the summary table
    summary_rows = []
    
    for sample in samples:
        group_id = sample['group_id']
        subject_id = sample['subject_id']
        
        # Resolve sample IDs from samplesheet
        tumor_sample_id, normal_sample_id = get_sample_ids_from_samplesheet(args.samplesheet, group_id)
        if not tumor_sample_id:
            tumor_sample_id = sample.get('sample_id')
        sample_key = tumor_sample_id or group_id
        
        # Check potential purple output directories:
        # Typically <output_dir>/<group_id>/purple/
        purple_dir = os.path.join(args.output_dir, group_id, 'purple')
        if not os.path.exists(purple_dir):
            # Try fallback to <output_dir>/purple
            purple_dir = os.path.join(args.output_dir, 'purple')

        pave_dir = os.path.join(args.output_dir, group_id, 'pave')
        if not os.path.exists(pave_dir):
            pave_dir = os.path.join(args.output_dir, 'pave')

        orange_dir = os.path.join(args.output_dir, group_id, 'orange')
        if not os.path.exists(orange_dir):
            orange_dir = os.path.join(args.output_dir, 'orange')
            
        if not os.path.exists(purple_dir):
            print(f"Warning: purple directory not found for sample {group_id} under {args.output_dir}", file=sys.stderr)
            # We can still add a row with NA if files aren't found yet
            row_data = {
                'Patient': subject_id or group_id,
                'Sample': sample_key,
                'Tumour.Fraction': 'NA',
                'Purity.Source': 'NA',
                'Purity.PURPLE': 'NA',
                'Purity.sWGS': 'NA',
                'Purity.Panel_ichor': 'NA',
                'Purity.VAF_Implied': 'NA',
                'MSI.Status': 'NA',
                'MSI.Indels.Per.Mb': 'NA',
                'MSI.Unstable.Sites.Proportion': 'NA',
                'Clinical.Classification': 'NA',
                'Clinical.Confidence': 'NA',
                'Clinical.Notes': 'NA',
                'Total.Reads': 'NA',
                'On.Target.Reads': 'NA',
                'On.Target.Coverage': 'NA',
                'Off.Target.Reads': 'NA',
                'Off.Target.Coverage': 'NA',
                'Panel.ichorCNA.MAD': 'NA',
                'sWGS.Coverage': 'NA',
                'sWGS.ichorCNA.MAD': 'NA'
            }
            for gene in GENES_OF_INTEREST:
                row_data[gene] = ''
            summary_rows.append(row_data)
            continue
            
        # Resolve sWGS and panel ichorCNA directories
        ichor_dir = resolve_swgs_ichor_dir(args.ichor_dir, args.output_dir, group_id, tumor_sample_id=tumor_sample_id)
        panel_ichor_dir = resolve_panel_ichor_dir(args.panel_ichor_dir, args.output_dir, group_id, tumor_sample_id=tumor_sample_id)

        # Resolve sample purity across platforms
        orthogonal_data = parse_orthogonal_metadata(args.orthogonal_tsv, group_id, subject_id) if args.orthogonal_tsv else None
        purity_dict = compute_all_purity_sources(purple_dir, ichor_dir, panel_ichor_dir, group_id, subject_id=subject_id, tumor_sample_id=tumor_sample_id, orthogonal_data=orthogonal_data)
        purity_float, ploidy_float, selected_source, purity_notice = resolve_sample_purity(purity_dict, requested_source=args.purity_source)
        
        purity = str(round(purity_float, 5)) if purity_float is not None else 'NA'
        msi_status = get_msi_status(purple_dir, group_id, sample_id=tumor_sample_id)
        cnv_loh = parse_copynumber_and_loh(purple_dir, group_id, sample_id=tumor_sample_id)
        
        # Apply copy number re-scaling ONLY if the selected purity source is NOT PURPLE
        # When PURPLE is used (default), preserve PURPLE's original calibrated copy numbers.
        if 'purple' not in selected_source and purity_float is not None:
            purple_pur = purity_dict.get('purple')
            if purple_pur is not None and purple_pur > 0.0 and abs(purple_pur - purity_float) > 0.0001:
                try:
                    for gene in cnv_loh:
                        raw_cn = cnv_loh[gene]['min_cn']
                        max_cn = cnv_loh[gene]['max_cn']
                        if raw_cn is not None:
                            # Back-calculate raw coverage ratio R
                            R = (raw_cn * purple_pur) / 2.0 + 1.0 - purple_pur
                            new_cn = 2.0 * (R - 1.0) / purity_float + 2.0
                            cnv_loh[gene]['min_cn'] = new_cn
                            
                            if max_cn is not None:
                                R_max = (max_cn * purple_pur) / 2.0 + 1.0 - purple_pur
                                new_max = 2.0 * (R_max - 1.0) / purity_float + 2.0
                                cnv_loh[gene]['max_cn'] = new_max
                            else:
                                new_max = new_cn
                                
                            # Re-classify and update copynumber status string based on un-capped values
                            if new_cn < 0.5:
                                cnv_loh[gene]['copynumber'] = 'Homozygous Deletion'
                            elif new_cn < 1.7:
                                cnv_loh[gene]['copynumber'] = 'Heterozygous Loss'
                            elif new_max > 2.8:
                                cnv_loh[gene]['copynumber'] = f'Gain (CN={max(0.0, new_max):.2f})'
                            else:
                                cnv_loh[gene]['copynumber'] = 'Normal'
                except Exception:
                    pass

        # Find PAVE directory
        pave_dir = os.path.join(args.output_dir, group_id, 'pave')
        if not os.path.exists(pave_dir):
            pave_dir = os.path.join(args.output_dir, 'pave')
            
        somatic_muts = parse_somatic_catalog(purple_dir, group_id, sample_id=tumor_sample_id)
        germline_muts = parse_pave_vcf(pave_dir, group_id, 'germline', sample_id=tumor_sample_id)
        germline_deletions = parse_germline_deletions(purple_dir, group_id, args.genome, sample_id=tumor_sample_id)
        germline_catalog = parse_germline_catalog(purple_dir, group_id, sample_id=tumor_sample_id)
        sample_meth_data = find_methylation_data(args.methylation_data, args.output_dir, group_id)
        methylation = parse_methylation(sample_meth_data, group_id, subject_id, args.samplesheet, args.genome)

        somatic_vafs_dict = {}
        for gene in GENES_OF_INTEREST:
            vars_list = get_detailed_variants(purple_dir, pave_dir, group_id, gene, subject_id, args.samplesheet, args.genome, sample_meth_data, min_somatic_vaf=args.min_somatic_vaf, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id)
            somatic_vafs = [v.get('vaf') for v in vars_list if v.get('var_type') == 'somatic' and v.get('vaf') is not None]
            somatic_vafs_dict[gene] = somatic_vafs

        high_quality_deletions = germline_deletions
        bias_dict = check_germline_deletion_bias(purple_dir, group_id, purity_float, args.genome, cnv_loh, somatic_vafs_dict, high_quality_deletions=None, sample_id=tumor_sample_id)
        
        # Cap all cnv_loh copy numbers at 0.0 now that the bias check has completed
        for gene in cnv_loh:
            if cnv_loh[gene]['min_cn'] is not None:
                cnv_loh[gene]['min_cn'] = max(0.0, cnv_loh[gene]['min_cn'])
            if cnv_loh[gene]['max_cn'] is not None:
                cnv_loh[gene]['max_cn'] = max(0.0, cnv_loh[gene]['max_cn'])
        
        # Fallback to ichorCNA copy number segments to fill gaps (e.g. PMS2 mapping issues)
        if ichor_dir:
            ichor_segments = parse_ichorcna_segments(ichor_dir, group_id, subject_id, args.genome)
            for gene in GENES_OF_INTEREST:
                if cnv_loh[gene]['copynumber'] is None and gene in ichor_segments and ichor_segments[gene] is not None:
                    ichor_cn = ichor_segments[gene]
                    if ichor_cn < 0.5:
                        cnv_loh[gene]['copynumber'] = 'Homozygous Deletion'
                        cnv_loh[gene]['is_ichor_cn'] = True
                    elif ichor_cn <= 1.7:
                        cnv_loh[gene]['copynumber'] = 'Heterozygous Loss'
                        cnv_loh[gene]['is_ichor_cn'] = True
                    elif ichor_cn >= 2.8:
                        cnv_loh[gene]['copynumber'] = f'Gain (CN={ichor_cn:.2f})'
                        cnv_loh[gene]['is_ichor_cn'] = True
                    else:
                        cnv_loh[gene]['copynumber'] = 'Normal'
                        cnv_loh[gene]['is_ichor_cn'] = True

        for gene in GENES_OF_INTEREST:
            is_ichor = cnv_loh[gene].get('is_ichor_cn', False)
            bias = bias_dict.get(gene)
            if bias and not is_ichor:
                corrected_val = bias['corrected_cn']
                if corrected_val < 0.5:
                    cnv_loh[gene]['copynumber'] = 'Homozygous Deletion'
                elif corrected_val < 1.7:
                    cnv_loh[gene]['copynumber'] = 'Heterozygous Loss'
                elif corrected_val > 2.8:
                    cnv_loh[gene]['copynumber'] = f'Gain (CN={corrected_val:.2f})'
                else:
                    cnv_loh[gene]['copynumber'] = 'Normal'
                    
        # Get msIndelsPerMb from purity file
        msi_indels_per_mb = "NA"
        purity_file = find_purple_file(purple_dir, "purity.tsv", group_id, sample_id=tumor_sample_id)
        if purity_file and os.path.exists(purity_file):
            try:
                with open(purity_file, mode='r', encoding='utf-8') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for r in reader:
                        val = r.get('msIndelsPerMb')
                        if val is not None:
                            msi_indels_per_mb = str(round(float(val), 5))
            except Exception as e:
                pass

        # Get MSI unstable sites proportion from purple.command.log
        unstable_prop = "NA"
        msi_indels_count, total_msi_sites, prop = get_msi_unstable_sites_proportion(purple_dir, group_id, sample_id=tumor_sample_id)
        if prop is not None:
            unstable_prop = str(round(prop, 5))

        row_data = {
            'Patient': subject_id or group_id,
            'Sample': sample_key,
            'Tumour.Fraction': purity,
            'Purity.Source': selected_source,
            'Purity.PURPLE': str(round(purity_dict['purple'], 5)) if purity_dict.get('purple') is not None else 'NA',
            'Purity.sWGS': str(round(purity_dict['swgs'], 5)) if purity_dict.get('swgs') is not None else 'NA',
            'Purity.Panel_ichor': str(round(purity_dict['panel_ichor'], 5)) if purity_dict.get('panel_ichor') is not None else 'NA',
            'Purity.VAF_Implied': str(round(purity_dict['vaf'], 5)) if purity_dict.get('vaf') is not None else 'NA',
            'MSI.Status': msi_status,
            'MSI.Indels.Per.Mb': msi_indels_per_mb,
            'MSI.Unstable.Sites.Proportion': unstable_prop,
            'Clinical.Classification': 'NA',
            'Clinical.Confidence': 'NA',
            'Clinical.Notes': 'NA',
            'Total.Reads': 'NA',
            'On.Target.Reads': 'NA',
            'On.Target.Coverage': 'NA',
            'Off.Target.Reads': 'NA',
            'Off.Target.Coverage': 'NA',
            'Panel.ichorCNA.MAD': 'NA',
            'sWGS.Coverage': 'NA',
            'sWGS.ichorCNA.MAD': 'NA'
        }
        
        for gene in GENES_OF_INTEREST:
            findings = []
            
            # 1. germline
            if gene in germline_muts:
                findings.append('germline')
            if gene in germline_deletions:
                findings.append('germline.deletion')
            # 2. somatic
            if gene in somatic_muts:
                findings.append('somatic')
            # 3. methylation
            if gene in methylation and methylation[gene].get('status') != 'unmethylated':
                findings.append('methylation')
            # 4. copynumber (cn.gain or cn.loss)
            cn_status = cnv_loh[gene]['copynumber']
            if cn_status:
                findings.append(cn_status)
            # 5. LOH
            if cnv_loh[gene]['LOH']:
                findings.append('LOH')
                
            row_data[gene] = '/'.join(findings) if findings else ''
            
        summary_rows.append(row_data)
        
    # 3. Create the summary directory and write the table
    summary_dir = os.path.join(args.output_dir, 'summary')
    os.makedirs(summary_dir, exist_ok=True)
    
    # Generate clinical MMR status report for each sample
    for sample in samples:
        group_id = sample['group_id']
        subject_id = sample['subject_id']
        tumor_sample_id, normal_sample_id = get_sample_ids_from_samplesheet(args.samplesheet, group_id)
        if not tumor_sample_id:
            tumor_sample_id = sample.get('sample_id')
        sample_key = tumor_sample_id or group_id

        purple_dir = os.path.join(args.output_dir, group_id, 'purple')
        if not os.path.exists(purple_dir):
            purple_dir = os.path.join(args.output_dir, 'purple')
        pave_dir = os.path.join(args.output_dir, group_id, 'pave')
        if not os.path.exists(pave_dir):
            pave_dir = os.path.join(args.output_dir, 'pave')
        orange_dir = os.path.join(args.output_dir, group_id, 'orange')
        if not os.path.exists(orange_dir):
            orange_dir = os.path.join(args.output_dir, 'orange')
            
        if os.path.exists(purple_dir):
            # Resolve sWGS and panel ichorCNA directories
            ichor_dir = resolve_swgs_ichor_dir(args.ichor_dir, args.output_dir, group_id, tumor_sample_id=tumor_sample_id)
            panel_ichor_dir = resolve_panel_ichor_dir(args.panel_ichor_dir, args.output_dir, group_id, tumor_sample_id=tumor_sample_id)
                    
            print(f"Generating clinical genomics MMR report for {group_id}...")
            clin_info = generate_clinical_report(purple_dir, pave_dir, group_id, subject_id, args.samplesheet, args.genome, args.methylation_data, summary_dir, gtf_path=args.gtf, ihc_tsv=args.ihc_tsv, ihc_file=args.ihc_file, cmd_age=args.age, cmd_tissue=args.tissue, ichor_dir=ichor_dir, min_somatic_vaf=args.min_somatic_vaf, panel_ichor_dir=panel_ichor_dir, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id, orange_dir=orange_dir, orthogonal_tsv=args.orthogonal_tsv, purity_source=args.purity_source)
            
            # Find the row in summary_rows and update it
            row = next((r for r in summary_rows if r['Sample'] in (sample_key, group_id)), None)

            if row and clin_info:
                row['Clinical.Classification'] = clin_info['classification']
                row['Clinical.Confidence'] = clin_info['confidence']
                row['Clinical.Notes'] = clin_info['notes']
                row['Mean.Driver.VAF.Diploid'] = clin_info.get('mean_driver_vaf_diploid', 'NA')
                row['Mean.Driver.VAF.All'] = clin_info.get('mean_driver_vaf_all', 'NA')
                row['mean driver VAF: diploid'] = clin_info.get('mean_driver_vaf_diploid', 'NA')
                row['mean driver VAF: all'] = clin_info.get('mean_driver_vaf_all', 'NA')
                qc_stats = clin_info.get('qc_stats')
                if qc_stats:
                    row['Total.Reads'] = str(qc_stats.get('total_reads', 'NA'))
                    row['On.Target.Reads'] = str(qc_stats.get('on_target_reads', 'NA'))
                    cov_on = qc_stats.get('on_target_coverage')
                    row['On.Target.Coverage'] = f"{float(cov_on):.2f}" if cov_on is not None and cov_on != 'NA' else 'NA'
                    row['Off.Target.Reads'] = str(qc_stats.get('off_target_reads', 'NA'))
                    cov_off = qc_stats.get('off_target_coverage')
                    row['Off.Target.Coverage'] = f"{float(cov_off):.2f}" if cov_off is not None and cov_off != 'NA' else 'NA'
                    p_mad = qc_stats.get('mad')
                    if p_mad is not None and p_mad != 'NA':
                        row['Panel.ichorCNA.MAD'] = f"{float(p_mad):.4f}"
                    else:
                        row['Panel.ichorCNA.MAD'] = 'NA'
                swgs_qc = clin_info.get('swgs_qc_stats')
                if swgs_qc:
                    c_val = swgs_qc.get('coverage') or swgs_qc.get('genome_coverage') or swgs_qc.get('off_target_coverage')
                    if c_val is not None and c_val != 'NA':
                        try:
                            if float(c_val) > 0.0:
                                row['sWGS.Coverage'] = f"{float(c_val):.2f}"
                        except (ValueError, TypeError):
                            pass
                    s_mad = swgs_qc.get('mad')
                    if s_mad is not None and s_mad != 'NA':
                        row['sWGS.ichorCNA.MAD'] = f"{float(s_mad):.4f}"
            
    output_file = os.path.join(summary_dir, f"{args.project}_final_summary.tsv")
    
    # Also write a copy directly to root of output-dir as final_summary.tsv for ease of access
    root_output_file = os.path.join(args.output_dir, "final_summary.tsv")
    
    purity_headers = ['Purity.Source', 'Purity.PURPLE', 'Purity.sWGS', 'Purity.Panel_ichor', 'Purity.VAF_Implied']
    qc_headers = ['Total.Reads', 'On.Target.Reads', 'On.Target.Coverage', 'Off.Target.Reads', 'Off.Target.Coverage', 'Panel.ichorCNA.MAD', 'sWGS.Coverage', 'sWGS.ichorCNA.MAD', 'Mean.Driver.VAF.Diploid', 'Mean.Driver.VAF.All', 'mean driver VAF: diploid', 'mean driver VAF: all']
    headers = ['Patient', 'Sample', 'Tumour.Fraction'] + purity_headers + ['MSI.Status', 'MSI.Indels.Per.Mb', 'MSI.Unstable.Sites.Proportion', 'Clinical.Classification', 'Clinical.Confidence', 'Clinical.Notes'] + qc_headers + GENES_OF_INTEREST
    
    for outfile in [output_file, root_output_file]:
        try:
            with open(outfile, mode='w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
                writer.writeheader()
                for row in summary_rows:
                    writer.writerow(row)
            print(f"Summary table successfully written to {outfile}")
        except Exception as e:
            print(f"Error: Failed to write output file {outfile}: {e}", file=sys.stderr)

    # 4. Generate BAF/LOH plots for genes with copy number, LOH, or mutation changes
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_script = os.path.join(script_dir, "plot_gene_baf.py")
    
    if os.path.exists(plot_script):
        print("Generating BAF/LOH plots for reported mutations/CNA/LOH genes...")
        for sample in samples:
            group_id = sample['group_id']
            subject_id = sample['subject_id']
            tumor_sample_id, normal_sample_id = get_sample_ids_from_samplesheet(args.samplesheet, group_id)
            if not tumor_sample_id:
                tumor_sample_id = sample.get('sample_id')
            sample_key = tumor_sample_id or group_id
            
            # Find input directories
            purple_dir = os.path.join(args.output_dir, group_id, 'purple')
            if not os.path.exists(purple_dir):
                purple_dir = os.path.join(args.output_dir, 'purple')
                
            pave_dir = os.path.join(args.output_dir, group_id, 'pave')
            if not os.path.exists(pave_dir):
                pave_dir = os.path.join(args.output_dir, 'pave')
                
            amber_dir = os.path.join(args.output_dir, group_id, 'amber')
            if not os.path.exists(amber_dir):
                amber_dir = os.path.join(args.output_dir, 'amber')
                
            # Locate input files
            amber_baf = None
            for p in filter(None, [tumor_sample_id, group_id]):
                cand = os.path.join(amber_dir, f"{p}.amber.baf.tsv.gz")
                if os.path.exists(cand):
                    amber_baf = cand
                    break
            if not amber_baf:
                import glob
                baf_matches = [f for f in glob.glob(os.path.join(amber_dir, "*.amber.baf.tsv.gz")) if not os.path.basename(f).startswith("._")]
                if baf_matches:
                    amber_baf = baf_matches[0]
                else:
                    amber_baf = os.path.join(amber_dir, f"{group_id}.amber.baf.tsv.gz")

            purple_segments = find_purple_file(purple_dir, "segment.tsv", group_id, sample_id=tumor_sample_id)
            if not purple_segments:
                purple_segments = os.path.join(purple_dir, f"{group_id}.purple.segment.tsv")
            
            if not os.path.exists(amber_baf):
                print(f"Warning: AMBER BAF file not found for sample {group_id} at {amber_baf}, skipping plots.", file=sys.stderr)
                continue
                
            row = next((r for r in summary_rows if r['Sample'] in (sample_key, group_id)), None)
            if not row:
                continue
                
            # Extract purity for expected raw BAF calculations
            purity_val = 1.0
            purity_str = row.get('Tumour.Fraction', '1.0')
            if purity_str != 'NA':
                try:
                    purity_val = float(purity_str)
                except ValueError:
                    pass
                    
            # Parse methylation data for CpG-level plotting
            sample_meth_data = find_methylation_data(args.methylation_data, args.output_dir, group_id)
            methylation = parse_methylation(sample_meth_data, group_id, subject_id, args.samplesheet, args.genome)
            ichor_dir = resolve_swgs_ichor_dir(args.ichor_dir, args.output_dir, group_id, tumor_sample_id=tumor_sample_id)
            
            for gene in GENES_OF_INTEREST:
                status = row.get(gene, '')
                # Generate plot if there is somatic/germline mutation, LOH, copy number gain/loss, or methylation
                has_cn_loh_mut = any(x in status for x in ['somatic', 'germline', 'LOH', 'cn.loss', 'cn.gain', 'Loss', 'Deletion', 'Gain', 'methylation', 'germline.deletion'])
                
                if has_cn_loh_mut:
                    variants = get_detailed_variants(purple_dir, pave_dir, group_id, gene, subject_id, args.samplesheet, args.genome, sample_meth_data, tumor_sample_id=tumor_sample_id, normal_sample_id=normal_sample_id)
                    variants_json = json.dumps(variants)
                    
                    # Get CpG methylation details
                    meth_data = methylation.get(gene)
                    meth_json = json.dumps(meth_data) if meth_data else None
                    
                    # Output plot paths (both clean PNG and interactive HTML)
                    plot_out = os.path.join(summary_dir, 'plots', f"{group_id}_{gene}_baf.png")
                    plot_html = os.path.join(summary_dir, 'plots', f"{group_id}_{gene}_baf.html")
                    
                    cmd = [
                        sys.executable, plot_script,
                        "--amber-baf", amber_baf,
                        "--gene", gene,
                        "--output", plot_out,
                        "--output-html", plot_html,
                        "--genome", args.genome,
                        "--purity", str(purity_val)
                    ]
                    if purple_segments and os.path.exists(purple_segments):
                        cmd.extend(["--purple-segments", purple_segments])
                    
                    # Target region capture probe resolution
                    target_cn_file = os.path.join(purple_dir, f"{group_id}.purple.target_region_cn.tsv")
                    if not os.path.exists(target_cn_file) and tumor_sample_id:
                        target_cn_file = os.path.join(purple_dir, f"{tumor_sample_id}.purple.target_region_cn.tsv")
                    if os.path.exists(target_cn_file):
                        cmd.extend(["--target-region-cn", target_cn_file])
                        
                    if variants:
                        cmd.extend(["--variants-json", variants_json])
                    if meth_json:
                        cmd.extend(["--methylation-json", meth_json])
                    if ichor_dir:
                        cmd.extend(["--ichor-dir", ichor_dir, "--sample-id", group_id, "--patient-id", subject_id])
                        
                    try:
                        print(f"  Plotting BAF for {group_id} - {gene}...")
                        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    except subprocess.CalledProcessError as e:
                        print(f"Warning: BAF plot failed for {group_id} - {gene}: {e.stderr.decode()}", file=sys.stderr)
                    except Exception as e:
                        print(f"Warning: BAF plot failed for {group_id} - {gene}: {e}", file=sys.stderr)
                        
            # 5. Generate IGV BED and HTML report for the sample
            generate_igv_bed_and_report(
                group_id=group_id,
                subject_id=subject_id,
                purple_dir=purple_dir,
                pave_dir=pave_dir,
                samplesheet_path=args.samplesheet,
                genome=args.genome,
                methylation_data=args.methylation_data,
                summary_dir=summary_dir,
                fasta_path=args.fasta,
                gtf_path=args.gtf,
                igv_env=args.igv_reports_env,
                tumor_sample_id=tumor_sample_id,
                normal_sample_id=normal_sample_id
            )
    else:
        print(f"Warning: Plotting script not found at {plot_script}, skipping plots.", file=sys.stderr)

if __name__ == '__main__':
    main()
