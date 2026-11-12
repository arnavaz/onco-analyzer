import os
import csv
import glob
import subprocess
import re
import argparse
import shlex
from multiprocessing import Pool

# Global paths
source_repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_base = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results"
uhn_ichor_dir = "/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/sWGS/IchorCNA/"
otb_ichor_dir = "/cluster/projects/pughlab/pughlab/projects/MultiMMR/OTB/sWGS/IchorCNA/"
methylation_data = "/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/EMSeq/MethylDackel/2025-10-10_MultiMMR_CpG_methylation_matrix.RData"
ihc_tsv = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v3/DNASeq/cohort_ihc.tsv"
panel_name = "multimmr_hmf"

endometrial_list = ["163243", "39096", "76452", "41463", "714363", "714364", "714576", "714823", "714834", "714960", "715035", "715042", "715465", "715620", "715628", "715779", "715900", "715916", "716627", "716666", "716757"]

def run_cmd_with_modules(cmd_list, stdout=None):
    cmd_str = " ".join(shlex.quote(arg) for arg in cmd_list)
    full_cmd = (
        "for init in \"$MODULESHOME/init/bash\" \"/usr/share/lmod/lmod/init/bash\" \"/usr/share/Modules/init/bash\" \"/etc/profile.d/modules.sh\" \"/etc/profile.d/z00_lmod.sh\" \"/etc/profile.d/z01_lmod.sh\"; do "
        "  if [ -f \"$init\" ]; then source \"$init\"; break; fi; "
        "done; "
        "if ! command -v readCounter &>/dev/null; then module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true; fi; "
        "if ! command -v Rscript &>/dev/null; then module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R 2>/dev/null || true; fi; "
        "if ! command -v samtools &>/dev/null; then module load samtools/1.19 2>/dev/null || module load samtools 2>/dev/null || true; fi; "
        f"{cmd_str}"
    )
    wrapper_cmd = ["bash", "-c", full_cmd]
    if stdout:
        return subprocess.run(wrapper_cmd, stdout=stdout, check=True)
    else:
        return subprocess.run(wrapper_cmd, check=True)

def check_output_with_modules(cmd_list):
    cmd_str = " ".join(shlex.quote(arg) for arg in cmd_list)
    full_cmd = (
        "for init in \"$MODULESHOME/init/bash\" \"/usr/share/lmod/lmod/init/bash\" \"/usr/share/Modules/init/bash\" \"/etc/profile.d/modules.sh\" \"/etc/profile.d/z00_lmod.sh\" \"/etc/profile.d/z01_lmod.sh\"; do "
        "  if [ -f \"$init\" ]; then source \"$init\"; break; fi; "
        "done; "
        "if ! command -v readCounter &>/dev/null; then module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true; fi; "
        "if ! command -v Rscript &>/dev/null; then module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R 2>/dev/null || true; fi; "
        "if ! command -v samtools &>/dev/null; then module load samtools/1.19 2>/dev/null || module load samtools 2>/dev/null || true; fi; "
        f"{cmd_str}"
    )
    wrapper_cmd = ["bash", "-c", full_cmd]
    return subprocess.check_output(wrapper_cmd)

def run_off_target_ichorcna(sample_dir, sheet, enable_subclones=False):
    """Re-run panel off-target ichorCNA using matched normal WIG and tuned parameters."""
    # 1. Find all redux.bam files in sample_dir
    bam_files = glob.glob(os.path.join(sample_dir, "**", "*.redux.bam"), recursive=True)
    bam_map = {}
    for bf in bam_files:
        sample_id = os.path.basename(bf).replace(".redux.bam", "")
        bam_map[sample_id] = bf
        
    # 2. Parse samplesheet to find tumor-normal pairing
    tumor_id = None
    tumor_path = None
    normal_id = None
    normal_path = None
    
    try:
        with open(sheet, 'r') as f:
            lines = [line for line in f if not line.strip().startswith('#')]
            reader = csv.DictReader(lines)
            samples = list(reader)
    except Exception as e:
        print(f"Error reading samplesheet {sheet}: {e}")
        return
        
    for row in samples:
        if row.get('sample_type', '').lower() == 'tumor' and row.get('sequence_type', '').lower() == 'dna':
            t_id = row['sample_id']
            # Match bam
            matched_key = next((k for k in bam_map.keys() if t_id in k or k in t_id), None)
            if matched_key:
                tumor_path = bam_map[matched_key]
                tumor_id = matched_key
                
            subject_id = row.get('subject_id')
            group_id = row.get('group_id')
            for nrow in samples:
                if nrow.get('sample_type', '').lower() == 'normal' and nrow.get('sequence_type', '').lower() == 'dna':
                    if nrow.get('subject_id') == subject_id or (group_id and nrow.get('group_id') == group_id):
                        n_id = nrow['sample_id']
                        matched_norm_key = next((k for k in bam_map.keys() if n_id in k or k in n_id), None)
                        if matched_norm_key:
                            normal_path = bam_map[matched_norm_key]
                            normal_id = matched_norm_key
                            break
            break
            
    if not tumor_path:
        print(f"Warning: No tumor redux BAM found for {sample_dir}, skipping ichorCNA.")
        return
        
    sample_ichor_dir = os.path.join(sample_dir, "ichorCNA", tumor_id)
    os.makedirs(sample_ichor_dir, exist_ok=True)
    
    # 3. Setup paths and parameters
    exclusion_bed = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_exclusion.bed"
    ichor_pon = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_panelOfNormals.rds"
    
    chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    genome_build = "hg38"
    chr_norm = "paste0('chr', 1:22)"
    chr_train = "paste0('chr', 1:22)"
    chr_analyze = "c(paste0('chr', 1:22),'chrX')"
    genome_style = "UCSC"
    
    # 4. Generate tumor WIG
    tumor_raw_wig = os.path.join(sample_ichor_dir, f"{tumor_id}.raw.wig")
    try:
        # idxstats to get total reads
        idx_out = check_output_with_modules(["samtools", "idxstats", tumor_path]).decode()
        total_reads = sum(int(line.split('\t')[2]) for line in idx_out.strip().split('\n') if line)
        
        # on-target reads
        on_target = int(check_output_with_modules(["samtools", "view", "-c", "-L", exclusion_bed, tumor_path]).decode().strip())
        off_target_reads = max(0, total_reads - on_target)
        off_target_cov = (off_target_reads * 150) / 3100000000
        
        downsample_bam = tumor_path
        is_temp_bam = False
        if off_target_cov > 3.0:
            fraction = 3.0 / off_target_cov
            s_param = f"{100 + fraction:.6f}"
            temp_bam = os.path.join(sample_ichor_dir, f"{tumor_id}.downsampled.bam")
            run_cmd_with_modules(["samtools", "view", "-b", "-s", s_param, tumor_path], stdout=open(temp_bam, 'wb'))
            run_cmd_with_modules(["samtools", "index", temp_bam])
            downsample_bam = temp_bam
            is_temp_bam = True
            
        run_cmd_with_modules(["readCounter", "--window", "1000000", "--quality", "20", "--chromosome", chrs, downsample_bam], stdout=open(tumor_raw_wig, 'wb'))
        if is_temp_bam:
            os.remove(downsample_bam)
            if os.path.exists(downsample_bam + ".bai"):
                os.remove(downsample_bam + ".bai")
    except Exception as e:
        print(f"Error counting tumor reads for {tumor_id}: {e}")
        return
        
    tumor_wig = os.path.join(sample_ichor_dir, f"{tumor_id}.wig")
    run_cmd_with_modules(["python3", os.path.join(source_repo, "scripts/filter_wig.py"), tumor_raw_wig, exclusion_bed, tumor_wig])
    if os.path.exists(tumor_raw_wig):
        os.remove(tumor_raw_wig)
        
    # Calculate target regions length from exclusion bed
    target_length = 0
    try:
        import gzip
        if exclusion_bed.endswith('.gz'):
            with gzip.open(exclusion_bed, 'rt') as f_bed:
                for line in f_bed:
                    if not line.strip() or line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        target_length += int(parts[2]) - int(parts[1])
        else:
            with open(exclusion_bed, 'r') as f_bed:
                for line in f_bed:
                    if not line.strip() or line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        target_length += int(parts[2]) - int(parts[1])
    except Exception as e:
        print(f"Warning: Failed to parse exclusion bed size: {e}")
        target_length = 0

    if target_length > 0:
        on_target_cov = (on_target * 150) / target_length
        off_target_cov = (off_target_reads * 150) / (3088286401 - target_length)
    else:
        on_target_cov = 0.0
        off_target_cov = (off_target_reads * 150) / 3088286401

    import json
    qc_stats = {
        "total_reads": total_reads,
        "on_target_reads": on_target,
        "on_target_coverage": on_target_cov,
        "off_target_reads": off_target_reads,
        "off_target_coverage": off_target_cov
    }
    with open(os.path.join(sample_ichor_dir, "qc_stats.json"), "w") as f_json:
        json.dump(qc_stats, f_json, indent=2)
    sub_qc_dir = os.path.join(sample_ichor_dir, tumor_id)
    os.makedirs(sub_qc_dir, exist_ok=True)
    with open(os.path.join(sub_qc_dir, "qc_stats.json"), "w") as f_json:
        json.dump(qc_stats, f_json, indent=2)

    # 5. Generate normal WIG if available
    normal_wig = None
    if normal_path and normal_id:
        normal_raw_wig = os.path.join(sample_ichor_dir, f"{normal_id}.raw.wig")
        try:
            n_idx_out = check_output_with_modules(["samtools", "idxstats", normal_path]).decode()
            n_total_reads = sum(int(line.split('\t')[2]) for line in n_idx_out.strip().split('\n') if line)
            
            n_on_target = int(check_output_with_modules(["samtools", "view", "-c", "-L", exclusion_bed, normal_path]).decode().strip())
            n_off_target_reads = max(0, n_total_reads - n_on_target)
            n_off_target_cov = (n_off_target_reads * 150) / 3100000000
            
            n_downsample_bam = normal_path
            n_is_temp_bam = False
            if n_off_target_cov > 3.0:
                n_fraction = 3.0 / n_off_target_cov
                n_s_param = f"{100 + n_fraction:.6f}"
                n_temp_bam = os.path.join(sample_ichor_dir, f"{normal_id}.downsampled.bam")
                run_cmd_with_modules(["samtools", "view", "-b", "-s", n_s_param, normal_path], stdout=open(n_temp_bam, 'wb'))
                run_cmd_with_modules(["samtools", "index", n_temp_bam])
                n_downsample_bam = n_temp_bam
                n_is_temp_bam = True
                
            run_cmd_with_modules(["readCounter", "--window", "1000000", "--quality", "20", "--chromosome", chrs, n_downsample_bam], stdout=open(normal_raw_wig, 'wb'))
            if n_is_temp_bam:
                os.remove(n_downsample_bam)
                if os.path.exists(n_downsample_bam + ".bai"):
                    os.remove(n_downsample_bam + ".bai")
        except Exception as e:
            print(f"Error counting normal reads for {normal_id}: {e}")
            return
            
        normal_wig = os.path.join(sample_ichor_dir, f"{normal_id}.wig")
        run_cmd_with_modules(["python3", os.path.join(source_repo, "scripts/filter_wig.py"), normal_raw_wig, exclusion_bed, normal_wig])
        if os.path.exists(normal_raw_wig):
            os.remove(normal_raw_wig)
            
    # 6. Resolve R references dynamically
    try:
        gc_wig = check_output_with_modules(["Rscript", "-e", "cat(system.file('extdata', 'gc_hg38_1000kb.wig', package='ichorCNA'))"]).decode().strip()
        map_wig = check_output_with_modules(["Rscript", "-e", "cat(system.file('extdata', 'map_hg38_1000kb.wig', package='ichorCNA'))"]).decode().strip()
        centromere = check_output_with_modules(["Rscript", "-e", "cat(system.file('extdata', 'GRCh38.GCA_000001405.2_centromere_acen.txt', package='ichorCNA'))"]).decode().strip()
    except Exception as e:
        print(f"Error resolving R references for {tumor_id}: {e}")
        return
        
    sc_states = "c(1,3)" if enable_subclones else "NULL"
    
    cmd_r = [
        "Rscript", os.path.join(source_repo, "resources/runIchorCNA.R"),
        "--WIG", tumor_wig,
        "--id", tumor_id,
        "--genomeStyle", genome_style,
        "--genomeBuild", genome_build,
        "--gcWig", gc_wig,
        "--mapWig", map_wig,
        "--centromere", centromere,
        "--chrNormalize", chr_norm,
        "--chrTrain", chr_train,
        "--chrs", chr_analyze,
        "--normal", "c(0.5,0.6,0.7,0.8,0.9)",
        "--ploidy", "c(2,3)",
        "--scStates", sc_states,
        "--outDir", sample_ichor_dir
    ]
    if normal_wig:
        cmd_r.extend(["--NORMWIG", normal_wig])
    cmd_r.extend(["--normalPanel", ichor_pon])
        
    print(f"  Executing runIchorCNA.R for {tumor_id}...")
    try:
        run_cmd_with_modules(cmd_r)
    except Exception as e:
        print(f"Error executing runIchorCNA.R: {e}")
    
    # Cleanup wigs
    if os.path.exists(tumor_wig):
        os.remove(tumor_wig)
    if normal_wig and os.path.exists(normal_wig):
        os.remove(normal_wig)

def process_sheet_wrapper(args_tuple):
    sheet, run_ichor, subclones = args_tuple
    sample_dir = os.path.dirname(sheet)
    sample_name = os.path.basename(sample_dir)
    
    if run_ichor:
        print(f"Re-running off-target ichorCNA for {sample_name}...")
        try:
            run_off_target_ichorcna(sample_dir, sheet, enable_subclones=subclones)
        except Exception as e:
            print(f"Error re-running ichorCNA for {sample_name}: {e}")
            
    clean_name = sample_name.replace("DYNA-", "").replace("DYNA_", "")
    m = re.search(r"[A-Z]+-[0-9]+", clean_name)
    patient_id = m.group(0) if m else clean_name
    
    ichor_dir = None
    for parent_dir in [uhn_ichor_dir, otb_ichor_dir]:
        if os.path.exists(parent_dir):
            candidates = [
                os.path.join(parent_dir, patient_id, f"{patient_id}-T"),
                os.path.join(parent_dir, patient_id),
                os.path.join(parent_dir, clean_name, f"{clean_name}-T"),
                os.path.join(parent_dir, sample_name)
            ]
            for c in candidates:
                if os.path.isdir(c):
                    ichor_dir = c
                    break
            if ichor_dir:
                break
                
            # Fuzzy match 2 levels deep
            if not ichor_dir:
                for root, dirs, files in os.walk(parent_dir, followlinks=True):
                    depth = root[len(parent_dir):].count(os.sep)
                    if depth > 1:
                        continue
                    for d in dirs:
                        if patient_id in d and d.endswith("-T"):
                            ichor_dir = os.path.join(root, d)
                            break
                    if ichor_dir:
                        break
    
    patient_code = sample_name.split("-")[-1]
    tissue = "endometrial" if patient_code in endometrial_list else "colorectal"
    
    cmd = [
        "python3", os.path.join(source_repo, "scripts/summarize_findings.py"),
        "--output-dir", sample_dir,
        "--samplesheet", sheet,
        "--methylation-data", methylation_data,
        "--ihc-tsv", ihc_tsv,
        "--tissue", tissue,
        "--project", panel_name,
        "--genome", "GRCh38_hmf",
        "--gtf", "/cluster/projects/pughlab/references/RefSeq/SYMBOL_MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"
    ]
    if ichor_dir:
        cmd += ["--ichor-dir", ichor_dir]
        
    print(f"Running summarize_findings for {sample_name}...")
    subprocess.run(cmd)

def main():
    parser = argparse.ArgumentParser(description='Regenerate reports in parallel.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of parallel threads/processes (default: 1)')
    parser.add_argument('--run-ichor', action='store_true', help='Re-run panel off-target ichorCNA on the redux BAMs')
    parser.add_argument('--subclones', action='store_true', help='Enable subclonal fraction estimation in ichorCNA')
    args = parser.parse_args()

    sample_sheets = glob.glob(os.path.join(results_base, "*", "samplesheet.csv"))
    print(f"Found {len(sample_sheets)} samplesheets to process. Using {args.threads} threads.")
    
    tasks = [(sheet, args.run_ichor, args.subclones) for sheet in sample_sheets]
    
    if args.threads > 1:
        with Pool(args.threads) as pool:
            pool.map(process_sheet_wrapper, tasks)
    else:
        for task in tasks:
            process_sheet_wrapper(task)

if __name__ == '__main__':
    main()
