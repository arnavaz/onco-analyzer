#!/usr/bin/env python3
"""
benchmark_ichor_bin_sizes.py
-----------------------------
Benchmarking script to evaluate ichorCNA copy number noise (MAD score),
tumor fraction estimates, and segmentation across alternate bin sizes and normalization modes:
  - Paired (Matched Normal WIG) + Custom PoN (1Mb)
  - Paired (Matched Normal WIG) Alone (1Mb)
  - Tumor-Only + Custom PoN (500kb, 1Mb, 2Mb)
  - Tumor-Only + Default/No PoN (500kb, 1Mb, 2Mb)

Target Sample: DYNA-OTB-715257-FFPE-T-DNA-U_S214
Matched Normal: DYNA-OTB-715257-B-N-DNA-U_S269
Matches exact runIchorCNA.R invocation used in run_oncoanalyser.sh.
"""

import os
import sys
import glob
import shlex
import subprocess

# Default Paths on the Cluster
SAMPLE_ID = "DYNA-OTB-715257-FFPE-T-DNA-U_S214"
NORMAL_SAMPLE_ID = "DYNA-OTB-715257-B-N-DNA-U_S269"

BASE_DIR = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/DYNA-OTB-715257"
BAM_PATH = os.path.join(BASE_DIR, SAMPLE_ID, "alignments/dna", f"{SAMPLE_ID}.redux.bam")
CRAM_PATH = os.path.join(BASE_DIR, SAMPLE_ID, "alignments/dna", f"{SAMPLE_ID}.redux.cram")

NORMAL_BAM_PATH = os.path.join(BASE_DIR, NORMAL_SAMPLE_ID, "alignments/dna", f"{NORMAL_SAMPLE_ID}.redux.bam")
NORMAL_CRAM_PATH = os.path.join(BASE_DIR, NORMAL_SAMPLE_ID, "alignments/dna", f"{NORMAL_SAMPLE_ID}.redux.cram")

CUSTOM_PON = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_panelOfNormals.rds"
OUT_DIR = os.path.join(BASE_DIR, "bin_size_benchmark")

# Chromosomes to include for readCounter
CHRS = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def run_bash_cmd(cmd_str, silent=False):
    """Executes bash command with module loading support for cluster environments."""
    full_cmd = (
        'for init in "$MODULESHOME/init/bash" "/usr/share/lmod/lmod/init/bash" "/usr/share/Modules/init/bash" "/etc/profile.d/modules.sh" "/etc/profile.d/z00_lmod.sh" "/etc/profile.d/z01_lmod.sh"; do '
        '  if [ -f "$init" ]; then source "$init"; break; fi; '
        'done; '
        'module load hmmcopy_utils/170718 2>/dev/null || module load hmmcopy_utils 2>/dev/null || true; '
        'module load R/4.1.0 2>/dev/null || module load R/4.2.0 2>/dev/null || module load R 2>/dev/null || true; '
        'module load samtools/1.19 2>/dev/null || module load samtools/1.20 2>/dev/null || module load samtools 2>/dev/null || true; '
        f'{cmd_str}'
    )
    if not silent:
        print(f"\n[RUNNING]: {cmd_str}")
    res = subprocess.run(["bash", "-c", full_cmd], capture_output=True, text=True)
    if not silent and res.returncode != 0:
        print(f"[WARNING/ERROR]: Exit code {res.returncode}")
        if res.stdout:
            print(f"STDOUT:\n{res.stdout}")
        if res.stderr:
            print(f"STDERR:\n{res.stderr}")
    return res

def parse_params_summary(params_file):
    """Parses key stats from ichorCNA params.txt output."""
    stats = {'mad': 'N/A', 'purity': 'N/A', 'ploidy': 'N/A'}
    if not os.path.exists(params_file):
        return stats
    with open(params_file, 'r') as f:
        for line in f:
            line_s = line.strip()
            if line_s.startswith("GC-Map correction MAD:"):
                stats['mad'] = line_s.split(":")[-1].strip()
            elif line_s.startswith("Tumor Fraction:"):
                stats['purity'] = line_s.split(":")[-1].strip()
            elif line_s.startswith("Ploidy:"):
                stats['ploidy'] = line_s.split(":")[-1].strip()
    return stats

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.normpath(os.path.join(script_dir, ".."))
    run_ichor_script = os.path.join(repo_root, "resources/runIchorCNA.R")
    
    if not os.path.exists(run_ichor_script):
        print(f"ERROR: Could not locate runIchorCNA.R at {run_ichor_script}", file=sys.stderr)
        sys.exit(1)
        
    # Recursively find clean tumor and normal redux BAM/CRAM files in BASE_DIR (preferring .bam over .cram)
    all_bams_crams = [f for f in glob.glob(os.path.join(BASE_DIR, "**", "*.redux.bam"), recursive=True) +
                              glob.glob(os.path.join(BASE_DIR, "**", "*.redux.cram"), recursive=True)
                   if not os.path.basename(f).startswith("00") and "esvee" not in f]

    target_align_file = None
    target_normal_align_file = None

    # Prefer BAM files first
    for bf in sorted(all_bams_crams, key=lambda x: (x.endswith(".cram"), x)):
        fname = os.path.basename(bf)
        if SAMPLE_ID in fname and target_align_file is None:
            target_align_file = bf
        elif (NORMAL_SAMPLE_ID in fname or "-N-" in fname or "-B-N-" in fname) and target_normal_align_file is None:
            target_normal_align_file = bf

    if not target_align_file or not os.path.exists(target_align_file):
        target_align_file = BAM_PATH if os.path.exists(BAM_PATH) else CRAM_PATH

    if not target_normal_align_file or not os.path.exists(target_normal_align_file):
        target_normal_align_file = NORMAL_BAM_PATH if os.path.exists(NORMAL_BAM_PATH) else NORMAL_CRAM_PATH

    has_normal = target_normal_align_file is not None and os.path.exists(target_normal_align_file)

    if not target_align_file or not os.path.exists(target_align_file):
        print(f"ERROR: Alignment file not found for sample {SAMPLE_ID} at {BAM_PATH} or {CRAM_PATH}", file=sys.stderr)
        sys.exit(1)

    print(f"=== ichorCNA Bin Size Benchmarking ===")
    print(f"Tumor Sample ID:    {SAMPLE_ID}")
    print(f"Alignment File:     {target_align_file}")
    print(f"Matched Normal ID:  {NORMAL_SAMPLE_ID}")
    print(f"Matched Normal File:{target_normal_align_file} (Exists: {has_normal})")
    print(f"runIchorCNA.R:      {run_ichor_script}")
    print(f"Custom PoN:         {CUSTOM_PON} (Exists: {os.path.exists(CUSTOM_PON)})")
    print(f"Benchmark Dir:      {OUT_DIR}\n")
    
    bin_configs = [
        {"name": "1Mb_paired_custom_pon", "window": 1000000, "wig_name": "1000kb", "bin_factor": 1, "pon": CUSTOM_PON, "use_normal": True},
        {"name": "1Mb_paired_no_pon", "window": 1000000, "wig_name": "1000kb", "bin_factor": 1, "pon": None, "use_normal": True},
        {"name": "1Mb_custom_pon", "window": 1000000, "wig_name": "1000kb", "bin_factor": 1, "pon": CUSTOM_PON, "use_normal": False},
        {"name": "1Mb_default_pon", "window": 1000000, "wig_name": "1000kb", "bin_factor": 1, "pon": None, "use_normal": False},
        {"name": "500kb_custom_pon", "window": 500000, "wig_name": "500kb", "bin_factor": 1, "pon": CUSTOM_PON, "use_normal": False},
        {"name": "500kb_default_pon", "window": 500000, "wig_name": "500kb", "bin_factor": 1, "pon": None, "use_normal": False},
        {"name": "2Mb_custom_pon", "window": 2000000, "wig_name": "2000kb", "bin_factor": 2, "pon": CUSTOM_PON, "use_normal": False},
        {"name": "2Mb_default_pon", "window": 2000000, "wig_name": "2000kb", "bin_factor": 2, "pon": None, "use_normal": False},
    ]
    
    results = []
    
    for cfg in bin_configs:
        cfg_name = cfg["name"]
        window_size = cfg["window"]
        wig_label = cfg["wig_name"]
        bin_factor = cfg["bin_factor"]
        use_pon = cfg["pon"]
        use_normal = cfg.get("use_normal", False)
        
        if use_normal and not has_normal:
            print(f"Skipping {cfg_name} because matched normal alignment file was not found.")
            continue

        cfg_out_dir = os.path.join(OUT_DIR, cfg_name)
        os.makedirs(cfg_out_dir, exist_ok=True)
        
        wig_file = os.path.join(OUT_DIR, f"{SAMPLE_ID}_{window_size}.wig")
        normal_wig_file = os.path.join(OUT_DIR, f"{NORMAL_SAMPLE_ID}_{window_size}.wig")
        
        # Step A: Generate Tumor WIG if not already present
        if os.path.exists(wig_file) and os.path.getsize(wig_file) == 0:
            os.remove(wig_file)

        if not os.path.exists(wig_file):
            print(f"Generating WIG file for bin size {window_size} bp...")
            align_in = target_align_file
            is_temp_tumor_bam = False
            if target_align_file.endswith(".cram"):
                temp_bam = os.path.join(OUT_DIR, f"{SAMPLE_ID}.temp.bam")
                print(f"Converting tumor CRAM to BAM for readCounter: {temp_bam}...")
                run_bash_cmd(f"samtools view -b -o {shlex.quote(temp_bam)} {shlex.quote(target_align_file)} && samtools index {shlex.quote(temp_bam)}")
                align_in = temp_bam
                is_temp_tumor_bam = True

            rc_cmd = f"readCounter --window {window_size} --quality 20 --chromosome {shlex.quote(CHRS)} {shlex.quote(align_in)} > {shlex.quote(wig_file)}"
            run_bash_cmd(rc_cmd)
            
            if is_temp_tumor_bam and os.path.exists(align_in):
                os.remove(align_in)
                if os.path.exists(align_in + ".bai"):
                    os.remove(align_in + ".bai")
        else:
            print(f"Reusing existing tumor WIG file: {wig_file}")
            
        # Step B: Generate Normal WIG if required & not already present
        if use_normal:
            if os.path.exists(normal_wig_file) and os.path.getsize(normal_wig_file) == 0:
                os.remove(normal_wig_file)

            if not os.path.exists(normal_wig_file):
                print(f"Generating Matched Normal WIG file for bin size {window_size} bp...")
                norm_align_in = target_normal_align_file
                is_temp_norm_bam = False
                if target_normal_align_file.endswith(".cram"):
                    temp_norm_bam = os.path.join(OUT_DIR, f"{NORMAL_SAMPLE_ID}.temp.bam")
                    print(f"Converting normal CRAM to BAM for readCounter: {temp_norm_bam}...")
                    run_bash_cmd(f"samtools view -b -o {shlex.quote(temp_norm_bam)} {shlex.quote(target_normal_align_file)} && samtools index {shlex.quote(temp_norm_bam)}")
                    norm_align_in = temp_norm_bam
                    is_temp_norm_bam = True

                rc_norm_cmd = f"readCounter --window {window_size} --quality 20 --chromosome {shlex.quote(CHRS)} {shlex.quote(norm_align_in)} > {shlex.quote(normal_wig_file)}"
                run_bash_cmd(rc_norm_cmd)

                if is_temp_norm_bam and os.path.exists(norm_align_in):
                    os.remove(norm_align_in)
                    if os.path.exists(norm_align_in + ".bai"):
                        os.remove(norm_align_in + ".bai")
            else:
                print(f"Reusing existing normal WIG file: {normal_wig_file}")

        # Step C: Build bash pre-command for GC, Map, and Centromere resolution on cluster
        if bin_factor == 1:
            ref_setup_cmd = (
                f"GC_WIG=$(Rscript -e \"cat(system.file('extdata', 'gc_hg38_{wig_label}.wig', package='ichorCNA'))\"); "
                f"MAP_WIG=$(Rscript -e \"cat(system.file('extdata', 'map_hg38_{wig_label}.wig', package='ichorCNA'))\"); "
            )
        else:
            gc_out = os.path.join(OUT_DIR, f"gc_hg38_{wig_label}.wig")
            map_out = os.path.join(OUT_DIR, f"map_hg38_{wig_label}.wig")

            # Remove previous corrupted files containing scientific notation if any
            for f_check in [gc_out, map_out]:
                if os.path.exists(f_check):
                    try:
                        with open(f_check, 'r') as f_in:
                            header_sample = f_in.read(500)
                            if 'e+' in header_sample or 'e-' in header_sample:
                                os.remove(f_check)
                    except Exception:
                        pass

            r_agg_code = f"""
            aggregateWig <- function(inFile, outFile, factor) {{
              if (file.exists(outFile) && file.info(outFile)$size > 0) return(outFile)
              lines <- readLines(inFile)
              header_idx <- grep('^fixedStep|^variableStep', lines)
              out_lines <- c()
              num_headers <- length(header_idx)
              new_step_str <- sprintf('step=%d', as.integer(1000000 * factor))
              new_span_str <- sprintf('span=%d', as.integer(1000000 * factor))
              for (h in 1:num_headers) {{
                start_i <- header_idx[h] + 1
                end_i <- if (h < num_headers) header_idx[h+1] - 1 else length(lines)
                if (start_i > end_i) next
                hdr <- lines[header_idx[h]]
                hdr <- gsub('step=1000000', new_step_str, hdr)
                hdr <- gsub('span=1000000', new_span_str, hdr)
                out_lines <- c(out_lines, hdr)
                vals <- as.numeric(lines[start_i:end_i])
                vals <- vals[!is.na(vals)]
                if (length(vals) == 0) next
                n <- length(vals)
                n_chunks <- ceiling(n / factor)
                for (c in 1:n_chunks) {{
                  chunk <- vals[((c-1)*factor + 1):min(c*factor, n)]
                  out_lines <- c(out_lines, sprintf('%.6f', mean(chunk)))
                }}
              }}
              writeLines(out_lines, outFile)
            }}
            gc1 <- system.file('extdata', 'gc_hg38_1000kb.wig', package='ichorCNA')
            map1 <- system.file('extdata', 'map_hg38_1000kb.wig', package='ichorCNA')
            aggregateWig(gc1, '{gc_out}', {bin_factor})
            aggregateWig(map1, '{map_out}', {bin_factor})
            """
            r_script_path = os.path.join(OUT_DIR, f"agg_{wig_label}.R")
            with open(r_script_path, 'w') as f:
                f.write(r_agg_code)
                
            ref_setup_cmd = (
                f"Rscript {shlex.quote(r_script_path)}; "
                f"GC_WIG={shlex.quote(gc_out)}; "
                f"MAP_WIG={shlex.quote(map_out)}; "
            )

        ref_setup_cmd += "CENTROMERE=$(Rscript -e \"cat(system.file('extdata', 'GRCh38.GCA_000001405.2_centromere_acen.txt', package='ichorCNA'))\"); "

        pon_arg = f"--normalPanel {shlex.quote(use_pon)}" if use_pon and os.path.exists(use_pon) else ""
        norm_wig_arg = f"--NORMWIG {shlex.quote(normal_wig_file)}" if use_normal else ""
        
        # Match exact runIchorCNA.R invocation from run_oncoanalyser.sh
        ichor_cmd = (
            f"cd {shlex.quote(cfg_out_dir)} && "
            f"{ref_setup_cmd} "
            f"Rscript {shlex.quote(run_ichor_script)} "
            f"--WIG {shlex.quote(wig_file)} "
            f"{norm_wig_arg} "
            f"--id {shlex.quote(SAMPLE_ID)} "
            f"--genomeStyle UCSC "
            f"--genomeBuild hg38 "
            f"--gcWig \"$GC_WIG\" "
            f"--mapWig \"$MAP_WIG\" "
            f"--centromere \"$CENTROMERE\" "
            f"--chrNormalize \"c(1:22)\" "
            f"--chrTrain \"c(1:22)\" "
            f"--chrs \"c(1:22, 'X')\" "
            f"--normal \"c(0.5,0.6,0.7,0.8,0.9)\" "
            f"--ploidy \"c(2,3)\" "
            f"--scStates \"NULL\" "
            f"{pon_arg}"
        )
        
        run_bash_cmd(ichor_cmd)
        
        # Step D: Parse & Store Results
        params_file = os.path.join(cfg_out_dir, f"{SAMPLE_ID}.params.txt")
        stats = parse_params_summary(params_file)
        results.append({
            'Config': cfg_name,
            'Bin_Size': f"{window_size // 1000}kb" if window_size < 1000000 else f"{window_size // 1000000}Mb",
            'Mode': "Paired + PoN" if (use_normal and use_pon) else ("Paired Only" if use_normal else ("Tumor-Only PoN" if use_pon else "Tumor-Only Default")),
            'PoN': "Custom 1Mb" if use_pon else "Default/None",
            'MAD_Score': stats['mad'],
            'Tumor_Fraction': stats['purity'],
            'Ploidy': stats['ploidy'],
            'Params_Path': params_file
        })

    # Summary Table Output
    print("\n" + "="*95)
    print(f" BENCHMARK RESULTS FOR {SAMPLE_ID}")
    print("="*95)
    header = f"{'Config':<22} | {'Bin Size':<10} | {'Mode':<18} | {'PoN Used':<14} | {'MAD Score':<12} | {'Tumor Fraction':<15}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(f"{r['Config']:<22} | {r['Bin_Size']:<10} | {r['Mode']:<18} | {r['PoN']:<14} | {r['MAD_Score']:<12} | {r['Tumor_Fraction']:<15}")
    print("="*95)
    
    # Write TSV report
    summary_tsv = os.path.join(OUT_DIR, "bin_size_benchmark_summary.tsv")
    with open(summary_tsv, 'w') as f_out:
        f_out.write("Config\tBin_Size\tMode\tPoN_Used\tMAD_Score\tTumor_Fraction\tPloidy\tParams_Path\n")
        for r in results:
            f_out.write(f"{r['Config']}\t{r['Bin_Size']}\t{r['Mode']}\t{r['PoN']}\t{r['MAD_Score']}\t{r['Tumor_Fraction']}\t{r['Ploidy']}\t{r['Params_Path']}\n")
    print(f"\nSummary table saved to: {summary_tsv}")

if __name__ == "__main__":
    main()
