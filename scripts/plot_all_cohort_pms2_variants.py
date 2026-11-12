#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Default directories
DEFAULT_RESULTS_DIR = "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results"

def parse_args():
    parser = argparse.ArgumentParser(description="Auto-discover all PMS2 variants in the cohort and generate coverage/homology plots.")
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR, help="Base cohort results directory.")
    parser.add_argument("--output-dir", help="Directory to save generated plots (defaults to each sample's summary/plots/ directory).")
    return parser.parse_args()

def extract_depths(format_str, sample_str):
    fields = format_str.split(":")
    vals = sample_str.split(":")
    if "AD" in fields:
        ad_idx = fields.index("AD")
        if ad_idx < len(vals):
            ad_val = vals[ad_idx]
            if "," in ad_val:
                ref_str, alt_str = ad_val.split(",", 1)
                if "," in alt_str:
                    alt_str = alt_str.split(",")[0]
                try:
                    return int(ref_str), int(alt_str)
                except ValueError:
                    pass
    return None, None

def parse_pms2_variants(vcf_path, var_type):
    if not os.path.exists(vcf_path):
        return []
    variants = []
    try:
        with gzip.open(vcf_path, "rt", encoding="utf-8") as f:
            tumor_col_idx = 10
            normal_col_idx = 9
            for line in f:
                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    for idx, col in enumerate(cols[9:], start=9):
                        if "-T-" in col or "T-DNA" in col or col.endswith("-T"):
                            tumor_col_idx = idx
                        elif "-N-" in col or "N-DNA" in col or col.endswith("-N"):
                            normal_col_idx = idx
                    continue
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                filt = parts[6].strip()
                if filt != "PASS":
                    continue
                pos = int(parts[1])
                chrom = parts[0]
                
                # Check PMS2 coordinates (chr7:5970925-6009130)
                if chrom not in ["chr7", "7"] or not (5970925 <= pos <= 6009130):
                    continue
                
                # Parse INFO
                info_str = parts[7]
                info = {}
                for item in info_str.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info[k] = v
                    else:
                        info[item] = True
                
                # Check gene
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
                if gene != "PMS2":
                    continue
                    
                # Get HGVS change
                change = ""
                if 'PAVE_TI' in info:
                    transcripts = info['PAVE_TI'].split(',')
                    for tx in transcripts:
                        tx_parts = tx.split('|')
                        if len(tx_parts) > 6:
                            is_canonical = tx_parts[4].strip() in ['true', '1']
                            if is_canonical:
                                change = tx_parts[6].strip() or tx_parts[5].strip()
                                if change:
                                    break
                    if not change:
                        for tx in transcripts:
                            tx_parts = tx.split('|')
                            if len(tx_parts) > 6:
                                change = tx_parts[6].strip() or tx_parts[5].strip()
                                if change:
                                    break
                
                # Get VAF
                vaf = None
                format_fields = parts[8].split(':')
                if 'AF' in format_fields and len(parts) > tumor_col_idx:
                    af_idx = format_fields.index('AF')
                    tumor_vals = parts[tumor_col_idx].split(':')
                    if len(tumor_vals) > af_idx:
                        try:
                            vaf = float(tumor_vals[af_idx].split(',')[0])
                        except ValueError:
                            pass
                
                # Extract read depths
                t_ref, t_alt = None, None
                n_ref, n_alt = None, None
                if len(parts) > tumor_col_idx:
                    t_ref, t_alt = extract_depths(parts[8], parts[tumor_col_idx])
                if len(parts) > normal_col_idx:
                    n_ref, n_alt = extract_depths(parts[8], parts[normal_col_idx])
                
                variants.append({
                    "pos": pos,
                    "change": change if change else f"g.{pos}",
                    "vaf": vaf,
                    "t_ref": t_ref,
                    "t_alt": t_alt,
                    "n_ref": n_ref,
                    "n_alt": n_alt,
                    "var_type": var_type
                })
    except Exception as e:
        print(f"Error parsing VCF {vcf_path}: {e}")
    return variants

def main():
    args = parse_args()
    results_dir = os.path.abspath(args.results_dir)
    
    if not os.path.exists(results_dir):
        print(f"Error: Results directory {results_dir} does not exist.")
        sys.exit(1)
        
    script_dir = os.path.dirname(os.path.abspath(__file__))
    homology_tsv = os.path.join(script_dir, "pms2_homology_profile.tsv")
    if os.path.exists(homology_tsv):
        print(f"Loading homology profile from {homology_tsv}...")
        df_homology = pd.read_csv(homology_tsv, sep="\t")
    else:
        print(f"Warning: Homology profile TSV not found at {homology_tsv}. Homology tracks will revert to static representation.")
        df_homology = None
        
    chrom = "chr7"
    start = 5970000
    end = 6015000
    
    exons = [
        {"num": 1, "start": 6008997, "end": 6009049},
        {"num": 2, "start": 6005892, "end": 6006031},
        {"num": 3, "start": 6003690, "end": 6003792},
        {"num": 4, "start": 6003972, "end": 6004058},
        {"num": 5, "start": 6002453, "end": 6002636},
        {"num": 6, "start": 5998987, "end": 5999275},
        {"num": 7, "start": 5997211, "end": 5997423},
        {"num": 8, "start": 5995534, "end": 5995633},
        {"num": 9, "start": 5991973, "end": 5992057},
        {"num": 10, "start": 5989800, "end": 5989955},
        {"num": 11, "start": 5986759, "end": 5987620},
        {"num": 12, "start": 5982824, "end": 5982991},
        {"num": 13, "start": 5978596, "end": 5978696},
        {"num": 14, "start": 5977588, "end": 5977757},
        {"num": 15, "start": 5973399, "end": 5973542}
    ]
    
    # Discover samples and variants
    discovered_variants = []
    print("Scanning results directories for PMS2 variants...")
    
    for item in sorted(os.listdir(results_dir)):
        item_path = os.path.join(results_dir, item)
        if not os.path.isdir(item_path):
            continue
        if not (item.startswith("DYNA-") or item.startswith("Seraseq-")):
            continue
            
        subfolders = [d for d in os.listdir(item_path) if os.path.isdir(os.path.join(item_path, d)) and d.startswith(item)]
        if not subfolders:
            continue
        sample_subfolder = os.path.join(item_path, subfolders[0])
        group_id = subfolders[0]
        
        # Check BAM exists
        bam_dir = os.path.join(sample_subfolder, "alignments", "dna")
        if not os.path.exists(bam_dir):
            continue
        bams = [f for f in os.listdir(bam_dir) if f.endswith(".redux.bam")]
        if not bams:
            continue
        bam_path = os.path.join(bam_dir, bams[0])
        
        # Check VCF files
        pave_dir = os.path.join(sample_subfolder, "pave")
        somatic_vcf = os.path.join(pave_dir, f"{group_id}.pave.somatic.vcf.gz")
        germline_vcf = os.path.join(pave_dir, f"{group_id}.pave.germline.vcf.gz")
        
        sample_vars = []
        sample_vars.extend(parse_pms2_variants(somatic_vcf, "somatic"))
        sample_vars.extend(parse_pms2_variants(germline_vcf, "germline"))
        
        for sv in sample_vars:
            discovered_variants.append({
                "patient": item,
                "group_id": group_id,
                "bam": bam_path,
                "pos": sv["pos"],
                "change": sv["change"],
                "vaf": sv["vaf"],
                "t_ref": sv["t_ref"],
                "t_alt": sv["t_alt"],
                "n_ref": sv["n_ref"],
                "n_alt": sv["n_alt"],
                "var_type": sv["var_type"],
                "summary_dir": os.path.join(item_path, "summary")
            })
            
    print(f"Found {len(discovered_variants)} PMS2 variants to plot.")
    
    # Process and plot each variant
    for i, var in enumerate(discovered_variants, 1):
        patient = var["patient"]
        group_id = var["group_id"]
        pos = var["pos"]
        change = var["change"]
        var_type = var["var_type"]
        bam_path = var["bam"]
        
        print(f"[{i}/{len(discovered_variants)}] Plotting {change} in sample {group_id}...")
        
        # Define output directory
        out_dir = args.output_dir if args.output_dir else os.path.join(var["summary_dir"], "plots")
        os.makedirs(out_dir, exist_ok=True)
        
        # Get raw and MAPQ>0 depth
        raw_cmd = f"samtools depth -a -r {chrom}:{start}-{end} {bam_path}"
        raw_output = subprocess.check_output(raw_cmd, shell=True).decode("utf-8")
        
        mapq_cmd = f"samtools depth -a -Q 1 -r {chrom}:{start}-{end} {bam_path}"
        mapq_output = subprocess.check_output(mapq_cmd, shell=True).decode("utf-8")
        
        def parse_depth(output_str):
            positions, depths = [], []
            for line in output_str.strip().split("\n"):
                if not line: continue
                parts = line.split("\t")
                positions.append(int(parts[1]))
                depths.append(int(parts[2]))
            return pd.DataFrame({"Position": positions, "Depth": depths})
            
        df_raw = parse_depth(raw_output)
        df_mapq = parse_depth(mapq_output)
        df_merged = pd.merge(df_raw, df_mapq, on="Position", suffixes=("_Raw", "_MAPQ"))
        
        # Smooth data
        bin_size = 100
        df_merged["Bin"] = (df_merged["Position"] - start) // bin_size
        df_binned = df_merged.groupby("Bin").agg({
            "Position": "mean",
            "Depth_Raw": "mean",
            "Depth_MAPQ": "mean"
        }).reset_index()
        
        # Create labels
        var_name = f"{var_type.capitalize()} PMS2 {change}"
        vaf_val = var["vaf"] * 100 if var["vaf"] else 0
        cov_tum = (var["t_ref"] + var["t_alt"]) if (var["t_ref"] is not None and var["t_alt"] is not None) else 0
        cov_norm = (var["n_ref"] + var["n_alt"]) if (var["n_ref"] is not None and var["n_alt"] is not None) else 0
        
        label_text = f"{var_name}\n"
        if vaf_val > 0:
            label_text += f"VAF: {vaf_val:.1f}%\n"
        label_text += f"Filtered Exon Coverage:\nTumour: {cov_tum} | Normal: {cov_norm}"
        
        color = "#7f1d1d" if var_type == "germline" else "#1e3a8a"
        
        # Generate Log and Linear Plots
        for scale in ["linear", "log"]:
            plot_df = df_binned.copy()
            if scale == "log":
                plot_df["Depth_Raw"] = plot_df["Depth_Raw"].clip(lower=1)
                plot_df["Depth_MAPQ"] = plot_df["Depth_MAPQ"].clip(lower=1)
                
            max_raw_depth = plot_df["Depth_Raw"].max()
            
            if scale == "log":
                # Continuous plot setup
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(24, 7), sharex=True, gridspec_kw={'height_ratios': [4, 1.2]})
                
                # 1. Coverage Plot (ax1)
                ax1.fill_between(plot_df["Position"]/1e6, plot_df["Depth_Raw"], color="#fca5a5", alpha=0.6, label="Raw BAM Coverage")
                ax1.fill_between(plot_df["Position"]/1e6, plot_df["Depth_MAPQ"], color="#dc2626", alpha=0.8, label="Unique BAM Coverage (MAPQ > 0)")
                ax1.set_yscale('log')
                ax1.legend(loc="upper right", fontsize=10)
                ax1.set_ylabel("Read Depth (Coverage, Log Scale)", fontsize=11, fontweight="bold")
                ax1.grid(True, which="both", linestyle=":", alpha=0.5)
                ax1.set_ylim(1, max_raw_depth * 3)
                
                # Annotate variant
                var_x = pos / 1e6
                ax1.axvline(var_x, color=color, linestyle="--", linewidth=1.2, alpha=0.9, zorder=3)
                y_pos = 10 ** (0.6 * np.log10(max_raw_depth))
                ax1.annotate(label_text, xy=(var_x, y_pos), xytext=(var_x + 0.001, y_pos),
                             fontsize=8.5, fontweight="bold",
                             arrowprops=dict(facecolor=color, arrowstyle="-|>", lw=1.0, alpha=0.8),
                             bbox=dict(boxstyle="round,pad=0.4", fc="#f8fafc", ec=color, lw=1.0, alpha=0.9),
                             zorder=4)
                             
                # 2. Exon & Homology Tracks (ax2)
                ax2.plot([start/1e6, end/1e6], [0.65, 0.65], color="gray", linewidth=1.2, zorder=1)
                for exon in exons:
                    e_start = exon["start"] / 1e6
                    e_end = exon["end"] / 1e6
                    e_w = max(e_end - e_start, 0.0003)
                    e_color = "#dc2626" if exon["num"] >= 11 else "#2563eb"
                    rect = patches.Rectangle((e_start, 0.55), e_w, 0.2, edgecolor="black", facecolor=e_color, linewidth=0.5, zorder=2)
                    ax2.add_patch(rect)
                    y_text = 0.82 if exon["num"] % 2 == 0 else 0.43
                    ax2.text(e_start + e_w/2, y_text, f"E{exon['num']}", ha="center", va="center", fontsize=7.5, fontweight="bold", color="#111111")
                    
                ax2.plot([start/1e6, end/1e6], [0.15, 0.15], color="gray", linewidth=1.0, linestyle=":", zorder=1)
                if df_homology is not None:
                    x_hom = df_homology["Position"] / 1e6
                    y_hom = 0.08 + 0.14 * df_homology["Homology_Score"]
                    ax2.fill_between(x_hom, 0.08, y_hom, color="#fca5a5", alpha=0.8, edgecolor="#dc2626", linewidth=0.8, zorder=2)
                    ax2.text((5.970 + 5.989)/2, 0.28, "PMS2CL Homology Profile", ha="center", va="center", fontsize=9, fontweight="bold", color="#dc2626")
                else:
                    hom_start = 5.970
                    hom_end = 5.989
                    hom_rect = patches.Rectangle((hom_start, 0.08), hom_end - hom_start, 0.14, edgecolor="black", facecolor="#fca5a5", hatch="//", alpha=0.7, zorder=2)
                    ax2.add_patch(hom_rect)
                    ax2.text((hom_start + hom_end)/2, 0.28, "PMS2CL Pseudogene Homology (Exons 11-15)", ha="center", va="center", fontsize=9, fontweight="bold", color="#dc2626")
                    
                ax2.set_ylim(-0.05, 1.05)
                ax2.set_yticks([0.15, 0.65])
                ax2.set_yticklabels(["Homology", "PMS2 Exons"], fontsize=9.5, fontweight="bold")
                ax2.set_xlabel("Genomic Position on Chromosome 7 (Mb)", fontsize=11, fontweight="bold")
                ax2.set_xlim(start/1e6, end/1e6)
                ax2.grid(True, axis='x', linestyle=":", alpha=0.6)
                
                plt.subplots_adjust(top=0.92, bottom=0.10, left=0.06, right=0.98, hspace=0.25)
                fig.suptitle(f"PMS2 Homology & Coverage: {group_id} - Continuous {var_name} (Log Scale)", fontsize=14, fontweight="bold")
                
            else:
                # Split x-axis linear scale plot
                exons_sorted = sorted(exons, key=lambda x: x["start"])
                num_exons = len(exons_sorted)
                
                padding = 150
                widths = []
                for exon in exons_sorted:
                    w = (exon["end"] - exon["start"]) + 2 * padding
                    widths.append(w)
                    
                fig, (axes_cov, axes_exon) = plt.subplots(
                    2, num_exons,
                    sharey="row",
                    sharex="col",
                    figsize=(24, 7.5),
                    gridspec_kw={'width_ratios': widths, 'height_ratios': [4, 1.2]}
                )
                
                raw_handle, mapq_handle = None, None
                for idx, exon in enumerate(exons_sorted):
                    ax_cov = axes_cov[idx]
                    ax_ex = axes_exon[idx]
                    
                    exon_start_bp = exon["start"] - padding
                    exon_end_bp = exon["end"] + padding
                    
                    # Slice data
                    df_slice = df_merged[(df_merged["Position"] >= exon_start_bp) & (df_merged["Position"] <= exon_end_bp)]
                    x_mb = df_slice["Position"] / 1e6
                    
                    h_raw = ax_cov.fill_between(x_mb, df_slice["Depth_Raw"], color="#fca5a5", alpha=0.6)
                    h_mapq = ax_cov.fill_between(x_mb, df_slice["Depth_MAPQ"], color="#dc2626", alpha=0.8)
                    if idx == 0:
                        raw_handle, mapq_handle = h_raw, h_mapq
                        
                    ax_cov.set_xlim(exon_start_bp/1e6, exon_end_bp/1e6)
                    ax_cov.set_ylim(0, max_raw_depth * 1.1)
                    
                    # Annotate variant if it falls within this exon's boundaries
                    if exon_start_bp <= pos <= exon_end_bp:
                        var_x = pos / 1e6
                        ax_cov.axvline(var_x, color=color, linestyle="--", linewidth=1.2, alpha=0.9, zorder=3)
                        ax_cov.annotate(label_text, xy=(var_x, max_raw_depth * 0.5), xytext=(var_x + 0.0001, max_raw_depth * 0.5),
                                         fontsize=8.0, fontweight="bold",
                                         arrowprops=dict(facecolor=color, arrowstyle="-|>", lw=1.0, alpha=0.8),
                                         bbox=dict(boxstyle="round,pad=0.25", fc="#f8fafc", ec=color, lw=1.0, alpha=0.95),
                                         zorder=4)
                                         
                    ax_cov.grid(True, linestyle=":", alpha=0.6)
                    if idx > 0:
                        ax_cov.tick_params(axis='y', which='both', left=False, labelleft=False)
                        ax_cov.spines['left'].set_visible(False)
                    else:
                        ax_cov.set_ylabel("Read Depth (Coverage)", fontsize=11, fontweight="bold")
                        
                    # Exon structure subplot
                    xmin, xmax = exon_start_bp / 1e6, exon_end_bp / 1e6
                    ax_ex.plot([xmin, xmax], [0.65, 0.65], color="gray", linewidth=1.2, zorder=1)
                    
                    e_start, e_end = exon["start"] / 1e6, exon["end"] / 1e6
                    e_w = e_end - e_start
                    e_color = "#dc2626" if exon["num"] >= 11 else "#2563eb"
                    rect = patches.Rectangle((e_start, 0.55), e_w, 0.2, edgecolor="black", facecolor=e_color, linewidth=0.5, zorder=2)
                    ax_ex.add_patch(rect)
                    ax_ex.text(e_start + e_w/2, 0.82, f"E{exon['num']}", ha="center", va="center", fontsize=8, fontweight="bold", color="#111111")
                    
                    # Homology profile track
                    ax_ex.plot([xmin, xmax], [0.15, 0.15], color="gray", linewidth=1.0, linestyle=":", zorder=1)
                    if df_homology is not None:
                        df_hom_slice = df_homology[(df_homology["Position"] >= exon_start_bp) & (df_homology["Position"] <= exon_end_bp)]
                        if not df_hom_slice.empty:
                            x_hom = df_hom_slice["Position"] / 1e6
                            y_hom = 0.08 + 0.14 * df_hom_slice["Homology_Score"]
                            ax_ex.fill_between(x_hom, 0.08, y_hom, color="#fca5a5", alpha=0.8, edgecolor="#dc2626", linewidth=0.8, zorder=2)
                        if exon["num"] == 13:
                            ax_ex.text((xmin + xmax)/2, 0.3, "PMS2CL Homology", ha="center", va="center", fontsize=8, fontweight="bold", color="#dc2626")
                    else:
                        if exon["num"] >= 11:
                            hom_rect = patches.Rectangle((xmin, 0.08), xmax - xmin, 0.14, edgecolor="black", facecolor="#fca5a5", hatch="//", alpha=0.7, zorder=2)
                            ax_ex.add_patch(hom_rect)
                            if exon["num"] == 13:
                                ax_ex.text((xmin + xmax)/2, 0.3, "PMS2CL Homology", ha="center", va="center", fontsize=8, fontweight="bold", color="#dc2626")
                                
                    ax_ex.set_ylim(-0.05, 1.05)
                    if idx > 0:
                        ax_ex.tick_params(axis='y', which='both', left=False, labelleft=False)
                        ax_ex.spines['left'].set_visible(False)
                    else:
                        ax_ex.set_yticks([0.15, 0.65])
                        ax_ex.set_yticklabels(["Homology", "PMS2 Exons"], fontsize=9.5, fontweight="bold")
                        
                    ax_ex.tick_params(axis='x', labelsize=8, rotation=45)
                    ax_ex.xaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))
                    ax_ex.grid(True, linestyle=":", alpha=0.6)
                    
                if raw_handle and mapq_handle:
                    axes_cov[-1].legend([raw_handle, mapq_handle], ["Raw BAM Coverage", "Unique (MAPQ > 0)"], loc="upper right", fontsize=9)
                    
                plt.subplots_adjust(top=0.92, bottom=0.15, left=0.06, right=0.98, hspace=0.15, wspace=0.15)
                fig.suptitle(f"PMS2 Homology & Coverage: {group_id} - Split Exons {var_name} (Linear Scale)", fontsize=14, fontweight="bold")
                
            # Save files
            clean_change = change.replace(">", "to").replace("*", "star").replace("/", "_")
            out_prefix = os.path.join(out_dir, f"pms2_coverage_homology_plot_{group_id}_{clean_change}_{scale}")
            
            plt.savefig(f"{out_prefix}.png", dpi=200)
            plt.savefig(f"{out_prefix}.pdf", format="pdf")
            plt.close()
            print(f"  -> Saved {scale} plots to: {out_prefix}.png/.pdf")

    print("\nAll tasks completed successfully!")

if __name__ == "__main__":
    main()
