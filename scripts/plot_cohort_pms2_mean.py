#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from concurrent.futures import ThreadPoolExecutor

# Cohort Results Dir
results_dir = "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results"

# Exon coordinates (exons 1 to 15, transcribed on minus strand)
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

# Variants to annotate
variants = [
    {
        "name": "p.Ile853Met (c.2559C>G) [VUS]\nUHN-27682 (Exon 15)\nFiltered: Tumour 74 | Normal 55",
        "pos": 5973429,
        "color": "#1e3a8a",
        "y_offset_linear": 0.4,
        "y_offset_log": 0.25,
        "side": "left"
    },
    {
        "name": "p.Glu836fs (c.2506delG) [Pathogenic]\nUHN-105729 (Exon 15)\nFiltered: Tumour 2,863 | Normal 2,528",
        "pos": 5973481,
        "color": "#7f1d1d",
        "y_offset_linear": 0.75,
        "y_offset_log": 0.65,
        "side": "right"
    },
    {
        "name": "p.Leu822= (c.2466T>C) [Benign]\nUHN-56254 (Exon 15)\nFiltered: Tumour 49 | Normal 56",
        "pos": 5973522,
        "color": "#b45309",
        "y_offset_linear": 0.22,
        "y_offset_log": 0.12,
        "side": "right"
    },
    {
        "name": "p.Gln374fs (c.1119_1122delTCAG) [Pathogenic]\nUHN-96803 (Exon 10)\nFiltered: Tumour 610 | Normal 748",
        "pos": 5989821,
        "color": "#065f46",
        "y_offset_linear": 0.55,
        "y_offset_log": 0.45,
        "side": "right"
    }
]

chrom = "chr7"
start = 5970000
end = 6015000

# 1. Identify all tumor BAMs
print("Identifying tumor BAM files...")
tumor_bams = []
# Blacklist the single BAM that timed out during diagnostics to prevent system hang
blacklist = ["DYNA-OTB-716535-P-T-DNA-U_S295"]

for p in os.listdir(results_dir):
    if p.startswith("DYNA-"):
        p_path = os.path.join(results_dir, p)
        if os.path.isdir(p_path):
            for s in os.listdir(p_path):
                if s in blacklist:
                    continue
                if "-T-" in s or "-T_S" in s or "-T-DNA" in s or "-P-T-" in s:
                    bam_file = os.path.join(p_path, s, "alignments", "dna", f"{s}.redux.bam")
                    if os.path.exists(bam_file):
                        tumor_bams.append((s, bam_file))

print(f"Found {len(tumor_bams)} tumor BAM files to query.")

# 2. Parallel BAM Coverage Extraction with strict timeout
import sys

def extract_single_bam_depth(item):
    sample_id, bam_path = item
    try:
        raw_cmd = f"samtools depth -a -r {chrom}:{start}-{end} {bam_path}"
        # Use timeout of 45 seconds to allow slow network reads to complete
        raw_out = subprocess.check_output(raw_cmd, shell=True, timeout=45).decode("utf-8")
        
        mapq_cmd = f"samtools depth -a -Q 1 -r {chrom}:{start}-{end} {bam_path}"
        mapq_out = subprocess.check_output(mapq_cmd, shell=True, timeout=45).decode("utf-8")
        
        raw_lines = raw_out.strip().split("\n")
        mapq_lines = mapq_out.strip().split("\n")
        
        if len(raw_lines) != len(mapq_lines):
            return None
            
        raw_depth = []
        mapq_depth = []
        for rl, ml in zip(raw_lines, mapq_lines):
            if not rl or not ml: continue
            raw_depth.append(int(rl.split("\t")[2]))
            mapq_depth.append(int(ml.split("\t")[2]))
            
        print(f"  {sample_id}: OK", flush=True)
        return np.array(raw_depth), np.array(mapq_depth)
    except subprocess.TimeoutExpired:
        print(f"  {sample_id}: TIMEOUT (skipped)", flush=True)
        return None
    except Exception as e:
        print(f"  {sample_id}: ERROR ({e})", flush=True)
        return None

# Extract coverages in parallel (4 workers to avoid thrashing the mount)
print("Querying coverages across cohort (max 4 concurrent workers)...", flush=True)
sum_raw = np.zeros(end - start + 1)
sum_mapq = np.zeros(end - start + 1)
count = 0

with ThreadPoolExecutor(max_workers=4) as executor:
    results = executor.map(extract_single_bam_depth, tumor_bams)
    for res in results:
        if res is not None:
            raw_arr, mapq_arr = res
            if len(raw_arr) == len(sum_raw):
                sum_raw += raw_arr
                sum_mapq += mapq_arr
                count += 1

print(f"Successfully aggregated coverage for {count} out of {len(tumor_bams)} samples.", flush=True)
if count == 0:
    print("Error: No coverages were successfully aggregated.", flush=True)
    exit(1)

mean_raw = sum_raw / count
mean_mapq = sum_mapq / count
positions_raw = np.arange(start, end + 1)



# Build dataframe
df = pd.DataFrame({
    "Position": positions_raw,
    "Depth_Raw": mean_raw,
    "Depth_MAPQ": mean_mapq
})

# Bin data into 100bp windows for smoothing
bin_size = 100
df["Bin"] = (df["Position"] - start) // bin_size
df_binned = df.groupby("Bin").agg({
    "Position": "mean",
    "Depth_Raw": "mean",
    "Depth_MAPQ": "mean"
}).reset_index()

# Load homology profile data
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
homology_tsv = os.path.join(script_dir, "pms2_homology_profile.tsv")
if os.path.exists(homology_tsv):
    print(f"Loading homology profile from {homology_tsv}...")
    df_homology = pd.read_csv(homology_tsv, sep="\t")
else:
    print(f"Warning: Homology profile TSV not found at {homology_tsv}")
    df_homology = None

# 4. Generate Linear and Log Plots
# 4. Generate Linear and Log Plots
for scale in ["linear", "log"]:
    print(f"Generating cohort plot on {scale} scale...")
    plot_df = df_binned.copy()
    
    if scale == "log":
        plot_df["Depth_Raw"] = plot_df["Depth_Raw"].clip(lower=1)
        plot_df["Depth_MAPQ"] = plot_df["Depth_MAPQ"].clip(lower=1)

        # Plot setup - 24 inches wide for cohort track log-scale
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(24, 6.5), sharex=True, 
                                       gridspec_kw={'height_ratios': [4, 1.2]})

        x_mb = plot_df["Position"] / 1e6

        # 1. Coverage Track
        ax1.fill_between(x_mb, plot_df["Depth_Raw"], label="Cohort Mean Raw Coverage (All Reads)", 
                         color="#fca5a5", alpha=0.6)
        ax1.fill_between(x_mb, plot_df["Depth_MAPQ"], label="Cohort Mean Unique Coverage (MAPQ > 0)", 
                         color="#dc2626", alpha=0.8)

        ax1.set_ylabel("Mean Read Depth (Coverage, Log Scale)", fontsize=11, fontweight="bold")
        ax1.set_title(f"PMS2 Cohort-Wide Mean Coverage Drop & Variant Annotations (N={count} Tumours) - Log Scale", 
                      fontsize=13, fontweight="bold", pad=15)
        ax1.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none")
        
        max_raw_depth = max(plot_df["Depth_Raw"])
        ax1.grid(True, which="both", linestyle=":", alpha=0.6)
        ax1.set_yscale("log")
        ax1.set_ylim(1, max_raw_depth * 3)

        # Annotate all variants
        for var in variants:
            var_x = var["pos"] / 1e6
            ax1.axvline(var_x, color=var["color"], linestyle="--", linewidth=1.1, alpha=0.9, zorder=3)
            y_pos = 10 ** (var["y_offset_log"] * np.log10(max_raw_depth))
            x_text_offset = -0.0055 if var["side"] == "left" else 0.0008
            
            ax1.annotate(var["name"], 
                         xy=(var_x, y_pos), 
                         xytext=(var_x + x_text_offset, y_pos),
                         fontsize=8.0, fontweight="bold",
                         arrowprops=dict(facecolor=var["color"], arrowstyle="-|>", lw=0.9, alpha=0.8),
                         bbox=dict(boxstyle="round,pad=0.35", fc="#f8fafc", ec=var["color"], lw=0.9, alpha=0.95),
                         zorder=4)

        # 2. Exon & Homology Tracks (ax2)
        ax2.plot([start/1e6, end/1e6], [0.65, 0.65], color="gray", linewidth=1.2, zorder=1)

        # Draw exons
        for exon in exons:
            e_start = exon["start"] / 1e6
            e_end = exon["end"] / 1e6
            width = e_end - e_start
            visual_width = max(width, 0.0003) 
            
            color = "#dc2626" if exon["num"] >= 11 else "#2563eb"
            rect = patches.Rectangle((e_start, 0.55), visual_width, 0.2, 
                                     edgecolor="black", facecolor=color, linewidth=0.5, zorder=2)
            ax2.add_patch(rect)
            
            y_text = 0.82 if exon["num"] % 2 == 0 else 0.43
            ax2.text(e_start + visual_width/2, y_text, f"E{exon['num']}", 
                     ha="center", va="center", fontsize=7.5, fontweight="bold", color="#111111")

        # Draw PMS2CL Homology Track underneath
        ax2.plot([start/1e6, end/1e6], [0.15, 0.15], color="gray", linewidth=1.0, linestyle=":", zorder=1)
        if df_homology is not None:
            x_hom = df_homology["Position"] / 1e6
            y_hom = 0.08 + 0.14 * df_homology["Homology_Score"]
            ax2.fill_between(x_hom, 0.08, y_hom, color="#fca5a5", alpha=0.8, edgecolor="#dc2626", linewidth=0.8, zorder=2)
            ax2.text((5.970 + 5.989)/2, 0.28, "PMS2CL Homology Profile", 
                     ha="center", va="center", fontsize=9, fontweight="bold", color="#dc2626")
        else:
            hom_start = 5.970
            hom_end = 5.989
            hom_rect = patches.Rectangle((hom_start, 0.08), hom_end - hom_start, 0.14, 
                                         edgecolor="black", facecolor="#fca5a5", hatch="//", alpha=0.7, zorder=2)
            ax2.add_patch(hom_rect)
            ax2.text((hom_start + hom_end)/2, 0.28, "PMS2CL Pseudogene Homology (Exons 11-15)", 
                     ha="center", va="center", fontsize=9, fontweight="bold", color="#dc2626")

        ax2.set_ylim(-0.05, 1.05)
        ax2.set_yticks([0.15, 0.65])
        ax2.set_yticklabels(["Homology", "PMS2 Exons"], fontsize=9.5, fontweight="bold")
        ax2.set_xlabel("Genomic Position on Chromosome 7 (Mb)", fontsize=11, fontweight="bold")
        ax2.set_xlim(start/1e6, end/1e6)
        ax2.grid(True, axis='x', linestyle=":", alpha=0.6)

        plt.tight_layout()

    else:
        # Split X-axis plot setup for cohort linear-scale
        exons_sorted = sorted(exons, key=lambda x: x["start"])
        num_exons = len(exons_sorted)
        
        padding = 150 # bp
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
        
        max_raw_depth = df["Depth_Raw"].max()
        raw_handle = None
        mapq_handle = None
        
        for idx, exon in enumerate(exons_sorted):
            ax_cov = axes_cov[idx]
            ax_ex = axes_exon[idx]
            
            exon_start_bp = exon["start"] - padding
            exon_end_bp = exon["end"] + padding
            
            # Slice data for this exon
            df_slice = df[(df["Position"] >= exon_start_bp) & (df["Position"] <= exon_end_bp)]
            x_mb = df_slice["Position"] / 1e6
            
            # Plot coverage
            h_raw = ax_cov.fill_between(x_mb, df_slice["Depth_Raw"], color="#fca5a5", alpha=0.6)
            h_mapq = ax_cov.fill_between(x_mb, df_slice["Depth_MAPQ"], color="#dc2626", alpha=0.8)
            
            if idx == 0:
                raw_handle = h_raw
                mapq_handle = h_mapq
                
            ax_cov.set_ylim(-max_raw_depth * 0.05, max_raw_depth * 1.15)
            
            # Annotate all variants that fall into this exon (with layout offsets to prevent overlap in E15)
            for var in variants:
                if exon_start_bp <= var["pos"] <= exon_end_bp:
                    var_x = var["pos"] / 1e6
                    ax_cov.axvline(var_x, color=var["color"], linestyle="--", linewidth=1.2, alpha=0.9, zorder=3)
                    
                    if exon["num"] == 15:
                        if "Ile853Met" in var["name"]:
                            xytext = (0.02, 1.12)
                            ha = 'left'
                        elif "Glu836fs" in var["name"]:
                            xytext = (0.5, 1.72)
                            ha = 'center'
                        else:  # Leu822=
                            xytext = (0.98, 1.42)
                            ha = 'right'
                    else:
                        xytext = (0.5, 1.12)
                        ha = 'center'
                        
                    ax_cov.annotate(
                        var["name"],
                        xy=(var_x, max_raw_depth * 0.95),
                        xytext=xytext,
                        textcoords='axes fraction',
                        fontsize=7.0, fontweight="bold",
                        ha=ha, va='bottom',
                        arrowprops=dict(facecolor=var["color"], arrowstyle="-|>", lw=1.0, alpha=0.8),
                        bbox=dict(boxstyle="round,pad=0.25", fc="#f8fafc", ec=var["color"], lw=1.0, alpha=0.95),
                        zorder=4
                    )
            
            ax_cov.grid(True, linestyle=":", alpha=0.6)
            if idx > 0:
                ax_cov.tick_params(axis='y', which='both', left=False, labelleft=False)
                ax_cov.spines['left'].set_visible(False)
            else:
                ax_cov.set_ylabel("Mean Read Depth (Coverage)", fontsize=11, fontweight="bold")
                
            # Exons track subplot
            xmin = exon_start_bp / 1e6
            xmax = exon_end_bp / 1e6
            ax_ex.plot([xmin, xmax], [0.65, 0.65], color="gray", linewidth=1.2, zorder=1)
            
            e_start = exon["start"] / 1e6
            e_end = exon["end"] / 1e6
            e_w = e_end - e_start
            
            e_color = "#dc2626" if exon["num"] >= 11 else "#2563eb"
            rect = patches.Rectangle((e_start, 0.55), e_w, 0.2, edgecolor="black", facecolor=e_color, linewidth=0.5, zorder=2)
            ax_ex.add_patch(rect)
            
            ax_ex.text(e_start + e_w/2, 0.82, f"E{exon['num']}", ha="center", va="center", fontsize=8, fontweight="bold", color="#111111")
            
            # Pseudogene homology track
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
            
        fig.suptitle(f"PMS2 Cohort-Wide Mean Coverage Profile & Variant Annotations (N={count} Tumours) - Split X-Axis (Linear Scale)", fontsize=14, fontweight="bold", y=0.98)
        fig.legend(handles=[raw_handle, mapq_handle], labels=["Raw BAM Coverage", "Unique BAM Coverage (MAPQ > 0)"], loc="upper right", bbox_to_anchor=(0.99, 0.95))
        
        # Increase top margin to fit the three stacked labels above E15
        plt.subplots_adjust(left=0.04, right=0.99, top=0.66, bottom=0.12, wspace=0.08)

    out_png = f"/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/pms2_cohort_mean_coverage_plot_{scale}.png"
    out_pdf = f"/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/pms2_cohort_mean_coverage_plot_{scale}.pdf"
    
    print(f"Saving cohort plot to: {out_png}")
    plt.savefig(out_png, dpi=300)
    print(f"Saving cohort plot to: {out_pdf}")
    plt.savefig(out_pdf)
    plt.close()

print("Cohort mean coverage plotting complete!")
