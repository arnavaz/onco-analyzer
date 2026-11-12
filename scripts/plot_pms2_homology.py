#!/usr/bin/env python3
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Configuration for all 4 samples with called PMS2 variants
samples = {
    "UHN-27682": {
        "bam": "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/DYNA-UHN-27682/DYNA-UHN-27682-FFPE-T-DNA-U_S53/alignments/dna/DYNA-UHN-27682-FFPE-T-DNA-U_S53.redux.bam",
        "var_name": "p.Ile853Met (c.2559C>G) [VUS]",
        "pos": 5973429,
        "cov_tum": 74,
        "cov_norm": 55,
        "color": "#1e3a8a",
        "y_offset_linear": 0.4,
        "y_offset_log": 0.3
    },
    "UHN-105729": {
        "bam": "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/DYNA-UHN-105729/DYNA-UHN-105729-FFPE-T-DNA-U_S31/alignments/dna/DYNA-UHN-105729-FFPE-T-DNA-U_S31.redux.bam",
        "var_name": "p.Glu836fs (c.2506delG) [Pathogenic]",
        "pos": 5973481,
        "cov_tum": 2863,
        "cov_norm": 2528,
        "color": "#7f1d1d",
        "y_offset_linear": 0.65,
        "y_offset_log": 0.6
    },
    "UHN-96803": {
        "bam": "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/DYNA-UHN-96803/DYNA-UHN-96803-FFPE-T-DNA-U_S43/alignments/dna/DYNA-UHN-96803-FFPE-T-DNA-U_S43.redux.bam",
        "var_name": "p.Gln374fs (c.1119_1122delTCAG) [Pathogenic]",
        "pos": 5989821,
        "cov_tum": 610,
        "cov_norm": 748,
        "color": "#065f46",
        "y_offset_linear": 0.55,
        "y_offset_log": 0.5
    },
    "UHN-56254": {
        "bam": "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/DYNA-UHN-56254/DYNA-UHN-56254-FFPE-T-DNA-U_S55/alignments/dna/DYNA-UHN-56254-FFPE-T-DNA-U_S55.redux.bam",
        "var_name": "p.Leu822= (c.2466T>C) [Benign]",
        "pos": 5973522,
        "cov_tum": 49,
        "cov_norm": 56,
        "color": "#b45309",
        "y_offset_linear": 0.4,
        "y_offset_log": 0.3
    }
}

chrom = "chr7"
start = 5970000
end = 6015000

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

# Load homology profile data
script_dir = os.path.dirname(os.path.abspath(__file__))
homology_tsv = os.path.join(script_dir, "pms2_homology_profile.tsv")
if os.path.exists(homology_tsv):
    print(f"Loading homology profile from {homology_tsv}...")
    df_homology = pd.read_csv(homology_tsv, sep="\t")
else:
    print(f"Warning: Homology profile TSV not found at {homology_tsv}")
    df_homology = None

# Run plotting for each sample
for sample_name, info in samples.items():
    bam_path = info["bam"]
    if not os.path.exists(bam_path):
        print(f"Warning: BAM not found for {sample_name}: {bam_path}")
        continue
        
    print(f"Processing data for {sample_name}...")
    
    # Extract raw coverage depth
    raw_cmd = f"samtools depth -a -r {chrom}:{start}-{end} {bam_path}"
    raw_output = subprocess.check_output(raw_cmd, shell=True).decode("utf-8")

    # Extract MAPQ > 0 coverage depth
    mapq_cmd = f"samtools depth -a -Q 1 -r {chrom}:{start}-{end} {bam_path}"
    mapq_output = subprocess.check_output(mapq_cmd, shell=True).decode("utf-8")

    # Parse depths
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

    # Merge
    df_merged = pd.merge(df_raw, df_mapq, on="Position", suffixes=("_Raw", "_MAPQ"))

    # Bin data into 100bp windows for smoothing
    bin_size = 100
    df_merged["Bin"] = (df_merged["Position"] - start) // bin_size
    df_binned = df_merged.groupby("Bin").agg({
        "Position": "mean",
        "Depth_Raw": "mean",
        "Depth_MAPQ": "mean"
    }).reset_index()

    # Generate both linear and log scale plots
    for scale in ["linear", "log"]:
        print(f"  Plotting {sample_name} on {scale} scale...")
        
        # Copy data so we can modify it
        plot_df = df_binned.copy()
        
        if scale == "log":
            # Clip depths to a minimum of 1 for robust log-scale plotting
            plot_df["Depth_Raw"] = plot_df["Depth_Raw"].clip(lower=1)
            plot_df["Depth_MAPQ"] = plot_df["Depth_MAPQ"].clip(lower=1)

            # Continuous axis plot setup for log-scale
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(24, 6.5), sharex=True, 
                                           gridspec_kw={'height_ratios': [4, 1.2]})

            x_mb = plot_df["Position"] / 1e6

            # 1. Coverage Track
            ax1.fill_between(x_mb, plot_df["Depth_Raw"], label="Raw BAM Coverage (All Reads)", 
                             color="#fca5a5", alpha=0.6)
            ax1.fill_between(x_mb, plot_df["Depth_MAPQ"], label="Unique BAM Coverage (MAPQ > 0)", 
                             color="#dc2626", alpha=0.8)

            ax1.set_ylabel("Read Depth (Coverage, Log Scale)", fontsize=11, fontweight="bold")
            ax1.set_title(f"PMS2 Coverage & Homology: {sample_name} (Tumour) - Log Scale", 
                          fontsize=13, fontweight="bold", pad=15)
            ax1.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none")
            
            max_raw_depth = max(plot_df["Depth_Raw"])
            ax1.grid(True, which="both", linestyle=":", alpha=0.6)
            ax1.set_yscale("log")
            ax1.set_ylim(1, max_raw_depth * 3)
            
            # Annotate specific variant for this patient
            var_x = info["pos"] / 1e6
            ax1.axvline(var_x, color=info["color"], linestyle="--", linewidth=1.2, alpha=0.9, zorder=3)
            
            label_text = f"{info['var_name']}\nFiltered Exon Coverage:\nTumour: {info['cov_tum']} | Normal: {info['cov_norm']}"
            y_pos = 10 ** (info["y_offset_log"] * np.log10(max_raw_depth))
            
            ax1.annotate(label_text, 
                         xy=(var_x, y_pos), 
                         xytext=(var_x + 0.001, y_pos),
                         fontsize=8.5, fontweight="bold",
                         arrowprops=dict(facecolor=info["color"], arrowstyle="-|>", lw=1.0, alpha=0.8),
                         bbox=dict(boxstyle="round,pad=0.4", fc="#f8fafc", ec=info["color"], lw=1.0, alpha=0.9),
                         zorder=4)

            # 2. Exon & Homology Tracks (ax2)
            ax2.plot([start/1e6, end/1e6], [0.65, 0.65], color="gray", linewidth=1.2, zorder=1)

            # Draw exons
            for exon in exons:
                e_start = exon["start"] / 1e6
                e_end = exon["end"] / 1e6
                width = e_end - e_start
                visual_width = max(width, 0.0003) 
                
                # Color homology exons (11-15) in red, unique exons (1-10) in blue
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
            # Split X-axis plot setup for linear-scale
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
                figsize=(24, 7), 
                gridspec_kw={'width_ratios': widths, 'height_ratios': [4, 1.2]}
            )
            
            max_raw_depth = df_merged["Depth_Raw"].max()
            raw_handle = None
            mapq_handle = None
            
            for idx, exon in enumerate(exons_sorted):
                ax_cov = axes_cov[idx]
                ax_ex = axes_exon[idx]
                
                exon_start_bp = exon["start"] - padding
                exon_end_bp = exon["end"] + padding
                
                # Slice data for this exon
                df_slice = df_merged[(df_merged["Position"] >= exon_start_bp) & (df_merged["Position"] <= exon_end_bp)]
                x_mb = df_slice["Position"] / 1e6
                
                # Plot coverage
                h_raw = ax_cov.fill_between(x_mb, df_slice["Depth_Raw"], color="#fca5a5", alpha=0.6)
                h_mapq = ax_cov.fill_between(x_mb, df_slice["Depth_MAPQ"], color="#dc2626", alpha=0.8)
                
                if idx == 0:
                    raw_handle = h_raw
                    mapq_handle = h_mapq
                    
                ax_cov.set_ylim(-max_raw_depth * 0.05, max_raw_depth * 1.15)
                
                # Plot vertical line and annotation if variant is here
                if exon_start_bp <= info["pos"] <= exon_end_bp:
                    var_x = info["pos"] / 1e6
                    ax_cov.axvline(var_x, color=info["color"], linestyle="--", linewidth=1.5, alpha=0.9, zorder=3)
                    
                    label_text = f"{info['var_name']}\nFiltered Exon Coverage:\nTumour: {info['cov_tum']} | Normal: {info['cov_norm']}"
                    ax_cov.annotate(
                        label_text,
                        xy=(var_x, max_raw_depth * 0.95),
                        xytext=(0.1, 1.12),
                        textcoords='axes fraction',
                        fontsize=8.0, fontweight="bold",
                        ha='left', va='bottom',
                        arrowprops=dict(facecolor=info["color"], arrowstyle="-|>", lw=1.2, alpha=0.8),
                        bbox=dict(boxstyle="round,pad=0.3", fc="#f8fafc", ec=info["color"], lw=1.2, alpha=0.9),
                        zorder=4
                    )
                
                ax_cov.grid(True, linestyle=":", alpha=0.6)
                if idx > 0:
                    ax_cov.tick_params(axis='y', which='both', left=False, labelleft=False)
                    ax_cov.spines['left'].set_visible(False)
                else:
                    ax_cov.set_ylabel("Read Depth (Coverage)", fontsize=11, fontweight="bold")
                    
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
                
            fig.suptitle(f"PMS2 Exon-Specific Coverage: {sample_name} (Tumour) - Split X-Axis (Linear Scale)", fontsize=14, fontweight="bold", y=0.98)
            fig.legend(handles=[raw_handle, mapq_handle], labels=["Raw BAM Coverage", "Unique BAM Coverage (MAPQ > 0)"], loc="upper right", bbox_to_anchor=(0.99, 0.95))
            
            plt.subplots_adjust(left=0.04, right=0.99, top=0.84, bottom=0.12, wspace=0.08)

        out_png_path = f"/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/pms2_coverage_homology_plot_{sample_name}_{scale}.png"
        out_pdf_path = f"/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/pms2_coverage_homology_plot_{sample_name}_{scale}.pdf"
        
        print(f"Saving PNG plot to: {out_png_path}")
        plt.savefig(out_png_path, dpi=300)
        print(f"Saving PDF plot to: {out_pdf_path}")
        plt.savefig(out_pdf_path)
        plt.close()

print("All plots generated successfully!")
