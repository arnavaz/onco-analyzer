#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # Safe for headless cluster runs
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker
import argparse
import os
import sys
import urllib.request
import json
import gzip

def get_offline_gene_structure(gtf_path, gene_symbol):
    """Parse a gzipped GTF file to extract canonical transcript structure for a gene offline."""
    gene_info = None
    transcripts = {} # transcript_id -> list of exons
    transcript_info = {} # transcript_id -> metadata dict
    
    try:
        with gzip.open(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = 1 if parts[6] == '+' else -1
                attributes = parts[8]
                
                # Normalize chromosome to 'chr1', 'chrX', etc.
                chrom_str = str(chrom)
                if not chrom_str.lower().startswith('chr'):
                    chrom_str = 'chr' + chrom_str
                
                # Quick check if it contains the gene symbol
                if f'gene_name "{gene_symbol}"' not in attributes:
                    continue
                
                # Parse attributes
                attr_dict = {}
                for item in attributes.split(';'):
                    item = item.strip()
                    if not item:
                        continue
                    if ' ' in item:
                        k, v = item.split(' ', 1)
                        attr_dict[k] = v.replace('"', '')
                
                if attr_dict.get('gene_name') != gene_symbol:
                    continue
                
                if feature == 'gene':
                    gene_info = {
                        'gene_start': start,
                        'gene_end': end,
                        'strand': strand,
                        'chrom': chrom_str
                    }
                    
                elif feature == 'transcript':
                    tx_id = attr_dict.get('transcript_id')
                    if tx_id:
                        is_canonical = ('tag' in attr_dict and attr_dict['tag'] == 'canonical') or ('ensembl_canonical' in attr_dict) or ('tag "canonical"' in attributes)
                        is_mane = ('tag' in attr_dict and attr_dict['tag'] == 'MANE_Select') or ('tag "MANE_Select"' in attributes)
                        
                        transcript_info[tx_id] = {
                            'gene_start': start,
                            'gene_end': end,
                            'strand': strand,
                            'chrom': chrom_str,
                            'is_canonical': is_canonical,
                            'is_mane': is_mane,
                            'exon_count': 0
                        }
                        
                elif feature == 'exon':
                    tx_id = attr_dict.get('transcript_id')
                    if tx_id:
                        if tx_id not in transcripts:
                            transcripts[tx_id] = []
                        transcripts[tx_id].append({
                            'start': start,
                            'end': end,
                            'id': attr_dict.get('exon_id', f"{tx_id}_exon")
                        })
    except Exception as e:
        print(f"Warning: Failed to parse GTF file offline for {gene_symbol}: {e}", file=sys.stderr)
        return None
        
    if not transcript_info:
        return None
        
    # Selection priority: MANE Select > Canonical > transcript with the most exons
    best_tx = None
    # 1. MANE Select
    for tx_id, info in transcript_info.items():
        if info['is_mane']:
            best_tx = tx_id
            break
    # 2. Canonical
    if not best_tx:
        for tx_id, info in transcript_info.items():
            if info['is_canonical']:
                best_tx = tx_id
                break
    # 3. Fallback to transcript with most exons
    if not best_tx:
        for tx_id, exons in transcripts.items():
            if tx_id in transcript_info:
                transcript_info[tx_id]['exon_count'] = len(exons)
        if transcript_info:
            best_tx = max(transcript_info.keys(), key=lambda k: transcript_info[k].get('exon_count', 0))
            
    if best_tx and best_tx in transcripts:
        exons = sorted(transcripts[best_tx], key=lambda x: x['start'])
        info = transcript_info[best_tx]
        return {
            'gene_start': gene_info['gene_start'] if gene_info else info['gene_start'],
            'gene_end': gene_info['gene_end'] if gene_info else info['gene_end'],
            'strand': info['strand'],
            'chrom': info['chrom'],
            'transcript_id': best_tx,
            'exons': exons
        }
        
    return None

def get_ensembl_canonical_exons(gene_symbol):
    """Fetch gene coordinates, strand, and exon structure for the canonical transcript from Ensembl."""
    url = f"https://rest.ensembl.org/lookup/symbol/human/{gene_symbol}?expand=1;mane=1"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    try:
        # Ensembl REST API rate-limiting compliance (up to 15 QPS)
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
            
            # Find the canonical transcript
            transcripts = data.get("Transcript", [])
            canonical = None
            for t in transcripts:
                if t.get("is_canonical") == 1:
                    canonical = t
                    break
            
            # Fallback to first transcript if no canonical flag is explicit
            if not canonical and transcripts:
                canonical = transcripts[0]
                
            if canonical:
                exons = []
                for ex in canonical.get("Exon", []):
                    exons.append({
                        'start': int(ex['start']),
                        'end': int(ex['end']),
                        'id': ex['id']
                    })
                # Sort exons by genomic start
                exons = sorted(exons, key=lambda x: x['start'])
                return {
                    'gene_start': int(data['start']),
                    'gene_end': int(data['end']),
                    'strand': int(data['strand']),
                    'chrom': str(data['seq_region_name']),
                    'transcript_id': canonical['id'],
                    'exons': exons
                }
    except Exception as e:
        print(f"Warning: Failed to fetch Ensembl exons for {gene_symbol}: {e}", file=sys.stderr)
    return None

def main():
    parser = argparse.ArgumentParser(description='Generate a gene-focused BAF plot from MAF')
    parser.add_argument('--maf', required=True, help='Path to MAF file')
    parser.add_argument('--genes-file', required=True, help='Path to MMR_drivers.tsv or similar gene file')
    parser.add_argument('--output', required=True, help='Path to output base path/folder')
    parser.add_argument('--title', help='Title of the plot')
    parser.add_argument('--purple-segments', help='Path to purple.segment.tsv file to draw clean segment lines')
    parser.add_argument('--min-depth', type=int, default=100, help='Minimum depth (t_depth) for SNPs')
    parser.add_argument('--vaf-min', type=float, default=0.15, help='Minimum VAF to filter out somatic/mapping noise')
    parser.add_argument('--vaf-max', type=float, default=0.85, help='Maximum VAF to filter out homozygous germline SNPs')
    parser.add_argument('--gtf', help='Path to gzipped Ensembl GTF file (.gtf.gz) for offline gene track annotations')
    
    args = parser.parse_args()
    
    # Load panel genes
    if not os.path.exists(args.genes_file):
        print(f"Error: Genes file not found: {args.genes_file}", file=sys.stderr)
        sys.exit(1)
        
    drivers_df = pd.read_csv(args.genes_file, sep='\t')
    panel_genes = sorted(list(set(drivers_df.iloc[:, 0].dropna().unique())))
    
    # Load MAF
    if not os.path.exists(args.maf):
        print(f"Error: MAF file not found: {args.maf}", file=sys.stderr)
        sys.exit(1)
        
    df = pd.read_csv(args.maf, sep='\t', comment='#', low_memory=False)
    
    # Check required columns
    required_cols = ['Hugo_Symbol', 'dbSNP_RS', 't_ref_count', 't_alt_count', 'Chromosome', 'Start_Position']
    for col in required_cols:
        if col not in df.columns:
            print(f"Error: Missing required column in MAF: {col}", file=sys.stderr)
            sys.exit(1)
            
    # Determine tumor sample ID from the MAF
    tumor_id = "UNKNOWN_TUMOR"
    if 'Tumor_Sample_Barcode' in df.columns and not df.empty:
        tumor_id = str(df['Tumor_Sample_Barcode'].iloc[0])
        
    # Filter for dbSNP sites (starts with 'rs')
    df = df[df['dbSNP_RS'].astype(str).str.startswith('rs', na=False)]
    
    # Calculate BAF and VAF
    df['total_count'] = df['t_ref_count'] + df['t_alt_count']
    
    # Filter out low depth to decrease noise (default 100x for panel)
    df = df[df['total_count'] >= args.min_depth]
    
    # Calculate BAF
    df['BAF'] = df['t_alt_count'] / df['total_count']
    
    # Filter out homozygous SNPs and low VAF somatic/mapping noise
    df = df[(df['BAF'] >= args.vaf_min) & (df['BAF'] <= args.vaf_max)]

    # Load PURPLE segments if provided
    segments_df = None
    if args.purple_segments and os.path.exists(args.purple_segments):
        try:
            segments_df = pd.read_csv(args.purple_segments, sep='\t')
            # Normalize segment chromosomes
            segments_df['norm_chrom'] = segments_df['chromosome'].astype(str).str.upper().str.replace('CHR', '')
        except Exception as e:
            print(f"Warning: Failed to load purple segments: {e}", file=sys.stderr)
            
    # Output directory resolution
    if args.output.endswith('.png') or args.output.endswith('.pdf'):
        output_dir = os.path.dirname(args.output)
    else:
        output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate one plot per gene
    for gene in panel_genes:
        print(f"Processing BAF plot for gene: {gene}")
        
        gene_df = df[df['Hugo_Symbol'] == gene]
        
        # Fetch exons and gene structure
        if args.gtf and os.path.exists(args.gtf):
            gene_info = get_offline_gene_structure(args.gtf, gene)
        else:
            gene_info = get_ensembl_canonical_exons(gene)
        
        # If no variants and no ensembl info, skip this gene
        if gene_df.empty and not gene_info:
            print(f"  No SNPs and no Ensembl info for {gene}. Skipping plot.")
            continue
            
        # Determine x-limits for plotting
        if gene_info:
            gene_start = gene_info['gene_start']
            gene_end = gene_info['gene_end']
            chrom = gene_info['chrom']
            x_min = gene_start - 5000
            x_max = gene_end + 5000
        else:
            # Fallback if Ensembl query failed
            chrom = str(gene_df['Chromosome'].iloc[0]) if not gene_df.empty else "UNK"
            v_coords = gene_df['Start_Position'].values
            if len(v_coords) > 0:
                gene_start = int(min(v_coords))
                gene_end = int(max(v_coords))
            else:
                gene_start, gene_end = 0, 1000
            x_min = gene_start - 5000
            x_max = gene_end + 5000
            
        # Setup figure with shared x-axis
        fig, (ax_plot, ax_track) = plt.subplots(
            2, 1, 
            figsize=(14, 8), 
            sharex=True, 
            gridspec_kw={'height_ratios': [5, 1.2]}
        )
        
        # 1. BAF Plot (Top Panel)
        ax_plot.axhline(y=0.5, color='#94a3b8', linestyle='--', linewidth=1.2, alpha=0.8, label='Diploid Heterozygous (0.5)')
        ax_plot.axhspan(0, 0.2, color='#fca5a5', alpha=0.15, label='LOH / Deletion Bounds')
        ax_plot.axhspan(0.8, 1.0, color='#fca5a5', alpha=0.15)
        
        # Plot individual variants
        if not gene_df.empty:
            ax_plot.scatter(
                gene_df['Start_Position'], gene_df['BAF'],
                c='#3b82f6', alpha=0.8, edgecolors='#1d4ed8',
                linewidths=0.5, s=60, zorder=3,
                label=f'SNP BAF (n={len(gene_df)})'
            )
            
        # Plot segment lines from PURPLE if available
        plotted_seg_legend = False
        if segments_df is not None:
            c = chrom.upper().replace('CHR', '')
            # Find segments overlapping with our x range on the correct chromosome
            overlapping_segs = segments_df[
                (segments_df['norm_chrom'] == c) &
                (segments_df['start'] <= x_max) &
                (segments_df['end'] >= x_min) &
                (segments_df['bafCount'] > 0)
            ]
            
            for idx, seg in overlapping_segs.iterrows():
                seg_baf = seg.get('tumorBAF')
                if pd.notna(seg_baf):
                    # Clamp segment coordinates to x range for visual tidiness
                    seg_start = max(int(seg['start']), x_min)
                    seg_end = min(int(seg['end']), x_max)
                    
                    # Major allele segment line
                    ax_plot.hlines(
                        y=seg_baf, xmin=seg_start, xmax=seg_end,
                        colors='#10b981', linewidth=4.0, zorder=4,
                        label='PURPLE Segment BAF' if not plotted_seg_legend else ""
                    )
                    # Minor allele segment line
                    ax_plot.hlines(
                        y=1 - seg_baf, xmin=seg_start, xmax=seg_end,
                        colors='#10b981', linewidth=4.0, zorder=4
                    )
                    plotted_seg_legend = True
                    
        ax_plot.set_ylim(-0.02, 1.02)
        ax_plot.set_ylabel('Tumor B-Allele Frequency (BAF)', fontsize=12, fontweight='bold')
        ax_plot.set_yticks(np.arange(0, 1.1, 0.1))
        ax_plot.tick_params(axis='y', labelsize=10)
        ax_plot.grid(axis='y', linestyle=':', color='#cbd5e1', alpha=0.7)
        ax_plot.legend(loc='lower left', frameon=True, facecolor='white', edgecolor='none', fontsize=9)
        
        # Set top plot title
        title_text = args.title or f"{tumor_id} - {gene} BAF Profile"
        if not gene_df.empty:
            title_text += f" (Min Depth {args.min_depth}X)"
        ax_plot.set_title(title_text, fontsize=14, fontweight='bold', pad=10)
        
        # Clean spines for top plot
        ax_plot.spines['top'].set_visible(False)
        ax_plot.spines['right'].set_visible(False)
        ax_plot.spines['left'].set_color('#64748b')
        ax_plot.spines['bottom'].set_color('#64748b')
        
        # 2. Gene Body Track (Bottom Panel)
        if gene_info:
            exons = gene_info['exons']
            # Draw intron line
            ax_track.plot([gene_start, gene_end], [0, 0], color='#64748b', linewidth=2.0, zorder=1)
            
            # Draw transcription direction arrows on introns
            strand = gene_info['strand']
            arrow_step = (gene_end - gene_start) / 10
            if arrow_step > 0:
                arrow_positions = np.arange(gene_start + arrow_step, gene_end, arrow_step)
                for pos in arrow_positions:
                    direction = '>' if strand == 1 else '<'
                    ax_track.text(
                        pos, 0, direction,
                        ha='center', va='center', color='#64748b',
                        fontsize=9, fontweight='bold', zorder=2
                    )
            
            # Draw exons
            for idx, ex in enumerate(exons, 1):
                estart = ex['start']
                eend = ex['end']
                width = eend - estart + 1
                
                rect = Rectangle(
                    (estart, -0.4), width, 0.8,
                    facecolor='#3b82f6', edgecolor='#1d4ed8',
                    linewidth=1.0, alpha=0.85, zorder=3
                )
                ax_track.add_patch(rect)
                
                # Draw exon number if it fits
                if width > (gene_end - gene_start) * 0.02:
                    mid_x = estart + width / 2
                    ax_track.text(
                        mid_x, 0, str(idx),
                        ha='center', va='center', color='white',
                        fontsize=8, fontweight='bold', zorder=4
                    )
            ax_track.text(
                gene_start, 0.6, f"5' ENST: {gene_info['transcript_id']}" if strand == 1 else "3'",
                fontsize=8, color='#64748b', ha='left', va='bottom'
            )
            ax_track.text(
                gene_end, 0.6, "3'" if strand == 1 else f"5' ENST: {gene_info['transcript_id']}",
                fontsize=8, color='#64748b', ha='right', va='bottom'
            )
        else:
            # Fallback message
            ax_track.text(
                (x_min + x_max)/2, 0, "Gene Body Track Unavailable (Ensembl query offline/failed)",
                ha='center', va='center', color='gray', fontsize=10, style='italic'
            )
            ax_track.plot([x_min, x_max], [0, 0], color='#cbd5e1', linestyle='--', zorder=1)
            
        ax_track.set_ylim(-1.0, 1.0)
        ax_track.set_yticks([]) # Hide y ticks
        
        # Format X axis ticks as Megabases
        @ticker.FuncFormatter
        def mb_formatter(x, pos):
            return f"{x / 1e6:.3f} Mb"
            
        ax_track.xaxis.set_major_formatter(mb_formatter)
        ax_track.set_xlabel(f"Genomic Coordinates on Chromosome {chrom} (GRCh38)", fontsize=11, fontweight='bold')
        ax_track.tick_params(axis='x', labelsize=10)
        
        # Clean spines for gene track
        ax_track.spines['top'].set_visible(False)
        ax_track.spines['right'].set_visible(False)
        ax_track.spines['left'].set_visible(False)
        ax_track.spines['bottom'].set_color('#64748b')
        
        # Apply strict x range limits to keep plots clean and aligned
        plt.xlim(x_min, x_max)
        
        plt.tight_layout()
        
        # Save plot
        gene_output_path = os.path.join(output_dir, f"{tumor_id}_{gene}_baf.png")
        plt.savefig(gene_output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved plot to {gene_output_path}")

if __name__ == '__main__':
    main()
