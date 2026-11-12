#!/usr/bin/env python3
"""
plot_purple_linear.py — Sequenza-style horizontal genome-wide plotter for HMF PURPLE outputs.

Renders a 3-track linear profile across chr1-chrY from PURPLE copy-number and purity outputs,
with optional overlays of raw AMBER BAF SNPs and binned COBALT depth ratios:
  1. Total & Allele-Specific Copy Number (Total, Major, and Minor allele CN segments)
  2. Read Depth Ratio (log2 tumor/normal GC-normalized COBALT ratio + fitted segments)
  3. B-Allele Frequency (AMBER tumorBAF scatter + PURPLE segment BAF / mirrored 1 - BAF)

Works on both complete and preliminary PURPLE runs (no somatic/SV VCFs required).
"""

import argparse
import csv
import gzip
import os
import sys
import tempfile
from pathlib import Path

# Ensure Matplotlib uses a writable config/cache directory on compute nodes
if not os.environ.get("MPLCONFIGDIR"):
    os.environ["MPLCONFIGDIR"] = os.path.join(tempfile.gettempdir(), f"mpl_cache_{os.getuid()}")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


CHROM_LENGTHS_38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257522, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

CHROM_LENGTHS_37 = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
}


def normalize_chrom(series: pd.Series) -> pd.Series:
    """Normalize chromosome names to canonical 'chr1'..'chr22','chrX','chrY'."""
    s = series.astype(str).str.strip()
    s = s.where(s.str.startswith("chr"), "chr" + s)
    return s


def safe_float(val, default=None):
    if val is None or val == "" or str(val).lower() in ("na", "nan", "none"):
        return default
    try:
        v = float(val)
        return v if np.isfinite(v) else default
    except (ValueError, TypeError):
        return default


def build_genome_coords(genome: str):
    """Return ordered chromosome lengths, cumulative offsets, centers, and total length."""
    chrom_lengths = CHROM_LENGTHS_37 if "37" in genome or "19" in genome else CHROM_LENGTHS_38
    offsets = {}
    centers = {}
    boundaries = []
    cur = 0
    for chrom, length in chrom_lengths.items():
        offsets[chrom] = cur
        centers[chrom] = cur + length / 2.0
        cur += length
        boundaries.append(cur)
    return chrom_lengths, offsets, centers, boundaries[:-1], cur


def _first_match(directory: Path, patterns: list[str]) -> Path | None:
    if not directory or not directory.is_dir():
        return None
    for pat in patterns:
        matches = sorted(
            p for p in directory.glob(pat)
            if p.is_file() and not p.name.startswith("._")
        )
        if matches:
            return matches[0]
    return None


def discover_inputs(args) -> dict:
    """Resolve PURPLE, AMBER, and COBALT input files from CLI flags or directory layout."""
    purple_dir = Path(args.purple_dir).resolve() if args.purple_dir else None
    sample_id = args.sample_id

    purity_tsv = Path(args.purity_tsv) if args.purity_tsv else None
    cnv_somatic = Path(args.cnv_somatic) if args.cnv_somatic else None
    segment_tsv = Path(args.segment_tsv) if args.segment_tsv else None
    qc_file = None

    if purple_dir:
        prefix = f"{sample_id}." if sample_id else "*."
        if not purity_tsv:
            purity_tsv = _first_match(purple_dir, [f"{prefix}purple.purity.tsv", "*.purple.purity.tsv"])
        if not cnv_somatic:
            cnv_somatic = _first_match(purple_dir, [f"{prefix}purple.cnv.somatic.tsv", "*.purple.cnv.somatic.tsv"])
        if not segment_tsv:
            segment_tsv = _first_match(purple_dir, [f"{prefix}purple.segment.tsv", "*.purple.segment.tsv"])
        qc_file = _first_match(purple_dir, [f"{prefix}purple.qc", "*.purple.qc"])

    if not sample_id:
        for p, suffix in [
            (purity_tsv, ".purple.purity.tsv"),
            (cnv_somatic, ".purple.cnv.somatic.tsv"),
            (segment_tsv, ".purple.segment.tsv"),
        ]:
            if p and p.name.endswith(suffix):
                sample_id = p.name[: -len(suffix)]
                break
        if not sample_id and purple_dir:
            sample_id = purple_dir.parent.name if purple_dir.name == "purple" else purple_dir.name

    # Resolve optional AMBER BAF
    amber_baf = Path(args.amber_baf) if args.amber_baf else None
    if not amber_baf:
        cand_amber_dirs = []
        if args.amber_dir:
            cand_amber_dirs.append(Path(args.amber_dir))
        if purple_dir:
            cand_amber_dirs.append(purple_dir.parent / "amber")
            # Also check standard sibling sample folder when purple lives in a prelim dir
            if sample_id:
                cand_amber_dirs.extend(sorted(purple_dir.parent.parent.parent.glob(f"*/{sample_id}/amber")))
        for ad in cand_amber_dirs:
            amber_baf = _first_match(ad, [f"{sample_id}.amber.baf.tsv.gz", "*.amber.baf.tsv.gz", "*.amber.baf.tsv"])
            if amber_baf:
                break

    # Resolve optional COBALT ratio
    cobalt_ratio = Path(args.cobalt_ratio) if args.cobalt_ratio else None
    if not cobalt_ratio:
        cand_cobalt_dirs = []
        if args.cobalt_dir:
            cand_cobalt_dirs.append(Path(args.cobalt_dir))
        if purple_dir:
            cand_cobalt_dirs.append(purple_dir.parent / "cobalt")
            if sample_id:
                cand_cobalt_dirs.extend(sorted(purple_dir.parent.parent.parent.glob(f"*/{sample_id}/cobalt")))
        for cd in cand_cobalt_dirs:
            cobalt_ratio = _first_match(
                cd,
                [
                    f"{sample_id}.cobalt.ratio.tsv.gz",
                    f"{sample_id}.cobalt.ratio.tsv",
                    "*.cobalt.ratio.tsv.gz",
                    "*.cobalt.ratio.tsv",
                ],
            )
            if cobalt_ratio:
                break

    return {
        "sample_id": sample_id or "Sample",
        "purple_dir": purple_dir,
        "purity_tsv": purity_tsv,
        "cnv_somatic": cnv_somatic,
        "segment_tsv": segment_tsv,
        "qc_file": qc_file,
        "amber_baf": amber_baf,
        "cobalt_ratio": cobalt_ratio,
    }


def parse_purity_info(purity_tsv: Path | None, qc_file: Path | None) -> dict:
    """Extract purity, ploidy, WGD, gender, and QC status from PURPLE summary files."""
    info = {
        "purity": None,
        "ploidy": None,
        "norm_factor": None,
        "wgd": None,
        "gender": "",
        "status": "",
    }
    if purity_tsv and purity_tsv.is_file():
        with open(purity_tsv, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                info["purity"] = safe_float(row.get("purity"))
                info["ploidy"] = safe_float(row.get("ploidy"))
                info["norm_factor"] = safe_float(row.get("normFactor"))
                wgd_raw = str(row.get("wholeGenomeDuplication", "")).strip().lower()
                if wgd_raw in ("true", "false"):
                    info["wgd"] = wgd_raw == "true"
                info["gender"] = (row.get("gender") or "").strip()
                info["status"] = (row.get("status") or "").strip()
                break

    if not info["status"] and qc_file and qc_file.is_file():
        with open(qc_file, "r", encoding="utf-8-sig") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0].strip() == "QCStatus":
                    info["status"] = parts[1].strip()
                    break
    return info


def load_purple_segments(cnv_somatic: Path | None, segment_tsv: Path | None, offsets: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load PURPLE CNV somatic segments and raw PURPLE segments (for observedTumorRatio)."""
    df_cnv = pd.DataFrame()
    if cnv_somatic and cnv_somatic.is_file():
        df_cnv = pd.read_csv(cnv_somatic, sep="\t")
    elif segment_tsv and segment_tsv.is_file():
        df_cnv = pd.read_csv(segment_tsv, sep="\t")

    if not df_cnv.empty:
        df_cnv["chrom_norm"] = normalize_chrom(df_cnv["chromosome"])
        df_cnv = df_cnv[df_cnv["chrom_norm"].isin(offsets)].copy()
        df_cnv["x0"] = df_cnv["chrom_norm"].map(offsets) + pd.to_numeric(df_cnv["start"], errors="coerce")
        df_cnv["x1"] = df_cnv["chrom_norm"].map(offsets) + pd.to_numeric(df_cnv["end"], errors="coerce")
        if "copyNumber" not in df_cnv.columns and "tumorCopyNumber" in df_cnv.columns:
            df_cnv["copyNumber"] = df_cnv["tumorCopyNumber"]
        for col in ["copyNumber", "minorAlleleCopyNumber", "majorAlleleCopyNumber", "baf", "observedBAF", "bafCount", "depthWindowCount"]:
            if col in df_cnv.columns:
                df_cnv[col] = pd.to_numeric(df_cnv[col], errors="coerce")
        if "majorAlleleCopyNumber" not in df_cnv.columns and "copyNumber" in df_cnv.columns and "minorAlleleCopyNumber" in df_cnv.columns:
            df_cnv["majorAlleleCopyNumber"] = df_cnv["copyNumber"] - df_cnv["minorAlleleCopyNumber"]
        df_cnv = df_cnv.dropna(subset=["x0", "x1"]).copy()

    df_seg = pd.DataFrame()
    if segment_tsv and segment_tsv.is_file():
        df_seg = pd.read_csv(segment_tsv, sep="\t")
        df_seg["chrom_norm"] = normalize_chrom(df_seg["chromosome"])
        df_seg = df_seg[df_seg["chrom_norm"].isin(offsets)].copy()
        df_seg["x0"] = df_seg["chrom_norm"].map(offsets) + pd.to_numeric(df_seg["start"], errors="coerce")
        df_seg["x1"] = df_seg["chrom_norm"].map(offsets) + pd.to_numeric(df_seg["end"], errors="coerce")
        for col in ["observedTumorRatio", "depthWindowCount", "tumorCopyNumber", "copyNumber"]:
            if col in df_seg.columns:
                df_seg[col] = pd.to_numeric(df_seg[col], errors="coerce")
        if "germlineStatus" in df_seg.columns:
            df_seg = df_seg[~df_seg["germlineStatus"].isin(["CENTROMETIC", "CENTROMERIC", "EXCLUDED", "NOISE"])].copy()
        df_seg = df_seg.dropna(subset=["x0", "x1"]).copy()

    return df_cnv, df_seg


def load_cobalt_binned(cobalt_path: Path | None, offsets: dict, bin_size: int) -> pd.DataFrame:
    """Read and median-bin COBALT tumorGCRatio values into genomic bins."""
    if not cobalt_path or not cobalt_path.is_file():
        return pd.DataFrame()

    # Peek header so we only read the required columns
    opener = gzip.open if str(cobalt_path).endswith(".gz") else open
    with opener(cobalt_path, "rt", encoding="utf-8-sig") as f:
        header = f.readline().rstrip("\n").split("\t")

    wanted = [c for c in ["chromosome", "position", "referenceReadDepth", "tumorGCRatio"] if c in header]
    if not {"chromosome", "position", "tumorGCRatio"}.issubset(wanted):
        return pd.DataFrame()

    df = pd.read_csv(cobalt_path, sep="\t", usecols=wanted, low_memory=False)
    df["chrom_norm"] = normalize_chrom(df["chromosome"])
    df = df[df["chrom_norm"].isin(offsets)].copy()
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["tumorGCRatio"] = pd.to_numeric(df["tumorGCRatio"], errors="coerce")
    mask = (df["tumorGCRatio"] > 0) & df["position"].notna()
    if "referenceReadDepth" in df.columns:
        df["referenceReadDepth"] = pd.to_numeric(df["referenceReadDepth"], errors="coerce")
        mask &= df["referenceReadDepth"] > 0
    df = df[mask].copy()
    if df.empty:
        return pd.DataFrame()

    df["bin"] = (df["position"] // bin_size) * bin_size + (bin_size // 2)
    binned = df.groupby(["chrom_norm", "bin"], as_index=False)["tumorGCRatio"].median()
    binned["x"] = binned["chrom_norm"].map(offsets) + binned["bin"]
    binned["log2_ratio"] = np.log2(np.clip(binned["tumorGCRatio"].values, 0.05, 8.0))
    return binned


def load_amber_baf(amber_path: Path | None, offsets: dict, max_points: int) -> pd.DataFrame:
    """Read and downsample AMBER tumorBAF SNP points across standard chromosomes."""
    if not amber_path or not amber_path.is_file():
        return pd.DataFrame()

    opener = gzip.open if str(amber_path).endswith(".gz") else open
    with opener(amber_path, "rt", encoding="utf-8-sig") as f:
        header = f.readline().rstrip("\n").split("\t")

    wanted = [c for c in ["chromosome", "position", "tumorBAF"] if c in header]
    if len(wanted) < 3:
        return pd.DataFrame()

    df = pd.read_csv(amber_path, sep="\t", usecols=wanted, low_memory=False)
    df["chrom_norm"] = normalize_chrom(df["chromosome"])
    df = df[df["chrom_norm"].isin(offsets)].copy()
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["tumorBAF"] = pd.to_numeric(df["tumorBAF"], errors="coerce")
    df = df.dropna(subset=["position", "tumorBAF"])
    df = df[(df["tumorBAF"] >= 0.0) & (df["tumorBAF"] <= 1.0)].copy()

    if len(df) > max_points:
        df = df.sample(n=max_points, random_state=42)
    df["x"] = df["chrom_norm"].map(offsets) + df["position"]
    return df


def _add_segment_collection(ax, x0, x1, y, color, lw=2.2, alpha=0.95, linestyle="-", zorder=4):
    """Efficiently render horizontal segments using LineCollection."""
    if len(x0) == 0:
        return
    segs = [[(float(a), float(v)), (float(b), float(v))] for a, b, v in zip(x0, x1, y)]
    lc = LineCollection(segs, colors=color, linewidths=lw, linestyles=linestyle, alpha=alpha, zorder=zorder)
    ax.add_collection(lc)


def _style_chrom_axis(ax, chrom_lengths, offsets, centers, boundaries, total_len, ylabel, show_xticks=False):
    ax.set_xlim(0, total_len)
    ax.set_ylabel(ylabel, fontsize=10.5, fontweight="bold", color="#0f172a")
    ax.tick_params(axis="y", labelsize=9, colors="#334155")
    for spine_name, spine in ax.spines.items():
        spine.set_color("#cbd5e1" if spine_name in ("top", "right") else "#64748b")

    for i, (chrom, length) in enumerate(chrom_lengths.items()):
        s = offsets[chrom]
        if i % 2 == 1:
            ax.axvspan(s, s + length, color="#f1f5f9", alpha=0.7, zorder=0)

    for b in boundaries:
        ax.axvline(b, color="#cbd5e1", linestyle="--", linewidth=0.7, zorder=1)

    tick_pos = [centers[c] for c in chrom_lengths]
    tick_lbl = [c.replace("chr", "") for c in chrom_lengths]
    ax.set_xticks(tick_pos)
    if show_xticks:
        ax.set_xticklabels(tick_lbl, fontsize=9, fontweight="bold", color="#1e293b")
        ax.set_xlabel("Chromosome", fontsize=11, fontweight="bold", color="#0f172a", labelpad=6)
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)


def compute_expected_log2_ratio(cn_vals: np.ndarray, purity: float | None, ploidy: float | None) -> np.ndarray:
    """Convert PURPLE segment copy number into expected log2 tumor/normal read-depth ratio."""
    p = purity if (purity is not None and 0.0 < purity <= 1.0) else 1.0
    pl = ploidy if (ploidy is not None and ploidy > 0.0) else 2.0
    sample_mean_ploidy = p * pl + (1.0 - p) * 2.0
    expected_ratio = ((1.0 - p) * 2.0 + p * np.maximum(0.0, cn_vals)) / max(sample_mean_ploidy, 1e-6)
    return np.log2(np.clip(expected_ratio, 0.05, 8.0))


def write_igv_seg(
    df_cnv: pd.DataFrame,
    df_seg: pd.DataFrame,
    sample_id: str,
    purity_info: dict,
    output_seg: Path,
    seg_mode: str = "ploidy",
    genome: str = "GRCh38",
):
    """Export PURPLE segments in standard 6-column IGV .seg format."""
    use_ratio = seg_mode == "ratio" and not df_seg.empty and "observedTumorRatio" in df_seg.columns
    src = df_seg.copy() if use_ratio else df_cnv.copy()

    if not use_ratio and "depthWindowCount" in src.columns:
        has_cov = src["depthWindowCount"].fillna(1) > 0
        if has_cov.any():
            src = src[has_cov].copy()

    src["start_i"] = pd.to_numeric(src["start"], errors="coerce")
    src["end_i"] = pd.to_numeric(src["end"], errors="coerce")
    src = src.dropna(subset=["start_i", "end_i"]).copy()
    src["start_i"] = src["start_i"].astype(int)
    src["end_i"] = src["end_i"].astype(int)
    src = src[src["end_i"] > src["start_i"]].copy()

    if src.empty:
        raise RuntimeError("No valid segments to export to IGV .seg")

    if "depthWindowCount" in src.columns:
        marks = pd.to_numeric(src["depthWindowCount"], errors="coerce").fillna(1).astype(int)
    elif "bafCount" in src.columns:
        marks = pd.to_numeric(src["bafCount"], errors="coerce").fillna(1).astype(int)
    else:
        marks = pd.Series(1, index=src.index, dtype=int)
    marks = marks.clip(lower=1)

    if use_ratio:
        vals = pd.to_numeric(src["observedTumorRatio"], errors="coerce").fillna(1.0).values
        seg_mean = np.maximum(-4.0, np.log2(np.maximum(vals, 1e-4)))
    elif seg_mode == "ratio":
        cn_vals = pd.to_numeric(src["copyNumber"], errors="coerce").fillna(2.0).values
        seg_mean = np.maximum(-4.0, compute_expected_log2_ratio(cn_vals, purity_info.get("purity"), purity_info.get("ploidy")))
    else:
        baseline = (
            purity_info["ploidy"]
            if (seg_mode == "ploidy" and purity_info.get("ploidy") and purity_info["ploidy"] > 0)
            else 2.0
        )
        cn_vals = pd.to_numeric(src["copyNumber"], errors="coerce").fillna(baseline).values
        seg_mean = np.maximum(-4.0, np.log2(np.maximum(cn_vals, 1e-4) / baseline))

    is_grch37 = "37" in genome or "19" in genome
    chrom_col = src["chrom_norm"].str.replace("^chr", "", regex=True) if is_grch37 else src["chrom_norm"]

    out_df = pd.DataFrame({
        "ID": sample_id,
        "chrom": chrom_col.values,
        "loc.start": src["start_i"].values,
        "loc.end": src["end_i"].values,
        "num.mark": marks.values,
        "seg.mean": np.round(seg_mean, 4),
    })

    output_seg.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_seg, sep="\t", index=False)
    print(f"[OK] Saved IGV SEG ({seg_mode}): {output_seg} ({len(out_df)} segments)")


def plot_purple_linear(inputs: dict, output_png: Path, output_pdf: Path | None, output_seg: Path | None, seg_mode: str, seg_only: bool, genome: str, max_cn_override: float | None, bin_size: int, max_baf_points: int, title_prefix: str | None):
    chrom_lengths, offsets, centers, boundaries, total_len = build_genome_coords(genome)
    purity_info = parse_purity_info(inputs["purity_tsv"], inputs["qc_file"])
    df_cnv, df_seg = load_purple_segments(inputs["cnv_somatic"], inputs["segment_tsv"], offsets)

    if df_cnv.empty:
        raise RuntimeError("No valid PURPLE segments found in *.purple.cnv.somatic.tsv or *.purple.segment.tsv")

    if output_seg:
        write_igv_seg(
            df_cnv=df_cnv,
            df_seg=df_seg,
            sample_id=inputs["sample_id"],
            purity_info=purity_info,
            output_seg=output_seg,
            seg_mode=seg_mode,
            genome=genome,
        )
        if seg_only:
            return

    df_cobalt = load_cobalt_binned(inputs["cobalt_ratio"], offsets, bin_size)
    df_amber = load_amber_baf(inputs["amber_baf"], offsets, max_baf_points)

    # Filter non-empty segments for plotting
    valid_cn = df_cnv.copy()
    if "depthWindowCount" in valid_cn.columns:
        has_cov = valid_cn["depthWindowCount"].fillna(1) > 0
        if has_cov.any():
            valid_cn = valid_cn[has_cov].copy()

    # Determine upper copy-number axis limit without squishing diploid regions on high-level focal amps
    cn_vals = valid_cn["copyNumber"].dropna().values if "copyNumber" in valid_cn.columns else np.array([2.0])
    raw_max_cn = float(np.nanmax(cn_vals)) if len(cn_vals) else 2.0
    if max_cn_override is not None:
        y_max_cn = max_cn_override
    else:
        p99 = float(np.nanpercentile(cn_vals, 99)) if len(cn_vals) else 4.0
        ploidy_val = purity_info["ploidy"] or 2.0
        y_max_cn = float(np.clip(np.ceil(max(6.0, p99 + 1.0, ploidy_val + 2.0)), 6.0, 12.0))

    fig, axes = plt.subplots(
        3, 1, figsize=(18, 9.2), dpi=200,
        gridspec_kw={"height_ratios": [1.35, 1.15, 1.15], "hspace": 0.14},
    )
    fig.patch.set_facecolor("#ffffff")
    ax_cn, ax_rat, ax_baf = axes

    # Header summary banner
    sample_label = title_prefix or inputs["sample_id"]
    meta_parts = []
    if purity_info["purity"] is not None:
        meta_parts.append(f"Purity: {purity_info['purity'] * 100:.1f}%")
    if purity_info["ploidy"] is not None:
        meta_parts.append(f"Ploidy: {purity_info['ploidy']:.2f}")
    if purity_info["wgd"] is not None:
        meta_parts.append(f"WGD: {'Yes' if purity_info['wgd'] else 'No'}")
    if purity_info["gender"]:
        meta_parts.append(f"Sex: {purity_info['gender']}")
    if purity_info["status"]:
        meta_parts.append(f"QC: {purity_info['status']}")
    meta_parts.append(f"Segments: {len(valid_cn)}")
    if raw_max_cn > y_max_cn:
        meta_parts.append(f"Max CN clipped at {y_max_cn:.0f} (peak {raw_max_cn:.1f})")

    fig.suptitle(f"{sample_label} — PURPLE Genome-Wide Linear Profile", fontsize=14.5, fontweight="bold", color="#0f172a", y=0.975)
    fig.text(0.5, 0.942, "   |   ".join(meta_parts), ha="center", va="top", fontsize=10, color="#334155", fontweight="medium")

    # ------------------------------------------------------------------
    # TRACK 1: Total & Allele-Specific Copy Number (Sequenza-style Major/Minor)
    # ------------------------------------------------------------------
    _style_chrom_axis(ax_cn, chrom_lengths, offsets, centers, boundaries, total_len, "Copy Number", show_xticks=False)
    ax_cn.set_ylim(-0.2, y_max_cn + 0.3)
    ax_cn.axhline(2.0, color="#94a3b8", linestyle=":", linewidth=1.0, zorder=1)
    if purity_info["ploidy"] is not None:
        ax_cn.axhline(purity_info["ploidy"], color="#0284c7", linestyle="--", linewidth=1.0, alpha=0.75, zorder=2)

    tot_cn = np.clip(valid_cn["copyNumber"].values, 0.0, y_max_cn)
    _add_segment_collection(ax_cn, valid_cn["x0"].values, valid_cn["x1"].values, tot_cn, color="#0f172a", lw=2.3, alpha=0.9, zorder=3)

    if "majorAlleleCopyNumber" in valid_cn.columns and "minorAlleleCopyNumber" in valid_cn.columns:
        sub = valid_cn.dropna(subset=["majorAlleleCopyNumber", "minorAlleleCopyNumber"])
        maj = np.clip(sub["majorAlleleCopyNumber"].values, 0.0, y_max_cn)
        min_a = np.clip(sub["minorAlleleCopyNumber"].values, 0.0, y_max_cn)
        # Apply tiny symmetric offset when major and minor coincide so both colors stay visible
        close_mask = np.abs(maj - min_a) < 0.08
        maj_plot = np.where(close_mask, maj + 0.045, maj)
        min_plot = np.where(close_mask, np.maximum(0.0, min_a - 0.045), min_a)
        _add_segment_collection(ax_cn, sub["x0"].values, sub["x1"].values, maj_plot, color="#dc2626", lw=2.0, alpha=0.9, zorder=4)
        _add_segment_collection(ax_cn, sub["x0"].values, sub["x1"].values, min_plot, color="#2563eb", lw=2.0, alpha=0.9, zorder=5)

    cn_legend = [
        Line2D([0], [0], color="#0f172a", lw=2.3, label="Total CN"),
        Line2D([0], [0], color="#dc2626", lw=2.0, label="Major Allele CN"),
        Line2D([0], [0], color="#2563eb", lw=2.0, label="Minor Allele CN"),
        Line2D([0], [0], color="#0284c7", lw=1.0, ls="--", label=f"Ploidy ({purity_info['ploidy']:.2f})" if purity_info["ploidy"] else "Ploidy"),
    ]
    ax_cn.legend(handles=cn_legend, loc="upper right", ncol=4, fontsize=8.5, frameon=True, facecolor="#ffffff", edgecolor="#cbd5e1")
    ax_cn.set_title("1. PURPLE Total & Allele-Specific Copy Number", loc="left", fontsize=10.5, fontweight="bold", color="#1e293b", pad=4)

    # ------------------------------------------------------------------
    # TRACK 2: Read Depth Ratio (COBALT log2 Tumor/Normal + Fitted/Observed Segments)
    # ------------------------------------------------------------------
    _style_chrom_axis(ax_rat, chrom_lengths, offsets, centers, boundaries, total_len, "Log2 Depth Ratio", show_xticks=False)
    ax_rat.set_ylim(-2.0, 2.0)
    ax_rat.set_yticks([-2.0, -1.0, 0.0, 1.0, 2.0])
    ax_rat.axhline(0.0, color="#64748b", linestyle="--", linewidth=0.9, alpha=0.8, zorder=1)

    rat_legend = []
    if not df_cobalt.empty:
        ax_rat.scatter(df_cobalt["x"].values, df_cobalt["log2_ratio"].values, s=4.5, c="#64748b", alpha=0.45, edgecolors="none", rasterized=True, zorder=2)
        rat_legend.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#64748b", markersize=5, label=f"COBALT {bin_size // 1000} kb median"))

    if not df_seg.empty and "observedTumorRatio" in df_seg.columns:
        seg_rat = df_seg[df_seg["observedTumorRatio"] > 0].copy()
        if not seg_rat.empty:
            y_obs = np.clip(np.log2(seg_rat["observedTumorRatio"].values), -2.0, 2.0)
            _add_segment_collection(ax_rat, seg_rat["x0"].values, seg_rat["x1"].values, y_obs, color="#0f172a", lw=1.8, alpha=0.85, zorder=3)
            rat_legend.append(Line2D([0], [0], color="#0f172a", lw=1.8, label="Observed Tumor Ratio"))

    expected_log2 = np.clip(compute_expected_log2_ratio(valid_cn["copyNumber"].values, purity_info["purity"], purity_info["ploidy"]), -2.0, 2.0)
    _add_segment_collection(ax_rat, valid_cn["x0"].values, valid_cn["x1"].values, expected_log2, color="#0d9488", lw=2.0, alpha=0.9, zorder=4)
    rat_legend.append(Line2D([0], [0], color="#0d9488", lw=2.0, label="PURPLE Expected Ratio"))
    ax_rat.legend(handles=rat_legend, loc="upper right", ncol=3, fontsize=8.5, frameon=True, facecolor="#ffffff", edgecolor="#cbd5e1")
    ax_rat.set_title("2. GC-Normalized Read Depth Ratio (Log2)", loc="left", fontsize=10.5, fontweight="bold", color="#1e293b", pad=4)

    # ------------------------------------------------------------------
    # TRACK 3: B-Allele Frequency (AMBER Het SNPs + PURPLE Segment BAF)
    # ------------------------------------------------------------------
    _style_chrom_axis(ax_baf, chrom_lengths, offsets, centers, boundaries, total_len, "Tumor BAF", show_xticks=True)
    ax_baf.set_ylim(-0.03, 1.03)
    ax_baf.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_baf.axhline(0.5, color="#64748b", linestyle="--", linewidth=0.9, alpha=0.8, zorder=1)

    baf_legend = []
    if not df_amber.empty:
        ax_baf.scatter(df_amber["x"].values, df_amber["tumorBAF"].values, s=2.0, c="#7c3aed", alpha=0.25, edgecolors="none", rasterized=True, zorder=2)
        baf_legend.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#7c3aed", markersize=5, label=f"AMBER SNPs (n={len(df_amber):,})"))

    baf_col = "baf" if "baf" in valid_cn.columns else ("observedBAF" if "observedBAF" in valid_cn.columns else None)
    if baf_col:
        sub_baf = valid_cn.dropna(subset=[baf_col]).copy()
        if "bafCount" in sub_baf.columns:
            sub_baf = sub_baf[sub_baf["bafCount"].fillna(0) > 0]
        sub_baf = sub_baf[(sub_baf[baf_col] >= 0.0) & (sub_baf[baf_col] <= 1.0)]
        if not sub_baf.empty:
            bvals = sub_baf[baf_col].values
            _add_segment_collection(ax_baf, sub_baf["x0"].values, sub_baf["x1"].values, bvals, color="#b45309", lw=2.1, alpha=0.95, zorder=4)
            _add_segment_collection(ax_baf, sub_baf["x0"].values, sub_baf["x1"].values, 1.0 - bvals, color="#b45309", lw=2.1, alpha=0.95, zorder=4)
            baf_legend.append(Line2D([0], [0], color="#b45309", lw=2.1, label="PURPLE Segment BAF / 1-BAF"))

    if baf_legend:
        ax_baf.legend(handles=baf_legend, loc="upper right", ncol=2, fontsize=8.5, frameon=True, facecolor="#ffffff", edgecolor="#cbd5e1")
    ax_baf.set_title("3. B-Allele Frequency (BAF)", loc="left", fontsize=10.5, fontweight="bold", color="#1e293b", pad=4)

    fig.subplots_adjust(left=0.055, right=0.985, top=0.895, bottom=0.07)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=200)
    print(f"[OK] Saved PNG: {output_png}")
    if output_pdf:
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_pdf)
        print(f"[OK] Saved PDF: {output_pdf}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot HMF PURPLE outputs as a horizontal Sequenza-style genome-wide profile."
    )
    parser.add_argument("--purple-dir", help="Path to sample purple/ directory")
    parser.add_argument("--sample-id", help="Sample ID override (auto-detected from PURPLE filenames)")
    parser.add_argument("--purity-tsv", help="Explicit path to *.purple.purity.tsv")
    parser.add_argument("--cnv-somatic", help="Explicit path to *.purple.cnv.somatic.tsv")
    parser.add_argument("--segment-tsv", help="Explicit path to *.purple.segment.tsv")
    parser.add_argument("--amber-dir", help="Path to sample amber/ directory (optional)")
    parser.add_argument("--amber-baf", help="Explicit path to *.amber.baf.tsv.gz (optional)")
    parser.add_argument("--cobalt-dir", help="Path to sample cobalt/ directory (optional)")
    parser.add_argument("--cobalt-ratio", help="Explicit path to *.cobalt.ratio.tsv[.gz] (optional)")
    parser.add_argument("-o", "--output", help="Output PNG path (default: <purple_dir>/<sample_id>.purple.linear.png)")
    parser.add_argument("--output-pdf", help="Optional output PDF path")
    parser.add_argument("--pdf", action="store_true", help="Also save a companion PDF alongside --output")
    parser.add_argument("--output-seg", help="Output IGV .seg path (default: <output_dir>/<sample_id>.purple.seg)")
    parser.add_argument("--no-seg", action="store_true", help="Skip generating the IGV .seg file")
    parser.add_argument("--seg-only", action="store_true", help="Only generate the IGV .seg file without plotting")
    parser.add_argument(
        "--seg-mode",
        choices=["ploidy", "diploid", "ratio"],
        default="ploidy",
        help="IGV seg.mean value: log2(copyNumber/ploidy) [default], log2(copyNumber/2), or log2(observedTumorRatio)",
    )
    parser.add_argument("--genome", default="GRCh38", help="Genome build: GRCh38 (default) or GRCh37")
    parser.add_argument("--max-cn", type=float, default=None, help="Upper Y-axis cap for copy number track")
    parser.add_argument("--bin-size", type=int, default=250000, help="COBALT binning window size in bp (default: 250000)")
    parser.add_argument("--max-baf-points", type=int, default=120000, help="Max AMBER SNP points to plot (default: 120000)")
    parser.add_argument("--title", help="Custom title prefix")
    args = parser.parse_args()

    if not args.purple_dir and not (args.cnv_somatic or args.segment_tsv):
        parser.error("Provide --purple-dir or --cnv-somatic / --segment-tsv")

    inputs = discover_inputs(args)
    if args.output:
        output_png = Path(args.output)
    elif inputs["purple_dir"]:
        output_png = inputs["purple_dir"] / f"{inputs['sample_id']}.purple.linear.png"
    else:
        output_png = Path(f"{inputs['sample_id']}.purple.linear.png")

    output_pdf = Path(args.output_pdf) if args.output_pdf else (output_png.with_suffix(".pdf") if args.pdf else None)
    if args.no_seg:
        output_seg = None
    elif args.output_seg:
        output_seg = Path(args.output_seg)
    else:
        output_seg = output_png.with_name(f"{inputs['sample_id']}.purple.seg")

    print(f"[INFO] Sample:        {inputs['sample_id']}")
    print(f"[INFO] Purity TSV:    {inputs['purity_tsv']}")
    print(f"[INFO] CNV Somatic:   {inputs['cnv_somatic']}")
    print(f"[INFO] Segment TSV:   {inputs['segment_tsv']}")
    print(f"[INFO] AMBER BAF:     {inputs['amber_baf'] or 'None (segment BAF only)'}")
    print(f"[INFO] COBALT Ratio:  {inputs['cobalt_ratio'] or 'None (segment ratio only)'}")

    plot_purple_linear(
        inputs=inputs,
        output_png=output_png,
        output_pdf=output_pdf,
        output_seg=output_seg,
        seg_mode=args.seg_mode,
        seg_only=args.seg_only,
        genome=args.genome,
        max_cn_override=args.max_cn,
        bin_size=args.bin_size,
        max_baf_points=args.max_baf_points,
        title_prefix=args.title,
    )


if __name__ == "__main__":
    main()
