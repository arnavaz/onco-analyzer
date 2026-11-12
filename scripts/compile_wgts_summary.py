#!/usr/bin/env python3
"""
compile_wgts_summary.py
=======================
Generalized cohort-level compilation and visualization tool for nf-core/oncoanalyser WGTS studies.
Supports WGS (and upcoming RNA-Seq) data, featuring:
  1. Cohort Landing Page with Overview KPIs & Samples/Patients Summary
  2. Quality Control (QC) Plots: Purity vs Ploidy, Sequencing Depth, TMB, MSI, ichorCNA
  3. Genomic Summary: Interactive Oncoprint of Top Mutated Genes & User Genes of Interest
  4. Minimal Residual Disease (MRD) Tab: Longitudinal plasma ctDNA trajectories with log-scale y-axis
     and greyed-out assay Limit of Detection (LOD) region

Generates:
  - cohort_wgts_summary.tsv          - Master patient/tumor metrics table
  - cohort_wgts_samples.tsv          - Tidy sample-level table (Tumor, Normal, Plasma)
  - cohort_wgts_oncoprint.tsv        - Gene x Sample alteration matrix
  - cohort_wgts_mrd_summary.tsv      - Detailed plasma ctDNA / WISP metrics (if available)
  - cohort_wgts_patient_matrix.tsv   - Longitudinal MRD matrix per patient
  - cohort_wgts_summary.png / .pdf   - Multi-panel publication figure
  - patient_plots/                   - Individual per-patient log-scale MRD trajectory plots
  - cohort_wgts_dashboard.html       - Standalone interactive HTML dashboard

Usage:
  python3 compile_wgts_summary.py --results-dir /path/to/Oncoanalyser/output
  python3 compile_wgts_summary.py --results-dir ... --cohort-name BFIERCE --genes-of-interest BAP1,NF2,SETD2,CDKN2A,TP53
"""

import os
import sys
import csv
import json
import gzip
import math
import glob
import re
import argparse
import shutil
import base64
from pathlib import Path
from collections import defaultdict
from datetime import datetime, date
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

# Set default writable MPLCONFIGDIR for HPC clusters and multiprocessing
if "MPLCONFIGDIR" not in os.environ:
    _user = os.environ.get("USER", "cluster_user")
    os.environ["MPLCONFIGDIR"] = f"/tmp/{_user}_matplotlib"

# Headless matplotlib backend for HPC clusters
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.collections import LineCollection
    from matplotlib.gridspec import GridSpec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import pandas as pd
    import numpy as np
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


# Human chromosome lengths (GRCh38)
CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895
}

def get_chrom_offsets():
    offsets = {}
    cum = 0
    for chrom, length in CHROM_LENGTHS.items():
        offsets[chrom] = cum
        offsets[chrom.replace("chr", "")] = cum
        cum += length
    return offsets, cum


DEFAULT_GENES_OF_INTEREST = [
    "TP53", "KRAS", "PIK3CA", "ARID1A", "NF2", "BAP1", "CDKN2A", "SETD2",
    "SMARCA4", "PBRM1", "PTEN", "BRCA1", "BRCA2", "EGFR", "BRAF", "ERBB2",
    "MYC", "CCND1", "RB1", "APC", "LATS1", "LATS2", "NF1", "CDK12"
]

# Known Pharmacogenomic (PGx) genes tracked by PEACH / PURPLE that are not oncogenic cancer drivers
PHARMACOGENOMIC_GENES = {
    "CYP1A2", "CYP2B6", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "CYP3A5",
    "DPYD", "UGT1A1", "UGT1A10", "TPMT", "NUDT15", "SLCO1B1", "VKORC1"
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile and visualize WGTS genomic and MRD results across an Oncoanalyser cohort."
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Base Oncoanalyser output directory containing patient folders (e.g. .../Oncoanalyser/).",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to save summaries and plots. Defaults to [results-dir]/wgts_summary/.",
    )
    parser.add_argument(
        "--cohort-name",
        help="Cohort display name (e.g. BFIERCE, TRACE_OC). Inferred from results directory if omitted.",
    )
    parser.add_argument(
        "--genes-of-interest",
        help="Comma-separated list or path to a file of gene symbols to prioritize in the Oncoprint.",
    )
    parser.add_argument(
        "--sample-sheets-dir",
        help="Path to sample_sheets directory to look up timepoint / sample aliases.",
    )
    parser.add_argument(
        "--metadata-tsv",
        help="Optional mapping/config TSV containing sample annotations or timepoints.",
    )
    parser.add_argument(
        "--timepoints-tsv",
        help="Optional mapping TSV/CSV containing sample timepoint labels, collection dates, or elapsed days.",
    )
    parser.add_argument(
        "--format",
        nargs="+",
        default=["png", "pdf"],
        help="Output image formats for static plots (default: png pdf).",
    )
    parser.add_argument(
        "--purity-overrides",
        help="Optional path to purity_overrides.tsv generated by curate_purity_fits.py to override purity and ploidy.",
    )
    parser.add_argument(
        "--skip-plasma-plots",
        action="store_true",
        help="Skip per-plasma multi-modal WISP decomposition plots to speed up summary generation.",
    )
    parser.add_argument(
        "--skip-fingerprint",
        "--skip-fingerprints",
        action="store_true",
        help="Skip Amber BAF sample fingerprinting and identity evaluation.",
    )
    parser.add_argument(
        "--fingerprint-stride",
        type=int,
        default=20,
        help="Locus sampling stride for Amber BAF fingerprinting (default: 20 for WGS, sampling ~10k loci across all autosomes; 1 disables downsampling).",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Number of parallel worker processes/threads for BAF fingerprinting and static figure generation (default: 8).",
    )
    parser.add_argument(
        "--fingerprint-threads",
        type=int,
        default=None,
        help="Legacy alias for --threads.",
    )
    parser.add_argument(
        "--embed-all-plots",
        action="store_true",
        help="Embed all high-resolution per-plasma multi-modal WISP decomposition plots as base64 in the interactive HTML dashboard (warning: increases HTML file size significantly for large cohorts).",
    )
    return parser.parse_args()


# ==============================================================================
# 1. HELPERS & SAFE PARSERS
# ==============================================================================

def safe_float(val, default=0.0):
    if val is None or val == "" or str(val).strip().upper() in {"NA", "NAN", "NULL"}:
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=0):
    if val is None or val == "" or str(val).strip().upper() in {"NA", "NAN", "NULL"}:
        return default
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return default


def safe_bool(val):
    if val is None:
        return False
    return str(val).strip().upper() in {"TRUE", "T", "1", "YES", "Y"}


def fmt_val(val, fmt=".2f", default="NA", mult=1.0):
    """Safely format a numerical or string value with format specifiers, handling None, NA, and invalid types cleanly."""
    if val is None or val == "" or str(val).strip().upper() in {"NA", "NAN", "NULL", "NONE"}:
        return default
    try:
        v = float(val) * mult
        return f"{v:{fmt}}"
    except Exception:
        return default


def parse_genes_arg(genes_arg):
    """Parse comma-separated gene list or text file into a set of uppercase gene symbols."""
    genes = []
    if not genes_arg:
        return DEFAULT_GENES_OF_INTEREST

    if os.path.isfile(genes_arg):
        try:
            with open(genes_arg, "r", encoding="utf-8") as f:
                for line in f:
                    for g in re.split(r"[,\t\s]+", line.strip()):
                        g_clean = g.strip().upper()
                        if g_clean:
                            genes.append(g_clean)
        except Exception as e:
            print(f"Warning: Failed to read genes-of-interest file {genes_arg}: {e}", file=sys.stderr)
    else:
        for g in re.split(r"[,\t\s]+", str(genes_arg).strip()):
            g_clean = g.strip().upper()
            if g_clean:
                genes.append(g_clean)

    # Deduplicate preserving order
    unique_genes = list(dict.fromkeys(genes))
    return unique_genes if unique_genes else DEFAULT_GENES_OF_INTEREST


# ==============================================================================
# 2. DISCOVERY & PARSING
# ==============================================================================

def find_candidate_patient_dirs(results_dir):
    """
    Search shallowly for patient directories under results_dir,
    strictly excluding massive work/ or .nextflow/ directories.
    """
    results_path = Path(results_dir).resolve()
    EXCLUDE_DIRS = {
        "work", ".nextflow", "logs", "work_dir", "scratch", ".git",
        "pipeline_info", "sample_sheets", "all_fastq_links", "realigned_fastqs",
        "wisp_summary", "wgts_summary", "cohort_summary", "metadata", "qc",
        "reports", "configs", "docs", "scripts", "tmp", "temp", "figures", "plots"
    }

    try:
        candidates = [
            d for d in results_path.iterdir()
            if d.is_dir() and d.name.lower() not in EXCLUDE_DIRS and not d.name.startswith(".")
        ]
    except Exception as e:
        print(f"Error accessing {results_dir}: {e}", file=sys.stderr)
        return []

    return sorted(candidates)


def parse_ichor_params_file(params_path):
    """Safely parse Tumor Fraction from an ichorCNA *params.txt file."""
    if not params_path or not os.path.isfile(params_path):
        return None
    try:
        with open(params_path, "r", encoding="utf-8-sig") as f:
            lines = [l.strip() for l in f if l.strip()]
        if not lines:
            return None
        # Check if line 0 is TSV header: Sample \t Tumor Fraction \t Ploidy
        parts0 = lines[0].split("\t")
        if len(lines) >= 2 and len(parts0) >= 2 and parts0[0].lower() == "sample":
            parts1 = lines[1].split("\t")
            if len(parts1) >= 2:
                v = safe_float(parts1[1])
                if v is not None:
                    return v
        # Alternatively search for 'Tumor Fraction:' label
        for line in lines:
            if line.lower().startswith("tumor fraction:"):
                val = line.split(":", 1)[-1].strip()
                v = safe_float(val)
                if v is not None:
                    return v
    except Exception:
        pass
    return None


def discover_cohort_files(results_dir):
    """
    Recursively and safely discover all relevant WGTS and WISP files across results_dir,
    strictly pruning huge scratch trees (work/, .nextflow/, logs/, etc.) in-place.
    Returns:
      - purple_files: list of Path to *.purple.purity.tsv
      - wisp_files: list of Path to *.wisp.summary.tsv or *.wisp.tsv
      - orange_map: dict mapping tumor_id/stem -> Path to *.orange.json
      - ichor_map: dict mapping sample_id/stem -> float tumor fraction
    """
    results_path = Path(results_dir).resolve()
    EXCLUDE_DIRS = {
        "work", ".nextflow", "logs", "work_dir", "scratch", ".git",
        "pipeline_info", "sample_sheets", "all_fastq_links", "realigned_fastqs",
        "wisp_summary", "wgts_summary", "cohort_summary", "metadata", "qc",
        "reports", "configs", "docs", "scripts", "tmp", "temp", "figures", "plots"
    }

    purple_files = []
    wisp_files = []
    orange_map = {}
    ichor_map = {}
    ichor_tumor_map = {}
    ichor_plasma_map = {}
    ichor_short_map = {}
    ichor_standard_map = {}
    ichor_prio = {}

    ichor_pdf_map = {}
    ichor_tumor_pdf_map = {}
    ichor_plasma_pdf_map = {}
    ichor_short_pdf_map = {}
    ichor_standard_pdf_map = {}

    for root, dirs, files in os.walk(results_path, followlinks=True):
        # In-place directory pruning to avoid traversing deep/massive scratch trees
        dirs[:] = [
            d for d in dirs
            if d.lower() not in EXCLUDE_DIRS and not d.startswith(".")
        ]

        root_path = Path(root)
        for f in files:
            if f.startswith("."):
                continue
            f_lower = f.lower()
            if f.endswith(".purple.purity.tsv"):
                purple_files.append(root_path / f)
            elif f.endswith(".wisp.summary.tsv") or (f.endswith(".wisp.tsv") and not any(f.endswith(x) for x in [".snv.wisp.tsv", ".cnv.wisp.tsv", ".loh.wisp.tsv"])):
                wisp_files.append(root_path / f)
            elif f.endswith(".orange.json"):
                stem = f.replace(".orange.json", "")
                try:
                    rel_p = root_path.relative_to(results_path).parts
                    if rel_p:
                        p_token = rel_p[0]
                        orange_map[(p_token, stem)] = root_path / f
                except Exception:
                    pass
                if stem not in orange_map:
                    orange_map[stem] = root_path / f
                elif orange_map[stem] is not None and orange_map[stem] != root_path / f:
                    orange_map[stem] = None
            elif f_lower.endswith(".pdf") and ("genomewide" in f_lower or "ichor" in str(root_path).lower()):
                file_stem = f[:-4]
                # Strip all genomeWide and short/standard variants to find base sample_id
                for sfx in [
                    "_genomeWide_all_sols", "_genomeWide", ".genomeWide", "_all_sols",
                    "_short_genomeWide", "_standard_genomeWide",
                    ".short_genomeWide", ".standard_genomeWide",
                    "_short.genomeWide", "_standard.genomeWide",
                ]:
                    if file_stem.endswith(sfx):
                        file_stem = file_stem[:-len(sfx)]
                        break
                for sub_sfx in ["_short", "_standard", ".short", ".standard"]:
                    if file_stem.endswith(sub_sfx):
                        file_stem = file_stem[:-len(sub_sfx)]
                        break

                sids = [file_stem]
                # Also include directory names if they look like sample identifiers
                for d_cand in [root_path.name, root_path.parent.name]:
                    if d_cand.lower() not in {"short", "standard", "ichorcna", "ichorcna_plasma_only", "ichorcna_tumor_only", "plasma", "tumor", "wisp", "amber"}:
                        if d_cand not in sids:
                            sids.append(d_cand)

                p_tokens = []
                try:
                    rel_p = root_path.relative_to(results_path).parts
                    if rel_p:
                        base_pt = rel_p[0]
                        p_tokens.append(base_pt)
                        if "_" in base_pt: p_tokens.append(base_pt.replace("_", "-"))
                        if "-" in base_pt: p_tokens.append(base_pt.replace("-", "_"))
                        m_pad = re.search(r'^(.*?)([-_])0*(\d+)$', base_pt)
                        if m_pad:
                            pref, sep, num = m_pad.groups()
                            for s_sep in ["-", "_"]:
                                p_tokens.append(f"{pref}{s_sep}{num}")
                                p_tokens.append(f"{pref}{s_sep}{int(num):04d}")
                                p_tokens.append(f"{pref}{s_sep}{int(num):02d}")
                except Exception:
                    pass

                prio = 2 if "_all_sols" not in f_lower else 1
                r_str = str(root_path).lower()
                if "tumor_only" in r_str or "plasma_only" in r_str:
                    prio += 2

                is_short = ("short" in f_lower or "short" in r_str or "size_select" in f_lower or "size_select" in r_str or "90_150" in f_lower or "90_150" in r_str)
                is_plasma = ("plasma" in r_str or "mrd" in r_str or "cfdna" in r_str)
                is_tumor = ("tumor" in r_str or "t_vs_n" in r_str or "tumour" in r_str)
                if is_short and is_plasma:
                    prio += 4

                for s in sids:
                    sample_variants = [s, s.replace("-", "_"), s.replace("_", "-")]
                    all_keys = list(sample_variants)
                    for pt in p_tokens:
                        clean_s = extract_clean_sample_name(s, pt)
                        if clean_s:
                            sample_variants.extend([clean_s, clean_s.replace("-", "_"), clean_s.replace("_", "-")])
                        for sv in sample_variants:
                            all_keys.extend([
                                (pt, sv),
                                f"{pt}_{sv}",
                                f"{pt}-{sv}",
                            ])
                    for k in all_keys:
                        curr_prio = ichor_prio.get(f"pdf_{k}", 0)
                        if prio >= curr_prio:
                            ichor_prio[f"pdf_{k}"] = prio
                            ichor_pdf_map[k] = root_path / f
                            if is_short:
                                ichor_short_pdf_map[k] = root_path / f
                            else:
                                ichor_standard_pdf_map[k] = root_path / f
                            if is_plasma:
                                ichor_plasma_pdf_map[k] = root_path / f
                            if is_tumor:
                                ichor_tumor_pdf_map[k] = root_path / f
            elif f_lower.endswith("params.txt") and ("ichor" in root_path.name.lower() or "ichorcna" in str(root_path).lower()):
                file_stem = f[:-10] if f_lower.endswith("params.txt") else f
                for sub_sfx in [".short", ".standard", "_short", "_standard"]:
                    if file_stem.endswith(sub_sfx):
                        file_stem = file_stem[:-len(sub_sfx)]
                        break
                sids = [file_stem]
                for d_cand in [root_path.name, root_path.parent.name]:
                    if d_cand.lower() not in {"short", "standard", "ichorcna", "ichorcna_plasma_only", "ichorcna_tumor_only", "plasma", "tumor", "wisp", "amber"}:
                        if d_cand not in sids:
                            sids.append(d_cand)

                p_tokens = []
                try:
                    rel_p = root_path.relative_to(results_path).parts
                    if rel_p:
                        base_pt = rel_p[0]
                        p_tokens.append(base_pt)
                        if "_" in base_pt: p_tokens.append(base_pt.replace("_", "-"))
                        if "-" in base_pt: p_tokens.append(base_pt.replace("-", "_"))
                        m_pad = re.search(r'^(.*?)([-_])0*(\d+)$', base_pt)
                        if m_pad:
                            pref, sep, num = m_pad.groups()
                            for s_sep in ["-", "_"]:
                                p_tokens.append(f"{pref}{s_sep}{num}")
                                p_tokens.append(f"{pref}{s_sep}{int(num):04d}")
                                p_tokens.append(f"{pref}{s_sep}{int(num):02d}")
                except Exception:
                    pass

                tf_val = parse_ichor_params_file(root_path / f)
                if tf_val is not None:
                    r_str = str(root_path).lower()
                    prio = 2 if ("tumor_only" in r_str or "plasma_only" in r_str) else 1
                    is_short = ("short" in f_lower or "short" in r_str or "size_select" in f_lower or "size_select" in r_str or "90_150" in f_lower or "90_150" in r_str)
                    is_plasma = ("plasma" in r_str or "mrd" in r_str or "cfdna" in r_str)
                    is_tumor = ("tumor" in r_str or "t_vs_n" in r_str or "tumour" in r_str)
                    if is_short and is_plasma:
                        prio += 4

                    for s in sids:
                        sample_variants = [s, s.replace("-", "_"), s.replace("_", "-")]
                        all_keys = list(sample_variants)
                        for pt in p_tokens:
                            clean_s = extract_clean_sample_name(s, pt)
                            if clean_s:
                                sample_variants.extend([clean_s, clean_s.replace("-", "_"), clean_s.replace("_", "-")])
                            for sv in sample_variants:
                                all_keys.extend([
                                    (pt, sv),
                                    f"{pt}_{sv}",
                                    f"{pt}-{sv}",
                                ])
                        for k in all_keys:
                            curr_prio = ichor_prio.get(k, 0)
                            if prio >= curr_prio:
                                ichor_prio[k] = prio
                                ichor_map[k] = tf_val
                                if is_short:
                                    ichor_short_map[k] = tf_val
                                else:
                                    ichor_standard_map[k] = tf_val
                                if is_plasma:
                                    ichor_plasma_map[k] = tf_val
                                if is_tumor:
                                    ichor_tumor_map[k] = tf_val

                # Look for companion genomeWide PDF in same directory or child directory
                for cand_stem in [file_stem, root_path.name]:
                    for cand_pdf_name in [f"{cand_stem}_genomeWide.pdf", f"{cand_stem}_genomeWide_all_sols.pdf"]:
                        cand_pdf = root_path / cand_pdf_name
                        if not cand_pdf.is_file():
                            cand_pdf = root_path / cand_stem / cand_pdf_name
                        if cand_pdf.is_file():
                            for s in sids:
                                clean_s = extract_clean_sample_name(s, p_token)
                                all_keys = [s, clean_s]
                                if p_token:
                                    all_keys.extend([(p_token, s), (p_token, clean_s), f"{p_token}_{s}", f"{p_token}_{clean_s}"])
                                for k in all_keys:
                                    if k not in ichor_pdf_map:
                                        ichor_pdf_map[k] = cand_pdf
                            break
            elif f == "tumor_fraction_comparison.tsv":
                p_token = None
                try:
                    rel_p = root_path.relative_to(results_path).parts
                    if rel_p:
                        p_token = rel_p[0]
                except Exception:
                    pass
                try:
                    with open(root_path / f, "r", encoding="utf-8-sig") as cf:
                        reader = csv.DictReader(cf, delimiter="\t")
                        for r in reader:
                            sid = r.get("sample_id")
                            tf = r.get("panel_tumor_fraction") or r.get("provided_tumor_fraction")
                            stype = (r.get("sample_type") or "").lower()
                            is_plasma = ("plasma" in stype or "cfdna" in stype)
                            is_tumor = ("tumor" in stype or "tumour" in stype or stype == "t")
                            if sid and tf and tf != "NA":
                                val = safe_float(tf)
                                if val is not None:
                                    prio = 1  # Fallback priority so direct PoN params.txt (prio=2) takes precedence
                                    clean_s = extract_clean_sample_name(sid, p_token)
                                    all_keys = [sid, clean_s]
                                    if p_token:
                                        all_keys.extend([
                                            (p_token, sid),
                                            (p_token, clean_s),
                                            f"{p_token}_{sid}",
                                            f"{p_token}_{clean_s}",
                                            (p_token, f"{p_token}_{sid}"),
                                            (p_token, f"{p_token}_{clean_s}"),
                                        ])
                                    for k in all_keys:
                                        curr_prio = ichor_prio.get(k, 0)
                                        if prio >= curr_prio:
                                            ichor_prio[k] = prio
                                            ichor_map[k] = val
                                            if is_plasma:
                                                ichor_plasma_map[k] = val
                                            if is_tumor:
                                                ichor_tumor_map[k] = val
                except Exception:
                    pass

    # Deduplicate wisp_files preferring .wisp.summary.tsv over .wisp.tsv
    wisp_by_stem = {}
    for wf in wisp_files:
        stem = wf.name.replace(".wisp.summary.tsv", "").replace(".wisp.tsv", "")
        if stem not in wisp_by_stem or wf.name.endswith(".wisp.summary.tsv"):
            wisp_by_stem[stem] = wf

    return {
        "purple_files": sorted(purple_files),
        "wisp_files": sorted(list(wisp_by_stem.values())),
        "orange_map": orange_map,
        "ichor_map": ichor_map,
        "ichor_tumor_map": ichor_tumor_map,
        "ichor_plasma_map": ichor_plasma_map,
        "ichor_short_map": ichor_short_map,
        "ichor_standard_map": ichor_standard_map,
        "ichor_pdf_map": ichor_pdf_map,
        "ichor_tumor_pdf_map": ichor_tumor_pdf_map,
        "ichor_plasma_pdf_map": ichor_plasma_pdf_map,
        "ichor_short_pdf_map": ichor_short_pdf_map,
        "ichor_standard_pdf_map": ichor_standard_pdf_map,
    }


def parse_purple_purity_tsv(purple_tsv_path):
    """Parse purity, ploidy, gender, and MSI from *.purple.purity.tsv."""
    if not os.path.isfile(purple_tsv_path):
        return {}
    try:
        with open(purple_tsv_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                return {
                    "purity": safe_float(r.get("purity")),
                    "norm_factor": safe_float(r.get("normFactor")),
                    "ploidy": safe_float(r.get("ploidy")),
                    "gender": (r.get("gender") or "").strip(),
                    "status": (r.get("status") or "").strip(),
                    "min_purity": safe_float(r.get("minPurity")),
                    "max_purity": safe_float(r.get("maxPurity")),
                    "ms_status": (r.get("msStatus") or "").strip(),
                }
    except Exception as e:
        print(f"Warning: Failed parsing {purple_tsv_path}: {e}", file=sys.stderr)
    return {}


def parse_purple_qc(qc_path):
    """Parse QC status, Amber mean depth, and contamination metrics from *.purple.qc."""
    if not os.path.isfile(qc_path):
        return {}
    res = {}
    try:
        with open(qc_path, "r", encoding="utf-8-sig") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    k, v = parts[0].strip(), parts[1].strip()
                    res[k] = v
    except Exception:
        pass
    return res


def parse_bamtools_summary(tsv_path):
    """Parse MeanCoverage, MedianCoverage, DuplicatePercent, etc. from *.bam_metric.summary.tsv."""
    if not tsv_path or not os.path.isfile(tsv_path):
        return {}
    try:
        with open(tsv_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                return {
                    "mean_coverage": safe_float(r.get("MeanCoverage")),
                    "median_coverage": safe_float(r.get("MedianCoverage")),
                    "duplicate_rate": safe_float(r.get("DuplicatePercent")),
                    "sd_coverage": safe_float(r.get("StdDevCoverage")),
                    "mad_coverage": safe_float(r.get("MadCoverage")),
                    "unmapped_rate": safe_float(r.get("UnmappedPercent")),
                }
    except Exception:
        pass
    return {}


def resolve_bamtools_coverage(tumor_dir, patient_dir=None, tumor_ids=None, ref_ids=None):
    """
    Locate and parse bamtools summary TSVs for tumor and matched normal.
    Checks:
      - tumor_dir / "bamtools" / <sample> / *.bam_metric.summary.tsv
      - patient_dir / ... / bamtools / ...
    """
    res = {"tumor": {}, "normal": {}}
    tumor_id_set = {str(x).strip() for x in (tumor_ids or []) if x and str(x).strip() != "None"}
    ref_id_set = {str(x).strip() for x in (ref_ids or []) if x and str(x).strip() != "None"}

    candidate_files = []
    if tumor_dir and Path(tumor_dir).is_dir():
        candidate_files.extend(list(Path(tumor_dir).glob("bamtools/*/*.bam_metric.summary.tsv")))
    if patient_dir and Path(patient_dir).is_dir():
        candidate_files.extend(list(Path(patient_dir).glob("*/bamtools/*/*.bam_metric.summary.tsv")))

    seen_files = set()
    unique_candidates = []
    for cf in candidate_files:
        p_str = str(cf.resolve())
        if p_str not in seen_files:
            seen_files.add(p_str)
            unique_candidates.append(cf)

    unmatched = []
    for s in unique_candidates:
        sid = s.parent.name
        metrics = parse_bamtools_summary(s)
        if not metrics:
            continue
        if sid in tumor_id_set:
            res["tumor"] = metrics
        elif sid in ref_id_set:
            res["normal"] = metrics
        else:
            unmatched.append((sid, metrics))

    for sid, metrics in unmatched:
        if not res["tumor"] and ("tumor" in sid.lower() or not res["normal"]):
            res["tumor"] = metrics
        elif not res["normal"]:
            res["normal"] = metrics

    return res


def find_orange_json(tumor_id, purple_dir, patient_dir=None, orange_map=None, patient_id=None):
    """
    Search for *.orange.json corresponding to a tumor sample.
    Checks sibling orange directories, patient directories, and scoped search before any global map fallback.
    """
    tumor_dir = purple_dir.parent if purple_dir else None

    # 1. Local Sibling Directory Checks (Highest Priority, guaranteed to be same sample)
    if tumor_dir and tumor_dir.is_dir():
        c = tumor_dir / "orange" / f"{tumor_id}.orange.json"
        if c.is_file():
            return c

        orange_dir = tumor_dir / "orange"
        if orange_dir.is_dir():
            jsons = [f for f in orange_dir.glob("*.orange.json") if not f.name.startswith(".")]
            if jsons:
                return jsons[0]

    # 2. Scoped Patient Directory Checks (Guaranteed to be same patient)
    if patient_dir and isinstance(patient_dir, Path) and patient_dir.is_dir():
        c = patient_dir / tumor_id / "orange" / f"{tumor_id}.orange.json"
        if c.is_file():
            return c

        c = patient_dir / "orange" / f"{tumor_id}.orange.json"
        if c.is_file():
            return c

        p_orange_dir = patient_dir / "orange"
        if p_orange_dir.is_dir():
            jsons = list(p_orange_dir.glob(f"*{tumor_id}*.orange.json"))
            if not jsons:
                jsons = list(p_orange_dir.glob("*.orange.json"))
            if jsons:
                return jsons[0]

        found = [f for f in patient_dir.glob(f"**/{tumor_id}.orange.json") if "/work/" not in str(f) and "/tmp/" not in str(f)]
        if found:
            return found[0]

        found_any = [f for f in patient_dir.glob("**/*.orange.json") if "/work/" not in str(f) and "/tmp/" not in str(f)]
        if found_any:
            return found_any[0]

    # 3. Global Map Fallback (Only if scoped to patient_id or unique)
    if orange_map:
        if patient_id and (patient_id, tumor_id) in orange_map:
            return orange_map[(patient_id, tumor_id)]
        if tumor_id in orange_map and orange_map[tumor_id] is not None:
            cand = orange_map[tumor_id]
            if not patient_id or patient_id in str(cand):
                return cand

    return None


def find_all_ichor_pdfs(sample_id, patient_dir=None, discovered=None, patient_id=None, raw_sample_id=None, is_plasma=None):
    """
    Search for both standard (unfiltered 1X) and short-fragment (90-150bp) ichorCNA plots.
    Robustly handles hyphen/underscore variations, zero-padding, nested subdirectories (short/standard),
    and scoped patient discovery.
    Returns dict:
      {
         "standard": Path or None,
         "short": Path or None,
         "primary": Path or None  # preferred plot (short for plasma if available; standard for tumor)
      }
    """
    patient_tokens = set()
    if patient_id:
        patient_tokens.update([patient_id, patient_id.replace("-", "_"), patient_id.replace("_", "-")])
        m_pad = re.search(r'^(.*?)([-_])0*(\d+)$', patient_id)
        if m_pad:
            pref, sep, num = m_pad.groups()
            for s_sep in ["-", "_"]:
                patient_tokens.add(f"{pref}{s_sep}{num}")
                patient_tokens.add(f"{pref}{s_sep}{int(num):04d}")
                patient_tokens.add(f"{pref}{s_sep}{int(num):02d}")

    raw_tokens = [s for s in [sample_id, raw_sample_id] if s]
    sample_tokens = set()
    for s in raw_tokens:
        sample_tokens.add(s)
        sample_tokens.add(s.replace("-", "_"))
        sample_tokens.add(s.replace("_", "-"))
        for pt in list(patient_tokens):
            clean = extract_clean_sample_name(s, pt)
            if clean:
                sample_tokens.add(clean)
                sample_tokens.add(clean.replace("-", "_"))
                sample_tokens.add(clean.replace("_", "-"))
            s_clean = re.sub(rf"^{re.escape(pt)}[-_]?", "", s, flags=re.IGNORECASE)
            if s_clean and s_clean != s:
                sample_tokens.add(s_clean)
                sample_tokens.add(s_clean.replace("-", "_"))
                sample_tokens.add(s_clean.replace("_", "-"))

    std_cand = None
    sh_cand = None

    # 1. Search scoped patient directories
    cand_pdirs = []
    if patient_dir and isinstance(patient_dir, (str, Path)):
        p_path = Path(patient_dir)
        cand_pdirs.append(p_path)
        p_parent = p_path.parent
        for pt in patient_tokens:
            cand_pdirs.append(p_parent / pt)
            cand_pdirs.append(p_path / pt)

    cand_pdirs = [d for d in cand_pdirs if d.is_dir()]

    if is_plasma is True:
        std_base_names = ["ichorCNA_plasma_vs_n", "ichorCNA_plasma_only", "ichorCNA", "plasma"]
        sh_base_names = ["ichorCNA_plasma_vs_n_short", "ichorCNA_plasma_only_short", "ichorCNA_short", "ichorCNA_size_select", "short"]
    elif is_plasma is False:
        std_base_names = ["ichorCNA_t_vs_n", "ichorCNA_tumor_only", "ichorCNA", "tumor"]
        sh_base_names = ["ichorCNA_t_vs_n_short", "ichorCNA_tumor_only_short", "ichorCNA_short"]
    else:
        std_base_names = ["ichorCNA_plasma_vs_n", "ichorCNA_t_vs_n", "ichorCNA_plasma_only", "ichorCNA_tumor_only", "ichorCNA", "plasma", "tumor"]
        sh_base_names = ["ichorCNA_plasma_vs_n_short", "ichorCNA_t_vs_n_short", "ichorCNA_plasma_only_short", "ichorCNA_short", "short"]

    for p_path in cand_pdirs:
        # Search direct candidate directory combinations
        for bn in (sh_base_names + std_base_names):
            for s in sample_tokens:
                sub_checks = [
                    p_path / bn / s,
                    p_path / bn / s / "short",
                    p_path / bn / s / "standard",
                    p_path / s / bn,
                    p_path / s / bn / "short",
                    p_path / s / bn / "standard",
                    p_path / "plasma" / s / bn,
                    p_path / "plasma" / s,
                    p_path / s,
                ]
                for cand_d in sub_checks:
                    if cand_d.is_dir():
                        for f in cand_d.glob("*genomeWide*.pdf"):
                            if f.name.startswith("._"): continue
                            fn_l = f.name.lower()
                            r_str = str(f).lower()
                            if "short" in fn_l or "short" in r_str or "size_select" in fn_l or "90_150" in fn_l:
                                if not sh_cand: sh_cand = f
                            else:
                                if not std_cand: std_cand = f

        # If still missing either, recursively glob within this patient's directory
        if not sh_cand or not std_cand:
            for f in p_path.glob("**/*genomeWide*.pdf"):
                if f.name.startswith("._"): continue
                fn_l = f.name.lower()
                r_str = str(f).lower()
                if any(tok.lower() in fn_l or tok.lower() in r_str for tok in sample_tokens):
                    if "short" in fn_l or "short" in r_str or "size_select" in fn_l or "90_150" in fn_l:
                        if not sh_cand: sh_cand = f
                    else:
                        if not std_cand: std_cand = f

    # 2. Check maps from discovered
    if discovered and isinstance(discovered, dict):
        sh_map = discovered.get("ichor_short_pdf_map", {})
        std_map = discovered.get("ichor_standard_pdf_map", {})
        gen_map = discovered.get("ichor_pdf_map", {})

        def lookup_map(the_map):
            if not the_map: return None
            for pt in patient_tokens:
                for s in sample_tokens:
                    cand = the_map.get((pt, s)) or the_map.get(f"{pt}_{s}") or the_map.get(f"{pt}-{s}")
                    if cand and os.path.isfile(cand):
                        return Path(cand)
            for s in sample_tokens:
                cand = the_map.get(s)
                if cand and os.path.isfile(cand):
                    return Path(cand)
            return None

        if not sh_cand:
            sh_cand = lookup_map(sh_map)
        if not std_cand:
            std_cand = lookup_map(std_map)
        if not std_cand and not sh_cand:
            gen_cand = lookup_map(gen_map)
            if gen_cand:
                gc_str = str(gen_cand).lower()
                if "short" in gc_str or "size_select" in gc_str:
                    sh_cand = gen_cand
                else:
                    std_cand = gen_cand

    # 3. Sibling recovery: if standard is found but short is missing (or vice versa), check sibling directory
    if sh_cand and not std_cand:
        std_siblings = [
            sh_cand.parent.parent / "standard" / sh_cand.name.replace("_short_", "_standard_").replace(".short.", ".standard."),
            sh_cand.parent.parent / "standard" / sh_cand.name.replace("_short_", "_"),
            sh_cand.parent / sh_cand.name.replace("_short_", "_standard_"),
            sh_cand.parent / sh_cand.name.replace("_short_", "_"),
        ]
        for sib in std_siblings:
            if sib.is_file():
                std_cand = sib
                break
    elif std_cand and not sh_cand:
        sh_siblings = [
            std_cand.parent.parent / "short" / std_cand.name.replace("_standard_", "_short_").replace(".standard.", ".short."),
            std_cand.parent.parent / "short" / f"{std_cand.stem.replace('_standard', '')}_short_genomeWide.pdf",
            std_cand.parent / std_cand.name.replace("_standard_", "_short_"),
        ]
        for sib in sh_siblings:
            if sib.is_file():
                sh_cand = sib
                break

    # Determine primary plot
    if is_plasma:
        primary = sh_cand or std_cand
    else:
        primary = std_cand or sh_cand

    return {
        "standard": std_cand,
        "short": sh_cand,
        "primary": primary
    }


def find_ichor_pdf(sample_id, patient_dir=None, ichor_pdf_map=None, patient_id=None, raw_sample_id=None, is_plasma=None):
    """Backwards-compatible wrapper returning primary ichorCNA plot."""
    disc = {"ichor_pdf_map": ichor_pdf_map} if ichor_pdf_map else None
    res = find_all_ichor_pdfs(sample_id, patient_dir=patient_dir, discovered=disc, patient_id=patient_id, raw_sample_id=raw_sample_id, is_plasma=is_plasma)
    return res.get("primary")


def stage_and_encode_pdf(pdf_path, out_dir, subdir_name, dest_filename):
    """
    Helper to copy a PDF into out_dir/subdir_name and optionally base64-encode it.
    Returns (rel_path, b64_string).
    """
    if not pdf_path or not os.path.isfile(pdf_path):
        return "", ""
    rel_path = ""
    b64_str = ""
    try:
        if os.path.getsize(pdf_path) < 15 * 1024 * 1024:
            with open(pdf_path, "rb") as f_pdf:
                b64_str = base64.b64encode(f_pdf.read()).decode("ascii")
    except Exception as e:
        print(f"Warning: Failed reading PDF {pdf_path} for base64: {e}", file=sys.stderr)

    if out_dir:
        target_dir = Path(out_dir) / subdir_name
        try:
            target_dir.mkdir(parents=True, exist_ok=True)
            target_file = target_dir / dest_filename
            if not target_file.exists() or target_file.stat().st_size != os.path.getsize(pdf_path):
                shutil.copy2(pdf_path, target_file)
            rel_path = f"{subdir_name}/{dest_filename}"
        except Exception as e:
            print(f"Warning: Failed staging PDF {pdf_path}: {e}", file=sys.stderr)
            rel_path = os.path.relpath(pdf_path, out_dir)
    else:
        rel_path = str(pdf_path)

    return rel_path, b64_str


def find_cuppa_tsv(tumor_id, purple_dir, patient_dir=None, patient_id=None):
    """
    Search for *.cuppa.pred_summ.tsv corresponding to a tumor sample.
    Checks sibling cuppa directories and patient directories.
    """
    tumor_dir = purple_dir.parent if purple_dir else None

    # 1. Local Sibling Directory Checks
    if tumor_dir and tumor_dir.is_dir():
        c = tumor_dir / "cuppa" / f"{tumor_id}.cuppa.pred_summ.tsv"
        if c.is_file():
            return c

        cuppa_dir = tumor_dir / "cuppa"
        if cuppa_dir.is_dir():
            tsvs = [f for f in cuppa_dir.glob("*.cuppa.pred_summ.tsv") if not f.name.startswith(".")]
            if tsvs:
                return tsvs[0]

    # 2. Scoped Patient Directory Checks
    if patient_dir and isinstance(patient_dir, Path) and patient_dir.is_dir():
        c = patient_dir / tumor_id / "cuppa" / f"{tumor_id}.cuppa.pred_summ.tsv"
        if c.is_file():
            return c

        c = patient_dir / "cuppa" / f"{tumor_id}.cuppa.pred_summ.tsv"
        if c.is_file():
            return c

        p_cuppa_dir = patient_dir / "cuppa"
        if p_cuppa_dir.is_dir():
            tsvs = list(p_cuppa_dir.glob(f"*{tumor_id}*.cuppa.pred_summ.tsv"))
            if not tsvs:
                tsvs = list(p_cuppa_dir.glob("*.cuppa.pred_summ.tsv"))
            if tsvs:
                return tsvs[0]

        found = [f for f in patient_dir.glob(f"**/{tumor_id}.cuppa.pred_summ.tsv") if "/work/" not in str(f) and "/tmp/" not in str(f)]
        if found:
            return found[0]

        found_any = [f for f in patient_dir.glob("**/*.cuppa.pred_summ.tsv") if "/work/" not in str(f) and "/tmp/" not in str(f)]
        if found_any:
            return found_any[0]

    return None


def parse_cuppa_pred_summ(cuppa_tsv_path):
    """
    Extract CUPPA cancer type predictions from *.cuppa.pred_summ.tsv.
    Returns:
      {
        "cancer_type": str,
        "likelihood": float,
        "predictions": list of dicts [{"cancerType": str, "likelihood": float}]
      }
    """
    if not cuppa_tsv_path or not os.path.isfile(cuppa_tsv_path):
        return {}

    try:
        with open(cuppa_tsv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        target_row = None
        for r in rows:
            if r.get("clf_name") == "dna_combined":
                target_row = r
                break
        if not target_row and rows:
            target_row = rows[0]

        if target_row:
            cancer_type = target_row.get("pred_class_1", "")
            likelihood = safe_float(target_row.get("pred_prob_1"))
            preds = []
            for i in (1, 2, 3):
                pc = target_row.get(f"pred_class_{i}")
                pp = safe_float(target_row.get(f"pred_prob_{i}"))
                if pc:
                    preds.append({"cancerType": pc, "likelihood": pp})
            return {
                "cancer_type": cancer_type,
                "likelihood": likelihood,
                "predictions": preds,
            }
    except Exception as e:
        print(f"Warning: Failed parsing CUPPA TSV {cuppa_tsv_path}: {e}", file=sys.stderr)
    return {}


def map_effect_to_cbio(effect_str, variant_type=""):
    """
    Map Orange/SnpEff/Ensembl variant consequence or driver call
    to standard cBioPortal mutation categories.
    """
    eff = (effect_str or "").upper()
    vtype = (variant_type or "").upper()
    if any(x in eff for x in ["NONSENSE", "FRAMESHIFT", "STOP"]):
        return "Truncating", "Frameshift" if "FRAMESHIFT" in eff else "Nonsense"
    if "SPLICE" in eff:
        return "Splice", "Splice Site"
    if "INFRAME" in eff or "IN_FRAME" in eff:
        return "In-frame InDel", "In-frame InDel"
    if "MISSENSE" in eff:
        return "Missense", "Missense"
    if "PROMOTER" in eff or "PROMOTER" in vtype:
        return "Promoter", "Promoter"
    if "AMP" in vtype or "GAIN" in vtype:
        return "Amplification", "Amplification"
    if "DEL" in vtype or "HOMDEL" in vtype:
        return "Deep Deletion", "Deep Deletion"
    if "FUSION" in vtype:
        return "Fusion", "Fusion"
    return "Mutation", eff.capitalize() if eff else "Mutation"


def extract_orange_variant(v, is_germline=False):
    """
    Extract a structured variant record with HGVS protein, HGVS coding,
    and VAF from an Orange JSON variant object.
    """
    gene = (v.get("gene") or "").strip().upper()
    if not gene:
        return None
    ci = v.get("canonicalImpact") or {}
    eff = v.get("worstCodingEffect") or ci.get("codingEffect") or ""
    alt_type, subtype = map_effect_to_cbio(eff, v.get("type", ""))

    hgvs_p = (ci.get("hgvsProteinImpact") or "").strip()
    hgvs_c = (ci.get("hgvsCodingImpact") or "").strip()

    if not hgvs_p:
        for oi in v.get("otherImpacts") or []:
            if oi.get("hgvsProteinImpact"):
                hgvs_p = oi.get("hgvsProteinImpact").strip()
                break
    if not hgvs_c:
        for oi in v.get("otherImpacts") or []:
            if oi.get("hgvsCodingImpact"):
                hgvs_c = oi.get("hgvsCodingImpact").strip()
                break

    if not hgvs_p:
        hgvs_p = (v.get("canonicalHgvsProtein") or v.get("hgvsProtein") or "").strip()
    if not hgvs_c:
        hgvs_c = (v.get("canonicalHgvsCoding") or v.get("hgvsCoding") or "").strip()

    # Normal depth and status unknown check
    nd = v.get("normalDepth") or {}
    norm_dp = safe_int(nd.get("totalReadCount")) if nd else safe_int(v.get("normalTotalReadCount"))
    norm_ad = safe_int(nd.get("alleleReadCount")) if nd else safe_int(v.get("normalAlleleReadCount"))
    is_status_unknown = (norm_dp is not None and norm_dp < 10)
    status_flag = "Low normal coverage (<10x)" if is_status_unknown else ""

    # Prioritize raw sequencing read depth for accurate VAF determination
    vaf = None
    tot = None
    alt = None
    td = v.get("tumorDepth") or {}
    if td:
        tot = safe_int(td.get("totalReadCount"))
        alt = safe_int(td.get("alleleReadCount"))
        if tot and tot > 0 and alt is not None:
            vaf = alt / tot
    if vaf is None and "totalReadCount" in v and "alleleReadCount" in v:
        tot = safe_int(v.get("totalReadCount"))
        alt = safe_int(v.get("alleleReadCount"))
        if tot and tot > 0 and alt is not None:
            vaf = alt / tot

    if vaf is None:
        adj_vaf = safe_float(v.get("adjustedVAF"))
        if adj_vaf is not None:
            if adj_vaf <= 1.0:
                vaf = adj_vaf
            elif adj_vaf <= 100.0:
                vaf = adj_vaf / 100.0

    return {
        "gene": gene,
        "alt_type": alt_type,
        "subtype": subtype,
        "driver": "MUTATION",
        "likelihood": 1.0,
        "hgvs_protein": hgvs_p,
        "hgvs_coding": hgvs_c,
        "vaf": vaf if vaf is not None else (0.5 if is_germline else 0.0),
        "depth_total": tot,
        "depth_alt": alt,
        "normal_dp": norm_dp,
        "normal_ad": norm_ad,
        "is_status_unknown": is_status_unknown,
        "status_flag": status_flag,
        "variant_copy_number": safe_float(v.get("variantCopyNumber")),
        "adjusted_vaf": safe_float(v.get("adjustedVAF")),
        "biallelic": v.get("biallelic", False),
        "coding_effect": eff,
        "chromosome": (v.get("chromosome") or "").strip(),
        "position": safe_int(v.get("position")),
        "ref": (v.get("ref") or "").strip(),
        "alt": (v.get("alt") or "").strip(),
        "transcript": (ci.get("transcript") or "").strip(),
        "is_germline": is_germline,
        "source": "orange_reportable",
    }


def parse_orange_json(orange_path):
    """
    Extract comprehensive clinical, QC, and genomic metrics from *.orange.json:
      - TMB (tumorMutationalBurdenPerMb, tumorMutationalBurdenStatus, tumorMutationalLoad)
      - MSI (microsatelliteStatus, microsatelliteIndelsPerMb)
      - Whole Genome Duplication (wholeGenomeDuplication)
      - Purity & Ploidy & Fitted Method
      - Tumor Coverage (meanCoverage, medianCoverage, sdCoverage, pctExcDupe)
      - Reference/Normal Coverage (meanCoverage, medianCoverage, pctExcDupe)
      - Quality & Flags (amberMeanDepth, contamination, status flags)
      - Driver Mutations & Reportable Variants (somatic & germline with HGVS protein, coding, VAF)
      - Driver Fusions, Disruptions, and CNAs
      - Complete set of all Orange driver genes
    """
    if not orange_path or not os.path.isfile(orange_path):
        return {}

    try:
        with open(orange_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        purple = data.get("purple", {})
        fit = purple.get("fit", {})
        fit_qc = fit.get("qc", {})
        chars = purple.get("characteristics", {})
        t_sample = data.get("tumorSample", {})
        t_metrics = t_sample.get("metrics", {})
        t_flag = t_sample.get("flagstat", {})
        r_sample = data.get("refSample", {})
        r_metrics = r_sample.get("metrics", {})
        r_flag = r_sample.get("flagstat", {})

        # CUPPA cancer type predictions
        cuppa = data.get("cuppa", {})
        cuppa_best = cuppa.get("bestPrediction", {}) if isinstance(cuppa, dict) else {}
        cuppa_cancer_type = cuppa_best.get("cancerType", "")
        cuppa_likelihood = safe_float(cuppa_best.get("likelihood"))
        cuppa_predictions = cuppa.get("predictions", []) if isinstance(cuppa, dict) else []

        qc_status = fit_qc.get("status", [])
        if isinstance(qc_status, list):
            qc_status_str = "; ".join(qc_status) if qc_status else "PASS"
        else:
            qc_status_str = str(qc_status) if qc_status else "PASS"

        # Extract all driver mutations, reportable variants, fusions, and CNAs from Orange report
        somatic_vars = []
        germline_vars = []
        all_driver_genes = set()

        # 1. Somatic variants (reportable + full somaticVariants list with protein changes)
        seen_somatic = set()
        all_somatic_input = purple.get("reportableSomaticVariants", []) + purple.get("somaticVariants", [])
        for v in all_somatic_input:
            extracted = extract_orange_variant(v, is_germline=False)
            if extracted:
                g = extracted["gene"].upper()
                v_sig = (g, extracted.get("chromosome"), extracted.get("position"), extracted.get("alt"))
                if v_sig not in seen_somatic and g not in PHARMACOGENOMIC_GENES:
                    seen_somatic.add(v_sig)
                    somatic_vars.append(extracted)
                    all_driver_genes.add(g)

        # 2. Germline variants (reportable + full germlineVariants list)
        seen_germline = set()
        all_germline_input = purple.get("reportableGermlineVariants", []) + purple.get("germlineVariants", [])
        for v in all_germline_input:
            extracted = extract_orange_variant(v, is_germline=True)
            if extracted:
                g = extracted["gene"].upper()
                v_sig = (g, extracted.get("chromosome"), extracted.get("position"), extracted.get("alt"))
                if v_sig not in seen_germline and g not in PHARMACOGENOMIC_GENES:
                    seen_germline.add(v_sig)
                    germline_vars.append(extracted)
                    all_driver_genes.add(g)

        # 3. Somatic and germline driver catalog items
        somatic_drivers = purple.get("somaticDrivers", [])
        for sd in somatic_drivers:
            g = (sd.get("gene") or "").strip().upper()
            if g and g not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(g)

        germline_drivers = purple.get("germlineDrivers", [])
        for gd in germline_drivers:
            g = (gd.get("gene") or "").strip().upper()
            if g and g not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(g)

        # 4. CNA driver lists
        cna_drivers = (
            purple.get("reportableSomaticGainsDels", [])
            + purple.get("nearReportableSomaticGains", [])
            + purple.get("reportableGermlineFullDels", [])
        )
        for cd in cna_drivers:
            g = (cd.get("gene") or "").strip().upper()
            if g and g not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(g)

        # 5. LINX structural variants, disruptions, and fusions
        linx = data.get("linx", {})
        for ld in linx.get("somaticDrivers", []):
            g = (ld.get("gene") or "").strip().upper()
            if g and g not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(g)

        # Disruptions and breakends (support both orange v5 somaticBreakends and legacy disruption keys)
        disruptions = []
        raw_breakends = (
            linx.get("somaticBreakends", [])
            + linx.get("somaticHomozygousDisruptions", [])
            + linx.get("germlineHomozygousDisruptions", [])
            + linx.get("reportableSomaticBreakends", [])
        )
        for dis in raw_breakends:
            g = (dis.get("gene") or "").strip().upper()
            is_disruptive = (
                dis.get("disruptive", False)
                or dis.get("driverType") == "DISRUPTION"
                or dis.get("reportedStatus") in {"REPORTED", "DRIVER"}
                or "DISRUPT" in str(dis.get("type", "")).upper()
            )
            if g and is_disruptive and g not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(g)
                disruptions.append({
                    "gene": g,
                    "type": dis.get("type", "SV"),
                    "driver_type": dis.get("driverType", "DISRUPTION"),
                    "reported_status": dis.get("reportedStatus", "REPORTED"),
                    "transcript": dis.get("transcript", ""),
                    "chromosome": dis.get("chromosome", ""),
                    "position": dis.get("position"),
                    "orientation": dis.get("orientation"),
                    "junction_copy_number": safe_float(dis.get("junctionCopyNumber")),
                    "region_type": dis.get("regionType", ""),
                    "exon_up": dis.get("exonUp"),
                    "exon_down": dis.get("exonDown"),
                    "sv_id": str(dis.get("svId")) if dis.get("svId") is not None else "",
                    "cluster_id": str(dis.get("clusterId")) if dis.get("clusterId") is not None else "",
                    "plot_filename": dis.get("plotFilename", ""),
                    "is_germline": False,
                })

        for gb in linx.get("germlineBreakends", []):
            g = (gb.get("gene") or "").strip().upper()
            if g and (gb.get("disruptive", False) or gb.get("reportedStatus") in {"REPORTED", "DRIVER"}) and g not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(g)
                disruptions.append({
                    "gene": g,
                    "type": gb.get("type", "SV"),
                    "driver_type": gb.get("driverType", "DISRUPTION"),
                    "reported_status": gb.get("reportedStatus", "REPORTED"),
                    "transcript": gb.get("transcript", ""),
                    "chromosome": gb.get("chromosome", ""),
                    "position": gb.get("position"),
                    "orientation": gb.get("orientation"),
                    "junction_copy_number": safe_float(gb.get("junctionCopyNumber")),
                    "region_type": gb.get("regionType", ""),
                    "exon_up": gb.get("exonUp"),
                    "exon_down": gb.get("exonDown"),
                    "sv_id": str(gb.get("svId")) if gb.get("svId") is not None else "",
                    "cluster_id": str(gb.get("clusterId")) if gb.get("clusterId") is not None else "",
                    "plot_filename": gb.get("plotFilename", ""),
                    "is_germline": True,
                })

        # Fusions (support orange v5 fusions with geneUp/geneDown as well as geneStart/geneEnd)
        fusions = []
        raw_fusions = (
            linx.get("fusions", [])
            + linx.get("reportableSomaticFusions", [])
            + linx.get("additionalViableSomaticFusions", [])
        )
        for fus in raw_fusions:
            gs = (fus.get("geneUp") or fus.get("geneStart") or "").strip().upper()
            ge = (fus.get("geneDown") or fus.get("geneEnd") or "").strip().upper()
            fname = f"{gs}_{ge}" if (gs and ge) else (gs or ge)
            if gs and gs not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(gs)
            if ge and ge not in PHARMACOGENOMIC_GENES:
                all_driver_genes.add(ge)
            if gs or ge:
                fusions.append({
                    "name": fname,
                    "type": fus.get("reportedType", "FUSION"),
                    "driver_interpretation": fus.get("driverInterpretation", "HIGH"),
                    "gene_start": gs,
                    "gene_end": ge,
                    "gene_up": gs,
                    "gene_down": ge,
                })

        return {
            "tmb": safe_float(chars.get("tumorMutationalBurdenPerMb")),
            "tmb_status": chars.get("tumorMutationalBurdenStatus", ""),
            "tmb_load": safe_int(chars.get("tumorMutationalLoad")),
            "msi_status": chars.get("microsatelliteStatus", "MSS"),
            "msi_score": safe_float(chars.get("microsatelliteIndelsPerMb")),
            "wgd": chars.get("wholeGenomeDuplication", False),
            "purity": safe_float(fit.get("purity")),
            "ploidy": safe_float(fit.get("ploidy")),
            "fitted_purity_method": fit.get("fittedPurityMethod", ""),
            "amber_mean_depth": safe_float(fit_qc.get("amberMeanDepth")),
            "contamination": safe_float(fit_qc.get("contamination")),
            "qc_status": qc_status_str,
            "deleted_genes_count": safe_int(fit_qc.get("deletedGenes")),
            # Tumor coverage & QC
            "tumor_mean_coverage": safe_float(t_metrics.get("meanCoverage")),
            "tumor_median_coverage": safe_float(t_metrics.get("medianCoverage")),
            "tumor_sd_coverage": safe_float(t_metrics.get("sdCoverage")),
            "tumor_duplication_rate": safe_float(t_metrics.get("pctExcDupe")),
            "tumor_mapped_proportion": safe_float(t_flag.get("mappedProportion")),
            "tumor_unique_reads": safe_int(t_flag.get("uniqueReadCount")),
            # Normal coverage & QC
            "normal_mean_coverage": safe_float(r_metrics.get("meanCoverage")),
            "normal_median_coverage": safe_float(r_metrics.get("medianCoverage")),
            "normal_duplication_rate": safe_float(r_metrics.get("pctExcDupe")),
            "normal_mapped_proportion": safe_float(r_flag.get("mappedProportion")),
            "normal_unique_reads": safe_int(r_flag.get("uniqueReadCount")),
            # Driver mutations, fusions, CNAs, and driver genes
            "sample_id": data.get("sampleId"),
            "ref_sample_id": data.get("referenceId"),
            "somatic_variants": somatic_vars,
            "germline_variants": germline_vars,
            "somatic_drivers": somatic_drivers,
            "germline_drivers": germline_drivers,
            "cna_drivers": cna_drivers,
            "fusions": fusions,
            "disruptions": disruptions,
            "all_driver_genes": sorted(list(all_driver_genes)),
            # CUPPA cancer type predictions
            "cuppa_cancer_type": cuppa_cancer_type,
            "cuppa_likelihood": cuppa_likelihood,
            "cuppa_predictions": cuppa_predictions,
        }
    except Exception as e:
        print(f"Warning: Failed parsing orange JSON {orange_path}: {e}", file=sys.stderr)
        return {}


def parse_somatic_vcf_normal_depth(purple_dir, raw_tumor_id, tumor_id):
    """
    Look up normal depth (DP) and allele depth (AD) from purple.somatic.vcf.gz
    or sage.somatic.vcf.gz to detect variants located in zero/low-coverage normal regions.
    """
    vcf_candidates = [
        purple_dir / f"{raw_tumor_id}.purple.somatic.vcf.gz",
        purple_dir / f"{tumor_id}.purple.somatic.vcf.gz",
        purple_dir.parent / "sage" / f"{raw_tumor_id}.sage.somatic.vcf.gz",
        purple_dir.parent / "sage" / f"{tumor_id}.sage.somatic.vcf.gz",
    ]
    vcf_path = None
    for cand in vcf_candidates:
        if cand.is_file():
            vcf_path = cand
            break
    if not vcf_path:
        for cand in [
            purple_dir / f"{raw_tumor_id}.purple.somatic.vcf",
            purple_dir / f"{tumor_id}.purple.somatic.vcf",
            purple_dir.parent / "sage" / f"{raw_tumor_id}.sage.somatic.vcf",
            purple_dir.parent / "sage" / f"{tumor_id}.sage.somatic.vcf",
        ]:
            if cand.is_file():
                vcf_path = cand
                break
    if not vcf_path:
        return {}

    depth_map = {}
    try:
        opener = gzip.open if str(vcf_path).endswith(".gz") else open
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
            normal_col_idx = 9  # Standard SAGE/PURPLE convention: sample 0 is Reference/Normal
            format_col_idx = 8
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) <= normal_col_idx:
                    continue
                chrom = parts[0].replace("chr", "")
                pos = safe_int(parts[1])
                fmt = parts[format_col_idx].split(":")
                if "DP" not in fmt:
                    continue
                dp_idx = fmt.index("DP")
                ad_idx = fmt.index("AD") if "AD" in fmt else None
                norm_sample = parts[normal_col_idx].split(":")
                norm_dp = safe_int(norm_sample[dp_idx]) if len(norm_sample) > dp_idx else None
                norm_ad = safe_int(norm_sample[ad_idx]) if ad_idx is not None and len(norm_sample) > ad_idx else 0

                gene_val = ""
                if len(parts) > 7:
                    for item in parts[7].split(";"):
                        if "=" in item:
                            ik, iv = item.split("=", 1)
                            if ik in {"GENE", "Gene", "SYMBOL"}:
                                gene_val = iv.strip().upper()
                                break

                v_info = {
                    "normal_dp": norm_dp,
                    "normal_ad": norm_ad,
                    "gene": gene_val,
                }
                depth_map[(chrom, pos)] = v_info
                if gene_val:
                    depth_map[(gene_val, pos)] = v_info
                    if gene_val not in depth_map:
                        depth_map[gene_val] = v_info
    except Exception as e:
        print(f"Warning: Failed reading normal depth from {vcf_path}: {e}", file=sys.stderr)
    return depth_map


def parse_somatic_driver_catalog(catalog_path):
    """
    Parse purple.driver.catalog.somatic.tsv to extract somatic driver mutations
    and high-impact coding variants.
    """
    if not os.path.isfile(catalog_path):
        return []

    drivers = []
    try:
        with open(catalog_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("gene") or "").strip().upper()
                if not gene:
                    continue

                driver_type = (r.get("driver") or "").strip()
                driver_lik = safe_float(r.get("driverLikelihood"))
                missense = safe_int(r.get("missense"))
                nonsense = safe_int(r.get("nonsense"))
                splice = safe_int(r.get("splice"))
                inframe = safe_int(r.get("inframe"))
                frameshift = safe_int(r.get("frameshift"))

                total_muts = missense + nonsense + splice + inframe + frameshift
                reported_status = (r.get("reportedStatus") or "").strip().upper()

                # In HMF PURPLE purple.driver.catalog.somatic.tsv:
                # reportedStatus is REPORTED for genuine driver events and NOT_REPORTED for candidate/panel genes with non-driver LOH or uncalled events.
                if reported_status in {"NOT_REPORTED", "UNREPORTED", "NON_DRIVER"}:
                    is_driver = False
                elif reported_status in {"REPORTED", "DRIVER", "PASS"}:
                    is_driver = (total_muts > 0 or driver_lik >= 0.1 or driver_type in {"DEL", "AMP", "PARTIAL_AMP"} or "PROMOTER" in driver_type.upper())
                else:
                    is_driver = (driver_lik >= 0.5 and total_muts > 0) or (driver_type in {"DEL", "AMP", "PARTIAL_AMP"} and driver_lik >= 0.5)

                # Strict check: LOH without coding mutations is NEVER a somatic sequence mutation
                if driver_type == "LOH" and total_muts == 0:
                    is_driver = False

                if not is_driver:
                    continue

                # Determine primary alteration classification according to cBioPortal taxonomy
                alt_type = "Mutation"
                subtype = "Mutation"
                if frameshift > 0 or nonsense > 0:
                    alt_type = "Truncating"
                    subtype = "Frameshift" if frameshift > 0 else "Nonsense"
                elif splice > 0:
                    alt_type = "Splice"
                    subtype = "Splice Site"
                elif inframe > 0:
                    alt_type = "In-frame InDel"
                    subtype = "In-frame InDel"
                elif missense > 0:
                    alt_type = "Missense"
                    subtype = "Missense"
                elif "PROMOTER" in driver_type.upper() or "PROMOTER" in (r.get("category") or "").upper():
                    alt_type = "Promoter"
                    subtype = "Promoter"
                elif "AMP" in driver_type.upper():
                    alt_type = "Amplification"
                    subtype = "Amplification"
                elif "DEL" in driver_type.upper():
                    alt_type = "Deep Deletion"
                    subtype = "Deep Deletion"
                else:
                    if total_muts == 0:
                        continue
                    alt_type = "Mutation"
                    subtype = "Mutation"

                vaf = safe_float(r.get("vaf") or r.get("tumorVaf"))
                hgvs_p = (r.get("canonicalHgvsProtein") or r.get("hgvsProtein") or "").strip()
                hgvs_c = (r.get("canonicalHgvsCoding") or r.get("hgvsCoding") or "").strip()

                drivers.append({
                        "gene": gene,
                        "alt_type": alt_type,
                        "subtype": subtype,
                        "driver": driver_type,
                        "likelihood": driver_lik,
                        "vaf": vaf,
                        "hgvs_protein": hgvs_p,
                        "hgvs_coding": hgvs_c,
                        "transcript": (r.get("canonicalTranscript") or "").strip(),
                        "chromosome": (r.get("chromosome") or r.get("chr") or "").strip(),
                        "start": safe_int(r.get("start")),
                        "is_germline": False,
                    })
    except Exception as e:
        print(f"Warning: Failed parsing driver catalog {catalog_path}: {e}", file=sys.stderr)

    return drivers


def parse_germline_driver_catalog(catalog_path):
    """
    Parse purple.driver.catalog.germline.tsv to extract germline driver mutations.
    """
    if not os.path.isfile(catalog_path):
        return []

    drivers = []
    try:
        with open(catalog_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("gene") or "").strip().upper()
                if not gene:
                    continue

                driver_type = (r.get("driver") or "").strip()
                driver_lik = safe_float(r.get("driverLikelihood"))
                missense = safe_int(r.get("missense"))
                nonsense = safe_int(r.get("nonsense"))
                splice = safe_int(r.get("splice"))
                inframe = safe_int(r.get("inframe"))
                frameshift = safe_int(r.get("frameshift"))

                total_muts = missense + nonsense + splice + inframe + frameshift
                reported_status = (r.get("reportedStatus") or "").strip().upper()

                # For germline driver catalog, PURPLE explicitly marks non-drivers with reportedStatus == "NOT_REPORTED".
                # A true germline driver mutation must be reported, or have high driverLikelihood (>= 0.5) with non-zero mutations.
                if reported_status in {"NOT_REPORTED", "UNREPORTED", "NON_DRIVER"}:
                    is_driver = False
                elif reported_status in {"REPORTED", "DRIVER", "PASS"}:
                    is_driver = (total_muts > 0 or (driver_lik is not None and driver_lik >= 0.1) or any(x in driver_type for x in ["DEL", "AMP", "DELETION"]))
                else:
                    is_driver = ((driver_lik is not None and driver_lik >= 0.5) and total_muts > 0)

                is_cnv = False
                alt_type = "Mutation"
                subtype = "Mutation"
                if frameshift > 0 or nonsense > 0:
                    alt_type = "Truncating"
                    subtype = "Frameshift" if frameshift > 0 else "Nonsense"
                elif splice > 0:
                    alt_type = "Splice"
                    subtype = "Splice Site"
                elif inframe > 0:
                    alt_type = "In-frame InDel"
                    subtype = "In-frame InDel"
                elif missense > 0:
                    alt_type = "Missense"
                    subtype = "Missense"
                elif "PROMOTER" in driver_type.upper() or "PROMOTER" in (r.get("category") or "").upper():
                    alt_type = "Promoter"
                    subtype = "Promoter"
                elif "DEL" in driver_type.upper():
                    alt_type = "Deep Deletion"
                    subtype = "Germline Deletion"
                    is_cnv = True
                elif "AMP" in driver_type.upper():
                    alt_type = "Amplification"
                    subtype = "Germline Amplification"
                    is_cnv = True

                vaf = None
                if not is_cnv:
                    vaf = safe_float(r.get("vaf") or r.get("tumorVaf") or 0.5)
                hgvs_p = (r.get("canonicalHgvsProtein") or r.get("hgvsProtein") or "").strip()
                hgvs_c = (r.get("canonicalHgvsCoding") or r.get("hgvsCoding") or "").strip()

                if is_driver:
                    drivers.append({
                        "gene": gene,
                        "alt_type": alt_type,
                        "subtype": subtype,
                        "driver": driver_type,
                        "likelihood": driver_lik,
                        "vaf": vaf,
                        "hgvs_protein": hgvs_p,
                        "hgvs_coding": hgvs_c,
                        "transcript": (r.get("canonicalTranscript") or "").strip(),
                        "chromosome": (r.get("chromosome") or r.get("chr") or "").strip(),
                        "start": safe_int(r.get("start")),
                        "is_germline": True,
                        "is_cnv": is_cnv,
                        "total_muts": total_muts,
                    })
    except Exception as e:
        print(f"Warning: Failed parsing germline driver catalog {catalog_path}: {e}", file=sys.stderr)

    return drivers


def parse_linx_svs(svs_path):
    """
    Parse linx.svs.tsv to index structural variant records by svId.
    Provides breakpoint coordinates, orientation, junction copy number, and gene partner annotations.
    """
    if not svs_path or not os.path.isfile(svs_path):
        return {}

    sv_map = {}
    try:
        with open(svs_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                sv_id = (r.get("svId") or "").strip()
                if sv_id:
                    sv_map[sv_id] = {
                        "sv_id": sv_id,
                        "type": (r.get("type") or "").strip(),
                        "coords_start": (r.get("coordsStart") or "").strip(),
                        "coords_end": (r.get("coordsEnd") or "").strip(),
                        "gene_start": (r.get("geneStart") or "").strip().upper(),
                        "gene_end": (r.get("geneEnd") or "").strip().upper(),
                        "jcn_min": safe_float(r.get("junctionCopyNumberMin")),
                        "jcn_max": safe_float(r.get("junctionCopyNumberMax")),
                        "cluster_id": (r.get("clusterId") or "").strip(),
                    }
    except Exception as e:
        print(f"Warning: Failed parsing LINX SVs {svs_path}: {e}", file=sys.stderr)

    return sv_map


def parse_linx_driver_catalog(catalog_path):
    """
    Parse linx.driver.catalog.tsv (or linx.drivers.tsv) to extract driver
    structural variants, gene disruptions, and complex SV alterations.
    """
    if not os.path.isfile(catalog_path):
        return []

    drivers = []
    try:
        with open(catalog_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("gene") or "").strip().upper()
                if not gene:
                    continue

                driver_type = (r.get("driver") or r.get("eventType") or "").strip().upper()
                reported_status = (r.get("reportedStatus") or "").strip().upper()
                driver_lik = safe_float(r.get("driverLikelihood"))

                if reported_status in {"NOT_REPORTED", "UNREPORTED", "NON_DRIVER"}:
                    is_driver = False
                elif reported_status in {"REPORTED", "DRIVER", "PASS"}:
                    is_driver = True
                else:
                    is_driver = (driver_lik is not None and driver_lik >= 0.5) or driver_type in {"DISRUPTION", "DEL", "AMP", "PARTIAL_AMP"}

                if not is_driver:
                    continue

                # CRITICAL: linx.driver.catalog.tsv is a merged driver catalog containing MUTATIONS, CNVs, and SVs!
                # Small mutations (point mutations, indels) and CNVs (AMP, DEL, LOH) are ALREADY handled by PURPLE!
                if driver_type in {"MUTATION", "POINT_MUTATION", "INDEL", "SNV", "LOH"}:
                    continue
                if driver_type in {"AMP", "PARTIAL_AMP", "DEL"}:
                    continue

                alt_type = "Structural Variant"
                subtype = "Gene Disruption"
                if "DISRUPT" in driver_type:
                    alt_type = "Structural Variant"
                    subtype = "Gene Disruption"
                elif "FUSION" in driver_type:
                    alt_type = "Fusion"
                    subtype = "Fusion"
                else:
                    alt_type = "Structural Variant"
                    subtype = driver_type

                drivers.append({
                    "gene": gene,
                    "alt_type": alt_type,
                    "subtype": subtype,
                    "driver": driver_type,
                    "likelihood": driver_lik if driver_lik is not None else 1.0,
                    "vaf": None,
                    "hgvs_protein": f"SV ({subtype})",
                    "hgvs_coding": "",
                    "transcript": (r.get("transcript") or "").strip(),
                    "chromosome": (r.get("chromosome") or r.get("chr") or "").strip(),
                    "start": None,
                    "is_germline": False,
                })
    except Exception as e:
        print(f"Warning: Failed parsing LINX driver catalog {catalog_path}: {e}", file=sys.stderr)

    return drivers


def parse_wisp_somatic_variants_tsv(snv_file, primary_driver_map=None):
    """
    Parse a .wisp.somatic_variants.tsv file to extract somatic SNVs/InDels
    evaluated in plasma ctDNA.
    """
    if not os.path.isfile(snv_file):
        return []

    primary_driver_map = primary_driver_map or {}
    results = []
    try:
        with open(snv_file, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("Gene") or r.get("gene") or r.get("gene_name") or r.get("SYMBOL") or "").strip().upper()
                chrom = (r.get("Chromosome") or r.get("chr") or r.get("#CHROM") or "").strip()
                pos = safe_int(r.get("Position") or r.get("pos") or r.get("POS"))

                # If gene not directly in column, check if matches a primary driver by chr/pos
                if not gene and primary_driver_map:
                    for pg, pdrv in primary_driver_map.items():
                        if pdrv.get("chromosome") == chrom and pdrv.get("start") == pos:
                            gene = pg
                            break

                if not gene or gene == "NAN":
                    continue

                sample_ad = safe_int(r.get("SampleAD") or r.get("sample_ad") or r.get("AltCount") or 0)
                sample_dp = safe_int(r.get("SampleDP") or r.get("sample_dp") or r.get("Depth") or 0)
                tumor_ad = safe_int(r.get("TumorAD") or 0)
                tumor_dp = safe_int(r.get("TumorDP") or 0)

                vaf = 0.0
                if sample_dp > 0:
                    vaf = sample_ad / sample_dp
                elif r.get("SampleVAF") or r.get("SampleAF") or r.get("PlasmaVAF"):
                    vaf = safe_float(r.get("SampleVAF") or r.get("SampleAF") or r.get("PlasmaVAF"))
                    if vaf > 1.0:
                        vaf = vaf / 100.0

                # Determine mutation type
                alt_type = "Mutation"
                subtype = "Mutation"
                is_status_unknown = False
                status_flag = ""
                normal_dp = None
                normal_ad = 0
                if gene in primary_driver_map:
                    pdrv = primary_driver_map[gene]
                    alt_type = pdrv.get("alt_type", "Mutation")
                    subtype = pdrv.get("subtype", alt_type)
                    is_status_unknown = pdrv.get("is_status_unknown", False)
                    status_flag = pdrv.get("status_flag", "")
                    normal_dp = pdrv.get("normal_dp")
                    normal_ad = pdrv.get("normal_ad", 0)
                else:
                    vtype = (r.get("Type") or r.get("VariantType") or r.get("Consequence") or "").upper()
                    if "INDEL" in vtype or "FRAME" in vtype:
                        alt_type = "Truncating" if "FRAME" in vtype else "In-frame InDel"
                    elif "MISSENSE" in vtype:
                        alt_type = "Missense"
                    elif "SPLICE" in vtype:
                        alt_type = "Splice"

                results.append({
                    "gene": gene,
                    "alt_type": alt_type,
                    "subtype": subtype,
                    "driver": "SOMATIC",
                    "likelihood": 1.0,
                    "vaf": vaf,
                    "sample_ad": sample_ad,
                    "sample_dp": sample_dp,
                    "tumor_ad": tumor_ad,
                    "tumor_dp": tumor_dp,
                    "is_detected": sample_ad > 0,
                    "hgvs_protein": (r.get("canonicalHgvsProtein") or r.get("hgvsProtein") or "").strip(),
                    "hgvs_coding": (r.get("canonicalHgvsCoding") or r.get("hgvsCoding") or "").strip(),
                    "chromosome": chrom,
                    "start": pos,
                    "is_germline": False,
                    "is_status_unknown": is_status_unknown,
                    "status_flag": status_flag,
                    "normal_dp": normal_dp,
                    "normal_ad": normal_ad,
                })
    except Exception as e:
        print(f"Warning: Failed parsing WISP somatic variants {snv_file}: {e}", file=sys.stderr)

    return results


def parse_linx_fusions(fusion_path):
    """
    Parse linx.fusion.tsv to extract high-confidence somatic gene fusions.
    """
    if not os.path.isfile(fusion_path):
        return {}

    fusions_by_gene = defaultdict(list)
    try:
        with open(fusion_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                g_start = (r.get("geneStart") or "").strip().upper()
                g_end = (r.get("geneEnd") or "").strip().upper()
                name = (r.get("name") or "").strip()
                reported = (r.get("reported") or "false").strip().lower() == "true"
                rep_type = (r.get("reportedType") or "").strip().upper()
                phased = (r.get("phased") or "false").strip().lower() == "true"
                reasons = (r.get("reportableReasons") or "").strip()

                is_valid = reported or (rep_type and rep_type != "NONE") or phased or ("INFRAME" in reasons.upper())
                if is_valid and (g_start or g_end):
                    fusion_entry = {
                        "name": name or f"{g_start}_{g_end}",
                        "gene_5p": g_start,
                        "gene_3p": g_end,
                        "reported_type": rep_type,
                    }
                    if g_start:
                        fusions_by_gene[g_start].append(fusion_entry)
                    if g_end and g_end != g_start:
                        fusions_by_gene[g_end].append(fusion_entry)
    except Exception as e:
        print(f"Warning: Failed parsing LINX fusions {fusion_path}: {e}", file=sys.stderr)

    return dict(fusions_by_gene)


# Canonical Cancer Driver & Census Genes for reporting focal CNVs and driver alterations
CANONICAL_CANCER_GENES = {
    "ABL1", "AKT1", "AKT2", "AKT3", "ALK", "APC", "AR", "ARAF", "ARID1A", "ARID1B", "ARID2", "ASXL1", "ATM", "ATRX",
    "B2M", "BAP1", "BCL2", "BCL6", "BCOR", "BRAF", "BRCA1", "BRCA2", "BRIP1", "BTK",
    "CARD11", "CASP8", "CBL", "CCND1", "CCND2", "CCND3", "CCNE1", "CD274", "CDH1", "CDK12", "CDK4", "CDK6",
    "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "CEBPA", "CHEK1", "CHEK2", "CREBBP", "CRLF2", "CSF1R",
    "CTCF", "CTNNA1", "CTNNB1", "DAXX", "DDR2", "DDX3X", "DNMT3A", "DOT1L",
    "EGFR", "EP300", "EPCAM", "EPHA2", "EPHA3", "ERBB2", "ERBB3", "ERBB4", "ERCC2", "ERG", "ESR1", "ETV1", "ETV4", "ETV6", "EZH2",
    "FANCA", "FANCC", "FAS", "FBXW7", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT3", "FOXL2", "FUBP1",
    "GATA1", "GATA2", "GATA3", "GID4", "GNA11", "GNAQ", "GNAS", "GREM1",
    "H3F3A", "HNF1A", "HRAS", "IDH1", "IDH2", "IGF1R", "IKZF1", "INPP4B", "JAK1", "JAK2", "JAK3", "JUN",
    "KDM5C", "KDM6A", "KDR", "KEAP1", "KIT", "KLF4", "KMT2A", "KMT2C", "KMT2D", "KRAS",
    "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAPK1", "MAX", "MCL1", "MDM2", "MDM4", "MED12", "MEN1", "MET", "MITF",
    "MLH1", "MLLT3", "MPL", "MRE11", "MSH2", "MSH3", "MSH6", "MTOR", "MUTYH", "MYC", "MYCL", "MYCN", "MYD88",
    "NBN", "NCOR1", "NF1", "NF2", "NFE2L2", "NOTCH1", "NOTCH2", "NOTCH3", "NPM1", "NRAS", "NSD1",
    "PALB2", "PARP1", "PAX5", "PBRM1", "PDCD1LG2", "PDGFRA", "PDGFRB", "PIK3C2B", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG",
    "PIK3R1", "PIK3R2", "PIM1", "PMS2", "POLE", "POLD1", "POT1", "PPARG", "PPM1D", "PPP2R1A", "PPP2R2A", "PRDM1", "PRKAR1A",
    "PTCH1", "PTEN", "PTPN11", "RAC1", "RAD21", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAF1", "RARA", "RB1",
    "RET", "RHEB", "RHOA", "RICTOR", "RIT1", "RNF43", "ROS1", "RUNX1",
    "SDHA", "SDHB", "SDHC", "SDHD", "SETD2", "SF3B1", "SH2B3", "SLX4", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMARCC2",
    "SMO", "SOCS1", "SOX17", "SOX2", "SOX9", "SPEN", "SPOP", "SRC", "STAG2", "STAT3", "STAT5B", "STK11", "SUFU",
    "TBX3", "TEK", "TERT", "TET2", "TGFBR2", "TIPARP", "TNFAIP3", "TNFRSF14", "TOP1", "TP53", "TSC1", "TSC2",
    "U2AF1", "VHL", "WT1", "XPO1", "ZFHX3", "ZFP36L1", "ZRSR2"
}


def is_driver_cna_candidate(gene: str, user_genes: set = None) -> bool:
    """
    Filter out passenger non-coding RNAs, olfactory receptors, and pseudogenes.
    Returns True only for user genes of interest or recognized canonical cancer driver genes.
    """
    g = gene.strip().upper()
    if user_genes and g in user_genes:
        return True
    if g in CANONICAL_CANCER_GENES:
        return True
    # Exclude pseudogenes, non-coding RNAs, and uncharacterized olfactory arrays
    if g.startswith(("MIR", "SNORD", "RNU", "LINC", "RNA5", "RN7", "Y_RNA", "OR")):
        return False
    if re.search(r"P\d+$", g) or g.endswith("P"):
        return False
    if any(tag in g for tag in ("-AS", "-DT", "-OT")):
        return False
    return False


def parse_cnv_gene_tsv(cnv_path, ploidy=2.0, amp_ratio=3.0, gender="", genes_of_interest=None):
    """
    Parse purple.cnv.gene.tsv for copy number status matching cBioPortal criteria:
      - Amplification (cBioPortal '2'): maxCopyNumber > (amp_ratio * ploidy)
      - Deep Deletion (cBioPortal '-2'): minCopyNumber < 0.5
      - Gain (cBioPortal '1'): maxCopyNumber > (2.0 * ploidy)
      - Loss (cBioPortal '-1'): minMinorAlleleCopyNumber < 0.5
    Only canonical transcript rows (isCanonical == true) are processed.
    ChrY absence in female patients is excluded to avoid false deep deletion calls.
    Focal driver candidacy (is_driver) is flagged to prevent passenger whole-genome loss noise.
    """
    if not os.path.isfile(cnv_path):
        return {}

    cnv_calls = {}
    effective_ploidy = float(ploidy) if ploidy and safe_float(ploidy) > 0 else 2.0
    amp_cutoff = amp_ratio * effective_ploidy
    gain_cutoff = 2.0 * effective_ploidy
    is_female = gender.strip().upper() in {"FEMALE", "F"}
    user_gene_set = {g.strip().upper() for g in (genes_of_interest or []) if g.strip()}

    try:
        with open(cnv_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("gene") or "").strip().upper()
                if not gene:
                    continue

                chrom = (r.get("chromosome") or r.get("chr") or "").strip().replace("chr", "").upper()
                if chrom == "Y" and is_female:
                    continue

                is_canonical = (r.get("isCanonical") or "true").strip().lower() == "true"
                if not is_canonical:
                    continue

                min_cn = safe_float(r.get("minCopyNumber"), default=2.0)
                max_cn = safe_float(r.get("maxCopyNumber"), default=2.0)
                min_minor_cn = safe_float(r.get("minMinorAlleleCopyNumber"), default=1.0)

                state = "Diploid"
                cbio_val = 0
                if min_cn < 0.5:
                    state = "Deep Deletion"
                    cbio_val = -2
                elif max_cn > amp_cutoff:
                    state = "Amplification"
                    cbio_val = 2
                elif max_cn > gain_cutoff:
                    state = "Gain"
                    cbio_val = 1
                elif min_minor_cn < 0.5:
                    state = "Loss"
                    cbio_val = -1

                if state != "Diploid":
                    is_driver = is_driver_cna_candidate(gene, user_gene_set)
                    cnv_calls[gene] = {
                        "state": state,
                        "cbio_val": cbio_val,
                        "min_cn": min_cn,
                        "max_cn": max_cn,
                        "is_driver": is_driver,
                    }
    except Exception as e:
        print(f"Warning: Failed parsing CNV gene file {cnv_path}: {e}", file=sys.stderr)

    return cnv_calls


def parse_ichor_metrics(patient_dir, is_plasma=None):
    """Extract ichorCNA tumor fractions for samples under a patient directory."""
    results = {}
    p_path = Path(patient_dir)
    p_name = p_path.name
    if is_plasma is True:
        ichor_dirs = [
            p_path / "ichorCNA_plasma_vs_n",
            p_path / "ichorCNA_plasma_vs_n_short",
            p_path / "ichorCNA_plasma_only",
            p_path / "ichorCNA_plasma_only_short",
            p_path / "ichorCNA_short",
        ]
        # Discover plasma-specific ichor directories: e.g. <patient_dir>/plasma/<plasma_id>/ichorCNA*
        for cand in p_path.glob("plasma/*/ichorCNA*"):
            if cand.is_dir() and cand not in ichor_dirs:
                ichor_dirs.append(cand)
    elif is_plasma is False:
        ichor_dirs = [
            p_path / "ichorCNA_t_vs_n",
            p_path / "ichorCNA_tumor_only",
            p_path / "ichorCNA_normal_only",
            p_path / "ichorCNA",
        ]
    else:
        ichor_dirs = [
            p_path / "ichorCNA_plasma_vs_n",
            p_path / "ichorCNA_t_vs_n",
            p_path / "ichorCNA_plasma_only",
            p_path / "ichorCNA_tumor_only",
            p_path / "ichorCNA_normal_only",
            p_path / "ichorCNA",
        ]
        for cand in p_path.glob("plasma/*/ichorCNA*"):
            if cand.is_dir() and cand not in ichor_dirs:
                ichor_dirs.append(cand)
        for cand in p_path.glob("*/ichorCNA*"):
            if cand.is_dir() and cand not in ichor_dirs:
                ichor_dirs.append(cand)

    def _register(sid, val):
        if not sid or val is None:
            return
        clean_sid = extract_clean_sample_name(sid, p_name)
        keys_to_set = [sid, clean_sid]
        if p_name:
            keys_to_set.extend([
                f"{p_name}_{sid}",
                f"{p_name}_{clean_sid}",
            ])
        for k in keys_to_set:
            if k not in results or results[k] is None:
                results[k] = val

    for ichor_base in ichor_dirs:
        if not ichor_base.is_dir():
            continue

        # 1. Check sample-specific params files first (direct PoN / params.txt takes priority; short takes priority for plasma)
        try:
            if is_plasma:
                for p_file in ichor_base.glob("**/short/*params.txt"):
                    sid = p_file.parent.parent.name if p_file.parent.name == "short" else p_file.parent.name
                    v = parse_ichor_params_file(p_file)
                    if v is not None:
                        _register(sid, v)
                for p_file in ichor_base.glob("**/*short*params.txt"):
                    sid = p_file.parent.name
                    file_sid = p_file.name.replace(".short.params.txt", "").replace("_short.params.txt", "").replace(".params.txt", "")
                    v = parse_ichor_params_file(p_file)
                    if v is not None:
                        _register(sid, v)
                        if file_sid and file_sid != sid:
                            _register(file_sid, v)
            for p_file in ichor_base.glob("*/*params.txt"):
                sid = p_file.parent.name
                file_sid = p_file.name.replace(".params.txt", "").replace(".params", "")
                v = parse_ichor_params_file(p_file)
                if v is not None:
                    _register(sid, v)
                    if file_sid and file_sid != sid:
                        _register(file_sid, v)
            for p_file in ichor_base.glob("*params.txt"):
                sid = p_file.stem.replace(".params", "")
                v = parse_ichor_params_file(p_file)
                if v is not None:
                    _register(sid, v)
        except Exception:
            pass

        # 2. Check comparison TSV as fallback if present
        comp_tsv = ichor_base / "tumor_fraction_comparison.tsv"
        if comp_tsv.is_file():
            try:
                with open(comp_tsv, "r", encoding="utf-8-sig") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for r in reader:
                        sid = r.get("sample_id")
                        tf = r.get("panel_tumor_fraction") or r.get("provided_tumor_fraction")
                        if sid and tf and tf != "NA":
                            v = safe_float(tf)
                            if v is not None:
                                _register(sid, v)
            except Exception:
                pass

    return results


def parse_wisp_summary_tsv(tsv_path):
    """Parse a single .wisp.summary.tsv file into a dict of typed values."""
    if not os.path.isfile(tsv_path):
        return None
    try:
        with open(tsv_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            if rows:
                return rows[0]
    except Exception:
        pass
    return None


def extract_metadata_from_wisp_path(fp, base_dir, patient_id):
    """Extract tumor_id and longitudinal_sample_id from WISP filename & parent directories."""
    fname = fp.name
    parts = list(fp.parts)
    tumor_id = "Primary_Tumor"

    try:
        w_idx = parts.index("wisp")
        if w_idx > 0:
            cand = parts[w_idx - 1]
            if cand not in {"purity", "output", "results", patient_id}:
                if not re.search(r"[-_](t[0-9]+|fu[0-9]+|plasma|cfdna)$", cand, re.IGNORECASE):
                    tumor_id = cand
    except ValueError:
        pass

    stem = fname.replace(".wisp.summary.tsv", "").replace(".wisp.tsv", "")
    sample_id = stem
    if patient_id and stem.startswith(f"{patient_id}_"):
        sample_id = stem[len(patient_id) + 1:]

    return tumor_id, sample_id


def parse_date_safe(date_str):
    """
    Parses a date string into a datetime object and ISO string (YYYY-MM-DD).
    Supports common date formats (YYYY-MM-DD, YYYY/MM/DD, MM/DD/YYYY, DD/MM/YYYY, etc.).
    """
    if not date_str:
        return None, ""
    s = str(date_str).strip()
    if s.lower() in {"na", "nan", "null", "none", ""}:
        return None, ""

    formats = [
        "%Y-%m-%d",
        "%Y/%m/%d",
        "%m/%d/%Y",
        "%d/%m/%Y",
        "%m-%d-%Y",
        "%d-%m-%Y",
        "%Y%m%d",
        "%b %d, %Y",
        "%b %d %Y",
        "%d-%b-%Y",
        "%d-%b-%y",
        "%Y-%m-%d %H:%M:%S",
    ]
    for fmt in formats:
        try:
            dt = datetime.strptime(s, fmt)
            return dt, dt.strftime("%Y-%m-%d")
        except ValueError:
            pass
    return None, s


def extract_clean_sample_name(sample_id, patient_id=None):
    """
    Extracts a clean, human-readable sample name/token from a raw sample ID or file stem.
    Strips file suffixes (.redux.bam, .bam, _wisp, etc.) and cohort/patient prefixes.
    E.g.:
      TRACEOC_0002_03_LB01-01 -> 03_LB01-01
      OV-P-02-T1 -> T1
      OV-P-02-03 -> 03
      BFIERCE_0001_04 -> 04
      DYNA-P01-PL1 -> PL1
      sample_2 -> sample_2
    """
    if not sample_id:
        return ""
    s = str(sample_id).strip()
    for ext in [".redux.bam", ".bam", ".fastq.gz", ".wisp.summary.tsv", ".wisp.tsv", "_wisp", ".wisp"]:
        if s.endswith(ext):
            s = s[:-len(ext)]

    # Strip patient prefix if present
    if patient_id:
        p_esc = re.escape(str(patient_id).strip())
        s = re.sub(rf"^{p_esc}[-_]?", "", s, flags=re.IGNORECASE)
        d_under = re.escape(str(patient_id).replace("-", "_").strip())
        s = re.sub(rf"^{d_under}[-_]?", "", s, flags=re.IGNORECASE)

    # Strip standard cohort prefixes: e.g. HCC-B-001_, TRACEOC_0002_, BFIERCE_0001_, etc.
    s = re.sub(
        r"^(HCC[_-]B[_-]\d+|TRACEOC[_-]\d+|OV[_-]P[_-]\d+|BFIERCE[_-]\d+|DYNA[_-]P\d+|MOHCCN[_-][A-Za-z0-9]+)[_-]?",
        "",
        s,
        flags=re.IGNORECASE
    )

    return s


def natural_sort_key(s):
    """
    Split string into numeric and non-numeric chunks for natural numerical sorting.
    Returns a tuple of 3-tuples (type_flag, int_val, str_val) so Python 3 never compares int and str.
    Guarantees that 'sample_2' < 'sample_10', '03' < '04' < '18', etc.
    """
    chunks = re.split(r'(\d+)', str(s or ""))
    key = []
    for c in chunks:
        if not c:
            continue
        if c.isdigit():
            key.append((0, int(c), ""))
        else:
            key.append((1, 0, c.lower()))
    return tuple(key)


def register_metadata_record(metadata_map, aliases, rec, extra_keys=None):
    """
    Registers all aliases, stems, clean names, and patient-prefixed variations
    of a sample into metadata_map and aliases. Merges non-empty fields into existing records.
    """
    sample_id = rec.get("sample_id") or ""
    donor_id = rec.get("patient_id") or ""
    clean_sname = rec.get("clean_sample_name") or ""

    base_keys = set()
    for k in [sample_id, clean_sname] + (extra_keys or []):
        if k:
            base_keys.add(k)
            # without aliquot/replicate suffix (e.g. -01, _01)
            k_no_ali = re.sub(r"[-_](01|02|03|\d{1,2})$", "", k)
            if k_no_ali:
                base_keys.add(k_no_ali)

    all_keys = set()
    for k in base_keys:
        all_keys.add(k)
        all_keys.add(k.lower())
        all_keys.add(f"{k}_wisp")
        all_keys.add(f"{k.lower()}_wisp")
        k_no_ext = k.replace(".redux.bam", "").replace(".bam", "").replace(".fastq.gz", "").replace("_wisp", "")
        all_keys.add(k_no_ext)
        all_keys.add(k_no_ext.lower())

        if donor_id:
            all_keys.add(f"{donor_id}_{k}")
            all_keys.add(f"{donor_id}-{k}")
            all_keys.add(f"{donor_id}_{k}".lower())
            all_keys.add(f"{donor_id}-{k}".lower())
            d_u = donor_id.replace("-", "_")
            all_keys.add(f"{d_u}_{k}")
            all_keys.add(f"{d_u}_{k}".lower())

    for k in all_keys:
        if k in metadata_map:
            existing = metadata_map[k]
            for field in ["patient_id", "timepoint", "date", "date_dt", "days", "tissue_type", "tissue_origin", "clean_sample_name"]:
                v_new = rec.get(field)
                v_old = existing.get(field)
                if v_new and (not v_old or v_old in {"Plasma ctDNA", "Plasma"} or (field == "timepoint" and v_old in {existing.get("sample_id"), existing.get("clean_sample_name")})):
                    existing[field] = v_new
        else:
            metadata_map[k] = dict(rec)

    if donor_id and sample_id:
        aliases[sample_id] = donor_id
        aliases[sample_id.lower()] = donor_id


def load_cohort_metadata(results_path, sample_sheets_dir=None, metadata_tsv=None, timepoints_tsv=None):
    """
    Loads metadata mappings across candidate metadata TSV/CSV files, timepoint tables, and samplesheets.
    Maps sample IDs, BAM stems, and patient tokens to:
      {
        "patient_id": donor/patient ID,
        "sample_id": cleaned sample ID,
        "timepoint": explicit or inferred timepoint,
        "date": formatted collection date,
        "date_dt": parsed datetime object,
        "days": numeric days from baseline,
        "tissue_type": tissue type,
        "clean_sample_name": clean sample token
      }
    """
    metadata_map = {}
    aliases = {}

    candidate_files = []
    if timepoints_tsv and os.path.isfile(timepoints_tsv):
        candidate_files.append(Path(timepoints_tsv))
    if metadata_tsv and os.path.isfile(metadata_tsv):
        candidate_files.append(Path(metadata_tsv))
    if sample_sheets_dir and os.path.isdir(sample_sheets_dir):
        candidate_files.extend(list(Path(sample_sheets_dir).glob("**/*.csv")))
        candidate_files.extend(list(Path(sample_sheets_dir).glob("**/*.tsv")))

    # Auto-discover metadata and timepoint files in results_path, its parent, and repo runs/metadata/
    script_dir = Path(__file__).resolve().parent
    repo_meta_dir = script_dir.parent / "runs" / "metadata"
    if repo_meta_dir.is_dir():
        for mf in repo_meta_dir.glob("*.tsv"):
            if not mf.name.startswith(".") and mf not in candidate_files:
                candidate_files.append(mf)
        for mf in repo_meta_dir.glob("*.csv"):
            if not mf.name.startswith(".") and mf not in candidate_files:
                candidate_files.append(mf)

    for parent_dir in [results_path, results_path.parent]:
        for pattern in ["*metadata*.tsv", "*metadata*.csv", "*config*.tsv", "*config*.csv", "*timepoint*.tsv", "*timepoint*.csv", "*timepoints*.tsv", "*timepoints*.csv"]:
            for mf in parent_dir.glob(pattern):
                if not mf.name.startswith(".") and mf not in candidate_files:
                    candidate_files.append(mf)

    for cf in candidate_files:
        if not cf.is_file() or cf.name.startswith("."):
            continue
        try:
            delim = "\t" if cf.suffix.lower() in {".tsv", ".txt"} else ","
            with open(cf, "r", encoding="utf-8-sig", errors="replace") as f:
                sample_lines = [f.readline() for _ in range(5)]
                f.seek(0)
                if sample_lines and "\t" in sample_lines[0]:
                    delim = "\t"
                elif sample_lines and "," in sample_lines[0]:
                    delim = ","

                reader = csv.reader(f, delimiter=delim)
                header_row = next(reader, None)
                if not header_row:
                    continue

                headers_clean = [h.strip().lower() for h in header_row]

                has_sample_header = any(h in headers_clean for h in [
                    "sample id", "sample_id", "sample", "sample_name", "sample name",
                    "plasma_sample_id", "plasma_sample", "plasma_id", "tumor_id", "normal_id",
                    "fastq file name", "bam file name", "filename", "file", "id"
                ])

                if has_sample_header:
                    dict_reader = csv.DictReader(f, fieldnames=header_row, delimiter=delim)
                    for r in dict_reader:
                        r_clean = {k.strip().lower(): (v or "").strip() for k, v in r.items() if k}

                        donor_id = (
                            r_clean.get("donor_id") or r_clean.get("donor") or
                            r_clean.get("patient_id") or r_clean.get("patient") or
                            r_clean.get("subject_id") or r_clean.get("subject") or
                            r_clean.get("individual_id") or ""
                        )

                        sample_id = (
                            r_clean.get("sample id") or r_clean.get("sample_id") or
                            r_clean.get("sample") or r_clean.get("sample_name") or
                            r_clean.get("sample name") or r_clean.get("plasma_sample_id") or
                            r_clean.get("plasma_sample") or r_clean.get("plasma_id") or
                            r_clean.get("tumor_id") or r_clean.get("normal_id") or
                            r_clean.get("id") or ""
                        )

                        bam_fn = (
                            r_clean.get("fastq file name") or r_clean.get("bam file name") or
                            r_clean.get("filename") or r_clean.get("file") or ""
                        )
                        bam_stem = Path(bam_fn).name.replace(".redux.bam", "").replace(".bam", "").replace(".fastq.gz", "") if bam_fn else ""

                        # Explicit timepoint
                        raw_tp = (
                            r_clean.get("timepoint") or r_clean.get("time_point") or
                            r_clean.get("time point") or r_clean.get("tp") or
                            r_clean.get("visit") or r_clean.get("timepoint_label") or
                            r_clean.get("timepoint_name") or ""
                        )

                        # Explicit collection date
                        date_raw = (
                            r_clean.get("date") or r_clean.get("collection_date") or
                            r_clean.get("collection date") or r_clean.get("draw_date") or
                            r_clean.get("draw date") or r_clean.get("sample_date") or
                            r_clean.get("sample date") or r_clean.get("sampling_date") or
                            r_clean.get("sampling date") or r_clean.get("blood_draw_date") or
                            r_clean.get("blood draw date") or r_clean.get("date_collected") or
                            r_clean.get("date collected") or ""
                        )
                        date_dt, date_str = parse_date_safe(date_raw)

                        # Explicit days from baseline
                        days_raw = (
                            r_clean.get("days") or r_clean.get("day") or
                            r_clean.get("days_from_baseline") or r_clean.get("days from baseline") or
                            r_clean.get("days_post_op") or r_clean.get("elapsed_days") or
                            r_clean.get("elapsed days") or r_clean.get("cycle_day") or
                            r_clean.get("cycle day") or ""
                        )
                        days_val = safe_float(days_raw, default=None) if days_raw != "" else None

                        tissue_origin = (r_clean.get("tissue origin") or r_clean.get("tissue_origin") or "").strip()
                        tissue_type = (r_clean.get("tissue type") or r_clean.get("type") or "").strip()

                        # If donor_id is present, strip it from raw_tp
                        if donor_id:
                            p_esc = re.escape(donor_id)
                            raw_tp_clean = re.sub(rf"^{p_esc}[-_]?", "", raw_tp, flags=re.IGNORECASE)
                            d_under = re.escape(donor_id.replace("-", "_"))
                            raw_tp_clean = re.sub(rf"^{d_under}[-_]?", "", raw_tp_clean, flags=re.IGNORECASE)
                        else:
                            raw_tp_clean = raw_tp

                        tp = raw_tp_clean.strip()
                        if not tp:
                            if date_str:
                                tp = date_str
                            elif days_val is not None:
                                tp = f"Day {int(days_val)}"
                            elif tissue_origin.lower() in {"primary", "tumor"}:
                                tp = "Baseline (Tumor)"
                            elif tissue_origin.lower() in {"pl", "plasma"}:
                                m_tp = re.search(r"[-_](T[0-9]+|FU[0-9]+)\b", sample_id, re.IGNORECASE)
                                tp = m_tp.group(1).upper() if m_tp else ""
                        elif tp.lower() in {"tumor", "primary"}:
                            tp = "Baseline (Tumor)"
                            if not tissue_origin:
                                tissue_origin = "Tumor"
                        elif re.match(r"^T\d+$", tp, re.IGNORECASE):
                            tp = tp.upper()
                            if not tissue_origin:
                                tissue_origin = "Plasma"
                                tissue_type = "Plasma ctDNA"

                        clean_sname = extract_clean_sample_name(sample_id or bam_stem, donor_id)

                        rec = {
                            "patient_id": donor_id,
                            "sample_id": sample_id or bam_stem,
                            "timepoint": tp,
                            "date": date_str,
                            "date_dt": date_dt,
                            "days": days_val,
                            "tissue_type": tissue_type,
                            "tissue_origin": tissue_origin,
                            "clean_sample_name": clean_sname,
                        }
                        register_metadata_record(metadata_map, aliases, rec, extra_keys=[bam_stem, raw_tp])
                else:
                    all_rows = [header_row] + list(reader)
                    for row in all_rows:
                        if not row or len(row) < 2:
                            continue
                        col0 = row[0].strip()
                        col1 = row[1].strip()
                        col2 = row[2].strip() if len(row) > 2 else ""

                        if col0.startswith("#") or col0.lower() in {"sample", "sample_id", "sample id"}:
                            continue

                        # Check if col1 is a date
                        dt_col1, date_str_col1 = parse_date_safe(col1)

                        if len(row) >= 3 and dt_col1 is None:
                            # 3 columns without date in col1:
                            # Col 0: patient_id (e.g. HCC-B-171)
                            # Col 1: sample_id  (e.g. HCC-B-171_HCCCFD2_0022_Pl_P_HCC-B-171-T4-P-DNA)
                            # Col 2: timepoint/date/days (e.g. T4)
                            donor_id = col0.replace("_", "-")
                            sample_id = col1
                            time_or_date = col2
                            third_col = ""
                        else:
                            sample_id = col0
                            time_or_date = col1
                            third_col = col2
                            donor_id = ""

                        # Extract donor from time_or_date or sample_id if not yet set
                        if not donor_id:
                            m_p1 = re.search(r"^(OV[-_]P[-_]\d+|TRACEOC[-_]\d+|DYNA[-_]P\d+|BFIERCE[-_]\d+|MOHCCN[-_][A-Za-z0-9]+)", time_or_date, re.IGNORECASE)
                            if m_p1:
                                donor_id = m_p1.group(1).replace("_", "-")
                            else:
                                for s in [time_or_date, sample_id]:
                                    m_p0 = re.search(r"\b(TRACEOC|OV[-_]P|DYNA[-_]P|BFIERCE|MOHCCN)[-_]?(\d+|[A-Za-z0-9]+)\b", s, re.IGNORECASE)
                                    if m_p0:
                                        donor_id = m_p0.group(0).replace("_", "-")
                                        break
                                    # Generalized donor pattern (e.g. HCC-B-171, PAT-01, DONOR-12)
                                    m_gen = re.search(r"^([A-Za-z0-9]+(?:[-_][A-Za-z0-9]+)+)[-_]", s)
                                    if m_gen:
                                        cand = m_gen.group(1)
                                        if not re.search(r"^(lb\d+|cfd\d*|pl|dna|rna|wgs|wes)", cand, re.IGNORECASE):
                                            donor_id = cand
                                            break

                        # Check date parsing
                        dt_obj, date_str = parse_date_safe(time_or_date)
                        dt_obj2, date_str2 = parse_date_safe(third_col) if third_col else (None, "")

                        tp = ""
                        days_val = None
                        tissue_origin = ""
                        tissue_type = ""

                        if dt_obj:
                            tp = date_str
                            days_val = safe_float(third_col, default=None) if third_col else None
                        else:
                            raw_val = time_or_date
                            if donor_id:
                                p_esc = re.escape(donor_id)
                                raw_val = re.sub(rf"^{p_esc}[-_]?", "", raw_val, flags=re.IGNORECASE)
                                d_under = re.escape(donor_id.replace("-", "_"))
                                raw_val = re.sub(rf"^{d_under}[-_]?", "", raw_val, flags=re.IGNORECASE)
                            raw_val = raw_val.strip()

                            if raw_val.lower() in {"tumor", "primary"}:
                                tp = "Baseline (Tumor)"
                                tissue_origin = "Tumor"
                                tissue_type = "Primary Tumor"
                            elif re.match(r"^T\d+$", raw_val, re.IGNORECASE) or raw_val.lower() in {"baseline", "pre-op", "pre", "plasma"}:
                                tissue_origin = "Plasma"
                                tissue_type = "Plasma ctDNA"
                                tp = "T0" if raw_val.lower() in {"baseline", "pre-op", "pre"} else raw_val.upper()
                            elif raw_val.lower() in {"normal", "germline", "bld", "blood"}:
                                tp = "Normal"
                                tissue_origin = "Normal"
                                tissue_type = "Germline Normal"
                            else:
                                tp = raw_val

                            if dt_obj2:
                                date_str = date_str2
                                dt_obj = dt_obj2
                            elif third_col:
                                days_val = safe_float(third_col, default=None)

                        clean_sname = extract_clean_sample_name(sample_id, donor_id)
                        rec = {
                            "patient_id": donor_id,
                            "sample_id": sample_id,
                            "timepoint": tp,
                            "date": date_str,
                            "date_dt": dt_obj,
                            "days": days_val,
                            "tissue_type": tissue_type,
                            "tissue_origin": tissue_origin,
                            "clean_sample_name": clean_sname,
                        }
                        register_metadata_record(metadata_map, aliases, rec, extra_keys=[time_or_date])
        except Exception as e:
            print(f"Warning: Failed parsing metadata file {cf}: {e}", file=sys.stderr)

    return {
        "metadata_map": metadata_map,
        "aliases": aliases,
    }


def infer_timepoint_from_name(name, is_plasma=False):
    """
    Infer human-readable timepoint label (e.g. Baseline, T0, T1, FU-1).
    Returns empty string if no explicit timepoint token is found so clean sample name can be used.
    """
    if not is_plasma:
        return "Baseline (Tumor)"

    n_lower = name.lower()
    # TRACEOC timepoint tokens (_03_ -> T0, _04_ -> T1, _05_ -> T2, _06_ -> T3, _18_ -> T4)
    if "_03_" in name or "-03-" in name:
        return "T0 (Baseline Plasma)"
    if "_04_" in name or "-04-" in name:
        return "T1"
    if "_05_" in name or "-05-" in name:
        return "T2"
    if "_06_" in name or "-06-" in name:
        return "T3"
    if "_18_" in name or "-18-" in name:
        return "T4"

    # Explicit T0, T1, T2...
    m_t = re.search(r"\b(T[0-9]+)\b", name, re.IGNORECASE)
    if m_t:
        return m_t.group(1).upper()

    # Explicit FU-1, FU-2, FU-002...
    m_fu = re.search(r"\b(FU[-_]?[0-9]+)\b", name, re.IGNORECASE)
    if m_fu:
        return m_fu.group(1).upper()

    if "baseline" in n_lower or "pre" in n_lower:
        return "T0 (Baseline Plasma)"

    return ""


def resolve_patient_and_sample(file_path, results_path, metadata_bundle, is_plasma=False):
    """
    Intelligently extracts patient_id, sample_id, and timepoint from file path,
    directory hierarchy, metadata annotations, and standard cohort regex patterns.
    Guarantees that tumor and plasma samples from the same patient group together!
    Returns: (patient_id, sample_id, timepoint, meta_info)
    """
    metadata_map = metadata_bundle.get("metadata_map", {}) if metadata_bundle else {}

    fname = file_path.name
    stem = fname
    for ext in [".purple.purity.tsv", ".wisp.summary.tsv", ".wisp.tsv", ".orange.json", ".purple.qc"]:
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
            break

    try:
        rel_parts = file_path.relative_to(results_path).parts
    except Exception:
        rel_parts = file_path.parts

    ancestors = [p for p in rel_parts[:-1] if p.lower() not in {
        "purple", "wisp", "purity", "plasma", "orange", "linx", "somatic_annotations", "results", "output"
    }]

    # Ground-truth Patient Directory Check:
    # In standard Oncoanalyser output trees, the top-level directory directly under results_path defines patient identity!
    dir_patient_id = None
    inferred_tp_from_dir = None
    if rel_parts and len(rel_parts) > 1:
        top_folder = rel_parts[0]
        # In flat directory layouts, top_folder may be a sample directory like "HCC-B-007-T0", "HCC-B-007-T1", or "HCC-B-007"
        # Check if top_folder has a trailing timepoint or sample token suffix
        m_tp_suf = re.search(r"[-_](t[0-9]+|fu[0-9]+|baseline|plasma|cfdna|normal|tumor|tumour|primary|bld|blood|germline)$", top_folder, re.IGNORECASE)
        clean_top_folder = top_folder
        if m_tp_suf:
            clean_top_folder = top_folder[:m_tp_suf.start()]
            suf_token = m_tp_suf.group(1).upper()
            if re.match(r"^T[0-9]+$", suf_token):
                inferred_tp_from_dir = suf_token
            elif re.match(r"^FU[0-9]+$", suf_token):
                inferred_tp_from_dir = suf_token
            elif suf_token in {"BASELINE", "PRE", "PRE-OP", "PREOP"}:
                inferred_tp_from_dir = "T0" if is_plasma else "Baseline (Tumor)"
            elif suf_token in {"TUMOR", "TUMOUR", "PRIMARY"}:
                inferred_tp_from_dir = "Baseline (Tumor)"

        cohort_patterns = [
            (r"HCC[_-]B[_-](\d+)", lambda m: f"HCC-B-{int(m.group(1)):03d}"),
            (r"BFIERCE[_-](\d+)", lambda m: f"BFIERCE_{int(m.group(1)):04d}"),
            (r"TRACEOC[_-](\d+)", lambda m: f"OV-P-{int(m.group(1)):02d}"),
            (r"OV[_-]P[_-](\d+)", lambda m: f"OV-P-{int(m.group(1)):02d}"),
            (r"DYNA[_-]P(\d+)", lambda m: f"DYNA-P{int(m.group(1)):02d}"),
            (r"MOHCCN[_-]([A-Za-z0-9]+)", lambda m: f"MOHCCN_{m.group(1)}"),
            (r"\b(P|PT|PAT|DONOR)[_-]?(\d{2,5})\b", lambda m: f"{m.group(1).upper()}-{int(m.group(2)):02d}"),
        ]
        for pat, fmt in cohort_patterns:
            m = re.search(pat, clean_top_folder, re.IGNORECASE)
            if m:
                dir_patient_id = fmt(m)
                break
        if not dir_patient_id and clean_top_folder.lower() not in {"tmp", "temp", "work", "sample_sheets", "metadata", "logs"}:
            if not re.search(r"(_wisp|_0[1-9]_|[-_]tumor|[-_]plasma|[-_]normal)", clean_top_folder, re.IGNORECASE):
                if re.match(r"^[A-Za-z0-9_-]+$", clean_top_folder) and len(clean_top_folder) >= 3:
                    dir_patient_id = clean_top_folder

    # 1. Direct metadata lookup by patient-qualified keys, stem, file name, or ancestors
    search_items = []
    if dir_patient_id:
        search_items.extend([f"{dir_patient_id}_{stem}", f"{dir_patient_id}-{stem}"])
    search_items.extend([stem, fname] + ancestors)
    cand_keys = []
    for item in search_items:
        if not item:
            continue
        cand_keys.append(item)
        cand_keys.append(item.lower())

        item_no_wisp = re.sub(r"[-_]wisp$", "", item, flags=re.IGNORECASE)
        cand_keys.append(item_no_wisp)
        cand_keys.append(item_no_wisp.lower())

        item_no_pid = re.sub(r"^(OV[-_]P[-_]\d+|TRACEOC[-_]\d+|BFIERCE[-_]\d+|DYNA[-_]P\d+|MOHCCN[-_][A-Za-z0-9]+)[-_]", "", item_no_wisp, flags=re.IGNORECASE)
        cand_keys.append(item_no_pid)
        cand_keys.append(item_no_pid.lower())

        item_no_ali = re.sub(r"[-_]\d{1,2}$", "", item_no_pid)
        cand_keys.append(item_no_ali)
        cand_keys.append(item_no_ali.lower())

        clean_tok = extract_clean_sample_name(item)
        if clean_tok:
            cand_keys.append(clean_tok)
            cand_keys.append(clean_tok.lower())

    found_m = None
    for cand_key in cand_keys:
        if cand_key in metadata_map:
            m_rec = metadata_map[cand_key]
            # Guard against cross-patient contamination on generic keys (e.g. A1, 1A)
            meta_pid = m_rec.get("patient_id")
            if dir_patient_id and meta_pid and meta_pid != dir_patient_id:
                if not (dir_patient_id.startswith(meta_pid) or meta_pid.startswith(dir_patient_id) or dir_patient_id.replace("_", "-") == meta_pid.replace("_", "-")):
                    continue
            found_m = m_rec
            break

    patient_id = dir_patient_id
    sample_id = stem
    if is_plasma and sample_id.endswith("_wisp"):
        sample_id = sample_id[:-5]

    tp = None
    meta_date = ""
    meta_date_dt = None
    meta_days = None

    if found_m:
        if found_m.get("patient_id"):
            patient_id = found_m["patient_id"]
        elif not patient_id:
            patient_id = dir_patient_id
        sample_id = found_m.get("sample_id") or sample_id
        tp = found_m.get("timepoint")
        meta_date = found_m.get("date", "")
        meta_date_dt = found_m.get("date_dt")
        meta_days = found_m.get("days")

    # 2. Extract Patient Token via Known Cohort Patterns if not from directory or metadata:
    if not patient_id:
        search_strings = ancestors + [stem]
        for s in search_strings:
            m_tr = re.search(r"TRACEOC[_-](\d+)", s, re.IGNORECASE)
            if m_tr:
                num = int(m_tr.group(1))
                patient_id = f"OV-P-{num:02d}"
                break
            m_ov = re.search(r"OV[_-]P[_-](\d+)", s, re.IGNORECASE)
            if m_ov:
                patient_id = f"OV-P-{int(m_ov.group(1)):02d}"
                break
            m_bf = re.search(r"BFIERCE[_-](\d+)", s, re.IGNORECASE)
            if m_bf:
                patient_id = f"BFIERCE_{int(m_bf.group(1)):04d}"
                break
            m_dy = re.search(r"DYNA[_-]P(\d+)", s, re.IGNORECASE)
            if m_dy:
                patient_id = f"DYNA-P{int(m_dy.group(1)):02d}"
                break
            m_moh = re.search(r"MOHCCN[_-]([A-Za-z0-9]+)", s, re.IGNORECASE)
            if m_moh:
                patient_id = f"MOHCCN_{m_moh.group(1)}"
                break
            m_gen = re.search(r"\b(P|PT|PAT|DONOR)[_-]?(\d{2,5})\b", s, re.IGNORECASE)
            if m_gen:
                patient_id = f"{m_gen.group(1).upper()}-{int(m_gen.group(2)):02d}"
                break

    # 3. Check if ancestor directory is a pure patient folder
    if not patient_id and ancestors:
        top_ancestor = ancestors[0]
        if not re.search(r"(_wisp|_0[1-9]_|[-_]tumor|[-_]plasma|[-_]normal)", top_ancestor, re.IGNORECASE):
            patient_id = top_ancestor

    # 4. Fallback: Suffix Stripping on stem
    if not patient_id:
        cleaned = stem
        for suf in [r"_wisp$", r"-wisp$", r"_0[1-9]_.*$", r"[-_](tumor|tumour|primary|plasma|cfdna|normal|blood|germline).*$", r"[-_](t[0-9]+|fu[0-9]+|baseline).*$", r"[-_]LB[0-9]+[-_][0-9]+$"]:
            cleaned = re.sub(suf, "", cleaned, flags=re.IGNORECASE)
        patient_id = cleaned if cleaned else (ancestors[0] if ancestors else stem)

    clean_sample_name = extract_clean_sample_name(sample_id, patient_id)

    # Resolve timepoint label
    if not is_plasma:
        timepoint = "Baseline (Tumor)"
    else:
        if tp and tp not in {"Plasma ctDNA", "Plasma"}:
            timepoint = tp
        elif meta_date:
            timepoint = meta_date
        elif meta_days is not None:
            timepoint = f"Day {int(meta_days)}"
        elif inferred_tp_from_dir:
            timepoint = inferred_tp_from_dir
        else:
            inferred = infer_timepoint_from_name(sample_id, is_plasma=True)
            if inferred and inferred not in {"Plasma ctDNA", "Plasma"}:
                timepoint = inferred
            else:
                timepoint = clean_sample_name or sample_id

    meta_info = {
        "patient_id": patient_id,
        "sample_id": sample_id,
        "timepoint": timepoint,
        "date": meta_date,
        "date_dt": meta_date_dt,
        "days": meta_days,
        "clean_sample_name": clean_sample_name,
    }

    return patient_id, sample_id, timepoint, meta_info


def is_plasma_sample(sample_name, file_path=None, metadata_bundle=None):
    """Detect whether a sample is a longitudinal plasma/cfDNA sample rather than a primary tumor."""
    s = str(sample_name)
    s_lower = s.lower()

    # 1. Metadata check
    if metadata_bundle and "metadata_map" in metadata_bundle:
        m_map = metadata_bundle["metadata_map"]
        cand_keys = [s, s_lower]
        if file_path:
            cand_keys.append(Path(file_path).name)
            cand_keys.append(Path(file_path).stem)
        for k in cand_keys:
            if k in m_map:
                tt = (m_map[k].get("tissue_type") or m_map[k].get("tissue_origin") or m_map[k].get("tissue origin") or "").lower()
                if "plasma" in tt or "cfdna" in tt or tt == "pl":
                    return True
                if "primary" in tt or "tumor" in tt:
                    return False

    # 2. TRACEOC sample numbering: _01_ is primary tumor, _02_ is normal, _03_+ are plasma timepoints
    m_tr = re.search(r"TRACEOC_\d+_(\d+)_", s, re.IGNORECASE)
    if m_tr:
        num = int(m_tr.group(1))
        if num >= 3:
            return True
        elif num == 1:
            return False

    # 3. Keyword / token checks
    if any(tok in s_lower for tok in ["plasma", "cfdna", "pl_p"]):
        return True
    if re.search(r"[-_](t[0-9]+|fu[0-9]+)$", s_lower):
        return True
    if re.search(r"[-_]r$", s_lower) or re.search(r"[-_](resection|tumor|tumour|prim|primary)\b", s_lower):
        return False
    if re.search(r"[-_](t[0-9]+|fu[0-9]+)\b", s_lower) and not re.search(r"[-_]r\b", s_lower):
        return True

    return False


def timepoint_sort_key(pl):
    """
    Sort plasma records chronologically by:
    1. Calendar date (if available): parsed datetime timestamp
    2. Elapsed days from baseline (if available): numerical days
    3. Explicit timepoint tokens: Baseline/T0/pre-op -> T1 -> T2 -> ...
    4. Natural numerical sort of clean sample name / sample ID
    """
    sid = str(pl.get("sample_id") or "").strip()
    clean_sid = pl.get("clean_sample_name") or extract_clean_sample_name(sid, pl.get("patient_id"))

    # 1. Calendar date
    dt_obj = pl.get("date_dt")
    if dt_obj and isinstance(dt_obj, datetime):
        return (0, dt_obj.timestamp(), natural_sort_key(clean_sid or sid))
    date_str = pl.get("date") or pl.get("sampling_date")
    if date_str:
        dt_parsed, _ = parse_date_safe(date_str)
        if dt_parsed:
            return (0, dt_parsed.timestamp(), natural_sort_key(clean_sid or sid))

    # 2. Elapsed days from baseline
    days_val = pl.get("days")
    if days_val is not None:
        try:
            return (1, float(days_val), natural_sort_key(clean_sid or sid))
        except (ValueError, TypeError):
            pass

    # 3. Explicit timepoint token matching
    tp = str(pl.get("timepoint") or "").strip().lower()
    combined = f"{tp} {sid.lower()}"

    # Tumor / baseline check
    if "01" in sid.lower() and any(w in sid.lower() for w in ["tumor", "primary"]):
        return (2, 0, natural_sort_key(clean_sid or sid))
    if any(tok in combined for tok in ["baseline", "pre-op", "preop", "pre", "t0", "t00", "day 0", "day0"]):
        return (2, 1, natural_sort_key(clean_sid or sid))

    m_t = re.search(r"\bt(\d+)\b", tp)
    if m_t:
        return (2, 10 + int(m_t.group(1)), natural_sort_key(clean_sid or sid))

    m_fu = re.search(r"\bfu[-_]?(\d+)\b", tp)
    if m_fu:
        return (2, 100 + int(m_fu.group(1)), natural_sort_key(clean_sid or sid))

    m_day = re.search(r"\b(?:day|d|cycle|c)[-_]?(\d+)\b", tp)
    if m_day:
        return (2, 1000 + int(m_day.group(1)), natural_sort_key(clean_sid or sid))

    # 4. Fallback: Natural numerical order of clean sample name / sample ID
    return (3, 0, natural_sort_key(clean_sid or sid))


def load_optional_annotations(sample_sheets_dir, metadata_tsv, timepoints_tsv=None):
    """Compatibility wrapper for loading timepoints."""
    bundle = load_cohort_metadata(Path.cwd(), sample_sheets_dir, metadata_tsv, timepoints_tsv)
    return {k: v.get("timepoint", "") for k, v in bundle.get("metadata_map", {}).items() if v.get("timepoint")}


# ==============================================================================
# 3. MASTER COHORT COMPILATION
# ==============================================================================

def load_purity_overrides(overrides_path):
    """
    Loads curated purity and ploidy overrides from purity_overrides.tsv.
    Supports keyed lookup by (patient_id, tumor_id), patient_id, or tumor_id.
    """
    overrides = {}
    if not overrides_path or not os.path.isfile(overrides_path):
        return overrides
    try:
        with open(overrides_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                pid = (r.get("patient_id") or "").strip()
                tid = (r.get("tumor_id") or "").strip()
                pur_val = safe_float(r.get("curated_purity"))
                pl_val = safe_float(r.get("curated_ploidy"))
                meth = (r.get("curation_method") or "").strip()
                notes = (r.get("notes") or "").strip()

                if pur_val is not None and pur_val > 1.0:
                    pur_val = pur_val / 100.0

                rec = {
                    "patient_id": pid,
                    "tumor_id": tid,
                    "purity": pur_val,
                    "ploidy": pl_val,
                    "method": meth,
                    "notes": notes
                }
                if pid and tid:
                    overrides[(pid, tid)] = rec
                if pid:
                    overrides[pid] = rec
                if tid:
                    overrides[tid] = rec
        print(f"Loaded {len(overrides)} purity/ploidy overrides from {overrides_path}")
    except Exception as e:
        print(f"Warning: Failed to load purity overrides from {overrides_path}: {e}", file=sys.stderr)
    return overrides


def compile_cohort(results_dir, genes_of_interest, sample_sheets_dir=None, metadata_tsv=None, timepoints_tsv=None, purity_overrides_tsv=None, output_dir=None):
    """
    Main compilation routine:
    Uses global recursive discovery across results_dir to discover all PURPLE purity/ploidy,
    Orange clinical/genomic metrics, somatic driver catalogs, CNVs, and longitudinal WISP ctDNA results.
    Intelligently associates tumor and plasma samples to the same patient regardless of flat or nested
    directory layouts.
    """
    results_path = Path(results_dir).resolve()
    out_path = Path(output_dir).resolve() if output_dir else (results_path / "wgts_summary")
    metadata_bundle = load_cohort_metadata(results_path, sample_sheets_dir, metadata_tsv, timepoints_tsv)
    discovered = discover_cohort_files(results_path)
    purity_overrides = load_purity_overrides(purity_overrides_tsv)

    purple_files = discovered["purple_files"]
    wisp_files = discovered["wisp_files"]
    orange_map = discovered["orange_map"]
    global_ichor_map = discovered["ichor_map"]
    global_ichor_pdf_map = discovered.get("ichor_pdf_map", {})

    cohort = {
        "cohort_name": results_path.name,
        "patients": {},
        "all_tumors": [],
        "all_plasmas": [],
        "all_normals": [],
        "gene_frequencies": defaultdict(int),
        "total_analyzed_patients": 0,
        "has_mrd_data": False,
        "all_orange_driver_genes": set(),
        "orange_reports_b64": {},
        "ichor_reports_b64": {},
    }

    def get_patient_record(pid, p_dir_path=None):
        def _find_matching_pdir(target_id):
            cand_dirs = [
                results_path / target_id,
                results_path / target_id.replace("-", "_"),
                results_path / target_id.replace("_", "-"),
            ]
            m_pad = re.search(r'^(.*?)([-_])0*(\d+)$', target_id)
            if m_pad:
                pref, sep, num = m_pad.groups()
                for s_sep in ["-", "_"]:
                    cand_dirs.append(results_path / f"{pref}{s_sep}{num}")
                    cand_dirs.append(results_path / f"{pref}{s_sep}{int(num):04d}")
            for d in cand_dirs:
                if d.is_dir():
                    return d
            return None

        if pid not in cohort["patients"]:
            matched_d = _find_matching_pdir(pid)
            canonical_pdir = matched_d if matched_d else (Path(p_dir_path) if p_dir_path else results_path)
            cohort["patients"][pid] = {
                "patient_id": pid,
                "patient_dir": str(canonical_pdir),
                "tumors": [],
                "normals": [],
                "plasmas": [],
                "ichor": {},
            }
        else:
            if cohort["patients"][pid]["patient_dir"] == str(results_path):
                matched_d = _find_matching_pdir(pid)
                if matched_d:
                    cohort["patients"][pid]["patient_dir"] = str(matched_d)
        return cohort["patients"][pid]

    # ------------------------------------------------------------------
    # 1. Discover and Process Primary Tumor & Purple Results
    # ------------------------------------------------------------------
    plasma_purple_files = {}
    for pf in purple_files:
        raw_tumor_id = pf.name.replace(".purple.purity.tsv", "")
        # Do not classify longitudinal plasma / cfDNA samples as Primary Tumors
        if is_plasma_sample(raw_tumor_id, pf, metadata_bundle):
            plasma_purple_files[raw_tumor_id] = pf
            continue

        p_stats = parse_purple_purity_tsv(pf)
        if not p_stats:
            continue

        purple_dir = pf.parent
        patient_id, tumor_id, tp, meta_info = resolve_patient_and_sample(pf, results_path, metadata_bundle, is_plasma=False)
        p_dir = (results_path / patient_id) if (results_path / patient_id).is_dir() else purple_dir.parent.parent
        p_record = get_patient_record(patient_id, p_dir)

        # Companion QC
        qc_file = purple_dir / f"{raw_tumor_id}.purple.qc"
        if not qc_file.is_file():
            qc_file = purple_dir / f"{tumor_id}.purple.qc"
        qc_data = parse_purple_qc(qc_file)

        # Companion Orange JSON
        orange_file = find_orange_json(raw_tumor_id, purple_dir, Path(p_record["patient_dir"]), orange_map, patient_id=patient_id)
        if not orange_file and tumor_id != raw_tumor_id:
            orange_file = find_orange_json(tumor_id, purple_dir, Path(p_record["patient_dir"]), orange_map, patient_id=patient_id)
        orange_data = parse_orange_json(orange_file) if orange_file else {}

        # CUPPA predicted cancer type and likelihood (with TSV fallback)
        cuppa_cancer_type = orange_data.get("cuppa_cancer_type")
        cuppa_likelihood = orange_data.get("cuppa_likelihood")
        cuppa_predictions = orange_data.get("cuppa_predictions") or []

        if not cuppa_cancer_type:
            cuppa_tsv = find_cuppa_tsv(raw_tumor_id, purple_dir, Path(p_record["patient_dir"]), patient_id=patient_id)
            if not cuppa_tsv and tumor_id != raw_tumor_id:
                cuppa_tsv = find_cuppa_tsv(tumor_id, purple_dir, Path(p_record["patient_dir"]), patient_id=patient_id)
            if cuppa_tsv:
                c_info = parse_cuppa_pred_summ(cuppa_tsv)
                cuppa_cancer_type = c_info.get("cancer_type", "")
                cuppa_likelihood = c_info.get("likelihood")
                cuppa_predictions = c_info.get("predictions", [])

        # TMB & MSI
        tmb = orange_data.get("tmb")
        if tmb is None:
            tmb = safe_float(qc_data.get("TMB"))

        msi_status = orange_data.get("msi_status") or p_stats.get("ms_status") or qc_data.get("MSI", "MSS")
        if not msi_status or str(msi_status).upper() in {"NA", "UNKNOWN", ""}:
            msi_status = "MSS"

        # Resolve coverage metrics from bamtools, orange_data, or amber QC fallback
        t_ids = {raw_tumor_id, tumor_id, orange_data.get("sample_id")}
        r_ids = {orange_data.get("ref_sample_id"), meta_info.get("normal_id")}
        bam_cov = resolve_bamtools_coverage(purple_dir.parent, patient_dir=p_dir, tumor_ids=t_ids, ref_ids=r_ids)

        tumor_mean_cov = bam_cov["tumor"].get("mean_coverage") or orange_data.get("tumor_mean_coverage")
        if tumor_mean_cov is None:
            tumor_mean_cov = safe_float(qc_data.get("AmberMeanDepth")) or orange_data.get("amber_mean_depth")

        tumor_median_cov = bam_cov["tumor"].get("median_coverage") or orange_data.get("tumor_median_coverage")
        normal_mean_cov = bam_cov["normal"].get("mean_coverage") or orange_data.get("normal_mean_coverage")
        normal_median_cov = bam_cov["normal"].get("median_coverage") or orange_data.get("normal_median_coverage")

        wgd = orange_data.get("wgd")
        if wgd is None:
            wgd = p_stats.get("ploidy", 2.0) >= 2.8

        qc_status = orange_data.get("qc_status") or qc_data.get("QCStatus", "PASS")

        # Parse Somatic Driver Catalog
        catalog_file = purple_dir / f"{raw_tumor_id}.purple.driver.catalog.somatic.tsv"
        if not catalog_file.is_file():
            catalog_file = purple_dir / f"{tumor_id}.purple.driver.catalog.somatic.tsv"
        drivers = [d for d in parse_somatic_driver_catalog(catalog_file) if d.get("gene") not in PHARMACOGENOMIC_GENES] if catalog_file.is_file() else []

        # Parse normal read depth from somatic VCF to annotate low-normal coverage status
        normal_depth_map = parse_somatic_vcf_normal_depth(purple_dir, raw_tumor_id, tumor_id)
        for d in drivers:
            c = str(d.get("chromosome", "")).replace("chr", "")
            p = d.get("start") or d.get("position")
            g = d.get("gene", "")
            nd = normal_depth_map.get((c, p)) or (normal_depth_map.get((c, p - 1)) if p else None) or (normal_depth_map.get((c, p + 1)) if p else None) or normal_depth_map.get((g, p)) or normal_depth_map.get(g)
            if nd:
                d["normal_dp"] = nd.get("normal_dp")
                d["normal_ad"] = nd.get("normal_ad")
                if nd.get("normal_dp") is not None and nd["normal_dp"] < 10:
                    d["status_flag"] = "somatic/germline status unknown due to low normal coverage"
                    d["is_status_unknown"] = True

        # Enrich and merge somatic mutations from Orange JSON
        for ov in orange_data.get("somatic_variants", []):
            g = ov["gene"]
            if g in PHARMACOGENOMIC_GENES:
                continue
            c = str(ov.get("chromosome", "")).replace("chr", "")
            p = ov.get("position")
            nd = normal_depth_map.get((c, p)) or (normal_depth_map.get((c, p - 1)) if p else None) or (normal_depth_map.get((c, p + 1)) if p else None) or normal_depth_map.get((g, p)) or normal_depth_map.get(g)
            if nd:
                ov["normal_dp"] = nd.get("normal_dp")
                ov["normal_ad"] = nd.get("normal_ad")
                if nd.get("normal_dp") is not None and nd["normal_dp"] < 10:
                    ov["status_flag"] = "somatic/germline status unknown due to low normal coverage"
                    ov["is_status_unknown"] = True
            matched = False
            for d in drivers:
                if d["gene"] == g:
                    if not d.get("hgvs_protein") or d.get("hgvs_protein") == ov.get("hgvs_protein"):
                        matched = True
                        if not d.get("hgvs_protein") and ov.get("hgvs_protein"):
                            d["hgvs_protein"] = ov["hgvs_protein"]
                        if not d.get("hgvs_coding") and ov.get("hgvs_coding"):
                            d["hgvs_coding"] = ov["hgvs_coding"]
                        if (d.get("vaf") is None or d.get("vaf") == 0) and ov.get("vaf") is not None:
                            d["vaf"] = ov["vaf"]
                        if ov.get("alt_type") and d.get("alt_type") == "Mutation":
                            d["alt_type"] = ov["alt_type"]
                            d["subtype"] = ov.get("subtype", ov["alt_type"])
                        if ov.get("is_status_unknown"):
                            d["is_status_unknown"] = True
                            d["status_flag"] = ov.get("status_flag")
                            d["normal_dp"] = ov.get("normal_dp")
                            d["normal_ad"] = ov.get("normal_ad")
                        break
            if not matched:
                drivers.append(ov)

        for sd in orange_data.get("somatic_drivers", []):
            g = (sd.get("gene") or "").strip().upper()
            if g and g not in PHARMACOGENOMIC_GENES and not any(d["gene"] == g for d in drivers):
                drv_t = sd.get("type", "MUTATION").upper()
                if "AMP" in drv_t:
                    alt_t, sub_t = "Amplification", "Amplification"
                elif "DEL" in drv_t:
                    alt_t, sub_t = "Deep Deletion", "Deep Deletion"
                else:
                    alt_t, sub_t = "Mutation", "Hotspot Driver" if sd.get("likelihoodMethod") == "HOTSPOT" else "Driver"
                drivers.append({
                    "gene": g,
                    "alt_type": alt_t,
                    "subtype": sub_t,
                    "driver": drv_t,
                    "likelihood": safe_float(sd.get("driverLikelihood", 1.0)),
                    "vaf": None,
                    "hgvs_protein": "",
                    "hgvs_coding": "",
                    "transcript": sd.get("transcript", ""),
                    "chromosome": "",
                    "start": None,
                    "is_germline": False,
                })

        # Parse Germline Driver Catalog
        germline_file = purple_dir / f"{raw_tumor_id}.purple.driver.catalog.germline.tsv"
        if not germline_file.is_file():
            germline_file = purple_dir / f"{tumor_id}.purple.driver.catalog.germline.tsv"
        germline_drivers = [d for d in parse_germline_driver_catalog(germline_file) if d.get("gene") not in PHARMACOGENOMIC_GENES] if germline_file.is_file() else []

        for ogv in orange_data.get("germline_variants", []):
            g = ogv["gene"]
            if g in PHARMACOGENOMIC_GENES:
                continue
            matched = False
            for gd in germline_drivers:
                if gd["gene"] == g:
                    matched = True
                    if not gd.get("hgvs_protein") and ogv.get("hgvs_protein"):
                        gd["hgvs_protein"] = ogv["hgvs_protein"]
                    if not gd.get("hgvs_coding") and ogv.get("hgvs_coding"):
                        gd["hgvs_coding"] = ogv["hgvs_coding"]
                    if (gd.get("vaf") is None or gd.get("vaf") == 0) and ogv.get("vaf") is not None:
                        gd["vaf"] = ogv["vaf"]
            if not matched:
                germline_drivers.append(ogv)

        # LINX Driver Catalog (somatic disruptions & SV drivers on disk)
        linx_cat_file = purple_dir.parent / "linx" / "somatic_annotations" / f"{raw_tumor_id}.linx.driver.catalog.tsv"
        if not linx_cat_file.is_file():
            linx_cat_file = purple_dir.parent / "linx" / "somatic_annotations" / f"{tumor_id}.linx.driver.catalog.tsv"
        if not linx_cat_file.is_file():
            linx_cat_file = purple_dir.parent / "linx" / f"{raw_tumor_id}.linx.driver.catalog.tsv"
        if not linx_cat_file.is_file():
            linx_cat_file = purple_dir.parent / "linx" / f"{tumor_id}.linx.driver.catalog.tsv"
        if linx_cat_file.is_file():
            for ld in parse_linx_driver_catalog(linx_cat_file):
                g = ld["gene"]
                if g and g not in PHARMACOGENOMIC_GENES:
                    if not any(d["gene"] == g and d.get("alt_type") == ld["alt_type"] for d in drivers):
                        drivers.append(ld)
                        cohort["all_orange_driver_genes"].add(g)

        # LINX disruptions from Orange JSON
        for dis in orange_data.get("disruptions", []):
            g = (dis.get("gene") or "").strip().upper()
            if g and g not in PHARMACOGENOMIC_GENES:
                if not any(d["gene"] == g and d.get("alt_type") == "Structural Variant" for d in drivers):
                    drivers.append({
                        "gene": g,
                        "alt_type": "Structural Variant",
                        "subtype": "Gene Disruption",
                        "driver": "DISRUPTION",
                        "likelihood": 1.0,
                        "vaf": None,
                        "hgvs_protein": f"Disruption ({dis.get('type', 'SV')})",
                        "hgvs_coding": "",
                        "transcript": dis.get("transcript", ""),
                        "chromosome": dis.get("chromosome", ""),
                        "start": dis.get("position"),
                        "is_germline": dis.get("is_germline", False),
                    })
                    cohort["all_orange_driver_genes"].add(g)

        # LINX Fusions
        linx_fusion_file = purple_dir.parent / "linx" / "somatic_annotations" / f"{raw_tumor_id}.linx.fusion.tsv"
        if not linx_fusion_file.is_file():
            linx_fusion_file = purple_dir.parent / "linx" / "somatic_annotations" / f"{tumor_id}.linx.fusion.tsv"
        if not linx_fusion_file.is_file():
            linx_fusion_file = purple_dir.parent / "linx" / f"{raw_tumor_id}.linx.fusion.tsv"
        if not linx_fusion_file.is_file():
            linx_fusion_file = purple_dir.parent / "linx" / f"{tumor_id}.linx.fusion.tsv"
        fusions = parse_linx_fusions(linx_fusion_file) if linx_fusion_file.is_file() else {}

        for ofus in orange_data.get("fusions", []):
            gs = (ofus.get("gene_start") or ofus.get("gene_up") or ofus.get("geneStart") or "").strip().upper()
            ge = (ofus.get("gene_end") or ofus.get("gene_down") or ofus.get("geneEnd") or "").strip().upper()
            fname = ofus.get("name") or (f"{gs}_{ge}" if gs and ge else (gs or ge))
            f_entry = {
                "name": fname,
                "type": ofus.get("type") or ofus.get("reportedType", "FUSION"),
                "driver_interpretation": ofus.get("driver_interpretation", "HIGH"),
                "gene_start": gs,
                "gene_end": ge,
            }
            if gs and gs not in PHARMACOGENOMIC_GENES:
                if gs not in fusions: fusions[gs] = []
                if not any(x["name"] == fname for x in fusions[gs]): fusions[gs].append(f_entry)
                cohort["all_orange_driver_genes"].add(gs)
            if ge and ge not in PHARMACOGENOMIC_GENES:
                if ge not in fusions: fusions[ge] = []
                if not any(x["name"] == fname for x in fusions[ge]): fusions[ge].append(f_entry)
                cohort["all_orange_driver_genes"].add(ge)

        for og in orange_data.get("all_driver_genes", []):
            if og not in PHARMACOGENOMIC_GENES:
                cohort["all_orange_driver_genes"].add(og)

        # Parse LINX SVs to annotate breakends with partners, coordinates, and junction copy number
        linx_svs_file = purple_dir.parent / "linx" / "somatic_annotations" / f"{raw_tumor_id}.linx.svs.tsv"
        if not linx_svs_file.is_file():
            linx_svs_file = purple_dir.parent / "linx" / "somatic_annotations" / f"{tumor_id}.linx.svs.tsv"
        if not linx_svs_file.is_file():
            linx_svs_file = purple_dir.parent / "linx" / f"{raw_tumor_id}.linx.svs.tsv"
        if not linx_svs_file.is_file():
            linx_svs_file = purple_dir.parent / "linx" / f"{tumor_id}.linx.svs.tsv"
        linx_svs_map = parse_linx_svs(linx_svs_file) if linx_svs_file.is_file() else {}

        # Resolve detailed structural variants / disruptions for all affected genes
        structural_variants_by_gene = defaultdict(list)
        for dis in orange_data.get("disruptions", []):
            g = (dis.get("gene") or "").strip().upper()
            if not g or g in PHARMACOGENOMIC_GENES:
                continue
            sv_id = str(dis.get("sv_id") or "")
            sv_info = linx_svs_map.get(sv_id, {})
            partner_gene = sv_info.get("gene_end") if sv_info.get("gene_start") == g else sv_info.get("gene_start")
            coords_start = sv_info.get("coords_start", "")
            coords_end = sv_info.get("coords_end", "")
            pos_str = str(dis.get("position") or "")
            if pos_str and pos_str in coords_start:
                bk_coords = coords_start
                partner_coords = coords_end
            elif pos_str and pos_str in coords_end:
                bk_coords = coords_end
                partner_coords = coords_start
            else:
                ori = f":{dis['orientation']}" if dis.get("orientation") is not None else ""
                bk_coords = f"{dis.get('chromosome')}:{dis.get('position')}{ori}" if dis.get("position") else ""
                partner_coords = coords_end or coords_start

            exon_str = ""
            e_up = dis.get("exon_up")
            e_down = dis.get("exon_down")
            if e_up is not None and e_down is not None:
                exon_str = f"Exons {e_up}-{e_down}" if e_up != e_down else f"Exon {e_up}"
            elif e_up is not None:
                exon_str = f"Exon {e_up}"

            rtype = dis.get("region_type", "")
            jcn = dis.get("junction_copy_number")
            if jcn is None and sv_info:
                jcn = sv_info.get("jcn_min")

            sv_name = f"{g} Disruption ({dis.get('type', 'SV')})"
            if partner_gene and partner_gene != "NONE":
                sv_name = f"{g}::{partner_gene} ({dis.get('type', 'SV')})"

            sv_entry = {
                "gene": g,
                "type": dis.get("type", "BND"),
                "event": dis.get("driver_type", "DISRUPTION"),
                "name": sv_name,
                "breakpoint": bk_coords,
                "partner": partner_gene or "Intergenic",
                "partner_coords": partner_coords,
                "exons": exon_str,
                "region_type": rtype,
                "jcn": jcn,
                "is_germline": dis.get("is_germline", False),
                "cluster_id": dis.get("cluster_id") or sv_info.get("cluster_id", ""),
                "plot_filename": dis.get("plot_filename", ""),
            }
            if not any(x.get("breakpoint") == bk_coords and x.get("type") == sv_entry["type"] for x in structural_variants_by_gene[g]):
                structural_variants_by_gene[g].append(sv_entry)

        for fg, flist in fusions.items():
            for f_item in flist:
                f_name = f_item.get("name", f"Fusion ({fg})")
                structural_variants_by_gene[fg].append({
                    "gene": fg,
                    "type": f_item.get("type", "FUSION"),
                    "event": "Gene Fusion",
                    "name": f_name,
                    "breakpoint": f_item.get("five_prime_coords", ""),
                    "partner": f_item.get("gene_3p") if f_item.get("gene_5p") == fg else f_item.get("gene_5p"),
                    "partner_coords": f_item.get("three_prime_coords", ""),
                    "exons": f"Fused Exons {f_item.get('fused_exon_up', '')} / {f_item.get('fused_exon_down', '')}".strip(" /"),
                    "region_type": "FUSION",
                    "jcn": None,
                    "is_germline": False,
                    "cluster_id": "",
                    "plot_filename": "",
                })

        cnv_file = purple_dir / f"{raw_tumor_id}.purple.cnv.gene.tsv"
        if not cnv_file.is_file():
            cnv_file = purple_dir / f"{tumor_id}.purple.cnv.gene.tsv"
        ploidy_val = p_stats.get("ploidy", 2.0)
        cnvs = parse_cnv_gene_tsv(cnv_file, ploidy=ploidy_val, amp_ratio=3.0, gender=p_stats.get("gender", ""), genes_of_interest=genes_of_interest)

        user_gset = {g.strip().upper() for g in (genes_of_interest or []) if g.strip()}
        altered_genes_in_tumor = set()
        for d in drivers:
            if d.get("gene") and d["gene"] not in PHARMACOGENOMIC_GENES:
                altered_genes_in_tumor.add(d["gene"])
        for gd in germline_drivers:
            if gd.get("gene") and gd["gene"] not in PHARMACOGENOMIC_GENES:
                altered_genes_in_tumor.add(gd["gene"])
        for g, cnv_info in cnvs.items():
            if g not in PHARMACOGENOMIC_GENES and cnv_info.get("is_driver") and cnv_info["state"] in {"Amplification", "Deep Deletion"}:
                altered_genes_in_tumor.add(g)
        for fg in fusions.keys():
            if fg not in PHARMACOGENOMIC_GENES and (is_driver_cna_candidate(fg, user_gset) or fg in cohort["all_orange_driver_genes"]):
                altered_genes_in_tumor.add(fg)
        for svg in structural_variants_by_gene.keys():
            altered_genes_in_tumor.add(svg)

        for g in altered_genes_in_tumor:
            cohort["gene_frequencies"][g] += 1

        ichor_tf = (
            discovered.get("ichor_tumor_map", {}).get((patient_id, raw_tumor_id))
            or discovered.get("ichor_tumor_map", {}).get((patient_id, tumor_id))
            or global_ichor_map.get((patient_id, raw_tumor_id))
            or global_ichor_map.get((patient_id, tumor_id))
        )
        if ichor_tf is None and p_record.get("patient_dir"):
            if "ichor_tumor" not in p_record:
                p_record["ichor_tumor"] = parse_ichor_metrics(p_record["patient_dir"], is_plasma=False)
            local_ichor = p_record["ichor_tumor"]
            ichor_tf = (
                local_ichor.get(raw_tumor_id)
                or local_ichor.get(tumor_id)
                or local_ichor.get(f"{patient_id}_{raw_tumor_id}")
                or local_ichor.get(f"{patient_id}_{tumor_id}")
            )
        if ichor_tf is None:
            ichor_tf = (
                discovered.get("ichor_tumor_map", {}).get(raw_tumor_id)
                or discovered.get("ichor_tumor_map", {}).get(tumor_id)
                or global_ichor_map.get(raw_tumor_id)
                or global_ichor_map.get(tumor_id)
            )

        # Check for curated purity / ploidy override
        ov = None
        if purity_overrides:
            ov = (
                purity_overrides.get((patient_id, raw_tumor_id))
                or purity_overrides.get((patient_id, tumor_id))
                or purity_overrides.get(patient_id)
                or purity_overrides.get(raw_tumor_id)
                or purity_overrides.get(tumor_id)
            )

        if ov and ov.get("purity") is not None:
            active_purity = ov["purity"]
            active_ploidy = ov["ploidy"] if ov.get("ploidy") is not None else ploidy_val
            curated_purity = ov["purity"]
            curated_ploidy = ov["ploidy"]
            curation_method = ov.get("method", "Curated")
            curation_notes = ov.get("notes", "")
        else:
            active_purity = p_stats.get("purity", 0.0)
            active_ploidy = ploidy_val
            curated_purity = None
            curated_ploidy = None
            curation_method = ""
            curation_notes = ""

        # Discover companion native HMF Orange PDF report
        orange_pdf = None
        if orange_file:
            cand_pdf = orange_file.with_name(orange_file.name.replace(".orange.json", ".orange.pdf").replace(".json", ".pdf"))
            if cand_pdf.is_file():
                orange_pdf = cand_pdf
            elif orange_file.with_suffix(".pdf").is_file():
                orange_pdf = orange_file.with_suffix(".pdf")
        if not orange_pdf:
            cand_odirs = []
            if p_record.get("patient_dir"):
                cand_odirs.append(Path(p_record["patient_dir"]) / tumor_id / "orange")
                cand_odirs.append(Path(p_record["patient_dir"]) / "orange")
            cand_odirs.append(purple_dir.parent / "orange")
            for cod in cand_odirs:
                if cod.is_dir():
                    c_pdf = cod / f"{tumor_id}.orange.pdf"
                    if c_pdf.is_file():
                        orange_pdf = c_pdf
                        break
                    c_list = list(cod.glob("*.orange.pdf"))
                    if c_list:
                        orange_pdf = c_list[0]
                        break

        if not orange_pdf and out_path:
            cand_staged = out_path / "orange_reports" / f"{tumor_id}.orange.pdf"
            if cand_staged.is_file():
                orange_pdf = cand_staged
            elif (out_path / "orange_reports" / f"{raw_tumor_id}.orange.pdf").is_file():
                orange_pdf = out_path / "orange_reports" / f"{raw_tumor_id}.orange.pdf"

        orange_pdf_file = str(orange_pdf) if orange_pdf else ""
        orange_pdf_rel = ""
        orange_pdf_b64 = ""
        if orange_pdf and os.path.isfile(orange_pdf):
            try:
                with open(orange_pdf, "rb") as f_pdf:
                    orange_pdf_b64 = base64.b64encode(f_pdf.read()).decode("ascii")
            except Exception as e:
                print(f"Warning: Failed reading orange PDF {orange_pdf} for base64: {e}", file=sys.stderr)

        if orange_pdf_b64:
            cohort.setdefault("orange_reports_b64", {})[tumor_id] = orange_pdf_b64

        if orange_pdf and out_path:
            top_orange_dir = out_path / "orange_reports"
            try:
                top_orange_dir.mkdir(parents=True, exist_ok=True)
                dest_name = f"{tumor_id}.orange.pdf"
                target_pdf = top_orange_dir / dest_name
                if not target_pdf.exists() or target_pdf.stat().st_size != os.path.getsize(orange_pdf):
                    shutil.copy2(orange_pdf, target_pdf)
                orange_pdf_rel = f"orange_reports/{dest_name}"
            except Exception as e:
                print(f"Warning: Failed staging orange PDF {orange_pdf}: {e}", file=sys.stderr)
                orange_pdf_rel = os.path.relpath(orange_pdf, out_path)
        elif orange_pdf:
            orange_pdf_rel = str(orange_pdf)

        # Discover companion ichorCNA genomeWide PDF reports (both standard and short)
        ichor_pdfs = find_all_ichor_pdfs(tumor_id, patient_dir=Path(p_record["patient_dir"]) if p_record.get("patient_dir") else None,
                                        discovered=discovered, patient_id=patient_id, raw_sample_id=raw_tumor_id, is_plasma=False)
        std_pdf = ichor_pdfs.get("standard") or ichor_pdfs.get("primary")
        std_rel, std_b64 = stage_and_encode_pdf(std_pdf, out_path, "ichor_reports", f"{patient_id}_{tumor_id}_genomeWide.pdf") if std_pdf else ("", "")
        sh_pdf = ichor_pdfs.get("short")
        sh_rel, sh_b64 = stage_and_encode_pdf(sh_pdf, out_path, "ichor_reports", f"{patient_id}_{tumor_id}_short_genomeWide.pdf") if sh_pdf else ("", "")

        prim_b64 = std_b64 or sh_b64
        prim_rel = std_rel or sh_rel
        if prim_b64 or std_b64 or sh_b64:
            p_variants = {patient_id, patient_id.replace("-", "_"), patient_id.replace("_", "-")}
            t_variants = {tumor_id, raw_tumor_id} if raw_tumor_id else {tumor_id}
            for p_tok in p_variants:
                for t_tok in t_variants:
                    if prim_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{t_tok}"] = prim_b64
                        cohort["ichor_reports_b64"][t_tok] = prim_b64
                    if std_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{t_tok}_standard"] = std_b64
                        cohort["ichor_reports_b64"][f"{t_tok}_standard"] = std_b64
                    if sh_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{t_tok}_short"] = sh_b64
                        cohort["ichor_reports_b64"][f"{t_tok}_short"] = sh_b64

        t_record = {
            "patient_id": patient_id,
            "tumor_id": tumor_id,
            "raw_tumor_id": raw_tumor_id,
            "purity": active_purity,
            "purity_pct": active_purity * 100.0,
            "ploidy": active_ploidy,
            "curated_purity": curated_purity,
            "curated_ploidy": curated_ploidy,
            "curation_method": curation_method,
            "curation_notes": curation_notes,
            "gender": p_stats.get("gender", ""),
            "status": p_stats.get("status", "NORMAL"),
            "tmb": tmb,
            "tmb_status": orange_data.get("tmb_status", ""),
            "msi_status": msi_status,
            "msi_score": orange_data.get("msi_score"),
            "wgd": bool(wgd),
            "tumor_mean_coverage": tumor_mean_cov,
            "tumor_median_coverage": tumor_median_cov,
            "normal_mean_coverage": normal_mean_cov,
            "normal_median_coverage": normal_median_cov,
            "qc_status": qc_status,
            "contamination": orange_data.get("contamination") or safe_float(qc_data.get("Contamination")),
            "amber_mean_depth": orange_data.get("amber_mean_depth") or safe_float(qc_data.get("AmberMeanDepth")),
            "ichor_tf": ichor_tf,
            "ichor_purity_pct": (ichor_tf * 100.0) if ichor_tf is not None else None,
            "ichor_mrd": bool(ichor_tf is not None and ichor_tf >= 0.03),
            "ichor_lod_pct": 3.0,
            "ichor_pdf_file": str(ichor_pdfs.get("primary") or "") if ichor_pdfs.get("primary") else "",
            "ichor_pdf_rel": prim_rel,
            "ichor_pdf_b64": prim_b64,
            "ichor_standard_pdf_rel": std_rel,
            "ichor_standard_pdf_b64": std_b64,
            "ichor_short_pdf_rel": sh_rel,
            "ichor_short_pdf_b64": sh_b64,
            "has_ichor_pdf": bool(prim_b64 or prim_rel),
            "has_ichor_standard": bool(std_b64 or std_rel),
            "has_ichor_short": bool(sh_b64 or sh_rel),
            "drivers": drivers,
            "germline_drivers": germline_drivers,
            "fusions": fusions,
            "structural_variants": dict(structural_variants_by_gene),
            "cnvs": cnvs,
            "altered_genes": sorted(list(altered_genes_in_tumor)),
            "purple_purity_file": str(pf),
            "orange_json_file": str(orange_file) if orange_file else "",
            "orange_pdf_file": orange_pdf_file,
            "orange_pdf_rel": orange_pdf_rel,
            "has_orange_pdf": bool(orange_pdf_b64 or orange_pdf),
            "cuppa_cancer_type": cuppa_cancer_type or "",
            "cuppa_likelihood": cuppa_likelihood,
            "cuppa_predictions": cuppa_predictions,
        }
        p_record["tumors"].append(t_record)
        cohort["all_tumors"].append(t_record)

    # ------------------------------------------------------------------
    # 2. Discover and Process Longitudinal Plasma / WISP MRD Results
    # ------------------------------------------------------------------
    for wf in wisp_files:
        raw_w = parse_wisp_summary_tsv(wf)
        if not raw_w:
            continue

        patient_id, plasma_id, tp_label, meta_info = resolve_patient_and_sample(wf, results_path, metadata_bundle, is_plasma=True)
        p_record = get_patient_record(patient_id, wf.parent.parent)

        # Match to primary tumor if available for this patient
        matched_tumor_id = "Primary_Tumor"
        cand_tumor, _ = extract_metadata_from_wisp_path(wf, results_path, patient_id)
        if p_record["tumors"]:
            matched_tumor = None
            parent_name = wf.parent.parent.name
            for t in p_record["tumors"]:
                tid = t.get("tumor_id")
                if tid and (tid == cand_tumor or tid in parent_name or (cand_tumor and cand_tumor in tid)):
                    matched_tumor = tid
                    break
            if matched_tumor:
                matched_tumor_id = matched_tumor
            else:
                matched_tumor_id = p_record["tumors"][0]["tumor_id"]
        elif cand_tumor and cand_tumor != "Primary_Tumor":
            matched_tumor_id = cand_tumor

        # 1. SNV Modality
        snv_mrd = safe_bool(raw_w.get("SNV_MRD"))
        snv_purity = safe_float(raw_w.get("SNVPurity"))
        snv_purity_low = safe_float(raw_w.get("SNVPurityLow"))
        snv_purity_high = safe_float(raw_w.get("SNVPurityHigh"))
        snv_pvalue = safe_float(raw_w.get("SNVPValue"))
        snv_lod = safe_float(raw_w.get("SNVLod"))

        # 2. LOH Modality
        loh_mrd = safe_bool(raw_w.get("LOH_MRD"))
        loh_purity = safe_float(raw_w.get("LOHPurity"))
        loh_purity_low = safe_float(raw_w.get("LOHPurityLow"))
        loh_purity_high = safe_float(raw_w.get("LOHPurityHigh"))
        loh_pvalue = safe_float(raw_w.get("LOHPValue"))
        loh_lod = safe_float(raw_w.get("LOHLod"))

        # 3. Copy Number / CNV Modality
        cnv_mrd = safe_bool(raw_w.get("CopyNumber_MRD") or raw_w.get("CNV_MRD"))
        cnv_purity = safe_float(raw_w.get("CopyNumberPurity") or raw_w.get("CNVPurity"))
        cnv_purity_low = safe_float(raw_w.get("CopyNumberPurityLow") or raw_w.get("CNVPurityLow"))
        cnv_purity_high = safe_float(raw_w.get("CopyNumberPurityHigh") or raw_w.get("CNVPurityHigh"))
        cnv_pvalue = safe_float(raw_w.get("CopyNumberPValue") or raw_w.get("CNVPValue"))
        cnv_lod = safe_float(raw_w.get("CopyNumberLod") or raw_w.get("CNVLod"))

        # 4. Integrated WISP Call
        wisp_mrd = safe_bool(raw_w.get("WispMRD"))
        wisp_purity = safe_float(raw_w.get("WispPurity"))
        wisp_method = (raw_w.get("WispMethod") or "").strip()
        wisp_pvalue = safe_float(raw_w.get("WispPValue"))

        is_mrd_positive = snv_mrd or loh_mrd or cnv_mrd or wisp_mrd
        mrd_status = "MRD Positive" if is_mrd_positive else "MRD Negative"

        # Check ichor tumor fraction across all possible candidate sample tokens
        wf_stem = wf.name.replace(".wisp.summary.tsv", "").replace(".wisp.tsv", "")
        clean_pname = meta_info.get("clean_sample_name") or extract_clean_sample_name(plasma_id, patient_id)
        folder_token = None
        for idx, part in enumerate(wf.parts):
            if part.lower() == "plasma" and idx + 1 < len(wf.parts):
                folder_token = wf.parts[idx + 1]
                break

        plasma_cand_keys = [
            plasma_id,
            clean_pname,
            folder_token,
            wf_stem,
            extract_clean_sample_name(wf_stem, patient_id),
            re.sub(rf"^{patient_id}[-_]?", "", plasma_id, flags=re.IGNORECASE) if patient_id else None,
            re.sub(rf"^{patient_id}[-_]?", "", wf_stem, flags=re.IGNORECASE) if patient_id else None,
        ]
        plasma_cand_keys = [k for k in plasma_cand_keys if k]

        ichor_tf = None
        # 1. First priority for plasma ctDNA: short-fragment size-selected ichorCNA
        for k in plasma_cand_keys:
            ichor_tf = (
                discovered.get("ichor_short_map", {}).get((patient_id, k))
                or discovered.get("ichor_short_map", {}).get(k)
            )
            if ichor_tf is not None:
                break

        # 2. Second priority: general plasma ichor map
        if ichor_tf is None:
            for k in plasma_cand_keys:
                ichor_tf = (
                    discovered.get("ichor_plasma_map", {}).get((patient_id, k))
                    or discovered.get("ichor_plasma_map", {}).get(k)
                    or global_ichor_map.get((patient_id, k))
                    or global_ichor_map.get(k)
                )
                if ichor_tf is not None:
                    break

        if ichor_tf is None and p_record.get("patient_dir"):
            if "ichor_plasma" not in p_record:
                p_record["ichor_plasma"] = parse_ichor_metrics(p_record["patient_dir"], is_plasma=True)
            local_ichor = p_record["ichor_plasma"]
            for k in plasma_cand_keys:
                ichor_tf = local_ichor.get(k) or local_ichor.get(f"{patient_id}_{k}")
                if ichor_tf is not None:
                    break

        plasma_record = {
            "patient_id": patient_id,
            "tumor_id": matched_tumor_id,
            "sample_id": plasma_id,
            "sample_type": "Plasma ctDNA",
            "timepoint": tp_label,
            "date": meta_info.get("date", ""),
            "date_dt": meta_info.get("date_dt"),
            "days": meta_info.get("days"),
            "clean_sample_name": meta_info.get("clean_sample_name") or plasma_id,
            "is_plasma": True,
            "mrd_status": mrd_status,
            "is_mrd_positive": is_mrd_positive,
            "drivers": [],
            "cnvs": {},
            "altered_genes": [],
            "tmb": None,
            "tmb_status": "NA",
            # SNV
            "snv_mrd": snv_mrd,
            "snv_purity": snv_purity,
            "snv_purity_pct": snv_purity * 100.0,
            "snv_purity_low_pct": snv_purity_low * 100.0,
            "snv_purity_high_pct": snv_purity_high * 100.0,
            "snv_pvalue": snv_pvalue,
            "snv_lod": snv_lod,
            "snv_lod_pct": snv_lod * 100.0,
            # LOH
            "loh_mrd": loh_mrd,
            "loh_purity": loh_purity,
            "loh_purity_pct": loh_purity * 100.0,
            "loh_purity_low": loh_purity_low,
            "loh_purity_low_pct": (loh_purity_low * 100.0) if loh_purity_low is not None else None,
            "loh_purity_high": loh_purity_high,
            "loh_purity_high_pct": (loh_purity_high * 100.0) if loh_purity_high is not None else None,
            "loh_pvalue": loh_pvalue,
            "loh_lod": loh_lod,
            "loh_lod_pct": loh_lod * 100.0,
            # Copy Number / CNV
            "cnv_mrd": cnv_mrd,
            "cnv_purity": cnv_purity,
            "cnv_purity_pct": cnv_purity * 100.0,
            "cnv_purity_low": cnv_purity_low,
            "cnv_purity_low_pct": (cnv_purity_low * 100.0) if cnv_purity_low is not None else None,
            "cnv_purity_high": cnv_purity_high,
            "cnv_purity_high_pct": (cnv_purity_high * 100.0) if cnv_purity_high is not None else None,
            "cnv_pvalue": cnv_pvalue,
            "cnv_lod": cnv_lod,
            "cnv_lod_pct": cnv_lod * 100.0,
            # WISP Integrated
            "wisp_mrd": wisp_mrd,
            "wisp_purity": wisp_purity,
            "wisp_purity_pct": wisp_purity * 100.0,
            "wisp_method": wisp_method,
            "wisp_pvalue": wisp_pvalue,
            # Technical QC & Coverage
            "calc_variants": safe_int(raw_w.get("CalcVariants")),
            "total_variants": safe_int(raw_w.get("TotalVariants")),
            "weighted_depth": safe_float(raw_w.get("WeightedAvgDepth")),
            "plasma_coverage_depth": safe_float(raw_w.get("WeightedAvgDepth")),
            "tumor_baseline_purity": safe_float(raw_w.get("TumorPurity")),
            "tumor_baseline_ploidy": safe_float(raw_w.get("TumorPloidy")),
            "ichor_tf": ichor_tf,
            "ichor_purity_pct": (ichor_tf * 100.0) if ichor_tf is not None else None,
            "ichor_mrd": bool(ichor_tf is not None and ichor_tf >= 0.03),
            "ichor_lod_pct": 3.0,
            "wisp_file": str(wf),
        }

        # Discover companion ichorCNA genomeWide PDF reports (both short and standard)
        ichor_pdfs = find_all_ichor_pdfs(plasma_id, patient_dir=Path(p_record["patient_dir"]) if p_record.get("patient_dir") else None,
                                        discovered=discovered, patient_id=patient_id, raw_sample_id=clean_pname, is_plasma=True)
        sh_pdf = ichor_pdfs.get("short")
        sh_rel, sh_b64 = stage_and_encode_pdf(sh_pdf, out_path, "ichor_reports", f"{patient_id}_{plasma_id}_short_genomeWide.pdf") if sh_pdf else ("", "")
        std_pdf = ichor_pdfs.get("standard")
        std_rel, std_b64 = stage_and_encode_pdf(std_pdf, out_path, "ichor_reports", f"{patient_id}_{plasma_id}_standard_genomeWide.pdf") if std_pdf else ("", "")

        prim_b64 = sh_b64 or std_b64
        prim_rel = sh_rel or std_rel
        if prim_b64 or std_b64 or sh_b64:
            p_variants = {patient_id, patient_id.replace("-", "_"), patient_id.replace("_", "-")}
            pl_variants = {plasma_id, clean_pname} if clean_pname else {plasma_id}
            for p_tok in p_variants:
                for pl_tok in pl_variants:
                    if prim_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{pl_tok}"] = prim_b64
                        cohort["ichor_reports_b64"][pl_tok] = prim_b64
                    if sh_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{pl_tok}_short"] = sh_b64
                        cohort["ichor_reports_b64"][f"{pl_tok}_short"] = sh_b64
                    if std_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{pl_tok}_standard"] = std_b64
                        cohort["ichor_reports_b64"][f"{pl_tok}_standard"] = std_b64

        plasma_record["ichor_pdf_file"] = str(ichor_pdfs.get("primary") or "") if ichor_pdfs.get("primary") else ""
        plasma_record["ichor_pdf_rel"] = prim_rel
        plasma_record["ichor_pdf_b64"] = prim_b64
        plasma_record["ichor_standard_pdf_rel"] = std_rel
        plasma_record["ichor_standard_pdf_b64"] = std_b64
        plasma_record["ichor_short_pdf_rel"] = sh_rel
        plasma_record["ichor_short_pdf_b64"] = sh_b64
        plasma_record["has_ichor_pdf"] = bool(prim_b64 or prim_rel)
        plasma_record["has_ichor_standard"] = bool(std_b64 or std_rel)
        plasma_record["has_ichor_short"] = bool(sh_b64 or sh_rel)

        # Check for companion PURPLE results for this plasma sample
        cand_pf = plasma_purple_files.get(plasma_id) or plasma_purple_files.get(raw_w.get("SampleId")) or plasma_purple_files.get(wf.name.replace(".wisp.summary.tsv", "").replace(".wisp.tsv", ""))
        if not cand_pf:
            for k, pf_cand in plasma_purple_files.items():
                if k == plasma_id or k in plasma_id or plasma_id in k:
                    cand_pf = pf_cand
                    break

        if cand_pf and cand_pf.is_file():
            plasma_record["purple_purity_file"] = str(cand_pf)
            p_stats = parse_purple_purity_tsv(cand_pf)
            if p_stats:
                plasma_record["purple_purity_pct"] = p_stats.get("purity", 0.0) * 100.0
                plasma_record["purple_ploidy"] = p_stats.get("ploidy", 2.0)

            p_cat = cand_pf.parent / f"{cand_pf.name.replace('.purple.purity.tsv', '')}.purple.driver.catalog.somatic.tsv"
            if not p_cat.is_file():
                for cf in cand_pf.parent.glob("*.purple.driver.catalog.somatic.tsv"):
                    p_cat = cf
                    break
            if p_cat.is_file():
                plasma_record["purple_drivers"] = parse_somatic_driver_catalog(p_cat)

            p_cnv = cand_pf.parent / f"{cand_pf.name.replace('.purple.purity.tsv', '')}.purple.cnv.gene.tsv"
            if not p_cnv.is_file():
                for cf in cand_pf.parent.glob("*.purple.cnv.gene.tsv"):
                    p_cnv = cf
                    break
            if p_cnv.is_file():
                plasma_record["purple_cnvs"] = parse_cnv_gene_tsv(p_cnv, ploidy=p_stats.get("ploidy", 2.0) if p_stats else 2.0, amp_ratio=3.0, genes_of_interest=genes_of_interest)

        p_record["plasmas"].append(plasma_record)
        cohort["all_plasmas"].append(plasma_record)
        cohort["has_mrd_data"] = True

    # Discover any plasma samples that had PURPLE run but no WISP summary file
    for raw_plasma_id, pf in plasma_purple_files.items():
        if any(pl["sample_id"] == raw_plasma_id or pl.get("purple_purity_file") == str(pf) for pl in cohort["all_plasmas"]):
            continue

        patient_id, plasma_id, tp_label, meta_info = resolve_patient_and_sample(pf, results_path, metadata_bundle, is_plasma=True)
        p_record = get_patient_record(patient_id, pf.parent.parent)
        matched_tumor_id = p_record["tumors"][0]["tumor_id"] if p_record["tumors"] else "Primary_Tumor"

        p_stats = parse_purple_purity_tsv(pf) or {}
        p_cat = pf.parent / f"{raw_plasma_id}.purple.driver.catalog.somatic.tsv"
        purple_drivers = parse_somatic_driver_catalog(p_cat) if p_cat.is_file() else []
        p_cnv = pf.parent / f"{raw_plasma_id}.purple.cnv.gene.tsv"
        purple_cnvs = parse_cnv_gene_tsv(p_cnv, ploidy=p_stats.get("ploidy", 2.0), amp_ratio=3.0, genes_of_interest=genes_of_interest) if p_cnv.is_file() else {}

        clean_pname = meta_info.get("clean_sample_name") or extract_clean_sample_name(plasma_id, patient_id)
        folder_token = None
        for idx, part in enumerate(pf.parts):
            if part.lower() == "plasma" and idx + 1 < len(pf.parts):
                folder_token = pf.parts[idx + 1]
                break

        plasma_cand_keys = [
            plasma_id,
            raw_plasma_id,
            clean_pname,
            folder_token,
            extract_clean_sample_name(raw_plasma_id, patient_id),
            re.sub(rf"^{patient_id}[-_]?", "", plasma_id, flags=re.IGNORECASE) if patient_id else None,
            re.sub(rf"^{patient_id}[-_]?", "", raw_plasma_id, flags=re.IGNORECASE) if patient_id else None,
        ]
        plasma_cand_keys = [k for k in plasma_cand_keys if k]

        ichor_tf = None
        # 1. First priority for plasma ctDNA: short-fragment size-selected ichorCNA
        for k in plasma_cand_keys:
            ichor_tf = (
                discovered.get("ichor_short_map", {}).get((patient_id, k))
                or discovered.get("ichor_short_map", {}).get(k)
            )
            if ichor_tf is not None:
                break

        # 2. Second priority: general plasma ichor map
        if ichor_tf is None:
            for k in plasma_cand_keys:
                ichor_tf = (
                    discovered.get("ichor_plasma_map", {}).get((patient_id, k))
                    or discovered.get("ichor_plasma_map", {}).get(k)
                    or global_ichor_map.get((patient_id, k))
                    or global_ichor_map.get(k)
                )
                if ichor_tf is not None:
                    break

        if ichor_tf is None and p_record.get("patient_dir"):
            if "ichor_plasma" not in p_record:
                p_record["ichor_plasma"] = parse_ichor_metrics(p_record["patient_dir"], is_plasma=True)
            local_ichor = p_record["ichor_plasma"]
            for k in plasma_cand_keys:
                ichor_tf = local_ichor.get(k) or local_ichor.get(f"{patient_id}_{k}")
                if ichor_tf is not None:
                    break

        pl_rec = {
            "patient_id": patient_id,
            "tumor_id": matched_tumor_id,
            "sample_id": plasma_id,
            "sample_type": "Plasma ctDNA",
            "timepoint": tp_label,
            "date": meta_info.get("date", ""),
            "date_dt": meta_info.get("date_dt"),
            "days": meta_info.get("days"),
            "clean_sample_name": meta_info.get("clean_sample_name") or plasma_id,
            "is_plasma": True,
            "mrd_status": "MRD Positive" if (p_stats.get("purity", 0.0) > 0.01 or bool(purple_drivers)) else "MRD Negative",
            "is_mrd_positive": (p_stats.get("purity", 0.0) > 0.01 or bool(purple_drivers)),
            "drivers": [],
            "cnvs": {},
            "altered_genes": [],
            "tmb": None,
            "tmb_status": "NA",
            "purple_purity_file": str(pf),
            "purple_purity_pct": p_stats.get("purity", 0.0) * 100.0,
            "purple_ploidy": p_stats.get("ploidy", 2.0),
            "purple_drivers": purple_drivers,
            "purple_cnvs": purple_cnvs,
            "snv_mrd": bool(purple_drivers),
            "snv_purity": p_stats.get("purity"),
            "snv_purity_pct": p_stats.get("purity", 0.0) * 100.0,
            "snv_pvalue": None,
            "snv_lod": None,
            "snv_lod_pct": None,
            "loh_mrd": False,
            "loh_purity": None,
            "loh_purity_pct": None,
            "cnv_mrd": bool(purple_cnvs),
            "cnv_purity": None,
            "cnv_purity_pct": None,
            "wisp_mrd": (p_stats.get("purity", 0.0) > 0.01 or bool(purple_drivers)),
            "wisp_purity": p_stats.get("purity"),
            "wisp_purity_pct": p_stats.get("purity", 0.0) * 100.0,
            "wisp_method": "PURPLE",
            "wisp_pvalue": None,
            "calc_variants": 0,
            "total_variants": len(purple_drivers),
            "weighted_depth": None,
            "plasma_coverage_depth": None,
            "tumor_baseline_purity": None,
            "tumor_baseline_ploidy": None,
            "ichor_tf": ichor_tf,
            "ichor_purity_pct": (ichor_tf * 100.0) if ichor_tf is not None else None,
            "ichor_mrd": bool(ichor_tf is not None and ichor_tf >= 0.03),
            "ichor_lod_pct": 3.0,
            "wisp_file": "",
        }

        # Discover companion ichorCNA genomeWide PDF reports (both short and standard)
        ichor_pdfs = find_all_ichor_pdfs(plasma_id, patient_dir=Path(p_record["patient_dir"]) if p_record.get("patient_dir") else None,
                                        discovered=discovered, patient_id=patient_id, raw_sample_id=clean_pname, is_plasma=True)
        sh_pdf = ichor_pdfs.get("short")
        sh_rel, sh_b64 = stage_and_encode_pdf(sh_pdf, out_path, "ichor_reports", f"{patient_id}_{plasma_id}_short_genomeWide.pdf") if sh_pdf else ("", "")
        std_pdf = ichor_pdfs.get("standard")
        std_rel, std_b64 = stage_and_encode_pdf(std_pdf, out_path, "ichor_reports", f"{patient_id}_{plasma_id}_standard_genomeWide.pdf") if std_pdf else ("", "")

        prim_b64 = sh_b64 or std_b64
        prim_rel = sh_rel or std_rel
        if prim_b64 or std_b64 or sh_b64:
            p_variants = {patient_id, patient_id.replace("-", "_"), patient_id.replace("_", "-")}
            pl_variants = {plasma_id, clean_pname} if clean_pname else {plasma_id}
            for p_tok in p_variants:
                for pl_tok in pl_variants:
                    if prim_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{pl_tok}"] = prim_b64
                        cohort["ichor_reports_b64"][pl_tok] = prim_b64
                    if sh_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{pl_tok}_short"] = sh_b64
                        cohort["ichor_reports_b64"][f"{pl_tok}_short"] = sh_b64
                    if std_b64:
                        cohort.setdefault("ichor_reports_b64", {})[f"{p_tok}_{pl_tok}_standard"] = std_b64
                        cohort["ichor_reports_b64"][f"{pl_tok}_standard"] = std_b64

        pl_rec["ichor_pdf_file"] = str(ichor_pdfs.get("primary") or "") if ichor_pdfs.get("primary") else ""
        pl_rec["ichor_pdf_rel"] = prim_rel
        pl_rec["ichor_pdf_b64"] = prim_b64
        pl_rec["ichor_standard_pdf_rel"] = std_rel
        pl_rec["ichor_standard_pdf_b64"] = std_b64
        pl_rec["ichor_short_pdf_rel"] = sh_rel
        pl_rec["ichor_short_pdf_b64"] = sh_b64
        pl_rec["has_ichor_pdf"] = bool(prim_b64 or prim_rel)
        pl_rec["has_ichor_standard"] = bool(std_b64 or std_rel)
        pl_rec["has_ichor_short"] = bool(sh_b64 or sh_rel)
        p_record["plasmas"].append(pl_rec)
        cohort["all_plasmas"].append(pl_rec)
        cohort["has_mrd_data"] = True

    # ------------------------------------------------------------------
    # 3. Chronological Sorting & Post-processing
    # ------------------------------------------------------------------
    for pid, prec in list(cohort["patients"].items()):
        for pl in prec["plasmas"]:
            if pl.get("timepoint") in {"Plasma ctDNA", "Plasma", ""}:
                if pl.get("date"):
                    pl["timepoint"] = pl["date"]
                elif pl.get("days") is not None:
                    pl["timepoint"] = f"Day {int(pl['days'])}"
                else:
                    pl["timepoint"] = pl.get("clean_sample_name") or pl["sample_id"]
        prec["plasmas"].sort(key=timepoint_sort_key)

    cohort["all_plasmas"].sort(key=lambda p: (natural_sort_key(p["patient_id"]), timepoint_sort_key(p)))

    # Apply strict clinical mutation filtering to plasma samples
    filter_and_assign_plasma_mutations(cohort)

    # Also register candidate patient dirs that had only normals or ichor if any
    for cand_p in find_candidate_patient_dirs(results_path):
        if cand_p.name not in cohort["patients"]:
            local_ichor = parse_ichor_metrics(cand_p)
            if local_ichor:
                cohort["patients"][cand_p.name] = {
                    "patient_id": cand_p.name,
                    "patient_dir": str(cand_p),
                    "tumors": [],
                    "normals": [],
                    "plasmas": [],
                    "ichor": local_ichor,
                }

    cohort["total_analyzed_patients"] = len(cohort["patients"])
    cohort["all_orange_driver_genes"] = sorted(list(cohort.get("all_orange_driver_genes", set())))
    return cohort


def filter_and_assign_plasma_mutations(cohort):
    """
    Applies strict user-defined mutation filtering for plasma ctDNA:
      - Primary tumor: Only mutations called in primary (somatic drivers).
      - Plasma samples:
          1. Mutations in genes called in the primary tumor (tracked from baseline), OR
          2. Novel/de novo mutations showing up in plasma with VAF > 5% (0.05).
      - Excludes background variants not called in primary with VAF <= 5%.
      - Clumps samples by patient with primary tumor followed by chronological plasma ctDNA timepoints.
    """
    for pid, prec in cohort["patients"].items():
        primary_mutated_genes = set()
        primary_driver_map = {}
        for t in prec.get("tumors", []):
            for d in t.get("drivers", []):
                g = (d.get("gene") or "").strip().upper()
                if g and g not in PHARMACOGENOMIC_GENES:
                    primary_mutated_genes.add(g)
                    if g not in primary_driver_map:
                        primary_driver_map[g] = d

        for pl in prec.get("plasmas", []):
            raw_candidates = []

            # 1. PURPLE candidate drivers if available
            if "purple_drivers" in pl:
                for pd in pl["purple_drivers"]:
                    raw_candidates.append(pd)

            # 2. WISP candidate somatic variants if available
            wisp_file = pl.get("wisp_file")
            if wisp_file and os.path.exists(wisp_file):
                wisp_dir = Path(wisp_file).parent
                sample_id = pl.get("sample_id", "")
                snv_path = None
                if wisp_dir.is_dir():
                    for f in wisp_dir.iterdir():
                        fname = f.name
                        if sample_id and sample_id in fname and fname.endswith(".wisp.somatic_variants.tsv"):
                            snv_path = f
                            break
                    if not snv_path:
                        for f in wisp_dir.glob("*.wisp.somatic_variants.tsv"):
                            snv_path = f
                            break
                if snv_path and snv_path.is_file():
                    wisp_muts = parse_wisp_somatic_variants_tsv(snv_path, primary_driver_map)
                    raw_candidates.extend(wisp_muts)

            filtered_drivers = []
            seen_genes = set()
            for cand in raw_candidates:
                g = (cand.get("gene") or "").strip().upper()
                if not g or g in PHARMACOGENOMIC_GENES or g in seen_genes:
                    continue

                vaf = cand.get("vaf", 0.0)
                if vaf is None:
                    vaf = 0.0
                vaf_frac = vaf / 100.0 if vaf > 1.0 else vaf
                vaf_pct = vaf if vaf > 1.0 else vaf * 100.0

                is_called_in_primary = (g in primary_mutated_genes)
                is_high_vaf = (vaf_pct > 5.0)

                keep = False
                if is_called_in_primary:
                    # In primary: keep if detected in plasma (VAF > 0 or SampleAD > 0)
                    sample_ad = cand.get("sample_ad", 0)
                    if vaf_frac > 0.0 or sample_ad > 0 or cand.get("is_detected", False):
                        keep = True
                        prim = primary_driver_map[g]
                        if not cand.get("alt_type") or cand.get("alt_type") == "Mutation":
                            cand["alt_type"] = prim.get("alt_type", "Mutation")
                            cand["subtype"] = prim.get("subtype", cand["alt_type"])
                        if not cand.get("hgvs_protein") and prim.get("hgvs_protein"):
                            cand["hgvs_protein"] = prim.get("hgvs_protein")
                        if not cand.get("hgvs_coding") and prim.get("hgvs_coding"):
                            cand["hgvs_coding"] = prim.get("hgvs_coding")
                        if prim.get("is_status_unknown"):
                            cand["is_status_unknown"] = True
                            cand["status_flag"] = prim.get("status_flag")
                            cand["normal_dp"] = prim.get("normal_dp")
                            cand["normal_ad"] = prim.get("normal_ad")
                elif is_high_vaf:
                    # Novel / de novo in plasma with VAF > 5%
                    keep = True

                if keep:
                    filtered_drivers.append(cand)
                    seen_genes.add(g)

            pl["drivers"] = filtered_drivers
            pl["altered_genes"] = sorted(list(seen_genes))

            # CNVs: keep only high confidence driver CNVs or those in primary driver CNVs
            filtered_cnvs = {}
            if "purple_cnvs" in pl:
                for cg, cinfo in pl["purple_cnvs"].items():
                    if cinfo.get("is_driver") or cinfo.get("state") in {"Amplification", "Deep Deletion"}:
                        filtered_cnvs[cg] = cinfo
            pl["cnvs"] = filtered_cnvs


# ==============================================================================
# 4. ONCOPRINT MATRIX GENERATION
# ==============================================================================

def generate_oncoprint_matrix(cohort, genes_of_interest, top_n_cohort=20):
    """
    Builds the Gene x Sample alteration matrix for the Oncoprint:
    Clumps samples from the same patient together (Primary Tumor followed by chronological Plasma ctDNA),
    with a larger visual gap between patient groups.
    Adheres strictly to cBioPortal Visual Layering Hierarchy & Glyph Taxonomy:
      - Layer 0: Base Canvas / Wild-Type track (light grey bar)
      - Layer 1: CNA (AMP, HOMDEL, GAIN, HETLOSS) & Expression (UP, DOWN)
      - Layer 2: Mutation (Missense, Truncating, In-frame InDel, Promoter, Splice)
      - Layer 3: Structural Variant / Fusion (vertical purple bar)
      - Layer 4: Germline (white horizontal stripe overlay across mutation)
    """
    gene_freqs = cohort["gene_frequencies"]
    sorted_top = sorted(gene_freqs.keys(), key=lambda g: (gene_freqs[g], g), reverse=True)
    top_cohort_genes = sorted_top[:top_n_cohort]

    # Combine prioritized genes of interest, all Orange report driver genes, and top cohort genes
    final_genes = []
    seen = set()
    for g in genes_of_interest:
        g_clean = g.strip().upper()
        if g_clean and g_clean not in seen and g_clean not in PHARMACOGENOMIC_GENES:
            seen.add(g_clean)
            final_genes.append(g_clean)

    for g in cohort.get("all_orange_driver_genes", []):
        g_clean = g.strip().upper()
        if g_clean and g_clean not in seen and g_clean not in PHARMACOGENOMIC_GENES:
            seen.add(g_clean)
            final_genes.append(g_clean)

    for g in top_cohort_genes:
        if g not in seen and g not in PHARMACOGENOMIC_GENES:
            seen.add(g)
            final_genes.append(g)

    matrix = {g: {} for g in final_genes}
    cell_details = {g: {} for g in final_genes}
    gene_summary = {g: {"total_altered": 0, "pct": 0.0} for g in final_genes}

    # Clump samples by patient: Primary Tumor(s) first, then chronological Plasma ctDNA timepoint(s)
    ordered_samples = []
    patient_boundaries = set()
    patient_blocks = []

    # Detect duplicate sample IDs across different patients to prevent key/column collisions
    sample_patient_counts = defaultdict(set)
    for pid, prec in cohort["patients"].items():
        for t in prec.get("tumors", []):
            sample_patient_counts[t["tumor_id"]].add(pid)
            clean_sname = t.get("clean_sample_name") or extract_clean_sample_name(t["tumor_id"], pid) or t["tumor_id"]
            sample_patient_counts[clean_sname].add(pid)
        for pl in prec.get("plasmas", []):
            sample_patient_counts[pl["sample_id"]].add(pid)
            clean_sname = pl.get("clean_sample_name") or extract_clean_sample_name(pl["sample_id"], pid) or pl["sample_id"]
            sample_patient_counts[clean_sname].add(pid)

    for pid, prec in cohort["patients"].items():
        p_samples = []
        # 1. Primary tumor(s)
        for t in prec.get("tumors", []):
            tid = t["tumor_id"]
            clean_sname = t.get("clean_sample_name") or extract_clean_sample_name(tid, pid) or tid
            if len(sample_patient_counts[tid]) > 1 or len(sample_patient_counts[clean_sname]) > 1:
                col_id = f"{pid}_{tid}"
            else:
                col_id = tid

            p_samples.append({
                "sample_id": col_id,
                "raw_sample_id": tid,
                "clean_sample_name": clean_sname,
                "patient_id": pid,
                "sample_type": "Primary Tumor",
                "timepoint": "Primary",
                "is_plasma": False,
                "record": t,
                "purity_pct": t.get("purity_pct"),
                "tmb": t.get("tmb"),
                "msi_status": t.get("msi_status", "MSS"),
                "cuppa_cancer_type": t.get("cuppa_cancer_type", ""),
                "cuppa_likelihood": t.get("cuppa_likelihood"),
                "cuppa_predictions": t.get("cuppa_predictions", []),
            })
        # 2. Plasma ctDNA timepoints (already sorted chronologically)
        for pl in prec.get("plasmas", []):
            psid = pl["sample_id"]
            clean_sname = pl.get("clean_sample_name") or extract_clean_sample_name(psid, pid) or psid
            if len(sample_patient_counts[psid]) > 1 or len(sample_patient_counts[clean_sname]) > 1:
                col_id = f"{pid}_{psid}"
            else:
                col_id = psid

            first_tumor = prec["tumors"][0] if prec.get("tumors") else {}
            p_samples.append({
                "sample_id": col_id,
                "raw_sample_id": psid,
                "clean_sample_name": clean_sname,
                "patient_id": pid,
                "sample_type": "Plasma ctDNA",
                "timepoint": pl.get("timepoint", "Plasma"),
                "is_plasma": True,
                "record": pl,
                "wisp_purity_pct": pl.get("wisp_purity_pct"),
                "purple_purity_pct": pl.get("purple_purity_pct"),
                "is_mrd_positive": pl.get("is_mrd_positive"),
                "mrd_status": pl.get("mrd_status"),
                "cuppa_cancer_type": first_tumor.get("cuppa_cancer_type", ""),
                "cuppa_likelihood": first_tumor.get("cuppa_likelihood"),
                "cuppa_predictions": first_tumor.get("cuppa_predictions", []),
            })

        if p_samples:
            patient_boundaries.add(p_samples[-1]["sample_id"])
            patient_blocks.append({
                "patient_id": pid,
                "count": len(p_samples),
                "samples": [s["sample_id"] for s in p_samples]
            })
            ordered_samples.extend(p_samples)

    total_primary_tumors = len(cohort["all_tumors"]) if cohort["all_tumors"] else 1

    for s_info in ordered_samples:
        sid = s_info["sample_id"]
        pid = s_info["patient_id"]
        stype = s_info["sample_type"]
        tp = s_info["timepoint"]
        is_plasma = s_info["is_plasma"]
        rec = s_info["record"]
        is_boundary = (sid in patient_boundaries)

        for g in final_genes:
            # 1. CNA Layer (discrete cBioPortal calls)
            cna_entry = rec.get("cnvs", {}).get(g)
            cna_type = None
            cna_label = None
            if cna_entry:
                cbio_val = cna_entry.get("cbio_val", 0)
                if cbio_val == 2:
                    cna_type = "AMP"
                    cna_label = "Amplification"
                elif cbio_val == -2:
                    cna_type = "HOMDEL"
                    cna_label = "Deep Deletion"
                elif cbio_val == 1:
                    cna_type = "GAIN"
                    cna_label = "Gain"
                elif cbio_val == -1:
                    cna_type = "HETLOSS"
                    cna_label = "Shallow Deletion"

            # 2. Mutation Layer (Somatic + Germline)
            mut_list = []
            for d in rec.get("drivers", []):
                if d.get("gene") == g and d.get("alt_type") != "Structural Variant":
                    mut_list.append({
                        "type": d["alt_type"],
                        "subtype": d.get("subtype", d["alt_type"]),
                        "hgvs_protein": d.get("hgvs_protein", ""),
                        "hgvs_coding": d.get("hgvs_coding", ""),
                        "vaf": d.get("vaf", 0.0),
                        "depth_total": d.get("depth_total") or d.get("tumor_dp") or d.get("sample_dp"),
                        "depth_alt": d.get("depth_alt") or d.get("tumor_ad") or d.get("sample_ad"),
                        "variant_copy_number": d.get("variant_copy_number"),
                        "biallelic": d.get("biallelic"),
                        "is_germline": False,
                        "is_status_unknown": d.get("is_status_unknown", False),
                        "status_flag": d.get("status_flag", ""),
                        "normal_dp": d.get("normal_dp"),
                        "normal_ad": d.get("normal_ad"),
                    })
            germline_cnv = None
            for gd in rec.get("germline_drivers", []):
                if gd.get("gene") == g:
                    if gd.get("is_cnv"):
                        germline_cnv = gd.get("subtype") or "Germline Deletion"
                    elif gd.get("hgvs_protein") or gd.get("hgvs_coding") or gd.get("total_muts", 0) > 0:
                        mut_list.append({
                            "type": gd["alt_type"],
                            "subtype": gd.get("subtype", gd["alt_type"]),
                            "hgvs_protein": gd.get("hgvs_protein", ""),
                            "hgvs_coding": gd.get("hgvs_coding", ""),
                            "vaf": gd.get("vaf", 0.5),
                            "depth_total": gd.get("depth_total"),
                            "depth_alt": gd.get("depth_alt"),
                            "variant_copy_number": gd.get("variant_copy_number"),
                            "biallelic": gd.get("biallelic"),
                            "is_germline": True,
                            "is_status_unknown": False,
                        })

            # 3. Structural Variant / Fusion Layer
            fusions_list = rec.get("fusions", {}).get(g, [])
            fusion_name = fusions_list[0]["name"] if fusions_list else ""
            gene_svs = rec.get("structural_variants", {}).get(g, [])
            has_sv = bool(fusion_name) or bool(gene_svs)
            sv_desc = gene_svs[0]["name"] if gene_svs else (fusion_name or "")

            # 4. Expression Layer
            expr = rec.get("expression", {}).get(g, None)

            # Build readable summary tags for TSV export and tooltips
            summary_tags = []
            if cna_label:
                summary_tags.append(cna_label)
            if germline_cnv and germline_cnv not in summary_tags:
                summary_tags.append(germline_cnv)
            for m in mut_list:
                if m["type"] not in {"Structural Variant", "Deep Deletion", "Amplification"}:
                    parts = []
                    if m.get("is_status_unknown"):
                        parts.append("Status Unknown")
                    elif m.get("is_germline"):
                        parts.append("Germline")
                    else:
                        parts.append("Somatic")
                    if m.get("vaf") is not None and m.get("vaf") > 0:
                        vaf_pct = m["vaf"] if m["vaf"] > 1.0 else m["vaf"] * 100.0
                        depth_info = ""
                        if m.get("depth_alt") is not None and m.get("depth_total"):
                            depth_info = f" ({m['depth_alt']}/{m['depth_total']} reads)"
                        parts.append(f"VAF: {vaf_pct:.1f}%{depth_info}")
                    suffix = f" ({', '.join(parts)})" if parts else ""
                    if m.get("hgvs_protein"):
                        m_tag = f"{m['hgvs_protein']} ({m['type']}{', ' + ', '.join(parts) if parts else ''})"
                    elif m.get("hgvs_coding"):
                        m_tag = f"{m['hgvs_coding']} ({m['type']}{', ' + ', '.join(parts) if parts else ''})"
                    else:
                        m_tag = f"{m['type']}{suffix}"
                    if m_tag not in summary_tags:
                        summary_tags.append(m_tag)
            if sv_desc:
                summary_tags.append(f"SV ({sv_desc})")
            if expr:
                summary_tags.append(f"mRNA {expr}")

            # Assess alteration status
            is_altered = (
                (cna_type in {"AMP", "HOMDEL"} or (cna_type in {"GAIN", "HETLOSS"} and cna_entry.get("is_driver")))
                or bool(germline_cnv)
                or any(m["type"] not in {"Structural Variant", "Deep Deletion", "Amplification"} for m in mut_list)
                or has_sv
                or bool(expr)
            )

            # Alteration frequency % in cohort header represents primary tumor/patient frequency
            if is_altered and not is_plasma:
                gene_summary[g]["total_altered"] += 1

            matrix[g][sid] = ";".join(summary_tags)

            cell_details[g][sid] = {
                "gene": g,
                "tumor_id": sid,
                "sample_id": sid,
                "raw_sample_id": s_info.get("raw_sample_id", sid),
                "clean_sample_name": s_info.get("clean_sample_name", sid),
                "patient_id": pid,
                "sample_type": stype,
                "timepoint": tp,
                "is_plasma": is_plasma,
                "is_last_in_patient": is_boundary,
                "cna": cna_type,
                "cna_label": cna_label,
                "germline_cnv": germline_cnv,
                "mutations": [m for m in mut_list if m["type"] not in {"Structural Variant", "Deep Deletion", "Amplification"}],
                "is_germline": any(m.get("is_germline") for m in mut_list) or bool(germline_cnv),
                "fusion": bool(fusion_name),
                "fusion_name": fusion_name,
                "has_sv": has_sv,
                "structural_variant": bool(gene_svs),
                "structural_variants": gene_svs,
                "sv_name": sv_desc,
                "expression": expr,
                "is_profiled": True,
                "summary": "; ".join(summary_tags) if summary_tags else "Wild Type",
                "cna_min_cn": cna_entry.get("min_cn") if cna_entry else None,
                "cna_max_cn": cna_entry.get("max_cn") if cna_entry else None,
            }

    for g in final_genes:
        gene_summary[g]["pct"] = (gene_summary[g]["total_altered"] / total_primary_tumors) * 100.0

    # Sort final_genes: Altered frequency descending, then Orange drivers / genes of interest, then alpha
    final_genes.sort(
        key=lambda g: (
            gene_summary[g]["pct"],
            g in cohort.get("all_orange_driver_genes", []) or g in genes_of_interest,
            g
        ),
        reverse=True
    )

    return {
        "genes": final_genes,
        "matrix": matrix,
        "cell_details": cell_details,
        "summary": gene_summary,
        "tumors": [s["sample_id"] for s in ordered_samples],
        "samples": [s["sample_id"] for s in ordered_samples],
        "patient_blocks": patient_blocks,
        "patient_boundaries": list(patient_boundaries),
        "sample_metadata": {s["sample_id"]: {k: v for k, v in s.items() if k != "record"} for s in ordered_samples},
    }



# ==============================================================================
# 4b. SAMPLE FINGERPRINTING & IDENTITY CONFIRMATION
# ==============================================================================

def evaluate_cohort_fingerprints(cohort, output_dir=None, stride=20, max_workers=8):
    """
    Fingerprinting confirmation routine:
    Compares polymorphic germline heterozygous SNP loci across Normal, Tumor, and serial Plasma samples
    using Overlap Concordance (|A ∩ B| / min(|A|, |B|)) to confirm that all samples associated with the same
    patient originate from the same biological individual, regardless of sequencing platform (Illumina vs Ultima)
    or depth differences.
    """
    # Check for existing fingerprint cache in output_dir
    if output_dir:
        cache_file = Path(output_dir) / "cohort_fingerprints.json"
        if cache_file.is_file():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    cached_data = json.load(f)
                cached_patients = set(cached_data.get("patient_summaries", {}).keys())
                cohort_patients = set(cohort.get("patients", {}).keys())
                # If cached patients cover all cohort patients and samples were fingerprinted
                if cohort_patients and cohort_patients.issubset(cached_patients) and cached_data.get("total_samples", 0) > 0:
                    print(f"Fingerprinting already complete for all {len(cohort_patients)} cohort patient(s) ({cached_data.get('total_samples')} samples).")
                    print(f"Reusing cached fingerprint results from: {cache_file}")
                    return cached_data
            except Exception:
                pass
    def find_sample_amber_baf(amber_dirs, sample_tokens):
        for d in amber_dirs:
            if not d or not d.is_dir():
                continue
            for token in sample_tokens:
                if not token:
                    continue
                cand = d / f"{token}.amber.baf.tsv.gz"
                if cand.is_file():
                    return cand
            try:
                files = [f for f in os.listdir(d) if f.endswith(".amber.baf.tsv.gz")]
                for token in sample_tokens:
                    if not token:
                        continue
                    for f in files:
                        f_stem = f.replace(".amber.baf.tsv.gz", "")
                        if token == f_stem or token in f_stem or f_stem in token:
                            return d / f
                if len(files) == 1:
                    return d / files[0]
            except Exception:
                pass
        return None

    def read_baf_loci(baf_file, is_plasma=False):
        t_sites = set()
        n_sites = set()
        min_dp = 6 if is_plasma else 10
        try:
            file_stride = stride
            try:
                if os.path.getsize(baf_file) < 500_000:
                    file_stride = 1
            except Exception:
                pass

            with gzip.open(baf_file, "rt") as f:
                header = f.readline().strip().split("\t")
                t_dp_idx = header.index("tumorDepth") if "tumorDepth" in header else 4
                n_dp_idx = header.index("normalDepth") if "normalDepth" in header else 7
                for line in f:
                    parts = line.split("\t")
                    if file_stride > 1 and int(parts[1]) % file_stride != 0:
                        continue
                    pos = f"{parts[0]}:{parts[1]}"
                    if n_dp_idx is not None and int(parts[n_dp_idx]) >= 10:
                        n_sites.add(pos)
                    if t_dp_idx is not None and int(parts[t_dp_idx]) >= min_dp:
                        t_sites.add(pos)
        except Exception as e:
            print(f"Warning reading amber BAF {baf_file}: {e}")
        return t_sites, n_sites

    # 1. Discover all candidate BAF files to read
    baf_file_tasks = {}
    tumor_targets = []
    plasma_targets = []

    for p_id, p_rec in cohort.get("patients", {}).items():
        pdir = Path(p_rec.get("patient_dir", ""))
        tumors = p_rec.get("tumors", [])
        plasmas = p_rec.get("plasmas", [])

        # Process Tumors
        for t in tumors:
            t_id = t["tumor_id"]
            t_clean = t.get("clean_sample_name") or t_id
            amber_dirs = []
            if t.get("purple_purity_file"):
                pf = Path(t["purple_purity_file"])
                amber_dirs.append(pf.parent.parent / "amber")
            if pdir.is_dir():
                amber_dirs.append(pdir / t_id / "amber")
                amber_dirs.append(pdir / "amber")
                if pdir.parent.is_dir():
                    amber_dirs.append(pdir.parent / t_id / "amber")
                    amber_dirs.append(pdir.parent / p_id / "amber")

            target_amber = next((d for d in amber_dirs if d.is_dir()), None)
            if target_amber:
                norm_id = "Normal"
                try:
                    vcfs = [f for f in os.listdir(target_amber) if f.endswith(".amber.snp.vcf.gz")]
                    if vcfs:
                        norm_id = vcfs[0].replace(".amber.snp.vcf.gz", "")
                    elif t.get("normal_id"):
                        norm_id = t["normal_id"]
                except Exception:
                    if t.get("normal_id"):
                        norm_id = t["normal_id"]

                baf_file = find_sample_amber_baf([target_amber], [t_id, t.get("raw_tumor_id"), t_clean])
                if baf_file and baf_file.is_file():
                    baf_path_str = str(baf_file)
                    baf_file_tasks[baf_path_str] = False
                    tumor_targets.append({
                        "patient_id": p_id,
                        "tumor_id": t_id,
                        "tumor_clean": t_clean,
                        "normal_id": norm_id,
                        "date": t.get("date") or "",
                        "baf_path": baf_path_str
                    })

        # Process Plasmas
        for pl_idx, pl in enumerate(plasmas):
            pl_id = pl["sample_id"]
            clean_pname = pl.get("clean_sample_name") or pl_id
            tp_label = pl.get("timepoint") or f"T{pl_idx}"
            matched_tumor = pl.get("tumor_id", "Primary_Tumor")

            amber_cand_dirs = []
            if pl.get("wisp_file"):
                wf = Path(pl["wisp_file"])
                amber_cand_dirs.append(wf.parent.parent / "amber")
                amber_cand_dirs.append(wf.parent / "amber")
                amber_cand_dirs.append(wf.parent.parent.parent / "amber")
            if pdir.is_dir():
                amber_cand_dirs.append(pdir / "plasma" / pl_id / matched_tumor / "amber")
                amber_cand_dirs.append(pdir / "plasma" / pl_id / "amber")
                amber_cand_dirs.append(pdir / pl_id / "amber")
                amber_cand_dirs.append(pdir / "amber")
                if pdir.parent.is_dir():
                    amber_cand_dirs.append(pdir.parent / pl_id / "amber")
                    amber_cand_dirs.append(pdir.parent / f"{p_id}-{tp_label}" / "amber")

            pl_tokens = [pl_id, clean_pname]
            p_strip = re.sub(rf"^{p_id}[-_]?", "", pl_id, flags=re.IGNORECASE)
            if p_strip:
                pl_tokens.append(p_strip)
            c_strip = re.sub(rf"^{p_id}[-_]?", "", clean_pname, flags=re.IGNORECASE)
            if c_strip:
                pl_tokens.append(c_strip)

            pl_baf = find_sample_amber_baf(amber_cand_dirs, pl_tokens)
            if pl_baf and pl_baf.is_file():
                baf_path_str = str(pl_baf)
                baf_file_tasks[baf_path_str] = True
                plasma_targets.append({
                    "patient_id": p_id,
                    "sample_id": pl_id,
                    "clean_sample_name": clean_pname,
                    "timepoint": tp_label,
                    "date": pl.get("date") or "",
                    "days": pl.get("days"),
                    "chron_idx": 10 + pl_idx,
                    "baf_path": baf_path_str
                })

    # 2. Read unique BAF files in parallel
    def _read_worker(item):
        path_str, is_pl = item
        return path_str, read_baf_loci(path_str, is_plasma=is_pl)

    num_threads = min(max_workers or 8, max(1, os.cpu_count() or 4))
    baf_data = {}
    if baf_file_tasks:
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            for baf_path, (t_sites, n_sites) in executor.map(_read_worker, baf_file_tasks.items()):
                baf_data[baf_path] = (t_sites, n_sites)

    # 3. Assemble samples list
    samples = []

    for t_item in tumor_targets:
        p_id = t_item["patient_id"]
        t_id = t_item["tumor_id"]
        t_clean = t_item["tumor_clean"]
        norm_id = t_item["normal_id"]
        baf_path = t_item["baf_path"]

        if baf_path in baf_data:
            t_sites, n_sites = baf_data[baf_path]

            # Add Normal if not already present
            if not any(s["patient_id"] == p_id and s["sample_type"] == "Normal" for s in samples):
                if n_sites:
                    samples.append({
                        "patient_id": p_id,
                        "sample_id": norm_id,
                        "clean_sample_name": norm_id,
                        "sample_type": "Normal",
                        "timepoint": "Normal",
                        "date": "",
                        "chron_idx": 0,
                        "label": f"{p_id} | Normal ({norm_id})",
                        "sites": n_sites
                    })

            # Add Tumor
            if t_sites:
                samples.append({
                    "patient_id": p_id,
                    "sample_id": t_id,
                    "clean_sample_name": t_clean,
                    "sample_type": "Tumor",
                    "timepoint": "Baseline (Tumor)",
                    "date": t_item["date"],
                    "chron_idx": 1,
                    "label": f"{p_id} | Tumor ({t_clean})",
                    "sites": t_sites
                })

    for pl_item in plasma_targets:
        p_id = pl_item["patient_id"]
        pl_id = pl_item["sample_id"]
        clean_pname = pl_item["clean_sample_name"]
        tp_label = pl_item["timepoint"]
        baf_path = pl_item["baf_path"]

        if baf_path in baf_data:
            p_sites, _ = baf_data[baf_path]
            if p_sites:
                samples.append({
                    "patient_id": p_id,
                    "sample_id": pl_id,
                    "clean_sample_name": clean_pname,
                    "sample_type": "Plasma",
                    "timepoint": tp_label,
                    "date": pl_item["date"],
                    "days": pl_item["days"],
                    "chron_idx": pl_item["chron_idx"],
                    "label": f"{p_id} | Plasma {tp_label} ({clean_pname})",
                    "sites": p_sites
                })

    # Sort samples: Patient grouping -> Normal (0) -> Tumor (1) -> Chronological Plasmas (10+)
    def sample_sort_key(s):
        pid_sort = natural_sort_key(s["patient_id"])
        type_rank = {"Normal": 0, "Tumor": 1, "Plasma": 2}.get(s["sample_type"], 9)
        chron_idx = s.get("chron_idx", 0)
        return (pid_sort, type_rank, chron_idx)

    samples.sort(key=sample_sort_key)

    n_s = len(samples)
    matrix = [[0.0] * n_s for _ in range(n_s)]
    dist_matrix = [[0.0] * n_s for _ in range(n_s)]
    pairwise_list = []
    all_intra_concordances = []
    all_inter_concordances = []
    all_intra_dists = []
    all_inter_dists = []

    for i in range(n_s):
        matrix[i][i] = 100.0
        dist_matrix[i][i] = 0.0

    for i in range(n_s):
        s_i = samples[i]
        for j in range(i + 1, n_s):
            s_j = samples[j]
            inter = len(s_i["sites"] & s_j["sites"])
            min_len = min(len(s_i["sites"]), len(s_j["sites"]))
            conc = round((inter / min_len * 100.0), 2) if min_len > 0 else 0.0
            dist = round(1.0 - (conc / 100.0), 4)
            matrix[i][j] = matrix[j][i] = conc
            dist_matrix[i][j] = dist_matrix[j][i] = dist

            is_same = (s_i["patient_id"] == s_j["patient_id"])
            if is_same:
                all_intra_concordances.append(conc)
                all_intra_dists.append(dist)
                status = "CONFIRMED_MATCH" if conc >= 85.0 else "MISMATCH_SUSPECTED_SWAP"
            else:
                all_inter_concordances.append(conc)
                all_inter_dists.append(dist)
                status = "DISTINCT_DONOR" if conc < 50.0 else "UNEXPECTED_SIMILARITY_WARNING"

            pairwise_list.append({
                "patient_1": s_i["patient_id"],
                "sample_1": s_i["sample_id"],
                "clean_name_1": s_i.get("clean_sample_name") or s_i["sample_id"],
                "type_1": s_i["sample_type"],
                "timepoint_1": s_i.get("timepoint", ""),
                "patient_2": s_j["patient_id"],
                "sample_2": s_j["sample_id"],
                "clean_name_2": s_j.get("clean_sample_name") or s_j["sample_id"],
                "type_2": s_j["sample_type"],
                "timepoint_2": s_j.get("timepoint", ""),
                "is_same_patient": is_same,
                "shared_loci": inter,
                "min_loci": min_len,
                "total_loci": len(s_i["sites"] | s_j["sites"]),
                "concordance_pct": conc,
                "distance": dist,
                "status": status
            })

    patient_summaries = {}
    for p_id in cohort.get("patients", {}).keys():
        p_s = [s for s in samples if s["patient_id"] == p_id]
        norm = next((s for s in p_s if s["sample_type"] == "Normal"), None)
        tum = next((s for s in p_s if s["sample_type"] == "Tumor"), None)
        p_plasmas = [s for s in p_s if s["sample_type"] == "Plasma"]

        t_n_conc = None
        t_n_dist = None
        t_n_shared = 0
        if norm and tum:
            inter = len(norm["sites"] & tum["sites"])
            min_len = min(len(norm["sites"]), len(tum["sites"]))
            t_n_conc = round((inter / min_len * 100.0), 2) if min_len else 0.0
            t_n_dist = round(1.0 - (t_n_conc / 100.0), 4)
            t_n_shared = inter

        pl_comparisons = []
        intra_concordances = []
        intra_dists = []
        if t_n_conc is not None:
            intra_concordances.append(t_n_conc)
            intra_dists.append(t_n_dist)

        for pl in p_plasmas:
            pl_n_conc = None
            pl_n_dist = None
            pl_n_shared = 0
            if norm:
                inter = len(norm["sites"] & pl["sites"])
                min_len = min(len(norm["sites"]), len(pl["sites"]))
                pl_n_conc = round((inter / min_len * 100.0), 2) if min_len else 0.0
                pl_n_dist = round(1.0 - (pl_n_conc / 100.0), 4)
                pl_n_shared = inter

            pl_t_conc = None
            pl_t_dist = None
            pl_t_shared = 0
            if tum:
                inter = len(tum["sites"] & pl["sites"])
                min_len = min(len(tum["sites"]), len(pl["sites"]))
                pl_t_conc = round((inter / min_len * 100.0), 2) if min_len else 0.0
                pl_t_dist = round(1.0 - (pl_t_conc / 100.0), 4)
                pl_t_shared = inter

            if pl_n_conc is not None:
                intra_concordances.append(pl_n_conc)
                intra_dists.append(pl_n_dist)
            if pl_t_conc is not None:
                intra_concordances.append(pl_t_conc)
                intra_dists.append(pl_t_dist)

            is_pl_match = ((pl_n_conc is not None and pl_n_conc >= 85.0) or (pl_t_conc is not None and pl_t_conc >= 85.0))
            pl_comparisons.append({
                "plasma_id": pl["sample_id"],
                "timepoint": pl.get("timepoint") or "",
                "clean_sample_name": pl.get("clean_sample_name") or pl["sample_id"],
                "date": pl.get("date") or "",
                "days": pl.get("days"),
                "vs_normal_concordance": pl_n_conc,
                "vs_normal_dist": pl_n_dist,
                "vs_normal_shared": pl_n_shared,
                "vs_tumor_concordance": pl_t_conc,
                "vs_tumor_dist": pl_t_dist,
                "vs_tumor_shared": pl_t_shared,
                "is_match": is_pl_match
            })

        min_intra_c = min(intra_concordances) if intra_concordances else None
        max_intra_d = max(intra_dists) if intra_dists else None

        is_confirmed = all(c >= 85.0 for c in intra_concordances) if intra_concordances else (True if len(p_s) <= 1 else False)
        status_str = "CONFIRMED_SAME_DONOR" if is_confirmed else "MISMATCH_SUSPECTED_SWAP"
        if not intra_concordances:
            status_str = "SINGLE_SAMPLE_NO_COMPARATOR"

        patient_summaries[p_id] = {
            "patient_id": p_id,
            "normal_id": norm["sample_id"] if norm else "NA",
            "tumor_id": tum["sample_id"] if tum else "NA",
            "tumor_clean_name": (tum.get("clean_sample_name") or tum["sample_id"]) if tum else "NA",
            "tumor_vs_normal_concordance": t_n_conc,
            "tumor_vs_normal_dist": t_n_dist,
            "tumor_vs_normal_shared": t_n_shared,
            "plasma_comparisons": pl_comparisons,
            "is_confirmed": is_confirmed,
            "status": status_str,
            "min_intra_concordance": min_intra_c,
            "max_intra_dist": max_intra_d,
            "intra_concordances": intra_concordances,
        }

    mean_intra_c = round(sum(all_intra_concordances) / len(all_intra_concordances), 2) if all_intra_concordances else 0.0
    mean_inter_c = round(sum(all_inter_concordances) / len(all_inter_concordances), 2) if all_inter_concordances else 0.0
    mean_intra_d = round(sum(all_intra_dists) / len(all_intra_dists), 4) if all_intra_dists else 0.0
    mean_inter_d = round(sum(all_inter_dists) / len(all_inter_dists), 4) if all_inter_dists else 0.0

    result = {
        "samples": [
            {
                "sample_id": s["sample_id"],
                "clean_sample_name": s.get("clean_sample_name", s["sample_id"]),
                "patient_id": s["patient_id"],
                "sample_type": s["sample_type"],
                "timepoint": s.get("timepoint", ""),
                "date": s.get("date", ""),
                "label": s["label"],
                "n_loci": len(s["sites"])
            }
            for s in samples
        ],
        "matrix": matrix,
        "distance_matrix": dist_matrix,
        "pairwise_list": pairwise_list,
        "patient_summaries": patient_summaries,
        "total_samples": len(samples),
        "total_patients": len(patient_summaries),
        "confirmed_patients": sum(1 for ps in patient_summaries.values() if ps["is_confirmed"]),
        "mismatched_patients": sum(1 for ps in patient_summaries.values() if not ps["is_confirmed"] and ps["status"] != "SINGLE_SAMPLE_NO_COMPARATOR"),
        "mean_intra_concordance": mean_intra_c,
        "mean_inter_concordance": mean_inter_c,
        "mean_intra_dist": mean_intra_d,
        "mean_inter_dist": mean_inter_d,
    }

    if output_dir:
        try:
            cache_file = Path(output_dir) / "cohort_fingerprints.json"
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)
        except Exception:
            pass

    return result


# ==============================================================================
# 5. TSV EXPORTS
# ==============================================================================

def write_tsv_reports(cohort, oncoprint_data, output_dir):
    """Write master TSVs for WGTS cohort, samples, oncoprint, and MRD."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. cohort_wgts_summary.tsv (patient/tumor level)
    sum_file = out_path / "cohort_wgts_summary.tsv"
    fps = cohort.get("fingerprints", {})
    with open(sum_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "patient_id", "tumor_id", "cuppa_predicted_cancer_type", "cuppa_likelihood", "purity_pct", "ploidy", "curated_purity_pct", "curated_ploidy",
            "curation_method", "curation_notes", "gender", "status",
            "tumor_mean_coverage", "plasma_mean_coverage", "normal_mean_coverage", "tmb_mut_per_mb", "tmb_status",
            "msi_status", "whole_genome_duplication", "qc_status", "ichor_tumor_fraction",
            "fingerprint_status", "fingerprint_min_intra_concordance_pct", "fingerprint_max_intra_dist",
            "altered_driver_genes_count", "altered_genes", "plasma_samples_count",
            "latest_mrd_status", "purple_purity_file", "orange_json_file", "orange_pdf_file"
        ])
        for p_id, p_rec in cohort["patients"].items():
            tumors = p_rec.get("tumors", [])
            plasmas = p_rec.get("plasmas", [])
            latest_mrd = plasmas[-1]["mrd_status"] if plasmas else "No Plasma"
            plasma_cov = fmt_val(plasmas[0].get("weighted_depth"), ".2f") if plasmas else "NA"

            ps = fps.get("patient_summaries", {}).get(p_id, {})
            fp_status = ps.get("status", "NA")
            fp_min_c = fmt_val(ps.get("min_intra_concordance"), ".2f")
            fp_max_d = fmt_val(ps.get("max_intra_dist"), ".4f")

            if tumors:
                for t in tumors:
                    cuppa_lk = (fmt_val(t.get("cuppa_likelihood"), ".1f", mult=100.0) + "%") if t.get("cuppa_likelihood") is not None else "NA"
                    writer.writerow([
                        p_id,
                        t["tumor_id"],
                        t.get("cuppa_cancer_type", ""),
                        cuppa_lk,
                        fmt_val(t.get("purity_pct"), ".2f"),
                        fmt_val(t.get("ploidy"), ".2f"),
                        fmt_val(t.get("curated_purity"), ".2f", mult=100.0),
                        fmt_val(t.get("curated_ploidy"), ".2f"),
                        t.get("curation_method", ""),
                        t.get("curation_notes", ""),
                        t.get("gender", ""),
                        t.get("status", ""),
                        fmt_val(t.get("tumor_mean_coverage"), ".2f"),
                        plasma_cov,
                        fmt_val(t.get("normal_mean_coverage"), ".2f"),
                        fmt_val(t.get("tmb"), ".2f"),
                        t.get("tmb_status", ""),
                        t.get("msi_status", ""),
                        "Yes" if t.get("wgd") else "No",
                        t.get("qc_status", "PASS"),
                        fmt_val(t.get("ichor_tf"), ".4f"),
                        fp_status,
                        fp_min_c,
                        fp_max_d,
                        len(t.get("altered_genes", [])),
                        ";".join(t.get("altered_genes", [])),
                        len(plasmas),
                        latest_mrd,
                        t.get("purple_purity_file", ""),
                        t.get("orange_json_file", ""),
                        t.get("orange_pdf_file", ""),
                    ])
            else:
                writer.writerow([
                    p_id, "NA", "", "NA", "NA", "NA", "NA", "NA", "", "", "NA", "NA", "NA", plasma_cov, "NA", "NA", "", "NA", "NA", "NA", "NA", fp_status, fp_min_c, fp_max_d, 0, "", len(plasmas), latest_mrd, "", "", ""
                ])

    # 2. cohort_wgts_samples.tsv (tidy sample level)
    samples_file = out_path / "cohort_wgts_samples.tsv"
    with open(samples_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "patient_id", "sample_id", "sample_type", "timepoint", "cuppa_predicted_cancer_type", "cuppa_likelihood", "purity_pct",
            "ploidy", "curated_purity_pct", "curated_ploidy", "curation_method", "curation_notes",
            "tumor_coverage", "plasma_coverage", "normal_coverage", "mrd_status", "ichor_tf",
            "fingerprint_status",
            "tmb", "msi", "wgd", "qc_status", "orange_pdf_file"
        ])
        for t in cohort.get("all_tumors", []):
            t_ps = fps.get("patient_summaries", {}).get(t["patient_id"], {})
            t_fp_status = t_ps.get("status", "NA")
            cuppa_lk = (fmt_val(t.get("cuppa_likelihood"), ".1f", mult=100.0) + "%") if t.get("cuppa_likelihood") is not None else "NA"
            writer.writerow([
                t["patient_id"], t["tumor_id"], "Primary Tumor", "Baseline",
                t.get("cuppa_cancer_type", ""),
                cuppa_lk,
                fmt_val(t.get("purity_pct"), ".2f"), fmt_val(t.get("ploidy"), ".2f"),
                fmt_val(t.get("curated_purity"), ".2f", mult=100.0),
                fmt_val(t.get("curated_ploidy"), ".2f"),
                t.get("curation_method", ""),
                t.get("curation_notes", ""),
                fmt_val(t.get("tumor_mean_coverage"), ".2f"),
                "NA",
                fmt_val(t.get("normal_mean_coverage"), ".2f"),
                "Baseline (Tumor)",
                fmt_val(t.get("ichor_tf"), ".4f"),
                t_fp_status,
                fmt_val(t.get("tmb"), ".2f"),
                t.get("msi_status", ""),
                "Yes" if t.get("wgd") else "No",
                t.get("qc_status", "PASS"),
                t.get("orange_pdf_file", "")
            ])
        for pl in cohort.get("all_plasmas", []):
            p_match = cohort.get("patients", {}).get(pl["patient_id"], {})
            p_tumors = p_match.get("tumors", [])
            cuppa_c = p_tumors[0].get("cuppa_cancer_type", "") if p_tumors else ""
            cuppa_lk = (fmt_val(p_tumors[0].get("cuppa_likelihood"), ".1f", mult=100.0) + "%") if (p_tumors and p_tumors[0].get("cuppa_likelihood") is not None) else "NA"
            pl_ps = fps.get("patient_summaries", {}).get(pl["patient_id"], {})
            pl_fp_status = pl_ps.get("status", "NA")
            writer.writerow([
                pl["patient_id"], pl["sample_id"], "Plasma ctDNA", pl.get("timepoint", ""),
                cuppa_c,
                cuppa_lk,
                fmt_val(pl.get("snv_purity_pct"), ".4f"), fmt_val(pl.get("tumor_baseline_ploidy"), ".2f"),
                "NA", "NA", "", "",
                "NA",
                fmt_val(pl.get("weighted_depth"), ".2f"),
                "NA",
                pl.get("mrd_status", ""),
                fmt_val(pl.get("ichor_tf"), ".4f"),
                pl_fp_status,
                "NA", "NA", "NA", "PASS", ""
            ])

    # 3. cohort_wgts_oncoprint.tsv
    onco_file = out_path / "cohort_wgts_oncoprint.tsv"
    with open(onco_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        tumor_ids = oncoprint_data.get("tumors", [])
        writer.writerow(["gene", "altered_count", "frequency_pct"] + tumor_ids)
        for g in oncoprint_data.get("genes", []):
            g_sum = oncoprint_data.get("summary", {}).get(g, {})
            pct_str = (fmt_val(g_sum.get("pct"), ".1f") + "%") if g_sum.get("pct") is not None else "0.0%"
            row = [g, g_sum.get("total_altered", 0), pct_str]
            for tid in tumor_ids:
                row.append(oncoprint_data.get("matrix", {}).get(g, {}).get(tid, ""))
            writer.writerow(row)

    # 4. MRD Tables (if plasma available)
    if cohort.get("has_mrd_data"):
        mrd_file = out_path / "cohort_wgts_mrd_summary.tsv"
        with open(mrd_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "patient_id", "tumor_id", "plasma_sample_id", "timepoint", "overall_mrd_status",
                "tumor_baseline_purity_pct",
                "snv_mrd", "plasma_snv_purity_pct", "plasma_snv_ci_low_pct", "plasma_snv_ci_high_pct", "snv_lod_pct", "snv_pvalue",
                "loh_mrd", "plasma_loh_purity_pct", "loh_lod_pct", "loh_pvalue",
                "cnv_mrd", "plasma_cnv_purity_pct", "cnv_lod_pct", "cnv_pvalue",
                "wisp_mrd", "wisp_purity_pct", "wisp_method", "wisp_pvalue",
                "calc_variants", "total_variants", "plasma_coverage_depth", "ichor_tumor_fraction"
            ])
            for pl in cohort.get("all_plasmas", []):
                writer.writerow([
                    pl["patient_id"], pl.get("tumor_id", "NA"), pl["sample_id"], pl.get("timepoint", ""), pl.get("mrd_status", "NA"),
                    fmt_val(pl.get("tumor_baseline_purity"), ".2f", mult=100.0),
                    pl.get("snv_mrd", False),
                    fmt_val(pl.get("snv_purity_pct"), ".4f"),
                    fmt_val(pl.get("snv_purity_low_pct"), ".4f"),
                    fmt_val(pl.get("snv_purity_high_pct"), ".4f"),
                    fmt_val(pl.get("snv_lod_pct"), ".4f"),
                    fmt_val(pl.get("snv_pvalue"), ".2e"),
                    pl.get("loh_mrd", False),
                    fmt_val(pl.get("loh_purity_pct"), ".4f"),
                    fmt_val(pl.get("loh_lod_pct"), ".4f"),
                    fmt_val(pl.get("loh_pvalue"), ".2e"),
                    pl.get("cnv_mrd", False),
                    fmt_val(pl.get("cnv_purity_pct"), ".4f"),
                    fmt_val(pl.get("cnv_lod_pct"), ".4f"),
                    fmt_val(pl.get("cnv_pvalue"), ".2e"),
                    pl.get("wisp_mrd", False),
                    fmt_val(pl.get("wisp_purity_pct"), ".4f"),
                    pl.get("wisp_method") or "NA",
                    fmt_val(pl.get("wisp_pvalue"), ".2e"),
                    pl.get("calc_variants", 0),
                    pl.get("total_variants", 0),
                    fmt_val(pl.get("weighted_depth"), ".1f"),
                    fmt_val(pl.get("ichor_tf"), ".4f"),
                ])

        # Backward compatibility symlink/duplicate
        compat_file = out_path / "cohort_wisp_summary.tsv"
        try:
            import shutil
            shutil.copyfile(mrd_file, compat_file)
        except Exception:
            pass

    # 5. Fingerprinting TSV Reports
    if fps:
        fp_sum_file = out_path / "cohort_fingerprint_summary.tsv"
        with open(fp_sum_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "patient_id", "normal_id", "tumor_id", "tumor_clean_name",
                "tumor_vs_normal_concordance_pct", "tumor_vs_normal_shared_loci",
                "plasma_id", "plasma_timepoint", "plasma_clean_name", "plasma_date",
                "plasma_vs_normal_concordance_pct", "plasma_vs_normal_shared_loci",
                "plasma_vs_tumor_concordance_pct", "plasma_vs_tumor_shared_loci",
                "min_intra_patient_concordance_pct", "identity_status", "confirmation_call"
            ])
            for pid, ps in fps.get("patient_summaries", {}).items():
                plasmas = ps.get("plasma_comparisons", [])
                t_n_conc = fmt_val(ps.get("tumor_vs_normal_concordance"), ".2f")
                t_n_shared = ps.get("tumor_vs_normal_shared", 0)
                min_c = fmt_val(ps.get("min_intra_concordance"), ".2f")
                call = "CONFIRMED" if ps.get("is_confirmed") else "MISMATCH"
                if plasmas:
                    for plc in plasmas:
                        pl_n_conc = fmt_val(plc.get("vs_normal_concordance"), ".2f")
                        pl_n_shared = plc.get("vs_normal_shared", 0)
                        pl_t_conc = fmt_val(plc.get("vs_tumor_concordance"), ".2f")
                        pl_t_shared = plc.get("vs_tumor_shared", 0)
                        writer.writerow([
                            pid, ps.get("normal_id", "NA"), ps.get("tumor_id", "NA"), ps.get("tumor_clean_name", ps.get("tumor_id", "NA")),
                            t_n_conc, t_n_shared,
                            plc["plasma_id"], plc.get("timepoint", ""), plc.get("clean_sample_name", plc["plasma_id"]), plc.get("date", ""),
                            pl_n_conc, pl_n_shared,
                            pl_t_conc, pl_t_shared,
                            min_c, ps.get("status", "NA"), call
                        ])
                else:
                    writer.writerow([
                        pid, ps.get("normal_id", "NA"), ps.get("tumor_id", "NA"), ps.get("tumor_clean_name", ps.get("tumor_id", "NA")),
                        t_n_conc, t_n_shared,
                        "NA", "NA", "NA", "NA",
                        "NA", 0, "NA", 0,
                        min_c, ps.get("status", "NA"), call
                    ])

        fp_mat_file = out_path / "cohort_fingerprint_matrix.tsv"
        with open(fp_mat_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            sample_labels = [s["label"] for s in fps.get("samples", [])]
            writer.writerow(["sample"] + sample_labels)
            for i, row in enumerate(fps.get("matrix", [])):
                writer.writerow([sample_labels[i]] + [fmt_val(val, ".2f", default="0.00") for val in row])

    print(f"TSV summaries written to: {out_path}")


# ==============================================================================
# 6. STATIC PUBLICATION PLOTS (MATPLOTLIB)
# ==============================================================================

def generate_plasma_wisp_multimodal_plot(pl, tumor_rec, out_img):
    """
    Generates a publication-grade 6-track stacked comparison figure
    comparing Primary Tumor vs Plasma ctDNA across CNV, LOH, and SNV modalities,
    with 3 diagnostic scatter/summary panels on the right.
    """
    if not HAS_MATPLOTLIB or not HAS_PANDAS:
        return False

    wisp_file = pl.get("wisp_file")
    if not wisp_file or not os.path.exists(wisp_file):
        return False

    wisp_dir = Path(wisp_file).parent
    sample_id = pl.get("sample_id", "")
    patient_id = pl.get("patient_id", "Patient")
    tumor_id = pl.get("tumor_id", "Tumor")

    # 1. Discover WISP modality files
    cn_path = None
    loh_path = None
    snv_path = None
    calc_path = None

    if wisp_dir.is_dir():
        for f in wisp_dir.iterdir():
            fname = f.name
            if sample_id and sample_id in fname:
                if fname.endswith(".wisp.cn_segments.tsv"): cn_path = str(f)
                elif fname.endswith(".wisp.amber_loh.tsv"): loh_path = str(f)
                elif fname.endswith(".wisp.somatic_variants.tsv"): snv_path = str(f)
                elif fname.endswith(".wisp.cn_plot_calcs.tsv"): calc_path = str(f)
        for f in wisp_dir.iterdir():
            fname = f.name
            if not cn_path and fname.endswith(".wisp.cn_segments.tsv"): cn_path = str(f)
            if not loh_path and fname.endswith(".wisp.amber_loh.tsv"): loh_path = str(f)
            if not snv_path and fname.endswith(".wisp.somatic_variants.tsv"): snv_path = str(f)
            if not calc_path and fname.endswith(".wisp.cn_plot_calcs.tsv"): calc_path = str(f)

    # 2. Discover Tumor PURPLE Somatic CNV file
    tumor_cnv_file = None
    if tumor_rec and tumor_rec.get("purple_purity_file"):
        purple_dir = Path(tumor_rec["purple_purity_file"]).parent
        cand = purple_dir / f"{tumor_id}.purple.cnv.somatic.tsv"
        if cand.exists():
            tumor_cnv_file = str(cand)
        else:
            for f in purple_dir.glob("*.purple.cnv.somatic.tsv"):
                tumor_cnv_file = str(f)
                break

    # Read DataFrames
    try:
        df_summary = pd.read_csv(wisp_file, sep="\t") if os.path.exists(wisp_file) else None
    except Exception:
        df_summary = None

    try:
        df_cn = pd.read_csv(cn_path, sep="\t") if cn_path and os.path.exists(cn_path) else pd.DataFrame()
    except Exception:
        df_cn = pd.DataFrame()

    is_amber_fallback = False
    try:
        df_loh = pd.read_csv(loh_path, sep="\t") if loh_path and os.path.exists(loh_path) else pd.DataFrame()
    except Exception:
        df_loh = pd.DataFrame()

    if df_loh.empty:
        cand_amber_dirs = [
            wisp_dir.parent / "amber",
            wisp_dir.parent.parent / "amber",
        ]
        if pl.get("patient_dir"):
            p_pdir = Path(pl["patient_dir"])
            clean_sname = pl.get("clean_sample_name") or pl.get("sample_id")
            for c_d in [p_pdir / "plasma" / clean_sname / "amber", p_pdir / "plasma" / clean_sname / tumor_id / "amber"]:
                if c_d.is_dir() and c_d not in cand_amber_dirs:
                    cand_amber_dirs.append(c_d)
        amber_pcf = None
        for ad in cand_amber_dirs:
            if ad.is_dir():
                for pf in ad.glob("*.amber.baf.pcf"):
                    amber_pcf = pf
                    break
            if amber_pcf:
                break
        
        if amber_pcf and amber_pcf.is_file():
            try:
                raw_pcf = pd.read_csv(amber_pcf, sep="\t")
                if "Chromosome" in raw_pcf.columns and "Start" in raw_pcf.columns and "MeanRatio" in raw_pcf.columns:
                    df_loh = pd.DataFrame({
                        "Chromosome": raw_pcf["Chromosome"].astype(str),
                        "CnSegmentStart": pd.to_numeric(raw_pcf["Start"], errors="coerce"),
                        "CnSegmentEnd": pd.to_numeric(raw_pcf["End"], errors="coerce"),
                        "AvgAF": pd.to_numeric(raw_pcf["MeanRatio"], errors="coerce"),
                        "Sites": 50,
                    }).dropna()
                    is_amber_fallback = True
            except Exception:
                df_loh = pd.DataFrame()

    try:
        df_snv = pd.read_csv(snv_path, sep="\t") if snv_path and os.path.exists(snv_path) else pd.DataFrame()
    except Exception:
        df_snv = pd.DataFrame()

    try:
        df_calc = pd.read_csv(calc_path, sep="\t") if calc_path and os.path.exists(calc_path) else pd.DataFrame()
    except Exception:
        df_calc = pd.DataFrame()

    try:
        df_tumor_cn = pd.read_csv(tumor_cnv_file, sep="\t") if tumor_cnv_file and os.path.exists(tumor_cnv_file) else pd.DataFrame()
    except Exception:
        df_tumor_cn = pd.DataFrame()

    offsets, total_len = get_chrom_offsets()

    # Set up figure
    plt.rcParams["font.sans-serif"] = "DejaVu Sans"
    fig = plt.figure(figsize=(22, 17), facecolor="#0b0f14")

    gs = GridSpec(6, 4, width_ratios=[3.3, 3.3, 3.3, 2.5], 
                  height_ratios=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                  hspace=0.28, wspace=0.22, left=0.06, right=0.96, top=0.93, bottom=0.04)

    is_pos = pl.get("is_mrd_positive", False)
    status_str = "MRD POSITIVE" if is_pos else "MRD NEGATIVE"

    tp_str = f" ({pl.get('timepoint')})" if pl.get('timepoint') else ""
    fig.suptitle(f"PATIENT: {patient_id}   |   SAMPLE: {sample_id}{tp_str}\nTumor Baseline: {tumor_id}   |   WISP Decomposition & Integrated Call: {status_str}",
                 color="#f0f6fc", fontsize=16.5, fontweight="bold", y=0.975)

    def style_chrom_axis(ax, y_label, y_limits=None, baseline=None, show_x_labels=False):
        ax.set_facecolor("#161b22")
        ax.set_xlim(0, total_len)
        if y_limits: ax.set_ylim(y_limits)
        ax.set_ylabel(y_label, color="#c9d1d9", fontsize=8.5, fontweight="600")
        ax.tick_params(colors="#8b949e", labelsize=8)
        for spine in ax.spines.values(): spine.set_color("#30363d")
        
        chrom_ticks = []
        chrom_labels = []
        for i, (c, l) in enumerate(CHROM_LENGTHS.items()):
            start = offsets[c]
            end = start + l
            mid = (start + end) / 2
            chrom_ticks.append(mid)
            lbl = c.replace("chr", "")
            chrom_labels.append(lbl if (i % 2 == 0 or l > 80000000) else "")
            if i % 2 == 1:
                ax.axvspan(start, end, color="#21262d", alpha=0.35, zorder=0)
            else:
                ax.axvspan(start, end, color="#161b22", alpha=0.35, zorder=0)
            ax.axvline(end, color="#30363d", lw=0.6, alpha=0.5)

        ax.set_xticks(chrom_ticks)
        if show_x_labels:
            ax.set_xticklabels(chrom_labels, color="#8b949e", fontsize=7.5)
        else:
            ax.set_xticklabels([])

        if baseline is not None:
            ax.axhline(baseline, color="#8b949e", ls="--", lw=0.8, alpha=0.7)

    # TRACK 1A: TUMOR COPY NUMBER
    ax_cn_t = fig.add_subplot(gs[0, 0:3])
    style_chrom_axis(ax_cn_t, "Tumor Copy Number", y_limits=(0, 6.5), baseline=2.0, show_x_labels=False)
    ax_cn_t.set_title("1. COPY NUMBER (CNV) — Primary Tumor Copy Number vs. Plasma GC-Corrected Depth Ratio", color="#58a6ff", fontsize=11.5, fontweight="bold", loc="left", pad=6)

    src_cn = df_tumor_cn if not df_tumor_cn.empty else df_cn
    if not src_cn.empty and ("Chromosome" in src_cn or "chromosome" in src_cn) and ("start" in src_cn or "SegmentStart" in src_cn or "startPos" in src_cn):
        chr_col = "chromosome" if "chromosome" in src_cn.columns else "Chromosome"
        start_col = "start" if "start" in src_cn.columns else ("SegmentStart" if "SegmentStart" in src_cn.columns else "startPos")
        end_col = "end" if "end" in src_cn.columns else ("SegmentEnd" if "SegmentEnd" in src_cn.columns else "endPos")
        cn_col = "copyNumber" if "copyNumber" in src_cn.columns else "CopyNumber"

        valid_cn = src_cn[src_cn[chr_col].astype(str).isin(offsets)].copy()
        if not valid_cn.empty:
            c_offsets = valid_cn[chr_col].astype(str).map(offsets)
            x0 = (c_offsets + valid_cn[start_col]).values
            x1 = (c_offsets + valid_cn[end_col]).values
            cn = valid_cn[cn_col].values

            gain_mask = cn > 2.6
            loss_mask = cn < 1.4
            neut_mask = ~gain_mask & ~loss_mask

            for m, col in [(gain_mask, "#f85149"), (loss_mask, "#388bfd"), (neut_mask, "#8b949e")]:
                if np.any(m):
                    lines = [((x0[i], cn[i]), (x1[i], cn[i])) for i in np.where(m)[0]]
                    ax_cn_t.add_collection(LineCollection(lines, colors=col, linewidths=2.2))

    # TRACK 1B: PLASMA GC RATIO
    ax_cn_p = fig.add_subplot(gs[1, 0:3])
    style_chrom_axis(ax_cn_p, "Plasma GC-Ratio", y_limits=(0.75, 1.35), baseline=1.0, show_x_labels=True)

    if not df_cn.empty and "Chromosome" in df_cn and "GcRatioMean" in df_cn and "SegmentStart" in df_cn:
        valid_df_cn = df_cn[df_cn["Chromosome"].astype(str).isin(offsets)].copy()
        if not valid_df_cn.empty:
            c_offsets = valid_df_cn["Chromosome"].astype(str).map(offsets)
            x0 = (c_offsets + valid_df_cn["SegmentStart"]).values
            x1 = (c_offsets + valid_df_cn["SegmentEnd"]).values
            gc = valid_df_cn["GcRatioMean"].values
            tcn = valid_df_cn.get("CopyNumber", pd.Series(2.0, index=valid_df_cn.index)).values

            gain_mask = tcn > 2.6
            loss_mask = tcn < 1.4
            neut_mask = ~gain_mask & ~loss_mask

            for m, col in [(gain_mask, "#f85149"), (loss_mask, "#388bfd"), (neut_mask, "#8b949e")]:
                if np.any(m):
                    lines = [((x0[i], gc[i]), (x1[i], gc[i])) for i in np.where(m)[0]]
                    ax_cn_p.add_collection(LineCollection(lines, colors=col, linewidths=1.8, alpha=0.85))

    # RIGHT INSET 1: WISP GC-ratio vs Tumor CN scatter
    ax_cn_fit = fig.add_subplot(gs[0:2, 3])
    ax_cn_fit.set_facecolor("#161b22")
    ax_cn_fit.tick_params(colors="#8b949e", labelsize=8)
    for spine in ax_cn_fit.spines.values(): spine.set_color("#30363d")
    ax_cn_fit.set_title("WISP CNV Regression Fit", color="#58a6ff", fontsize=10.5, fontweight="bold", pad=8)
    ax_cn_fit.set_xlabel("Tumor Copy Number", color="#c9d1d9", fontsize=9)
    ax_cn_fit.set_ylabel("Plasma GC Ratio", color="#c9d1d9", fontsize=9)

    cn_pur = safe_float(pl.get("cnv_purity_pct"), 0.0)
    cn_lod = safe_float(pl.get("cnv_lod_pct"), 0.0)
    cn_mrd = bool(pl.get("cnv_mrd", False))
    cn_call = "POSITIVE" if cn_mrd else "Negative (below LOD)"
    call_col = "#3fb950" if cn_mrd else "#8b949e"

    if not df_cn.empty and "CopyNumber" in df_cn and "GcRatioMean" in df_cn:
        t_cn = df_cn["CopyNumber"].values
        p_gc = df_cn["GcRatioMean"].values
        ax_cn_fit.scatter(t_cn, p_gc, c="#d29922", s=14, alpha=0.6, edgecolors="none")
        ax_cn_fit.axhline(1.0, color="#8b949e", ls="--", lw=0.8)
        ax_cn_fit.axvline(2.0, color="#8b949e", ls="--", lw=0.8)

        if not df_calc.empty and "FitCoefficient" in df_calc and "FitIntercept" in df_calc:
            slope = df_calc["FitCoefficient"].values[0]
            intercept = df_calc["FitIntercept"].values[0]
            x_line = np.linspace(max(0, min(t_cn)), min(8, max(t_cn)), 50)
            y_line = slope * (x_line - 2.0) + intercept
            ax_cn_fit.plot(x_line, y_line, color="#f85149", lw=2, label="Fit Line")

    ax_cn_fit.text(0.05, 0.90, f"CNV Purity: {cn_pur:.2f}%\nCNV LOD:   {cn_lod:.2f}%\nStatus:    {cn_call}",
                   transform=ax_cn_fit.transAxes, color=call_col, fontsize=9, fontfamily="monospace", fontweight="bold",
                   bbox=dict(boxstyle="round,pad=0.4", facecolor="#0d1117", edgecolor="#30363d", alpha=0.9))

    # TRACK 2A: TUMOR MINOR ALLELE COPY NUMBER (LOH BASELINE)
    ax_loh_t = fig.add_subplot(gs[2, 0:3])
    style_chrom_axis(ax_loh_t, "Tumor Minor CN", y_limits=(-0.1, 2.2), baseline=1.0, show_x_labels=False)
    ax_loh_t.set_title("2. LOSS OF HETEROZYGOSITY (LOH) — Tumor Allele Loss vs. Phased Plasma SNP Allelic Imbalance", color="#bc8cff", fontsize=11.5, fontweight="bold", loc="left", pad=6)

    if not df_tumor_cn.empty and "minorAlleleCopyNumber" in df_tumor_cn and "chromosome" in df_tumor_cn:
        valid_t_cn = df_tumor_cn[df_tumor_cn["chromosome"].astype(str).isin(offsets)].copy()
        if not valid_t_cn.empty:
            c_offsets = valid_t_cn["chromosome"].astype(str).map(offsets)
            x0 = (c_offsets + valid_t_cn["start"]).values
            x1 = (c_offsets + valid_t_cn["end"]).values
            mcn = valid_t_cn["minorAlleleCopyNumber"].values
            is_loh = mcn < 0.2

            for m, col in [(is_loh, "#a371f7"), (~is_loh, "#8b949e")]:
                if np.any(m):
                    lines = [((x0[i], mcn[i]), (x1[i], mcn[i])) for i in np.where(m)[0]]
                    ax_loh_t.add_collection(LineCollection(lines, colors=col, linewidths=2.2, alpha=0.9))

            for i in np.where(is_loh)[0]:
                ax_loh_t.fill_between([x0[i], x1[i]], 0, mcn[i], color="#a371f7", alpha=0.25)

    # TRACK 2B: PLASMA ALLELE FREQUENCY AT LOH REGIONS
    ax_loh_p = fig.add_subplot(gs[3, 0:3])
    style_chrom_axis(ax_loh_p, "Plasma Allele Freq", y_limits=(0.43, 0.57), baseline=0.50, show_x_labels=True)

    if not df_loh.empty and "Chromosome" in df_loh and "AvgAF" in df_loh and "CnSegmentStart" in df_loh:
        valid_loh = df_loh[df_loh["Chromosome"].astype(str).isin(offsets)].copy()
        if not valid_loh.empty:
            c_offsets = valid_loh["Chromosome"].astype(str).map(offsets)
            x0 = (c_offsets + valid_loh["CnSegmentStart"]).values
            x1 = (c_offsets + valid_loh["CnSegmentEnd"]).values
            mid_x = (x0 + x1) / 2.0
            af = valid_loh["AvgAF"].values
            sites = valid_loh.get("Sites", pd.Series(1, index=valid_loh.index)).values
            size_pt = np.clip(np.sqrt(sites) * 2.2, 14, 130)

            is_imb = np.abs(af - 0.50) > 0.005
            for m, col in [(is_imb, "#a371f7"), (~is_imb, "#8b949e")]:
                if np.any(m):
                    lines = [((x0[i], af[i]), (x1[i], af[i])) for i in np.where(m)[0]]
                    ax_loh_p.add_collection(LineCollection(lines, colors=col, linewidths=1.8, alpha=0.8))
                    ax_loh_p.scatter(mid_x[m], af[m], s=size_pt[m], c=col, alpha=0.8, edgecolors="#ffffff", linewidths=0.4)

    # RIGHT INSET 2: LOH Summary Card
    ax_loh_sum = fig.add_subplot(gs[2:4, 3])
    ax_loh_sum.set_facecolor("#161b22")
    ax_loh_sum.set_xticks([])
    ax_loh_sum.set_yticks([])
    for spine in ax_loh_sum.spines.values(): spine.set_color("#30363d")
    ax_loh_sum.set_title("WISP LOH Aggregation Power", color="#bc8cff", fontsize=10.5, fontweight="bold", pad=8)

    loh_pur = safe_float(pl.get("loh_purity_pct"), 0.0)
    loh_pval = safe_float(pl.get("loh_pvalue"), 1.0)
    loh_lod = safe_float(pl.get("loh_lod_pct"), 0.0)
    loh_mrd = bool(pl.get("loh_mrd", False))
    loh_call = "POSITIVE" if loh_mrd else "Negative (below LOD)"

    if is_amber_fallback:
        loh_regions = len(df_loh)
        med_imb = float(np.median(np.abs(df_loh["AvgAF"] - 0.50))) if not df_loh.empty else 0.0
        sum_text = (
            f"LOH Methodology (AMBER BAF):\n"
            f"  Track 2B displays AMBER\n"
            f"  segmented BAF. Native WISP\n"
            f"  LOH pending Step 2 re-run.\n"
            f"---------------------------\n"
            f"AMBER Segments:  {loh_regions:,}\n"
            f"Median Imbalance:{med_imb:.4f}\n"
            f"---------------------------\n"
            f"WISP LOH Status: Pending Re-run\n"
            f"(Step 2 patch applied)"
        )
        ax_loh_sum.text(0.06, 0.50, sum_text, color="#f0f6fc", fontsize=9.2, fontfamily="monospace", va="center",
                        bbox=dict(boxstyle="round,pad=0.5", facecolor="#0d1117", edgecolor="#bc8cff", alpha=0.9))
    else:
        loh_sites = 0
        loh_regions = 0
        loh_frags = 0
        if df_summary is not None and not df_summary.empty and "LOHSiteCount" in df_summary:
            loh_sites = df_summary["LOHSiteCount"].values[0]
            loh_regions = df_summary["LOHRegionCount"].values[0]
            loh_frags = df_summary["LOHFragments"].values[0]
        elif not df_loh.empty:
            loh_regions = len(df_loh)
            loh_sites = df_loh["Sites"].sum() if "Sites" in df_loh else 0
            loh_frags = df_loh["Depth"].sum() if "Depth" in df_loh else 0

        sum_text = (
            f"LOH Methodology:\n"
            f"  Tumor LOH establishes\n"
            f"  which allele is deleted.\n"
            f"  Plasma tests if germline\n"
            f"  SNPs deviate from 50/50.\n"
            f"---------------------------\n"
            f"LOH Segments:    {loh_regions:,}\n"
            f"Phased SNPs:     {loh_sites:,}\n"
            f"Aggregated Reads:{loh_frags:,}\n"
            f"---------------------------\n"
            f"LOH Purity:      {loh_pur:.2f}%\n"
            f"p-value:         {loh_pval:.2e}\n"
            f"LOH LOD:         {loh_lod:.2f}%\n"
            f"Status:          {loh_call}"
        )
        ax_loh_sum.text(0.06, 0.50, sum_text, color="#f0f6fc", fontsize=9.2, fontfamily="monospace", va="center",
                        bbox=dict(boxstyle="round,pad=0.5", facecolor="#0d1117", edgecolor="#bc8cff", alpha=0.9))

    # TRACK 3A: TUMOR SOMATIC VARIANTS (VAF)
    ax_snv_t = fig.add_subplot(gs[4, 0:3])
    style_chrom_axis(ax_snv_t, "Tumor VAF (%)", y_limits=(0, 100), baseline=50.0, show_x_labels=False)
    ax_snv_t.set_title("3. SOMATIC SNVs & INDELs — Tumor Variant Allele Frequency vs. Plasma Detection", color="#3fb950", fontsize=11.5, fontweight="bold", loc="left", pad=6)

    # TRACK 3B: PLASMA SOMATIC VARIANTS (VAF)
    ax_snv_p = fig.add_subplot(gs[5, 0:3])
    style_chrom_axis(ax_snv_p, "Plasma VAF (%)", y_limits=(0, 25), baseline=None, show_x_labels=True)

    if not df_snv.empty and "Chromosome" in df_snv and "TumorAD" in df_snv and "TumorDP" in df_snv:
        valid_snv = df_snv[df_snv["Chromosome"].astype(str).isin(offsets)].copy()
        if not valid_snv.empty:
            valid_snv["TumorVAF"] = (valid_snv["TumorAD"] / valid_snv["TumorDP"].replace(0, 1)) * 100.0
            valid_snv["PlasmaVAF"] = (valid_snv["SampleAD"] / valid_snv["SampleDP"].replace(0, 1)) * 100.0

            c_offsets = valid_snv["Chromosome"].astype(str).map(offsets)
            valid_snv["gx"] = c_offsets + valid_snv["Position"]

            det_mask = valid_snv["SampleAD"] > 0
            df_det = valid_snv[det_mask]
            df_undet = valid_snv[~det_mask]

            # Downsample background undetected SNVs to max 15,000 for fast vector rendering
            if len(df_undet) > 15000:
                df_undet_plot = df_undet.sample(n=15000, random_state=42)
            else:
                df_undet_plot = df_undet

            # TRACK 3A: TUMOR SOMATIC VARIANTS
            if not df_undet_plot.empty:
                ax_snv_t.scatter(df_undet_plot["gx"], df_undet_plot["TumorVAF"], color="#8b949e", s=8, alpha=0.3, edgecolors="none")
            if not df_det.empty:
                ax_snv_t.scatter(df_det["gx"], df_det["TumorVAF"], color="#3fb950", s=20, alpha=0.85, edgecolors="none")

            # TRACK 3B: PLASMA SOMATIC VARIANTS
            if not df_undet_plot.empty:
                ax_snv_p.scatter(df_undet_plot["gx"], np.zeros(len(df_undet_plot)), color="#8b949e", s=6, alpha=0.25, zorder=2, edgecolors="none")
            if not df_det.empty:
                ax_snv_p.scatter(df_det["gx"], df_det["PlasmaVAF"], color="#3fb950", s=28, zorder=5, edgecolors="#ffffff", lw=0.4)

                # Annotate top detected driver genes
                high_det = df_det[(df_det["PlasmaVAF"] > 6.0) & df_det["Gene"].notna() & (df_det["Gene"] != "nan")]
                for _, r in high_det.iterrows():
                    gene = str(r["Gene"]).strip()
                    if gene and gene != "nan":
                        ax_snv_p.annotate(gene, (r["gx"], r["PlasmaVAF"]), textcoords="offset points", xytext=(0, 5),
                                         ha="center", fontsize=7.5, color="#58a6ff", fontweight="bold")

    # RIGHT INSET 3: Tumor VAF vs Plasma VAF scatter
    ax_snv_scat = fig.add_subplot(gs[4:6, 3])
    ax_snv_scat.set_facecolor("#161b22")
    ax_snv_scat.tick_params(colors="#8b949e", labelsize=8)
    for spine in ax_snv_scat.spines.values(): spine.set_color("#30363d")
    ax_snv_scat.set_title("Tumor VAF vs. Plasma VAF", color="#3fb950", fontsize=10.5, fontweight="bold", pad=8)
    ax_snv_scat.set_xlabel("Tumor VAF (%)", color="#c9d1d9", fontsize=9)
    ax_snv_scat.set_ylabel("Plasma VAF (%)", color="#c9d1d9", fontsize=9)

    snv_pur = safe_float(pl.get("snv_purity_pct"), 0.0)
    snv_pval = safe_float(pl.get("snv_pvalue"), 1.0)
    snv_lod = safe_float(pl.get("snv_lod_pct"), 0.0)
    snv_mrd = bool(pl.get("snv_mrd", False))
    snv_call = "POSITIVE" if snv_mrd else "Negative (below LOD)"
    scall_col = "#3fb950" if snv_mrd else "#8b949e"

    if not df_snv.empty and "PlasmaVAF" in df_snv:
        det_mask = df_snv["SampleAD"] > 0
        ax_snv_scat.scatter(df_snv.loc[~det_mask, "TumorVAF"], df_snv.loc[~det_mask, "PlasmaVAF"],
                           c="#8b949e", s=10, alpha=0.25, label="0 reads")
        ax_snv_scat.scatter(df_snv.loc[det_mask, "TumorVAF"], df_snv.loc[det_mask, "PlasmaVAF"],
                           c="#3fb950", s=26, alpha=0.85, edgecolors="#ffffff", lw=0.4, label="Detected (AD>0)")

    ax_snv_scat.text(0.05, 0.90, f"SNV Purity: {snv_pur:.2f}%\np-value:    {snv_pval:.3f}\nSNV LOD:    {snv_lod:.2f}%\nStatus:     {snv_call}",
                     transform=ax_snv_scat.transAxes, color=scall_col, fontsize=9, fontfamily="monospace", fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.4", facecolor="#0d1117", edgecolor="#30363d", alpha=0.9))

    os.makedirs(os.path.dirname(out_img), exist_ok=True)
    plt.savefig(out_img, dpi=200, facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)
    return True


def render_mrd_trajectory_on_axis(ax_p, plasmas, modality="SNV", patient_id="Patient",
                                  x_labels=None, show_legend=True, show_title=True,
                                  is_bottom_track=True):
    """
    Renders a clinical ctDNA trajectory on a matplotlib Axes matching publication standards:
      - Log-scale y-axis from 0.001% to 100% labeled 'Tumor Fraction (%)'
      - Dotted vertical line at x=0 for 'WGS Baseline Timestamp'
      - Continuous dashed cutoff line with shaded grey region below for Limit of Detection (NO scatter points)
      - Solid trajectory line with circular markers for detected ctDNA
      - Dashed trajectory line to grey circular markers for below-LOD measured values
      - Smooth continuous ribbon for 95% Confidence Interval (NO per-point boxes)
      - Exact same visual layout and legend across all modalities (SNV / CNV / LOH)
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    mod_configs = {
        "SNV": {
            "metric_key": "snv_purity_pct",
            "fallback_metric": "wisp_purity_pct",
            "mrd_key": "snv_mrd",
            "lod_key": "snv_lod_pct",
            "fallback_lod": "wisp_lod_pct",
            "low_key": "snv_purity_low_pct",
            "high_key": "snv_purity_high_pct",
            "color": "#ef4444",
            "edge_color": "#b91c1c",
            "ci_color": "#fee2e2",
            "title_mod": "SNV Modality",
            "label_mod": "SNV",
            "y_label": "Tumor Fraction (%) — SNV ctDNA",
        },
        "CNV": {
            "metric_key": "cnv_purity_pct",
            "fallback_metric": None,
            "mrd_key": "cnv_mrd",
            "lod_key": "cnv_lod_pct",
            "fallback_lod": None,
            "low_key": "cnv_purity_low_pct",
            "high_key": "cnv_purity_high_pct",
            "color": "#a855f7",
            "edge_color": "#7e22ce",
            "ci_color": "#f3e8ff",
            "title_mod": "CNV Modality",
            "label_mod": "CNV",
            "y_label": "Tumor Fraction (%) — CNV ctDNA",
        },
        "LOH": {
            "metric_key": "loh_purity_pct",
            "fallback_metric": None,
            "mrd_key": "loh_mrd",
            "lod_key": "loh_lod_pct",
            "fallback_lod": None,
            "low_key": "loh_purity_low_pct",
            "high_key": "loh_purity_high_pct",
            "color": "#f59e0b",
            "edge_color": "#b45309",
            "ci_color": "#fef3c7",
            "title_mod": "LOH Modality",
            "label_mod": "LOH",
            "y_label": "Tumor Fraction (%) — LOH ctDNA",
            "marker": "o",
        },
        "ICHOR": {
            "metric_key": "ichor_purity_pct",
            "fallback_metric": None,
            "mrd_key": "ichor_mrd",
            "lod_key": "ichor_lod_pct",
            "fallback_lod": None,
            "low_key": None,
            "high_key": None,
            "color": "#d97706",
            "edge_color": "#92400e",
            "ci_color": "#fef3c7",
            "title_mod": "1X ichorCNA Modality (Independent Copy Number)",
            "label_mod": "ichorCNA",
            "y_label": "Tumor Fraction (%) — ichorCNA",
            "marker": "D",
        },
    }

    cfg = mod_configs.get(modality.upper(), mod_configs["SNV"])
    num_pts = len(plasmas)
    x_coords = np.arange(num_pts)

    y_min = 0.0005
    y_max = 105.0
    ax_p.set_yscale("log")
    ax_p.set_ylim(y_min, y_max)
    ax_p.set_xlim(-0.6, max(num_pts - 0.4, 0.6))

    # Y-axis ticks and labels
    ax_p.set_yticks([0.001, 0.01, 0.1, 1.0, 10.0, 100.0])
    ax_p.set_yticklabels(["0.001%", "0.01%", "0.1%", "1%", "10%", "100%"], fontsize=9.5)
    ax_p.set_ylabel(cfg["y_label"], fontsize=10.5, fontweight="bold")

    if show_title:
        ax_p.set_title(f"{patient_id} — Longitudinal ctDNA Trajectory ({cfg['title_mod']})", fontsize=11.5, fontweight="bold", pad=10)

    # Clean grid
    ax_p.grid(True, which="both", linestyle="-", color="#e5e7eb", linewidth=0.8, alpha=0.8, zorder=1)

    # WGS Baseline Timestamp line (at index 0 or baseline)
    ax_p.axvline(0, color="black", linestyle=":", linewidth=1.5, zorder=2)

    # Baseline horizontal dashed reference line (0.05%)
    ax_p.axhline(0.05, color="#9ca3af", linestyle="--", linewidth=0.8, alpha=0.7, zorder=2)

    # Midpoint vertical line if multiple points
    if num_pts > 4:
        ax_p.axvline(num_pts // 2, color="#9ca3af", linestyle="-.", linewidth=0.8, alpha=0.6, zorder=2)

    # LOD array
    lods = []
    for pl in plasmas:
        lod_val = pl.get(cfg["lod_key"]) or (pl.get(cfg["fallback_lod"]) if cfg["fallback_lod"] else None) or 0.002
        if lod_val <= 0:
            lod_val = 0.002
        lods.append(lod_val)
    lods = np.array(lods)

    # Below Limit of Detection: shaded grey region
    ax_p.fill_between(x_coords, y_min, lods, color="#e2e8f0", alpha=0.9, zorder=2)

    # Limit of Detection stepped dashed cutoff line (NO points, clean line only!)
    ax_p.plot(x_coords, lods, color="#94a3b8", linestyle="--", linewidth=1.2, zorder=2.5)

    # Tumor fractions, CI, and marker partitioning
    all_plot_y = []
    detected_indices = []
    detected_y = []
    below_lod_indices = []
    below_lod_y = []
    ci_lows = []
    ci_highs = []

    for idx, pl in enumerate(plasmas):
        is_det = pl.get(cfg["mrd_key"]) if pl.get(cfg["mrd_key"]) is not None else pl.get("is_mrd_positive", False)
        val = pl.get(cfg["metric_key"])
        if val is None and cfg["fallback_metric"]:
            val = pl.get(cfg["fallback_metric"])
        val = val or 0.0

        c_low = pl.get(cfg["low_key"])
        c_high = pl.get(cfg["high_key"])
        lod_val = lods[idx]

        is_detected = bool(is_det and (val > lod_val))
        actual_val = max(val, 0.0)
        plot_y = max(actual_val, y_min * 1.15) if actual_val > 0 else (y_min * 1.15)
        all_plot_y.append(plot_y)

        if is_detected:
            detected_indices.append(idx)
            detected_y.append(plot_y)
            c_low_val = c_low if (c_low is not None and c_low > 0) else max(plot_y * 0.75, lod_val)
            c_high_val = c_high if (c_high is not None and c_high > 0) else min(plot_y * 1.25, 100.0)
        else:
            below_lod_indices.append(idx)
            below_lod_y.append(plot_y)
            # Below-LOD / negative sample confidence interval:
            # Upper bound: use provided c_high if > 0, otherwise the sample Limit of Detection (LOD)
            c_high_val = c_high if (c_high is not None and c_high > 0) else lod_val
            c_high_val = max(c_high_val, lod_val)
            # Lower bound: use provided c_low if > 0, otherwise the plot baseline floor
            c_low_val = c_low if (c_low is not None and c_low > 0) else (y_min * 1.15)
            c_low_val = min(c_low_val, c_high_val)

        ci_lows.append(c_low_val)
        ci_highs.append(c_high_val)

    # 1. Draw Confidence Interval ribbon across all timepoints (including negative/below-LOD samples)
    if num_pts >= 2:
        for i in range(num_pts - 1):
            ax_p.fill_between([x_coords[i], x_coords[i + 1]],
                              [ci_lows[i], ci_lows[i + 1]],
                              [ci_highs[i], ci_highs[i + 1]],
                              color=cfg["ci_color"], alpha=0.75, zorder=3)
    elif num_pts == 1:
        # Single timepoint: draw vertical confidence interval bar with caps
        ax_p.plot([x_coords[0], x_coords[0]], [ci_lows[0], ci_highs[0]],
                  color=cfg["color"], linewidth=2.5, alpha=0.6, zorder=3)
        cap_w = 0.08
        ax_p.plot([x_coords[0] - cap_w, x_coords[0] + cap_w], [ci_highs[0], ci_highs[0]],
                  color=cfg["color"], linewidth=2.0, alpha=0.7, zorder=3)
        ax_p.plot([x_coords[0] - cap_w, x_coords[0] + cap_w], [ci_lows[0], ci_lows[0]],
                  color=cfg["color"], linewidth=2.0, alpha=0.7, zorder=3)

    # 2. Draw trajectory lines connecting timepoints
    for i in range(num_pts - 1):
        x_seg = [x_coords[i], x_coords[i + 1]]
        y_seg = [all_plot_y[i], all_plot_y[i + 1]]
        if i in detected_indices and (i + 1) in detected_indices:
            # Both detected: solid trajectory line
            ax_p.plot(x_seg, y_seg, color=cfg["color"], linewidth=2.0, zorder=4)
        else:
            # Transition to or between below-LOD timepoints: dashed line to actual measured level
            ax_p.plot(x_seg, y_seg, color="#9ca3af", linestyle="--", linewidth=1.5, zorder=3.5)

    marker_sym = cfg.get("marker", "o")
    # 3. Draw detected markers (solid circular/diamond points with dark border)
    if detected_indices:
        ax_p.scatter([x_coords[i] for i in detected_indices], detected_y,
                     color=cfg["color"], edgecolors=cfg["edge_color"], s=60, marker=marker_sym, linewidths=1.5, zorder=5)

    # 4. Draw below-LOD markers (grey circular/diamond points at actual measured ctDNA level!)
    if below_lod_indices:
        ax_p.scatter([x_coords[i] for i in below_lod_indices], below_lod_y,
                     color="#6b7280", edgecolors="#374151", s=60, marker=marker_sym, linewidths=1.5, zorder=5)

    # X-axis ticks
    ax_p.set_xticks(x_coords)
    if is_bottom_track and x_labels:
        max_label_len = max(len(str(lbl).split("\n")[0]) for lbl in x_labels) if x_labels else 0
        needs_rotation = max_label_len > 4 or num_pts > 6
        rot_angle = 35 if needs_rotation else 0
        ha_align = "right" if needs_rotation else "center"
        ax_p.set_xticklabels(x_labels, fontsize=9.0, rotation=rot_angle, ha=ha_align)
        ax_p.set_xlabel("Sampling Date / Timepoint", fontsize=10.5, fontweight="bold")
    elif not is_bottom_track:
        ax_p.set_xticklabels([])

    # Legend matching clinical standard
    if show_legend:
        legend_elements = [
            Line2D([0], [0], color="black", linestyle=":", linewidth=1.5, label="WGS Baseline Timestamp"),
            Line2D([0], [0], color=cfg["color"], marker=marker_sym, markersize=7, linewidth=2.0, label=f"cfDNA Detected ({cfg['label_mod']})"),
            Patch(facecolor=cfg["ci_color"], edgecolor="none", alpha=0.9, label="Confidence Interval (95% CI)"),
            Line2D([0], [0], color="#6b7280", marker=marker_sym, markersize=7, linestyle="--", linewidth=1.5, label=f"Below Limit of Detection (Actual {cfg['label_mod']} Level)"),
            Line2D([0], [0], color="#94a3b8", linestyle="--", linewidth=1.2, label="Limit of Detection (LOD)"),
            Patch(facecolor="#e2e8f0", edgecolor="none", alpha=0.9, label="Below Limit of Detection Region"),
        ]
        ax_p.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 0.98),
                    ncol=2, frameon=True, framealpha=0.95, edgecolor="#cbd5e1", fontsize=8.5)


def generate_patient_mrd_trajectory_plot(patient_id, p_rec, plots_dir, formats=["png", "pdf"]):
    """
    Generates a publication-quality longitudinal MRD ctDNA trajectory plot for a patient,
    matching clinical cfDNA tracking standards for SNV modality.
    """
    if not HAS_MATPLOTLIB:
        return None, ""

    plasmas = sorted(p_rec.get("plasmas", []), key=timepoint_sort_key)
    if not plasmas:
        return None, ""

    # Determine X labels (sampling dates if available, otherwise timepoint labels)
    x_labels = []
    for idx, pl in enumerate(plasmas):
        dt = pl.get("sampling_date") or pl.get("date")
        tp = pl.get("timepoint")
        if tp and dt and tp != dt:
            lbl = f"{tp}\n{dt}"
        elif tp and tp not in {"Plasma ctDNA", "Plasma"}:
            lbl = tp
        elif dt:
            lbl = dt
        else:
            lbl = pl.get("clean_sample_name") or pl.get("sample_id") or f"T{idx}"
        x_labels.append(lbl)

    num_pts = len(plasmas)
    fig_w = min(36.0, max(10.5, num_pts * 1.1 + 2.0))
    fig_p, ax_p = plt.subplots(figsize=(fig_w, 5.8), dpi=200 if num_pts > 40 else 300)

    render_mrd_trajectory_on_axis(ax_p, plasmas, modality="SNV", patient_id=patient_id,
                                  x_labels=x_labels, show_legend=True, show_title=True,
                                  is_bottom_track=True)

    fig_p.tight_layout()

    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    out_png = plots_dir / f"{patient_id}_mrd_trajectory.png"
    out_pdf = plots_dir / f"{patient_id}_mrd_trajectory.pdf"

    fig_p.savefig(out_png, dpi=150, bbox_inches="tight")
    if "pdf" in formats:
        fig_p.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig_p)

    b64_str = ""
    try:
        import base64
        with open(out_png, "rb") as f_img:
            b64_str = f"data:image/png;base64,{base64.b64encode(f_img.read()).decode('ascii')}"
    except Exception:
        pass

    return str(out_png), b64_str


def generate_patient_mrd_multimodal_plot(patient_id, p_rec, plots_dir, formats=["png", "pdf"], embed_b64=False):
    """
    Generates a publication-quality 3-track longitudinal MRD ctDNA trajectory breakdown plot
    for a patient, rendering SNV, CNV, and LOH modalities with the EXACT SAME layout, styling,
    ribbons, and LOD shading as the standard trajectory plots.
    """
    if not HAS_MATPLOTLIB:
        return None, ""

    plasmas = sorted(p_rec.get("plasmas", []), key=timepoint_sort_key)
    if not plasmas:
        return None, ""

    x_labels = []
    for idx, pl in enumerate(plasmas):
        dt = pl.get("sampling_date") or pl.get("date")
        tp = pl.get("timepoint")
        if tp and dt and tp != dt:
            lbl = f"{tp}\n{dt}"
        elif tp and tp not in {"Plasma ctDNA", "Plasma"}:
            lbl = tp
        elif dt:
            lbl = dt
        else:
            lbl = pl.get("clean_sample_name") or pl.get("sample_id") or f"T{idx}"
        x_labels.append(lbl)

    num_pts = len(plasmas)
    fig_w = min(36.0, max(11.0, num_pts * 1.1 + 2.5))
    patient_dpi = 200 if num_pts > 40 else 300
    has_ichor_plasmas = any(pl.get("ichor_tf") is not None for pl in plasmas)
    if has_ichor_plasmas:
        fig_multi, (ax_s, ax_c, ax_l, ax_i) = plt.subplots(4, 1, figsize=(fig_w, 17.5), sharex=True, dpi=patient_dpi)
        render_mrd_trajectory_on_axis(ax_s, plasmas, modality="SNV", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=False)
        render_mrd_trajectory_on_axis(ax_c, plasmas, modality="CNV", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=False)
        render_mrd_trajectory_on_axis(ax_l, plasmas, modality="LOH", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=False)
        render_mrd_trajectory_on_axis(ax_i, plasmas, modality="ICHOR", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=True)
        fig_multi.suptitle(f"{patient_id} — Multi-Modal ctDNA Trajectory Breakdown (SNV / CNV / LOH / 1X ichorCNA)",
                           fontsize=13.5, fontweight="bold", y=0.995)
    else:
        fig_multi, (ax_s, ax_c, ax_l) = plt.subplots(3, 1, figsize=(fig_w, 14.0), sharex=True, dpi=patient_dpi)
        render_mrd_trajectory_on_axis(ax_s, plasmas, modality="SNV", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=False)
        render_mrd_trajectory_on_axis(ax_c, plasmas, modality="CNV", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=False)
        render_mrd_trajectory_on_axis(ax_l, plasmas, modality="LOH", patient_id=patient_id,
                                      x_labels=x_labels, show_legend=True, show_title=True,
                                      is_bottom_track=True)
        fig_multi.suptitle(f"{patient_id} — Multi-Modal ctDNA Trajectory Breakdown (SNV / CNV / LOH)",
                           fontsize=13.5, fontweight="bold", y=0.995)
    fig_multi.tight_layout()

    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    out_png = plots_dir / f"{patient_id}_mrd_multimodal.png"
    out_pdf = plots_dir / f"{patient_id}_mrd_multimodal.pdf"

    fig_multi.savefig(out_png, dpi=150, bbox_inches="tight")
    if "pdf" in formats:
        fig_multi.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig_multi)

    b64_str = ""
    if embed_b64:
        try:
            import base64
            with open(out_png, "rb") as f_img:
                b64_str = f"data:image/png;base64,{base64.b64encode(f_img.read()).decode('ascii')}"
        except Exception:
            pass

    return str(out_png), b64_str


def _render_patient_plots_worker(args_tuple):
    """
    Worker function to render patient trajectory and multimodal plots in parallel.
    args_tuple: (patient_id, p_rec, plots_dir_str, formats[, embed_all_plots])
    Returns: (patient_id, traj_png, traj_b64, multi_png, multi_b64)
    """
    if len(args_tuple) >= 5:
        pid, p_rec, plots_dir_str, formats, embed_all_plots = args_tuple[:5]
    else:
        pid, p_rec, plots_dir_str, formats = args_tuple[:4]
        embed_all_plots = False

    plots_dir = Path(plots_dir_str)
    traj_png, traj_b64 = None, ""
    multi_png, multi_b64 = None, ""
    try:
        traj_png, traj_b64 = generate_patient_mrd_trajectory_plot(pid, p_rec, plots_dir, formats=formats)
    except Exception as e:
        print(f"Warning: Failed generating trajectory plot for patient {pid}: {e}", file=sys.stderr)

    try:
        multi_png, multi_b64 = generate_patient_mrd_multimodal_plot(pid, p_rec, plots_dir, formats=formats, embed_b64=embed_all_plots)
    except Exception as e:
        print(f"Warning: Failed generating multi-modal plot for patient {pid}: {e}", file=sys.stderr)

    return pid, traj_png, traj_b64, multi_png, multi_b64


def _render_plasma_plot_worker(args_tuple):
    """
    Worker function to render a single plasma multimodal decomposition plot in parallel.
    args_tuple: (pid, sid, pl, tumor_rec, out_img_str, plot_fname[, embed_all_plots])
    Returns: (pid, sid, rel_path, b64_str, ok, err_msg)
    """
    if len(args_tuple) >= 7:
        pid, sid, pl, tumor_rec, out_img_str, plot_fname, embed_all_plots = args_tuple[:7]
    else:
        pid, sid, pl, tumor_rec, out_img_str, plot_fname = args_tuple[:6]
        embed_all_plots = False

    try:
        out_img = Path(out_img_str)
        ok = generate_plasma_wisp_multimodal_plot(pl, tumor_rec, out_img_str)
        b64_str = ""
        if embed_all_plots and ok and out_img.exists():
            try:
                import base64
                with open(out_img, "rb") as f_img:
                    b64_str = f"data:image/png;base64,{base64.b64encode(f_img.read()).decode('ascii')}"
            except Exception:
                pass
        return pid, sid, f"plasma_plots/{plot_fname}", b64_str, True, None
    except Exception as e:
        return pid, sid, "", "", False, str(e)


def generate_static_plots(cohort, oncoprint_data, output_dir, formats=["png", "pdf"], skip_plasma_plots=False, threads=8, embed_all_plots=False):
    """
    Generate cohort-level multi-panel publication plot and per-patient longitudinal
    trajectories with log-scale y-axis and greyed-out Limit of Detection (LOD).
    Parallelizes per-patient and per-plasma figure rendering across multiple workers.
    """
    if not HAS_MATPLOTLIB:
        print("Notice: matplotlib not installed; skipping static image generation.", file=sys.stderr)
        return

    out_path = Path(output_dir)
    plots_dir = out_path / "patient_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------------------------
    # Plot 1: Multi-panel Cohort Summary Figure
    # --------------------------------------------------------------------------
    fig = plt.figure(figsize=(18, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)

    # Panel A: Purity vs Ploidy
    ax1 = fig.add_subplot(gs[0, 0])
    tumors = cohort.get("all_tumors", [])
    valid_tumors = [t for t in tumors if t.get("purity_pct") is not None and t.get("ploidy") is not None]
    if valid_tumors:
        purities = [t["purity_pct"] for t in valid_tumors]
        ploidies = [t["ploidy"] for t in valid_tumors]
        t_labels = [t["tumor_id"] for t in valid_tumors]
        scatter = ax1.scatter(ploidies, purities, c="#1f77b4", s=90, edgecolors="#0f3d61", alpha=0.85, zorder=3)
        ax1.set_xlabel("Tumor Ploidy", fontsize=11, fontweight="bold")
        ax1.set_ylabel("Tumor Purity (%)", fontsize=11, fontweight="bold")
        ax1.set_title("A. Primary Tumor Purity vs. Ploidy", fontsize=12, fontweight="bold", loc="left")
        ax1.grid(True, linestyle="--", alpha=0.4, zorder=1)
        ax1.set_ylim(0, 105)
        for i, lbl in enumerate(t_labels):
            if i < 8:  # annotate top few to avoid clutter
                ax1.annotate(lbl, (ploidies[i], purities[i]), textcoords="offset points", xytext=(5, 5), fontsize=8, alpha=0.8)
    else:
        ax1.text(0.5, 0.5, "No primary tumors detected", ha="center", va="center", transform=ax1.transAxes)

    # Panel B: Top Mutated / Driver Genes Frequency
    ax2 = fig.add_subplot(gs[0, 1])
    top_genes = oncoprint_data["genes"][:12]
    if top_genes:
        counts = [oncoprint_data["summary"][g]["total_altered"] for g in top_genes]
        pcts = [oncoprint_data["summary"][g]["pct"] for g in top_genes]
        y_pos = range(len(top_genes))

        ax2.barh(y_pos, pcts, color="#2a9d8f", edgecolor="#1a5f57", alpha=0.85, height=0.65, zorder=3)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(top_genes, fontsize=10, fontweight="bold")
        ax2.invert_yaxis()
        ax2.set_xlabel("% of Analyzed Tumors Altered", fontsize=11, fontweight="bold")
        ax2.set_title("B. Frequent Genomic Alterations & Drivers", fontsize=12, fontweight="bold", loc="left")
        ax2.grid(True, linestyle="--", alpha=0.4, axis="x", zorder=1)
        for i, pct in enumerate(pcts):
            ax2.text(pct + 0.8, i, f"{pct:.1f}% ({counts[i]})", va="center", fontsize=9, fontweight="bold", color="#264653")
    else:
        ax2.text(0.5, 0.5, "No genomic alterations detected", ha="center", va="center", transform=ax2.transAxes)

    # Panel C: TMB Distribution by MSI Status
    ax3 = fig.add_subplot(gs[1, 0])
    tmbs = [t["tmb"] for t in tumors if t["tmb"] is not None and t["tmb"] > 0]
    if tmbs:
        mss_tmbs = [t["tmb"] for t in tumors if t["tmb"] is not None and t["tmb"] > 0 and t["msi_status"] == "MSS"]
        msi_tmbs = [t["tmb"] for t in tumors if t["tmb"] is not None and t["tmb"] > 0 and t["msi_status"] != "MSS"]

        boxes = []
        labels = []
        if mss_tmbs:
            boxes.append(mss_tmbs)
            labels.append(f"MSS (n={len(mss_tmbs)})")
        if msi_tmbs:
            boxes.append(msi_tmbs)
            labels.append(f"MSI-H (n={len(msi_tmbs)})")

        if boxes:
            try:
                ax3.boxplot(boxes, tick_labels=labels, patch_artist=True, boxprops=dict(facecolor="#e76f51", alpha=0.7), medianprops=dict(color="black", linewidth=1.5))
            except TypeError:
                ax3.boxplot(boxes, labels=labels, patch_artist=True, boxprops=dict(facecolor="#e76f51", alpha=0.7), medianprops=dict(color="black", linewidth=1.5))
        ax3.set_yscale("log")
        ax3.set_ylabel("TMB (Mutations / Mb - Log Scale)", fontsize=11, fontweight="bold")
        ax3.set_title("C. Tumor Mutational Burden (TMB)", fontsize=12, fontweight="bold", loc="left")
        ax3.grid(True, linestyle="--", alpha=0.4, zorder=1)
    else:
        ax3.text(0.5, 0.5, "No TMB metrics detected", ha="center", va="center", transform=ax3.transAxes)

    # Panel D: Cohort Plasma ctDNA Purity by Modality (SNV, LOH, CNV) & Sample LODs
    ax4 = fig.add_subplot(gs[1, 1])
    if cohort["has_mrd_data"]:
        ax4.set_yscale("log")
        y_min = 0.005  # 0.005% log floor
        y_max = 100.0  # 100% max

        def mrd_sort_key(pl):
            if pl.get("snv_mrd"):
                return (0, -(safe_float(pl.get("snv_purity_pct"), 0.0)))
            if pl.get("cnv_mrd"):
                return (1, -(safe_float(pl.get("cnv_purity_pct"), 0.0)))
            if pl.get("loh_mrd"):
                return (2, -(safe_float(pl.get("loh_purity_pct"), 0.0)))
            return (3, pl.get("patient_id", ""))

        plasmas = sorted(cohort["all_plasmas"], key=mrd_sort_key)
        all_lods = [safe_float(pl.get("snv_lod_pct"), 0.0) for pl in plasmas if safe_float(pl.get("snv_lod_pct"), 0.0) > 0]
        cohort_lod = (sum(all_lods) / len(all_lods)) if all_lods else 0.03

        # Greyed-out sub-analytical LOD region at bottom
        ax4.axhspan(y_min, cohort_lod, color="#ced4da", alpha=0.5, zorder=1, label=f"Sub-LOD Floor (< {cohort_lod:.3f}%)")
        ax4.axhline(cohort_lod, color="#6c757d", linestyle="--", linewidth=1.0, zorder=2)

        x_indices = list(range(len(plasmas)))
        x_labels = []
        for pl in plasmas:
            p_id = pl.get("patient_id", "")
            s_name = pl.get("clean_sample_name") or pl.get("sample_id", "")
            tp = pl.get("timepoint", "")
            if s_name and s_name != p_id:
                lbl = f"{p_id}\n{s_name}"
            else:
                lbl = f"{p_id}"
            if tp and tp not in lbl:
                lbl += f" ({tp})"
            x_labels.append(lbl)

        for idx, pl in enumerate(plasmas):
            # 1. SNV (offset -0.2)
            snv_p = safe_float(pl.get("snv_purity_pct"), 0.0)
            snv_l = safe_float(pl.get("snv_lod_pct"), 0.0)
            snv_y = max(snv_p, y_min * 1.2) if snv_p > 0 else (y_min * 1.1)
            snv_face = "#0ea5e9" if pl.get("snv_mrd") else "none"
            ax4.scatter([idx - 0.2], [snv_y], color=snv_face, edgecolors="#0369a1", s=75, marker="o", linewidths=1.5, zorder=4)
            if snv_l > 0:
                ax4.scatter([idx - 0.2], [max(snv_l, y_min)], color="#0284c7", marker="_", s=100, linewidths=2.0, zorder=3)

            # 2. LOH (offset 0.0)
            loh_p = safe_float(pl.get("loh_purity_pct"), 0.0)
            loh_l = safe_float(pl.get("loh_lod_pct"), 0.0)
            loh_y = max(loh_p, y_min * 1.2) if loh_p > 0 else (y_min * 1.1)
            loh_face = "#f59e0b" if pl.get("loh_mrd") else "none"
            ax4.scatter([idx], [loh_y], color=loh_face, edgecolors="#b45309", s=70, marker="D", linewidths=1.5, zorder=4)
            if loh_l > 0:
                ax4.scatter([idx], [max(loh_l, y_min)], color="#b45309", marker="_", s=100, linewidths=2.0, zorder=3)

            # 3. CNV (offset +0.2)
            cnv_p = safe_float(pl.get("cnv_purity_pct"), 0.0)
            cnv_l = safe_float(pl.get("cnv_lod_pct"), 0.0)
            cnv_y = max(cnv_p, y_min * 1.2) if cnv_p > 0 else (y_min * 1.1)
            cnv_face = "#a855f7" if pl.get("cnv_mrd") else "none"
            ax4.scatter([idx + 0.2], [cnv_y], color=cnv_face, edgecolors="#7e22ce", s=70, marker="s", linewidths=1.5, zorder=4)
            if cnv_l > 0:
                ax4.scatter([idx + 0.2], [max(cnv_l, y_min)], color="#7e22ce", marker="_", s=100, linewidths=2.0, zorder=3)

        ax4.set_xticks(x_indices)
        ax4.set_xticklabels(x_labels, fontsize=8.0, rotation=35, ha="right")
        ax4.set_ylim(y_min, y_max)
        ax4.set_ylabel("Plasma ctDNA Purity (%) - Log Scale", fontsize=11, fontweight="bold")
        ax4.set_title("D. Cohort Plasma ctDNA (Sorted SNV > CNV > LOH) & LODs", fontsize=12, fontweight="bold", loc="left")
        ax4.grid(True, linestyle="--", alpha=0.35, which="both", zorder=1)

        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#0ea5e9", markeredgecolor="#0369a1", markersize=8, label="SNV ctDNA"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#a855f7", markeredgecolor="#7e22ce", markersize=7, label="CNV ctDNA"),
            Line2D([0], [0], marker="D", color="w", markerfacecolor="#f59e0b", markeredgecolor="#b45309", markersize=7, label="LOH ctDNA"),
            Line2D([0], [0], marker="_", color="#475569", markersize=10, markeredgewidth=2, linestyle="None", label="Sample-specific LOD"),
            Patch(facecolor="#ced4da", edgecolor="#6c757d", alpha=0.5, label=f"Sub-LOD Floor (< {cohort_lod:.3f}%)"),
        ]
        ax4.legend(handles=legend_elements, loc="upper right", fontsize=8, framealpha=0.9)
    else:
        ax4.text(0.5, 0.5, "No plasma MRD (WISP) data available in cohort", ha="center", va="center", transform=ax4.transAxes)

    plt.suptitle(f"Oncoanalyser WGTS Multi-Omics Cohort Summary: {cohort['cohort_name']}", fontsize=15, fontweight="bold", y=0.98)

    for fmt in formats:
        out_fig = out_path / f"cohort_wgts_summary.{fmt}"
        plt.savefig(out_fig, format=fmt, bbox_inches="tight")
        print(f"Cohort plot saved: {out_fig}")
    plt.close(fig)

    # Dedicated Plot: Stacked Multi-Modal MRD Publication Figure (SNV / CNV / LOH / 1X ichorCNA)
    if cohort["has_mrd_data"] and cohort["all_plasmas"]:
        # Cap figure width to prevent matplotlib Agg backend 65536 pixel limit overflow in large cohorts
        fig_w = min(36.0, max(14.0, len(plasmas) * 0.16 + 4.0))
        mrd_dpi = 150 if len(plasmas) > 40 else 200
        has_ichor_cohort = any(pl.get("ichor_tf") is not None for pl in plasmas)
        y_min = 0.0005
        y_max = 105.0

        if has_ichor_cohort:
            fig_mrd, (ax_s, ax_c, ax_l, ax_i) = plt.subplots(4, 1, figsize=(fig_w, 15.0), sharex=True, dpi=mrd_dpi)
            track_configs = [
                (ax_s, "snv_purity_pct", "snv_mrd", "snv_lod_pct", "#ef4444", "#b91c1c", "Track 1: Cohort SNV ctDNA Purity (%) & Sample-Specific SNV LOD", "SNV", "o", "cfDNA Detected (SNV)"),
                (ax_c, "cnv_purity_pct", "cnv_mrd", "cnv_lod_pct", "#a855f7", "#7e22ce", "Track 2: Cohort Copy Number (CNV) ctDNA Purity (%) & Sample-Specific CNV LOD", "CNV", "o", "cfDNA Detected (CNV)"),
                (ax_l, "loh_purity_pct", "loh_mrd", "loh_lod_pct", "#f59e0b", "#b45309", "Track 3: Cohort Loss of Heterozygosity (LOH) ctDNA Purity (%) & Sample-Specific LOH LOD", "LOH", "o", "cfDNA Detected (LOH)"),
                (ax_i, "ichor_purity_pct", "ichor_mrd", "ichor_lod_pct", "#d97706", "#92400e", "Track 4: Cohort 1X ichorCNA Tumor Fraction (%) [LoD Floor: 3.0%]", "1X ichorCNA", "D", "Detected (TF >= 3%)"),
            ]
            bottom_ax = ax_i
        else:
            fig_mrd, (ax_s, ax_c, ax_l) = plt.subplots(3, 1, figsize=(fig_w, 11.5), sharex=True, dpi=mrd_dpi)
            track_configs = [
                (ax_s, "snv_purity_pct", "snv_mrd", "snv_lod_pct", "#ef4444", "#b91c1c", "Track 1: Cohort SNV ctDNA Purity (%) & Sample-Specific SNV LOD", "SNV", "o", "cfDNA Detected (SNV)"),
                (ax_c, "cnv_purity_pct", "cnv_mrd", "cnv_lod_pct", "#a855f7", "#7e22ce", "Track 2: Cohort Copy Number (CNV) ctDNA Purity (%) & Sample-Specific CNV LOD", "CNV", "o", "cfDNA Detected (CNV)"),
                (ax_l, "loh_purity_pct", "loh_mrd", "loh_lod_pct", "#f59e0b", "#b45309", "Track 3: Cohort Loss of Heterozygosity (LOH) ctDNA Purity (%) & Sample-Specific LOH LOD", "LOH", "o", "cfDNA Detected (LOH)"),
            ]
            bottom_ax = ax_l

        for ax, metric_k, mrd_k, lod_k, color, edge_color, title, mod_lbl, marker, det_lbl in track_configs:
            ax.set_yscale("log")
            ax.set_ylim(y_min, y_max)
            ax.set_yticks([0.001, 0.01, 0.1, 1.0, 10.0, 100.0])
            ax.set_yticklabels(["0.001%", "0.01%", "0.1%", "1%", "10%", "100%"], fontsize=9.0)
            ax.grid(True, which="both", linestyle="-", color="#e5e7eb", linewidth=0.8, alpha=0.8, zorder=1)
            ax.set_ylabel("Tumor Fraction (%)", fontsize=10, fontweight="bold")
            ax.set_title(title, fontsize=11, fontweight="bold", loc="left")

            if mod_lbl == "1X ichorCNA":
                cohort_lods = [3.0 for _ in plasmas]
                below_lod_lbl = "Below 3% LoD (Actual TF Level)"
                lod_leg_lbl = "Assay LoD Floor (3.0%)"
            else:
                cohort_lods = [max(pl.get(lod_k) or 0.002, y_min) for pl in plasmas]
                below_lod_lbl = f"Below LOD (Actual {mod_lbl} Level)"
                lod_leg_lbl = "Limit of Detection (LOD)"

            # Continuous shaded grey region below LOD
            ax.fill_between(x_indices, y_min, cohort_lods, color="#e2e8f0", alpha=0.9, zorder=2)
            # Continuous dashed cutoff line (NO square points!)
            ax.plot(x_indices, cohort_lods, color="#94a3b8", linestyle="--", linewidth=1.2, zorder=2.5)

            for idx, pl in enumerate(plasmas):
                val = pl.get(metric_k)
                is_det = pl.get(mrd_k)
                if val is not None and val > 0:
                    y_pt = max(val, y_min * 1.15)
                else:
                    y_pt = (y_min * 1.15)

                if is_det:
                    ax.scatter([idx], [y_pt], color=color, edgecolors=edge_color, s=65, marker=marker, linewidths=1.5, zorder=4)
                else:
                    ax.scatter([idx], [y_pt], color="#6b7280" if val is not None else "#9ca3af", edgecolors="#374151", s=65, marker=marker, linewidths=1.5, zorder=4)

            leg_elems = [
                Line2D([0], [0], color=color, marker=marker, markersize=7, linewidth=0, label=det_lbl),
                Line2D([0], [0], color="#6b7280", marker=marker, markersize=7, linewidth=0, label=below_lod_lbl),
                Line2D([0], [0], color="#94a3b8", linestyle="--", linewidth=1.2, label=lod_leg_lbl),
                Patch(facecolor="#e2e8f0", edgecolor="none", alpha=0.9, label="Below Limit of Detection Region"),
            ]
            ax.legend(handles=leg_elems, loc="upper right", fontsize=8.5, framealpha=0.95, edgecolor="#cbd5e1")

        bottom_ax.set_xticks(x_indices)
        tick_font_size = 8.5 if len(plasmas) <= 40 else (7.0 if len(plasmas) <= 80 else 5.5)
        if len(plasmas) > 120:
            step = max(1, len(plasmas) // 60)
            thinned_labels = [lbl if (i % step == 0 or i == len(plasmas) - 1) else "" for i, lbl in enumerate(x_labels)]
            bottom_ax.set_xticklabels(thinned_labels, fontsize=6.5, rotation=45, ha="right")
        else:
            bottom_ax.set_xticklabels(x_labels, fontsize=tick_font_size, rotation=35, ha="right")

        bottom_ax.set_xlabel("Analyzed Plasma Samples (Sorted: SNV MRD+ > CNV MRD+ > LOH MRD+ > Neg)", fontsize=10, fontweight="bold")
        fig_mrd.suptitle(f"Cohort Plasma ctDNA Multi-Modal Profile: {cohort['cohort_name']}", fontsize=13, fontweight="bold", y=0.98)
        fig_mrd.tight_layout()

        for fmt in formats:
            out_mrd_fig = out_path / f"cohort_wgts_mrd_multimodal.{fmt}"
            fig_mrd.savefig(out_mrd_fig, format=fmt, dpi=mrd_dpi, bbox_inches="tight")
            print(f"Multi-modal MRD plot saved: {out_mrd_fig}")

        cohort["cohort_mrd_multimodal_plot"] = "cohort_wgts_mrd_multimodal.png"
        out_mrd_png = out_path / "cohort_wgts_mrd_multimodal.png"
        if out_mrd_png.is_file():
            try:
                import base64
                with open(out_mrd_png, "rb") as f_img:
                    cohort["cohort_mrd_multimodal_plot_b64"] = f"data:image/png;base64,{base64.b64encode(f_img.read()).decode('ascii')}"
            except Exception:
                pass

        plt.close(fig_mrd)

    # --------------------------------------------------------------------------
    # Plot 2: Per-Patient Log-Scale Plots (SNV, LOH, CNV, WISP + Sample LODs)
    # --------------------------------------------------------------------------
    if cohort.get("has_mrd_data"):
        patient_tasks = []
        for p_id, p_rec in cohort["patients"].items():
            plasmas = p_rec.get("plasmas", [])
            if not plasmas:
                continue
            patient_tasks.append((p_id, p_rec, str(plots_dir), formats, embed_all_plots))

        if patient_tasks:
            worker_count = min(threads, len(patient_tasks)) if threads > 1 else 1
            print(f"Generating MRD plots for {len(patient_tasks)} patient(s) (using {worker_count} worker(s))...")

            def apply_patient_result(pid, traj_png, traj_b64, multi_png, multi_b64):
                p_rec = cohort["patients"][pid]
                plasmas = p_rec.get("plasmas", [])
                if traj_png:
                    rel_plot = f"patient_plots/{Path(traj_png).name}"
                    p_rec["trajectory_plot"] = rel_plot
                    p_rec["trajectory_plot_b64"] = traj_b64
                    for pl in plasmas:
                        pl["trajectory_plot"] = rel_plot
                    print(f"Patient MRD trajectory plot saved: {traj_png}")
                if multi_png:
                    rel_multi = f"patient_plots/{Path(multi_png).name}"
                    p_rec["multimodal_plot"] = rel_multi
                    p_rec["multimodal_plot_b64"] = multi_b64
                    for pl in plasmas:
                        pl["multimodal_plot"] = rel_multi
                    print(f"Patient MRD multi-modal plot saved: {multi_png}")

            if worker_count > 1:
                with ProcessPoolExecutor(max_workers=worker_count) as executor:
                    for res in executor.map(_render_patient_plots_worker, patient_tasks):
                        apply_patient_result(*res)
            else:
                for task in patient_tasks:
                    apply_patient_result(*_render_patient_plots_worker(task))

    # --------------------------------------------------------------------------
    # Plot 3: Per-Plasma ctDNA Multi-Modal Decomposition (CNV, LOH, SNV)
    # --------------------------------------------------------------------------
    if not skip_plasma_plots and cohort.get("has_mrd_data") and HAS_MATPLOTLIB and HAS_PANDAS:
        plasma_plots_dir = out_path / "plasma_plots"
        plasma_plots_dir.mkdir(parents=True, exist_ok=True)

        plasma_tasks = []
        for pl in cohort.get("all_plasmas", []):
            pid = pl.get("patient_id")
            sid = pl.get("sample_id")
            prec = cohort.get("patients", {}).get(pid, {})
            tumors = prec.get("tumors", [])
            tumor_rec = None
            for t in tumors:
                if t.get("tumor_id") == pl.get("tumor_id"):
                    tumor_rec = t
                    break
            if not tumor_rec and tumors:
                tumor_rec = tumors[0]

            plot_fname = f"{pid}_{sid}_wisp_multimodal.png"
            out_img = str(plasma_plots_dir / plot_fname)
            plasma_tasks.append((pid, sid, pl, tumor_rec, out_img, plot_fname, embed_all_plots))

        if plasma_tasks:
            worker_count = min(threads, len(plasma_tasks)) if threads > 1 else 1
            print(f"Generating {len(plasma_tasks)} plasma ctDNA decomposition plot(s) (using {worker_count} worker(s))...")

            plasmas_by_key = {(pl.get("patient_id"), pl.get("sample_id")): pl for pl in cohort.get("all_plasmas", [])}

            def apply_plasma_result(pid, sid, rel_path, b64_str, ok, err_msg):
                if ok:
                    target_pl = plasmas_by_key.get((pid, sid))
                    if target_pl:
                        target_pl["multimodal_plot"] = rel_path
                        target_pl["multimodal_plot_b64"] = b64_str
                    print(f"Plasma multi-modal plot saved: {plasma_plots_dir / f'{pid}_{sid}_wisp_multimodal.png'}")
                elif err_msg:
                    print(f"Warning: Failed to generate multi-modal plot for {pid} {sid}: {err_msg}", file=sys.stderr)

            if worker_count > 1:
                with ProcessPoolExecutor(max_workers=worker_count) as executor:
                    for res in executor.map(_render_plasma_plot_worker, plasma_tasks):
                        apply_plasma_result(*res)
            else:
                for task in plasma_tasks:
                    apply_plasma_result(*_render_plasma_plot_worker(task))


# ==============================================================================
# 7. INTERACTIVE HTML DASHBOARD (LANDING PAGE + QC + ONCOPRINT + LOG-SCALE MRD)
# ==============================================================================

def generate_interactive_dashboard(cohort, oncoprint_data, output_dir, embed_all_plots=False):
    """
    Renders a modern, single-page interactive application containing:
      - Landing Page: KPIs & Patient/Sample Master Table
      - Tab 1: QC Plots (Purity, Depth, TMB, MSI, ichorCNA)
      - Tab 2: Genomic Summary with interactive Oncoprint & live gene filter
      - Tab 3: MRD (Log-scaled longitudinal ctDNA trajectories + greyed-out LOD region)
    """
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    html_file = out_path / "cohort_wgts_dashboard.html"

    # Prepare pruned, lightweight copies of sample records to prevent HTML file ballooning
    clean_all_tumors = []
    for t in cohort.get("all_tumors", []):
        t_c = dict(t)
        # Deduplicate per-sample ichor b64 (available in ichor_reports_b64)
        t_c["ichor_pdf_b64"] = ""
        t_c["ichor_standard_pdf_b64"] = ""
        t_c["ichor_short_pdf_b64"] = ""
        clean_all_tumors.append(t_c)

    clean_all_plasmas = []
    for pl in cohort.get("all_plasmas", []):
        pl_c = dict(pl)
        if not embed_all_plots:
            pl_c["multimodal_plot_b64"] = ""
        pl_c["ichor_pdf_b64"] = ""
        pl_c["ichor_standard_pdf_b64"] = ""
        pl_c["ichor_short_pdf_b64"] = ""
        clean_all_plasmas.append(pl_c)

    payload_patients = {}
    for pid, prec in cohort.get("patients", {}).items():
        p_tumors = []
        for t in prec.get("tumors", []):
            t_c = dict(t)
            t_c["ichor_pdf_b64"] = ""
            t_c["ichor_standard_pdf_b64"] = ""
            t_c["ichor_short_pdf_b64"] = ""
            p_tumors.append(t_c)

        p_plasmas = []
        for pl in prec.get("plasmas", []):
            pl_c = dict(pl)
            if not embed_all_plots:
                pl_c["multimodal_plot_b64"] = ""
            pl_c["ichor_pdf_b64"] = ""
            pl_c["ichor_standard_pdf_b64"] = ""
            pl_c["ichor_short_pdf_b64"] = ""
            p_plasmas.append(pl_c)

        payload_patients[pid] = {
            "patient_id": pid,
            "tumors": p_tumors,
            "plasmas": p_plasmas,
            "ichor": prec["ichor"],
            "trajectory_plot": prec.get("trajectory_plot", ""),
            "trajectory_plot_b64": prec.get("trajectory_plot_b64", ""),
            "multimodal_plot": prec.get("multimodal_plot", ""),
            "multimodal_plot_b64": prec.get("multimodal_plot_b64", "") if embed_all_plots else "",
        }

    # Prepare JSON payloads
    payload = {
        "cohort_name": cohort["cohort_name"],
        "total_patients": cohort["total_analyzed_patients"],
        "total_tumors": len(clean_all_tumors),
        "total_plasmas": len(clean_all_plasmas),
        "has_mrd": cohort["has_mrd_data"],
        "cohort_mrd_multimodal_plot": cohort.get("cohort_mrd_multimodal_plot", ""),
        "cohort_mrd_multimodal_plot_b64": cohort.get("cohort_mrd_multimodal_plot_b64", ""),
        "orange_reports_b64": cohort.get("orange_reports_b64", {}),
        "ichor_reports_b64": cohort.get("ichor_reports_b64", {}),
        "all_tumors": clean_all_tumors,
        "all_plasmas": clean_all_plasmas,
        "all_orange_driver_genes": cohort.get("all_orange_driver_genes", []),
        "fingerprints": cohort.get("fingerprints", {}),
        "oncoprint": {
            "genes": oncoprint_data["genes"],
            "tumors": oncoprint_data["tumors"],
            "matrix": oncoprint_data["matrix"],
            "cell_details": oncoprint_data.get("cell_details", {}),
            "summary": oncoprint_data["summary"],
            "patient_blocks": oncoprint_data.get("patient_blocks", []),
            "patient_boundaries": oncoprint_data.get("patient_boundaries", []),
            "sample_metadata": oncoprint_data.get("sample_metadata", {}),
        },
        "patients": payload_patients
    }

    def _sanitize_for_json(obj):
        if isinstance(obj, dict):
            res = {}
            for k, v in obj.items():
                if isinstance(k, tuple):
                    str_k = "_".join(str(part) for part in k if part is not None)
                else:
                    str_k = str(k)
                res[str_k] = _sanitize_for_json(v)
            return res
        elif isinstance(obj, (list, tuple)):
            return [_sanitize_for_json(x) for x in obj]
        elif isinstance(obj, set):
            return sorted([_sanitize_for_json(x) for x in obj], key=lambda x: str(x))
        elif isinstance(obj, (datetime, date)) or (hasattr(obj, "isoformat") and callable(obj.isoformat)):
            try:
                return obj.isoformat()
            except Exception:
                return str(obj)
        elif isinstance(obj, Path):
            return str(obj)
        elif hasattr(obj, "item"):
            try:
                return obj.item()
            except Exception:
                return str(obj)
        return obj

    clean_payload = _sanitize_for_json(payload)
    payload_json = json.dumps(clean_payload, separators=(',', ':'), default=str)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{cohort['cohort_name']} - WGTS Multi-Omics & MRD Dashboard</title>
  <!-- Google Fonts: Inter & Outfit -->
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=Outfit:wght@400;600;700;800&display=swap" rel="stylesheet">
  <!-- Chart.js for responsive interactive plotting -->
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
  <!-- PDF.js for in-browser PDF canvas rendering (compatible with OneDrive / sandboxed frames) -->
  <script src="https://cdnjs.cloudflare.com/ajax/libs/pdf.js/3.11.174/pdf.min.js"></script>

  <style>
    :root {{
      --bg-primary: #f8fafc;
      --bg-secondary: #ffffff;
      --bg-card: #ffffff;
      --bg-subtle: #f1f5f9;
      --bg-glass: rgba(15, 23, 42, 0.03);
      --border-color: #e2e8f0;
      --border-strong: #cbd5e1;
      --text-primary: #0f172a;
      --text-secondary: #475569;
      --text-muted: #94a3b8;
      --accent-blue: #0284c7;
      --accent-teal: #0d9488;
      --accent-crimson: #dc2626;
      --accent-amber: #d97706;
      --accent-purple: #7c3aed;
      --badge-pos: #dc2626;
      --badge-neg: #059669;
      --badge-pos-bg: #fee2e2;
      --badge-neg-bg: #dcfce7;
      --font-body: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      --font-head: 'Outfit', 'Inter', -apple-system, sans-serif;
    }}

    * {{ box-sizing: border-box; margin: 0; padding: 0; }}

    body {{
      font-family: var(--font-body);
      background-color: var(--bg-primary);
      color: var(--text-primary);
      line-height: 1.5;
      padding: 0;
      margin: 0;
      min-height: 100vh;
      -webkit-font-smoothing: antialiased;
    }}

    /* Header Nav */
    header {{
      background: var(--bg-secondary);
      border-bottom: 1px solid var(--border-color);
      padding: 1rem 2rem;
      display: flex;
      justify-content: space-between;
      align-items: center;
      position: sticky;
      top: 0;
      z-index: 100;
      box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    }}

    .logo-group h1 {{
      font-family: var(--font-head);
      font-size: 1.35rem;
      font-weight: 700;
      color: var(--text-primary);
      display: flex;
      align-items: center;
      gap: 0.5rem;
      letter-spacing: -0.01em;
    }}

    .logo-group p {{
      font-size: 0.8rem;
      color: var(--text-secondary);
      margin-top: 2px;
    }}

    .badge {{
      display: inline-flex;
      align-items: center;
      gap: 4px;
      padding: 0.2rem 0.55rem;
      border-radius: 4px;
      font-size: 0.72rem;
      font-weight: 600;
      letter-spacing: 0.02em;
      white-space: nowrap;
    }}

    .badge-primary {{ background: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd; }}
    .badge-mrd-pos {{ background: #fee2e2; color: #991b1b; border: 1px solid #fecaca; }}
    .badge-mrd-neg {{ background: #dcfce7; color: #166534; border: 1px solid #bbf7d0; }}
    .badge-normal {{ background: #f1f5f9; color: #334155; border: 1px solid #e2e8f0; }}
    .badge-cuppa {{ background: #fef3c7; color: #92400e; border: 1px solid #fde68a; }}

    /* Tab Navigation: Quiet neutral segmented pill */
    .tabs-nav {{
      display: flex;
      gap: 3px;
      background: var(--bg-subtle);
      padding: 4px;
      border-radius: 8px;
      border: 1px solid var(--border-color);
    }}

    .tab-btn {{
      background: none;
      border: none;
      color: var(--text-secondary);
      font-family: var(--font-body);
      font-size: 0.85rem;
      font-weight: 500;
      padding: 0.45rem 1.1rem;
      border-radius: 6px;
      cursor: pointer;
      transition: all 0.15s ease;
      display: flex;
      align-items: center;
      gap: 0.4rem;
    }}

    .tab-btn:hover {{
      color: var(--text-primary);
      background: rgba(0, 0, 0, 0.04);
    }}

    .tab-btn.active {{
      background: #ffffff;
      color: var(--text-primary);
      font-weight: 600;
      border: 1px solid var(--border-color);
      box-shadow: 0 1px 2px rgba(0, 0, 0, 0.05);
    }}

    /* Main Container */
    main {{
      max-width: 1560px;
      margin: 0 auto;
      padding: 1.75rem 2rem;
    }}

    .tab-pane {{
      display: none;
      animation: fadeIn 0.2s ease-in-out;
    }}

    .tab-pane.active {{
      display: block;
    }}

    @keyframes fadeIn {{
      from {{ opacity: 0; transform: translateY(3px); }}
      to {{ opacity: 1; transform: translateY(0); }}
    }}

    /* KPI Grid - Flat Surfaces, Typography First */
    .kpi-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
      gap: 1rem;
      margin-bottom: 1.75rem;
    }}

    .kpi-card {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 1.15rem 1.25rem;
      box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    }}

    .kpi-card h3 {{
      font-size: 0.72rem;
      text-transform: uppercase;
      font-weight: 600;
      letter-spacing: 0.05em;
      color: var(--text-secondary);
      margin-bottom: 0.35rem;
    }}

    .kpi-card .value {{
      font-family: var(--font-body);
      font-size: 1.75rem;
      font-weight: 700;
      color: var(--text-primary);
      font-variant-numeric: tabular-nums;
      line-height: 1.2;
    }}

    .kpi-card .subtitle {{
      font-size: 0.72rem;
      color: var(--text-muted);
      margin-top: 0.3rem;
    }}

    /* Panels & Cards */
    .panel {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 1.5rem;
      margin-bottom: 1.75rem;
      box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    }}

    .panel-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 1.25rem;
      border-bottom: 1px solid var(--border-color);
      padding-bottom: 0.85rem;
    }}

    .panel-header h2 {{
      font-family: var(--font-body);
      font-size: 1.15rem;
      font-weight: 700;
      color: var(--text-primary);
      letter-spacing: -0.01em;
    }}

    /* Data Tables: Clean Rules and Right Alignment */
    .table-container {{
      overflow-x: auto;
      max-height: 560px;
      border: 1px solid var(--border-color);
      border-radius: 6px;
    }}

    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.825rem;
      text-align: left;
    }}

    th {{
      background: var(--bg-subtle);
      color: var(--text-secondary);
      font-weight: 600;
      font-size: 0.72rem;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      padding: 0.65rem 0.85rem;
      position: sticky;
      top: 0;
      z-index: 10;
      border-bottom: 1px solid var(--border-strong);
    }}

    td {{
      padding: 0.65rem 0.85rem;
      border-bottom: 1px solid var(--border-color);
      color: var(--text-primary);
      font-size: 0.825rem;
    }}

    tr:hover td {{
      background: #f8fafc;
    }}

    /* Search & Filter inputs */
    .toolbar {{
      display: flex;
      gap: 0.75rem;
      margin-bottom: 1.25rem;
      flex-wrap: wrap;
      align-items: center;
    }}

    .search-input {{
      background: #ffffff;
      border: 1px solid var(--border-strong);
      color: var(--text-primary);
      padding: 0.5rem 0.85rem;
      border-radius: 6px;
      font-size: 0.85rem;
      min-width: 280px;
    }}

    .search-input:focus {{
      outline: none;
      border-color: var(--accent-blue);
      box-shadow: 0 0 0 2px rgba(2, 132, 199, 0.15);
    }}

    .btn {{
      background: #ffffff;
      border: 1px solid var(--border-strong);
      color: var(--text-primary);
      padding: 0.45rem 0.9rem;
      border-radius: 6px;
      font-size: 0.825rem;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.15s ease;
      display: inline-flex;
      align-items: center;
      gap: 4px;
    }}

    .btn:hover {{
      background: #f8fafc;
      border-color: #94a3b8;
    }}

    .btn.active, .btn-primary {{
      background: var(--text-primary);
      color: #ffffff;
      border-color: var(--text-primary);
    }}

    .btn.active:hover, .btn-primary:hover {{
      background: #1e293b;
      color: #ffffff;
    }}

    .btn-sm {{
      padding: 0.3rem 0.65rem;
      font-size: 0.75rem;
    }}

    /* Oncoprint Styles conforming to cBioPortal Hierarchy & Taxonomy */
    .oncoprint-wrapper {{
      overflow-x: auto;
      padding: 1.25rem 0.75rem;
      background: #ffffff;
      border: 1px solid var(--border-color);
      border-radius: 8px;
    }}

    .oncoprint-table {{
      width: max-content !important;
      border-collapse: separate;
      border-spacing: 2px 2px;
      table-layout: auto;
      user-select: none;
    }}

    .oncoprint-table td {{
      padding: 0;
      vertical-align: middle;
    }}

    .oncoprint-table th.gene-header {{
      width: 110px;
      min-width: 110px;
      text-align: right;
      padding: 4px 8px 4px 0;
      font-size: 0.72rem;
      font-weight: 700;
      color: var(--text-secondary);
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }}

    .oncoprint-table th.sample-header {{
      width: 28px;
      min-width: 28px;
      max-width: 28px;
      height: 130px;
      vertical-align: bottom;
      padding: 4px 1px;
      text-align: center;
    }}

    .oncoprint-table th.patient-boundary,
    .oncoprint-table td.patient-boundary {{
      padding-right: 14px !important;
    }}

    .oncoprint-table th.patient-group-boundary {{
      padding-right: 14px !important;
    }}

    .oncoprint-table th.patient-group-header {{
      font-size: 0.72rem;
      font-weight: 700;
      color: var(--text-primary);
      background: var(--bg-subtle);
      border: 1px solid var(--border-strong);
      border-radius: 4px;
      padding: 3px 6px;
      text-align: center;
      letter-spacing: 0.5px;
      white-space: nowrap;
    }}

    .sample-type-badge {{
      display: inline-block;
      font-size: 0.55rem;
      font-weight: 700;
      padding: 1px 3px;
      border-radius: 3px;
      margin-top: 2px;
      white-space: nowrap;
      text-transform: uppercase;
      letter-spacing: 0.2px;
    }}

    .sample-type-badge.type-primary {{
      background: #e0f2fe;
      color: #0369a1;
      border: 1px solid #bae6fd;
    }}

    .sample-type-badge.type-plasma {{
      background: #f3e8ff;
      color: #6b21a8;
      border: 1px solid #e9d5ff;
    }}

    .oncoprint-table th.sample-header .sample-label {{
      writing-mode: vertical-rl;
      transform: rotate(180deg);
      font-size: 0.68rem;
      font-weight: 600;
      color: #334155;
      white-space: nowrap;
      letter-spacing: 0.3px;
      display: inline-block;
      max-height: 125px;
      overflow: hidden;
      text-overflow: ellipsis;
    }}

    .oncoprint-table th.freq-header {{
      width: 48px;
      min-width: 48px;
      text-align: left;
      font-size: 0.72rem;
      font-weight: 700;
      color: var(--text-secondary);
      padding-left: 8px;
    }}

    .oncoprint-table td.gene-name {{
      font-weight: 700;
      font-size: 0.78rem;
      color: var(--text-primary);
      padding: 0 8px 0 0;
      white-space: nowrap;
      text-align: right;
    }}

    .oncoprint-table td.freq-val {{
      text-align: left;
      font-size: 0.72rem;
      font-weight: 600;
      color: var(--text-secondary);
      padding-left: 8px;
      font-variant-numeric: tabular-nums;
    }}

    /* Track Header Rows (CUPPA, Purity, TMB, MSI) */
    .onco-track-row td {{
      padding: 0;
    }}

    .onco-track-label {{
      font-size: 0.7rem;
      font-weight: 600;
      color: var(--text-secondary);
      padding: 0 8px 0 0;
      white-space: nowrap;
      text-align: right;
    }}

    .onco-track-cell {{
      width: 28px;
      min-width: 28px;
      max-width: 28px;
      height: 18px;
      position: relative;
      background: var(--bg-subtle);
      border: 1px solid var(--border-color);
      border-radius: 2px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-size: 0.52rem;
      font-weight: 600;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      padding: 0 1px;
      box-sizing: border-box;
      color: var(--text-primary);
    }}

    /* ---------------------------------------------------------
       cBioPortal Visual Layering Hierarchy (Z-Order)
       Base Canvas: Crisp light-grey wild-type cell
       Cell: 28px width x 32px height
       --------------------------------------------------------- */
    .onco-cell-container {{
      width: 28px;
      min-width: 28px;
      max-width: 28px;
      height: 32px;
      position: relative;
      background: #e2e8f0;
      border-radius: 1px;
      cursor: pointer;
      box-sizing: border-box;
      transition: transform 0.1s ease;
    }}

    .onco-cell-container:hover {{
      outline: 2px solid var(--accent-blue);
      outline-offset: 1px;
      z-index: 50 !important;
      transform: scale(1.15);
    }}

    /* Floating Interactive Oncoprint Hover Tooltip */
    .onco-floating-tooltip {{
      position: fixed;
      z-index: 99999;
      pointer-events: none;
      background: #ffffff;
      border: 1px solid var(--border-strong);
      border-radius: 8px;
      padding: 10px 14px;
      box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.15), 0 8px 10px -6px rgba(0, 0, 0, 0.1);
      color: var(--text-primary);
      font-size: 0.82rem;
      line-height: 1.4;
      max-width: 340px;
      transition: opacity 0.08s ease, transform 0.08s ease;
      opacity: 0;
      transform: translateY(4px);
    }}

    .onco-floating-tooltip.visible {{
      opacity: 1;
      transform: translateY(0);
    }}

    .onco-tooltip-gene {{
      font-size: 1.05rem;
      font-weight: 700;
      color: var(--accent-blue);
      margin-bottom: 2px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 8px;
    }}

    .onco-tooltip-sample {{
      font-size: 0.76rem;
      color: var(--text-secondary);
      margin-bottom: 6px;
      padding-bottom: 6px;
      border-bottom: 1px solid var(--border-color);
    }}

    .onco-tooltip-mut-item {{
      margin-top: 4px;
      padding: 4px 8px;
      background: #f8fafc;
      border-radius: 4px;
      border: 1px solid var(--border-color);
      border-left: 3px solid #16a34a;
    }}

    .onco-tooltip-mut-item.truncating {{
      border-left-color: #0f172a;
    }}

    .onco-tooltip-mut-item.germline {{
      border-left-color: #7c3aed;
    }}

    .onco-tooltip-aa {{
      font-weight: 700;
      color: #0284c7;
      font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
      font-size: 0.88rem;
    }}

    .onco-tooltip-tag {{
      font-size: 0.7rem;
      padding: 1px 5px;
      border-radius: 3px;
      background: #f1f5f9;
      color: #334155;
      display: inline-block;
      margin-left: 4px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.2px;
      border: 1px solid #e2e8f0;
    }}

    .onco-layer-not-profiled {{
      position: absolute;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background: repeating-linear-gradient(
        45deg,
        rgba(0, 0, 0, 0.04),
        rgba(0, 0, 0, 0.04) 2px,
        transparent 2px,
        transparent 6px
      );
      z-index: 0;
    }}

    /* Layer 1: Copy Number Alterations (Full-height filled rectangle) */
    .onco-layer-cna {{
      position: absolute;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      z-index: 1;
      border-radius: 1px;
    }}

    .cna-amp {{ background-color: #DC2626 !important; }}       /* High-level Gain: Red */
    .cna-homdel {{ background-color: #2563EB !important; }}    /* Homozygous Loss: Blue */
    .cna-gain {{ background-color: #FCA5A5 !important; }}      /* Low-level Gain: Light Red */
    .cna-hetloss {{ background-color: #93C5FD !important; }}   /* Heterozygous Loss: Light Blue */

    /* Layer 1 Modifier: Expression Outlines */
    .expr-up {{
      box-shadow: inset 0 0 0 2px #DC2626;
      z-index: 1;
    }}

    .expr-down {{
      box-shadow: inset 0 0 0 2px #2563EB;
      z-index: 1;
    }}

    /* Layer 2: Mutation Layer (Centered sub-height rectangle ~37.5% of cell: 12px / 32px) */
    .onco-layer-mut {{
      position: absolute;
      top: 50%;
      left: 0;
      width: 100%;
      height: 12px;
      transform: translateY(-50%);
      z-index: 2;
      border-radius: 1px;
    }}

    .mut-missense {{ background-color: #16A34A !important; }}    /* Solid Green */
    .mut-truncating {{ background-color: #0F172A !important; }}  /* Solid Black / Dark Charcoal */
    .mut-inframe {{ background-color: #B45309 !important; }}     /* Ochre / Brown */
    .mut-promoter {{ background-color: #0D9488 !important; }}    /* Teal */
    .mut-splice {{ background-color: #EA580C !important; }}      /* Orange */

    /* Layer 3: Structural Variant / Fusion Layer */
    .onco-layer-fusion {{
      position: absolute;
      top: 0;
      left: 50%;
      width: 6px;
      height: 100%;
      transform: translateX(-50%);
      background-color: #9333EA;                               /* Magenta / Purple vertical bar */
      z-index: 3;
      border-radius: 1px;
    }}

    /* Layer 4: Germline / Modifier Layer (Centered horizontal white overlay stripe) */
    .onco-layer-germline {{
      position: absolute;
      top: 50%;
      left: 0;
      width: 100%;
      height: 2px;
      transform: translateY(-50%);
      background-color: #FFFFFF !important;
      z-index: 4;
      pointer-events: none;
    }}

    /* cBioPortal Oncoprint Legend Box */
    .oncoprint-legend-box {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
      gap: 1.25rem;
      margin-top: 1.25rem;
      padding: 1.15rem;
      background: #ffffff;
      border: 1px solid var(--border-color);
      border-radius: 8px;
      font-size: 0.8rem;
    }}

    .legend-group-title {{
      font-weight: 700;
      color: var(--text-secondary);
      text-transform: uppercase;
      font-size: 0.72rem;
      letter-spacing: 0.5px;
      margin-bottom: 0.6rem;
    }}

    .legend-items-list {{
      display: flex;
      flex-direction: column;
      gap: 0.45rem;
    }}

    .onco-legend-item {{
      display: flex;
      align-items: center;
      gap: 0.6rem;
      color: var(--text-primary);
    }}

    .legend-glyph-box {{
      width: 16px;
      height: 24px;
      position: relative;
      background: #e2e8f0;
      border: 1px solid var(--border-strong);
      border-radius: 1px;
      flex-shrink: 0;
      display: flex;
      align-items: center;
      justify-content: center;
      box-sizing: border-box;
    }}

    .legend-glyph-mut {{
      position: absolute;
      top: 50%;
      left: 0;
      width: 100%;
      height: 9px;
      transform: translateY(-50%);
      border-radius: 1px;
    }}

    .legend-glyph-fusion {{
      position: absolute;
      top: 0;
      left: 50%;
      width: 4px;
      height: 100%;
      transform: translateX(-50%);
      background-color: #9333EA;
      border-radius: 1px;
    }}

    .legend-glyph-germline {{
      position: absolute;
      top: 50%;
      left: 0;
      width: 100%;
      height: 2px;
      transform: translateY(-50%);
      background-color: #FFFFFF !important;
      z-index: 4;
    }}

    /* Chart Grid */
    .charts-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(460px, 1fr));
      gap: 1.5rem;
      margin-bottom: 2rem;
    }}

    .chart-box {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 1.5rem;
      min-height: 380px;
      box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    }}

    .chart-box h3 {{
      font-family: var(--font-body);
      font-size: 0.95rem;
      font-weight: 700;
      margin-bottom: 1rem;
      color: var(--text-primary);
    }}
  </style>
</head>
<body>

  <header>
    <div class="logo-group">
      <h1>
        <span>{cohort['cohort_name']}</span>
        <span class="badge badge-primary">WGTS Multi-Omics</span>
      </h1>
      <p>nf-core/oncoanalyser Cohort-Level Summary & MRD Trajectory Dashboard</p>
    </div>

    <nav class="tabs-nav">
      <button class="tab-btn active" onclick="switchTab('tab-overview')">Overview</button>
      <button class="tab-btn" onclick="switchTab('tab-qc')">QC Plots</button>
      <button class="tab-btn" onclick="switchTab('tab-genomics')">Genomic Summary</button>
      <button class="tab-btn" onclick="switchTab('tab-mrd')">MRD / WISP</button>
      <button class="tab-btn" onclick="switchTab('tab-fingerprint')">🧬 Sample Fingerprinting</button>
    </nav>
  </header>

  <main>

    <!-- =================================================================== -->
    <!-- TAB 1: COHORT OVERVIEW (LANDING PAGE) -->
    <!-- =================================================================== -->
    <div id="tab-overview" class="tab-pane active">
      <div class="kpi-grid">
        <div class="kpi-card">
          <h3>Total Patients</h3>
          <div class="value">{cohort['total_analyzed_patients']}</div>
          <div class="subtitle">Cohort Size</div>
        </div>
        <div class="kpi-card">
          <h3>Primary Tumors</h3>
          <div class="value">{len(cohort['all_tumors'])}</div>
          <div class="subtitle">WGS DNA Genomes</div>
        </div>
        <div class="kpi-card">
          <h3>Longitudinal Plasmas</h3>
          <div class="value">{len(cohort['all_plasmas'])}</div>
          <div class="subtitle">ctDNA Samples (WISP)</div>
        </div>
        <div class="kpi-card">
          <h3>Median Purity</h3>
          <div class="value" id="kpi-purity">--%</div>
          <div class="subtitle">PURPLE fit estimate</div>
        </div>
        <div class="kpi-card">
          <h3>Median Ploidy</h3>
          <div class="value" id="kpi-ploidy">--</div>
          <div class="subtitle">Genomic copy number</div>
        </div>
        <div class="kpi-card">
          <h3>Median Coverage</h3>
          <div class="value" id="kpi-coverage">--x</div>
          <div class="subtitle">Tumor Mean Depth</div>
        </div>
        <div class="kpi-card">
          <h3>Median TMB</h3>
          <div class="value" id="kpi-tmb">--</div>
          <div class="subtitle">Mutations / Mb</div>
        </div>
      </div>

      <div class="panel">
        <div class="panel-header">
          <h2>Samples & Patient Summary</h2>
          <button class="btn" onclick="exportTableToCSV('overview-table', 'cohort_wgts_overview.csv')">Export CSV</button>
        </div>

        <div class="toolbar">
          <input type="text" id="overview-search" class="search-input" placeholder="Search by patient ID, sample ID, status..." onkeyup="filterOverviewTable()">
        </div>

        <div class="table-container">
          <table id="overview-table">
            <thead>
              <tr>
                <th>Patient ID</th>
                <th>Primary Tumor</th>
                <th>CUPPA Cancer Type</th>
                <th>Purity (%)</th>
                <th>Ploidy</th>
                <th>Coverage (Tumor / Normal)</th>
                <th>TMB (mut/Mb)</th>
                <th>MSI Status</th>
                <th>WGD</th>
                <th>ichorCNA TF</th>
                <th>Fingerprint Identity</th>
                <th>Plasmas</th>
                <th>Latest MRD</th>
                <th>MRD Trajectory</th>
                <th>Key Altered Genes</th>
              </tr>
            </thead>
            <tbody id="overview-tbody">
              <!-- Rendered via JS -->
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- TAB 2: QC PLOTS -->
    <!-- =================================================================== -->
    <div id="tab-qc" class="tab-pane">
      <div class="charts-grid">
        <div class="chart-box">
          <h3>Primary Tumor Purity vs. Ploidy</h3>
          <div style="height: 320px; position: relative;">
            <canvas id="chart-purity-ploidy"></canvas>
          </div>
        </div>

        <div class="chart-box">
          <h3>Tumor Mutational Burden (TMB) by MSI Status</h3>
          <div style="height: 320px; position: relative;">
            <canvas id="chart-tmb"></canvas>
          </div>
        </div>

        <div class="chart-box">
          <h3>WGS Sequencing Coverage (Tumor, Plasma ctDNA & Normal)</h3>
          <div style="height: 320px; position: relative;">
            <canvas id="chart-coverage"></canvas>
          </div>
        </div>

        <div class="chart-box">
          <h3>Variant Consequence Breakdown (Somatic Drivers & Focal CNAs)</h3>
          <div style="height: 320px; position: relative;">
            <canvas id="chart-consequences"></canvas>
          </div>
        </div>

        <div class="chart-box">
          <h3>sWGS ichorCNA vs PURPLE Purity Concordance</h3>
          <div style="height: 320px; position: relative;">
            <canvas id="chart-ichor-concordance"></canvas>
          </div>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- TAB 3: GENOMIC SUMMARY (ONCOPRINT) -->
    <!-- =================================================================== -->
    <div id="tab-genomics" class="tab-pane">
      <div class="panel">
        <div class="panel-header">
          <h2>Genomic Landscape & Oncoprint</h2>
          <div style="display: flex; gap: 0.5rem;">
            <button class="btn" onclick="renderOncoprint(10)">Top 10</button>
            <button class="btn" onclick="renderOncoprint(20)">Top 20</button>
            <button class="btn" onclick="renderOncoprint(999)">All Prioritized</button>
          </div>
        </div>

        <!-- CUPPA Cancer of Unknown Primary Cohort Distribution Bar -->
        <div id="cuppa-distribution-panel" style="margin-bottom: 1.25rem; padding: 0.85rem 1.15rem; background: var(--bg-subtle); border: 1px solid var(--border-color); border-radius: 6px;">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem; flex-wrap: wrap; gap: 0.5rem;">
            <div style="display: flex; align-items: center; gap: 0.5rem;">
              <span style="font-size: 0.75rem; font-weight: 700; color: var(--text-secondary); text-transform: uppercase; letter-spacing: 0.05em;">CUPPA Predicted Cancer Types (HMF Classifier):</span>
              <span id="cuppa-cohort-total" style="font-size: 0.72rem; color: var(--text-muted);"></span>
            </div>
            <span style="font-size: 0.7rem; color: var(--text-muted); font-style: italic;">Primary tumor whole-genome classification</span>
          </div>
          <div id="cuppa-distribution-pills" style="display: flex; flex-wrap: wrap; gap: 0.5rem; align-items: center;">
            <!-- Dynamically populated via JS -->
          </div>
        </div>

        <div class="toolbar" style="flex-direction: column; align-items: stretch; gap: 0.75rem;">
          <div style="display: flex; gap: 0.5rem; align-items: center;">
            <input type="text" id="gene-add-input" list="available-genes-list" class="search-input" style="flex: 1; max-width: 420px;" placeholder="Type or select gene symbol to add (e.g. PIK3CA, PTEN, ARAF, CTNNB1)..." onkeydown="if(event.key==='Enter') addCustomGene()">
            <datalist id="available-genes-list"></datalist>
            <button class="btn btn-primary" onclick="addCustomGene()">Add Gene</button>
          </div>
          <div id="orange-driver-pills-bar" style="display: flex; flex-wrap: wrap; gap: 0.4rem; align-items: center; padding-top: 0.25rem;">
            <span style="font-size: 0.75rem; color: var(--text-secondary); text-transform: uppercase; font-weight: 600; margin-right: 0.25rem;">Orange Report Drivers:</span>
            <div id="orange-driver-pills" style="display: inline-flex; flex-wrap: wrap; gap: 0.35rem;"></div>
          </div>
        </div>

        <div class="oncoprint-wrapper">
          <div id="oncoprint-container">
            <!-- Dynamically populated Oncoprint Grid -->
          </div>
        </div>

        <div class="oncoprint-legend-box">
          <div class="legend-group">
            <div class="legend-group-title">Mutations (Small Variants & Indels)</div>
            <div class="legend-items-list">
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-mut mut-missense"></div></div>
                <span>Missense Mutation</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-mut mut-truncating"></div></div>
                <span>Truncating (Nonsense, Frameshift)</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-mut mut-inframe"></div></div>
                <span>In-Frame InDel / Insertion / Deletion</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-mut mut-splice"></div></div>
                <span>Splice Mutation</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-mut mut-promoter"></div></div>
                <span>Promoter / Non-coding Mutation</span>
              </div>
            </div>
          </div>

          <div class="legend-group">
            <div class="legend-group-title">Copy Number Alterations (CNA)</div>
            <div class="legend-items-list">
              <div class="onco-legend-item">
                <div class="legend-glyph-box cna-amp"></div>
                <span>Amplification (High-level Gain, '2')</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box cna-homdel"></div>
                <span>Deep Deletion (Homozygous Loss, '-2')</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box cna-gain"></div>
                <span>Gain (Low-level Gain, '1')</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box cna-hetloss"></div>
                <span>Shallow Deletion (Heterozygous Loss, '-1')</span>
              </div>
            </div>
          </div>

          <div class="legend-group">
            <div class="legend-group-title">Structural & Modifiers</div>
            <div class="legend-items-list">
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-fusion"></div></div>
                <span>SV / Gene Fusion</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box"><div class="legend-glyph-mut mut-missense"></div><div class="legend-glyph-germline"></div></div>
                <span>Germline Mutation (White Stripe)</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box expr-up"></div>
                <span>mRNA / Protein Upregulation</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box expr-down"></div>
                <span>mRNA / Protein Downregulation</span>
              </div>
              <div class="onco-legend-item">
                <div class="legend-glyph-box"></div>
                <span>Wild Type / Unaltered</span>
              </div>
            </div>
          </div>

          <div class="legend-group">
            <div class="legend-group-title">CUPPA Predicted Cancer Types</div>
            <div class="legend-items-list">
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #ffe4e6; border: 1px solid #f43f5e;"></div>
                <span>Breast</span>
              </div>
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #ede9fe; border: 1px solid #8b5cf6;"></div>
                <span>Endometrium / Uterus</span>
              </div>
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #fae8ff; border: 1px solid #d946ef;"></div>
                <span>Ovary / Fallopian Tube</span>
              </div>
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #fef3c7; border: 1px solid #f59e0b;"></div>
                <span>Bone / Soft Tissue</span>
              </div>
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #d1fae5; border: 1px solid #10b981;"></div>
                <span>Head and Neck</span>
              </div>
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #e0f2fe; border: 1px solid #0ea5e9;"></div>
                <span>Colorectum</span>
              </div>
              <div class="onco-legend-item">
                <div style="width: 14px; height: 14px; border-radius: 3px; background: #f1f5f9; border: 1px solid #94a3b8;"></div>
                <span>Other / Low Confidence</span>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div class="panel">
        <div class="panel-header">
          <h2>Detected Somatic Drivers & High-Confidence Alterations</h2>
          <button class="btn" onclick="exportTableToCSV('drivers-table', 'cohort_drivers.csv')">Export CSV</button>
        </div>
        <div class="table-container">
          <table id="drivers-table">
            <thead>
              <tr>
                <th>Patient ID</th>
                <th>Tumor ID</th>
                <th>Gene</th>
                <th>Classification</th>
                <th>Protein Change</th>
                <th>Coding Change</th>
                <th>VAF</th>
                <th>Driver Status</th>
              </tr>
            </thead>
            <tbody id="drivers-tbody">
              <!-- Rendered via JS -->
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- TAB 4: MRD / WISP (LOG-SCALED WITH LIMIT OF DETECTION SHADING) -->
    <!-- =================================================================== -->
    <div id="tab-mrd" class="tab-pane">
      <div id="mrd-empty-notice" style="display: none; padding: 3rem; text-align: center; color: var(--text-secondary);">
        <h3>No plasma ctDNA / WISP samples detected in this cohort</h3>
        <p>This tab activates automatically when longitudinal plasma genomes are processed.</p>
      </div>

      <div id="mrd-content">
        <div class="panel">
          <div class="panel-header" style="border-bottom: 1px solid rgba(255,255,255,0.08); padding-bottom: 0.85rem; margin-bottom: 1rem;">
            <div>
              <h2 style="margin:0 0 4px 0;">Longitudinal ctDNA Trajectory (Clinical Standard)</h2>
              <p style="margin:0; font-size:0.825rem; color:var(--text-secondary);">
                Log-scale tumor fraction tracking with baseline WGS timestamp, limit of detection floor, and confidence intervals.
              </p>
            </div>
            <div style="display: flex; gap: 0.75rem; align-items: center; flex-wrap: wrap;">
              <div style="display: flex; gap: 0.35rem;">
                <button class="btn btn-sm active" id="btn-mrd-view-clinical" onclick="setMrdViewMode('CLINICAL')">📈 Clinical Trajectory (Standard)</button>
                <button class="btn btn-sm" id="btn-mrd-view-tracks" onclick="setMrdViewMode('TRACKS')">📊 Multi-Modal Breakdown (SNV/CNV/LOH)</button>
              </div>
              <div style="display: flex; gap: 0.4rem; align-items: center;">
                <label style="font-size: 0.85rem; color: var(--text-secondary);">Patient:</label>
                <select id="mrd-patient-select" class="search-input" style="min-width: 180px;" onchange="updateMrdTrajectoryPlot()">
                  <option value="ALL">All Cohort Patients (Overview)</option>
                </select>
                <a id="btn-download-mrd-png" href="#" download class="btn btn-sm" style="display:none; color:#34d399; border-color:#059669; font-weight:600; white-space:nowrap; gap:4px; align-items:center; text-decoration:none;">⬇ PNG</a>
                <a id="btn-download-mrd-pdf" href="#" download class="btn btn-sm" style="display:none; color:#60a5fa; border-color:#2563eb; font-weight:600; white-space:nowrap; gap:4px; align-items:center; text-decoration:none;">⬇ PDF</a>
                <button id="btn-view-patient-trajectory" class="btn btn-sm" onclick="viewSelectedPatientTrajectory()" style="display:none; color:#34d399; border-color:#059669; font-weight:600; white-space:nowrap; gap:4px; align-items:center;">🔍 Full Screen</button>
              </div>
            </div>
          </div>

          <!-- Mode 1: Clinical Trajectory (Standard Attached Image Style) -->
          <div id="mrd-clinical-view">
            <!-- Single Patient View -->
            <div id="mrd-clinical-single-patient" style="display:none; text-align:center;">
              <div id="mrd-clinical-img-container" style="max-width:1100px; margin:0 auto; background:#ffffff; border-radius:8px; padding:12px; border:1px solid #30363d; box-shadow:0 4px 16px rgba(0,0,0,0.35);">
                <img id="mrd-clinical-img" src="" alt="Clinical ctDNA Trajectory" style="max-width:100%; height:auto; display:block; margin:0 auto; border-radius:4px;" />
              </div>
            </div>

            <!-- Cohort Overview Grid when 'ALL' is selected -->
            <div id="mrd-clinical-cohort-grid" style="display:grid; grid-template-columns:repeat(auto-fit, minmax(500px, 1fr)); gap:1.25rem;">
            </div>
          </div>

          <!-- Mode 2: Multi-Modal Tracks View (SNV / CNV / LOH) -->
          <div id="mrd-tracks-view" style="display:none;">
            <!-- Single Patient / Cohort Multi-Modal Publication Figure View -->
            <div id="mrd-tracks-single-patient" style="text-align:center; margin-bottom:1.5rem;">
              <div id="mrd-tracks-img-container" style="max-width:1100px; margin:0 auto; background:#ffffff; border-radius:8px; padding:12px; border:1px solid #30363d; box-shadow:0 4px 16px rgba(0,0,0,0.35);">
                <img id="mrd-tracks-img" src="" alt="Multi-Modal ctDNA Trajectory" style="max-width:100%; height:auto; display:block; margin:0 auto; border-radius:4px;" />
              </div>
            </div>

            <div style="background: rgba(255, 255, 255, 0.02); border-radius: 8px; padding: 0.75rem 1rem; margin-bottom: 1rem; font-size: 0.825rem; color: var(--text-secondary); display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.5rem;">
              <div>
                <strong>Interactive Modality Breakdown:</strong> Longitudinal trajectory across <strong>SNV</strong> (<span style="color:#ef4444;">● Red</span>), <strong>Copy Number</strong> (<span style="color:#a855f7;">● Purple</span>), <strong>LOH</strong> (<span style="color:#f59e0b;">● Amber</span>), and <strong>1X ichorCNA</strong> (<span style="color:#d97706;">◆ Gold Diamond</span>) formatted identically to the trajectory plots with continuous LOD floor shading and CI ribbons.
              </div>
              <div style="display:flex; gap:0.4rem; align-items:center;">
                <div style="display: flex; gap: 0.35rem;">
                  <button class="btn btn-sm active" id="btn-mrd-mod-all" onclick="setMrdModality('ALL')">All Tracks</button>
                  <button class="btn btn-sm" id="btn-mrd-mod-snv" onclick="setMrdModality('SNV')">SNV</button>
                  <button class="btn btn-sm" id="btn-mrd-mod-cnv" onclick="setMrdModality('CNV')">CNV</button>
                  <button class="btn btn-sm" id="btn-mrd-mod-loh" onclick="setMrdModality('LOH')">LOH</button>
                  <button class="btn btn-sm" id="btn-mrd-mod-ichor" onclick="setMrdModality('ICHOR')">ichorCNA</button>
                </div>
                <label style="font-size:0.8rem; color:var(--text-primary); cursor:pointer; display:flex; align-items:center; gap:0.35rem; margin-left:0.5rem;">
                  <input type="checkbox" id="toggle-sample-lods" checked onchange="updateMrdTrajectoryPlot()"> Show LOD Shading
                </label>
              </div>
            </div>

            <!-- Stacked Multi-Modal Tracks (Interactive Canvases) -->
            <div id="track-card-snv" style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:14px 18px 10px 18px; margin-bottom:16px; box-shadow:0 2px 8px rgba(0,0,0,0.06);">
              <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px; flex-wrap: wrap; gap: 8px;">
                <div style="font-weight: 700; font-size: 0.95rem; color: #0f172a; display: flex; align-items: center; gap: 8px;">
                  <span style="color:#ef4444;">● Track 1: Longitudinal ctDNA Trajectory (SNV Modality)</span>
                  <span class="badge" style="background: rgba(239, 68, 68, 0.1); color: #ef4444; border: 1px solid rgba(239, 68, 68, 0.3);">High Sensitivity</span>
                </div>
                <div style="display:flex; gap:12px; font-size:0.75rem; color:#475569; align-items:center; flex-wrap:wrap;">
                  <span><span style="display:inline-block; width:12px; height:0; border-top:2px dotted #000; vertical-align:middle; margin-right:3px;"></span>WGS Baseline</span>
                  <span><span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:#ef4444; vertical-align:middle; margin-right:3px;"></span>Detected (SNV)</span>
                  <span><span style="display:inline-block; width:12px; height:8px; background:#fee2e2; vertical-align:middle; margin-right:3px;"></span>95% CI Ribbon</span>
                  <span><span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:#6b7280; vertical-align:middle; margin-right:3px;"></span>Below LOD</span>
                  <span><span style="display:inline-block; width:14px; height:0; border-top:2px dashed #94a3b8; vertical-align:middle; margin-right:3px;"></span>LOD Line</span>
                </div>
              </div>
              <div style="height: 220px; position: relative;">
                <canvas id="chart-mrd-snv"></canvas>
              </div>
            </div>

            <div id="track-card-cnv" style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:14px 18px 10px 18px; margin-bottom:16px; box-shadow:0 2px 8px rgba(0,0,0,0.06);">
              <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px; flex-wrap: wrap; gap: 8px;">
                <div style="font-weight: 700; font-size: 0.95rem; color: #0f172a; display: flex; align-items: center; gap: 8px;">
                  <span style="color:#a855f7;">● Track 2: Longitudinal ctDNA Trajectory (CNV Modality)</span>
                  <span class="badge" style="background: rgba(168, 85, 247, 0.1); color: #a855f7; border: 1px solid rgba(168, 85, 247, 0.3);">Broad Genome</span>
                </div>
                <div style="display:flex; gap:12px; font-size:0.75rem; color:#475569; align-items:center; flex-wrap:wrap;">
                  <span><span style="display:inline-block; width:12px; height:0; border-top:2px dotted #000; vertical-align:middle; margin-right:3px;"></span>WGS Baseline</span>
                  <span><span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:#a855f7; vertical-align:middle; margin-right:3px;"></span>Detected (CNV)</span>
                  <span><span style="display:inline-block; width:12px; height:8px; background:#f3e8ff; vertical-align:middle; margin-right:3px;"></span>95% CI Ribbon</span>
                  <span><span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:#6b7280; vertical-align:middle; margin-right:3px;"></span>Below LOD</span>
                  <span><span style="display:inline-block; width:14px; height:0; border-top:2px dashed #94a3b8; vertical-align:middle; margin-right:3px;"></span>LOD Line</span>
                </div>
              </div>
              <div style="height: 220px; position: relative;">
                <canvas id="chart-mrd-cnv"></canvas>
              </div>
            </div>

            <div id="track-card-loh" style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:14px 18px 10px 18px; margin-bottom:16px; box-shadow:0 2px 8px rgba(0,0,0,0.06);">
              <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px; flex-wrap: wrap; gap: 8px;">
                <div style="font-weight: 700; font-size: 0.95rem; color: #0f172a; display: flex; align-items: center; gap: 8px;">
                  <span style="color:#f59e0b;">● Track 3: Longitudinal ctDNA Trajectory (LOH Modality)</span>
                  <span class="badge" style="background: rgba(245, 158, 11, 0.1); color: #f59e0b; border: 1px solid rgba(245, 158, 11, 0.3);">Intermediate Sensitivity</span>
                </div>
                <div style="display:flex; gap:12px; font-size:0.75rem; color:#475569; align-items:center; flex-wrap:wrap;">
                  <span><span style="display:inline-block; width:12px; height:0; border-top:2px dotted #000; vertical-align:middle; margin-right:3px;"></span>WGS Baseline</span>
                  <span><span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:#f59e0b; vertical-align:middle; margin-right:3px;"></span>Detected (LOH)</span>
                  <span><span style="display:inline-block; width:12px; height:8px; background:#fef3c7; vertical-align:middle; margin-right:3px;"></span>95% CI Ribbon</span>
                  <span><span style="display:inline-block; width:8px; height:8px; border-radius:50%; background:#6b7280; vertical-align:middle; margin-right:3px;"></span>Below LOD</span>
                  <span><span style="display:inline-block; width:14px; height:0; border-top:2px dashed #94a3b8; vertical-align:middle; margin-right:3px;"></span>LOD Line</span>
                </div>
              </div>
              <div style="height: 240px; position: relative;">
                <canvas id="chart-mrd-loh"></canvas>
              </div>
            </div>

            <div id="track-card-ichor" style="background:#ffffff; border:1px solid #cbd5e1; border-radius:8px; padding:14px 18px 10px 18px; margin-bottom:16px; box-shadow:0 2px 8px rgba(0,0,0,0.06);">
              <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px; flex-wrap: wrap; gap: 8px;">
                <div style="font-weight: 700; font-size: 0.95rem; color: #0f172a; display: flex; align-items: center; gap: 8px;">
                  <span style="color:#d97706;">◆ Track 4: 1X ichorCNA (Short Fragments 90–150 bp)</span>
                  <span class="badge" style="background: rgba(217, 119, 6, 0.1); color: #d97706; border: 1px solid rgba(217, 119, 6, 0.3);">Downsampled 1.0X</span>
                </div>
                <div style="display:flex; gap:12px; font-size:0.75rem; color:#475569; align-items:center; flex-wrap:wrap;">
                  <span><span style="display:inline-block; width:12px; height:0; border-top:2px dotted #000; vertical-align:middle; margin-right:3px;"></span>WGS Baseline</span>
                  <span><span style="display:inline-block; width:8px; height:8px; background:#d97706; transform:rotate(45deg); vertical-align:middle; margin-right:3px;"></span>Detected (TF &ge; 3%)</span>
                  <span><span style="display:inline-block; width:8px; height:8px; background:#6b7280; transform:rotate(45deg); vertical-align:middle; margin-right:3px;"></span>Below 3% LoD</span>
                  <span><span style="display:inline-block; width:14px; height:0; border-top:2px dashed #94a3b8; vertical-align:middle; margin-right:3px;"></span>3% LoD Floor</span>
                </div>
              </div>
              <div style="height: 240px; position: relative;">
                <canvas id="chart-mrd-ichor"></canvas>
              </div>
            </div>
          </div>
        </div>

        <div class="panel">
          <div class="panel-header">
            <h2>Plasma MRD Metrics & Multi-Modal LODs</h2>
            <button class="btn" onclick="exportTableToCSV('mrd-table', 'cohort_mrd_summary.csv')">Export CSV</button>
          </div>
          <div class="table-container">
            <table id="mrd-table">
              <thead>
                <tr>
                  <th>Patient ID</th>
                  <th>Sample ID</th>
                  <th>Timepoint</th>
                  <th>Overall MRD</th>
                  <th>SNV Purity (%)</th>
                  <th>SNV LOD (%)</th>
                  <th>CNV Purity (%)</th>
                  <th>CNV LOD (%)</th>
                  <th>LOH Purity (%)</th>
                  <th>LOH LOD (%)</th>
                  <th>WISP Status</th>
                  <th>ichorCNA TF (%)</th>
                  <th>Plasma Coverage (Depth)</th>
                  <th>Decomposition & Trajectory</th>
                </tr>
              </thead>
              <tbody id="mrd-tbody">
                <!-- Rendered via JS -->
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- TAB 5: SAMPLE FINGERPRINTING & IDENTITY CONFIRMATION -->
    <!-- =================================================================== -->
    <div id="tab-fingerprint" class="tab-pane">
      <div style="background:var(--bg-secondary); border:1px solid var(--border-color); border-radius:8px; padding:1.25rem 1.5rem; margin-bottom:1.5rem; box-shadow:0 1px 2px rgba(0,0,0,0.03);">
        <div style="display:flex; justify-content:space-between; align-items:flex-start; flex-wrap:wrap; gap:1rem;">
          <div>
            <h2 style="font-family:var(--font-head); font-size:1.3rem; font-weight:700; color:var(--text-primary); margin-bottom:0.3rem; display:flex; align-items:center; gap:0.5rem;">
              <span>🧬 Sample Identity & Germline Fingerprinting</span>
              <span class="badge" style="background:#ecfdf5; color:#065f46; border:1px solid #a7f3d0; font-size:0.75rem;">AMBER SNP BAF Concordance</span>
            </h2>
            <p style="font-size:0.875rem; color:var(--text-secondary); max-width:960px; line-height:1.5;">
              Pairwise germline heterozygous SNP concordance (|A &cap; B| / min(|A|, |B|)) computed across Normal, Primary Tumor, and serial Plasma ctDNA samples. Intra-patient concordance &ge; 85% (typically &gt; 98%) confirms identical biological donor origin across sequencing platforms (Illumina and Ultima) and excludes sample mislabeling or contamination.
            </p>
          </div>
          <div style="display:flex; gap:0.5rem;">
            <button class="btn" onclick="exportTableToCSV('fingerprint-table', 'cohort_fingerprint_summary.csv')">Export Concordance CSV</button>
            <button class="btn" onclick="exportFingerprintMatrixCSV()">Export Matrix CSV</button>
          </div>
        </div>
      </div>

      <div class="kpi-grid">
        <div class="kpi-card">
          <h3>Donor Concordance</h3>
          <div class="value" id="fp-kpi-rate" style="color:#059669;">--%</div>
          <div class="subtitle" id="fp-kpi-rate-sub">Confirmed Same Donor</div>
        </div>
        <div class="kpi-card">
          <h3>Fingerprinted Samples</h3>
          <div class="value" id="fp-kpi-samples">--</div>
          <div class="subtitle">Normal, Tumor & Plasma</div>
        </div>
        <div class="kpi-card">
          <h3>Mean Intra-Donor Concordance</h3>
          <div class="value" id="fp-kpi-intra" style="color:#0284c7;">--</div>
          <div class="subtitle">Expected &gt; 90% (Identical)</div>
        </div>
        <div class="kpi-card">
          <h3>Mean Inter-Donor Concordance</h3>
          <div class="value" id="fp-kpi-inter" style="color:#475569;">--</div>
          <div class="subtitle">Expected &lt; 35% (Distinct)</div>
        </div>
      </div>

      <!-- Panel 1: Intra-Patient Sample Concordance Table -->
      <div class="panel">
        <div class="panel-header">
          <h2>Intra-Patient Sample Concordance</h2>
          <span style="font-size:0.8rem; color:var(--text-muted);">Threshold: Intra-donor Concordance &ge; 85.0% confirms match</span>
        </div>
        <div class="toolbar" style="margin-bottom:1rem; display:flex; gap:1rem; align-items:center;">
          <input type="text" id="fp-search" class="search-input" placeholder="Search by patient ID, sample ID, status..." onkeyup="filterFingerprintTable()" style="max-width:380px;">
        </div>
        <div class="table-container">
          <table id="fingerprint-table">
            <thead>
              <tr>
                <th>Patient ID</th>
                <th>Normal Sample</th>
                <th>Tumor Sample</th>
                <th>Tumor vs Normal Concordance</th>
                <th>Plasma Timepoint / Sample</th>
                <th>Plasma vs Normal Concordance</th>
                <th>Plasma vs Tumor Concordance</th>
                <th>Min Intra Concordance</th>
                <th>Concordance Call</th>
              </tr>
            </thead>
            <tbody id="fingerprint-tbody">
              <!-- Rendered via JS -->
            </tbody>
          </table>
        </div>
      </div>

      <!-- Panel 2: Interactive Pairwise Concordance Heatmap Matrix -->
      <div class="panel">
        <div class="panel-header">
          <h2>Pairwise Concordance Heatmap Matrix</h2>
          <div style="display:flex; align-items:center; gap:1.25rem;">
            <div style="display:flex; align-items:center; gap:0.5rem; font-size:0.8rem; color:var(--text-secondary);">
              <span style="display:inline-block; width:12px; height:12px; background:#dcfce7; border:1px solid #86efac; border-radius:2px;"></span> Intra-Donor Match (&ge; 85%)
              <span style="display:inline-block; width:12px; height:12px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:2px; margin-left:8px;"></span> Unrelated Donors (&lt; 35%)
              <span style="display:inline-block; width:12px; height:12px; background:#fee2e2; border:1px solid #fca5a5; border-radius:2px; margin-left:8px;"></span> Mismatch Alert (&lt; 85% intra)
            </div>
            <select id="fp-matrix-filter" onchange="renderFingerprintMatrix()" style="padding:0.35rem 0.75rem; font-size:0.8rem; border-radius:6px; border:1px solid var(--border-color); background:white; color:var(--text-primary); cursor:pointer;">
              <option value="all">All Samples (Normal + Tumor + Plasma)</option>
              <option value="tn">Normal & Tumor Only</option>
              <option value="plasma">Plasma ctDNA Only</option>
            </select>
          </div>
        </div>
        <div class="table-container" style="max-height:640px; overflow:auto; position:relative;" id="fp-matrix-scroll-wrapper">
          <div id="fp-matrix-container">
            <!-- Heatmap rendered via JS -->
          </div>
        </div>
      </div>
    </div>

    <!-- GLOBAL APPLICATION MODALS & OVERLAYS (OUTSIDE ANY TAB-PANE) -->

    <!-- Oncoprint Cell Details Modal -->
    <div id="onco-cell-modal" onclick="if(event.target === this) closeOncoCellModal()" style="display:none; position:fixed; z-index:9999; left:0; top:0; width:100%; height:100%; background:rgba(15,23,42,0.6); backdrop-filter:blur(3px); align-items:center; justify-content:center;">
      <div style="background:#ffffff; border:1px solid var(--border-strong); border-radius:10px; width:90%; max-width:560px; padding:1.75rem; position:relative; box-shadow:0 16px 40px rgba(0,0,0,0.12); max-height:90vh; overflow-y:auto;">
        <button onclick="closeOncoCellModal()" style="position:absolute; top:1rem; right:1rem; background:transparent; border:none; color:var(--text-muted); font-size:1.4rem; cursor:pointer; line-height:1;">&times;</button>
        <div id="onco-modal-content"></div>
      </div>
    </div>

    <!-- Floating Interactive Oncoprint Hover Tooltip -->
    <div id="oncoprint-floating-tooltip" class="onco-floating-tooltip" style="display:none;"></div>

    <!-- Plasma ctDNA Multi-Modal Decomposition Modal -->
    <div id="plasma-modality-modal" onclick="if(event.target === this) closePlasmaModalityModal()" style="display:none; position:fixed; z-index:9999; left:0; top:0; width:100%; height:100%; background:rgba(15,23,42,0.6); backdrop-filter:blur(4px); align-items:center; justify-content:center; padding:1.5rem;">
      <div style="background:#ffffff; border:1px solid var(--border-strong); border-radius:10px; width:96%; max-width:1400px; max-height:92vh; display:flex; flex-direction:column; position:relative; box-shadow:0 20px 50px rgba(0,0,0,0.15); overflow:hidden;">
        <div style="padding:1rem 1.5rem; border-bottom:1px solid var(--border-color); display:flex; justify-content:space-between; align-items:center; background:var(--bg-subtle);">
          <div id="plasma-modal-header-info"></div>
          <div style="display:flex; gap:12px; align-items:center;">
            <a id="plasma-modal-download-btn" href="#" download class="btn btn-sm btn-primary" style="text-decoration:none; display:inline-flex; align-items:center; gap:6px;">⬇ Download Full-Res PNG</a>
            <button onclick="closePlasmaModalityModal()" style="background:transparent; border:none; color:var(--text-muted); font-size:1.6rem; cursor:pointer; line-height:1; padding:0 4px;">&times;</button>
          </div>
        </div>
        <div id="plasma-modal-body" style="padding:1.5rem; overflow-y:auto; flex:1;"></div>
      </div>
    </div>

    <!-- Patient Longitudinal MRD Trajectory Modal -->
    <div id="patient-trajectory-modal" onclick="if(event.target === this) closePatientTrajectoryModal()" style="display:none; position:fixed; z-index:9999; left:0; top:0; width:100%; height:100%; background:rgba(15,23,42,0.6); backdrop-filter:blur(4px); align-items:center; justify-content:center; padding:1.5rem;">
      <div style="background:#ffffff; border:1px solid var(--border-strong); border-radius:10px; width:96%; max-width:1400px; max-height:92vh; display:flex; flex-direction:column; position:relative; box-shadow:0 20px 50px rgba(0,0,0,0.15); overflow:hidden;">
        <div style="padding:1rem 1.5rem; border-bottom:1px solid var(--border-color); display:flex; justify-content:space-between; align-items:center; background:var(--bg-subtle);">
          <div id="patient-trajectory-header-info"></div>
          <div style="display:flex; gap:10px; align-items:center;">
            <a id="patient-trajectory-download-png" href="#" download class="btn btn-sm btn-primary" style="text-decoration:none; display:inline-flex; align-items:center; gap:6px;">⬇ Download PNG</a>
            <a id="patient-trajectory-download-pdf" href="#" download class="btn btn-sm" style="text-decoration:none; display:inline-flex; align-items:center; gap:6px;">⬇ Download PDF</a>
            <button onclick="closePatientTrajectoryModal()" style="background:transparent; border:none; color:var(--text-muted); font-size:1.6rem; cursor:pointer; line-height:1; padding:0 4px;">&times;</button>
          </div>
        </div>
        <div id="patient-trajectory-body" style="padding:1.5rem; overflow-y:auto; flex:1;"></div>
      </div>
    </div>

    <!-- HMF Orange Clinical PDF Report Modal -->
    <div id="orange-pdf-modal" onclick="if(event.target === this) closeOrangePdfModal()" style="display:none; position:fixed; z-index:10000; left:0; top:0; width:100%; height:100%; background:rgba(15,23,42,0.7); backdrop-filter:blur(4px); align-items:center; justify-content:center; padding:1rem;">
      <div style="background:#ffffff; border:1px solid var(--border-strong); border-radius:10px; width:98%; max-width:1450px; height:94vh; display:flex; flex-direction:column; position:relative; box-shadow:0 25px 50px -12px rgba(0,0,0,0.25); overflow:hidden;">
        <div style="padding:0.75rem 1.25rem; border-bottom:1px solid var(--border-color); display:flex; justify-content:space-between; align-items:center; background:var(--bg-subtle);">
          <div id="orange-pdf-modal-title" style="display:flex; align-items:center; gap:10px; flex-wrap:wrap;">
            <span style="font-weight:700; font-size:1.05rem; color:#0f172a;">📄 HMF Orange Clinical Molecular Report</span>
          </div>
          <div style="display:flex; gap:8px; align-items:center;">
            <button id="orange-pdf-zoom-out" type="button" onclick="zoomOrangePdf(-0.2)" class="btn btn-sm" style="font-size:0.75rem; padding:4px 8px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; cursor:pointer;" title="Zoom Out">🔍 -</button>
            <span id="orange-pdf-zoom-level" style="font-size:0.75rem; color:#64748b; font-weight:600; min-width:40px; text-align:center;">135%</span>
            <button id="orange-pdf-zoom-in" type="button" onclick="zoomOrangePdf(0.2)" class="btn btn-sm" style="font-size:0.75rem; padding:4px 8px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; cursor:pointer;" title="Zoom In">🔍 +</button>
            <button id="orange-pdf-newtab-btn" class="btn btn-sm" style="font-size:0.75rem; padding:4px 10px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; color:#0f172a; cursor:pointer; font-weight:600; display:inline-flex; align-items:center; gap:4px;" title="Open PDF in a full dedicated browser tab">
              ↗ Open in New Tab
            </button>
            <a id="orange-pdf-download-btn" href="#" download class="btn btn-sm btn-primary" style="font-size:0.75rem; padding:4px 10px; text-decoration:none; display:inline-flex; align-items:center; gap:4px;" title="Save PDF to local disk">
              ⬇ Download PDF
            </a>
            <button onclick="closeOrangePdfModal()" style="background:transparent; border:none; color:var(--text-muted); font-size:1.6rem; cursor:pointer; line-height:1; padding:0 4px; margin-left:6px;" title="Close modal">&times;</button>
          </div>
        </div>
        <div id="orange-pdf-scroll-area" style="flex:1; background:#525659; position:relative; overflow-y:auto; overflow-x:auto; display:flex; flex-direction:column; align-items:center; padding:1.5rem 1rem;">
          <iframe id="orange-pdf-frame" style="width:100%; height:100%; border:none; display:none;" src=""></iframe>
          <div id="orange-pdf-canvas-container" style="display:flex; flex-direction:column; align-items:center; width:100%;"></div>
        </div>
      </div>
    </div>

    <!-- 1X ichorCNA Genome-Wide PDF Report Modal -->
    <div id="ichor-pdf-modal" onclick="if(event.target === this) closeIchorPdfModal()" style="display:none; position:fixed; z-index:10000; left:0; top:0; width:100%; height:100%; background:rgba(15,23,42,0.7); backdrop-filter:blur(4px); align-items:center; justify-content:center; padding:1rem;">
      <div style="background:#ffffff; border:1px solid var(--border-strong); border-radius:10px; width:98%; max-width:1450px; height:94vh; display:flex; flex-direction:column; position:relative; box-shadow:0 25px 50px -12px rgba(0,0,0,0.25); overflow:hidden;">
        <div style="padding:0.75rem 1.25rem; border-bottom:1px solid var(--border-color); display:flex; justify-content:space-between; align-items:center; background:var(--bg-subtle);">
          <div id="ichor-pdf-modal-title" style="display:flex; align-items:center; gap:10px; flex-wrap:wrap;">
            <span style="font-weight:700; font-size:1.05rem; color:#0f172a;">📊 1X ichorCNA Genome-Wide Solution Plot</span>
            <span id="ichor-pdf-modal-sample-badge" class="badge badge-normal" style="font-weight:600; display:none;"></span>
            <span id="ichor-pdf-modal-patient-badge" class="badge badge-primary" style="font-size:0.8rem; display:none;"></span>
            <div id="ichor-modal-tab-group" style="display:none; align-items:center; gap:4px; margin-left:10px; background:#e2e8f0; padding:2px; border-radius:6px;">
              <button id="btn-ichor-tab-short" type="button" onclick="switchIchorModalTab('short')" style="font-size:0.75rem; padding:3px 10px; border-radius:4px; border:none; background:transparent; color:#475569; font-weight:600; cursor:pointer;" title="View in-silico short fragment size-selected (90-150 bp) 1X solution">⚡ Short Fragments (90–150 bp)</button>
              <button id="btn-ichor-tab-standard" type="button" onclick="switchIchorModalTab('standard')" style="font-size:0.75rem; padding:3px 10px; border-radius:4px; border:none; background:transparent; color:#475569; font-weight:600; cursor:pointer;" title="View standard unfiltered 1X solution">Standard 1X (All Reads)</button>
            </div>
            <span id="ichor-pdf-spinner" style="font-size:0.8rem; color:#64748b; margin-left:8px; display:inline-flex; align-items:center; gap:4px;">⏳ Rendering PDF document...</span>
          </div>
          <div style="display:flex; gap:8px; align-items:center;">
            <button id="ichor-pdf-zoom-out" type="button" onclick="zoomIchorPdf(-0.2)" class="btn btn-sm" style="font-size:0.75rem; padding:4px 8px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; cursor:pointer;" title="Zoom Out">🔍 -</button>
            <span id="ichor-pdf-zoom-level" style="font-size:0.75rem; color:#64748b; font-weight:600; min-width:40px; text-align:center;">135%</span>
            <button id="ichor-pdf-zoom-in" type="button" onclick="zoomIchorPdf(0.2)" class="btn btn-sm" style="font-size:0.75rem; padding:4px 8px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; cursor:pointer;" title="Zoom In">🔍 +</button>
            <button id="ichor-pdf-newtab-btn" class="btn btn-sm" style="font-size:0.75rem; padding:4px 10px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; color:#0f172a; cursor:pointer; font-weight:600; display:inline-flex; align-items:center; gap:4px;" title="Open PDF in a full dedicated browser tab">
              ↗ Open in New Tab
            </button>
            <a id="ichor-pdf-download-btn" href="#" download class="btn btn-sm btn-primary" style="font-size:0.75rem; padding:4px 10px; text-decoration:none; display:inline-flex; align-items:center; gap:4px; background:#d97706; border-color:#b45309;" title="Save PDF to local disk">
              ⬇ Download PDF
            </a>
            <button onclick="closeIchorPdfModal()" style="background:transparent; border:none; color:var(--text-muted); font-size:1.6rem; cursor:pointer; line-height:1; padding:0 4px; margin-left:6px;" title="Close modal">&times;</button>
          </div>
        </div>
        <div id="ichor-pdf-scroll-area" style="flex:1; background:#525659; position:relative; overflow-y:auto; overflow-x:auto; display:flex; flex-direction:column; align-items:center; padding:1.5rem 1rem;">
          <iframe id="ichor-pdf-frame" style="width:100%; height:100%; border:none; display:none;" src=""></iframe>
          <div id="ichor-pdf-canvas-container" style="display:flex; flex-direction:column; align-items:center; width:100%;"></div>
        </div>
      </div>
    </div>

  </main>

  <!-- EMBEDDED COHORT JSON PAYLOAD -->
  <script>
    const COHORT_DATA = {payload_json};
  </script>

  <script>
    // Tab Switching
    function switchTab(tabId) {{
      document.querySelectorAll('.tab-pane').forEach(el => el.classList.remove('active'));
      document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));

      const target = document.getElementById(tabId);
      if (target) target.classList.add('active');

      const btn = Array.from(document.querySelectorAll('.tab-btn')).find(b => b.getAttribute('onclick').includes(tabId));
      if (btn) btn.classList.add('active');

      if (tabId === 'tab-qc') renderQcCharts();
      if (tabId === 'tab-genomics') renderOncoprint();
      if (tabId === 'tab-mrd') renderMrdSection();
      if (tabId === 'tab-fingerprint') renderFingerprintSection();
    }}

    // CUPPA Cancer Type Color Scheme (Clean Pastel Palette for Clinical Taxonomy)
    function getCuppaStyle(cancerType) {{
      if (!cancerType || cancerType === 'Unknown' || cancerType === '--') {{
        return {{ bg: '#f1f5f9', border: '#cbd5e1', text: '#64748b', short: 'NA' }};
      }}
      const ct = cancerType.toLowerCase();
      if (ct.includes('breast')) {{
        return {{ bg: '#ffe4e6', border: '#f43f5e', text: '#881337', short: 'Breast' }};
      }} else if (ct.includes('endometri') || ct.includes('uter')) {{
        return {{ bg: '#ede9fe', border: '#8b5cf6', text: '#5b21b6', short: 'Endometrium' }};
      }} else if (ct.includes('ovar') || ct.includes('fallopian')) {{
        return {{ bg: '#fae8ff', border: '#d946ef', text: '#86198f', short: 'Ovary' }};
      }} else if (ct.includes('bone') || ct.includes('cartilag') || ct.includes('sarcoma') || ct.includes('soft tissue')) {{
        return {{ bg: '#fef3c7', border: '#f59e0b', text: '#92400e', short: 'Bone/ST' }};
      }} else if (ct.includes('head') || ct.includes('neck') || ct.includes('salivary') || ct.includes('adenoid')) {{
        return {{ bg: '#d1fae5', border: '#10b981', text: '#065f46', short: 'Head/Neck' }};
      }} else if (ct.includes('colorect') || ct.includes('colon') || ct.includes('rectum')) {{
        return {{ bg: '#e0f2fe', border: '#0ea5e9', text: '#0369a1', short: 'Colorectum' }};
      }} else if (ct.includes('lung')) {{
        return {{ bg: '#e2e8f0', border: '#64748b', text: '#1e293b', short: 'Lung' }};
      }} else if (ct.includes('prostate')) {{
        return {{ bg: '#ccfbf1', border: '#14b8a6', text: '#115e59', short: 'Prostate' }};
      }} else if (ct.includes('melanoma') || ct.includes('skin')) {{
        return {{ bg: '#ffedd5', border: '#f97316', text: '#9a3412', short: 'Melanoma' }};
      }} else if (ct.includes('kidney') || ct.includes('renal')) {{
        return {{ bg: '#e0e7ff', border: '#6366f1', text: '#3730a3', short: 'Kidney' }};
      }} else if (ct.includes('pancrea')) {{
        return {{ bg: '#fdf4ff', border: '#c026d3', text: '#701a75', short: 'Pancreas' }};
      }} else {{
        const words = cancerType.split(/[\\s_-]+/);
        const shortLabel = words[0].length > 10 ? words[0].substring(0, 8) + '..' : words[0];
        return {{ bg: '#f1f5f9', border: '#94a3b8', text: '#334155', short: shortLabel }};
      }}
    }}

    // Master Table & Overview
    function populateOverview() {{
      const tbody = document.getElementById('overview-tbody');
      tbody.innerHTML = '';

      let purities = [];
      let ploidies = [];
      let tmbs = [];
      let coverages = [];

      Object.values(COHORT_DATA.patients).forEach(p => {{
        const tumors = p.tumors || [];
        const plasmas = p.plasmas || [];
        const latestMrd = plasmas.length > 0 ? plasmas[plasmas.length - 1].mrd_status : 'No Plasma';

        if (tumors.length > 0) {{
          tumors.forEach(t => {{
            if (t.purity) purities.push(t.purity_pct);
            if (t.ploidy) ploidies.push(t.ploidy);
            if (t.tmb !== null && t.tmb !== undefined) tmbs.push(t.tmb);
            if (t.tumor_mean_coverage) coverages.push(t.tumor_mean_coverage);

            const mrdBadge = latestMrd === 'MRD Positive' ? '<span class="badge badge-mrd-pos">Positive</span>' :
                             (latestMrd === 'MRD Negative' ? '<span class="badge badge-mrd-neg">Negative</span>' : '<span style="color:var(--text-muted)">--</span>');

            const covStr = (t.tumor_mean_coverage ? t.tumor_mean_coverage.toFixed(1) + 'x (T)' : '--') +
                           ' / ' +
                           (t.normal_mean_coverage ? t.normal_mean_coverage.toFixed(1) + 'x (N)' : '--');

            const msiBadge = t.msi_status === 'MSI' || t.msi_status === 'MSI-H' ?
                             `<span class="badge badge-mrd-pos">${{t.msi_status}}</span>` :
                             `<span class="badge badge-normal">${{t.msi_status || 'MSS'}}</span>`;

            const wgdBadge = t.wgd ?
                             '<span class="badge badge-mrd-pos">Yes</span>' :
                             '<span class="badge badge-normal">No</span>';

            const trajBtn = plasmas.length > 0
              ? `<button class="btn btn-sm" onclick="showPatientTrajectoryModal('${{p.patient_id}}')" style="font-size:0.75rem; padding:3px 8px; border-radius:4px; border:1px solid #10b981; background:#ffffff; color:#059669; cursor:pointer; font-weight:600; white-space:nowrap; display:inline-flex; align-items:center; gap:3px;">📈 Trajectory</button>`
              : '<span style="color:var(--text-muted)">--</span>';

            const hasOrange = Boolean(
              (COHORT_DATA.orange_reports_b64 && (COHORT_DATA.orange_reports_b64[t.tumor_id] || COHORT_DATA.orange_reports_b64[t.raw_tumor_id])) ||
              t.orange_pdf_b64 ||
              t.orange_pdf_rel ||
              t.has_orange_pdf
            );

            const orangeBtn = hasOrange
              ? `<button type="button" onclick="showOrangePdfModal('${{t.tumor_id}}', '${{p.patient_id}}')" class="badge" style="background:#fef3c7; color:#92400e; border:1px solid #f59e0b; cursor:pointer; margin-left:6px; font-size:0.68rem; padding:2px 7px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:3px;" title="View native HMF Orange Clinical PDF Report directly in browser">📄 Report</button>`
              : '';

            const hasIchor = Boolean(
              (COHORT_DATA.ichor_reports_b64 && (
                COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{t.tumor_id}}`] ||
                COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{t.tumor_id}}_standard`] ||
                COHORT_DATA.ichor_reports_b64[t.tumor_id] ||
                COHORT_DATA.ichor_reports_b64[t.raw_tumor_id]
              )) ||
              t.ichor_pdf_b64 ||
              t.ichor_pdf_rel ||
              t.has_ichor_pdf ||
              t.has_ichor_standard ||
              t.has_ichor_short
            );

            const ichorBtn = hasIchor
              ? `<button type="button" onclick="showIchorPdfModal('${{t.tumor_id}}', '${{p.patient_id}}')" class="badge" style="background:#fef3c7; color:#b45309; border:1px solid #f59e0b; cursor:pointer; margin-left:4px; font-size:0.68rem; padding:2px 6px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:2px;" title="View 1X ichorCNA Genome-Wide Solution Plot in browser">📊 Plot</button>`
              : '';

            const cType = t.cuppa_cancer_type || '--';
            const lhood = t.cuppa_likelihood;
            const cStyle = getCuppaStyle(cType);
            const lhoodStr = (lhood !== null && lhood !== undefined && cType !== '--') ? ` <span style="opacity:0.8; font-weight:normal; font-size:0.7rem;">(${{(lhood * 100).toFixed(0)}}%)</span>` : '';
            const cuppaBadge = cType !== '--'
              ? `<span class="badge" style="background:${{cStyle.bg}}; color:${{cStyle.text}}; border:1px solid ${{cStyle.border}}; font-weight:600;" title="CUPPA Predicted: ${{cType}}${{lhood !== null ? ' (' + (lhood * 100).toFixed(1) + '%)' : ''}}">${{cType}}${{lhoodStr}}</span>`
              : '<span style="color:var(--text-muted)">--</span>';

            let purityBadge = `${{t.purity_pct.toFixed(1)}}%`;
            if (t.curated_purity !== null && t.curated_purity !== undefined) {{
              purityBadge += ` <span class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #bae6fd; font-size:0.65rem;" title="Curated via ${{t.curation_method}}: ${{t.curation_notes}}">CURATED</span>`;
            }} else if (t.status === 'NORMAL' || t.status === 'WARN_LOW_PURITY') {{
              purityBadge += ` <span class="badge badge-normal" style="font-size:0.65rem;" title="Low cellularity baseline: copy-neutral genome / sub-threshold CNV">LOW-CELL</span>`;
            }}

            const fpData = COHORT_DATA.fingerprints || {{}};
            const ps = (fpData.patient_summaries || {{}})[p.patient_id];
            let fpBadge = '<span style="color:var(--text-muted)">--</span>';
            if (ps) {{
              if (ps.is_confirmed) {{
                const minC = (ps.min_intra_concordance !== null && ps.min_intra_concordance !== undefined) ? ps.min_intra_concordance.toFixed(1) + '%' : '100%';
                fpBadge = `<span class="badge" style="background:#ecfdf5; color:#065f46; border:1px solid #a7f3d0;" title="Normal, Tumor, and Plasma(s) confirmed same individual (${{minC}} intra-donor concordance)">✓ Confirmed (${{minC}})</span>`;
              }} else if (ps.status === 'SINGLE_SAMPLE_NO_COMPARATOR') {{
                fpBadge = `<span class="badge badge-normal" title="Single sample available, no comparator">Single Sample</span>`;
              }} else {{
                const minC = (ps.min_intra_concordance !== null && ps.min_intra_concordance !== undefined) ? ps.min_intra_concordance.toFixed(1) + '%' : '--';
                fpBadge = `<span class="badge" style="background:#fee2e2; color:#991b1b; border:1px solid #fecaca;" title="Sample mismatch suspected! Min intra-donor concordance: ${{minC}}">⚠ Mismatch (${{minC}})</span>`;
              }}
            }}

            // 1. Primary tumor ichor button and label
            const hasTumorIchor = Boolean(
              t.ichor_pdf_b64 || t.ichor_pdf_rel || t.has_ichor_pdf ||
              t.has_ichor_standard || t.has_ichor_short ||
              (COHORT_DATA.ichor_reports_b64 && (
                COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{t.tumor_id}}`] ||
                COHORT_DATA.ichor_reports_b64[t.tumor_id] ||
                (t.raw_tumor_id && COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{t.raw_tumor_id}}`]) ||
                (t.raw_tumor_id && COHORT_DATA.ichor_reports_b64[t.raw_tumor_id])
              ))
            );
            const tumorTfVal = (t.ichor_tf !== null && t.ichor_tf !== undefined) ? `${{(t.ichor_tf * 100).toFixed(1)}}%` : '';
            const tumorIchorBtn = hasTumorIchor
              ? `<button type="button" onclick="showIchorPdfModal('${{t.tumor_id}}', '${{p.patient_id}}')" class="badge" style="background:#fef3c7; color:#b45309; border:1px solid #f59e0b; cursor:pointer; margin-left:3px; font-size:0.68rem; padding:2px 6px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:2px;" title="View Primary Tumor 1X ichorCNA Solution Plot">📊 Tumor Plot</button>`
              : '';

            let tumorHtml = '';
            if (tumorTfVal || tumorIchorBtn) {{
              tumorHtml = `<div style="white-space:nowrap; margin-bottom:2px;"><strong>Tumor:</strong> ${{tumorTfVal}} ${{tumorIchorBtn}}</div>`;
            }}

            // 2. Plasma ichor buttons and labels (all timepoints)
            const plIchorRows = [];
            plasmas.forEach(pl => {{
              const hasPlIchor = Boolean(
                pl.ichor_pdf_b64 || pl.ichor_pdf_rel || pl.has_ichor_pdf ||
                pl.has_ichor_short || pl.has_ichor_standard ||
                (COHORT_DATA.ichor_reports_b64 && (
                  COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{pl.sample_id}}`] ||
                  COHORT_DATA.ichor_reports_b64[pl.sample_id] ||
                  (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{pl.clean_sample_name}}`]) ||
                  (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[pl.clean_sample_name])
                ))
              );
              const plTfVal = (pl.ichor_tf !== null && pl.ichor_tf !== undefined) ? `${{(pl.ichor_tf * 100).toFixed(1)}}%` : '';
              const plIchorBtn = hasPlIchor
                ? `<button type="button" onclick="showIchorPdfModal('${{pl.sample_id}}', '${{p.patient_id}}')" class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc; cursor:pointer; margin-left:3px; font-size:0.68rem; padding:2px 6px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:2px;" title="View Plasma 1X ichorCNA Solution Plot (${{pl.sample_id}})">📊 Plasma Plot</button>`
                : '';
              if (plTfVal || plIchorBtn) {{
                const plLabel = pl.timepoint || pl.clean_sample_name || pl.sample_id;
                plIchorRows.push(`<div style="font-size:0.75rem; color:#0369a1; margin-top:2px; white-space:nowrap;"><strong>Plasma (${{plLabel}}):</strong> ${{plTfVal}} ${{plIchorBtn}}</div>`);
              }}
            }});

            let ichorCellHtml = '<span style="color:var(--text-muted)">--</span>';
            if (tumorHtml || plIchorRows.length > 0) {{
              ichorCellHtml = `<div>${{tumorHtml}}${{plIchorRows.join('')}}</div>`;
            }}

            const row = document.createElement('tr');
            row.innerHTML = `
              <td><strong>${{p.patient_id}}</strong></td>
              <td><span class="badge badge-normal">${{t.tumor_id}}</span>${{orangeBtn}}</td>
              <td>${{cuppaBadge}}</td>
              <td>${{purityBadge}}</td>
              <td>${{t.ploidy.toFixed(2)}}</td>
              <td>${{covStr}}</td>
              <td>${{t.tmb !== null && t.tmb !== undefined ? t.tmb.toFixed(2) : '--'}}</td>
              <td>${{msiBadge}}</td>
              <td>${{wgdBadge}}</td>
              <td>${{ichorCellHtml}}</td>
              <td>${{fpBadge}}</td>
              <td>${{plasmas.length}}</td>
              <td>${{mrdBadge}}</td>
              <td>${{trajBtn}}</td>
              <td style="font-size:0.8rem; color:var(--accent-teal);">${{t.altered_genes.slice(0, 5).join(', ') || '--'}}</td>
            `;
            tbody.appendChild(row);
          }});
        }} else {{
          const trajBtn = plasmas.length > 0
            ? `<button class="btn btn-sm" onclick="showPatientTrajectoryModal('${{p.patient_id}}')" style="font-size:0.75rem; padding:3px 8px; border-radius:4px; border:1px solid #10b981; background:#ffffff; color:#059669; cursor:pointer; font-weight:600; white-space:nowrap; display:inline-flex; align-items:center; gap:3px;">📈 Trajectory</button>`
            : '<span style="color:var(--text-muted)">--</span>';

          const fpData = COHORT_DATA.fingerprints || {{}};
          const ps = (fpData.patient_summaries || {{}})[p.patient_id];
          let fpBadge = '<span style="color:var(--text-muted)">--</span>';
          if (ps) {{
            if (ps.is_confirmed) {{
              const minC = (ps.min_intra_concordance !== null && ps.min_intra_concordance !== undefined) ? ps.min_intra_concordance.toFixed(1) + '%' : '100%';
              fpBadge = `<span class="badge" style="background:#ecfdf5; color:#065f46; border:1px solid #a7f3d0;">✓ Confirmed (${{minC}})</span>`;
            }} else {{
              const minC = (ps.min_intra_concordance !== null && ps.min_intra_concordance !== undefined) ? ps.min_intra_concordance.toFixed(1) + '%' : '--';
              fpBadge = `<span class="badge" style="background:#fee2e2; color:#991b1b; border:1px solid #fecaca;">⚠ Mismatch (${{minC}})</span>`;
            }}
          }}

          const plIchorRows = [];
          plasmas.forEach(pl => {{
            const hasPlIchor = Boolean(
              pl.ichor_pdf_b64 || pl.ichor_pdf_rel || pl.has_ichor_pdf ||
              pl.has_ichor_short || pl.has_ichor_standard ||
              (COHORT_DATA.ichor_reports_b64 && (
                COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{pl.sample_id}}`] ||
                COHORT_DATA.ichor_reports_b64[pl.sample_id] ||
                (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[`${{p.patient_id}}_${{pl.clean_sample_name}}`]) ||
                (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[pl.clean_sample_name])
              ))
            );
            const plTfVal = (pl.ichor_tf !== null && pl.ichor_tf !== undefined) ? `${{(pl.ichor_tf * 100).toFixed(1)}}%` : '';
            const plIchorBtn = hasPlIchor
              ? `<button type="button" onclick="showIchorPdfModal('${{pl.sample_id}}', '${{p.patient_id}}')" class="badge" style="background:#e0f2fe; color:#0369a1; border:1px solid #7dd3fc; cursor:pointer; margin-left:3px; font-size:0.68rem; padding:2px 6px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:2px;" title="View Plasma 1X ichorCNA Solution Plot (${{pl.sample_id}})">📊 Plasma Plot</button>`
              : '';
            if (plTfVal || plIchorBtn) {{
              const plLabel = pl.timepoint || pl.clean_sample_name || pl.sample_id;
              plIchorRows.push(`<div style="font-size:0.75rem; color:#0369a1; margin-top:2px; white-space:nowrap;"><strong>Plasma (${{plLabel}}):</strong> ${{plTfVal}} ${{plIchorBtn}}</div>`);
            }}
          }});
          const noTumorIchorHtml = plIchorRows.length > 0 ? `<div>${{plIchorRows.join('')}}</div>` : '<span style="color:var(--text-muted)">--</span>';

          const row = document.createElement('tr');
          row.innerHTML = `
            <td><strong>${{p.patient_id}}</strong></td>
            <td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td>
            <td>${{noTumorIchorHtml}}</td>
            <td>${{fpBadge}}</td>
            <td>${{plasmas.length}}</td>
            <td>--</td>
            <td>${{trajBtn}}</td>
            <td>--</td>
          `;
          tbody.appendChild(row);
        }}
      }});

      // Update KPIs
      if (purities.length > 0) {{
        purities.sort((a,b)=>a-b);
        document.getElementById('kpi-purity').innerText = purities[Math.floor(purities.length/2)].toFixed(1) + '%';
      }}
      if (ploidies.length > 0) {{
        ploidies.sort((a,b)=>a-b);
        document.getElementById('kpi-ploidy').innerText = ploidies[Math.floor(ploidies.length/2)].toFixed(2);
      }}
      if (coverages.length > 0) {{
        coverages.sort((a,b)=>a-b);
        document.getElementById('kpi-coverage').innerText = coverages[Math.floor(coverages.length/2)].toFixed(1) + 'x';
      }}
      if (tmbs.length > 0) {{
        tmbs.sort((a,b)=>a-b);
        document.getElementById('kpi-tmb').innerText = tmbs[Math.floor(tmbs.length/2)].toFixed(2);
      }}
    }}

    function filterOverviewTable() {{
      const q = document.getElementById('overview-search').value.toLowerCase();
      document.querySelectorAll('#overview-tbody tr').forEach(r => {{
        r.style.display = r.innerText.toLowerCase().includes(q) ? '' : 'none';
      }});
    }}

    // QC Charts
    let qcCharts = {{}};
    function renderQcCharts() {{
      if (qcCharts.purityPloidy) return; // already rendered

      // 1. Purity vs Ploidy
      const ctx1 = document.getElementById('chart-purity-ploidy').getContext('2d');
      const pts1 = COHORT_DATA.all_tumors.map(t => ({{ x: t.ploidy, y: t.purity_pct, label: t.tumor_id }}));
      qcCharts.purityPloidy = new Chart(ctx1, {{
        type: 'scatter',
        data: {{
          datasets: [{{
            label: 'Primary Tumors',
            data: pts1,
            backgroundColor: '#0284c7',
            borderColor: '#0369a1',
            pointRadius: 6,
            pointHoverRadius: 9
          }}]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          scales: {{
            x: {{ title: {{ display: true, text: 'Tumor Ploidy', color: '#475569', font: {{ weight: '600' }} }}, grid: {{ color: '#f1f5f9' }}, ticks: {{ color: '#64748b' }} }},
            y: {{ title: {{ display: true, text: 'Tumor Purity (%)', color: '#475569', font: {{ weight: '600' }} }}, min: 0, max: 100, grid: {{ color: '#f1f5f9' }}, ticks: {{ color: '#64748b' }} }}
          }}
        }}
      }});

      // 2. TMB distribution
      const ctx2 = document.getElementById('chart-tmb').getContext('2d');
      const tmbs = COHORT_DATA.all_tumors.filter(t => t.tmb !== null && t.tmb !== undefined).map(t => ({{
        x: t.tumor_id,
        y: Math.max(t.tmb, 0.01),
        raw: t.tmb,
        msi: t.msi_status || 'MSS'
      }}));

      qcCharts.tmb = new Chart(ctx2, {{
        type: 'bar',
        data: {{
          labels: tmbs.map(t => t.x),
          datasets: [{{
            label: 'TMB (mut/Mb)',
            data: tmbs.map(t => t.y),
            backgroundColor: tmbs.map(t => (t.msi === 'MSI' || t.msi === 'MSI-H') ? '#a855f7' : '#f43f5e'),
            borderRadius: 4
          }}]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          plugins: {{
            tooltip: {{
              callbacks: {{
                label: function(ctx) {{
                  const pt = tmbs[ctx.dataIndex];
                  return `TMB: ${{pt.raw.toFixed(2)}} mut/Mb (${{pt.msi}})`;
                }}
              }}
            }}
          }},
          scales: {{
            y: {{
              type: 'logarithmic',
              title: {{ display: true, text: 'TMB (Mut/Mb) - Log Scale', color: '#475569', font: {{ weight: '600' }} }},
              grid: {{ color: '#f1f5f9' }},
              ticks: {{ color: '#64748b' }}
            }},
            x: {{
              grid: {{ display: false }},
              ticks: {{ color: '#64748b', font: {{ size: 10 }}, maxRotation: 45 }}
            }}
          }}
        }}
      }});

      // 3. WGS Sequencing Coverage (Tumor, Plasma ctDNA & Matched Normal)
      const ctxCov = document.getElementById('chart-coverage').getContext('2d');
      const covTumors = COHORT_DATA.all_tumors.filter(t => t.tumor_mean_coverage || t.normal_mean_coverage);
      const covLabels = covTumors.map(t => t.tumor_id);
      const tumorCovData = covTumors.map(t => t.tumor_mean_coverage || 0);
      const normalCovData = covTumors.map(t => t.normal_mean_coverage || 0);
      const plasmaCovData = covTumors.map(t => {{
        const p = COHORT_DATA.patients[t.patient_id];
        if (p && p.plasmas && p.plasmas.length > 0) {{
          return p.plasmas[0].weighted_depth || 0;
        }}
        return 0;
      }});

      const hasPlasmaCov = plasmaCovData.some(v => v > 0);
      const covDatasets = [
        {{
          label: 'Primary Tumor Coverage (Mean Depth)',
          data: tumorCovData,
          backgroundColor: '#0ea5e9',
          borderRadius: 4
        }}
      ];
      if (hasPlasmaCov) {{
        covDatasets.push({{
          label: 'Plasma ctDNA Coverage (Weighted Depth)',
          data: plasmaCovData,
          backgroundColor: '#a855f7',
          borderRadius: 4
        }});
      }}
      covDatasets.push({{
        label: 'Matched Normal Coverage (Mean Depth)',
        data: normalCovData,
        backgroundColor: '#94a3b8',
        borderRadius: 4
      }});

      qcCharts.coverage = new Chart(ctxCov, {{
        type: 'bar',
        data: {{
          labels: covLabels,
          datasets: covDatasets
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          plugins: {{
            legend: {{ labels: {{ color: '#334155' }} }},
            tooltip: {{
              callbacks: {{
                label: function(ctx) {{
                  return `${{ctx.dataset.label}}: ${{ctx.raw ? ctx.raw.toFixed(1) + 'x' : 'NA'}}`;
                }}
              }}
            }}
          }},
          scales: {{
            y: {{
              title: {{ display: true, text: 'Mean Coverage Depth (X)', color: '#475569', font: {{ weight: '600' }} }},
              grid: {{ color: '#f1f5f9' }},
              ticks: {{ color: '#64748b' }},
              beginAtZero: true
            }},
            x: {{
              grid: {{ display: false }},
              ticks: {{ color: '#64748b', font: {{ size: 10 }}, maxRotation: 45 }}
            }}
          }}
        }}
      }});

      // 4. Consequence breakdown (Somatic Drivers & High-Confidence Focal CNAs)
      const ctx3 = document.getElementById('chart-consequences').getContext('2d');
      const counts = {{ 'Missense': 0, 'Frameshift': 0, 'Nonsense': 0, 'Splice': 0, 'In-frame InDel': 0, 'Amplification': 0, 'Deep Deletion': 0 }};
      COHORT_DATA.all_tumors.forEach(t => {{
        const countedAlts = new Set();
        (t.drivers || []).forEach(d => {{
          if (counts[d.alt_type] !== undefined) {{
            counts[d.alt_type]++;
            countedAlts.add(`${{d.gene}}:${{d.alt_type}}`);
          }}
        }});
        Object.entries(t.cnvs || {{}}).forEach(([gene, c]) => {{
          if (c.is_driver && (c.state === 'Amplification' || c.state === 'Deep Deletion') && counts[c.state] !== undefined) {{
            if (!countedAlts.has(`${{gene}}:${{c.state}}`)) {{
              counts[c.state]++;
              countedAlts.add(`${{gene}}:${{c.state}}`);
            }}
          }}
        }});
      }});

      qcCharts.consequences = new Chart(ctx3, {{
        type: 'doughnut',
        data: {{
          labels: Object.keys(counts),
          datasets: [{{
            data: Object.values(counts),
            backgroundColor: ['#0d9488', '#9333ea', '#e11d48', '#d97706', '#0284c7', '#dc2626', '#2563eb']
          }}]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          plugins: {{
            legend: {{ position: 'right', labels: {{ color: '#334155', font: {{ size: 11 }} }} }},
            tooltip: {{
              callbacks: {{
                label: function(ctx) {{
                  const val = ctx.raw || 0;
                  const total = Object.values(counts).reduce((a, b) => a + b, 0);
                  const pct = total > 0 ? ((val / total) * 100).toFixed(1) : 0;
                  return `${{ctx.label}}: ${{val}} (${{pct}}%)`;
                }}
              }}
            }}
          }}
        }}
      }});

      // 5. ichorCNA concordance (Tumor PURPLE & Plasma WISP/PURPLE)
      const ctx4 = document.getElementById('chart-ichor-concordance').getContext('2d');
      const tumorConcPts = COHORT_DATA.all_tumors
        .filter(t => t.ichor_tf !== null && t.ichor_tf !== undefined)
        .map(t => ({{ x: t.purity_pct, y: t.ichor_tf * 100, label: `${{t.patient_id}} (${{t.tumor_id}})`, sample_id: t.tumor_id, patient_id: t.patient_id }}));
      const plasmaConcPts = (COHORT_DATA.all_plasmas || [])
        .filter(p => p.ichor_tf !== null && p.ichor_tf !== undefined)
        .map(p => {{
          const xVal = (p.wisp_purity_pct !== null && p.wisp_purity_pct !== undefined && p.wisp_purity_pct > 0)
            ? p.wisp_purity_pct
            : (p.purple_purity_pct !== null && p.purple_purity_pct !== undefined ? p.purple_purity_pct : (p.snv_purity_pct || 0));
          return {{ x: xVal, y: p.ichor_tf * 100, label: `${{p.patient_id}} ${{p.timepoint}} (${{p.sample_id}})`, sample_id: p.sample_id, patient_id: p.patient_id }};
        }});

      qcCharts.ichor = new Chart(ctx4, {{
        type: 'scatter',
        data: {{
          datasets: [
            {{
              label: 'Primary Tumors (PURPLE vs ichorCNA)',
              data: tumorConcPts,
              backgroundColor: '#0284c7',
              borderColor: '#0369a1',
              pointRadius: 6,
              pointHoverRadius: 8
            }},
            {{
              label: 'Plasma ctDNA (WISP/PURPLE vs ichorCNA)',
              data: plasmaConcPts,
              backgroundColor: '#9333ea',
              borderColor: '#7e22ce',
              pointRadius: 6,
              pointHoverRadius: 8
            }}
          ]
        }},
        options: {{
          responsive: true,
          maintainAspectRatio: false,
          onClick: function(evt, activeEls) {{
            if (activeEls && activeEls.length > 0) {{
              const el = activeEls[0];
              const dataset = qcCharts.ichor.data.datasets[el.datasetIndex];
              const pt = dataset && dataset.data ? dataset.data[el.index] : null;
              if (pt && pt.sample_id) {{
                showIchorPdfModal(pt.sample_id, pt.patient_id);
              }}
            }}
          }},
          onHover: function(evt, activeEls) {{
            if (evt && evt.native && evt.native.target) {{
              evt.native.target.style.cursor = (activeEls && activeEls.length > 0) ? 'pointer' : 'default';
            }}
          }},
          plugins: {{
            legend: {{ labels: {{ color: '#334155' }} }},
            tooltip: {{
              callbacks: {{
                label: function(ctx) {{
                  const pt = ctx.raw;
                  return `${{pt.label || ''}}: Expected ${{ctx.parsed.x.toFixed(2)}}%, ichorCNA ${{ctx.parsed.y.toFixed(2)}}%`;
                }},
                afterLabel: function(ctx) {{
                  return '💡 Click dot to view 1X ichorCNA solution plot';
                }}
              }}
            }}
          }},
          scales: {{
            x: {{ title: {{ display: true, text: 'Expected Tumor Fraction / Purity (%)', color: '#475569', font: {{ weight: '600' }} }}, grid: {{ color: '#f1f5f9' }}, ticks: {{ color: '#64748b' }} }},
            y: {{ title: {{ display: true, text: '1X ichorCNA Tumor Fraction (%)', color: '#475569', font: {{ weight: '600' }} }}, grid: {{ color: '#f1f5f9' }}, ticks: {{ color: '#64748b' }} }}
          }}
        }}
      }});
    }}

    // Oncoprint conforming to cBioPortal Visual Hierarchy & Taxonomy
    let oncoGenes = COHORT_DATA.oncoprint.genes;

    function getDynamicCellDetails(gene, tumorId) {{
      const sm = (COHORT_DATA.oncoprint && COHORT_DATA.oncoprint.sample_metadata && COHORT_DATA.oncoprint.sample_metadata[tumorId]) ? COHORT_DATA.oncoprint.sample_metadata[tumorId] : {{}};
      const t = (sm.patient_id ? COHORT_DATA.all_tumors.find(x => x.patient_id === sm.patient_id && (x.tumor_id === sm.raw_sample_id || x.tumor_id === tumorId)) : null) ||
                COHORT_DATA.all_tumors.find(x => x.tumor_id === tumorId) ||
                (sm.patient_id ? (COHORT_DATA.all_plasmas || []).find(x => x.patient_id === sm.patient_id && (x.sample_id === sm.raw_sample_id || x.sample_id === tumorId)) : null) ||
                (COHORT_DATA.all_plasmas || []).find(x => x.sample_id === tumorId) || {{}};
      const cnv = (t.cnvs && t.cnvs[gene]) ? t.cnvs[gene] : null;
      let cna_type = null;
      let cna_label = null;
      if (cnv) {{
        if (cnv.cbio_val === 2) {{ cna_type = 'AMP'; cna_label = 'Amplification'; }}
        else if (cnv.cbio_val === -2) {{ cna_type = 'HOMDEL'; cna_label = 'Deep Deletion'; }}
        else if (cnv.cbio_val === 1) {{ cna_type = 'GAIN'; cna_label = 'Gain'; }}
        else if (cnv.cbio_val === -1) {{ cna_type = 'HETLOSS'; cna_label = 'Shallow Deletion'; }}
      }}

      const muts = [];
      (t.drivers || []).filter(d => d.gene === gene).forEach(d => {{
        muts.push({{
          type: d.alt_type,
          subtype: d.subtype || d.alt_type,
          hgvs_protein: d.hgvs_protein,
          hgvs_coding: d.hgvs_coding,
          vaf: d.vaf,
          is_germline: false
        }});
      }});
      (t.germline_drivers || []).filter(gd => gd.gene === gene).forEach(gd => {{
        muts.push({{
          type: gd.alt_type,
          subtype: gd.subtype || gd.alt_type,
          hgvs_protein: gd.hgvs_protein,
          hgvs_coding: gd.hgvs_coding,
          vaf: gd.vaf,
          is_germline: true
        }});
      }});

      const fusions = (t.fusions && t.fusions[gene]) ? t.fusions[gene] : [];
      const fusion_name = fusions.length > 0 ? fusions[0].name : '';

      const summaryTags = [];
      if (cna_label) summaryTags.push(cna_label);
      muts.forEach(m => {{
        let desc = m.type;
        const parts = [];
        if (m.is_germline) parts.push('Germline');
        if (m.vaf !== undefined && m.vaf !== null && m.vaf > 0) {{
          const vafPct = m.vaf > 1.0 ? m.vaf : m.vaf * 100.0;
          parts.push('VAF: ' + vafPct.toFixed(1) + '%');
        }}
        const suffix = parts.length > 0 ? (' (' + parts.join(', ') + ')') : '';
        if (m.hgvs_protein) {{
          desc = m.hgvs_protein + ' (' + m.type + (parts.length > 0 ? ', ' + parts.join(', ') : '') + ')';
        }} else if (m.hgvs_coding) {{
          desc = m.hgvs_coding + ' (' + m.type + (parts.length > 0 ? ', ' + parts.join(', ') : '') + ')';
        }} else {{
          desc = m.type + suffix;
        }}
        summaryTags.push(desc);
      }});
      if (fusion_name) summaryTags.push('Fusion (' + fusion_name + ')');

      return {{
        gene: gene,
        tumor_id: tumorId,
        patient_id: t.patient_id || sm.patient_id,
        cna: cna_type,
        cna_label: cna_label,
        mutations: muts,
        is_germline: muts.some(m => m.is_germline),
        fusion: Boolean(fusion_name),
        fusion_name: fusion_name,
        expression: null,
        is_profiled: true,
        summary: summaryTags.length > 0 ? summaryTags.join('; ') : 'Wild Type',
        cna_min_cn: cnv ? cnv.min_cn : null,
        cna_max_cn: cnv ? cnv.max_cn : null,
      }};
    }}

    function renderOncoprint(limit = 20) {{
      const container = document.getElementById('oncoprint-container');
      const tumors = COHORT_DATA.oncoprint.tumors;
      const cellDetails = COHORT_DATA.oncoprint.cell_details || {{}};
      const patientBlocks = COHORT_DATA.oncoprint.patient_blocks || [];
      const patientBoundaries = COHORT_DATA.oncoprint.patient_boundaries || [];
      const sampleMeta = COHORT_DATA.oncoprint.sample_metadata || {{}};

      // Populate CUPPA Cohort Distribution Summary Pills
      const cuppaCounts = {{}};
      let totalCuppaTumors = 0;
      (COHORT_DATA.all_tumors || []).forEach(t => {{
        if (t.cuppa_cancer_type) {{
          cuppaCounts[t.cuppa_cancer_type] = (cuppaCounts[t.cuppa_cancer_type] || 0) + 1;
          totalCuppaTumors++;
        }}
      }});
      const cuppaTotalEl = document.getElementById('cuppa-cohort-total');
      if (cuppaTotalEl) {{
        cuppaTotalEl.textContent = `${{totalCuppaTumors}} Primary Tumor${{totalCuppaTumors !== 1 ? 's' : ''}} Profiled`;
      }}
      const cuppaPillsEl = document.getElementById('cuppa-distribution-pills');
      if (cuppaPillsEl) {{
        cuppaPillsEl.innerHTML = '';
        const sortedCuppa = Object.entries(cuppaCounts).sort((a, b) => b[1] - a[1]);
        if (sortedCuppa.length === 0) {{
          cuppaPillsEl.innerHTML = '<span style="color:var(--text-muted); font-size:0.8rem;">No CUPPA predictions available</span>';
        }} else {{
          sortedCuppa.forEach(([cType, count]) => {{
            const cStyle = getCuppaStyle(cType);
            const pct = totalCuppaTumors > 0 ? ((count / totalCuppaTumors) * 100).toFixed(0) : 0;
            const pill = document.createElement('div');
            pill.className = 'badge';
            pill.style.cssText = `background:${{cStyle.bg}}; color:${{cStyle.text}}; border:1px solid ${{cStyle.border}}; padding:3px 9px; border-radius:12px; font-size:0.75rem; font-weight:600; display:inline-flex; align-items:center; gap:6px;`;
            pill.innerHTML = `<span>${{cType}}</span><strong style="background:rgba(0,0,0,0.06); padding:1px 5px; border-radius:8px; font-size:0.7rem;">${{count}} (${{pct}}%)</strong>`;
            cuppaPillsEl.appendChild(pill);
          }});
        }}
      }}

      const genesToShow = oncoGenes.slice(0, limit);

      let html = `<table class="oncoprint-table"><thead>`;

      // 1. Patient Group Header Row
      if (patientBlocks && patientBlocks.length > 0) {{
        html += `<tr class="patient-header-row"><th class="gene-header" style="background:var(--bg-subtle); border-bottom:2px solid var(--border-color); color:var(--text-secondary);">Patient</th>`;
        patientBlocks.forEach((pb, idx) => {{
          const isLast = (idx === patientBlocks.length - 1);
          const bCls = isLast ? '' : ' patient-group-boundary';
          html += `<th colspan="${{pb.count}}" class="patient-group-header${{bCls}}">${{pb.patient_id}}</th>`;
        }});
        html += `<th class="freq-header" style="background:var(--bg-subtle); border-bottom:2px solid var(--border-color);"></th></tr>`;
      }}

      // 2. Sample Label Header Row
      html += `<tr class="sample-header-row"><th class="gene-header">Sample</th>`;
      tumors.forEach(sid => {{
        const sm = sampleMeta[sid] || {{}};
        const isBoundary = patientBoundaries.includes(sid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const isPlasma = sm.is_plasma;
        const tp = sm.timepoint || (isPlasma ? 'Plasma' : 'Primary');
        const badgeCls = isPlasma ? 'type-plasma' : 'type-primary';
        const displayId = sm.clean_sample_name || sid;
        html += `<th class="sample-header${{bCls}}">` +
                `<div class="sample-label" title="${{sid}}">${{displayId}}</div>` +
                `<div class="sample-type-badge ${{badgeCls}}">${{tp}}</div>` +
                `</th>`;
      }});
      html += `<th class="freq-header">Freq</th></tr></thead><tbody>`;

      // 3. Clinical Track: Sample Type
      html += `<tr class="onco-track-row"><td class="onco-track-label">Sample Type</td>`;
      tumors.forEach(sid => {{
        const sm = sampleMeta[sid] || {{}};
        const isBoundary = patientBoundaries.includes(sid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const isPlasma = sm.is_plasma;
        const shortLabel = isPlasma ? 'Plasma' : 'Tumor';
        const typeBg = isPlasma ? '#f3e8ff' : '#e0f2fe';
        const typeColor = isPlasma ? '#6b21a8' : '#0369a1';
        const typeBorder = isPlasma ? '#e9d5ff' : '#bae6fd';
        html += `<td class="${{bCls}}"><div class="onco-track-cell" style="background:${{typeBg}}; color:${{typeColor}}; border:1px solid ${{typeBorder}};" title="${{sid}}: ${{sm.sample_type || shortLabel}}">${{shortLabel}}</div></td>`;
      }});
      html += `<td></td></tr>`;

      // 4. Clinical Track: CUPPA Predicted Cancer Type
      html += `<tr class="onco-track-row"><td class="onco-track-label">CUPPA Cancer Type</td>`;
      tumors.forEach(sid => {{
        const sm = sampleMeta[sid] || {{}};
        const isBoundary = patientBoundaries.includes(sid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const cType = sm.cuppa_cancer_type || '';
        const lhood = sm.cuppa_likelihood;
        const cStyle = getCuppaStyle(cType);
        const lhoodPct = (lhood !== null && lhood !== undefined && cType) ? ' (' + (lhood * 100).toFixed(0) + '%)' : '';
        const titleStr = cType ? `${{sid}}: CUPPA ${{cType}}${{lhoodPct}}` : `${{sid}}: CUPPA Not Profiled`;
        html += `<td class="${{bCls}}"><div class="onco-track-cell" style="background:${{cStyle.bg}}; color:${{cStyle.text}}; border:1px solid ${{cStyle.border}}; font-size:0.56rem; font-weight:700;" title="${{titleStr}}">${{cStyle.short}}</div></td>`;
      }});
      html += `<td></td></tr>`;

      // 5. Clinical Track: Purity
      html += `<tr class="onco-track-row"><td class="onco-track-label">Purity</td>`;
      tumors.forEach(sid => {{
        const sm = sampleMeta[sid] || {{}};
        const isBoundary = patientBoundaries.includes(sid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const isPlasma = sm.is_plasma;

        let pur = 0;
        let purDisplay = '--';
        let bg = '#f1f5f9';
        let fg = '#64748b';
        let border = '#e2e8f0';

        if (!isPlasma) {{
          const tObj = (sm.patient_id ? COHORT_DATA.all_tumors.find(x => x.patient_id === sm.patient_id && (x.tumor_id === sm.raw_sample_id || x.tumor_id === sid)) : null) ||
                       COHORT_DATA.all_tumors.find(x => x.tumor_id === sid) || {{}};
          pur = (sm.purity_pct !== undefined && sm.purity_pct !== null) ? sm.purity_pct :
                ((tObj.purity_pct !== undefined && tObj.purity_pct !== null) ? tObj.purity_pct : 0);
          const alpha = Math.min(1.0, Math.max(0.18, pur / 100.0));
          purDisplay = pur > 0 ? (pur >= 99.5 ? '100%' : Math.round(pur) + '%') : '--';
          bg = `rgba(14, 165, 233, ${{alpha}})`;
          fg = (alpha > 0.55) ? '#ffffff' : '#0369a1';
          border = 'rgba(14, 165, 233, 0.4)';
        }} else {{
          const plObj = (sm.patient_id ? (COHORT_DATA.all_plasmas || []).find(x => x.patient_id === sm.patient_id && (x.sample_id === sm.raw_sample_id || x.sample_id === sid)) : null) ||
                        (COHORT_DATA.all_plasmas || []).find(x => x.sample_id === sid) || {{}};
          pur = (sm.wisp_purity_pct !== undefined && sm.wisp_purity_pct !== null) ? sm.wisp_purity_pct :
                ((sm.purple_purity_pct !== undefined && sm.purple_purity_pct !== null) ? sm.purple_purity_pct :
                ((plObj.wisp_purity_pct !== undefined && plObj.wisp_purity_pct !== null) ? plObj.wisp_purity_pct :
                ((plObj.purple_purity_pct !== undefined && plObj.purple_purity_pct !== null) ? plObj.purple_purity_pct : 0)));
          if (pur > 0) {{
            purDisplay = pur >= 10 ? Math.round(pur) + '%' : pur.toFixed(1) + '%';
            const alpha = Math.min(1.0, Math.max(0.25, pur / 20.0));
            bg = `rgba(168, 85, 247, ${{alpha}})`;
            fg = (alpha > 0.55) ? '#ffffff' : '#6b21a8';
            border = 'rgba(168, 85, 247, 0.4)';
          }}
        }}
        html += `<td class="${{bCls}}"><div class="onco-track-cell" style="background:${{bg}}; color:${{fg}}; border:1px solid ${{border}};" title="${{sid}} Purity: ${{pur > 0 ? pur.toFixed(1) + '%' : 'NA'}}">${{purDisplay}}</div></td>`;
      }});
      html += `<td></td></tr>`;

      // 6. Clinical Track: MRD Status
      html += `<tr class="onco-track-row"><td class="onco-track-label">MRD Status</td>`;
      tumors.forEach(sid => {{
        const sm = sampleMeta[sid] || {{}};
        const isBoundary = patientBoundaries.includes(sid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const isPlasma = sm.is_plasma;

        let mrdText = 'Baseline';
        let mrdBg = '#f1f5f9';
        let mrdColor = '#64748b';
        let mrdBorder = '#e2e8f0';

        if (isPlasma) {{
          const plObj = (sm.patient_id ? (COHORT_DATA.all_plasmas || []).find(x => x.patient_id === sm.patient_id && (x.sample_id === sm.raw_sample_id || x.sample_id === sid)) : null) ||
                        (COHORT_DATA.all_plasmas || []).find(x => x.sample_id === sid) || {{}};
          const isPos = sm.is_mrd_positive !== undefined ? sm.is_mrd_positive : plObj.is_mrd_positive;
          const status = sm.mrd_status || plObj.mrd_status;
          if (isPos) {{
            mrdText = 'MRD Pos';
            mrdBg = '#fee2e2';
            mrdColor = '#991b1b';
            mrdBorder = '#fecaca';
          }} else if (status === 'MRD Negative') {{
            mrdText = 'MRD Neg';
            mrdBg = '#dcfce7';
            mrdColor = '#166534';
            mrdBorder = '#bbf7d0';
          }} else {{
            mrdText = '--';
            mrdBg = '#f1f5f9';
            mrdColor = '#94a3b8';
            mrdBorder = '#e2e8f0';
          }}
        }}
        html += `<td class="${{bCls}}"><div class="onco-track-cell" style="background:${{mrdBg}}; color:${{mrdColor}}; border:1px solid ${{mrdBorder}};" title="${{sid}} MRD: ${{mrdText}}">${{mrdText}}</div></td>`;
      }});
      html += `<td></td></tr>`;

      // 7. Clinical Track: TMB Track
      html += `<tr class="onco-track-row"><td class="onco-track-label">TMB (mut/Mb)</td>`;
      tumors.forEach(tid => {{
        const sm = sampleMeta[tid] || {{}};
        const isBoundary = patientBoundaries.includes(tid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const tObj = (sm.patient_id ? COHORT_DATA.all_tumors.find(x => x.patient_id === sm.patient_id && (x.tumor_id === sm.raw_sample_id || x.tumor_id === tid)) : null) ||
                     COHORT_DATA.all_tumors.find(x => x.tumor_id === tid) || {{}};
        const tmb = (sm.tmb !== undefined && sm.tmb !== null) ? sm.tmb :
                    ((tObj.tmb !== undefined && tObj.tmb !== null) ? tObj.tmb : null);
        const tmbStr = tmb !== null ? tmb.toFixed(1) : '--';
        const isHigh = tmb !== null && tmb >= 10;
        const bg = isHigh ? '#fee2e2' : '#f8fafc';
        const fg = isHigh ? '#991b1b' : '#334155';
        const border = isHigh ? '#fecaca' : '#e2e8f0';
        html += `<td class="${{bCls}}"><div class="onco-track-cell" style="background:${{bg}}; color:${{fg}}; border:1px solid ${{border}}; font-weight:${{isHigh ? '700' : '500'}};" title="${{tid}} TMB: ${{tmbStr}} mut/Mb">${{tmbStr}}</div></td>`;
      }});
      html += `<td></td></tr>`;

      // 8. Clinical Track: MSI Track
      html += `<tr class="onco-track-row"><td class="onco-track-label">MSI Status</td>`;
      tumors.forEach(tid => {{
        const sm = sampleMeta[tid] || {{}};
        const isBoundary = patientBoundaries.includes(tid);
        const bCls = isBoundary ? ' patient-boundary' : '';
        const tObj = (sm.patient_id ? COHORT_DATA.all_tumors.find(x => x.patient_id === sm.patient_id && (x.tumor_id === sm.raw_sample_id || x.tumor_id === tid)) : null) ||
                     COHORT_DATA.all_tumors.find(x => x.tumor_id === tid) || {{}};
        const msi = sm.msi_status || tObj.msi_status || (sm.is_plasma ? '--' : 'MSS');
        const isMsiH = msi === 'MSI' || msi === 'MSI-H';
        const bg = isMsiH ? '#fee2e2' : (msi === '--' ? '#f1f5f9' : '#dcfce7');
        const fg = isMsiH ? '#991b1b' : (msi === '--' ? '#94a3b8' : '#166534');
        const border = isMsiH ? '#fecaca' : (msi === '--' ? '#e2e8f0' : '#bbf7d0');
        html += `<td class="${{bCls}}"><div class="onco-track-cell" style="background:${{bg}}; color:${{fg}}; border:1px solid ${{border}}; font-weight:${{isMsiH ? '700' : '500'}};" title="${{tid}} MSI: ${{msi}}">${{isMsiH ? 'MSI' : (msi === '--' ? '--' : 'MSS')}}</div></td>`;
      }});
      html += `<td></td></tr>`;

      // Gene Rows
      genesToShow.forEach(g => {{
        let gInfo = (COHORT_DATA.oncoprint.summary && COHORT_DATA.oncoprint.summary[g]) ? COHORT_DATA.oncoprint.summary[g] : null;
        if (!gInfo) {{
          let alteredCount = 0;
          tumors.forEach(tid => {{
            let cell = (cellDetails[g] && cellDetails[g][tid]) ? cellDetails[g][tid] : getDynamicCellDetails(g, tid);
            if (cell && (cell.cna || (cell.mutations && cell.mutations.length > 0) || cell.fusion || cell.has_sv || cell.structural_variant)) {{
              alteredCount++;
            }}
          }});
          gInfo = {{ total_altered: alteredCount, pct: tumors.length > 0 ? (alteredCount / tumors.length) * 100.0 : 0 }};
        }}
        html += `<tr><td class="gene-name">${{g}}</td>`;

        tumors.forEach(tid => {{
          const isBoundary = patientBoundaries.includes(tid);
          const bCls = isBoundary ? ' patient-boundary' : '';

          let cell = (cellDetails[g] && cellDetails[g][tid]) ? cellDetails[g][tid] : null;
          if (!cell) {{
            cell = getDynamicCellDetails(g, tid);
          }}

          let cellTooltipTitle = `${{g}} | ${{cell.summary}} | Sample: ${{tid}} (Patient: ${{cell.patient_id || 'NA'}})`;
          let cellHtml = `<div class="onco-cell-container" onclick="showOncoCellModal('${{g}}', '${{tid}}')" onmouseenter="showOncoTooltip(event, '${{g}}', '${{tid}}')" onmousemove="moveOncoTooltip(event)" onmouseleave="hideOncoTooltip()" title="${{cellTooltipTitle}}">`;

          // Layer 0: Base / Wild-Type Canvas (cell container background is solid grey)
          if (!cell.is_profiled) {{
            cellHtml += `<div class="onco-layer-not-profiled"></div>`;
          }}

          // Layer 1: CNA (Full-height filled rectangle) & Expression
          if (cell.cna === 'AMP') {{
            cellHtml += `<div class="onco-layer-cna cna-amp"></div>`;
          }} else if (cell.cna === 'HOMDEL') {{
            cellHtml += `<div class="onco-layer-cna cna-homdel"></div>`;
          }} else if (cell.cna === 'GAIN') {{
            cellHtml += `<div class="onco-layer-cna cna-gain"></div>`;
          }} else if (cell.cna === 'HETLOSS') {{
            cellHtml += `<div class="onco-layer-cna cna-hetloss"></div>`;
          }}

          // Layer 1 Modifier: Expression Outline
          if (cell.expression === 'UP') {{
            cellHtml += `<div class="expr-up" style="position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;"></div>`;
          }} else if (cell.expression === 'DOWN') {{
            cellHtml += `<div class="expr-down" style="position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;"></div>`;
          }}

          // Layer 2: Mutation Layer (Centered sub-height rectangle ~37.5% of cell: small sequence variants)
          const smallMuts = (cell.mutations || []).filter(m => m.type !== 'Structural Variant' && m.type !== 'SV');
          if (smallMuts.length > 0) {{
            let mutCls = '';
            let customStyle = '';
            if (smallMuts.length === 1) {{
              const m = smallMuts[0];
              if (m.type === 'Missense') mutCls = 'mut-missense';
              else if (m.type === 'Truncating') mutCls = 'mut-truncating';
              else if (m.type === 'In-frame InDel') mutCls = 'mut-inframe';
              else if (m.type === 'Splice') mutCls = 'mut-splice';
              else if (m.type === 'Promoter') mutCls = 'mut-promoter';
              else mutCls = 'mut-missense';
            }} else {{
              // Concurrent multi-hit mutations: split gradient
              const colorMap = {{
                'Missense': '#008000',
                'Truncating': '#000000',
                'In-frame InDel': '#996633',
                'Splice': '#FF7F00',
                'Promoter': '#008080'
              }};
              const c1 = colorMap[smallMuts[0].type] || '#008000';
              const c2 = colorMap[smallMuts[1].type] || '#000000';
              customStyle = `background: linear-gradient(135deg, ${{c1}} 50%, ${{c2}} 50%);`;
            }}
            cellHtml += `<div class="onco-layer-mut ${{mutCls}}" style="${{customStyle}}"></div>`;
          }}

          // Layer 3: Structural Variant / Fusion Layer (Vertical purple/magenta bar)
          if (cell.fusion || cell.has_sv || cell.structural_variant) {{
            cellHtml += `<div class="onco-layer-fusion" title="Structural Variant / Fusion"></div>`;
          }}

          // Layer 4: Germline / Modifier Layer (Centered horizontal white overlay stripe)
          if (cell.is_germline) {{
            cellHtml += `<div class="onco-layer-germline"></div>`;
          }}

          cellHtml += `</div>`;
          html += `<td class="${{bCls}}">${{cellHtml}}</td>`;
        }});

        html += `<td class="freq-val">${{gInfo.pct.toFixed(0)}}%</td></tr>`;
      }});

      html += `</tbody></table>`;
      container.innerHTML = html;

      // Populate Driver Alterations Table
      populateDriversTable();
    }}

    function showOncoCellModal(gene, tumorId) {{
      const cell = (COHORT_DATA.oncoprint.cell_details && COHORT_DATA.oncoprint.cell_details[gene] && COHORT_DATA.oncoprint.cell_details[gene][tumorId])
        ? COHORT_DATA.oncoprint.cell_details[gene][tumorId]
        : getDynamicCellDetails(gene, tumorId);

      const sm = (COHORT_DATA.oncoprint.sample_metadata && COHORT_DATA.oncoprint.sample_metadata[tumorId]) || {{}};
      const typeLabel = sm.is_plasma ? ('Plasma ctDNA (' + (sm.timepoint || 'Plasma') + ')') : 'Primary Tumor';

      const modal = document.getElementById('onco-cell-modal');
      const content = document.getElementById('onco-modal-content');
      if (!modal || !content) return;

      let cnaInfo = 'None (Diploid / Unaltered)';
      if (cell.cna) {{
        cnaInfo = `<strong>${{cell.cna_label}}</strong> (${{cell.cna}})`;
        if (cell.cna_min_cn !== null && cell.cna_max_cn !== null) {{
          cnaInfo += ` [Copy Number: ${{cell.cna_min_cn.toFixed(2)}} - ${{cell.cna_max_cn.toFixed(2)}}]`;
        }}
      }}

      let mutsHtml = '<p style="color:var(--text-muted);">No sequence variants detected in this sample.</p>';
      if (cell.mutations && cell.mutations.length > 0) {{
        mutsHtml = '<ul style="margin: 0.5rem 0 0 1.2rem; padding: 0;">';
        cell.mutations.forEach(m => {{
          let origin = m.is_germline ? '<span class="badge" style="background:#8b5cf6; color:#fff;">GERMLINE</span>' : '<span class="badge badge-success">SOMATIC</span>';
          if (m.is_status_unknown) {{
            origin = '<span class="badge" style="background:#d97706; color:#fff; font-weight:700;">SOMATIC / GERMLINE STATUS UNKNOWN</span>';
          }}
          const pChange = m.hgvs_protein || 'Protein change unavailable';
          const cChange = m.hgvs_coding ? ('(' + m.hgvs_coding + ')') : '';
          const vafVal = (m.vaf !== undefined && m.vaf !== null) ? (m.vaf > 1.0 ? m.vaf : m.vaf * 100.0) : null;
          let depthStr = '';
          if (m.depth_alt !== undefined && m.depth_alt !== null && m.depth_total) {{
            depthStr = ` (${{m.depth_alt}}/${{m.depth_total}} reads)`;
          }} else if (m.sample_ad !== undefined && m.sample_ad !== null && m.sample_dp) {{
            depthStr = ` (${{m.sample_ad}}/${{m.sample_dp}} reads)`;
          }}
          let vcnStr = '';
          if (m.variant_copy_number !== undefined && m.variant_copy_number !== null) {{
            vcnStr = ` [VCN: ${{m.variant_copy_number.toFixed(1)}}]`;
          }}
          const vafStr = vafVal !== null ? ('VAF: ' + vafVal.toFixed(1) + '%' + depthStr + vcnStr) : '';
          const biallelicBadge = m.biallelic ? ' <span class="badge" style="background:#dc2626; color:#fff; font-size:0.68rem; margin-left:0.25rem;">BIALLELIC</span>' : '';
          let noteStr = '';
          if (m.is_status_unknown) {{
            const ndpStr = (m.normal_dp !== undefined && m.normal_dp !== null) ? m.normal_dp + 'x' : '0x';
            noteStr = `<div style="color:#b45309; font-size:0.78rem; margin-top:0.25rem; font-style:italic;">⚠️ Somatic / germline status unknown due to low normal coverage (${{ndpStr}} &lt; 10x)</div>`;
          }}
          mutsHtml += '<li style="margin-bottom: 0.45rem;"><strong>' + m.type + '</strong> (' + m.subtype + ') ' + origin + biallelicBadge + '<br><span style="color:var(--accent-blue); font-weight:600;">' + pChange + '</span> ' + cChange + ' <span style="color:var(--text-muted); font-size:0.8rem; margin-left:0.35rem;">' + vafStr + '</span>' + noteStr + '</li>';
        }});
        mutsHtml += '</ul>';
      }}

      let germlineCnvBadge = cell.germline_cnv ? ` <span class="badge" style="background:#8b5cf6; color:#fff; margin-left:0.4rem;">${{cell.germline_cnv.toUpperCase()}}</span>` : '';

      let svInfo = '';
      if (cell.structural_variants && cell.structural_variants.length > 0) {{
        svInfo = `
          <div style="margin-bottom:1rem; padding:0.75rem; background:var(--bg-subtle); border-radius:6px; border:1px solid #c084fc;">
            <h4 style="margin:0 0 0.5rem 0; font-size:0.8rem; color:#7c3aed; text-transform:uppercase; display:flex; align-items:center; gap:0.5rem;">
              <span>Structural Variant / Disruption / Fusion</span>
              <span class="badge" style="background:#7c3aed; color:#fff; font-size:0.7rem;">${{cell.structural_variants.length}} EVENT${{cell.structural_variants.length > 1 ? 'S' : ''}}</span>
            </h4>
            <div style="display:flex; flex-direction:column; gap:0.5rem;">
        `;
        cell.structural_variants.forEach(sv => {{
          const typeBadge = `<span class="badge" style="background:#8b5cf6; color:#fff;">${{sv.type || 'SV'}}</span>`;
          const eventLabel = sv.event ? `<span class="badge" style="background:#3b82f6; color:#fff;">${{sv.event}}</span>` : '';
          const originBadge = sv.is_germline ? `<span class="badge" style="background:#ec4899; color:#fff;">GERMLINE</span>` : `<span class="badge badge-success">SOMATIC</span>`;
          const bkStr = sv.breakpoint ? `<div><span style="color:var(--text-muted); font-size:0.8rem;">Breakpoint:</span> <code style="color:var(--accent-blue); font-size:0.85rem;">${{sv.breakpoint}}</code></div>` : '';
          const partnerStr = sv.partner && sv.partner !== 'Intergenic'
            ? `<div><span style="color:var(--text-muted); font-size:0.8rem;">Partner Gene:</span> <strong style="color:#d97706;">${{sv.partner}}</strong> <code style="color:var(--text-muted); font-size:0.8rem;">${{sv.partner_coords ? '(' + sv.partner_coords + ')' : ''}}</code></div>`
            : (sv.partner_coords ? `<div><span style="color:var(--text-muted); font-size:0.8rem;">Partner Breakpoint:</span> <code style="color:var(--text-muted); font-size:0.8rem;">${{sv.partner_coords}}</code> (Intergenic)</div>` : '');
          const exonStr = sv.exons ? `<div><span style="color:var(--text-muted); font-size:0.8rem;">Affected Region:</span> <span style="color:var(--text-primary); font-size:0.85rem;">${{sv.exons}}</span> ${{sv.region_type ? ('<span style="color:var(--text-muted); font-size:0.75rem;">(' + sv.region_type + ')</span>') : ''}}</div>` : '';
          const jcnStr = (sv.jcn !== null && sv.jcn !== undefined) ? `<div><span style="color:var(--text-muted); font-size:0.8rem;">Junction Copy Number (JCN):</span> <strong style="color:#059669;">${{typeof sv.jcn === 'number' ? sv.jcn.toFixed(2) : sv.jcn}}</strong></div>` : '';
          
          svInfo += `
            <div style="background:#ffffff; padding:0.5rem 0.65rem; border-radius:4px; border:1px solid var(--border-color); font-size:0.85rem;">
              <div style="display:flex; align-items:center; gap:0.4rem; margin-bottom:0.35rem; flex-wrap:wrap;">
                ${{typeBadge}} ${{eventLabel}} ${{originBadge}} <strong style="color:var(--text-primary);">${{sv.name || sv.gene}}</strong>
              </div>
              <div style="display:grid; grid-template-columns: 1fr; gap:0.2rem; margin-left:0.25rem;">
                ${{bkStr}}
                ${{partnerStr}}
                ${{exonStr}}
                ${{jcnStr}}
              </div>
            </div>
          `;
        }});
        svInfo += `
            </div>
          </div>
        `;
      }} else if (cell.fusion || cell.has_sv) {{
        const svName = cell.sv_name || cell.fusion_name || 'Driver Structural Variant / Breakend';
        svInfo = `
          <div style="margin-bottom:1rem; padding:0.75rem; background:var(--bg-subtle); border-radius:6px; border:1px solid #c084fc;">
            <h4 style="margin:0 0 0.35rem 0; font-size:0.8rem; color:#7c3aed; text-transform:uppercase;">Structural Variant / Fusion</h4>
            <div style="color:var(--text-primary); font-weight:600; font-size:0.9rem;">${{svName}}</div>
            <div style="color:var(--text-muted); font-size:0.78rem; margin-top:0.25rem;">Detected via LINX / Orange Report driver breakend annotations.</div>
          </div>
        `;
      }}
      let exprInfo = cell.expression ? `<p style="color:#0284c7;"><strong>Expression:</strong> mRNA ${{cell.expression}}</p>` : '';

      let cuppaInfo = '';
      if (sm.cuppa_cancer_type) {{
        const lhood = sm.cuppa_likelihood;
        const lhoodStr = (lhood !== null && lhood !== undefined) ? ` (${{(lhood * 100).toFixed(1)}}% confidence)` : '';
        const cStyle = getCuppaStyle(sm.cuppa_cancer_type);
        let predsHtml = '';
        if (sm.cuppa_predictions && sm.cuppa_predictions.length > 0) {{
          predsHtml = '<div style="margin-top:0.4rem; font-size:0.75rem; display:flex; flex-wrap:wrap; gap:4px;">';
          sm.cuppa_predictions.slice(0, 4).forEach(p => {{
            predsHtml += `<span class="badge" style="background:#ffffff; color:var(--text-primary); border:1px solid var(--border-color);">${{p.cancer_type}}: ${{(p.likelihood * 100).toFixed(1)}}%</span>`;
          }});
          predsHtml += '</div>';
        }}
        cuppaInfo = `
          <div style="margin-bottom:1rem; padding:0.75rem; background:var(--bg-subtle); border-radius:6px; border:1px solid var(--border-color);">
            <h4 style="margin:0 0 0.35rem 0; font-size:0.8rem; color:var(--text-secondary); text-transform:uppercase;">CUPPA Cancer of Unknown Primary Prediction</h4>
            <div style="display:flex; align-items:center; gap:8px;">
              <span class="badge" style="background:${{cStyle.bg}}; color:${{cStyle.text}}; border:1px solid ${{cStyle.border}}; font-weight:700; font-size:0.85rem; padding:3px 8px;">${{sm.cuppa_cancer_type}}</span>
              <span style="font-size:0.82rem; color:var(--text-muted);">${{lhoodStr}}</span>
            </div>
            ${{predsHtml}}
          </div>
        `;
      }}

      content.innerHTML = `
        <h3 style="margin-top:0; color:var(--accent-blue); font-size:1.15rem;">${{gene}} in ${{tumorId}}</h3>
        <p style="color:var(--text-secondary); font-size:0.85rem; margin-bottom:1rem;">Patient ID: <strong>${{cell.patient_id}}</strong> | Sample: <strong>${{tumorId}}</strong> (${{typeLabel}})</p>
        ${{cuppaInfo}}
        <div style="margin-bottom:1rem; padding:0.75rem; background:var(--bg-subtle); border-radius:6px; border:1px solid var(--border-color);">
          <h4 style="margin:0 0 0.35rem 0; font-size:0.8rem; color:var(--text-secondary); text-transform:uppercase;">Copy Number Alteration (CNA)</h4>
          <div>${{cnaInfo}}${{germlineCnvBadge}}</div>
        </div>
        <div style="margin-bottom:1rem; padding:0.75rem; background:var(--bg-subtle); border-radius:6px; border:1px solid var(--border-color);">
          <h4 style="margin:0 0 0.35rem 0; font-size:0.8rem; color:var(--text-secondary); text-transform:uppercase;">Mutations (Small Variants)</h4>
          ${{mutsHtml}}
        </div>
        ${{svInfo}}
        ${{exprInfo}}
      `;

      modal.style.display = 'flex';
    }}

    function closeOncoCellModal() {{
      const modal = document.getElementById('onco-cell-modal');
      if (modal) modal.style.display = 'none';
    }}

    function showOncoTooltip(e, gene, tumorId) {{
      const tip = document.getElementById('oncoprint-floating-tooltip');
      if (!tip) return;
      const cell = (COHORT_DATA.oncoprint.cell_details && COHORT_DATA.oncoprint.cell_details[gene] && COHORT_DATA.oncoprint.cell_details[gene][tumorId])
        ? COHORT_DATA.oncoprint.cell_details[gene][tumorId]
        : getDynamicCellDetails(gene, tumorId);

      const sm = (COHORT_DATA.oncoprint.sample_metadata && COHORT_DATA.oncoprint.sample_metadata[tumorId]) || {{}};
      const typeLabel = sm.is_plasma ? ('Plasma ctDNA (' + (sm.timepoint || 'Plasma') + ')') : 'Primary Tumor';
      const patientId = cell.patient_id || sm.patient_id || 'Unknown';

      let contentHtml = `
        <div class="onco-tooltip-gene">
          <span>${{gene}}</span>
          <span style="font-size:0.75rem; font-weight:normal; color:var(--text-muted);">${{cell.cna_label || (cell.mutations && cell.mutations.length > 0 ? cell.mutations[0].type : 'Wild Type')}}</span>
        </div>
        <div class="onco-tooltip-sample">
          <strong>${{patientId}}</strong> · ${{tumorId}} <span class="onco-tooltip-tag">${{typeLabel}}</span>
        </div>
      `;

      if (sm.cuppa_cancer_type) {{
        const lhood = sm.cuppa_likelihood;
        const lhoodStr = (lhood !== null && lhood !== undefined) ? ` (${{(lhood * 100).toFixed(0)}}%)` : '';
        const cStyle = getCuppaStyle(sm.cuppa_cancer_type);
        contentHtml += `
          <div style="font-size:0.74rem; margin-bottom:5px; display:flex; align-items:center; gap:4px;">
            <span style="color:var(--text-secondary); font-weight:600;">CUPPA:</span>
            <span class="badge" style="background:${{cStyle.bg}}; color:${{cStyle.text}}; border:1px solid ${{cStyle.border}}; font-size:0.68rem; padding:1px 5px;">${{sm.cuppa_cancer_type}}${{lhoodStr}}</span>
          </div>
        `;
      }}

      if (cell.mutations && cell.mutations.length > 0) {{
        contentHtml += `<div style="font-size:0.72rem; text-transform:uppercase; color:#8b949e; margin-top:4px; font-weight:600;">Sequence Alteration(s):</div>`;
        cell.mutations.forEach(m => {{
          const isGermline = m.is_germline;
          const isTrunc = (m.type === 'Truncating');
          const itemCls = isGermline ? 'germline' : (isTrunc ? 'truncating' : '');
          const aa = m.hgvs_protein || m.hgvs_coding || m.type;
          const coding = m.hgvs_coding && m.hgvs_protein ? ('(' + m.hgvs_coding + ')') : '';
          const vafVal = (m.vaf !== undefined && m.vaf !== null) ? (m.vaf > 1.0 ? m.vaf : m.vaf * 100.0) : null;
          let depthStr = '';
          if (m.depth_alt !== undefined && m.depth_alt !== null && m.depth_total) {{
            depthStr = ` (${{m.depth_alt}}/${{m.depth_total}})`;
          }} else if (m.sample_ad !== undefined && m.sample_ad !== null && m.sample_dp) {{
            depthStr = ` (${{m.sample_ad}}/${{m.sample_dp}})`;
          }}
          const vafStr = vafVal !== null ? (' · VAF ' + vafVal.toFixed(1) + '%' + depthStr) : '';
          let originBadge = isGermline ? '<span class="onco-tooltip-tag" style="background:#8b5cf6; color:#fff;">Germline</span>' : '<span class="onco-tooltip-tag" style="background:#238636; color:#fff;">Somatic</span>';
          if (m.is_status_unknown) {{
            originBadge = '<span class="onco-tooltip-tag" style="background:#d97706; color:#fff;">Status Unknown</span>';
          }}
          
          contentHtml += `
            <div class="onco-tooltip-mut-item ${{itemCls}}">
              <div style="display:flex; justify-content:space-between; align-items:center;">
                <span class="onco-tooltip-aa">${{aa}}</span>
                ${{originBadge}}
              </div>
              <div style="font-size:0.74rem; color:#8b949e; margin-top:1px;">
                ${{m.type}} ${{coding}}${{vafStr}}
              </div>
            </div>
          `;
        }});
      }}

      if (cell.cna) {{
        let cnRange = '';
        if (cell.cna_min_cn !== null && cell.cna_max_cn !== null) {{
          cnRange = ` (CN: ${{cell.cna_min_cn.toFixed(1)}}-${{cell.cna_max_cn.toFixed(1)}})`;
        }}
        contentHtml += `
          <div style="margin-top:5px; padding:3px 6px; background:#f1f5f9; border-radius:4px; border:1px solid #e2e8f0;">
            <strong style="color:#d97706;">CNA:</strong> ${{cell.cna_label}} (${{cell.cna}})${{cnRange}}
          </div>
        `;
      }}

      if (cell.fusion || cell.has_sv || (cell.structural_variants && cell.structural_variants.length > 0)) {{
        const svName = cell.sv_name || cell.fusion_name || (cell.structural_variants && cell.structural_variants[0] ? cell.structural_variants[0].name : 'Driver Breakend');
        contentHtml += `
          <div style="margin-top:5px; padding:3px 6px; background:#f1f5f9; border-radius:4px; border-left:3px solid #8b5cf6; border-top:1px solid #e2e8f0; border-right:1px solid #e2e8f0; border-bottom:1px solid #e2e8f0;">
            <strong style="color:#7c3aed;">SV / Fusion:</strong> ${{svName}}
          </div>
        `;
      }}

      if (!cell.cna && (!cell.mutations || cell.mutations.length === 0) && !cell.fusion && !cell.has_sv) {{
        contentHtml += `<div style="color:#64748b; font-style:italic; margin-top:2px;">Wild Type (Diploid / No Drivers Detected)</div>`;
      }}

      contentHtml += `<div style="font-size:0.68rem; color:#94a3b8; margin-top:6px; text-align:right;">Click cell to inspect full alteration profile</div>`;

      tip.innerHTML = contentHtml;
      tip.style.display = 'block';
      tip.classList.add('visible');
      positionOncoTooltip(e);
    }}

    function moveOncoTooltip(e) {{
      positionOncoTooltip(e);
    }}

    function positionOncoTooltip(e) {{
      const tip = document.getElementById('oncoprint-floating-tooltip');
      if (!tip || tip.style.display === 'none') return;
      const pad = 14;
      let left = e.clientX + pad;
      let top = e.clientY + pad;

      const rect = tip.getBoundingClientRect();
      if (left + rect.width > window.innerWidth - 10) {{
        left = e.clientX - rect.width - pad;
      }}
      if (top + rect.height > window.innerHeight - 10) {{
        top = e.clientY - rect.height - pad;
      }}
      tip.style.left = Math.max(8, left) + 'px';
      tip.style.top = Math.max(8, top) + 'px';
    }}

    function hideOncoTooltip() {{
      const tip = document.getElementById('oncoprint-floating-tooltip');
      if (!tip) return;
      tip.classList.remove('visible');
      tip.style.display = 'none';
    }}

    function showPlasmaModalityModal(sampleId, patientId) {{
      let pl = null;
      if (patientId) {{
        pl = (COHORT_DATA.all_plasmas || []).find(p => p.patient_id === patientId && (p.sample_id === sampleId || p.clean_sample_name === sampleId));
      }}
      if (!pl) {{
        pl = (COHORT_DATA.all_plasmas || []).find(p => p.sample_id === sampleId || p.clean_sample_name === sampleId);
      }}
      if (!pl) return;

      const modal = document.getElementById('plasma-modality-modal');
      const headerInfo = document.getElementById('plasma-modal-header-info');
      const downloadBtn = document.getElementById('plasma-modal-download-btn');
      const body = document.getElementById('plasma-modal-body');
      if (!modal || !headerInfo || !body) return;

      const statusBadge = pl.is_mrd_positive 
        ? '<span class="badge badge-mrd-pos" style="font-size:0.85rem; padding:4px 10px;">MRD POSITIVE</span>'
        : '<span class="badge badge-mrd-neg" style="font-size:0.85rem; padding:4px 10px;">MRD NEGATIVE</span>';

      const sName = pl.clean_sample_name || pl.sample_id;
      headerInfo.innerHTML = `
        <div style="display:flex; align-items:center; gap:12px; flex-wrap:wrap;">
          <h3 style="margin:0; font-size:1.15rem; color:#0f172a; font-weight:700;">
            <span style="color:#0284c7;">Patient: ${{pl.patient_id}}</span> &nbsp;|&nbsp; 
            <span style="color:#0f172a;">Sample: ${{sName}}</span>
          </h3>
          <span class="badge badge-primary" style="font-size:0.8rem;">Timepoint: ${{pl.timepoint || '--'}}</span>
          <span style="color:#64748b; font-size:0.85rem;">Tumor Baseline: <strong style="color:#0f172a;">${{pl.tumor_id || 'Tumor'}}</strong></span>
          ${{statusBadge}}
        </div>
      `;

      const plotUrl = pl.multimodal_plot_b64 || pl.multimodal_plot || `plasma_plots/${{pl.patient_id}}_${{pl.sample_id}}_wisp_multimodal.png`;
      downloadBtn.href = plotUrl;
      downloadBtn.download = `${{pl.patient_id}}_${{pl.sample_id}}_wisp_multimodal.png`;

      const snvCall = pl.snv_mrd ? 'POSITIVE' : 'Negative (&lt; LOD)';
      const snvBadgeClass = pl.snv_mrd ? 'badge-mrd-pos' : 'badge-mrd-neg';
      const cnvCall = pl.cnv_mrd ? 'POSITIVE' : 'Negative (&lt; LOD)';
      const cnvBadgeClass = pl.cnv_mrd ? 'badge-mrd-pos' : 'badge-mrd-neg';
      const lohCall = pl.loh_mrd ? 'POSITIVE' : 'Negative (&lt; LOD)';
      const lohBadgeClass = pl.loh_mrd ? 'badge-mrd-pos' : 'badge-mrd-neg';

      const hasIchor = Boolean(
        (COHORT_DATA.ichor_reports_b64 && (COHORT_DATA.ichor_reports_b64[pl.sample_id] || (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[pl.clean_sample_name]))) ||
        pl.ichor_pdf_b64 ||
        pl.ichor_pdf_rel ||
        pl.has_ichor_pdf
      );
      const ichorCardBtn = hasIchor
        ? `<button type="button" onclick="showIchorPdfModal('${{pl.sample_id}}', '${{pl.patient_id}}')" class="btn btn-sm" style="width:100%; font-size:0.75rem; padding:4px 10px; border-radius:4px; border:1px solid #f59e0b; background:#fef3c7; color:#92400e; cursor:pointer; font-weight:600; display:inline-flex; align-items:center; justify-content:center; gap:4px;" title="View native ichorCNA Genome-Wide Solution Plot in browser">📊 View Genome-Wide Plot (PDF)</button>`
        : '';

      body.innerHTML = `
        <!-- Main Plot Image -->
        <div style="text-align:center; margin-bottom:1.5rem; background:#f8fafc; padding:0.75rem; border-radius:8px; border:1px solid #e2e8f0;">
          <img src="${{plotUrl}}" alt="Multi-Modal WISP Decomposition" 
               style="width:100%; max-height:68vh; object-fit:contain; border-radius:4px; display:block; margin:0 auto; cursor:zoom-in;"
               onclick="window.open('${{plotUrl}}', '_blank')"
               onerror="this.onerror=null; this.parentElement.innerHTML='<div style=\\'padding:3rem; color:#64748b;\\'>Plot image not found at ${{plotUrl}}. High-resolution figures are generated when running compile_wgts_summary.py.</div>';" />
          <div style="font-size:0.75rem; color:#64748b; margin-top:6px; text-align:right;">
            <em>Click image to view in new tab or download using the button above</em>
          </div>
        </div>

        <!-- 4-Column Modality Breakdown Cards -->
        <div style="display:grid; grid-template-columns:repeat(auto-fit, minmax(280px, 1fr)); gap:1rem;">
          <!-- Card 1: CNV -->
          <div style="background:#ffffff; border:1px solid #e2e8f0; border-top:3px solid #0284c7; border-radius:6px; padding:1rem; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.5rem;">
              <h4 style="margin:0; font-size:0.95rem; color:#0284c7; font-weight:700;">1. Copy Number (CNV) Modality</h4>
              <span class="badge ${{cnvBadgeClass}}">${{cnvCall}}</span>
            </div>
            <p style="font-size:0.78rem; color:#64748b; margin:0 0 0.75rem 0; line-height:1.4;">
              Compares GC-corrected plasma depth ratio against tumor segment copy numbers using linear regression.
            </p>
            <div style="font-size:0.82rem; font-family:monospace; background:#f8fafc; padding:0.6rem; border-radius:4px; border:1px solid #e2e8f0; color:#0f172a;">
              <div>Purity:    <strong style="color:#0284c7;">${{(pl.cnv_purity_pct || 0).toFixed(4)}}%</strong></div>
              <div>Assay LOD: <span style="color:#7c3aed;">${{(pl.cnv_lod_pct || 0).toFixed(4)}}%</span></div>
              <div>p-value:   <span style="color:#0f172a;">${{pl.cnv_pvalue ? pl.cnv_pvalue.toExponential(2) : 'NA'}}</span></div>
            </div>
          </div>

          <!-- Card 2: LOH -->
          <div style="background:#ffffff; border:1px solid #e2e8f0; border-top:3px solid #8b5cf6; border-radius:6px; padding:1rem; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.5rem;">
              <h4 style="margin:0; font-size:0.95rem; color:#7c3aed; font-weight:700;">2. Loss of Heterozygosity (LOH)</h4>
              <span class="badge ${{lohBadgeClass}}">${{lohCall}}</span>
            </div>
            <p style="font-size:0.78rem; color:#64748b; margin:0 0 0.75rem 0; line-height:1.4;">
              Aggregates phased germline SNPs across tumor LOH segments to measure deviation from the 50/50 baseline.
            </p>
            <div style="font-size:0.82rem; font-family:monospace; background:#f8fafc; padding:0.6rem; border-radius:4px; border:1px solid #e2e8f0; color:#0f172a;">
              <div>Purity:    <strong style="color:#7c3aed;">${{(pl.loh_purity_pct || 0).toFixed(4)}}%</strong></div>
              <div>Assay LOD: <span style="color:#d97706;">${{(pl.loh_lod_pct || 0).toFixed(4)}}%</span></div>
              <div>p-value:   <span style="color:#0f172a;">${{pl.loh_pvalue ? pl.loh_pvalue.toExponential(2) : 'NA'}}</span></div>
            </div>
          </div>

          <!-- Card 3: SNV -->
          <div style="background:#ffffff; border:1px solid #e2e8f0; border-top:3px solid #16a34a; border-radius:6px; padding:1rem; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.5rem;">
              <h4 style="margin:0; font-size:0.95rem; color:#16a34a; font-weight:700;">3. Somatic Mutations (SNV)</h4>
              <span class="badge ${{snvBadgeClass}}">${{snvCall}}</span>
            </div>
            <p style="font-size:0.78rem; color:#64748b; margin:0 0 0.75rem 0; line-height:1.4;">
              Tracks tumor somatic SNVs & indels in plasma cfDNA against Hartwig BQR background sequencing error.
            </p>
            <div style="font-size:0.82rem; font-family:monospace; background:#f8fafc; padding:0.6rem; border-radius:4px; border:1px solid #e2e8f0; color:#0f172a;">
              <div>Purity:    <strong style="color:#16a34a;">${{(pl.snv_purity_pct || 0).toFixed(4)}}%</strong></div>
              <div>Assay LOD: <span style="color:#0284c7;">${{(pl.snv_lod_pct || 0).toFixed(4)}}%</span></div>
              <div>p-value:   <span style="color:#0f172a;">${{pl.snv_pvalue ? pl.snv_pvalue.toExponential(2) : 'NA'}}</span></div>
            </div>
          </div>

          <!-- Card 4: 1X ichorCNA (Short Fragments) -->
          <div style="background:#ffffff; border:1px solid #e2e8f0; border-top:3px solid #d97706; border-radius:6px; padding:1rem; box-shadow:0 1px 3px rgba(0,0,0,0.05); display:flex; flex-direction:column; justify-content:space-between;">
            <div>
              <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.5rem;">
                <h4 style="margin:0; font-size:0.95rem; color:#d97706; font-weight:700;">4. 1X ichorCNA (Short Fragments)</h4>
                <span class="badge ${{pl.ichor_tf !== null && pl.ichor_tf >= 0.03 ? 'badge-mrd-pos' : 'badge-mrd-neg'}}">${{pl.ichor_tf !== null ? (pl.ichor_tf >= 0.03 ? 'POSITIVE' : 'Negative (&lt; 3%)') : 'NA'}}</span>
              </div>
              <p style="font-size:0.78rem; color:#64748b; margin:0 0 0.75rem 0; line-height:1.4;">
                Downsampled ~1.0X whole-genome profile with in-silico 90–150 bp size selection and 1 Mb HMM binning.
              </p>
              <div style="font-size:0.82rem; font-family:monospace; background:#f8fafc; padding:0.6rem; border-radius:4px; border:1px solid #e2e8f0; color:#0f172a; margin-bottom:0.75rem;">
                <div>Tumor Fraction: <strong style="color:#d97706;">${{(pl.ichor_tf !== null && pl.ichor_tf !== undefined) ? (pl.ichor_tf * 100).toFixed(2) + '%' : 'NA'}}</strong></div>
                <div>Assay LoD Floor: <span style="color:#94a3b8;">3.00%</span></div>
                <div>Companion Plot:  <span style="color:#0f172a;">${{hasIchor ? 'Available' : 'None'}}</span></div>
              </div>
            </div>
            ${{ichorCardBtn}}
          </div>
        </div>
      `;

      modal.style.display = 'flex';
    }}

    function closePlasmaModalityModal() {{
      const modal = document.getElementById('plasma-modality-modal');
      if (modal) modal.style.display = 'none';
    }}

    function showPatientTrajectoryModal(patientId) {{
      const prec = COHORT_DATA.patients[patientId];
      if (!prec) return;

      const modal = document.getElementById('patient-trajectory-modal');
      const headerInfo = document.getElementById('patient-trajectory-header-info');
      const downloadPng = document.getElementById('patient-trajectory-download-png');
      const downloadPdf = document.getElementById('patient-trajectory-download-pdf');
      const body = document.getElementById('patient-trajectory-body');
      if (!modal || !headerInfo || !body) return;

      const plasmas = prec.plasmas || [];
      const tumors = prec.tumors || [];
      const latestPl = plasmas.length > 0 ? plasmas[plasmas.length - 1] : null;
      const latestMrd = latestPl ? latestPl.mrd_status : 'Unknown';
      const isPos = latestPl ? latestPl.is_mrd_positive : false;
      const statusBadge = isPos 
        ? '<span class="badge badge-mrd-pos" style="font-size:0.85rem; padding:4px 10px;">LATEST: MRD POSITIVE</span>'
        : '<span class="badge badge-mrd-neg" style="font-size:0.85rem; padding:4px 10px;">LATEST: MRD NEGATIVE</span>';

      const baseTumor = tumors.length > 0 ? tumors[0] : null;
      const tumorStr = baseTumor ? `${{baseTumor.tumor_id}} (${{baseTumor.purity_pct.toFixed(1)}}% purity)` : 'None';
      const hasOrangeBase = baseTumor && Boolean(
        (COHORT_DATA.orange_reports_b64 && (COHORT_DATA.orange_reports_b64[baseTumor.tumor_id] || COHORT_DATA.orange_reports_b64[baseTumor.raw_tumor_id])) ||
        baseTumor.orange_pdf_b64 ||
        baseTumor.orange_pdf_rel ||
        baseTumor.has_orange_pdf
      );
      const baseOrangeBtn = hasOrangeBase
        ? `<button type="button" onclick="showOrangePdfModal('${{baseTumor.tumor_id}}', '${{patientId}}')" class="btn btn-sm" style="font-size:0.72rem; padding:2px 8px; border-radius:4px; border:1px solid #f59e0b; background:#fef3c7; color:#92400e; cursor:pointer; font-weight:600; display:inline-flex; align-items:center; gap:3px; margin-left:6px;" title="View native HMF Orange Clinical PDF Report">📄 Orange Report</button>`
        : '';

      const isTracks = (currentMrdViewMode === 'TRACKS');
      headerInfo.innerHTML = `
        <div style="display:flex; align-items:center; gap:12px; flex-wrap:wrap;">
          <h3 style="margin:0; font-size:1.15rem; color:#0f172a; font-weight:700;">${{patientId}} — Longitudinal Plasma MRD Trajectory ${{isTracks ? '(Multi-Modal Breakdown)' : '(Clinical Standard)'}}</h3>
          <span class="badge badge-primary" style="font-size:0.8rem;">${{plasmas.length}} Plasma Timepoints</span>
          <span style="color:#64748b; font-size:0.85rem;">Baseline Tumor: <strong style="color:#0f172a;">${{tumorStr}}</strong>${{baseOrangeBtn}}</span>
          ${{statusBadge}}
        </div>
      `;

      const plotUrl = isTracks
        ? (prec.multimodal_plot_b64 || prec.multimodal_plot || `patient_plots/${{patientId}}_mrd_multimodal.png`)
        : (prec.trajectory_plot_b64 || prec.trajectory_plot || `patient_plots/${{patientId}}_mrd_trajectory.png`);
      downloadPng.href = plotUrl;
      downloadPng.download = isTracks ? `${{patientId}}_mrd_multimodal.png` : `${{patientId}}_mrd_trajectory.png`;

      const pdfUrl = isTracks ? `patient_plots/${{patientId}}_mrd_multimodal.pdf` : `patient_plots/${{patientId}}_mrd_trajectory.pdf`;
      downloadPdf.href = pdfUrl;
      downloadPdf.download = isTracks ? `${{patientId}}_mrd_multimodal.pdf` : `${{patientId}}_mrd_trajectory.pdf`;

      let timepointRows = '';
      plasmas.forEach(pl => {{
        const badge = pl.is_mrd_positive ? '<span class="badge badge-mrd-pos">Positive</span>' : '<span class="badge badge-mrd-neg">Negative</span>';
        const snvStr = pl.snv_purity_pct > 0 ? `<strong>${{pl.snv_purity_pct.toFixed(4)}}%</strong>` : '<span style="color:#64748b;">0.00% (&lt; LOD)</span>';
        const cnvStr = pl.cnv_purity_pct > 0 ? `${{pl.cnv_purity_pct.toFixed(4)}}%` : '<span style="color:#64748b;">0.00% (&lt; LOD)</span>';
        const lohStr = pl.loh_purity_pct > 0 ? `${{pl.loh_purity_pct.toFixed(4)}}%` : '<span style="color:#64748b;">0.00% (&lt; LOD)</span>';
        const wispStr = pl.wisp_purity_pct > 0 ? `<strong>${{pl.wisp_purity_pct.toFixed(4)}}%</strong>` : '<span style="color:#64748b;">0.00% (&lt; LOD)</span>';
        const subDecompBtn = `<button class="btn btn-sm" onclick="showPlasmaModalityModal('${{pl.sample_id}}', '${{patientId}}')" style="font-size:0.75rem; padding:2px 8px; border-radius:4px; border:1px solid #cbd5e1; background:#ffffff; color:#0284c7; cursor:pointer;">📊 Decomp</button>`;

        const hasIchor = Boolean(
          (COHORT_DATA.ichor_reports_b64 && (COHORT_DATA.ichor_reports_b64[pl.sample_id] || (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[pl.clean_sample_name]))) ||
          pl.ichor_pdf_b64 ||
          pl.ichor_pdf_rel ||
          pl.has_ichor_pdf
        );
        const ichorBtn = hasIchor
          ? `<button type="button" onclick="showIchorPdfModal('${{pl.sample_id}}', '${{patientId}}')" class="badge" style="background:#fef3c7; color:#b45309; border:1px solid #f59e0b; cursor:pointer; margin-left:4px; font-size:0.68rem; padding:2px 6px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:2px;" title="View 1X ichorCNA Genome-Wide Solution Plot in browser">📊 Plot</button>`
          : '';
        const ichorStr = (pl.ichor_tf !== null && pl.ichor_tf !== undefined) ? `<strong style="color:#d97706;">${{(pl.ichor_tf * 100).toFixed(2)}}%</strong>${{ichorBtn}}` : '<span style="color:#64748b;">--</span>';

        timepointRows += `
          <tr>
            <td><span class="badge badge-primary">${{pl.timepoint}}</span></td>
            <td><code>${{pl.sample_id}}</code></td>
            <td>${{badge}}</td>
            <td style="color:#0284c7;">${{snvStr}}</td>
            <td style="color:#64748b; font-size:0.8rem;">${{pl.snv_lod_pct ? pl.snv_lod_pct.toFixed(4) + '%' : '--'}}</td>
            <td style="color:#7c3aed;">${{cnvStr}}</td>
            <td style="color:#64748b; font-size:0.8rem;">${{pl.cnv_lod_pct ? pl.cnv_lod_pct.toFixed(4) + '%' : '--'}}</td>
            <td style="color:#d97706;">${{lohStr}}</td>
            <td style="color:#64748b; font-size:0.8rem;">${{pl.loh_lod_pct ? pl.loh_lod_pct.toFixed(4) + '%' : '--'}}</td>
            <td style="color:#dc2626;">${{wispStr}}</td>
            <td style="color:#d97706;">${{ichorStr}}</td>
            <td>${{pl.weighted_depth ? pl.weighted_depth.toFixed(1) + 'x' : '--'}}</td>
            <td>${{subDecompBtn}}</td>
          </tr>
        `;
      }});

      body.innerHTML = `
        <div style="text-align:center; margin-bottom:1.5rem; background:#f8fafc; padding:0.75rem; border-radius:8px; border:1px solid #e2e8f0;">
          <img src="${{plotUrl}}" alt="MRD Trajectory: ${{patientId}}" 
               style="width:100%; max-height:68vh; object-fit:contain; border-radius:4px; display:block; margin:0 auto; cursor:zoom-in;"
               onclick="window.open('${{plotUrl}}', '_blank')"
               onerror="this.onerror=null; this.parentElement.innerHTML='<div style=\\'padding:3rem; color:#64748b;\\'>Trajectory plot image not found at ${{plotUrl}}. High-resolution figures are generated when running compile_wgts_summary.py with matplotlib.</div>';" />
          <div style="font-size:0.75rem; color:#64748b; margin-top:6px; text-align:right;">
            <em>Click image to view in new tab or download using the buttons above</em>
          </div>
        </div>

        <div style="background:#ffffff; border:1px solid #e2e8f0; border-radius:8px; padding:1rem; box-shadow:0 1px 3px rgba(0,0,0,0.05);">
          <h4 style="margin:0 0 0.75rem 0; font-size:0.95rem; color:#0f172a; font-weight:700;">Longitudinal Timepoint Summary</h4>
          <div class="table-container" style="max-height:260px; overflow-y:auto;">
            <table>
              <thead>
                <tr>
                  <th>Timepoint</th>
                  <th>Sample ID</th>
                  <th>MRD Status</th>
                  <th>SNV ctDNA</th>
                  <th>SNV LOD</th>
                  <th>CNV ctDNA</th>
                  <th>CNV LOD</th>
                  <th>LOH ctDNA</th>
                  <th>LOH LOD</th>
                  <th>WISP ctDNA</th>
                  <th>1X ichorCNA TF</th>
                  <th>Depth</th>
                  <th>Decomposition</th>
                </tr>
              </thead>
              <tbody>
                ${{timepointRows}}
              </tbody>
            </table>
          </div>
        </div>
      `;

      modal.style.display = 'flex';
    }}

    function closePatientTrajectoryModal() {{
      const modal = document.getElementById('patient-trajectory-modal');
      if (modal) modal.style.display = 'none';
    }}

    // Cache for decoded PDF Blob URLs and raw bytes
    const orangeBlobUrlCache = {{}};
    const orangeBytesCache = {{}};
    let currentPdfDoc = null;
    let currentPdfScale = 1.35;

    function b64ToUint8Array(b64Data) {{
      const cleanB64 = b64Data.indexOf(',') !== -1 ? b64Data.split(',')[1] : b64Data;
      const byteChars = atob(cleanB64.trim());
      const bytes = new Uint8Array(byteChars.length);
      for (let i = 0; i < byteChars.length; i++) {{
        bytes[i] = byteChars.charCodeAt(i);
      }}
      return bytes;
    }}

    function getOrangePdfBase64(tumorId) {{
      let b64 = (COHORT_DATA.orange_reports_b64 && COHORT_DATA.orange_reports_b64[tumorId]) || null;
      if (!b64) {{
        const t = (COHORT_DATA.all_tumors || []).find(x => x.tumor_id === tumorId || x.raw_tumor_id === tumorId);
        if (t) {{
          b64 = (COHORT_DATA.orange_reports_b64 && (
            COHORT_DATA.orange_reports_b64[t.tumor_id] || 
            (t.raw_tumor_id && COHORT_DATA.orange_reports_b64[t.raw_tumor_id])
          )) || t.orange_pdf_b64;
        }}
      }}
      return b64;
    }}

    function getOrangePdfBlobUrl(tumorId) {{
      if (orangeBlobUrlCache[tumorId]) {{
        return orangeBlobUrlCache[tumorId];
      }}

      const b64 = getOrangePdfBase64(tumorId);
      if (b64) {{
        try {{
          const bytes = orangeBytesCache[tumorId] || b64ToUint8Array(b64);
          orangeBytesCache[tumorId] = bytes;
          const blob = new Blob([bytes], {{ type: 'application/pdf' }});
          const blobUrl = URL.createObjectURL(blob);
          orangeBlobUrlCache[tumorId] = blobUrl;
          return blobUrl;
        }} catch (err) {{
          console.warn('Blob URL creation not supported, falling back to data URL:', err);
          return 'data:application/pdf;base64,' + b64;
        }}
      }}

      // Fallback to relative URL if base64 not available
      const t = (COHORT_DATA.all_tumors || []).find(x => x.tumor_id === tumorId || x.raw_tumor_id === tumorId);
      if (t && t.orange_pdf_rel) {{
        return t.orange_pdf_rel;
      }}

      return null;
    }}

    async function renderOrangePdfPages(pdfDoc, scale) {{
      const container = document.getElementById('orange-pdf-canvas-container');
      if (!container || !pdfDoc) return;
      container.innerHTML = '';

      const outputScale = window.devicePixelRatio || 1;

      for (let pageNum = 1; pageNum <= pdfDoc.numPages; pageNum++) {{
        const page = await pdfDoc.getPage(pageNum);
        const viewport = page.getViewport({{ scale: scale }});

        const pageWrapper = document.createElement('div');
        pageWrapper.style.cssText = 'margin-bottom: 24px; background: #ffffff; box-shadow: 0 6px 24px rgba(0,0,0,0.35); border-radius: 4px; overflow: hidden; display: flex; flex-direction: column; align-items: center; max-width: 98%;';

        const canvas = document.createElement('canvas');
        canvas.width = Math.floor(viewport.width * outputScale);
        canvas.height = Math.floor(viewport.height * outputScale);
        canvas.style.width = Math.floor(viewport.width) + 'px';
        canvas.style.maxWidth = '100%';
        canvas.style.height = 'auto';
        canvas.style.display = 'block';

        const ctx = canvas.getContext('2d');
        const transform = outputScale !== 1 ? [outputScale, 0, 0, outputScale, 0, 0] : null;

        const footer = document.createElement('div');
        footer.style.cssText = 'padding: 6px 14px; font-size: 0.72rem; color: #64748b; background: #f8fafc; width: 100%; border-top: 1px solid #e2e8f0; display: flex; justify-content: space-between; align-items: center; box-sizing: border-box;';
        footer.innerHTML = `<span>Page ${{pageNum}} of ${{pdfDoc.numPages}}</span><span style="font-size:0.68rem; color:#94a3b8;">HMF Orange Report</span>`;

        pageWrapper.appendChild(canvas);
        pageWrapper.appendChild(footer);
        container.appendChild(pageWrapper);

        await page.render({{
          canvasContext: ctx,
          viewport: viewport,
          transform: transform
        }}).promise;
      }}
    }}

    function zoomOrangePdf(delta) {{
      if (!currentPdfDoc) return;
      currentPdfScale = Math.max(0.6, Math.min(2.5, currentPdfScale + delta));
      const zoomLevelElem = document.getElementById('orange-pdf-zoom-level');
      if (zoomLevelElem) zoomLevelElem.textContent = Math.round(currentPdfScale * 100) + '%';
      renderOrangePdfPages(currentPdfDoc, currentPdfScale);
    }}

    function showOrangePdfModal(tumorId, patientId) {{
      const modal = document.getElementById('orange-pdf-modal');
      const titleElem = document.getElementById('orange-pdf-modal-title');
      const canvasContainer = document.getElementById('orange-pdf-canvas-container');
      const frame = document.getElementById('orange-pdf-frame');
      const downloadBtn = document.getElementById('orange-pdf-download-btn');
      const newtabBtn = document.getElementById('orange-pdf-newtab-btn');
      const zoomLevelElem = document.getElementById('orange-pdf-zoom-level');
      if (!modal) return;

      const pStr = patientId ? `<span class="badge badge-primary" style="font-size:0.8rem;">Patient: ${{patientId}}</span>` : '';
      titleElem.innerHTML = `
        <span style="font-weight:700; font-size:1.05rem; color:#0f172a;">📄 HMF Orange Clinical Molecular Report</span>
        <span class="badge badge-normal" style="font-weight:600;">${{tumorId}}</span>
        ${{pStr}}
        <span id="orange-pdf-spinner" style="font-size:0.8rem; color:#64748b; margin-left:8px; display:inline-flex; align-items:center; gap:4px;">⏳ Rendering PDF document...</span>
      `;

      if (canvasContainer) canvasContainer.innerHTML = '';
      if (frame) {{
        frame.style.display = 'none';
        try {{ frame.src = 'about:blank'; }} catch(e) {{}}
      }}
      modal.style.display = 'flex';

      setTimeout(async () => {{
        const spinner = document.getElementById('orange-pdf-spinner');
        try {{
          const b64 = getOrangePdfBase64(tumorId);
          let uint8Bytes = null;
          let pdfUrl = null;

          if (b64) {{
            uint8Bytes = orangeBytesCache[tumorId] || b64ToUint8Array(b64);
            orangeBytesCache[tumorId] = uint8Bytes;
            pdfUrl = getOrangePdfBlobUrl(tumorId);
          }} else {{
            const t = (COHORT_DATA.all_tumors || []).find(x => x.tumor_id === tumorId || x.raw_tumor_id === tumorId);
            if (t && t.orange_pdf_rel) {{
              pdfUrl = t.orange_pdf_rel;
            }}
          }}

          if (!pdfUrl && !uint8Bytes) {{
            if (spinner) spinner.innerHTML = '<span style="color:#ef4444;">⚠️ Orange PDF report not available</span>';
            alert('Orange PDF report not available for ' + tumorId);
            return;
          }}

          // Wire up download button safely
          if (downloadBtn && pdfUrl) {{
            downloadBtn.href = pdfUrl;
            downloadBtn.download = `${{tumorId}}.orange.pdf`;
          }}

          // Wire up open-in-new-tab button safely
          if (newtabBtn && pdfUrl) {{
            newtabBtn.onclick = function() {{
              try {{
                window.open(pdfUrl, '_blank');
              }} catch (openErr) {{
                console.error('Failed opening window:', openErr);
              }}
            }};
          }}

          // 1. PRIMARY RENDERER: PDF.js to HTML5 Canvas (100% works in OneDrive / sandboxed iframes)
          if (window.pdfjsLib && uint8Bytes) {{
            try {{
              // Disable worker in sandboxed / about:srcdoc contexts to avoid origin restrictions
              pdfjsLib.GlobalWorkerOptions.workerSrc = '';
              const loadingTask = pdfjsLib.getDocument({{ data: uint8Bytes }});
              currentPdfDoc = await loadingTask.promise;
              currentPdfScale = 1.35;
              if (zoomLevelElem) zoomLevelElem.textContent = Math.round(currentPdfScale * 100) + '%';
              await renderOrangePdfPages(currentPdfDoc, currentPdfScale);
              if (spinner) spinner.style.display = 'none';
              return; // Rendered successfully with PDF.js!
            }} catch (pdfjsErr) {{
              console.warn('PDF.js canvas rendering failed, trying iframe fallback:', pdfjsErr);
            }}
          }}

          // 2. SECONDARY FALLBACK: native iframe (works on local file:// or open web)
          if (frame && pdfUrl) {{
            try {{
              frame.style.display = 'block';
              frame.src = pdfUrl;
              if (spinner) spinner.style.display = 'none';
            }} catch (frameErr) {{
              console.warn('Native iframe embedding blocked by environment:', frameErr);
              frame.style.display = 'none';
              if (canvasContainer) {{
                canvasContainer.innerHTML = `
                  <div style="background:#ffffff; padding:2rem; border-radius:8px; text-align:center; max-width:600px; margin-top:2rem;">
                    <div style="font-size:2.5rem; margin-bottom:1rem;">📄</div>
                    <h3 style="margin-bottom:0.5rem; color:#0f172a;">PDF Preview Restricted by Cloud Viewer</h3>
                    <p style="color:#64748b; font-size:0.9rem; margin-bottom:1.5rem;">
                      The OneDrive / SharePoint security sandbox does not permit embedded PDF plugins.
                    </p>
                    <div style="display:flex; justify-content:center; gap:12px;">
                      <a href="${{pdfUrl}}" target="_blank" class="btn btn-primary" style="text-decoration:none; padding:8px 16px; border-radius:6px; font-weight:600;">↗ Open PDF in Dedicated Tab</a>
                      <a href="${{pdfUrl}}" download="${{tumorId}}.orange.pdf" class="btn btn-secondary" style="text-decoration:none; padding:8px 16px; border-radius:6px; font-weight:600; border:1px solid #cbd5e1; background:#f8fafc; color:#0f172a;">⬇ Download PDF</a>
                    </div>
                  </div>
                `;
              }}
              if (spinner) spinner.style.display = 'none';
            }}
          }}
        }} catch (err) {{
          console.error('Error displaying Orange PDF:', err);
          if (spinner) spinner.innerHTML = `<span style="color:#ef4444;">⚠️ Error loading PDF: ${{err.message}}</span>`;
        }}
      }}, 30);
    }}

    function closeOrangePdfModal() {{
      const modal = document.getElementById('orange-pdf-modal');
      const canvasContainer = document.getElementById('orange-pdf-canvas-container');
      const frame = document.getElementById('orange-pdf-frame');
      if (modal) modal.style.display = 'none';
      if (frame) {{
        frame.style.display = 'none';
        try {{ frame.src = 'about:blank'; }} catch(e) {{}}
      }}
      if (canvasContainer) canvasContainer.innerHTML = '';
      currentPdfDoc = null;
    }}

    // Cache for decoded ichorCNA PDF Blob URLs and raw bytes
    const ichorBlobUrlCache = {{}};
    const ichorBytesCache = {{}};
    let currentIchorPdfDoc = null;
    let currentIchorPdfScale = 1.35;
    let currentIchorModalSampleId = null;
    let currentIchorModalPatientId = null;
    let currentIchorModalTab = null;

    function getIchorPdfBase64(sampleId, patientId, plotType) {{
      let pl = null;
      if (patientId) {{
        pl = (COHORT_DATA.all_plasmas || []).find(x => x.patient_id === patientId && (x.sample_id === sampleId || x.clean_sample_name === sampleId));
      }}
      if (!pl) {{
        pl = (COHORT_DATA.all_plasmas || []).find(x => x.sample_id === sampleId || x.clean_sample_name === sampleId);
      }}

      let t = null;
      if (patientId) {{
        t = (COHORT_DATA.all_tumors || []).find(x => x.patient_id === patientId && (x.tumor_id === sampleId || x.raw_tumor_id === sampleId));
      }}
      if (!t) {{
        t = (COHORT_DATA.all_tumors || []).find(x => x.tumor_id === sampleId || x.raw_tumor_id === sampleId);
      }}

      const rec = pl || t;
      if (rec) {{
        if (plotType === 'short' && rec.ichor_short_pdf_b64) return rec.ichor_short_pdf_b64;
        if (plotType === 'standard' && rec.ichor_standard_pdf_b64) return rec.ichor_standard_pdf_b64;
        if (plotType === 'short' && !rec.ichor_short_pdf_b64 && rec.ichor_standard_pdf_b64) return rec.ichor_standard_pdf_b64;
        if (plotType === 'standard' && !rec.ichor_standard_pdf_b64 && rec.ichor_short_pdf_b64) return rec.ichor_short_pdf_b64;
        if (rec.ichor_pdf_b64) return rec.ichor_pdf_b64;
      }}

      if (COHORT_DATA.ichor_reports_b64) {{
        if (patientId) {{
          if (plotType === 'short' && COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}_short`]) return COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}_short`];
          if (plotType === 'standard' && COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}_standard`]) return COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}_standard`];
          if (COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}`]) return COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}`];
        }}
        if (plotType === 'short' && COHORT_DATA.ichor_reports_b64[`${{sampleId}}_short`]) return COHORT_DATA.ichor_reports_b64[`${{sampleId}}_short`];
        if (plotType === 'standard' && COHORT_DATA.ichor_reports_b64[`${{sampleId}}_standard`]) return COHORT_DATA.ichor_reports_b64[`${{sampleId}}_standard`];
        if (COHORT_DATA.ichor_reports_b64[sampleId]) return COHORT_DATA.ichor_reports_b64[sampleId];
      }}
      return null;
    }}

    function getIchorPdfBlobUrl(sampleId, patientId, plotType) {{
      const cacheKey = `${{patientId || ''}}_${{sampleId}}_${{plotType || 'auto'}}`;
      if (ichorBlobUrlCache[cacheKey]) {{
        return ichorBlobUrlCache[cacheKey];
      }}

      const b64 = getIchorPdfBase64(sampleId, patientId, plotType);
      if (b64) {{
        try {{
          const bytes = ichorBytesCache[cacheKey] || b64ToUint8Array(b64);
          ichorBytesCache[cacheKey] = bytes;
          const blob = new Blob([bytes], {{ type: 'application/pdf' }});
          const blobUrl = URL.createObjectURL(blob);
          ichorBlobUrlCache[cacheKey] = blobUrl;
          return blobUrl;
        }} catch (err) {{
          console.warn('Blob URL creation not supported, falling back to data URL:', err);
          return 'data:application/pdf;base64,' + b64;
        }}
      }}

      // Fallback to relative URL
      let pl = null;
      if (patientId) pl = (COHORT_DATA.all_plasmas || []).find(x => x.patient_id === patientId && (x.sample_id === sampleId || x.clean_sample_name === sampleId));
      if (!pl) pl = (COHORT_DATA.all_plasmas || []).find(x => x.sample_id === sampleId || x.clean_sample_name === sampleId);

      let t = null;
      if (patientId) t = (COHORT_DATA.all_tumors || []).find(x => x.patient_id === patientId && (x.tumor_id === sampleId || x.raw_tumor_id === sampleId));
      if (!t) t = (COHORT_DATA.all_tumors || []).find(x => x.tumor_id === sampleId || x.raw_tumor_id === sampleId);

      const rec = pl || t;
      if (rec) {{
        if (plotType === 'short' && rec.ichor_short_pdf_rel) return rec.ichor_short_pdf_rel;
        if (plotType === 'standard' && rec.ichor_standard_pdf_rel) return rec.ichor_standard_pdf_rel;
        if (rec.ichor_pdf_rel) return rec.ichor_pdf_rel;
      }}
      return null;
    }}

    async function renderIchorPdfPages(pdfDoc, scale) {{
      const container = document.getElementById('ichor-pdf-canvas-container');
      if (!container || !pdfDoc) return;
      container.innerHTML = '';

      const outputScale = window.devicePixelRatio || 1;

      for (let pageNum = 1; pageNum <= pdfDoc.numPages; pageNum++) {{
        const page = await pdfDoc.getPage(pageNum);
        const viewport = page.getViewport({{ scale: scale }});

        const pageWrapper = document.createElement('div');
        pageWrapper.style.cssText = 'margin-bottom: 24px; background: #ffffff; box-shadow: 0 6px 24px rgba(0,0,0,0.35); border-radius: 4px; overflow: hidden; display: flex; flex-direction: column; align-items: center; max-width: 98%;';

        const canvas = document.createElement('canvas');
        canvas.width = Math.floor(viewport.width * outputScale);
        canvas.height = Math.floor(viewport.height * outputScale);
        canvas.style.width = Math.floor(viewport.width) + 'px';
        canvas.style.maxWidth = '100%';
        canvas.style.height = 'auto';
        canvas.style.display = 'block';

        const ctx = canvas.getContext('2d');
        const transform = outputScale !== 1 ? [outputScale, 0, 0, outputScale, 0, 0] : null;

        const footer = document.createElement('div');
        footer.style.cssText = 'padding: 6px 14px; font-size: 0.72rem; color: #64748b; background: #f8fafc; width: 100%; border-top: 1px solid #e2e8f0; display: flex; justify-content: space-between; align-items: center; box-sizing: border-box;';
        footer.innerHTML = `<span>Page ${{pageNum}} of ${{pdfDoc.numPages}}</span><span style="font-size:0.68rem; color:#94a3b8;">1X ichorCNA Genome-Wide Plot</span>`;

        pageWrapper.appendChild(canvas);
        pageWrapper.appendChild(footer);
        container.appendChild(pageWrapper);

        await page.render({{
          canvasContext: ctx,
          viewport: viewport,
          transform: transform
        }}).promise;
      }}
    }}

    function zoomIchorPdf(delta) {{
      if (!currentIchorPdfDoc) return;
      currentIchorPdfScale = Math.max(0.6, Math.min(2.5, currentIchorPdfScale + delta));
      const zoomLevelElem = document.getElementById('ichor-pdf-zoom-level');
      if (zoomLevelElem) zoomLevelElem.textContent = Math.round(currentIchorPdfScale * 100) + '%';
      renderIchorPdfPages(currentIchorPdfDoc, currentIchorPdfScale);
    }}

    function showIchorPdfModal(sampleId, patientId, preferredTab) {{
      const modal = document.getElementById('ichor-pdf-modal');
      const sampleBadge = document.getElementById('ichor-pdf-modal-sample-badge');
      const patientBadge = document.getElementById('ichor-pdf-modal-patient-badge');
      const tabGroup = document.getElementById('ichor-modal-tab-group');
      const tabShort = document.getElementById('btn-ichor-tab-short');
      const tabStandard = document.getElementById('btn-ichor-tab-standard');
      const canvasContainer = document.getElementById('ichor-pdf-canvas-container');
      const frame = document.getElementById('ichor-pdf-frame');
      if (!modal) return;

      currentIchorModalSampleId = sampleId;
      currentIchorModalPatientId = patientId || null;

      // Check available plots
      let pl = null;
      if (patientId) pl = (COHORT_DATA.all_plasmas || []).find(x => x.patient_id === patientId && (x.sample_id === sampleId || x.clean_sample_name === sampleId));
      if (!pl) pl = (COHORT_DATA.all_plasmas || []).find(x => x.sample_id === sampleId || x.clean_sample_name === sampleId);

      let t = null;
      if (patientId) t = (COHORT_DATA.all_tumors || []).find(x => x.patient_id === patientId && (x.tumor_id === sampleId || x.raw_tumor_id === sampleId));
      if (!t) t = (COHORT_DATA.all_tumors || []).find(x => x.tumor_id === sampleId || x.raw_tumor_id === sampleId);

      if (sampleBadge) {{
        const isTumor = Boolean(t);
        const typeLabel = isTumor ? 'Tumor' : (pl ? 'Plasma' : 'Sample');
        sampleBadge.textContent = `${{typeLabel}}: ${{sampleId}}`;
        sampleBadge.style.display = 'inline-block';
      }}
      if (patientBadge) {{
        if (patientId) {{
          patientBadge.textContent = `Patient: ${{patientId}}`;
          patientBadge.style.display = 'inline-block';
        }} else {{
          patientBadge.style.display = 'none';
        }}
      }}

      const rec = pl || t;
      const hasShort = Boolean(
        (rec && (rec.has_ichor_short || rec.ichor_short_pdf_b64 || rec.ichor_short_pdf_rel)) ||
        (COHORT_DATA.ichor_reports_b64 && (
          (patientId && COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}_short`]) ||
          COHORT_DATA.ichor_reports_b64[`${{sampleId}}_short`]
        ))
      );
      const hasStandard = Boolean(
        (rec && (rec.has_ichor_standard || rec.ichor_standard_pdf_b64 || rec.ichor_standard_pdf_rel)) ||
        (COHORT_DATA.ichor_reports_b64 && (
          (patientId && COHORT_DATA.ichor_reports_b64[`${{patientId}}_${{sampleId}}_standard`]) ||
          COHORT_DATA.ichor_reports_b64[`${{sampleId}}_standard`]
        ))
      );

      // Show/hide tab group
      if (tabGroup) {{
        if (hasShort && hasStandard) {{
          tabGroup.style.display = 'inline-flex';
          if (tabShort) tabShort.style.display = 'inline-block';
          if (tabStandard) tabStandard.style.display = 'inline-block';
        }} else {{
          tabGroup.style.display = 'none';
        }}
      }}

      let activeTab = preferredTab;
      if (!activeTab) {{
        if (hasShort) {{
          activeTab = 'short';
        }} else {{
          activeTab = 'standard';
        }}
      }} else if (activeTab === 'short' && !hasShort && hasStandard) {{
        activeTab = 'standard';
      }} else if (activeTab === 'standard' && !hasStandard && hasShort) {{
        activeTab = 'short';
      }}

      currentIchorModalTab = activeTab;

      if (canvasContainer) canvasContainer.innerHTML = '';
      if (frame) {{
        frame.style.display = 'none';
        try {{ frame.src = 'about:blank'; }} catch(e) {{}}
      }}
      modal.style.display = 'flex';

      loadIchorModalPdf(sampleId, patientId, activeTab);
    }}

    function switchIchorModalTab(tabType) {{
      if (!currentIchorModalSampleId) return;
      currentIchorModalTab = tabType;
      loadIchorModalPdf(currentIchorModalSampleId, currentIchorModalPatientId, tabType);
    }}

    async function loadIchorModalPdf(sampleId, patientId, plotType) {{
      const tabShort = document.getElementById('btn-ichor-tab-short');
      const tabStandard = document.getElementById('btn-ichor-tab-standard');
      const spinner = document.getElementById('ichor-pdf-spinner');
      const canvasContainer = document.getElementById('ichor-pdf-canvas-container');
      const frame = document.getElementById('ichor-pdf-frame');
      const downloadBtn = document.getElementById('ichor-pdf-download-btn');
      const newtabBtn = document.getElementById('ichor-pdf-newtab-btn');
      const zoomLevelElem = document.getElementById('ichor-pdf-zoom-level');

      if (spinner) {{
        spinner.style.display = 'inline-flex';
        spinner.innerHTML = '⏳ Rendering PDF document...';
      }}

      // Update tab active classes/styles
      if (tabShort) {{
        const isShort = (plotType === 'short');
        tabShort.style.background = isShort ? '#ffffff' : 'transparent';
        tabShort.style.color = isShort ? '#0f172a' : '#475569';
        tabShort.style.boxShadow = isShort ? '0 1px 3px rgba(0,0,0,0.1)' : 'none';
      }}
      if (tabStandard) {{
        const isStd = (plotType === 'standard');
        tabStandard.style.background = isStd ? '#ffffff' : 'transparent';
        tabStandard.style.color = isStd ? '#0f172a' : '#475569';
        tabStandard.style.boxShadow = isStd ? '0 1px 3px rgba(0,0,0,0.1)' : 'none';
      }}

      try {{
        const cacheKey = `${{patientId || ''}}_${{sampleId}}_${{plotType || 'auto'}}`;
        const b64 = getIchorPdfBase64(sampleId, patientId, plotType);
        let uint8Bytes = null;
        let pdfUrl = null;

        if (b64) {{
          uint8Bytes = ichorBytesCache[cacheKey] || b64ToUint8Array(b64);
          ichorBytesCache[cacheKey] = uint8Bytes;
          pdfUrl = getIchorPdfBlobUrl(sampleId, patientId, plotType);
        }} else {{
          pdfUrl = getIchorPdfBlobUrl(sampleId, patientId, plotType);
        }}

        if (!pdfUrl && !uint8Bytes) {{
          if (spinner) spinner.innerHTML = `<span style="color:#ef4444;">⚠️ ichorCNA ${{plotType || ''}} PDF plot not available</span>`;
          return;
        }}

        const plotSuffix = plotType === 'short' ? '_short' : '';
        const dlFilename = `${{patientId ? patientId + '_' : ''}}${{sampleId}}${{plotSuffix}}_genomeWide.pdf`;

        if (downloadBtn && pdfUrl) {{
          downloadBtn.href = pdfUrl;
          downloadBtn.download = dlFilename;
        }}

        if (newtabBtn && pdfUrl) {{
          newtabBtn.onclick = function() {{
            try {{
              window.open(pdfUrl, '_blank');
            }} catch (openErr) {{
              console.error('Failed opening window:', openErr);
            }}
          }};
        }}

        if (window.pdfjsLib && uint8Bytes) {{
          try {{
            pdfjsLib.GlobalWorkerOptions.workerSrc = '';
            const loadingTask = pdfjsLib.getDocument({{ data: uint8Bytes }});
            currentIchorPdfDoc = await loadingTask.promise;
            currentIchorPdfScale = 1.35;
            if (zoomLevelElem) zoomLevelElem.textContent = Math.round(currentIchorPdfScale * 100) + '%';
            await renderIchorPdfPages(currentIchorPdfDoc, currentIchorPdfScale);
            if (spinner) spinner.style.display = 'none';
            return;
          }} catch (pdfjsErr) {{
            console.warn('PDF.js canvas rendering failed, trying iframe fallback:', pdfjsErr);
          }}
        }}

        if (frame && pdfUrl) {{
          try {{
            frame.style.display = 'block';
            frame.src = pdfUrl;
            if (spinner) spinner.style.display = 'none';
          }} catch (frameErr) {{
            console.warn('Native iframe embedding blocked by environment:', frameErr);
            frame.style.display = 'none';
            if (canvasContainer) {{
              canvasContainer.innerHTML = `
                <div style="background:#ffffff; padding:2rem; border-radius:8px; text-align:center; max-width:600px; margin-top:2rem;">
                  <div style="font-size:2.5rem; margin-bottom:1rem;">📊</div>
                  <h3 style="margin-bottom:0.5rem; color:#0f172a;">PDF Preview Restricted by Cloud Viewer</h3>
                  <p style="color:#64748b; font-size:0.9rem; margin-bottom:1.5rem;">
                    The OneDrive / SharePoint security sandbox does not permit embedded PDF plugins.
                  </p>
                  <div style="display:flex; justify-content:center; gap:12px;">
                    <a href="${{pdfUrl}}" target="_blank" class="btn btn-primary" style="text-decoration:none; padding:8px 16px; border-radius:6px; font-weight:600; background:#d97706; border-color:#b45309;">↗ Open PDF in Dedicated Tab</a>
                    <a href="${{pdfUrl}}" download="${{dlFilename}}" class="btn btn-secondary" style="text-decoration:none; padding:8px 16px; border-radius:6px; font-weight:600; border:1px solid #cbd5e1; background:#f8fafc; color:#0f172a;">⬇ Download PDF</a>
                  </div>
                </div>
              `;
            }}
            if (spinner) spinner.style.display = 'none';
          }}
        }}
      }} catch (err) {{
        console.error('Error displaying ichorCNA PDF:', err);
        if (spinner) spinner.innerHTML = `<span style="color:#ef4444;">⚠️ Error loading PDF: ${{err.message}}</span>`;
      }}
    }}

    function closeIchorPdfModal() {{
      const modal = document.getElementById('ichor-pdf-modal');
      const canvasContainer = document.getElementById('ichor-pdf-canvas-container');
      const frame = document.getElementById('ichor-pdf-frame');
      if (modal) modal.style.display = 'none';
      if (frame) {{
        frame.style.display = 'none';
        try {{ frame.src = 'about:blank'; }} catch(e) {{}}
      }}
      if (canvasContainer) canvasContainer.innerHTML = '';
      currentIchorPdfDoc = null;
    }}

    function viewSelectedPatientTrajectory() {{
      const p = document.getElementById('mrd-patient-select').value;
      if (p && p !== 'ALL') {{
        showPatientTrajectoryModal(p);
      }}
    }}

    function populateDriversTable() {{
      const dtbody = document.getElementById('drivers-tbody');
      if (!dtbody) return;
      dtbody.innerHTML = '';

      COHORT_DATA.all_tumors.forEach(t => {{
        (t.drivers || []).forEach(d => {{
          const row = document.createElement('tr');
          row.innerHTML = `
            <td><strong>${{t.patient_id}}</strong></td>
            <td>${{t.tumor_id}}</td>
            <td style="color:#58a6ff; font-weight:bold;">${{d.gene}}</td>
            <td><span class="badge ${{d.alt_type === 'Truncating' ? 'badge-danger' : 'badge-primary'}}">${{d.alt_type}}</span></td>
            <td>${{d.hgvs_protein || '--'}}</td>
            <td>${{d.hgvs_coding || '--'}}</td>
            <td>${{d.vaf ? (d.vaf * 100).toFixed(1) + '%' : '--'}}</td>
            <td><span class="badge badge-success">SOMATIC</span></td>
          `;
          dtbody.appendChild(row);
        }});

        (t.germline_drivers || []).forEach(gd => {{
          const row = document.createElement('tr');
          row.innerHTML = `
            <td><strong>${{t.patient_id}}</strong></td>
            <td>${{t.tumor_id}}</td>
            <td style="color:#58a6ff; font-weight:bold;">${{gd.gene}}</td>
            <td><span class="badge badge-warning">${{gd.alt_type}}</span></td>
            <td>${{gd.hgvs_protein || '--'}}</td>
            <td>${{gd.hgvs_coding || '--'}}</td>
            <td>${{gd.vaf ? (gd.vaf * 100).toFixed(1) + '%' : '50.0%'}}</td>
            <td><span class="badge" style="background:#8b5cf6; color:#fff;">GERMLINE</span></td>
          `;
          dtbody.appendChild(row);
        }});
      }});
    }}

    function addCustomGene(customSymbol) {{
      const input = document.getElementById('gene-add-input');
      const gene = (customSymbol || (input ? input.value : '')).trim().toUpperCase();
      if (!gene) return;

      if (!oncoGenes.includes(gene)) {{
        oncoGenes.unshift(gene);
      }} else {{
        // Move to the top if already present so it is immediately visible
        oncoGenes = [gene, ...oncoGenes.filter(g => g !== gene)];
      }}

      if (input) input.value = '';
      renderOncoprint(999);

      // Temporary visual highlight and scroll into view
      setTimeout(() => {{
        const rows = document.querySelectorAll('.oncoprint-table tbody tr');
        rows.forEach(r => {{
          const nameTd = r.querySelector('.gene-name');
          if (nameTd && nameTd.textContent.trim().toUpperCase() === gene) {{
            r.style.outline = '2px solid #58a6ff';
            r.scrollIntoView({{ behavior: 'smooth', block: 'nearest' }});
            setTimeout(() => {{ r.style.outline = 'none'; }}, 2200);
          }}
        }});
      }}, 100);
    }}

    // MRD Tab & Multi-Modal Log-Scale Plotting with Sample-Specific LODs
    let mrdCharts = {{ snv: null, cnv: null, loh: null, ichor: null }};
    let currentMrdModality = 'ALL';

    function setMrdModality(mod) {{
      currentMrdModality = mod;
      ['all', 'snv', 'cnv', 'loh', 'ichor'].forEach(m => {{
        const btn = document.getElementById('btn-mrd-mod-' + m);
        if (btn) {{
          if (m === mod.toLowerCase()) btn.classList.add('active');
          else btn.classList.remove('active');
        }}
      }});
      const tSnv = document.getElementById('track-card-snv');
      const tCnv = document.getElementById('track-card-cnv');
      const tLoh = document.getElementById('track-card-loh');
      const tIchor = document.getElementById('track-card-ichor');
      if (tSnv) tSnv.style.display = (mod === 'ALL' || mod === 'SNV') ? 'block' : 'none';
      if (tCnv) tCnv.style.display = (mod === 'ALL' || mod === 'CNV') ? 'block' : 'none';
      if (tLoh) tLoh.style.display = (mod === 'ALL' || mod === 'LOH') ? 'block' : 'none';
      if (tIchor) tIchor.style.display = (mod === 'ALL' || mod === 'ICHOR') ? 'block' : 'none';
      updateMrdTrajectoryPlot();
    }}

    // Sort cohort plasmas hierarchically: SNV MRD+ -> CNV MRD+ -> LOH MRD+ -> Triple Negative
    function getSortedPlasmas(plasmas) {{
      return [...plasmas].sort((a, b) => {{
        // 1. SNV MRD positive (highest purity first)
        if (a.snv_mrd !== b.snv_mrd) return a.snv_mrd ? -1 : 1;
        if (a.snv_mrd && b.snv_mrd) return (b.snv_purity_pct || 0) - (a.snv_purity_pct || 0);

        // 2. CNV MRD positive (highest purity first)
        if (a.cnv_mrd !== b.cnv_mrd) return a.cnv_mrd ? -1 : 1;
        if (a.cnv_mrd && b.cnv_mrd) return (b.cnv_purity_pct || 0) - (a.cnv_purity_pct || 0);

        // 3. LOH MRD positive (highest purity first)
        if (a.loh_mrd !== b.loh_mrd) return a.loh_mrd ? -1 : 1;
        if (a.loh_mrd && b.loh_mrd) return (b.loh_purity_pct || 0) - (a.loh_purity_pct || 0);

        // 4. Triple negative: sort by patient_id
        return (a.patient_id || '').localeCompare(b.patient_id || '');
      }});
    }}

    // Custom Chart.js Plugin to match Trajectory Plot standards:
    // - Continuous shaded grey region below LOD
    // - Continuous dashed LOD line (NO points!)
    // - Dotted vertical WGS baseline timestamp line
    // - Dashed baseline horizontal reference line at 0.05%
    // - Continuous 95% Confidence Interval ribbon between detected points (NO per-point boxes!)
    function createSampleLodPlugin(lodKey, metricKey, mrdKey, lowKey, highKey, ciColor, currentPlasmas, showLods) {{
      return {{
        id: 'sampleLodShading_' + lodKey,
        beforeDatasetsDraw: (chart) => {{
          const {{ ctx, chartArea, scales: {{ x, y }} }} = chart;
          if (!chartArea) return;
          const {{ bottom, top, left, right }} = chartArea;
          const n = currentPlasmas.length;
          if (n === 0) return;

          ctx.save();

          // 1. Below Limit of Detection: shaded grey region
          if (showLods) {{
            ctx.beginPath();
            ctx.moveTo(x.getPixelForTick(0), bottom);
            currentPlasmas.forEach((pl, idx) => {{
              const rawLod = pl[lodKey] || 0.002;
              const lodClamped = Math.min(Math.max(rawLod, 0.0005), 100.0);
              ctx.lineTo(x.getPixelForTick(idx), y.getPixelForValue(lodClamped));
            }});
            ctx.lineTo(x.getPixelForTick(n - 1), bottom);
            ctx.closePath();
            ctx.fillStyle = 'rgba(226, 232, 240, 0.85)';
            ctx.fill();

            // 2. Limit of Detection continuous dashed line (NO points, clean line only!)
            ctx.beginPath();
            currentPlasmas.forEach((pl, idx) => {{
              const rawLod = pl[lodKey] || 0.002;
              const lodClamped = Math.min(Math.max(rawLod, 0.0005), 100.0);
              const xPos = x.getPixelForTick(idx);
              const yPos = y.getPixelForValue(lodClamped);
              if (idx === 0) ctx.moveTo(xPos, yPos);
              else ctx.lineTo(xPos, yPos);
            }});
            ctx.strokeStyle = '#94a3b8';
            ctx.lineWidth = 1.3;
            ctx.setLineDash([5, 4]);
            ctx.stroke();
          }}

          // 3. WGS Baseline Timestamp line (dotted vertical line at index 0)
          ctx.beginPath();
          const x0 = x.getPixelForTick(0);
          ctx.moveTo(x0, top);
          ctx.lineTo(x0, bottom);
          ctx.strokeStyle = '#000000';
          ctx.lineWidth = 1.5;
          ctx.setLineDash([2, 3]);
          ctx.stroke();

          // 4. Baseline horizontal dashed reference line (0.05%)
          ctx.beginPath();
          const yRef = y.getPixelForValue(0.05);
          ctx.moveTo(left, yRef);
          ctx.lineTo(right, yRef);
          ctx.strokeStyle = '#9ca3af';
          ctx.lineWidth = 0.8;
          ctx.setLineDash([4, 4]);
          ctx.stroke();

          // 5. Confidence Interval Ribbon across all timepoints (including negative/below-LOD samples)
          if (metricKey && currentPlasmas.length >= 2) {{
            const ciValues = currentPlasmas.map((pl, idx) => {{
              const isDet = pl[mrdKey];
              const val = pl[metricKey] || 0;
              const lod = pl[lodKey] || 0.002;
              const isDetected = Boolean(isDet && (val > lod));
              const cLow = pl[lowKey];
              const cHigh = pl[highKey];

              let lowVal, highVal;
              if (isDetected) {{
                lowVal = (cLow !== null && cLow !== undefined && cLow > 0) ? cLow : Math.max(val * 0.75, lod);
                highVal = (cHigh !== null && cHigh !== undefined && cHigh > 0) ? cHigh : Math.min(val * 1.25, 100.0);
              }} else {{
                highVal = (cHigh !== null && cHigh !== undefined && cHigh > 0) ? Math.max(cHigh, lod) : lod;
                lowVal = (cLow !== null && cLow !== undefined && cLow > 0) ? cLow : 0.000575;
                if (lowVal > highVal) lowVal = highVal * 0.5;
              }}
              return {{ low: lowVal, high: highVal }};
            }});

            for (let i = 0; i < currentPlasmas.length - 1; i++) {{
              const ci1 = ciValues[i];
              const ci2 = ciValues[i + 1];

              const x1 = x.getPixelForTick(i);
              const x2 = x.getPixelForTick(i + 1);
              const yHigh1 = y.getPixelForValue(ci1.high);
              const yHigh2 = y.getPixelForValue(ci2.high);
              const yLow2 = Math.min(y.getPixelForValue(ci2.low), bottom);
              const yLow1 = Math.min(y.getPixelForValue(ci1.low), bottom);

              ctx.beginPath();
              ctx.moveTo(x1, yHigh1);
              ctx.lineTo(x2, yHigh2);
              ctx.lineTo(x2, yLow2);
              ctx.lineTo(x1, yLow1);
              ctx.closePath();
              ctx.fillStyle = ciColor || 'rgba(254, 226, 226, 0.75)';
              ctx.fill();
            }}
          }}

          ctx.restore();
        }}
      }};
    }}

    function renderMrdSection() {{
      if (!COHORT_DATA.has_mrd) {{
        document.getElementById('mrd-empty-notice').style.display = 'block';
        document.getElementById('mrd-content').style.display = 'none';
        return;
      }}

      // Populate Patient Dropdown
      const pSelect = document.getElementById('mrd-patient-select');
      if (pSelect.children.length === 1) {{
        Object.keys(COHORT_DATA.patients).forEach(pid => {{
          if (COHORT_DATA.patients[pid].plasmas && COHORT_DATA.patients[pid].plasmas.length > 0) {{
            const opt = document.createElement('option');
            opt.value = pid;
            opt.innerText = pid;
            pSelect.appendChild(opt);
          }}
        }});
      }}

      // Populate MRD Table with all modalities & sample-specific LODs (sorted SNV > CNV > LOH)
      const mtbody = document.getElementById('mrd-tbody');
      mtbody.innerHTML = '';
      const sortedPlasmas = getSortedPlasmas(COHORT_DATA.all_plasmas);
      sortedPlasmas.forEach(pl => {{
        const row = document.createElement('tr');
        const badge = pl.is_mrd_positive ? '<span class="badge badge-mrd-pos">Positive</span>' : '<span class="badge badge-mrd-neg">Negative</span>';
        const snvPurStr = pl.snv_purity_pct > 0 ? `<strong>${{pl.snv_purity_pct.toFixed(4)}}%</strong>` : '<span style="color:#8b949e;">0.00% (&lt; LOD)</span>';
        const cnvPurStr = pl.cnv_purity_pct > 0 ? `${{pl.cnv_purity_pct.toFixed(4)}}%` : '<span style="color:#8b949e;">0.00% (&lt; LOD)</span>';
        const lohPurStr = pl.loh_purity_pct > 0 ? `${{pl.loh_purity_pct.toFixed(4)}}%` : '<span style="color:#8b949e;">0.00% (&lt; LOD)</span>';
        const decompBtn = `<button class="btn btn-sm" onclick="showPlasmaModalityModal('${{pl.sample_id}}', '${{pl.patient_id}}')" style="font-size:0.75rem; padding:3px 8px; border-radius:4px; border:1px solid #388bfd; background:#161b22; color:#58a6ff; cursor:pointer; font-weight:600; white-space:nowrap; display:inline-flex; align-items:center; gap:3px;">📊 Decomp</button>`;
        const trajBtn = `<button class="btn btn-sm" onclick="showPatientTrajectoryModal('${{pl.patient_id}}')" style="font-size:0.75rem; padding:3px 8px; border-radius:4px; border:1px solid #10b981; background:#161b22; color:#34d399; cursor:pointer; font-weight:600; white-space:nowrap; display:inline-flex; align-items:center; gap:3px;">📈 Trajectory</button>`;

        const hasIchor = Boolean(
          (COHORT_DATA.ichor_reports_b64 && (COHORT_DATA.ichor_reports_b64[pl.sample_id] || (pl.clean_sample_name && COHORT_DATA.ichor_reports_b64[pl.clean_sample_name]))) ||
          pl.ichor_pdf_b64 ||
          pl.ichor_pdf_rel ||
          pl.has_ichor_pdf
        );
        const ichorBtn = hasIchor
          ? `<button type="button" onclick="showIchorPdfModal('${{pl.sample_id}}', '${{pl.patient_id}}')" class="badge" style="background:#fef3c7; color:#b45309; border:1px solid #f59e0b; cursor:pointer; margin-left:4px; font-size:0.68rem; padding:2px 6px; border-radius:4px; font-weight:600; display:inline-flex; align-items:center; gap:2px;" title="View 1X ichorCNA Genome-Wide Solution Plot in browser">📊 Plot</button>`
          : '';

        const ichorStr = (pl.ichor_tf !== null && pl.ichor_tf !== undefined) ? `<strong style="color:#d97706;">${{(pl.ichor_tf * 100).toFixed(2)}}%</strong>${{ichorBtn}}` : '<span style="color:#8b949e;">--</span>';

        row.innerHTML = `
          <td><strong>${{pl.patient_id}}</strong></td>
          <td>${{pl.sample_id}}</td>
          <td><span class="badge badge-primary">${{pl.timepoint}}</span></td>
          <td>${{badge}}</td>
          <td>${{snvPurStr}}</td>
          <td style="color:#0ea5e9; font-weight:600;">${{pl.snv_lod_pct ? pl.snv_lod_pct.toFixed(4) + '%' : '--'}}</td>
          <td>${{cnvPurStr}}</td>
          <td style="color:#a855f7; font-weight:600;">${{pl.cnv_lod_pct ? pl.cnv_lod_pct.toFixed(4) + '%' : '--'}}</td>
          <td>${{lohPurStr}}</td>
          <td style="color:#f59e0b; font-weight:600;">${{pl.loh_lod_pct ? pl.loh_lod_pct.toFixed(4) + '%' : '--'}}</td>
          <td><span class="badge badge-primary">${{pl.wisp_method || (pl.is_mrd_positive ? 'WISP POS' : 'WISP NEG')}}</span></td>
          <td>${{ichorStr}}</td>
          <td>${{pl.weighted_depth ? pl.weighted_depth.toFixed(1) + 'x' : '--'}}</td>
          <td style="white-space:nowrap;"><div style="display:inline-flex; gap:6px;">${{decompBtn}}${{trajBtn}}</div></td>
        `;
        mtbody.appendChild(row);
      }});

      updateMrdTrajectoryPlot();
    }}

    let currentMrdViewMode = 'CLINICAL';

    function setMrdViewMode(mode) {{
      currentMrdViewMode = mode;
      const btnClin = document.getElementById('btn-mrd-view-clinical');
      const btnTrk = document.getElementById('btn-mrd-view-tracks');
      if (btnClin) btnClin.classList.toggle('active', mode === 'CLINICAL');
      if (btnTrk) btnTrk.classList.toggle('active', mode === 'TRACKS');
      const clinView = document.getElementById('mrd-clinical-view');
      const trkView = document.getElementById('mrd-tracks-view');
      if (clinView) clinView.style.display = (mode === 'CLINICAL') ? 'block' : 'none';
      if (trkView) trkView.style.display = (mode === 'TRACKS') ? 'block' : 'none';
      updateMrdTrajectoryPlot();
    }}

    function updateMrdTrajectoryPlot() {{
      const selectedPatient = document.getElementById('mrd-patient-select').value;
      const showLods = document.getElementById('toggle-sample-lods') ? document.getElementById('toggle-sample-lods').checked : true;

      const trajBtn = document.getElementById('btn-view-patient-trajectory');
      const dlPng = document.getElementById('btn-download-mrd-png');
      const dlPdf = document.getElementById('btn-download-mrd-pdf');

      // Update Clinical Trajectory Mode vs Multi-Modal Tracks Mode
      if (currentMrdViewMode === 'CLINICAL') {{
        const clinicalSingle = document.getElementById('mrd-clinical-single-patient');
        const clinicalGrid = document.getElementById('mrd-clinical-cohort-grid');

        if (selectedPatient === 'ALL') {{
          if (clinicalGrid) clinicalGrid.style.display = 'none';
          if (clinicalSingle) {{
            clinicalSingle.style.display = 'block';
            const imgEl = document.getElementById('mrd-clinical-img');
            const cohortPlotUrl = COHORT_DATA.cohort_mrd_multimodal_plot_b64 || COHORT_DATA.cohort_mrd_multimodal_plot || 'cohort_wgts_mrd_multimodal.png';
            if (imgEl) {{
              imgEl.src = cohortPlotUrl;
              imgEl.alt = 'Cohort Multi-Modal MRD Publication Figure';
            }}
          }}
          if (trajBtn) trajBtn.style.display = 'none';
          if (dlPng) {{
            dlPng.href = COHORT_DATA.cohort_mrd_multimodal_plot || 'cohort_wgts_mrd_multimodal.png';
            dlPng.download = 'cohort_wgts_mrd_multimodal.png';
            dlPng.style.display = 'inline-flex';
          }}
          if (dlPdf) {{
            dlPdf.href = 'cohort_wgts_mrd_multimodal.pdf';
            dlPdf.download = 'cohort_wgts_mrd_multimodal.pdf';
            dlPdf.style.display = 'inline-flex';
          }}
        }} else {{
          const prec = COHORT_DATA.patients[selectedPatient];
          if (clinicalGrid) clinicalGrid.style.display = 'none';
          if (clinicalSingle) {{
            clinicalSingle.style.display = 'block';
            const imgEl = document.getElementById('mrd-clinical-img');
            const plotUrl = prec ? (prec.trajectory_plot_b64 || prec.trajectory_plot || `patient_plots/${{selectedPatient}}_mrd_trajectory.png`) : '';
            if (imgEl) {{
              imgEl.src = plotUrl;
              imgEl.alt = `${{selectedPatient}} Longitudinal Trajectory`;
            }}
          }}
          if (trajBtn) trajBtn.style.display = 'inline-flex';
          if (dlPng && prec) {{
            dlPng.href = prec.trajectory_plot || `patient_plots/${{selectedPatient}}_mrd_trajectory.png`;
            dlPng.download = `${{selectedPatient}}_mrd_trajectory.png`;
            dlPng.style.display = 'inline-flex';
          }}
          if (dlPdf && prec) {{
            dlPdf.href = `patient_plots/${{selectedPatient}}_mrd_trajectory.pdf`;
            dlPdf.download = `${{selectedPatient}}_mrd_trajectory.pdf`;
            dlPdf.style.display = 'inline-flex';
          }}
        }}
      }} else {{
        // Multi-Modal Breakdown Mode
        const tracksSingle = document.getElementById('mrd-tracks-single-patient');
        if (selectedPatient === 'ALL') {{
          if (tracksSingle) {{
            tracksSingle.style.display = 'block';
            const imgEl = document.getElementById('mrd-tracks-img');
            const cohortPlotUrl = COHORT_DATA.cohort_mrd_multimodal_plot_b64 || COHORT_DATA.cohort_mrd_multimodal_plot || 'cohort_wgts_mrd_multimodal.png';
            if (imgEl) {{
              imgEl.src = cohortPlotUrl;
              imgEl.alt = 'Cohort Multi-Modal ctDNA Profile';
            }}
          }}
          if (trajBtn) trajBtn.style.display = 'none';
          if (dlPng) {{
            dlPng.href = COHORT_DATA.cohort_mrd_multimodal_plot || 'cohort_wgts_mrd_multimodal.png';
            dlPng.download = 'cohort_wgts_mrd_multimodal.png';
            dlPng.style.display = 'inline-flex';
          }}
          if (dlPdf) {{
            dlPdf.href = 'cohort_wgts_mrd_multimodal.pdf';
            dlPdf.download = 'cohort_wgts_mrd_multimodal.pdf';
            dlPdf.style.display = 'inline-flex';
          }}
        }} else {{
          const prec = COHORT_DATA.patients[selectedPatient];
          if (tracksSingle) {{
            tracksSingle.style.display = 'block';
            const imgEl = document.getElementById('mrd-tracks-img');
            const plotUrl = prec ? (prec.multimodal_plot_b64 || prec.multimodal_plot || `patient_plots/${{selectedPatient}}_mrd_multimodal.png`) : '';
            if (imgEl) {{
              imgEl.src = plotUrl;
              imgEl.alt = `${{selectedPatient}} Multi-Modal ctDNA Trajectory Breakdown`;
            }}
          }}
          if (trajBtn) trajBtn.style.display = 'inline-flex';
          if (dlPng && prec) {{
            dlPng.href = prec.multimodal_plot || `patient_plots/${{selectedPatient}}_mrd_multimodal.png`;
            dlPng.download = `${{selectedPatient}}_mrd_multimodal.png`;
            dlPng.style.display = 'inline-flex';
          }}
          if (dlPdf && prec) {{
            dlPdf.href = `patient_plots/${{selectedPatient}}_mrd_multimodal.pdf`;
            dlPdf.download = `${{selectedPatient}}_mrd_multimodal.pdf`;
            dlPdf.style.display = 'inline-flex';
          }}
        }}
      }}

      const yMinFloor = 0.0005; // 0.0005% log baseline matching trajectory plot
      let xLabels = [];
      let plasmasToPlot = [];

      if (selectedPatient === 'ALL') {{
        plasmasToPlot = getSortedPlasmas(COHORT_DATA.all_plasmas);
        xLabels = plasmasToPlot.map(pl => `${{pl.patient_id}} (${{pl.timepoint}})`);
      }} else {{
        const prec = COHORT_DATA.patients[selectedPatient];
        plasmasToPlot = prec ? (prec.plasmas || []) : [];
        xLabels = plasmasToPlot.map((pl, idx) => {{
          const dt = pl.sampling_date || pl.date;
          const tp = pl.timepoint;
          if (tp && dt && tp !== dt) return `${{tp}}\\n${{dt}}`;
          if (tp && tp !== 'Plasma ctDNA' && tp !== 'Plasma') return tp;
          if (dt) return dt;
          return pl.clean_sample_name || pl.sample_id || `T${{idx}}`;
        }});
      }}

      // Destroy existing charts
      Object.keys(mrdCharts).forEach(k => {{
        if (mrdCharts[k]) {{
          mrdCharts[k].destroy();
          mrdCharts[k] = null;
        }}
      }});

      if (plasmasToPlot.length === 0) return;

      const activeTracks = [];
      if (currentMrdModality === 'ALL' || currentMrdModality === 'SNV') activeTracks.push('snv');
      if (currentMrdModality === 'ALL' || currentMrdModality === 'CNV') activeTracks.push('cnv');
      if (currentMrdModality === 'ALL' || currentMrdModality === 'LOH') activeTracks.push('loh');
      if (currentMrdModality === 'ALL' || currentMrdModality === 'ICHOR') activeTracks.push('ichor');

      const trackConfigs = {{
        snv: {{
          canvasId: 'chart-mrd-snv',
          title: 'Tumor Fraction (%) — SNV ctDNA',
          metricKey: 'snv_purity_pct',
          mrdKey: 'snv_mrd',
          lodKey: 'snv_lod_pct',
          lowKey: 'snv_purity_low_pct',
          highKey: 'snv_purity_high_pct',
          color: {{ fill: '#ef4444', border: '#b91c1c' }},
          ciColor: 'rgba(254, 226, 226, 0.75)',
          pointStyle: 'circle'
        }},
        cnv: {{
          canvasId: 'chart-mrd-cnv',
          title: 'Tumor Fraction (%) — CNV ctDNA',
          metricKey: 'cnv_purity_pct',
          mrdKey: 'cnv_mrd',
          lodKey: 'cnv_lod_pct',
          lowKey: 'cnv_purity_low_pct',
          highKey: 'cnv_purity_high_pct',
          color: {{ fill: '#a855f7', border: '#7e22ce' }},
          ciColor: 'rgba(243, 232, 255, 0.75)',
          pointStyle: 'circle'
        }},
        loh: {{
          canvasId: 'chart-mrd-loh',
          title: 'Tumor Fraction (%) — LOH ctDNA',
          metricKey: 'loh_purity_pct',
          mrdKey: 'loh_mrd',
          lodKey: 'loh_lod_pct',
          lowKey: 'loh_purity_low_pct',
          highKey: 'loh_purity_high_pct',
          color: {{ fill: '#f59e0b', border: '#b45309' }},
          ciColor: 'rgba(254, 243, 199, 0.75)',
          pointStyle: 'circle'
        }},
        ichor: {{
          canvasId: 'chart-mrd-ichor',
          title: 'Tumor Fraction (%) — 1X ichorCNA (Short Fragments)',
          metricKey: 'ichor_purity_pct',
          mrdKey: 'ichor_mrd',
          lodKey: 'ichor_lod_pct',
          lowKey: null,
          highKey: null,
          color: {{ fill: '#d97706', border: '#92400e' }},
          ciColor: 'rgba(254, 243, 199, 0.5)',
          pointStyle: 'rectRot'
        }}
      }};

      const lastTrackId = activeTracks[activeTracks.length - 1];

      activeTracks.forEach(trackKey => {{
        const cfg = trackConfigs[trackKey];
        const canvas = document.getElementById(cfg.canvasId);
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const isBottomTrack = (trackKey === lastTrackId);

        const points = plasmasToPlot.map((pl, idx) => {{
          const val = pl[cfg.metricKey];
          const isDet = pl[cfg.mrdKey];
          const lod = pl[cfg.lodKey] || 0.002;
          const isDetected = Boolean(isDet && (val > lod));
          const actualVal = Math.max(val || 0.0, 0.0);
          const plotY = (actualVal > 0) ? Math.max(actualVal, yMinFloor * 1.15) : (yMinFloor * 1.15);
          return {{
            x: xLabels[idx],
            y: plotY,
            rawVal: val,
            lod: lod,
            isDetected: isDetected,
            _pl: pl
          }};
        }});

        const datasets = [
          {{
            label: cfg.title,
            data: points,
            borderColor: cfg.color.fill,
            backgroundColor: cfg.color.fill,
            pointBackgroundColor: points.map(p => p.isDetected ? cfg.color.fill : '#6b7280'),
            pointBorderColor: points.map(p => p.isDetected ? cfg.color.border : '#374151'),
            pointBorderWidth: 1.5,
            pointStyle: cfg.pointStyle || 'circle',
            pointRadius: 6,
            pointHoverRadius: 8.5,
            fill: false,
            tension: 0,
            showLine: true,
            segment: {{
              borderColor: ctxSeg => {{
                const p0 = points[ctxSeg.p0DataIndex];
                const p1 = points[ctxSeg.p1DataIndex];
                if (p0 && p1 && p0.isDetected && p1.isDetected) {{
                  return cfg.color.fill;
                }}
                return '#9ca3af';
              }},
              borderDash: ctxSeg => {{
                const p0 = points[ctxSeg.p0DataIndex];
                const p1 = points[ctxSeg.p1DataIndex];
                if (p0 && p1 && p0.isDetected && p1.isDetected) {{
                  return [];
                }}
                return [5, 4];
              }},
              borderWidth: ctxSeg => {{
                const p0 = points[ctxSeg.p0DataIndex];
                const p1 = points[ctxSeg.p1DataIndex];
                return (p0 && p1 && p0.isDetected && p1.isDetected) ? 2.0 : 1.5;
              }}
            }}
          }}
        ];

        mrdCharts[trackKey] = new Chart(ctx, {{
          type: 'line',
          data: {{
            labels: xLabels,
            datasets: datasets
          }},
          options: {{
            responsive: true,
            maintainAspectRatio: false,
            onClick: function(evt, activeEls) {{
              if (activeEls && activeEls.length > 0) {{
                const el = activeEls[0];
                const pt = points[el.index];
                if (pt && pt._pl) {{
                  const sId = pt._pl.sample_id || pt._pl.clean_sample_name;
                  const pId = pt._pl.patient_id;
                  showIchorPdfModal(sId, pId);
                }}
              }}
            }},
            onHover: function(evt, activeEls) {{
              if (evt && evt.native && evt.native.target) {{
                evt.native.target.style.cursor = (activeEls && activeEls.length > 0) ? 'pointer' : 'default';
              }}
            }},
            scales: {{
              y: {{
                type: 'logarithmic',
                min: yMinFloor,
                max: 105.0,
                title: {{
                  display: true,
                  text: cfg.title,
                  color: '#334155',
                  font: {{ size: 10.5, weight: 'bold' }}
                }},
                ticks: {{
                  color: '#64748b',
                  font: {{ size: 9.5 }},
                  callback: function(val) {{
                    if (val === 100 || val === 10 || val === 1 || val === 0.1 || val === 0.01 || val === 0.001) return val + '%';
                    return null;
                  }}
                }},
                grid: {{ color: '#e5e7eb' }}
              }},
              x: {{
                display: isBottomTrack,
                title: {{
                  display: isBottomTrack,
                  text: selectedPatient === 'ALL' ? 'Analyzed Plasma Samples (Sorted: SNV MRD+ > CNV MRD+ > LOH MRD+ > Neg)' : 'Sampling Date / Timepoint',
                  color: '#334155',
                  font: {{ size: 11, weight: 'bold' }}
                }},
                ticks: {{
                  color: '#64748b',
                  font: {{ size: selectedPatient === 'ALL' ? 9.5 : 10.5 }},
                  maxRotation: selectedPatient === 'ALL' ? 35 : (xLabels.length > 6 ? 35 : 0)
                }},
                grid: {{ color: '#e5e7eb' }}
              }}
            }},
            plugins: {{
              legend: {{
                display: false
              }},
              tooltip: {{
                callbacks: {{
                  label: function(ctxTip) {{
                    const raw = ctxTip.raw;
                    if (!raw) return '';
                    const lbl = ctxTip.dataset.label;
                    const valStr = (raw.rawVal && raw.rawVal > 0) ? `${{raw.rawVal.toFixed(4)}}%` : 'Not Detected (< LOD)';
                    const statusStr = raw.isDetected ? 'MRD POSITIVE' : 'MRD NEGATIVE';
                    const lodVal = raw.lod ? raw.lod.toFixed(4) + '%' : 'NA';
                    return `${{lbl}}: ${{valStr}} [${{statusStr}}] (Sample LOD: ${{lodVal}})`;
                  }},
                  afterLabel: function(ctxTip) {{
                    return '💡 Click dot to view 1X ichorCNA solution plot';
                  }}
                }}
              }}
            }}
          }},
          plugins: [
            createSampleLodPlugin(cfg.lodKey, cfg.metricKey, cfg.mrdKey, cfg.lowKey, cfg.highKey, cfg.ciColor, plasmasToPlot, showLods)
          ]
        }});
      }});
    }}

    // Export Table to CSV
    function exportTableToCSV(tableId, filename) {{
      const table = document.getElementById(tableId);
      let csv = [];
      for (let r of table.rows) {{
        let row = [];
        for (let c of r.cells) {{
          row.push('"' + c.innerText.replace(/"/g, '""') + '"');
        }}
        csv.push(row.join(','));
      }}
      const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.click();
    }}

    function initGeneAddControls() {{
      const dlist = document.getElementById('available-genes-list');
      const pillsContainer = document.getElementById('orange-driver-pills');
      
      const geneSet = new Set();
      (COHORT_DATA.all_orange_driver_genes || []).forEach(g => geneSet.add(g));
      (COHORT_DATA.oncoprint.genes || []).forEach(g => geneSet.add(g));
      (COHORT_DATA.all_tumors || []).forEach(t => {{
        (t.drivers || []).forEach(d => {{ if (d.gene) geneSet.add(d.gene); }});
        (t.germline_drivers || []).forEach(gd => {{ if (gd.gene) geneSet.add(gd.gene); }});
        Object.keys(t.cnvs || {{}}).forEach(g => {{
          if (t.cnvs[g] && t.cnvs[g].is_driver) geneSet.add(g);
        }});
      }});

      if (dlist) {{
        dlist.innerHTML = '';
        Array.from(geneSet).sort().forEach(g => {{
          const opt = document.createElement('option');
          opt.value = g;
          dlist.appendChild(opt);
        }});
      }}

      if (pillsContainer) {{
        pillsContainer.innerHTML = '';
        const orangeGenes = (COHORT_DATA.all_orange_driver_genes && COHORT_DATA.all_orange_driver_genes.length > 0)
          ? COHORT_DATA.all_orange_driver_genes
          : Array.from(geneSet).slice(0, 20);

        orangeGenes.forEach(g => {{
          const pill = document.createElement('button');
          pill.className = 'badge';
          pill.style.cssText = 'background: #161b22; border: 1px solid #30363d; color: #58a6ff; cursor: pointer; padding: 2px 7px; border-radius: 12px; font-size: 0.72rem; font-family: monospace; font-weight: 600; transition: all 0.15s ease;';
          pill.textContent = '+ ' + g;
          pill.title = 'Click to add ' + g + ' to Oncoprint';
          pill.onmouseover = () => {{ pill.style.borderColor = '#58a6ff'; pill.style.background = 'rgba(88, 166, 255, 0.15)'; }};
          pill.onmouseout = () => {{ pill.style.borderColor = '#30363d'; pill.style.background = '#161b22'; }};
          pill.onclick = () => {{ addCustomGene(g); }};
          pillsContainer.appendChild(pill);
        }});
      }}
    }}

    // ==============================================================================
    // FINGERPRINTING & IDENTITY VERIFICATION SECTION
    // ==============================================================================
    function renderFingerprintSection() {{
      const fps = COHORT_DATA.fingerprints;
      if (!fps || !fps.patient_summaries) return;

      // Update KPIs
      const kpiRate = document.getElementById('fp-kpi-rate');
      const kpiRateSub = document.getElementById('fp-kpi-rate-sub');
      const kpiSamples = document.getElementById('fp-kpi-samples');
      const kpiIntra = document.getElementById('fp-kpi-intra');
      const kpiInter = document.getElementById('fp-kpi-inter');

      if (kpiRate) {{
        const rate = fps.total_patients > 0 ? (fps.confirmed_patients / fps.total_patients * 100).toFixed(0) : 0;
        kpiRate.textContent = rate + '%';
        if (kpiRateSub) kpiRateSub.textContent = `${{fps.confirmed_patients}} / ${{fps.total_patients}} Patients Verified`;
      }}
      if (kpiSamples) kpiSamples.textContent = fps.total_samples;
      if (kpiIntra) kpiIntra.textContent = (fps.mean_intra_concordance !== undefined && fps.mean_intra_concordance !== null) ? fps.mean_intra_concordance.toFixed(1) + '%' : (fps.mean_intra_dist !== undefined ? ((1 - fps.mean_intra_dist) * 100).toFixed(1) + '%' : '--');
      if (kpiInter) kpiInter.textContent = (fps.mean_inter_concordance !== undefined && fps.mean_inter_concordance !== null) ? fps.mean_inter_concordance.toFixed(1) + '%' : (fps.mean_inter_dist !== undefined ? ((1 - fps.mean_inter_dist) * 100).toFixed(1) + '%' : '--');

      // Render Concordance Table
      const tbody = document.getElementById('fingerprint-tbody');
      if (tbody) {{
        tbody.innerHTML = '';
        Object.values(fps.patient_summaries).forEach(ps => {{
          const isConfirmed = ps.is_confirmed;
          const statusBadge = isConfirmed
            ? '<span class="badge" style="background:#ecfdf5; color:#065f46; border:1px solid #a7f3d0; font-weight:700;">✓ CONFIRMED SAME DONOR</span>'
            : (ps.status === 'SINGLE_SAMPLE_NO_COMPARATOR'
               ? '<span class="badge badge-normal">SINGLE SAMPLE</span>'
               : '<span class="badge" style="background:#fee2e2; color:#991b1b; border:1px solid #fecaca; font-weight:700;">⚠ MISMATCH SUSPECTED</span>');

          const tnConcVal = ps.tumor_vs_normal_concordance;
          const tnShared = ps.tumor_vs_normal_shared;
          const tnConcStr = (tnConcVal !== null && tnConcVal !== undefined)
            ? `<strong style="color:${{tnConcVal >= 85.0 ? '#059669' : '#dc2626'}};">${{tnConcVal.toFixed(1)}}%</strong> <span style="font-size:0.75rem; color:var(--text-muted);">(${{tnShared ? tnShared.toLocaleString() : ''}} shared)</span>`
            : '<span style="color:var(--text-muted)">--</span>';

          const minConcVal = ps.min_intra_concordance;
          const minConcStr = (minConcVal !== null && minConcVal !== undefined)
            ? `<strong style="color:${{minConcVal >= 85.0 ? '#059669' : '#dc2626'}}; font-size:0.95rem;">${{minConcVal.toFixed(1)}}%</strong>`
            : '<span style="color:var(--text-muted)">--</span>';

          const tumName = ps.tumor_clean_name || ps.tumor_id;
          const tumSubText = (ps.tumor_clean_name && ps.tumor_id && ps.tumor_clean_name !== ps.tumor_id)
            ? `<div style="font-size:0.72rem; color:var(--text-muted); font-family:monospace;">${{ps.tumor_id}}</div>`
            : '';
          const tumCell = `<span class="badge badge-normal">${{tumName}}</span>${{tumSubText}}`;

          const plasmas = ps.plasma_comparisons || [];
          if (plasmas.length === 0) {{
            const row = document.createElement('tr');
            row.innerHTML = `
              <td><strong>${{ps.patient_id}}</strong></td>
              <td><span class="badge badge-normal">${{ps.normal_id}}</span></td>
              <td>${{tumCell}}</td>
              <td>${{tnConcStr}}</td>
              <td style="color:var(--text-muted)">No Plasma</td>
              <td style="color:var(--text-muted)">--</td>
              <td style="color:var(--text-muted)">--</td>
              <td>${{minConcStr}}</td>
              <td>${{statusBadge}}</td>
            `;
            tbody.appendChild(row);
          }} else {{
            plasmas.forEach((pl, idx) => {{
              const row = document.createElement('tr');
              const plNConc = pl.vs_normal_concordance;
              const plNShared = pl.vs_normal_shared;
              const plNStr = (plNConc !== null && plNConc !== undefined)
                ? `<strong style="color:${{plNConc >= 85.0 ? '#059669' : '#dc2626'}};">${{plNConc.toFixed(1)}}%</strong> <span style="font-size:0.75rem; color:var(--text-muted);">(${{plNShared ? plNShared.toLocaleString() : ''}} shared)</span>`
                : '<span style="color:var(--text-muted)">--</span>';

              const plTConc = pl.vs_tumor_concordance;
              const plTShared = pl.vs_tumor_shared;
              const plTStr = (plTConc !== null && plTConc !== undefined)
                ? `<strong style="color:${{plTConc >= 85.0 ? '#059669' : '#dc2626'}};">${{plTConc.toFixed(1)}}%</strong> <span style="font-size:0.75rem; color:var(--text-muted);">(${{plTShared ? plTShared.toLocaleString() : ''}} shared)</span>`
                : '<span style="color:var(--text-muted)">--</span>';

              const plName = pl.clean_sample_name || pl.plasma_id;
              const plTpBadge = pl.timepoint ? `<span class="badge badge-primary" style="font-size:0.75rem; font-weight:700;">${{pl.timepoint}}</span>` : '';
              const plDateStr = pl.date ? `<span style="font-size:0.75rem; color:var(--text-muted); margin-left:6px;">(${{pl.date}})</span>` : (pl.days !== null && pl.days !== undefined ? `<span style="font-size:0.75rem; color:var(--text-muted); margin-left:6px;">(Day ${{pl.days}})</span>` : '');
              const plSubText = (pl.clean_sample_name && pl.plasma_id && pl.clean_sample_name !== pl.plasma_id)
                ? `<div style="font-size:0.72rem; color:var(--text-muted); margin-top:2px; font-family:monospace;">${{pl.plasma_id}}</div>`
                : '';

              const plCell = `
                <div style="display:flex; align-items:center; gap:6px; flex-wrap:wrap;">
                  ${{plTpBadge}}
                  <span style="font-weight:600; color:var(--text-primary);">${{plName}}</span>
                  ${{plDateStr}}
                </div>
                ${{plSubText}}
              `;

              row.innerHTML = `
                <td><strong>${{idx === 0 ? ps.patient_id : ''}}</strong></td>
                <td>${{idx === 0 ? `<span class="badge badge-normal">${{ps.normal_id}}</span>` : ''}}</td>
                <td>${{idx === 0 ? tumCell : ''}}</td>
                <td>${{idx === 0 ? tnConcStr : ''}}</td>
                <td>${{plCell}}</td>
                <td>${{plNStr}}</td>
                <td>${{plTStr}}</td>
                <td>${{idx === 0 ? minConcStr : ''}}</td>
                <td>${{idx === 0 ? statusBadge : ''}}</td>
              `;
              tbody.appendChild(row);
            }});
          }}
        }});
      }}

      renderFingerprintMatrix();
    }}

    function renderFingerprintMatrix() {{
      const fps = COHORT_DATA.fingerprints;
      if (!fps || !fps.samples || !fps.matrix) return;

      const filterVal = document.getElementById('fp-matrix-filter') ? document.getElementById('fp-matrix-filter').value : 'all';
      const container = document.getElementById('fp-matrix-container');
      if (!container) return;

      let indices = [];
      fps.samples.forEach((s, idx) => {{
        if (filterVal === 'all') indices.push(idx);
        else if (filterVal === 'tn' && (s.sample_type === 'Normal' || s.sample_type === 'Tumor')) indices.push(idx);
        else if (filterVal === 'plasma' && s.sample_type === 'Plasma') indices.push(idx);
      }});

      if (indices.length === 0) {{
        container.innerHTML = '<p style="padding:1rem; color:var(--text-muted);">No samples match the selected filter.</p>';
        return;
      }}

      let html = '<table class="matrix-table" style="border-collapse:collapse; font-size:0.75rem; font-variant-numeric:tabular-nums; width:100%;">';
      html += '<thead><tr><th style="position:sticky; top:0; left:0; z-index:20; background:#f8fafc; border:1px solid #e2e8f0; padding:6px 8px; font-weight:600; min-width:140px; text-align:left;">Sample</th>';
      indices.forEach(idx => {{
        const s = fps.samples[idx];
        const typeBadge = s.sample_type === 'Normal' ? 'N' : (s.sample_type === 'Tumor' ? 'T' : (s.timepoint || 'P'));
        const color = s.sample_type === 'Normal' ? '#0284c7' : (s.sample_type === 'Tumor' ? '#7c3aed' : '#059669');
        const shortName = s.clean_sample_name || s.sample_id;
        html += `<th style="position:sticky; top:0; z-index:10; background:#f8fafc; border:1px solid #e2e8f0; padding:6px 4px; font-weight:600; text-align:center; min-width:48px; white-space:nowrap;" title="${{s.label}}">
          <div style="font-size:0.65rem; color:${{color}}; font-weight:700;">${{typeBadge}}</div>
          <div style="font-size:0.7rem; font-weight:600; max-width:54px; overflow:hidden; text-overflow:ellipsis;">${{shortName}}</div>
        </th>`;
      }});
      html += '</tr></thead><tbody>';

      indices.forEach(i => {{
        const s_i = fps.samples[i];
        const typeColor = s_i.sample_type === 'Normal' ? '#0284c7' : (s_i.sample_type === 'Tumor' ? '#7c3aed' : '#059669');
        const typeLetter = s_i.sample_type === 'Normal' ? 'N' : (s_i.sample_type === 'Tumor' ? 'T' : (s_i.timepoint || 'P'));
        const rowName = s_i.clean_sample_name || s_i.sample_id;
        html += `<tr><th style="position:sticky; left:0; z-index:5; background:#f8fafc; border:1px solid #e2e8f0; padding:4px 8px; text-align:left; font-weight:600; white-space:nowrap;">
          <span style="color:${{typeColor}}; font-size:0.7rem; font-weight:700; margin-right:4px;">[${{typeLetter}}]</span>
          <span style="color:var(--text-primary); font-size:0.75rem;">${{rowName}}</span>
          <span style="color:var(--text-muted); font-size:0.7rem; margin-left:4px;">${{s_i.patient_id}}</span>
        </th>`;

        indices.forEach(j => {{
          const val = fps.matrix[i][j];
          const s_j = fps.samples[j];
          const isSamePatient = (s_i.patient_id === s_j.patient_id);

          let bg = '#ffffff';
          let textColor = '#64748b';
          let weight = 'normal';
          let border = '1px solid #f1f5f9';

          if (i === j) {{
            bg = '#f1f5f9';
            textColor = '#94a3b8';
            weight = '600';
          }} else if (isSamePatient) {{
            if (val >= 85.0) {{
              bg = '#dcfce7';
              textColor = '#166534';
              weight = '700';
              border = '1px solid #bbf7d0';
            }} else {{
              bg = '#fee2e2';
              textColor = '#991b1b';
              weight = '700';
              border = '1px solid #fecaca';
            }}
          }} else {{
            if (val >= 50.0) {{
              bg = '#fef3c7';
              textColor = '#92400e';
              weight = '600';
            }} else {{
              bg = '#ffffff';
              textColor = '#94a3b8';
            }}
          }}

          const cellText = (i === j) ? '100%' : (val !== null && val !== undefined ? val.toFixed(1) + '%' : '--');
          const tooltip = `${{s_i.label}} vs ${{s_j.label}}\\nConcordance: ${{cellText}}\\nStatus: ${{isSamePatient ? (val >= 85.0 ? 'Intra-Patient Match (Confirmed Same Donor)' : 'MISMATCH WARNING') : 'Unrelated Donors'}}`;

          html += `<td style="background:${{bg}}; color:${{textColor}}; font-weight:${{weight}}; border:${{border}}; text-align:center; padding:4px 2px; cursor:default;" title="${{tooltip}}">
            ${{cellText}}
          </td>`;
        }});
        html += '</tr>';
      }});

      html += '</tbody></table>';
      container.innerHTML = html;
    }}

    function exportFingerprintMatrixCSV() {{
      const fps = COHORT_DATA.fingerprints;
      if (!fps || !fps.samples || !fps.matrix) return;
      let csvContent = 'data:text/csv;charset=utf-8,';
      const header = ['Sample', ...fps.samples.map(s => s.label)];
      csvContent += header.map(h => `"${{h}}"`).join(',') + String.fromCharCode(10);
      fps.matrix.forEach((row, i) => {{
        const r = [fps.samples[i].label, ...row.map(v => (v !== null && v !== undefined ? v.toFixed(2) + '%' : 'NA'))];
        csvContent += r.map(c => `"${{c}}"`).join(',') + String.fromCharCode(10);
      }});
      const encodedUri = encodeURI(csvContent);
      const link = document.createElement('a');
      link.setAttribute('href', encodedUri);
      link.setAttribute('download', 'cohort_fingerprint_concordance_matrix.csv');
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }}

    function filterFingerprintTable() {{
      const query = (document.getElementById('fp-search').value || '').toLowerCase().trim();
      const rows = document.querySelectorAll('#fingerprint-tbody tr');
      rows.forEach(r => {{
        const text = r.innerText.toLowerCase();
        r.style.display = text.includes(query) ? '' : 'none';
      }});
    }}

    // Modal keydown & click-outside handlers
    window.addEventListener('keydown', (e) => {{
      if (e.key === 'Escape') {{
        closeOrangePdfModal();
        closeIchorPdfModal();
        closePatientTrajectoryModal();
        closePlasmaModalityModal();
      }}
    }});

    // Init
    window.addEventListener('DOMContentLoaded', () => {{
      populateOverview();
      initGeneAddControls();
      renderFingerprintSection();

      const orangeModalElem = document.getElementById('orange-pdf-modal');
      if (orangeModalElem) {{
        orangeModalElem.addEventListener('click', (e) => {{
          if (e.target === orangeModalElem) closeOrangePdfModal();
        }});
      }}

      const ichorModalElem = document.getElementById('ichor-pdf-modal');
      if (ichorModalElem) {{
        ichorModalElem.addEventListener('click', (e) => {{
          if (e.target === ichorModalElem) closeIchorPdfModal();
        }});
      }}
    }});
  </script>
</body>
</html>
"""

    with open(html_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    # Symlink or duplicate to index.html and cohort_wisp_dashboard.html
    try:
        index_file = out_path / "index.html"
        with open(index_file, "w", encoding="utf-8") as f:
            f.write(html_content)
        compat_dashboard = out_path / "cohort_wisp_dashboard.html"
        with open(compat_dashboard, "w", encoding="utf-8") as f:
            f.write(html_content)
    except Exception:
        pass

    html_size_mb = html_file.stat().st_size / (1024 * 1024) if html_file.exists() else 0.0
    print(f"Interactive dashboard generated: {html_file} ({html_size_mb:.2f} MB)")


# ==============================================================================
# 8. MAIN ENTRYPOINT
# ==============================================================================

def main():
    args = parse_args()
    results_dir = os.path.abspath(args.results_dir)

    if not os.path.isdir(results_dir):
        print(f"Error: Results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = args.output_dir or os.path.join(results_dir, "wgts_summary")
    genes_of_interest = parse_genes_arg(args.genes_of_interest)

    print(f"==================================================================")
    print(f"  compile_wgts_summary.py - WGTS Multi-Omics & MRD Reporting")
    print(f"  Results Directory:   {results_dir}")
    print(f"  Output Directory:    {output_dir}")
    print(f"  Prioritized Genes:   {len(genes_of_interest)} genes")
    print(f"==================================================================")

    # 1. Compile Cohort
    cohort = compile_cohort(
        results_dir=results_dir,
        genes_of_interest=genes_of_interest,
        sample_sheets_dir=args.sample_sheets_dir,
        metadata_tsv=args.metadata_tsv,
        timepoints_tsv=args.timepoints_tsv,
        purity_overrides_tsv=args.purity_overrides,
        output_dir=output_dir,
    )
    if args.cohort_name:
        cohort["cohort_name"] = args.cohort_name

    print(f"Discovered {cohort['total_analyzed_patients']} patient(s), "
          f"{len(cohort['all_tumors'])} primary tumor(s), "
          f"{len(cohort['all_plasmas'])} plasma ctDNA sample(s).")

    # 2. Build Oncoprint Matrix
    oncoprint_data = generate_oncoprint_matrix(cohort, genes_of_interest, top_n_cohort=20)

    # Determine effective thread/worker count
    threads = getattr(args, "threads", 8)
    if getattr(args, "fingerprint_threads", None) is not None:
        threads = args.fingerprint_threads

    # 3. Evaluate Sample Fingerprinting & Identity Confirmation
    if getattr(args, "skip_fingerprint", False):
        print("Skipping sample fingerprinting (--skip-fingerprint specified).")
        cohort["fingerprints"] = {
            "samples": [], "matrix": [], "distance_matrix": [], "pairwise_list": [],
            "patient_summaries": {}, "total_samples": 0, "total_patients": 0,
            "confirmed_patients": 0, "mismatched_patients": 0,
            "mean_intra_concordance": 0.0, "mean_inter_concordance": 0.0,
            "mean_intra_dist": 0.0, "mean_inter_dist": 0.0
        }
    else:
        print("Evaluating sample fingerprinting and donor concordance (Amber BAF)...")
        cohort["fingerprints"] = evaluate_cohort_fingerprints(
            cohort,
            output_dir=output_dir,
            stride=getattr(args, "fingerprint_stride", 20),
            max_workers=threads
        )
    if cohort["fingerprints"]["total_samples"] > 0:
        print(f"Fingerprinted {cohort['fingerprints']['total_samples']} sample(s): "
              f"{cohort['fingerprints']['confirmed_patients']}/{cohort['fingerprints']['total_patients']} "
              f"patients confirmed same individual (mean intra-conc: {cohort['fingerprints']['mean_intra_concordance']:.1f}%, "
              f"mean inter-conc: {cohort['fingerprints']['mean_inter_concordance']:.1f}%).")

    # 4. Write Master TSVs
    write_tsv_reports(cohort, oncoprint_data, output_dir)

    embed_all = getattr(args, 'embed_all_plots', False)

    # 5. Generate Static Figures (parallelized across worker processes)
    generate_static_plots(
        cohort,
        oncoprint_data,
        output_dir,
        formats=args.format,
        skip_plasma_plots=getattr(args, 'skip_plasma_plots', False),
        threads=threads,
        embed_all_plots=embed_all
    )

    # 6. Generate Standalone HTML Interactive Dashboard
    generate_interactive_dashboard(cohort, oncoprint_data, output_dir, embed_all_plots=embed_all)

    print(f"==================================================================")
    print(f"  Successfully compiled WGTS summary & dashboard in:")
    print(f"  {output_dir}")
    print(f"==================================================================")


if __name__ == "__main__":
    main()
