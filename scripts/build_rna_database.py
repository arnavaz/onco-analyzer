#!/usr/bin/env python3
"""
Pan-Cancer RNA-Seq Ingestion & Cross-Cohort Normalization Engine
===============================================================
Parses multi-cohort HMF Isofox RNA-Seq outputs (gene-level expression,
transcript isoforms, novel splice junctions, and chimeric fusions)
alongside cBioPortal / CanDIG clinical metadata into an analytical
database (DuckDB or SQLite) with global TMM normalization and
empirical Bayes (ComBat) batch effect adjustment.

Usage:
  python3 scripts/build_rna_database.py \
      --results-dir /path/to/oncoanalyser_outputs \
      --cbio-dir /path/to/cbio_import_dir \
      --output-db hmf_rnaseq_cohort.duckdb \
      --gene-distribution /path/to/isofox.hmf_38_151_2600.gene_distribution.csv.gz
"""

import os
import sys
import re
import csv
import gzip
import json
import math
import bisect
import argparse
import datetime
import gc
import concurrent.futures
from array import array
from pathlib import Path
from collections import defaultdict
from contextlib import contextmanager

# Numerical packages
try:
    import numpy as np
except ImportError:
    print("Error: numpy is required for normalization. Please install via: pip install numpy", file=sys.stderr)
    sys.exit(1)

# Check for DuckDB or SQLite backend
HAS_DUCKDB = False
try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False

import sqlite3

# Known candidate locations for HMF gene distribution reference
DEFAULT_GENE_DIST_CANDIDATES = [
    Path("/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0/hmf_pipeline_resources.38_v3.1.0--1/rna/isofox.hmf_38_151_2600.gene_distribution.csv.gz"),
    Path("/Users/jbruce/Desktop/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0/hmf_pipeline_resources.38_v3.1.0--1/rna/isofox.hmf_38_151_2600.gene_distribution.csv.gz"),
]

# Canonical mapping from clinical diagnosis / cancer type strings to HMF CancerType categories
HMF_CANCER_TYPE_MAP = [
    (re.compile(r"glioblastoma|astrocytoma|glioma|cns|brain|medulloblastoma", re.I), "Nervous system"),
    (re.compile(r"colorectal|colon|rectum|bowel|gi|gastrointestinal", re.I), "Colorectum"),
    (re.compile(r"pancreas|pancreatic", re.I), "Pancreas"),
    (re.compile(r"breast", re.I), "Breast"),
    (re.compile(r"lung|non-small|small cell", re.I), "Lung"),
    (re.compile(r"prostate", re.I), "Prostate"),
    (re.compile(r"ovary|ovarian", re.I), "Ovary"),
    (re.compile(r"urothelial|bladder|urinary|kidney|renal", re.I), "Urothelial tract"),
    (re.compile(r"melanoma|skin", re.I), "Skin"),
    (re.compile(r"head and neck|larynx|oral|pharynx", re.I), "Head and neck"),
    (re.compile(r"esophagus|stomach|gastric", re.I), "Esophagus"),
    (re.compile(r"liver|hepatocellular|biliary|cholangio", re.I), "Biliary"),
    (re.compile(r"sarcoma|bone|soft tissue", re.I), "Bone/Soft tissue"),
    (re.compile(r"uterus|uterine|endometrial|cervix", re.I), "Uterus"),
]


def resolve_gene_dist_file(custom_path=None):
    """Find the Isofox gene distribution reference file."""
    if custom_path and Path(custom_path).is_file():
        return Path(custom_path)
    for cand in DEFAULT_GENE_DIST_CANDIDATES:
        if cand.is_file():
            return cand
    return None


def map_to_hmf_cancer_type(cancer_type_str, diagnosis_str=""):
    """Map any cancer type or diagnosis string to an HMF baseline cancer type category."""
    combined = f"{cancer_type_str or ''} {diagnosis_str or ''}".strip()
    for pattern, hmf_cat in HMF_CANCER_TYPE_MAP:
        if pattern.search(combined):
            return hmf_cat
    return "ALL"


# ---------------------------------------------------------------------------
# File Utilities & Smart Text Opener
# ---------------------------------------------------------------------------

@contextmanager
def open_smart_text(filepath):
    """
    Open a file as text, automatically detecting if it is gzip-compressed
    (via magic bytes b'\\x1f\\x8b') regardless of file extension, and using
    utf-8 with replacement for non-decodable bytes.
    Yields: (file_handle, detected_delimiter)
    """
    with open(filepath, "rb") as f:
        magic = f.read(2)

    is_gz = (magic == b"\x1f\x8b")
    if is_gz:
        handle = gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    else:
        handle = open(filepath, "r", encoding="utf-8", errors="replace")

    try:
        first_line = handle.readline()
        sep = "\t" if "\t" in first_line else ","
        handle.seek(0)
        yield handle, sep
    finally:
        handle.close()


def filter_isofox_files(file_list):
    """Filter out hidden/AppleDouble (._*) and companion index files (.tbi, .idx, .csi, etc.)."""
    valid = []
    seen = set()
    for f in file_list:
        if f.name.startswith((".", "._")):
            continue
        name_lower = f.name.lower()
        if any(name_lower.endswith(ext) for ext in [".tbi", ".idx", ".csi", ".tmp", ".md5", ".sha256", ".bai"]):
            continue
        if str(f) not in seen:
            seen.add(str(f))
            valid.append(f)
    return valid


# ---------------------------------------------------------------------------
# Isofox Output Parsers
# ---------------------------------------------------------------------------

def find_sample_dirs(results_dirs):
    """
    Locate all sample directories within one or more results directories.
    A valid sample directory contains an 'isofox' subdirectory or *.isf.gene_data.* files.
    """
    sample_dirs = []
    seen_samples = set()

    for r_dir in results_dirs:
        r_path = Path(r_dir).resolve()
        if not r_path.is_dir():
            print(f"Warning: Directory not found: {r_path}", file=sys.stderr)
            continue

        # Check direct subdirectories
        for child in sorted(r_path.iterdir()):
            if not child.is_dir() or child.name.startswith((".", "work", "pipeline_info")):
                continue

            # Check if this child is a sample dir
            is_sample = False
            if (child / "isofox").is_dir():
                is_sample = True
            elif list(child.glob("*.isf.gene_data.*")) or list(child.glob("*.gene_data.*")):
                is_sample = True
            elif (child / "purple").is_dir() or (child / "linx").is_dir():
                # Check nested or sibling isofox
                if (child / "isofox").is_dir() or list(child.glob("*.isf.*")):
                    is_sample = True

            if is_sample and child.name not in seen_samples:
                seen_samples.add(child.name)
                sample_dirs.append(child)
            elif child.is_dir():
                # Check grandchild directories (e.g. gcp_download/COHORT/PATIENT/SAMPLE)
                try:
                    for grandchild in sorted(child.iterdir()):
                        if not grandchild.is_dir() or grandchild.name.startswith((".", "work")):
                            continue
                        if (grandchild / "isofox").is_dir() or list(grandchild.glob("*.isf.gene_data.*")):
                            if grandchild.name not in seen_samples:
                                seen_samples.add(grandchild.name)
                                sample_dirs.append(grandchild)
                except OSError:
                    pass

    return sample_dirs


def parse_isofox_gene_data(sample_dir):
    """
    Extract gene-level expression data from *.isf.gene_data.* in sample_dir.
    Returns: list of dicts with gene metrics.
    """
    isofox_dir = sample_dir / "isofox" if (sample_dir / "isofox").is_dir() else sample_dir
    raw_files = list(isofox_dir.glob("*.isf.gene_data.*")) + list(isofox_dir.glob("*.gene_data.*"))
    files = filter_isofox_files(raw_files)
    if not files:
        return []

    target = sorted(files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]

    genes = []
    try:
        with open_smart_text(target) as (f, sep):
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                gene_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                gene_id = (row.get("GeneId") or row.get("gene_id") or "").strip()
                if not gene_name and not gene_id:
                    continue
                if not gene_name:
                    gene_name = gene_id

                adj_tpm = _safe_float(row.get("AdjTPM") or row.get("AdjustedTpm") or row.get("TPM") or row.get("Tpm") or row.get("RawTPM"), 0.0)
                raw_tpm = _safe_float(row.get("RawTPM") or row.get("TPM") or row.get("Tpm") or row.get("AdjTPM"), 0.0)
                spliced_frags = _safe_float(row.get("SplicedFragments") or row.get("spliced_fragments") or row.get("FittedFragments"), 0.0)
                unspliced_frags = _safe_float(row.get("UnsplicedFragments") or row.get("unspliced_fragments"), 0.0)
                total_fitted = spliced_frags + unspliced_frags

                genes.append({
                    "gene_name": gene_name,
                    "ensembl_gene_id": gene_id,
                    "chromosome": (row.get("Chromosome") or row.get("chromosome") or "").strip(),
                    "gene_length": _safe_int(row.get("GeneLength") or row.get("gene_length"), 0),
                    "spliced_fragments": spliced_frags,
                    "unspliced_fragments": unspliced_frags,
                    "fitted_fragments": total_fitted,
                    "adj_tpm": adj_tpm,
                    "raw_tpm": raw_tpm,
                    "log2_adj_tpm": round(math.log2(adj_tpm + 1.0), 4),
                    "percentile_cancer": _safe_float(row.get("PercentileCancer"), -1.0),
                    "median_tpm_cancer": _safe_float(row.get("MedianTpmCancer"), -1.0),
                    "percentile_cohort": _safe_float(row.get("PercentileCohort"), -1.0),
                    "median_tpm_cohort": _safe_float(row.get("MedianTpmCohort"), -1.0),
                })
    except Exception as e:
        print(f"Warning: Failed reading gene data in {sample_dir}: {e}", file=sys.stderr)

    return genes


def parse_isofox_trans_data(sample_dir):
    """
    Extract transcript isoform expression data from *.isf.trans_data.* in sample_dir.
    Returns: list of dicts with transcript metrics.
    """
    isofox_dir = sample_dir / "isofox" if (sample_dir / "isofox").is_dir() else sample_dir
    raw_files = list(isofox_dir.glob("*.isf.trans_data.*")) + list(isofox_dir.glob("*.trans_data.*"))
    files = filter_isofox_files(raw_files)
    if not files:
        return []

    target = sorted(files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]

    transcripts = []
    try:
        with open_smart_text(target) as (f, sep):
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                gene_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                trans_id = (row.get("TranscriptId") or row.get("transcript_id") or "").strip()
                if not trans_id:
                    continue

                canonical_val = str(row.get("Canonical") or row.get("is_canonical") or "").strip().lower()
                is_canonical = 1 if canonical_val in ["true", "1", "t", "yes"] else 0

                tpm = _safe_float(row.get("TPM") or row.get("tpm") or row.get("AdjTPM"), 0.0)
                fitted_frags = _safe_float(row.get("FittedFragments") or row.get("fitted_fragments"), 0.0)

                transcripts.append({
                    "gene_name": gene_name,
                    "transcript_id": trans_id,
                    "transcript_name": (row.get("TranscriptName") or row.get("transcript_name") or trans_id).strip(),
                    "is_canonical": is_canonical,
                    "exon_count": _safe_int(row.get("ExonCount") or row.get("exon_count"), 0),
                    "transcript_length": _safe_int(row.get("TranscriptLength") or row.get("transcript_length"), 0),
                    "effective_length": _safe_float(row.get("EffectiveLength") or row.get("effective_length"), 0.0),
                    "fitted_fragments": fitted_frags,
                    "tpm": tpm,
                    "log2_tpm": round(math.log2(tpm + 1.0), 4),
                })
    except Exception as e:
        print(f"Warning: Failed reading transcript data in {sample_dir}: {e}", file=sys.stderr)

    return transcripts


def parse_isofox_alt_splice(sample_dir):
    """
    Extract novel splice junctions from *.isf.alt_splice_junc.* or *.alt_splice_junctions.*.
    """
    isofox_dir = sample_dir / "isofox" if (sample_dir / "isofox").is_dir() else sample_dir
    raw_files = (
        list(isofox_dir.glob("*.isf.alt_splice_junc.*")) +
        list(isofox_dir.glob("*.alt_splice_junctions.*")) +
        list(isofox_dir.glob("*.isf.alt_splice_junctions.*"))
    )
    files = filter_isofox_files(raw_files)
    if not files:
        return []

    target = sorted(files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]

    junctions = []
    try:
        with open_smart_text(target) as (f, sep):
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                gene_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                chrom = (row.get("Chromosome") or row.get("chromosome") or "").strip()
                start = _safe_int(row.get("Start") or row.get("start") or row.get("JunctionStart"), 0)
                end = _safe_int(row.get("End") or row.get("end") or row.get("JunctionEnd"), 0)
                if not chrom or start <= 0 or end <= 0:
                    continue

                junctions.append({
                    "gene_name": gene_name,
                    "chromosome": chrom,
                    "start_pos": start,
                    "end_pos": end,
                    "strand": (row.get("Strand") or row.get("strand") or "+").strip(),
                    "junction_type": (row.get("Type") or row.get("type") or row.get("JunctionType") or "NOVEL_SPLICE").strip(),
                    "fragment_support": _safe_int(row.get("FragmentCount") or row.get("fragment_count") or row.get("SupportCount"), 0),
                    "cohort_frequency": _safe_int(row.get("CohortFrequency") or row.get("cohort_frequency"), 0),
                })
    except Exception as e:
        print(f"Warning: Failed reading alternate splice junctions in {sample_dir}: {e}", file=sys.stderr)

    return junctions


def parse_isofox_fusions(sample_dir):
    """
    Extract chimeric RNA fusions from *.isf.fusions.* or *.pass_fusions.*.
    """
    isofox_dir = sample_dir / "isofox" if (sample_dir / "isofox").is_dir() else sample_dir
    raw_files = (
        list(isofox_dir.glob("*.pass_fusions.*")) +
        list(isofox_dir.glob("*.isf.pass_fusions.*")) +
        list(isofox_dir.glob("*.fusions.*")) +
        list(isofox_dir.glob("*.isf.fusions.*"))
    )
    files = filter_isofox_files(raw_files)
    if not files:
        return []

    target = sorted(files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]

    fusions = []
    try:
        with open_smart_text(target) as (f, sep):
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                g_up = (row.get("GeneUp") or row.get("gene_up") or row.get("FivePrimeGene") or "").strip()
                g_down = (row.get("GeneDown") or row.get("gene_down") or row.get("ThreePrimeGene") or "").strip()
                if not g_up and not g_down:
                    continue

                in_frame_val = str(row.get("InFrame") or row.get("in_frame") or "").strip().lower()
                in_frame = 1 if in_frame_val in ["true", "1", "yes"] else 0

                split_cnt = _safe_int(row.get("SplitFragments") or row.get("split_fragments") or row.get("SplitReads"), 0)
                realigned_cnt = _safe_int(row.get("RealignedFragments") or row.get("realigned_fragments"), 0)
                discordant_cnt = _safe_int(row.get("DiscordantFragments") or row.get("discordant_fragments"), 0)
                total_support = split_cnt + realigned_cnt + discordant_cnt

                fusions.append({
                    "upstream_gene": g_up,
                    "downstream_gene": g_down,
                    "in_frame": in_frame,
                    "split_reads": split_cnt + realigned_cnt,
                    "discordant_reads": discordant_cnt,
                    "total_support": total_support,
                    "fusion_type": (row.get("Type") or row.get("type") or "FUSION").strip(),
                })
    except Exception as e:
        print(f"Warning: Failed reading fusions in {sample_dir}: {e}", file=sys.stderr)

    return fusions


# ---------------------------------------------------------------------------
# Parallel Workers (Pass 1 & Pass 2 Streaming)
# ---------------------------------------------------------------------------

def extract_sample_counts_worker(s_dir_str):
    """
    Pass 1 Worker: Extracts only (sample_id, {gene_name: fitted_fragments})
    to construct global count matrix for TMM/ComBat with near-zero memory footprint.
    """
    s_dir = Path(s_dir_str)
    s_id = s_dir.name
    isofox_dir = s_dir / "isofox" if (s_dir / "isofox").is_dir() else s_dir
    raw_files = list(isofox_dir.glob("*.isf.gene_data.*")) + list(isofox_dir.glob("*.gene_data.*"))
    files = filter_isofox_files(raw_files)
    if not files:
        return s_id, None

    target = sorted(files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]
    counts = {}
    try:
        with open_smart_text(target) as (f, sep):
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                g_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                if not g_name:
                    g_name = (row.get("GeneId") or row.get("gene_id") or "").strip()
                if not g_name:
                    continue
                spliced = _safe_float(row.get("SplicedFragments") or row.get("spliced_fragments") or row.get("FittedFragments"), 0.0)
                unspliced = _safe_float(row.get("UnsplicedFragments") or row.get("unspliced_fragments"), 0.0)
                counts[g_name] = spliced + unspliced
    except Exception as e:
        print(f"Warning: Failed reading gene counts for {s_id}: {e}", file=sys.stderr)
        return s_id, None

    return s_id, counts


def parse_and_format_sample_worker(args):
    """
    Pass 2 Worker: Streams all records (genes, transcripts, splice junctions, fusions)
    for a single sample and pre-formats them into database insertion tuples.
    args: (sample_dir_str, norm_factor, sample_batch_adj_dict)
    """
    s_dir_str, norm_factor, sample_batch_adj = args
    s_dir = Path(s_dir_str)
    s_id = s_dir.name

    gene_rows = []
    trans_rows = []
    splice_rows = []
    fusion_rows = []

    # 1. Gene expression
    isofox_dir = s_dir / "isofox" if (s_dir / "isofox").is_dir() else s_dir
    raw_gene_files = list(isofox_dir.glob("*.isf.gene_data.*")) + list(isofox_dir.glob("*.gene_data.*"))
    gene_files = filter_isofox_files(raw_gene_files)
    if gene_files:
        target = sorted(gene_files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]
        try:
            with open_smart_text(target) as (f, sep):
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    gene_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                    gene_id = (row.get("GeneId") or row.get("gene_id") or "").strip()
                    if not gene_name and not gene_id:
                        continue
                    if not gene_name:
                        gene_name = gene_id

                    adj_tpm = _safe_float(row.get("AdjTPM") or row.get("AdjustedTpm") or row.get("TPM") or row.get("Tpm") or row.get("RawTPM"), 0.0)
                    raw_tpm = _safe_float(row.get("RawTPM") or row.get("TPM") or row.get("Tpm") or row.get("AdjTPM"), 0.0)
                    spliced_frags = _safe_float(row.get("SplicedFragments") or row.get("spliced_fragments") or row.get("FittedFragments"), 0.0)
                    unspliced_frags = _safe_float(row.get("UnsplicedFragments") or row.get("unspliced_fragments"), 0.0)
                    total_fitted = spliced_frags + unspliced_frags

                    norm_cnt = total_fitted * norm_factor
                    combat_expr = sample_batch_adj.get(gene_name, 0.0) if sample_batch_adj is not None else norm_cnt

                    gene_rows.append((
                        s_id,
                        gene_name,
                        gene_id,
                        (row.get("Chromosome") or row.get("chromosome") or "").strip(),
                        _safe_int(row.get("GeneLength") or row.get("gene_length"), 0),
                        total_fitted,
                        round(norm_cnt, 2),
                        round(raw_tpm, 4),
                        round(adj_tpm, 4),
                        round(math.log2(adj_tpm + 1.0), 4),
                        round(math.log2(norm_cnt + 1.0), 4),
                        round(combat_expr, 4),
                        _safe_float(row.get("PercentileCancer"), -1.0),
                        _safe_float(row.get("MedianTpmCancer"), -1.0),
                        _safe_float(row.get("PercentileCohort"), -1.0),
                        _safe_float(row.get("MedianTpmCohort"), -1.0),
                    ))
        except Exception as e:
            print(f"Warning: Failed reading gene data in {s_dir}: {e}", file=sys.stderr)

    # 2. Transcript Isoforms
    raw_trans_files = list(isofox_dir.glob("*.isf.trans_data.*")) + list(isofox_dir.glob("*.trans_data.*"))
    trans_files = filter_isofox_files(raw_trans_files)
    if trans_files:
        target = sorted(trans_files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]
        try:
            raw_transcripts = []
            gene_total_tpm = defaultdict(float)
            with open_smart_text(target) as (f, sep):
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    trans_id = (row.get("TranscriptId") or row.get("transcript_id") or "").strip()
                    if not trans_id:
                        continue
                    gene_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                    tpm = _safe_float(row.get("TPM") or row.get("tpm") or row.get("AdjTPM"), 0.0)
                    gene_total_tpm[gene_name] += tpm

                    canonical_val = str(row.get("Canonical") or row.get("is_canonical") or "").strip().lower()
                    is_canonical = 1 if canonical_val in ["true", "1", "t", "yes"] else 0
                    raw_transcripts.append((
                        gene_name,
                        trans_id,
                        (row.get("TranscriptName") or row.get("transcript_name") or trans_id).strip(),
                        is_canonical,
                        _safe_int(row.get("ExonCount") or row.get("exon_count"), 0),
                        _safe_int(row.get("TranscriptLength") or row.get("transcript_length"), 0),
                        _safe_float(row.get("EffectiveLength") or row.get("effective_length"), 0.0),
                        _safe_float(row.get("FittedFragments") or row.get("fitted_fragments"), 0.0),
                        tpm
                    ))

            for (g_name, tr_id, tr_name, is_can, ex_cnt, tr_len, eff_len, fit_fr, tpm) in raw_transcripts:
                tot_tpm = gene_total_tpm[g_name]
                iso_frac = round(tpm / tot_tpm, 4) if tot_tpm > 0 else 0.0
                trans_rows.append((
                    s_id, g_name, tr_id, tr_name, is_can, ex_cnt, tr_len, eff_len, fit_fr,
                    round(tpm, 4), round(math.log2(tpm + 1.0), 4), iso_frac
                ))
        except Exception as e:
            print(f"Warning: Failed reading transcript data in {s_dir}: {e}", file=sys.stderr)

    # 3. Alternate Splice Junctions
    raw_splice_files = (
        list(isofox_dir.glob("*.isf.alt_splice_junc.*")) +
        list(isofox_dir.glob("*.alt_splice_junctions.*")) +
        list(isofox_dir.glob("*.isf.alt_splice_junctions.*"))
    )
    splice_files = filter_isofox_files(raw_splice_files)
    if splice_files:
        target = sorted(splice_files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]
        try:
            with open_smart_text(target) as (f, sep):
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    gene_name = (row.get("GeneName") or row.get("gene_name") or "").strip()
                    chrom = (row.get("Chromosome") or row.get("chromosome") or "").strip()
                    start = _safe_int(row.get("Start") or row.get("start") or row.get("JunctionStart"), 0)
                    end = _safe_int(row.get("End") or row.get("end") or row.get("JunctionEnd"), 0)
                    if not chrom or start <= 0 or end <= 0:
                        continue
                    splice_rows.append((
                        s_id,
                        gene_name,
                        chrom,
                        start,
                        end,
                        (row.get("Strand") or row.get("strand") or "+").strip(),
                        (row.get("Type") or row.get("type") or row.get("JunctionType") or "NOVEL_SPLICE").strip(),
                        _safe_int(row.get("FragmentCount") or row.get("fragment_count") or row.get("SupportCount"), 0),
                        _safe_int(row.get("CohortFrequency") or row.get("cohort_frequency"), 0)
                    ))
        except Exception as e:
            print(f"Warning: Failed reading splice junctions in {s_dir}: {e}", file=sys.stderr)

    # 4. Fusions
    raw_fusion_files = (
        list(isofox_dir.glob("*.pass_fusions.*")) +
        list(isofox_dir.glob("*.isf.pass_fusions.*")) +
        list(isofox_dir.glob("*.fusions.*")) +
        list(isofox_dir.glob("*.isf.fusions.*"))
    )
    fusion_files = filter_isofox_files(raw_fusion_files)
    if fusion_files:
        target = sorted(fusion_files, key=lambda x: (not any(x.name.endswith(e) for e in ['.csv', '.tsv', '.csv.gz', '.tsv.gz']), x.name))[0]
        try:
            with open_smart_text(target) as (f, sep):
                reader = csv.DictReader(f, delimiter=sep)
                for row in reader:
                    g_up = (row.get("GeneUp") or row.get("gene_up") or row.get("FivePrimeGene") or "").strip()
                    g_down = (row.get("GeneDown") or row.get("gene_down") or row.get("ThreePrimeGene") or "").strip()
                    if not g_up and not g_down:
                        continue
                    in_frame_val = str(row.get("InFrame") or row.get("in_frame") or "").strip().lower()
                    in_frame = 1 if in_frame_val in ["true", "1", "yes"] else 0
                    split_cnt = _safe_int(row.get("SplitFragments") or row.get("split_fragments") or row.get("SplitReads"), 0)
                    realigned_cnt = _safe_int(row.get("RealignedFragments") or row.get("realigned_fragments"), 0)
                    discordant_cnt = _safe_int(row.get("DiscordantFragments") or row.get("discordant_fragments"), 0)
                    fusion_rows.append((
                        s_id,
                        g_up,
                        g_down,
                        in_frame,
                        split_cnt + realigned_cnt,
                        discordant_cnt,
                        split_cnt + realigned_cnt + discordant_cnt,
                        (row.get("Type") or row.get("type") or "FUSION").strip()
                    ))
        except Exception as e:
            print(f"Warning: Failed reading fusions in {s_dir}: {e}", file=sys.stderr)

    return s_id, gene_rows, trans_rows, splice_rows, fusion_rows


# ---------------------------------------------------------------------------
# Clinical Metadata Parser (cBioPortal & CanDIG Integration)
# ---------------------------------------------------------------------------

def parse_clinical_metadata(cbio_dirs=None, metadata_tsv=None):
    """
    Parse cBioPortal clinical sample and patient tables across one or more cohort directories,
    or a flat metadata TSV.
    Returns: dict mapping sample_id -> dict of clinical attributes.
    """
    clinical_data = defaultdict(dict)

    # 1. Read cBioPortal data_clinical_sample.txt and data_clinical_patient.txt across all dirs
    if cbio_dirs:
        if isinstance(cbio_dirs, (str, Path)):
            cbio_dirs = [cbio_dirs]

        for cbio_dir in cbio_dirs:
            c_dir = Path(cbio_dir).resolve()
            if not c_dir.exists():
                print(f"Warning: Clinical directory not found: {c_dir}", file=sys.stderr)
                continue

            # Look for clinical files directly (fast path), falling back to glob only if needed
            candidates = []
            if (c_dir / "import_folder" / "data_clinical_sample.txt").is_file():
                candidates.append(c_dir / "import_folder" / "data_clinical_sample.txt")
            elif (c_dir / "data_clinical_sample.txt").is_file():
                candidates.append(c_dir / "data_clinical_sample.txt")
            else:
                candidates = [f for f in c_dir.glob("**/data_clinical_sample.txt") if not f.name.startswith((".", "._"))]
            sample_files = [f for f in candidates if not f.name.startswith((".", "._"))]
            
            for sample_file in sample_files:
                patient_file = sample_file.parent / "data_clinical_patient.txt"
                patient_map = {}
                if patient_file.is_file():
                    try:
                        with open(patient_file, "r", encoding="utf-8", errors="replace") as f:
                            headers = None
                            for line in f:
                                if line.startswith("#") or not line.strip():
                                    continue
                                parts = [p.strip('"') for p in line.rstrip("\r\n").split("\t")]
                                if headers is None:
                                    headers = parts
                                    continue
                                if len(parts) > 0:
                                    p_id = parts[0]
                                    patient_map[p_id] = {headers[i]: parts[i] if i < len(parts) else "NA" for i in range(len(headers))}
                    except Exception as e:
                        print(f"Warning: Failed reading patient clinical file {patient_file}: {e}", file=sys.stderr)

                if sample_file.is_file():
                    try:
                        with open(sample_file, "r", encoding="utf-8", errors="replace") as f:
                            headers = None
                            for line in f:
                                if line.startswith("#") or not line.strip():
                                    continue
                                parts = [p.strip('"') for p in line.rstrip("\r\n").split("\t")]
                                if headers is None:
                                    headers = parts
                                    continue
                                if len(parts) > 1:
                                    s_idx = headers.index("SAMPLE_ID") if "SAMPLE_ID" in headers else (1 if len(parts) > 1 else 0)
                                    s_id = parts[s_idx]
                                    row_dict = {headers[i]: parts[i] if i < len(parts) else "NA" for i in range(len(headers))}
                                    
                                    p_id = row_dict.get("PATIENT_ID") or parts[0]
                                    if p_id in patient_map:
                                        for pk, pv in patient_map[p_id].items():
                                            if pk not in row_dict:
                                                row_dict[pk] = pv

                                    # Index by BOTH sample ID and patient ID so folders named by patient match!
                                    clinical_data[s_id].update(row_dict)
                                    if p_id:
                                        clinical_data[p_id].update(row_dict)
                    except Exception as e:
                        print(f"Warning: Failed reading sample clinical file {sample_file}: {e}", file=sys.stderr)

    # 2. Read explicit metadata TSV if provided
    if metadata_tsv:
        m_path = Path(metadata_tsv)
        if m_path.is_file():
            try:
                sep = "\t" if m_path.suffix in [".tsv", ".txt"] else ","
                with open(m_path, "r", encoding="utf-8", errors="replace") as f:
                    reader = csv.DictReader(f, delimiter=sep)
                    for row in reader:
                        s_id = row.get("SAMPLE_ID") or row.get("SampleId") or row.get("sample_id") or row.get("Sample_ID") or row.get("ID")
                        if not s_id:
                            continue
                        s_id = s_id.strip()
                        clinical_data[s_id].update(row)
            except Exception as e:
                print(f"Warning: Failed reading metadata file {metadata_tsv}: {e}", file=sys.stderr)

    return clinical_data


# ---------------------------------------------------------------------------
# Cross-Cohort Normalization Engine (TMM + ComBat)
# ---------------------------------------------------------------------------

def compute_tmm_factors(counts_matrix):
    """
    Compute Robinson & Oshlack (2010) Trimmed Mean of M-values (TMM) scaling factors.
    counts_matrix: 2D numpy array (genes x samples) of fitted fragment counts.
    Returns: 1D numpy array of scaling factors (length = num_samples).
    """
    n_genes, n_samples = counts_matrix.shape
    if n_samples <= 1:
        return np.ones(n_samples, dtype=float)

    # Library sizes
    lib_sizes = np.sum(counts_matrix, axis=0)
    # Avoid zero division
    lib_sizes[lib_sizes == 0] = 1.0

    # Pick reference sample (sample whose 75th percentile is closest to the mean 75th percentile)
    q75 = np.percentile(counts_matrix, 75, axis=0) / lib_sizes
    ref_idx = int(np.argmin(np.abs(q75 - np.mean(q75))))
    ref_counts = counts_matrix[:, ref_idx]
    ref_lib = lib_sizes[ref_idx]

    tmm_factors = np.ones(n_samples, dtype=float)

    for j in range(n_samples):
        if j == ref_idx:
            continue

        obs_counts = counts_matrix[:, j]
        obs_lib = lib_sizes[j]

        # Valid genes with non-zero counts in both samples
        valid = (obs_counts > 0) & (ref_counts > 0)
        if np.sum(valid) < 50:
            continue

        # Log fold changes (M) and average absolute levels (A)
        m_vals = np.log2((obs_counts[valid] / obs_lib) / (ref_counts[valid] / ref_lib))
        a_vals = 0.5 * (np.log2(obs_counts[valid] / obs_lib) + np.log2(ref_counts[valid] / ref_lib))

        # Trim M (30% on each tail) and A (5% on each tail)
        m_low, m_high = np.percentile(m_vals, [30, 70])
        a_low, a_high = np.percentile(a_vals, [5, 95])

        keep = (m_vals >= m_low) & (m_vals <= m_high) & (a_vals >= a_low) & (a_vals <= a_high)
        if np.sum(keep) < 20:
            continue

        # Approximate inverse variance weights
        w = (obs_lib - obs_counts[valid][keep]) / (obs_lib * obs_counts[valid][keep]) + \
            (ref_lib - ref_counts[valid][keep]) / (ref_lib * ref_counts[valid][keep])
        weights = 1.0 / np.maximum(w, 1e-6)

        tmm_log = np.sum(weights * m_vals[keep]) / np.sum(weights)
        tmm_factors[j] = 2.0 ** tmm_log

    # Normalize factors so their geometric mean is 1
    geo_mean = np.exp(np.mean(np.log(np.maximum(tmm_factors, 1e-6))))
    return tmm_factors / geo_mean


def run_combat_adjustment(data_matrix, batch_ids):
    """
    Parametric Empirical Bayes Batch Effect Adjustment (Johnson, Li, & Rabinovic, 2007).
    data_matrix: 2D numpy array (genes x samples) of log-expression values.
    batch_ids: list or array of batch labels for each sample.
    Returns: batch-adjusted 2D numpy array.
    """
    unique_batches, batch_indices = np.unique(batch_ids, return_inverse=True)
    n_batches = len(unique_batches)

    # Need at least 2 batches with >= 2 samples each to perform ComBat
    batch_counts = np.bincount(batch_indices)
    if n_batches < 2 or np.all(batch_counts < 2):
        return data_matrix

    n_genes, n_samples = data_matrix.shape

    # Overall gene means and variances
    grand_mean = np.mean(data_matrix, axis=1, keepdims=True)
    var_pooled = np.mean((data_matrix - grand_mean) ** 2, axis=1, keepdims=True)
    var_pooled[var_pooled <= 1e-8] = 1.0
    std_pooled = np.sqrt(var_pooled)

    # Standardize data
    standardized = (data_matrix - grand_mean) / std_pooled

    # Batch effects: gamma (additive shift) and delta (multiplicative variance)
    adjusted = np.zeros_like(data_matrix)

    for b in range(n_batches):
        idx = np.where(batch_indices == b)[0]
        n_b = len(idx)
        if n_b < 2:
            adjusted[:, idx] = data_matrix[:, idx]
            continue

        batch_data = standardized[:, idx]
        gamma_hat = np.mean(batch_data, axis=1, keepdims=True)
        delta_hat = np.var(batch_data, axis=1, ddof=1, keepdims=True)
        delta_hat[delta_hat <= 1e-8] = 1.0

        # Empirical Bayes shrinkage towards prior means
        gamma_bar = np.mean(gamma_hat)
        tau_bar = np.var(gamma_hat, ddof=1) if n_genes > 1 else 1.0
        delta_bar = np.mean(delta_hat)

        # Shrinkage weight
        weight = tau_bar / (tau_bar + delta_hat / n_b)
        gamma_star = weight * gamma_hat + (1.0 - weight) * gamma_bar

        # Adjust standardized data: (standardized - gamma*) / sqrt(delta*)
        batch_adj = (batch_data - gamma_star) / np.sqrt(np.maximum(delta_hat, 1e-4))
        # Project back to original scale
        adjusted[:, idx] = batch_adj * std_pooled + grand_mean

    # Non-negative expression values
    return np.maximum(adjusted, 0.0)


# ---------------------------------------------------------------------------
# Database Storage Layer (DuckDB / SQLite)
# ---------------------------------------------------------------------------

def initialize_database(db_path, use_duckdb=True):
    """
    Create and initialize tables in DuckDB or SQLite with uniform schema.
    """
    if use_duckdb and HAS_DUCKDB:
        conn = duckdb.connect(str(db_path))
    else:
        conn = sqlite3.connect(str(db_path))

    cursor = conn.cursor()

    # 1. Clinical Samples Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS clinical_samples (
        sample_id VARCHAR PRIMARY KEY,
        patient_id VARCHAR,
        cancer_type VARCHAR,
        cancer_type_detailed VARCHAR,
        cuppa_predicted_cancer_type VARCHAR,
        cuppa_likelihood DOUBLE,
        purity DOUBLE,
        ploidy DOUBLE,
        tmb_per_mb DOUBLE,
        msi_status VARCHAR,
        ecdna_status VARCHAR,
        mean_telomere_length_tumour_bp DOUBLE,
        batch_id VARCHAR,
        clinical_json VARCHAR
    )
    """)

    # 2. Gene Expression Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS gene_expression (
        sample_id VARCHAR,
        gene_name VARCHAR,
        ensembl_gene_id VARCHAR,
        chromosome VARCHAR,
        gene_length INTEGER,
        fitted_fragments DOUBLE,
        norm_counts DOUBLE,
        raw_tpm DOUBLE,
        adj_tpm DOUBLE,
        log2_adj_tpm DOUBLE,
        log2_norm_counts DOUBLE,
        batch_corrected_expr DOUBLE,
        percentile_cancer DOUBLE,
        median_tpm_cancer DOUBLE,
        percentile_cohort DOUBLE,
        median_tpm_cohort DOUBLE,
        PRIMARY KEY (sample_id, gene_name)
    )
    """)

    # 3. Transcript Isoforms Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS transcript_expression (
        sample_id VARCHAR,
        gene_name VARCHAR,
        transcript_id VARCHAR,
        transcript_name VARCHAR,
        is_canonical INTEGER,
        exon_count INTEGER,
        transcript_length INTEGER,
        effective_length DOUBLE,
        fitted_fragments DOUBLE,
        tpm DOUBLE,
        log2_tpm DOUBLE,
        isoform_fraction DOUBLE,
        PRIMARY KEY (sample_id, transcript_id)
    )
    """)

    # 4. Alternate Splice Junctions Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS alt_splice_junctions (
        sample_id VARCHAR,
        gene_name VARCHAR,
        chromosome VARCHAR,
        start_pos INTEGER,
        end_pos INTEGER,
        strand VARCHAR,
        junction_type VARCHAR,
        fragment_support INTEGER,
        cohort_frequency INTEGER
    )
    """)

    # 5. RNA Fusions Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS rna_fusions (
        sample_id VARCHAR,
        upstream_gene VARCHAR,
        downstream_gene VARCHAR,
        in_frame INTEGER,
        split_reads INTEGER,
        discordant_reads INTEGER,
        total_support INTEGER,
        fusion_type VARCHAR
    )
    """)

    # 6. Database Metadata Table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS database_metadata (
        key VARCHAR PRIMARY KEY,
        value VARCHAR
    )
    """)

    # 7. Cohort Landscape Coordinates Table (t-SNE, PCA, UMAP)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS cohort_umap_coordinates (
        sample_id VARCHAR PRIMARY KEY,
        tsne_1 DOUBLE,
        tsne_2 DOUBLE,
        pca_1 DOUBLE,
        pca_2 DOUBLE,
        pca_3 DOUBLE,
        umap_1 DOUBLE,
        umap_2 DOUBLE,
        specimen_category VARCHAR,
        expr_erbb2 DOUBLE,
        expr_tacstd2 DOUBLE,
        expr_cldn18 DOUBLE,
        expr_folr1 DOUBLE,
        expr_vtcn1 DOUBLE
    )
    """)

    conn.commit()
    return conn



def _fast_tsne_embedding(X_pca, perplexity=30.0, n_iter=450, early_exag=12.0, lr=200.0, seed=42):
    """NumPy implementation of t-SNE with Student-t heavy tails & early exaggeration on PCA space."""
    np.random.seed(seed)
    N, D = X_pca.shape
    if N <= 3:
        return np.zeros((N, 2))

    sq_norms = np.sum(X_pca**2, axis=1)
    dists_sq = np.maximum(0, sq_norms[:, None] + sq_norms[None, :] - 2.0 * np.dot(X_pca, X_pca.T))
    np.fill_diagonal(dists_sq, 0)

    target_entropy = np.log(perplexity)
    P = np.zeros((N, N), dtype=np.float32)
    for i in range(N):
        d_i = dists_sq[i]
        mask = np.ones(N, dtype=bool)
        mask[i] = False
        d_i_nz = d_i[mask]

        beta_min, beta_max = -np.inf, np.inf
        beta = 1.0
        for _ in range(40):
            p_i = np.exp(-d_i_nz * beta)
            sum_p = np.sum(p_i)
            if sum_p == 0:
                h = 0
            else:
                h = np.log(sum_p) + beta * np.sum(d_i_nz * p_i) / sum_p
            h_diff = h - target_entropy
            if abs(h_diff) < 1e-4:
                break
            if h_diff > 0:
                beta_min = beta
                beta = beta * 2.0 if beta_max == np.inf else (beta + beta_max) / 2.0
            else:
                beta_max = beta
                beta = beta / 2.0 if beta_min == -np.inf else (beta + beta_min) / 2.0

        P[i, mask] = p_i / max(1e-12, sum_p)

    P = (P + P.T) / (2.0 * N)
    P = np.maximum(P, 1e-12)

    Y = X_pca[:, :2] / (np.std(X_pca[:, 0]) * 1e4)
    Y_prev = Y.copy()
    gains = np.ones_like(Y)

    for it in range(n_iter):
        P_eff = P * early_exag if it < 150 else P

        diff = Y[:, None, :] - Y[None, :, :]
        sq_dist = np.sum(diff**2, axis=-1)
        inv_dist = 1.0 / (1.0 + sq_dist)
        np.fill_diagonal(inv_dist, 0)

        sum_q = np.sum(inv_dist)
        Q = inv_dist / max(1e-12, sum_q)

        grad = 4.0 * np.sum(((P_eff - Q) * inv_dist)[:, :, None] * diff, axis=1)
        momentum = 0.5 if it < 150 else 0.8

        gains = np.where(np.sign(grad) != np.sign(Y - Y_prev), gains + 0.2, gains * 0.8)
        gains = np.maximum(gains, 0.01)

        inc = momentum * (Y - Y_prev) - lr * (gains * grad)
        Y_prev = Y.copy()
        Y = Y + inc

    return Y


def compute_and_store_cohort_umap(conn, active_samples, sorted_genes, use_duckdb=True):
    """Compute 2D coordinates (t-SNE, PCA PC1/2, PC1/3) across cohort using top 1,500 variable genes."""
    cursor = conn.cursor()
    sample_to_idx = {s: i for i, s in enumerate(active_samples)}
    N = len(active_samples)
    if N < 5:
        return

    # 1. Top variable genes (mean > 0.5, preserving biological scale)
    if use_duckdb:
        q_var = """
            SELECT gene_name, STDDEV(log2_adj_tpm) as var_score
            FROM gene_expression
            GROUP BY gene_name
            HAVING AVG(log2_adj_tpm) > 0.5 AND STDDEV(log2_adj_tpm) IS NOT NULL
            ORDER BY var_score DESC
            LIMIT 1500
        """
    else:
        q_var = """
            SELECT gene_name, (AVG(log2_adj_tpm * log2_adj_tpm) - AVG(log2_adj_tpm) * AVG(log2_adj_tpm)) as var_score
            FROM gene_expression
            GROUP BY gene_name
            HAVING AVG(log2_adj_tpm) > 0.5
            ORDER BY var_score DESC
            LIMIT 1500
        """
    cursor.execute(q_var)
    top_genes = [r[0] for r in cursor.fetchall()]
    gene_to_idx = {g: i for i, g in enumerate(top_genes)}
    G = len(top_genes)

    # 2. Build expression matrix
    mat = np.zeros((N, G), dtype=np.float32)
    batch_size = 400
    for b_start in range(0, G, batch_size):
        b_genes = top_genes[b_start:b_start + batch_size]
        quoted = ",".join(f"'{g.replace(chr(39), chr(39)+chr(39))}'" for g in b_genes)
        q_mat = f"SELECT sample_id, gene_name, log2_adj_tpm FROM gene_expression WHERE gene_name IN ({quoted})"
        cursor.execute(q_mat)
        for s, g, val in cursor.fetchall():
            if s in sample_to_idx and g in gene_to_idx:
                mat[sample_to_idx[s], gene_to_idx[g]] = val or 0.0

    # 3. Mean-centered SVD/PCA (do NOT divide by std; preserves dynamic lineage range)
    mat_centered = mat - np.mean(mat, axis=0)
    U, S, Vt = np.linalg.svd(mat_centered, full_matrices=False)

    def norm_span(pts):
        p_min, p_max = np.min(pts, axis=0), np.max(pts, axis=0)
        span = np.maximum(1e-5, p_max - p_min)
        return ((pts - p_min) / span * 20.0) - 10.0

    pca_raw_1_2 = U[:, :2] * S[:2]
    pca_norm_1_2 = norm_span(pca_raw_1_2)

    pca_raw_1_3 = np.column_stack([U[:, 0] * S[0], U[:, 2] * S[2]])
    pca_norm_1_3 = norm_span(pca_raw_1_3)

    # 4. Run t-SNE with early exaggeration on top 30 PCs
    pca_30 = U[:, :30] * S[:30]
    tsne_raw = _fast_tsne_embedding(pca_30, perplexity=30.0, n_iter=450)
    tsne_norm = norm_span(tsne_raw)

    # 5. Extract specimen categories
    cursor.execute("SELECT sample_id, clinical_json FROM clinical_samples")
    specimen_cat_map = {}
    for s_id, c_json in cursor.fetchall():
        cat = "Unknown"
        if c_json:
            try:
                d = json.loads(c_json)
                st = d.get("SPECIMEN_TYPE", "")
                if "primary" in st.lower():
                    cat = "Primary"
                elif "metastat" in st.lower():
                    cat = "Metastasis"
                elif "recurrent" in st.lower():
                    cat = "Local Recurrence"
                elif st:
                    cat = "Other"
            except Exception:
                pass
        specimen_cat_map[s_id] = cat

    # 6. Horizon-Beacon target expressions
    horizon_genes = ['ERBB2', 'TACSTD2', 'CLDN18', 'FOLR1', 'VTCN1']
    hb_expr = {g: {} for g in horizon_genes}
    for g in horizon_genes:
        cursor.execute("SELECT sample_id, log2_adj_tpm FROM gene_expression WHERE gene_name = ?", [g])
        for s, val in cursor.fetchall():
            hb_expr[g][s] = round(val or 0.0, 3)

    # 7. Store rows
    insert_rows = []
    for s in active_samples:
        i = sample_to_idx[s]
        t1 = round(float(tsne_norm[i, 0]), 2)
        t2 = round(float(tsne_norm[i, 1]), 2)
        p1 = round(float(pca_norm_1_2[i, 0]), 2)
        p2 = round(float(pca_norm_1_2[i, 1]), 2)
        p3 = round(float(pca_norm_1_3[i, 1]), 2)
        u1 = round(float(pca_norm_1_2[i, 0]), 2)
        u2 = round(float(pca_norm_1_2[i, 1]), 2)
        cat = specimen_cat_map.get(s, "Unknown")
        e_erbb2 = hb_expr['ERBB2'].get(s, 0.0)
        e_tacstd2 = hb_expr['TACSTD2'].get(s, 0.0)
        e_cldn18 = hb_expr['CLDN18'].get(s, 0.0)
        e_folr1 = hb_expr['FOLR1'].get(s, 0.0)
        e_vtcn1 = hb_expr['VTCN1'].get(s, 0.0)
        insert_rows.append((s, t1, t2, p1, p2, p3, u1, u2, cat, e_erbb2, e_tacstd2, e_cldn18, e_folr1, e_vtcn1))

    cursor.executemany("INSERT INTO cohort_umap_coordinates VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", insert_rows)
    conn.commit()
    print(f"  [Landscape] Successfully generated and stored t-SNE & PCA coordinates for {len(active_samples)} samples.")


def create_database_indexes(conn, use_duckdb=True):
    """Create indexes for high-speed exploratory querying."""
    cursor = conn.cursor()
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene_name ON gene_expression (gene_name)")
    if not use_duckdb:
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene_name_nocase ON gene_expression (gene_name COLLATE NOCASE)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_sample_id ON gene_expression (sample_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_trans_gene ON transcript_expression (gene_name)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_trans_sample ON transcript_expression (sample_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_clin_cancer ON clinical_samples (cancer_type)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_alt_gene ON alt_splice_junctions (gene_name)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_fus_up ON rna_fusions (upstream_gene)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_fus_down ON rna_fusions (downstream_gene)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_umap_sample ON cohort_umap_coordinates (sample_id)")
    conn.commit()


# ---------------------------------------------------------------------------
# Main ETL Workflow
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build analytical database from HMF Isofox RNA-Seq and cBioPortal clinical outputs."
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        nargs="+",
        help="One or more directories containing sample Oncoanalyser / Isofox output folders."
    )
    parser.add_argument(
        "--cbio-dir",
        default=None,
        nargs="+",
        help="One or more directories containing cBioPortal data_clinical_sample.txt and data_clinical_patient.txt."
    )
    parser.add_argument(
        "--metadata-tsv",
        default=None,
        help="Optional clinical metadata TSV/CSV."
    )
    parser.add_argument(
        "--gene-distribution",
        default=None,
        help="Path to HMF isofox.hmf_38_151_2600.gene_distribution.csv.gz (auto-discovered if omitted)."
    )
    parser.add_argument(
        "--output-db",
        default="hmf_rnaseq_cohort.duckdb",
        help="Output database file path (e.g. hmf_rnaseq_cohort.duckdb or hmf_rnaseq_cohort.sqlite3)."
    )
    parser.add_argument(
        "--batch-column",
        default="BATCH_ID",
        help="Clinical attribute column to use for ComBat batch adjustment (default: BATCH_ID)."
    )
    parser.add_argument(
        "--skip-combat",
        action="store_true",
        help="Skip ComBat batch adjustment and only compute TMM library size normalization."
    )
    parser.add_argument(
        "--threads", "--cores", "-t",
        dest="threads",
        type=int,
        default=min(8, os.cpu_count() or 4),
        help="Number of parallel worker processes to use (default: min(8, cpu_count))."
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "duckdb", "sqlite"],
        default="auto",
        help="Database backend to use (default: auto, chooses duckdb if available else sqlite)."
    )

    args = parser.parse_args()

    # Determine backend
    use_duckdb = False
    if args.backend == "duckdb":
        if not HAS_DUCKDB:
            print("Error: DuckDB requested but 'duckdb' python package is not installed. Run: pip install duckdb", file=sys.stderr)
            sys.exit(1)
        use_duckdb = True
    elif args.backend == "sqlite":
        use_duckdb = False
    else:
        # Auto
        use_duckdb = HAS_DUCKDB and (not args.output_db.endswith(".sqlite") and not args.output_db.endswith(".sqlite3"))

    backend_label = "DuckDB" if use_duckdb else "SQLite3"
    print(f"=== Starting Pan-Cancer RNA-Seq Database Builder ({backend_label}) ===")

    # 1. Discover Sample Directories
    sample_dirs = find_sample_dirs(args.results_dir)
    print(f"Discovered {len(sample_dirs)} valid sample directories across results paths.")
    if not sample_dirs:
        print("Error: No sample directories with Isofox outputs found.", file=sys.stderr)
        sys.exit(1)

    # 2. Parse Clinical Metadata
    clinical_map = parse_clinical_metadata(cbio_dirs=args.cbio_dir, metadata_tsv=args.metadata_tsv)
    print(f"Loaded clinical attributes for {len(clinical_map)} samples.")

    # 3. Resolve HMF Gene Distribution Reference
    gene_dist_file = resolve_gene_dist_file(args.gene_distribution)
    if gene_dist_file:
        print(f"Using HMF gene distribution reference: {gene_dist_file}")
    else:
        print("Note: HMF gene distribution reference file not provided/found. Will use pre-annotated percentiles or cohort percentiles.")

    # 4. Pass 1: Parallel extraction of gene counts for TMM / ComBat
    threads = max(1, args.threads)
    print(f"Pass 1: Extracting fitted fragment counts across {len(sample_dirs)} samples using {threads} parallel worker(s)...")

    counts_by_sample = {}
    unique_genes = set()

    if threads > 1 and len(sample_dirs) > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            dir_iter = iter(sample_dirs)
            futures = set()
            for _ in range(min(threads * 3, len(sample_dirs))):
                d = next(dir_iter, None)
                if d:
                    futures.add(executor.submit(extract_sample_counts_worker, str(d)))

            done_count = 0
            while futures:
                done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    done_count += 1
                    s_id, c_dict = fut.result()
                    if c_dict:
                        counts_by_sample[s_id] = c_dict
                        unique_genes.update(c_dict.keys())
                    if done_count % 25 == 0 or done_count == len(sample_dirs):
                        print(f"  [Pass 1] Processed {done_count}/{len(sample_dirs)} samples ({len(unique_genes)} unique genes)...")

                    d = next(dir_iter, None)
                    if d:
                        futures.add(executor.submit(extract_sample_counts_worker, str(d)))
    else:
        for idx, s_dir in enumerate(sample_dirs, 1):
            s_id, c_dict = extract_sample_counts_worker(str(s_dir))
            if c_dict:
                counts_by_sample[s_id] = c_dict
                unique_genes.update(c_dict.keys())
            if idx % 25 == 0 or idx == len(sample_dirs):
                print(f"  [Pass 1] Processed {idx}/{len(sample_dirs)} samples ({len(unique_genes)} unique genes)...")

    active_samples = sorted(counts_by_sample.keys())
    sorted_genes = sorted(unique_genes)
    if not active_samples:
        print("Error: No gene expression data could be parsed from any of the discovered sample directories.", file=sys.stderr)
        sys.exit(1)
    print(f"Successfully discovered {len(active_samples)} valid RNA-Seq samples across {len(sorted_genes)} unique genes.")

    # 5. Build Lightweight Count Matrix & Run Global TMM Normalization
    print("Computing global cross-cohort TMM normalization...")
    gene_idx_map = {g: i for i, g in enumerate(sorted_genes)}
    sample_idx_map = {s: j for j, s in enumerate(active_samples)}

    # Matrix: genes x samples (fitted fragments) in float32 to conserve RAM
    counts_mat = np.zeros((len(sorted_genes), len(active_samples)), dtype=np.float32)
    for s_id in active_samples:
        j = sample_idx_map[s_id]
        c_dict = counts_by_sample[s_id]
        for g_name, cnt in c_dict.items():
            i = gene_idx_map.get(g_name)
            if i is not None:
                counts_mat[i, j] = cnt

    del counts_by_sample
    gc.collect()

    tmm_factors = compute_tmm_factors(counts_mat)
    lib_sizes = np.sum(counts_mat, axis=0)
    effective_lib_sizes = lib_sizes * tmm_factors
    median_eff_lib = float(np.median(effective_lib_sizes)) if len(effective_lib_sizes) > 0 else 1.0

    # Normalization scaling factor for each sample: median_eff_lib / effective_lib_sizes[j]
    norm_factors = {}
    for j, s_id in enumerate(active_samples):
        eff_lib = float(effective_lib_sizes[j])
        norm_factors[s_id] = (median_eff_lib / eff_lib) if eff_lib > 0 else 1.0

    print(f"TMM scaling computed across {len(active_samples)} samples (Median effective library: {median_eff_lib:,.0f} fragments).")

    # 6. ComBat Batch Effect Adjustment (Optional)
    batch_adj_lookup = {}
    combat_applied = False
    if not args.skip_combat:
        batch_ids = []
        for s_id in active_samples:
            b_val = clinical_map.get(s_id, {}).get(args.batch_column) or clinical_map.get(s_id, {}).get("BATCH") or clinical_map.get(s_id, {}).get("CENTER") or "BATCH_1"
            batch_ids.append(str(b_val))

        unique_b = np.unique(batch_ids)
        if len(unique_b) > 1:
            print(f"Running ComBat batch effect adjustment across {len(unique_b)} batches ({args.batch_column})...")
            norm_counts_mat = np.zeros_like(counts_mat)
            for j in range(len(active_samples)):
                norm_counts_mat[:, j] = counts_mat[:, j] * norm_factors[active_samples[j]]

            log_norm_mat = np.log2(norm_counts_mat + 1.0)
            log_adj_mat = run_combat_adjustment(log_norm_mat, batch_ids)
            combat_mat = np.maximum((2.0 ** log_adj_mat) - 1.0, 0.0)

            for j, s_id in enumerate(active_samples):
                sample_dict = {}
                for i, g_name in enumerate(sorted_genes):
                    val = float(combat_mat[i, j])
                    if val > 0.01:
                        sample_dict[g_name] = round(val, 2)
                batch_adj_lookup[s_id] = sample_dict

            del norm_counts_mat, log_norm_mat, log_adj_mat, combat_mat
            gc.collect()
            combat_applied = True
            print("ComBat batch effect adjustment complete.")
        else:
            print(f"Only 1 unique batch found for column '{args.batch_column}'. Skipping ComBat adjustment.")

    del counts_mat
    gc.collect()

    # 7. Initialize Database & Insert Records
    out_path = Path(args.output_db).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()

    conn = initialize_database(out_path, use_duckdb=use_duckdb)
    cursor = conn.cursor()

    if not use_duckdb:
        cursor.execute("PRAGMA journal_mode = WAL")
        cursor.execute("PRAGMA synchronous = OFF")
        cursor.execute("PRAGMA cache_size = -64000")
        cursor.execute("PRAGMA temp_store = MEMORY")

    print(f"Writing records to {backend_label} database at {out_path}...")

    # A. Insert Clinical Samples
    clin_rows = []
    for s_id in active_samples:
        c_dict = clinical_map.get(s_id, {})
        c_type = c_dict.get("CANCER_TYPE") or c_dict.get("Cancer_Type") or c_dict.get("CANCER_TYPE_DETAILED") or "UNKNOWN"
        c_det = c_dict.get("CANCER_TYPE_DETAILED") or c_dict.get("DIAGNOSIS") or c_type
        cuppa_type = c_dict.get("CUPPA_PREDICTED_CANCER_TYPE") or "NA"
        cuppa_lh = _safe_float(c_dict.get("CUPPA_LIKELIHOOD"), -1.0)
        purity = _safe_float(c_dict.get("PURITY"), -1.0)
        ploidy = _safe_float(c_dict.get("PLOIDY"), -1.0)
        tmb = _safe_float(c_dict.get("TMBPERMB") or c_dict.get("TMB_PER_MB"), -1.0)
        msi = c_dict.get("MSISTATUS") or c_dict.get("MSI_STATUS") or "NA"
        ecdna = c_dict.get("ECDNA_STATUS") or "NA"
        tel_len = _safe_float(c_dict.get("MEAN_TELOMERE_LENGTH_TUMOUR_BP") or c_dict.get("MEAN_TELOMERE_LENGTH_TUMOR_BP"), -1.0)
        b_id = c_dict.get(args.batch_column) or c_dict.get("BATCH_ID") or c_dict.get("BATCH") or "BATCH_1"

        clin_rows.append((
            s_id,
            c_dict.get("PATIENT_ID") or s_id,
            c_type,
            c_det,
            cuppa_type,
            cuppa_lh,
            purity,
            ploidy,
            tmb,
            msi,
            ecdna,
            tel_len,
            b_id,
            json.dumps(c_dict)
        ))

    cursor.executemany("""
    INSERT INTO clinical_samples VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, clin_rows)
    conn.commit()

    # B. Pass 2: Streaming ingestion of gene, transcript, splice, and fusion records
    print(f"Pass 2: Parsing & streaming records into {backend_label} database using {threads} parallel worker(s)...")

    sample_dir_map = {s.name: s for s in sample_dirs if s.name in active_samples}
    tasks = []
    for s_id in active_samples:
        s_dir = sample_dir_map[s_id]
        n_fac = norm_factors.get(s_id, 1.0)
        b_adj = batch_adj_lookup.get(s_id)
        tasks.append((str(s_dir), n_fac, b_adj))

    total_samples = len(tasks)
    inserted_samples = 0
    total_genes_inserted = 0
    total_trans_inserted = 0

    def insert_sample_records(s_id, gene_rows, trans_rows, splice_rows, fusion_rows):
        nonlocal total_genes_inserted, total_trans_inserted
        if gene_rows:
            cursor.executemany("INSERT INTO gene_expression VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", gene_rows)
            total_genes_inserted += len(gene_rows)
        if trans_rows:
            cursor.executemany("INSERT INTO transcript_expression VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", trans_rows)
            total_trans_inserted += len(trans_rows)
        if splice_rows:
            cursor.executemany("INSERT INTO alt_splice_junctions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", splice_rows)
        if fusion_rows:
            cursor.executemany("INSERT INTO rna_fusions VALUES (?, ?, ?, ?, ?, ?, ?, ?)", fusion_rows)

    if threads > 1 and total_samples > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            task_iter = iter(tasks)
            futures = set()
            # Bounded in-flight queue: at most threads * 2 samples simultaneously in memory
            for _ in range(min(threads * 2, total_samples)):
                t = next(task_iter, None)
                if t:
                    futures.add(executor.submit(parse_and_format_sample_worker, t))

            while futures:
                done, futures = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                for fut in done:
                    s_id, g_rows, t_rows, sp_rows, f_rows = fut.result()
                    insert_sample_records(s_id, g_rows, t_rows, sp_rows, f_rows)
                    del g_rows, t_rows, sp_rows, f_rows
                    inserted_samples += 1

                    if inserted_samples % 10 == 0:
                        conn.commit()

                    if inserted_samples % 10 == 0 or inserted_samples == total_samples:
                        pct = (inserted_samples / total_samples) * 100
                        print(f"  [Pass 2] Inserted {inserted_samples}/{total_samples} samples ({pct:.1f}%) | {total_genes_inserted:,} gene rows...")

                    t = next(task_iter, None)
                    if t:
                        futures.add(executor.submit(parse_and_format_sample_worker, t))
    else:
        for t in tasks:
            s_id, g_rows, t_rows, sp_rows, f_rows = parse_and_format_sample_worker(t)
            insert_sample_records(s_id, g_rows, t_rows, sp_rows, f_rows)
            del g_rows, t_rows, sp_rows, f_rows
            inserted_samples += 1

            if inserted_samples % 10 == 0:
                conn.commit()

            if inserted_samples % 10 == 0 or inserted_samples == total_samples:
                pct = (inserted_samples / total_samples) * 100
                print(f"  [Pass 2] Inserted {inserted_samples}/{total_samples} samples ({pct:.1f}%) | {total_genes_inserted:,} gene rows...")

    conn.commit()

    # C. Insert Metadata
    meta_rows = [
        ("created_at", datetime.datetime.now().isoformat()),
        ("backend", backend_label),
        ("total_samples", str(len(active_samples))),
        ("total_genes", str(len(sorted_genes))),
        ("batch_column", args.batch_column),
        ("combat_applied", "True" if combat_applied else "False"),
    ]
    cursor.executemany("INSERT INTO database_metadata VALUES (?, ?)", meta_rows)
    conn.commit()

    # 7.5 Compute and store UMAP landscape
    print("Computing whole-cohort 2D UMAP embedding & target expressions...")
    compute_and_store_cohort_umap(conn, active_samples, sorted_genes, use_duckdb=use_duckdb)

    # 8. Create Indexes
    print("Building high-performance database indexes...")
    create_database_indexes(conn, use_duckdb=use_duckdb)

    if not use_duckdb:
        cursor.execute("PRAGMA optimize")
        cursor.execute("PRAGMA synchronous = NORMAL")

    conn.close()

    db_size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"\nSUCCESS: Database generated at: {out_path} ({db_size_mb:.2f} MB)")
    print(f"Summary: {len(active_samples)} samples, {len(sorted_genes)} genes, global TMM normalization applied.")


def _safe_float(val, default):
    try:
        if val is None or str(val).upper() in ["NA", "NAN", "NULL", ""]:
            return default
        return float(val)
    except (ValueError, TypeError):
        return default


def _safe_int(val, default):
    try:
        if val is None or str(val).upper() in ["NA", "NAN", "NULL", ""]:
            return default
        return int(float(val))
    except (ValueError, TypeError):
        return default


if __name__ == "__main__":
    main()

