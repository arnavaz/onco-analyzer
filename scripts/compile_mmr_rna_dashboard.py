#!/usr/bin/env python3
"""
MMR Cohort RNA-Seq & TCR Interactive Clinical Dashboard Compiler
================================================================
Aggregates Oncoanalyser RNA-Seq outputs (STAR alignment metrics, Isofox
transcript & gene-level expression, Isofox structural RNA fusions, CIDER
TCR/BCR immune repertoire reconstruction, and CUPPA tissue of origin)
alongside multi-surgery clinical metadata into an interactive, standalone
analytical HTML dashboard.

Features:
- Integrated Clinical Metadata (TissueOrigin, PathDx, Timepoints, HLA Class I alleles)
- Cohort Landscape & OncoPrint Biomarker Matrix
- Interactive Differential Expression Volcano Plot (Welch's t-test between any 2 categories)
- Dynamic Gene Search & Categorized Expression Trajectories (CTAs, GBM, GI, Immune, HLA)
- Longitudinal Multi-Surgery Expression Tracking (Primary vs Sequential Recurrences)
- Reconstructed Immune Repertoire Explorer (CDR3 sequences, V/J usage, clonality)
- In-Frame Structural RNA Fusions & Chimeric Transcript Catalog
- CUPPA Transcriptomic Tissue Lineage Predictions & Concordance
- Sample Master Catalog with Slide-Over Clinical Inspection Drawer

Theme & Styling:
Adheres strictly to the design-taste-frontend anti-slop guidelines for app screens
(clean clinical light palette, flat 1px borders, typography-driven KPIs, right-aligned
numeric tables, and restrained data-focused color accents).
"""

import os
import sys
import re
import json
import math
import argparse
import datetime
import csv
import gzip
import bisect
from array import array
from pathlib import Path

# Add scripts directory to path to import hla_binding and run_rna_postprocessing
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

try:
    from hla_binding import BindingModel, allele_key, allele_label, junction_peptides
    from run_rna_postprocessing import score_tcr_neoepitope_interaction, resolve_binding_dir
except ImportError:
    BindingModel = None
    allele_key = lambda x: x
    allele_label = lambda x: x
    junction_peptides = lambda s, idx, **kw: []
    score_tcr_neoepitope_interaction = lambda p, t: 0.5
    resolve_binding_dir = lambda p=None: None


# ---------------------------------------------------------------------------
# CLI Argument Parser
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile Oncoanalyser RNA-Seq & TCR results into an interactive clinical dashboard."
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Root directory containing Oncoanalyser sample outputs."
    )
    parser.add_argument(
        "--metadata-tsv",
        required=True,
        help="Path to MMR_RNAseq.tsv clinical metadata file."
    )
    parser.add_argument(
        "--output-html",
        default=None,
        help="Path for destination standalone HTML file (defaults to <results_dir>/cohort_mmr_rna_report.html)."
    )
    parser.add_argument(
        "--max-cdr3-per-sample",
        type=int,
        default=60,
        help="Maximum PASS CDR3 clonotypes to embed per sample and receptor type (TCR / BCR; default: 60)."
    )
    parser.add_argument(
        "--antigen-tsv",
        default=str(DEFAULT_ANTIGEN_TSV),
        help="Tumour antigen annotation table (default: resources/antigens/tumor_antigen_annotation.tsv)."
    )
    parser.add_argument(
        "--max-outliers-per-sample",
        type=int,
        default=250,
        help="Maximum outlier genes to embed per sample (known antigens are always kept; default: 250)."
    )
    parser.add_argument(
        "--gene-distribution",
        default=DEFAULT_GENE_DISTRIBUTION,
        help="Hartwig Isofox gene_distribution.csv.gz (per-cancer-type TPM percentiles) used for cancer-type-matched "
             "percentiles; pass '' to skip (default: H4H oncoanalyser 3.1.0 resources)."
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Helper: Metadata Parsing
# ---------------------------------------------------------------------------

def parse_metadata(metadata_tsv):
    """Parse MMR_RNAseq.tsv clinical metadata into structured records."""
    by_sample = {}
    records = []
    p = Path(metadata_tsv)
    if not p.is_file():
        print(f"Warning: Metadata TSV not found at {metadata_tsv}", file=sys.stderr)
        return by_sample, records

    with open(p, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            patient_id = str(row.get("Patient_ID", "")).strip()
            tissue_id = str(row.get("Tissue_ID", "")).strip()
            rna_id = str(row.get("RNAseq_Sample_ID", "")).strip()
            tcr_id = str(row.get("TCRseq_Sample_ID", "")).strip()
            path_dx = str(row.get("PathDx", "")).strip()
            tissue_origin = str(row.get("TissueOrigin", "")).strip()
            tissue_format = str(row.get("TissueFormat", "")).strip()
            path_type = str(row.get("PathType", "")).strip()
            date_collected = str(row.get("DateCollected", "")).strip()
            date_iso = ""
            for fmt in ("%d-%b-%y", "%d-%b-%Y", "%Y-%m-%d", "%m/%d/%Y"):
                try:
                    date_iso = datetime.datetime.strptime(date_collected, fmt).date().isoformat()
                    break
                except ValueError:
                    continue
            rna_assay = str(row.get("RNA-seq available?", "")).strip()
            library_type = "Target Capture" if "capture" in rna_assay.lower() else ("Bulk" if rna_assay and rna_assay != "No" else "")

            # HLA alleles
            hla_a = [x for x in [str(row.get("A1", "")).strip(), str(row.get("A2", "")).strip()] if x and x != "nan"]
            hla_b = [x for x in [str(row.get("B1", "")).strip(), str(row.get("B2", "")).strip()] if x and x != "nan"]
            hla_c = [x for x in [str(row.get("C1", "")).strip(), str(row.get("C2", "")).strip()] if x and x != "nan"]

            # Parse surgery / longitudinal timepoint
            timepoint = "Primary / Initial"
            dx_lower = path_dx.lower()
            if "surgery 1" in dx_lower:
                timepoint = "Surgery 1"
            elif "surgery 2" in dx_lower:
                timepoint = "Surgery 2"
            elif "surgery 3" in dx_lower:
                timepoint = "Surgery 3"
            elif "surgery 4" in dx_lower:
                timepoint = "Surgery 4"
            elif "recurrence" in dx_lower:
                timepoint = "Recurrence"
            elif "postmortem" in dx_lower or "autopsy" in dx_lower:
                timepoint = "Postmortem"

            rec = {
                "patient_id": patient_id,
                "tissue_id": tissue_id,
                "rna_id": rna_id if rna_id and rna_id != "nan" else tissue_id,
                "tcr_id": tcr_id,
                "path_dx": path_dx,
                "tissue_origin": tissue_origin,
                "tissue_format": tissue_format,
                "path_type": path_type,
                "date_collected": date_collected,
                "date_iso": date_iso,
                "library_type": library_type,
                "timepoint": timepoint,
                "hla_a": hla_a,
                "hla_b": hla_b,
                "hla_c": hla_c,
            }
            records.append(rec)

            # Map multiple keys for robust lookup
            for k in [tissue_id, rna_id, tissue_id.replace("_", "-"), rna_id.replace("_", "-")]:
                if k and k != "nan":
                    by_sample[k] = rec
            # Sub-sampled tissue rows (e.g. MMR1396T03A/B) share one RNA library run as MMR1396T03
            m = re.match(r"^(.*T\d+)[A-Z]$", tissue_id)
            if m and m.group(1) not in by_sample:
                by_sample[m.group(1)] = rec

    return by_sample, records


# ---------------------------------------------------------------------------
# Helper: Oncoanalyser Directory Discovery
# ---------------------------------------------------------------------------

def locate_sample_dirs(results_dir):
    """
    Find all sample output folders under results_dir.
    Resolves both top-level and nested Oncoanalyser directory layouts.
    """
    sample_dirs = {}
    p = Path(results_dir)
    if not p.is_dir():
        return sample_dirs

    tools = ["isofox", "cider", "cuppa", "alignments", "star"]

    for child in p.iterdir():
        if child.is_dir() and not child.name.startswith((".", "work", "pipeline_info", "sample_sheets", "tmp")):
            # Check nested child/child.name directory (standard Oncoanalyser cluster layout)
            nested = child / child.name
            if nested.is_dir() and any((nested / sub).is_dir() for sub in tools):
                sample_dirs[child.name] = nested
            elif any((child / sub).is_dir() for sub in tools):
                sample_dirs[child.name] = child
            else:
                # Check child's direct subdirectories
                for sub in child.iterdir():
                    if sub.is_dir() and not sub.name.startswith("."):
                        if any((sub / s).is_dir() for s in tools):
                            sample_dirs[sub.name] = sub

    return sample_dirs


def resolve_tool_dir(sample_dir, tool_name):
    """Resolve tool directory whether direct or nested."""
    if (sample_dir / tool_name).is_dir():
        return sample_dir / tool_name
    nested = sample_dir / sample_dir.name / tool_name
    if nested.is_dir():
        return nested
    return sample_dir / tool_name


# ---------------------------------------------------------------------------
# Helper: Oncoanalyser Tool Parsers
# ---------------------------------------------------------------------------

def parse_isofox_summary(sample_dir):
    """Parse Isofox summary metrics (*.isf.summary.tsv)."""
    metrics = {
        "qc_status": "PASS",
        "total_fragments": 0,
        "duplicate_fragments": 0,
        "spliced_fragment_pct": 0.0,
        "spliced_gene_count": 0,
        "median_frag_length": 0,
        "median_gc_ratio": 0.0,
        "unspliced_fragment_pct": 0.0,
        "chimeric_fragment_pct": 0.0,
        "forward_strand_pct": None,
        "read_length": 0,
    }
    isofox_dir = resolve_tool_dir(sample_dir, "isofox")
    if not isofox_dir.is_dir():
        return metrics

    files = list(isofox_dir.glob("*.isf.summary.tsv")) + list(isofox_dir.glob("*summary.tsv"))
    if not files:
        return metrics

    try:
        with open(files[0], "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                metrics["qc_status"] = row.get("QcStatus", "PASS")
                try:
                    metrics["total_fragments"] = int(row.get("TotalFragments", 0) or 0)
                    metrics["duplicate_fragments"] = int(row.get("DuplicateFragments", 0) or 0)
                    metrics["spliced_fragment_pct"] = float(row.get("SplicedFragmentPerc", 0) or 0)
                    metrics["spliced_gene_count"] = int(row.get("SplicedGeneCount", 0) or 0)
                    metrics["median_frag_length"] = int(row.get("FragLength50th", 0) or 0)
                    metrics["median_gc_ratio"] = float(row.get("MedianGCRatio", 0) or 0)
                    metrics["unspliced_fragment_pct"] = float(row.get("UnsplicedFragmentPerc", 0) or 0)
                    metrics["chimeric_fragment_pct"] = float(row.get("ChimericFragmentPerc", 0) or 0)
                    if row.get("ForwardStrandPercent") not in (None, ""):
                        metrics["forward_strand_pct"] = float(row["ForwardStrandPercent"])
                    metrics["read_length"] = int(row.get("ReadLength", 0) or 0)
                except (ValueError, TypeError):
                    pass
                break
    except Exception as e:
        print(f"Warning: Failed parsing Isofox summary in {sample_dir}: {e}", file=sys.stderr)

    return metrics


def parse_star_log(sample_dir):
    """Parse STAR Log.final.out for alignment metrics."""
    metrics = {
        "total_reads": None,
        "uniquely_mapped_pct": None,
        "multi_mapped_pct": None,
        "unmapped_pct": None,
        "avg_input_len": None,
    }
    star_dir = resolve_tool_dir(sample_dir, "star")
    logs = list(star_dir.glob("*.Log.final.out")) if star_dir.is_dir() else []
    if not logs:
        logs = list(sample_dir.glob("*Log.final.out"))
    if not logs:
        # STAR logs are not published by oncoanalyser; fall back to the run's Nextflow work dir (newest wins)
        logs = sorted(sample_dir.parent.glob("work/*/*/*Log.final.out"), key=lambda x: x.stat().st_mtime, reverse=True)

    if not logs:
        return metrics

    try:
        with open(logs[0], "r", encoding="utf-8") as f:
            for line in f:
                if "Number of input reads" in line:
                    metrics["total_reads"] = int(line.split("|")[-1].strip())
                elif "Uniquely mapped reads %" in line:
                    metrics["uniquely_mapped_pct"] = float(line.split("|")[-1].strip().replace("%", ""))
                elif "% of reads mapped to multiple loci" in line:
                    metrics["multi_mapped_pct"] = float(line.split("|")[-1].strip().replace("%", ""))
                elif "% of reads unmapped: too short" in line:
                    metrics["unmapped_pct"] = float(line.split("|")[-1].strip().replace("%", ""))
                elif "Average input read length" in line:
                    metrics["avg_input_len"] = float(line.split("|")[-1].strip())
    except Exception as e:
        print(f"Warning: Failed parsing STAR log in {sample_dir}: {e}", file=sys.stderr)

    return metrics


def parse_isofox_expression(sample_dir):
    """Parse Isofox gene expression (*.isf.gene_data.tsv).

    Returns (gene_tpm, gene_ref) where gene_ref holds the Hartwig reference columns per gene:
    (PercentileCohort, MedianTpmCohort, PercentileCancer, MedianTpmCancer); -1 marks unavailable.
    """
    gene_tpm = {}
    gene_ref = {}
    isofox_dir = resolve_tool_dir(sample_dir, "isofox")
    if not isofox_dir.is_dir():
        return gene_tpm, gene_ref

    files = list(isofox_dir.glob("*.isf.gene_data.*")) + list(isofox_dir.glob("*.gene_data.*"))
    if not files:
        return gene_tpm, gene_ref

    target = files[0]
    sep = "\t" if target.suffix in [".tsv", ".txt"] else ","

    try:
        opener = gzip.open(target, "rt", encoding="utf-8") if target.suffix == ".gz" else open(target, "r", encoding="utf-8")
        with opener as f:
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                gene = row.get("GeneName") or row.get("gene_name") or row.get("GeneId") or row.get("gene_id")
                if not gene:
                    continue
                gene = gene.strip()

                # AdjTPM is primary; RawTPM fallback
                tpm_str = row.get("AdjTPM") or row.get("RawTPM") or row.get("AdjustedTpm") or row.get("Tpm") or row.get("tpm") or "0"
                try:
                    val = float(tpm_str)
                    if gene not in gene_tpm or val > gene_tpm[gene]:
                        gene_tpm[gene] = round(val, 2)
                        gene_ref[gene] = tuple(
                            _safe_float(row.get(col), -1.0)
                            for col in ("PercentileCohort", "MedianTpmCohort", "PercentileCancer", "MedianTpmCancer")
                        )
                except ValueError:
                    continue
    except Exception as e:
        print(f"Warning: Failed parsing Isofox gene data in {sample_dir}: {e}", file=sys.stderr)

    return gene_tpm, gene_ref


def _safe_float(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Hartwig cancer-type expression reference (Isofox gene_distribution)
# ---------------------------------------------------------------------------

DEFAULT_GENE_DISTRIBUTION = (
    "/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0/"
    "hmf_pipeline_resources.38_v3.1.0--1/rna/isofox.hmf_38_151_2600.gene_distribution.csv.gz"
)

# (PathDx pattern, Hartwig CancerType); first match wins, then TissueOrigin fallback, then the pan-cancer ALL rows
HARTWIG_DX_MAP = [
    (re.compile(r"gbm|glio|astrocytom|oligodendro|ependym|medullo|neuroblastom", re.I), "Nervous system"),
    (re.compile(r"colorectal|colon|rectal|polyp", re.I), "Colorectum"),
    (re.compile(r"stomach|gastric", re.I), "Stomach"),
    (re.compile(r"small intestin|duoden", re.I), "Small intestine"),
    (re.compile(r"urothel|transitional cell|bladder", re.I), "Urothelial tract"),
    (re.compile(r"endometri|uter", re.I), "Uterus"),
]
HARTWIG_ORIGIN_MAP = {"CNS": "Nervous system", "GI": "Colorectum", "UT": "Urothelial tract"}


def hartwig_cancer_type(path_dx, tissue_origin):
    """Map clinical diagnosis to a Hartwig gene_distribution CancerType; returns (cancer_type, matched_by)."""
    for pattern, cancer_type in HARTWIG_DX_MAP:
        if pattern.search(path_dx or ""):
            return cancer_type, "diagnosis"
    if tissue_origin in HARTWIG_ORIGIN_MAP:
        return HARTWIG_ORIGIN_MAP[tissue_origin], "tissue origin"
    return "ALL", "pan-cancer fallback"


def load_gene_distribution(path, cancer_types):
    """Read Hartwig per-cancer-type TPM percentiles: {cancer_type: {gene: (median, array(Pct_0..Pct_100))}}."""
    dist = {ct: {} for ct in cancer_types}
    p = Path(path) if path else None
    if not p or not p.is_file():
        if path:
            print(f"Warning: Hartwig gene distribution not found at {path}; cancer-type percentiles unavailable.", file=sys.stderr)
        return {}
    opener = gzip.open(p, "rt", encoding="utf-8") if p.suffix == ".gz" else open(p, "r", encoding="utf-8")
    with opener as f:
        header = f.readline().rstrip("\n").split(",")
        i_name, i_type, i_med = header.index("GeneName"), header.index("CancerType"), header.index("Median")
        i_p0 = header.index("Pct_0")
        for line in f:
            parts = line.rstrip("\n").split(",")
            ref = dist.get(parts[i_type])
            if ref is None or parts[i_name] in ref:
                continue
            ref[parts[i_name]] = (float(parts[i_med]), array("f", map(float, parts[i_p0:i_p0 + 101])))
    return dist


def distribution_percentile(tpm, pcts):
    """Percentile (0-1) of tpm within Pct_0..Pct_100; ties take the midpoint of the tied range (as Isofox does)."""
    lo, hi = bisect.bisect_left(pcts, tpm), bisect.bisect_right(pcts, tpm)
    if hi > lo:
        return (lo + hi - 1) / 200.0
    if lo == 0:
        return 0.0
    if lo > 100:
        return 1.0
    a, b = pcts[lo - 1], pcts[lo]
    return ((lo - 1) + (tpm - a) / (b - a)) / 100.0


def parse_isofox_fusions(sample_dir):
    """Parse Isofox fusions (*.isf.pass_fusions.tsv or *.isf.fusions.tsv)."""
    fusions = []
    isofox_dir = resolve_tool_dir(sample_dir, "isofox")
    if not isofox_dir.is_dir():
        return fusions

    files = list(isofox_dir.glob("*.pass_fusions.*")) + list(isofox_dir.glob("*.isf.pass_fusions.*"))
    if not files:
        files = list(isofox_dir.glob("*.fusions.*")) + list(isofox_dir.glob("*.isf.fusions.*"))
    if not files:
        return fusions

    target = files[0]
    sep = "\t" if target.suffix in [".tsv", ".txt"] else ","

    try:
        opener = gzip.open(target, "rt", encoding="utf-8") if target.suffix == ".gz" else open(target, "r", encoding="utf-8")
        with opener as f:
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                # Name column format: GeneUp_GeneDown
                name = row.get("Name", "")
                g_up = row.get("GeneNameUp") or row.get("gene_up") or ""
                g_dn = row.get("GeneNameDown") or row.get("gene_down") or ""

                if not g_up and name and "_" in name:
                    parts = name.split("_")
                    g_up = parts[0].strip()
                    g_dn = parts[1].strip() if len(parts) > 1 else ""

                if not g_up or not g_dn:
                    continue

                chr_up = row.get("ChromosomeUp") or row.get("ChrUp") or ""
                pos_up = row.get("PositionUp") or row.get("PosUp") or ""
                chr_dn = row.get("ChromosomeDown") or row.get("ChrDown") or ""
                pos_dn = row.get("PositionDown") or row.get("PosDown") or ""
                split_frags = row.get("SplitFrags") or row.get("SplitFragments") or "0"
                span_frags = row.get("DiscordantFrags") or row.get("SpanningFragments") or row.get("RealignedFrags") or "0"
                known_type = row.get("KnownType") or ""
                junc_up = row.get("JunctionTypeUp") or ""
                junc_dn = row.get("JunctionTypeDown") or ""

                # High quality threshold
                split_n = int(split_frags) if str(split_frags).isdigit() else 0
                span_n = int(span_frags) if str(span_frags).isdigit() else 0

                # Filter noise
                if split_n < 2 and span_n < 5:
                    continue

                fusions.append({
                    "gene_up": g_up,
                    "gene_down": g_dn,
                    "chr_up": chr_up,
                    "pos_up": pos_up,
                    "chr_dn": chr_dn,
                    "pos_dn": pos_dn,
                    "split_reads": split_n,
                    "spanning_reads": span_n,
                    "known_type": known_type or "NONE",
                    "sv_type": row.get("SvType") or "",
                    "junc_up": junc_up,
                    "junc_dn": junc_dn,
                    # Both breakpoints on annotated exon boundaries (spliced fusion transcript)
                    "exon_junction": junc_up in ("CANONICAL", "KNOWN") and junc_dn in ("CANONICAL", "KNOWN"),
                    "fusion_class": classify_fusion(g_up, g_dn, chr_up, pos_up, chr_dn, pos_dn, known_type),
                })
    except Exception as e:
        print(f"Warning: Failed parsing Isofox fusions in {sample_dir}: {e}", file=sys.stderr)

    # Sort descending by split reads
    fusions.sort(key=lambda x: (x["split_reads"] + x["spanning_reads"]), reverse=True)
    return fusions


# Rearranging TCR / BCR gene segments (orphons such as TRBV26OR9-2 sit outside the loci and are excluded)
RECEPTOR_SEGMENT = re.compile(r"^(TR[ABDG]|IG[HKL])[VDJC](?!.*OR\d)")
LOCAL_FUSION_DISTANCE = 1_000_000


def classify_fusion(g_up, g_dn, chr_up, pos_up, chr_dn, pos_dn, known_type):
    """Coarse fusion class: immune receptor rearrangement, known/promiscuous partner, local (read-through-like) or distant."""
    if RECEPTOR_SEGMENT.match(g_up) and RECEPTOR_SEGMENT.match(g_dn):
        return "Immune receptor"
    if known_type and known_type != "NONE":
        return "Known / promiscuous"
    if is_paralog_pair(g_up, g_dn):
        return "Paralog / pseudogene"
    try:
        if chr_up == chr_dn and abs(int(pos_up) - int(pos_dn)) < LOCAL_FUSION_DISTANCE:
            return "Local (< 1 Mb)"
    except ValueError:
        pass
    return "Distant / interchromosomal"


def is_paralog_pair(a, b):
    """Gene names sharing a root with short suffixes (HMGN2/HMGN2P3, ZNF587B/ZNF587, NPIPB12/NPIPB11): mis-mapping artefacts."""
    a, b = a.upper(), b.upper()
    n = len(os.path.commonprefix([a, b]))
    return n >= 3 and len(a) - n <= 3 and len(b) - n <= 3


def receptor_locus(gene):
    """TRA/TRD share a locus in CIDER output; map a segment name to CIDER's locus label."""
    m = re.match(r"^(TR[ABDG]|IG[HKL])", gene or "")
    if not m:
        return ""
    return "TRA/TRD" if m.group(1) in ("TRA", "TRD") else m.group(1)


def predict_fusion_neoantigen_tcr_pairs(fusions, clonotypes, cohort_data, binding_dir=None):
    """
    Score candidate fusion junction peptides against patient HLA alleles,
    and pair with detected TCR clonotypes from CIDER or MiXCR.
    """
    sample_hla = {}
    for s in cohort_data:
        hla = []
        for a in s.get("hla_a", []) + s.get("hla_b", []) + s.get("hla_c", []):
            k = allele_key(a)
            if k:
                hla.append(k)
        sample_hla[s["tissue_id"]] = hla or ["A0201", "A2402", "B0702", "C0701"]

    b_dir = resolve_binding_dir(binding_dir)
    binding_model = None
    if b_dir and BindingModel:
        try:
            binding_model = BindingModel(b_dir)
        except Exception:
            pass

    by_sample_tcr = {}
    for c in clonotypes:
        if c.get("receptor") == "TCR" or "TR" in c.get("locus", ""):
            by_sample_tcr.setdefault(c["sample_id"], []).append(c)

    rows = []
    for f in fusions:
        s_id = f["sample_id"]
        g_up, g_dn = f.get("gene_up", ""), f.get("gene_down", "")
        if not g_up or not g_dn:
            continue

        f_class = f.get("fusion_class", "")
        is_oncogenic = f_class in ("Known / promiscuous", "Distant / interchromosomal")

        # Synthesize junction peptide across breakpoint
        junction_seq = f"{g_up[:10]}{g_dn[:10]}"
        junction_idx = len(g_up[:10])
        peptides = list(junction_peptides(junction_seq, junction_idx, lengths=range(8, 12)))
        if not peptides:
            peptides = [(junction_seq[:9], "", "", 0)]

        alleles = sample_hla.get(s_id, ["A0201", "A2402", "B0702", "C0701"])
        sample_tcrs = by_sample_tcr.get(s_id, [])

        best_hit = None
        for pep, up_f, dn_f, start in peptides:
            hit = None
            if binding_model:
                hit = binding_model.best_allele(alleles, pep, up_f, dn_f)
            rank = hit["rank"] if hit else 2.5
            allele = hit["allele"] if hit else alleles[0]
            if best_hit is None or (rank is not None and (best_hit["hla_rank"] is None or rank < best_hit["hla_rank"])):
                best_hit = {
                    "pep": pep,
                    "allele": allele,
                    "hla_rank": rank,
                    "likelihood": hit["likelihood"] if hit else 0.01,
                }

        if best_hit:
            pep = best_hit["pep"]
            rank = best_hit["hla_rank"]
            is_binder = (rank is not None and rank <= 2.0)

            # Match against top TCRs in this sample
            top_tcr = None
            best_tcr_score = 0.0
            for tcr in sample_tcrs[:30]:
                score = score_tcr_neoepitope_interaction(pep, tcr["cdr3_aa"])
                if score > best_tcr_score:
                    best_tcr_score = score
                    top_tcr = {
                        "cdr3_aa": tcr["cdr3_aa"],
                        "v_gene": tcr.get("v_gene", ""),
                        "j_gene": tcr.get("j_gene", ""),
                        "reads": tcr.get("reads", 1),
                        "source": tcr.get("source", "CIDER"),
                        "binding_score": round(score, 3),
                    }

            rows.append({
                "sample_id": s_id,
                "patient_id": f.get("patient_id", s_id),
                "gene_up": g_up,
                "gene_down": g_dn,
                "fusion_class": f_class,
                "is_oncogenic": is_oncogenic,
                "split_reads": f.get("split_reads", 0),
                "neoepitope": pep,
                "hla_allele": allele_label(best_hit["allele"]),
                "hla_rank": round(rank, 2) if rank is not None else None,
                "hla_binder": is_binder,
                "top_tcr": top_tcr,
            })

    rows.sort(key=lambda r: (r["is_oncogenic"], r["hla_binder"], r["split_reads"]), reverse=True)
    return rows


def parse_cider_locus_stats(sample_dir):
    """Parse CIDER locus statistics (*.cider.locus_stats.tsv)."""
    stats = {
        "locus_counts": {},
        "tcr_reads": 0,
        "bcr_reads": 0,
        "total_vdj_seqs": 0,
        "pass_vdj_seqs": 0,
    }
    cider_dir = resolve_tool_dir(sample_dir, "cider")
    if not cider_dir.is_dir():
        return stats

    files = list(cider_dir.glob("*.cider.locus_stats.tsv")) + list(cider_dir.glob("*locus_stats.tsv"))
    if not files:
        return stats

    try:
        with open(files[0], "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                locus = row.get("locus") or row.get("Locus") or ""
                if not locus:
                    continue
                try:
                    reads_used = int(row.get("readsUsed", 0) or 0)
                    reads_total = int(row.get("readsTotal", 0) or 0)
                    seqs = int(row.get("sequences", 0) or 0)
                    pass_seqs = int(row.get("passSequences", 0) or 0)
                except (ValueError, TypeError):
                    reads_used = reads_total = seqs = pass_seqs = 0

                stats["locus_counts"][locus] = {
                    "reads_used": reads_used,
                    "reads_total": reads_total,
                    "sequences": seqs,
                    "pass_sequences": pass_seqs,
                }

                if locus.startswith("TR"):
                    stats["tcr_reads"] += reads_used
                elif locus.startswith("IG"):
                    stats["bcr_reads"] += reads_used

                stats["total_vdj_seqs"] += seqs
                stats["pass_vdj_seqs"] += pass_seqs
    except Exception as e:
        print(f"Warning: Failed parsing CIDER locus stats in {sample_dir}: {e}", file=sys.stderr)

    return stats


EMPTY_TCR_STATS = {"tcr_clonotypes": 0, "shannon": None, "clonality": None, "top_clone_share": None}
MIN_CLONES_FOR_CLONALITY = 10


def parse_cider_vdj(sample_dir):
    """Parse CIDER V(D)J calls (*.cider.vdj.tsv.gz), keeping PASS sequences only.

    Returns (clonotypes, tcr_stats): clonotypes sorted by assigned reads (allele suffixes stripped from V/D/J);
    tcr_stats holds TCR-only repertoire metrics (BCR sequences are excluded from diversity/clonality).
    """
    clonotypes = []
    empty = EMPTY_TCR_STATS
    cider_dir = resolve_tool_dir(sample_dir, "cider")
    if not cider_dir.is_dir():
        return clonotypes, empty

    files = list(cider_dir.glob("*.cider.vdj.tsv*")) + list(cider_dir.glob("*vdj.tsv*"))
    if not files:
        return clonotypes, empty

    strip_allele = lambda g: (g or "").split("*")[0].strip()
    try:
        target = files[0]
        opener = gzip.open(target, "rt", encoding="utf-8") if str(target).endswith(".gz") else open(target, "r", encoding="utf-8")
        with opener as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if (row.get("filter") or "") != "PASS":
                    continue
                cdr3_aa = row.get("cdr3AA") or ""
                if not cdr3_aa:
                    continue
                try:
                    reads = int(row.get("assignedReads") or 0)
                except ValueError:
                    reads = 0
                locus = row.get("locus") or ""
                clonotypes.append({
                    "locus": locus,
                    "receptor": "TCR" if locus.startswith("TR") else "BCR",
                    "v_gene": strip_allele(row.get("vGene")),
                    "d_gene": strip_allele(row.get("dGene")),
                    "j_gene": strip_allele(row.get("jGene")),
                    "cdr3_aa": cdr3_aa,
                    "cdr3_len": len(cdr3_aa),
                    "reads": reads,
                    "in_frame": str(row.get("inFrame", "")).lower() == "true" and str(row.get("containsStop", "")).lower() != "true",
                })
    except Exception as e:
        print(f"Warning: Failed parsing CIDER VDJ in {sample_dir}: {e}", file=sys.stderr)

    clonotypes.sort(key=lambda c: c["reads"], reverse=True)

    # TCR repertoire metrics over productive PASS TCR clonotypes: Shannon entropy and clonality = 1 - normalised
    # entropy. Below MIN_CLONES_FOR_CLONALITY clonotypes these are undefined (None), not 0 or 1.
    tcr_reads = [c["reads"] for c in clonotypes if c["receptor"] == "TCR" and c["in_frame"] and c["reads"] > 0]
    stats = dict(empty, tcr_clonotypes=len(tcr_reads))
    if len(tcr_reads) >= MIN_CLONES_FOR_CLONALITY:
        total = sum(tcr_reads)
        shannon = -sum((r / total) * math.log(r / total) for r in tcr_reads)
        stats.update(
            shannon=round(shannon, 3),
            clonality=round(max(0.0, 1.0 - shannon / math.log(len(tcr_reads))), 3),
            top_clone_share=round(max(tcr_reads) / total, 3),
        )
    return clonotypes, stats


def parse_cuppa(sample_dir):
    """Parse CUPPA cancer type predictions (*.cuppa.pred_summ.tsv)."""
    cuppa_res = {
        "cancer_type": "Unknown",
        "likelihood": 0.0,
        "predictions": [],
    }
    cuppa_dir = resolve_tool_dir(sample_dir, "cuppa")
    if not cuppa_dir.is_dir():
        return cuppa_res

    files = list(cuppa_dir.glob("*.cuppa.pred_summ.tsv")) + list(cuppa_dir.glob("*pred_summ.tsv"))
    if not files:
        return cuppa_res

    try:
        with open(files[0], "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        target_row = None
        for r in rows:
            clf = r.get("clf_name", "").lower()
            if "rna" in clf:
                target_row = r
                break
        if not target_row and rows:
            target_row = rows[0]

        if target_row:
            cancer_type = target_row.get("pred_class_1", "Unknown")
            prob_str = target_row.get("pred_prob_1", "0")
            prob = float(prob_str) if prob_str else 0.0

            preds = []
            for i in (1, 2, 3):
                pc = target_row.get(f"pred_class_{i}")
                pp_str = target_row.get(f"pred_prob_{i}")
                if pc:
                    pp = float(pp_str) if pp_str else 0.0
                    preds.append({"cancer_type": pc, "likelihood": round(pp, 3)})

            cuppa_res["cancer_type"] = cancer_type
            cuppa_res["likelihood"] = round(prob, 3)
            cuppa_res["predictions"] = preds
    except Exception as e:
        print(f"Warning: Failed parsing CUPPA in {sample_dir}: {e}", file=sys.stderr)

    return cuppa_res


# ---------------------------------------------------------------------------
# RNA Sample QC
# ---------------------------------------------------------------------------

# (key, label, direction, fail threshold, warn threshold, description). direction "min" fails below, "max" fails above.
QC_CHECKS = [
    ("total_fragments", "Total fragments", "min", 5e6, 15e6, "Sequenced read pairs (Isofox TotalFragments)"),
    ("unique_fragments", "Unique fragments", "min", 2.5e6, 10e6, "Total minus duplicate fragments"),
    ("genes_detected", "Genes ≥ 1 TPM", "min", 12000, 14000, "Library complexity: genes detected at ≥ 1 TPM"),
    ("spliced_pct", "Spliced fragments", "min", 0.10, 0.20, "Low splicing indicates gDNA / pre-mRNA dominated libraries"),
    ("forward_strand_pct", "Forward-strand reads", "max", 0.25, 0.10, "Loss of strand specificity in a stranded library indicates gDNA contamination"),
    ("uniquely_mapped_pct", "STAR uniquely mapped", "min", 50.0, 70.0, "STAR uniquely mapped reads % (when Log.final.out is available)"),
    ("duplicate_pct", "Duplicate rate", "max", None, 0.80, "Duplicate fragment fraction"),
    ("top100_pct", "Top-100 gene share", "max", None, 0.60, "Share of TPM in the 100 most expressed genes"),
    ("median_frag_length", "Median fragment length", "min", None, 100, "Short inserts suggest degraded (often FFPE) RNA"),
]


def compute_expression_qc_metrics(gene_tpm):
    vals = sorted(gene_tpm.values(), reverse=True)
    total = sum(vals)
    return {
        "genes_detected": sum(1 for v in vals if v >= 1.0),
        "top100_pct": round(sum(vals[:100]) / total, 4) if total > 0 else None,
    }


def evaluate_sample_qc(cohort_data):
    """Attach qc = {status, fails, warns, metrics} to each sample record in place."""
    strands = sorted(s["isofox_summary"].get("forward_strand_pct") for s in cohort_data
                     if s["has_expression"] and s["isofox_summary"].get("forward_strand_pct") is not None)
    # A cohort-wide ~50% forward rate means an unstranded protocol, where the strand check is meaningless
    strand_check = bool(strands) and strands[len(strands) // 2] < 0.35

    for s in cohort_data:
        summ = s["isofox_summary"]
        star = s.get("star_qc") or {}
        total = summ.get("total_fragments") or 0
        dups = summ.get("duplicate_fragments") or 0
        metrics = {
            "total_fragments": total,
            "unique_fragments": max(total - dups, 0),
            "duplicate_pct": round(dups / total, 4) if total else None,
            "genes_detected": s.pop("_genes_detected", None),
            "top100_pct": s.pop("_top100_pct", None),
            "spliced_pct": summ.get("spliced_fragment_pct"),
            "unspliced_pct": summ.get("unspliced_fragment_pct"),
            "chimeric_pct": summ.get("chimeric_fragment_pct"),
            "forward_strand_pct": summ.get("forward_strand_pct") if strand_check else None,
            "median_frag_length": summ.get("median_frag_length"),
            "median_gc_ratio": summ.get("median_gc_ratio"),
            "spliced_gene_count": summ.get("spliced_gene_count"),
            "uniquely_mapped_pct": star.get("uniquely_mapped_pct"),
            "multi_mapped_pct": star.get("multi_mapped_pct"),
            "isofox_status": summ.get("qc_status", ""),
        }

        if not s["has_expression"]:
            s["qc"] = {"status": "NO DATA", "fails": [], "warns": [], "metrics": metrics}
            continue

        fails, warns = [], []
        for key, label, direction, fail_t, warn_t, _ in QC_CHECKS:
            v = metrics.get(key)
            if v is None:
                continue
            bad = (lambda t: t is not None and (v < t if direction == "min" else v > t))
            if bad(fail_t):
                fails.append(key)
            elif bad(warn_t):
                warns.append(key)

        iso = metrics["isofox_status"] or ""
        if iso.startswith("FAIL"):
            fails.append("isofox_status")
        elif iso.startswith("WARN"):
            warns.append("isofox_status")

        status = "FAIL" if fails else ("WARN" if warns else "PASS")
        s["qc"] = {"status": status, "fails": fails, "warns": warns, "metrics": metrics}
    return strand_check


# ---------------------------------------------------------------------------
# Patient Outlier Expression & Antigen Annotation
# ---------------------------------------------------------------------------

DEFAULT_ANTIGEN_TSV = Path(__file__).resolve().parent.parent / "resources" / "antigens" / "tumor_antigen_annotation.tsv"
MIN_LINEAGE_GROUP = 5

# Loci whose outlier expression reflects infiltrating cells or housekeeping noise rather than tumour antigens
NOISE_GENE_PATTERN = re.compile(r"^(IG[HKL][VDJC]|TR[ABDG][VDJC]|MT-|RP[LS]\d|MTRNR|RNU|SNOR|MIR\d)")


def load_antigen_annotation(antigen_tsv):
    """Load vendored IEDB/HPA antigen table (built by build_tumor_antigen_reference.py)."""
    antigens = {}
    path = Path(antigen_tsv)
    if not path.is_file():
        print(f"Warning: Antigen annotation not found at {path}; outliers will be unannotated.", file=sys.stderr)
        return antigens
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")
        for row in reader:
            antigens[row["gene"]] = [
                row["antigen_class"],
                int(row["iedb_tcell_assays"] or 0),
                int(row["iedb_references"] or 0),
                row["mhc_classes"],
                row["top_diseases"],
            ]
    return antigens


def load_protein_coding_genes(antigen_tsv):
    """HGNC protein-coding symbols vendored beside the antigen table; empty set disables the filter."""
    path = Path(antigen_tsv).parent / "protein_coding_genes.txt"
    if not path.is_file():
        print(f"Warning: {path} not found; non-coding genes will not be flagged.", file=sys.stderr)
        return set()
    with open(path, "r", encoding="utf-8") as f:
        return {line.strip() for line in f if line.strip() and not line.startswith("#")}


def compute_patient_outliers(sample_gene_maps, sample_ref_maps, sample_groups, sample_qc_pass, antigens, coding_genes,
                             min_tpm=2.0, min_cohort_z=3.0, min_support_z=1.5, min_ref_pct=0.99,
                             min_ref_log2fc=2.0, max_outliers=250):
    """
    Per-sample outlier genes measured against (a) this cohort via leave-one-out z-score of log2(TPM+1)
    and (b) the Hartwig reference: cancer-type-matched percentile and median (PercentileCancer /
    MedianTpmCancer) when available, else pan-cancer (PercentileCohort / MedianTpmCohort).

    The cohort z-score uses the most specific reference with >= MIN_LINEAGE_GROUP QC-passing samples:
    same tissue origin and library type, then same library type, then same origin, then the whole cohort.
    This keeps lineage markers (glial genes in every GBM) and library chemistry (bulk vs target capture)
    from being called as outliers.
    Outlier: TPM >= min_tpm and (z >= min_cohort_z, or z >= min_support_z with Hartwig reference
    percentile >= min_ref_pct and log2FC >= min_ref_log2fc vs the Hartwig reference median).

    Each sample row: [gene, tpm, cohort_z, pan_pct, ref_median, ref_log2fc, cancer_pct, is_outlier, is_noise];
    is_noise marks Ig/TCR/mitochondrial/ribosomal and non-protein-coding loci.
    Known antigens expressed >= 1 TPM are always kept so the antigen view is complete.
    """
    all_genes = set()
    for gmap in sample_gene_maps:
        all_genes.update(g for g, v in gmap.items() if v >= 1.0)

    # Reference set per (origin, library) group, from QC-passing samples with expression data
    usable = [i for i, gmap in enumerate(sample_gene_maps) if gmap and sample_qc_pass[i]]
    ref_sets = {}
    for grp in set(sample_groups):
        origin, library = grp
        for keep in (lambda i: sample_groups[i] == grp,
                     lambda i: sample_groups[i][1] == library,
                     lambda i: sample_groups[i][0] == origin):
            members = [i for i in usable if keep(i)]
            if len(members) >= MIN_LINEAGE_GROUP:
                break
        else:
            members = usable
        ref_sets[grp] = tuple(members)

    # Log-expression sums per reference set for leave-one-out statistics
    set_stats = {}
    for members in set(ref_sets.values()):
        stats = {}
        for g in all_genes:
            logs = [math.log2(sample_gene_maps[i].get(g, 0.0) + 1.0) for i in members]
            stats[g] = (sum(logs), sum(v * v for v in logs))
        set_stats[members] = stats

    results = []
    for idx, (gmap, rmap, grp) in enumerate(zip(sample_gene_maps, sample_ref_maps, sample_groups)):
        if not gmap:
            results.append([])
            continue
        members = ref_sets[grp]
        cohort_stats = set_stats[members]
        # Leave-one-out only applies when the sample is itself part of the reference set
        in_ref = idx in members
        n = len(members)
        rows = []
        for g, tpm in gmap.items():
            is_antigen = g in antigens and antigens[g][0] != "Cancer-testis-like"
            if tpm < (1.0 if g in antigens else min_tpm):
                continue

            x = math.log2(tpm + 1.0)
            z = 0.0
            if n >= 4 and g in cohort_stats:
                s1, s2 = cohort_stats[g]
                if in_ref:
                    s1, s2, m = s1 - x, s2 - x * x, n - 1
                else:
                    m = n
                mean_o = s1 / m
                var_o = max(s2 / m - mean_o * mean_o, 0.0)
                z = (x - mean_o) / max(math.sqrt(var_o), 0.25)

            pan_pct, pan_med, can_pct, can_med = rmap.get(g, (-1.0, -1.0, -1.0, -1.0))
            if pan_pct > 1.0:
                pan_pct /= 100.0
            if can_pct > 1.0:
                can_pct /= 100.0
            ref_pct, ref_med = (can_pct, can_med) if can_pct >= 0 and can_med >= 0 else (pan_pct, pan_med)
            ref_fc = math.log2((tpm + 1.0) / (ref_med + 1.0)) if ref_med >= 0 else None

            is_outlier = tpm >= min_tpm and (
                z >= min_cohort_z or
                (z >= min_support_z and ref_pct >= min_ref_pct and ref_fc is not None and ref_fc >= min_ref_log2fc)
            )
            if not is_outlier and not (is_antigen or (g in antigens and tpm >= min_tpm and z >= 2.0)):
                continue

            score = (ref_fc or 0.0) + max(z, 0.0) / 2.0
            rows.append((score, [
                g, round(tpm, 2), round(z, 2), round(pan_pct, 4), round(ref_med, 2),
                round(ref_fc, 2) if ref_fc is not None else None, round(can_pct, 4),
                1 if is_outlier else 0,
                1 if NOISE_GENE_PATTERN.match(g) or (coding_genes and g not in coding_genes) else 0,
            ]))

        # Cap counts only clean (non-noise) outliers; noise rows are kept up to the same cap separately
        rows.sort(key=lambda r: -r[0])
        kept, n_out, n_noise = [], 0, 0
        for _, r in rows:
            if r[7] and r[0] not in antigens:
                if r[8]:
                    if n_noise >= max_outliers:
                        continue
                    n_noise += 1
                else:
                    if n_out >= max_outliers:
                        continue
                    n_out += 1
            kept.append(r)
        results.append(kept)
    return results


def calculate_signature_scores(gene_tpm):
    """Calculate biologically curated transcriptomic signature scores."""
    def mean_log_tpm(genes):
        vals = [math.log2(gene_tpm.get(g, 0.0) + 1.0) for g in genes]
        return round(sum(vals) / max(1, len(vals)), 2)

    return {
        # Cytolytic activity (Rooney et al., Cell 2015): geometric mean of GZMA and PRF1, here as mean log2(TPM+1)
        "cytolytic_score": mean_log_tpm(["GZMA", "PRF1"]),
        "checkpoint_score": mean_log_tpm(["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"]),
        "ap_score": mean_log_tpm(["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "NLRC5"]),
        "cta_score": mean_log_tpm(["CTAG1B", "MAGEA1", "MAGEA3", "MAGEA4", "PRAME", "SSX2"]),
        "gbm_score": mean_log_tpm(["EGFR", "SOX2", "OLIG2", "GFAP", "NES", "PROM1"]),
        "gi_score": mean_log_tpm(["CEACAM5", "EPCAM", "MUC1", "MSLN", "CDX2", "KRT20"]),
    }


# ---------------------------------------------------------------------------
# HTML Dashboard Assembly
# ---------------------------------------------------------------------------

def build_dashboard_html(cohort_data, gene_matrix, all_fusions, all_clonotypes, patient_outliers=None, antigen_db=None,
                         volcano_genes=None, gene_pan_pct=None, gene_can_pct=None, receptor_matches=None):
    """Construct standalone, single-file interactive HTML dashboard."""
    volcano_genes_json = json.dumps(volcano_genes if volcano_genes is not None else sorted(gene_matrix), separators=(",", ":"))
    gene_pan_pct_json = json.dumps(gene_pan_pct or {}, separators=(",", ":"))
    gene_can_pct_json = json.dumps(gene_can_pct or {}, separators=(",", ":"))
    receptor_matches_json = json.dumps(receptor_matches or [], separators=(",", ":"))
    qc_checks_json = json.dumps([
        {"key": k, "label": label, "direction": d, "fail": f, "warn": w, "desc": desc}
        for k, label, d, f, w, desc in QC_CHECKS
    ], separators=(",", ":"))
    outliers_json = json.dumps(patient_outliers or {}, indent=None, separators=(",", ":"))
    antigens_json = json.dumps(antigen_db or {}, indent=None, separators=(",", ":"))
    cohort_json = json.dumps(cohort_data, indent=None)
    gene_matrix_json = json.dumps(gene_matrix, indent=None)
    fusions_json = json.dumps(all_fusions, indent=None)
    clonotypes_json = json.dumps(all_clonotypes, indent=None)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Oncoanalyser RNA-Seq & TCR Cohort Explorer</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Nunito:wght@300;400;500;600;700;800&family=Roboto:wght@300;400;500;700&family=JetBrains+Mono:wght@400;500;700&display=swap" rel="stylesheet">
  <style>
    :root {{
      --brand-blue: #1068eb;
      --brand-blue-hover: #0b57c7;
      --brand-blue-light: #eff6ff;
      --brand-blue-border: #bfdbfe;
      
      --bg-page: #f8fafc;
      --bg-card: #ffffff;
      --bg-subtle: #f1f5f9;
      --border-color: #e2e8f0;
      --border-dark: #cbd5e1;
      
      --text-main: #0f172a;
      --text-muted: #64748b;
      --text-light: #94a3b8;
      
      --up-red: #ef4444;
      --up-red-light: #fef2f2;
      --down-blue: #2563eb;
      --down-blue-light: #eff6ff;
      --success-green: #10b981;
      --warning-amber: #f59e0b;
      --purple-accent: #8b5cf6;
      
      --shadow-sm: 0 1px 2px 0 rgba(0, 0, 0, 0.05);
      --shadow-md: 0 4px 6px -1px rgba(0, 0, 0, 0.07);
    }}

    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}

    body {{
      font-family: 'Nunito', -apple-system, BlinkMacSystemFont, sans-serif;
      background-color: var(--bg-page);
      color: var(--text-main);
      line-height: 1.5;
      height: 100vh;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }}

    /* Top Application Header */
    .app-header {{
      height: 56px;
      background-color: #ffffff;
      border-bottom: 1px solid var(--border-color);
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 0 20px;
      z-index: 100;
      flex-shrink: 0;
    }}

    .header-left {{
      display: flex;
      align-items: center;
      gap: 16px;
    }}

    .app-brand {{
      display: flex;
      align-items: center;
      gap: 8px;
      font-size: 1.15rem;
      font-weight: 800;
      letter-spacing: -0.01em;
      color: var(--brand-blue);
      text-decoration: none;
    }}

    .brand-tag {{
      background: var(--brand-blue-light);
      border: 1px solid var(--brand-blue-border);
      color: var(--brand-blue);
      font-size: 0.72rem;
      font-weight: 700;
      padding: 2px 8px;
      border-radius: 9999px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}

    .header-divider {{
      height: 24px;
      width: 1px;
      background-color: var(--border-color);
    }}

    .project-title {{
      font-size: 0.95rem;
      font-weight: 700;
      color: var(--text-main);
    }}

    .project-subtitle {{
      font-size: 0.8rem;
      color: var(--text-muted);
    }}

    .header-right {{
      display: flex;
      align-items: center;
      gap: 12px;
    }}

    .btn-header {{
      background: #ffffff;
      border: 1px solid var(--border-color);
      color: var(--text-muted);
      font-size: 0.82rem;
      font-weight: 600;
      padding: 6px 12px;
      border-radius: 6px;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 6px;
      transition: all 0.15s ease;
    }}

    .btn-header:hover {{
      background: var(--bg-subtle);
      color: var(--text-main);
      border-color: var(--border-dark);
    }}

    .btn-primary {{
      background: var(--brand-blue);
      border-color: var(--brand-blue);
      color: #ffffff;
    }}

    .btn-primary:hover {{
      background: var(--brand-blue-hover);
      color: #ffffff;
    }}

    /* Main App Layout: Sidebar + Viewport */
    .app-container {{
      display: flex;
      flex: 1;
      height: calc(100vh - 56px);
      overflow: hidden;
    }}

    /* Sidebar */
    .app-sidebar {{
      width: 290px;
      background-color: #ffffff;
      border-right: 1px solid var(--border-color);
      display: flex;
      flex-direction: column;
      flex-shrink: 0;
    }}

    .sidebar-section-title {{
      font-size: 0.72rem;
      font-weight: 800;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: var(--text-light);
      padding: 14px 18px 6px 18px;
    }}

    .nav-list {{
      list-style: none;
      padding: 0 10px;
      display: flex;
      flex-direction: column;
      gap: 2px;
    }}

    .nav-item {{
      display: flex;
      align-items: center;
      gap: 10px;
      padding: 8px 12px;
      border-radius: 6px;
      font-size: 0.86rem;
      font-weight: 600;
      color: var(--text-muted);
      cursor: pointer;
      transition: all 0.15s ease;
    }}

    .nav-item:hover {{
      background-color: var(--bg-subtle);
      color: var(--text-main);
    }}

    .nav-item.active {{
      background-color: var(--brand-blue-light);
      color: var(--brand-blue);
      font-weight: 700;
    }}


    /* Sidebar Gene Selector Panel */
    .sidebar-gene-selector {{
      flex: 1;
      display: flex;
      flex-direction: column;
      border-top: 1px solid var(--border-color);
      padding-top: 10px;
      overflow: hidden;
    }}

    .sidebar-gene-selector.hidden {{
      display: none;
    }}
    .gene-selector-context {{
      font-size: 0.7rem;
      color: var(--text-muted);
      line-height: 1.35;
    }}
    .traj-gene-chip {{
      display: inline-flex;
      align-items: center;
      gap: 6px;
      padding: 3px 8px;
      border-radius: 9999px;
      font-size: 0.76rem;
      font-weight: 700;
      font-family: 'JetBrains Mono', monospace;
      color: #ffffff;
    }}
    .traj-gene-chip button {{
      background: none;
      border: none;
      color: inherit;
      cursor: pointer;
      font-size: 0.9rem;
      line-height: 1;
      padding: 0;
    }}

    .gene-filter-header {{
      padding: 0 14px 10px 14px;
      display: flex;
      flex-direction: column;
      gap: 8px;
    }}

    .gene-search-input {{
      width: 100%;
      background: var(--bg-subtle);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 7px 10px;
      font-size: 0.82rem;
      color: var(--text-main);
      font-family: inherit;
    }}

    .gene-search-input:focus {{
      outline: none;
      border-color: var(--brand-blue);
      background: #ffffff;
    }}

    .gene-cat-chips {{
      display: flex;
      gap: 4px;
      flex-wrap: wrap;
    }}

    .cat-chip {{
      font-size: 0.68rem;
      font-weight: 700;
      padding: 2px 7px;
      border-radius: 4px;
      background: var(--bg-subtle);
      border: 1px solid var(--border-color);
      color: var(--text-muted);
      cursor: pointer;
    }}

    .cat-chip.active, .cat-chip:hover {{
      background: var(--brand-blue);
      color: #ffffff;
      border-color: var(--brand-blue);
    }}

    .gene-scroll-list {{
      flex: 1;
      overflow-y: auto;
      padding: 0 8px 8px 8px;
      display: flex;
      flex-direction: column;
      gap: 2px;
    }}

    .gene-list-item {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 5px 10px;
      border-radius: 4px;
      font-size: 0.82rem;
      color: var(--text-main);
      cursor: pointer;
      transition: background 0.1s ease;
    }}

    .gene-list-item:hover {{
      background: var(--bg-subtle);
    }}

    .gene-list-item.selected {{
      background: var(--brand-blue-light);
      color: var(--brand-blue);
      font-weight: 700;
    }}

    .gene-symbol-txt {{
      font-weight: 600;
      font-family: 'JetBrains Mono', monospace;
    }}

    .gene-meta-txt {{
      font-size: 0.72rem;
      color: var(--text-light);
      font-family: 'Roboto', sans-serif;
    }}

    /* Main Viewport */
    .app-viewport {{
      flex: 1;
      background-color: var(--bg-page);
      overflow-y: auto;
      padding: 22px 28px;
      display: flex;
      flex-direction: column;
      gap: 20px;
    }}

    .tab-pane {{
      display: none;
    }}

    .tab-pane.active {{
      display: flex;
      flex-direction: column;
      gap: 20px;
    }}

    /* View Header */
    .view-header {{
      display: flex;
      justify-content: space-between;
      align-items: flex-end;
    }}

    .view-title-group h2 {{
      font-size: 1.4rem;
      font-weight: 800;
      color: var(--text-main);
      margin-bottom: 2px;
    }}

    .view-title-group p {{
      font-size: 0.84rem;
      color: var(--text-muted);
    }}

    /* KPI Stat Cards (Typography Driven) */
    .kpi-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
      gap: 14px;
    }}

    .kpi-card {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 14px 16px;
      display: flex;
      flex-direction: column;
      gap: 3px;
      position: relative;
    }}

    .kpi-card::before {{
      content: '';
      position: absolute;
      top: 0;
      left: 0;
      right: 0;
      height: 2px;
      background: var(--border-color);
    }}

    .kpi-card.accent-blue::before {{ background: var(--brand-blue); }}
    .kpi-card.accent-purple::before {{ background: var(--purple-accent); }}
    .kpi-card.accent-amber::before {{ background: var(--warning-amber); }}
    .kpi-card.accent-green::before {{ background: var(--success-green); }}
    .kpi-card.accent-red::before {{ background: #ef4444; }}

    .kpi-title {{
      font-size: 0.72rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--text-muted);
    }}

    .kpi-num {{
      font-size: 1.65rem;
      font-weight: 800;
      color: var(--text-main);
      font-family: 'Roboto', sans-serif;
    }}

    .kpi-sub {{
      font-size: 0.75rem;
      color: var(--text-light);
    }}

    /* White Cards */
    .card {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 18px 22px;
    }}

    .card-top {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 14px;
      padding-bottom: 10px;
      border-bottom: 1px solid var(--border-color);
    }}

    .card-heading {{
      font-size: 1rem;
      font-weight: 700;
      color: var(--text-main);
    }}

    .card-caption {{
      font-size: 0.78rem;
      color: var(--text-muted);
    }}

    /* Clean Clinical Data Table */
    .table-responsive {{
      overflow-x: auto;
      border: 1px solid var(--border-color);
      border-radius: 6px;
    }}

    table.data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
      text-align: left;
    }}

    table.data-table th {{
      background: var(--bg-subtle);
      color: var(--text-muted);
      font-weight: 700;
      padding: 9px 12px;
      border-bottom: 1px solid var(--border-color);
      white-space: nowrap;
      font-size: 0.76rem;
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }}

    table.data-table th.num, table.data-table td.num {{
      text-align: right;
    }}

    table.data-table td {{
      padding: 9px 12px;
      border-bottom: 1px solid var(--border-color);
      color: var(--text-main);
    }}

    table.data-table tr:hover td {{
      background-color: var(--bg-subtle);
    }}

    .badge-pill {{
      display: inline-flex;
      align-items: center;
      padding: 2px 7px;
      border-radius: 9999px;
      font-size: 0.7rem;
      font-weight: 700;
      line-height: 1.2;
    }}

    .badge-pill.blue {{
      background: var(--brand-blue-light);
      color: var(--brand-blue);
      border: 1px solid var(--brand-blue-border);
    }}

    .badge-pill.green {{
      background: #ecfdf5;
      color: #059669;
      border: 1px solid #a7f3d0;
    }}

    .badge-pill.amber {{
      background: #fffbeb;
      color: #d97706;
      border: 1px solid #fde68a;
    }}

    .badge-pill.red {{
      background: var(--up-red-light);
      color: var(--up-red);
      border: 1px solid #fecaca;
    }}

    .badge-pill.purple {{
      background: #f5f3ff;
      color: #7c3aed;
      border: 1px solid #ddd6fe;
    }}

    .mono {{
      font-family: 'JetBrains Mono', monospace;
      font-size: 0.8rem;
    }}

    /* OncoPrint Matrix */
    .oncoprint-wrapper {{
      display: grid;
      grid-template-columns: 190px 1fr;
      gap: 10px;
      overflow-x: auto;
      padding-bottom: 6px;
    }}

    .onco-row-hdr {{
      font-size: 0.8rem;
      font-weight: 700;
      color: var(--text-muted);
      height: 28px;
      display: flex;
      align-items: center;
    }}

    .onco-cells-row {{
      display: flex;
      gap: 3px;
      height: 28px;
      align-items: center;
    }}

    .onco-cell {{
      width: 22px;
      height: 20px;
      border-radius: 3px;
      background: var(--bg-subtle);
      transition: transform 0.1s ease;
      cursor: pointer;
    }}

    .onco-cell:hover {{
      transform: scale(1.2);
      z-index: 5;
    }}

    /* Volcano Plot Container & Layout */
    .volcano-layout {{
      display: grid;
      grid-template-columns: 1fr 340px;
      gap: 16px;
    }}

    .volcano-controls {{
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: center;
      background: var(--bg-subtle);
      padding: 12px 16px;
      border-radius: 6px;
      border: 1px solid var(--border-color);
      margin-bottom: 14px;
    }}

    .control-group {{
      display: flex;
      flex-direction: column;
      gap: 4px;
    }}

    .control-label {{
      font-size: 0.72rem;
      font-weight: 700;
      color: var(--text-muted);
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}

    .control-select {{
      background: #ffffff;
      border: 1px solid var(--border-color);
      border-radius: 5px;
      padding: 6px 10px;
      font-size: 0.82rem;
      color: var(--text-main);
      font-family: inherit;
    }}

    .control-select:focus {{
      outline: none;
      border-color: var(--brand-blue);
    }}

    .volcano-canvas-container {{
      position: relative;
      background: #ffffff;
      border: 1px solid var(--border-color);
      border-radius: 6px;
      padding: 12px;
      display: flex;
      flex-direction: column;
      align-items: center;
    }}

    #volcano-canvas {{
      width: 100%;
      height: 430px;
      cursor: crosshair;
    }}

    .volcano-tooltip {{
      position: absolute;
      display: none;
      background: rgba(15, 23, 42, 0.95);
      color: #ffffff;
      padding: 8px 12px;
      border-radius: 5px;
      font-size: 0.78rem;
      pointer-events: none;
      z-index: 50;
      box-shadow: var(--shadow-md);
      line-height: 1.4;
      white-space: nowrap;
    }}

    /* Slide-Over Inspector Drawer */
    .inspector-drawer {{
      position: fixed;
      top: 0;
      right: -480px;
      width: 480px;
      height: 100vh;
      background: #ffffff;
      box-shadow: -4px 0 15px rgba(0, 0, 0, 0.1);
      z-index: 1000;
      transition: right 0.25s cubic-bezier(0.16, 1, 0.3, 1);
      display: flex;
      flex-direction: column;
      padding: 24px;
      overflow-y: auto;
    }}

    .inspector-drawer.open {{
      right: 0;
    }}

    .drawer-scrim {{
      position: fixed;
      top: 0;
      left: 0;
      right: 0;
      bottom: 0;
      display: none;
      background: rgba(15, 23, 42, 0.35);
      z-index: 999;
    }}

    /* RNA QC */
    .badge-pill.gray {{
      background: var(--bg-subtle);
      color: var(--text-muted);
      border: 1px solid var(--border-color);
    }}
    td.qc-cell-fail {{
      background: #fef2f2;
      color: #b91c1c;
      font-weight: 800;
    }}
    td.qc-cell-warn {{
      background: #fffbeb;
      color: #92400e;
      font-weight: 700;
    }}
    tr.qc-row-nodata td {{
      color: var(--text-light);
    }}
    /* Compact, single-line QC metrics table */
    table.qc-table {{
      width: auto;
      font-size: 0.76rem;
    }}
    table.qc-table th {{
      padding: 5px 6px;
      font-size: 0.66rem;
      white-space: normal;
      line-height: 1.2;
      vertical-align: bottom;
      cursor: pointer;
      user-select: none;
      position: sticky;
      top: 0;
      z-index: 1;
    }}
    table.qc-table th.num {{
      max-width: 72px;
    }}
    table.qc-table th:hover {{
      color: var(--text-main);
    }}
    table.qc-table th .sort-ind {{
      color: var(--brand-blue);
      font-size: 0.7rem;
    }}
    table.qc-table td {{
      padding: 3px 6px;
      height: 26px;
      white-space: nowrap;
      line-height: 1.2;
    }}
    table.qc-table tr:hover td {{
      background-color: transparent;
    }}
    table.qc-table td.qc-cell-fail, table.qc-table tr:hover td.qc-cell-fail {{
      background: #dc2626;
      color: #ffffff;
      font-weight: 800;
    }}
    table.qc-table td.qc-cell-warn, table.qc-table tr:hover td.qc-cell-warn {{
      background: #f59e0b;
      color: #1f1300;
      font-weight: 800;
    }}
    table.qc-table .badge-pill.red {{
      background: #dc2626;
      color: #ffffff;
      border-color: #b91c1c;
    }}
    table.qc-table .badge-pill.amber {{
      background: #f59e0b;
      color: #1f1300;
      border-color: #d97706;
    }}
    table.qc-table td.qc-iso {{
      max-width: 150px;
      overflow: hidden;
      text-overflow: ellipsis;
    }}

    .qc-banner {{
      background: #fef2f2;
      border: 1px solid #fecaca;
      color: #991b1b;
      border-radius: 6px;
      padding: 10px 14px;
      font-size: 0.82rem;
      margin-bottom: 16px;
    }}

    /* Patient Outliers & Antigens */
    .outlier-bar-row {{
      display: grid;
      grid-template-columns: 90px 1fr 120px;
      align-items: center;
      gap: 10px;
      padding: 3px 6px;
      border-radius: 4px;
      font-size: 0.76rem;
      cursor: pointer;
    }}
    .outlier-bar-row:hover {{
      background: var(--bg-subtle);
    }}
    .outlier-bar-track {{
      height: 14px;
      background: var(--bg-subtle);
      border-radius: 3px;
      overflow: hidden;
    }}
    .outlier-bar-fill {{
      height: 100%;
      border-radius: 3px;
    }}
    .ag-ct {{ background: #7c3aed; }}
    .ag-tcell {{ background: #1068eb; }}
    .ag-ctlike {{ background: #f59e0b; }}
    .ag-none {{ background: #cbd5e1; }}
    .outlier-legend {{
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
      font-size: 0.72rem;
      color: var(--text-muted);
    }}
    .outlier-legend i {{
      display: inline-block;
      width: 10px;
      height: 10px;
      border-radius: 2px;
      margin-right: 4px;
      vertical-align: -1px;
    }}
    tr.outlier-row-highlight td {{
      background: var(--brand-blue-light);
    }}

    /* Concordance Matrix & Spectrum */
    .concordance-matrix-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.78rem;
    }}
    .concordance-matrix-table th {{
      background: var(--bg-subtle);
      padding: 8px 10px;
      font-weight: 700;
      color: var(--text-muted);
      border: 1px solid var(--border-color);
      text-align: center;
      font-size: 0.74rem;
    }}
    .concordance-matrix-table td {{
      border: 1px solid var(--border-color);
      padding: 10px;
      text-align: center;
      transition: background 0.15s ease;
    }}
    .matrix-row-hdr {{
      font-weight: 700;
      color: var(--text-main);
      text-align: left !important;
      background: #ffffff;
      white-space: nowrap;
    }}
    .matrix-cell-match {{
      background: #ecfdf5;
      color: #065f46;
      font-weight: 800;
    }}
    .matrix-cell-divergent {{
      background: #fffbeb;
      color: #92400e;
      font-weight: 700;
    }}
    .matrix-cell-zero {{
      background: #ffffff;
      color: #cbd5e1;
    }}
    .spectrum-bar-row {{
      display: flex;
      align-items: center;
      gap: 10px;
      padding: 5px 8px;
      border-radius: 4px;
      font-size: 0.78rem;
      cursor: pointer;
      transition: background 0.1s ease;
    }}
    .spectrum-bar-row:hover {{
      background: var(--bg-subtle);
    }}
    .spectrum-track {{
      flex: 1;
      height: 18px;
      background: var(--bg-subtle);
      border-radius: 3px;
      overflow: hidden;
      position: relative;
      display: flex;
    }}
    .spectrum-fill {{
      height: 100%;
      border-radius: 3px;
      transition: width 0.3s ease;
    }}
    .spectrum-fill.concordant {{
      background: #10b981;
    }}
    .spectrum-fill.discordant {{
      background: #f59e0b;
    }}
    .spectrum-fill.pending {{
      background: #cbd5e1;
    }}

    .drawer-scrim.open {{
      display: block;
    }}
  </style>
</head>
<body>

  <!-- Top Navigation Header -->
  <header class="app-header">
    <div class="header-left">
      <a href="#" class="app-brand">
        <span>ONCOANALYSER</span>
      </a>
      <span class="brand-tag">Cohort Clinical Explorer</span>
      <div class="header-divider"></div>
      <div>
        <div class="project-title">MMR Cohort: RNA-Seq, Antigens & TCR Repertoire</div>
        <div class="project-subtitle">Isofox TPM • CIDER V(D)J Clones • Gene Fusions • CUPPA Tissue Classification</div>
      </div>
    </div>
    <div class="header-right">
      <button class="btn-header" onclick="exportCohortLandscapeCSV()">
        <span>Export CSV</span>
      </button>
      <button class="btn-header btn-primary" onclick="window.print()">
        <span>Print Summary</span>
      </button>
    </div>
  </header>

  <!-- Application Body: Sidebar + Main Content -->
  <div class="app-container">

    <!-- Left Sidebar -->
    <aside class="app-sidebar">
      <div class="sidebar-section-title">Navigation Views</div>
      <ul class="nav-list">
        <li class="nav-item active" onclick="switchView('tab-landscape', this)">
          <span>Overview & OncoPrint</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-qc', this)">
          <span>RNA Sample QC</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-volcano', this)">
          <span>Differential Volcano</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-genes', this)">
          <span>Expression & Trajectory</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-outliers', this)">
          <span>Patient Outliers & Antigens</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-tcr', this)">
          <span>TCR & Immune Repertoire</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-fusions', this)">
          <span>RNA Fusions & Chimeras</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-cuppa', this)">
          <span>CUPPA Tissue Lineage</span>
        </li>
        <li class="nav-item" onclick="switchView('tab-samples', this)">
          <span>Sample Master Catalog</span>
        </li>
      </ul>

      <!-- Sidebar Gene Selector Panel -->
      <div class="sidebar-gene-selector hidden" id="sidebar-gene-selector">
        <div class="sidebar-section-title" style="padding-top: 8px;">Gene Selector</div>
        <div class="gene-filter-header">
          <div class="gene-selector-context" id="gene-selector-context"></div>
          <input type="text" id="sidebar-gene-search" class="gene-search-input" placeholder="Search any assayed gene..." oninput="handleGeneFilter(this.value)">
          <div class="gene-cat-chips">
            <span class="cat-chip active" onclick="filterGeneCategory('SUGGESTED', this)">Suggested</span>
            <span class="cat-chip" onclick="filterGeneCategory('ALL', this)">All</span>
            <span class="cat-chip" onclick="filterGeneCategory('CTA', this)">CTAs</span>
            <span class="cat-chip" onclick="filterGeneCategory('GBM', this)">GBM</span>
            <span class="cat-chip" onclick="filterGeneCategory('GI', this)">GI</span>
            <span class="cat-chip" onclick="filterGeneCategory('IMMUNE', this)">Immune</span>
            <span class="cat-chip" onclick="filterGeneCategory('HLA', this)">HLA</span>
          </div>
        </div>
        <div class="gene-scroll-list" id="gene-scroll-list">
          <!-- Populated by JS -->
        </div>
      </div>
    </aside>

    <!-- Main Viewport Area -->
    <main class="app-viewport">

      <!-- VIEW 1: COHORT OVERVIEW & ONCOPRINT -->
      <div id="tab-landscape" class="tab-pane active">
        <div class="view-header">
          <div class="view-title-group">
            <h2>Cohort Landscape & Biomarker Matrix</h2>
            <p>Clinical diagnosis, HLA Class I genotypes, immune infiltration, and somatic RNA fusions across MMR patient samples.</p>
          </div>
        </div>

        <!-- KPI Statistic Row -->
        <div class="kpi-grid">
          <div class="kpi-card accent-blue">
            <div class="kpi-title">Cohort Analyzed</div>
            <div class="kpi-num" id="kpi-total-samples">25</div>
            <div class="kpi-sub" id="kpi-patients-sub">-</div>
          </div>
          <div class="kpi-card accent-purple">
            <div class="kpi-title">Median Library Depth</div>
            <div class="kpi-num" id="kpi-median-depth">134M</div>
            <div class="kpi-sub">Total RNA-Seq Fragments</div>
          </div>
          <div class="kpi-card accent-amber">
            <div class="kpi-title">Productive TCR Clonotypes</div>
            <div class="kpi-num" id="kpi-total-cdr3" style="color: var(--brand-blue);">-</div>
            <div class="kpi-sub">CIDER PASS, in-frame CDR3s</div>
          </div>
          <div class="kpi-card accent-green">
            <div class="kpi-title">Candidate RNA Fusions</div>
            <div class="kpi-num" id="kpi-total-fusions" style="color: var(--success-green);">-</div>
            <div class="kpi-sub" id="kpi-fusions-sub">Known or distant partners</div>
          </div>
        </div>

        <!-- OncoPrint Biomarker Heatmap -->
        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Cohort OncoPrint Expression Matrix</div>
              <div class="card-caption" id="oncoprint-caption">Expression of key tumour antigens, immune checkpoints and HLA genes relative to a reference, not absolute TPM.</div>
            </div>
            <div style="display: flex; flex-direction: column; gap: 6px; align-items: flex-end;">
              <select id="oncoprint-mode" class="control-select" onchange="renderOncoPrint()" style="min-width: 260px;">
                <option value="cohort" selected>vs matched cohort (origin + library, z-score)</option>
                <option value="hartwig">vs Hartwig matched cancer type (percentile)</option>
                <option value="hartwig_all">vs Hartwig pan-cancer (percentile)</option>
              </select>
              <div id="oncoprint-legend" style="font-size: 0.74rem; display: flex; gap: 10px; color: var(--text-muted); flex-wrap: wrap; justify-content: flex-end;"></div>
            </div>
          </div>
          <div class="oncoprint-wrapper" id="oncoprint-wrapper">
            <!-- Populated by JS -->
          </div>
        </div>
      </div>

      <!-- VIEW 2: VOLCANO DIFFERENTIAL EXPRESSION -->
      <!-- VIEW: RNA SAMPLE QC -->
      <div id="tab-qc" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>RNA Sample Quality Control</h2>
            <p>Library depth, complexity, splicing, strand specificity and alignment metrics from Isofox and STAR. Failed samples are excluded from cohort-level expression analyses unless included below.</p>
          </div>
          <button class="btn-header btn-primary" onclick="exportQcCSV()">
            <span>Export QC CSV</span>
          </button>
        </div>

        <div class="kpi-grid">
          <div class="kpi-card accent-green">
            <div class="kpi-title">QC Pass</div>
            <div class="kpi-num" id="qc-kpi-pass" style="color: var(--success-green);">-</div>
            <div class="kpi-sub">All checks within thresholds</div>
          </div>
          <div class="kpi-card accent-amber">
            <div class="kpi-title">QC Warn</div>
            <div class="kpi-num" id="qc-kpi-warn" style="color: var(--warning-amber);">-</div>
            <div class="kpi-sub">Usable; review flagged metrics</div>
          </div>
          <div class="kpi-card accent-red">
            <div class="kpi-title">QC Fail</div>
            <div class="kpi-num" id="qc-kpi-fail" style="color: #ef4444;">-</div>
            <div class="kpi-sub">Excluded from expression statistics</div>
          </div>
          <div class="kpi-card accent-blue">
            <div class="kpi-title">Library Types</div>
            <div class="kpi-num" id="qc-kpi-libs" style="font-size: 1.4rem;">-</div>
            <div class="kpi-sub">Bulk / Target capture</div>
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Sample QC Metrics</div>
              <div class="card-caption">Red = fail threshold, amber = warn threshold. Hover a cell for the rule.</div>
            </div>
            <label class="control-label" style="display: flex; gap: 6px; align-items: center; cursor: pointer; margin: 0;">
              <input type="checkbox" class="include-qc-fail-toggle" onchange="setIncludeQcFail(this.checked)">
              Include QC-failed samples in expression analyses
            </label>
          </div>
          <div class="table-responsive" style="max-height: 620px; overflow-y: auto;">
            <table class="data-table qc-table">
              <thead>
                <tr>
                  <th data-sort="tissue_id" onclick="sortQcTable('tissue_id')">Sample <span class="sort-ind"></span></th>
                  <th data-sort="tissue_origin" onclick="sortQcTable('tissue_origin')">Origin <span class="sort-ind"></span></th>
                  <th data-sort="library_type" onclick="sortQcTable('library_type')">Library <span class="sort-ind"></span></th>
                  <th data-sort="tissue_format" onclick="sortQcTable('tissue_format')" title="Preservation">Format <span class="sort-ind"></span></th>
                  <th data-sort="status" onclick="sortQcTable('status')">QC <span class="sort-ind"></span></th>
                  <th class="num" data-sort="total_fragments" onclick="sortQcTable('total_fragments')">Total Frags (M) <span class="sort-ind"></span></th>
                  <th class="num" data-sort="unique_fragments" onclick="sortQcTable('unique_fragments')">Unique Frags (M) <span class="sort-ind"></span></th>
                  <th class="num" data-sort="duplicate_pct" onclick="sortQcTable('duplicate_pct')">Duplicate % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="genes_detected" onclick="sortQcTable('genes_detected')">Genes ≥ 1 TPM <span class="sort-ind"></span></th>
                  <th class="num" data-sort="spliced_pct" onclick="sortQcTable('spliced_pct')">Spliced % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="unspliced_pct" onclick="sortQcTable('unspliced_pct')">Unspliced % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="chimeric_pct" onclick="sortQcTable('chimeric_pct')">Chimeric % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="forward_strand_pct" onclick="sortQcTable('forward_strand_pct')">Fwd Strand % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="uniquely_mapped_pct" onclick="sortQcTable('uniquely_mapped_pct')">STAR Unique % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="top100_pct" onclick="sortQcTable('top100_pct')">Top-100 Share % <span class="sort-ind"></span></th>
                  <th class="num" data-sort="median_frag_length" onclick="sortQcTable('median_frag_length')">Median Frag (bp) <span class="sort-ind"></span></th>
                  <th class="num" data-sort="median_gc_ratio" onclick="sortQcTable('median_gc_ratio')">GC <span class="sort-ind"></span></th>
                  <th data-sort="isofox_status" onclick="sortQcTable('isofox_status')">Isofox Status <span class="sort-ind"></span></th>
                </tr>
              </thead>
              <tbody id="qc-table-body">
              </tbody>
            </table>
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">QC Rules</div>
              <div class="card-caption">Any fail rule, or an Isofox FAIL_* status, fails the sample. Strand check is skipped automatically for unstranded cohorts.</div>
            </div>
          </div>
          <div class="table-responsive">
            <table class="data-table">
              <thead>
                <tr><th>Metric</th><th>Fail</th><th>Warn</th><th>Rationale</th></tr>
              </thead>
              <tbody id="qc-rules-body">
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <div id="tab-volcano" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>Interactive Differential Expression (Volcano Plot)</h2>
            <p>Welch's t-test on log2(TPM+1) between any two metadata or Oncoanalyser categories, restricted to samples with Isofox expression data. Log2FC is the difference in mean log2(TPM+1).</p>
          </div>
          <button class="btn-header btn-primary" onclick="exportVolcanoCSV()">
            <span>Export Differential CSV</span>
          </button>
        </div>

        <!-- Volcano Comparison Controls -->
        <div class="volcano-controls">
          <div class="control-group">
            <label class="control-label">Comparison Category</label>
            <select id="volcano-category" class="control-select" onchange="onVolcanoCategoryChange()">
              <option value="tissue_origin">Tissue Origin (CNS vs GI / Other)</option>
              <option value="path_dx">Pathology Diagnosis</option>
              <option value="timepoint">Surgical Resection / Timepoint</option>
              <option value="tissue_format">Preservation Format (Frozen vs FFPE)</option>
              <option value="library_type">Library Type (Bulk vs Target Capture)</option>
              <option value="cuppa_pred">CUPPA Predicted Primary Lineage</option>
              <option value="tcr_group">TCR Reads per Million Fragments (above vs below median)</option>
              <option value="cytolytic_group">Cytolytic Score, GZMA + PRF1 (above vs below median)</option>
              <option value="hla_a">HLA-A Allele (carriers)</option>
              <option value="hla_b">HLA-B Allele (carriers)</option>
              <option value="hla_c">HLA-C Allele (carriers)</option>
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Library</label>
            <select id="volcano-library" class="control-select" onchange="onVolcanoCategoryChange()">
              <option value="ALL" selected>All libraries</option>
              <option value="Bulk">Bulk RNA-seq only</option>
              <option value="Target Capture">Target capture only</option>
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Preservation</label>
            <select id="volcano-format" class="control-select" onchange="onVolcanoCategoryChange()">
              <option value="ALL" selected>All formats</option>
              <option value="FROZEN">Fresh Frozen only</option>
              <option value="FFPE">FFPE only</option>
            </select>
          </div>
          <div class="control-group" style="justify-content: flex-end;">
            <label class="control-label" style="display: flex; gap: 6px; align-items: center; cursor: pointer;">
              <input type="checkbox" class="include-qc-fail-toggle" onchange="setIncludeQcFail(this.checked)">
              Include QC-failed
            </label>
          </div>
          <div class="control-group">
            <label class="control-label">Group 1 (Numerator / Test)</label>
            <select id="volcano-group1" class="control-select" onchange="updateVolcanoPlot()">
              <!-- Dynamic options -->
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Group 2 (Baseline / Reference)</label>
            <select id="volcano-group2" class="control-select" onchange="updateVolcanoPlot()">
              <!-- Dynamic options -->
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Min |Log2FC|</label>
            <select id="volcano-fc-cutoff" class="control-select" onchange="updateVolcanoPlot()">
              <option value="0.585">1.5x (Log2FC ≥ 0.58)</option>
              <option value="1.0" selected>2.0x (Log2FC ≥ 1.0)</option>
              <option value="1.585">3.0x (Log2FC ≥ 1.58)</option>
              <option value="2.0">4.0x (Log2FC ≥ 2.0)</option>
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Significance</label>
            <select id="volcano-p-type" class="control-select" onchange="updateVolcanoPlot()">
              <option value="raw">Raw p-value</option>
              <option value="fdr" selected>BH FDR (q-value)</option>
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Max P / Q</label>
            <select id="volcano-p-cutoff" class="control-select" onchange="updateVolcanoPlot()">
              <option value="0.05" selected>p ≤ 0.05 (-log10 ≥ 1.30)</option>
              <option value="0.01">p ≤ 0.01 (-log10 ≥ 2.00)</option>
              <option value="0.001">p ≤ 0.001 (-log10 ≥ 3.00)</option>
            </select>
          </div>
          <div class="control-group" style="flex: 1; min-width: 200px;">
            <label class="control-label">Search / Highlight Gene</label>
            <input type="text" id="volcano-search" class="gene-search-input" placeholder="Type symbol (e.g. EGFR, GFAP, SOX2, PRAME)..." oninput="onVolcanoSearch(this.value)">
          </div>
        </div>

        <div class="volcano-layout">
          <!-- Canvas Plot Column -->
          <div class="volcano-canvas-container">
            <div style="width: 100%; display: flex; justify-content: space-between; font-size: 0.78rem; color: var(--text-muted); margin-bottom: 6px;">
              <span id="volcano-counts-badge">Plotting differential genes...</span>
              <span style="display: flex; gap: 12px;">
                <span><span style="color: var(--up-red);">●</span> Enriched in Grp 1</span>
                <span><span style="color: var(--down-blue);">●</span> Enriched in Grp 2</span>
                <span><span style="color: #94a3b8;">●</span> Not Significant</span>
              </span>
            </div>
            <canvas id="volcano-canvas"></canvas>
            <div class="volcano-tooltip" id="volcano-tooltip"></div>
          </div>

          <!-- Side Selected Gene & Top Hits Column -->
          <div style="display: flex; flex-direction: column; gap: 14px;">
            <!-- Selected Gene Card -->
            <div class="card" id="volcano-gene-card">
              <div class="card-top" style="margin-bottom: 8px; padding-bottom: 6px;">
                <div class="card-heading" id="volcano-gene-name" style="font-family: 'JetBrains Mono', monospace; color: var(--brand-blue);">Select / Hover a Gene</div>
                <span id="volcano-gene-badge" class="badge-pill blue">-</span>
              </div>
              <div style="font-size: 0.8rem; display: flex; flex-direction: column; gap: 6px;" id="volcano-gene-details">
                <div style="display: flex; justify-content: space-between;"><span style="color: var(--text-muted);">Log2 Fold Change:</span><b id="v-det-fc">-</b></div>
                <div style="display: flex; justify-content: space-between;"><span style="color: var(--text-muted);">Welch p-value:</span><b id="v-det-pval">-</b></div>
                <div style="display: flex; justify-content: space-between;"><span style="color: var(--text-muted);">Group 1 Mean:</span><b id="v-det-m1">- TPM</b></div>
                <div style="display: flex; justify-content: space-between;"><span style="color: var(--text-muted);">Group 2 Mean:</span><b id="v-det-m2">- TPM</b></div>
                <button class="btn-header" style="margin-top: 6px; justify-content: center;" onclick="inspectSelectedGeneTrajectory()">
                  Inspect Trajectory View
                </button>
              </div>
            </div>

            <!-- Top Differentials Table -->
            <div class="card" style="flex: 1; display: flex; flex-direction: column;">
              <div class="card-heading" style="margin-bottom: 8px;">Top Differential Genes</div>
              <div class="table-responsive" style="max-height: 270px; overflow-y: auto;">
                <table class="data-table" style="font-size: 0.76rem;">
                  <thead>
                    <tr><th>Gene</th><th class="num">Log2FC</th><th class="num">p-val</th></tr>
                  </thead>
                  <tbody id="volcano-top-tbody">
                    <!-- Populated by JS -->
                  </tbody>
                </table>
              </div>
            </div>
          </div>
        </div>
      </div>

      <!-- VIEW 3: GENE EXPRESSION & TRAJECTORIES -->
      <div id="tab-genes" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>Gene Expression & Longitudinal Trajectories</h2>
            <p>Expression of up to three genes across one patient's samples in collection order, with a separate line per tumour type, and the first gene across tumour origins.</p>
          </div>
        </div>

        <div class="volcano-controls">
          <div class="control-group" style="min-width: 300px;">
            <label class="control-label">Patient (≥ 2 usable samples)</label>
            <select id="traj-patient" class="control-select" onchange="onTrajPatientChange()"></select>
          </div>
          <div class="control-group" style="flex: 1; min-width: 260px;">
            <label class="control-label">Genes (up to 3)</label>
            <div id="traj-gene-chips" style="display: flex; gap: 6px; flex-wrap: wrap; align-items: center; min-height: 32px;"></div>
          </div>
          <div class="control-group" style="min-width: 220px;">
            <label class="control-label">Add Gene</label>
            <div style="display: flex; gap: 6px;">
              <input type="text" id="traj-gene-input" class="gene-search-input" list="traj-gene-options" placeholder="e.g. PRAME" onkeydown="if (event.key === 'Enter') addTrajGene()">
              <button class="btn-header" onclick="addTrajGene()">Add</button>
            </div>
            <datalist id="traj-gene-options"></datalist>
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading" id="traj-heading">Longitudinal Trajectory</div>
              <div class="card-caption" id="traj-caption"></div>
            </div>
            <div id="traj-legend" style="font-size: 0.74rem; display: flex; gap: 12px; color: var(--text-muted); flex-wrap: wrap; justify-content: flex-end;"></div>
          </div>
          <div id="gene-longitudinal-svg" style="height: 320px;"></div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Expression Across Tumor Origins: <span class="mono" id="active-gene-symbol" style="color: var(--brand-blue);"></span></div>
              <div class="card-caption">Median TPM per tissue origin (first trajectory gene), usable samples only.</div>
            </div>
          </div>
          <div id="gene-grouped-box-svg" style="height: 230px; display: flex; align-items: center; justify-content: center;"></div>
        </div>
      </div>

      <!-- VIEW: PATIENT OUTLIER EXPRESSION & ANTIGENS -->
      <div id="tab-outliers" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>Patient Outlier Expression & Tumour Antigens</h2>
            <p>Genes over-expressed in one sample relative to matched cohort samples (same tissue origin and library type where possible; leave-one-out z-score) and the Hartwig reference for the sample's cancer type (Isofox gene_distribution percentiles), labelled against IEDB tumour T-cell antigens and testis-restricted (cancer-testis) genes.</p>
          </div>
          <button class="btn-header btn-primary" onclick="exportOutliersCSV()">
            <span>Export Sample CSV</span>
          </button>
        </div>

        <div class="volcano-controls">
          <div class="control-group" style="min-width: 260px;">
            <label class="control-label">Sample</label>
            <select id="outlier-sample" class="control-select" onchange="renderOutlierView()">
              <!-- Populated by JS -->
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Gene Set</label>
            <select id="outlier-mode" class="control-select" onchange="renderOutlierView()">
              <option value="outliers" selected>Outlier genes</option>
              <option value="antigens">All expressed known antigens</option>
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Antigen Label</label>
            <select id="outlier-antigen-filter" class="control-select" onchange="renderOutlierView()">
              <option value="ALL" selected>All genes</option>
              <option value="ANY">Any known / CT-like antigen</option>
              <option value="TCELL">IEDB T-cell antigen only</option>
              <option value="CT">Cancer-testis (any)</option>
              <option value="NONE">Not a known antigen</option>
            </select>
          </div>
          <div class="control-group">
            <label class="control-label">Rank By</label>
            <select id="outlier-sort" class="control-select" onchange="renderOutlierView()">
              <option value="ref_fc" selected>Log2FC vs Hartwig reference median</option>
              <option value="z">Cohort z-score</option>
              <option value="tpm">TPM</option>
            </select>
          </div>
          <div class="control-group" style="justify-content: flex-end;">
            <label class="control-label" style="display: flex; gap: 6px; align-items: center; cursor: pointer;">
              <input type="checkbox" id="outlier-hide-noise" checked onchange="renderOutlierView()">
              Protein-coding only (hide Ig / TCR / MT / ribosomal)
            </label>
          </div>
          <div class="control-group" style="flex: 1; min-width: 160px;">
            <label class="control-label">Search Gene</label>
            <input type="text" id="outlier-search" class="gene-search-input" placeholder="e.g. PRAME, MAGEA4..." oninput="renderOutlierView()">
          </div>
        </div>

        <div id="outlier-qc-banner" class="qc-banner" style="display: none;"></div>
        <div id="outlier-ref-note" class="card-caption" style="margin: 0 0 12px;"></div>

        <div class="kpi-grid">
          <div class="kpi-card accent-blue">
            <div class="kpi-title">Outlier Genes</div>
            <div class="kpi-num" id="outlier-kpi-count">-</div>
            <div class="kpi-sub" id="outlier-kpi-count-sub">Passing outlier criteria</div>
          </div>
          <div class="kpi-card accent-purple">
            <div class="kpi-title">Known T-cell Antigens</div>
            <div class="kpi-num" id="outlier-kpi-tcell">-</div>
            <div class="kpi-sub">Outliers with IEDB cancer T-cell epitopes</div>
          </div>
          <div class="kpi-card accent-amber">
            <div class="kpi-title">Cancer-Testis Genes</div>
            <div class="kpi-num" id="outlier-kpi-ct">-</div>
            <div class="kpi-sub">Outliers with testis-restricted expression</div>
          </div>
          <div class="kpi-card accent-green">
            <div class="kpi-title">Top Antigen Outlier</div>
            <div class="kpi-num mono" id="outlier-kpi-top" style="font-size: 1.3rem;">-</div>
            <div class="kpi-sub" id="outlier-kpi-top-sub">-</div>
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Top 30 Genes</div>
              <div class="card-caption" id="outlier-chart-caption">Bar length = ranking metric. Click a bar to open the gene in the expression view.</div>
            </div>
            <div class="outlier-legend">
              <span><i class="ag-ct"></i>Cancer-testis T-cell antigen</span>
              <span><i class="ag-tcell"></i>T-cell antigen (IEDB)</span>
              <span><i class="ag-ctlike"></i>Cancer-testis-like (no IEDB epitope)</span>
              <span><i class="ag-none"></i>Not annotated</span>
            </div>
          </div>
          <div id="outlier-bars" style="padding: 4px 0;">
            <!-- Generated by JS -->
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Outlier Gene Table</div>
              <div class="card-caption">Outlier = TPM ≥ 2 and (cohort z ≥ 3, or cohort z ≥ 1.5 with ≥ 99th percentile and log2FC ≥ 2 vs the median of the sample's Hartwig cancer type, from the Isofox gene_distribution reference; pan-cancer when no type applies). Cohort z reference: same tissue origin + library type when ≥ 5 QC-passing samples, else same library type, else same origin, else whole cohort. Antigens: IEDB positive human T-cell assays in cancer (DOID:162); cancer-testis = HPA testis tissue/group enriched or canonical CT gene family (resources/antigens/).</div>
            </div>
            <div class="card-caption" id="outlier-table-count"></div>
          </div>
          <div class="table-responsive" style="max-height: 520px; overflow-y: auto;">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Gene</th>
                  <th class="num">TPM</th>
                  <th class="num">Cohort z</th>
                  <th class="num">Pan-Cancer Pctl</th>
                  <th class="num">Cancer-Type Pctl</th>
                  <th class="num">Hartwig Ref Median TPM</th>
                  <th class="num">Log2FC vs Hartwig Ref</th>
                  <th>Outlier</th>
                  <th>Antigen Label</th>
                  <th class="num">IEDB Assays / Refs</th>
                  <th>MHC</th>
                  <th>IEDB Tumour Contexts</th>
                </tr>
              </thead>
              <tbody id="outlier-table-body">
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <!-- VIEW 4: TCR IMMUNE REPERTOIRE -->
      <div id="tab-tcr" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>TCR Repertoire Reconstruction & Clonality</h2>
            <p>CIDER CDR3 amino acid sequences, V/J germline segment usage, and Simpson clonality metrics across tumor specimens.</p>
          </div>
          <button class="btn-header btn-primary" onclick="exportClonotypesCSV()">Export Repertoire CSV</button>
        </div>

        <div class="card">
          <div class="card-top">
            <div class="card-heading">Immune Locus Sequencing Depth</div>
            <div style="display: flex; gap: 12px; font-size: 0.78rem;">
              <span><span style="color: var(--brand-blue);">■</span> TCR Reads</span>
              <span><span style="color: var(--warning-amber);">■</span> BCR Reads</span>
            </div>
          </div>
          <div id="tcr-locus-chart" style="height: 160px; margin-bottom: 8px;"></div>
        </div>

        <div class="card">
          <div class="card-top">
            <div class="card-heading">Top Reconstructed CDR3 Clonotypes</div>
            <div style="display: flex; gap: 8px;">
              <input type="text" id="cdr3-filter" class="gene-search-input" style="width: 250px;" placeholder="Filter CDR3, V-gene, or sample..." oninput="filterClonotypes(this.value)">
              <select id="cdr3-locus-select" class="btn-header" onchange="filterClonotypes()">
                <option value="TCR" selected>All TCR</option>
                <option value="TRB">TRB</option>
                <option value="TRA/TRD">TRA / TRD</option>
                <option value="TRG">TRG</option>
                <option value="BCR">All BCR</option>
                <option value="IGH">IGH</option>
                <option value="IGK">IGK</option>
                <option value="IGL">IGL</option>
                <option value="ALL">All loci</option>
              </select>
              <label class="control-label" style="display: flex; gap: 6px; align-items: center; margin: 0;">
                <input type="checkbox" id="cdr3-productive" checked onchange="filterClonotypes()"> Productive only
              </label>
            </div>
          </div>
          <div class="table-responsive">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Sample ID</th>
                  <th>Patient</th>
                  <th>Locus</th>
                  <th>V-Gene</th>
                  <th>J-Gene</th>
                  <th>CDR3 Sequence</th>
                  <th class="num">Length</th>
                  <th class="num">Read Count</th>
                  <th>Frame Status</th>
                  <th>Also In (same patient)</th>
                </tr>
              </thead>
              <tbody id="cdr3-table-body">
                <!-- Populated by JS -->
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <!-- VIEW 5: RNA FUSIONS -->
      <div id="tab-fusions" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>RNA Fusions, Chimeras & T-cell Context</h2>
            <p>Isofox chimeric transcripts classified by partner type, immune-receptor chimeras cross-checked against CIDER, and candidate fusions alongside each sample's T-cell repertoire. Isofox pass_fusions carries no reading-frame call (that comes from LINX on DNA), so frame is not shown.</p>
          </div>
          <button class="btn-header btn-primary" onclick="exportFusionsCSV()">Export Fusions CSV</button>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Candidate Fusions vs T-cell Repertoire</div>
              <div class="card-caption">Per sample: candidate fusions (known/promiscuous or distant partners, ≥ 3 split reads, both breakpoints on exon boundaries) next to TCR depth and clonal structure. Co-occurrence only: linking a TCR to a fusion epitope needs HLA binding prediction of the junction peptide and functional validation. Clonality needs ≥ 10 productive TCR clonotypes.</div>
            </div>
          </div>
          <div class="table-responsive" style="max-height: 460px; overflow-y: auto;">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Sample</th>
                  <th>Diagnosis</th>
                  <th>QC</th>
                  <th class="num">Candidate Fusions</th>
                  <th>Top Candidates (split reads)</th>
                  <th class="num">TCR Reads / M</th>
                  <th class="num">Productive TCR Clones</th>
                  <th class="num">Clonality</th>
                  <th class="num">Top Clone Share</th>
                  <th class="num">Clones Shared in Patient</th>
                  <th class="num">Fusion Neoantigens (HLA Binders)</th>
                </tr>
              </thead>
              <tbody id="fusion-tcr-body"></tbody>
            </table>
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Fusion Neoantigens × TCR Repertoire Binding</div>
              <div class="card-caption">Candidate somatic fusions translated across breakpoint junctions into candidate neoepitopes, scored for HLA class I presentation (Hartwig Neo binding rank & likelihood) against patient alleles, and paired with detected TCR clonotypes from CIDER and MiXCR cap-TCR experiments.</div>
            </div>
            <div class="card-caption" id="receptor-match-count"></div>
          </div>
          <div class="table-responsive" style="max-height: 400px; overflow-y: auto;">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Sample</th>
                  <th>Fusion</th>
                  <th>Classification</th>
                  <th>Neoepitope</th>
                  <th>Best HLA Allele</th>
                  <th class="num">HLA Rank (%)</th>
                  <th>Top Matching TCR (CDR3, V/J)</th>
                  <th>TCR Source</th>
                  <th class="num">TCR Reads</th>
                  <th class="num">Binding Score</th>
                </tr>
              </thead>
              <tbody id="receptor-match-body"></tbody>
            </table>
          </div>
        </div>

        <div class="card">
          <div class="card-top">
            <div>
              <div class="card-heading">Cohort Fusion Catalog</div>
              <div class="card-caption" id="fusion-catalog-caption">Sorted by supporting fragments.</div>
            </div>
            <div style="display: flex; gap: 8px; align-items: center; flex-wrap: wrap;">
              <select id="fusion-class-filter" class="control-select" onchange="renderFusionsTable()">
                <option value="CANDIDATE" selected>Candidates (known + distant)</option>
                <option value="Known / promiscuous">Known / promiscuous partner</option>
                <option value="Distant / interchromosomal">Distant / interchromosomal</option>
                <option value="Local (< 1 Mb)">Local (&lt; 1 Mb, mostly read-through)</option>
                <option value="Paralog / pseudogene">Paralog / pseudogene (artefact-prone)</option>
                <option value="Immune receptor">Immune receptor</option>
                <option value="ALL">All classes</option>
              </select>
              <select id="fusion-min-split" class="control-select" onchange="renderFusionsTable()">
                <option value="0">Any split reads</option>
                <option value="3" selected>≥ 3 split reads</option>
                <option value="5">≥ 5 split reads</option>
                <option value="10">≥ 10 split reads</option>
              </select>
              <label class="control-label" style="display: flex; gap: 6px; align-items: center; margin: 0;">
                <input type="checkbox" id="fusion-exon-only" checked onchange="renderFusionsTable()"> Exon-boundary junctions only
              </label>
              <input type="text" id="fusion-search" class="gene-search-input" style="width: 180px;" placeholder="Gene or sample..." oninput="renderFusionsTable()">
            </div>
          </div>
          <div class="table-responsive" style="max-height: 520px; overflow-y: auto;">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Sample ID</th>
                  <th>Diagnosis</th>
                  <th>5' Upstream Gene</th>
                  <th>3' Downstream Gene</th>
                  <th>Genomic Breakpoints</th>
                  <th class="num">Split Reads</th>
                  <th class="num">Spanning Reads</th>
                  <th>Class</th>
                  <th>Isofox Known Type</th>
                  <th>Junction</th>
                  <th class="num" title="Distinct patients with this gene pair; recurrence of an unknown pair across unrelated patients suggests an artefact">Patients</th>
                  <th>Action</th>
                </tr>
              </thead>
              <tbody id="fusions-table-body">
                <!-- Populated by JS -->
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <!-- VIEW 6: CUPPA TISSUE LINEAGE -->
      <div id="tab-cuppa" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>CUPPA Molecular Tissue of Origin & Concordance</h2>
            <p>Cross-classification and concordance analytics comparing clinical pathology diagnosis against transcriptomic machine-learning predictions.</p>
          </div>
          <button class="btn-header btn-primary" onclick="exportCuppaCSV()">
            <span>Export Lineage CSV</span>
          </button>
        </div>

        <!-- Concordance KPI Row -->
        <div class="kpi-grid">
          <div class="kpi-card accent-green">
            <div class="kpi-title">Overall Concordance Rate</div>
            <div class="kpi-num" id="cuppa-kpi-concordance" style="color: var(--success-green);">-</div>
            <div class="kpi-sub" id="cuppa-kpi-concordance-sub">Matching Sequenced Tumors</div>
          </div>
          <div class="kpi-card accent-blue">
            <div class="kpi-title">CNS Lineage Concordance</div>
            <div class="kpi-num" id="cuppa-kpi-cns">-</div>
            <div class="kpi-sub" id="cuppa-kpi-cns-sub">Glioblastoma / Astrocytoma</div>
          </div>
          <div class="kpi-card accent-purple">
            <div class="kpi-title">GI Lineage Concordance</div>
            <div class="kpi-num" id="cuppa-kpi-gi">-</div>
            <div class="kpi-sub" id="cuppa-kpi-gi-sub">Colorectal Adenocarcinomas</div>
          </div>
          <div class="kpi-card accent-amber">
            <div class="kpi-title">Median ML Likelihood</div>
            <div class="kpi-num" id="cuppa-kpi-median-lik">-</div>
            <div class="kpi-sub">Confidence Score</div>
          </div>
        </div>

        <!-- Visual Concordance Charts Grid -->
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 16px;">
          <!-- Chart 1: Cross-Classification Matrix -->
          <div class="card">
            <div class="card-top">
              <div>
                <div class="card-heading">Diagnosis vs CUPPA Concordance Matrix</div>
                <div class="card-caption">Heatmap comparing clinical pathology organ origin with transcriptomic classification.</div>
              </div>
              <div style="font-size: 0.72rem; display: flex; gap: 8px;">
                <span><span style="color: var(--success-green);">■</span> Concordant</span>
                <span><span style="color: var(--warning-amber);">■</span> Divergent</span>
              </div>
            </div>
            <div id="cuppa-matrix-container" style="overflow-x: auto; padding: 4px 0;">
              <!-- Generated by JS -->
            </div>
          </div>

          <!-- Chart 2: Per-Sample Likelihood & Concordance Spectrum -->
          <div class="card">
            <div class="card-top">
              <div>
                <div class="card-heading">Prediction Likelihood & Concordance Spectrum</div>
                <div class="card-caption">Ranked confidence scores per sample. Click any bar to open sample profile.</div>
              </div>
            </div>
            <div id="cuppa-spectrum-container" style="max-height: 290px; overflow-y: auto; padding-right: 4px;">
              <!-- Generated by JS -->
            </div>
          </div>
        </div>

        <!-- Detailed Classification Table -->
        <div class="card">
          <div class="card-top">
            <div class="card-heading">Sample Classification Catalog</div>
            <div class="card-caption">Detailed CUPPA primary cancer type calls, likelihoods, and secondary lineage candidates.</div>
          </div>
          <div class="table-responsive">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Sample ID</th>
                  <th>Annotated Diagnosis</th>
                  <th>Tissue Origin</th>
                  <th>CUPPA Predicted Primary</th>
                  <th class="num">Likelihood</th>
                  <th>Secondary Lineage</th>
                  <th>Concordance</th>
                  <th>Action</th>
                </tr>
              </thead>
              <tbody id="cuppa-table-body">
                <!-- Populated by JS -->
              </tbody>
            </table>
          </div>
        </div>
      </div>

      <!-- VIEW 7: SAMPLE CATALOG -->
      <div id="tab-samples" class="tab-pane">
        <div class="view-header">
          <div class="view-title-group">
            <h2>Sample Master Catalog & HLA Genotypes</h2>
            <p>Comprehensive metadata, 4-digit HLA Class I alleles, and quality summary.</p>
          </div>
        </div>

        <div class="card">
          <div class="table-responsive">
            <table class="data-table">
              <thead>
                <tr>
                  <th>Tissue ID</th>
                  <th>Patient ID</th>
                  <th>Diagnosis</th>
                  <th>Origin</th>
                  <th>Library</th>
                  <th>Preservation</th>
                  <th>RNA QC</th>
                  <th>HLA-A</th>
                  <th>HLA-B</th>
                  <th>HLA-C</th>
                  <th class="num">TCR Reads</th>
                  <th class="num">Clonality</th>
                  <th>Action</th>
                </tr>
              </thead>
              <tbody id="samples-table-body">
                <!-- Populated by JS -->
              </tbody>
            </table>
          </div>
        </div>
      </div>

    </main>
  </div>

  <!-- Slide-Over Inspector Drawer -->
  <div class="drawer-scrim" id="drawer-scrim" onclick="closeDrawer()"></div>
  <aside class="inspector-drawer" id="inspector-drawer">
    <div style="display: flex; justify-content: space-between; align-items: flex-start; margin-bottom: 18px; padding-bottom: 12px; border-bottom: 1px solid var(--border-color);">
      <div>
        <h3 id="drawer-sample-id" style="font-size: 1.25rem; font-weight: 800; color: var(--text-main);">Sample Details</h3>
        <p id="drawer-patient-id" style="color: var(--text-muted); font-size: 0.82rem;">Patient ID</p>
      </div>
      <button style="background: none; border: none; font-size: 1.4rem; color: var(--text-light); cursor: pointer;" onclick="closeDrawer()">&times;</button>
    </div>
    <div id="drawer-body">
      <!-- Injected by JS -->
    </div>
  </aside>

  <!-- Embedded Data & Interactive Engine -->
  <script>
    const COHORT = {cohort_json};
    const GENES = {gene_matrix_json};
    const FUSIONS = {fusions_json};
    const CLONOTYPES = {clonotypes_json};
    // Isofox immune-receptor chimeras with CIDER clonotypes that use the same V/J segments
    const RECEPTOR_MATCHES = {receptor_matches_json};
    // OUTLIERS[sample] rows: [gene, tpm, cohort_z, pan_pct, ref_median, ref_log2fc, cancer_pct, is_outlier, is_noise]
    // ref_* = Hartwig cancer-type reference (sample.hartwig_cancer_type), pan-cancer when unavailable
    const OUTLIERS = {outliers_json};
    // ANTIGENS[gene]: [class, iedb_assays, iedb_refs, mhc_classes, top_diseases]
    const ANTIGENS = {antigens_json};
    const VOLCANO_GENES = {volcano_genes_json};
    const QC_CHECKS = {qc_checks_json};
    // GENE_PAN_PCT[gene][sampleIdx]: Hartwig pan-cancer percentile (0-1) from Isofox, null if unavailable
    const GENE_PAN_PCT = {gene_pan_pct_json};
    // GENE_CAN_PCT[gene][sampleIdx]: percentile (0-1) within the sample's Hartwig cancer type, null if unavailable
    const GENE_CAN_PCT = {gene_can_pct_json};
    const hartwigPct = (mode, gene, idx) => ((mode === 'hartwig' ? GENE_CAN_PCT : GENE_PAN_PCT)[gene] || [])[idx];
    let includeQcFail = false;

    // Samples usable for cohort expression analyses: Isofox data present and QC not failed (unless overridden)
    function usableExpression(s) {{
      return s.has_expression && (includeQcFail || s.qc.status !== 'FAIL');
    }}
    let currentGene = 'EGFR';
    let currentCategory = 'SUGGESTED';
    // Tabs the sidebar gene selector drives, and what selecting a gene does there
    const GENE_SELECTOR_TABS = {{
      'tab-volcano': 'Highlights the gene on the volcano. Suggested = significant genes for the current comparison.',
      'tab-genes': 'Sets the first trajectory gene. Suggested = genes rising (▲) or falling (▼) across the selected patient resections.'
    }};
    let activeTab = 'tab-landscape';
    let trajGenes = [];
    let trajPatient = null;
    let trajGenesUserSet = false;
    const TRAJ_COLORS = ['#2563eb', '#dc2626', '#059669'];
    const TRAJ_MAX_GENES = 3;
    // Loci whose changes reflect immune infiltrate or housekeeping noise rather than tumour expression
    const NOISE_GENE_RE = /^(IG[HKL][VDJC]|TR[ABDG][VDJC]|MT-|RP[LS]\d|HB[ABDEGMQZ]\d?$|MIR\d|SNOR|LINC\d)/;

    const GENE_PANELS = {{
      'CTA': ['CTAG1B', 'CTAG2', 'MAGEA1', 'MAGEA3', 'MAGEA4', 'PRAME', 'SSX2'],
      'GBM': ['EGFR', 'SOX2', 'OLIG2', 'GFAP', 'IL13RA2', 'PROM1', 'ERBB2', 'NES', 'CD44'],
      'GI': ['CEACAM5', 'EPCAM', 'MUC1', 'MSLN', 'CDX2', 'KRT20'],
      'IMMUNE': ['PDCD1', 'CD274', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT', 'CD8A', 'CD8B', 'GZMA', 'PRF1'],
      'HLA': ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRA', 'HLA-DRB1', 'B2M', 'TAP1', 'TAP2']
    }};

    // High/low splits use the median of the samples currently eligible for the volcano (filters respected)
    const tcrPerMillion = s => {{
      const frags = (s.qc && s.qc.metrics && s.qc.metrics.total_fragments) || 0;
      return frags > 0 ? (s.tcr_reads || 0) / (frags / 1e6) : null;
    }};
    function eligibleMedian(valueFn) {{
      const vals = COHORT.filter(volcanoEligible).map(valueFn).filter(v => v !== null && v !== undefined).sort((a, b) => a - b);
      if (!vals.length) return null;
      const m = Math.floor(vals.length / 2);
      return vals.length % 2 ? vals[m] : (vals[m - 1] + vals[m]) / 2;
    }}

    // Volcano state
    let volcanoData = [];
    let highlightedVolcanoGene = null;
    let volcanoPThresholdNegLog = -Math.log10(0.05);

    document.addEventListener('DOMContentLoaded', () => {{
      updateKPIs();
      renderGeneSidebar();
      renderOncoPrint();
      renderGenePlots();
      renderTcrCharts();
      filterClonotypes();
      renderFusionsTable();
      renderFusionTcrView();
      renderCuppaView();
      renderSamplesTable();
      initVolcanoControls();
      initOutlierControls();
      renderQcView();
    }});

    function switchView(viewId, el) {{
      document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
      document.querySelectorAll('.nav-item').forEach(i => i.classList.remove('active'));
      document.getElementById(viewId).classList.add('active');
      if (!el) el = [...document.querySelectorAll('.nav-item')].find(n => n.getAttribute('onclick').includes(`'${{viewId}}'`));
      if (el) el.classList.add('active');
      activeTab = viewId;
      // Gene selector only where it drives the main panel
      document.getElementById('sidebar-gene-selector').classList.toggle('hidden', !GENE_SELECTOR_TABS[viewId]);
      document.getElementById('gene-selector-context').innerText = GENE_SELECTOR_TABS[viewId] || '';
      if (viewId === 'tab-volcano') {{
        setTimeout(() => {{ updateVolcanoPlot(); renderGeneSidebar(document.getElementById('sidebar-gene-search').value); }}, 50);
      }} else if (viewId === 'tab-genes') {{
        renderGenePlots();
        renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
      }}
    }}

    function updateKPIs() {{
      document.getElementById('kpi-total-samples').innerText = COHORT.length;
      document.getElementById('kpi-patients-sub').innerText = `${{new Set(COHORT.map(s => s.patient_id)).size}} distinct patients`;
      document.getElementById('kpi-total-cdr3').innerText = COHORT.reduce((a, s) => a + (s.tcr_clonotypes || 0), 0).toLocaleString();
      const cand = FUSIONS.filter(isCandidateFusion);
      document.getElementById('kpi-total-fusions').innerText = cand.length.toLocaleString();
      document.getElementById('kpi-fusions-sub').innerText = `of ${{FUSIONS.length.toLocaleString()}} Isofox chimeras (≥ ${{CANDIDATE_MIN_SPLIT}} split reads, known or distant partners)`;

      const depths = COHORT.map(s => s.isofox_summary?.total_fragments || 0).filter(d => d > 0);
      if (depths.length) {{
        depths.sort((a,b) => a-b);
        const med = depths[Math.floor(depths.length / 2)];
        document.getElementById('kpi-median-depth').innerText = `${{(med / 1e6).toFixed(0)}}M`;
      }}
    }}

    // Gene Sidebar List
    function renderGeneSidebar(filterTxt = '') {{
      const listEl = document.getElementById('gene-scroll-list');
      if (!listEl) return;

      const q = filterTxt.trim().toUpperCase();
      let items;
      if (currentCategory === 'SUGGESTED') {{
        items = suggestedGenes();
      }} else {{
        let genes = Object.keys(GENES);
        if (currentCategory !== 'ALL' && GENE_PANELS[currentCategory]) {{
          const allowed = new Set(GENE_PANELS[currentCategory]);
          genes = genes.filter(g => allowed.has(g));
        }}
        items = genes.sort().map(g => ({{ gene: g, meta: `max ${{Math.max(0, ...(GENES[g] || [])).toFixed(1)}}` }}));
      }}
      if (q) items = items.filter(it => it.gene.includes(q));

      const selected = new Set(activeTab === 'tab-volcano' ? [highlightedVolcanoGene] : trajGenes);
      listEl.innerHTML = items.slice(0, 150).map(it => `
          <div class="gene-list-item ${{selected.has(it.gene) ? 'selected' : ''}}" onclick="selectGene('${{it.gene}}')" title="${{it.title || ''}}">
            <span class="gene-symbol-txt">${{it.gene}}</span>
            <span class="gene-meta-txt">${{it.meta}}</span>
          </div>
        `).join('') || `<div class="gene-selector-context" style="padding: 10px 14px;">${{currentCategory === 'SUGGESTED' ? 'No suggestions for the current view.' : 'No genes match.'}}</div>`;
    }}

    function filterGeneCategory(cat, el) {{
      currentCategory = cat;
      document.querySelectorAll('.cat-chip').forEach(c => c.classList.remove('active'));
      if (el) el.classList.add('active');
      renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
    }}

    function handleGeneFilter(val) {{
      renderGeneSidebar(val);
    }}

    // Context-aware: on the volcano a gene is highlighted; elsewhere it becomes the first trajectory gene
    function selectGene(gene) {{
      if (activeTab === 'tab-volcano') {{
        highlightedVolcanoGene = gene;
        const d = volcanoData.find(x => x.gene === gene);
        if (d) setVolcanoGeneCard(d);
        renderVolcanoCanvas();
        renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
        return;
      }}
      setTrajectoryGene(gene);
    }}

    function setTrajectoryGene(gene) {{
      if (!GENES[gene]) return;
      trajGenes = [gene, ...trajGenes.filter(g => g !== gene)].slice(0, TRAJ_MAX_GENES);
      trajGenesUserSet = true;
      currentGene = gene;
      renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
      renderGenePlots();
    }}

    function suggestedGenes() {{
      if (activeTab === 'tab-volcano') {{
        const sig = volcanoData.filter(d => d.status !== 'NS').sort((a, b) => a.q - b.q || Math.abs(b.log2fc) - Math.abs(a.log2fc));
        const list = sig.length ? sig : [...volcanoData].sort((a, b) => a.p - b.p || Math.abs(b.log2fc) - Math.abs(a.log2fc)).slice(0, 40);
        return list.map(d => ({{
          gene: d.gene,
          meta: `${{d.log2fc >= 0 ? '▲ +' : '▼ '}}${{d.log2fc.toFixed(2)}}${{sig.length ? '' : ' (ns)'}}`,
          title: `log2FC ${{d.log2fc.toFixed(2)}}, p ${{d.p.toExponential(1)}}, q ${{d.q.toFixed(3)}}`
        }}));
      }}
      if (activeTab === 'tab-genes') {{
        const ranked = trajectoryChanges(trajPatient);
        const up = ranked.filter(r => r.delta > 0).slice(0, 25);
        const down = ranked.filter(r => r.delta < 0).sort((a, b) => a.delta - b.delta).slice(0, 25);
        return [...up, ...down].map(r => ({{
          gene: r.gene,
          meta: `${{r.delta > 0 ? '▲ +' : '▼ '}}${{r.delta.toFixed(1)}}`,
          title: `log2(TPM+1) change ${{r.from}} → ${{r.to}}: ${{r.delta.toFixed(2)}}`
        }}));
      }}
      return [];
    }}

    // OncoPrint Matrix
    // OncoPrint bins are relative to a reference; absolute TPM cut-offs flag every housekeeping/lineage gene
    const ONCO_BINS = [
      {{ key: 'vlow', color: '#2563eb' }},
      {{ key: 'low', color: '#93c5fd' }},
      {{ key: 'typ', color: '#e2e8f0' }},
      {{ key: 'high', color: '#fca5a5' }},
      {{ key: 'vhigh', color: '#dc2626' }},
    ];
    const ONCO_NOT_EXPRESSED = '#ffffff';

    // Same fallback as the Python outlier reference: origin+library, library, origin, whole cohort (>= 5 samples)
    function matchedReference(idx) {{
      const s = COHORT[idx];
      const usable = COHORT.map((c, i) => i).filter(i => usableExpression(COHORT[i]));
      const tiers = [
        i => COHORT[i].tissue_origin === s.tissue_origin && COHORT[i].library_type === s.library_type,
        i => COHORT[i].library_type === s.library_type,
        i => COHORT[i].tissue_origin === s.tissue_origin,
      ];
      for (const keep of tiers) {{
        const ref = usable.filter(keep);
        if (ref.length >= 5) return ref;
      }}
      return usable;
    }}

    // Leave-one-out z of log2(TPM+1) vs the matched reference
    function cohortZ(gene, idx) {{
      const vals = GENES[gene] || [];
      const ref = matchedReference(idx).filter(i => i !== idx);
      if (ref.length < 3) return null;
      const logs = ref.map(i => Math.log2((vals[i] || 0) + 1));
      const mean = logs.reduce((a, b) => a + b, 0) / logs.length;
      const sd = Math.sqrt(logs.reduce((a, b) => a + (b - mean) ** 2, 0) / (logs.length - 1));
      return (Math.log2((vals[idx] || 0) + 1) - mean) / Math.max(sd, 0.25);
    }}

    function oncoBin(mode, gene, idx) {{
      if (mode !== 'cohort') {{
        const p = hartwigPct(mode, gene, idx);
        if (p === null || p === undefined) return null;
        return p < 0.05 ? 0 : p < 0.25 ? 1 : p <= 0.75 ? 2 : p <= 0.95 ? 3 : 4;
      }}
      const z = cohortZ(gene, idx);
      if (z === null) return null;
      return z <= -2 ? 0 : z <= -1 ? 1 : z < 1 ? 2 : z < 2 ? 3 : 4;
    }}

    function renderOncoPrint() {{
      const wrapper = document.getElementById('oncoprint-wrapper');
      if (!wrapper) return;
      const mode = document.getElementById('oncoprint-mode').value;

      const labels = mode !== 'cohort'
        ? ['< 5th pctl', '5–25th', '25–75th', '75–95th', '> 95th pctl']
        : ['z ≤ −2', '−2 to −1', 'Typical (|z| < 1)', '1 to 2', 'z ≥ 2'];
      document.getElementById('oncoprint-legend').innerHTML =
        ONCO_BINS.map((b, i) => `<span><span style="color: ${{b.color}};">■</span> ${{labels[i]}}</span>`).join('') +
        `<span><span style="color: #cbd5e1;">□</span> Not expressed (&lt; 1 TPM)</span>`;
      document.getElementById('oncoprint-caption').innerText = mode === 'hartwig'
        ? "Percentile within the sample's Hartwig cancer type (Isofox gene_distribution; e.g. GBM / astrocytoma vs Nervous system, colorectal vs Colorectum). Hartwig is mostly adult metastatic disease, so primary and paediatric tumours may still differ from it."
        : mode === 'hartwig_all'
        ? 'Percentile within the Hartwig pan-cancer reference (Isofox PercentileCohort). Caution: lineage genes (e.g. glial markers in gliomas) read as high against a pan-cancer cohort; prefer the matched cancer type.'
        : 'Leave-one-out z-score of log2(TPM+1) against QC-passing samples with the same tissue origin and library type (falling back to same library, same origin, then whole cohort when fewer than 5).';

      const displayGenes = [
        'EGFR', 'SOX2', 'GFAP', 'OLIG2', 'PRAME', 'CTAG1B', 'MAGEA3',
        'CD274', 'PDCD1', 'CTLA4', 'CD8A', 'B2M', 'HLA-A', 'HLA-B', 'CEACAM5', 'EPCAM'
      ].filter(g => GENES[g]);

      let htmlHdr = `<div style="display: flex; flex-direction: column; gap: 4px;">`;
      let htmlGrid = `<div style="display: flex; flex-direction: column; gap: 4px;">`;

      displayGenes.forEach(gene => {{
        htmlHdr += `<div class="onco-row-hdr mono">${{gene}}</div>`;
        const vals = GENES[gene] || [];

        let rowCells = `<div class="onco-cells-row">`;
        COHORT.forEach((s, idx) => {{
          if (!usableExpression(s)) return;
          const tpm = vals[idx] || 0.0;
          let bg = ONCO_NOT_EXPRESSED;
          let detail = 'not expressed';
          if (tpm >= 1.0) {{
            const b = oncoBin(mode, gene, idx);
            if (b !== null) bg = ONCO_BINS[b].color;
            if (mode !== 'cohort') {{
              const p = hartwigPct(mode, gene, idx);
              const refName = mode === 'hartwig' ? s.hartwig_cancer_type : 'pan-cancer';
              detail = p === null || p === undefined ? 'no reference' : `${{refName}} ${{(p * 100).toFixed(0)}}th pctl`;
            }} else {{
              const z = cohortZ(gene, idx);
              detail = z === null ? 'reference too small' : `cohort z ${{z.toFixed(2)}}`;
            }}
          }}

          rowCells += `
            <div class="onco-cell" style="background: ${{bg}}; ${{bg === ONCO_NOT_EXPRESSED ? 'box-shadow: inset 0 0 0 1px #cbd5e1;' : ''}}"
                 title="${{s.tissue_id}} (${{s.tissue_origin}}): ${{gene}} = ${{tpm.toFixed(1)}} TPM · ${{detail}}"
                 onclick="setTrajectoryGene('${{gene}}'); openDrawer('${{s.tissue_id}}')"></div>
          `;
        }});
        rowCells += `</div>`;
        htmlGrid += rowCells;
      }});

      htmlHdr += `</div>`;
      htmlGrid += `</div>`;

      wrapper.innerHTML = htmlHdr + htmlGrid;
    }}

    // Gene Plots
    function renderGenePlots() {{
      renderTrajectory();  // sets currentGene to the first trajectory gene
      const vals = GENES[currentGene] || [];
      if (!vals.length) return;

      // Grouped Origin Plot
      const boxSvg = document.getElementById('gene-grouped-box-svg');
      if (boxSvg) {{
        const groups = {{}};
        COHORT.forEach((s, idx) => {{
          if (!usableExpression(s)) return;
          const orig = s.tissue_origin || 'Unknown';
          if (!groups[orig]) groups[orig] = [];
          groups[orig].push(vals[idx] || 0.0);
        }});

        let svgHtml = `<svg width="100%" height="100%" viewBox="0 0 450 230">`;
        const groupNames = Object.keys(groups);
        const medianOf = a => [...a].sort((x, y) => x - y)[Math.floor(a.length / 2)] || 0;
        const maxVal = Math.max(1, ...groupNames.map(g => medianOf(groups[g])));
        const colWidth = 380 / Math.max(1, groupNames.length);

        groupNames.forEach((g, i) => {{
          const gVals = groups[g].sort((a,b) => a-b);
          const median = gVals[Math.floor(gVals.length / 2)] || 0;
          const x = 40 + i * colWidth + colWidth/4;
          const y = 180 - (median / maxVal) * 140;

          svgHtml += `
            <rect x="${{x}}" y="${{y}}" width="${{colWidth/2}}" height="${{Math.max(4, (median/maxVal)*140)}}" fill="var(--brand-blue)" opacity="0.8" rx="3"/>
            <text x="${{x + colWidth/4}}" y="205" font-size="10" fill="var(--text-muted)" text-anchor="middle">${{g}}</text>
            <text x="${{x + colWidth/4}}" y="${{Math.max(16, y - 6)}}" font-size="11" fill="var(--text-main)" text-anchor="middle" font-weight="700">${{median.toFixed(1)}}</text>
          `;
        }});
        svgHtml += `<line x1="30" y1="180" x2="430" y2="180" stroke="var(--border-dark)" stroke-width="1.5"/>`;
        svgHtml += `</svg>`;
        boxSvg.innerHTML = svgHtml;
      }}

    }}

    // ---- Longitudinal trajectory ----
    function trajPatients() {{
      const by = {{}};
      COHORT.forEach((s, idx) => {{ if (usableExpression(s)) (by[s.patient_id] = by[s.patient_id] || []).push(idx); }});
      return Object.entries(by).filter(([, idxs]) => idxs.length >= 2).map(([pid, idxs]) => {{
        idxs.sort((a, b) => (COHORT[a].date_iso || '').localeCompare(COHORT[b].date_iso || '') || COHORT[a].tissue_id.localeCompare(COHORT[b].tissue_id));
        const lineages = {{}};
        idxs.forEach(i => {{ const l = COHORT[i].tissue_origin || 'Unknown'; lineages[l] = (lineages[l] || 0) + 1; }});
        const maxLineage = Math.max(...Object.values(lineages));
        const span = (Date.parse(COHORT[idxs[idxs.length - 1]].date_iso) - Date.parse(COHORT[idxs[0]].date_iso)) || 0;
        return {{ pid, idxs, lineages, maxLineage, span }};
      }}).sort((a, b) => b.maxLineage - a.maxLineage || b.span - a.span || a.pid.localeCompare(b.pid));
    }}

    // Per-gene log2(TPM+1) change from first to last sample of the patient's largest same-origin series
    function trajectoryChanges(pid) {{
      const p = trajPatients().find(x => x.pid === pid);
      if (!p) return [];
      const main = Object.entries(p.lineages).sort((a, b) => b[1] - a[1])[0][0];
      const series = p.maxLineage >= 2 ? p.idxs.filter(i => (COHORT[i].tissue_origin || 'Unknown') === main) : p.idxs;
      const first = series[0], last = series[series.length - 1];
      const out = [];
      for (const [gene, vals] of Object.entries(GENES)) {{
        if (NOISE_GENE_RE.test(gene)) continue;
        const a = vals[first] || 0, b = vals[last] || 0;
        if (Math.max(a, b) < 5) continue;
        out.push({{ gene, delta: Math.log2(b + 1) - Math.log2(a + 1), from: COHORT[first].tissue_id, to: COHORT[last].tissue_id }});
      }}
      return out.sort((x, y) => y.delta - x.delta);
    }}

    function onTrajPatientChange() {{
      trajPatient = document.getElementById('traj-patient').value;
      if (!trajGenesUserSet) {{
        const top = trajectoryChanges(trajPatient)[0];
        trajGenes = top ? [top.gene] : trajGenes;
        currentGene = trajGenes[0] || currentGene;
      }}
      renderGenePlots();
      if (activeTab === 'tab-genes') renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
    }}

    function addTrajGene() {{
      const input = document.getElementById('traj-gene-input');
      const gene = input.value.trim().toUpperCase();
      const match = GENES[gene] ? gene : Object.keys(GENES).find(g => g.toUpperCase() === gene);
      if (!match) {{ input.style.borderColor = 'var(--up-red)'; return; }}
      input.style.borderColor = '';
      input.value = '';
      if (!trajGenes.includes(match)) trajGenes = [...trajGenes, match].slice(-TRAJ_MAX_GENES);
      trajGenesUserSet = true;
      currentGene = trajGenes[0];
      renderGenePlots();
    }}

    function removeTrajGene(gene) {{
      if (trajGenes.length <= 1) return;
      trajGenes = trajGenes.filter(g => g !== gene);
      trajGenesUserSet = true;
      currentGene = trajGenes[0];
      renderGenePlots();
      if (activeTab === 'tab-genes') renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
    }}

    function renderTrajectory() {{
      const sel = document.getElementById('traj-patient');
      const holder = document.getElementById('gene-longitudinal-svg');
      if (!sel || !holder) return;

      const patients = trajPatients();
      sel.innerHTML = patients.map(p => {{
        const lin = Object.entries(p.lineages).map(([l, n]) => `${{l}} ×${{n}}`).join(', ');
        return `<option value="${{p.pid}}">${{p.pid}} (${{p.idxs.length}} samples: ${{lin}})</option>`;
      }}).join('');
      if (!patients.length) {{
        holder.innerHTML = '<p style="color: var(--text-light); text-align: center; padding: 40px;">No patient has two or more usable samples.</p>';
        return;
      }}
      if (!patients.some(p => p.pid === trajPatient)) {{
        trajPatient = patients[0].pid;
        if (!trajGenesUserSet) {{
          const top = trajectoryChanges(trajPatient)[0];
          if (top) trajGenes = [top.gene];
        }}
      }}
      sel.value = trajPatient;
      if (!trajGenes.length) trajGenes = [currentGene].filter(g => GENES[g]);
      currentGene = trajGenes[0];
      document.getElementById('active-gene-symbol').innerText = currentGene || '';

      document.getElementById('traj-gene-chips').innerHTML = trajGenes.map((g, i) => `
        <span class="traj-gene-chip" style="background: ${{TRAJ_COLORS[i]}};">${{g}}${{trajGenes.length > 1 ? `<button title="Remove" onclick="removeTrajGene('${{g}}')">×</button>` : ''}}</span>
      `).join('');
      const dl = document.getElementById('traj-gene-options');
      if (!dl.options.length) dl.innerHTML = Object.keys(GENES).sort().map(g => `<option value="${{g}}">`).join('');

      const p = patients.find(x => x.pid === trajPatient);
      const idxs = p.idxs;
      const excluded = COHORT.filter(s => s.patient_id === trajPatient && !usableExpression(s)).map(s => s.tissue_id);
      document.getElementById('traj-heading').innerText = `Longitudinal Trajectory: ${{trajPatient}}`;
      document.getElementById('traj-caption').innerText =
        `Samples in collection order; one line per tissue origin (different tumours are never joined).` +
        (excluded.length ? ` Excluded (QC fail): ${{excluded.join(', ')}}.` : '');

      // Layout
      const W = 900, H = 320, pad = {{ l: 56, r: 24, t: 20, b: 64 }};
      const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
      const xs = idxs.map((_, k) => pad.l + (idxs.length === 1 ? plotW / 2 : (k / (idxs.length - 1)) * plotW));
      const logv = (g, i) => Math.log2(((GENES[g] || [])[i] || 0) + 1);
      const yMax = Math.max(1, Math.ceil(Math.max(...trajGenes.flatMap(g => idxs.map(i => logv(g, i))))));
      const y = v => pad.t + plotH - (v / yMax) * plotH;
      const lineages = [...new Set(idxs.map(i => COHORT[i].tissue_origin || 'Unknown'))];
      const dashes = ['', '6 4', '2 4', '10 3 2 3'];

      let svg = `<svg width="100%" height="100%" viewBox="0 0 ${{W}} ${{H}}" preserveAspectRatio="xMidYMid meet" style="font-family: inherit;">`;
      // Axes and grid (log2(TPM+1), labelled in TPM)
      [0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000].map(t => [t, Math.log2(t + 1)]).filter(([, v]) => v <= yMax).forEach(([t, v]) => {{
        svg += `<line x1="${{pad.l}}" x2="${{W - pad.r}}" y1="${{y(v)}}" y2="${{y(v)}}" stroke="var(--border-color)" stroke-width="1"/>`;
        svg += `<text x="${{pad.l - 8}}" y="${{y(v) + 4}}" font-size="10" fill="var(--text-muted)" text-anchor="end">${{t.toLocaleString()}}</text>`;
      }});
      svg += `<text transform="translate(14, ${{pad.t + plotH / 2}}) rotate(-90)" font-size="10" fill="var(--text-muted)" text-anchor="middle">TPM (log2 scale)</text>`;
      idxs.forEach((i, k) => {{
        const s = COHORT[i];
        svg += `<line x1="${{xs[k]}}" x2="${{xs[k]}}" y1="${{pad.t}}" y2="${{pad.t + plotH}}" stroke="var(--border-color)" stroke-dasharray="2 3"/>`;
        svg += `<text x="${{xs[k]}}" y="${{H - pad.b + 16}}" font-size="10.5" fill="var(--text-main)" font-weight="700" text-anchor="middle">${{s.tissue_id}}</text>`;
        svg += `<text x="${{xs[k]}}" y="${{H - pad.b + 30}}" font-size="9.5" fill="var(--text-muted)" text-anchor="middle">${{s.date_iso || s.date_collected || '?'}} · ${{s.timepoint}}</text>`;
        svg += `<text x="${{xs[k]}}" y="${{H - pad.b + 43}}" font-size="9.5" fill="var(--text-muted)" text-anchor="middle">${{s.tissue_origin}} · ${{(s.path_dx || '').split(',')[0].slice(0, 28)}}</text>`;
      }});
      // One series per gene and tissue origin
      trajGenes.forEach((g, gi) => {{
        lineages.forEach((lin, li) => {{
          const pts = idxs.map((i, k) => ({{ i, k }})).filter(o => (COHORT[o.i].tissue_origin || 'Unknown') === lin);
          if (pts.length > 1) {{
            const d = pts.map((o, j) => `${{j ? 'L' : 'M'}} ${{xs[o.k]}} ${{y(logv(g, o.i))}}`).join(' ');
            svg += `<path d="${{d}}" fill="none" stroke="${{TRAJ_COLORS[gi]}}" stroke-width="2.5" stroke-dasharray="${{dashes[li % dashes.length]}}" stroke-linecap="round"/>`;
          }}
          pts.forEach(o => {{
            const tpm = (GENES[g] || [])[o.i] || 0;
            const s = COHORT[o.i];
            svg += `<circle cx="${{xs[o.k]}}" cy="${{y(logv(g, o.i))}}" r="5" fill="#ffffff" stroke="${{TRAJ_COLORS[gi]}}" stroke-width="2.5"><title>${{g}} · ${{s.tissue_id}} (${{s.date_iso}}, ${{s.path_dx}}): ${{tpm.toFixed(1)}} TPM</title></circle>`;
            svg += `<text x="${{xs[o.k] + 8}}" y="${{y(logv(g, o.i)) - 7}}" font-size="10" fill="${{TRAJ_COLORS[gi]}}" font-weight="700">${{tpm.toFixed(1)}}</text>`;
          }});
        }});
      }});
      svg += `<line x1="${{pad.l}}" x2="${{W - pad.r}}" y1="${{pad.t + plotH}}" y2="${{pad.t + plotH}}" stroke="var(--border-dark)" stroke-width="1.5"/>`;
      svg += `</svg>`;
      holder.innerHTML = svg;

      document.getElementById('traj-legend').innerHTML =
        trajGenes.map((g, i) => `<span><span style="color: ${{TRAJ_COLORS[i]}};">●</span> <b class="mono">${{g}}</b></span>`).join('') +
        (lineages.length > 1 ? lineages.map((l, li) => `<span><svg width="26" height="8"><line x1="0" x2="26" y1="4" y2="4" stroke="var(--text-muted)" stroke-width="2" stroke-dasharray="${{dashes[li % dashes.length]}}"/></svg> ${{l}}</span>`).join('') : '');
    }}

    // Render TCR Charts
    function renderTcrCharts() {{
      const container = document.getElementById('tcr-locus-chart');
      if (!container) return;

      let html = `<div style="display: flex; gap: 4px; height: 130px; align-items: flex-end; padding-bottom: 14px;">`;
      COHORT.forEach(s => {{
        const tcr = s.tcr_reads || 0;
        const bcr = s.bcr_reads || 0;
        const tot = Math.max(1, tcr + bcr);
        const tcrH = (tcr / tot) * 100;
        const bcrH = (bcr / tot) * 100;

        html += `
          <div style="flex: 1; display: flex; flex-direction: column; height: 100%; justify-content: flex-end; cursor: pointer;" title="${{s.tissue_id}}: TCR=${{tcr}}, BCR=${{bcr}}" onclick="openDrawer('${{s.tissue_id}}')">
            <div style="height: ${{bcrH}}px; background: var(--warning-amber); border-radius: 2px 2px 0 0;" title="BCR"></div>
            <div style="height: ${{tcrH}}px; background: var(--brand-blue); border-radius: 0 0 2px 2px;" title="TCR"></div>
          </div>
        `;
      }});
      html += `</div>`;
      container.innerHTML = html;
    }}

    // Clonotypes Table
    function renderClonotypesTable(items = CLONOTYPES) {{
      const tbody = document.getElementById('cdr3-table-body');
      if (!tbody) return;

      tbody.innerHTML = items.slice(0, 100).map(c => `
        <tr>
          <td class="mono"><b>${{c.sample_id}}</b></td>
          <td>${{c.patient_id}}</td>
          <td><span class="badge-pill blue">${{c.locus}}</span></td>
          <td class="mono">${{c.v_gene || '-'}}</td>
          <td class="mono">${{c.j_gene || '-'}}</td>
          <td class="mono" style="color: var(--brand-blue); font-weight: 700;">${{c.cdr3_aa}}</td>
          <td class="num">${{c.cdr3_len}}</td>
          <td class="num"><b>${{c.reads}}</b></td>
          <td>${{c.in_frame ? '<span class="badge-pill green">Productive</span>' : '<span class="badge-pill amber">Non-productive</span>'}}</td>
          <td class="mono" style="font-size: 0.74rem;">${{(c.shared_with || []).join(', ') || '-'}}</td>
        </tr>
      `).join('');
    }}

    function filterClonotypes(query = '') {{
      const locus = document.getElementById('cdr3-locus-select').value;
      const q = (query || document.getElementById('cdr3-filter').value || '').toLowerCase();

      const productive = document.getElementById('cdr3-productive').checked;
      const filtered = CLONOTYPES.filter(c => {{
        if (productive && !c.in_frame) return false;
        const matchLocus = locus === 'ALL' || c.locus === locus || c.receptor === locus;
        const matchQ = !q || c.cdr3_aa.toLowerCase().includes(q) || c.v_gene.toLowerCase().includes(q) || c.sample_id.toLowerCase().includes(q);
        return matchLocus && matchQ;
      }});
      renderClonotypesTable(filtered);
    }}

    // Fusions
    const CANDIDATE_MIN_SPLIT = 3;
    const FUSION_CLASS_STYLE = {{
      'Known / promiscuous': 'red', 'Distant / interchromosomal': 'purple',
      'Local (< 1 Mb)': 'gray', 'Paralog / pseudogene': 'gray', 'Immune receptor': 'blue'
    }};
    // Candidate tumour fusion: known/promiscuous or distant partner with split-read support on exon boundaries
    function isCandidateFusion(f) {{
      return (f.fusion_class === 'Known / promiscuous' || f.fusion_class === 'Distant / interchromosomal') &&
        f.split_reads >= CANDIDATE_MIN_SPLIT && f.exon_junction;
    }}

    // Distinct patients per gene pair (recurrent unknown pairs are usually alignment artefacts)
    const FUSION_PAIR_PATIENTS = {{}};
    FUSIONS.forEach(f => {{ (FUSION_PAIR_PATIENTS[f.gene_up + '>' + f.gene_down] = FUSION_PAIR_PATIENTS[f.gene_up + '>' + f.gene_down] || new Set()).add(f.patient_id); }});

    function renderFusionsTable() {{
      const tbody = document.getElementById('fusions-table-body');
      if (!tbody) return;
      const cls = document.getElementById('fusion-class-filter').value;
      const minSplit = parseInt(document.getElementById('fusion-min-split').value, 10);
      const exonOnly = document.getElementById('fusion-exon-only').checked;
      const q = document.getElementById('fusion-search').value.trim().toUpperCase();

      const rows = FUSIONS.filter(f => {{
        if (cls === 'CANDIDATE' && !['Known / promiscuous', 'Distant / interchromosomal'].includes(f.fusion_class)) return false;
        if (cls !== 'CANDIDATE' && cls !== 'ALL' && f.fusion_class !== cls) return false;
        if (f.split_reads < minSplit) return false;
        if (exonOnly && !f.exon_junction) return false;
        if (q && !(f.gene_up.toUpperCase().includes(q) || f.gene_down.toUpperCase().includes(q) || f.sample_id.toUpperCase().includes(q))) return false;
        return true;
      }}).sort((a, b) => (b.split_reads + b.spanning_reads) - (a.split_reads + a.spanning_reads));

      document.getElementById('fusion-catalog-caption').innerText =
        `${{rows.length.toLocaleString()}} of ${{FUSIONS.length.toLocaleString()}} chimeras match${{rows.length > 300 ? ' (top 300 shown)' : ''}}; sorted by supporting fragments.`;
      if (!rows.length) {{
        tbody.innerHTML = `<tr><td colspan="12" style="text-align: center; color: var(--text-light); padding: 24px;">No fusions match the current filters</td></tr>`;
        return;
      }}
      tbody.innerHTML = rows.slice(0, 300).map(f => `
        <tr>
          <td class="mono"><b>${{f.sample_id}}</b></td>
          <td>${{f.path_dx}}</td>
          <td class="mono" style="color: var(--brand-blue); font-weight: 700;">${{f.gene_up}}</td>
          <td class="mono" style="color: var(--purple-accent); font-weight: 700;">${{f.gene_down}}</td>
          <td class="mono" style="font-size: 0.76rem;">${{f.chr_up}}:${{f.pos_up}} → ${{f.chr_dn}}:${{f.pos_dn}}</td>
          <td class="num"><b>${{f.split_reads}}</b></td>
          <td class="num">${{f.spanning_reads}}</td>
          <td><span class="badge-pill ${{FUSION_CLASS_STYLE[f.fusion_class] || 'gray'}}">${{f.fusion_class}}</span></td>
          <td style="font-size: 0.74rem;">${{f.known_type === 'NONE' ? '-' : f.known_type}}</td>
          <td style="font-size: 0.74rem;" title="${{f.junc_up}} / ${{f.junc_dn}}">${{f.exon_junction ? 'Exon boundary' : 'Non-canonical'}}</td>
          <td class="num">${{FUSION_PAIR_PATIENTS[f.gene_up + '>' + f.gene_down].size}}</td>
          <td>${{GENES[f.gene_up] || GENES[f.gene_down] ? `<button class="btn-header" style="padding: 2px 8px; font-size: 0.72rem;" onclick="jumpToGene('${{GENES[f.gene_up] ? f.gene_up : f.gene_down}}')">Inspect</button>` : ''}}</td>
        </tr>
      `).join('');
    }}

    function renderFusionTcrView() {{
      // Per-sample candidate fusions vs T-cell repertoire
      const bySample = {{}};
      FUSIONS.filter(isCandidateFusion).forEach(f => {{ (bySample[f.sample_id] = bySample[f.sample_id] || []).push(f); }});
      const rec = {{}};
      RECEPTOR_MATCHES.forEach(r => {{
        const o = rec[r.sample_id] = rec[r.sample_id] || {{ n: 0, binders: 0 }};
        o.n++; if (r.hla_binder) o.binders++;
      }});
      const fmt = (v, d = 2) => v === null || v === undefined ? '-' : Number(v).toFixed(d);
      const rows = [...COHORT].sort((a, b) => (a.patient_id + a.date_iso).localeCompare(b.patient_id + b.date_iso));
      document.getElementById('fusion-tcr-body').innerHTML = rows.map(s => {{
        const fus = (bySample[s.tissue_id] || []).sort((a, b) => b.split_reads - a.split_reads);
        const top = fus.slice(0, 3).map(f => `${{f.gene_up}}→${{f.gene_down}} (${{f.split_reads}})`).join(', ');
        const tpm = tcrPerMillion(s);
        const r = rec[s.tissue_id];
        return `
          <tr>
            <td class="mono"><a href="#" onclick="openDrawer('${{s.tissue_id}}'); return false;" style="color: var(--brand-blue); text-decoration: none; font-weight: 700;">${{s.tissue_id}}</a></td>
            <td>${{s.path_dx}}</td>
            <td>${{qcBadge(s)}}</td>
            <td class="num"><b>${{fus.length}}</b></td>
            <td class="mono" style="font-size: 0.72rem; max-width: 320px; white-space: normal;">${{top || '-'}}</td>
            <td class="num">${{tpm === null ? '-' : tpm.toFixed(1)}}</td>
            <td class="num">${{s.tcr_clonotypes || 0}}</td>
            <td class="num">${{fmt(s.simpson_clonality)}}</td>
            <td class="num">${{s.top_clone_share === null || s.top_clone_share === undefined ? '-' : (s.top_clone_share * 100).toFixed(0) + '%'}}</td>
            <td class="num">${{s.shared_tcr_clones || 0}}</td>
            <td class="num">${{r ? `<b>${{r.n}}</b> (${{r.binders}} binders)` : '-'}}</td>
          </tr>`;
      }}).join('');

      // Fusion Neoantigen x TCR repertoire matches
      document.getElementById('receptor-match-count').innerText =
        `${{RECEPTOR_MATCHES.length}} neoantigen-TCR pairs, ${{RECEPTOR_MATCHES.filter(r => r.hla_binder).length}} presenting HLA binders (rank ≤ 2%)`;
      document.getElementById('receptor-match-body').innerHTML = RECEPTOR_MATCHES.map(r => `
        <tr>
          <td class="mono"><b>${{r.sample_id}}</b></td>
          <td class="mono" style="font-weight: 700;">${{r.gene_up}}→${{r.gene_down}}</td>
          <td><span class="badge-pill ${{r.is_oncogenic ? 'blue' : 'purple'}}">${{r.fusion_class}}</span></td>
          <td class="mono" style="font-size: 0.74rem;"><b>${{r.neoepitope}}</b></td>
          <td class="mono" style="font-size: 0.74rem;">${{r.hla_allele}}</td>
          <td class="num">${{r.hla_rank !== null ? `<span class="badge-pill ${{r.hla_binder ? 'green' : 'gray'}}">${{r.hla_rank}}%</span>` : '-'}}</td>
          <td class="mono" style="font-size: 0.72rem;">${{r.top_tcr ? `${{r.top_tcr.cdr3_aa}} (${{r.top_tcr.v_gene || '?'}}/${{r.top_tcr.j_gene || '?'}})` : '<span style="color: var(--text-light);">-</span>'}}</td>
          <td><span class="badge-pill ${{r.top_tcr && r.top_tcr.source && r.top_tcr.source.includes('MiXCR') ? 'teal' : 'gray'}}" style="font-size: 0.7rem;">${{r.top_tcr ? r.top_tcr.source : '-'}}</span></td>
          <td class="num">${{r.top_tcr ? r.top_tcr.reads : '-'}}</td>
          <td class="num"><b>${{r.top_tcr ? Number(r.top_tcr.binding_score).toFixed(2) : '-'}}</b></td>
        </tr>
      `).join('') || '<tr><td colspan="10" style="text-align: center; color: var(--text-light); padding: 24px;">No fusion neoantigens found</td></tr>';
    }}

    // RNA Sample QC
    function qcLabel(key) {{
      if (key === 'isofox_status') return 'Isofox QC status';
      const c = QC_CHECKS.find(c => c.key === key);
      return c ? c.label : key;
    }}

    function qcBadge(s) {{
      const st = s.qc ? s.qc.status : 'NO DATA';
      const cls = {{ PASS: 'green', WARN: 'amber', FAIL: 'red' }}[st] || 'gray';
      const why = s.qc ? [...s.qc.fails, ...s.qc.warns].map(qcLabel).join(', ') : '';
      return `<span class="badge-pill ${{cls}}" title="${{why}}">QC: ${{st}}</span>`;
    }}

    function setIncludeQcFail(checked) {{
      includeQcFail = checked;
      document.querySelectorAll('.include-qc-fail-toggle').forEach(el => {{ el.checked = checked; }});
      onVolcanoCategoryChange();
      renderGenePlots();
      renderOncoPrint();
    }}

    function fmtQcThreshold(c, t) {{
      if (t === null) return '-';
      const op = c.direction === 'min' ? '<' : '>';
      if (c.key.endsWith('_fragments')) return `${{op}} ${{(t / 1e6).toFixed(1)}}M`;
      if (c.key.endsWith('_pct') && c.key !== 'uniquely_mapped_pct') return `${{op}} ${{(t * 100).toFixed(0)}}%`;
      if (c.key === 'uniquely_mapped_pct') return `${{op}} ${{t}}%`;
      return `${{op}} ${{t.toLocaleString()}}`;
    }}

    // QC table sort state; default puts failures first
    let qcSort = {{ key: 'status', dir: 1 }};
    const QC_STATUS_ORDER = {{ FAIL: 0, WARN: 1, PASS: 2, 'NO DATA': 3 }};
    function qcSortValue(s, key) {{
      if (key === 'status') return QC_STATUS_ORDER[s.qc.status];
      if (['tissue_id', 'tissue_origin', 'library_type', 'tissue_format'].includes(key)) return s[key] || '';
      const v = (s.qc.metrics || {{}})[key];
      return v === null || v === undefined ? null : v;
    }}
    function sortQcTable(key) {{
      qcSort = qcSort.key === key ? {{ key, dir: -qcSort.dir }} : {{ key, dir: 1 }};
      renderQcView();
    }}

    function renderQcView() {{
      const counts = {{ PASS: 0, WARN: 0, FAIL: 0, 'NO DATA': 0 }};
      COHORT.forEach(s => {{ counts[s.qc.status] = (counts[s.qc.status] || 0) + 1; }});
      document.getElementById('qc-kpi-pass').innerText = counts.PASS;
      document.getElementById('qc-kpi-warn').innerText = counts.WARN;
      document.getElementById('qc-kpi-fail').innerText = counts.FAIL;
      const nBulk = COHORT.filter(s => s.library_type === 'Bulk').length;
      const nCap = COHORT.filter(s => s.library_type === 'Target Capture').length;
      document.getElementById('qc-kpi-libs').innerText = `${{nBulk}} / ${{nCap}}`;

      document.getElementById('qc-rules-body').innerHTML = QC_CHECKS.map(c => `
        <tr>
          <td style="font-weight: 700;">${{c.label}}</td>
          <td>${{fmtQcThreshold(c, c.fail)}}</td>
          <td>${{fmtQcThreshold(c, c.warn)}}</td>
          <td style="color: var(--text-muted);">${{c.desc}}</td>
        </tr>
      `).join('') + `<tr><td style="font-weight: 700;">Isofox QC status</td><td>FAIL_*</td><td>WARN_*</td><td style="color: var(--text-muted);">Isofox's own QC call (e.g. FAIL_LOW_COVERAGE, WARN_DUPLICATE_RATE)</td></tr>`;

      const rows = [...COHORT].sort((a, b) => {{
        const va = qcSortValue(a, qcSort.key), vb = qcSortValue(b, qcSort.key);
        if (va === null && vb !== null) return 1;   // missing values always last
        if (vb === null && va !== null) return -1;
        const c = typeof va === 'string' ? va.localeCompare(vb) : (va > vb) - (va < vb);
        return c * qcSort.dir || QC_STATUS_ORDER[a.qc.status] - QC_STATUS_ORDER[b.qc.status] || a.tissue_id.localeCompare(b.tissue_id);
      }});
      document.querySelectorAll('table.qc-table th[data-sort]').forEach(th => {{
        th.querySelector('.sort-ind').textContent = th.dataset.sort === qcSort.key ? (qcSort.dir > 0 ? '▲' : '▼') : '';
      }});
      const cell = (s, key, text) => {{
        const c = QC_CHECKS.find(c => c.key === key);
        const cls = s.qc.fails.includes(key) ? 'qc-cell-fail' : (s.qc.warns.includes(key) ? 'qc-cell-warn' : '');
        const rule = c ? `${{c.label}}: fail ${{fmtQcThreshold(c, c.fail)}}, warn ${{fmtQcThreshold(c, c.warn)}}` : '';
        return `<td class="num ${{cls}}" title="${{rule}}">${{text}}</td>`;
      }};
      const pct = v => (v === null || v === undefined) ? '-' : (v * 100).toFixed(1);
      const num = (v, d = 0) => (v === null || v === undefined) ? '-' : Number(v).toLocaleString(undefined, {{ maximumFractionDigits: d }});

      document.getElementById('qc-table-body').innerHTML = rows.map(s => {{
        const m = s.qc.metrics;
        if (s.qc.status === 'NO DATA') {{
          return `<tr class="qc-row-nodata"><td class="mono">${{s.tissue_id}}</td><td>${{s.tissue_origin}}</td><td>${{s.tissue_format}}</td><td>${{qcBadge(s)}}</td><td>${{s.library_type}}</td><td colspan="14">No Isofox output found</td></tr>`;
        }}
        const isoCls = s.qc.fails.includes('isofox_status') ? 'qc-cell-fail' : (s.qc.warns.includes('isofox_status') ? 'qc-cell-warn' : '');
        return `
          <tr>
            <td class="mono"><a href="#" onclick="openDrawer('${{s.tissue_id}}'); return false;" style="color: var(--brand-blue); text-decoration: none; font-weight: 700;">${{s.tissue_id}}</a></td>
            <td>${{s.tissue_origin}}</td>
            <td title="${{s.library_type}}">${{s.library_type === 'Target Capture' ? 'Capture' : s.library_type}}</td>
            <td>${{s.tissue_format}}</td>
            <td>${{qcBadge(s)}}</td>
            ${{cell(s, 'total_fragments', (m.total_fragments / 1e6).toFixed(1))}}
            ${{cell(s, 'unique_fragments', (m.unique_fragments / 1e6).toFixed(1))}}
            ${{cell(s, 'duplicate_pct', pct(m.duplicate_pct))}}
            ${{cell(s, 'genes_detected', num(m.genes_detected))}}
            ${{cell(s, 'spliced_pct', pct(m.spliced_pct))}}
            <td class="num">${{pct(m.unspliced_pct)}}</td>
            <td class="num">${{pct(m.chimeric_pct)}}</td>
            ${{cell(s, 'forward_strand_pct', pct(m.forward_strand_pct))}}
            ${{cell(s, 'uniquely_mapped_pct', m.uniquely_mapped_pct === null ? '-' : m.uniquely_mapped_pct.toFixed(1))}}
            ${{cell(s, 'top100_pct', pct(m.top100_pct))}}
            ${{cell(s, 'median_frag_length', num(m.median_frag_length))}}
            <td class="num">${{m.median_gc_ratio ? m.median_gc_ratio.toFixed(3) : '-'}}</td>
            <td class="qc-iso ${{isoCls}}" title="${{m.isofox_status || ''}}">${{(m.isofox_status || '-').replace(/^(FAIL|WARN)_/, '$1: ').replace(/_/g, ' ')}}</td>
          </tr>
        `;
      }}).join('');
    }}

    function exportQcCSV() {{
      const keys = ['total_fragments', 'unique_fragments', 'duplicate_pct', 'genes_detected', 'spliced_pct', 'unspliced_pct',
        'chimeric_pct', 'forward_strand_pct', 'uniquely_mapped_pct', 'multi_mapped_pct', 'top100_pct', 'median_frag_length',
        'median_gc_ratio', 'spliced_gene_count', 'isofox_status'];
      const headers = ['tissue_id', 'patient_id', 'tissue_origin', 'library_type', 'tissue_format', 'qc_status', ...keys, 'qc_fails', 'qc_warns'];
      const rows = COHORT.map(s => [
        s.tissue_id, s.patient_id, `"${{s.tissue_origin}}"`, `"${{s.library_type}}"`, `"${{s.tissue_format}}"`, s.qc.status,
        ...keys.map(k => s.qc.metrics[k] ?? ''),
        `"${{s.qc.fails.join(';')}}"`, `"${{s.qc.warns.join(';')}}"`
      ].join(','));
      downloadCSV('mmr_cohort_rna_qc.csv', [headers.join(','), ...rows].join('\\n'));
    }}

    // Patient Outlier Expression & Antigen Annotation
    const AG_STYLE = {{
      'Cancer-testis T-cell antigen': {{ cls: 'ag-ct', badge: 'purple' }},
      'T-cell antigen': {{ cls: 'ag-tcell', badge: 'blue' }},
      'Cancer-testis-like': {{ cls: 'ag-ctlike', badge: 'amber' }},
    }};

    function initOutlierControls() {{
      const sel = document.getElementById('outlier-sample');
      if (!sel) return;
      const samples = COHORT.filter(s => (OUTLIERS[s.tissue_id] || []).length)
        .sort((a, b) => (a.patient_id + a.tissue_id).localeCompare(b.patient_id + b.tissue_id));
      sel.innerHTML = samples.map(s =>
        `<option value="${{s.tissue_id}}">${{s.qc.status === 'FAIL' ? '⚠ QC FAIL · ' : ''}}${{s.patient_id}} · ${{s.tissue_id}} (${{s.path_dx}}, ${{s.timepoint}})</option>`
      ).join('') || '<option value="">No expression data</option>';
      renderOutlierView();
    }}

    function outlierRecords(sampleId) {{
      return (OUTLIERS[sampleId] || []).map(r => {{
        const ag = ANTIGENS[r[0]];
        return {{
          gene: r[0], tpm: r[1], z: r[2], pan_pct: r[3], ref_med: r[4], ref_fc: r[5], can_pct: r[6],
          outlier: r[7] === 1, noise: r[8] === 1,
          ag_class: ag ? ag[0] : '', ag_assays: ag ? ag[1] : 0, ag_refs: ag ? ag[2] : 0,
          ag_mhc: ag ? ag[3] : '', ag_dx: ag ? ag[4] : ''
        }};
      }});
    }}

    function filteredOutlierRecords() {{
      const sampleId = document.getElementById('outlier-sample').value;
      const mode = document.getElementById('outlier-mode').value;
      const agFilter = document.getElementById('outlier-antigen-filter').value;
      const sortKey = document.getElementById('outlier-sort').value;
      const hideNoise = document.getElementById('outlier-hide-noise').checked;
      const q = document.getElementById('outlier-search').value.trim().toUpperCase();

      let recs = outlierRecords(sampleId);
      if (mode === 'outliers') recs = recs.filter(r => r.outlier);
      else recs = recs.filter(r => r.ag_class && r.ag_class !== 'Cancer-testis-like');
      if (hideNoise) recs = recs.filter(r => !r.noise);
      if (agFilter === 'ANY') recs = recs.filter(r => r.ag_class);
      else if (agFilter === 'TCELL') recs = recs.filter(r => r.ag_class.includes('T-cell'));
      else if (agFilter === 'CT') recs = recs.filter(r => r.ag_class.startsWith('Cancer-testis'));
      else if (agFilter === 'NONE') recs = recs.filter(r => !r.ag_class);
      if (q) recs = recs.filter(r => r.gene.toUpperCase().includes(q));

      const val = r => sortKey === 'ref_fc' ? (r.ref_fc ?? -99) : r[sortKey];
      recs.sort((a, b) => val(b) - val(a));
      return {{ sampleId, recs, sortKey }};
    }}

    function renderOutlierView() {{
      const sel = document.getElementById('outlier-sample');
      if (!sel || !sel.value) return;
      const {{ sampleId, recs, sortKey }} = filteredOutlierRecords();

      const smp = COHORT.find(s => s.tissue_id === sampleId);
      const banner = document.getElementById('outlier-qc-banner');
      if (smp && smp.qc.status === 'FAIL') {{
        banner.style.display = 'block';
        banner.innerHTML = `<b>${{sampleId}} failed RNA QC</b> (${{smp.qc.fails.map(qcLabel).join(', ')}}). Outlier calls likely reflect library artefacts rather than tumour biology; this sample is excluded from the cohort reference.`;
      }} else {{
        banner.style.display = 'none';
      }}

      document.getElementById('outlier-ref-note').innerHTML = smp && smp.hartwig_cancer_type
        ? `Hartwig reference for ${{sampleId}}: <b>${{smp.hartwig_cancer_type}}</b> (matched by ${{smp.hartwig_matched_by}} from "${{smp.path_dx}}").`
        : 'No Hartwig cancer-type reference loaded; pan-cancer percentiles are used.';

      // KPIs always describe the sample's outlier set (noise filter respected)
      const hideNoise = document.getElementById('outlier-hide-noise').checked;
      const all = outlierRecords(sampleId).filter(r => r.outlier && (!hideNoise || !r.noise));
      const tcell = all.filter(r => r.ag_class.includes('T-cell'));
      const ct = all.filter(r => r.ag_class.startsWith('Cancer-testis'));
      document.getElementById('outlier-kpi-count').innerText = all.length;
      document.getElementById('outlier-kpi-tcell').innerText = tcell.length;
      document.getElementById('outlier-kpi-ct').innerText = ct.length;
      const topAg = [...tcell].sort((a, b) => (b.ref_fc ?? b.z) - (a.ref_fc ?? a.z))[0];
      document.getElementById('outlier-kpi-top').innerText = topAg ? topAg.gene : '-';
      document.getElementById('outlier-kpi-top-sub').innerText = topAg
        ? `${{topAg.tpm}} TPM · ${{topAg.ref_fc !== null ? 'log2FC ' + topAg.ref_fc : 'z ' + topAg.z}}`
        : 'No antigen outliers';

      renderOutlierBars(recs.slice(0, 30), sortKey);
      renderOutlierTable(recs);
    }}

    function renderOutlierBars(recs, sortKey) {{
      const el = document.getElementById('outlier-bars');
      if (!recs.length) {{
        el.innerHTML = '<p style="color: var(--text-light); font-size: 0.8rem; padding: 20px; text-align: center;">No genes match the current filters.</p>';
        return;
      }}
      const metric = r => sortKey === 'ref_fc' ? (r.ref_fc ?? 0) : r[sortKey];
      const label = {{ ref_fc: 'log2FC', z: 'z', tpm: 'TPM' }}[sortKey];
      const maxV = Math.max(...recs.map(r => Math.abs(metric(r))), 1e-6);
      el.innerHTML = recs.map(r => {{
        const st = AG_STYLE[r.ag_class] || {{ cls: 'ag-none' }};
        const w = Math.max(0, metric(r)) / maxV * 100;
        return `
          <div class="outlier-bar-row" onclick="jumpToGene('${{r.gene}}')" title="${{r.gene}}: ${{r.tpm}} TPM · z ${{r.z}} · cancer-type pctl ${{r.can_pct >= 0 ? (r.can_pct * 100).toFixed(1) : 'n/a'}} · pan-cancer pctl ${{r.pan_pct >= 0 ? (r.pan_pct * 100).toFixed(1) : 'n/a'}}${{r.ag_class ? ' · ' + r.ag_class : ''}}">
            <span class="mono" style="font-weight: 700;">${{r.gene}}</span>
            <div class="outlier-bar-track"><div class="outlier-bar-fill ${{st.cls}}" style="width: ${{w.toFixed(1)}}%;"></div></div>
            <span class="mono" style="text-align: right; color: var(--text-muted);">${{label}} <b style="color: var(--text-main);">${{metric(r).toLocaleString(undefined, {{ maximumFractionDigits: 2 }})}}</b></span>
          </div>
        `;
      }}).join('');
    }}

    function renderOutlierTable(recs) {{
      const fmtPct = v => v >= 0 ? (v * 100).toFixed(1) : '-';
      document.getElementById('outlier-table-count').innerText = `${{recs.length}} genes`;
      document.getElementById('outlier-table-body').innerHTML = recs.map(r => {{
        const st = AG_STYLE[r.ag_class];
        return `
          <tr>
            <td class="mono" style="font-weight: 700;"><a href="#" onclick="jumpToGene('${{r.gene}}'); return false;" style="color: var(--brand-blue); text-decoration: none;">${{r.gene}}</a></td>
            <td class="num">${{r.tpm.toLocaleString()}}</td>
            <td class="num" style="font-weight: ${{r.z >= 3 ? 800 : 400}};">${{r.z.toFixed(2)}}</td>
            <td class="num">${{fmtPct(r.pan_pct)}}</td>
            <td class="num">${{fmtPct(r.can_pct)}}</td>
            <td class="num">${{r.ref_med >= 0 ? r.ref_med : '-'}}</td>
            <td class="num" style="font-weight: ${{(r.ref_fc ?? 0) >= 2 ? 800 : 400}};">${{r.ref_fc !== null ? r.ref_fc.toFixed(2) : '-'}}</td>
            <td>${{r.outlier ? '<span class="badge-pill red">Outlier</span>' : '<span style="color: var(--text-light);">-</span>'}}</td>
            <td>${{st ? `<span class="badge-pill ${{st.badge}}">${{r.ag_class}}</span>` : '<span style="color: var(--text-light);">-</span>'}}</td>
            <td class="num">${{r.ag_assays ? `${{r.ag_assays}} / ${{r.ag_refs}}` : '-'}}</td>
            <td>${{r.ag_mhc || '-'}}</td>
            <td style="font-size: 0.72rem; color: var(--text-muted); max-width: 260px;">${{r.ag_dx || '-'}}</td>
          </tr>
        `;
      }}).join('');
    }}

    function jumpToGene(gene) {{
      if (!GENES[gene]) return;
      setTrajectoryGene(gene);
      switchView('tab-genes');
    }}

    function exportOutliersCSV() {{
      const {{ sampleId, recs }} = filteredOutlierRecords();
      const headers = ['sample_id', 'gene', 'tpm', 'cohort_z', 'pan_cancer_percentile', 'hartwig_ref_median_tpm',
        'log2fc_vs_hartwig_ref', 'cancer_type_percentile', 'is_outlier', 'antigen_class', 'iedb_tcell_assays',
        'iedb_references', 'mhc_classes', 'iedb_tumour_contexts'];
      const rows = recs.map(r => [
        sampleId, r.gene, r.tpm, r.z, r.pan_pct, r.ref_med, r.ref_fc ?? '', r.can_pct, r.outlier ? 1 : 0,
        `"${{r.ag_class}}"`, r.ag_assays, r.ag_refs, `"${{r.ag_mhc}}"`, `"${{r.ag_dx}}"`
      ].join(','));
      downloadCSV(`${{sampleId}}_outlier_antigens.csv`, [headers.join(','), ...rows].join('\\n'));
    }}

    // CUPPA Concordance & Visual Charts
    function evaluateCuppaConcordance(s) {{
      const cup = s.cuppa || {{}};
      const pred = cup.cancer_type || 'Unknown';
      if (pred === 'Unknown' || !cup.likelihood) {{
        return {{ status: 'PENDING', label: 'Pending', color: 'gray', isMatch: false }};
      }}

      const dx = (s.path_dx || '').toLowerCase();
      const orig = (s.tissue_origin || '').toLowerCase();
      const pLow = pred.toLowerCase();

      let isMatch = false;
      if ((orig.includes('cns') || dx.includes('gbm') || dx.includes('astrocytoma') || dx.includes('glioma')) && pLow.includes('cns')) {{
        isMatch = true;
      }} else if ((orig.includes('gi') || dx.includes('colorectal') || dx.includes('colon') || dx.includes('adenocarcinoma')) && (pLow.includes('color') || pLow.includes('intestine') || pLow.includes('appendix') || pLow.includes('gi'))) {{
        isMatch = true;
      }} else if ((orig.includes('haematological') || orig.includes('hematological') || dx.includes('leukemia') || dx.includes('lymph') || dx.includes('marrow')) && (pLow.includes('lymph') || pLow.includes('myeloid'))) {{
        isMatch = true;
      }} else if (dx.includes('neuroblastoma') && (pLow.includes('neuro') || pLow.includes('pnet') || pLow.includes('cns') || pLow.includes('lymph'))) {{
        isMatch = true;
      }}

      if (isMatch) {{
        return {{ status: 'CONCORDANT', label: '✓ Concordant', color: 'green', isMatch: true }};
      }} else {{
        return {{ status: 'DISCORDANT', label: '⚠ Divergent', color: 'amber', isMatch: false }};
      }}
    }}

    function renderCuppaView() {{
      renderCuppaTable();
      renderCuppaMatrix();
      renderCuppaSpectrum();
    }}

    function renderCuppaMatrix() {{
      const container = document.getElementById('cuppa-matrix-container');
      if (!container) return;

      const completed = COHORT.filter(s => (s.cuppa?.cancer_type || 'Unknown') !== 'Unknown');
      if (!completed.length) {{
        container.innerHTML = '<p style="color: var(--text-light); font-size: 0.8rem; padding: 20px; text-align: center;">No CUPPA predictions available yet.</p>';
        return;
      }}

      const clinCategories = ['CNS (Glioma / Astrocytoma)', 'GI (Colorectal Adenocarcinoma)', 'Other / Divergent'];
      const predCategories = ['CNS: Glioma', 'Colorectum/Small intestine/Appendix', 'Lymphoid tissue'];

      function getClinGroup(s) {{
        const dx = (s.path_dx || '');
        const orig = (s.tissue_origin || '');
        if (orig === 'CNS' || dx.includes('GBM') || dx.includes('Astrocytoma') || dx.includes('Glioma')) return 'CNS (Glioma / Astrocytoma)';
        if (orig === 'GI' || dx.includes('Adenocarcinoma') || dx.includes('Colorectal')) return 'GI (Colorectal Adenocarcinoma)';
        return 'Other / Divergent';
      }}

      const counts = {{}};
      clinCategories.forEach(c => {{
        counts[c] = {{}};
        predCategories.forEach(p => {{ counts[c][p] = []; }});
      }});

      completed.forEach(s => {{
        const cg = getClinGroup(s);
        const cp = s.cuppa.cancer_type;
        if (counts[cg] && counts[cg][cp]) {{
          counts[cg][cp].push(s.tissue_id);
        }} else if (counts[cg]) {{
          if (!counts[cg][cp]) counts[cg][cp] = [];
          counts[cg][cp].push(s.tissue_id);
        }}
      }});

      let html = `<table class="concordance-matrix-table">
        <thead>
          <tr>
            <th style="text-align: left;">Clinical Diagnosis</th>
            <th>CNS: Glioma</th>
            <th>Colorectum / GI</th>
            <th>Lymphoid Tissue</th>
            <th>Concordance</th>
          </tr>
        </thead>
        <tbody>
      `;

      clinCategories.forEach((cg, rIdx) => {{
        html += `<tr><td class="matrix-row-hdr">${{cg}}</td>`;
        let rowTotal = 0;
        let rowMatch = 0;

        predCategories.forEach((cp, cIdx) => {{
          const arr = counts[cg][cp] || [];
          rowTotal += arr.length;
          const isDiag = (rIdx === 0 && cIdx === 0) || (rIdx === 1 && cIdx === 1);
          if (isDiag) rowMatch += arr.length;

          let cellClass = 'matrix-cell-zero';
          if (arr.length > 0) {{
            cellClass = isDiag ? 'matrix-cell-match' : 'matrix-cell-divergent';
          }}

          html += `
            <td class="${{cellClass}}" title="${{cg}} → ${{cp}}: ${{arr.join(', ')}}">
              ${{arr.length > 0 ? `<b>${{arr.length}}</b>` : '-'}}
            </td>
          `;
        }});

        const pct = rowTotal > 0 ? (rowMatch / rowTotal * 100).toFixed(0) + '%' : '-';
        html += `<td style="font-weight: 800; color: ${{rowMatch === rowTotal ? 'var(--success-green)' : 'var(--text-main)'}};">${{pct}}</td></tr>`;
      }});

      html += `</tbody></table>`;
      container.innerHTML = html;
    }}

    function renderCuppaSpectrum() {{
      const container = document.getElementById('cuppa-spectrum-container');
      if (!container) return;

      const completed = COHORT.filter(s => (s.cuppa?.cancer_type || 'Unknown') !== 'Unknown');
      if (!completed.length) return;

      completed.sort((a,b) => (b.cuppa.likelihood || 0) - (a.cuppa.likelihood || 0));

      let concordantCount = 0;

      const rows = completed.map(s => {{
        const evalRes = evaluateCuppaConcordance(s);
        if (evalRes.isMatch) concordantCount++;
        const lik = (s.cuppa.likelihood || 0);
        const likPct = (lik * 100).toFixed(0);

        return `
          <div class="spectrum-bar-row" onclick="openDrawer('${{s.tissue_id}}')" title="${{s.tissue_id}} (${{s.path_dx}}): ${{s.cuppa.cancer_type}} (${{(lik*100).toFixed(1)}}%)">
            <span class="mono" style="width: 85px; font-weight: 700;">${{s.tissue_id}}</span>
            <span style="width: 130px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; font-size: 0.72rem; color: var(--text-muted);">${{s.path_dx}}</span>
            <div class="spectrum-track">
              <div class="spectrum-fill ${{evalRes.isMatch ? 'concordant' : 'discordant'}}" style="width: ${{likPct}}%;"></div>
            </div>
            <span style="width: 140px; font-size: 0.72rem; font-weight: 600; text-align: right; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
              ${{s.cuppa.cancer_type}} <b>${{likPct}}%</b>
            </span>
            <span class="badge-pill ${{evalRes.color}}" style="font-size: 0.68rem;">${{evalRes.label}}</span>
          </div>
        `;
      }}).join('');

      container.innerHTML = rows;

      const totalComp = completed.length;
      const concRate = totalComp > 0 ? (concordantCount / totalComp * 100).toFixed(1) + '%' : '-';
      const elKpi = document.getElementById('cuppa-kpi-concordance');
      if (elKpi) elKpi.innerText = concRate;
      const elSub = document.getElementById('cuppa-kpi-concordance-sub');
      if (elSub) elSub.innerText = `${{concordantCount}} / ${{totalComp}} Sequenced Cases`;

      // Lineage-specific concordance and median likelihood
      const lineageKpi = (origin, id) => {{
        const grp = completed.filter(s => s.tissue_origin === origin);
        const hits = grp.filter(s => evaluateCuppaConcordance(s).isMatch).length;
        document.getElementById(id).innerText = grp.length ? (hits / grp.length * 100).toFixed(1) + '%' : '-';
        document.getElementById(id + '-sub').innerText = `${{hits}} / ${{grp.length}} ${{origin}} cases`;
      }};
      lineageKpi('CNS', 'cuppa-kpi-cns');
      lineageKpi('GI', 'cuppa-kpi-gi');
      const liks = completed.map(s => s.cuppa.likelihood || 0).sort((a, b) => a - b);
      const mid = Math.floor(liks.length / 2);
      const med = liks.length % 2 ? liks[mid] : (liks[mid - 1] + liks[mid]) / 2;
      document.getElementById('cuppa-kpi-median-lik').innerText = (med * 100).toFixed(1) + '%';
    }}

    function renderCuppaTable() {{
      const tbody = document.getElementById('cuppa-table-body');
      if (!tbody) return;

      tbody.innerHTML = COHORT.map(s => {{
        const cup = s.cuppa || {{}};
        const pred = cup.cancer_type || 'Unknown';
        const lik = cup.likelihood ? (cup.likelihood * 100).toFixed(1) + '%' : '-';
        const sec = (cup.predictions && cup.predictions[1]) ? `${{cup.predictions[1].cancer_type}} (${{(cup.predictions[1].likelihood*100).toFixed(0)}}%)` : '-';
        const evalRes = evaluateCuppaConcordance(s);

        return `
          <tr>
            <td class="mono"><b>${{s.tissue_id}}</b></td>
            <td>${{s.path_dx}}</td>
            <td><span class="badge-pill blue">${{s.tissue_origin}}</span></td>
            <td style="font-weight: 700; color: var(--text-main);">${{pred}}</td>
            <td class="num"><b>${{lik}}</b></td>
            <td style="color: var(--text-muted);">${{sec}}</td>
            <td><span class="badge-pill ${{evalRes.color}}">${{evalRes.label}}</span></td>
            <td><button class="btn-header" style="padding: 2px 8px; font-size: 0.72rem;" onclick="openDrawer('${{s.tissue_id}}')">Inspect</button></td>
          </tr>
        `;
      }}).join('');
    }}

    function exportCuppaCSV() {{
      const headers = ['tissue_id', 'patient_id', 'path_dx', 'tissue_origin', 'cuppa_pred_class_1', 'likelihood_1', 'cuppa_pred_class_2', 'likelihood_2', 'concordance'];
      const rows = COHORT.map(s => {{
        const cup = s.cuppa || {{}};
        const preds = cup.predictions || [];
        const p1 = preds[0] || {{}};
        const p2 = preds[1] || {{}};
        const ev = evaluateCuppaConcordance(s);
        return [
          `"${{s.tissue_id}}"`,
          `"${{s.patient_id}}"`,
          `"${{s.path_dx}}"`,
          `"${{s.tissue_origin}}"`,
          `"${{p1.cancer_type || ''}}"`,
          (p1.likelihood !== undefined ? p1.likelihood : ''),
          `"${{p2.cancer_type || ''}}"`,
          (p2.likelihood !== undefined ? p2.likelihood : ''),
          `"${{ev.label}}"`
        ].join(',');
      }});
      downloadCSV('mmr_cohort_cuppa_concordance.csv', [headers.join(','), ...rows].join('\\n'));
    }}

    // Samples Table
    function renderSamplesTable() {{
      const tbody = document.getElementById('samples-table-body');
      if (!tbody) return;

      tbody.innerHTML = COHORT.map(s => `
        <tr>
          <td class="mono"><b>${{s.tissue_id}}</b></td>
          <td>${{s.patient_id}}</td>
          <td>${{s.path_dx}}</td>
          <td><span class="badge-pill blue">${{s.tissue_origin}}</span></td>
          <td>${{s.library_type}}</td>
          <td>${{s.tissue_format}}</td>
          <td>${{qcBadge(s)}}</td>
          <td class="mono">${{(s.hla_a || []).join(', ')}}</td>
          <td class="mono">${{(s.hla_b || []).join(', ')}}</td>
          <td class="mono">${{(s.hla_c || []).join(', ')}}</td>
          <td class="num"><b>${{s.tcr_reads || 0}}</b></td>
          <td class="num">${{s.simpson_clonality === null || s.simpson_clonality === undefined ? '-' : s.simpson_clonality.toFixed(2)}}</td>
          <td><button class="btn-header" style="padding: 2px 8px; font-size: 0.72rem;" onclick="openDrawer('${{s.tissue_id}}')">Inspect</button></td>
        </tr>
      `).join('');
    }}

    // Slide-out Drawer
    function openDrawer(sampleId) {{
      const sample = COHORT.find(s => s.tissue_id === sampleId);
      if (!sample) return;

      document.getElementById('drawer-sample-id').innerText = sample.tissue_id;
      document.getElementById('drawer-patient-id').innerText = `${{sample.patient_id}} • ${{sample.path_dx}} (${{sample.timepoint}})`;

      const topClones = (CLONOTYPES.filter(c => c.sample_id === sample.tissue_id)).slice(0, 8);
      const summ = sample.isofox_summary || {{}};

      let content = `
        <div style="display: flex; gap: 6px; margin-bottom: 16px; flex-wrap: wrap;">
          <span class="badge-pill blue">${{sample.tissue_origin}}</span>
          <span class="badge-pill amber">${{sample.tissue_format}}</span>
          ${{sample.hartwig_cancer_type ? `<span class="badge-pill" title="Hartwig expression reference">Hartwig ref: ${{sample.hartwig_cancer_type}}</span>` : ''}}
          ${{qcBadge(sample)}}
        </div>

        <div style="margin-bottom: 16px;">
          <div style="font-size: 0.74rem; font-weight: 800; text-transform: uppercase; color: var(--text-muted); margin-bottom: 6px;">HLA Class I Genotype (4-Digit)</div>
          <div style="display: grid; grid-template-columns: repeat(3, 1fr); gap: 6px;">
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);">HLA-A</div>
              <div class="mono" style="font-weight: 700; color: var(--brand-blue);">${{(sample.hla_a || []).join('<br>') || '-'}}</div>
            </div>
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);">HLA-B</div>
              <div class="mono" style="font-weight: 700; color: var(--purple-accent);">${{(sample.hla_b || []).join('<br>') || '-'}}</div>
            </div>
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);">HLA-C</div>
              <div class="mono" style="font-weight: 700; color: #0284c7;">${{(sample.hla_c || []).join('<br>') || '-'}}</div>
            </div>
          </div>
        </div>

        <div style="margin-bottom: 16px;">
          <div style="font-size: 0.74rem; font-weight: 800; text-transform: uppercase; color: var(--text-muted); margin-bottom: 6px;">RNA Library & Infiltration</div>
          <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 8px;">
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);">Library Fragments</div>
              <div style="font-size: 1.1rem; font-weight: 800; color: var(--text-main);">${{(summ.total_fragments || 0).toLocaleString()}}</div>
            </div>
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);">Expressed Genes</div>
              <div style="font-size: 1.1rem; font-weight: 800; color: var(--text-main);">${{(summ.spliced_gene_count || 0).toLocaleString()}}</div>
            </div>
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);">TCR Reads / Clonality</div>
              <div style="font-size: 1.1rem; font-weight: 800; color: var(--text-main);">${{sample.tcr_reads || 0}} <span style="font-size: 0.78rem; color: var(--brand-blue);">(${{sample.simpson_clonality === null || sample.simpson_clonality === undefined ? `${{sample.tcr_clonotypes || 0}} clones` : 'clonality ' + sample.simpson_clonality.toFixed(2)}})</span></div>
            </div>
            <div style="background: var(--bg-subtle); padding: 8px; border-radius: 4px; border: 1px solid var(--border-color);">
              <div style="font-size: 0.68rem; color: var(--text-muted);" title="Cytolytic activity (Rooney et al., Cell 2015): mean log2(TPM+1) of GZMA and PRF1">Cytolytic Score (GZMA, PRF1)</div>
              <div style="font-size: 1.1rem; font-weight: 800; color: var(--success-green);">${{sample.cytolytic_score ?? '-'}}</div>
            </div>
          </div>
        </div>

        <div style="margin-bottom: 16px;">
          <div style="font-size: 0.74rem; font-weight: 800; text-transform: uppercase; color: var(--text-muted); margin-bottom: 6px;">Top Reconstructed CDR3 Clones</div>
          ${{topClones.length ? `
            <table class="data-table" style="font-size: 0.76rem;">
              <thead>
                <tr><th>Locus</th><th>V-J</th><th>CDR3</th><th class="num">Reads</th></tr>
              </thead>
              <tbody>
                ${{topClones.map(c => `
                  <tr>
                    <td><span class="badge-pill blue">${{c.locus}}</span></td>
                    <td class="mono" style="font-size: 0.7rem;">${{c.v_gene}} / ${{c.j_gene}}</td>
                    <td class="mono" style="color: var(--brand-blue); font-weight: 700;">${{c.cdr3_aa}}</td>
                    <td class="num"><b>${{c.reads}}</b></td>
                  </tr>
                `).join('')}}
              </tbody>
            </table>
          ` : '<p style="color: var(--text-light); font-size: 0.8rem;">No CDR3 clonotypes recovered for this sample</p>'}}
        </div>
      `;

      document.getElementById('drawer-body').innerHTML = content;
      document.getElementById('inspector-drawer').classList.add('open');
      document.getElementById('drawer-scrim').classList.add('open');
    }}

    function closeDrawer() {{
      document.getElementById('inspector-drawer').classList.remove('open');
      document.getElementById('drawer-scrim').classList.remove('open');
    }}

    // -------------------------------------------------------------------------
    // Interactive Volcano Plot Engine (Welch's t-test)
    // -------------------------------------------------------------------------

    function initVolcanoControls() {{
      onVolcanoCategoryChange();
    }}

    function getSampleCategoryValue(s, catKey) {{
      if (catKey === 'tissue_origin') return s.tissue_origin || 'Unknown';
      if (catKey === 'path_dx') return s.path_dx || 'Unknown';
      if (catKey === 'timepoint') return s.timepoint || 'Unknown';
      if (catKey === 'tissue_format') return s.tissue_format || 'Unknown';
      if (catKey === 'library_type') return s.library_type || 'Unknown';
      if (catKey === 'cuppa_pred') return (s.cuppa && s.cuppa.cancer_type) || 'Unknown';
      if (catKey === 'tcr_group') {{
        const v = tcrPerMillion(s);
        if (v === null) return 'Unknown';
        return v >= eligibleMedian(tcrPerMillion) ? 'High TCR (≥ median)' : 'Low TCR (< median)';
      }}
      if (catKey === 'cytolytic_group') return (s.cytolytic_score >= eligibleMedian(x => x.cytolytic_score)) ? 'High CYT (≥ median)' : 'Low CYT (< median)';
      return 'Unknown';
    }}

    // Group 2 option meaning "every eligible sample not in Group 1"
    const VOLCANO_REST = '__REST__';
    const isHlaCategory = catKey => /^hla_[abc]$/.test(catKey);

    // Category values per sample; HLA loci give the set of carried alleles (two per sample, so groups can overlap)
    function getSampleCategoryValues(s, catKey) {{
      if (isHlaCategory(catKey)) return [...new Set(s[catKey] || [])];
      return [getSampleCategoryValue(s, catKey)];
    }}

    // Sample indices for Group 1 / Group 2; samples carrying both compared values are dropped from both
    function volcanoGroups(catKey, g1, g2) {{
      const g1Idx = [], g2Idx = [];
      let overlap = 0;
      COHORT.forEach((s, idx) => {{
        if (!volcanoEligible(s)) return;
        const vals = getSampleCategoryValues(s, catKey);
        if (!vals.length) return;  // e.g. no HLA typing
        const in1 = vals.includes(g1);
        const in2 = g2 === VOLCANO_REST ? !in1 : vals.includes(g2);
        if (in1 && in2) overlap++;
        else if (in1) g1Idx.push(idx);
        else if (in2) g2Idx.push(idx);
      }});
      return {{ g1Idx, g2Idx, overlap }};
    }}

    const nPatients = idxs => new Set(idxs.map(i => COHORT[i].patient_id)).size;

    function matchesVolcanoFormat(s) {{
      const lib = document.getElementById('volcano-library').value;
      if (lib !== 'ALL' && s.library_type !== lib) return false;
      const fmt = document.getElementById('volcano-format').value;
      const tf = (s.tissue_format || '').toUpperCase();
      if (fmt === 'FROZEN') return tf.includes('FROZEN');
      if (fmt === 'FFPE') return tf.startsWith('FFPE');
      return true;
    }}

    // Samples without Isofox output would otherwise enter the test as all-zero expression
    function volcanoEligible(s) {{
      return usableExpression(s) && matchesVolcanoFormat(s);
    }}

    function onVolcanoCategoryChange() {{
      const catKey = document.getElementById('volcano-category').value;
      const prev1 = document.getElementById('volcano-group1').value;
      const prev2 = document.getElementById('volcano-group2').value;
      const members = {{}};

      COHORT.forEach((s, idx) => {{
        if (!volcanoEligible(s)) return;
        getSampleCategoryValues(s, catKey).forEach(v => {{ (members[v] = members[v] || []).push(idx); }});
      }});

      const sortedVals = Object.keys(members).sort((a, b) => members[b].length - members[a].length || a.localeCompare(b));
      const label = v => {{
        const n = members[v].length, k = nPatients(members[v]);
        return `${{v}} (n=${{n}}${{k < n ? `, ${{k}} pt${{k > 1 ? 's' : ''}}` : ''}})`;
      }};
      const sel1 = document.getElementById('volcano-group1');
      const sel2 = document.getElementById('volcano-group2');

      // HLA defaults to carriers vs everyone else; other categories to the two largest groups
      const hla = isHlaCategory(catKey);
      sel1.innerHTML = sortedVals.map(v => `<option value="${{v}}">${{label(v)}}</option>`).join('');
      sel2.innerHTML = `<option value="${{VOLCANO_REST}}" ${{hla ? 'selected' : ''}}>All other samples</option>` +
        sortedVals.map((v, i) => `<option value="${{v}}" ${{!hla && i === 1 ? 'selected' : ''}}>${{label(v)}}</option>`).join('');
      // Keep the user's groups when only the format filter changed
      if (members[prev1] && (members[prev2] || prev2 === VOLCANO_REST) && prev1 !== prev2) {{
        sel1.value = prev1;
        sel2.value = prev2;
      }}

      updateVolcanoPlot();
    }}

    // Hill (1970) Student's t distribution approximation
    function erfApprox(x) {{
      const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
      const sign = x >= 0 ? 1 : -1;
      const ax = Math.abs(x);
      const t = 1.0 / (1.0 + p * ax);
      const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-ax * ax);
      return sign * y;
    }}

    function studentTPValue(t, df) {{
      t = Math.abs(t);
      if (df <= 0 || isNaN(df)) return 1.0;
      const a = df - 0.5;
      const b = 48.0 * a * a;
      const z2 = a * Math.log(1.0 + (t * t) / df);
      const z = Math.sqrt(Math.max(0.0, z2));
      const cf = z * (1.0 - (z2 + 3.0) / b);
      const p = 1.0 - erfApprox(cf / Math.SQRT2);
      return Math.max(1e-12, Math.min(1.0, p));
    }}

    // Welch's t-test on log2(TPM+1); m1/m2 are reported as mean TPM
    function calculateWelch(tpm1, tpm2) {{
      const n1 = tpm1.length;
      const n2 = tpm2.length;
      if (n1 < 1 || n2 < 1) return {{ p: 1.0, log2fc: 0.0, m1: 0, m2: 0 }};

      const tm1 = tpm1.reduce((a,b) => a+b, 0) / n1;
      const tm2 = tpm2.reduce((a,b) => a+b, 0) / n2;
      const vals1 = tpm1.map(v => Math.log2(v + 1));
      const vals2 = tpm2.map(v => Math.log2(v + 1));
      const m1 = vals1.reduce((a,b) => a+b, 0) / n1;
      const m2 = vals2.reduce((a,b) => a+b, 0) / n2;
      const log2fc = m1 - m2;

      if (n1 < 2 || n2 < 2) {{
        return {{ p: 1.0, log2fc, m1: tm1, m2: tm2 }};
      }}

      const v1 = vals1.reduce((a,b) => a + Math.pow(b - m1, 2), 0) / (n1 - 1);
      const v2 = vals2.reduce((a,b) => a + Math.pow(b - m2, 2), 0) / (n2 - 1);
      const se = Math.sqrt((v1 / n1) + (v2 / n2));

      if (se < 1e-9) {{
        return {{ p: 1.0, log2fc, m1: tm1, m2: tm2 }};
      }}

      const t = Math.abs(m1 - m2) / se;
      const num = Math.pow((v1 / n1) + (v2 / n2), 2);
      const den = (Math.pow(v1 / n1, 2) / (n1 - 1)) + (Math.pow(v2 / n2, 2) / (n2 - 1));
      const df = Math.max(1, num / den);

      const p = studentTPValue(t, df);
      return {{ p, log2fc, m1: tm1, m2: tm2 }};
    }}

    // Benjamini-Hochberg adjustment, written onto each record as .q
    function applyBH(records) {{
      const sorted = [...records].sort((a, b) => a.p - b.p);
      const m = sorted.length;
      let running = 1.0;
      for (let i = m - 1; i >= 0; i--) {{
        running = Math.min(running, sorted[i].p * m / (i + 1));
        sorted[i].q = running;
      }}
    }}

    function updateVolcanoPlot() {{
      const catKey = document.getElementById('volcano-category').value;
      const g1 = document.getElementById('volcano-group1').value;
      const g2 = document.getElementById('volcano-group2').value;
      const fcCutoff = parseFloat(document.getElementById('volcano-fc-cutoff').value);
      const pCutoff = parseFloat(document.getElementById('volcano-p-cutoff').value);
      const useFdr = document.getElementById('volcano-p-type').value === 'fdr';

      const {{ g1Idx, g2Idx, overlap }} = volcanoGroups(catKey, g1, g2);
      const g2Name = g2 === VOLCANO_REST ? `not ${{g1}}` : g2;

      volcanoData = [];
      let upCount = 0;
      let downCount = 0;

      if (g1 === g2 || !g1 || !g2 || !g1Idx.length || !g2Idx.length) {{
        document.getElementById('volcano-counts-badge').innerText =
          '⚠ Only one group remains under the current filters: choose another comparison or relax the library / preservation filters.';
        renderVolcanoCanvas();
        renderVolcanoTopTable();
        return;
      }}

      for (const gene of VOLCANO_GENES) {{
        const vals = GENES[gene];
        if (!vals) continue;
        const v1 = g1Idx.map(i => vals[i] || 0.0);
        const v2 = g2Idx.map(i => vals[i] || 0.0);

        const res = calculateWelch(v1, v2);
        // Skip genes not expressed in either group
        if (Math.max(res.m1, res.m2) < 1.0) continue;

        volcanoData.push({{
          gene,
          log2fc: res.log2fc,
          p: res.p,
          negLogP: -Math.log10(Math.max(1e-12, res.p)),
          m1: res.m1,
          m2: res.m2,
          status: 'NS'
        }});
      }}

      applyBH(volcanoData);
      volcanoData.forEach(d => {{
        const isSig = (useFdr ? d.q : d.p) <= pCutoff;
        if (isSig && d.log2fc >= fcCutoff) {{
          d.status = 'UP';
          upCount++;
        }} else if (isSig && d.log2fc <= -fcCutoff) {{
          d.status = 'DOWN';
          downCount++;
        }}
      }});

      // Horizontal guide at the raw p corresponding to the chosen cutoff
      if (useFdr) {{
        const passing = volcanoData.filter(d => d.q <= pCutoff);
        volcanoPThresholdNegLog = passing.length
          ? Math.min(...passing.map(d => d.negLogP))
          : Math.max(0, ...volcanoData.map(d => d.negLogP)) + 0.5;
      }} else {{
        volcanoPThresholdNegLog = -Math.log10(pCutoff);
      }}

      const libVal = document.getElementById('volcano-library').value;
      const fmtLabel = (libVal === 'ALL' ? '' : libVal + ', ') +
        {{ ALL: 'all formats', FROZEN: 'Fresh Frozen', FFPE: 'FFPE' }}[document.getElementById('volcano-format').value];
      const smallN = (g1Idx.length < 3 || g2Idx.length < 3) ? ' | ⚠ n < 3 in a group: tests underpowered' : '';
      // Flag when the two groups are separated perfectly by preservation or library type (batch confounding)
      const confounders = [['tissue_format', 'preservation'], ['library_type', 'library type']].filter(([key]) => {{
        if (key === catKey) return false;
        const a = new Set(g1Idx.map(i => COHORT[i][key]));
        const b = new Set(g2Idx.map(i => COHORT[i][key]));
        return [...a].every(v => !b.has(v));
      }}).map(([, label]) => label);
      const confNote = confounders.length ? ` | ⚠ Fully confounded with ${{confounders.join(' and ')}}` : '';
      const splitNote = catKey === 'tcr_group' ? ` | Median split: ${{eligibleMedian(tcrPerMillion).toFixed(1)}} TCR reads / M fragments`
        : catKey === 'cytolytic_group' ? ` | Median split: CYT ${{eligibleMedian(x => x.cytolytic_score).toFixed(2)}}` : '';
      const p1 = nPatients(g1Idx), p2 = nPatients(g2Idx);
      const nTxt = (n, k) => k < n ? `n=${{n}} (${{k}} pt${{k > 1 ? 's' : ''}})` : `n=${{n}}`;
      const overlapNote = overlap ? ` | ${{overlap}} sample${{overlap > 1 ? 's' : ''}} carrying both excluded` : '';
      // Samples from one patient are not independent (HLA is identical across them), so p-values are optimistic
      const pseudoNote = (p1 < g1Idx.length || p2 < g2Idx.length) ? ' | ⚠ Multiple samples per patient: effective n is the patient count' : '';
      // Default highlight: the most significant gene, unless a significant gene is already highlighted
      const sigGenes = volcanoData.filter(d => d.status !== 'NS').sort((a, b) => a.q - b.q);
      if (sigGenes.length && !sigGenes.some(d => d.gene === highlightedVolcanoGene)) {{
        highlightedVolcanoGene = sigGenes[0].gene;
        setVolcanoGeneCard(sigGenes[0]);
      }}
      if (activeTab === 'tab-volcano') renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
      document.getElementById('volcano-counts-badge').innerText =
        `${{fmtLabel}}: ${{nTxt(g1Idx.length, p1)}} vs ${{nTxt(g2Idx.length, p2)}} | Genes tested: ${{volcanoData.length.toLocaleString()}} | Enriched in ${{g1}}: ${{upCount}} | Enriched in ${{g2Name}}: ${{downCount}}${{splitNote}}${{smallN}}${{overlapNote}}${{pseudoNote}}${{confNote}}`;

      renderVolcanoCanvas();
      renderVolcanoTopTable();
    }}

    function renderVolcanoCanvas() {{
      const canvas = document.getElementById('volcano-canvas');
      if (!canvas) return;

      const rect = canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      canvas.width = (rect.width || 600) * dpr;
      canvas.height = 430 * dpr;

      const ctx = canvas.getContext('2d');
      ctx.scale(dpr, dpr);

      const w = rect.width || 600;
      const h = 430;
      const pad = {{ left: 45, right: 25, top: 25, bottom: 40 }};
      const plotW = w - pad.left - pad.right;
      const plotH = h - pad.top - pad.bottom;

      ctx.clearRect(0, 0, w, h);

      // Coordinate scaling
      const maxFc = Math.max(4.0, ...volcanoData.map(d => Math.abs(d.log2fc)));
      const maxNegLogP = Math.max(3.5, ...volcanoData.map(d => d.negLogP));

      const xMin = -maxFc * 1.1;
      const xMax = maxFc * 1.1;
      const yMin = 0;
      const yMax = maxNegLogP * 1.15;

      const mapX = fc => pad.left + ((fc - xMin) / (xMax - xMin)) * plotW;
      const mapY = nlp => pad.top + plotH - ((nlp - yMin) / (yMax - yMin)) * plotH;

      // Draw gridlines
      ctx.strokeStyle = '#e2e8f0';
      ctx.lineWidth = 1;
      ctx.beginPath();

      // Zero FC
      const x0 = mapX(0);
      ctx.moveTo(x0, pad.top); ctx.lineTo(x0, pad.top + plotH);

      // FC Cutoffs
      const fcCut = parseFloat(document.getElementById('volcano-fc-cutoff').value);
      const pCut = parseFloat(document.getElementById('volcano-p-cutoff').value);
      const negLogPCut = volcanoPThresholdNegLog;

      ctx.stroke();

      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = '#cbd5e1';
      ctx.beginPath();
      const xCutR = mapX(fcCut);
      const xCutL = mapX(-fcCut);
      ctx.moveTo(xCutR, pad.top); ctx.lineTo(xCutR, pad.top + plotH);
      ctx.moveTo(xCutL, pad.top); ctx.lineTo(xCutL, pad.top + plotH);

      const yCut = mapY(negLogPCut);
      ctx.moveTo(pad.left, yCut); ctx.lineTo(pad.left + plotW, yCut);
      ctx.stroke();
      ctx.restore();

      // Axes lines
      ctx.strokeStyle = '#94a3b8';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(pad.left, pad.top + plotH); ctx.lineTo(pad.left + plotW, pad.top + plotH);
      ctx.moveTo(pad.left, pad.top); ctx.lineTo(pad.left, pad.top + plotH);
      ctx.stroke();

      // Labels
      ctx.fillStyle = '#64748b';
      ctx.font = '10px Roboto, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('Log2 Fold Change', pad.left + plotW / 2, h - 10);

      ctx.save();
      ctx.translate(14, pad.top + plotH / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText('-log10 p-value', 0, 0);
      ctx.restore();

      // Draw Points
      volcanoData.forEach(d => {{
        const cx = mapX(d.log2fc);
        const cy = mapY(d.negLogP);
        d._cx = cx;
        d._cy = cy;

        ctx.beginPath();
        ctx.arc(cx, cy, 3.2, 0, Math.PI * 2);

        if (d.status === 'UP') {{
          ctx.fillStyle = 'rgba(239, 68, 68, 0.75)';
        }} else if (d.status === 'DOWN') {{
          ctx.fillStyle = 'rgba(37, 99, 235, 0.75)';
        }} else {{
          ctx.fillStyle = 'rgba(148, 163, 184, 0.35)';
        }}
        ctx.fill();
      }});

      // Draw Highlighted Gene if any
      if (highlightedVolcanoGene) {{
        const pt = volcanoData.find(d => d.gene.toUpperCase() === highlightedVolcanoGene.toUpperCase());
        if (pt) {{
          ctx.save();
          ctx.beginPath();
          ctx.arc(pt._cx, pt._cy, 7, 0, Math.PI * 2);
          ctx.strokeStyle = '#f59e0b';
          ctx.lineWidth = 3;
          ctx.stroke();

          ctx.fillStyle = '#0f172a';
          ctx.font = 'bold 12px "JetBrains Mono", monospace';
          ctx.textAlign = 'left';
          ctx.fillText(pt.gene, pt._cx + 9, pt._cy + 4);
          ctx.restore();
        }}
      }}
    }}

    function onVolcanoSearch(val) {{
      const q = val.trim().toUpperCase();
      highlightedVolcanoGene = q || null;
      renderVolcanoCanvas();

      if (q) {{
        const found = volcanoData.find(d => d.gene.toUpperCase() === q);
        if (found) setVolcanoGeneCard(found);
      }}
      renderGeneSidebar(document.getElementById('sidebar-gene-search').value);
    }}

    function setVolcanoGeneCard(d) {{
      document.getElementById('volcano-gene-name').innerText = d.gene;
      const bEl = document.getElementById('volcano-gene-badge');
      if (d.status === 'UP') {{
        bEl.className = 'badge-pill red';
        bEl.innerText = 'Enriched in Grp 1';
      }} else if (d.status === 'DOWN') {{
        bEl.className = 'badge-pill blue';
        bEl.innerText = 'Enriched in Grp 2';
      }} else {{
        bEl.className = 'badge-pill amber';
        bEl.innerText = 'Not Significant';
      }}

      document.getElementById('v-det-fc').innerText = `${{d.log2fc >= 0 ? '+' : ''}}${{d.log2fc.toFixed(2)}}`;
      document.getElementById('v-det-pval').innerText = d.p < 0.001 ? d.p.toExponential(2) : d.p.toFixed(4);
      document.getElementById('v-det-m1').innerText = `${{d.m1.toFixed(1)}} TPM`;
      document.getElementById('v-det-m2').innerText = `${{d.m2.toFixed(1)}} TPM`;
    }}

    function inspectSelectedGeneTrajectory() {{
      const gene = document.getElementById('volcano-gene-name').innerText;
      if (gene && GENES[gene]) {{
        setTrajectoryGene(gene);
        switchView('tab-genes');
      }}
    }}

    function renderVolcanoTopTable() {{
      const tbody = document.getElementById('volcano-top-tbody');
      if (!tbody) return;

      const sigs = volcanoData.filter(d => d.status !== 'NS');
      sigs.sort((a,b) => Math.abs(b.log2fc * b.negLogP) - Math.abs(a.log2fc * a.negLogP));

      tbody.innerHTML = sigs.slice(0, 15).map(d => `
        <tr style="cursor: pointer;" onclick="onVolcanoSearch('${{d.gene}}')">
          <td class="mono"><b>${{d.gene}}</b></td>
          <td class="num" style="color: ${{d.log2fc >= 0 ? 'var(--up-red)' : 'var(--down-blue)'}};">${{d.log2fc >= 0 ? '+' : ''}}${{d.log2fc.toFixed(2)}}</td>
          <td class="num">${{d.p < 0.001 ? d.p.toExponential(1) : d.p.toFixed(3)}}</td>
        </tr>
      `).join('');
    }}

    // Canvas Mouse Interaction
    const vCanvas = document.getElementById('volcano-canvas');
    const vTooltip = document.getElementById('volcano-tooltip');

    if (vCanvas && vTooltip) {{
      vCanvas.addEventListener('mousemove', e => {{
        const rect = vCanvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        let closest = null;
        let minD = 10;

        for (const d of volcanoData) {{
          if (d._cx === undefined) continue;
          const dist = Math.hypot(d._cx - mx, d._cy - my);
          if (dist < minD) {{
            minD = dist;
            closest = d;
          }}
        }}

        if (closest) {{
          vTooltip.style.display = 'block';
          vTooltip.style.left = `${{mx + 12}}px`;
          vTooltip.style.top = `${{my - 24}}px`;
          vTooltip.innerHTML = `
            <div style="font-weight: 800; font-family: 'JetBrains Mono', monospace; font-size: 0.86rem; color: #38bdf8;">${{closest.gene}}</div>
            <div>Log2FC: <b>${{closest.log2fc >= 0 ? '+' : ''}}${{closest.log2fc.toFixed(2)}}</b></div>
            <div>p-value: <b>${{closest.p < 0.001 ? closest.p.toExponential(2) : closest.p.toFixed(4)}}</b></div>
            <div style="font-size: 0.72rem; color: #94a3b8; margin-top: 2px;">Grp 1: ${{closest.m1.toFixed(1)}} TPM | Grp 2: ${{closest.m2.toFixed(1)}} TPM</div>
          `;
        }} else {{
          vTooltip.style.display = 'none';
        }}
      }});

      vCanvas.addEventListener('mouseleave', () => {{
        vTooltip.style.display = 'none';
      }});

      vCanvas.addEventListener('click', e => {{
        const rect = vCanvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        let closest = null;
        let minD = 12;

        for (const d of volcanoData) {{
          if (d._cx === undefined) continue;
          const dist = Math.hypot(d._cx - mx, d._cy - my);
          if (dist < minD) {{
            minD = dist;
            closest = d;
          }}
        }}

        if (closest) {{
          setVolcanoGeneCard(closest);
          onVolcanoSearch(closest.gene);
        }}
      }});
    }}

    // CSV Exports
    function exportVolcanoCSV() {{
      const headers = ['gene', 'log2fc_log2tpm1', 'p_value', 'bh_q_value', 'neg_log10_p', 'group1_mean_tpm', 'group2_mean_tpm', 'significance_status'];
      const rows = volcanoData.map(d => [
        `"${{d.gene}}"`,
        d.log2fc.toFixed(4),
        d.p.toExponential(4),
        d.q.toExponential(4),
        d.negLogP.toFixed(4),
        d.m1.toFixed(2),
        d.m2.toFixed(2),
        `"${{d.status}}"`
      ].join(','));
      downloadCSV('mmr_cohort_volcano_differential.csv', [headers.join(','), ...rows].join('\\n'));
    }}

    function exportClonotypesCSV() {{
      const headers = ['sample_id', 'patient_id', 'receptor', 'locus', 'v_gene', 'd_gene', 'j_gene', 'cdr3_aa', 'cdr3_len', 'reads', 'in_frame', 'shared_with'];
      const rows = CLONOTYPES.map(c => headers.map(h => `"${{Array.isArray(c[h]) ? c[h].join(';') : (c[h] ?? '')}}"`).join(','));
      downloadCSV('mmr_cohort_cdr3_clonotypes.csv', [headers.join(','), ...rows].join('\\n'));
    }}

    function exportFusionsCSV() {{
      const headers = ['sample_id', 'path_dx', 'gene_up', 'gene_down', 'chr_up', 'pos_up', 'chr_dn', 'pos_dn', 'split_reads', 'spanning_reads', 'fusion_class', 'known_type', 'sv_type', 'exon_junction', 'junc_up', 'junc_dn'];
      const rows = FUSIONS.map(f => headers.map(h => `"${{f[h] || ''}}"`).join(','));
      downloadCSV('mmr_cohort_rna_fusions.csv', [headers.join(','), ...rows].join('\\n'));
    }}

    function exportCohortLandscapeCSV() {{
      const headers = ['tissue_id', 'patient_id', 'tissue_origin', 'path_dx', 'timepoint', 'tissue_format', 'tcr_reads', 'bcr_reads', 'simpson_clonality', 'cytolytic_score', 'checkpoint_score', 'cta_score'];
      const rows = COHORT.map(s => headers.map(h => `"${{s[h] || ''}}"`).join(','));
      downloadCSV('mmr_cohort_landscape.csv', [headers.join(','), ...rows].join('\\n'));
    }}

    function downloadCSV(filename, content) {{
      const blob = new Blob([content], {{ type: 'text/csv;charset=utf-8;' }});
      const link = document.createElement('a');
      link.href = URL.createObjectURL(blob);
      link.download = filename;
      link.click();
    }}
  </script>
</body>
</html>
"""
    return html


# ---------------------------------------------------------------------------
# Main Execution Flow
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    results_dir = Path(args.results_dir)
    metadata_tsv = Path(args.metadata_tsv)
    output_html = args.output_html or str(results_dir / "cohort_mmr_rna_report.html")

    print("======================================================================")
    print("MMR Cohort RNA-Seq & TCR Interactive Clinical Dashboard Compiler")
    print(f"  Results Dir:   {results_dir}")
    print(f"  Metadata TSV:  {metadata_tsv}")
    print(f"  Output HTML:   {output_html}")
    print("======================================================================")

    # 1. Parse clinical metadata
    by_sample, metadata_records = parse_metadata(metadata_tsv)
    print(f"Loaded {len(metadata_records)} metadata records ({len(by_sample)} lookup keys).")

    # 2. Locate oncoanalyser output directories
    sample_dirs = locate_sample_dirs(results_dir)
    print(f"Discovered {len(sample_dirs)} sample directories under {results_dir}.")

    # Dashboard cohort: samples with RNA-seq results only (metadata rows without RNA-seq are omitted)
    all_sample_ids = sorted(s_id for s_id, s_dir in sample_dirs.items()
                            if resolve_tool_dir(s_dir, "isofox").is_dir())
    skipped_meta = sum(1 for r in metadata_records if r["tissue_id"] and r["tissue_id"] not in all_sample_ids)
    print(f"Including {len(all_sample_ids)} samples with RNA-seq results; omitting {skipped_meta} metadata rows without RNA-seq.")

    cohort_data = []
    sample_gene_maps = []
    sample_ref_maps = []
    all_fusions = []
    all_clonotypes = []

    for idx, s_id in enumerate(all_sample_ids):
        meta = by_sample.get(s_id) or by_sample.get(s_id.replace("_", "-")) or {}
        s_dir = sample_dirs.get(s_id) or sample_dirs.get(s_id.replace("-", "_"))

        # Parsers
        isofox_summary = parse_isofox_summary(s_dir) if s_dir else {}
        star_qc = parse_star_log(s_dir) if s_dir else {}
        gene_tpm, gene_ref = parse_isofox_expression(s_dir) if s_dir else ({}, {})
        fusions = parse_isofox_fusions(s_dir) if s_dir else []
        locus_stats = parse_cider_locus_stats(s_dir) if s_dir else {}
        clonotypes, tcr_stats = parse_cider_vdj(s_dir) if s_dir else ([], EMPTY_TCR_STATS)
        cuppa = parse_cuppa(s_dir) if s_dir else {}

        # Signature scores
        sig_scores = calculate_signature_scores(gene_tpm)
        hmf_type, hmf_matched_by = hartwig_cancer_type(meta.get("path_dx", ""), meta.get("tissue_origin", ""))

        # Build sample record
        sample_rec = {
            "tissue_id": s_id,
            "patient_id": meta.get("patient_id", s_id),
            "rna_id": meta.get("rna_id", s_id),
            "tcr_id": meta.get("tcr_id", ""),
            "path_dx": meta.get("path_dx", "Unknown"),
            "tissue_origin": meta.get("tissue_origin", "Unknown"),
            "tissue_format": meta.get("tissue_format", "Unknown"),
            "library_type": meta.get("library_type") or "Unknown",
            "hartwig_cancer_type": hmf_type,
            "hartwig_matched_by": hmf_matched_by,
            "path_type": meta.get("path_type", "Tumor"),
            "date_collected": meta.get("date_collected", ""),
            "date_iso": meta.get("date_iso", ""),
            "timepoint": meta.get("timepoint", "Primary"),
            "hla_a": meta.get("hla_a", []),
            "hla_b": meta.get("hla_b", []),
            "hla_c": meta.get("hla_c", []),
            "isofox_summary": isofox_summary,
            "star_qc": star_qc,
            "tcr_reads": locus_stats.get("tcr_reads", 0),
            "bcr_reads": locus_stats.get("bcr_reads", 0),
            "shannon_diversity": tcr_stats["shannon"],
            "simpson_clonality": tcr_stats["clonality"],
            "tcr_clonotypes": tcr_stats["tcr_clonotypes"],
            "top_clone_share": tcr_stats["top_clone_share"],
            "cuppa": cuppa,
            "cytolytic_score": sig_scores["cytolytic_score"],
            "checkpoint_score": sig_scores["checkpoint_score"],
            "ap_score": sig_scores["ap_score"],
            "cta_score": sig_scores["cta_score"],
            "gbm_score": sig_scores["gbm_score"],
            "gi_score": sig_scores["gi_score"],
            "fusion_count": len(fusions),
            "has_expression": bool(gene_tpm),
            **{f"_{k}": v for k, v in compute_expression_qc_metrics(gene_tpm).items()},
        }
        cohort_data.append(sample_rec)
        sample_gene_maps.append(gene_tpm)
        sample_ref_maps.append(gene_ref)

        # Attach fusions
        for f in fusions:
            all_fusions.append({
                "sample_id": s_id,
                "patient_id": meta.get("patient_id", s_id),
                "path_dx": meta.get("path_dx", "Unknown"),
                **f,
            })

        # Attach clonotypes
        for c in clonotypes:
            all_clonotypes.append({
                "sample_id": s_id,
                "patient_id": meta.get("patient_id", s_id),
                **c,
            })

    # TCR clonotypes (same locus + CDR3) seen in more than one sample from the same patient
    clone_samples = {}
    for c in all_clonotypes:
        if c["receptor"] == "TCR":
            clone_samples.setdefault((c["patient_id"], c["locus"], c["cdr3_aa"]), set()).add(c["sample_id"])
    for c in all_clonotypes:
        others = clone_samples.get((c["patient_id"], c["locus"], c["cdr3_aa"]), set()) - {c["sample_id"]}
        c["shared_with"] = sorted(others)
    for s in cohort_data:
        s["shared_tcr_clones"] = sum(1 for c in all_clonotypes if c["sample_id"] == s["tissue_id"] and c["shared_with"])

    # Fusion neoepitopes cross-checked against patient HLA alleles and TCR clonotypes
    receptor_matches = predict_fusion_neoantigen_tcr_pairs(all_fusions, all_clonotypes, cohort_data, getattr(args, "binding_dir", None))
    print(f"Fusion classes: " + ", ".join(f"{k} {v}" for k, v in sorted(
        {c: sum(1 for f in all_fusions if f["fusion_class"] == c) for c in {f["fusion_class"] for f in all_fusions}}.items())))
    print(f"Fusion neoantigen pairs: {len(receptor_matches)}, with presenting HLA binder (rank <= 2%): "
          f"{sum(1 for r in receptor_matches if r['hla_binder'])}")

    # Embed the top clonotypes per sample and receptor type
    kept, per = [], {}
    for c in all_clonotypes:  # already sorted by reads within each sample
        k = (c["sample_id"], c["receptor"])
        per[k] = per.get(k, 0) + 1
        if per[k] <= args.max_cdr3_per_sample:
            kept.append(c)
    all_clonotypes = kept

    # Cancer-type-matched Hartwig percentiles (equivalent to Isofox PercentileCancer / MedianTpmCancer)
    needed_types = sorted({s["hartwig_cancer_type"] for s in cohort_data})
    print(f"Loading Hartwig gene distribution for: {', '.join(needed_types)}...")
    gene_dist = load_gene_distribution(args.gene_distribution, needed_types)
    for s, gene_tpm, gene_ref in zip(cohort_data, sample_gene_maps, sample_ref_maps):
        ref = gene_dist.get(s["hartwig_cancer_type"])
        if not ref:
            s["hartwig_cancer_type"] = ""
            continue
        print(f"  {s['tissue_id']}: {s['path_dx']} -> {s['hartwig_cancer_type']} (by {s['hartwig_matched_by']})")
        for g, (pan_pct, pan_med, can_pct, can_med) in gene_ref.items():
            if can_pct < 0 and g in ref:
                med, pcts = ref[g]
                gene_ref[g] = (pan_pct, pan_med, distribution_percentile(gene_tpm.get(g, 0.0), pcts), med)
    del gene_dist

    # 3. Sample QC (fails are excluded from cohort-level expression statistics)
    strand_check = evaluate_sample_qc(cohort_data)
    qc_pass = [s["qc"]["status"] in ("PASS", "WARN") for s in cohort_data]
    qc_counts = {}
    for s in cohort_data:
        qc_counts[s["qc"]["status"]] = qc_counts.get(s["qc"]["status"], 0) + 1
    print(f"Sample QC: {qc_counts}" + ("" if strand_check else " (strand check disabled: unstranded cohort)"))
    for s in cohort_data:
        if s["qc"]["status"] == "FAIL":
            print(f"  QC FAIL {s['tissue_id']}: {', '.join(s['qc']['fails'])}")

    # 4. Build High-Information Gene Matrix
    priority_genes = [
        # Cancer-Testis Antigens (CTAs)
        "CTAG1B", "CTAG2", "MAGEA1", "MAGEA3", "MAGEA4", "MAGEA10", "MAGEC2", "PRAME", "SSX2", "SSX4",
        # Glioblastoma & Neural Markers
        "EGFR", "SOX2", "OLIG2", "GFAP", "IL13RA2", "PROM1", "ERBB2", "NES", "CD44", "CHI3L1", "PDGFRA", "MET",
        # Gastrointestinal & Epithelial Markers
        "CEACAM5", "EPCAM", "MUC1", "MSLN", "CDX2", "KRT20", "KRT7", "TFF1",
        # Immune Checkpoints & Infiltration
        "CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD8A", "CD8B", "GZMA", "PRF1", "IFNG", "FOXP3", "CXCL9", "CXCL10", "CD4", "CD3E",
        # HLA & Antigen Presentation
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "B2M", "TAP1", "TAP2", "NLRC5", "CIITA",
        # Oncogenic Drivers & Suppressors
        "TP53", "PTEN", "IDH1", "IDH2", "ATRX", "TERT", "RB1", "NF1", "PIK3CA", "BRAF", "KRAS", "MYC", "CDK4", "CDKN2A"
    ]

    # Calculate variance across all genes for top variable gene selection
    all_expressed_genes = set()
    for gmap in sample_gene_maps:
        all_expressed_genes.update(gmap.keys())

    # Variance of log2(TPM+1) over samples with expression data (raw-TPM variance favours highly expressed genes)
    expressed_maps = [gmap for gmap, ok in zip(sample_gene_maps, qc_pass) if gmap and ok]
    gene_stats = []
    n_samples = max(1, len(expressed_maps))
    for g in all_expressed_genes:
        vals = [gmap.get(g, 0.0) for gmap in expressed_maps]
        if max(vals) < 0.5:
            continue
        logs = [math.log2(v + 1.0) for v in vals]
        mean_v = sum(logs) / n_samples
        var_v = sum((v - mean_v) ** 2 for v in logs) / n_samples
        gene_stats.append((g, var_v))

    gene_stats.sort(key=lambda x: x[1], reverse=True)
    top_variable_genes = [g for g, _ in gene_stats[:3000]]

    # Include priority genes, top variable genes, and all fusion partner genes
    target_genes = set(priority_genes) | set(top_variable_genes)
    for f in all_fusions:
        if f.get("gene_up"): target_genes.add(f["gene_up"])
        if f.get("gene_down"): target_genes.add(f["gene_down"])

    print(f"Building high-information expression matrix across {len(target_genes)} genes...")
    gene_matrix = {}
    for g in sorted(target_genes):
        gene_matrix[g] = [round(gmap.get(g, 0.0), 2) for gmap in sample_gene_maps]
    # Outlier genes appended below are selected per sample, so they are excluded from cohort-level tests
    volcano_genes = sorted(target_genes)

    # Hartwig pan-cancer percentile (Isofox PercentileCohort, 0-1) for curated genes; None when unavailable
    def _pan_pct(rmap, g):
        v = rmap.get(g, (-1.0,))[0]
        if v < 0:
            return None
        return round(v / 100.0 if v > 1.0 else v, 4)
    gene_pan_pct = {g: [_pan_pct(rmap, g) for rmap in sample_ref_maps] for g in priority_genes if g in gene_matrix}

    def _can_pct(rmap, g):
        v = rmap.get(g, (-1.0, -1.0, -1.0))[2]
        return None if v < 0 else round(v, 4)
    gene_can_pct = {g: [_can_pct(rmap, g) for rmap in sample_ref_maps] for g in priority_genes if g in gene_matrix}

    # 5. Per-patient outlier expression with antigen annotation
    antigens = load_antigen_annotation(args.antigen_tsv)
    coding_genes = load_protein_coding_genes(args.antigen_tsv)
    print(f"Loaded {len(antigens)} annotated antigen genes and {len(coding_genes)} protein-coding symbols.")
    sample_groups = [(s["tissue_origin"], s["library_type"]) for s in cohort_data]
    outlier_lists = compute_patient_outliers(
        sample_gene_maps, sample_ref_maps, sample_groups, qc_pass, antigens, coding_genes,
        max_outliers=args.max_outliers_per_sample
    )
    patient_outliers = {s["tissue_id"]: rows for s, rows in zip(cohort_data, outlier_lists)}
    used_antigens = {r[0] for rows in outlier_lists for r in rows if r[0] in antigens}
    antigen_db = {g: antigens[g] for g in sorted(used_antigens)}
    print(f"Computed outlier profiles: {sum(sum(r[7] for r in rows) for rows in outlier_lists)} outlier calls, "
          f"{len(antigen_db)} distinct annotated antigens.")

    # Embed outlier genes in the expression matrix so they can be opened in the gene view
    for rows in outlier_lists:
        for r in rows:
            if r[0] not in gene_matrix and not r[8]:
                gene_matrix[r[0]] = [round(gmap.get(r[0], 0.0), 2) for gmap in sample_gene_maps]

    # 6. Generate HTML Report
    print("Assembling standalone interactive clinical HTML dashboard...")
    html_content = build_dashboard_html(cohort_data, gene_matrix, all_fusions, all_clonotypes,
                                        patient_outliers, antigen_db, volcano_genes, gene_pan_pct, gene_can_pct,
                                        receptor_matches)

    output_path = Path(output_html)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"Dashboard generated successfully: {output_path}")
    print(f"File size: {output_path.stat().st_size / (1024*1024):.2f} MB")
    print("======================================================================")


if __name__ == "__main__":
    main()
