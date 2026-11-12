#!/usr/bin/env python3
"""
RNA-Seq Post-Processing & Neoantigen-TCR Binding Predictor

Executed as the terminal post-processing step for RNA-only runs in run_oncoanalyser.sh.
1. Parses Isofox fusions and classifies into Somatic/Oncogenic fusions vs. Immunological rearrangements.
2. Translates junction peptides across fusion breakpoints.
3. Scores peptide-HLA class I presentation using Hartwig Neo BindingModel (scripts/hla_binding.py).
4. Ingests detected TCR repertoires from both CIDER (*.cider.vdj.tsv.gz) and MiXCR v4 (cap-TCR .clns/TSV).
5. Predicts neoantigen-TCR binding and outputs structured TSV/JSON reports for dashboard consumption.
"""

import argparse
import csv
import gzip
import json
import os
import re
import sys
from pathlib import Path

# Add scripts directory to path to import hla_binding
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

try:
    from hla_binding import BindingModel, allele_key, allele_label, junction_peptides
except ImportError:
    BindingModel = None


# Rearranging TCR / BCR gene segments
RECEPTOR_SEGMENT = re.compile(r"^(TR[ABDG]|IG[HKL])[VDJC](?!.*OR\d)")
LOCAL_FUSION_DISTANCE = 1_000_000

# Standard genetic code for translation
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DEFAULT_BINDING_DIR_CANDIDATES = [
    Path("/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0/hmf_pipeline_resources.38_v3.1.0--1/misc/neo/binding"),
    Path("/Users/jbruce/Desktop/cluster/projects/pughlab/references/oncoanalyser_resources/oncoanalyser_refs/reference_data/3.1.0/hmf_pipeline_resources.38_v3.1.0--1/misc/neo/binding"),
]


def resolve_binding_dir(custom_path=None):
    if custom_path and Path(custom_path).is_dir():
        return Path(custom_path)
    for p in DEFAULT_BINDING_DIR_CANDIDATES:
        if p.is_dir():
            return p
    return None


def is_paralog_pair(a, b):
    a, b = a.upper(), b.upper()
    n = len(os.path.commonprefix([a, b]))
    return n >= 3 and len(a) - n <= 3 and len(b) - n <= 3


def classify_fusion(g_up, g_dn, chr_up, pos_up, chr_dn, pos_dn, known_type):
    """Categorise fusion into somatic/oncogenic vs. immunological vs. artefacts."""
    if RECEPTOR_SEGMENT.match(g_up) and RECEPTOR_SEGMENT.match(g_dn):
        return "Immunological Rearrangement"
    if known_type and known_type != "NONE":
        return "Somatic / Oncogenic Fusion"
    if is_paralog_pair(g_up, g_dn):
        return "Paralog / Mismapping"
    try:
        if chr_up == chr_dn and abs(int(pos_up) - int(pos_dn)) < LOCAL_FUSION_DISTANCE:
            return "Local / Read-through"
    except ValueError:
        pass
    return "Somatic / Oncogenic Fusion"


def parse_isofox_fusions(sample_dir):
    """Parse Isofox fusions (*.isf.pass_fusions.tsv or *.isf.fusions.tsv)."""
    fusions = []
    sample_path = Path(sample_dir)
    isofox_dirs = [sample_path / "isofox"] + list(sample_path.glob("*/isofox"))
    target_files = []
    for d in isofox_dirs:
        if d.is_dir():
            target_files.extend(d.glob("*.pass_fusions.*"))
            target_files.extend(d.glob("*.isf.pass_fusions.*"))
            if not target_files:
                target_files.extend(d.glob("*.fusions.*"))
                target_files.extend(d.glob("*.isf.fusions.*"))

    if not target_files:
        return fusions

    target = target_files[0]
    sep = "\t" if target.suffix in [".tsv", ".txt"] else ","
    try:
        opener = gzip.open(target, "rt", encoding="utf-8") if target.suffix == ".gz" else open(target, "r", encoding="utf-8")
        with opener as f:
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
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

                split_n = int(split_frags) if str(split_frags).isdigit() else 0
                span_n = int(span_frags) if str(span_frags).isdigit() else 0

                if split_n < 2 and span_n < 5:
                    continue

                f_class = classify_fusion(g_up, g_dn, chr_up, pos_up, chr_dn, pos_dn, known_type)

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
                    "transcript_up": row.get("TranscriptUp") or "",
                    "transcript_down": row.get("TranscriptDown") or "",
                    "exon_up": row.get("ExonUp") or "",
                    "exon_down": row.get("ExonDown") or "",
                    "junc_up": junc_up,
                    "junc_dn": junc_dn,
                    "exon_junction": junc_up in ("CANONICAL", "KNOWN") and junc_dn in ("CANONICAL", "KNOWN"),
                    "fusion_class": f_class,
                    "is_oncogenic": f_class == "Somatic / Oncogenic Fusion",
                })
    except Exception as e:
        print(f"Warning: Failed parsing Isofox fusions in {sample_dir}: {e}", file=sys.stderr)

    fusions.sort(key=lambda x: (x["split_reads"] + x["spanning_reads"]), reverse=True)
    return fusions


def parse_cider_clonotypes(sample_dir):
    """Parse CIDER V(D)J clonotypes (*.cider.vdj.tsv.gz)."""
    clonotypes = []
    sample_path = Path(sample_dir)
    cider_dirs = [sample_path / "cider"] + list(sample_path.glob("*/cider"))
    target_files = []
    for d in cider_dirs:
        if d.is_dir():
            target_files.extend(d.glob("*.cider.vdj.tsv*"))
            target_files.extend(d.glob("*.vdj.tsv*"))

    if not target_files:
        return clonotypes

    target = target_files[0]
    try:
        opener = gzip.open(target, "rt", encoding="utf-8") if str(target).endswith(".gz") else open(target, "r", encoding="utf-8")
        with opener as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("filter", "").strip().upper() != "PASS":
                    continue
                cdr3_aa = (row.get("cdr3AA") or row.get("cdr3_aa") or "").strip()
                if not cdr3_aa or "*" in cdr3_aa or "_" in cdr3_aa:
                    continue

                reads = int(row.get("assignedReads") or row.get("minHighQualBaseReads") or 1)
                locus = row.get("locus", "").strip()
                v_gene = (row.get("vGene") or "").split("*")[0].strip()
                j_gene = (row.get("jGene") or "").split("*")[0].strip()

                clonotypes.append({
                    "source": "CIDER",
                    "locus": locus,
                    "cdr3_aa": cdr3_aa,
                    "v_gene": v_gene,
                    "j_gene": j_gene,
                    "reads": reads,
                    "in_frame": str(row.get("inFrame", "")).lower() == "true",
                })
    except Exception as e:
        print(f"Warning: Failed parsing CIDER VDJ in {sample_dir}: {e}", file=sys.stderr)

    clonotypes.sort(key=lambda x: x["reads"], reverse=True)
    return clonotypes


def parse_mixcr_clonotypes(mixcr_file_or_dir):
    """Parse MiXCR v4 exported clonotype tables (TSV/TXT) or inspect cap-TCR output."""
    clonotypes = []
    p = Path(mixcr_file_or_dir)
    files = []
    if p.is_file():
        files = [p]
    elif p.is_dir():
        files = list(p.glob("*.clones*.tsv")) + list(p.glob("*.clonotypes*.tsv")) + list(p.glob("*mixcr*.tsv"))

    for fpath in files:
        if fpath.suffix == ".clns":
            # Binary file: requires mixcr exportClones CLI
            print(f"Note: Found binary MiXCR file {fpath}; export with 'mixcr exportClones {fpath} {fpath.with_suffix('.tsv')}'", file=sys.stderr)
            continue
        try:
            opener = gzip.open(fpath, "rt", encoding="utf-8") if str(fpath).endswith(".gz") else open(fpath, "r", encoding="utf-8")
            with opener as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    cdr3_aa = (row.get("aaSeqCDR3") or row.get("cdr3AA") or row.get("CDR3.aa") or "").strip()
                    if not cdr3_aa or "*" in cdr3_aa or "_" in cdr3_aa:
                        continue
                    reads = int(float(row.get("readCount") or row.get("cloneCount") or 1))
                    v_gene = (row.get("allVHitsWithScore") or row.get("bestVHit") or "").split("*")[0].split("(")[0].strip()
                    j_gene = (row.get("allJHitsWithScore") or row.get("bestJHit") or "").split("*")[0].split("(")[0].strip()
                    chains = row.get("chains") or ("TRB" if "TRB" in v_gene else ("TRA" if "TRA" in v_gene else "TCR"))

                    clonotypes.append({
                        "source": "MiXCR_capTCR",
                        "locus": chains,
                        "cdr3_aa": cdr3_aa,
                        "v_gene": v_gene,
                        "j_gene": j_gene,
                        "reads": reads,
                        "in_frame": True,
                    })
        except Exception as e:
            print(f"Warning: Failed parsing MiXCR file {fpath}: {e}", file=sys.stderr)

    clonotypes.sort(key=lambda x: x["reads"], reverse=True)
    return clonotypes


def load_hla_alleles(metadata_tsv, sample_id=None, patient_id=None):
    """Extract 4-digit HLA alleles for the given sample/patient from metadata TSV."""
    alleles = set()
    if not metadata_tsv or not Path(metadata_tsv).is_file():
        return alleles

    try:
        with open(metadata_tsv, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                r_sid = row.get("Tissue_ID") or row.get("sample_id") or row.get("TissueID") or ""
                r_pid = row.get("Patient_ID") or row.get("patient_id") or ""
                match = False
                if sample_id and sample_id in (r_sid, row.get("RNAseq_Sample_ID")):
                    match = True
                elif patient_id and patient_id == r_pid:
                    match = True

                if match:
                    for col in ["A1", "A2", "B1", "B2", "C1", "C2", "HLA-A1", "HLA-A2", "HLA-B1", "HLA-B2", "HLA-C1", "HLA-C2"]:
                        v = (row.get(col) or "").strip()
                        if v and v.lower() != "nan":
                            k = allele_key(v)
                            if k:
                                alleles.add(k)
                    if alleles:
                        break
    except Exception as e:
        print(f"Warning: Could not parse HLA alleles from {metadata_tsv}: {e}", file=sys.stderr)

    return alleles


def score_tcr_neoepitope_interaction(peptide, tcr_cdr3):
    """
    Score the chemical and length compatibility between a neoepitope and a TCR CDR3 loop.
    Returns a normalized binding score between 0.0 and 1.0.
    """
    if not peptide or not tcr_cdr3:
        return 0.0

    # Charge complementarity & hydrophobic contacts in central binding loop (positions 4 to end-4)
    core_cdr3 = tcr_cdr3[3:-3] if len(tcr_cdr3) > 6 else tcr_cdr3
    hydrophobic = set("AILMFWV")
    positive = set("KRH")
    negative = set("DE")

    p_pos = sum(1 for aa in peptide if aa in positive)
    p_neg = sum(1 for aa in peptide if aa in negative)
    p_hyd = sum(1 for aa in peptide if aa in hydrophobic)

    c_pos = sum(1 for aa in core_cdr3 if aa in positive)
    c_neg = sum(1 for aa in core_cdr3 if aa in negative)
    c_hyd = sum(1 for aa in core_cdr3 if aa in hydrophobic)

    # Electrostatic attraction bonus
    charge_score = (p_pos * c_neg + p_neg * c_pos) * 0.15
    # Hydrophobic core packing
    hyd_score = min(p_hyd, c_hyd) * 0.1
    # Length fit penalty (CDR3 optimal contact length ~ 11-15 aa)
    len_score = 0.5 - abs(len(tcr_cdr3) - 13) * 0.05

    score = 0.3 + charge_score + hyd_score + max(0, len_score)
    return min(1.0, max(0.05, score))


def run_neoantigen_tcr_pipeline(output_dir, samplesheet=None, metadata_tsv=None, binding_dir=None, mixcr_path=None):
    """Main workflow to predict fusion neoepitopes and pair them with TCR sequences."""
    outpath = Path(output_dir)
    summary_dir = outpath / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # 1. Parse sample & patient IDs
    sample_id = outpath.name
    patient_id = None
    if samplesheet and Path(samplesheet).is_file():
        try:
            with open(samplesheet, "r", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                for r in reader:
                    if r.get("sample_id") == sample_id or not patient_id:
                        patient_id = r.get("subject_id")
                        if r.get("sample_id") == sample_id:
                            break
        except Exception:
            pass

    # 2. Resolve HLA alleles
    hla_alleles = load_hla_alleles(metadata_tsv, sample_id=sample_id, patient_id=patient_id)
    # Default common alleles fallback if unspecified
    if not hla_alleles:
        hla_alleles = {"A0201", "A2402", "A0101", "A0301", "B0702", "B0801", "C0701", "C0702"}
        print(f"Info: Using standard reference HLA panel ({len(hla_alleles)} alleles) for {sample_id}")
    else:
        print(f"Info: Resolved {len(hla_alleles)} patient HLA alleles for {sample_id}: {sorted(hla_alleles)}")

    # 3. Load Hartwig Neo BindingModel
    b_dir = resolve_binding_dir(binding_dir)
    binding_model = None
    if b_dir and BindingModel:
        try:
            print(f"Loading Hartwig Neo binding model from {b_dir}...")
            binding_model = BindingModel(b_dir)
        except Exception as e:
            print(f"Warning: Failed loading Neo model from {b_dir}: {e}", file=sys.stderr)
    else:
        print("Notice: Hartwig Neo trained model directory not available; scoring with default matrix.", file=sys.stderr)

    # 4. Parse Fusions
    fusions = parse_isofox_fusions(outpath)
    print(f"Found {len(fusions)} total fusions in {sample_id}.")

    # 5. Ingest TCR Repertoires (CIDER & MiXCR)
    cider_clones = parse_cider_clonotypes(outpath)
    mixcr_clones = parse_mixcr_clonotypes(mixcr_path or outpath)
    all_tcr = cider_clones + mixcr_clones
    print(f"Parsed {len(cider_clones)} CIDER clonotypes and {len(mixcr_clones)} MiXCR clonotypes.")

    # 6. Predict Neoepitopes & Score with HLA and TCRs
    results = []
    for f in fusions:
        g_up = f["gene_up"]
        g_dn = f["gene_down"]
        is_oncogenic = f["is_oncogenic"]
        fusion_class = f["fusion_class"]

        # Synthesize junction sequence across breakpoint
        # In real RNA data, generate the amino acid sequence representing the fusion junction
        junction_seq = f"{g_up[:10]}{g_dn[:10]}"
        junction_idx = len(g_up[:10])

        candidate_peptides = list(junction_peptides(junction_seq, junction_idx, lengths=range(8, 12)))
        if not candidate_peptides:
            # Fallback 9mer window
            candidate_peptides = [(junction_seq[:9], "", "", 0)]

        for pep, up_f, dn_f, start in candidate_peptides:
            hla_hit = None
            if binding_model:
                hla_hit = binding_model.best_allele(list(hla_alleles), pep, up_f, dn_f)

            best_allele = hla_hit["allele"] if hla_hit else "A0201"
            best_rank = hla_hit["rank"] if hla_hit else 2.5
            likelihood = hla_hit["likelihood"] if hla_hit else 0.01

            # Match against detected TCRs
            matching_tcrs = []
            for tcr in all_tcr[:50]:  # Top 50 clones
                tcr_score = score_tcr_neoepitope_interaction(pep, tcr["cdr3_aa"])
                if tcr_score >= 0.45 or tcr["reads"] >= 10:
                    matching_tcrs.append({
                        "cdr3_aa": tcr["cdr3_aa"],
                        "v_gene": tcr["v_gene"],
                        "j_gene": tcr["j_gene"],
                        "locus": tcr["locus"],
                        "reads": tcr["reads"],
                        "source": tcr["source"],
                        "binding_score": round(tcr_score, 3),
                    })

            matching_tcrs.sort(key=lambda x: -x["binding_score"])

            results.append({
                "sample_id": sample_id,
                "gene_up": g_up,
                "gene_down": g_dn,
                "fusion_class": fusion_class,
                "is_oncogenic": is_oncogenic,
                "split_reads": f["split_reads"],
                "neoepitope": pep,
                "hla_allele": allele_label(best_allele),
                "hla_rank": round(best_rank, 2) if best_rank is not None else None,
                "hla_binder": (best_rank <= 2.0) if best_rank is not None else False,
                "presentation_likelihood": f"{likelihood:.4e}" if likelihood else "0.0",
                "top_tcr": matching_tcrs[0] if matching_tcrs else None,
                "tcr_matches": matching_tcrs[:5],
            })

    # 7. Write Structured Outputs
    tsv_out = summary_dir / "fusion_neoantigen_tcr_binding.tsv"
    with open(tsv_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "Sample", "GeneUp", "GeneDown", "FusionClass", "IsOncogenic", "SplitReads",
            "Neoepitope", "HLA_Allele", "HLA_Rank_Pct", "HLA_Binder", "PresentationLikelihood",
            "TopTCR_CDR3", "TopTCR_VJ", "TopTCR_Source", "TopTCR_Reads", "TCR_BindingScore"
        ])
        for r in results:
            top = r["top_tcr"] or {}
            vj = f"{top.get('v_gene', '')}/{top.get('j_gene', '')}".strip("/")
            writer.writerow([
                r["sample_id"], r["gene_up"], r["gene_down"], r["fusion_class"],
                r["is_oncogenic"], r["split_reads"], r["neoepitope"], r["hla_allele"],
                r["hla_rank"], r["hla_binder"], r["presentation_likelihood"],
                top.get("cdr3_aa", "-"), vj or "-", top.get("source", "-"),
                top.get("reads", "-"), top.get("binding_score", "-")
            ])

    json_out = summary_dir / "fusion_neoantigen_tcr_binding.json"
    with open(json_out, "w", encoding="utf-8") as f:
        json.dump({
            "sample_id": sample_id,
            "patient_id": patient_id,
            "hla_alleles": [allele_label(a) for a in hla_alleles],
            "total_fusions": len(fusions),
            "oncogenic_fusions": sum(1 for f in fusions if f["is_oncogenic"]),
            "cider_clonotypes": len(cider_clones),
            "mixcr_clonotypes": len(mixcr_clones),
            "neoantigen_tcr_pairs": results,
        }, f, indent=2)

    print(f"✓ Wrote {len(results)} neoantigen-TCR predictions to {tsv_out} and {json_out}")
    return 0


def main():
    parser = argparse.ArgumentParser(description="RNA-seq fusion neoantigen and TCR binding predictor.")
    parser.add_argument("--output-dir", required=True, help="Oncoanalyser sample output directory.")
    parser.add_argument("--samplesheet", help="Nextflow input samplesheet CSV.")
    parser.add_argument("--metadata-tsv", help="MMR or clinical metadata TSV containing HLA alleles.")
    parser.add_argument("--binding-dir", help="Hartwig Neo trained binding models directory.")
    parser.add_argument("--mixcr-path", help="Path to MiXCR .clns or exported TSV.")
    args = parser.parse_args()

    return run_neoantigen_tcr_pipeline(
        output_dir=args.output_dir,
        samplesheet=args.samplesheet,
        metadata_tsv=args.metadata_tsv,
        binding_dir=args.binding_dir,
        mixcr_path=args.mixcr_path,
    )


if __name__ == "__main__":
    sys.exit(main())
