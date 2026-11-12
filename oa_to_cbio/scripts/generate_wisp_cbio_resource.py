#!/usr/bin/env python3
"""
generate_wisp_cbio_resource.py
==============================
Generates cBioPortal Resource Framework files and interactive HTML panels
for plasma ctDNA / WISP multi-modal decomposition and longitudinal MRD monitoring.

Creates:
  - <import_dir>/reports/wisp/<patient_id>_wisp.html (interactive patient view report)
  - <import_dir>/reports/wisp/plots/<patient_id>_<sample_id>_wisp_multimodal.png
  - <import_dir>/meta_resource_definition.txt
  - <import_dir>/data_resource_definition.txt
  - <import_dir>/meta_resource_patient.txt
  - <import_dir>/data_resource_patient.txt
  - <import_dir>/meta_resource_sample.txt
  - <import_dir>/data_resource_sample.txt
"""

import os
import sys
import re
import csv
import json
import base64
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# Add root repo to path to allow importing compile_wgts_summary utilities
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:
    from scripts import compile_wgts_summary as cws
    HAS_CWS = True
except Exception as e:
    HAS_CWS = False
    print(f"Warning: compile_wgts_summary could not be imported ({e}); running with native fallbacks.", file=sys.stderr)


def safe_float(val, default=0.0):
    if val is None or val == "" or str(val).lower() in {"na", "nan", "null"}:
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=0):
    if val is None or val == "" or str(val).lower() in {"na", "nan", "null"}:
        return default
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return default


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


def parse_wisp_somatic_variants(snv_file):
    """Parse .wisp.somatic_variants.tsv file to extract somatic SNVs/InDels."""
    if not os.path.isfile(snv_file):
        return []
    variants = []
    try:
        with open(snv_file, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                gene = (r.get("Gene") or r.get("gene") or r.get("gene_name") or r.get("SYMBOL") or "").strip().upper()
                chrom = (r.get("Chromosome") or r.get("chr") or r.get("#CHROM") or "").strip()
                pos = safe_int(r.get("Position") or r.get("pos") or r.get("POS"))
                if not gene or gene == "NAN":
                    continue
                sample_ad = safe_int(r.get("SampleAD") or r.get("sample_ad") or r.get("AltCount") or 0)
                sample_dp = safe_int(r.get("SampleDP") or r.get("sample_dp") or r.get("Depth") or 0)
                tumor_ad = safe_int(r.get("TumorAD") or 0)
                tumor_dp = safe_int(r.get("TumorDP") or 0)
                vaf = (sample_ad / sample_dp) if sample_dp > 0 else 0.0
                vtype = (r.get("Type") or r.get("VariantType") or r.get("Consequence") or "Mutation").upper()
                variants.append({
                    "gene": gene,
                    "chrom": chrom,
                    "pos": pos,
                    "type": vtype,
                    "tumor_ad": tumor_ad,
                    "tumor_dp": tumor_dp,
                    "sample_ad": sample_ad,
                    "sample_dp": sample_dp,
                    "vaf": vaf,
                    "is_detected": sample_ad > 0,
                    "hgvs_protein": (r.get("canonicalHgvsProtein") or r.get("hgvsProtein") or "").strip(),
                    "hgvs_coding": (r.get("canonicalHgvsCoding") or r.get("hgvsCoding") or "").strip(),
                })
    except Exception as e:
        print(f"Warning: Failed parsing WISP somatic variants from {snv_file}: {e}", file=sys.stderr)
    return variants


def discover_wisp_data(in_dir, import_dir=None):
    """
    Scans in_dir for WISP summary, CNV, LOH, SNV, and matched tumor PURPLE files.
    Groups results by patient.
    """
    in_path = Path(in_dir).resolve()
    if not in_path.is_dir():
        return {}

    # 0. Load sample-to-patient mapping from clinical data if available
    sample_to_patient = {}
    if import_dir:
        clin_sample_file = Path(import_dir) / "data_clinical_sample.txt"
        if clin_sample_file.is_file():
            try:
                with open(clin_sample_file, "r", encoding="utf-8-sig") as f:
                    for line in f:
                        if line.startswith("#"): continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 2 and parts[0] != "SAMPLE_ID" and parts[1] != "PATIENT_ID":
                            # In cBioPortal, Col 1 is often SAMPLE_ID or PATIENT_ID
                            # Let's inspect headers if first row
                            pass
                # Standard reader
                with open(clin_sample_file, "r", encoding="utf-8-sig") as f:
                    lines = [l for l in f if not l.startswith("#")]
                    if lines:
                        reader = csv.DictReader(lines, delimiter="\t")
                        for r in reader:
                            s = r.get("SAMPLE_ID") or r.get("SAMPLE")
                            p = r.get("PATIENT_ID") or r.get("PATIENT")
                            if s and p:
                                sample_to_patient[s] = p
                                sample_to_patient[s.lower()] = p
            except Exception:
                pass

    EXCLUDE_DIRS = {
        "work", ".nextflow", "logs", "work_dir", "scratch", ".git",
        "pipeline_info", "sample_sheets", "all_fastq_links", "realigned_fastqs",
        "tmp", "temp"
    }

    wisp_summaries = []
    purple_purities = {}
    purple_cnvs = {}
    orange_reports = {}
    ichor_map = {}

    for root, dirs, files in os.walk(in_path, followlinks=True):
        dirs[:] = [d for d in dirs if d.lower() not in EXCLUDE_DIRS and not d.startswith(".")]
        root_p = Path(root)
        for f in files:
            if f.startswith("."):
                continue
            f_lower = f.lower()
            if f.endswith(".wisp.summary.tsv") or (f.endswith(".wisp.tsv") and not any(f.endswith(x) for x in [".snv.wisp.tsv", ".cnv.wisp.tsv", ".loh.wisp.tsv"])):
                wisp_summaries.append(root_p / f)
            elif f.endswith(".purple.purity.tsv"):
                stem = f.replace(".purple.purity.tsv", "")
                purple_purities[stem] = root_p / f
            elif f.endswith(".purple.cnv.somatic.tsv"):
                stem = f.replace(".purple.cnv.somatic.tsv", "")
                purple_cnvs[stem] = root_p / f
            elif f.endswith(".orange.pdf") or f.endswith(".orange.html"):
                stem = f.replace(".orange.pdf", "").replace(".orange.html", "")
                orange_reports[stem] = root_p / f
            elif f_lower.endswith("params.txt") and ("ichor" in str(root_p).lower() or "cna" in str(root_p).lower()):
                try:
                    with open(root_p / f, "r") as pf:
                        for line in pf:
                            if "TF:" in line or "Tumor Fraction:" in line:
                                tf_val = safe_float(line.split(":")[-1].strip())
                                sid = f.replace(".params.txt", "").replace("_params.txt", "")
                                ichor_map[sid] = tf_val
                                break
                except Exception:
                    pass

    if not wisp_summaries:
        return {}

    # Group by patient
    patients_data = defaultdict(lambda: {"plasmas": [], "tumors": []})

    for wf in wisp_summaries:
        raw_row = parse_wisp_summary_tsv(str(wf))
        if not raw_row:
            continue

        wisp_dir = wf.parent
        fname = wf.name
        stem = fname.replace(".wisp.summary.tsv", "").replace(".wisp.tsv", "")

        # Infer patient ID and sample ID
        patient_id = None
        if stem in sample_to_patient:
            patient_id = sample_to_patient[stem]
        elif stem.lower() in sample_to_patient:
            patient_id = sample_to_patient[stem.lower()]

        parts = wf.parts
        if not patient_id:
            for p in parts:
                if re.match(r"^(OV[-_]P[-_]\d+|TRACEOC[-_]\d+|BFIERCE[-_]\d+|DYNA[-_]P\d+|MOHCCN[-_][A-Za-z0-9]+|OCT[-_]\d+[-_]\d+|CHARM[-_]\d+)", p, re.IGNORECASE):
                    patient_id = p
                    break

        if not patient_id:
            # Nested check: patient_id / sample_id / wisp / file.tsv
            if len(parts) >= 4 and parts[-2] == "wisp":
                cand_sample = parts[-3]
                cand_patient = parts[-4]
                if any(k in cand_sample.lower() for k in ["plasma", "cfdna", "_t", "_fu", "-t", "cycle"]):
                    patient_id = cand_patient
                else:
                    patient_id = cand_sample
            elif len(parts) >= 3 and parts[-2] == "wisp":
                patient_id = parts[-3]
            elif len(parts) >= 2:
                patient_id = parts[-2]

        # Suffix strip check from stem (e.g. TEST-PAT-001_Plasma_T0 -> TEST-PAT-001)
        if not patient_id or any(k in str(patient_id).lower() for k in ["_plasma", "-plasma", "_cfdna", "-cfdna", "_t0", "_t1", "_t2"]):
            m_stem = re.match(r"^(.*?)([-_](?:plasma|cfdna|t[0-9]+|fu[0-9]+|baseline|cycle[-_]?[0-9]+).*)$", stem, re.IGNORECASE)
            if m_stem:
                patient_id = m_stem.group(1)

        patient_id = patient_id or "Patient"

        sample_id = stem
        if patient_id and stem.startswith(f"{patient_id}_"):
            sample_id = stem[len(patient_id) + 1:]

        # Extract timepoint / date
        timepoint = raw_row.get("Timepoint") or raw_row.get("timepoint") or ""
        sampling_date = raw_row.get("Date") or raw_row.get("date") or raw_row.get("SamplingDate") or ""
        if not timepoint:
            m_tp = re.search(r"[-_](t[0-9]+|fu[0-9]+|baseline|cycle[-_]?[0-9]+|post[-_]?op|progression)", sample_id, re.IGNORECASE)
            timepoint = m_tp.group(1).upper() if m_tp else "Plasma ctDNA"

        # Modality stats
        wisp_purity_pct = safe_float(raw_row.get("Purity") or raw_row.get("purity") or raw_row.get("WispPurity")) * 100.0
        wisp_lod_pct = safe_float(raw_row.get("LOD") or raw_row.get("lod") or raw_row.get("WispLOD")) * 100.0
        wisp_pval = safe_float(raw_row.get("PValue") or raw_row.get("p_value") or raw_row.get("WispPValue"))
        is_mrd_pos = raw_row.get("IsMRDPositive") or raw_row.get("is_mrd_positive") or raw_row.get("MRD_Status")
        if isinstance(is_mrd_pos, str):
            is_mrd_pos = is_mrd_pos.strip().lower() in {"true", "1", "positive", "pos"}
        else:
            is_mrd_pos = bool(is_mrd_pos) or (wisp_purity_pct > wisp_lod_pct and wisp_purity_pct > 0.05)

        snv_purity_pct = safe_float(raw_row.get("SnvPurity") or raw_row.get("snv_purity")) * 100.0
        snv_lod_pct = safe_float(raw_row.get("SnvLOD") or raw_row.get("snv_lod")) * 100.0
        snv_pval = safe_float(raw_row.get("SnvPValue") or raw_row.get("snv_pvalue"))
        snv_mrd = bool(raw_row.get("SnvMRD") or (snv_purity_pct > snv_lod_pct and snv_purity_pct > 0))

        cnv_purity_pct = safe_float(raw_row.get("CnvPurity") or raw_row.get("cnv_purity")) * 100.0
        cnv_lod_pct = safe_float(raw_row.get("CnvLOD") or raw_row.get("cnv_lod")) * 100.0
        cnv_pval = safe_float(raw_row.get("CnvPValue") or raw_row.get("cnv_pvalue"))
        cnv_mrd = bool(raw_row.get("CnvMRD") or (cnv_purity_pct > cnv_lod_pct and cnv_purity_pct > 0))

        loh_purity_pct = safe_float(raw_row.get("LohPurity") or raw_row.get("loh_purity")) * 100.0
        loh_lod_pct = safe_float(raw_row.get("LohLOD") or raw_row.get("loh_lod")) * 100.0
        loh_pval = safe_float(raw_row.get("LohPValue") or raw_row.get("loh_pvalue"))
        loh_mrd = bool(raw_row.get("LohMRD") or (loh_purity_pct > loh_lod_pct and loh_purity_pct > 0))

        # Somatic variants file
        snv_file = None
        for sf in [wisp_dir / f"{stem}.wisp.somatic_variants.tsv", wisp_dir / f"{sample_id}.wisp.somatic_variants.tsv"]:
            if sf.is_file():
                snv_file = sf
                break
        if not snv_file:
            for f in wisp_dir.glob("*.wisp.somatic_variants.tsv"):
                snv_file = f
                break

        somatic_vars = parse_wisp_somatic_variants(str(snv_file)) if snv_file else []

        # Matched ichorCNA
        ichor_tf = ichor_map.get(sample_id) or ichor_map.get(stem)

        plasma_record = {
            "sample_id": sample_id,
            "patient_id": patient_id,
            "stem": stem,
            "wisp_file": str(wf),
            "snv_file": str(snv_file) if snv_file else None,
            "timepoint": timepoint,
            "sampling_date": sampling_date,
            "wisp_purity_pct": wisp_purity_pct,
            "wisp_lod_pct": wisp_lod_pct,
            "wisp_pvalue": wisp_pval,
            "is_mrd_positive": is_mrd_pos,
            "snv_purity_pct": snv_purity_pct,
            "snv_lod_pct": snv_lod_pct,
            "snv_pvalue": snv_pval,
            "snv_mrd": snv_mrd,
            "cnv_purity_pct": cnv_purity_pct,
            "cnv_lod_pct": cnv_lod_pct,
            "cnv_pvalue": cnv_pval,
            "cnv_mrd": cnv_mrd,
            "loh_purity_pct": loh_purity_pct,
            "loh_lod_pct": loh_lod_pct,
            "loh_pvalue": loh_pval,
            "loh_mrd": loh_mrd,
            "ichor_tf": ichor_tf,
            "somatic_variants": somatic_vars,
            "multimodal_plot": None,
            "multimodal_plot_b64": None,
        }

        patients_data[patient_id]["plasmas"].append(plasma_record)

    # Attach matched tumor records if found
    for pid, pdata in patients_data.items():
        matched_tumor = None
        for stem, pp_path in purple_purities.items():
            if pid in stem or pid.replace("-", "_") in stem:
                matched_tumor = {
                    "tumor_id": stem,
                    "purple_purity_file": str(pp_path),
                    "purple_cnv_file": str(purple_cnvs.get(stem)),
                }
                break
        if matched_tumor:
            pdata["tumors"].append(matched_tumor)

        # Attach Orange reports for this patient
        pdata["orange_reports"] = []
        for stem, orange_p in orange_reports.items():
            if pid in stem or pid.replace("-", "_") in stem or (sample_to_patient.get(stem) == pid):
                pdata["orange_reports"].append({
                    "sample_id": stem,
                    "path": str(orange_p),
                    "ext": orange_p.suffix
                })

    # Also capture patients who only have Orange reports (no WISP)
    for stem, orange_p in orange_reports.items():
        pid = sample_to_patient.get(stem) or sample_to_patient.get(stem.lower())
        if not pid:
            m_stem = re.match(r"^(.*?)([-_](?:tumor|tumour|t|primary).*)$", stem, re.IGNORECASE)
            pid = m_stem.group(1) if m_stem else stem
        if pid not in patients_data:
            patients_data[pid]["orange_reports"].append({
                "sample_id": stem,
                "path": str(orange_p),
                "ext": orange_p.suffix
            })

    return patients_data


def generate_interactive_patient_wisp_html(patient_id, p_rec, out_html_path):
    """
    Renders a modern, self-contained interactive patient view HTML dashboard
    matching the cBioPortal 'Summary' tab styling, with:
      - Top Clinical & MRD KPI Strip
      - Interactive Longitudinal Trajectory Timeline Chart (Log10/Linear toggle, LOD shading, hover tooltips)
      - Timepoint Selector Tabs / Pills
      - 4-Card Multi-Modal Breakdown (CNV, LOH, SNV, 1X ichorCNA)
      - Multi-Modal Decomposition Figure Viewer with zoom modal
      - Interactive Somatic Variant Table with live search and sorting
      - Zero external internet dependencies
    """
    plasmas = sorted(p_rec.get("plasmas", []), key=lambda p: p.get("sampling_date") or p.get("timepoint") or p.get("sample_id"))
    if not plasmas:
        return

    # Calculate Patient-level longitudinal summary metrics
    latest_pl = plasmas[-1]
    is_mrd_pos = any(p.get("is_mrd_positive", False) for p in plasmas)
    latest_mrd_pos = latest_pl.get("is_mrd_positive", False)
    all_purities = [p.get("wisp_purity_pct", 0.0) for p in plasmas]
    peak_purity = max(all_purities) if all_purities else 0.0
    latest_purity = latest_pl.get("wisp_purity_pct", 0.0)
    nadir_purity = min(all_purities) if all_purities else 0.0

    tumors = p_rec.get("tumors", [])
    tumor_id = tumors[0]["tumor_id"] if tumors else "Primary Tumor"

    # Convert payload to JSON for embedded offline execution
    plasmas_json = json.dumps(plasmas)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Liquid Biopsy (WISP) — {patient_id}</title>
  <style>
    :root {{
      --primary: #0284c7;
      --primary-dark: #0369a1;
      --text: #0f172a;
      --text-muted: #64748b;
      --bg: #f8fafc;
      --card-bg: #ffffff;
      --border: #e2e8f0;
      --border-dark: #cbd5e1;
      --danger: #dc2626;
      --danger-bg: #fef2f2;
      --success: #16a34a;
      --success-bg: #f0fdf4;
      --warning: #d97706;
      --warning-bg: #fffbeb;
      --purple: #7c3aed;
      --purple-bg: #f5f3ff;
      --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    }}

    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: var(--font-sans);
      background-color: var(--bg);
      color: var(--text);
      line-height: 1.5;
      padding: 1.25rem;
      -webkit-font-smoothing: antialiased;
    }}

    .container {{
      max-width: 1400px;
      margin: 0 auto;
    }}

    /* Top Patient Header */
    .patient-header {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1.25rem 1.5rem;
      margin-bottom: 1.25rem;
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 1rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    }}
    .patient-title-group {{
      display: flex;
      align-items: center;
      gap: 12px;
      flex-wrap: wrap;
    }}
    .patient-title {{
      font-size: 1.35rem;
      font-weight: 700;
      color: var(--text);
      letter-spacing: -0.01em;
    }}
    .badge {{
      display: inline-flex;
      align-items: center;
      padding: 0.25rem 0.65rem;
      border-radius: 9999px;
      font-size: 0.75rem;
      font-weight: 700;
      letter-spacing: 0.03em;
      text-transform: uppercase;
    }}
    .badge-pos {{ background: var(--danger-bg); color: var(--danger); border: 1px solid #fecaca; }}
    .badge-neg {{ background: var(--success-bg); color: var(--success); border: 1px solid #bbf7d0; }}
    .badge-info {{ background: #eff6ff; color: var(--primary); border: 1px solid #bfdbfe; }}

    /* KPI Summary Strip */
    .kpi-strip {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 1rem;
      margin-bottom: 1.25rem;
    }}
    .kpi-card {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1rem;
      box-shadow: 0 1px 2px rgba(0,0,0,0.03);
    }}
    .kpi-label {{
      font-size: 0.75rem;
      font-weight: 600;
      color: var(--text-muted);
      text-transform: uppercase;
      letter-spacing: 0.05em;
      margin-bottom: 0.25rem;
    }}
    .kpi-value {{
      font-size: 1.45rem;
      font-weight: 700;
      color: var(--text);
    }}
    .kpi-subtext {{
      font-size: 0.75rem;
      color: var(--text-muted);
      margin-top: 0.2rem;
    }}

    /* Main Chart Card */
    .card {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1.25rem;
      margin-bottom: 1.25rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    }}
    .card-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 1rem;
      padding-bottom: 0.75rem;
      border-bottom: 1px solid var(--border);
      flex-wrap: wrap;
      gap: 0.75rem;
    }}
    .card-title {{
      font-size: 1.05rem;
      font-weight: 700;
      color: var(--text);
      display: flex;
      align-items: center;
      gap: 8px;
    }}

    /* Chart Controls */
    .btn-group {{
      display: inline-flex;
      border-radius: 6px;
      overflow: hidden;
      border: 1px solid var(--border-dark);
    }}
    .btn {{
      background: #ffffff;
      border: none;
      padding: 0.35rem 0.75rem;
      font-size: 0.8rem;
      font-weight: 600;
      color: var(--text-muted);
      cursor: pointer;
      transition: all 0.15s ease;
    }}
    .btn:not(:last-child) {{ border-right: 1px solid var(--border-dark); }}
    .btn:hover {{ background: #f1f5f9; color: var(--text); }}
    .btn.active {{ background: var(--primary); color: #ffffff; }}

    /* SVG Chart Styles */
    .chart-container {{
      position: relative;
      width: 100%;
      height: 320px;
      background: #ffffff;
    }}
    svg.trajectory-chart {{
      width: 100%;
      height: 100%;
      overflow: visible;
    }}
    .chart-grid line {{ stroke: #f1f5f9; stroke-width: 1; }}
    .chart-axis text {{ fill: #64748b; font-size: 11px; font-family: var(--font-sans); }}
    .lod-zone {{ fill: #f8fafc; opacity: 0.9; }}
    .lod-line {{ stroke: #cbd5e1; stroke-dasharray: 4 4; stroke-width: 1.5; }}
    .chart-line {{ fill: none; stroke-width: 2.5; stroke-linecap: round; stroke-linejoin: round; }}
    .chart-dot {{
      stroke-width: 2;
      cursor: pointer;
      transition: r 0.15s ease, stroke-width 0.15s ease;
    }}
    .chart-dot:hover {{ r: 7.5; stroke-width: 3.5; }}
    .chart-dot.active {{ r: 8.5; stroke: #0f172a !important; stroke-width: 3; }}

    /* Tooltip */
    .chart-tooltip {{
      position: absolute;
      background: rgba(15, 23, 42, 0.95);
      color: #ffffff;
      padding: 8px 12px;
      border-radius: 6px;
      font-size: 0.75rem;
      pointer-events: none;
      display: none;
      z-index: 100;
      box-shadow: 0 4px 6px -1px rgba(0,0,0,0.2);
      line-height: 1.4;
      white-space: nowrap;
    }}

    /* Timepoint Pill Tabs */
    .timepoint-nav {{
      display: flex;
      gap: 0.5rem;
      overflow-x: auto;
      padding-bottom: 0.25rem;
      margin-bottom: 1.25rem;
    }}
    .timepoint-pill {{
      background: var(--card-bg);
      border: 1px solid var(--border-dark);
      border-radius: 6px;
      padding: 0.5rem 0.9rem;
      font-size: 0.82rem;
      font-weight: 600;
      color: var(--text-muted);
      cursor: pointer;
      white-space: nowrap;
      display: flex;
      align-items: center;
      gap: 8px;
      transition: all 0.15s ease;
    }}
    .timepoint-pill:hover {{ background: #f8fafc; color: var(--text); border-color: var(--primary); }}
    .timepoint-pill.active {{
      background: #eff6ff;
      border-color: var(--primary);
      color: var(--primary-dark);
      box-shadow: 0 1px 2px rgba(2, 132, 199, 0.15);
    }}

    /* 4-Card Modality Grid */
    .modality-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 1rem;
      margin-bottom: 1.25rem;
    }}
    .modality-card {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1rem;
      box-shadow: 0 1px 2px rgba(0,0,0,0.03);
    }}
    .modality-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 0.5rem;
    }}
    .modality-title {{ font-size: 0.92rem; font-weight: 700; }}
    .modality-desc {{ font-size: 0.75rem; color: var(--text-muted); margin-bottom: 0.75rem; line-height: 1.35; }}
    .modality-stats {{
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 0.6rem 0.75rem;
      font-family: var(--font-mono);
      font-size: 0.78rem;
      display: grid;
      gap: 0.25rem;
    }}

    /* Image figure container */
    .plot-figure-box {{
      text-align: center;
      background: #0b0f14;
      padding: 0.75rem;
      border-radius: 8px;
      border: 1px solid #1f2937;
      margin-bottom: 1.25rem;
    }}
    .plot-figure-img {{
      width: 100%;
      max-height: 65vh;
      object-fit: contain;
      border-radius: 4px;
      display: block;
      margin: 0 auto;
      cursor: zoom-in;
    }}

    /* Somatic Variants Table */
    .table-container {{
      overflow-x: auto;
      border: 1px solid var(--border);
      border-radius: 6px;
    }}
    table.data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.82rem;
      text-align: left;
    }}
    table.data-table th {{
      background: #f8fafc;
      color: var(--text);
      font-weight: 700;
      padding: 0.65rem 0.85rem;
      border-bottom: 1px solid var(--border);
      cursor: pointer;
      user-select: none;
    }}
    table.data-table td {{
      padding: 0.6rem 0.85rem;
      border-bottom: 1px solid #f1f5f9;
      color: var(--text);
    }}
    table.data-table tr:hover {{ background: #f8fafc; }}

    /* Zoom Modal */
    .modal-overlay {{
      display: none;
      position: fixed;
      inset: 0;
      background: rgba(0, 0, 0, 0.85);
      z-index: 1000;
      align-items: center;
      justify-content: center;
      padding: 2rem;
    }}
    .modal-overlay.open {{ display: flex; }}
    .modal-img {{
      max-width: 95vw;
      max-height: 92vh;
      object-fit: contain;
      border-radius: 6px;
      box-shadow: 0 20px 25px -5px rgba(0,0,0,0.5);
    }}
    .modal-close {{
      position: absolute;
      top: 1.5rem;
      right: 2rem;
      color: #ffffff;
      font-size: 2rem;
      font-weight: 700;
      cursor: pointer;
    }}
  </style>
</head>
<body>

<div class="container">

  <!-- 1. Patient Header -->
  <div class="patient-header">
    <div class="patient-title-group">
      <span style="font-size:1.5rem;">🩸</span>
      <div>
        <h1 class="patient-title">{patient_id} — Liquid Biopsy (WISP ctDNA)</h1>
        <div style="font-size:0.8rem; color:var(--text-muted); margin-top:2px;">
          Primary Baseline: <strong style="color:var(--text);">{tumor_id}</strong> &bull; Total Plasma Timepoints: <strong>{len(plasmas)}</strong>
        </div>
      </div>
      <span class="badge {'badge-pos' if latest_mrd_pos else 'badge-neg'}" style="margin-left:8px; font-size:0.82rem; padding:0.35rem 0.85rem;">
        {'MRD POSITIVE' if latest_mrd_pos else 'MRD NEGATIVE'}
      </span>
    </div>
    <div style="display:flex; gap:8px;">
      <button class="btn" onclick="exportDataTSV()" style="border:1px solid var(--border-dark); border-radius:6px;">📥 Export TSV</button>
      <button class="btn" onclick="window.print()" style="border:1px solid var(--border-dark); border-radius:6px;">🖨️ Print View</button>
    </div>
  </div>

  <!-- 2. KPI Summary Strip -->
  <div class="kpi-strip">
    <div class="kpi-card" style="border-top:3px solid var(--primary);">
      <div class="kpi-label">Latest Tumor Fraction</div>
      <div class="kpi-value" style="color:var(--primary);">{latest_purity:.3f}%</div>
      <div class="kpi-subtext">{latest_pl.get('timepoint', 'Latest')} ({latest_pl.get('sampling_date') or 'Sample'})</div>
    </div>
    <div class="kpi-card" style="border-top:3px solid var(--danger);">
      <div class="kpi-label">Peak ctDNA Fraction</div>
      <div class="kpi-value" style="color:var(--danger);">{peak_purity:.3f}%</div>
      <div class="kpi-subtext">Across {len(plasmas)} monitored samples</div>
    </div>
    <div class="kpi-card" style="border-top:3px solid var(--success);">
      <div class="kpi-label">Nadir ctDNA Fraction</div>
      <div class="kpi-value" style="color:var(--success);">{nadir_purity:.3f}%</div>
      <div class="kpi-subtext">Lowest observed ctDNA level</div>
    </div>
    <div class="kpi-card" style="border-top:3px solid var(--purple);">
      <div class="kpi-label">Assay Modalities</div>
      <div class="kpi-value" style="color:var(--purple); font-size:1.15rem;">SNV + CNV + LOH</div>
      <div class="kpi-subtext">Multi-modal WISP decomposition</div>
    </div>
  </div>

  <!-- 3. Longitudinal MRD Trajectory Card -->
  <div class="card">
    <div class="card-header">
      <div class="card-title">
        <span>📈</span> Longitudinal ctDNA Monitoring Trajectory
      </div>
      <div style="display:flex; gap:12px; align-items:center;">
        <div class="btn-group">
          <button class="btn active" id="btn-scale-log" onclick="setScale('log')">Log10 Scale</button>
          <button class="btn" id="btn-scale-lin" onclick="setScale('linear')">Linear Scale</button>
        </div>
      </div>
    </div>
    <div class="chart-container">
      <div id="chart-tooltip" class="chart-tooltip"></div>
      <svg id="trajectory-chart" class="trajectory-chart"></svg>
    </div>
    <div style="display:flex; justify-content:space-between; align-items:center; margin-top:8px; font-size:0.75rem; color:var(--text-muted);">
      <div>
        <span style="display:inline-block; width:10px; height:10px; background:#0284c7; border-radius:50%; margin-right:4px;"></span> WISP ctDNA (%)
        <span style="display:inline-block; width:10px; height:10px; background:#16a34a; border-radius:50%; margin-left:12px; margin-right:4px;"></span> Somatic SNV
        <span style="display:inline-block; width:10px; height:10px; background:#7c3aed; border-radius:50%; margin-left:12px; margin-right:4px;"></span> CNV Ratio
        <span style="display:inline-block; width:12px; height:6px; background:#cbd5e1; margin-left:12px; margin-right:4px;"></span> Grey Shading: Below Limit of Detection (&lt; LOD)
      </div>
      <div><em>Click any dot to inspect timepoint details below</em></div>
    </div>
  </div>

  <!-- 4. Timepoint Selector Nav -->
  <div style="margin-bottom:0.5rem; font-size:0.85rem; font-weight:700; color:var(--text);">Select Longitudinal Timepoint:</div>
  <div class="timepoint-nav" id="timepoint-nav"></div>

  <!-- 5. Active Timepoint Decomposition Details -->
  <div class="card" id="timepoint-detail-card">
    <div class="card-header">
      <div class="card-title" id="active-timepoint-title">
        <span>🔬</span> Multi-Modal Decomposition Breakdown
      </div>
      <div id="active-timepoint-badge"></div>
    </div>

    <!-- 4-Card Modality Grid -->
    <div class="modality-grid">
      <!-- Card 1: CNV -->
      <div class="modality-card" style="border-top:3px solid #0284c7;">
        <div class="modality-header">
          <span class="modality-title" style="color:#0284c7;">1. Copy Number (CNV)</span>
          <span class="badge" id="m-cnv-badge"></span>
        </div>
        <div class="modality-desc">Compares GC-corrected plasma depth ratio against tumor segment copy numbers via regression.</div>
        <div class="modality-stats">
          <div>Purity:    <strong id="m-cnv-purity" style="color:#0284c7;">-</strong></div>
          <div>Assay LOD: <span id="m-cnv-lod" style="color:#7c3aed;">-</span></div>
          <div>p-value:   <span id="m-cnv-pval">-</span></div>
        </div>
      </div>

      <!-- Card 2: LOH -->
      <div class="modality-card" style="border-top:3px solid #8b5cf6;">
        <div class="modality-header">
          <span class="modality-title" style="color:#7c3aed;">2. Loss of Heterozygosity</span>
          <span class="badge" id="m-loh-badge"></span>
        </div>
        <div class="modality-desc">Aggregates phased germline SNPs across tumor LOH segments to measure deviation from 50/50 baseline.</div>
        <div class="modality-stats">
          <div>Purity:    <strong id="m-loh-purity" style="color:#7c3aed;">-</strong></div>
          <div>Assay LOD: <span id="m-loh-lod" style="color:#d97706;">-</span></div>
          <div>p-value:   <span id="m-loh-pval">-</span></div>
        </div>
      </div>

      <!-- Card 3: SNV -->
      <div class="modality-card" style="border-top:3px solid #16a34a;">
        <div class="modality-header">
          <span class="modality-title" style="color:#16a34a;">3. Somatic SNVs/InDels</span>
          <span class="badge" id="m-snv-badge"></span>
        </div>
        <div class="modality-desc">Tracks tumor somatic mutations in plasma cfDNA against background sequencing error.</div>
        <div class="modality-stats">
          <div>Purity:    <strong id="m-snv-purity" style="color:#16a34a;">-</strong></div>
          <div>Assay LOD: <span id="m-snv-lod" style="color:#0284c7;">-</span></div>
          <div>p-value:   <span id="m-snv-pval">-</span></div>
        </div>
      </div>

      <!-- Card 4: 1X ichorCNA -->
      <div class="modality-card" style="border-top:3px solid #d97706;">
        <div class="modality-header">
          <span class="modality-title" style="color:#d97706;">4. 1X ichorCNA</span>
          <span class="badge" id="m-ichor-badge"></span>
        </div>
        <div class="modality-desc">Ultra-low pass whole-genome sequencing copy-number tumor fraction estimation.</div>
        <div class="modality-stats">
          <div>Tumor Fraction: <strong id="m-ichor-tf" style="color:#d97706;">-</strong></div>
          <div>Assay Method:   <span>sWGS (1X)</span></div>
          <div>Concordance:    <span id="m-ichor-status">-</span></div>
        </div>
      </div>
    </div>

    <!-- Multi-Modal Plot Figure Box -->
    <div style="font-size:0.9rem; font-weight:700; color:var(--text); margin-bottom:6px;">
      Genome-Wide Multi-Modal Decomposition Figure:
    </div>
    <div class="plot-figure-box" id="plot-figure-container">
      <img id="active-plot-img" class="plot-figure-img" alt="Multi-Modal WISP Decomposition Plot" onclick="openZoomModal(this.src)" />
      <div style="font-size:0.75rem; color:#9ca3af; margin-top:8px; text-align:right;">
        <em>Click figure to view in full-screen zoom modal</em>
      </div>
    </div>

    <!-- Somatic Variants Table -->
    <div style="display:flex; justify-content:space-between; align-items:center; margin-top:1.25rem; margin-bottom:0.6rem; flex-wrap:wrap; gap:8px;">
      <div style="font-size:0.95rem; font-weight:700; color:var(--text);">
        Tracked Somatic Variants in Plasma (<span id="variant-count">0</span> variants)
      </div>
      <div>
        <input type="text" id="variant-search" placeholder="🔍 Search gene, locus, type..." 
               style="padding:0.35rem 0.65rem; font-size:0.8rem; border:1px solid var(--border-dark); border-radius:6px; width:220px;"
               onkeyup="filterVariantTable()" />
      </div>
    </div>
    <div class="table-container">
      <table class="data-table" id="variants-table">
        <thead>
          <tr>
            <th onclick="sortTable(0)">Gene</th>
            <th onclick="sortTable(1)">Locus</th>
            <th onclick="sortTable(2)">Type / Consequence</th>
            <th onclick="sortTable(3)">Tumor Reads (VAF)</th>
            <th onclick="sortTable(4)">Plasma Reads (AD/DP)</th>
            <th onclick="sortTable(5)">Plasma VAF</th>
            <th onclick="sortTable(6)">Detection Call</th>
          </tr>
        </thead>
        <tbody id="variants-tbody">
          <!-- Populated by JS -->
        </tbody>
      </table>
    </div>

  </div>

</div>

<!-- Full Screen Modal -->
<div class="modal-overlay" id="zoom-modal" onclick="closeZoomModal()">
  <span class="modal-close">&times;</span>
  <img id="modal-img" class="modal-img" alt="Expanded WISP Plot" />
</div>

<script>
  const PLASMAS = {plasmas_json};
  let activeIndex = 0;
  let currentScale = 'log';

  function init() {{
    renderTimepointNav();
    renderChart();
    selectTimepoint(0);
  }}

  function setScale(scale) {{
    currentScale = scale;
    document.getElementById('btn-scale-log').classList.toggle('active', scale === 'log');
    document.getElementById('btn-scale-lin').classList.toggle('active', scale === 'linear');
    renderChart();
  }}

  function renderTimepointNav() {{
    const nav = document.getElementById('timepoint-nav');
    nav.innerHTML = '';
    PLASMAS.forEach((pl, idx) => {{
      const pill = document.createElement('div');
      pill.className = `timepoint-pill ${{idx === activeIndex ? 'active' : ''}}`;
      pill.id = `tp-pill-${{idx}}`;
      const dateStr = pl.sampling_date ? ` (${{pl.sampling_date}})` : '';
      const badge = pl.is_mrd_positive 
        ? '<span class="badge badge-pos" style="font-size:0.65rem; padding:1px 5px;">POS</span>'
        : '<span class="badge badge-neg" style="font-size:0.65rem; padding:1px 5px;">&lt; LOD</span>';
      pill.innerHTML = `<span>📅 ${{pl.timepoint || pl.sample_id}}${{dateStr}}</span> ${{badge}}`;
      pill.onclick = () => selectTimepoint(idx);
      nav.appendChild(pill);
    }});
  }}

  function selectTimepoint(idx) {{
    activeIndex = idx;
    document.querySelectorAll('.timepoint-pill').forEach((p, i) => {{
      p.classList.toggle('active', i === idx);
    }});

    // Update active dot in chart
    document.querySelectorAll('.chart-dot').forEach((dot, i) => {{
      dot.classList.toggle('active', i === idx);
    }});

    const pl = PLASMAS[idx];
    if (!pl) return;

    // Header info
    document.getElementById('active-timepoint-title').innerHTML = `
      <span>🔬</span> Timepoint: ${{pl.timepoint || pl.sample_id}} <span style="font-weight:normal; font-size:0.85rem; color:var(--text-muted); margin-left:6px;">(Sample ID: ${{pl.sample_id}})</span>
    `;
    document.getElementById('active-timepoint-badge').innerHTML = pl.is_mrd_positive
      ? '<span class="badge badge-pos">MRD POSITIVE</span>'
      : '<span class="badge badge-neg">MRD NEGATIVE (&lt; LOD)</span>';

    // Cards
    // 1. CNV
    document.getElementById('m-cnv-badge').className = `badge ${{pl.cnv_mrd ? 'badge-pos' : 'badge-neg'}}`;
    document.getElementById('m-cnv-badge').textContent = pl.cnv_mrd ? 'POSITIVE' : 'Negative';
    document.getElementById('m-cnv-purity').textContent = `${{pl.cnv_purity_pct.toFixed(4)}}%`;
    document.getElementById('m-cnv-lod').textContent = `${{pl.cnv_lod_pct.toFixed(4)}}%`;
    document.getElementById('m-cnv-pval').textContent = pl.cnv_pvalue ? pl.cnv_pvalue.toExponential(2) : 'NA';

    // 2. LOH
    document.getElementById('m-loh-badge').className = `badge ${{pl.loh_mrd ? 'badge-pos' : 'badge-neg'}}`;
    document.getElementById('m-loh-badge').textContent = pl.loh_mrd ? 'POSITIVE' : 'Negative';
    document.getElementById('m-loh-purity').textContent = `${{pl.loh_purity_pct.toFixed(4)}}%`;
    document.getElementById('m-loh-lod').textContent = `${{pl.loh_lod_pct.toFixed(4)}}%`;
    document.getElementById('m-loh-pval').textContent = pl.loh_pvalue ? pl.loh_pvalue.toExponential(2) : 'NA';

    // 3. SNV
    document.getElementById('m-snv-badge').className = `badge ${{pl.snv_mrd ? 'badge-pos' : 'badge-neg'}}`;
    document.getElementById('m-snv-badge').textContent = pl.snv_mrd ? 'POSITIVE' : 'Negative';
    document.getElementById('m-snv-purity').textContent = `${{pl.snv_purity_pct.toFixed(4)}}%`;
    document.getElementById('m-snv-lod').textContent = `${{pl.snv_lod_pct.toFixed(4)}}%`;
    document.getElementById('m-snv-pval').textContent = pl.snv_pvalue ? pl.snv_pvalue.toExponential(2) : 'NA';

    // 4. ichorCNA
    if (pl.ichor_tf !== null && pl.ichor_tf !== undefined) {{
      const tfPct = pl.ichor_tf * 100.0;
      document.getElementById('m-ichor-badge').className = `badge ${{tfPct >= 3.0 ? 'badge-pos' : 'badge-neg'}}`;
      document.getElementById('m-ichor-badge').textContent = tfPct >= 3.0 ? 'DETECTED' : 'LOW/SUB';
      document.getElementById('m-ichor-tf').textContent = `${{tfPct.toFixed(2)}}%`;
      document.getElementById('m-ichor-status').textContent = 'Concordant';
    }} else {{
      document.getElementById('m-ichor-badge').className = 'badge badge-info';
      document.getElementById('m-ichor-badge').textContent = 'NA';
      document.getElementById('m-ichor-tf').textContent = 'Not run / NA';
      document.getElementById('m-ichor-status').textContent = 'NA';
    }}

    // Image Figure
    const imgEl = document.getElementById('active-plot-img');
    const plotSrc = pl.multimodal_plot_b64 || pl.multimodal_plot || `plots/${{pl.patient_id}}_${{pl.sample_id}}_wisp_multimodal.png`;
    imgEl.src = plotSrc;
    imgEl.onerror = () => {{
      imgEl.onerror = null;
      imgEl.parentElement.innerHTML = '<div style="padding:2.5rem; color:#9ca3af; font-size:0.85rem;">Figure plot image will render here when generated via matplotlib.</div>';
    }};

    // Somatic Variants Table
    renderVariantsTable(pl.somatic_variants || []);
  }}

  function renderVariantsTable(variants) {{
    const tbody = document.getElementById('variants-tbody');
    document.getElementById('variant-count').textContent = variants.length;
    tbody.innerHTML = '';
    if (variants.length === 0) {{
      tbody.innerHTML = '<tr><td colspan="7" style="text-align:center; color:var(--text-muted); padding:1.5rem;">No somatic mutations tracked for this sample.</td></tr>';
      return;
    }}
    variants.forEach(v => {{
      const tr = document.createElement('tr');
      const vafPct = (v.vaf * 100).toFixed(2);
      const isDet = v.is_detected;
      tr.innerHTML = `
        <td><strong style="color:var(--primary);">${{v.gene}}</strong></td>
        <td style="font-family:var(--font-mono); font-size:0.75rem;">${{v.chrom}}:${{v.pos}}</td>
        <td><span style="font-size:0.75rem;">${{v.type}}</span></td>
        <td style="font-family:var(--font-mono);">${{v.tumor_ad}} / ${{v.tumor_dp}}</td>
        <td style="font-family:var(--font-mono);">${{v.sample_ad}} / ${{v.sample_dp}}</td>
        <td><strong style="color:${{isDet ? 'var(--primary)' : 'var(--text-muted)'}};">${{vafPct}}%</strong></td>
        <td>
          <span class="badge ${{isDet ? 'badge-pos' : 'badge-neg'}}" style="font-size:0.7rem; padding:1px 6px;">
            ${{isDet ? 'DETECTED' : '&lt; LOD'}}
          </span>
        </td>
      `;
      tbody.appendChild(tr);
    }});
  }}

  function filterVariantTable() {{
    const q = document.getElementById('variant-search').value.toLowerCase();
    const rows = document.querySelectorAll('#variants-tbody tr');
    rows.forEach(r => {{
      const txt = r.textContent.toLowerCase();
      r.style.display = txt.includes(q) ? '' : 'none';
    }});
  }}

  function renderChart() {{
    const svg = document.getElementById('trajectory-chart');
    svg.innerHTML = '';
    const W = svg.clientWidth || 900;
    const H = 300;
    const padL = 60, padR = 40, padT = 25, padB = 40;
    const chartW = W - padL - padR;
    const chartH = H - padT - padB;

    const N = PLASMAS.length;
    if (N === 0) return;

    // Y mapping
    const minY_log = 0.01;
    const maxY_log = Math.max(100.0, ...PLASMAS.map(p => p.wisp_purity_pct * 1.5 || 10.0));
    const maxY_lin = Math.max(10.0, ...PLASMAS.map(p => p.wisp_purity_pct * 1.2 || 10.0));

    function getY(val) {{
      if (currentScale === 'log') {{
        const safeVal = Math.max(val || minY_log, minY_log);
        const logMin = Math.log10(minY_log);
        const logMax = Math.log10(maxY_log);
        const norm = (Math.log10(safeVal) - logMin) / (logMax - logMin);
        return padT + chartH * (1 - norm);
      }} else {{
        const safeVal = Math.max(0, val || 0);
        const norm = safeVal / maxY_lin;
        return padT + chartH * (1 - norm);
      }}
    }}

    function getX(idx) {{
      if (N === 1) return padL + chartW / 2;
      return padL + (idx / (N - 1)) * chartW;
    }}

    // 1. Shaded LOD region (< 0.05%)
    const lodCutoff = 0.05;
    const lodY = getY(lodCutoff);
    const bottomY = padT + chartH;

    const lodRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    lodRect.setAttribute('x', padL);
    lodRect.setAttribute('y', lodY);
    lodRect.setAttribute('width', chartW);
    lodRect.setAttribute('height', Math.max(0, bottomY - lodY));
    lodRect.setAttribute('class', 'lod-zone');
    svg.appendChild(lodRect);

    // LOD dashed line
    const lodLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    lodLine.setAttribute('x1', padL);
    lodLine.setAttribute('x2', padL + chartW);
    lodLine.setAttribute('y1', lodY);
    lodLine.setAttribute('y2', lodY);
    lodLine.setAttribute('class', 'lod-line');
    svg.appendChild(lodLine);

    // LOD text label
    const lodText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    lodText.setAttribute('x', padL + 8);
    lodText.setAttribute('y', lodY - 5);
    lodText.setAttribute('fill', '#94a3b8');
    lodText.setAttribute('font-size', '10px');
    lodText.setAttribute('font-weight', '600');
    lodText.textContent = 'Assay Limit of Detection (LOD ~0.05%)';
    svg.appendChild(lodText);

    // 2. Grid & Y-Ticks
    const yTicks = currentScale === 'log' ? [0.01, 0.1, 1.0, 10.0, 100.0] : [0, 2, 5, 10, 20, 50].filter(v => v <= maxY_lin);
    yTicks.forEach(tickVal => {{
      const y = getY(tickVal);
      const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
      line.setAttribute('x1', padL);
      line.setAttribute('x2', padL + chartW);
      line.setAttribute('y1', y);
      line.setAttribute('y2', y);
      line.setAttribute('stroke', '#f1f5f9');
      svg.appendChild(line);

      const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      label.setAttribute('x', padL - 8);
      label.setAttribute('y', y + 3);
      label.setAttribute('text-anchor', 'end');
      label.setAttribute('fill', '#64748b');
      label.setAttribute('font-size', '10px');
      label.textContent = tickVal >= 1 ? `${{tickVal}}%` : `${{tickVal}}%`;
      svg.appendChild(label);
    }});

    // 3. Connect line for WISP ctDNA
    let dPath = '';
    PLASMAS.forEach((pl, idx) => {{
      const x = getX(idx);
      const y = getY(pl.wisp_purity_pct);
      dPath += (idx === 0 ? `M ${{x}} ${{y}}` : ` L ${{x}} ${{y}}`);
    }});

    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    path.setAttribute('d', dPath);
    path.setAttribute('class', 'chart-line');
    path.setAttribute('stroke', '#0284c7');
    svg.appendChild(path);

    // 4. Data Points & Tooltips
    const tip = document.getElementById('chart-tooltip');
    PLASMAS.forEach((pl, idx) => {{
      const x = getX(idx);
      const y = getY(pl.wisp_purity_pct);

      // X label
      const xLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      xLabel.setAttribute('x', x);
      xLabel.setAttribute('y', padT + chartH + 20);
      xLabel.setAttribute('text-anchor', 'middle');
      xLabel.setAttribute('fill', '#475569');
      xLabel.setAttribute('font-size', '11px');
      xLabel.setAttribute('font-weight', '600');
      xLabel.textContent = pl.timepoint || `T${{idx}}`;
      svg.appendChild(xLabel);

      if (pl.sampling_date) {{
        const dLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        dLabel.setAttribute('x', x);
        dLabel.setAttribute('y', padT + chartH + 34);
        dLabel.setAttribute('text-anchor', 'middle');
        dLabel.setAttribute('fill', '#94a3b8');
        dLabel.setAttribute('font-size', '9.5px');
        dLabel.textContent = pl.sampling_date;
        svg.appendChild(dLabel);
      }}

      // Dot
      const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
      circle.setAttribute('cx', x);
      circle.setAttribute('cy', y);
      circle.setAttribute('r', 5.5);
      circle.setAttribute('fill', pl.is_mrd_positive ? '#0284c7' : '#16a34a');
      circle.setAttribute('stroke', '#ffffff');
      circle.setAttribute('class', `chart-dot ${{idx === activeIndex ? 'active' : ''}}`);

      circle.onclick = () => selectTimepoint(idx);

      circle.onmouseenter = (e) => {{
        tip.style.display = 'block';
        tip.style.left = `${{e.pageX + 10}}px`;
        tip.style.top = `${{e.pageY - 25}}px`;
        tip.innerHTML = `
          <strong>${{pl.timepoint || pl.sample_id}}</strong> (${{pl.sampling_date || 'Plasma'}})<br/>
          WISP ctDNA: <strong>${{pl.wisp_purity_pct.toFixed(4)}}%</strong><br/>
          SNV Purity: ${{pl.snv_purity_pct.toFixed(4)}}% &bull; LOD: ${{pl.snv_lod_pct.toFixed(4)}}%<br/>
          CNV Purity: ${{pl.cnv_purity_pct.toFixed(4)}}% &bull; LOH: ${{pl.loh_purity_pct.toFixed(4)}}%<br/>
          Status: <strong>${{pl.is_mrd_positive ? 'MRD POSITIVE' : 'MRD NEGATIVE'}}</strong>
        `;
      }};
      circle.onmousemove = (e) => {{
        tip.style.left = `${{e.pageX + 10}}px`;
        tip.style.top = `${{e.pageY - 25}}px`;
      }};
      circle.onmouseleave = () => {{
        tip.style.display = 'none';
      }};

      svg.appendChild(circle);
    }});
  }}

  function openZoomModal(src) {{
    const modal = document.getElementById('zoom-modal');
    const modalImg = document.getElementById('modal-img');
    modalImg.src = src;
    modal.classList.add('open');
  }}

  function closeZoomModal() {{
    document.getElementById('zoom-modal').classList.remove('open');
  }}

  function exportDataTSV() {{
    let tsv = 'PATIENT_ID\\tSAMPLE_ID\\tTIMEPOINT\\tSAMPLING_DATE\\tWISP_PURITY_PCT\\tWISP_LOD_PCT\\tIS_MRD_POSITIVE\\tSNV_PURITY_PCT\\tCNV_PURITY_PCT\\tLOH_PURITY_PCT\\tICHOR_TF\\n';
    PLASMAS.forEach(p => {{
      tsv += `${{p.patient_id}}\\t${{p.sample_id}}\\t${{p.timepoint}}\\t${{p.sampling_date}}\\t${{p.wisp_purity_pct}}\\t${{p.wisp_lod_pct}}\\t${{p.is_mrd_positive}}\\t${{p.snv_purity_pct}}\\t${{p.cnv_purity_pct}}\\t${{p.loh_purity_pct}}\\t${{p.ichor_tf || 'NA'}}\\n`;
    }});
    const blob = new Blob([tsv], {{ type: 'text/tab-separated-values' }});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${{PLASMAS[0]?.patient_id || 'patient'}}_wisp_ctdna_summary.tsv`;
    a.click();
    URL.revokeObjectURL(url);
  }}

  window.addEventListener('resize', renderChart);
  window.addEventListener('DOMContentLoaded', init);
</script>

</body>
</html>
"""
    out_html_path = Path(out_html_path)
    out_html_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"  Generated interactive WISP patient view: {out_html_path}")


def write_cbio_resource_files(import_dir, study_id, patients_data, relative_url_prefix="reports/wisp/"):
    """
    Writes cBioPortal Resource Definition and Mapping files:
      - meta_resource_definition.txt
      - data_resource_definition.txt
      - meta_resource_patient.txt
      - data_resource_patient.txt
      - meta_resource_sample.txt
      - data_resource_sample.txt
    """
    import_path = Path(import_dir)
    import_path.mkdir(parents=True, exist_ok=True)

    has_wisp = any(len(pdata.get("plasmas", [])) > 0 for pdata in patients_data.values())
    has_orange = any(len(pdata.get("orange_reports", [])) > 0 for pdata in patients_data.values())

    # Stage Orange reports
    if has_orange:
        orange_out_dir = import_path / "reports" / "orange"
        orange_out_dir.mkdir(parents=True, exist_ok=True)
        import shutil
        for pid, pdata in patients_data.items():
            for o_item in pdata.get("orange_reports", []):
                src_path = Path(o_item["path"])
                dest_file = orange_out_dir / f"{o_item['sample_id']}_orange{o_item['ext']}"
                if src_path.is_file() and not dest_file.exists():
                    try:
                        shutil.copy2(src_path, dest_file)
                    except Exception as e:
                        print(f"Warning: Failed staging orange report {src_path}: {e}", file=sys.stderr)
                o_item["staged_rel"] = f"reports/orange/{dest_file.name}"

    # 1. meta_resource_definition.txt
    meta_def_file = import_path / "meta_resource_definition.txt"
    with open(meta_def_file, "w") as f:
        f.write(f"cancer_study_identifier: {study_id}\n")
        f.write("resource_type: DEFINITION\n")
        f.write("data_filename: data_resource_definition.txt\n")

    # 2. data_resource_definition.txt
    # RESOURCE_ID, DISPLAY_NAME, RESOURCE_TYPE, DESCRIPTION, OPEN_BY_DEFAULT, PRIORITY
    data_def_file = import_path / "data_resource_definition.txt"
    with open(data_def_file, "w") as f:
        f.write("RESOURCE_ID\tDISPLAY_NAME\tRESOURCE_TYPE\tDESCRIPTION\tOPEN_BY_DEFAULT\tPRIORITY\n")
        if has_wisp:
            f.write("WISP_PLOTS\tLiquid Biopsy (WISP)\tPATIENT\tPlasma ctDNA multi-modal WISP decomposition and longitudinal MRD monitoring\tFALSE\t1\n")
        if has_orange:
            f.write("ORANGE_REPORT\tOrange Report\tSAMPLE\tComprehensive HMFTools Orange whole-genome clinical interpretation report\tFALSE\t2\n")

    # 3. meta_resource_patient.txt
    meta_pat_file = import_path / "meta_resource_patient.txt"
    with open(meta_pat_file, "w") as f:
        f.write(f"cancer_study_identifier: {study_id}\n")
        f.write("resource_type: PATIENT\n")
        f.write("data_filename: data_resource_patient.txt\n")

    # 4. data_resource_patient.txt
    # PATIENT_ID, RESOURCE_ID, URL
    data_pat_file = import_path / "data_resource_patient.txt"
    with open(data_pat_file, "w") as f:
        f.write("PATIENT_ID\tRESOURCE_ID\tURL\n")
        for pid, pdata in sorted(patients_data.items()):
            if pdata.get("plasmas"):
                rel_url = f"{relative_url_prefix}{pid}_wisp.html"
                f.write(f"{pid}\tWISP_PLOTS\t{rel_url}\n")
            if pdata.get("orange_reports"):
                first_orange = pdata["orange_reports"][0]
                o_url = first_orange.get("staged_rel") or f"reports/orange/{first_orange['sample_id']}_orange{first_orange['ext']}"
                f.write(f"{pid}\tORANGE_REPORT\t{o_url}\n")

    # 5. meta_resource_sample.txt
    meta_smp_file = import_path / "meta_resource_sample.txt"
    with open(meta_smp_file, "w") as f:
        f.write(f"cancer_study_identifier: {study_id}\n")
        f.write("resource_type: SAMPLE\n")
        f.write("data_filename: data_resource_sample.txt\n")

    # 6. data_resource_sample.txt
    # PATIENT_ID, SAMPLE_ID, RESOURCE_ID, URL
    data_smp_file = import_path / "data_resource_sample.txt"
    with open(data_smp_file, "w") as f:
        f.write("PATIENT_ID\tSAMPLE_ID\tRESOURCE_ID\tURL\n")
        for pid, pdata in sorted(patients_data.items()):
            for pl in pdata.get("plasmas", []):
                sid = pl["sample_id"]
                rel_url = f"{relative_url_prefix}{pid}_wisp.html#{sid}"
                f.write(f"{pid}\t{sid}\tWISP_PLOTS\t{rel_url}\n")
            for o_item in pdata.get("orange_reports", []):
                sid = o_item["sample_id"]
                o_url = o_item.get("staged_rel") or f"reports/orange/{sid}_orange{o_item['ext']}"
                f.write(f"{pid}\t{sid}\tORANGE_REPORT\t{o_url}\n")

    print(f"  Wrote cBioPortal resource definition and mapping files to: {import_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate cBioPortal WISP Liquid Biopsy Resource Panel & Interactive HTML Views")
    parser.add_argument("-i", "--in-dir", required=True, help="Input directory containing Oncoanalyser / WISP outputs")
    parser.add_argument("-o", "--import-dir", required=True, help="cBioPortal study import directory (where data_resource_*.txt will be placed)")
    parser.add_argument("-s", "--study-id", required=True, help="cBioPortal Cancer Study Identifier")
    parser.add_argument("--relative-url-prefix", default="reports/wisp/", help="Relative URL prefix for HTML reports (default: reports/wisp/)")
    parser.add_argument("--skip-plots", action="store_true", help="Skip matplotlib high-resolution image generation")

    args = parser.parse_args()

    in_dir = Path(args.in_dir).resolve()
    import_dir = Path(args.import_dir).resolve()
    study_id = args.study_id

    print(f"[WISP Resource] Scanning {in_dir} for plasma WGS and WISP files...")
    patients_data = discover_wisp_data(in_dir, import_dir=import_dir)

    has_wisp = any(len(p.get("plasmas", [])) > 0 for p in patients_data.values()) if patients_data else False
    has_orange = any(len(p.get("orange_reports", [])) > 0 for p in patients_data.values()) if patients_data else False

    if not has_wisp and not has_orange:
        print(f"[Resource Generator] No plasma WGS / WISP files or Orange reports discovered in {in_dir}. Skipping resource panel generation.")
        sys.exit(0)

    total_patients = len(patients_data)
    total_plasmas = sum(len(p.get("plasmas", [])) for p in patients_data.values())
    total_orange = sum(len(p.get("orange_reports", [])) for p in patients_data.values())
    print(f"[Resource Generator] Discovered {total_plasmas} plasma samples (WISP) and {total_orange} Orange reports across {total_patients} patients.")

    reports_dir = import_dir / "reports" / "wisp"
    plots_dir = reports_dir / "plots"
    if has_wisp:
        reports_dir.mkdir(parents=True, exist_ok=True)
        plots_dir.mkdir(parents=True, exist_ok=True)

    # Generate multi-modal figures and interactive HTML for each patient with plasma data
    for pid, pdata in patients_data.items():
        plasmas = pdata.get("plasmas", [])
        if not plasmas:
            continue
        tumors = pdata.get("tumors", [])
        tumor_rec = tumors[0] if tumors else None

        for pl in plasmas:
            sid = pl["sample_id"]
            plot_fname = f"{pid}_{sid}_wisp_multimodal.png"
            out_img = plots_dir / plot_fname

            if not args.skip_plots and HAS_CWS:
                try:
                    ok = cws.generate_plasma_wisp_multimodal_plot(pl, tumor_rec, str(out_img))
                    if ok and out_img.is_file():
                        pl["multimodal_plot"] = f"plots/{plot_fname}"
                        try:
                            with open(out_img, "rb") as f_img:
                                pl["multimodal_plot_b64"] = f"data:image/png;base64,{base64.b64encode(f_img.read()).decode('ascii')}"
                        except Exception:
                            pass
                except Exception as e:
                    print(f"Warning: Failed generating multi-modal plot for {pid} {sid}: {e}", file=sys.stderr)

        # Generate interactive patient HTML
        html_out = reports_dir / f"{pid}_wisp.html"
        generate_interactive_patient_wisp_html(pid, pdata, html_out)

    # Write cBioPortal Resource Manifests
    write_cbio_resource_files(import_dir, study_id, patients_data, relative_url_prefix=args.relative_url_prefix)
    print(f"[WISP Resource] Successfully created WISP panel resource definition for cBioPortal study '{study_id}'.")


if __name__ == "__main__":
    main()
