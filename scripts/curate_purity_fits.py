#!/usr/bin/env python3
"""
curate_purity_fits.py - Interactive Purity & Ploidy Curation Workflow

Generates an interactive HTML console and helper CLI to systematically review:
  1. PURPLE 2D Purity/Ploidy contour plots (*.purity.range.png)
  2. High-quality somatic driver variant VAFs & implied clonal purity (2 * VAF)
  3. Alternative local minima valleys from *.purple.purity.range.tsv
  4. Native HMF ORANGE clinical reports (*.orange.pdf)
  5. Native sWGS ichorCNA multi-solution genome-wide CNV plots (*.params.txt, *_genomeWide_*.pdf/png)
  6. Exports/manages a standard purity_overrides.tsv consumed by compile_wgts_summary.py
"""

import os
import sys
import json
import csv
import argparse
import base64
from pathlib import Path


def safe_float(v):
    if v is None:
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def safe_int(v):
    if v is None:
        return None
    try:
        return int(float(v))
    except (ValueError, TypeError):
        return None


def parse_purple_purity_tsv(path):
    if not path or not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                return {
                    "purity": safe_float(r.get("purity")),
                    "ploidy": safe_float(r.get("ploidy")),
                    "status": (r.get("status") or "").strip(),
                    "score": safe_float(r.get("score")),
                    "wgd": str(r.get("wholeGenomeDuplication", "")).lower() == "true",
                    "tmb": safe_float(r.get("tmbPerMb")),
                }
    except Exception:
        pass
    return {}


def parse_purity_range_tsv(path, top_n=6):
    if not path or not os.path.isfile(path):
        return []
    rows = []
    try:
        with open(path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                sc = safe_float(r.get("score"))
                pur = safe_float(r.get("purity"))
                pl = safe_float(r.get("ploidy"))
                sp = safe_float(r.get("somaticPenalty", 0.0))
                dip = safe_float(r.get("diploidProportion", 0.0))
                if sc is not None and pur is not None and pl is not None:
                    rows.append((sc, pur, pl, sp, dip))
        rows.sort(key=lambda x: x[0])
        
        candidates = []
        seen = set()
        for sc, pur, pl, sp, dip in rows:
            # Cluster by purity and ploidy bins to return distinct valleys
            key = (round(pur, 2), round(pl, 1))
            if key not in seen:
                seen.add(key)
                candidates.append({
                    "score": sc,
                    "purity": pur,
                    "ploidy": pl,
                    "somatic_penalty": sp,
                    "diploid_prop": dip
                })
                if len(candidates) >= top_n:
                    break
        return candidates
    except Exception:
        return []


def parse_ichorcna_params(path):
    """
    Parses ichorCNA *.params.txt file.
    Extracts top selected tumor fraction, ploidy, MAD, and all candidate solutions.
    """
    res = {
        "optimal_tumor_fraction": None,
        "optimal_ploidy": None,
        "mad": None,
        "gender": None,
        "solutions": []
    }
    if not path or not os.path.isfile(path):
        return res
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = [l.strip() for l in f if l.strip()]
        
        in_table = False
        headers = []
        for line in lines:
            if line.startswith("Gender:"):
                res["gender"] = line.split(":", 1)[1].strip()
            elif line.startswith("Tumor Fraction:"):
                try:
                    res["optimal_tumor_fraction"] = float(line.split(":", 1)[1].strip())
                except Exception:
                    pass
            elif line.startswith("Ploidy:"):
                try:
                    res["optimal_ploidy"] = float(line.split(":", 1)[1].strip())
                except Exception:
                    pass
            elif line.startswith("GC-Map correction MAD:"):
                try:
                    res["mad"] = float(line.split(":", 1)[1].strip())
                except Exception:
                    pass
            elif line.startswith("init\t") or line.startswith("init "):
                in_table = True
                headers = [h.strip() for h in line.split()]
            elif in_table:
                parts = line.split()
                if len(parts) >= len(headers):
                    row = dict(zip(headers, parts))
                    init_id = row.get("init")
                    n_est = None
                    phi = None
                    loglik = None
                    subc = None
                    try:
                        n_est = float(row.get("n_est"))
                    except Exception:
                        pass
                    try:
                        phi = float(row.get("phi_est"))
                    except Exception:
                        pass
                    try:
                        loglik = float(row.get("loglik"))
                    except Exception:
                        pass
                    try:
                        subc = float(row.get("Frac_genome_subclonal"))
                    except Exception:
                        pass
                    
                    tf = (1.0 - n_est) if n_est is not None else None
                    res["solutions"].append({
                        "init": init_id,
                        "normal_est": n_est,
                        "tumor_fraction": tf,
                        "ploidy": phi,
                        "loglik": loglik,
                        "subclonal_frac": subc,
                        "is_optimal": False
                    })
        if res["solutions"]:
            # Sort descending by log-likelihood
            res["solutions"].sort(key=lambda s: -(s["loglik"] if s["loglik"] is not None else -999999))
            for i, s in enumerate(res["solutions"]):
                s["is_optimal"] = (i == 0)
    except Exception:
        pass
    return res


def discover_ichor_data(pat_dir, tumor_id, out_dir=None, embed_images=False):
    """
    Discovers all ichorCNA metrics, params, and plots for a given tumor.
    """
    ichor_data = {
        "has_ichor": False,
        "optimal_tumor_fraction": None,
        "optimal_ploidy": None,
        "mad": None,
        "gender": None,
        "all_sols_pdf": "",
        "all_sols_rel": "",
        "optimal_pdf": "",
        "optimal_pdf_rel": "",
        "solutions": []
    }
    
    # Locate ichorCNA directory
    cand_dirs = [
        pat_dir / "ichorCNA_t_vs_n" / tumor_id,
        pat_dir / "ichorCNA_plasma_vs_n" / tumor_id,
        pat_dir / "ichorCNA_tumor_only" / tumor_id,
        pat_dir / "ichorCNA_plasma_only" / tumor_id,
        pat_dir / "ichorCNA" / tumor_id,
        pat_dir / "ichorCNA_normal_only" / tumor_id,
        pat_dir / "plasma" / tumor_id / "ichorCNA" / tumor_id,
    ]
    ichor_dir = None
    for cd in cand_dirs:
        if cd.is_dir():
            ichor_dir = cd
            break
            
    if not ichor_dir:
        # Check any subdirectory inside known ichor directory types
        for b_name in ["ichorCNA_t_vs_n", "ichorCNA_plasma_vs_n", "ichorCNA_tumor_only", "ichorCNA_plasma_only", "ichorCNA"]:
            base_ichor = pat_dir / b_name
            if base_ichor.is_dir():
                sub_dirs = [d for d in base_ichor.iterdir() if d.is_dir() and not d.name.startswith(".")]
                if len(sub_dirs) == 1:
                    ichor_dir = sub_dirs[0]
                    break

    if not ichor_dir or not ichor_dir.is_dir():
        return ichor_data

    params_file = ichor_dir / f"{tumor_id}.params.txt"
    if not params_file.is_file():
        cand_p = list(ichor_dir.glob("*.params.txt"))
        if cand_p:
            params_file = cand_p[0]

    p_parsed = parse_ichorcna_params(params_file)
    ichor_data["has_ichor"] = True
    ichor_data["optimal_tumor_fraction"] = p_parsed.get("optimal_tumor_fraction")
    ichor_data["optimal_ploidy"] = p_parsed.get("optimal_ploidy")
    ichor_data["mad"] = p_parsed.get("mad")
    ichor_data["gender"] = p_parsed.get("gender")

    nested = ichor_dir / tumor_id
    if not nested.is_dir():
        # Check if plots are in ichor_dir itself or subfolder
        sub_nested = [d for d in ichor_dir.iterdir() if d.is_dir() and not d.name.startswith(".")]
        nested = sub_nested[0] if sub_nested else ichor_dir

    all_sols_pdf = nested / f"{tumor_id}_genomeWide_all_sols.pdf"
    if not all_sols_pdf.is_file():
        cand = list(nested.glob("*genomeWide_all_sols.pdf"))
        if cand:
            all_sols_pdf = cand[0]

    top_pdf = nested / f"{tumor_id}_genomeWide.pdf"
    if not top_pdf.is_file():
        cand = list(nested.glob("*genomeWide.pdf"))
        if cand:
            top_pdf = cand[0]

    if all_sols_pdf.is_file():
        ichor_data["all_sols_pdf"] = str(all_sols_pdf)
        ichor_data["all_sols_rel"] = os.path.relpath(all_sols_pdf, out_dir) if out_dir else str(all_sols_pdf)

    if top_pdf.is_file():
        ichor_data["optimal_pdf"] = str(top_pdf)
        ichor_data["optimal_pdf_rel"] = os.path.relpath(top_pdf, out_dir) if out_dir else str(top_pdf)

    # Attach plots for each candidate solution
    solutions = []
    for sol in p_parsed.get("solutions", []):
        init_id = sol["init"]
        sol_pdf = nested / f"{tumor_id}_genomeWide_{init_id}.pdf"
        sol_png = nested / f"{tumor_id}_genomeWide_{init_id}.png"
        
        pdf_path = str(sol_pdf) if sol_pdf.is_file() else ""
        pdf_rel = os.path.relpath(sol_pdf, out_dir) if (out_dir and sol_pdf.is_file()) else pdf_path
        png_path = str(sol_png) if sol_png.is_file() else ""
        png_rel = os.path.relpath(sol_png, out_dir) if (out_dir and sol_png.is_file()) else png_path
        
        png_b64 = ""
        if embed_images and sol_png.is_file():
            try:
                with open(sol_png, "rb") as pf:
                    png_b64 = "data:image/png;base64," + base64.b64encode(pf.read()).decode("ascii")
            except Exception:
                pass

        solutions.append({
            "init": init_id,
            "normal_est": sol["normal_est"],
            "tumor_fraction": sol["tumor_fraction"],
            "ploidy": sol["ploidy"],
            "loglik": sol["loglik"],
            "subclonal_frac": sol["subclonal_frac"],
            "is_optimal": sol["is_optimal"],
            "pdf_path": pdf_path,
            "pdf_rel": pdf_rel,
            "png_path": png_path,
            "png_rel": png_rel,
            "png_base64": png_b64
        })

    ichor_data["solutions"] = solutions
    return ichor_data


def discover_orange_report(pat_dir, tumor_id, p_dir=None, out_dir=None):
    """
    Locates native HMF Orange clinical PDF and JSON reports.
    """
    orange_data = {
        "has_orange": False,
        "json_path": "",
        "pdf_path": "",
        "pdf_rel": ""
    }
    
    cand_dirs = []
    if p_dir:
        cand_dirs.append(p_dir.parent / "orange")
    cand_dirs.append(pat_dir / tumor_id / "orange")
    cand_dirs.append(pat_dir / "orange")
    
    for od in cand_dirs:
        if od.is_dir():
            pdf_file = od / f"{tumor_id}.orange.pdf"
            if not pdf_file.is_file():
                cand_pdfs = list(od.glob("*.orange.pdf"))
                if cand_pdfs:
                    pdf_file = cand_pdfs[0]
            
            json_file = od / f"{tumor_id}.orange.json"
            if not json_file.is_file():
                cand_jsons = list(od.glob("*.orange.json"))
                if cand_jsons:
                    json_file = cand_jsons[0]
                    
            if pdf_file.is_file():
                orange_data["has_orange"] = True
                orange_data["pdf_path"] = str(pdf_file)
                orange_data["pdf_rel"] = os.path.relpath(pdf_file, out_dir) if out_dir else str(pdf_file)
            if json_file.is_file():
                orange_data["json_path"] = str(json_file)
            if orange_data["has_orange"]:
                break
                
    return orange_data


def parse_drivers(orange_path, cat_somatic_path):
    drivers = []
    seen = set()
    
    # 1. Orange JSON (primary source for detailed depth and clinical annotations)
    if orange_path and os.path.isfile(orange_path):
        try:
            with open(orange_path, "r", encoding="utf-8-sig") as f:
                oj = json.load(f)
                p_data = oj.get("purple", {})
                
                # Somatic variants
                for v in p_data.get("somaticVariants", []):
                    g = (v.get("gene") or "").strip().upper()
                    ci = v.get("canonicalImpact") or {}
                    p = (ci.get("hgvsProteinImpact") or v.get("canonicalHgvsProtein") or "").strip()
                    c = (ci.get("hgvsCodingImpact") or v.get("canonicalHgvsCoding") or "").strip()
                    td = v.get("tumorDepth") or {}
                    nd = v.get("normalDepth") or {}
                    tot = safe_int(td.get("totalReadCount") or v.get("totalReadCount"))
                    alt = safe_int(td.get("alleleReadCount") or v.get("alleleReadCount"))
                    norm_tot = safe_int(nd.get("totalReadCount"))
                    norm_alt = safe_int(nd.get("alleleReadCount"))
                    vaf = (alt / tot) if (tot and alt is not None and tot > 0) else safe_float(v.get("adjustedVAF"))
                    
                    # Flag potential germline / zero-normal variants
                    is_zero_normal = (norm_tot is not None and norm_tot == 0) or (norm_alt is not None and norm_alt > 0)
                    
                    key = (g, p, c)
                    if g and key not in seen:
                        seen.add(key)
                        drivers.append({
                            "gene": g,
                            "protein": p,
                            "coding": c,
                            "alt": alt,
                            "total": tot,
                            "norm_alt": norm_alt,
                            "norm_total": norm_tot,
                            "is_zero_normal": is_zero_normal,
                            "vaf": vaf,
                            "implied_purity": min(1.0, 2.0 * vaf) if vaf is not None else None,
                            "origin": "Somatic",
                            "driver_type": v.get("worstCodingEffect", "Mutation"),
                            "is_hotspot": v.get("hotspot") == "HOTSPOT"
                        })
                        
                # Germline variants
                for v in p_data.get("germlineVariants", []):
                    g = (v.get("gene") or "").strip().upper()
                    ci = v.get("canonicalImpact") or {}
                    p = (ci.get("hgvsProteinImpact") or v.get("canonicalHgvsProtein") or "").strip()
                    c = (ci.get("hgvsCodingImpact") or v.get("canonicalHgvsCoding") or "").strip()
                    td = v.get("tumorDepth") or {}
                    tot = safe_int(td.get("totalReadCount") or v.get("totalReadCount"))
                    alt = safe_int(td.get("alleleReadCount") or v.get("alleleReadCount"))
                    vaf = (alt / tot) if (tot and alt is not None and tot > 0) else None
                    
                    key = (g, p, c)
                    if g and key not in seen:
                        seen.add(key)
                        drivers.append({
                            "gene": g,
                            "protein": p,
                            "coding": c,
                            "alt": alt,
                            "total": tot,
                            "norm_alt": None,
                            "norm_total": None,
                            "is_zero_normal": False,
                            "vaf": vaf,
                            "implied_purity": None,
                            "origin": "Germline",
                            "driver_type": v.get("worstCodingEffect", "Germline"),
                            "is_hotspot": False
                        })
        except Exception:
            pass

    # 2. Driver Catalog TSV fallback
    if cat_somatic_path and os.path.isfile(cat_somatic_path):
        try:
            with open(cat_somatic_path, "r", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    g = (r.get("gene") or "").strip().upper()
                    p = (r.get("canonicalHgvsProtein") or r.get("hgvsProtein") or "").strip()
                    c = (r.get("canonicalHgvsCoding") or r.get("hgvsCoding") or "").strip()
                    vaf = safe_float(r.get("vaf"))
                    key = (g, p, c)
                    if g and key not in seen:
                        seen.add(key)
                        drivers.append({
                            "gene": g,
                            "protein": p,
                            "coding": c,
                            "alt": None,
                            "total": None,
                            "norm_alt": None,
                            "norm_total": None,
                            "is_zero_normal": False,
                            "vaf": vaf,
                            "implied_purity": min(1.0, 2.0 * vaf) if vaf is not None else None,
                            "origin": "Somatic",
                            "driver_type": r.get("driver", "Mutation"),
                            "is_hotspot": "HOTSPOT" in str(r.get("driver", "")).upper()
                        })
        except Exception:
            pass

    # Sort somatic mutations by VAF descending (clonal drivers on top)
    drivers.sort(key=lambda x: (x["origin"] != "Somatic", -(x["vaf"] or 0.0)))
    return drivers


def build_curation_dashboard(cases, output_html, existing_overrides=None):
    if existing_overrides is None:
        existing_overrides = {}

    cases_json = json.dumps(cases)
    overrides_json = json.dumps(existing_overrides)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>WGTS Tumor Purity & Ploidy Curation Console</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Outfit:wght@600;700;800&display=swap" rel="stylesheet">
  <style>
    :root {{
      --bg: #0d1117;
      --card-bg: #161b22;
      --card-inner: #1f242c;
      --border: #30363d;
      --text: #c9d1d9;
      --text-muted: #8b949e;
      --accent: #58a6ff;
      --accent-green: #238636;
      --accent-amber: #d29922;
      --accent-purple: #a855f7;
      --accent-orange: #f59e0b;
      --accent-red: #f85149;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      background: var(--bg);
      color: var(--text);
      margin: 0;
      padding: 1.5rem;
    }}
    .header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      border-bottom: 1px solid var(--border);
      padding-bottom: 1rem;
      margin-bottom: 1.25rem;
    }}
    h1 {{
      font-family: 'Outfit', sans-serif;
      margin: 0;
      font-size: 1.6rem;
      font-weight: 800;
      color: #fff;
      letter-spacing: -0.02em;
    }}
    .subtitle {{ color: var(--text-muted); font-size: 0.85rem; margin-top: 0.25rem; }}
    .actions {{ display: flex; gap: 0.75rem; }}
    .btn {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      color: var(--text);
      padding: 0.45rem 0.9rem;
      border-radius: 6px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 5px;
      text-decoration: none;
      transition: all 0.15s ease;
    }}
    .btn:hover {{ background: #21262d; border-color: #8b949e; color: #fff; }}
    .btn-primary {{ background: var(--accent-green); color: #fff; border-color: rgba(240,246,252,0.1); }}
    .btn-primary:hover {{ background: #2ea043; }}
    .btn-orange {{ background: #92400e; color: #fef3c7; border-color: #d97706; }}
    .btn-orange:hover {{ background: #b45309; color: #fff; }}
    .btn-purple {{ background: #581c87; color: #f3e8ff; border-color: #9333ea; }}
    .btn-purple:hover {{ background: #6b21a8; color: #fff; }}
    
    .badge {{
      display: inline-block;
      padding: 0.2rem 0.5rem;
      border-radius: 12px;
      font-size: 0.72rem;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }}
    .badge-normal {{ background: rgba(210, 153, 34, 0.15); color: #e3b341; border: 1px solid rgba(210, 153, 34, 0.4); }}
    .badge-somatic {{ background: rgba(35, 134, 54, 0.2); color: #3fb950; border: 1px solid rgba(35, 134, 54, 0.4); }}
    .badge-curated {{ background: rgba(88, 166, 255, 0.2); color: #79c0ff; border: 1px solid rgba(88, 166, 255, 0.4); }}
    .badge-optimal {{ background: rgba(35, 134, 54, 0.25); color: #4ade80; border: 1px solid rgba(74, 222, 128, 0.4); }}
    .badge-warn {{ background: rgba(239, 68, 68, 0.2); color: #f87171; border: 1px solid rgba(239, 68, 68, 0.4); }}
    
    .grid-container {{
      display: grid;
      grid-template-columns: 290px 1fr;
      gap: 1.25rem;
      min-height: calc(100vh - 120px);
    }}
    .sample-list {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 0.5rem;
      overflow-y: auto;
      max-height: calc(100vh - 130px);
    }}
    .sample-item {{
      padding: 0.75rem;
      border-radius: 6px;
      cursor: pointer;
      border: 1px solid transparent;
      margin-bottom: 0.35rem;
      transition: all 0.15s ease;
    }}
    .sample-item:hover {{ background: #21262d; }}
    .sample-item.active {{
      background: #1f2937;
      border-color: var(--accent);
    }}
    .sample-item-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 0.25rem;
    }}
    .sample-item-id {{ font-weight: 700; color: #fff; font-size: 0.9rem; }}
    .sample-item-metrics {{ font-size: 0.78rem; color: var(--text-muted); }}

    .main-panel {{
      display: flex;
      flex-direction: column;
      gap: 1.25rem;
    }}
    .card {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 1.25rem;
    }}
    .card-title {{
      font-size: 1rem;
      font-weight: 700;
      color: #fff;
      margin: 0 0 0.75rem 0;
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 0.5rem;
    }}
    
    .two-col {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 1.25rem;
    }}
    
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.8rem;
    }}
    th, td {{
      padding: 0.5rem 0.6rem;
      text-align: left;
      border-bottom: 1px solid #21262d;
    }}
    th {{
      color: var(--text-muted);
      font-weight: 600;
      background: rgba(0,0,0,0.2);
    }}
    tr:hover td {{ background: rgba(255,255,255,0.02); }}
    
    /* Tabs for Visual Plots */
    .tabs-nav {{
      display: flex;
      gap: 0.5rem;
      border-bottom: 1px solid var(--border);
      padding-bottom: 0.6rem;
      margin-bottom: 0.75rem;
      align-items: center;
    }}
    .tab-btn {{
      background: var(--card-inner);
      border: 1px solid var(--border);
      color: var(--text-muted);
      padding: 0.4rem 0.85rem;
      border-radius: 6px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.15s ease;
    }}
    .tab-btn:hover {{ color: #fff; border-color: var(--text-muted); }}
    .tab-btn.active {{
      background: #1f2937;
      color: #fff;
      border-color: var(--accent);
    }}
    
    .tab-pane {{
      display: none;
    }}
    .tab-pane.active {{
      display: block;
    }}
    
    .contour-img-wrap {{
      text-align: center;
      background: #000;
      padding: 0.5rem;
      border-radius: 6px;
      border: 1px solid var(--border);
      min-height: 480px;
      display: flex;
      align-items: center;
      justify-content: center;
    }}
    .contour-img {{
      max-width: 100%;
      max-height: 540px;
      height: auto;
      border-radius: 4px;
      display: block;
      margin: 0 auto;
    }}
    
    /* ichorCNA Controls & Viewer */
    .ichor-controls {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 0.75rem;
      margin-bottom: 0.75rem;
      flex-wrap: wrap;
    }}
    .ichor-viewer-wrap {{
      border: 1px solid var(--border);
      border-radius: 6px;
      overflow: hidden;
      background: #ffffff;
      min-height: 520px;
      display: flex;
      flex-direction: column;
    }}
    .ichor-iframe {{
      width: 100%;
      height: 550px;
      border: none;
      background: #fff;
    }}
    .ichor-png {{
      width: 100%;
      max-height: 550px;
      object-fit: contain;
      background: #fff;
      display: block;
    }}
    
    .curation-box {{
      background: #1c2128;
      border: 1px solid #388bfd33;
      border-radius: 8px;
      padding: 1rem;
      display: flex;
      align-items: center;
      gap: 1.25rem;
      margin-top: 0.5rem;
      flex-wrap: wrap;
    }}
    .input-group {{
      display: flex;
      flex-direction: column;
      gap: 0.25rem;
    }}
    .input-group label {{
      font-size: 0.72rem;
      font-weight: 700;
      color: var(--text-muted);
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}
    .input-group input, .input-group select {{
      background: #0d1117;
      border: 1px solid var(--border);
      color: #fff;
      padding: 0.45rem 0.65rem;
      border-radius: 6px;
      font-size: 0.88rem;
      min-width: 110px;
    }}
    .input-group input:focus, .input-group select:focus {{
      border-color: var(--accent);
      outline: none;
    }}
    .input-group.notes {{
      flex: 1;
      min-width: 240px;
    }}
    .input-group.notes input {{
      width: 100%;
    }}
  </style>
</head>
<body>

  <div class="header">
    <div>
      <h1>WGTS Tumor Purity & Ploidy Curation Console</h1>
      <div class="subtitle">Systematic multi-tier review: PURPLE 2D contour penalty maps, sWGS ichorCNA candidate solutions, and driver variant anchors</div>
    </div>
    <div class="actions">
      <button class="btn" onclick="exportOverridesTSV()">📥 Export overrides.tsv</button>
      <button class="btn btn-primary" onclick="copyOverridesTSV()">📋 Copy to Clipboard</button>
    </div>
  </div>

  <div class="grid-container">
    <div class="sample-list" id="sampleList"></div>

    <div class="main-panel">
      <!-- Active Case Header & Curation Controls -->
      <div class="card">
        <div class="card-title">
          <div>
            <span id="caseHeaderTitle" style="font-size:1.15rem; color:#fff;">Select a patient</span>
            <span id="caseStatusBadge" style="margin-left: 0.5rem;"></span>
          </div>
          <div style="display: flex; gap: 0.5rem; align-items: center;">
            <a id="btnOrangeReport" href="#" target="_blank" class="btn btn-orange" style="display: none;" title="Open native HMF Orange Clinical PDF Report">
              📄 View ORANGE Report (PDF) ↗
            </a>
            <a id="btnIchorReport" href="#" target="_blank" class="btn btn-purple" style="display: none;" title="Open ichorCNA all-solutions genome-wide CNV plot">
              📊 ichorCNA All Solutions (PDF) ↗
            </a>
          </div>
        </div>
        
        <div class="curation-box">
          <div class="input-group">
            <label>Curated Purity (%)</label>
            <input type="number" id="inpPurity" step="0.5" min="1" max="100" onchange="updateActiveCase()">
          </div>
          <div class="input-group">
            <label>Curated Ploidy</label>
            <input type="number" id="inpPloidy" step="0.05" min="1" max="8" onchange="updateActiveCase()">
          </div>
          <div class="input-group">
            <label>Curation Method</label>
            <select id="selMethod" onchange="updateActiveCase()">
              <option value="sWGS ichorCNA concordant">sWGS ichorCNA concordant</option>
              <option value="Low Cellularity Baseline (Sub-threshold CNV)">Low Cellularity Baseline (Sub-threshold CNV)</option>
              <option value="PURPLE 2D contour valley">PURPLE 2D contour valley</option>
              <option value="Clonal VAF anchor">Clonal VAF anchor</option>
              <option value="Pathology cellularity">Pathology cellularity</option>
              <option value="Manual estimate">Manual estimate</option>
            </select>
          </div>
          <div class="input-group notes">
            <label>Rationale & Notes</label>
            <input type="text" id="inpNotes" placeholder="e.g. sWGS ichorCNA 15.0% TF + copy-neutral quiet genome" oninput="updateActiveCase()">
          </div>
          <button class="btn btn-primary" style="align-self: flex-end;" onclick="saveAndNext()">Save & Next &rarr;</button>
        </div>
      </div>

      <!-- Two-column view: Visual Plot Tabs on left, Drivers & Minima on right -->
      <div class="two-col">
        <!-- Visual Plot Container with Tabs -->
        <div class="card">
          <div class="tabs-nav">
            <button id="tabBtnPurple" class="tab-btn active" onclick="switchVisualTab('purple')">PURPLE 2D Contour Map</button>
            <button id="tabBtnIchor" class="tab-btn" onclick="switchVisualTab('ichor')">ichorCNA Genome-Wide CNV</button>
            <div style="flex:1;"></div>
            <span id="tabVisualBadge" style="font-size:0.75rem; color:var(--text-muted);"></span>
          </div>

          <!-- Tab 1: PURPLE Contour -->
          <div id="tabPanePurple" class="tab-pane active">
            <div class="contour-img-wrap">
              <img id="contourPlot" class="contour-img" src="" alt="PURPLE 2D Contour Plot">
            </div>
            <div style="margin-top:0.5rem; font-size:0.75rem; color:var(--text-muted); display:flex; justify-content:space-between;">
              <span>Red cross indicates PURPLE selected minimum. Blue contours show penalty surface.</span>
              <span id="lblPurpleScore"></span>
            </div>
          </div>

          <!-- Tab 2: ichorCNA Solutions -->
          <div id="tabPaneIchor" class="tab-pane">
            <div class="ichor-controls">
              <div style="display:flex; align-items:center; gap:0.5rem; flex:1;">
                <label style="font-size:0.75rem; font-weight:700; color:var(--text-muted); text-transform:uppercase;">Solution:</label>
                <select id="selIchorSol" style="background:#0d1117; border:1px solid var(--border); color:#fff; padding:0.35rem 0.6rem; border-radius:6px; font-size:0.82rem; flex:1;" onchange="renderSelectedIchorSolution()">
                  <!-- Options filled dynamically -->
                </select>
              </div>
              <div style="display:flex; gap:0.4rem;">
                <button class="btn btn-purple" onclick="adoptActiveIchorSolution()" title="Adopt selected solution into curation form">
                  Adopt Solution
                </button>
                <a id="linkIchorSolPdf" target="_blank" href="#" class="btn" style="text-decoration:none;" title="Open active solution plot in new tab">
                  ↗ Pop-out
                </a>
              </div>
            </div>

            <!-- Comparison Table of all Evaluated ichorCNA Solutions -->
            <div style="max-height:165px; overflow-y:auto; border:1px solid var(--border); border-radius:6px; margin-bottom:0.75rem;">
              <table>
                <thead>
                  <tr>
                    <th>Init</th>
                    <th>Tumor Fraction</th>
                    <th>Ploidy</th>
                    <th>Log-Likelihood</th>
                    <th>Subclonal</th>
                    <th>Action</th>
                  </tr>
                </thead>
                <tbody id="ichorTableBody"></tbody>
              </table>
            </div>

            <!-- ichorCNA Plot Viewer -->
            <div class="ichor-viewer-wrap">
              <img id="ichorPlotImg" class="ichor-png" style="display:none;" src="" alt="ichorCNA Genome-Wide Plot">
              <iframe id="ichorPlotFrame" class="ichor-iframe" src=""></iframe>
            </div>
          </div>
        </div>

        <!-- Drivers & Alternative Valleys -->
        <div style="display: flex; flex-direction: column; gap: 1.25rem;">
          <!-- High-Quality Driver Variant Table -->
          <div class="card">
            <div class="card-title">
              <span>Reportable Driver Variants</span>
              <span style="font-size: 0.72rem; color: var(--text-muted);">Clonal VAF Anchors (Heuristic Reference)</span>
            </div>
            <table>
              <thead>
                <tr>
                  <th>Gene</th>
                  <th>Protein / Coding</th>
                  <th>Raw Depth (T / N)</th>
                  <th>Raw VAF</th>
                  <th>Implied Purity (2&times;VAF) *</th>
                  <th>Action</th>
                </tr>
              </thead>
              <tbody id="driverTableBody"></tbody>
            </table>
            <div style="font-size:0.72rem; color:var(--text-muted); margin-top:0.5rem; line-height:1.3;">
              * Note: Single-VAF implied purity assumes an unphased clonal heterozygous mutation in a copy-neutral region. If copy number is flat and BAF is balanced, whole-genome evidence and ichorCNA take precedence over single somatic/germline outliers.
            </div>
          </div>

          <!-- Top Alternative Local Minima from purity.range.tsv -->
          <div class="card">
            <div class="card-title">
              <span>Top PURPLE Valleys in Solution Grid</span>
              <span style="font-size: 0.72rem; color: var(--text-muted);">From purity.range.tsv</span>
            </div>
            <table>
              <thead>
                <tr>
                  <th>Purity</th>
                  <th>Ploidy</th>
                  <th>Score</th>
                  <th>Somatic Penalty</th>
                  <th>Diploid Prop</th>
                  <th>Action</th>
                </tr>
              </thead>
              <tbody id="candidatesTableBody"></tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  </div>

  <script>
    const CASES = {cases_json};
    const OVERRIDES = {overrides_json};
    let activeIdx = 0;
    let activeVisualTab = 'purple';

    function init() {{
      renderSampleList();
      selectCase(0);
    }}

    function renderSampleList() {{
      const listEl = document.getElementById('sampleList');
      listEl.innerHTML = '';
      CASES.forEach((c, idx) => {{
        const item = document.createElement('div');
        item.className = 'sample-item' + (idx === activeIdx ? ' active' : '');
        item.onclick = () => selectCase(idx);
        
        const ov = OVERRIDES[c.patient];
        const isCurated = !!ov;
        const curPur = isCurated ? ov.purity : (c.current_purity ? (c.current_purity * 100) : 0);
        const curPl = isCurated ? ov.ploidy : (c.current_ploidy || 2.0);
        
        let badgeHtml = '';
        if (isCurated) {{
          badgeHtml = '<span class="badge badge-curated">CURATED</span>';
        }} else if (c.status === 'NORMAL') {{
          badgeHtml = '<span class="badge badge-normal">LOW-CELL</span>';
        }} else {{
          badgeHtml = '<span class="badge badge-somatic">SOMATIC</span>';
        }}

        item.innerHTML = `
          <div class="sample-item-header">
            <span class="sample-item-id">${{c.patient}}</span>
            ${{badgeHtml}}
          </div>
          <div class="sample-item-metrics">
            Purity: <strong>${{curPur.toFixed(1)}}%</strong> · Ploidy: <strong>${{curPl.toFixed(2)}}</strong>
          </div>
        `;
        listEl.appendChild(item);
      }});
    }}

    function switchVisualTab(tab) {{
      activeVisualTab = tab;
      const btnPurple = document.getElementById('tabBtnPurple');
      const btnIchor = document.getElementById('tabBtnIchor');
      const panePurple = document.getElementById('tabPanePurple');
      const paneIchor = document.getElementById('tabPaneIchor');
      const lblBadge = document.getElementById('tabVisualBadge');

      if (tab === 'purple') {{
        btnPurple.classList.add('active');
        btnIchor.classList.remove('active');
        panePurple.classList.add('active');
        paneIchor.classList.remove('active');
        lblBadge.innerText = 'PURPLE Solution Surface';
      }} else {{
        btnPurple.classList.remove('active');
        btnIchor.classList.add('active');
        panePurple.classList.remove('active');
        paneIchor.classList.add('active');
        lblBadge.innerText = 'sWGS ichorCNA EM Log-Likelihood Models';
        renderSelectedIchorSolution();
      }}
    }}

    function selectCase(idx) {{
      activeIdx = idx;
      renderSampleList();
      const c = CASES[idx];
      
      document.getElementById('caseHeaderTitle').innerText = `${{c.patient}} (${{c.tumor}})`;
      
      const ov = OVERRIDES[c.patient];
      const isCurated = !!ov;
      
      let badgeHtml = '';
      if (isCurated) {{
        badgeHtml = '<span class="badge badge-curated">CURATED (' + ov.method + ')</span>';
      }} else if (c.status === 'NORMAL') {{
        badgeHtml = '<span class="badge badge-normal">LOW CELLULARITY BASELINE (' + c.status + ')</span>';
      }} else {{
        badgeHtml = '<span class="badge badge-somatic">SOMATIC FIT (' + c.status + ')</span>';
      }}
      document.getElementById('caseStatusBadge').innerHTML = badgeHtml;

      // Orange Report Link
      const btnOrange = document.getElementById('btnOrangeReport');
      if (c.orange && (c.orange.pdf_rel || c.orange.pdf_path)) {{
        btnOrange.href = c.orange.pdf_rel || c.orange.pdf_path;
        btnOrange.style.display = 'inline-flex';
      }} else {{
        btnOrange.style.display = 'none';
      }}

      // ichorCNA All Solutions Link
      const btnIchor = document.getElementById('btnIchorReport');
      if (c.ichor && (c.ichor.all_sols_rel || c.ichor.all_sols_pdf)) {{
        btnIchor.href = c.ichor.all_sols_rel || c.ichor.all_sols_pdf;
        btnIchor.style.display = 'inline-flex';
      }} else {{
        btnIchor.style.display = 'none';
      }}

      // Inputs
      const defaultPur = isCurated ? ov.purity : (c.current_purity ? (c.current_purity * 100) : 15.0);
      const defaultPl = isCurated ? ov.ploidy : (c.current_ploidy || 2.0);
      const defaultMeth = isCurated ? ov.method : (c.status === 'NORMAL' ? 'Low Cellularity Baseline (Sub-threshold CNV)' : 'sWGS ichorCNA concordant');
      const defaultNotes = isCurated ? ov.notes : '';

      document.getElementById('inpPurity').value = defaultPur.toFixed(1);
      document.getElementById('inpPloidy').value = defaultPl.toFixed(2);
      document.getElementById('selMethod').value = defaultMeth;
      document.getElementById('inpNotes').value = defaultNotes;

      // Tab 1: PURPLE Contour Image
      document.getElementById('contourPlot').src = c.plot_base64 || c.plot_path;
      document.getElementById('lblPurpleScore').innerText = c.candidates && c.candidates[0] ? `Optimal Score: ${{c.candidates[0].score.toFixed(4)}}` : '';

      // Tab 2: ichorCNA Setup
      setupIchorTab(c);

      // Drivers Table
      const dtBody = document.getElementById('driverTableBody');
      dtBody.innerHTML = '';
      if (c.drivers && c.drivers.length > 0) {{
        c.drivers.forEach(d => {{
          const tr = document.createElement('tr');
          const tDepth = (d.alt !== null && d.total) ? `${{d.alt}}/${{d.total}}` : 'NA';
          const nDepth = (d.norm_alt !== null && d.norm_total) ? `${{d.norm_alt}}/${{d.norm_total}}` : (d.is_zero_normal ? '<span class="badge badge-warn" style="font-size:0.65rem;">0x Normal</span>' : 'NA');
          const depthStr = `${{tDepth}} <span style="color:#8b949e; font-size:0.75rem;">(N: ${{nDepth}})</span>`;
          const vafStr = d.vaf !== null ? `${{(d.vaf * 100).toFixed(1)}}%` : 'NA';
          const impliedStr = d.implied_purity !== null ? `${{(d.implied_purity * 100).toFixed(1)}}%` : '—';
          
          let btnHtml = '';
          if (d.implied_purity !== null) {{
            const targetPur = (d.implied_purity * 100).toFixed(1);
            btnHtml = `<button class="btn" style="padding:0.2rem 0.5rem; font-size:0.72rem;" onclick="applyPreset(${{targetPur}}, null, 'Clonal VAF anchor (${{d.gene}} ${{vafStr}})')">Set Purity</button>`;
          }}
          
          let warnBadge = '';
          if (d.is_zero_normal || (d.vaf && d.vaf > 0.45 && d.vaf < 0.60)) {{
            warnBadge = ' <span class="badge badge-normal" style="font-size:0.62rem;" title="VAF ~50% with flat copy-number: verify whether germline before forcing high purity override">Likely Germline / Quiet CN</span>';
          }}
          
          tr.innerHTML = `
            <td><strong>${{d.gene}}</strong>${{warnBadge}}</td>
            <td><span style="color:#58a6ff;">${{d.protein || '—'}}</span> <span style="color:#8b949e; font-size:0.72rem;">${{d.coding || ''}}</span></td>
            <td>${{depthStr}}</td>
            <td><strong>${{vafStr}}</strong></td>
            <td style="color:#3fb950; font-weight:700;">${{impliedStr}}</td>
            <td>${{btnHtml}}</td>
          `;
          dtBody.appendChild(tr);
        }});
      }} else {{
        dtBody.innerHTML = '<tr><td colspan="6" style="text-align:center; color:#8b949e;">No reportable driver variants identified</td></tr>';
      }}

      // Candidates Table
      const ctBody = document.getElementById('candidatesTableBody');
      ctBody.innerHTML = '';
      if (c.candidates && c.candidates.length > 0) {{
        c.candidates.forEach(cand => {{
          const tr = document.createElement('tr');
          const pVal = (cand.purity * 100).toFixed(1);
          const plVal = cand.ploidy.toFixed(2);
          const btnHtml = `<button class="btn" style="padding:0.2rem 0.5rem; font-size:0.72rem;" onclick="applyPreset(${{pVal}}, ${{plVal}}, 'PURPLE 2D contour valley')">Adopt Fit</button>`;
          tr.innerHTML = `
            <td><strong>${{pVal}}%</strong></td>
            <td>${{plVal}}</td>
            <td>${{cand.score.toFixed(4)}}</td>
            <td>${{cand.somatic_penalty.toFixed(4)}}</td>
            <td>${{(cand.diploid_prop * 100).toFixed(1)}}%</td>
            <td>${{btnHtml}}</td>
          `;
          ctBody.appendChild(tr);
        }});
      }} else {{
        ctBody.innerHTML = '<tr><td colspan="6" style="text-align:center; color:#8b949e;">No alternative candidate valleys available</td></tr>';
      }}
    }}

    function setupIchorTab(c) {{
      const sel = document.getElementById('selIchorSol');
      const tbody = document.getElementById('ichorTableBody');
      sel.innerHTML = '';
      tbody.innerHTML = '';

      if (!c.ichor || !c.ichor.has_ichor || !c.ichor.solutions || c.ichor.solutions.length === 0) {{
        sel.innerHTML = '<option value="">No ichorCNA data available</option>';
        tbody.innerHTML = '<tr><td colspan="6" style="text-align:center; color:#8b949e;">No ichorCNA results found</td></tr>';
        document.getElementById('ichorPlotFrame').src = '';
        return;
      }}

      // Populate options
      c.ichor.solutions.forEach((sol, sIdx) => {{
        const opt = document.createElement('option');
        opt.value = sIdx;
        const optPfx = sol.is_optimal ? '★ [Optimal] ' : '';
        const tfStr = sol.tumor_fraction !== null ? (sol.tumor_fraction * 100).toFixed(1) + '%' : 'NA';
        const plStr = sol.ploidy !== null ? sol.ploidy.toFixed(2) : 'NA';
        const llStr = sol.loglik !== null ? sol.loglik.toFixed(0) : 'NA';
        opt.textContent = `${{optPfx}}${{sol.init}} — TF: ${{tfStr}}, Ploidy: ${{plStr}}, LogLik: ${{llStr}}`;
        sel.appendChild(opt);

        // Populate Table
        const tr = document.createElement('tr');
        if (sol.is_optimal) tr.style.background = 'rgba(35, 134, 54, 0.12)';
        const btnAdopt = `<button class="btn" style="padding:0.2rem 0.5rem; font-size:0.72rem;" onclick="adoptIchorByIndex(${{sIdx}})">Adopt Fit</button>`;
        const btnView = `<button class="btn" style="padding:0.2rem 0.5rem; font-size:0.72rem;" onclick="selectIchorSolutionByIndex(${{sIdx}})">View</button>`;
        const badgeOpt = sol.is_optimal ? ' <span class="badge badge-optimal" style="font-size:0.6rem;">TOP</span>' : '';
        
        tr.innerHTML = `
          <td><strong>${{sol.init}}</strong>${{badgeOpt}}</td>
          <td>${{tfStr}}</td>
          <td>${{plStr}}</td>
          <td>${{llStr}}</td>
          <td>${{sol.subclonal_frac !== null ? sol.subclonal_frac.toFixed(2) : '0'}}</td>
          <td>${{btnAdopt}} ${{btnView}}</td>
        `;
        tbody.appendChild(tr);
      }});

      // Add "All Solutions Overview" option
      if (c.ichor.all_sols_rel || c.ichor.all_sols_pdf) {{
        const optAll = document.createElement('option');
        optAll.value = 'ALL';
        optAll.textContent = '📑 All Solutions Multi-Page Overview PDF';
        sel.appendChild(optAll);
      }}

      renderSelectedIchorSolution();
    }}

    function selectIchorSolutionByIndex(sIdx) {{
      document.getElementById('selIchorSol').value = sIdx;
      renderSelectedIchorSolution();
    }}

    function renderSelectedIchorSolution() {{
      const c = CASES[activeIdx];
      if (!c.ichor || !c.ichor.solutions) return;

      const selVal = document.getElementById('selIchorSol').value;
      const iframe = document.getElementById('ichorPlotFrame');
      const img = document.getElementById('ichorPlotImg');
      const linkPopout = document.getElementById('linkIchorSolPdf');

      if (selVal === 'ALL') {{
        const pdfSrc = c.ichor.all_sols_rel || c.ichor.all_sols_pdf;
        iframe.src = pdfSrc;
        iframe.style.display = 'block';
        img.style.display = 'none';
        linkPopout.href = pdfSrc;
        return;
      }}

      const sIdx = parseInt(selVal);
      const sol = c.ichor.solutions[sIdx];
      if (!sol) return;

      const pdfSrc = sol.pdf_rel || sol.pdf_path || c.ichor.optimal_pdf_rel || c.ichor.optimal_pdf;
      linkPopout.href = pdfSrc;

      if (sol.png_base64 || sol.png_rel || sol.png_path) {{
        img.src = sol.png_base64 || sol.png_rel || sol.png_path;
        img.style.display = 'block';
        iframe.style.display = 'none';
      }} else if (pdfSrc) {{
        iframe.src = pdfSrc;
        iframe.style.display = 'block';
        img.style.display = 'none';
      }} else {{
        iframe.style.display = 'none';
        img.style.display = 'none';
      }}
    }}

    function adoptActiveIchorSolution() {{
      const selVal = document.getElementById('selIchorSol').value;
      if (selVal === 'ALL') return;
      adoptIchorByIndex(parseInt(selVal));
    }}

    function adoptIchorByIndex(sIdx) {{
      const c = CASES[activeIdx];
      if (!c.ichor || !c.ichor.solutions) return;
      const sol = c.ichor.solutions[sIdx];
      if (!sol) return;

      const tfPct = sol.tumor_fraction !== null ? (sol.tumor_fraction * 100).toFixed(1) : '15.0';
      const pl = sol.ploidy !== null ? sol.ploidy.toFixed(2) : '2.00';
      const note = `sWGS ichorCNA concordant: ${{sol.init}} (TF ${{tfPct}}%, loglik ${{sol.loglik}})`;
      
      applyPreset(tfPct, pl, note);
      document.getElementById('selMethod').value = 'sWGS ichorCNA concordant';
      updateActiveCase();
    }}

    function applyPreset(pur, pl, note) {{
      if (pur !== null) document.getElementById('inpPurity').value = parseFloat(pur).toFixed(1);
      if (pl !== null) document.getElementById('inpPloidy').value = parseFloat(pl).toFixed(2);
      if (note) {{
        document.getElementById('inpNotes').value = note;
        if (note.includes('VAF')) document.getElementById('selMethod').value = 'Clonal VAF anchor';
        else if (note.includes('ichorCNA')) document.getElementById('selMethod').value = 'sWGS ichorCNA concordant';
        else if (note.includes('contour')) document.getElementById('selMethod').value = 'PURPLE 2D contour valley';
      }}
      updateActiveCase();
    }}

    function updateActiveCase() {{
      const c = CASES[activeIdx];
      const pur = parseFloat(document.getElementById('inpPurity').value);
      const pl = parseFloat(document.getElementById('inpPloidy').value);
      const meth = document.getElementById('selMethod').value;
      const notes = document.getElementById('inpNotes').value;

      OVERRIDES[c.patient] = {{
        patient_id: c.patient,
        tumor_id: c.tumor,
        purity: pur,
        ploidy: pl,
        method: meth,
        notes: notes
      }};
      renderSampleList();
    }}

    function saveAndNext() {{
      updateActiveCase();
      if (activeIdx < CASES.length - 1) {{
        selectCase(activeIdx + 1);
      }}
    }}

    function generateTSV() {{
      let tsv = 'patient_id\\ttumor_id\\tcurated_purity\\tcurated_ploidy\\tcuration_method\\tnotes\\n';
      CASES.forEach(c => {{
        const ov = OVERRIDES[c.patient];
        if (ov) {{
          tsv += `${{ov.patient_id}}\\t${{ov.tumor_id}}\\t${{(ov.purity / 100.0).toFixed(4)}}\\t${{ov.ploidy.toFixed(2)}}\\t${{ov.method}}\\t${{ov.notes}}\\n`;
        }}
      }});
      return tsv;
    }}

    function exportOverridesTSV() {{
      const tsv = generateTSV();
      const blob = new Blob([tsv], {{ type: 'text/tab-separated-values' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'purity_overrides.tsv';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }}

    function copyOverridesTSV() {{
      const tsv = generateTSV();
      navigator.clipboard.writeText(tsv).then(() => {{
        alert('purity_overrides.tsv copied to clipboard!');
      }}).catch(err => {{
        alert('Failed to copy: ' + err);
      }});
    }}

    window.onload = init;
  </script>
</body>
</html>
"""
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Interactive curation console successfully generated: {output_html}")


def main():
    parser = argparse.ArgumentParser(description="Generate Purity Curation Console and overrides.tsv")
    parser.add_argument("--results-dir", required=True, help="Directory containing Oncoanalyser output cases")
    parser.add_argument("--output-html", default="purity_curation_dashboard.html", help="Path to output HTML curation console")
    parser.add_argument("--overrides-file", help="Existing purity_overrides.tsv to load")
    parser.add_argument("--embed-images", action="store_true", help="Embed contour and CNV images as base64 in HTML for offline portability")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        sys.exit(f"Error: {results_dir} is not a directory")

    out_html_path = Path(args.output_html).resolve()
    out_dir = out_html_path.parent

    # Load existing overrides if provided
    overrides = {}
    if args.overrides_file and os.path.isfile(args.overrides_file):
        with open(args.overrides_file, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                pid = r.get("patient_id")
                if pid:
                    overrides[pid] = {
                        "patient_id": pid,
                        "tumor_id": r.get("tumor_id"),
                        "purity": safe_float(r.get("curated_purity")) * 100.0 if safe_float(r.get("curated_purity")) <= 1.0 else safe_float(r.get("curated_purity")),
                        "ploidy": safe_float(r.get("curated_ploidy")),
                        "method": r.get("curation_method", "Manual"),
                        "notes": r.get("notes", "")
                    }

    # Discover cases
    cases = []
    for pat_dir in sorted(results_dir.iterdir()):
        if not pat_dir.is_dir() or pat_dir.name.startswith((".", "sample_sheets", "tmp", "logs", "wgts_summary")):
            continue
        patient_id = pat_dir.name
        
        # Look for purple directories directly at pat_dir / tumor_id / purple or pat_dir / purple
        cand_dirs = []
        if (pat_dir / "purple").is_dir():
            cand_dirs.append(pat_dir / "purple")
        for sub in pat_dir.iterdir():
            if sub.is_dir() and not sub.name.startswith((".", "logs", "tmp", "plasma", "purity", "ichorCNA")):
                if (sub / "purple").is_dir():
                    cand_dirs.append(sub / "purple")

        for p_dir in cand_dirs:
            tumor_id = p_dir.parent.name
            
            purity_tsv = p_dir / f"{tumor_id}.purple.purity.tsv"
            if not purity_tsv.is_file():
                # Try matching any *.purple.purity.tsv
                cand = list(p_dir.glob("*.purple.purity.tsv"))
                if cand:
                    purity_tsv = cand[0]
                    tumor_id = purity_tsv.name.replace(".purple.purity.tsv", "")
            
            p_stats = parse_purple_purity_tsv(purity_tsv)
            range_tsv = p_dir / f"{tumor_id}.purple.purity.range.tsv"
            candidates = parse_purity_range_tsv(range_tsv)
            
            # Orange Report
            orange_data = discover_orange_report(pat_dir, tumor_id, p_dir=p_dir, out_dir=out_dir)
            
            # Drivers
            cat_somatic = p_dir / f"{tumor_id}.purple.driver.catalog.somatic.tsv"
            drivers = parse_drivers(orange_data.get("json_path"), cat_somatic)
            
            # ichorCNA
            ichor_data = discover_ichor_data(pat_dir, tumor_id, out_dir=out_dir, embed_images=args.embed_images)
            
            # Contour Plot
            plot_png = p_dir / "plot" / f"{tumor_id}.purity.range.png"
            plot_b64 = ""
            plot_rel = os.path.relpath(plot_png, out_dir) if plot_png.is_file() else ""
            if args.embed_images and plot_png.is_file():
                try:
                    with open(plot_png, "rb") as img_f:
                        plot_b64 = "data:image/png;base64," + base64.b64encode(img_f.read()).decode("ascii")
                except Exception:
                    pass

            cases.append({
                "patient": patient_id,
                "tumor": tumor_id,
                "current_purity": p_stats.get("purity"),
                "current_ploidy": p_stats.get("ploidy"),
                "status": p_stats.get("status", "NORMAL"),
                "wgd": p_stats.get("wgd", False),
                "tmb": p_stats.get("tmb"),
                "drivers": drivers,
                "candidates": candidates,
                "orange": orange_data,
                "ichor": ichor_data,
                "plot_path": str(plot_png),
                "plot_rel": plot_rel,
                "plot_base64": plot_b64
            })

    print(f"Discovered {len(cases)} cases across {results_dir}")
    build_curation_dashboard(cases, str(out_html_path), existing_overrides=overrides)


if __name__ == "__main__":
    main()
