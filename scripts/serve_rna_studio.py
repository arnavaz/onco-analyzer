#!/usr/bin/env python3
"""
Pan-Cancer RNA-Seq Clinical Analytics Studio Server & Interface
==============================================================
Serves an exploratory, anti-slop analytical clinical web application
backed by DuckDB / SQLite RNA-Seq database.

Supports:
1. Cross-cancer gene distribution queries (boxplots, strip plots).
2. Patient over-expression threshold matrix benchmarked against HMF baseline.
3. Dynamic Group vs. Group differential expression (Volcano plot & FDR table).
4. Whole-cohort dimensionality reduction landscape (UMAP, t-SNE, PCA).
5. Tumour Mutation Burden (TMB) & somatic mutation load analysis.

Usage:
  python3 scripts/serve_rna_studio.py --db hmf_rnaseq_cohort.duckdb --port 8080
"""

import os
import sys
import json
import math
import sqlite3
import argparse
import threading
import urllib.parse
from collections import defaultdict
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from pathlib import Path

# DuckDB support if available
HAS_DUCKDB = False
try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False


class StudioDataStore:
    """Analytical data access layer over SQLite3 & DuckDB with thread-local connections."""

    def __init__(self, db_path):
        self.db_path = Path(db_path).resolve()
        if not self.db_path.is_file():
            raise FileNotFoundError(f"Database not found: {self.db_path}")

        # Inspect binary magic bytes to determine actual database format
        with open(self.db_path, "rb") as f:
            hdr = f.read(16)
            self.is_sqlite = hdr.startswith(b"SQLite format 3")

        self.use_duckdb = HAS_DUCKDB
        self._local = threading.local()
        self._diff_expr_cache = {}

        # One-time index verification for SQLite
        if self.is_sqlite:
            init_conn = sqlite3.connect(str(self.db_path))
            init_conn.execute("CREATE INDEX IF NOT EXISTS idx_gene_name_nocase ON gene_expression (gene_name COLLATE NOCASE)")
            init_conn.execute("CREATE INDEX IF NOT EXISTS idx_trans_gene_nocase ON transcript_expression (gene_name COLLATE NOCASE)")
            init_conn.execute("CREATE INDEX IF NOT EXISTS idx_alt_gene_nocase ON alt_splice_junctions (gene_name COLLATE NOCASE)")
            init_conn.commit()
            init_conn.close()

        # Fast startup: no heavy background table scans on 77M rows
        pass

    def _get_conn(self):
        if not hasattr(self._local, "conn") or self._local.conn is None:
            if not self.is_sqlite and self.use_duckdb:
                self._local.conn = duckdb.connect(str(self.db_path), read_only=True)
            else:
                c = sqlite3.connect(str(self.db_path), check_same_thread=False)
                c.row_factory = sqlite3.Row
                c.execute("PRAGMA cache_size = -64000")
                c.execute("PRAGMA mmap_size = 2147483648")
                self._local.conn = c
        return self._local.conn

    def get_cursor(self):
        return self._get_conn().cursor()

    def get_cohort_overview(self):
        """Get summary statistics, cancer types, and available clinical attributes."""
        cur = self.get_cursor()

        # Metadata
        cur.execute("SELECT key, value FROM database_metadata")
        meta = {row[0]: row[1] for row in cur.fetchall()}

        # Total samples and genes (fast lookup without scanning 77M rows)
        cur.execute("SELECT COUNT(*) FROM clinical_samples")
        n_samples = cur.fetchone()[0]
        n_genes = int(meta.get("total_genes") or 39500)

        # Cancer types distribution
        cur.execute("""
            SELECT cancer_type, COUNT(*) as count
            FROM clinical_samples
            GROUP BY cancer_type
            ORDER BY count DESC
        """)
        cancer_types = [{"cancer_type": row[0], "count": row[1]} for row in cur.fetchall()]

        # Metadata
        cur.execute("SELECT key, value FROM database_metadata")
        meta = {row[0]: row[1] for row in cur.fetchall()}

        # Sample clinical attributes available for grouping
        cur.execute("SELECT * FROM clinical_samples LIMIT 1")
        sample_row = cur.fetchone()
        clin_cols = []
        if sample_row:
            if hasattr(sample_row, "keys"):
                clin_cols = [k for k in sample_row.keys() if k != "clinical_json"]
            else:
                desc = cur.description
                clin_cols = [d[0] for d in desc if d[0] != "clinical_json"]

        return {
            "total_samples": n_samples,
            "total_genes": n_genes,
            "cancer_types": cancer_types,
            "clinical_columns": clin_cols,
            "metadata": meta,
        }

    def get_gene_distribution(self, gene_name, metric="log2_adj_tpm", group_by="cancer_type"):
        """Get gene expression distribution across cancer types or another clinical group."""
        cur = self.get_cursor()

        valid_metrics = {
            "log2_adj_tpm": "ge.log2_adj_tpm",
            "adj_tpm": "ge.adj_tpm",
            "norm_counts": "ge.norm_counts",
            "log2_norm_counts": "ge.log2_norm_counts",
            "batch_corrected": "ge.batch_corrected_expr",
            "raw_tpm": "ge.raw_tpm",
        }
        metric_col = valid_metrics.get(metric, "ge.log2_adj_tpm")

        query = f"""
            SELECT 
                cs.{group_by} as group_val,
                ge.sample_id,
                {metric_col} as expr_val,
                ge.adj_tpm,
                ge.percentile_cancer,
                ge.median_tpm_cancer,
                cs.cancer_type,
                cs.purity,
                cs.cuppa_predicted_cancer_type,
                cs.ecdna_status
            FROM gene_expression ge
            JOIN clinical_samples cs ON ge.sample_id = cs.sample_id
            WHERE ge.gene_name = ? COLLATE NOCASE
            ORDER BY cs.{group_by}, expr_val ASC
        """
        cur.execute(query, [gene_name])
        rows = cur.fetchall()

        # Group data points for boxplot calculations
        groups = {}
        all_points = []
        for r in rows:
            g_val = str(r[0] or "Unknown")
            val = float(r[2]) if r[2] is not None else 0.0
            pt = {
                "sample_id": r[1],
                "expr": round(val, 2),
                "adj_tpm": round(float(r[3] or 0.0), 2),
                "pct_cancer": round(float(r[4] or -1.0), 1),
                "med_cancer": round(float(r[5] or -1.0), 2),
                "cancer_type": r[6],
                "purity": round(float(r[7] or -1.0), 2),
                "cuppa": r[8],
                "ecdna": r[9],
            }
            if g_val not in groups:
                groups[g_val] = []
            groups[g_val].append(pt)
            all_points.append(pt)

        # Compute summary stats per group
        group_summaries = []
        for g_name, pts in groups.items():
            vals = [p["expr"] for p in pts]
            vals.sort()
            n = len(vals)
            if n == 0:
                continue

            q25 = vals[int(0.25 * n)]
            q50 = vals[int(0.50 * n)]
            q75 = vals[int(0.75 * n)]
            iqr = q75 - q25
            low_whisker = max(vals[0], q25 - 1.5 * iqr)
            high_whisker = min(vals[-1], q75 + 1.5 * iqr)

            group_summaries.append({
                "group": g_name,
                "count": n,
                "median": round(q50, 2),
                "q25": round(q25, 2),
                "q75": round(q75, 2),
                "min": round(low_whisker, 2),
                "max": round(high_whisker, 2),
                "points": pts,
            })

        group_summaries.sort(key=lambda x: x["median"], reverse=True)

        return {
            "gene": gene_name.upper(),
            "metric": metric,
            "group_by": group_by,
            "total_points": len(all_points),
            "groups": group_summaries,
        }

    def get_overexpression_matrix(self, gene_name, threshold_pct=90.0):
        """Compute proportion of patients with over-expression relative to HMF benchmark."""
        cur = self.get_cursor()

        # Normalize cutoff: e.g. 90 -> 0.90
        cutoff_val = (threshold_pct / 100.0) if threshold_pct > 1.0 else threshold_pct

        query = """
            SELECT 
                cs.cancer_type,
                ge.sample_id,
                ge.adj_tpm,
                ge.percentile_cancer,
                ge.median_tpm_cancer,
                ge.percentile_cohort
            FROM gene_expression ge
            JOIN clinical_samples cs ON ge.sample_id = cs.sample_id
            WHERE ge.gene_name = ? COLLATE NOCASE
        """
        cur.execute(query, [gene_name])
        rows = cur.fetchall()

        by_cancer = {}
        for r in rows:
            c_type = str(r[0] or "Unknown")
            if c_type not in by_cancer:
                by_cancer[c_type] = {"total": 0, "elevated": 0, "samples": []}

            by_cancer[c_type]["total"] += 1
            pct_cancer = float(r[3] or -1.0)
            pct_cohort = float(r[5] or -1.0)

            # Check against tumor-matched cancer baseline if available, otherwise pan-cancer cohort baseline
            ref_pct = pct_cancer if pct_cancer >= 0 else pct_cohort
            is_elevated = (ref_pct >= cutoff_val) if ref_pct >= 0 else False

            if is_elevated:
                by_cancer[c_type]["elevated"] += 1

            by_cancer[c_type]["samples"].append({
                "sample_id": r[1],
                "adj_tpm": round(float(r[2] or 0.0), 2),
                "percentile_cancer": pct_cancer,
                "percentile_cohort": pct_cohort,
                "is_elevated": is_elevated,
            })

        matrix = []
        for c_type, stats in by_cancer.items():
            tot = stats["total"]
            el = stats["elevated"]
            matrix.append({
                "cancer_type": c_type,
                "total_patients": tot,
                "count_elevated": el,
                "percent_elevated": round((el / tot) * 100.0, 1) if tot > 0 else 0.0,
            })

        matrix.sort(key=lambda x: (x["percent_elevated"], x["count_elevated"]), reverse=True)

        return {
            "gene": gene_name.upper(),
            "threshold_percentile": threshold_pct,
            "cancer_types": matrix,
        }

    def compute_differential_expression(self, contrast_col, group_a_val, group_b_val, min_samples=3):
        """
        Run Welch's t-test and compute log2 fold changes between Group A and Group B
        using in-engine SQL aggregation (instant, avoids pulling 60M raw tuples into Python).
        """
        contrast_col = contrast_col.lower()
        cache_key = (contrast_col, group_a_val, group_b_val)
        if cache_key in self._diff_expr_cache:
            return self._diff_expr_cache[cache_key]

        cur = self.get_cursor()

        # 1. Fetch sample counts in Group A and Group B
        cur.execute(f"SELECT COUNT(*) FROM clinical_samples WHERE {contrast_col} = ?", [group_a_val])
        n_a = cur.fetchone()[0]

        cur.execute(f"SELECT COUNT(*) FROM clinical_samples WHERE {contrast_col} = ?", [group_b_val])
        n_b = cur.fetchone()[0]

        if n_a < min_samples or n_b < min_samples:
            return {
                "error": f"Insufficient sample size: Group '{group_a_val}' has {n_a} samples, Group '{group_b_val}' has {n_b} samples (minimum {min_samples} required)."
            }

        # 2. Pushdown aggregation: Use DuckDB if available for 10x parallel speedup
        if self.use_duckdb:
            try:
                d_con = duckdb.connect()
                d_con.load_extension("sqlite")
                d_con.execute(f"ATTACH '{self.db_path}' AS hmf (TYPE SQLITE)")
                query = f"""
                    SELECT 
                        ge.gene_name,
                        AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm END) as mean_a,
                        AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm * ge.log2_adj_tpm END) as sq_a,
                        AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm END) as mean_b,
                        AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm * ge.log2_adj_tpm END) as sq_b
                    FROM hmf.gene_expression ge
                    JOIN hmf.clinical_samples cs ON ge.sample_id = cs.sample_id
                    WHERE cs.{contrast_col} IN (?, ?)
                    GROUP BY ge.gene_name
                    HAVING mean_a IS NOT NULL AND mean_b IS NOT NULL
                """
                rows = d_con.execute(query, [group_a_val, group_a_val, group_b_val, group_b_val, group_a_val, group_b_val]).fetchall()
                d_con.close()
            except Exception:
                # Fallback to SQLite cursor if DuckDB attach encounters any issue
                rows = []
        else:
            rows = []

        if not rows:
            query = f"""
                SELECT 
                    ge.gene_name,
                    AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm END) as mean_a,
                    AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm * ge.log2_adj_tpm END) as sq_a,
                    AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm END) as mean_b,
                    AVG(CASE WHEN cs.{contrast_col} = ? THEN ge.log2_adj_tpm * ge.log2_adj_tpm END) as sq_b
                FROM gene_expression ge
                JOIN clinical_samples cs ON ge.sample_id = cs.sample_id
                WHERE cs.{contrast_col} IN (?, ?)
                GROUP BY ge.gene_name
                HAVING mean_a IS NOT NULL AND mean_b IS NOT NULL
            """
            cur.execute(query, [group_a_val, group_a_val, group_b_val, group_b_val, group_a_val, group_b_val])
            rows = cur.fetchall()

        results = []
        for r in rows:
            gene = r[0]
            ma = float(r[1] or 0.0)
            sqa = float(r[2] or 0.0)
            mb = float(r[3] or 0.0)
            sqb = float(r[4] or 0.0)

            # Filter out low-expressed noise (mean < 0.2 in both groups)
            if ma < 0.2 and mb < 0.2:
                continue

            fc = round(ma - mb, 3)

            # Sample variance with moderate variance floor (avoids artificial p=1e-30 from zero variance)
            vara = max(0.05, (sqa - ma * ma) * (n_a / max(1, n_a - 1)))
            varb = max(0.05, (sqb - mb * mb) * (n_b / max(1, n_b - 1)))

            denom = math.sqrt(vara / n_a + varb / n_b + 1e-6)
            t_stat = (ma - mb) / denom
            z = abs(t_stat)
            p_val = math.erfc(z / math.sqrt(2.0))
            p_val = max(p_val, 1e-20)

            results.append({
                "gene": gene,
                "log2_fc": fc,
                "mean_a": round(ma, 2),
                "mean_b": round(mb, 2),
                "p_value": p_val,
                "neg_log10_p": round(-math.log10(p_val), 2),
            })

        # Benjamini-Hochberg FDR
        results.sort(key=lambda x: x["p_value"])
        m = len(results)
        for i, res in enumerate(results, 1):
            fdr = min(1.0, res["p_value"] * m / i)
            res["fdr"] = fdr
            res["neg_log10_fdr"] = round(-math.log10(max(fdr, 1e-30)), 2)

        # Balanced volcano point selection: UP + DOWN + representative background
        sig_up = [r for r in results if r["log2_fc"] >= 0.5 and r["fdr"] <= 0.05][:300]
        sig_down = [r for r in results if r["log2_fc"] <= -0.5 and r["fdr"] <= 0.05][:300]
        other = [r for r in results if r not in sig_up and r not in sig_down]
        step = max(1, len(other) // 400)
        bg_sample = other[::step][:400]

        data = {
            "contrast_col": contrast_col,
            "group_a": group_a_val,
            "group_b": group_b_val,
            "n_a": n_a,
            "n_b": n_b,
            "total_tested": m,
            "top_genes": results[:50],
            "volcano_points": [
                {
                    "gene": r["gene"],
                    "x": r["log2_fc"],
                    "y": r["neg_log10_p"],
                    "log2_fc": r["log2_fc"],
                    "neg_log10_p": r["neg_log10_p"],
                    "fdr": r["fdr"]
                }
                for r in (sig_up + sig_down + bg_sample)
            ],
        }
        self._diff_expr_cache[cache_key] = data
        return data

    def get_multi_gene_overexpression(
        self, 
        gene_list, 
        threshold_pct=90.0, 
        selected_cancers=None,
        min_purity=0.0,
        weight_mode="bulk",
        purity_floor=0.15
    ):
        """
        Compute multi-gene over-expression counts, cancer breakdown, and co-expression
        with optional tumor type filtering, purity threshold filtering, and tumor-cell / TME weighting.
        """
        if not gene_list:
            return {"genes": [], "total_cohort": 0, "matrix": [], "gene_summary": []}

        # Keep unique in order
        seen = set()
        clean_genes = []
        for g in gene_list:
            gu = g.strip().upper()
            if gu and gu not in seen:
                clean_genes.append(gu)
                seen.add(gu)

        cutoff_val = (threshold_pct / 100.0) if threshold_pct > 1.0 else threshold_pct
        min_purity = float(min_purity) if min_purity is not None else 0.0
        purity_floor = float(purity_floor) if purity_floor is not None else 0.15
        weight_mode = str(weight_mode or "bulk").lower()

        cur = self.get_cursor()

        # Clean selected cancers
        clean_cancers = []
        if selected_cancers:
            if isinstance(selected_cancers, str):
                selected_cancers = selected_cancers.split(",")
            for c in selected_cancers:
                c_str = c.strip()
                if c_str and c_str.upper() != "ALL":
                    clean_cancers.append(c_str)

        # Build clinical filter clause
        clin_where = []
        clin_params = []

        if min_purity > 0:
            clin_where.append("purity >= ?")
            clin_params.append(min_purity)

        if clean_cancers:
            c_placeholders = ','.join(['?'] * len(clean_cancers))
            clin_where.append(f"cancer_type IN ({c_placeholders})")
            clin_params.extend(clean_cancers)

        clin_sql = ("WHERE " + " AND ".join(clin_where)) if clin_where else ""

        # Fetch clinical cancer counts
        cur.execute(
            f"SELECT cancer_type, COUNT(*) FROM clinical_samples {clin_sql} GROUP BY cancer_type ORDER BY COUNT(*) DESC",
            clin_params
        )
        cancer_counts = {r[0]: r[1] for r in cur.fetchall()}
        total_cohort = sum(cancer_counts.values())

        if total_cohort == 0:
            return {
                "genes": clean_genes,
                "selected_cancers": clean_cancers,
                "min_purity": min_purity,
                "weight_mode": weight_mode,
                "purity_floor": purity_floor,
                "below_floor_count": 0,
                "uncalled_purity_count": 0,
                "threshold_percentile": threshold_pct,
                "total_cohort": 0,
                "any_target_count": 0,
                "any_target_pct": 0.0,
                "co_elevated_count": 0,
                "co_elevated_pct": 0.0,
                "gene_summary": [{"gene": g, "elevated_count": 0, "total_cohort": 0, "percent_elevated": 0.0} for g in clean_genes],
                "cancer_matrix": []
            }

        # Fetch expression for all requested genes
        gene_placeholders = ','.join(['?'] * len(clean_genes))
        expr_where = [f"ge.gene_name IN ({gene_placeholders}) COLLATE NOCASE"]
        expr_params = list(clean_genes)

        if min_purity > 0:
            expr_where.append("cs.purity >= ?")
            expr_params.append(min_purity)

        if clean_cancers:
            c_placeholders = ','.join(['?'] * len(clean_cancers))
            expr_where.append(f"cs.cancer_type IN ({c_placeholders})")
            expr_params.extend(clean_cancers)

        expr_sql = "WHERE " + " AND ".join(expr_where)

        query = f"""
            SELECT 
                ge.gene_name,
                ge.sample_id,
                cs.cancer_type,
                ge.adj_tpm,
                cs.purity,
                ge.percentile_cancer,
                ge.percentile_cohort
            FROM gene_expression ge
            JOIN clinical_samples cs ON ge.sample_id = cs.sample_id
            {expr_sql}
        """
        cur.execute(query, expr_params)
        rows = cur.fetchall()

        # Track purity stats across all evaluated samples
        sample_purity_map = {}
        sample_cancer_map = {}
        for r in rows:
            sid = r[1]
            ct = str(r[2] or "Unknown")
            p = float(r[4]) if (r[4] is not None) else -1.0
            sample_purity_map[sid] = p
            sample_cancer_map[sid] = ct

        total_evaluated_samples = len(sample_purity_map)
        below_floor_samples = set(sid for sid, p in sample_purity_map.items() if 0 < p < purity_floor)
        uncalled_purity_samples = set(sid for sid, p in sample_purity_map.items() if p <= 0)

        # Organize rows by gene
        gene_rows = {g: [] for g in clean_genes}
        for r in rows:
            g = str(r[0]).upper()
            if g in gene_rows:
                gene_rows[g].append(r)

        # Compute elevation per gene
        import numpy as np
        gene_stats = {g: {"gene": g, "elevated_count": 0, "below_floor_elevated": 0} for g in clean_genes}
        sample_hits = {}

        cancer_matrix = {
            ct: {
                "cancer_type": ct,
                "total_patients": total_n,
                "gene_counts": {g: 0 for g in clean_genes},
                "below_floor_count": sum(1 for sid in below_floor_samples if sample_cancer_map.get(sid) == ct),
                "uncalled_purity_count": sum(1 for sid in uncalled_purity_samples if sample_cancer_map.get(sid) == ct),
                "_sample_hits": set()
            }
            for ct, total_n in cancer_counts.items()
        }

        for g in clean_genes:
            g_data = gene_rows[g]
            if not g_data:
                continue

            if weight_mode == "tumor_adjusted":
                # Adjusted expression = AdjTPM / max(purity, floor)
                # Uncalled samples assigned cohort median 0.40
                vals = []
                for r in g_data:
                    tpm = float(r[3] or 0.0)
                    p = float(r[4]) if (r[4] is not None) else -1.0
                    eff_p = max(p, purity_floor) if p > 0 else 0.40
                    vals.append(tpm / eff_p)
                th_val = float(np.percentile(vals, threshold_pct)) if vals else 0.0
                is_el_list = [v >= th_val for v in vals]

            elif weight_mode == "tme_weighted":
                # TME expression = AdjTPM * (1 - clamped_purity)
                vals = []
                for r in g_data:
                    tpm = float(r[3] or 0.0)
                    p = float(r[4]) if (r[4] is not None) else -1.0
                    eff_p = min(max(p, 0.0), 1.0) if p > 0 else 0.40
                    vals.append(tpm * (1.0 - eff_p))
                th_val = float(np.percentile(vals, threshold_pct)) if vals else 0.0
                is_el_list = [v >= th_val for v in vals]

            else:
                # Bulk precomputed reference percentile
                is_el_list = []
                for r in g_data:
                    pct_c = float(r[5] or -1.0)
                    pct_coh = float(r[6] or -1.0)
                    ref_pct = pct_c if pct_c >= 0 else pct_coh
                    is_el_list.append(ref_pct >= cutoff_val if ref_pct >= 0 else False)

            # Record elevated samples
            for r, is_el in zip(g_data, is_el_list):
                if is_el:
                    sid = r[1]
                    ct = str(r[2] or "Unknown")
                    gene_stats[g]["elevated_count"] += 1
                    if sid in below_floor_samples:
                        gene_stats[g]["below_floor_elevated"] += 1

                    if ct in cancer_matrix:
                        cancer_matrix[ct]["gene_counts"][g] += 1
                        cancer_matrix[ct]["_sample_hits"].add(sid)

                    sample_hits.setdefault(sid, set()).add(g)

        # Finalize gene summary list
        gene_summary = []
        for g in clean_genes:
            cnt = gene_stats[g]["elevated_count"]
            tot = total_cohort
            pct = round((cnt / tot * 100.0), 1) if tot > 0 else 0.0
            gene_summary.append({
                "gene": g,
                "elevated_count": cnt,
                "total_cohort": tot,
                "percent_elevated": pct,
                "below_floor_elevated": gene_stats[g]["below_floor_elevated"]
            })

        # Finalize cancer breakdown table
        matrix_rows = []
        for ct, data in cancer_matrix.items():
            tot = data["total_patients"]
            any_c = len(data["_sample_hits"])
            any_pct = round((any_c / tot * 100.0), 1) if tot > 0 else 0.0
            row = {
                "cancer_type": ct,
                "total_patients": tot,
                "any_elevated_count": any_c,
                "any_elevated_pct": any_pct,
                "below_floor_count": data["below_floor_count"],
                "uncalled_purity_count": data["uncalled_purity_count"],
                "gene_counts": {g: data["gene_counts"].get(g, 0) for g in clean_genes},
                "gene_pcts": {g: round((data["gene_counts"].get(g, 0) / tot * 100.0), 1) if tot > 0 else 0.0 for g in clean_genes}
            }
            matrix_rows.append(row)

        matrix_rows.sort(key=lambda x: (x["any_elevated_pct"], x["any_elevated_count"]), reverse=True)

        co_elevated_count = sum(1 for tgts in sample_hits.values() if len(tgts) >= 2)
        any_elevated_count = len(sample_hits)

        return {
            "genes": clean_genes,
            "selected_cancers": clean_cancers,
            "min_purity": min_purity,
            "weight_mode": weight_mode,
            "purity_floor": purity_floor,
            "below_floor_count": len(below_floor_samples),
            "below_floor_pct": round(len(below_floor_samples) / max(1, total_evaluated_samples) * 100.0, 1),
            "uncalled_purity_count": len(uncalled_purity_samples),
            "uncalled_purity_pct": round(len(uncalled_purity_samples) / max(1, total_evaluated_samples) * 100.0, 1),
            "threshold_percentile": threshold_pct,
            "total_cohort": total_cohort,
            "any_target_count": any_elevated_count,
            "any_target_pct": round((any_elevated_count / total_cohort * 100.0), 1) if total_cohort > 0 else 0.0,
            "co_elevated_count": co_elevated_count,
            "co_elevated_pct": round((co_elevated_count / total_cohort * 100.0), 1) if total_cohort > 0 else 0.0,
            "gene_summary": gene_summary,
            "cancer_matrix": matrix_rows
        }
# ---------------------------------------------------------------------------
# HTTP Request Handler & Frontend Application
# ---------------------------------------------------------------------------

    def get_cohort_umap(self):
        """Fetch precomputed coordinates (t-SNE, PCA, UMAP) and sample metadata for immediate plotting."""
        cur = self.get_cursor()
        query = """
            SELECT 
                u.sample_id,
                u.tsne_1,
                u.tsne_2,
                u.pca_1,
                u.pca_2,
                u.pca_3,
                u.umap_1,
                u.umap_2,
                cs.cancer_type,
                cs.cancer_type_detailed,
                cs.cuppa_predicted_cancer_type,
                cs.cuppa_likelihood,
                cs.purity,
                u.specimen_category,
                u.expr_erbb2,
                u.expr_tacstd2,
                u.expr_cldn18,
                u.expr_folr1,
                u.expr_vtcn1,
                cs.clinical_json
            FROM cohort_umap_coordinates u
            JOIN clinical_samples cs ON u.sample_id = cs.sample_id
            ORDER BY u.sample_id
        """
        rows = cur.execute(query).fetchall()
        
        samples = []
        cancer_counts = {}
        cuppa_counts = {}
        specimen_counts = {}
        
        for r in rows:
            cancer = r[8] or "Unknown"
            cuppa = r[10] or "Unclassified"
            specimen = r[13] or "Unknown"
            cuppa_top = ""
            if r[19]:
                try:
                    c_dict = json.loads(r[19])
                    cuppa_top = c_dict.get("CUPPA_TOP_PREDICTIONS", "")
                except Exception:
                    pass
            
            cancer_counts[cancer] = cancer_counts.get(cancer, 0) + 1
            if cuppa != "NA":
                cuppa_counts[cuppa] = cuppa_counts.get(cuppa, 0) + 1
            specimen_counts[specimen] = specimen_counts.get(specimen, 0) + 1
            
            samples.append({
                "id": r[0],
                "tsne_x": round(float(r[1]), 2),
                "tsne_y": round(float(r[2]), 2),
                "pca12_x": round(float(r[3]), 2),
                "pca12_y": round(float(r[4]), 2),
                "pca13_x": round(float(r[3]), 2),
                "pca13_y": round(float(r[5]), 2),
                "umap_x": round(float(r[6]), 2),
                "umap_y": round(float(r[7]), 2),
                "cancer": cancer,
                "cancer_det": r[9] or "",
                "cuppa": cuppa,
                "cuppa_p": round(float(r[11] or 0.0), 3),
                "cuppa_top": cuppa_top,
                "purity": round(float(r[12]), 2) if r[12] is not None else None,
                "specimen": specimen,
                "erbb2": round(float(r[14] or 0.0), 2),
                "tacstd2": round(float(r[15] or 0.0), 2),
                "cldn18": round(float(r[16] or 0.0), 2),
                "folr1": round(float(r[17] or 0.0), 2),
                "vtcn1": round(float(r[18] or 0.0), 2)
            })
            
        return {
            "total_samples": len(samples),
            "cancer_types": sorted(list(cancer_counts.keys())),
            "cancer_counts": cancer_counts,
            "cuppa_types": sorted(list(cuppa_counts.keys())),
            "cuppa_counts": cuppa_counts,
            "specimen_types": ["Primary", "Metastasis", "Local Recurrence", "Other"],
            "specimen_counts": specimen_counts,
            "genes": ["ERBB2", "TACSTD2", "CLDN18", "FOLR1", "VTCN1"],
            "samples": samples
        }

    def get_tmb_data(self):
        """Fetch TMB, mutation load, purity, MSI, POLE signature, and clinical cohort for all biopsies."""
        cur = self.get_cursor()
        query = """
            SELECT 
                sample_id,
                cancer_type,
                cancer_type_detailed,
                cuppa_predicted_cancer_type,
                purity,
                ploidy,
                tmb_per_mb,
                msi_status,
                clinical_json
            FROM clinical_samples
            WHERE tmb_per_mb IS NOT NULL AND tmb_per_mb >= 0
            ORDER BY cancer_type, sample_id
        """
        rows = cur.execute(query).fetchall()
        samples = []
        for r in rows:
            pole_val = 0.0
            specimen_type = "Unknown"
            if r[8]:
                try:
                    d = json.loads(r[8])
                    v = d.get('SIG10_RECURRENT_POLE_SOMATIC_MUTATIONS')
                    if v:
                        pole_val = float(v)
                    st = d.get('SPECIMEN_TYPE', '')
                    if 'primary' in st.lower(): specimen_type = 'Primary'
                    elif 'metastat' in st.lower(): specimen_type = 'Metastasis'
                    elif 'recurrent' in st.lower(): specimen_type = 'Local Recurrence'
                    elif st: specimen_type = 'Other'
                except Exception:
                    pass
            
            tmb_val = round(float(r[6]), 3)
            samples.append({
                "id": r[0],
                "cancer": r[1] or "Unknown",
                "cancer_det": r[2] or "",
                "cuppa": r[3] or "Unclassified",
                "purity": round(float(r[4]), 2) if r[4] is not None and r[4] >= 0 else None,
                "ploidy": round(float(r[5]), 2) if r[5] is not None and r[5] >= 0 else None,
                "tmb": tmb_val,
                "mut_count": int(round(tmb_val * 2850)),
                "coding_mut": int(round(tmb_val * 32)),
                "msi": r[7] or "NA",
                "specimen": specimen_type,
                "pole_sig": round(pole_val, 2),
                "is_pole_mt": pole_val > 0.05
            })
        return {"total_samples": len(samples), "samples": samples}

class StudioRequestHandler(BaseHTTPRequestHandler):

    store = None

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        params = urllib.parse.parse_qs(parsed.query)

        if path == "/" or path == "/index.html":
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(get_studio_html().encode("utf-8"))
            return

        elif path == "/api/cohort_overview":
            data = self.store.get_cohort_overview()
            self._send_json(data)
            return

        elif path == "/api/cohort_umap":
            data = self.store.get_cohort_umap()
            self._send_json(data)
            return

        elif path == "/api/tmb_data":
            data = self.store.get_tmb_data()
            self._send_json(data)
            return

        elif path.startswith("/api/genes/") and "/distribution" in path:
            gene = path.split("/api/genes/")[1].split("/distribution")[0].strip()
            metric = params.get("metric", ["log2_adj_tpm"])[0]
            group_by = params.get("group_by", ["cancer_type"])[0]
            data = self.store.get_gene_distribution(gene, metric=metric, group_by=group_by)
            self._send_json(data)
            return

        elif path == "/api/genes/multi_overexpression":
            genes_raw = params.get("genes", ["ERBB2,TACSTD2,NECTIN4,FOLR1"])[0]
            gene_list = [g.strip().upper() for g in genes_raw.replace(";", ",").replace(" ", ",").split(",") if g.strip()]
            th = float(params.get("threshold", [90.0])[0])
            cancers_raw = params.get("cancers", [""])[0]
            selected_cancers = [c.strip() for c in cancers_raw.split(",") if c.strip() and c.strip().upper() != "ALL"]
            min_p = float(params.get("min_purity", [0.0])[0])
            w_mode = params.get("weight_mode", ["bulk"])[0]
            p_floor = float(params.get("purity_floor", [0.15])[0])
            data = self.store.get_multi_gene_overexpression(
                gene_list, 
                threshold_pct=th, 
                selected_cancers=selected_cancers,
                min_purity=min_p,
                weight_mode=w_mode,
                purity_floor=p_floor
            )
            self._send_json(data)
            return

        elif path.startswith("/api/genes/") and "/overexpression" in path:
            gene = path.split("/api/genes/")[1].split("/overexpression")[0].strip()
            th = float(params.get("threshold", [90.0])[0])
            data = self.store.get_overexpression_matrix(gene, threshold_pct=th)
            self._send_json(data)
            return

        elif path == "/api/diff_expr":
            col = params.get("col", ["ecdna_status"])[0]
            grp_a = params.get("group_a", ["Yes"])[0]
            grp_b = params.get("group_b", ["No"])[0]
            data = self.store.compute_differential_expression(col, grp_a, grp_b)
            self._send_json(data)
            return

        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(b"Not Found")

    def _send_json(self, payload):
        self.send_response(200)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(json.dumps(payload).encode("utf-8"))

    def log_message(self, format, *args):
        # Quiet logger for clean CLI output
        pass


def get_studio_html():
    """Returns the standalone publication-grade Clinical Analytics Studio HTML interface with export capabilities."""
    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pan-Cancer RNA-Seq Clinical Studio</title>
<style>
  :root {
    --bg-page: #f8fafc;
    --bg-surface: #ffffff;
    --border: #e2e8f0;
    --border-light: #f1f5f9;
    --text-primary: #0f172a;
    --text-secondary: #475569;
    --text-muted: #94a3b8;
    --accent: #2563eb;
    --accent-hover: #1d4ed8;
    --accent-light: #eff6ff;
    --up-reg: #059669;
    --down-reg: #e11d48;
    --sidebar-w: 240px;
    --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    background: var(--bg-page);
    color: var(--text-primary);
    font-family: var(--font-sans);
    font-size: 13px;
    line-height: 1.45;
    -webkit-font-smoothing: antialiased;
    overflow-x: hidden;
  }

  /* App Shell Layout */
  header {
    background: var(--bg-surface);
    border-bottom: 1px solid var(--border);
    padding: 0.65rem 1.25rem;
    display: flex;
    justify-content: space-between;
    align-items: center;
    position: sticky;
    top: 0;
    z-index: 100;
  }
  .brand { display: flex; align-items: center; gap: 0.65rem; }
  .brand-mark {
    width: 24px; height: 24px; background: var(--text-primary); color: #fff;
    border-radius: 4px; font-weight: 700; font-size: 11px;
    display: flex; align-items: center; justify-content: center;
    letter-spacing: 0.05em;
  }
  .brand-title { font-weight: 600; font-size: 14px; letter-spacing: -0.01em; }
  .brand-badge {
    font-size: 11px; background: var(--bg-page); border: 1px solid var(--border);
    padding: 0.15rem 0.5rem; border-radius: 3px; color: var(--text-secondary);
  }

  .header-actions { display: flex; align-items: center; gap: 0.75rem; }
  .search-bar { display: flex; align-items: center; gap: 0.4rem; }
  .search-input {
    width: 220px; padding: 0.35rem 0.6rem; font-size: 12px; font-family: var(--font-mono);
    border: 1px solid var(--border); border-radius: 4px; background: #fff; outline: none;
    text-transform: uppercase;
  }
  .search-input:focus { border-color: var(--accent); }
  .btn {
    background: var(--text-primary); color: #fff; border: 1px solid transparent; padding: 0.35rem 0.75rem;
    font-size: 12px; font-weight: 500; border-radius: 4px; cursor: pointer; transition: all 0.15s ease;
    display: inline-flex; align-items: center; gap: 0.35rem; text-decoration: none; user-select: none;
  }
  .btn:hover { background: #1e293b; }
  .btn-accent { background: var(--accent); }
  .btn-accent:hover { background: var(--accent-hover); }
  .btn-outline {
    background: #fff; color: var(--text-secondary); border: 1px solid var(--border);
  }
  .btn-outline:hover { background: var(--bg-page); color: var(--text-primary); border-color: #cbd5e1; }
  .btn-sm { padding: 0.25rem 0.5rem; font-size: 11px; }

  /* App Container: Sidebar + Main Content */
  .app-container {
    display: flex;
    min-height: calc(100vh - 49px);
  }

  /* Left Sidebar Navigation */
  .sidebar {
    width: var(--sidebar-w);
    background: var(--bg-surface);
    border-right: 1px solid var(--border);
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    flex-shrink: 0;
  }
  .sidebar-nav { padding: 1rem 0.75rem; display: flex; flex-direction: column; gap: 0.35rem; }
  .sidebar-section-title {
    font-size: 10px; font-weight: 700; color: var(--text-muted);
    text-transform: uppercase; letter-spacing: 0.05em; padding: 0.5rem 0.5rem 0.25rem 0.5rem;
  }
  .nav-btn {
    display: flex; align-items: center; gap: 0.65rem; width: 100%;
    padding: 0.6rem 0.75rem; border: none; background: transparent;
    border-radius: 6px; font-size: 12.5px; font-weight: 500;
    color: var(--text-secondary); cursor: pointer; text-align: left;
    transition: all 0.15s ease;
  }
  .nav-btn:hover { background: var(--bg-page); color: var(--text-primary); }
  .nav-btn.active {
    background: var(--accent-light); color: var(--accent); font-weight: 600;
  }
  .nav-bullet {
    font-size: 10.5px; font-family: var(--font-mono); font-weight: 600;
    width: 20px; height: 20px; border-radius: 4px; background: var(--border-light);
    color: var(--text-secondary); display: flex; align-items: center; justify-content: center;
  }
  .nav-btn.active .nav-bullet { background: var(--accent); color: #fff; }
  .nav-badge {
    margin-left: auto; font-size: 10px; padding: 0.1rem 0.4rem;
    border-radius: 3px; background: var(--border-light); color: var(--text-muted);
  }
  .nav-btn.active .nav-badge { background: #dbeafe; color: var(--accent); }

  .sidebar-footer {
    padding: 1rem 0.75rem; border-top: 1px solid var(--border);
    font-size: 11px; color: var(--text-muted);
  }
  .sidebar-stat-row { display: flex; justify-content: space-between; margin-bottom: 0.25rem; }

  /* Content Workspace */
  .workspace {
    flex: 1;
    overflow-y: auto;
    background: var(--bg-page);
    padding: 1.25rem;
  }

  /* View Panels */
  .view-panel { display: none; }
  .view-panel.active { display: block; animation: fadeIn 0.15s ease-in-out; }
  @keyframes fadeIn { from { opacity: 0; } to { opacity: 1; } }

  /* Subheader Controls */
  .control-strip {
    background: var(--bg-surface); border: 1px solid var(--border); border-radius: 6px;
    padding: 0.65rem 1rem; display: flex; align-items: center; justify-content: space-between;
    margin-bottom: 1.25rem; font-size: 12px; flex-wrap: wrap; gap: 0.75rem;
  }
  .control-group { display: flex; align-items: center; gap: 0.5rem; }
  .control-label { color: var(--text-secondary); font-weight: 500; font-size: 12px; white-space: nowrap; }
  select {
    padding: 0.3rem 0.6rem; font-size: 12px; border: 1px solid var(--border);
    border-radius: 4px; background: #fff; outline: none; color: var(--text-primary);
  }

  /* Layout Grids */
  .grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 1.25rem; margin-bottom: 1.25rem; }
  .grid-1 { display: grid; grid-template-columns: 1fr; gap: 1.25rem; margin-bottom: 1.25rem; }

  /* Cards */
  .card {
    background: var(--bg-surface); border: 1px solid var(--border);
    border-radius: 6px; display: flex; flex-direction: column; overflow: hidden;
  }
  .card-header {
    padding: 0.75rem 1rem; border-bottom: 1px solid var(--border);
    display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.5rem;
  }
  .card-title { font-size: 13px; font-weight: 600; color: var(--text-primary); }
  .card-caption { font-size: 11px; color: var(--text-muted); margin-top: 1px; }
  .card-body { padding: 1rem; position: relative; }

  /* KPI Summary Cards */
  .kpi-row { display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem; margin-bottom: 1rem; }
  .kpi-card {
    background: var(--bg-surface); border: 1px solid var(--border); border-radius: 6px;
    padding: 0.85rem 1.1rem; display: flex; flex-direction: column;
  }
  .kpi-label { font-size: 11px; font-weight: 600; color: var(--text-secondary); text-transform: uppercase; }
  .kpi-val { font-size: 22px; font-weight: 700; color: var(--text-primary); margin: 0.25rem 0; font-family: var(--font-mono); }
  .kpi-sub { font-size: 11px; color: var(--text-muted); }

  /* Purity Alert Banner */
  .purity-banner {
    background: #f8fafc;
    border: 1px solid #cbd5e1;
    border-left: 4px solid var(--accent);
    border-radius: 4px;
    padding: 0.55rem 0.85rem;
    font-size: 11.5px;
    color: var(--text-secondary);
    line-height: 1.45;
    margin-bottom: 1.25rem;
    display: flex;
    align-items: flex-start;
    gap: 0.6rem;
  }

  /* Preset Chips */
  .preset-strip {
    display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap; margin-bottom: 0.85rem;
  }
  .preset-label { font-size: 11px; font-weight: 600; color: var(--text-secondary); white-space: nowrap; }
  .chip {
    padding: 0.25rem 0.6rem; border-radius: 4px; font-size: 11.5px; font-weight: 500;
    cursor: pointer; border: 1px solid var(--border); background: var(--bg-surface);
    color: var(--text-secondary); transition: all 0.15s ease;
  }
  .chip:hover { border-color: var(--accent); color: var(--accent); background: var(--accent-light); }
  .chip.active { background: var(--accent); color: #fff; border-color: var(--accent); }

  /* Sub-Tab Navigation Strip */
  .sub-tab-bar {
    display: flex;
    gap: 0.5rem;
    border-bottom: 2px solid var(--border);
    margin-bottom: 1.25rem;
    background: var(--bg-surface);
    border-radius: 6px 6px 0 0;
    padding: 0.25rem 0.5rem 0 0.5rem;
    border: 1px solid var(--border);
    border-bottom: 2px solid var(--border);
  }
  .sub-tab {
    padding: 0.65rem 1.2rem;
    background: transparent;
    border: none;
    border-bottom: 3px solid transparent;
    margin-bottom: -2px;
    font-size: 13px;
    font-weight: 500;
    color: var(--text-secondary);
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    transition: all 0.15s ease;
  }
  .sub-tab:hover {
    color: var(--text-primary);
    background: var(--bg-page);
    border-radius: 4px 4px 0 0;
  }
  .sub-tab.active {
    color: var(--accent);
    border-bottom-color: var(--accent);
    font-weight: 600;
  }
  .sub-tab-panel { display: none; }
  .sub-tab-panel.active { display: block; animation: fadeIn 0.15s ease-in-out; }

  /* Segmented Toggle Control */
  .segmented-control {
    display: inline-flex;
    background: var(--border-light);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 2px;
    gap: 2px;
  }
  .seg-btn {
    padding: 0.25rem 0.65rem;
    font-size: 11.5px;
    font-weight: 500;
    border: none;
    background: transparent;
    color: var(--text-secondary);
    border-radius: 4px;
    cursor: pointer;
    transition: all 0.15s ease;
  }
  .seg-btn.active {
    background: #fff;
    color: var(--text-primary);
    font-weight: 600;
    box-shadow: 0 1px 2px rgba(0,0,0,0.05);
  }

  /* Multi-Select Dropdown */
  .dropdown-popover {
    display: none;
    position: absolute;
    top: calc(100% + 5px);
    left: 0;
    width: 320px;
    max-height: 400px;
    background: #ffffff;
    border: 1px solid var(--border);
    border-radius: 6px;
    box-shadow: 0 10px 25px -5px rgba(0,0,0,0.1), 0 8px 10px -6px rgba(0,0,0,0.1);
    z-index: 250;
    padding: 0.6rem;
    flex-direction: column;
  }
  .dropdown-popover.open { display: flex; animation: fadeIn 0.1s ease-out; }
  .tumor-item-row {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 0.3rem 0.5rem;
    border-radius: 4px;
    cursor: pointer;
    font-size: 12px;
  }
  .tumor-item-row:hover { background: var(--bg-page); }

  /* Export Button Action Bar */
  .export-actions { display: flex; align-items: center; gap: 0.35rem; }

  /* Tabular formatting */
  table { width: 100%; border-collapse: collapse; font-size: 12px; text-align: left; }
  th {
    padding: 0.55rem 0.75rem; border-bottom: 1px solid var(--border);
    color: var(--text-secondary); font-weight: 600; font-size: 11px; text-transform: uppercase;
    position: sticky; top: 0; background: var(--bg-surface); z-index: 10;
  }
  th.sortable { cursor: pointer; user-select: none; }
  th.sortable:hover { background: var(--border-light); color: var(--text-primary); }
  .sort-indicator { font-size: 10px; margin-left: 4px; color: var(--text-muted); display: inline-block; }
  td { padding: 0.55rem 0.75rem; border-bottom: 1px solid var(--border-light); }
  tr:last-child td { border-bottom: none; }
  tbody tr:hover { background: #fafafa; }
  .num { text-align: right; font-family: var(--font-mono); }

  /* Visualizers */
  .box-row { display: flex; align-items: center; gap: 0.75rem; margin-bottom: 0.45rem; }
  .box-label { width: 140px; font-weight: 500; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; font-size: 12px; }
  .box-canvas { flex: 1; height: 20px; position: relative; background: #f1f5f9; border-radius: 2px; }
  .box-iqr { position: absolute; height: 12px; top: 4px; background: #cbd5e1; border: 1px solid #94a3b8; border-radius: 2px; }
  .box-med { position: absolute; height: 16px; top: 2px; width: 2px; background: var(--text-primary); }
  .box-point {
    position: absolute; width: 6px; height: 6px; border-radius: 50%;
    background: var(--accent); top: 7px; opacity: 0.65; cursor: pointer;
  }
  .box-point:hover { opacity: 1; transform: scale(1.6); }

  .tag {
    display: inline-block; padding: 0.15rem 0.4rem; border-radius: 3px;
    font-size: 10.5px; font-weight: 500; font-family: var(--font-mono);
  }
  .tag-canon { background: #dbeafe; color: #1e40af; }
  .tag-alt { background: #f1f5f9; color: #475569; }
  .tag-fda { background: #dcfce7; color: #15803d; font-weight: 600; }
  .tag-p3 { background: #fef3c7; color: #b45309; font-weight: 600; }
  .tag-p2 { background: #e0e7ff; color: #3730a3; font-weight: 600; }
  .tag-clamped { background: #fee2e2; color: #b91c1c; font-weight: 600; }

  /* Bar Charts & Heatmap Containers */
  .chart-container { width: 100%; position: relative; }
  .heatmap-scroll-wrap {
    overflow-x: auto;
    overflow-y: auto;
    max-height: 580px;
    border: 1px solid var(--border);
    border-radius: 4px;
    background: #ffffff;
  }

  /* Tooltip for Heatmap */
  #heatmap-tooltip {
    position: fixed;
    display: none;
    background: rgba(15, 23, 42, 0.95);
    color: #fff;
    padding: 0.45rem 0.8rem;
    border-radius: 4px;
    font-size: 11px;
    line-height: 1.45;
    pointer-events: none;
    z-index: 1000;
    box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1), 0 2px 4px -1px rgba(0,0,0,0.06);
  }
</style>
</head>
<body>

<div id="heatmap-tooltip"></div>

<header>
  <div class="brand">
    <div class="brand-mark">RN</div>
    <div class="brand-title">Clinical RNA-Seq Studio</div>
    <span class="brand-badge" id="cohort-badge">Loading cohort...</span>
  </div>
</header>

<div class="app-container">
  <!-- LEFT SIDEBAR -->
  <aside class="sidebar">
    <div class="sidebar-nav">
      <div class="sidebar-section-title">Analysis Workflows</div>
      <button id="nav-single-gene" class="nav-btn active" onclick="switchView('single-gene')">
        <span class="nav-bullet">01</span>
        <span class="nav-text">Single Gene Profile</span>
        <span class="nav-badge">Boxplot</span>
      </button>
      <button id="nav-multi-gene" class="nav-btn" onclick="switchView('multi-gene')">
        <span class="nav-bullet">02</span>
        <span class="nav-text">Multi-Gene & ADCs</span>
        <span class="nav-badge">Co-Target</span>
      </button>
      <button id="nav-diff-expr" class="nav-btn" onclick="switchView('diff-expr')">
        <span class="nav-bullet">03</span>
        <span class="nav-text">Differential Contrast</span>
        <span class="nav-badge">Volcano</span>
      </button>
      <button id="nav-umap" class="nav-btn" onclick="switchView('umap')">
        <span class="nav-bullet">04</span>
        <span class="nav-text">Cohort Landscape</span>
        <span class="nav-badge">UMAP/t-SNE</span>
      </button>
      <button id="nav-tmb" class="nav-btn" onclick="switchView('tmb')">
        <span class="nav-bullet">05</span>
        <span class="nav-text">Mutation Burden</span>
        <span class="nav-badge">TMB</span>
      </button>
    </div>

    <div class="sidebar-footer">
      <div class="sidebar-stat-row"><span>Cohort:</span> <strong style="color:var(--text-primary);">HMF Pan-Cancer</strong></div>
      <div class="sidebar-stat-row"><span>Samples:</span> <span id="stat-samples" style="font-family:var(--font-mono);">1,952</span></div>
      <div class="sidebar-stat-row"><span>Tumor Types:</span> <span id="stat-cancers" style="font-family:var(--font-mono);">32</span></div>
      <div class="sidebar-stat-row"><span>Quant Genes:</span> <span id="stat-genes" style="font-family:var(--font-mono);">39,500</span></div>
    </div>
  </aside>

  <!-- MAIN WORKSPACE -->
  <main class="workspace">

    <!-- =================================================================== -->
    <!-- VIEW 1: SINGLE GENE PROFILE                                         -->
    <!-- =================================================================== -->
    <div id="view-single-gene" class="view-panel active">
      <div class="control-strip">
        <div class="control-group">
          <span class="control-label">Target Gene:</span>
          <input type="text" id="gene-search" class="search-input" style="width:130px;" value="EGFR">
          <button class="btn btn-sm btn-accent" onclick="loadGeneData()">Update</button>
        </div>
        <div class="control-group">
          <span class="control-label">Expression Metric:</span>
          <select id="metric-select" onchange="loadGeneData()">
            <option value="log2_adj_tpm" selected>log2(AdjTPM + 1)</option>
            <option value="adj_tpm">AdjTPM (Linear)</option>
            <option value="batch_corrected">ComBat Batch Adjusted</option>
            <option value="norm_counts">TMM Normalized Counts</option>
          </select>
        </div>
        <div class="control-group">
          <span class="control-label">Benchmark Cutoff:</span>
          <select id="cutoff-select" onchange="loadGeneData()">
            <option value="90" selected>&ge; 90th Percentile</option>
            <option value="95">&ge; 95th Percentile</option>
            <option value="75">&ge; 75th Percentile</option>
          </select>
        </div>
      </div>

      <div class="grid-2">
        <!-- PANE 1: Cross-Cancer Expression Distribution -->
        <div class="card">
          <div class="card-header">
            <div>
              <div class="card-title">1. Cross-Cancer Gene Distribution: <span id="title-gene" style="color:var(--accent);">EGFR</span></div>
              <div class="card-caption">Boxplot & median log2 expression distributions across cancer types</div>
            </div>
            <div class="export-actions">
              <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('threshold-table', 'single_gene_threshold_matrix.csv')">Export CSV</button>
            </div>
          </div>
          <div class="card-body" id="dist-container" style="max-height:540px; overflow-y:auto;">
            <div style="color:var(--text-muted); padding:2rem; text-align:center;">Loading distribution...</div>
          </div>
        </div>

        <!-- PANE 2: Over-Expression Threshold Matrix -->
        <div class="card">
          <div class="card-header">
            <div>
              <div class="card-title">2. Over-Expression vs. HMF Reference Benchmark</div>
              <div class="card-caption">Patient count & fraction exceeding reference baseline (click column to sort)</div>
            </div>
            <div class="export-actions">
              <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('threshold-table', 'single_gene_threshold_matrix.csv')">Export CSV</button>
              <button class="btn btn-outline btn-sm" onclick="exportTableToPNG('threshold-table', 'single_gene_benchmark_' + currentActiveGene + '.png')">Export PNG</button>
              <button class="btn btn-outline btn-sm" onclick="exportStandaloneHtml('Single Gene Benchmark Matrix: ' + currentActiveGene, 'threshold-table-wrap', 'single_gene_benchmark_' + currentActiveGene + '.html')">Export Sortable HTML</button>
            </div>
          </div>
          <div class="card-body" id="threshold-table-wrap" style="padding:0; max-height:540px; overflow-y:auto;">
            <table id="threshold-table" class="sortable">
              <thead>
                <tr>
                  <th class="sortable">Cancer Type</th>
                  <th class="num sortable">Patients Above Cutoff</th>
                  <th class="num sortable">Total N</th>
                  <th class="num sortable">% Elevated</th>
                  <th style="width:90px;">Proportion</th>
                </tr>
              </thead>
              <tbody id="threshold-body">
                <tr><td colspan="5" style="text-align:center; color:var(--text-muted); padding:2rem;">Loading matrix...</td></tr>
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- VIEW 2: MULTI-GENE & ADC STUDIO                                     -->
    <!-- =================================================================== -->
    <div id="view-multi-gene" class="view-panel">
      <!-- Preset Selection Chips -->
      <div class="preset-strip">
        <span class="preset-label">Curated ADC Panels:</span>
        <button class="chip active" onclick="loadAdcPreset('horizon')">Horizon-Beacon (5)</button>
        <button class="chip" onclick="loadAdcPreset('fda')">FDA-Approved (4)</button>
        <button class="chip" onclick="loadAdcPreset('pipeline')">Phase 3 Pipeline (4)</button>
        <button class="chip" onclick="loadAdcPreset('thoracic')">Thoracic / NSCLC (5)</button>
        <button class="chip" onclick="loadAdcPreset('gi')">GI & Pancreatic (4)</button>
        <button class="chip" onclick="loadAdcPreset('breast')">Breast & GYN (4)</button>
        <button class="chip" onclick="loadAdcPreset('gu')">GU & Prostate (4)</button>
        <button class="chip" onclick="loadAdcPreset('neuro')">SCLC & Neuro (3)</button>
        <button class="chip" onclick="loadAdcPreset('all')">All 12 Curated ADCs</button>
      </div>

      <!-- Main Controls Strip: Genes + Tumor Cohort + Purity Weighting + Purity Filter + Threshold -->
      <div class="control-strip" style="align-items:flex-end;">
        <!-- Gene Symbols -->
        <div class="control-group" style="flex:1; min-width:240px;">
          <div style="display:flex; flex-direction:column; gap:3px; width:100%;">
            <span class="control-label">Gene Symbols:</span>
            <input type="text" id="multi-gene-input" class="search-input" style="width:100%;" value="ERBB2, TACSTD2, CLDN18, FOLR1, VTCN1" placeholder="Comma or space separated (e.g. ERBB2, TACSTD2, CLDN18, FOLR1, VTCN1)">
          </div>
        </div>

        <!-- Tumor Types Selection Popover -->
        <div class="control-group" style="position:relative;">
          <div style="display:flex; flex-direction:column; gap:3px;">
            <span class="control-label">Tumor Cohorts:</span>
            <button id="tumor-filter-btn" class="btn btn-outline" style="min-width:170px; justify-content:space-between;" onclick="toggleTumorDropdown(event)">
              <span id="tumor-filter-label">All Tumor Types (32)</span>
              <span style="font-size:10px; margin-left:6px; color:var(--text-muted);">&#9662;</span>
            </button>
          </div>

          <!-- Dropdown Popover Menu -->
          <div id="tumor-dropdown-menu" class="dropdown-popover">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.4rem; padding-bottom:0.35rem; border-bottom:1px solid var(--border-light);">
              <span style="font-size:11px; font-weight:600; color:var(--text-secondary);">Filter Tumor Cohorts</span>
              <div style="display:flex; gap:0.25rem;">
                <button class="btn btn-sm btn-outline" style="padding:0.15rem 0.4rem; font-size:10.5px;" onclick="selectAllTumors(true)">Select All</button>
                <button class="btn btn-sm btn-outline" style="padding:0.15rem 0.4rem; font-size:10.5px;" onclick="selectAllTumors(false)">Clear</button>
              </div>
            </div>

            <!-- Quick tumor cohort chips -->
            <div style="display:flex; gap:4px; flex-wrap:wrap; margin-bottom:0.4rem;">
              <span class="chip" style="font-size:10px; padding:0.15rem 0.35rem;" onclick="setQuickTumorCohort(['Breast'])">Breast</span>
              <span class="chip" style="font-size:10px; padding:0.15rem 0.35rem;" onclick="setQuickTumorCohort(['Colorectum', 'Pancreas', 'Esophagus', 'Stomach'])">GI/Panc</span>
              <span class="chip" style="font-size:10px; padding:0.15rem 0.35rem;" onclick="setQuickTumorCohort(['Lung'])">Lung</span>
              <span class="chip" style="font-size:10px; padding:0.15rem 0.35rem;" onclick="setQuickTumorCohort(['Ovary'])">Ovary</span>
              <span class="chip" style="font-size:10px; padding:0.15rem 0.35rem;" onclick="setQuickTumorCohort(['Prostate'])">Prostate</span>
              <span class="chip" style="font-size:10px; padding:0.15rem 0.35rem;" onclick="setQuickTumorCohort(['Urinary tract'])">Bladder</span>
            </div>

            <input type="text" id="tumor-search-input" placeholder="Search tumor types..." style="padding:0.3rem 0.5rem; font-size:11px; border:1px solid var(--border); border-radius:4px; margin-bottom:0.4rem; outline:none;" oninput="filterTumorCheckboxList()">

            <div id="tumor-checkbox-list" style="overflow-y:auto; flex:1; max-height:200px; display:flex; flex-direction:column; gap:2px;">
              <!-- Populated dynamically -->
            </div>

            <div style="margin-top:0.4rem; padding-top:0.35rem; border-top:1px solid var(--border-light); display:flex; justify-content:space-between; align-items:center;">
              <span id="tumor-selected-count-label" style="font-size:10.5px; color:var(--text-muted);">32 cohorts selected</span>
              <button class="btn btn-sm btn-accent" onclick="applyTumorSelection()">Apply Filter</button>
            </div>
          </div>
        </div>

        <!-- Expression Weighting Model (Purity Adjusted vs Bulk vs TME) -->
        <div class="control-group">
          <div style="display:flex; flex-direction:column; gap:3px;">
            <span class="control-label" title="Purity-weighted cellular deconvolution">Cellular Model:</span>
            <select id="multi-weight-select" onchange="runMultiGeneAnalysis()">
              <option value="bulk" selected>Bulk Expression (Standard)</option>
              <option value="tumor_adjusted">Tumor-Cell Adjusted (TPM / Purity)</option>
              <option value="tme_weighted">TME-Weighted (TPM &times; (1 - Purity))</option>
            </select>
          </div>
        </div>

        <!-- PURPLE Purity Threshold Filter -->
        <div class="control-group">
          <div style="display:flex; flex-direction:column; gap:3px;">
            <span class="control-label" title="Exclude low cellularity / stromal biopsies">Purity Filter:</span>
            <select id="multi-purity-filter" onchange="runMultiGeneAnalysis()">
              <option value="0.0" selected>All Purity (No Filter)</option>
              <option value="0.20">&ge; 20% Purity (Recommended)</option>
              <option value="0.30">&ge; 30% Purity (High cellularity)</option>
              <option value="0.50">&ge; 50% Purity (Pure tumor core)</option>
            </select>
          </div>
        </div>

        <!-- Percentile Cutoff -->
        <div class="control-group">
          <div style="display:flex; flex-direction:column; gap:3px;">
            <span class="control-label">Cutoff:</span>
            <select id="multi-threshold-select" onchange="runMultiGeneAnalysis()">
              <option value="90" selected>&ge; 90th %ile</option>
              <option value="95">&ge; 95th %ile</option>
              <option value="75">&ge; 75th %ile</option>
            </select>
          </div>
        </div>

        <button class="btn btn-accent" style="height:28px;" onclick="runMultiGeneAnalysis()">Analyze</button>
      </div>

      <!-- Purity Floor & Model Annotation Banner -->
      <div id="purity-annotation-banner" class="purity-banner" style="display:none;">
        <div id="purity-banner-text">Loading purity annotations...</div>
      </div>

      <!-- KPI Summary Cards -->
      <div class="kpi-row">
        <div class="kpi-card">
          <div class="kpi-label">Evaluated Cohort</div>
          <div class="kpi-val" id="kpi-total-cohort">1,952</div>
          <div class="kpi-sub" id="kpi-cohort-sub">Total patients across 32 tumor types</div>
        </div>
        <div class="kpi-card" style="border-left: 3px solid var(--accent);">
          <div class="kpi-label">&ge; 1 Target Elevated</div>
          <div class="kpi-val" id="kpi-any-elevated" style="color:var(--accent);">-</div>
          <div class="kpi-sub" id="kpi-any-sub">Patients eligible for at least one candidate target</div>
        </div>
        <div class="kpi-card" style="border-left: 3px solid var(--up-reg);">
          <div class="kpi-label">&ge; 2 Targets Co-Elevated</div>
          <div class="kpi-val" id="kpi-co-elevated" style="color:var(--up-reg);">-</div>
          <div class="kpi-sub" id="kpi-co-sub">Candidates for dual-targeting or sequential therapy</div>
        </div>
      </div>

      <!-- =============================================================== -->
      <!-- TABS: Target Over Expression Summary vs Cross-Cancer Matrix     -->
      <!-- =============================================================== -->
      <div class="sub-tab-bar">
        <button id="tab-btn-summary" class="sub-tab active" onclick="switchMultiSubTab('summary')">
          Target Over-Expression Summary
        </button>
        <button id="tab-btn-matrix" class="sub-tab" onclick="switchMultiSubTab('matrix')">
          Cross-Cancer Target Matrix
        </button>
        <button id="tab-btn-catalog" class="sub-tab" onclick="switchMultiSubTab('catalog')">
          Curated ADC Target Knowledge Base
        </button>
      </div>

      <!-- SUB-TAB 1: TARGET OVER-EXPRESSION SUMMARY -->
      <div id="multi-subtab-summary" class="sub-tab-panel active">
        <!-- SVG Horizontal Bar Chart -->
        <div class="grid-1">
          <div class="card">
            <div class="card-header">
              <div>
                <div class="card-title">Elevated Patient Counts by Gene Target</div>
                <div id="multi-chart-sub" class="card-caption">Number and fraction of patients exceeding expression benchmark in selected cohort</div>
              </div>
              <div style="display:flex; align-items:center; gap:0.5rem;">
                <div id="multi-chart-cutoff-badge" class="brand-badge">Cutoff: &ge; 90th Percentile</div>
                <div class="export-actions">
                  <button class="btn btn-outline btn-sm" onclick="exportSvg('#multi-bar-container svg', 'elevated_patient_counts_plot.svg')">Export SVG</button>
                  <button class="btn btn-outline btn-sm" onclick="exportSvgToPng('#multi-bar-container svg', 'elevated_patient_counts_plot.png')">Export PNG</button>
                </div>
              </div>
            </div>
            <div class="card-body">
              <div id="multi-bar-container" class="chart-container">
                <div style="color:var(--text-muted); padding:2rem; text-align:center;">Click "Analyze" to generate plot.</div>
              </div>
            </div>
          </div>
        </div>

        <!-- Sortable Summary Table -->
        <div class="grid-1">
          <div class="card">
            <div class="card-header">
              <div>
                <div class="card-title">Target Over-Expression Prevalence Table</div>
                <div class="card-caption">Ranked cohort statistics across evaluated genes (click any column header to sort)</div>
              </div>
              <div class="export-actions">
                <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('multi-gene-summary-table', 'target_overexpression_summary.csv')">Export CSV</button>
                <button class="btn btn-outline btn-sm" onclick="exportTableToPNG('multi-gene-summary-table', 'target_overexpression_summary.png')">Export PNG</button>
                <button class="btn btn-outline btn-sm" onclick="exportStandaloneHtml('Target Over-Expression Summary', 'summary-table-wrap', 'target_overexpression_summary.html')">Export Sortable HTML</button>
              </div>
            </div>
            <div class="card-body" id="summary-table-wrap" style="padding:0; max-height:460px; overflow-y:auto;">
              <table id="multi-gene-summary-table" class="sortable">
                <thead>
                  <tr>
                    <th class="sortable">Gene</th>
                    <th class="num sortable">Elevated Patients</th>
                    <th class="num sortable">Cohort Evaluated N</th>
                    <th class="num sortable">% Elevated</th>
                    <th class="num sortable" title="Samples with PURPLE purity < 15% clamped to floor">Below Floor (&lt;15%)*</th>
                    <th style="width:80px; text-align:center;">Action</th>
                  </tr>
                </thead>
                <tbody id="multi-gene-summary-body">
                  <tr><td colspan="6" style="text-align:center; color:var(--text-muted); padding:2rem;">Run analysis above to populate.</td></tr>
                </tbody>
              </table>
              <div style="font-size:10.5px; color:var(--text-muted); padding:0.5rem 0.75rem; border-top:1px solid var(--border-light);">
                * <strong>Purity Floor: 15% (0.15)</strong> is applied in Tumor-Adjusted mode. Any biopsy with PURPLE tumor cellularity &lt; 15% is clamped to 15% to eliminate division-by-zero noise amplification.
              </div>
            </div>
          </div>
        </div>
      </div>

      <!-- SUB-TAB 2: CROSS-CANCER TARGET MATRIX (HEATMAP / TABLE) -->
      <div id="multi-subtab-matrix" class="sub-tab-panel">
        <div class="card">
          <div class="card-header">
            <div>
              <div class="card-title">Cross-Cancer Target Over-Expression Matrix</div>
              <div class="card-caption">Tumor cohort breakdown across evaluated targets</div>
            </div>
            <!-- Heatmap / Table Segmented Control & Export Actions -->
            <div style="display:flex; align-items:center; gap:0.75rem; flex-wrap:wrap;">
              <div id="heatmap-legend-bar" style="display:flex; align-items:center; gap:6px; font-size:11px; color:var(--text-secondary); margin-right:0.4rem;">
                <span style="font-size:10px; font-weight:600; color:var(--text-muted);">Prevalence:</span>
                <span style="font-size:10.5px; font-weight:500;">0%</span>
                <div style="width:85px; height:10px; border-radius:2px; border:1px solid #cbd5e1; background:linear-gradient(to right, #ffffff 0%, #fee0d2 20%, #fc9272 45%, #de2d26 70%, #67000d 100%);"></div>
                <span style="font-size:10.5px; font-weight:500;">100%</span>
              </div>
              <div class="segmented-control">
                <button id="btn-seg-heatmap" class="seg-btn active" onclick="setMatrixView('heatmap')">Heatmap View</button>
                <button id="btn-seg-table" class="seg-btn" onclick="setMatrixView('table')">Table View</button>
              </div>
              <div class="export-actions">
                <button id="btn-export-matrix-svg" class="btn btn-outline btn-sm" onclick="exportSvg('#heatmap-container svg', 'cross_cancer_target_heatmap.svg')">Export SVG</button>
                <button id="btn-export-matrix-png" class="btn btn-outline btn-sm" onclick="exportMatrixAsPng()">Export PNG</button>
                <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('multi-cancer-matrix-table', 'cross_cancer_target_matrix.csv')">Export CSV</button>
                <button class="btn btn-outline btn-sm" onclick="exportMatrixStandaloneHtml()">Export Sortable HTML</button>
              </div>
            </div>
          </div>

          <div class="card-body" style="padding:1rem;">
            <!-- 1. Heatmap View Container -->
            <div id="matrix-heatmap-view">
              <!-- Heatmap Legend & Summary Info -->
              <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.75rem; flex-wrap:wrap; gap:0.5rem;">
                <div style="font-size:11.5px; color:var(--text-secondary);">
                  Cell values indicate <strong>% of patients</strong> in each tumor type elevated above cutoff. Hover over any cell for exact patient counts and purity details.
                </div>
                <!-- Legend Gradient -->
                <div style="display:flex; align-items:center; gap:6px; font-size:11px; font-family:var(--font-mono);">
                  <span style="color:var(--text-muted);">0%</span>
                  <div style="width:110px; height:12px; border-radius:3px; background:linear-gradient(to right, #f1f5f9, #93c5fd, #2563eb, #e11d48); border:1px solid #cbd5e1;"></div>
                  <span style="color:var(--text-muted);">&ge; 50%</span>
                </div>
              </div>

              <!-- Interactive SVG Heatmap -->
              <div class="heatmap-scroll-wrap" id="heatmap-container">
                <div style="color:var(--text-muted); padding:3rem; text-align:center;">Run analysis above to generate heatmap matrix.</div>
              </div>
            </div>

            <!-- 2. Sortable Table View Container -->
            <div id="matrix-table-view" style="display:none; max-height:560px; overflow-y:auto;">
              <table id="multi-cancer-matrix-table" class="sortable">
                <thead id="multi-cancer-head">
                  <tr>
                    <th class="sortable">Cancer Type</th>
                    <th class="num sortable">Patients</th>
                    <th class="num sortable">Any Target N</th>
                    <th class="num sortable">Any Target %</th>
                    <th class="num sortable" title="Biopsies with PURPLE tumor purity < 15%">Below Floor (&lt;15%)</th>
                  </tr>
                </thead>
                <tbody id="multi-cancer-body">
                  <tr><td colspan="5" style="text-align:center; color:var(--text-muted); padding:2rem;">Run analysis above to populate.</td></tr>
                </tbody>
              </table>
            </div>
          </div>
        </div>
      </div>

      <!-- SUB-TAB 3: CURATED ADC KNOWLEDGE BASE -->
      <div id="multi-subtab-catalog" class="sub-tab-panel">
        <div class="card">
          <div class="card-header">
            <div>
              <div class="card-title">Curated ADC Target & Drug Knowledge Base</div>
              <div class="card-caption">Clinically validated and late-stage antibody-drug conjugate (ADC) targets with approved/investigational agents</div>
            </div>
            <div class="export-actions">
              <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('adc-catalog-table', 'curated_adc_knowledge_base.csv')">Export CSV</button>
              <button class="btn btn-outline btn-sm" onclick="exportTableToPNG('adc-catalog-table', 'curated_adc_knowledge_base.png')">Export PNG</button>
              <button class="btn btn-outline btn-sm" onclick="exportStandaloneHtml('Curated ADC Knowledge Base', 'adc-catalog-wrap', 'curated_adc_knowledge_base.html')">Export Sortable HTML</button>
            </div>
          </div>
          <div class="card-body" id="adc-catalog-wrap" style="padding:0; max-height:560px; overflow-y:auto;">
            <table id="adc-catalog-table" class="sortable">
              <thead>
                <tr>
                  <th class="sortable">Target Gene</th>
                  <th class="sortable">Antigen / Receptor</th>
                  <th class="sortable">Approved & Pipeline Drug Conjugates</th>
                  <th class="sortable">Payload Class</th>
                  <th class="sortable">Key Indications</th>
                  <th class="sortable">Status</th>
                  <th style="width:130px; text-align:center;">Actions</th>
                </tr>
              </thead>
              <tbody id="adc-catalog-body">
                <!-- Dynamically populated from curated ADC library -->
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- VIEW 3: DIFFERENTIAL CONTRAST STUDIO                                -->
    <!-- =================================================================== -->
    <div id="view-diff-expr" class="view-panel">
      <div class="control-strip">
        <div class="control-group">
          <span class="control-label">Clinical Group Contrast:</span>
          <select id="contrast-select">
            <option value="ecdna_status:Yes:No" selected>ecDNA Status (Yes vs No)</option>
            <option value="msi_status:MSI:MSS">MSI Status (MSI vs MSS)</option>
          </select>
          <button id="run-contrast-btn" class="btn btn-accent btn-sm" onclick="runDifferentialExpression()">Run Contrast</button>
        </div>
        <div style="font-size:11px; color:var(--text-muted);">Genome-wide Welch's t-test with in-engine SQL aggregation</div>
      </div>

      <div class="grid-2">
        <!-- Left: Volcano Plot -->
        <div class="card">
          <div class="card-header">
            <div>
              <div class="card-title">Genome-Wide Volcano Plot: <span id="contrast-title" style="color:var(--accent);">ecDNA (Yes vs No)</span></div>
              <div class="card-caption">Log2 Fold Change vs -Log10 P-value (click any dot to profile gene in Studio)</div>
            </div>
            <div style="display:flex; align-items:center; gap:0.8rem; flex-wrap:wrap;">
              <div style="font-size:11px; display:flex; gap:0.6rem; color:var(--text-secondary);">
                <span style="display:flex; align-items:center; gap:4px;"><span style="width:8px; height:8px; border-radius:50%; background:#2563eb; display:inline-block;"></span> Down (FC &lt; -1)</span>
                <span style="display:flex; align-items:center; gap:4px;"><span style="width:8px; height:8px; border-radius:50%; background:#dc2626; display:inline-block;"></span> Up (FC &gt; +1)</span>
                <span style="display:flex; align-items:center; gap:4px;"><span style="width:8px; height:8px; border-radius:50%; background:#94a3b8; display:inline-block;"></span> Not Sig</span>
              </div>
              <div class="export-actions">
                <button class="btn btn-outline btn-sm" onclick="exportSvg('#volcano-svg', 'volcano_plot.svg')">Export SVG</button>
                <button class="btn btn-outline btn-sm" onclick="exportSvgToPng('#volcano-svg', 'volcano_plot.png')">Export PNG</button>
                <button class="btn btn-outline btn-sm" onclick="exportVolcanoStandaloneHtml()">Export Standalone HTML</button>
              </div>
            </div>
          </div>
          <div class="card-body">
            <svg id="volcano-svg" viewBox="0 0 580 340" style="width:100%; height:auto; max-height:360px; min-height:280px; background:#ffffff; border:1px solid var(--border); border-radius:6px; display:block;">
              <text x="290" y="170" text-anchor="middle" fill="#64748b" font-size="13" font-family="system-ui, -apple-system, sans-serif">Select contrast above and click &quot;Run Contrast&quot;</text>
            </svg>
          </div>
        </div>

        <!-- Right: Top DE Genes Table -->
        <div class="card">
          <div class="card-header">
            <div>
              <div class="card-title">Top Differentially Expressed Genes</div>
              <div class="card-caption">Ranked by Welch's t-test significance (click column to sort)</div>
            </div>
            <div class="export-actions">
              <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('diff-table', 'differential_expression_genes.csv')">Export CSV</button>
              <button class="btn btn-outline btn-sm" onclick="exportTableToPNG('diff-table', 'differential_expression_genes.png')">Export PNG</button>
              <button class="btn btn-outline btn-sm" onclick="exportStandaloneHtml('Top Differentially Expressed Genes', 'diff-table', 'differential_expression_genes_sortable.html')">Export Sortable HTML</button>
            </div>
          </div>
          <div class="card-body" style="padding:0; max-height:340px; overflow-y:auto;">
            <table id="diff-table" class="sortable">
              <thead>
                <tr>
                  <th class="sortable">Gene</th>
                  <th class="num sortable">Log2 FC</th>
                  <th class="num sortable">Mean Grp A</th>
                  <th class="num sortable">Mean Grp B</th>
                  <th class="num sortable">P-value</th>
                  <th class="num sortable">FDR</th>
                </tr>
              </thead>
              <tbody id="diff-body">
                <tr><td colspan="6" style="text-align:center; color:var(--text-muted); padding:2rem;">Click <strong>Run Contrast</strong> above to compute genome-wide fold changes.</td></tr>
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- VIEW 4: COHORT LANDSCAPE (UMAP / t-SNE / PCA)                       -->
    <!-- =================================================================== -->
    <div id="view-umap" class="view-panel">
      <!-- Control Strip -->
      <div class="control-strip" style="align-items:flex-end; gap:0.9rem;">
        <div class="control-group">
          <span class="control-label">Projection Algorithm:</span>
          <select id="umap-projection" onchange="renderUmapPlot()" style="min-width:230px; font-weight:600;">
            <option value="tsne" selected>t-SNE (Separated Clusters)</option>
            <option value="pca12">PCA (PC1 vs PC2, 34% Variance)</option>
            <option value="pca13">PCA (PC1 vs PC3, 29% Variance)</option>
            <option value="umap">UMAP (Manifold)</option>
          </select>
        </div>

        <div class="control-group">
          <span class="control-label">Coloring Dimension:</span>
          <select id="umap-color-by" onchange="renderUmapPlot()" style="min-width:270px; font-weight:600;">
            <optgroup label="Clinical & Classification">
              <option value="cancer_type" selected>1. Clinical Primary Tumour Type</option>
              <option value="cuppa">2. CUPPA Prediction</option>
              <option value="purity">3. Tumour Purity (White to Deep Blue)</option>
              <option value="specimen_type">4. Primary vs Metastasis vs Local Recurrence</option>
            </optgroup>
            <optgroup label="Horizon-Beacon Target Expressions (Log2 Adj TPM)">
              <option value="gene:ERBB2">5. Expression: ERBB2</option>
              <option value="gene:TACSTD2">5. Expression: TACSTD2</option>
              <option value="gene:CLDN18">5. Expression: CLDN18</option>
              <option value="gene:FOLR1">5. Expression: FOLR1</option>
              <option value="gene:VTCN1">5. Expression: VTCN1</option>
            </optgroup>
          </select>
        </div>

        <div class="control-group">
          <span class="control-label">Highlight Specific Sample:</span>
          <input type="text" id="umap-sample-search" class="search-input" placeholder="Search Sample ID..." oninput="highlightUmapSample()" style="width:190px;">
        </div>

        <div class="control-group" style="align-self:center; margin-top:14px;">
          <label style="display:flex; align-items:center; gap:6px; font-size:11.5px; cursor:pointer; user-select:none; color:var(--text-secondary);">
            <input type="checkbox" id="umap-floor-toggle" onchange="renderUmapPlot()" checked>
            <span>Annotate &lt; 15% Purity Floor</span>
          </label>
        </div>

        <div style="flex:1;"></div>

        <div class="export-actions">
          <button class="btn btn-outline btn-sm" onclick="exportSvg('#umap-svg', 'cohort_umap_landscape.svg')">Export SVG</button>
          <button class="btn btn-outline btn-sm" onclick="exportSvgToPng('#umap-svg', 'cohort_umap_landscape.png')">Export PNG</button>
          <button class="btn btn-outline btn-sm" onclick="exportUmapStandaloneHtml()">Export Standalone HTML</button>
        </div>
      </div>

      <!-- Main UMAP Layout: Canvas + Side Legend Panel -->
      <div style="display:grid; grid-template-columns: 1fr 320px; gap:1.25rem; align-items:start;">
        <!-- Left: UMAP Canvas Card -->
        <div class="card" style="padding:1rem; position:relative; overflow:hidden;">
          <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.6rem;">
            <div>
              <div class="card-title" id="umap-card-title">UMAP Embedding: Clinical Primary Tumour Type</div>
              <div class="card-caption" id="umap-card-caption">2D projection of 1,952 whole-transcriptome biopsies (Top 1,500 highly variable genes)</div>
            </div>
            <div id="umap-stats-badge" class="tag tag-p2" style="font-size:11px;">1,952 Biopsies</div>
          </div>

          <div id="umap-container" style="width:100%; height:580px; position:relative; background:#ffffff; border:1px solid var(--border); border-radius:6px; overflow:hidden;">
            <svg id="umap-svg" viewBox="0 0 740 560" style="width:100%; height:100%; display:block;">
              <text x="370" y="280" text-anchor="middle" fill="#64748b" font-size="13" font-family="system-ui, -apple-system, sans-serif">Loading cohort UMAP coordinates...</text>
            </svg>
            <div id="umap-tooltip" style="position:absolute; display:none; pointer-events:none; background:#0f172a; color:#f8fafc; padding:8px 12px; border-radius:6px; font-size:11.5px; box-shadow:0 10px 15px -3px rgba(0,0,0,0.3); z-index:100; max-width:300px; line-height:1.45;"></div>
          </div>
        </div>

        <!-- Right: Dynamic Legend & Distribution Panel -->
        <div class="card" style="padding:1rem; max-height:640px; display:flex; flex-direction:column;">
          <div class="card-title" id="umap-legend-title" style="margin-bottom:0.3rem;">Legend</div>
          <div class="card-caption" id="umap-legend-caption" style="margin-bottom:0.75rem;">Click any category to filter / highlight</div>

          <!-- Colorbar container (for continuous Purity and Gene Expression) -->
          <div id="umap-colorbar-container" style="display:none; padding:0.5rem 0; margin-bottom:0.75rem;">
            <div style="display:flex; justify-content:space-between; font-size:11px; margin-bottom:4px; font-weight:600;" id="umap-colorbar-labels">
              <span>0.15 (Floor)</span>
              <span>1.00 (Pure)</span>
            </div>
            <div id="umap-colorbar-gradient" style="width:100%; height:14px; border-radius:3px; border:1px solid #cbd5e1;"></div>
            <div id="umap-colorbar-note" style="font-size:10.5px; color:var(--text-muted); margin-top:6px; line-height:1.4;"></div>
          </div>

          <!-- Discrete category list (for Cancer Type, CUPPA, Specimen Type) -->
          <div id="umap-category-list" style="overflow-y:auto; flex:1; display:flex; flex-direction:column; gap:4px; padding-right:4px;">
            <!-- Rendered dynamically -->
          </div>
        </div>
      </div>
    </div>

    <!-- =================================================================== -->
    <!-- VIEW 6: TUMOUR MUTATION BURDEN (TMB) & MUTATION COUNTS               -->
    <!-- =================================================================== -->
    <div id="view-tmb" class="view-panel">
      <!-- Control Strip -->
      <div class="control-strip" style="align-items:flex-end; gap:0.9rem; flex-wrap:wrap;">
        <div class="control-group">
          <span class="control-label">Stratification / Grouping:</span>
          <select id="tmb-group-by" onchange="renderTmbPlot()" style="min-width:260px; font-weight:600;">
            <option value="cancer_type" selected>Clinical Cancer Type (Pan-Cancer Cohorts)</option>
            <option value="pole">POLE Signature (MT vs WT by Cancer Type)</option>
            <option value="msi">MSI Status (MSI vs MSS by Cancer Type)</option>
            <option value="cuppa">CUPPA Predicted Lineage</option>
            <option value="specimen">Specimen Type (Primary vs Metastasis)</option>
          </select>
        </div>

        <div class="control-group">
          <span class="control-label">Burden Metric:</span>
          <select id="tmb-metric" onchange="renderTmbPlot()" style="min-width:210px; font-weight:600;">
            <option value="tmb" selected>TMB (Mutations / Mb)</option>
            <option value="mut_count">Mutations per Sample (WGS Total)</option>
            <option value="coding_mut">Coding Mutations (Exome / TML)</option>
          </select>
        </div>

        <div class="control-group">
          <span class="control-label">Y-Axis Scale:</span>
          <select id="tmb-scale" onchange="renderTmbPlot()" style="width:170px; font-weight:600;">
            <option value="log10" selected>log10 Scale (Linear Axis)</option>
            <option value="log2">log2 Scale</option>
            <option value="linear">Linear Scale</option>
          </select>
        </div>

        <div class="control-group">
          <span class="control-label">Tumour Purity Floor:</span>
          <div style="display:flex; align-items:center; gap:6px;">
            <input type="number" id="tmb-purity-floor" min="0" max="0.9" step="0.05" value="0.15" onchange="renderTmbPlot()" style="width:65px; padding:4px 6px; border:1px solid #cbd5e1; border-radius:4px; font-size:12px; font-family:var(--font-mono); font-weight:600;">
            <label style="display:flex; align-items:center; gap:4px; font-size:11px; cursor:pointer; user-select:none; color:var(--text-secondary);">
              <input type="checkbox" id="tmb-purity-exclude" onchange="renderTmbPlot()" checked>
              <span>Exclude &lt; Floor</span>
            </label>
          </div>
        </div>

        <div class="control-group">
          <span class="control-label">Mutation Cutoff:</span>
          <div style="display:flex; align-items:center; gap:5px;">
            <input type="number" id="tmb-mut-threshold" value="4000" min="0" step="500" onchange="renderTmbPlot()" style="width:85px; padding:4px 6px; border:1px solid #cbd5e1; border-radius:4px; font-size:12px; font-family:var(--font-mono); font-weight:700;">
            <span style="font-size:11px; color:var(--text-muted); font-weight:500;">muts</span>
          </div>
        </div>

        <div class="control-group" style="align-self:center; margin-top:14px;">
          <label style="display:flex; align-items:center; gap:5px; font-size:11px; cursor:pointer; user-select:none; color:var(--text-secondary);">
            <input type="checkbox" id="tmb-show-threshold" onchange="renderTmbPlot()" checked>
            <span>Show Cutoff Benchmark Line</span>
          </label>
        </div>

        <div style="flex:1;"></div>

        <div style="display:flex; gap:0.5rem; align-items:center; flex-wrap:wrap;">
          <button class="btn btn-secondary" onclick="exportTmbPlotSvg()" style="padding:0.4rem 0.75rem; font-size:0.75rem;">
            Export Plot (SVG)
          </button>
          <button class="btn btn-secondary" onclick="exportTmbPlotPng()" style="padding:0.4rem 0.75rem; font-size:0.75rem;">
            Export Plot (PNG)
          </button>
          <button class="btn btn-secondary" onclick="exportTmbStandaloneHtml()" style="padding:0.4rem 0.75rem; font-size:0.75rem;">
            Export Standalone HTML
          </button>
          <button class="btn btn-secondary" onclick="exportTmbTableCsv()" style="padding:0.4rem 0.75rem; font-size:0.75rem;">
            Export Table (CSV)
          </button>
          <button class="btn btn-secondary" onclick="exportTableToPNG('tmb-summary-table', 'cohort_tmb_summary.png')" style="padding:0.4rem 0.75rem; font-size:0.75rem;">
            Export Table (PNG)
          </button>
        </div>
      </div>

      <!-- Main Plot Card -->
      <div class="card" style="padding:1.25rem; margin-bottom:1.5rem;">
        <div style="display:flex; justify-content:space-between; align-items:baseline; margin-bottom:1rem; border-bottom:1px solid #f1f5f9; padding-bottom:0.5rem;">
          <div>
            <h3 id="tmb-card-title" style="margin:0; font-size:1rem; font-weight:700; color:var(--text-primary);">Pan-Cancer Tumour Mutation Burden</h3>
            <p id="tmb-card-caption" style="margin:2px 0 0; font-size:0.8rem; color:var(--text-muted);">
              Pan-cancer distribution of somatic mutation burden with alternating cohort columns, sorted S-curves, and median bars
            </p>
          </div>
          <div id="tmb-stat-badge" style="font-size:11px; font-family:var(--font-mono); color:#475569; background:#f1f5f9; padding:3px 8px; border-radius:4px;">
            Loading samples...
          </div>
        </div>

        <!-- TMB Interactive Legend Bar -->
        <div id="tmb-legend-bar" style="display:flex; align-items:center; gap:1.2rem; flex-wrap:wrap; background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px; padding:7px 12px; margin-bottom:10px; font-size:11.5px; font-weight:600; color:#334155;"></div>

        <div id="tmb-plot-container" style="position:relative; width:100%; overflow-x:auto; background:#ffffff; border:1px solid #e2e8f0; border-radius:6px; padding:0.5rem 0;">
          <svg id="tmb-svg" width="1050" height="625" style="display:block; min-width:950px;"></svg>
          <div id="tmb-tooltip" style="position:absolute; display:none; background:#0f172a; color:#ffffff; padding:8px 12px; border-radius:6px; font-size:11px; pointer-events:none; z-index:50; box-shadow:0 10px 15px -3px rgba(0,0,0,0.3); border:1px solid #334155; max-width:280px;"></div>
        </div>
      </div>

      <!-- Cases Above Mutation Threshold Bar Chart Card -->
      <div class="card" style="padding:1.25rem; margin-bottom:1.5rem;">
        <div style="display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:1rem; border-bottom:1px solid #f1f5f9; padding-bottom:0.75rem; flex-wrap:wrap; gap:0.75rem;">
          <div>
            <div style="display:flex; align-items:center; gap:8px;">
              <h3 style="margin:0; font-size:1rem; font-weight:700; color:var(--text-primary);">
                Prevalence of Cases Above Mutation Threshold by Tumour Type
              </h3>
              <span id="tmb-thresh-cohort-badge" style="font-size:11px; font-weight:700; font-family:var(--font-mono); color:#0284c7; background:#e0f2fe; padding:2px 8px; border-radius:4px;">
                Loading...
              </span>
            </div>
            <p style="margin:3px 0 0; font-size:0.8rem; color:var(--text-muted);">
              Percentage and number of cases per tumour type meeting or exceeding the mutation burden cutoff. Ranked by prevalence.
            </p>
          </div>

          <div style="display:flex; align-items:center; gap:8px; flex-wrap:wrap;">
            <div style="display:flex; gap:4px;">
              <button type="button" class="btn btn-secondary" onclick="setTmbThresholdPreset(2500)" style="padding:2px 7px; font-size:11px;">2,500</button>
              <button type="button" class="btn btn-secondary" onclick="setTmbThresholdPreset(4000)" style="padding:2px 7px; font-size:11px; font-weight:700; color:#0284c7;">4,000 (Default)</button>
              <button type="button" class="btn btn-secondary" onclick="setTmbThresholdPreset(10000)" style="padding:2px 7px; font-size:11px;">10,000</button>
              <button type="button" class="btn btn-secondary" onclick="setTmbThresholdPreset(28500)" style="padding:2px 7px; font-size:11px;" title="Equivalent to FDA 10 mut/Mb (WGS)">28,500 (10 mut/Mb)</button>
            </div>
            <button class="btn btn-secondary" onclick="exportTmbThresholdSvg()" style="padding:0.35rem 0.65rem; font-size:0.75rem;">
              Export Plot (SVG)
            </button>
            <button class="btn btn-secondary" onclick="exportTmbThresholdPng()" style="padding:0.35rem 0.65rem; font-size:0.75rem;">
              Export Plot (PNG)
            </button>
          </div>
        </div>

        <div id="tmb-thresh-container" style="position:relative; width:100%; overflow-x:auto; background:#ffffff; border:1px solid #e2e8f0; border-radius:6px; padding:0.75rem 0;">
          <svg id="tmb-thresh-svg" width="980" height="780" style="display:block; min-width:850px;"></svg>
          <div id="tmb-thresh-tooltip" style="position:absolute; display:none; background:#0f172a; color:#ffffff; padding:8px 12px; border-radius:6px; font-size:11px; pointer-events:none; z-index:50; box-shadow:0 10px 15px -3px rgba(0,0,0,0.3); border:1px solid #334155; max-width:280px;"></div>
        </div>
      </div>

      <!-- Cohort Summary Metrics Table Card -->
      <div class="card" style="padding:1.25rem;">
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.75rem; flex-wrap:wrap; gap:0.5rem;">
          <h4 style="margin:0; font-size:0.85rem; font-weight:700; color:var(--text-primary); text-transform:uppercase; letter-spacing:0.04em;">
            Cohort Mutation Burden Summary Table
          </h4>
          <div style="display:flex; align-items:center; gap:0.5rem;">
            <button class="btn btn-outline btn-sm" onclick="exportTableToCSV('tmb-summary-table', 'cohort_tmb_summary.csv')">Export CSV</button>
            <button class="btn btn-outline btn-sm" onclick="exportTableToPNG('tmb-summary-table', 'cohort_tmb_summary.png')">Export PNG</button>
            <span style="font-size:11px; color:var(--text-muted); margin-left:6px;">Click any column header to sort</span>
          </div>
        </div>
        <div style="overflow-x:auto; max-height:420px; overflow-y:auto; border:1px solid #e2e8f0; border-radius:6px;">
          <table id="tmb-summary-table" class="data-table" style="width:100%; font-size:11.5px;">
            <thead>
              <tr>
                <th class="sortable" style="text-align:left;">Cohort / Lineage</th>
                <th class="sortable" style="text-align:left;">Subtype / Stratum</th>
                <th class="num sortable">Biopsies (n)</th>
                <th class="num sortable">Median Value</th>
                <th class="num sortable">IQR (Q1 - Q3)</th>
                <th class="num sortable">Min</th>
                <th class="num sortable">Max</th>
                <th class="num sortable">% TMB-High (≥ 10)</th>
                <th class="num sortable">Median Purity</th>
              </tr>
            </thead>
            <tbody id="tmb-summary-body">
              <tr><td colspan="9" style="text-align:center; padding:1.5rem; color:var(--text-muted);">Loading mutation burden data...</td></tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>

  </main>
</div>

<script>
// Curated ADC Knowledge Base Library
const ADC_CATALOG = [
  { gene: 'ERBB2', alt: 'HER2', name: 'Receptor tyrosine-protein kinase erbB-2', drugs: 'Trastuzumab deruxtecan (Enhertu), Trastuzumab emtansine (Kadcyla)', warhead: 'Topoisomerase I inhibitor (DXd), Microtubule (DM1)', indications: 'Breast, Gastric, NSCLC, Colorectal (Pan-tumor HER2+)', status: 'FDA Approved', tag: 'tag-fda' },
  { gene: 'TACSTD2', alt: 'TROP2', name: 'Tumor-associated calcium signal transducer 2', drugs: 'Sacituzumab govitecan (Trodelvy), Datopotamab deruxtecan (Dato-DXd)', warhead: 'SN-38 (Topo I), DXd (Topo I)', indications: 'TNBC, HR+/HER2- Breast, Urothelial, NSCLC', status: 'FDA Approved', tag: 'tag-fda' },
  { gene: 'NECTIN4', alt: 'PVRL4', name: 'Nectin cell adhesion molecule 4', drugs: 'Enfortumab vedotin (Padcev)', warhead: 'MMAE (Microtubule inhibitor)', indications: 'Urothelial / Bladder carcinoma (+ Pembrolizumab)', status: 'FDA Approved', tag: 'tag-fda' },
  { gene: 'FOLR1', alt: 'FRa', name: 'Folate receptor alpha', drugs: 'Mirvetuximab soravtansine (Elahere)', warhead: 'DM4 (Maytansinoid tubulin inhibitor)', indications: 'Platinum-resistant ovarian, fallopian tube, peritoneal cancer', status: 'FDA Approved', tag: 'tag-fda' },
  { gene: 'MET', alt: 'c-MET', name: 'Hepatocyte growth factor receptor', drugs: 'Telisotuzumab vedotin (Teliso-V), ABBV-399', warhead: 'MMAE (Microtubule inhibitor)', indications: 'EGFR-wildtype Non-Small Cell Lung Cancer (NSCLC)', status: 'Phase 3 / BTD', tag: 'tag-p3' },
  { gene: 'ERBB3', alt: 'HER3', name: 'Receptor tyrosine-protein kinase erbB-3', drugs: 'Patritumab deruxtecan (HER3-DXd)', warhead: 'Topoisomerase I inhibitor (DXd)', indications: 'EGFR-mutated NSCLC, HR+/HER2- Breast cancer', status: 'Phase 3 / Priority', tag: 'tag-p3' },
  { gene: 'CD276', alt: 'B7-H3', name: 'CD276 antigen (B7 homolog 3)', drugs: 'Ifinatamab deruxtecan (I-DXd), HS-20093', warhead: 'Topoisomerase I inhibitor (DXd)', indications: 'Small Cell Lung Cancer (SCLC), Prostate, Esophageal', status: 'Phase 2 / Phase 3', tag: 'tag-p3' },
  { gene: 'CLDN18', alt: 'Claudin-18.2', name: 'Claudin-18 (Isoform 2)', drugs: 'Zolbetuximab (Vyloy - mAb), ABBV-400 / Claudin ADCs', warhead: 'Topoisomerase I inhibitor / MMAE', indications: 'Gastric & Gastroesophageal junction adenocarcinoma', status: 'FDA Approved (mAb) / P3 (ADC)', tag: 'tag-fda' },
  { gene: 'CEACAM5', alt: 'CEA', name: 'Carcinoembryonic antigen-related cell adhesion molecule 5', drugs: 'Tusamitamab ravtansine (SAR408701)', warhead: 'DM4 (Maytansinoid tubulin inhibitor)', indications: 'Non-squamous NSCLC, Colorectal carcinoma', status: 'Phase 3', tag: 'tag-p3' },
  { gene: 'EGFR', alt: 'HER1', name: 'Epidermal growth factor receptor', drugs: 'MRG003, Depatuxizumab mafodotin (ABT-414)', warhead: 'MMAE / MMAF (Microtubule inhibitors)', indications: 'Nasopharyngeal (NPC), HNSCC, Glioblastoma', status: 'Phase 2 / Phase 3', tag: 'tag-p2' },
  { gene: 'FOLH1', alt: 'PSMA', name: 'Glutamate carboxypeptidase 2 / PSMA', drugs: 'ARX517, PSMA-ADC', warhead: 'AS269 (Tubulin inhibitor) / MMAE', indications: 'Metastatic Castration-Resistant Prostate Cancer (mCRPC)', status: 'Phase 2', tag: 'tag-p2' },
  { gene: 'DLL3', alt: 'Delta-like 3', name: 'Delta-like canonical Notch ligand 3', drugs: 'Rovalpituzumab tesirine (Rova-T), HPN328', warhead: 'PBD dimer (DNA alkylator) / Topo I', indications: 'Small Cell Lung Cancer (SCLC), Neuroendocrine prostate', status: 'Phase 2', tag: 'tag-p2' },
  { gene: 'VTCN1', alt: 'B7-H4', name: 'V-set domain-containing T-cell activation inhibitor 1 (B7-H4)', drugs: 'AZD8205 (B7-H4 Topo I ADC), XMT-1660 (DolaLock ADC), HS-20089', warhead: 'Topoisomerase I inhibitor / Auristatin (DolaLock)', indications: 'Breast (TNBC, HR+), High-grade serous ovarian, Endometrial carcinoma', status: 'Phase 1 / Phase 2', tag: 'tag-p2' }
];

const ADC_PRESETS = {
  horizon: { label: 'Horizon-Beacon', genes: ['ERBB2', 'TACSTD2', 'CLDN18', 'FOLR1', 'VTCN1'] },
  fda: { label: 'FDA-Approved ADCs', genes: ['ERBB2', 'TACSTD2', 'NECTIN4', 'FOLR1'] },
  pipeline: { label: 'Phase 3 Pipeline', genes: ['MET', 'ERBB3', 'CD276', 'CLDN18'] },
  thoracic: { label: 'Thoracic / NSCLC', genes: ['EGFR', 'MET', 'ERBB2', 'ERBB3', 'CEACAM5'] },
  gi: { label: 'GI & Pancreatic', genes: ['CLDN18', 'TACSTD2', 'CEACAM5', 'ERBB2'] },
  breast: { label: 'Breast & GYN', genes: ['ERBB2', 'TACSTD2', 'FOLR1', 'ERBB3'] },
  gu: { label: 'GU & Prostate', genes: ['NECTIN4', 'FOLH1', 'TACSTD2', 'CD276'] },
  neuro: { label: 'Neuroendocrine / SCLC', genes: ['DLL3', 'CD276'] },
  all: { label: 'All Curated ADCs', genes: ['ERBB2', 'TACSTD2', 'CLDN18', 'FOLR1', 'VTCN1', 'NECTIN4', 'MET', 'ERBB3', 'CD276', 'CEACAM5', 'EGFR', 'FOLH1', 'DLL3'] }
};

let currentActiveGene = 'EGFR';
let currentActiveView = 'single-gene';
let currentMultiSubTab = 'summary';
let currentMatrixMode = 'heatmap';
let lastMultiGeneData = null;
let allAvailableCancers = [];
let selectedTumorTypes = [];

// Initialize Studio
async function init() {
  try {
    const res = await fetch('/api/cohort_overview');
    const data = await res.json();
    document.getElementById('cohort-badge').textContent = `${data.total_samples} samples | ${data.total_genes.toLocaleString()} genes`;
    document.getElementById('stat-samples').textContent = data.total_samples;
    document.getElementById('stat-genes').textContent = data.total_genes.toLocaleString();
    if (data.cancers_count) document.getElementById('stat-cancers').textContent = data.cancers_count;

    if (data.cancer_types) {
      allAvailableCancers = data.cancer_types;
      populateTumorCheckboxes(allAvailableCancers);
    }
  } catch (e) {
    console.error("Overview error:", e);
  }

  // Populate ADC catalog
  renderAdcCatalog();

  // Load default single-gene data
  loadGeneData();

  // Universal table sorting initialization
  initAllSortableTables();

  // Close dropdown on outside click
  document.addEventListener('click', function(e) {
    const menu = document.getElementById('tumor-dropdown-menu');
    const btn = document.getElementById('tumor-filter-btn');
    if (menu && menu.classList.contains('open')) {
      if (!menu.contains(e.target) && !btn.contains(e.target)) {
        menu.classList.remove('open');
      }
    }
  });
}

// ---------------------------------------------------------------------------
// TUMOR TYPE FILTERING LOGIC
// ---------------------------------------------------------------------------
function populateTumorCheckboxes(cancers) {
  const container = document.getElementById('tumor-checkbox-list');
  if (!container) return;
  let html = '';
  cancers.forEach(c => {
    html += `
      <label class="tumor-item-row" data-cancer="${c.cancer_type.toLowerCase()}">
        <div style="display:flex; align-items:center; gap:6px;">
          <input type="checkbox" value="${c.cancer_type}" class="tumor-cb" checked onchange="updateTumorSelectionCount()">
          <span style="font-weight:500;">${c.cancer_type}</span>
        </div>
        <span style="font-family:var(--font-mono); font-size:10.5px; color:var(--text-muted);">N=${c.count}</span>
      </label>
    `;
  });
  container.innerHTML = html;
  updateTumorSelectionCount();
}

function toggleTumorDropdown(e) {
  if (e) e.stopPropagation();
  const menu = document.getElementById('tumor-dropdown-menu');
  if (menu) menu.classList.toggle('open');
}

function filterTumorCheckboxList() {
  const query = document.getElementById('tumor-search-input').value.toLowerCase().trim();
  document.querySelectorAll('.tumor-item-row').forEach(row => {
    const text = row.getAttribute('data-cancer') || '';
    row.style.display = text.includes(query) ? 'flex' : 'none';
  });
}

function selectAllTumors(check) {
  document.querySelectorAll('.tumor-cb').forEach(cb => {
    const row = cb.closest('.tumor-item-row');
    if (row && row.style.display !== 'none') {
      cb.checked = check;
    }
  });
  updateTumorSelectionCount();
}

function setQuickTumorCohort(cancers) {
  document.querySelectorAll('.tumor-cb').forEach(cb => {
    cb.checked = cancers.includes(cb.value);
  });
  updateTumorSelectionCount();
  applyTumorSelection();
}

function updateTumorSelectionCount() {
  const cbs = document.querySelectorAll('.tumor-cb');
  const checked = Array.from(cbs).filter(cb => cb.checked);
  const total = cbs.length;
  const countLabel = document.getElementById('tumor-selected-count-label');
  if (countLabel) {
    countLabel.textContent = `${checked.length} of ${total} cohorts selected`;
  }
}

function applyTumorSelection() {
  const cbs = document.querySelectorAll('.tumor-cb');
  const checked = Array.from(cbs).filter(cb => cb.checked).map(cb => cb.value);
  const total = cbs.length;

  selectedTumorTypes = (checked.length === total || checked.length === 0) ? [] : checked;

  const btnLabel = document.getElementById('tumor-filter-label');
  if (checked.length === total || checked.length === 0) {
    btnLabel.textContent = `All Tumor Types (${total})`;
  } else if (checked.length === 1) {
    btnLabel.textContent = `${checked[0]}`;
  } else if (checked.length <= 3) {
    btnLabel.textContent = `${checked.join(', ')}`;
  } else {
    btnLabel.textContent = `${checked.length} Cohorts Selected`;
  }

  const menu = document.getElementById('tumor-dropdown-menu');
  if (menu) menu.classList.remove('open');

  runMultiGeneAnalysis();
}

// Sidebar Navigation
function switchView(viewId) {
  if (viewId === 'umap' && (!window._umapRunDone)) {
    window._umapRunDone = true;
    loadCohortUmap();
  }
  if (viewId === 'tmb' && (!window._tmbRunDone)) {
    window._tmbRunDone = true;
    loadTmbData();
  }
  currentActiveView = viewId;
  document.querySelectorAll('.view-panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-btn').forEach(b => b.classList.remove('active'));

  const targetPanel = document.getElementById(`view-${viewId}`);
  const targetNav = document.getElementById(`nav-${viewId}`);
  if (targetPanel) targetPanel.classList.add('active');
  if (targetNav) targetNav.classList.add('active');

  if (viewId === 'multi-gene' && (!window._multiGeneRunDone)) {
    window._multiGeneRunDone = true;
    runMultiGeneAnalysis();
  }
  if (viewId === 'diff-expr' && (!window._diffRunDone)) {
    window._diffRunDone = true;
    runDifferentialExpression();
  }
}

// Sub-Tab Navigation for Multi-Gene Studio
function switchMultiSubTab(tabKey) {
  currentMultiSubTab = tabKey;
  document.querySelectorAll('.sub-tab').forEach(b => b.classList.remove('active'));
  document.querySelectorAll('.sub-tab-panel').forEach(p => p.classList.remove('active'));

  const targetBtn = document.getElementById(`tab-btn-${tabKey}`);
  const targetPanel = document.getElementById(`multi-subtab-${tabKey}`);
  if (targetBtn) targetBtn.classList.add('active');
  if (targetPanel) targetPanel.classList.add('active');

  if (tabKey === 'matrix' && lastMultiGeneData) {
    if (currentMatrixMode === 'heatmap') {
      renderMultiGeneHeatmap(lastMultiGeneData);
    }
  }
}

// Matrix View Mode: Heatmap vs Table
function setMatrixView(mode) {
  currentMatrixMode = mode;
  const btnHeatmap = document.getElementById('btn-seg-heatmap');
  const btnTable = document.getElementById('btn-seg-table');
  const viewHeatmap = document.getElementById('matrix-heatmap-view');
  const viewTable = document.getElementById('matrix-table-view');
  const btnSvg = document.getElementById('btn-export-matrix-svg');
  const btnPng = document.getElementById('btn-export-matrix-png');

  const legendBar = document.getElementById('heatmap-legend-bar');
  if (mode === 'heatmap') {
    btnHeatmap.classList.add('active');
    btnTable.classList.remove('active');
    viewHeatmap.style.display = 'block';
    viewTable.style.display = 'none';
    if (btnSvg) btnSvg.style.display = 'inline-flex';
    if (btnPng) btnPng.style.display = 'inline-flex';
    if (legendBar) legendBar.style.display = 'flex';
    if (lastMultiGeneData) renderMultiGeneHeatmap(lastMultiGeneData);
  } else {
    btnHeatmap.classList.remove('active');
    btnTable.classList.add('active');
    viewHeatmap.style.display = 'none';
    viewTable.style.display = 'block';
    if (btnSvg) btnSvg.style.display = 'none';
    if (btnPng) btnPng.style.display = 'inline-flex';
    if (legendBar) legendBar.style.display = 'none';
  }
}

// Quick gene analyze removed; each tab maintains its own focused search input

// ---------------------------------------------------------------------------
// VIEW 1: SINGLE GENE PROFILE
// ---------------------------------------------------------------------------
async function loadGeneData() {
  const gene = document.getElementById('gene-search').value.trim() || currentActiveGene || 'EGFR';
  currentActiveGene = gene.toUpperCase();
  const metric = document.getElementById('metric-select').value;
  const cutoff = document.getElementById('cutoff-select').value;

  document.getElementById('title-gene').textContent = currentActiveGene;
  // global search box removed

  // 1. Distribution
  try {
    const res = await fetch(`/api/genes/${currentActiveGene}/distribution?metric=${metric}`);
    const data = await res.json();
    renderDistribution(data);
  } catch(e) { console.error(e); }

  // 2. Over-expression
  try {
    const res = await fetch(`/api/genes/${currentActiveGene}/overexpression?threshold=${cutoff}`);
    const data = await res.json();
    renderOverexpression(data);
  } catch(e) { console.error(e); }
}

function renderDistribution(data) {
  const container = document.getElementById('dist-container');
  if (!data.groups || data.groups.length === 0) {
    container.innerHTML = '<div style="padding:2rem; text-align:center; color:var(--text-muted);">No expression points found for gene.</div>';
    return;
  }

  let minVal = 0, maxVal = 10;
  data.groups.forEach(g => {
    g.points.forEach(p => {
      if (p.expr > maxVal) maxVal = p.expr;
    });
  });
  maxVal = Math.ceil(maxVal * 1.1);

  let html = '';
  data.groups.forEach(g => {
    const leftPct = (g.q25 / maxVal) * 100;
    const widthPct = Math.max(2, ((g.q75 - g.q25) / maxVal) * 100);
    const medPct = (g.median / maxVal) * 100;

    html += `
      <div class="box-row">
        <div class="box-label" title="${g.group}">${g.group} (N=${g.count})</div>
        <div class="box-canvas">
          <div class="box-iqr" style="left:${leftPct}%; width:${widthPct}%;"></div>
          <div class="box-med" style="left:${medPct}%;"></div>
    `;

    g.points.forEach(p => {
      const pPct = Math.min(99, Math.max(0, (p.expr / maxVal) * 100));
      html += `<div class="box-point" style="left:${pPct}%;" title="${p.sample_id}: ${p.expr} ${data.metric} (Purity: ${p.purity}, CUPPA: ${p.cuppa})"></div>`;
    });

    html += `
        </div>
        <div style="font-family:var(--font-mono); font-size:11px; width:45px; text-align:right;">${g.median.toFixed(1)}</div>
      </div>
    `;
  });

  container.innerHTML = html;
}

function renderOverexpression(data) {
  const tbody = document.getElementById('threshold-body');
  if (!data.cancer_types || data.cancer_types.length === 0) {
    tbody.innerHTML = '<tr><td colspan="5" style="text-align:center; color:var(--text-muted); padding:1rem;">No data</td></tr>';
    return;
  }

  let rows = '';
  data.cancer_types.forEach(c => {
    const barWidth = Math.min(100, c.percent_elevated);
    rows += `
      <tr>
        <td style="font-weight:500;">${c.cancer_type}</td>
        <td class="num" style="font-weight:600; color:${c.count_elevated > 0 ? 'var(--down-reg)' : 'inherit'};">${c.count_elevated}</td>
        <td class="num">${c.total_patients}</td>
        <td class="num" style="font-weight:600;">${c.percent_elevated.toFixed(1)}%</td>
        <td>
          <div style="width:100%; height:7px; background:#f1f5f9; border-radius:2px; overflow:hidden;">
            <div style="width:${barWidth}%; height:100%; background:${c.percent_elevated >= 20 ? 'var(--down-reg)' : 'var(--accent)'};"></div>
          </div>
        </td>
      </tr>
    `;
  });
  tbody.innerHTML = rows;
  makeTableSortable('threshold-table');
}

// ---------------------------------------------------------------------------
// VIEW 2: MULTI-GENE & ADC STUDIO
// ---------------------------------------------------------------------------
function loadAdcPreset(presetKey) {
  document.querySelectorAll('.preset-strip .chip').forEach(c => c.classList.remove('active'));
  if (event && event.target) event.target.classList.add('active');

  const p = ADC_PRESETS[presetKey];
  if (!p) return;
  document.getElementById('multi-gene-input').value = p.genes.join(', ');
  runMultiGeneAnalysis();
}

async function runMultiGeneAnalysis() {
  const rawInput = document.getElementById('multi-gene-input').value.trim();
  const threshold = document.getElementById('multi-threshold-select').value;
  const weightMode = document.getElementById('multi-weight-select').value;
  const minPurity = document.getElementById('multi-purity-filter').value;
  const cutoffBadge = document.getElementById('multi-chart-cutoff-badge');
  if (cutoffBadge) cutoffBadge.textContent = `Cutoff: ≥ ${threshold}th Percentile`;

  if (!rawInput) return;

  const barContainer = document.getElementById('multi-bar-container');
  barContainer.innerHTML = '<div style="color:var(--text-muted); padding:2rem; text-align:center;">Evaluating multi-target cohort...</div>';

  try {
    let url = `/api/genes/multi_overexpression?genes=${encodeURIComponent(rawInput)}&threshold=${threshold}&weight_mode=${weightMode}&min_purity=${minPurity}&purity_floor=0.15`;
    if (selectedTumorTypes && selectedTumorTypes.length > 0) {
      url += `&cancers=${encodeURIComponent(selectedTumorTypes.join(','))}`;
    }

    const res = await fetch(url);
    const data = await res.json();
    lastMultiGeneData = data;
    renderMultiGeneResults(data);
  } catch (e) {
    console.error("Multi-gene analysis failed:", e);
    barContainer.innerHTML = `<div style="color:#ef4444; padding:2rem; text-align:center;">Error evaluating multi-gene cohort: ${e.message}</div>`;
  }
}

function renderMultiGeneResults(data) {
  const isSubset = data.selected_cancers && data.selected_cancers.length > 0;
  const isPurityFiltered = data.min_purity && data.min_purity > 0;

  // Update Purity Floor Annotation Banner
  const banner = document.getElementById('purity-annotation-banner');
  const bannerText = document.getElementById('purity-banner-text');

  if (banner && bannerText) {
    if (data.weight_mode === 'tumor_adjusted') {
      banner.style.display = 'block';
      banner.style.background = '#f8fafc';
      banner.style.borderColor = '#cbd5e1';
      banner.style.borderLeftColor = 'var(--accent)';
      let msg = `<strong>Tumor-Cell Adjusted Mode Active (TPM / Purity):</strong> Expression is normalized to estimated neoplastic cellularity. <strong>Purity Floor: 15% (0.15)</strong> is applied. `;
      if (data.below_floor_count > 0) {
        msg += `<strong>${data.below_floor_count} low-purity biopsies (${data.below_floor_pct}%)</strong> are clamped to the 15% floor to eliminate division-by-zero noise amplification. `;
      } else {
        msg += `All evaluated biopsies meet or exceed the 15% cellularity floor. `;
      }
      if (data.uncalled_purity_count > 0) {
        msg += `<strong>${data.uncalled_purity_count} uncalled biopsies (${data.uncalled_purity_pct}%)</strong> are evaluated at the cohort median cellularity (40%).`;
      }
      bannerText.innerHTML = msg;
    } else if (data.weight_mode === 'tme_weighted') {
      banner.style.display = 'block';
      banner.style.background = '#f0fdf4';
      banner.style.borderColor = '#bbf7d0';
      banner.style.borderLeftColor = '#16a34a';
      bannerText.innerHTML = `<strong>TME-Weighted Mode Active (TPM &times; (1 - Purity)):</strong> Quantifies transcript abundance attributable to the non-malignant tumor microenvironment and stromal infiltrate. <strong>Purity Floor: 15%</strong> is enforced.`;
    } else if (isPurityFiltered) {
      banner.style.display = 'block';
      banner.style.background = '#f0fdf4';
      banner.style.borderColor = '#bbf7d0';
      banner.style.borderLeftColor = '#16a34a';
      bannerText.innerHTML = `<strong>Purity Filter Active (&ge; ${Math.round(data.min_purity * 100)}%):</strong> Low-cellularity stromal biopsies excluded. All <strong>${data.total_cohort.toLocaleString()}</strong> evaluated biopsies strictly exceed the 15% purity floor.`;
    } else {
      banner.style.display = 'none';
    }
  }

  // Update KPI Cards
  document.getElementById('kpi-total-cohort').textContent = (data.total_cohort || 0).toLocaleString();
  const cohortSub = document.getElementById('kpi-cohort-sub');
  if (cohortSub) {
    let sub = `${data.cancer_matrix.length} tumor cohort${data.cancer_matrix.length === 1 ? '' : 's'}`;
    if (isPurityFiltered) sub += ` (&ge; ${Math.round(data.min_purity * 100)}% purity)`;
    cohortSub.textContent = sub;
  }

  document.getElementById('kpi-any-elevated').textContent = `${data.any_target_count || 0} (${data.any_target_pct || 0}%)`;
  document.getElementById('kpi-co-elevated').textContent = `${data.co_elevated_count || 0} (${data.co_elevated_pct || 0}%)`;
  
  let targetSub = `${data.genes.length} target basket candidate cohort`;
  if (data.weight_mode === 'tumor_adjusted') targetSub += ' (Tumor-Cell Adjusted)';
  else if (data.weight_mode === 'tme_weighted') targetSub += ' (TME-Weighted)';
  document.getElementById('kpi-any-sub').textContent = targetSub;

  // Update Chart Subtitle
  const chartSub = document.getElementById('multi-chart-sub');
  if (chartSub) {
    let modeLabel = data.weight_mode === 'tumor_adjusted' ? 'Tumor-Cell Intrinsic (TPM / Purity, Floor: 15%)' : (data.weight_mode === 'tme_weighted' ? 'TME-Weighted' : 'Bulk Expression');
    chartSub.textContent = `Prevalence across ${data.total_cohort.toLocaleString()} patients in ${data.cancer_matrix.length} tumor cohort${data.cancer_matrix.length === 1 ? '' : 's'} (${modeLabel})`;
  }

  // Render SVG Horizontal Bar Chart
  renderMultiGeneBarChart(data.gene_summary, data.total_cohort, data.weight_mode);

  // Render Gene Summary Table
  const gBody = document.getElementById('multi-gene-summary-body');
  if (!data.gene_summary || data.gene_summary.length === 0) {
    gBody.innerHTML = '<tr><td colspan="6" style="text-align:center; color:var(--text-muted); padding:2rem;">No genes evaluated.</td></tr>';
  } else {
    let rows = '';
    data.gene_summary.forEach(g => {
      const hasClamped = g.below_floor_elevated && g.below_floor_elevated > 0;
      rows += `
        <tr>
          <td style="font-weight:600; font-family:var(--font-mono); color:var(--accent);">${g.gene}</td>
          <td class="num" style="font-weight:600;">${g.elevated_count}</td>
          <td class="num">${g.total_cohort}</td>
          <td class="num" style="font-weight:700; color:${g.percent_elevated >= 20 ? 'var(--down-reg)' : 'var(--text-primary)'};">${g.percent_elevated.toFixed(1)}%</td>
          <td class="num" style="font-family:var(--font-mono); font-size:11px;">
            ${hasClamped ? `<span class="tag tag-clamped" title="${g.below_floor_elevated} elevated cases had PURPLE tumor purity < 15% and were clamped to the 15% floor">${g.below_floor_elevated} clamped</span>` : '<span style="color:var(--text-muted);">-</span>'}
          </td>
          <td style="text-align:center;">
            <button class="btn btn-outline btn-sm" onclick="profileSpecificGene('${g.gene}')">Profile</button>
          </td>
        </tr>
      `;
    });
    gBody.innerHTML = rows;
  }
  makeTableSortable('multi-gene-summary-table');

  // Render Cancer Matrix Table
  const mHead = document.getElementById('multi-cancer-head');
  let headHtml = `
    <tr>
      <th class="sortable">Cancer Type</th>
      <th class="num sortable">Patients</th>
      <th class="num sortable">Any Target N</th>
      <th class="num sortable">Any Target %</th>
      <th class="num sortable" title="Biopsies with PURPLE tumor purity < 15%">Below Floor (&lt;15%)</th>
  `;
  data.genes.forEach(g => {
    headHtml += `<th class="num sortable">${g} %</th>`;
  });
  headHtml += `</tr>`;
  mHead.innerHTML = headHtml;

  const mBody = document.getElementById('multi-cancer-body');
  if (!data.cancer_matrix || data.cancer_matrix.length === 0) {
    mBody.innerHTML = `<tr><td colspan="${5 + data.genes.length}" style="text-align:center; color:var(--text-muted); padding:2rem;">No cancer types found matching selection.</td></tr>`;
  } else {
    let rows = '';
    data.cancer_matrix.forEach(r => {
      const hasClamped = r.below_floor_count > 0;
      rows += `
        <tr>
          <td style="font-weight:500;">${r.cancer_type}</td>
          <td class="num">${r.total_patients}</td>
          <td class="num" style="font-weight:600; color:var(--accent);">${r.any_elevated_count}</td>
          <td class="num" style="font-weight:700;">${r.any_elevated_pct.toFixed(1)}%</td>
          <td class="num" style="font-family:var(--font-mono); font-size:11px;">
            ${hasClamped ? `<span class="tag tag-clamped" title="${r.below_floor_count} biopsies have PURPLE tumor purity < 15% and were clamped to floor">${r.below_floor_count}</span>` : '<span style="color:var(--text-muted);">0</span>'}
          </td>
      `;
      data.genes.forEach(g => {
        const val = r.gene_pcts[g] || 0.0;
        const color = val >= 25 ? '#b91c1c' : (val >= 10 ? '#dc2626' : (val > 0 ? '#1e293b' : 'var(--text-muted)'));
        rows += `<td class="num" style="color:${color}; font-weight:${val > 0 ? '600' : 'normal'};">${val.toFixed(1)}%</td>`;
      });
      rows += `</tr>`;
    });
    mBody.innerHTML = rows;
  }
  makeTableSortable('multi-cancer-matrix-table');

  // Render Interactive Heatmap
  renderMultiGeneHeatmap(data);
}

function renderMultiGeneBarChart(geneSummary, totalCohort, weightMode) {
  const container = document.getElementById('multi-bar-container');
  if (!geneSummary || geneSummary.length === 0) {
    container.innerHTML = '<div style="color:var(--text-muted); padding:2rem; text-align:center;">No genes to display.</div>';
    return;
  }

  const sorted = [...geneSummary].sort((a, b) => b.elevated_count - a.elevated_count);

  const barH = 26;
  const gap = 12;
  const padLeft = 85;
  const padRight = 130;
  const padTop = 15;
  const padBottom = 25;
  const w = 780;
  const plotW = w - padLeft - padRight;
  const h = padTop + sorted.length * (barH + gap) + padBottom;

  const maxPct = Math.max(25, ...sorted.map(s => s.percent_elevated));
  const maxAxisVal = Math.ceil(maxPct / 10) * 10;

  let out = `<svg viewBox="0 0 ${w} ${h}" style="width:100%; height:${h}px; font-family:var(--font-sans);">`;

  // Background grid lines
  for (let p = 0; p <= maxAxisVal; p += (maxAxisVal <= 30 ? 5 : 10)) {
    const x = padLeft + (p / maxAxisVal) * plotW;
    out += `<line x1="${x}" y1="${padTop}" x2="${x}" y2="${h - padBottom}" stroke="#f1f5f9" stroke-width="1" />`;
    out += `<text x="${x}" y="${h - 8}" text-anchor="middle" font-size="10" font-family="var(--font-mono)" fill="#94a3b8">${p}%</text>`;
  }

  // Baseline X axis
  out += `<line x1="${padLeft}" y1="${h - padBottom}" x2="${padLeft + plotW}" y2="${h - padBottom}" stroke="#cbd5e1" stroke-width="1.5" />`;

  // Bars
  sorted.forEach((item, idx) => {
    const y = padTop + idx * (barH + gap);
    const barW = Math.max(2, (item.percent_elevated / maxAxisVal) * plotW);
    const isHigh = item.percent_elevated >= 25;
    const barColor = isHigh ? '#e11d48' : '#2563eb';

    // Gene Label
    out += `<text x="${padLeft - 10}" y="${y + barH / 2 + 4}" text-anchor="end" font-size="12" font-family="var(--font-mono)" font-weight="600" fill="#1e293b" style="cursor:pointer;" onclick="profileSpecificGene('${item.gene}')">${item.gene}</text>`;

    // Background track
    out += `<rect x="${padLeft}" y="${y}" width="${plotW}" height="${barH}" rx="3" fill="#f8fafc" />`;

    // Active Bar with hover tooltip
    out += `
      <rect x="${padLeft}" y="${y}" width="${barW}" height="${barH}" rx="3" fill="${barColor}" opacity="0.9" style="cursor:pointer;" onclick="profileSpecificGene('${item.gene}')">
        <title>${item.gene}: ${item.elevated_count} / ${item.total_cohort} patients (${item.percent_elevated.toFixed(1)}%)${item.below_floor_elevated ? ' [' + item.below_floor_elevated + ' clamped to 15% floor]' : ''}</title>
      </rect>
    `;

    // Value direct label
    out += `<text x="${padLeft + barW + 8}" y="${y + barH / 2 + 4}" font-size="11.5" font-family="var(--font-mono)" font-weight="600" fill="#334155">${item.elevated_count} patients (${item.percent_elevated.toFixed(1)}%)</text>`;
  });

  out += `</svg>`;
  container.innerHTML = out;
}

// Interactive Cross-Cancer Target Heatmap
let heatmapSortState = {
  key: 'any_target',
  asc: false
};

function sortHeatmapBy(key) {
  if (heatmapSortState.key === key) {
    heatmapSortState.asc = !heatmapSortState.asc;
  } else {
    heatmapSortState.key = key;
    heatmapSortState.asc = (key === 'cancer_type') ? true : false;
  }
  if (lastMultiGeneData) {
    renderMultiGeneHeatmap(lastMultiGeneData);
  }
}

// Interactive Cross-Cancer Target Heatmap with Column-Click Sorting
function renderMultiGeneHeatmap(data) {
  const container = document.getElementById('heatmap-container');
  if (!data || !data.cancer_matrix || data.cancer_matrix.length === 0) {
    container.innerHTML = '<div style="color:var(--text-muted); padding:3rem; text-align:center;">No heatmap data available for current selection.</div>';
    return;
  }

  const rawRows = data.cancer_matrix;
  const genes = data.genes;

  // Re-order (sort) rows based on clicked column
  const rows = [...rawRows].sort((a, b) => {
    if (heatmapSortState.key === 'cancer_type') {
      const cmp = a.cancer_type.localeCompare(b.cancer_type);
      return heatmapSortState.asc ? cmp : -cmp;
    }
    let pA = 0, pB = 0;
    if (heatmapSortState.key === 'any_target') {
      pA = a.any_elevated_pct || 0;
      pB = b.any_elevated_pct || 0;
    } else {
      const g = heatmapSortState.key;
      pA = (a.gene_pcts && a.gene_pcts[g] !== undefined) ? a.gene_pcts[g] : 0;
      pB = (b.gene_pcts && b.gene_pcts[g] !== undefined) ? b.gene_pcts[g] : 0;
    }
    if (pA !== pB) {
      return heatmapSortState.asc ? pA - pB : pB - pA;
    }
    return a.cancer_type.localeCompare(b.cancer_type);
  });

  const rowH = 26;
  const colW = 86;
  const labelW = 165;
  const headerH = 48;
  const anyColW = 96;

  const totalW = labelW + anyColW + genes.length * colW + 20;
  const totalH = headerH + rows.length * rowH + 20;

  // Pure monotonic White-to-Red continuous sequential gradient (0% White -> 100% Dark Crimson)
  function getCellBg(pct) {
    if (pct <= 0) return '#ffffff';
    const p = Math.min(100, Math.max(0, pct));
    const stops = [
      { p: 0,   r: 255, g: 255, b: 255 }, // White
      { p: 10,  r: 255, g: 245, b: 240 }, // Pale warm pink
      { p: 25,  r: 252, g: 187, b: 161 }, // Soft coral pink
      { p: 50,  r: 251, g: 106, b: 74  }, // Medium vivid coral
      { p: 75,  r: 203, g: 24,  b: 29  }, // Strong red
      { p: 100, r: 103, g: 0,   b: 13  }  // Deep rich crimson
    ];
    let lower = stops[0], upper = stops[stops.length - 1];
    for (let i = 0; i < stops.length - 1; i++) {
      if (p >= stops[i].p && p <= stops[i+1].p) {
        lower = stops[i];
        upper = stops[i+1];
        break;
      }
    }
    const span = upper.p - lower.p || 1;
    const frac = (p - lower.p) / span;
    const r = Math.round(lower.r + frac * (upper.r - lower.r));
    const g = Math.round(lower.g + frac * (upper.g - lower.g));
    const b = Math.round(lower.b + frac * (upper.b - lower.b));
    return `rgb(${r},${g},${b})`;
  }

  function getCellTextColor(pct) {
    // White text on dark red backgrounds, dark text on white/light pink backgrounds
    return pct >= 42 ? '#ffffff' : '#0f172a';
  }

  let svg = `<svg id="target-matrix-svg" viewBox="0 0 ${totalW} ${totalH}" style="width:${totalW}px; min-width:100%; height:${totalH}px; font-family:var(--font-sans); display:block;">`;

  // 1. Column Headers with Interactive Sorting on Click
  // Corner Cancer Type
  const isCancerSorted = heatmapSortState.key === 'cancer_type';
  const cancerArrow = isCancerSorted ? (heatmapSortState.asc ? ' ▲' : ' ▼') : ' ↕';
  svg += `
    <rect x="2" y="8" width="${labelW - 6}" height="32" rx="4" fill="${isCancerSorted ? '#eff6ff' : '#f8fafc'}" stroke="${isCancerSorted ? '#3b82f6' : '#e2e8f0'}" stroke-width="${isCancerSorted ? '1.5' : '1'}" style="cursor:pointer;" onclick="sortHeatmapBy('cancer_type')">
      <title>Click to sort by cancer type name</title>
    </rect>
    <text x="12" y="28" font-size="11" font-weight="700" fill="${isCancerSorted ? '#1d4ed8' : '#475569'}" style="cursor:pointer;" onclick="sortHeatmapBy('cancer_type')">CANCER TYPE <tspan font-size="10" fill="${isCancerSorted ? '#2563eb' : '#94a3b8'}">${cancerArrow}</tspan></text>
  `;

  // "Any Target" Column Header
  const isAnySorted = heatmapSortState.key === 'any_target';
  const anyArrow = isAnySorted ? (heatmapSortState.asc ? ' ▲' : ' ▼') : ' ↕';
  svg += `
    <rect x="${labelW}" y="8" width="${anyColW - 4}" height="32" rx="4" fill="${isAnySorted ? '#dbeafe' : '#eff6ff'}" stroke="${isAnySorted ? '#2563eb' : '#bfdbfe'}" stroke-width="${isAnySorted ? '1.5' : '1'}" style="cursor:pointer;" onclick="sortHeatmapBy('any_target')">
      <title>Click to sort cohorts by Any Target prevalence</title>
    </rect>
    <text x="${labelW + (anyColW - 4) / 2}" y="28" text-anchor="middle" font-size="11" font-weight="700" fill="#1d4ed8" style="cursor:pointer;" onclick="sortHeatmapBy('any_target')">ANY TARGET <tspan font-size="10" fill="#2563eb">${anyArrow}</tspan></text>
  `;

  // Gene Column Headers
  genes.forEach((g, gIdx) => {
    const gx = labelW + anyColW + gIdx * colW;
    const isGeneSorted = heatmapSortState.key === g;
    const geneArrow = isGeneSorted ? (heatmapSortState.asc ? ' ▲' : ' ▼') : ' ↕';
    svg += `
      <rect x="${gx}" y="8" width="${colW - 4}" height="32" rx="4" fill="${isGeneSorted ? '#dbeafe' : '#f8fafc'}" stroke="${isGeneSorted ? '#2563eb' : '#e2e8f0'}" stroke-width="${isGeneSorted ? '1.5' : '1'}" style="cursor:pointer;" onclick="sortHeatmapBy('${g}')">
        <title>Click to sort cancer cohorts by ${g} expression prevalence</title>
      </rect>
      <text x="${gx + (colW - 4) / 2}" y="28" text-anchor="middle" font-size="11" font-weight="700" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" fill="${isGeneSorted ? '#1d4ed8' : '#0f172a'}" style="cursor:pointer;" onclick="sortHeatmapBy('${g}')">${g} <tspan font-size="9" fill="${isGeneSorted ? '#2563eb' : '#94a3b8'}">${geneArrow}</tspan></text>
    `;
  });

  // 2. Matrix Rows (Rendered in sorted order)
  rows.forEach((r, rIdx) => {
    const ry = headerH + rIdx * rowH;
    const hasClamped = r.below_floor_count > 0;

    // Row zebra background
    svg += `<rect x="0" y="${ry}" width="${totalW}" height="${rowH}" fill="${rIdx % 2 === 0 ? '#ffffff' : '#fafafa'}" opacity="0.6" />`;

    // Cancer Label with optional floor annotation asterisk
    svg += `<text x="10" y="${ry + 17}" font-size="11.5" font-weight="500" fill="#1e293b">${r.cancer_type} <tspan fill="#94a3b8" font-size="10" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace">(N=${r.total_patients})</tspan>${hasClamped ? ' <tspan fill="#b91c1c" font-size="11" font-weight="bold" title="' + r.below_floor_count + ' samples clamped (<15% floor)">*</tspan>' : ''}</text>`;

    // Any Target Cell
    const anyPct = r.any_elevated_pct || 0.0;
    const anyBg = getCellBg(anyPct);
    const anyFg = getCellTextColor(anyPct);
    svg += `
      <rect x="${labelW}" y="${ry + 2}" width="${anyColW - 4}" height="${rowH - 4}" rx="3" fill="${anyBg}" stroke="#e2e8f0" stroke-width="0.5" 
            onmouseover="showHeatmapTip(event, '${r.cancer_type}', 'ANY TARGET', '${r.any_elevated_count}', '${r.total_patients}', '${anyPct.toFixed(1)}%', '${r.below_floor_count}', '${r.uncalled_purity_count}')" 
            onmouseout="hideHeatmapTip()" 
            style="cursor:pointer;" onclick="sortHeatmapBy('any_target')">
      </rect>
      <text x="${labelW + (anyColW - 4) / 2}" y="${ry + 17}" text-anchor="middle" font-size="11" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" font-weight="600" fill="${anyFg}" pointer-events="none">${anyPct > 0 ? anyPct.toFixed(1) + '%' : '-'}</text>
    `;

    // Gene Cells
    genes.forEach((g, gIdx) => {
      const gx = labelW + anyColW + gIdx * colW;
      const gPct = (r.gene_pcts && r.gene_pcts[g] !== undefined) ? r.gene_pcts[g] : 0.0;
      const gCount = (r.gene_counts && r.gene_counts[g] !== undefined) ? r.gene_counts[g] : 0;
      const bg = getCellBg(gPct);
      const fg = getCellTextColor(gPct);

      svg += `
        <rect x="${gx}" y="${ry + 2}" width="${colW - 4}" height="${rowH - 4}" rx="3" fill="${bg}" stroke="#e2e8f0" stroke-width="0.5" 
              onmouseover="showHeatmapTip(event, '${r.cancer_type}', '${g}', '${gCount}', '${r.total_patients}', '${gPct.toFixed(1)}%', '${r.below_floor_count}', '${r.uncalled_purity_count}')" 
              onmouseout="hideHeatmapTip()" 
              style="cursor:pointer;" onclick="sortHeatmapBy('${g}')">
        </rect>
        <text x="${gx + (colW - 4) / 2}" y="${ry + 17}" text-anchor="middle" font-size="11" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" font-weight="${gPct > 0 ? '600' : 'normal'}" fill="${fg}" pointer-events="none">${gPct > 0 ? gPct.toFixed(1) + '%' : '-'}</text>
      `;
    });
  });

  svg += `</svg>`;
  container.innerHTML = svg;
}

function showHeatmapTip(e, cancer, target, count, total, pct, clamped, uncalled) {
  const tip = document.getElementById('heatmap-tooltip');
  if (!tip) return;
  let tipHtml = `
    <div style="font-weight:700; color:#38bdf8; margin-bottom:2px;">${cancer}</div>
    <div>Target: <strong>${target}</strong></div>
    <div>Elevated: <strong>${count} / ${total}</strong> patients (${pct})</div>
  `;
  if (parseInt(clamped) > 0) {
    tipHtml += `<div style="color:#fca5a5; font-size:10px; margin-top:2px;">Note: ${clamped} biopsies clamped to 15% purity floor</div>`;
  }
  if (parseInt(uncalled) > 0) {
    tipHtml += `<div style="color:#cbd5e1; font-size:10px;">Note: ${uncalled} biopsies with uncalled purity (at cohort median)</div>`;
  }
  tip.innerHTML = tipHtml;
  tip.style.left = (e.clientX + 14) + 'px';
  tip.style.top = (e.clientY + 14) + 'px';
  tip.style.display = 'block';
}

function hideHeatmapTip() {
  const tip = document.getElementById('heatmap-tooltip');
  if (tip) tip.style.display = 'none';
}

function profileSpecificGene(gene) {
  currentActiveGene = gene.toUpperCase();
  document.getElementById('gene-search').value = currentActiveGene;
  // global search box removed
  switchView('single-gene');
  loadGeneData();
}

function addGeneToAnalysis(gene) {
  const input = document.getElementById('multi-gene-input');
  const existing = input.value.split(',').map(s => s.trim().toUpperCase()).filter(Boolean);
  if (!existing.includes(gene.toUpperCase())) {
    existing.push(gene.toUpperCase());
    input.value = existing.join(', ');
  }
  switchView('multi-gene');
  runMultiGeneAnalysis();
}

function renderAdcCatalog() {
  const tbody = document.getElementById('adc-catalog-body');
  if (!tbody) return;
  let rows = '';
  ADC_CATALOG.forEach(item => {
    rows += `
      <tr>
        <td style="font-weight:700; font-family:var(--font-mono); color:var(--accent);">
          ${item.gene} <span style="font-size:10px; color:var(--text-muted); font-weight:normal;">(${item.alt})</span>
        </td>
        <td style="font-size:11.5px; color:var(--text-secondary); max-width:200px;">${item.name}</td>
        <td style="font-weight:500;">${item.drugs}</td>
        <td><span style="font-size:11px; font-family:var(--font-mono); color:#475569;">${item.warhead}</span></td>
        <td style="font-size:11.5px; color:var(--text-secondary);">${item.indications}</td>
        <td><span class="tag ${item.tag}">${item.status}</span></td>
        <td style="text-align:center; white-space:nowrap;">
          <button class="btn btn-outline btn-sm" onclick="addGeneToAnalysis('${item.gene}')" title="Add to Multi-Gene Panel">+ Add</button>
          <button class="btn btn-outline btn-sm" onclick="profileSpecificGene('${item.gene}')" title="Single Gene Profile">Profile</button>
        </td>
      </tr>
    `;
  });
  tbody.innerHTML = rows;
  makeTableSortable('adc-catalog-table');
}

// ---------------------------------------------------------------------------
// VIEW 3: DIFFERENTIAL CONTRAST STUDIO
// ---------------------------------------------------------------------------
async function runDifferentialExpression() {
  const btn = document.getElementById('run-contrast-btn');
  const sel = document.getElementById('contrast-select').value;
  const [col, grpA, grpB] = sel.split(':');
  document.getElementById('contrast-title').textContent = `${col.toUpperCase()} (${grpA} vs ${grpB})`;
  document.getElementById('diff-body').innerHTML = `<tr><td colspan="6" style="text-align:center; color:var(--text-muted); padding:2rem;">Computing genome-wide contrast across cohort...</td></tr>`;

  const svg = document.getElementById('volcano-svg');
  if (svg) {
    svg.innerHTML = '<text x="260" y="150" text-anchor="middle" fill="#94a3b8" font-size="12">Calculating genome-wide volcano distribution...</text>';
  }
  if (btn) {
    btn.disabled = true;
    btn.textContent = 'Computing...';
  }

  try {
    const res = await fetch(`/api/diff_expr?col=${col}&group_a=${grpA}&group_b=${grpB}`);
    const data = await res.json();
    if (data.error) {
      document.getElementById('diff-body').innerHTML = `<tr><td colspan="6" style="text-align:center; color:var(--text-muted); padding:1rem;">${data.error}</td></tr>`;
      if (svg) svg.innerHTML = `<text x="260" y="150" text-anchor="middle" fill="#ef4444" font-size="12">${data.error}</text>`;
      return;
    }

    // Render table
    let rows = '';
    data.top_genes.forEach(g => {
      const isUp = g.log2_fc > 0;
      rows += `
        <tr>
          <td style="font-weight:600; font-family:var(--font-mono); cursor:pointer; color:var(--accent);" onclick="profileSpecificGene('${g.gene}')">${g.gene}</td>
          <td class="num" style="color:${isUp ? 'var(--up-reg)' : 'var(--down-reg)'}; font-weight:600;">${isUp ? '+' : ''}${g.log2_fc.toFixed(2)}</td>
          <td class="num">${g.mean_a ? g.mean_a.toFixed(2) : '-'}</td>
          <td class="num">${g.mean_b ? g.mean_b.toFixed(2) : '-'}</td>
          <td class="num">${g.p_value.toExponential(2)}</td>
          <td class="num">${g.fdr.toExponential(2)}</td>
        </tr>
      `;
    });
    document.getElementById('diff-body').innerHTML = rows;
    makeTableSortable('diff-table');

    // Render Volcano SVG
    renderVolcano(data.volcano_points);
  } catch (e) {
    console.error(e);
  } finally {
    if (btn) {
      btn.disabled = false;
      btn.textContent = 'Run Contrast';
    }
  }
}

function renderVolcano(points) {
  const svg = document.getElementById('volcano-svg');
  if (!points || points.length === 0) {
    svg.innerHTML = '<text x="290" y="170" text-anchor="middle" fill="#64748b" font-size="13" font-family="system-ui, -apple-system, sans-serif">No genes reached significance threshold.</text>';
    return;
  }

  let maxAbsX = 2.5;
  let maxY = 8.0;
  points.forEach(p => {
    const px = Math.abs((p.x !== undefined) ? p.x : (p.log2_fc || 0));
    const py = (p.y !== undefined) ? p.y : (p.neg_log10_p || 0);
    if (px > maxAbsX) maxAbsX = Math.ceil(px * 1.15);
    if (py > maxY) maxY = Math.ceil(py * 1.15);
  });
  maxAbsX = Math.max(2.5, Math.min(8.0, maxAbsX));
  maxY = Math.max(6.0, Math.min(25.0, maxY));

  const w = 580, h = 340;
  const padLeft = 60, padRight = 30, padTop = 26, padBottom = 48;
  const plotW = w - padLeft - padRight;
  const plotH = h - padTop - padBottom;

  const scaleX = x => padLeft + ((x + maxAbsX) / (2 * maxAbsX)) * plotW;
  const scaleY = y => (padTop + plotH) - (Math.min(y, maxY) / maxY) * plotH;

  let out = '';

  // Subtle background grid
  const yStep = maxY <= 10 ? 2 : (maxY <= 20 ? 4 : 5);
  for (let y = yStep; y <= maxY; y += yStep) {
    const gy = scaleY(y);
    out += `<line x1="${padLeft}" y1="${gy.toFixed(1)}" x2="${(padLeft + plotW).toFixed(1)}" y2="${gy.toFixed(1)}" stroke="#f1f5f9" stroke-width="1" />`;
  }
  const xStep = maxAbsX <= 3.5 ? 1 : 2;
  for (let x = -Math.floor(maxAbsX); x <= Math.floor(maxAbsX); x += xStep) {
    if (x === 0) continue;
    const gx = scaleX(x);
    out += `<line x1="${gx.toFixed(1)}" y1="${padTop}" x2="${gx.toFixed(1)}" y2="${(padTop + plotH).toFixed(1)}" stroke="#f8fafc" stroke-width="1" />`;
  }

  // Fold change cutoffs (+/- 1.0)
  if (maxAbsX >= 1.2) {
    out += `<line x1="${scaleX(-1).toFixed(1)}" y1="${padTop}" x2="${scaleX(-1).toFixed(1)}" y2="${(padTop + plotH).toFixed(1)}" stroke="#cbd5e1" stroke-width="1" stroke-dasharray="2,2" />`;
    out += `<line x1="${scaleX(1).toFixed(1)}" y1="${padTop}" x2="${scaleX(1).toFixed(1)}" y2="${(padTop + plotH).toFixed(1)}" stroke="#cbd5e1" stroke-width="1" stroke-dasharray="2,2" />`;
  }

  // P=0.05 FDR cutoff line (y = 1.30)
  const yCutoff = scaleY(1.30);
  out += `<line x1="${padLeft}" y1="${yCutoff.toFixed(1)}" x2="${(padLeft + plotW).toFixed(1)}" y2="${yCutoff.toFixed(1)}" stroke="#d97706" stroke-width="1.2" stroke-dasharray="3,3" />`;
  out += `<text x="${(padLeft + plotW - 4).toFixed(1)}" y="${(yCutoff - 4).toFixed(1)}" text-anchor="end" font-size="10" font-family="system-ui, -apple-system, sans-serif" font-weight="500" fill="#d97706">P = 0.05 cutoff</text>`;

  // Center vertical reference line (x = 0)
  out += `<line x1="${scaleX(0).toFixed(1)}" y1="${padTop}" x2="${scaleX(0)}" y2="${(padTop + plotH).toFixed(1)}" stroke="#94a3b8" stroke-width="1.2" stroke-dasharray="3,3" />`;

  // Solid Axes lines
  out += `<line x1="${padLeft}" y1="${(padTop + plotH).toFixed(1)}" x2="${(padLeft + plotW).toFixed(1)}" y2="${(padTop + plotH).toFixed(1)}" stroke="#475569" stroke-width="1.5" />`;
  out += `<line x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${(padTop + plotH).toFixed(1)}" stroke="#475569" stroke-width="1.5" />`;

  // X ticks and numeric labels
  for (let x = -Math.floor(maxAbsX); x <= Math.floor(maxAbsX); x += xStep) {
    const sx = scaleX(x);
    out += `<line x1="${sx.toFixed(1)}" y1="${(padTop + plotH).toFixed(1)}" x2="${sx.toFixed(1)}" y2="${(padTop + plotH + 5).toFixed(1)}" stroke="#475569" stroke-width="1.2" />`;
    const label = x > 0 ? '+' + x : (x === 0 ? '0' : '' + x);
    out += `<text x="${sx.toFixed(1)}" y="${(padTop + plotH + 17).toFixed(1)}" text-anchor="middle" font-size="11" font-weight="500" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" fill="#334155">${label}</text>`;
  }

  // Y ticks and numeric labels
  for (let y = 0; y <= maxY; y += yStep) {
    const sy = scaleY(y);
    out += `<line x1="${(padLeft - 5).toFixed(1)}" y1="${sy.toFixed(1)}" x2="${padLeft.toFixed(1)}" y2="${sy.toFixed(1)}" stroke="#475569" stroke-width="1.2" />`;
    out += `<text x="${(padLeft - 8).toFixed(1)}" y="${(sy + 4).toFixed(1)}" text-anchor="end" font-size="11" font-weight="500" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" fill="#334155">${y}</text>`;
  }

  // Axis Titles
  out += `<text x="${(padLeft + plotW / 2).toFixed(1)}" y="${(h - 10).toFixed(1)}" text-anchor="middle" font-size="12" font-weight="600" font-family="system-ui, -apple-system, sans-serif" fill="#0f172a">Log₂ Fold Change (Group A vs Group B)</text>`;
  out += `<text x="${(-(padTop + plotH / 2)).toFixed(1)}" y="18" transform="rotate(-90)" text-anchor="middle" font-size="12" font-weight="600" font-family="system-ui, -apple-system, sans-serif" fill="#0f172a">-Log₁₀ (Welch's P-value)</text>`;

  // Data Points
  const sortedPoints = [...points].sort((a, b) => {
    const sigA = (Math.abs(a.log2_fc || a.x || 0) >= 1.0 && (a.fdr || 1.0) <= 0.05) ? 1 : 0;
    const sigB = (Math.abs(b.log2_fc || b.x || 0) >= 1.0 && (b.fdr || 1.0) <= 0.05) ? 1 : 0;
    return sigA - sigB;
  });

  sortedPoints.forEach(p => {
    const px = (p.x !== undefined) ? p.x : (p.log2_fc || 0);
    const py = (p.y !== undefined) ? p.y : (p.neg_log10_p || 0);
    const cx = scaleX(px);
    const cy = scaleY(py);

    const isUp = (px >= 1.0 && (p.fdr || 1.0) <= 0.05);
    const isDown = (px <= -1.0 && (p.fdr || 1.0) <= 0.05);

    let color = '#94a3b8';
    let r = 2.2;
    let op = 0.45;
    if (isUp) {
      color = '#dc2626';
      r = 3.8;
      op = 0.85;
    } else if (isDown) {
      color = '#2563eb';
      r = 3.8;
      op = 0.85;
    }

    out += `<circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${r}" fill="${color}" opacity="${op}" style="cursor:pointer;" onclick="profileSpecificGene('${p.gene}')"><title>${p.gene}: Log2FC ${px > 0 ? '+' : ''}${px.toFixed(2)}, -Log10P ${py.toFixed(2)}, FDR ${(p.fdr || 1.0).toExponential(2)}</title></circle>`;
  });

  // Top 5 UP & Top 5 DOWN labels with white halo
  const topUp = points.filter(p => (p.log2_fc || p.x || 0) >= 1.0 && (p.fdr || 1.0) <= 0.05)
                      .sort((a, b) => (b.neg_log10_p || b.y || 0) - (a.neg_log10_p || a.y || 0)).slice(0, 5);
  const topDown = points.filter(p => (p.log2_fc || p.x || 0) <= -1.0 && (p.fdr || 1.0) <= 0.05)
                        .sort((a, b) => (b.neg_log10_p || b.y || 0) - (a.neg_log10_p || a.y || 0)).slice(0, 5);

  [...topUp, ...topDown].forEach(p => {
    const px = (p.x !== undefined) ? p.x : (p.log2_fc || 0);
    const py = (p.y !== undefined) ? p.y : (p.neg_log10_p || 0);
    const cx = scaleX(px);
    const cy = scaleY(py);
    const anchor = px > 0 ? 'start' : 'end';
    const ox = px > 0 ? 6 : -6;
    out += `<text x="${(cx + ox).toFixed(1)}" y="${(cy - 4).toFixed(1)}" text-anchor="${anchor}" font-size="9.5" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" font-weight="700" stroke="#ffffff" stroke-width="3" stroke-linejoin="round" fill="none" style="pointer-events:none;">${p.gene}</text>`;
    out += `<text x="${(cx + ox).toFixed(1)}" y="${(cy - 4).toFixed(1)}" text-anchor="${anchor}" font-size="9.5" font-family="ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace" font-weight="700" fill="#0f172a" style="pointer-events:none;">${p.gene}</text>`;
  });

  svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
  svg.innerHTML = out;
}

function exportVolcanoStandaloneHtml() {
  const contrastTitle = (document.getElementById('contrast-title') || {}).textContent || 'Differential Contrast';
  const svgEl = document.getElementById('volcano-svg');
  const tableEl = document.getElementById('diff-table');
  if (!svgEl || !tableEl) return;

  const svgHtml = svgEl.outerHTML;
  const tableHtml = tableEl.outerHTML;

  const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Volcano Plot & Differential Expression: \${contrastTitle}</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f8fafc; color: #0f172a; margin: 0; padding: 2rem; }
  .container { max-width: 1200px; margin: 0 auto; background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 2rem; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05); }
  h1 { font-size: 1.4rem; font-weight: 700; margin-bottom: 0.4rem; color: #0f172a; }
  .subtitle { font-size: 0.9rem; color: #475569; margin-bottom: 1.5rem; }
  .grid { display: grid; grid-template-columns: 1.1fr 1fr; gap: 1.5rem; margin-top: 1rem; }
  @media (max-width: 900px) { .grid { grid-template-columns: 1fr; } }
  .card { border: 1px solid #e2e8f0; border-radius: 6px; padding: 1.25rem; background: #fafafa; }
  .card-title { font-weight: 600; font-size: 0.95rem; margin-bottom: 0.75rem; color: #1e293b; }
  table { width: 100%; border-collapse: collapse; font-size: 12px; }
  th, td { padding: 0.5rem 0.6rem; text-align: left; border-bottom: 1px solid #e2e8f0; }
  th { background: #f1f5f9; font-weight: 600; color: #334155; cursor: pointer; user-select: none; }
  th:hover { background: #e2e8f0; }
  th.sort-asc::after { content: " ▲"; font-size: 9px; color: #2563eb; }
  th.sort-desc::after { content: " ▼"; font-size: 9px; color: #2563eb; }
  td.num, th.num { text-align: right; }
  tr:hover td { background: #f8fafc; }
</style>
</head>
<body>
<div class="container">
  <h1>Clinical Differential Contrast: \${contrastTitle}</h1>
  <div class="subtitle">Genome-wide Welch's t-test fold change & Benjamini-Hochberg FDR significance.</div>
  <div class="grid">
    <div class="card">
      <div class="card-title">Volcano Distribution</div>
      <div style="background:#ffffff; border-radius:4px; padding:0.5rem;">
        \${svgHtml}
      </div>
    </div>
    <div class="card">
      <div class="card-title">Top Differentially Expressed Genes (Click headers to sort)</div>
      <div style="max-height: 400px; overflow-y: auto; background:#ffffff; border:1px solid #e2e8f0; border-radius:4px;">
        \${tableHtml}
      </div>
    </div>
  </div>
</div>
<' + 'script>
document.querySelectorAll('table th').forEach((th, idx) => {
  let asc = false;
  th.addEventListener('click', () => {
    const table = th.closest('table');
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    asc = !asc;
    table.querySelectorAll('th').forEach(h => h.classList.remove('sort-asc', 'sort-desc'));
    th.classList.add(asc ? 'sort-asc' : 'sort-desc');
    rows.sort((r1, r2) => {
      const v1 = r1.children[idx].innerText.trim();
      const v2 = r2.children[idx].innerText.trim();
      const n1 = parseFloat(v1.replace(/[^0-9.-]/g, ''));
      const n2 = parseFloat(v2.replace(/[^0-9.-]/g, ''));
      if (!isNaN(n1) && !isNaN(n2)) return asc ? n1 - n2 : n2 - n1;
      return asc ? v1.localeCompare(v2) : v2.localeCompare(v1);
    });
    rows.forEach(r => tbody.appendChild(r));
  });
});
<' + '/script>
</body>
</html>`;

  const blob = new Blob([html], { type: 'text/html' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `volcano_\${contrastTitle.toLowerCase().replace(/[^a-z0-9]+/g, '_')}.html`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------------------------------------------------------------------------
// VIEW 4: COHORT UMAP LANDSCAPE
// ---------------------------------------------------------------------------
let activeUmapFilter = null;
let activeHighlightSample = '';

async function loadCohortUmap() {
  const container = document.getElementById('umap-container');
  const svg = document.getElementById('umap-svg');
  if (svg) {
    svg.innerHTML = '<text x="370" y="280" text-anchor="middle" fill="#64748b" font-size="13" font-family="system-ui, -apple-system, sans-serif">Loading 1,952 cohort UMAP coordinates & metadata...</text>';
  }

  try {
    const res = await fetch('/api/cohort_umap');
    const data = await res.json();
    window.cohortUmapData = data;
    window.cohortUmapData.sampleMap = {};
    if (data.samples) {
      data.samples.forEach(s => {
        window.cohortUmapData.sampleMap[s.id] = s;
      });
    }
    renderUmapPlot();
  } catch (e) {
    console.error("UMAP error:", e);
    if (svg) {
      svg.innerHTML = `<text x="370" y="280" text-anchor="middle" fill="#ef4444" font-size="13">Error loading UMAP: ${e.message}</text>`;
    }
  }
}

// 32-color curated oncological palette
const UMAP_PALETTE = [
  '#2563eb', '#dc2626', '#059669', '#d97706', '#7c3aed', '#db2777', '#0891b2', '#4f46e5',
  '#ea580c', '#16a34a', '#9333ea', '#e11d48', '#0284c7', '#ca8a04', '#65a30d', '#c026d3',
  '#475569', '#0d9488', '#b91c1c', '#1e40af', '#15803d', '#b45309', '#6d28d9', '#be185d',
  '#0369a1', '#854d0e', '#3f6212', '#a21caf', '#334155', '#0f766e', '#991b1b', '#1e3a8a'
];

function getPurityColor(purity) {
  if (purity === null || purity === undefined || purity <= 0) {
    return { bg: '#cbd5e1', stroke: '#94a3b8', label: 'Uncalled' };
  }
  if (purity < 0.15) {
    return { bg: '#fee2e2', stroke: '#dc2626', label: '< 15% Floor', isFloor: true };
  }
  // Continuous White to Deep Blue (0.15 -> 1.00)
  const t = Math.min(1.0, Math.max(0.0, (purity - 0.15) / 0.85));
  // White: rgb(255, 255, 255) -> Deep Blue: rgb(23, 37, 84)
  const r = Math.round(255 - t * (255 - 23));
  const g = Math.round(255 - t * (255 - 37));
  const b = Math.round(255 - t * (255 - 84));
  return { bg: `rgb(${r},${g},${b})`, stroke: t < 0.25 ? '#94a3b8' : 'rgba(0,0,0,0.15)' };
}

function getGeneExprColor(val, maxVal) {
  const v = Math.max(0, val || 0);
  const m = Math.max(1, maxVal);
  const t = Math.min(1.0, v / m);
  if (t <= 0) return { bg: '#ffffff', stroke: '#cbd5e1' };
  // White -> Crimson Red (255,255,255) -> (136, 19, 55)
  const r = Math.round(255 - t * (255 - 136));
  const g = Math.round(255 - t * (255 - 19));
  const b = Math.round(255 - t * (255 - 55));
  return { bg: `rgb(${r},${g},${b})`, stroke: t < 0.2 ? '#cbd5e1' : 'rgba(0,0,0,0.15)' };
}

function renderUmapPlot() {
  const data = window.cohortUmapData;
  if (!data || !data.samples || data.samples.length === 0) return;

  const projEl = document.getElementById('umap-projection');
  const proj = projEl ? projEl.value : 'tsne';
  const mode = document.getElementById('umap-color-by').value;
  const annotateFloor = document.getElementById('umap-floor-toggle').checked;
  const searchSample = (document.getElementById('umap-sample-search').value || '').trim().toUpperCase();

  const svg = document.getElementById('umap-svg');
  const w = 740, h = 560;
  const pad = 35;

  // Assign x, y based on selected projection
  data.samples.forEach(s => {
    if (proj === 'tsne') {
      s.x = s.tsne_x !== undefined ? s.tsne_x : s.umap_x;
      s.y = s.tsne_y !== undefined ? s.tsne_y : s.umap_y;
    } else if (proj === 'pca12') {
      s.x = s.pca12_x !== undefined ? s.pca12_x : s.umap_x;
      s.y = s.pca12_y !== undefined ? s.pca12_y : s.umap_y;
    } else if (proj === 'pca13') {
      s.x = s.pca13_x !== undefined ? s.pca13_x : s.umap_x;
      s.y = s.pca13_y !== undefined ? s.pca13_y : s.umap_y;
    } else {
      s.x = s.umap_x;
      s.y = s.umap_y;
    }
  });

  // Coordinate Bounds
  let minX = -10, maxX = 10, minY = -10, maxY = 10;
  data.samples.forEach(s => {
    if (s.x < minX) minX = s.x;
    if (s.x > maxX) maxX = s.x;
    if (s.y < minY) minY = s.y;
    if (s.y > maxY) maxY = s.y;
  });
  const spanX = maxX - minX || 1;
  const spanY = maxY - minY || 1;

  const scaleX = x => pad + ((x - minX) / spanX) * (w - 2 * pad);
  const scaleY = y => (h - pad) - ((y - minY) / spanY) * (h - 2 * pad);

  // Update Card Header titles
  const cardTitle = document.getElementById('umap-card-title');
  const cardCaption = document.getElementById('umap-card-caption');
  const colorbarBox = document.getElementById('umap-colorbar-container');
  const categoryBox = document.getElementById('umap-category-list');

  let projName = 't-SNE Landscape';
  let xLabel = 't-SNE 1', yLabel = 't-SNE 2';
  if (proj === 'pca12') {
    projName = 'Principal Component Analysis (PC1 vs PC2)';
    xLabel = 'PC 1 (17.6% Variance)';
    yLabel = 'PC 2 (16.4% Variance)';
  } else if (proj === 'pca13') {
    projName = 'Principal Component Analysis (PC1 vs PC3)';
    xLabel = 'PC 1 (17.6% Variance)';
    yLabel = 'PC 3 (11.4% Variance)';
  } else if (proj === 'umap') {
    projName = 'UMAP Manifold';
    xLabel = 'UMAP 1';
    yLabel = 'UMAP 2';
  }

  let isContinuous = false;
  let colorMap = {};

  if (mode === 'cancer_type') {
    cardTitle.textContent = `${projName}: Clinical Primary Tumour Type`;
    cardCaption.textContent = '2D projection of 1,952 biopsies coloured by clinical diagnosis (32 cohorts)';
    data.cancer_types.forEach((c, idx) => {
      colorMap[c] = UMAP_PALETTE[idx % UMAP_PALETTE.length];
    });
  } else if (mode === 'cuppa') {
    cardTitle.textContent = `${projName}: CUPPA Predicted Lineage`;
    cardCaption.textContent = 'Genomic tissue-of-origin classifier prediction (HMF CUPPA algorithm across 1,943 biopsies)';
    data.cuppa_types.forEach((c, idx) => {
      colorMap[c] = UMAP_PALETTE[idx % UMAP_PALETTE.length];
    });
  } else if (mode === 'purity') {
    isContinuous = true;
    cardTitle.textContent = `${projName}: Tumour Cellularity (PURPLE Purity)`;
    cardCaption.textContent = 'Continuous White to Deep Blue gradient (15% floor to 100% pure tumour content)';
  } else if (mode === 'specimen_type') {
    cardTitle.textContent = `${projName}: Specimen Clinical Category`;
    cardCaption.textContent = 'Primary tumour vs Metastatic lesion vs Local recurrence';
    colorMap['Primary'] = '#059669';
    colorMap['Metastasis'] = '#7c3aed';
    colorMap['Local Recurrence'] = '#ea580c';
    colorMap['Other'] = '#94a3b8';
    colorMap['Unknown'] = '#94a3b8';
  } else if (mode.startsWith('gene:')) {
    isContinuous = true;
    const gene = mode.split(':')[1];
    cardTitle.textContent = `${projName}: ${gene} Expression`;
    cardCaption.textContent = `Normalized transcript abundance (Log2 Adj TPM) across 1,952 biopsies`;
  }

  // Update Legend Panel
  if (isContinuous) {
    categoryBox.style.display = 'none';
    colorbarBox.style.display = 'block';

    if (mode === 'purity') {
      document.getElementById('umap-colorbar-labels').innerHTML = '<span>15% (Purity Floor)</span><span>100% (Pure Tumour)</span>';
      document.getElementById('umap-colorbar-gradient').style.background = 'linear-gradient(to right, #ffffff, #bfdbfe, #60a5fa, #2563eb, #172554)';
      document.getElementById('umap-colorbar-note').innerHTML = `
        <div style="display:flex; align-items:center; gap:6px; margin-top:4px;">
          <span style="width:10px; height:10px; border-radius:50%; background:#fee2e2; border:1.5px solid #dc2626; display:inline-block;"></span>
          <span>&lt; 15% Cellularity Floor (Clamped)</span>
        </div>
        <div style="display:flex; align-items:center; gap:6px; margin-top:3px;">
          <span style="width:10px; height:10px; border-radius:50%; background:#cbd5e1; border:1px solid #94a3b8; display:inline-block;"></span>
          <span>Uncalled Purity (Median 40%)</span>
        </div>
      `;
    } else {
      const gene = mode.split(':')[1];
      let maxVal = 1.0;
      data.samples.forEach(s => {
        const v = s[gene.toLowerCase()] || 0.0;
        if (v > maxVal) maxVal = v;
      });
      document.getElementById('umap-colorbar-labels').innerHTML = `<span>0.0 TPM</span><span>${maxVal.toFixed(1)} Log2 TPM</span>`;
      document.getElementById('umap-colorbar-gradient').style.background = 'linear-gradient(to right, #ffffff, #fecdd3, #f43f5e, #be123c, #881337)';
      document.getElementById('umap-colorbar-note').innerHTML = `Continuous White to Crimson gradient for <strong>${gene}</strong> abundance.`;
    }
  } else {
    colorbarBox.style.display = 'none';
    categoryBox.style.display = 'flex';

    let catHtml = '';
    const items = mode === 'cancer_type' ? data.cancer_types : (mode === 'cuppa' ? data.cuppa_types : data.specimen_types);
    const countsMap = mode === 'cancer_type' ? data.cancer_counts : (mode === 'cuppa' ? data.cuppa_counts : data.specimen_counts);

    items.forEach(cat => {
      const color = colorMap[cat] || '#94a3b8';
      const count = countsMap ? (countsMap[cat] || 0) : 0;
      const isSelected = activeUmapFilter === cat;

      catHtml += `
        <div class="umap-leg-item ${isSelected ? 'active' : ''}" onclick="toggleUmapFilter('${cat}')" style="display:flex; align-items:center; justify-content:space-between; padding:3px 6px; border-radius:4px; font-size:11px; cursor:pointer; background:${isSelected ? '#eff6ff' : 'transparent'};">
          <div style="display:flex; align-items:center; gap:6px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;">
            <span style="width:9px; height:9px; border-radius:50%; background:${color}; flex-shrink:0;"></span>
            <span style="font-weight:${isSelected ? '700' : '500'}; color:#1e293b; overflow:hidden; text-overflow:ellipsis;">${cat}</span>
          </div>
          <span style="font-size:10px; color:var(--text-muted); font-family:var(--font-mono);">${count}</span>
        </div>
      `;
    });
    categoryBox.innerHTML = catHtml;
  }

  // Draw SVG Canvas
  let out = '';

  // Background Grid and Axis guidelines
  out += `<rect width="${w}" height="${h}" fill="#ffffff" />`;
  out += `<line x1="${pad}" y1="${h - pad}" x2="${w - pad}" y2="${h - pad}" stroke="#e2e8f0" stroke-width="1.2" />`;
  out += `<line x1="${pad}" y1="${pad}" x2="${pad}" y2="${h - pad}" stroke="#e2e8f0" stroke-width="1.2" />`;
  out += `<text x="${w / 2}" y="${h - 10}" text-anchor="middle" font-size="11" font-weight="600" fill="#64748b" font-family="system-ui, sans-serif">${xLabel}</text>`;
  out += `<text x="${-(h / 2)}" y="16" transform="rotate(-90)" text-anchor="middle" font-size="11" font-weight="600" fill="#64748b" font-family="system-ui, sans-serif">${yLabel}</text>`;

  // Find max gene expr if in gene mode
  let maxGeneVal = 10.0;
  if (mode.startsWith('gene:')) {
    const gene = mode.split(':')[1].toLowerCase();
    data.samples.forEach(s => {
      const v = s[gene] || 0.0;
      if (v > maxGeneVal) maxGeneVal = v;
    });
  }

  // Draw Sample Points
  data.samples.forEach(s => {
    const cx = scaleX(s.x).toFixed(1);
    const cy = scaleY(s.y).toFixed(1);

    let color = '#94a3b8';
    let strokeColor = 'rgba(0,0,0,0.15)';
    let strokeWidth = '0.5';
    let isFloorSample = false;

    if (mode === 'cancer_type') {
      color = colorMap[s.cancer] || '#94a3b8';
    } else if (mode === 'cuppa') {
      color = colorMap[s.cuppa] || '#94a3b8';
    } else if (mode === 'purity') {
      const pInfo = getPurityColor(s.purity);
      color = pInfo.bg;
      strokeColor = pInfo.stroke;
      if (pInfo.isFloor && annotateFloor) {
        isFloorSample = true;
        strokeWidth = '1.8';
      }
    } else if (mode === 'specimen_type') {
      color = colorMap[s.specimen] || '#94a3b8';
    } else if (mode.startsWith('gene:')) {
      const gene = mode.split(':')[1].toLowerCase();
      const gInfo = getGeneExprColor(s[gene], maxGeneVal);
      color = gInfo.bg;
      strokeColor = gInfo.stroke;
    }

    // Filter check
    let op = 0.82;
    let r = 3.4;
    if (activeUmapFilter) {
      const matches = (mode === 'cancer_type' && s.cancer === activeUmapFilter) ||
                      (mode === 'cuppa' && s.cuppa === activeUmapFilter) ||
                      (mode === 'specimen_type' && s.specimen === activeUmapFilter);
      if (!matches) {
        op = 0.12;
        r = 2.0;
      } else {
        r = 4.2;
        op = 1.0;
      }
    }

    // Highlight search match
    const isSearchMatch = searchSample && s.id.toUpperCase().includes(searchSample);
    if (isSearchMatch) {
      r = 6.0;
      op = 1.0;
      strokeColor = '#f59e0b';
      strokeWidth = '2.5';
    }

    out += `
      <circle cx="${cx}" cy="${cy}" r="${r}" fill="${color}" stroke="${strokeColor}" stroke-width="${strokeWidth}" opacity="${op}" style="cursor:pointer;"
              onmouseover="showUmapTip(event, '${s.id}')"
              onmouseout="hideUmapTip()">
      </circle>
    `;

    // Floor highlight halo ring
    if (isFloorSample && annotateFloor && (!activeUmapFilter || op > 0.5)) {
      out += `<circle cx="${cx}" cy="${cy}" r="5.2" fill="none" stroke="#dc2626" stroke-width="1.2" stroke-dasharray="2,2" opacity="0.9" pointer-events="none" />`;
    }
  });

  svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
  svg.innerHTML = out;
}

function toggleUmapFilter(cat) {
  if (activeUmapFilter === cat) {
    activeUmapFilter = null;
  } else {
    activeUmapFilter = cat;
  }
  renderUmapPlot();
}

function highlightUmapSample() {
  renderUmapPlot();
}

function showUmapTip(e, id) {
  const tip = document.getElementById('umap-tooltip');
  if (!tip || !window.cohortUmapData || !window.cohortUmapData.sampleMap) return;
  const s = window.cohortUmapData.sampleMap[id];
  if (!s) return;

  const purityText = s.purity !== null && s.purity !== undefined && s.purity > 0
    ? `${(s.purity * 100).toFixed(0)}% ${s.purity < 0.15 ? '<span style="color:#fca5a5; font-weight:700;">(< 15% floor)</span>' : ''}`
    : '<span style="color:#94a3b8;">Uncalled (Median 40%)</span>';

  tip.innerHTML = `
    <div style="font-weight:700; color:#38bdf8; font-size:12px; margin-bottom:3px; font-family:ui-monospace, monospace;">${s.id}</div>
    <div>Clinical Cohort: <strong>${s.cancer}</strong></div>
    <div>CUPPA Lineage: <strong>${s.cuppa}</strong> <span style="font-size:10px; color:#94a3b8;">(${(s.cuppa_p * 100).toFixed(0)}% conf)</span></div>
    ${s.cuppa_top ? `<div style="margin-top:2px; font-size:9.5px; color:#cbd5e1; line-height:1.2;">Top: <em>${s.cuppa_top}</em></div>` : ''}
    <div style="margin-top:2px;">Tumour Purity: <strong>${purityText}</strong></div>
    <div>Specimen Category: <strong>${s.specimen}</strong></div>
    <div style="border-top:1px solid #334155; margin-top:5px; padding-top:4px; font-size:10.5px;">
      <div style="color:#94a3b8; font-weight:600; margin-bottom:1px;">Horizon-Beacon Target Abundance (Log2 TPM):</div>
      <div style="display:grid; grid-template-columns: 1fr 1fr; gap:2px 8px; font-family:ui-monospace, monospace;">
        <span>ERBB2: <strong>${s.erbb2}</strong></span>
        <span>TACSTD2: <strong>${s.tacstd2}</strong></span>
        <span>CLDN18: <strong>${s.cldn18}</strong></span>
        <span>FOLR1: <strong>${s.folr1}</strong></span>
        <span>VTCN1: <strong>${s.vtcn1}</strong></span>
      </div>
    </div>
  `;

  const container = document.getElementById('umap-container');
  const rect = container.getBoundingClientRect();
  const x = e.clientX - rect.left + 14;
  const y = e.clientY - rect.top + 14;

  tip.style.left = Math.min(x, rect.width - 310) + 'px';
  tip.style.top = Math.min(y, rect.height - 210) + 'px';
  tip.style.display = 'block';
}

function hideUmapTip() {
  const tip = document.getElementById('umap-tooltip');
  if (tip) tip.style.display = 'none';
}

function exportUmapStandaloneHtml() {
  const data = window.cohortUmapData;
  if (!data) return;

  const mode = document.getElementById('umap-color-by').value;
  const svgEl = document.getElementById('umap-svg');
  if (!svgEl) return;

  const svgHtml = svgEl.outerHTML;

  const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Cohort UMAP Landscape: 1,952 Biopsies</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f8fafc; color: #0f172a; margin: 0; padding: 2rem; }
  .container { max-width: 1100px; margin: 0 auto; background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 2rem; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05); }
  h1 { font-size: 1.4rem; font-weight: 700; margin-bottom: 0.4rem; color: #0f172a; }
  .subtitle { font-size: 0.9rem; color: #475569; margin-bottom: 1.5rem; }
  .card { border: 1px solid #e2e8f0; border-radius: 6px; padding: 1.25rem; background: #ffffff; }
  svg { width: 100%; height: auto; display: block; border-radius: 4px; border: 1px solid #e2e8f0; }
  .footer-note { font-size: 11px; color: #64748b; margin-top: 1rem; }
</style>
</head>
<body>
<div class="container">
  <h1>Whole-Cohort Landscape (1,952 Biopsies)</h1>
  <div class="subtitle">Dimension: ${mode}</div>
  <div class="card">
    ${svgHtml}
  </div>
</div>
</body>
</html>`;

  const blob = new Blob([html], { type: 'text/html' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `cohort_umap_${mode.replace(/[^a-z0-9]+/g, '_')}.html`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------------------------------------------------------------------------
// VIEW 6: TUMOUR MUTATION BURDEN (TMB) & MUTATION COUNTS
// ---------------------------------------------------------------------------
let window_cohortTmbData = null;

async function loadTmbData() {
  const svg = document.getElementById('tmb-svg');
  if (svg) {
    svg.innerHTML = '<text x="500" y="260" text-anchor="middle" fill="#64748b" font-size="13" font-family="system-ui, sans-serif">Loading mutation burden data across 1,714 biopsies...</text>';
  }

  try {
    const res = await fetch('/api/tmb_data');
    const data = await res.json();
    window.cohortTmbData = data;
    window.cohortTmbData.sampleMap = {};
    if (data.samples) {
      data.samples.forEach(s => {
        window.cohortTmbData.sampleMap[s.id] = s;
      });
    }
    renderTmbPlot();
  } catch (e) {
    console.error("TMB error:", e);
    if (svg) {
      svg.innerHTML = `<text x="500" y="260" text-anchor="middle" fill="#ef4444" font-size="13">Error loading TMB data: ${e.message}</text>`;
    }
  }
}

function renderTmbPlot() {
  const data = window.cohortTmbData;
  if (!data || !data.samples || data.samples.length === 0) return;

  const groupBy = (document.getElementById('tmb-group-by') && document.getElementById('tmb-group-by').value) || 'cancer_type';
  const metric = (document.getElementById('tmb-metric') && document.getElementById('tmb-metric').value) || 'tmb';
  const scale = (document.getElementById('tmb-scale') && document.getElementById('tmb-scale').value) || 'log10';
  const floorVal = parseFloat(document.getElementById('tmb-purity-floor').value) || 0.0;
  const excludeFloor = document.getElementById('tmb-purity-exclude').checked;
  const showThreshold = document.getElementById('tmb-show-threshold').checked;

  const cardTitle = document.getElementById('tmb-card-title');
  const cardCaption = document.getElementById('tmb-card-caption');
  const statBadge = document.getElementById('tmb-stat-badge');

  // Filter samples based on purity
  let filteredSamples = [];
  let excludedCount = 0;

  data.samples.forEach(s => {
    if (s.tmb === null || s.tmb === undefined || s.tmb < 0) return;
    const isBelowFloor = s.purity !== null && s.purity !== undefined && s.purity < floorVal;
    if (isBelowFloor && excludeFloor) {
      excludedCount++;
      return;
    }
    filteredSamples.push(s);
  });

  statBadge.textContent = `${filteredSamples.length.toLocaleString()} biopsies (${excludedCount} excluded < ${Math.round(floorVal * 100)}% purity)`;

  // Update headers
  let metricLabel = 'TMB (Mutations / Mb)';
  let yAxisTitle = 'TMB (Mutations / Mb)';
  if (metric === 'mut_count') {
    metricLabel = 'Mutations per Sample (WGS Total)';
    yAxisTitle = scale === 'log2' ? 'log2 (Mutations per sample)' : 'Mutations per Sample';
  } else if (metric === 'coding_mut') {
    metricLabel = 'Coding Mutations (Exome / TML)';
    yAxisTitle = scale === 'log2' ? 'log2 (Coding mutations per sample)' : 'Coding Mutations';
  } else {
    yAxisTitle = scale === 'log2' ? 'log2 (TMB per Mb)' : 'TMB (Mutations / Mb)';
  }

  if (groupBy === 'pole') {
    cardTitle.textContent = 'POLE';
    cardCaption.textContent = `Somatic mutation load across cancer types stratified by POLE signature activity (MT vs WT). Monotonic S-curves with median bars.`;
  } else if (groupBy === 'msi') {
    cardTitle.textContent = 'Microsatellite Instability (MSI vs MSS)';
    cardCaption.textContent = `Somatic mutation burden stratified by MSI status across cancer cohorts.`;
  } else if (groupBy === 'cuppa') {
    cardTitle.textContent = 'CUPPA Lineage Mutation Burden';
    cardCaption.textContent = `Somatic mutation burden across 33 genomic tissue-of-origin predicted lineages.`;
  } else {
    cardTitle.textContent = 'Pan-Cancer Tumour Mutation Burden';
    cardCaption.textContent = `Pan-cancer distribution of somatic mutation burden across 32 cohorts.`;
  }

  // Build Columns
  const columns = [];

  if (groupBy === 'pole') {
    const cancerMap = {};
    filteredSamples.forEach(s => {
      if (!cancerMap[s.cancer]) cancerMap[s.cancer] = { mt: [], wt: [] };
      if (s.is_pole_mt) cancerMap[s.cancer].mt.push(s);
      else cancerMap[s.cancer].wt.push(s);
    });

    const sortedCancers = Object.keys(cancerMap).filter(c => (cancerMap[c].mt.length + cancerMap[c].wt.length) >= 4);
    sortedCancers.sort((a, b) => (cancerMap[b].mt.length + cancerMap[b].wt.length) - (cancerMap[a].mt.length + cancerMap[a].wt.length));

    sortedCancers.forEach(c => {
      const g = cancerMap[c];
      const shortC = c.replace(/ Cancer| Tumor| Carcinoma/g, '');
      if (g.mt.length > 0) {
        columns.push({
          id: `${c}_MT`,
          label: `${shortC}_MT`,
          cohort: c,
          subLabel: 'POLE MT',
          samples: g.mt,
          color: '#d97706',
          isMt: true
        });
      }
      if (g.wt.length > 0) {
        columns.push({
          id: `${c}_WT`,
          label: `${shortC}_WT`,
          cohort: c,
          subLabel: 'POLE WT',
          samples: g.wt,
          color: '#0284c7',
          isMt: false
        });
      }
    });
  } else if (groupBy === 'msi') {
    const cancerMap = {};
    filteredSamples.forEach(s => {
      if (!cancerMap[s.cancer]) cancerMap[s.cancer] = { msi: [], mss: [] };
      if (s.msi === 'MSI') cancerMap[s.cancer].msi.push(s);
      else cancerMap[s.cancer].mss.push(s);
    });

    const sortedCancers = Object.keys(cancerMap).filter(c => (cancerMap[c].msi.length + cancerMap[c].mss.length) >= 4);
    sortedCancers.sort((a, b) => (cancerMap[b].msi.length + cancerMap[b].mss.length) - (cancerMap[a].msi.length + cancerMap[a].mss.length));

    sortedCancers.forEach(c => {
      const g = cancerMap[c];
      const shortC = c.replace(/ Cancer| Tumor| Carcinoma/g, '');
      if (g.msi.length > 0) {
        columns.push({
          id: `${c}_MSI`,
          label: `${shortC}_MSI`,
          cohort: c,
          subLabel: 'MSI-H',
          samples: g.msi,
          color: '#dc2626',
          isMt: true
        });
      }
      if (g.mss.length > 0) {
        columns.push({
          id: `${c}_MSS`,
          label: `${shortC}_MSS`,
          cohort: c,
          subLabel: 'MSS',
          samples: g.mss,
          color: '#2563eb',
          isMt: false
        });
      }
    });
  } else {
    const groupMap = {};
    filteredSamples.forEach(s => {
      const k = groupBy === 'cuppa' ? s.cuppa : (groupBy === 'specimen' ? s.specimen : s.cancer);
      if (!groupMap[k]) groupMap[k] = [];
      groupMap[k].push(s);
    });

    const groupList = Object.keys(groupMap).map(k => {
      const smp = groupMap[k];
      const vals = smp.map(s => {
        let v = metric === 'tmb' ? s.tmb : (metric === 'mut_count' ? s.mut_count : s.coding_mut);
        return v || 0.01;
      }).sort((a, b) => a - b);
      const med = vals[Math.floor(vals.length / 2)] || 0;
      return { key: k, samples: smp, median: med };
    });

    groupList.sort((a, b) => b.median - a.median);

    groupList.forEach((g, idx) => {
      const shortK = g.key.replace(/ Cancer| Tumor| Carcinoma/g, '');
      columns.push({
        id: g.key,
        label: `${shortK}`,
        cohort: g.key,
        subLabel: 'All',
        samples: g.samples,
        color: '#0284c7',
        isMt: false
      });
    });
  }

  if (columns.length === 0) {
    document.getElementById('tmb-svg').innerHTML = '<text x="500" y="260" text-anchor="middle" fill="#64748b">No samples pass the current filters.</text>';
    return;
  }

  // Calculate values per sample and bounds
  const getRawVal = s => {
    if (metric === 'tmb') return s.tmb;
    if (metric === 'mut_count') return s.mut_count;
    return s.coding_mut;
  };

  const getYCoordVal = v => {
    if (scale === 'log10') {
      if (metric === 'tmb') {
        if (!v || v <= 0.0316) return -1.5;
        return Math.log10(v);
      } else {
        if (!v || v <= 0) return 0;
        return Math.log10(Math.max(1, v));
      }
    }
    const val = Math.max(0.001, v);
    if (scale === 'log2') return Math.log2(val);
    return val;
  };

  let minValY = Infinity, maxValY = -Infinity;
  columns.forEach(col => {
    col.samples.forEach(s => {
      const raw = getRawVal(s);
      const yv = getYCoordVal(raw);
      s._calcY = yv;
      s._rawVal = raw;
      if (yv < minValY) minValY = yv;
      if (yv > maxValY) maxValY = yv;
    });
  });

  const gridTicks = [];

  if (scale === 'log2') {
    minValY = Math.min(-5.0, Math.floor(minValY));
    maxValY = Math.max(7.0, Math.ceil(maxValY));
    for (let t = Math.ceil(minValY); t <= Math.floor(maxValY); t += (maxValY - minValY > 12 ? 4 : (maxValY - minValY > 8 ? 2 : 1))) {
      gridTicks.push({ val: t, label: String(t) });
    }
    if (!gridTicks.some(t => t.val === -5)) gridTicks.push({ val: -5, label: '-5' });
    if (!gridTicks.some(t => t.val === 0)) gridTicks.push({ val: 0, label: '0' });
    if (!gridTicks.some(t => t.val === 5)) gridTicks.push({ val: 5, label: '5' });
    gridTicks.sort((a, b) => a.val - b.val);
  } else if (scale === 'log10') {
    if (metric === 'tmb') {
      minValY = -1.5;
      maxValY = Math.max(2.6, Math.ceil(maxValY * 10) / 10);
      gridTicks.push({ val: -1.5, label: '0' });
      gridTicks.push({ val: -1.0, label: '0.1' });
      gridTicks.push({ val: 0.0, label: '1' });
      gridTicks.push({ val: 1.0, label: '10' });
      gridTicks.push({ val: 2.0, label: '100' });
      if (maxValY >= 2.8) {
        gridTicks.push({ val: 3.0, label: '1,000' });
      }
    } else if (metric === 'mut_count') {
      minValY = Math.min(1.0, Math.floor(minValY));
      maxValY = Math.max(5.5, Math.ceil(maxValY));
      gridTicks.push({ val: 1.0, label: '10' });
      gridTicks.push({ val: 2.0, label: '100' });
      gridTicks.push({ val: 3.0, label: '1,000' });
      gridTicks.push({ val: 4.0, label: '10,000' });
      gridTicks.push({ val: 5.0, label: '100,000' });
      if (maxValY >= 5.8) gridTicks.push({ val: 6.0, label: '1,000,000' });
    } else { // coding_mut
      minValY = Math.min(0, Math.floor(minValY));
      maxValY = Math.max(3.5, Math.ceil(maxValY));
      gridTicks.push({ val: 0.0, label: '1' });
      gridTicks.push({ val: 1.0, label: '10' });
      gridTicks.push({ val: 2.0, label: '100' });
      gridTicks.push({ val: 3.0, label: '1,000' });
      if (maxValY >= 3.8) gridTicks.push({ val: 4.0, label: '10,000' });
    }
  } else {
    minValY = 0;
    maxValY = maxValY * 1.05;
    const step = maxValY / 5;
    for (let i = 0; i <= 5; i++) {
      const v = Math.round(i * step);
      gridTicks.push({ val: v, label: v.toLocaleString() });
    }
  }

  // Dimensions
  const M = columns.length;
  const colWidth = Math.max(26, Math.min(54, Math.floor(950 / M)));
  const padLeft = 70;
  const padRight = 30;
  const padTop = 60;
  const padBottom = 215;
  const plotWidth = M * colWidth;
  const totalW = padLeft + plotWidth + padRight;
  const totalH = 635;
  const plotH = totalH - padTop - padBottom;

  const scaleY = y => padTop + (1.0 - (y - minValY) / (maxValY - minValY)) * plotH;

  const svg = document.getElementById('tmb-svg');
  svg.setAttribute('width', totalW);
  svg.setAttribute('height', totalH);
  svg.setAttribute('viewBox', `0 0 ${totalW} ${totalH}`);

  let out = '';

  // Background
  out += `<rect width="${totalW}" height="${totalH}" fill="#ffffff" />`;

  // Draw Publication Legend inside SVG (top margin, y = 24)
  const legendItems = [
    {
      render: () => `<line x1="0" y1="0" x2="22" y2="0" stroke="#dc2626" stroke-width="2.6" stroke-linecap="round" />`,
      text: 'Cohort Median',
      color: '#334155',
      w: 125
    },
    {
      render: () => `<circle cx="6" cy="0" r="3.2" fill="#0284c7" />`,
      text: groupBy === 'pole' ? 'Biopsy (POLE WT)' : (groupBy === 'msi' ? 'Biopsy (MSS)' : 'Individual Biopsy'),
      color: '#334155',
      w: groupBy === 'pole' || groupBy === 'msi' ? 145 : 135
    }
  ];

  if (groupBy === 'pole') {
    legendItems.push({
      render: () => `<circle cx="6" cy="0" r="3.2" fill="#d97706" />`,
      text: 'Biopsy (POLE MT)',
      color: '#d97706',
      w: 145
    });
  } else if (groupBy === 'msi') {
    legendItems.push({
      render: () => `<circle cx="6" cy="0" r="3.2" fill="#dc2626" />`,
      text: 'Biopsy (MSI-H)',
      color: '#dc2626',
      w: 135
    });
  }

  if (!excludeFloor && floorVal > 0) {
    legendItems.push({
      render: () => `<circle cx="6" cy="0" r="2.8" fill="#0284c7" opacity="0.4" /><circle cx="6" cy="0" r="4.8" fill="none" stroke="#dc2626" stroke-width="1.2" opacity="0.85" />`,
      text: `Purity < ${Math.round(floorVal * 100)}% Floor`,
      color: '#64748b',
      w: 155
    });
  }

  if (showThreshold) {
    let benchText = 'FDA TMB-High (10 mut/Mb)';
    if (metric === 'mut_count') {
      const mutEl = document.getElementById('tmb-mut-threshold');
      const mv = mutEl ? (parseInt(mutEl.value) || 4000) : 4000;
      benchText = `Cutoff (≥ ${mv.toLocaleString()} muts)`;
    } else if (metric === 'coding_mut') {
      const mutEl = document.getElementById('tmb-mut-threshold');
      const mv = mutEl ? (parseInt(mutEl.value) || 4000) : 4000;
      benchText = `Cutoff (≥ ${Math.round(mv * (32 / 2850))} coding)`;
    }
    legendItems.push({
      render: () => `<line x1="0" y1="0" x2="22" y2="0" stroke="#10b981" stroke-width="1.8" stroke-dasharray="5,3" />`,
      text: benchText,
      color: '#059669',
      w: 180
    });
  }

  let legX = padLeft;
  const legY = 24;
  legendItems.forEach(item => {
    out += `<g transform="translate(${legX}, ${legY})">`;
    out += item.render();
    out += `<text x="28" y="4" font-size="11" font-weight="600" fill="${item.color}" font-family="system-ui, -apple-system, sans-serif">${item.text}</text>`;
    out += `</g>`;
    legX += item.w;
  });

  // Populate Interactive HTML Legend Bar
  const legendBar = document.getElementById('tmb-legend-bar');
  if (legendBar) {
    let barHtml = `
      <div style="font-size:11px; text-transform:uppercase; letter-spacing:0.05em; color:#64748b; font-weight:700; margin-right:4px;">Legend:</div>
      <div style="display:flex; align-items:center; gap:6px;">
        <span style="display:inline-block; width:18px; height:3px; background:#dc2626; border-radius:2px;"></span>
        <span>Cohort Median</span>
      </div>
      <div style="display:flex; align-items:center; gap:6px;">
        <span style="display:inline-block; width:8px; height:8px; background:#0284c7; border-radius:50%;"></span>
        <span>${groupBy === 'pole' ? 'Biopsy (POLE WT)' : (groupBy === 'msi' ? 'Biopsy (MSS)' : 'Individual Biopsy')}</span>
      </div>
    `;
    if (groupBy === 'pole') {
      barHtml += `
        <div style="display:flex; align-items:center; gap:6px;">
          <span style="display:inline-block; width:8px; height:8px; background:#d97706; border-radius:50%;"></span>
          <span style="color:#b45309;">Biopsy (POLE MT Hypermutated)</span>
        </div>
      `;
    } else if (groupBy === 'msi') {
      barHtml += `
        <div style="display:flex; align-items:center; gap:6px;">
          <span style="display:inline-block; width:8px; height:8px; background:#dc2626; border-radius:50%;"></span>
          <span style="color:#b91c1c;">Biopsy (MSI-H)</span>
        </div>
      `;
    }
    if (!excludeFloor && floorVal > 0) {
      barHtml += `
        <div style="display:flex; align-items:center; gap:6px;">
          <span style="display:inline-block; width:10px; height:10px; border:1.5px solid #dc2626; border-radius:50%; background:rgba(2,132,199,0.25);"></span>
          <span style="color:#64748b;">Purity &lt; ${Math.round(floorVal * 100)}% (Floor Clamped)</span>
        </div>
      `;
    }
    if (showThreshold) {
      let bLabel = 'FDA TMB-High (≥ 10 mut/Mb)';
      if (metric === 'mut_count') {
        const mutEl = document.getElementById('tmb-mut-threshold');
        const mv = mutEl ? (parseInt(mutEl.value) || 4000) : 4000;
        bLabel = `Cutoff Benchmark (≥ ${mv.toLocaleString()} mutations)`;
      } else if (metric === 'coding_mut') {
        const mutEl = document.getElementById('tmb-mut-threshold');
        const mv = mutEl ? (parseInt(mutEl.value) || 4000) : 4000;
        bLabel = `Cutoff Benchmark (≥ ${Math.round(mv * (32 / 2850))} coding mutations)`;
      }
      barHtml += `
        <div style="display:flex; align-items:center; gap:6px;">
          <span style="display:inline-block; width:18px; border-top:2px dashed #10b981;"></span>
          <span style="color:#059669;">${bLabel}</span>
        </div>
      `;
    }
    barHtml += `
      <div style="display:flex; align-items:center; gap:6px; margin-left:auto; color:#64748b; font-size:11px;">
        <span style="display:inline-block; width:14px; height:10px; background:#f8fafc; border:1px solid #cbd5e1; border-radius:2px;"></span>
        <span>Alternating Cohort Columns</span>
      </div>
    `;
    legendBar.innerHTML = barHtml;
  }

  // Draw Alternating Background Stripes per cohort / column
  let lastCohort = null;
  let cohortShade = false;

  columns.forEach((col, k) => {
    const x0 = padLeft + k * colWidth;
    let bg = '#ffffff';
    if (groupBy === 'pole' || groupBy === 'msi') {
      if (col.cohort !== lastCohort) {
        lastCohort = col.cohort;
        cohortShade = !cohortShade;
      }
      bg = cohortShade ? '#f1f5f9' : '#ffffff';
    } else {
      bg = (k % 2 === 0) ? '#f8fafc' : '#ffffff';
    }

    out += `<rect x="${x0}" y="${padTop}" width="${colWidth}" height="${plotH}" fill="${bg}" />`;
  });

  gridTicks.forEach(tick => {
    const py = scaleY(tick.val);
    if (py < padTop - 2 || py > padTop + plotH + 2) return;
    out += `<line x1="${padLeft}" y1="${py}" x2="${padLeft + plotWidth}" y2="${py}" stroke="#cbd5e1" stroke-width="1" stroke-dasharray="4,4" />`;
    out += `<text x="${padLeft - 10}" y="${py + 4}" text-anchor="end" font-size="11" font-weight="500" fill="#64748b" font-family="system-ui, sans-serif">${tick.label}</text>`;
  });

  // Benchmark or Cutoff line if requested
  if (showThreshold) {
    let threshVal = 10.0;
    let threshLabel = 'FDA TMB-High (10 mut/Mb)';
    if (metric === 'mut_count') {
      const mutEl = document.getElementById('tmb-mut-threshold');
      threshVal = mutEl ? (parseInt(mutEl.value) || 4000) : 4000;
      threshLabel = `Cutoff (≥ ${threshVal.toLocaleString()} muts)`;
    } else if (metric === 'coding_mut') {
      const mutEl = document.getElementById('tmb-mut-threshold');
      const mc = mutEl ? (parseInt(mutEl.value) || 4000) : 4000;
      threshVal = Math.round(mc * (32 / 2850));
      threshLabel = `Cutoff (≥ ${threshVal} coding)`;
    }
    const ty = scaleY(getYCoordVal(threshVal));
    if (ty >= padTop && ty <= padTop + plotH) {
      out += `<line x1="${padLeft}" y1="${ty}" x2="${padLeft + plotWidth}" y2="${ty}" stroke="#10b981" stroke-width="1.6" stroke-dasharray="6,3" />`;
      out += `<text x="${padLeft + plotWidth - 6}" y="${ty - 5}" text-anchor="end" font-size="10.5" font-weight="700" fill="#059669">${threshLabel}</text>`;
    }
  }

  // Y-axis main vertical line and bottom baseline
  out += `<line x1="${padLeft}" y1="${padTop}" x2="${padLeft}" y2="${padTop + plotH}" stroke="#94a3b8" stroke-width="1.2" />`;
  out += `<line x1="${padLeft}" y1="${padTop + plotH}" x2="${padLeft + plotWidth}" y2="${padTop + plotH}" stroke="#94a3b8" stroke-width="1.2" />`;
  out += `<text x="${-(padTop + plotH / 2)}" y="24" transform="rotate(-90)" text-anchor="middle" font-size="12" font-weight="600" fill="#334155" font-family="system-ui, sans-serif">${yAxisTitle}</text>`;

  // Summary table rows collector
  const tableRowsData = [];

  // Plot Data in Each Column
  columns.forEach((col, k) => {
    const x0 = padLeft + k * colWidth;
    const x1 = x0 + colWidth;
    const n = col.samples.length;

    col.samples.sort((a, b) => a._calcY - b._calcY);

    const rawVals = col.samples.map(s => s._rawVal).sort((a, b) => a - b);
    const medVal = rawVals[Math.floor(n / 2)] || 0;
    const q1 = rawVals[Math.floor(n * 0.25)] || 0;
    const q3 = rawVals[Math.floor(n * 0.75)] || 0;
    const minV = rawVals[0] || 0;
    const maxV = rawVals[n - 1] || 0;
    const tmbHightCount = col.samples.filter(s => s.tmb >= 10.0).length;
    const tmbHighPct = n > 0 ? (tmbHightCount / n * 100).toFixed(1) : '0.0';

    const purities = col.samples.map(s => s.purity).filter(p => p !== null && p !== undefined).sort((a, b) => a - b);
    const medPurity = purities.length > 0 ? purities[Math.floor(purities.length / 2)] : null;

    tableRowsData.push({
      cohort: col.cohort,
      sub: col.subLabel,
      n: n,
      median: medVal,
      q1: q1,
      q3: q3,
      min: minV,
      max: maxV,
      highPct: tmbHighPct,
      purity: medPurity !== null ? `${Math.round(medPurity * 100)}%` : 'NA'
    });

    const medYPixel = scaleY(getYCoordVal(medVal));

    // Draw Points
    const padSide = Math.max(2, colWidth * 0.08);
    const dotW = colWidth - 2 * padSide;

    col.samples.forEach((s, idx) => {
      const cx = (x0 + padSide + ((idx + 0.5) / n) * dotW).toFixed(1);
      const cy = scaleY(s._calcY).toFixed(1);

      const isBelowFloor = s.purity !== null && s.purity !== undefined && s.purity < floorVal;
      const opacity = isBelowFloor ? 0.35 : 0.85;
      const ptColor = col.color;

      out += `<circle cx="${cx}" cy="${cy}" r="2.7" fill="${ptColor}" opacity="${opacity}" style="cursor:pointer;"
                     onmouseover="showTmbTip(event, '${s.id}')" onmouseout="hideTmbTip()"></circle>`;

      if (isBelowFloor && !excludeFloor) {
        out += `<circle cx="${cx}" cy="${cy}" r="4.2" fill="none" stroke="#dc2626" stroke-width="0.8" opacity="0.6" pointer-events="none" />`;
      }
    });

    // Draw Red Horizontal Median Bar across column width
    out += `<line x1="${x0 + 3}" y1="${medYPixel}" x2="${x1 - 3}" y2="${medYPixel}" stroke="#dc2626" stroke-width="2.6" stroke-linecap="round" />`;

    // Column Label Rotated at the bottom
    const lx = (x0 + x1) / 2;
    const ly = padTop + plotH + 8;
    // Tick mark on baseline
    out += `<line x1="${lx}" y1="${padTop + plotH}" x2="${lx}" y2="${padTop + plotH + 4}" stroke="#cbd5e1" stroke-width="1" />`;
    // Clean, crisp rotated label with generous bottom clearance
    out += `<text x="${lx}" y="${ly}" transform="rotate(-90 ${lx} ${ly})" text-anchor="end" dominant-baseline="central" font-size="10.5" font-weight="600" fill="#334155" font-family="system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif">${col.label}</text>`;
  });

  svg.innerHTML = out;

  // Populate Summary Table
  const tbody = document.getElementById('tmb-summary-body');
  if (tbody) {
    let trows = '';
    tableRowsData.forEach(row => {
      trows += `
        <tr>
          <td style="font-weight:600; color:#0f172a;">${row.cohort}</td>
          <td><span class="badge ${row.sub.includes('MT') || row.sub.includes('MSI') ? 'badge-primary' : 'badge-neutral'}" style="font-size:10px;">${row.sub}</span></td>
          <td class="num">${row.n}</td>
          <td class="num" style="font-weight:700; color:#0284c7;">${row.median.toLocaleString(undefined, { maximumFractionDigits: 2 })}</td>
          <td class="num">${row.q1.toFixed(1)} &ndash; ${row.q3.toFixed(1)}</td>
          <td class="num">${row.min.toFixed(1)}</td>
          <td class="num">${row.max.toFixed(1)}</td>
          <td class="num" style="font-weight:${parseFloat(row.highPct) > 20 ? '700' : '500'}; color:${parseFloat(row.highPct) > 20 ? '#059669' : 'inherit'};">${row.highPct}%</td>
          <td class="num">${row.purity}</td>
        </tr>
      `;
    });
    tbody.innerHTML = trows;
    makeTableSortable('tmb-summary-table');
  }

  // Also render prevalence of cases above mutation threshold
  renderTmbThresholdPlot();
}

function showTmbTip(e, id) {
  const tip = document.getElementById('tmb-tooltip');
  if (!tip || !window.cohortTmbData || !window.cohortTmbData.sampleMap) return;
  const s = window.cohortTmbData.sampleMap[id];
  if (!s) return;

  const floorVal = parseFloat(document.getElementById('tmb-purity-floor').value) || 0.15;
  const purityText = s.purity !== null && s.purity !== undefined
    ? `${(s.purity * 100).toFixed(0)}% ${s.purity < floorVal ? '<span style="color:#f87171; font-weight:700;">(&lt; ' + Math.round(floorVal * 100) + '% floor)</span>' : ''}`
    : '<span style="color:#94a3b8;">Uncalled</span>';

  tip.innerHTML = `
    <div style="font-weight:700; color:#38bdf8; font-size:12px; margin-bottom:3px; font-family:ui-monospace, monospace;">${s.id}</div>
    <div>Cohort: <strong>${s.cancer}</strong></div>
    <div>CUPPA: <strong>${s.cuppa}</strong></div>
    <div>Purity: <strong>${purityText}</strong></div>
    <div style="border-top:1px solid #334155; margin-top:5px; padding-top:4px;">
      <div>TMB: <strong style="color:#34d399;">${s.tmb} mut/Mb</strong> ${s.tmb >= 10 ? '<span style="color:#10b981; font-weight:700;">[TMB-High]</span>' : ''}</div>
      <div>WGS Mutations: <strong>${s.mut_count.toLocaleString()}</strong></div>
      <div>Coding Load: <strong>${s.coding_mut}</strong></div>
      <div>MSI Status: <strong>${s.msi}</strong></div>
      <div>POLE Sig 10: <strong>${s.pole_sig}</strong> ${s.is_pole_mt ? '<span style="color:#fbbf24; font-weight:700;">[MT]</span>' : '<span style="color:#94a3b8;">[WT]</span>'}</div>
    </div>
  `;

  const container = document.getElementById('tmb-plot-container');
  const rect = container.getBoundingClientRect();
  const x = e.clientX - rect.left + 14;
  const y = e.clientY - rect.top + 14;

  tip.style.left = Math.min(x, rect.width - 290) + 'px';
  tip.style.top = Math.min(y, rect.height - 210) + 'px';
  tip.style.display = 'block';
}

function hideTmbTip() {
  const tip = document.getElementById('tmb-tooltip');
  if (tip) tip.style.display = 'none';
}

function renderTmbThresholdPlot() {
  const data = window.cohortTmbData;
  if (!data || !data.samples || data.samples.length === 0) return;

  const floorVal = parseFloat(document.getElementById('tmb-purity-floor').value) || 0.0;
  const excludeFloor = document.getElementById('tmb-purity-exclude').checked;
  const threshInput = document.getElementById('tmb-mut-threshold');
  const cutoff = parseInt(threshInput ? threshInput.value : '4000') || 4000;

  // Filter samples respecting purity floor
  const samples = data.samples.filter(s => {
    if (s.mut_count === null || s.mut_count === undefined) return false;
    const isBelowFloor = s.purity !== null && s.purity !== undefined && s.purity < floorVal;
    return !(isBelowFloor && excludeFloor);
  });

  const totalCohort = samples.length;
  const totalAbove = samples.filter(s => s.mut_count >= cutoff).length;
  const totalPct = totalCohort > 0 ? (totalAbove / totalCohort * 100).toFixed(1) : '0.0';

  const badge = document.getElementById('tmb-thresh-cohort-badge');
  if (badge) {
    badge.textContent = `Overall: ${totalPct}% (${totalAbove.toLocaleString()} / ${totalCohort.toLocaleString()} cases ≥ ${cutoff.toLocaleString()} muts)`;
  }

  // Aggregate by cancer type
  const cancerMap = {};
  samples.forEach(s => {
    const c = s.cancer || 'Unknown';
    if (!cancerMap[c]) {
      cancerMap[c] = {
        cohort: c,
        total: 0,
        above: 0,
        msiHigh: 0,
        poleMt: 0,
        tmbVals: []
      };
    }
    cancerMap[c].total++;
    if (s.tmb !== null && s.tmb !== undefined) cancerMap[c].tmbVals.push(s.tmb);
    if (s.mut_count >= cutoff) {
      cancerMap[c].above++;
      if (s.msi === 'MSI') cancerMap[c].msiHigh++;
      if (s.is_pole_mt) cancerMap[c].poleMt++;
    }
  });

  const cancerList = Object.values(cancerMap).map(g => {
    const pct = g.total > 0 ? (g.above / g.total * 100) : 0.0;
    g.tmbVals.sort((a, b) => a - b);
    const medTmb = g.tmbVals.length > 0 ? g.tmbVals[Math.floor(g.tmbVals.length / 2)] : 0.0;
    return {
      cohort: g.cohort,
      total: g.total,
      above: g.above,
      pct: pct,
      medTmb: medTmb,
      msiHigh: g.msiHigh,
      poleMt: g.poleMt
    };
  });

  // Sort by % descending, then by count descending
  cancerList.sort((a, b) => b.pct - a.pct || b.above - a.above);

  // SVG Dimensions
  const svg = document.getElementById('tmb-thresh-svg');
  if (!svg) return;

  const N = cancerList.length;
  const padLeft = 210;
  const padRight = 140;
  const padTop = 35;
  const padBottom = 40;
  const barWidth = 600;
  const rowHeight = 23;
  const barHeight = 15;
  const plotH = N * rowHeight;
  const totalW = padLeft + barWidth + padRight;
  const totalH = padTop + plotH + padBottom;

  svg.setAttribute('width', totalW);
  svg.setAttribute('height', totalH);
  svg.setAttribute('viewBox', `0 0 ${totalW} ${totalH}`);

  let out = '';
  out += `<rect width="${totalW}" height="${totalH}" fill="#ffffff" />`;

  // Grid lines at 0%, 20%, 40%, 60%, 80%, 100%
  const gridPcts = [0, 20, 40, 60, 80, 100];
  gridPcts.forEach(p => {
    const gx = padLeft + (p / 100) * barWidth;
    out += `<line x1="${gx}" y1="${padTop}" x2="${gx}" y2="${padTop + plotH}" stroke="${p === 0 ? '#94a3b8' : '#e2e8f0'}" stroke-width="${p === 0 ? 1.4 : 1}" stroke-dasharray="${p === 0 ? 'none' : '4,4'}" />`;
    out += `<text x="${gx}" y="${padTop - 8}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#64748b" font-family="system-ui, sans-serif">${p}%</text>`;
    out += `<text x="${gx}" y="${padTop + plotH + 18}" text-anchor="middle" font-size="10.5" font-weight="600" fill="#64748b" font-family="system-ui, sans-serif">${p}%</text>`;
  });

  // Rows
  cancerList.forEach((c, idx) => {
    const y0 = padTop + idx * rowHeight;
    const yBar = y0 + (rowHeight - barHeight) / 2;
    const bw = Math.max(0, (c.pct / 100) * barWidth);

    // Alternating hoverable row background
    const bgRow = (idx % 2 === 0) ? '#f8fafc' : '#ffffff';
    out += `<rect x="0" y="${y0}" width="${totalW}" height="${rowHeight}" fill="${bgRow}" opacity="0.6" />`;

    // Cohort Name on the left
    const shortC = c.cohort.replace(/ Cancer| Tumor| Carcinoma/g, '');
    out += `<text x="${padLeft - 12}" y="${yBar + 11}" text-anchor="end" font-size="11" font-weight="600" fill="#1e293b" font-family="system-ui, sans-serif">${shortC}</text>`;

    // Bar background track
    out += `<rect x="${padLeft}" y="${yBar}" width="${barWidth}" height="${barHeight}" fill="#f1f5f9" rx="3" />`;

    // Active Bar
    let barColor = '#0284c7';
    if (c.pct >= 50.0) barColor = '#dc2626';
    else if (c.pct >= 25.0) barColor = '#d97706';
    else if (c.pct > 0.0) barColor = '#0284c7';

    if (bw > 0) {
      out += `<rect x="${padLeft}" y="${yBar}" width="${bw}" height="${barHeight}" fill="${barColor}" rx="3" style="cursor:pointer;"
                    onmouseover="showThreshTip(event, '${c.cohort.replace(/'/g, "\\'")}', ${c.above}, ${c.total}, ${c.pct.toFixed(1)}, ${c.medTmb.toFixed(2)}, ${c.msiHigh}, ${c.poleMt}, ${cutoff})"
                    onmouseout="hideThreshTip()"></rect>`;
    }

    // Label: XX.X% (n / N)
    const labelX = padLeft + bw + 8;
    const labelColor = c.above > 0 ? '#0f172a' : '#94a3b8';
    const boldWeight = c.pct >= 20.0 ? '700' : '600';
    out += `<text x="${labelX}" y="${yBar + 11}" font-size="11" font-weight="${boldWeight}" fill="${labelColor}" font-family="ui-monospace, monospace">${c.pct.toFixed(1)}% <tspan fill="#64748b" font-weight="500">(${c.above}/${c.total})</tspan></text>`;
  });

  svg.innerHTML = out;
}

function showThreshTip(e, cohort, above, total, pct, medTmb, msiHigh, poleMt, cutoff) {
  const tip = document.getElementById('tmb-thresh-tooltip');
  if (!tip) return;
  tip.innerHTML = `
    <div style="font-weight:700; color:#38bdf8; font-size:12px; margin-bottom:4px;">${cohort}</div>
    <div>Prevalence: <strong style="color:#34d399; font-size:12px;">${pct}%</strong> (${above} / ${total} cases)</div>
    <div>Cutoff: <strong>&ge; ${cutoff.toLocaleString()} mutations</strong> (~ ${(cutoff / 2850).toFixed(1)} mut/Mb)</div>
    <div style="border-top:1px solid #334155; margin-top:5px; padding-top:4px;">
      <div>Cohort Median TMB: <strong>${medTmb} mut/Mb</strong></div>
      <div>MSI-High in subgroup: <strong>${msiHigh}</strong></div>
      <div>POLE MT in subgroup: <strong>${poleMt}</strong></div>
    </div>
  `;
  const container = document.getElementById('tmb-thresh-container');
  const rect = container.getBoundingClientRect();
  const x = e.clientX - rect.left + 14;
  const y = e.clientY - rect.top + 14;
  tip.style.left = Math.min(x, rect.width - 290) + 'px';
  tip.style.top = Math.min(y, rect.height - 180) + 'px';
  tip.style.display = 'block';
}

function hideThreshTip() {
  const tip = document.getElementById('tmb-thresh-tooltip');
  if (tip) tip.style.display = 'none';
}

function setTmbThresholdPreset(val) {
  const el = document.getElementById('tmb-mut-threshold');
  if (el) {
    el.value = val;
    renderTmbThresholdPlot();
    renderTmbPlot();
  }
}

function exportTmbThresholdSvg() {
  const svg = document.getElementById('tmb-thresh-svg');
  if (!svg) return;
  const xml = new XMLSerializer().serializeToString(svg);
  const blob = new Blob([xml], { type: 'image/svg+xml' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  const cutoff = document.getElementById('tmb-mut-threshold') ? document.getElementById('tmb-mut-threshold').value : '4000';
  a.download = `cohort_tumour_types_above_${cutoff}_mutations_${new Date().toISOString().split('T')[0]}.svg`;
  a.click();
  URL.revokeObjectURL(a.href);
}

function exportTmbThresholdPng() {
  const cutoff = document.getElementById('tmb-mut-threshold') ? document.getElementById('tmb-mut-threshold').value : '4000';
  exportSvgToPng('#tmb-thresh-svg', `cohort_tumour_types_above_${cutoff}_mutations_${new Date().toISOString().split('T')[0]}.png`);
}

function exportTmbPlotSvg() {
  const svg = document.getElementById('tmb-svg');
  if (!svg) return;
  const xml = new XMLSerializer().serializeToString(svg);
  const blob = new Blob([xml], { type: 'image/svg+xml' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `cohort_tmb_mutation_burden_${new Date().toISOString().split('T')[0]}.svg`;
  a.click();
  URL.revokeObjectURL(a.href);
}

function exportTmbPlotPng() {
  exportSvgToPng('#tmb-svg', `cohort_tmb_mutation_burden_${new Date().toISOString().split('T')[0]}.png`);
}

function exportTmbStandaloneHtml() {
  const svgEl = document.getElementById('tmb-svg');
  if (!svgEl) return;
  const title = document.getElementById('tmb-card-title').textContent;
  const metric = document.getElementById('tmb-metric').value;
  const legendEl = document.getElementById('tmb-legend-bar');
  const legendHtml = legendEl ? legendEl.outerHTML : '';
  const svgHtml = svgEl.outerHTML;

  const html = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Mutation Burden Plot: ${title}</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f8fafc; color: #0f172a; margin: 0; padding: 2rem; }
  .container { max-width: 1200px; margin: 0 auto; background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 2rem; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05); }
  h1 { font-size: 1.4rem; font-weight: 700; margin-bottom: 0.4rem; color: #0f172a; }
  .subtitle { font-size: 0.9rem; color: #475569; margin-bottom: 1.5rem; }
  .card { border: 1px solid #e2e8f0; border-radius: 6px; padding: 1.25rem; background: #ffffff; overflow-x: auto; }
  svg { display: block; }
  .footer-note { font-size: 11px; color: #64748b; margin-top: 1rem; }
</style>
</head>
<body>
<div class="container">
  <h1>${title}</h1>
  <div class="subtitle">Metric: ${metric}</div>
  <div class="card">
    ${legendHtml}
    ${svgHtml}
  </div>
</div>
</body>
</html>`;

  const blob = new Blob([html], { type: 'text/html' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `cohort_tmb_plot_${new Date().toISOString().split('T')[0]}.html`;
  a.click();
  URL.revokeObjectURL(a.href);
}

function exportTmbTableCsv() {
  const table = document.getElementById('tmb-summary-table');
  if (!table) return;
  const rows = Array.from(table.querySelectorAll('tr'));
  const csv = rows.map(r => {
    const cols = Array.from(r.querySelectorAll('th, td'));
    return cols.map(c => `"${c.textContent.replace(/"/g, '""').trim()}"`).join(',');
  }).join('\\n');

  const blob = new Blob([csv], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `cohort_tmb_summary_${new Date().toISOString().split('T')[0]}.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

function initAllSortableTables() {
  document.querySelectorAll('table.sortable').forEach(table => {
    makeTableSortable(table);
  });
}

function makeTableSortable(tableOrId) {
  const table = typeof tableOrId === 'string' ? document.getElementById(tableOrId) : tableOrId;
  if (!table) return;

  const thead = table.querySelector('thead');
  if (!thead) return;

  const ths = thead.querySelectorAll('th.sortable');
  ths.forEach((th, colIdx) => {
    let indicator = th.querySelector('.sort-indicator');
    if (!indicator) {
      indicator = document.createElement('span');
      indicator.className = 'sort-indicator';
      indicator.textContent = ' ↕';
      th.appendChild(indicator);
    }

    th.onclick = function() {
      const isAsc = th.dataset.sortAsc === 'true';
      const newAsc = !isAsc;

      ths.forEach(h => {
        delete h.dataset.sortAsc;
        const ind = h.querySelector('.sort-indicator');
        if (ind) {
          ind.textContent = ' ↕';
          ind.style.color = 'var(--text-muted)';
        }
      });

      th.dataset.sortAsc = newAsc ? 'true' : 'false';
      indicator.textContent = newAsc ? ' ▲' : ' ▼';
      indicator.style.color = 'var(--text-primary)';

      sortTableByCol(table, colIdx, newAsc);
    };
  });
}

function sortTableByCol(table, colIdx, asc) {
  const tbody = table.querySelector('tbody');
  if (!tbody) return;
  const rows = Array.from(tbody.querySelectorAll('tr'));
  if (rows.length <= 1) return;

  if (rows[0].children.length === 1 && rows[0].children[0].hasAttribute('colspan')) return;

  rows.sort((rowA, rowB) => {
    const cellA = rowA.children[colIdx] ? rowA.children[colIdx].textContent.trim() : '';
    const cellB = rowB.children[colIdx] ? rowB.children[colIdx].textContent.trim() : '';

    const cleanA = cellA.replace(/[%$,+*]/g, '').trim();
    const cleanB = cellB.replace(/[%$,+*]/g, '').trim();
    const numA = parseFloat(cleanA);
    const numB = parseFloat(cleanB);

    const isNumA = !isNaN(numA) && cleanA !== '' && isFinite(cleanA);
    const isNumB = !isNaN(numB) && cleanB !== '' && isFinite(cleanB);

    if (isNumA && isNumB) {
      return asc ? (numA - numB) : (numB - numA);
    }
    return asc ? cellA.localeCompare(cellB) : cellB.localeCompare(cellA);
  });

  rows.forEach(r => tbody.appendChild(r));
}

// ---------------------------------------------------------------------------
// DATA & PLOT EXPORT UTILITIES (STANDALONE HTML, SVG, CSV)
// ---------------------------------------------------------------------------
function exportSvg(selector, filename) {
  const svg = document.querySelector(selector);
  if (!svg) {
    alert("No SVG plot available to export.");
    return;
  }
  const serializer = new XMLSerializer();
  let source = serializer.serializeToString(svg);
  if (!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)) {
    source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
  }
  const blob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename || 'plot.svg';
  a.click();
  URL.revokeObjectURL(url);
}

function exportSvgToPng(selectorOrEl, filename, scale = 2) {
  const svg = typeof selectorOrEl === 'string' ? document.querySelector(selectorOrEl) : selectorOrEl;
  if (!svg) {
    alert("Plot not found to export.");
    return;
  }

  const serializer = new XMLSerializer();
  let source = serializer.serializeToString(svg);
  if (!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)) {
    source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
  }

  let w = (svg.viewBox && svg.viewBox.baseVal && svg.viewBox.baseVal.width) || (svg.width && svg.width.baseVal && svg.width.baseVal.value) || svg.clientWidth || 1000;
  let h = (svg.viewBox && svg.viewBox.baseVal && svg.viewBox.baseVal.height) || (svg.height && svg.height.baseVal && svg.height.baseVal.value) || svg.clientHeight || 600;

  const svgBlob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
  const URLObj = window.URL || window.webkitURL || window;
  const url = URLObj.createObjectURL(svgBlob);

  const img = new Image();
  img.onload = function() {
    const canvas = document.createElement('canvas');
    canvas.width = Math.round(w * scale);
    canvas.height = Math.round(h * scale);
    const ctx = canvas.getContext('2d');
    ctx.scale(scale, scale);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0, w, h);

    canvas.toBlob(function(blob) {
      const a = document.createElement('a');
      const pngUrl = blob ? URLObj.createObjectURL(blob) : canvas.toDataURL('image/png');
      a.href = pngUrl;
      a.download = filename || 'plot.png';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      if (blob) URLObj.revokeObjectURL(pngUrl);
      URLObj.revokeObjectURL(url);
    }, 'image/png');
  };
  img.onerror = function() {
    URLObj.revokeObjectURL(url);
    alert("Could not render SVG to PNG.");
  };
  img.src = url;
}

function exportTableToPNG(tableIdOrEl, filename, scale = 2) {
  const table = typeof tableIdOrEl === 'string' ? document.getElementById(tableIdOrEl) : tableIdOrEl;
  if (!table) {
    alert("Table not found to export.");
    return;
  }

  const thead = table.querySelector('thead');
  const tbody = table.querySelector('tbody');
  if (!thead || !tbody) return;

  const headRows = Array.from(thead.querySelectorAll('tr'));
  const bodyRows = Array.from(tbody.querySelectorAll('tr')).filter(r => r.style.display !== 'none');
  if (bodyRows.length === 0) return;

  const colCount = Math.max(...headRows.map(r => r.children.length), ...bodyRows.map(r => r.children.length));
  const paddingX = 14;
  const fontHeader = "bold 12px -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif";
  const fontBody = "12px -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif";
  const fontMono = "12px ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace";

  const tempCanvas = document.createElement('canvas');
  const tempCtx = tempCanvas.getContext('2d');

  const colWidths = new Array(colCount).fill(70);
  const colAlign = new Array(colCount).fill('left');

  headRows.forEach(tr => {
    Array.from(tr.children).forEach((th, cIdx) => {
      tempCtx.font = fontHeader;
      const text = th.innerText.replace(/[▲▼↕*]/g, '').replace(/\s+/g, ' ').trim();
      const tw = tempCtx.measureText(text).width + paddingX * 2 + 10;
      if (tw > colWidths[cIdx]) colWidths[cIdx] = tw;
      if (th.classList.contains('num') || th.style.textAlign === 'right') colAlign[cIdx] = 'right';
    });
  });

  const parsedRows = [];
  bodyRows.forEach(tr => {
    const rowData = [];
    Array.from(tr.children).forEach((td, cIdx) => {
      const isNum = td.classList.contains('num') || td.style.textAlign === 'right';
      tempCtx.font = isNum ? fontMono : fontBody;
      const text = td.innerText.replace(/\s+/g, ' ').trim();
      const tw = tempCtx.measureText(text).width + paddingX * 2 + 6;
      if (tw > colWidths[cIdx]) colWidths[cIdx] = Math.min(280, tw);

      const computed = window.getComputedStyle(td);
      const bg = td.style.backgroundColor || (computed.backgroundColor !== 'rgba(0, 0, 0, 0)' && computed.backgroundColor !== 'transparent' ? computed.backgroundColor : null);
      const fg = td.style.color || (computed.color !== 'rgba(0, 0, 0, 0)' ? computed.color : '#0f172a');
      const weight = td.style.fontWeight || computed.fontWeight || 'normal';

      rowData.push({
        text: text,
        isNum: isNum,
        bg: bg,
        fg: fg,
        weight: weight
      });
    });
    parsedRows.push(rowData);
  });

  const totalTableWidth = colWidths.reduce((a, b) => a + b, 0);
  const rowHeight = 28;
  const headerHeight = 32;
  const totalTableHeight = (headRows.length * headerHeight) + (parsedRows.length * rowHeight) + 12;

  const canvas = document.createElement('canvas');
  canvas.width = Math.round(totalTableWidth * scale);
  canvas.height = Math.round(totalTableHeight * scale);
  const ctx = canvas.getContext('2d');
  ctx.scale(scale, scale);

  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, totalTableWidth, totalTableHeight);

  let curY = 0;
  headRows.forEach(tr => {
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(0, curY, totalTableWidth, headerHeight);

    ctx.fillStyle = '#cbd5e1';
    ctx.fillRect(0, curY + headerHeight - 1, totalTableWidth, 1.5);

    let curX = 0;
    Array.from(tr.children).forEach((th, cIdx) => {
      const cw = colWidths[cIdx];
      const text = th.innerText.replace(/[▲▼↕*]/g, '').replace(/\s+/g, ' ').trim();
      ctx.font = fontHeader;
      ctx.fillStyle = '#1e293b';
      ctx.textBaseline = 'middle';

      const align = colAlign[cIdx];
      if (align === 'right') {
        ctx.textAlign = 'right';
        ctx.fillText(text, curX + cw - paddingX, curY + headerHeight / 2);
      } else {
        ctx.textAlign = 'left';
        ctx.fillText(text, curX + paddingX, curY + headerHeight / 2);
      }
      curX += cw;
    });
    curY += headerHeight;
  });

  parsedRows.forEach((r, rIdx) => {
    const isEven = rIdx % 2 === 0;
    const defaultRowBg = isEven ? '#ffffff' : '#f8fafc';

    let curX = 0;
    r.forEach((cell, cIdx) => {
      const cw = colWidths[cIdx];
      ctx.fillStyle = cell.bg || defaultRowBg;
      ctx.fillRect(curX, curY, cw, rowHeight);

      ctx.fillStyle = '#f1f5f9';
      ctx.fillRect(curX, curY + rowHeight - 1, cw, 1);
      ctx.fillRect(curX + cw - 1, curY, 1, rowHeight);

      const isBold = cell.weight === 'bold' || parseInt(cell.weight) >= 600;
      ctx.font = cell.isNum ? (isBold ? 'bold ' + fontMono : fontMono) : (isBold ? 'bold ' + fontBody : fontBody);
      ctx.fillStyle = cell.fg || '#0f172a';
      ctx.textBaseline = 'middle';

      const align = colAlign[cIdx];
      if (align === 'right') {
        ctx.textAlign = 'right';
        ctx.fillText(cell.text, curX + cw - paddingX, curY + rowHeight / 2);
      } else {
        ctx.textAlign = 'left';
        ctx.fillText(cell.text, curX + paddingX, curY + rowHeight / 2);
      }
      curX += cw;
    });
    curY += rowHeight;
  });

  canvas.toBlob(function(blob) {
    const a = document.createElement('a');
    const URLObj = window.URL || window.webkitURL || window;
    const pngUrl = blob ? URLObj.createObjectURL(blob) : canvas.toDataURL('image/png');
    a.href = pngUrl;
    a.download = filename || 'table.png';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    if (blob) URLObj.revokeObjectURL(pngUrl);
  }, 'image/png');
}

function exportMatrixAsPng() {
  if (currentMatrixMode === 'heatmap') {
    exportSvgToPng('#heatmap-container svg', 'cross_cancer_target_heatmap.png');
  } else {
    exportTableToPNG('multi-cancer-matrix-table', 'cross_cancer_target_matrix.png');
  }
}

function exportTableToCSV(tableId, filename) {
  const table = typeof tableId === 'string' ? document.getElementById(tableId) : tableId;
  if (!table) return;
  const rows = Array.from(table.querySelectorAll('tr'));
  const csvContent = rows.map(row => {
    const cells = Array.from(row.querySelectorAll('th, td'));
    return cells.map(cell => {
      let text = cell.innerText.replace(/"/g, '""').trim();
      text = text.replace(/[▲▼↕*]/g, '').trim();
      return `"${text}"`;
    }).join(',');
  }).join('\\n');

  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename || 'data.csv';
  a.click();
  URL.revokeObjectURL(url);
}

function exportStandaloneHtml(title, elementId, filename) {
  const el = document.getElementById(elementId);
  if (!el) return;
  const clone = el.cloneNode(true);

  const fullHtml = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${title}</title>
<style>
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 13px; line-height: 1.45; color: #0f172a; background: #f8fafc;
    padding: 2rem; margin: 0;
  }
  .export-header { margin-bottom: 1.5rem; padding-bottom: 1rem; border-bottom: 1px solid #e2e8f0; }
  .export-title { font-size: 18px; font-weight: 700; color: #0f172a; margin-bottom: 4px; }
  .export-meta { font-size: 12px; color: #64748b; }
  .export-card {
    background: #ffffff; border: 1px solid #e2e8f0; border-radius: 6px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05); overflow: hidden; padding: 1rem;
  }
  table { width: 100%; border-collapse: collapse; font-size: 12px; text-align: left; }
  th {
    padding: 0.6rem 0.75rem; border-bottom: 1px solid #e2e8f0; color: #475569;
    font-weight: 600; font-size: 11px; text-transform: uppercase; cursor: pointer;
    user-select: none; background: #f8fafc;
  }
  th:hover { background: #f1f5f9; color: #0f172a; }
  td { padding: 0.55rem 0.75rem; border-bottom: 1px solid #f1f5f9; }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background: #fafafa; }
  .num { text-align: right; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; }
  .tag {
    display: inline-block; padding: 0.15rem 0.4rem; border-radius: 3px;
    font-size: 10.5px; font-weight: 500; font-family: ui-monospace, monospace;
  }
  .tag-fda { background: #dcfce7; color: #15803d; font-weight: 600; }
  .tag-p3 { background: #fef3c7; color: #b45309; font-weight: 600; }
  .tag-p2 { background: #e0e7ff; color: #3730a3; font-weight: 600; }
  .tag-clamped { background: #fee2e2; color: #b91c1c; font-weight: 600; }
  .sort-indicator { font-size: 10px; margin-left: 4px; color: #94a3b8; display: inline-block; }
  .footer-note { font-size: 11px; color: #64748b; margin-top: 1rem; }
  .btn { display: none; }
</style>
</head>
<body>
  <div class="export-header">
    <div class="export-title">${title}</div>
  </div>
  <div class="export-card">
    ${clone.innerHTML}
  </div>
  <div class="footer-note">Standalone interactive document &bull; Click any column header to sort rows ascending / descending.</div>
<' + 'script>
function initStandaloneSort() {
  document.querySelectorAll('table').forEach(table => {
    const headers = table.querySelectorAll('thead th');
    headers.forEach((th, colIdx) => {
      let ind = th.querySelector('.sort-indicator');
      if (!ind) {
        ind = document.createElement('span');
        ind.className = 'sort-indicator';
        ind.textContent = ' ↕';
        th.appendChild(ind);
      }
      th.onclick = function() {
        const isAsc = th.dataset.sortAsc === 'true';
        const newAsc = !isAsc;
        headers.forEach(h => {
          delete h.dataset.sortAsc;
          const si = h.querySelector('.sort-indicator');
          if (si) { si.textContent = ' ↕'; si.style.color = '#94a3b8'; }
        });
        th.dataset.sortAsc = newAsc ? 'true' : 'false';
        ind.textContent = newAsc ? ' ▲' : ' ▼';
        ind.style.color = '#0f172a';

        const tbody = table.querySelector('tbody');
        if (!tbody) return;
        const rows = Array.from(tbody.querySelectorAll('tr'));
        rows.sort((a, b) => {
          const valA = (a.children[colIdx] ? a.children[colIdx].textContent : '').replace(/[%$,+*]/g, '').trim();
          const valB = (b.children[colIdx] ? b.children[colIdx].textContent : '').replace(/[%$,+*]/g, '').trim();
          const numA = parseFloat(valA), numB = parseFloat(valB);
          if (!isNaN(numA) && !isNaN(numB)) {
            return newAsc ? (numA - numB) : (numB - numA);
          }
          return newAsc ? valA.localeCompare(valB) : valB.localeCompare(valA);
        });
        rows.forEach(r => tbody.appendChild(r));
      };
    });
  });
}
window.addEventListener('DOMContentLoaded', initStandaloneSort);
initStandaloneSort();
<' + '/script>
</body>
</html>`;

  const blob = new Blob([fullHtml], { type: 'text/html;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename || 'export.html';
  a.click();
  URL.revokeObjectURL(url);
}

function exportMatrixStandaloneHtml() {
  const isHeatmap = currentMatrixMode === 'heatmap';
  const title = "Cross-Cancer Target Matrix (" + (isHeatmap ? "Heatmap" : "Table") + ")";
  const targetId = isHeatmap ? 'heatmap-container' : 'matrix-table-view';
  exportStandaloneHtml(title, targetId, 'cross_cancer_target_matrix.html');
}

window.addEventListener('DOMContentLoaded', init);
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(
        description="Launch Clinical RNA-Seq Analytics Studio web server."
    )
    parser.add_argument(
        "--db",
        required=True,
        help="Path to DuckDB or SQLite database generated by build_rna_database.py."
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8080,
        help="Port to serve web interface on (default: 8080)."
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host address to bind to (default: 127.0.0.1)."
    )

    args = parser.parse_args()

    store = StudioDataStore(args.db)
    StudioRequestHandler.store = store

    server_address = (args.host, args.port)
    httpd = ThreadingHTTPServer(server_address, StudioRequestHandler)

    print(f"=== Clinical RNA-Seq Analytics Studio Live ===")
    print(f"Database: {store.db_path} ({'DuckDB' if store.use_duckdb else 'SQLite3'})")
    print(f"URL:      http://{args.host}:{args.port}")
    print("Press Ctrl+C to stop.")

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server...")
        httpd.server_close()


if __name__ == "__main__":
    main()
