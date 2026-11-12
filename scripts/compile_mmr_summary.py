#!/usr/bin/env python3
import os
import csv
import sys
import argparse
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Compile individual MMR reports into a single cohort HTML dashboard.")
    parser.add_argument("--results-dir", required=True, help="Base results directory containing sample folders.")
    parser.add_argument("--output", help="Path to write the output HTML dashboard. Defaults to [results-dir]/cohort_mmr_dashboard.html")
    return parser.parse_args()

def get_risk_classification(alterations_str):
    """Determine risk category and color styling from TSV cell contents."""
    if not alterations_str:
        return {
            "label": "Wildtype",
            "class": "risk-wildtype",
            "score": 0
        }
        
    alts = [a.strip().lower() for a in alterations_str.split('/')]
    
    # Biallelic: LOH + (somatic or germline or methylation)
    has_loh = 'loh' in alts
    has_mut_or_meth = ('somatic' in alts) or ('germline' in alts) or ('methylation' in alts)
    
    # Or multiple mutations (e.g. compound heterozygous)
    mut_count = alts.count('somatic') + alts.count('germline')
    
    if (has_loh and has_mut_or_meth) or (mut_count >= 2):
        return {
            "label": "Biallelic Inactivation",
            "class": "risk-biallelic",
            "score": 3,
            "details": alterations_str
        }
    elif has_mut_or_meth or ('cn.loss' in alts) or (has_loh):
        return {
            "label": "Monoallelic Alteration",
            "class": "risk-monoallelic",
            "score": 2,
            "details": alterations_str
        }
    elif 'cn.gain' in alts:
        return {
            "label": "Somatic Gain",
            "class": "risk-low",
            "score": 1,
            "details": alterations_str
        }
    else:
        return {
            "label": alterations_str,
            "class": "risk-unknown",
            "score": 0,
            "details": alterations_str
        }

def get_class_for_classification(classification_str):
    if not classification_str:
        return "risk-unknown"
    cls_l = classification_str.lower()
    if "proficient" in cls_l:
        return "risk-wildtype"
    elif "hereditary" in cls_l or "lynch syndrome" in cls_l:
        return "risk-biallelic"
    elif "double-hit" in cls_l or "lynch-like" in cls_l:
        return "risk-monoallelic"
    elif "sporadic" in cls_l:
        return "risk-low"
    else:
        return "risk-unknown"

def find_clinical_report(results_dir, sample_id, patient):
    """Locate the clinical HTML report for a sample relative to the results dir."""
    keys = list(dict.fromkeys(filter(None, [sample_id, patient])))
    candidates = []
    for k in keys:
        candidates.extend([
            os.path.join(results_dir, patient, "summary", f"{k}_clinical_mmr_report.html"),
            os.path.join(results_dir, patient, patient, "summary", f"{k}_clinical_mmr_report.html"),
            os.path.join(results_dir, sample_id, "summary", f"{k}_clinical_mmr_report.html"),
            os.path.join(results_dir, patient, sample_id, "summary", f"{k}_clinical_mmr_report.html"),
            os.path.join(results_dir, "summary", f"{k}_clinical_mmr_report.html"),
        ])
        
    for cand in candidates:
        if os.path.exists(cand):
            return os.path.relpath(cand, results_dir)
            
    # Shallow search in patient / sample dirs
    for parent in filter(None, [patient, sample_id]):
        for sub in ["summary", os.path.join(parent, "summary")]:
            s_dir = os.path.join(results_dir, parent, sub)
            if os.path.isdir(s_dir):
                matches = [f for f in glob.glob(os.path.join(s_dir, "*_clinical_mmr_report.html")) if not os.path.basename(f).startswith('.')]
                if matches:
                    return os.path.relpath(matches[0], results_dir)
                    
    return None

def find_orange_report(results_dir, sample_id, patient):
    """Locate the Orange PDF report for a sample relative to the results dir."""
    keys = list(dict.fromkeys(filter(None, [sample_id, patient])))
    candidates = []
    for k in keys:
        candidates.extend([
            os.path.join(results_dir, patient, patient, "orange", f"{k}.orange.pdf"),
            os.path.join(results_dir, patient, "orange", f"{k}.orange.pdf"),
            os.path.join(results_dir, sample_id, "orange", f"{k}.orange.pdf"),
            os.path.join(results_dir, patient, sample_id, "orange", f"{k}.orange.pdf"),
            os.path.join(results_dir, "orange", f"{k}.orange.pdf"),
        ])
    for cand in candidates:
        if os.path.exists(cand):
            return os.path.relpath(cand, results_dir)
            
    # Shallow glob in orange directories
    for parent in filter(None, [patient, sample_id]):
        for orange_sub in [
            os.path.join(results_dir, parent, parent, "orange"),
            os.path.join(results_dir, parent, "orange"),
            os.path.join(results_dir, parent, "orange", "output")
        ]:
            if os.path.isdir(orange_sub):
                matches = [f for f in glob.glob(os.path.join(orange_sub, "*.orange.pdf")) if not os.path.basename(f).startswith('.')]
                if matches:
                    return os.path.relpath(matches[0], results_dir)
                    
    return None

def compile_tumor_fraction_comparisons(results_dir):
    """Scan all sample directories, read tumor_fraction_comparison.tsv, and merge into a cohort-wide TSV."""
    output_path = os.path.join(results_dir, "cohort_tumor_fraction_comparison.tsv")
    
    # Locate all individual comparison TSV files (support both ichorCNA and ichorCNA_panel)
    patterns = [
        os.path.join(results_dir, "*", "ichorCNA", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "*", "ichorCNA", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "ichorCNA", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "ichorCNA_t_vs_n", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "ichorCNA_plasma_vs_n", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "ichorCNA_tumor_only", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "ichorCNA_plasma_only", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "ichorCNA_panel", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "*", "*", "ichorCNA_panel", "tumor_fraction_comparison.tsv"),
        os.path.join(results_dir, "ichorCNA_panel", "tumor_fraction_comparison.tsv"),
    ]
    matches = []
    for pat in patterns:
        matches.extend(glob.glob(pat))
    matches = list(dict.fromkeys(matches))
    
    if not matches:
        return
        
    compiled_rows = []
    for match in matches:
        try:
            with open(match, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    compiled_rows.append(row)
        except Exception as e:
            print(f"Warning: Failed to read {match}: {e}", file=sys.stderr)
            
    if not compiled_rows:
        return
        
    try:
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ["sample_id", "sample_type", "panel_tumor_fraction", "provided_tumor_fraction", "hmf_tumor_fraction", "delfi_tumor_fraction", "delfi_pa_score"]
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore', delimiter='\t')
            writer.writeheader()
            for r in compiled_rows:
                writer.writerow(r)
        print(f"Cohort tumor fraction comparison TSV written to: {output_path}")
    except Exception as e:
        print(f"Warning: Failed to write cohort comparison TSV {output_path}: {e}", file=sys.stderr)

def main():
    args = parse_args()
    results_dir = os.path.abspath(args.results_dir)
    
    if not os.path.exists(results_dir):
        print(f"Error: Results directory does not exist: {results_dir}", file=sys.stderr)
        sys.exit(1)
        
    output_html = args.output
    if not output_html:
        output_html = os.path.join(results_dir, "cohort_mmr_dashboard.html")
        
    # Find summary TSVs
    # Check for individual sample summaries first across cohort subdirectories
    individual_tsvs = glob.glob(os.path.join(results_dir, "*", "final_summary.tsv"))
    if not individual_tsvs:
        individual_tsvs = glob.glob(os.path.join(results_dir, "*", "summary", "*final_summary.tsv"))
    individual_tsvs = [p for p in individual_tsvs if not any(part.startswith('.') for part in p.split(os.sep))]
    individual_tsvs = list(dict.fromkeys(individual_tsvs))
    
    if len(individual_tsvs) > 1:
        print(f"Found {len(individual_tsvs)} individual summaries. Merging...", file=sys.stderr)
        rows = []
        headers = []
        for path in sorted(individual_tsvs):
            with open(path, mode='r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                if not headers:
                    headers = reader.fieldnames
                for r in reader:
                    r['_source_path'] = path
                    rows.append(r)
    else:
        # First check root final_summary.tsv
        root_tsv = os.path.join(results_dir, "final_summary.tsv")
        tsv_path = None
        if os.path.exists(root_tsv):
            tsv_path = root_tsv
        else:
            pattern = os.path.join(results_dir, "summary", "*final_summary.tsv")
            matches = glob.glob(pattern)
            if matches:
                tsv_path = matches[0]
            elif individual_tsvs:
                tsv_path = individual_tsvs[0]
                
        if not tsv_path:
            print("Error: Could not locate any final_summary.tsv files.", file=sys.stderr)
            sys.exit(1)
            
        print(f"Reading summary from: {tsv_path}", file=sys.stderr)
        rows = []
        with open(tsv_path, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            headers = reader.fieldnames
            for r in reader:
                r['_source_path'] = tsv_path
                rows.append(r)
                
    # Deduplicate rows by (Patient, Sample), prioritizing rows with Purity.Source populated
    dedup_dict = {}
    for r in rows:
        key = (r.get('Patient'), r.get('Sample'))
        if key not in dedup_dict:
            dedup_dict[key] = r
        else:
            existing = dedup_dict[key]
            if r.get('Purity.Source') and r.get('Purity.Source') != 'NA' and (not existing.get('Purity.Source') or existing.get('Purity.Source') == 'NA'):
                dedup_dict[key] = r
    rows = list(dedup_dict.values())
                
    if not rows:
        print("Error: No data rows found to compile.", file=sys.stderr)
        sys.exit(1)
        
    # Filter out normal baseline rows
    filtered_rows = []
    for r in rows:
        classification = str(r.get('Clinical.Classification', '')).strip().lower()
        sample_name = str(r.get('Sample', '')).strip()
        if classification == 'normal baseline' or sample_name.endswith('_Normal'):
            continue
        filtered_rows.append(r)
    rows = filtered_rows
    
    if not rows:
        print("Error: No tumor data rows remaining after filtering normal baselines.", file=sys.stderr)
        sys.exit(1)
        
    # Apply cohort frequency filter for germline deletions (technical panel artifacts)
    n_total_samples = len(rows)
    genes_to_strip = []
    core_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2']
    other_headers = [h for h in headers if h not in ['Patient', 'Sample', 'Tumour.Fraction', 'MSI.Status', 'MSI.Indels.Per.Mb', 'MSI.Unstable.Sites.Proportion', 'Clinical.Classification', 'Clinical.Confidence', 'Clinical.Notes'] + core_genes]
    all_genes_list = core_genes + other_headers

    if n_total_samples > 0:
        for gene in all_genes_list:
            if gene not in headers:
                continue
            del_count = 0
            for r in rows:
                val = r.get(gene, '')
                if val and 'germline.deletion' in str(val):
                    del_count += 1
            freq = del_count / n_total_samples
            if freq > 0.25:
                print(f"Cohort frequency filter: stripping germline.deletion for {gene} (called in {del_count}/{n_total_samples} samples, {freq:.1%})", file=sys.stderr)
                genes_to_strip.append(gene)

    if genes_to_strip:
        for r in rows:
            # 1. Strip from gene columns
            for gene in genes_to_strip:
                val = r.get(gene, '')
                if val and 'germline.deletion' in str(val):
                    parts = [p.strip() for p in str(val).split('/') if p.strip() != 'germline.deletion']
                    new_val = '/'.join(parts) if parts else 'Normal'
                    r[gene] = new_val
            
            # 2. Update Clinical.Notes / Clinical.Classification to reflect the stripped genes
            notes = r.get('Clinical.Notes', '')
            if notes and "Pathogenic germline variant detected in:" in str(notes):
                try:
                    prefix, gene_list_str = str(notes).split("Pathogenic germline variant detected in:", 1)
                    import re
                    words = re.findall(r'\b[A-Z0-9]+\b', gene_list_str)
                    kept_genes = [w for w in words if w not in genes_to_strip]
                    
                    if kept_genes:
                        new_gene_list = ", ".join(kept_genes) + "."
                        r['Clinical.Notes'] = f"{prefix}Pathogenic germline variant detected in: {new_gene_list}"
                    else:
                        parts = prefix.split(" | ", 1)
                        r['Clinical.Notes'] = parts[1] if len(parts) > 1 else "Tumor displays intact mismatch repair mechanisms."
                except Exception:
                    pass

        # Rewrite modified individual summaries to keep them in sync
        if 'individual_tsvs' in locals() and individual_tsvs:
            rows_by_path = {}
            for r in rows:
                path = r.get('_source_path')
                if path:
                    if path not in rows_by_path:
                        rows_by_path[path] = []
                    rows_by_path[path].append(r)
            
            for path, path_rows in rows_by_path.items():
                try:
                    clean_headers = [h for h in headers if h != '_source_path']
                    with open(path, mode='w', newline='', encoding='utf-8') as f:
                        writer = csv.DictWriter(f, fieldnames=clean_headers, delimiter='\t')
                        writer.writeheader()
                        for pr in path_rows:
                            clean_pr = {k: pr.get(k, '') for k in clean_headers}
                            writer.writerow(clean_pr)
                except Exception as e:
                    print(f"Warning: Failed to update individual summary {path}: {e}", file=sys.stderr)
        
    # Core MMR genes we display first, followed by others
    core_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2']
    other_headers = [h for h in headers if h not in ['Patient', 'Sample', 'Tumour.Fraction', 'MSI.Status', 'MSI.Indels.Per.Mb', 'MSI.Unstable.Sites.Proportion', 'Clinical.Classification', 'Clinical.Confidence', 'Clinical.Notes'] + core_genes]
    display_genes = core_genes + other_headers
    # Preprocess rows to determine Biallelic Inactivation status, alteration counts, and MSI status sort keys
    core_mmr_genes = ['MLH1', 'MSH2', 'MSH6', 'PMS2']
    for row in rows:
        has_biallelic_mmr = False
        num_alts = 0
        for gene in display_genes:
            cell_val = row.get(gene, '')
            if cell_val and cell_val.strip() and cell_val.strip().lower() != 'wildtype':
                num_alts += 1
                if gene in core_mmr_genes:
                    risk = get_risk_classification(cell_val)
                    if risk['class'] == 'risk-biallelic':
                        has_biallelic_mmr = True
                        
        row['_has_biallelic_mmr'] = has_biallelic_mmr
        row['MMR_Biallelic_Inactivation'] = 'Yes' if has_biallelic_mmr else 'No'
        row['_num_alterations'] = num_alts
        
        # MSI sorting value: MSI (0) -> MSS (1) -> other/missing (2)
        msi = row.get('MSI.Status', '').strip().upper()
        if msi == 'MSI':
            row['_msi_sort'] = 0
        elif msi == 'MSS':
            row['_msi_sort'] = 1
        else:
            row['_msi_sort'] = 2

    # Sort: Biallelic Inactivation (Yes/True first, so 0 first), then MSI status, then number of alterations (descending)
    rows.sort(key=lambda r: (
        0 if r['_has_biallelic_mmr'] else 1,
        r['_msi_sort'],
        -r['_num_alterations']
    ))
    
    # Generate HTML content
    html_lines = []
    html_lines.append("""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MMR & Lynch Syndrome Cohort Summary Dashboard</title>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&family=JetBrains+Mono:wght@400;700&display=swap" rel="stylesheet">
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --border-color: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            
            --color-biallelic: #fecaca;
            --text-biallelic: #991b1b;
            --border-biallelic: #f87171;
            
            --color-monoallelic: #ffedd5;
            --text-monoallelic: #9a3412;
            --border-monoallelic: #fb923c;
            
            --color-low: #e0f2fe;
            --text-low: #0369a1;
            --border-low: #38bdf8;
            
            --color-wildtype: #d1fae5;
            --text-wildtype: #065f46;
            --border-wildtype: #34d399;
        }
        
        body {
            font-family: 'Outfit', sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            margin: 0;
            padding: 24px;
        }
        
        .header {
            margin-bottom: 32px;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 24px;
        }
        
        h1 {
            font-size: 2.2rem;
            font-weight: 700;
            margin: 0 0 8px 0;
            background: linear-gradient(135deg, #38bdf8, #818cf8);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        
        .subtitle {
            color: var(--text-secondary);
            font-size: 1.1rem;
            margin: 0;
        }
        
        .dashboard-container {
            background-color: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            overflow: hidden;
            box-shadow: 0 10px 30px rgba(0, 0, 0, 0.25);
        }
        
        .table-responsive {
            overflow-x: auto;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            text-align: left;
        }
        
        th {
            background-color: #0f172a;
            color: var(--text-secondary);
            font-weight: 600;
            font-size: 0.95rem;
            padding: 16px 20px;
            border-bottom: 1px solid var(--border-color);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            cursor: pointer;
            user-select: none;
            transition: background-color 0.2s ease, color 0.2s ease;
        }
        
        th:hover {
            background-color: #1e293b !important;
            color: #38bdf8 !important;
        }
        
        td {
            padding: 16px 20px;
            border-bottom: 1px solid var(--border-color);
            font-size: 1rem;
            vertical-align: middle;
        }
        
        tr:hover td {
            background-color: #273549;
        }
        
        .sample-link {
            color: #38bdf8;
            text-decoration: none;
            font-weight: 600;
            border-bottom: 1px dashed transparent;
            transition: all 0.2s ease;
        }
        
        .sample-link:hover {
            color: #7dd3fc;
            border-bottom-color: #7dd3fc;
        }
        
        .badge {
            display: inline-block;
            padding: 6px 12px;
            border-radius: 9999px;
            font-size: 0.85rem;
            font-weight: 700;
            letter-spacing: 0.02em;
        }
        
        .badge-box {
            display: inline-block;
            padding: 6px 12px;
            border-radius: 6px;
            font-size: 0.85rem;
            font-weight: 700;
            letter-spacing: 0.02em;
        }
        
        .badge-msi {
            background-color: #fef08a;
            color: #854d0e;
        }
        
        .badge-mss {
            background-color: #e2e8f0;
            color: #475569;
        }
        
        /* Cell Risk Styles */
        .cell-risk {
            border-radius: 8px;
            padding: 10px 14px;
            border: 1px solid transparent;
            font-size: 0.9rem;
            line-height: 1.3;
        }
        
        .risk-biallelic {
            background-color: var(--color-biallelic);
            color: var(--text-biallelic);
            border-color: var(--border-biallelic);
        }
        
        .risk-monoallelic {
            background-color: var(--color-monoallelic);
            color: var(--text-monoallelic);
            border-color: var(--border-monoallelic);
        }
        
        .risk-low {
            background-color: var(--color-low);
            color: var(--text-low);
            border-color: var(--border-low);
        }
        
        .risk-wildtype {
            background-color: var(--color-wildtype);
            color: var(--text-wildtype);
            border-color: var(--border-wildtype);
        }
        
        .details-sub {
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.75rem;
            opacity: 0.75;
            margin-top: 4px;
            display: block;
        }
        
        .legend {
            display: flex;
            gap: 16px;
            margin-top: 24px;
            flex-wrap: wrap;
            background-color: var(--card-bg);
            border: 1px solid var(--border-color);
            padding: 16px 24px;
            border-radius: 12px;
        }
        
        .legend-item {
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 0.9rem;
            color: var(--text-secondary);
        }
        
        .legend-color {
            width: 16px;
            height: 16px;
            border-radius: 4px;
            border: 1px solid transparent;
        }
    </style>
    <script>
    document.addEventListener('DOMContentLoaded', () => {
        const table = document.querySelector('table');
        if (!table) return;
        const headers = table.querySelectorAll('th');
        
        headers.forEach((header, index) => {
            header.addEventListener('click', () => {
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                const isAscending = header.dataset.sortDirection !== 'asc';
                
                // Clear sorting icons and states on other headers
                headers.forEach(h => {
                    if (h !== header) {
                        delete h.dataset.sortDirection;
                        h.innerHTML = h.innerHTML.replace(/ ▲| ▼/g, '');
                    }
                });
                
                rows.sort((rowA, rowB) => {
                    const cellA = rowA.cells[index].innerText.trim();
                    const cellB = rowB.cells[index].innerText.trim();
                    
                    // Extract numeric values if cells contain numbers
                    const cleanA = cellA.replace(/[^0-9.-]/g, '');
                    const cleanB = cellB.replace(/[^0-9.-]/g, '');
                    const valA = parseFloat(cleanA);
                    const valB = parseFloat(cleanB);
                    
                    if (cleanA !== '' && cleanB !== '' && !isNaN(valA) && !isNaN(valB)) {
                        return isAscending ? valA - valB : valB - valA;
                    }
                    
                    return isAscending 
                        ? cellA.localeCompare(cellB, undefined, {numeric: true, sensitivity: 'base'}) 
                        : cellB.localeCompare(cellA, undefined, {numeric: true, sensitivity: 'base'});
                });
                
                header.dataset.sortDirection = isAscending ? 'asc' : 'desc';
                header.innerHTML = header.innerHTML.replace(/ ▲| ▼/g, '') + (isAscending ? ' ▲' : ' ▼');
                
                tbody.innerHTML = '';
                rows.forEach(row => tbody.appendChild(row));
            });
        });
    });
    </script>
</head>
<body>
    <div class="header">
        <h1>MMR & Lynch Syndrome Cohort Summary Dashboard</h1>
        <p class="subtitle">Compiled clinical findings and mismatch repair pathway status across the cohort</p>
    </div>
    
    <div class="dashboard-container">
        <div class="table-responsive">
            <table>
                <thead>
                    <tr>
                        <th>Patient</th>
                        <th>Sample</th>
                        <th>Orange Report</th>
                        <th>Tumour Fraction</th>
                        <th>MSI Status</th>
                        <th>MMR gene Biallelic Inactivation</th>
                        <th>Clinical Classification</th>
                        <th>Confidence</th>
""")
    
    # Table headers for genes
    for gene in display_genes:
        html_lines.append(f"                        <th>{gene}</th>\n")
        
    html_lines.append("""                    </tr>
                </thead>
                <tbody>
""")
    
    # Table rows
    for row in rows:
        patient = row.get('Patient', '')
        sample = row.get('Sample', '')
        tf = row.get('Tumour.Fraction', '')
        msi = row.get('MSI.Status', '')
        
        # Find clickable link
        link = find_clinical_report(results_dir, sample, patient)
        if link:
            sample_cell = f'<a class="sample-link" href="{link}" target="_blank">{sample}</a>'
        else:
            sample_cell = sample
            
        # Orange Report badge
        orange_link = find_orange_report(results_dir, sample, patient)
        if orange_link:
            orange_cell = f'<a href="{orange_link}" target="_blank" class="badge" style="background-color: #ffedd5; color: #c2410c; border: 1px solid #fdba74; text-decoration: none; font-weight: 600; display: inline-flex; align-items: center; gap: 4px;">📄 Orange PDF</a>'
        else:
            orange_cell = '<span style="color: var(--text-secondary); font-size: 0.8rem;">N/A</span>'

        # MSI badge
        msi_class = "badge-msi" if msi.upper() == "MSI" else "badge-mss"
        msi_cell = f'<span class="badge {msi_class}">{msi}</span>' if msi else "N/A"
        
        # Biallelic Inactivation badge
        bi_class = "risk-biallelic" if row['_has_biallelic_mmr'] else "badge-mss"
        bi_cell = f'<span class="badge {bi_class}">{row["MMR_Biallelic_Inactivation"]}</span>'
        
        # Clinical classification badges
        clin_class = row.get('Clinical.Classification', 'N/A')
        clin_badge_class = get_class_for_classification(clin_class)
        clin_cell = f'<span class="badge-box {clin_badge_class}">{clin_class}</span>'
        conf_cell = row.get('Clinical.Confidence', 'N/A')
        if str(conf_cell).startswith("Low"):
            conf_cell_styled = f'<span class="badge" style="background-color: #fef2f2; color: #991b1b; border: 1px solid #fee2e2;">⚠️ {conf_cell}</span>'
        else:
            conf_cell_styled = conf_cell
            
        tf_cell = tf if tf else 'N/A'
        if tf and tf != 'NA':
            try:
                tf_val = float(tf)
                if tf_val < 0.3:
                    tf_cell = f'<span class="badge" style="background-color: #fffbeb; color: #92400e; border: 1px solid #fef3c7; font-weight: bold;">⚠️ {tf} (Low)</span>'
            except ValueError:
                pass
        
        html_lines.append(f"""                    <tr>
                        <td>{patient}</td>
                        <td>{sample_cell}</td>
                        <td>{orange_cell}</td>
                        <td>{tf_cell}</td>
                        <td>{msi_cell}</td>
                        <td>{bi_cell}</td>
                        <td>{clin_cell}</td>
                        <td>{conf_cell_styled}</td>
""")
        
        # Gene status cells
        for gene in display_genes:
            cell_val = row.get(gene, '')
            risk = get_risk_classification(cell_val)
            
            if risk['label'] == "Wildtype":
                html_lines.append(f"""                        <td>
                            <div class="cell-risk {risk['class']}">
                                <strong>Wildtype</strong>
                            </div>
                        </td>
""")
            else:
                details_html = f'<span class="details-sub">{risk["details"]}</span>' if "details" in risk else ""
                html_lines.append(f"""                        <td>
                            <div class="cell-risk {risk['class']}">
                                <strong>{risk['label']}</strong>
                                {details_html}
                            </div>
                        </td>
""")
                
        html_lines.append("                    </tr>\n")
        
    html_lines.append("""                </tbody>
            </table>
        </div>
    </div>
    
    <div class="legend">
        <div class="legend-item">
            <div class="legend-color risk-biallelic"></div>
            <span>Biallelic Inactivation (High Risk)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color risk-monoallelic"></div>
            <span>Monoallelic Alteration / CN Loss (Moderate Risk)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color risk-low"></div>
            <span>Somatic Gain (Low Risk)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color risk-wildtype"></div>
            <span>Wildtype / No Alterations (Normal)</span>
        </div>
    </div>
</body>
</html>
""")
    
    # Write output
    with open(output_html, 'w', encoding='utf-8') as f:
        f.writelines(html_lines)
        
    print(f"Cohort MMR dashboard HTML written to: {output_html}")

    # Write merged TSV file
    output_tsv = os.path.join(results_dir, "cohort_final_summary.tsv")
    try:
        with open(output_tsv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            for row in rows:
                clean_row = {k: row.get(k, '') for k in headers}
                writer.writerow(clean_row)
        print(f"Cohort merged TSV written to: {output_tsv}")
        
        # Merge individual sample tumor fraction comparison TSV files
        compile_tumor_fraction_comparisons(results_dir)
    except Exception as e:
        print(f"Warning: Failed to write merged TSV {output_tsv}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
