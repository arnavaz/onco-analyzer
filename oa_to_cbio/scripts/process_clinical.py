import json
import sys
import os
import re
import argparse
import csv

def sanitize_attribute_name(name):
    # Convert to uppercase
    name = name.upper()
    # Replace non-alphanumeric characters with underscore
    name = re.sub(r'[^A-Z0-9_]', '_', name)
    # Collapse multiple underscores
    name = re.sub(r'_+', '_', name)
    # Strip leading/trailing underscores
    name = name.strip('_')
    return name

def parse_teal_file(teal_path):
    if not os.path.exists(teal_path):
        return {}
    try:
        with open(teal_path, 'r') as f:
            headers = [h.strip() for h in f.readline().rstrip('\r\n').split('\t')]
            values = [v.strip() for v in f.readline().rstrip('\r\n').split('\t')]
            data = {}
            for h, v in zip(headers, values):
                if v == "" or v.upper() == "NA" or v.upper() == "NAN":
                    data[h] = "NA"
                else:
                    data[h] = v
            return data
    except Exception as e:
        print(f"Warning: Failed to parse TEAL file {teal_path}: {e}", file=sys.stderr)
        return {}

def find_cuppa_tsv(tumor_id, sample_dir):
    if not sample_dir or not os.path.isdir(sample_dir):
        return None
    # 1. Direct path in cuppa/ subfolder
    cand = os.path.join(sample_dir, "cuppa", f"{tumor_id}.cuppa.pred_summ.tsv")
    if os.path.exists(cand):
        return cand
    # 2. Any .cuppa.pred_summ.tsv in cuppa/ subfolder
    cuppa_dir = os.path.join(sample_dir, "cuppa")
    if os.path.isdir(cuppa_dir):
        try:
            for f in sorted(os.listdir(cuppa_dir)):
                if f.endswith(".cuppa.pred_summ.tsv") and not f.startswith("."):
                    return os.path.join(cuppa_dir, f)
        except OSError:
            pass
    # 3. Direct in sample_dir
    cand = os.path.join(sample_dir, f"{tumor_id}.cuppa.pred_summ.tsv")
    if os.path.exists(cand):
        return cand
    # 4. Sibling cuppa/ or parent search
    parent_dir = os.path.dirname(sample_dir.rstrip("/"))
    if os.path.isdir(parent_dir):
        cand = os.path.join(parent_dir, "cuppa", f"{tumor_id}.cuppa.pred_summ.tsv")
        if os.path.exists(cand):
            return cand
        cand = os.path.join(parent_dir, tumor_id, "cuppa", f"{tumor_id}.cuppa.pred_summ.tsv")
        if os.path.exists(cand):
            return cand
    return None

def parse_cuppa_pred_summ(cuppa_tsv_path):
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
            likelihood = target_row.get("pred_prob_1")
            preds = []
            for i in range(1, 5):
                pc = target_row.get(f"pred_class_{i}")
                pp = target_row.get(f"pred_prob_{i}")
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

def format_numeric(val, decimals=None):
    if val is None or str(val).upper() in ["NA", "NAN", "NULL", ""]:
        return "NA"
    try:
        num = float(val)
        if decimals is None or decimals == 0:
            return str(round(num))
        else:
            return f"{num:.{decimals}f}"
    except ValueError:
        return "NA"

def get_value(d, *keys):
    for k in keys:
        if d is None:
            return "NA"
        d = d.get(k)
    return d if d is not None else "NA"

def resolve_sample_dir(tumor_id, json_in):
    # Try environment variable first
    s_dir = os.environ.get("SAMPLE_DIR")
    if s_dir:
        s_dir = os.path.abspath(s_dir)
        # Check if the folder contains purple or orange or linx subdirectories directly
        if os.path.exists(os.path.join(s_dir, "purple")) or os.path.exists(os.path.join(s_dir, "linx")):
            return s_dir
        
        # Sibling directory resolving: if not direct, check if sibling/child named after tumor_id exists
        parent = os.path.dirname(s_dir.rstrip("/"))
        sibling = os.path.join(parent, tumor_id)
        if os.path.exists(os.path.join(sibling, "purple")) or os.path.exists(os.path.join(sibling, "linx")):
            return sibling
            
        child = os.path.join(s_dir, tumor_id)
        if os.path.exists(os.path.join(child, "purple")) or os.path.exists(os.path.join(child, "linx")):
            return child
            
        return s_dir
        
    # Fallback to json path
    if json_in and os.path.exists(json_in) and json_in.upper() != "NA":
        s_dir = os.path.dirname(os.path.dirname(os.path.abspath(json_in)))
        if os.path.exists(os.path.join(s_dir, "purple")) or os.path.exists(os.path.join(s_dir, "linx")):
            return s_dir
            
    return os.path.abspath(os.getcwd())

def main():
    parser = argparse.ArgumentParser(description="Process clinical data for cBioPortal.")
    parser.add_argument("sample_id", help="Sample ID (Tumor or Normal)")
    parser.add_argument("orange_json", help="Path to orange JSON file (or empty/NA)")
    parser.add_argument("output_tsv", help="Path to output clinical TSV file")
    parser.add_argument("--normal-id", default=None, help="Matched normal Sample ID")
    parser.add_argument("--is-normal", action="store_true", help="Flag indicating the sample is a normal control")

    args = parser.parse_args()

    tumor_id = args.sample_id
    json_in = args.orange_json
    out_tsv = args.output_tsv
    normal_id = args.normal_id
    is_normal = args.is_normal

    # These cohorts encode the donor ID in the first three sample ID components.
    patient_id = tumor_id
    if tumor_id.startswith(("PANCAN-", "DYNA-", "MOHCCN-")):
        parts = tumor_id.split('-')
        if len(parts) >= 3:
            patient_id = "-".join(parts[:3])

    metrics = {
        "SAMPLE_ID": tumor_id,
        "PATIENT_ID": patient_id,
        "IS_NORMAL": "Yes" if is_normal else "No"
    }

    cuppa_cancer_type = None
    cuppa_likelihood = None
    cuppa_predictions = []

    if json_in and os.path.exists(json_in) and json_in.upper() != "NA":
        try:
            with open(json_in, 'r') as f:
                data = json.load(f)
                
                # Use the true sample ID from the JSON if available
                true_sample_id = data.get("sampleId")
                if true_sample_id:
                    metrics["SAMPLE_ID"] = true_sample_id
                
                # Core metrics
                metrics["PURITY"] = get_value(data, "purple", "fit", "purity")
                metrics["PLOIDY"] = get_value(data, "purple", "fit", "ploidy")
                metrics["MSISTATUS"] = get_value(data, "purple", "characteristics", "microsatelliteStatus")
                metrics["MSISCORE"] = get_value(data, "purple", "characteristics", "microsatelliteIndelsPerMb")
                metrics["TMBPERMB"] = get_value(data, "purple", "characteristics", "tumorMutationalBurdenPerMb")
                
                # 1. HLA Alleles
                alleles = get_value(data, "lilac", "alleles")
                if isinstance(alleles, list):
                    counts = {}
                    for a in alleles:
                        allele_str = a.get('allele', '')
                        if allele_str:
                            hla_type = allele_str.split('*')[0] if '*' in allele_str else "UNKNOWN"
                            counts[hla_type] = counts.get(hla_type, 0) + 1
                            metrics[f"HLA_{hla_type}_{counts[hla_type]}"] = allele_str

                # 2. Immune Escape
                escape_flags = data.get("immuneEscape", {})
                escapes = []
                if escape_flags.get("hasHlaEscape"): escapes.append("HLA Escape")
                if escape_flags.get("hasAntigenPresentationPathwayEscape"): escapes.append("Antigen Presentation Escape")
                if escape_flags.get("hasIFNGammaPathwayEscape"): escapes.append("IFN-Gamma Escape")
                if escape_flags.get("hasPDL1OverexpressionEscape"): escapes.append("PDL1 Overexpression")
                if escape_flags.get("hasCD58InactivationEscape"): escapes.append("CD58 Inactivation")
                if escape_flags.get("hasEpigeneticSETDB1Escape"): escapes.append("Epigenetic SETDB1 Escape")
                metrics["ESCAPE_MECHANISM"] = ", ".join(escapes) if escapes else "None"
                
                # 3. Coverage Stats
                cov = data.get("tumorSample", {}).get("metrics", {})
                if cov:
                    metrics["MEAN_COVERAGE"] = cov.get("meanCoverage", "NA")
                    metrics["MEDIAN_COVERAGE"] = cov.get("medianCoverage", "NA")
                    
                # 4. Signatures
                sigs = data.get("sigAllocations", [])
                for sig in sigs:
                    s_name = sig.get("signature", "")
                    s_et = sig.get("etiology", "")
                    s_pct = sig.get("percent")
                    if s_name and s_pct is not None:
                        col_name = f"{s_name}({s_et})" if s_et else s_name
                        metrics[col_name] = f"{float(s_pct)*100:.2f}"
                        
                # 5 & 6. Viruses
                vi = data.get("virusInterpreter", {})
                rep_v = vi.get("reportableViruses", [])
                all_v = vi.get("allViruses", [])
                
                if not all_v and rep_v:
                    all_v = rep_v
                    
                if all_v:
                    for idx, v in enumerate(all_v, 1):
                        prefix = f"VIRUS_{idx}"
                        if isinstance(v, dict):
                            v_name = v.get("virus", v.get("name", "Unknown"))
                            metrics[f"{prefix}_ORGANISM"] = v_name
                            
                            is_driver = "No"
                            for rv in rep_v:
                                if isinstance(rv, dict) and rv.get("virus", rv.get("name")) == v_name:
                                    is_driver = "Yes"
                                elif rv == v_name:
                                    is_driver = "Yes"
                            metrics[f"{prefix}_IS_DRIVER"] = is_driver
                            
                            for k, val in v.items():
                                if k not in ["virus", "name"]:
                                    metrics[f"{prefix}_{k.upper()}"] = str(val)
                        else:
                            metrics[f"{prefix}_ORGANISM"] = str(v)
                            metrics[f"{prefix}_IS_DRIVER"] = "Yes" if v in rep_v else "No"
                else:
                    metrics["VIRUS_1_ORGANISM"] = "None"

                # 6. CUPPA cancer type predictions
                cuppa = data.get("cuppa", {})
                if isinstance(cuppa, dict):
                    cuppa_best = cuppa.get("bestPrediction", {}) if isinstance(cuppa.get("bestPrediction"), dict) else {}
                    cuppa_cancer_type = cuppa_best.get("cancerType")
                    cuppa_likelihood = cuppa_best.get("likelihood")
                    cuppa_predictions = cuppa.get("predictions", [])

        except Exception as e:
            print(f"Warning: Failed to parse {json_in}: {e}", file=sys.stderr)

    # 7. TEAL Telomere Length Metrics Extraction
    sample_dir = resolve_sample_dir(tumor_id, json_in)

    teal_file = os.path.join(sample_dir, "teal", f"{tumor_id}.teal.tellength.tsv")
    teal_data = parse_teal_file(teal_file) if os.path.exists(teal_file) else {}

    if teal_data:
        metrics["MEAN_TELOMERE_LENGTH_TUMOUR_BP"] = format_numeric(teal_data.get("finalTelomereLength"), 1)
        metrics["MEAN_TELOMERE_LENGTH_TUMOUR_RAW_BP"] = format_numeric(teal_data.get("rawTelomereLength"), 1)
        metrics["DUPLICATE_READ_PROPORTION_TUMOUR"] = format_numeric(teal_data.get("duplicateProportion"), 3)

    # 8. Paired Normal Baseline & Delta Metrics
    if not is_normal:
        normal_teal_file = None
        candidates = []
        teal_dir = os.path.join(sample_dir, "teal")

        if normal_id and normal_id != "UNKNOWN_NORMAL":
            # 1. Primary: Same teal directory as tumor (standard oncoanalyser paired output)
            candidates.append(os.path.join(teal_dir, f"{normal_id}.teal.tellength.tsv"))
            # 2. Sibling normal directory (if flat or separate sample output directories)
            parent_dir = os.path.dirname(sample_dir.rstrip("/"))
            candidates.append(os.path.join(parent_dir, normal_id, "teal", f"{normal_id}.teal.tellength.tsv"))
            # 3. Fuzzy search for normal_id in teal_dir
            if os.path.isdir(teal_dir):
                try:
                    for f in sorted(os.listdir(teal_dir)):
                        if normal_id in f and f.endswith(".teal.tellength.tsv") and not f.startswith(f"{tumor_id}."):
                            candidates.append(os.path.join(teal_dir, f))
                except OSError:
                    pass

        # 4. Auto-discovery fallback: any .teal.tellength.tsv in teal_dir that is NOT the tumor file
        if os.path.isdir(teal_dir):
            try:
                for f in sorted(os.listdir(teal_dir)):
                    if f.endswith(".teal.tellength.tsv") and not f.startswith(f"{tumor_id}."):
                        candidates.append(os.path.join(teal_dir, f))
            except OSError:
                pass

        # 5. Sibling directory fallback: look in sibling directories matching normal naming
        parent_dir = os.path.dirname(sample_dir.rstrip("/"))
        if os.path.isdir(parent_dir):
            try:
                for sibling in sorted(os.listdir(parent_dir)):
                    if sibling != tumor_id and any(k in sibling.lower() for k in ["normal", "-n", "_n", "germline"]):
                        sib_teal = os.path.join(parent_dir, sibling, "teal")
                        if os.path.isdir(sib_teal):
                            for f in sorted(os.listdir(sib_teal)):
                                if f.endswith(".teal.tellength.tsv"):
                                    candidates.append(os.path.join(sib_teal, f))
            except OSError:
                pass

        for cand in candidates:
            if cand and os.path.exists(cand):
                normal_teal_file = cand
                break

        normal_teal_data = parse_teal_file(normal_teal_file) if normal_teal_file else {}
        
        if normal_teal_data:
            metrics["MEAN_TELOMERE_LENGTH_NORMAL_BP"] = format_numeric(normal_teal_data.get("finalTelomereLength"), 1)
            metrics["DUPLICATE_READ_PROPORTION_NORMAL"] = format_numeric(normal_teal_data.get("duplicateProportion"), 3)
            
            # Calculate Delta: finalTelomereLength_Tumor - finalTelomereLength_Normal
            t_len = teal_data.get("finalTelomereLength") if teal_data else None
            n_len = normal_teal_data.get("finalTelomereLength")
            if t_len and n_len and str(t_len).upper() != "NA" and str(n_len).upper() != "NA":
                try:
                    delta = float(t_len) - float(n_len)
                    metrics["MEAN_TELOMERE_LENGTH_DELTA_BP"] = format_numeric(delta, 1)
                except ValueError:
                    metrics["MEAN_TELOMERE_LENGTH_DELTA_BP"] = "NA"
            else:
                metrics["MEAN_TELOMERE_LENGTH_DELTA_BP"] = "NA"
        else:
            if teal_data:
                metrics["MEAN_TELOMERE_LENGTH_NORMAL_BP"] = "NA"
                metrics["DUPLICATE_READ_PROPORTION_NORMAL"] = "NA"
                metrics["MEAN_TELOMERE_LENGTH_DELTA_BP"] = "NA"

    # 9. ecDNA / Double Minute Metrics Extraction
    clusters_file = os.path.join(sample_dir, "linx", "somatic_annotations", f"{tumor_id}.linx.clusters.tsv")
    links_file = os.path.join(sample_dir, "linx", "somatic_annotations", f"{tumor_id}.linx.links.tsv")
    drivers_file = os.path.join(sample_dir, "linx", "somatic_annotations", f"{tumor_id}.linx.drivers.tsv")
    svs_file = os.path.join(sample_dir, "linx", "somatic_annotations", f"{tumor_id}.linx.svs.tsv")
    breakend_file = os.path.join(sample_dir, "linx", "somatic_annotations", f"{tumor_id}.linx.breakend.tsv")

    if is_normal or not os.path.exists(clusters_file):
        metrics["ECDNA_STATUS"] = "NA"
        metrics["ECDNA_CLUSTER_COUNT"] = "NA"
        metrics["ECDNA_AMPLICON_GENES"] = "NA"
        metrics["ECDNA_MAX_COPY_NUMBER"] = "NA"
        metrics["ECDNA_CLOSED_SEG_LENGTH_MB"] = "NA"
    else:
        # Default values if no DM clusters are found (LINX evaluated, 0 ecDNA clusters)
        metrics["ECDNA_STATUS"] = "No"
        metrics["ECDNA_CLUSTER_COUNT"] = "0"
        metrics["ECDNA_AMPLICON_GENES"] = "None"
        metrics["ECDNA_MAX_COPY_NUMBER"] = "NA"
        metrics["ECDNA_CLOSED_SEG_LENGTH_MB"] = "NA"
        
        if os.path.exists(clusters_file):
            try:
                # Find all Double Minute (DM) or ecDNA-linked clusters
                dm_resolved_cluster_ids = set()
                dm_cluster_ids = set()
                
                # 1. Read clusters file
                with open(clusters_file, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        r_type = row.get("resolvedType", "").strip().upper()
                        if r_type == "DM" or r_type == "DOUBLE_MINUTE":
                            c_id = row.get("clusterId", "").strip()
                            if c_id:
                                dm_resolved_cluster_ids.add(c_id)
                                dm_cluster_ids.add(c_id)
                                
                # 2. Read links file to capture complex clusters containing ecDNA segments
                if os.path.exists(links_file):
                    with open(links_file, 'r', encoding='utf-8') as f:
                        reader = csv.DictReader(f, delimiter='\t')
                        for row in reader:
                            is_ec_link = row.get("ecDna", "").strip().lower() == "true"
                            if is_ec_link:
                                c_id = row.get("clusterId", "").strip()
                                if c_id:
                                    dm_cluster_ids.add(c_id)
                
                dm_count = len(dm_cluster_ids)
                if dm_count > 0:
                    metrics["ECDNA_STATUS"] = "Yes"
                    metrics["ECDNA_CLUSTER_COUNT"] = str(dm_count)
                    
                    # Compute segment lengths and copy number from links
                    max_cn = 0.0
                    total_closed_len = 0.0
                    if os.path.exists(links_file):
                        with open(links_file, 'r', encoding='utf-8') as f:
                            reader = csv.DictReader(f, delimiter='\t')
                            for row in reader:
                                c_id = row.get("clusterId", "").strip()
                                is_ec_link = row.get("ecDna", "").strip().lower() == "true"
                                is_dm_cluster = c_id in dm_resolved_cluster_ids
                                
                                # Only count the link if it is flagged as ecDna, or if the cluster is a resolved double minute
                                if c_id in dm_cluster_ids and (is_ec_link or is_dm_cluster):
                                    try:
                                        length = float(row.get("length", 0.0))
                                        total_closed_len += length
                                    except ValueError:
                                        pass
                                    try:
                                        cn = float(row.get("junctionCopyNumber", 0.0))
                                        if cn > max_cn:
                                            max_cn = cn
                                    except ValueError:
                                        pass
                    
                    metrics["ECDNA_MAX_COPY_NUMBER"] = f"{max_cn:.2f}"
                    metrics["ECDNA_CLOSED_SEG_LENGTH_MB"] = f"{(total_closed_len / 1000000.0):.3f}"
                    
                    # Identify amplified genes associated with the double minute cluster
                    all_genes = set()
                    
                    # 1. Driver genes from drivers file
                    if os.path.exists(drivers_file):
                        with open(drivers_file, 'r', encoding='utf-8') as f:
                            reader = csv.DictReader(f, delimiter='\t')
                            for row in reader:
                                c_id = row.get("clusterId", "").strip()
                                e_type = row.get("eventType", "").strip().upper()
                                if c_id in dm_cluster_ids and e_type == "GAIN":
                                    gene = row.get("gene", "").strip()
                                    if gene:
                                        all_genes.add(gene)
                    
                    # 2. Genes directly disrupted by DM SVs
                    if os.path.exists(svs_file) and os.path.exists(breakend_file):
                        # Map svId to clusterId
                        sv_to_cluster = {}
                        with open(svs_file, 'r', encoding='utf-8') as f:
                            reader = csv.DictReader(f, delimiter='\t')
                            for row in reader:
                                s_id = row.get("svId", "").strip()
                                c_id = row.get("clusterId", "").strip()
                                if s_id and c_id:
                                    sv_to_cluster[s_id] = c_id
                        
                        # Find overlapping genes from breakends
                        with open(breakend_file, 'r', encoding='utf-8') as f:
                            reader = csv.DictReader(f, delimiter='\t')
                            for row in reader:
                                s_id = row.get("SvId", "").strip()
                                if s_id in sv_to_cluster and sv_to_cluster[s_id] in dm_cluster_ids:
                                    gene = row.get("Gene", "").strip()
                                    if gene and gene.upper() != "NA":
                                        all_genes.add(gene)
                    
                    if all_genes:
                        metrics["ECDNA_AMPLICON_GENES"] = "; ".join(sorted(all_genes))
                    else:
                        metrics["ECDNA_AMPLICON_GENES"] = "None"
            except Exception as e:
                print(f"Warning: Failed to parse ecDNA metrics from clusters/links: {e}", file=sys.stderr)

    # 10. CUPPA Cancer Type Predictions
    if is_normal:
        metrics["CUPPA_PREDICTED_CANCER_TYPE"] = "NA"
        metrics["CUPPA_LIKELIHOOD"] = "NA"
        metrics["CUPPA_TOP_PREDICTIONS"] = "NA"
    else:
        if not cuppa_cancer_type:
            cuppa_tsv = find_cuppa_tsv(tumor_id, sample_dir)
            if cuppa_tsv:
                c_info = parse_cuppa_pred_summ(cuppa_tsv)
                cuppa_cancer_type = c_info.get("cancer_type")
                cuppa_likelihood = c_info.get("likelihood")
                cuppa_predictions = c_info.get("predictions", [])

        if cuppa_cancer_type:
            metrics["CUPPA_PREDICTED_CANCER_TYPE"] = str(cuppa_cancer_type).strip()
            metrics["CUPPA_LIKELIHOOD"] = format_numeric(cuppa_likelihood, 4)

            top_preds_list = []
            if isinstance(cuppa_predictions, list):
                for p in cuppa_predictions[:4]:
                    if isinstance(p, dict):
                        ct = p.get("cancerType")
                        lh = p.get("likelihood")
                        if ct:
                            if lh is not None and str(lh).upper() not in ["NA", "NAN", ""]:
                                try:
                                    lh_pct = f"{float(lh)*100.0:.1f}%"
                                    top_preds_list.append(f"{ct} ({lh_pct})")
                                except ValueError:
                                    top_preds_list.append(str(ct))
                            else:
                                top_preds_list.append(str(ct))

            metrics["CUPPA_TOP_PREDICTIONS"] = "; ".join(top_preds_list) if top_preds_list else str(cuppa_cancer_type)
        else:
            metrics["CUPPA_PREDICTED_CANCER_TYPE"] = "NA"
            metrics["CUPPA_LIKELIHOOD"] = "NA"
            metrics["CUPPA_TOP_PREDICTIONS"] = "NA"

    # Write output to TSV
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    with open(out_tsv, 'w') as fout:
        for k, v in metrics.items():
            sanitized_k = sanitize_attribute_name(k)
            fout.write(f"{sanitized_k}\t{v}\n")

if __name__ == "__main__":
    main()
