#!/usr/bin/env python3
import os
import sys
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract tumor fraction comparisons from panel ichorCNA and sWGS results.")
    parser.add_argument("--results-dir", required=True, help="Path to the cohort results directory.")
    parser.add_argument("--ichor-parent-dir", nargs="+", default=[
        "/cluster/projects/pughlab/pughlab/projects/MultiMMR/UHNBB/sWGS/IchorCNA",
        "/cluster/projects/pughlab/pughlab/projects/MultiMMR/OTB/sWGS/IchorCNA"
    ], help="Path(s) to the sWGS ichorCNA parent directories.")
    return parser.parse_args()

def parse_tf_from_file(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if 'tumor fraction' in line.lower() and ':' in line:
                    return line.split(':', 1)[1].strip()
    except Exception:
        pass
    return "NA"

def find_tf_in_dir(directory, patient_id=None):
    """
    Finds a params.txt or metrics.txt file in directory and extracts the TF.
    Prioritizes files containing patient_id and prioritizes params.txt over metrics.txt.
    """
    if not os.path.isdir(directory):
        return "NA"
        
    candidates = []
    
    # 1. Prioritize patient_id + params.txt
    if patient_id:
        for root, dirs, files in os.walk(directory):
            for f in files:
                if patient_id in f and 'params.txt' in f:
                    candidates.append(os.path.join(root, f))
                    
        # 2. Prioritize patient_id + metrics.txt
        for root, dirs, files in os.walk(directory):
            for f in files:
                if patient_id in f and 'metrics.txt' in f:
                    candidates.append(os.path.join(root, f))
                    
    # 3. Fallback to general params.txt
    for root, dirs, files in os.walk(directory):
        for f in files:
            if 'params.txt' in f:
                candidates.append(os.path.join(root, f))
                
    # 4. Fallback to general metrics.txt
    for root, dirs, files in os.walk(directory):
        for f in files:
            if 'metrics.txt' in f:
                candidates.append(os.path.join(root, f))
                
    # Evaluate candidates in order and return the first valid TF
    for path in candidates:
        tf = parse_tf_from_file(path)
        if tf != "NA":
            return tf
            
    return "NA"

def main():
    args = parse_args()
    results_dir = os.path.abspath(args.results_dir)
    ichor_parents = [os.path.abspath(d) for d in args.ichor_parent_dir]

    if not os.path.exists(results_dir):
        print(f"Error: Results directory does not exist: {results_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning results directory: {results_dir}")
    print(f"Using sWGS parent directories: {ichor_parents}")

    processed_count = 0
    for patient_name in os.listdir(results_dir):
        patient_path = os.path.join(results_dir, patient_name)
        if not os.path.isdir(patient_path):
            continue
            
        ichor_bases = [
            os.path.join(patient_path, 'ichorCNA_t_vs_n'),
            os.path.join(patient_path, 'ichorCNA_plasma_vs_n'),
            os.path.join(patient_path, 'ichorCNA_tumor_only'),
            os.path.join(patient_path, 'ichorCNA_plasma_only'),
            os.path.join(patient_path, 'ichorCNA'),
        ]
        
        seen_samples = set()
        for ichor_base in ichor_bases:
            if not os.path.exists(ichor_base):
                continue
                
            # Find all sample subdirectories under this ichorCNA directory
            for sample_name in os.listdir(ichor_base):
                if sample_name in seen_samples:
                    continue
                sample_path = os.path.join(ichor_base, sample_name)
                if not os.path.isdir(sample_path):
                    continue
                # Skip normal control directories
                if sample_name.endswith('-N') or '-N-' in sample_name or 'normal' in sample_name.lower() or 'control' in sample_name.lower():
                    continue
                seen_samples.add(sample_name)
                
            # Extract panel TF
            panel_tf = find_tf_in_dir(sample_path)
                    
            # Resolve original/provided sWGS folder across parent dirs
            provided_tf = "NA"
            clean_name = sample_name
            for prefix in ['DYNA-', 'DYNA_']:
                if clean_name.startswith(prefix):
                    clean_name = clean_name[len(prefix):]
            m = re.search(r'[A-Z]+-[0-9]+', clean_name)
            patient_id = m.group(0) if m else clean_name
            
            s_dir = ""
            for parent_dir in ichor_parents:
                if not os.path.exists(parent_dir):
                    continue
                temp_dir = os.path.join(parent_dir, patient_id, f"{patient_id}-T")
                if not os.path.exists(temp_dir):
                    temp_dir = os.path.join(parent_dir, patient_id)
                if not os.path.exists(temp_dir):
                    candidates = []
                    try:
                        for d1 in os.listdir(parent_dir):
                            p1 = os.path.join(parent_dir, d1)
                            if os.path.isdir(p1):
                                if patient_id in d1 and d1.endswith('-T'):
                                    candidates.append(p1)
                                for d2 in os.listdir(p1):
                                    p2 = os.path.join(p1, d2)
                                    if os.path.isdir(p2):
                                        if patient_id in d2 and d2.endswith('-T'):
                                            candidates.append(p2)
                    except Exception:
                        pass
                    if candidates:
                        temp_dir = candidates[0]
                
                if os.path.exists(temp_dir):
                    s_dir = temp_dir
                    break
                    
            if os.path.exists(s_dir):
                provided_tf = find_tf_in_dir(s_dir, patient_id=patient_id)
                
            # Extract HMF pipeline tumor fraction (PURPLE purity)
            hmf_tf = "NA"
            purple_dir = os.path.join(patient_path, sample_name, "purple")
            if os.path.isdir(purple_dir):
                for root, dirs, files in os.walk(purple_dir):
                    for f in files:
                        if f.endswith('.purple.purity.tsv'):
                            try:
                                with open(os.path.join(root, f), 'r', encoding='utf-8') as purity_f:
                                    lines = purity_f.readlines()
                                    if len(lines) >= 2:
                                        hmf_tf = lines[1].split('\t')[0].strip()
                            except Exception:
                                pass
                            break
                        
            comparison_file = os.path.join(ichor_base, 'tumor_fraction_comparison.tsv')
            try:
                with open(comparison_file, 'w', encoding='utf-8') as out_f:
                    out_f.write("sample_id\tpanel_tumor_fraction\tprovided_tumor_fraction\thmf_tumor_fraction\n")
                    out_f.write(f"{sample_name}\t{panel_tf}\t{provided_tf}\t{hmf_tf}\n")
                print(f"  [{patient_name}] Generated comparison TSV: panel={panel_tf} vs provided={provided_tf} vs hmf={hmf_tf}")
                legacy_dir = os.path.join(patient_path, 'ichorCNA')
                if os.path.isdir(legacy_dir) and legacy_dir != ichor_base:
                    try:
                        import shutil
                        shutil.copy2(comparison_file, os.path.join(legacy_dir, 'tumor_fraction_comparison.tsv'))
                    except Exception:
                        pass
                processed_count += 1
            except Exception as e:
                print(f"Warning: Failed to write comparison file {comparison_file}: {e}", file=sys.stderr)
                
    print(f"\nDone! Successfully processed {processed_count} sample(s).")

if __name__ == "__main__":
    main()
