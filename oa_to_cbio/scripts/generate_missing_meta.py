#!/usr/bin/env python3
import os
import sys

if len(sys.argv) < 3:
    print("Usage: python3 generate_missing_meta.py <import_dir> <study_id>")
    sys.exit(1)

import_dir = sys.argv[1]
study_id = sys.argv[2]

# Standard cBioPortal metadata templates based on filename patterns
METADATA_TEMPLATES = {
    # Timelines (e.g. data_clinical_timeline_specimen.txt)
    "clinical_timeline_": {
        "genetic_alteration_type": "CLINICAL",
        "datatype": "TIMELINE",
    },
    # Clinical Data
    "clinical_patient": {
        "genetic_alteration_type": "CLINICAL",
        "datatype": "PATIENT_ATTRIBUTES",
    },
    "clinical_sample": {
        "genetic_alteration_type": "CLINICAL",
        "datatype": "SAMPLE_ATTRIBUTES",
    },
    # Mutations
    "mutations": {
        "genetic_alteration_type": "MUTATION_EXTENDED",
        "datatype": "MAF",
        "stable_id": "mutations",
        "show_profile_in_analysis_tab": "true",
        "profile_name": "Mutations",
        "profile_description": "Mutation data.",
    },
    # Copy Number Alteration
    "CNA": {
        "genetic_alteration_type": "COPY_NUMBER_ALTERATION",
        "datatype": "DISCRETE",
        "stable_id": "cna",
        "show_profile_in_analysis_tab": "true",
        "profile_name": "Putative copy-number alterations",
        "profile_description": "Putative copy-number alterations. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.",
    },
    "segment": {
        "genetic_alteration_type": "COPY_NUMBER_ALTERATION",
        "datatype": "SEG",
        "reference_genome_id": "hg38",
        "description": "Copy number segment data.",
    },
    # Structural Variants
    "structural_variants": {
        "genetic_alteration_type": "STRUCTURAL_VARIANT",
        "datatype": "SV",
        "stable_id": "structural_variants",
        "show_profile_in_analysis_tab": "true",
        "profile_name": "Structural variants",
        "profile_description": "Structural variant data.",
    },
    # Expression / mRNA
    "mrna_seq_v2_rsem": {
        "genetic_alteration_type": "MRNA_EXPRESSION",
        "datatype": "CONTINUOUS",
        "stable_id": "rna_seq_v2_mrna",
        "show_profile_in_analysis_tab": "true",
        "profile_name": "mRNA expression",
        "profile_description": "mRNA expression data.",
    },
    # Resource Data (WISP / Liquid Biopsy / Custom Tabs)
    "resource_definition": {
        "resource_type": "DEFINITION",
    },
    "resource_patient": {
        "resource_type": "PATIENT",
    },
    "resource_sample": {
        "resource_type": "SAMPLE",
    },
    "resource_study": {
        "resource_type": "STUDY",
    },
}

def sanitize_timeline_file(filepath):
    print(f"  Sanitizing timeline file: {filepath}")
    headers = None
    rows = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('#'):
                rows.append((True, line))
                continue
            parts = line.split('\t')
            parts = [p.strip('"') for p in parts]
            if headers is None:
                headers = parts
            rows.append((False, parts))
            
    if not headers:
        return
        
    start_date_idx = None
    try:
        start_date_idx = headers.index("START_DATE")
    except ValueError:
        pass
        
    with open(filepath, 'w', encoding='utf-8') as f:
        for is_comment, parts in rows:
            if is_comment:
                f.write(parts + "\n")
                continue
            if parts == headers:
                f.write("\t".join(parts) + "\n")
                continue
            
            # Skip the row entirely if START_DATE is missing or NA
            if start_date_idx is not None and start_date_idx < len(parts):
                sd_val = parts[start_date_idx]
                if not sd_val or sd_val.upper() == "NA":
                    continue
            
            new_parts = []
            for idx, val in enumerate(parts):
                col = headers[idx] if idx < len(headers) else ""
                
                # Rule 2: No "NA" values in STOP_DATE; it needs to be a blank
                if col == "STOP_DATE":
                    if not val or val.upper() == "NA":
                        new_parts.append("")
                        continue
                
                # Rule 3: No floating values for DATE elements
                if col.endswith("_DATE"):
                    if not val or val.upper() == "NA":
                        new_parts.append("" if col == "STOP_DATE" else "NA")
                        continue
                    try:
                        val = str(int(round(float(val))))
                    except ValueError:
                        pass
                
                new_parts.append(val)
            f.write("\t".join(new_parts) + "\n")


if not os.path.exists(import_dir):
    print(f"Error: Import directory {import_dir} does not exist.")
    sys.exit(1)

print(f"Scanning {import_dir} for missing metadata files...")

for filename in os.listdir(import_dir):
    if filename.startswith("data_") and filename.endswith(".txt"):
        base_name = filename[5:-4]  # Strip "data_" and ".txt"
        filepath = os.path.join(import_dir, filename)
        
        if "clinical_timeline_" in base_name:
            sanitize_timeline_file(filepath)
            
        meta_filename = f"meta_{base_name}.txt"
        meta_filepath = os.path.join(import_dir, meta_filename)
        
        if not os.path.exists(meta_filepath):
            # Find matching template based on substring match of filename
            template = None
            for key, t in METADATA_TEMPLATES.items():
                if key in base_name:
                    template = t
                    break
            
            if template:
                print(f"  Missing metadata file for '{filename}' -> generating '{meta_filename}'")
                with open(meta_filepath, 'w') as f:
                    f.write(f"cancer_study_identifier: {study_id}\n")
                    for k, v in template.items():
                        f.write(f"{k}: {v}\n")
                    f.write(f"data_filename: {filename}\n")
            else:
                print(f"  Warning: Could not infer metadata fields for '{filename}'. Skip generating.")
