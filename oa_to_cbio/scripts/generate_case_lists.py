import os
import sys

if len(sys.argv) < 3:
    print("Usage: python3 generate_case_lists.py <import_dir> <study_id>")
    sys.exit(1)

import_dir = sys.argv[1]
study_id = sys.argv[2]

case_lists_dir = os.path.join(import_dir, "case_lists")
os.makedirs(case_lists_dir, exist_ok=True)

def write_case_list(suffix, name, desc, samples, category=None):
    if not samples:
        return
    stable_id = f"{study_id}_{suffix}"
    out_file = os.path.join(case_lists_dir, f"cases_{suffix}.txt")
    with open(out_file, 'w') as f:
        f.write(f"cancer_study_identifier: {study_id}\n")
        f.write(f"stable_id: {stable_id}\n")
        f.write(f"case_list_name: {name}\n")
        f.write(f"case_list_description: {desc} ({len(samples)} samples)\n")
        ids_str = '\t'.join(samples)
        f.write(f"case_list_ids: {ids_str}\n")
        if category:
            f.write(f"case_list_category: {category}\n")

# 1. Parse samples from each file type
samples_mut = set()
mut_file = os.path.join(import_dir, "data_mutations.txt")
if os.path.exists(mut_file):
    with open(mut_file, 'r') as f:
        header = None
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip('\n').split('\t')
            if not header:
                header = parts
                tsb_idx = header.index("Tumor_Sample_Barcode") if "Tumor_Sample_Barcode" in header else -1
            elif tsb_idx >= 0 and len(parts) > tsb_idx:
                samples_mut.add(parts[tsb_idx])

samples_cna = set()
cna_file = os.path.join(import_dir, "data_CNA.txt")
if os.path.exists(cna_file):
    with open(cna_file, 'r') as f:
        for line in f:
            if line.startswith('Hugo_Symbol'):
                parts = line.strip('\n').split('\t')
                samples_cna.update(parts[1:])
                break

samples_sv = set()
sv_file = os.path.join(import_dir, "data_structural_variants.txt")
if os.path.exists(sv_file):
    with open(sv_file, 'r') as f:
        header = None
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip('\n').split('\t')
            if not header:
                header = parts
                tsb_idx = header.index("Sample_Id") if "Sample_Id" in header else -1
            elif tsb_idx >= 0 and len(parts) > tsb_idx:
                samples_sv.add(parts[tsb_idx])

samples_rna = set()
rna_file = os.path.join(import_dir, "data_mrna_seq_v2_rsem.txt")
if os.path.exists(rna_file):
    with open(rna_file, 'r') as f:
        for line in f:
            if line.startswith('Hugo_Symbol'):
                parts = line.strip('\n').split('\t')
                if "Entrez_Gene_Id" in parts:
                    samples_rna.update(parts[2:]) # Skip Hugo_Symbol and Entrez_Gene_Id
                else:
                    samples_rna.update(parts[1:]) # Skip Hugo_Symbol only
                break

# Intersections
samples_cnaseq = samples_mut.intersection(samples_cna)
samples_3way = samples_cnaseq.intersection(samples_rna)

# Write Case Lists
write_case_list(
    suffix="sequenced",
    name="Sequenced Tumors",
    desc="All sequenced tumors",
    samples=sorted(list(samples_mut)),
    category="all_cases_with_mutation_data"
)

write_case_list(
    suffix="cna",
    name="CNA Tumors",
    desc="All tumors with CNA data",
    samples=sorted(list(samples_cna)),
    category="all_cases_with_cna_data"
)

write_case_list(
    suffix="cnaseq",
    name="Sequenced and CNA Tumors",
    desc="All tumors with mutation and CNA data",
    samples=sorted(list(samples_cnaseq)),
    category="all_cases_with_mutation_and_cna_data"
)

write_case_list(
    suffix="rna_seq_v2_mrna",
    name="mRNA (RNA-SeqV2) Tumors",
    desc="All tumors with mRNA expression data",
    samples=sorted(list(samples_rna)),
    category="all_cases_with_mrna_rnaseq_data"
)

write_case_list(
    suffix="sv",
    name="SV Tumors",
    desc="All tumors with Structural Variant data",
    samples=sorted(list(samples_sv)),
    category="all_cases_with_sv_data"
)

write_case_list(
    suffix="3way_complete",
    name="Complete Tumors",
    desc="All tumors with mutation, CNA, and mRNA data",
    samples=sorted(list(samples_3way)),
    category="all_cases_with_mutation_and_cna_and_mrna_data"
)

print("Generated case lists in", case_lists_dir)
