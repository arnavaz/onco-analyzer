#!/usr/bin/env python3
"""Build per-sample nf-core/oncoanalyser samplesheets for RNA-Seq only samples.

Reads a list of FASTQ paths (or a FASTQ directory) and a metadata TSV mapping
samples/tissues to patients. For each RNA-Seq sample, writes a valid Oncoanalyser
samplesheet to:
    <out_dir>/<sample_name>/samplesheet.csv
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build oncoanalyser samplesheets for RNA-Seq only runs"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--fastq-list",
        help="Path to a text file containing one FASTQ filepath per line",
    )
    group.add_argument(
        "--fastq-dir",
        help="Directory containing FASTQ files",
    )
    parser.add_argument(
        "--metadata-tsv",
        required=True,
        help="Path to sample metadata TSV (e.g. MMR_RNAseq.tsv)",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output directory where per-sample folders with samplesheet.csv will be created",
    )
    parser.add_argument(
        "--sample-type",
        default="tumor",
        choices=["tumor", "normal"],
        help="sample_type for Oncoanalyser samplesheet (default: tumor; oncoanalyser requires tumor for RNA quantification)",
    )
    return parser.parse_args()


def collect_fastq_pairs(fastq_list_path=None, fastq_dir=None):
    if fastq_list_path:
        with open(fastq_list_path, "r", encoding="utf-8") as f:
            fastqs = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]
    elif fastq_dir:
        p = Path(fastq_dir)
        fastqs = [str(f) for f in p.glob("**/*") if f.is_file() and f.name.endswith((".fastq.gz", ".fq.gz"))]
    else:
        fastqs = []

    # Identify R1s and find their matching R2s
    r1_files = sorted([f for f in fastqs if "_R1" in Path(f).name])
    pairs = {}
    for r1 in r1_files:
        r1_name = Path(r1).name
        # Match standard R1/R2 replacements
        if "_R1.fastq.gz" in r1:
            r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")
            prefix = r1_name.replace("_R1.fastq.gz", "")
        elif "_R1.fq.gz" in r1:
            r2 = r1.replace("_R1.fq.gz", "_R2.fq.gz")
            prefix = r1_name.replace("_R1.fq.gz", "")
        elif "_R1_001.fastq.gz" in r1:
            r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
            prefix = r1_name.replace("_R1_001.fastq.gz", "")
        else:
            r2 = re.sub(r"_R1\b", "_R2", r1)
            prefix = re.sub(r"_R1\b.*", "", r1_name)

        pairs[prefix] = (r1, r2)
    return pairs


def generate_samplesheets(pairs, metadata_tsv_path, out_dir, default_sample_type="tumor"):
    df = pd.read_csv(metadata_tsv_path, sep="\t")

    # Filter out rows explicitly annotated as having no RNA-seq
    if "RNA-seq available?" in df.columns:
        df_rna = df[df["RNA-seq available?"].astype(str).str.strip().str.lower() != "no"].copy()
    else:
        df_rna = df.copy()

    records = []
    unmapped = []

    for prefix, (r1, r2) in sorted(pairs.items()):
        norm_p = prefix.replace("_", "-")
        matched_rows = []

        for _, row in df_rna.iterrows():
            sid = str(row.get("RNAseq_Sample_ID", "")).strip()
            tid = str(row.get("Tissue_ID", "")).strip()
            pid = str(row.get("Patient_ID", "")).strip()

            norm_sid = sid.replace("_", "-")
            clean_sid = norm_sid.replace("-R1", "").replace("-R2", "")

            # Matching criteria
            if (
                prefix == sid
                or norm_p == norm_sid
                or norm_p == clean_sid
                or prefix in sid
                or sid in prefix
                or prefix == tid
                or (prefix == "JG" and "JG" in sid)
            ):
                matched_rows.append(row)

        if not matched_rows:
            unmapped.append((prefix, r1, r2))
            continue

        # When multiple metadata rows match (e.g. multiple tumor biopsies or normal/tumor pairs)
        if len(matched_rows) > 1:
            tumor_rows = [r for r in matched_rows if str(r.get("PathType", "")).strip().lower() == "tumor"]
            selected_row = tumor_rows[0] if tumor_rows else matched_rows[0]
        else:
            selected_row = matched_rows[0]

        patient_id = str(selected_row.get("Patient_ID", "")).strip()
        tissue_id = str(selected_row.get("Tissue_ID", "")).strip()

        # If tissue_id ends in a sub-aliquot letter like A/B for a single shared sequencing run, normalize
        sample_name = tissue_id
        if prefix == "19_16409" and sample_name.endswith(("A", "B")):
            sample_name = sample_name[:-1]

        records.append({
            "sample_name": sample_name,
            "patient_id": patient_id,
            "prefix": prefix,
            "r1": r1,
            "r2": r2,
        })

    os.makedirs(out_dir, exist_ok=True)
    created_count = 0

    print(f"Generating samplesheets for {len(records)} RNA samples in: {out_dir}")
    for rec in records:
        sample_dir = os.path.join(out_dir, rec["sample_name"])
        os.makedirs(sample_dir, exist_ok=True)
        sheet_path = os.path.join(sample_dir, "samplesheet.csv")

        # Oncoanalyser schema fields:
        # group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath
        with open(sheet_path, "w", encoding="utf-8") as f:
            f.write("group_id,subject_id,sample_id,sample_type,sequence_type,filetype,info,filepath\n")
            f.write(
                f"{rec['sample_name']},{rec['patient_id']},{rec['sample_name']},"
                f"{default_sample_type},rna,fastq,library_id:1;lane:1,"
                f"{rec['r1']};{rec['r2']}\n"
            )
        created_count += 1
        print(f"  [{rec['sample_name']}] Patient={rec['patient_id']} (FASTQ: {rec['prefix']}) -> {sheet_path}")

    if unmapped:
        print(f"\nWARNING: {len(unmapped)} FASTQ prefix(es) could not be mapped to metadata:")
        for prefix, r1, _ in unmapped:
            print(f"  - {prefix} ({r1})")

    print(f"\nDone. Successfully generated {created_count} samplesheet(s).")


def main():
    args = parse_args()
    pairs = collect_fastq_pairs(args.fastq_list, args.fastq_dir)
    if not pairs:
        print("ERROR: No FASTQ pairs found.", file=sys.stderr)
        sys.exit(1)
    generate_samplesheets(pairs, args.metadata_tsv, args.out, args.sample_type)


if __name__ == "__main__":
    main()
