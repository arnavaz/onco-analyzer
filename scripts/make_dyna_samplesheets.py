#!/usr/bin/env python3
"""
Generate per-patient samplesheets for Dynacare MultiMMR Panel v4 runs.

Expects paired tumour and normal FASTQs named:
  Tumour: T{num}-*_L00{lane}_R[12]_001.fastq.gz
  Normal: N{num}-*_L00{lane}_R[12]_001.fastq.gz

Creates:
  <out>/<prefix><num>/samplesheet.csv (e.g. <out>/DYNA-P1/samplesheet.csv)
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Generate per-patient Dynacare Panel v4 samplesheets")
    parser.add_argument(
        "--fastq-dir",
        required=True,
        help="Directory containing Dynacare FASTQ files"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Root results directory where per-patient subdirectories will be created"
    )
    parser.add_argument(
        "--patient-prefix",
        default="DYNA-P",
        help="Prefix for subject/group IDs and directories (default: DYNA-P)"
    )
    return parser.parse_args()


def normalize_cluster_path(p: str) -> str:
    """Ensure path uses /cluster/projects/... even if executed on macOS where /cluster is a symlink."""
    p_str = str(p)
    mac_cluster_prefix = "/Users/jbruce/Desktop/cluster"
    if p_str.startswith(mac_cluster_prefix):
        return "/cluster" + p_str[len(mac_cluster_prefix):]
    return p_str


def main():
    args = parse_args()
    fastq_dir = Path(args.fastq_dir)
    out_root = Path(args.out)

    if not fastq_dir.is_dir():
        print(f"ERROR: FASTQ directory not found: {fastq_dir}", file=sys.stderr)
        return 1

    # Regexp to extract Dynacare FASTQ metadata:
    # e.g. T1-19-00888-TE_S1_L001_R1_001.fastq.gz
    # Group 1: 'T' or 'N'
    # Group 2: Patient number (e.g. '1', '2')
    # Group 3: Prefix / Sample ID before _L00X (e.g. 'T1-19-00888-TE_S1')
    # Group 4: Lane (e.g. '001')
    # Group 5: Read ('1' or '2')
    pattern = re.compile(
        r"^(?P<type_char>[TN])(?P<pat_num>\d+)-(?P<rest>.+?)_L(?P<lane>\d{3})_R(?P<read>[12])(?:_\d+)?\.(?:fastq|fq)\.gz$"
    )

    # Structure: patients[pat_num][sample_type][sample_id][lane][read] = path
    patients = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))

    all_files = sorted(fastq_dir.glob("*.f*q.gz"))
    if not all_files:
        print(f"ERROR: No FASTQ files found in {fastq_dir}", file=sys.stderr)
        return 1

    matched_count = 0
    for f in all_files:
        m = pattern.match(f.name)
        if not m:
            continue
        type_char = m.group("type_char")
        pat_num = int(m.group("pat_num"))
        lane = m.group("lane")
        read = m.group("read")

        sample_type = "tumor" if type_char == "T" else "normal"
        prefix = f.name.split(f"_L{lane}_")[0]

        norm_path = normalize_cluster_path(str(f))
        patients[pat_num][sample_type][prefix][lane][read] = norm_path
        matched_count += 1

    print(f"Found {matched_count} matching FASTQ files across {len(patients)} patient(s).")

    fieldnames = [
        "group_id",
        "subject_id",
        "sample_id",
        "sample_type",
        "sequence_type",
        "filetype",
        "info",
        "filepath",
    ]

    for pat_num in sorted(patients.keys()):
        pat_id = f"{args.patient_prefix}{pat_num}"
        pat_data = patients[pat_num]

        pat_dir = out_root / pat_id
        pat_dir.mkdir(parents=True, exist_ok=True)
        sheet_path = pat_dir / "samplesheet.csv"

        rows = []
        for stype in ["tumor", "normal"]:
            for sample_id in sorted(pat_data[stype].keys()):
                for lane in sorted(pat_data[stype][sample_id].keys()):
                    pair = pat_data[stype][sample_id][lane]
                    r1 = pair.get("1")
                    r2 = pair.get("2")
                    if not r1 or not r2:
                        print(f"WARNING: Incomplete pair for {sample_id} lane {lane} in patient {pat_id}", file=sys.stderr)
                        continue

                    rows.append({
                        "group_id": pat_id,
                        "subject_id": pat_id,
                        "sample_id": sample_id,
                        "sample_type": stype,
                        "sequence_type": "dna",
                        "filetype": "fastq",
                        "info": f"library_id:{sample_id};lane:{lane}",
                        "filepath": f"{r1};{r2}",
                    })

        with open(sheet_path, "w", newline="", encoding="utf-8") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        print(f"  {pat_id}: Wrote {len(rows)} rows to {sheet_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
