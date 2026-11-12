#!/usr/bin/env python3
"""Build per-patient nf-core/oncoanalyser samplesheets from a directory of lane-level FASTQs.

For each patient (subject_id) found in the FASTQs, one samplesheet is written to:
  <out_dir>/<subject_id>/samplesheet.csv

Expected FASTQ naming convention:
  <prefix>_L00X_R1_001.fastq.gz
  <prefix>_L00X_R2_001.fastq.gz

If a mapping TSV is provided (e.g. FULL_PANCAN-H_CHANGE_injection_sheet.tsv),
the script will infer sample_id, sample_type, and sequence_type from that file.

For FASTQ inputs, nf-core/oncoanalyser expects sequencing metadata in the `info`
column, e.g. `library_id:S1;lane:001`.

The `library_id` must be unique per biological sample + library combination. This
script derives it from the FASTQ prefix so multiple libraries for the same sample
can safely reuse lane numbers without colliding.
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

FASTQ_RE = re.compile(
    r"^(?P<prefix>.+)_L(?P<lane>\d{3})_R(?P<read>[12])(?:_\d+)?\.(?:fastq|fq)\.gz$"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create oncoanalyser samplesheet from FASTQs or BAMs")
    p.add_argument(
        "--fastq-dir",
        default=None,
        help="Directory containing FASTQ files (for single/combined directory search)"
    )
    p.add_argument(
        "--tumor-fastq-dir",
        default=None,
        help="Directory containing tumor FASTQ files (for paired directory search)"
    )
    p.add_argument(
        "--normal-fastq-dir",
        default=None,
        help="Directory containing normal/reference FASTQ files (for paired directory search)"
    )
    p.add_argument(
        "--bam-dir",
        default=None,
        help="Directory containing BAM/CRAM files (for single/combined directory search)"
    )
    p.add_argument(
        "--tumor-bam-dir",
        default=None,
        help="Directory containing tumor BAM/CRAM files (for paired directory search)"
    )
    p.add_argument(
        "--normal-bam-dir",
        default=None,
        help="Directory containing normal/reference BAM/CRAM files (for paired directory search)"
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output directory (cohort root). Each patient gets a subdirectory: <out>/<subject_id>/samplesheet.csv",
    )
    p.add_argument(
        "--mapping-tsv",
        default=None,
        help="Optional injection/mapping TSV with columns including Fastq File Name and Sample ID",
    )
    p.add_argument(
        "--default-sequence-type",
        default="dna",
        choices=["dna", "rna"],
        help="sequence_type used when mapping TSV is not provided",
    )
    p.add_argument(
        "--default-sample-type",
        default="tumor",
        choices=["tumor", "normal"],
        help="sample_type used when mapping TSV is not provided",
    )
    p.add_argument(
        "--emit-normal-rna-sheet",
        action="store_true",
        help="Also write per-patient normal RNA FASTQ sheet as <out>/<subject_id>/samplesheet.normal_rna.csv",
    )
    p.add_argument(
        "--config-tsv",
        default=None,
        help="Configuration TSV with explicit filepaths (BAM/CRAM/FASTQ) and sample annotations (e.g. redux_config.tsv)",
    )
    p.add_argument(
        "--oncoanalyser-outdir",
        default=None,
        help="Root directory where oncoanalyser output is written (used to reference Step 1 purple/amber/normal BAM in samplesheet.purity.csv). Defaults to parent of --out.",
    )
    p.add_argument(
        "--use-existing-redux-bams",
        default=None,
        help="Path to previous oncoanalyser results directory. When specified, searches for existing *.redux.bam files under <dir>/<subject_id>/*/alignments/*/<sample_id>.redux.bam and reuses them with filetype 'bam_redux' and info 'generate_redux_tsvs_only' (skipping BWA-MEM2 alignment).",
    )
    p.add_argument(
        "--generate-redux-tsvs-only",
        action="store_true",
        help="When supplying BAM inputs, set filetype to 'bam_redux' and add 'generate_redux_tsvs_only' to info field so oncoanalyser v3.1 skips alignment and runs lightweight Redux TSV generation.",
    )
    p.add_argument(
        "--purity-plasma-from-fastq",
        action="store_true",
        help="When using --use-existing-redux-bams, force longitudinal plasma ctDNA in Step 2 purity sheets to use FASTQs instead of existing plasma BAMs.",
    )
    p.add_argument(
        "--purity-reference-dir",
        default=None,
        help="Directory to reference for primary tumor purple/amber outputs in Step 2 purity sheets (defaults to --oncoanalyser-outdir).",
    )
    return p.parse_args()


def norm_header(h: str) -> str:
    return re.sub(r"\s+", " ", h.strip().strip('"').lower())


def is_plasma_meta(meta: Dict[str, str] | None) -> bool:
    if not meta:
        return False
    to = (meta.get("tissue_origin") or "").strip().lower()
    if to in ["pl", "plasma", "cfdna"]:
        return True
    tt = (meta.get("tissue_type") or "").strip().lower()
    if tt in ["pl", "plasma", "cfdna"]:
        return True
    lt = (meta.get("library_type") or "").strip().lower()
    if lt == "pg":
        return True
    return False


def find_redux_bam(base_dir: str | Path, subject_id: str, sample_id: str) -> Path | None:
    """Find existing .redux.bam (or .redux.cram / .bam) for subject and sample.
    Supports both native paths and cluster-to-Mac mount translation.
    """
    bdir_str = str(base_dir).rstrip("/")
    dirs_to_check = [bdir_str]
    if bdir_str.startswith("/cluster/"):
        mac_mount = bdir_str.replace("/cluster", "/Users/jbruce/Desktop/cluster")
        dirs_to_check.append(mac_mount)
    elif bdir_str.startswith("/Users/jbruce/Desktop/cluster"):
        dirs_to_check.append(bdir_str.replace("/Users/jbruce/Desktop/cluster", "/cluster"))

    search_dirs = [d for d in dirs_to_check if Path(d).exists()]
    if not search_dirs:
        search_dirs = [bdir_str]

    for d in search_dirs:
        patterns = [
            f"{d}/{subject_id}/*/alignments/dna/{sample_id}.redux.bam",
            f"{d}/{subject_id}/*/alignments/{sample_id}/{sample_id}.redux.bam",
            f"{d}/{subject_id}/alignments/dna/{sample_id}.redux.bam",
            f"{d}/{subject_id}/alignments/{sample_id}/{sample_id}.redux.bam",
            f"{d}/{subject_id}/*/{sample_id}.redux.bam",
            f"{d}/{subject_id}/{sample_id}.redux.bam",
            f"{d}/{sample_id}.redux.bam",
            f"{d}/{subject_id}/*/alignments/dna/{sample_id}.bam",
            f"{d}/{subject_id}/*/{sample_id}.bam",
            f"{d}/{sample_id}.bam",
        ]
        for pat in patterns:
            matches = glob.glob(pat)
            if matches:
                found_path = matches[0]
                if bdir_str.startswith("/cluster/") and found_path.startswith("/Users/jbruce/Desktop/cluster"):
                    found_path = found_path.replace("/Users/jbruce/Desktop/cluster", "/cluster")
                return Path(found_path)
    return None


def load_mapping(path: Path) -> Tuple[Dict[str, Dict[str, str]], Dict[Tuple[str, str], Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in mapping TSV: {path}")

        header_map = {norm_header(h): h for h in reader.fieldnames}

        def get_col(candidates: List[str], required: bool = True) -> str | None:
            for c in candidates:
                norm_c = norm_header(c)
                if norm_c in header_map:
                    return header_map[norm_c]
            if required:
                raise ValueError(
                    f"Required mapping column not found from candidates {candidates}. Found: {reader.fieldnames}"
                )
            return None

        fastq_col = get_col(["fastq file name", "files"])
        sample_col = get_col(["sample id", "sample_id"])
        tissue_col = get_col(["tissue type (eg. primary tumor, reference)", "tissue type", "tissue_type"])
        lib_col = get_col(["library type (wg:whole genome, wt:whole transcriptome, pg:plasma, ex:exome)", "library type", "library_type"])
        subject_col = get_col(["donor_id", "subject_id", "subject id", "oicr_donor"], required=False)
        tissue_origin_col = get_col(["tissue origin", "tissue_origin"], required=False)

        mapping: Dict[str, Dict[str, str]] = {}
        sample_meta: Dict[Tuple[str, str], Dict[str, str]] = {}
        for row in reader:
            fastq_val = (row.get(fastq_col) or "").strip().strip('"')
            if not fastq_val:
                continue

            fastq_names = [f.strip() for f in fastq_val.split(";")]

            raw_sample_id = (row.get(sample_col) or "").strip().strip('"')
            clean_sample_id = re.sub(r"\s+", "-", raw_sample_id)
            row_data = {
                "sample_id": clean_sample_id,
                "tissue_type": (row.get(tissue_col) or "").strip().strip('"'),
                "library_type": (row.get(lib_col) or "").strip().strip('"'),
                "tissue_origin": (row.get(tissue_origin_col) or "").strip().strip('"') if tissue_origin_col else "",
            }
            if subject_col:
                raw_subject_id = (row.get(subject_col) or "").strip().strip('"')
                row_data["subject_id"] = re.sub(r"\s+", "-", raw_subject_id)
            else:
                row_data["subject_id"] = infer_subject_id(clean_sample_id)

            sub_id = row_data["subject_id"]
            sample_meta[(sub_id, clean_sample_id)] = row_data

            mapping[(sub_id, clean_sample_id)] = row_data
            mapping[f"{sub_id}_{clean_sample_id}"] = row_data
            if clean_sample_id not in mapping:
                mapping[clean_sample_id] = row_data
            elif mapping[clean_sample_id] is not None and mapping[clean_sample_id].get("subject_id") != sub_id:
                mapping[clean_sample_id] = None

            for name in fastq_names:
                if name:
                    mapping[name] = row_data
        return mapping, sample_meta


def infer_subject_id(sample_id: str) -> str:
    # Pugh Lab / OICR pattern: e.g., NL82-MMR1277T01-F -> NL82-MMR1277
    # e.g., DYNA-UHN-1194-FFPE-T-DNA-U_S36 -> DYNA-UHN-1194
    parts = sample_id.split("-")
    if len(parts) >= 2:
        # If the second part has a T\d+ or B\d+ suffix, strip it
        if re.search(r"(?i)[TB]\d+", parts[1]):
            clean_pat = re.sub(r"(?i)[TB]\d+.*$", "", parts[1])
            return f"{parts[0]}-{clean_pat}"
        # If there are more parts and the first part is DYNA
        if parts[0].upper() == "DYNA" and len(parts) >= 3:
            return f"{parts[0]}-{parts[1]}-{parts[2]}"

    # PANCAN-style IDs, e.g. PANCAN-H-0845-01T-01 / PANCAN-H-0845-01N-01
    # -> PANCAN-H-0845
    m = re.match(r"^(.*)-\d{2}[TN]-\d{2}$", sample_id, flags=re.IGNORECASE)
    if m:
        return m.group(1)
    m = re.match(r"^(.*)-\d{2}[TN]$", sample_id, flags=re.IGNORECASE)
    if m:
        return m.group(1)

    m = re.match(r"^(.*)-\d{4}-[TN]$", sample_id)
    if m:
        return m.group(1)
    m = re.match(r"^(.*)-[TN]$", sample_id)
    if m:
        return m.group(1)
    return sample_id


def infer_sample_type(tissue_type: str, sample_id: str, default_value: str, tissue_origin: str = "") -> str:
    if tissue_origin:
        to = tissue_origin.lower().strip()
        if to in ["pl", "plasma", "cfdna"]:
            return "tumor"
        if to in ["blood", "normal", "b"]:
            return "normal"
    if tissue_type:
        tt = tissue_type.lower().strip()
        if tt in ["pl", "plasma", "cfdna"]:
            return "tumor"
        elif tt in ["r", "ref", "reference", "normal", "norm", "n", "blood"]:
            return "normal"
        elif tt in ["t", "tumor", "primary tumor", "primary", "relapse", "metastasis"]:
            return "tumor"
    sid = sample_id.upper()
    # Tissue annotation is the most reliable source when available.
    if sid.endswith("-N") or sid.endswith("-NORMAL"):
        return "normal"
    if sid.endswith("-T") or sid.endswith("-TUMOR"):
        return "tumor"
    return default_value


def infer_sequence_type(library_type: str, default_value: str) -> str:
    lt = library_type.strip().upper()
    if lt == "WT":
        return "rna"
    if lt in {"WG", "PG", "EX"}:
        return "dna"
    return default_value


def infer_tumor_anchor(sample_id: str) -> str:
    """Return a stable tumor anchor used to form one analysis group per tumor specimen."""
    sid = sample_id.strip()
    # Normalize common RNA sample suffixes so DNA/RNA from same tumor co-group.
    sid = re.sub(r"(?i)([-_]RNA)$", "", sid)
    return sid


def main() -> int:
    args = parse_args()

    if args.fastq_dir and (args.tumor_fastq_dir or args.normal_fastq_dir):
        print("ERROR: Cannot combine --fastq-dir with --tumor-fastq-dir or --normal-fastq-dir", file=sys.stderr)
        return 1

    if args.bam_dir and (args.tumor_bam_dir or args.normal_bam_dir):
        print("ERROR: Cannot combine --bam-dir with --tumor-bam-dir or --normal-bam-dir", file=sys.stderr)
        return 1

    if not (args.config_tsv or args.fastq_dir or args.tumor_fastq_dir or args.normal_fastq_dir or args.bam_dir or args.tumor_bam_dir or args.normal_bam_dir or args.use_existing_redux_bams):
        print("ERROR: Must provide at least one input source (--config-tsv, --fastq-dir, --tumor-fastq-dir/--normal-fastq-dir, --bam-dir, --tumor-bam-dir/--normal-bam-dir, or --use-existing-redux-bams)", file=sys.stderr)
        return 1

    if (args.tumor_fastq_dir or args.normal_fastq_dir) and not (args.tumor_fastq_dir and args.normal_fastq_dir):
        print("ERROR: Both --tumor-fastq-dir and --normal-fastq-dir must be provided for paired FASTQ search", file=sys.stderr)
        return 1

    if (args.tumor_bam_dir or args.normal_bam_dir) and not (args.tumor_bam_dir and args.normal_bam_dir):
        print("ERROR: Both --tumor-bam-dir and --normal-bam-dir must be provided for paired BAM search", file=sys.stderr)
        return 1

    mapping: Dict[str, Dict[str, str]] = {}
    sample_meta: Dict[Tuple[str, str], Dict[str, str]] = {}
    if args.mapping_tsv:
        mapping, sample_meta = load_mapping(Path(args.mapping_tsv))

    redux_bam_rows: List[Dict[str, str]] = []
    satisfied_by_redux: set = set()

    if args.use_existing_redux_bams:
        if sample_meta:
            for (sub_id, s_id), s_info in sorted(sample_meta.items()):
                seq_type = infer_sequence_type(s_info.get("library_type", ""), args.default_sequence_type)
                if seq_type != "dna":
                    continue
                is_plasma = is_plasma_meta(s_info)
                if is_plasma and args.purity_plasma_from_fastq:
                    continue
                b = find_redux_bam(args.use_existing_redux_bams, sub_id, s_id)
                if b:
                    stype = infer_sample_type(
                        s_info.get("tissue_type", ""), s_id, args.default_sample_type, tissue_origin=s_info.get("tissue_origin", "")
                    )
                    is_cram = str(b).lower().endswith((".cram", ".redux.cram"))
                    ftype = "cram_redux" if is_cram else "bam_redux"
                    info_p = []
                    if is_plasma:
                        info_p.append("longitudinal_sample")
                    info_p.append("generate_redux_tsvs_only")
                    redux_bam_rows.append({
                        "subject_id": sub_id,
                        "sample_id": s_id,
                        "sample_type": stype,
                        "sequence_type": "dna",
                        "filetype": ftype,
                        "info": ";".join(info_p),
                        "filepath": str(b),
                        "is_plasma": is_plasma,
                    })
                    satisfied_by_redux.add((sub_id, s_id))
                    print(f"INFO: Reusing existing REDUX alignment for {sub_id} {s_id}: {b}")
                else:
                    print(f"WARN: Existing REDUX alignment not found for {sub_id} {s_id} under {args.use_existing_redux_bams}")
        else:
            print("WARN: --use-existing-redux-bams requires --mapping-tsv to resolve sample annotations", file=sys.stderr)

    def scan_dir(d: Path) -> Dict[Tuple[str, str], Dict[str, Path]]:
        reads = defaultdict(dict)
        for fp in sorted(d.glob("*.f*q.gz")):
            m = FASTQ_RE.match(fp.name)
            if not m:
                m = re.match(r"^(?P<prefix>.+?)_(?P<lane>\d+)_(?P<barcode>[ACGT-]+)_R(?P<read>[12])\.(?:fastq|fq)\.gz$", fp.name)
            if not m:
                continue
            prefix = m.group("prefix")
            lane = m.group("lane")
            read = m.group("read")
            reads[(prefix, lane)][read] = fp.resolve()
        return reads

    rows: List[Dict[str, str]] = []
    normal_rna_rows: List[Dict[str, str]] = []
    skipped = 0
    skipped_normal_rna = 0
    skipped_orphan_tumor_rna = 0
    skipped_subjects_no_tumor_dna = 0

    def process_lane_reads(
        reads: Dict[Tuple[str, str], Dict[str, Path]], 
        dir_default_sample_type: str
    ) -> List[Dict[str, str]]:
        nonlocal skipped, skipped_normal_rna
        local_rows = []
        for (prefix, lane), pair in sorted(reads.items()):
            if "1" not in pair or "2" not in pair:
                skipped += 1
                print(
                    f"WARN: Skipping {prefix} lane {lane} because paired R1/R2 is incomplete",
                    file=sys.stderr,
                )
                continue

            r1 = pair["1"]
            r2 = pair["2"]

            sample_id = prefix
            sample_type = dir_default_sample_type
            sequence_type = args.default_sequence_type
            is_plasma = False

            if mapping:
                meta = mapping.get(r1.name)
                if not meta:
                    skipped += 1
                    print(
                        f"WARN: {r1.name} not found in mapping TSV; skipping lane {lane}",
                        file=sys.stderr,
                    )
                    continue
                sample_id = meta["sample_id"] or prefix
                tissue_origin = meta.get("tissue_origin", "")
                sample_type = infer_sample_type(
                    meta["tissue_type"], sample_id, dir_default_sample_type, tissue_origin=tissue_origin
                )
                sequence_type = infer_sequence_type(meta["library_type"], args.default_sequence_type)
                is_plasma = is_plasma_meta(meta)

            subject_id = infer_subject_id(sample_id)
            if mapping and r1.name in mapping and mapping[r1.name].get("subject_id"):
                subject_id = mapping[r1.name]["subject_id"]

            if (subject_id, sample_id) in satisfied_by_redux:
                # Biological sample already has an existing REDUX BAM aligned!
                continue

            # oncoanalyser does not use normal RNA FASTQ inputs in WGTS/TARGETED workflows.
            if sample_type == "normal" and sequence_type == "rna":
                skipped_normal_rna += 1
                subject_id = infer_subject_id(sample_id)
                if mapping and r1.name in mapping and mapping[r1.name].get("subject_id"):
                    subject_id = mapping[r1.name]["subject_id"]
                library_id = prefix
                normal_rna_rows.append(
                    {
                        "subject_id": subject_id,
                        "sample_id": sample_id,
                        "sample_type": sample_type,
                        "sequence_type": sequence_type,
                        "filetype": "fastq",
                        "info": f"library_id:{library_id};lane:{lane}",
                        "filepath": f"{r1};{r2}",
                        "is_plasma": False,
                    }
                )
                continue

            subject_id = infer_subject_id(sample_id)
            if mapping and r1.name in mapping and mapping[r1.name].get("subject_id"):
                subject_id = mapping[r1.name]["subject_id"]
            library_id = prefix

            info_parts = [f"library_id:{library_id}", f"lane:{lane}"]
            if is_plasma:
                info_parts.append("longitudinal_sample")
            info_str = ";".join(info_parts)

            local_rows.append(
                {
                    "subject_id": subject_id,
                    "sample_id": sample_id,
                    "sample_type": sample_type,
                    "sequence_type": sequence_type,
                    "filetype": "fastq",
                    "info": info_str,
                    "filepath": f"{r1};{r2}",
                    "is_plasma": is_plasma,
                }
            )
        return local_rows

    def scan_bam_dir(d: Path) -> List[Path]:
        files = []
        for ext in ["*.bam", "*.cram"]:
            files.extend(d.glob(ext))
            files.extend(d.glob(f"*/*/{ext}"))
            files.extend(d.glob(f"*/*/*/{ext}"))
        return sorted([fp.resolve() for fp in files if fp.is_file()])

    def process_bam_files(
        bam_files: List[Path],
        dir_default_sample_type: str
    ) -> List[Dict[str, str]]:
        nonlocal skipped, skipped_normal_rna
        local_rows = []
        for bam in bam_files:
            sample_id = bam.stem
            for suffix in ["_realigned_recalibrated", ".realigned.recalibrated", ".redux"]:
                if sample_id.endswith(suffix):
                    sample_id = sample_id[:-len(suffix)]

            sample_type = dir_default_sample_type
            sequence_type = args.default_sequence_type
            is_plasma = False

            # Infer subject from directory path if present (e.g. .../<subject_id>/.../sample.bam)
            dir_subject = None
            for part in reversed(bam.parts[:-1]):
                if any(k in part.upper() for k in ["BFIERCE", "TRACEOC", "OV-P", "DYNA", "MOHCCN"]) or re.search(r"\b(P|PT|PAT|DONOR)[_-]?\d{2,5}\b", part, re.IGNORECASE):
                    dir_subject = part
                    break

            if mapping:
                meta = None
                if dir_subject:
                    meta = mapping.get((dir_subject, sample_id)) or mapping.get(f"{dir_subject}_{sample_id}")
                if not meta:
                    meta = mapping.get(bam.name) or mapping.get(bam.stem) or mapping.get(sample_id)
                if not meta:
                    skipped += 1
                    print(
                        f"WARN: {bam.name} not found in mapping TSV; skipping file",
                        file=sys.stderr,
                    )
                    continue
                sample_id = meta["sample_id"] or sample_id
                tissue_origin = meta.get("tissue_origin", "")
                sample_type = infer_sample_type(
                    meta["tissue_type"], sample_id, dir_default_sample_type, tissue_origin=tissue_origin
                )
                sequence_type = infer_sequence_type(meta["library_type"], args.default_sequence_type)
                is_plasma = is_plasma_meta(meta)

            is_redux = bam.name.lower().endswith((".redux.bam", ".redux.cram")) or args.generate_redux_tsvs_only
            ftype = "cram_redux" if (bam.suffix.lower() == ".cram" and is_redux) else (
                "bam_redux" if is_redux else ("cram" if bam.suffix.lower() == ".cram" else "bam")
            )

            subject_id = dir_subject or infer_subject_id(sample_id)
            if mapping and meta and meta.get("subject_id"):
                subject_id = meta["subject_id"]

            if sample_type == "normal" and sequence_type == "rna":
                skipped_normal_rna += 1
                normal_rna_rows.append(
                    {
                        "subject_id": subject_id,
                        "sample_id": sample_id,
                        "sample_type": sample_type,
                        "sequence_type": sequence_type,
                        "filetype": "cram" if bam.suffix.lower() == ".cram" else "bam",
                        "info": "",
                        "filepath": str(bam),
                        "is_plasma": False,
                    }
                )
                continue

            info_parts = []
            if is_plasma:
                info_parts.append("longitudinal_sample")
            if args.generate_redux_tsvs_only or is_redux:
                info_parts.append("generate_redux_tsvs_only")
            info_str = ";".join(info_parts)

            local_rows.append(
                {
                    "subject_id": subject_id,
                    "sample_id": sample_id,
                    "sample_type": sample_type,
                    "sequence_type": sequence_type,
                    "filetype": ftype,
                    "info": info_str,
                    "filepath": str(bam),
                    "is_plasma": is_plasma,
                }
            )
        return local_rows

    # Gather FASTQ rows
    fastq_rows: List[Dict[str, str]] = []
    fastq_normal_candidates: List[Dict[str, str]] = []

    if args.fastq_dir:
        fastq_dir = Path(args.fastq_dir)
        if not fastq_dir.is_dir():
            print(f"ERROR: --fastq-dir is not a directory: {fastq_dir}", file=sys.stderr)
            return 1
        lane_reads = scan_dir(fastq_dir)
        if not lane_reads:
            print("ERROR: No FASTQ files matching *_L###_R[12]_*.fastq.gz found", file=sys.stderr)
            return 1
        fastq_rows = process_lane_reads(lane_reads, args.default_sample_type)
    elif args.tumor_fastq_dir and args.normal_fastq_dir:
        tumor_dir = Path(args.tumor_fastq_dir)
        normal_dir = Path(args.normal_fastq_dir)
        if not tumor_dir.is_dir():
            print(f"ERROR: --tumor-fastq-dir is not a directory: {tumor_dir}", file=sys.stderr)
            return 1
        if not normal_dir.is_dir():
            print(f"ERROR: --normal-fastq-dir is not a directory: {normal_dir}", file=sys.stderr)
            return 1

        tumor_lane_reads = scan_dir(tumor_dir)
        normal_lane_reads = scan_dir(normal_dir)

        if not tumor_lane_reads and not normal_lane_reads:
            print("ERROR: No FASTQ files matching *_L###_R[12]_*.fastq.gz found in paired FASTQ directories", file=sys.stderr)
            return 1

        if tumor_lane_reads:
            fastq_rows = process_lane_reads(tumor_lane_reads, "tumor")
        if normal_lane_reads:
            fastq_normal_candidates = process_lane_reads(normal_lane_reads, "normal")

    # Gather BAM/CRAM rows
    bam_rows: List[Dict[str, str]] = []
    bam_normal_candidates: List[Dict[str, str]] = []

    if args.bam_dir:
        bam_dir = Path(args.bam_dir)
        if not bam_dir.is_dir():
            print(f"ERROR: --bam-dir is not a directory: {bam_dir}", file=sys.stderr)
            return 1
        bam_files = scan_bam_dir(bam_dir)
        if not bam_files:
            print("ERROR: No BAM/CRAM files found in BAM directory", file=sys.stderr)
            return 1
        bam_rows = process_bam_files(bam_files, args.default_sample_type)
    elif args.tumor_bam_dir and args.normal_bam_dir:
        tumor_b_dir = Path(args.tumor_bam_dir)
        normal_b_dir = Path(args.normal_bam_dir)
        if not tumor_b_dir.is_dir():
            print(f"ERROR: --tumor-bam-dir is not a directory: {tumor_b_dir}", file=sys.stderr)
            return 1
        if not normal_b_dir.is_dir():
            print(f"ERROR: --normal-bam-dir is not a directory: {normal_b_dir}", file=sys.stderr)
            return 1

        tumor_bam_files = scan_bam_dir(tumor_b_dir)
        normal_bam_files = scan_bam_dir(normal_b_dir)

        if not tumor_bam_files and not normal_bam_files:
            print("ERROR: No BAM/CRAM files found in paired BAM directories", file=sys.stderr)
            return 1

        if tumor_bam_files:
            bam_rows = process_bam_files(tumor_bam_files, "tumor")
        if normal_bam_files:
            bam_normal_candidates = process_bam_files(normal_bam_files, "normal")

    # Determine tumor subjects for paired normal matching
    all_tumor_rows = [r for r in (redux_bam_rows + fastq_rows + bam_rows) if r["sample_type"] == "tumor"]
    tumor_subjects = {r["subject_id"] for r in all_tumor_rows}

    matched_normal_rows: List[Dict[str, str]] = []
    if args.tumor_fastq_dir:
        for r in fastq_normal_candidates:
            if r["subject_id"] in tumor_subjects:
                matched_normal_rows.append(r)
        print(f"Matched {len(matched_normal_rows)} out of {len(fastq_normal_candidates)} normal/reference FASTQ lanes to tumor subjects.")

    matched_normal_bams: List[Dict[str, str]] = []
    if args.tumor_bam_dir:
        for r in bam_normal_candidates:
            if r["subject_id"] in tumor_subjects:
                matched_normal_bams.append(r)
        print(f"Matched {len(matched_normal_bams)} out of {len(bam_normal_candidates)} normal/reference BAM files to tumor subjects.")

    config_rows: List[Dict[str, str]] = []
    if args.config_tsv:
        config_path = Path(args.config_tsv)
        if not config_path.is_file():
            print(f"ERROR: --config-tsv file not found: {config_path}", file=sys.stderr)
            return 1
        with config_path.open("r", newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"No header in config TSV: {config_path}")
            header_map = {norm_header(h): h for h in reader.fieldnames}

            def get_col_cfg(candidates: List[str], required: bool = True) -> str | None:
                for c in candidates:
                    norm_c = norm_header(c)
                    if norm_c in header_map:
                        return header_map[norm_c]
                if required:
                    raise ValueError(f"Required column not found from {candidates}. Found: {reader.fieldnames}")
                return None

            file_col = get_col_cfg(["fastq file name", "file", "filepath", "files", "bam"])
            sample_col = get_col_cfg(["sample id", "sample_id", "sample"])
            tissue_col = get_col_cfg(["tissue type (eg. primary tumor, reference)", "tissue type", "tissue_type"], required=False)
            lib_col = get_col_cfg(["library type (wg:whole genome, wt:whole transcriptome, pg:plasma, ex:exome)", "library type", "library_type"], required=False)
            subject_col = get_col_cfg(["donor_id", "subject_id", "subject id", "oicr_donor"], required=False)
            tissue_origin_col = get_col_cfg(["tissue origin", "tissue_origin"], required=False)

            for r in reader:
                fp_val = (r.get(file_col) or "").strip().strip('"')
                if not fp_val:
                    continue
                raw_sample_id = (r.get(sample_col) or "").strip().strip('"')
                sample_id = re.sub(r"\s+", "-", raw_sample_id)
                tissue_type = (r.get(tissue_col) or "").strip().strip('"') if tissue_col else ""
                library_type = (r.get(lib_col) or "").strip().strip('"') if lib_col else ""
                tissue_origin = (r.get(tissue_origin_col) or "").strip().strip('"') if tissue_origin_col else ""
                if subject_col and r.get(subject_col):
                    subject_id = re.sub(r"\s+", "-", r[subject_col].strip().strip('"'))
                else:
                    subject_id = infer_subject_id(sample_id)

                meta_dict = {
                    "sample_id": sample_id,
                    "tissue_type": tissue_type,
                    "library_type": library_type,
                    "tissue_origin": tissue_origin,
                    "subject_id": subject_id,
                }
                is_plasma = is_plasma_meta(meta_dict)
                sample_type = infer_sample_type(tissue_type, sample_id, args.default_sample_type, tissue_origin=tissue_origin)
                sequence_type = infer_sequence_type(library_type, args.default_sequence_type)

                fp_lower = fp_val.lower()
                is_redux_bam = fp_lower.endswith(".redux.bam")
                is_redux_cram = fp_lower.endswith(".redux.cram")

                if is_redux_bam:
                    ftype = "bam_redux"
                elif fp_lower.endswith(".bam"):
                    ftype = "bam"
                elif is_redux_cram:
                    ftype = "cram_redux"
                elif fp_lower.endswith(".cram"):
                    ftype = "cram"
                elif fp_lower.endswith((".fastq.gz", ".fq.gz")):
                    ftype = "fastq"
                else:
                    ftype = "bam"

                info_parts = []
                if is_plasma:
                    info_parts.append("longitudinal_sample")
                if ftype in {"bam_redux", "cram_redux"} or args.generate_redux_tsvs_only:
                    info_parts.append("generate_redux_tsvs_only")
                info_str = ";".join(info_parts)

                base_row = {
                    "subject_id": subject_id,
                    "sample_id": sample_id,
                    "sample_type": sample_type,
                    "sequence_type": sequence_type,
                    "filetype": ftype,
                    "info": info_str,
                    "filepath": fp_val,
                    "is_plasma": is_plasma,
                }
                config_rows.append(base_row)

    # Combine everything
    rows = []
    rows.extend(redux_bam_rows)
    rows.extend(config_rows)
    rows.extend(fastq_rows)
    rows.extend(bam_rows)
    rows.extend(matched_normal_rows)
    rows.extend(matched_normal_bams)

    # Filter normal RNA rows based on matched tumor subjects
    normal_rna_rows = [r for r in normal_rna_rows if r["subject_id"] in tumor_subjects]

    if not rows:
        print("ERROR: No rows generated. Check FASTQ/BAM names and mapping TSV.", file=sys.stderr)
        return 1

    # Resolve oncoanalyser output directory for referencing Step 1 outputs in Step 2 sheets
    if args.oncoanalyser_outdir:
        onco_outdir = str(args.oncoanalyser_outdir).rstrip("/")
    else:
        onco_outdir = str(Path(args.out).parent).rstrip("/")

    # Build analysis groups per subject.
    # Step 1: WGTS mode receives primary tumor tissue DNA (+ matched normal DNA/RNA).
    # Step 2: Purity Estimate mode receives longitudinal plasma ctDNA + primary normal BAM + primary purple/amber dirs.
    grouped_rows: List[Dict[str, str]] = []
    purity_sheets_by_subject: Dict[str, Dict[str, List[Dict[str, str]]]] = defaultdict(dict)
    rows_by_subject: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        rows_by_subject[row["subject_id"]].append(row)

    for subject_id, subject_rows in rows_by_subject.items():
        primary_tumor_dna_rows = [
            r for r in subject_rows
            if r["sample_type"] == "tumor" and r["sequence_type"] == "dna" and not r.get("is_plasma", False)
        ]
        plasma_tumor_dna_rows = [
            r for r in subject_rows
            if r["sample_type"] == "tumor" and r["sequence_type"] == "dna" and r.get("is_plasma", False)
        ]
        normal_dna_rows = [
            r for r in subject_rows
            if r["sample_type"] == "normal" and r["sequence_type"] == "dna"
        ]

        if primary_tumor_dna_rows:
            tumor_dna_anchors = sorted(
                {infer_tumor_anchor(r["sample_id"]) for r in primary_tumor_dna_rows}
            )
        elif plasma_tumor_dna_rows:
            tumor_dna_anchors = sorted(
                {infer_tumor_anchor(r["sample_id"]) for r in plasma_tumor_dna_rows}
            )
        else:
            skipped_subjects_no_tumor_dna += 1
            continue

        tumor_rna_anchors = {
            infer_tumor_anchor(r["sample_id"])
            for r in subject_rows
            if r["sample_type"] == "tumor" and r["sequence_type"] == "rna"
        }
        skipped_orphan_tumor_rna += len(tumor_rna_anchors - set(tumor_dna_anchors))

        distinct_plasma_sids = sorted(set(r["sample_id"] for r in plasma_tumor_dna_rows))

        # --- Step 1: WGTS Candidate Rows ---
        # When primary tissue tumor DNA is present, exclude plasma rows so Nextflow only aligns the genuine primary tumor.
        for anchor in tumor_dna_anchors:
            candidate_rows: List[Dict[str, str]] = []
            for r in subject_rows:
                if r["sample_type"] == "normal" and r["sequence_type"] == "dna":
                    candidate_rows.append(dict(r))
                elif (
                    r["sample_type"] == "tumor"
                    and r["sequence_type"] in {"dna", "rna"}
                ):
                    if primary_tumor_dna_rows:
                        if not r.get("is_plasma", False) and infer_tumor_anchor(r["sample_id"]) == anchor:
                            candidate_rows.append(dict(r))
                    else:
                        # Fallback for plasma-only cohorts
                        if infer_tumor_anchor(r["sample_id"]) == anchor:
                            candidate_rows.append(dict(r))

            has_tumor_dna = any(
                r["sample_type"] == "tumor" and r["sequence_type"] == "dna"
                for r in candidate_rows
            )
            if not has_tumor_dna:
                continue

            tumor_dna_sample_ids = {
                r["sample_id"] for r in candidate_rows
                if r["sample_type"] == "tumor" and r["sequence_type"] == "dna" and not r.get("is_plasma", False)
            }

            for r in candidate_rows:
                r_copy = dict(r)
                r_copy["group_id"] = anchor
                if (
                    r_copy["sample_type"] == "tumor"
                    and r_copy["sequence_type"] == "rna"
                    and r_copy["sample_id"] in tumor_dna_sample_ids
                ):
                    r_copy["sample_id"] = f"{r_copy['sample_id']}-RNA"
                grouped_rows.append(r_copy)

        # --- Step 2: Purity Estimate (MRD / WISP) Candidate Rows ---
        # If both primary tumor and longitudinal plasma ctDNA exist, build the purity estimate sheet(s)
        if primary_tumor_dna_rows and plasma_tumor_dna_rows:
            primary_anchor = tumor_dna_anchors[0]
            primary_sid = primary_tumor_dna_rows[0]["sample_id"]
            normal_sid = normal_dna_rows[0]["sample_id"] if normal_dna_rows else None

            ref_outdir = args.purity_reference_dir or onco_outdir
            for p_sid in distinct_plasma_sids:
                p_rows: List[Dict[str, str]] = []

                # 1. Plasma sample rows for this plasma sample
                for r in plasma_tumor_dna_rows:
                    if r["sample_id"] == p_sid:
                        r_copy = dict(r)
                        r_copy["group_id"] = primary_anchor
                        info_parts = [p for p in (r_copy.get("info") or "").split(";") if p]
                        if "longitudinal_sample" not in info_parts:
                            info_parts.append("longitudinal_sample")
                        if r_copy.get("filetype") in {"bam_redux", "cram_redux"} and "generate_redux_tsvs_only" not in info_parts:
                            info_parts.append("generate_redux_tsvs_only")
                        r_copy["info"] = ";".join(info_parts)
                        p_rows.append(r_copy)

                # 2. Matched Normal DNA alignment (redux_dir from Step 1)
                # Nextflow natively supports filetype=redux_dir for normal alignments in purity_estimate mode.
                # Pointing to the Step 1 alignment folder causes Nextflow to skip REDUX on the normal,
                # avoiding redundant re-processing and preventing a duplicate 30-280GB normal BAM/CRAM from
                # being written into the plasma output directory.
                if normal_sid:
                    normal_redux_path = None
                    if normal_dna_rows and normal_dna_rows[0].get("filetype") == "redux_dir":
                        normal_redux_path = normal_dna_rows[0]["filepath"]
                    else:
                        candidate_dirs = [
                            Path(ref_outdir) / subject_id / primary_anchor / "alignments" / normal_sid,
                            Path(ref_outdir) / subject_id / primary_anchor / "alignments" / "dna" / normal_sid,
                            Path(onco_outdir) / subject_id / primary_anchor / "alignments" / normal_sid,
                            Path(onco_outdir) / subject_id / primary_anchor / "alignments" / "dna" / normal_sid,
                        ]
                        for cd in candidate_dirs:
                            if cd.is_dir():
                                normal_redux_path = str(cd)
                                break
                        if not normal_redux_path:
                            normal_redux_path = f"{ref_outdir}/{subject_id}/{primary_anchor}/alignments/{normal_sid}"

                    p_rows.append({
                        "group_id": primary_anchor,
                        "subject_id": subject_id,
                        "sample_id": normal_sid,
                        "sample_type": "normal",
                        "sequence_type": "dna",
                        "filetype": "redux_dir",
                        "info": "",
                        "filepath": normal_redux_path,
                    })

                # 3. Primary tumor PURPLE directory from Step 1
                purple_dir_path = f"{ref_outdir}/{subject_id}/{primary_anchor}/purple"
                p_rows.append({
                    "group_id": primary_anchor,
                    "subject_id": subject_id,
                    "sample_id": primary_sid,
                    "sample_type": "tumor",
                    "sequence_type": "dna",
                    "filetype": "purple_dir",
                    "info": "",
                    "filepath": purple_dir_path,
                })

                # 4. Primary tumor AMBER directory from Step 1
                amber_dir_path = f"{ref_outdir}/{subject_id}/{primary_anchor}/amber"
                p_rows.append({
                    "group_id": primary_anchor,
                    "subject_id": subject_id,
                    "sample_id": primary_sid,
                    "sample_type": "tumor",
                    "sequence_type": "dna",
                    "filetype": "amber_dir",
                    "info": "",
                    "filepath": amber_dir_path,
                })

                purity_sheets_by_subject[subject_id][p_sid] = p_rows

    if not grouped_rows:
        print(
            "ERROR: No oncoanalyser-compatible rows were generated. "
            "Ensure there is at least one tumor DNA sample per patient.",
            file=sys.stderr,
        )
        return 1

    # Group rows by subject_id and write one samplesheet per patient.
    out_root = Path(args.out)
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

    patients: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in grouped_rows:
        patients[row["subject_id"]].append(row)

    normal_rna_patients: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in normal_rna_rows:
        normal_rna_patients[row["subject_id"]].append(row)

    written_total = 0
    for subject_id, patient_rows in sorted(patients.items()):
        patient_dir = out_root / subject_id
        patient_dir.mkdir(parents=True, exist_ok=True)
        out_path = patient_dir / "samplesheet.csv"
        with out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(patient_rows)
        print(f"  {subject_id}: {len(patient_rows)} WGTS rows -> {out_path}")
        written_total += len(patient_rows)

        if subject_id in purity_sheets_by_subject:
            purity_dict = purity_sheets_by_subject[subject_id]
            for p_sid, p_rows in purity_dict.items():
                purity_out_path = patient_dir / f"samplesheet.purity.{p_sid}.csv"
                with purity_out_path.open("w", newline="", encoding="utf-8") as fh:
                    writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
                    writer.writeheader()
                    writer.writerows(p_rows)
                print(f"    wrote {len(p_rows)} purity (MRD/WISP) rows -> {purity_out_path}")
                primary_purity_link = patient_dir / "samplesheet.purity.csv"
                if not primary_purity_link.exists():
                    try:
                        primary_purity_link.symlink_to(purity_out_path.name)
                    except OSError:
                        pass

        if args.emit_normal_rna_sheet and subject_id in normal_rna_patients:
            normal_rows = []
            for r in normal_rna_patients[subject_id]:
                row = dict(r)
                row["group_id"] = f"{subject_id}_NORMAL_RNA"
                normal_rows.append(row)
            normal_out_path = patient_dir / "samplesheet.normal_rna.csv"
            with normal_out_path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(normal_rows)
            print(f"    wrote {len(normal_rows)} normal RNA rows -> {normal_out_path}")

    print(f"\nWrote {written_total} total rows across {len(patients)} patient(s) under {out_root}")
    if skipped:
        print(f"Skipped {skipped} lane/file entries (see warnings)")
    if skipped_normal_rna:
        print(f"Skipped {skipped_normal_rna} normal RNA lane/file entries from main oncoanalyser samplesheet")
    if skipped_orphan_tumor_rna:
        print(f"Skipped {skipped_orphan_tumor_rna} tumor RNA lane/file entries without matching tumor DNA anchor")
    if skipped_subjects_no_tumor_dna:
        print(f"Skipped {skipped_subjects_no_tumor_dna} subject(s) with no tumor DNA inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
