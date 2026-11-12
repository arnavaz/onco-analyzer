#!/usr/bin/env python3
"""Automated Cohort Intake & Runner Generator for Oncoanalyser.

Ingests FASTQ map files across multiple release batches and arbitrary schemas,
resolves duplicate FASTQs across archive/transfer paths, normalizes sample metadata,
and produces:
  1. A unified, deduplicated mapping TSV: runs/metadata/<COHORT>.fastqs.map.tsv
  2. A deduplicated FASTQ paths list:   runs/metadata/<COHORT>_fastq_list.dedup.txt
  3. A production-ready SLURM runner:     runs/run_<COHORT>.sh

Usage:
  python3 Setup/prepare_cohort_run.py \
    --cohort-name NETMRD \
    --map-files metadata/*.fastqs.map.tsv \
    --fastq-list /path/to/fastq_list.txt \
    --exclude-duplicate-pattern "OICR_data/fastqs" \
    --account pughlab \
    --partition pughlab \
    --outdir /cluster/projects/rskimlab/TEMP/NETMRD/Oncoanalyser \
    --target-link-dir /cluster/projects/rskimlab/TEMP/NETMRD/all_fastq_links
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


STANDARD_FIELDS = [
    "oicr_case",
    "oicr_donor",
    "donor_id",
    "library_id",
    "library_type",
    "tissue_type",
    "tissue_origin",
    "run",
    "barcode",
    "sample_id",
    "files",
]


def norm_h(h: str) -> str:
    return re.sub(r"\s+", " ", h.strip().strip('"').lower())


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Cohort Intake & Oncoanalyser Runner Generator")
    p.add_argument("--cohort-name", required=True, help="Cohort identifier (e.g. NETMRD, BFIERCE)")
    p.add_argument("--map-files", nargs="+", required=True, help="One or more fastq map TSVs (in release order)")
    p.add_argument("--fastq-list", default=None, help="Text file containing full cluster paths to fastq.gz files")
    p.add_argument("--fastq-dir", default=None, help="Directory containing fastq.gz files")
    p.add_argument("--exclude-duplicate-pattern", default="", help="Substring pattern to exclude when duplicate fastqs exist (e.g. 'OICR_data/fastqs')")
    p.add_argument("--account", default="pughlab", help="SLURM account (default: pughlab)")
    p.add_argument("--partition", default=None, help="SLURM partition (default: auto based on account)")
    p.add_argument("--outdir", default=None, help="Root Oncoanalyser output directory")
    p.add_argument("--target-link-dir", default=None, help="Directory where FASTQs will be symlinked on cluster")
    p.add_argument("--source-repo", default="/cluster/projects/pughlab/src/oncoanalyser_pughlab", help="Pipeline source repository on cluster")
    p.add_argument("--output-script", default=None, help="Output path for the generated bash runner script")
    p.add_argument("--output-map", default=None, help="Output path for the unified map TSV")
    p.add_argument("--output-fastq-list", default=None, help="Output path for the deduplicated fastq list")
    return p.parse_args()


def load_and_harmonize_maps(map_paths: List[Path]) -> Tuple[List[Dict[str, str]], Set[str]]:
    merged_rows: Dict[Tuple[str, ...], Dict[str, str]] = {}
    all_fastqs: Set[str] = set()

    for mp in map_paths:
        if not mp.is_file():
            print(f"WARN: Mapping file not found: {mp}", file=sys.stderr)
            continue
        with mp.open("r", newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                continue

            for r in reader:
                files_raw = r.get("files") or r.get("fastq file name") or r.get("filepath") or ""
                if not files_raw.strip():
                    continue

                filenames = [f.strip() for f in files_raw.split(";") if f.strip()]
                all_fastqs.update(filenames)
                f_key = tuple(sorted(filenames))

                # Flexible header resolution across release schemas
                oicr_case = r.get("oicr_case") or r.get("case_id") or r.get("sample") or ""
                oicr_donor = r.get("oicr_donor") or r.get("sample") or r.get("case_id") or ""
                
                # In 2024 schema, clinical donor was in 'sample_id'; in 2025/2026 schema, in 'donor_id'
                donor_id = r.get("donor_id") or (r.get("sample_id") if "group_id" in r else "")
                if not donor_id:
                    donor_id = oicr_donor or oicr_case

                library_id = r.get("library_id") or r.get("library") or ""
                library_type = r.get("library_type") or r.get("library_source") or "WG"
                
                # In 2024 schema, sample ID was in 'group_id'; in 2025/2026, in 'sample_id'
                sample_id = r.get("sample_id") if "group_id" not in r else r.get("group_id")
                if not sample_id:
                    sample_id = library_id

                tissue_type = r.get("tissue_type") or ""
                tissue_origin = r.get("tissue_origin") or ""
                run = r.get("run") or ""
                barcode = r.get("barcode") or ""

                merged_rows[f_key] = {
                    "oicr_case": oicr_case,
                    "oicr_donor": oicr_donor,
                    "donor_id": donor_id,
                    "library_id": library_id,
                    "library_type": library_type,
                    "tissue_type": tissue_type,
                    "tissue_origin": tissue_origin,
                    "run": run,
                    "barcode": barcode,
                    "sample_id": sample_id,
                    "files": ";".join(filenames),
                }

    sorted_rows = sorted(
        merged_rows.values(),
        key=lambda x: (x["donor_id"], x["library_id"], x["run"])
    )
    return sorted_rows, all_fastqs


def deduplicate_fastqs(paths: List[str], exclude_pattern: str) -> Tuple[List[str], List[str]]:
    by_base = defaultdict(list)
    for p in paths:
        by_base[os.path.basename(p)].append(p)

    kept: List[str] = []
    excluded: List[str] = []

    for base, p_list in sorted(by_base.items()):
        if len(p_list) == 1:
            kept.append(p_list[0])
        else:
            if exclude_pattern:
                non_excluded = [p for p in p_list if exclude_pattern not in p]
                has_excluded = [p for p in p_list if exclude_pattern in p]
                if non_excluded and has_excluded:
                    kept.extend(non_excluded)
                    excluded.extend(has_excluded)
                    continue
            # Fallback: keep the first path
            kept.append(p_list[0])
            excluded.extend(p_list[1:])

    return sorted(kept), sorted(excluded)


def generate_runner_script(
    cohort_name: str,
    target_link_dir: str,
    fastq_list_cluster: str,
    mapping_tsv_cluster: str,
    outdir_cluster: str,
    account: str,
    partition: str,
    source_repo: str,
) -> str:
    return f"""#!/bin/bash

# ==============================================================================
# {cohort_name} Workflow Runner: Step 1 (WGTS) + Step 2 (Longitudinal Plasma WISP/MRD)
# Generated automatically by Setup/prepare_cohort_run.py
# ==============================================================================

TARGET_LINK_DIR="{target_link_dir}"
mkdir -p "${{TARGET_LINK_DIR}}"

# 1. Symlink using the deduplicated FASTQ list
SOURCE_REPO="{source_repo}"
FASTQ_LIST="{fastq_list_cluster}"
[[ -f "${{FASTQ_LIST}}" ]] || FASTQ_LIST="${{SOURCE_REPO}}/runs/metadata/{cohort_name}_fastq_list.dedup.txt"

if [[ -f "${{FASTQ_LIST}}" ]]; then
    while IFS= read -r source_path; do
        [[ -z "$source_path" || "$source_path" =~ ^# ]] && continue
        filename=$(basename "$source_path")
        ln -sf "$source_path" "${{TARGET_LINK_DIR}}/${{filename}}"
    done < "${{FASTQ_LIST}}"
fi

# 2. Configuration & Sample Sheets
MAPPING_TSV="{mapping_tsv_cluster}"
[[ -f "${{MAPPING_TSV}}" ]] || MAPPING_TSV="${{SOURCE_REPO}}/runs/metadata/{cohort_name}.fastqs.map.tsv"

OUTDIR="${{OUTDIR:-{outdir_cluster}}}"
SHEETS_DIR="${{OUTDIR}}/sample_sheets"

# SLURM Account & Partition configuration:
ACCOUNT="${{ACCOUNT:-{account}}}"
PARTITION="${{PARTITION:-{partition}}}"

# Guard against known invalid account/partition combinations on H4H
if [[ "${{ACCOUNT}}" == "deperrotgroup" && "${{PARTITION}}" == "pughlab" ]]; then
    echo "WARNING: Partition 'pughlab' is private to account 'pughlab'. Switching partition to 'all'."
    PARTITION="all"
fi

if [[ -z "${{WGTS_TIME:-}}" ]]; then
    if [[ "${{PARTITION}}" == "pughlab" || "${{PARTITION}}" == "long" ]]; then
        WGTS_TIME="14-00:00:00"
    elif [[ "${{PARTITION}}" == "himem" ]]; then
        WGTS_TIME="7-00:00:00"
    else
        WGTS_TIME="5-00:00:00"
    fi
fi

module load python3/3.12.11

echo "Generating Oncoanalyser sample sheets for {cohort_name}..."
python3 "${{SOURCE_REPO}}/Setup/make_oncoanalyser_samplesheet.py" \\
  --fastq-dir "${{TARGET_LINK_DIR}}" \\
  --mapping-tsv "${{MAPPING_TSV}}" \\
  --out "${{SHEETS_DIR}}" \\
  --oncoanalyser-outdir "${{OUTDIR}}"

export SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_SINGULARITY_CACHEDIR=/cluster/projects/pughlab/references/oncoanalyser_resources/singularity_cache
export NXF_OFFLINE=true

# Setup user-controlled temporary directory to prevent compute node /tmp exhaustion
export TMPDIR="${{OUTDIR}}/tmp"
export TEMP="${{OUTDIR}}/tmp"
export TMP="${{OUTDIR}}/tmp"
export SINGULARITY_TMPDIR="${{OUTDIR}}/tmp"
export APPTAINER_TMPDIR="${{OUTDIR}}/tmp"
export SINGULARITYENV_TMPDIR="${{OUTDIR}}/tmp"
export SINGULARITYENV_TEMP="${{OUTDIR}}/tmp"
export SINGULARITYENV_TMP="${{OUTDIR}}/tmp"
mkdir -p "${{TMPDIR}}"

echo ""
echo "=========================================================="
echo "Submitting {cohort_name} Oncoanalyser pipeline jobs..."
echo "Account:   ${{ACCOUNT}}"
echo "Partition: ${{PARTITION}}"
echo "Outdir:    ${{OUTDIR}}"
echo "Time Limit:${{WGTS_TIME}}"
echo "=========================================================="

for sample_sheet in "${{SHEETS_DIR}}"/*/samplesheet.csv; do
  [[ -f "$sample_sheet" ]] || continue
  patient=$(basename "$(dirname "$sample_sheet")")
  patient_outdir="${{OUTDIR}}/${{patient}}"
  patient_tmpdir="${{patient_outdir}}/tmp"
  
  mkdir -p "${{patient_outdir}}/logs" "${{patient_tmpdir}}"
  
  # Determine tumor anchor and check Step 1 (WGTS) status
  tumor_anchor=$(grep -m 1 ",tumor," "$sample_sheet" | cut -d',' -f1 || true)
  purple_dir="${{patient_outdir}}/${{tumor_anchor}}/purple"
  amber_dir="${{patient_outdir}}/${{tumor_anchor}}/amber"

  # Check if Step 1 (WGTS) is currently running or completed
  running_wgts=$(squeue -u "$USER" --name="onco_wgts_${{patient}}" -h -o "%i" | head -n 1)

  has_purple="false"
  if compgen -G "${{patient_outdir}}/*/purple/*.purple.purity.tsv" > /dev/null || [[ -d "$purple_dir" && -d "$amber_dir" ]]; then
    has_purple="true"
  fi

  wgts_job_id=""
  if [[ -n "$running_wgts" ]]; then
    echo "[${{patient}}] Step 1 (WGTS) is currently running in SLURM queue (Job ID: ${{running_wgts}})"
    wgts_job_id="$running_wgts"
  elif [[ "${{has_purple}}" == "true" ]]; then
    echo "[${{patient}}] Step 1 (WGTS) outputs already exist (${{purple_dir}}); skipping Step 1 submission."
  else
    echo "[${{patient}}] Submitting Step 1 (WGTS mode)..."
    sbatch_out=$(sbatch --parsable \\
      --partition="${{PARTITION}}" \\
      --time="${{WGTS_TIME}}" \\
      --export=ALL,SOURCE_REPO=${{SOURCE_REPO}},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${{patient_tmpdir}},TEMP=${{patient_tmpdir}},TMP=${{patient_tmpdir}},SINGULARITY_TMPDIR=${{patient_tmpdir}},APPTAINER_TMPDIR=${{patient_tmpdir}},SINGULARITYENV_TMPDIR=${{patient_tmpdir}},SINGULARITYENV_TEMP=${{patient_tmpdir}},SINGULARITYENV_TMP=${{patient_tmpdir}} \\
      -A "$ACCOUNT" \\
      --job-name="onco_wgts_${{patient}}" \\
      --output="${{patient_outdir}}/logs/onco_wgts_${{patient}}_%j.out" \\
      --error="${{patient_outdir}}/logs/onco_wgts_${{patient}}_%j.err" \\
      "${{SOURCE_REPO}}/run_oncoanalyser.sh" \\
      "$sample_sheet" \\
      "$patient_outdir" \\
      wgts \\
      GRCh38_hmf \\
      "$ACCOUNT" \\
      --processes_exclude neo,virusinterpreter 2>&1) || true

    if [[ "$sbatch_out" =~ ^[0-9]+ ]]; then
      wgts_job_id=$(echo "$sbatch_out" | cut -d';' -f1)
      echo "  Submitted WGTS job ID: ${{wgts_job_id}}"
    else
      echo "  ERROR: Failed to submit Step 1 (WGTS) for ${{patient}}: ${{sbatch_out}}"
    fi
  fi

  # Step 2: Longitudinal Plasma Purity Estimate (MRD / WISP & ichorCNA)
  purity_sheets=()
  for ps in "${{SHEETS_DIR}}/${{patient}}"/samplesheet.purity.*.csv; do
    [[ -f "$ps" ]] && purity_sheets+=("$ps")
  done
  if [[ ${{#purity_sheets[@]}} -eq 0 && -f "${{SHEETS_DIR}}/${{patient}}/samplesheet.purity.csv" ]]; then
    purity_sheets=("${{SHEETS_DIR}}/${{patient}}/samplesheet.purity.csv")
  fi

  if [[ ${{#purity_sheets[@]}} -gt 0 ]]; then
    for purity_sheet in "${{purity_sheets[@]}}"; do
      [[ -f "$purity_sheet" ]] || continue
      p_name=$(basename "$purity_sheet" .csv)
      p_name="${{p_name#samplesheet.purity.}}"
      p_name="${{p_name#samplesheet.plasma.}}"
      p_name="${{p_name#samplesheet.}}"
      if [[ -z "$p_name" || "$p_name" == "purity" ]]; then
        p_name=$(awk -F',' '$7 ~ /longitudinal_sample/ || $4 ~ /^[pP]lasma/ {{print $3}}' "$purity_sheet" | head -n 1)
        p_name="${{p_name:-plasma}}"
      fi

      plasma_outdir="${{patient_outdir}}/plasma/${{p_name}}"
      timepoint_tmpdir="${{plasma_outdir}}/tmp"
      mkdir -p "${{plasma_outdir}}" "${{patient_outdir}}/logs" "${{timepoint_tmpdir}}"

      # 1. Skip if both WISP and ichorCNA are already complete
      has_wisp=$(compgen -G "${{plasma_outdir}}/*/wisp/*.wisp.summary.tsv" "${{plasma_outdir}}/*/wisp/*.wisp.tsv" > /dev/null && echo "true" || echo "false")
      has_ichor=$(compgen -G "${{patient_outdir}}/ichorCNA_plasma_vs_n/${{p_name}}/*.params.txt" "${{patient_outdir}}/ichorCNA_plasma_only/${{p_name}}/*.params.txt" "${{plasma_outdir}}/ichorCNA*/*.params.txt" > /dev/null && echo "true" || echo "false")
      if [[ "${{has_wisp}}" == "true" && "${{has_ichor}}" == "true" && "${{FORCE_ICHOR:-false}}" != "true" ]]; then
        echo "  [${{patient}}] Plasma (${{p_name}}) already completed (both WISP and ichorCNA exist); skipping."
        continue
      fi

      # 2. Check if currently running in SLURM queue
      running_plasma=$(squeue -u "$USER" --name="onco_plasma_${{patient}}_${{p_name}}" -h -o "%i" | head -n 1)
      if [[ -n "$running_plasma" ]]; then
        plasma_reason=$(squeue -j "$running_plasma" -h -o "%r" 2>/dev/null || true)
        if [[ "$plasma_reason" == "DependencyNeverSatisfied" ]]; then
          echo "  [${{patient}}] Plasma (${{p_name}}) Job ${{running_plasma}} is stuck in DependencyNeverSatisfied; cancelling..."
          scancel "$running_plasma" 2>/dev/null || true
          running_plasma=""
        else
          echo "  [${{patient}}] Plasma (${{p_name}}) already running in SLURM (Job ID: ${{running_plasma}}); skipping."
          continue
        fi
      fi

      # 3. Dependency gating: Only launch Step 2 if Step 1 is running or already complete
      walltime="24:00:00"
      extra_flags=()
      dep_args=()
      if [[ "${{has_wisp}}" == "true" ]]; then
        echo "  [${{patient}}] WISP already complete for ${{p_name}} -> Running ichorCNA post-processing only (--skip-nextflow)."
        walltime="12:00:00"
        extra_flags=(--skip-nextflow)
      elif [[ -n "$wgts_job_id" ]]; then
        dep_args=(--dependency=afterok:"${{wgts_job_id}}")
        echo "  [${{patient}}] Queuing Step 2 (${{p_name}}) with dependency on WGTS Job ${{wgts_job_id}}..."
      elif [[ "${{has_purple}}" == "true" ]]; then
        echo "  [${{patient}}] Launching Step 2 (${{p_name}}) in parallel (WGTS outputs exist)..."
      else
        echo "  [${{patient}}] SKIPPING Step 2 (${{p_name}}): Step 1 (WGTS) has not completed and is not running."
        continue
      fi

      sbatch "${{dep_args[@]}}" \\
        --partition="${{PARTITION}}" \\
        --time="${{walltime}}" \\
        --cpus-per-task=8 \\
        --mem=32G \\
        --export=ALL,SOURCE_REPO=${{SOURCE_REPO}},NXF_OFFLINE=true,RESUME=true,CONVERT_FINAL_BAM_TO_CRAM=false,TMPDIR=${{timepoint_tmpdir}},TEMP=${{timepoint_tmpdir}},TMP=${{timepoint_tmpdir}},SINGULARITY_TMPDIR=${{timepoint_tmpdir}},APPTAINER_TMPDIR=${{timepoint_tmpdir}},SINGULARITYENV_TMPDIR=${{timepoint_tmpdir}},SINGULARITYENV_TEMP=${{timepoint_tmpdir}},SINGULARITYENV_TMP=${{timepoint_tmpdir}} \\
        -A "$ACCOUNT" \\
        --job-name="onco_plasma_${{patient}}_${{p_name}}" \\
        --output="${{patient_outdir}}/logs/onco_plasma_${{patient}}_${{p_name}}_%j.out" \\
        --error="${{patient_outdir}}/logs/onco_plasma_${{patient}}_${{p_name}}_%j.err" \\
        "${{SOURCE_REPO}}/run_oncoanalyser.sh" \\
        "$purity_sheet" \\
        "$plasma_outdir" \\
        purity_estimate \\
        GRCh38_hmf \\
        "$ACCOUNT" \\
        --purity_estimate_mode wgts \\
        --processes_exclude neo,virusinterpreter \\
        "${{extra_flags[@]}}"
    done
  else
    echo "  No purity/MRD sample sheet found for ${{patient}}; skipping Step 2."
  fi
done
"""


def main() -> int:
    args = parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    cohort = args.cohort_name.strip()

    # Expand map files (supports globs)
    expanded_map_paths: List[Path] = []
    for mf in args.map_files:
        matches = glob.glob(mf)
        if matches:
            expanded_map_paths.extend(Path(m) for m in matches)
        else:
            expanded_map_paths.append(Path(mf))

    print(f"Reading {len(expanded_map_paths)} mapping TSV file(s)...")
    harmonized_rows, map_fastqs = load_and_harmonize_maps(expanded_map_paths)
    print(f"  Loaded {len(harmonized_rows)} unique lane entries across {len(map_fastqs)} FASTQs.")

    # FASTQ list discovery and deduplication
    raw_fastq_paths: List[str] = []
    if args.fastq_list:
        fl_path = Path(args.fastq_list)
        if fl_path.is_file():
            with fl_path.open("r", encoding="utf-8") as fh:
                raw_fastq_paths = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        else:
            print(f"WARN: --fastq-list file not found: {fl_path}", file=sys.stderr)
    elif args.fastq_dir:
        fd_path = Path(args.fastq_dir)
        if fd_path.is_dir():
            raw_fastq_paths = [str(p) for p in fd_path.glob("*.f*q.gz")]

    dedup_paths: List[str] = []
    if raw_fastq_paths:
        dedup_paths, excluded_paths = deduplicate_fastqs(raw_fastq_paths, args.exclude_duplicate_pattern)
        print(f"  FASTQs: {len(raw_fastq_paths)} paths provided -> {len(dedup_paths)} deduplicated, {len(excluded_paths)} excluded.")
        
        # Verify coverage
        kept_basenames = {os.path.basename(p) for p in dedup_paths}
        missing = [fn for fn in map_fastqs if fn not in kept_basenames]
        if missing:
            print(f"  WARN: {len(missing)} FASTQs from mapping files were not found in the deduplicated path list!", file=sys.stderr)
            print(f"  Example missing: {missing[:5]}", file=sys.stderr)
        else:
            print(f"  Concordance: 100% of mapping FASTQs ({len(map_fastqs)}) verified in deduplicated path list.")

    # Paths for output files
    out_map = Path(args.output_map or (repo_root / "runs" / "metadata" / f"{cohort}.fastqs.map.tsv"))
    out_map.parent.mkdir(parents=True, exist_ok=True)
    with out_map.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=STANDARD_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(harmonized_rows)
    print(f"  Wrote harmonized mapping TSV: {out_map}")

    out_fastq_list = Path(args.output_fastq_list or (repo_root / "runs" / "metadata" / f"{cohort}_fastq_list.dedup.txt"))
    if dedup_paths:
        out_fastq_list.parent.mkdir(parents=True, exist_ok=True)
        with out_fastq_list.open("w", encoding="utf-8") as fh:
            for p in dedup_paths:
                fh.write(f"{p}\n")
        print(f"  Wrote deduplicated FASTQ list: {out_fastq_list}")

    # Output runner script
    account = args.account
    partition = args.partition
    if not partition:
        partition = "pughlab" if account == "pughlab" else "all"

    target_link_dir = args.target_link_dir or f"/cluster/projects/pughlab/projects/{cohort}/all_fastq_links"
    outdir_cluster = args.outdir or f"/cluster/projects/pughlab/projects/{cohort}/Oncoanalyser"
    fastq_list_cluster = f"/cluster/projects/pughlab/src/oncoanalyser_pughlab/runs/metadata/{cohort}_fastq_list.dedup.txt"
    mapping_tsv_cluster = f"/cluster/projects/pughlab/src/oncoanalyser_pughlab/runs/metadata/{cohort}.fastqs.map.tsv"

    script_content = generate_runner_script(
        cohort_name=cohort,
        target_link_dir=target_link_dir,
        fastq_list_cluster=fastq_list_cluster,
        mapping_tsv_cluster=mapping_tsv_cluster,
        outdir_cluster=outdir_cluster,
        account=account,
        partition=partition,
        source_repo=args.source_repo,
    )

    out_script = Path(args.output_script or (repo_root / "runs" / f"run_{cohort}.sh"))
    out_script.parent.mkdir(parents=True, exist_ok=True)
    with out_script.open("w", encoding="utf-8") as fh:
        fh.write(script_content)
    os.chmod(out_script, 0o755)
    print(f"  Wrote executable SLURM runner script: {out_script}")

    print("\nCohort intake successfully completed!")
    print(f"To deploy: git add {out_map} {out_fastq_list} {out_script} && git commit -m 'Setup {cohort}' && git pushall")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
