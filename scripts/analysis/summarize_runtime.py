#!/usr/bin/env python3
"""
summarize_runtime.py

Summarizes wall time and per-task runtime across an oncoanalyser cohort,
grouped by SLURM account.

Reads:
  - <results_base>/<patient>/logs/nf_*.out          (account, start/end date)
    - OR <results_base>/logs/nf_*.out                 (shared logs directory)
  - <results_base>/<patient>/pipeline_info/execution_trace.txt  (per-task timing)

Outputs (in <results_base>/runtime_summary/ unless --output-dir is set):
  - patient_summary.tsv          per-patient wall time and account
  - task_stats.tsv               per-task mean/median/max across cohort
  - patient_wall_time_hist.png   wall time distribution
  - wall_time_by_account.png     wall time grouped by account (if >1 account)
  - top_tasks_mean_duration.png  top N tasks by mean duration

Usage:
  module load python3/3.12.11
  python3 scripts/analysis/summarize_runtime.py <results_base> [options]

Options:
  --output-dir <path>    Output directory (default: <results_base>/runtime_summary)
  --top-tasks <N>        Top N tasks to plot (default: 20)
  --format png|pdf       Plot format (default: png)
"""

import argparse
import re
import sys
import subprocess
import os
import shlex
from pathlib import Path
from datetime import datetime

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as e:
    sys.exit(
        f"Missing dependency: {e}\n"
        "Install with: pip install pandas matplotlib"
    )


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

# Nextflow duration strings: e.g. "1d 2h 3m 4.5s", "23m 5s", "45.3s"
_NF_DUR_RE = re.compile(
    r"(?:(\d+)\s*d\s*)?(?:(\d+)\s*h\s*)?(?:(\d+)\s*m\s*)?(?:([\d.]+)\s*s)?"
)

def _parse_nf_duration(s: str) -> float:
    """Return duration in seconds, or NaN on failure."""
    if not s or s.strip() in ("-", ""):
        return float("nan")
    m = _NF_DUR_RE.fullmatch(s.strip())
    if not m or not any(m.groups()):
        return float("nan")
    d, h, mi, sec = m.groups(default="0")
    return int(d) * 86400 + int(h) * 3600 + int(mi) * 60 + float(sec or 0)


# Date line in the wrapper header/footer:  "  Date:        Thu May 08 10:22:01 EDT 2025"
_DATE_RE    = re.compile(r"Date:\s+(.+)")
_COMPLETE_RE = re.compile(r"Pipeline complete:\s+(.+)")
_ACCOUNT_RE = re.compile(r"Account:\s+(\S+)")
_OUTPUT_RE  = re.compile(r"Output:\s+(.+)")
_MODE_RE    = re.compile(r"Mode:\s+(\S+)")
_GENOME_RE  = re.compile(r"Genome:\s+(\S+)")
_SHEET_RE   = re.compile(r"Samplesheet:\s+(.+)")
_DONE_RE    = re.compile(r"Pipeline complete:")

_DATE_FMTS = [
    "%a %b %d %H:%M:%S %Z %Y",   # "Thu May  8 10:22:01 EDT 2025"
    "%a %b  %d %H:%M:%S %Z %Y",  # single-digit day with double space
]

def _try_parse_date(s: str):
    for fmt in _DATE_FMTS:
        try:
            return datetime.strptime(s.strip(), fmt)
        except ValueError:
            pass
    # Fallback parser handles additional date variants and timezone tokens.
    try:
        dt = pd.to_datetime(s.strip(), errors="coerce")
        if pd.isna(dt):
            return None
        return dt.to_pydatetime()
    except Exception:
        pass
    return None


def _parse_nf_log(log_path: Path) -> dict:
    """Extract start time, end time, wall seconds, and SLURM account from nf_*.out."""
    start_dt = None
    end_dt = None
    account = None
    output_dir = None
    completed = False
    mode = None
    genome = None
    samplesheet = None
    for line in log_path.read_text(errors="replace").splitlines():
        m = _DATE_RE.search(line)
        if m:
            dt = _try_parse_date(m.group(1))
            if dt:
                start_dt = dt
        c = _COMPLETE_RE.search(line)
        if c:
            dt = _try_parse_date(c.group(1))
            if dt:
                end_dt = dt
        a = _ACCOUNT_RE.search(line)
        if a:
            account = a.group(1)
        o = _OUTPUT_RE.search(line)
        if o:
            output_dir = o.group(1).strip()
        m = _MODE_RE.search(line)
        if m:
            mode = m.group(1).strip()
        g = _GENOME_RE.search(line)
        if g:
            genome = g.group(1).strip()
        s = _SHEET_RE.search(line)
        if s:
            samplesheet = s.group(1).strip()
        if _DONE_RE.search(line):
            completed = True

    wall_seconds = None
    if start_dt and end_dt:
        wall_seconds = (end_dt - start_dt).total_seconds()
        if wall_seconds < 0:
            wall_seconds = None

    return {
        "account": account,
        "wall_seconds": wall_seconds,
        "output_dir": output_dir,
        "completed": completed,
        "mode": mode,
        "genome": genome,
        "samplesheet": samplesheet,
    }


def _build_global_log_map(global_log_dir: Path) -> dict:
    """Map patient name -> best matching log info from shared logs directory."""
    mapping = {}
    if not global_log_dir.is_dir():
        return mapping

    logs = sorted(global_log_dir.glob("nf_*.out"))
    for log_path in logs:
        info = _parse_nf_log(log_path)
        out = info.get("output_dir")
        if not out:
            continue

        patient = Path(out).name
        prev = mapping.get(patient)

        # Prefer completed runs, then runs with wall-time, then newest filename sort.
        score = (
            2 if info.get("completed") else 0,
            1 if info.get("wall_seconds") else 0,
            str(log_path),
        )
        if prev is None or score > prev["_score"]:
            mapping[patient] = {
                "account": info.get("account"),
                "wall_seconds": info.get("wall_seconds"),
                "log_path": str(log_path),
                "completed": info.get("completed", False),
                "mode": info.get("mode"),
                "genome": info.get("genome"),
                "samplesheet": info.get("samplesheet"),
                "_score": score,
            }

    # Drop internal score before returning.
    for patient in list(mapping.keys()):
        mapping[patient].pop("_score", None)
    return mapping


def _parse_trace(trace_path: Path) -> pd.DataFrame:
    """Return a DataFrame with columns [task_name, status, duration_s, realtime_s]."""
    try:
        df = pd.read_csv(trace_path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame()

    if "name" not in df.columns or "duration" not in df.columns:
        return pd.DataFrame()

    df["duration_s"] = df["duration"].apply(_parse_nf_duration)
    df["realtime_s"] = (
        df["realtime"].apply(_parse_nf_duration)
        if "realtime" in df.columns
        else float("nan")
    )
    # Strip process index e.g. "NFCORE_ONCOANALYSER:SAGE (PATIENT1)" → "NFCORE_ONCOANALYSER:SAGE"
    df["task_name"] = (
        df["name"]
        .str.replace(r"\s*\([^)]*\)$", "", regex=True)
        .str.strip()
    )
    return df[["task_name", "status", "duration_s", "realtime_s"]]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Summarize oncoanalyser cohort wall time and per-task runtimes."
    )
    parser.add_argument("results_base", help="Cohort results base directory")
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory (default: <results_base>/runtime_summary)",
    )
    parser.add_argument(
        "--top-tasks", type=int, default=20,
        help="Number of top tasks to show in bar chart (default: 20)",
    )
    parser.add_argument(
        "--format", choices=["png", "pdf"], default="png",
        help="Plot file format (default: png)",
    )
    parser.add_argument(
        "--global-logs-dir", default=None,
        help="Shared logs directory containing nf_*.out (default: <results_base>/logs)",
    )
    parser.add_argument(
        "--patient-prefix", default="PANCAN-H-",
        help="Only include top-level directories with this prefix (default: PANCAN-H-). Use '' to disable.",
    )
    parser.add_argument(
        "--wrapper-script",
        default="/cluster/projects/pughlab/src/oncoanalyser_pughlab/run_oncoanalyser.sh",
        help="Path to run_oncoanalyser.sh used when printing resume commands.",
    )
    args = parser.parse_args()

    results_base = Path(args.results_base).resolve()
    if not results_base.is_dir():
        sys.exit(f"ERROR: Directory not found: {results_base}")

    output_dir = Path(args.output_dir) if args.output_dir else results_base / "runtime_summary"
    output_dir.mkdir(parents=True, exist_ok=True)

    global_logs_dir = Path(args.global_logs_dir) if args.global_logs_dir else results_base / "logs"
    global_log_map = _build_global_log_map(global_logs_dir)

    fmt = args.format
    top_n = args.top_tasks

    # ── Collect data ──────────────────────────────────────────────────────────
    patient_rows = []
    all_task_frames = []

    for patient_dir in sorted(results_base.iterdir()):
        if not patient_dir.is_dir() or patient_dir.name in {"runtime_summary", "logs", "fingerprint_metrics"}:
            continue
        if patient_dir.name.startswith(".") or patient_dir.name.startswith("._"):
            continue
        if args.patient_prefix and not patient_dir.name.startswith(args.patient_prefix):
            continue

        patient = patient_dir.name
        log_dir = patient_dir / "logs"
        trace_path = patient_dir / "pipeline_info" / "execution_trace.txt"

        # Newest log file for this patient
        wall_s = None
        account = None
        log_source = "none"
        log_path = None
        completed = False
        mode = None
        genome = None
        samplesheet = None
        if log_dir.is_dir():
            logs = sorted(log_dir.glob("nf_*.out"))
            if logs:
                info = _parse_nf_log(logs[-1])
                wall_s = info["wall_seconds"]
                account = info["account"]
                log_source = "patient_logs"
                log_path = str(logs[-1])
                completed = info.get("completed", False)
                mode = info.get("mode")
                genome = info.get("genome")
                samplesheet = info.get("samplesheet")

        # Fallback to shared logs directory, mapped by wrapper Output path basename.
        if wall_s is None:
            g = global_log_map.get(patient)
            if g:
                wall_s = g.get("wall_seconds")
                account = g.get("account")
                log_source = "global_logs"
                log_path = g.get("log_path")
                completed = g.get("completed", False)
                mode = g.get("mode")
                genome = g.get("genome")
                samplesheet = g.get("samplesheet")

        patient_rows.append({
            "patient":    patient,
            "account":    account if account and account != "auto" else "unknown",
            "wall_hours": round(wall_s / 3600, 3) if wall_s is not None else None,
            "log_source": log_source,
            "log_path": log_path,
            "log_completed": completed,
            "mode": mode,
            "genome": genome,
            "samplesheet": samplesheet,
        })

        if trace_path.exists():
            df = _parse_trace(trace_path)
            if not df.empty:
                df["patient"] = patient
                all_task_frames.append(df)

    if not patient_rows:
        sys.exit("No patient directories with log files found.")

    # ── Patient summary ───────────────────────────────────────────────────────
    pat_df = pd.DataFrame(patient_rows)
    pat_df.to_csv(output_dir / "patient_summary.tsv", sep="\t", index=False)

    complete = pat_df["wall_hours"].notna()
    mapped = pat_df["log_source"].ne("none")
    n_complete = complete.sum()
    n_mapped = mapped.sum()
    print(f"\nPatients found:        {len(pat_df)}")
    print(f"Patients with logs:    {n_mapped}")
    print(f"Patients with walltime:{n_complete}")
    if "log_source" in pat_df.columns:
        src_counts = pat_df["log_source"].value_counts().to_string()
        print(f"Log source breakdown:\n{src_counts}\n")

    if n_mapped and n_complete < n_mapped:
        missing_wall = pat_df.loc[mapped & ~complete, ["patient", "log_source", "log_path", "log_completed"]]
        missing_wall.to_csv(output_dir / "patients_missing_walltime.tsv", sep="\t", index=False)
        print(f"Wrote diagnostics: {output_dir / 'patients_missing_walltime.tsv'}\n")

    # ── Incomplete run section with resume commands ─────────────────────────
    incomplete = pat_df.loc[mapped & (~pat_df["log_completed"])].copy()
    if not incomplete.empty:
        resume_rows = []
        print("Incomplete/failed cases (no 'Pipeline complete' marker):")
        for _, row in incomplete.sort_values("patient").iterrows():
            patient = row["patient"]
            mode = row["mode"] if pd.notna(row["mode"]) else None
            genome = row["genome"] if pd.notna(row["genome"]) else None
            samplesheet = row["samplesheet"] if pd.notna(row["samplesheet"]) else None
            account = row["account"] if pd.notna(row["account"]) and row["account"] != "unknown" else "pughlab"
            outdir = str((results_base / patient).resolve())
            resume_args = None

            if samplesheet and mode and genome:
                resume_args = [
                    "sbatch",
                    "--parsable",
                    f"--export=ALL,SOURCE_REPO={Path(args.wrapper_script).resolve().parent},NXF_OFFLINE=true",
                    "-A",
                    account,
                    args.wrapper_script,
                    samplesheet,
                    outdir,
                    mode,
                    genome,
                    account,
                ]
                cmd = "RESUME=true " + " ".join(shlex.quote(x) for x in resume_args)
                resume_ready = True
            else:
                cmd = (
                    f"# Missing mode/genome/samplesheet in log for {patient}.\n"
                    f"# Fill placeholders and run:\n"
                    f"RESUME=true sbatch -A {account} {args.wrapper_script} "
                    f"/path/to/{patient}/samplesheet.csv {outdir} <mode> <genome> {account}"
                )
                resume_ready = False

            print(f"- {patient}")
            print(f"  log: {row['log_path']}")
            print(f"  resume: {cmd}")

            resume_rows.append({
                "patient": patient,
                "log_path": row["log_path"],
                "account": account,
                "mode": mode,
                "genome": genome,
                "samplesheet": samplesheet,
                "resume_ready": resume_ready,
                "resume_command": cmd,
                "resume_args": resume_args,
            })

        pd.DataFrame(resume_rows).to_csv(output_dir / "resume_incomplete_cases.tsv", sep="\t", index=False)
        print(f"\nWrote resume commands: {output_dir / 'resume_incomplete_cases.tsv'}\n")

        # Optional interactive execution of generated resume commands.
        run_ready = [r for r in resume_rows if r["resume_ready"]]
        if run_ready:
            prompt = f"Run resume commands for all [{len(run_ready)}] samples now? (y/n): "
            try:
                answer = input(prompt).strip().lower()
            except EOFError:
                answer = "n"

            if answer == "y":
                print("Submitting resume jobs...")
                submitted = 0
                failed = 0
                submission_rows = []
                for row in run_ready:
                    patient = row["patient"]
                    run_env = dict(os.environ)
                    run_env["RESUME"] = "true"
                    result = subprocess.run(
                        row["resume_args"],
                        text=True,
                        capture_output=True,
                        env=run_env,
                    )
                    if result.returncode == 0:
                        submitted += 1
                        stdout = result.stdout.strip()
                        # For --parsable, first token is job id (sometimes suffixed by ';cluster').
                        job_id = stdout.split(";")[0].splitlines()[-1].strip() if stdout else ""
                        msg = f"job_id={job_id}" if job_id else "submitted (job id not returned)"
                        print(f"  [OK] {patient}: {msg}")
                        submission_rows.append({
                            "patient": patient,
                            "status": "OK",
                            "job_id": job_id,
                            "stdout": result.stdout.strip(),
                            "stderr": result.stderr.strip(),
                            "command": row["resume_command"],
                        })
                    else:
                        failed += 1
                        err = result.stderr.strip() or result.stdout.strip() or "unknown error"
                        print(f"  [FAILED] {patient}: {err}")
                        submission_rows.append({
                            "patient": patient,
                            "status": "FAILED",
                            "job_id": "",
                            "stdout": result.stdout.strip(),
                            "stderr": result.stderr.strip(),
                            "command": row["resume_command"],
                        })

                if submission_rows:
                    sub_out = output_dir / "resume_submission_results.tsv"
                    pd.DataFrame(submission_rows).to_csv(sub_out, sep="\t", index=False)
                    print(f"Wrote submission results: {sub_out}")
                print(f"Resume submissions complete: {submitted} submitted, {failed} failed.")
            else:
                print("Resume submission skipped.")

    if n_complete:
        wh = pat_df.loc[complete, "wall_hours"]
        print(f"Wall time (hours)")
        print(f"  mean:   {wh.mean():.2f}")
        print(f"  median: {wh.median():.2f}")
        print(f"  min:    {wh.min():.2f}")
        print(f"  max:    {wh.max():.2f}")
        print()
        acc_counts = pat_df.loc[complete, "account"].value_counts().to_string()
        print(f"Accounts:\n{acc_counts}\n")

    # ── Task summary ─────────────────────────────────────────────────────────
    task_stats = None
    if all_task_frames:
        tasks_df = pd.concat(all_task_frames, ignore_index=True)
        task_stats = (
            tasks_df.groupby("task_name")["duration_s"]
            .agg(n_runs="count", mean_s="mean", median_s="median", max_s="max")
            .sort_values("mean_s", ascending=False)
            .reset_index()
        )
        for col in ("mean_s", "median_s", "max_s"):
            task_stats[col.replace("_s", "_h")] = task_stats[col] / 3600
        task_stats.to_csv(output_dir / "task_stats.tsv", sep="\t", index=False)

    # ── Plots ─────────────────────────────────────────────────────────────────
    wh_vals = pat_df["wall_hours"].dropna()

    # 1. Patient wall time histogram
    if not wh_vals.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(wh_vals, bins="auto", edgecolor="white", color="#4C72B0")
        ax.set_xlabel("Wall Time (hours)")
        ax.set_ylabel("Number of Patients")
        ax.set_title(f"Patient Wall Time Distribution (n={len(wh_vals)})")
        fig.tight_layout()
        fig.savefig(output_dir / f"patient_wall_time_hist.{fmt}", dpi=150)
        plt.close(fig)

    # 2. Wall time by SLURM account (box plot, only if >1 account)
    accounts = pat_df.loc[complete, "account"].unique()
    if len(accounts) > 1 and not wh_vals.empty:
        groups = [
            pat_df.loc[(pat_df["account"] == a) & complete, "wall_hours"].dropna().tolist()
            for a in accounts
        ]
        fig, ax = plt.subplots(figsize=(max(6, len(accounts) * 1.8), 4))
        ax.boxplot(groups, tick_labels=accounts)
        ax.set_xlabel("SLURM Account")
        ax.set_ylabel("Wall Time (hours)")
        ax.set_title("Wall Time by SLURM Account")
        fig.tight_layout()
        fig.savefig(output_dir / f"wall_time_by_account.{fmt}", dpi=150)
        plt.close(fig)

    # 3. Top N tasks by mean duration
    if task_stats is not None and not task_stats.empty:
        top = task_stats.head(top_n)
        fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.4)))
        colors = ["#DD8452" if h > 1 else "#4C72B0" for h in top["mean_h"]]
        ax.barh(top["task_name"][::-1], top["mean_h"][::-1], color=colors[::-1])
        ax.set_xlabel("Mean Duration (hours)")
        ax.set_title(f"Top {len(top)} Tasks by Mean Duration")
        fig.tight_layout()
        fig.savefig(output_dir / f"top_tasks_mean_duration.{fmt}", dpi=150)
        plt.close(fig)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"Outputs written to: {output_dir}/")
    print("  patient_summary.tsv")
    if task_stats is not None:
        print("  task_stats.tsv")
    print(f"  *.{fmt}  (plots)")


if __name__ == "__main__":
    main()
