# Oncoanalyser Analysis Scripts

Reusable Python and Bash utilities for cohort-level analysis of oncoanalyser outputs.

## Structure

- `analysis/` — analysis and comparison scripts
- `shared/` — shared utilities and helpers

## Available Scripts

### Fingerprint Comparison

**Script:** `analysis/compare_fingerprints.sh`

Compare germline variants across samples in a cohort to identify related or identical individuals.

Uses GATK4 `GenotypeConcordance` for all-vs-all sample comparisons on SAGE germline VCFs (including both tumor and normal samples in each VCF).

```bash
bash scripts/analysis/compare_fingerprints.sh /path/to/cohort/results
```

Outputs:

- `sample_inventory.tsv`: discovered samples and source VCFs
- `pairwise_concordance.tsv`: all-vs-all comparison run status
- `pairs/*.genotype_concordance_*_metrics`: GATK summary/detail/contingency metrics for each sample pair

Stepwise testing examples:

```bash
# 1) Build inventory only (no GATK pair runs)
bash scripts/analysis/compare_fingerprints.sh /path/to/cohort --discover-only

# 2) Test only first 5 pairwise comparisons
bash scripts/analysis/compare_fingerprints.sh /path/to/cohort --pair-limit 5

# 3) Use local env file (recommended for development)
cp scripts/dev/compare_fingerprints.env.example scripts/dev/compare_fingerprints.env
bash scripts/analysis/compare_fingerprints.sh --env-file scripts/dev/compare_fingerprints.env
```

Recommended development pattern:

- Keep local settings in `scripts/dev/*.env` (ignored by git)
- Commit only `scripts/dev/*.env.example`
- Prefer `--discover-only` and `--pair-limit` when debugging

See `analysis/compare_fingerprints.sh --help` for full usage and output details.

### Runtime Summary

**Script:** `analysis/summarize_runtime.py`

Summarizes wall time and per-task runtimes across a cohort, broken down by SLURM account.

Reads the Nextflow SLURM wrapper logs (`logs/nf_*.out`) for job-level wall time and account, and
`pipeline_info/execution_trace.txt` for per-task timing.

```bash
module load python3/3.12.11
pip install pandas matplotlib  # one-time

python3 scripts/analysis/summarize_runtime.py /path/to/results_base
```

Outputs (written to `<results_base>/runtime_summary/` by default):

- `patient_summary.tsv`: per-patient wall time and account
- `task_stats.tsv`: per-task mean/median/max duration across cohort
- `patient_wall_time_hist.png`: wall time distribution
- `wall_time_by_account.png`: wall time grouped by SLURM account (if >1 account)
- `top_tasks_mean_duration.png`: top N tasks by mean duration

Options: `--output-dir`, `--top-tasks N`, `--format png|pdf`

## Adding New Scripts

1. Create a new file in `analysis/` (bash or Python)
2. Add a header comment block explaining:
   - Purpose
   - Input assumptions (paths, file formats)
   - Output structure
   - Required modules/dependencies
3. Use utilities from `shared/` for common tasks (path handling, logging, etc.)
4. Update this README with a brief description

## Notes

- Scripts assume oncoanalyser output structure (as documented in main README)
- All paths are relative to the cohort results directory unless otherwise noted
- Bash scripts should be idempotent where practical (safe to re-run)
