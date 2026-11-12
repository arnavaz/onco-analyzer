#!/usr/bin/env python3
"""
HLA class I binding model, re-implemented from the Hartwig Neo trained model files.

The reference bundle ships the trained model under
  <hmf_pipeline_resources>/misc/neo/binding/
    neo_train_pos_weight.csv         per allele x peptide length: position weight matrix
    neo_train_flank_pos_weight.csv   amino-acid weights for the 3 up / 3 down flank residues
    neo_train_rand_dist.csv          score -> rank percentile, per allele x length
    neo_train_likelihood.csv         rank percentile -> binding likelihood, per allele x length
    neo_train_expression_dist.csv    TPM bucket -> observed binding rate
    neo_train_recognition.csv        peptides with measured T-cell recognition, with allele

Scoring one peptide against one allele is a sum of position weights, converted to a rank
percentile against the allele's random-peptide distribution (lower rank = better binder;
<= 2% is the conventional binder threshold, <= 0.5% strong). This mirrors
com.hartwig.hmftools.neo.bind.BindScorer so it can run offline with no Java, and
`validate_against_jar` in predict_neoepitope_binding.py checks it against the real scorer.
"""

import csv
import gzip
import math
import re
from bisect import bisect_left
from pathlib import Path

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MIN_PEPTIDE_LENGTH = 8
MAX_PEPTIDE_LENGTH = 12
BINDER_RANK = 2.0
STRONG_BINDER_RANK = 0.5

# A*24:02 / HLA-A*24:02 / A24:02 / A2402 all mean allele key A2402
ALLELE_RE = re.compile(r"^(?:HLA[-_]?)?([ABC])\*?(\d{2}):?(\d{2,3})")


def allele_key(text):
    """Normalise an HLA class I allele to the Hartwig key (A*24:02 -> A2402)."""
    m = ALLELE_RE.match(str(text).strip())
    if not m:
        return None
    return f"{m.group(1)}{m.group(2)}{m.group(3)}"


def allele_label(key):
    """A2402 -> A*24:02"""
    return f"{key[0]}*{key[1:3]}:{key[3:]}" if key and len(key) >= 5 else key


def _open(path):
    path = str(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


class BindingModel:
    def __init__(self, binding_dir):
        d = Path(binding_dir)
        self.pos_weights = {}      # (allele, length) -> {aa: [w per position]}
        self.flank_weights = {}    # aa -> [U3, U2, U1, D1, D2, D3]
        self.rand_dist = {}        # (allele, length) -> ([scores desc], [rank percentiles asc])
        self.likelihood = {}       # (allele, length) -> ([ranks asc], [likelihoods])
        self.expression_dist = []  # [(tpm bucket, binding rate)]
        self.recognition = {}      # (peptide, allele) -> source string
        self._load_pos_weights(d / "neo_train_pos_weight.csv")
        self._load_flank_weights(d / "neo_train_flank_pos_weight.csv")
        self._load_rand_dist(d / "neo_train_rand_dist.csv")
        self._load_likelihood(d / "neo_train_likelihood.csv")
        self._load_expression_dist(d / "neo_train_expression_dist.csv")
        self._load_recognition(d / "neo_train_recognition.csv")

    # ---------- loading ----------

    def _load_pos_weights(self, path):
        with _open(path) as f:
            for row in csv.DictReader(f):
                key = (row["Allele"], int(row["PeptideLength"]))
                weights = [float(row[f"P{i}"]) for i in range(12) if f"P{i}" in row]
                self.pos_weights.setdefault(key, {})[row["AminoAcid"]] = weights

    def _load_flank_weights(self, path):
        with _open(path) as f:
            for row in csv.DictReader(f):
                self.flank_weights[row["AminoAcid"]] = [
                    float(row[c]) for c in ("U3", "U2", "U1", "D1", "D2", "D3")
                ]

    def _load_rand_dist(self, path):
        # Rows are grouped by allele+length, ordered by increasing rank bucket
        # (ScoreBucket) and so by decreasing score.
        raw = {}
        with _open(path) as f:
            for row in csv.DictReader(f):
                key = (row["Allele"], int(row["PeptideLength"]))
                raw.setdefault(key, []).append((float(row["Score"]), float(row["ScoreBucket"]) * 100.0))
        for key, pairs in raw.items():
            pairs.sort(key=lambda x: -x[0])
            self.rand_dist[key] = ([p[0] for p in pairs], [p[1] for p in pairs])

    def _load_likelihood(self, path):
        with _open(path) as f:
            reader = csv.DictReader(f)
            rank_cols = [(c, float(c[2:]) * 100.0) for c in (reader.fieldnames or []) if c.startswith("P_")]
            rank_cols.sort(key=lambda x: x[1])
            for row in reader:
                key = (row["Allele"], int(row["PeptideLength"]))
                self.likelihood[key] = (
                    [r for _, r in rank_cols],
                    [float(row[c]) for c, _ in rank_cols],
                )

    def _load_expression_dist(self, path):
        with _open(path) as f:
            for row in csv.DictReader(f):
                if row.get("DataType") == "TpmLevel":
                    self.expression_dist.append((float(row["Bucket"]), float(row["BindingRate"])))
        self.expression_dist.sort()

    def _load_recognition(self, path):
        if not Path(path).is_file():
            return
        with _open(path) as f:
            for row in csv.DictReader(f):
                if str(row.get("Immunogenic", "")).upper() != "TRUE":
                    continue
                key = (row["Peptide"].strip(), row["Allele"].strip())
                self.recognition[key] = (row.get("Source") or "").strip()

    # ---------- scoring ----------

    @property
    def alleles(self):
        return sorted({a for a, _ in self.pos_weights})

    def supports(self, allele, length):
        return (allele, length) in self.pos_weights

    def raw_score(self, allele, peptide):
        """Position-weight score, or None when the allele/length is not in the model."""
        matrix = self.pos_weights.get((allele, len(peptide)))
        if matrix is None:
            return None
        total = 0.0
        for i, aa in enumerate(peptide):
            weights = matrix.get(aa)
            if weights is None:  # X, U, *, ...
                return None
            total += weights[i]
        return total

    def flank_score(self, up_flank, down_flank):
        """Weight contributed by up to 3 residues either side of the peptide."""
        total = 0.0
        for offset, aa in enumerate(reversed((up_flank or "")[-3:])):
            w = self.flank_weights.get(aa)
            if w:
                total += w[2 - offset]
        for offset, aa in enumerate((down_flank or "")[:3]):
            w = self.flank_weights.get(aa)
            if w:
                total += w[3 + offset]
        return total

    def score_rank(self, allele, length, score):
        """Percentile rank of a score against the allele's random peptide distribution."""
        dist = self.rand_dist.get((allele, length))
        if dist is None or score is None:
            return None
        scores, ranks = dist
        # scores descend; find the first entry at or below this score
        lo, hi = 0, len(scores)
        while lo < hi:
            mid = (lo + hi) // 2
            if scores[mid] > score:
                lo = mid + 1
            else:
                hi = mid
        if lo == 0:
            return ranks[0]
        if lo >= len(scores):
            return ranks[-1]
        s_hi, s_lo = scores[lo - 1], scores[lo]
        r_lo, r_hi = ranks[lo - 1], ranks[lo]
        if s_hi == s_lo:
            return r_hi
        return r_lo + (r_hi - r_lo) * (s_hi - score) / (s_hi - s_lo)

    def rank_likelihood(self, allele, length, rank):
        """Trained likelihood that a peptide at this rank is presented."""
        entry = self.likelihood.get((allele, length))
        if entry is None or rank is None:
            return None
        ranks, values = entry
        if rank <= ranks[0]:
            return values[0]
        if rank >= ranks[-1]:
            return values[-1]
        i = bisect_left(ranks, rank)
        r0, r1 = ranks[i - 1], ranks[i]
        v0, v1 = values[i - 1], values[i]
        if r1 == r0:
            return v1
        # likelihoods fall over orders of magnitude, so interpolate in log space
        if v0 > 0 and v1 > 0:
            frac = (rank - r0) / (r1 - r0)
            return math.exp(math.log(v0) + frac * (math.log(v1) - math.log(v0)))
        return v0 + (rank - r0) / (r1 - r0) * (v1 - v0)

    def expression_rate(self, tpm):
        """Observed binding rate for peptides from genes at this expression level."""
        if tpm is None or not self.expression_dist:
            return None
        rate = self.expression_dist[0][1]
        for bucket, value in self.expression_dist:
            if tpm < bucket:
                break
            rate = value
        return rate

    def score(self, allele, peptide, up_flank="", down_flank="", tpm=None):
        """Full result for one peptide/allele, or None when it cannot be scored."""
        raw = self.raw_score(allele, peptide)
        if raw is None:
            return None
        flank = self.flank_score(up_flank, down_flank)
        rank = self.score_rank(allele, len(peptide), raw)
        return {
            "allele": allele,
            "peptide": peptide,
            "score": raw,
            "flank_score": flank,
            "rank": rank,
            "likelihood": self.rank_likelihood(allele, len(peptide), rank),
            "expression_rate": self.expression_rate(tpm),
            "recognition": self.recognition.get((peptide, allele), ""),
        }

    def best_allele(self, alleles, peptide, up_flank="", down_flank="", tpm=None):
        """Lowest-rank (best binding) result across a patient's alleles."""
        results = [self.score(a, peptide, up_flank, down_flank, tpm) for a in alleles]
        results = [r for r in results if r and r["rank"] is not None]
        return min(results, key=lambda r: r["rank"]) if results else None


def tile_peptides(sequence, lengths=range(MIN_PEPTIDE_LENGTH, MAX_PEPTIDE_LENGTH + 1), flank=3):
    """Yield (peptide, up_flank, down_flank, start) for every window of each length."""
    n = len(sequence)
    for length in lengths:
        for start in range(0, n - length + 1):
            end = start + length
            yield (
                sequence[start:end],
                sequence[max(0, start - flank):start],
                sequence[end:end + flank],
                start,
            )


def junction_peptides(sequence, junction_index, lengths=range(MIN_PEPTIDE_LENGTH, MAX_PEPTIDE_LENGTH + 1), flank=3):
    """Windows that span a junction: they must cover residues on both sides.

    junction_index is the index of the first residue contributed by the downstream
    partner, so a window qualifies when start < junction_index < start + length.
    """
    for peptide, up_flank, down_flank, start in tile_peptides(sequence, lengths, flank):
        if start < junction_index < start + len(peptide):
            yield peptide, up_flank, down_flank, start
