#!/usr/bin/env python3
"""
Build the vendored tumour antigen annotation table used by compile_mmr_rna_dashboard.py.

Cluster nodes are offline, so this is run on a machine with internet access and the
resulting TSV is committed to resources/antigens/.

Sources:
  * IEDB (https://www.iedb.org) T-cell assays: human source protein, human host,
    positive outcome, disease under DOID:162 (cancer). Aggregated per UniProt accession.
  * Human Protein Atlas RNA tissue specificity: genes "Tissue enriched" or "Group enriched"
    in testis, plus canonical cancer-testis gene families (CTdatabase nomenclature; HPA calls
    some of these, e.g. CTAG2, only "Tissue enhanced"), used as the cancer-testis class.
  * UniProt REST: accession -> primary gene symbol.
  * HGNC protein-coding gene set (approved, alias and previous symbols), written to
    protein_coding_genes.txt beside the antigen table so outlier views can exclude non-coding loci.

Usage:
  python3 scripts/build_tumor_antigen_reference.py \
      --output resources/antigens/tumor_antigen_annotation.tsv
"""

import sys
import csv
import re
import json
import argparse
import datetime
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

IEDB_URL = "https://query-api.iedb.org/tcell_search"
HPA_URL = "https://www.proteinatlas.org/api/search_download.php"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/accessions"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_groups/protein-coding_gene.txt"
PAGE_SIZE = 5000

# Canonical cancer-testis antigen families (CTdatabase nomenclature)
CT_FAMILY_PATTERN = re.compile(
    r"^(MAGE[ABC]\d+[A-Z]?|GAGE\d+[A-Z]?|PAGE\d+|XAGE\d+[A-Z]?|SSX\d+|CTAG1[AB]|CTAG2|CT45A\d+|CT47[AB]\d*|"
    r"SPANX[A-Z]\d*|CSAG\d|SAGE1|NXF2B?|BAGE\d*|CTCFL|SYCP1|BRDT|PRAME|DDX43|CT83|TPTE|ADAM2|SPA17|"
    r"CXorf67|TEX15|TDRD1|PIWIL[12]|LUZP4|CTNNA2|DSCR8|POTE[A-Z]|VCX\d*[A-Z]?|VCY\d*[A-Z]?)$"
)


def parse_args():
    parser = argparse.ArgumentParser(description="Build tumour antigen annotation TSV from IEDB + HPA.")
    parser.add_argument(
        "--output",
        default=str(Path(__file__).resolve().parent.parent / "resources" / "antigens" / "tumor_antigen_annotation.tsv"),
        help="Destination TSV path."
    )
    return parser.parse_args()


def http_get(url, headers=None, timeout=120):
    req = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


def fetch_iedb_cancer_tcell():
    """Return {uniprot_acc: {"name", "assays", "refs", "diseases"}} for positive human cancer T-cell assays."""
    params = {
        "select": "parent_source_antigen_name,parent_source_antigen_iri,pubmed_id,reference_id,disease_names,mhc_class",
        "source_organism_iri": "eq.NCBITaxon:9606",
        "host_organism_iri": "eq.NCBITaxon:9606",
        "disease_iri_search": "cs.{DOID:162}",
        "qualitative_measure": "like.Positive*",
        "order": "tcell_id",
    }
    base = IEDB_URL + "?" + urllib.parse.urlencode(params, safe="{}:*.,")
    antigens = defaultdict(lambda: {"name": "", "assays": 0, "refs": set(), "diseases": defaultdict(int), "mhc": set()})

    offset = 0
    while True:
        rows = json.loads(http_get(base, headers={"Range": f"{offset}-{offset + PAGE_SIZE - 1}"}))
        for r in rows:
            iri = r.get("parent_source_antigen_iri") or ""
            if not iri.startswith("UNIPROT:"):
                continue
            acc = iri.split(":", 1)[1]
            rec = antigens[acc]
            rec["name"] = (r.get("parent_source_antigen_name") or "").split(" (UniProt")[0]
            rec["assays"] += 1
            if r.get("reference_id"):
                rec["refs"].add(r["reference_id"])
            if r.get("mhc_class"):
                rec["mhc"].add(r["mhc_class"])
            for d in r.get("disease_names") or []:
                rec["diseases"][d] += 1
        print(f"  IEDB: fetched {offset + len(rows)} assays", file=sys.stderr)
        if len(rows) < PAGE_SIZE:
            break
        offset += PAGE_SIZE
    return antigens


def map_uniprot_to_gene(accessions):
    """Map UniProt accessions to primary gene symbols (UniProt REST, batched).

    Proteins encoded by duplicated loci (e.g. NY-ESO-1 = CTAG1A; CTAG1B) map to every gene.
    """
    mapping = {}
    accs = sorted(accessions)
    for i in range(0, len(accs), 500):
        chunk = accs[i:i + 500]
        q = urllib.parse.urlencode({"accessions": ",".join(chunk), "fields": "accession,gene_primary", "format": "tsv"})
        text = http_get(f"{UNIPROT_URL}?{q}")
        for row in csv.DictReader(text.splitlines(), delimiter="\t"):
            genes = [g.strip() for g in (row.get("Gene Names (primary)") or "").split(";") if g.strip()]
            if genes:
                mapping[row["Entry"]] = genes
    return mapping


def fetch_hpa_testis_enriched():
    genes = set()
    for category in ("Tissue enriched", "Group enriched"):
        q = urllib.parse.urlencode({
            "search": f"tissue_category_rna:testis;{category}",
            "format": "tsv",
            "columns": "g,eg",
            "compress": "no",
        })
        text = http_get(f"{HPA_URL}?{q}")
        genes |= {row["Gene"].strip() for row in csv.DictReader(text.splitlines(), delimiter="\t") if row.get("Gene")}
    return genes


def fetch_hgnc_protein_coding():
    symbols = set()
    for row in csv.DictReader(http_get(HGNC_URL, timeout=300).splitlines(), delimiter="\t"):
        if row.get("status") != "Approved":
            continue
        symbols.add(row["symbol"])
        for col in ("alias_symbol", "prev_symbol"):
            symbols.update(x.strip() for x in (row.get(col) or "").split("|") if x.strip())
    return symbols


def main():
    args = parse_args()

    print("Fetching IEDB cancer T-cell epitope assays...", file=sys.stderr)
    iedb = fetch_iedb_cancer_tcell()
    print(f"  {len(iedb)} distinct UniProt source antigens", file=sys.stderr)

    print("Mapping UniProt accessions to gene symbols...", file=sys.stderr)
    acc_to_gene = map_uniprot_to_gene(iedb.keys())

    print("Fetching HPA testis-enriched (tissue + group) genes...", file=sys.stderr)
    testis = fetch_hpa_testis_enriched()
    print(f"  {len(testis)} testis-enriched genes", file=sys.stderr)

    # Collapse to gene level (several accessions/isoforms can share a gene)
    genes = {}
    for acc, rec in iedb.items():
        for gene in acc_to_gene.get(acc, []):
            g = genes.setdefault(gene, {"accs": set(), "name": rec["name"], "assays": 0, "refs": set(),
                                        "diseases": defaultdict(int), "mhc": set()})
            g["accs"].add(acc)
            g["assays"] += rec["assays"]
            g["refs"] |= rec["refs"]
            g["mhc"] |= rec["mhc"]
            for d, n in rec["diseases"].items():
                g["diseases"][d] += n

    # IEDB antigens from canonical CT families count as cancer-testis even when HPA is less strict
    testis |= {g for g in genes if CT_FAMILY_PATTERN.match(g)}

    for gene in testis:
        genes.setdefault(gene, {"accs": set(), "name": "", "assays": 0, "refs": set(),
                                "diseases": defaultdict(int), "mhc": set()})

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        f.write(f"# Built {datetime.date.today().isoformat()} by scripts/build_tumor_antigen_reference.py\n")
        f.write("# Sources: IEDB T-cell assays (human, positive, DOID:162 cancer); HPA testis tissue/group enriched + CT gene families; UniProt gene mapping\n")
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "antigen_class", "iedb_tcell_assays", "iedb_references", "mhc_classes",
                    "testis_restricted", "top_diseases", "uniprot", "protein_name"])
        for gene in sorted(genes):
            g = genes[gene]
            in_iedb = g["assays"] > 0
            in_testis = gene in testis
            if in_iedb and in_testis:
                cls = "Cancer-testis T-cell antigen"
            elif in_iedb:
                cls = "T-cell antigen"
            else:
                cls = "Cancer-testis-like"
            top_dx = sorted(g["diseases"].items(), key=lambda x: -x[1])[:3]
            w.writerow([
                gene, cls, g["assays"], len(g["refs"]), ",".join(sorted(g["mhc"])),
                "yes" if in_testis else "no",
                "; ".join(d for d, _ in top_dx),
                ",".join(sorted(g["accs"])), g["name"],
            ])
    print(f"Wrote {len(genes)} genes to {out_path}", file=sys.stderr)

    print("Fetching HGNC protein-coding gene symbols...", file=sys.stderr)
    coding = fetch_hgnc_protein_coding()
    coding_path = out_path.parent / "protein_coding_genes.txt"
    with open(coding_path, "w", encoding="utf-8") as f:
        f.write(f"# Built {datetime.date.today().isoformat()} from HGNC protein-coding_gene.txt (approved + alias + previous symbols)\n")
        f.write("\n".join(sorted(coding)) + "\n")
    print(f"Wrote {len(coding)} protein-coding symbols to {coding_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
