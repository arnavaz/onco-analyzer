#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
import urllib.request
import json
import time
import ssl

ssl_context = ssl._create_unverified_context()

# Standard 12 MMR/driver genes to target by default
DEFAULT_GENES = ["MUTYH", "PMS2", "TP53", "MLH1", "MSH2", "MSH6", "EPCAM", "MLH3", "POLD1", "POLE", "APC", "BRAF"]

def ensembl_lookup_gene(gene_symbol, assembly="GRCh38"):
    """Lookup gene coordinates from Ensembl REST API."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, context=ssl_context) as response:
            data = json.loads(response.read().decode())
            chrom = str(data.get("seq_region_name"))
            start = int(data.get("start"))
            end = int(data.get("end"))
            return chrom, start, end
    except Exception as e:
        print(f"Warning: Failed to look up coordinates for gene {gene_symbol} via Ensembl: {e}", file=sys.stderr)
        return None

def _query_gnomad_variants_region_single(chrom, start, end, genome_version="GRCh38"):
    """Query gnomAD region API for a single region interval."""
    url = "https://gnomad.broadinstitute.org/api"
    query = """
    query($chrom: String!, $start: Int!, $stop: Int!, $dataset: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
      region(chrom: $chrom, start: $start, stop: $stop, reference_genome: $referenceGenome) {
        variants(dataset: $dataset) {
          variant_id
          exome { af }
          genome { af }
        }
      }
    }
    """
    dataset = "gnomad_r4" if genome_version == "GRCh38" else "gnomad_r2_1"
    variables = {
        "chrom": str(chrom),
        "start": int(start),
        "stop": int(end),
        "dataset": dataset,
        "referenceGenome": genome_version
    }
    req_data = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    req = urllib.request.Request(url, data=req_data, headers={"Content-Type": "application/json"})
    
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, context=ssl_context) as response:
                res_data = json.loads(response.read().decode())
                if "errors" in res_data:
                    raise Exception(f"GraphQL Errors: {res_data['errors']}")
                res = res_data.get("data", {}).get("region", {}).get("variants", [])
                return res if res is not None else []
        except KeyboardInterrupt:
            print("\nOperation cancelled by user.", file=sys.stderr)
            sys.exit(130)
        except Exception as e:
            print(f"Warning: gnomAD region query attempt {attempt+1} failed for {chrom}:{start}-{end}: {e}", file=sys.stderr)
            backoff = 5 * (attempt + 1)
            if "429" in str(e):
                backoff = 30 * (attempt + 1)
                print(f"  Rate limit (429) hit. Backing off for {backoff} seconds...", file=sys.stderr)
            time.sleep(backoff)
    print(f"Warning: Skipping region query for {chrom}:{start}-{end} after 6 failed attempts.", file=sys.stderr)
    return []

def query_gnomad_variants_region(chrom, start, end, genome_version="GRCh38", max_chunk_size=50000):
    """Query gnomAD region API, automatically chunking large regions into max_chunk_size intervals."""
    if (end - start) > max_chunk_size:
        all_variants = []
        curr = start
        while curr < end:
            chunk_end = min(curr + max_chunk_size, end)
            chunk_vars = _query_gnomad_variants_region_single(chrom, curr, chunk_end, genome_version)
            all_variants.extend(chunk_vars)
            curr = chunk_end
            time.sleep(0.3)
        return all_variants
    else:
        return _query_gnomad_variants_region_single(chrom, start, end, genome_version)

def query_gnomad_variants_gene(gene_symbol, genome_version="GRCh38"):
    """Query gnomAD gene API for all variants in a specific gene."""
    url = "https://gnomad.broadinstitute.org/api"
    query = """
    query($geneSymbol: String!, $dataset: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
      gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
        variants(dataset: $dataset) {
          variant_id
          exome { af }
          genome { af }
        }
      }
    }
    """
    dataset = "gnomad_r4" if genome_version == "GRCh38" else "gnomad_r2_1"
    variables = {
        "geneSymbol": gene_symbol,
        "dataset": dataset,
        "referenceGenome": genome_version
    }
    req_data = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    req = urllib.request.Request(url, data=req_data, headers={"Content-Type": "application/json"})
    
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, context=ssl_context) as response:
                res_data = json.loads(response.read().decode())
                if "errors" in res_data:
                    raise Exception(f"GraphQL Errors: {res_data['errors']}")
                res = res_data.get("data", {}).get("gene", {}).get("variants", [])
                return res if res is not None else []
        except Exception as e:
            print(f"Warning: gnomAD gene query attempt {attempt+1} failed for {gene_symbol}: {e}", file=sys.stderr)
            backoff = 5 * (attempt + 1)
            if "429" in str(e):
                backoff = 30 * (attempt + 1)
                print(f"  Rate limit (429) hit. Backing off for {backoff} seconds...", file=sys.stderr)
            time.sleep(backoff)
    print(f"Warning: Skipping gene query for {gene_symbol} after 6 failed attempts.", file=sys.stderr)
    return []

def merge_intervals(intervals, max_gap=10000):
    """Merge intervals that are within max_gap of each other to minimize region queries."""
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start - prev_end <= max_gap:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return merged

def main():
    parser = argparse.ArgumentParser(description="Expand AmberGermlineSites list with high-frequency panel-targeted SNPs")
    parser.add_argument("--target-regions-bed", required=True, help="Path to panel targeted regions BED file (.bed or .bed.gz)")
    parser.add_argument("--default-loci", required=True, help="Path to standard default AmberGermlineSites.tsv.gz file")
    parser.add_argument("--output", required=True, help="Path to write the output expanded, gzipped loci file")
    parser.add_argument("--genome", default="GRCh38_hmf", help="Genome version (GRCh38_hmf or GRCh37_hmf)")
    parser.add_argument("--maf", type=float, default=0.001, help="Minimum global Minor Allele Frequency (MAF) threshold to include (default: 0.001)")
    parser.add_argument("--genes", default=",".join(DEFAULT_GENES), help="Comma-separated list of genes to expand (or 'all_bed' to fetch SNPs for every bed region, warning: very slow)")
    
    parser.add_argument("--no-cache", "--force", action="store_true", help="Ignore any existing checkpoint cache file and re-query gnomAD")
    
    args = parser.parse_args()
    
    is_grch37 = "37" in args.genome
    ref_genome = "GRCh37" if is_grch37 else "GRCh38"
    
    # 1. Load BED regions
    print("Loading targeted panel BED intervals...")
    bed_regions = {}
    opener = gzip.open if args.target_regions_bed.endswith(".gz") else open
    with opener(args.target_regions_bed, "rt") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom = parts[0].replace("chr", "")
                start, end = int(parts[1]), int(parts[2])
                if chrom not in bed_regions:
                    bed_regions[chrom] = []
                bed_regions[chrom].append((start, end))
                
    # 2. Load standard heterozygous sites first
    print(f"Loading default loci file from {args.default_loci}...")
    default_df_rows = []
    seen_loci = set()
    with gzip.open(args.default_loci, "rt") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                default_df_rows.append(parts)
                
    # Parse header columns mapping
    c_idx = header.index("Chromosome")
    p_idx = header.index("Position")
    r_idx = header.index("Ref")
    a_idx = header.index("Alt")
    s_idx = header.index("SnpCheck") if "SnpCheck" in header else -1
    
    print("Populating seen loci from default set...")
    for row in default_df_rows:
        chrom = row[c_idx]
        pos = int(row[p_idx])
        chrom_norm = f"chr{chrom.replace('chr', '')}"
        seen_loci.add((chrom_norm, pos))

    # 3. Retrieve variants from gnomAD and directly append new ones to default list
    final_rows = list(default_df_rows)
    extra_snps_count = 0
    
    if args.genes.lower() == "auto":
        print("Auto-detecting genes from the 4th column of the BED file...")
        detected_genes = set()
        opener = gzip.open if args.target_regions_bed.endswith(".gz") else open
        with opener(args.target_regions_bed, "rt") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split("\t")
                    if len(parts) > 3:
                        # Extract gene name, e.g. "MSH2" or "MSH2_exon1"
                        gene_name = parts[3].split("_")[0].split("-")[0].split(".")[0].strip()
                        if gene_name and gene_name.isalnum() and not gene_name.isdigit():
                            detected_genes.add(gene_name)
        if detected_genes:
            print(f"Auto-detected {len(detected_genes)} genes: {sorted(list(detected_genes))}")
            args.genes = ",".join(sorted(list(detected_genes)))
        else:
            print("Warning: No genes detected in the 4th column. Defaulting to 'all_bed' approach.")
            args.genes = "all_bed"
            
    # Checkpoint / Cache setup for graceful restart
    cache_file = os.path.abspath(args.output) + ".cache.json"
    cached_queries = set()
    cached_snps_data = []
    if os.path.exists(cache_file) and not args.no_cache:
        try:
            with open(cache_file, "r") as cf:
                cdata = json.load(cf)
                cached_maf = cdata.get("maf")
                if cached_maf is not None and float(cached_maf) != float(args.maf):
                    print(f"Notice: Cache MAF ({cached_maf}) does not match current MAF ({args.maf}). Invalidating cache and re-querying.")
                    cached_queries = set()
                    cached_snps_data = []
                else:
                    cached_queries = set(cdata.get("completed_queries", []))
                    cached_snps_data = cdata.get("cached_snps", [])
                    print(f"Resuming from checkpoint file ({cache_file}): {len(cached_queries)} queries previously completed, {len(cached_snps_data)} SNPs loaded.")
                    for row in cached_snps_data:
                        chrom_norm = f"chr{row[c_idx].replace('chr', '')}"
                        pos = int(row[p_idx])
                        seen_loci.add((chrom_norm, pos))
                        final_rows.append(row)
                        extra_snps_count += 1
        except Exception as e:
            print(f"Warning: Failed to load checkpoint file {cache_file}: {e}. Starting fresh.", file=sys.stderr)
            cached_queries = set()
            cached_snps_data = []

    def save_checkpoint(query_key, new_snps_rows):
        cached_queries.add(query_key)
        cached_snps_data.extend(new_snps_rows)
        try:
            with open(cache_file, "w") as cf:
                json.dump({
                    "maf": args.maf,
                    "completed_queries": list(cached_queries),
                    "cached_snps": cached_snps_data
                }, cf)
        except Exception as e:
            print(f"Warning: Failed to save checkpoint: {e}", file=sys.stderr)

    if args.genes.lower() == "all_bed":
        print("Using all BED intervals to query variations from gnomAD...")
        # Merge overlapping/nearby intervals to minimize API calls
        query_regions = []
        for chrom, intervals in bed_regions.items():
            merged = merge_intervals(intervals, max_gap=10000)
            for start, end in merged:
                query_regions.append((chrom, start, end))
                
        print(f"Merged {sum(len(v) for v in bed_regions.values())} BED intervals into {len(query_regions)} query regions.")
        
        for idx, (chrom, start, end) in enumerate(query_regions):
            q_key = f"region_{chrom}:{start}-{end}"
            if q_key in cached_queries:
                print(f"[{idx+1}/{len(query_regions)}] Region ({chrom}:{start}-{end}) already completed (cached). Skipping.")
                continue

            print(f"[{idx+1}/{len(query_regions)}] Querying region ({chrom}:{start}-{end})...")
            variants = query_gnomad_variants_region(chrom, start, end, genome_version=ref_genome)
            time.sleep(0.5) # rate limit politeness
            
            region_added_count = 0
            new_region_rows = []
            for v in variants:
                vid = v.get("variant_id")
                if not vid: continue
                v_chrom, pos_str, ref, alt = vid.split("-")
                pos = int(pos_str)
                if len(ref) != 1 or len(alt) != 1: continue
                
                # Check if variant falls inside any of the original BED intervals for this chromosome
                v_chrom_norm = v_chrom.replace("chr", "")
                if v_chrom_norm not in bed_regions: continue
                if not any(b_start <= pos <= b_end for b_start, b_end in bed_regions[v_chrom_norm]):
                    continue
                
                exome_af = v.get("exome", {}).get("af", 0) if v.get("exome") else 0
                genome_af = v.get("genome", {}).get("af", 0) if v.get("genome") else 0
                max_af = max(exome_af or 0, genome_af or 0)
                
                if max_af >= args.maf:
                    chrom_norm = f"chr{v_chrom_norm}"
                    locus_key = (chrom_norm, pos)
                    if locus_key not in seen_loci:
                        seen_loci.add(locus_key)
                        
                        # Add to final list matching original chromosome formatting
                        out_chrom = f"chr{v_chrom_norm}" if not is_grch37 else v_chrom_norm
                        new_row = [""] * len(header)
                        new_row[c_idx] = out_chrom
                        new_row[p_idx] = str(pos)
                        new_row[r_idx] = ref
                        new_row[a_idx] = alt
                        if s_idx != -1:
                            new_row[s_idx] = "false"
                        final_rows.append(new_row)
                        new_region_rows.append(new_row)
                        extra_snps_count += 1
                        region_added_count += 1

            save_checkpoint(q_key, new_region_rows)
            if region_added_count > 0:
                print(f"  Added {region_added_count} new common SNPs from region {chrom}:{start}-{end}")
    else:
        genes_list = [g.strip() for g in args.genes.split(",") if g.strip()]
        print(f"Querying gnomAD by gene for {len(genes_list)} genes...")
        for idx, gene in enumerate(genes_list):
            print(f"[{idx+1}/{len(genes_list)}] Querying variants for gene: {gene}...")
            variants = query_gnomad_variants_gene(gene, genome_version=ref_genome)
            time.sleep(10.0) # rate limit politeness
            print(f"  Retrieved {len(variants)} variants from gnomAD for {gene}.")
            
            gene_added_count = 0
            for v in variants:
                vid = v.get("variant_id")
                if not vid: continue
                v_chrom, pos_str, ref, alt = vid.split("-")
                pos = int(pos_str)
                
                # Check if variant falls inside targeted BED regions
                if v_chrom not in bed_regions: continue
                if not any(b_start <= pos <= b_end for b_start, b_end in bed_regions[v_chrom]):
                    continue
                    
                # Must be SNV
                if len(ref) != 1 or len(alt) != 1: continue
                
                # Filter by MAF
                exome_af = v.get("exome", {}).get("af", 0) if v.get("exome") else 0
                genome_af = v.get("genome", {}).get("af", 0) if v.get("genome") else 0
                max_af = max(exome_af or 0, genome_af or 0)
                
                if max_af >= args.maf:
                    chrom_norm = f"chr{v_chrom.replace('chr', '')}"
                    locus_key = (chrom_norm, pos)
                    if locus_key not in seen_loci:
                        seen_loci.add(locus_key)
                        
                        out_chrom = f"chr{v_chrom.replace('chr', '')}" if not is_grch37 else v_chrom.replace("chr", "")
                        new_row = [""] * len(header)
                        new_row[c_idx] = out_chrom
                        new_row[p_idx] = str(pos)
                        new_row[r_idx] = ref
                        new_row[a_idx] = alt
                        if s_idx != -1:
                            new_row[s_idx] = "false"
                        final_rows.append(new_row)
                        extra_snps_count += 1
                        gene_added_count += 1
            print(f"  Added {gene_added_count} new common SNPs in targeted BED regions of {gene}")
                        
    print(f"Added total of {extra_snps_count} new common panel-targeted SNPs to the reference set.")
    
    # Sort final rows by Chromosome and Position
    def sort_key(row):
        chr_str = row[c_idx].replace("chr", "")
        if chr_str == "X": chr_num = 100
        elif chr_str == "Y": chr_num = 101
        elif chr_str == "M" or chr_str == "MT": chr_num = 102
        else:
            try: chr_num = int(chr_str)
            except ValueError: chr_num = 999
        return (chr_num, int(row[p_idx]))
        
    final_rows.sort(key=sort_key)
    
    # 4. Write expanded file to output
    print(f"Writing expanded loci file to {args.output}...")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with gzip.open(args.output, "wt") as f:
        f.write("\t".join(header) + "\n")
        for row in final_rows:
            f.write("\t".join(row) + "\n")
            
    print(f"Done! Created expanded loci file with {len(final_rows)} sites (added {extra_snps_count} targeted region SNPs).")

if __name__ == "__main__":
    main()
