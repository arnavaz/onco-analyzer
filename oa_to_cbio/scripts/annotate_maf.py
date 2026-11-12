import json
import sys
import os

if len(sys.argv) < 5:
    print("Usage: python3 annotate_maf.py <input_maf> <orange_json> <output_maf> <mutation_status> [--unfiltered]")
    sys.exit(1)

in_maf = sys.argv[1]
orange_json = sys.argv[2]
out_maf = sys.argv[3]
mutation_status = sys.argv[4]

unfiltered = "--unfiltered" in sys.argv
non_germline_filter = "--non-germline-filter" in sys.argv

# cBioPortal default exclusion list
EXCLUDED_CLASSES = {
    'Silent', 'Intron', 'IGR', "3'UTR", "5'UTR", "3'Flank", "5'Flank", 'RNA'
}

# Coordinate matching helper function
def match_variant(maf_chrom, maf_start, maf_ref, maf_alt, json_chrom, json_pos, json_ref, json_alt):
    # Normalize chromosomes
    m_chrom = str(maf_chrom).upper().replace('CHR', '')
    j_chrom = str(json_chrom).upper().replace('CHR', '')
    if m_chrom != j_chrom:
        return False
        
    try:
        m_start = int(maf_start)
        j_pos = int(json_pos)
    except ValueError:
        return False
    
    # 1. Exact match
    if m_start == j_pos and maf_ref == json_ref and maf_alt == json_alt:
        return True
        
    # 2. Insertion match (VCF to MAF coords)
    if len(json_ref) < len(json_alt):
        converted_start = j_pos + 1
        converted_ref = '-'
        converted_alt = json_alt[1:]
        if m_start == converted_start and maf_ref == converted_ref and maf_alt == converted_alt:
            return True
            
    # 3. Deletion match (VCF to MAF coords)
    if len(json_ref) > len(json_alt):
        converted_start = j_pos + 1
        converted_ref = json_ref[1:]
        converted_alt = '-'
        if m_start == converted_start and maf_ref == converted_ref and maf_alt == converted_alt:
            return True
            
    return False

# Function to search reported variants
def find_matching_reported_variant(reported_vars, hugo, chrom, start, ref, alt, transcript, hgvsc, hgvsp):
    candidates = reported_vars.get(hugo, [])
    if not candidates:
        return None
        
    for v in candidates:
        # Check coordinate match
        if match_variant(chrom, start, ref, alt, v['chromosome'], v['position'], v['ref'], v['alt']):
            return v
            
        # Fallback to transcript and HGVS match
        v_tx = v.get('transcript')
        if v_tx and transcript and v_tx == transcript:
            v_hgvsc = v.get('hgvsc')
            v_hgvsp = v.get('hgvsp')
            if v_hgvsc and hgvsc and v_hgvsc == hgvsc:
                return v
            if v_hgvsp and hgvsp and v_hgvsp == hgvsp:
                return v
                
    return None

# Parse JSON for driver variants
somatic_reported_vars = {}
germline_reported_vars = {}

if os.path.exists(orange_json):
    try:
        with open(orange_json, 'r') as f:
            data = json.load(f)
            
            # Somatic drivers
            somatic_drivers = {}
            for d in data.get('purple', {}).get('somaticDrivers', []):
                if d.get('type') == 'MUTATION':
                    gene = d.get('gene')
                    tx = d.get('transcript')
                    dl = d.get('driverLikelihood')
                    if gene and dl is not None:
                        somatic_drivers[(gene, tx)] = float(dl)
            
            # Somatic reported/suspect variants
            for key in ['reportableSomaticVariants', 'additionalSuspectSomaticVariants']:
                for v in data.get('purple', {}).get(key, []):
                    gene = v.get('gene')
                    if not gene:
                        continue
                    impact = v.get('canonicalImpact', {})
                    tx = impact.get('transcript')
                    hgvsc = impact.get('hgvsCodingImpact')
                    hgvsp = impact.get('hgvsProteinImpact')
                    
                    dl = somatic_drivers.get((gene, tx))
                    if dl is None:
                        # Gene-level fallback
                        for (dg, dt), dval in somatic_drivers.items():
                            if dg == gene:
                                dl = dval
                                break
                    
                    var_info = {
                        'chromosome': v.get('chromosome'),
                        'position': v.get('position'),
                        'ref': v.get('ref'),
                        'alt': v.get('alt'),
                        'transcript': tx,
                        'hgvsc': hgvsc,
                        'hgvsp': hgvsp,
                        'driverLikelihood': dl,
                        'status': 'Reportable' if key == 'reportableSomaticVariants' else 'Suspect'
                    }
                    if gene not in somatic_reported_vars:
                        somatic_reported_vars[gene] = []
                    somatic_reported_vars[gene].append(var_info)
                    
            # Germline drivers
            germline_drivers = {}
            for d in data.get('purple', {}).get('germlineDrivers', []):
                if d.get('type') == 'MUTATION':
                    gene = d.get('gene')
                    tx = d.get('transcript')
                    dl = d.get('driverLikelihood')
                    if gene and dl is not None:
                        germline_drivers[(gene, tx)] = float(dl)
            
            # Germline reported/suspect variants
            for key in ['reportableGermlineVariants', 'additionalSuspectGermlineVariants']:
                for v in data.get('purple', {}).get(key, []):
                    gene = v.get('gene')
                    if not gene:
                        continue
                    impact = v.get('canonicalImpact', {})
                    tx = impact.get('transcript')
                    hgvsc = impact.get('hgvsCodingImpact')
                    hgvsp = impact.get('hgvsProteinImpact')
                    
                    dl = germline_drivers.get((gene, tx))
                    if dl is None:
                        # Gene-level fallback
                        for (dg, dt), dval in germline_drivers.items():
                            if dg == gene:
                                dl = dval
                                break
                                
                    status_text = "Reportable Germline Variant" if key == 'reportableGermlineVariants' else "Suspect Germline Variant"
                    var_info = {
                        'chromosome': v.get('chromosome'),
                        'position': v.get('position'),
                        'ref': v.get('ref'),
                        'alt': v.get('alt'),
                        'transcript': tx,
                        'hgvsc': hgvsc,
                        'hgvsp': hgvsp,
                        'driverLikelihood': dl,
                        'status': status_text
                    }
                    if gene not in germline_reported_vars:
                        germline_reported_vars[gene] = []
                    germline_reported_vars[gene].append(var_info)
                    
    except Exception as e:
        print(f"Warning: Failed to parse {orange_json}: {e}", file=sys.stderr)

# 3. Process MAF
if not os.path.exists(in_maf):
    print(f"Warning: Input MAF not found: {in_maf}", file=sys.stderr)
    sys.exit(0)

with open(in_maf, 'r') as fin, open(out_maf, 'w') as fout:
    header = None
    cbp_cols = ['cbp_driver', 'cbp_driver_annotation', 'cbp_driver_tiers', 'cbp_driver_tiers_annotation']
    
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        
        parts = line.rstrip('\n').split('\t')
        
        if header is None:
            header = parts
            # Add CBP columns if not present
            for col in cbp_cols:
                if col not in header:
                    header.append(col)
            if 'Mutation_Status' not in header:
                header.append('Mutation_Status')
            fout.write('\t'.join(header) + '\n')
            
            # Map indices
            hugo_idx = header.index('Hugo_Symbol') if 'Hugo_Symbol' in header else -1
            status_idx = header.index('Mutation_Status')
            class_idx = header.index('Variant_Classification') if 'Variant_Classification' in header else -1
            
            chrom_idx = header.index('Chromosome') if 'Chromosome' in header else -1
            start_idx = header.index('Start_Position') if 'Start_Position' in header else -1
            ref_idx = header.index('Reference_Allele') if 'Reference_Allele' in header else -1
            alt_idx = header.index('Tumor_Seq_Allele2') if 'Tumor_Seq_Allele2' in header else -1
            
            transcript_idx = header.index('Transcript_ID') if 'Transcript_ID' in header else -1
            hgvsc_idx = header.index('HGVSc') if 'HGVSc' in header else -1
            hgvsp_idx = header.index('HGVSp_Short') if 'HGVSp_Short' in header else (header.index('HGVSp') if 'HGVSp' in header else -1)
            continue
            
        # Pad row to match header length
        while len(parts) < len(header):
            parts.append('')
            
        parts[status_idx] = mutation_status
        hugo = parts[hugo_idx] if hugo_idx >= 0 else ''
        var_class = parts[class_idx] if class_idx >= 0 else ''
        
        chrom = parts[chrom_idx] if chrom_idx >= 0 else ''
        start = parts[start_idx] if start_idx >= 0 else ''
        ref = parts[ref_idx] if ref_idx >= 0 else ''
        alt = parts[alt_idx] if alt_idx >= 0 else ''
        
        transcript = parts[transcript_idx] if transcript_idx >= 0 else ''
        hgvsc = parts[hgvsc_idx] if hgvsc_idx >= 0 else ''
        hgvsp = parts[hgvsp_idx] if hgvsp_idx >= 0 else ''
        
        # Match variant
        matched_var = None
        if mutation_status == 'Somatic':
            matched_var = find_matching_reported_variant(
                somatic_reported_vars, hugo, chrom, start, ref, alt, transcript, hgvsc, hgvsp
            )
        elif mutation_status == 'Germline':
            matched_var = find_matching_reported_variant(
                germline_reported_vars, hugo, chrom, start, ref, alt, transcript, hgvsc, hgvsp
            )
            # If this is a germline MAF, strict filter to only include genes/variants annotated in orange.json
            if matched_var is None and not non_germline_filter:
                continue
        
        # Apply filter if not unfiltered
        if not unfiltered and var_class in EXCLUDED_CLASSES:
            # Keep TERT promoter variants (5'Flank)
            is_tert = (hugo == 'TERT' and var_class == "5'Flank")
            # Keep any variant that has driver annotations
            is_driver = (matched_var is not None)
            
            if not (is_tert or is_driver):
                continue
        
        if matched_var is not None:
            dl = matched_var.get('driverLikelihood')
            if dl is not None:
                # Map Driver Likelihood to cBioPortal fields
                if dl >= 0.8:
                    cbp_driver = 'Putative_Driver'
                    cbp_driver_annotation = 'High Likelihood Driver'
                elif dl > 0.1:
                    cbp_driver = 'Unknown'
                    cbp_driver_annotation = 'VUS / Moderate Likelihood'
                else:
                    cbp_driver = 'Putative_Passenger'
                    cbp_driver_annotation = 'Low Likelihood'
                    
                cbp_driver_tiers = f"DL {dl*100:.0f}%"
                cbp_driver_tiers_annotation = "HMFTools Driver Likelihood"
                
                parts[header.index('cbp_driver')] = cbp_driver
                parts[header.index('cbp_driver_annotation')] = cbp_driver_annotation
                parts[header.index('cbp_driver_tiers')] = cbp_driver_tiers
                parts[header.index('cbp_driver_tiers_annotation')] = cbp_driver_tiers_annotation
            else:
                # Fallback for variants that are reportable/suspect but lack a driver likelihood score
                status_text = matched_var.get('status', 'Reportable Variant')
                parts[header.index('cbp_driver')] = 'Unknown' if 'Suspect' in status_text else 'Putative_Driver'
                parts[header.index('cbp_driver_annotation')] = status_text
                parts[header.index('cbp_driver_tiers')] = 'N/A'
                parts[header.index('cbp_driver_tiers_annotation')] = 'No Driver Likelihood Provided'
                
        fout.write('\t'.join(parts) + '\n')
