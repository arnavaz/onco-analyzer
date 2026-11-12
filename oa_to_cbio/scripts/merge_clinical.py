import argparse
import os
import glob
import re
from collections import OrderedDict
import sys
from icd10_translator import get_oncotree_name, translate_icd10

def parse_cbio_clinical(filepath, id_column):
    headers = []
    metadata_lines = []
    data = OrderedDict()
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('#') or line.startswith('"#'):
                metadata_lines.append([p.strip('"') for p in line.split('\t')])
                continue
            
            parts = [p.strip('"') for p in line.split('\t')]
            if not headers:
                headers = parts
                id_idx = headers.index(id_column) if id_column in headers else 0
                continue
            
            # Data row
            if len(parts) > 0:
                row_id = parts[id_idx] if id_idx < len(parts) else ""
                row_id = row_id.strip()
                
                # Skip rows where the ID is blank, NA, or the first column is blank/NA
                if not row_id or row_id.upper() == "NA":
                    continue
                if not parts[0].strip() or parts[0].strip().upper() == "NA":
                    continue
                    
                row_data = {headers[i]: parts[i] if i < len(parts) else "NA" for i in range(len(headers))}
                data[row_id] = row_data
                
    return headers, metadata_lines, data

def write_cbio_clinical(filepath, headers, metadata_lines, data):
    ACRONYM_TOKENS = {"DNA", "HR", "NTHL1", "QC"}
    HLA_DISPLAY_NAMES = {
        "HLA_DQA1_1": "HLA DQA1 1",
        "HLA_DPB1_2": "HLA DPB1 2",
        "HLA_DQB1_2": "HLA DQB1 2",
        "HLA_DPA1_1": "HLA DPA1 1",
        "HLA_A_1": "HLA A 1",
        "HLA_B_1": "HLA B 1",
        "HLA_C_1": "HLA C 1",
    }
    ACRONYM_DISPLAY_NAMES = {
        "TMBPERMB": "TMB per Mb",
        "TMB_PER_MB": "TMB per Mb",
        "MISALLOC": "MIS Alloc",
        "MIS_ALLOC": "MIS Alloc",
        "MSISCORE": "MSI Score",
        "MSI_SCORE": "MSI Score",
        "MSISTATUS": "MSI Status",
        "MSI_STATUS": "MSI Status",
        "CUPPA_PREDICTED_CANCER_TYPE": "CUPPA Predicted Cancer Type",
        "CUPPA_LIKELIHOOD": "CUPPA Likelihood",
        "CUPPA_TOP_PREDICTIONS": "CUPPA Top Predictions",
    }
    VIRUS_SUFFIX_DISPLAY_NAMES = {
        "PERCENTAGECOVERED": "Percentage Covered",
        "MEANCOVERAGE": "Mean Coverage",
        "EXPECTEDCLONALCOVERAGE": "Expected Clonal Coverage",
        "DRIVERINTERPRETATION": "Driver Interpretation",
    }

    def format_header_display(header):
        if header.upper() in HLA_DISPLAY_NAMES:
            return HLA_DISPLAY_NAMES[header.upper()]

        header_upper = header.upper()
        if header_upper.startswith("HLA_"):
            parts = [p for p in header.split("_") if p]
            if len(parts) >= 3:
                locus = parts[1].upper()
                suffix = parts[-1]
                if suffix.isdigit():
                    return f"HLA {locus} {suffix}"
                return f"HLA {' '.join(parts[1:]).upper()}"

        virus_match = re.fullmatch(r"VIRUS_(\d+)_(.+)", header_upper)
        if virus_match:
            virus_num = virus_match.group(1)
            suffix = virus_match.group(2)
            suffix_display = VIRUS_SUFFIX_DISPLAY_NAMES.get(suffix)
            if suffix_display:
                return f"Virus {virus_num} {suffix_display}"

        lowercase_words = {"a", "an", "and", "as", "at", "but", "by", "for", "in", "nor", "of", "on", "or", "the", "to", "via"}
        display_tokens = []
        for token in header.split("_"):
            token_upper = token.upper()
            if token_upper in ACRONYM_TOKENS or re.fullmatch(r"(SIG|SBS|DBS|ID)\d+", token_upper):
                display_tokens.append(token_upper)
            elif token_upper == "QCSTATUS":
                display_tokens.extend(["QC", "Status"])
            else:
                display_tokens.append(token.lower() if token.lower() in lowercase_words else token.title())
        if display_tokens and display_tokens[0].lower() in lowercase_words:
            display_tokens[0] = display_tokens[0].title()
        return " ".join(display_tokens)

    CLINICAL_DATATYPE_OVERRIDES = {
        "OS_MONTHS": "NUMBER",
        "AGE": "NUMBER",
        "DFS_MONTHS": "NUMBER",
        "DFS_STATUS": "STRING",
        "OS_STATUS": "STRING",
        "PURITY": "NUMBER",
        "PLOIDY": "NUMBER",
        "MSISCORE": "NUMBER",
        "TMBPERMB": "NUMBER",
        "MEAN_COVERAGE": "NUMBER",
        "MEDIAN_COVERAGE": "NUMBER",
        "MEAN_TELOMERE_LENGTH_TUMOUR_BP": "NUMBER",
        "MEAN_TELOMERE_LENGTH_TUMOUR_RAW_BP": "NUMBER",
        "MEAN_TELOMERE_LENGTH_NORMAL_BP": "NUMBER",
        "MEAN_TELOMERE_LENGTH_DELTA_BP": "NUMBER",
        "MEAN_TELOMERE_LENGTH_TUMOR_BP": "NUMBER",
        "MEAN_TELOMERE_LENGTH_RAW_BP": "NUMBER",
        "DUPLICATE_READ_PROPORTION_TUMOUR": "NUMBER",
        "DUPLICATE_READ_PROPORTION_NORMAL": "NUMBER",
        "ECDNA_STATUS": "STRING",
        "ECDNA_CLUSTER_COUNT": "NUMBER",
        "ECDNA_AMPLICON_GENES": "STRING",
        "ECDNA_MAX_COPY_NUMBER": "NUMBER",
        "ECDNA_CLOSED_SEG_LENGTH_MB": "NUMBER",
        "CUPPA_PREDICTED_CANCER_TYPE": "STRING",
        "CUPPA_LIKELIHOOD": "NUMBER",
        "CUPPA_TOP_PREDICTIONS": "STRING",
    }

    for idx, header in enumerate(headers):
        if idx == 0:
            continue
        acronym_display = ACRONYM_DISPLAY_NAMES.get(header.upper()) or HLA_DISPLAY_NAMES.get(header.upper())
        if acronym_display:
            metadata_lines[0][idx] = acronym_display
            if metadata_lines[1][idx] == header or "_" in metadata_lines[1][idx]:
                metadata_lines[1][idx] = acronym_display
        elif metadata_lines[0][idx] == header or "_" in metadata_lines[0][idx]:
            metadata_lines[0][idx] = format_header_display(header)
            metadata_lines[1][idx] = metadata_lines[0][idx]
    
    # Correct any incorrect datatypes in the metadata datatype line (index 2)
    if len(metadata_lines) >= 3:
        for idx, h in enumerate(headers):
            h_upper = h.upper()
            if h_upper in CLINICAL_DATATYPE_OVERRIDES:
                expected_type = CLINICAL_DATATYPE_OVERRIDES[h_upper]
                if idx == 0:
                    metadata_lines[2][idx] = f"#{expected_type}"
                else:
                    metadata_lines[2][idx] = expected_type
            elif re.match(r"^(SIG\d+|SBS\d+|DBS\d+|ID\d+|MISALLOC|MIS_ALLOC)(_|$)", h_upper):
                if idx == 0:
                    metadata_lines[2][idx] = "#NUMBER"
                else:
                    metadata_lines[2][idx] = "NUMBER"

    with open(filepath, 'w') as f:
        for meta in metadata_lines:
            f.write("\t".join(meta) + "\n")
        f.write("\t".join(headers) + "\n")
        for row_id, row_dict in data.items():
            f.write("\t".join([row_dict.get(h, "NA") for h in headers]) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp-dir', required=True)
    parser.add_argument('--out-sample', required=True)
    parser.add_argument('--out-patient', required=True)
    parser.add_argument('--pre-sample', default=None)
    parser.add_argument('--pre-patient', default=None)
    args = parser.parse_args()

    # 1. Read existing sample data if provided
    sample_headers = ["PATIENT_ID", "SAMPLE_ID"]
    sample_meta = [
        ["#Patient Identifier", "Sample Identifier"],
        ["#Patient Identifier", "Sample Identifier"],
        ["#STRING", "STRING"],
        ["#1", "1"]
    ]
    sample_data = OrderedDict()
    
    if args.pre_sample and os.path.exists(args.pre_sample):
        p_headers, p_meta, p_data = parse_cbio_clinical(args.pre_sample, "SAMPLE_ID")
        if p_headers: sample_headers = p_headers
        if len(p_meta) >= 4: sample_meta = p_meta
        elif len(p_meta) > 0:
            while len(p_meta) < 4: p_meta.append(["#"] * len(sample_headers))
            sample_meta = p_meta
        sample_data = p_data

    # 2. Read new sample data from TSVs
    tsv_files = glob.glob(os.path.join(args.tmp_dir, "*.clinical.tsv"))
    new_keys = []
    
    for tsv in tsv_files:
        s_id = None
        s_data = {}
        with open(tsv, 'r') as f:
            for line in f:
                parts = line.rstrip('\r\n').split('\t', 1)
                if len(parts) == 2:
                    s_data[parts[0]] = parts[1]
        if s_data.get("IS_NORMAL") == "Yes" or len(s_data) <= 3:
            continue
            
        for k in s_data:
            if k not in new_keys and k not in ["PATIENT_ID", "SAMPLE_ID", "IS_NORMAL"]:
                new_keys.append(k)
                
        s_id = s_data.get("SAMPLE_ID")
        if s_id:
            # These cohorts encode the donor ID in the first three sample ID components.
            cohort_patient_id = None
            if s_id.startswith(("PANCAN-", "DYNA-", "MOHCCN-")):
                parts = s_id.split('-')
                if len(parts) >= 3:
                    cohort_patient_id = "-".join(parts[:3])
                    s_data["PATIENT_ID"] = cohort_patient_id
            
            # Try to find a matching existing sample ID (e.g. handle -R suffixes)
            target_id = s_id
            if target_id not in sample_data and f"{target_id}-R" in sample_data:
                target_id = f"{target_id}-R"
                
            if target_id in sample_data:
                # Correct stale cohort IDs; otherwise preserve the existing source ID.
                if cohort_patient_id:
                    sample_data[target_id]["PATIENT_ID"] = cohort_patient_id
                elif "PATIENT_ID" in sample_data[target_id]:
                    s_data.pop("PATIENT_ID", None)
            else:
                sample_data[target_id] = {}
                
            # Update the existing dictionary with new/overwritten values
            sample_data[target_id].update(s_data)
            # Ensure SAMPLE_ID matches the key
            sample_data[target_id]["SAMPLE_ID"] = target_id

    # Ensure all existing metadata lines are padded to the length of sample_headers
    for idx in range(len(sample_meta)):
        while len(sample_meta[idx]) < len(sample_headers):
            if idx == 0 or idx == 1:
                col_name = sample_headers[len(sample_meta[idx])]
                sample_meta[idx].append(col_name)
            elif idx == 2:
                sample_meta[idx].append("STRING")
            elif idx == 3:
                sample_meta[idx].append("1")
            else:
                sample_meta[idx].append("NA")

    TELOMERE_METADATA = {
        "MEAN_TELOMERE_LENGTH_TUMOUR_BP": {
            "display": "Mean Telomere Length (Tumour Corrected, bp)",
            "desc": "Purity- and ploidy-corrected mean telomere length per chromosome end in base pairs",
            "type": "NUMBER",
            "priority": "1"
        },
        "MEAN_TELOMERE_LENGTH_TUMOUR_RAW_BP": {
            "display": "Mean Telomere Length (Tumour Uncorrected, bp)",
            "desc": "Uncorrected bulk mean telomere length per chromosome end in base pairs",
            "type": "NUMBER",
            "priority": "2"
        },
        "MEAN_TELOMERE_LENGTH_NORMAL_BP": {
            "display": "Mean Telomere Length (Matched Normal, bp)",
            "desc": "Mean telomere length per chromosome end in base pairs for matched normal sample",
            "type": "NUMBER",
            "priority": "2"
        },
        "MEAN_TELOMERE_LENGTH_DELTA_BP": {
            "display": "Mean Telomere Length Delta (Tumour - Normal, bp)",
            "desc": "Difference in mean telomere length between tumour and matched normal sample",
            "type": "NUMBER",
            "priority": "2"
        },
        "DUPLICATE_READ_PROPORTION_TUMOUR": {
            "display": "Duplicate Read Fraction (Tumour)",
            "desc": "Fraction of marked duplicate reads in tumour sequencing BAM",
            "type": "NUMBER",
            "priority": "3"
        },
        "DUPLICATE_READ_PROPORTION_NORMAL": {
            "display": "Duplicate Read Fraction (Matched Normal)",
            "desc": "Fraction of marked duplicate reads in matched normal sequencing BAM",
            "type": "NUMBER",
            "priority": "4"
        },
        # Legacy/fallback definitions:
        "MEAN_TELOMERE_LENGTH_TUMOR_BP": {
            "display": "Mean Telomere Length (Tumour Corrected, bp)",
            "desc": "Purity- and ploidy-corrected mean telomere length per chromosome end in base pairs",
            "type": "NUMBER",
            "priority": "1"
        },
        "MEAN_TELOMERE_LENGTH_RAW_BP": {
            "display": "Mean Telomere Length (Tumour Uncorrected, bp)",
            "desc": "Uncorrected bulk mean telomere length per chromosome end in base pairs",
            "type": "NUMBER",
            "priority": "2"
        },
        "TELOMERE_LENGTH_DELTA_BP": {
            "display": "Telomere Length Delta (bp)",
            "desc": "Difference in corrected telomere length between tumour and matched normal sample",
            "type": "NUMBER",
            "priority": "2"
        },
        "ECDNA_STATUS": {
            "display": "ecDNA Status",
            "desc": "Indicates if circular extrachromosomal DNA (Double Minute) was detected in the sample (Yes/No)",
            "type": "STRING",
            "priority": "1"
        },
        "ECDNA_CLUSTER_COUNT": {
            "display": "ecDNA Cluster Count",
            "desc": "Total number of circular ecDNA/Double Minute clusters detected in the sample by LINX",
            "type": "NUMBER",
            "priority": "2"
        },
        "ECDNA_AMPLICON_GENES": {
            "display": "ecDNA Amplified Genes",
            "desc": "Semicolon-separated list of amplified genes residing on the detected ecDNA segments",
            "type": "STRING",
            "priority": "2"
        },
        "ECDNA_MAX_COPY_NUMBER": {
            "display": "ecDNA Max Copy Number",
            "desc": "Maximum copy number observed among the circular ecDNA amplicons in the sample",
            "type": "NUMBER",
            "priority": "2"
        },
        "ECDNA_CLOSED_SEG_LENGTH_MB": {
            "display": "ecDNA Closed Segment Length (Mb)",
            "desc": "Total genomic length of circularized closed segments in Megabases",
            "type": "NUMBER",
            "priority": "2"
        }
    }

    CUPPA_METADATA = {
        "CUPPA_PREDICTED_CANCER_TYPE": {
            "display": "CUPPA Predicted Cancer Type",
            "desc": "Cancer of Unknown Primary (CUPPA) predicted primary cancer type",
            "type": "STRING",
            "priority": "1"
        },
        "CUPPA_LIKELIHOOD": {
            "display": "CUPPA Likelihood",
            "desc": "CUPPA classifier prediction confidence likelihood score (0 to 1)",
            "type": "NUMBER",
            "priority": "2"
        },
        "CUPPA_TOP_PREDICTIONS": {
            "display": "CUPPA Top Predictions",
            "desc": "Top candidate cancer types and associated likelihood scores predicted by CUPPA",
            "type": "STRING",
            "priority": "2"
        }
    }

    # Add any new columns to headers and metadata
    for k in new_keys:
        if k not in sample_headers:
            sample_headers.append(k)
            
            # Check if this is a predefined metadata column
            if k in TELOMERE_METADATA:
                meta_item = TELOMERE_METADATA[k]
                display = meta_item["display"]
                desc = meta_item["desc"]
                datatype = meta_item["type"]
                priority = meta_item["priority"]
            elif k in CUPPA_METADATA:
                meta_item = CUPPA_METADATA[k]
                display = meta_item["display"]
                desc = meta_item["desc"]
                datatype = meta_item["type"]
                priority = meta_item["priority"]
            else:
                display = k
                desc = k
                # Determine datatype based on values across all samples
                is_numeric = True
                has_value = False
                for s_id, s_dict in sample_data.items():
                    val = s_dict.get(k)
                    if val and val != "NA":
                        has_value = True
                        try:
                            float(val)
                        except ValueError:
                            is_numeric = False
                            break
                datatype = "NUMBER" if (is_numeric and has_value) else "STRING"
                priority = "1"
            
            # Pad all metadata lines in sample_meta
            for idx in range(len(sample_meta)):
                if idx == 0:
                    sample_meta[idx].append(display)
                elif idx == 1:
                    sample_meta[idx].append(desc)
                elif idx == 2:
                    sample_meta[idx].append(datatype)
                elif idx == 3:
                    sample_meta[idx].append(priority)
                else:
                    sample_meta[idx].append("NA")

    # Apply ICD-10 to Cancer Type translation if diagnosis fields contain ICD-10 codes
    for s_id, s_dict in sample_data.items():
        for field in ["CANCER_TYPE", "PRIMARY_DIAGNOSIS", "DIAGNOSIS", "DIAGNOSIS_CODE"]:
            val = s_dict.get(field)
            if val and str(val).strip().upper() != "NA":
                broad, detailed, oncotree = translate_icd10(val)
                if broad != val:
                    s_dict["CANCER_TYPE"] = broad
                    curr_det = s_dict.get("CANCER_TYPE_DETAILED", "")
                    if not curr_det or curr_det.upper() in ["NA", ""] or curr_det == val or len(curr_det) <= 15:
                        s_dict["CANCER_TYPE_DETAILED"] = detailed
                    if oncotree != "NA" and ("ONCOTREE_CODE" not in s_dict or s_dict.get("ONCOTREE_CODE") in ["NA", ""]):
                        s_dict["ONCOTREE_CODE"] = oncotree
                    break

    # Prefer the canonical OncoTree name over a free-text diagnosis description.
    for s_dict in sample_data.values():
        oncotree_name = get_oncotree_name(s_dict.get("ONCOTREE_CODE"))
        if oncotree_name:
            s_dict["CANCER_TYPE_DETAILED"] = oncotree_name

    # Ensure CANCER_TYPE, CANCER_TYPE_DETAILED, and ONCOTREE_CODE are in sample_headers if present in data
    for col in ["CANCER_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE"]:
        if col not in sample_headers and any(col in s for s in sample_data.values()):
            sample_headers.append(col)
            for idx in range(len(sample_meta)):
                if idx == 0: sample_meta[idx].append("Cancer Type" if col == "CANCER_TYPE" else ("Cancer Type Detailed" if col == "CANCER_TYPE_DETAILED" else "OncoTree Code"))
                elif idx == 1: sample_meta[idx].append("Cancer Type" if col == "CANCER_TYPE" else ("Detailed Cancer Type" if col == "CANCER_TYPE_DETAILED" else "OncoTree Code"))
                elif idx == 2: sample_meta[idx].append("STRING")
                elif idx == 3: sample_meta[idx].append("1")
                else: sample_meta[idx].append("NA")

    # Ensure ECDNA_AMPLICON_GENES is "None" if ECDNA_STATUS is "No"
    for s_dict in sample_data.values():
        if s_dict.get("ECDNA_STATUS") == "No":
            if s_dict.get("ECDNA_AMPLICON_GENES") in (None, "", "NA"):
                s_dict["ECDNA_AMPLICON_GENES"] = "None"

    # Flatten unallocated/sub-threshold signatures to "0.00" for samples where signatures were evaluated
    sig_cols = [
        h for h in sample_headers
        if re.match(r"^(SIG\d+|SBS\d+|DBS\d+|ID\d+|MISALLOC|MIS_ALLOC)(_|$)", h.upper())
    ]
    if sig_cols:
        for s_dict in sample_data.values():
            has_sig_analysis = any(
                col in s_dict and s_dict[col] not in (None, "", "NA")
                for col in sig_cols
            ) or (s_dict.get("TMBPERMB") not in (None, "", "NA")) or (s_dict.get("PURITY") not in (None, "", "NA"))

            if has_sig_analysis:
                for col in sig_cols:
                    if col not in s_dict or s_dict[col] in (None, "", "NA"):
                        s_dict[col] = "0.00"

    write_cbio_clinical(args.out_sample, sample_headers, sample_meta, sample_data)

    # 3. Read existing patient data if provided
    patient_headers = ["PATIENT_ID", "SUBMITTER_DONOR_ID", "SEX"]
    patient_meta = [
        ["#Patient Identifier", "Submitter Donor ID", "Sex"],
        ["#Patient Identifier", "Submitter Donor ID", "Sex"],
        ["#STRING", "STRING", "STRING"],
        ["#1", "1", "1"]
    ]
    patient_data = OrderedDict()
    
    if args.pre_patient and os.path.exists(args.pre_patient):
        p_headers, p_meta, p_data = parse_cbio_clinical(args.pre_patient, "PATIENT_ID")
        if p_headers: patient_headers = p_headers
        if len(p_meta) >= 4: patient_meta = p_meta
        elif len(p_meta) > 0:
            while len(p_meta) < 4: p_meta.append(["#"] * len(patient_headers))
            patient_meta = p_meta
        patient_data = p_data

    PATIENT_COLUMN_DISPLAY_NAMES = {
        "PATIENT_ID": "Patient Identifier",
        "SUBMITTER_DONOR_ID": "Submitter Donor ID",
        "PROGRAM_ID": "Program ID",
        "COHORT_ID": "Cohort ID",
        "GENDER": "Gender",
        "SEX": "Sex",
        "SEX_AT_BIRTH": "Sex at Birth",
        "IS_DECEASED": "Is Deceased",
        "CAUSE_OF_DEATH": "Cause of Death",
        "LOST_TO_FOLLOWUP_REASON": "Lost to Follow-up Reason",
        "CANCER_TYPE_CODE": "Cancer Type Code",
        "BASIS_OF_DIAGNOSIS": "Basis of Diagnosis",
        "LATERALITY": "Laterality",
        "CLINICAL_TUMOUR_STAGING_SYSTEM": "Clinical Tumour Staging System",
        "CLINICAL_T_CATEGORY": "Clinical T Category",
        "CLINICAL_N_CATEGORY": "Clinical N Category",
        "CLINICAL_M_CATEGORY": "Clinical M Category",
        "CLINICAL_STAGE_GROUP": "Clinical Stage Group",
        "PATHOLOGICAL_TUMOUR_STAGING_SYSTEM": "Pathological Tumour Staging System",
        "PATHOLOGICAL_T_CATEGORY": "Pathological T Category",
        "PATHOLOGICAL_N_CATEGORY": "Pathological N Category",
        "PATHOLOGICAL_M_CATEGORY": "Pathological M Category",
        "PATHOLOGICAL_STAGE_GROUP": "Pathological Stage Group",
        "TOBACCO_SMOKING_STATUS": "Tobacco Smoking Status",
        "TOBACCO_TYPE": "Tobacco Type",
        "PACK_YEARS_SMOKED": "Pack Years Smoked",
        "OS_MONTHS": "Overall Survival (Months)",
        "OS_STATUS": "Overall Survival Status",
        "DFS_MONTHS": "Disease Free (Months)",
        "DFS_STATUS": "Disease Free Status",
        "AGE": "Age at Diagnosis",
    }

    # Ensure all existing metadata lines are padded to the length of patient_headers
    for idx in range(len(patient_meta)):
        while len(patient_meta[idx]) < len(patient_headers):
            col_name = patient_headers[len(patient_meta[idx])]
            if idx == 0 or idx == 1:
                display = PATIENT_COLUMN_DISPLAY_NAMES.get(col_name) or col_name.replace("_", " ").title()
                patient_meta[idx].append(display)
            elif idx == 2:
                datatype = "NUMBER" if col_name in ["OS_MONTHS", "DFS_MONTHS", "AGE", "PACK_YEARS_SMOKED"] or col_name.startswith("AGE_") else "STRING"
                patient_meta[idx].append(datatype)
            elif idx == 3:
                patient_meta[idx].append("1")
            else:
                patient_meta[idx].append("NA")

    # 4. Gather all PATIENT_IDs from the final sample data
    for s_id, s_dict in sample_data.items():
        p_id = s_dict.get("PATIENT_ID")
        if p_id:
            p_id = p_id.strip()
            if p_id and p_id.upper() != "NA" and p_id not in patient_data:
                patient_data[p_id] = {"PATIENT_ID": p_id}

    # Apply ICD-10 to Cancer Type translation if patient diagnosis fields contain ICD-10 codes
    for p_id, p_dict in patient_data.items():
        for field in ["CANCER_TYPE", "PRIMARY_DIAGNOSIS", "DIAGNOSIS", "DIAGNOSIS_CODE"]:
            val = p_dict.get(field)
            if val and str(val).strip().upper() != "NA":
                broad, detailed, oncotree = translate_icd10(val)
                if broad != val:
                    p_dict["CANCER_TYPE"] = broad
                    curr_det = p_dict.get("CANCER_TYPE_DETAILED", "")
                    if not curr_det or curr_det.upper() in ["NA", ""] or curr_det == val or len(curr_det) <= 15:
                        p_dict["CANCER_TYPE_DETAILED"] = detailed
                    if oncotree != "NA" and ("ONCOTREE_CODE" not in p_dict or p_dict.get("ONCOTREE_CODE") in ["NA", ""]):
                        p_dict["ONCOTREE_CODE"] = oncotree
                    break

    # Prefer the canonical OncoTree name over a free-text diagnosis description.
    for p_dict in patient_data.values():
        oncotree_name = get_oncotree_name(p_dict.get("ONCOTREE_CODE"))
        if oncotree_name:
            p_dict["CANCER_TYPE_DETAILED"] = oncotree_name

    # Ensure CANCER_TYPE, CANCER_TYPE_DETAILED, and ONCOTREE_CODE are in patient_headers if present in patient_data
    for col in ["CANCER_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE"]:
        if col not in patient_headers and any(col in p for p in patient_data.values()):
            patient_headers.append(col)
            for idx in range(len(patient_meta)):
                if idx == 0: patient_meta[idx].append("Cancer Type" if col == "CANCER_TYPE" else ("Cancer Type Detailed" if col == "CANCER_TYPE_DETAILED" else "OncoTree Code"))
                elif idx == 1: patient_meta[idx].append("Cancer Type" if col == "CANCER_TYPE" else ("Detailed Cancer Type" if col == "CANCER_TYPE_DETAILED" else "OncoTree Code"))
                elif idx == 2: patient_meta[idx].append("STRING")
                elif idx == 3: patient_meta[idx].append("1")
                else: patient_meta[idx].append("NA")

    if patient_data:
        write_cbio_clinical(args.out_patient, patient_headers, patient_meta, patient_data)
    else:
        print("No patient IDs found. Skipping patient clinical file generation.")
        if os.path.exists(args.out_patient):
            os.remove(args.out_patient)

if __name__ == "__main__":
    main()
