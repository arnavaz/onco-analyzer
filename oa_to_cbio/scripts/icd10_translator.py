import os
import json
import re

# Load official NIH NLM ICD-10-CM codes dictionary
JSON_PATH = os.path.join(os.path.dirname(__file__), 'icd10_codes.json')
OFFICIAL_ICD10_CODES = {}

if os.path.exists(JSON_PATH):
    try:
        with open(JSON_PATH, 'r') as f:
            OFFICIAL_ICD10_CODES = json.load(f)
    except Exception as e:
        pass

ONCOTREE_NAMES_PATH = os.path.join(os.path.dirname(__file__), 'oncotree_code_names.json')
ONCOTREE_NAME_MAP = {}
if os.path.exists(ONCOTREE_NAMES_PATH):
    try:
        with open(ONCOTREE_NAMES_PATH, 'r') as f:
            ONCOTREE_NAME_MAP = json.load(f)
    except Exception:
        pass


def get_oncotree_name(code):
    """Return the canonical OncoTree name for a code, or None if unknown."""
    if not code:
        return None
    return ONCOTREE_NAME_MAP.get(str(code).strip().upper())

# Official WHO ICD-10 3-character category to Broad cBioPortal Cancer Type mapping
CATEGORY_BROAD_MAP = {
    'C00': 'Head and Neck Cancer', 'C01': 'Head and Neck Cancer', 'C02': 'Head and Neck Cancer', 'C03': 'Head and Neck Cancer',
    'C04': 'Head and Neck Cancer', 'C05': 'Head and Neck Cancer', 'C06': 'Head and Neck Cancer', 'C07': 'Head and Neck Cancer',
    'C08': 'Head and Neck Cancer', 'C09': 'Head and Neck Cancer', 'C10': 'Head and Neck Cancer', 'C11': 'Head and Neck Cancer',
    'C12': 'Head and Neck Cancer', 'C13': 'Head and Neck Cancer', 'C14': 'Head and Neck Cancer',
    'C15': 'Esophagogastric Cancer', 'C16': 'Esophagogastric Cancer',
    'C17': 'Gastrointestinal Cancer', 'C18': 'Colorectal Cancer', 'C19': 'Colorectal Cancer', 'C20': 'Colorectal Cancer',
    'C21': 'Anal Cancer', 'C22': 'Hepatobiliary Cancer', 'C23': 'Hepatobiliary Cancer', 'C24': 'Hepatobiliary Cancer',
    'C25': 'Pancreatic Cancer', 'C26': 'Gastrointestinal Cancer',
    'C30': 'Head and Neck Cancer', 'C31': 'Head and Neck Cancer', 'C32': 'Head and Neck Cancer',
    'C33': 'Non-Small Cell Lung Cancer', 'C34': 'Non-Small Cell Lung Cancer', 'C37': 'Thymic Tumor',
    'C38': 'Pleural Mesothelioma', 'C39': 'Thoracic Cancer',
    'C40': 'Bone Cancer', 'C41': 'Bone Cancer',
    'C43': 'Melanoma', 'C44': 'Skin Cancer, Non-Melanoma',
    'C45': 'Mesothelioma', 'C46': 'Sarcoma', 'C47': 'Nerve Sheath Tumor', 'C48': 'Peritoneal Cancer', 'C49': 'Soft Tissue Sarcoma',
    'C50': 'Breast Cancer',
    'C51': 'Gynecologic Cancer', 'C52': 'Gynecologic Cancer', 'C53': 'Cervical Cancer', 'C54': 'Uterine Cancer', 'C55': 'Uterine Cancer',
    'C56': 'Ovarian Cancer', 'C57': 'Gynecologic Cancer', 'C58': 'Gynecologic Cancer',
    'C60': 'Penile Cancer', 'C61': 'Prostate Cancer', 'C62': 'Testicular Cancer', 'C63': 'Male Genital Cancer',
    'C64': 'Renal Cell Carcinoma', 'C65': 'Urothelial Cancer', 'C66': 'Urothelial Cancer', 'C67': 'Bladder Cancer', 'C68': 'Urothelial Cancer',
    'C69': 'Ocular Melanoma', 'C70': 'CNS/Brain Cancer', 'C71': 'Glioma', 'C72': 'CNS/Brain Cancer',
    'C73': 'Thyroid Cancer', 'C74': 'Adrenocortical Carcinoma', 'C75': 'Endocrine Cancer',
    'C7A': 'Neuroendocrine Tumor', 'C7B': 'Neuroendocrine Tumor',
    'C76': 'Cancer of Unknown Primary', 'C77': 'Cancer of Unknown Primary', 'C78': 'Cancer of Unknown Primary',
    'C79': 'Cancer of Unknown Primary', 'C80': 'Cancer of Unknown Primary',
    'C81': 'Hodgkin Lymphoma', 'C82': 'Non-Hodgkin Lymphoma', 'C83': 'Non-Hodgkin Lymphoma', 'C84': 'Mature T/NK-Cell Lymphomas',
    'C85': 'Non-Hodgkin Lymphoma', 'C86': 'Mature T/NK-Cell Lymphomas', 'C88': 'Lymphoproliferative Disorder',
    'C90': 'Multiple Myeloma', 'C91': 'Lymphoid Leukemia', 'C92': 'Myeloid Leukemia', 'C93': 'Leukemia', 'C94': 'Leukemia',
    'C95': 'Leukemia', 'C96': 'Hematologic Malignancy',
    'D00': 'Carcinoma in Situ', 'D01': 'Carcinoma in Situ', 'D02': 'Carcinoma in Situ', 'D03': 'Melanoma in Situ',
    'D04': 'Carcinoma in Situ', 'D05': 'Breast Carcinoma in Situ', 'D06': 'Cervical Carcinoma in Situ',
    'D07': 'Carcinoma in Situ', 'D09': 'Carcinoma in Situ'
}

# Standard MSK OncoTree Codes for ICD-10 categories
ONCOTREE_CODE_MAP = {
    'C00': 'HN', 'C01': 'HN', 'C02': 'HN', 'C03': 'HN', 'C04': 'HN', 'C05': 'HN', 'C06': 'HN', 'C07': 'HN',
    'C08': 'HN', 'C09': 'HN', 'C10': 'HN', 'C11': 'HN', 'C12': 'HN', 'C13': 'HN', 'C14': 'HN',
    'C15': 'ESCA', 'C16': 'STAD', 'C17': 'SBN', 'C18': 'COAD', 'C19': 'READ', 'C20': 'READ',
    'C21': 'ANAL', 'C22': 'LIHC', 'C23': 'GBC', 'C24': 'IHCH', 'C25': 'PAAD', 'C26': 'GI',
    'C30': 'HN', 'C31': 'HN', 'C32': 'HN', 'C33': 'NSCLC', 'C34': 'NSCLC', 'C37': 'THYM',
    'C38': 'MESO', 'C39': 'THOR', 'C40': 'BONE', 'C41': 'BONE', 'C43': 'MEL', 'C44': 'SKIN',
    'C45': 'MESO', 'C46': 'KAPOSI', 'C47': 'PNST', 'C48': 'PERI', 'C49': 'STS', 'C50': 'BRCA',
    'C51': 'VULVA', 'C52': 'VAGINA', 'C53': 'CESC', 'C54': 'UCEC', 'C55': 'UTEC', 'C56': 'OV',
    'C57': 'GYN', 'C58': 'GTD', 'C60': 'PENIS', 'C61': 'PRAD', 'C62': 'TEST', 'C63': 'MGU',
    'C64': 'KIRC', 'C65': 'UTUC', 'C66': 'UTUC', 'C67': 'BLCA', 'C68': 'UTUC', 'C69': 'OM',
    'C70': 'CNS', 'C71': 'GLIOMA', 'C72': 'CNS', 'C73': 'THCA', 'C74': 'ACC', 'C75': 'ENDOCRINE',
    'C7A': 'NET', 'C7B': 'NET', 'C76': 'CUP', 'C77': 'CUP', 'C78': 'CUP', 'C79': 'CUP', 'C80': 'CUP',
    'C81': 'HD', 'C82': 'NHL', 'C83': 'NHL', 'C84': 'PTCL', 'C85': 'NHL', 'C86': 'PTCL', 'C88': 'LPD',
    'C90': 'MM', 'C91': 'ALL', 'C92': 'AML', 'C93': 'LEUK', 'C94': 'LEUK', 'C95': 'LEUK', 'C96': 'HEM'
}

def translate_icd10(val):
    """
    Translates an ICD-10-CM diagnosis code string (e.g. 'C18.9' or 'C34.90')
    into a tuple of (broad_cancer_type, official_detailed_description, oncotree_code).
    If val is not an ICD-10 code, returns (val, val, 'NA').
    """
    if not val or str(val).strip().upper() in ['NA', 'UNKNOWN', 'NULL', '']:
        return 'Cancer of Unknown Primary', 'Cancer of Unknown Primary', 'CUP'
    
    val_clean = str(val).strip()
    match = re.search(r'([A-Z][0-9]{2}(\.[0-9A-Z]+)?)', val_clean)
    if not match:
        return val_clean, val_clean, 'NA'
        
    full_code = match.group(1)
    cat_code = full_code.split('.')[0]
    
    # Retrieve official CDC/NIH NLM description if available
    official_desc = OFFICIAL_ICD10_CODES.get(full_code, OFFICIAL_ICD10_CODES.get(cat_code, val_clean))
    
    # Retrieve broad cBioPortal cancer type & OncoTree code
    broad_type = CATEGORY_BROAD_MAP.get(cat_code, val_clean)
    oncotree_code = ONCOTREE_CODE_MAP.get(cat_code, 'NA')
    
    return broad_type, official_desc, oncotree_code
