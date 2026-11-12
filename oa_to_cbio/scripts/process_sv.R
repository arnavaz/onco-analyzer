#!/usr/bin/env Rscript

# Convert Oncoanalyser LINX SV output to cBioPortal SV format.
# Arguments:
# 1. input_linx_vis_sv_data_tsv
# 2. input_linx_svs_tsv
# 3. output_cbio_sv_tsv
# 4. tumor_sample_id

library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: process_sv.R <input_vis_sv_data_tsv> <input_svs_tsv> <input_bnd_tsv> <input_fusion_csv> <input_sv_vcf> <output_sv_file> <tumor_sample_id>")
}

input_vis_file <- args[1]
input_svs_file <- args[2]
input_bnd_file <- args[3]
input_fusion_file <- args[4]
input_sv_vcf_file <- args[5]
output_file <- args[6]
sample_name <- args[7]

## TEST
if (FALSE) {
  sample_name <- "PANCAN-H-0835-02T-01"
  input_sv_vcf_file <- "/Users/jbruce/Desktop/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/PANCAN-H-0835/PANCAN-H-0835-02T-01/purple/PANCAN-H-0835-02T-01.purple.sv.vcf.gz" 
  input_fusion_file <- "/Users/jbruce/Desktop/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/PANCAN-H-0835/PANCAN-H-0835-02T-01/isofox/PANCAN-H-0835-02T-01-RNA.isf.fusions.csv"
  input_vis_file <- "/Users/jbruce/Desktop/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/PANCAN-H-0835/PANCAN-H-0835-02T-01/linx/somatic_annotations/PANCAN-H-0835-02T-01.linx.vis_sv_data.tsv"
  input_svs_file <- "/Users/jbruce/Desktop/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/PANCAN-H-0835/PANCAN-H-0835-02T-01/linx/somatic_annotations/PANCAN-H-0835-02T-01.linx.svs.tsv"
  input_bnd_file <- "/Users/jbruce/Desktop/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/Oncoanalyser/PANCAN-H-0835/PANCAN-H-0835-02T-01/linx/somatic_annotations/PANCAN-H-0835-02T-01.linx.breakend.tsv"
  output_file <- "/Users/jbruce/Desktop/cluster/projects/mohccn/external_PANCAN/PANCAN-H_CHANGE/cbioportal_input/sv_tmp/PANCAN-H-0835-02T-01.sv.tsv"
}

cat("Processing SVs for", sample_name, "\n")

if (!file.exists(input_vis_file) || !file.exists(input_svs_file) || !file.exists(input_bnd_file)) {
  cat("Warning: Input file does not exist, generating empty SV file.\n")
  data_sv <- data.frame(
    Sample_Id = character(), SV_Status = character(),
    Site1_Hugo_Symbol = character(), Site1_Chromosome = character(), Site1_Contig = character(), Site1_Position = character(),
    Site2_Hugo_Symbol = character(), Site2_Chromosome = character(), Site2_Contig = character(), Site2_Position = character(),
    SV_Class = character()
  )
  write.table(data_sv, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

vis_data <- read.table(input_vis_file, sep = "\t", header = TRUE, comment.char = "", stringsAsFactors = FALSE)
svs_data <- read.table(input_svs_file, sep = "\t", header = TRUE, comment.char = "", stringsAsFactors = FALSE)
bnd_data <- read.table(input_bnd_file, sep = "\t", header = TRUE, comment.char = "", stringsAsFactors = FALSE)

if (nrow(vis_data) == 0 || nrow(svs_data) == 0 || nrow(bnd_data) == 0) {
  data_sv <- data.frame(
    Sample_Id = character(), SV_Status = character(),
    Site1_Hugo_Symbol = character(), Site1_Chromosome = character(), Site1_Contig = character(), Site1_Position = character(),
    Site2_Hugo_Symbol = character(), Site2_Chromosome = character(), Site2_Contig = character(), Site2_Position = character(),
    SV_Class = character()
  )
  write.table(data_sv, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

# Merge them on SvId and svId
# Note: In vis_data it is SvId, in svs_data it is svId
merged <- merge(vis_data, svs_data, by.x = "SvId", by.y = "svId", all.x = TRUE)

# Process bnd_data into wide format (one row per svId) to prevent row explosion
bnd_start <- bnd_data %>% 
  filter(tolower(isStart) == "true") %>% 
  group_by(svId) %>% slice(1) %>% ungroup() %>%
  select(svId, transcriptStart = transcriptId, exonStart = exonUp, junctionContextStart = regionType)

bnd_end <- bnd_data %>% 
  filter(tolower(isStart) == "false") %>% 
  group_by(svId) %>% slice(1) %>% ungroup() %>%
  select(svId, transcriptEnd = transcriptId, exonEnd = exonUp, junctionContextEnd = regionType)

bnd_wide <- merge(bnd_start, bnd_end, by = "svId", all = TRUE)

merged <- merge(merged, bnd_wide, by.x = "SvId", by.y = "svId", all.x = TRUE)

get_col <- function(df, colname) {
  if (colname %in% colnames(df)) {
    return(df[[colname]])
  }
  return(rep(NA, nrow(df)))
}

## Expand rows with multiple genes, spearated by ";" in geneStart or geneEnd
merged_expanded <- merged %>%
  mutate(
    geneStart = str_split(geneStart, ";"),
    geneEnd   = str_split(geneEnd, ";")
  ) %>%
  mutate(
    n = pmax(map_int(geneStart, length), map_int(geneEnd, length)),
    geneStart = map2(geneStart, n, ~ {
      x <- .x
      length(x) <- .y
      x
    }),
    geneEnd = map2(geneEnd, n, ~ {
      x <- .x
      length(x) <- .y
      x
    })
  ) %>%
  unnest(c(geneStart, geneEnd)) %>%
  mutate(across(c(geneStart, geneEnd), ~ na_if(str_trim(.), ""))) %>%
  select(-n)

# ---------------------------------------------------------
# Parse purple.sv.vcf.gz for read support
# ---------------------------------------------------------
if (file.exists(input_sv_vcf_file)) {
  # Read VCF ignoring header lines
  vcf_lines <- readLines(input_sv_vcf_file)
  vcf_lines <- vcf_lines[!grepl("^##", vcf_lines)]
  
  if (length(vcf_lines) > 1) {
    vcf_df <- read.table(text = vcf_lines, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE)
    
    tumor_col <- vcf_df[[ncol(vcf_df)]]
    format_col <- vcf_df$FORMAT
    
    sf_counts <- rep(NA, nrow(vcf_df))
    vf_counts <- rep(NA, nrow(vcf_df))
    
    for (i in seq_along(format_col)) {
      fmt <- unlist(strsplit(format_col[i], ":"))
      val <- unlist(strsplit(tumor_col[i], ":"))
      sf_idx <- match("SF", fmt)
      vf_idx <- match("VF", fmt)
      
      if (!is.na(sf_idx) && length(val) >= sf_idx && val[sf_idx] != ".") sf_counts[i] <- as.numeric(val[sf_idx])
      if (!is.na(vf_idx) && length(val) >= vf_idx && val[vf_idx] != ".") vf_counts[i] <- as.numeric(val[vf_idx])
    }
    
    # Paired-End = Total Variant Fragments (VF) - Split Fragments (SF)
    pe_counts <- vf_counts - sf_counts
    # Protect against negative paired end counts if any artifacts exist
    pe_counts <- ifelse(pe_counts < 0, 0, pe_counts)
    
    vcf_counts <- data.frame(
      vcfId = vcf_df$ID,
      Tumor_Split_Read_Count = sf_counts,
      Tumor_Paired_End_Read_Count = pe_counts,
      stringsAsFactors = FALSE
    )
    
    # Merge the read counts into merged_expanded
    # (merged_expanded has vcfId from svs_data)
    merged_expanded <- merge(merged_expanded, vcf_counts, by = "vcfId", all.x = TRUE)
  }
}

# Safely initialize columns to NA if the VCF was missing or empty
if (!("Tumor_Split_Read_Count" %in% colnames(merged_expanded))) {
  merged_expanded$Tumor_Split_Read_Count <- NA
}
if (!("Tumor_Paired_End_Read_Count" %in% colnames(merged_expanded))) {
  merged_expanded$Tumor_Paired_End_Read_Count <- NA
}

# ---------------------------------------------------------
# Build the final cBioPortal data_sv dataframe
# ---------------------------------------------------------
geneStart <- get_col(merged_expanded, "geneStart")
geneEnd <- get_col(merged_expanded, "geneEnd")
chrStart <- get_col(merged_expanded, "ChrStart")
posStart <- get_col(merged_expanded, "PosStart")
chrEnd <- get_col(merged_expanded, "ChrEnd")
posEnd <- get_col(merged_expanded, "PosEnd")
sv_type <- get_col(merged_expanded, "Type")

# Additional fields for Mavis-like output
transcriptStart <- get_col(merged_expanded, "transcriptStart")
transcriptEnd <- get_col(merged_expanded, "transcriptEnd")
exonStart <- get_col(merged_expanded, "exonStart")
exonEnd <- get_col(merged_expanded, "exonEnd")
juncContextStart <- get_col(merged_expanded, "junctionContextStart")
juncContextEnd <- get_col(merged_expanded, "junctionContextEnd")
orientStart <- get_col(merged_expanded, "OrientStart")
orientEnd <- get_col(merged_expanded, "OrientEnd")
resolvedType <- get_col(merged_expanded, "ResolvedType")

geneStart[is.na(geneStart) | geneStart == "" | geneStart == "None"] <- "Unknown"
geneEnd[is.na(geneEnd) | geneEnd == "" | geneEnd == "None"] <- "Unknown"

data_sv <- data.frame(
  Sample_Id = rep(sample_name, nrow(merged_expanded)),
  SV_Status = rep("SOMATIC", nrow(merged_expanded)),
  
  Site1_Ensembl_Transcript_Id = transcriptStart,
  Site1_Hugo_Symbol = geneStart,
  Site1_Region_Number = exonStart,
  Site1_Region = juncContextStart,
  Site1_Chromosome = chrStart,
  Site1_Contig = NA_character_,
  Site1_Position = posStart,
  Site1_Description = NA_character_,
  
  Site2_Ensembl_Transcript_Id = transcriptEnd,
  Site2_Hugo_Symbol = geneEnd,
  Site2_Region_Number = exonEnd,
  Site2_Region = juncContextEnd,
  Site2_Chromosome = chrEnd,
  Site2_Contig = NA_character_,
  Site2_Position = posEnd,
  Site2_Description = NA_character_,
  
  Site2_Effect_On_Frame = NA_character_,
  SV_Class = toupper(sv_type),
  
  Tumor_Split_Read_Count = merged_expanded$Tumor_Split_Read_Count,
  Tumor_Paired_End_Read_Count = merged_expanded$Tumor_Paired_End_Read_Count,
  
  Event_Info = paste0(sv_type, " event"),
  Connection_Type = paste0(ifelse(orientStart == 1, "3", "5"), "to", ifelse(orientEnd == 1, "5", "3")),
  Breakpoint_Type = ifelse(grepl("INF", resolvedType), "IMPRECISE", "PRECISE"),
  DNA_Support = "Yes",
  SV_Length = ifelse(chrStart == chrEnd & chrEnd != "-1", abs(as.numeric(posEnd) - as.numeric(posStart)), 0),
  NCBI_Build = rep("GRCh38", nrow(merged_expanded))
)
# Replace Unknowns to match your Mavis script style if needed.
# ... feel free to plug the rest of the logic from your Mavis script here!

# Check RNA Support using Fusion File
data_sv$RNA_Support <- "No"

if (file.exists(input_fusion_file)) {
  delim <- if (grepl("\\.tsv$", input_fusion_file, ignore.case = TRUE)) "\t" else ","
  fusion_data <- read.table(input_fusion_file, sep = delim, header = TRUE, comment.char = "", stringsAsFactors = FALSE)
  if (nrow(fusion_data) > 0) {
    
    # Prepare fusion pairs, keeping JuncType info for the frame
    fus_pairs <- fusion_data %>%
      select(g1 = GeneNameUp, g2 = GeneNameDown, JuncTypeUp, JuncTypeDown) %>%
      filter(g1 != "Unknown", g1 != "None", g1 != "", !is.na(g1)) %>%
      filter(g2 != "Unknown", g2 != "None", g2 != "", !is.na(g2)) %>%
      separate_rows(g1, sep = ";") %>%
      separate_rows(g2, sep = ";") %>%
      group_by(g1, g2) %>%
      summarise(
        JuncTypeUp = paste(unique(JuncTypeUp), collapse = ","),
        JuncTypeDown = paste(unique(JuncTypeDown), collapse = ","),
        .groups = "drop"
      ) %>%
      distinct()
      
    # Add an index to track rows
    data_sv <- data_sv %>% mutate(.sv_idx = row_number())
    
    # Direction 1: SV Site1 matches Fusion g1, Site2 matches Fusion g2
    match1 <- data_sv %>%
      inner_join(fus_pairs, by = c("Site1_Hugo_Symbol" = "g1", "Site2_Hugo_Symbol" = "g2")) %>%
      mutate(JuncInfo = paste(JuncTypeUp, JuncTypeDown, sep = "-")) %>%
      select(.sv_idx, JuncInfo)
      
    # Direction 2: SV Site1 matches Fusion g2, Site2 matches Fusion g1
    match2 <- data_sv %>%
      inner_join(fus_pairs, by = c("Site1_Hugo_Symbol" = "g2", "Site2_Hugo_Symbol" = "g1")) %>%
      mutate(JuncInfo = paste(JuncTypeDown, JuncTypeUp, sep = "-")) %>%
      select(.sv_idx, JuncInfo)
      
    matched_df <- bind_rows(match1, match2) %>%
      group_by(.sv_idx) %>%
      summarise(JuncInfo = paste(unique(JuncInfo), collapse = "; "), .groups = "drop")
    
    # Update RNA_Support and map frame labels
    data_sv <- data_sv %>%
      left_join(matched_df, by = ".sv_idx") %>%
      mutate(
        RNA_Support = if_else(!is.na(JuncInfo), "Yes", "No"),
        # Map Isofox junction types to cBioPortal accepted frame labels
        # KNOWN implies an annotated fusion transcript, so we can reasonably infer in-frame
        # CANONICAL just means canonical splice sites (GT-AG), so the frame is actually unknown
        Site2_Effect_On_Frame = case_when(
          is.na(JuncInfo) ~ NA_character_,
          grepl("KNOWN", JuncInfo) ~ "in-frame",
          TRUE ~ NA_character_
        )
      ) %>%
      select(-.sv_idx, -JuncInfo)
  }
}

# Filter out structural variants where both endpoints are Unknown genes
data_sv <- data_sv[!(data_sv$Site1_Hugo_Symbol == "Unknown" & data_sv$Site2_Hugo_Symbol == "Unknown"), ]

# Remove duplicates based on standard cBioPortal SV key columns
key_cols <- c("Sample_Id", "Site1_Hugo_Symbol", "Site1_Chromosome", "Site1_Position", "Site2_Hugo_Symbol", "Site2_Chromosome", "Site2_Position", "SV_Class")
data_sv <- data_sv[!duplicated(data_sv[, key_cols]), ]

write.table(data_sv, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote SV output to", output_file, "\n")
