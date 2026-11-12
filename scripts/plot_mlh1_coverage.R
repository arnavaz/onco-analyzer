#!/usr/bin/env Rscript

# plot_mlh1_coverage.R
# Usage: Rscript plot_mlh1_coverage.R [results_dir] [output_base_path]

# Helper to check and load packages
load_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "is required but not installed. Attempting to install..."))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

load_library("ggplot2")

args <- commandArgs(trailingOnly = TRUE)
results_dir <- "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/"
output_base_path <- "/Users/jbruce/Desktop/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/mlh1_binned_coverage_plot"

if (length(args) >= 1) {
  results_dir <- args[1]
}
if (length(args) >= 2) {
  output_base_path <- args[2]
}

sample_files <- list(
  "UHN-27682" = list(
    "cn" = file.path(results_dir, "DYNA-UHN-27682/DYNA-UHN-27682-FFPE-T-DNA-U_S53/purple/DYNA-UHN-27682-FFPE-T-DNA-U_S53.purple.target_region_cn.tsv"),
    "cobalt" = file.path(results_dir, "DYNA-UHN-27682/DYNA-UHN-27682-FFPE-T-DNA-U_S53/cobalt/DYNA-UHN-27682-FFPE-T-DNA-U_S53.cobalt.ratio.tsv.gz")
  ),
  "UHN-31882" = list(
    "cn" = file.path(results_dir, "DYNA-UHN-31882/DYNA-UHN-31882-FFPE-T-DNA-U_S44/purple/DYNA-UHN-31882-FFPE-T-DNA-U_S44.purple.target_region_cn.tsv"),
    "cobalt" = file.path(results_dir, "DYNA-UHN-31882/DYNA-UHN-31882-FFPE-T-DNA-U_S44/cobalt/DYNA-UHN-31882-FFPE-T-DNA-U_S44.cobalt.ratio.tsv.gz")
  )
)

gene_starts <- c("MLH1" = 36993226, "BAP1" = 52402288, "MITF" = 69739598)

# Parse exon coordinates from bedRegions to draw the gene structure track at the bottom
exons_list <- list()

for (sample_name in names(sample_files)) {
  cn_path <- sample_files[[sample_name]]$cn
  if (!file.exists(cn_path)) next
  
  df <- read.delim(cn_path, sep = "\t", stringsAsFactors = FALSE)
  df_chr3 <- df[df$chromosome == "chr3", ]
  if (nrow(df_chr3) == 0) next
  
  df_chr3$Gene <- NA
  df_chr3$Gene[grepl("MLH1", df_chr3$bedRegions, ignore.case = TRUE)] <- "MLH1"
  df_chr3$Gene[grepl("BAP1", df_chr3$bedRegions, ignore.case = TRUE)] <- "BAP1"
  df_chr3$Gene[grepl("MITF", df_chr3$bedRegions, ignore.case = TRUE)] <- "MITF"
  
  df_chr3 <- df_chr3[!is.na(df_chr3$Gene), ]
  if (nrow(df_chr3) == 0) next
  
  for (i in 1:nrow(df_chr3)) {
    region_str <- df_chr3$bedRegions[i]
    gene_name <- df_chr3$Gene[i]
    
    parts <- strsplit(region_str, ";")[[1]]
    for (part in parts) {
      if (grepl(":", part)) {
        coord_part <- strsplit(part, ":")[[1]][2]
        start_end <- strsplit(coord_part, "-")[[1]]
        start_val <- as.numeric(start_end[1])
        end_val <- as.numeric(start_end[2])
        
        exons_list <- rbind(exons_list, data.frame(
          Gene = gene_name,
          Start = start_val,
          End = end_val
        ))
      }
    }
  }
}

exons_df <- unique(exons_list)
exons_df$xmin <- (exons_df$Start - gene_starts[exons_df$Gene]) / 1000
exons_df$xmax <- (exons_df$End - gene_starts[exons_df$Gene]) / 1000
exons_df$ymin <- -0.4
exons_df$ymax <- -0.2

# Intron lines
introns_list <- list()
for (g in unique(exons_df$Gene)) {
  g_exons <- exons_df[exons_df$Gene == g, ]
  introns_list <- rbind(introns_list, data.frame(
    Gene = g,
    xmin = min(g_exons$xmin),
    xmax = max(g_exons$xmax),
    y = -0.3
  ))
}
introns_df <- data.frame(introns_list)

# Process each sample separately
for (sample_name in names(sample_files)) {
  cn_path <- sample_files[[sample_name]]$cn
  cobalt_path <- sample_files[[sample_name]]$cobalt
  
  if (!file.exists(cn_path) || !file.exists(cobalt_path)) {
    warning(paste("Warning: Files missing for", sample_name))
    next
  }
  
  message(paste("Processing copy number data for", sample_name))
  
  # Read COBALT ratios and filter for targeted exon bins on chr3
  df_cobalt <- read.delim(gzfile(cobalt_path), sep = "\t", stringsAsFactors = FALSE)
  df_cobalt <- df_cobalt[df_cobalt$chromosome == "chr3", ]
  df_cobalt <- df_cobalt[df_cobalt$referenceGCRatio > 0 & df_cobalt$referenceGCDiploidRatio > 0, ]
  
  # Map COBALT bins to genes by coordinates
  df_cobalt$Gene <- NA
  df_cobalt$Gene[df_cobalt$position >= 36990000 & df_cobalt$position <= 37051000] <- "MLH1"
  df_cobalt$Gene[df_cobalt$position >= 52402000 & df_cobalt$position <= 52410000] <- "BAP1"
  df_cobalt$Gene[df_cobalt$position >= 69730000 & df_cobalt$position <= 69970000] <- "MITF"
  df_cobalt <- df_cobalt[!is.na(df_cobalt$Gene), ]
  
  # Calculate binned copy numbers using COBALT GC-diploid ratio directly
  df_cobalt$Normal_CN_Binned <- df_cobalt$referenceGCDiploidRatio * 2
  df_cobalt$Tumour_CN_Binned <- (df_cobalt$tumorReadDepth / df_cobalt$referenceReadDepth) * df_cobalt$referenceGCDiploidRatio * 2
  
  # Convert position to relative kb from gene start
  df_cobalt$Relative_Pos_Kb <- (df_cobalt$position - gene_starts[df_cobalt$Gene]) / 1000
  df_cobalt$Gene <- factor(df_cobalt$Gene, levels = c("MLH1", "BAP1", "MITF"))
  
  # Reshape data to long format for clean plotting of both normal and tumor binned points
  plot_df <- data.frame(
    Relative_Pos_Kb = rep(df_cobalt$Relative_Pos_Kb, 2),
    Gene = rep(df_cobalt$Gene, 2),
    Copy_Number = c(df_cobalt$Normal_CN_Binned, df_cobalt$Tumour_CN_Binned),
    Sample_Type = c(rep("Normal (Germline)", nrow(df_cobalt)), rep("Tumour (Somatic Observed)", nrow(df_cobalt)))
  )
  plot_df$Sample_Type <- factor(plot_df$Sample_Type, levels = c("Normal (Germline)", "Tumour (Somatic Observed)"))
  
  exons_sub <- exons_df[exons_df$Gene %in% c("MLH1", "BAP1", "MITF"), ]
  exons_sub$Gene <- factor(exons_sub$Gene, levels = c("MLH1", "BAP1", "MITF"))
  introns_sub <- introns_df[introns_df$Gene %in% c("MLH1", "BAP1", "MITF"), ]
  introns_sub$Gene <- factor(introns_sub$Gene, levels = c("MLH1", "BAP1", "MITF"))
  
  # Create plot
  p <- ggplot() +
    # Gene track at bottom
    geom_segment(data = introns_sub, aes(x = xmin, xend = xmax, y = y, yend = y), color = "gray20", linewidth = 0.5) +
    geom_rect(data = exons_sub, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#475569", color = "gray10", size = 0.2) +
    
    # Normal CN (Germline) and Tumour CN (Somatic) binned points
    geom_point(data = plot_df, aes(x = Relative_Pos_Kb, y = Copy_Number, color = Sample_Type), size = 2.5, alpha = 0.8) +
    
    facet_grid(. ~ Gene, scales = "free_x") +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_hline(yintercept = 1.0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0.0, linetype = "solid", color = "red3", linewidth = 0.4) +
    
    scale_color_manual(values = c(
      "Normal (Germline)" = "#377EB8",           # Blue points
      "Tumour (Somatic Observed)" = "#E41A1C"    # Red points
    )) +
    
    scale_y_continuous(limits = c(-0.5, 4.0), breaks = c(-0.3, 0, 1, 2, 3, 4), labels = c("Gene Track", "0", "1", "2", "3", "4")) +
    labs(
      title = paste("Exon-level Copy Number Profile:", sample_name),
      subtitle = "COBALT GC-diploid binned normal and observed tumor copy numbers",
      x = "Relative Position from Gene Start (kb)",
      y = "Copy Number (GC-Diploid Normalized)",
      color = "Sample Channel"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray30", size = 10),
      strip.background = element_rect(fill = "#1e293b", color = NA),
      strip.text = element_text(color = "white", face = "bold", size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
    
  plot_png_path <- paste0(output_base_path, "_", sample_name, ".png")
  plot_pdf_path <- paste0(output_base_path, "_", sample_name, ".pdf")
  
  message(paste("Saving PNG to:", plot_png_path))
  ggsave(plot_png_path, plot = p, width = 11, height = 5.5, dpi = 300)
  
  message(paste("Saving PDF to:", plot_pdf_path))
  ggsave(plot_pdf_path, plot = p, width = 11, height = 5.5)
}

message("Individual binned coverage plots complete!")
