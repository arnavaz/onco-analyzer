#!/usr/bin/env Rscript

# plot_msi_unstable_sites.R
# Usage: Rscript plot_msi_unstable_sites.R [path_to_final_summary.tsv] [output_plot_path]

# Helper to check and load packages
load_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "is required but not installed. Attempting to install..."))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

load_library("ggplot2")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

tsv_path <- NULL
output_path <- NULL

if (length(args) >= 1) {
  tsv_path <- args[1]
}
if (length(args) >= 2) {
  output_path <- args[2]
}

# Auto-locate final_summary.tsv if not specified
if (is.null(tsv_path)) {
  candidates <- c(
    "final_summary.tsv",
    "summary/final_summary.tsv",
    "../final_summary.tsv"
  )
  for (cand in candidates) {
    if (file.exists(cand)) {
      tsv_path <- cand
      break
    }
  }
}

if (is.null(tsv_path) || !file.exists(tsv_path)) {
  stop("Error: Input summary TSV file not found. Please provide the path to final_summary.tsv.")
}

if (is.null(output_path)) {
  output_path <- "msi_unstable_sites_plot.png"
}

message(paste("Reading summary TSV from:", tsv_path))
data <- read.delim(tsv_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Check required columns
req_cols <- c("Sample", "Tumour.Fraction", "MSI.Unstable.Sites.Proportion")
missing_cols <- setdiff(req_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("Error: Missing required columns in TSV:", paste(missing_cols, collapse = ", ")))
}

# Convert data types and handle NAs
data$Tumour.Fraction <- as.numeric(as.character(data$Tumour.Fraction))
data$MSI.Unstable.Sites.Proportion <- as.numeric(as.character(data$MSI.Unstable.Sites.Proportion))

# Filter out rows that do not have valid MSI proportions
original_rows <- nrow(data)
data <- data[!is.na(data$MSI.Unstable.Sites.Proportion), ]
message(paste("Filtered out", original_rows - nrow(data), "rows with missing MSI proportion data."))

if (nrow(data) == 0) {
  stop("Error: No samples with valid MSI proportion data to plot.")
}

# Group samples dynamically: Normals, cfDNA, Tumours
data$Group <- "Tumours"

# Check for cell-free DNA
cfdna_indices <- grepl("cfDNA|liquid|cf-DNA|-P-T-|-P-", data$Sample, ignore.case = TRUE)
data$Group[cfdna_indices] <- "cfDNA"

# Check for normal background samples
normal_indices <- grepl("normal|buffy|control|-N-", data$Sample, ignore.case = TRUE)
data$Group[normal_indices] <- "Normals"

data$Group <- factor(data$Group, levels = c("Tumours", "cfDNA", "Normals"))

data$MSI.Status <- toupper(trimws(as.character(data$MSI.Status)))
data$MSI.Status[is.na(data$MSI.Status) | data$MSI.Status == ""] <- "MSS"
data$MSI.Status <- factor(data$MSI.Status, levels = c("MSI", "MSS"))

# Sort samples: descending order of unstable sites within each Group
data <- data[order(data$Group, -data$MSI.Unstable.Sites.Proportion), ]
data$Sample <- factor(data$Sample, levels = unique(data$Sample))

# Create the plot
p <- ggplot(data, aes(x = Sample, y = MSI.Unstable.Sites.Proportion, fill = MSI.Status)) +
  geom_col(color = "black", size = 0.3, width = 0.75) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("MSI" = "#E41A1C", "MSS" = "#94a3b8")) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, 0.2), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Proportion of MS sites that are unstable",
    fill = "MSI Status"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    strip.background = element_rect(fill = "#377EB8", color = NA),
    strip.text = element_text(color = "white", face = "bold", size = 11),
    legend.position = "right",
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )

# If there are many samples, hide the x-axis sample IDs to keep the layout clean
if (nrow(data) > 25) {
  message("Large cohort detected (> 25 samples). Hiding individual x-axis sample labels for readability.")
  p <- p + theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
}

# Save the plot
message(paste("Saving plot to:", output_path))
ggsave(output_path, plot = p, width = 14, height = 5)
message("Plot generation complete!")
