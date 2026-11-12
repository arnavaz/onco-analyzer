#!/usr/bin/env Rscript

# Usage:
#   Rscript plot_linearized_input.R <cobalt_ratio_tsv> <cobalt_pcf> <amber_baf_tsv_gz> <output_png> [sample_name]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript plot_linearized_input.R <cobalt_ratio_tsv> <cobalt_pcf> <amber_baf_tsv_gz> <output_png> [sample_name]\n")
  quit(status = 1)
}

cobalt_ratio_file <- args[1]
cobalt_pcf_file   <- args[2]
amber_baf_file    <- args[3]
output_png        <- args[4]
sample_name       <- ifelse(length(args) >= 5, args[5], "Sample")


cobalt_ratio_file <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238/S25-22052/cobalt/S25-22052.cobalt.ratio.tsv"
cobalt_pcf_file   <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238/S25-22052/cobalt/S25-22052.cobalt.ratio.pcf"
amber_baf_file    <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238/S25-22052/amber/S25-22052.amber.baf.tsv.gz"
output_png        <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238_S25-22052.linearized.png"
sample_name       <- "CHM2-01-0238_S25-22052"


library(ggplot2)
library(dplyr)
library(grid)

# Define GRCh38 chromosome lengths and offsets
chr_lengths <- c(
  "chr1" = 248956422, "chr2" = 242193529, "chr3" = 198295559, "chr4" = 190214555,
  "chr5" = 181538259, "chr6" = 170805979, "chr7" = 159345973, "chr8" = 145138636,
  "chr9" = 138394717, "chr10" = 133797422, "chr11" = 135086622, "chr12" = 133275309,
  "chr13" = 114364328, "chr14" = 107043718, "chr15" = 101991189, "chr16" = 90338345,
  "chr17" = 83257522, "chr18" = 80373285, "chr19" = 58617616, "chr20" = 64444167,
  "chr21" = 46709983, "chr22" = 50818468, "chrX" = 156040895, "chrY" = 57227415
)
chr_offsets <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
names(chr_offsets) <- names(chr_lengths)

boundaries <- chr_offsets + chr_lengths
boundaries <- boundaries[-length(boundaries)]

chrom_info <- data.frame(
  chrom = names(chr_lengths),
  length = chr_lengths,
  offset = chr_offsets,
  stringsAsFactors = FALSE
) %>%
  mutate(
    center = offset + length / 2,
    label = gsub("chr", "", chrom)
  )

# --- 1. Load and bin COBALT ratios ---
cat("Loading and binning COBALT ratio data...\n")
cobalt <- read.table(cobalt_ratio_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Filter for standard chromosomes
cobalt <- cobalt %>% filter(chromosome %in% names(chr_lengths))

# GC-normalized ratio is in tumorGCRatio column. Filter out unnormalized bins (-1) or zero depth
cobalt_filtered <- cobalt %>% filter(tumorGCRatio >= 0, referenceReadDepth > 0)

# Bin to 100kb windows to keep plotting fast and clean
bin_size <- 100000000
cobalt_binned <- cobalt_filtered %>%
  mutate(bin = floor(position / bin_size) * bin_size) %>%
  group_by(chromosome, bin) %>%
  summarise(
    tumorGCRatio = median(tumorGCRatio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(position = bin) %>%
  mutate(
    pos_cum = position + chr_offsets[chromosome],
    log2_ratio = log2(pmax(0.001, tumorGCRatio))
  )

# --- 2. Load COBALT PCF segments ---
cat("Loading COBALT segments...\n")
segments <- read.table(cobalt_pcf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
segments <- segments %>% 
  filter(chrom %in% names(chr_lengths)) %>%
  mutate(
    start_cum = start.pos + chr_offsets[chrom],
    end_cum = end.pos + chr_offsets[chrom]
  )

# --- 3. Load and downsample AMBER BAF ---
cat("Loading and downsampling AMBER BAF data...\n")
# Open gzipped file
con <- gzfile(amber_baf_file)
baf <- read.table(con, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# close(con)

baf <- baf %>% filter(chromosome %in% names(chr_lengths))

# Downsample BAF to 100k points for performance
set.seed(42)
if (nrow(baf) > 100000) {
  baf <- baf[sample(1:nrow(baf), 100000), ]
}

baf <- baf %>%
  mutate(pos_cum = position + chr_offsets[chromosome])

# --- 4. Plotting ---
cat("Generating plots...\n")

# A. COBALT Ratio Plot
p1 <- ggplot() +
  # Background vertical gridlines for chromosomes
  geom_vline(xintercept = boundaries, color = "grey85", linetype = "dashed") +
  # Ratio points
  geom_point(data = cobalt_binned, aes(x = pos_cum, y = log2_ratio), color = "grey50", alpha = 0.5, size = 0.4) +
  # PCF Segment lines
  geom_segment(data = segments, aes(x = start_cum, xend = end_cum, y = mean, yend = mean), color = "darkblue", size = 1.0) +
  scale_x_continuous(breaks = chrom_info$center, labels = chrom_info$label, expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(-2.0, 2.0), breaks = seq(-2, 2, 1)) +
  labs(
    title = paste0(sample_name, " - Linearized Genomic Profile"),
    y = "Log2 Ratio (Tumor/Ref)",
    x = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

# B. AMBER BAF Plot
p2 <- ggplot() +
  geom_vline(xintercept = boundaries, color = "grey85", linetype = "dashed") +
  geom_point(data = baf, aes(x = pos_cum, y = tumorBAF), color = "darkorchid", alpha = 0.4, size = 0.3) +
  scale_x_continuous(breaks = chrom_info$center, labels = chrom_info$label, expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.25)) +
  labs(
    y = "Tumor BAF",
    x = "Chromosome"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save combined layout

output_png        <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238_S25-22052.linearized_bin_100MB.png"


cat("Saving combined plot to", output_png, "...\n")
png(output_png, width = 1600, height = 900, res = 150)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 1), "null"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

cat("Done!\n")

output_pdf        <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/CHARM_WGTS/Oncoanalyser/CHM2-01-0238_S25-22052.linearized_bin_100MB_ratio.pdf"

cat("Saving just log plot to", output_pdf, "...\n")
pdf(output_pdf, width = 10, height = 5)
print(p1 +scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.1)))
dev.off()

cat("Done!\n")

