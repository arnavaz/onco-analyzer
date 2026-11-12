library(HMMcopy)
library(GenomicRanges)
library(optparse)
library(ichorCNA)

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

option_list <- list(
  make_option(c("-f", "--filelist"), type = "character", help = "List of wig files."),
  make_option(c("-o", "--outfile"), type = "character", help = "Output file (base path, will append .rds)."),
  make_option(c("-g", "--genomeBuild"), type = "character", default = "hg38", help = "hg19 or hg38. [Default: %default]"),
  make_option(c("-w", "--windowSize"), type = "character", default = "1000kb", help = "Window size (e.g. 1000kb). [Default: %default]"),
  make_option(c("--genomeStyle"), type = "character", default = "UCSC", help = "UCSC or NCBI. [Default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Resolve package file names based on genome build
build_for_wig <- opt$genomeBuild
if (build_for_wig == "GRCh37") {
  build_for_wig <- "hg19"
} else if (build_for_wig == "GRCh38") {
  build_for_wig <- "hg38"
}

centromere_file <- if (build_for_wig == "hg38") {
  "GRCh38.GCA_000001405.2_centromere_acen.txt"
} else if (build_for_wig == "hg19") {
  "GRCh37.p13_centromere_UCSC-gapTable.txt"
} else {
  paste0(opt$genomeBuild, "_centromeres.txt")
}

# Dynamically find resources within the installed ichorCNA R package
gcWig <- system.file("extdata", paste0("gc_", build_for_wig, "_", opt$windowSize, ".wig"), package = "ichorCNA")
mapWig <- system.file("extdata", paste0("map_", build_for_wig, "_", opt$windowSize, ".wig"), package = "ichorCNA")
centromere <- system.file("extdata", centromere_file, package = "ichorCNA")

if (gcWig == "") {
  stop(paste0("Could not find GC Wig file for build: ", build_for_wig, ", window size: ", opt$windowSize))
}
if (mapWig == "") {
  stop(paste0("Could not find Map Wig file for build: ", build_for_wig, ", window size: ", opt$windowSize))
}
if (centromere == "") {
  stop(paste0("Could not find Centromere file: ", centromere_file))
}

message("Using GC Wig: ", gcWig)
message("Using Map Wig: ", mapWig)
message("Using Centromere file: ", centromere)

chrs <- if (opt$genomeBuild == "hg38") {
  c(paste0("chr", 1:22), "chrX")
} else {
  c(1:22, "X")
}
chrNormalize <- if (opt$genomeBuild == "hg38") {
  paste0("chr", 1:22)
} else {
  1:22
}

seqlevelsStyle(chrs) <- opt$genomeStyle
seqlevelsStyle(chrNormalize) <- opt$genomeStyle

centromere_df <- read.delim(centromere, header=TRUE, stringsAsFactors=FALSE, sep="\t")
seqlevelsStyle(centromere_df$Chr) <- opt$genomeStyle

files <- read.delim(opt$filelist, header = FALSE, stringsAsFactors=FALSE, sep ="\t")[, 1]

# Load GC and Map reference WIGs once outside the loop
gc <- wigToGRanges(gcWig)
map <- wigToGRanges(mapWig)

# Align seqlevels style manually for gc and map
if (opt$genomeStyle == "UCSC") {
  # gc
  levels_gc <- seqlevels(gc)
  if (!any(grepl("^chr", levels_gc))) {
    new_levels <- paste0("chr", levels_gc)
    new_levels[new_levels == "chrMT"] <- "chrM"
    names(new_levels) <- levels_gc
    gc <- renameSeqlevels(gc, new_levels)
  }
  # map
  levels_map <- seqlevels(map)
  if (!any(grepl("^chr", levels_map))) {
    new_levels <- paste0("chr", levels_map)
    new_levels[new_levels == "chrMT"] <- "chrM"
    names(new_levels) <- levels_map
    map <- renameSeqlevels(map, new_levels)
  }
} else if (opt$genomeStyle == "NCBI") {
  # gc
  levels_gc <- seqlevels(gc)
  if (any(grepl("^chr", levels_gc))) {
    new_levels <- gsub("^chr", "", levels_gc)
    new_levels[new_levels == "M"] <- "MT"
    names(new_levels) <- levels_gc
    gc <- renameSeqlevels(gc, new_levels)
  }
  # map
  levels_map <- seqlevels(map)
  if (any(grepl("^chr", levels_map))) {
    new_levels <- gsub("^chr", "", levels_map)
    new_levels[new_levels == "M"] <- "MT"
    names(new_levels) <- levels_map
    map <- renameSeqlevels(map, new_levels)
  }
}

# Ensure gc and map have all chrs in seqlevels
seqlevels(gc) <- unique(c(seqlevels(gc), chrs))
seqlevels(map) <- unique(c(seqlevels(map), chrs))

normalGR <- NULL

for (i in 1:length(files)){
  sid <- gsub(".wig", "", basename(files[i]))
  message("Loading and correcting normal file: ", files[i])
  normal_reads <- wigToGRanges(files[i])
  
  # Align seqlevels style manually for normal_reads
  if (opt$genomeStyle == "UCSC") {
    levels_reads <- seqlevels(normal_reads)
    if (!any(grepl("^chr", levels_reads))) {
      new_levels <- paste0("chr", levels_reads)
      new_levels[new_levels == "chrMT"] <- "chrM"
      names(new_levels) <- levels_reads
      normal_reads <- renameSeqlevels(normal_reads, new_levels)
    }
  } else if (opt$genomeStyle == "NCBI") {
    levels_reads <- seqlevels(normal_reads)
    if (any(grepl("^chr", levels_reads))) {
      new_levels <- gsub("^chr", "", levels_reads)
      new_levels[new_levels == "M"] <- "MT"
      names(new_levels) <- levels_reads
      normal_reads <- renameSeqlevels(normal_reads, new_levels)
    }
  }
  
  # Crucial safety check: ensure all elements in chrs are present in the GRanges' seqlevels
  seqlevels(normal_reads) <- unique(c(seqlevels(normal_reads), chrs))
  
  # Debug prints
  message("seqlevels(normal_reads): ", paste(seqlevels(normal_reads), collapse=", "))
  message("seqlevels(gc): ", paste(seqlevels(gc), collapse=", "))
  message("seqlevels(map): ", paste(seqlevels(map), collapse=", "))
  message("First 3 elements of normal_reads:")
  print(head(normal_reads, 3))
  message("First 3 elements of gc:")
  print(head(gc, 3))
  
  normal_counts <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, 
                                        centromere=centromere_df, targetedSequences=NULL,
                                        genomeStyle=opt$genomeStyle, chrNormalize=chrNormalize)
  
  if (is.null(normalGR)){
    normalGR <- normal_counts$counts
    values(normalGR) <- values(normalGR)$copy
    colnames(values(normalGR)) <- sid
  } else {
    values(normalGR)[[sid]] <- normal_counts$counts$copy
  }
  
  gender <- normal_counts$gender
  chrXMedian <- gender$chrXMedian
  chrXStr <- grep("X", chrs, value = TRUE)
  chrXInd <- as.character(seqnames(normalGR)) == chrXStr
  values(normalGR)[[sid]][chrXInd] <- values(normalGR)[[sid]][chrXInd] - chrXMedian
}

medianVal <- apply(values(normalGR), 1, median, na.rm = TRUE)
values(normalGR)[["Median"]] <- medianVal

# Save GR object directly to user-specified path
saveRDS(normalGR, file = opt$outfile)
message("Panel of Normals successfully written to: ", opt$outfile)
