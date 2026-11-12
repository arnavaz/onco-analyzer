#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ichorCNA))
suppressPackageStartupMessages(library(GenomicRanges))

sample_id <- "DYNA-OTB-715257-FFPE-T-DNA-U_S214"
base_dir <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/Dyna_v2/DNASeq/results/DYNA-OTB-715257"
wig_file <- file.path(base_dir, "bin_size_benchmark", paste0(sample_id, "_1000000.wig"))
pon_file <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_panelOfNormals.rds"
excl_file <- "/cluster/projects/mohccn/Sandbox/JEFF_TMP/MultiMMR_HMF/panel/v2/ichor_resources/custom_exclusion.bed"

cat("=== FUNCTION SIGNATURES ===\n")
cat("loadReadCountsFromWig arguments:\n")
print(args(loadReadCountsFromWig))
cat("\nnormalizeByPanelOrMatchedNormal arguments:\n")
print(args(normalizeByPanelOrMatchedNormal))

cat("\n=== normalizeByPanelOrMatchedNormal BODY ===\n")
print(body(normalizeByPanelOrMatchedNormal))

cat("\n=== 1. WIG FILE CHECK ===\n")
cat("WIG file:", wig_file, "\n")
cat("WIG file exists:", file.exists(wig_file), "\n")

cat("\n=== 2. PON CHECK ===\n")
cat("PoN file:", pon_file, "\n")
cat("PoN file exists:", file.exists(pon_file), "\n")
if (file.exists(pon_file)) {
  pon <- readRDS(pon_file)
  cat("PoN class:", class(pon), "\n")
  if (inherits(pon, "GRanges")) {
    cat("PoN seqlevels:", head(seqlevels(pon)), "\n")
    cat("PoN total bins:", length(pon), "\n")
  }
}

cat("\n=== 3. EXCLUSION BED CHECK ===\n")
cat("Exclusion BED:", excl_file, "\n")
cat("Exclusion BED exists:", file.exists(excl_file), "\n")
excl_df <- NULL
if (file.exists(excl_file)) {
  excl_raw <- read.delim(excl_file, header=FALSE)
  if (ncol(excl_raw) >= 3) {
    colnames(excl_raw)[1:3] <- c("chr", "start", "end")
  }
  excl_df <- excl_raw
}

cat("\n=== 4. CALL loadReadCountsFromWig USING FORMAL ARGS ===\n")
gcWig <- system.file('extdata', 'gc_hg38_1000kb.wig', package='ichorCNA')
mapWig <- system.file('extdata', 'map_hg38_1000kb.wig', package='ichorCNA')
centromere <- system.file('extdata', 'GRCh38.GCA_000001405.2_centromere_acen.txt', package='ichorCNA')

chrs <- paste0("chr", c(1:22, "X"))
wig_df <- data.frame(patientID=sample_id, tumour_file=wig_file, stringsAsFactors=FALSE)

# Extract formal argument names dynamically
formal_args <- names(formals(loadReadCountsFromWig))
cat("Formal args of loadReadCountsFromWig:", paste(formal_args, collapse=", "), "\n")

# Construct call dynamically based on formal argument names
arg_list <- list()
if ("wigFiles" %in% formal_args) arg_list[["wigFiles"]] <- wig_df
if ("countsfile" %in% formal_args) arg_list[["countsfile"]] <- wig_df
if ("chrs" %in% formal_args) arg_list[["chrs"]] <- chrs
if ("gcWig" %in% formal_args) arg_list[["gcWig"]] <- gcWig
if ("gc" %in% formal_args) arg_list[["gc"]] <- gcWig
if ("mapWig" %in% formal_args) arg_list[["mapWig"]] <- mapWig
if ("map" %in% formal_args) arg_list[["map"]] <- mapWig
if ("centromere" %in% formal_args) arg_list[["centromere"]] <- centromere
if ("targetedSequences" %in% formal_args) arg_list[["targetedSequences"]] <- excl_df
if ("genomeStyle" %in% formal_args) arg_list[["genomeStyle"]] <- "UCSC"
if ("genomeBuild" %in% formal_args) arg_list[["genomeBuild"]] <- "hg38"

counts <- tryCatch({
  do.call(loadReadCountsFromWig, arg_list)
}, error = function(e) {
  cat("do.call loadReadCountsFromWig failed:", e$message, "\n")
  return(NULL)
})

if (!is.null(counts)) {
  tumour_copy <- counts[[1]]
  cat("Total bins in tumour_copy:", length(tumour_copy), "\n")
  cat("Valid bins count:", sum(tumour_copy$valid, na.rm=TRUE), "\n")
  cat("NA copy count:", sum(is.na(tumour_copy$copy)), "\n")

  cat("\n=== 5. CHECK PON NORMALIZATION ===\n")
  if (file.exists(pon_file)) {
    norm_res <- tryCatch({
      normalizeByPanelOrMatchedNormal(tumour_copy, chrs = chrs, normal_panel = pon_file, gender = "male")
    }, error = function(e) {
      cat("ERROR in normalizeByPanelOrMatchedNormal:", e$message, "\n")
      return(NULL)
    })
    if (!is.null(norm_res)) {
      cat("After PoN norm - Total bins:", length(norm_res), "\n")
      cat("After PoN norm - Valid bins:", sum(norm_res$valid, na.rm=TRUE), "\n")
      cat("After PoN norm - NA copy count:", sum(is.na(norm_res$copy)), "\n")
      cat("Summary of copy values after PoN norm:\n")
      print(summary(norm_res$copy))
    }
  }

  cat("\n=== 6. CHECK PARAMETERS & RUNEM CALL ===\n")
  chrTrain <- chrs[chrs != "chrX"]
  chrInd <- as.character(seqnames(tumour_copy)) %in% chrTrain
  valid <- tumour_copy$valid
  
  cat("chrInd sum:", sum(chrInd), "\n")
  cat("valid & chrInd sum:", sum(valid & chrInd, na.rm=TRUE), "\n")
  
  logR <- as.data.frame(tumour_copy$copy)
  cat("logR rows:", nrow(logR), "\n")
  cat("logR[valid & chrInd, , drop=F] rows:", nrow(logR[valid & chrInd, , drop=F]), "\n")
  
  param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = 5, includeHOMD = FALSE, ct.sc = NULL, ploidy = 2)
  cat("getDefaultParameters output structure:\n")
  str(param)
}
