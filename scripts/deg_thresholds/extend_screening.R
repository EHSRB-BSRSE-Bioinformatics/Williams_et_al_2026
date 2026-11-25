#!/usr/bin/env Rscript
# Extend the screening endpoint using stored BH values from the frozen trials.
source("scripts/deg_thresholds/analyze.R")
suppressPackageStartupMessages(library(DESeq2)) # S4 deserialization only

args <- commandArgs(TRUE)
if (length(args) > 3L) {
  stop("Usage: Rscript scripts/deg_thresholds/extend_screening.R [archive_root] [output_tsv] [cell]")
}
archive_root <- if (length(args)) args[[1L]] else "../2023_Meier_FDR"
output <- if (length(args) >= 2L) args[[2L]] else
  "data/derived/deg_thresholds/threshold_screening_p010.tsv"
frozen <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
index <- trial_index(names(frozen))
check_ids(index$trial, index$trial, "Frozen metadata")
manifest <- read.delim("data/reference/deg_source_archives.tsv")
archive_files <- c(TK6 = "tk6.RDS", `MCF-7` = "mcf7.RDS", HepaRG = "heparg.RDS",
                   `Liver spheroids` = "spheroids.RDS")
cells <- names(archive_files)
if (length(args) == 3L) {
  check(args[[3L]] %in% cells, "Unknown cell model")
  cells <- args[[3L]]
}
grid <- subset(threshold_grid(), padj_cutoff == .05)
grid$padj_cutoff <- .10
grid$threshold_id <- paste0("p0.1_fc", grid$fc_label)
fc_cutoffs <- log2(grid$linear_fc[-1L])

count_screening <- function(id, entry) {
  expected <- frozen[[id]]
  check(identical(as.character(entry$exp_metadata$sample_id), as.character(expected$sample_id)) &&
        identical(as.character(entry$exp_metadata$group), as.character(expected$group)),
        paste(id, "archive assignments differ"))
  results <- entry$DESeqRes
  p <- results$padj
  selected <- which(p < .10)
  abs_lfc <- abs(results$log2FoldChange[selected])
  c(length(selected), vapply(fc_cutoffs, function(fc) sum(abs_lfc > fc, na.rm = TRUE), integer(1)))
}

rates <- lapply(cells, function(cell) {
  source <- manifest[match(archive_files[[cell]], manifest$archive_file), ]
  archive_path <- file.path(archive_root, source$source_directory, source$archive_file)
  strata <- index[index$cell == cell, ]
  message("Reading ", cell, " archive")
  archive <- readRDS(archive_path)
  check_ids(names(archive), strata$trial, paste(cell, "archive"))
  counts <- vapply(strata$trial, function(id) count_screening(id, archive[[id]]), integer(nrow(grid)))
  rm(archive)
  gc()
  message("Verifying ", cell, " archive checksum")
  check(unname(sha256(archive_path)) == source$sha256, "Archive fingerprint differs")
  do.call(rbind, lapply(3:10, function(n) {
    selected <- strata$n == n
    values <- counts[, selected, drop = FALSE]
    check(ncol(values) == 1000L && setequal(strata$iteration[selected], 1:1000),
          "Incomplete frozen stratum")
    positive <- rowSums(values > 0L)
    cbind(data.frame(endpoint = "screening", cell = cell, n = n), grid,
          data.frame(positive_trials = positive, trials = ncol(values),
                     rate = positive / ncol(values), total_degs = rowSums(values),
                     median_degs = apply(values, 1, median), max_degs = apply(values, 1, max)))
  }))
})
write_tsv(do.call(rbind, rates), output)
