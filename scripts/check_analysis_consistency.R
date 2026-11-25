# Deterministic revision regression checks; no model fits or archive writes.
script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}
suppressPackageStartupMessages(source("R/config.R"))
suppressPackageStartupMessages(source("R/summaries.R"))
suppressPackageStartupMessages(source("R/figures.R"))

# Rate is P(at least one call), including groups with no/all positive trials.
flags <- data.frame(cell_line = rep(c("none", "some", "all"), each = 4),
                    n = 3L, positive = c(rep(FALSE, 4), TRUE, FALSE, TRUE, FALSE, rep(TRUE, 4)))
rates <- false_positive_rates(flags, "positive", "n", expected_trials = 4L)
stopifnot(all(rates$n_trials_observed == 4L),
          identical(rates$false_positive_rate[match(c("none", "some", "all"), rates$cell_line)], c(0, .5, 1)))
warned <- FALSE
rates <- withCallingHandlers(false_positive_rates(flags, "positive", "n", expected_trials = 1000L),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
stopifnot(warned, all(rates$n_trials_observed == 4L), all(rates$expected_trials == 1000L),
          rates$false_positive_rate[rates$cell_line == "some"] == .5)

# Do not count the NA-index placeholder rows of the historical data.frame code.
res <- data.frame(padj = c(.01, .05, NA, .01), log2FoldChange = c(1, 1, 1, log2(1.5)),
                  row.names = c("pass", "p_boundary", "missing", "fc_boundary"))
trial <- list(DESeqRes = res, DESeqResDEGs = res[FALSE, ], DESeqResFiltered = res)
summary <- summarise_degs_from_trials(list(TK6_3_3_1 = trial))
stopifnot(summary$nDEGs[summary$Filter == "P-value filter only"] == 2,
          summary$nDEGs[summary$Filter == "P-value + fold change"] == 1,
          summary$nDEGs[summary$Filter == "Final DEGs"] == 0,
          nrow(res[res$padj < .05, ]) == 3)

# The archive records excluded significant IDs; they reconstruct final calls
# independently of the recycled intermediate. Check every bundled reference.
reference <- readRDS("data/reference/deg_reference.rds")
for (x in reference) {
  r <- x$DESeqRes
  sig <- r[!is.na(r$padj) & r$padj < .05, , drop = FALSE]
  excluded <- unique(unlist(lapply(x$DESeqjs_filtered_out[2:3], rownames)))
  pass <- !rownames(sig) %in% excluded
  final <- if (nrow(sig)) sig[pass & abs(sig$log2FoldChange) > log2(1.5), , drop = FALSE] else data.frame()
  recycled <- if (nrow(sig)) r[pass, , drop = FALSE] else data.frame()
  stopifnot(identical(final, x$DESeqResDEGs), identical(recycled, x$DESeqResFiltered))
}

# A real archive example: 3 significant rows become a one-row legacy A count,
# while final calls correctly contain 2 features.
x <- reference[["Liver spheroids_3_3_334"]]
stopifnot(sum(x$DESeqRes$padj < .05, na.rm = TRUE) == 3,
          nrow(x$DESeqRes[x$DESeqResFiltered$padj < .05, ]) == 1,
          nrow(x$DESeqResDEGs) == 2)

# Figure 2 plots the two screening stages. Validate the archived final calls separately.
frozen <- factor_cell_lines(as.data.frame(readRDS(derived_files$deg)))
corrected <- read_deg_figure_summary()
stopifnot(inherits(tryCatch(build_deg_figure(frozen), error = identity), "error"))
plot <- build_deg_figure(corrected)
heatmap_a <- plot[[1]][[1]]$data
heatmap_b <- plot[[1]][[2]]$data
final <- filter(corrected, Filter == "Final DEGs")
stopifnot(all(as.character(heatmap_a$Filter) == "P-value filter only"),
          all(as.character(heatmap_b$Filter) == "P-value + fold-change"),
          sum(final$nDEGs) == sum(frozen$nDEGs[frozen$Filter %in% "Final DEGs"]))

# Every plotted condition has 1000 unique frozen assignments, no missing flags.
metadata <- readRDS(derived_files$trial_metadata)
stopifnot(length(metadata) == 32000L, !anyDuplicated(names(metadata)))
for (nm in names(reference)) {
  stopifnot(identical(as.character(reference[[nm]]$exp_metadata$sample_id),
                      as.character(metadata[[nm]]$sample_id)),
            identical(as.character(reference[[nm]]$exp_metadata$group),
                      as.character(metadata[[nm]]$group)))
}
for (dat in list(heatmap_a, heatmap_b, final, as.data.frame(readRDS(derived_files$ora)),
                 as.data.frame(readRDS(derived_files$gsea)))) {
  id <- if ("run" %in% names(dat)) dat$run else dat$name
  stopifnot(setequal(id, names(metadata)), !anyDuplicated(id))
  n <- if ("n_per_group" %in% names(dat)) dat$n_per_group else dat$n
  stopifnot(all(table(dat$cell_line, n) == 1000L))
  flags <- names(dat)[startsWith(names(dat), "contains_")]
  stopifnot(!anyNA(dat[flags]))
}

# Size filtering changes the BH family: use pval, not archived full-family padj.
gsea <- data.frame(size = c(14L, 15L, 500L, 501L), pval = c(.00001, .02, .8, .00001),
                   padj = 0, NES = c(3, 2, -3, 3))
g <- summarise_gsea_from_trials(list(TK6_3_3_1 = gsea))
stopifnot(g$num_rows == 2L, g$p_lt_0.05 == 1L, g$p_lt_0.05_ES_gt_2 == 0L)

# Optional pathway dependencies: show why pre-selection Entrez collapsing is an
# analytical change, even though repeated annotation rows alone are harmless.
if (requireNamespace("data.table", quietly = TRUE)) {
  source("R/deg_trials.R")
  source("R/pathway_trials.R")
  r <- data.frame(log2FoldChange = c(1, 2, -3, 4, -1, 5, NA),
                  padj = c(.01, .2, .01, .4, .01, .01, .01),
                  row.names = c("ALPHA_p1", "ALPHA_p2", "BETA_p1", "GAMMA_p1",
                                "DELTA_p1", "UNMAPPED_p1", "EPSILON_p1"))
  annotation <- data.frame(ALIAS = c("ALPHA", "BETA", "GAMMA", "DELTA", "DELTA", "DELTA"),
                           ENTREZID = c("1", "2", "2", "3", "4", "3"))
  collapsed <- collapse_deseq_results_to_entrez(list(DESeqRes = r), annotation)
  calls <- collapsed$ENTREZID[collapsed$padj < .05 & abs(collapsed$linearFoldChange) > 1.5]
  stopifnot(setequal(collapsed$ENTREZID, as.character(1:4)), !anyDuplicated(collapsed$ENTREZID),
            collapsed$log2FoldChange[collapsed$ENTREZID == "1"] == 2,
            setequal(calls, c("3", "4")))
  # Original ORA lines 50-70: choose a probe per alias, then threshold the
  # expanded annotation join without choosing a representative per Entrez ID.
  historical <- data.table::data.table(r, keep.rownames = "rn")
  historical <- historical[, .(log2FoldChange, rn, padj)]
  historical[, rn := sub("_(.*)", "", rn)]
  historical <- historical[, .SD[which.max(abs(log2FoldChange))], by = rn]
  historical <- dplyr::left_join(as.data.frame(historical), annotation, by = c("rn" = "ALIAS"), relationship = "many-to-many")
  old_calls <- with(historical, unique(ENTREZID[!is.na(ENTREZID) & padj < .05 & abs(log2FoldChange) > log2(1.5)]))
  stopifnot(setequal(old_calls, c("2", "3", "4")))
  message("Pathway alias/probe/Entrez collapse checks passed (synthetic counterexample).")
} else {
  message("Optional pathway collapse check skipped: data.table unavailable.")
}
message("Analysis consistency regression checks passed (17 archived DEG references; all plotted trial IDs).")
