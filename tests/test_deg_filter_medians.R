source("scripts/deg_thresholds/plot_filter_medians.R")
z <- trial_count_statistics(c(0, 0, 2, 10))
stopifnot(z$median_all_trials == 1, z$median_positive_trials == 6,
          z$mean_degs == 3, z$trials == 4, z$positive_trials == 2, z$total_degs == 12)
z <- trial_count_statistics(c(0, 0, 0, 100))
stopifnot(z$median_all_trials == 0, z$median_positive_trials == 100)
z <- trial_count_statistics(c(0, 0))
stopifnot(z$median_all_trials == 0, is.na(z$median_positive_trials), z$positive_trials == 0)
for (bad in list(numeric(), NA_real_, -1, 1.5, Inf)) {
  stopifnot(inherits(tryCatch(trial_count_statistics(bad), error = identity), "error"))
}
x <- read_deg_figure_summary()
if (file.exists("data/derived/deg_thresholds/filter_medians.rds")) {
  results <- readRDS("data/derived/deg_thresholds/filter_medians.rds")
  stopifnot(nrow(results) == 96L, all(results$trials == 1000L), all(results$median_all_trials == 0),
    all(results$mean_degs == results$total_degs / results$trials))
  for (i in seq_len(nrow(results))) {
    d <- results[i, ]
    raw <- x$nDEGs[x$cell_line == d$cell & x$n_per_group == d$n & x$Filter == d$stage]
    expected <- trial_count_statistics(raw)
    stopifnot(isTRUE(all.equal(as.list(d[names(expected)]), as.list(expected), check.attributes = FALSE)))
  }
}
cat("Arithmetic means, zero-inclusive and conditional medians, empty-positive sets, invalid counts and saved stratum statistics verified.\n")
