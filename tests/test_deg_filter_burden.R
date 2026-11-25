source("scripts/deg_thresholds/plot_filter_burden.R")
# Partial removal and loss of the last call must be counted separately.
z <- summarise_filter_burden(c(10, 2, 0, 4), c(5, 2, 0, 0), c(3, 0, 0, 0))
stopifnot(identical(unname(z$stages$total_degs), c(16, 7, 3)),
  identical(z$transitions$degs_removed, c(9, 4)),
  identical(z$transitions$reduced_still_positive, c(1L, 1L)),
  identical(z$transitions$reduced_to_zero, c(1L, 1L)),
  all(z$transitions$trials_with_reduction == 2L))
zero <- summarise_filter_burden(c(0, 0), c(0, 0), c(0, 0))
stopifnot(all(is.na(zero$transitions$fraction_degs_removed)))
for (bad in list(list(1, 2, 0), list(NA, 0, 0), list(1.5, 1, 0), list(c(1, 2), 1, 0))) {
  stopifnot(inherits(tryCatch(do.call(summarise_filter_burden, bad), error = identity), "error"))
}
cat("Paired count reductions, trial-state decomposition, zero denominators and invalid inputs verified.\n")

if (file.exists("data/derived/deg_thresholds/filter_burden.rds")) {
  actual <- readRDS("data/derived/deg_thresholds/filter_burden.rds")
  stopifnot(identical(unname(actual$pooled$stages$total_degs), c(925785, 187166, 185673)),
    identical(unname(actual$pooled$stages$positive_trials), c(5316, 2231, 2202)),
    actual$pooled$transitions$reduced_still_positive[2] == 308L,
    actual$pooled$transitions$reduced_to_zero[2] == 29L,
    all(actual$stages$trials == 1000L))
  for (stage in actual$pooled$stages$stage) {
    stopifnot(sum(actual$stages$total_degs[actual$stages$stage == stage]) ==
                actual$pooled$stages$total_degs[actual$pooled$stages$stage == stage])
  }
  cat("Saved pooled counts and condition-level totals reconcile.\n")
}
