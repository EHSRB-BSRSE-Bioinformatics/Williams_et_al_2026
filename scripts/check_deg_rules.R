# Small boundary checks for the rules most easily changed during cleanup.
script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/deg_trials.R")
groups <- factor(rep(c("A", "B"), each = 4))
cpm <- rbind(
  three_at_threshold = c(0.5, 0.5, 0.5, 0, 0, 0, 0, 0),
  only_two = c(0.5, 0.5, 0, 0, 0, 0, 0, 0),
  other_group = c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0),
  below_threshold = rep(0.499, 8)
)
stopifnot(identical(unname(rodaf_relevance_filter(cpm, groups, 0.5)),
                    c(TRUE, FALSE, TRUE, FALSE)))
groups <- factor(rep(c("A", "B"), each = 3))
counts <- rbind(
  separated = c(10, 10, 10, 1, 1, 1),
  tied = rep(10, 6),
  one_spike = c(0, 0, 100, 1, 1, 1),
  two_spikes = c(0, 0, 100, 0, 0, 100),
  zero_group = c(0, 0, 0, 0, 0, 100)
)
filters <- rodaf_post_test_filters(counts, groups)
stopifnot(identical(as.integer(filters[, "quantile"]), c(1L, 0L, 0L, 0L, 0L)),
          identical(as.integer(filters[, "spike"]), c(1L, 1L, 1L, 0L, 1L)))
message("DEG rule boundary checks passed.")
