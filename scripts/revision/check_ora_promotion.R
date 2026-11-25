source("R/config.R")
source("R/summaries.R")
ora <- readRDS(derived_files$ora)
trials <- readRDS(derived_files$reactome_ora)
joined <- attach_reactome_ora(ora, trials)
stopifnot(identical(as.data.frame(joined)[names(ora)], as.data.frame(ora)),
          identical(joined$run, ora$run), nrow(joined) == 32000L,
          setequal(joined$run, names(readRDS(derived_files$trial_metadata))))
reordered <- attach_reactome_ora(ora, trials[rev(seq_len(nrow(trials))), ])
stopifnot(identical(joined, reordered))
rates <- false_positive_rates(joined, "contains_reactome_false_positive", "n")
expected <- readRDS("docs/revision/reactome_ora_v2/summary.rds")
idx <- match(paste(rates$cell_line, rates$n), paste(expected$cell_line, expected$n))
stopifnot(!anyNA(idx), nrow(rates) == 32L,
          all(rates$n_trials_observed == 1000L),
          all(rates$false_positive_rate == expected$ora_rate[idx]),
          all(joined$reactome[joined$reactome_status == "insufficient_input"] == 0),
          sum(joined$contains_reactome_false_positive) == 557L)
expect_error <- function(value) stopifnot(inherits(tryCatch(value, error = identity), "error"))
bad <- trials; bad$run[2] <- bad$run[1]
expect_error(attach_reactome_ora(ora, bad))
expect_error(attach_reactome_ora(ora, trials[-1, ]))
bad <- trials; bad$status[1] <- "failed"
expect_error(attach_reactome_ora(ora, bad))
bad <- trials; bad$n_significant[1] <- NA_integer_
expect_error(attach_reactome_ora(ora, bad))
bad <- trials; bad$n[1] <- 99L
expect_error(attach_reactome_ora(ora, bad))
bad <- trials; i <- which(bad$status == "insufficient_input")[1]
bad$n_significant[i] <- 1L
expect_error(attach_reactome_ora(ora, bad))
cat("PASS: all historical columns/ranks retained; joins invariant to input order.\n",
    "PASS: exact frozen IDs, all 32 denominators, 557 Reactome-positive trials and zero-input handling.\n",
    "PASS: duplicate/missing IDs, failed trials, missing counts, wrong metadata and inconsistent statuses rejected.\n")
