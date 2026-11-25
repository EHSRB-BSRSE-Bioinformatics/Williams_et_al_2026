# Summarise saved ORA universes and apply alternative export cutoffs.
source("R/config.R")
suppressPackageStartupMessages(source("R/summaries.R"))
meta <- parse_trial_names(names(readRDS(derived_files$trial_metadata)))
frozen <- readRDS(derived_files$ora)
reactome <- readRDS(derived_files$reactome_ora)
stems <- c("TK6" = "tk6", "MCF-7" = "mcf7", "HepaRG" = "heparg", "Liver spheroids" = "spheroids")
count_export <- function(obj, p = attr(obj, "pvalueCutoff"), q = attr(obj, "qvalueCutoff")) {
  if (!isS4(obj)) return(0L)
  tab <- attr(obj, "result")
  take <- tab$pvalue <= p & tab$p.adjust <= p
  tab <- tab[take, , drop = FALSE]
  if (!anyNA(tab$qvalue)) tab <- tab[tab$qvalue <= q, , drop = FALSE]
  nrow(na.omit(tab))
}
universe_size <- function(obj) if (isS4(obj)) length(attr(obj, "universe")) else NA_integer_
results <- parallel::mclapply(names(stems), function(model) {
  ids <- meta[meta$cell_line == model, ]
  expected <- frozen[match(ids$run, frozen$run), ]
  x <- readRDS(file.path("../2023_Meier_FDR/bootstrap_ORA/2023-11-09",
                         paste0(stems[[model]], "_ora.RDS")))
  stopifnot(length(x) == nrow(ids), nrow(ids) == 8000L,
            is.null(names(x)) || identical(names(x), ids$run), identical(expected$run, ids$run))
  go <- lapply(x, `[[`, "enrich_go_bp")
  wp <- lapply(x, `[[`, "enricher_wikipathways")
  go_count <- vapply(go, count_export, integer(1))
  wp_count <- vapply(wp, count_export, integer(1))
  stopifnot(identical(unname(go_count), expected$go_bp),
            identical(unname(wp_count), expected$enricher_wikipathways))
  paired <- match(ids$run, reactome$run)
  stopifnot(!anyNA(paired), all(reactome$cell_line[paired] == model),
            all(reactome$n[paired] == ids$n_per_group))
  cat("Validated archived counts and read saved universes for", model, "\n")
  data.frame(cell_model = model, n_per_group = ids$n_per_group,
             go_universe = vapply(go, universe_size, integer(1)),
             wp_universe = vapply(wp, universe_size, integer(1)),
             reactome_universe = reactome$n_background_annotated[paired],
             go_positive = go_count > 0, wp_positive = wp_count > 0,
             wp_stricter_p = vapply(wp, count_export, integer(1), p = 0.01, q = 0.20) > 0,
             wp_go_cutoffs = vapply(wp, count_export, integer(1), p = 0.01, q = 0.05) > 0)
}, mc.cores = 4L)
stopifnot(all(vapply(results, is.data.frame, logical(1))))
raw <- bind_rows(results)
summary <- raw |>
  group_by(cell_model, n_per_group) |>
  group_modify(function(d, keys) {
    paired <- !is.na(d$go_universe) & !is.na(d$wp_universe)
    stopifnot(nrow(d) == 1000L, any(paired))
    data.frame(trials = nrow(d), trials_with_both_universes = sum(paired),
               go_median_universe_on_paired = median(d$go_universe[paired]),
               reactome_median_universe_on_paired = median(d$reactome_universe[paired]),
               wikipathways_median_universe_on_paired = median(d$wp_universe[paired]),
               go_eFPR_percent = 100 * mean(d$go_positive),
               wikipathways_eFPR_percent = 100 * mean(d$wp_positive),
               wikipathways_eFPR_percent_at_p_0.01 = 100 * mean(d$wp_stricter_p),
               wikipathways_eFPR_percent_at_GO_cutoffs = 100 * mean(d$wp_go_cutoffs))
  }) |>
  ungroup()
out <- "docs/revision/efpr_by_stratum/ora_universes_and_thresholds.tsv"
write.table(summary, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(summary[summary$n_per_group == 3L, ]), row.names = FALSE)
cat("Saved universe and cutoff diagnostics for the 32 matched strata.\n")
