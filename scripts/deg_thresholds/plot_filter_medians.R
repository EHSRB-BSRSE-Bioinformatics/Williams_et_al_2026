# Describe per-trial DEG counts, retaining zero-call trials explicitly.
source("R/config.R")
source("R/figure_registry.R")
suppressPackageStartupMessages(source("R/summaries.R"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

trial_count_statistics <- function(counts) {
  if (
    !length(counts) ||
      anyNA(counts) ||
      any(!is.finite(counts)) ||
      any(counts < 0) ||
      any(counts != floor(counts))
  ) {
    stop("Invalid trial counts.")
  }
  positive <- counts[counts > 0]
  data.frame(
    trials = length(counts),
    positive_trials = length(positive),
    total_degs = sum(counts),
    mean_degs = mean(counts),
    median_all_trials = median(counts),
    q25_all_trials = unname(quantile(counts, .25)),
    q75_all_trials = unname(quantile(counts, .75)),
    median_positive_trials = if (length(positive)) {
      median(positive)
    } else {
      NA_real_
    },
    q25_positive_trials = if (length(positive)) {
      unname(quantile(positive, .25))
    } else {
      NA_real_
    },
    q75_positive_trials = if (length(positive)) {
      unname(quantile(positive, .75))
    } else {
      NA_real_
    },
    max_degs = max(counts)
  )
}

run_filter_medians <- function() {
  x <- read_deg_figure_summary()
  groups <- split(
    seq_len(nrow(x)),
    interaction(x$cell_line, x$n_per_group, x$Filter, drop = TRUE)
  )
  stats <- do.call(
    rbind,
    lapply(groups, function(i) {
      cbind(
        cell = as.character(x$cell_line[i[1]]),
        n = x$n_per_group[i[1]],
        stage = as.character(x$Filter[i[1]]),
        trial_count_statistics(x$nDEGs[i])
      )
    })
  )
  stopifnot(nrow(stats) == 96L, all(stats$trials == 1000L))
  out <- fdr_path("data", "derived", "deg_thresholds")
  write.table(
    stats,
    file.path(out, "filter_medians.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  saveRDS(stats, file.path(out, "filter_medians.rds"), compress = "xz")
  stats$cell <- factor(stats$cell, levels = cell_line_levels)
  stats$stage <- factor(
    stats$stage,
    levels = levels(x$Filter),
    labels = c("BH p < 0.05", "+ FC > 1.5", "+ quantile/spike")
  )
  shared_limits <- c(
    1,
    max(c(stats$mean_degs, stats$median_positive_trials), na.rm = TRUE) + 1
  )
  make_panel <- function(metric, title, subtitle, conditional = FALSE) {
    d <- stats
    d$value <- d[[metric]]
    d$label <- ifelse(is.na(d$value), "NA", format(d$value, trim = TRUE))
    if (conditional) {
      d$label <- paste0(d$label, "\n(n=", d$positive_trials, ")")
    }
    ggplot(d, aes(cell, factor(n), fill = value + 1)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = label), size = 3) +
      scale_fill_gradient(
        low = "#fff5eb",
        high = "#ef8a62",
        trans = "log10",
        limits = shared_limits,
        na.value = "grey85",
        breaks = c(1, 2, 6, 11, 21, 51, 101),
        labels = c(0, 1, 5, 10, 20, 50, 100),
        name = "DEGs per trial\n(log scale)"
      ) +
      facet_wrap(~stage, nrow = 1) +
      labs(
        title = title,
        subtitle = subtitle,
        x = NULL,
        y = "Sample size per group"
      )
  }
  top <- make_panel(
    "mean_degs",
    "Mean DEGs per trial: all frozen trials",
    "Arithmetic mean = total calls / 1,000 trials, including zero-call trials."
  )
  bottom <- make_panel(
    "median_positive_trials",
    "Median DEGs among trials with at least one DEG",
    "Tile labels: median (n = positive trials). Conditioning changes the trial subset at each stage.",
    TRUE
  )
  plot <- (top / bottom) +
    plot_layout(guides = "collect") &
    theme_minimal(base_size = 11) &
    theme(panel.grid = element_blank())
  figures <- figure_path("deg_filter_medians_supplementary", c("pdf", "png"))
  ggsave(figures[1], plot, width = 13, height = 10, device = cairo_pdf)
  ggsave(figures[2], plot, width = 13, height = 10, dpi = 300)
  paths <- c(
    fdr_path("data", "derived", "deg_thresholds", "threshold_sensitivity.rds"),
    fdr_path("scripts", "deg_thresholds", "plot_filter_medians.R"),
    fdr_path("R", "figure_registry.R"),
    fdr_path("tests", "test_deg_filter_medians.R"),
    fdr_path("R", "summaries.R"),
    file.path(out, c("filter_medians.tsv", "filter_medians.rds")),
    figures
  )
  hashes <- vapply(
    paths,
    function(p) {
      strsplit(system2("sha256sum", shQuote(p), stdout = TRUE), " +")[[1]][1]
    },
    ""
  )
  write.table(
    data.frame(path = paths, sha256 = hashes),
    file.path(out, "filter_medians_provenance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  print(stats[stats$cell == "TK6" & stats$n == 10, ])
  cat(
    "Nonzero all-trial medians:",
    sum(stats$median_all_trials > 0),
    "of",
    nrow(stats),
    "conditions/stages.\n"
  )
}
if (sys.nframe() == 0L) {
  run_filter_medians()
}
