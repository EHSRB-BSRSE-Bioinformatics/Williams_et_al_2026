# Summarise paired filtering stages from validated, frozen trial counts.
source("R/config.R")
source("R/figure_registry.R")
suppressPackageStartupMessages(source("R/summaries.R"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

summarise_filter_burden <- function(a, b, c) {
  if (length(a) != length(b) || length(a) != length(c)) {
    stop("Unpaired count lengths.")
  }
  counts <- cbind(a, b, c)
  if (
    !nrow(counts) ||
      anyNA(counts) ||
      any(!is.finite(counts)) ||
      any(counts < 0) ||
      any(counts != floor(counts)) ||
      any(b > a) ||
      any(c > b)
  ) {
    stop("Invalid or nonnested paired counts.")
  }
  total <- colSums(counts)
  stages <- data.frame(
    stage = c("BH p < 0.05", "+ FC > 1.5", "+ quantile/spike"),
    trials = nrow(counts),
    total_degs = total,
    mean_degs_per_trial = total / nrow(counts),
    positive_trials = colSums(counts > 0)
  )
  transitions <- do.call(
    rbind,
    lapply(1:2, function(i) {
      before <- counts[, i]
      after <- counts[, i + 1L]
      removed <- before - after
      data.frame(
        transition = c("Adding FC > 1.5", "Adding quantile/spike")[i],
        trials = length(before),
        degs_before = sum(before),
        degs_after = sum(after),
        degs_removed = sum(removed),
        fraction_degs_removed = if (sum(before)) {
          sum(removed) / sum(before)
        } else {
          NA_real_
        },
        positive_before = sum(before > 0),
        positive_after = sum(after > 0),
        trials_with_reduction = sum(removed > 0),
        reduced_still_positive = sum(removed > 0 & after > 0),
        reduced_to_zero = sum(before > 0 & after == 0),
        unchanged_positive = sum(before > 0 & removed == 0),
        unchanged_zero = sum(before == 0)
      )
    })
  )
  list(stages = stages, transitions = transitions)
}

run_filter_burden <- function() {
  x <- read_deg_figure_summary()
  labels <- levels(x$Filter)
  a <- x[x$Filter == labels[1], ]
  b <- x[x$Filter == labels[2], ]
  b <- b[match(a$name, b$name), ]
  c <- x[x$Filter == labels[3], ]
  c <- c[match(a$name, c$name), ]
  stopifnot(identical(a$name, b$name), identical(a$name, c$name))
  groups <- split(
    seq_len(nrow(a)),
    interaction(a$cell_line, a$n_per_group, drop = TRUE)
  )
  results <- lapply(groups, function(i) {
    z <- summarise_filter_burden(a$nDEGs[i], b$nDEGs[i], c$nDEGs[i])
    lapply(z, function(d) {
      cbind(cell = as.character(a$cell_line[i[1]]), n = a$n_per_group[i[1]], d)
    })
  })
  stages <- do.call(rbind, lapply(results, `[[`, "stages"))
  transitions <- do.call(rbind, lapply(results, `[[`, "transitions"))
  pooled <- summarise_filter_burden(a$nDEGs, b$nDEGs, c$nDEGs)
  stopifnot(
    nrow(stages) == 96L,
    nrow(transitions) == 64L,
    all(stages$trials == 1000L),
    all(transitions$trials == 1000L),
    all(
      transitions$trials_with_reduction ==
        transitions$reduced_still_positive + transitions$reduced_to_zero
    ),
    all(
      transitions$positive_before ==
        transitions$positive_after + transitions$reduced_to_zero
    )
  )
  out <- fdr_path("data", "derived", "deg_thresholds")
  tables <- list(
    filter_burden_stages = stages,
    filter_burden_transitions = transitions,
    filter_burden_pooled_stages = pooled$stages,
    filter_burden_pooled_transitions = pooled$transitions
  )
  for (nm in names(tables)) {
    write.table(
      tables[[nm]],
      file.path(out, paste0(nm, ".tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
  saveRDS(
    list(stages = stages, transitions = transitions, pooled = pooled),
    file.path(out, "filter_burden.rds"),
    compress = "xz"
  )
  stages$stage <- factor(stages$stage, levels = pooled$stages$stage)
  transitions$transition <- factor(
    transitions$transition,
    levels = pooled$transitions$transition
  )
  for (nm in c("stages", "transitions")) {
    d <- get(nm)
    d$cell <- factor(d$cell, levels = cell_line_levels)
    assign(nm, d)
  }
  top <- ggplot(stages, aes(cell, factor(n), fill = total_degs + 1)) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = format(total_degs, big.mark = ",", trim = TRUE)),
      size = 3
    ) +
    scale_fill_gradient(
      low = "#fff5eb",
      high = "#ef8a62",
      trans = "log10",
      name = "DEG calls + 1\n(log scale)"
    ) +
    facet_wrap(~stage, nrow = 1) +
    labs(
      title = "Total DEG calls across 1,000 frozen trials per condition",
      subtitle = "The same gene can contribute once in each trial; totals are not unique-gene counts.",
      x = NULL,
      y = "Sample size per group"
    )
  bottom <- ggplot(
    transitions,
    aes(cell, factor(n), fill = fraction_degs_removed)
  ) +
    geom_tile(colour = "white") +
    geom_text(
      aes(
        label = paste0(
          round(100 * fraction_degs_removed, 1),
          "%\n(",
          format(degs_removed, big.mark = ",", trim = TRUE),
          ")"
        )
      ),
      size = 3
    ) +
    scale_fill_gradient(
      low = "#fff5eb",
      high = "#d94801",
      limits = c(0, 1),
      labels = function(v) paste0(100 * v, "%"),
      name = "DEG calls\nremoved"
    ) +
    facet_wrap(~transition, nrow = 1) +
    labs(
      title = "Reduction from the preceding filtering stage",
      subtitle = "Tile labels: percent of preceding DEG calls removed (number removed).",
      x = NULL,
      y = "Sample size per group"
    )
  plot <- (top / bottom) +
    plot_layout(heights = c(1, 1.15)) &
    theme_minimal(base_size = 11) &
    theme(panel.grid = element_blank())
  for (ext in c("pdf", "png")) {
    path <- figure_path("deg_filter_burden_supplementary", ext)
    if (ext == "pdf") {
      ggsave(path, plot, width = 13, height = 10, device = cairo_pdf)
    } else {
      ggsave(path, plot, width = 13, height = 10, dpi = 300)
    }
  }
  paths <- c(
    fdr_path("data", "derived", "deg_thresholds", "threshold_sensitivity.rds"),
    fdr_path("scripts", "deg_thresholds", "plot_filter_burden.R"),
    fdr_path("R", "figure_registry.R"),
    fdr_path("R", "summaries.R"),
    fdr_path("tests", "test_deg_filter_burden.R"),
    file.path(out, paste0(names(tables), ".tsv")),
    file.path(out, "filter_burden.rds"),
    figure_path("deg_filter_burden_supplementary", c("pdf", "png"))
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
    file.path(out, "filter_burden_provenance.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  print(pooled)
}

if (sys.nframe() == 0L) {
  run_filter_burden()
}
