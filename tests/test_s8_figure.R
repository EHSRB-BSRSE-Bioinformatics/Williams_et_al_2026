source("R/config.R")
suppressPackageStartupMessages(source("R/summaries.R"))
suppressPackageStartupMessages(source("R/figures.R"))

summary <- read_s8_figure_summary()
data <- summary$data

stopifnot(
    nrow(data) > 0L,
    setequal(unique(as.character(data$cell_line)), cell_line_levels),
    all(is.finite(data$mean_log2cpm)),
    all(is.finite(data$sd_log2cpm)),
    all(data$sd_log2cpm >= 0),
    all(data$count[!is.na(data$count)] >= 1L),
    all(data$count[!is.na(data$count)] <= n_trials_default),
    all(data$n[!is.na(data$n)] %in% sample_size_levels),
    all(is.na(data$n) == is.na(data$count)),
    any(is.na(data$n)),
    !anyDuplicated(data[c("cell_line", "feature", "n")])
)

for (cell_line in cell_line_levels) {
    cell_data <- data[data$cell_line == cell_line, , drop = FALSE]
    stopifnot(all(
        sample_size_levels %in% unique(cell_data$n[!is.na(cell_data$n)])
    ))
    stopifnot(any(is.na(cell_data$n)))
}

plot <- build_s8_mean_variance_figure(summary)
smoother <- plot$layers[[2L]]
stopifnot(
    inherits(smoother$stat, "StatSmooth"),
    identical(smoother$stat_params$method, "loess"),
    identical(smoother$stat_params$span, 0.75),
    identical(smoother$stat_params$se, FALSE),
    isFALSE(smoother$inherit.aes)
)

cat(
    "Figure S8 summary and descriptive per-panel LOESS specification verified.\n"
)
