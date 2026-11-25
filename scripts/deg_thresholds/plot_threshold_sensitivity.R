#!/usr/bin/env Rscript
args <- if (sys.nframe() == 0L) {
  commandArgs(trailingOnly = TRUE)
} else {
  character()
}
if (length(args) && !identical(args, "--only=heatmap") && !identical(args, "--include-drafts")) {
  stop(
    "Usage: Rscript scripts/deg_thresholds/plot_threshold_sensitivity.R [--only=heatmap | --include-drafts]",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(ggplot2))
source("R/config.R")

padj_cutoffs <- c(.01, .05, .10)
padj_labels <- sprintf("%.2f", padj_cutoffs)
fc_labels <- c("none", "1.25", "1.5", "2.0")
fc_plot_labels <- c("None", fc_labels[-1L])
data_dir <- fdr_path("data", "derived", "deg_thresholds")
archive <- read.delim(file.path(data_dir, "threshold_rates.tsv"))
extension <- read.delim(file.path(data_dir, "threshold_screening_p010.tsv"))
screening_rates <- rbind(
  subset(archive, endpoint == "screening" & padj_cutoff %in% padj_cutoffs),
  extension
)
screening <- screening_rates
stopifnot(
  all(extension$padj_cutoff == .10),
  all(screening$endpoint == "screening"),
  all(screening$trials == n_trials_default),
  all(screening$rate == screening$positive_trials / screening$trials)
)

screening$cell <- factor(screening$cell, levels = cell_line_levels)
screening$n <- factor(screening$n, levels = sample_size_levels)
screening$fc <- factor(screening$fc_label, levels = fc_labels)
screening$padj <- factor(
  screening$padj_cutoff,
  levels = padj_cutoffs,
  labels = padj_labels
)
coverage <- with(screening, table(cell, n, fc, padj))
stopifnot(sum(coverage) == nrow(screening), all(coverage == 1L))
rates <- with(screening, tapply(rate, list(cell, n, fc, padj), identity))
stopifnot(
  all(rates[,,, 1L] <= rates[,,, 2L]),
  all(rates[,,, 2L] <= rates[,,, 3L]),
  all(rates[,, 1:3, ] >= rates[,, 2:4, ])
)

rate_max <- ceiling(10 * max(screening$rate)) / 10
published <- subset(screening, padj_cutoff == .05 & linear_fc == 1.5)
rate_label <- "Empirical false-positive rate"
fc_label <- "Absolute linear fold-change cutoff"
cutoff_label <- "BH-adjusted p-value cutoff"
published_label <- "the nominal cutoffs used in the main study (adjusted p < 0.05 and fold change > 1.5)."
outline <- "#000000"

panels <- ggplot(screening, aes(x = fc)) +
  facet_grid(
    rows = vars(n),
    cols = vars(cell),
    switch = "y",
    labeller = labeller(n = function(n) paste("n =", n))
  ) +
  labs(
    title = "Impact of filtering thresholds on DEG false positive rates",
    subtitle = "1,000 mock trials per cell model and sample size; numbers show percentages of experiments with ≥1 DEG."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(10, "pt"),
    panel.spacing.y = grid::unit(6, "pt"),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", colour = "black"),
    strip.text.y.left = element_text(angle = 0),
    axis.text = element_text(colour = "grey20"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(10, 12, 10, 10)
  )

fold_change_axis <- list(
  scale_x_discrete(labels = fc_plot_labels, expand = expansion(add = 0.5)),
  labs(x = fc_label)
)

# Main DEG, ORA and GSEA heatmaps all use reversed Rocket (viridis option F).
heatmap <- panels +
  fold_change_axis +
  aes(y = padj) +
  geom_tile(aes(fill = rate), colour = "white", linewidth = 0.35) +
  geom_tile(
    data = published,
    width = 0.96,
    height = 0.91,
    fill = NA,
    colour = "white",
    linewidth = 0.8
  ) +
  geom_tile(
    data = published,
    width = 0.96,
    height = 0.91,
    fill = NA,
    colour = outline,
    linewidth = 0.42
  ) +
  geom_text(
    aes(label = sprintf("%.1f", 100 * rate), colour = rate >= rate_max * 0.42),
    size = 3.1
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "black", `TRUE` = "white"),
    guide = "none"
  ) +
  scale_fill_viridis_c(
    option = "F",
    direction = -1,
    limits = c(0, rate_max),
    breaks = seq.int(0, round(10 * rate_max)) / 10,
    labels = scales::label_percent(),
    name = rate_label,
    guide = guide_colourbar(
      title.position = "top",
      barwidth = grid::unit(2, "in")
    )
  ) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(
    y = cutoff_label,
    caption = paste(
      "Black outlines mark",
      published_label
    )
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title.x.bottom = element_text(
      margin = margin(t = 12)
    )
  )


rate_panels <- panels +
  aes(y = rate) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  geom_point(
    data = published,
    shape = 1,
    size = 3,
    stroke = 0.5,
    colour = outline,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, rate_max),
    breaks = seq(0, rate_max, length.out = 3L),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(y = rate_label, caption = paste("Blue circles mark", published_label)) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    panel.spacing.y = grid::unit(9, "pt"),
    legend.key.width = grid::unit(24, "pt")
  )

line_plot <- rate_panels +
  aes(group = padj, colour = padj, linetype = padj, shape = padj) +
  fold_change_axis +
  scale_colour_manual(
    name = cutoff_label,
    values = viridisLite::viridis(9, option = "F")[c(1, 3, 6)]
  ) +
  scale_linetype_manual(
    name = cutoff_label,
    values = c("dotted", "solid", "longdash")
  ) +
  scale_shape_manual(name = cutoff_label, values = c(16, 17, 15))

swapped_line_plot <- line_plot +
  coord_flip() +
  theme(
    panel.grid.major.x = element_line(),
    panel.grid.major.y = element_blank()
  )

padj_axis_line_plot <- rate_panels +
  aes(x = padj, group = fc, colour = fc, linetype = fc, shape = fc) +
  scale_x_discrete(
    name = cutoff_label,
    limits = rev(padj_labels),
    expand = expansion(add = 0.5)
  ) +
  scale_colour_manual(
    name = fc_label,
    labels = fc_plot_labels,
    values = viridisLite::viridis(13, option = "F")[c(1, 4, 7, 9)]
  ) +
  scale_linetype_manual(
    name = fc_label,
    labels = fc_plot_labels,
    values = c("dotted", "dotdash", "solid", "longdash")
  ) +
  scale_shape_manual(
    name = fc_label,
    labels = fc_plot_labels,
    values = c(16, 17, 15, 18)
  )

render_first_draft <- function(path, data) {
  colours <- hcl.colors(101, "YlOrRd", rev = TRUE)
  pdf(path, width = 13, height = 11)
  on.exit(dev.off())
  for (endpoint in c("screening", "full_rodaf")) {
    layout(matrix(1:32, nrow = 8, ncol = 4, byrow = TRUE))
    par(oma = c(4, 4.4, 5, .5), mar = c(.75, 1, .5, .35))
    for (n in sample_size_levels) {
      for (cell in cell_line_levels) {
        d <- data[data$endpoint == endpoint & data$cell == cell & data$n == n, ]
        plot.new()
        plot.window(
          xlim = c(.5, 4.5),
          ylim = c(.5, 3.5),
          xaxs = "i",
          yaxs = "i"
        )
        if (endpoint == "full_rodaf") {
          rect(.5, 2.5, 4.5, 3.5, col = "grey94", border = "white")
          text(2.5, 3, "Unavailable", cex = .75, col = "grey35")
        }
        x <- match(d$fc_label, fc_labels)
        y <- match(d$padj_cutoff, padj_cutoffs)
        proportion <- d$rate / rate_max
        rect(
          x - .5,
          y - .5,
          x + .5,
          y + .5,
          col = colours[1 + round(100 * proportion)],
          border = "white"
        )
        text(
          x,
          y,
          sprintf("%.1f", 100 * d$rate),
          cex = .85,
          col = ifelse(proportion > .6, "white", "#222222")
        )
        marked <- which(d$padj_cutoff == .05 & d$linear_fc == 1.5)
        rect(
          x[marked] - .47,
          y[marked] - .47,
          x[marked] + .47,
          y[marked] + .47,
          border = outline,
          lwd = 1.5
        )
        if (n == min(sample_size_levels)) {
          mtext(cell, side = 3, line = .15, cex = .9)
        }
        if (cell == cell_line_levels[1L]) {
          axis(
            2,
            at = 1:3,
            labels = padj_labels,
            las = 1,
            tick = FALSE,
            cex.axis = .65,
            line = -.5
          )
          mtext(paste0("n=", n), side = 2, line = 1.8, cex = .65)
        }
        if (n == max(sample_size_levels)) {
          axis(
            1,
            at = 1:4,
            labels = fc_labels,
            tick = FALSE,
            cex.axis = .75,
            line = -.5
          )
        }
      }
    }
    title <- if (endpoint == "screening") {
      "A. Stored BH p + fold-change screening"
    } else {
      "B. Stored BH p + fold change + quantile/spike rules"
    }
    subtitle <- if (endpoint == "screening") {
      sprintf(
        "Trials with at least one DEG (%% of 1,000); colour scale 0 to %.0f%%",
        100 * rate_max
      )
    } else {
      "Archived quantile and spike rule outcomes cover adjusted p < 0.05."
    }
    mtext(title, side = 3, outer = TRUE, line = 3.2, cex = 1.1)
    mtext(subtitle, side = 3, outer = TRUE, line = 1.6, cex = .8)
    mtext(fc_label, side = 1, outer = TRUE, line = 1.3, cex = .85)
    mtext(
      paste("Blue outlines show", published_label),
      side = 1,
      outer = TRUE,
      line = 2.55,
      cex = .75
    )
    mtext(cutoff_label, side = 2, outer = TRUE, line = 2.3, cex = .85)
  }
}

save_figure <- function(plot, stem, dpi) {
  ggsave(
    paste0(stem, ".pdf"),
    plot,
    width = 11,
    height = 10.5,
    device = cairo_pdf
  )
  ggsave(
    paste0(stem, ".png"),
    plot,
    width = 11,
    height = 10.5,
    dpi = dpi,
    bg = "white"
  )
}

source("R/figure_registry.R")
save_figure(
  heatmap,
  figure_path("s4_threshold_heatmap"),
  dpi = 600
)
if (identical(args, "--include-drafts")) {
  save_figure(
    line_plot,
    figure_path("supplementary_threshold_lines_draft"),
    dpi = 300
  )
  save_figure(
    swapped_line_plot,
    figure_path("supplementary_threshold_lines_swapped_draft"),
    dpi = 300
  )
  save_figure(
    padj_axis_line_plot,
    figure_path("supplementary_threshold_lines_padj_axis_draft"),
    dpi = 300
  )
  full_rodaf <- subset(
    archive,
    endpoint == "full_rodaf" & padj_cutoff %in% padj_cutoffs
  )
  render_first_draft(
    figure_path("supplementary_threshold_heatmap_first_draft", "pdf"),
    rbind(screening_rates, full_rodaf)
  )
}
