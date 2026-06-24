library(dplyr)
library(ggplot2)
library(patchwork)

heatmap_base <- list(
  geom_tile(width = 1, height = 1, colour = NA, linewidth = 0),
  scale_y_reverse(
    limits = c(n_trials_default + 0.5, 0.5),
    expand = c(0, 0)
  ),
  scale_x_continuous(
    breaks = sample_size_levels,
    limits = c(min(sample_size_levels) - 0.5, max(sample_size_levels) + 0.5),
    expand = c(0, 0)
  ),
  theme_minimal(),
  xlab("Sample size per group"),
  ylab("Monte Carlo trial rank"),
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
)

rate_base <- list(
  geom_point(),
  geom_line(),
  geom_hline(yintercept = 0.05, col = "red", lty = "dashed"),
  theme_bw(),
  xlab("Sample size per group"),
  ylab("Empirical false-positive rate"),
  scale_x_continuous(breaks = sample_size_levels),
  scale_y_continuous(limits = c(0, 1)),
  theme(
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank()
  )
)

count_fill_scale <- function(title) {
  scale_fill_viridis_c(
    option = "F",
    name = title,
    direction = -1,
    trans = "log10",
    breaks = c(1, 11, 101, 1001, 10001),
    labels = c("0", "10", "100", "1,000", "10,000")
  )
}

build_deg_figure <- function(deg_summary) {
  deg_summary <- deg_summary %>%
    mutate(
      Filter = gsub("R-ODAF \\+ P-Value", "P-value filter only", Filter),
      Filter = factor(Filter, levels = c("P-value filter only", "Final DEGs"))
    ) %>%
    filter(!is.na(Filter))

  deg_scale <- count_fill_scale("Number\nof DEGs\nin trial")

  heatmap_unfiltered <- ggplot(
    filter(deg_summary, Filter == "P-value filter only"),
    aes(x = n_per_group, y = rank, fill = nDEGs + 1)
  ) +
    heatmap_base +
    deg_scale +
    facet_grid(~cell_line) +
    ggtitle("DEGs filtered on adjusted p-value")

  heatmap_filtered <- ggplot(
    filter(deg_summary, Filter == "Final DEGs"),
    aes(x = n_per_group, y = rank, fill = nDEGs + 1)
  ) +
    heatmap_base +
    deg_scale +
    facet_grid(~cell_line) +
    theme(legend.position = "none") +
    ggtitle("DEGs filtered on adjusted p-value and fold change")

  rates <- false_positive_rates(
    deg_summary,
    flag_col = "contains_false_positive",
    n_col = "n_per_group",
    group_cols = c("cell_line", "Filter")
  )

  rate_unfiltered <- ggplot(
    filter(rates, Filter == "P-value filter only"),
    aes(x = n_per_group, y = false_positive_rate, group = cell_line)
  ) +
    rate_base +
    facet_grid(~cell_line)

  rate_filtered <- ggplot(
    filter(rates, Filter == "Final DEGs"),
    aes(x = n_per_group, y = false_positive_rate, group = cell_line)
  ) +
    rate_base +
    facet_grid(~cell_line)

  (heatmap_unfiltered | heatmap_filtered) /
    (rate_unfiltered | rate_filtered) +
    plot_layout(guides = "collect", heights = c(1, 0.75)) +
    plot_annotation(tag_levels = list(c("A", "B", "", "")))
}

build_ora_figure <- function(ora_summary) {
  ora_scale <- count_fill_scale("Enriched\npathways\nin trial")

  go_heatmap <- ggplot(
    ora_summary,
    aes(x = n, y = go_bp_rank, fill = go_bp + 1)
  ) +
    heatmap_base +
    ora_scale +
    facet_grid(~cell_line) +
    theme(legend.position = "none") +
    ggtitle("GO over-representation analysis")

  wiki_heatmap <- ggplot(
    ora_summary,
    aes(x = n, y = wiki_rank, fill = enricher_wikipathways + 1)
  ) +
    heatmap_base +
    ora_scale +
    facet_grid(~cell_line) +
    ggtitle("WikiPathways over-representation analysis")

  go_rates <- false_positive_rates(
    ora_summary,
    flag_col = "contains_GO_false_positive",
    n_col = "n"
  )
  wiki_rates <- false_positive_rates(
    ora_summary,
    flag_col = "contains_wiki_false_positive",
    n_col = "n"
  )

  go_rate_plot <- ggplot(
    go_rates,
    aes(x = n, y = false_positive_rate, group = cell_line)
  ) +
    rate_base +
    facet_grid(~cell_line)
  wiki_rate_plot <- ggplot(
    wiki_rates,
    aes(x = n, y = false_positive_rate, group = cell_line)
  ) +
    rate_base +
    facet_grid(~cell_line)

  (go_heatmap | wiki_heatmap) /
    (go_rate_plot | wiki_rate_plot) +
    plot_layout(guides = "collect", heights = c(1, 0.75)) +
    plot_annotation(tag_levels = list(c("A", "B", "", "")))
}

build_gsea_figure <- function(gsea_summary,
                              min_size = gsea_min_size_default,
                              max_size = gsea_max_size_default) {
  if (!all(c("min_size", "max_size") %in% names(gsea_summary))) {
    stop(
      "The GSEA summary does not record min_size and max_size. Rebuild it ",
      "from trial-level results before rendering the publication figure."
    )
  }

  gsea_summary <- gsea_summary %>%
    filter(.data$min_size == min_size, .data$max_size == max_size)
  if (nrow(gsea_summary) == 0L) {
    stop(
      "The GSEA summary has no results for gene-set sizes ", min_size,
      "–", max_size, "."
    )
  }

  gsea_scale <- count_fill_scale("Enriched\npathways\nin trial")

  pval_heatmap <- ggplot(
    gsea_summary,
    aes(x = n, y = rank, fill = p_lt_0.05 + 1)
  ) +
    heatmap_base +
    gsea_scale +
    facet_grid(~cell_line) +
    theme(legend.position = "none") +
    ggtitle("GSEA, adjusted p-value < 0.05")

  nes_heatmap <- ggplot(
    gsea_summary,
    aes(x = n, y = rank_filtered, fill = p_lt_0.05_ES_gt_2 + 1)
  ) +
    heatmap_base +
    gsea_scale +
    facet_grid(~cell_line) +
    ggtitle("GSEA, adjusted p-value < 0.05 and absolute NES > 2")

  pval_rates <- false_positive_rates(
    gsea_summary,
    flag_col = "contains_false_positive_pval_only",
    n_col = "n"
  )
  nes_rates <- false_positive_rates(
    gsea_summary,
    flag_col = "contains_false_positive_NES_pval_filtered",
    n_col = "n"
  )

  pval_rate_plot <- ggplot(
    pval_rates,
    aes(x = n, y = false_positive_rate, group = cell_line)
  ) +
    rate_base +
    facet_grid(~cell_line)
  nes_rate_plot <- ggplot(
    nes_rates,
    aes(x = n, y = false_positive_rate, group = cell_line)
  ) +
    rate_base +
    facet_grid(~cell_line)

  (pval_heatmap | nes_heatmap) /
    (pval_rate_plot | nes_rate_plot) +
    plot_layout(guides = "collect", heights = c(1, 0.75)) +
    plot_annotation(
      tag_levels = list(c("A", "B", "", "")),
      subtitle = paste0(
        "Gene-set size: ", min_size, "–", max_size, " genes"
      )
    )
}

build_tgx_ddi_figure <- function(tgx_summary) {
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("The ggridges package is required for TGx-DDI ridge plots.")
  }

  tgx_summary <- tgx_summary %>%
    mutate(n_per_group = factor(as.integer(as.character(n_per_group))))

  pc1_plot <- ggplot(
    filter(tgx_summary, tgx_ddi_test == "pc1_res", tgx_ddi_value >= 0),
    aes(x = tgx_ddi_value, y = n_per_group, group = n_per_group)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    ggridges::geom_density_ridges(
      alpha = 1,
      point_alpha = 0.01,
      jittered_points = TRUE
    ) +
    geom_point(
      data = filter(tgx_summary, tgx_ddi_test == "pc1_res", tgx_ddi_value < 0),
      color = "red",
      size = 0.8
    ) +
    facet_grid(~cell_line) +
    xlab("PC1 value") +
    ylab("Sample size per group") +
    theme_minimal() +
    ggtitle("TGx-DDI first principal component distribution")

  prob_plot <- ggplot(
    filter(tgx_summary, tgx_ddi_test == "probG_res", tgx_ddi_value >= 0),
    aes(x = tgx_ddi_value, y = n_per_group, group = n_per_group)
  ) +
    geom_vline(xintercept = 0.9, linetype = "dashed", color = "blue") +
    ggridges::geom_density_ridges(
      alpha = 1,
      point_alpha = 0.01,
      jittered_points = TRUE
    ) +
    geom_point(
      data = filter(
        tgx_summary,
        tgx_ddi_test == "probG_res",
        tgx_ddi_value >= 0.9
      ),
      color = "red",
      size = 0.8
    ) +
    facet_grid(~cell_line) +
    scale_x_log10() +
    xlab("Probability analysis score") +
    ylab("Sample size per group") +
    theme_minimal() +
    ggtitle("TGx-DDI probability analysis distribution")

  pc1_plot / prob_plot + plot_annotation(tag_levels = "A")
}
