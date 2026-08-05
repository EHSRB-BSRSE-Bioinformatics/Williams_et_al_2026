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

tpod_design_labels <- c(
  d4n6 = "4d x 6n",
  d4n12 = "4d x 12n",
  d4n18 = "4d x 18n",
  d6n4 = "6d x 4n",
  d6n8 = "6d x 8n",
  d6n12 = "6d x 12n",
  d8n3 = "8d x 3n",
  d8n6 = "8d x 6n",
  d8n9 = "8d x 9n",
  d12n2 = "12d x 2n",
  d12n4 = "12d x 4n",
  d12n6 = "12d x 6n"
)

tpod_compact_design_labels <- c(
  d4n6 = "4 × 6",
  d4n12 = "4 × 12",
  d4n18 = "4 × 18",
  d6n4 = "6 × 4",
  d6n8 = "6 × 8",
  d6n12 = "6 × 12",
  d8n3 = "8 × 3",
  d8n6 = "8 × 6",
  d8n9 = "8 × 9",
  d12n2 = "12 × 2",
  d12n4 = "12 × 4",
  d12n6 = "12 × 6"
)

tpod_cell_line_levels <- c("TK6", "MCF-7", "HepaRG", "Liver Spheroids")
tpod_analysis_levels <- c(
  "Gene (rank 25)",
  "Pathway-KEGG",
  "Pathway-Reactome",
  "Pathway-Hallmark",
  "Pathway-GO"
)

build_tpod_fdr_panel <- function(tpod_summary,
                                 cell_line,
                                 parameter_levels,
                                 show_x_labels = FALSE) {
  panel_data <- tpod_summary %>%
    filter(CellLine == cell_line) %>%
    mutate(
      Design = factor(
        Design,
        levels = names(tpod_design_labels),
        labels = unname(tpod_design_labels)
      ),
      tPoD = factor(tPoD, levels = tpod_analysis_levels),
      Param = factor(Param, levels = rev(parameter_levels))
    ) %>%
    droplevels()

  ggplot(panel_data, aes(x = Design, y = Param, fill = FPR)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(FPR > 0, sprintf("%.2f", FPR), "")),
      size = 3,
      color = "black",
      fontface = "bold"
    ) +
    facet_wrap(~tPoD, nrow = 1) +
    scale_fill_gradient(
      low = "#FFFFFF",
      high = "#CC0000",
      limits = c(0, 1),
      name = "FDR"
    ) +
    labs(
      title = cell_line,
      x = if (show_x_labels) "Experimental design" else NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      axis.text.x = if (show_x_labels) {
        element_text(angle = 45, hjust = 1, color = "black")
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_labels) element_line() else element_blank(),
      axis.text.y = element_text(color = "black"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

build_tpod_fdr_figure <- function(tpod_summary) {
  required_columns <- c(
    "R_Squared_Filter",
    "Fold_Change_Cutoff",
    "P_Value_Cutoff",
    "tPoD",
    "FPR",
    "Design",
    "CellLine"
  )
  missing_columns <- setdiff(required_columns, names(tpod_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "The tPoD FDR summary is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  tpod_summary <- tpod_summary %>%
    mutate(
      Param = paste0(
        ifelse(R_Squared_Filter == "R-Sq > 0.6", "R2>=0.6", "No R2"),
        " | FC", Fold_Change_Cutoff,
        " | p", P_Value_Cutoff
      )
    )

  parameter_levels <- tpod_summary %>%
    arrange(
      P_Value_Cutoff,
      Fold_Change_Cutoff,
      desc(R_Squared_Filter == "R-Sq > 0.6")
    ) %>%
    pull(Param) %>%
    unique()

  panels <- Map(
    function(cell_line, show_x_labels) {
      build_tpod_fdr_panel(
        tpod_summary = tpod_summary,
        cell_line = cell_line,
        parameter_levels = parameter_levels,
        show_x_labels = show_x_labels
      )
    },
    tpod_cell_line_levels,
    tpod_cell_line_levels == tail(tpod_cell_line_levels, 1L)
  )

  wrap_plots(panels, ncol = 1, heights = c(1, 1, 1, 1.3), guides = "collect") +
    plot_annotation(
      title = "False discovery rate by experimental design and analysis parameters",
      subtitle = "Values are shown inside tiles where FDR > 0",
      caption = paste(
        "d = dose groups; n = animals per dose group |",
        "FC = fold-change cutoff; p = p-value cutoff"
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", hjust = 0.5),
        plot.caption = element_text(color = "grey50", hjust = 0.5)
      )
    ) &
    theme(legend.position = "right")
}

build_tpod_robustness_figure <- function(tpod_summary, fdr_limit = 0.05) {
  required_columns <- c("tPoD", "FPR", "Design", "CellLine")
  missing_columns <- setdiff(required_columns, names(tpod_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "The tPoD FDR summary is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (length(fdr_limit) != 1L || !is.finite(fdr_limit) ||
      fdr_limit < 0 || fdr_limit > 1) {
    stop("fdr_limit must be one finite value between 0 and 1.")
  }

  robustness_data <- tpod_summary %>%
    mutate(
      Design = factor(
        Design,
        levels = names(tpod_compact_design_labels),
        labels = unname(tpod_compact_design_labels)
      ),
      CellLine = factor(CellLine, levels = tpod_cell_line_levels),
      tPoD = factor(tPoD, levels = rev(tpod_analysis_levels))
    ) %>%
    group_by(CellLine, tPoD, Design) %>%
    summarise(
      n_settings = sum(!is.na(FPR)),
      n_exceed = sum(FPR > fdr_limit, na.rm = TRUE),
      max_fdr = if (all(is.na(FPR))) NA_real_ else max(FPR, na.rm = TRUE),
      .groups = "drop"
    )

  expected_settings <- max(robustness_data$n_settings, na.rm = TRUE)
  if (!is.finite(expected_settings) || expected_settings < 1L) {
    stop("The tPoD FDR summary contains no non-missing FDR values.")
  }

  ggplot(robustness_data, aes(x = Design, y = tPoD)) +
    geom_point(
      data = filter(robustness_data, n_exceed > 0L),
      aes(size = n_exceed, fill = max_fdr),
      shape = 21,
      colour = "grey20",
      stroke = 0.3
    ) +
    geom_vline(
      xintercept = c(3.5, 6.5, 9.5),
      colour = "grey75",
      linewidth = 0.3
    ) +
    facet_wrap(~CellLine, ncol = 1, strip.position = "left") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_size_area(
      limits = c(0, expected_settings),
      breaks = unique(c(1, 2, 4, 6, expected_settings)),
      max_size = 8,
      name = paste0(
        "Parameter settings\nwith FDR > ",
        format(fdr_limit, trim = TRUE)
      )
    ) +
    scale_fill_gradient(
      low = "#FEE8C8",
      high = "#B30000",
      limits = c(fdr_limit, 1),
      name = "Maximum FDR"
    ) +
    labs(
      x = "Experimental design: dose groups × replicates per group",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold",
        hjust = 1,
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

build_tpod_parameter_setting_figure <- function(tpod_summary) {
  required_columns <- c(
    "R_Squared_Filter",
    "Fold_Change_Cutoff",
    "P_Value_Cutoff",
    "tPoD",
    "FPR",
    "Design",
    "CellLine"
  )
  missing_columns <- setdiff(required_columns, names(tpod_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "The tPoD FDR summary is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  setting_data <- tpod_summary %>%
    mutate(
      Param = paste0(
        ifelse(R_Squared_Filter == "R-Sq > 0.6", "R² ≥ 0.6", "No R²"),
        " | FC ", Fold_Change_Cutoff,
        " | p ", P_Value_Cutoff
      )
    )

  parameter_levels <- setting_data %>%
    arrange(
      P_Value_Cutoff,
      Fold_Change_Cutoff,
      desc(R_Squared_Filter == "R-Sq > 0.6")
    ) %>%
    pull(Param) %>%
    unique()

  setting_data <- setting_data %>%
    mutate(
      Design = factor(
        Design,
        levels = names(tpod_compact_design_labels),
        labels = unname(tpod_compact_design_labels)
      ),
      CellLine = factor(CellLine, levels = tpod_cell_line_levels),
      tPoD = factor(tPoD, levels = rev(tpod_analysis_levels)),
      Param = factor(Param, levels = parameter_levels),
      design_x = as.integer(Design) + case_when(
        P_Value_Cutoff == 0.01 & Fold_Change_Cutoff == 1.5 ~ -0.27,
        P_Value_Cutoff == 0.01 & Fold_Change_Cutoff == 2.0 ~ -0.09,
        P_Value_Cutoff == 0.05 & Fold_Change_Cutoff == 1.5 ~ 0.09,
        P_Value_Cutoff == 0.05 & Fold_Change_Cutoff == 2.0 ~ 0.27,
        TRUE ~ 0
      ),
      tpod_y = as.integer(tPoD) + if_else(
        R_Squared_Filter == "R-Sq > 0.6",
        0.13,
        -0.13
      )
    )

  setting_colours <- c(
    "#0072B2", "#56B4E9", "#009E73", "#7FC97F",
    "#D55E00", "#E69F00", "#CC79A7", "#000000"
  )
  names(setting_colours) <- parameter_levels

  ggplot(
    setting_data,
    aes(x = design_x, y = tpod_y, size = FPR, fill = Param)
  ) +
    geom_point(
      data = filter(setting_data, FPR > 0),
      shape = 21,
      colour = "grey20",
      stroke = 0.2,
      alpha = 0.9,
      show.legend = TRUE
    ) +
    geom_vline(
      xintercept = c(3.5, 6.5, 9.5),
      colour = "grey65",
      linewidth = 0.35
    ) +
    facet_wrap(~CellLine, ncol = 1, strip.position = "left") +
    scale_x_continuous(
      breaks = seq_along(tpod_compact_design_labels),
      labels = unname(tpod_compact_design_labels),
      limits = c(0.5, length(tpod_compact_design_labels) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(tpod_analysis_levels),
      labels = rev(tpod_analysis_levels),
      limits = c(0.5, length(tpod_analysis_levels) + 0.5),
      expand = c(0, 0)
    ) +
    scale_size_area(
      limits = c(0, 1),
      breaks = c(0.05, 0.10, 0.25, 0.50, 1),
      max_size = 3.5,
      name = "FDR"
    ) +
    scale_fill_manual(
      values = setting_colours,
      drop = FALSE,
      name = "Analysis parameter setting"
    ) +
    labs(
      x = "Experimental design: dose groups × replicates per group",
      y = NULL
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(size = 4, alpha = 1),
        ncol = 1
      ),
      size = guide_legend(
        order = 2,
        override.aes = list(fill = "grey60", alpha = 1)
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold",
        hjust = 1,
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

build_tpod_fdr_supplementary_figure <- function(tpod_summary) {
  required_columns <- c(
    "R_Squared_Filter",
    "Fold_Change_Cutoff",
    "P_Value_Cutoff",
    "tPoD",
    "FPR",
    "Design",
    "CellLine"
  )
  missing_columns <- setdiff(required_columns, names(tpod_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "The tPoD FDR summary is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  supplementary_data <- tpod_summary %>%
    mutate(
      Param = paste0(
        ifelse(R_Squared_Filter == "R-Sq > 0.6", "R2>=0.6", "No R2"),
        " | FC", Fold_Change_Cutoff,
        " | p", P_Value_Cutoff
      ),
      Design = factor(
        Design,
        levels = names(tpod_compact_design_labels),
        labels = unname(tpod_compact_design_labels)
      ),
      CellLine = factor(CellLine, levels = tpod_cell_line_levels),
      tPoD = factor(tPoD, levels = tpod_analysis_levels),
      FDR_bin = cut(
        FPR,
        breaks = c(-Inf, 0, 0.05, 0.10, 0.25, Inf),
        labels = c("0", ">0–0.05", ">0.05–0.10", ">0.10–0.25", ">0.25")
      )
    )

  parameter_levels <- supplementary_data %>%
    arrange(
      P_Value_Cutoff,
      Fold_Change_Cutoff,
      desc(R_Squared_Filter == "R-Sq > 0.6")
    ) %>%
    pull(Param) %>%
    unique()

  supplementary_data <- supplementary_data %>%
    mutate(Param = factor(Param, levels = rev(parameter_levels)))

  ggplot(supplementary_data, aes(x = Design, y = Param, fill = FDR_bin)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(CellLine ~ tPoD) +
    scale_fill_manual(
      values = c(
        "0" = "#F2F2F2",
        ">0–0.05" = "#FEE8C8",
        ">0.05–0.10" = "#FDBB84",
        ">0.10–0.25" = "#E34A33",
        ">0.25" = "#B30000"
      ),
      drop = FALSE,
      name = "FDR"
    ) +
    labs(
      x = "Experimental design: dose groups × replicates per group",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}
