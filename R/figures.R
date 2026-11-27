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
  y_scale = scale_y_continuous(limits = c(0, 1)),
  theme(
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank()
  )
)

count_fill_scale <- function(title, limits = NULL) {
  scale_fill_viridis_c(
    option = "F",
    name = title,
    direction = -1,
    trans = "log10",
    limits = limits,
    breaks = c(1, 11, 101, 1001, 10001),
    labels = c("0", "10", "100", "1,000", "10,000")
  )
}

build_deg_figure <- function(deg_summary) {
  # Reject legacy summaries rather than relabeling recycled intermediate counts.
  sources <- c(
    "archived_DESeqRes_padj",
    "archived_DESeqRes_padj_lfc",
    "verified_archived_final"
  )
  labels <- c("P-value filter only", "P-value + fold-change", "Final DEGs")
  if (
    !"count_source" %in% names(deg_summary) ||
      !setequal(unique(as.character(deg_summary$Filter)), labels) ||
      anyNA(deg_summary$count_source) ||
      any(
        deg_summary$count_source !=
          sources[match(as.character(deg_summary$Filter), labels)]
      )
  ) {
    stop(
      "Use read_deg_figure_summary() for validated canonical DEG plotting data."
    )
  }
  plot_labels <- labels[1:2]
  deg_summary <- deg_summary %>%
    filter(Filter %in% plot_labels) %>%
    mutate(Filter = factor(Filter, levels = plot_labels))

  deg_scale <- count_fill_scale(
    "Number\nof DEGs\nin trial",
    limits = c(1, max(1, deg_summary$nDEGs) + 1)
  )
  titles <- c("DEGs: BH-adjusted p < 0.05", "DEGs: BH p < 0.05 and FC > 1.5")
  heatmaps <- lapply(seq_along(plot_labels), function(i) {
    ggplot(
      filter(deg_summary, Filter == plot_labels[i]),
      aes(x = n_per_group, y = rank, fill = nDEGs + 1)
    ) +
      heatmap_base +
      deg_scale +
      facet_grid(~cell_line) +
      ggtitle(titles[i])
  })
  rates <- false_positive_rates(
    deg_summary,
    flag_col = "contains_false_positive",
    n_col = "n_per_group",
    group_cols = c("cell_line", "Filter")
  )
  rate_plots <- lapply(plot_labels, function(label) {
    ggplot(
      filter(rates, Filter == label),
      aes(x = n_per_group, y = false_positive_rate, group = cell_line)
    ) +
      rate_base +
      facet_grid(~cell_line)
  })

  (heatmaps[[1]] | heatmaps[[2]]) /
    (rate_plots[[1]] | rate_plots[[2]]) +
    plot_layout(guides = "collect", heights = c(1, 0.75)) +
    plot_annotation(tag_levels = list(c("A", "B", "", "")))
}

build_ora_figure <- function(
  ora_summary,
  collection = c("Reactome", "WikiPathways"),
  include_go = TRUE
) {
  collection <- match.arg(collection, several.ok = TRUE)
  panels <- data.frame(
    row.names = c("GO", "Reactome", "WikiPathways"),
    count = c("go_bp", "reactome", "enricher_wikipathways"),
    rank = c("go_bp_rank", "reactome_rank", "wiki_rank"),
    flag = c(
      "contains_GO_false_positive",
      "contains_reactome_false_positive",
      "contains_wiki_false_positive"
    )
  )
  panels <- panels[c(if (include_go) "GO", unique(collection)), , drop = FALSE]
  required <- unlist(panels, use.names = FALSE)
  missing <- setdiff(required, names(ora_summary))
  if (length(missing)) {
    stop(
      "Missing ORA columns: ",
      paste(missing, collapse = ", "),
      ". Load inputs with read_summary_inputs()."
    )
  }

  count_max <- max(1, unlist(ora_summary[panels$count], use.names = FALSE))
  ora_scale <- count_fill_scale(
    "Enriched\npathways\nin trial",
    limits = c(1, count_max + 1)
  )
  heatmaps <- lapply(seq_len(nrow(panels)), function(i) {
    count <- panels$count[i]
    trial_rank <- panels$rank[i]
    ggplot(
      ora_summary,
      aes(x = n, y = .data[[trial_rank]], fill = .data[[count]] + 1)
    ) +
      heatmap_base +
      ora_scale +
      facet_grid(~cell_line) +
      ggtitle(paste(rownames(panels)[i], "over-representation analysis")) +
      theme(legend.position = if (i == nrow(panels)) "right" else "none")
  })

  rates <- lapply(panels$flag, function(flag) {
    false_positive_rates(ora_summary, flag, "n")
  })
  rate_style <- rate_base
  if (nrow(panels) == 3L) {
    rate_style$y_scale <- scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent(accuracy = 1)
    )
  }
  rate_plots <- lapply(rates, function(data) {
    ggplot(data, aes(x = n, y = false_positive_rate, group = cell_line)) +
      rate_style +
      facet_grid(~cell_line)
  })

  wrap_plots(heatmaps, nrow = 1) /
    wrap_plots(rate_plots, nrow = 1) +
    plot_layout(guides = "collect", heights = c(1, 0.75)) +
    plot_annotation(
      tag_levels = list(c(
        LETTERS[seq_len(nrow(panels))],
        rep("", nrow(panels))
      ))
    )
}

build_gsea_figure <- function(
  gsea_summary,
  min_size = gsea_min_size_default,
  max_size = gsea_max_size_default
) {
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
      "The GSEA summary has no results for gene-set sizes ",
      min_size,
      "–",
      max_size,
      "."
    )
  }

  gsea_scale <- count_fill_scale(
    "Enriched\npathways\nin trial",
    limits = c(
      1,
      max(1, gsea_summary$p_lt_0.05, gsea_summary$p_lt_0.05_ES_gt_2) + 1
    )
  )

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
        "Gene-set size: ",
        min_size,
        "–",
        max_size,
        " genes"
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

build_tpod_fdr_panel <- function(
  tpod_summary,
  cell_line,
  parameter_levels,
  show_x_labels = FALSE
) {
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
        " | FC",
        Fold_Change_Cutoff,
        " | p",
        P_Value_Cutoff
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
        "d = dose groups; n = samples per dose group |",
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

tpod_matrix_key <- function() {
  grid::grobTree(
    grid::textGrob(
      c("1.5", "2.0"),
      x = c(0.56, 0.79),
      y = 0.90,
      gp = grid::gpar(fontsize = 9)
    ),
    grid::textGrob(
      c("No R²", "R² > 0.6"),
      x = 0.32,
      y = c(0.64, 0.29),
      just = "right",
      gp = grid::gpar(fontsize = 9)
    ),
    grid::rectGrob(
      x = rep(c(0.56, 0.79), 2),
      y = rep(c(0.64, 0.29), each = 2),
      width = grid::unit(0.62, "cm"),
      height = grid::unit(0.62, "cm"),
      gp = grid::gpar(fill = "grey92", col = "grey50", lwd = 0.6)
    )
  )
}

build_tpod_fdr_main_figure <- function(
  tpod_summary,
  p_main = 0.05,
  b_summary = c("mean", "max"),
  n_mock = 100L
) {
  b_summary <- match.arg(b_summary)
  required <- c(
    "R_Squared_Filter",
    "Fold_Change_Cutoff",
    "P_Value_Cutoff",
    "tPoD",
    "FPR",
    "Design",
    "CellLine"
  )
  missing <- setdiff(required, names(tpod_summary))
  if (length(missing)) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  # Normalise potentially factor/character cutoff columns.
  dat <- tpod_summary %>%
    dplyr::mutate(
      P_Value_Cutoff = as.numeric(as.character(P_Value_Cutoff)),
      Fold_Change_Cutoff = as.numeric(as.character(Fold_Change_Cutoff)),
      FPR = as.numeric(as.character(FPR)),
      R_Squared_Filter = as.character(R_Squared_Filter),
      Design = as.character(Design),
      tPoD = as.character(tPoD),
      CellLine = as.character(CellLine)
    )

  expected_fc <- c(1.5, 2)
  expected_r2 <- c("No R-Sq", "R-Sq > 0.6")
  if (
    anyNA(dat[, setdiff(required, "FPR")]) ||
      any(!is.na(dat$FPR) & !is.finite(dat$FPR)) ||
      any(dat$FPR < 0 | dat$FPR > 1, na.rm = TRUE) ||
      !setequal(unique(dat$Fold_Change_Cutoff), expected_fc) ||
      !setequal(unique(dat$R_Squared_Filter), expected_r2) ||
      !setequal(unique(dat$Design), names(tpod_compact_design_labels)) ||
      !setequal(unique(dat$tPoD), tpod_analysis_levels) ||
      !setequal(unique(dat$CellLine), tpod_cell_line_levels)
  ) {
    stop(
      "Invalid or unexpected tPoD design, endpoint, cell line, filter, or FPR value."
    )
  }
  if (
    anyDuplicated(dat[, c(
      "CellLine",
      "Design",
      "tPoD",
      "R_Squared_Filter",
      "Fold_Change_Cutoff",
      "P_Value_Cutoff"
    )])
  ) {
    stop(
      "Expected one FPR per cell line, design, endpoint and parameter setting."
    )
  }
  if (!all(c(0.01, 0.05) %in% dat$P_Value_Cutoff)) {
    stop("Panel A requires both p-value cutoffs (0.01 and 0.05).")
  }

  # A: descriptive averages across the same endpoint/design/FC/R-squared grid.
  # No uncertainty bars: settings re-use mock experiments and are correlated.
  overview <- dat %>%
    dplyr::filter(P_Value_Cutoff %in% c(0.01, 0.05)) %>%
    dplyr::group_by(CellLine, P_Value_Cutoff) %>%
    dplyr::summarise(
      mean_fpr = if (all(is.na(FPR))) NA_real_ else mean(FPR, na.rm = TRUE),
      n_nonzero = sum(FPR > 0, na.rm = TRUE),
      n_settings = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      CellLine = factor(CellLine, levels = tpod_cell_line_levels),
      p_label = factor(
        sprintf("p = %.2f", P_Value_Cutoff),
        levels = c("p = 0.01", "p = 0.05")
      ),
      text_label = ifelse(
        is.na(mean_fpr),
        "Missing",
        ifelse(
          n_nonzero == 0L,
          "None\nobserved",
          scales::percent(mean_fpr, accuracy = 0.1)
        )
      ),
      label_y = ifelse(n_nonzero == 0L & !is.na(mean_fpr), 0.004, mean_fpr)
    )

  panel_a <- ggplot2::ggplot(
    overview,
    ggplot2::aes(x = CellLine, y = mean_fpr)
  ) +
    ggplot2::geom_col(width = 0.66, fill = "#C44E46") +
    ggplot2::geom_text(
      ggplot2::aes(y = label_y, label = text_label),
      vjust = -0.45,
      size = 3.3
    ) +
    ggplot2::facet_wrap(~p_label, nrow = 1) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::labs(
      title = "A   Effect of the Williams trend-test threshold",
      subtitle = "Mean observed false-positive rate across designs, endpoints and other filtering settings",
      x = NULL,
      y = "Mean observed FPR"
    ) +
    ggplot2::theme_bw(base_size = 11.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      strip.text = ggplot2::element_text(face = "bold", size = 11.5),
      strip.background = ggplot2::element_rect(fill = "grey93", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = "black", size = 11)
    )

  # B: show each FC x R-squared setting at the chosen p cutoff, without averaging.
  setting_data <- dat %>%
    dplyr::filter(P_Value_Cutoff == p_main)
  if (!nrow(setting_data)) {
    stop("No data at p_main = ", p_main)
  }

  # Complete the full grid explicitly: absent records remain NA, never zero.
  matrix_grid <- expand.grid(
    CellLine = tpod_cell_line_levels,
    Design = names(tpod_compact_design_labels),
    tPoD = tpod_analysis_levels,
    Fold_Change_Cutoff = expected_fc,
    R_Squared_Filter = expected_r2,
    stringsAsFactors = FALSE
  )
  key_cols <- names(matrix_grid)
  key <- function(x) do.call(paste, c(x[key_cols], sep = "\r"))
  matched <- match(key(matrix_grid), key(setting_data))
  if (anyNA(matched) || anyNA(setting_data$FPR)) {
    warning(
      sum(is.na(matched)),
      " missing Panel B combinations and ",
      sum(is.na(setting_data$FPR)),
      " missing FPR values; shown in grey."
    )
  }
  detail <- matrix_grid
  detail$FPR <- setting_data$FPR[matched]

  detail <- detail %>%
    dplyr::mutate(
      CellLine = factor(CellLine, levels = tpod_cell_line_levels),
      design_x = match(Design, names(tpod_compact_design_labels)),
      endpoint_y = match(tPoD, rev(tpod_analysis_levels)),
      tile_x = design_x + ifelse(Fold_Change_Cutoff == 1.5, -0.23, 0.23),
      tile_y = endpoint_y + ifelse(R_Squared_Filter == "No R-Sq", 0.23, -0.23),
      fpr_label = ifelse(is.na(FPR), "NA", sprintf("%.2f", FPR)),
      label_colour = ifelse(is.na(FPR) | FPR >= 0.25, "white", "grey15")
    )

  panel_b <- ggplot2::ggplot(
    detail,
    ggplot2::aes(x = tile_x, y = tile_y, fill = FPR)
  ) +
    ggplot2::geom_tile(
      width = 0.44,
      height = 0.44,
      colour = "white",
      linewidth = 0.15
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = fpr_label, colour = label_colour),
      size = 2.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_tile(
      data = dplyr::distinct(detail, CellLine, design_x, endpoint_y),
      ggplot2::aes(x = design_x, y = endpoint_y),
      width = 0.94,
      height = 0.94,
      fill = NA,
      colour = "grey30",
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(CellLine ~ ., switch = "y", axes = "all_x") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(tpod_compact_design_labels),
      labels = unname(tpod_compact_design_labels),
      limits = c(0.5, length(tpod_compact_design_labels) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(tpod_analysis_levels),
      labels = rev(tpod_analysis_levels),
      limits = c(0.5, length(tpod_analysis_levels) + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#F1F2F3", "#FEE8C8", "#FB6A4A", "#B30000"),
      values = scales::rescale(c(0, 0.05, 0.25, 1)),
      limits = c(0, 1),
      breaks = c(0, 0.10, 0.25, 0.50, 1),
      labels = scales::label_percent(accuracy = 1),
      na.value = "#718096",
      name = "Observed FPR"
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::guides(
      custom = ggplot2::guide_custom(
        tpod_matrix_key(),
        title = "Matrix key (FC × R²)",
        width = grid::unit(4, "cm"),
        height = grid::unit(2.2, "cm"),
        order = 1
      ),
      fill = ggplot2::guide_colourbar(order = 2)
    ) +
    ggplot2::labs(
      title = "B   Experimental design and tPOD metric (FC × R² matrix)",
      subtitle = sprintf(
        "p = %.2f; each outlined cell shows four observed FPRs",
        p_main
      ),
      x = "Experimental design (dose groups × samples per group)",
      y = NULL
    ) +
    ggplot2::coord_fixed(ratio = 1.15) +
    ggplot2::theme_bw(base_size = 11.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(
        angle = 0,
        face = "bold",
        size = 14,
        hjust = 1,
        margin = ggplot2::margin(r = 8)
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        colour = "black",
        size = 10
      ),
      axis.text.y = ggplot2::element_text(colour = "black", size = 10.5),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(0.8, "lines")
    )

  patchwork::free(panel_a, type = "label", side = "l") /
    panel_b +
    patchwork::plot_layout(heights = c(0.28, 2.8))
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
  if (
    length(fdr_limit) != 1L ||
      !is.finite(fdr_limit) ||
      fdr_limit < 0 ||
      fdr_limit > 1
  ) {
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
        " | FC ",
        Fold_Change_Cutoff,
        " | p ",
        P_Value_Cutoff
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
      design_x = as.integer(Design) +
        case_when(
          P_Value_Cutoff == 0.01 & Fold_Change_Cutoff == 1.5 ~ -0.27,
          P_Value_Cutoff == 0.01 & Fold_Change_Cutoff == 2.0 ~ -0.09,
          P_Value_Cutoff == 0.05 & Fold_Change_Cutoff == 1.5 ~ 0.09,
          P_Value_Cutoff == 0.05 & Fold_Change_Cutoff == 2.0 ~ 0.27,
          TRUE ~ 0
        ),
      tpod_y = as.integer(tPoD) +
        if_else(
          R_Squared_Filter == "R-Sq > 0.6",
          0.13,
          -0.13
        )
    )

  setting_colours <- c(
    "#0072B2",
    "#56B4E9",
    "#009E73",
    "#7FC97F",
    "#D55E00",
    "#E69F00",
    "#CC79A7",
    "#000000"
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
        " | FC",
        Fold_Change_Cutoff,
        " | p",
        P_Value_Cutoff
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

build_s8_mean_variance_figure <- function(s8_summary) {
  plot_data <- if (is.list(s8_summary) && is.data.frame(s8_summary$data)) {
    s8_summary$data
  } else {
    s8_summary
  }
  required_columns <- c(
    "cell_line",
    "mean_log2cpm",
    "sd_log2cpm",
    "n",
    "count"
  )
  missing_columns <- setdiff(required_columns, names(plot_data))
  if (length(missing_columns) > 0L) {
    stop(
      "Figure S8 data are missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  plot_data <- plot_data %>%
    mutate(
      cell_line = factor(
        as.character(cell_line),
        levels = c("HepaRG", "Liver spheroids", "MCF-7", "TK6")
      ),
      n_panel = factor(
        ifelse(is.na(n), "Never", as.character(n)),
        levels = c(as.character(sample_size_levels), "Never")
      )
    )

  ggplot(plot_data, aes(x = mean_log2cpm, y = sd_log2cpm, colour = count)) +
    geom_point(alpha = 0.12, size = 0.7) +
    geom_smooth(
      aes(x = mean_log2cpm, y = sd_log2cpm, group = 1),
      inherit.aes = FALSE,
      method = "loess",
      formula = y ~ x,
      span = 0.75,
      se = FALSE,
      colour = "white",
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    scale_colour_gradient(
      low = "#132B43",
      high = "#56B1F7",
      na.value = "grey55",
      name = "Number of trials\ncalled DEG"
    ) +
    facet_grid(cell_line ~ n_panel, drop = FALSE) +
    labs(
      title = "Baseline mean expression and variability in DMSO controls",
      subtitle = paste(
        "LOESS curves are shown in white (span = 0.75)."
      ),
      x = "Mean log2 CPM",
      y = "Standard deviation of log2 CPM",
      caption = paste(
        "Columns n = 3–10 show probes called DEG in at least one of 1,000 trials;",
        "Never shows probes with no calls across the tested sample sizes."
      )
    ) +
    theme_grey(base_size = 10.5) +
    theme(
      panel.spacing = grid::unit(0.35, "lines"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey30", size = 9.5),
      plot.caption = element_text(colour = "grey30", size = 8, hjust = 0),
      strip.text = element_text(size = 9),
      legend.position = "right"
    )
}
