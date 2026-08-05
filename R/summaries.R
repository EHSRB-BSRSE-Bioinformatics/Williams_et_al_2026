library(dplyr)
library(tidyr)

parse_trial_names <- function(trial_names) {
  if (is.null(trial_names) || anyNA(trial_names) || any(!nzchar(trial_names))) {
    stop("Trial names must be non-missing, non-empty strings.")
  }

  parts <- strsplit(trial_names, "_", fixed = TRUE)
  if (any(lengths(parts) != 4L)) {
    stop("Trial names must follow '<cell_line>_<n_A>_<n_B>_<trial_id>'.")
  }

  parsed <- data.frame(
    run = trial_names,
    cell_line = vapply(parts, `[[`, character(1), 1L),
    n_per_group = suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 2L))),
    n_group_b = suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 3L))),
    trial_id = suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 4L))),
    stringsAsFactors = FALSE
  )
  if (anyNA(parsed[c("n_per_group", "n_group_b", "trial_id")])) {
    stop("Trial sample sizes and IDs must be integers.")
  }
  parsed
}

attach_trial_names <- function(trial_list, trial_metadata, cell_line) {
  if (is.null(names(trial_metadata))) {
    stop("The frozen trial metadata must be a named list.")
  }

  metadata_index <- parse_trial_names(names(trial_metadata))
  expected_names <- metadata_index$run[metadata_index$cell_line == cell_line]

  if (length(expected_names) == 0L) {
    stop("No frozen trial metadata found for cell line '", cell_line, "'.")
  }
  if (length(trial_list) != length(expected_names)) {
    stop(
      "The GSEA list for '", cell_line, "' has ", length(trial_list),
      " results, but the frozen metadata has ", length(expected_names),
      " trials. Refusing to attach trial names by position."
    )
  }

  existing_names <- names(trial_list)
  has_existing_names <- !is.null(existing_names) &&
    all(!is.na(existing_names) & nzchar(existing_names))
  if (!is.null(existing_names) && !has_existing_names) {
    stop(
      "The GSEA list for '", cell_line,
      "' has incomplete trial names; expected either all names or none."
    )
  }
  if (has_existing_names && !identical(existing_names, expected_names)) {
    stop(
      "Existing GSEA trial names for '", cell_line,
      "' do not match the frozen metadata order."
    )
  }

  names(trial_list) <- expected_names
  trial_list
}

factor_cell_lines <- function(data, column = "cell_line") {
  data[[column]] <- factor(data[[column]], levels = cell_line_levels)
  data
}

add_rank <- function(data, value_col, rank_col = "rank", group_cols = c("cell_line", "n")) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    mutate("{rank_col}" := rank(.data[[value_col]], ties.method = "first")) %>%
    ungroup()
}

false_positive_rates <- function(data,
                                 flag_col,
                                 n_col,
                                 group_cols = "cell_line",
                                 expected_trials = n_trials_default,
                                 rate_col = "false_positive_rate") {
  grouping_cols <- c(group_cols, n_col)

  rates <- data %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      n_trials_observed = n(),
      false_positive_trials = sum(.data[[flag_col]] %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      expected_trials = expected_trials,
      "{rate_col}" := false_positive_trials / n_trials_observed,
      FDR = .data[[rate_col]]
    )

  unexpected_trial_counts <- rates %>%
    filter(n_trials_observed != expected_trials)

  if (nrow(unexpected_trial_counts) > 0) {
    warning(
      "Some groups do not contain the expected number of Monte Carlo trials. ",
      "Rates were calculated using the observed trial count."
    )
  }

  rates %>%
    rename(trial_count = false_positive_trials) %>%
    ungroup() %>%
    arrange(across(all_of(grouping_cols)))
}

summarise_degs_from_trials <- function(trial_list) {
  final_degs <- lapply(trial_list, function(x) nrow(x$DESeqResDEGs))
  fdr_p_only <- lapply(trial_list, function(x) {
    res <- as.data.frame(x$DESeqRes)
    sum(!is.na(res$padj) & res$padj < 0.05)
  })
  p_and_fc <- lapply(trial_list, function(x) {
    res <- as.data.frame(x$DESeqRes)
    sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5))
  })
  rodaf_only <- lapply(trial_list, function(x) nrow(x$DESeqResFiltered))

  out <- bind_rows(
    data.frame(Filter = "Final DEGs", nDEGs = unlist(final_degs)),
    data.frame(Filter = "P-value filter only", nDEGs = unlist(fdr_p_only)),
    data.frame(Filter = "P-value + fold change", nDEGs = unlist(p_and_fc)),
    data.frame(Filter = "R-ODAF only", nDEGs = unlist(rodaf_only))
  )

  trial_info <- parse_trial_names(rep(names(trial_list), times = 4))
  bind_cols(out, trial_info) %>%
    mutate(contains_false_positive = nDEGs > 0) %>%
    group_by(cell_line, n_per_group, Filter) %>%
    mutate(rank = rank(nDEGs, ties.method = "first")) %>%
    ungroup() %>%
    factor_cell_lines()
}

summarise_gsea_from_trials <- function(gsea_list,
                                       min_size = gsea_min_size_default,
                                       max_size = gsea_max_size_default,
                                       alpha = 0.05,
                                       nes_cutoff = 2) {
  if (is.null(names(gsea_list)) || any(!nzchar(names(gsea_list)))) {
    stop(
      "GSEA results must have trial names. Use attach_trial_names() with ",
      "the frozen Monte Carlo trial metadata before summarising."
    )
  }
  if (min_size < 1L || max_size < min_size) {
    stop("GSEA size limits must satisfy 1 <= min_size <= max_size.")
  }

  rows <- Map(function(x, trial_name) {
    x <- as.data.frame(x)
    required_columns <- c("size", "pval", "NES")
    missing_columns <- setdiff(required_columns, names(x))
    if (length(missing_columns) > 0L) {
      stop(
        "GSEA result '", trial_name, "' is missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    x_filt <- x %>%
      dplyr::filter(
        !is.na(size),
        size >= min_size,
        size <= max_size,
        !is.na(pval)
      ) %>%
      dplyr::mutate(
        padj_recomputed = p.adjust(pval, method = "BH")
      )

    data.frame(
      num_rows_raw = nrow(x),
      num_rows = nrow(x_filt),
      min_size = min_size,
      max_size = max_size,
      p_lt_0.05 = sum(x_filt$padj_recomputed < alpha, na.rm = TRUE),
      p_lt_0.05_ES_gt_2 = sum(
        abs(x_filt$NES) > nes_cutoff &
          x_filt$padj_recomputed < alpha,
        na.rm = TRUE
      )
    )
  }, gsea_list, names(gsea_list))

  bind_rows(rows) %>%
    bind_cols(parse_trial_names(names(gsea_list))) %>%
    mutate(
      n = n_per_group,
      contains_false_positive_pval_only = p_lt_0.05 > 0,
      contains_false_positive_NES_pval_filtered = p_lt_0.05_ES_gt_2 > 0
    ) %>%
    add_rank("p_lt_0.05", "rank", c("cell_line", "n", "min_size")) %>%
    add_rank("p_lt_0.05_ES_gt_2", "rank_filtered", c("cell_line", "n", "min_size")) %>%
    factor_cell_lines()
}

summarise_gsea_files <- function(paths, trial_metadata) {
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("GSEA paths must be named with their publication cell-line labels.")
  }

  summaries <- Map(function(path, cell_line) {
    message("Summarising GSEA results for ", cell_line, ": ", path)
    gsea_list <- readRDS(path)
    gsea_list <- attach_trial_names(gsea_list, trial_metadata, cell_line)
    result <- summarise_gsea_from_trials(gsea_list)
    rm(gsea_list)
    gc()
    result
  }, unname(paths), names(paths))

  bind_rows(summaries) %>% factor_cell_lines()
}

summarise_ora_from_trials <- function(ora_list) {
  rows <- lapply(ora_list, function(x) {
    data.frame(
      go_bp = nrow(na.omit(as.data.frame(x$enrich_go_bp))),
      enricher_wikipathways = nrow(na.omit(as.data.frame(x$enricher_wikipathways)))
    )
  })

  bind_rows(rows) %>%
    bind_cols(parse_trial_names(names(ora_list))) %>%
    mutate(
      n = n_per_group,
      contains_GO_false_positive = go_bp > 0,
      contains_wiki_false_positive = enricher_wikipathways > 0
    ) %>%
    add_rank("go_bp", "go_bp_rank", c("cell_line", "n")) %>%
    add_rank("enricher_wikipathways", "wiki_rank", c("cell_line", "n")) %>%
    factor_cell_lines()
}

read_summary_inputs <- function() {
  list(
    deg = readRDS(derived_files$deg) %>% factor_cell_lines(),
    gsea = readRDS(derived_files$gsea) %>% factor_cell_lines(),
    ora = readRDS(derived_files$ora) %>% factor_cell_lines(),
    tpod_fdr = readRDS(derived_files$tpod_fdr),
    tgx_ddi = readRDS(derived_files$tgx_ddi) %>% factor_cell_lines()
  )
}
