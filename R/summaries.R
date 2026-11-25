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
    n_per_group = suppressWarnings(as.integer(vapply(
      parts,
      `[[`,
      character(1),
      2L
    ))),
    n_group_b = suppressWarnings(as.integer(vapply(
      parts,
      `[[`,
      character(1),
      3L
    ))),
    trial_id = suppressWarnings(as.integer(vapply(
      parts,
      `[[`,
      character(1),
      4L
    ))),
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
      "The GSEA list for '",
      cell_line,
      "' has ",
      length(trial_list),
      " results, but the frozen metadata has ",
      length(expected_names),
      " trials. Refusing to attach trial names by position."
    )
  }

  existing_names <- names(trial_list)
  has_existing_names <- !is.null(existing_names) &&
    all(!is.na(existing_names) & nzchar(existing_names))
  if (!is.null(existing_names) && !has_existing_names) {
    stop(
      "The GSEA list for '",
      cell_line,
      "' has incomplete trial names; expected either all names or none."
    )
  }
  if (has_existing_names && !identical(existing_names, expected_names)) {
    stop(
      "Existing GSEA trial names for '",
      cell_line,
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

add_rank <- function(
  data,
  value_col,
  rank_col = "rank",
  group_cols = c("cell_line", "n")
) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    mutate("{rank_col}" := rank(.data[[value_col]], ties.method = "first")) %>%
    ungroup()
}

false_positive_rates <- function(
  data,
  flag_col,
  n_col,
  group_cols = "cell_line",
  expected_trials = n_trials_default,
  rate_col = "false_positive_rate"
) {
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
    sum(
      !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5)
    )
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

summarise_gsea_from_trials <- function(
  gsea_list,
  min_size = gsea_min_size_default,
  max_size = gsea_max_size_default,
  alpha = 0.05,
  nes_cutoff = 2
) {
  if (is.null(names(gsea_list)) || any(!nzchar(names(gsea_list)))) {
    stop(
      "GSEA results must have trial names. Use attach_trial_names() with ",
      "the frozen Monte Carlo trial metadata before summarising."
    )
  }
  if (min_size < 1L || max_size < min_size) {
    stop("GSEA size limits must satisfy 1 <= min_size <= max_size.")
  }

  rows <- Map(
    function(x, trial_name) {
      x <- as.data.frame(x)
      required_columns <- c("size", "pval", "NES")
      missing_columns <- setdiff(required_columns, names(x))
      if (length(missing_columns) > 0L) {
        stop(
          "GSEA result '",
          trial_name,
          "' is missing columns: ",
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
    },
    gsea_list,
    names(gsea_list)
  )

  bind_rows(rows) %>%
    bind_cols(parse_trial_names(names(gsea_list))) %>%
    mutate(
      n = n_per_group,
      contains_false_positive_pval_only = p_lt_0.05 > 0,
      contains_false_positive_NES_pval_filtered = p_lt_0.05_ES_gt_2 > 0
    ) %>%
    add_rank("p_lt_0.05", "rank", c("cell_line", "n", "min_size")) %>%
    add_rank(
      "p_lt_0.05_ES_gt_2",
      "rank_filtered",
      c("cell_line", "n", "min_size")
    ) %>%
    factor_cell_lines()
}

summarise_gsea_files <- function(paths, trial_metadata) {
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("GSEA paths must be named with their publication cell-line labels.")
  }

  summaries <- Map(
    function(path, cell_line) {
      message("Summarising GSEA results for ", cell_line, ": ", path)
      gsea_list <- readRDS(path)
      gsea_list <- attach_trial_names(gsea_list, trial_metadata, cell_line)
      result <- summarise_gsea_from_trials(gsea_list)
      rm(gsea_list)
      gc()
      result
    },
    unname(paths),
    names(paths)
  )

  bind_rows(summaries) %>% factor_cell_lines()
}

summarise_ora_from_trials <- function(ora_list) {
  rows <- lapply(ora_list, function(x) {
    data.frame(
      go_bp = nrow(na.omit(as.data.frame(x$enrich_go_bp))),
      enricher_wikipathways = nrow(na.omit(as.data.frame(
        x$enricher_wikipathways
      )))
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

attach_reactome_ora <- function(ora_summary, reactome_trials) {
  required <- c("run", "cell_line", "n", "n_significant", "status")
  if (
    !all(required %in% names(reactome_trials)) ||
      anyNA(reactome_trials$run) ||
      anyDuplicated(reactome_trials$run) ||
      anyNA(ora_summary$run) ||
      anyDuplicated(ora_summary$run) ||
      !setequal(ora_summary$run, reactome_trials$run)
  ) {
    stop(
      "Reactome ORA results must have unique trial IDs matching the ORA summary."
    )
  }
  reactome_trials <- reactome_trials[
    match(ora_summary$run, reactome_trials$run),
  ]
  if (
    anyNA(reactome_trials$n_significant) ||
      any(!is.finite(reactome_trials$n_significant)) ||
      any(reactome_trials$n_significant < 0) ||
      any(
        reactome_trials$n_significant != floor(reactome_trials$n_significant)
      ) ||
      anyNA(reactome_trials$status) ||
      any(
        !reactome_trials$status %in%
          c("significant", "none", "insufficient_input")
      )
  ) {
    stop("Reactome ORA contains invalid counts or unresolved trial statuses.")
  }
  if (
    !identical(
      as.character(ora_summary$cell_line),
      as.character(reactome_trials$cell_line)
    ) ||
      !identical(as.integer(ora_summary$n), as.integer(reactome_trials$n)) ||
      any(
        (reactome_trials$status == "significant") !=
          (reactome_trials$n_significant > 0)
      )
  ) {
    stop(
      "Reactome ORA trial metadata or statuses disagree with the trial counts."
    )
  }
  ora_summary %>%
    ungroup() %>%
    mutate(
      reactome = reactome_trials$n_significant,
      reactome_status = reactome_trials$status,
      contains_reactome_false_positive = reactome > 0
    ) %>%
    add_rank("reactome", "reactome_rank", c("cell_line", "n"))
}

# Canonical plotting data from feature-aligned archived counts at fixed cutoffs.
read_deg_figure_summary <- function(
  path = fdr_path(
    "data",
    "derived",
    "deg_thresholds",
    "threshold_sensitivity.rds"
  )
) {
  if (!file.exists(path)) {
    stop(
      "Missing validated DEG counts: ",
      path,
      ". Run scripts/deg_thresholds/analyze.R; no model refitting is needed."
    )
  }
  fixed <- readRDS(path)
  counts <- fixed$trial_counts
  columns <- c(
    "screening__p0.05_fcnone",
    "screening__p0.05_fc1.5",
    "full_rodaf__p0.05_fc1.5"
  )
  frozen <- readRDS(derived_files$trial_metadata)
  ids <- names(frozen)
  if (
    length(ids) != 32000L ||
      anyDuplicated(ids) ||
      is.null(rownames(counts)) ||
      anyDuplicated(rownames(counts)) ||
      !setequal(rownames(counts), ids) ||
      !all(columns %in% colnames(counts)) ||
      is.null(fixed$validation$passed) ||
      !all(fixed$validation$passed)
  ) {
    stop("DEG figure counts lack complete validated frozen trial coverage.")
  }
  counts <- counts[match(ids, rownames(counts)), columns, drop = FALSE]
  if (
    anyNA(counts) ||
      any(!is.finite(counts)) ||
      any(counts < 0) ||
      any(counts != floor(counts)) ||
      any(counts[, 2] > counts[, 1]) ||
      any(counts[, 3] > counts[, 2])
  ) {
    stop("Invalid DEG figure counts or final calls outside p-only selection.")
  }
  index <- parse_trial_names(ids)
  if (
    any(index$n_per_group != index$n_group_b) ||
      !setequal(unique(index$cell_line), cell_line_levels) ||
      !setequal(unique(index$n_per_group), sample_size_levels) ||
      any(table(index$cell_line, index$n_per_group) != 1000L)
  ) {
    stop("DEG figure strata do not contain the frozen 1000-trial denominators.")
  }
  archived <- as.data.frame(readRDS(derived_files$deg))
  archived <- archived[
    !is.na(archived$Filter) & archived$Filter == "Final DEGs",
  ]
  if (
    anyDuplicated(archived$name) ||
      !setequal(archived$name, ids) ||
      !identical(
        as.integer(archived$nDEGs[match(ids, archived$name)]),
        as.integer(counts[, 3])
      )
  ) {
    stop("DEG figure final counts differ from the frozen publication summary.")
  }
  labels <- c("P-value filter only", "P-value + fold-change", "Final DEGs")
  sources <- c(
    "archived_DESeqRes_padj",
    "archived_DESeqRes_padj_lfc",
    "verified_archived_final"
  )
  rows <- lapply(seq_along(columns), function(i) {
    data.frame(
      name = ids,
      n_per_group = index$n_per_group,
      cell_line = index$cell_line,
      Filter = labels[i],
      nDEGs = counts[, i],
      contains_false_positive = counts[, i] > 0,
      count_source = sources[i]
    )
  })
  bind_rows(rows) %>%
    mutate(Filter = factor(Filter, levels = labels)) %>%
    factor_cell_lines() %>%
    group_by(cell_line, Filter, n_per_group) %>%
    mutate(rank = rank(nDEGs, ties.method = "first")) %>%
    ungroup()
}

read_summary_inputs <- function() {
  gsea_archived <- readRDS(derived_files$gsea) %>% factor_cell_lines()
  gsea <- readRDS(derived_files$gsea_unique) %>% factor_cell_lines()
  archived_rows <- match(gsea$run, gsea_archived$run)
  if (
    nrow(gsea) != nrow(gsea_archived) ||
      anyDuplicated(gsea$run) ||
      anyNA(archived_rows) ||
      !identical(gsea$cell_line, gsea_archived$cell_line[archived_rows]) ||
      !identical(
        as.integer(gsea$n),
        as.integer(gsea_archived$n[archived_rows])
      ) ||
      anyNA(gsea$p_lt_0.05) ||
      anyNA(gsea$p_lt_0.05_ES_gt_2)
  ) {
    stop(
      "Corrected GSEA summary does not match the validated frozen trial groups."
    )
  }
  list(
    deg = read_deg_figure_summary(),
    gsea = gsea,
    gsea_archived = gsea_archived,
    ora = attach_reactome_ora(
      readRDS(derived_files$ora),
      readRDS(derived_files$reactome_ora)
    ) %>%
      factor_cell_lines(),
    tpod_fdr = readRDS(derived_files$tpod_fdr),
    tgx_ddi = readRDS(derived_files$tgx_ddi) %>% factor_cell_lines()
  )
}

s8_cache_stems <- c(
  TK6 = "tk6",
  `MCF-7` = "mcf7",
  HepaRG = "heparg",
  `Liver spheroids` = "spheroids"
)

# Reconstruct Figure S8 from the historical p-value-only recurrence lists and
# the QC-applied DMSO count matrix. Missing recurrence rows intentionally remain
# missing, matching the original sparse join rather than imputing zero counts.
build_s8_figure_summary <- function(
  count_path = fdr_path("data", "input", "count_table.tsv"),
  metadata_path = fdr_path("data", "input", "metadata.QC_applied.txt"),
  trial_metadata_path = derived_files$trial_metadata,
  cache_dir = Sys.getenv(
    "FDR_S8_RECURRENCE_DIR",
    unset = fdr_path("data", "large", "gene_recurrence")
  )
) {
  required_paths <- c(count_path, metadata_path, trial_metadata_path)
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "Missing Figure S8 input files: ",
      paste(missing_paths, collapse = ", ")
    )
  }

  metadata <- read.delim(
    metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_metadata <- c("Sample_ID", "chemical", "cell_line")
  missing_columns <- setdiff(required_metadata, names(metadata))
  if (length(missing_columns) > 0L || anyDuplicated(metadata$Sample_ID)) {
    stop(
      "Figure S8 metadata must have unique sample IDs and columns: ",
      paste(required_metadata, collapse = ", ")
    )
  }

  counts <- as.matrix(read.delim(
    count_path,
    check.names = FALSE,
    row.names = 1
  ))
  storage.mode(counts) <- "numeric"
  if (
    is.null(rownames(counts)) ||
      is.null(colnames(counts)) ||
      anyDuplicated(rownames(counts)) ||
      anyDuplicated(colnames(counts)) ||
      anyNA(counts) ||
      any(!is.finite(counts)) ||
      any(counts < 0)
  ) {
    stop(
      "Figure S8 counts must have unique feature/sample IDs and finite non-negative values."
    )
  }

  frozen <- readRDS(trial_metadata_path)
  if (is.null(names(frozen)) || anyDuplicated(names(frozen))) {
    stop("Figure S8 requires unique names in the frozen trial metadata.")
  }
  trial_index <- parse_trial_names(names(frozen))
  if (
    nrow(trial_index) !=
      length(cell_line_levels) *
        length(sample_size_levels) *
        n_trials_default ||
      any(trial_index$n_per_group != trial_index$n_group_b) ||
      !setequal(unique(trial_index$cell_line), cell_line_levels) ||
      !setequal(unique(trial_index$n_per_group), sample_size_levels) ||
      any(
        table(trial_index$cell_line, trial_index$n_per_group) !=
          n_trials_default
      )
  ) {
    stop(
      "Figure S8 trial metadata does not contain the complete frozen design."
    )
  }

  expression <- lapply(cell_line_levels, function(cell_line) {
    sample_ids <- metadata$Sample_ID[
      !is.na(metadata$chemical) &
        metadata$chemical == "DMSO" &
        metadata$cell_line == cell_line
    ]
    sample_columns <- match(sample_ids, colnames(counts))
    if (
      length(sample_ids) < 2L ||
        anyNA(sample_columns) ||
        anyDuplicated(sample_ids)
    ) {
      stop("Missing or invalid DMSO count columns for ", cell_line, ".")
    }

    cell_counts <- counts[, sample_columns, drop = FALSE]
    library_sizes <- colSums(cell_counts)
    if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
      stop("DMSO samples have invalid library sizes for ", cell_line, ".")
    }

    # CPM after DESeq2 size-factor normalization is identical to raw-count CPM:
    # the per-sample size factor cancels when CPM divides by library size.
    log2_cpm <- log2(sweep(cell_counts, 2, library_sizes / 1e6, "/") + 1)
    data.frame(
      cell_line = cell_line,
      feature = rownames(counts),
      mean_log2cpm = rowMeans(log2_cpm),
      sd_log2cpm = apply(log2_cpm, 1L, stats::sd),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  recurrence <- lapply(cell_line_levels, function(cell_line) {
    cache_path <- file.path(
      cache_dir,
      paste0(s8_cache_stems[[cell_line]], ".rds")
    )
    if (!file.exists(cache_path)) {
      stop(
        "Missing Figure S8 recurrence cache: ",
        cache_path,
        ". Rebuild the legacy recurrence cache before rebuilding the summary."
      )
    }
    cache <- readRDS(cache_path)
    expected_ids <- trial_index$run[trial_index$cell_line == cell_line]
    if (
      !identical(cache$cell, cell_line) ||
        is.null(cache$screen) ||
        is.null(names(cache$screen)) ||
        anyDuplicated(names(cache$screen)) ||
        !setequal(names(cache$screen), expected_ids) ||
        length(cache$screen) != length(expected_ids) ||
        length(cache$archive) != 1L ||
        is.na(cache$archive) ||
        length(cache$sha256) != 1L ||
        is.na(cache$sha256)
    ) {
      stop(
        "Figure S8 recurrence cache does not match frozen trials for ",
        cell_line,
        "."
      )
    }

    screens <- cache$screen[expected_ids]
    invalid_screen <- vapply(
      screens,
      function(features) {
        !is.character(features) ||
          anyNA(features) ||
          any(!nzchar(features)) ||
          anyDuplicated(features) > 0L
      },
      logical(1)
    )
    if (any(invalid_screen)) {
      stop(
        "Invalid or duplicated DEG feature IDs in the ",
        cell_line,
        " recurrence cache."
      )
    }
    unknown_features <- setdiff(
      unique(unlist(screens, use.names = FALSE)),
      rownames(counts)
    )
    if (length(unknown_features) > 0L) {
      stop(
        "Recurrence cache contains features absent from the count matrix for ",
        cell_line,
        "."
      )
    }

    by_sample_size <- lapply(sample_size_levels, function(sample_size) {
      trial_ids <- trial_index$run[
        trial_index$cell_line == cell_line &
          trial_index$n_per_group == sample_size
      ]
      calls <- unlist(screens[trial_ids], use.names = FALSE)
      if (!length(calls)) {
        stop(
          "No p-value-only DEG calls were found for ",
          cell_line,
          ", n = ",
          sample_size,
          "."
        )
      }
      counts_by_feature <- table(calls)
      data.frame(
        cell_line = cell_line,
        n = as.integer(sample_size),
        feature = names(counts_by_feature),
        count = as.integer(counts_by_feature),
        stringsAsFactors = FALSE
      )
    }) %>%
      bind_rows()

    if (
      any(by_sample_size$count < 1L) ||
        any(by_sample_size$count > n_trials_default) ||
        anyDuplicated(by_sample_size[c("cell_line", "n", "feature")])
    ) {
      stop("Invalid per-probe DEG recurrence counts for ", cell_line, ".")
    }
    attr(by_sample_size, "source_archive") <- data.frame(
      cell_line = cell_line,
      path = cache$archive,
      sha256 = cache$sha256,
      stringsAsFactors = FALSE
    )
    by_sample_size
  })

  source_archives <- bind_rows(lapply(
    recurrence,
    attr,
    which = "source_archive"
  ))
  recurrence <- bind_rows(recurrence)
  figure_data <- dplyr::left_join(
    expression,
    recurrence,
    by = c("cell_line", "feature")
  )
  if (
    anyDuplicated(figure_data[c("cell_line", "feature", "n")]) ||
      any(!is.finite(figure_data$mean_log2cpm)) ||
      any(!is.finite(figure_data$sd_log2cpm))
  ) {
    stop(
      "Figure S8 data contain duplicate probe panels or invalid expression summaries."
    )
  }

  list(
    data = figure_data,
    source_archives = source_archives,
    method = list(
      expression = paste(
        "Mean and sample SD of log2(CPM + 1) across all QC-applied DMSO wells",
        "within each cell model."
      ),
      recurrence = paste(
        "Number of archived DMSO-versus-DMSO trials with BH-adjusted p < 0.05,",
        "stratified by per-group sample size; no fold-change filter."
      ),
      join = paste(
        "Matches the original sparse join: only positive recurrence counts are",
        "shown in n = 3–10 panels; probes never called DEG are shown once with n = NA."
      )
    )
  )
}

read_s8_figure_summary <- function(path = derived_files$s8_mean_variance) {
  if (!file.exists(path)) {
    stop(
      "Missing Figure S8 summary: ",
      path,
      ". Run scripts/rebuild_s8_figure_data.R when the recurrence caches are available."
    )
  }
  summary <- readRDS(path)
  required_columns <- c(
    "cell_line",
    "feature",
    "mean_log2cpm",
    "sd_log2cpm",
    "n",
    "count"
  )
  if (!is.list(summary) || !is.data.frame(summary$data)) {
    stop("Figure S8 summary must contain a data frame named 'data'.")
  }
  missing_columns <- setdiff(required_columns, names(summary$data))
  data <- summary$data
  if (
    length(missing_columns) > 0L ||
      nrow(data) == 0L ||
      anyNA(data$cell_line) ||
      anyNA(data$feature) ||
      any(!data$cell_line %in% cell_line_levels) ||
      !all(cell_line_levels %in% unique(as.character(data$cell_line))) ||
      any(!is.finite(data$mean_log2cpm)) ||
      any(!is.finite(data$sd_log2cpm)) ||
      any(data$sd_log2cpm < 0) ||
      any(!is.na(data$n) & !data$n %in% sample_size_levels) ||
      any(is.na(data$n) != is.na(data$count)) ||
      any(
        !is.na(data$count) &
          (!is.finite(data$count) |
            data$count < 1L |
            data$count > n_trials_default |
            data$count != floor(data$count))
      ) ||
      anyDuplicated(data[c("cell_line", "feature", "n")])
  ) {
    stop("Figure S8 summary has missing, invalid, or duplicate plot data.")
  }
  if (
    !is.data.frame(summary$source_archives) ||
      !all(
        c("cell_line", "path", "sha256") %in% names(summary$source_archives)
      ) ||
      !setequal(summary$source_archives$cell_line, cell_line_levels)
  ) {
    stop("Figure S8 summary is missing per-cell-model archive provenance.")
  }
  summary
}
