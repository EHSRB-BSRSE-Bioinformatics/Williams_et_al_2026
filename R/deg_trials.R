# Study-specific DESeq2/R-ODAF analysis. See docs/deseq2_methods.md for the
# execution order, historical indexing behavior, and validation scope.

make_deg_params <- function(alpha = 0.05,
                            linear_fc_filter = 1.5,
                            nmr_threshold = 100000,
                            cooks = FALSE,
                            design_formula = stats::as.formula("~ group"),
                            min_cpm = 0.5) {
  list(alpha = alpha, linear_fc_filter = linear_fc_filter,
       nmr_threshold = nmr_threshold, cooks = cooks,
       design_formula = design_formula, min_cpm = min_cpm)
}

linear_fold_change <- function(log2_fold_change) {
  ifelse(log2_fold_change >= 0, 2 ^ log2_fold_change,
         -1 / (2 ^ log2_fold_change))
}

require_deg_packages <- function() {
  packages <- c("DESeq2", "edgeR", "ashr")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing DEG packages: ", paste(missing, collapse = ", "),
         ". Use the optional Pixi environment: pixi run -e deg-study ...")
  }
}

rodaf_relevance_filter <- function(cpm, groups, min_cpm) {
  passes <- vapply(levels(groups), function(group) {
    selected <- groups == group
    rowSums(cpm[, selected, drop = FALSE] >= min_cpm) >= 0.75 * sum(selected)
  }, logical(nrow(cpm)))
  rowSums(matrix(passes, nrow = nrow(cpm))) > 0
}

rodaf_post_test_filters <- function(normalized_counts, groups) {
  # Only significant genes enter these two rules in the study implementation.
  group_columns <- lapply(levels(groups), function(group) which(groups == group))
  filters <- matrix(1, nrow = nrow(normalized_counts), ncol = 3,
                    dimnames = list(rownames(normalized_counts),
                                    c("Low", "quantile", "spike")))
  for (gene in seq_len(nrow(normalized_counts))) {
    group_counts <- lapply(group_columns, function(columns) {
      as.numeric(normalized_counts[gene, columns])
    })
    filters[gene, "quantile"] <- as.integer(
      stats::median(group_counts[[1L]]) > stats::quantile(group_counts[[2L]], 0.75)[[1L]] ||
      stats::median(group_counts[[2L]]) > stats::quantile(group_counts[[1L]], 0.75)[[1L]]
    )
    spikes <- vapply(group_counts, function(counts) {
      if (max(counts) == 0) return(FALSE)
      max(counts) / sum(counts) >= 1.4 * length(counts) ^ (-0.66)
    }, logical(1))
    # The historical implementation rejects a gene only if BOTH groups spike.
    filters[gene, "spike"] <- as.integer(sum(spikes) <= 1)
  }
  filters
}

run_deg_trial <- function(trial_name,
                          trial_metadata,
                          count_data,
                          params = make_deg_params()) {
  require_deg_packages()
  trial <- trial_metadata[[trial_name]]
  if (is.null(trial)) stop("Unknown trial: ", trial_name)

  # Historically, samples below the read-depth cutoff were dropped from the
  # trial rather than raising the missing-sample error of prepare_trial_counts().
  count_data <- count_data[, colSums(count_data) > params$nmr_threshold, drop = FALSE]
  trial <- trial[trial$sample_id %in% colnames(count_data), , drop = FALSE]
  trial$original_names <- as.character(trial$sample_id)
  trial$group <- factor(as.character(trial$group),
                        levels = sort(unique(as.character(trial$group))))
  count_subset <- count_data[, trial$original_names, drop = FALSE]
  stopifnot(ncol(count_subset) > 0, nrow(trial) > 0)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_subset), colData = as.data.frame(trial),
    design = params$design_formula
  )
  # This was unconditional, despite filter_gene_counts = FALSE in the old script.
  dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]
  dds <- DESeq2::DESeq(dds, parallel = FALSE)
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  cpm <- edgeR::cpm(normalized_counts)
  relevant <- rodaf_relevance_filter(cpm, trial$group, params$min_cpm)
  relevant_counts <- normalized_counts[relevant, , drop = FALSE]

  contrast <- c("group", "A", "B")
  res <- DESeq2::results(
    dds[rownames(relevant_counts), ], contrast = contrast,
    alpha = params$alpha, pAdjustMethod = "fdr",
    cooksCutoff = params$cooks, parallel = FALSE
  )
  res <- DESeq2::lfcShrink(dds, contrast = contrast, res = res,
                          type = "ashr", parallel = FALSE)
  significant <- subset(res, res$padj < params$alpha)

  filtered <- final_degs <- data.frame()
  filtered_out <- list()
  filter_table <- data.frame()
  merged_degs <- NULL
  if (nrow(significant) > 0L) {
    significant_counts <- relevant_counts[
      rownames(relevant_counts) %in% rownames(significant), , drop = FALSE
    ]
    filters <- rodaf_post_test_filters(significant_counts, trial$group)
    passes <- rowSums(filters) == 3

    # Preserve the original positional logical subsetting. The logical vector
    # is sized to significant genes; R recycles it across the full result table.
    # Replacing this with matching by feature ID would change archived results.
    filtered <- res[passes, ]
    final_degs <- significant[
      passes & !is.na(significant$padj) &
        abs(significant$log2FoldChange) > log2(params$linear_fc_filter),
    ]
    no_quantile <- significant[filters[, "quantile"] == 0, ]
    spikes <- significant[filters[, "spike"] == 0, ]
    filtered_out <- list(normalized_counts[!relevant, , drop = FALSE],
                         no_quantile, spikes)
    merged_degs <- rownames(final_degs)
    filter_table <- data.frame(
      facet = NA, contrast = "A_vs_B", initial = nrow(normalized_counts),
      relevance_filtered = sum(!relevant), quantile_filtered = nrow(no_quantile),
      spike_filtered = nrow(spikes), passed_all_filters = nrow(final_degs)
    )
  }

  list(
    DESeqRes = res,
    DESeqResFiltered = filtered,
    DESeqResDEGs = final_degs,
    DESeqmergedDEGs = merged_degs,
    DESeqfiltered_table = filter_table,
    DESeqjs_filtered_out = filtered_out,
    exp_metadata = trial,
    contrasts = data.frame(group_A = "A", group_B = "B"),
    intgroup = "group",
    num_rows_metadata = nrow(trial),
    trial_name = trial_name,
    params = params
  )
}

run_deg_trials <- function(trial_metadata,
                           count_data,
                           params = make_deg_params(),
                           cores = 1L) {
  validate_trial_metadata(trial_metadata)
  require_deg_packages()
  trial_names <- names(trial_metadata)
  out <- run_trials_parallel(
    trial_names, run_deg_trial, trial_metadata = trial_metadata,
    count_data = count_data, params = params, cores = cores
  )
  names(out) <- trial_names
  failed <- vapply(out, inherits, logical(1), "try-error")
  if (any(failed)) stop("DEG trial(s) failed: ", paste(trial_names[failed], collapse = ", "))
  out
}
