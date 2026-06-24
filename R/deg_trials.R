# Trial-level DESeq2 analysis for DMSO-versus-DMSO Monte Carlo simulations.

make_deg_params <- function(alpha = 0.05,
                            linear_fc_filter = 1.5,
                            nmr_threshold = 100000,
                            cooks = FALSE,
                            design_formula = stats::as.formula("~ group")) {
  list(
    alpha = alpha,
    linear_fc_filter = linear_fc_filter,
    nmr_threshold = nmr_threshold,
    cooks = cooks,
    design_formula = design_formula
  )
}

linear_fold_change <- function(log2_fold_change) {
  ifelse(
    log2_fold_change >= 0,
    2 ^ log2_fold_change,
    -1 / (2 ^ log2_fold_change)
  )
}

run_deg_trial <- function(trial_name,
                          trial_metadata,
                          count_data,
                          params = make_deg_params()) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("The DESeq2 package is required to run DEG trials.")
  }

  trial <- trial_metadata[[trial_name]]
  prepared <- prepare_trial_counts(
    trial = trial,
    count_data = count_data,
    nmr_threshold = params$nmr_threshold
  )

  count_subset <- round(prepared$count_data)
  storage.mode(count_subset) <- "integer"
  exp_metadata <- prepared$exp_metadata

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_subset,
    colData = exp_metadata,
    design = params$design_formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  res <- DESeq2::results(
    dds,
    contrast = c("group", "B", "A"),
    alpha = params$alpha,
    cooksCutoff = params$cooks
  )

  res_df <- as.data.frame(res)
  res_df$linearFoldChange <- linear_fold_change(res_df$log2FoldChange)
  res_df$pass_p_value <- !is.na(res_df$padj) & res_df$padj < params$alpha
  res_df$pass_fold_change <- !is.na(res_df$log2FoldChange) &
    abs(res_df$log2FoldChange) > log2(params$linear_fc_filter)
  res_df$pass_final_filter <- res_df$pass_p_value & res_df$pass_fold_change

  # The original exploratory workflow referred to an R-ODAF-filtered table, but
  # the private R-ODAF implementation is not present in the material supplied to
  # this repository. For the reproducible public rerun, DESeqResFiltered records
  # the analysable rows with finite adjusted p-values; DESeqResDEGs records the
  # publication figure's p-value + fold-change false-positive calls.
  filtered <- res_df[!is.na(res_df$padj), , drop = FALSE]
  final_degs <- res_df[res_df$pass_final_filter %in% TRUE, , drop = FALSE]

  list(
    DESeqRes = res,
    DESeqResFiltered = filtered,
    DESeqResDEGs = final_degs,
    exp_metadata = exp_metadata,
    contrasts = data.frame(group_A = "A", group_B = "B"),
    num_rows_metadata = nrow(exp_metadata),
    trial_name = trial_name,
    params = params
  )
}

run_deg_trials <- function(trial_metadata,
                           count_data,
                           params = make_deg_params(),
                           cores = 1L) {
  validate_trial_metadata(trial_metadata)
  trial_names <- names(trial_metadata)
  out <- run_trials_parallel(
    trial_names,
    run_deg_trial,
    trial_metadata = trial_metadata,
    count_data = count_data,
    params = params,
    cores = cores
  )
  names(out) <- trial_names
  out
}
