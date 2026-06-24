# Trial-level TGx-DDI biomarker calls for DMSO-versus-DMSO simulations.

read_tgx_ddi_classifier <- function(path) {
  if (!file.exists(path)) {
    stop("TGx-DDI classifier file does not exist: ", path)
  }
  classifier <- read.delim(path, check.names = FALSE)
  required <- c(
    "GENE_SYMBOL",
    "Genotoxic.score",
    "Non.Genotoxic.score",
    "stdev",
    "pc.mn",
    "pc.std",
    "PC1"
  )
  missing <- setdiff(required, names(classifier))
  if (length(missing) > 0L) {
    stop("TGx-DDI classifier is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (any(classifier$stdev <= 0 | classifier$pc.std <= 0, na.rm = TRUE)) {
    stop("TGx-DDI classifier stdev and pc.std columns must be positive.")
  }
  classifier
}

make_tgx_ddi_input <- function(count_subset, classifier) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR is required to calculate CPM values for TGx-DDI trials.")
  }

  gene <- sub("_.*", "", rownames(count_subset))
  collapsed <- rowsum(as.matrix(count_subset), group = gene, reorder = FALSE)
  log_cpm <- log2(edgeR::cpm(collapsed) + 1)
  log_cpm <- log_cpm[rownames(log_cpm) %in% classifier$GENE_SYMBOL, , drop = FALSE]

  data.frame(
    SampleID = rownames(log_cpm),
    as.data.frame(log_cpm, check.names = FALSE),
    check.names = FALSE,
    row.names = NULL
  )
}

prepare_tgx_ddi_trial <- function(trial_name,
                                  trial_metadata,
                                  count_data,
                                  classifier,
                                  nmr_threshold = 100000) {
  trial <- trial_metadata[[trial_name]]
  prepared <- prepare_trial_counts(
    trial = trial,
    count_data = count_data,
    nmr_threshold = nmr_threshold
  )

  biomarker_input <- make_tgx_ddi_input(prepared$count_data, classifier)
  list(
    trial_name = trial_name,
    biomarker_input = biomarker_input,
    exp_metadata = prepared$exp_metadata,
    num_rows_metadata = nrow(prepared$exp_metadata)
  )
}

tgx_ddi_hca_call <- function(z) {
  if (ncol(z) < 40L) {
    return(NA_character_)
  }

  training_cols <- 11:40
  hca_data <- z[, training_cols, drop = FALSE]
  numeric_cols <- vapply(hca_data, is.numeric, logical(1))
  if (!all(numeric_cols)) {
    return(NA_character_)
  }

  dendro <- stats::dist(t(hca_data))
  hc <- stats::hclust(dendro, method = "average")
  clus <- stats::cutree(hc, k = 2)

  # This preserves the Li et al. TGx-DDI classifier convention used in the
  # exploratory script: columns 11:23 are DDI references, 24:39 are non-DDI
  # references, and column 40 is the trial profile appended as lfc.
  clus <- clus[clus[-30L] == clus[30L]]
  if (length(clus) == 13L && identical(names(clus), colnames(z)[11:23])) {
    return("DDI")
  }
  if (length(clus) == 16L && identical(names(clus), colnames(z)[24:39])) {
    return("NonDDI")
  }
  "unknown"
}

classify_tgx_ddi_trial <- function(prepared_trial,
                                   classifier,
                                   probability_threshold = 0.9,
                                   pc1_threshold = 0,
                                   prior_prob_genotoxic = 0.5) {
  trial <- prepared_trial$trial_name
  exp_metadata <- prepared_trial$exp_metadata
  biomarker_data <- prepared_trial$biomarker_input

  a_names <- exp_metadata$sample_id[exp_metadata$group == "A"]
  b_names <- exp_metadata$sample_id[exp_metadata$group == "B"]
  if (length(a_names) == 0L || length(b_names) == 0L) {
    stop("Trial must contain both A and B groups: ", trial)
  }
  missing_biomarker_samples <- setdiff(c(a_names, b_names), names(biomarker_data))
  if (length(missing_biomarker_samples) > 0L) {
    stop(
      "TGx-DDI biomarker input is missing trial samples: ",
      paste(missing_biomarker_samples, collapse = ", ")
    )
  }

  biomarker_matrix <- as.matrix(biomarker_data[, setdiff(names(biomarker_data), "SampleID"), drop = FALSE])
  rownames(biomarker_matrix) <- biomarker_data$SampleID

  trial_fcs <- data.frame(
    ID = rownames(biomarker_matrix),
    lfc = rowMeans(biomarker_matrix[, b_names, drop = FALSE]) -
      rowMeans(biomarker_matrix[, a_names, drop = FALSE]),
    stringsAsFactors = FALSE
  )
  trial_fcs$lfc <- log10(2 ^ trial_fcs$lfc)

  z <- merge(classifier, trial_fcs, by.x = "GENE_SYMBOL", by.y = "ID")
  if (nrow(z) == 0L) {
    stop("No TGx-DDI classifier genes were present in trial: ", trial)
  }

  g_score <- sum(((z$lfc - z$Genotoxic.score) / z$stdev) ^ 2) -
    2 * log(prior_prob_genotoxic)
  ng_score <- sum(((z$lfc - z$Non.Genotoxic.score) / z$stdev) ^ 2) -
    2 * log(1 - prior_prob_genotoxic)

  prob_g <- 1 / (1 + exp(g_score / 2 - ng_score / 2))
  prob_ng <- 1 - prob_g

  pc_x <- as.matrix((z$lfc - z$pc.mn) / z$pc.std)
  pc1 <- as.numeric(as.vector(z$PC1) %*% pc_x)
  hca <- tgx_ddi_hca_call(z)

  list(
    trial_res = trial,
    probG_res = prob_g,
    probNG_res = prob_ng,
    pc1_res = pc1,
    hca = hca,
    n_genes_used = nrow(z),
    positive_probability = prob_g > probability_threshold,
    positive_pc1 = pc1 < pc1_threshold
  )
}

run_tgx_ddi_trials <- function(trial_metadata,
                               count_data,
                               classifier,
                               nmr_threshold = 100000,
                               cores = 1L,
                               probability_threshold = 0.9,
                               pc1_threshold = 0) {
  validate_trial_metadata(trial_metadata)
  trial_names <- names(trial_metadata)

  out <- run_trials_parallel(
    trial_names,
    function(trial_name) {
      prepared <- prepare_tgx_ddi_trial(
        trial_name = trial_name,
        trial_metadata = trial_metadata,
        count_data = count_data,
        classifier = classifier,
        nmr_threshold = nmr_threshold
      )
      classify_tgx_ddi_trial(
        prepared_trial = prepared,
        classifier = classifier,
        probability_threshold = probability_threshold,
        pc1_threshold = pc1_threshold
      )
    },
    cores = cores
  )
  names(out) <- trial_names
  out
}

summarise_tgx_ddi_results <- function(results) {
  if (!is.list(results) || length(results) == 0L) {
    stop("TGx-DDI results must be a non-empty list.")
  }

  rows <- do.call(rbind, lapply(results, function(entry) {
    data.frame(
      trial_res = as.character(entry$trial_res),
      probG_res = as.numeric(entry$probG_res),
      probNG_res = as.numeric(entry$probNG_res),
      pc1_res = as.numeric(entry$pc1_res),
      hca = as.character(entry$hca),
      n_genes_used = as.integer(entry$n_genes_used),
      positive_probability = as.logical(entry$positive_probability),
      positive_pc1 = as.logical(entry$positive_pc1),
      stringsAsFactors = FALSE
    )
  }))

  index <- parse_trial_index(rows$trial_res)
  rows$cell_line <- index$cell_line
  rows$n_per_group <- index$n_per_group
  rows$n_group_b <- index$n_group_b
  rows$trial_id <- index$trial_id
  rows$cell_line <- factor(rows$cell_line, levels = cell_line_levels)

  tidyr::pivot_longer(
    rows,
    cols = c("probG_res", "pc1_res"),
    names_to = "tgx_ddi_test",
    values_to = "tgx_ddi_value"
  )
}
