script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")

issues <- character()
add_issue <- function(...) {
  issues <<- c(issues, paste0(...))
}

required_files <- c(
  derived_files$deg,
  derived_files$gsea,
  derived_files$ora,
  derived_files$tgx_ddi,
  derived_files$trial_metadata,
  fdr_path("data", "input", "count_table.tsv"),
  fdr_path("data", "input", "metadata.QC_applied.txt"),
  fdr_path("data", "external", "wikipathways-20231010-gmt-Homo_sapiens.gmt")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  add_issue("Missing required files: ", paste(missing_files, collapse = ", "))
}

required_packages <- c("dplyr", "tidyr", "ggplot2", "patchwork", "ggridges")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  add_issue("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

read_rds_checked <- function(path, label) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    readRDS(path),
    error = function(error) {
      add_issue(label, " could not be read: ", conditionMessage(error))
      NULL
    }
  )
}

check_columns <- function(data, required, label) {
  if (is.null(data)) {
    return(invisible(NULL))
  }
  if (!is.data.frame(data)) {
    add_issue(label, " must be a data frame.")
    return(invisible(NULL))
  }
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    add_issue(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

trial_metadata <- read_rds_checked(derived_files$trial_metadata, "Trial metadata")
expected_trial_count <- length(cell_line_levels) * length(sample_size_levels) * n_trials_default
if (!is.null(trial_metadata)) {
  if (!is.list(trial_metadata) || length(trial_metadata) != expected_trial_count) {
    add_issue("Trial metadata must contain ", expected_trial_count, " trial definitions.")
  }
  if (is.null(names(trial_metadata)) || anyNA(names(trial_metadata)) ||
      any(!nzchar(names(trial_metadata))) || anyDuplicated(names(trial_metadata))) {
    add_issue("Trial metadata must have complete, unique trial names.")
  }
}

deg <- read_rds_checked(derived_files$deg, "DEG summary")
gsea <- read_rds_checked(derived_files$gsea, "GSEA summary")
ora <- read_rds_checked(derived_files$ora, "ORA summary")
tgx_ddi <- read_rds_checked(derived_files$tgx_ddi, "TGx-DDI summary")

check_columns(
  deg,
  c("cell_line", "n_per_group", "Filter", "nDEGs", "rank", "contains_false_positive"),
  "DEG summary"
)
check_columns(
  gsea,
  c(
    "cell_line", "n", "min_size", "max_size", "p_lt_0.05",
    "p_lt_0.05_ES_gt_2", "rank", "rank_filtered"
  ),
  "GSEA summary"
)
check_columns(
  ora,
  c("cell_line", "n", "go_bp", "enricher_wikipathways", "go_bp_rank", "wiki_rank"),
  "ORA summary"
)
check_columns(
  tgx_ddi,
  c("cell_line", "n_per_group", "tgx_ddi_test", "tgx_ddi_value"),
  "TGx-DDI summary"
)

if (!is.null(gsea) && is.data.frame(gsea)) {
  valid_gsea_range <- all(c("min_size", "max_size") %in% names(gsea)) &&
    any(
      gsea$min_size == gsea_min_size_default &
        gsea$max_size == gsea_max_size_default,
      na.rm = TRUE
    )
  if (!valid_gsea_range) {
    add_issue(
      "GSEA summary does not contain the publication gene-set range ",
      gsea_min_size_default, "--", gsea_max_size_default, "."
    )
  }
}

expected_rows <- c(gsea = expected_trial_count, ora = expected_trial_count, tgx_ddi = 2L * expected_trial_count)
summary_objects <- list(gsea = gsea, ora = ora, tgx_ddi = tgx_ddi)
for (label in names(expected_rows)) {
  object <- summary_objects[[label]]
  if (!is.null(object) && is.data.frame(object) && nrow(object) != expected_rows[[label]]) {
    add_issue(
      label, " summary has ", nrow(object), " rows; expected ",
      expected_rows[[label]], "."
    )
  }
}

metadata_path <- fdr_path("data", "input", "metadata.QC_applied.txt")
count_path <- fdr_path("data", "input", "count_table.tsv")
if (file.exists(metadata_path) && file.exists(count_path)) {
  metadata <- tryCatch(
    read.delim(metadata_path, check.names = FALSE),
    error = function(error) {
      add_issue("Sample metadata could not be read: ", conditionMessage(error))
      NULL
    }
  )
  count_header <- strsplit(readLines(count_path, n = 1L), "\t", fixed = TRUE)[[1L]]

  if (!is.null(metadata)) {
    required_metadata_columns <- c("Sample_ID", "chemical", "cell_line")
    missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
    if (length(missing_metadata_columns) > 0L) {
      add_issue(
        "Sample metadata is missing columns: ",
        paste(missing_metadata_columns, collapse = ", ")
      )
    } else {
      if (anyDuplicated(metadata$Sample_ID)) {
        add_issue("Sample metadata contains duplicate Sample_ID values.")
      }
      missing_count_samples <- setdiff(metadata$Sample_ID, count_header)
      if (length(missing_count_samples) > 0L) {
        add_issue(
          "Count matrix is missing metadata samples: ",
          paste(missing_count_samples, collapse = ", ")
        )
      }

      dmso_counts <- table(metadata$cell_line[metadata$chemical == "DMSO"])
      insufficient <- cell_line_levels[
        !cell_line_levels %in% names(dmso_counts) |
          dmso_counts[cell_line_levels] < 2L * max(sample_size_levels)
      ]
      if (length(insufficient) > 0L) {
        add_issue(
          "Insufficient DMSO samples for the maximum trial size in: ",
          paste(insufficient, collapse = ", ")
        )
      }
    }
  }
}

if (length(issues) > 0L) {
  message("Repository check failed:")
  message(paste(" -", issues, collapse = "\n"))
  quit(status = 1L)
}

message("Repository check passed.")
