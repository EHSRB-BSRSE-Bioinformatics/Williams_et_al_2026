script_file <- sub(
  "--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/simulation_utils.R")

issues <- character()
add_issue <- function(...) {
  issues <<- c(issues, paste0(...))
}

required_packages <- c(
  "dplyr",
  "tidyr",
  "data.table",
  "DESeq2",
  "edgeR",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "reactome.db",
  "clusterProfiler",
  "fgsea"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  add_issue("Missing simulation packages: ", paste(missing_packages, collapse = ", "))
}

required_files <- c(
  derived_files$trial_metadata,
  fdr_path("data", "input", "count_table.tsv"),
  fdr_path("data", "input", "metadata.QC_applied.txt"),
  fdr_path("data", "external", "wikipathways-20231010-gmt-Homo_sapiens.gmt")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  add_issue("Missing required simulation files: ", paste(missing_files, collapse = ", "))
}

classifier_path <- Sys.getenv("FDR_TGX_DDI_CLASSIFIER", unset = "")
if (!nzchar(classifier_path)) {
  add_issue(
    "FDR_TGX_DDI_CLASSIFIER is not set. Set it to the TGx-DDI classifier TSV ",
    "before running scripts/run_tgx_ddi_trials.R."
  )
} else if (!file.exists(classifier_path)) {
  add_issue("TGx-DDI classifier file does not exist: ", classifier_path)
}

if (file.exists(derived_files$trial_metadata)) {
  trial_metadata <- tryCatch(
    readRDS(derived_files$trial_metadata),
    error = function(error) {
      add_issue("Trial metadata could not be read: ", conditionMessage(error))
      NULL
    }
  )
  if (!is.null(trial_metadata)) {
    tryCatch(
      validate_trial_metadata(trial_metadata),
      error = function(error) add_issue(conditionMessage(error))
    )
    expected_trial_count <- length(cell_line_levels) *
      length(sample_size_levels) *
      n_trials_default
    if (length(trial_metadata) != expected_trial_count) {
      add_issue(
        "Trial metadata has ", length(trial_metadata),
        " entries; expected ", expected_trial_count, "."
      )
    }
  }
}

if (length(issues) > 0L) {
  message("Simulation setup check failed:")
  message(paste(" -", issues, collapse = "\n"))
  quit(status = 1L)
}

message("Simulation setup check passed.")
