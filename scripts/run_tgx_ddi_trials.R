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
source("R/tgx_ddi_trials.R")

args <- commandArgs(trailingOnly = TRUE)
cell_lines <- parse_cell_lines_arg(cli_arg(
  args,
  "--cell-lines",
  Sys.getenv("FDR_CELL_LINES", unset = "all")
))
cores <- parse_int(cli_arg(args, "--cores", Sys.getenv("FDR_CORES", unset = "1")), 1L)
overwrite <- as_logical_flag(cli_arg(args, "--overwrite", Sys.getenv("FDR_OVERWRITE", unset = "false")))

count_path <- cli_arg(args, "--count-table", Sys.getenv(
  "FDR_COUNT_TABLE",
  unset = fdr_path("data", "input", "count_table.tsv")
))
trial_metadata_path <- cli_arg(args, "--trial-metadata", Sys.getenv(
  "FDR_TRIAL_METADATA",
  unset = derived_files$trial_metadata
))
classifier_path <- cli_arg(args, "--classifier", Sys.getenv("FDR_TGX_DDI_CLASSIFIER", unset = ""))
if (!nzchar(classifier_path)) {
  stop("Set --classifier or FDR_TGX_DDI_CLASSIFIER to the TGx-DDI classifier TSV.")
}
output_dir <- cli_arg(args, "--output-dir", Sys.getenv(
  "FDR_TGX_DDI_DIR",
  unset = fdr_path("data", "large", "TGx_DDI_simulations")
))
summary_path <- cli_arg(args, "--summary-output", Sys.getenv(
  "FDR_TGX_DDI_SUMMARY",
  unset = derived_files$tgx_ddi
))

nmr_threshold <- parse_number(Sys.getenv("FDR_NMR_THRESHOLD", unset = "100000"), 100000)
probability_threshold <- parse_number(Sys.getenv("FDR_TGX_DDI_PROB_THRESHOLD", unset = "0.9"), 0.9)
pc1_threshold <- parse_number(Sys.getenv("FDR_TGX_DDI_PC1_THRESHOLD", unset = "0"), 0)

message("Reading count table: ", count_path)
count_data <- read_count_table(count_path)
message("Reading trial metadata: ", trial_metadata_path)
trial_metadata <- readRDS(trial_metadata_path)
validate_trial_metadata(trial_metadata)
message("Reading TGx-DDI classifier: ", classifier_path)
classifier <- read_tgx_ddi_classifier(classifier_path)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
combined <- list()
for (cell_line in cell_lines) {
  message("Running TGx-DDI trials for ", cell_line)
  cell_trials <- subset_trial_metadata(trial_metadata, cell_line)
  results <- run_tgx_ddi_trials(
    trial_metadata = cell_trials,
    count_data = count_data,
    classifier = classifier,
    nmr_threshold = nmr_threshold,
    cores = cores,
    probability_threshold = probability_threshold,
    pc1_threshold = pc1_threshold
  )
  combined <- c(combined, results)
  write_rds_checked(
    results,
    file.path(output_dir, paste0(cell_line_file_stem(cell_line), "_tgx_ddi_results.RDS")),
    overwrite = overwrite
  )
  rm(results)
  gc()
}

summary <- summarise_tgx_ddi_results(combined)
write_rds_checked(summary, summary_path, overwrite = overwrite)
message("TGx-DDI trial workflow complete. Outputs are in ", output_dir)
