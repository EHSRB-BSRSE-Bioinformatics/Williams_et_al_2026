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
source("R/deg_trials.R")

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
output_dir <- cli_arg(args, "--output-dir", Sys.getenv(
  "FDR_DEG_DIR",
  unset = fdr_path("data", "large", "DEGs_simulations")
))

params <- make_deg_params(
  alpha = parse_number(Sys.getenv("FDR_DEG_ALPHA", unset = "0.05"), 0.05),
  linear_fc_filter = parse_number(Sys.getenv("FDR_DEG_LINEAR_FC", unset = "1.5"), 1.5),
  nmr_threshold = parse_number(Sys.getenv("FDR_NMR_THRESHOLD", unset = "100000"), 100000),
  cooks = as_logical_flag(Sys.getenv("FDR_DEG_COOKS", unset = "false"))
)

message("Reading count table: ", count_path)
count_data <- read_count_table(count_path)
message("Reading trial metadata: ", trial_metadata_path)
trial_metadata <- readRDS(trial_metadata_path)
validate_trial_metadata(trial_metadata)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (cell_line in cell_lines) {
  message("Running DEG trials for ", cell_line)
  cell_trials <- subset_trial_metadata(trial_metadata, cell_line)
  results <- run_deg_trials(
    trial_metadata = cell_trials,
    count_data = count_data,
    params = params,
    cores = cores
  )
  out_path <- file.path(output_dir, paste0(cell_line_file_stem(cell_line), ".RDS"))
  write_rds_checked(results, out_path, overwrite = overwrite)
  rm(results)
  gc()
}

message("DEG trial workflow complete. Outputs are in ", output_dir)
