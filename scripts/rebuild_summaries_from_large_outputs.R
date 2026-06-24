script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/summaries.R")

large_root <- Sys.getenv("FDR_LARGE_ROOT", unset = fdr_path("data", "large"))

deg_dir <- file.path(large_root, "DEGs_simulations")
gsea_dir <- Sys.getenv(
  "FDR_GSEA_DIR",
  unset = file.path(large_root, "GSEA_simulations")
)
trial_metadata_path <- Sys.getenv(
  "FDR_TRIAL_METADATA",
  unset = derived_files$trial_metadata
)
ora_dir <- file.path(large_root, "ORA_simulations", "2023-11-09")

load_named_cell_line_files <- function(paths) {
  objects <- lapply(paths, readRDS)
  unlist(objects, recursive = FALSE)
}

deg_paths <- file.path(deg_dir, c("spheroids.RDS", "tk6.RDS", "heparg.RDS", "mcf7.RDS"))
gsea_paths <- setNames(
  file.path(gsea_dir, c("spheroids_gsea.RDS", "tk6_gsea.RDS", "heparg_gsea.RDS", "mcf7_gsea.RDS")),
  c("Liver spheroids", "TK6", "HepaRG", "MCF-7")
)
ora_paths <- file.path(ora_dir, c("spheroids_ora.RDS", "tk6_ora.RDS", "heparg_ora.RDS", "mcf7_ora.RDS"))

required_files <- c(deg_paths, gsea_paths, ora_paths, trial_metadata_path)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing large input RDS files:\n", paste(missing, collapse = "\n"))
}

deg_summary <- summarise_degs_from_trials(load_named_cell_line_files(deg_paths))
trial_metadata <- readRDS(trial_metadata_path)
gsea_summary <- summarise_gsea_files(gsea_paths, trial_metadata)
ora_summary <- summarise_ora_from_trials(load_named_cell_line_files(ora_paths))

saveRDS(deg_summary, derived_files$deg)
saveRDS(gsea_summary, derived_files$gsea)
saveRDS(ora_summary, derived_files$ora)

message("Rebuilt DEG, GSEA, and ORA compact summaries.")
