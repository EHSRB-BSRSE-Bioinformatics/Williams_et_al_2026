script_file <- sub(
  "--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/summaries.R")
source("R/simulation_utils.R")
source("R/tgx_ddi_trials.R")

large_root <- Sys.getenv("FDR_LARGE_ROOT", unset = fdr_path("data", "large"))

deg_dir <- Sys.getenv("FDR_DEG_DIR", unset = file.path(large_root, "DEGs_simulations"))
gsea_dir <- Sys.getenv(
  "FDR_GSEA_DIR",
  unset = file.path(large_root, "GSEA_simulations")
)
trial_metadata_path <- Sys.getenv(
  "FDR_TRIAL_METADATA",
  unset = derived_files$trial_metadata
)
ora_dir <- Sys.getenv(
  "FDR_ORA_DIR",
  unset = file.path(large_root, "ORA_simulations")
)
tgx_ddi_dir <- Sys.getenv(
  "FDR_TGX_DDI_DIR",
  unset = file.path(large_root, "TGx_DDI_simulations")
)

load_named_cell_line_files <- function(paths) {
  objects <- lapply(paths, readRDS)
  unlist(objects, recursive = FALSE)
}

cell_line_stems <- vapply(cell_line_levels, cell_line_file_stem, character(1))
names(cell_line_stems) <- cell_line_levels

deg_paths <- file.path(deg_dir, paste0(cell_line_stems, ".RDS"))
gsea_paths <- setNames(
  file.path(gsea_dir, paste0(cell_line_stems, "_gsea.RDS")),
  names(cell_line_stems)
)
ora_paths <- file.path(ora_dir, paste0(cell_line_stems, "_ora.RDS"))
tgx_ddi_paths <- setNames(
  file.path(tgx_ddi_dir, paste0(cell_line_stems, "_tgx_ddi_results.RDS")),
  names(cell_line_stems)
)

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

if (all(file.exists(tgx_ddi_paths))) {
  tgx_ddi_results <- load_named_cell_line_files(tgx_ddi_paths)
  tgx_ddi_summary <- summarise_tgx_ddi_results(tgx_ddi_results)
  saveRDS(tgx_ddi_summary, derived_files$tgx_ddi)
  message("Rebuilt DEG, GSEA, ORA, and TGx-DDI compact summaries.")
} else {
  missing_tgx <- tgx_ddi_paths[!file.exists(tgx_ddi_paths)]
  message(
    "Rebuilt DEG, GSEA, and ORA compact summaries. TGx-DDI large results ",
    "were not all present, so the existing compact TGx-DDI summary was left ",
    "unchanged. Missing TGx-DDI files:\n",
    paste(missing_tgx, collapse = "\n")
  )
}
