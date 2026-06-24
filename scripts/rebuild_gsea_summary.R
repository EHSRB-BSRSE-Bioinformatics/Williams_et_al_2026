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

large_root <- Sys.getenv("FDR_LARGE_ROOT", unset = fdr_path("data", "large"))
gsea_dir <- Sys.getenv(
  "FDR_GSEA_DIR",
  unset = file.path(large_root, "GSEA_simulations")
)
trial_metadata_path <- Sys.getenv(
  "FDR_TRIAL_METADATA",
  unset = derived_files$trial_metadata
)
gsea_paths <- setNames(
  file.path(
    gsea_dir,
    c(
      "spheroids_gsea.RDS",
      "tk6_gsea.RDS",
      "heparg_gsea.RDS",
      "mcf7_gsea.RDS"
    )
  ),
  c("Liver spheroids", "TK6", "HepaRG", "MCF-7")
)

required_files <- c(gsea_paths, trial_metadata_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing GSEA rebuild inputs:\n", paste(missing_files, collapse = "\n"))
}

trial_metadata <- readRDS(trial_metadata_path)
gsea_summary <- summarise_gsea_files(gsea_paths, trial_metadata)

dir.create(dirname(derived_files$gsea), recursive = TRUE, showWarnings = FALSE)
saveRDS(gsea_summary, derived_files$gsea)
message(
  "Wrote ", nrow(gsea_summary), " GSEA trial summaries with min_size = ",
  gsea_min_size_default, " to ", derived_files$gsea
)
