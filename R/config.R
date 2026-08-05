project_root <- function() {
  env_root <- Sys.getenv("FDR_CODE_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  candidate <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    is_project_root <- file.exists(file.path(candidate, "pixi.toml")) &&
      dir.exists(file.path(candidate, "R")) &&
      dir.exists(file.path(candidate, "scripts"))
    if (is_project_root) {
      return(candidate)
    }

    parent <- dirname(candidate)
    if (identical(parent, candidate)) {
      break
    }
    candidate <- parent
  }

  stop(
    "Could not locate the repository root. Run from inside the repository or ",
    "set FDR_CODE_ROOT."
  )
}

fdr_path <- function(...) {
  file.path(project_root(), ...)
}

cell_line_levels <- c("TK6", "MCF-7", "HepaRG", "Liver spheroids")
n_trials_default <- 1000L
sample_size_levels <- 3:10
gsea_min_size_default <- 15L
gsea_max_size_default <- 500L

derived_files <- list(
  deg = fdr_path("data", "derived", "deg_monte_carlo_summary.rds"),
  gsea = fdr_path("data", "derived", "gsea_monte_carlo_summary.rds"),
  ora = fdr_path("data", "derived", "ora_monte_carlo_summary.rds"),
  tpod_fdr = fdr_path("data", "derived", "tpod_fdr_grid_summary.rds"),
  tgx_ddi = fdr_path("data", "derived", "tgx_ddi_monte_carlo_summary.rds"),
  trial_metadata = fdr_path("data", "derived", "monte_carlo_trial_metadata.rds")
)
