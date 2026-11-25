script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/monte_carlo_trials.R")

args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(Sys.getenv("FDR_MONTE_CARLO_SEED", unset = "1"))
output <- fdr_path("data", "derived", paste0("monte_carlo_trial_metadata_seed", seed, ".rds"))

if ("--overwrite" %in% args && file.exists(output)) {
  unlink(output)
}

if (file.exists(output)) {
  stop("Refusing to overwrite existing trial metadata: ", output)
}

metadata <- read.delim(fdr_path("data", "input", "metadata.QC_applied.txt"))
trial_metadata <- make_trial_metadata(metadata, seed = seed)
saveRDS(trial_metadata, output)
message("Wrote ", length(trial_metadata), " Monte Carlo trial definitions to ", output)
