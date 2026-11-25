#!/usr/bin/env Rscript
script_file <- sub(
    "--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
    setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}
source("R/config.R")
source("R/figure_registry.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "--dry-run")) {
    stop("Usage: Rscript scripts/organize_figures.R [--dry-run]")
}
dry_run <- identical(args, "--dry-run")
plan <- organize_figure_outputs(dry_run = dry_run)
print(plan, row.names = FALSE)
if (!dry_run) {
    write_figure_manifest()
    message("Organized figures and wrote docs/figure_manifest.tsv.")
}
