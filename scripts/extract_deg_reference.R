# Extract a small, unchanged reference set from the archived per-cell-line RDS
# files. Reads one large file at a time; never fits a model or rewrites archives.
script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/simulation_utils.R")
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Use pixi run -e deg.")
args <- commandArgs(trailingOnly = TRUE)
archive_dir <- cli_arg(args, "--archive-dir", "")
if (!nzchar(archive_dir)) stop("Supply --archive-dir containing the four DEG RDS files.")
output_dir <- cli_arg(args, "--output-dir", "data/reference")
if (any(file.exists(file.path(output_dir, c("deg_reference.rds",
    "deg_reference_metadata.rds", "deg_reference_trials.tsv"))))) {
  stop("Reference files already exist; use a fresh --output-dir to re-extract.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
metadata <- readRDS(derived_files$trial_metadata)
reference <- list()
selection <- list()
versions <- list()
for (cell_line in cell_line_levels) {
  archive_path <- file.path(archive_dir, paste0(cell_line_file_stem(cell_line), ".RDS"))
  message("Reading archive: ", archive_path)
  archive <- readRDS(archive_path)
  index <- parse_trial_index(names(archive))
  selected <- c(which(index$n_per_group == 3L)[1L],
                which(index$n_per_group == 10L)[1L])
  reasons <- c("first n=3 trial", "first n=10 trial")
  # Exercise significant calls, quantile rejection, and spike rejection when
  # present in this archive. Selection is deterministic and uses no refitted data.
  targets <- list(
    "nonzero final DEGs" = function(x) nrow(x$DESeqResDEGs) > 0L,
    "quantile rejection" = function(x) nrow(x$DESeqfiltered_table) > 0L &&
      x$DESeqfiltered_table$quantile_filtered[1L] > 0L,
    "spike rejection" = function(x) nrow(x$DESeqfiltered_table) > 0L &&
      x$DESeqfiltered_table$spike_filtered[1L] > 0L
  )
  for (reason in names(targets)) {
    matches <- which(vapply(archive, targets[[reason]], logical(1)))
    if (length(matches) && !matches[1L] %in% selected) {
      selected <- c(selected, matches[1L])
      reasons <- c(reasons, reason)
    }
  }
  for (j in seq_along(selected)) {
    name <- names(archive)[selected[j]]
    trial <- archive[[selected[j]]]
    stopifnot(identical(as.character(trial$exp_metadata$sample_id),
                        as.character(metadata[[name]]$sample_id)),
              identical(as.character(trial$exp_metadata$group),
                        as.character(metadata[[name]]$group)))
    reference[[name]] <- list(
      DESeqRes = as.data.frame(trial$DESeqRes),
      DESeqResFiltered = as.data.frame(trial$DESeqResFiltered),
      DESeqResDEGs = as.data.frame(trial$DESeqResDEGs),
      DESeqfiltered_table = trial$DESeqfiltered_table,
      DESeqjs_filtered_out = lapply(trial$DESeqjs_filtered_out, as.data.frame),
      exp_metadata = trial$exp_metadata,
      num_rows_metadata = trial$num_rows_metadata
    )
    selection[[name]] <- data.frame(
      trial_name = name, reason = reasons[j],
      archive_file = basename(archive_path),
      tested_features = nrow(trial$DESeqRes),
      significant_features = sum(trial$DESeqRes$padj < 0.05, na.rm = TRUE),
      final_degs = nrow(trial$DESeqResDEGs)
    )
    versions[[name]] <- list(metadata = S4Vectors::metadata(trial$DESeqRes),
                             priorInfo = trial$DESeqRes@priorInfo)
  }
  rm(archive)
  gc()
}
saveRDS(reference, file.path(output_dir, "deg_reference.rds"), compress = "xz")
saveRDS(versions, file.path(output_dir, "deg_reference_metadata.rds"), compress = "xz")
write.table(do.call(rbind, selection), file.path(output_dir, "deg_reference_trials.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Extracted ", length(reference), " reference trials.")
