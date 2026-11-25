# Shared helpers for the trial-level simulation workflow.
# These functions intentionally avoid loading heavyweight Bioconductor packages;
# analysis-specific dependencies are loaded in the DEG, pathway, and TGx-DDI
# modules.

trim_ws <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

cli_arg <- function(args, name, default = NULL) {
  exact <- which(args == name)
  if (length(exact) > 0L) {
    if (exact[[1L]] == length(args)) {
      stop("Missing value after ", name, ".")
    }
    return(args[[exact[[1L]] + 1L]])
  }

  prefix <- paste0(name, "=")
  inline <- grep(paste0("^", gsub("-", "\\\\-", prefix)), args, value = TRUE)
  if (length(inline) > 0L) {
    return(sub(prefix, "", inline[[1L]], fixed = TRUE))
  }

  default
}

as_logical_flag <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y", "overwrite")
}

parse_int <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  out <- suppressWarnings(as.integer(x))
  if (is.na(out)) {
    stop("Expected an integer, got: ", x)
  }
  out
}

parse_number <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  out <- suppressWarnings(as.numeric(x))
  if (is.na(out)) {
    stop("Expected a number, got: ", x)
  }
  out
}

parse_cell_lines_arg <- function(value, allowed = cell_line_levels) {
  if (is.null(value) || !nzchar(value) || identical(tolower(value), "all")) {
    return(allowed)
  }
  requested <- trim_ws(unlist(strsplit(value, ",", fixed = TRUE)))
  requested <- requested[nzchar(requested)]
  if (length(requested) == 0L) {
    return(allowed)
  }
  invalid <- setdiff(requested, allowed)
  if (length(invalid) > 0L) {
    stop(
      "Unknown cell line(s): ", paste(invalid, collapse = ", "),
      ". Valid values are: ", paste(allowed, collapse = ", ")
    )
  }
  requested
}

cell_line_file_stem <- function(cell_line) {
  stems <- c(
    "TK6" = "tk6",
    "MCF-7" = "mcf7",
    "HepaRG" = "heparg",
    "Liver spheroids" = "spheroids"
  )
  out <- unname(stems[cell_line])
  if (is.na(out)) {
    stop("No file stem is defined for cell line: ", cell_line)
  }
  out
}

read_count_table <- function(path) {
  if (!file.exists(path)) {
    stop("Count table does not exist: ", path)
  }
  counts <- read.delim(path, check.names = FALSE, row.names = 1)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  if (anyNA(counts)) {
    stop("Count table contains NA values: ", path)
  }
  if (any(counts < 0)) {
    stop("Count table contains negative values: ", path)
  }
  counts
}

parse_trial_index <- function(trial_names) {
  if (is.null(trial_names) || anyNA(trial_names) || any(!nzchar(trial_names))) {
    stop("Trial names must be complete, non-empty strings.")
  }

  parts <- strsplit(trial_names, "_", fixed = TRUE)
  if (any(lengths(parts) != 4L)) {
    stop("Trial names must follow '<cell_line>_<n_A>_<n_B>_<trial_id>'.")
  }

  out <- data.frame(
    run = trial_names,
    cell_line = vapply(parts, `[[`, character(1), 1L),
    n_per_group = suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 2L))),
    n_group_b = suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 3L))),
    trial_id = suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 4L))),
    stringsAsFactors = FALSE
  )

  if (anyNA(out[c("n_per_group", "n_group_b", "trial_id")])) {
    stop("Trial sample sizes and trial IDs must be integers.")
  }
  out
}

validate_trial_metadata <- function(trial_metadata) {
  if (!is.list(trial_metadata) || is.null(names(trial_metadata))) {
    stop("Trial metadata must be a named list.")
  }
  if (anyNA(names(trial_metadata)) || any(!nzchar(names(trial_metadata)))) {
    stop("Trial metadata list names must be complete.")
  }
  if (anyDuplicated(names(trial_metadata))) {
    stop("Trial metadata list names must be unique.")
  }

  required <- c("sample_id", "group")
  bad <- vapply(trial_metadata, function(x) {
    !is.data.frame(x) || !all(required %in% names(x))
  }, logical(1))
  if (any(bad)) {
    stop(
      "Every trial metadata entry must be a data frame with columns: ",
      paste(required, collapse = ", "),
      ". First invalid entry: ", names(trial_metadata)[which(bad)[[1L]]]
    )
  }

  invisible(trial_metadata)
}

subset_trial_metadata <- function(trial_metadata, cell_line = NULL) {
  validate_trial_metadata(trial_metadata)
  if (is.null(cell_line)) {
    return(trial_metadata)
  }

  index <- parse_trial_index(names(trial_metadata))
  keep <- index$cell_line == cell_line
  if (!any(keep)) {
    stop("No trials found for cell line: ", cell_line)
  }
  trial_metadata[keep]
}

prepare_trial_counts <- function(trial,
                                 count_data,
                                 nmr_threshold = 100000,
                                 group_levels = c("A", "B")) {
  if (!all(c("sample_id", "group") %in% names(trial))) {
    stop("Trial metadata must contain sample_id and group columns.")
  }

  trial$sample_id <- as.character(trial$sample_id)
  trial$group <- as.character(trial$group)

  total_reads <- colSums(count_data)
  count_data <- count_data[, total_reads > nmr_threshold, drop = FALSE]

  missing_samples <- setdiff(trial$sample_id, colnames(count_data))
  if (length(missing_samples) > 0L) {
    stop(
      "Trial contains samples that are absent from the count table or fail the ",
      "read-depth threshold: ", paste(missing_samples, collapse = ", ")
    )
  }

  count_subset <- count_data[, trial$sample_id, drop = FALSE]
  trial$group <- factor(as.character(trial$group), levels = group_levels)
  if (anyNA(trial$group)) {
    stop("Trial contains group labels outside expected levels: ", paste(group_levels, collapse = ", "))
  }
  rownames(trial) <- trial$sample_id

  list(count_data = count_subset, exp_metadata = trial)
}

run_trials_parallel <- function(x, fun, cores = 1L, ..., mc_preschedule = TRUE) {
  cores <- max(1L, as.integer(cores))
  if (cores > 1L && .Platform$OS.type == "unix") {
    return(parallel::mclapply(
      X = x,
      FUN = fun,
      ...,
      mc.cores = cores,
      mc.preschedule = mc_preschedule
    ))
  }
  if (cores > 1L && .Platform$OS.type != "unix") {
    warning("Parallel fork execution is unavailable on this platform; using serial lapply().")
  }
  lapply(x, fun, ...)
}

write_rds_checked <- function(object, path, overwrite = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !overwrite) {
    stop("Refusing to overwrite existing file: ", path)
  }
  saveRDS(object, path)
  message("Wrote ", path)
  invisible(path)
}
