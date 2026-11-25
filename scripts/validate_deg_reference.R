# Bounded regression check using full gene matrices and selected frozen trials.
script_file <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/simulation_utils.R")
source("R/deg_trials.R")
args <- commandArgs(trailingOnly = TRUE)
reference_path <- cli_arg(args, "--reference", "data/reference/deg_reference.rds")
output_dir <- cli_arg(args, "--output-dir", "data/large/deg_validation")
original_path <- cli_arg(args, "--original-source", "")
strict_numerics <- as_logical_flag(cli_arg(args, "--strict-numerics", "true"))
results_dir <- cli_arg(args, "--results-dir", "")
tolerance <- parse_number(cli_arg(args, "--tolerance", "1e-6"), 1e-6)
require_deg_packages()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
reference <- readRDS(reference_path)
metadata <- readRDS(derived_files$trial_metadata)
counts <- if (!nzchar(results_dir) || nzchar(original_path)) {
  read_count_table("data/input/count_table.tsv")
} else NULL
cached <- NULL
if (nzchar(results_dir)) {
  files <- list.files(results_dir, pattern = "[.]RDS$", full.names = TRUE)
  cached <- unlist(lapply(files, readRDS), recursive = FALSE)
  if (anyDuplicated(names(cached)) || !all(names(reference) %in% names(cached))) {
    stop("--results-dir must contain every reference trial exactly once.")
  }
  for (name in names(reference)) {
    if (!isTRUE(all.equal(cached[[name]]$params, make_deg_params()))) {
      stop("Cached trial parameters differ from study defaults: ", name)
    }
  }
}

# Optionally evaluate ONLY the original function definition. Never source its
# top-level loops, hard-coded paths, or 32,000-trial execution statements.
original <- NULL
if (nzchar(original_path)) {
  expressions <- parse(original_path)
  definitions <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name("js_fun_everything")), as.list(expressions))
  stopifnot(length(definitions) == 1L)
  original_env <- new.env(parent = globalenv())
  for (package in c("DESeq2", "edgeR")) {
    suppressPackageStartupMessages(library(package, character.only = TRUE))
  }
  # Original groups are exclusively A/B; mixedsort and sort agree for these.
  original_env$mixedsort <- sort
  original_env$js_predictor_variable <- "group"
  original_env$res.design <- "group"
  eval(definitions[[1L]], envir = original_env)
  original <- original_env$js_fun_everything
}

records <- list()
add_record <- function(trial, comparison, field, passed, detail) {
  records[[length(records) + 1L]] <<- data.frame(
    trial_name = trial, comparison = comparison, field = field,
    passed = passed, detail = detail
  )
}
compare_table <- function(actual, expected, trial, comparison, field) {
  actual <- as.data.frame(actual)
  expected <- as.data.frame(expected)
  # Assert membership separately; compare numbers after matching feature IDs.
  same_ids <- setequal(rownames(actual), rownames(expected))
  same_columns <- identical(names(actual), names(expected))
  add_record(trial, comparison, paste0(field, ":feature_ids"), same_ids,
             paste(nrow(actual), "observed;", nrow(expected), "expected"))
  add_record(trial, comparison, paste0(field, ":columns"), same_columns,
             paste(names(actual), collapse = ","))
  if (!same_ids || !same_columns || nrow(expected) == 0L) return(invisible(NULL))
  actual <- actual[match(rownames(expected), rownames(actual)), , drop = FALSE]
  for (column in names(expected)) {
    a <- actual[[column]]
    b <- expected[[column]]
    if (is.numeric(a) && is.numeric(b)) {
      nonfinite <- !is.finite(a) | !is.finite(b)
      same_missing <- identical(is.na(a), is.na(b)) &&
        identical(a[nonfinite], b[nonfinite])
      finite <- !nonfinite
      difference <- if (any(finite)) max(abs(a[finite] - b[finite])) else 0
      relative <- if (any(finite)) max(abs(a[finite] - b[finite]) /
                                       pmax(1, abs(b[finite]))) else 0
      passed <- same_missing && relative <= tolerance
      detail <- paste0("max_absolute_difference=", signif(difference, 8),
                       "; max_scaled_difference=", signif(relative, 8))
    } else {
      passed <- identical(a, b)
      detail <- "exact nonnumeric comparison"
    }
    add_record(trial, comparison, paste0(field, ":", column), passed, detail)
  }
}
compare_trial <- function(actual, expected, name, comparison) {
  for (field in c("DESeqRes", "DESeqResFiltered", "DESeqResDEGs",
                  "DESeqfiltered_table")) {
    compare_table(actual[[field]], expected[[field]], name, comparison, field)
  }
  add_record(name, comparison, "sample_ids",
             identical(as.character(actual$exp_metadata$sample_id),
                       as.character(expected$exp_metadata$sample_id)), "ordered sample IDs")
  add_record(name, comparison, "group_labels",
             identical(as.character(actual$exp_metadata$group),
                       as.character(expected$exp_metadata$group)), "ordered group labels")
  add_record(name, comparison, "num_rows_metadata",
             identical(actual$num_rows_metadata, expected$num_rows_metadata), "retained sample count")
  a <- actual$DESeqjs_filtered_out
  b <- expected$DESeqjs_filtered_out
  add_record(name, comparison, "filtered_out:length", length(a) == length(b),
             paste(length(a), "observed;", length(b), "expected"))
  for (j in seq_len(min(length(a), length(b)))) {
    compare_table(a[[j]], b[[j]], name, comparison, paste0("filtered_out_", j))
  }
}

# Compare feature IDs before values. The original-function comparison is stricter
# than the historical archive comparison: no numerical drift is expected when
# both functions run with the same installed dependencies.
for (name in names(reference)) {
  message("Validating ", name)
  clean <- if (is.null(cached)) run_deg_trial(name, metadata, counts) else cached[[name]]
  compare_trial(clean, reference[[name]], name, "archive")
  if (!is.null(original)) {
    legacy_params <- make_deg_params()
    legacy_params$MinCount <- legacy_params$min_cpm
    legacy_params$filter_gene_counts <- FALSE
    legacy_params$parallel <- FALSE
    legacy <- original(
      i = 1L, ls_metadata = metadata[name], count_data = counts,
      params = legacy_params, intgroup = "group",
      contrasts = data.frame(experimental = "A", control = "B"),
      js_model_formula = ~ group, ls_run_name = name
    )
    archive_tolerance <- tolerance
    tolerance <- 0
    compare_trial(clean, legacy, name, "original_function")
    tolerance <- archive_tolerance
  }
  rm(clean)
  gc()
}
report <- do.call(rbind, records)
# A relaxed archive check still requires exact feature/call membership and all
# structural/metadata checks. Only archive numeric-column differences may be
# informational. Original-function comparisons always require exact numerics.
report$required <- strict_numerics | report$comparison == "original_function" |
  startsWith(report$field, "DESeqfiltered_table:") |
  !startsWith(report$detail, "max_absolute_difference=")
write.table(report, file.path(output_dir, "comparisons.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
writeLines(c(paste("Validation time (UTC):", format(Sys.time(), tz = "UTC", usetz = TRUE)),
             paste("Scaled numerical tolerance:", tolerance),
             paste("Strict archive numerics:", strict_numerics),
             paste("Cached result directory:", if (nzchar(results_dir)) results_dir else "none; refitted"),
             "Parameters:", capture.output(dput(make_deg_params())),
             capture.output(sessionInfo())), file.path(output_dir, "session_info.txt"))
summary <- aggregate(passed ~ comparison, report,
                     function(x) paste(sum(x), "of", length(x), "checks passed"))
print(summary, row.names = FALSE)
if (any(!report$passed & report$required)) {
  print(report[!report$passed, ], row.names = FALSE)
  stop("Reference differences found; see ", output_dir, "/comparisons.tsv")
}
if (any(!report$passed)) {
  message(sum(!report$passed), " archive numeric-column checks exceeded tolerance; ",
          "all required feature/call, filter-count, and metadata checks passed.")
} else {
  message("All selected-trial comparisons passed.")
}
message("Report: ", output_dir)
