source("R/config.R")
source("R/simulation_utils.R")
source("R/tpod.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/tpod/summarize_fpr.R --export TSV --annotation CSV\n",
      "  --cell-line NAME --design DESIGN --output CSV\n",
      "  [--pathways TSV] [--expected-trials 100]\n", sep = "")
  quit(save = "no")
}

export_path <- cli_arg(args, "--export", "")
annotation_path <- cli_arg(args, "--annotation", "")
output_path <- cli_arg(args, "--output", "")
cell_line <- cli_arg(args, "--cell-line", "")
design <- cli_arg(args, "--design", "")
pathway_path <- cli_arg(args, "--pathways")
expected_trials <- parse_int(cli_arg(args, "--expected-trials"), 100L)

if (any(!nzchar(c(export_path, annotation_path, output_path))) ||
    !cell_line %in% cell_line_levels || !design %in% names(tpod_designs())) {
  stop("Supply --export, --annotation, --output, a valid --cell-line and --design. See --help.")
}
if (file.exists(output_path)) stop("Refusing to overwrite: ", output_path)
if (expected_trials < 1L) stop("--expected-trials must be positive.")

export <- read.delim(export_path, check.names = FALSE, colClasses = "character",
                     comment.char = "")
annotation <- read.csv(annotation_path, check.names = FALSE, colClasses = "character")
pathways <- if (!is.null(pathway_path)) {
  read.delim(pathway_path, colClasses = "character")
}
summary <- summarize_tpod_fpr(export, annotation, pathways, expected_trials)
summary$Design <- design
# Match the spelling used by the bundled summary and figure code.
summary$CellLine <- if (cell_line == "Liver spheroids") "Liver Spheroids" else cell_line
summary <- summary[c(setdiff(names(summary), "Param"), "Param")]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, output_path, row.names = FALSE)
message("Wrote ", output_path)
