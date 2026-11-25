source("R/config.R")
source("R/simulation_utils.R")
source("R/tpod.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Usage: Rscript scripts/tpod/prepare_mock_doses.R --output-dir DIR\n",
      "  [--cell-lines LIST] [--designs LIST] [--n-trials 100] [--seed 90214]\n",
      "  [--count-table TSV] [--metadata TSV]\n",
      "Cell lines and designs are comma-separated and default to all.\n", sep = "")
  quit(save = "no")
}

output_dir <- cli_arg(args, "--output-dir")
if (is.null(output_dir) || !nzchar(output_dir)) stop("Supply --output-dir.")
if (file.exists(output_dir)) stop("Choose a new output directory: ", output_dir)

cell_lines <- unique(parse_cell_lines_arg(cli_arg(args, "--cell-lines")))
design_arg <- cli_arg(args, "--designs", "all")
designs <- if (design_arg == "all") {
  names(tpod_designs())
} else {
  trimws(strsplit(design_arg, ",", fixed = TRUE)[[1L]])
}
n_trials <- parse_int(cli_arg(args, "--n-trials"), 100L)
seed <- parse_int(cli_arg(args, "--seed"), 90214L)
count_path <- cli_arg(args, "--count-table", fdr_path("data", "input", "count_table.tsv"))
metadata_path <- cli_arg(args, "--metadata", fdr_path("data", "input", "metadata.QC_applied.txt"))

counts <- read_count_table(count_path)
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE)
prepared <- prepare_tpod_expression(counts, metadata)
assignments <- lapply(cell_lines, function(cell_line) {
  samples <- colnames(prepared$expression)[which(prepared$cell_line == cell_line)]
  make_tpod_assignments(samples, designs, n_trials, seed)
})
for (i in seq_along(cell_lines)) {
  stem <- cell_line_file_stem(cell_lines[[i]])
  for (design in designs) {
    design_dir <- file.path(output_dir, stem, design)
    write_tpod_inputs(prepared$expression, assignments[[i]][[design]], design_dir)
  }
}
message("Wrote BMDExpress inputs under ", output_dir)
