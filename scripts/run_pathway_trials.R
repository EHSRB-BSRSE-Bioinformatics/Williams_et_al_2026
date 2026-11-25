script_file <- sub(
  "--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/simulation_utils.R")
source("R/deg_trials.R")
source("R/pathway_trials.R")

args <- commandArgs(trailingOnly = TRUE)
analysis <- cli_arg(args, "--analysis", Sys.getenv("FDR_PATHWAY_ANALYSIS", unset = "both"))
analysis <- tolower(analysis)
if (!analysis %in% c("ora", "gsea", "both")) {
  stop("--analysis must be one of: ora, gsea, both")
}
analyses <- if (identical(analysis, "both")) c("ora", "gsea") else analysis

cell_lines <- parse_cell_lines_arg(cli_arg(
  args,
  "--cell-lines",
  Sys.getenv("FDR_CELL_LINES", unset = "all")
))
cores <- parse_int(cli_arg(args, "--cores", Sys.getenv("FDR_CORES", unset = "1")), 1L)
overwrite <- as_logical_flag(cli_arg(args, "--overwrite", Sys.getenv("FDR_OVERWRITE", unset = "false")))

deg_dir <- cli_arg(args, "--deg-dir", Sys.getenv(
  "FDR_DEG_DIR",
  unset = fdr_path("data", "large", "DEGs_simulations")
))
gsea_dir <- cli_arg(args, "--gsea-dir", Sys.getenv(
  "FDR_GSEA_DIR",
  unset = fdr_path("data", "large", "GSEA_simulations")
))
ora_dir <- cli_arg(args, "--ora-dir", Sys.getenv(
  "FDR_ORA_DIR",
  unset = fdr_path("data", "large", "ORA_simulations", as.character(Sys.Date()))
))
wikipathways_path <- cli_arg(args, "--wikipathways", Sys.getenv(
  "FDR_WIKIPATHWAYS_GMT",
  unset = fdr_path("data", "external", "wikipathways-20231010-gmt-Homo_sapiens.gmt")
))

message("Loading human gene alias annotation")
annotation <- load_human_alias_annotation()

wikipathways <- NULL
if ("ora" %in% analyses) {
  message("Loading WikiPathways GMT: ", wikipathways_path)
  wikipathways <- read_wikipathways_gmt(wikipathways_path)
}

reactome_pathways <- NULL
if ("gsea" %in% analyses) {
  message("Loading Reactome pathways")
  reactome_pathways <- load_reactome_pathways()
}

for (cell_line in cell_lines) {
  stem <- cell_line_file_stem(cell_line)
  deg_path <- file.path(deg_dir, paste0(stem, ".RDS"))
  if (!file.exists(deg_path)) {
    stop("Missing DEG trial file for ", cell_line, ": ", deg_path)
  }

  message("Reading DEG trials for ", cell_line, ": ", deg_path)
  deg_trials <- readRDS(deg_path)

  if ("ora" %in% analyses) {
    message("Running ORA trials for ", cell_line)
    ora_results <- run_pathway_trials(
      deg_trials = deg_trials,
      analysis = "ora",
      annotation = annotation,
      wikipathways = wikipathways,
      cores = cores
    )
    write_rds_checked(
      ora_results,
      file.path(ora_dir, paste0(stem, "_ora.RDS")),
      overwrite = overwrite
    )
    rm(ora_results)
    gc()
  }

  if ("gsea" %in% analyses) {
    message("Running GSEA trials for ", cell_line)
    gsea_results <- run_pathway_trials(
      deg_trials = deg_trials,
      analysis = "gsea",
      annotation = annotation,
      reactome_pathways = reactome_pathways,
      cores = cores,
      min_size = 1L,
      max_size = gsea_max_size_default,
      n_perm_simple = parse_int(Sys.getenv("FDR_GSEA_N_PERM_SIMPLE", unset = "500"), 500L)
    )
    write_rds_checked(
      gsea_results,
      file.path(gsea_dir, paste0(stem, "_gsea.RDS")),
      overwrite = overwrite
    )
    rm(gsea_results)
    gc()
  }

  rm(deg_trials)
  gc()
}

message("Pathway trial workflow complete.")
