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
required_packages <- c("dplyr", "tidyr", "ggplot2", "patchwork", "ggridges")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages required to render figures: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before rerunning scripts/render_figures.R."
  )
}

source("R/summaries.R")
source("R/figures.R")
source("R/recurrence_figures.R")

dir.create(fdr_path("figures"), showWarnings = FALSE, recursive = TRUE)

figure_specs <- list(
  s5_recurrence = list(
    plot = function() build_s5_recurrence_figure(summaries$s56),
    pdf_size = c(14, 8),
    png_size = c(14, 8)
  ),
  s6_expression = list(
    plot = function() build_s6_expression_figure(summaries$s56),
    pdf_size = c(16, 12),
    png_size = c(16, 12)
  ),
  deg = list(
    plot = function() build_deg_figure(summaries$deg),
    pdf_size = c(12, 7),
    png_size = c(10.5, 6)
  ),
  ora = list(
    plot = function() build_ora_figure(summaries$ora),
    pdf_size = c(18, 7),
    png_size = c(16, 6)
  ),
  ora_wikipathways_supplementary = list(
    plot = function() {
      build_ora_figure(
        summaries$ora,
        collection = "WikiPathways",
        include_go = FALSE
      )
    },
    pdf_size = c(6, 7),
    png_size = c(5.25, 6)
  ),
  gsea = list(
    plot = function() build_gsea_figure(summaries$gsea),
    pdf_size = c(12, 7),
    png_size = c(10.5, 6)
  ),
  gsea_archived_supplementary = list(
    plot = function() build_gsea_figure(summaries$gsea_archived),
    pdf_size = c(12, 7),
    png_size = c(10.5, 6)
  ),
  tgx_ddi = list(
    plot = function() build_tgx_ddi_figure(summaries$tgx_ddi),
    pdf_size = c(12, 8),
    png_size = c(10.5, 7)
  ),
  tpod_fdr = list(
    plot = function() build_tpod_fdr_figure(summaries$tpod_fdr),
    pdf_size = c(24, 13),
    png_size = c(24, 13)
  ),
  tpod_fdr_main = list(
    plot = function() build_tpod_fdr_main_figure(summaries$tpod_fdr),
    pdf_size = c(12, 17.5),
    png_size = c(12, 17.5)
  ),
  tpod_robustness = list(
    plot = function() build_tpod_robustness_figure(summaries$tpod_fdr),
    pdf_size = c(11, 10),
    png_size = c(11, 10)
  ),
  tpod_parameter_settings = list(
    plot = function() build_tpod_parameter_setting_figure(summaries$tpod_fdr),
    pdf_size = c(12, 10),
    png_size = c(12, 10)
  ),
  tpod_fdr_supplementary = list(
    plot = function() build_tpod_fdr_supplementary_figure(summaries$tpod_fdr),
    pdf_size = c(24, 13),
    png_size = c(24, 13)
  ),
  s8_mean_variance = list(
    plot = function() {
      build_s8_mean_variance_figure(summaries$s8_mean_variance)
    },
    pdf_size = c(11, 5.5),
    png_size = c(11, 5.5)
  )
)

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--publication")) {
  registry <- figure_registry()
  publication_keys <- registry$key[registry$category %in% c("main", "supplementary")]
  figure_specs <- figure_specs[names(figure_specs) %in% publication_keys]
} else if (length(args)) {
  if (length(args) != 1L || !startsWith(args, "--only=")) {
    stop("Usage: Rscript scripts/render_figures.R [--publication | --only=figure1,figure2]")
  }
  selected <- strsplit(sub("^--only=", "", args), ",", fixed = TRUE)[[1]]
  if (!length(selected) || any(!selected %in% names(figure_specs))) {
    stop(
      "Unknown figure selection. Available: ",
      paste(names(figure_specs), collapse = ", ")
    )
  }
  figure_specs <- figure_specs[unique(selected)]
}

summaries <- if (
  all(names(figure_specs) %in% c("s5_recurrence", "s6_expression"))
) {
  list(s56 = read_s56_summary())
} else if (identical(names(figure_specs), "s8_mean_variance")) {
  list(s8_mean_variance = read_s8_figure_summary())
} else if (identical(names(figure_specs), "deg")) {
  list(deg = read_deg_figure_summary())
} else {
  read_summary_inputs()
}
if (
  any(c("s5_recurrence", "s6_expression") %in% names(figure_specs)) &&
    !"s56" %in% names(summaries)
) {
  summaries$s56 <- read_s56_summary()
}
if (
  "s8_mean_variance" %in%
    names(figure_specs) &&
    !"s8_mean_variance" %in% names(summaries)
) {
  summaries$s8_mean_variance <- read_s8_figure_summary()
}
set.seed(1L)

if ("deg" %in% names(figure_specs)) {
  saveRDS(
    summaries$deg,
    fdr_path("data", "derived", "deg_figure_summary.rds"),
    compress = "xz"
  )
  deg_rates <- false_positive_rates(
    summaries$deg,
    "contains_false_positive",
    "n_per_group",
    c("cell_line", "Filter")
  )
  deg_rates$FDR <- NULL
  write.table(
    deg_rates,
    fdr_path("data", "derived", "deg_figure_rates.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

for (name in names(figure_specs)) {
  spec <- figure_specs[[name]]
  plot <- spec$plot()
  print(paste0("Writing ", name))
  ggplot2::ggsave(
    filename = figure_path(name, "pdf"),
    plot = plot,
    width = spec$pdf_size[[1L]],
    height = spec$pdf_size[[2L]],
    units = "in",
    device = cairo_pdf
  )
  ggplot2::ggsave(
    filename = figure_path(name, "png"),
    plot = plot,
    width = spec$png_size[[1L]],
    height = spec$png_size[[2L]],
    units = "in",
    dpi = 600,
    device = png
  )
}

if ("deg" %in% names(figure_specs)) {
  sources <- c(
    fdr_path("data", "derived", "deg_thresholds", "threshold_sensitivity.rds"),
    derived_files$trial_metadata,
    derived_files$deg,
    fdr_path("R", "summaries.R"),
    fdr_path("R", "figures.R"),
    fdr_path("R", "figure_registry.R"),
    fdr_path("scripts", "render_figures.R"),
    fdr_path("tests", "test_deg_figure.R"),
    fdr_path("scripts", "check_analysis_consistency.R")
  )
  outputs <- c(
    figure_path("deg", c("pdf", "png")),
    fdr_path("data", "derived", "deg_figure_summary.rds"),
    fdr_path("data", "derived", "deg_figure_rates.tsv")
  )
  paths <- c(sources, outputs)
  if (!nzchar(Sys.which("sha256sum"))) {
    stop("sha256sum is required for figure provenance.")
  }
  hashes <- vapply(
    paths,
    function(path) {
      strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " +")[[1L]][
        1L
      ]
    },
    character(1)
  )
  write.table(
    data.frame(
      role = c(rep("input", length(sources)), rep("output", length(outputs))),
      path = substring(paths, nchar(project_root()) + 2L),
      sha256 = unname(hashes)
    ),
    fdr_path("data", "derived", "deg_figure_provenance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if ("s8_mean_variance" %in% names(figure_specs)) {
  sources <- c(
    fdr_path("data", "input", "count_table.tsv"),
    fdr_path("data", "input", "metadata.QC_applied.txt"),
    fdr_path("data", "derived", "s8_mean_variance_summary.rds"),
    fdr_path("R", "config.R"),
    fdr_path("R", "summaries.R"),
    fdr_path("R", "figures.R"),
    fdr_path("R", "figure_registry.R"),
    fdr_path("scripts", "rebuild_s8_figure_data.R"),
    fdr_path("scripts", "render_figures.R")
  )
  outputs <- c(
    figure_path("s8_mean_variance", c("pdf", "png"))
  )
  paths <- c(sources, outputs)
  if (!nzchar(Sys.which("sha256sum"))) {
    stop("sha256sum is required for Figure S8 provenance.")
  }
  hashes <- vapply(
    paths,
    function(path) {
      strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " +")[[1L]][
        1L
      ]
    },
    character(1)
  )
  write.table(
    data.frame(
      role = c(rep("input", length(sources)), rep("output", length(outputs))),
      path = substring(paths, nchar(project_root()) + 2L),
      sha256 = unname(hashes)
    ),
    fdr_path("data", "derived", "s8_figure_provenance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if (any(c("s5_recurrence", "s6_expression") %in% names(figure_specs))) {
  paths <- c(
    "data/derived/s56_figure_summary.rds",
    "R/config.R",
    "R/recurrence_figures.R",
    "R/figure_registry.R",
    "scripts/render_figures.R"
  )
  rendered <- intersect(
    c("s5_recurrence", "s6_expression"),
    names(figure_specs)
  )
  outputs <- unlist(
    lapply(rendered, function(name) {
      figure_filename(name, c("pdf", "png"))
    }),
    use.names = FALSE
  )
  s56_write_tsv(
    data.frame(
      role = c(
        rep("input_or_code", length(paths)),
        rep("output", length(outputs))
      ),
      path = c(paths, outputs),
      sha256 = vapply(c(paths, outputs), s56_sha256, "")
    ),
    "data/derived/s56_figure_provenance.tsv"
  )
}
message("Wrote figures to ", fdr_path("figures"))
