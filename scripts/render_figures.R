script_file <- sub(
  "--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
  setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
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

dir.create(fdr_path("figures"), showWarnings = FALSE, recursive = TRUE)

set.seed(1L)
summaries <- read_summary_inputs()
set.seed(1L)

figure_specs <- list(
  deg = list(
    plot = build_deg_figure(summaries$deg),
    filename = "deg_monte_carlo_false_positive_rates",
    pdf_size = c(12, 7),
    png_size = c(10.5, 6)
  ),
  ora = list(
    plot = build_ora_figure(summaries$ora),
    filename = "ora_monte_carlo_false_positive_rates",
    pdf_size = c(12, 7),
    png_size = c(10.5, 6)
  ),
  gsea = list(
    plot = build_gsea_figure(summaries$gsea),
    filename = "gsea_monte_carlo_false_positive_rates",
    pdf_size = c(12, 7),
    png_size = c(10.5, 6)
  ),
  tgx_ddi = list(
    plot = build_tgx_ddi_figure(summaries$tgx_ddi),
    filename = "tgx_ddi_monte_carlo_false_positive_rates",
    pdf_size = c(12, 8),
    png_size = c(10.5, 7)
  ),
  tpod_fdr = list(
    plot = build_tpod_fdr_figure(summaries$tpod_fdr),
    filename = "tpod_fdr_grid",
    pdf_size = c(24, 13),
    png_size = c(24, 13)
  ),
  tpod_robustness = list(
    plot = build_tpod_robustness_figure(summaries$tpod_fdr),
    filename = "tpod_fdr_robustness",
    pdf_size = c(11, 10),
    png_size = c(11, 10)
  ),
  tpod_parameter_settings = list(
    plot = build_tpod_parameter_setting_figure(summaries$tpod_fdr),
    filename = "tpod_fdr_by_parameter_setting",
    pdf_size = c(12, 10),
    png_size = c(12, 10)
  ),
  tpod_fdr_supplementary = list(
    plot = build_tpod_fdr_supplementary_figure(summaries$tpod_fdr),
    filename = "tpod_fdr_grid_supplementary",
    pdf_size = c(24, 13),
    png_size = c(24, 13)
  )
)

for (name in names(figure_specs)) {
  spec <- figure_specs[[name]]
  print(paste0("Writing ", name))
  ggplot2::ggsave(
    filename = fdr_path("figures", paste0(spec$filename, ".pdf")),
    plot = spec$plot,
    width = spec$pdf_size[[1L]],
    height = spec$pdf_size[[2L]],
    units = "in",
    device = cairo_pdf
  )
  ggplot2::ggsave(
    filename = fdr_path("figures", paste0(spec$filename, ".png")),
    plot = spec$plot,
    width = spec$png_size[[1L]],
    height = spec$png_size[[2L]],
    units = "in",
    dpi = 600,
    device = png
  )
}

message("Wrote figures to ", fdr_path("figures"))
