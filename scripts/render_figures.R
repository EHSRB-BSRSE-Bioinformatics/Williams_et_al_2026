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

figures <- list(
  deg = build_deg_figure(summaries$deg),
  ora = build_ora_figure(summaries$ora),
  gsea = build_gsea_figure(summaries$gsea),
  tgx_ddi = build_tgx_ddi_figure(summaries$tgx_ddi)
)

for (name in names(figures)) {
  print(paste0("Writing ", name))
  ggplot2::ggsave(
    filename = fdr_path(
      "figures",
      paste0(name, "_monte_carlo_false_positive_rates.pdf")
    ),
    plot = figures[[name]],
    width = 12,
    height = if (name == "tgx_ddi") 8 else 7,
    units = "in",
    device = cairo_pdf
  )
  ggplot2::ggsave(
    filename = fdr_path(
      "figures",
      paste0(name, "_monte_carlo_false_positive_rates.png")
    ),
    plot = figures[[name]],
    width = 10.5,
    height = if (name == "tgx_ddi") 7 else 6,
    units = "in",
    dpi = 300,
    device = png
  )
}

message("Wrote figures to ", fdr_path("figures"))
