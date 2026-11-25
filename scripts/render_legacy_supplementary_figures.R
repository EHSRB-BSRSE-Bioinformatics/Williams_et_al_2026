script <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/figure_registry.R")
source("R/legacy_supplementary_figures.R")
legacy_assert(length(commandArgs(TRUE)) == 0L, "Render takes no arguments")
x <- readRDS(legacy_summary_path)
plots <- legacy_plots(x)
outputs <- character()
for (key in names(plots)) {
    width <- if (key %in% c("S7", "S9")) {
        10
    } else if (key == "S12") {
        15
    } else {
        13
    }
    height <- if (key == "S12") {
        max(8, 2 + length(x$pathways$order) * 0.24)
    } else {
        6
    }
    for (extension in c("png", "pdf")) {
        path <- figure_filename(tolower(key), extension)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(
            path,
            plots[[key]],
            width = width,
            height = height,
            units = "in",
            dpi = 300,
            bg = "white",
            device = if (extension == "pdf") grDevices::cairo_pdf else "png",
            limitsize = FALSE
        )
        outputs <- c(outputs, path)
        message("Wrote ", path)
    }
}
legacy_write_tsv(
    rbind(
        legacy_provenance(c(
            legacy_summary_path,
            "R/legacy_supplementary_figures.R",
            "R/figure_registry.R",
            "scripts/render_legacy_supplementary_figures.R"
        )),
        legacy_provenance(outputs, "figure_output")
    ),
    "data/derived/legacy_supplementary_figure_provenance.tsv"
)
