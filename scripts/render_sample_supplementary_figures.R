script <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/figure_registry.R")
source("R/sample_supplementary_figures.R")
args <- commandArgs(TRUE)
sample_supp_assert(
    length(args) == 0L ||
        identical(args, "--s3-raw") ||
        identical(args, "--s3-normalized"),
    "Usage: Rscript scripts/render_sample_supplementary_figures.R [--s3-normalized | --s3-raw]"
)
scale <- if (identical(args, "--s3-normalized")) "normalized" else "raw"
summary_path <- "data/derived/sample_supplementary_summary.rds"
summary <- readRDS(summary_path)
sample_supp_assert(
    identical(summary$schema_version, 3L),
    "Unsupported summary schema"
)
plots <- sample_supp_s3_plots(summary, scale)
render <- function(key, width, height, draw) {
    for (extension in c("png", "pdf")) {
        path <- figure_filename(key, extension)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        if (extension == "png") {
            grDevices::png(
                path,
                width = width,
                height = height,
                units = "in",
                res = 300,
                bg = "white",
                type = "cairo"
            )
        } else {
            grDevices::pdf(
                path,
                width = width,
                height = height,
                bg = "white",
                useDingbats = FALSE
            )
        }
        tryCatch(draw(), finally = grDevices::dev.off())
    }
}
render("s1", 10, 10, function() sample_supp_draw_s1(summary))
render("s3", 11, 12, function() sample_supp_draw_s3(plots))
extensions <- c("png", "pdf")
outputs <- c(figure_filename("s1", extensions), figure_filename("s3", extensions))
sample_supp_write_tsv(
    data.frame(
        path = outputs,
        sha256 = vapply(outputs, sample_supp_sha256, ""),
        summary_sha256 = sample_supp_sha256(summary_path),
        s3_scale = scale,
        background = "white",
        png_dpi = 300,
        R_version = as.character(getRversion()),
        ggplot2_version = as.character(utils::packageVersion("ggplot2"))
    ),
    "data/derived/sample_supplementary_figure_provenance.tsv"
)
cat("Rendered Fig_S1 and Fig_S3 PNG/PDF; S3 scale:", scale, "\n")
