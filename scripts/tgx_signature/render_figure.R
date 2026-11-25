#!/usr/bin/env Rscript
# Render the TGx-DDI signature supplementary figure from the compact analysis RDS.
script_file <- sub(
    "--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
    setwd(dirname(dirname(dirname(normalizePath(
        script_file,
        mustWork = TRUE
    )))))
}

source("R/config.R")
source("R/figure_registry.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
    i <- match(key, args)
    if (is.na(i)) default else args[[i + 1L]]
}
if (
    length(args) %% 2L != 0L ||
        (length(args) > 0L &&
            any(!grepl("^--", args[seq(1L, length(args), 2L)])))
) {
    stop(
        "Usage: Rscript scripts/tgx_signature/render_figure.R [--input PATH] [--output-dir PATH]"
    )
}

input_path <- get_arg(
    "--input",
    "data/derived/tgx_signature/tgx_signature_null.rds"
)
output_dir <- get_arg("--output-dir", fdr_path("figures"))
if (!file.exists(input_path)) {
    stop("Missing TGx-DDI signature analysis RDS: ", input_path)
}

result <- readRDS(input_path)
if (!is.list(result) || !is.data.frame(result$summary)) {
    stop("Input RDS must contain a summary data frame.")
}
required <- c(
    "cell",
    "n",
    "stage",
    "trials",
    "signature_trial_rate",
    "whole_trial_rate"
)
if (!all(required %in% names(result$summary))) {
    stop("TGx-DDI summary is missing required plotting columns.")
}
final <- result$summary[result$summary$stage == "final", , drop = FALSE]
cell_levels <- c("TK6", "MCF-7", "HepaRG", "Liver spheroids")
if (
    nrow(final) != length(cell_levels) * 8L ||
        any(final$trials != 1000L) ||
        anyNA(final$whole_trial_rate) ||
        any(!is.finite(final$signature_trial_rate)) ||
        any(!is.finite(final$whole_trial_rate)) ||
        any(!final$cell %in% cell_levels) ||
        !setequal(final$n, 3:10) ||
        anyDuplicated(final[c("cell", "n")])
) {
    stop("Final-call summary is incomplete or invalid for figure rendering.")
}

# --output-dir overrides the figures root, retaining the registered category.
spec <- figure_spec("s13_tgx_signature")
output_dir <- file.path(output_dir, basename(spec$directory))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
pdf_path <- file.path(output_dir, basename(figure_filename(spec$key, "pdf")))
png_path <- file.path(output_dir, basename(figure_filename(spec$key, "png")))

draw_figure <- function() {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(old_par), add = TRUE)
    for (cell in cell_levels) {
        d <- final[final$cell == cell, , drop = FALSE]
        d <- d[order(d$n), , drop = FALSE]
        whole_rate <- 100 * d$whole_trial_rate
        signature_rate <- 100 * d$signature_trial_rate
        ylim <- c(0, max(1, whole_rate, signature_rate))
        plot(
            d$n,
            whole_rate,
            type = "b",
            pch = 16,
            col = "#606060",
            ylim = ylim,
            xlim = c(3, 10),
            xlab = "Samples per group",
            ylab = "Trials with >=1 final DEG (%)",
            main = cell
        )
        lines(d$n, signature_rate, type = "b", pch = 17, col = "#007C91")
        legend(
            "topright",
            c("Whole transcriptome", "TGx-DDI genes"),
            col = c("#606060", "#007C91"),
            pch = c(16, 17),
            lty = 1,
            bty = "n",
            cex = 0.8
        )
    }
}

grDevices::pdf(pdf_path, width = 9, height = 7, useDingbats = FALSE)
tryCatch(draw_figure(), finally = grDevices::dev.off())
grDevices::png(
    png_path,
    width = 9,
    height = 7,
    units = "in",
    res = 600,
    type = "cairo"
)
tryCatch(draw_figure(), finally = grDevices::dev.off())

if (!nzchar(Sys.which("sha256sum"))) {
    stop("sha256sum is required to record figure provenance.")
}
provenance_paths <- c(
    input_path,
    "scripts/tgx_signature/render_figure.R",
    "R/figure_registry.R",
    pdf_path,
    png_path
)
hashes <- vapply(
    provenance_paths,
    function(path) {
        strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " +")[[
            1L
        ]][1L]
    },
    character(1)
)
absolute_paths <- normalizePath(provenance_paths, mustWork = TRUE)
root_prefix <- paste0(project_root(), "/")
local <- startsWith(absolute_paths, root_prefix)
absolute_paths[local] <- substring(absolute_paths[local], nchar(root_prefix) + 1L)
write.table(
    data.frame(
        role = c("input", "source", "source", "output", "output"),
        path = absolute_paths,
        sha256 = unname(hashes)
    ),
    file.path(dirname(input_path), "figure_provenance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
message("Wrote ", pdf_path, " and ", png_path)
