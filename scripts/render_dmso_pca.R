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

required_packages <- c("DESeq2", "ggplot2")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
    stop(
        "Missing packages required for the DMSO PCA: ",
        paste(missing_packages, collapse = ", "),
        ". Run this script in the optional 'deg' Pixi environment."
    )
}

metadata <- read.delim(
    fdr_path("data", "input", "metadata.QC_applied.txt"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
counts <- as.matrix(read.delim(
    fdr_path("data", "input", "count_table.tsv"),
    check.names = FALSE,
    row.names = 1
))

required_metadata <- c(
    "Sample_ID",
    "chemical",
    "cell_line",
    "technical_control",
    "reference_rna"
)
missing_columns <- setdiff(required_metadata, names(metadata))
if (length(missing_columns) > 0L) {
    stop(
        "Sample metadata is missing: ",
        paste(missing_columns, collapse = ", ")
    )
}
if (anyDuplicated(metadata$Sample_ID)) {
    stop("QC-applied metadata contains duplicate Sample_ID values.")
}
if (
    !is.logical(metadata$technical_control) ||
        !is.logical(metadata$reference_rna)
) {
    stop("Technical-control flags must be logical TRUE/FALSE values.")
}

# The matrix also contains media and unannotated control columns. Use only the
# post-QC DMSO records, excluding any technical or reference-RNA controls.
dmso_metadata <- metadata[
    !is.na(metadata$chemical) & metadata$chemical == "DMSO",
    ,
    drop = FALSE
]
if (
    anyNA(dmso_metadata$technical_control) || anyNA(dmso_metadata$reference_rna)
) {
    stop("DMSO metadata has missing technical-control flags.")
}
dmso_metadata <- dmso_metadata[
    !dmso_metadata$technical_control & !dmso_metadata$reference_rna,
    ,
    drop = FALSE
]
if (!nrow(dmso_metadata)) {
    stop(
        "No QC-applied DMSO samples remain after excluding technical controls."
    )
}

unknown_cell_lines <- setdiff(unique(dmso_metadata$cell_line), cell_line_levels)
missing_cell_lines <- setdiff(cell_line_levels, unique(dmso_metadata$cell_line))
if (length(unknown_cell_lines) > 0L || length(missing_cell_lines) > 0L) {
    stop("DMSO samples do not contain exactly the expected four cell models.")
}
sample_ids <- as.character(dmso_metadata$Sample_ID)
sample_columns <- match(sample_ids, colnames(counts))
if (anyNA(sample_columns)) {
    missing_samples <- sample_ids[is.na(sample_columns)]
    stop(
        "Count matrix is missing QC-applied DMSO samples: ",
        paste(missing_samples, collapse = ", ")
    )
}
if (anyDuplicated(sample_columns)) {
    stop("Selected DMSO Sample_ID values do not map to unique count columns.")
}

counts <- counts[, sample_columns, drop = FALSE]
storage.mode(counts) <- "numeric"
if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("Selected DMSO counts must be finite and non-negative.")
}
counts <- round(counts)
storage.mode(counts) <- "integer"

dmso_metadata$cell_line <- factor(
    dmso_metadata$cell_line,
    levels = cell_line_levels
)
rownames(dmso_metadata) <- dmso_metadata$Sample_ID

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = dmso_metadata,
    design = ~cell_line
)
dds <- dds[rowSums(DESeq2::counts(dds)) > 1L, ]
if (nrow(dds) < 2L) {
    stop("Too few genes remain for PCA after the historical count filter.")
}

# Match the historical report workflow: variance-stabilizing transformation,
# then PCA across all retained genes rather than only the default top 500.
vsd <- DESeq2::vst(dds, blind = FALSE)
pca_data <- DESeq2::plotPCA(
    vsd,
    intgroup = "cell_line",
    ntop = nrow(dds),
    returnData = TRUE
)
percent_variance <- round(100 * attr(pca_data, "percentVar"))
pca_data$cell_line <- factor(pca_data$cell_line, levels = cell_line_levels)

cell_line_colours <- c(
    "TK6" = "#0072B2",
    "MCF-7" = "#D55E00",
    "HepaRG" = "#009E73",
    "Liver spheroids" = "#CC79A7"
)

p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = PC1, y = PC2, colour = cell_line)
) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(shape = 16, size = 2.2, alpha = 0.8) +
    ggplot2::scale_colour_manual(
        values = cell_line_colours,
        breaks = cell_line_levels,
        name = "Cell model"
    ) +
    ggplot2::labs(
        title = "PCA of DMSO treated samples from four cell models",
        x = paste0("PC1: ", percent_variance[[1L]], "% variance"),
        y = paste0("PC2: ", percent_variance[[2L]], "% variance")
    ) +
    ggplot2::guides(
        colour = ggplot2::guide_legend(
            override.aes = list(size = 3, alpha = 1)
        )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "right",
        plot.title = ggplot2::element_text(face = "bold"),
        plot.caption = ggplot2::element_text(hjust = 0)
    )

output_prefix <- figure_path("s2_pca_dmso")
dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)
ggplot2::ggsave(
    figure_path("s2_pca_dmso", "pdf"),
    p,
    width = 7.5,
    height = 5.1,
    units = "in",
    device = grDevices::cairo_pdf
)
ggplot2::ggsave(
    figure_path("s2_pca_dmso", "png"),
    p,
    width = 7.5,
    height = 5.1,
    units = "in",
    dpi = 600
)

message(
    "Wrote PCA for ",
    ncol(dds),
    " QC-applied DMSO samples across ",
    length(cell_line_levels),
    " cell models."
)
message("Output: ", output_prefix, ".pdf and .png")
