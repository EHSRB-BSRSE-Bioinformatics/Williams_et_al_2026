script <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/sample_supplementary_figures.R")
summary <- readRDS("data/derived/sample_supplementary_summary.rds")
inputs <- sample_supp_read_inputs()
build <- function(
    distance = summary$distance,
    annotation = summary$annotation,
    selection = summary$selection,
    data = inputs
) {
    sample_supp_build(
        distance,
        annotation,
        selection,
        data,
        summary$s3_decisions
    )
}
rebuilt <- build()
for (name in names(rebuilt)) {
    stopifnot(identical(rebuilt[[name]], summary[[name]]))
}
stopifnot(
    identical(dim(summary$distance), c(365L, 365L)),
    identical(
        as.integer(table(factor(
            summary$annotation$chemical,
            levels = c("DMSO", "Media")
        ))),
        c(301L, 64L)
    ),
    identical(summary$tree$method, "complete"),
    identical(summary$tree$labels, rownames(summary$distance)),
    identical(summary$panels$feature, sample_supp_panels$feature),
    identical(summary$panels$cell_line, sample_supp_panels$cell_line)
)
# A Euclidean distance must satisfy the triangle inequality (all pairs, ten anchors).
d <- summary$distance
for (anchor in seq_len(10)) {
    stopifnot(all(d <= outer(d[, anchor], d[anchor, ], "+") + 1e-8))
}
expression <- summary$expression
frozen <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
evidence <- read.delim("data/reference/sample_supplementary/s3_raster_evidence.tsv")
recovered <- sample_supp_select_s3(evidence, frozen, inputs$counts)
stopifnot(
    identical(recovered$selection, summary$selection),
    isTRUE(all.equal(recovered$decisions, summary$s3_decisions)),
    identical(summary$s3_decisions$trial_name, c(
        "Liver spheroids_4_4_712", "HepaRG_4_4_582", "MCF-7_4_4_59", "TK6_4_4_136"
    )),
    identical(recovered$decisions$trials_scanned, rep(8000L, 4L)),
    identical(recovered$decisions$pair_score, c(6, 6, 9, 4)),
    identical(recovered$decisions$runner_up_score, c(65, 61, 60, 39))
)
for (trial in unique(summary$selection$trial_name)) {
    x <- summary$selection[summary$selection$trial_name == trial, ]
    stopifnot(
        identical(x$sample_id, frozen[[trial]]$sample_id),
        identical(x$group, as.character(frozen[[trial]]$group))
    )
}
for (i in seq_len(8)) {
    x <- expression[expression$panel == LETTERS[i], ]
    wells <- inputs$metadata$Sample_ID[
        inputs$metadata$cell_line == sample_supp_panels$cell_line[i] &
            inputs$metadata$chemical == "DMSO"
    ]
    stopifnot(
        setequal(x$sample_id[x$group == "All DMSO"], wells),
        sum(x$group == "A") == 4L,
        sum(x$group == "B") == 4L,
        !anyDuplicated(x[c("sample_id", "group")]),
        all(x$normalized_count == x$raw_count / x$size_factor),
        all(
            x$raw_count ==
                inputs$counts[cbind(
                    match(x$input_feature, rownames(inputs$counts)),
                    match(x$sample_id, colnames(inputs$counts))
                )]
        )
    )
    stopifnot(identical(unique(x$trial_name), summary$s3_decisions$trial_name[
        summary$s3_decisions$cell_line == sample_supp_panels$cell_line[i]
    ]))
}
fails <- function(expr) {
    stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
}
bad <- d
bad[1, 2] <- bad[1, 2] + 1
fails(sample_supp_validate_distance(bad, summary$annotation))
bad <- summary$annotation
bad$sample_id[1] <- bad$sample_id[2]
fails(sample_supp_validate_distance(d, bad))
bad <- summary$selection
bad$group[1] <- "B"
fails(build(selection = bad))
bad <- inputs
bad$factors$size_factor[1] <- 0
fails(build(data = bad))
bad <- inputs
bad$metadata$chemical[match(
    summary$selection$sample_id[1],
    bad$metadata$Sample_ID
)] <- "Media"
fails(build(data = bad))
aliases <- inputs
rownames(aliases$counts)[
    rownames(aliases$counts) == "H2AC17_33569"
] <- "HIST1H2AI_33569"
alias_summary <- build(data = aliases)
stopifnot(
    identical(alias_summary$expression$raw_count, expression$raw_count),
    identical(
        alias_summary$expression$normalized_count,
        expression$normalized_count
    ),
    all(
        alias_summary$expression$input_feature[
            alias_summary$expression$panel == "E"
        ] ==
            "HIST1H2AI_33569"
    )
)
bad <- inputs
bad$factors <- bad$factors[-1L, ]
fails(build(data = bad))
fails(sample_supp_select_s3(evidence[-1L, ], frozen, inputs$counts))
fails(sample_supp_select_s3(evidence, frozen[-1L], inputs$counts))
plots <- sample_supp_s3_plots(summary)
stopifnot(length(plots) == 8L)
for (plot in plots) {
    stopifnot(
        identical(levels(plot$data$group), c("All DMSO", "B", "A")),
        identical(plot$data$display_count, plot$data$raw_count),
        all(is.finite(ggplot2::ggplot_build(plot)$data[[1]]$y))
    )
}
# Preserve zeros explicitly rather than silently dropping them on a log axis.
zero <- summary
zero$expression$raw_count[1] <- 0
stopifnot(
    ggplot2::ggplot_build(sample_supp_s3_plots(zero)[[1]])$data[[1]]$y[1] == 0
)
seed_before <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    .Random.seed
} else {
    NULL
}
invisible(sample_supp_s3_plots(summary))
seed_after <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    .Random.seed
} else {
    NULL
}
stopifnot(identical(seed_before, seed_after))
provenance <- read.delim("data/derived/sample_supplementary_provenance.tsv")
stopifnot(all(
    vapply(provenance$path, sample_supp_sha256, "") == provenance$sha256
))
figure_provenance <- "data/derived/sample_supplementary_figure_provenance.tsv"
if (file.exists(figure_provenance)) {
    p <- read.delim(figure_provenance)
    stopifnot(
        nrow(p) == 4L,
        all(p$background == "white"),
        all(file.info(p$path)$size > 1000),
        all(vapply(p$path, sample_supp_sha256, "") == p$sha256),
        all(
            p$summary_sha256 ==
                sample_supp_sha256(
                    "data/derived/sample_supplementary_summary.rds"
                )
        )
    )
}
cat(
    "S1/S3 tests passed: frozen distances, four unique n=4 raster matches, raw counts and provenance.\n"
)
