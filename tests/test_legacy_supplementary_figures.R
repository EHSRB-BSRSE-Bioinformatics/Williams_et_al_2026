script <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/legacy_supplementary_figures.R")
x <- readRDS(legacy_summary_path)
legacy_validate(x)
ref <- legacy_reference
screen <- readRDS(file.path(ref, "screening_counts.rds"))
if (file.exists("data/reference/s56/screening_counts.rds")) {
    original_screen <- readRDS("data/reference/s56/screening_counts.rds")
    columns <- c("cell_line", "n", "feature", "count")
    canonical <- function(d) {
        d <- d[order(d$cell_line, d$n, d$feature), columns]
        rownames(d) <- NULL
        d
    }
    stopifnot(identical(canonical(screen), canonical(original_screen)))
}
expr <- readRDS(file.path(ref, "expression_reference.rds"))
stopifnot(
    length(expr$features) == 22460L,
    nrow(expr$factors) == 301L,
    nrow(x$expression) == 4L * 22460L,
    nrow(x$scatter) == 220809L
)
counts <- as.matrix(read.delim(
    "data/input/count_table.tsv",
    row.names = 1L,
    check.names = FALSE
))
metadata <- read.delim("data/input/metadata.QC_applied.txt")
trials <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
final <- readRDS(file.path(ref, "final_trial_counts.rds"))
compiled <- as.data.frame(readRDS("data/derived/deg_monte_carlo_summary.rds"))
compiled <- compiled[
    !is.na(compiled$Filter) & as.character(compiled$Filter) == "Final DEGs",
]
stopifnot(
    !anyDuplicated(compiled$name),
    setequal(final$trial_id, compiled$name),
    identical(
        final$n_final,
        as.integer(compiled$nDEGs[match(final$trial_id, compiled$name)])
    )
)
# Independent normalization spot-check for every model and 15 distributed probes.
features <- expr$features[unique(round(seq(
    1,
    length(expr$features),
    length.out = 15
)))]
for (model in legacy_models) {
    wells <- expr$factors$sample_id[
        metadata$cell_line[match(expr$factors$sample_id, metadata$Sample_ID)] ==
            model
    ]
    sf <- expr$factors$size_factor[match(wells, expr$factors$sample_id)]
    z <- log2(sweep(counts[features, wells, drop = FALSE], 2L, sf, "/") + 1)
    observed <- x$expression[x$expression$cell_line == model, ]
    observed <- observed[match(features, observed$feature), ]
    stopifnot(
        isTRUE(all.equal(
            unname(rowMeans(z)),
            observed$mean_expression,
            tolerance = 1e-12
        )),
        isTRUE(all.equal(
            unname(apply(z, 1L, var)),
            observed$variance,
            tolerance = 1e-12
        ))
    )
}
# Literal historical left_join without n must equal exported scatter multiplicities.
literal <- merge(
    x$expression,
    screen[c("cell_line", "feature", "count")],
    by = c("cell_line", "feature"),
    all.x = TRUE
)
literal$count[is.na(literal$count)] <- 0L
key <- function(d) sort(paste(d$cell_line, d$feature, d$count))
stopifnot(identical(key(literal), key(x$scatter)))
# Independent scan for one actual sample/n per model and both randomized groups.
for (model in legacy_models) {
    s <- x$samples[x$samples$cell_line == model & x$samples$n == 3L, ][1L, ]
    ids <- names(trials)[startsWith(names(trials), paste0(model, "_3_3_"))]
    contains <- vapply(
        trials[ids],
        function(t) s$sample_id %in% t$sample_id,
        logical(1)
    )
    hit <- final$n_final[match(ids, final$trial_id)] > 0
    stopifnot(
        sum(contains) == s$total_experiments,
        sum(contains & hit) == s$n_high_DEG
    )
}
# Plate is the historical last drawn row, not an average across n.
keys <- paste(x$samples$cell_line, x$samples$row, x$samples$column)
for (i in seq_len(nrow(x$plate))) {
    p <- x$plate[i, ]
    j <- tail(which(keys == paste(p$cell_line, p$row, p$column)), 1L)
    stopifnot(identical(
        p$high_DEG_proportion,
        x$samples$high_DEG_proportion[j]
    ))
}
# Boundary ties are included; grouping is model/n and selection is by name union.
toy <- data.frame(
    cell_line = "TK6",
    n = 3L,
    pathway = paste0("P", 1:22),
    count = c(100:82, 80L, 80L, 1L)
)
anno <- data.frame(pathway = toy$pathway, path_name = toy$pathway)
test <- legacy_pathways(toy, anno)
stopifnot(length(test$selected) == 21L, !"P22" %in% test$selected)
# No silent pooling of n or silently replacing zero/unknown heatmap cells.
toy_expr <- data.frame(
    cell_line = "TK6",
    feature = c("a", "b"),
    mean_expression = c(1, 2),
    variance = c(1, 2)
)
toy_screen <- data.frame(
    cell_line = "TK6",
    feature = "a",
    n = c(3L, 4L),
    count = c(4L, 8L)
)
joined <- legacy_scatter(toy_expr, toy_screen)
stopifnot(nrow(joined) == 3L, identical(sort(joined$count), c(0L, 4L, 8L)))
# References and public inputs are cryptographically pinned.
manifest <- read.delim(file.path(ref, "source_provenance.tsv"))
for (i in which(manifest$role %in% c("reference_output", "public_input"))) {
    stopifnot(identical(legacy_sha256(manifest$path[i]), manifest$sha256[i]))
}
for (path in c(
    "data/derived/legacy_supplementary_data_provenance.tsv",
    "data/derived/legacy_supplementary_figure_provenance.tsv"
)) {
    if (!file.exists(path)) {
        next
    }
    p <- read.delim(path)
    stopifnot(!any(grepl("^/|[A-Za-z]:[/\\\\]|\\.\\./", p$path)))
    for (i in seq_len(nrow(p))) {
        stopifnot(identical(legacy_sha256(p$path[i]), p$sha256[i]))
    }
}
plots <- legacy_plots(x)
stopifnot(identical(names(plots), c("S7", "S9", "S10", "S11", "S12")))
for (p in plots) {
    invisible(ggplot2::ggplot_build(p))
}
message(
    "PASS: historical definitions, normalization, joins, participation, overplot, ties, hashes and five plot builds"
)
