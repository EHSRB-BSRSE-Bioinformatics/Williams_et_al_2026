source("R/config.R")
source("R/recurrence_figures.R")
summary <- read_s56_summary()
expected <- list(
    TK6 = c(
        "CD55",
        "CFL2",
        "CNOT10",
        "FKBP4",
        "LBHD1",
        "PHF21A",
        "RPL26",
        "RRM2",
        "SCAPER",
        "SPECC1L",
        "ZNF782"
    ),
    `MCF-7` = c(
        "COX8A",
        "H2AC13",
        "H2AC17",
        "H2AC8",
        "H3C7",
        "MALAT1",
        "POTEG",
        "PSMD6",
        "RALGPS2",
        "TMSB10"
    ),
    HepaRG = c(
        "AGXT",
        "APOA1",
        "BNIP3",
        "ERO1A",
        "HAO2",
        "LDHA",
        "NXPH4",
        "PGK1",
        "SCP2",
        "SERPIND1",
        "TMEM45A",
        "TMSB10"
    ),
    `Liver spheroids` = c(
        "CFHR2",
        "FAH",
        "FTH1",
        "GDA",
        "HPN",
        "NBPF8",
        "QPRT",
        "SERPIND1",
        "TMEM123",
        "UGT2B7"
    )
)
for (model in names(expected)) {
    stopifnot(setequal(
        summary$selected$symbol[summary$selected$cell_line == model],
        expected[[model]]
    ))
}
stopifnot(
    length(unique(summary$selected$feature)) == 41L,
    nrow(summary$selected) == 43L,
    length(unique(summary$expression$sample_id)) == 301L
)
screening <- readRDS(file.path(s56_reference_dir(), "screening_counts.rds"))
stopifnot(isTRUE(all.equal(
    s56_select(screening),
    summary$selected,
    check.attributes = FALSE
)))
counts <- as.matrix(read.delim(
    "data/input/count_table.tsv",
    row.names = 1L,
    check.names = FALSE
))
data <- summary$expression
raw <- counts[cbind(
    match(data$feature, rownames(counts)),
    match(data$sample_id, colnames(counts))
)]
sf <- summary$size_factors$size_factor[match(
    data$sample_id,
    summary$size_factors$sample_id
)]
stopifnot(isTRUE(all.equal(
    data$log2_normalized_count,
    log2(raw / sf + 1),
    tolerance = 1e-12
)))

# Boundary ties and repeated n-specific features must not become pooled top ten.
toy <- do.call(
    rbind,
    lapply(cell_line_levels, function(model) {
        grid <- expand.grid(
            feature = paste0("GENE", 1:12, "_1"),
            n = 3:10,
            stringsAsFactors = FALSE
        )
        grid$count <- 1L
        grid$count[grid$n == 3L] <- c(20:12, 11L, 11L, 1L)
        grid$cell_line <- model
        grid
    })
)
stopifnot(all(table(s56_select(toy)$cell_line) == 11L))
bad <- summary
bad$endpoint <- "final_DEGs"
stopifnot(inherits(try(validate_s56_summary(bad), silent = TRUE), "try-error"))
bad <- summary
bad$heatmap$count[1L] <- -1L
stopifnot(inherits(try(validate_s56_summary(bad), silent = TRUE), "try-error"))
bad <- summary
bad$stars <- bad$stars[-1L, ]
stopifnot(inherits(try(validate_s56_summary(bad), silent = TRUE), "try-error"))
correlations <- s56_baseline_correlations(summary)
stopifnot(
    nrow(correlations$pairs) == 211L,
    all(abs(correlations$pairs$rho) <= 1),
    all(correlations$summary$pairs == c(55L, 45L, 66L, 45L))
)
for (model in cell_line_levels) {
    row <- correlations$pairs[correlations$pairs$cell_line == model, ][1L, ]
    d <- data[data$cell_line == model, ]
    a <- d[d$feature == row$feature1, ]
    b <- d[d$feature == row$feature2, ]
    b <- b[match(a$sample_id, b$sample_id), ]
    stopifnot(isTRUE(all.equal(
        row$rho,
        cor(
            a$log2_normalized_count,
            b$log2_normalized_count,
            method = "spearman"
        )
    )))
}
s5 <- build_s5_recurrence_figure(summary)
s6 <- build_s6_expression_figure(summary)
stopifnot(
    inherits(s5, "ggplot"),
    inherits(s6, "ggplot"),
    identical(s6$layers[[2L]]$position$seed, 1L),
    nrow(s6$layers[[3L]]$data) == nrow(summary$stars)
)
stopifnot(
    any(!summary$heatmap$eligible_for_test),
    any(summary$heatmap$eligible_for_test & summary$heatmap$count == 0L),
    all(is.na(s5$data$plot_count) == !s5$data$eligible_for_test)
)
fill <- s5$scales$get_scales("fill")
stopifnot(
    identical(fill$map(NA_real_), "grey75"),
    identical(toupper(fill$map(0)), "#000000")
)
invisible(ggplot2::ggplot_build(s5))
invisible(ggplot2::ggplot_build(s6))
cat(
    "S5/S6 historical selection, archived normalization, stars, correlations and plot builds verified.\n"
)
