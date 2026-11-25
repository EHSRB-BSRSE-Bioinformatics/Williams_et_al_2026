# Historical S5/S6 screening endpoint; distinct from pooled FINAL DEG rankings.
s56_endpoint <- "archived_DESeqRes_padj_lt_0.05"
s56_stems <- c(
    TK6 = "tk6",
    `MCF-7` = "mcf7",
    HepaRG = "heparg",
    `Liver spheroids` = "spheroids"
)
s56_reference_dir <- function() fdr_path("data", "reference", "s56")
s56_summary_path <- function() {
    fdr_path("data", "derived", "s56_figure_summary.rds")
}

s56_assert <- function(value, message) {
    if (!isTRUE(value)) stop(message, call. = FALSE)
}

s56_sha256 <- function(path) {
    s56_assert(file.exists(path), paste("Missing provenance input:", path))
    result <- system2("sha256sum", shQuote(path), stdout = TRUE)
    s56_assert(
        length(result) == 1L && grepl("^[[:xdigit:]]{64} ", result),
        "SHA-256 failed."
    )
    sub(" .*", "", result)
}

s56_write_tsv <- function(data, path) {
    write.table(
        data,
        path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        na = "NA"
    )
}

s56_select <- function(screening) {
    s56_assert(
        all(c("cell_line", "n", "feature", "count") %in% names(screening)),
        "Invalid screening columns."
    )
    s56_assert(
        !anyNA(screening) &&
            !anyDuplicated(screening[c("cell_line", "n", "feature")]),
        "Missing or duplicate screening rows."
    )
    s56_assert(
        setequal(unique(screening$cell_line), cell_line_levels) &&
            all(screening$n %in% sample_size_levels) &&
            all(screening$count > 0) &&
            all(screening$count <= n_trials_default) &&
            all(screening$count == as.integer(screening$count)),
        "Invalid screening counts or strata."
    )
    do.call(
        rbind,
        lapply(cell_line_levels, function(model) {
            data <- screening[screening$cell_line == model, ]
            s56_assert(
                setequal(unique(data$n), sample_size_levels),
                "Incomplete screening sample sizes."
            )
            # Literal historical top_n(10) on feature-by-n rows, NOT pooled sums.
            cutoff <- sort(data$count, decreasing = TRUE)[min(10L, nrow(data))]
            features <- sort(unique(data$feature[data$count >= cutoff]))
            data.frame(
                cell_line = model,
                feature = features,
                symbol = sub("_[0-9]+$", "", features),
                selection_cutoff = cutoff,
                stringsAsFactors = FALSE
            )
        })
    )
}

build_s56_summary <- function(
    screening,
    size_factors,
    counts,
    metadata,
    eligibility
) {
    selected <- s56_select(screening)
    features <- sort(unique(selected$feature))
    s56_assert(
        all(c("sample_id", "size_factor") %in% names(size_factors)) &&
            !anyNA(size_factors) &&
            !anyDuplicated(size_factors$sample_id) &&
            all(is.finite(size_factors$size_factor)) &&
            all(size_factors$size_factor > 0),
        "Invalid archived S6 size factors."
    )
    s56_assert(
        !anyDuplicated(metadata$Sample_ID) &&
            !anyDuplicated(colnames(counts)) &&
            !anyDuplicated(rownames(counts)),
        "Duplicate input identifiers."
    )
    wells <- sort(size_factors$sample_id)
    md <- metadata[match(wells, metadata$Sample_ID), ]
    s56_assert(
        all(features %in% rownames(counts)) &&
            all(wells %in% colnames(counts)) &&
            !anyNA(md$Sample_ID) &&
            all(md$chemical == "DMSO") &&
            setequal(unique(md$cell_line), cell_line_levels),
        "S6 features/wells do not match DMSO inputs."
    )
    factors <- size_factors$size_factor[match(wells, size_factors$sample_id)]
    raw <- counts[features, wells, drop = FALSE]
    s56_assert(
        all(is.finite(raw)) && all(raw >= 0),
        "Invalid expression counts."
    )
    expression <- log2(sweep(raw, 2L, factors, "/") + 1)
    grid <- expand.grid(
        feature = features,
        n = sample_size_levels,
        cell_line = cell_line_levels,
        stringsAsFactors = FALSE
    )
    heatmap <- merge(
        grid,
        screening,
        by = c("cell_line", "n", "feature"),
        all.x = TRUE,
        sort = FALSE
    )
    heatmap$count[is.na(heatmap$count)] <- 0L
    s56_assert(
        all(
            c("cell_line", "n", "feature", "tested_trials", "bh_trials") %in%
                names(eligibility)
        ) &&
            !anyNA(eligibility) &&
            !anyDuplicated(eligibility[c("cell_line", "n", "feature")]),
        "Missing or duplicate screening eligibility rows."
    )
    heatmap <- merge(
        heatmap,
        eligibility,
        by = c("cell_line", "n", "feature"),
        all.x = TRUE,
        sort = FALSE
    )
    s56_assert(
        !anyNA(heatmap$tested_trials) && !anyNA(heatmap$bh_trials),
        "Incomplete screening eligibility for selected probes."
    )
    heatmap$eligible_for_test <- heatmap$tested_trials > 0L
    heatmap$symbol <- sub("_[0-9]+$", "", heatmap$feature)
    nonzero <- heatmap[heatmap$count > 0L, ]
    # Original reorder(gene, count) uses mean over the observed nonzero rows.
    means <- tapply(nonzero$count, nonzero$symbol, mean)
    gene_order <- names(means)[order(means, names(means))]
    expression_long <- data.frame(
        feature = rep(features, times = length(wells)),
        sample_id = rep(wells, each = length(features)),
        cell_line = rep(md$cell_line, each = length(features)),
        log2_normalized_count = as.vector(expression),
        stringsAsFactors = FALSE
    )
    expression_long$symbol <- sub("_[0-9]+$", "", expression_long$feature)
    # One star per probe/model if any padj-only screening call at any n.
    stars <- aggregate(count ~ cell_line + feature + symbol, heatmap, sum)
    stars <- stars[stars$count > 0L, c("cell_line", "feature", "symbol")]
    # A per-facet star height replaces the original repeated/global-max star.
    limits <- aggregate(log2_normalized_count ~ feature, expression_long, max)
    stars$y <- limits$log2_normalized_count[match(
        stars$feature,
        limits$feature
    )] +
        0.5
    list(
        schema_version = 2L,
        endpoint = s56_endpoint,
        selection_rule = "top_10_feature_by_n_rows_per_model_with_boundary_ties",
        expression_rule = "log2(raw_count / archived_dds_all_size_factor + 1)",
        trials_per_stratum = n_trials_default,
        selected = selected,
        heatmap = heatmap,
        expression = expression_long,
        stars = stars,
        gene_order = gene_order,
        size_factors = size_factors
    )
}

validate_s56_summary <- function(summary) {
    s56_assert(
        identical(summary$schema_version, 2L) &&
            identical(summary$endpoint, s56_endpoint) &&
            identical(
                summary$selection_rule,
                "top_10_feature_by_n_rows_per_model_with_boundary_ties"
            ) &&
            identical(
                summary$expression_rule,
                "log2(raw_count / archived_dds_all_size_factor + 1)"
            ),
        "Wrong S5/S6 endpoint or normalization."
    )
    selected <- summary$selected
    features <- unique(selected$feature)
    heatmap <- summary$heatmap
    expression <- summary$expression
    s56_assert(
        !anyNA(selected) &&
            !anyDuplicated(selected[c("cell_line", "feature")]) &&
            !anyDuplicated(sub("_[0-9]+$", "", features)),
        "Invalid selected probe identities."
    )
    s56_assert(
        nrow(heatmap) == length(features) * 32L &&
            !anyNA(heatmap) &&
            !anyDuplicated(heatmap[c("cell_line", "n", "feature")]) &&
            all(heatmap$count >= 0 & heatmap$count <= n_trials_default) &&
            setequal(heatmap$feature, features) &&
            setequal(heatmap$cell_line, cell_line_levels) &&
            setequal(heatmap$n, sample_size_levels),
        "Invalid S5 grid."
    )
    s56_assert(
        all(
            heatmap$tested_trials >= 0 &
                heatmap$tested_trials <= n_trials_default
        ) &&
            all(
                heatmap$bh_trials >= 0 &
                    heatmap$bh_trials <= heatmap$tested_trials
            ) &&
            all(heatmap$count <= heatmap$bh_trials) &&
            all(heatmap$eligible_for_test == (heatmap$tested_trials > 0L)),
        "Invalid S5 testing eligibility."
    )
    s56_assert(
        !anyNA(expression) &&
            !anyDuplicated(expression[c("feature", "sample_id")]) &&
            all(is.finite(expression$log2_normalized_count)) &&
            all(expression$log2_normalized_count >= 0) &&
            setequal(expression$feature, features) &&
            nrow(expression) == length(features) * nrow(summary$size_factors) &&
            setequal(expression$sample_id, summary$size_factors$sample_id) &&
            setequal(summary$gene_order, selected$symbol),
        "Invalid S6 expression data."
    )
    expected_stars <- unique(heatmap[
        heatmap$count > 0L,
        c("cell_line", "feature")
    ])
    s56_assert(
        setequal(
            paste(expected_stars$cell_line, expected_stars$feature),
            paste(summary$stars$cell_line, summary$stars$feature)
        ) &&
            !anyDuplicated(summary$stars[c("cell_line", "feature")]),
        "Incorrect S6 screening stars."
    )
    invisible(summary)
}

read_s56_summary <- function(path = s56_summary_path()) {
    summary <- readRDS(path)
    validate_s56_summary(summary)
    summary
}

s56_baseline_correlations <- function(summary) {
    validate_s56_summary(summary)
    pairs <- do.call(
        rbind,
        lapply(cell_line_levels, function(model) {
            features <- summary$selected$feature[
                summary$selected$cell_line == model
            ]
            data <- summary$expression[summary$expression$cell_line == model, ]
            wells <- sort(unique(data$sample_id))
            mat <- matrix(
                NA_real_,
                length(features),
                length(wells),
                dimnames = list(features, wells)
            )
            ii <- cbind(
                match(data$feature, features),
                match(data$sample_id, wells)
            )
            keep <- !is.na(ii[, 1L])
            mat[ii[keep, , drop = FALSE]] <- data$log2_normalized_count[keep]
            s56_assert(
                !anyNA(mat) && all(apply(mat, 1L, sd) > 0),
                "Incomplete or constant correlation features."
            )
            rho <- cor(t(mat), method = "spearman")
            index <- t(combn(seq_along(features), 2L))
            data.frame(
                cell_line = model,
                feature1 = features[index[, 1L]],
                feature2 = features[index[, 2L]],
                rho = rho[index],
                wells = length(wells)
            )
        })
    )
    context <- do.call(
        rbind,
        lapply(cell_line_levels, function(model) {
            data <- pairs[pairs$cell_line == model, ]
            data.frame(
                cell_line = model,
                selected_probes = sum(summary$selected$cell_line == model),
                wells = unique(data$wells),
                pairs = nrow(data),
                median_rho = median(data$rho),
                median_abs_rho = median(abs(data$rho)),
                min_rho = min(data$rho),
                max_rho = max(data$rho)
            )
        })
    )
    list(pairs = pairs, summary = context)
}

build_s5_recurrence_figure <- function(summary) {
    validate_s56_summary(summary)
    data <- summary$heatmap
    data$cell_line <- factor(data$cell_line, levels = cell_line_levels)
    data$symbol <- factor(data$symbol, levels = summary$gene_order)
    data$plot_count <- ifelse(data$eligible_for_test, data$count, NA_real_)
    ggplot2::ggplot(
        data,
        ggplot2::aes(x = factor(n), y = symbol, fill = plot_count)
    ) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(~cell_line) +
        ggplot2::scale_fill_gradientn(
            colours = c(
                "black",
                "#1c1044",
                "#4f127b",
                "#812581",
                "#b5367a",
                "#e55063",
                "#fb8761",
                "#fec287",
                "#fcfdbf"
            ),
            na.value = "grey75",
            limits = c(0, max(data$count)),
            name = "Trials with\na screening call"
        ) +
        ggplot2::labs(
            title = "Recurrent false-positive DEG calls",
            subtitle = "BH-adjusted p < 0.05",
            x = "Sample size per group",
            y = "Gene symbol (one selected assay probe)",
            caption = "Grey: no eligible DESeq2 test at that sample size. Black: tested, zero screening calls."
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 8)
        )
}

build_s6_expression_figure <- function(summary) {
    validate_s56_summary(summary)
    data <- summary$expression
    stars <- summary$stars
    levels <- sort(unique(data$symbol))
    data$cell_line <- factor(data$cell_line, levels = cell_line_levels)
    stars$cell_line <- factor(stars$cell_line, levels = cell_line_levels)
    data$symbol <- factor(data$symbol, levels = levels)
    stars$symbol <- factor(stars$symbol, levels = levels)
    ggplot2::ggplot(
        data,
        ggplot2::aes(x = cell_line, y = log2_normalized_count, fill = cell_line)
    ) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_point(
            position = ggplot2::position_jitter(
                width = 0.2,
                height = 0,
                seed = 1L
            ),
            alpha = 0.2,
            size = 0.7,
            colour = "black"
        ) +
        ggplot2::geom_point(
            data = stars,
            ggplot2::aes(x = cell_line, y = y),
            inherit.aes = FALSE,
            shape = 8,
            colour = "red",
            size = 2
        ) +
        ggplot2::facet_wrap(~symbol, ncol = 7, scales = "free_y") +
        ggplot2::labs(
            title = "Baseline expression of probes identified as recurrent DEGs",
            subtitle = "Red asterisk: at least one BH-adjusted-p-only screening call in that model",
            x = NULL,
            y = "log2(DESeq2-normalized count + 1)"
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                hjust = 1,
                size = 7
            ),
            strip.text = ggplot2::element_text(size = 8),
            legend.position = "none"
        )
}
