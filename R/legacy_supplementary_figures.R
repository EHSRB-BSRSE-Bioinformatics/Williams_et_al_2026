# Standalone historical S7/S9-S12 port. No dependency on the shared figure registry.
legacy_models <- c("TK6", "MCF-7", "HepaRG", "Liver spheroids")
legacy_stems <- c("tk6", "mcf7", "heparg", "spheroids")
legacy_reference <- "data/reference/legacy_supplementary"
legacy_summary_path <- "data/derived/legacy_supplementary_summary.rds"
legacy_definitions <- c(
    probe = "archived DESeqRes: non-missing padj < 0.05; probe IDs, not FINAL genes",
    expression = "mean and sample variance of log2(raw / archived ddsList$all size factor + 1), DMSO only",
    scatter = "left join recurrence by probe/model after dropping n; positive counts remain separate by n; never-recurrent probes have one zero row",
    sample = "sample/model/n: trials with historical Final DEGs > 0 divided by all trials containing sample, either randomized group",
    histogram = "30 bins; each sample/n is an observation (not one pooled observation per sample)",
    plate = "historical tile overplot: descending n_high_DEG, stable sample/n/model tie order; final drawn value per model/row/column",
    gsea = "archived GSEA padj < 0.05; no post-hoc size, NES or unique-gene-ID filter; top 20 rows per model/n with ties, union of pathway names"
)
legacy_assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
}
legacy_sha256 <- function(path) {
    legacy_assert(file.exists(path), paste("Missing input:", path))
    x <- system2("sha256sum", shQuote(path), stdout = TRUE)
    legacy_assert(
        length(x) == 1L && grepl("^[[:xdigit:]]{64} ", x),
        "SHA-256 failed"
    )
    sub(" .*", "", x)
}
legacy_write_tsv <- function(x, path) {
    old <- options(digits = 17)
    on.exit(options(old))
    write.table(
        x,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}
legacy_provenance <- function(paths, role = "input_or_code") {
    data.frame(
        role = rep(role, length(paths)),
        path = paths,
        sha256 = unname(vapply(paths, legacy_sha256, ""))
    )
}
legacy_extract_expression <- function(source_dir, counts, metadata) {
    relative <- paste0(
        "analysis.cell_lines_comparison/DEG_lists/",
        "TempO-Seq_FDR_Project_all_05-05-2023.10.12-Per_sample_normalized_counts.txt"
    )
    path <- file.path(source_dir, relative)
    legacy_assert(
        file.exists(path),
        "Missing historical dds-all normalized-count export"
    )
    archived <- as.matrix(read.delim(path, row.names = 1L, check.names = FALSE))
    factor_path <- "data/reference/s56/size_factors.tsv"
    manifest <- read.delim("data/reference/s56/source_provenance.tsv")
    expected <- manifest$sha256[manifest$path == factor_path]
    legacy_assert(
        length(expected) == 1L &&
            identical(expected, legacy_sha256(factor_path)),
        "Archived size-factor checksum mismatch"
    )
    factors <- read.delim(factor_path)
    features <- rownames(archived)
    legacy_assert(
        nrow(archived) == 22460L &&
            ncol(archived) == 365L &&
            all(features %in% rownames(counts)) &&
            all(factors$sample_id %in% colnames(archived)),
        "Unexpected historical expression universe"
    )
    observed <- sweep(
        counts[features, factors$sample_id, drop = FALSE],
        2L,
        factors$size_factor,
        "/"
    )
    expected <- archived[, factors$sample_id, drop = FALSE]
    # Historical text export rounded normalized counts to five decimal places.
    legacy_assert(
        max(abs(expected - observed)) < 1e-4,
        "Archived normalized-count export differs from public reconstruction"
    )
    list(
        reference = list(features = features, factors = factors),
        provenance = data.frame(
            role = "historical_expression_universe_and_values",
            path = relative,
            sha256 = legacy_sha256(path)
        )
    )
}
legacy_trial_info <- function(ids) {
    pieces <- strsplit(ids, "_", fixed = TRUE)
    legacy_assert(all(lengths(pieces) == 4L), "Malformed trial identifiers")
    x <- data.frame(
        trial_id = ids,
        cell_line = vapply(pieces, `[`, "", 1L),
        n = as.integer(vapply(pieces, `[`, "", 2L)),
        n_b = as.integer(vapply(pieces, `[`, "", 3L)),
        iteration = as.integer(vapply(pieces, `[`, "", 4L))
    )
    legacy_assert(
        !anyNA(x) &&
            !anyDuplicated(ids) &&
            all(x$n == x$n_b) &&
            all(x$n %in% 3:10) &&
            all(x$cell_line %in% legacy_models),
        "Invalid trial design"
    )
    x
}
legacy_validate_counts <- function(x, id) {
    legacy_assert(
        all(c("cell_line", "n", id, "count") %in% names(x)) &&
            !anyNA(x) &&
            !anyDuplicated(x[c("cell_line", "n", id)]) &&
            all(x$count > 0 & x$count <= 1000) &&
            all(x$count == as.integer(x$count)) &&
            all(x$cell_line %in% legacy_models) &&
            all(x$n %in% 3:10),
        paste("Invalid recurrence table:", id)
    )
}
legacy_expression <- function(counts, metadata, factors, features) {
    legacy_assert(
        !anyDuplicated(rownames(counts)) &&
            !anyDuplicated(colnames(counts)) &&
            !anyDuplicated(metadata$Sample_ID) &&
            !anyDuplicated(factors$sample_id) &&
            !anyDuplicated(features) &&
            all(features %in% rownames(counts)) &&
            all(factors$sample_id %in% colnames(counts)) &&
            all(is.finite(factors$size_factor) & factors$size_factor > 0),
        "Invalid expression identifiers or factors"
    )
    md <- metadata[match(factors$sample_id, metadata$Sample_ID), ]
    legacy_assert(
        !anyNA(md$Sample_ID) &&
            all(md$chemical == "DMSO") &&
            setequal(md$cell_line, legacy_models),
        "Expression sample mismatch"
    )
    raw <- counts[features, factors$sample_id, drop = FALSE]
    legacy_assert(all(is.finite(raw) & raw >= 0), "Invalid public counts")
    z <- log2(sweep(raw, 2L, factors$size_factor, "/") + 1)
    do.call(
        rbind,
        lapply(legacy_models, function(model) {
            v <- z[, md$cell_line == model, drop = FALSE]
            legacy_assert(ncol(v) > 1L, "Too few DMSO samples")
            data.frame(
                feature = features,
                cell_line = model,
                mean_expression = rowMeans(v),
                variance = apply(v, 1L, stats::var)
            )
        })
    )
}
legacy_scatter <- function(expression, screening) {
    legacy_validate_counts(screening, "feature")
    # Keep n for auditability, though the historical join selected it away.
    x <- merge(
        expression,
        screening,
        by = c("cell_line", "feature"),
        all.x = TRUE,
        sort = FALSE
    )
    x$count[is.na(x$count)] <- 0L
    x <- x[
        order(
            match(x$cell_line, legacy_models),
            x$feature,
            x$n,
            na.last = TRUE
        ),
    ]
    rownames(x) <- NULL
    x
}
legacy_sample_participation <- function(trials, final, metadata) {
    info <- legacy_trial_info(names(trials))
    legacy_assert(
        nrow(info) == 32000L &&
            all(table(info$cell_line, info$n) == 1000L) &&
            !anyDuplicated(final$trial_id) &&
            setequal(final$trial_id, info$trial_id) &&
            !anyNA(final$n_final) &&
            all(final$n_final >= 0) &&
            all(final$n_final == as.integer(final$n_final)),
        "Incomplete FINAL trial design"
    )
    final <- final[match(info$trial_id, final$trial_id), ]
    parts <- lapply(seq_along(trials), function(i) {
        t <- trials[[i]]
        legacy_assert(
            all(c("sample_id", "group") %in% names(t)) &&
                !anyNA(t$sample_id) &&
                !anyDuplicated(t$sample_id) &&
                nrow(t) == 2L * info$n[i] &&
                setequal(as.character(t$group), c("A", "B")) &&
                all(table(t$group) == info$n[i]),
            "Invalid trial samples/groups"
        )
        data.frame(
            sample_id = t$sample_id,
            cell_line = info$cell_line[i],
            n = info$n[i],
            n_high_DEG = as.integer(final$n_final[i] > 0),
            total_experiments = 1L
        )
    })
    x <- do.call(rbind, parts)
    md <- metadata[match(x$sample_id, metadata$Sample_ID), ]
    legacy_assert(
        !anyDuplicated(metadata$Sample_ID) &&
            !anyNA(md$Sample_ID) &&
            all(md$chemical == "DMSO") &&
            all(md$cell_line == x$cell_line),
        "Trial samples do not match public DMSO metadata"
    )
    x <- aggregate(
        cbind(n_high_DEG, total_experiments) ~ sample_id + n + cell_line,
        x,
        sum
    )
    # Historical group_by sorts keys, then arrange(desc(n_high_DEG)) is stable.
    x <- x[order(x$sample_id, x$n, x$cell_line), ]
    x <- x[order(-x$n_high_DEG, method = "radix"), ]
    x$high_DEG_proportion <- x$n_high_DEG / x$total_experiments
    md <- metadata[match(x$sample_id, metadata$Sample_ID), ]
    x$row <- as.character(md$row)
    x$column <- as.integer(md$column)
    x$batch <- as.character(md$batch)
    legacy_assert(
        !anyNA(x) && all(x$row %in% LETTERS[1:8]) && all(x$column %in% 1:12),
        "Missing or invalid plate positions"
    )
    rownames(x) <- NULL
    x
}
legacy_plate <- function(samples) {
    key <- paste(samples$cell_line, samples$row, samples$column)
    x <- samples[!duplicated(key, fromLast = TRUE), ]
    x$overplotted_rows <- as.integer(table(key)[key[
        !duplicated(key, fromLast = TRUE)
    ]])
    x <- x[order(match(x$cell_line, legacy_models), x$row, x$column), ]
    rownames(x) <- NULL
    x
}
legacy_pathways <- function(counts, annotation) {
    legacy_validate_counts(counts, "pathway")
    legacy_assert(
        all(c("pathway", "path_name") %in% names(annotation)) &&
            !anyDuplicated(annotation$pathway),
        "Invalid pathway annotation"
    )
    x <- merge(counts, annotation, by = "pathway", all.x = TRUE)
    legacy_assert(!anyNA(x$path_name), "Missing historical pathway names")
    x$path_name <- sub("^Homo sapiens: ", "", x$path_name)
    # The Rmd title says top 10, but executable selection is top_n(20).
    selected <- unique(unlist(
        lapply(split(x, list(x$cell_line, x$n), drop = TRUE), function(d) {
            cutoff <- sort(d$count, decreasing = TRUE)[min(20L, nrow(d))]
            d$path_name[d$count >= cutoff]
        }),
        use.names = FALSE
    ))
    observed <- x[x$path_name %in% selected, ]
    means <- tapply(observed$count, observed$path_name, mean)
    path_order <- names(means)[order(means, names(means))]
    # Grey for cells absent from historical positive-count table, not inferred eligibility.
    grid <- expand.grid(
        path_name = path_order,
        n = 3:10,
        cell_line = legacy_models,
        stringsAsFactors = FALSE
    )
    legacy_assert(
        !anyDuplicated(observed[c("cell_line", "n", "path_name")]),
        "Multiple pathway IDs share a selected name; do not silently combine"
    )
    grid <- merge(
        grid,
        observed,
        by = c("cell_line", "n", "path_name"),
        all.x = TRUE
    )
    list(
        selected = selected,
        order = path_order,
        observed = observed,
        grid = grid
    )
}
legacy_build <- function(counts, metadata, trials, refs) {
    expression <- legacy_expression(
        counts,
        metadata,
        refs$factors,
        refs$features
    )
    samples <- legacy_sample_participation(trials, refs$final, metadata)
    list(
        schema_version = 1L,
        definitions = legacy_definitions,
        expression = expression,
        scatter = legacy_scatter(expression, refs$screening),
        samples = samples,
        plate = legacy_plate(samples),
        pathways = legacy_pathways(refs$gsea, refs$annotation)
    )
}
legacy_validate <- function(x) {
    legacy_assert(
        identical(x$schema_version, 1L) &&
            identical(x$definitions, legacy_definitions),
        "Wrong legacy schema/definitions"
    )
    legacy_assert(
        !anyNA(x$expression) &&
            all(is.finite(x$expression$variance)) &&
            all(x$expression$variance >= 0) &&
            nrow(x$expression) == 4L * length(unique(x$expression$feature)),
        "Invalid expression summary"
    )
    legacy_assert(
        nrow(x$scatter) >= nrow(x$expression) &&
            !anyNA(x$scatter$count) &&
            all(x$scatter$count %in% 0:1000) &&
            all(is.na(x$scatter$n) == (x$scatter$count == 0)),
        "Invalid historical scatter join"
    )
    s <- x$samples
    legacy_assert(
        !anyDuplicated(s[c("sample_id", "cell_line", "n")]) &&
            all(s$high_DEG_proportion >= 0 & s$high_DEG_proportion <= 1) &&
            all(s$n_high_DEG <= s$total_experiments) &&
            sum(s$total_experiments) == 4L * 1000L * sum(2L * (3:10)) &&
            nrow(s) == 8L * length(unique(s$sample_id)),
        "Invalid sample participation"
    )
    legacy_assert(
        identical(x$plate, legacy_plate(s)),
        "Wrong historical plate overplot"
    )
    g <- x$pathways$grid
    legacy_assert(
        nrow(g) == 32L * length(x$pathways$order) &&
            all(g$count[!is.na(g$count)] %in% 1:1000),
        "Invalid pathway heatmap"
    )
    invisible(TRUE)
}
legacy_plots <- function(x) {
    legacy_validate(x)
    suppressPackageStartupMessages(library(ggplot2))
    colors <- setNames(
        c("#2166AC", "#B2182B", "#1B7837", "#762A83"),
        legacy_models
    )
    theme <- theme_minimal(base_size = 12) +
        theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            strip.text = element_text(face = "bold"),
            plot.subtitle = element_text(size = 10),
            legend.position = "bottom"
        )
    s <- x$scatter
    s$cell_line <- factor(s$cell_line, levels = legacy_models)
    scatter <- function(column, label, number) {
        ggplot(s, aes(x = .data[[column]], y = count)) +
            geom_point(aes(color = cell_line), alpha = 0.2, size = 1) +
            geom_smooth(
                method = "lm",
                formula = y ~ x,
                se = TRUE,
                color = "blue"
            ) +
            facet_wrap(~cell_line, ncol = 2) +
            scale_color_manual(values = colors, guide = "none") +
            labs(
                title = paste(
                    number,
                    if (column == "variance") {
                        "Expression variance and false positives"
                    } else {
                        "Expression level and false positives"
                    }
                ),
                subtitle = "Historical padj-only probe recurrence; separate counts per sample size (not pooled FINAL)",
                x = label,
                y = "Times probe was a false positive"
            ) +
            theme
    }
    samples <- x$samples
    samples$cell_line <- factor(samples$cell_line, levels = legacy_models)
    plate <- x$plate
    plate$cell_line <- factor(plate$cell_line, levels = legacy_models)
    grid <- x$pathways$grid
    grid$cell_line <- factor(grid$cell_line, levels = legacy_models)
    grid$path_name <- factor(grid$path_name, levels = x$pathways$order)
    list(
        S7 = scatter("mean_expression", "Mean log2 normalized count", "S7."),
        S9 = scatter("variance", "Variance of log2 normalized count", "S9."),
        S10 = ggplot(samples, aes(high_DEG_proportion)) +
            geom_histogram(
                bins = 30,
                fill = "#2166AC",
                color = "black",
                linewidth = 0.2
            ) +
            facet_grid(~cell_line) +
            labs(
                title = "S10. Sample participation in trials with non-zero FINAL DEGs",
                subtitle = "Each observation is a sample / sample-size combination; denominator: all trials containing that sample",
                x = "Proportion of participating trials with non-zero FINAL DEGs",
                y = "Sample / sample-size frequency"
            ) +
            theme,
        S11 = ggplot(
            plate,
            aes(row, factor(column), fill = high_DEG_proportion)
        ) +
            geom_tile() +
            facet_grid(~cell_line) +
            scale_fill_viridis_c(name = "Proportion", na.value = "grey85") +
            labs(
                title = "S11. Plate position and non-zero FINAL DEG trials",
                subtitle = "Historical last-drawn sample-size value per well (overplot reproduced; not a pooled plate effect)",
                x = "Row",
                y = "Column"
            ) +
            theme +
            theme(panel.grid = element_blank()),
        S12 = ggplot(grid, aes(factor(n), path_name, fill = count)) +
            geom_tile() +
            facet_grid(~cell_line) +
            scale_fill_viridis_c(
                option = "F",
                na.value = "grey90",
                name = "Trials enriched"
            ) +
            scale_y_discrete(labels = function(z) {
                vapply(
                    z,
                    function(s) paste(strwrap(s, 62), collapse = "\n"),
                    ""
                )
            }) +
            labs(
                title = "S12. Top recurrent historical GSEA pathways",
                subtitle = "Archived padj < 0.05; top 20 per model / sample size with ties; no revised size, NES or unique-ID filter",
                x = "Sample size per group",
                y = "Reactome pathway"
            ) +
            theme +
            theme(
                panel.grid = element_blank(),
                axis.text.y = element_text(size = 9)
            )
    )
}
