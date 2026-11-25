# Standalone S1/S3 implementation: no shared configuration or renderer required.
sample_supp_models <- c("Liver spheroids", "HepaRG", "MCF-7", "TK6")
sample_supp_panels <- data.frame(
    panel = LETTERS[1:8],
    feature = c(
        "KRT7_10561",
        "FLNC_89124",
        "GSTA2_2817",
        "ABCB11_82408",
        "H2AC17_33569",
        "H4C3_87348",
        "FOXD4L3_33551",
        "H3-3A_92932"
    ),
    cell_line = rep(sample_supp_models, each = 2L),
    stringsAsFactors = FALSE
)
sample_supp_assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
}
sample_supp_sha256 <- function(path) {
    sample_supp_assert(file.exists(path), paste("Missing input:", path))
    value <- system2("sha256sum", shQuote(path), stdout = TRUE)
    sample_supp_assert(
        length(value) == 1L &&
            grepl("^[[:xdigit:]]{64} ", value),
        "SHA-256 calculation failed"
    )
    sub(" .*", "", value)
}
sample_supp_write_tsv <- function(x, path) {
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
sample_supp_read_inputs <- function() {
    list(
        counts = as.matrix(read.delim(
            "data/input/count_table.tsv",
            row.names = 1L,
            check.names = FALSE
        )),
        metadata = read.delim(
            "data/input/metadata.QC_applied.txt",
            check.names = FALSE
        ),
        factors = read.delim("data/reference/s56/size_factors.tsv")
    )
}
sample_supp_validate_distance <- function(distance, metadata) {
    sample_supp_assert(
        is.matrix(distance) &&
            is.numeric(distance) &&
            nrow(distance) > 1L &&
            nrow(distance) == ncol(distance) &&
            !is.null(rownames(distance)) &&
            identical(rownames(distance), colnames(distance)) &&
            !anyDuplicated(rownames(distance)) &&
            all(is.finite(distance)) &&
            all(distance >= 0) &&
            max(abs(distance - t(distance))) < 1e-10 &&
            max(abs(diag(distance))) < 1e-10,
        "Invalid frozen distance matrix"
    )
    sample_supp_assert(
        all(
            c("sample_id", "cell_line", "chemical") %in%
                names(metadata)
        ) &&
            !anyNA(metadata) &&
            !anyDuplicated(metadata$sample_id) &&
            identical(metadata$sample_id, rownames(distance)) &&
            setequal(unique(metadata$cell_line), sample_supp_models),
        "S1 distance/annotation identifiers do not match"
    )
    invisible(TRUE)
}
sample_supp_select_s3 <- function(evidence, frozen, counts) {
    sample_supp_assert(
        identical(evidence$feature, sample_supp_panels$feature) &&
            !anyNA(evidence) && !anyDuplicated(names(frozen)),
        "Invalid S3 raster evidence or trial metadata"
    )
    parse_rows <- function(rows) {
        as.numeric(strsplit(rows, ",", fixed = TRUE)[[1L]])
    }
    panel_score <- function(panel, assignments) {
        ticks <- parse_rows(panel$tick_counts)
        tick_rows <- parse_rows(panel$tick_rows)
        log_counts <- log10(ticks)
        fit <- stats::lm(tick_rows ~ log_counts)
        raster <- list(A = parse_rows(panel$a_rows), B = parse_rows(panel$b_rows))
        slope <- unname(stats::coef(fit)[2L])
        intercept <- unname(stats::coef(fit)[1L])
        suffix <- sub(".*_", "_", panel$feature)
        feature <- rownames(counts)[endsWith(rownames(counts), suffix)]
        sample_supp_assert(length(feature) == 1L, "Missing or ambiguous S3 probe suffix")
        raw <- counts[feature, ]
        positive <- raw > 0
        centers <- round(intercept + slope * log10(raw[positive]), 8)
        # At this raster resolution the thresholded markers extend about 3 px.
        # A union of rows retains the shape of merged or overlapping markers.
        projected <- setNames(vector("list", length(raw)), names(raw))
        projected[positive] <- lapply(centers, function(y) {
            seq.int(ceiling(y - 3), floor(y + 3))
        })
        vapply(assignments, function(groups) {
            sum(vapply(names(raster), function(group) {
                rows <- unique(unlist(projected[groups[[group]]], use.names = FALSE))
                observed <- raster[[group]]
                length(rows) + length(observed) - 2L * sum(rows %in% observed)
            }, integer(1)))
        }, integer(1))
    }
    selections <- decisions <- vector("list", length(sample_supp_models))
    for (i in seq_along(sample_supp_models)) {
        model <- sample_supp_models[i]
        ids <- as.vector(outer(3:10, seq_len(1000L), function(n, run) {
            paste(model, n, n, run, sep = "_")
        }))
        sample_supp_assert(all(ids %in% names(frozen)), "Incomplete S3 trial universe")
        trials <- frozen[ids]
        assignments <- lapply(trials, function(trial) {
            split(trial$sample_id, trial$group)
        })
        panels <- which(sample_supp_panels$cell_line == model)
        scores <- vapply(panels, function(panel) {
            panel_score(evidence[panel, ], assignments)
        }, integer(length(ids)))
        total <- rowSums(scores)
        ranked <- order(total)
        best <- ranked[1L]
        second <- ranked[2L]
        trial <- trials[[best]]
        sample_supp_assert(
            all(table(trial$group) == 4L) && length(trial$sample_id) == 8L &&
                all(scores[best, ] <= 5L) && total[second] - total[best] >= 30L,
            paste("S3 raster does not identify a clear n=4 match for", model)
        )
        decisions[[i]] <- data.frame(
            cell_line = model,
            trials_scanned = length(ids),
            trial_name = ids[best],
            left_panel_score = scores[best, 1L],
            right_panel_score = scores[best, 2L],
            pair_score = total[best],
            runner_up_score = total[second],
            row.names = NULL
        )
        selections[[i]] <- data.frame(
            cell_line = model,
            trial_name = ids[best],
            sample_id = trial$sample_id,
            group = as.character(trial$group)
        )
    }
    list(selection = do.call(rbind, selections), decisions = do.call(rbind, decisions))
}
sample_supp_build <- function(distance, annotation, selection, inputs, decisions) {
    sample_supp_validate_distance(distance, annotation)
    counts <- inputs$counts
    md <- inputs$metadata
    sf <- inputs$factors
    sample_supp_assert(
        !anyDuplicated(rownames(counts)) &&
            !anyDuplicated(colnames(counts)) &&
            !anyDuplicated(md$Sample_ID) &&
            !anyDuplicated(sf$sample_id) &&
            all(is.finite(sf$size_factor)) &&
            all(sf$size_factor > 0),
        "Invalid public identifiers or size factors"
    )
    sample_supp_assert(
        all(annotation$sample_id %in% colnames(counts)) &&
            identical(
                as.character(md$cell_line[match(
                    annotation$sample_id,
                    md$Sample_ID
                )]),
                annotation$cell_line
            ) &&
            identical(
                as.character(md$chemical[match(
                    annotation$sample_id,
                    md$Sample_ID
                )]),
                annotation$chemical
            ),
        "Public metadata differ from frozen S1 annotations"
    )
    sample_supp_assert(
        all(
            c("sample_id", "group", "cell_line", "trial_name") %in%
                names(selection)
        ) &&
            !anyNA(selection) &&
            nrow(selection) == 8L * length(sample_supp_models) &&
            !anyDuplicated(selection$sample_id) &&
            setequal(selection$cell_line, sample_supp_models) &&
            all(selection$group %in% c("A", "B")),
        "Invalid S3 selection"
    )
    for (model in sample_supp_models) {
        sel <- selection[selection$cell_line == model, ]
        trial <- decisions$trial_name[decisions$cell_line == model]
        sample_supp_assert(
            identical(as.integer(table(factor(sel$group, levels = c("A", "B")))), c(4L, 4L)) &&
                identical(unique(sel$trial_name), trial),
            "S3 requires the raster-matched n=4 trial with four wells in each group"
        )
        matched <- match(sel$sample_id, md$Sample_ID)
        sample_supp_assert(
            !anyNA(matched) &&
                all(md$chemical[matched] == "DMSO") &&
                all(md$cell_line[matched] == model),
            "S3 selected wells are not model-matched DMSO"
        )
    }
    rows <- lapply(seq_len(nrow(sample_supp_panels)), function(i) {
        def <- sample_supp_panels[i, ]
        # Probe suffix is stable across historical histone-symbol updates.
        suffix <- sub(".*_", "", def$feature)
        feature <- rownames(counts)[grepl(
            paste0("_", suffix, "$"),
            rownames(counts)
        )]
        sample_supp_assert(
            length(feature) == 1L,
            "Missing or ambiguous S3 probe suffix"
        )
        wells <- sort(md$Sample_ID[
            md$cell_line == def$cell_line & md$chemical == "DMSO"
        ])
        factors <- sf$size_factor[match(wells, sf$sample_id)]
        raw <- as.numeric(counts[feature, wells])
        sample_supp_assert(
            !anyNA(factors) &&
                all(wells %in% colnames(counts)) &&
                all(is.finite(raw)) &&
                all(raw >= 0),
            "Incomplete S3 counts/normalization"
        )
        all <- data.frame(
            panel = def$panel,
            feature = def$feature,
            input_feature = feature,
            cell_line = def$cell_line,
            sample_id = wells,
            group = "All DMSO",
            raw_count = raw,
            size_factor = factors,
            normalized_count = raw / factors,
            stringsAsFactors = FALSE
        )
        sel <- selection[selection$cell_line == def$cell_line, ]
        chosen <- all[match(sel$sample_id, all$sample_id), ]
        chosen$group <- sel$group
        chosen$trial_name <- sel$trial_name
        all$trial_name <- sel$trial_name[1L]
        rbind(all, chosen) # Selected wells also remain in the full DMSO population.
    })
    expression <- do.call(rbind, rows)
    rownames(expression) <- NULL
    list(
        schema_version = 3L,
        distance = distance,
        annotation = annotation,
        tree = stats::hclust(stats::as.dist(distance), method = "complete"),
        panels = sample_supp_panels,
        selection = selection,
        expression = expression,
        s3_decisions = decisions,
        s1_rule = "Euclidean distance of frozen rldList$all assay (VST producer)",
        s3_selection = "Archived raster matched to frozen n=4 assignments over both probes per model",
        s3_rule = "Raw counts by default; display log10(count + 1) retains zeros"
    )
}
sample_supp_draw_s1 <- function(summary) {
    sample_supp_validate_distance(summary$distance, summary$annotation)
    tree <- summary$tree
    order <- tree$order
    n <- length(order)
    matrix <- summary$distance[order, order]
    model_colors <- setNames(
        c("#009E73", "#E69F00", "#0072B2", "#CC79A7"),
        sample_supp_models
    )
    # Reversed Blues: low within-model distances are blue, high distances white.
    colors <- grDevices::colorRampPalette(c(
        "#08306B",
        "#2171B5",
        "#6BAED6",
        "#C6DBEF",
        "#F7FBFF"
    ))(255L)
    old <- par(no.readonly = TRUE)
    on.exit({
        layout(1)
        par(old)
    })
    layout(
        matrix(c(0, 1, 0, 0, 2, 0, 3, 4, 5), nrow = 3, byrow = TRUE),
        widths = c(1.3, 7, 1.7),
        heights = c(1.3, .3, 7)
    )
    par(mar = c(0, 0, 3, 0))
    plot(
        as.dendrogram(tree),
        leaflab = "none",
        axes = FALSE,
        main = "S1  Sample-to-sample distances",
        cex.main = 1.1
    )
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(.5, n + .5), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(
        seq_len(n) - .5,
        0,
        seq_len(n) + .5,
        1,
        col = model_colors[summary$annotation$cell_line[order]],
        border = NA
    )
    par(mar = c(0, 0, 0, 0))
    plot(
        as.dendrogram(tree),
        horiz = TRUE,
        leaflab = "none",
        axes = FALSE,
        xlim = c(max(tree$height), 0),
        ylim = c(n + .5, .5)
    )
    par(mar = c(2.5, 0, 0, 0))
    image(
        seq_len(n),
        seq_len(n),
        t(matrix[n:1L, , drop = FALSE]),
        col = colors,
        axes = FALSE,
        xlab = "",
        ylab = "",
        useRaster = TRUE
    )
    box(col = "#888888")
    mtext(
        "365 QC samples: 301 DMSO + 64 Media; frozen VST-assay distances",
        side = 1,
        line = .6,
        cex = .65
    )
    par(mar = c(3, 1, 0, 4))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend(
        "top",
        legend = sample_supp_models,
        fill = model_colors,
        bty = "n",
        cex = .8,
        title = "Cell model"
    )
    maximum <- max(matrix)
    rect(
        .1,
        seq(0, .48, length.out = 256L)[1:255],
        .3,
        seq(0, .48, length.out = 256L)[2:256],
        col = colors,
        border = NA
    )
    ticks <- pretty(c(0, maximum), n = 4)
    ticks <- ticks[ticks >= 0 & ticks <= maximum]
    text(.4, ticks / maximum * .48, labels = round(ticks), adj = 0, cex = .8)
    text(.1, .53, "Euclidean\ndistance", adj = 0, cex = .8)
}
sample_supp_s3_plots <- function(summary, scale = "raw") {
    sample_supp_assert(scale %in% c("normalized", "raw"), "Unknown S3 scale")
    sample_supp_assert(
        requireNamespace("ggplot2", quietly = TRUE),
        "S3 needs ggplot2"
    )
    lapply(seq_len(nrow(summary$panels)), function(i) {
        def <- summary$panels[i, ]
        data <- summary$expression[summary$expression$panel == def$panel, ]
        data$group <- factor(data$group, levels = c("All DMSO", "B", "A"))
        # Deterministic offsets, independent of RNG and device/render order.
        data$x <- ave(seq_len(nrow(data)), data$group, FUN = function(x) {
            seq(-.05, .05, length.out = length(x))
        })
        data$display_count <- if (scale == "raw") {
            data$raw_count
        } else {
            data$normalized_count
        }
        ggplot2::ggplot(
            data,
            ggplot2::aes(x, display_count + 1, color = group)
        ) +
            ggplot2::geom_point(size = 1.75, alpha = .8) +
            ggplot2::facet_grid(
                cols = ggplot2::vars(group),
                drop = FALSE
            ) +
            ggplot2::scale_x_continuous(
                limits = c(-.24, .24),
                expand = ggplot2::expansion(mult = 0)
            ) +
            ggplot2::scale_y_log10() +
            ggplot2::scale_color_manual(
                values = c("All DMSO" = "#F8766D", B = "#00BA38", A = "#619CFF")
            ) +
            ggplot2::labs(
                title = paste0(
                    def$panel,
                    "  ",
                    def$feature,
                    ", ",
                    def$cell_line
                ),
                x = NULL,
                y = if (scale == "raw") {
                    "Raw count + 1 (log10)"
                } else {
                    "Normalized count + 1 (log10)"
                }
            ) +
            ggplot2::theme_bw(base_size = 9) +
            ggplot2::theme(
                legend.position = "none",
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(size = 10),
                panel.grid.minor = ggplot2::element_blank()
            )
    })
}
sample_supp_draw_s3 <- function(plots) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        layout = grid::grid.layout(
            5,
            2,
            heights = grid::unit(c(.18, 1, 1, 1, 1), "null")
        )
    ))
    grid::grid.text(
        "S3  DMSO mock comparisons with four wells per group",
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2),
        gp = grid::gpar(fontsize = 11)
    )
    for (i in seq_along(plots)) {
        print(
            plots[[i]],
            newpage = FALSE,
            vp = grid::viewport(
                layout.pos.row = 1L + ceiling(i / 2),
                layout.pos.col = 1L + (i - 1L) %% 2L
            )
        )
    }
    grid::popViewport()
}
