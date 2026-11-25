# Single source of truth for manuscript output names. No plotting dependencies.
figure_registry <- function() {
    main_script <- "scripts/render_figures.R"
    sample_script <- "scripts/render_sample_supplementary_figures.R"
    legacy_script <- "scripts/render_legacy_supplementary_figures.R"
    keys <- c(
        "deg",
        "ora",
        "gsea",
        "tgx_ddi",
        "tpod_fdr_main",
        "s1",
        "s2_pca_dmso",
        "s3",
        "s4_threshold_heatmap",
        "s5_recurrence",
        "s6_expression",
        "s7",
        "s8_mean_variance",
        "s9",
        "s10",
        "s11",
        "s12",
        "s13_tgx_signature",
        "deg_filter_burden_supplementary",
        "deg_filter_medians_supplementary",
        "gsea_archived_supplementary",
        "ora_wikipathways_supplementary",
        "tpod_parameter_settings",
        "tpod_fdr",
        "tpod_fdr_supplementary",
        "tpod_robustness"
    )
    legacy <- c(
        "deg_monte_carlo_false_positive_rates",
        "ora_monte_carlo_false_positive_rates",
        "gsea_monte_carlo_false_positive_rates",
        "tgx_ddi_monte_carlo_false_positive_rates",
        "tpod_fdr_main",
        "supplementary_figure_s1",
        "supplementary_figure_s2_pca_dmso",
        "supplementary_figure_s3",
        "supplementary_threshold_heatmap",
        "supplementary_figure_s5_recurrence",
        "supplementary_figure_s6_expression",
        "supplementary_figure_s7",
        "supplementary_figure_s8_mean_variance",
        "supplementary_figure_s9",
        "supplementary_figure_s10",
        "supplementary_figure_s11",
        "supplementary_figure_s12",
        "tgx_ddi_signature_supplementary",
        "deg_filter_burden_supplementary",
        "deg_filter_medians_supplementary",
        "gsea_archived_duplicate_ids_supplementary",
        "ora_wikipathways_supplementary",
        "tpod_fdr_by_parameter_setting",
        "tpod_fdr_grid",
        "tpod_fdr_grid_supplementary",
        "tpod_fdr_robustness"
    )
    number <- c(paste0("Fig_", 2:6), paste0("Fig_S", 1:13), rep("", 8))
    category <- c(
        rep("main", 5),
        rep("supplementary", 13),
        rep("exploratory", 8)
    )
    stem <- ifelse(nzchar(number), number, legacy)
    directory <- c(
        rep("figures", 5),
        rep("figures/supplementary", 13),
        rep("figures/exploratory", 8)
    )
    scripts <- c(
        rep(main_script, 5),
        sample_script,
        "scripts/render_dmso_pca.R",
        sample_script,
        "scripts/deg_thresholds/plot_threshold_sensitivity.R",
        main_script,
        main_script,
        legacy_script,
        main_script,
        rep(legacy_script, 4),
        "scripts/tgx_signature/render_figure.R",
        "scripts/deg_thresholds/plot_filter_burden.R",
        "scripts/deg_thresholds/plot_filter_medians.R",
        rep(main_script, 6)
    )
    note <- rep("", length(keys))
    note[keys == "s1"] <- "Frozen archived VST sample distances."
    note[keys == "s3"] <-
        "Four archived n=4 DMSO trials matched to the historical manuscript raster."
    note[keys == "s7"] <- "Expression level versus archived adjusted-p-only probe recurrence."
    note[keys == "s8_mean_variance"] <-
        "Baseline mean-variance figure. Archived variance versus false positives is S9."
    note[keys == "s9"] <- "Expression variance versus archived adjusted-p-only false positives."
    note[keys == "s10"] <- "Frequency of nonzero historical FINAL DEG trial events."
    note[keys == "s11"] <- "Historical last-drawn plate tiles from nonzero FINAL DEG trial events."
    note[keys == "s12"] <- "Archived top-20 GSEA selection per stratum, including ties."
    registry <- data.frame(
        key = keys,
        number = number,
        category = category,
        stem = stem,
        directory = directory,
        script = scripts,
        selection = keys,
        legacy_stem = legacy,
        note = note,
        stringsAsFactors = FALSE
    )
    drafts <- c(
        "supplementary_threshold_lines_draft",
        "supplementary_threshold_lines_swapped_draft",
        "supplementary_threshold_lines_padj_axis_draft",
        "supplementary_threshold_heatmap_first_draft"
    )
    registry <- rbind(
        registry,
        data.frame(
            key = drafts,
            number = "",
            category = "exploratory",
            stem = drafts,
            directory = "figures/exploratory",
            script = "scripts/deg_thresholds/plot_threshold_sensitivity.R",
            selection = drafts,
            legacy_stem = drafts,
            note = "Unselected threshold draft; render with the threshold-drafts task.",
            stringsAsFactors = FALSE
        )
    )
    registry$extensions <- rep(list(c("png", "pdf")), nrow(registry))
    registry$extensions[[match("supplementary_threshold_heatmap_first_draft", registry$key)]] <- "pdf"
    registry
}

figure_spec <- function(key) {
    registry <- figure_registry()
    if (!is.character(key) || length(key) != 1L || is.na(key)) {
        stop("Figure key must be one non-missing string.")
    }
    # Number and old descriptive stem are accepted for standalone producers.
    i <- which(
        registry$key == key |
            registry$number == key |
            registry$legacy_stem == key
    )
    if (length(i) != 1L) {
        stop("Unknown or ambiguous figure key: ", key)
    }
    registry[i, , drop = FALSE]
}

figure_filename <- function(key, extension = NULL) {
    spec <- figure_spec(key)
    if (is.null(extension)) {
        return(file.path(spec$directory, spec$stem))
    }
    if (
        !is.character(extension) ||
            !length(extension) ||
            anyNA(extension) ||
            any(!extension %in% c("png", "pdf"))
    ) {
        stop("Figure extension must be png or pdf (without a leading dot).")
    }
    file.path(spec$directory, paste0(spec$stem, ".", extension))
}

figure_path <- function(
    key,
    extension = NULL,
    root = project_root(),
    create_dir = TRUE
) {
    path <- file.path(root, figure_filename(key, extension))
    if (create_dir) {
        directory <- unique(dirname(path))
        for (d in directory) {
            if (!dir.exists(d) && !dir.create(d, recursive = TRUE)) {
                stop("Could not create figure directory: ", d)
            }
        }
    }
    path
}

# Preflight all conflicts before moving anything. Identical duplicates may be
# removed, but a different destination is never replaced, even on repeat runs.
organize_figure_outputs <- function(root = project_root(), dry_run = FALSE) {
    registry <- figure_registry()
    files <- list.files(
        file.path(root, "figures"),
        pattern = "\\.(png|pdf)$",
        full.names = TRUE,
        recursive = FALSE
    )
    plan <- data.frame(
        source = files,
        destination = files,
        action = "unchanged",
        stringsAsFactors = FALSE
    )
    for (i in seq_along(files)) {
        ext <- tools::file_ext(files[i])
        stem <- tools::file_path_sans_ext(basename(files[i]))
        matched <- which(registry$legacy_stem == stem | registry$stem == stem)
        destination <- if (length(matched) == 1L) {
            figure_path(registry$key[matched], ext, root, create_dir = FALSE)
        } else {
            file.path(root, "figures", "exploratory", basename(files[i]))
        }
        plan$destination[i] <- destination
        if (identical(files[i], destination)) {
            next
        }
        plan$action[i] <- "move"
        if (file.exists(destination)) {
            hashes <- unname(tools::md5sum(c(files[i], destination)))
            if (anyNA(hashes) || hashes[1] != hashes[2]) {
                stop(
                    "Conflicting figure destination; no files moved: ",
                    destination
                )
            }
            plan$action[i] <- "remove_identical_duplicate"
        }
    }
    targets <- plan$destination[plan$action != "unchanged"]
    if (anyDuplicated(targets)) {
        stop(
            "Multiple legacy figures target the same destination; no files moved."
        )
    }
    if (dry_run) {
        return(plan)
    }
    for (i in seq_len(nrow(plan))) {
        if (plan$action[i] == "unchanged") {
            next
        }
        dir.create(
            dirname(plan$destination[i]),
            recursive = TRUE,
            showWarnings = FALSE
        )
        # Recheck immediately before the operation as well as during preflight.
        if (file.exists(plan$destination[i])) {
            hashes <- unname(tools::md5sum(c(
                plan$source[i],
                plan$destination[i]
            )))
            if (anyNA(hashes) || hashes[1] != hashes[2]) {
                stop("Conflicting figure destination: ", plan$destination[i])
            }
            if (!file.remove(plan$source[i])) {
                stop("Could not remove identical duplicate: ", plan$source[i])
            }
        } else if (!file.rename(plan$source[i], plan$destination[i])) {
            stop("Could not move figure: ", plan$source[i])
        }
    }
    plan
}

write_figure_manifest <- function(root = project_root()) {
    registry <- figure_registry()
    rows <- do.call(
        rbind,
        lapply(seq_len(nrow(registry)), function(i) {
            spec <- registry[i, ]
            filenames <- figure_filename(spec$key, spec$extensions[[1L]])
            data.frame(
                number = spec$number,
                category = spec$category,
                filename = filenames,
                script = spec$script,
                selection = spec$selection,
                status = ifelse(
                    file.exists(file.path(root, filenames)),
                    "present",
                    "pending"
                ),
                script_status = ifelse(
                    file.exists(file.path(root, spec$script)),
                    "available",
                    "pending"
                ),
                note = spec$note,
                stringsAsFactors = FALSE
            )
        })
    )
    rows <- rbind(
        data.frame(
            number = "Fig_1",
            category = "main",
            filename = "",
            script = "",
            selection = "",
            status = "not_provided",
            script_status = "not_provided",
            note = "Figure and source were not supplied."
        ),
        rows
    )
    # Retain visibility of unregistered root outputs moved to exploratory.
    extra <- list.files(
        file.path(root, "figures", "exploratory"),
        pattern = "\\.(png|pdf)$",
        full.names = FALSE
    )
    extra <- setdiff(paste0("figures/exploratory/", extra), rows$filename)
    if (length(extra)) {
        rows <- rbind(
            rows,
            data.frame(
                number = "",
                category = "exploratory",
                filename = extra,
                script = "",
                selection = "",
                status = "present",
                script_status = "unknown",
                note = "Unregistered legacy output; producer and manuscript selection unverified."
            )
        )
    }
    path <- file.path(root, "docs", "figure_manifest.tsv")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(
        rows,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
    )
    invisible(rows)
}
