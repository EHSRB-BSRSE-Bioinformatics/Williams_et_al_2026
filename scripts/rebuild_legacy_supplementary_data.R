script <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/legacy_supplementary_figures.R")
args <- commandArgs(TRUE)
extract <- "--extract-reference" %in% args
expression_only <- "--extract-expression-reference" %in% args
source_arg <- grep("^--source-dir=", args, value = TRUE)
if ("--source-dir" %in% args) {
    i <- match("--source-dir", args)
    legacy_assert(i < length(args), "--source-dir needs a directory")
    source_arg <- paste0("--source-dir=", args[i + 1L])
    args <- args[-c(i, i + 1L)]
}
legacy_assert(
    all(
        args %in%
            c("--extract-reference", "--extract-expression-reference") |
            grepl("^--source-dir=", args)
    ) &&
        length(source_arg) <= 1L,
    "Usage: Rscript scripts/rebuild_legacy_supplementary_data.R [--extract-reference|--extract-expression-reference --source-dir DIR]"
)
legacy_assert(!extract || !expression_only, "Choose one extraction mode")
legacy_assert(
    extract || expression_only || length(source_arg) == 0L,
    "--source-dir is only used for extraction"
)
ref <- legacy_reference
paths <- file.path(
    ref,
    c(
        "screening_counts.rds",
        "expression_reference.rds",
        "final_trial_counts.rds",
        "gsea_counts.rds",
        "pathway_names.tsv"
    )
)
inputs <- c(
    "data/input/count_table.tsv",
    "data/input/metadata.QC_applied.txt",
    "data/derived/monte_carlo_trial_metadata.rds"
)
counts <- as.matrix(read.delim(inputs[1L], row.names = 1L, check.names = FALSE))
metadata <- read.delim(inputs[2L], check.names = FALSE)
trials <- readRDS(inputs[3L])
source_dir <- if (length(source_arg)) {
    sub("^--source-dir=", "", source_arg)
} else {
    "data/large"
}
if (extract) {
    provenance <- legacy_provenance(inputs, "public_input")
    record_source <- function(path, relative, role) {
        # Source paths are relative to the supplied archive root, never machine paths.
        provenance <<- rbind(
            provenance,
            data.frame(
                role = role,
                path = relative,
                sha256 = legacy_sha256(path)
            )
        )
    }
    dir.create(ref, recursive = TRUE, showWarnings = FALSE)
    screen <- finals <- list()
    expression_features <- NULL
    for (i in seq_along(legacy_models)) {
        model <- legacy_models[i]
        ids <- names(trials)[startsWith(names(trials), paste0(model, "_"))]
        info <- legacy_trial_info(ids)
        cache_path <- file.path(
            "data/large/gene_recurrence",
            paste0(legacy_stems[i], ".rds")
        )
        relative <- file.path(
            "bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10",
            paste0(legacy_stems[i], ".RDS")
        )
        archive <- file.path(source_dir, relative)
        if (file.exists(cache_path)) {
            cache <- readRDS(cache_path)
            legacy_assert(
                cache$cell == model &&
                    setequal(names(cache$screen), ids) &&
                    setequal(names(cache$final), ids),
                "Recurrence cache trial mismatch"
            )
            if (file.exists(archive)) {
                legacy_assert(
                    identical(legacy_sha256(archive), cache$sha256),
                    "DEG archive/cache SHA mismatch"
                )
            }
            provenance <- rbind(
                provenance,
                legacy_provenance(cache_path, "current_DEG_cache"),
                data.frame(
                    role = "historical_DEG_archive",
                    path = relative,
                    sha256 = cache$sha256
                )
            )
            screening <- cache$screen[ids]
            final <- lengths(cache$final[ids])
            # All four archived DESeqRes universes must agree with ddsList$all.
            if (is.null(expression_features)) {
                expression_features <- rownames(cache$tested)
            }
            legacy_assert(
                setequal(expression_features, rownames(cache$tested)),
                "Different DEG feature universes"
            )
            rm(cache)
        } else {
            legacy_assert(
                file.exists(archive),
                "No current DEG cache or historical archive; supply --source-dir"
            )
            raw <- readRDS(archive)
            legacy_assert(
                length(raw) == 8000L,
                "DEG archive must contain 8000 ordered trials"
            )
            if (!is.null(names(raw))) {
                legacy_assert(
                    identical(names(raw), ids),
                    "DEG archive order differs from frozen design"
                )
            }
            screening <- lapply(raw, function(t) {
                d <- as.data.frame(t$DESeqRes)
                rownames(d)[!is.na(d$padj) & d$padj < 0.05]
            })
            # Historical count definitions are authoritative; do not substitute a new FINAL rule.
            final <- NULL
            if (is.null(expression_features)) {
                expression_features <- rownames(raw[[1L]]$DESeqRes)
            }
            record_source(archive, relative, "historical_DEG_archive")
            rm(raw)
        }
        screen[[i]] <- do.call(
            rbind,
            lapply(3:10, function(n) {
                tab <- table(unlist(screening[info$n == n], use.names = FALSE))
                data.frame(
                    cell_line = model,
                    n = n,
                    feature = names(tab),
                    count = as.integer(tab)
                )
            })
        )
        compiled_path <- "data/derived/deg_monte_carlo_summary.rds"
        compiled <- as.data.frame(readRDS(compiled_path))
        keep <- !is.na(compiled$Filter) &
            as.character(compiled$Filter) == "Final DEGs"
        compiled <- compiled[keep, ]
        legacy_assert(
            !anyDuplicated(compiled$name) && all(ids %in% compiled$name),
            "Missing historical FINAL counts"
        )
        observed <- compiled$nDEGs[match(ids, compiled$name)]
        if (!is.null(final)) {
            legacy_assert(
                identical(as.integer(final), as.integer(observed)),
                "FINAL cache/compiled mismatch"
            )
        }
        finals[[i]] <- data.frame(
            trial_id = ids,
            n_final = as.integer(observed)
        )
        rm(screening)
        gc(FALSE)
    }
    provenance <- rbind(
        provenance,
        legacy_provenance(compiled_path, "historical_FINAL_summary")
    )
    saveRDS(do.call(rbind, screen), paths[1L], compress = "xz")
    factor_path <- "data/reference/s56/size_factors.tsv"
    # Existing export has already been verified against archived ddsList$all.
    factor_manifest <- read.delim("data/reference/s56/source_provenance.tsv")
    expected <- factor_manifest$sha256[factor_manifest$path == factor_path]
    legacy_assert(
        length(expected) == 1L &&
            identical(expected, legacy_sha256(factor_path)),
        "Archived size-factor checksum mismatch"
    )
    expression_export <- legacy_extract_expression(source_dir, counts, metadata)
    saveRDS(expression_export$reference, paths[2L], compress = "xz")
    provenance <- rbind(provenance, expression_export$provenance)
    saveRDS(do.call(rbind, finals), paths[3L], compress = "xz")
    provenance <- rbind(
        provenance,
        legacy_provenance(
            c(factor_path, "data/reference/s56/source_provenance.tsv"),
            "archived_expression_export"
        )
    )
    workspace <- factor_manifest[
        factor_manifest$role == "source_expression_workspace",
    ]
    if (nrow(workspace)) {
        workspace$path <- "analysis.cell_lines_comparison/RData/FDR Project_DEG_data.RData"
        provenance <- rbind(provenance, workspace)
    }
    annotation_path <- "data/large/reactome_ora_v1/workspace_provenance/p.anno.rds"
    annotation <- readRDS(annotation_path)
    annotation <- unique(data.frame(
        pathway = annotation$DB_ID,
        path_name = annotation$path_name
    ))
    provenance <- rbind(
        provenance,
        legacy_provenance(annotation_path, "historical_workspace_pathway_names")
    )
    gsea <- vector("list", 4L)
    for (i in seq_along(legacy_models)) {
        relative <- file.path(
            "bootstrap_GSEA/2023-11-08",
            paste0(legacy_stems[i], "_gsea.RDS")
        )
        archive <- file.path(source_dir, relative)
        legacy_assert(
            file.exists(archive),
            "Missing unfiltered historical GSEA archive; supply --source-dir"
        )
        message("Extracting historical GSEA: ", legacy_models[i])
        raw <- readRDS(archive)
        ids <- names(trials)[startsWith(
            names(trials),
            paste0(legacy_models[i], "_")
        )]
        legacy_assert(
            length(raw) == 8000L,
            "GSEA archive must contain 8000 ordered trials"
        )
        if (!is.null(names(raw))) {
            legacy_assert(
                identical(names(raw), ids),
                "GSEA archive order mismatch"
            )
        }
        calls <- lapply(raw, function(d) {
            d <- as.data.frame(d)
            legacy_assert(
                all(c("pathway", "padj", "NES") %in% names(d)) &&
                    !anyDuplicated(d$pathway),
                "Invalid raw GSEA schema"
            )
            d$pathway[!is.na(d$padj) & d$padj < 0.05]
        })
        gsea[[i]] <- do.call(
            rbind,
            lapply(3:10, function(n) {
                tab <- table(unlist(
                    calls[((n - 3L) * 1000L + 1L):((n - 2L) * 1000L)],
                    use.names = FALSE
                ))
                data.frame(
                    cell_line = legacy_models[i],
                    n = n,
                    pathway = names(tab),
                    count = as.integer(tab)
                )
            })
        )
        record_source(archive, relative, "unfiltered_historical_GSEA_archive")
        rm(raw, calls)
        gc(FALSE)
    }
    gsea <- do.call(rbind, gsea)
    annotation <- annotation[annotation$pathway %in% gsea$pathway, ]
    legacy_assert(
        all(gsea$pathway %in% annotation$pathway),
        "Historical GSEA annotation incomplete"
    )
    saveRDS(gsea, paths[4L], compress = "xz")
    legacy_write_tsv(annotation, paths[5L])
    for (rmd in c(
        "Process bootstrapped results.Rmd",
        "Process bootstrapped results GSEA.Rmd"
    )) {
        relative <- file.path("Rmd", rmd)
        if (file.exists(file.path(source_dir, relative))) {
            record_source(
                file.path(source_dir, relative),
                relative,
                "historical_plot_source"
            )
        }
    }
    provenance <- rbind(
        provenance,
        legacy_provenance(paths, "reference_output")
    )
    legacy_write_tsv(provenance, file.path(ref, "source_provenance.tsv"))
}
manifest_path <- file.path(ref, "source_provenance.tsv")
if (expression_only) {
    # Refresh just the independently sourced expression universe, without loading GSEA again.
    manifest <- read.delim(manifest_path)
    expression_export <- legacy_extract_expression(source_dir, counts, metadata)
    saveRDS(expression_export$reference, paths[2L], compress = "xz")
    manifest <- manifest[
        manifest$path != paths[2L] &
            manifest$role != "historical_expression_universe_and_values",
    ]
    manifest <- rbind(
        manifest,
        expression_export$provenance,
        legacy_provenance(paths[2L], "reference_output")
    )
    manifest$path[manifest$role == "source_expression_workspace"] <-
        "analysis.cell_lines_comparison/RData/FDR Project_DEG_data.RData"
    legacy_write_tsv(manifest, manifest_path)
}
manifest <- read.delim(manifest_path)
for (path in c(inputs, paths)) {
    expected <- manifest$sha256[manifest$path == path]
    legacy_assert(
        length(expected) == 1L && identical(expected, legacy_sha256(path)),
        paste("Reference checksum mismatch:", path)
    )
}
expr <- readRDS(paths[2L])
refs <- list(
    screening = readRDS(paths[1L]),
    features = expr$features,
    factors = expr$factors,
    final = readRDS(paths[3L]),
    gsea = readRDS(paths[4L]),
    annotation = read.delim(paths[5L])
)
summary <- legacy_build(counts, metadata, trials, refs)
legacy_validate(summary)
saveRDS(summary, legacy_summary_path, compress = "xz")
outputs <- c(
    legacy_summary_path,
    "data/derived/legacy_supplementary_samples.tsv",
    "data/derived/legacy_supplementary_plate.tsv",
    "data/derived/legacy_supplementary_gsea_selected.tsv",
    "data/derived/legacy_supplementary_definitions.tsv"
)
legacy_write_tsv(summary$samples, outputs[2L])
legacy_write_tsv(summary$plate, outputs[3L])
legacy_write_tsv(summary$pathways$observed, outputs[4L])
legacy_write_tsv(
    data.frame(
        endpoint = names(legacy_definitions),
        definition = unname(legacy_definitions)
    ),
    outputs[5L]
)
legacy_write_tsv(
    rbind(
        legacy_provenance(c(
            inputs,
            paths,
            manifest_path,
            "R/legacy_supplementary_figures.R",
            "scripts/rebuild_legacy_supplementary_data.R"
        )),
        legacy_provenance(outputs, "output")
    ),
    "data/derived/legacy_supplementary_data_provenance.tsv"
)
message(
    "Rebuilt S7/S9-S12: ",
    nrow(summary$scatter),
    " scatter rows; ",
    nrow(summary$samples),
    " sample/n rows; ",
    length(summary$pathways$order),
    " selected pathway names"
)
