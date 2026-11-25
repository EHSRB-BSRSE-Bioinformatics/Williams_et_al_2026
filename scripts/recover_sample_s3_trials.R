# Audit the earlier n=3 investigation. Publication S3 uses the separate draft raster evidence.
script <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/sample_supplementary_figures.R")
original <- Sys.getenv("FDR_ORIGINAL_ANALYSIS", "../2023_Meier_FDR")
ref <- "data/reference/sample_supplementary"
inputs <- sample_supp_read_inputs()
frozen <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
archive_dir <- "bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10"
archived_md <- readRDS(file.path(original, archive_dir, "ls_metadata.RDS"))
stems <- c("spheroids", "heparg", "mcf7", "tk6")
tables <- sources <- vector("list", 4L)
for (i in seq_along(stems)) {
    model <- sample_supp_models[i]
    path <- file.path(original, archive_dir, paste0(stems[i], ".RDS"))
    cache_path <- file.path(
        "data/large/gene_recurrence",
        paste0(stems[i], ".rds")
    )
    cache <- readRDS(cache_path)
    sha <- sample_supp_sha256(path)
    sample_supp_assert(
        identical(cache$sha256, sha),
        "Screening cache/archive hash mismatch"
    )
    archive <- readRDS(path)
    ids <- names(frozen)[startsWith(names(frozen), paste0(model, "_3_3_"))]
    ids <- ids[order(as.integer(sub(".*_", "", ids)))]
    sample_supp_assert(
        length(ids) == 1000L && all(ids %in% names(archive)),
        "Incomplete frozen n=3 archive"
    )
    defs <- sample_supp_panels[sample_supp_panels$cell_line == model, ]
    tables[[i]] <- do.call(
        rbind,
        lapply(ids, function(id) {
            x <- archive[[id]]
            md <- frozen[[id]]
            sample_supp_assert(
                identical(md$sample_id, archived_md[[id]]$sample_id) &&
                    identical(
                        as.character(md$group),
                        as.character(archived_md[[id]]$group)
                    ) &&
                    identical(md$sample_id, x$exp_metadata$sample_id) &&
                    identical(
                        as.character(md$group),
                        as.character(x$exp_metadata$group)
                    ),
                "Archived result/frozen sample assignment mismatch"
            )
            # Unfiltered DESeq2 results, not FINAL or fold-change-filtered calls.
            res <- as.data.frame(x$DESeqRes)
            sample_supp_assert(
                all(c("padj", "log2FoldChange") %in% names(res)),
                paste("Missing unfiltered results:", id)
            )
            suffix <- sub(".*_", "", defs$feature)
            row <- match(suffix, sub(".*_", "", rownames(res)))
            raw_row <- match(suffix, sub(".*_", "", rownames(inputs$counts)))
            sample_supp_assert(
                !anyNA(raw_row),
                "Missing requested raw-count probe"
            )
            raw <- inputs$counts[raw_row, md$sample_id, drop = FALSE]
            padj <- res$padj[row]
            called <- !is.na(padj) & padj < .05
            cached <- suffix %in% sub(".*_", "", cache$screen[[id]])
            sample_supp_assert(
                identical(called, cached),
                "Numeric result/cache calls differ"
            )
            data.frame(
                cell_line = model,
                trial_name = id,
                panel = defs$panel,
                feature = defs$feature,
                archive_feature = rownames(res)[row],
                padj = padj,
                log2FoldChange = res$log2FoldChange[row],
                mean_A = rowMeans(raw[, md$group == "A", drop = FALSE]),
                mean_B = rowMeans(raw[, md$group == "B", drop = FALSE]),
                min_A = apply(raw[, md$group == "A", drop = FALSE], 1, min),
                max_A = apply(raw[, md$group == "A", drop = FALSE], 1, max),
                min_B = apply(raw[, md$group == "B", drop = FALSE], 1, min),
                max_B = apply(raw[, md$group == "B", drop = FALSE], 1, max),
                called = called,
                stringsAsFactors = FALSE
            )
        })
    )
    sources[[i]] <- data.frame(
        role = "s3_recovery_source",
        path = paste0("legacy/", archive_dir, "/", stems[i], ".RDS"),
        sha256 = sha
    )
    cat(model, ": verified 1000 frozen n=3 trials\n")
    rm(archive, cache)
    gc()
}
scan <- do.call(rbind, tables)
sample_supp_write_tsv(scan, file.path(ref, "s3_n3_probe_results.tsv"))
evidence_paths <- c(".Rhistory", ".RData", "scripts/show_FP_DEGs.R")
sources <- rbind(
    do.call(rbind, sources),
    data.frame(
        role = "s3_recovery_source",
        path = paste0("legacy/", evidence_paths),
        sha256 = vapply(
            file.path(original, evidence_paths),
            sample_supp_sha256,
            ""
        )
    )
)
sample_supp_write_tsv(sources, file.path(ref, "s3_recovery_sources.tsv"))
cat(
    "Exported the historical n=3 probe scan for comparison with the recovered n=4 examples.\n"
)
