script <- sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(script) == 1L) {
    setwd(dirname(dirname(normalizePath(script))))
}
source("R/sample_supplementary_figures.R")
args <- commandArgs(TRUE)
sample_supp_assert(
    length(args) == 0L || identical(args, "--extract-reference"),
    "Usage: Rscript scripts/rebuild_sample_supplementary_data.R [--extract-reference]"
)
ref <- "data/reference/sample_supplementary"
dir.create(ref, recursive = TRUE, showWarnings = FALSE)
inputs <- sample_supp_read_inputs()
distance_path <- file.path(ref, "s1_distances.rds")
annotation_path <- file.path(ref, "s1_samples.tsv")
selection_path <- file.path(ref, "s3_trial_selection.tsv")
panel_path <- file.path(ref, "s3_panels.tsv")
evidence_path <- file.path(ref, "s3_raster_evidence.tsv")
raster_sources_path <- file.path(ref, "s3_raster_sources.tsv")
decision_path <- file.path(ref, "s3_selection_decisions.tsv")
frozen_path <- "data/derived/monte_carlo_trial_metadata.rds"
frozen <- readRDS(frozen_path)
evidence <- read.delim(evidence_path)
recovered <- sample_supp_select_s3(evidence, frozen, inputs$counts)
if (identical(args, "--extract-reference")) {
    sample_supp_assert(
        requireNamespace("DESeq2", quietly = TRUE),
        "Reference extraction requires the pixi deg or deg-study environment"
    )
    original <- Sys.getenv("FDR_ORIGINAL_ANALYSIS", "../2023_Meier_FDR")
    archive_relative <- "analysis.cell_lines_comparison/RData/FDR Project_DEG_data.RData"
    archive <- file.path(original, archive_relative)
    env <- new.env()
    load(archive, envir = env)
    transformed <- env$rldList$all
    dds <- env$ddsList$all
    assay <- SummarizedExperiment::assay(transformed)
    sample_supp_assert(
        is.matrix(assay) &&
            all(is.finite(assay)) &&
            identical(dim(assay), c(22460L, 365L)) &&
            identical(colnames(assay), colnames(dds)),
        "Unexpected archived transformed assay"
    )
    public_raw <- inputs$counts[rownames(dds), colnames(dds), drop = FALSE]
    sample_supp_assert(
        isTRUE(all.equal(
            unname(DESeq2::counts(dds)),
            unname(public_raw),
            tolerance = 0
        )),
        "Archive counts differ from public counts"
    )
    sample_supp_assert(
        !isTRUE(all.equal(unname(assay), unname(public_raw))) &&
            !isTRUE(all.equal(unname(assay), unname(log2(public_raw + 1)))),
        "Archive does not contain the expected transformed assay"
    )
    ids <- sort(colnames(assay))
    distance <- as.matrix(stats::dist(t(assay[, ids, drop = FALSE])))
    md <- inputs$metadata[match(ids, inputs$metadata$Sample_ID), ]
    annotation <- data.frame(
        sample_id = ids,
        cell_line = md$cell_line,
        chemical = md$chemical,
        stringsAsFactors = FALSE
    )
    archived_md <- as.data.frame(SummarizedExperiment::colData(transformed))
    sample_supp_assert(
        identical(
            as.character(archived_md[ids, "cell_line"]),
            annotation$cell_line
        ) &&
            identical(
                as.character(archived_md[ids, "chemical"]),
                annotation$chemical
            ),
        "Frozen/public annotations differ"
    )
    sample_supp_validate_distance(distance, annotation)
    saveRDS(distance, distance_path, compress = "xz", version = 2)
    sample_supp_write_tsv(annotation, annotation_path)
    trial_archive_relative <- "bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/ls_metadata.RDS"
    sample_supp_write_tsv(recovered$selection, selection_path)
    sample_supp_write_tsv(recovered$decisions, decision_path)
    sample_supp_write_tsv(sample_supp_panels, panel_path)
    # Hash the full extracted assay without shipping its ~65 MB matrix.
    assay_tmp <- tempfile(fileext = ".rds")
    saveRDS(assay, assay_tmp, compress = FALSE, version = 2)
    assay_sha <- sample_supp_sha256(assay_tmp)
    unlink(assay_tmp)
    transformation <- data.frame(
        object = "rldList$all",
        features = nrow(assay),
        samples = ncol(assay),
        serialized_assay_sha256 = assay_sha,
        serialization = "R saveRDS version=2 compress=FALSE, original row/column order",
        producer = "scripts/DESeq_functions.R regularize_data: VST blind=FALSE, optional batch removal",
        fidelity = "Exact archived assay distances; not a fresh transform or raw-count proxy",
        deseq2_extraction_version = as.character(utils::packageVersion(
            "DESeq2"
        ))
    )
    sample_supp_write_tsv(
        transformation,
        file.path(ref, "s1_transformation.tsv")
    )
    source_paths <- c(
        archive_relative,
        trial_archive_relative,
        "Rmd/stats_report.Rmd",
        "scripts/DESeq_functions.R",
        "scripts/load-facet-data.R",
        "scripts/run_DESeq2.R",
        "scripts/show_FP_DEGs.R",
        "Rmd/Genes of interest.Rmd",
        "analysis.cell_lines_comparison/DEG_reports/data_explorer_TempO-Seq_FDR_Project_05-05-2023.11.24.html"
    )
    sources <- data.frame(
        role = "legacy_source",
        path = paste0("legacy/", source_paths),
        sha256 = vapply(
            file.path(original, source_paths),
            sample_supp_sha256,
            ""
        )
    )
    sources <- rbind(
        sources,
        read.delim(raster_sources_path)
    )
    public_paths <- c(
        "data/input/count_table.tsv",
        "data/input/metadata.QC_applied.txt",
        "data/reference/s56/size_factors.tsv",
        frozen_path,
        distance_path,
        annotation_path,
        selection_path,
        panel_path,
        evidence_path,
        decision_path,
        raster_sources_path,
        file.path(ref, "s1_transformation.tsv")
    )
    sources <- rbind(
        sources,
        data.frame(
            role = "public_input_or_reference",
            path = public_paths,
            sha256 = vapply(public_paths, sample_supp_sha256, "")
        )
    )
    sample_supp_write_tsv(sources, file.path(ref, "source_provenance.tsv"))
    rm(env, dds, transformed, assay, public_raw)
    gc()
}
provenance <- read.delim(
    file.path(ref, "source_provenance.tsv"),
    check.names = FALSE
)
local <- provenance[provenance$role == "public_input_or_reference", ]
sample_supp_assert(
    all(vapply(local$path, sample_supp_sha256, "") == local$sha256),
    "Public inputs/references changed: review provenance before re-extraction"
)
distance <- readRDS(distance_path)
annotation <- read.delim(annotation_path)
selection <- read.delim(selection_path)
decisions <- read.delim(decision_path)
sample_supp_assert(
    identical(read.delim(panel_path), sample_supp_panels),
    "Compact S3 panel manifest differs from requested panels"
)
sample_supp_assert(
    isTRUE(all.equal(
        recovered$selection,
        selection,
        check.attributes = FALSE
    )) &&
        isTRUE(all.equal(
            recovered$decisions,
            decisions,
            check.attributes = FALSE
        )),
    "S3 selection differs from the deterministic historical raster match"
)
for (trial in unique(selection$trial_name)) {
    selected <- selection[selection$trial_name == trial, ]
    sample_supp_assert(
        identical(selected$sample_id, frozen[[trial]]$sample_id) &&
            identical(selected$group, as.character(frozen[[trial]]$group)),
        "S3 compact selection differs from the actual frozen trial"
    )
}
summary <- sample_supp_build(
    distance,
    annotation,
    selection,
    inputs,
    decisions
)
summary$source_provenance <- provenance
saveRDS(
    summary,
    "data/derived/sample_supplementary_summary.rds",
    compress = "xz",
    version = 2
)
counts_path <- "data/derived/sample_supplementary_s3_counts.tsv"
summary_path <- "data/derived/sample_supplementary_summary.rds"
sample_supp_write_tsv(summary$expression, counts_path)
sample_supp_write_tsv(
    data.frame(
        role = c(rep("rebuild_input", nrow(local)), "derived_summary", "derived_counts"),
        path = c(local$path, summary_path, counts_path),
        sha256 = c(
            local$sha256,
            sample_supp_sha256(summary_path),
            sample_supp_sha256(counts_path)
        )
    ),
    "data/derived/sample_supplementary_provenance.tsv"
)
cat(
    "Rebuilt S1 and S3 from frozen distances and four raster-matched n=4 trials.\n"
)
