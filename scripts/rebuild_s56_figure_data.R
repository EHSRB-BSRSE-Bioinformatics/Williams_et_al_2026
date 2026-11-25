script_file <- sub(
    "--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L) {
    setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}
source("R/config.R")
source("R/recurrence_figures.R")

args <- commandArgs(TRUE)
s56_assert(
    length(args) == 0L || identical(args, "--extract-reference"),
    "Usage: Rscript scripts/rebuild_s56_figure_data.R [--extract-reference]"
)
# Keep published provenance paths repository-relative, including reference files.
reference <- file.path("data", "reference", "s56")
dir.create(reference, recursive = TRUE, showWarnings = FALSE)
screen_path <- file.path(reference, "screening_counts.rds")
factor_path <- file.path(reference, "size_factors.tsv")
eligibility_path <- file.path(reference, "screening_eligibility.rds")
input_paths <- c(
    "data/input/count_table.tsv",
    "data/input/metadata.QC_applied.txt"
)
counts <- as.matrix(read.delim(
    input_paths[1L],
    row.names = 1L,
    check.names = FALSE
))
metadata <- read.delim(input_paths[2L], check.names = FALSE)

if (identical(args, "--extract-reference")) {
    cache_dir <- Sys.getenv(
        "FDR_S56_RECURRENCE_DIR",
        unset = "data/large/gene_recurrence"
    )
    original_dir <- Sys.getenv(
        "FDR_ORIGINAL_ANALYSIS",
        unset = "../2023_Meier_FDR"
    )
    workspace <- Sys.getenv(
        "FDR_S56_EXPRESSION_WORKSPACE",
        unset = file.path(
            original_dir,
            "analysis.cell_lines_comparison",
            "RData",
            "FDR Project_DEG_data.RData"
        )
    )
    frozen <- readRDS(derived_files$trial_metadata)
    tables <- provenance <- vector("list", length(s56_stems))
    eligibility_caches <- vector("list", length(s56_stems))
    for (i in seq_along(s56_stems)) {
        model <- names(s56_stems)[i]
        cache_path <- file.path(cache_dir, paste0(s56_stems[i], ".rds"))
        cache <- readRDS(cache_path)
        ids <- names(frozen)[startsWith(names(frozen), paste0(model, "_"))]
        s56_assert(
            identical(cache$cell, model) &&
                length(ids) == 8000L &&
                setequal(names(cache$screen), ids),
            "Cache does not match the frozen 8,000 trials per model."
        )
        tables[[i]] <- do.call(
            rbind,
            lapply(sample_size_levels, function(n) {
                trial_ids <- ids[
                    as.integer(vapply(
                        strsplit(ids, "_", fixed = TRUE),
                        `[`,
                        "",
                        2L
                    )) ==
                        n
                ]
                s56_assert(
                    length(trial_ids) == n_trials_default,
                    "Incomplete screening stratum."
                )
                counts_n <- table(unlist(
                    cache$screen[trial_ids],
                    use.names = FALSE
                ))
                data.frame(
                    cell_line = model,
                    n = n,
                    feature = names(counts_n),
                    count = as.integer(counts_n)
                )
            })
        )
        provenance[[i]] <- data.frame(
            role = c("recurrence_cache", "source_DEG_archive"),
            path = c(
                file.path("recurrence_cache", basename(cache_path)),
                file.path("source_DEG_archive", basename(cache$archive))
            ),
            sha256 = c(s56_sha256(cache_path), cache$sha256)
        )
        eligibility_caches[[i]] <- list(tested = cache$tested, bh = cache$bh)
    }
    screening <- do.call(rbind, tables)
    screening <- screening[
        order(
            match(screening$cell_line, cell_line_levels),
            screening$n,
            screening$feature
        ),
    ]
    rownames(screening) <- NULL
    saveRDS(screening, screen_path, compress = "xz")
    selected_features <- sort(unique(s56_select(screening)$feature))
    eligibility <- do.call(
        rbind,
        lapply(seq_along(s56_stems), function(i) {
            grid <- expand.grid(
                feature = selected_features,
                n = sample_size_levels,
                stringsAsFactors = FALSE
            )
            cache <- eligibility_caches[[i]]
            index <- cbind(
                match(grid$feature, rownames(cache$tested)),
                match(as.character(grid$n), colnames(cache$tested))
            )
            s56_assert(!anyNA(index), "Missing cached testing eligibility.")
            grid$cell_line <- names(s56_stems)[i]
            grid$tested_trials <- cache$tested[index]
            grid$bh_trials <- cache$bh[index]
            grid
        })
    )
    rownames(eligibility) <- NULL
    saveRDS(eligibility, eligibility_path, compress = "xz")

    s56_assert(
        requireNamespace("DESeq2", quietly = TRUE),
        "DESeq2 is needed only for reference extraction; use the deg-study environment."
    )
    env <- new.env()
    load(workspace, envir = env)
    dds <- env$ddsList$all
    s56_assert(
        !is.null(dds) && is.null(DESeq2::normalizationFactors(dds)),
        "Unsupported archived normalization."
    )
    wells <- sort(intersect(
        colnames(dds),
        metadata$Sample_ID[metadata$chemical == "DMSO"]
    ))
    archived_raw <- as.matrix(DESeq2::counts(dds)[, wells, drop = FALSE])
    s56_assert(
        identical(
            unname(archived_raw),
            unname(counts[rownames(dds), wells, drop = FALSE])
        ),
        "Archived S6 counts differ from public inputs."
    )
    factors <- data.frame(
        sample_id = wells,
        size_factor = unname(DESeq2::sizeFactors(dds)[wells])
    )
    # Preserve full precision for exact public reconstruction without DESeq2.
    old_options <- options(digits = 17)
    s56_write_tsv(factors, factor_path)
    options(old_options)
    selected <- unique(s56_select(screening)$feature)
    expected <- log2(
        DESeq2::counts(dds, normalized = TRUE)[selected, wells, drop = FALSE] +
            1
    )
    observed <- log2(
        sweep(
            counts[selected, wells, drop = FALSE],
            2L,
            factors$size_factor,
            "/"
        ) +
            1
    )
    s56_assert(
        isTRUE(all.equal(expected, observed, tolerance = 1e-12)),
        "S6 archived expression reconstruction failed."
    )
    provenance <- do.call(rbind, provenance)
    provenance <- rbind(
        provenance,
        data.frame(
            role = "source_expression_workspace",
            path = basename(workspace),
            sha256 = s56_sha256(workspace)
        )
    )
    original_source <- file.path(
        original_dir,
        "Rmd",
        "Process bootstrapped results.Rmd"
    )
    if (file.exists(original_source)) {
        provenance <- rbind(
            provenance,
            data.frame(
                role = "original_plot_source",
                path = "Rmd/Process bootstrapped results.Rmd",
                sha256 = s56_sha256(original_source)
            )
        )
    }
    for (path in c(
        input_paths,
        "data/derived/monte_carlo_trial_metadata.rds",
        screen_path,
        factor_path,
        eligibility_path
    )) {
        provenance <- rbind(
            provenance,
            data.frame(
                role = "reference_input_or_output",
                path = path,
                sha256 = s56_sha256(path)
            )
        )
    }
    s56_write_tsv(provenance, file.path(reference, "source_provenance.tsv"))
    rm(env, dds, archived_raw, expected, observed)
}

source_provenance <- read.delim(file.path(reference, "source_provenance.tsv"))
for (path in c(input_paths, screen_path, factor_path, eligibility_path)) {
    expected <- source_provenance$sha256[source_provenance$path == path]
    s56_assert(
        length(expected) == 1L && identical(s56_sha256(path), expected),
        paste("Reference checksum mismatch:", path)
    )
}
screening <- readRDS(screen_path)
factors <- read.delim(factor_path)
eligibility <- readRDS(eligibility_path)
summary <- build_s56_summary(screening, factors, counts, metadata, eligibility)
validate_s56_summary(summary)
dir.create("data/derived", showWarnings = FALSE)
saveRDS(summary, s56_summary_path(), compress = "xz")
s56_write_tsv(summary$selected, "data/derived/s56_selected_probes.tsv")
correlations <- s56_baseline_correlations(summary)
s56_write_tsv(correlations$pairs, "data/derived/s56_baseline_correlations.tsv")
s56_write_tsv(
    correlations$summary,
    "data/derived/s56_baseline_correlation_summary.tsv"
)
paths <- c(
    input_paths,
    screen_path,
    factor_path,
    eligibility_path,
    file.path(reference, "source_provenance.tsv"),
    "R/config.R",
    "R/recurrence_figures.R",
    "scripts/rebuild_s56_figure_data.R",
    "data/derived/s56_figure_summary.rds",
    "data/derived/s56_selected_probes.tsv",
    "data/derived/s56_baseline_correlations.tsv",
    "data/derived/s56_baseline_correlation_summary.tsv"
)
s56_write_tsv(
    data.frame(
        role = c(rep("input_or_code", 9L), rep("output", 4L)),
        path = paths,
        sha256 = vapply(paths, s56_sha256, "")
    ),
    "data/derived/s56_data_provenance.tsv"
)
message(
    "S5/S6 summary rebuilt; selected probes and descriptive baseline correlations written."
)
print(correlations$summary, row.names = FALSE)
