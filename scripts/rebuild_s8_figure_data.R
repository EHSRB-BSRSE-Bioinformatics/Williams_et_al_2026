script_file <- sub(
    "--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)
)
if (length(script_file) == 1L && nzchar(script_file)) {
    setwd(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
}

source("R/config.R")
source("R/summaries.R")

summary <- build_s8_figure_summary()
output_path <- derived_files$s8_mean_variance
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(summary, output_path, compress = "xz")

sha256_file <- function(path) {
    if (!file.exists(path)) {
        return(NA_character_)
    }
    result <- system2("sha256sum", shQuote(path), stdout = TRUE)
    if (length(result) != 1L || !grepl("^[[:xdigit:]]{64} ", result)) {
        stop("Could not calculate SHA-256 for ", path)
    }
    sub(" .*", "", result)
}

relative_path <- function(path) {
    root <- normalizePath(project_root(), mustWork = TRUE)
    absolute <- normalizePath(path, mustWork = FALSE)
    prefix <- paste0(root, .Platform$file.sep)
    if (startsWith(absolute, prefix)) {
        substring(absolute, nchar(prefix) + 1L)
    } else {
        basename(path)
    }
}

cache_dir <- Sys.getenv(
    "FDR_S8_RECURRENCE_DIR",
    unset = fdr_path("data", "large", "gene_recurrence")
)
cache_paths <- file.path(
    cache_dir,
    paste0(unname(s8_cache_stems), ".rds")
)
cache_hashes <- vapply(cache_paths, sha256_file, character(1))
if (anyNA(cache_hashes)) {
    stop("A recurrence cache disappeared while building Figure S8 provenance.")
}

archive_paths <- file.path(
    project_root(),
    summary$source_archives$path
)
archive_hashes <- summary$source_archives$sha256
for (i in seq_along(archive_paths)) {
    if (file.exists(archive_paths[[i]])) {
        observed_hash <- sha256_file(archive_paths[[i]])
        if (!identical(observed_hash, archive_hashes[[i]])) {
            stop(
                "Source archive checksum differs from the recurrence cache: ",
                archive_paths[[i]]
            )
        }
    }
}

original_rmd <- fdr_path(
    "..",
    "2023_Meier_FDR",
    "Rmd",
    "Process bootstrapped results.Rmd"
)
source_paths <- c(
    fdr_path("data", "input", "count_table.tsv"),
    fdr_path("data", "input", "metadata.QC_applied.txt"),
    derived_files$trial_metadata,
    fdr_path("R", "config.R"),
    fdr_path("R", "summaries.R"),
    fdr_path("scripts", "rebuild_s8_figure_data.R")
)
source_roles <- c(
    "input",
    "input",
    "frozen_trial_design",
    "code",
    "code",
    "code"
)
provenance <- data.frame(
    role = source_roles,
    path = vapply(source_paths, relative_path, character(1)),
    sha256 = vapply(source_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
)
provenance <- rbind(
    provenance,
    data.frame(
        role = "recurrence_cache",
        path = vapply(cache_paths, relative_path, character(1)),
        sha256 = unname(cache_hashes),
        stringsAsFactors = FALSE
    ),
    data.frame(
        role = "source_DEG_archive",
        path = summary$source_archives$path,
        sha256 = archive_hashes,
        stringsAsFactors = FALSE
    )
)
if (file.exists(original_rmd)) {
    provenance <- rbind(
        provenance,
        data.frame(
            role = "original_plot_source",
            path = file.path("Rmd", basename(original_rmd)),
            sha256 = sha256_file(original_rmd),
            stringsAsFactors = FALSE
        )
    )
}
provenance <- rbind(
    provenance,
    data.frame(
        role = "output",
        path = relative_path(output_path),
        sha256 = sha256_file(output_path),
        stringsAsFactors = FALSE
    )
)
write.table(
    provenance,
    fdr_path("data", "derived", "s8_mean_variance_provenance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

message("Wrote Figure S8 summary to ", output_path)
message(
    "Wrote provenance to ",
    fdr_path("data", "derived", "s8_mean_variance_provenance.tsv")
)
