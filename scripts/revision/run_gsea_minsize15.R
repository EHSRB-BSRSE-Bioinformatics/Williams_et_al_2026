# Sensitivity analysis using the existing unique-Entrez ranking checkpoints.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("--full", "--pilot")) {
  stop("Usage: Rscript scripts/revision/run_gsea_minsize15.R --full|--pilot")
}
pilot_only <- args == "--pilot"

source("R/revision/gsea_unique_ids.R")
source("R/revision/reactome_ora.R")
gsea_environment()

cores <- suppressWarnings(as.integer(Sys.getenv("GSEA_UNIQUE_CORES", "16")))
if (is.na(cores) || cores < 1L) stop("GSEA_UNIQUE_CORES must be a positive integer.")

pilot <- readRDS(file.path(gsea_root, "pilot.rds"))
stopifnot(isTRUE(pilot$validated), length(pilot$rows) == 33L)
frozen_inputs <- do.call(c, unname(pilot$provenance))
if (!identical(tools::md5sum(names(frozen_inputs)), frozen_inputs)) {
  stop("The frozen ranking inputs differ from the validated unique-ID run.")
}

annotation <- gsea_annotation()
pathways <- names(annotation$pathways)
ids <- names(readRDS("data/derived/monte_carlo_trial_metadata.rds"))
trial_info <- trial_identity(ids)
stopifnot(length(ids) == 32000L, setequal(trial_info$cell_line, names(gsea_stems)),
          setequal(trial_info$n, 3:10), all(table(trial_info$cell_line, trial_info$n) == 1000L),
          all(names(pilot$rows) %in% ids))
seeds <- setNames(20261116L + seq_along(ids), ids)
archives <- read.delim(file.path(gsea_docs, "source_archives.tsv"), colClasses = "character")
posthoc <- readRDS(file.path(gsea_docs, "trial_results.rds"))
stopifnot(!anyDuplicated(posthoc$run), setequal(posthoc$run, ids))
posthoc <- posthoc[match(ids, posthoc$run), ]
stopifnot(identical(posthoc$seed, unname(seeds)),
          !anyNA(posthoc[c("n_significant", "n_significant_NES2")]))

out <- file.path("data/large/gsea_unique_ids_minsize15_v1",
                 if (pilot_only) "pilot" else "full")
target_ids <- if (pilot_only) ids[ids %in% names(pilot$rows)] else ids
selected_info <- trial_info[match(target_ids, ids), ]
blocks <- split(target_ids, sprintf("%s_%02d", gsea_stems[selected_info$cell_line], selected_info$n))
rank_files <- file.path(gsea_root, "rankings", paste0(names(blocks), ".rds"))
if (any(!file.exists(rank_files))) {
  stop("Missing frozen ranking checkpoints: ", paste(rank_files[!file.exists(rank_files)], collapse = ", "))
}
provenance <- list(
  min_size = 15L, max_size = 500L, frozen_inputs = frozen_inputs,
  runner_md5 = tools::md5sum("scripts/revision/run_gsea_minsize15.R"),
  versions = vapply(c("fgsea", "BiocParallel", "data.table"), function(package) {
    as.character(packageVersion(package))
  }, character(1))
)

evaluate <- function(ranking, id) {
  stopifnot(length(ranking) > 0L, length(names(ranking)) == length(ranking),
            !anyNA(names(ranking)), all(nzchar(names(ranking))),
            !anyDuplicated(names(ranking)), all(is.finite(ranking)))
  set.seed(seeds[[id]])
  warnings <- character()
  raw <- withCallingHandlers(
    fgsea::fgseaMultilevel(annotation$pathways, ranking, minSize = 15, maxSize = 500,
                           nPermSimple = 500, sampleSize = 101, eps = 1e-50,
                           scoreType = "std", gseaParam = 1, nproc = 0,
                           BPPARAM = BiocParallel::SerialParam()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  result <- gsea_postprocess(as.data.frame(raw))
  stopifnot(result$n_eligible == nrow(raw))
  result$seed <- seeds[[id]]
  result$warnings <- unique(warnings)
  result
}

summarise_block <- function(results, model, n) {
  trial_counts <- do.call(rbind, lapply(names(results), function(id) {
    z <- results[[id]]
    stopifnot(z$seed == seeds[[id]], z$status %in% c("significant", "none", "untestable"),
              z$n_eligible == nrow(z$raw), all(z$raw$size >= 15L & z$raw$size <= 500L))
    data.frame(run = id, cell_line = model, n = n, seed = z$seed, status = z$status,
               n_eligible = z$n_eligible, n_tested = z$n_tested, n_undefined = z$n_undefined,
               direct_significant = z$n_significant, direct_significant_NES2 = z$n_significant_NES2)
  }))
  counts <- matrix(0L, nrow = length(pathways), ncol = 4L)
  for (z in results) {
    idx <- match(z$raw$pathway, pathways)
    tested <- !is.na(z$raw$padj_size15)
    significant <- tested & z$raw$padj_size15 < .05
    stopifnot(!anyNA(idx), sum(tested) == z$n_tested, sum(significant) == z$n_significant,
              sum(significant & abs(z$raw$NES) > 2) == z$n_significant_NES2)
    counts[idx, ] <- counts[idx, ] + cbind(TRUE, tested, significant, significant & abs(z$raw$NES) > 2)
  }
  pathway_counts <- data.frame(cell_line = model, n = n, pathway = pathways,
                               n_trials = nrow(trial_counts), n_eligible = counts[, 1L],
                               n_tested = counts[, 2L], n_significant = counts[, 3L],
                               n_significant_NES2 = counts[, 4L])
  list(trials = trial_counts, pathways = pathway_counts[pathway_counts$n_eligible > 0L, ])
}

run_block <- function(block, selected, rank_file) {
  checkpoint <- file.path(out, paste0(block, ".rds"))
  ranking_md5 <- tools::md5sum(rank_file)
  first <- match(selected[[1L]], ids)
  model <- trial_info$cell_line[first]
  n <- trial_info$n[first]
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    if (!identical(saved$provenance, provenance) || !identical(saved$ranking_md5, ranking_md5) ||
        !identical(names(saved$rows), selected)) {
      stop("Checkpoint inputs differ from the current run: ", checkpoint)
    }
    message("Verified, reusing ", checkpoint)
    return(summarise_block(saved$rows, model, n))
  }

  rankings <- readRDS(rank_file)
  expected <- trial_info$run[trial_info$cell_line == model & trial_info$n == n]
  source_md5 <- archives$md5[match(paste0(gsea_stems[[model]], ".RDS"), archives$archive)]
  stopifnot(identical(rankings$provenance, pilot$provenance), identical(names(rankings$rows), expected),
            identical(unname(rankings$archive_md5), source_md5))
  message("Running ", block, " (", length(selected), " trials)")
  results <- parallel::mclapply(selected, function(id) {
    tryCatch(evaluate(rankings$rows[[id]]$unique, id), error = identity)
  }, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = TRUE)
  failed <- which(!vapply(results, function(z) is.list(z) && !inherits(z, "error"), logical(1)))
  if (length(failed)) {
    i <- failed[[1L]]
    stop("fgsea failed in ", selected[[i]], ": ", paste(as.character(results[[i]]), collapse = " "))
  }
  names(results) <- selected
  summaries <- summarise_block(results, model, n)
  gsea_atomic_save(list(provenance = provenance, ranking_md5 = ranking_md5,
                        rows = results, session = sessionInfo()), checkpoint)
  message("Saved ", checkpoint)
  summaries
}

summaries <- Map(run_block, names(blocks), blocks, rank_files)
trials <- do.call(rbind, lapply(summaries, `[[`, "trials"))
stopifnot(!anyDuplicated(trials$run), setequal(trials$run, target_ids))
trials <- trials[match(target_ids, trials$run), ]
baseline <- match(trials$run, posthoc$run)
trials$posthoc_significant <- posthoc$n_significant[baseline]
trials$posthoc_significant_NES2 <- posthoc$n_significant_NES2[baseline]

summarise_rates <- function(z, model, n) {
  direct <- z$direct_significant > 0L
  prior <- z$posthoc_significant > 0L
  data.frame(cell_line = model, n = n, n_trials = nrow(z),
             direct_positive = sum(direct), posthoc_positive = sum(prior),
             direct_rate = mean(direct), posthoc_rate = mean(prior),
             direct_NES2_positive = sum(z$direct_significant_NES2 > 0L),
             posthoc_NES2_positive = sum(z$posthoc_significant_NES2 > 0L),
             lost_positive = sum(prior & !direct), gained_positive = sum(!prior & direct),
             undefined_pathways = sum(z$n_undefined))
}
groups <- split(trials, list(trials$cell_line, trials$n), drop = TRUE)
rates <- do.call(rbind, c(lapply(groups, function(z) summarise_rates(z, z$cell_line[1L], z$n[1L])),
                         list(summarise_rates(trials, "All", NA_integer_))))
write_tsv <- function(data, filename) {
  path <- file.path(out, filename)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  write.table(data, temporary, sep = "\t", quote = FALSE, row.names = FALSE)
  stopifnot(file.rename(temporary, path))
}
write_tsv(trials, "trial_results.tsv")
write_tsv(rates, "summary.tsv")
write_tsv(do.call(rbind, lapply(summaries, `[[`, "pathways")), "pathway_counts.tsv")
print(tail(rates, 1L), row.names = FALSE)
message("Validated ", nrow(trials), " paired GSEA trials. Results: ", out)
