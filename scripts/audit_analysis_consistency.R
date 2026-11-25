# Read-only revision audit. Never refits trials or writes publication summaries.
# Run from the repository root; --large also reads the archived ORA/GSEA lists.
args <- commandArgs(TRUE)
root_arg <- grep("^--historical-root=", args, value = TRUE)
historical_root <- if (length(root_arg)) sub("^--historical-root=", "", root_arg) else "../2023_Meier_FDR"
suppressPackageStartupMessages(source("R/config.R"))
suppressPackageStartupMessages(source("R/summaries.R"))
source("R/deg_trials.R")

cat("Analysis consistency audit; read-only\n")
cat("R:", as.character(getRversion()), " historical root:", historical_root, "\n")
deg <- as.data.frame(readRDS(derived_files$deg))
ora <- as.data.frame(readRDS(derived_files$ora))
gsea <- as.data.frame(readRDS(derived_files$gsea))
metadata <- readRDS(derived_files$trial_metadata)
stopifnot(length(metadata) == 32000L, !anyDuplicated(names(metadata)))
read_old <- function(name) as.data.frame(readRDS(file.path(historical_root, name)))
stopifnot(identical(deg, read_old("compiled_summary_boostrapped_DEGs_FDR_data.wrangled.RDS")),
          identical(ora, read_old("compiled_summary_boostrapped_ora_FDR_data.wrangled.RDS")))
cat("DEG and lowercase-ora historical compact objects: identical\n")
cat("Uppercase-ORA historical compact identical:",
    identical(ora, read_old("compiled_summary_boostrapped_ORA_FDR_data.wrangled.RDS")), "\n")
wide <- tidyr::pivot_wider(deg[!is.na(deg$Filter), c("name", "cell_line", "n_per_group", "Filter", "nDEGs")],
                           names_from = Filter, values_from = nDEGs)
compare_counts <- function(a, b, label) {
  cat(label, ": count differences=", sum(a != b), "; positive-status differences=",
      sum((a > 0) != (b > 0)), "\n", sep = "")
}
compare_counts(wide[["R-ODAF + P-Value"]], wide[["P-Value Only"]], "Legacy A versus historical direct-p summary")
compare_counts(wide[["Final DEGs"]], wide[["P-Value + Fold Change"]], "Final versus historical p+FC summary")
print(as.data.frame(wide %>% group_by(cell_line) %>% summarise(
  legacy_A = sum(`R-ODAF + P-Value` > 0), historical_p = sum(`P-Value Only` > 0),
  final = sum(`Final DEGs` > 0), historical_p_fc = sum(`P-Value + Fold Change` > 0))))

reference <- readRDS("data/reference/deg_reference.rds")
for (nm in names(reference)) {
  x <- reference[[nm]]
  res <- x$DESeqRes
  sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
  excluded <- unique(unlist(lapply(x$DESeqjs_filtered_out[2:3], rownames)))
  passes <- !rownames(sig) %in% excluded
  final <- if (nrow(sig)) sig[passes & abs(sig$log2FoldChange) > log2(1.5), , drop = FALSE] else data.frame()
  intermediate <- if (nrow(sig)) res[passes, , drop = FALSE] else data.frame()
  stopifnot(identical(final, x$DESeqResDEGs), identical(intermediate, x$DESeqResFiltered))
  historical <- deg[deg$name == nm & !is.na(deg$Filter), ]
  # Reproduce the original Rmd literally, including NA-index placeholder rows.
  legacy_A <- res[as.data.frame(x$DESeqResFiltered)$padj < 0.05, , drop = FALSE]
  legacy_p <- res[res$padj < 0.05, , drop = FALSE]
  legacy_fc <- res[res$padj < 0.05 & abs(res$log2FoldChange) > log2(1.5), , drop = FALSE]
  expected <- c("R-ODAF + P-Value" = nrow(legacy_A), "P-Value Only" = nrow(legacy_p),
                "P-Value + Fold Change" = nrow(legacy_fc), "Final DEGs" = nrow(final))
  stopifnot(all(historical$nDEGs[match(names(expected), historical$Filter)] == expected))
  cat(nm, ": sig=", nrow(sig), "; legacy_A=", nrow(legacy_A), "; legacy_p=", nrow(legacy_p),
      "; final=", nrow(final), "; p_fc=", sum(!is.na(res$padj) & res$padj < .05 & abs(res$log2FoldChange) > log2(1.5)),
      "; intermediate=", nrow(intermediate), "; NA_padj=", sum(is.na(res$padj)), "; exact=TRUE\n", sep = "")
}
old_gsea <- read_old("compiled_summary_boostrapped_GSEA_FDR_data.wrangled.RDS")
old_gsea <- old_gsea[match(gsea$run, old_gsea$run), ]
stopifnot(!anyNA(old_gsea$run))
compare_counts(gsea$p_lt_0.05, old_gsea$p_lt_0.05, "GSEA current versus historical A")
compare_counts(gsea$p_lt_0.05_ES_gt_2, old_gsea$p_lt_0.05_ES_gt_2, "GSEA current versus historical B")

# Count exactly the table exported by DOSE 3.24.0 get_enriched() followed by
# na.omit(), as in the original Rmd. The archived slots record all cutoffs.
# Direct attribute access avoids requiring the old S4 class/package to be loaded.
archived_ora_count <- function(object) {
  if (is.null(object) || !isS4(object)) return(0L)
  tab <- attr(object, "result")
  p <- attr(object, "pvalueCutoff")
  q <- attr(object, "qvalueCutoff")
  if (length(p)) {
    tab <- tab[tab$pvalue <= p, , drop = FALSE]
    tab <- tab[tab$p.adjust <= p, , drop = FALSE]
  }
  if (length(q) && !any(is.na(tab$qvalue))) tab <- tab[tab$qvalue <= q, , drop = FALSE]
  nrow(na.omit(tab))
}

if ("--large" %in% args) {
  for (cell in cell_line_levels) {
    stem <- c("TK6" = "tk6", "MCF-7" = "mcf7", "HepaRG" = "heparg", "Liver spheroids" = "spheroids")[[cell]]
    expected_names <- names(metadata)[parse_trial_names(names(metadata))$cell_line == cell]
    cat("Reading ORA:", cell, "\n"); flush.console()
    objects <- readRDS(file.path(historical_root, "bootstrap_ORA/2023-11-09", paste0(stem, "_ora.RDS")))
    stopifnot(length(objects) == length(expected_names))
    if (!is.null(names(objects))) stopifnot(identical(names(objects), expected_names))
    counts <- t(vapply(objects, function(x) c(archived_ora_count(x$enrich_go_bp), archived_ora_count(x$enricher_wikipathways)), integer(2)))
    frozen <- ora[match(expected_names, ora$run), ]
    stopifnot(identical(as.integer(counts[, 1]), frozen$go_bp),
              identical(as.integer(counts[, 2]), frozen$enricher_wikipathways))
    for (field in c("enrich_go_bp", "enricher_wikipathways")) {
      s4 <- lapply(objects, `[[`, field)
      s4 <- s4[vapply(s4, isS4, logical(1))]
      settings <- unique(vapply(s4, function(x) paste(attr(x, "pvalueCutoff"), attr(x, "qvalueCutoff"), attr(x, "pAdjustMethod")), character(1)))
      cat("ORA", cell, field, "cutoffs/method:", settings, "; S4 trials=", length(s4), "\n")
      cat("ORA", cell, field, "unique gene/universe slots:",
          all(vapply(s4, function(x) !anyDuplicated(attr(x, "gene")) && !anyDuplicated(attr(x, "universe")), logical(1))), "\n")
    }
    cat("ORA", cell, ": all 8000 trial counts match\n"); flush.console()
    rm(objects, s4); gc()
    cat("Reading GSEA:", cell, "\n"); flush.console()
    objects <- readRDS(file.path(historical_root, "bootstrap_GSEA/2023-11-08", paste0(stem, "_gsea.RDS")))
    objects <- attach_trial_names(objects, metadata, cell)
    fresh <- summarise_gsea_from_trials(objects)
    frozen <- gsea[match(fresh$run, gsea$run), ]
    fields <- c("num_rows_raw", "num_rows", "p_lt_0.05", "p_lt_0.05_ES_gt_2", "rank", "rank_filtered")
    stopifnot(all(vapply(fields, function(field) identical(fresh[[field]], frozen[[field]]), logical(1))))
    raw_counts <- t(vapply(objects, function(x) c(sum(x$padj < .05, na.rm = TRUE),
                                               sum(x$padj < .05 & abs(x$NES) > 2, na.rm = TRUE)), integer(2)))
    historical <- old_gsea[match(fresh$run, old_gsea$run), ]
    stopifnot(all(raw_counts[, 1] == historical$p_lt_0.05), all(raw_counts[, 2] == historical$p_lt_0.05_ES_gt_2))
    cat("GSEA", cell, ": all 8000 current and historical trial counts match their respective pipelines\n"); flush.console()
    rm(objects); gc()
  }
}
cat("Read-only consistency checks passed.\n")
