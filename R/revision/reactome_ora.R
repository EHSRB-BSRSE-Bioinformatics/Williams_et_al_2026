# Historical-compatible ORA inputs; intentionally separate from public helpers.
historical_ora_input <- function(trial, alias_map) {
  x <- as.data.frame(trial$DESeqRes)
  stopifnot(all(c("log2FoldChange", "padj") %in% names(x)), !is.null(rownames(x)))
  alias <- sub("_(.*)", "", rownames(x))
  # which.max ignores NA and uses the first maximum, in original probe order.
  ord <- order(alias, -abs(x$log2FoldChange), seq_len(nrow(x)), na.last = TRUE)
  idx <- ord[!duplicated(alias[ord]) & !is.na(x$log2FoldChange[ord])]
  idx <- idx[order(match(alias[idx], unique(alias)))]
  ids <- alias_map[alias[idx]]
  lens <- lengths(ids)
  mapped <- unlist(ids, use.names = FALSE)
  expanded <- rep(idx, lens)
  fc <- x$log2FoldChange[expanded]
  linear <- ifelse(fc > 0, 2^fc, -1/(2^fc))
  keep <- !is.na(x$padj[expanded]) & x$padj[expanded] < .05 & abs(linear) > 1.5
  list(selected = unique(mapped[keep]), background = unique(mapped),
       selected_expanded = sum(keep), n_features = nrow(x),
       n_aliases = length(idx), n_mapped_aliases = sum(lens > 0L))
}

# Exact upper-tail hypergeometric ORA. The BH family includes EVERY size-eligible
# term, including zero-overlap terms (p=1); q-values are not used or estimated.
# This explicit family is more comparable to GSEA than enricher's hit-only family.
compile_reactome <- function(pathways) {
  genes <- unique(unlist(pathways, use.names = FALSE))
  mat <- Matrix::sparseMatrix(i = rep(seq_along(pathways), lengths(pathways)),
       j = match(unlist(pathways, use.names = FALSE), genes), x = 1,
       dims = c(length(pathways), length(genes)))
  mat@x[] <- 1
  list(genes = genes, matrix = mat, names = names(pathways))
}

reactome_hypergeom <- function(input, pathways, family = names(pathways), compiled = NULL) {
  if (!is.null(compiled)) {
    universe <- intersect(input$background, compiled$genes)
    genes <- intersect(input$selected, universe)
    sizes <- as.numeric(compiled$matrix %*% as.integer(compiled$genes %in% universe))
    overlaps <- as.numeric(compiled$matrix %*% as.integer(compiled$genes %in% genes))
    eligible <- sizes >= 15L & sizes <= 500L & compiled$names %in% family
    sizes <- sizes[eligible]; overlaps <- overlaps[eligible]
    term_names <- compiled$names[eligible]
  } else {
  annotated <- unique(unlist(pathways, use.names = FALSE))
  universe <- intersect(input$background, annotated)
  genes <- intersect(input$selected, universe)
  effective <- lapply(pathways[family], intersect, y = universe)
  sizes <- lengths(effective)
  effective <- effective[sizes >= 15L & sizes <= 500L]
  sizes <- lengths(effective)
  overlaps <- vapply(effective, function(z) sum(z %in% genes), integer(1))
  term_names <- names(effective)
  }
  n <- length(genes); N <- length(universe)
  p <- if (length(sizes) && N > 0L) phyper(overlaps - 1L, sizes, N - sizes, n, lower.tail = FALSE) else numeric()
  terms <- data.frame(pathway = term_names, size = sizes, overlap = overlaps,
                      pvalue = p, padj = p.adjust(p, "BH"), row.names = NULL)
  status <- if (n < 2L || !nrow(terms)) "insufficient_input" else if (any(terms$padj < .05)) "significant" else "none"
  # Prespecified minimum of two unique annotated selected genes.
  if (status == "insufficient_input" && nrow(terms)) terms$padj <- p.adjust(terms$pvalue, "BH")
  list(status = status, n_selected = length(input$selected), n_background = length(input$background),
       n_selected_annotated = n, n_background_annotated = N, n_tested = nrow(terms),
       n_significant = if (status == "insufficient_input") 0L else sum(terms$padj < .05),
       terms = terms)
}

trial_identity <- function(ids) {
  p <- strsplit(ids, "_", fixed = TRUE)
  stopifnot(!anyDuplicated(ids), all(lengths(p) == 4L))
  out <- data.frame(run = ids, cell_line = vapply(p, `[`, "", 1L),
                    n = as.integer(vapply(p, `[`, "", 2L)))
  stopifnot(all(vapply(p, `[`, "", 2L) == vapply(p, `[`, "", 3L)), !anyNA(out$n))
  out
}

check_trial_coverage <- function(results, expected) {
  stopifnot(!anyDuplicated(results$run), setequal(results$run, expected),
            !anyNA(results$status), !any(results$status == "failed"),
            all(table(results$cell_line, results$n) == 1000L))
  invisible(TRUE)
}
