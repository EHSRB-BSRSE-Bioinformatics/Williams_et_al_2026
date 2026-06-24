# Trial-level ORA and GSEA functions for DESeq2 Monte Carlo outputs.

load_human_alias_annotation <- function() {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("AnnotationDbi and org.Hs.eg.db are required for pathway trials.")
  }

  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  keys <- AnnotationDbi::keys(orgdb, keytype = "ENTREZID")
  anno <- AnnotationDbi::select(
    orgdb,
    keys = keys,
    keytype = "ENTREZID",
    columns = c("SYMBOL", "ENTREZID", "ALIAS", "ENSEMBL")
  )

  anno <- anno[!is.na(anno$ALIAS) & !is.na(anno$ENTREZID), , drop = FALSE]
  unique(anno[, c("ENTREZID", "ALIAS")])
}

collapse_deseq_results_to_entrez <- function(current_trial, annotation) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required for pathway trials.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The dplyr package is required for pathway trials.")
  }

  if (!"DESeqRes" %in% names(current_trial)) {
    stop("DEG trial object is missing DESeqRes.")
  }

  degs <- data.table::data.table(
    as.data.frame(current_trial$DESeqRes),
    keep.rownames = "rn"
  )
  required <- c("rn", "log2FoldChange", "padj")
  missing <- setdiff(required, names(degs))
  if (length(missing) > 0L) {
    stop("DESeqRes is missing required columns: ", paste(missing, collapse = ", "))
  }

  degs <- degs[!is.na(log2FoldChange)]
  degs[, rn := sub("_(.*)", "", rn)]
  degs <- degs[, .SD[which.max(abs(log2FoldChange))], by = rn]

  joined <- dplyr::left_join(
    as.data.frame(degs),
    annotation,
    by = c("rn" = "ALIAS")
  )
  joined <- joined[!is.na(joined$ENTREZID), , drop = FALSE]
  joined <- joined[order(-abs(joined$log2FoldChange)), , drop = FALSE]
  joined <- joined[!duplicated(joined$ENTREZID), , drop = FALSE]
  joined$linearFoldChange <- linear_fold_change(joined$log2FoldChange)
  joined
}

read_wikipathways_gmt <- function(path) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required to read the WikiPathways GMT file.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr and tidyr are required to process WikiPathways annotations.")
  }
  if (!file.exists(path)) {
    stop("WikiPathways GMT file does not exist: ", path)
  }

  wp2gene <- clusterProfiler::read.gmt(path)
  wp2gene <- tidyr::separate(
    wp2gene,
    col = 1,
    into = c("name", "version", "wpid", "org"),
    sep = "%",
    remove = FALSE,
    fill = "right"
  )

  list(
    term2gene = dplyr::distinct(wp2gene[, c("wpid", "gene")]),
    term2name = dplyr::distinct(wp2gene[, c("wpid", "name")])
  )
}

run_ora_trial <- function(current_trial,
                          annotation,
                          wikipathways,
                          alpha = 0.05,
                          linear_fc_filter = 1.5,
                          go_pvalue_cutoff = 0.01,
                          go_qvalue_cutoff = 0.05) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required for ORA trials.")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required for GO ORA trials.")
  }

  collapsed <- collapse_deseq_results_to_entrez(current_trial, annotation)
  entrez_degs <- collapsed$ENTREZID[
    !is.na(collapsed$padj) &
      collapsed$padj < alpha &
      abs(collapsed$linearFoldChange) > linear_fc_filter
  ]
  background <- collapsed$ENTREZID

  results <- list(enrich_go_bp = NA, enricher_wikipathways = NA)
  if (length(entrez_degs) > 1L) {
    results$enrich_go_bp <- clusterProfiler::enrichGO(
      gene = entrez_degs,
      universe = background,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      readable = TRUE,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = go_pvalue_cutoff,
      qvalueCutoff = go_qvalue_cutoff
    )
  }

  matches <- wikipathways$term2gene[wikipathways$term2gene$gene %in% entrez_degs, , drop = FALSE]
  if (nrow(matches) > 1L) {
    results$enricher_wikipathways <- clusterProfiler::enricher(
      entrez_degs,
      universe = background,
      pAdjustMethod = "fdr",
      pvalueCutoff = alpha,
      TERM2GENE = wikipathways$term2gene,
      TERM2NAME = wikipathways$term2name
    )
  }

  results
}

load_reactome_pathways <- function() {
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    stop("reactome.db is required for Reactome GSEA trials.")
  }
  pathways <- as.list(reactome.db::reactomePATHID2EXTID)
  pathways[grep("R-HSA-", names(pathways))]
}

run_gsea_trial <- function(current_trial,
                           annotation,
                           pathways,
                           min_size = 1L,
                           max_size = gsea_max_size_default,
                           n_perm_simple = 500L) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea is required for GSEA trials.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required for GSEA trials.")
  }

  collapsed <- collapse_deseq_results_to_entrez(current_trial, annotation)
  ranked <- collapsed[!is.na(collapsed$log2FoldChange), , drop = FALSE]
  ranked <- ranked[order(ranked$log2FoldChange), , drop = FALSE]
  ranked_genes <- ranked$log2FoldChange
  names(ranked_genes) <- ranked$ENTREZID
  ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranked_genes,
    minSize = min_size,
    maxSize = max_size,
    nPermSimple = n_perm_simple
  )
}

run_pathway_trials <- function(deg_trials,
                               analysis = c("ora", "gsea"),
                               annotation,
                               wikipathways = NULL,
                               reactome_pathways = NULL,
                               cores = 1L,
                               ...) {
  analysis <- match.arg(analysis)
  extra_args <- list(...)

  fun <- switch(
    analysis,
    ora = function(current_trial) {
      do.call(
        run_ora_trial,
        c(
          list(
            current_trial = current_trial,
            annotation = annotation,
            wikipathways = wikipathways
          ),
          extra_args
        )
      )
    },
    gsea = function(current_trial) {
      do.call(
        run_gsea_trial,
        c(
          list(
            current_trial = current_trial,
            annotation = annotation,
            pathways = reactome_pathways
          ),
          extra_args
        )
      )
    }
  )

  out <- run_trials_parallel(deg_trials, fun, cores = cores)
  names(out) <- names(deg_trials)
  out
}
