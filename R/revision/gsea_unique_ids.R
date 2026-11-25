# Controlled deduplication of historical GSEA inputs. No DESeq2 fitting.
gsea_root <- "data/large/gsea_unique_ids_v1"
gsea_docs <- "docs/revision/gsea_unique_ids_v1"
gsea_stems <- c(TK6="tk6", `MCF-7`="mcf7", HepaRG="heparg", `Liver spheroids`="spheroids")

gsea_environment <- function() {
  .libPaths(c(file.path(gsea_root,"library"), .libPaths()))
  stopifnot(getRversion() == "4.2.3", packageVersion("fgsea") == "1.26.0")
  # Load namespaces once in the parent, before forking per-trial workers.
  loadNamespace("fgsea")
  data.table::setDTthreads(1L)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  invisible(TRUE)
}

gsea_annotation <- function() {
  a <- readRDS("docs/revision/reactome_ora_v2/frozen_annotation.rds")
  # Restore the exact saved pathway order as well as memberships.
  p <- readRDS("data/large/reactome_ora_v1/workspace_provenance/pathways.rds")
  p <- p[grep("R-HSA-",names(p))]
  stopifnot(setequal(names(p),names(a$pathways)),
    all(vapply(names(p),function(id)setequal(p[[id]],a$pathways[[id]]),logical(1))))
  a$pathways <- p
  a
}

gsea_ranking <- function(trial, alias_map) {
  x <- as.data.frame(trial$DESeqRes)
  alias <- sub("_(.*)", "", rownames(x))
  stopifnot("log2FoldChange" %in% names(x), !is.null(rownames(x)))
  ord <- order(alias, -abs(x$log2FoldChange), seq_len(nrow(x)), na.last=TRUE)
  idx <- ord[!duplicated(alias[ord]) & !is.na(x$log2FoldChange[ord])]
  idx <- idx[order(match(alias[idx],unique(alias)))]
  maps <- alias_map[alias[idx]]
  score <- rep(x$log2FoldChange[idx], lengths(maps))
  id <- unlist(maps,use.names=FALSE)
  stopifnot(length(id)==length(score), !anyNA(id), all(nzchar(id)), all(is.finite(score)))
  historical <- setNames(score,id)
  historical <- historical[order(score,seq_along(score))]
  # First eliminate identical pairs, then choose max |score| per Entrez.
  # Exact opposite-sign or same-sign ties retain first historical mapping row.
  pairs <- !duplicated(data.frame(id=id,score=score))
  ii <- which(pairs)
  distinct <- table(id[ii])
  choose <- ii[order(-abs(score[ii]),ii)]
  choose <- choose[!duplicated(id[choose])]
  choose <- choose[order(score[choose],choose)]
  unique_stats <- setNames(score[choose],id[choose])
  stopifnot(!anyDuplicated(names(unique_stats)),setequal(names(unique_stats),id))
  list(historical=historical, unique=unique_stats,
       diagnostics=list(n_features=nrow(x), n_aliases=length(idx),
         n_historical=length(score), n_unique=length(choose),
         repeated_pairs=sum(!pairs), multi_score_ids=sum(distinct>1L),
         discarded_distinct_scores=sum(pairs)-length(choose),
         tied_nonzero_historical=sum(duplicated(score[score!=0])),
         tied_nonzero_unique=sum(duplicated(unique_stats[unique_stats!=0])),
         zero_unique=sum(unique_stats==0)))
}

gsea_postprocess <- function(raw) {
  stopifnot(!anyDuplicated(raw$pathway),all(is.na(raw$pval) | (is.finite(raw$pval) & raw$pval>=0 & raw$pval<=1)))
  eligible <- raw$size>=15L & raw$size<=500L
  family <- eligible & !is.na(raw$pval)
  raw$padj_size15 <- NA_real_
  raw$padj_size15[family] <- p.adjust(raw$pval[family],"BH")
  stopifnot(all(is.finite(raw$NES[family])))
  sig <- raw$padj_size15[family] < .05
  list(raw=raw, n_eligible=sum(eligible), n_undefined=sum(eligible & is.na(raw$pval)),
       n_undefined_all=sum(is.na(raw$pval)), n_tested=sum(family),
       n_significant=sum(sig), n_significant_NES2=sum(sig & abs(raw$NES[family])>2),
       status=if(!any(family)) "untestable" else if(any(sig)) "significant" else "none")
}

gsea_evaluate <- function(ranking, pathways, seed) {
  stopifnot(all(is.finite(ranking)),!anyNA(names(ranking)))
  set.seed(seed)
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  raw <- withCallingHandlers(fgsea::fgseaMultilevel(pathways,ranking,
    minSize=1, maxSize=500, nPermSimple=500, sampleSize=101, eps=1e-50,
    scoreType="std", gseaParam=1, nproc=0,
    BPPARAM=BiocParallel::SerialParam()),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
    })
  result <- gsea_postprocess(as.data.frame(raw))
  result$warnings <- unique(warnings)
  result$seed <- seed
  result$elapsed <- proc.time()[["elapsed"]]-start
  result
}

gsea_atomic_save <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temp <- paste0(path,".tmp-",Sys.getpid())
  saveRDS(x,temp,compress="gzip")
  stopifnot(file.rename(temp,path))
}
