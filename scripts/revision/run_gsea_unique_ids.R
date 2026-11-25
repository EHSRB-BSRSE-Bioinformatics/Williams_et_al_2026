source("R/revision/gsea_unique_ids.R")
source("R/revision/reactome_ora.R")
gsea_environment()
args <- commandArgs(TRUE)
mode <- if("--prepare" %in% args) "prepare" else if("--full" %in% args) "full" else "pilot"
if(mode=="prepare") suppressPackageStartupMessages(library(DESeq2))
cores <- as.integer(Sys.getenv("GSEA_UNIQUE_CORES","16"))
stopifnot(is.finite(cores),cores>=1L)
a <- gsea_annotation()
ids <- names(readRDS("data/derived/monte_carlo_trial_metadata.rds"))
identity <- trial_identity(ids)
seeds <- setNames(20261116L+seq_along(ids),ids)
sources <- c("R/revision/gsea_unique_ids.R","scripts/revision/run_gsea_unique_ids.R")
provenance <- list(code_md5=tools::md5sum(sources), annotation_md5=tools::md5sum(c(
  "docs/revision/reactome_ora_v2/frozen_annotation.rds",
  "data/large/reactome_ora_v1/workspace_provenance/pathways.rds")),
  metadata_md5=tools::md5sum("data/derived/monte_carlo_trial_metadata.rds"),
  fgsea_source_md5=tools::md5sum("data/large/reactome_ora_v1/annotation/fgsea_1.26.0.tar.gz"),
  environment_md5=tools::md5sum(file.path(gsea_docs,"conda-explicit.txt")))
read_checked <- function(path,expected=NULL) {
  z <- readRDS(path)
  stopifnot(identical(z$provenance,provenance))
  if(!is.null(expected)) stopifnot(identical(names(z$rows),expected))
  z
}
parallel_map <- function(x,f) {
  z <- parallel::mclapply(x,f,mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE)
  stopifnot(all(vapply(z,is.list,logical(1))))
  names(z) <- x
  z
}
refs <- readRDS("data/reference/deg_reference.rds")
ora <- readRDS("docs/revision/reactome_ora_v2/trial_results.rds")
# Fixed before enrichment: reference set plus min/max feature count at both
# endpoint sample sizes for every cell model. Ties use frozen metadata order.
pilot_ids <- names(refs)
for(cell in names(gsea_stems)) for(n in c(3L,10L)) {
  z <- ora[match(identity$run[identity$cell_line==cell & identity$n==n],ora$run),]
  pilot_ids <- c(pilot_ids,z$run[c(which.min(z$n_features),which.max(z$n_features))])
}
pilot_ids <- ids[ids %in% pilot_ids]
write.table(data.frame(run=pilot_ids,seed=seeds[pilot_ids]),file.path(gsea_docs,"pilot_selection.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(gsea_docs,"analysis_sessionInfo.txt"))
if(mode=="full") {
  gate <- read_checked(file.path(gsea_root,"pilot.rds"),pilot_ids)
  stopifnot(isTRUE(gate$validated))
}
for(cell in names(gsea_stems)) {
  stem <- gsea_stems[[cell]]
  cell_ids <- identity$run[identity$cell_line==cell]
  if(mode=="prepare") {
    paths <- file.path(gsea_root,"rankings",sprintf("%s_%02d.rds",stem,3:10))
    if(all(file.exists(paths))) {
      for(n in 3:10) read_checked(paths[n-2L],cell_ids[identity$n[match(cell_ids,ids)]==n])
      next
    }
    message(Sys.time()," reading frozen DEG archive: ",cell)
    archive <- file.path("../2023_Meier_FDR/bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10",paste0(stem,".RDS"))
    trials <- readRDS(archive)
    message(Sys.time()," archive loaded: ",cell)
    archive_md5 <- tools::md5sum(archive)
    stopifnot(length(trials)==length(cell_ids))
    if(!is.null(names(trials))) stopifnot(identical(names(trials),cell_ids))
    names(trials) <- cell_ids
    for(n in 3:10) {
      current <- cell_ids[identity$n[match(cell_ids,ids)]==n]
      path <- paths[n-2L]
      if(file.exists(path)) {read_checked(path,current); next}
      rows <- parallel_map(current,function(id) gsea_ranking(trials[[id]],a$alias_map))
      stopifnot(all(vapply(rows,function(r)!anyDuplicated(names(r$unique)),logical(1))))
      gsea_atomic_save(list(rows=rows,provenance=provenance,archive_md5=archive_md5),path)
      message(Sys.time()," ranking checkpoint ",path)
    }
    rm(trials); gc()
  } else if(mode=="full") {
    for(n in 3:10) {
      current <- cell_ids[identity$n[match(cell_ids,ids)]==n]
      ranking_path <- file.path(gsea_root,"rankings",sprintf("%s_%02d.rds",stem,n))
      path <- file.path(gsea_root,"full",sprintf("%s_%02d.rds",stem,n))
      if(file.exists(path)) {
        z <- read_checked(path,current)
        stopifnot(identical(z$ranking_md5,tools::md5sum(ranking_path)),
          !any(vapply(z$rows,function(r) r$result$status=="failed",logical(1))))
        message("Verified, skipping ",path); next
      }
      ranks <- read_checked(ranking_path,current)$rows
      start <- proc.time()[["elapsed"]]
      rows <- parallel_map(current,function(id) {
        tryCatch(list(run=id,diagnostics=ranks[[id]]$diagnostics,
          result=gsea_evaluate(ranks[[id]]$unique,a$pathways,seeds[[id]])),
          error=function(e)list(run=id,result=list(status="failed",error=conditionMessage(e))))
      })
      gsea_atomic_save(list(rows=rows,provenance=provenance,ranking_md5=tools::md5sum(ranking_path),
        cores=cores,elapsed=proc.time()[["elapsed"]]-start,session=sessionInfo()),path)
      stopifnot(!any(vapply(rows,function(r)r$result$status=="failed",logical(1))))
      message(Sys.time()," GSEA checkpoint ",path)
      rm(ranks,rows); gc()
    }
  }
}
if(mode=="pilot") {
  # --references-only is an early timing diagnostic; cannot authorize full run.
  references_only <- "--references-only" %in% args
  selected <- if(references_only) intersect(pilot_ids,names(refs)) else pilot_ids
  ranks <- list()
  for(id in selected) {
    if(id %in% names(refs)) ranks[[id]] <- gsea_ranking(refs[[id]],a$alias_map) else {
      meta <- identity[match(id,ids),]
      p <- file.path(gsea_root,"rankings",sprintf("%s_%02d.rds",gsea_stems[[meta$cell_line]],meta$n))
      ranks[[id]] <- read_checked(p)$rows[[id]]
    }
  }
  rows <- parallel_map(selected,function(id) {
    r <- ranks[[id]]
    u <- gsea_evaluate(r$unique,a$pathways,seeds[[id]])
    d <- gsea_evaluate(r$historical,a$pathways,seeds[[id]])
    repeat_u <- gsea_evaluate(r$unique,a$pathways,seeds[[id]])
    repeat_d <- gsea_evaluate(r$historical,a$pathways,seeds[[id]])
    stopifnot(identical(u$raw,repeat_u$raw),identical(d$raw,repeat_d$raw))
    list(run=id,diagnostics=r$diagnostics,unique=u,duplicated=d)
  })
  # Also verify outside fork workers to catch scheduling/environment differences.
  for(id in selected[seq_len(min(4L,length(selected)))]) {
    z <- gsea_evaluate(ranks[[id]]$unique,a$pathways,seeds[[id]])
    stopifnot(identical(z$raw,rows[[id]]$unique$raw))
  }
  gsea_atomic_save(list(rows=rows,provenance=provenance,validated=!references_only,
    cores=cores,session=sessionInfo()),file.path(gsea_root,if(references_only) "reference_pilot.rds" else "pilot.rds"))
  stats <- do.call(rbind,lapply(rows,function(z)data.frame(run=z$run,
    unique_significant=z$unique$n_significant,duplicate_significant=z$duplicated$n_significant,
    unique_undefined=z$unique$n_undefined,duplicate_undefined=z$duplicated$n_undefined,
    unique_elapsed=z$unique$elapsed,duplicate_elapsed=z$duplicated$elapsed)))
  write.table(stats,file.path(gsea_docs,if(references_only) "reference_pilot.tsv" else "pilot.tsv"),
    sep="\t",quote=FALSE,row.names=FALSE)
  print(stats)
  cat("PASS: repeated unique and duplicate runs identical; serial/fork checks identical\n")
}
