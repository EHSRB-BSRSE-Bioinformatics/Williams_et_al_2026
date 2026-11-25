# Primary v2: the archive-compatible Reactome release, using saved original
# DESeqRes-derived inputs from v1. Never fits or rereads DESeq2 results.
source("R/revision/reactome_ora.R")
args <- commandArgs(TRUE)
root <- "data/large/reactome_ora_v2"
dir.create(file.path(root,"full"),recursive=TRUE,showWarnings=FALSE)
annotation_file <- file.path(root,"annotation/frozen_annotation.rds")
a <- readRDS(annotation_file)
compiled <- compile_reactome(a$pathways)
annotation_md5 <- unname(tools::md5sum(annotation_file))
source_md5 <- unname(tools::md5sum(c("R/revision/reactome_ora.R","scripts/revision/retest_reactome_collection.R")))
cores <- as.integer(Sys.getenv("REACTOME_ORA_CORES","1"))
stopifnot(!is.na(cores),cores>=1L)
metadata <- trial_identity(names(readRDS("data/derived/monte_carlo_trial_metadata.rds")))
if("--pilot" %in% args) {
  refs <- readRDS("data/reference/deg_reference.rds")
  rows <- lapply(names(refs),function(id) {
    input <- historical_ora_input(refs[[id]],a$alias_map)
    direct <- reactome_hypergeom(input,a$pathways)
    fast <- reactome_hypergeom(input,a$pathways,compiled=compiled)
    stopifnot(isTRUE(all.equal(direct,fast)),identical(fast,reactome_hypergeom(input,a$pathways,compiled=compiled)))
    list(run=id,input=input,result=fast)
  })
  names(rows) <- names(refs)
  saveRDS(rows,file.path(root,"pilot.rds"))
  cat("Reactome 1.84.0: 17 pilot direct/compiled/repeat checks passed.\n")
  quit(status=0L)
}
stems <- c(TK6="tk6",`MCF-7`="mcf7",HepaRG="heparg",`Liver spheroids`="spheroids")
if(length(args)) stems <- stems[stems %in% args]
stopifnot(length(stems)>0L)
for(cell in names(stems)) for(n in 3:10) {
  filename <- sprintf("%s_%02d.rds",stems[[cell]],n)
  target <- file.path(root,"full",filename)
  ids <- metadata$run[metadata$cell_line==cell&metadata$n==n]
  if(file.exists(target)) {
    old <- readRDS(target)
    stopifnot(identical(names(old$rows),ids),identical(old$source_md5,source_md5),identical(old$annotation_md5,annotation_md5),
      !any(vapply(old$rows,function(z)z$result$status=="failed",logical(1))))
    next
  }
  input_file <- file.path("data/large/reactome_ora_v1/full",filename)
  if(!file.exists(input_file)) stop("Prepared-input checkpoint not yet available: ",input_file)
  inputs <- readRDS(input_file)$rows
  stopifnot(identical(names(inputs),ids),!any(vapply(inputs,function(z)z$result$status=="failed",logical(1))))
  start <- proc.time()
  rows <- parallel::mclapply(inputs,function(z) {
    z$result <- tryCatch(reactome_hypergeom(z$input,a$pathways,compiled=compiled),
      error=function(e)list(status="failed",error=conditionMessage(e)))
    z
  },mc.cores=cores,mc.preschedule=TRUE,mc.set.seed=FALSE)
  stopifnot(all(vapply(rows,is.list,logical(1))))
  payload <- list(schema="reactome_ora_v2",rows=rows,elapsed=(proc.time()-start)[["elapsed"]],
    session=sessionInfo(),annotation_md5=annotation_md5,source_md5=source_md5,
    prepared_input_md5=unname(tools::md5sum(input_file)),cores=cores)
  tmp <- paste0(target,".tmp"); con <- gzfile(tmp,"wb",compression=1L)
  saveRDS(payload,con,version=3); close(con)
  stopifnot(file.rename(tmp,target))
  message(Sys.time()," checkpoint ",target," elapsed=",payload$elapsed)
  rm(inputs,rows,payload); gc()
}
