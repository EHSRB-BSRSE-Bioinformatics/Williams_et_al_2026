# Sensitivity analysis only: reuses saved p-values; never calls GSEA.
source("R/revision/reactome_ora.R")
root <- "data/large/reactome_ora_v1"
out <- "docs/revision/reactome_ora_v1"
metadata <- trial_identity(names(readRDS("data/derived/monte_carlo_trial_metadata.rds")))
stems <- c(TK6="tk6", `MCF-7`="mcf7", HepaRG="heparg", `Liver spheroids`="spheroids")
all_stems <- stems
args <- commandArgs(TRUE)
if (length(args)) stems <- stems[stems %in% args]
stopifnot(length(stems)>0L)
source_md5 <- unname(tools::md5sum("scripts/revision/reactome_common_family.R"))
for (cell in names(stems)) {
  target <- file.path(root,paste0("common_",stems[[cell]],".rds"))
  if (file.exists(target)) {
    previous <- readRDS(target)
    stopifnot(identical(attr(previous,"source_md5"),source_md5),setequal(previous$run,metadata$run[metadata$cell_line==cell]))
    next
  }
  message(Sys.time()," reading archived GSEA for common-family comparison: ",cell)
  g <- readRDS(file.path("../2023_Meier_FDR/bootstrap_GSEA/2023-11-08",paste0(stems[[cell]],"_gsea.RDS")))
  ids <- metadata$run[metadata$cell_line==cell]
  stopifnot(length(g)==length(ids))
  if (!is.null(names(g))) stopifnot(identical(names(g),ids))
  names(g) <- ids
  rows <- unlist(lapply(3:10,function(n) readRDS(file.path(root,"full",sprintf("%s_%02d.rds",stems[[cell]],n)))$rows),recursive=FALSE)
  stopifnot(setequal(names(rows),ids))
  result <- lapply(ids,function(id) {
    o <- rows[[id]]; terms <- o$result$terms
    z <- as.data.frame(g[[id]])
    z <- z[!is.na(z$pval)&z$size>=15L&z$size<=500L,]
    common <- intersect(terms$pathway,z$pathway)
    oi <- match(common,terms$pathway); gi <- match(common,z$pathway)
    op <- p.adjust(terms$pvalue[oi],"BH"); gp <- p.adjust(z$pval[gi],"BH")
    data.frame(run=id,n_common=length(common),ora_family=nrow(terms),gsea_family=nrow(z),
      family_jaccard=length(common)/length(union(terms$pathway,z$pathway)),
      size_disagreements=sum(terms$size[oi]!=z$size[gi]),
      ora_common=if(o$result$status=="insufficient_input") 0L else sum(op<.05),
      gsea_common=sum(gp<.05))
  })
  result <- do.call(rbind,result)
  attr(result,"source_md5") <- source_md5
  temp <- paste0(target,".tmp"); saveRDS(result,temp)
  stopifnot(file.rename(temp,target))
  rm(g,rows); gc()
}
if (!all(file.exists(file.path(root,paste0("common_",all_stems,".rds"))))) {
  message("Requested cell(s) complete; remaining common-family cells have not run.")
  quit(status=0L)
}
results <- do.call(rbind,lapply(all_stems,function(s) readRDS(file.path(root,paste0("common_",s,".rds")))))
results <- merge(metadata,results,by="run",sort=FALSE)
stopifnot(nrow(results)==32000L,!anyDuplicated(results$run),all(table(results$cell_line,results$n)==1000L))
saveRDS(results,file.path(out,"common_family_trial_results.rds"))
write.table(results,file.path(out,"common_family_trial_results.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
library(dplyr)
s <- results %>% group_by(cell_line,n) %>% summarise(n_trials=n(),ora_rate=mean(ora_common>0),gsea_rate=mean(gsea_common>0),
 difference=ora_rate-gsea_rate,n_common_median=median(n_common),jaccard_median=median(family_jaccard),
 trials_with_size_disagreement=sum(size_disagreements>0),.groups="drop")
saveRDS(s,file.path(out,"common_family_summary.rds"))
write.table(s,file.path(out,"common_family_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
print(s,n=32)
