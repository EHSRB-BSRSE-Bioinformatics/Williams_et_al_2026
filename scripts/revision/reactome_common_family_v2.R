# Primary v2 common-family sensitivity, using extracted archived raw p-values.
source("R/revision/reactome_ora.R")
root <- "data/large/reactome_ora_v2"; out <- "docs/revision/reactome_ora_v2"
metadata <- trial_identity(names(readRDS("data/derived/monte_carlo_trial_metadata.rds")))
all_stems <- c(TK6="tk6",`MCF-7`="mcf7",HepaRG="heparg",`Liver spheroids`="spheroids")
stems <- all_stems; args <- commandArgs(TRUE)
if(length(args)) stems <- stems[stems %in% args]
source_md5 <- unname(tools::md5sum("scripts/revision/reactome_common_family_v2.R"))
annotation_md5 <- unname(tools::md5sum(file.path(root,"annotation/frozen_annotation.rds")))
for(cell in names(stems)) {
  target <- file.path(root,paste0("common_",stems[[cell]],".rds"))
  if(file.exists(target)) {
    z <- readRDS(target)
    stopifnot(identical(attr(z,"source_md5"),source_md5),identical(attr(z,"annotation_md5"),annotation_md5))
    next
  }
  g <- readRDS(file.path("data/large/reactome_ora_v1",paste0("gsea_family_",stems[[cell]],".rds")))
  ids <- metadata$run[metadata$cell_line==cell]
  stopifnot(identical(names(g),ids))
  rows <- list()
  for(n in 3:10) {
    ora <- readRDS(file.path(root,"full",sprintf("%s_%02d.rds",stems[[cell]],n)))$rows
    for(id in names(ora)) {
      o <- ora[[id]]; terms <- o$result$terms; z <- g[[id]]
      common <- intersect(terms$pathway,z$pathway)
      oi <- match(common,terms$pathway); gi <- match(common,z$pathway)
      op <- p.adjust(terms$pvalue[oi],"BH"); gp <- p.adjust(z$pval[gi],"BH")
      rows[[id]] <- data.frame(run=id,n_common=length(common),ora_family=nrow(terms),gsea_family=nrow(z),
        family_jaccard=length(common)/length(union(terms$pathway,z$pathway)),
        size_disagreements=sum(terms$size[oi]!=z$size[gi]),
        ora_common=if(o$result$status=="insufficient_input") 0L else sum(op<.05),gsea_common=sum(gp<.05))
    }
  }
  result <- do.call(rbind,rows); rownames(result) <- NULL
  stopifnot(setequal(result$run,ids),!anyDuplicated(result$run))
  attr(result,"source_md5") <- source_md5; attr(result,"annotation_md5") <- annotation_md5
  tmp <- paste0(target,".tmp"); saveRDS(result,tmp); stopifnot(file.rename(tmp,target))
  message("Common family complete: ",cell)
  rm(g,ora,rows); gc()
}
if(!all(file.exists(file.path(root,paste0("common_",all_stems,".rds"))))) quit(status=0L)
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
