# Targeted ranking/collection documentation only; no fgsea calls.
source("R/revision/reactome_ora.R")
a <- readRDS("data/large/reactome_ora_v1/annotation/frozen_annotation_org3.17.rds")
ref <- readRDS("data/reference/deg_reference.rds")
out <- lapply(names(ref),function(id) {
  r <- as.data.frame(ref[[id]]$DESeqRes)
  alias <- sub("_(.*)","",rownames(r))
  # Original group order and first maximum; expand every annotation match.
  idx <- unlist(lapply(split(seq_len(nrow(r)),factor(alias,levels=unique(alias))),function(i) i[which.max(abs(r$log2FoldChange[i]))]),use.names=FALSE)
  maps <- a$alias_map[alias[idx]]
  stat <- rep(r$log2FoldChange[idx],lengths(maps))
  names(stat) <- unlist(maps,use.names=FALSE)
  stat <- stat[order(stat)]
  data.frame(run=id,n_ranked=length(stat),n_unique_entrez=length(unique(names(stat))),
             duplicated_names=sum(duplicated(names(stat))),
             tied_nonzero=sum(duplicated(stat[stat!=0])),zero_stats=sum(stat==0),
             missing_stats=sum(is.na(stat)))
})
write.table(do.call(rbind,out),"docs/revision/reactome_ora_v1/gsea_reference_rankings.tsv",sep="\t",row.names=FALSE,quote=FALSE)
print(do.call(rbind,out))
