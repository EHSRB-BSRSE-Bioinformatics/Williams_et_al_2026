source("R/revision/gsea_unique_ids.R")
gsea_environment()
a <- gsea_annotation()
refs <- readRDS("data/reference/deg_reference.rds")
out <- lapply(names(refs),function(id) {
  x <- as.data.frame(refs[[id]]$DESeqRes)
  # Independent literal historical data.table + dplyr reconstruction.
  d <- data.table::data.table(x,keep.rownames=TRUE)[,.(log2FoldChange,rn)]
  d[,rn:=sub("_(.*)","",rn)]
  d <- d[,.SD[which.max(abs(log2FoldChange))],by=rn]
  d <- d[d$rn %in% a$short_annotation$ALIAS,]
  joined <- merge(as.data.frame(d),transform(a$short_annotation,.map_order=seq_len(nrow(a$short_annotation))),
    by.x="rn",by.y="ALIAS",sort=FALSE)
  joined <- joined[order(match(joined$rn,d$rn),joined$.map_order),]
  expected <- setNames(joined$log2FoldChange,joined$ENTREZID)
  expected <- expected[order(expected)]
  r <- gsea_ranking(refs[[id]],a$alias_map)
  stopifnot(identical(r$historical,expected), identical(r,gsea_ranking(refs[[id]],a$alias_map)))
  maxima <- tapply(abs(expected),names(expected),max)
  stopifnot(all(abs(r$unique)==maxima[names(r$unique)]))
  sizes <- lengths(lapply(a$pathways,intersect,y=names(r$unique)))
  old_sizes <- lengths(lapply(a$pathways,intersect,y=names(expected)))
  stopifnot(identical(sizes,old_sizes))
  data.frame(run=id,as.data.frame(r$diagnostics),size_changes=sum(sizes!=old_sizes),
    eligible_pathways=sum(sizes>=15 & sizes<=500))
})
# Probe-first selection, opposite-sign Entrez ties, repeated annotation pairs,
# unmapped aliases, all-NA alias groups and deterministic equal-score order.
toy <- list(DESeqRes=data.frame(log2FoldChange=c(1,-2,2,NA,3),
  row.names=c("A_1","A_2","B","C","UNMAPPED")))
z <- gsea_ranking(toy,list(A=c("2","2","3"),B=c("2","1"),C="4"))
stopifnot(identical(z$unique,c(`2`=-2,`3`=-2,`1`=2)),z$diagnostics$multi_score_ids==1L,
          z$diagnostics$repeated_pairs==1L)
raw <- data.frame(pathway=letters[1:5],size=c(14,15,500,501,20),
                  pval=c(.0001,.025,.05,.0001,NA),NES=rep(3,5))
p <- gsea_postprocess(raw)
stopifnot(p$n_tested==2L,p$n_undefined==1L,p$n_significant==0L,
          identical(p$raw$padj_size15[c(2,3)],c(.05,.05)))
boundary <- data.frame(pathway=c("positive","negative"),size=c(15,500),pval=c(.001,.001),NES=c(2,-2))
stopifnot(gsea_postprocess(boundary)$n_significant==2L,
          gsea_postprocess(boundary)$n_significant_NES2==0L)
boundary$NES <- c(2.0001,-2.0001)
stopifnot(gsea_postprocess(boundary)$n_significant_NES2==2L)
bad <- try(gsea_ranking(toy,list(A=NA_character_,B="1")),silent=TRUE)
stopifnot(inherits(bad,"try-error"))
write.table(do.call(rbind,out),file.path(gsea_docs,"reference_rankings.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(do.call(rbind,out))
cat("PASS: 17 historical vector reconstructions, deduplication, deterministic ties and BH boundaries\n")
