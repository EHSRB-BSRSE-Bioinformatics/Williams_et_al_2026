# Resolve the pilot discrepancy using an archived adjacent annotation release.
library(AnnotationDbi)
source("R/revision/reactome_ora.R")
root <- "data/large/reactome_ora_v1/annotation"
a <- readRDS(file.path(root,"frozen_annotation.rds"))
db <- loadDb(file.path(root,"candidate_3.17/org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite"))
anno <- select(db,keys=keys(db),columns=c("SYMBOL","ENTREZID","ALIAS","ENSEMBL"))
short <- anno[!is.na(anno$ALIAS)&!is.na(anno$ENTREZID),c("ENTREZID","ALIAS")]
b <- a; b$short_annotation <- short
b$alias_map <- split(as.character(short$ENTREZID),short$ALIAS)
b$org_metadata <- metadata(db); b$versions["org.Hs.eg.db"] <- "3.17.0"
b$mapping_provenance <- "org.Hs.eg.db 3.17.0 source SQLite loaded using AnnotationDbi 1.60.0; candidate chosen by pilot slot concordance, not a recovered runtime record"
saveRDS(b,file.path(root,"frozen_annotation_org3.17.rds"),version=3)
ref <- readRDS("data/reference/deg_reference.rds")
old <- read.delim("docs/revision/reactome_ora_v1/pilot_archive_mapping.tsv",check.names=FALSE,colClasses="character",na.strings="NA")
wp <- unique(unlist(lapply(strsplit(readLines("data/external/wikipathways-20231010-gmt-Homo_sapiens.gmt"),"\t"),function(x)x[-c(1,2)])))
out <- lapply(seq_len(nrow(old)),function(i) {
  id <- old$run[i]
  za <- historical_ora_input(ref[[id]],a$alias_map)
  zb <- historical_ora_input(ref[[id]],b$alias_map)
  expected <- union(setdiff(za$selected,strsplit(old$selected_added[i],";",fixed=TRUE)[[1]]),strsplit(old$selected_missing[i],";",fixed=TRUE)[[1]])
  expected <- setdiff(expected,"")
  data.frame(run=id,db=old$db[i],selected_equal=setequal(zb$selected,expected),
    wp_universe_equal=if(old$db[i]=="enricher_wikipathways") setequal(intersect(za$background,wp),intersect(zb$background,wp)) else NA,
    n_selected=length(zb$selected),selected_missing=paste(setdiff(expected,zb$selected),collapse=";"),
    selected_added=paste(setdiff(zb$selected,expected),collapse=";"))
})
out <- do.call(rbind,out)
write.table(out,"docs/revision/reactome_ora_v1/pilot_annotation_candidate.tsv",sep="\t",row.names=FALSE,quote=FALSE)
print(out); print(b$org_metadata)
stopifnot(all(out$selected_equal),all(is.na(out$wp_universe_equal)|out$wp_universe_equal))
cat("Candidate org.Hs.eg.db 3.17.0 matches all 20 available archived selected-gene slots and all 10 checked WikiPathways universes.\n")
