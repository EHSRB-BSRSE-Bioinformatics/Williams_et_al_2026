# Bounded collection compatibility: compare frozen membership-derived sizes with
# a few archived GSEA reference tables, without invoking enrichment software.
library(DBI)
library(RSQLite)
source("R/revision/reactome_ora.R")
dir <- "data/large/reactome_ora_v1/annotation"
a <- readRDS(file.path(dir,"frozen_annotation_org3.17.rds"))
con <- dbConnect(SQLite(),file.path(dir,"candidate_reactome1.84/reactome.db/inst/extdata/reactome.sqlite"))
links <- dbGetQuery(con,"SELECT DB_ID,gene_id FROM pathway2gene WHERE DB_ID LIKE 'R-HSA-%'")
b <- a; b$pathways <- split(as.character(links$gene_id),links$DB_ID)
b$reactome_metadata <- dbGetQuery(con,"SELECT * FROM metadata")
b$versions["reactome.db"] <- "1.84.0"
b$path_names <- dbGetQuery(con,"SELECT path_name,DB_ID FROM pathway2name")
dbDisconnect(con)
b$reactome_provenance <- "reactome.db 1.84.0 source SQLite pathway2gene, human R-HSA identifiers; bounded compatibility check against archived GSEA effective sizes"
candidate <- file.path(dir,"frozen_annotation_reactome1.84_org3.17.rds")
if(file.exists(candidate)) stopifnot(identical(readRDS(candidate),b)) else saveRDS(b,candidate,version=3)
refs <- readRDS("data/reference/deg_reference.rds")
files <- list.files("data/large/reactome_ora_v1",pattern="^gsea_reference_.*[.]rds$",full.names=TRUE)
g <- unlist(lapply(files,readRDS),recursive=FALSE)
result <- list()
for(version in c("1.77.0","1.84.0")) {
  paths <- if(version=="1.77.0") a$pathways else b$pathways
  for(id in names(g)) {
    input <- historical_ora_input(refs[[id]],a$alias_map)
    sizes <- vapply(paths,function(z)length(intersect(z,input$background)),integer(1))
    observed <- as.data.frame(g[[id]])
    predicted <- names(sizes)[sizes>=1L&sizes<=500L]
    common <- intersect(predicted,observed$pathway)
    eligible <- names(sizes)[sizes>=15L&sizes<=500L]
    observed_eligible <- observed$pathway[observed$size>=15L&observed$size<=500L]
    result[[length(result)+1L]] <- data.frame(version=version,run=id,
      predicted_raw=length(predicted),observed_raw=nrow(observed),
      missing_paths=length(setdiff(observed$pathway,predicted)),added_paths=length(setdiff(predicted,observed$pathway)),
      shared=length(common),size_disagreements=sum(sizes[common]!=observed$size[match(common,observed$pathway)]),
      leading_edge_violations=sum(vapply(seq_len(nrow(observed)),function(j) !all(observed$leadingEdge[[j]] %in% paths[[observed$pathway[j]]]),logical(1))),
      eligible_jaccard=length(intersect(eligible,observed_eligible))/length(union(eligible,observed_eligible)))
  }
}
result <- do.call(rbind,result)
write.table(result,"docs/revision/reactome_ora_v1/collection_reference_compatibility.tsv",sep="\t",row.names=FALSE,quote=FALSE)
print(result,row.names=FALSE); print(b$reactome_metadata)
cat("Candidate human pathways:",length(b$pathways),"\n")
