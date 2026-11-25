# Compare extracted plain historical data; no enrichment or saved function calls.
root <- "data/large/reactome_ora_v1/workspace_provenance"
a <- readRDS("data/large/reactome_ora_v2/annotation/frozen_annotation.rds")
aliases <- readRDS(file.path(root,"short_anno.rds"))
pathways <- readRDS(file.path(root,"pathways.rds"))
pathways <- pathways[grep("R-HSA-",names(pathways))]
stopifnot(identical(as.data.frame(aliases),a$short_annotation),setequal(names(pathways),names(a$pathways)))
stopifnot(all(vapply(names(pathways),function(id)setequal(pathways[[id]],a$pathways[[id]]),logical(1))))
cat("Saved alias rows/order identical to frozen: TRUE; rows=",nrow(aliases),"\n",sep="")
cat("Saved/frozen human pathway IDs and memberships identical: TRUE; sets=",length(pathways),"\n",sep="")
cat("Pathway-list order identical:",identical(names(pathways),names(a$pathways)),"(irrelevant to hypergeometric/BH decisions)\n")
assign(".Random.seed",readRDS(file.path(root,"Random.seed.rds")),envir=.GlobalEnv)
cat("Saved-state RNG kinds:",paste(RNGkind(),collapse=" / "),"\n")
cat("The saved state does not identify initial GSEA seeds or parallel worker streams.\n")
