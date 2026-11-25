# Cache only the archived Figure 4A family p-values needed by the new sensitivity.
# This is an extraction, with no enrichment calculation or full-summary audit.
root <- "data/large/reactome_ora_v1"
metadata <- names(readRDS("data/derived/monte_carlo_trial_metadata.rds"))
stems <- c(TK6="tk6",`MCF-7`="mcf7",HepaRG="heparg",`Liver spheroids`="spheroids")
args <- commandArgs(TRUE)
if(length(args)) stems <- stems[stems %in% args]
for(cell in names(stems)) {
  target <- file.path(root,paste0("gsea_family_",stems[[cell]],".rds"))
  if(file.exists(target)) next
  message(Sys.time()," extracting archived GSEA comparison inputs: ",cell)
  g <- readRDS(file.path("../2023_Meier_FDR/bootstrap_GSEA/2023-11-08",paste0(stems[[cell]],"_gsea.RDS")))
  ids <- metadata[startsWith(metadata,paste0(cell,"_"))]
  stopifnot(length(g)==length(ids))
  if(!is.null(names(g))) stopifnot(identical(names(g),ids))
  names(g) <- ids
  refids <- intersect(ids,names(readRDS("data/reference/deg_reference.rds")))
  saveRDS(lapply(g[refids],as.data.frame),file.path(root,paste0("gsea_reference_",stems[[cell]],".rds")))
  slim <- lapply(g,function(z) {
    z <- as.data.frame(z)
    z[!is.na(z$pval)&z$size>=15L&z$size<=500L,c("pathway","size","pval")]
  })
  tmp <- paste0(target,".tmp"); saveRDS(slim,tmp); stopifnot(file.rename(tmp,target))
  rm(g,slim); gc()
}
