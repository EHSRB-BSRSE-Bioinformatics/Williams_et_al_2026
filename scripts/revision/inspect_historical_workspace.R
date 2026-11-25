# Read-only provenance inspection; never evaluate saved analysis functions.
file <- "../2023_Meier_FDR/.RData"
out <- "data/large/reactome_ora_v1/workspace_provenance"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
e <- new.env(parent=baseenv())
message(Sys.time()," inspecting saved 2023-11-09 workspace; bytes=",file.info(file)$size)
problem <- tryCatch({load(file,envir=e); NULL},error=function(error)conditionMessage(error))
if(!is.null(problem)) cat("Workspace load limitation:",problem,"\n")
names <- ls(e,all.names=TRUE)
cat("Relevant saved object names:\n")
print(names[grepl("anno|pathway|session|version|seed|RNG",names,ignore.case=TRUE)])
targets <- intersect(c("pathways","anno","short_anno","p.anno",".Random.seed","sessionInfo","session_info","sessioninfo"),names)
for(nm in targets) {
  x <- get(nm,envir=e,inherits=FALSE)
  cat(nm,": class=",paste(class(x),collapse=",")," length=",length(x),"\n",sep="")
  if(is.data.frame(x)||is.list(x)||is.atomic(x)) saveRDS(x,file.path(out,paste0(gsub("^\\.","",nm),".rds")))
}
if("pathways" %in% targets) {
  original <- get("pathways",e)
  a <- readRDS("data/large/reactome_ora_v2/annotation/frozen_annotation.rds")
  original <- original[grep("R-HSA-",names(original))]
  common <- intersect(names(original),names(a$pathways))
  cat("Saved human pathways:",length(original),"; frozen:",length(a$pathways),
      "; missing:",length(setdiff(names(original),names(a$pathways))),
      "; added:",length(setdiff(names(a$pathways),names(original))),
      "; membership disagreements:",sum(!vapply(common,function(id)setequal(original[[id]],a$pathways[[id]]),logical(1))),"\n")
}
if("short_anno" %in% targets) {
  original <- as.data.frame(get("short_anno",e))
  a <- readRDS("data/large/reactome_ora_v2/annotation/frozen_annotation.rds")
  pairs <- function(z) unique(paste(z$ALIAS,z$ENTREZID,sep="\t"))
  cat("Saved alias/Entrez unique pairs match frozen:",setequal(pairs(original),pairs(a$short_annotation)),"\n")
  cat("Saved annotation rows:",nrow(original),"; frozen:",nrow(a$short_annotation),"\n")
}
if(".Random.seed" %in% targets) {
  seed <- get(".Random.seed",e)
  cat("Saved RNG state length:",length(seed),"; kind code:",seed[1],". This is a saved end-state, not a documented GSEA initialization/parallel stream.\n")
}
cat("No saved functions or enrichment calculations were executed.\n")
