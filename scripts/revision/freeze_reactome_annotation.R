dir <- "data/large/reactome_ora_v1/annotation"
lib <- file.path(dir, "library")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("org.Hs.eg.db_3.16.0", "reactome.db_1.77.0")) {
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", paste0("--library=", lib), file.path(dir, paste0(pkg, ".tar.gz"))))
  stopifnot(status == 0L)
}
.libPaths(c(normalizePath(lib), .libPaths()))
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reactome.db)
keys <- keys(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, keys = keys, columns = c("SYMBOL", "ENTREZID", "ALIAS", "ENSEMBL"))
# Retain expanded rows, matching the original source (including ENSEMBL repeats).
short <- anno[!is.na(anno$ALIAS) & !is.na(anno$ENTREZID), c("ENTREZID", "ALIAS")]
pathways <- as.list(reactomePATHID2EXTID)
pathways <- pathways[grep("R-HSA-", names(pathways))]
frozen <- list(alias_map = split(as.character(short$ENTREZID), short$ALIAS),
               short_annotation = short, pathways = pathways,
               path_names = as.data.frame(reactomePATHNAME2ID),
               org_metadata = metadata(org.Hs.eg.db), reactome_metadata = metadata(reactome.db),
               versions = sapply(c("AnnotationDbi", "org.Hs.eg.db", "reactome.db"), function(p) as.character(packageVersion(p))))
saveRDS(frozen, file.path(dir, "frozen_annotation.rds"), version = 3)
capture.output(sessionInfo(), file = file.path(dir, "freeze_sessionInfo.txt"))
print(frozen$versions)
print(frozen$reactome_metadata)
cat("Human pathways:", length(pathways), "\n")
