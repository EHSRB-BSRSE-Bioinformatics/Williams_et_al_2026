source("R/revision/reactome_ora.R")
root <- "data/large/reactome_ora_v1"
compatible <- "--org3.17" %in% commandArgs(TRUE)
a <- readRDS(file.path(root, if(compatible) "annotation/frozen_annotation_org3.17.rds" else "annotation/frozen_annotation.rds"))
data.table::setDTthreads(1L)
ref <- readRDS("data/reference/deg_reference.rds")
compiled <- compile_reactome(a$pathways)
stems <- c(TK6="tk6", `MCF-7`="mcf7", HepaRG="heparg", `Liver spheroids`="spheroids")
metadata <- trial_identity(names(readRDS("data/derived/monte_carlo_trial_metadata.rds")))
wp <- readLines("data/external/wikipathways-20231010-gmt-Homo_sapiens.gmt")
wpgenes <- unique(unlist(lapply(strsplit(wp,"\t"),function(x) x[-c(1,2)])))
gogenes <- NULL
if (compatible) {
  con <- DBI::dbConnect(RSQLite::SQLite(),file.path(root,"annotation/candidate_3.17/org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite"))
  gogenes <- DBI::dbGetQuery(con,"SELECT DISTINCT genes.gene_id FROM genes JOIN go_bp_all ON genes._id=go_bp_all._id")$gene_id
  DBI::dbDisconnect(con)
}
rows <- list()
for (cell in names(stems)) {
  ids <- metadata$run[metadata$cell_line == cell]
  ora <- readRDS(file.path("../2023_Meier_FDR/bootstrap_ORA/2023-11-09",paste0(stems[[cell]],"_ora.RDS")))
  for (id in intersect(ids,names(ref))) {
    z <- historical_ora_input(ref[[id]],a$alias_map)
    # Literal original data.table selection independently verifies representatives.
    d <- data.table::data.table(as.data.frame(ref[[id]]$DESeqRes),keep.rownames="rn")
    d <- d[, .(log2FoldChange,rn,padj)]
    d[,rn:=sub("_(.*)","",rn)]
    d <- d[, .SD[which.max(abs(log2FoldChange))],by=rn]
    d <- d[d$rn %in% a$short_annotation$ALIAS,]
    x <- merge(as.data.frame(d),a$short_annotation,by.x="rn",by.y="ALIAS",all.x=TRUE,sort=FALSE)
    linear <- ifelse(x$log2FoldChange>0,2^x$log2FoldChange,-1/(2^x$log2FoldChange))
    sel <- unique(x$ENTREZID[!is.na(x$padj) & x$padj < .05 & abs(linear)>1.5])
    stopifnot(setequal(sel,z$selected),setequal(x$ENTREZID,z$background))
    direct <- reactome_hypergeom(z,a$pathways)
    fast <- reactome_hypergeom(z,a$pathways,compiled=compiled)
    stopifnot(isTRUE(all.equal(direct,fast)),identical(fast,reactome_hypergeom(z,a$pathways,compiled=compiled)))
    old <- ora[[match(id,ids)]]
    for (db in c("enrich_go_bp","enricher_wikipathways")) {
      obj <- old[[db]]
      if (isS4(obj)) {
        bg <- methods::slot(obj,"universe")
        genes <- methods::slot(obj,"gene")
        rows[[length(rows)+1L]] <- data.frame(run=id,db=db,
          n_selected=length(z$selected), n_archived=length(genes),
          selected_equal=setequal(z$selected,genes),
          background_archived_subset=all(bg %in% z$background),
          background_equal=if(db=="enricher_wikipathways") setequal(intersect(z$background,wpgenes),bg) else if(compatible) setequal(intersect(z$background,gogenes),bg) else NA,
          selected_missing=paste(setdiff(genes,z$selected),collapse=";"),
          selected_added=paste(setdiff(z$selected,genes),collapse=";"))
      }
    }
    cat(id," selected=",length(z$selected)," background=",length(z$background)," annotated=",fast$n_background_annotated," family=",fast$n_tested," status=",fast$status,"\n",sep="")
  }
  rm(ora); gc()
}
out <- do.call(rbind,rows)
write.table(out,if(compatible) "docs/revision/reactome_ora_v1/pilot_archive_mapping_org3.17.tsv" else "docs/revision/reactome_ora_v1/pilot_archive_mapping.tsv",sep="\t",row.names=FALSE,quote=FALSE)
print(out)
stopifnot(all(out$background_archived_subset),all(is.na(out$background_equal)|out$background_equal))
if ("--allow-known-annotation-drift" %in% commandArgs(TRUE)) {
  stopifnot(all(out$selected_missing==""),all(out$selected_added %in% c("","84215")),sum(!out$selected_equal)==4L)
  cat("Initial org3.16 pilot: known single-ID annotation discrepancy in two trials; candidate release must be checked separately.\n")
} else stopifnot(all(out$selected_equal))
