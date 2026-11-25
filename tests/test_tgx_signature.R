source("scripts/tgx_signature/analyze.R")
fails<-function(expr)inherits(tryCatch({force(expr);NULL},error=identity),"error")
stopifnot(fails(check_ids(c("a","a"),c("a","b"),"test")),
          fails(check_ids("a",c("a","b"),"test")),
          fails(parse_ids("TK6_3_4_1")),
          parse_ids("Liver spheroids_10_10_1000")$n==10L)
r<-data.frame(pvalue=c(.001,.001,.001,.001,.001,NA,.8),
              padj=c(.049,.05,.001,.001,.001,NA,NA),
              log2FoldChange=c(log2(1.5),2,2,2,2,2,2),row.names=letters[1:7])
e<-list(DESeqRes=r,DESeqResDEGs=r["e",,drop=FALSE],
        DESeqjs_filtered_out=list(data.frame(),r["c",,drop=FALSE],r["d",,drop=FALSE]),
        DESeqResFiltered=r[rep(TRUE,7),,drop=FALSE])
s<-select_stages(e)
stopifnot(identical(s$final,"e"),setequal(s$significant,c("a","c","d","e")),
          !"a" %in% s$significant_fc,!"b" %in% s$significant,
          "g" %in% s$tested,!"g" %in% s$bh_eligible,
          identical(unname(gene_members(c("p2","p3"),list(c("p1","p2"),"p3",character()))),c(TRUE,TRUE,FALSE)))
e$DESeqResFiltered<-NULL;stopifnot(identical(select_stages(e),s))
e$DESeqResDEGs<-r[c("c","e"),,drop=FALSE];stopifnot(fails(select_stages(e)))
# A trial with no significant genes legitimately has no exclusion tables.
e<-list(DESeqRes=r["b",,drop=FALSE],DESeqResDEGs=data.frame(),DESeqjs_filtered_out=list())
stopifnot(length(select_stages(e)$final)==0L)
cat("All TGx signature unit tests passed.\n")

reference<-readRDS("data/reference/deg_reference.rds")
for(e in reference) {
  stages<-select_stages(e)
  stopifnot(setequal(stages$final,rownames(e$DESeqResDEGs)))
}
cat("All 17 portable archived reference trials match reconstructed final calls.\n")

ref<-read.delim("data/reference/tgx_ddi_64_gene_map.tsv")
probe_sets<-lapply(ref$reference_probes,function(p)if(is.na(p))character() else strsplit(p,";",fixed=TRUE)[[1]])
for(e in reference)for(features in select_stages(e)) {
  slow<-vapply(probe_sets,function(p)any(p %in% features),logical(1))
  stopifnot(identical(gene_members(features,probe_sets),slow))
}
cat("Optimized membership agrees with direct ID matching for all reference stages.\n")
