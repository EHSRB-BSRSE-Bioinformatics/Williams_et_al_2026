source("R/revision/reactome_ora.R")
out <- Sys.getenv("REACTOME_ORA_DOCS","docs/revision/reactome_ora_v2")
root <- Sys.getenv("REACTOME_ORA_ROOT","data/large/reactome_ora_v2")
d <- readRDS(file.path(out,"trial_results.rds"))
s <- as.data.frame(readRDS(file.path(out,"summary.rds")))
check_trial_coverage(d,names(readRDS("data/derived/monte_carlo_trial_metadata.rds")))
stopifnot(nrow(d)==32000L,nrow(s)==32L,all(s$n_trials==1000L),
 all(d$n_significant[d$status=="insufficient_input"]==0L),
 all(s$both+s$ora_only+s$gsea_only+s$neither==1000L),
 all(abs(s$difference-(s$ora_rate-s$gsea_rate))<1e-12),
 all(abs(s$both_proportion+s$ora_only_proportion+s$gsea_only_proportion+s$neither_proportion-1)<1e-12),
 all((d$status=="significant")== (d$n_significant>0)),
 all(d$n_selected_annotated<=d$n_selected),all(d$n_background_annotated<=d$n_background))
for (i in seq_len(nrow(s))) {
  z <- d[d$cell_line==s$cell_line[i]&d$n==s$n[i],]
  stopifnot(mean(z$n_significant>0)==s$ora_rate[i],mean(z$gsea_significant>0)==s$gsea_rate[i],
    sum(z$concordance=="Both")==s$both[i])
}
# Reconstruct every significant-count endpoint from saved unfiltered p-values.
for (p in list.files(file.path(root,"full"),pattern="\\.rds$",full.names=TRUE)) {
  x <- readRDS(p)
  for (z in x$rows) {
    r <- z$result
    stopifnot(identical(r$terms$padj,p.adjust(r$terms$pvalue,"BH")),
      !anyDuplicated(r$terms$pathway),all(is.finite(r$terms$pvalue)),all(r$terms$pvalue>=0&r$terms$pvalue<=1),
      all(z$input$selected %in% z$input$background),
      r$n_tested==nrow(r$terms),all(r$terms$size>=15L&r$terms$size<=500L),
      r$n_significant==if(r$status=="insufficient_input") 0L else sum(r$terms$padj<.05),
      d$n_significant[match(z$run,d$run)]==r$n_significant)
  }
}
common_path <- file.path(out,"common_family_trial_results.rds")
if(file.exists(common_path)) {
  c <- readRDS(common_path)
  cs <- readRDS(file.path(out,"common_family_summary.rds"))
  stopifnot(setequal(c$run,d$run),!anyDuplicated(c$run),all(table(c$cell_line,c$n)==1000L),
    all(c$n_common<=c$ora_family&c$n_common<=c$gsea_family),
    all(c$ora_common>=0&c$ora_common<=c$n_common),all(c$gsea_common>=0&c$gsea_common<=c$n_common))
  for(i in seq_len(nrow(cs))) {
    z <- c[c$cell_line==cs$cell_line[i]&c$n==cs$n[i],]
    stopifnot(mean(z$ora_common>0)==cs$ora_rate[i],mean(z$gsea_common>0)==cs$gsea_rate[i])
  }
  cat("Common-family identity, denominators, bounds and summary rates verified.\n")
}
cat("32,000 trial IDs, all 32 denominators, concordance partitions and all saved p-value summaries verified.\n")
