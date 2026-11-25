source("R/revision/reactome_ora.R")
library(dplyr)
root <- Sys.getenv("REACTOME_ORA_ROOT","data/large/reactome_ora_v2")
out <- Sys.getenv("REACTOME_ORA_DOCS","docs/revision/reactome_ora_v2")
metadata <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
paths <- list.files(file.path(root,"full"),pattern="\\.rds$",full.names=TRUE)
stopifnot(length(paths)==32L)
rows <- unlist(lapply(paths,function(p) readRDS(p)$rows),recursive=FALSE)
compact <- do.call(rbind,lapply(rows,function(z) {
  r <- z$result
  if (r$status=="failed") stop(z$run,": ",r$error)
  data.frame(run=z$run,status=r$status,n_significant=r$n_significant,
    n_selected=r$n_selected,n_background=r$n_background,
    n_selected_annotated=r$n_selected_annotated,n_background_annotated=r$n_background_annotated,
    n_tested=r$n_tested,n_features=z$input$n_features,n_aliases=z$input$n_aliases,
    n_mapped_aliases=z$input$n_mapped_aliases,n_selected_expanded=z$input$selected_expanded,
    single_gene_raw_BH_positive=r$n_selected_annotated==1L && any(r$terms$padj<.05))
}))
compact <- merge(trial_identity(names(metadata)),compact,by="run",sort=FALSE)
compact <- compact[match(names(metadata),compact$run),]
check_trial_coverage(compact,names(metadata))
g <- as.data.frame(readRDS("data/derived/gsea_monte_carlo_summary.rds"))
o <- as.data.frame(readRDS("data/derived/ora_monte_carlo_summary.rds"))
stopifnot(!anyDuplicated(g$run),!anyDuplicated(o$run),setequal(compact$run,g$run),setequal(compact$run,o$run))
g <- g[match(compact$run,g$run),]; o <- o[match(compact$run,o$run),]
compact$gsea_significant <- g$p_lt_0.05
compact$gsea_significant_NES2 <- g$p_lt_0.05_ES_gt_2
compact$gsea_n_tested <- g$num_rows
compact$go_significant <- o$go_bp
compact$wp_significant <- o$enricher_wikipathways
compact$concordance <- with(compact,ifelse(n_significant>0,ifelse(gsea_significant>0,"Both","ORA only"),ifelse(gsea_significant>0,"GSEA only","Neither")))
saveRDS(compact,file.path(out,"trial_results.rds"))
write.table(compact,file.path(out,"trial_results.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# Fixed paired-record bootstrap. Multinomial draws are exactly the ordinary
# bootstrap for the three-valued paired difference (-1,0,+1).
set.seed(20261116)
summary <- compact %>% group_by(cell_line,n) %>% group_modify(function(d,k) {
  delta <- as.integer(d$n_significant>0)-as.integer(d$gsea_significant>0)
  boot <- rmultinom(10000,nrow(d),tabulate(delta+2L,nbins=3L)/nrow(d))
  diffs <- (boot[3,]-boot[1,])/nrow(d)
  data.frame(n_trials=nrow(d),ora_rate=mean(d$n_significant>0),gsea_rate=mean(d$gsea_significant>0),
    difference=mean(delta),ci_lower=unname(quantile(diffs,.025)),ci_upper=unname(quantile(diffs,.975)),
    both=sum(d$concordance=="Both"),ora_only=sum(d$concordance=="ORA only"),
    gsea_only=sum(d$concordance=="GSEA only"),neither=sum(d$concordance=="Neither"),
    insufficient=sum(d$status=="insufficient_input"),go_rate=mean(d$go_significant>0),wp_rate=mean(d$wp_significant>0),
    single_gene_raw_BH_positive=sum(d$single_gene_raw_BH_positive),
    ora_count_mean=mean(d$n_significant),ora_count_median=median(d$n_significant),ora_count_q95=unname(quantile(d$n_significant,.95)),ora_count_max=max(d$n_significant),
    gsea_count_mean=mean(d$gsea_significant),gsea_count_median=median(d$gsea_significant),gsea_count_q95=unname(quantile(d$gsea_significant,.95)),gsea_count_max=max(d$gsea_significant),
    selected_median=median(d$n_selected),selected_max=max(d$n_selected),
    background_median=median(d$n_background),annotated_background_median=median(d$n_background_annotated),
    ora_family_median=median(d$n_tested),gsea_family_median=median(d$gsea_n_tested))
}) %>% ungroup() %>% mutate(both_proportion=both/n_trials,ora_only_proportion=ora_only/n_trials,gsea_only_proportion=gsea_only/n_trials,neither_proportion=neither/n_trials)
saveRDS(summary,file.path(out,"summary.rds"))
write.table(summary,file.path(out,"summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
counts <- bind_rows(
  compact %>% transmute(cell_line,n,method="Reactome ORA",significant_count=n_significant),
  compact %>% transmute(cell_line,n,method="Reactome GSEA",significant_count=gsea_significant)) %>%
  count(cell_line,n,method,significant_count,name="n_trials") %>% mutate(proportion=n_trials/1000)
write.table(counts,file.path(out,"count_distribution.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
diagnostics <- compact %>% mutate(selected_bin=cut(n_selected_annotated,breaks=c(-1,0,1,9,99,Inf),
  labels=c("0","1","2-9","10-99","100+"))) %>% group_by(cell_line,n,selected_bin) %>%
  summarise(n_trials=n(),ora_rate=mean(n_significant>0),gsea_rate=mean(gsea_significant>0),
    mapped_background_median=median(n_background),reactome_background_median=median(n_background_annotated),
    .groups="drop")
write.table(diagnostics,file.path(out,"input_diagnostics.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
print(as.data.frame(summary)[,1:12],row.names=FALSE)
cat("Overall paired trial rates: ORA=",mean(compact$n_significant>0)," GSEA=",mean(compact$gsea_significant>0),"\n")
