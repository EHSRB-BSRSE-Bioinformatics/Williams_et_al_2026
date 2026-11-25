# Reconstruct all decisions from saved raw p-values; pair exact frozen IDs.
source("R/config.R")
source("R/summaries.R")
source("R/revision/reactome_ora.R")
source("R/revision/gsea_unique_ids.R")
ids <- names(readRDS(derived_files$trial_metadata))
identity <- trial_identity(ids)
archived <- readRDS(derived_files$gsea)
ora <- readRDS("docs/revision/reactome_ora_v2/trial_results.rds")
source_archives <- read.delim(file.path(gsea_docs,"source_archives.tsv"))
pilot <- readRDS(file.path(gsea_root,"pilot.rds"))
stopifnot(isTRUE(pilot$validated))
stopifnot(!anyDuplicated(archived$run),!anyDuplicated(ora$run),
  setequal(ids,archived$run),setequal(ids,ora$run))
rows <- list(); compact <- list()
manifest <- list()
for(cell in names(gsea_stems)) for(n in 3:10) {
  stem <- gsea_stems[[cell]]
  current <- identity$run[identity$cell_line==cell & identity$n==n]
  file <- file.path(gsea_root,"full",sprintf("%s_%02d.rds",stem,n))
  rank_file <- file.path(gsea_root,"rankings",sprintf("%s_%02d.rds",stem,n))
  z <- readRDS(file)
  ranks <- readRDS(rank_file)
  old <- readRDS(file.path("data/large/reactome_ora_v2/full",sprintf("%s_%02d.rds",stem,n)))$rows
  stopifnot(identical(names(z$rows),current),identical(names(ranks$rows),current),
    unname(ranks$archive_md5)==source_archives$md5[match(paste0(stem,".RDS"),source_archives$archive)],
    identical(z$provenance$code_md5,tools::md5sum(c("R/revision/gsea_unique_ids.R","scripts/revision/run_gsea_unique_ids.R"))),
    identical(z$provenance,ranks$provenance),identical(z$ranking_md5,tools::md5sum(rank_file)))
  manifest[[length(manifest)+1L]] <- data.frame(file=c(file,rank_file),md5=unname(tools::md5sum(c(file,rank_file))))
  for(id in current) {
    r <- z$rows[[id]]$result
    if(id %in% names(pilot$rows)) stopifnot(identical(r$raw,pilot$rows[[id]]$unique$raw))
    stopifnot(r$status %in% c("significant","none","untestable"))
    rebuilt <- gsea_postprocess(r$raw)
    for(field in c("raw","status","n_tested","n_eligible","n_undefined","n_undefined_all","n_significant","n_significant_NES2"))
      stopifnot(identical(rebuilt[[field]],r[[field]]))
    rank <- ranks$rows[[id]]$unique
    stopifnot(!anyDuplicated(names(rank)),!anyNA(names(rank)),all(nzchar(names(rank))),all(is.finite(rank)),
      setequal(names(rank),old[[id]]$input$background),
      length(rank)==old[[id]]$result$n_background,
      identical(z$rows[[id]]$diagnostics,ranks$rows[[id]]$diagnostics),
      r$seed==20261116L+match(id,ids))
    ot <- old[[id]]$result$terms
    gt <- r$raw[!is.na(r$raw$padj_size15),]
    common <- intersect(ot$pathway,gt$pathway)
    oi <- match(common,ot$pathway); gi <- match(common,gt$pathway)
    og <- old[[id]]$result$status!="insufficient_input"
    old_count <- ora$n_significant[match(id,ora$run)]
    stopifnot(old_count==old[[id]]$result$n_significant)
    rows[[id]] <- data.frame(run=id,status=r$status,as.data.frame(z$rows[[id]]$diagnostics),
      n_eligible=r$n_eligible,n_undefined=r$n_undefined,n_undefined_all=r$n_undefined_all,
      n_tested=r$n_tested,n_significant=r$n_significant,n_significant_NES2=r$n_significant_NES2,
      elapsed=r$elapsed,seed=r$seed,n_warnings=length(r$warnings),
      n_common=length(common),family_jaccard=length(common)/length(union(ot$pathway,gt$pathway)),
      size_disagreements=sum(ot$size[oi]!=gt$size[gi]),
      ora_common=if(og)sum(p.adjust(ot$pvalue[oi],"BH")<.05) else 0L,
      gsea_common=sum(p.adjust(gt$pval[gi],"BH")<.05),
      ora_significant=old_count,
      archived_significant=archived$p_lt_0.05[match(id,archived$run)],
      archived_significant_NES2=archived$p_lt_0.05_ES_gt_2[match(id,archived$run)],
      archived_n_tested=archived$num_rows[match(id,archived$run)])
  }
  # Independent public summary builder verifies both plotted endpoints.
  s <- summarise_gsea_from_trials(lapply(z$rows,function(r)r$result$raw))
  stopifnot(identical(as.integer(s$p_lt_0.05),unname(vapply(z$rows,function(r)r$result$n_significant,integer(1)))),
    identical(as.integer(s$p_lt_0.05_ES_gt_2),unname(vapply(z$rows,function(r)r$result$n_significant_NES2,integer(1)))))
  compact[[paste(cell,n)]] <- s
  message("Validated raw p-values, rankings and summary: ",cell," n=",n)
}
results <- do.call(rbind,rows)
results <- merge(identity,results,by="run",sort=FALSE)
results <- results[match(ids,results$run),]
stopifnot(identical(results$run,ids),nrow(results)==32000L,!anyDuplicated(results$run),
  all(table(results$cell_line,results$n)==1000L),!any(results$size_disagreements))
results$changed_positive <- (results$n_significant>0)!=(results$archived_significant>0)
results$lost_positive <- results$n_significant==0 & results$archived_significant>0
results$gained_positive <- results$n_significant>0 & results$archived_significant==0
results$concordance <- with(results,ifelse(ora_significant>0,ifelse(n_significant>0,"Both","ORA only"),
  ifelse(n_significant>0,"GSEA only","Neither")))
saveRDS(results,file.path(gsea_docs,"trial_results.rds"))
write.table(results,file.path(gsea_docs,"trial_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
pilot_comparison <- do.call(rbind,lapply(pilot$rows,function(p) {
  archived_count <- results$archived_significant[match(p$run,results$run)]
  data.frame(run=p$run,unique_significant=p$unique$n_significant,
    duplicate_same_env_significant=p$duplicated$n_significant,archived_significant=archived_count,
    input_delta=p$unique$n_significant-p$duplicated$n_significant,
    environment_rng_delta=p$duplicated$n_significant-archived_count)
}))
write.table(pilot_comparison,file.path(gsea_docs,"pilot_archive_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
compact <- bind_rows(compact)
compact <- compact[match(ids,compact$run),]
saveRDS(compact,file.path(gsea_docs,"gsea_summary.rds"))
set.seed(20261116L)
summary <- results %>% group_by(cell_line,n) %>% group_modify(function(d,k) {
  delta <- as.integer(d$ora_significant>0)-as.integer(d$n_significant>0)
  boot <- rmultinom(10000,nrow(d),tabulate(delta+2L,nbins=3L)/nrow(d))
  ci <- quantile((boot[3,]-boot[1,])/nrow(d),c(.025,.975),names=FALSE)
  data.frame(n_trials=nrow(d),gsea_rate=mean(d$n_significant>0),
    gsea_NES2_rate=mean(d$n_significant_NES2>0),archived_rate=mean(d$archived_significant>0),
    ora_rate=mean(d$ora_significant>0),ora_minus_gsea=mean(delta),ci_lower=ci[1],ci_upper=ci[2],
    lost_positive=sum(d$lost_positive),gained_positive=sum(d$gained_positive),
    both=sum(d$concordance=="Both"),ora_only=sum(d$concordance=="ORA only"),
    gsea_only=sum(d$concordance=="GSEA only"),neither=sum(d$concordance=="Neither"),
    untestable=sum(d$status=="untestable"),trials_undefined=sum(d$n_undefined>0),
    undefined_pathways=sum(d$n_undefined),family_median=median(d$n_tested),
    family_loss_median=median(d$archived_n_tested-d$n_tested),
    count_mean=mean(d$n_significant),count_median=median(d$n_significant),
    count_q95=unname(quantile(d$n_significant,.95)),count_max=max(d$n_significant),
    ora_common_rate=mean(d$ora_common>0),gsea_common_rate=mean(d$gsea_common>0),
    common_median=median(d$n_common),family_jaccard_median=median(d$family_jaccard),
    unique_ids_median=median(d$n_unique),repeated_pairs_median=median(d$repeated_pairs),
    multi_score_ids_median=median(d$multi_score_ids),multi_score_trials=sum(d$multi_score_ids>0),
    tied_nonzero_unique_median=median(d$tied_nonzero_unique),
    tied_nonzero_historical_median=median(d$tied_nonzero_historical))
}) %>% ungroup()
saveRDS(summary,file.path(gsea_docs,"summary.rds"))
write.table(summary,file.path(gsea_docs,"summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
distribution <- bind_rows(
  results %>% transmute(cell_line,n,method="Unique-ID GSEA",count=n_significant),
  results %>% transmute(cell_line,n,method="Archived GSEA",count=archived_significant),
  results %>% transmute(cell_line,n,method="Reactome ORA",count=ora_significant)) %>%
  count(cell_line,n,method,count,name="n_trials") %>% mutate(proportion=n_trials/1000)
write.table(distribution,file.path(gsea_docs,"count_distribution.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,manifest),file.path(gsea_docs,"checkpoint_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(as.data.frame(summary),row.names=FALSE)
cat("PASS: 32,000 frozen unique IDs, every ranking unique/finite, raw BH/count reconstruction, common families and 1,000/group\n")
cat("Overall: unique GSEA=",sum(results$n_significant>0)," archived=",sum(results$archived_significant>0),
  " ORA=",sum(results$ora_significant>0)," lost=",sum(results$lost_positive)," gained=",sum(results$gained_positive),"\n")
