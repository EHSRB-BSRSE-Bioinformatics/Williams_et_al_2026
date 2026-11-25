source("R/config.R")
source("R/summaries.R")
out <- "docs/revision/gsea_unique_ids_v1"
trials <- readRDS(file.path(out,"trial_results.rds"))
expected <- readRDS(file.path(out,"summary.rds"))
stopifnot(nrow(trials)==32000L,!anyDuplicated(trials$run),
  setequal(trials$run,names(readRDS(derived_files$trial_metadata))),
  all(table(trials$cell_line,trials$n)==1000L),
  all(trials$status %in% c("significant","none","untestable")))
x <- read_summary_inputs()
ix <- match(x$gsea$run,trials$run)
stopifnot(identical(as.integer(x$gsea$p_lt_0.05),as.integer(trials$n_significant[ix])),
  identical(as.integer(x$gsea$p_lt_0.05_ES_gt_2),as.integer(trials$n_significant_NES2[ix])))
rates <- x$gsea %>% group_by(cell_line,n) %>% summarise(n_trials=n(),
  rate=mean(contains_false_positive_pval_only),rate_NES2=mean(contains_false_positive_NES_pval_filtered),.groups="drop")
ei <- match(paste(rates$cell_line,rates$n),paste(expected$cell_line,expected$n))
stopifnot(all(rates$n_trials==1000L),all(rates$rate==expected$gsea_rate[ei]),
  all(rates$rate_NES2==expected$gsea_NES2_rate[ei]))
# Frozen publication summaries and unrelated figures must remain byte-identical.
preserved <- c(list.files("data/derived",pattern="[.]rds$",full.names=TRUE),
  list.files("figures",pattern="[.](pdf|png)$",full.names=TRUE))
preserved <- preserved[!grepl("/gsea_(monte_carlo_false_positive_rates|archived_duplicate_ids_supplementary)[.]",preserved)]
before <- tools::md5sum(preserved)
status <- system2(file.path(R.home("bin"),"Rscript"),
  c("scripts/render_figures.R","--only=gsea,gsea_archived_supplementary"))
stopifnot(status==0L,identical(before,tools::md5sum(preserved)))
write.table(rates,file.path(out,"figure_rates.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"promotion_sessionInfo.txt"))
cat("PASS: 64 corrected plotted rates, both raw-derived count endpoints, archived summary and unrelated artifacts preserved\n")
