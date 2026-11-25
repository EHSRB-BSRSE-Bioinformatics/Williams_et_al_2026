source("R/config.R")
suppressPackageStartupMessages(source("R/summaries.R"))
suppressPackageStartupMessages(source("R/figures.R"))
x<-read_deg_figure_summary()
frozen<-readRDS(derived_files$trial_metadata)
stopifnot(nrow(x)==96000L,!anyNA(x$nDEGs),!anyNA(x$contains_false_positive))
for(label in c("P-value filter only","P-value + fold-change","Final DEGs")) {
  d<-x[x$Filter==label,]
  stopifnot(nrow(d)==32000L,!anyDuplicated(d$name),setequal(d$name,names(frozen)),
            all(table(d$cell_line,d$n_per_group)==1000L))
  for(i in split(seq_len(nrow(d)),interaction(d$cell_line,d$n_per_group,drop=TRUE))) {
    stopifnot(setequal(d$rank[i],1:1000),all(diff(d$nDEGs[i[order(d$rank[i])]])>=0))
  }
}
a<-x[x$Filter=="P-value filter only",];b<-x[x$Filter=="Final DEGs",]
stopifnot(sum(a$contains_false_positive)==5316L,sum(b$contains_false_positive)==2202L,
          a$nDEGs[a$name=="Liver spheroids_3_3_334"]==3L,
          b$nDEGs[b$name=="Liver spheroids_3_3_334"]==2L,
          a$nDEGs[a$name=="HepaRG_3_3_5"]==23L)
middle<-x[x$Filter=="P-value + fold-change",]
stopifnot(sum(middle$contains_false_positive)==2231L,
          all(middle$nDEGs<=a$nDEGs),all(b$nDEGs<=middle$nDEGs))
plot<-build_deg_figure(x)
stopifnot(length(plot[[1]]$patches$plots)==1L,length(plot[[2]]$patches$plots)==1L)
fills<-lapply(1:2,function(i)ggplot_build(plot[[1]][[i]])$plot$scales$get_scales("fill")$get_limits())
stopifnot(identical(fills[[1]],fills[[2]]))
fixed<-readRDS("data/derived/deg_thresholds/threshold_sensitivity.rds")
for(j in 1:2) {
  panel<-plot[[2]][[j]]
  expected<-fixed$summary[fixed$summary$endpoint=="screening" & fixed$summary$padj_cutoff==.05 &
    fixed$summary$fc_label==c("none","1.5")[j],]
  actual<-panel$data
  key<-paste(actual$cell_line,actual$n_per_group)
  i<-match(key,paste(expected$cell,expected$n))
  stopifnot(!anyNA(i),all(actual$n_trials_observed==1000L),
            identical(actual$false_positive_rate,expected$rate[i]))
}
legacy<-readRDS(derived_files$deg)
stopifnot(inherits(tryCatch(build_deg_figure(legacy),error=identity),"error"))
cat("All three DEG stages, both plotted panels, trial IDs, denominators and ranks verified. Legacy relabeling rejected.\n")
