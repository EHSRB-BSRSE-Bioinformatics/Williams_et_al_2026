source("scripts/deg_thresholds/analyze.R")
fails<-function(expr)inherits(tryCatch({force(expr);NULL},error=identity),"error")
p<-c(.009,.01,.024,.025,.049,.05,NA,NaN,.001,.001,.001,.001)
l<-c(-log2(1.5),log2(1.5),-1,1,0,1,1,1,NA,Inf,-Inf,NaN)
stopifnot(identical(select_threshold(p,l,.01),c(1L,9L,10L,11L,12L)),
          identical(select_threshold(p,l,.025,1.5),c(3L,10L,11L)),
          identical(select_threshold(p,l,.05,1.5),c(3L,4L,10L,11L)),
          5L %in% select_threshold(p,l,.05),
          !5L %in% select_threshold(p,l,.05,1.25),
          fails(select_threshold(p,l,.1)),fails(select_threshold(p,l,.05,1)),
          fails(check_ids(c("a","a"),c("a","b"),"test")))
r<-data.frame(padj=c(.001,.02,.04,.049,.05,NA),log2FoldChange=c(-2,1,.1,NA,2,2),row.names=letters[1:6])
e<-list(DESeqRes=r,DESeqResDEGs=r["a",,drop=FALSE],
        DESeqjs_filtered_out=list(data.frame(),r["b",,drop=FALSE],r["c",,drop=FALSE]),
        DESeqfiltered_table=data.frame(quantile_filtered=1L,spike_filtered=1L),
        DESeqResFiltered=r)
x<-extract_envelope(e);g<-threshold_grid();v<-count_grid(x$features,g)
stopifnot(v[g$padj_cutoff==.05 & g$fc_label=="none","screening"]==4L,
          v[g$padj_cutoff==.05 & g$fc_label=="none","full_rodaf"]==2L,
          v[g$padj_cutoff==.05 & g$fc_label=="1.25","screening"]==2L,
          v[g$padj_cutoff==.05 & g$fc_label=="1.25","full_rodaf"]==1L)
e$DESeqResFiltered<-NULL;stopifnot(identical(x,extract_envelope(e)))
e$DESeqResDEGs<-r["b",,drop=FALSE];stopifnot(fails(extract_envelope(e)))
refs<-readRDS("data/reference/deg_reference.rds")
for(e in refs) {
  x<-extract_envelope(e);v<-count_grid(x$features,g);k<-which(g$padj_cutoff==.05 & g$fc_label=="1.5")
  stopifnot(v[k,"full_rodaf"]==nrow(e$DESeqResDEGs),all(v[,"full_rodaf"]<=v[,"screening"]))
}
cat("Boundary, sign, NA, no-FC, cached post-rule coverage, and all 17 archive-reference checks passed.\n")
if(file.exists("data/derived/deg_thresholds/threshold_sensitivity.rds")) {
  x<-readRDS("data/derived/deg_thresholds/threshold_sensitivity.rds")
  stopifnot(nrow(x$summary)==2L*4L*8L*12L,all(x$summary$trials==1000L),
            nrow(x$trial_counts)==32000L,ncol(x$trial_counts)==24L,
            all(x$validation$passed))
  m<-x$trial_counts
  for(ep in c("screening","full_rodaf"))for(fc in c("none","1.25","1.5","2.0")) {
    z<-m[,paste0(ep,"__p",c(.01,.025,.05),"_fc",fc)]
    stopifnot(all(z[,1]<=z[,2]),all(z[,2]<=z[,3]))
  }
  for(ep in c("screening","full_rodaf"))for(p in c(.01,.025,.05)) {
    z<-m[,paste0(ep,"__p",p,"_fc",c("none","1.25","1.5","2.0"))]
    stopifnot(all(z[,1]>=z[,2]),all(z[,2]>=z[,3]),all(z[,3]>=z[,4]))
  }
  stopifnot(all(m[,1:12]>=m[,13:24]))
  den<-tapply(x$distribution$trials,interaction(x$distribution$endpoint,x$distribution$cell,x$distribution$n,x$distribution$threshold_id,drop=TRUE),sum)
  stopifnot(all(den==1000L))
  cat("All 32000-trial completeness, endpoint inclusion, distribution and monotonicity checks passed.\n")
}
