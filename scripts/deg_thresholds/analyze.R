#!/usr/bin/env Rscript
# Post hoc thresholds on fixed archived results; no model fits or BH adjustment.
check <- function(ok, message) if (!isTRUE(ok)) stop(message, call.=FALSE)
write_tsv <- function(x, path) write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
check_ids <- function(x, expected, label) {
  check(!anyNA(x) && !anyDuplicated(x) && setequal(x,expected),paste(label,"trial IDs differ or are duplicated"))
}
trial_index <- function(ids) {
  parts<-strsplit(ids,"_",fixed=TRUE)
  check(all(lengths(parts)==4L),"Malformed trial IDs")
  z<-do.call(rbind,parts)
  check(all(z[,2]==z[,3]),"Unequal group sizes")
  data.frame(trial=ids,cell=z[,1],n=as.integer(z[,2]),iteration=as.integer(z[,4]))
}
threshold_grid <- function() {
  g<-expand.grid(padj_cutoff=c(.01,.025,.05),linear_fc=c(NA_real_,1.25,1.5,2),
                 KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  g$fc_label<-ifelse(is.na(g$linear_fc),"none",format(g$linear_fc,trim=TRUE))
  # Numeric labels are fixed explicitly to avoid platform-dependent formatting.
  g$fc_label<-rep(c("none","1.25","1.5","2.0"),each=3)
  g$threshold_id<-paste0("p",g$padj_cutoff,"_fc",g$fc_label)
  g
}
select_threshold <- function(padj,lfc,p,fc=NA_real_) {
  check(p<=.05 && p>0,"Unsupported p cutoff outside archived post-rule coverage")
  keep<-!is.na(padj)&padj<p
  if(!is.na(fc)) {
    check(fc>1,"Use NA for no fold-change cutoff")
    keep<-keep & !is.na(lfc) & abs(lfc)>log2(fc)
  }
  which(keep)
}
extract_envelope <- function(entry) {
  r<-as.data.frame(entry$DESeqRes)
  check(all(c("padj","log2FoldChange") %in% names(r)),"Missing archived columns")
  check(!anyDuplicated(rownames(r)),"Duplicate result feature IDs")
  check(all(is.na(r$padj) | (r$padj>=0 & r$padj<=1)),"Invalid stored adjusted p-value")
  sig<-select_threshold(r$padj,r$log2FoldChange,.05)
  s<-r[sig,,drop=FALSE]
  qbad<-sbad<-character()
  if(nrow(s)) {
    check(length(entry$DESeqjs_filtered_out)==3L,"Missing post-rule exclusions")
    qbad<-rownames(entry$DESeqjs_filtered_out[[2]])
    sbad<-rownames(entry$DESeqjs_filtered_out[[3]])
    check(!anyDuplicated(qbad) && !anyDuplicated(sbad) && all(c(qbad,sbad) %in% rownames(s)),"Post-rule IDs outside significant envelope")
    check(entry$DESeqfiltered_table$quantile_filtered==length(qbad) &&
          entry$DESeqfiltered_table$spike_filtered==length(sbad),"Exclusion tally differs")
  }
  s$quantile_pass<-!rownames(s) %in% qbad
  s$spike_pass<-!rownames(s) %in% sbad
  original<-select_threshold(s$padj,s$log2FoldChange,.05,1.5)
  final<-rownames(s)[original[s$quantile_pass[original]&s$spike_pass[original]]]
  check(setequal(final,rownames(entry$DESeqResDEGs)),"Original final feature set does not match archive")
  list(features=data.frame(feature=rownames(s),padj=s$padj,log2FoldChange=s$log2FoldChange,
                          quantile_pass=s$quantile_pass,spike_pass=s$spike_pass),
       qc=c(relevant=nrow(r),missing_padj=sum(is.na(r$padj)),significant_005=nrow(s),
            quantile_rejections=length(qbad),spike_rejections=length(sbad),archived_final=length(final)))
}
count_grid <- function(features,grid=threshold_grid()) {
  screen<-full<-integer(nrow(grid))
  for(k in seq_len(nrow(grid))) {
    i<-select_threshold(features$padj,features$log2FoldChange,grid$padj_cutoff[k],grid$linear_fc[k])
    screen[k]<-length(i)
    full[k]<-sum(features$quantile_pass[i]&features$spike_pass[i])
  }
  cbind(screening=screen,full_rodaf=full)
}
sha256 <- function(paths) {
  check(nzchar(Sys.which("sha256sum")),"sha256sum is required")
  vapply(paths,function(p)strsplit(system2("sha256sum",shQuote(p),stdout=TRUE)," +")[[1]][1],"")
}
assemble <- function(out,frozen,idx,grid,inputs) {
  stems<-c(TK6="tk6",`MCF-7`="mcf7",HepaRG="heparg",`Liver spheroids`="spheroids")
  parts<-lapply(stems,function(stem)readRDS(file.path(out,paste0("counts_",stem,".rds"))))
  ids<-unlist(lapply(parts,function(x)x$index$trial),use.names=FALSE)
  check_ids(ids,names(frozen),"Combined counts")
  trial_counts<-do.call(rbind,lapply(parts,`[[`,"counts"));rownames(trial_counts)<-ids
  qc<-do.call(rbind,lapply(parts,`[[`,"qc"));qc$trial<-ids
  summary<-list();distributions<-list()
  for(endpoint in c("screening","full_rodaf"))for(cell in names(stems))for(n in 3:10)for(k in seq_len(nrow(grid))) {
    selected<-which(idx$cell==cell&idx$n==n)
    values<-trial_counts[match(idx$trial[selected],ids),paste(endpoint,grid$threshold_id[k],sep="__")]
    check(length(values)==1000 && !anyNA(values),"Incomplete stratum")
    key<-paste(endpoint,cell,n,k)
    summary[[key]]<-cbind(data.frame(endpoint=endpoint,cell=cell,n=n),grid[k,,drop=FALSE],
      data.frame(positive_trials=sum(values>0),trials=length(values),rate=mean(values>0),
                 total_degs=sum(values),median_degs=median(values),max_degs=max(values)))
    tab<-table(values)
    distributions[[key]]<-data.frame(endpoint=endpoint,cell=cell,n=n,threshold_id=grid$threshold_id[k],
                                    n_degs=as.integer(names(tab)),trials=as.integer(tab))
  }
  summary<-do.call(rbind,summary);distribution<-do.call(rbind,distributions)
  # Compare literal compact labels and the actual source selected by Figure 2.
  old<-readRDS("data/derived/deg_monte_carlo_summary.rds")
  audit<-list();mismatch<-list()
  comparisons<-data.frame(label=c("Final DEGs","P-Value Only","P-Value + Fold Change","R-ODAF + P-Value"),
     endpoint=c("full_rodaf","screening","screening","screening"),fc=c("1.5","none","1.5","none"))
  for(i in seq_len(nrow(comparisons))) {
    spec<-comparisons[i,];a<-old[!is.na(old$Filter)&old$Filter==spec$label,]
    check_ids(a$name,names(frozen),spec$label)
    col<-paste0(spec$endpoint,"__p0.05_fc",spec$fc)
    actual<-trial_counts[match(a$name,ids),col]
    bad<-which(a$nDEGs!=actual)
    if(spec$label=="Final DEGs")check(!length(bad),"Compact published final counts differ")
    for(cell in names(stems))for(n in 3:10) {
      sel<-which(as.character(a$cell_line)==cell&a$n_per_group==n)
      audit[[paste(i,cell,n)]]<-data.frame(publication_label=spec$label,endpoint=spec$endpoint,fc_label=spec$fc,
        cell=cell,n=n,trials=length(sel),discrepant_count_trials=sum(a$nDEGs[sel]!=actual[sel]),
        publication_positive_trials=sum(a$nDEGs[sel]>0),archive_positive_trials=sum(actual[sel]>0))
    }
    if(length(bad))mismatch[[as.character(i)]]<-data.frame(trial=a$name[bad],publication_label=spec$label,
      endpoint=spec$endpoint,fc_label=spec$fc,publication_count=a$nDEGs[bad],archive_count=actual[bad])
  }
  audit<-do.call(rbind,audit);mismatch<-do.call(rbind,mismatch)
  validation<-data.frame(check=c("trial_IDs_and_assignments","original_final_feature_membership","original_final_compact_counts"),
                         trials=32000L,passed=TRUE)
  result<-list(summary=summary,trial_index=idx[match(ids,idx$trial),],trial_counts=trial_counts,
    grid=grid,distribution=distribution,trial_qc=qc,publication_comparison=audit,
    validation=validation,provenance=inputs,session=public_session_info(),
    estimand="Proportion of frozen trials with >=1 retained feature; stored BH and shrinkage fixed",
    post_rule_basis="Feature-aligned archived exclusions covering all padj<0.05 features, before FC filtering")
  saveRDS(result,file.path(out,"threshold_sensitivity.rds"),compress="xz")
  write_tsv(summary,file.path(out,"threshold_rates.tsv"))
  write_tsv(distribution,file.path(out,"deg_count_distributions.tsv"))
  write_tsv(qc,file.path(out,"trial_qc.tsv"))
  write_tsv(audit,file.path(out,"publication_cutoff_comparison.tsv"))
  con<-gzfile(file.path(out,"publication_count_mismatches.tsv.gz"),"wt");write_tsv(mismatch,con);close(con)
  con<-gzfile(file.path(out,"trial_counts.tsv.gz"),"wt")
  write_tsv(cbind(result$trial_index,as.data.frame(trial_counts)),con);close(con)
  write_tsv(inputs,file.path(out,"provenance.tsv"));write_tsv(validation,file.path(out,"validation.tsv"))
  writeLines(result$session,file.path(out,"session_info.txt"))
  message("Completed: ",out)
}
main <- function(args=commandArgs(TRUE)) {
  arg<-function(key,default){i<-match(key,args);if(is.na(i))default else args[i+1L]}
  source("R/provenance.R")
  parent<-normalizePath(arg("--source","../2023_Meier_FDR"))
  out<-arg("--output","data/derived/deg_thresholds")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  suppressPackageStartupMessages(library(DESeq2)) # S4 deserialization only
  frozen<-readRDS("data/derived/monte_carlo_trial_metadata.rds")
  check_ids(names(frozen),unique(names(frozen)),"Frozen metadata")
  idx<-trial_index(names(frozen));check(nrow(idx)==32000L,"Expected 32000 trials")
  for(s in split(seq_len(nrow(idx)),interaction(idx$cell,idx$n,drop=TRUE))) {
    check(length(s)==1000L && setequal(idx$iteration[s],1:1000),"Incomplete frozen stratum")
    for(i in s)check(identical(names(table(as.character(frozen[[idx$trial[i]]]$group))),c("A","B")) &&
                     all(table(frozen[[idx$trial[i]]]$group)==idx$n[i]),"Frozen group sizes differ from ID")
  }
  counts<-as.matrix(read.delim("data/input/count_table.tsv",row.names=1,check.names=FALSE))
  check(!anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)),"Duplicate input identifiers")
  totals<-colSums(counts);check(all(is.finite(totals)) && all(counts>=0),"Invalid input counts")
  grid<-threshold_grid()
  stems<-c(TK6="tk6",`MCF-7`="mcf7",HepaRG="heparg",`Liver spheroids`="spheroids")
  sources<-c("data/derived/monte_carlo_trial_metadata.rds","data/input/count_table.tsv",
    "data/derived/deg_monte_carlo_summary.rds","data/reference/deg_source_archives.tsv",
    "docs/deseq2_methods.md","R/deg_trials.R","R/summaries.R","R/figures.R",
    "scripts/deg_thresholds/analyze.R","R/provenance.R","tests/test_deg_thresholds.R")
  archives<-file.path(parent,"bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10",paste0(stems,".RDS"))
  manifest<-read.delim("data/reference/deg_source_archives.tsv")
  input_fingerprint<-sha256(sources)
  for(cell in names(stems)) {
    stem<-stems[[cell]];path<-archives[match(cell,names(stems))]
    cache<-file.path(out,paste0("counts_",stem,".rds"))
    if("--resume" %in% args && file.exists(cache)) {
      cached<-readRDS(cache)
      check(identical(cached$input_fingerprint,input_fingerprint) && identical(cached$grid,grid),"Cached input fingerprints or grid differ")
      check(unname(sha256(path))==cached$archive_sha256,"Cached archive fingerprint differs")
      message("Reusing checked cache for ",cell);next
    }
    message("Reading ",cell," archive")
    archive<-readRDS(path);expected<-idx$trial[idx$cell==cell]
    check_ids(names(archive),expected,paste(cell,"archive"))
    values<-matrix(0L,length(expected),nrow(grid)*2L,
       dimnames=list(expected,c(paste0("screening__",grid$threshold_id),paste0("full_rodaf__",grid$threshold_id))))
    qc<-vector("list",length(expected))
    for(j in seq_along(expected)) {
      id<-expected[j];e<-archive[[id]];f<-frozen[[id]]
      check(identical(as.character(e$exp_metadata$sample_id),as.character(f$sample_id)) &&
            identical(as.character(e$exp_metadata$group),as.character(f$group)),paste(id,"archive assignments differ"))
      check(!anyDuplicated(f$sample_id) && all(f$sample_id %in% colnames(counts)) &&
            all(totals[f$sample_id]>100000),paste(id,"invalid or low-depth samples"))
      check(all(rownames(e$DESeqRes) %in% rownames(counts)),paste(id,"unknown feature IDs"))
      envelope<-extract_envelope(e)
      values[j,]<-as.vector(count_grid(envelope$features,grid))
      qc[[j]]<-envelope$qc
      if(j%%1000L==0L)message(cell," ",j,"/",length(expected))
    }
    rm(archive);gc()
    message("Verifying ",cell," archive checksum")
    hash<-unname(sha256(path));check(hash==manifest$sha256[match(basename(path),manifest$archive_file)],"Archive fingerprint differs")
    saveRDS(list(index=idx[match(expected,idx$trial),],counts=values,qc=as.data.frame(do.call(rbind,qc)),
                 archive_path=archive_relative_path(path,parent),archive_sha256=hash,grid=grid,input_fingerprint=input_fingerprint),cache,compress="xz")
    rm(values,qc);gc()
  }
  message("Fingerprinting inputs")
  check(identical(input_fingerprint,sha256(sources)),"Inputs changed during extraction")
  inputs<-data.frame(path=sources,sha256=unname(input_fingerprint))
  for(cell in names(stems)) {
    cache<-readRDS(file.path(out,paste0("counts_",stems[[cell]],".rds")))
    check(identical(cache$grid,grid),"Cached grid differs")
    inputs<-rbind(inputs,data.frame(path=archive_relative_path(cache$archive_path,parent),sha256=cache$archive_sha256))
  }
  assemble(out,frozen,idx,grid,inputs)
}
if(sys.nframe()==0L)main()
