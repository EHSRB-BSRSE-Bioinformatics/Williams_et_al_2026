#!/usr/bin/env Rscript
# Extract signature-specific null behavior from frozen assignments and archives.
# No model fitting, biomarker training, or subset-specific p-value adjustment.
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
check_ids <- function(observed, expected, label) {
  assert(!anyNA(observed) && !anyDuplicated(observed), paste(label,"missing/duplicate IDs"))
  assert(setequal(observed, expected), paste(label,"ID set differs from frozen assignments"))
}
parse_ids <- function(ids) {
  m <- regexec("^(.*)_([0-9]+)_([0-9]+)_([0-9]+)$", ids)
  z <- regmatches(ids,m)
  assert(all(lengths(z)==5L),"Invalid trial ID")
  z <- do.call(rbind,z)
  assert(all(z[,3]==z[,4]),"Unequal group sizes")
  data.frame(trial=ids,cell=z[,2],n=as.integer(z[,3]),iteration=as.integer(z[,5]))
}
select_stages <- function(entry) {
  r <- as.data.frame(entry$DESeqRes)
  assert(!anyDuplicated(rownames(r)),"Duplicate archived features")
  sig <- rownames(r)[!is.na(r$padj) & r$padj < .05]
  fc <- rownames(r)[!is.na(r$padj) & r$padj < .05 & !is.na(r$log2FoldChange) & abs(r$log2FoldChange)>log2(1.5)]
  qbad <- sbad <- character()
  if(length(sig)) {
    assert(length(entry$DESeqjs_filtered_out)==3L,"Missing exclusion tables")
    qbad <- rownames(entry$DESeqjs_filtered_out[[2]])
    sbad <- rownames(entry$DESeqjs_filtered_out[[3]])
    assert(all(c(qbad,sbad) %in% sig),"Exclusions are not significant features")
  }
  q <- setdiff(sig,qbad); s <- setdiff(sig,sbad); both <- intersect(q,s)
  final <- intersect(both,fc)
  assert(setequal(final,rownames(entry$DESeqResDEGs)),"Final DEG archive disagrees with reconstructed rules")
  list(relevance=rownames(r),tested=rownames(r)[!is.na(r$pvalue)],
       bh_eligible=rownames(r)[!is.na(r$padj)],significant=sig,
       significant_fc=fc,quantile=q,spike=s,quantile_spike=both,final=final)
}
gene_members <- function(features, probe_sets) {
  # Restrict the large tested-feature vector once, before checking each gene.
  hits<-intersect(features,unlist(probe_sets,use.names=FALSE))
  vapply(probe_sets,function(p) any(p %in% hits),logical(1))
}
main <- function(args=commandArgs(TRUE)) {
  arg <- function(key,default) {i<-match(key,args);if(is.na(i))default else args[i+1L]}
  root <- normalizePath(arg("--root","."));setwd(root)
  parent <- normalizePath(arg("--source","../2023_Meier_FDR"))
  out <- arg("--output","data/derived/tgx_signature")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  suppressPackageStartupMessages(library(DESeq2)) # deserialize historical S4 tables only
  refpath <- "data/reference/tgx_ddi_64_gene_map.tsv"
  ref <- read.delim(refpath,check.names=FALSE,na.strings="NA")
  assert(nrow(ref)==64L && !anyDuplicated(ref$published_symbol) && !anyDuplicated(ref$measured_symbol),"Signature must contain 64 unique genes")
  counts <- as.matrix(read.delim("data/input/count_table.tsv",row.names=1,check.names=FALSE))
  assert(all(is.finite(counts)) && all(counts>=0),"Invalid count values")
  assert(!anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)),"Duplicate count IDs")
  symbols <- sub("_[^_]*$","",rownames(counts))
  probes <- lapply(ref$measured_symbol,function(s)rownames(counts)[symbols==s])
  names(probes)<-ref$published_symbol
  ref$n_measured_probes <- lengths(probes)
  ref$measured_probes <- vapply(probes,paste,collapse=";",FUN.VALUE="")
  ref$multiple_probes <- lengths(probes)>1L
  # A multi-probe gene is retained once if any of its unchanged probe tests passes.
  allprobes <- unlist(probes,use.names=FALSE)
  assert(!anyDuplicated(allprobes),"Probe maps to multiple signature genes")
  ref$mapping_status <- vapply(seq_len(64),function(i) {
    p<-probes[[i]];r<-if(is.na(ref$reference_probes[i]))character() else strsplit(ref$reference_probes[i],";",fixed=TRUE)[[1]]
    if(!length(p))"unmeasured" else if(setequal(p,r))"exact_probe_reference" else if(!length(r))"symbol_only_no_reference_probe" else "probe_reference_mismatch"
  },"")
  assert(!any(ref$mapping_status=="probe_reference_mismatch"),"Unresolved probe-reference mismatch")
  write_tsv(ref,file.path(out,"gene_map.tsv"))
  frozen <- readRDS("data/derived/monte_carlo_trial_metadata.rds")
  idx <- parse_ids(names(frozen));check_ids(names(frozen),unique(names(frozen)),"Frozen")
  assert(length(frozen)==32000L,"Expected 32000 frozen trials")
  for(key in split(seq_len(nrow(idx)),interaction(idx$cell,idx$n,drop=TRUE))) {
    assert(setequal(idx$iteration[key],1:1000),"Frozen iteration IDs incomplete")
    for(i in key) {
      groups<-table(as.character(frozen[[idx$trial[i]]]$group))
      assert(identical(names(groups),c("A","B")) && all(groups==idx$n[i]),"Frozen n/group count mismatch")
    }
  }
  den <- as.data.frame(table(idx$cell,idx$n));assert(nrow(den)==32L && all(den$Freq==1000L),"Incomplete cell/n denominators")
  paths <- c(refpath,"data/input/count_table.tsv","data/derived/monte_carlo_trial_metadata.rds","data/derived/tgx_ddi_monte_carlo_summary.rds","data/reference/tgx_ddi_mapping_sources.tsv")
  rows <- list(); genestage <- list(); sparse <- list(); qc <- list(); audit <- list()
  stems<-c(TK6="tk6",`MCF-7`="mcf7",HepaRG="heparg",`Liver spheroids`="spheroids")
  if("--resume" %in% args && file.exists(file.path(out,"extraction_checkpoint.rds"))) {
    checkpoint<-readRDS(file.path(out,"extraction_checkpoint.rds"))
    rows<-checkpoint$rows;sparse<-checkpoint$stages;genestage<-checkpoint$genes
    qc<-checkpoint$qc;audit<-checkpoint$calls;paths<-checkpoint$paths
    message("Resuming completed cells: ",paste(names(rows),collapse=", "))
  }
  stage_names <- c("detected","initial","relevance","tested","bh_eligible","significant","significant_fc","quantile","spike","quantile_spike","final")
  selected_cells<-if("--cells" %in% args)strsplit(arg("--cells",""),",",fixed=TRUE)[[1]] else names(stems)
  assert(all(selected_cells %in% names(stems)),"Unknown selected cell")
  for(cell in setdiff(selected_cells,names(rows))) {
    path <- file.path(parent,"bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10",paste0(stems[[cell]],".RDS"))
    message("Reading ",cell,": ",path)
    archive<-readRDS(path); paths<-c(paths,path)
    expected<-idx$trial[idx$cell==cell];check_ids(names(archive),expected,paste(cell,"DEG archive"))
    cube<-array(FALSE,c(length(expected),64L,length(stage_names)),dimnames=list(expected,ref$published_symbol,stage_names))
    cellrows<-vector("list",length(expected))
    for(j in seq_along(expected)) {
      id<-expected[j]; e<-archive[[id]]; f<-frozen[[id]]
      assert(identical(as.character(e$exp_metadata$sample_id),as.character(f$sample_id)) && identical(as.character(e$exp_metadata$group),as.character(f$group)),paste(id,"sample assignments differ"))
      assert(!anyDuplicated(f$sample_id) && all(f$sample_id %in% colnames(counts)),paste(id,"sample IDs"))
      assert(all(colSums(counts[,f$sample_id,drop=FALSE])>100000),paste(id,"sample fails depth cutoff"))
      if(!is.null(e$trial_name)) assert(identical(e$trial_name,id),paste(id,"embedded trial ID differs"))
      stages<-select_stages(e)
      assert(all(stages$relevance %in% rownames(counts)),paste(id,"unmapped archive features"))
      subset<-counts[allprobes,f$sample_id,drop=FALSE]
      initial<-rownames(subset)[rowSums(round(subset))>1]
      stages<-c(list(detected=rownames(subset)[rowSums(subset)>0],initial=initial),stages)
      for(k in seq_along(stages))cube[j,,k]<-gene_members(stages[[k]],probes)
      assert(all(stages$relevance[stages$relevance %in% allprobes] %in% initial),paste(id,"relevance/initial mismatch"))
      if(length(stages$significant)) {
        assert(e$DESeqfiltered_table$passed_all_filters==length(stages$final),"Final filter tally mismatch")
      }
      cellrows[[j]]<-data.frame(trial=id,cell=cell,n=idx$n[match(id,idx$trial)],stage=stage_names,
        signature_genes=as.integer(colSums(cube[j,,,drop=FALSE],dims=2)),
        signature_genes_strict_probe=as.integer(colSums(cube[j,ref$mapping_status=="exact_probe_reference",,drop=FALSE],dims=2)),
        signature_features=vapply(stages,function(z)sum(allprobes %in% z),integer(1)),
        whole_features=c(NA_integer_,NA_integer_,vapply(stages[-c(1,2)],length,integer(1))))
      if(j%%1000L==0L)message(cell," ",j,"/",length(expected))
    }
    cr<-do.call(rbind,cellrows);rows[[cell]]<-cr;sparse[[cell]]<-cube
    for(n in 3:10) {
      sel<-idx$n[match(expected,idx$trial)]==n
      for(k in seq_along(stage_names)) {
        mat<-matrix(cube[sel,,k],ncol=64)
        genestage[[paste(cell,n,k)]]<-data.frame(cell=cell,n=n,published_symbol=ref$published_symbol,stage=stage_names[k],trials_retained=colSums(mat),trials= sum(sel),proportion_trials=colMeans(mat))
      }
    }
    ids<-unique(unlist(lapply(frozen[expected],function(z)z$sample_id)))
    refcell<-ref;refcell$cell<-cell
    refcell$positive_count_samples<-vapply(probes,function(p)if(length(p))sum(colSums(counts[p,ids,drop=FALSE])>0)else 0L,integer(1))
    refcell$total_samples<-length(ids);qc[[cell]]<-refcell
    # Classifier records contain their own IDs; list names are not trusted.
    cp<-file.path(parent,"bootstrap_TGx-DDI/output/1000_iterations/nobatch_2023-09-01",paste0(stems[[cell]],"_tgx_ddi_results.RDS"))
    calls<-readRDS(cp);paths<-c(paths,cp)
    callids<-vapply(calls,function(z)as.character(z$trial_res),"")
    check_ids(callids,expected,paste(cell,"classifier archive"))
    audit[[cell]]<-do.call(rbind,lapply(calls,function(z)data.frame(trial=as.character(z$trial_res),probG=as.numeric(z$probG_res),pc1=as.numeric(z$pc1_res),hca=if(length(z$hca))as.character(z$hca) else NA_character_)))
    rm(archive,cube,cellrows);gc()
    saveRDS(list(rows=rows,stages=sparse,genes=genestage,qc=qc,calls=audit,paths=paths),file.path(out,"extraction_checkpoint.rds"),compress=FALSE)
  }
  if("--extract-only" %in% args)return(invisible(NULL))
  trials<-do.call(rbind,rows);calls<-do.call(rbind,audit);check_ids(calls$trial,idx$trial,"Combined classifier")
  classifier_input_rows<-list()
  for(cell in names(stems)) {
    bp<-file.path(parent,"bootstrap_TGx-DDI/output/1000_iterations/nobatch_2023-09-01",paste0(stems[[cell]],"-biomarker.RDS"))
    input<-readRDS(bp);paths<-c(paths,bp)
    expected<-idx$trial[idx$cell==cell];check_ids(names(input),expected,paste(cell,"classifier input"))
    for(id in expected) {
      e<-input[[id]];f<-frozen[[id]]
      assert(identical(as.character(e$exp_metadata$sample_id),as.character(f$sample_id)) &&
             identical(as.character(e$exp_metadata$group),as.character(f$group)),paste(id,"classifier input assignments differ"))
      genes<-as.character(e[[1]][[1]])
      assert(!anyDuplicated(genes) && setequal(genes,ref$measured_symbol),paste(id,"classifier input gene set differs"))
    }
    classifier_input_rows[[cell]]<-data.frame(cell=cell,trials=length(expected),assignment_matches=length(expected),classifier_genes=64L)
    rm(input);gc()
  }
  write_tsv(do.call(rbind,classifier_input_rows),file.path(out,"classifier_input_qc.tsv"))
  compact<-readRDS("data/derived/tgx_ddi_monte_carlo_summary.rds")
  hc<-unique(compact[,c("trial_res","hca")]);check_ids(hc$trial_res,idx$trial,"Frozen HCA summary")
  saved_hca<-as.character(hc$hca[match(calls$trial,hc$trial_res)])
  calls$raw_hca_available<-!is.na(calls$hca)
  assert(all(is.na(calls$hca)|calls$hca==saved_hca),"Raw/summary HCA mismatch")
  calls$hca<-saved_hca
  calls$hca_source<-ifelse(calls$raw_hca_available,"raw_and_frozen_summary","frozen_summary_only")
  assert(all(is.finite(calls$probG)) && all(calls$probG>=0 & calls$probG<=1) && all(is.finite(calls$pc1)) && !anyNA(saved_hca),"Incomplete classifier values")
  calls$positive_probability<-calls$probG>.9;calls$positive_pc1<-calls$pc1<0;calls$positive_hca<-calls$hca=="DDI"
  calls$positive_any_component<-calls$positive_probability|calls$positive_pc1|calls$positive_hca
  for(test in c("probG_res","pc1_res")) {
    csub<-compact[compact$tgx_ddi_test==test,];check_ids(csub$trial_res,idx$trial,paste("Classifier summary",test))
    target<-if(test=="probG_res")calls$probG else calls$pc1
    assert(isTRUE(all.equal(csub$tgx_ddi_value,target[match(csub$trial_res,calls$trial)],check.attributes=FALSE)),"Classifier summary numeric mismatch")
  }
  trials<-merge(trials,calls,by="trial",sort=FALSE,all.x=TRUE)
  assert(!anyNA(trials$positive_any_component),"Missing classifier outcomes")
  summary<-do.call(rbind,lapply(split(trials,interaction(trials$cell,trials$n,trials$stage,drop=TRUE)),function(d) {
    data.frame(cell=d$cell[1],n=d$n[1],stage=d$stage[1],trials=nrow(d),
      signature_positive_trials=sum(d$signature_genes>0),signature_trial_rate=mean(d$signature_genes>0),
      strict_probe_trial_rate=mean(d$signature_genes_strict_probe>0),
      whole_positive_trials=if(all(is.na(d$whole_features)))NA_integer_ else sum(d$whole_features>0),
      whole_trial_rate=if(all(is.na(d$whole_features)))NA_real_ else mean(d$whole_features>0),
      genes_min=min(d$signature_genes),genes_median=median(d$signature_genes),genes_mean=mean(d$signature_genes),genes_max=max(d$signature_genes),
      classifier_positive_trials=sum(d$positive_any_component),
      joint_signature_classifier=sum(d$signature_genes>0 & d$positive_any_component),
      signature_only=sum(d$signature_genes>0 & !d$positive_any_component),
      classifier_only=sum(d$signature_genes==0 & d$positive_any_component),
      neither=sum(d$signature_genes==0 & !d$positive_any_component))
  }))
  assert(all(summary$trials==1000L),"Output denominator incomplete")
  distribution<-as.data.frame(table(trials$cell,trials$n,trials$stage,factor(trials$signature_genes,levels=0:64)))
  names(distribution)<-c("cell","n","stage","signature_genes","trials")
  distribution<-distribution[distribution$trials>0,];distribution$proportion_trials<-distribution$trials/1000
  genes<-do.call(rbind,genestage);mapqc<-do.call(rbind,qc)
  paths<-unique(c(paths,"scripts/tgx_signature/analyze.R","tests/test_tgx_signature.R",
                  "tests/check_tgx_signature_outputs.R","R/provenance.R","data/reference/deg_source_archives.tsv",
                  "data/derived/deg_monte_carlo_summary.rds"))
  assert(nzchar(Sys.which("sha256sum")),"sha256sum is required for input verification")
  message("Fingerprinting inputs and checking archive SHA-256 values")
  hashes<-vapply(paths,function(p)strsplit(system2("sha256sum",shQuote(p),stdout=TRUE)," +")[[1]][1],"")
  source("R/provenance.R")
  provenance<-data.frame(path=archive_relative_path(paths,parent),sha256=unname(hashes))
  manifest<-read.delim("data/reference/deg_source_archives.tsv")
  for(i in seq_len(nrow(manifest))) {
    hit<-which(basename(paths)==manifest$archive_file[i])
    assert(length(hit)==1L && hashes[hit]==manifest$sha256[i],"Archive SHA-256 does not match frozen source manifest")
  }
  result<-list(summary=summary,trials=trials,distribution=distribution,gene_frequencies=genes,gene_map=mapqc,
               membership=sparse,provenance=provenance,session=public_session_info())
  saveRDS(result,file.path(out,"tgx_signature_null.rds"),compress="xz")
  for(name in c("summary","trials","distribution","gene_frequencies","gene_map","provenance")) {
    value<-switch(name,gene_frequencies=genes,gene_map=mapqc,get(name))
    write_tsv(value,file.path(out,paste0(name,".tsv")))
  }
  final<-summary[summary$stage=="final",];final<-final[order(final$cell,final$n),]
  write_tsv(final,file.path(out,"supplementary_table.tsv"))
  pdf(file.path(out,"supplementary_figure.pdf"),width=9,height=7)
  par(mfrow=c(2,2),mar=c(4,4,3,1))
  for(cell in names(stems)) {
    d<-final[final$cell==cell,]
    plot(d$n,100*d$whole_trial_rate,type="b",pch=16,col="#606060",ylim=c(0,max(1,100*d$whole_trial_rate)),xlab="Samples per group",ylab="Trials with >=1 final DEG (%)",main=cell)
    lines(d$n,100*d$signature_trial_rate,type="b",pch=17,col="#007C91")
    legend("topright",c("Whole transcriptome","TGx-DDI genes"),col=c("#606060","#007C91"),pch=c(16,17),lty=1,bty="n",cex=.8)
  };dev.off()
  writeLines(result$session,file.path(out,"session_info.txt"))
  unlink(file.path(out,"extraction_checkpoint.rds"))
  finish_outputs(out)
  message("Completed: ",out)
}
finish_outputs <- function(out) {
  x<-readRDS(file.path(out,"tgx_signature_null.rds"))
  detail<-merge(x$gene_frequencies,x$gene_map,by=c("cell","published_symbol"),sort=FALSE)
  write_tsv(detail,file.path(out,"gene_map_filter_qc.tsv"))
  t<-x$trials
  pair<-do.call(rbind,lapply(c("positive_probability","positive_pc1","positive_hca","positive_any_component"),function(component) {
    do.call(rbind,lapply(split(t,interaction(t$cell,t$n,t$stage,drop=TRUE)),function(d) {
      p<-d[[component]];a<-d$signature_genes>0
      data.frame(cell=d$cell[1],n=d$n[1],stage=d$stage[1],classifier_component=component,
        trials=nrow(d),both=sum(a&p),signature_only=sum(a&!p),classifier_only=sum(!a&p),neither=sum(!a&!p))
    }))
  }))
  write_tsv(pair,file.path(out,"classifier_pairings.tsv"))
  # Compact trial data retain every stage count, with one classifier record per trial.
  keys<-unique(t[,c("trial","cell","n","probG","pc1","hca","positive_probability","positive_pc1","positive_hca","positive_any_component","raw_hca_available","hca_source")])
  for(stage in unique(t$stage)) {
    z<-t[t$stage==stage,];i<-match(keys$trial,z$trial)
    keys[[paste0("signature_",stage)]]<-z$signature_genes[i]
    keys[[paste0("signature_strict_probe_",stage)]]<-z$signature_genes_strict_probe[i]
    keys[[paste0("signature_features_",stage)]]<-z$signature_features[i]
    keys[[paste0("whole_",stage)]]<-z$whole_features[i]
  }
  con<-gzfile(file.path(out,"trials_wide.tsv.gz"),"wt");write_tsv(keys,con);close(con)
  unlink(file.path(out,"trials.tsv"))
}
if(sys.nframe()==0L)main()
