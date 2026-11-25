x<-readRDS("data/derived/tgx_signature/tgx_signature_null.rds")
s<-x$summary;t<-x$trials;g<-x$gene_frequencies
stopifnot(nrow(s)==32*11,all(s$trials==1000),nrow(t)==32000*11,
          all(table(t$trial)==11),nrow(g)==32*64*11,
          all(s$signature_trial_rate<=s$whole_trial_rate | is.na(s$whole_trial_rate)),
          all(s$signature_positive_trials==s$joint_signature_classifier+s$signature_only),
          all(s$classifier_positive_trials==s$joint_signature_classifier+s$classifier_only),
          all(s$trials==s$joint_signature_classifier+s$signature_only+s$classifier_only+s$neither),
          all(t$signature_genes<=64),all(t$signature_genes_strict_probe<=t$signature_genes))
key<-interaction(x$distribution$cell,x$distribution$n,x$distribution$stage,drop=TRUE)
stopifnot(all(tapply(x$distribution$trials,key,sum)==1000))
for(cube in x$membership) {
  stopifnot(all(cube[,,"final"]<=cube[,,"quantile_spike"]),
            all(cube[,,"quantile_spike"]<=cube[,,"quantile"]),
            all(cube[,,"quantile_spike"]<=cube[,,"spike"]),
            all(cube[,,"significant"]<=cube[,,"bh_eligible"]),
            all(cube[,,"bh_eligible"]<=cube[,,"tested"]),
            all(cube[,,"tested"]<=cube[,,"relevance"]),
            all(cube[,,"relevance"]<=cube[,,"initial"]),
            all(cube[,,"initial"]<=cube[,,"detected"]))
}
cat("All 32000 trials: denominator, pairing, distribution and stage checks passed.\n")

frozen<-readRDS("data/derived/monte_carlo_trial_metadata.rds")
source("scripts/tgx_signature/analyze.R")
idx<-parse_ids(names(frozen))
for(key in split(seq_len(nrow(idx)),interaction(idx$cell,idx$n,drop=TRUE))) {
  stopifnot(setequal(idx$iteration[key],1:1000))
  for(i in key)stopifnot(all(table(frozen[[idx$trial[i]]]$group)==idx$n[i]))
}
old<-readRDS("data/derived/deg_monte_carlo_summary.rds")
audit<-list();mismatches<-list()
for(label in c("Final DEGs","P-Value Only","P-Value + Fold Change")) {
  stage<-switch(label,"Final DEGs"="final","P-Value Only"="significant","P-Value + Fold Change"="significant_fc")
  a<-old[!is.na(old$Filter)&old$Filter==label,];b<-t[t$stage==stage,]
  check_ids(a$name,b$trial,label)
  counts<-b$whole_features[match(a$name,b$trial)]
  delta<-as.integer(a$nDEGs)-as.integer(counts)
  bad<-which(delta!=0L)
  audit[[label]]<-data.frame(publication_label=label,archive_stage=stage,trials=nrow(a),discrepant_count_trials=length(bad))
  if(length(bad))mismatches[[label]]<-data.frame(trial=a$name[bad],publication_label=label,
    archive_stage=stage,publication_count=a$nDEGs[bad],archive_rule_count=counts[bad])
  if(label=="Final DEGs")stopifnot(length(bad)==0L)
}
write_tsv(do.call(rbind,audit),"data/derived/tgx_signature/publication_summary_audit.tsv")
con<-gzfile("data/derived/tgx_signature/publication_summary_mismatches.tsv.gz","wt")
write_tsv(do.call(rbind,mismatches),con);close(con)
cat("Frozen iteration/group completeness and all publication final DEG counts match.\n")
print(do.call(rbind,audit),row.names=FALSE)
cat("Intermediate publication-label differences are diagnostic; archived feature tables define the new stage counts.\n")
