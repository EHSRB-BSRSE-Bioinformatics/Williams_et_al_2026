# Exercise the post-processing with two synthetic BMDExpress analysis sets.
source("R/tpod.R")

probes <- paste0("probe_", seq_len(30L))
annotation <- data.frame(Probe_Name = probes, Enzrez_ID = as.character(seq_along(probes)))
pathways <- data.frame(
  collection = rep(tpod_pathway_collections, each = length(probes)),
  pathway = "example_pathway",
  entrez_id = annotation$Enzrez_ID
)
export <- data.frame(
  Analysis = rep(c("Exp 1", "Exp 2"), each = length(probes)),
  `Probe ID` = probes,
  `Best BMD` = rep(c(1, 100), each = length(probes)),
  `Best BMDL` = 0.5,
  `Best BMDU` = 2,
  `Best RSquared` = 0.6,
  `Best fitPValue` = 0.5,
  `Prefilter P-Value` = 0.01,
  `Max Fold Change Absolute Value` = 1.5,
  check.names = FALSE
)
summary <- summarize_tpod_fpr(export, annotation, pathways, expected_trials = 2L)
expected_positive <- summary$Fold_Change_Cutoff == 1.5 & summary$P_Value_Cutoff == 0.05
stopifnot(nrow(summary) == 40L, all(summary$n_total == 2L),
          all(summary$FPR == ifelse(expected_positive, 0.5, 0)))
print(summary, row.names = FALSE)
