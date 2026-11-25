# Consistency audit of the submitted analysis

Run `Rscript scripts/audit_analysis_consistency.R` in the locked default
environment from the repository root. The script reads the original compact
objects from `../2023_Meier_FDR` by default. Pass `--historical-root=PATH`
to use another copy and `--large` to check the original pathway trial lists.

- **DEG.** The submitted p-only panel used the intermediate
  `R-ODAF + P-Value` count under a new label. The archived intermediate
  recycles a logical index across the full feature table. Independently,
  the alternate direct-p summary counts placeholder rows from missing
  adjusted p-values. The audit reproduces both effects in the 17 bundled
  reference trials and reconstructs final calls from feature identities.
- **ORA.** The frozen GO/WikiPathways summary matches the lowercase
  `compiled_summary_boostrapped_ora_FDR_data.wrangled.RDS`. The archived
  workflow selects the largest absolute shrunken fold change per alias,
  maps to Entrez IDs, and applies stored BH p < 0.05 and absolute fold
  change > 1.5. Its background is the mapped result table, intersected
  with the collection by the enrichment software. The quantile and spike
  rules used for final DEGs are separate. The current trial-generation
  helper also collapses by Entrez before screening, which can change the
  selected IDs. The archived ORA lists are unnamed, so a rebuild needs to
  restore the frozen trial IDs before combining them.
- **GSEA.** The original call used `fgseaMultilevel` with
  `nPermSimple=500` for the preliminary simple-permutation stage.
  The checked-in compact archive applies its own 15–500 effective-size
  restriction and recalculates BH from defined raw p-values per trial.
  These historical postprocessing differences are separate from the
  later unique-ID rerun.

Use the fraction of frozen mock trials with at least one call as the
experiment-level false-positive rate. The procedure that adjusts feature or
pathway p-values is a separate part of the analysis. Trial assignments reuse
the same control wells.
