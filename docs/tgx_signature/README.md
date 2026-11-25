# TGx-DDI signature-specific null behavior

Run from the repository root in the locked study environment:

```bash
pixi run --locked -e deg-study Rscript tests/test_tgx_signature.R
pixi run --locked -e deg-study Rscript scripts/tgx_signature/analyze.R \
  --source ../2023_Meier_FDR --output data/derived/tgx_signature
pixi run --locked -e deg-study Rscript tests/check_tgx_signature_outputs.R
pixi run --locked -e deg-study Rscript scripts/tgx_signature/render_figure.R
```

The extraction requires the four original DEG archives (19 GB compressed in total),
the four small classifier-call archives, the frozen trial metadata, and the count
table. It reads one DEG archive at a time and needs substantial RAM. It does not
fit DESeq2 models, train the classifier, change frozen assignments, or adjust
p-values within the signature subset. Results are written only to the new output
directory. Input checksums and the execution session are recorded with results.

## Signature identity and mapping

The fixed panel is the **64 distinct genes**, rather than the original 65
microarray probes of TGx-28.65. The published list is Appendix I of the
[FDA letter of support, 24 October 2017](https://www.fda.gov/files/drugs/published/LOS-HESI-TGx-DDI.pdf).
The underlying classifier derives from Li et al. (2015); the later annotation
consolidation to 64 genes is described by
[Buick et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6588084/).
The published symbols were checked against the 64 `ID` values in the retained
classifier file. No genes were selected using their null-trial behavior.

`data/reference/tgx_ddi_64_gene_map.tsv` preserves all 64 published symbols,
updated measured symbols, Entrez and Ensembl identifiers, and the explicit
TempO-Seq probe references. It was transcribed from these source files.

- `Biomarker_TGx-DDI_Li_2015_Tk6_cells.txt`: original symbol, Ensembl ID,
  Entrez ID, and fixed centroid scores (64 rows).
- `Biomarker TGx-DDI Li 2015 Tk6 cells revised Gene Symbol for TempoSeq v2 pc comps and data for heatmap log10.txt`:
  original-to-measured symbol crosswalk and fixed classifier (64 rows).
- `classifier TGx-DDI v2 TempOSeqv2 anno rev.txt`: explicit probe crosswalk
  (72 rows including the blank USP41 probe entry).

The original and revised files agree exactly on original symbol order, Entrez IDs,
and genotoxic centroid scores. The measured-symbol updates are C12orf5→TIGAR,
FAM123B→AMER1, HIST1H1E→H1-4, HIST1H2BB→H2BC3, HIST1H2BC→H2BC4,
HIST1H2BG→H2BC8, HIST1H2BI→H2BC10, HIST1H2BM→H2BC14,
HIST1H2BN→H2BC15, and HIST1H3D→H3C4. PDF extraction contains symbol OCR artifacts; original classifier IDs resolve
these against the published panel.

Measured IDs are parsed at the final underscore. Exact symbol matches are checked
against explicit probe IDs. USP41_33775 is a unique measured-symbol match without
a reference probe ID; it is flagged, and a sensitivity result excludes it. No
fuzzy matching or additional alias expansion is performed. No probe maps to two
signature genes. Unmeasured genes and missing reference probes remain explicit.

Six symbols map to multiple measured probes: TIGAR, CBLB, IKBIP, MDM2, PHLDA3,
TP53I3. These are one-to-many probe mappings, not unresolved distinct-gene
identities. A gene counts once if any mapped probe passes. Probe tests retain their
original whole-transcriptome p-values and are never collapsed or retested. Both
gene and feature counts are supplied. The historical classifier separately sums
counts by symbol before CPM calculation.

## Definitions

Every cell model and n=3–10 uses the same 1,000 frozen trials, including zero-call
trials. Signature experiment-level rate means the fraction containing at least
one retained signature gene. The whole-transcriptome comparator uses at least one
retained feature under the same rule. These proportions differ from gene counts
and the proportion of the 64 genes called significant.

- `detected`: any positive input count in selected samples.
- `initial`: total rounded counts >1 in selected samples.
- `relevance`: membership in archived `DESeqRes` after the study CPM rule.
- `tested`: nonmissing raw p-value.
- `bh_eligible`: nonmissing adjusted p-value after independent filtering.
- `significant`: archived BH adjusted p-value <0.05.
- `significant_fc`: significant, absolute shrunken log2 fold change >log2(1.5).
- `quantile`, `spike`: significant probes surviving each rule separately.
- `quantile_spike`: significant probes surviving both rules.
- `final`: both rules plus fold-change cutoff; membership must exactly match
  archived `DESeqResDEGs` in every trial.

Quantile and spike stages are separate, not cumulative. The cumulative path is
significant → quantile → quantile_spike → final. Detection and initial stages
describe availability and have no whole-transcriptome inferential rate. Missing
p-values and adjusted p-values remain visible. `DESeqResFiltered` is never read.
See `docs/deseq2_methods.md` and `R/deg_trials.R` for the historical rules.

Pairing uses embedded classifier `trial_res` IDs, with exact ID-set checks. Saved
probability (>0.9), PC1 (<0), and HCA (DDI) calls remain separate. Their union is
explicitly called `positive_any_component`, a descriptive union of saved calls.
Probability and PC1 values must agree with the frozen compact classifier summary.
The classifier uses B relative to A, whereas DESeq2 uses A relative to B; absolute
DEG fold-change membership requires no reversal of classifier calls.

Rates describe empirical null behavior conditional on these samples and frozen
assignments. They are not conventional BH FDR estimates or independent biological
replicates. Under a complete null, an experiment with calls has false discovery
proportion one; overlapping Monte Carlo draws do not establish general FDR
control. No independence-based confidence intervals or tests are used. Per-gene
frequencies use all 1,000 trials, with filtered or absent genes contributing zero
and availability reported separately.


## Classifier provenance and output files

The raw call files lack HCA fields. HCA is therefore joined by exact trial ID from
the existing frozen `tgx_ddi_monte_carlo_summary.rds`; `raw_hca_available` and
`hca_source` expose this distinction. Raw probability and PC1 records are checked
against that same summary. No classifier component is recomputed.

The primary reproducible entry point is `scripts/tgx_signature/analyze.R`; tests
are separate. Optional `--resume` continues completed-cell checkpoints, while
`--cells` and `--extract-only` support bounded extraction to distinct directories.
The completed analysis preserves the same 32,000 assignments regardless of how
archives are scheduled. The mapping-source SHA-256 manifest is stored in
`data/reference/tgx_ddi_mapping_sources.tsv`.

- `tgx_signature_null.rds`: compact complete results, including per-trial,
  per-gene Boolean stage membership.
- `summary.tsv`: all eleven stage counts/rates and gene-count ranges for each
  cell model and n; strict-probe sensitivity and classifier joint counts.
- `trials_wide.tsv.gz`: all 32,000 trial records, every stage count, and saved
  classifier outcomes with provenance.
- `distribution.tsv`: observed gene-count bins, frequencies and trial proportions;
  omitted bins have frequency zero. Zero-DEG trials are retained.
- `gene_frequencies.tsv`: descriptive per-gene retention frequencies for every
  stage, using the full 1,000-trial denominators.
- `gene_map.tsv`: mapping of all 64 genes to measured probes.
- `gene_map_filter_qc.tsv`: detailed mapping plus per-cell/n/stage retention,
  including all unmapped or filtered entries.
- `classifier_pairings.tsv`: full 2×2 counts for every stage and each classifier
  component, plus their explicitly labelled union.
- `supplementary_table.tsv`: final-call comparisons and exact rates.
- `figures/supplementary/Fig_S13.{pdf,png}`: final-call comparison;
  render from the compact RDS with `scripts/tgx_signature/render_figure.R`.
- `figure_provenance.tsv`: checksums for the compact input, renderer, and both
  figure outputs.
- `provenance.tsv`, `session_info.txt`: exact input fingerprints and environment.

Mapping checks establish consistency with the retained probe annotation. They do
not assess sequence-level probe specificity or resolve cross-hybridization.


Pooled detection across all frozen samples differs from measurement: TK6 detects
63/64 (USP41 absent); MCF-7 detects 64/64; HepaRG detects 61/64 (HIST1H2BI/H2BC10,
LCE1E and USP41 absent); liver spheroids detect 63/64 (USP41 absent). Trial-specific
detection, initial retention, testing and final retention remain separate outputs.


All 32,000 archived classifier input sample/group assignments also match the
frozen metadata, and each input contains the same 64 measured symbols. This is
checked before pairing outcomes; `classifier_input_qc.tsv` records completeness.
The final extraction verifies SHA-256 values of all four DEG archives against
`data/reference/deg_source_archives.tsv`, in addition to fingerprinting every
analysis input and the executable/test scripts.


`tests/check_tgx_signature_outputs.R` additionally writes
`publication_summary_audit.tsv` and `publication_summary_mismatches.tsv.gz`.
Only final DEG counts are required to match the compact publication summary.
Intermediate label differences are disclosed as diagnostics; the archived feature
results and explicit historical rules remain the source for this focused analysis.
See `findings.md` and `validation.md` for the observed outcomes.
