# Post hoc DEG threshold sensitivity

The analysis includes all 32,000 frozen mock trials (1,000 in each cell-model/n
stratum), 12 threshold combinations, and two separate endpoints. All numbers
below use stored BH values and shrunken log2 fold changes. No BH recalculation,
DESeq2 results extraction, model fitting or shrinkage was performed.

For a compact descriptive comparison, pooling the eight equally sized n strata
within each cell model gives the following **full final-rule** rates. Each entry
uses 8,000 trials. The complete per-n numerators and denominators are in
`data/derived/deg_thresholds/threshold_rates.tsv`. The two-page heatmap draft
shows the per-n rates at the revised cutoffs.

| Cell model | Stored BH p<0.05, no FC cutoff | Published p<0.05, FC>1.5 | Stored p<0.01, FC>2.0 |
|---|---:|---:|---:|
| TK6 | 1,056/8,000 (13.20%) | 706/8,000 (8.825%) | 449/8,000 (5.6125%) |
| MCF-7 | 1,086/8,000 (13.575%) | 388/8,000 (4.85%) | 183/8,000 (2.2875%) |
| HepaRG | 1,209/8,000 (15.1125%) | 248/8,000 (3.10%) | 69/8,000 (0.8625%) |
| Liver spheroids | 1,915/8,000 (23.9375%) | 860/8,000 (10.75%) | 527/8,000 (6.5875%) |
| All models | 5,266/32,000 (16.45625%) | 2,202/32,000 (6.88125%) | 1,228/32,000 (3.8375%) |

The corresponding **p-plus-FC screening** pooled rates across all models are
5,316/32,000 (16.6125%) at p<0.05/no FC cutoff, 2,231/32,000 (6.971875%) at
p<0.05/FC>1.5, and 1,248/32,000 (3.90%) at p<0.01/FC>2.0. Screening and full
final-rule calls are distinct even when their experiment-level rates are close.
Feature counts and complete count distributions are also supplied. “No FC cutoff”
omits that condition, rather than requiring a nonzero fold change.

Every individual trial has nonincreasing counts as the p or FC cutoff becomes
stricter. This is a nested-selection check, not an independence-based hypothesis
test. Rates still depend on n and cell model; these null-only data do not measure
power or justify selecting new primary thresholds. Pooling here is descriptive,
with the same equal weights across n; it does not replace the per-n results.

## Published-cutoff and Figure 2 audit

All 32,000 original full final feature sets match archived `DESeqResDEGs` exactly.
Their counts also match the compact publication summary exactly. At the two
literal intermediate compact labels, direct selections differ in 2,151 trial
counts for `P-Value Only` and 622 for `P-Value + Fold Change`.

The submitted figure scripts
rename `R-ODAF + P-Value` to `P-value filter only`. This actual p-value-panel source
has **2,286 trial feature-count discrepancies** against direct selection of
archived `DESeqRes` rows at stored BH p<0.05 with no FC cutoff. Of those,
**57 trials change the >=1-call status**: the legacy panel source has 5,259 positive
trials, while direct screening has 5,316. The maximum cell-model/n rate difference
is 1.0 percentage point (TK6, n=9: 11.3% versus 12.3%). This is a reproducible
panel-source/label inconsistency, separate from the threshold sensitivity and
from the unaffected final-call panel. The source trace identifies NA-index placeholder counting and two positional
subsetting errors.

The primary fitted inputs and archived summaries are retained. Correct per-trial
p-only counts are supplied in `trial_counts.tsv.gz`, column
`screening__p0.05_fcnone`. The three fixed-cutoff stages are the inputs for
the corrected revised Figure 2.
`publication_cutoff_comparison.tsv` gives all per-cell/n comparison numerators,
and `publication_count_mismatches.tsv.gz` lists every affected trial ID and both
counts for each compared label.

## Reproducible deliverables

- `scripts/deg_thresholds/analyze.R` produces the original threshold extraction.
- `scripts/deg_thresholds/extend_screening.R` adds 0.10 screening rates in
  `threshold_screening_p010.tsv`.
- `scripts/deg_thresholds/plot_threshold_sensitivity.R` renders the single-panel
  screening heatmap as Figure S4. Run `pixi run --locked threshold-drafts`
  to generate the alternate layouts under `figures/exploratory/`.
  The line-chart drafts vary fold change along the horizontal axis, flip that
  layout, or vary adjusted p-value along the horizontal axis with a line for
  each fold-change cutoff. All figures use adjusted p-value cutoffs 0.01, 0.05
  and 0.10.
  The two-page draft marks 0.10 full-rule rates as unavailable.
- `tests/test_deg_thresholds.R`: boundary/reference and full-result tests.
- `threshold_sensitivity.rds`: complete summary, trial counts, grid, QC,
  distributions, publication audit and SHA-256 provenance.
- `threshold_rates.tsv`: all 768 endpoint/cell/n/threshold results, with exact
  numerators, denominators, total DEG counts, medians and maxima.
- `trial_counts.tsv.gz`: all 32,000 IDs and 24 endpoint/threshold counts.
- `deg_count_distributions.tsv`, `trial_qc.tsv`: count distributions, fixed
  availability and original post-rule exclusion tallies.
- `provenance.tsv`, `validation.tsv`, `session_info.txt`: provenance and checks.

See `methods.md` for exact selection rules, why original post-rule outcomes can
be carried forward for this grid, and the fixed-independent-filtering limitation.
