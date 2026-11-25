# Post hoc sensitivity to DEG decision thresholds

## Reproduction

From the repository root, render just the final supplementary screening heatmap
from the bundled rate tables. This writes the PDF and PNG to `figures/`.

```bash
pixi run --locked Rscript scripts/deg_thresholds/plot_threshold_sensitivity.R --only=heatmap
```

To render the final heatmap and all drafts, omit the option.

```bash
pixi run --locked Rscript scripts/deg_thresholds/plot_threshold_sensitivity.R
```

To rebuild the 0.10 screening rates from the original archives, run
`pixi run --locked -e deg-study Rscript scripts/deg_thresholds/extend_screening.R`.
This checks the archived trial IDs and sample assignments against the frozen
metadata and verifies the SHA-256 hash of each archive. The 0.01 and 0.05 rates
come from the original extraction, which can be rebuilt with the following commands.

```bash
pixi run --locked -e deg-study Rscript tests/test_deg_thresholds.R
pixi run --locked -e deg-study Rscript scripts/deg_thresholds/analyze.R \
  --source ../2023_Meier_FDR --output data/derived/deg_thresholds
pixi run --locked -e deg-study Rscript tests/test_deg_thresholds.R
```

This reads one large per-cell-model source archive at a time and releases it
before the next. Small `counts_*.rds` files checkpoint each completed cell model.
`--resume` verifies the input fingerprints, grid and archive checksum before
reusing a checkpoint. No DESeq2 fit, results extraction, shrinkage, classifier
training or BH recalculation is performed. The study environment loads DESeq2
only to deserialize historical S4 tables. Source archives remain read-only.
The source archive manifest, count matrix, frozen assignments, compact publication
summary, methods/code and test script are fingerprinted with SHA-256. Archive
checksums must match `data/reference/deg_source_archives.tsv`.

## Frozen inputs and thresholds

The source tables are the published fitted `DESeqRes` outputs, which have already
passed the relevance filter and contain shrunken log2 fold changes and BH adjusted
p-values extracted using `results(alpha=0.05)` with independent filtering enabled.
The original extraction uses strict adjusted-p-value cutoffs 0.01, 0.025 and 0.05.
The updated figures use 0.01, 0.05 and 0.10 for screening. Both cross the p-value
cutoffs with absolute linear fold-change cutoffs none, 1.25, 1.5 and 2.0. For a finite cutoff
c, selection requires `abs(shrunken log2FoldChange) > log2(c)`: increases and
decreases are treated symmetrically. No exponentiation or signed-linear-FC
approximation is needed. Equality fails both strict cutoffs. Missing adjusted
p-values never pass. A finite fold-change cutoff also excludes missing fold
changes; “none” omits the fold-change condition entirely, retaining zero or
missing fold changes if the other conditions pass. “None” is not a cutoff of 1.

Two endpoints are reported separately:

1. **screening**: adjusted p-value plus optional fold change on the archived
   relevance-retained `DESeqRes` rows.
2. **full_rodaf**: the same selection plus both original quantile and spike rules.

Feature IDs are kept at their original resolution; multiple probes are not
collapsed or retested. `DESeqResFiltered` is never read or treated as final calls.

## Why the final-rule endpoint is identifiable without refitting

The original sample-level counts and frozen assignments are available, and exact
sample/group identity and depth eligibility are checked in every archived trial.
The original normalized counts/size factors are not reconstructed from raw counts.
Instead, the archives preserve the exact feature-aligned quantile and spike
rejection sets in `DESeqjs_filtered_out[[2]]` and `[[3]]`. These rules were evaluated
for **every** feature with stored adjusted p-value <0.05, before fold-change
filtering. The original quantile rule compares group medians with the opposite
group 75th percentile. The spike rule rejects only if both groups satisfy
`max/sum >= 1.4*n^(-0.66)`; an all-zero group does not contribute a spike.

For p cutoffs at or below 0.05, the candidate set is a subset of the original
evaluated set, even when the fold-change cutoff is relaxed or omitted.
Quantile/spike outcomes depend only on the fixed trial-level normalized counts
and groups, not the p or fold-change cutoffs. Thus the complete post hoc final
selection is exactly the intersection of the new stored-p/FC selection and the
complement of the saved rejection sets. This reuses the original sample-level
rule outcomes by feature ID; it does not recalculate normalization, apply rules
to raw counts, or infer them from the problematic recycled intermediate table.
Every original final feature set at p<0.05 and FC>1.5 must match `DESeqResDEGs`.

At 0.10, the archived adjusted p-values and fold changes support screening.
Additional candidates lack saved quantile and spike outcomes, so the backup
heatmap marks the corresponding full-rule rates as unavailable. Extending that
endpoint would require original fitted normalized counts or a separately
justified reconstruction.

## Estimand, checks and limitations

For each endpoint, threshold pair, cell model and n=3–10, the numerator is the
number of frozen null trials with at least one retained feature; the denominator
is all 1,000 frozen trials, including those with zero calls. DEG counts, total
calls, medians, maxima and count distributions are also supplied. No trial is
silently dropped. All 32,000 IDs, iteration coverage, group sizes, exact archived
sample assignments and source fingerprints are checked. Published-cutoff final
feature sets and final compact-summary counts must agree exactly. Comparisons to
literal intermediate compact labels and the actual Figure 2 p-value-panel source
are diagnostic audits rather than alternative sources for the new calculations.
The test suite checks threshold/sign/NA boundaries, no-FC behavior, saved-rule
coverage, portable reference trials, all output denominators and per-trial
monotonicity over p and fold-change cutoffs.

These are **post hoc thresholds on fixed stored outputs**. Changing `alpha` in
`DESeq2::results()` can change the independent-filter threshold, the features
eligible for adjustment and the resulting BH values. All original NA adjusted p-values, the original testing
universe, fitted models, normalization, outlier handling and shrinkage remain
fixed. No alternative BH analysis is included. Fold-change cutoffs are decision
filters, not tests of a nonzero `lfcThreshold`. Accordingly, these results do not
estimate the performance of a completely rerun pipeline at each nominal alpha.

These are empirical experiment-level null rates conditional on the study samples
and frozen, overlapping Monte Carlo assignments. They are not proportions of
genes, general proofs of BH FDR control, or independent biological replicates.
No independence-based tests or confidence intervals are applied. The grid is
exploratory and is not used to select new primary analysis thresholds. Null-only
trials do not assess loss of power or identify an optimal cutoff for detecting
true treatment effects.

## Supplementary figure caption

Empirical false-positive rates after DEG screening at three stored BH-adjusted
p-value cutoffs (0.01, 0.05 and 0.10) and four absolute linear fold-change
cutoffs (none, 1.25, 1.5 and 2.0). Columns show the four cell models and rows show
n = 3–10 samples per group. Tile labels give the percentage of 1,000 frozen null
trials with at least one feature passing both cutoffs. “None” omits the
fold-change condition. Blue outlines mark the published screening cutoffs
(adjusted p < 0.05 and fold change > 1.5). Adjusted p-values and shrunken fold
changes come from the original fits with independent filtering at alpha = 0.05.
