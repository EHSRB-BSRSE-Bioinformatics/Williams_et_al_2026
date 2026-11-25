# Validation record

Completed on 2026-11-17 on branch `analysis/deg-posthoc-thresholds`, using the
locked `deg-study` environment. The archive extraction read one cell model at a
time. The four large archive SHA-256 values match the frozen source manifest;
count-data, trial-assignment, code/method and summary fingerprints are recorded.
Inputs were fingerprinted before and after extraction and checked for changes.

Unit checks pass for strict p-value and FC boundaries, symmetric positive and
negative fold changes, NA/NaN handling, zero and missing fold changes under no-FC
selection, rejection of unsupported p cutoffs, duplicate-ID rejection, correct
cached quantile/spike application and independence from `DESeqResFiltered`.
All 17 portable archived reference trials reproduce their original final counts.

The complete extraction validates all 32,000 frozen IDs and exact sample/group
assignments, group sizes and iteration coverage, depth eligibility and measured
feature IDs. Original final feature membership at stored p<0.05/FC>1.5 matches the
archive in every trial. Final compact-summary counts match in every trial.

The full output suite passes all 768 stratum denominators (1,000 each), all
32,000-trial per-threshold count comparisons for endpoint inclusion and
monotonicity, and all distribution totals. Original screening and final cutoff
counts are explicitly separated from legacy intermediate-label comparisons.
The Figure 2 source audit finds 2,286 count discrepancies and 57 positive-status
discrepancies in the p-value panel; the final panel agrees exactly. The sensitivity extraction preserved primary figures and analyses.

Both supplementary heatmap pages were rendered for visual inspection. A
rendering-only change removed a PDF character-conversion warning and separated
axis/header labels. All non-rendering function bodies were checked identical
before and after that change; the source fingerprint and compact caches were
updated accordingly. The final rendering completes with warnings treated as
errors and no warnings. Numeric tables were unchanged by this layout adjustment.
