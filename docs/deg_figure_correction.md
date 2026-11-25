# Corrected DEG Monte Carlo figure

The files `figures/Fig_2.pdf` and `.png` now use
feature-aligned archived counts at the original decision thresholds. This is a
correction to counting and display, with unchanged fitted results, BH values,
independent filtering, shrinkage and final DEG rules. The DEG figure uses separately derived plotting inputs. Original figure files and pre-correction source snapshots are kept
under the gitignored `data/large/deg_figure_correction/pre_correction/` directory.

## Where the discrepancies originated

The retained `../2023_Meier_FDR/Rmd/Process bootstrapped results.Rmd` reveals two
bookkeeping errors and a later plotting-label substitution.

First, the literal p-only summary selected rows with:

```r
res[res$padj < 0.05, ]
```

When an adjusted p-value is missing, that comparison is NA. R data-frame indexing
creates a placeholder row for an NA selector; it does not silently discard the
row. The old code then counted row names, including those placeholders. Across
all 32,000 trials, its literal `P-Value Only` count equals the actual significant
count **plus the number of missing adjusted p-values**, exactly. For example,
`HepaRG_3_3_5` has 23 significant features and 12,435 missing adjusted p-values;
the legacy count is 12,458. Those missing values are not significant genes.
The p-plus-FC summary has the analogous NA-indexing problem whenever the combined
logical selector is NA. Re-evaluating the original expressions reproduces all
three audited legacy labels in all 17 portable reference trials.

Second, the historical `DESeqResFiltered` intermediate applies a logical vector
sized to significant genes to the full result table. R recycles that shorter
vector across unrelated rows. The downstream summary then applies another mask
from that intermediate to the full table:

```r
res[filtered$padj < 0.05, ]
```

These rows are selected by position rather than matched feature IDs. The result
is not a valid p-only count or a valid final-rule count. A concrete reference
trial, `Liver spheroids_3_3_334`, has 3 BH-significant features and 2 final DEGs,
but the legacy p-panel count is 1.

Finally, the figure code renamed `R-ODAF + P-Value` to `P-value filter only`, so
this positional selection was displayed as if it were direct p-only screening.
That displayed source differs from correct direct screening in 2,286 trial
feature counts and 57 >=1-call statuses. Its positive-trial total was 5,259;
the correct total is 5,316. The largest rate change is 1 percentage point
(TK6 n=9: 11.3% to 12.3%). Count-color changes can be much larger than rate
changes, because thousands of erroneous placeholder rows still contribute only
one positive experiment.

The final DEG table uses the significant subset correctly and is separate from
these intermediate summaries. Its feature identities and counts agree with the
fixed-rule reconstruction in every one of the 32,000 trials; the total remains
2,202 positive trials.

## Corrected data path and reproduction

`R/summaries.R::read_deg_figure_summary()` constructs three analysis stages
from the validated per-trial archive extraction in
`data/derived/deg_thresholds/threshold_sensitivity.rds`:

- Panel A: `screening__p0.05_fcnone`, direct nonmissing stored BH p<0.05 on the
  archived relevance-retained `DESeqRes` rows, without an FC cutoff.
- Panel B: `screening__p0.05_fc1.5`, direct screening of the same rows with
  nonmissing stored shrunken absolute log2 fold-change > log2(1.5).
- Final DEG counts in the saved data: `full_rodaf__p0.05_fc1.5`, verified
  original full final calls, including the quantile/spike rules.

The positive-trial totals are 5,316, 2,231 and 2,202 of 32,000, respectively
(16.61%, 6.97% and 6.88%). These are proportions of mock trials with at least
one call, not proportions of genes or conventional BH FDR estimates. The
additional quantile/spike rules remove all remaining calls in 29 trials
after p-plus-FC screening. The overlapping Monte Carlo trials are descriptive
resamples; these totals do not assume independent experiments.

All IDs and all 1,000-trial strata are validated. Heatmap ranks are recomputed
within each cell/model, sample size and endpoint from these counts. The figure
builder rejects the unvalidated legacy summary. Both heatmaps use the same
count-color scale. Titles state the fixed thresholds.

```bash
pixi run --locked Rscript tests/test_deg_figure.R
pixi run --locked Rscript scripts/render_figures.R --only=deg
```

The targeted rendering command loads only DEG plotting inputs. The standard full
render also uses the corrected reader. The original
`data/derived/deg_monte_carlo_summary.rds` remains unchanged for provenance;
corrected counts and rates for all three stages are exported to
`data/derived/deg_figure_summary.rds` and `deg_figure_rates.tsv`.

The dedicated tests verify all three stage counts and per-trial nesting, every ID and denominator,
within-stratum ranks, all 64 plotted rates, shared color limits, rejection of
legacy relabeling, and concrete historical counterexamples. The existing
`scripts/check_analysis_consistency.R` checks also pass with the corrected figure
path. Figure input/output fingerprints are recorded in
`data/derived/deg_figure_provenance.tsv`.
