# DEG count reductions across successive filters

The companion `figures/deg_filter_burden_supplementary.pdf` (also PNG) shows
model-by-sample-size totals across each condition's 1,000 frozen trials.
The upper panels use one common logarithmic color scale and label the exact
number of DEG calls at each stage: stored BH p<0.05; additionally absolute
shrunken log2 fold-change>log2(1.5); additionally the archived quantile/spike
rules. Lower panels show the percentage of preceding-stage calls removed,
with exact numbers removed in parentheses. Stages are nested feature sets;
the display order does not imply that the original pipeline applied FC first.

A DEG call is one feature in one trial. A feature appearing in multiple trials
contributes repeatedly. These are neither unique-gene totals nor conventional
BH FDR estimates. Trial assignments overlap; no independence-based inference
or confidence intervals are applied. Large-call trials can dominate totals;
use this figure alongside the existing trial-level count distributions and
>=1-DEG rates, rather than as a replacement for them.

Across 32,000 trials, totals decrease from 925,785 to 187,166 to 185,673 calls.
FC screening removes 738,619 calls (79.78%); the additional quantile/spike rules
remove 1,493 calls (0.80% of calls remaining after FC). Quantile/spike filtering
reduces counts in 337 trials: 308 retain at least one DEG and 29 lose their last
DEG. Thus the experiment-level positive-trial rate hides most affected trials,
but the overall extra reduction in DEG calls is also small. Most quantile/spike
removals occur in TK6 (1,449 calls); MCF-7 contributes 1, HepaRG 0, and liver
spheroids 43. Some condition-specific percentages reflect very small call
denominators; exact counts are displayed to make this visible.

The TSV transition table provides before/after totals, calls removed, positive
trial counts, and a paired trial-state decomposition: reduced but still
positive; reduced to zero; unchanged positive; unchanged zero. The proportion
of DEG calls removed uses calls as its denominator. Trial counts use all
1,000 frozen trials per condition; they are a distinct unit of measurement.

Reproduce from the repository root, using the existing validated compact counts:

```bash
pixi run --locked Rscript scripts/deg_thresholds/plot_filter_burden.R
pixi run --locked Rscript tests/test_deg_filter_burden.R
```

The script reads through `read_deg_figure_summary()`, which checks all trial IDs,
denominators, nested counts and agreement with archived final calls. No model
fitting, BH recalculation or changes to the primary figure are performed.
Tables and a compact RDS are saved as `data/derived/deg_thresholds/filter_burden*`;
the provenance TSV contains SHA-256 hashes of source counts, scripts and outputs.
Tests check a controlled partial-removal example, complete removal, invalid
counts, zero denominators and agreement of saved condition and pooled totals.

## Per-trial means and conditional medians

`figures/deg_filter_medians_supplementary.pdf` (also PNG) provides a separate
per-trial figure. The upper row shows the arithmetic mean across all 1,000
trials per condition, including zeros: total calls divided by 1,000. This
replaces the uninformative all-trial median, which is zero at every stage.
The lower row conditions on
trials with at least one DEG at the displayed stage, with each tile showing its
positive-trial denominator. Both rows share the same log(count+1) color scale;
legend labels give the untransformed DEG count. Grey/NA would indicate no
positive trials, rather than a conditional median of zero.

The conditional subset can change with filtering. Its median can increase as
low-count trials lose their last call, even though every paired trial's count
is nonincreasing. Conditional medians are descriptive counts, not paired
reductions, experiment-level rates or BH FDR estimates. Median reductions must
not be inferred by subtracting these conditional medians.

For TK6 at sample size 10, the successive totals are 48,759, 9,453 and 9,209
across 1,000 trials, with 109, 73 and 70 positive trials. The corresponding
all-trial medians are 0, 0 and 0; conditional medians are 6, 14 and 27. Thus the
final 9,209 total corresponds to a mean of 9.209 calls per trial and a conditional
median of 27, rather than 9,209 DEGs in a typical trial. The maximum final count
in this condition is 1,047; quartiles and maxima are included in the TSV.

Reproduce independently of the totals figure:

```bash
pixi run --locked Rscript scripts/deg_thresholds/plot_filter_medians.R
pixi run --locked Rscript tests/test_deg_filter_medians.R
```

Results are saved in `data/derived/deg_thresholds/filter_medians.tsv` and `.rds`,
with hashes in `filter_medians_provenance.tsv`. Tests verify arithmetic means equal totals divided by trial counts, inclusion of zero
trials, conditional denominators, empty-positive handling, fractional medians,
and exact agreement of all 96 saved summaries with the validated trial counts.

The arithmetic mean directly measures average calls per mock trial and can be
influenced by trials with many DEGs. A conventional geometric mean is zero
whenever any count is zero. The shifted geometric mean exp(mean(log(count+1)))-1
would be a different statistic dependent on the added pseudocount; it is not
used in this figure. The existing median filenames are retained for continuity.
The previous all-zero top-panel version is backed up under the gitignored
`data/large/deg_figure_correction/pre_mean_panel/` directory.
