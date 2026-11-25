# Supplementary Figure S8: baseline mean–variance pattern

## Data and original analysis

The original plot was produced in the sibling `2023_Meier_FDR` analysis from
`Rmd/Process bootstrapped results.Rmd`. It summarized baseline expression in
DMSO control wells and joined each probe's mean and standard deviation to the
number of historical DMSO-versus-DMSO trials in which that probe passed the
BH-adjusted p-value threshold. The recurrence endpoint is `padj < 0.05` without
a fold-change filter; it is not the later FINAL DEG endpoint.

For each probe and cell model, expression is `log2(CPM + 1)` across all
QC-applied DMSO wells. The original code calculated CPM from DESeq2-normalized
counts. The port uses raw counts and divides by each sample's raw library size;
this gives the same CPM because the DESeq2 size factor is a constant multiplier
within a sample and cancels in the CPM denominator. The plotted SD is the sample
standard deviation (`stats::sd`).

Each `n = 3`–`10` panel shows probes with at least one DEG call among the 1,000
trials at that sample size, colored by the number of such trials. As in the
original left join, probes with no call at any sample size appear in the final
“Never” column; zero recurrence counts are not imputed into each individual
sample-size panel.

## Revised visualization and interpretation

Each panel has one LOESS smoother (`stats::loess`, span 0.75). The confidence
band is omitted. The smoother is descriptive: no correlation coefficient,
p-value, or inferential test is calculated for this figure. This is intentional
because the mean–variance pattern is curved and non-monotonic, so a single
correlation coefficient would not summarize it appropriately.

Figure definition

> **Figure S8. Baseline mean expression and variability in DMSO controls.** Each
> point is an assay probe. Mean log2 CPM and sample standard deviation of log2
> CPM were calculated across QC-applied DMSO wells separately for each cell
> model. Panels `n = 3`–`10` show probes called differentially expressed
> (BH-adjusted p < 0.05, without a fold-change filter) in at least one of 1,000
> DMSO-versus-DMSO trials at that per-group sample size; color indicates the
> number of trials. “Never” denotes probes with no DEG calls across the tested
> sample sizes. The solid line is a descriptive LOESS smoother (span = 0.75); no
> correlation coefficient or inferential test was calculated.

## Reproduction

The compact plotting data are stored in
`data/derived/s8_mean_variance_summary.rds`, so the figure can be rendered from
the tracked repository inputs with:

```bash
pixi run s8-figure
```

The full figure suite also renders it with `pixi run figures`. To rebuild the
compact summary, run `pixi run s8-data`; that rebuild requires the four local
recurrence caches in `data/large/gene_recurrence/` (or set
`FDR_S8_RECURRENCE_DIR`). Those caches retain the original trial IDs and source
archive checksums. The resulting provenance is written to
`data/derived/s8_mean_variance_provenance.tsv`; figure-render inputs and output
checksums are recorded in `data/derived/s8_figure_provenance.tsv`.
