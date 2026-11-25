# Mock-dose benchmark concentration analysis

The tPoD analysis used four cell models, 12 mock-dose designs, and 100 DMSO
experiments per model and design, for a total of 4,800 experiments. Each
experiment assigned DMSO samples to nominal doses, including dose zero.
Designs combined 4, 6, 8, or 12 dose groups with 24, 48, or 72 samples in
total. For example, `d6n4` denotes six groups with four replicates each.

The repository entry points use base R and share one design definition and
one parameter grid across cell models. To summarize exported fits, supply
the probe annotation and dated pathway memberships used for those fits.

## Prepare mock-dose inputs

Run from the repository root with the bundled count table and post-QC metadata.
This command writes one import file for HepaRG to a new output directory.

```bash
Rscript scripts/tpod/prepare_mock_doses.R \
  --cell-lines HepaRG --designs d6n4 --n-trials 1 \
  --output-dir data/large/tpod/example_inputs
```

The file `heparg/d6n4/Exp 1.txt` starts with a `Dose` row followed by probe
expression rows. The preparation script screens libraries below 10% of the
median library size in the supplied count table, matches sample IDs to the
supplied metadata, and selects eligible DMSO samples. It computes
`log2(1e6 * (count + 0.5) / (library_size + 1))` and floors the result at
`log2(5)`. It resets seed 90214 for each cell model and draws samples without
replacement within each experiment, using the original design order. With
the same number of trials, selecting a subset of designs preserves their
position in the random number stream.

Omit the selection flags to generate 100 inputs per model and design. Pass
`--count-table` and `--metadata` to use other source files. The bundled
metadata is the post-QC copy, while the supplied script refers to a local
`metadata.txt`. The original sample assignments are unavailable here, so
the generated inputs illustrate the design with the supplied data.

In the supplied `d6n8` code, positions 31 and 32 are assigned to dose 1 and
then overwritten by dose 10. The generator follows those assignments, giving
six samples at dose 1 and ten at dose 10. The bundled summary alone cannot
establish which assignments went into the historical BMDExpress fits.

## Summarize BMDExpress exports

Fit the mock-dose inputs in BMDExpress and export a tab-separated table for
each cell model and design. The R summarizer takes one such export and applies
all eight settings of the Williams trend-test p-value, absolute fold change,
and R² grid.

```bash
Rscript scripts/tpod/summarize_fpr.R \
  --export data/large/tpod/Export_HepaRG_d6n4.txt \
  --annotation data/large/tpod/probe_manifest.csv \
  --pathways data/large/tpod/pathway_memberships.tsv \
  --cell-line HepaRG --design d6n4 \
  --output data/large/tpod/heparg_d6n4_fpr.csv
```

- The export needs `Analysis`, `Probe ID`, `Best BMD`, `Best BMDL`,
  `Best BMDU`, `Best RSquared`, and `Best fitPValue`. Columns matching
  `Prefilter.*P-Value` and `Max Fold Change Absolute Value` supply the two
  prefilters. As in the collaborator's code, an absent prefilter field
  skips that criterion. The entry point warns when this happens.
- The probe manifest uses the original `Probe_Name` and `Enzrez_ID`
  column names. Pathway memberships use `collection`, `pathway`, and
  `entrez_id`, with one row per membership across KEGG, Reactome, Hallmark,
  and GO. Omit `--pathways` to report the gene endpoint alone.
- The summarizer checks for 100 distinct analysis sets by default. Use
  `--expected-trials` for a deliberately smaller export. Include an analysis
  row for every experiment, including experiments with no qualifying fits,
  so they contribute to the denominator.

The trend-test cutoffs are strict (`p < 0.05` or `p < 0.01`), and absolute
fold change must meet or exceed 1.5 or 2.0. The fit rules require
`0 < BMD < 100`, `fit p > 0`, `BMDL > 0`, `BMD/BMDL < 20`,
and `BMDU/BMDL < 40`. The R² cutoffs are inclusive (`>= 0.6` or `>= 0`).
The CSV retains the historical labels `R-Sq > 0.6` and `No R-Sq`.

An experiment counts as positive for the gene tPoD when at least 25 unique
Entrez genes pass the filters. A pathway endpoint counts as positive when
any pathway in that collection has at least 30 annotated genes, at least
three passing genes, and at least 5% of its annotated genes passing. The
reported FPR is the fraction of experiments with a positive endpoint.

For a runnable post-processing example, use two synthetic analysis sets.
At the 1.5-fold, `p < 0.05` settings, one experiment passes each of the five
endpoints, yielding an FPR of 0.5. No experiment passes at the stricter trend
p-value or fold-change cutoff.

```bash
Rscript scripts/tpod/example.R
```

## Render Figure 6

Figure 6 uses the bundled
[`tpod_fdr_grid_summary.rds`](../data/derived/tpod_fdr_grid_summary.rds).
Rendering uses the repository's bundled figure inputs and standard figure environment.
The supplied R scripts cover input preparation and post-processing. The original
BMDExpress project, fit configuration, exports, probe manifest, and dated
pathway snapshots are unavailable in the repository, so an end-to-end match
to the bundled summary has not been verified.

```bash
pixi run --locked Rscript scripts/render_figures.R --only=tpod_fdr_main
```
