# Study-specific DESeq2 / R-ODAF implementation

The Monte Carlo analysis uses a study-specific implementation of the R-ODAF
rules around DESeq2, rather than a separate statistical model or a modified
DESeq2 package. The executable implementation is `R/deg_trials.R`. This document
describes the executed code, including historical behavior that must be retained
to compare with archived results. It does not claim that every implementation
detail is required by the general R-ODAF framework.

## Analysis sequence and parameters

Each frozen trial assigns DMSO samples to groups A and B. There are four cell
models, eight nominal sample sizes per group (3–10), and 1,000 trials per
combination, totaling 32,000 trials. The assignments in
`data/derived/monte_carlo_trial_metadata.rds` are the source of record.

| Stage | Study implementation |
| --- | --- |
| Input | QC-filtered counts in `data/input/count_table.tsv`; preserve the full feature matrix for each selected trial. |
| Sample inclusion | Retain samples with total counts strictly greater than 100,000, calculated on the input matrix. Remove failing samples from trial metadata before fitting. None of the frozen trial samples fail this cutoff (minimum total: 631,284). |
| Counts | Round counts with R's `round()` before creating the DESeq2 object. |
| Model | Negative binomial model, design `~ group`, without additional covariates; contrast `c("group", "A", "B")`, i.e. A relative to B. |
| Initial feature filter | After object creation, retain genes with total rounded counts strictly greater than 1 across retained trial samples. This is unconditional in the historical script. |
| Fit | `DESeq2::DESeq(dds, parallel = FALSE)` with other arguments left at package defaults: median-ratio size factors (`sfType = "ratio"`), parametric dispersion fitting, Wald tests, and automatic outlier replacement with `minReplicatesForReplace = 7`. |
| Relevance | Calculate DESeq2 size-factor-normalized counts, then `edgeR::cpm()` on that matrix. Retain genes with CPM **≥0.5** in **≥75%** of samples in at least one group. For n samples, this requires at least `ceiling(0.75*n)` samples. |
| Results extraction | Apply the relevance subset to the already fitted object, then call `results()` with the A/B contrast, `alpha = 0.05`, `pAdjustMethod = "fdr"` (BH), and `cooksCutoff = FALSE`. Independent filtering remains enabled (`independentFiltering = TRUE`), and testing uses the default zero log2-fold-change threshold (`lfcThreshold = 0`). Thus the relevance filter precedes results extraction and multiple-testing adjustment, but follows model fitting. |
| Fold-change shrinkage | Call `lfcShrink()` with the fitted object, the extracted results, the same contrast, and `type = "ashr"`. Subsequent fold-change thresholds use the shrunken estimates. |
| Significance | Retain rows with adjusted p-value strictly less than 0.05. Missing adjusted p-values do not pass. Only these genes enter the quantile and spike rules. |
| Quantile rule | Using size-factor-normalized counts, require the median in either group to be strictly greater than the other group's 75th percentile. R's default `quantile()` interpolation is used. |
| Spike rule | Flag a group if `max(counts)/sum(counts) >= 1.4*n^(-0.66)`. An all-zero group does not contribute a spike. Reject a gene only when **both** groups are flagged. |
| Final DEG calls | Significant genes must pass the quantile and spike rules and satisfy `abs(shrunken log2FoldChange) > log2(1.5)`. |

Disabling Cook's distance filtering in `results()` does not explicitly disable
the outlier replacement behavior of `DESeq()`. The implementation preserves
these separate settings rather than equating them.

The rows in `DESeqRes` have already passed the relevance rule and have undergone
shrinkage. Consequently, the publication's “P-Value Only” and “P-Value + Fold
Change” labels describe additional selections from that result table; they do
not denote completely unfiltered DESeq2 analyses.

## Preserved historical behavior

`DESeqResDEGs` contains the final DEG calls described above. Trials with no
significant genes return empty data frames for both filtered output tables.

The intermediate `DESeqResFiltered` table retains a positional subsetting
operation from the historical implementation: a logical vector for significant
genes is applied to the full `DESeqRes` table. R recycles this vector when it is
shorter than the table. Therefore, this field must not be interpreted as a
feature-ID-aligned application of the quantile and spike rules to every tested
gene. It is retained for compatibility with archived intermediate objects and
legacy summaries. The final DEG table applies the vector to the significant
genes themselves. Changing this indexing would be an analysis correction,
rather than a behavior-preserving cleanup.

The old report template describes a different spike formula. The formula in the
table above is the one actually implemented in the Monte Carlo script.

## Provenance

The cleaned implementation was extracted from the locally retained
`scripts/run_DESeq2_bootstrapped.R` in the parent study repository. The original
script remains an internal provenance artifact and is not part of the public
code presentation. A local copy is retained under the Git-ignored
`data/large/provenance/` directory.

- Original script SHA-256:
  `5e58aadac419fe54728ae84f31c33f198f1fa23d3efba21830ac09c9a96c74c6`.
- Parent checkout HEAD at extraction:
  `c94a25d39d66074b31be33a0cd1a295b15854b78`.
- The script had uncommitted modifications. The checksum identifies the extracted
  source; the checkout commit alone does **not** identify that file's contents.
- Source archives: `bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/`
  in the parent study repository, with one RDS file per cell model.

The cleanup separates model fitting, relevance filtering, and post-test rules;
removes plotting setup, global scratch variables, and hard-coded execution
loops; and exposes selected trials through a command-line entry point. It
preserves analytical operations and historical subsetting. Console messages,
scratch globals, and incidental output-list bookkeeping are not preserved.
The simplified DEG implementation previously supplied here is superseded.

## Software versions

The source workflow environment (`workflow/envs/reports.yml`) declares DESeq2
1.38.0, edgeR 3.40.0, and ashr 2.2.54. The selected archived results independently
record ashr 2.2.54 in their shrinkage metadata. A complete original `sessionInfo()`
is not available in these objects, so the declared environment is distinguished
from a verified historical runtime. The environment file SHA-256 is
`30e27da8fefa4f4d6d9e4f38fb4ff4010b88c402e6dcc550fb8c89da7c94b969`;
a local copy is retained in the Git-ignored provenance directory.

Two optional environments are supplied: `deg-study` fixes these three declared
package versions with a compatible R 4.2 series; `deg` uses the current versions
compatible with the figure environment. Both are resolved in `pixi.lock`; the legacy GenomeInfoDbData package additionally uses its checksum-verified installation hook.
The former reconstructs the declared analysis dependencies, not the complete
historical machine or every transitive package version. The latter supports
portability checks and same-environment comparisons of the old and cleaned
functions. Current package versions can change shrunken estimates even when the
analysis code is unchanged.

## Runnable example and validation

The optional locked study environment supplies DESeq2, edgeR, and ashr without
adding them to the default figure environment:

```bash
pixi install --locked -e deg-study --run-post-link-scripts
pixi run --locked -e deg-study Rscript scripts/run_deg_trials.R \
  --trial-ids data/reference/deg_reference_trials.tsv \
  --output-dir data/large/deg_example --cores 1
pixi run --locked -e deg-study Rscript scripts/validate_deg_reference.R --strict-numerics false
```

The install flag enables the historical Bioconductor package's download and
checksum verification of GenomeInfoDbData. If the environment was installed
without this flag and DESeq2 cannot load, complete that dependency with:

```bash
pixi reinstall --locked -e deg-study --run-post-link-scripts bioconductor-genomeinfodbdata
```

The example retains all genes and runs only the selected frozen assignments.
Omitting `--trial-ids` invokes the expensive full trial workflow; this is not
needed for figure reproduction or the bounded validation.

`data/reference/deg_reference.rds` contains unchanged selected archived output
tables converted to ordinary data frames for portability. Selection includes
the first n=3 and n=10 trial per cell model, supplemented by the first trial
with final DEGs, quantile rejection, and spike rejection when available and not
already selected. The TSV records the selection reasons. These are targeted
regression examples, not a random sample for estimating overall concordance.
The current reference set contains 17 trials (224,375 result rows in total).
No gene-level input subset is used.

Validation checks feature membership, numeric result columns, final DEG calls,
filter counts, excluded-feature tables, and sample assignments. Numeric
comparisons allow `abs(observed - expected) / max(1, abs(expected)) <= 1e-6`;
missing/nonfinite values must match. Feature IDs and calls must match exactly.
Results and `sessionInfo()` are written to `data/large/deg_validation/`.
The example above requires exact feature/call membership, filter counts, sample
metadata, and schema agreement, while reporting archive numeric-column
differences separately. Omitting `--strict-numerics false` requires every numeric
comparison to pass as well and returns a nonzero status on numerical differences.
The source-function comparison always requires exact numeric agreement,
regardless of this flag. The numerical tolerance remains 1e-6 in both modes.

For internal verification against the original source in the **same** current
software environment, add:

```bash
pixi run --locked -e deg Rscript scripts/validate_deg_reference.R --strict-numerics false \
  --original-source data/large/provenance/run_DESeq2_bootstrapped.R
```

The validator parses and evaluates only the original trial function. It never
sources the top-level script or launches its full simulation loops. A separate
comparison with archived output detects historical software or source drift;
same-environment agreement alone cannot establish agreement with the archived
analysis. The same-environment numeric comparison requires exact agreement
(tolerance zero).
Neither check establishes equivalence for every one of 32,000 trials.

To reconstruct the small reference set from the complete archives:

```bash
pixi run --locked -e deg-study Rscript scripts/extract_deg_reference.R \
  --archive-dir /path/to/DEGs_simulations
```

This reads the large archived objects one cell model at a time and does not
refit models. The checked-in report `docs/deg_validation.md` records the actual
validation outcome and version limitations.
