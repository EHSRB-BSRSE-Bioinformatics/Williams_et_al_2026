# Supplementary Figures S5 and S6: recurrent screening calls

## Scope

These figures reproduce the historical **BH-adjusted p < 0.05 only** screening
analysis. They do not use the full FINAL DEG rule, a fold-change cutoff, or the
pooled FINAL top-ten gene lists. The latter lists were a separate exploratory
analysis and must not be substituted into these captions or their Results text.
Historical recurrence notes called these figures S4/S5; the current manuscript
calls them S5/S6.

## Selection and recurrence (S5)

Each of 32,000 frozen DMSO-versus-DMSO trials contributes the archived probe IDs
with a non-missing BH-adjusted p < 0.05. Appearances are counted separately for
each model and per-group sample size n = 3–10, with 1,000 trials per stratum.
No p-values or DEG models are refitted.

The historical algorithm applies `top_n(10, wt = count)` to **probe-by-n rows**
grouped by cell model, retains boundary ties, and then takes distinct probe IDs.
It does not sum counts across sample sizes before ranking, and does not select
ten distinct genes at every sample size. This selects 11 probes in TK6, 10 in
MCF-7, 12 in HepaRG, and 10 in liver spheroids: 43 model/probe selections and 41
unique probes. Each of these 41 probes has a distinct displayed gene symbol.

The union of selected probes is displayed across all four models. A probe can
therefore appear in a model's heatmap without having been selected by that
model. Gene symbols are display labels, not the units of statistical testing.
The machine-readable selections are in
[s56_selected_probes.tsv](../data/derived/s56_selected_probes.tsv).

The new plot retains the original mean-over-nonzero-rows y ordering. Unlike the
historical sparse raster, it distinguishes **black = tested, zero screening
calls** from **grey = no eligible DESeq2 test at that sample size**. Eligibility
means at least one trial with a non-missing archived unadjusted p-value. The
number of tested and BH-eligible trials is retained for every selected
probe/model/n in
[screening_eligibility.rds](../data/reference/s56/screening_eligibility.rds).
Grey does not establish biological absence: every selected probe has at least
some raw reads in every model, including low-count background observations.
No new expression threshold is imposed. Model order is fixed to TK6, MCF-7,
HepaRG, liver spheroids.

## Baseline expression (S6)

S6 uses the same selected probes and 301 QC-applied DMSO wells: 77 TK6, 77 MCF-7,
73 HepaRG and 74 liver spheroids. Each well appears once per displayed probe.
The plotted values are **log2(DESeq2-normalized count + 1)**, not log2 CPM.

The historical plot loaded `ddsList$all` from an archived expression workspace.
This object contains 365 samples; its size factors were already estimated before
restricting the plot to DMSO wells. They must not be re-estimated within each
cell model or from the selected probes. The archive has no gene-specific
normalization-factor matrix. Its DMSO raw counts were verified to match the
public count matrix exactly. Consequently, the historical selected-probe
expression values can be reconstructed from public raw counts and the archived
DMSO size factors in [size_factors.tsv](../data/reference/s56/size_factors.tsv).
Extraction verifies agreement with the archived normalized values at tolerance
1e-12; normal rebuilding needs neither DESeq2 nor the archived workspace.

Boxes show medians, quartiles and 1.5-IQR whiskers; individual wells are shown
with deterministic horizontal jitter. A red star means at least one screening
call for that probe/model at any n. It does **not** mark a treatment effect,
correlation, enrichment, or selection into that model's top set. There is one
star per probe/model, placed above that probe's expression range rather than
repeated once per well at a global maximum as in the historical code.

The underlying selections, counts and normalized values are preserved. Layout,
eligibility/zero display, deterministic jitter and star placement are corrected for public
rendering; the PNG files are not byte-identical copies of the Word images.

## Descriptive baseline correlations among S5-selected probes

Baseline pairwise Spearman correlations are computed separately within each
model among **that model's historically selected probes** (11/10/12/10), not
among the pooled FINAL leaders and not across different cell models. Each
eligible DMSO well is used once. The expression scale is the same archived
normalized-count scale as S6. No new figure or pathway enrichment analysis is
introduced.

| Model | Probes | Wells | Pairs | Median absolute rho | Signed rho range |
|---|---:|---:|---:|---:|---|
| TK6 | 11 | 77 | 55 | 0.754 | 0.641 to 0.873 |
| MCF-7 | 10 | 77 | 45 | 0.929 | -0.968 to 0.979 |
| HepaRG | 12 | 73 | 66 | 0.843 | -0.891 to 0.940 |
| Liver spheroids | 10 | 74 | 45 | 0.880 | -0.924 to 0.933 |

Values are generated in
[s56_baseline_correlation_summary.tsv](../data/derived/s56_baseline_correlation_summary.tsv);
the 211 signed pairwise values are in
[s56_baseline_correlations.tsv](../data/derived/s56_baseline_correlations.tsv).
These are **unadjusted descriptive correlations**, conditional on selecting
recurrent probes in the same finite dataset. Normalization, library properties,
plate effects and selection can contribute to the observed covariance. No
matched-background comparison, significance test, causal attribution, independent
replication or coordinated treatment-response claim is made. A large absolute
correlation can be negative: it must not be described as uniform co-upregulation.

S6's expression distributions alone do not measure these pairwise correlations.
Functional examples such as histone genes are not formal pathway enrichment of
this screening list. The separate FINAL-list Reactome and adjusted-correlation
results must not be attached to these probes.

## Reproduction

From the repository root:

```bash
pixi run --locked s56-data
pixi run --locked s56-check
pixi run --locked s56-figures
```

`s56-data` rebuilds the summary and correlation tables using tracked raw counts,
metadata, complete screening-count reference, selected-probe testing eligibility,
and archived size factors.
`s56-figures` validates the summary and writes PDF and 600-dpi PNG files.
`pixi run figures` also includes both figures.

One-time reference extraction, if recovering or intentionally replacing the
reference inputs, requires the original archives' four recurrence caches and
the expression workspace:

```bash
pixi run --locked -e deg-study Rscript scripts/rebuild_s56_figure_data.R --extract-reference
```

Set `FDR_S56_RECURRENCE_DIR`, `FDR_ORIGINAL_ANALYSIS` or
`FDR_S56_EXPRESSION_WORKSPACE` to relocate those inputs. Normal public rendering
and rebuilding do not depend on them. Extraction reconstructs screening counts
from the cached exact probe memberships and checks frozen trial IDs. The cache
creation and original archive hashes are already documented by the retained
recurrence/S8 provenance; this command does not re-extract multi-gigabyte DEG
archives or refit simulations.

Reference source hashes are in
[source_provenance.tsv](../data/reference/s56/source_provenance.tsv).
Data-rebuild and render hashes are in
[s56_data_provenance.tsv](../data/derived/s56_data_provenance.tsv) and
[s56_figure_provenance.tsv](../data/derived/s56_figure_provenance.tsv).
Tests verify exact selected identities, boundary-tie behavior, the full heatmap
grid, normalization, star membership, independently recomputed correlations,
and plot builds.

## Figure interpretation

**Supplementary Figure S5. Recurrent false-positive screening calls in
DMSO-versus-DMSO mock experiments.** Recurrence is the number of trials in which
an assay probe had BH-adjusted p < 0.05, without a fold-change requirement.
Counts are shown separately for each cell model and per-group sample size
(n = 3–10; 1,000 trials per stratum). Within each model, probes associated with
the ten highest probe-by-sample-size recurrence rows were selected, retaining
boundary ties, before taking distinct probe IDs. This selected 11 probes in TK6,
10 in MCF-7, 12 in HepaRG and 10 in liver spheroids, comprising 41 unique probes.
Their counts are displayed across all four models; gene symbols label the
individual probes. Black denotes zero screening calls when the probe was
eligible for testing in at least one trial at that sample size; grey denotes no
eligible test (no non-missing archived unadjusted p-value). Grey is not proof of
biological absence. This historical
screening selection is distinct from a pooled ranking under the full final DEG
rule.

**Supplementary Figure S6. Baseline expression of the probes selected in
Supplementary Figure S5.** Each point represents one QC-applied DMSO well
(77 TK6, 77 MCF-7, 73 HepaRG and 74 liver spheroids). Expression is
log2(DESeq2-normalized count + 1), using the archived study-wide size factors.
Boxes show the median and interquartile range; whiskers extend to the most
extreme observations within 1.5 times the interquartile range. Panels have
probe-specific y scales. A red asterisk indicates that the probe had at least
one BH-adjusted-p-only screening call in that model across the tested sample
sizes; it does not indicate a treatment response or a gene-gene correlation.
