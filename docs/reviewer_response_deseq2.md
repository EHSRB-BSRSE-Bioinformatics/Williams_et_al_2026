# Suggested manuscript Methods text

Differential expression was assessed with DESeq2 using a study-specific
implementation of the R-ODAF filtering rules. The source analysis environment
specified DESeq2 1.38.0, edgeR 3.40.0, and ashr 2.2.54; the archived results
independently record ashr 2.2.54. Each frozen DMSO-versus-DMSO trial was fitted
with design `~ group`, using group A relative to group B. Samples required more
than 100,000 total counts, and features required more than one total rounded
count across the selected samples before fitting. DESeq2 model fitting otherwise
used package defaults. Size-factor-normalized counts were converted to CPM with
edgeR, and features were retained for results extraction if at least 75% of
samples in either group had CPM ≥0.5. Results were extracted with alpha 0.05,
Benjamini–Hochberg adjustment, and Cook's distance filtering disabled;
independent filtering and model-stage outlier replacement retained their
package defaults. Log2 fold changes were shrunk using ashr. Final DEGs required
adjusted p <0.05, absolute shrunken log2 fold change >log2(1.5), and passage of
the quantile and spike rules. The quantile rule required the median normalized
count in either group to exceed the other group's 75th percentile. A group was
flagged for a spike when its maximum normalized count divided by its summed
normalized counts was ≥1.4n^(-0.66); features were excluded by this rule only
when both groups were flagged. Executable code, detailed filtering order,
historical implementation notes, selected archived trial outputs, and validation
results are provided in the accompanying repository.

Use the manuscript's existing DESeq2 and R-ODAF citations and add the ashr
citation if it is not already present. The version wording intentionally
separates declared dependencies from a complete historical session record.

The historical `DESeqResFiltered` intermediate has positional logical recycling;
it is not a feature-aligned application of all rules to every tested feature.
See `deseq2_methods.md` before describing legacy “R-ODAF Only” intermediate
summaries in a manuscript. The final DEG table uses the significant-gene rows.

# Suggested reviewer-response wording

We have added a cleaned, runnable reference implementation of the study-specific
DESeq2/R-ODAF workflow, together with explicit analysis settings, filtering order,
software-version provenance, and a bounded validation set of 17 frozen Monte
Carlo trials covering all four cell models. The original script is retained
internally for provenance. The public code preserves the historical analysis
behavior, including documented intermediate-table indexing. Validation compares
the cleaned and original trial functions in the same software environment and
separately compares selected reruns with archived trial-level outputs. All 703
same-environment source-function checks passed with exact numerical agreement. Selected reruns preserved archived feature membership, final DEG
membership, and filter counts in both current and declared study-version
environments. Some numerical estimates differed from the archived values;
these residual differences are documented rather than presented as exact
historical replay. Detailed outcomes are recorded in the accompanying
validation report. The complete 32,000-trial analysis was not rerun.
