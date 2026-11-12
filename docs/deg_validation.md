# Selected-trial DEG validation

Validation was performed on 12 November 2026 using the 17 frozen trials listed
in `data/reference/deg_reference_trials.tsv`. Each fit retained the full input
feature matrix. The reference set covers all four cell models, n=3 and n=10,
zero-significance cases, nonzero final DEG calls, and quantile/spike rejection.
There are 224,375 archived `DESeqRes` rows across the selected trials. These
targeted examples do not estimate concordance across all 32,000 trials.

## Code equivalence in the current environment

The cleaned function and the locally retained original trial function were each
run on the same selected inputs using R 4.5.3, DESeq2 1.50.2, edgeR 4.10.7, and
ashr 2.2.63. **All 703 comparisons passed**, with zero numerical tolerance.
Checked fields included the main results, historical filtered intermediate,
final DEG table, filter counts, excluded-feature tables, sample IDs/groups,
and retained sample counts. Numeric columns and feature membership matched
exactly. This supports behavior preservation for these selected trials.

The original source was parsed to evaluate only its trial function. The script's
full top-level workflow was never executed. The source checksum and checkout
status are documented in `deseq2_methods.md`.

## Current software compared with the archived analysis

Using the same current environment, **684 of 703 archive comparisons passed**.
All feature-membership comparisons, final DEG memberships, and filter counts
matched for all 17 trials. Nineteen numeric-column checks exceeded the archive
tolerance in three trials: `MCF-7_3_3_1`, `HepaRG_3_3_1`, and
`Liver spheroids_3_3_1`. Differences affected shrunken estimates/standard errors
and, in two trials, p-values or adjusted p-values. The largest absolute
difference among failed checks was about 0.0190 in `lfcSE` for HepaRG.

The archive tolerance was `abs(observed - expected) / max(1, abs(expected)) <=
1e-6`, with matching missing/nonfinite values required. This strict numerical comparison exited with status 1.
Current software is therefore not presented as an exact replay of archived
numerical estimates. Agreement between old and cleaned functions in this same
current environment separates this drift from the cleanup itself.

## Declared study-version replay

The `deg-study` environment fixes DESeq2 1.38.0, edgeR 3.40.0, and ashr 2.2.54,
as declared in the source workflow environment. The selected archived objects
also record ashr 2.2.54. The resolved R version is 4.2.3. This reconstructs the
three declared analysis packages, rather than claiming recovery of the complete
historical session or all transitive dependencies.

**683 of 703 numerical/structural comparisons passed at tolerance 1e-6.**
All feature memberships, final DEG memberships, and filter counts matched across
all 17 trials. Twenty numeric-column comparisons exceeded tolerance in four
trials: `MCF-7_3_3_1`, `MCF-7_10_10_1`, `HepaRG_3_3_1`, and
`Liver spheroids_3_3_1`. The largest failed absolute difference was about
0.000325 in p-values for `MCF-7_3_3_1`; the HepaRG shrinkage differences were
reduced to approximately 8.3e-6 or smaller. The strict check returned status 1.

Thus the declared package versions improve some agreement but do not establish
an exact historical runtime. Remaining differences may reflect transitive
software, numerical libraries, or historical source/runtime details that were
not captured. Their cause has not been isolated. The old and cleaned functions
agree exactly when run in the same tested current environment, which is the
specific evidence for preserving analysis semantics during cleanup.

## Decision checks and strict numerical checks

For the reviewer's requested methods transparency, the documented command uses
`--strict-numerics false`. It requires matching feature membership, final calls,
filter counts, schemas, ordered sample IDs/groups, and retained sample counts,
while retaining numeric estimate differences in the report. This mode passed
using the 17 saved selected-trial CLI outputs, without fitting them again.

Strict numeric checking remains the script default; omit that flag to require
every archived numeric column to agree within the stated tolerance. Both strict
archive runs above failed and are reported as such. The original-function
comparison always requires exact numerics, even with the relaxed archive flag.
The numerical tolerance remains 1e-6 in both modes; exact historical numerical
replay has not been established.

## Other verification

- The reference extractor completed successfully against the complete archived
  files. A repeat extraction reproduced all three reference files byte for byte.
  Source-archive and reference-file SHA-256 checksums are recorded in
  `data/reference/deg_source_archives.tsv` and `docs/data_manifest.tsv`.
- The selected-trial CLI completed using two worker processes and produced
  exactly 17 trial outputs across four per-cell-model files.
- The CPM, quantile, and spike boundary checks passed; they are also included
  in the existing lightweight CI workflow.
- The compact repository setup check passed, all R scripts parsed, and
  `git diff --check` passed. The default figure environment's lock entries were
  unchanged. Frozen trial assignments, publication summaries, and figures were
  not regenerated.
- The complete 32,000-trial analysis was not rerun.

## Evidence and commands

- `validation/deg_current_comparisons.tsv`: all current-environment comparisons,
  including both original-function and archive checks.
- `validation/deg_current_session_info.txt`: current validation parameters and
  complete session information.
- `validation/deg_study_comparisons.tsv` and
  `validation/deg_study_session_info.txt`: strict declared-version replay.
- `validation/deg_decision_comparisons.tsv` and
  `validation/deg_decision_session_info.txt`: successful decision check using
  saved current-environment CLI outputs.

```bash
# Archive comparison using the declared study versions.
pixi install --locked -e deg-study --run-post-link-scripts
pixi run --locked -e deg-study Rscript scripts/validate_deg_reference.R --strict-numerics false

# Internal exact code-equivalence check plus archive decision checks
# under current software. Original source is retained locally, outside Git.
pixi run --locked -e deg Rscript scripts/validate_deg_reference.R --strict-numerics false \
  --original-source data/large/provenance/run_DESeq2_bootstrapped.R \
  --output-dir data/large/deg_validation_current

Rscript scripts/check_deg_rules.R
```
