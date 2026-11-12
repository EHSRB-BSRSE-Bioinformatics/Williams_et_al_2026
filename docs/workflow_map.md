# Publication workflow map

This repository packages the false-positive analysis in Williams et al. (2026)
as two reproducibility levels: a compact figure workflow that runs entirely from
Git, and an expensive trial-level workflow that regenerates the large simulation
objects before rebuilding compact summaries.

## Compact figure workflow

```text
data/derived/*_monte_carlo_summary.rds
                    |
                    v
        scripts/render_figures.R
                    |
                    v
             figures/*.pdf
             figures/*.png
```

Run `pixi run check` and then `pixi run figures`. This is the supported route
for reproducing the publication plots from the files stored directly in Git.

## Source inputs and frozen trial design

- `data/input/count_table.tsv` is the QC-filtered count matrix used by the
  DESeq2 Monte Carlo trials.
- `data/input/metadata.QC_applied.txt` is the post-QC sample metadata. Trials
  select samples where `chemical == "DMSO"`.
- `data/external/wikipathways-20231010-gmt-Homo_sapiens.gmt` is the frozen
  WikiPathways release used by ORA.
- `data/derived/monte_carlo_trial_metadata.rds` contains 32,000 frozen trial
  definitions: four cell lines, eight per-group sample sizes, and 1,000 trials
  per combination.

`scripts/01_make_trial_metadata.R` regenerates trial definitions with seed 1 by
default and writes a separate seed-labelled file. The frozen object remains the
source of record for the publication analysis.

## Trial-level simulation workflow

```text
data/input/count_table.tsv
                 +
data/derived/monte_carlo_trial_metadata.rds
                 |
                 v
        scripts/run_deg_trials.R
                 |
                 v
data/large/DEGs_simulations/*.RDS
       |                         \
       |                          \
       v                           v
scripts/run_pathway_trials.R   scripts/run_tgx_ddi_trials.R
       |                           |
       v                           v
data/large/ORA_simulations/     data/large/TGx_DDI_simulations/
data/large/GSEA_simulations/              |
       \___________________________       |
                                   \      |
                                    v     v
                scripts/rebuild_summaries_from_large_outputs.R
                                      |
                                      v
                    data/derived/*_monte_carlo_summary.rds
                                      |
                                      v
                            scripts/render_figures.R
```

The full workflow is documented in `docs/simulation_workflow.md`. The exact
study-specific DESeq2/R-ODAF sequence and a bounded validation path are documented
in `docs/deseq2_methods.md`. It is
computationally expensive and creates multi-gigabyte intermediate objects, so the
large outputs are intentionally kept out of Git and distributed separately.

## Expensive re-summarization workflow

```text
data/large/DEGs_simulations/*.RDS -------------------+
data/large/GSEA_simulations/*_gsea.RDS --------------+--> compact summaries
data/large/ORA_simulations/*_ora.RDS ----------------+          |
data/large/TGx_DDI_simulations/*_tgx_ddi_results.RDS-+          v
                                                        figure workflow
```

Run `pixi run rebuild-summaries` or
`Rscript scripts/rebuild_summaries_from_large_outputs.R` after extracting the
large archive. Set `FDR_ORA_DIR` when rebuilding from a newly dated ORA run.
TGx-DDI large outputs are optional for re-summarization; when absent, the
existing compact TGx-DDI summary is left unchanged.

For GSEA, parallel results may be unnamed. The rebuild code restores trial names
from the frozen trial design, validates existing names when they are present,
and refuses to continue when result counts do not match. The publication GSEA
summary uses gene sets containing 15--500 genes.

## False-positive rate definition

For each cell line and sample size, the empirical false-positive rate is:

```text
trials with at least one false positive / observed trials
```

The code records both the observed and expected trial counts and warns when an
analysis group does not contain 1,000 trials.

## Integrity metadata

`docs/data_manifest.tsv` records the role and SHA-256 checksum of each included
data object. Recalculate the affected checksum whenever an included data file is
deliberately replaced; unexplained checksum changes should block release.
