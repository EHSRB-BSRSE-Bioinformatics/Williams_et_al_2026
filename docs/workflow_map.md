# Publication workflow map

This repository packages the false-positive analysis in Williams et al. (2026)
as two reproducibility levels: a compact figure workflow that runs entirely
from Git, and an expensive re-summarization workflow that also requires the
separate large-file archive.

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
for reproducing the publication plots.

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

## Expensive re-summarization workflow

```text
data/large/DEGs_simulations/*.RDS ------------------+
data/large/GSEA_simulations/*_gsea.RDS -------------+--> compact summaries
data/large/ORA_simulations/2023-11-09/*_ora.RDS ----+          |
                                                               v
                                                        figure workflow
```

The large trial-level objects will be distributed through Zenodo. Their exact
directory layout and filenames are documented in `data/large/README.md`. Run
`pixi run rebuild-summaries` after extracting them. The TGx-DDI compact summary
is included in Git; its large provenance objects are archived for traceability
but are not consumed by the current rebuild script.

For GSEA, parallel `foreach` results may be unnamed. The rebuild code restores
trial names from the frozen trial design, validates existing names when they
are present, and refuses to continue when result counts do not match. The
publication GSEA summary uses gene sets containing 15--500 genes.

## False-positive rate definition

For each cell line and sample size, the empirical false-positive rate is:

```text
trials with at least one false positive / observed trials
```

The code records both the observed and expected trial counts and warns when an
analysis group does not contain 1,000 trials.

## Integrity metadata

`docs/data_manifest.tsv` records the role and SHA-256 checksum of each included
data object. Recalculate the affected checksum whenever an included data file
is deliberately replaced; unexplained checksum changes should block release.

