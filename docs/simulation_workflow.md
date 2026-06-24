# Trial-level simulation workflow

The publication repository has two reproducibility layers:

1. **Compact figure workflow**: rebuild the published figures from the compact
   `data/derived/*_monte_carlo_summary.rds` files.
2. **Trial-level simulation workflow**: regenerate the large DESeq2, ORA, GSEA,
   and TGx-DDI trial objects from the frozen Monte Carlo trial assignments.

The compact workflow remains the default because the full simulation workflow is
computationally expensive and produces multi-gigabyte intermediate objects.

## Required inputs

The full simulation workflow expects these files to be present in the repository
layout already documented for the compact workflow:

- `data/input/count_table.tsv`
- `data/input/metadata.QC_applied.txt`
- `data/derived/monte_carlo_trial_metadata.rds`
- `data/external/wikipathways-20231010-gmt-Homo_sapiens.gmt`

TGx-DDI reruns additionally require the classifier table used for the
publication analysis. Because that classifier is not currently bundled in the
repository, pass it explicitly:

```bash
export FDR_TGX_DDI_CLASSIFIER="/path/to/tgx-ddi-classifier.tsv"
```

Run the simulation setup check before launching the expensive jobs:

```bash
Rscript scripts/check_simulation_setup.R
```

## 1. Regenerate DESeq2 trial-level objects

```bash
Rscript scripts/run_deg_trials.R --cores 40 --overwrite true
```

By default this writes:

```text
data/large/DEGs_simulations/tk6.RDS
data/large/DEGs_simulations/mcf7.RDS
data/large/DEGs_simulations/heparg.RDS
data/large/DEGs_simulations/spheroids.RDS
```

Useful options:

```bash
Rscript scripts/run_deg_trials.R --cell-lines "TK6,HepaRG" --cores 20
```

Environment variables are also supported:

- `FDR_DEG_DIR`
- `FDR_COUNT_TABLE`
- `FDR_TRIAL_METADATA`
- `FDR_DEG_ALPHA` default `0.05`
- `FDR_DEG_LINEAR_FC` default `1.5`
- `FDR_NMR_THRESHOLD` default `100000`

Note: the old exploratory script referred to an `R-ODAF`-filtered table, but the
private implementation of that filter was not present in the material used to
build this public workflow. The reproducible DESeq2 rerun therefore records
`DESeqResFiltered` as rows with finite adjusted p-values and records
`DESeqResDEGs` as the p-value + fold-change calls used by the publication figure.

## 2. Regenerate ORA and GSEA trial-level objects

```bash
Rscript scripts/run_pathway_trials.R --analysis both --cores 40 --overwrite true
```

By default this writes GSEA outputs directly under:

```text
data/large/GSEA_simulations/
```

ORA outputs are dated by default:

```text
data/large/ORA_simulations/<YYYY-MM-DD>/
```

To rebuild the compact summaries from a new ORA run, pass the matching directory:

```bash
export FDR_ORA_DIR="data/large/ORA_simulations/2026-06-24"
Rscript scripts/rebuild_summaries_from_large_outputs.R
```

## 3. Regenerate TGx-DDI trial calls

```bash
export FDR_TGX_DDI_CLASSIFIER="/path/to/tgx-ddi-classifier.tsv"
Rscript scripts/run_tgx_ddi_trials.R --cores 40 --overwrite true
```

By default this writes per-cell-line trial calls under:

```text
data/large/TGx_DDI_simulations/
```

It also rebuilds:

```text
data/derived/tgx_ddi_monte_carlo_summary.rds
```

## 4. Rebuild compact summaries and figures

After DESeq2, ORA, GSEA, and optionally TGx-DDI trial-level outputs are present:

```bash
Rscript scripts/rebuild_summaries_from_large_outputs.R
Rscript scripts/render_figures.R
```

If TGx-DDI large outputs are absent, `rebuild_summaries_from_large_outputs.R`
leaves the existing compact TGx-DDI summary unchanged.
