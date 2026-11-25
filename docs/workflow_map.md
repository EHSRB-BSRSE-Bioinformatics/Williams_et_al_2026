# Publication workflow map

This repository packages the false-positive analysis in Williams et al. (2026)
as two reproducibility levels: a compact figure workflow that runs entirely from
Git, and an expensive trial-level workflow that regenerates the large simulation
objects before rebuilding compact summaries.

## Compact figure workflow

Run `pixi run --locked figures` from the repository root. The task validates
tracked inputs, runs the summary and standalone producers, and calls the S2
PCA producer in the `deg` environment. It writes the publication PNG/PDF pairs
to `figures/` and `figures/supplementary/`. The optional `threshold-drafts`
task also writes alternate threshold plots under `figures/exploratory/`.
The figure manifest records the paths and producers. [Figure manifest](figure_manifest.tsv)

S2 uses QC-passing DMSO samples and colors the PCA by cell model. The same
treatment applies to all samples and batch is confounded with cell model in
this subset. To render S2 alone, run `pixi run --locked -e deg dmso-pca`.

S1/S3 and S7/S9–S12 each have a compact rebuild from frozen references.
Use `sample-data` or `legacy-data` to rebuild, followed by the corresponding
`sample-figures` or `legacy-figures` task to render. [S1/S3 methods](sample_supplementary_figures.md)
[Historical figure methods](legacy_supplementary_figures.md)

Supplementary Figure S8 follows the same compact-render pattern. Its tracked
summary contains the baseline expression statistics and sparse per-probe DEG
recurrence counts; `scripts/render_figures.R` adds the per-panel descriptive
LOESS smoothers. Run `pixi run s8-figure` to render only S8.

```text
data/input/count_table.tsv + data/input/metadata.QC_applied.txt
                 + data/large/gene_recurrence/*.rds
                                  |
                                  v
                    scripts/rebuild_s8_figure_data.R
                                  |
                                  v
         data/derived/s8_mean_variance_summary.rds (tracked)
                                  |
                                  v
                        scripts/render_figures.R
                                  |
                                  v
             figures/supplementary/Fig_S8.*
```

The summary rebuild requires the four local recurrence caches. It records the
source DEG archive hashes and analysis inputs in
`data/derived/s8_mean_variance_provenance.tsv`; normal figure rendering uses the
tracked summary and does not need those caches.

Supplementary Figures S5/S6 now have a fully tracked compact rebuild path:

```text
data/reference/s56/screening_counts.rds + size_factors.tsv
       + data/input/count_table.tsv + metadata.QC_applied.txt
          |
          v
        scripts/rebuild_s56_figure_data.R
          |
          v
       data/derived/s56_figure_summary.rds
       selected-probe and baseline-correlation TSVs
          |
          v
      scripts/render_figures.R
          |
          v
  figures/supplementary/Fig_S5.*
  figures/supplementary/Fig_S6.*
```

Run `pixi run s56-data`, `pixi run s56-check`, and `pixi run s56-figures`.
Unlike S8's cache-dependent rebuild, the S5/S6 rebuild uses tracked complete
screening-count references and archived DMSO size factors. Optional reference
re-extraction needs the original recurrence caches and expression workspace.
See [S5/S6 methods](supplementary_figures_s5_s6.md) for endpoint, normalization,
selection, and descriptive-correlation limitations.

## Source inputs and frozen trial design

- `data/input/count_table.tsv` is the QC-filtered count matrix used by the
  DESeq2 Monte Carlo trials.
- `data/input/metadata.QC_applied.txt` is the post-QC sample metadata. Trials
  select samples where `chemical == "DMSO"`.
- `data/external/wikipathways-20231010-gmt-Homo_sapiens.gmt` is the frozen
  WikiPathways release used by ORA.
- `data/derived/monte_carlo_trial_metadata.rds` contains 32,000 frozen trial
  definitions: four cell models, eight per-group sample sizes, and 1,000 trials
  per combination.

`scripts/01_make_trial_metadata.R` regenerates trial definitions with seed 1 by
default and writes a separate seed-labelled file. The frozen object remains the
source of record for the publication analysis.

## Trial-level simulation workflow

The diagram follows the generic simulation entry points. The publication
Reactome ORA and GSEA recipes pin annotation and trial inputs separately.
[Reactome ORA recipe](revision/reactome_ora_v2/REPRODUCE.md) ·
[Unique-Entrez GSEA recipe](revision/gsea_unique_ids_v1/REPRODUCE.md)

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

## Mock-dose benchmark concentration workflow

Four cell models, twelve mock-dose designs, and 100 experiments per design
yield 4,800 experiments. The two-group trial design above has its own sample
assignments. The mock-dose workflow starts from the DMSO controls and uses
BMDExpress for benchmark concentration modeling.

```text
Count table + sample metadata
  -> scripts/tpod/prepare_mock_doses.R
  -> mock-dose input files
  -> BMDExpress imports, fits, and exports (external)
  -> scripts/tpod/summarize_fpr.R
  -> per-design tPoD false-positive rate summaries

Bundled data/derived/tpod_fdr_grid_summary.rds
  -> scripts/render_figures.R --only=tpod_fdr_main
  -> Figure 6
```

The mock-dose methods page documents the original scripts, the inputs needed
to rerun the analysis, and two runnable examples. [Mock-dose workflow](tpod_workflow.md)

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

For each cell line and design, the empirical false-positive rate is:

```text
trials with at least one false positive / observed trials
```

The two-group summary code records both the observed and expected trial counts
and warns when an analysis group does not contain 1,000 trials. The tPoD export
summarizer expects 100 experiments per cell model and design by default and
checks the number of analysis sets before computing rates.

## Integrity metadata

`docs/data_manifest.tsv` records the role and SHA-256 checksum of each included
data object. Recalculate the affected checksum whenever an included data file is
deliberately replaced. Run `python3 scripts/check_data_manifest.py` to check
the tracked release files. Archived extraction provenance preserves the hashes
of the code and data used at extraction time. Review the separate figure
provenance tables for hashes of the code and outputs used for rendering.
