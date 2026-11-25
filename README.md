# Empirical false-positive rates in transcriptomics

This repository is the publication companion for Williams et al. (2026),
*Empirically determined false discovery rates for differentially expressed
genes, benchmark concentration modeling, transcriptomic biomarker activation,
and pathway analysis in four human cell models*.

It contains the code, frozen Monte Carlo trial assignments, compact derived
results, source count data, and rendered figures needed to inspect and reproduce
the false-positive analyses based on DMSO controls. The study used 32,000
two-group comparisons for differential expression, pathway analysis, and
TGx-DDI, plus 4,800 mock-dose experiments for benchmark concentration
modeling and transcriptomic points of departure (tPoDs).

## Reproduce the figures

The locked [Pixi](https://pixi.sh/) environment currently supports Linux x86-64.
From the repository root, run:

```bash
pixi run --locked figures
```

The task validates the compact inputs and writes PNG/PDF pairs for Figures 2–6
under `figures/` and Figures S1–S13 under `figures/supplementary/`. It runs
the DMSO PCA (S2) in the `deg` environment. Run `pixi run --locked threshold-drafts`
for alternate threshold layouts or `pixi run --locked exploratory-figures`
for the full set of optional plots. The figure manifest lists every output and
its producer. Figure 1 was not supplied. [Figure manifest](docs/figure_manifest.tsv)

| Figures | Targeted render | Compact rebuild or methods |
| --- | --- | --- |
| Main figures and summary-based plots | `pixi run --locked publication-summary-figures` | [Workflow map](docs/workflow_map.md) |
| Figure 3 (GO, Reactome and WikiPathways ORA) | `pixi run --locked Rscript scripts/render_figures.R --only=ora` | [Frozen Reactome ORA workflow](docs/revision/reactome_ora_v2/REPRODUCE.md) |
| Figure 4 (unique-Entrez Reactome GSEA) | `pixi run --locked Rscript scripts/render_figures.R --only=gsea` | [Frozen ranking and GSEA workflow](docs/revision/gsea_unique_ids_v1/REPRODUCE.md) |
| Figure 6 (mock-dose tPoDs) | `pixi run --locked Rscript scripts/render_figures.R --only=tpod_fdr_main` | [BMDExpress workflow and runnable examples](docs/tpod_workflow.md) |
| S1 and S3 | `pixi run --locked sample-figures` | `pixi run --locked sample-data` · [Methods](docs/sample_supplementary_figures.md) |
| S2 | `pixi run --locked -e deg dmso-pca` | [PCA inputs](docs/workflow_map.md#compact-figure-workflow) |
| S4 | `pixi run --locked Rscript scripts/deg_thresholds/plot_threshold_sensitivity.R --only=heatmap` | [Workflow map](docs/workflow_map.md) |
| S5 and S6 | `pixi run --locked s56-figures` | `pixi run --locked s56-data` · [Methods](docs/supplementary_figures_s5_s6.md) |
| S7 and S9–S12 | `pixi run --locked legacy-figures` | `pixi run --locked legacy-data` · [Historical definitions](docs/legacy_supplementary_figures.md) |
| S8 | `pixi run --locked s8-figure` | `pixi run --locked s8-data` with local recurrence caches · [Methods](docs/supplementary_figure_s8.md) |
| S13 | `pixi run --locked tgx-signature-figure` | [TGx-DDI signature](docs/tgx_signature/README.md) |

S3 uses the frozen four-well-per-group DMSO comparisons matched to the archived
draft figure for all four cell models. [S3 selection and provenance](docs/sample_supplementary_figures.md#recovery-evidence-and-deterministic-selection)

For a subset of the main plots, run:

```bash
pixi run --locked Rscript scripts/render_figures.R --only=deg,ora,gsea
```

If Pixi is unavailable, install the plotting packages listed in `DESCRIPTION`
and DESeq2 for S2, then run the producer scripts in the figure manifest from
the repository root.

## Repository contents

- `R/`: reusable analysis and plotting functions.
- `scripts/`: validation, figure rendering, and summary-rebuild entry points.
- `data/input/`: QC-filtered count matrix and sample metadata.
- `data/derived/`: frozen trial assignments and compact analysis summaries.
- `data/external/`: the dated WikiPathways gene-set release used for ORA.
- `data/large/`: documented location for the separately archived trial-level
  objects; downloaded objects are intentionally ignored by Git.
- `figures/`: main, supplementary, and exploratory outputs.
- `docs/figure_manifest.tsv`: canonical filenames and producer scripts.
- `docs/data_manifest.tsv`: provenance, inclusion status, and SHA-256 checksums.
- `docs/workflow_map.md`: data flow and rerun policy.

## Reproducibility scope

The compact summaries are sufficient to regenerate the included figures. The
exact 32,000 Monte Carlo sample assignments are frozen in
`data/derived/monte_carlo_trial_metadata.rds`; they can also be regenerated with
seed 1 using:

```bash
pixi run trial-metadata
```

The full DESeq2, GSEA, ORA, and TGx-DDI trial-level objects are multi-gigabyte
intermediates and are not stored in Git. They will be deposited in Zenodo. Until
that deposit is complete, the compact summaries in this repository are the
supported reproduction path. After downloading the archive, extract it under
`data/large/` using the layout in `data/large/README.md`.

The figure commands read bundled results. Selected summaries can also be rebuilt
from tracked references, as indicated in the table. The archive-dependent
trial recipes require the original inputs and are computationally expensive.
Figure 3 combines archived GO and WikiPathways results with frozen Reactome ORA.
Figure 4 uses frozen unique-Entrez Reactome rankings and a 15–500 gene-set
testing family. The generic pathway runner provides a separate simulation route.

The mock-dose workflow uses its own DMSO assignments. The repository includes
scripts to generate illustrative BMDExpress inputs and summarize exported fits.
Figure 6 renders from the bundled tPoD summary. [Mock-dose workflow](docs/tpod_workflow.md)

The study-specific DESeq2/R-ODAF implementation, parameters, filtering order,
and provenance are documented in [docs/deseq2_methods.md](docs/deseq2_methods.md).
A small set of archived trial outputs is bundled in `data/reference/` for bounded
validation with the full gene matrix:

```bash
pixi install --locked -e deg-study --run-post-link-scripts
pixi run --locked -e deg-study Rscript scripts/validate_deg_reference.R --strict-numerics false
```

To run just these selected frozen trials and save their outputs:

```bash
pixi run --locked -e deg-study Rscript scripts/run_deg_trials.R \
  --trial-ids data/reference/deg_reference_trials.tsv \
  --output-dir data/large/deg_example --cores 1
```

The optional `deg-study` environment fixes the study's declared DESeq2, edgeR,
and ashr versions. It is separate from the default figure environment.
See [docs/deg_validation.md](docs/deg_validation.md) for the measured validation
outcome and historical-version limitations. The full original pipeline and its
report-generation dependencies are not bundled.

The approximately 23 MB count table is deliberately stored directly in Git so
the compact figure workflow remains self-contained. The multi-gigabyte
trial-level objects are the only files delegated to Zenodo.

## Data and citation status

The Git repository is the versioned source for code and compact files. The
Zenodo DOI and the final publication citation will be added after those records
exist. At that point, cite both the paper and the archived repository release.
See the [publication companion repository](https://github.com/EHSRB-BSRSE-Bioinformatics/Williams_et_al_2026)
for the versioned releases.

The code, documentation, and study-generated data are distributed under the
[MIT License](LICENSE). Terms for bundled third-party material are recorded in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Contributing

This is an archive-style research repository. Reproducibility corrections are
welcome; see `CONTRIBUTING.md` before opening an issue or pull request.
