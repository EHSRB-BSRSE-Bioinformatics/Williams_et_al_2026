# Empirical false-positive rates in transcriptomics

This repository is the publication companion for Williams et al. (2026),
*Empirically determined false discovery rates for differentially expressed
genes, benchmark concentration modeling, transcriptomic biomarker activation,
and pathway analysis in four human cell models*.

It contains the code, frozen Monte Carlo trial assignments, compact derived
results, source count data, and rendered figures needed to inspect and reproduce
the false-positive analyses based on DMSO-versus-DMSO trials.

## Reproduce the figures

The locked [Pixi](https://pixi.sh/) environment currently supports Linux x86-64.
From the repository root:

```bash
pixi run check
pixi run figures
```

`pixi run check` validates the required files, R packages, sample identifiers,
frozen trial design, and compact-summary schemas. `pixi run figures` writes PDF
and PNG versions of all four figures to `figures/`.

If Pixi is unavailable, install the R packages listed in `DESCRIPTION`, then run:

```bash
Rscript scripts/check_setup.R
Rscript scripts/render_figures.R
```

## Repository contents

- `R/`: reusable analysis and plotting functions.
- `scripts/`: validation, figure rendering, and summary-rebuild entry points.
- `data/input/`: QC-filtered count matrix and sample metadata.
- `data/derived/`: frozen trial assignments and compact analysis summaries.
- `data/external/`: the dated WikiPathways gene-set release used for ORA.
- `data/large/`: documented location for the separately archived trial-level
  objects; downloaded objects are intentionally ignored by Git.
- `figures/`: versioned publication outputs generated from compact summaries.
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

Full reruns are computationally expensive and are only needed when the frozen
trial design, filtering thresholds, or annotation versions change. The
publication GSEA analysis uses gene sets containing 15–500 genes.

The upstream workflow that originally generated the trial-level DESeq2,
enrichment, and biomarker objects is not included here. Its additional software
dependencies are therefore intentionally excluded from the locked environment.

The approximately 23 MB count table is deliberately stored directly in Git so
the compact figure workflow remains self-contained. The multi-gigabyte
trial-level objects are the only files delegated to Zenodo.

## Data and citation status

The Git repository is the versioned source for code and compact files. The
Zenodo DOI and the final publication citation will be added after those records
exist. At that point, cite both the paper and the archived repository release.
The remaining publication and archive steps are tracked in
`RELEASE_CHECKLIST.md`.

The repository's code and documentation are licensed under the MIT License; see
`LICENSE`. The bundled WikiPathways data has separate terms recorded in
`THIRD_PARTY_NOTICES.md`. Confirm the distribution terms for the study-generated
data with the publication authors before the public release.

## Contributing

This is an archive-style research repository. Reproducibility corrections are
welcome; see `CONTRIBUTING.md` before opening an issue or pull request.
