# Reproduce the controlled GSEA calculation

Run from the repository root. Large original DEG archives and the recovered
historical pathway list are required at the paths specified by the runner. This workflow
does not use DESeqResDEGs, refit models or regenerate GO/WikiPathways.

```bash
CONDA_PACKAGE_CACHE=/path/to/conda/pkgs python scripts/revision/setup_gsea_unique_env.py
.pixi/envs/deg-study/bin/Rscript scripts/revision/inspect_gsea_unique_env.R
python scripts/revision/verify_gsea_unique_sources.py
.pixi/envs/deg-study/bin/Rscript scripts/revision/check_gsea_unique_rankings.R
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GSEA_UNIQUE_CORES=24 \
  .pixi/envs/deg-study/bin/Rscript scripts/revision/run_gsea_unique_ids.R --prepare
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GSEA_UNIQUE_CORES=32 \
  .pixi/envs/deg-study/bin/Rscript scripts/revision/run_gsea_unique_ids.R
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GSEA_UNIQUE_CORES=32 \
  .pixi/envs/deg-study/bin/Rscript scripts/revision/run_gsea_unique_ids.R --full
.pixi/envs/default/bin/Rscript scripts/revision/summarise_gsea_unique_ids.R
```

Before running the commands, recover the ordered `pathways.rds` and
`short_anno.rds` using `scripts/revision/inspect_historical_workspace.R`.
The GSEA summarizer also consumes the Reactome ORA v2 raw checkpoints.
Place the [fgsea 1.26.0 source tarball](https://bioconductor.org/packages/3.17/bioc/src/contrib/fgsea_1.26.0.tar.gz)
under `data/large/reactome_ora_v1/annotation/` and check it against the
bundled `fgsea_source.sha256`.

The setup helper consumes cached dependencies and builds fgsea 1.26.0 into
`data/large/gsea_unique_ids_v1/library`. For reconstruction on another machine,
create a dedicated environment from `conda-explicit.txt` and build the preserved
fgsea source into its R library; the `deg-study` R executable path may be
adjusted to that environment. The explicit artifact lock pins all dependencies;
the preserved source tarball checksum pins fgsea. Parent R versions are checked
by the runner. The runner records source and output hashes in its checkpoint records.

Rankings and raw results are in `data/large/gsea_unique_ids_v1/{rankings,full}`;
the repeated same-environment pilot is `pilot.rds`. Compact summaries, trial
diagnostics, comparisons and logs are retained in this directory. The runner
refuses full enrichment until a validated complete pilot is present. A
`--references-only` timing run cannot satisfy this gate. Repeating `--prepare`
or `--full` checks and reuses matching completed checkpoints.

## Direct minimum-size sensitivity

With the original ranking checkpoints and pinned GSEA environment already
available, run from the repository root:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GSEA_UNIQUE_CORES=32 \
  flock -n data/large/gsea_unique_ids_minsize15_v1.lock \
  .pixi/envs/deg-study/bin/Rscript scripts/revision/run_gsea_minsize15.R --full
```

This passes `minSize = 15` and `maxSize = 500` directly to
`fgseaMultilevel()` on the 32,000 frozen unique-Entrez rankings. It reuses
the recovered Reactome pathways and the original per-trial seeds. BH
adjustment uses the pathways with defined raw p-values in each trial.

The run writes under `data/large/gsea_unique_ids_minsize15_v1/full/` and
checks and reuses completed checkpoints when rerun with the same command.
The raw block RDS files contain pathway results and warnings.
`trial_results.tsv` pairs the direct results with the post hoc Figure 4
counts by trial ID. `summary.tsv` reports the positive-trial comparisons
by model and sample size and overall. `pathway_counts.tsv` gives counts
of eligible, tested and significant trials per pathway, model and sample
size. The outputs are for sensitivity analysis and comparison with Figure 4.

Use `--pilot` in place of `--full` to run the same procedure on the
33 frozen pilot trials. Its outputs are under
`data/large/gsea_unique_ids_minsize15_v1/pilot/`.
