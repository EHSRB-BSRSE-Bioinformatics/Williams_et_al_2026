# Reactome ORA reproduction

The primary analysis uses archive-compatible Reactome 1.84.0 memberships and
org.Hs.eg.db 3.17.0 alias mappings. The companion v1 analysis retains the
Reactome 1.77.0 database-version sensitivity. See the
[v1 instructions](../reactome_ora_v1/REPRODUCE.md) for the original DEG/GSEA
archive locations and locked Pixi environments.

Run from the repository root:

```bash
mkdir -p data/large/reactome_ora_v1/annotation data/large/reactome_ora_v2/annotation
cp docs/revision/reactome_ora_v1/frozen_annotation_org3.17.rds data/large/reactome_ora_v1/annotation/
cp docs/revision/reactome_ora_v2/frozen_annotation.rds data/large/reactome_ora_v2/annotation/
.pixi/envs/deg-study/bin/Rscript scripts/revision/check_reactome_ora.R
OPENBLAS_NUM_THREADS=1 .pixi/envs/deg-study/bin/Rscript scripts/revision/run_reactome_ora.R --full
OPENBLAS_NUM_THREADS=1 .pixi/envs/deg-study/bin/Rscript scripts/revision/retest_reactome_collection.R
.pixi/envs/default/bin/Rscript scripts/revision/summarise_reactome_ora.R
.pixi/envs/default/bin/Rscript scripts/revision/check_reactome_summary.R
.pixi/envs/default/bin/Rscript scripts/revision/plot_reactome_ora.R
.pixi/envs/default/bin/Rscript scripts/revision/cache_gsea_comparison_inputs.R
.pixi/envs/default/bin/Rscript scripts/revision/reactome_common_family_v2.R
```

V2 reads the inputs prepared by v1. The analysis reuses the frozen DESeq2
results and trial assignments. It accounts for all 32,000 trials and retains
insufficient-input trials in the denominator. There were 557 positive ORA
trials (1.74%) in 32 groups of 1,000. The paired columns in these compact
results refer to the archived GSEA analysis. Compare matched Reactome
workflows with their input selection and eligibility rules in view.

`trial_results.tsv` and `trial_results.rds` contain the trial IDs, statuses,
counts and paired diagnostics. `summary.tsv` contains stratum-level rates.
The separately adjusted common-family sensitivity has its own trial and
summary tables. Raw checkpoint records remain under `data/large/`.
See [the analysis specification](../reactome_ora_methods.md) for test
families, the input gate and the interpretation of conditional intervals.

Figure 3 shows GO biological processes in panel A, Reactome in panel B, and
WikiPathways in panel C. The three panels share a color scale and a 0–100%
rate axis, with a dashed reference line at 5%.
Render the figure from the bundled summaries with
`pixi run Rscript scripts/render_figures.R --only=ora`.
