# Reproducing Reactome ORA revision v1

This directory retains the manuscript-declared Reactome 1.77.0 sensitivity and
the annotation investigation. The **primary** archive-compatible 1.84.0
analysis is documented in `../reactome_ora_v2/REPRODUCE.md`.

Run from the publication repository root. Original DEG archives must remain at
`../2023_Meier_FDR/bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/`.
For the common-family sensitivity, GSEA archives must remain at
`../2023_Meier_FDR/bootstrap_GSEA/2023-11-08/`. Source paths are explicit in the
scripts and no archived results are written. Each primary checkpoint stores
code/annotation hashes, software session, elapsed time, IDs, complete prepared
input lists/backgrounds, and raw and adjusted pathway p-values.

## Primary reproduction using the frozen annotation

The checked-in frozen annotation avoids future database or download drift.
Install the locked environments if needed:

```bash
pixi install --locked -e deg-study --run-post-link-scripts
pixi install --locked
mkdir -p data/large/reactome_ora_v1/annotation
cp docs/revision/reactome_ora_v1/frozen_annotation_org3.17.rds data/large/reactome_ora_v1/annotation/
.pixi/envs/deg-study/bin/Rscript scripts/revision/check_reactome_ora.R
env OPENBLAS_NUM_THREADS=1 .pixi/envs/deg-study/bin/Rscript scripts/revision/run_reactome_ora.R
env REACTOME_ORA_CORES=8 OPENBLAS_NUM_THREADS=1 .pixi/envs/deg-study/bin/Rscript scripts/revision/run_reactome_ora.R --pilot-parallel
env REACTOME_ORA_CORES=8 OPENBLAS_NUM_THREADS=1 .pixi/envs/deg-study/bin/Rscript scripts/revision/run_reactome_ora.R --full
env REACTOME_ORA_ROOT=data/large/reactome_ora_v1 REACTOME_ORA_DOCS=docs/revision/reactome_ora_v1 .pixi/envs/default/bin/Rscript scripts/revision/summarise_reactome_ora.R
env REACTOME_ORA_ROOT=data/large/reactome_ora_v1 REACTOME_ORA_DOCS=docs/revision/reactome_ora_v1 .pixi/envs/default/bin/Rscript scripts/revision/check_reactome_summary.R
env REACTOME_ORA_DOCS=docs/revision/reactome_ora_v1 .pixi/envs/default/bin/Rscript scripts/revision/plot_reactome_ora.R
.pixi/envs/default/bin/Rscript scripts/revision/reactome_common_family.R
.pixi/envs/deg-study/bin/Rscript scripts/revision/verify_gsea_details.R
```

The full run resumes automatically. Completed checkpoints must have matching
IDs, code and annotation hashes and no failures; otherwise it stops. For a
changed analysis, create a new versioned output root instead of bypassing the
hash check. Atomic `.tmp` files do not count as completed checkpoints.
The direct Rscript paths select the existing locked Pixi environments without
changing the publication's dependency solves. The default environment is used
for summaries/plots; deg-study supplies R 4.2.3, Matrix and DESeq2 classes for
reading original S4 result tables. No model-fitting function is called.

`reactome_common_family.R` optionally accepts one or more archive stems
(`tk6`, `mcf7`, `heparg`, `spheroids`) for per-cell execution. The final invocation
without arguments assembles all four outputs. This is new restricted-family
postprocessing, not another reproduction audit of the submitted GSEA summary.

## Reconstructing annotation provenance and the pilot investigation

The release URLs and SHA-256 digests are recorded in the Methods and hash
files. To recreate the initial declared-environment annotation and then the
compatible candidate:

```bash
curl -L --fail https://bioconductor.org/packages/3.14/data/annotation/src/contrib/reactome.db_1.77.0.tar.gz -o data/large/reactome_ora_v1/annotation/reactome.db_1.77.0.tar.gz
curl -L --fail https://bioconductor.org/packages/3.16/data/annotation/src/contrib/org.Hs.eg.db_3.16.0.tar.gz -o data/large/reactome_ora_v1/annotation/org.Hs.eg.db_3.16.0.tar.gz
.pixi/envs/deg-study/bin/Rscript scripts/revision/freeze_reactome_annotation.R
curl -L --fail https://bioconductor.org/packages/3.17/data/annotation/src/contrib/org.Hs.eg.db_3.17.0.tar.gz -o data/large/reactome_ora_v1/annotation/org.Hs.eg.db_3.17.0.tar.gz
mkdir -p data/large/reactome_ora_v1/annotation/candidate_3.17
tar -xzf data/large/reactome_ora_v1/annotation/org.Hs.eg.db_3.17.0.tar.gz -C data/large/reactome_ora_v1/annotation/candidate_3.17
.pixi/envs/deg-study/bin/Rscript scripts/revision/check_annotation_candidate.R
```

The candidate check uses the retained initial pilot mapping TSV to reconstruct
the actually inspected saved gene sets, then compares candidate inputs.
For a fresh direct comparison with the archived S4 slots, run the pilot
validator. The execution here used cached data.table 1.14.6 from the original
Conda package cache, without installing it into or altering locked environments:

```bash
env OPENBLAS_NUM_THREADS=1 R_LIBS="$DATA_TABLE_R_LIBRARY" .pixi/envs/deg-study/bin/Rscript scripts/revision/validate_reactome_pilot.R --allow-known-annotation-drift
env OPENBLAS_NUM_THREADS=1 R_LIBS="$DATA_TABLE_R_LIBRARY" .pixi/envs/deg-study/bin/Rscript scripts/revision/validate_reactome_pilot.R --org3.17
```

Set `DATA_TABLE_R_LIBRARY` to an R-4.2-compatible library containing data.table 1.14.6. Production ORA, summary reconstruction
and its deterministic tests do not require data.table. Initial 3.16.0 pilot
discrepancies are preserved as evidence rather than erased. The final validator
checks the chosen 3.17.0 annotation directly against the archived slots.
