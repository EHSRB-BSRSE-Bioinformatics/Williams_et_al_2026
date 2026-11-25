# Controlled unique-Entrez GSEA

The analysis uses frozen relevance-filtered, ashr-shrunken `DESeqRes` from
the original four cell-model archives, without refitting DESeq2 or selecting
significant genes. Trial identity and order come from the frozen Monte Carlo
metadata. Reactome memberships and list order come from the recovered saved
historical workspace; the alias map is the identical recovered mapping frozen
with the Reactome ORA revision. These contents match Reactome 1.84.0 and
org.Hs.eg.db 3.17.0. This is evidence of data identity, not a recovered original
software session.

We remove probe suffixes beginning at the first underscore and retain the
first probe having maximum absolute shrunken log2 fold change per alias, before
expanding every alias mapping row in its historical order. We remove identical
Entrez/score pairs, then retain the maximum absolute score for each Entrez ID.
Exact ties, including opposite-sign ties, retain the first historical mapping
occurrence. This second selection can favor extreme fold changes; its frequency
is reported separately from repeated annotation rows. Final rankings retain
every unique mapped Entrez ID and require nonmissing names and finite scores.
They are ordered ascending by score, with first historical mapping occurrence
breaking score ties; fgsea sorts descending internally. Fixed input order makes
remaining tie processing reproducible in the pinned runtime.
There is no random jitter. Remaining equal scores are retained and counted.

We run fgsea 1.26.0, built from the preserved source archive, with R 4.2.3.
The original manuscript declared fgsea 1.26.0; its executed historical runtime
has not been recovered. This new runtime is fully recorded. A dedicated local
package library overlays the locked `deg-study` environment; the explicit conda
artifact lock includes the base and added dependencies. The parent environment
is not changed. `pixi.lock`, `conda-explicit.txt`, and the checkpoint provenance identify this environment.

The primary call is `fgseaMultilevel(minSize=1, maxSize=500, nPermSimple=500)`.
Remaining algorithm defaults are made explicit: `sampleSize=101`, `eps=1e-50`,
`scoreType="std"`, `gseaParam=1`, `nproc=0`, and no deprecated `absEps` override.
`BPPARAM=SerialParam()` prevents nested parallelism. Data.table, OpenMP and
OpenBLAS each use one thread; trial-level fork workers provide parallelism.
Mersenne-Twister / Inversion / Rejection RNG is initialized separately for
every trial, with seed `20261116 + position_in_frozen_metadata` (1-based).
Seeds do not depend on worker scheduling. Original archive seeds are unknown.

All raw pathway results, including effective sizes, p-values, original package
adjustments, NES, ES, error estimates and leading edges are saved. For the
primary family, effective sizes 15–500 and nonmissing raw p-values are retained,
and BH is recalculated within that family for each trial. Significance is
strictly BH p<0.05. The Figure 4B sensitivity additionally requires |NES|>2.
Pathways with undefined p-values are excluded from the family and counted,
including separately the number among size-eligible pathways. fgsea's warning
identifies insufficient positive/negative preliminary permutations for some
pathways (`modeFraction<10`); these are scientific/numerical undefined tests,
not execution failures. We retain nPermSimple=500 for this controlled comparison
and do not replace these values by zero or one. Other missing precision estimates
(`log2err`) do not exclude a finite p-value. A trial with no defined eligible
tests is explicitly `untestable`; computational errors are `failed` and block
validation. No failed trial is counted as negative.

The pilot is selected before examining enrichment results: all 17 frozen
references plus minimum/maximum feature-count trials at n=3 and n=10 in every
cell model (ties follow frozen metadata order). Unique and historical duplicate
rankings use identical trial seeds and settings in the same environment. Each
is repeated; additional serial comparisons check scheduling independence.
This control separates a ranking change from the unknown software/RNG change
between the new runtime and unseeded archives. Small stochastic differences
near BH 0.05 cannot be assigned to deduplication using the archive comparison.

Full execution saves atomic checkpoints per sample size, one cell model at a
time, with code, annotation, input, environment and metadata hashes. Resume
rejects incompatible provenance and unresolved failures. The validation rebuilds
BH decisions and both Figure 4 counts from every raw result, checks all 32,000
unique frozen IDs and all 32 denominators of 1,000, and joins exact IDs to ORA
and archived GSEA. Common-family sensitivity intersects defined eligible
pathway IDs and recalculates BH from each method's raw p-values; the ORA
insufficient-input rule remains in force.

Paired ORA minus GSEA intervals use 10,000 paired-record bootstrap resamples
(multinomial equivalent; seed 20261116). Reused underlying biological samples
mean that these intervals describe conditional Monte Carlo variation and do
not represent independent biological replication. Matching pathways and BH
thresholds still leaves different inputs and selection gates for ORA and GSEA;
the experiment-level difference cannot be attributed solely to the algorithms.
