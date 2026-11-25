# Reactome ORA revision: analysis specification and provenance

Primary revision v2, 16 November 2026. Branch: `revision/reactome-ora-v1`.
Version 1 retains the analysis using the manuscript-declared Reactome 1.77.0;
version 2 uses the archive-compatible Reactome 1.84.0 collection.
This analysis is independent of the archived publication outputs.
No DESeq2 model or GSEA enrichment calculation is rerun here.

## Inputs and mapping

The input is the original per-cell-model DEG archive in
`../2023_Meier_FDR/bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/`:
`tk6.RDS`, `mcf7.RDS`, `heparg.RDS`, and `spheroids.RDS`.
These are the relevance-filtered, ashr-shrunken `DESeqRes` tables. The analysis
never reads `DESeqResFiltered` or `DESeqResDEGs`. It therefore does not depend
on the Figure 2A legacy indexing intermediate and does not apply the final
DEG quantile/spike rules.

`R/revision/reactome_ora.R::historical_ora_input()` removes suffixes beginning
at the first underscore and chooses the first maximum absolute shrunken log2
fold change per alias, ignoring missing estimates as `which.max()` does.
The selected probe supplies both fold change and adjusted p-value. Only then
is the alias mapped to all corresponding Entrez IDs. ENSEMBL-expanded repeated
annotation rows are retained in the frozen alias map. An Entrez ID is selected
if any mapped representative passes adjusted p < 0.05 and absolute signed
linear fold change > 1.5. Signed linear FC is `2^LFC` for positive LFC and
`-1/2^LFC` otherwise. The background contains all mapped IDs before selection,
including rows with missing adjusted p-values. Selected and background IDs
are subsequently made unique for hypergeometric testing.

There is no maximum-FC selection per Entrez ID before significance testing.
The existing public helper does perform that additional selection; it is
deliberately not reused or modified. Pilot comparisons against archived input
slots establish the achieved compatibility; annotation package declarations
alone are not treated as proof of original runtime identity.

## Frozen annotations

The **primary** Reactome collection is `reactome.db` **1.84.0**, with 2,596
human pathways, recorded source date 2023-04-05 and schema version 84. Human
Entrez memberships are read from its source SQLite `pathway2gene` table,
restricted to IDs beginning `R-HSA-`, the underlying collection exposed by
`reactomePATHID2EXTID` in the original GSEA source.

The manuscript-declared **1.77.0** release was evaluated first (v1; 2,520
human pathways; recorded source date 2021-10-01; schema 77). The common-family
diagnostic exposed effective-size mismatches with the archived GSEA. A bounded
comparison of all 17 archived reference trials across four cell models showed
that 1.84.0 exactly matched the raw pathway families, every effective size and
all leading-edge memberships. In the initial five TK6 references, 1.77.0
omitted 75–76 archived pathways and disagreed on 246–247 shared raw-family
sizes. Thus 1.84.0 is the closer
available historical match and supplies the primary v2 comparison. Exact
runtime version provenance still cannot be proved from size/membership
compatibility alone. The initial full v1 analysis is retained independently
as a database-version sensitivity, not substituted for v2.

The 9 November 2023 historical `.RData` workspace was subsequently inspected
read-only: all 2,596 human pathway IDs and complete gene memberships exactly
match the chosen 1.84.0 collection. Its entire 235,861-row alias/Entrez table
is identical to the frozen 3.17.0 table, including row order. Historical
annotation contents are therefore recovered directly, while original executed
package versions still lack a runtime record. Pathway-list ordering differs
between the SQLite split and saved list; set membership and ORA/BH results
are invariant to that ordering. No saved functions were evaluated.

The primary alias source is `org.Hs.eg.db` **3.17.0**, read from its source
package SQLite using AnnotationDbi **1.60.0** and the original `keys()`/`select()`
columns SYMBOL, ENTREZID, ALIAS, ENSEMBL. The original reports environment's
declared 3.16.0 release was tested first: it added Entrez 84215 to two TK6
reference lists relative to the archived slots. In 3.17.0, ZFP41 no longer maps
additionally to 84215. The 3.17.0 candidate matched all 20 available saved
GO/WikiPathways selected-gene slots (10 reference trials) and all 10 checked
WikiPathways effective universes. A subsequent direct validation also exactly
matched all 10 saved GO effective universes using the source SQLite BP-all
annotations. This motivated
the choice of 3.17.0; it does not establish an original runtime record or
all-trial mapping equivalence. The packages were retrieved from:

- <https://bioconductor.org/packages/3.14/data/annotation/src/contrib/reactome.db_1.77.0.tar.gz>
- <https://bioconductor.org/packages/3.17/data/annotation/src/contrib/reactome.db_1.84.0.tar.gz>
- <https://bioconductor.org/packages/3.16/data/annotation/src/contrib/org.Hs.eg.db_3.16.0.tar.gz>
- <https://bioconductor.org/packages/3.17/data/annotation/src/contrib/org.Hs.eg.db_3.17.0.tar.gz>

Downloads, initial annotation packages and candidate maps are under
`data/large/reactome_ora_v1/annotation/`. The primary v2 frozen map is
`data/large/reactome_ora_v2/annotation/frozen_annotation.rds`, also included in
`docs/revision/reactome_ora_v2/` for reproduction. Checksums and metadata are
retained with the revision. In v1, `frozen_annotation.rds` preserves the initial
3.16.0 pilot and `frozen_annotation_org3.17.rds` supplies its full 1.77.0 analysis.
The original GSEA runtime's package versions are
not recoverable from its result tables; the new analysis uses the declared
archive-compatible Reactome release, with compatibility diagnostics reported separately.

## Hypergeometric test and multiple testing

We implement the exact upper-tail hypergeometric test with R's `phyper()`.
For each trial, the effective universe is the supplied mapped background
intersected with the union of human Reactome IDs. The selected list is
intersected with that universe. For each pathway, let N be universe size,
n be selected-list size, K be pathway size within the universe, and k be
selected overlap. The raw p-value is
`phyper(k-1, K, N-K, n, lower.tail=FALSE)`.

The prespecified family comprises **all** pathways with effective K between
15 and 500 inclusive. Zero-overlap pathways remain in this family with p=1.
BH is calculated separately per trial across that complete family. A term
is significant only at BH-adjusted p < 0.05, strictly. There is no additional
raw-p cutoff and no q-value gate; q-values are neither estimated nor used.
This implements the requested hypergeometric ORA directly rather than relying
on `enricher()`'s conditional hit-containing family and exported-table gates.
Sparse-matrix multiplication accelerates size/overlap counting; direct
set-intersection and matrix implementations are checked for equivalence.

Trials with fewer than two unique selected Reactome-annotated genes or no
size-eligible pathways are classified as `insufficient_input`, contributing
zero positive results to the experiment-level denominator. The two-unique-gene
rule is explicit; it is not a reproduction of the historical pre-deduplication
vector-length execution gate. Other statuses are `significant`, `none`, and
`failed`. Failures are not counted as zero and block final summary validation.
Counts of selected IDs, mapped and effective backgrounds, annotated selected
IDs, tested pathways and significant terms are saved for every trial.

## Trial identity, execution and comparison

Frozen assignments in `data/derived/monte_carlo_trial_metadata.rds` define
32,000 unique IDs: four cell models, n=3–10 per group, 1,000 trials per group.
Names/order are checked against these blocks before attaching identifiers to
archived DEG results. This follows the documented historical positional
association; an undocumented historical reorder cannot be independently ruled out.
Execution proceeds one cell model at a time with atomic checkpoints per n.
Resume checks the IDs, annotation hash, code hashes and absence of failures.
Outputs never replace frozen summaries or figures.

V2 reuses the prepared Entrez lists/backgrounds saved from the original frozen
DESeqRes inputs by v1, records the source checkpoint digest, and recalculates
only Reactome hypergeometric tests with the compatible collection. It does
not refit DESeq2 or use regenerated differential-expression estimates.

The comparator is the unchanged compact GSEA summary underlying Figure 4A,
already verified in session 01 against every archived GSEA result. It uses
effective size 15–500, nonmissing raw p-values, and recalculated within-family
BH p < 0.05, without the optional |NES| > 2 filter. The Figure 4B counts are
retained as a descriptive additional column. Original GO and WikiPathways
counts are joined unchanged for context; their cutoffs and eligible families
differ and are not harmonized retrospectively.
GO used raw and BH-adjusted p≤0.01 with q≤0.05 when available; WikiPathways
used raw and BH-adjusted p≤0.05 with its recorded q cutoff 0.2 when available.
Both used the original effective-size defaults 10–500. Published counts retain
the historical exported-table/complete-row semantics documented in session 01.

Per-group results include ORA/GSEA positive fractions, ORA minus GSEA paired
differences, significant-count distributions, and both/ORA-only/GSEA-only/
neither counts and proportions. Conditional 95% percentile intervals use
10,000 ordinary paired-record bootstrap resamples (seed 20261116), implemented
equivalently by multinomial resampling the three-valued paired difference.
These quantify Monte Carlo record variation conditional on this data and
assignment pool. Reused biological samples mean trials are not independent
biological experiments; the intervals are not population-level biological
uncertainty estimates.

The common-family sensitivity intersects each trial's size-eligible ORA family
with the archived GSEA size-eligible, nonmissing-p-value family. It reuses raw
p-values and recalculates BH separately for **both** methods within exactly
that intersection, retaining the same mapped inputs and hypergeometric
background. This is distinct from simply counting pre-adjusted significant
terms after subsetting. Effective-size agreements and family overlap are
reported as diagnostics; no new GSEA statistics are generated.

Matching database, size limits and BH threshold reduces database/threshold
confounding. ORA uses a selected list and an annotation-intersected background;
GSEA uses the full historical ranking including duplicate IDs. Effective
eligibility and test families can therefore differ. The comparison does not
isolate the statistical algorithms perfectly, and pathway dependence also
precludes interpreting the 0.05 plot benchmark as guaranteed error control.
