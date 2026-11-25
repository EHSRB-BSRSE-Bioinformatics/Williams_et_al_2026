# Historical supplementary figures S7 and S9–S12

## Scope and numbering

This standalone port reconstructs historical plotting semantics, not a new
analysis of the revised endpoints. **S8 remains the existing baseline
mean–variance figure.** The historical variance-versus-false-positive screenshot
sometimes called S8 is rendered as **S9** here.

Implementation: [R/legacy_supplementary_figures.R](../R/legacy_supplementary_figures.R).
The figure registry provides the output filenames. Outputs are PNG (300 dpi)
and Cairo PDF, with white backgrounds, under `figures/supplementary/`.

## Public reconstruction

Run from the repository root in the existing default pixi environment:

```sh
pixi run --locked legacy-data
pixi run --locked legacy-figures
pixi run --locked legacy-check
```

Reconstruction uses only the compact references in
`data/reference/legacy_supplementary`, the public count table, QC-applied sample
metadata, and frozen Monte Carlo trial metadata. It does not consult a sibling
repository, large caches, Reactome packages, or an online service. DESeq2, edgeR
and enrichment computation are not needed for this public step. `sha256sum` must
be available. Plotting uses the existing ggplot2 dependency.

## Re-extraction from scientific sources

```sh
pixi run Rscript --vanilla scripts/rebuild_legacy_supplementary_data.R \
  --extract-reference --source-dir /path/to/historical-analysis-or-archive-root
```

The directory argument accepts both `--source-dir DIR` and `--source-dir=DIR`.
The archive root must contain these historical relative paths:

- `bootstrap_DEGs/output/1000_iterations/nobatch_2023-05-10/{tk6,mcf7,heparg,spheroids}.RDS`
- `bootstrap_GSEA/2023-11-08/{tk6,mcf7,heparg,spheroids}_gsea.RDS`
- optionally, the two historical processing notebooks under `Rmd`, to hash the
  exact plot definitions.
- the May 2023 `ddsList$all` normalized-count export under
  `analysis.cell_lines_comparison/DEG_lists` (the 2023.10.12 timestamped file).

Without a directory argument extraction looks in the current large-data directory.
The current installation has DEG recurrence caches there, but raw historical GSEA
files are in the archives, not revised GSEA outputs. Do not substitute the revised
size-filtered GSEA family cache or unique-gene-ID results for these raw archives.

For DEG extraction the existing current recurrence caches are used when available,
checking model, all 8,000 trial IDs, FINAL counts, and archive SHA-256 if the archive
is supplied. If a cache is unavailable, its raw DEG archive is read directly.
Historical compiled FINAL counts are exported and compared with cached FINAL
probe-list lengths. Missing `Filter` labels in the compiled table are explicitly
excluded; they must not create NA pseudo-trials.

Two already-extracted current provenance resources are also required for
re-extraction: the verified S5/S6 archived size-factor export and the historical
workspace pathway-name table under `data/large/reactome_ora_v1/workspace_provenance`.
These are frozen scientific inputs, not regenerated using a contemporary database.
The compact port exports its own copies, so neither is needed for public rebuild.

The full public count table and trial caches contain 22,533 probes, but the
historical `ddsList$all` expression object contains **22,460**. Extraction takes
that exact feature universe from the archived normalized-count export and checks
all 301 DMSO columns against raw counts divided by the frozen size factors,
allowing only the text export's five-decimal rounding. S7/S9 therefore contain
220,809 historical join rows, not an extra 292 zero rows from the broader trial
universe. The full screening reference still retains every recurrent probe.

To independently refresh only this expression reference without reading the
GSEA archives again, use `--extract-expression-reference --source-dir DIR`.
This mode also rebuilds data and refreshes hashes; it requires the remaining
reference files to have already been extracted.

### Sufficient reference exports

| Reference | Purpose |
|---|---|
| `screening_counts.rds` | All positive padj-only probe recurrence counts by model and n; not only selected probes |
| `expression_reference.rds` | Archived feature universe and DMSO `ddsList$all` size factors, preserving numeric precision |
| `final_trial_counts.rds` | Historical FINAL DEG count for each of the 32,000 trial IDs |
| `gsea_counts.rds` | All positive raw historical GSEA recurrence counts by pathway ID/model/n, sufficient to rerun selection |
| `pathway_names.tsv` | Frozen historical Reactome names for all recurrent pathway IDs |
| `source_provenance.tsv` | SHA-256 of public inputs, references, caches, archives, normalization/name exports and historical notebooks |

Counts and metadata are not duplicated in the reference directory. No individual
trial DEG/GSEA result matrices, full expression matrices, saved analysis functions,
local absolute paths, or credentials are exported.

## Exact definitions and historical source anchors

The DEG processing notebook `Rmd/Process bootstrapped results.Rmd` defines the
padj-only recurrence in its “Never DEGs, always DEGs” section (approximately
lines 809–885), expression transformation and join around 926–1034, S9 variance
scatter around 1036–1050, S7 mean-expression scatter around 1137–1150, and S10/S11
sample participation around 1246–1314. The GSEA notebook
`Rmd/Process bootstrapped results GSEA.Rmd`, approximately lines 180–330, supplies
S12 selection and plotting. Source hashes, rather than line numbers alone, pin
the exact versions. The source notebooks are read-only.

### S7: expression versus false-positive count

- Unit of expression is an archived **probe**, retaining the full identifier,
  not a symbol-level aggregated gene.
- For each model, compute mean of `log2(raw_count / archived_size_factor + 1)`
  over the QC-applied DMSO wells. These are DESeq2 normalized counts, not CPM,
  VST, or log of the mean count.
- Recurrence is non-missing archived `DESeqRes$padj < 0.05`, counted separately
  in each n=3,…,10 block of 1,000 trials. No FINAL filter is applied.
- The original join drops n before joining on probe/model. Thus a recurrent
  probe contributes one point for **each n with a positive recurrence count**;
  a probe never recurrent at any n contributes one zero-count point. It does
  **not** sum eight sample sizes or make an eightfold zero-completed grid.
  The exported scatter retains n (NA only for the singleton zero row) for audit.
- All points are included, with a blue within-model ordinary least-squares line
  and its default 95% confidence band, as in the original. These descriptive
  lines are not independent-trial causal evidence or corrected inferential tests.

### S9: expression variance versus false-positive count

Same probes, recurrence table, historical join, four model panels and linear
fit as S7. The x axis is the **sample variance** of log2 normalized counts within
each model's DMSO wells (denominator number of wells minus one), not variance of
raw counts, CPM variance, standard deviation, or the existing S8 mean–variance
data. S7/S9 use common axis scales across model panels.

### S10: sample participation histogram

For each sample/model/n, the numerator is the number of participating trials
whose historical archived **FINAL** DEG count exceeds zero. The denominator is
all trials in that stratum containing that sample, in either randomized group.
A trial contains distinct samples and contributes once to each participating
sample. It is not divided by all 1,000 stratum trials and does not count multiple
DEGs multiple times. Historical `n_final_DEGs` is the length of archived
`DESeqResDEGs`, mirrored by the current compiled `Final DEGs` label.

The histogram has 30 bins and four model columns. Just as in the source, each
sample/n proportion is an observation: a sample used at all eight n values
appears eight times, **not** as a single pooled proportion. The y label makes
this historical multiplicity explicit. These FINAL trial events must not be
described as the padj-only probe endpoint used in S7/S9.

### S11: plate position heatmap

The original `geom_tile` maps plate row (A–H) to x, plate column to y, model to
facet, and the S10 participation proportion to fill. It does not facet or
aggregate by n, so eight values per sample can occupy the same coordinate.
The preceding descending sort on the numerator determines the last visible
tile. This port selects the **last-drawn row** deterministically using the
historical stable sample/n/model key order followed by descending numerator
order; the exported table retains selected n, numerator, denominator, sample,
batch and number of overplotted rows.

This is a faithful historical rendering, **not** a pooled estimate of a plate
effect. Across different sample sizes the denominators differ, so selecting a
visible tile by numerator is not selecting the smallest proportion. Missing
coordinates remain blank. Like the historical plot, the shared color scale spans
the observed visible proportions, rather than washing out differences by forcing
0–1 limits; the values themselves remain raw proportions in [0,1].
A future pooled or n-faceted plate analysis should be separately named and
specified rather than silently replacing this endpoint.

### S12: top GSEA pathway recurrence

- Read raw archived November 2023 GSEA results; count pathway calls with
  non-missing `padj < 0.05` by model and n. Do not impose an additional NES
  threshold, restrict sizes after testing, recompute BH, or rerun ranking with
  revised unique gene IDs. “Unfiltered” means **no revised post-hoc filters**,
  not removal of whatever settings were used to produce the original archives.
- The source sorts on absolute NES but does not threshold it for this figure.
  Sorting does not change counts. The heading in intermediate source plots
  says “top 10”, but the executable selection is **top 20 per model/n with
  boundary ties**. The union of selected **pathway names** is then retained
  across all observed model/n combinations.
- Use the saved historical workspace human Reactome names, strip the
  `Homo sapiens: ` prefix, and order names by their mean observed positive count
  as in `reorder(path_name, count)`. Duplicate IDs sharing a selected name
  cause a deliberate error instead of silent aggregation.
- Four model columns and n=3,…,10 on x; height adapts to the selected list for
  legible pathway labels. Grey denotes **absence from the positive-count
  recurrence table**, not established biological absence or testing eligibility.
  Positive tiles use counts, not pooled counts or rates.

## Provenance, validation and limitations

The data and figure provenance tables under `data/derived/legacy_supplementary*`
hash public inputs, code, reference files and outputs. The source manifest
uses archive-root-relative paths; no absolute machine paths are saved. Public
rebuild checks the source-manifest checksums of every reference and public input
before computing data. Dataset definitions are also exported as a readable TSV.

Tests cover independent normalization/variance spot checks, literal historical
join multiplicities, real frozen-trial sample denominators and FINAL events,
last-drawn plate coordinates, top-20 boundary ties, no silent n pooling,
reference/output hashes and construction of all five ggplots. The summary's
schema and endpoint definitions are validated before rendering. Repeating public
rebuild should reproduce data hashes; PDF byte identity is not promised across
different Cairo/font/toolchain installations.

The executed extraction and independent checks produced 220,809 S7/S9 scatter
rows, 2,408 sample/n rows from 301 DMSO samples, and 57 selected GSEA pathway
names. The compact reference exports occupy approximately 676 KiB. All five
PNG/PDF pairs were rendered and checked via reduced-size previews. A repeated
public rebuild reproduced the SHA-256 hashes of the summary and every scientific
TSV output. Tests passed before and after rendering. The local Cairo/font stack
emitted a non-fatal Fontconfig initialization warning; figure files and tests
were nevertheless successfully generated.

These are computationally reproduced historical definitions, not pixel-identical
images. Typography, point size, facet order, white backgrounds and explanatory
labels are modernized. Archived GSEA statistics are extracted, not rerun: random
seeds, duplicate-gene-ID ranking problems and revised enrichment validity are not
resolved by this figure port. The old scatter weighting and plate overplot are
explicitly preserved for comparison, not endorsed as optimal analysis choices.

The full `pixi run --locked figures` workflow renders S7 and S9–S12 directly
from their tracked compact summary.
