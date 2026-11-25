# Supplementary Figures S1 and S3: sample distances and gene examples

## Scientific status

S1 uses the frozen archived all-gene sample distances. S3 uses the four
frozen DMSO mock comparisons recovered from the gene-count figure embedded
in two historical manuscript drafts. Each group contains four wells.

## S1: transformed assay, population, and clustering

The legacy `Rmd/stats_report.Rmd`, chunk `sample-to-sample_heatmaps`, calculates
`stats::dist(t(assay(current_rld)))`, sorts the distance representation, and uses
the same sample distances for both pheatmap dendrograms. Its palette is reversed
RColorBrewer Blues. The all-gene tab calls this function on `rld`.
`scripts/load-facet-data.R` obtains that object from `rldList[['all']]` for the
all-sample comparison. The recovered workspace is identified by its relative
source path and SHA-256 in the reference provenance table.

The archived object is **not rlog despite its `rld` name**. Its producer,
`scripts/DESeq_functions.R::regularize_data`, calls DESeq2 VST with `blind=FALSE`
and supports optional limma batch-effect removal. The exact archived assay is
used without altering it. Configuration-dependent batch handling is not inferred
from the object's name. This port does not claim to have re-estimated the exact
original transformation from scratch; preserving the frozen assay distances
avoids changing dispersion fits, DESeq2 versions, or optional batch handling.

The frozen matrix has **22,460 probes and 365 QC-passing samples: 301 DMSO and
64 Media**, in four cell models. S1 includes all 365 samples, as does this archived
all-gene report; it is not a DMSO-only distance analysis. Extraction checks every
archived raw count against the public input and checks both cell-model and
chemical annotations. The recovered distance matrix has all 365 sample IDs.
The uncompressed version-2 RDS serialization of the full original-order assay
is SHA-256-hashed during extraction, but not shipped as a large redundant matrix.
The compact distance RDS itself is independently hashed.

Rendering uses **one global complete-linkage tree** computed from these distances,
not four independently clustered blocks and not distances between distance
profiles. Rows and columns share its order. Low distances are dark blue and high
distances pale/white; a four-color model annotation, global top/left trees, and
Euclidean-distance legend are included. Dense individual sample labels are
omitted; IDs remain in the compact reference. Base R replaces pheatmap, so branch
orientation and layout need not be pixel-identical to pheatmap/dendextend output.
No raw counts, CPM, log2 counts, or normalized-count proxy are substituted for VST.

## S3: exact panels, actual frozen samples, and units

Panels are ordered as four rows and two columns:

| Panels | Model | Left probe | Right probe |
|---|---|---|---|
| A, B | Liver spheroids | KRT7_10561 | FLNC_89124 |
| C, D | HepaRG | GSTA2_2817 | ABCB11_82408 |
| E, F | MCF-7 | H2AC17_33569 | H4C3_87348 |
| G, H | TK6 | FOXD4L3_33551 | H3-3A_92932 |

Probe suffixes, not historical histone-symbol spelling, are used to resolve
public input rows. Resolution must be unique; the requested modern symbol and
the actual input feature are both retained in the expression export. The current
public count table contains all eight requested spellings.

### Recovery evidence and deterministic selection

Two drafts in the read-only historical project archive contain byte-identical
709 × 911 PNGs of the eight requested panels. The image entries are
`word/media/image13.png` in `FDR draft v3.1.docx` and
`word/media/image11.png` in `FDR draft v2 aw.docx`, both under
`Paper 1 - FDRs`. The embedded PNG SHA-256 and source locations are recorded
in [s3_raster_sources.tsv](../data/reference/sample_supplementary/s3_raster_sources.tsv).
The original chat screenshot was unavailable for direct comparison.

The compact raster evidence records log10-axis tick values and their vertical
pixel centers plus the occupied blue (A) and green (B) pixel rows at each
facet's center. Pixel coordinates are zero-based. Within an 11-pixel window
centered on the group, a row is marked for B when green exceeds red by 100 and
blue by 80, and for A when blue exceeds red by 100 and green exceeds red by
30. This threshold excludes the background and panel grid. Fitting pixel row
against log10 tick count calibrates each panel independently. The PNG windows
below use an inclusive start row and an exclusive end row.

| Panel | A center column | B center column | Rows |
|---|---:|---:|---|
| A | 187 | 133 | 63–188 |
| B | 546 | 492 | 63–188 |
| C | 187 | 133 | 300–424 |
| D | 545 | 492 | 300–424 |
| E | 188 | 138 | 537–662 |
| F | 548 | 495 | 537–660 |
| G | 187 | 133 | 774–900 |
| H | 546 | 492 | 774–900 |

The rebuild projects each positive raw count to a pixel row, expands that
center by three pixels to represent the dot, and takes the union of occupied
rows within each group. It compares those rows with the observed rows using
the size of the symmetric difference. Overlapping dots can share pixel rows.
Scores are summed across A/B and both probes of the same cell model. The
search covers all 8,000 frozen trials per model, across n=3–10, with original
group labels. Each winner uses four wells per group and fits both probes.
The script checks for a clear gap to the runner-up and validates the exact
frozen sample IDs against the selection table.

| Model | Matched trial | Best score | Next score |
|---|---|---:|---:|
| Liver spheroids | `Liver spheroids_4_4_712` | 6 | 65 |
| HepaRG | `HepaRG_4_4_582` | 6 | 61 |
| MCF-7 | `MCF-7_4_4_59` | 9 | 60 |
| TK6 | `TK6_4_4_136` | 4 | 39 |

[Raster evidence](../data/reference/sample_supplementary/s3_raster_evidence.tsv)
· [Match scores](../data/reference/sample_supplementary/s3_selection_decisions.tsv)
· [Frozen wells and A/B assignments](../data/reference/sample_supplementary/s3_trial_selection.tsv)

Each panel shows all QC-passing DMSO wells for its model and repeats the four
selected wells in each of B and A. The historical red population was also
checked against the public counts. Its plotted pixel rows supported the
DMSO-only population, with several panels excluding Media counts that would
extend or add values to the observed distribution. In HepaRG, B has the higher
raw-count mean for both genes, matching the draft A/B labels.

The draft places raw counts on a log10 axis. This renderer defaults to raw
counts and uses a labelled `log10(count + 1)` display to retain zero counts.
Optional `--s3-normalized` uses raw count divided by the archived full-precision
S5/S6 DESeq2 size factor. `--s3-raw` is an alias for the default. Display
settings are recorded in render provenance. Horizontal offsets are
deterministic.

The matched A/B wells are historical illustrations of separation under mock
labels. The available draft raster and frozen counts do not provide archived
DESeq2 adjusted p-values for these n=4 trials. The renderer does not annotate
trial-level significance or imply that the panels pass the FINAL DEG filters.
The previous n=3 probe scan and its source hashes remain in the repository as
records of the earlier investigation.

## Files and public workflow

- [R/sample_supplementary_figures.R](../R/sample_supplementary_figures.R): validation,
  summary building, suffix resolution, base-R S1 and ggplot2/grid S3 rendering.
- [scripts/rebuild_sample_supplementary_data.R](../scripts/rebuild_sample_supplementary_data.R):
  public rebuild and opt-in frozen reference extraction.
- [scripts/render_sample_supplementary_figures.R](../scripts/render_sample_supplementary_figures.R):
  standalone S1/S3 PNG/PDF rendering.
- [tests/test_sample_supplementary_figures.R](../tests/test_sample_supplementary_figures.R):
  numerical/identifier/trial/group/population/normalization/provenance tests.

Run from the repository root:

```sh
pixi run --locked sample-data
pixi run --locked sample-figures
pixi run --locked sample-check
```

Routine rebuilding requires only shipped raw counts, metadata, the frozen trial
metadata, the S5/S6 size-factor reference, and the compact sample-supplementary
references. SHA-256 checks fail closed if any scientific input changes. Legacy
files are not read by this workflow. A second rebuild must produce the same
summary and count-table hashes.

The optional S1 archive reference refresh command is:

```sh
pixi run -e deg-study Rscript scripts/rebuild_sample_supplementary_data.R --extract-reference
```

The default read-only legacy root is the sibling `../2023_Meier_FDR`; an alternate
root can be supplied through `FDR_ORIGINAL_ANALYSIS`. Only relative source labels
and hashes are written into public provenance, never the machine's original
absolute location. Extraction needs DESeq2 and SummarizedExperiment; routine
rebuild/testing uses base R and rendering/testing uses existing ggplot2 and grid.
Linux `sha256sum` is required. No package installation or shared environment
change is needed; **pheatmap is not required** (not installed in the inspected
environments). Existing pixi default, deg, and deg-study environments suffice for
their respective jobs.

The compact references contain the frozen S1 distance matrix and annotation,
the S3 panel and raster evidence tables, and the recovered A/B selections.
Derived outputs include the rebuilt summary, S3 count table, and rebuild and
render provenance. The registered outputs are
`figures/supplementary/Fig_S1.{png,pdf}` and
`figures/supplementary/Fig_S3.{png,pdf}`. PNGs are 300 dpi.

The full `pixi run --locked figures` workflow renders both figures directly
from the tracked summary. PDF creation dates and graphics-device versions may
change artifact hashes across renders.
