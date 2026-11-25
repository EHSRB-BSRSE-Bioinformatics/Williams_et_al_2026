# Selected archived DEG results

These 17 frozen trials are small regression references extracted from the
original per-cell-model DEG archives, without refitting. All feature rows and
numeric columns in each selected result table are retained; S4 results were
converted to ordinary data frames. The reference RDS is approximately 7.13 MiB.
The input gene matrix is not reduced for validation.

- `deg_reference_trials.tsv`: trial IDs, deterministic selection reasons, source
  archive filenames, and original tested/significant/final feature counts.
- `deg_reference.rds`: main result tables, filter summaries, excluded-feature
  tables, retained sample metadata, and retained sample counts.
- `deg_reference_metadata.rds`: original results metadata and ashr prior metadata,
  including the recorded ashr version (2.2.54).
- `deg_source_archives.tsv`: SHA-256 checksums of the complete source archives.

See `docs/deseq2_methods.md` for extraction and validation commands and
`docs/deg_validation.md` for measured outcomes. The reference set is deliberately
chosen to exercise boundary sample sizes, nonzero DEG calls, and quantile/spike
filter rejections; it is not a statistical sample of the 32,000 trials.

The original complete R session was not stored in these result tables. The
historical DESeq2, edgeR, and R versions cannot be recovered from this reference
metadata alone. The source environment file separately declares DESeq2 1.38.0
and edgeR 3.40.0; `deg-study` reconstructs those declared versions. The locked optional environment records the software used for
the current validation; it is not presented as the historical environment.
