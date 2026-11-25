# Validation record

Completed on branch `analysis/tgx-signature-null-behavior` on 2026-11-17.
The extraction used the locked `deg-study` R environment only to deserialize
historical S4 results; no DESeq2 model or classifier was fitted.

The unit suite passes strict adjusted-p-value and fold-change boundary checks,
missing-p-value handling, duplicate/missing trial-ID rejection, unequal-group ID
rejection, multiple-probe deduplication, exclusion-set reconstruction, and
independence from `DESeqResFiltered`. All 17 portable archived reference trials
match the reconstructed final calls. Optimized set membership agrees with direct
feature-ID matching for every reference stage.

The complete extraction verifies every trial ID and exact sample/group assignments
against frozen metadata, rejects duplicate feature IDs, checks all relevant
archive features against the count matrix, and checks sample depth and initial
signature availability. All 32,000 reconstructed final feature sets exactly match
archived `DESeqResDEGs`. All four archive SHA-256 values agree with the recorded
source manifest. Classifier input assignments match all 32,000 frozen trials;
each contains the same 64 measured gene symbols. Raw probability and PC1 values
agree with the frozen compact classifier summary. HCA is explicitly sourced from
that summary where absent in the raw calls.

The output suite passes denominator completeness (1,000 trials in each of 32
strata, complete iteration IDs and group sizes), all stage membership inclusion
checks, gene-count distribution totals, all classifier pairing marginal totals,
whole-transcriptome versus signature rate bounds, and strict-probe sensitivity
bounds. Final whole-transcriptome counts match the compact publication summary
for all 32,000 trials. Intermediate label comparisons are retained as diagnostics:
2,151 `P-Value Only` counts and 622 `P-Value + Fold Change` counts differ from the
direct archived-rule selections. The new stage results do not use those legacy
counts. Exact IDs and counts are in `publication_summary_mismatches.tsv.gz`.

The supplementary figure was rendered from the compact analysis RDS to
`figures/supplementary/Fig_S13.pdf` and the 600-dpi
`figures/supplementary/Fig_S13.png`; both outputs were visually
checked. `figure_provenance.tsv` records the input, renderer, and output hashes.
Results, mapping, frequencies, all 32,000 trial records and per-gene Boolean memberships occupy
approximately 6 MB, including an approximately 0.8 MB RDS and a 1.1 MB compressed
wide trial TSV. Original publication inputs and outputs were not modified.
Existing unrelated workspace changes were left untouched.
