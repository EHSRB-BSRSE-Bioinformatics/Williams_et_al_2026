# Unique-ID GSEA validation

Unique-ID GSEA was positive in 31,078/32,000 trials (97.118750%), compared with 30,661/32,000 (95.815625%) in the archived duplicate-ID analysis and 557/32,000 (1.740625%) for Reactome ORA. 598 archived positive calls became negative and 1,015 archived negative calls became positive.

The high experiment-level positive frequency persists in every cell-model/sample-size group after requiring unique Entrez ranking identifiers. Primary group rates range from 92.9–99.1%; the additional |NES|>2 sensitivity
is positive in 30,079/32,000 (93.996875%), with group rates 89.1–98.4%.
The reusable `gsea_summary.rds` supplies the corrected figure inputs.

## Coverage and reconstruction

All 32,000 exact frozen trial IDs are accounted for, with 1,000 per group, zero
unresolved execution failures and 0 explicitly untestable trials.
Every evaluated ranking has unique, nonmissing Entrez names and finite scores.
The full rerun exactly reproduced every selected pilot raw table with its
recorded trial seed, including reference-derived and archive-derived rankings.
Every raw result was checked by rebuilding size-restricted BH and both displayed
count endpoints, including an independent reconstruction through the public
summary builder. Exact common-family sizes agree with ORA in every trial.
The full runner checks hashes and IDs before reusing checkpoints.

All 17 reference historical vectors exactly reconstruct the literal historical
mapping procedure. Every unique mapped ID is preserved. Effective pathway sizes
are unchanged in the references. Repeated Entrez/score pairs range from
1,494 to 1,848;
352–545
Entrez IDs per reference require selection among distinct scores. The maximum
absolute-score rule can prefer extreme estimates; this is a declared ranking
choice, not a claim that all aliases estimate an identical quantity.

Across the full run, selection among distinct scores was needed in
32,000/32,000 trials, totaling
15,228,499 trial–Entrez instances.
Unique ranking lengths range from 10,142
to 13,588. Remaining nonzero tied-score
occurrences range from 437
to 611; these are kept in
deterministic input order, with no jitter. Every full-run ranking's identifier
set also exactly matches the independently prepared ORA background.

## Pilot and undefined tests

The 33-trial same-environment pilot changed positive/negative status in one trial between duplicate and unique rankings. All repeated unique and duplicate raw tables were identical; the serial/fork comparisons also agreed. Significant pathway counts changed in 31 pilot trials, including trials that remained positive. The same-environment duplicate control differed from the archived positive call in any pilot trial. The historical archive has unrecoverable seeds and no complete runtime record, so its differences from the new results are not a pure estimate of a deduplication effect.

Among size-eligible tests, 6,012/34,669,235 p-values are undefined,
affecting 503 trials. These tests are retained in raw outputs,
excluded before BH, and reported through `n_eligible`, `n_undefined` and
`n_tested`. Thus undefined tests reduce the BH family rather than being counted
as negative pathway tests. The primary nPermSimple=500 setting is preserved;
the fgsea unbalanced-score warnings and other warnings are retained per trial.
Undefined NES associated with undefined p-values is excluded with that test.
Missing log2err alone does not exclude a finite p-value. Computational errors
would block the final summary, and are not treated as zero-positive trials.

## Paired Reactome comparison

The common-defined-eligible-family sensitivity is positive for unique-ID GSEA
in 31,078/32,000 (97.118750%) and for ORA in
557/32,000 (1.740625%). Both methods' BH values are
recalculated from raw p-values within each trial's intersection, rather than
subsetting previously adjusted discoveries. `summary.tsv` includes paired ORA
minus GSEA differences and 95% conditional paired-record bootstrap intervals,
four-way concordance, count distributions and ranking/family diagnostics.

Trials reuse underlying biological samples; these intervals describe conditional
Monte Carlo record variation and do not represent independent biological
replication. ORA retains its sparse-input selection and insufficient-input
gate, whereas GSEA uses all mapped fold changes. The contrast cannot be
attributed solely to the enrichment algorithms.

| Cell model | n/group | Unique GSEA | Archived GSEA | ORA | Lost positives | Gained positives |
|---|---:|---:|---:|---:|---:|---:|
| HepaRG | 3 | 0.991 | 0.982 | 0.024 | 5 | 14 |
| HepaRG | 4 | 0.988 | 0.970 | 0.006 | 9 | 27 |
| HepaRG | 5 | 0.985 | 0.979 | 0.003 | 12 | 18 |
| HepaRG | 6 | 0.982 | 0.965 | 0.001 | 13 | 30 |
| HepaRG | 7 | 0.980 | 0.974 | 0.001 | 20 | 26 |
| HepaRG | 8 | 0.977 | 0.960 | 0.000 | 18 | 35 |
| HepaRG | 9 | 0.975 | 0.962 | 0.001 | 18 | 31 |
| HepaRG | 10 | 0.972 | 0.962 | 0.000 | 25 | 35 |
| Liver spheroids | 3 | 0.929 | 0.923 | 0.042 | 35 | 41 |
| Liver spheroids | 4 | 0.936 | 0.939 | 0.027 | 29 | 26 |
| Liver spheroids | 5 | 0.954 | 0.947 | 0.019 | 18 | 25 |
| Liver spheroids | 6 | 0.952 | 0.956 | 0.022 | 29 | 25 |
| Liver spheroids | 7 | 0.962 | 0.955 | 0.006 | 18 | 25 |
| Liver spheroids | 8 | 0.961 | 0.969 | 0.006 | 27 | 19 |
| Liver spheroids | 9 | 0.958 | 0.972 | 0.006 | 32 | 18 |
| Liver spheroids | 10 | 0.959 | 0.965 | 0.003 | 27 | 21 |
| MCF-7 | 3 | 0.991 | 0.974 | 0.056 | 3 | 20 |
| MCF-7 | 4 | 0.989 | 0.976 | 0.022 | 10 | 23 |
| MCF-7 | 5 | 0.983 | 0.973 | 0.019 | 15 | 25 |
| MCF-7 | 6 | 0.980 | 0.971 | 0.008 | 16 | 25 |
| MCF-7 | 7 | 0.971 | 0.954 | 0.003 | 21 | 38 |
| MCF-7 | 8 | 0.984 | 0.969 | 0.003 | 13 | 28 |
| MCF-7 | 9 | 0.972 | 0.952 | 0.001 | 22 | 42 |
| MCF-7 | 10 | 0.984 | 0.951 | 0.002 | 13 | 46 |
| TK6 | 3 | 0.990 | 0.976 | 0.060 | 6 | 20 |
| TK6 | 4 | 0.978 | 0.971 | 0.040 | 19 | 26 |
| TK6 | 5 | 0.979 | 0.961 | 0.035 | 14 | 32 |
| TK6 | 6 | 0.975 | 0.944 | 0.026 | 13 | 44 |
| TK6 | 7 | 0.967 | 0.938 | 0.032 | 21 | 50 |
| TK6 | 8 | 0.966 | 0.938 | 0.023 | 23 | 51 |
| TK6 | 9 | 0.955 | 0.918 | 0.033 | 27 | 64 |
| TK6 | 10 | 0.953 | 0.915 | 0.027 | 27 | 65 |
