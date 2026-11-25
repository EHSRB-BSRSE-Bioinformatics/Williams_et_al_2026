# Signature-specific null results

All 32,000 frozen trials are present: four cell models × eight sample sizes
(n=3–10 per group) × 1,000 trials. Exact sample/group assignments match both
archived DEG results and archived classifier inputs. No fits or classifier
training were performed. Final DEG identities match the archives in every trial.

All 64 published genes map to 72 measured features. Six symbols have multiple
probes (TIGAR, CBLB, IKBIP, MDM2, PHLDA3, TP53I3); no measured feature maps to two
signature genes. Ten historical symbols require the explicit reference renamings.
The probe reference verifies 63 genes/71 probes. USP41_33775 matches uniquely by
symbol but lacks a reference probe ID; excluding it changes none of the reported
experiment-positive rates, at any stage.

The following ranges cover n=3–10; availability ranges additionally cover all
individual trials, while rate ranges compare the eight 1,000-trial strata.

| Cell model | Pooled detected /64 | Detected per trial | Initial per trial | Relevant/tested per trial | Nonmissing BH per trial | Final signature trial rate | Whole final trial rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| TK6 | 63 | 58–63 | 58–62 | 53–59 | 6–59 | 0.9–3.3% | 6.3–14.9% |
| MCF-7 | 64 | 59–64 | 58–61 | 52–60 | 4–60 | 0–4.2% | 0.9–15.4% |
| HepaRG | 61 | 56–61 | 55–60 | 53–57 | 2–57 | 0–0.1% | 0.9–8.7% |
| Liver spheroids | 63 | 54–63 | 54–63 | 44–56 | 2–56 | 0–0.9% | 3.3–24.4% |

USP41 has no counts in TK6, HepaRG or spheroids. HIST1H2BI/H2BC10 and LCE1E also
have no counts in HepaRG. These genes remain explicit in all mapping and frequency
outputs. Low-count detection is distinct from initial retention and testing.
`summary.tsv` reports all eleven stages by cell and n, including significance,
fold change, quantile, spike, their intersection and final calls; the detailed
QC table reports every gene, including zero-retention entries.

Using BH significance alone, 1,081/32,000 trials (3.38%) contain a signature gene:
2.9–4.9% in TK6, 3.0–9.6% in MCF-7, 0.7–2.6% in HepaRG, and 1.4–4.8% in
spheroids. These BH values come from the archived whole-transcriptome feature
tests; they were not recalculated within the 64-gene subset.

With all final rules, 245/32,000 trials (0.765625%) contain a signature DEG;
2,202/32,000 (6.88125%) contain a whole-transcriptome DEG. Of the signature
trials, 140 are TK6, 91 MCF-7, 1 HepaRG and 13 spheroids. The gene-count
distribution has 31,755 zero-call trials, 118 with one gene and 127 with two or
more genes; the maximum is 18 and every cell/n stratum has median zero.
The complete distribution by cell and n is in `distribution.tsv`.

Saved classifier component calls give the following final-rule pairing:

| | Any positive saved classifier component | No positive saved component |
|---|---:|---:|
| At least one final signature DEG | 1 | 244 |
| No final signature DEG | 16 | 31,739 |

Probability >0.9 is positive in 17 trials (one TK6, 16 spheroids); PC1 <0 in four
spheroid trials, all included in those probability positives; HCA is never DDI.
Raw call files omit HCA, so saved HCA calls are joined from the frozen classifier
summary, with source fields exposed. Each component and the descriptive union
have separate pairing tables. Individual signature DEGs and classifier positivity
are different endpoints.

These are trial proportions conditional on these samples and frozen overlapping
Monte Carlo draws. They are neither proportions of the 64 genes nor evidence of
general BH FDR control. No independence-based uncertainty intervals are reported.
Per-gene frequencies are descriptive fractions of the same full 1,000-trial
strata; stage availability is supplied alongside them.

## Publication-summary comparison

Final whole-transcriptome counts agree with the compact publication summary for
all 32,000 trials. Literal comparisons to the compact labels `P-Value Only` and
`P-Value + Fold Change` differ in 2,151 and 622 trials, respectively. The focused
analysis derives those stages directly from archived adjusted p-values and
shrunken fold changes, as specified in the methods; it does not infer their
meaning from legacy labels or use `DESeqResFiltered`. The diagnostic audit lists
all affected IDs and both counts. Publication outputs remain unchanged; no cause
is inferred from these label-level discrepancies alone.
