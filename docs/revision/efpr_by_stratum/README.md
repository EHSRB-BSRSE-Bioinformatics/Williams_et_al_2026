# Experiment-level false-positive rates by cell model and sample size

[Cell-model percentages](efpr_by_cell_model.tsv) pool the eight equally
represented sample sizes within each model, giving 8,000 mock trials per
endpoint. [Cell-model and sample-size percentages](efpr_by_cell_model_and_n.tsv)
use 1,000 trials per endpoint and report the sample size per group. A positive
trial has at least one qualifying feature or pathway. All percentages use the
full trial denominator for the corresponding group.

The tables use the [DEG screening and final-call rates](../../../data/derived/deg_figure_rates.tsv),
the [three Figure 3 ORA rates](../reactome_ora_v2/ora_figure_rates.tsv),
and the [unique-ID GSEA results](../gsea_unique_ids_v1/summary.tsv). The direct
`minSize = 15` column uses
`data/large/gsea_unique_ids_minsize15_v1/full/summary.tsv` from the
[completed sensitivity workflow](../gsea_unique_ids_v1/REPRODUCE.md). Figure 4B
uses BH-adjusted p < 0.05 and |NES| > 2.

[ORA universe and threshold diagnostics](ora_universes_and_thresholds.tsv)
reports median effective gene-universe sizes over the same trials with saved
GO and WikiPathways universes. The eFPR columns use all 1,000 trials in each
stratum. The WikiPathways cutoff checks reapply the stricter raw and
BH-adjusted p-value cutoff of 0.01 to the saved WikiPathways results, first
with the original q-value cutoff of 0.20 and then with GO's cutoff of 0.05.
The gene universe, pathway memberships, and existing BH adjustments are
held fixed. The audit checks that the original export rules reconstruct all
published GO and WikiPathways counts.

From the repository root, with the completed direct GSEA run and the original
ORA archives in the sibling `2023_Meier_FDR` directory:

```bash
python scripts/revision/summarise_stratified_efpr.py
pixi run Rscript scripts/revision/audit_ora_universes.R
```
