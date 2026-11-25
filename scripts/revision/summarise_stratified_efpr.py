"""Write descriptive eFPR tables from the validated 32,000-trial summaries."""
import csv
from collections import defaultdict
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path


def read_table(path):
    with open(path, newline="") as table:
        return list(csv.DictReader(table, delimiter="\t"))


def percentage(positive, trials, digits):
    rate = Decimal(positive) * 100 / trials
    return rate.quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP)


def write_table(path, headers, rows):
    with path.open("w", newline="") as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


models = ("TK6", "MCF-7", "HepaRG", "Liver spheroids")
endpoints = (
    "DEG adjusted p only", "DEG adjusted p and fold change", "Final DEGs",
    "ORA GO BP", "ORA Reactome", "ORA WikiPathways",
    "GSEA Figure 4A", "GSEA Figure 4B", "GSEA direct minSize 15",
)
grid = defaultdict(dict)


def add_count(row, size_column, endpoint, positive, total, rate):
    model, n = row["cell_line"], int(row[size_column])
    positive, total = int(positive), int(total)
    assert model in models and 3 <= n <= 10 and 0 <= positive <= total == 1000
    assert endpoint not in grid[model, n]
    assert Decimal(rate) == Decimal(positive) / total
    grid[model, n][endpoint] = positive


deg_stages = dict(zip(("P-value filter only", "P-value + fold-change", "Final DEGs"), endpoints[:3]))
for row in read_table("data/derived/deg_figure_rates.tsv"):
    add_count(row, "n_per_group", deg_stages[row["Filter"]], row["trial_count"],
              row["n_trials_observed"], row["false_positive_rate"])

ora_collections = {"GO": "ORA GO BP", "Reactome": "ORA Reactome", "WikiPathways": "ORA WikiPathways"}
for row in read_table("docs/revision/reactome_ora_v2/ora_figure_rates.tsv"):
    add_count(row, "n", ora_collections[row["collection"]], row["n_positive"],
              row["n_trials"], row["rate"])

for row in read_table("docs/revision/gsea_unique_ids_v1/summary.tsv"):
    for column, endpoint in (("gsea_rate", "GSEA Figure 4A"), ("gsea_NES2_rate", "GSEA Figure 4B")):
        positive = Decimal(row[column]) * int(row["n_trials"])
        add_count(row, "n", endpoint, positive, row["n_trials"], row[column])

direct_path = "data/large/gsea_unique_ids_minsize15_v1/full/summary.tsv"
for row in read_table(direct_path):
    if row["cell_line"] == "All":
        continue
    stratum = grid[row["cell_line"], int(row["n"])]
    assert int(row["posthoc_positive"]) == stratum["GSEA Figure 4A"]
    add_count(row, "n", "GSEA direct minSize 15", row["direct_positive"],
              row["n_trials"], row["direct_rate"])

assert set(grid) == {(model, n) for model in models for n in range(3, 11)}
assert all(set(counts) == set(endpoints) for counts in grid.values())
out = Path("docs/revision/efpr_by_stratum")
headers = [endpoint + " (%)" for endpoint in endpoints]
by_size = (
    (model, n, 1000, *(percentage(grid[model, n][endpoint], 1000, 1) for endpoint in endpoints))
    for model in models for n in range(3, 11)
)
write_table(out / "efpr_by_cell_model_and_n.tsv",
            ["cell_model", "n_per_group", "trials_per_endpoint", *headers], by_size)
by_model = (
    (model, 8000, *(percentage(sum(grid[model, n][endpoint] for n in range(3, 11)), 8000, 2)
                    for endpoint in endpoints))
    for model in models
)
write_table(out / "efpr_by_cell_model.tsv", ["cell_model", "trials_per_endpoint", *headers], by_model)
