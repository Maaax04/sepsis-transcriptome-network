"""Performs differential expression analysis on the cleaned RNA-seq data."""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = (PROJECT_ROOT/"data")
RESULTS_DIR = (PROJECT_ROOT/"results"/"tables")

def main() -> None:

    expr = pd.read_csv(DATA_DIR/"GSE154918_counts_clean.csv.gz", index_col=0)
    metadata = pd.read_csv(DATA_DIR/"GSE154918_metadata_clean.csv.gz")

    expr_healthy = expr[metadata[metadata["condition"] == "healthy"]["sample_title"]]
    expr_sepsis = expr[metadata[metadata["condition"] == "sepsis"]["sample_title"]]

    log2FoldChanges = expr_sepsis.mean(axis=1) - expr_healthy.mean(axis=1)

    t_test = stats.ttest_ind(expr_sepsis, expr_healthy, axis=1, equal_var=False)
    t_stat, p_values = t_test.statistic, t_test.pvalue

    _, p_values_adj, _, _ = multipletests(p_values, method="fdr_bh")

    results = pd.DataFrame({
        "mean_sepsis": expr_sepsis.mean(axis=1),
        "mean_healthy": expr_healthy.mean(axis=1),
        "log2FoldChange": log2FoldChanges,
        "t_statistic": t_stat,
        "p_value": p_values,
        "p_value_adj": p_values_adj
    }).sort_values(by="p_value_adj")

    results["regulation"] = np.select(
        condlist=[
            (results["log2FoldChange"] > 1.0) & (results["p_value_adj"] < 0.05),
            (results["log2FoldChange"] < -1.0) & (results["p_value_adj"] < 0.05)
        ],
        choicelist=["upregulated", "downregulated"],
        default="not_significant"
    )

    results.to_csv(RESULTS_DIR/"GSE154918_differential_expression_results.csv", index=True)
    print(f"Exported {RESULTS_DIR/'GSE154918_differential_expression_results.csv'} {results.shape}")

if __name__ == "__main__":
    main()