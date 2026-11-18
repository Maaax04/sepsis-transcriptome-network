"""Visualizes the differential expression results."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from adjustText import adjust_text

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = (PROJECT_ROOT/"data")
TABLES_DIR = (PROJECT_ROOT/"results"/"tables")
FIGURES_DIR = (PROJECT_ROOT/"results"/"figures")

def main() -> None:
    results = pd.read_csv(TABLES_DIR/"GSE154918_differential_expression_results.csv", index_col="gene_symbol")

    results["neg_log10_padj"] = -np.log10(np.clip(results["p_value_adj"], 1e-300, None))

    fig, ax = plt.subplots(figsize=(10,7))

    groups = {
        "not_significant": {"color": "lightgray", "label": "Not Significant"},
        "upregulated": {"color": "crimson", "label": "Upregulated"},
        "downregulated": {"color": "royalblue", "label": "Downregulated"}
    }

    for group, props in groups.items():
        subset = results[results["regulation"] == group]
        ax.scatter(subset["log2FoldChange"], subset["neg_log10_padj"],
                   color=props["color"], label=props["label"], alpha=0.7, edgecolors="none", s=12)
    ax.legend()
    ax.axvline(x=1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(x=-1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(y=-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.6)

    top_genes =results[results["regulation"] != "not_significant"].head(10)
    texts = []
    for gene, row in top_genes.iterrows():
        text = ax.text(row["log2FoldChange"], row["neg_log10_padj"], gene, fontsize=8)
        texts.append(text)

    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.6)
    )

    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("-Log10 Adjusted P-value")
    ax.set_title("Volcano Plot of Differential Expression Results")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(FIGURES_DIR/"volcano_plot.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()