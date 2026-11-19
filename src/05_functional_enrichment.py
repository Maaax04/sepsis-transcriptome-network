"""Functional Enrichment Analysis"""
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from matplotlib import pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = (PROJECT_ROOT/"data")
RESULTS_DIR = (PROJECT_ROOT/"results")
TABLE_DIR = (RESULTS_DIR/"tables").resolve()
FIGURE_DIR = (RESULTS_DIR/"figures").resolve()
API_URL = "https://string-db.org/api/tsv/enrichment"

def fetch_enrichment(gene_list: list[str], direction: str) -> pd.DataFrame:
    genes = "\r".join(gene_list)

    payload = {
        "identifiers": genes,
        "species": 9606,
        "caller_identity": "sepsis-transcriptome-network"
    }

    response = requests.post(API_URL, data=payload, timeout=30)
    response.raise_for_status()

    enrichment_results = pd.read_csv(StringIO(response.text), sep="\t")
    enrichment_results = enrichment_results[enrichment_results["fdr"] < 0.05]
    enrichment_results["direction"] = direction

    return enrichment_results

def main() -> None:
    de_results = pd.read_csv(TABLE_DIR/"GSE154918_differential_expression_results.csv", index_col="gene_symbol")

    upregulated_genes = de_results[de_results["regulation"] == "upregulated"].head(400).index.to_list()
    downregulated_genes = de_results[de_results["regulation"] == "downregulated"].head(400).index.to_list()

    up_enrichment = fetch_enrichment(upregulated_genes, "upregulated")
    down_enrichment = fetch_enrichment(downregulated_genes, "downregulated")

    enrichment_results = pd.concat([up_enrichment, down_enrichment], ignore_index=True)

    enrichment_results.to_csv(TABLE_DIR/"GSE154918_functional_enrichment_results.csv", index=False)
    print(f"Exported {TABLE_DIR/'GSE154918_functional_enrichment_results.csv'} {enrichment_results.shape}")

    plot_enrichment(enrichment_results, FIGURE_DIR/"GSE154918_functional_enrichment_plot.png")

def plot_enrichment(enrichment_results: pd.DataFrame, output_path: Path) -> None:
    neg_log10_fdr = -np.log10(np.clip(enrichment_results["fdr"], 1e-300, None))
    enrichment_results = enrichment_results.assign(neg_log10_fdr=neg_log10_fdr)

    fig, axes = plt.subplots(2, 2, figsize=(20, 12))

    configs = [
        (axes[0, 0], "Process", "upregulated", "GO Biological Processes (Upregulated)", "crimson"),
        (axes[0, 1], "Process", "downregulated", "GO Biological Processes (Downregulated)", "royalblue"),
        (axes[1, 0], "RCTM", "upregulated", "Reactome Processes (Upregulated)", "crimson"),
        (axes[1, 1], "RCTM", "downregulated", "Reactome Processes (Downregulated)", "royalblue")
    ]

    for ax, category, direction, title, color in configs:
        subset = enrichment_results[(enrichment_results["category"] == category) & (enrichment_results["direction"] == direction)].sort_values("fdr", ascending=True).head(10)
        ax.barh(subset["description"], subset["neg_log10_fdr"], color=color)
        ax.invert_yaxis()
        ax.axvline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel("-log10(FDR)")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    main()