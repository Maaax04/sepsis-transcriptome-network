"""Constructs a protein-protein interaction network from the differential expression results."""

from pathlib import Path
from io import StringIO
from matplotlib import pyplot as plt
from adjustText import adjust_text

import pandas as pd
import networkx as nx
import requests


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = (PROJECT_ROOT/"data")
RESULTS_DIR = (PROJECT_ROOT/"results")
TABLE_DIR = (RESULTS_DIR/"tables").resolve()
FIGURE_DIR = (RESULTS_DIR/"figures").resolve()
URL = "https://string-db.org/api/tsv/network"

def main() -> None:

    de_results = pd.read_csv(TABLE_DIR/"GSE154918_differential_expression_results.csv", index_col="gene_symbol")

    significant_genes = de_results[de_results["regulation"] != "not_significant"].head(400).index.to_list()

    genes = "\r".join(significant_genes)

    payload = {
        "identifiers": genes,
        "species": 9606,
        "required_score": 700,
        "network_type": "functional",
        "caller_identity": "sepsis-transcriptome-network"
    }

    response = requests.post(URL, data=payload, timeout=30)
    response.raise_for_status()

    interactions = pd.read_csv(StringIO(response.text), sep="\t")

    print(f"Fetched {len(interactions)} interactions for {len(significant_genes)} significant genes.")


    G = nx.Graph()
    for _, row in interactions.iterrows():
        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

    print(f"Constructed PPI network with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

    connected_components = nx.connected_components(G)

    largest_connected_component = max(connected_components, key=len)
    G_largest = G.subgraph(largest_connected_component).copy()
    print(f"Constructed largest connected component with {G_largest.number_of_nodes()} nodes and {G_largest.number_of_edges()} edges.")

    degrees = dict(G_largest.degree())
    betweenness_centrality = nx.betweenness_centrality(G_largest)
    eigenvector_centrality = nx.eigenvector_centrality(G_largest)

    hub_df = pd.DataFrame({
        "degree": pd.Series(degrees),
        "betweenness_centrality": pd.Series(betweenness_centrality),
        "eigenvector_centrality": pd.Series(eigenvector_centrality)
    })
    hub_df.index.name = "gene_symbol"

    hub_df = hub_df.join(de_results[["log2FoldChange", "p_value_adj", "regulation"]], how="left")

    hub_df = hub_df.sort_values(by=["degree", "betweenness_centrality", "eigenvector_centrality"], ascending=False)

    hub_df.to_csv(TABLE_DIR/"GSE154918_ppi_network_hub_genes.csv", index=True)

    plot_ppi_network(G_largest, hub_df, FIGURE_DIR/"GSE154918_ppi_network.png")

def plot_ppi_network(G: nx.Graph, hub_df: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 12))

    color_map = {
        "upregulated": "crimson",
        "downregulated": "royalblue",
    }

    node_colors = []
    node_sizes = []

    for node in G.nodes():
        if node in hub_df.index:
            regulation = hub_df.loc[node, "regulation"]
            node_colors.append(color_map.get(regulation, "lightgray"))
            size = hub_df.loc[node, "degree"] * 60 + 80
            node_sizes.append(size)
        else:
            node_colors.append("lightgray")
            node_sizes.append(80)

    pos = nx.spring_layout(G, seed=42, k=0.75)
    nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color="slategray")
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.85, edgecolors="white", linewidths=0.8)

    top_genes = hub_df.head(20).index.to_list()

    x_coords = [pos[node][0] for node in G.nodes()]
    y_coords = [pos[node][1] for node in G.nodes()]

    texts = []
    for gene in top_genes:
        if gene in pos:
            x, y = pos[gene]
            text = plt.text(x, y - 0.035, gene, fontsize=10, fontweight="bold", color="black", ha="center", va="top")
            texts.append(text)
    adjust_text(
        texts,
        x=x_coords,
        y=y_coords,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
    )

    plt.title("Protein-Protein Interaction Network of Significant Genes", fontsize=16)
    plt.axis("off")
    ax.scatter([], [], c="crimson", label="Upregulated", s=100)
    ax.scatter([], [], c="royalblue", label="Downregulated", s=100)
    plt.legend(scatterpoints=1, frameon=True, fontsize=12, loc="lower left", title="Gene Regulation")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()