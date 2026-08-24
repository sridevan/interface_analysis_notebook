"""Matplotlib figures for the notebook: similarity heatmap, dendrogram,
residue-pair frequency heatmap.

Each function draws its figure and returns None. Returning the Figure would
render it twice in a notebook, once from `plt.show()` and again from the cell's
display hook when the call is the last expression. Kept out of the notebook so
the narrative stays readable and the plotting logic can be tested.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, leaves_list

from pdbe_interfaces import outputs
from pdbe_interfaces.representation import InterfaceRecord
from pdbe_interfaces.similarity import ClusterResult


def similarity_heatmap(
    similarity_matrix: np.ndarray,
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    distance_cut: float,
) -> None:
    """Jaccard similarity heatmap, rows and columns in dendrogram leaf order.

    Leaf ordering places cluster members adjacent, so clusters appear as blocks
    on the diagonal rather than being scattered in retrieval order. White lines
    mark the boundaries between clusters.

    Does nothing when there are fewer than two interfaces to compare.
    """
    n = len(records)
    if n < 2:
        print("Only one interface; the similarity matrix is trivial.")
        return

    labels = [outputs.label_for_record(r) for r in records]
    leaf_order = leaves_list(cluster_result.linkage)
    reordered = similarity_matrix[leaf_order, :][:, leaf_order]
    labels_reordered = [labels[i] for i in leaf_order]

    flat = cluster_result.flat_assignment
    boundaries = [
        i for i in range(1, len(leaf_order))
        if flat[leaf_order[i]] != flat[leaf_order[i - 1]]
    ]

    _, ax = plt.subplots(figsize=(max(6, 0.5 * n + 3), max(5, 0.45 * n + 2)))
    sns.heatmap(
        reordered, ax=ax, cmap="viridis", square=True,
        cbar_kws={"label": "Jaccard similarity (typed)"},
        xticklabels=labels_reordered, yticklabels=labels_reordered,
        vmin=0.0, vmax=1.0,
    )
    for b in boundaries:
        ax.axhline(b, color="white", linewidth=0.5, alpha=0.6)
        ax.axvline(b, color="white", linewidth=0.5, alpha=0.6)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_title(
        f"Jaccard similarity in dendrogram leaf order "
        f"({len(set(flat))} clusters at cut {distance_cut})"
    )
    plt.tight_layout()
    plt.show()


def cluster_dendrogram(
    cluster_result: ClusterResult,
    records: list[InterfaceRecord],
    distance_cut: float,
) -> None:
    """Average-linkage dendrogram with the current cut drawn as a red line."""
    n = len(records)
    if n < 2:
        print("At least two interfaces are required for clustering.")
        return

    labels = [outputs.label_for_record(r) for r in records]
    _, ax = plt.subplots(figsize=(max(8, 0.5 * n + 3), 5))
    dendrogram(cluster_result.linkage, labels=labels, leaf_rotation=90, ax=ax)
    ax.axhline(
        distance_cut, color="red", linestyle="--", linewidth=1,
        label=f"distance cut = {distance_cut}",
    )
    ax.set_ylabel("1 - Jaccard")
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()
    print(
        f"Clusters at distance cut {distance_cut}: "
        f"{len(set(cluster_result.flat_assignment))}"
    )


def pair_frequency_heatmap(freq: dict, top_n: int, scope: str) -> None:
    """Residue-pair frequency heatmap, partner 1 residues by partner 2 residues.

    `freq` is the dict returned by `outputs.interface_frequency_summary`. Cell
    values are the fraction of interfaces in scope containing that pair. Each
    axis shows the top `top_n` residues by frequency, ordered by ascending
    UniProt residue number so positions read in sequence order.
    """
    matrix = freq["pair_matrix"]
    p1_lab, p2_lab = freq["partner_1_labels"], freq["partner_2_labels"]
    n_int = freq["n_interfaces"]
    p1_show, p2_show = min(len(p1_lab), top_n), min(len(p2_lab), top_n)
    if p1_show == 0 or p2_show == 0:
        print("No residue pairs to plot.")
        return

    p1_table, p2_table = freq["partner_1_residues"], freq["partner_2_residues"]
    p1_idx = sorted(range(p1_show), key=lambda i: int(p1_table.iloc[i]["position"]))
    p2_idx = sorted(range(p2_show), key=lambda i: int(p2_table.iloc[i]["position"]))
    matrix_show = matrix[p1_idx, :][:, p2_idx] / max(n_int, 1)

    _, ax = plt.subplots(
        figsize=(max(5, 0.32 * p2_show), max(4, 0.3 * p1_show))
    )
    sns.heatmap(
        matrix_show, ax=ax, cmap="viridis",
        xticklabels=[p2_lab[i] for i in p2_idx],
        yticklabels=[p1_lab[i] for i in p1_idx],
        vmin=0.0, vmax=1.0,
        cbar_kws={"label": f"fraction of {n_int}"},
        linewidths=0.3, linecolor="white",
    )
    ax.set_xlabel("Partner 2 residue")
    ax.set_ylabel("Partner 1 residue")
    ax.set_title(f"Top {p1_show} by {p2_show} pairs, {scope}")
    plt.tight_layout()
    plt.show()
