"""Jaccard similarity matrix and hierarchical clustering."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


@dataclass
class ClusterResult:
    linkage: np.ndarray
    flat_assignment: np.ndarray
    distance_cut: float


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def jaccard_similarity_matrix(sets: list[set]) -> np.ndarray:
    n = len(sets)
    m = np.eye(n, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            v = jaccard(sets[i], sets[j])
            m[i, j] = v
            m[j, i] = v
    return m


def untyped(pairs: Iterable[tuple]) -> set:
    """Collapse bond_type from interaction pair tuples (residue_1, residue_2, bond)."""
    return {(p[0], p[1]) for p in pairs}


DEFAULT_SWEEP_CUTS: tuple[float, ...] = (0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70)


def sweep_cluster_cuts(
    linkage_matrix: np.ndarray,
    cuts: Sequence[float] = DEFAULT_SWEEP_CUTS,
) -> pd.DataFrame:
    """Return a table of cluster counts and size distributions at each cut.

    Helps the user pick `cluster_distance_cut` from data: a stable plateau
    (same cluster count over a range of cuts) is a defensible default; a
    fragmentation gradient calls for inspecting the dendrogram and density
    of singletons before deciding.

    Columns: cut, n_clusters, largest_5_sizes (list[int]), n_singletons.
    """
    rows = []
    for cut in cuts:
        flat = fcluster(linkage_matrix, t=cut, criterion="distance")
        sizes = sorted(Counter(flat.tolist()).values(), reverse=True)
        n_singletons = sum(1 for s in sizes if s == 1)
        rows.append({
            "cut": float(cut),
            "n_clusters": len(sizes),
            "largest_5_sizes": sizes[:5],
            "n_singletons": int(n_singletons),
        })
    return pd.DataFrame(rows)


def cluster_interfaces(
    similarity_matrix: np.ndarray,
    distance_cut: float,
    method: str = "average",
) -> ClusterResult:
    n = similarity_matrix.shape[0]
    if n < 2:
        return ClusterResult(
            linkage=np.empty((0, 4)),
            flat_assignment=np.ones(n, dtype=int),
            distance_cut=distance_cut,
        )
    distance = 1.0 - similarity_matrix
    np.fill_diagonal(distance, 0.0)
    distance = np.clip(distance, 0.0, None)
    condensed = squareform(distance, checks=False)
    Z = linkage(condensed, method=method)
    assignment = fcluster(Z, t=distance_cut, criterion="distance")
    return ClusterResult(linkage=Z, flat_assignment=assignment, distance_cut=distance_cut)
