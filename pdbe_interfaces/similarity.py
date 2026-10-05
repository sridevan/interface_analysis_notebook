# Copyright 2026 EMBL - European Bioinformatics Institute
# Author: Sri Devan Appasamy
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
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


def describe_cut_sensitivity(linkage_matrix: np.ndarray, distance_cut: float) -> str:
    """Say how close `distance_cut` is to the linkage distances around it.

    Distances are on the dendrogram's scale, `1 - Jaccard`. A linkage close to
    the cut means a small change in the cut changes the number of clusters.
    Clusters are not named: their ids are labels, not stable identities.
    """
    heights = np.asarray(linkage_matrix, dtype=float).reshape(-1, 4)[:, 2]
    if heights.size == 0:
        return "Fewer than two interface instances: there is no clustering cut to assess."
    n_groups = len(set(fcluster(linkage_matrix, t=distance_cut, criterion="distance").tolist()))
    at_cut = np.isclose(heights, distance_cut, rtol=0.0, atol=1e-9)
    below = heights[(heights < distance_cut) & ~at_cut]
    above = heights[(heights > distance_cut) & ~at_cut]

    def fmt(value: float) -> str:
        # Two decimals, unless that would print a different distance as the cut itself.
        text = f"{value:.2f}"
        return text if text != f"{distance_cut:.2f}" or value == distance_cut else f"{value:.3f}"

    plural = "group" if n_groups == 1 else "groups"
    first = f"At cut {distance_cut:.2f}, the analysis gives {n_groups} interaction {plural}."
    parts = []
    if below.size:
        parts.append(f"the nearest linkage below the cut is at {fmt(below.max())}")
    else:
        parts.append("no linkage lies below the cut")
    if above.size:
        parts.append(f"the next linkage above the cut is at {fmt(above.min())}")
    else:
        parts.append("no linkage lies above the cut")
    second = (parts[0][0].upper() + parts[0][1:]) + ", and " + parts[1] + "."
    if at_cut.any():
        second = (
            f"A linkage occurs exactly at the cut ({distance_cut:.2f}), so the number of "
            f"groups changes with any small change in the cut; " + parts[0] + ", and "
            + parts[1] + "."
        )
    return f"{first} {second}"


def describe_typed_untyped_divergence(
    sim_typed: np.ndarray, sim_untyped: np.ndarray,
) -> str:
    """Compare typed and untyped Jaccard matrices.

    Substantial divergence, above roughly 0.2, indicates that bond-type
    assignment contributes more to the clustering than the residue-pair
    topology does. Bond types are geometry-derived and resolution-dependent,
    so in that situation the clustering is more sensitive to data quality.
    """
    diff = float(np.abs(sim_typed - sim_untyped).max())
    return f"Maximum absolute difference between typed and untyped Jaccard: {diff:.3f}"
