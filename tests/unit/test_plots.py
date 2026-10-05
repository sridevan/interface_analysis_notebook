"""Plot options in `pdbe_interfaces.plots` that must agree with the clustering."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster

from pdbe_interfaces import plots
from pdbe_interfaces.representation import InterfaceRecord
from pdbe_interfaces.similarity import ClusterResult

# 0+1 join at 0.20, 2+3 at 0.30, everything at 1.00. scipy's default colour
# threshold is 70% of the tree height (0.70), far from the cuts used below.
TREE = np.array([
    [0, 1, 0.20, 2],
    [2, 3, 0.30, 2],
    [4, 5, 1.00, 4],
], dtype=float)


def _colours_agree_with_clusters(cut, **options):
    """Adjacent leaves share a non-neutral colour exactly when they share a cluster."""
    drawn = dendrogram(TREE, no_plot=True, **options)
    cluster = fcluster(TREE, t=cut, criterion="distance")
    leaves, colours = drawn["leaves"], drawn["leaves_color_list"]
    for i in range(len(leaves) - 1):
        same_cluster = cluster[leaves[i]] == cluster[leaves[i + 1]]
        same_colour = colours[i] == colours[i + 1] and colours[i] != plots.ABOVE_CUT_COLOUR
        if same_cluster != same_colour:
            return False
    return True


def test_dendrogram_colours_follow_the_clustering_cut():
    # Cut 0.25: clusters {0, 1}, {2}, {3}. Cut 0.30: the 2+3 link is exactly at
    # the cut, which fcluster includes in a cluster. Cut 0.50: {0, 1}, {2, 3}.
    for cut in (0.10, 0.25, 0.30, 0.50, 1.00):
        assert _colours_agree_with_clusters(cut, **plots.dendrogram_colour_options(cut)), cut


def test_default_dendrogram_colours_would_disagree_with_the_cut():
    # The bug this guards against: with scipy's default threshold, leaves 2 and
    # 3 share a colour although a cut of 0.25 puts them in different clusters.
    assert not _colours_agree_with_clusters(0.25, above_threshold_color=plots.ABOVE_CUT_COLOUR)


def test_cluster_dendrogram_passes_the_cut_to_the_plot(monkeypatch):
    seen = {}

    def fake_dendrogram(linkage, **kwargs):
        seen.update(kwargs)
        return dendrogram(linkage, **kwargs)

    monkeypatch.setattr(plots, "dendrogram", fake_dendrogram)
    monkeypatch.setattr(plots.plt, "show", lambda: None)
    records = [
        InterfaceRecord(pdb_id=f"{i}abc", assembly_id="1", interface_id=1, interface_info={},
                        unp_accession_1="P1", unp_accession_2="P2")
        for i in range(4)
    ]
    result = ClusterResult(linkage=TREE, flat_assignment=fcluster(TREE, 0.25, "distance"),
                           distance_cut=0.25)
    plots.cluster_dendrogram(result, records, 0.25)
    assert seen["color_threshold"] > 0.25 and seen["color_threshold"] < 0.25 + 1e-9
    assert seen["above_threshold_color"] == plots.ABOVE_CUT_COLOUR
