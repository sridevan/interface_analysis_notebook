"""Jaccard similarity and hierarchical clustering: `pdbe_interfaces.similarity`."""

from __future__ import annotations

import numpy as np
import pytest

from pdbe_interfaces import similarity


@pytest.mark.parametrize(
    "a, b, expected",
    [
        ({1, 2, 3}, {1, 2, 3}, 1.0),
        ({1, 2}, {3, 4}, 0.0),
        ({1, 2, 3}, {2, 3, 4}, 0.5),      # intersection 2, union 4
        # Documented behaviour: two interfaces with no mapped contacts score as
        # identical, so they would cluster together.
        (set(), set(), 1.0),
        (set(), {1}, 0.0),
    ],
    ids=["identical", "disjoint", "partial", "empty-vs-empty", "empty-vs-nonempty"],
)
def test_jaccard(a, b, expected):
    assert similarity.jaccard(a, b) == pytest.approx(expected)
    assert similarity.jaccard(b, a) == pytest.approx(expected)


def test_similarity_matrix_is_symmetric_with_unit_diagonal():
    sets = [{1, 2, 3}, {2, 3, 4}, {9}]
    m = similarity.jaccard_similarity_matrix(sets)

    assert m.shape == (3, 3)
    assert np.array_equal(m, m.T)
    assert np.array_equal(np.diag(m), np.ones(3))
    assert m[0, 1] == pytest.approx(0.5)
    assert m[0, 2] == 0.0
    assert m[1, 2] == 0.0


def test_clustering_groups_identical_and_separates_disjoint():
    sets = [{"c1", "c2", "c3"}, {"c1", "c2", "c3"}, {"x1", "x2"}]
    sim = similarity.jaccard_similarity_matrix(sets)
    labels = similarity.cluster_interfaces(sim, distance_cut=0.5).flat_assignment

    assert len(labels) == 3
    assert labels[0] == labels[1]
    assert labels[0] != labels[2]


def test_cut_is_in_one_minus_jaccard_units():
    """`cluster_distance_cut` merges pairs whose distance is at or below it."""
    sets = [{1, 2, 3, 4}, {3, 4, 5, 6}, {9}]      # jaccard(A, B) = 2/6
    sim = similarity.jaccard_similarity_matrix(sets)
    d_ab = 1 - 2 / 6
    tight = similarity.cluster_interfaces(sim, distance_cut=d_ab - 0.05).flat_assignment
    loose = similarity.cluster_interfaces(sim, distance_cut=d_ab + 0.05).flat_assignment

    assert tight[0] != tight[1]
    assert loose[0] == loose[1]
    assert loose[0] != loose[2]


def test_single_interface_is_one_cluster():
    sim = similarity.jaccard_similarity_matrix([{1, 2}])
    result = similarity.cluster_interfaces(sim, distance_cut=0.6)

    assert result.flat_assignment.tolist() == [1]
    assert result.linkage.shape == (0, 4)


def test_typed_and_untyped_diverge_only_on_bond_type():
    typed = [{("r1", "r2", "hydrogen_bond")}, {("r1", "r2", "salt_bridge")}]
    untyped = [similarity.untyped(p) for p in typed]

    assert untyped == [{("r1", "r2")}, {("r1", "r2")}]
    assert similarity.jaccard_similarity_matrix(typed)[0, 1] == 0.0
    assert similarity.jaccard_similarity_matrix(untyped)[0, 1] == 1.0


# --- how close the cut is to a linkage ---------------------------------------

# Four interfaces: 0+1 join at 0.20, 2+3 at 0.60, everything at 0.90.
TREE = np.array([
    [0, 1, 0.20, 2],
    [2, 3, 0.60, 2],
    [4, 5, 0.90, 4],
], dtype=float)


def test_cut_sensitivity_between_two_linkages():
    assert similarity.describe_cut_sensitivity(TREE, 0.40) == (
        "At cut 0.40, the analysis gives 3 interaction groups. The nearest linkage below "
        "the cut is at 0.20, and the next linkage above the cut is at 0.60."
    )


def test_cut_sensitivity_states_a_linkage_exactly_at_the_cut():
    text = similarity.describe_cut_sensitivity(TREE, 0.60)
    # fcluster joins at distance <= cut, so the 0.60 linkage is inside a group.
    assert text.startswith("At cut 0.60, the analysis gives 2 interaction groups.")
    assert "A linkage occurs exactly at the cut (0.60)" in text
    assert "the nearest linkage below the cut is at 0.20" in text
    assert "the next linkage above the cut is at 0.90" in text


def test_cut_sensitivity_boundaries():
    low = similarity.describe_cut_sensitivity(TREE, 0.10)
    assert low.startswith("At cut 0.10, the analysis gives 4 interaction groups.")
    assert "No linkage lies below the cut" in low and "next linkage above the cut is at 0.20" in low

    high = similarity.describe_cut_sensitivity(TREE, 0.95)
    assert high.startswith("At cut 0.95, the analysis gives 1 interaction group.")
    assert "nearest linkage below the cut is at 0.90" in high and "no linkage lies above the cut" in high

    assert "no clustering cut" in similarity.describe_cut_sensitivity(np.empty((0, 4)), 0.6)


def test_cut_sensitivity_does_not_print_a_near_linkage_as_the_cut():
    near = TREE.copy()
    near[1, 2] = 0.604
    assert "next linkage above the cut is at 0.604" in similarity.describe_cut_sensitivity(near, 0.60)


def test_cut_sensitivity_group_count_matches_cluster_interfaces():
    sets = [{1, 2, 3}, {1, 2, 3}, {1, 2, 4}, {7, 8}, {7, 9}]
    result = similarity.cluster_interfaces(similarity.jaccard_similarity_matrix(sets), 0.6)
    n_groups = len(set(result.flat_assignment.tolist()))
    assert f"gives {n_groups} interaction groups" in similarity.describe_cut_sensitivity(
        result.linkage, 0.6)
