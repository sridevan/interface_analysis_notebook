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
