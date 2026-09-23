"""Frequency and conservation arithmetic: `pdbe_interfaces.outputs`.

Denominators are numbers of interfaces, a residue counts once per interface
however many contacts it makes, and thresholds are inclusive.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest

from pdbe_interfaces import outputs
from pdbe_interfaces.representation import InterfaceRecord
from pdbe_interfaces.similarity import ClusterResult

P1 = "P00001"
P2 = "P00002"
HB = "hydrogen_bond"
SB = "salt_bridge"


def _record(pdb_id, pairs):
    """Interface from (partner-1 position, partner-2 position, bond) tuples."""
    rec = InterfaceRecord(
        pdb_id=pdb_id, assembly_id="1", interface_id=1, interface_info={},
        unp_accession_1=P1, unp_accession_2=P2,
    )
    for (pos1, pos2, bond) in pairs:
        rec.uniprot_pairs.add(((P1, pos1, 1), (P2, pos2, 2), bond))
        rec.residue_identity[(P1, pos1, 1)] = "A"
        rec.residue_identity[(P2, pos2, 2)] = "G"
    return rec


@pytest.fixture
def records():
    """Partner-1 residue presence: 10 in 3/3, 20 in 2/3, 30 in 1/3, 40 in 1/3.

    Interface 1 has residue 10 in two contacts, to test once-per-interface counting.
    """
    return [
        _record("1aaa", [(10, 5, HB), (10, 6, HB), (20, 5, HB), (30, 5, HB)]),
        _record("2bbb", [(10, 5, HB), (20, 5, HB)]),
        _record("3ccc", [(10, 5, HB), (40, 5, HB)]),
    ]


def test_residue_and_pair_frequencies_use_interfaces_as_denominator(records):
    freq = outputs.interface_frequency_summary(records)
    assert freq["n_interfaces"] == 3

    p1 = freq["partner_1_residues"].set_index("position")
    assert p1["n_interfaces"].to_dict() == {10: 3, 20: 2, 30: 1, 40: 1}   # 10 is not 4
    assert p1.loc[10, "fraction"] == pytest.approx(1.0)
    assert p1.loc[20, "fraction"] == pytest.approx(2 / 3, abs=1e-3)
    p2 = freq["partner_2_residues"].set_index("position")
    assert p2["n_interfaces"].to_dict() == {5: 3, 6: 1}

    pairs = freq["pairs"].set_index("contact")["n_interfaces"]
    assert pairs["P00001:A10 – P00002:G5"] == 3
    assert pairs["P00001:A20 – P00002:G5"] == 2
    assert pairs["P00001:A10 – P00002:G6"] == 1

    matrix, rows, cols = freq["pair_matrix"], freq["partner_1_labels"], freq["partner_2_labels"]
    assert matrix[rows.index("P00001:A10"), cols.index("P00002:G5")] == 3
    assert matrix[rows.index("P00001:A40"), cols.index("P00002:G6")] == 0
    assert matrix.sum() == sum(pairs)


@pytest.mark.parametrize(
    "threshold, expected_partner_1",
    [
        (1.0, {10}),              # all interfaces
        (2 / 3, {10, 20}),        # exact equality is inclusive
        (0.67, {10}),             # just above 2/3 excludes residue 20
        (0.0, {10, 20, 30, 40}),  # everything seen at least once
    ],
    ids=["all", "inclusive-equality", "just-above", "zero"],
)
def test_conserved_residues_threshold_is_inclusive(records, threshold, expected_partner_1):
    conserved = outputs.conserved_residues(records, threshold=threshold)
    assert {pos for (acc, pos, role) in conserved if role == 1} == expected_partner_1
    assert (P2, 5, 2) in conserved


def test_conserved_pairs_are_typed(records):
    assert outputs.conserved_interaction_pairs(records, threshold=1.0) == {
        ((P1, 10, 1), (P2, 5, 2), HB),
    }
    # Same residues, different bond type in one interface: no longer in all three.
    records[2].uniprot_pairs.remove(((P1, 10, 1), (P2, 5, 2), HB))
    records[2].uniprot_pairs.add(((P1, 10, 1), (P2, 5, 2), SB))
    assert outputs.conserved_interaction_pairs(records, threshold=1.0) == set()
    assert outputs.conserved_interaction_pairs(records, threshold=2 / 3) == {
        ((P1, 10, 1), (P2, 5, 2), HB),
        ((P1, 20, 1), (P2, 5, 2), HB),
    }


def test_export_json_frequencies_levels_and_typed_vs_untyped(records, tmp_path):
    # Give interface 3 a second bond type on the (10, 5) pair.
    records[2].uniprot_pairs.add(((P1, 10, 1), (P2, 5, 2), SB))
    config = SimpleNamespace(output_dir=str(tmp_path))

    typed = outputs.export_interface_frequency_json(records, complex_id="T", config=config)
    residues = {(r["role"], r["unp_residue_number"]): r for r in typed["residue_frequencies"]}
    assert residues[(1, 10)]["n_interfaces"] == 3
    assert residues[(1, 10)]["conservation_level"] == "strong"      # 1.0
    assert residues[(1, 20)]["conservation_level"] == "medium"      # 0.667
    assert residues[(1, 30)]["conservation_level"] == "weak"        # 0.333
    assert residues[(1, 10)]["unp_residue_label"] == "A10"
    assert all(c["partner_1"]["role"] == 1 and c["partner_2"]["role"] == 2
               for c in typed["contact_frequencies"])
    typed_10_5 = {(c["bond_type"], c["n_interfaces"]) for c in typed["contact_frequencies"]
                  if (c["partner_1"]["unp_residue_number"], c["partner_2"]["unp_residue_number"]) == (10, 5)}
    assert typed_10_5 == {(HB, 3), (SB, 1)}

    untyped = outputs.export_interface_frequency_json(
        records, complex_id="U", config=config, typed_contacts=False,
    )
    untyped_10_5 = [c for c in untyped["contact_frequencies"]
                    if (c["partner_1"]["unp_residue_number"], c["partner_2"]["unp_residue_number"]) == (10, 5)]
    assert len(untyped_10_5) == 1
    assert untyped_10_5[0]["n_interfaces"] == 3     # once per interface, not 4

    written = json.loads((tmp_path / "T.json").read_text())
    assert written["metadata"]["n_interfaces"] == 3
    assert written["residue_frequencies"] == typed["residue_frequencies"]


# --- cluster-level frequencies -------------------------------------------------


@pytest.fixture
def clustered():
    """Two clusters of two. c0 everywhere, c1 only in cluster A, c2 only in B."""
    records = [
        _record("1aaa", [(1, 1, HB), (2, 2, HB)]),
        _record("2bbb", [(1, 1, HB), (2, 2, HB)]),
        _record("3ccc", [(1, 1, HB), (3, 3, SB)]),
        _record("4ddd", [(1, 1, HB), (3, 3, SB)]),
    ]
    result = ClusterResult(
        linkage=np.empty((0, 4)), flat_assignment=np.array([1, 1, 2, 2]), distance_cut=0.6,
    )
    return records, result


def test_rewiring_labels_from_cluster_fractions(clustered):
    records, result = clustered
    table = outputs.compare_cluster_contacts(records, result, cluster_a=1, cluster_b=2)
    by_contact = table.set_index("contact")

    shared = by_contact.loc["P00001:A1-P00002:G1 hydrogen_bond"]
    assert (shared.cluster_A_fraction, shared.cluster_B_fraction) == (1.0, 1.0)
    assert shared.contact_direction == "shared core"

    a_only = by_contact.loc["P00001:A2-P00002:G2 hydrogen_bond"]
    assert (a_only.cluster_A_count, a_only.cluster_B_count) == (2, 0)
    assert a_only.contact_direction == "higher in A"

    b_only = by_contact.loc["P00001:A3-P00002:G3 salt_bridge"]
    assert (b_only.cluster_A_count, b_only.cluster_B_count) == (0, 2)
    assert b_only.contact_direction == "higher in B"


def test_cluster_report_in_and_out_counts(clustered):
    records, result = clustered
    report = outputs.cluster_interpretation_report(records, result, outputs.AnnotationOverlap())

    assert report.cluster_size.tolist() == [2, 2]
    row_a = report[report.cluster_id == 1].iloc[0]
    assert row_a.member_pdb_ids == "1aaa,2bbb"
    # Contact in every member and absent from the rest of the dataset.
    assert "P00001:A2-P00002:G2 hydrogen_bond (2/2 interfaces = 100%, 2/2 entries; rest 0/2 = 0%)" \
        in row_a.cluster_contacts
    # Contact present everywhere.
    assert "P00001:A1-P00002:G1 hydrogen_bond (2/2 interfaces = 100%, 2/2 entries; rest 2/2 = 100%)" \
        in row_a.cluster_contacts
    assert row_a.cluster_contacts_total == 2
    assert not row_a.has_annotations
