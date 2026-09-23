"""Scientific regression: heterodimer partner correspondence, Spike RBD with ACE2.

Runs the notebook's workflow offline on the frozen fixture in
`tests/fixtures/recorded/spike_ace2_reversal/` (captured 2026-09-23), a
deliberately small subset of `PDB-CPX-140195`: five interface instances from
three entries, chosen so that the same biological interface is reported by PISA
in both partner orderings.

STING (tranche 2A) tested the full scientific state-analysis path on a
homodimer, where the partner-consistency check has nothing to act on. This
fixture tests the behaviour STING cannot exercise: a heterodimer whose partner
order differs between structural instances, which must be normalised so that
equivalent contacts become comparable and annotations stay attached to the
correct biological partner.

Biological identities were established from PDBe entity metadata and UniProt,
not from the workflow's own output:

    Q9BYF1  Angiotensin-converting enzyme 2 (human)      role 1
    P0DTC2  Spike glycoprotein (SARS-CoV-2)              role 2

    6m0j  chain A = ACE2, chain E = Spike
    7p19  chains A, B = ACE2; chains C, E = Spike
    7rpv  chains A-D = ACE2;  chains E-H = Spike

This is not a clustering study; the Spike/ACE2 conformational landscape is out
of scope here.
"""

from __future__ import annotations

import pytest

from pdbe_interfaces import representation, similarity

ACE2 = "Q9BYF1"
SPIKE = "P0DTC2"
HB = "hydrogen_bond"
SB = "salt_bridge"

CANONICAL = [("6m0j", "1", 1), ("7p19", "1", 1), ("7rpv", "1", 1)]
REVERSED = [("7p19", "2", 1), ("7rpv", "4", 1)]


def _contact(ace2_pos, spike_pos, bond):
    """Contact in the normalised representation: ACE2 as role 1, Spike as role 2."""
    return ((ACE2, ace2_pos, 1), (SPIKE, spike_pos, 2), bond)


def _by_label(records):
    return {r.label(): r for r in records}


@pytest.fixture(scope="module")
def raw_records(spike_ace2):
    """Records built from the same frozen input but WITHOUT partner normalisation."""
    return _by_label(representation.build_interface_records(spike_ace2.selection.interfaces))


@pytest.fixture(scope="module")
def norm_records(spike_ace2):
    """Records after the production partner-consistency check."""
    return _by_label(spike_ace2.all_records)


def test_spike_ace2_recorded_normalises_partner_orientation(spike_ace2, raw_records, norm_records):
    run = spike_ace2

    # Fixture integrity: every request answered from the frozen capture.
    assert len(run.recorded.urls) == run.recorded.metadata["n_requests_recorded"] == 22
    assert run.complex_id == "PDB-CPX-140195"
    assert run.details["oligomeric_state"] == "Heterodimer"

    # Canonical biological roles, asserted rather than assumed.
    assert run.partner_map == {(ACE2, 1): "ACE2", (SPIKE, 2): "S"}

    # The scoped cohort: five instances from three entries, nothing excluded.
    assert run.selection.n_interfaces_before == 5
    assert run.selection.pdb_ids == ["6m0j", "7p19", "7rpv"]
    assert [r.key for r in run.records] == CANONICAL[:1] + [("7p19", "1", 1), ("7p19", "2", 1),
                                                            ("7rpv", "1", 1), ("7rpv", "4", 1)]
    assert run.comparable.excluded == []

    # Before normalisation PISA reports two instances with the partners swapped.
    reversed_labels = {f"{p} a{a} i{i}" for p, a, i in REVERSED}
    for label, rec in raw_records.items():
        expected = (SPIKE, ACE2) if label in reversed_labels else (ACE2, SPIKE)
        assert (rec.unp_accession_1, rec.unp_accession_2) == expected, label

    # After normalisation every instance carries ACE2 in role 1 and Spike in role 2.
    for label, rec in norm_records.items():
        assert (rec.unp_accession_1, rec.unp_accession_2) == (ACE2, SPIKE), label


def test_spike_ace2_recorded_reverses_every_view_of_a_record_together(norm_records, raw_records):
    """Reversal must move contacts, author mappings and residue identities as one.

    `7rpv` assembly 4 is the cleanest case available: assemblies 1-3 of the same
    entry are reported canonically and assembly 4 with the partners swapped. Two
    assemblies of one deposition represent the same biological ACE2-RBD
    interface, which is a strong control for partner-order reversal, while still
    allowing genuine differences in the contacts detected in each copy.
    """
    a4, a1 = norm_records["7rpv a4 i1"], norm_records["7rpv a1 i1"]

    # Author chains: D is an ACE2 chain and H a Spike chain (PDBe entity metadata).
    assert {k[1] for k in a4.author_to_uniprot} == {"D", "H"}
    for (_, chain, _, _), (acc, _, role) in a4.author_to_uniprot.items():
        assert (acc, role) == ((ACE2, 1) if chain == "D" else (SPIKE, 2)), chain

    # Residue identities follow the new roles and agree with the canonical assembly.
    for key, aa in (((ACE2, 30, 1), "D"), ((ACE2, 41, 1), "Y"),
                    ((SPIKE, 417, 2), "K"), ((SPIKE, 500, 2), "T")):
        assert a4.residue_identity[key] == aa == a1.residue_identity[key]

    # Author-space pairs stay internally coherent: ACE2 chain first in every pair.
    assert all(first[1] == "D" and second[1] == "H" for first, second, _ in a4.author_pairs)

    # Every UniProt-keyed contact is oriented ACE2 role 1 -> Spike role 2.
    for k1, k2, _ in a4.uniprot_pairs:
        assert (k1[0], k1[2]) == (ACE2, 1) and (k2[0], k2[2]) == (SPIKE, 2)

    # Nothing was lost or invented by the reversal.
    assert len(a4.uniprot_pairs) == len(raw_records["7rpv a4 i1"].uniprot_pairs) == 15

    # Tranche-1.5 rule: per-residue mappings survive and stay self-consistent.
    assert set(a4.author_to_uniprot.values()) == set(a4.residue_identity)
    assert a4.n_residues_dropped_no_uniprot == 0
    assert a4.n_microheterogeneity_collisions == 0
    assert a4.author_to_uniprot[("7rpv", "D", 325, None)] == (ACE2, 325, 1)


def test_spike_ace2_recorded_preserves_equivalent_contacts_after_reversal(raw_records, norm_records):
    """Real ACE2-RBD contacts appear in opposite order before, and the same order after."""
    # Established ACE2-RBD interface contacts, verified present in this fixture.
    equivalents = [
        _contact(30, 417, SB),    # ACE2 Asp30 with Spike Lys417
        _contact(38, 498, HB),    # ACE2 Asp38 with Spike Gln498
        _contact(41, 500, HB),    # ACE2 Tyr41 with Spike Thr500
    ]
    a1_raw, a4_raw = raw_records["7rpv a1 i1"], raw_records["7rpv a4 i1"]
    a1, a4 = norm_records["7rpv a1 i1"], norm_records["7rpv a4 i1"]

    for contact in equivalents:
        (acc1, pos1, _), (acc2, pos2, _), bond = contact
        flipped = ((acc2, pos2, 1), (acc1, pos1, 2), bond)
        # Before: the canonical instance has it one way round, the reversed the other.
        assert contact in a1_raw.uniprot_pairs
        assert contact not in a4_raw.uniprot_pairs
        assert flipped in a4_raw.uniprot_pairs
        # After: both instances express it identically, ACE2 first.
        assert contact in a1.uniprot_pairs and contact in a4.uniprot_pairs

    # The two assemblies genuinely share a substantial contact core once comparable.
    assert len(a1.uniprot_pairs & a4.uniprot_pairs) == 11


def test_spike_ace2_recorded_similarity_is_not_broken_by_partner_order(raw_records, norm_records):
    """Partner order alone must not make equivalent interfaces look unrelated."""
    for canonical, reversed_ in (("7rpv a1 i1", "7rpv a4 i1"), ("7p19 a1 i1", "7p19 a2 i1")):
        before = similarity.jaccard(raw_records[canonical].uniprot_pairs,
                                    raw_records[reversed_].uniprot_pairs)
        after = similarity.jaccard(norm_records[canonical].uniprot_pairs,
                                   norm_records[reversed_].uniprot_pairs)
        # Reported in opposite orders, the two share no tuple at all.
        assert before == 0.0, canonical
        assert after > before
        # Not 1.0: the assemblies really do differ in a few detected contacts.
        assert after == pytest.approx({"7rpv a1 i1": 0.55, "7p19 a1 i1": 0.6923}[canonical], abs=1e-4)

    # Whole-cohort consequence: no instance is left isolated by partner order.
    sets = [r.uniprot_pairs for r in norm_records.values()]
    matrix = similarity.jaccard_similarity_matrix(sets)
    n = len(sets)
    assert min(matrix[i, j] for i in range(n) for j in range(n) if i != j) > 0.0


def test_spike_ace2_recorded_preserves_annotation_partner_identity(spike_ace2, norm_records):
    """An annotation on a reversed instance must stay on the right biological partner."""
    muts = spike_ace2.overlap.mutations
    assert spike_ace2.overlap.modifications.empty
    assert spike_ace2.overlap.ligands.empty      # no ligand reaches these interfaces

    rows = {(r.pdb_id, r.assembly_id): r for r in muts.itertuples()}
    assert set(rows) == {("7p19", "1"), ("7rpv", "1"), ("7rpv", "4")}

    # 7rpv assembly 4 is REVERSED. Q325Y sits on author chain D, an ACE2 chain,
    # and must resolve to ACE2 Gln325 in role 1 -- not to Spike residue 325.
    rev = rows[("7rpv", "4")]
    assert (rev.auth_asym_id, rev.auth_seq_id) == ("D", 325)
    assert (rev.unp_acc, rev.unp_seq_id, rev.role) == (ACE2, 325, 1)
    assert rev.mutation_label == "Q325Y"
    assert norm_records["7rpv a4 i1"].residue_identity[(ACE2, 325, 1)] == "Q"

    # The canonical assembly of the same entry carries the same mutation on a
    # different author chain and must resolve to the same biological identity.
    can = rows[("7rpv", "1")]
    assert (can.auth_asym_id, can.auth_seq_id) == ("A", 325)
    assert (can.unp_acc, can.unp_seq_id, can.role) == (ACE2, 325, 1)

    # The other partner: Q498Y on 7p19 chain E, a Spike chain, resolves to role 2.
    spike_mut = rows[("7p19", "1")]
    assert (spike_mut.auth_asym_id, spike_mut.auth_seq_id) == ("E", 498)
    assert (spike_mut.unp_acc, spike_mut.unp_seq_id, spike_mut.role) == (SPIKE, 498, 2)
    assert spike_mut.mutation_label == "Q498Y"
