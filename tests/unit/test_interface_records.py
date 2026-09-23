"""Residue and contact correspondence: `representation.build_interface_records`.

The invariant under test: equivalent structural residues in different PDB
entries must land on the same UniProt-keyed tuple, while author-space keys stay
available for the annotation join.
"""

from __future__ import annotations

import pytest

from pdbe_interfaces import representation, similarity

X = "P00001"
Y = "P00002"


def _build(*interfaces):
    return representation.build_interface_records(list(interfaces))


def test_one_contact_yields_author_and_uniprot_keys_with_roles(make_contact, make_interface):
    contact = make_contact("A", 10, "B", 20, unp1=110, unp2=220, aa1="K", aa2="D")
    [rec] = _build(make_interface("1abc", [contact]))

    assert (rec.pdb_id, rec.assembly_id, rec.interface_id) == ("1abc", "1", 1)
    assert (rec.unp_accession_1, rec.unp_accession_2) == (X, Y)
    assert rec.author_pairs == {
        (("1abc", "A", 10, None), ("1abc", "B", 20, None), "hydrogen_bond"),
    }
    assert rec.uniprot_pairs == {((X, 110, 1), (Y, 220, 2), "hydrogen_bond")}
    assert rec.author_to_uniprot == {
        ("1abc", "A", 10, None): (X, 110, 1),
        ("1abc", "B", 20, None): (Y, 220, 2),
    }
    assert rec.residue_identity == {(X, 110, 1): "K", (Y, 220, 2): "D"}
    assert rec.n_residues_dropped_no_uniprot == 0
    assert rec.n_microheterogeneity_collisions == 0


def test_comparison_space_is_keyed_by_uniprot_not_author_numbering(make_contact, make_interface):
    """Different chains and author numbers, same UniProt positions: identical.
    Same chains and author numbers, different UniProt positions: unrelated."""
    entry_a = make_interface("1abc", [
        make_contact("A", 10, "B", 20, unp1=110, unp2=220),
        make_contact("A", 11, "B", 21, unp1=111, unp2=221, bond="salt_bridge"),
    ])
    entry_b = make_interface("2xyz", [
        make_contact("H", 310, "L", 520, unp1=110, unp2=220),
        make_contact("H", 311, "L", 521, unp1=111, unp2=221, bond="salt_bridge"),
    ])
    entry_c = make_interface("3qqq", [
        make_contact("A", 10, "B", 20, unp1=310, unp2=220),
    ])
    rec_a, rec_b, rec_c = _build(entry_a, entry_b, entry_c)

    assert rec_a.uniprot_pairs == rec_b.uniprot_pairs
    assert rec_a.author_pairs.isdisjoint(rec_b.author_pairs)
    assert similarity.jaccard(rec_a.uniprot_pairs, rec_b.uniprot_pairs) == 1.0
    assert rec_a.uniprot_pairs.isdisjoint(rec_c.uniprot_pairs)


def test_role_follows_pisa_side_not_accession_order(make_contact, make_interface):
    contact = make_contact("A", 10, "B", 20, acc1="Q99999", acc2="A00001")
    [rec] = _build(make_interface("1abc", [contact]))

    assert (rec.unp_accession_1, rec.unp_accession_2) == ("Q99999", "A00001")
    assert rec.uniprot_pairs == {(("Q99999", 10, 1), ("A00001", 20, 2), "hydrogen_bond")}


@pytest.mark.parametrize(
    "raw, expected",
    [(None, None), ("", None), (" ", None), (" B ", "B")],
    ids=["none", "empty", "blank", "padded-code"],
)
def test_insertion_code_normalisation(raw, expected, make_contact, make_interface):
    [rec] = _build(make_interface("1abc", [make_contact("A", 10, "B", 20, ins1=raw)]))
    assert {a for (a, _, _) in rec.author_pairs} == {("1abc", "A", 10, expected)}


def test_duplicate_observations_collapse_but_bond_types_stay_distinct(make_contact, make_interface):
    """Atom-level repeats count once; a second bond type is a second typed tuple."""
    hbond = make_contact("A", 10, "B", 20, bond="hydrogen_bond")
    salt = make_contact("A", 10, "B", 20, bond="salt_bridge")
    [rec] = _build(make_interface("1abc", [hbond, dict(hbond), dict(hbond), salt]))

    assert len(rec.author_pairs) == 2
    assert len(rec.uniprot_pairs) == 2
    assert len(similarity.untyped(rec.uniprot_pairs)) == 1


def test_contact_missing_uniprot_on_one_side(make_contact, make_interface):
    """The pair leaves the comparison space; the mapped residue keeps its mapping.

    Residue correspondence and pair comparability are different things: a
    residue with an accession and a position is mapped whether or not its
    contact partner is, but the pair is only comparable if both sides map.
    """
    mapped = make_contact("A", 10, "B", 20)
    unmapped_pos = make_contact("A", 11, "B", 21, unp2=None)
    unmapped_acc = make_contact("A", 12, "B", 22, acc1=None)
    [rec] = _build(make_interface("1abc", [mapped, unmapped_pos, unmapped_acc]))

    assert len(rec.author_pairs) == 3
    assert rec.uniprot_pairs == {((X, 10, 1), (Y, 20, 2), "hydrogen_bond")}
    assert rec.n_residues_dropped_no_uniprot == 2
    assert rec.author_to_uniprot == {
        ("1abc", "A", 10, None): (X, 10, 1),
        ("1abc", "B", 20, None): (Y, 20, 2),
        ("1abc", "A", 11, None): (X, 11, 1),     # mapped side of a half-mapped contact
        ("1abc", "B", 22, None): (Y, 22, 2),     # mapped side of a half-mapped contact
    }
    assert rec.residue_identity.keys() == set(rec.author_to_uniprot.values())


def test_mapped_residue_is_consistent_across_partial_and_full_contacts(make_contact, make_interface):
    contacts = [
        make_contact("A", 10, "B", 21, unp2=None),   # A10 seen first in a half-mapped contact
        make_contact("A", 10, "B", 20),               # then in a fully mapped one
    ]
    [rec] = _build(make_interface("1abc", contacts))

    assert rec.author_to_uniprot[("1abc", "A", 10, None)] == (X, 10, 1)
    assert rec.uniprot_pairs == {((X, 10, 1), (Y, 20, 2), "hydrogen_bond")}
    assert rec.n_microheterogeneity_collisions == 0


@pytest.mark.parametrize(
    "contacts_spec, expected",
    [
        # No accession anywhere: nothing to infer.
        ([dict(acc1=None, acc2=None), dict(acc1=None, acc2=None)], ("", "")),
        # Accessions present but positions missing: accessions are still taken,
        # although the contact itself is dropped from the comparison space.
        ([dict(unp1=None, unp2=None)], (X, Y)),
        # First contact lacks one accession; the next fully mapped one is used.
        ([dict(acc2=None), dict()], (X, Y)),
    ],
    ids=["no-accessions", "accessions-without-positions", "first-mapped-contact"],
)
def test_partner_accessions_come_from_first_contact_carrying_both(
    contacts_spec, expected, make_contact, make_interface,
):
    contacts = [make_contact("A", 10 + i, "B", 20 + i, **spec) for i, spec in enumerate(contacts_spec)]
    [rec] = _build(make_interface("1abc", contacts))
    assert (rec.unp_accession_1, rec.unp_accession_2) == expected


def test_microheterogeneity_keeps_first_seen_author_key(make_contact, make_interface):
    """Two author residues (A10 and A10A) mapping to one UniProt position."""
    contacts = [
        make_contact("A", 10, "B", 20, unp1=110, unp2=220),
        make_contact("A", 10, "B", 20, unp1=110, unp2=220, ins1="A"),
    ]
    [rec] = _build(make_interface("1abc", contacts))

    assert rec.n_microheterogeneity_collisions == 1
    assert len(rec.author_pairs) == 2
    assert rec.uniprot_pairs == {((X, 110, 1), (Y, 220, 2), "hydrogen_bond")}
    assert rec.author_to_uniprot[("1abc", "A", 10, None)] == (X, 110, 1)
    assert ("1abc", "A", 10, "A") not in rec.author_to_uniprot


def test_one_record_per_interface_in_input_order(make_contact, make_interface):
    """Cluster labels are matched to records by position, so order must hold."""
    items = [
        make_interface("1abc", [make_contact("A", 10, "B", 20)], assembly_id=1, interface_id=1),
        make_interface("1abc", [make_contact("C", 10, "D", 20)], assembly_id=2, interface_id=1),
        make_interface("2xyz", [make_contact("A", 10, "B", 20)], assembly_id=1, interface_id=3),
    ]
    records = _build(*items)
    assert [r.key for r in records] == [("1abc", "1", 1), ("1abc", "2", 1), ("2xyz", "1", 3)]


# --- comparable-record selection -------------------------------------------


def test_select_comparable_records_excludes_zero_pair_interfaces(make_contact, make_interface):
    """Eligibility for contact-based comparison needs >= 1 UniProt pair.

    An interface with PISA contacts but no comparable pair carries insufficient
    comparable information; it is not "identical to another empty interface".
    """
    comparable = make_interface("1aaa", [make_contact("A", 10, "B", 20)])
    empty = make_interface("2bbb", [
        make_contact("A", 10, "B", 20, unp2=None),
        make_contact("A", 11, "B", 21, acc1=None),
    ])
    partial = make_interface("3ccc", [
        make_contact("A", 10, "B", 20, unp2=None),
        make_contact("A", 11, "B", 21),            # one usable pair is enough
    ])
    records = _build(comparable, empty, partial)

    selection = representation.select_comparable_records(records)

    assert [r.pdb_id for r in selection.records] == ["1aaa", "3ccc"]
    assert [r.pdb_id for r in selection.excluded] == ["2bbb"]
    assert selection.reasons == {("2bbb", "1", 1): representation.NO_COMPARABLE_CONTACTS}
    assert selection.excluded[0].author_pairs            # retained for provenance
    assert "2bbb" in selection.summary()


def test_excluded_record_does_not_reach_similarity_clustering_or_frequencies(
    make_contact, make_interface,
):
    """Three comparable interfaces share contact X; a fourth has no comparable pair.

    The frequency of X is 3/3, not 3/4: a missing mapping is not biological
    variation. The matrix and the cluster labels cover the three only.
    """
    from pdbe_interfaces import outputs, similarity

    items = [make_interface(pdb, [make_contact("A", 10, "B", 20)]) for pdb in ("1aaa", "2bbb", "3ccc")]
    items.append(make_interface("4ddd", [make_contact("A", 10, "B", 20, unp2=None)]))
    selection = representation.select_comparable_records(_build(*items))
    records = selection.records

    sim = similarity.jaccard_similarity_matrix([r.uniprot_pairs for r in records])
    labels = similarity.cluster_interfaces(sim, distance_cut=0.6).flat_assignment
    assert sim.shape == (3, 3)
    assert len(labels) == 3
    assert len(set(labels.tolist())) == 1

    freq = outputs.interface_frequency_summary(records)
    assert freq["n_interfaces"] == 3
    assert freq["pairs"]["fraction"].tolist() == [1.0]
    assert outputs.conserved_residues(records, threshold=1.0) == {(X, 10, 1), (Y, 20, 2)}
    assert outputs.conserved_interaction_pairs(records, threshold=1.0) == {
        ((X, 10, 1), (Y, 20, 2), "hydrogen_bond"),
    }
