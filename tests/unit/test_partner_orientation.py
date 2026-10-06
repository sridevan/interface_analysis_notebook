"""Partner orientation: `representation.check_partner_consistency`.

Heterodimer: the majority `(unp_accession_1, unp_accession_2)` ordering is
canonical and minority records have their roles reversed. Homodimer: both
accessions are equal, so the check cannot act and has nothing to do; the two
copies are instead made interchangeable in the contact key itself, by
`representation.canonical_uniprot_pair`, which these tests cover.
"""

from __future__ import annotations

import pytest

from pdbe_interfaces import representation, similarity

X = "P0000X"
Y = "P0000Y"
H = "P0HOMO"


def _canonical(make_contact, make_interface, pdb_id, c1, c2):
    """X (role 1) against Y (role 2); chain and author numbers vary per entry."""
    return make_interface(pdb_id, [
        make_contact(c1, 10, c2, 20, acc1=X, acc2=Y, unp1=10, unp2=20, aa1="K", aa2="D"),
        make_contact(c1, 11, c2, 21, acc1=X, acc2=Y, unp1=11, unp2=21, aa1="R", aa2="E",
                     bond="salt_bridge"),
    ])


def _reversed(make_contact, make_interface, pdb_id="2bbb"):
    """The same interface with PISA listing Y first: Y (role 1) against X (role 2)."""
    return make_interface(pdb_id, [
        make_contact("L", 220, "H", 110, acc1=Y, acc2=X, unp1=20, unp2=10, aa1="D", aa2="K"),
        make_contact("L", 221, "H", 111, acc1=Y, acc2=X, unp1=21, unp2=11, aa1="E", aa2="R",
                     bond="salt_bridge"),
    ])


def _check(items):
    return representation.check_partner_consistency(
        representation.build_interface_records(items)
    )


def test_minority_ordering_is_reversed_onto_the_canonical_roles(make_contact, make_interface):
    records = representation.build_interface_records([
        _canonical(make_contact, make_interface, "1aaa", "A", "B"),
        _reversed(make_contact, make_interface, "2bbb"),
        _canonical(make_contact, make_interface, "3ccc", "C", "D"),
    ])
    canonical, reversed_, other = records
    canonical_pairs_before = set(canonical.uniprot_pairs)
    canonical_a2u_before = dict(canonical.author_to_uniprot)

    # Before the check the reversed record shares nothing with its equivalents.
    assert similarity.jaccard(canonical.uniprot_pairs, reversed_.uniprot_pairs) == 0.0

    representation.check_partner_consistency(records)

    # Same accession in the same role everywhere.
    assert [(r.unp_accession_1, r.unp_accession_2) for r in records] == [(X, Y)] * 3
    # Contact tuples now identical, so equivalent interfaces are identical.
    assert reversed_.uniprot_pairs == canonical.uniprot_pairs == other.uniprot_pairs
    assert similarity.jaccard(canonical.uniprot_pairs, reversed_.uniprot_pairs) == 1.0
    # Author-space joins follow the new roles.
    assert reversed_.author_to_uniprot == {
        ("2bbb", "H", 110, None): (X, 10, 1),
        ("2bbb", "H", 111, None): (X, 11, 1),
        ("2bbb", "L", 220, None): (Y, 20, 2),
        ("2bbb", "L", 221, None): (Y, 21, 2),
    }
    assert reversed_.author_pairs == {
        (("2bbb", "H", 110, None), ("2bbb", "L", 220, None), "hydrogen_bond"),
        (("2bbb", "H", 111, None), ("2bbb", "L", 221, None), "salt_bridge"),
    }
    # Residue identities follow the new roles.
    assert reversed_.residue_identity == canonical.residue_identity == {
        (X, 10, 1): "K", (X, 11, 1): "R", (Y, 20, 2): "D", (Y, 21, 2): "E",
    }
    # Majority records are untouched.
    assert canonical.uniprot_pairs == canonical_pairs_before
    assert canonical.author_to_uniprot == canonical_a2u_before


def test_two_record_tie_keeps_the_first_seen_ordering(make_contact, make_interface):
    """Observed: with a 1:1 split the first record's ordering becomes canonical."""
    records = _check([
        _reversed(make_contact, make_interface, "2bbb"),
        _canonical(make_contact, make_interface, "1aaa", "A", "B"),
    ])
    assert [(r.unp_accession_1, r.unp_accession_2) for r in records] == [(Y, X), (Y, X)]
    assert records[0].uniprot_pairs == records[1].uniprot_pairs


@pytest.mark.parametrize(
    "odd_contact, expected_accessions, expected_pairs",
    [
        # Third accession: matches neither ordering, left as-is with a warning.
        (dict(acc1="P0000Z", acc2=X), ("P0000Z", X),
         {(("P0000Z", 10, 1), (X, 10, 2), "hydrogen_bond")}),
        # No accessions at all: nothing to compare, left as-is.
        (dict(acc1=None, acc2=None), ("", ""), set()),
    ],
    ids=["unrelated-accession-pair", "no-accessions"],
)
def test_records_matching_neither_ordering_are_left_alone(
    odd_contact, expected_accessions, expected_pairs, make_contact, make_interface,
):
    records = _check([
        _canonical(make_contact, make_interface, "1aaa", "A", "B"),
        _canonical(make_contact, make_interface, "3ccc", "C", "D"),
        make_interface("4zzz", [make_contact("A", 10, "B", 20, unp1=10, unp2=10, **odd_contact)]),
    ])
    odd = records[2]
    assert (odd.unp_accession_1, odd.unp_accession_2) == expected_accessions
    assert odd.uniprot_pairs == expected_pairs


# --- Homodimer: the two copies are interchangeable ---------------------------


def _homo_records(make_contact, make_interface, pairs_1, pairs_2):
    """Two homodimer interfaces from (side-1 residue, side-2 residue) tuples."""
    items = [
        make_interface(pdb_id, [
            make_contact("A", r1, "B", r2, acc1=H, acc2=H, unp1=r1, unp2=r2)
            for (r1, r2) in pairs
        ])
        for pdb_id, pairs in (("1aaa", pairs_1), ("2bbb", pairs_2))
    ]
    return _check(items)


def _canonical_pairs(pairs):
    """The contact keys `pairs` should produce: one role, ordered by position."""
    return {
        ((H, min(r1, r2), 1), (H, max(r1, r2), 1), "hydrogen_bond")
        for (r1, r2) in pairs
    }


@pytest.mark.parametrize(
    "pairs_1, pairs_2",
    [
        # The same interface with the chains listed in the opposite order.
        ([(10, 20), (11, 21)], [(20, 10), (21, 11)]),
        # With a self-symmetric contact among them.
        ([(10, 20), (11, 21), (10, 10)], [(20, 10), (21, 11), (10, 10)]),
        # PISA reporting both directions of every contact.
        ([(10, 20), (20, 10)], [(20, 10), (10, 20)]),
    ],
    ids=["flipped-asymmetric", "flipped-partially-symmetric", "both-directions-reported"],
)
def test_homodimer_similarity_is_independent_of_chain_order(
    pairs_1, pairs_2, make_contact, make_interface,
):
    """Which copy PISA lists first is an artefact of the deposition, so two
    recordings of the same interface must compare as identical."""
    rec1, rec2 = _homo_records(make_contact, make_interface, pairs_1, pairs_2)

    assert rec1.unp_accession_1 == rec1.unp_accession_2 == H
    assert rec1.uniprot_pairs == _canonical_pairs(pairs_1)
    assert rec2.uniprot_pairs == _canonical_pairs(pairs_2)
    assert rec1.uniprot_pairs == rec2.uniprot_pairs
    assert similarity.jaccard(rec1.uniprot_pairs, rec2.uniprot_pairs) == 1.0


def test_homodimer_flipped_chain_order_clusters_together(make_contact, make_interface):
    rec1, rec2 = _homo_records(
        make_contact, make_interface, [(10, 20), (11, 21)], [(20, 10), (21, 11)],
    )
    sim = similarity.jaccard_similarity_matrix([rec1.uniprot_pairs, rec2.uniprot_pairs])
    labels = similarity.cluster_interfaces(sim, distance_cut=0.6).flat_assignment
    assert labels[0] == labels[1]


def test_homodimer_keeps_the_deposited_chains_and_roles(make_contact, make_interface):
    """Canonicalisation applies to the comparison key only."""
    rec1, rec2 = _homo_records(
        make_contact, make_interface, [(10, 20)], [(20, 10)],
    )
    # Author-space pairs follow the deposition, so the two differ.
    assert rec1.author_pairs == {
        (("1aaa", "A", 10, None), ("1aaa", "B", 20, None), "hydrogen_bond")}
    assert rec2.author_pairs == {
        (("2bbb", "A", 20, None), ("2bbb", "B", 10, None), "hydrogen_bond")}
    # And so do the author-to-UniProt joins used by the annotation overlap.
    assert rec1.author_to_uniprot == {
        ("1aaa", "A", 10, None): (H, 10, 1), ("1aaa", "B", 20, None): (H, 20, 2)}
    assert rec2.author_to_uniprot == {
        ("2bbb", "A", 20, None): (H, 20, 1), ("2bbb", "B", 10, None): (H, 10, 2)}
