"""Homodimer contact symmetry: `representation.canonical_uniprot_pair`.

The two copies of a homodimer are interchangeable, so a contact between them
is an unordered pair of UniProt residues. The canonicalisation happens once,
where `uniprot_pairs` is built, and every analysis path inherits it; these
tests check the rule itself and then each surface that reads it.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pdbe_interfaces import outputs, representation, similarity
from pdbe_interfaces.similarity import ClusterResult

H = "P0HOMO"
X = "P0000X"
Y = "P0000Y"
HB = "hydrogen_bond"
SB = "salt_bridge"
ROLE = representation.HOMODIMER_ROLE


def _homo(make_contact, make_interface, pdb_id, pairs, chains=("A", "B")):
    """One homodimer interface from (side-1 residue, side-2 residue[, bond])."""
    c1, c2 = chains
    return make_interface(pdb_id, [
        make_contact(c1, p[0], c2, p[1], acc1=H, acc2=H, unp1=p[0], unp2=p[1],
                     bond=(p[2] if len(p) > 2 else HB), aa1="A", aa2="A")
        for p in pairs
    ])


def _record(make_contact, make_interface, pdb_id, pairs, chains=("A", "B")):
    return representation.build_interface_records(
        [_homo(make_contact, make_interface, pdb_id, pairs, chains)]
    )[0]


def _key(pos1, pos2, bond=HB):
    return ((H, pos1, ROLE), (H, pos2, ROLE), bond)


# --- the rule ----------------------------------------------------------------


def test_reciprocal_homodimer_contacts_give_one_key():
    forward = representation.canonical_uniprot_pair((H, 209, 1), (H, 232, 2), HB, True)
    reverse = representation.canonical_uniprot_pair((H, 232, 1), (H, 209, 2), HB, True)
    assert forward == reverse == _key(209, 232)


def test_homodimer_detection_uses_the_uniprot_accession():
    assert representation.is_homodimer(H, H) is True
    assert representation.is_homodimer(X, Y) is False
    assert representation.is_homodimer("", "") is False


def test_heterodimer_contact_orientation_is_untouched():
    pair = representation.canonical_uniprot_pair((X, 42, 1), (Y, 266, 2), HB, False)
    assert pair == ((X, 42, 1), (Y, 266, 2), HB)
    # The reverse listing stays a different key: the partners are different
    # proteins, so which side a residue is on is part of the contact.
    other = representation.canonical_uniprot_pair((Y, 266, 1), (X, 42, 2), HB, False)
    assert other != pair


# --- the fingerprint ---------------------------------------------------------


def test_reciprocal_duplicates_within_one_interface_are_counted_once(
    make_contact, make_interface,
):
    rec = _record(make_contact, make_interface, "1aaa", [(209, 232), (232, 209)])
    assert rec.uniprot_pairs == {_key(209, 232)}


def test_a_self_pair_stays_one_contact(make_contact, make_interface):
    rec = _record(make_contact, make_interface, "1aaa", [(238, 238)])
    assert rec.uniprot_pairs == {_key(238, 238)}


def test_the_same_residue_pair_with_two_bond_types_stays_two_typed_contacts(
    make_contact, make_interface,
):
    rec = _record(make_contact, make_interface, "1aaa",
                  [(209, 231, HB), (231, 209, SB)])
    assert rec.uniprot_pairs == {_key(209, 231, HB), _key(209, 231, SB)}
    # The untyped view drops the bond type, as it always did, leaving one pair.
    assert similarity.untyped(rec.uniprot_pairs) == {((H, 209, ROLE), (H, 231, ROLE))}


def test_role_reversed_recordings_of_one_interface_are_equivalent(
    make_contact, make_interface,
):
    a = _record(make_contact, make_interface, "1aaa", [(10, 20), (11, 21)])
    b = _record(make_contact, make_interface, "2bbb", [(20, 10), (21, 11)], chains=("B", "A"))
    assert a.uniprot_pairs == b.uniprot_pairs
    assert similarity.jaccard(a.uniprot_pairs, b.uniprot_pairs) == 1.0


def test_the_deposited_structural_mapping_is_preserved(make_contact, make_interface):
    """Annotations and the 3D views read the author side, which must not move."""
    rec = _record(make_contact, make_interface, "1aaa", [(232, 209)])
    assert rec.author_pairs == {
        (("1aaa", "A", 232, None), ("1aaa", "B", 209, None), HB)}
    # Author residues keep the role of the copy they sit on, so an annotation
    # on chain B is still attributed to chain B.
    assert rec.author_to_uniprot == {
        ("1aaa", "A", 232, None): (H, 232, 1),
        ("1aaa", "B", 209, None): (H, 209, 2),
    }
    # The canonical key the analysis uses is the unordered one.
    assert rec.uniprot_pairs == {_key(209, 232)}


# --- the surfaces that read the fingerprint ----------------------------------


@pytest.fixture
def homo_run(make_contact, make_interface):
    """Four homodimer interfaces: two listing the chains each way, in two groups."""
    records = [
        _record(make_contact, make_interface, "1aaa", [(10, 20), (30, 40)]),
        _record(make_contact, make_interface, "2bbb", [(20, 10), (40, 30)]),
        _record(make_contact, make_interface, "3ccc", [(10, 20), (50, 60)]),
        _record(make_contact, make_interface, "4ddd", [(60, 50), (20, 10)]),
    ]
    sim = similarity.jaccard_similarity_matrix([r.uniprot_pairs for r in records])
    return records, similarity.cluster_interfaces(sim, distance_cut=0.6)


def test_contact_frequencies_have_no_reciprocal_duplicate_rows(homo_run):
    records, _ = homo_run
    pairs = outputs.interface_frequency_summary(records)["pairs"]
    seen = set(zip(pairs.partner_1, pairs.partner_2))
    assert not {(a, b) for (a, b) in seen if (b, a) in seen and a != b}
    # The contact listed both ways across the four interfaces is one row of 4.
    row = pairs[(pairs.partner_1.str.endswith(":A10")) & (pairs.partner_2.str.endswith(":A20"))]
    assert len(row) == 1 and int(row.iloc[0].n_interfaces) == 4


def test_exported_contacts_have_no_reciprocal_duplicates(homo_run, tmp_path):
    records, _ = homo_run
    export = outputs.export_interface_frequency_json(
        records, complex_id="T", config=SimpleNamespace(output_dir=str(tmp_path)))
    keys = {(c["partner_1"]["unp_residue_number"], c["partner_2"]["unp_residue_number"],
             c["bond_type"]) for c in export["contact_frequencies"]}
    assert not {(a, b, t) for (a, b, t) in keys if (b, a, t) in keys and a != b}
    assert (10, 20, HB) in keys and (20, 10, HB) not in keys
    # The collapse is declared, and neither side claims a partner role.
    assert export["metadata"]["contact_symmetry"] == "unordered"
    assert export["metadata"]["homodimer_symmetry_collapsed"] is True
    assert all(c["partner_1"]["role"] is None and c["partner_2"]["role"] is None
               for c in export["contact_frequencies"])
    # The partners array still describes the deposited complex, both copies.
    assert {p["role"] for p in export["partners"]} == {1, 2}


def test_n_distinct_fingerprints_uses_the_canonical_fingerprint(homo_run):
    records, cluster_result = homo_run
    report = outputs.cluster_interpretation_report(
        records, cluster_result, outputs.AnnotationOverlap()).set_index("cluster_id")
    # 1aaa and 2bbb are the same interface listed each way: one fingerprint.
    group = report.loc[int(cluster_result.flat_assignment[0])]
    assert group.cluster_size == 2 and group.n_distinct_fingerprints == 1


def test_a_homodimer_residue_is_counted_once_per_interface(
    make_contact, make_interface, tmp_path,
):
    """Residue 10 carries a contact on both copies of the protein. It is one
    UniProt residue, so it counts once, not once per copy."""
    records = [_record(make_contact, make_interface, "1aaa", [(10, 20), (30, 10)])]
    freq = outputs.interface_frequency_summary(records)
    residues = {
        (row.accession, int(row.position)): int(row.n_interfaces)
        for table in ("partner_1_residues", "partner_2_residues")
        for row in freq[table].itertuples()
    }
    assert residues[(H, 10)] == 1
    assert outputs.conserved_residues(records, threshold=1.0) == {
        (H, 10, ROLE), (H, 20, ROLE), (H, 30, ROLE)}
    export = outputs.export_interface_frequency_json(
        records, complex_id="T", config=SimpleNamespace(output_dir=str(tmp_path)))
    assert [r["unp_residue_number"] for r in export["residue_frequencies"]].count(10) == 1
    # The row describes a UniProt position, not one copy of the protein.
    assert {r["role"] for r in export["residue_frequencies"]} == {None}


# --- heterodimer control -----------------------------------------------------


def test_heterodimer_fingerprints_and_similarity_are_unchanged(
    make_contact, make_interface,
):
    def hetero(pdb_id, pairs):
        return representation.build_interface_records([make_interface(pdb_id, [
            make_contact("A", p[0], "B", p[1], acc1=X, acc2=Y, unp1=p[0], unp2=p[1])
            for p in pairs
        ])])[0]

    a = hetero("1aaa", [(42, 266), (44, 295)])
    b = hetero("2bbb", [(266, 42), (295, 44)])
    # Partner order is meaningful, so these keep the orientation they were given.
    assert a.uniprot_pairs == {((X, 42, 1), (Y, 266, 2), HB), ((X, 44, 1), (Y, 295, 2), HB)}
    assert b.uniprot_pairs == {((X, 266, 1), (Y, 42, 2), HB), ((X, 295, 1), (Y, 44, 2), HB)}
    assert similarity.jaccard(a.uniprot_pairs, b.uniprot_pairs) == 0.0
    # And a heterodimer residue keeps its side.
    result = ClusterResult(linkage=np.empty((0, 4)), flat_assignment=np.array([1, 1]),
                           distance_cut=0.6)
    report = outputs.cluster_interpretation_report([a, a], result, outputs.AnnotationOverlap())
    assert report.n_distinct_fingerprints.tolist() == [1]


def test_heterodimer_export_keeps_partner_roles(make_contact, make_interface, tmp_path):
    records = representation.build_interface_records([make_interface("1aaa", [
        make_contact("A", 42, "B", 266, acc1=X, acc2=Y, unp1=42, unp2=266),
    ])])
    export = outputs.export_interface_frequency_json(
        records, complex_id="T", config=SimpleNamespace(output_dir=str(tmp_path)))
    assert export["metadata"]["contact_symmetry"] == "ordered"
    assert export["metadata"]["homodimer_symmetry_collapsed"] is False
    assert {r["role"] for r in export["residue_frequencies"]} == {1, 2}
    contact = export["contact_frequencies"][0]
    assert (contact["partner_1"]["role"], contact["partner_2"]["role"]) == (1, 2)
    assert (contact["partner_1"]["unp_accession"], contact["partner_2"]["unp_accession"]) == (X, Y)
