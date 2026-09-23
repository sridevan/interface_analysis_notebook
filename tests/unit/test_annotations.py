"""Annotation mapping: `annotations.overlap_annotations` and the ligand path.

Every annotation is joined by author key, so the dangerous failures are a match
on the wrong chain, a missed insertion code, and an off-interface residue being
reported as at the interface.
"""

from __future__ import annotations

import pandas as pd
import pytest

from pdbe_interfaces import annotations, api, representation

X = "P00001"
Y = "P00002"


def _mutation(chain, resnum, *, ins=None, frm="A", to="G", mtype="Engineered mutation"):
    return {
        "chain_id": chain,
        "author_residue_number": resnum,
        "author_insertion_code": ins,
        "mutation_details": {"from": frm, "to": to, "type": mtype},
    }


def _modification(chain, resnum, ccd="SEP"):
    return {
        "chain_id": chain, "author_residue_number": resnum, "author_insertion_code": None,
        "chem_comp_id": ccd, "chem_comp_name": f"{ccd} name",
    }


def _uniprot_key(row):
    """(unp_acc, unp_seq_id, role) with pandas missing values normalised to None."""
    return tuple(None if pd.isna(v) else v for v in (row.unp_acc, row.unp_seq_id, row.role))


@pytest.fixture
def records(make_contact, make_interface):
    """One interface of entry 1abc.

    Interface residues in author space: A10, A11, A12A, A13 and B20, B21, B22, B23.
    A13/B23 have no UniProt mapping on side 2, so they exist only in author space.
    """
    contacts = [
        make_contact("A", 10, "B", 20, unp1=110, unp2=220),
        make_contact("A", 11, "B", 21, unp1=111, unp2=221),
        make_contact("A", 12, "B", 22, unp1=112, unp2=222, ins1="A"),
        make_contact("A", 13, "B", 23, unp1=113, unp2=None),
    ]
    return representation.check_partner_consistency(
        representation.build_interface_records([make_interface("1abc", contacts)])
    )


def _overlap(records, mutations=None, modifications=None, ligands=None, **kwargs):
    return annotations.overlap_annotations(
        records,
        mutations_response=mutations or {},
        modifications_response=modifications or {},
        ligand_contacts_by_pdb=ligands or {},
        **kwargs,
    )


# --- mutations and modifications ---------------------------------------------


def test_annotations_map_to_uniprot_key_and_role(records):
    overlap = _overlap(
        records,
        mutations={"1abc": [_mutation("A", 10, frm="K", to="A")]},
        modifications={"1abc": [_modification("B", 20, "SEP")]},
    )

    assert len(overlap.mutations) == 1
    mut = overlap.mutations.iloc[0]
    assert (mut.pdb_id, mut.assembly_id, mut.interface_id) == ("1abc", "1", 1)
    assert (mut.auth_asym_id, mut.auth_seq_id, mut.ins_code) == ("A", 10, None)
    assert (mut.unp_acc, mut.unp_seq_id, mut.role) == (X, 110, 1)
    assert mut.mutation_label == "K10A"

    assert len(overlap.modifications) == 1
    mod = overlap.modifications.iloc[0]
    assert (mod.unp_acc, mod.unp_seq_id, mod.role) == (Y, 220, 2)
    assert mod.modification_chem_comp_id == "SEP"


@pytest.mark.parametrize(
    "mutation",
    [
        _mutation("B", 10),            # residue 10 is at the interface on chain A only
        _mutation("A", 99),            # not an interface residue
        _mutation("A", 12),            # A12 is not A12A
        _mutation("A", 12, ins=" "),   # blank code normalises to None: still not A12A
    ],
    ids=["wrong-chain", "off-interface", "missing-insertion-code", "blank-insertion-code"],
)
def test_annotation_on_a_non_interface_residue_is_dropped(records, mutation):
    assert _overlap(records, mutations={"1abc": [mutation]}).mutations.empty


def test_insertion_code_is_part_of_the_join(records):
    df = _overlap(records, mutations={"1abc": [_mutation("A", 12, ins="A")]}).mutations
    assert len(df) == 1
    assert (df.iloc[0].unp_seq_id, df.iloc[0].ins_code) == (112, "A")


def test_annotation_on_partially_mapped_contact(records):
    """A13 maps to UniProt; its contact partner B23 does not.

    The mapped residue resolves normally. The unmapped residue is still at the
    interface in author space, so its row is kept with the
    (None, None, None) sentinel: counted at the interface, not placeable in
    the cross-structure comparison space.
    """
    df = _overlap(records, mutations={"1abc": [_mutation("A", 13), _mutation("B", 23)]}).mutations
    assert len(df) == 2
    by_chain = df.set_index("auth_asym_id")
    assert _uniprot_key(by_chain.loc["A"]) == (X, 113, 1)
    assert _uniprot_key(by_chain.loc["B"]) == (None, None, None)


def test_conflict_mutations_excluded_by_default_and_included_on_request(records):
    response = {"1abc": [
        _mutation("A", 10, mtype="Engineered mutation"),
        _mutation("A", 11, mtype="Conflict"),
        _mutation("B", 20, mtype="Expression tag"),
    ]}
    default = _overlap(records, mutations=response).mutations
    assert default.mutation_type.tolist() == ["Engineered mutation"]

    widened = _overlap(
        records, mutations=response,
        mutation_type_filter=("Engineered mutation", "Conflict"),
    ).mutations
    assert sorted(widened.mutation_type.tolist()) == ["Conflict", "Engineered mutation"]


def test_empty_responses_give_empty_frames(records):
    overlap = _overlap(records)
    assert overlap.mutations.empty and overlap.modifications.empty and overlap.ligands.empty


def test_one_row_per_interface_carrying_the_residue(make_contact, make_interface):
    """An entry with two assemblies that both expose chain A at the interface."""
    items = [
        make_interface("1abc", [make_contact("A", 10, "B", 20)], assembly_id=1),
        make_interface("1abc", [make_contact("A", 10, "C", 20)], assembly_id=2),
    ]
    records = representation.build_interface_records(items)
    df = _overlap(records, mutations={"1abc": [_mutation("A", 10)]}).mutations
    assert sorted(df.assembly_id.tolist()) == ["1", "2"]


def test_mutation_on_reversed_record_lands_on_canonical_role(make_contact, make_interface):
    """Partner reversal must carry the annotation join with it."""
    xa, ya = "P0000X", "P0000Y"
    items = [
        make_interface("1aaa", [make_contact("A", 10, "B", 20, acc1=xa, acc2=ya)]),
        make_interface("3ccc", [make_contact("A", 10, "B", 20, acc1=xa, acc2=ya)]),
        # PISA lists Y first for this entry; X residue 10 sits on chain H.
        make_interface("2bbb", [
            make_contact("L", 220, "H", 110, acc1=ya, acc2=xa, unp1=20, unp2=10),
        ]),
    ]
    records = representation.check_partner_consistency(
        representation.build_interface_records(items)
    )
    df = _overlap(records, mutations={"2bbb": [_mutation("H", 110)]}).mutations

    assert len(df) == 1
    assert (df.iloc[0].unp_acc, df.iloc[0].unp_seq_id, df.iloc[0].role) == (xa, 10, 1)


# --- ligands -----------------------------------------------------------------


BOUND_MOLECULES = [{
    "composition": {"ligands": [
        {"chem_comp_id": "SO4", "chain_id": "A", "author_residue_number": 301,
         "author_insertion_code": "", "molecule_type": "Bound molecule"},
        {"chem_comp_id": "LIG", "chain_id": "A", "author_residue_number": 302,
         "author_insertion_code": " ", "molecule_type": "Bound molecule"},
        {"chem_comp_id": "NAG", "chain_id": "C", "author_residue_number": 1,
         "author_insertion_code": "", "molecule_type": "Carbohydrate-polymer"},
    ]},
}]

LIGAND_INTERACTIONS = [{
    "ligand": {"chem_comp_id": "LIG", "chain_id": "A", "author_residue_number": 302},
    "interactions": [
        # Two atom-level records to the same residue with the same label.
        {"end": {"chain_id": "A", "author_residue_number": 10, "author_insertion_code": ""},
         "interaction_details": ["hbond"]},
        {"end": {"chain_id": "A", "author_residue_number": 10, "author_insertion_code": " "},
         "interaction_details": ["hbond", "vdw"]},
        # Contact residue outside the interface.
        {"end": {"chain_id": "A", "author_residue_number": 99},
         "interaction_details": ["hbond"]},
    ],
}]


def test_filter_bound_molecules_blocklist_and_glycan_flag():
    kept = annotations.filter_bound_molecules(BOUND_MOLECULES, blocklist=frozenset({"SO4"}))
    assert [lig["chem_comp_id"] for lig in kept] == ["LIG", "NAG"]
    assert kept[0] == {"chain_id": "A", "author_residue_number": 302, "ins_code": None,
                       "chem_comp_id": "LIG"}

    kept = annotations.filter_bound_molecules(
        BOUND_MOLECULES, blocklist=frozenset({"SO4"}), drop_carbohydrate_polymers=True,
    )
    assert [lig["chem_comp_id"] for lig in kept] == ["LIG"]


@pytest.fixture
def ligand_contacts(monkeypatch):
    """Collect contacts for 1abc with the two ligand fetchers faked.

    Returns (contacts_by_pdb, requested ligand keys).
    """
    requested: list = []
    monkeypatch.setattr(
        api, "fetch_bound_molecules_many",
        lambda pdb_ids, max_workers=8: {pdb_id: BOUND_MOLECULES for pdb_id in pdb_ids},
    )

    def ligand_interactions_many(keys, max_workers=8):
        requested.extend(keys)
        return [LIGAND_INTERACTIONS if key[2] == 302 else [] for key in keys]

    monkeypatch.setattr(api, "fetch_ligand_interactions_many", ligand_interactions_many)
    contacts = annotations.collect_ligand_contacts_for_entries(
        ["1abc"], blocklist=frozenset({"SO4"}), drop_carbohydrate_polymers=True,
    )
    return contacts, requested


def test_ligand_contacts_collapse_atoms_and_only_fetch_surviving_ligands(ligand_contacts):
    contacts, requested = ligand_contacts

    assert requested == [("1abc", "A", 302)]      # SO4 blocked, NAG dropped
    tuples = {
        (c.interface_residue, c.ligand_chem_comp_id, c.ligand_chain_id,
         c.ligand_author_residue_number, c.contact_type)
        for c in contacts["1abc"]
    }
    assert tuples == {
        (("1abc", "A", 10, None), "LIG", "A", 302, "hbond"),
        (("1abc", "A", 10, None), "LIG", "A", 302, "vdw"),
        (("1abc", "A", 99, None), "LIG", "A", 302, "hbond"),
    }


def test_ligand_contact_maps_to_correct_interface_residue(records, ligand_contacts):
    contacts, _ = ligand_contacts
    df = _overlap(records, ligands=contacts).ligands

    # A99 is a ligand contact but not an interface residue, so it is dropped.
    assert sorted(df.contact_type.tolist()) == ["hbond", "vdw"]
    assert set(zip(df.auth_asym_id, df.auth_seq_id)) == {("A", 10)}
    assert set(zip(df.unp_acc, df.unp_seq_id, df.role)) == {(X, 110, 1)}
    assert set(df.ligand_chem_comp_id) == {"LIG"}


# --- interaction of comparability and residue mapping -------------------------


def test_half_mapped_interface_keeps_residue_mappings_but_is_not_comparable(
    make_contact, make_interface,
):
    """Every contact maps on one side only.

    The interface holds useful residue-level correspondence for annotations
    but has no comparable contact pair, so it must stay out of the
    similarity, clustering and frequency stages.
    """
    items = [
        make_interface("1aaa", [make_contact("A", 10, "B", 20)]),
        make_interface("2bbb", [
            make_contact("A", 50, "B", 106, unp2=None),
            make_contact("A", 51, "B", 107, unp2=None),
        ]),
    ]
    all_records = representation.check_partner_consistency(
        representation.build_interface_records(items)
    )
    half = all_records[1]
    assert half.author_to_uniprot == {
        ("2bbb", "A", 50, None): (X, 50, 1),
        ("2bbb", "A", 51, None): (X, 51, 1),
    }
    assert half.uniprot_pairs == set()

    selection = representation.select_comparable_records(all_records)
    assert [r.pdb_id for r in selection.records] == ["1aaa"]
    assert selection.reasons[half.key] == representation.NO_COMPARABLE_CONTACTS

    # Annotations are joined over all records, so the retained mappings are used.
    df = _overlap(all_records, mutations={"2bbb": [_mutation("A", 50), _mutation("B", 106)]}).mutations
    by_chain = df.set_index("auth_asym_id")
    assert _uniprot_key(by_chain.loc["A"]) == (X, 50, 1)
    assert _uniprot_key(by_chain.loc["B"]) == (None, None, None)
