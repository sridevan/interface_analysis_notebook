"""Scientific regression: the STING homodimer (`11gl` -> PDB-CPX-172174) on frozen data.

Runs the notebook's complete workflow on the recorded PDBe responses in
`tests/fixtures/recorded/PDB-CPX-172174/` (captured 2026-09-23) and checks that
it still reproduces the behaviour documented in `spec/new/demo_notes.md`: four
interface states of 9, 3, 1 and 1, the contacts that separate the two large
states, the ligands at known interfaces, the QC flags on the sparse singleton,
and a consistent export. Expected values were established from the demo notes
and verified against the fixture before being encoded here.

STING is a homodimer, so contacts are compared in the orientation PISA reports
(see the known limitation in final_spec.md). The tranche-1.5 diagnostic found a
maximum orientation delta of about 0.15 for this complex and no clustering
change, so the states below are the current production behaviour and do not
depend on orientation handling.

Cluster ids from `fcluster` are arbitrary; only relationships and size
multisets are asserted.
"""

from __future__ import annotations

import json

import numpy as np

from pdbe_interfaces import outputs, representation

STING = "Q3TBT3"
HB = "hydrogen_bond"
SB = "salt_bridge"


def _contact(pos1, pos2, bond):
    """Typed contact tuple in the production representation, role 1 then role 2."""
    return ((STING, pos1, 1), (STING, pos2, 2), bond)


def _members(run, cluster_id):
    return [r for i, r in enumerate(run.records)
            if int(run.cluster_result.flat_assignment[i]) == cluster_id]


def _fraction_with(records, contact):
    return sum(contact in r.uniprot_pairs for r in records) / len(records)


def test_sting_recorded_workflow_builds_expected_interface_states(sting):
    run = sting

    # Every request was answered from the frozen fixture: 32 recorded calls
    # (2 lookups, 1 interface list, 2 POSTs, 12 bound-molecule GETs, 15 ligand GETs).
    assert len(run.recorded.urls) == run.recorded.metadata["n_requests_recorded"] == 32

    # Identifier resolution, dimer check and partner roles.
    assert run.complex_id == "PDB-CPX-172174"
    assert run.details["oligomeric_state"] == "Homodimer"
    assert run.details["total_chains"] == 2
    assert set(run.partner_map) == {(STING, 1), (STING, 2)}

    # Interface cohort: fixture-specific exact counts.
    assert run.selection.n_interfaces_before == 14
    assert len(run.selection.interfaces) == 14
    assert dict(run.selection.dropped_by_reason) == {}
    assert run.selection.pdb_ids == [
        "11gl", "11gm", "11gn", "4jc5", "4kby", "4kc0",
        "4loj", "4lok", "4lol", "4yp1", "6xnn", "9ltf",
    ]
    assert len(run.all_records) == 14
    assert {(r.unp_accession_1, r.unp_accession_2) for r in run.all_records} == {(STING, STING)}
    assert all(r.n_residues_dropped_no_uniprot == 0 for r in run.all_records)
    assert all(r.n_microheterogeneity_collisions == 0 for r in run.all_records)

    # Comparable-record selection (tranche 1.5): every instance has mapped contacts.
    assert run.comparable.excluded == []
    assert len(run.records) == 14
    assert representation.bond_type_counts(run.records) == {HB: 169, SB: 35}

    # Similarity matrix sanity, end to end.
    assert run.sim_typed.shape == (14, 14)
    assert np.array_equal(run.sim_typed, run.sim_typed.T)
    assert np.array_equal(np.diag(run.sim_typed), np.ones(14))

    # Four states at the default cut of 0.6, sizes 9, 3, 1, 1.
    assert run.cluster_sizes() == [1, 1, 3, 9]
    assert run.cluster_of("11gl") == run.cluster_of("11gm") == run.cluster_of("11gm", "2")
    assert run.cluster_of("11gl") == run.cluster_of("4loj") == run.cluster_of("6xnn")
    assert run.cluster_of("4kby") == run.cluster_of("4kc0") == run.cluster_of("9ltf")
    assert run.cluster_of("11gl") != run.cluster_of("4kby")
    assert len(_members(run, run.cluster_of("4jc5"))) == 1
    assert len(_members(run, run.cluster_of("4lol"))) == 1
    assert len({run.cluster_of(p) for p in ("11gl", "4kby", "4jc5", "4lol")}) == 4

    # Structure table mirrors the comparable cohort and the labels.
    assert len(run.structure_table) == 14
    assert set(run.structure_table.partner_1) == {"Sting1 (role 1)"}
    assert set(run.structure_table.experimental_method) == {"X-ray diffraction"}


def test_sting_recorded_reproduces_known_rewiring(sting):
    run = sting
    large = _members(run, run.cluster_of("11gl"))
    small = _members(run, run.cluster_of("4kby"))
    assert (len(large), len(small)) == (9, 3)

    # Hydrogen bonds present in every member of the large state and absent from
    # the small one, in the production (role 1, role 2) representation.
    large_only = [_contact(232, 209, HB), _contact(232, 260, HB), _contact(209, 233, HB)]
    for contact in large_only:
        assert _fraction_with(large, contact) == 1.0
        assert _fraction_with(small, contact) == 0.0
    # Salt bridge in the opposite direction.
    small_only = _contact(273, 156, SB)
    assert _fraction_with(small, small_only) == 1.0
    assert _fraction_with(large, small_only) == 0.0

    # The rewiring table labels them the same way.
    table = outputs.compare_cluster_contacts(
        run.records, run.cluster_result,
        cluster_a=run.cluster_of("11gl"), cluster_b=run.cluster_of("4kby"),
        partner_map=run.partner_map, top_n=None,
    )
    direction = table.set_index("contact")["contact_direction"]
    assert direction["Sting1:A232-Sting1:D209 hydrogen_bond"] == "higher in A"
    assert direction["Sting1:A232-Sting1:Y260 hydrogen_bond"] == "higher in A"
    assert direction["Sting1:D209-Sting1:G233 hydrogen_bond"] == "higher in A"
    assert direction["Sting1:D273-Sting1:H156 salt_bridge"] == "higher in B"
    assert "shared core" not in set(direction)   # the two states share no core contact


def test_sting_recorded_maps_known_interface_ligands(sting):
    run = sting
    ligands = run.overlap.ligands
    assert run.overlap.mutations.empty          # no engineered mutation at any interface
    assert run.overlap.modifications.empty      # endpoint returned 404 for all entries

    def at(pdb_id, ccd):
        rows = ligands[(ligands.pdb_id == pdb_id) & (ligands.ligand_chem_comp_id == ccd)]
        return {(r.auth_asym_id, r.auth_seq_id, r.unp_acc, r.unp_seq_id, r.role) for r in rows.itertuples()}

    # 1YE binds at the 4lol interface on Ser161 of both copies.
    assert at("4lol", "1YE") == {("A", 161, STING, 161, 1), ("B", 161, STING, 161, 2)}
    assert set(ligands[ligands.pdb_id == "4lol"].ligand_chem_comp_id) == {"1YE"}

    # Cyclic di-GMP (C2E) contacts the 4kby interface, including Arg237 of both copies.
    c2e = at("4kby", "C2E")
    assert ("A", 237, STING, 237, 1) in c2e and ("B", 237, STING, 237, 2) in c2e
    assert {row[2] for row in c2e} == {STING}

    # Ligands sit in the states the demo notes describe.
    lig_by_pdb = ligands.groupby("pdb_id")["ligand_chem_comp_id"].agg(set)
    assert lig_by_pdb["4lok"] == {"1YD"} and run.cluster_of("4lok") == run.cluster_of("11gl")
    assert lig_by_pdb["9ltf"] == {"A1ELY"} and run.cluster_of("9ltf") == run.cluster_of("4kby")
    assert "4jc5" not in lig_by_pdb.index


def test_sting_recorded_reports_expected_qc_and_frequencies(sting):
    run = sting
    report = run.report.set_index("cluster_id")
    assert report.cluster_size.tolist() == [9, 3, 1, 1]      # largest state first

    # 4jc5: singleton with only four residue pairs, flagged on both counts.
    row = report.loc[run.cluster_of("4jc5")]
    assert row.member_pdb_ids == "4jc5"
    assert row.residue_pair_count_median == 4
    assert "singleton" in row.qc_warnings
    assert "sparse" in row.qc_warnings
    # 4lol: singleton with seven pairs, so only the singleton flag.
    row = report.loc[run.cluster_of("4lol")]
    assert "singleton" in row.qc_warnings and "sparse" not in row.qc_warnings
    assert "1YE" in row.cluster_ligands
    # The large state carries no QC warning.
    assert report.loc[run.cluster_of("11gl")].qc_warnings == ""

    # Aggregation runs over the 14 comparable instances; nothing reaches 0.8.
    freq = outputs.interface_frequency_summary(run.records, partner_map=run.partner_map)
    assert freq["n_interfaces"] == 14
    top = freq["partner_1_residues"].iloc[0]
    assert (top.position, top.n_interfaces) == (209, 10)
    assert top.fraction == round(10 / 14, 3)
    assert outputs.conserved_residues(run.records, threshold=0.8) == set()
    assert outputs.conserved_residues(run.records, threshold=10 / 14) >= {(STING, 209, 1), (STING, 235, 1)}
    # The large-state contact is present in exactly its nine members.
    pairs = freq["pairs"].set_index("contact")
    assert pairs.loc["Sting1:A232 – Sting1:D209", "n_interfaces"] == 9


def test_sting_recorded_export_is_consistent(sting):
    run = sting
    assert run.output_path.exists()
    written = json.loads(run.output_path.read_text())
    assert written["metadata"]["complex_id"] == "PDB-CPX-172174"
    assert written["metadata"]["oligomeric_state"] == "Homodimer"
    assert written["metadata"]["n_interfaces"] == 14
    assert written["residue_frequencies"] == run.export["residue_frequencies"]

    partners = {(p["role"], p["unp_accession"], p["gene_name"]) for p in written["partners"]}
    assert partners == {(1, STING, "Sting1"), (2, STING, "Sting1")}

    residues = {(r["role"], r["unp_residue_number"]): r for r in written["residue_frequencies"]}
    assert residues[(1, 209)]["n_interfaces"] == 10
    assert residues[(1, 209)]["frequency"] == round(10 / 14, 4)
    assert residues[(1, 209)]["conservation_level"] == "medium"
    assert residues[(1, 209)]["unp_residue_label"] == "D209"
    assert all(0 < r["n_interfaces"] <= 14 for r in written["residue_frequencies"])

    contacts = {
        (c["partner_1"]["unp_residue_number"], c["partner_2"]["unp_residue_number"], c["bond_type"]): c
        for c in written["contact_frequencies"]
    }
    assert contacts[(232, 209, HB)]["n_interfaces"] == 9
    assert contacts[(273, 156, SB)]["n_interfaces"] == 3
    assert all(c["partner_1"]["role"] == 1 and c["partner_2"]["role"] == 2
               for c in written["contact_frequencies"])
