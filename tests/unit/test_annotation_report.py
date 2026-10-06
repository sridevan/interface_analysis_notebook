"""Group-level annotation summary: `outputs.annotation_report`.

Counting is once per interface instance, with supporting PDB entries counted
separately. Ligand frequencies are taken over preferred-assembly interfaces
only, because PDBe computes Arpeggio protein-ligand interactions for the
preferred assembly; mutations and modifications come from per-entry endpoints
with no assembly dimension, so their denominator is the whole group.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pdbe_interfaces import outputs
from pdbe_interfaces.annotations import AnnotationOverlap
from pdbe_interfaces.representation import InterfaceRecord
from pdbe_interfaces.similarity import ClusterResult

P1, P2 = "P00001", "P00002"


def _record(pdb_id, assembly_id="1", area=None):
    info = {} if area is None else {"interface_area": area}
    return InterfaceRecord(
        pdb_id=pdb_id, assembly_id=assembly_id, interface_id=1, interface_info=info,
        unp_accession_1=P1, unp_accession_2=P2,
    )


def _cluster(n, assignment=None):
    flat = np.array(assignment if assignment is not None else [1] * n)
    return ClusterResult(linkage=np.empty((0, 4)), flat_assignment=flat, distance_cut=0.6)


def _ligand_rows(*triples, ccd="ATP", contact="hbond"):
    """One row per (pdb_id, assembly_id, interface residue)."""
    return pd.DataFrame([
        {"pdb_id": p, "assembly_id": a, "interface_id": 1, "auth_seq_id": res,
         "ligand_chem_comp_id": ccd, "contact_type": contact}
        for (p, a, res) in triples
    ])


def _meta(*entries):
    """(pdb_id, assembly_id, preferred) triples -> assembly metadata."""
    return {(p, a): {"experimental_method": "X-ray diffraction", "resolution": 2.0,
                     "preferred_assembly": pref} for (p, a, pref) in entries}


# --- interface area ----------------------------------------------------------


def test_median_interface_area_is_reported_per_group():
    records = [_record("1aaa", area=100.0), _record("2bbb", area=300.0),
               _record("3ccc", area=200.0)]
    report = outputs.cluster_interpretation_report(records, _cluster(3), AnnotationOverlap())
    assert report.interface_area_median.tolist() == [200.0]
    assert report.interface_area_min.tolist() == [100.0]
    assert report.interface_area_max.tolist() == [300.0]


def test_missing_interface_areas_are_skipped_not_zeroed():
    """A missing area shrinks the sample rather than pulling the median down,
    which is the behaviour the QC logic already relies on."""
    records = [_record("1aaa", area=100.0), _record("2bbb", area=None),
               _record("3ccc", area=200.0)]
    report = outputs.cluster_interpretation_report(records, _cluster(3), AnnotationOverlap())
    assert report.interface_area_median.tolist() == [150.0]
    assert "n=2" in report.interface_area_range.iloc[0]

    # No area at all leaves the statistics empty rather than zero.
    none_only = [_record("1aaa"), _record("2bbb")]
    empty = outputs.cluster_interpretation_report(none_only, _cluster(2), AnnotationOverlap())
    assert empty.interface_area_median.iloc[0] is None
    assert empty.interface_area_range.iloc[0] == ""


# --- counting unit -----------------------------------------------------------


def test_repeated_ligand_rows_in_one_interface_count_once():
    records = [_record("1aaa"), _record("2bbb")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(
        ("1aaa", "1", 10), ("1aaa", "1", 11), ("1aaa", "1", 12),   # one interface, three residues
    ))
    table = outputs.annotation_report(records, _cluster(2), overlap,
                                      assembly_metadata=_meta(("1aaa", "1", True), ("2bbb", "1", True)))
    row = table.iloc[0]
    assert row.n_interfaces_with_annotation == 1
    assert (row.n_eligible_interfaces, row.group_size) == (2, 2)
    assert row.frequency_among_eligible == pytest.approx(0.5)


def test_two_interfaces_from_one_entry_are_two_interfaces_but_one_entry():
    records = [_record("1aaa", "1"), _record("1aaa", "2"), _record("2bbb", "1")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "1", 10), ("1aaa", "2", 10)))
    table = outputs.annotation_report(
        records, _cluster(3), overlap,
        assembly_metadata=_meta(("1aaa", "1", True), ("1aaa", "2", True), ("2bbb", "1", True)))
    row = table.iloc[0]
    assert row.n_interfaces_with_annotation == 2
    assert row.n_pdb_entries_with_annotation == 1
    assert (row.n_eligible_interfaces, row.n_eligible_pdb_entries) == (3, 2)


# --- ligand eligibility ------------------------------------------------------


def test_non_preferred_interfaces_leave_the_ligand_denominator():
    """Three interfaces, one on a non-preferred assembly: the denominator is 2,
    and the non-preferred interface is not counted as ligand-negative."""
    records = [_record("1aaa", "1"), _record("2bbb", "1"), _record("2bbb", "2")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "1", 10)))
    meta = _meta(("1aaa", "1", True), ("2bbb", "1", False), ("2bbb", "2", True))
    row = outputs.annotation_report(records, _cluster(3), overlap, assembly_metadata=meta).iloc[0]

    assert (row.n_eligible_interfaces, row.group_size) == (2, 3)
    assert row.frequency_among_eligible == pytest.approx(0.5)      # 1/2, not 1/3
    assert row.eligibility_basis == "preferred_assembly"
    # Eligibility and occurrence are separate quantities.
    assert row.n_interfaces_with_annotation == 1 and row.n_eligible_interfaces == 2


def test_a_ligand_row_on_a_non_preferred_interface_does_not_count():
    records = [_record("1aaa", "1"), _record("1aaa", "2")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "2", 10)))   # non-preferred
    meta = _meta(("1aaa", "1", True), ("1aaa", "2", False))
    table = outputs.annotation_report(records, _cluster(2), overlap, assembly_metadata=meta)
    assert table.empty      # nothing eligible carried the ligand


def test_without_assembly_metadata_no_interface_is_ligand_eligible():
    """Eligibility is never inferred from whether ligand rows came back."""
    records = [_record("1aaa"), _record("2bbb")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "1", 10)))
    assert outputs.annotation_report(records, _cluster(2), overlap).empty


# --- mutations and modifications --------------------------------------------


def test_mutation_and_modification_denominators_ignore_the_ligand_rule():
    records = [_record("1aaa", "1"), _record("1aaa", "2"), _record("2bbb", "1")]
    overlap = AnnotationOverlap(
        mutations=pd.DataFrame([
            {"pdb_id": "1aaa", "assembly_id": "2", "interface_id": 1,
             "mutation_label": "D237N", "mutation_type": "Engineered mutation"},
        ]),
        modifications=pd.DataFrame([
            {"pdb_id": "1aaa", "assembly_id": "2", "interface_id": 1,
             "modification_chem_comp_id": "SEP",
             "modification_chem_comp_name": "PHOSPHOSERINE"},
        ]),
    )
    meta = _meta(("1aaa", "1", True), ("1aaa", "2", False), ("2bbb", "1", True))
    table = outputs.annotation_report(records, _cluster(3), overlap, assembly_metadata=meta)

    by_type = table.set_index("annotation_type")
    for kind in ("mutation", "modification"):
        row = by_type.loc[kind]
        # The annotation sits on the non-preferred assembly and still counts.
        assert row.n_interfaces_with_annotation == 1
        assert (row.n_eligible_interfaces, row.group_size) == (3, 3)
        assert row.eligibility_basis == "all_interfaces"
    assert by_type.loc["mutation"].annotation == "D237N"
    assert by_type.loc["mutation"].annotation_detail == "Engineered mutation"
    assert by_type.loc["modification"].annotation == "SEP"
    assert by_type.loc["modification"].annotation_detail == "PHOSPHOSERINE"


def test_distinct_annotations_are_not_collapsed():
    records = [_record("1aaa"), _record("2bbb")]
    overlap = AnnotationOverlap(ligands=pd.concat([
        _ligand_rows(("1aaa", "1", 10), ccd="ATP"),
        _ligand_rows(("2bbb", "1", 10), ccd="GTP"),
    ], ignore_index=True))
    table = outputs.annotation_report(records, _cluster(2), overlap,
                                      assembly_metadata=_meta(("1aaa", "1", True), ("2bbb", "1", True)))
    assert sorted(table.annotation) == ["ATP", "GTP"]
    assert table.n_interfaces_with_annotation.tolist() == [1, 1]


# --- the report does not disturb the grouping --------------------------------


def test_annotation_report_does_not_change_group_membership():
    records = [_record("1aaa"), _record("2bbb"), _record("3ccc")]
    result = _cluster(3, [1, 1, 2])
    before = result.flat_assignment.tolist()
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "1", 10)))
    meta = _meta(("1aaa", "1", True), ("2bbb", "1", True), ("3ccc", "1", True))

    table = outputs.annotation_report(records, result, overlap, assembly_metadata=meta)
    report = outputs.cluster_interpretation_report(records, result, overlap,
                                                   assembly_metadata=meta)
    assert result.flat_assignment.tolist() == before
    assert report.cluster_size.tolist() == [2, 1]
    assert set(table.cluster_id) == {1}


# --- presentation ------------------------------------------------------------


def test_formatted_view_reports_occurrence_and_eligibility_separately():
    records = [_record("1aaa", "1"), _record("2bbb", "1"), _record("2bbb", "2")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "1", 10)))
    meta = _meta(("1aaa", "1", True), ("2bbb", "1", False), ("2bbb", "2", True))
    table = outputs.annotation_report(records, _cluster(3), overlap, assembly_metadata=meta)
    view = outputs.format_annotation_report(table)

    assert list(view.columns) == ["Interaction group", "Type", "Annotation",
                                  "Interfaces with annotation",
                                  "PDB entries with annotation", "Eligible"]
    row = view.iloc[0]
    assert (row.Type, row.Annotation) == ("Ligand", "ATP")
    assert row["Interfaces with annotation"] == "1/2 (50%)"
    assert row.Eligible == "2/3 interfaces"
    assert table.equals(outputs.annotation_report(records, _cluster(3), overlap,
                                                  assembly_metadata=meta))


def test_missing_annotation_types_are_named_not_tabulated():
    records = [_record("1aaa")]
    overlap = AnnotationOverlap(ligands=_ligand_rows(("1aaa", "1", 10)))
    table = outputs.annotation_report(records, _cluster(1), overlap,
                                      assembly_metadata=_meta(("1aaa", "1", True)))
    assert outputs.describe_missing_annotation_types(table) == (
        "No mutations or modifications were observed at any interface in this dataset.")
    assert outputs.describe_missing_annotation_types(
        outputs.annotation_report(records, _cluster(1), AnnotationOverlap())
    ).startswith("No ligands, mutations or modifications")


def test_group_report_display_renames_without_touching_the_dataframe():
    records = [_record("1aaa", area=100.0), _record("2bbb", area=300.0)]
    report = outputs.cluster_interpretation_report(records, _cluster(2), AnnotationOverlap())
    before = report.copy()
    view = outputs.format_group_report(report)

    assert list(view.columns) == [
        "Interaction group", "Group size", "PDB entries", "Distinct fingerprints",
        "Member interfaces", "Methods", "Resolution", "Median area (Å²)", "Quality flags",
    ]
    assert view["Interaction group"].tolist() == report["cluster_id"].tolist()
    assert view["Group size"].tolist() == report["cluster_size"].tolist()
    assert view["Median area (Å²)"].tolist() == report["interface_area_median"].tolist()
    # The group-size-denominated annotation columns are not shown beside the
    # eligibility-aware annotation report.
    assert not {"cluster_ligands", "cluster_mutations", "cluster_modifications"} & set(view.columns)
    assert report.equals(before)
    # The dataframe keeps its own field names for downstream use.
    assert "cluster_id" in report.columns and "cluster_size" in report.columns
