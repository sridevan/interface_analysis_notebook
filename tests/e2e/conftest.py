"""End-to-end fixtures: run the notebook's workflow on frozen, recorded PDBe responses.

`RecordedPDBe` stands in for `requests.Session.request` and serves the JSON captured
under `tests/fixtures/recorded/<complex_id>/`, so the real `pdbe_interfaces.api`
layer (URL building, 404 policy, batched POST, thread-pool fan-out) runs unchanged.
Any URL that was not recorded raises, which keeps the test offline and makes an
incomplete fixture fail loudly. The root `block_network` fixture stays active for
everything outside the recorded run.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import requests

from pdbe_interfaces import Config, annotations, api, outputs, representation, similarity

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "recorded"


class RecordedResponse:
    def __init__(self, status_code: int, body):
        self.status_code = status_code
        self._body = body

    def json(self):
        return self._body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}", response=self)


class RecordedPDBe:
    """Route PDBe API v2 URLs to the recorded responses of one complex."""

    def __init__(self, fixture_dir: str):
        d = FIXTURES / fixture_dir
        load = lambda name: json.loads((d / name).read_text())  # noqa: E731
        self.fixture_dir = fixture_dir
        self.metadata = load("metadata.json")
        self.complex_id = self.metadata["complex_id"]
        # The entry-id lookup file is named after the starting PDB entry.
        [entry_file] = [f for f in d.glob("complex_details_*.json")]
        self.entry_lookup = load(entry_file.name)
        self.complex_details = load("complex_details.json")
        self.interface_interactions = load("interface_interactions.json")
        self.mutations_post = load("mutations_post.json")
        self.modifications_post = load("modifications_post.json")
        self.bound_molecules = load("bound_molecules.json")
        self.ligand_interactions = load("ligand_interactions.json")
        self.urls: list[str] = []

    def __call__(self, method: str, url: str, **kwargs) -> RecordedResponse:
        self.urls.append(url)
        rec = self._route(method, url, kwargs)
        return RecordedResponse(rec["status"], rec["body"])

    def _route(self, method: str, url: str, kwargs: dict) -> dict:
        if url == self.entry_lookup["url"]:
            return self.entry_lookup
        if url == self.complex_details["url"]:
            return self.complex_details
        if url == self.interface_interactions["url"]:
            return self.interface_interactions
        for post in (self.mutations_post, self.modifications_post):
            if url == post["url"]:
                assert method == "POST"
                sent = json.loads(kwargs["data"])
                if sent != post["request_body"]:
                    raise RuntimeError(f"Unrecorded POST body for {url}: {sent!r}")
                return post
        m = re.search(r"/pdb/bound_molecules/([^/?]+)$", url)
        if m and m.group(1) in self.bound_molecules:
            return self.bound_molecules[m.group(1)]
        m = re.search(r"/pdb/bound_ligand_interactions/([^?]+)\?preserve_case=false$", url)
        if m and m.group(1) in self.ligand_interactions:
            return self.ligand_interactions[m.group(1)]
        raise RuntimeError(f"Unrecorded request in offline e2e test: {method} {url}")


@dataclass
class WorkflowRun:
    """Everything the notebook computes, in the notebook's order."""

    config: Config
    complex_id: str
    details: dict
    partner_map: dict
    assembly_metadata: dict
    selection: representation.InterfaceSelection
    mutations_response: dict
    modifications_response: dict
    ligand_contacts_by_pdb: dict
    all_records: list
    comparable: representation.ComparableSelection
    records: list
    sim_typed: np.ndarray
    cluster_result: similarity.ClusterResult
    overlap: annotations.AnnotationOverlap
    structure_table: pd.DataFrame
    report: pd.DataFrame
    export: dict
    output_path: Path
    recorded: RecordedPDBe

    def cluster_of(self, pdb_id: str, assembly_id: str = "1", interface_id: int = 1) -> int:
        """Cluster label of one interface instance; labels are arbitrary integers."""
        for i, r in enumerate(self.records):
            if r.key == (pdb_id, assembly_id, interface_id):
                return int(self.cluster_result.flat_assignment[i])
        raise KeyError((pdb_id, assembly_id, interface_id))

    def cluster_sizes(self) -> list[int]:
        return sorted(np.bincount(self.cluster_result.flat_assignment)[1:].tolist())


def run_workflow(identifier: str, fixture_dir: str, output_dir: Path) -> WorkflowRun:
    """The notebook's Phase 1 to Phase 7 call sequence, unchanged, on recorded responses."""
    cfg = Config(identifier=identifier, output_dir=str(output_dir))
    recorded = RecordedPDBe(fixture_dir)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(requests.Session, "request", recorded)

        # Phase 1
        resolved_id, details = api.resolve_complex_id(cfg.identifier, require_dimer=True)
        partner_map = outputs.build_partner_map(details)
        assembly_metadata = outputs.extract_assembly_metadata(details)
        selection = representation.select_interfaces(
            api.fetch_interface_interactions(resolved_id), details,
            max_resolution=cfg.max_resolution, max_entries=cfg.max_entries,
        )
        pdb_ids = selection.pdb_ids
        mutations_response = api.fetch_mutations(pdb_ids)
        modifications_response = api.fetch_modifications(pdb_ids)
        ligand_contacts_by_pdb = annotations.collect_ligand_contacts_for_entries(
            pdb_ids, blocklist=cfg.ligand_blocklist,
            drop_carbohydrate_polymers=cfg.drop_carbohydrate_polymers,
            max_workers=cfg.max_workers,
        )

    # Phases 2 to 7 are pure computation; the recorded session is no longer needed.
    all_records = representation.check_partner_consistency(
        representation.build_interface_records(selection.interfaces)
    )
    comparable = representation.select_comparable_records(all_records)
    records = comparable.records
    sim_typed = similarity.jaccard_similarity_matrix([r.uniprot_pairs for r in records])
    cluster_result = similarity.cluster_interfaces(sim_typed, distance_cut=cfg.cluster_distance_cut)
    overlap = annotations.overlap_annotations(
        all_records, mutations_response=mutations_response,
        modifications_response=modifications_response,
        ligand_contacts_by_pdb=ligand_contacts_by_pdb,
        mutation_type_filter=cfg.mutation_type_filter,
    )
    structure_table = outputs.build_structure_table(
        records, cluster_result, overlap, partner_map, assembly_metadata=assembly_metadata,
    )
    report = outputs.cluster_interpretation_report(
        records, cluster_result, overlap, assembly_metadata=assembly_metadata,
        partner_map=partner_map, conservation_threshold=cfg.conservation_threshold,
    )
    export = outputs.export_interface_frequency_json(
        records, partner_map=partner_map, complex_id=resolved_id,
        complex_name=details.get("name"), oligomeric_state=details.get("oligomeric_state"),
        config=cfg, typed_contacts=True, complex_details=details,
    )
    return WorkflowRun(
        config=cfg, complex_id=resolved_id, details=details, partner_map=partner_map,
        assembly_metadata=assembly_metadata, selection=selection,
        mutations_response=mutations_response, modifications_response=modifications_response,
        ligand_contacts_by_pdb=ligand_contacts_by_pdb, all_records=all_records,
        comparable=comparable, records=records, sim_typed=sim_typed,
        cluster_result=cluster_result, overlap=overlap, structure_table=structure_table,
        report=report, export=export, output_path=output_dir / f"{resolved_id}.json",
        recorded=recorded,
    )


@pytest.fixture(scope="module")
def sting(tmp_path_factory) -> WorkflowRun:
    """The STING homodimer run (`11gl` -> PDB-CPX-172174) on the frozen fixture."""
    return run_workflow("11gl", "PDB-CPX-172174", tmp_path_factory.mktemp("sting"))


@pytest.fixture(scope="module")
def spike_ace2(tmp_path_factory) -> WorkflowRun:
    """The Spike RBD-ACE2 heterodimer run (`6m0j` -> PDB-CPX-140195), scoped to the
    five interface instances selected to exercise partner-order normalisation."""
    return run_workflow("6m0j", "spike_ace2_reversal", tmp_path_factory.mktemp("spike"))
