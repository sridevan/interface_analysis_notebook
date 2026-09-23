# Frozen regression fixture: STING homodimer, `11gl` / `PDB-CPX-172174`

Real PDBe API v2 responses for the notebook's working example, captured once on
2026-09-23 and frozen. Used by `tests/e2e/test_sting_recorded.py` to run the complete
scientific workflow offline and check that it still reproduces the STING interface
states, rewiring, annotations and export described in `spec/new/demo_notes.md`.

This fixture is **not refreshed automatically** and must stay frozen until it is
deliberately reviewed and re-captured. It is a regression fixture for our workflow,
not a test of the live PDBe database, which will have grown since the capture date.

## What was captured

| File | Endpoint | Notes |
|---|---|---|
| `complex_details_11gl.json` | `complex/details/11gl?id_type=pdb_id` | entry-id lookup that resolves to the complex |
| `complex_details.json` | `complex/details/PDB-CPX-172174?id_type=pdb_complex_id` | participants, assemblies, method, resolution |
| `interface_interactions.json` | `complex/interface_interactions/PDB-CPX-172174` | 14 PISA interfaces with residue-pair contacts and UniProt mappings |
| `mutations_post.json` | POST `pdb/entry/mutated_AA_or_NA` | one chunk for the 12 retained entries; `request_body` is the exact id string sent |
| `modifications_post.json` | POST `pdb/entry/modified_AA_or_NA` | returned 404 (no modifications on file), stored as such |
| `bound_molecules.json` | `pdb/bound_molecules/{pdb_id}` | one entry per retained PDB id; `4kc0` returned 404 |
| `ligand_interactions.json` | `pdb/bound_ligand_interactions/{pdb}/{chain}/{res}?preserve_case=false` | one entry per ligand instance surviving the default blocklist (15 instances) |
| `metadata.json` | | complex id, starting identifier, capture date, retained ids, config used, request count |

Each file stores the raw response body with its HTTP status, so the production
`pdbe_interfaces.api` layer parses it exactly as it would a live response.

## How it was generated

The notebook's Phase 1 call sequence was run once against the live API with the
default `Config(identifier="11gl")` while `requests.Session.request` was wrapped to
record every request and response. The 32 recorded responses were then grouped by
endpoint into the files above. The ligand-interaction file contains only the ligand
instances that survive the default `ligand_blocklist`, because the workflow never
requests the others; a run with a different blocklist would need a fresh capture.

## What it protects

Identifier resolution, dimer check, assembly filtering, batched annotation retrieval,
ligand retrieval and filtering, interface-record construction, partner consistency,
comparable-record selection, Jaccard similarity, clustering, annotation overlap,
conservation and frequency, the cluster interpretation report, rewiring, and the
JSON export, all on real data whose expected behaviour is documented.
