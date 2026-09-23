# Tranche 2A: frozen STING end-to-end regression

Completed 2026-09-23. First scientific end-to-end regression test: the STING
homodimer working example, `11gl` / `PDB-CPX-172174`, run offline on frozen real
PDBe responses through the same production functions the notebook calls.

The question it answers: given frozen real input for a biological example whose
behaviour is already understood, does the complete scientific workflow still
reproduce the expected interface states, rewiring, annotations and export? It is a
regression test of our workflow, not of the live PDBe database.

Not included, by instruction: Spike/ACE2, triosephosphate isomerase, live PDBe
integration tests, notebook smoke tests, CI, fixture-refresh automation.

---

## Why STING

It is the notebook's working example and the one whose behaviour the repository
documents in most detail (`spec/new/demo_notes.md`, 2026-08-24). It is small enough
to run in well under a second on recorded data, all X-ray, with no partner reversals
and no assembly exclusions, so it isolates the core pipeline from the filtering
logic. It still exercises every stage: four interface states, rewiring between the
two large ones, ligands in three states, and QC flags on a sparse singleton. The
tranche-1.5 orientation diagnostic showed a maximum orientation delta of about 0.15
and no clustering change for this complex, so the states are stable production
behaviour.

---

## Fixture provenance

Location: `tests/fixtures/recorded/PDB-CPX-172174/` (664 KB, nine files, README and
`metadata.json` included).

| File | Size | Endpoint |
|---|---|---|
| `complex_details_11gl.json` | 0.5 KB | `complex/details/11gl?id_type=pdb_id` |
| `complex_details.json` | 5.8 KB | `complex/details/PDB-CPX-172174?id_type=pdb_complex_id` |
| `interface_interactions.json` | 128 KB | `complex/interface_interactions/PDB-CPX-172174` |
| `mutations_post.json` | 8.3 KB | POST `pdb/entry/mutated_AA_or_NA`, one chunk of 12 ids |
| `modifications_post.json` | 0.2 KB | POST `pdb/entry/modified_AA_or_NA`, recorded 404 |
| `bound_molecules.json` | 6.6 KB | `pdb/bound_molecules/{pdb_id}` for 12 entries (`4kc0` 404) |
| `ligand_interactions.json` | 503 KB | `pdb/bound_ligand_interactions/...` for the 15 ligand instances surviving the default blocklist |
| `metadata.json` | 0.8 KB | complex id, starting id, capture date, retained ids, config, request count |

Captured once on 2026-09-23 by running the notebook's Phase 1 call sequence with the
default `Config(identifier="11gl")` while `requests.Session.request` was wrapped to
record every request and response; 32 responses were recorded and grouped by
endpoint. Each file keeps the raw body and HTTP status. The fixture is frozen and is
not refreshed automatically. The ligand file covers only the instances the default
blocklist lets through, because the workflow never requests the others.

---

## How the test runs the workflow

`tests/e2e/conftest.py` provides `RecordedPDBe`, a stand-in for
`requests.Session.request` that routes each URL to the recorded response and raises
on anything unrecorded. `run_workflow` executes the notebook's Phase 1 to Phase 7
call sequence unchanged, with the recorded session active only for Phase 1 (inside a
`MonkeyPatch.context()`), and returns a `WorkflowRun` holding every intermediate. A
module-scoped `sting` fixture runs it once per test module. The root `block_network`
guard remains active for everything else, so the test cannot reach PDBe.

Production stages exercised, in order: `resolve_complex_id` (entry-id lookup,
complex details, dimer check), `build_partner_map`, `extract_assembly_metadata`,
`fetch_interface_interactions`, `select_interfaces`, `fetch_mutations`,
`fetch_modifications` (404 path), `collect_ligand_contacts_for_entries` (bound
molecules incl. one 404, blocklist, per-ligand fetch on the thread pool),
`build_interface_records`, `check_partner_consistency`, `select_comparable_records`,
`jaccard_similarity_matrix`, `cluster_interfaces`, `overlap_annotations`,
`build_structure_table`, `cluster_interpretation_report`, `interface_frequency_summary`,
`conserved_residues`, `compare_cluster_contacts`, `export_interface_frequency_json`.
Not exercised: plots, the explorer widget, Mol* rendering.

---

## Scientific expectations protected

Expected values were taken from the demo notes and verified against the frozen
fixture before being encoded. Five tests:

**`test_sting_recorded_workflow_builds_expected_interface_states`**

- All 32 requests answered from the fixture.
- `11gl` resolves to `PDB-CPX-172174`; homodimer with two chains; partner map has
  `Q3TBT3` in both roles.
- 14 interface instances from 12 entries (exact, fixture-specific), no assembly
  exclusions, no dropped contacts, no microheterogeneity, no zero-comparable-contact
  exclusions; bond vocabulary `hydrogen_bond: 169, salt_bridge: 35`.
- Similarity matrix 14 by 14, symmetric, unit diagonal.
- Cluster size multiset `[1, 1, 3, 9]` at the default cut 0.6.
- Relationships: `11gl`, `11gm` (both assemblies), `4loj`, `6xnn` co-cluster;
  `4kby`, `4kc0`, `9ltf` co-cluster; `11gl` and `4kby` differ; `4jc5` and `4lol` are
  singletons; the four representatives fall in four distinct states.
- Structure table has 14 rows, role-1 label `Sting1 (role 1)`, all X-ray.

**`test_sting_recorded_reproduces_known_rewiring`**

- In the production `((Q3TBT3, pos, 1), (Q3TBT3, pos, 2), bond)` representation:
  A232–D209, A232–Y260 and D209–G233 hydrogen bonds are in 9/9 members of the
  large state and 0/3 of the small one; the D273–H156 salt bridge is in 3/3 of the
  small state and 0/9 of the large one.
- `compare_cluster_contacts` labels the first three `higher in A` and the salt
  bridge `higher in B`, and reports no `shared core` contact between these states.

**`test_sting_recorded_maps_known_interface_ligands`**

- No engineered mutation at any interface; modifications empty (404 path).
- `1YE` at the `4lol` interface on Ser161 of both copies, mapped to `Q3TBT3` 161
  with roles 1 and 2; `1YE` is the only ligand at that interface.
- `C2E` (cyclic di-GMP) at the `4kby` interface, including Arg237 of both copies.
- `1YD` on `4lok` in the large state; `A1ELY` on `9ltf` in the small state; no
  ligand on `4jc5`.

**`test_sting_recorded_reports_expected_qc_and_frequencies`**

- Report ordered largest state first (`[9, 3, 1, 1]`).
- `4jc5`: singleton, median residue-pair count 4, both the singleton and the sparse
  warning present. `4lol`: singleton warning only, `1YE` in its ligand column. The
  large state carries no warning.
- Frequency denominator is 14; the most frequent residue is D209 (role 1) in
  10/14; nothing reaches the 0.8 conservation threshold; at 10/14 the conserved
  set includes D209 and K235; the A232–D209 pair is in exactly 9 interfaces.

**`test_sting_recorded_export_is_consistent`**

- File written under `tmp_path`; complex id, oligomeric state and 14 interfaces in
  the metadata; residue frequencies on disk equal the returned dict.
- Both roles of `Q3TBT3` with gene name `Sting1` in `partners`.
- D209 role 1: 10 interfaces, frequency 0.7143, level `medium`, label `D209`.
- A232–D209 hydrogen bond in 9 interfaces, D273–H156 salt bridge in 3; every
  contact has `partner_1.role == 1` and `partner_2.role == 2`.

Warning categories are asserted by keyword (`singleton`, `sparse`), not by full
sentence. Cluster ids are never asserted numerically.

---

## Intentionally not asserted

- Numeric cluster ids.
- The complete cluster membership table (a few representative relationships stand
  for it), the complete rewiring table, the complete ligand list.
- Exact area, resolution or contact-density statistics.
- Full JSON snapshot or byte-for-byte file comparison.
- Warning wording character for character.
- `ZNT` at the `11gl` interface: present in the fixture but not listed in the demo
  notes, so left out of the assertions and recorded here as an observation only.
- Anything about the live database: the fixture is frozen; counts will differ live.

---

## Results

```text
tests before      87
new e2e tests      5
tests after       92
passed            92
failed             0
skipped            0
xfail              0
runtime           1.9 s for the suite; the five e2e tests take about 0.1 s
```

The STING tests passed on the first run, as expected for characterisation of
behaviour that was verified against the fixture beforehand. No test was made to
fail artificially.

Network: the e2e run is fully offline. Every URL is served by `RecordedPDBe`, which
raises on any unrecorded request, and the root `block_network` fixture is active for
the rest of each test. The request-count assertion confirms all 32 calls went
through the router.

---

## Production changes

None. `pdbe_interfaces/` and `notebook.ipynb` were not touched in this tranche. No
testability plumbing was needed: the notebook's call sequence could be reproduced in
test code as-is.

## Unexpected findings

None that affect the workflow. The frozen fixture reproduces every value in the demo
notes of 2026-08-24: 14 instances, 12 entries, four states of 9/3/1/1 with the same
membership, the same four state-separating contacts, the same ligands per state,
and the same QC flags on `4jc5`. The only difference from the demo notes is the
additional `ZNT` ligand contact at the `11gl` interface, which the notes did not
list; it is a data observation, not a workflow discrepancy, and is not asserted.
Raw mutation records exist for 8 of the 12 entries but none maps to an interface
residue, so the fixture offers no meaningful mutation assertion.
