# Tranche 2A final report: frozen STING end-to-end regression

Completed 2026-09-23. Summary of the first scientific end-to-end regression test,
the STING homodimer working example `11gl` / `PDB-CPX-172174`, run offline on
frozen real PDBe responses through the production functions the notebook calls. The
detailed tranche document is `testing_tranche2a_sting_e2e.md`; the fixture's own
provenance is in `tests/fixtures/recorded/PDB-CPX-172174/README.md`.

Not implemented, by instruction: Spike/ACE2 regression, triosephosphate isomerase
regression, live PDBe tests, notebook smoke tests, CI, fixture-refresh automation.

---

## Fixture

```text
location   tests/fixtures/recorded/PDB-CPX-172174/
files      complex_details_11gl.json, complex_details.json, interface_interactions.json,
           mutations_post.json, modifications_post.json (recorded 404), bound_molecules.json
           (12 entries, 4kc0 is a 404), ligand_interactions.json (15 instances),
           metadata.json, README.md
size       664 KB, of which ligand_interactions.json is 503 KB
source     live PDBe API v2, captured once on 2026-09-23 by wrapping requests.Session.request
           around the notebook's Phase 1 call sequence with the default Config("11gl")
requests   32 recorded; the fixture is frozen and not refreshed automatically
```

Each file keeps the raw response body with its HTTP status, so the production API
layer parses it exactly as it would live data. The ligand file covers only the
instances the default blocklist lets through, because the workflow never requests
the others. The README records what was captured, how, why it is frozen and what it
protects.

---

## Workflow coverage

`tests/e2e/conftest.py` holds a `RecordedPDBe` router that replaces
`requests.Session.request` for Phase 1 only and raises on any unrecorded URL, plus
`run_workflow`, which is the notebook's Phase 1 to Phase 7 call sequence unchanged.
A module-scoped `sting` fixture runs it once.

Exercised, in order:

- identifier resolution through both `complex/details` lookups and the dimer check;
- partner map and assembly metadata;
- interface retrieval and assembly selection;
- batched mutation POST; modification POST on the 404 path;
- bound-molecule retrieval including one 404, blocklist filtering, per-ligand fetches
  on the thread pool;
- interface-record construction, partner consistency, comparable-record selection;
- Jaccard similarity matrix and clustering;
- annotation overlap;
- structure table, cluster interpretation report, frequency summary, conservation,
  rewiring;
- JSON export to a temporary directory.

Not exercised: plots, the explorer widget, Mol* rendering.

---

## Scientific assertions

Expected values were taken from `spec/new/demo_notes.md` and verified against the
frozen fixture before being encoded. Cluster ids are never asserted numerically.

- **Interface cohort**: 14 instances from the 12 documented entries; no assembly
  exclusions, no dropped contacts, no microheterogeneity, no zero-comparable-contact
  exclusions; both roles `Q3TBT3`; bond vocabulary 169 hydrogen bonds and 35 salt
  bridges; matrix 14 by 14, symmetric, unit diagonal; all 32 requests served from
  the fixture.
- **Cluster relationships**: size multiset `[1, 1, 3, 9]` at cut 0.6; `11gl` with
  `11gm` (both assemblies), `4loj`, `6xnn`; `4kby` with `4kc0`, `9ltf`; `11gl` apart
  from `4kby`; `4jc5` and `4lol` singletons; the four representatives in four
  distinct states.
- **Rewiring**: A232–D209, A232–Y260 and D209–G233 hydrogen bonds in 9/9 of the
  large state and 0/3 of the small; the D273–H156 salt bridge the reverse. Asserted
  on the production `((Q3TBT3, pos, 1), (Q3TBT3, pos, 2), bond)` tuples and on the
  rewiring table's direction labels; no `shared core` contact between these states.
- **Annotations**: `1YE` on `4lol` Ser161 of both copies, mapped to `Q3TBT3` 161
  with roles 1 and 2, the only ligand at that interface; `C2E` on `4kby` including
  Arg237 of both copies; `1YD` on `4lok` in the large state; `A1ELY` on `9ltf` in
  the small state; nothing on `4jc5`; no interface mutations or modifications in
  the data.
- **Frequency and conservation**: denominator 14; D209 (role 1) most frequent at
  10/14; nothing conserved at 0.8; at 10/14 the conserved set includes D209 and
  K235; the A232–D209 pair in exactly 9 interfaces.
- **QC**: report ordered largest first; `4jc5` singleton with median 4 residue pairs
  and both the singleton and sparse flags; `4lol` singleton flag only with `1YE` in
  its ligand column; the large state unflagged. Flags asserted by keyword, not full
  sentence.
- **Export**: file written under `tmp_path`; correct complex id, oligomeric state and
  14 interfaces; both roles of `Q3TBT3` with gene name `Sting1`; D209 role 1 at 10
  interfaces, frequency 0.7143, level `medium`, label `D209`; A232–D209 hydrogen
  bond in 9 and D273–H156 salt bridge in 3; every contact role 1 then role 2; residue
  frequencies on disk equal the returned dict.

Intentionally not asserted: numeric cluster ids, complete membership, rewiring or
ligand tables, area and resolution statistics, full JSON snapshots, warning wording,
the `ZNT` contact at the `11gl` interface (present in the fixture, absent from the
demo notes), and anything about the live database.

---

## Results

```text
tests before      87
new E2E tests      5
tests after       92

passed            92
failed             0
skipped            0
xfail              0

runtime           1.9 s for the whole suite (2.7 s wall clock); e2e about 0.1 s
```

All five e2e tests passed on their first run, which is expected for
characterisation of behaviour that was verified against the fixture before being
encoded. No test was made to fail artificially.

Tests:

```text
tests/e2e/test_sting_recorded.py::test_sting_recorded_workflow_builds_expected_interface_states
tests/e2e/test_sting_recorded.py::test_sting_recorded_reproduces_known_rewiring
tests/e2e/test_sting_recorded.py::test_sting_recorded_maps_known_interface_ligands
tests/e2e/test_sting_recorded.py::test_sting_recorded_reports_expected_qc_and_frequencies
tests/e2e/test_sting_recorded.py::test_sting_recorded_export_is_consistent
```

---

## Network

Fully offline. Every request in the e2e run is answered by the recorded router,
which raises on anything unrecorded, and the root `block_network` guard stays active
outside the Phase 1 context. The first test asserts that exactly 32 calls went
through the router. A plain `pytest` makes no external requests.

---

## Production changes

None. `pdbe_interfaces/` and `notebook.ipynb` are unchanged in this tranche. No
testability plumbing was needed: the notebook's call sequence could be reproduced in
test code as-is.

Files added or changed in this tranche:

```text
tests/e2e/conftest.py                              recorded router, run_workflow, sting fixture
tests/e2e/test_sting_recorded.py                   the five regression tests
tests/fixtures/recorded/PDB-CPX-172174/            frozen fixture (nine files)
spec/new/testing_tranche2a_sting_e2e.md            detailed tranche document
spec/new/testing_tranche2a_final_report.md         this document
spec/new/testing_strategy.md                       one status note under Candidate A
```

---

## Unexpected findings

No discrepancy with the trusted STING analysis. The frozen data reproduces every
value in the demo notes of 2026-08-24: cohort, state sizes, membership, the four
state-separating contacts, ligands per state and the QC flags on `4jc5`.

Two data observations, neither a workflow issue:

- a `ZNT` ligand contact at the `11gl` interface that the demo notes did not list
  (not asserted);
- raw mutation records exist for 8 of the 12 entries but none maps to an interface
  residue, so no mutation assertion was possible.

Nothing is committed. Spike/ACE2, TIM, live tests, notebook smoke tests and CI
remain unstarted pending review of STING.
