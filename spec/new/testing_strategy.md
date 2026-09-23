# Testing strategy: Aggregated Interface Interaction Analysis in Dimers

Proposed automated testing strategy for the notebook and the `pdbe_interfaces/`
modules. Written from a read-through of the code, the notebook, the spec and the
demo notes on 2026-09-23. No tests have been written yet; this document is the
plan to review before implementation.

The tests are organised around what must work correctly for the scientific result
to be trusted, not around coverage. Four layers:

1. Core functions required by the scientific workflow.
2. Input / identifier validation.
3. API tests: mocked API behaviour plus a small number of live integration tests.
4. End-to-end scientific regression tests on known biological examples.

---

## 1. Workflow map

The notebook is a thin driver. Every step below is a call into `pdbe_interfaces/`,
and the arrows are the real data dependencies.

```text
Config.identifier  ("11gl" or "PDB-CPX-172174")
    ↓  api.resolve_complex_id(require_dimer=True)
       complex/details?id_type=pdb_id  →  pick complex (preferred assembly / only dimer)
       complex/details?id_type=pdb_complex_id  →  _require_dimer(total_chains == 2)
    ↓  outputs.build_partner_map / extract_assembly_metadata      (gene labels, method, resolution)
    ↓  api.fetch_interface_interactions(complex_id)                (one list, all PISA interfaces)
    ↓  representation.select_interfaces
       assembly_exclusions: bound_macromolecules (always), max_resolution (optional)
       drop interfaces not listed in complex/details; apply max_entries
    ↓  annotation retrieval for the retained pdb_ids
       api.fetch_mutations / fetch_modifications        (batched POST, 50 ids per chunk)
       annotations.collect_ligand_contacts_for_entries
         bound_molecules per entry → filter_bound_molecules (blocklist, glycans)
         bound_ligand_interactions per ligand instance → _parse_ligand_interactions
    ↓  representation.build_interface_records  (_build_one)
       author key  (pdb, chain, resnum, ins)   → author_pairs   (annotation join space)
       UniProt key (acc, unp_pos, role 1|2)    → uniprot_pairs  (comparison space)
       drop pairs lacking UniProt on either side; microheterogeneity keeps first mapping
    ↓  representation.check_partner_consistency  (_reverse_roles)
       majority (acc_1, acc_2) ordering is canonical; minority records have roles swapped
       (no effect on homodimers: both accessions are equal, ordered tuples are kept as-is)
    ↓  similarity.jaccard_similarity_matrix (typed and untyped sets)
    ↓  similarity.cluster_interfaces  (1 - Jaccard, average linkage, fcluster at distance cut)
       similarity.sweep_cluster_cuts   (advisory table)
    ↓  annotations.overlap_annotations
       author key → author_to_uniprot → (acc, pos, role); off-interface annotations dropped
    ↓  outputs.build_structure_table, describe_conservation
       conserved_residues / conserved_interaction_pairs  (fraction of interfaces >= threshold)
    ↓  outputs.cluster_interpretation_report
       _core_contacts, _contact_profile, _label_profile, _qc_warnings
    ↓  outputs.interface_frequency_summary (explorer), rewiring_table → compare_cluster_contacts
    ↓  visualize.visualize_clusters_grid  (Mol*, downloads CIF; not testable offline)
    ↓  outputs.export_interface_frequency_json  →  {output_dir}/{complex_id}.json
```

Two facts about the map matter for the tests.

- Everything after Phase 1 is pure Python over `InterfaceRecord` objects and dicts,
  so the entire scientific core runs offline once the six API responses are provided.
- The only place orientation is normalised is `check_partner_consistency`, and it acts
  on accession order only. There is no symmetrisation for homodimers. The spec states
  this is intended (symmetric contacts are kept as distinct ordered tuples), which
  means homodimer correctness rests entirely on PISA labelling chains consistently
  across entries. The tests should make that assumption visible rather than hide it.

---

## 2. Critical functions

| Function | What it does | Why a wrong answer corrupts the science | Test type |
|---|---|---|---|
| `api.resolve_complex_id` | Turns an entry id or complex id into `(complex_id, details)`, choosing among multiple complexes and enforcing the dimer rule | Picks the wrong complex, or accepts a non-dimer, and every downstream number describes the wrong thing | Unit with mocked `fetch_complexes_for_pdb_id` / `fetch_complex_details`: single complex, preferred-assembly dimer, only-dimer fallback, ambiguous, non-dimer, empty, casing and whitespace |
| `representation.assembly_exclusions` / `select_interfaces` | Drops assemblies with bound macromolecules and, optionally, poor or missing resolution; applies `max_entries`; records reasons | Including a Fab-bound instance mixes non-comparable interfaces into the clusters; excluding wrongly silently shrinks the dataset | Unit on hand-built `complex_details`: bound list present, resolution None with and without the filter, unparseable resolution, interface whose assembly is absent from details, nothing surviving raises |
| `representation._build_one` | Builds the author and UniProt pair sets, role labels, residue identities, drop and collision counts | This is the residue correspondence. Wrong role, wrong position, or a dropped pair changes the fingerprint that everything else compares | Unit on synthetic interaction lists: role assignment, insertion-code normalisation, missing UniProt on one side drops the pair but keeps the author pair, duplicate atom-level records collapse to one tuple, microheterogeneity keeps first mapping and does not add the colliding author key to `author_to_uniprot` |
| `representation.check_partner_consistency` + `_reverse_roles` | Makes role 1 the same accession in every record | Without it, a reversed heterodimer entry shares no tuples with the others and forms a spurious cluster. A wrong swap leaves `author_to_uniprot` or `residue_identity` on the old role, so annotations land on the wrong partner | Unit: reversed record ends with canonical accession in role 1, every tuple, `author_to_uniprot` and `residue_identity` flipped consistently, Jaccard against an unreversed twin is 1.0; homodimer records untouched; non-matching accession pair left alone; a 50/50 tie is deterministic |
| `similarity.jaccard` / `jaccard_similarity_matrix` | Pairwise set similarity | Everything about states comes from this matrix | Unit: identical, disjoint, partial, both empty, symmetry, diagonal 1.0, typed differs from untyped when only bond type differs |
| `similarity.cluster_interfaces` | Average-linkage clustering, cut at `distance_cut` | Wrong linkage or wrong criterion changes the number and membership of states | Unit on a designed 3-block matrix: known membership at a cut, one interface returns a single cluster, cut monotonicity (fewer clusters at higher cut), `flat_assignment` index order matches record order |
| `annotations.filter_bound_molecules` | Flattens bound-molecule records and applies blocklist and glycan filter | A missed blocklist entry floods the ligand overlay with sulfate; dropping a real ligand hides the correlate | Unit: blocklist hit, carbohydrate-polymer only when flag set, missing chain or residue number skipped, insertion code normalised |
| `annotations._parse_ligand_interactions` | Expands atom-level ligand contacts into `(residue, ligand instance, contact type)` | Contact residues off by one chain or number map the ligand to the wrong interface residue | Unit: multiple `interaction_details` labels expand, atom duplicates collapse, response ligand metadata overrides requested key, missing end residue skipped |
| `annotations.overlap_annotations` | Joins mutations, modifications, ligands onto interface residues by author key and attaches the UniProt key | This is "annotation mapped to the wrong residue". A join that ignores chain, insertion code, or the mutation type filter puts the mark on the wrong residue or the wrong partner | Unit: annotation on interface residue lands with the right `(acc, pos, role)`; same residue number on a non-interface chain is dropped; residue in the interface but without UniProt gets the `(None, None, None)` sentinel; `mutation_type_filter` excludes Conflict by default; one entry with two assemblies produces two rows |
| `outputs.conserved_residues` / `conserved_interaction_pairs` / `_core_contacts` | Frequency across interfaces against a threshold | A residue counted twice per interface, or a `>` instead of `>=` at the boundary, changes what is called conserved | Unit: residue appearing in several tuples of one interface counts once; exact-boundary case (e.g. 4 of 5 at 0.8) is included; threshold 1.0 requires all |
| `outputs.interface_frequency_summary` / `export_interface_frequency_json` | Per-residue and per-pair counts and fractions; `pair_matrix`; JSON with conservation levels | These are the numbers a downstream Mol* colouring consumes | Unit: fractions use interface count as denominator, `pair_matrix[i, j]` matches the `pairs` table, `partner_1` is always role 1, conservation level boundaries (0.80, 0.50, 0.20 inclusive), file written to `config.output_dir` under `tmp_path`, typed vs untyped row counts |
| `outputs.cluster_interpretation_report` | Per-cluster stats, core contacts, in/out counts, annotation profiles, QC warnings, size-descending order | Wrong in/out denominators or entry de-duplication misstates the evidence for a state | Unit on a 3-cluster synthetic set: `in_cluster_size` and `out_cluster_size` sum to N, entry counts de-duplicate assemblies of one entry, singleton and sparse warnings fire where designed, the dominant cluster never gets the range warning, rows ordered by size |
| `outputs.compare_cluster_contacts` / `rewiring_table` | Two-state contact comparison with `shared core` / `higher in A` / `higher in B` / `rare` labels | The rewiring narrative is read directly from these labels | Unit: known fractions produce known labels at the 0.8 and 0.5 boundaries, default picks the two largest non-singleton clusters, fewer than two non-singletons returns the explanatory message |
| `api._retry_request` / `_get_json` / `_batched_post` | Retry, 404 policy, chunked POST merge | A swallowed failure on a critical endpoint yields an incomplete dataset that still "runs" | Mocked API tests (section 4) |

Deliberately left out: `Config.describe`, `label_for_record`, `describe_complex`,
`_stats`, formatting helpers, plots, explorer and Mol* builders. They cannot change a
scientific result, only its presentation.

---

## 3. Identifier handling

### Accepted forms

`Config.identifier` is a string. If it matches `^PDB-CPX-\d+$` case-insensitively, it
is used directly after upper-casing. Anything else is stripped, lower-cased and treated
as a PDB entry id. There is no other input identifier; everything else is configuration.

### Where validation happens today

All of it is inside `api.resolve_complex_id` and the two fetchers it calls. `Config`
itself validates nothing.

### Syntactic checks testable offline

Exercised with the fetchers monkeypatched, no network.

- Empty or whitespace-only string raises `ValueError` before any request.
- `"pdb-cpx-172174"` and `" PDB-CPX-172174 "` normalise to `PDB-CPX-172174` and go to
  the complex-id endpoint.
- `"11GL"` and `" 11gl "` normalise to `11gl` and go to the entry-id endpoint.
- Malformed complex ids such as `PDB-CPX-`, `PDB-CPX-12A`, `CPX-2128` (the Complex
  Portal id, which users will try) fall through to the entry path. That is current
  behaviour and a test should pin it, because it means a typo in a complex id produces
  a confusing "did not resolve to any complex" message.
- There is no offline check that an entry id looks like one (four characters, leading
  digit, or the extended `pdb_0000xxxx` form). A test asserting that `"not-an-id"`
  raises without a network call would fail today. Whether to add that check is a
  decision for the maintainer; write the test and mark it `xfail` with a reason so the
  gap is recorded.

### Semantic checks that need PDBe

Live, marked integration, one call each.

- Valid entry that resolves to a dimer: `11gl` → `PDB-CPX-172174` (documented in the
  README and demo notes).
- Valid entry that is not a dimer: `4hhb` raises `ValueError` naming both the entry and
  the complex with its chain count (demo notes).
- Valid entry mapping to several complexes: the code handles it but no example is
  recorded in the repository, so do not assert one until a case has been found.
- Nonexistent entry id (e.g. `0xxx`): the entry endpoint returns 404, which
  `allow_404` turns into `ValueError`.
- Nonexistent complex id (e.g. `PDB-CPX-999999999`): `fetch_complex_details` does not
  pass `allow_404`, so a 404 surfaces as `requests.HTTPError`, not the `ValueError` the
  spec describes. A mocked unit test should pin whichever behaviour is wanted; the live
  test should confirm what PDBe actually returns for an unknown complex id, which
  cannot be known offline.
- A two-chain complex where one participant is not UniProt-mapped: `_require_dimer`
  only checks `total_chains == 2`, and `build_partner_map` then yields a single
  partner. The README says both components must map to UniProt, so this is another
  documented-but-unenforced rule to capture as a test with a synthetic
  `complex_details`.

---

## 4. API dependencies

All six endpoints are PDBe API v2 over HTTPS with no authentication. A seventh
external resource, the mmCIF download used by Mol*, is only touched by `visualize`
and should stay out of automated tests.

### Recommended mocking layer

Patch `requests.Session.request` at the class level with `monkeypatch`. That covers
the thread-local sessions used by the parallel fetchers, so no code change is needed.
Patch `time.sleep` and set `api.BACKOFF_BASE_S` to zero in the same fixture so retry
tests run in milliseconds. For tests of the representation and annotation layers,
patch the `api.fetch_*` functions instead and feed recorded JSON.

### Per-endpoint expectations

| Endpoint | What the workflow expects | Mocked cases worth testing | Live test |
|---|---|---|---|
| `complex/details/{pdb_id}?id_type=pdb_id` | `{pdb_id: [ {pdb_complex_id, assemblies:[{preferred_assembly}], ...} ]}` | 200 with one complex; 200 with two complexes and a preferred flag; 200 with empty list → `ValueError`; 404 → `ValueError`; 500 then 200 → succeeds after retry; five 503s → raises; `Timeout` five times → raises | Yes, one call: `11gl` resolves to `PDB-CPX-172174` and carries `pdb_complex_id` |
| `complex/details/{cpx}?id_type=pdb_complex_id` | `{cpx: [ {name, total_chains, oligomeric_state, participants:[{accession, accession_type, gene_name, stoichiometry}], assemblies:[{pdb_id, assembly_id, resolution, experimental_method, bound_macromolecules}] } ]}` | Full record; `total_chains` 4 → dimer error; missing `assemblies` key → no exclusions and no metadata rather than a crash; missing `participants` → empty partner map; 404 → currently `HTTPError` | Yes: field presence check for `PDB-CPX-172174` (`total_chains == 2`, participants contain `Q3TBT3`, every assembly has `pdb_id` and `assembly_id`). This is the test that catches an upstream schema change |
| `complex/interface_interactions/{cpx}` | `{cpx: [ {entry_id, assembly_id, interface_id, interface_info{interface_area,...}, interactions:[{auth_asym_id_1, auth_seq_id_1, unp_accession_1, unp_seq_id_1, unp_one_letter_code_1, ..._2, bond_type}] } ]}` | Normal; empty list → `ValueError`; key missing → `ValueError`; interaction missing `unp_seq_id_2` → dropped and counted; `auth_ins_code_1` vs `ins_code_1` fallback; whitespace insertion code; interface with no interactions produces an empty record | Yes: for `PDB-CPX-172174`, at least one interface for entry `11gl`, every interaction carries both UniProt accessions (the spec records that none lacked a mapping), bond types non-empty |
| `pdb/entry/mutated_AA_or_NA` (POST) | `{pdb_id: [ {chain_id, author_residue_number, author_insertion_code, mutation_details{from, to, type}} ]}` | Body is the JSON-encoded comma-joined id string; 120 ids produce three POSTs and a merged dict; 404 on one chunk → that chunk empty, others kept; 500 retried; record with `type: Conflict` filtered out later | Yes, cheap: POST for `7bh9` returns a record with `chain_id E`, residue 477, type `Engineered mutation` (from the saved notebook output). Guards the author-chain semantics the join depends on |
| `pdb/entry/modified_AA_or_NA` (POST) | Same shape with `chem_comp_id`, `chem_comp_name` | Same as above; 404 for all chunks yields `{}` and the workflow continues with zero modifications (the normal case for both examples) | Optional. The saved run shows 404 for every chunk on Spike/ACE2, so a live test would only confirm the empty path |
| `pdb/bound_molecules/{pdb_id}` | `{pdb_id: [ {composition:{ligands:[{chem_comp_id, chain_id, author_residue_number, author_insertion_code, molecule_type}]}} ]}` | 404 → `[]`; persistent 500 → `[]` with warning (tolerated); ligand missing `chain_id` skipped; `_parallel_map` preserves input order when the first item is slowest | Yes, one call for `4lol`: response parses and contains a non-blocklisted ligand (`1YE` per demo notes) |
| `pdb/bound_ligand_interactions/{pdb}/{chain}/{res}` | `{pdb_id: [ {ligand:{chem_comp_id, chain_id, author_residue_number}, interactions:[{end:{chain_id, author_residue_number, author_insertion_code}, interaction_details:[...]}]} ]}` | 404 → `[]`; malformed `end` skipped; several `interaction_details` expand; two atom records collapse | One call, for the `4lol` `1YE` instance, asserting at least one contact residue comes back. Slow endpoint, so keep it to one |

### Cross-cutting mocked tests

- `RETRY_STATUSES` (429, 500, 502, 503, 504) are retried; other 4xx are not.
- `MAX_ATTEMPTS` is respected and the last exception is re-raised.
- A `tolerate_failure` call never raises; a critical call does.
- `_batched_post` chunks at `POST_BATCH_SIZE` and merges; an empty id list makes no
  request.
- `_parallel_map` returns results in input order and falls back to serial for one item.

Live tests should be few, marked, and skipped by default. One per endpoint is enough,
and they should assert shape and one known value, never counts.

---

## 5. Scientific end-to-end regression test

### Candidate A: `11gl` / `PDB-CPX-172174`, STING homodimer

*Status (2026-09-23): implemented as the frozen offline regression
`tests/e2e/test_sting_recorded.py` on the fixture under
`tests/fixtures/recorded/PDB-CPX-172174/`; see `testing_tranche2a_sting_e2e.md`.
The live variant, Candidates B and C, and the notebook smoke test remain open.*

The notebook's working example. The demo notes record its behaviour in detail as of
2026-08-24. It is small (14 interfaces, 12 entries, all X-ray), runs in a few seconds,
exercises clustering, rewiring, ligand overlap and QC warnings, and has no partner
reversals or exclusions, so it isolates the core pipeline from the filtering logic.

Facts already demonstrated that can serve as expected results:

- Four states at the default cut, sizes 9, 3, 1, 1, with a plateau from cut 0.55 to 0.70.
- The nine-member state contains `11gl`, `11gm`, `11gn`, `4loj`, `4lok`, `4yp1`,
  `6xnn`; the three-member state contains `4kby`, `4kc0`, `9ltf`; `4jc5` and `4lol`
  are singletons.
- Rewiring between the two large states: `A232-D209`, `A232-Y260` and `D209-G233`
  hydrogen bonds at 100% in the large state and 0% in the small one; the `D273-H156`
  salt bridge at 0% versus 100%.
- Ligands at the interface: `1YD`, `2BA`, `V67` in the large state; `A1ELY`, `C2E` in
  the three-member state; `1YE` in `4lol`; none in `4jc5`.
- `4jc5` carries both the singleton and the sparse-fingerprint warning (4 residue pairs).
- Bond types are `hydrogen_bond` and `salt_bridge` only.

Stable assertions (robust to new depositions):

- The listed entries co-cluster; `4kby` and `11gl` are in different clusters.
- The four rewiring contacts have the stated direction between the clusters containing
  `11gl` and `4kby`.
- `1YE` maps to an interface residue of `4lol` and `C2E` to one of `4kby`.
- Every record has `unp_accession_1 == unp_accession_2 == Q3TBT3`.
- Bond types are a subset of `{hydrogen_bond, salt_bridge}`.
- The similarity matrix is symmetric with unit diagonal.
- The JSON export has both roles for `Q3TBT3` and `partner_1.role == 1` for every
  contact.

Brittle assertions (will break as PDBe grows): 14 interfaces, 12 entries, exactly four
clusters, cluster ids, the exact plateau boundaries, "zero exclusions", exact fractions
and medians, `hydrogen_bond: 169`, and any comparison against a complete JSON dump.

### Candidate B: `6m0j` / `PDB-CPX-140195`, Spike RBD with ACE2 heterodimer

*Status (2026-09-23): a scoped offline regression covering partner-order
normalisation and annotation identity is implemented as
`tests/e2e/test_spike_ace2_recorded.py` on a five-instance fixture; see
`testing_tranche2b_spike_ace2_e2e.md`. The full-aggregation and live variants
remain open.*

The saved notebook outputs (also 2026-08-24) demonstrate exactly the failure modes
the tests must guard:

- 36 interfaces had partner order reversed (`7a91`, `7a92`, `7mjn`, `7rpv` assembly 4,
  and others).
- 21 antibody-bearing assemblies were excluded (`7l0n`, `7tn0`, `7xcp`, and others).
- Engineered mutations were mapped: `7bh9` chain E 477 → `P0DTC2` 477 role 2 (S477N);
  `7dmu` chains A and C 35 → `Q9BYF1` 35 role 1 (E35K); `7rpv` all four assemblies
  Q325Y on role 1, including the reversed assembly 4.
- A glycan contact was mapped: `7sn0` assembly 2, chain D 417, `NAG`.

The strongest single assertion here is that after `check_partner_consistency` every
record has `Q9BYF1` in role 1 and `P0DTC2` in role 2, and that the Q325Y mutation on
the reversed `7rpv` assembly still lands on role 1. That directly tests partner
reversal and annotation mapping on reversed records, which no unit test with
synthetic data can do as convincingly.

Stable assertions: those above, plus `ACE2:Y41-S:T500` and `ACE2:K353-S:G502`
hydrogen bonds present in the large majority of interfaces (129 and 128 of 130 at the
time; assert a fraction above 0.9 rather than a count); `7l0n` absent from the
records; and the rewiring table between the two largest states labelling
`ACE2:D38-S:Q498 salt_bridge` and `ACE2:S19-S:S477 hydrogen_bond` as higher in one
state and `ACE2:E37-S:Y505 hydrogen_bond` and `ACE2:D30-S:K417 salt_bridge` as higher
in the other.

Brittle: 130/114/13 counts, the 36 and 21 figures, which cluster is state 6 or 8, the
exact mutation row count of 26, and the ordering of the rewiring table. It is also
slow (thousands of ligand calls), so this belongs only in the live tier.

### Candidate C: `1spq` / `PDB-CPX-130029`, triosephosphate isomerase (cheap symmetry check)

The README and spec state its interface is invariant and forms one cluster at every
cut. That is the one existing piece of evidence that homodimer orientation is handled
adequately, so a live test asserting a single cluster at the default cut and across
the sweep is worth its few seconds.

### Recommended shape

Record the six responses for `PDB-CPX-172174` once into `tests/fixtures/recorded/`
and run the whole pipeline offline against them. On that frozen fixture the brittle
assertions become legitimate and exact, and the test is deterministic and fast. Keep
a second, live version of Candidate A with only the stable assertions, and add
Candidates B and C as live tests marked slow.

---

## 6. Proposed test architecture

```text
tests/
    conftest.py                     shared fixtures: synthetic records, recorded-fixture loader,
                                    session patch with zero backoff, autouse network block
                                    for anything not marked integration
    fixtures/
        synthetic/                  tiny hand-written payloads (a 3-interface heterodimer,
                                    a 2-interface homodimer, one reversed record, one Fab-bound assembly)
        recorded/PDB-CPX-172174/    real responses captured once: complex_details, interface_interactions,
                                    mutations, modifications, bound_molecules/*.json, ligand_interactions/*.json
    unit/                           offline, pure functions, synthetic inputs
        test_identifier_resolution.py
        test_assembly_selection.py
        test_interface_records.py
        test_partner_orientation.py
        test_similarity_clustering.py
        test_ligand_filtering.py
        test_annotation_overlap.py
        test_conservation_frequencies.py
        test_cluster_report.py
        test_rewiring.py
        test_export_json.py
    integration/
        test_api_mocked.py          retry, 404 policy, batching, order preservation; offline
        test_api_live.py            @pytest.mark.integration, one call per endpoint
    e2e/
        test_sting_recorded.py      full pipeline on recorded fixtures; exact assertions allowed
        test_sting_live.py          @pytest.mark.integration; stable assertions only
        test_spike_ace2_live.py     @pytest.mark.integration + slow; reversal, exclusion, mutation mapping
        test_tim_live.py            @pytest.mark.integration; single cluster
```

What belongs where:

- **unit/**: pure functions with synthetic inputs. Always offline. Any test here that
  touches the network is a bug in the test.
- **integration/**: the API layer. The mocked file exercises retry, 404 and batching
  logic against a patched session. The live file makes one real call per endpoint and
  checks response shape plus one known value.
- **e2e/**: the full pipeline. The recorded variant is offline and exact. The live
  variants are marked `integration` (and `slow` where relevant) and assert only
  behaviour that survives new depositions.

Configuration goes in a new `pytest.ini` or `pyproject.toml`: register the
`integration` and `slow` markers, and set the default `addopts` to
`-m "not integration"` so a plain `pytest` never touches the network. CI runs the
default; a nightly or manual job runs `-m integration`. Only `pytest` needs adding as
a development dependency, in a `requirements-dev.txt`. No mocking library is
required, since `monkeypatch` on `requests.Session.request` covers every case.

The notebook itself can get one optional smoke test:
`PDBE_INTERFACES_STATIC=1 jupyter nbconvert --execute` on the default configuration,
marked integration and slow. It proves the narrative cells still match the module
API, which the module tests cannot. The Mol* cell downloads CIF files, so this stays
out of the default run.

---

## 7. Priorities

Ordered by how easily the code could run green while producing a wrong scientific
answer.

1. **Residue correspondence in `_build_one`.** Role assignment, UniProt position,
   insertion codes, dropped pairs, duplicate collapse, microheterogeneity. Everything
   else consumes this.
2. **Partner reversal in `check_partner_consistency` and `_reverse_roles`.** Reversed
   record equals its unreversed twin at Jaccard 1.0; `author_to_uniprot` and
   `residue_identity` flip together; homodimers untouched; tie behaviour pinned.
3. **Annotation join in `overlap_annotations`.** Right residue, right role,
   chain-sensitive, insertion-code-sensitive, type filter applied, sentinel for
   unmapped interface residues. Then the same for `_parse_ligand_interactions`.
4. **Frequency and conservation arithmetic.** Once-per-interface counting, inclusive
   threshold boundary, denominators in `interface_frequency_summary`,
   `_contact_profile`, `_label_profile` (entry de-duplication) and the JSON export
   levels.
5. **Clustering on a designed matrix.** Membership at a known cut, single-interface
   case, index alignment between `flat_assignment` and `records`, rewiring labels at
   their boundaries.
6. **Assembly selection and identifier resolution** with mocked fetchers, including
   the two documented-but-unenforced rules (non-UniProt participant, malformed ids) as
   `xfail` so the gap is recorded.
7. **Mocked API behaviour**: retry set, 404 policy per endpoint, batching, tolerated
   versus fatal failures, thread-pool order preservation.
8. **Recorded-fixture end-to-end on STING**, with exact assertions.
9. **Live tier**: one call per endpoint, STING stable assertions, TIM single cluster,
   then Spike/ACE2 for reversal and exclusion.

### Observations the tests should decide rather than assume

Three things the read-through turned up. None needs fixing before tests are written,
but each deserves a test that states the intended behaviour.

- Homodimer orientation is not normalised and rests on PISA's chain ordering. The spec
  calls this intended; `1spq` forming a single cluster is the empirical evidence it
  holds in practice.
- An unknown complex id raises `requests.HTTPError` from `fetch_complex_details`
  rather than the `ValueError` the spec promises.
- `max_entries=0` produces an empty selection without an error.
