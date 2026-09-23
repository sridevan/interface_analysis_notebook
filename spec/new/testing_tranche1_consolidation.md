# Testing tranche 1: suite consolidation report

Follow-up to `testing_tranche1_report.md`, completed 2026-09-23. The 128-test
offline suite from tranche 1 was audited and streamlined to a smaller,
high-signal suite that protects every scientifically important behaviour while
removing redundancy, repeated permutations, and tests coupled to implementation
detail. No production code was changed. No design question was resolved.

The retained suite runs offline in under two seconds. All 83 items pass.

---

## Audit of the 128 tests

Every test was classified before anything was changed.

| Category | Tests | Action |
|---|---|---|
| Scientifically essential | 41 | Kept, many merged into fewer functions with more assertions |
| Important software contract | 23 | Kept the distinct branches, dropped per-endpoint repeats |
| Redundant variant | 34 | Removed or folded into an existing test |
| Implementation-detail test | 12 | Removed |
| Candidate for parametrisation | 18 | Collapsed into 9 parametrised tests |

The main patterns found:

- seven separate partner-reversal tests that each checked one field of the same
  three-record scenario;
- six single-assertion Jaccard tests;
- per-endpoint 404 tests and per-status retry tests that all enter one branch;
- identifier normalisation variants that differed only by whitespace or case.

---

## Before / after

```text
                              before     after    (functions)
interface correspondence        21        14         9
partner orientation             16         8         5
similarity/clustering           14        10         6
annotations                     20        14        11
frequency/conservation          10         9         6
identifiers                     24        16         8
API                             23        12        12
TOTAL                          128        83        57
```

Counts are pytest items. Parametrised cases inflate the item count relative to
functions; 57 functions is the number a reader has to understand.

---

## What was consolidated

- **7 reversal tests → 1.** `test_minority_ordering_is_reversed_onto_the_canonical_roles`
  runs the three-record heterodimer scenario once and asserts, in order: Jaccard 0
  before the check, same accession in the same role, identical tuples and Jaccard
  1.0, `author_to_uniprot` on new roles, author pairs swapped, `residue_identity`
  on new roles, majority records untouched.
- **4 homodimer tests → 1 parametrised + 1.** Flipped asymmetric (0.0), flipped
  partially symmetric (0.2) and fully symmetric (1.0) are three cases of
  `test_homodimer_similarity_depends_on_chain_order`, whose docstring states the
  semantics are an open design question. A separate test keeps the observable
  consequence: the flipped pair separates at the default cut of 0.6.
- **6 Jaccard tests → 1 parametrised** (`test_jaccard`) with five cases, symmetry
  asserted inside it.
- **Mutation and modification join → 1 test** asserting a role-1 mutation and a
  role-2 modification from one fixture.
- **4 "does not match" cases → 1 parametrised**: wrong chain, off-interface,
  missing insertion code, blank insertion code.
- **3 partner-accession inference tests → 1 parametrised.**
- **4 conservation-threshold assertions → 1 parametrised** (all, exact equality,
  just above, zero).
- **2 export tests → 1** covering conservation levels, once-per-interface
  counting, typed versus untyped, and the file written to disk.
- **Identifier normalisation: 3 + 4 + 3 items → 3 + 3 + 2**, the two ambiguity
  outcomes into one parametrised test, and three "no complex" payloads into one
  parametrised test including the missing-complex-id case.
- **API: batched chunking, merge and 404-chunk tolerance → 1 test.** HTTP-status
  and exception exhaustion → 1 test. Timeout and connection-error retry → 1 test.
  404 and empty-result on the critical endpoint → 1 test.

---

## What was removed

### Implementation-detail tests

- `ins_code_1` fallback field; missing `bond_type` becomes `"unknown"`; contact
  with missing author fields skipped silently (recorded as design question 7 in
  the tranche-1 report).
- Duplicate complex ids collapsing to one in the resolver.
- `fetch_bound_molecules_many` keyed by id; parallel fetch with no items; entries
  without ligands still getting a key; unknown entry in the mutation response
  ignored; malformed mutation records skipped; empty blocklist keeps everything.
- Cluster report core contacts inclusive threshold (same `>=` rule already tested
  through `conserved_residues`); empty record list has no conservation.

### Redundant variants

- Whitespace insertion codes collapsing (covered by normalisation plus the
  duplicate-collapse test).
- Cluster labels follow input order (covered by the designed clustering test and
  the record-order test).
- Everything merges above cut 1.0 (a permutation of the cut-height test).
- Reversed record clusters with its equivalents (follows from Jaccard 1.0 plus
  the clustering test).
- Homodimer same orientation is identical (trivial).
- Off-interface modification (same branch as mutation).
- 500 as a second retry status; 404 on a second best-effort endpoint;
  interface-interactions success as a second success case; timeout exhaustion
  tolerated (same `except Exception` branch as HTTP exhaustion tolerated).

### Removed for policy reasons

- 404 on complex details raises `HTTPError`, and 404 on the entry lookup becomes
  `ValueError`. The first pinned an exception type listed as undecided; the
  second was already covered through the identifier tests. The design question
  stays in the tranche-1 report.
- Non-dimer accepted when `require_dimer=False`: the notebook always passes
  `True`.
- "Complex id that does not resolve": on inspection it exercised the test's own
  mock, not production code.

---

## Scientific protection

The retained suite still protects each of these, with the test that carries it:

- **Residue correspondence**:
  `test_one_contact_yields_author_and_uniprot_keys_with_roles`,
  `test_comparison_space_is_keyed_by_uniprot_not_author_numbering`,
  `test_role_follows_pisa_side_not_accession_order`,
  `test_microheterogeneity_keeps_first_seen_author_key`.
- **Heterodimer reversal**:
  `test_minority_ordering_is_reversed_onto_the_canonical_roles`, plus
  `test_mutation_on_reversed_record_lands_on_canonical_role` for annotation roles.
- **Homodimer orientation issue**:
  `test_homodimer_similarity_depends_on_chain_order` (flipped asymmetric gives
  Jaccard 0) and `test_homodimer_flipped_chain_order_separates_at_default_cut`.
- **Annotation role mapping**:
  `test_annotations_map_to_uniprot_key_and_role`,
  `test_annotation_on_a_non_interface_residue_is_dropped`,
  `test_ligand_contact_maps_to_correct_interface_residue`.
- **Missing mappings**:
  `test_contact_missing_uniprot_on_one_side` (mapped side not entered into
  `author_to_uniprot`),
  `test_mutation_on_interface_residue_without_uniprot_uses_sentinel`, and the
  empty-vs-empty case of `test_jaccard`.
- **Insertion codes**:
  `test_insertion_code_normalisation` and `test_insertion_code_is_part_of_the_join`.
- **Contact similarity**:
  `test_jaccard`, `test_similarity_matrix_is_symmetric_with_unit_diagonal`,
  `test_typed_and_untyped_diverge_only_on_bond_type`,
  `test_duplicate_observations_collapse_but_bond_types_stay_distinct`.
- **Clustering**:
  `test_clustering_groups_identical_and_separates_disjoint`,
  `test_cut_is_in_one_minus_jaccard_units`, `test_single_interface_is_one_cluster`.
- **Frequency denominators**:
  `test_residue_and_pair_frequencies_use_interfaces_as_denominator`,
  `test_export_json_frequencies_levels_and_typed_vs_untyped`,
  `test_cluster_report_in_and_out_counts`.
- **Conservation thresholds**:
  `test_conserved_residues_threshold_is_inclusive`, `test_conserved_pairs_are_typed`.
- **API retry and failure**:
  `test_transient_status_is_retried_then_succeeds`,
  `test_network_exception_is_retried_then_succeeds`,
  `test_exhausted_retries_raise_on_critical_endpoint`,
  `test_exhausted_retries_are_tolerated_on_best_effort_endpoint`,
  `test_non_retryable_client_error_raises_immediately`, the two 404-policy tests,
  and the batched POST tests. Parallel ordering was kept because
  `collect_ligand_contacts_for_entries` zips requests with responses by position.

The network-blocking fixture in `tests/conftest.py` is unchanged and
`test_network_is_blocked_by_default` still guards it.

---

## Results

```text
passed   83
failed   0
skipped  0
xfail    0
```

## Runtime

About 1.9 seconds for the tests, 2.6 seconds wall clock including interpreter and
import startup. The single warning is the pydantic deprecation from `molviewspec`
on import.

## Production files

None modified. `git status` shows only the untracked `pytest.ini`,
`requirements-dev.txt`, `tests/`, and the documents under `spec/new/`. Nothing is
committed. The end-to-end tranche has not been started.
