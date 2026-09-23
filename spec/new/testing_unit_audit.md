# Unit-suite audit before commit

Final review of `tests/unit/` before the testing milestone is committed.
Completed 2026-09-23.

Scope was `tests/unit/` only. Nothing was modified: not the unit tests, not
`tests/e2e/`, not `tests/smoke/`, not the recorded fixtures, not production
code, not `notebook.ipynb`, not scientific behaviour. This document is the
proposed consolidation report; no change follows from it.

The purpose was not to reach a smaller number.

**Current state: 62 test functions, 87 items, all passing.**

---

## Classification scheme

Each logical test group was classified as one of:

1. **scientifically essential** — protects behaviour whose failure would change
   the biological interpretation;
2. **software-contract essential** — protects API, error or input behaviour the
   workflow genuinely relies on;
3. **useful edge-case regression** — narrower, but guards a real hazard;
4. **redundant with another test**;
5. **implementation-detail / overly coupled**.

---

## `test_interface_records.py` — 11 functions, 17 items

| Test | Class | Contract |
|---|---|---|
| `one_contact_yields_author_and_uniprot_keys_with_roles` | 1 | the base record shape: author key, UniProt key, role, residue identity, counters |
| `comparison_space_is_keyed_by_uniprot_not_author_numbering` | 1 | equivalent residues in different entries converge; unrelated ones stay apart |
| `role_follows_pisa_side_not_accession_order` | 1 | role comes from the PISA side, not from alphabetical accession order |
| `insertion_code_normalisation` (4 params) | 3 | `strip() or None` across the forms endpoints emit |
| `duplicate_observations_collapse_but_bond_types_stay_distinct` | 1 | atom-level collapse and typed/untyped semantics at the **builder** layer |
| `contact_missing_uniprot_on_one_side` | 1 | per-residue mapping independence (tranche 1.5) |
| `mapped_residue_is_consistent_across_partial_and_full_contacts` | 3 | guards a hazard the tranche-1.5 per-side loop introduced: a residue first seen in a half-mapped contact must not later collide with itself |
| `partner_accessions_come_from_first_contact_carrying_both` (3 params) | 3 | accession inference across three distinct input shapes |
| `microheterogeneity_keeps_first_seen_author_key` | 1 | first-seen rule, now applied per residue |
| `select_comparable_records_excludes_zero_pair_interfaces` | 1 | the 0 versus ≥1 rule, the recorded reason, provenance retention |
| `excluded_record_does_not_reach_similarity_clustering_or_frequencies` | 1 | denominator integrity: a contact in three of three comparable interfaces reads 3/3, not 3/4 |
| `one_record_per_interface_in_input_order` | 2 | cluster labels are matched to records by position, so order must hold |

## `test_partner_orientation.py` — 5 functions, 8 items

All class 1.

| Test | Contract |
|---|---|
| `minority_ordering_is_reversed_onto_the_canonical_roles` | reversal moves contacts, `author_to_uniprot`, `residue_identity` and author pairs together; majority records untouched |
| `two_record_tie_keeps_the_first_seen_ordering` | documented tie behaviour |
| `records_matching_neither_ordering_are_left_alone` (2 params) | third accession, and no accessions at all |
| `homodimer_similarity_depends_on_chain_order` (3 params) | flipped asymmetric 0.0, partially symmetric 0.2, fully symmetric 1.0 |
| `homodimer_flipped_chain_order_separates_at_default_cut` | the observable consequence: equivalent homodimer interfaces separate into different states |

## `test_similarity.py` — 6 functions, 10 items

All class 1: `jaccard` (5 params, including the documented empty-vs-empty 1.0),
matrix symmetry and unit diagonal, clustering membership, the cut in
`1 − Jaccard` units, the single-interface case, and typed versus untyped at the
**similarity** layer.

## `test_annotations.py` — 12 functions, 14 items

Class 1 unless noted.

| Test | Contract |
|---|---|
| `annotations_map_to_uniprot_key_and_role` | mutation on role 1, modification on role 2 |
| `annotation_on_a_non_interface_residue_is_dropped` (4 params) | wrong chain, off-interface, missing insertion code, blank insertion code |
| `insertion_code_is_part_of_the_join` | the positive case |
| `annotation_on_partially_mapped_contact` | mapped residue resolves; unmapped partner keeps the `(None, None, None)` sentinel |
| `conflict_mutations_excluded_by_default_and_included_on_request` | mutation-type filtering |
| `mutation_on_reversed_record_lands_on_canonical_role` | partner reversal carries the annotation join with it |
| `filter_bound_molecules_blocklist_and_glycan_flag` | ligand blocklist and glycan flag |
| `ligand_contacts_collapse_atoms_and_only_fetch_surviving_ligands` | atom collapse; blocked ligands are never fetched |
| `ligand_contact_maps_to_correct_interface_residue` | ligand contact resolves to the right interface residue |
| `half_mapped_interface_keeps_residue_mappings_but_is_not_comparable` | **the interaction** (see below) |
| `empty_responses_give_empty_frames` | class 2: no annotations is not an error |
| `one_row_per_interface_carrying_the_residue` | class 2: per-interface duplication is intentional |

## `test_frequencies.py` — 6 functions, 9 items

All class 1: denominators and the pair matrix, inclusive threshold (4 params),
typed conservation, export levels with typed versus untyped, rewiring labels,
and cluster-report in/out counts.

## `test_identifiers.py` (8 functions, 16 items) and `test_api.py` (12 functions, 12 items)

All class 2, each covering a distinct branch: identifier normalisation for both
forms, empty input rejected before any request, dimer rejection, the three
multiple-complex resolution branches, the no-complex cases, the two 404
policies, status versus exception retry, exhaustion on tolerant versus critical
endpoints, non-retryable 4xx, batched POST chunking and merging, empty input
making no request, parallel-fetch ordering, and the network guard itself.

---

## Nothing in class 4 or class 5

No test was found that is redundant with another, and none asserts internal
mechanics without an externally meaningful contract.

---

## The tranche-1.5 overlap, examined explicitly

The tranche-1.5 additions do overlap at the **assertion** level, and one test
deserves explicit defence.

`test_half_mapped_interface_keeps_residue_mappings_but_is_not_comparable`
re-asserts pieces that three other tests already cover:

| Assertion it makes | Already covered by |
|---|---|
| mapped sides retained in `author_to_uniprot` | `test_contact_missing_uniprot_on_one_side` |
| `uniprot_pairs` empty, excluded with the reason | `test_select_comparable_records_excludes_zero_pair_interfaces` |
| annotation resolves, unmapped partner gets the sentinel | `test_annotation_on_partially_mapped_contact` |

It survives because it is the only test asserting the **conjunction**, which is
a distinct scientific contract: *an interface can carry useful residue-level
annotation information while being unusable for contact-based cross-instance
comparison.* Each of the other three asserts one half. Removing it would leave
that interaction unprotected, and the interaction was the explicit purpose of
Part E of tranche 1.5.

Classified **1, scientifically essential**, with overlapping assertions
acknowledged.

---

## Candidates considered and rejected

Four reductions are technically available, totalling four items (87 → 83). Each
is stated in full, with the reason for rejecting it.

### 1. `test_pdb_id_with_no_complex_raises`: merge two params

- **Current**: three params — a 404 (`_get_json` returns `None`), an empty list
  for the queried id, and a record with no `pdb_complex_id`.
- **Contract**: an entry that resolves to nothing raises `ValueError` before any
  complex-details call.
- **Duplication**: the 404 and empty-list params land in the same
  `if not data or pdb_id not in data or not data[pdb_id]` branch.
- **Remaining protection**: either param alone covers the branch; the third
  param covers a different one.
- **Change**: −1 item.
- **Rejected because**: the two params document two genuinely different PDBe
  behaviours, a 404 and an empty 200, which a reader would otherwise have to
  infer.

### 2. `test_insertion_code_normalisation`: merge two params

- **Current**: four params — `None`, `""`, `" "`, `" B "`.
- **Contract**: insertion codes normalise via `strip() or None`.
- **Duplication**: `""` and `" "` both resolve through the same expression.
- **Remaining protection**: either alone covers it; `None` and `" B "` cover the
  early return and the padded-code path.
- **Change**: −1 item.
- **Rejected because**: the spec records both as forms the endpoints actually
  emit, so both are documentation of real input.

### 3. `test_conserved_residues_threshold_is_inclusive`: drop the `0.0` param

- **Current**: four params — `1.0`, `2/3`, `0.67`, `0.0`.
- **Contract**: threshold comparison is inclusive at exact equality.
- **Duplication**: the discrimination is carried by `2/3` (included) against
  `0.67` (excluded); `1.0` and `0.0` are boundary sanity.
- **Remaining protection**: the middle two params.
- **Change**: −1 item.
- **Rejected because**: `threshold=0.0` returning everything seen at least once
  is a real configuration outcome, cheap to keep.

### 4. `test_empty_responses_give_empty_frames`: fold into another test

- **Current**: one item asserting all three annotation frames are empty when the
  responses are empty.
- **Contract**: absent annotations are not an error.
- **Duplication**: partly implied by the drop tests, which assert
  `.mutations.empty`.
- **Remaining protection**: `test_annotation_on_a_non_interface_residue_is_dropped`.
- **Change**: −1 item.
- **Rejected because**: the spec states this explicitly under failure handling
  ("No mutations or no ligands returned. Not an error"), and a named test is the
  clearest place for it. The drop tests exercise a different branch, one where
  the loop runs and matches nothing.

---

## Recommendation

**Leave `tests/unit/` unchanged at 87 items.**

The consolidation in tranche 1 already did this work: 128 items were reduced to
83 by removing 34 redundant variants and 12 implementation-detail tests, and by
collapsing 18 cases into 9 parametrised tests. The four items added since were
targeted at behaviours that did not exist at that point:

- per-residue mapping independence;
- the zero-comparable-contact exclusion rule;
- its downstream effect on frequency denominators;
- the interaction between exclusion and annotation usability.

All eleven contracts named for protection are covered, each by a test whose
failure would name the specific contract broken:

| Contract | Protected by |
|---|---|
| residue correspondence | `test_interface_records.py`, first three tests |
| partial mappings | `contact_missing_uniprot_on_one_side`, `mapped_residue_is_consistent_across_partial_and_full_contacts`, `annotation_on_partially_mapped_contact` |
| zero-comparable-contact exclusion | `select_comparable_records_excludes_zero_pair_interfaces`, `excluded_record_does_not_reach_similarity_clustering_or_frequencies` |
| heterodimer partner reversal | `minority_ordering_is_reversed_onto_the_canonical_roles`, `mutation_on_reversed_record_lands_on_canonical_role` |
| homodimer orientation behaviour | `homodimer_similarity_depends_on_chain_order`, `homodimer_flipped_chain_order_separates_at_default_cut` |
| typed vs untyped contact semantics | `duplicate_observations_collapse_but_bond_types_stay_distinct` (builder), `typed_and_untyped_diverge_only_on_bond_type` (similarity), `conserved_pairs_are_typed`, export test |
| annotation mapping | `test_annotations.py`, nine tests |
| frequency denominators | `residue_and_pair_frequencies_use_interfaces_as_denominator`, `cluster_report_in_and_out_counts` |
| conservation / rewiring arithmetic | `conserved_residues_threshold_is_inclusive`, `conserved_pairs_are_typed`, `rewiring_labels_from_cluster_fractions` |
| identifier resolution | `test_identifiers.py`, eight tests |
| API retry / failure behaviour | `test_api.py`, twelve tests |

No test was removed on the grounds that an end-to-end test also exercises the
code. Unit and end-to-end tests serve different purposes: the unit tests name
the broken contract, the end-to-end tests prove the assembled workflow still
reproduces known biology.

---

## Status

Nothing was changed by this audit. The suite stands at 98 tests (87 unit, 10
end-to-end, 1 smoke), all passing.

The testing milestone remains **uncommitted**. An earlier commit of it was made
and then reversed with `git reset HEAD~1` at the user's request, which preserved
every file; nothing was pushed, and `origin/main` is still at `630c891`.

This document is itself untracked and should be included in the milestone commit
when it is made.
