# Tranche 1.5: comparable-contact exclusion, partial mappings, homodimer diagnostic

Completed 2026-09-23, following `testing_tranche1_consolidation.md`. Three deliberately
small items:

1. Interface instances with zero comparable UniProt contact pairs are excluded from
   contact-based comparative analysis, with provenance retained.
2. Individually valid author-residue to UniProt mappings are preserved even when the
   other side of the contact is unmapped.
3. A diagnostic of homodimer orientation sensitivity on real data, without any change
   to homodimer behaviour.

Both behaviour changes were made test-first: the regression tests were written, run
and confirmed failing for the intended reason, then the smallest production change was
made. End-to-end regression tests were not started.

---

## Behaviour changes

### Zero-comparable-contact interfaces

A new `representation.select_comparable_records` splits records into those with at
least one UniProt contact pair and those with none. Records with none are held out of
the similarity matrix, clustering, conservation, frequency and rewiring stages, kept on
a `ComparableSelection.excluded` list, and tagged in `reasons` with
`no_comparable_uniprot_contacts`. There is no minimum count: one pair is enough.
`similarity.jaccard(set(), set())` is unchanged at 1.0.

The notebook's Phase 2 cell now performs the split, prints the exclusion summary, and
shows the per-interface summary table for all records. Phase 4 joins annotations over
all records, so excluded interfaces still receive mutations, modifications and
ligands. The structure table, cluster report, rewiring table and JSON export use the
comparable records only.

The scientific rule: an interface with PISA contacts but no comparable pair carries
insufficient comparable information. It is not "identical to another empty
interface", and a missing mapping must not look like biological variation.

### Partial UniProt mappings

In `_build_one`, residue mapping is now decided per side. A residue with both an
accession and a position enters `author_to_uniprot` and `residue_identity` whether or
not its contact partner maps. The pair still enters `uniprot_pairs` only when both
sides map, and the dropped-contact count is unchanged. Nothing is inferred: a side
missing either component stays unmapped, and annotations on it keep the
`(None, None, None)` sentinel. Microheterogeneity keeps its first-seen semantics, now
applied per side.

---

## Tests

```text
before             83 items / 57 functions
new tests added     4 items / 4 functions (plus 2 existing tests rewritten to the new rule)
final              87 items / 61 functions
passed             87
failed             0
skipped            0
xfail              0
runtime            1.8 s (2.6 s wall clock)
```

New or rewritten tests:

- `test_interface_records.py`
  - `test_contact_missing_uniprot_on_one_side` (rewritten): mapped side retained,
    unmapped side absent, pair excluded from comparison space.
  - `test_mapped_residue_is_consistent_across_partial_and_full_contacts` (new).
  - `test_select_comparable_records_excludes_zero_pair_interfaces` (new): comparable,
    empty and partially mapped records; reason recorded; provenance retained.
  - `test_excluded_record_does_not_reach_similarity_clustering_or_frequencies` (new):
    3 comparable records with contact X plus 1 empty record gives a 3 by 3 matrix,
    3 cluster labels and a frequency of 3/3.
- `test_annotations.py`
  - `test_annotation_on_partially_mapped_contact` (rewritten from the sentinel test):
    mutation on the mapped residue resolves to accession, position and role; mutation
    on the unmapped partner keeps the sentinel.
  - `test_half_mapped_interface_keeps_residue_mappings_but_is_not_comparable` (new,
    Part E): every contact half-mapped; residue mappings present, comparison pairs
    empty, record excluded, annotations still resolve.

All five RED tests were confirmed failing for the intended reason before the
production change. Two then failed on a pandas artefact (a frame mixing mapped and
unmapped rows stores the sentinel's numeric fields as NaN); the tests now normalise
NaN to None through a small helper. The network guard still passes and a plain
`pytest` makes no requests.

---

## Empty-interface example

Three interfaces carry the contact A50 to B106; a fourth has two PISA contacts whose
B side has no UniProt mapping.

```text
Excluded 1 of 4 interfaces from the contact-based comparison:
  4ddd a1 i1  no_comparable_uniprot_contacts  (2 author-space contacts, 2 residues mapped to UniProt)
Retained 3 interfaces for comparison.
similarity matrix shape: (3, 3)   cluster labels: [1, 1, 1]
contact P12345:K50 – P12345:D106   n_interfaces 3   fraction 1.0
```

Before this change the same input gave a 4 by 4 matrix, a fourth cluster label, and a
frequency of 0.75.

---

## Partial-mapping example

The excluded record above, with a mutation on each side of one contact:

```text
author_to_uniprot: {A50 -> (P12345, 50, role 1), A51 -> (P12345, 51, role 1)}
uniprot_pairs:     set()
excluded reason:   no_comparable_uniprot_contacts

auth_asym_id  auth_seq_id  unp_acc  unp_seq_id  role  mutation_label
A             50           P12345   50          1     K50A
B             106          None     NaN         NaN   D106A
```

The mapped residue is retained, the unmapped partner stays unmapped, the pair is
excluded from comparison, and the annotation on the mapped residue resolves correctly.

---

## Homodimer diagnostic

### Method

Script: `analysis/homodimer_orientation_diagnostic.py`. Outputs: `analysis/results/`
(`pairwise.csv`, `per_complex.csv`, `clustering_impact.csv`, `top_examples.csv`,
`skipped.csv`, `summary.json`). Live PDBe API calls; not part of the pytest suite;
no production behaviour touched.

For every pair of comparable interface instances (A, B) of a homodimer complex:

```text
direct  = Jaccard(A, B)
swapped = Jaccard(A, global_swap(B))     # copy 1 <-> copy 2 for the whole interface
delta   = swapped - direct
```

The swap is global across the interface; individual pairs are not canonicalised.
For complexes with any delta above 0.2, clustering at the notebook default cut of 0.6
was compared between the production similarity and a diagnostic orientation-aware
similarity `max(direct, swapped)`, which exists only in the script.

Sample: 30 homodimer complexes drawn from the PDBe-KB `assemblies_data.csv`
(composition string `{accession}_2`, 8 to 60 assemblies, seed 1), with STING
(`PDB-CPX-172174`) and triosephosphate isomerase (`PDB-CPX-130029`) forced in. Two
were skipped: one fully antibody-bound, one with a single comparable interface.
Reused existing functions throughout (`fetch_complex_details`,
`fetch_interface_interactions`, `select_interfaces`, `build_interface_records`,
`check_partner_consistency`, `select_comparable_records`, `jaccard`,
`cluster_interfaces`). Sampling needed no new infrastructure beyond one CSV read.

### Results

```text
complexes analysed                    28
interfaces                           474
pairwise comparisons                6479
fraction swapped > direct           4.4 %
fraction delta > 0.2                0.46 %   (30 pairs, 6 complexes)
fraction delta > 0.5                0.08 %   ( 5 pairs, 2 complexes)
complexes whose partition changes      3
complexes whose cluster count changes  2
complexes whose singletons decrease    2
```

The effect is concentrated in sparse fingerprints. Grouped by the smaller contact
count in each pair:

| contacts in sparser interface | pairs | delta > 0.2 | delta > 0.5 | max delta |
|---|---|---|---|---|
| 1 to 2 | 202 | 7.4 % | 2.5 % | 1.00 |
| 3 to 4 | 77 | 13.0 % | 0 | 0.50 |
| 5 to 9 | 600 | 0.17 % | 0 | 0.27 |
| 10 or more | 5600 | 0.07 % | 0 | 0.33 |

Clustering impact for the six complexes with any delta above 0.2:

| complex | interfaces | clusters current | clusters aware | singletons current | singletons aware | partition changed |
|---|---|---|---|---|---|---|
| `PDB-CPX-127157` Phosphite dehydrogenase | 21 | 1 | 1 | 0 | 0 | no |
| `PDB-CPX-118634` TetR-type regulator | 8 | 5 | 4 | 4 | 3 | yes |
| `PDB-CPX-151132` Anthranilate phosphoribosyltransferase | 11 | 1 | 1 | 0 | 0 | no |
| `PDB-CPX-172479` PARP15 | 51 | 2 | 2 | 1 | 1 | no |
| `PDB-CPX-152431` Motility protein B | 17 | 7 | 5 | 2 | 0 | yes |
| `PDB-CPX-138209` Tetracycline repressor D | 42 | 5 | 5 | 3 | 3 | yes (one interface moves) |

### Strongest examples

- **Motility protein B, `PDB-CPX-152431`**: one-contact interfaces flip between
  Jaccard 0 and 1 (`3imp` a3 vs a5, `3cyq` a1 vs a6, both same-entry assembly pairs).
  Clustering changes from 7 states with 2 singletons to 5 states with none. Median
  direct similarity is 0, so this complex already carries the sparse-fingerprint
  warning.
- **TetR-type regulator, `PDB-CPX-118634`**: same pattern (`5vlg` a3 vs a4); 5 states
  with 4 singletons becomes 4 with 3.
- **Tetracycline repressor D, `PDB-CPX-138209`**: `5fkn` vs `6yr2`, 17 and 19
  contacts, 0.38 direct to 0.71 swapped. One interface moves between clusters.
- **PARP15, `PDB-CPX-172479`**: `7pwk` and `7z2q` vs `8zv7`, 17 contacts each, 0.79
  direct to 1.00 swapped. A genuinely flipped identical interface with a rich
  fingerprint, though it does not change clustering at the 0.6 cut.
- **Phosphite dehydrogenase, `PDB-CPX-127157`**: `4e5n` vs `4e5p`, 22 contacts,
  0.63 to 0.83. No clustering change.
- **Anthranilate phosphoribosyltransferase, `PDB-CPX-151132`**: 10 of 55 pairs above
  0.2 (max 0.50); single cluster either way.
- **STING, `PDB-CPX-172174`** (the working example): max delta 0.15, no clustering
  change. **Triosephosphate isomerase**: max delta 0.07.

Six of the twelve pairs with delta of at least 0.5 are two assemblies of the same
entry, which is direct evidence that PISA does report the two copies in opposite
order within one deposition.

### Interpretation

Diagnostic finding. Orientation effect: **real but generally modest** for the clustering in this exploratory sample, and most consequential for sparse contact fingerprints. Swapping
increases similarity in about one pair in twenty, exceeds 0.2 in fewer than one in
two hundred, and moves cluster membership in three of 28 complexes. Where it changes
the state count, the interfaces involved have one or two contacts and would already
be flagged as sparse.

Two things temper that conclusion. The phenomenon is real for rich interfaces too
(PARP15, tetracycline repressor), and it can create or dissolve singleton states in
sparse-interface complexes. An interpretation note is therefore justified now; whether
to canonicalise orientation in production can be decided separately on these numbers.
Production homodimer behaviour is unchanged.

Answers to the five questions posed:

1. How often does a global swap increase similarity: about 4 percent of pairs.
2. How often substantially: under 0.5 percent of pairs above 0.2; under 0.1 percent
   above 0.5.
3. Does it change clustering materially: in 3 of 28 complexes, all with sparse
   interfaces among the affected members.
4. Are singleton states sometimes explained by orientation: yes, in two complexes,
   both sparse-fingerprint cases.
5. Is a note sufficient: on this evidence, yes, with the sparse-fingerprint caveat
   made explicit.

### Proposed homodimer note

Diagnostic finding, not production behaviour: the wording below was added to
`spec/new/final_spec.md` under *Known Limitation: Homodimer Orientation Dependence*,
with the numbers in a developer note that describes them as an exploratory sample.

> **Homodimer interfaces.** Contacts are compared using the partner orientation
> reported for each structural instance. Because the two partners have the same
> molecular identity, the partner-consistency check has nothing to act on, and
> equivalent homodimer interfaces reported with opposite partner orientations can
> appear less similar than they are. How much less depends on the interface: a fully
> symmetric contact set is unaffected, a partially symmetric one keeps the contacts
> that are their own mirror image, and an asymmetric one with a sparse fingerprint can
> lose all overlap. An exploratory diagnostic indicates that this has little effect on
> most interfaces but can influence clustering for sparse contact fingerprints.
> Singleton or small states in homodimers with few mapped contacts should therefore be
> interpreted with caution, and it is worth checking whether their members are
> assemblies of the same entry.

Production homodimer behaviour is unchanged: ordered contact tuples, the Jaccard
measure, the clustering method and the rewiring logic are as before. The figures above
describe this exploratory sample only and are not a PDB-wide prevalence estimate.

---

## Production files changed

- `pdbe_interfaces/representation.py`: per-side residue mapping in `_build_one`; new
  `NO_COMPARABLE_CONTACTS`, `ComparableSelection` and `select_comparable_records`.
- `notebook.ipynb`: Phase 2 guidance bullet and code cell (select comparable records,
  print the summary, summary table over all records); Phase 4 cell joins annotations
  over `all_records`. The stale saved outputs of those two cells were cleared, which
  accounts for most of the diff.

New, non-production: `analysis/homodimer_orientation_diagnostic.py`,
`analysis/results/` (six CSV and JSON files), and the four new tests listed above.

Not changed: Jaccard semantics, homodimer comparison, rewiring and annotation logic
for homodimers, identifier validation, `spec/new/final_spec.md`. Nothing is
committed. The end-to-end regression tranche has not been started.
