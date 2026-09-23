# Tranche 2B final report: STING close-out edits and Spike RBD–ACE2 regression

Completed 2026-09-23. Two sequential tasks: three interpretation and
documentation corrections to the STING material, then a focused frozen
Spike RBD–ACE2 end-to-end regression.

The detailed tranche document is `testing_tranche2b_spike_ace2_e2e.md`. This
document is the close-out summary.

97 tests pass offline. No production code was changed.

---

## Part A — STING close-out

### Aggregation-rule wording

`final_spec.md` did **not** contain the rule. Line 1167 addressed out-of-scope
sequence conservation ("this workflow uses structural instances of the same
complex, not sequence homologues"), which is a different point about what the
workflow does not do, not a statement about how complexes are aggregated.

A new paragraph was added to the domain-context section, immediately before
*"Why a single complex has many structures"*. It was first written in terms of
species and **superseded on review**, because a complex may legitimately span
organisms. The current wording, in `final_spec.md`, is stated in terms of mapped
component identity and stoichiometry rather than species; the correction is
recorded in `testing_notebook_smoke.md`.

Both STING literature documents state the rule and note that human comparisons
serve only to relate the mouse observations to the literature and to test
transferability, playing no part in producing the states.

### D273–H156 reframing

Three passages were rewritten. The contact is now presented as a valid,
reproducible state-associated feature *within* mouse STING which, because the
aggregation contains only mouse structures, cannot have arisen from a sequence
difference between orthologues. What fails to carry over is the
**interpretation**, since mouse D273 corresponds to human Y274.

The distinction maintained throughout is:

```text
variation among mouse STING structures   <- what the workflow measured; sound
cross-species transferability            <- absent here; do not generalise
```

Wording that could imply species heterogeneity confounded the clustering has
been removed; the phrase "no species heterogeneity inside the dataset" no longer
appears anywhere in `spec/new/`.

### Recommendation wording

The restrictive "reproduction rather than discovery" phrasing was replaced with:

> **STING is primarily a validation case study: the aggregated interface
> analysis recovers established conformational biology while expressing it as
> residue-level contact rewiring. It also highlights specific contact-level
> observations not identified explicitly in the literature searched, which
> require further validation before any novelty claim.**

A following paragraph preserves the four evidence classes and assigns the
observations to them: **established biology reproduced** (CDN/mimetic state
separation, SR-717 co-clustering); **mechanistic interpretation supported**
(state-defining contacts in the lid region, ligand contacts on canonical pocket
residues, c-di-GMP separating from other cyclic dinucleotides);
**workflow-derived observation** (state-dependent A232–D209, A232–Y260,
D209–G233; the `4lol` partial lid-contact network; D273–H156); and **requires
further validation**, applying to every item in the third class.

### "partial lid engagement" replaced

Changed to **"partial lid-contact network"** in all five occurrences across the
two documents, with an added note that the wording is deliberately about which
contacts are present and absent and asserts no particular global structural
mechanism.

### Files changed in Part A

```text
spec/new/final_spec.md
spec/new/testing_tranche2a_sting_literature.md
spec/new/testing_tranche2a_sting_literature_report.md
```

No code, tests or fixtures were touched.

---

## Part B — Spike RBD–ACE2 fixture

`6m0j` resolves to `PDB-CPX-140195`, a heterodimer of ACE2 and Spike. The full
aggregation holds 131 comparable interface instances across 115 entries, of
which **36 are reported by PISA with the two partners in the opposite order**.
Three entries carry both orientations across their own assemblies: `7p19`,
`7rpv` and `8xye`.

Five instances from three entries were selected, 200 KB in
`tests/fixtures/recorded/spike_ace2_reversal/`:

| Instance | Order | Why selected |
|---|---|---|
| `6m0j` a1 i1 | canonical | archetypal Spike RBD–ACE2 structure, the named starting entry |
| `7p19` a1 i1 | canonical | Spike Q498Y engineered mutation at the interface (chain E) |
| `7p19` a2 i1 | **reversed** | same entry as a1, opposite reported order |
| `7rpv` a1 i1 | canonical | same entry and same biological ACE2-RBD interface as the reversed assembly; a strong control for partner-order reversal, while allowing genuine differences in the contacts detected in each copy (a4) |
| `7rpv` a4 i1 | **reversed** | ACE2 Q325Y engineered mutation at the interface (chain D) |

Three canonical against two reversed gives `check_partner_consistency` an
unambiguous majority. A 2:2 split would have tested the first-seen tie-break
rather than the majority rule, which is why five instances were taken rather
than four.

Files: `complex_details_6m0j.json`, `complex_details.json` (verbatim, all 151
assemblies), `interface_interactions.json` (**subset** to the five instances,
recorded in the file and in metadata), `mutations_post.json`,
`modifications_post.json` (recorded 404), `bound_molecules.json`,
`ligand_interactions.json` (14 instances), `metadata.json`. Captured once on
2026-09-23; 22 requests recorded; frozen and not refreshed automatically.

---

## Canonical roles

Established from PDBe entity metadata and UniProt, not from the workflow's own
output, then asserted:

| Role | Biological partner | UniProt | Organism |
|---|---|---|---|
| 1 | Angiotensin-converting enzyme 2 (ACE2) | `Q9BYF1` | *Homo sapiens* |
| 2 | Spike glycoprotein (RBD) | `P0DTC2` | SARS-CoV-2 |

Chain assignment from the PDBe `molecules` endpoint:

```text
6m0j   chain A = ACE2      chain E = Spike
7p19   chains A, B = ACE2  chains C, E = Spike
7rpv   chains A-D = ACE2   chains E-H = Spike
```

So the reversed instance `7rpv` a4 uses chain D for ACE2 and chain H for Spike,
which is exactly what `author_to_uniprot` reproduces after normalisation.

---

## Contact regression

ACE2 Tyr41 with Spike Thr500, a hydrogen bond, in the normalised representation:

```text
before   7rpv a1 (canonical): ((Q9BYF1,41,1),(P0DTC2,500,2))   present
         7rpv a4 (reversed) : ((P0DTC2,500,1),(Q9BYF1,41,2))   present, flipped
after    both instances     : ((Q9BYF1,41,1),(P0DTC2,500,2))
```

Two further established ACE2–RBD contacts are asserted the same way: ACE2 Asp30
with Spike Lys417 (salt bridge) and ACE2 Asp38 with Spike Gln498 (hydrogen
bond). Residue identities are asserted as well (`D30`, `Y41`, `K417`, `T500`),
so the mapping is checked against residues rather than against numbers alone.
The two assemblies share 11 contacts once comparable.

### Similarity

| Pair | Jaccard before | Jaccard after |
|---|---|---|
| `7rpv` a1 vs a4 | 0.0000 | 0.5500 |
| `7p19` a1 vs a2 | 0.0000 | 0.6923 |

Not forced to 1.0: the assemblies genuinely differ in a few detected contacts,
as independent copies of a crystallographic assembly normally do. The
**before** value of exactly 0.0 is the scientifically important number, since it
shows that two views of the same biological interface share no contact tuple at
all when reported in opposite orders, and would cluster apart purely because of
a reporting convention.

---

## Annotation regression

| Entry / assembly | Orientation | Author chain, residue | Resolves to | Mutation |
|---|---|---|---|---|
| `7rpv` a4 | **reversed** | D 325 | `Q9BYF1` 325, role 1 (ACE2) | Q325Y |
| `7rpv` a1 | canonical | A 325 | `Q9BYF1` 325, role 1 (ACE2) | Q325Y |
| `7p19` a1 | canonical | E 498 | `P0DTC2` 498, role 2 (Spike) | Q498Y |

The central assertion is the first row: a mutation on a reversed instance,
sitting on author chain D, resolves to ACE2 Gln325 in role 1 and is not
attributed to Spike residue 325. The canonical assembly of the same entry is the
control. The canonical UniProt residue identity at `(Q9BYF1, 325, 1)` is `Q`,
consistent with the `Q325Y` label. The `7p19` row shows the other partner
mapping to role 2.

---

## Tests

```text
tests before      92
new tests          5
tests after       97

passed            97
failed             0
skipped            0
xfail              0

runtime           2.1 s for the suite (2.8 s wall clock); Spike/ACE2 about 0.05 s
```

Five tests in `tests/e2e/test_spike_ace2_recorded.py`:

```text
test_spike_ace2_recorded_normalises_partner_orientation
test_spike_ace2_recorded_reverses_every_view_of_a_record_together
test_spike_ace2_recorded_preserves_equivalent_contacts_after_reversal
test_spike_ace2_recorded_similarity_is_not_broken_by_partner_order
test_spike_ace2_recorded_preserves_annotation_partner_identity
```

Fully offline. All 22 requests go through the recorded router, which raises on
anything unrecorded, and the root `block_network` guard is unchanged.

---

## Production changes

None. `pdbe_interfaces/` and `notebook.ipynb` were not touched in this tranche.

One **test-plumbing** change: `tests/e2e/conftest.py` previously assumed the
fixture directory was named after the complex id and that the entry-lookup file
was `complex_details_11gl.json`. It now takes a fixture directory name, reads the
complex id from `metadata.json`, and discovers the entry-lookup file by glob.
This changes no behaviour for the STING tests, which continue to pass unchanged.

---

## Unexpected findings

- **`7p19` assembly 2 contributes no mutation row.** The Q498Y record exists for
  author chain C, the Spike chain in that assembly, but residue 498 is not an
  interface residue there, so the workflow correctly reports nothing. The
  equivalent residue *is* at the interface in assembly 1. This is a real
  difference between two assemblies of one deposition, not a mapping failure,
  and it is why the Spike-side annotation assertion uses the canonical assembly.
- **`7rpv` carries five engineered ACE2 mutations** (L79F, M82Y, Q325Y, H374A,
  H378A); only Q325Y lies at the interface and only it is reported. The
  H374A/H378A pair are the catalytic-site substitutions used to inactivate ACE2
  and sit far from the RBD interface, so the filtering behaves correctly.
- **The shared contact core is the textbook ACE2–RBD interface.** The 11
  contacts shared between `7rpv` a1 and a4 involve ACE2 residues 30, 35, 38, 41,
  42 and Spike residues 417, 493, 498, 500 among others. No literature review
  was performed for this tranche, but the correspondence is an informal sanity
  check on the mapping.

---

## Not implemented, by instruction

Full Spike/ACE2 clustering survey, live PDBe tests, notebook smoke tests, CI,
STING conformational-metric analysis, human STING analysis, and Spike/ACE2
literature grounding. Enough authoritative metadata was recorded to establish
the biological identities; whether a separate literature-grounding analysis is
useful can be decided separately.

---

## Document sequence

```text
spec/new/testing_strategy.md                            long-term plan
spec/new/testing_tranche1_report.md                     first offline suite, design questions
spec/new/testing_tranche1_consolidation.md              suite streamlined to 83 items
spec/new/testing_tranche1_5_report.md                   exclusion rule, partial mappings, homodimer diagnostic
spec/new/testing_tranche1_5_closeout.md                 spec and consistency update
spec/new/testing_tranche2a_sting_e2e.md                 frozen STING end-to-end regression
spec/new/testing_tranche2a_final_report.md              tranche 2A close-out
spec/new/testing_tranche2a_sting_literature.md          STING literature grounding (detailed)
spec/new/testing_tranche2a_sting_literature_report.md   STING literature grounding (summary)
spec/new/testing_tranche2b_spike_ace2_e2e.md            frozen Spike RBD-ACE2 regression
spec/new/testing_tranche2b_final_report.md              this document
```
