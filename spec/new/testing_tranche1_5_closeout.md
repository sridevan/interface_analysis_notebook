# Tranche 1.5 close-out: specification and consistency update

Completed 2026-09-23, following the review of `testing_tranche1_5_report.md`. This was
a documentation and consistency task only. The two production behaviour changes from
tranche 1.5 were accepted as implemented:

1. Interface instances with zero comparable UniProt contact pairs are excluded from
   contact-based comparative analysis while provenance is retained.
2. Individually valid author-residue to UniProt mappings are retained independently
   even when the interaction partner is unmapped.

The homodimer diagnostic was accepted as evidence that the orientation issue is real
but does not currently justify redesigning homodimer comparison.

No code changed in this task. Tranche 2 (end-to-end regression) has not been started.

---

## Specification changes

All in `spec/new/final_spec.md`, following its existing section structure
(77 insertions, 7 deletions).

- **Revision note**: three bullets pointing to the new rules and the new known
  limitation.
- **Warning conditions**: the "residues lacking a UniProt mapping" bullet now
  describes per-contact dropping with per-residue mapping retained, plus a new bullet
  for interfaces with no comparable pair.
- **Edge cases**: the homodimer bullet now states that the partner-consistency check
  cannot act on homodimers and points to the known limitation. The microheterogeneity
  bullet, which still said the workflow does not detect collisions, now describes the
  implemented first-seen rule, applied per residue.
- **Execution flow, Phase 2**: step 5 rewritten as independent residue mapping with
  both-sides-required pair comparability, and a new step 6a for the comparable split.
- **New subsection "Residue correspondence and contact-pair comparability"** after
  the aggregation rule: the two-level rule, the A50/B106 worked example, the Phase 4
  consequence, and the no-inference statement.
- **New subsection "Comparable-interface eligibility"**: the 0 versus at least 1 rule,
  no other minimum, the seven exclusion points, provenance retention, the
  methodological reason (missing correspondence is neither variation nor identity),
  and the note that `jaccard(set(), set())` stays 1.0 but such records never reach it.
- **Failure handling**: the stale "flagged in the structure table" bullet replaced
  with the actual behaviour.
- **Module `representation`**: `select_comparable_records` and the
  `ComparableSelection` dataclass added.
- **User-facing copy**: Phase 2 preamble and the mapping warning rewritten; a new
  exclusion-summary warning template added.
- **New section "Known Limitation: Homodimer Orientation Dependence"** after the
  disulfide limitation.

### Comparable-interface eligibility rule (as documented)

An interface instance is eligible for contact-based cross-instance comparison only
when at least one of its residue–residue contacts is represented in the common
UniProt coordinate system:

```text
0 comparable pairs   → excluded from the comparative analysis
>= 1 comparable pair → retained
```

No other minimum is applied. An excluded instance is retained for provenance with the
reason `no_comparable_uniprot_contacts`, may still hold author-space contacts, mapped
residues and annotations, and is excluded from the similarity matrix, clustering,
conservation on the comparison representation, contact frequencies (including the
JSON export) and rewiring. Missing correspondence must not be read as biological
variation, and two instances that both lack correspondence must not be read as
biologically identical. `jaccard(set(), set())` remains 1.0 as a mathematical
convention; records with an empty comparison set never reach that stage.

### Independent residue-mapping rule (as documented)

Residue-level UniProt correspondence is decided per residue: a residue is mapped
whenever both a UniProt accession and a UniProt sequence position are present, and
nothing is inferred. Contact-pair comparability is decided per contact: a contact
enters the comparison set only when both residues map. For a partially mapped contact:

```text
A50  → P12345:50 (role 1)        mapping retained in author_to_uniprot
B106 → no UniProt mapping        remains unmapped
A50–B106                         not added to the UniProt-keyed comparison set
```

Annotations on `A50` therefore resolve to `P12345:50`; annotations on `B106` are
reported at the interface with no UniProt key. Microheterogeneity keeps the first-seen
rule, applied per residue.

---

## Homodimer wording

User-level paragraph added to the spec:

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

It is followed by a "Production behaviour: unchanged" paragraph stating that the
evidence supports a note, not a redesign, and that canonicalisation or an
orientation-aware similarity would be a separate decision. The earlier absolute
sentence ("shares no ordered contacts with its twin") appears nowhere in `spec/new/`.

---

## Diagnostic interpretation

The numbers sit in a developer note under the known limitation, introduced as "an
exploratory diagnostic sample, not an unbiased PDB-wide prevalence estimate", with
the sampling method (homodimer complexes with 8 to 60 assemblies from the PDBe-KB
complexes listing) and the forced STING and triosephosphate isomerase examples stated.
Each figure is phrased as a fraction of pairwise comparisons "in this sample":

```text
28 homodimer complexes, 474 interfaces, 6,479 pairwise comparisons
swapped > direct        4.4 % of pairwise comparisons
delta > 0.2             0.46 %
delta > 0.5             0.08 %
partition changed       3 of 28 complexes
cluster count changed   2 of 28 complexes
```

The same-entry observation is recorded there: six of the twelve pairs with a delta of
at least 0.5 were two assemblies of one PDB entry, so opposite orientation can arise
from the reported chain order alone. It is explicitly marked as something to check by
hand, not an algorithmic rule.

The tranche report's interpretation and proposed-note sections were revised to the
same non-absolute wording and now label the content as a diagnostic finding distinct
from production behaviour. The documented conclusion: the homodimer orientation effect
is real but generally modest for clustering in the exploratory sample, and most
consequential for sparse contact fingerprints. It does not support redesigning
orientation handling, canonicalising contacts, or changing similarity or rewiring now.

---

## Provenance visibility

No notebook change was needed. Phase 1 prints the assembly-filter breakdown and
"Retained N interfaces across M PDB entries". Phase 2 then prints either
"All N interfaces have at least one comparable UniProt contact pair" or:

```text
Excluded k of N interfaces from the contact-based comparison:
  <label>  no_comparable_uniprot_contacts  (<a> author-space contacts, <m> residues mapped to UniProt)
Retained n interfaces for comparison.
```

The per-interface summary table below it lists every instance with `n_uniprot_pairs`.
Total, retained, excluded and reason are all visible.

---

## Tests

```text
passed   87
failed   0
skipped  0
xfail    0
runtime  2.5 s
```

The network guard test passes and no test was added or changed.

---

## Files changed in this task

- `spec/new/final_spec.md`: the sections listed above.
- `spec/new/testing_tranche1_5_report.md`: interpretation lead-in and proposed-note
  section revised.

No other file was touched. The working tree also carries the earlier, already-accepted
tranche 1.5 changes to `pdbe_interfaces/representation.py` and `notebook.ipynb`, plus
the untracked `tests/`, `analysis/`, `pytest.ini`, `requirements-dev.txt` and the
other `spec/new/` reports. No scientific algorithm, homodimer behaviour, Jaccard
semantics, clustering, rewiring, identifier validation or API behaviour was modified.
Nothing is committed.

## Document sequence

```text
spec/new/testing_strategy.md                 long-term plan
spec/new/testing_tranche1_report.md          first offline suite, design questions
spec/new/testing_tranche1_consolidation.md   suite streamlined to 83 items
spec/new/testing_tranche1_5_report.md        exclusion rule, partial mappings, homodimer diagnostic
spec/new/testing_tranche1_5_closeout.md      this document
```
