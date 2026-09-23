# Testing tranche 1: offline unit suite, results and observations

Report on the first implementation tranche of the testing strategy in
`testing_strategy.md`, completed 2026-09-23. Scope was deliberately limited to an
offline unit suite protecting the parts of the workflow most likely to produce a
scientifically wrong result while still executing successfully. End-to-end
regression tests, live API tests, notebook smoke tests and CI were left for the
next phase.

All tests run without network access. A plain `pytest` makes zero requests.

---

## Test inventory

| Area | File | Tests |
|---|---|---|
| Interface / residue correspondence | `tests/unit/test_interface_records.py` | 21 |
| Partner reversal and homodimer orientation | `tests/unit/test_partner_orientation.py` | 16 |
| Similarity and clustering | `tests/unit/test_similarity.py` | 14 |
| Annotation mapping (mutations, modifications, ligands) | `tests/unit/test_annotations.py` | 20 |
| Frequency, conservation, rewiring labels, cluster in/out counts | `tests/unit/test_frequencies.py` | 10 |
| Identifier resolution | `tests/unit/test_identifiers.py` | 24 |
| Mocked API behaviour | `tests/unit/test_api.py` | 23 |
| **Total** | | **128** |

### Infrastructure

- `tests/conftest.py`: an autouse fixture replaces `requests.Session.request` so
  any accidental network call fails loudly, plus two builders (`make_contact`,
  `make_interface`) for minimal `interface_interactions` payloads whose expected
  values can be worked out by inspection.
- `pytest.ini`: sets `testpaths`, `pythonpath = .` (the package is not installed,
  so a bare `pytest` could not import it otherwise), registers an `integration`
  marker and excludes it by default with `-m "not integration"`.
- `requirements-dev.txt`: lists `pytest` only. No mocking library was needed;
  `monkeypatch` covers every case.

### Running

```bash
pip install -r requirements-dev.txt
pytest            # offline suite, default
pytest -v         # verbose
```

---

## Results

```text
passed   128
failed   0
skipped  0
xfail    0
```

The single warning in the output is a pydantic deprecation raised by `molviewspec`
on import, not by this code.

Four tests failed on the first run. All four were mistakes in the test
expectations, not in production code, and were corrected in the tests:

- two tests removed UniProt positions but left accessions in place, and the record
  builder infers partner accessions from accessions alone (now recorded as an
  observation below);
- one expected set omitted a contact pair that legitimately occurs in two of three
  interfaces;
- one compared a numpy bool with `is False`.

---

## Scientifically important observations

### Homodimer orientation

Ordered homodimer contacts do affect cross-instance similarity. With the same
accession on both sides, `check_partner_consistency` cannot act, and tuples stay
ordered. The tests record three outcomes:

- Two identical interfaces with chains listed in opposite order score Jaccard 0.0
  and separate at the default cut of 0.6.
- A partially symmetric interface (two asymmetric contacts plus one self-symmetric
  contact) scores 0.2 after a flip, since only the self-symmetric contact survives.
- An interface where PISA reports both directions of every contact is
  orientation-invariant.

The spec (`final_spec.md`, "Edge cases: Homodimers") states that ordered tuples
are intended. What the code cannot guarantee is that PISA orders homodimer chains
consistently across entries. If it does not for some complex, equivalent
interfaces will form separate states and the rewiring table will report
"contacts" that are the same contact seen from the other chain. The `1spq`
single-cluster result is evidence for that complex only, not for the general case.

### Partner reversal

Works as designed for heterodimers. After the check, the reversed record has the
canonical accession in role 1, identical UniProt tuples to its equivalents,
`author_to_uniprot` and `residue_identity` on the new roles, author pairs swapped,
and Jaccard 1.0 against its equivalents. A mutation on the reversed entry's
X-chain residue lands on role 1.

Two further behaviours are recorded:

- A 1:1 tie makes the first record's ordering canonical (`Counter.most_common`
  preserves first-seen order among equal counts).
- A record whose accession pair matches neither ordering is left as-is with a
  warning.

### Missing UniProt mappings

- When one side of a contact lacks a mapping, neither residue of that contact
  enters `author_to_uniprot`, including the side that did carry a mapping. An
  annotation on that residue is then reported with the `(None, None, None)`
  sentinel unless the residue also appears in a fully mapped contact.
- Partner accessions are inferred from accessions alone. A contact with
  accessions but no positions still sets the record's partner accessions while
  contributing nothing to the comparison space.
- Two interfaces that lost every contact to missing mappings have Jaccard 1.0 by
  the documented empty-vs-empty rule, so they would cluster together as identical.
- Contacts missing author chain or residue number are skipped silently and are not
  counted in `n_residues_dropped_no_uniprot`.

### Insertion codes

Normalisation (`""`, `" "`, padded codes) and the fallback `ins_code_1` field
behave correctly. The annotation join is code-sensitive in both directions: A12
does not match A12A and A12A does not match A12. A whitespace insertion code on an
annotation normalises to None.

### Duplicate contacts

Repeated atom-level observations collapse to one residue pair in both spaces.
Different bond types on one residue pair are distinct typed tuples and one untyped
tuple. Duplicate atom-level ligand contacts collapse to one contact per
interaction label. A missing `bond_type` becomes `"unknown"`.

### Annotation role mapping

Correct accession, position and role for both partners. The same residue number on
another chain does not match. Off-interface annotations are dropped. `Conflict`
mutations are excluded by default and included on request. An entry with two
assemblies exposing the same chain produces one row per interface. Annotations for
an entry not in the record set are ignored. Ligand contacts to residues outside the
interface are dropped, and only ligands surviving the blocklist and glycan filter
are fetched.

### Frequency denominators

Residue and contact frequencies use the interface count as denominator and count a
residue once per interface regardless of how many contacts it makes. Thresholds are
inclusive at exact equality (2 of 3 at threshold 2/3 qualifies; 0.67 excludes it)
in `conserved_residues`, `conserved_interaction_pairs` and the report's core
contacts. Conserved pairs are typed. The untyped JSON export still counts a residue
pair once per interface when it has two bond types. Rewiring labels (`shared core`,
`higher in A`, `higher in B`) and the report's in-cluster / rest-of-dataset counts
were verified on a two-cluster example.

### API layer

Retry on 500/503, timeout and connection errors up to `MAX_ATTEMPTS`; no retry on
400 or 404. Tolerated endpoints (`bound_molecules`, `bound_ligand_interactions`)
return empty after exhaustion while critical endpoints raise. Batched POST chunks at
`POST_BATCH_SIZE` and merges, and a 404 on one chunk does not discard the others.
An empty id list makes no request. Parallel fetchers preserve input order even when
the first item responds last.

---

## Design questions

Recorded as observed behaviour in the tests, not asserted as required. Each needs
a maintainer decision before any production change.

1. **Homodimer chain ordering.** Should equivalent homodimer interfaces with
   flipped PISA chain order be treated as the same state? Options range from
   leaving it (current, spec-endorsed) to canonicalising each contact or each
   record. Evidence on how PISA actually orders homodimer chains is needed first.
2. **Sentinel for mapped residues in dropped contacts.** A residue with a UniProt
   mapping loses it if its only contact has an unmapped partner. Should the mapped
   side still be entered into `author_to_uniprot`?
3. **Empty-vs-empty Jaccard of 1.0.** Should interfaces with no mapped contacts be
   excluded from clustering rather than treated as identical? The spec's
   failure-handling section says such an interface should be "excluded from the
   similarity matrix"; the code does not do that.
4. **Unknown complex id.** A 404 from `complex/details?id_type=pdb_complex_id`
   raises `requests.HTTPError`; the spec says `ValueError`. Which is intended?
5. **Malformed entry ids.** Strings like `not-an-id` or `CPX-2128` reach the API
   as entry ids. No test was written, per instruction.
6. **Non-UniProt participant.** The dimer check tests `total_chains == 2` only; a
   two-chain complex with a non-UniProt participant passes. No test was written.
7. **Contacts missing author fields.** Silently skipped and not counted in
   `n_residues_dropped_no_uniprot`. Should they be counted or logged?
8. **Partner-ordering tie.** With an even split the first record wins. Acceptable,
   or should it raise like the ambiguous-complex case?
9. **`max_entries=0`.** Yields an empty selection with no error (from the earlier
   read-through; not tested here).

---

## Production changes

None. The working tree shows only new files:

```text
pytest.ini
requirements-dev.txt
spec/new/testing_strategy.md
spec/new/testing_tranche1_report.md
tests/
```

The `pythonpath = .` line in `pytest.ini` is the only testing-enablement
configuration and does not affect the notebook or the modules. Nothing is
committed.

---

## Not implemented in this tranche

Left for the next phase, per instruction: recorded STING end-to-end fixture, live
STING test, Spike/ACE2 and TIM end-to-end tests, live PDBe contract tests,
notebook execution smoke test, CI, coverage targets, `xfail` tests for undecided
behaviour, and any change to homodimer semantics.
