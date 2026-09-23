# Notebook smoke test: final report

Close-out summary for the two documentation corrections and the notebook smoke
test. Completed 2026-09-23.

The detailed document is `testing_notebook_smoke.md`. This document is the
summary report.

98 tests pass offline. No production code was changed. CI remains the next and
final testing task and was deliberately not started.

---

## Documentation corrections

### Component and species aggregation wording

The statement added during the tranche-2B close-out described the aggregation
rule in terms of species and said every run analyses variation "within one
species". Too broad: a complex may legitimately contain components from
different organisms, and one of the worked examples already does. The paragraph
in `final_spec.md` now reads:

> **Complexes are aggregated by mapped component identity and stoichiometry.**
> Assemblies grouped under a PDB complex identifier share the same mapped
> component composition, represented by UniProt accessions or, where a polymer is
> unmapped, by the corresponding unmapped-entity identifiers. Because orthologous
> proteins carry different component identities, replacing a component with its
> orthologue yields a different composition and therefore a different complex
> identifier, so variation between complexes containing different orthologues is
> never mixed into one aggregation. Every run of this workflow consequently
> compares structural instances of the same component composition, and the
> interface states and contact differences it reports cannot be produced by
> substituting one orthologue for another.
>
> An individual complex may nevertheless contain components from more than one
> organism: the Spike RBD with ACE2 complex pairs human ACE2 (`Q9BYF1`) with
> SARS-CoV-2 Spike (`P0DTC2`). The constraint is on component identity, not on
> the complex being single-organism. Where a result is compared with literature
> on an orthologue of one component, that comparison concerns transferability of
> the interpretation and requires explicit residue correspondence; it is not part
> of the aggregation.

Both STING literature documents were updated the same way: the general rule is
now stated in terms of component identity, with an explicit note that a complex
may span organisms, followed by the stronger case-specific conclusion, that this
complex's composition is a single accession at stoichiometry two, so **all
instances in the STING aggregation contain mouse STING `Q3TBT3` and its
interface states cannot be caused by differences between mouse and human STING
orthologues**.

The phrases "within one species" and "within a single species" no longer appear
as general claims anywhere in `spec/new/`.

### The `7rpv` control wording

`7rpv` assembly 1 against assembly 4 had been described as "same entry, so
orientation is the only difference". Too strong: their normalised Jaccard is
0.55, so the assemblies do differ in the contacts detected in each copy. The
wording is now, in the documents, the test docstring and the fixture metadata:

> Same entry and same biological ACE2-RBD interface as the reversed assembly; a
> strong control for partner-order reversal, while allowing genuine differences
> in the contacts detected in each copy.

The regression itself is unchanged. It already asserted the real values, 0.0
before normalisation and 0.55 and 0.6923 after, rather than an idealised
identity; only the prose overstated the control.

---

## Smoke-test approach

`tests/smoke/test_notebook_smoke.py` loads the real `notebook.ipynb` with
`nbformat` and executes it in a real IPython kernel with `nbclient`, so the
notebook's own cell sources run rather than a Python transcription of them.
Three cells are injected by the test and no notebook cell is edited:

1. **Setup**, first: selects the `Agg` matplotlib backend, sets
   `PDBE_INTERFACES_STATIC=1`, puts the repository and the e2e helpers on
   `sys.path`, and installs `RecordedPDBe` as `requests.Session.request`.
2. **Override**, immediately after the configuration cell: redirects
   `Config.output_dir` to a pytest temporary directory with
   `dataclasses.replace`, so the Phase 7 export never writes into the repository.
3. **Summary**, last: prints a small JSON blob of the values the test checks,
   parsed back out of the executed notebook's outputs.

`allow_errors=False`, so any cell that raises fails the test with that cell's
traceback.

> **The smoke test verifies notebook orchestration. Scientific correctness is
> protected separately by the STING and Spike/ACE2 end-to-end regressions.**

It catches what those tests structurally cannot, because they reproduce the
notebook's call sequence in Python: a cell that no longer matches the module
API, a renamed variable a later cell depends on, a cell-ordering mistake, or a
display call that raises headless.

---

## Cells and stages exercised

All **18** code cells execute.

| Phase | Reached |
|---|---|
| Configuration | yes |
| 1. Data retrieval (resolve, select, mutations, modifications, ligands) | yes |
| 2. Build representations, comparable selection, summary table | yes |
| 3. Similarity, clustering, heatmap, sweep, dendrogram | yes |
| 4. Annotation overlap and preview tables | yes |
| 5. Structure table, conservation, cluster interpretation report | yes |
| 5a. Key residue pairs (explorer, static mode) | yes |
| 5b. Interface rewiring | yes |
| 6. Mol* cluster representatives | yes |
| 7. JSON export and preview | yes |

---

## Exclusions

**None.** Every code cell runs; no cell is skipped, commented out or stubbed.

Two environment settings make the display-dependent cells behave headlessly, and
both are settings the project already documents for headless use rather than
test-only workarounds:

- `matplotlib.use("Agg")` renders figures to memory instead of opening windows.
- `PDBE_INTERFACES_STATIC=1` makes the Phase 5a explorer render the
  whole-dataset view once instead of building ipywidgets. The README records
  that the widget layout hangs `nbconvert --execute`, so this is the documented
  headless mode.

Mol* was expected to need exclusion and did not: `visualize_clusters_grid`
builds a MolViewSpec and emits a display object without fetching anything from
Python, so it runs cleanly in the kernel.

---

## Tests

```text
tests before       97
new smoke tests     1
tests after        98

passed             98
failed              0
skipped             0
xfail               0

full-suite runtime  7.3 s   (8.4 s wall clock)
smoke-test alone    5.0 s   (5.9 s wall clock)
```

Kernel startup dominates the smoke test; the suite ran in about 2.5 s before it,
and remains practical for routine development.

### Network and artefacts

The kernel is a separate process, so the suite's in-process `block_network`
guard does not reach it. Offline behaviour is enforced by the recorded router
installed inside the kernel, which raises on any URL not captured in the frozen
STING fixture. The test asserts that all **32** requests were served by the
router, the same count the STING end-to-end test asserts.

No artefacts remain: the export is redirected to a pytest temporary directory,
the test asserts no `interface_frequencies/` directory appears in the
repository, and `notebook.ipynb` is read but never written back (its git diff is
unchanged from tranche 1.5).

---

## Dependencies

`nbformat>=5.10` and `nbclient>=0.10` added to `requirements-dev.txt`. Both were
already installed as transitive dependencies of `jupyter` in
`requirements.txt`; they are declared because the test imports them directly.

No Jupyter server package was added, and no notebook-testing framework such as
`nbval` or `pytest-notebook`.

---

## Production changes

None. `pdbe_interfaces/` and `notebook.ipynb` were not modified. The notebook
executed on the first attempt, so there is no execution bug to report.

Files changed in this task:

```text
spec/new/final_spec.md                                  aggregation-rule wording
spec/new/testing_tranche2a_sting_literature.md          aggregation-rule wording
spec/new/testing_tranche2a_sting_literature_report.md   aggregation-rule wording
spec/new/testing_tranche2b_spike_ace2_e2e.md            7rpv control wording
spec/new/testing_tranche2b_final_report.md              7rpv control wording, superseded quote
tests/e2e/test_spike_ace2_recorded.py                   docstring only, no assertion changed
tests/fixtures/recorded/spike_ace2_reversal/metadata.json   selection rationale text only
tests/smoke/test_notebook_smoke.py                      new
requirements-dev.txt                                    nbformat, nbclient
spec/new/testing_notebook_smoke.md                      new
spec/new/testing_notebook_smoke_report.md               this document
```

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
spec/new/testing_tranche2b_final_report.md              tranche 2B close-out
spec/new/testing_notebook_smoke.md                      notebook smoke test and wording corrections
spec/new/testing_notebook_smoke_report.md               this document
```

Remaining: CI, the next and final testing task.
