# Notebook smoke test, and two documentation corrections

Completed 2026-09-23. One shallow offline smoke test that executes
`notebook.ipynb` itself, plus two wording corrections identified during review of
the tranche-2A and tranche-2B material.

> **The smoke test verifies notebook orchestration. Scientific correctness is
> protected separately by the STING and Spike/ACE2 end-to-end regressions.**

CI is deliberately not part of this task.

---

## Part A — Component and species aggregation wording

The statement added to `final_spec.md` during the tranche-2B close-out described
the aggregation rule in terms of species, and said that every run analyses
variation "within one species". That is too broad: a biological complex can
legitimately contain components from different organisms, and the immediate
counter-example is already one of the worked examples, human ACE2 (`Q9BYF1`)
bound to SARS-CoV-2 Spike (`P0DTC2`).

The rule is about mapped component identity, not organism count. The paragraph in
`final_spec.md` now reads:

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

Both STING literature documents were updated the same way: they now state the
general rule in terms of component identity, note explicitly that a complex may
span organisms, and then keep the stronger case-specific conclusion, that this
complex's composition is a single accession at stoichiometry two, so **all
instances in the STING aggregation contain mouse STING `Q3TBT3` and its interface
states cannot be caused by differences between mouse and human STING
orthologues**.

The phrases "within one species" and "within a single species" no longer appear
as general claims anywhere in `spec/new/`.

## Part B — The `7rpv` control wording

The tranche-2B documentation described `7rpv` assembly 1 against assembly 4 as
"same entry, so orientation is the only difference". That is too strong: their
normalised Jaccard is 0.55, so the two assemblies do differ in the contacts
detected in each copy.

The wording is now, in the documents, the test docstring and the fixture
metadata:

> Same entry and same biological ACE2-RBD interface as the reversed assembly; a
> strong control for partner-order reversal, while allowing genuine differences
> in the contacts detected in each copy.

The regression itself is unchanged. It already asserted the real similarity
values (0.0 before normalisation, 0.55 and 0.6923 after) rather than an
idealised identity; only the prose overstated the control.

---

## What the smoke test protects

It answers one question: can the notebook run its analysis cells to completion
and reach its principal outputs? It catches the class of failure the end-to-end
tests structurally cannot, because those reproduce the notebook's call sequence
in Python rather than running the notebook:

- a notebook cell that no longer matches the module API after a refactor;
- a renamed variable that a later cell depends on;
- a cell-ordering mistake;
- a display or plotting call that raises when run headless.

It is deliberately shallow. It does not re-assert STING cluster membership,
contacts, ligands or export contents; those belong to
`tests/e2e/test_sting_recorded.py`, and duplicating them here would create a
second regression to maintain.

## How the notebook is executed

`tests/smoke/test_notebook_smoke.py` loads the real `notebook.ipynb` with
`nbformat` and executes it in a real IPython kernel with `nbclient`, so the
notebook's own cell sources run. Three cells are injected by the test, and no
notebook cell is edited:

1. **Setup**, before everything: selects the `Agg` matplotlib backend, sets
   `PDBE_INTERFACES_STATIC=1` (the headless mode for the Phase 5a explorer that
   the README documents), puts the repository and the e2e test helpers on
   `sys.path`, and installs `RecordedPDBe` as `requests.Session.request`.
2. **Override**, immediately after the configuration cell: redirects
   `Config.output_dir` to a pytest temporary directory with
   `dataclasses.replace`, so the Phase 7 export never writes into the
   repository.
3. **Summary**, at the end: prints a small JSON blob of the values the test
   checks, which the test parses out of the executed notebook's outputs.

`allow_errors=False`, so any cell that raises fails the test with that cell's
traceback. The kernel's working directory is a temporary directory as well.

## Notebook stages exercised

All **18** code cells execute. In notebook phase terms:

| Phase | Cells | Reached |
|---|---|---|
| Configuration | 1 | yes |
| 1. Data retrieval (resolve, select, mutations, modifications, ligands) | 4 | yes |
| 2. Build representations, comparable selection, summary table | 1 | yes |
| 3. Similarity, clustering, heatmap, sweep, dendrogram | 2 | yes |
| 4. Annotation overlap and preview tables | 3 | yes |
| 5. Structure table, conservation, cluster interpretation report | 3 | yes |
| 5a. Key residue pairs (explorer, static mode) | 1 | yes |
| 5b. Interface rewiring | 1 | yes |
| 6. Mol* cluster representatives | 1 | yes |
| 7. JSON export and preview | 1 | yes |

## Exclusions

**None.** Every code cell runs, including the Phase 3 matplotlib figures, the
Phase 5a explorer and the Phase 6 Mol* rendering. No cell is skipped, commented
out or stubbed.

Two environment settings make the display-dependent cells behave headlessly,
and both are settings the project already documents for headless use rather than
test-only workarounds:

- `matplotlib.use("Agg")` renders figures to memory instead of opening windows;
  `plt.show()` becomes a no-op.
- `PDBE_INTERFACES_STATIC=1` makes `explorer.residue_pair_explorer` render the
  whole-dataset view once instead of building ipywidgets. The README records
  that the widget layout hangs `nbconvert --execute`, so this is the documented
  way to run the notebook headlessly.

Mol* was expected to need exclusion and did not: `visualize.visualize_clusters_grid`
builds a MolViewSpec and emits a display object without fetching anything from
Python, so it runs cleanly in the kernel. Only a browser would fetch the
structure files it references.

## Recorded fixture and network behaviour

The frozen STING fixture `tests/fixtures/recorded/PDB-CPX-172174/` is reused
unchanged, through the same `RecordedPDBe` router the end-to-end tests use. No
new fixture system was introduced.

The kernel is a separate process, so the suite's in-process `block_network`
guard does not reach it. What keeps this test offline is the recorded router
installed in the kernel, which raises on any URL that was not captured. The test
asserts that all **32** requests were served by the router, which is the same
count the STING end-to-end test asserts, so a live call for a recorded endpoint
would change the count and a call for an unrecorded one would raise.

The export is redirected to a pytest temporary directory, and the test asserts
both that the file landed there and that no `interface_frequencies/` directory
was created in the repository. `notebook.ipynb` is read but never written back.

## Assertions

Deliberately few:

- all 18 notebook code cells executed (implicit: `nbclient` raises otherwise);
- 32 requests served by the recorded router;
- the complex resolved to `PDB-CPX-172174`;
- comparable records exist, and comparable plus excluded equals all records;
- the similarity matrix is square with one row per comparable record;
- one cluster label per comparable record;
- one structure-table row per comparable record;
- the cluster report and rewiring table are non-empty, with at most one report
  row per record;
- the export produced residue and contact frequencies and wrote its file to the
  temporary directory.

## Runtime

```text
smoke test alone      5.4 s   (6.4 s wall clock, including kernel startup)
full suite            6.9 s   (8.0 s wall clock)
```

Kernel startup dominates. The suite remains practical for routine development;
before this test it ran in about 2.5 s.

## Dependencies

`nbformat>=5.10` and `nbclient>=0.10` added to `requirements-dev.txt`. Both were
already installed as transitive dependencies of `jupyter` in
`requirements.txt`; they are declared because the test imports them directly. No
Jupyter server package, and no notebook-testing framework such as `nbval` or
`pytest-notebook`, was added.

## Production changes

None. `pdbe_interfaces/` and `notebook.ipynb` were not modified. The notebook
executed on the first attempt with no bug to report.
