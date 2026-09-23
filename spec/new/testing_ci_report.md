# Continuous integration: final testing report

Completed 2026-09-23. Final task in the testing sequence: a minimal GitHub
Actions workflow that runs the existing offline pytest suite on pushes and pull
requests.

Nothing beyond `checkout → install → pytest` was added. No coverage, linting,
type checking, security or dependency scanning, Docker, test matrix, caching,
scheduled job, badge or deployment step.

---

## Workflow

```text
file            .github/workflows/tests.yml
name            tests
triggers        push to main; pull_request (any branch)
runner          ubuntu-latest
Python          3.11 (single version, no matrix)
install         python -m pip install --upgrade pip
                pip install -r requirements.txt -r requirements-dev.txt
test command    pytest -v
```

The repository had no `.github/` directory, so no existing workflow was
duplicated or replaced.

### Why pushes are restricted to `main`

The repository's primary branch is `main`. Restricting the `push` trigger to it
means work on a feature branch is checked once, by the `pull_request` event,
rather than twice. Every change still reaches CI before it can be merged.

### Why Python 3.11

`README.md` states "Tested on Python 3.11", `requirements.txt` says "Tested with
Python 3.11", and the local development interpreter is 3.11.9. There is no
`pyproject.toml`, `setup.py` or `tox.ini` declaring anything broader. This CI is
a regression guard, not a compatibility matrix, so one version matching the
documented environment is the right scope.

### Dependency installation

Both existing requirements files are used as-is. `requirements.txt` carries the
runtime stack including `jupyter`, `ipykernel`, `matplotlib` and `molviewspec`;
`requirements-dev.txt` carries `pytest`, plus `nbformat` and `nbclient` for the
notebook smoke test. No dependency management was restructured and no package
manager was introduced.

---

## Verification

The workflow was validated three ways before being finalised.

**YAML structure.** Parsed and the job, runner and step sequence inspected. Note
that a YAML 1.1 parser reads the key `on` as boolean `True`; GitHub's own parser
does not, and the file is written in the conventional form.

**Clean-environment install.** A fresh Python 3.11 virtual environment was
created and both requirements files installed, exactly as the workflow does. The
install succeeded with no missing declaration. Dependency resolution picked
newer versions than the local environment in several cases (pandas 3.0.6,
numpy 2.4.6, pytest 9.1.1, matplotlib 3.11.2, nbclient 0.11.0), which is what a
fresh CI runner will also do.

**Full suite in that environment.** All 98 tests passed. The one risk specific
to the notebook smoke test, whether the `python3` kernelspec is discoverable
after a plain `pip install`, was checked explicitly: `ipykernel` registers it in
the package's share directory and `jupyter_client` resolves it, both locally and
in the clean environment.

### Runtime

| Run | Result | Time |
|---|---|---|
| Clean venv, first (cold bytecode cache) | 98 passed | 39 s |
| Clean venv, warm | 98 passed | 6.9 s |
| Clean venv, unit + e2e only | 97 passed | 0.25 s |
| Clean venv, smoke only | 1 passed | 5.8 s |
| Local environment | 98 passed | 7.0 s |

The first clean run is slower because Python compiles the freshly installed
packages; CI will see something in that range on a cold runner. The notebook
smoke test dominates, almost entirely kernel startup. Either figure is
comfortable for a regression guard.

---

## Network behaviour

The pytest suite is fully offline and CI requires no secrets, tokens or
credentials. The workflow references none.

- The STING, Spike/ACE2 and notebook smoke tests are served from the frozen
  responses under `tests/fixtures/recorded/`, through the `RecordedPDBe` router,
  which raises on any URL that was not captured.
- `tests/conftest.py` installs an autouse guard that fails any in-process HTTP
  request, and `tests/unit/test_api.py::test_network_is_blocked_by_default`
  asserts the guard is active.
- `pytest.ini` sets `addopts = -m "not integration"`, so any future live test
  would be excluded from the default run that CI executes.
- The notebook smoke test runs in a separate kernel process, outside the reach of
  the in-process guard; the recorded router installed inside that kernel is what
  keeps it offline, and the test asserts that all 32 requests were served by it.

Package installation obviously needs network access. Test execution does not.

---

## Current result

```text
passed    98
failed     0
skipped    0
xfail      0
runtime    7.0 s locally, 6.9 s warm in a clean CI-like environment
```

---

## GitHub Actions status

**Workflow created but not yet executed on GitHub.** Nothing in this session has
been committed or pushed, by consistent intent throughout the testing work, so
no run has been triggered. The workflow will run on the first push to `main` or
the first pull request after these changes are committed.

---

## Files changed

```text
.github/workflows/tests.yml        new; the CI workflow
README.md                          new Tests section; tests/ added to the project layout
spec/new/testing_ci_report.md      new; this document
```

No production code was modified in this task. `pdbe_interfaces/`,
`notebook.ipynb`, the fixtures and every scientific regression assertion are
untouched; the modifications those files show in `git status` are the
already-accepted tranche-1.5 and tranche-2 changes.

---

## Testing project status

The suite now stands at **98 offline tests** across four layers:

| Layer | Tests | Protects |
|---|---|---|
| `tests/unit/` | 87 | residue correspondence, partner orientation, similarity and clustering, annotation mapping, frequency and conservation arithmetic, identifier resolution, API retry and failure handling |
| `tests/e2e/` STING | 5 | the full scientific state-analysis path on a homodimer: interface states, rewiring, annotations, QC, export |
| `tests/e2e/` Spike/ACE2 | 5 | heterodimer partner correspondence and annotation integrity through partner-order reversal |
| `tests/smoke/` | 1 | notebook orchestration: all 18 code cells execute against recorded data |

Running automatically on pushes to `main` and on pull requests.

> Core testing work is complete for now.

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
spec/new/testing_notebook_smoke_report.md               notebook smoke test close-out
spec/new/testing_ci_report.md                           this document
```
