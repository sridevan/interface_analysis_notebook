"""Smoke test: the notebook's analysis path executes end to end on recorded data.

This protects **orchestration**, not science. It answers one question: can
`notebook.ipynb` run its analysis cells to completion and reach its principal
outputs? Scientific correctness is protected separately and in detail by
`tests/e2e/test_sting_recorded.py` (STING) and
`tests/e2e/test_spike_ace2_recorded.py` (Spike RBD with ACE2), so the assertions
here are deliberately shallow and must not be extended into another regression.

What it catches that the end-to-end tests cannot: a notebook cell that no longer
matches the module API, a renamed variable a later cell depends on, a cell
ordering mistake, or a display call that raises headless.

How it runs. The real notebook is loaded with `nbformat` and executed in a real
IPython kernel with `nbclient`, so the notebook's own cell sources run, not a
Python transcription of them. Three cells are injected by the test:

1. a setup cell before everything, which selects the non-interactive matplotlib
   backend, sets `PDBE_INTERFACES_STATIC=1` (the headless mode the README
   documents for the Phase 5a explorer), puts the repository and the e2e test
   helpers on `sys.path`, and installs `RecordedPDBe` as
   `requests.Session.request` so every call is served from the frozen STING
   fixture;
2. an override cell immediately after the configuration cell, which redirects
   `Config.output_dir` to a pytest temporary directory so the Phase 7 export
   never writes into the repository;
3. a summary cell at the end, which prints the few values this test checks.

No notebook cell is skipped or neutralised. All eighteen code cells execute,
including the Phase 3 plots, the Phase 5a explorer and the Phase 6 Mol*
rendering. The kernel runs in its own process, so the suite's in-process network
guard does not reach it; the recorded router is what keeps this test offline,
and it raises on any URL that was not captured.
"""

from __future__ import annotations

import json
import pathlib
import tempfile

import nbformat
import pytest
from nbclient import NotebookClient

REPO = pathlib.Path(__file__).resolve().parents[2]
NOTEBOOK = REPO / "notebook.ipynb"
FIXTURE = "PDB-CPX-172174"
KERNEL_TIMEOUT_S = 300

SETUP_CELL = """
import os, sys, logging, warnings
warnings.filterwarnings("ignore")
os.environ["PDBE_INTERFACES_STATIC"] = "1"      # headless Phase 5a, per the README
import matplotlib; matplotlib.use("Agg")        # no interactive figure windows
sys.path.insert(0, {repo!r})
sys.path.insert(0, {helpers!r})
import requests
from conftest import RecordedPDBe
_RECORDED = RecordedPDBe({fixture!r})
requests.Session.request = _RECORDED            # offline: raises on any unrecorded URL
logging.disable(logging.WARNING)
"""

OVERRIDE_CELL = """
import dataclasses
cfg = dataclasses.replace(cfg, output_dir={output_dir!r})
"""

SUMMARY_CELL = """
import json as _json
print("SMOKE_SUMMARY=" + _json.dumps({
    "complex_id": complex_id,
    "n_all_records": len(all_records),
    "n_comparable": len(records),
    "n_excluded": len(comparable.excluded),
    "sim_rows": int(sim_typed.shape[0]),
    "sim_cols": int(sim_typed.shape[1]),
    "n_cluster_labels": int(len(cluster_result.flat_assignment)),
    "structure_table_rows": int(len(outputs.build_structure_table(
        records, cluster_result, overlap, partner_map))),
    "report_rows": int(len(report)),
    "rewiring_rows": int(len(rewiring)),
    "n_residue_freqs": len(export_data["residue_frequencies"]),
    "n_contact_freqs": len(export_data["contact_frequencies"]),
    "requests_served": len(_RECORDED.urls),
}))
"""


@pytest.fixture(scope="module")
def executed_notebook(tmp_path_factory):
    """Execute the real notebook once in a kernel; return (summary, output_dir)."""
    output_dir = tmp_path_factory.mktemp("notebook_export")
    nb = nbformat.read(NOTEBOOK, as_version=4)

    cells = [nbformat.v4.new_code_cell(SETUP_CELL.format(
        repo=str(REPO), helpers=str(REPO / "tests" / "e2e"), fixture=FIXTURE))]
    n_notebook_cells = 0
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        cells.append(cell)
        n_notebook_cells += 1
        if "cfg = Config(" in cell.source:
            cells.append(nbformat.v4.new_code_cell(
                OVERRIDE_CELL.format(output_dir=str(output_dir))))
    cells.append(nbformat.v4.new_code_cell(SUMMARY_CELL))
    nb.cells = cells

    # allow_errors=False: any cell raising fails the test with that cell's traceback.
    with tempfile.TemporaryDirectory() as run_dir:
        NotebookClient(
            nb, timeout=KERNEL_TIMEOUT_S, kernel_name="python3",
            allow_errors=False, resources={"metadata": {"path": run_dir}},
        ).execute()

    summary = None
    for cell in nb.cells:
        for out in cell.get("outputs", []):
            text = out.get("text") or ""
            if "SMOKE_SUMMARY=" in text:
                summary = json.loads(text.split("SMOKE_SUMMARY=", 1)[1].strip().splitlines()[0])
    assert summary is not None, "notebook ran but produced no summary output"
    summary["n_notebook_code_cells"] = n_notebook_cells
    return summary, output_dir


def test_notebook_analysis_path_executes_on_recorded_data(executed_notebook):
    summary, output_dir = executed_notebook

    # Every code cell of the notebook ran; nbclient would have raised otherwise.
    assert summary["n_notebook_code_cells"] == 18

    # Offline: every request was served by the recorded router.
    assert summary["requests_served"] == 32

    # The principal variables later cells depend on exist and are self-consistent.
    assert summary["complex_id"] == FIXTURE
    assert summary["n_comparable"] > 0
    assert summary["n_all_records"] == summary["n_comparable"] + summary["n_excluded"]
    assert summary["sim_rows"] == summary["sim_cols"] == summary["n_comparable"]
    assert summary["n_cluster_labels"] == summary["n_comparable"]
    assert summary["structure_table_rows"] == summary["n_comparable"]

    # The interpretation and rewiring stages produced output.
    assert 0 < summary["report_rows"] <= summary["n_comparable"]
    assert summary["rewiring_rows"] > 0

    # Phase 7 export completed, into the temporary directory rather than the repo.
    assert summary["n_residue_freqs"] > 0 and summary["n_contact_freqs"] > 0
    assert (output_dir / f"{FIXTURE}.json").is_file()
    assert not (REPO / "interface_frequencies").exists()
