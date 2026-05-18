# interface_analysis_notebook

Jupyter notebook implementation of the Aggregated Interface Interaction Analysis workflow.

Compares protein–protein interfaces across all deposited structures of a single PDBe-KB complex by interaction-pair similarity, then overlays mutations, modifications, and ligands. Cluster representatives are rendered in Mol*, and per-residue conservation plus residue–residue contact frequencies are exported as JSON for downstream visualisation (e.g. AFDB-style Mol* colouring).

See `spec/` for the full specification (`spec/new/final_spec.md` is the current spec; `spec/old/` holds the earlier draft and the implementation brief).

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter notebook notebook.ipynb
```

Tested on Python 3.11.

## Working example

`PDB-CPX-130306` (KRAS–RAF1 heterodimer). Edit the `Config` instance in the first cell of the notebook to target a different complex. `PDB-CPX-130029` (homodimer) is a useful second case for sanity-checking role symmetry.

Key `Config` fields:
- `complex_id` — PDBe-KB complex identifier.
- `max_resolution` — drop assemblies above this Å threshold (and those with no reported resolution); `None` keeps everything.
- `output_dir` — destination for the Phase 7 JSON export. Defaults to `interface_frequencies/` (created if missing); set to `None` to write to the current directory.

## Project layout

```
pdbe_interfaces/
  api.py             PDBe v2 API calls
  representation.py  InterfaceRecord, interaction-pair sets
  similarity.py      Jaccard, hierarchical clustering
  annotations.py     mutation / modification / ligand workflow
  outputs.py         structure table, conserved sets, cluster report, JSON export
  visualize.py       Mol* cluster-representative builders
notebook.ipynb       phase-by-phase narrative (Phases 1–7)
spec/
  new/               current spec (final_spec.md) + demo notes
  old/               earlier draft + implementation brief
```

## Notebook phases

1. Data retrieval (PDBe interfaces, complex details, annotations)
2. Build `InterfaceRecord`s (author- and UniProt-keyed interaction pairs)
3. Jaccard similarity + hierarchical clustering
4. Annotation overlap (mutations, modifications, ligands)
5. Summarisation & per-cluster interpretation
   - 5a. Key residue pairs (interactive widgets + contact heatmap, axes ordered by ascending UniProt position)
   - 5b. Interface rewiring between interaction states
6. Mol* visualisation of cluster representatives
7. JSON export of per-residue conservation and residue–residue contact frequencies (`{complex_id}.json`)

## Scope (v1)

Dimer complexes only (upstream API constraint). Single complex per run. Descriptive analysis at sample sizes ~10–30 interfaces — no statistical enrichment tests; the cluster-enrichment heuristic uses simple ratio thresholds (`min_count=2`, `min_enrichment=2.0`).
