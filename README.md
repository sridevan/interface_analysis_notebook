# Aggregated Interface Interaction Analysis in Dimers

Jupyter notebook implementation of the Aggregated Interface Interaction Analysis workflow,
for **dimeric complexes only**.

Compares the protein–protein interface across all deposited structures of a single dimeric
complex by interaction-pair similarity, then overlays mutations, modifications, and ligands.
Cluster representatives are rendered in Mol*, and per-residue conservation plus
residue–residue contact frequencies are exported as JSON for downstream visualisation
(e.g. AFDB-style Mol* colouring).

The complex is named by a PDB entry id or a PDBe-KB complex id, and both components must map
to UniProt. Complexes with more than two components are rejected, and assembly instances
carrying an additional bound macromolecule are excluded: chain correspondence cannot be
determined beyond two components, so those interfaces are not comparable across structures.

The notebook holds the narrative, the configuration and the guidance for reading each
output; all implementation lives in `pdbe_interfaces/`.

See `spec/` for the full specification (`spec/new/final_spec.md` is the current spec; `spec/old/` holds the earlier draft and the implementation brief).

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter notebook notebook.ipynb
```

Tested on Python 3.11.

## Interpreting the output

Scope and limitations are documented in the notebook, as guidance at the phase each applies
to. In short: the contact and annotation columns are counts rather than enrichment
statistics, since clusters are built from the same contacts those columns list; deposited
structures are not independent observations, so frequencies describe the deposition set
analysed; and the interaction types available at residue-pair level should be checked per
complex before interpreting any comparison.

## Working example

`11gl` (STING, complex `PDB-CPX-172174`, Complex Portal `CPX-2128`): a homodimer of 14
interfaces from 12 entries, all X-ray at 1.29 to 2.75 Å, separating into four states of 9, 3,
1 and 1, with ligands in three of them. Small enough to read end to end while still
exercising the clustering, rewiring and QC output. Edit the `Config` instance in the first
cell of the notebook to target a different complex.

Other useful cases: `1spq` (triosephosphate isomerase, `PDB-CPX-130029`) is a homodimer whose
interface is invariant, forming a single cluster at every cut; `6m0j` (Spike RBD with ACE2,
`PDB-CPX-140195`) is a 130-interface heterodimer; `PDB-CPX-130306` (KRAS with RAF1) is a
smaller heterodimer.

## Configuration

Every setting is written out in the notebook's configuration cell and documented field by
field in the `pdbe_interfaces.config.Config` docstring. The main ones:
- `identifier`: a PDB entry ID (`"11gl"`, resolved via `complex/details/{pdb_id}?id_type=pdb_id`) or a PDBe-KB complex ID (`"PDB-CPX-172174"`). Either way the analysis covers every deposited structure of the resolved complex. Dimers only.
- `max_entries`: restrict the analysis to the first N PDB entries alphabetically, for a quick check; `None` (default) uses the full dataset.
- `max_workers`: threads for the Phase 1 per-entry and per-ligand retrieval (default 8). The ligand endpoint is called once per ligand instance, so this determines retrieval time.
- `max_resolution`: exclude assemblies above this Å threshold and those with no reported resolution; `None` keeps everything.
- `mutation_type_filter`: mutation types kept from the annotation API; the default keeps engineered mutations only, add `"Conflict"` for natural variants.
- `ligand_blocklist`: components excluded from the ligand analysis (ions, buffers, cryoprotectants).
- `cluster_distance_cut`: dendrogram cut height in `1 - Jaccard` units, which sets how many interface interaction states are reported.
- `conservation_threshold`: fraction of interfaces at which a residue or contact counts as conserved.
- `output_dir`: destination for the Phase 7 JSON export. Defaults to `interface_frequencies/`; `None` writes to the working directory.

Retrieval is the slow step and runs concurrently: on the 114-entry Spike RBD with ACE2
complex the ligand fetches take 4.5 s with the default 8 workers, against 32 s sequentially.

## Project layout

```
pdbe_interfaces/
  api.py             PDBe v2 API calls
  config.py          Config dataclass and logging setup
  plots.py           similarity heatmap, dendrogram, pair-frequency heatmap
  explorer.py        interactive residue-pair explorer (ipywidgets)
  representation.py  InterfaceRecord, interaction-pair sets
  similarity.py      Jaccard, hierarchical clustering
  annotations.py     mutation / modification / ligand workflow
  outputs.py         structure table, conserved sets, cluster report, JSON export
  visualize.py       Mol* cluster-representative builders
notebook.ipynb       phase narrative, configuration and guidance (Phases 1–7)
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

Dimer complexes only, with both components mapped to UniProt: chain correspondence cannot
be determined for instances with more than two components, and the `interface_interactions`
endpoint returns dimers only. Single complex per run. The analysis is
descriptive throughout: contacts and annotations are reported as counts inside and outside
each cluster, with no enrichment statistic or significance threshold.

## Licence

Copyright 2026 EMBL - European Bioinformatics Institute. Author: Sri Devan Appasamy.

Licensed under the Apache License, Version 2.0. See [`LICENSE`](LICENSE) for the full text.
You may not use these files except in compliance with the License. Unless required by
applicable law or agreed to in writing, software distributed under the License is
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
express or implied.
