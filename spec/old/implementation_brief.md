# Implementation Brief

Companion to `aggregated_interface_analysis_spec.md`. The spec describes the system; this brief describes the build. Read the spec first.

## Tech stack

- Python 3.11+
- `requests` for HTTP (synchronous; v1 call volume is small enough not to need async)
- `pandas` for tables
- `numpy` for the similarity matrix
- `scipy.cluster.hierarchy` for clustering and dendrogram
- `matplotlib` + `seaborn` for the heatmap and dendrogram
- Standard library: `dataclasses`, `typing`, `logging`, `json`, `functools`

No additional dependencies unless a specific need arises during the build.

## Project structure

```
project_root/
├── notebook.ipynb              # the deliverable; five-phase narrative
├── pdbe_interfaces/            # helper package, importable from the notebook
│   ├── __init__.py
│   ├── api.py                  # all PDBe API calls; returns parsed JSON dicts
│   ├── representation.py       # builds residue keys, interaction-pair sets
│   ├── similarity.py           # Jaccard, distance matrix, clustering
│   ├── annotations.py          # mutations, modifications, ligand workflow
│   └── outputs.py              # structure table, conserved sets, cluster report
└── README.md
```

The notebook orchestrates and visualises. Logic lives in the package so it is testable and reusable. Each notebook cell is short — calls into the package, displays a result.

## Function boundaries

Suggested top-level functions (signatures are guidance, not contracts):

```python
# api.py
def fetch_complex_details(complex_id: str) -> dict
def fetch_interface_interactions(complex_id: str) -> dict
def fetch_mutations(pdb_ids: list[str]) -> dict
def fetch_modifications(pdb_ids: list[str]) -> dict
def fetch_bound_molecules(pdb_id: str) -> list[dict]
def fetch_ligand_interactions(pdb_id: str, chain_id: str, author_residue_number: int) -> list[dict]

# representation.py
def build_interface_records(interface_response: dict) -> list[InterfaceRecord]
def build_interaction_pair_sets(record: InterfaceRecord) -> tuple[set, set]   # author-keyed, uniprot-keyed
def check_partner_consistency(records: list[InterfaceRecord]) -> list[InterfaceRecord]   # warns and reverses if needed

# similarity.py
def jaccard_similarity_matrix(sets: list[set]) -> np.ndarray
def cluster_interfaces(distance_matrix: np.ndarray, method: str = "average") -> ClusterResult

# annotations.py
def filter_bound_molecules(molecules: list[dict], blocklist: set[str]) -> list[dict]
def overlap_annotations(records: list[InterfaceRecord], mutations: dict, modifications: dict, ligand_contacts: dict) -> AnnotationOverlap

# outputs.py
def build_structure_table(...) -> pd.DataFrame
def conserved_residues(records: list[InterfaceRecord], threshold: float = 1.0) -> set
def conserved_interaction_pairs(records: list[InterfaceRecord], threshold: float = 1.0) -> set
def cluster_interpretation_report(...) -> pd.DataFrame
```

`InterfaceRecord` is a dataclass holding the parsed per-interface data — `pdb_id`, `assembly_id`, `interface_id`, the `interface_info` PISA metrics, the two interaction-pair sets, the residue→UniProt mapping, the partner accessions and roles. Keep it dumb — no methods beyond construction.

## Configuration

A single `Config` dataclass at the top of the notebook, exposing every parameter the spec marks as user-overridable:

```python
@dataclass
class Config:
    complex_id: str = "PDB-CPX-130306"
    mutation_type_filter: tuple[str, ...] = ("Engineered mutation",)
    ligand_blocklist: frozenset[str] = frozenset({
        "SO4", "GOL", "HOH", "EDO", "PEG", "MPD",
        "CL", "NA", "MG", "ZN", "CA", "K", "ACT",
    })
    drop_carbohydrate_polymers: bool = True
    conservation_threshold: float = 0.8           # fraction of interfaces; tolerates 1-2 drop-outs at N≈10-15. Set to 1.0 for strict
    cluster_distance_cut: float = 0.5             # 1 - Jaccard
    log_level: str = "INFO"
```

Defaults match the spec. Users edit the dataclass instance, not the function arguments.

## Error handling

- API failures: bail on first error with a clear message. No retries, no partial-data path. v1 is exploratory; surfacing failures fast is more useful than masking them.
- Validation failure (input is not a dimer, complex ID does not resolve): raise `ValueError` with a one-line explanation.
- Logged warnings (per spec): residues dropped for missing UniProt mapping, partner-order inconsistency, ligand filter applied. Use `logging` at WARNING level so they appear in the notebook output without halting.
- Empty results (e.g. no mutations returned for a complex): not an error. Continue with empty annotation lists.

## Visualisation defaults

- Similarity matrix heatmap: `seaborn.heatmap`, `viridis` colormap, square cells, annotations off (Jaccard values readable from the colour scale at this matrix size). Row/column labels: `pdb_id` plus a short annotation suffix where present (e.g. `6vjj [G12C, sotorasib]`).
- Dendrogram: `scipy.cluster.hierarchy.dendrogram`, default orientation, no automatic colouring of clusters in v1 (the cluster cut is shown as a horizontal line at `cluster_distance_cut`).
- Structure table: rendered with `df.style` for the notebook display, no exotic formatting.
- All plots: white background, default fonts, no titles inside the figure (titles in the surrounding markdown).

## Acceptance checks

When run on the working example (`PDB-CPX-130306`), the notebook should produce:

- A structure table with at least 10 rows (~12 PDB entries × assemblies).
- A similarity matrix that is symmetric with 1.0 on the diagonal.
- A dendrogram that visually separates at least two clusters at the default cut.
- A non-empty mutation overlap (KRAS structures will have G12-something mutations at the interface for at least some entries).
- A non-empty ligand overlap (KRAS-G12C inhibitor structures will have sotorasib or similar bound, at the GTP-binding pocket; whether the inhibitor reaches the RAF1 interface depends on the entry).
- No unhandled exceptions.

If the workflow runs on `PDB-CPX-130306` end-to-end and produces these outputs, v1 is functionally complete. Quality of the analysis (i.e. whether the clusters are interpretable) is a separate question and is the user's call.

## Out of scope for the build

These are explicitly not implemented in v1, even if the spec mentions them:
- API response caching to disk.
- Retry/backoff on API failures.
- Async / parallel API calls.
- Tests beyond the acceptance checks above.
- CI, packaging, or distribution.

Add only if a concrete need arises.
