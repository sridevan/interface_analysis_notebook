# Aggregated Interface Interaction Analysis Notebook — Specification

## Overview

The Aggregated Interface Interaction Analysis Notebook is a Jupyter notebook that compares protein–protein interfaces across all available structural instances of the same complex. Given a single PDBe-KB complex identifier, it retrieves every deposited interface for that complex from the PDBe v2 APIs, computes interaction-pair similarity between them, clusters the interfaces, and overlays annotations (engineered mutations, post-translational modifications, bound ligands) so that structural variation can be related to the experimental conditions under which it was observed.

**Core problem.** Interface analyses are usually reported per-structure, even when a complex has many deposited PDB entries that differ in mutations, modifications, ligands, crystal forms, or experimental conditions. Treating each structure in isolation discards the comparative signal that exists across the dataset.

**Intended users.** Structural biologists analysing a known complex; medicinal chemists comparing apo and ligand-bound states; bioinformaticians curating reference datasets of interface variability for downstream pipelines.

**What the notebook does.** Fetches all interfaces for one complex from PDBe, builds a UniProt-keyed interaction-pair representation that is comparable across structures, computes a Jaccard similarity matrix and hierarchical clustering, and joins per-residue annotations onto the interface so each cluster can be inspected for an annotation correlate (e.g. "all G12C structures cluster here", "this cluster has sotorasib bound"). Each cluster is treated as a candidate **interface interaction state** — a recurring pattern of residue–residue contacts that distinguishes a subset of structures from the rest.

**Primary outcome.** A structure table, a similarity heatmap, a dendrogram, conserved-residue and conserved-interaction-pair sets, annotation overlap views, a cluster interpretation report (with per-state statistics, core contacts, enriched contacts and annotations, and QC warnings), a head-to-head **interface rewiring** comparison between any two interaction states, and a lightweight **interface frequency JSON export** (per-residue conservation and per-contact frequencies for downstream visualisation) — each surfaced as an inspectable artefact in the notebook.

**Broader context.** The interaction-pair representation is structure-source agnostic. Although v1 is restricted to PDBe-derived dimer interfaces, the same downstream machinery is intended to extend to non-dimer interfaces (when upstream APIs support them) and to predicted complexes (AlphaFold-Multimer, Boltz, Chai) by replacing only the data-retrieval layer. The notebook is an analysis tool — not a docking, prediction, or sequence-conservation tool.

The working example throughout development is **PDB-CPX-130306 (KRAS–RAF1 heterodimer)**: a clean heterodimer with multiple PDB structures, mutations, and ligands.

---

## Domain / Business / Scientific Context

Implementing this notebook correctly requires understanding several domain concepts that govern how PDBe data is shaped and what interface comparisons mean.

**Protein complexes and interfaces.** A protein complex is an assembly of two or more chains that contact each other in a defined geometry. The "interface" between two chains is the set of residues in physical contact, where contact is operationalised by an upstream pipeline (PISA in this workflow) as atoms within a specified distance, with a contact type (hydrogen bond, salt bridge, covalent bond, etc.). PISA reports contacts at atom level; this workflow aggregates atom-level contacts to residue-pair-and-bond-type elements (see §6 of the source spec; reproduced in Output Parsing).

**Why a single complex has many structures.** PDBe-KB groups PDB entries by complex composition. The same heterodimer can be deposited tens of times — different mutants, different bound ligands, different crystal forms, different resolutions, different research groups. Each deposition contributes one or more *assemblies*, and each assembly can contain one or more *interfaces*. The unit of analysis here is the interface, keyed by `(pdb_id, assembly_id, interface_id)`.

**Three coexisting numbering systems.**
- **Author numbering** — chain and residue labels chosen by the depositor, including an optional one-character insertion code. Used by both the interface API and all annotation APIs.
- **UniProt numbering** — canonical sequence position. Provided inline by the interface API for interface residues. *Not* provided by any annotation API.
- **Label numbering (SEQRES position)** — provided by some endpoints but not used in this workflow.

The notebook joins annotations to interfaces in author space (because annotation endpoints do not return UniProt mappings), then reads the UniProt key off the interface side for cross-structure analysis. This two-step join is the load-bearing mechanic of the workflow.

**Why UniProt-keyed comparison is necessary.** Author chain identifiers and residue numbers are inconsistent across PDB entries. A residue can be `A.123` in one structure and `B.123` (or `A.45`) in another, even when it is biologically the same position. UniProt accession plus UniProt sequence position is the only stable cross-structure identifier. The role label (1 or 2) disambiguates the two sides of homodimer interfaces.

**Why interaction pairs, not just interface residues.** Two interfaces can share their residue rosters but differ in *how* those residues contact each other. Comparing on (residue1, residue2, bond_type) tuples captures the contact geometry at residue resolution and is the right level of detail for the sample sizes involved (~10–30 interfaces per complex).

**Why Jaccard + hierarchical clustering and not ML.** Sample sizes are small. Jaccard on set membership is interpretable, parameter-free at this scale, and tolerates the heterogeneous interface sizes that arise across crystal forms. Hierarchical clustering with average linkage on `1 − Jaccard` distance produces a dendrogram the user can inspect directly rather than a black-box assignment.

**Annotations and what they mean.**
- *Engineered mutations* are deliberate substitutions introduced by the depositors relative to the canonical UniProt sequence — disease-associated (e.g. KRAS G12C), stability/solubility, catalytic knockouts, phosphomimics, crystallisation aids. The mutation API also returns categories that are *not* deliberate experimental choices (sequence conflicts, cloning artefacts, expression tags). The default filter retains only `"Engineered mutation"`; the filter is exposed as a parameter.
- *Modifications* are PTMs and chemical modifications encoded as alternative chem_comp_ids (e.g. `SEP` = phospho-serine). For protein-only complexes the data is generally clean.
- *Ligands* are bound molecules — drug-like inhibitors, cofactors, substrates — and also a long tail of buffer components, cryoprotectants, counter-ions. The default ligand blocklist removes the obvious noise (`SO4`, `GOL`, `HOH`, `EDO`, `PEG`, `MPD`, `CL`, `NA`, `MG`, `ZN`, `CA`, `K`, `ACT`); carbohydrate polymers are also dropped by default. Both are exposed as parameters because some users want catalytic metals or glycosylation retained.

**Caveats the notebook must surface, not hide.**
- Annotations that fall *outside* the interface are dropped. This is correct for an interface-comparison workflow but means a structure carrying a peripheral mutation will not show that mutation in its overlap row. The structure table records counts only for interface-resident annotations.
- A cluster having no annotation correlate is an expected outcome, not a failure. It can mean the structures genuinely differ for reasons the annotation pipeline does not capture (conformational state, crystal form, construct boundaries, refinement protocol), or that the cluster is methodological. The cluster interpretation report says so explicitly and does not speculate.
- At v1 sample sizes, claims about cluster–annotation correspondence are descriptive, not statistical. No formal enrichment tests are computed.

**Misconception the workflow must prevent: that absence of annotation overlap implies absence of effect.** Many residues that drive interface differences are not annotated in PDBe. The notebook surfaces what is annotated and is silent on what is not.

---

## Input Specification

### Primary input

A single PDBe-KB complex identifier, e.g. `PDB-CPX-130306`, supplied as a string in the notebook's `Config` dataclass.

### v1 scope restriction

v1 supports **dimer complexes only**. This is an upstream constraint: the `interface_interactions` endpoint currently returns data for dimers only. The notebook validates the input by calling `complex/details` first and aborts with a `ValueError` and a one-line explanation if:

- the complex ID does not resolve, or
- the complex is not a dimer.

### External APIs

All APIs are PDBe v2 endpoints over HTTPS. No authentication is required.

| Data | Endpoint | Method |
|---|---|---|
| Complex membership and partner identities | `https://www.ebi.ac.uk/pdbe/api/v2/complex/details/{complex_id}?id_type=pdb_complex_id` | GET |
| Interface interactions (PISA contacts with UniProt mappings and interface-level metrics) | `https://www.ebi.ac.uk/pdbe/api/v2/complex/interface_interactions/{complex_id}` | GET |
| Mutations | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/mutated_AA_or_NA` | POST (JSON-encoded comma-separated string of PDB IDs) |
| Modifications | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/modified_AA_or_NA` | POST (JSON-encoded comma-separated string of PDB IDs) |
| Bound molecules | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_molecules/{pdb_id}` | GET |
| Ligand interactions | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_ligand_interactions/{pdb_id}/{chain_id}/{author_residue_number}?preserve_case=false` | GET (one per surviving ligand instance) |

### Validation logic

1. **Complex resolution.** `complex/details` must return a 200 with a non-empty body containing component data. Empty or 4xx responses raise `ValueError("Complex ID {complex_id} did not resolve")`.
2. **Dimer check.** The component count from `complex/details` must equal 2. Otherwise raise `ValueError("Complex {complex_id} is not a dimer (n components = {n}); v1 supports dimer complexes only")`.
3. **Interface response non-empty.** If `interface_interactions` returns no interfaces for the complex, raise `ValueError("No interfaces returned for {complex_id}")`. v1 does not support empty-result rendering.
4. **Assembly filter non-empty.** If every assembly in `complex/details` has bound macromolecules (so the assembly filter drops every interface), raise `ValueError("All assemblies for {complex_id} have bound macromolecules; nothing to analyse after applying the assembly filter")`.

### Resolution filtering (optional)

`Config.max_resolution: float | None = None` is an optional cap on assembly resolution. When set (e.g. `2.0`), assemblies are dropped if:
- their `resolution` value (from `complex/details`) exceeds the threshold, OR
- their `resolution` value is `None` (NMR, unreported, structures without a crystallographic resolution).

The strict interpretation of "resolution N or better" is intentional — when the user filters by resolution they typically want only structures with explicit, sufficient resolution data. Default `None` disables the filter (all resolutions kept, including None).

Useful for: trimming low-resolution cryo-EM entries from datasets like Spike–ACE2 (which contains entries at 4.30 Å), or restricting analysis to high-resolution X-ray data when fine contact geometry matters.

### Assembly filtering

Interfaces from assemblies with bound macromolecules are excluded. The `complex/details` response provides a `bound_macromolecules` list per assembly; assemblies whose list is non-empty are dropped, along with every interface they would contribute.

The rationale is that a bound macromolecule (a co-chaperone, an antibody Fab, a third protein chain not part of the canonical dimer) can perturb the dimer interface either directly (by contacting interface residues) or allosterically (by stabilising a non-native conformation). Excluding these assemblies keeps the cross-structure comparison clean: every contributing interface is the canonical dimer in a state without third-party macromolecular perturbation.

Bound *small-molecule* ligands are handled separately by the ligand workflow (§8 of the source spec) and are not the subject of this filter. Many complexes carry catalytic cofactors or substrate analogues by design, so excluding any assembly with a small-molecule ligand would discard most of the dataset.

Implementation: a set of valid `(pdb_id, assembly_id)` tuples is built from `complex/details` in Phase 1. Interfaces from `interface_interactions` are filtered against this set before any downstream processing. Dropped interfaces are logged at WARNING with their `(pdb_id, assembly_id)` and the count.

The filter is exposed as a configuration parameter (`Config.exclude_assemblies_with_bound_macromolecules: bool = True`) so users who want to retain such assemblies can override it.

### Warning conditions (logged at WARNING, do not halt)

- Interface residues lacking a UniProt mapping in the response. Drop them, log a count, and record the count per interface in the structure table.
- `(unp_accession_1, unp_accession_2)` ordering inconsistent across interfaces of the same complex. Log a warning, reverse role assignment for the inconsistent entries, then proceed.
- Bound molecules dropped by the carbohydrate-polymer filter or the ligand blocklist. Log the dropped chem_comp_ids per PDB ID at INFO; this is an audit trail for the filter, not a problem.

### Preprocessing requirements

- **Field-name normalisation on ingest.** Annotation endpoints use `chain_id` and `author_residue_number`; the interface endpoint uses `auth_asym_id_{1,2}` and `auth_seq_id_{1,2}`. Both are normalised to a common author key shape `(pdb_id, auth_asym_id, auth_seq_id, ins_code)`.
- **Insertion-code normalisation.** Endpoints emit `""`, `" "`, or alphabetic codes. Apply `s.strip() or None` to every insertion-code value; whitespace-only is treated as "no insertion code."
- **Field rename.** The `interface_interactions` response calls the PDB identifier `entry_id`. Rename to `pdb_id` on ingest for consistency.
- **Role assignment from PISA.** Residues under `auth_asym_id_1` / `unp_accession_1` are role 1; residues under `auth_asym_id_2` / `unp_accession_2` are role 2. Applied uniformly to heterodimers and homodimers.

### Error handling and retry policy

The API layer distinguishes **transient** errors (where retrying is appropriate) from **persistent** errors (where bailing is appropriate):

- **Transient errors** — `ConnectionError`, request `Timeout`, and HTTP statuses **429 / 500 / 502 / 503 / 504** — trigger up to 5 total attempts (1 initial + 4 retries) with exponential backoff (1 s, 2 s, 4 s, 8 s, 16 s — total ~31 s wait). Each retry is logged at WARNING with the endpoint and attempt count. PDBe under load occasionally has 10–30 s hiccups; this policy makes most such cases recoverable.
- **Best-effort calls** — `bound_molecules` and `bound_ligand_interactions` are per-PDB-entry / per-ligand-instance calls (often 200+ in a single Phase 1). They use `tolerate_failure=True`: when retries exhaust on a transient error, the call returns an empty result with a logged WARNING rather than raising. Losing one ligand-instance fetch loses information about that specific ligand only; failing the whole analysis would be disproportionate. Critical calls (`complex/details`, `interface_interactions`, batch mutations / modifications) remain bail-fast.
- **Persistent errors** — 4xx other than 404, JSON decode failures, exhaustion of retries — raise the underlying `requests` exception with a message naming the failed endpoint and inputs.
- **404 as "no data"** is the documented behaviour for `bound_molecules`, `bound_ligand_interactions`, `mutated_AA_or_NA`, and `modified_AA_or_NA`. PDBe returns 404 (not an empty 200) when no entries have data of that type. These endpoints are called with `allow_404=True`, which converts 404 into an empty result and logs a WARNING; no retry is attempted.
- **404 on `complex/details` or `interface_interactions`** is fatal — these are required inputs.

### Batched POST endpoints

`mutated_AA_or_NA` and `modified_AA_or_NA` accept a comma-separated list of PDB IDs as their JSON body. To bound request size and protect against future server-side limits, the workflow chunks PDB IDs into blocks of **50 per request** (configurable as `POST_BATCH_SIZE` in `api.py`) and merges per-chunk results into a single `{pdb_id: [records]}` dict. Each chunk is logged at INFO with its index and total count. Empty chunks (no data for those entries) contribute nothing to the merged result.

This change makes the workflow safe at scale (large multimer complexes, when v2 supports them, may reach into hundreds of PDB IDs) without materially affecting the small-complex case (a single chunk of 50 covers most v1 dimer datasets).

### Edge cases

- **PDB entries with multiple assemblies, or multiple interfaces per assembly,** contribute multiple rows to the structure table. This is correct behaviour, not deduplication.
- **Homodimers.** Both partners share a UniProt accession; role distinguishes them. Symmetric contacts (e.g. `(chain1:25, chain2:87)` and `(chain1:87, chain2:25)`) are preserved as distinct ordered tuples in the interaction set, which is the correct physical behaviour.
- **No mutations or no ligands returned.** Not an error. Annotation lists are empty; structure table columns are zero.
- **Two distinct author residues mapping to the same UniProt position** (microheterogeneity, alternate residue identities). The workflow does not detect this — it is documented in §11 of the source spec as a load-bearing assumption that SIFTS provides a one-to-one mapping per `(pdb_id, auth_asym_id, auth_seq_id, ins_code)`.

---

## System Execution / Processing Logic

### Tech stack

- Python 3.11+
- `requests` for HTTP (synchronous; v1 call volume is small)
- `pandas` for tables
- `numpy` for the similarity matrix
- `scipy.cluster.hierarchy` for clustering and dendrogram
- `matplotlib` + `seaborn` for the heatmap and dendrogram
- Standard library: `dataclasses`, `typing`, `logging`, `json`, `functools`

No additional dependencies unless a concrete need arises during the build.

### Installation / setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install jupyter requests pandas numpy scipy matplotlib seaborn
```

The notebook expects the `pdbe_interfaces/` package (see Component Breakdown) to be importable from the working directory. No external data files are required; all inputs are fetched at runtime.

### Configuration

A single `Config` dataclass is defined at the top of the notebook. Users edit the dataclass instance, not function arguments. All defaults match the spec:

```python
@dataclass
class Config:
    complex_id: str = "PDB-CPX-130306"
    mutation_type_filter: tuple[str, ...] = ("Engineered mutation",)
    ligand_blocklist: frozenset[str] = frozenset({
        "CL", "NA", "MG", "ZN", "CA", "K",
        "SO4", "GOL", "HOH", "EDO", "PEG", "MPD", "ACT",
        "SGM", "DTT", "BME", "TRS",
    })
    drop_carbohydrate_polymers: bool = False
    exclude_assemblies_with_bound_macromolecules: bool = True
    max_resolution: float | None = None  # e.g. 2.0 to keep only ≤ 2 Å
    conservation_threshold: float = 0.8     # fraction of interfaces; 1.0 is strict
    cluster_distance_cut: float = 0.5       # 1 - Jaccard
    log_level: str = "INFO"
    output_dir: str | None = None           # JSON export destination; None = current dir
```

**Parameter rationale.**
- `mutation_type_filter`: deliberate experimental choices only; cloning artefacts and expression tags are noise for interface comparison.
- `ligand_blocklist`: removes buffers, cryoprotectants, common counter-ions. Users retaining catalytic metals (e.g. Mg²⁺ in kinases) override.
- `drop_carbohydrate_polymers`: glycosylation is biologically meaningful for some complexes but out of scope by default.
- `conservation_threshold = 0.8`: tolerates 1–2 drop-outs at N≈10–15. Set to 1.0 for strict conservation; 0.5 shifts the question from "conserved" to "majority-present."
- `cluster_distance_cut = 0.5`: a midpoint on `1 − Jaccard` that separates clearly different interfaces from clearly similar ones at the working-example sample size. Users override directly or read clusters off the dendrogram.
- `output_dir`: where the lightweight JSON export (Phase 7) is written. `None` falls back to the current working directory; missing directories are created with `pathlib.Path.mkdir(parents=True, exist_ok=True)`.

### Execution flow — phases

The notebook is organised as a sequence of phases, each mapping to one or two cells. Every cell ends with an inspectable artefact (a dataframe, a plot, a printed summary). Phases 1–5 form the analysis core; Phase 7 is a self-contained export with no effect on prior outputs.

**Phase 1 — Data retrieval.**
1. Validate the input by calling `complex/details`. Abort if invalid. Build the `partner_map` from the response. Build the set of valid `(pdb_id, assembly_id)` tuples — those whose `bound_macromolecules` is empty (see Assembly filtering above).
2. Call `interface_interactions` for the complex. Apply the assembly filter; drop interfaces from assemblies with bound macromolecules and log the count.
3. Call annotation endpoints in parallel (logical parallelism — the notebook executes sequentially in v1, but the four classes of fetch are independent and could be reordered): batched POSTs for mutations and modifications covering the PDB IDs from the surviving interfaces; per-entry GETs for `bound_molecules`; per-surviving-ligand GETs for `bound_ligand_interactions` after filtering.

**Phase 2 — Build representations.**
4. For each interface in the response, build the author-keyed and UniProt-keyed interaction-pair sets, including the atom-level → residue-pair-and-bond-type aggregation.
5. Drop interface residues lacking a UniProt mapping, log a warning, record the count per interface.
6. Run the cross-structure consistency check on `(unp_accession_1, unp_accession_2)`. Reverse role assignment for inconsistent interfaces.

**Phase 3 — Similarity and clustering.**
7. Compute pairwise Jaccard similarity between all interfaces using the UniProt-keyed, typed interaction-pair sets. Compute an untyped variant (collapsing across `bond_type`) as a robustness check.
8. Cluster by hierarchical clustering with average linkage on `1 − Jaccard` distance. Render the dendrogram. Default cut from `Config.cluster_distance_cut`.

**Phase 4 — Annotation overlap.**
9. For each interface, join annotations (mutations, modifications, ligand-contact residues) to interface residues using the author key. Drop annotations that fall outside the interface.
10. After the join, attach the UniProt key from the interface side for cross-structure comparison.

**Phase 5 — Summarisation and interpretation.**
11. Build the structure table.
12. Identify conserved residues and conserved interaction pairs at the configured threshold.
13. Build the cluster interpretation report — one row per interface interaction state, including per-interface distributions of residue-pair count, typed-tuple interaction count, and PISA interface area; median contact density; **core contacts** present in ≥ `conservation_threshold` of cluster members; enriched contacts and annotations; and QC warnings (singleton, sparse, tiny interface, mixed UniProt residue range vs the dominant cluster).
14. Build the **interface rewiring** table — a head-to-head comparison between two interaction states, labelling each contact as `shared core`, `A-enriched`, `B-enriched`, or `rare`. Defaults to the two largest non-singleton states; users can pass arbitrary cluster ids.

**Phase 7 — Export interface conservation JSON.**
15. Call `export_interface_frequency_json` to write a lightweight JSON containing UniProt residue-level interface conservation and residue–residue contact frequencies for downstream visualisation (Mol* colouring, frontend contact matrices/heatmaps, AFDB-style overlays). This step is export-only — it does not modify clustering, the heatmap, or the interpretation report. Output is written to `Config.output_dir` (created if missing) or the current directory if `output_dir` is `None`. See **Output schema — Interface frequency JSON export** below for the file structure.

### Caching

v1 has no caching. API responses are fetched fresh on each notebook run; the structure table records the fetch date so outputs remain interpretable across PDBe releases. Disk caching is out of scope for v1 (see Out of Scope).

### Performance expectations

For PDB-CPX-130306 (~12 PDB entries, ~10–20 interfaces total), end-to-end runtime is dominated by the per-entry `bound_molecules` and per-ligand `bound_ligand_interactions` GETs — typically tens of HTTP round trips, on the order of 30–120 seconds depending on PDBe latency. The Jaccard matrix and clustering are negligible at this scale.

### Invocation

The notebook is the deliverable. There is no CLI in v1. Users open `notebook.ipynb`, edit the `Config` instance in the first cell, and run all cells.

---

## Output Parsing and Interpretation

### Interaction-pair tuple schema (the central data structure)

**Author-numbered (within-structure, for annotation overlap):**
```
(
  (pdb_id, auth_asym_id_1, auth_seq_id_1, ins_code_1),
  (pdb_id, auth_asym_id_2, auth_seq_id_2, ins_code_2),
  bond_type
)
```

**UniProt-numbered (cross-structure, for similarity and clustering):**
```
(
  (unp_acc_1, unp_seq_id_1, 1),
  (unp_acc_2, unp_seq_id_2, 2),
  bond_type
)
```

`bond_type` takes values from the PISA vocabulary (e.g. `hydrogen_bond`, `salt_bridge`, `covalent_bond`, `disulfide_bond`, plus others as returned). Values are recorded as they appear in API responses; **no manual remapping**.

### Aggregation rule

Each row in the PISA response is an atom-level contact. Atom-level contacts sharing the same key (residue pair plus `bond_type`) collapse to a single set element. A residue pair forming contacts of different `bond_type` values appears as multiple elements. Atom-level detail is not retained.

Aggregation is per interface and per key. For homodimers, ordering preserves both directions of a symmetric contact as distinct elements (correct physical behaviour).

### Output schema — Structure table

One row per interface, keyed on `(pdb_id, assembly_id, interface_id)`. Columns:

- **Identification**: `pdb_id`, `assembly_id`, `interface_id`, `partner_1`, `partner_2` (gene names from `partner_map`).
- **Methodology**: `experimental_method` (string, e.g. "X-ray diffraction", "Electron microscopy"), `resolution` (float, Å). Both are taken per-assembly from `complex/details`. Surfaced alongside biological annotations so the user can distinguish methodological clusters (e.g. "all members are 3.5 Å cryo-EM") from biological ones (e.g. "all members are Omicron").
- **PISA metrics**: `interface_area`, `solvation_energy`, `stabilisation_energy`, `pisa_p_value`, `n_interface_residues`, `n_hydrogen_bonds`, `n_salt_bridges`, `n_covalent_bonds`, `n_disulfide_bonds`.
- **Annotation counts**: `n_engineered_mutations_at_interface`, `n_modifications_at_interface`, `n_distinct_ligands_at_interface`, `n_ligand_contact_residues_at_interface`.
- **Cluster assignment**: `cluster_id` at the default cut.
- **Quality**: `n_residues_dropped_no_uniprot`, `fetch_date` (ISO-8601 date).

### Output schema — Similarity matrix

A square `numpy` array of shape `(n_interfaces, n_interfaces)` of pairwise Jaccard similarities computed on the UniProt-keyed, typed interaction-pair sets. Symmetric. Diagonal is exactly 1.0. Index order matches the structure table row order.

A second matrix is computed alongside using the *untyped* sets (bond_type collapsed) as a robustness check. Substantial divergence between the two suggests bond-type assignment is doing more work than the residue-pair topology — informative but not actionable in v1.

### Output schema — Clusters

`scipy.cluster.hierarchy.linkage` output (an `(n−1, 4)` array) plus a flat assignment vector at `cluster_distance_cut`. The flat vector populates the `cluster_id` column of the structure table.

A **cluster-cut sweep helper** (`similarity.sweep_cluster_cuts`) returns cluster count and top-5 cluster sizes at standard cuts (0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70). The notebook displays this table alongside the dendrogram so the user can pick `cluster_distance_cut` from data: a stable plateau (same cluster count across a range of cuts) is a defensible default; a fragmentation gradient calls for inspecting the dendrogram. The default of 0.5 is a sensible starting point but is not guaranteed to be optimal for every complex.

### Output schema — Conserved residues and conserved interaction pairs

Two sets:
- Conserved residues: UniProt-keyed residues `(unp_acc, unp_seq_id, role)` present in at least `conservation_threshold` × N interfaces.
- Conserved interaction pairs: UniProt-keyed interaction-pair tuples present in at least `conservation_threshold` × N interfaces.

Default threshold 0.8.

### Output schema — Annotation overlap views

Three long-form dataframes (one per stream — mutations, modifications, ligands), each with columns:

- `pdb_id`, `assembly_id`, `interface_id`
- `auth_asym_id`, `auth_seq_id`, `ins_code`
- `unp_acc`, `unp_seq_id`, `role`
- Stream-specific metadata: `mutation_label` (e.g. `G12C`), or `modification_chem_comp_id` + `modification_chem_comp_name`, or `ligand_chem_comp_id` + `ligand_chem_comp_name` + `ligand_chain_id` + `ligand_author_residue_number` + `contact_type`.

These pivot cleanly into structure-by-residue heatmaps for visual inspection.

### Cross-cluster comparison helpers

Two helpers are exposed; both produce a per-contact diff between two clusters but with different framings:

**`outputs.compare_cluster_contacts(records, cluster_result, cluster_a=None, cluster_b=None, typed=True, top_n=20, partner_map=None, ...)`** is the **interface rewiring** helper. It returns a DataFrame with columns `contact`, `cluster_A_count`, `cluster_B_count`, `cluster_A_fraction`, `cluster_B_fraction`, `fraction_difference`, `enrichment_direction`, sorted by `|fraction_difference|` descending and capped at `top_n`. The direction label categorises each contact as:

- `shared core` — fraction ≥ `shared_core_threshold` (default 0.8) in *both* clusters; the conserved core that is unchanged between states;
- `A-enriched` — `fraction_difference ≥ a_enriched_threshold` (default 0.5); a contact strongly gained in state A relative to state B;
- `B-enriched` — `fraction_difference ≤ -a_enriched_threshold`; a contact strongly gained in state B;
- `rare` — neither state shows strong enrichment (the long tail).

When `cluster_a` and `cluster_b` are omitted, the two largest non-singleton clusters are selected automatically. When `typed=False`, contacts collapse across `bond_type` so a residue-pair shared via different bond types counts once. The chosen cluster ids and sizes are surfaced via `df.attrs["cluster_a"] / "cluster_b" / "size_a" / "size_b" / "typed"` for downstream display.

**`outputs.compare_clusters(records, cluster_result, cluster_a, cluster_b, partner_map=None, top_n=20, ...)`** is the original divergence helper, retained for backward compatibility. It returns a DataFrame ranked by `magnitude = |fraction_in_a - fraction_in_b|` with directional labels `A-only`, `B-only`, `biased_to_A`, `biased_to_B`, `tied`.

Both are the natural follow-up to the cluster interpretation report: the report tells you *what each interaction state is enriched for*; the comparison tells you *what specifically distinguishes one state from another*. Worked example on ACE2–Spike (cluster 13 = Omicron, cluster 25 = engineered ACE2):

- `ACE2:353-S:496 hydrogen_bond` is present in 22/23 cluster-25 members and 0/34 cluster-13 members (`B-enriched`, fraction_difference ≈ −0.96).
- `ACE2:19-S:477 hydrogen_bond` is present in 30/34 cluster-13 members and 0/23 cluster-25 members (`A-enriched`, fraction_difference ≈ +0.88) — the S477N variant signature.
- `ACE2:83-S:489 hydrogen_bond` is present in nearly every member of both clusters — flagged `shared core` because it is part of the conserved core that is unchanged between states.

The comparison directly surfaces the structural distinction between two affinity-enhancement strategies (variant mutation vs. ACE2 engineering): variants strengthen existing contacts, engineering introduces new contacts at different positions.

### Output schema — Cluster interpretation report

A pandas DataFrame, one row per **interface interaction state** (cluster), with columns:

- `cluster_id`, `cluster_size`, `member_pdb_ids`
- `experimental_methods` — comma-joined `"method (count)"`, e.g. `"X-ray diffraction (24), Electron microscopy (10)"`. Per-interface counts.
- `resolution_range` — `"min–max Å (median X.XX, n=N)"`. `n` may be less than `cluster_size` if some assemblies lack a resolution value.
- `interface_area_range` — `"min–max Å² (median X, n=N)"`. Sourced from PISA `interface_area` per interface.
- **Per-interface distributions** (one numeric column each — `min` / `median` / `mean` / `max`), computed across the cluster's members:
  - `residue_pair_count_*` — distinct UniProt-keyed residue–residue pairs per interface (collapsed across `bond_type`).
  - `interaction_count_*` — distinct typed-tuple interactions per interface (the input to Jaccard similarity).
  - `interface_area_*` — PISA `interface_area` per interface in Å² (None when missing).
- `contact_density_median` — `interaction_count_median / interface_area_median`, in interactions/Å². None when either median is missing or zero. A coarse density measure that distinguishes "compact, contact-rich" interfaces from "extended, sparse" ones.
- `core_contacts` — semicolon-joined `"PartnerA:X{pos}-PartnerB:Y{pos} bond_type (count/denominator, pct%)"` for typed contacts present in ≥ `conservation_threshold` of cluster members. Sorted by `count` descending, then by contact label. Capped at the top 10. Complements `enriched_contacts`: the core surface is *what is conserved within the state*; enriched contacts are *what distinguishes the state from the rest of the dataset*.
- `enriched_contacts` — semicolon-joined `"PartnerA:X{pos}-PartnerB:Y{pos} bond_type (count, Nx)"` where `X` / `Y` are the canonical UniProt one-letter amino-acid codes (e.g. `KRAS:D33-RAF1:K84 salt_bridge`, `ACE2:D38-S:N501 hydrogen_bond`). Identifies UniProt-keyed interaction-pair tuples enriched in this cluster vs the rest of the dataset (same enrichment heuristic as annotations). Top `top_n_contacts` per cluster (default 10). Partner labels come from `partner_map`; residue identities come from `unp_one_letter_code_{1,2}` in the interface API response, captured into `InterfaceRecord.residue_identity` during Phase 2. The one-letter code is the canonical UniProt residue, not the deposited residue — variants are surfaced separately via the mutation overlap. This surfaces interface-level cluster signatures directly: e.g. for ACE2–Spike cluster 13, `ACE2:D38-S:N501 hydrogen_bond (12, 8.5×)` corresponds exactly to the N501Y Omicron mutation. The contact enrichment frequently surfaces biological cluster signatures even when no mutation/modification/ligand annotation is recorded.
- `enriched_mutations`, `enriched_modifications`, `enriched_ligands` — semicolon-joined `"label (count, Nx)"` or `"label (count, exclusive)"`
- `interfaces_with_any_mutation`, `interfaces_with_any_modification`, `interfaces_with_any_ligand` — count of cluster-member interfaces carrying at least one at-interface annotation of that stream
- `annotation_correlate_present` (bool) — true if any stream has at least one enriched label
- `qc_warnings` — semicolon-joined warnings; empty string when none. See *Cluster QC warnings* below.
- `notes` (free-text narrative summary in `"interface interaction state {id} ({n} interfaces)"` form, with the methodology summary as a `[methods: ...; resolution: ...]` suffix and the per-state statistics, core contacts, and QC warnings appended on indented lines for printable display)

### Cluster QC warnings

The cluster interpretation report attaches automatic warnings to flag clusters that are likely to be misinterpreted. Each warning is included in the row's `qc_warnings` column (semicolon-joined) and surfaced in `notes`. The four warnings:

- **Singleton** (`"singleton cluster; interpret cautiously"`) — `cluster_size == 1`. A single interface cannot establish a recurring state on its own.
- **Sparse fingerprint** (`"sparse interaction fingerprint; clustering may reflect few contacts"`) — `residue_pair_count_median < sparse_pair_threshold` (default 5). Clusters built from very few residue pairs per interface are sensitive to single-bond differences and unreliable.
- **Tiny interface** (`"small interface area; possible weak or non-biological interface"`) — `interface_area_median < tiny_area_threshold_a2` (default 500 Å²). Below this is the regime of crystal contacts and weak transient interactions, not stable biological dimers.
- **Mixed UniProt residue range** (`"residue range poorly overlaps the dominant cluster; this may indicate a different domain, processed product, or polyprotein segment"`) — for each `(unp_accession, role)` in the cluster, the `(min, max)` UniProt position is compared against the dominant (largest) cluster's range; if the overlap fraction (overlap length / span of the larger interval) is below `range_overlap_min` (default 0.2) for any partner, the warning fires. This catches polyprotein / processed-product / domain mixups, e.g. HIV gag-pol where different processed proteins map to the same UniProt accession at very different residue ranges. The dominant cluster never receives this warning (it is the reference).

All thresholds are exposed as `cluster_interpretation_report` parameters (`sparse_pair_threshold`, `tiny_area_threshold_a2`, `range_overlap_min`) so they can be tuned per dataset.

### Output schema — Interface frequency JSON export

A self-contained JSON file written to `Config.output_dir` (or the current working directory if `output_dir` is `None`) by `export_interface_frequency_json`. The filename is `{complex_id}_interface_frequencies.json`. The file is consumed by frontend visualisation tooling (Mol* colouring by conservation, residue–residue contact matrices/heatmaps, AFDB-style overlays) and is independent of clustering. The export is **always** produced from the full set of retained interfaces — interaction states are not partitioned in this view.

**Top-level structure.**

```json
{
  "metadata": { ... },
  "partners": [ ... ],
  "residue_frequencies": [ ... ],
  "contact_frequencies": [ ... ]
}
```

**`metadata` block.**

| Field | Description |
|---|---|
| `complex_id` | input complex identifier (e.g. `"PDB-CPX-140195"`) |
| `complex_name` | human-readable complex name from `complex/details` |
| `oligomeric_state` | e.g. `"Heterodimer"` |
| `generated_on` | ISO date the export was written |
| `n_interfaces` | total retained interfaces (denominator for all frequencies) |
| `representation` | fixed string: `"UniProt residue-level interface interactions"` |
| `contact_definition` | fixed string: `"PISA-derived residue-residue interface interactions mapped to UniProt residue positions"` |
| `typed_contacts` | bool — true when bond type is part of the contact key |

**`partners` array.** One object per `(unp_accession, role)` pair. Roles observed in the records and in `partner_map` are unioned, so homodimers (one accession, two roles) and heterodimers both expose role 1 *and* role 2:

```json
{
  "role": 1,
  "unp_accession": "Q9BYF1",
  "gene_name": "ACE2",
  "name": "Angiotensin-converting enzyme 2",
  "label": "ACE2"
}
```

`gene_name` and `name` are populated from `complex_details["participants"]` when the caller passes `complex_details`; otherwise they are `null`. The `label` field follows a fallback rule: `gene_name` → `name` → `partner_map` value with `" (role 1)"` / `" (role 2)"` suffixes stripped → `unp_accession`. No additional API calls are made.

**`residue_frequencies` array.** One object per UniProt residue observed at any interface, with `frequency = n_interfaces_with_residue / metadata.n_interfaces`. Each residue is counted at most once per interface. Sorted deterministically by `(role, unp_accession, unp_residue_number)`:

```json
{
  "unp_accession": "Q9BYF1",
  "role": 1,
  "unp_residue_number": 38,
  "unp_residue_label": "D38",
  "n_interfaces": 120,
  "frequency": 0.9231,
  "conservation_level": "strong"
}
```

`unp_residue_label` is the canonical UniProt one-letter code at this position concatenated with the residue number (e.g. `"D38"`). `conservation_level` follows fixed cutoffs:

| Bucket | Frequency range |
|---|---|
| `strong` | `>= 0.80` |
| `medium` | `0.50 <= frequency < 0.80` |
| `weak`   | `0.20 <= frequency < 0.50` |
| `rare`   | `< 0.20` |

`conservation_level` is **not** included on `contact_frequencies` rows; conversely, `residue_frequencies` rows do **not** carry `contact_types` or bond information.

**`contact_frequencies` array.** One object per UniProt-keyed residue–residue contact observed at any interface, with `frequency = n_interfaces_with_contact / metadata.n_interfaces`. Each contact is counted at most once per interface. When `typed_contacts=True` (default) the contact key includes `bond_type`, so a residue pair seen with both a hydrogen bond and a salt bridge contributes two rows; when `typed_contacts=False`, bond types are collapsed and `bond_type` is omitted from the row. Sorted by `(partner_1.role, partner_1.unp_residue_number, partner_2.role, partner_2.unp_residue_number, bond_type)`:

```json
{
  "partner_1": {
    "unp_accession": "Q9BYF1",
    "role": 1,
    "unp_residue_number": 38,
    "unp_residue_label": "D38"
  },
  "partner_2": {
    "unp_accession": "P0DTC2",
    "role": 2,
    "unp_residue_number": 498,
    "unp_residue_label": "Q498"
  },
  "bond_type": "salt_bridge",
  "n_interfaces": 26,
  "frequency": 0.51
}
```

`partner_1` is always the role-1 partner and `partner_2` the role-2 partner; partners are not reordered alphabetically or by accession. For homodimers the role assignment chosen during representation building is preserved.

**Determinism.** Sort orders are explicit on both `residue_frequencies` and `contact_frequencies`; output is independent of Python set/dict iteration order. Preview tables shown in the notebook (top-N by frequency descending) sort the in-memory copy for readability — they do not affect the on-disk file.

**Constraints.**

- This is an export-only feature. It does not modify clustering, the heatmap, the cluster interpretation report, or any earlier output.
- No additional web/API calls are made to enrich partner names; the function consumes only what is already in `records`, `partner_map`, and the optionally passed `complex_details`.

### Enrichment heuristic

For each cluster and each annotation stream, a label is reported when:

- It appears on at least `min_count` distinct cluster-member interfaces (default 2) AND its proportion-of-interfaces in the cluster is at least `min_enrichment` × its proportion in the rest of the dataset (default 2.0); OR
- It is **exclusive** to this cluster — zero occurrences in the rest of the dataset — and appears on at least one member interface.

Counts are by distinct interface (de-duplicated on `(pdb_id, assembly_id, interface_id)`), so a ligand with multiple atom-level contacts on one interface counts once.

This replaces the earlier "must dominate ≥60% of cluster members" rule, which missed combinations of low-frequency markers. Variant-class signatures — where each individual mutation hits 1–3 cluster members but together define the cluster — are now correctly surfaced. Worked example (PDB-CPX-140195, ACE2–Spike): cluster 13 (n=34) is flagged as carrying S477N (3 interfaces, 8.5×), Q498R (3, 8.5×), Q493R (1, exclusive) — the canonical Omicron RBD signature.

Both `min_count` and `min_enrichment` are exposed as `cluster_interpretation_report` parameters; the defaults (2 and 2.0) are tuned for v1 sample sizes (~10–150 interfaces per complex). At larger N a stricter enrichment threshold may be appropriate.

### Confidence interpretation and prioritisation rules

- **Jaccard similarity** is interpreted as the fraction of interaction-pair-and-bond-type tuples shared between two interfaces. It is a descriptive measure, not a probability.
- **Cluster compactness** at the default cut is informative but not load-bearing. The dendrogram is the source of truth; the flat cut is a convenience.
- **Conservation threshold** is a knob for the user. Default 0.8. Setting to 1.0 produces the strict conserved core; setting to 0.5 produces the majority-present set, which is a different question.
- **Annotation correlate present** means the cluster has a distinctive feature (mutation, ligand class, modification) that distinguishes it from other clusters. Absence of correlate is *not* a failure; the report says so explicitly.

### Failure handling and empty-result behaviour

- Empty annotation streams (no mutations, no ligands) produce empty long-form dataframes and zero counts in the structure table — not errors.
- A cluster with no annotation correlate produces a row in the cluster interpretation report with `annotation_correlate_present = False` and a `notes` value such as `"No distinguishing mutation, modification, or ligand identified for this cluster's members."`
- An interface with all residues dropped for missing UniProt mapping (pathological) is excluded from the similarity matrix and flagged in the structure table with `n_residues_dropped_no_uniprot` equal to the original residue count and zero contribution to clustering.

### Ambiguity handling

- If two interfaces have identical UniProt-keyed interaction-pair sets but differ in author space (e.g. different chain labels), they are correctly identical at Jaccard 1.0 and will collapse to a single cluster at any positive distance threshold.
- If `(unp_accession_1, unp_accession_2)` ordering is inconsistent across interfaces, role assignment is reversed for the inconsistent entries before aggregation. This is not arbitrary — the partner under `unp_accession_1` for the *majority* of interfaces is taken as the canonical role-1 partner.

---

## Visualisations / UI Components

### Similarity heatmap

- **Purpose.** Quick visual read of which interfaces are similar and where the obvious clusters lie.
- **Data source.** The UniProt-keyed, typed Jaccard matrix.
- **Library.** `seaborn.heatmap`.
- **Row/column ordering.** Reordered by dendrogram leaf order (`scipy.cluster.hierarchy.leaves_list`) so cluster blocks appear as bright squares on the diagonal. Without this reordering, interfaces are in API-fetch order and the cluster structure is invisible.
- **Cluster boundaries.** Thin white lines drawn between neighbouring leaves whose flat-cluster ID changes, marking the block structure explicitly.
- **Defaults.** `viridis` colormap, square cells, value annotations off (Jaccard values are readable from the colour scale at this matrix size).
- **Row/column labels.** `pdb_id` plus a short annotation suffix where present, e.g. `6vjj [G12C, sotorasib]`. Where the suffix would exceed ~30 characters, truncate with an ellipsis.
- **Interaction.** Static. The notebook uses inline matplotlib output; no zoom/pan.
- **Fallback if data missing.** If the matrix is 1×1 (single interface), skip the heatmap and print `"Only one interface; similarity matrix is trivial."`

### Dendrogram

- **Purpose.** The primary visual for cluster structure. Lets the user choose a distance cut by eye.
- **Data source.** `scipy.cluster.hierarchy.linkage` output on `1 − Jaccard` distance.
- **Library.** `scipy.cluster.hierarchy.dendrogram` rendering through matplotlib.
- **Defaults.** Default orientation (top-down). No automatic colouring of clusters in v1. A horizontal line is drawn at `cluster_distance_cut` to show the cut.
- **Leaf labels.** Same as heatmap row labels.
- **Interaction.** Static.
- **Fallback if data missing.** If fewer than two interfaces, skip the dendrogram and print `"Need at least two interfaces to cluster."`

### Structure table

- **Purpose.** One-look summary of every interface: identification, PISA metrics, annotation counts, cluster assignment, quality.
- **Data source.** The structure-table dataframe (see Output Parsing).
- **Library.** `pandas.DataFrame.style` for the notebook display. No exotic formatting.
- **Defaults.** All numeric columns formatted to sensible precision (`interface_area` to 0 decimals, energies to 1 decimal, Jaccard-derived columns to 3 decimals). Cluster ID rendered as integer.
- **Highlighting.** Optional row-shading by cluster ID via `df.style.apply` to make cluster membership scannable. Off by default; enabled with one line if the user wants it.

### Mol* / MolViewSpec interface viewer

3D viewers for interface inspection, rendered via the `pdbe_interfaces.visualize` module:

- `visualize_interface(record, overlap=None)` — render one interface.
- `visualize_cluster_representative(cluster_id, records, cluster_result, overlap=None, assembly_metadata=None)` — pick a representative and render it.
- `visualize_clusters_grid(records, cluster_result, overlap=None, assembly_metadata=None, min_cluster_size=2, max_clusters=8)` — render one representative per cluster, stacked vertically, with a markdown header naming each. Defaults skip singletons and cap at 8 clusters; override to widen or narrow.

**Representative selection:** when `assembly_metadata` is provided, the representative is the cluster member with the **lowest (best) resolution**. Ties broken by `(pdb_id, assembly_id, interface_id)` so the choice is deterministic across runs. Members without a resolution value (NMR, unreported) are de-prioritised. When metadata is not provided, falls back to the first member by the same deterministic sort. This avoids the pitfall of "first-by-API-order" selection where a low-resolution outlier could become the visual representative of the cluster.

**Library:** `molviewspec` (Python builder for MolViewSpec) → embedded Mol* viewer.

**Encoding:** both partner chains as a semi-transparent molecular surface in a neutral colour (default `lightgray`, opacity 0.4) — a gray backdrop so the colored interface residues stand out. Interface residues drawn as ball-and-stick in the chain's colour (`cornflowerblue` / `lightcoral`); chain identity is read from the stick colour, not the surface colour. At-interface mutations overlaid red. `stick_size_factor` (default 0.5) scales the sticks; `surface_color` is exposed for users who want per-chain coloured surfaces back (pass any colour name or set per-chain by editing `PARTNER_COLOURS`).

**Camera:** the viewer auto-focuses on the union of all interface residues with a small padding factor (`zoom_radius_factor`, default 1.3) so the binding interface is filling the viewport on first render. Both `surface_opacity` and `zoom_radius_factor` are exposed on all three viewer helpers.

**Selectors:** `auth_asym_id` and `auth_seq_id` — same author numbering used elsewhere in the workflow, so no conversion is needed.

**Source structure:** PDBe coordinate file, `https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif`.

**Display:** side-effect via IPython `display()`. Each viewer embeds a few hundred kB of HTML/JS, so the workflow recommends visualising a small number of cluster representatives rather than every interface. The grid helper enforces this by default via `min_cluster_size` and `max_clusters`.

**Fallback:** if `molviewspec` is not installed, the import fails fast with a clear `ImportError`. The rest of the workflow does not depend on it.

### Residue-pair frequency heatmap and tables

`outputs.interface_frequency_summary(records, partner_map=None, residue_identity=None)` returns a dict containing:
- `pairs` — DataFrame of UniProt-keyed residue-pairs (collapsed across bond_type), sorted by how many interfaces contain each pair. Surfaces the "anchor" contacts that are universal across the dataset.
- `partner_1_residues` / `partner_2_residues` — DataFrames listing each unique residue per partner with its count of participating interfaces. Surfaces "key residues" — the ones that always appear at the interface.
- `pair_matrix`, `partner_1_labels`, `partner_2_labels` — 2D numpy matrix (partner-1 residues × partner-2 residues) with per-cell counts, plus axis labels. Suitable for a frequency heatmap.

The notebook wraps these outputs in an **interactive widget** (`ipywidgets.Dropdown` + `IntSlider`):
- Cluster dropdown: `"All clusters"` (full dataset) or a specific cluster ID. Re-runs `interface_frequency_summary` on the cluster's member subset and re-renders tables + heatmap on change.
- Top-N slider: caps rows displayed in each table and the heatmap's matrix dimension (default 15, range 5–40).

This gives an "all interfaces vs cluster X" comparison without re-running the workflow. Cell values in the heatmap are fractions in `[0, 1]`; 1.0 means every interface in the (filtered) selection contains that residue-pair contact.

### Annotation overlap heatmaps (optional, derived)

- **Purpose.** Show which residues carry which annotations across structures.
- **Data source.** Pivot of the annotation overlap long-form dataframes: rows = UniProt residue keys, columns = `pdb_id`, cell values = annotation labels (e.g. mutation substitution).
- **Library.** `seaborn.heatmap` with categorical encoding, or `pandas.DataFrame.style` with conditional colouring.
- **Fallback if data missing.** If the stream is empty, skip the heatmap and print `"No {stream} found at any interface for this complex."`

### General plot defaults

- White background, default fonts, no titles inside the figure (titles live in the surrounding markdown cells).
- All plots inline in the notebook. No separate output files in v1.

### Accessibility considerations

- Colormap choice (`viridis`) is colour-blind-safe.
- Long axis labels are rotated 90° on the heatmap to avoid overlap.
- The cluster interpretation report is a dataframe (text-readable) rather than a graphic, so the substantive findings are accessible without relying on the visualisations.

---

## Workflow / User Journey

The user is a structural biologist or computational scientist who wants to understand how a known protein–protein complex's interface varies across deposited structures.

### Stage 1 — Open the notebook and configure

- **User action.** Open `notebook.ipynb` in Jupyter. Edit the `Config` dataclass instance in the first cell to set `complex_id` to the target PDBe-KB complex. Optionally adjust thresholds, mutation filter, ligand blocklist.
- **System action.** None yet.
- **UI response.** The first cell renders the `Config` instance for confirmation.
- **Decision points.** User decides whether to override defaults. The defaults are tuned for the working example; for other complexes the user may need to adjust the ligand blocklist (e.g. retain Mg²⁺ for kinases) or the conservation threshold.

### Stage 2 — Run Phase 1: data retrieval

- **User action.** Execute the Phase 1 cells.
- **System action.** Calls `complex/details`, validates dimer status, builds `partner_map`. Calls `interface_interactions`. Calls mutation and modification batched POSTs. Calls per-entry `bound_molecules` and per-surviving-ligand `bound_ligand_interactions`.
- **UI response.** Each cell prints a short summary: complex name, partner labels, number of interfaces retrieved, number of PDB entries, count of mutations / modifications / ligands fetched. Logged warnings appear inline (residues without UniProt mapping; ligand blocklist applications).
- **Loading state.** The notebook cell shows a running indicator. For the working example, expect 30–120 seconds.
- **Error states.**
    - Invalid complex ID: `ValueError` at the validation step, with the message visible in the cell output. The user fixes `Config.complex_id` and re-runs.
    - Network error: the underlying `requests` exception is raised. The user retries (no in-notebook retry logic in v1).
    - Non-dimer complex: `ValueError` with explanation. Out of scope for v1.

### Stage 3 — Run Phase 2: build representations

- **User action.** Execute the Phase 2 cells.
- **System action.** Builds author-keyed and UniProt-keyed interaction-pair sets per interface. Drops residues lacking UniProt mappings. Runs the cross-structure consistency check.
- **UI response.** Each cell prints: total interfaces, total interaction-pair tuples per interface (summary stats), residues dropped per interface, and any partner-order corrections applied.

### Stage 4 — Run Phase 3: similarity and clustering

- **User action.** Execute the Phase 3 cells.
- **System action.** Computes Jaccard matrix (typed and untyped). Computes linkage. Renders heatmap and dendrogram.
- **UI response.** Heatmap and dendrogram inline. Brief printed summary: number of clusters at the default cut, range of Jaccard values.
- **Decision points.** The user inspects the dendrogram. If the default cut produces too few or too many clusters, the user adjusts `Config.cluster_distance_cut` and re-runs Phase 3 (Phase 1 and 2 do not need re-running — but in v1 the simplest workflow is to re-run all cells; caching is out of scope).

### Stage 5 — Run Phase 4: annotation overlap

- **User action.** Execute the Phase 4 cells.
- **System action.** Joins each annotation stream onto interface residues by author key. Attaches the UniProt key from the interface side. Builds the three long-form annotation overlap dataframes.
- **UI response.** A printed count of matched annotations per stream. Optional inline display of each long-form dataframe (head only, full dataframe available via variable reference).

### Stage 6 — Run Phase 5: summarisation and interpretation

- **User action.** Execute the Phase 5 cells.
- **System action.** Builds the structure table. Computes conserved residues and conserved interaction pairs at the configured threshold. Builds the cluster interpretation report.
- **UI response.** Structure table rendered with `df.style`. Conserved sets printed as lists. Cluster interpretation report rendered as a dataframe.
- **Decision points.** The user reads the cluster interpretation report and decides whether the clusters are interpretable. If most clusters report `annotation_correlate_present = False`, the user may want to broaden the mutation filter, narrow the ligand blocklist, or accept that the clusters are driven by features the annotation pipeline does not capture.

### Stage 7 — Iterate

- **User action.** Adjust `Config` parameters. Re-run.
- **System action.** Same as above.
- **Typical iterations.** Adjusting `cluster_distance_cut` after seeing the dendrogram; adjusting `ligand_blocklist` after seeing what was dropped; adjusting `conservation_threshold` to interrogate the strict-vs-majority distinction; adjusting `mutation_type_filter` to include `"Conflict"` if the user wants to see strain variants.

---

## Component / Module Breakdown

The notebook orchestrates and visualises. Logic lives in the `pdbe_interfaces/` package so it is testable and reusable. Each notebook cell is short — calls into the package, displays a result.

### Project structure

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

### Module: `pdbe_interfaces.api`

- **Responsibility.** All HTTP. Parses JSON. No domain logic.
- **Inputs.** Complex IDs, PDB IDs, chain/residue identifiers (per endpoint).
- **Outputs.** Parsed dicts/lists.
- **Dependencies.** `requests`, `logging`.
- **Functions.**
  ```python
  def fetch_complex_details(complex_id: str) -> dict
  def fetch_interface_interactions(complex_id: str) -> dict
  def fetch_mutations(pdb_ids: list[str]) -> dict
  def fetch_modifications(pdb_ids: list[str]) -> dict
  def fetch_bound_molecules(pdb_id: str) -> list[dict]
  def fetch_ligand_interactions(pdb_id: str, chain_id: str,
                                author_residue_number: int) -> list[dict]
  ```
- **Implementation notes.** No retries in v1. 404 on `fetch_ligand_interactions` returns `[]` and logs a warning; other 4xx/5xx raise.
- **Integration points.** Called from notebook Phase 1 and from `annotations.py`.

### Module: `pdbe_interfaces.representation`

- **Responsibility.** Build the `InterfaceRecord` dataclass per interface, including the two interaction-pair sets and the residue→UniProt mapping. Run partner-consistency check.
- **Inputs.** The `interface_interactions` JSON.
- **Outputs.** `list[InterfaceRecord]`.
- **Dependencies.** Standard library only.
- **Functions.**
  ```python
  def build_interface_records(interface_response: dict) -> list[InterfaceRecord]
  def build_interaction_pair_sets(record: InterfaceRecord) -> tuple[set, set]
  def check_partner_consistency(records: list[InterfaceRecord]) -> list[InterfaceRecord]
  ```
- **`InterfaceRecord` dataclass.** Holds `pdb_id`, `assembly_id`, `interface_id`, the `interface_info` PISA metrics (as a nested dict), the two interaction-pair sets, the residue→UniProt mapping, the partner accessions and roles. No methods beyond construction.
- **Integration points.** Outputs feed `similarity.py`, `annotations.py`, and `outputs.py`.

### Module: `pdbe_interfaces.similarity`

- **Responsibility.** Jaccard matrix and clustering.
- **Inputs.** A list of sets (one per interface).
- **Outputs.** `numpy.ndarray` similarity matrix; `ClusterResult` dataclass holding linkage and flat assignment.
- **Dependencies.** `numpy`, `scipy.cluster.hierarchy`.
- **Functions.**
  ```python
  def jaccard_similarity_matrix(sets: list[set]) -> np.ndarray
  def cluster_interfaces(distance_matrix: np.ndarray,
                         method: str = "average") -> ClusterResult
  ```
- **Implementation notes.** Symmetric matrix construction; diagonal forced to 1.0 for floating-point cleanliness.
- **Integration points.** Called from notebook Phase 3.

### Module: `pdbe_interfaces.annotations`

- **Responsibility.** Mutation filtering, modification ingestion, ligand workflow (filter → fetch interactions → join). Annotation→interface joins.
- **Inputs.** Mutation/modification responses, bound-molecule lists, interface records.
- **Outputs.** `AnnotationOverlap` dataclass holding the three long-form dataframes.
- **Dependencies.** `pandas`, `pdbe_interfaces.api`, `logging`.
- **Functions.**
  ```python
  def filter_bound_molecules(molecules: list[dict],
                             blocklist: set[str]) -> list[dict]
  def overlap_annotations(records: list[InterfaceRecord],
                          mutations: dict, modifications: dict,
                          ligand_contacts: dict) -> AnnotationOverlap
  ```
- **Integration points.** Called from notebook Phase 4.

### Module: `pdbe_interfaces.outputs`

- **Responsibility.** Build user-facing artefacts.
- **Inputs.** Interface records, similarity result, cluster result, annotation overlap.
- **Outputs.** Pandas DataFrames.
- **Dependencies.** `pandas`.
- **Functions.**
  ```python
  def build_structure_table(records, cluster_result, overlap, fetch_date) -> pd.DataFrame
  def conserved_residues(records: list[InterfaceRecord], threshold: float = 1.0) -> set
  def conserved_interaction_pairs(records: list[InterfaceRecord], threshold: float = 1.0) -> set
  def interface_frequency_summary(records, partner_map=None) -> dict   # for Phase 5a
  def cluster_interpretation_report(
      records, cluster_result, overlap,
      assembly_metadata=None, partner_map=None,
      min_count=2, min_enrichment=2.0, top_n_contacts=10,
      conservation_threshold=0.8,
      sparse_pair_threshold=5,
      tiny_area_threshold_a2=500.0,
      range_overlap_min=0.2,
  ) -> pd.DataFrame
  def compare_clusters(records, cluster_result, cluster_a, cluster_b, ...) -> pd.DataFrame
  def compare_cluster_contacts(
      records, cluster_result,
      cluster_a=None, cluster_b=None,
      typed=True, top_n=20,
      partner_map=None,
      a_enriched_threshold=0.5, shared_core_threshold=0.8,
  ) -> pd.DataFrame
  def export_interface_frequency_json(
      records,
      partner_map=None,
      complex_id=None, complex_name=None, oligomeric_state=None,
      config=None,
      typed_contacts=True,
      complex_details=None,
  ) -> dict
  ```
- **Integration points.** Called from notebook Phase 5 (cluster interpretation report, conserved sets), Phase 5a (interactive frequency summary), Phase 5b (interface rewiring), and Phase 7 (`export_interface_frequency_json` writes the per-residue and per-contact frequency JSON for downstream visualisation). `cluster_interpretation_report` carries the per-state statistics, core contacts, and QC warnings; `compare_cluster_contacts` produces the rewiring table; the JSON export reads only `records` (plus optional `partner_map` and `complex_details`) and is independent of clustering.

### Notebook (`notebook.ipynb`)

- **Responsibility.** Five-phase narrative. Configuration. Visualisation. Markdown commentary explaining each phase.
- **Implementation notes.** Each cell is short — calls into the package, displays a result. Markdown cells precede each code cell with a one-sentence explanation aimed at the user audience (see User-Facing Explanations).

---

## User-Facing Explanations and Copy

The notebook embeds the following copy verbatim. All text is for a technically literate audience (structural biologists, medicinal chemists, bioinformaticians). Tone is direct and avoids hedging.

### Title cell

> **Aggregated Interface Interaction Analysis**
>
> Compare protein–protein interfaces across all deposited structures of one complex. Cluster the interfaces by interaction-pair similarity, then overlay mutations, modifications, and ligands to interpret what each cluster represents.
>
> Working example: **PDB-CPX-130306 (KRAS–RAF1 heterodimer).**

### Configuration cell preamble

> Edit the `Config` instance below. The defaults are tuned for the working example. Common adjustments: change `complex_id` to your target complex; broaden `mutation_type_filter` to include `"Conflict"` if you want to see strain variants; remove a chem_comp_id from `ligand_blocklist` to retain a catalytic metal or buffer component.

### Phase 1 preamble

> **Phase 1 — Data retrieval.** We validate the complex ID, confirm it is a dimer (v1 limitation, set by the `interface_interactions` API), and fetch interfaces, mutations, modifications, and ligands. Annotations covering all PDB entries are batched. Per-entry calls are sequential. Expect 30–120 seconds for a complex like the working example.

### Phase 2 preamble

> **Phase 2 — Build representations.** Each interface is represented as a set of `(residue_1, residue_2, bond_type)` tuples. We build two versions: one keyed by author chain and residue number (used for joining annotations within a structure), one keyed by UniProt accession and sequence position with a role label (used for comparison across structures). Interface residues lacking a UniProt mapping are dropped — these are usually termini or expression-tag residues, and they are not comparable across structures anyway.

### Phase 3 preamble

> **Phase 3 — Similarity and clustering.** Jaccard similarity between two interfaces is the fraction of interaction-pair tuples they share. We cluster by hierarchical clustering with average linkage on `1 − Jaccard` distance. The dendrogram is the source of truth — adjust `cluster_distance_cut` after looking at it.

### Phase 4 preamble

> **Phase 4 — Annotation overlap.** Mutations, modifications, and ligand-contact residues are joined onto interface residues. Annotations falling outside the interface are dropped: this workflow is about what happens *at* the interface. A residue carrying a peripheral mutation will not show that mutation in its overlap row; the structure table reflects interface counts only.

### Phase 5 preamble

> **Phase 5 — Summarisation and interpretation.** The structure table consolidates one row per interface. Conserved residues and conserved interaction pairs are those present in at least `conservation_threshold` of interfaces (default 0.8). The cluster interpretation report reads each cluster as a candidate **interface interaction state** and surfaces, for each: per-interface distributions of residue-pair count, typed-interaction count, and PISA interface area; median contact density; **core contacts** present in ≥ `conservation_threshold` of cluster members; enriched contacts and annotations vs the rest of the dataset; and **QC warnings** for singleton clusters, sparse fingerprints, tiny interfaces, and clusters whose UniProt residue range poorly overlaps the dominant cluster. Clusters with no obvious correlate are reported as such, not speculated about.

### Phase 5b preamble — interface rewiring

> **Phase 5b — Interface rewiring between interaction states.** A head-to-head comparison of two states. Each typed contact is labelled `shared core` (≥ 80% in both states), `A-enriched` / `B-enriched` (fraction differs by ≥ 0.5), or `rare`. Defaults to the two largest non-singleton states; pass `cluster_a` / `cluster_b` to compare any pair from the report above.

### Tooltip / inline help — conservation threshold

> Default `0.8` tolerates 1–2 drop-outs at typical sample sizes (N≈10–15). Set to `1.0` to require strict conservation across every interface; set to `0.5` to surface "majority-present" residues and pairs rather than conserved ones.

### Tooltip / inline help — cluster distance cut

> The cut on `1 − Jaccard` distance that defines flat cluster assignments. Default `0.5` is a reasonable starting point for the working example. Read the dendrogram and set the cut where the visible structure of the data sits.

### Warning — residues dropped for missing UniProt mapping

> Dropped {n} residues at interface {pdb_id}/{assembly_id}/{interface_id} for missing UniProt mapping. These residues are not included in cross-structure comparison. Total drops are recorded in the structure table.

### Warning — partner-order inconsistency

> Partner ordering inconsistent across interfaces: `(unp_accession_1, unp_accession_2) = ({a}, {b})` for the majority, but `({b}, {a})` for {pdb_id}/{assembly_id}/{interface_id}. Reversing role assignment for the inconsistent entries before aggregation.

### Warning — ligand filter applied

> Filtered out {n} bound molecules across {m} PDB entries: {dropped_chem_comp_ids}. Adjust `ligand_blocklist` or `drop_carbohydrate_polymers` in `Config` to retain.

### Empty-state — no mutations at any interface

> No engineered mutations found at any interface for this complex. The mutation overlap table is empty; the structure table mutation count is zero everywhere. This is a valid result, not an error — it means none of the deposited structures carry an engineered mutation that lands on the interface (mutations elsewhere in the protein are dropped at the join step).

### Empty-state — no ligands at any interface

> No ligand contacts at any interface for this complex (after applying `ligand_blocklist` and the carbohydrate-polymer filter). The ligand overlap table is empty.

### Cluster report — annotation correlate present

> Cluster {cluster_id} ({n} interfaces): all members carry {feature}. {member_pdb_ids}.

### Cluster report — no annotation correlate

> Cluster {cluster_id} ({n} interfaces): no distinguishing mutation, modification, or ligand identified for this cluster's members. The cluster may reflect conformational state, crystal form, construct boundaries, refinement protocol, or an artefact of the similarity metric. Inspect the member structures directly to interpret further.

### Error — complex did not resolve

> Complex ID `{complex_id}` did not resolve. Check the identifier is a valid PDBe-KB complex ID (format `PDB-CPX-NNNNNN`).

### Error — non-dimer complex

> Complex `{complex_id}` has {n} components. v1 supports dimer complexes only — this is an upstream limitation of the `interface_interactions` API, not a design choice. Higher-order oligomers are a planned v2 extension.

### Error — no interfaces returned

> No interfaces returned for `{complex_id}`. The complex ID is valid but the interface API has no contacts on file. This is unexpected for a deposited dimer complex; report the case if it persists across reruns.

### Trust / safety

> This workflow is descriptive at v1 sample sizes (~10–30 interfaces per complex). It does not compute statistical enrichment tests, and claims about cluster–annotation correspondence should be read as observations, not inferences. Annotations not captured by the PDBe pipelines (conformational state, crystal form, refinement protocol) are invisible to this analysis; absence of annotation overlap does not imply absence of effect.

---

## Data, Logging, and Observability

### Data persistence

v1 persists nothing to disk by default. All artefacts live in the notebook's output cells and the in-memory dataframes. A user wanting to export takes a manual step: `df.to_csv(...)`, `plt.savefig(...)`. This is intentional — v1 is exploratory; persistence patterns are deferred until they have been used.

The structure table records `fetch_date` (ISO-8601 date) so that any exported CSV remains interpretable across PDBe releases.

### Logging

Use the standard library `logging` module, configured at the top of the notebook with `Config.log_level` (default `INFO`).

**INFO-level events.**
- Each API call: endpoint, inputs, response status, parsed-record count.
- Filter applications: count of dropped bound molecules per PDB entry, with chem_comp_ids.
- Phase boundaries: `"Phase 1 complete: {n} interfaces, {m} PDB entries"`.

**WARNING-level events.**
- Residues dropped for missing UniProt mapping (count and interface key).
- Partner-order inconsistency (the reversed entries).
- 404 on `bound_ligand_interactions` (the ligand instance skipped).

**ERROR-level events.**
- Any HTTP or JSON parsing failure that is about to be raised. Logged before raise so the message survives in notebook output even if the exception traceback is collapsed.

Logging configuration:
```python
logging.basicConfig(
    level=Config.log_level,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
```

### Metrics (in-notebook)

Surfaced as printed summaries at phase boundaries, not collected to an external store:
- Number of interfaces fetched.
- Number of PDB entries.
- Number of mutations / modifications / ligand contacts before and after filtering.
- Number of interface residues dropped for missing UniProt mapping.
- Number of partner-order corrections applied.
- Distribution of Jaccard values (min/median/max).
- Number of clusters at the default cut.

### Telemetry

None. v1 is a local notebook; no telemetry endpoint is wired up.

### Debugging

Each module function is testable in isolation. The `InterfaceRecord` dataclass is the central debugging surface — a list of `InterfaceRecord` objects is sufficient to reproduce any downstream computation. Encourage users debugging unexpected output to inspect `records[i]` directly.

### Auditability

The combination of `fetch_date` in the structure table and the inline log of filter applications gives a sufficient audit trail at v1: any export from the notebook is dated, and the filter list is recorded in `Config`.

### Privacy / security

PDBe data is public. No user credentials are required. No personally identifiable information is handled. The notebook makes outbound HTTPS requests to `www.ebi.ac.uk` only.

---

## Integration Points

### Upstream dependencies

- **PDBe v2 API** — `complex/details`, `complex/interface_interactions`, `pdb/entry/mutated_AA_or_NA`, `pdb/entry/modified_AA_or_NA`, `pdb/bound_molecules/{pdb_id}`, `pdb/bound_ligand_interactions/{pdb_id}/{chain_id}/{author_residue_number}`. The workflow trusts these endpoints' schemas as documented.
- **PDBe-KB complex catalogue** — the source of `PDB-CPX-NNNNNN` identifiers. Users obtain complex IDs by browsing PDBe-KB; the notebook does not include a search interface.
- **PISA pipeline (upstream of `interface_interactions`)** — the source of the contact-type vocabulary (`hydrogen_bond`, `salt_bridge`, etc.) and the interface-level metrics.
- **SIFTS (upstream of the inline UniProt mappings)** — the source of the residue-level UniProt mappings inlined in the interface response. The workflow relies on SIFTS providing a one-to-one mapping per `(pdb_id, auth_asym_id, auth_seq_id, ins_code)`.

### Downstream consumers

None defined for v1. The natural consumers are:
- **Other notebooks or scripts** that want the structure table, conserved interaction pairs, or cluster assignments for a complex. Manual export via `df.to_csv` or pickling of the `InterfaceRecord` list.
- **Reference dataset curation pipelines** that aggregate conserved interaction pairs across many complexes. The interaction-pair tuple format is stable and reusable.

### Export formats

Manual on demand:
- Structure table → CSV (`pd.DataFrame.to_csv`).
- Annotation overlap dataframes → CSV.
- Conserved sets → JSON or pickle.
- Plots → PNG via `plt.savefig`.

No automated export in v1.

### Workflow integrations and automation hooks

None in v1. The notebook is run interactively. Headless execution via `papermill` would work in principle (the `Config` dataclass is parameter-friendly), but is out of scope.

### Reusability conditions and promotion to production

The `pdbe_interfaces/` package is the reusable surface. Its API is the function signatures listed in Component Breakdown. Promotion to a production pipeline would require:
- Pinning dependency versions.
- Adding API response caching to disk (out of scope for v1).
- Adding retries with exponential backoff on HTTP failures (out of scope for v1).
- Test coverage on the representation, similarity, and annotation-join layers.
- Explicit handling of the assumptions in §11 of the source spec (currently load-bearing and unchecked).

### Interoperability expectations

The interaction-pair tuple format is designed to be source-agnostic. A future v3 extension that compares predicted complexes (AlphaFold-Multimer, Boltz, Chai) against experimental structures replaces `pdbe_interfaces.api` with a local interaction-pair extractor that produces the same `InterfaceRecord` shape. Downstream similarity, clustering, and annotation overlap machinery is unchanged.

---

## Non-Functional Requirements

### Performance

- **End-to-end runtime on the working example (PDB-CPX-130306).** Target: ≤120 seconds. Dominated by per-entry `bound_molecules` and per-ligand `bound_ligand_interactions` GETs.
- **Jaccard matrix computation.** Negligible at v1 scale (≤30 interfaces). No optimisation target.
- **Clustering and dendrogram rendering.** Negligible.

### Scalability

- **Per complex.** v1 targets ≤30 interfaces per complex. Behaviour at higher counts (e.g. 100+) is not tested but should remain correct; runtime scales linearly with the number of bound-molecule fetches.
- **Across complexes.** v1 handles one complex per notebook run. Multi-complex aggregation is a v2 extension.

### Reliability

- API failures bail fast with a clear message. No silent partial-data paths.
- Cross-structure consistency check (partner ordering) actively corrects an inconsistency rather than ignoring it.
- Documented load-bearing assumptions (§11 of the source spec) are listed so they are auditable, even though they are not actively checked.

### Latency

Not user-interactive in the sub-second sense; the notebook is exploratory. Per-cell latency targets:
- Phase 1 (data retrieval): up to 120 s.
- Phases 2–5: <5 s each.

### Maintainability

- Logic is in `pdbe_interfaces/`, not in the notebook, so it is testable.
- The `Config` dataclass centralises every parameter the spec marks as user-overridable.
- The `InterfaceRecord` dataclass is the canonical intermediate representation; downstream modules consume it.

### Portability

- Python 3.11+ on any OS supported by `numpy`/`scipy`/`pandas`.
- No native dependencies beyond what `numpy`/`scipy` carry.

### Usability

- Five-phase notebook narrative with markdown preambles for each phase.
- All parameters in one `Config` dataclass at the top of the notebook; no hunting through cells.
- Visualisations have sensible defaults; no per-plot tuning required for the working example.

### Accessibility

- Colour-blind-safe colormap (`viridis`) on the heatmap.
- Cluster interpretation report is text-readable.
- 90° axis-label rotation on the heatmap to avoid overlap.

### Security

- No credentials handled; all PDBe v2 endpoints are public.
- All outbound traffic is to `www.ebi.ac.uk` over HTTPS.
- No user-supplied code execution paths.

### Cost constraints

- Zero external service cost (PDBe is free).
- Local compute only.

---

## Out of Scope for v1

- **Higher-order oligomers (trimers, tetramers, larger assemblies).** Rationale: the `interface_interactions` API currently returns data for dimers only. A v2 extension when the API supports them.
- **Multiple complexes per notebook run.** Rationale: v1 establishes the per-complex workflow; cross-complex aggregation is a v2 extension on the same axis.
- **Typed ligand contacts feeding into similarity.** Rationale: changes the conceptual unit from "protein–protein interface" to "extended interaction surface." Useful for drug discovery but a v2 design choice.
- **Statistical enrichment tests for cluster–annotation correspondence.** Rationale: not meaningful at v1 sample sizes (~10–30 interfaces per complex). Becomes useful above ~30, which v1 does not assume.
- **Predicted complexes (AlphaFold-Multimer, Boltz, Chai) as input.** Rationale: requires a local interaction-pair extractor since PDBe APIs do not run on user-supplied models. The interaction-pair abstraction supports this; v3 extension.
- **Machine learning approaches to similarity, clustering, or annotation prediction.** Rationale: Jaccard + hierarchical clustering is appropriate for the sample sizes, is interpretable, and avoids the data-volume and validation problems of model-based methods at this scale.
- **Sequence-based conservation analysis.** Rationale: this workflow uses structural instances of the same complex, not sequence homologues. Cross-species or cross-paralogue conservation is a separate problem.
- **Structure prediction or docking.** Rationale: this is an analysis tool, not a generation tool.
- **Per-atom contact analysis.** Rationale: the aggregation rule collapses atom-level detail; specific donor–acceptor geometry is intentionally outside the workflow's scope.
- **API response caching to disk.** Rationale: v1 fetches fresh; caching is a build optimisation deferred until concrete need.
- **Retry / backoff on API failures.** Rationale: v1 bails fast; retries are a build optimisation deferred until concrete need.
- **Async / parallel API calls.** Rationale: v1 call volume does not justify async machinery.
- **Tests beyond the acceptance checks.** Rationale: v1 is exploratory; comprehensive tests are deferred until the API surface stabilises.
- **CI, packaging, distribution.** Rationale: v1 is a single-user notebook; no distribution surface.
- **Headless / `papermill` execution.** Rationale: parameters are dataclass-friendly so this would work, but it is not a v1 deliverable.
- **Automated export of artefacts.** Rationale: users export manually as needed; standard patterns are not yet established.

### Acceptance checks for v1 (definition of done)

When run on the working example (`PDB-CPX-130306`), the notebook produces:
- A structure table with at least 10 rows (~12 PDB entries × assemblies).
- A similarity matrix that is symmetric with 1.0 on the diagonal.
- A dendrogram that visually separates at least two clusters at the default cut.
- A non-empty ligand overlap. KRAS structures carry GNP (a GTP analogue) at the GTP-binding pocket; whether GNP also reaches the RAF1 interface depends on the entry, but some entries do.
- No unhandled exceptions.

**Note on mutation overlap.** An earlier version of this spec listed "non-empty mutation overlap" as an acceptance check on the assumption that KRAS oncogenic mutations would land at the interface. They do not — G12, G13, Q61 and C118 sit at the GTP-binding site and the activation loop, not at switch I (which is where RAF1 binds). The workflow correctly returns zero at-interface mutations for this complex. For a complex where interface mutations *are* expected, run on `PDB-CPX-140195` (ACE2–Spike) — RBD variant residues like S477N, N501Y, Q493R, Q498R land at the interface and are surfaced as enriched signatures in the cluster interpretation report.

If the workflow runs end-to-end on `PDB-CPX-130306` and produces these outputs, v1 is functionally complete. Quality of the analysis (whether the clusters are interpretable) is a separate question and is the user's call.

---

## Open Questions Requiring Resolution Before Build

### 1. Microheterogeneity in SIFTS mappings

- **Question.** Are there real PDB entries in the v1 working example (or near neighbours) where two distinct author residues map to the same UniProt position?
- **Why it matters.** The workflow assumes a one-to-one mapping per `(pdb_id, auth_asym_id, auth_seq_id, ins_code)`. If violated, distinct residues silently collapse in the UniProt-keyed set, producing wrong Jaccard values without any visible failure.
- **Proposed next step.** Add a one-shot check on the inline UniProt mappings during Phase 2: for each PDB entry, verify that no two distinct author keys map to the same `(unp_acc, unp_seq_id)`. If violated, log at WARNING and (open question) decide whether to drop the affected residues, keep one, or abort. Recommend: log and keep the first occurrence; revisit if observed in practice.

### 2. PISA partner-ordering inconsistency: directionality

- **Question.** When `(unp_accession_1, unp_accession_2)` ordering is inconsistent across interfaces, the workflow reverses role assignment for the *minority* ordering. Is the majority always the correct canonical direction?
- **Why it matters.** Picking the wrong canonical direction does not break correctness (Jaccard is invariant to a global swap), but it does affect the partner labels in plots and the cluster interpretation report. A user expecting role-1 to be KRAS will see the opposite if PDBe's majority ordering is RAF1-first.
- **Proposed next step.** Add a `Config.canonical_partner_accession: str | None = None`. If set, force role-1 assignment to that accession regardless of majority. If unset, use the majority. Document the override in the configuration cell preamble.

### 3. Interface-level metric semantics (PISA p-value)

- **Question.** The `interface_interactions` response includes a PISA p-value among the interface-level metrics. The structure table surfaces it. What does the user do with it?
- **Why it matters.** PISA's p-value is a measure of biological-vs-crystallographic interface significance. At v1 we surface it without filtering. If a user excludes interfaces above a p-value threshold, the cluster structure changes. We currently expose no such filter.
- **Proposed next step.** Add `Config.pisa_pvalue_max: float | None = None`. If set, filter interfaces above the threshold before clustering and record the filter application in the structure table. Default `None` (no filter) preserves v1 behaviour. Decide before build whether to expose the filter or keep the surface lean.

### 4. Untyped Jaccard variant: report-only or alternative

- **Question.** The spec calls for an untyped Jaccard variant as a robustness check. Is it report-only, or should the user be able to switch the primary clustering to use it?
- **Why it matters.** The bond-type vocabulary varies in granularity across PISA versions; a user wanting cross-release stability might prefer untyped clustering. The spec presents the typed version as primary.
- **Proposed next step.** Keep typed primary (per spec). Render both heatmaps. Print a summary of where they differ. Defer the "switch primary" knob until a user asks for it.

### 5. Cluster-distance-cut default

- **Question.** Is `0.5` (on `1 − Jaccard`) the right default for the working example? Is it sensible across the broader set of dimer complexes in PDBe?
- **Why it matters.** The default determines what users see first. A bad default produces cluster reports that look wrong before the user adjusts the cut.
- **Proposed next step.** During build, run the workflow on the working example and 2–3 other dimer complexes (e.g. an antibody–antigen pair, an enzyme–substrate pair) and confirm that `0.5` produces a sensible cluster count (>1, <N/2). Adjust if not. Document the empirical basis in the `Config` docstring.

### 6. Partner labelling for homodimers in the structure table

- **Question.** For homodimers, both partners share an accession and (often) the same gene name. The structure table's `partner_1_label` and `partner_2_label` will be identical. Should we append a chain identifier to disambiguate?
- **Why it matters.** Identical labels in both partner columns are visually unhelpful and risk confusion in plot axes.
- **Proposed next step.** For homodimers, append `(role 1)` and `(role 2)` to the gene name in the structure table and in plot labels. For heterodimers, leave gene names as-is. Spec already notes the user can post-process; this is the in-notebook treatment.

---
