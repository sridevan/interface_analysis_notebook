# Aggregated Interface Interaction Analysis Notebook - Spec

## 1. One line idea
A Jupyter notebook to compare protein–protein interfaces across multiple structural instances of the same complex, using interaction-pair similarity, clustering, and annotation enrichment.

---

## 2. Working example
Example complex:
Complex ID: PDB-CPX-130306 
Complex name: KRAS–RAF1 heterodimer

This example is used for v1 development because:
- clean protein–protein heterodimer
- multiple PDB structures
- contains mutations and ligands
- suitable for interface comparison

---

## 3. Why
Most interface analyses are reported per-structure, even when a complex has many deposited structures. This workflow treats a complex as a dataset of structural instances — varying in mutations, bound ligands, and experimental conditions — and supports:
- comparing interfaces across structures
- identifying structurally conserved versus variable interactions
- relating structural variation to mutations, modifications, and ligands

Useful for structural biologists analysing a known complex, medicinal chemists comparing apo and ligand-bound states, and bioinformaticians building reference datasets. The interaction-pair abstraction is structure-source agnostic, so the same approach extends naturally to comparing predicted complexes (e.g. AlphaFold-Multimer) against experimental structures, though this is out of scope for v1 (see §12). Not a prediction or docking tool, and not a sequence-conservation analysis.

---

## 4. Input

### Primary input
A PDBe-KB complex ID (e.g. `PDB-CPX-130306`).

### Scope
v1 supports **dimer complexes only**, because the `interface_interactions` endpoint currently returns data for dimers only (see §12). The notebook validates the input by calling `complex/details` first and aborts with a clear message if the complex ID does not resolve, or if the complex is not a dimer.

### Endpoints used

| Data | Endpoint | Method |
|---|---|---|
| Complex membership and partner identities | `https://www.ebi.ac.uk/pdbe/api/v2/complex/details/{complex_id}?id_type=pdb_complex_id` | GET |
| Interface interactions (PISA contacts with UniProt mappings and interface-level metrics) | `https://www.ebi.ac.uk/pdbe/api/v2/complex/interface_interactions/{complex_id}` | GET |
| Mutations | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/mutated_AA_or_NA` | POST (batched, JSON-encoded comma-separated string of PDB IDs) |
| Modifications | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/modified_AA_or_NA` | POST (batched, JSON-encoded comma-separated string of PDB IDs) |
| Bound molecules | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_molecules/{pdb_id}` | GET |
| Ligand interactions | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/bound_ligand_interactions/{pdb_id}/{chain_id}/{author_residue_number}?preserve_case=false` | GET (one per ligand instance) |

### Resolution filtering (optional)

`Config.max_resolution: float | None = None` is an optional cap on assembly resolution from `complex/details`. When set (e.g. `2.0`), assemblies are dropped if their resolution exceeds the threshold or is `None`. The strict interpretation matches typical user intent when filtering: keep only structures with explicit, sufficient resolution data. Default `None` disables the filter.

### Assembly filtering

Interfaces from assemblies with bound macromolecules are excluded from the analysis. The `complex/details` response provides a `bound_macromolecules` list per assembly; assemblies whose list is non-empty are dropped, along with every interface they contribute.

The rationale is that a bound macromolecule (a co-chaperone, an antibody Fab, a third protein chain not part of the canonical dimer) can perturb the dimer interface — either directly by contacting interface residues, or allosterically by stabilising a non-native conformation. Excluding these assemblies keeps the cross-structure comparison clean: every contributing interface is the canonical dimer in a state without third-party macromolecular perturbation.

Bound *small-molecule* ligands (the subject of §8) are handled separately and are not the subject of this filter — many complexes have catalytic cofactors or substrate analogues bound by design, and excluding those would discard most of the dataset.

The filter is implemented in Phase 1: a set of valid `(pdb_id, assembly_id)` tuples is built from `complex/details`, and interfaces from `interface_interactions` are filtered against it before any downstream processing. Dropped interfaces are logged at WARNING with the count.

### Notes
- The `interface_interactions` response provides per-residue UniProt mappings inline (`unp_accession`, `unp_seq_id`) alongside author numbering, so no separate SIFTS call is required. UniProt mappings are available for interface residues; annotations falling outside the interface are not mapped to UniProt in this workflow (see §8).
- The `interface_interactions` response also provides interface-level PISA metrics (interface area, solvation energy, stabilisation energy, p-value, hydrogen bond / salt bridge / covalent / disulfide bond counts, interface residue count) which are surfaced in §9.
- API responses are fetched fresh on each notebook run; the structure table records the fetch date so outputs remain interpretable across PDBe releases.
- API calls retry up to 5 times with exponential backoff (1 s / 2 s / 4 s / 8 s / 16 s) on transient errors (HTTP 429 / 500 / 502 / 503 / 504, `ConnectionError`, `Timeout`). Critical calls (`complex/details`, `interface_interactions`, batched mutations / modifications) raise after the budget is exhausted. Per-entry / per-ligand best-effort calls (`bound_molecules`, `bound_ligand_interactions`) return an empty result with a logged WARNING after the budget is exhausted, so a single failed call does not kill an analysis with hundreds of such fetches. The two batch-POST endpoints (`mutated_AA_or_NA`, `modified_AA_or_NA`) chunk the PDB ID list into blocks of 50 per request and merge results.
- The unit of analysis is the **interface**, keyed by `(pdb_id, assembly_id, interface_id)` from the `interface_interactions` response. (The response itself uses the field name `entry_id` for the PDB identifier; the notebook renames this to `pdb_id` on ingest for consistency with the rest of the workflow.) Every interface returned by the API is treated as an independent observation. No filtering by assembly is performed — `assembly_id` is recorded as metadata on each row of the structure table. PDB entries with multiple assemblies or multiple interfaces per assembly contribute multiple rows.

---

## 5. Residue numbering strategy

Three numbering systems coexist in the data:

- **Author numbering** — chain and residue labels assigned by depositors, including an optional insertion code. Used by both the interface API and all annotation APIs (mutations, modifications, bound molecules, ligand interactions).
- **UniProt numbering** — canonical sequence position. Provided inline by the interface API for interface residues; not provided by any annotation API.
- **Label numbering** (SEQRES position) — provided by some endpoints but not used in this workflow.

### Two keys

**Author key** (within-structure join layer):
`(pdb_id, auth_asym_id, auth_seq_id, ins_code)`

Used to join annotations onto interface residues within a single PDB entry. Both sides have author numbering, so the join is direct.

**UniProt key** (cross-structure comparison layer):
`(unp_acc, unp_seq_id, role)`

Used to compare equivalent residues across different structures. The `role` field (1 or 2) is taken from the PISA interface response — the partner under `unp_accession_1` is role 1, the partner under `unp_accession_2` is role 2. `role` disambiguates the two sides of homodimer interfaces, where accession alone is insufficient. UniProt mappings are available only for interface residues (provided inline by the interface API); annotations falling outside the interface are dropped at the overlap step (see §8).

### Field-name normalisation across endpoints

Different endpoints use different field names for what is logically the same author key. The notebook normalises them to a common shape on ingest:

| Source | Author chain field | Author residue field | Insertion code field |
|---|---|---|---|
| Interface (PISA) | `auth_asym_id_1` / `auth_asym_id_2` | `auth_seq_id_1` / `auth_seq_id_2` | (when present) |
| Mutations | `chain_id` | `author_residue_number` | `author_insertion_code` |
| Modifications | `chain_id` | `author_residue_number` | `author_insertion_code` |
| Bound molecules | `chain_id` | `author_residue_number` | `author_insertion_code` |
| Ligand interactions (`end` block) | `chain_id` | `author_residue_number` | `author_insertion_code` |

Insertion-code values vary across endpoints (`""` in mutations and modifications, `" "` in bound molecules and ligand interactions). The notebook treats any whitespace-only value as "no insertion code" via `s.strip() or None`.

### Order of operations

1. Within each interface response, build the author-keyed and UniProt-keyed interaction-pair sets (per §6).
2. Within each PDB entry, join annotations to interface residues using the author key.
3. After joining, read off the UniProt key from the interface side for cross-structure analysis.
4. Interface residues that lack a UniProt mapping in the API response are dropped from the analysis with a logged warning, and the count is recorded per interface in the structure table.

### Notes
- Insertion codes are retained in the author key for correctness. The workflow trusts SIFTS to provide a one-to-one mapping per `(pdb_id, auth_asym_id, auth_seq_id, ins_code)` and does not handle the pathological case of two author residues mapping to the same UniProt position. See §11 for the full statement of this assumption.

---

## 6. Data representation

### Unit of analysis

A single **interface**, keyed by:

`(pdb_id, assembly_id, interface_id)`

Each interface in the `interface_interactions` response contributes one row to the structure table and one set of interaction pairs to the analysis. PDB entries with multiple assemblies or multiple interfaces per assembly contribute multiple rows.

### Residue identifiers

**Author-numbered residue:**
`(pdb_id, auth_asym_id, auth_seq_id, ins_code)`

Used to join annotations onto interface residues within a single structure.

**UniProt-numbered residue:**
`(unp_acc, unp_seq_id, role)`

Used to compare equivalent residues across structures. `role ∈ {1, 2}` is taken directly from the PISA response — residues under `auth_asym_id_1` / `unp_accession_1` are role 1; residues under `auth_asym_id_2` / `unp_accession_2` are role 2. This rule applies uniformly to heterodimers and homodimers.

### Interaction-pair tuples

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

`bond_type` takes values from the PISA vocabulary (e.g. `hydrogen_bond`, `salt_bridge`, `covalent_bond`, `disulfide_bond`, and others as returned by the API). Values are recorded as they appear in responses; no manual remapping is performed.

### Aggregation rule

Each row in the PISA response is an atom-level contact between two residues. Atom-level contacts sharing the same key (residue pair plus `bond_type`) collapse to a single set element. A residue pair forming contacts of different `bond_type` values appears as multiple elements. Atom-level detail is not retained.

Aggregation is performed per interface and per key:
- The author-numbered interaction set aggregates by author key.
- The UniProt-numbered interaction set aggregates by UniProt key (which collapses any author-level distinctions that map to the same UniProt position).

Aggregation keys treat interaction-pair tuples as ordered: `((residue_role_1), (residue_role_2), bond_type)`. For heterodimers this is unambiguous. For homodimers, this preserves both directions of a symmetric contact (e.g. `(chain1:25, chain2:87)` and `(chain1:87, chain2:25)`) as distinct elements, which is the correct behaviour for physical analysis.

### Cross-structure consistency check

When building the cross-structure interaction set, the notebook verifies that `(unp_accession_1, unp_accession_2)` is the same for every interface in the complex. If a different ordering appears, the notebook logs a warning and reverses role assignment for the inconsistent entries before aggregation, so that the same biological partner is always role 1 in the cross-structure set. Within a single interface, PISA is trusted to be internally consistent and no per-interface deduplication is performed.

### Partner map

A partner_map is built once per complex from the complex/details response, keyed on (unp_accession, role) and mapping to the gene name returned for each component. Used only for output labelling (structure table, residue labels, plot axes). Interaction-pair tuples themselves remain keyed on UniProt accession and role.
For homodimers, both keys share an accession; the role distinguishes them in display, and both can map to the same gene name (the user can post-process to add chain identifiers if a finer distinction is needed).

### Source-agnosticism

The interaction-pair representation does not depend on the PISA response shape. It depends only on having: a list of contacts between two residues, each residue identified by UniProt accession and position, with a contact type label. Any source that produces this — including local interaction-pair extraction from predicted complexes (e.g. AlphaFold-Multimer) — can feed into the same downstream analysis. This is the integration point referenced in §3 for the v3 extension.

---

## 7. Method

The workflow proceeds in five phases. Each phase is intended to map to one or two cells in the notebook, with each cell ending in a display of an intermediate or final result so that the notebook reads as a sequence of phases each producing an inspectable artefact.

### Phase 1 — Data retrieval

1. Validate the input complex ID by calling `complex/details`. Abort if the ID does not resolve or the complex is not a dimer. Build the `partner_map` from the response. Build a set of valid `(pdb_id, assembly_id)` tuples — those whose `bound_macromolecules` list is empty (per the assembly filter described in §4).
2. Call `interface_interactions` for the complex. The response provides per-interface PISA contacts and metrics, with UniProt mappings inline. Filter the response against the valid-assembly set; drop interfaces from assemblies with bound macromolecules and log the count.
3. Call the annotation endpoints in parallel: a single batched POST each for mutations and modifications (covering the PDB IDs from the surviving interfaces), and per-entry GETs for `bound_molecules`. For each surviving ligand (after filtering, see §8), call `bound_ligand_interactions`.

### Phase 2 — Build representations

4. For each interface in the response, build the author-numbered and UniProt-numbered interaction-pair sets per §6, including aggregation of atom-level contacts to residue-pair-and-bond-type elements.
5. Drop interface residues that lack a UniProt mapping, with a logged warning. Record the drop count per interface.
6. Run the cross-structure consistency check on `(unp_accession_1, unp_accession_2)`. Flag and reverse role assignment for any inconsistent interfaces.

### Phase 3 — Similarity and clustering

7. Compute pairwise Jaccard similarity between all interfaces, using the UniProt-numbered interaction-pair sets (typed by `bond_type`) as the primary metric. Compute an untyped variant (collapsing across `bond_type`) as a robustness check.
8. Cluster the interfaces by hierarchical clustering with average linkage on `1 − Jaccard` distance. Expose the dendrogram for visual inspection. Display a **cluster-cut sweep table** alongside the dendrogram showing cluster count and size distribution at standard cuts (0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70). Stable plateaus (same cluster count across a range of cuts) indicate a defensible default; fragmentation gradients call for closer inspection. Default cut: 0.5; the user can override based on the sweep and dendrogram.

### Phase 4 — Annotation overlap

9. For each interface, join annotations (mutations, modifications, ligand-contact residues) to interface residues using the author key. Annotations falling outside the interface are not retained.
10. After the join, attach the UniProt key of each matched interface residue from the interface side, enabling cross-structure comparison of which annotations land on which positions.

### Phase 5 — Summarisation and interpretation

11. Build the structure table (one row per interface) with PISA metrics, annotation counts, and cluster assignment.
12. Identify conserved residues and conserved interaction pairs — those present in at least `conservation_threshold` of interfaces (default 0.8); the threshold is exposed as a parameter.
13. Build the **cluster interpretation report**, treating each cluster as a candidate **interface interaction state**. For each state, report: per-interface distributions of residue-pair count, typed-tuple interaction count, and PISA interface area (min / median / mean / max), plus median contact density (interactions / Å²); **core contacts** — typed contacts present in ≥ `conservation_threshold` of cluster members; **enriched contacts** and annotations (mutations / modifications / ligands) vs the rest of the dataset; and **QC warnings** (singleton, sparse fingerprint, tiny interface, mixed UniProt residue range vs the dominant cluster). Where no correspondence is apparent, report this explicitly rather than speculate (per §10).
14. Build the **interface rewiring** table — a head-to-head comparison between two interaction states. Each typed contact is labelled `shared core` (≥ 80% in both), `A-enriched` / `B-enriched` (fraction differs by ≥ 0.5), or `rare`. Defaults to the two largest non-singleton states.

### Phase 7 — Export interface conservation JSON

15. Call `export_interface_frequency_json` to write a lightweight JSON containing UniProt residue-level interface conservation (`residue_frequencies`) and residue–residue contact frequencies (`contact_frequencies`) across the full set of retained interfaces. The export is consumed by frontend visualisation tooling (Mol* colouring, contact matrices/heatmaps, AFDB-style overlays). Filename: `{complex_id}_interface_frequencies.json`. Output directory comes from `Config.output_dir` (created if missing) or defaults to the current working directory. The export is independent of clustering — it does not modify any prior output. See §9 for the schema.

---

## 8. Annotation enrichment

Three annotation streams are joined to interface residues using the author key (per §5). Annotations falling outside the interface are dropped at the join step and not retained. Annotations are used for cluster interpretation and the structure table; they do not contribute to the similarity metric (which is computed on PISA interface contacts only — see §7).

### Mutations

A single batched POST to `mutated_AA_or_NA` covering all PDB IDs in the complex returns mutations per entry. Each mutation record carries `mutation_details.type`, `mutation_details.from`, and `mutation_details.to`.

The vocabulary for `type` includes (non-exhaustively):
- `"Engineered mutation"` — a point mutation deliberately introduced by the depositors, relative to the canonical UniProt sequence. Covers a wide range of cases: disease-associated mutations (e.g. KRAS G12C), stability or solubility mutations, catalytic knockouts, phosphomimics, and crystallisation aids.
- `"Conflict"` — a discrepancy between deposited and reference sequence, often a sequencing error or strain variant.
- `"Cloning artifact"`, `"Expression tag"` — residues from cloning or tag sequences not fully removed before crystallisation; almost always at termini.

The default filter retains `"Engineered mutation"` only, on the basis that the other categories are not deliberate experimental choices and are usually noise for an interface-comparison workflow. The filter is exposed as a parameter so the user can broaden the set or narrow it further (e.g. by `mutation_details.from`/`to` for specific substitutions of interest). The substitution is carried through as a label (e.g. `C116S`) on each matched interface residue.

### Modifications

A single batched POST to `modified_AA_or_NA` covering all PDB IDs in the complex returns modifications per entry. Each modification record carries `chem_comp_id` (the modified residue's three-letter code, e.g. `SEP` for phospho-serine, `8AN`) and `chem_comp_name`. Both are carried through to each matched interface residue. No filtering is applied by default; for protein-only complexes the response is dominated by canonical PTMs and the data is generally clean.

### Ligands

The ligand workflow is per-entry and includes a filtering step:

1. GET `bound_molecules` for each PDB ID.
2. Filter the returned bound molecules:
   - Drop entries whose ligand `chem_comp_id` is in a blocklist of buffer / cryoprotectant / counter-ion residues. Default blocklist: `CL`, `NA`, `MG`, `ZN`, `CA`, `K`, `SO4`, `GOL`, `HOH`, `EDO`, `PEG`, `MPD`, `ACT`, `SGM`, `DTT`, `BME`, `TRS`. The blocklist is exposed as a parameter so the user can override it (e.g. to retain catalytic metals).
   - Optionally drop entries with `molecule_type == "Carbohydrate-polymer"` (typically N-glycans). The parameter `drop_carbohydrate_polymers` controls this and **defaults to `False`** so glycan contacts at the interface are surfaced. For glycoprotein complexes (Spike–ACE2, antibody–Fc, lectin–glycan, viral attachment receptors) glycans at the interface are biologically meaningful — e.g. ACE2 N322 and N90 glycans participate in RBD recognition. For cytoplasmic complexes (most kinases, GTPases, transcription factors) carbohydrate polymers are not modelled at all, so the parameter has no effect. Glycan contacts only appear in the at-interface ligand overlap report; they do not influence the protein–protein similarity matrix or clustering, so changing this parameter does not fragment clusters.
3. For each surviving ligand instance, GET `bound_ligand_interactions/{pdb_id}/{chain_id}/{author_residue_number}?preserve_case=false`.
4. Each interaction record describes an atom-level contact between a ligand atom and a protein residue, with a list of contact-type labels in `interaction_details`. The list is expanded so each contact-type label produces a separate element. Atom-level detail is then collapsed: multiple atom-level records sharing the same `(interface residue, ligand instance, contact_type)` key aggregate to a single element. The output unit per ligand contact is:
   ```
   (
     (pdb_id, auth_asym_id, auth_seq_id, ins_code),
     (ligand_chem_comp_id, ligand_chain_id, ligand_author_residue_number),
     contact_type
   )
   ```
   where `contact_type` is a single value from the `interaction_details` vocabulary (e.g. `hbond`, `weak_hbond`, `polar`, `weak_polar`, `vdw_clash`, `vdw`, `hydrophobic`). Values are recorded as they appear in responses; no manual remapping is performed.
5. Join ligand-contact residues to interface residues using the author key (the first element of the tuple above). Each matched interface residue carries through the ligand instance identifier, the contact type, and the ligand's `chem_comp_id` and `chem_comp_name` for display, so cluster interpretation can identify which ligand is at which interface and via which contacts.

The PISA `bond_type` vocabulary (used for protein–protein contacts in §6) and the ligand `interaction_details` vocabulary are produced by different pipelines with different cutoffs and are not reconciled.

### Field-name normalisation

The annotation endpoints use `chain_id` and `author_residue_number`; the interface endpoint uses `auth_asym_id` and `auth_seq_id`. Both are normalised to the common author key on ingest (per §5). Insertion-code values are normalised by stripping whitespace; any whitespace-only value is treated as "no insertion code."

### Output of this phase

For each interface, three lists of matched annotations (mutations, modifications, ligand contacts), each annotation tagged with its source-specific metadata (substitution label, modification chem_comp, or ligand identity and contact type) and the UniProt key of the matched interface residue. These feed directly into the structure table, the conserved-residue report, and cluster interpretation.

---

## 9. Outputs

The notebook produces the following outputs, each available as a dataframe or visualisation:

### Structure table

One row per interface, keyed on `(pdb_id, assembly_id, interface_id)`. Columns:

- **Identification**: PDB ID, assembly ID, interface ID, partner labels (from `partner_map`).
- **Methodology**: experimental method (X-ray diffraction, Electron microscopy, etc.), resolution in Å. Sourced from `complex/details` per assembly. Surfacing these alongside biological annotations lets the user disentangle methodological clusters (e.g. "all members are 3.5 Å cryo-EM") from biological ones (e.g. "all members are Omicron").
- **PISA metrics**: interface area, solvation energy, stabilisation energy, p-value, number of interface residues, number of hydrogen bonds, number of salt bridges, number of covalent bonds, number of disulfide bonds.
- **Annotation counts**: number of engineered mutations at the interface, number of modifications at the interface, number of distinct ligands contacting the interface, number of ligand-contact residues at the interface.
- **Cluster assignment**: integer cluster ID at the default cut.
- **Quality**: count of interface residues dropped for missing UniProt mapping, fetch date.

### Similarity matrix

A square matrix of pairwise Jaccard similarities between all interfaces, computed on the UniProt-numbered, typed interaction-pair sets. An untyped variant is computed alongside as a robustness check. Rendered as a heatmap.

### Dendrogram and cluster assignments

Hierarchical clustering with average linkage on `1 − Jaccard` distance. The dendrogram is the primary visual; cluster assignments at the default cut are written into the structure table.

### Conserved residues and conserved interaction pairs

The set of UniProt-numbered residues, and the set of UniProt-numbered interaction-pair elements, present in at least a chosen fraction of interfaces. Default threshold: 0.8 (present in at least 80% of interfaces). The threshold is exposed as a parameter; users wanting strict conservation can set it to 1.0, while a lower value (e.g. 0.5) shifts the question from "conserved" to "majority-present."

### Annotation overlap views

For each annotation stream (mutations, modifications, ligands), a long-form table of `(interface, interface residue, annotation metadata)` tuples — every place where an annotation lands on an interface residue. Used directly for cluster interpretation; also pivots cleanly into structure-by-residue heatmaps.

### Cluster interpretation report

One row per **interface interaction state** (cluster) with: cluster size, member structures, methodology summary (experimental-method distribution, resolution range, interface area range); **per-interface statistics** (residue-pair count, typed-interaction count, interface area — min / median / mean / max — and median contact density); **core contacts** (typed contacts present in ≥ `conservation_threshold` of cluster members); **enriched contacts** (UniProt-keyed interaction-pair tuples present in cluster members at significantly higher rate than rest of dataset); **enriched** mutations / modifications / ligands at the interface; per-stream annotation density; **QC warnings** (singleton, sparse fingerprint, tiny interface, mixed UniProt residue range vs the dominant cluster — the last catching polyprotein / processed-product / domain mixups within a single UniProt accession); and a narrative note.

Contact enrichment uses the same statistical heuristic as annotation enrichment: a tuple is reported when it appears on ≥ `min_count` cluster-member interfaces AND is at least `min_enrichment`× more frequent than in the rest of the dataset, OR when it is exclusive to the cluster. Tuples are formatted with gene-name labels via `partner_map` so they read directly as `Partner_A:pos-Partner_B:pos bond_type`. This surfaces interface-level cluster signatures even when no mutation/modification/ligand annotation is recorded — a cluster of structures sharing a distinctive contact pattern but no annotated mutation is now visible. For ACE2–Spike, the Omicron mutations show up both as enriched mutations (S477N, Q498R, Q493R) AND as the corresponding enriched contacts (ACE2:19-S:477, ACE2:38-S:501) — these match exactly, validating the mutation-to-contact link directly.

The methodology summary is included alongside annotation enrichment so a cluster that lacks an obvious biological correlate but is methodologically homogeneous (e.g. "all members are cryo-EM at 3.0–3.5 Å resolution") is visibly flagged as such. This addresses a known confound: clusters can reflect modelling pipeline rather than biology, and the prior cluster report did not surface the data needed to make that judgement.

The report uses an **enrichment-based** heuristic rather than a strict dominance threshold. For each cluster, an annotation label is reported when:
- it appears on ≥ `min_count` distinct cluster-member interfaces (default 2) AND its proportion-of-interfaces in the cluster is at least `min_enrichment` × the proportion in the rest of the dataset (default 2.0); OR
- it is **exclusive** to the cluster (zero occurrences elsewhere in the dataset) and appears on ≥ 1 member interface.

Counts are by *distinct interface*, de-duplicated on `(pdb_id, assembly_id, interface_id)`, so a ligand making multiple contact tuples on a single interface counts once.

This replaces an earlier "must dominate ≥ 60% of cluster members" heuristic. The earlier rule missed combinations of low-frequency markers — e.g. variant signatures where each individual mutation appears on 1–3 cluster members but together define the cluster. The enrichment rule surfaces both single dominant features and signature-style combinations.

Per-stream **annotation density** is also reported: the count of cluster-member interfaces that carry at least one at-interface mutation / modification / ligand contact. This gives a coarse "how mutated is this cluster" read independent of whether any individual annotation crosses the enrichment threshold.

A cluster is flagged as having an annotation correlate when at least one stream surfaces an enriched label. Where no enrichment is apparent the report says so explicitly (per §10) and does not speculate.

### Interface frequency JSON export

A self-contained JSON file written by `export_interface_frequency_json` (Phase 7) for downstream visualisation tooling. Filename: `{complex_id}_interface_frequencies.json`. Output directory: `Config.output_dir` (created if missing) or current working directory if `output_dir` is `None`. Top-level keys:

- `metadata` — `complex_id`, `complex_name`, `oligomeric_state`, `generated_on` (ISO date), `n_interfaces`, `representation`, `contact_definition`, `typed_contacts` (bool).
- `partners` — one object per `(unp_accession, role)` pair, with `role`, `unp_accession`, `gene_name`, `name`, `label`. `gene_name` and `name` are populated from `complex_details["participants"]` when supplied; otherwise `null`. `label` falls back through `gene_name` → `name` → `partner_map` value (with `(role 1)` / `(role 2)` suffixes stripped) → `unp_accession`. No additional API calls are made.
- `residue_frequencies` — one object per UniProt residue with `unp_accession`, `role`, `unp_residue_number`, `unp_residue_label` (canonical UniProt one-letter code + position), `n_interfaces`, `frequency` (= `n_interfaces / metadata.n_interfaces`), and `conservation_level` ∈ `{strong (≥0.80), medium (0.50–0.79), weak (0.20–0.49), rare (<0.20)}`. Each residue is counted at most once per interface. Sorted by `(role, unp_accession, unp_residue_number)`.
- `contact_frequencies` — one object per UniProt-keyed residue–residue contact with `partner_1` and `partner_2` (each `{unp_accession, role, unp_residue_number, unp_residue_label}`), `bond_type` (when `typed_contacts=True`), `n_interfaces`, `frequency`. `partner_1` is always role 1, `partner_2` always role 2; partners are not reordered. Each contact is counted at most once per interface; when `typed_contacts=False` the bond type is collapsed and `bond_type` is omitted from the row. Sorted by `(partner_1.role, partner_1.unp_residue_number, partner_2.role, partner_2.unp_residue_number, bond_type)`.

`conservation_level` is **only** present on `residue_frequencies`; `contact_types` are **only** considered as part of the `bond_type` field on `contact_frequencies`. Sort orders are explicit, so output is independent of Python set/dict iteration order. The export is independent of clustering — it is always built from the full set of retained interfaces.

---

## 10. Success criteria

The workflow is successful when it produces clusters that the user can interpret confidently — either by relating them to tracked annotations (mutations, modifications, ligands) or by acknowledging when no annotation correlate exists.

Specifically:

- **A cluster has an annotation correlate** when its member interfaces share a distinctive feature that distinguishes them from other clusters — e.g. a cluster of structures all bearing the same mutation, all with the same ligand class bound, or all carrying a specific modification at the interface. The cluster interpretation report names the correlate and identifies the interfaces it applies to.

- **A cluster may have no annotation correlate.** This is an expected outcome, not a failure. It can mean: (i) the structures genuinely differ in their interfaces for reasons the annotation pipeline does not capture (conformational state, crystal form, construct boundaries, resolution, refinement protocol), or (ii) the cluster is methodological — driven by structures from the same study using the same construct, or by an artefact of the similarity metric. The cluster interpretation report says so explicitly and does not speculate beyond what the data shows.

- **At v1 sample sizes** (typically ~10–30 interfaces per complex), claims about cluster–annotation correspondence are descriptive, not statistical. The workflow does not compute formal enrichment tests; it surfaces the composition of each cluster and lets the user judge.

- **The workflow is unsuccessful** if it fails to ingest the data correctly, produces clusters that contradict the input (e.g. identical interfaces ending up in different clusters), or generates an interpretation report that overstates what the annotations support.

---

## 11. Assumptions and limitations

This section lists the load-bearing assumptions the workflow makes about its data sources and the cases it does not handle. These are the assumptions that, if violated, would cause silently wrong results rather than visible failures.

- **PISA partner ordering is consistent within a single interface response.** The workflow does not deduplicate within an interface; if PISA were to emit the same atom-level contact in both column orders, both would appear in the set.
- **PISA partner ordering is consistent across the entries of a complex.** A sanity check verifies this and reverses role assignment for inconsistent entries; the assumption is that the check itself is reliable.
- **SIFTS provides a one-to-one mapping per `(pdb_id, auth_asym_id, auth_seq_id, ins_code)`.** The workflow does not handle two distinct author residues mapping to the same UniProt position (microheterogeneity, alternate residue identities); such cases would silently collapse in the UniProt-keyed set.
- **The `interface_interactions` endpoint provides per-residue UniProt mappings inline.** The workflow has no fallback to a separate SIFTS call; if mappings stop being inline, the cross-structure analysis breaks.
- **Insertion-code values across all endpoints are either empty, whitespace-only, or alphabetic.** Whitespace-stripping is the only normalisation; any non-standard convention (e.g. multi-character codes, special characters) is not handled.

Annotation defaults (engineered mutations only, ligand blocklist contents) are stated in §8 and exposed as parameters; they are policy choices, not assumptions about the data.

---

## 12. Triage

### v1 — current scope

- A single dimer complex per notebook run.
- Five-phase workflow per §7: data retrieval, representation building, similarity and clustering, annotation overlap, summarisation.
- Outputs per §9: structure table, similarity matrix and dendrogram, conserved residues and interaction pairs, annotation overlap views, cluster interpretation report.
- Dimer-only restriction follows from the `interface_interactions` endpoint, which currently returns data for dimers only.

### v2 — natural extensions on the same axis

- **Multiple complexes per run** — comparing the same workflow across complexes that share a partner, or aggregating findings across a complex family.
- **Higher-order oligomers** — contingent on the upstream API supporting non-dimer interfaces, since the current v1 dimer restriction is API-driven rather than design-driven.
- **Typed ligand contacts feeding into similarity** — extending the interaction-pair set to include ligand contacts as elements alongside protein–protein contacts. Changes the conceptual unit from "protein–protein interface" to "extended interaction surface" (see §8). Useful for drug discovery applications where the ligand pocket and the interface overlap.
- **Statistical enrichment tests** — Fisher's exact, hypergeometric, or similar, for cluster–annotation correspondence. Becomes meaningful at sample sizes above ~30 interfaces per complex, which v1 does not assume.

### v3 — different input source

- **Predicted complexes** (e.g. AlphaFold-Multimer, Boltz, Chai) compared against experimental structures using the same workflow. Requires a local interaction-pair extractor (since PDBe APIs do not run on user-supplied models) producing the same internal representation as defined in §6. The interaction-pair abstraction is structure-source agnostic, so the downstream similarity, clustering, and annotation overlap machinery is unchanged.

### Out of scope

- **Machine learning approaches** for similarity, clustering, or annotation prediction. The Jaccard + hierarchical clustering approach is appropriate for the sample sizes involved (~10–30 interfaces per complex), is interpretable, and avoids the data-volume and validation problems that come with model-based methods at this scale.
- **Sequence-based conservation analysis.** The workflow uses structural instances of the same complex, not sequence homologues; cross-species or cross-paralogue conservation is a separate problem (per §3).
- **Structure prediction or docking.** This is an analysis tool, not a generation tool.
- **Per-atom contact analysis.** The aggregation rule in §6 collapses atom-level detail; specific donor–acceptor geometry is intentionally outside the workflow's scope.
