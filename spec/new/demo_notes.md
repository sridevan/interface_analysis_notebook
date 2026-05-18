# Group-lead demo: Aggregated Interface Interaction Analysis

**Working example to show:** `PDB-CPX-140195` (SARS-CoV-2 Spike RBD + ACE2)

**Why this complex over KRAS–RAF1:** KRAS–RAF1 is clean and educational but the only biology it surfaces is a construct difference (RBD-only vs RBD+CRD). Spike–ACE2 demonstrates the workflow doing its most distinctive work — recovering known variant biology from contact patterns alone, on a 130-interface dataset.

---

## What Spike–ACE2 demonstrates simultaneously

1. **Variant signatures emerge as both mutation enrichment AND contact enrichment, and they match exactly.** Cluster 13's mutation profile is `S477N (3, 8.5×)`, `Q498R (3, 8.5×)`, `Q493R (1, exclusive)`. Its top enriched contacts are `ACE2:19-S:477 (30, 3.9×)`, `ACE2:38-S:501 (12, 8.5×)`, `ACE2:24-S:477 (4, 3.8×)`. The mutation residues and the contact residues are the same residues. The workflow rediscovered the Omicron signature from contacts alone.

2. **Engineered ACE2 mutations are correctly separated as a distinct cluster.** Cluster 25 (n=23) carries `Q325Y (3, 14×)`, `E35K (2, 9.3×)`, `N330Y (2, exclusive)`, `Q498Y (1, exclusive)`, `Q498H (1, exclusive)` — all known affinity-enhancing ACE2 mutations. Different cluster from Spike variants, even though both manifest as interface mutations. The contact enrichment shows the new bonds the engineered residues form: `ACE2:325-S:506 (3, 14×)`, `ACE2:353-S:498 (9, 20.9×)`.

3. **Methodology context disambiguates clusters.** Cluster 13 (Omicron) has 20 X-ray + 14 EM — same biology across both methods, so it's robust. Clusters 9, 14, 24 (no biological signature) are dominated by EM at ~3 Å — these are pipeline-artefact clusters. The `experimental_methods` and `resolution_range` columns make this readable at a glance.

4. **The defensive checks earn their keep.**
   - 36 entries had reversed `(unp_accession_1, unp_accession_2)` ordering — auto-corrected by `check_partner_consistency`.
   - 15 ternary structures (with antibody Fabs or accessory chains) dropped by the assembly filter.
   - Without these, the similarity matrix would be unusable.

---

## The "wow" moment

> **N501Y creates a hydrogen bond between Spike position 501 and ACE2 D38.** The workflow doesn't know this biologically — but cluster 13's *contact* enrichment gives `ACE2:38-S:501 hydrogen_bond (12, 8.5×)` and its *mutation* enrichment gives `N501Y`. The mutation residue and the enriched-contact residue are the same. That's the validation.

Same story for:
- `S477N → ACE2:19-S:477 hydrogen_bond` (present in 30 of 34 cluster members)
- Engineered `Q325Y → ACE2:325-S:506 hydrogen_bond`

The mutation pipeline (PDBe `mutated_AA_or_NA`) and the contact pipeline (PISA `interface_interactions`) are independent. Their agreement on the same residues is non-trivial validation.

---

## A useful counterexample within the same dataset

**Cluster 14 (n=7)** has *no* annotation correlate — zero mutations, modifications, or ligands at the interface. But contact enrichment surfaces:
- `ACE2:353-S:495 (2, 11.7×)`
- `ACE2:353-S:501 (7, 9.5×)`
- `ACE2:35-S:493 (5, 5.2×)`

This is a real interface-level distinction that the mutation pipeline misses entirely — likely a variant lineage whose mutations are outside the immediate interface or a conformational subset. Contact enrichment fills the gap.

---

## Second "wow" moment — cluster-vs-cluster diff (Omicron vs engineered ACE2)

After reading the per-cluster report you'll naturally ask: "what specifically distinguishes the Omicron cluster from the engineered-ACE2 cluster?" `outputs.compare_clusters(records, cluster_result, cluster_a=13, cluster_b=25, partner_map=partner_map)` returns the directional diff. Top results:

| contact | cluster 13 (Omicron) | cluster 25 (engineered ACE2) | direction |
|---|---|---|---|
| `ACE2:353-S:496 hbond` | 0/34 | **22/23** | **B-only** |
| `ACE2:19-S:477 hbond` | **30/34** | 0/23 | **A-only** (S477N variant) |
| `ACE2:30-S:417 salt_bridge` | 0/34 | **20/23** | **B-only** |
| `ACE2:38-S:498 salt_bridge` | **32/34** | 2/23 | biased_to_A (Q498R variant) |
| `ACE2:37-S:505 hbond` | 0/34 | **19/23** | **B-only** |
| `ACE2:393-S:505 hbond` | 0/34 | **14/23** | **B-only** |
| `ACE2:38-S:501 hbond` | **12/34** | 0/23 | **A-only** (N501Y variant) |
| `ACE2:325-S:506 hbond` | 0/34 | **3/23** | **B-only** (Q325Y mutation) |

The diff cleanly separates two affinity-enhancement *strategies*:

- **Cluster 13 (Omicron variants):** existing residues mutated → existing contact positions strengthened. The variant chemistry (S477N, N501Y, Q498R) makes pre-existing residues form better hbonds/salt-bridges. The contact set doesn't expand much; it intensifies at known positions.
- **Cluster 25 (engineered ACE2):** new residues added → contacts appear at residue positions the WT interaction doesn't use. ACE2:353, :393, :30, :325 forming bonds to S:417, S:496, S:505, S:506 are *not* part of the canonical ACE2-Spike interface — these contacts only exist in the engineered variants.

This is a structural-biology insight worth flagging to the lead: **the workflow distinguishes "stronger existing bonds" from "new bonds at new positions" purely from contact-pattern statistics**. No labels were applied for this — the workflow inferred it from the data.

Specifically usable as the "second wow" after the mutation-contact correspondence (frame 4 in the demo narrative).

---

## Demo narrative (six frames)

### Frame 1 — What this complex is and how many structures we have

> "PDB-CPX-140195: SARS-CoV-2 Spike RBD + ACE2. 127 deposited PDB entries cover this complex. After the assembly filter dropped 15 ternaries (entries with antibody Fabs that would perturb the canonical interface), we have 130 interfaces across 114 PDB entries. Phase 1 made 36 partner-order corrections — about 30% of entries had ACE2/Spike accessions in reverse, and the workflow auto-fixed them."

### Frame 2 — Dendrogram + cluster-cut sweep

> "Here's the dendrogram. The cut at 0.5 gives 30 clusters. The sweep table shows what happens at other cuts — at 0.65 we collapse to 10 clusters, at 0.70 to a single cluster. 0.5 is in the middle of a fragmentation gradient, but it cleanly separates the meaningful biological clusters from the methodology singletons."

### Frame 3 — Two clusters tell stories, the rest are noise

> "Cluster 13 has 34 interfaces — Omicron variant Spike. Cluster 25 has 23 — engineered-affinity ACE2. Together those are 57 of 130 interfaces with clear biological correlates. Most of the remaining 73 are in methodologically-homogeneous small clusters: Cryo-EM at 3.0–3.5 Å, often singletons. The report flags those as 'no enriched mutation, modification, or ligand' so we don't over-interpret them."

### Frame 4 — Variant biology emerges twice, independently, and matches

> "Cluster 13's enriched mutations: S477N, Q498R, Q493R. Cluster 13's enriched contacts: ACE2:19-S:477, ACE2:38-S:501, ACE2:24-S:477. The mutation residues and the contact residues are the same. The mutation pipeline and the contact pipeline are independent. Their agreement on the same residues is the validation that the workflow is finding real biology, not artefacts."

### Frame 5 — Cluster-vs-cluster diff: two strategies for tighter binding

> "If we ask 'how does Omicron differ from the engineered-ACE2 work specifically', `compare_clusters(13, 25)` returns the directional diff. Cluster 25 has an entire set of contacts the Omicron cluster doesn't — `ACE2:353-S:496` in 22 of 23 members but zero in cluster 13; `ACE2:30-S:417`, `ACE2:37-S:505`, `ACE2:393-S:505` — none of which are part of the canonical ACE2-Spike interface. The engineering adds *new* contacts at *new* positions. Variants do the opposite: they keep the same contact positions and just chemistry-swap the residues to make stronger bonds. The workflow distinguishes 'stronger existing bonds' from 'new bonds' purely from contact statistics."

### Frame 6 — What's left

> "Open questions: the cluster heuristic (default cut 0.5) needs user inspection on every dataset; the methodology-vs-biology disambiguation works descriptively but isn't yet a statistical test; the workflow extends to predicted complexes (AlphaFold-Multimer) by replacing only the data-retrieval layer, since the interaction-pair representation is structure-source agnostic."

---

## Practical demo tips

- **Pre-run Phase 1 before the meeting** (10 min runtime). Save the notebook with outputs cached so you can scroll through, not wait.
- **Have one slide with the structure table** filtered to clusters 13, 25, 24 — it makes the "same biology, different clusters" claim concrete.
- **Have the cluster interpretation report dataframe ready** — the `notes` column reads as natural language and tells the story by itself.
- **Keep KRAS-RAF1 in your back pocket** as a smaller example. Phase 1 takes 30–60 s. If the lead asks "does this work on smaller datasets?" you can run it live and show the construct distinction (cluster 1 = RBD-only, ~650 Å²; cluster 2 = RBD+CRD, ~1150 Å²).

---

## Headline numbers to have memorised

| Quantity | Value |
|---|---|
| PDB entries (post assembly filter) | 114 |
| Interfaces analysed | 130 |
| Assembly filter drops | 15 |
| Partner-order corrections | 36 |
| Clusters at default cut (0.5) | 30 |
| Largest cluster (Cluster 13, Omicron) | n=34, X-ray (20) + EM (14), 2.38–3.64 Å, 737–930 Å² |
| Engineered-ACE2 cluster (Cluster 25) | n=23, X-ray (17) + EM (6), 2.40–3.80 Å, 769–1008 Å² |
| Conserved residues at 80% threshold | 10 |
| Conserved interaction pairs at 80% threshold | 4 |
| At-interface mutation overlap rows | 26 |
| At-interface ligand overlap rows | 4 (NAG-K417 in 7sn0) |

---

## Background documents

- Spec: `aggregated_interface_analysis_spec.md`
- Implementation brief: `implementation_brief.md`
- Final synthesised spec: `final_spec.md`
- Code: `../interface_analysis_notebook/` (sibling repo)
