# Group demo: Aggregated Interface Interaction Analysis

**Working example:** `6m0j`, which resolves to `PDB-CPX-140195` (SARS-CoV-2 Spike RBD with ACE2).

**Why this complex:** it is the largest dataset the workflow has been run on, 130 interfaces from 114 PDB entries, so it exercises the assembly filters, the concurrency, the clustering and the QC warnings at a realistic scale. KRAS with RAF1 (`PDB-CPX-130306`) remains a smaller and cleaner teaching example.

All figures below were produced by the current notebook on 2026-08-24 and are reproducible by running it. They will drift as new structures are deposited.

---

## What to show

### 1. Identifier resolution

Type a PDB entry id, not a complex id. `6m0j` resolves to `PDB-CPX-140195` and the analysis then covers all 114 entries of that complex, not just the one typed. Entering a non-dimer, for example `4hhb`, produces a clear error naming both the entry and the complex it resolved to, with the chain count.

### 2. The assembly filter, with reasons

Of 151 assemblies listed for the complex, 21 are excluded, and the notebook states why rather than reporting a bare total:

```
Assembly filters excluded 21 of 151 interfaces:
    21  bound_macromolecules: antibody
```

Every one is a ternary complex with an antibody Fab bound alongside the dimer. Retaining them would mix Fab-perturbed interfaces into the comparison.

### 3. Defensive checks that earn their keep

- **36 interfaces** had reversed `(unp_accession_1, unp_accession_2)` ordering, corrected automatically by `check_partner_consistency`. Without this the similarity matrix would compare partner 1 of one structure against partner 2 of another.
- **0 residues** were dropped for a missing UniProt mapping in this complex, though the count is reported per interface because other complexes lose termini and expression tags here.

### 4. Retrieval speed

The ligand step issues one call per ligand instance across 114 entries. Sequentially this takes 32.0 s; on 8 workers it takes 4.5 s and returns identical contacts. All phases except the Mol* rendering complete in under 10 s.

### 5. Clustering and what it does not settle

Thirteen clusters at `cluster_distance_cut = 0.6`, of sizes 66, 44, 4, 3, 3, 2, 2 and six singletons. Two points to make explicitly during the demo:

- The number of states is a function of the cut. Show the dendrogram and the cut sweep, and change the cut to demonstrate that the count moves.
- The two large clusters have the same median resolution (3.11 Å) and both mix electron microscopy with X-ray diffraction, so neither is obviously a methodological artefact. That is a check the report supports, not a conclusion it delivers.

### 6. QC warnings doing their job

Six of the thirteen clusters carry warnings. Two singletons are flagged for a resolution gap, at 4.30 Å and 4.12 Å against 3.10 Å for the dominant cluster, with the note that missing contacts there may be undetected rather than absent. One singleton is additionally flagged for a sparse fingerprint (3 residue pairs) and a poorly overlapping residue range.

### 7. Counts, not enrichment

This is the part of the demo that has changed most, and it is worth being direct about why.

The report gives, for every contact and every annotation, the count inside the cluster and the count in the rest of the dataset, at interface and at entry level:

```
ACE2:D38-S:Q498 salt_bridge (2/2 interfaces = 100%, 1/1 entries; rest 1/128 = 1%)
```

The entry count is the important column. "2/2 interfaces" from "1/1 entries" is one deposition counted twice, which the interface count alone would hide.

No fold-enrichment, p-value or significance threshold is reported. For contacts this is forced by the method: clusters are defined by Jaccard similarity over the same contacts, so testing those contacts against the clustering is circular. Permuting cluster labels at random, preserving sizes, reproduced roughly half of the contacts that the previous enrichment rule reported. Adding Fisher's exact test with Benjamini-Hochberg control removes the small-count noise but not the circularity, so it was evaluated and not adopted. For annotations the comparison is not circular, but the structures are not independent observations, so counts are reported there too.

---

## A claim this demo previously made, and why it has been withdrawn

Earlier versions of these notes claimed that the workflow rediscovered the Omicron signature from contacts alone, citing enrichment figures such as `S477N (3, 8.5x)` alongside `ACE2:19-S:477 (30, 3.9x)`, and presented the correspondence between mutated residues and enriched contact residues as validation of the method.

That claim does not survive examination and should not be repeated.

**Mutation coverage is sparse.** In the current dataset, engineered mutation annotations at the interface reach very few structures: N501Y appears on 3 interfaces from 3 entries, S477N on 4 interfaces from 3 entries, Q498R on 4 interfaces from 3 entries. The enrichment multipliers quoted previously were ratios computed on counts of this size.

**The contacts involved are common across the dataset**, so their presence in mutation carriers is not distinctive:

| mutation | contact | carriers | non-carriers |
|---|---|---|---|
| N501Y | ACE2:D38-S:N501 hydrogen_bond | 2/3 | 14/127 |
| S477N | ACE2:S19-S:S477 hydrogen_bond | 4/4 | 48/126 |
| Q498R | ACE2:D38-S:Q498 salt_bridge | 4/4 | 59/126 |

Only the first shows a marked difference between carriers and non-carriers, and it rests on three interfaces. The other two contacts occur in roughly 38 and 47 percent of non-carriers, so observing them in every carrier says little.

What can honestly be said is narrower: the workflow surfaces which contacts occur in which structures, with the counts needed to judge the weight of the evidence, and in this dataset the annotation coverage is too thin to support a claim about variant biology. That is a more useful thing to show a group than a multiplier that turns three structures into "8.5x".

---

## Suggested running order

1. Configuration cell, changing `identifier` to show entry-id resolution.
2. Phase 1, reading out the exclusion breakdown and the partner-order corrections.
3. Phase 3, the heatmap and dendrogram, moving the cut to show the state count change.
4. Phase 5, one large cluster and one flagged singleton, drawing attention to the entry counts and the QC line.
5. Phase 5a, filtering the explorer to a single state.
6. Phase 6, one Mol* representative.
7. Close on the counts-not-enrichment point above, which is the main methodological change since the last demo.
