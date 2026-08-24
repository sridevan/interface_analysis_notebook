# Group demo: Aggregated Interface Interaction Analysis in Dimers

**Working example:** `11gl`, which resolves to `PDB-CPX-172174` (STING homodimer, Complex Portal `CPX-2128`, UniProt Q3TBT3).

**Why this complex:** small enough to read end to end, 14 interfaces from 12 PDB entries, while still producing everything the notebook is built to show: four interface interaction states, a rewiring comparison between two of them, QC warnings on the singletons, and ligands at the interface in three states. All 14 interfaces are X-ray at 1.29 to 2.75 Å, so a difference between states cannot be dismissed as a difference in method.

All figures below were produced by the current notebook on 2026-08-24 and are reproducible by running it. They will drift as new structures are deposited.

Two contrasting examples to mention: `1spq` (triosephosphate isomerase) has an invariant interface, one cluster at every cut, and `6m0j` (Spike RBD with ACE2) has 130 interfaces from 114 entries.

---

## Running order

### 1. Identifier resolution

Type a PDB entry id rather than a complex id. `11gl` resolves to `PDB-CPX-172174`, and the analysis then covers all 12 entries of that complex rather than the one typed. Entering a non-dimer, for example `4hhb`, gives an error naming both the entry and the complex it resolved to, with the chain count.

Worth mentioning: about 4 percent of PDB entries map to more than one complex, because assemblies of one entry can differ in composition. The resolver takes the complex holding the entry's preferred assembly.

### 2. Retrieval and filtering

For this complex, 14 of 14 assemblies survive, so the exclusion list is empty. Show it anyway and explain what it would report: instances carrying an antibody Fab, peptide, nucleic acid fragment or accessory subunit are excluded, because chain correspondence cannot be determined beyond two components. On `6m0j` this line reads `21  bound_macromolecules: antibody`.

Partner-order corrections are also zero here. On `6m0j`, 36 interfaces need reversing before they are comparable.

### 3. Contact vocabulary

Phase 2 prints `{'hydrogen_bond': 169, 'salt_bridge': 35}`. Make the point that the comparison sees only what the endpoint returns at residue-pair level, and that this should be checked per complex rather than assumed.

Interfaces here carry a median of 17 residue pairs (range 4 to 20) over a median 1,771 Å² interface, which is a reasonably rich fingerprint. Contrast with a complex where the median is 4 or 5, in which case the clustering rests on very few contacts and the QC warning says so.

### 4. Clustering, and the cut

The cut sweep is the important table:

| cut | states | sizes |
|---|---|---|
| 0.30 | 9 | 5, 2, 1, 1, 1 |
| 0.40 | 6 | 8, 2, 1, 1, 1 |
| 0.50 | 5 | 9, 2, 1, 1, 1 |
| 0.55 to 0.70 | 4 | 9, 3, 1, 1 |

The default cut of 0.6 sits on a plateau running from 0.55 to 0.70, which is the argument for accepting it here. Move the cut to 0.3 during the demo to show the state count change, and make the point that the number of states is a parameter, not a finding.

### 5. The four states

| state | n | entries | resolution (median) | pairs (median) | area (median) | ligands at interface |
|---|---|---|---|---|---|---|
| 2 | 9 | 11gl, 11gm, 11gn, 4loj, 4lok, 4yp1, 6xnn | 1.84 Å | 17 | 1,860 Å² | 1YD, 2BA, V67 |
| 1 | 3 | 4kby, 4kc0, 9ltf | 2.20 Å | 11 | 1,246 Å² | A1ELY, C2E |
| 3 | 1 | 4jc5 | 2.75 Å | 4 | 1,548 Å² | none |
| 4 | 1 | 4lol | 2.43 Å | 7 | 1,487 Å² | 1YE |

States are listed largest first, since the cluster ids come from `fcluster` and their numeric order means nothing.

The two singletons carry QC warnings. State 3 gets both a singleton warning and a sparse-fingerprint warning, because its 4 residue pairs are too few to establish a state. That is the notebook telling the reader not to interpret it, and it is worth showing rather than skipping.

### 6. Rewiring between states 2 and 1

Nine interfaces against three, all X-ray, medians 1.84 and 2.20 Å:

```
Sting1:A232-Sting1:D209 hydrogen_bond    state 2: 100%   state 1:   0%   higher in A
Sting1:A232-Sting1:Y260 hydrogen_bond    state 2: 100%   state 1:   0%   higher in A
Sting1:D209-Sting1:G233 hydrogen_bond    state 2: 100%   state 1:   0%   higher in A
Sting1:D273-Sting1:H156 salt_bridge      state 2:   0%   state 1: 100%   higher in B
```

State 1 also has a smaller interface, 1,246 versus 1,860 Å², and fewer pairs, 11 versus 17.

This is the point in the demo to be careful. The table says these two groups of structures differ in these contacts. It does not say why. Candidate explanations to check before claiming anything: the ligands differ between the states (C2E, cyclic di-GMP, and A1ELY in state 1; 1YD, 2BA and V67 in state 2), the constructs may differ, and although both are X-ray the resolutions are not identical. The notebook gives the counts and the method and resolution columns to support that check; it does not perform it.

### 7. Counts, not enrichment

The columns report, for every contact and annotation, the count inside the state and the count in the rest of the dataset, at interface and at entry level:

```
1YE (1/1 interfaces = 100%, 1/1 entries; rest 0/13 = 0%)
```

Read as: present in the single interface of state 4, from one PDB entry, and in none of the other 13 interfaces. The entry count is what stops "100 percent of the state" being mistaken for independent evidence.

No fold-enrichment, p-value or significance threshold is reported. For contacts this is forced by the method: interfaces are grouped by contact similarity, so a statistic scoring those contacts against the grouping restates the grouping. Permuting cluster labels at random on a larger dataset reproduced roughly half of the contacts an earlier enrichment rule reported. A Fisher exact test with Benjamini-Hochberg control was prototyped, and it does remove the small-count noise, but it addresses neither the circularity nor the fact that repeated assemblies and laboratory-correlated depositions break the independence assumption. For annotations the comparison is not circular, but the independence problem remains, so those are reported as counts too.

If asked what was lost: nothing that was trustworthy. The previous rule reported a contact seen in one interface and absent from the other 129 as infinite-fold enrichment.

---

## Suggested close

Runtime is about 2.5 seconds for the whole notebook on this complex, excluding the Mol* rendering. On `6m0j`, 114 entries, retrieval is concurrent and the ligand step takes 4.5 seconds against 32 seconds sequentially.

The methodological change since the last demo is the move from enrichment scores to counts, and the reasoning behind it is in the Phase 5 guidance in the notebook.
