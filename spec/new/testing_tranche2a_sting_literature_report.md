# STING literature grounding: final report

Close-out summary of the literature-grounding analysis of the biological
observations recovered by the aggregated interface workflow on the STING working
example, `11gl` / `PDB-CPX-172174`, as frozen in the tranche-2A regression
fixture. Completed 2026-09-23.

The detailed analysis is `testing_tranche2a_sting_literature.md`. This document
is the summary report.

No production code, no regression test, no fixture and no other specification
document was changed. The offline suite still passes at 92 tests.

---

## Papers reviewed

Eight primary structural or functional papers covering ten of the twelve
entries, plus five reviews and contextual sources, plus PDBe and UniProt records
for all twelve entries and ten chemical components. Four entries (`11gl`,
`11gm`, `11gn`, `9ltf`) carry no journal citation and are marked
to-be-published.

**An access limitation shapes every negative finding.** Gao 2013 in *Cell*
returned HTTP 403 and PubMed Central presented a bot check for Zhang 2013, so
those two rest on abstracts and search-engine extraction rather than a full
reading. Interface contact tables are exactly what tends to live in
supplementary material, so "not identified in the literature searched" is weaker
here than it would be after full-text review. This caveat applies to all
negative findings in the analysis document.

---

## Main established biology

- STING's C-terminal domain is a constitutive homodimer with an extensive
  hydrophobic interface, around 916 Å² buried per monomer.
- Cyclic dinucleotides bind at the dimer interface itself, which is why a ligand
  can reorganise inter-protomer contacts at all.
- Ligand binding drives an open-to-closed transition in which a four-stranded
  antiparallel β-sheet lid forms from residues 219–249 of each protomer (human
  numbering).
- Y167, R238, Y240, E260 and T263 are the characterised
  cyclic-dinucleotide-coordinating residues.
- SR-717 is a direct cGAMP mimetic inducing the same closed conformation; DMXAA
  and CMA are mouse-specific agonists.
- The words "open" and "closed" are used in two incompatible senses across
  papers, and apo mouse STING is reported as already closed. This is the single
  biggest interpretation hazard in the dataset.

---

## Exact rewiring contacts

**PDBe-KB Complexes aggregates equivalent complex instances by mapped component
identity and stoichiometry**, represented by UniProt accessions; orthologues
carry different accessions and so form different complexes. (The rule constrains
component identity, not organism count: a single complex may span organisms, as
Spike RBD with ACE2 does.) This complex's composition is one accession at
stoichiometry two, so **all twelve entries contain mouse STING Q3TBT3 by
construction, and its interface states cannot be caused by differences between
mouse and human STING orthologues.** Human comparisons serve only to relate
the mouse observations to the largely human-focused literature and to test
cross-species transferability of the interpretation.

Needleman-Wunsch alignment of mouse Q3TBT3 against human Q86WV6 gives a uniform
offset across the cyclic-dinucleotide-binding domain: **mouse residue *n*
corresponds to human *n*+1**, confirmed independently by the known mouse R231 /
human H232 pair.

| Contact | Exact pair described? | Residues individually implicated? | State-dependent rewiring described? | Confidence |
|---|---|---|---|---|
| A232–D209 | No | Partially; A232 sits in the lid, neither residue has a documented individual role | No | Moderate that it is unpublished; high that it is real in the frozen data |
| A232–Y260 | No | Partially, with a numbering trap: human E260 is mouse E259, **not** mouse Y260 (which is human Y261) | No | As above |
| D209–G233 | No | Yes for G233; human Gly234 appears in a documented lid contact, but paired with Tyr245 | No | As above |
| D273–H156 | No | No | No | Valid within mouse STING; **not transferable**, since mouse D273 corresponds to human **Y274** and the salt bridge cannot exist in human STING |

**Two further core contacts of the large state are documented at pair level.**
Zhang 2013 describes the Tyr245 side chain with the Gly234 main-chain carbonyl,
and the Ser243 side chain with the Lys236 main-chain amide (human numbering).
These are exactly the workflow's **G233–Y244** and **S242–K235** at conserved
positions. That is genuine external agreement on two of the six core contacts of
the nine-member state.

---

## Cluster interpretation

**The nine-member state is coherent and matches published structural work.**
Every member is a cyclic dinucleotide complex (2'3'-cUMP-AMP, 3'3'-cUMP-AMP,
c-di-AMP, 2'3'-cGAMP, 3'3'-cGAMP) or SR-717, which was independently shown to
induce the cGAMP-like closed conformation. Ligand class and interface state agree
across seven entries and four chemical scaffolds.

**The three-member state is mixed and needs qualification.** It groups a
c-di-GMP complex with an **apo** structure and an unpublished small-molecule
complex. Its coherence is the absence of the lid contact network, not a shared
published label. The presence of the apo structure argues against labelling it
by ligand status.

**The `4jc5` singleton is not interpretable.** Four mapped contacts at 2.75 Å,
the dataset's lowest resolution; the workflow flags it as both a singleton and a
sparse fingerprint, which is the desired behaviour.

**`4lol` is the interesting case.** Gao 2013 calls DMXAA-bound mouse STING
closed, yet the workflow separates it from the CDN state. Inspection shows why:
it retains S242–K235 but lacks the D209 / A232 / G233 / Y260 network. That is
a **partial lid-contact network** — a difference in what is being measured
(which specific contacts are present, versus a global conformational metric),
not a contradiction. The wording is deliberately about contact composition and
does not assert a particular global structural mechanism.

---

## Potentially interesting workflow-derived findings

Not identified in the literature searched, stated cautiously:

- The state-dependent presence of A232–D209, A232–Y260 and D209–G233 at 9/9
  versus 0/3.
- The partial lid-contact network of DMXAA (`4lol`), retaining one documented
  lid contact while lacking the rest of the network.
- D273–H156 marking the three-member state at 3/3 versus 0/9, which is
  mouse-specific at best because human STING carries Y274 at that position.

None of these is described as novel. Each requires further validation.

---

## ZNT

`ZNT` is **2'3'-cUMP-AMP**, a cyclic dinucleotide (C19H23N7O14P2, MW 635.4), not
a zinc species despite the resemblance to the blocklisted `ZN` code. It is the
ligand of `11gl` and contacts mouse residues 231, 237 and 262, whose human
equivalents (H232/R232, R238, T263) are among the best-characterised
cyclic-dinucleotide pocket residues.

It is a biologically meaningful agonist complex, and its absence from the demo
notes looks like an oversight in the notes rather than a workflow problem. No
blocklist change or production annotation change is needed.

---

## Recommendation

**STING is primarily a validation case study: the aggregated interface analysis
recovers established conformational biology while expressing it as residue-level
contact rewiring. It also highlights specific contact-level observations not
identified explicitly in the literature searched, which require further
validation before any novelty claim.**

The defensible claim is that an unbiased aggregation over deposited structures,
using no conformational metric and no structural superposition, recovers a
mechanism that took the field several dedicated crystallographic studies to
establish, and expresses it at the level of named residue pairs rather than
global angles. Two of those pairs independently match the primary literature,
which validates the representation.

The four evidence classes should be kept distinct when any of this is written up:
**established biology reproduced** (the CDN/mimetic state separation, SR-717
co-clustering); **mechanistic interpretation supported** (state-defining contacts
in the lid region, ligand contacts on canonical pocket residues, c-di-GMP
separating from other cyclic dinucleotides); **workflow-derived observation**
(the state-dependent A232–D209, A232–Y260 and D209–G233 frequencies, the `4lol`
partial lid-contact network, D273–H156); and **requires further validation**,
which applies to every item in the third class.

Three cautions should travel with it:

1. State the species and give both numbering systems, since D273 does not exist
   in human STING.
2. Avoid "open" and "closed" without defining the sense, given the conflicting
   usage and the report that apo mouse STING is already closed.
3. Do not present four states as four conformational states of STING; that
   number is a function of the chosen cut height.

Two things would strengthen it materially, neither proposed for implementation
now:

- Compute a global conformational metric such as α2-tip separation for each
  instance, to test directly whether the workflow states align with the
  published open/closed axis and to quantify where `4lol` sits.
- Run the same analysis on human STING complexes, where the literature is deeper
  and D273 does not exist, to test which contacts transfer across species.

---

## Sources

- [Gao *et al.* 2013, *Cell*](https://pubmed.ncbi.nlm.nih.gov/23910378/) — PDB 4LOJ, 4LOK, 4LOL
- [Chin KH *et al.* 2013, *Acta Cryst D*](https://journals.iucr.org/paper?S0907444912047269) — PDB 4KBY, 4KC0
- [Chin EN *et al.* 2020, *Science*](https://www.science.org/doi/10.1126/science.abb4255) — PDB 6XNN, 6XNP
- [Cavlar *et al.* 2013, *EMBO J*](https://pubmed.ncbi.nlm.nih.gov/23604073/) — PDB 4JC5
- [Ouyang *et al.* 2012, *Immunity*](https://pubmed.ncbi.nlm.nih.gov/22579474/)
- [Zhang *et al.* 2013, *Molecular Cell*](https://pmc.ncbi.nlm.nih.gov/articles/PMC3808999/)
- [Konno *et al.* 2014, *Cell Reports*](https://pubmed.ncbi.nlm.nih.gov/25199835/)
- [*Frontiers in Immunology* 2022 review](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.808607/full)
- [*Nature Communications* 2025](https://www.nature.com/articles/s41467-025-58641-5)
- [*Molecular Cell* 2023](https://www.cell.com/molecular-cell/fulltext/S1097-2765(23)00243-5)

---

## Document sequence

```text
spec/new/testing_strategy.md                        long-term plan
spec/new/testing_tranche1_report.md                 first offline suite, design questions
spec/new/testing_tranche1_consolidation.md          suite streamlined to 83 items
spec/new/testing_tranche1_5_report.md               exclusion rule, partial mappings, homodimer diagnostic
spec/new/testing_tranche1_5_closeout.md             spec and consistency update
spec/new/testing_tranche2a_sting_e2e.md             frozen STING end-to-end regression
spec/new/testing_tranche2a_final_report.md          tranche 2A close-out
spec/new/testing_tranche2a_sting_literature.md      STING literature grounding (detailed)
spec/new/testing_tranche2a_sting_literature_report.md   this document
```
