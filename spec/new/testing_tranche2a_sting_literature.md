# STING literature grounding

Literature-grounding analysis of the biological observations recovered by the
aggregated interface workflow on the STING working example, `11gl` /
`PDB-CPX-172174`, as frozen in the tranche-2A regression fixture.

Written 2026-09-23. This is an interpretation document only. No production code,
no regression test and no other specification document was changed while writing
it. The regression assertions protect reproducible workflow behaviour and are a
separate concern from the biological interpretation attached to that behaviour.

Throughout, the phrase **not identified in the literature searched** is used in
preference to "novel". Failure to find a published description is a statement
about the search, not about the world.

---

## Scope and search method

### Databases and tools

| Resource | Use | Access achieved |
|---|---|---|
| PDBe REST API (`www.ebi.ac.uk/pdbe/api`) | entry summaries, primary citations, source organism, chemical components | full, programmatic |
| UniProt REST (`rest.uniprot.org`) | Q3TBT3 (mouse) and Q86WV6 (human) sequences, domain and binding-site features | full, programmatic |
| RCSB PDB entry pages | cross-check of entry metadata | via search results |
| General web search | primary literature discovery | full |
| IUCr Journals | Chin *et al.* 2013 abstract | full abstract |
| Cell / ScienceDirect | Gao *et al.* 2013 full text | **blocked (HTTP 403)** |
| PubMed Central | Zhang *et al.* 2013, review PMC9841179 | **blocked (reCAPTCHA)** |
| Europe PMC REST | attempted metadata retrieval | returned version stub only |

### Principal search terms

STING dimer interface; cyclic dinucleotide binding; open versus closed
conformation; four-stranded beta-sheet lid; lid region residues; mSTING c-di-GMP
recognition modes; DMXAA species specificity; SR-717 cGAMP mimetic; combinations
of residue identifiers in one-letter and three-letter form (`D209`/`Asp209`,
`A232`/`Ala233`, `G233`/`Gly234`, `Y260`/`Tyr261`, `D273`, `H156`, and human
equivalents); `Tyr245 Gly234`, `Ser243 Lys236`; systematic comparison of STING
crystal structures; aggregated interface contact analysis.

Searches were run on 2026-09-23.

### Important limitation on full-text access

Two of the most directly relevant primary papers, Gao *et al.* Cell 2013 and
Zhang *et al.* Molecular Cell 2013, could not be retrieved as full text: the
publisher returned 403 and PubMed Central presented a bot check. Statements
attributed to them below rest on abstracts, publisher summaries and
search-engine extracts of their text, not on a full reading. Any claim of the
form "not identified in the literature searched" is therefore weaker than it
would be after a full-text review of those two papers plus their supplementary
material, which is where interface contact tables are most likely to appear.
**This caveat applies to every negative finding in this document.**

### A naming collision worth recording

Literature searches combining "STING" with "interface" and "contacts" are
confounded by an unrelated bioinformatics package also called STING (Sequence To
and INtegrated Graphics; "STING Contacts", "Blue Star STING"), which computes
residue contacts in arbitrary proteins. Searches were filtered accordingly. A
future search on this topic should anticipate the collision.

---

## Structural context of the frozen STING dataset

All fourteen interface instances in the frozen fixture come from twelve PDB
entries, and **every entry is mouse STING** (*Mus musculus*, UniProt Q3TBT3,
gene *Sting1*, formerly *Tmem173*). All are X-ray structures between 1.29 and
2.75 Å.

This is a property of the aggregation, not a coincidence of this dataset.
**PDBe-KB Complexes aggregates equivalent complex instances by mapped component
identity and stoichiometry**, represented by UniProt accessions. Because
orthologues carry different accessions, a complex containing human STING would
have a different composition and a different complex identifier. (The rule is
about component identity, not about a complex being single-organism: a complex
may legitimately span organisms, as Spike RBD with ACE2 does.)

For this complex the composition is a single accession at stoichiometry two, so
**all instances in this STING aggregation contain mouse STING Q3TBT3, and its
interface states cannot be caused by differences between mouse and human STING
orthologues.** The states and rewiring described below are structural variation
among mouse STING structures.
Comparisons with human STING are used here only to relate the mouse observations
to the broader, largely human-focused literature and to assess cross-species
transferability of the interpretation; they play no part in producing the
states. Every such comparison requires explicit residue correspondence, which is
established below.

| PDB | Ligand (CCD) | Chemical identity | Resolution | Primary citation |
|---|---|---|---|---|
| `11gl` | ZNT | 2'3'-cUMP-AMP | 1.98 Å | to be published (released 2026-07-29) |
| `11gm` | A1AEP | 3'3'-cUMP-AMP | 1.84 Å | to be published (released 2026-07-29) |
| `11gn` | 2BA | 3'3'-c-di-AMP | 1.29 Å | to be published (released 2026-07-29) |
| `4loj` | 1SY | 2'3'-cGAMP, c[G(2',5')pA(3',5')p] | 1.77 Å | Gao *et al.* Cell 2013 |
| `4lok` | 1YD | 3'3'-cGAMP, c[G(3',5')pA(3',5')p] | 2.07 Å | Gao *et al.* Cell 2013 |
| `4lol` | 1YE | DMXAA | 2.43 Å | Gao *et al.* Cell 2013 |
| `4yp1` | 2BA | c-di-AMP | 2.65 Å | Zhu *et al.* Biochemistry 2015 |
| `6xnn` | V67 | SR-717 | 2.49 Å | Chin EN *et al.* Science 2020 |
| `4kby` | C2E | c-di-GMP | 2.36 Å | Chin KH *et al.* Acta Cryst D 2013 |
| `4kc0` | none | apo | 2.20 Å | Chin KH *et al.* Acta Cryst D 2013 |
| `9ltf` | A1ELY | telatinib analogue | 1.92 Å | to be published (released 2026-02-11) |
| `4jc5` | 1K5 | CMA (10-carboxymethyl-9-acridanone) | 2.75 Å | Cavlar *et al.* EMBO J 2013 |

Four entries (`11gl`, `11gm`, `11gn`, `9ltf`) carry no journal citation in PDBe
and are marked to-be-published. Their structural interpretation therefore rests
on the deposition titles and on this analysis, not on peer-reviewed text.

One citation is worth flagging: the primary citation attached to `4yp1` is a
Biochemistry 2015 paper whose title concerns cyclic di-AMP binding to
*SaCpaA_RCK*, a different protein. The mouse STING/c-di-AMP structure appears in
that work as a comparison, so `4yp1` is not the subject of its own dedicated
structural paper.

### Species and residue numbering correspondence

Mouse Q3TBT3 (378 aa) and human Q86WV6 (379 aa) were aligned by
Needleman-Wunsch (identity scoring; 267 identities, 70.6% of the mouse
sequence). Across the whole cyclic-dinucleotide-binding domain the offset is
uniform: **mouse residue *n* corresponds to human residue *n*+1**. This is
independently confirmed by the well-known allele pair in the literature, mouse
R231 and human H232/R232, which the alignment reproduces exactly.

| Mouse (workflow numbering) | Human equivalent | Conserved? |
|---|---|---|
| H156 | H157 | yes |
| S161 | S162 | yes |
| D209 | D210 | yes |
| N210 | N211 | yes |
| R231 | H232 (R232 allele) | **no** (species/allele difference) |
| A232 | A233 | yes |
| G233 | G234 | yes |
| K235 | K236 | yes |
| R237 | R238 | yes |
| V238 | V239 | yes |
| Y239 | Y240 | yes |
| S242 | S243 | yes |
| Y244 | Y245 | yes |
| Y260 | Y261 | yes |
| T262 | T263 | yes |
| T266 | T267 | yes |
| **D273** | **Y274** | **no** |
| K275 | Q276 | **no** |

Two entries in this table carry real interpretive weight:

- **Mouse D273 corresponds to human Y274, not to an acidic residue.** The
  D273–H156 salt bridge recovered by the workflow *cannot* form in human STING,
  because the human protein has tyrosine at the equivalent position. This is a
  statement about **transferability, not about the state**: the contact is a
  valid, reproducible state-associated feature *within* mouse STING, and since
  the aggregation contains only mouse structures it cannot have arisen from a
  sequence difference between orthologues. What does not carry over is the
  *interpretation*, which should not be generalised to human STING or to the
  STING family.
- Mouse I229 corresponds to human G230. This is the natural species difference
  in the lid that underlies DMXAA specificity (see below).

---

## Published STING conformational biology

The following points are firmly established in the primary structural
literature and are not in dispute.

1. **STING's C-terminal ligand-binding domain is a constitutive homodimer.**
   The dimer is described as V-shaped or butterfly-like, with each protomer
   forming a wing and a crevice between them. The dimer interface is extensive
   and largely hydrophobic, reported at roughly 916 Å² buried per monomer, and
   is mediated principally by the first long kinked helix (LBDα1), with
   additional contributions from helices α5 and α7 (Ouyang *et al.* Immunity
   2012).

2. **Cyclic dinucleotides bind at the dimer interface**, in the cleft between
   the two protomers, rather than within a single protomer. This places the
   ligand and the protein-protein interface in direct structural competition,
   which is the reason a CDN can reorganise inter-protomer contacts at all.

3. **Ligand binding drives an open-to-closed transition in which a
   four-stranded antiparallel β-sheet "lid" forms over the ligand.** The lid and
   its connecting loops are formed by residues **219–249 of each protomer**
   (human numbering; Zhang *et al.* Molecular Cell 2013), which corresponds to
   mouse 218–248. On closure the two protomers rotate inward, deepening the
   pocket. Reported quantifications include the α2-tip separation falling into
   two narrow ranges, roughly 47–54 Å (open) and 34–35 Å (closed), and a large
   rotation of the ligand-binding domain relative to the transmembrane domain in
   the full-length protein.

4. **Specific CDN-coordinating residues are well characterised.** Y167, R238,
   Y240, N242, E260 and T263 (human) are repeatedly identified in CDN
   recognition, with Y167 and R238 conserved across STING homologues. Alanine
   substitution of Y167, Y240 or R238 abolishes ligand binding and the cellular
   interferon response. UniProt annotates mouse binding sites at 165–166,
   237–240 and 262, consistent with this.

5. **Non-nucleotide agonists can reproduce the CDN-bound closed state.** SR-717
   was shown by a 1.8 Å co-crystal structure to be a direct cGAMP mimetic that
   induces "the same closed conformation of STING" (Chin EN *et al.* Science
   2020); `6xnn` is the mouse structure from that work. Not all non-nucleotide
   agonists behave this way: diABZI compounds have been reported to retain the
   open conformation, so agonist class does not map onto conformation uniformly.

6. **DMXAA and CMA are mouse-specific agonists.** DMXAA binds mouse STING as two
   copies in the pocket and was reported to induce a closed conformation
   (Gao *et al.* Cell 2013). Human STING can be rendered DMXAA-responsive by
   engineered substitutions: S162A in the binding pocket, G230I in the lid, and
   Q266I, acting cooperatively. Note carefully that **S162A and Q266I are
   engineered substitutions, not humanisation-to-mouse changes** — mouse carries
   serine and glutamine at the equivalent positions (S161, Q265). Only G230I
   matches the natural mouse residue (I229). CMA is likewise mouse-selective
   (Cavlar *et al.* EMBO J 2013).

### A terminology hazard that must be handled explicitly

The words "open" and "closed" are used in at least two incompatible senses in
this literature.

- **Gao 2013 and most later work**: closed = the lid-formed, CDN-bound,
  activation-competent state, contrasted with the apo state.
- **Chin KH *et al.* Acta Cryst D 2013** (the source of `4kby` and `4kc0`):
  their *mouse* structures, including the **apo** form, are described as
  "novel closed-form structures" exhibiting considerable differences from
  "previously reported open-form human STING-CTD structures". Here "closed"
  denotes a more compact mouse dimer relative to open human structures.

Independent sources state that apo human STING is open and closes on cGAMP
binding, whereas **apo mouse STING is already in a closed conformation**. That
statement is directly relevant here, because the frozen dataset is entirely
mouse and includes an apo structure.

The consequence for this document is stated plainly: **the workflow's two large
states should not be labelled "open" and "closed" without qualification.** They
are distinguished by the identity of specific inter-protomer contacts, which is
a finer and less ambiguous discriminator than a global conformational label.

---

## Workflow states versus published structural states

Cluster identifiers from `fcluster` are arbitrary and are not used. States are
named by membership.

| Workflow state | Members | Ligands reaching interface residues | Published structural interpretation | Relationship to workflow state |
|---|---|---|---|---|
| **large 9-member state** | `11gl`, `11gm` (2 assemblies), `11gn` (2 assemblies), `4loj`, `4lok`, `4yp1`, `6xnn` | ZNT (2'3'-cUMP-AMP), 1YD (3'3'-cGAMP), 2BA (c-di-AMP), V67 (SR-717) | all are cyclic-dinucleotide or confirmed cGAMP-mimetic complexes; `4loj`/`4lok` are the Gao 2013 closed-conformation structures; `6xnn` is SR-717, shown to induce the cGAMP-like closed state | **Coherent.** Every member is a CDN or CDN-mimetic complex, and the state carries the lid contact network described below. This is the closest thing in the dataset to a recognisable published grouping. |
| **3-member state containing `4kby`** | `4kby`, `4kc0`, `9ltf` | C2E (c-di-GMP), A1ELY (telatinib analogue) | `4kby` and `4kc0` are the Chin 2013 mouse c-di-GMP and **apo** structures, called "closed-form" in their own paper's sense; `9ltf` is unpublished | **Mixed and requires qualification.** It groups an apo structure with a c-di-GMP complex and an unpublished small-molecule complex. Its coherence is that all three lack the lid contact network, not that they share a published label. The presence of the apo structure argues against calling this state "ligand-bound" anything. |
| **`4lol` singleton** | `4lol` | 1YE (DMXAA) | Gao 2013 describes DMXAA-bound mouse STING as adopting a closed conformation | **Apparent tension, see below.** Literature-"closed" but separated from the CDN closed state by contact identity. |
| **`4jc5` singleton** | `4jc5` | none reaching interface residues | Cavlar 2013, mouse-specific CMA complex | **Not interpretable as a state.** Four mapped residue pairs at 2.75 Å, the lowest resolution in the dataset; the workflow flags it as both a singleton and a sparse fingerprint. |

### The `4lol` tension, stated carefully

Gao 2013 reports DMXAA-bound mouse STING as closed, in the same sense as the
CDN complexes from the same paper. The workflow nonetheless separates `4lol`
from the nine-member CDN state. Inspection of the frozen contacts shows this is
not a contradiction but a difference in what is being measured:

- `4lol` **does** carry S242–K235, one of the two lid contacts documented in the
  literature (see next section), and also Y244–K235.
- `4lol` **lacks** D209–A232, A232–Y260, D209–G233 and G233–Y244, the rest of
  the network that all nine CDN structures share.

A global conformational descriptor such as α2-tip separation can score `4lol` as
closed while its specific inter-protomer contact set differs from the CDN
complexes. DMXAA is a 282 Da xanthenone present as two copies, against
635–690 Da cyclic dinucleotides present as one; a different pocket occupancy
producing a partially different contact set is a plausible reading, and the
frozen data is consistent with it. This should be described as a
**partial lid-contact network**, a statement about which contacts are present
and absent rather than about a global structural mechanism. It was **not
identified in the literature searched** as an explicit contact-level
statement.

---

## State-separating residues

Evidence is reported at two levels, kept strictly separate: whether the
individual residue has an established structural role, and whether the specific
residue-residue contact has been described.

| Mouse residue | Human equiv. | In lid (human 219–249)? | Established individual role in the literature searched |
|---|---|---|---|
| D209 | D210 | no (immediately N-terminal to the lid region) | No specific functional or structural role for Asp210 was identified in the literature searched. It lies in the α-helical region preceding the lid. |
| A232 | A233 | yes | No residue-specific role identified. Adjacent to the polymorphic/allelic position 232 (human) / 231 (mouse), which is well characterised: human R232 contacts the CDN α-phosphates and residue 232 is explicitly described as part of the β-sheet lid over the binding pocket. |
| G233 | G234 | yes | **Yes.** Gly234's main-chain carbonyl is named in the literature as a partner in an interdomain lid contact (see next section). |
| Y260 | Y261 | no (C-terminal to the lid) | Adjacent to human E260, which is repeatedly listed among key CDN-coordinating residues. Note the off-by-one hazard: **human E260 is mouse E259, not mouse Y260.** Mouse Y260 is human Y261. No specific published role for Tyr261 was identified in the literature searched. |
| D273 | **Y274** | no | **Not conserved.** No role identified; the human protein cannot present an acidic side chain here. |
| H156 | H157 | no | No specific role identified in the literature searched. |

The off-by-one hazard flagged for Y260 is worth emphasising, because a careless
search would "confirm" mouse Y260 using literature about human E260. They are
different residues one position apart in different species. The workflow's Y260
is human Y261.

Two further residues that appear in the large-state core contacts do have
established roles: mouse R231 (human H232/R232, the allelic CDN-contacting lid
residue) and mouse Y244 and K235 and S242 (human Y245, K236, S243), which appear
in the documented lid interaction set.

---

## Exact contact-pair evidence

The four contacts protected by the regression, in the production representation
(role 1 residue, role 2 residue, bond type), all on Q3TBT3, together with the
two additional core contacts that turned out to be literature-documented:

| Contact (mouse) | Human equivalent | 9-member state | 3-member state | `4jc5` | `4lol` |
|---|---|---|---|---|---|
| A232–D209 hydrogen bond | A233–D210 | 9/9 | 0/3 | 0/1 | 0/1 |
| A232–Y260 hydrogen bond | A233–Y261 | 9/9 | 0/3 | 0/1 | 0/1 |
| D209–G233 hydrogen bond | D210–G234 | 9/9 | 0/3 | 0/1 | 0/1 |
| D273–H156 salt bridge | *(Y274)*–H157 | 0/9 | 3/3 | 0/1 | 0/1 |
| G233–Y244 hydrogen bond | G234–Y245 | 9/9 | 0/3 | 0/1 | 0/1 |
| S242–K235 hydrogen bond | S243–K236 | 9/9 | 0/3 | 0/1 | **1/1** |

### What the literature does and does not say about these pairs

**Documented at pair level.** Zhang *et al.* Molecular Cell 2013 describe
interdomain interactions within the lid of the STING dimer including
"the side group of Tyr245 and the main-chain carbonyl oxygen atom of Gly234"
and "the side group of Ser243 and the main-chain amide nitrogen atom of
Lys236" (human numbering). These correspond exactly to the workflow's mouse
**G233–Y244** and **S242–K235** contacts. Both residue pairs are conserved
between mouse and human. This is genuine pair-level agreement: two of the
workflow's large-state core contacts reproduce inter-protomer lid interactions
that were described in the primary literature from an independent structure.

Note the limitation: this passage was recovered through search-engine extraction
of the paper's text, not a full reading, and the paper's structure is human
STING with 2'3'-cGAMP. The agreement is with a human structure, at conserved
positions, not with a mouse structure.

**Not identified at pair level.** For the four contacts the regression
specifically protects:

| Contact | A. both residues individually known to matter | B. residues in the relevant structural region | C. the exact residue-residue contact described | D. state-dependent gain/loss of that contact described |
|---|---|---|---|---|
| A232–D209 | partial: neither Ala233 nor Asp210 has a documented individual role; A233 is adjacent to the characterised position 232 | yes: A232 is in the lid, D209 immediately precedes it | **not identified in the literature searched** | **not identified in the literature searched** |
| A232–Y260 | partial, as above; Y260 (human Y261) has no documented role and must not be confused with human E260 | partially: A232 in the lid, Y260 outside it | **not identified in the literature searched** | **not identified in the literature searched** |
| D209–G233 | partial: Gly234 (human) is named in a documented lid contact, but with Tyr245, not with Asp210 | yes | **not identified in the literature searched** | **not identified in the literature searched** |
| D273–H156 | no: neither residue has a documented role, and D273 is not conserved in human | no: both outside the lid | **not identified in the literature searched** | **not identified in the literature searched** |

Applying the required standard to the clearest case: *Gly234 participates in a
documented inter-protomer lid interaction in the cGAMP-bound STING dimer, which
makes the observed D209–G233 rewiring structurally plausible; however, I did not
identify a publication explicitly describing the D209–G233 (human D210–G234)
contact or its state-dependent formation.*

### The one observation that needs the strongest qualification

The **D273–H156 salt bridge** is the contact that most cleanly marks the
3-member state (3/3 versus 0/9). Within mouse STING it is a genuine
state-associated contact: the aggregation is entirely mouse, so the difference
between the two states is structural variation among mouse structures and cannot
be an artefact of comparing orthologues.

It is, however, **the least transferable observation in the dataset**. Mouse
D273 corresponds to human **Y274**, so the same salt bridge is structurally
impossible in human STING. The distinction to maintain is between
*within-species structural variation*, which is what the workflow measured and
which is sound, and *cross-species transferability of the interpretation*, which
here is absent. It should be presented as a mouse STING interface feature and
never as a general property of the STING dimer. It remains a legitimate
regression property of the frozen analysis either way.

---

## Ligand-state relationships

### Do similar ligands cluster together?

Largely yes, and the exception is informative.

- Every member of the 9-member state is a cyclic dinucleotide complex
  (2'3'-cUMP-AMP, 3'3'-cUMP-AMP, c-di-AMP, 2'3'-cGAMP, 3'3'-cGAMP) or SR-717,
  a compound independently shown to be a direct cGAMP mimetic inducing the same
  closed conformation. Ligand class and interface state agree across seven
  distinct entries and four distinct chemical scaffolds.
- The 3-member state contains a c-di-GMP complex, an **apo** structure and an
  unpublished small-molecule complex. The presence of an apo structure means
  ligand identity cannot be the whole explanation for this state.
- c-di-GMP (`4kby`) sits in a different state from every other cyclic
  dinucleotide in the dataset. This is consistent with published reports that
  c-di-GMP and 2'3'-cGAMP produce distinct STING conformations, including the
  observation that the α2-helices form a larger angle in the V-shaped c-di-GMP
  complex than in the U-shaped cGAMP complex, and with dedicated work on the
  distinct dynamical and conformational features of human STING in response to
  2'3'-cGAMP versus c-di-GMP. It is also consistent with c-di-GMP being a lower
  affinity, bacterial-derived ligand relative to the endogenous 2'3'-cGAMP.

### Confounders that cannot be excluded

The frozen dataset cannot separate ligand identity from other explanations, and
the following should be stated whenever the ligand association is mentioned:

- **Crystallisation and construct.** `4kby` and `4kc0` come from the same study
  and are explicitly described as having been obtained "under different
  crystallization conditions" from earlier work, with construct mSTING(137–344).
  The other entries use varying construct boundaries (observed polypeptide
  lengths in the fixture range from 185 to 210 residues). Two structures from
  one laboratory, one crystal form and one construct sharing an interface state
  is not independent evidence of a ligand effect.
- **Deposition era.** The 3-member state mixes 2013 and 2026 depositions, and
  the 9-member state mixes 2013, 2015, 2020 and 2026, so era alone does not
  explain the split, which is reassuring.
- **Resolution.** The 9-member state has a median of 17 mapped residue pairs;
  the 3-member state 11; `4lol` 7; `4jc5` 4. `4jc5` is the lowest-resolution
  entry at 2.75 Å. Contact detection is geometry-derived and
  resolution-dependent, which the workflow's own QC warning states.

Accordingly the correct formulation is: cyclic-dinucleotide and cGAMP-mimetic
complexes are **associated with** the lid contact network, and this association
is **consistent with** the established ligand-induced lid closure mechanism.
The frozen data does not establish causation.

### Does the `4lol` singleton make structural sense?

Partially, and it is the most interesting single-structure observation in the
dataset. DMXAA is a mouse-specific non-nucleotide agonist that binds as two
copies. The workflow maps its interface contact onto mouse S161, whose human
equivalent S162 is precisely the residue whose S162A substitution renders human
STING DMXAA-sensitive. That the annotation pipeline places the DMXAA contact on
the known species-specificity determinant is a meaningful independent check on
the residue mapping. The singleton status reflects partial rather than absent
lid engagement, as set out above.

### Does the `4jc5` sparse singleton make sense?

It should not be over-interpreted, and the workflow says so itself. `4jc5` is
the CMA complex at 2.75 Å with four mapped residue pairs, and carries both the
singleton and sparse-fingerprint QC warnings. CMA is, like DMXAA, a
mouse-specific agonist, so a distinct interface state would not be surprising;
but with four contacts at the dataset's lowest resolution, the state cannot be
distinguished from a detection artefact. The correct reading is that the
workflow has flagged it as uninterpretable, which is the desired behaviour.

---

## ZNT observation

`ZNT` was absent from the earlier demo notes and was queried separately.

**`ZNT` is 2'3'-cUMP-AMP**, a cyclic dinucleotide: formula C19H23N7O14P2,
molecular weight 635.4, PDBe compound name "2'3'-cUA". It is not a zinc species;
the three-letter code is coincidentally similar to `ZN`, which is on the default
ligand blocklist, but `ZNT` is a distinct chemical component and is correctly
not blocklisted.

It is the ligand of `11gl`, whose deposition title is "Crystal structure of
mSTING in complex with 2'3'-cUMP-AMP" (released 2026-07-29, no journal citation
yet). The workflow maps its interface contacts onto mouse residues 231, 237 and
262, whose human equivalents H232/R232, R238 and T263 are among the best
characterised CDN-coordinating residues in the STING literature. In other words,
this ligand contacts exactly the canonical cyclic-dinucleotide pocket.

**Interpretation: `ZNT` is a biologically meaningful agonist complex, not a
crystallisation artefact**, and `11gl`'s membership of the CDN-bound 9-member
state is coherent. Its omission from the demo notes appears to be an oversight
in the notes rather than a workflow problem. No change to the blocklist or to
production annotation behaviour is proposed or required.

A related observation: the ligands of `11gm`, `11gn` and `4loj` do not appear in
the interface annotation output at all, because the workflow reports only
ligands contacting *interface* residues. This is expected behaviour and not a
defect, but it means the ligand column understates which structures are
ligand-bound. Worth remembering when reading the cluster report.

---

## Prior systematic interface analysis of STING

Searches for prior work that systematically aggregates the STING dimer interface
across many deposited structures did not identify such a study. What was found
falls into three groups:

1. **Single-structure or few-structure comparisons.** The dominant pattern.
   Papers solve one to three structures and compare them with a small number of
   previously published ones, typically apo versus one ligand
   (Ouyang 2012, Chin KH 2013, Gao 2013, Zhang 2013, Chin EN 2020).
2. **Pairwise or small-set conformational comparisons using global metrics.**
   α2-tip separation, α1–α1 distance and interhelical angle, domain rotation
   angle. These classify structures as open or closed but do not report
   contact-level interface composition across many structures.
3. **Simulation-based studies.** Molecular dynamics of human STING with
   different CDNs, examining dynamical and conformational differences. These are
   trajectory analyses of one or two systems, not aggregations of deposited
   experimental structures.

No study was identified that aggregates all deposited experimental STING dimer
structures, computes residue-pair contact frequencies across them, and clusters
the structures by contact similarity. Given the full-text access limitation
recorded above, this is stated as **not identified in the literature searched**,
not as an absence of prior work. A proper claim of priority would require a
full-text review including supplementary material, and ideally a structured
search of the structural-bioinformatics methods literature rather than the
STING-biology literature.

---

## What is established versus workflow-derived

| Workflow observation | Literature status | Evidence | Interpretation |
|---|---|---|---|
| STING analysed as a homodimer, both roles Q3TBT3 | **Established** | Constitutive homodimer with an extensive hydrophobic interface, Ouyang 2012 and all subsequent structural work | Correct by construction; confirms the workflow's dimer handling on a real case |
| Multiple interface-contact states recovered from one complex's deposited structures | **Consistent with established biology** | STING is the textbook case of a ligand-induced open-to-closed transition; multiple conformational states are expected | Recovering more than one state is the biologically correct outcome; the *number* of states is a function of the chosen cut and is not a biological claim |
| Four states of sizes 9 / 3 / 1 / 1 at cut 0.6 | **Workflow-derived** | frozen fixture; the state count changes with the cut (9 states at 0.30, 4 from 0.55 to 0.70) | A parameter-dependent description of this deposition set, not a claim that STING has four conformational states |
| Cyclic dinucleotide and cGAMP-mimetic complexes share one interface state | **Consistent with established biology** | Gao 2013 (2'3'- and 3'3'-cGAMP closed); Chin EN 2020 (SR-717 induces the same closed conformation) | Ligand-class coherence across seven entries and four scaffolds; strong agreement with the published mechanism |
| That state is defined by a lid-region inter-protomer contact network | **Partly established, partly workflow-derived** | lid = human 219–249 from both protomers (Zhang 2013); two of the six core contacts are documented (below) | Contact-level expression of a known mechanism |
| G233–Y244 present in 9/9 CDN structures, absent elsewhere | **Established at pair level, workflow-derived at state level** | Zhang 2013 describes Tyr245 side chain with Gly234 main-chain carbonyl (human) | Pair-level agreement; the state-dependent frequency was not identified in the literature searched |
| S242–K235 present in 9/9 CDN structures and in `4lol` | **Established at pair level, workflow-derived at state level** | Zhang 2013 describes Ser243 with Lys236 main-chain amide (human) | As above; its presence in `4lol` is what makes the `4lol` contact set partial rather than absent |
| A232–D209 rewiring, 9/9 vs 0/3 | **Not identified in the literature searched** | residues in/adjacent to the lid; no pair-level description found | Structurally plausible as part of lid closure; requires validation |
| A232–Y260 rewiring, 9/9 vs 0/3 | **Not identified in the literature searched** | A232 in the lid, Y260 outside; beware human E260 vs mouse Y260 confusion | As above, with an explicit numbering caution |
| D209–G233 rewiring, 9/9 vs 0/3 | **Not identified in the literature searched** | G233 documented in a different lid pair | As above |
| D273–H156 rewiring, 0/9 vs 3/3 | **Requires qualification** | mouse D273 corresponds to human **Y274**; not conserved | Mouse-specific at best; must never be generalised to STING as a family |
| `4lol` (DMXAA) singleton | **Requires qualification** | Gao 2013 calls DMXAA-bound mSTING closed; the workflow separates it on contact identity | Partial lid engagement; a real difference in what is measured, not a contradiction |
| `4jc5` sparse singleton with QC warnings | **Requires qualification** | 4 pairs at 2.75 Å, the dataset's lowest resolution | Correctly flagged as uninterpretable; may be a detection artefact |
| DMXAA and SR-717 interface contacts map to mouse S161 / S161+S240 | **Consistent with established biology** | human S162A confers DMXAA sensitivity; residue 241 region is in the lid | Independent check that residue mapping lands on the right residues |
| `ZNT` at the `11gl` interface | **Established as chemistry, workflow-derived as annotation** | PDBe compound record: 2'3'-cUMP-AMP; contacts mouse 231/237/262 | Genuine CDN complex at the canonical pocket |
| Apo `4kc0` groups with c-di-GMP `4kby` | **Requires qualification** | Chin 2013 solved both; apo mouse STING reported as already closed | Cautions against labelling this state by ligand status; supports labelling states by contact composition instead |

---

## Implications for the interface-rewiring project

### What STING demonstrates as a case study

Stated conservatively, the STING example demonstrates four things.

1. **The workflow reproduces known conformational biology from deposited data
   alone.** Without any conformational metric, structural superposition or prior
   knowledge, clustering on interface contact similarity separates the cyclic
   dinucleotide and cGAMP-mimetic complexes from the apo and c-di-GMP
   structures. That separation matches the established ligand-induced lid
   closure mechanism, and it places SR-717 with the cyclic dinucleotides exactly
   as the Science 2020 co-crystal analysis concluded.

2. **It provides a contact-level representation of a state change that the
   literature usually describes with global metrics.** Published descriptions
   rely on α2-tip separation, interhelical angle or domain rotation. The
   workflow instead names the specific inter-protomer residue pairs gained and
   lost, with frequencies across the deposition set. Two of those pairs match
   interactions independently described in the primary literature, which is a
   meaningful external check on the representation.

3. **It exposes structures that global conformational labels obscure.** `4lol`
   is called closed in the literature yet carries only part of the CDN contact
   network. Whether that is biologically important or an artefact of ligand size
   and pocket occupancy is exactly the kind of question a contact-level
   representation raises and a global label does not.

4. **Its QC behaves as intended on real data.** The sparse, low-resolution
   `4jc5` singleton is flagged rather than presented as a fifth conformational
   state.

### Candidate observations for the manuscript

Labelled by evidence class, as required.

| Observation | Label |
|---|---|
| Clustering on interface contacts alone separates CDN/mimetic complexes from apo and c-di-GMP structures, matching the published lid-closure mechanism | **established biology reproduced** |
| SR-717 (`6xnn`) co-clusters with the cyclic dinucleotides, consistent with its published characterisation as a direct cGAMP mimetic | **established biology reproduced** |
| The state-defining contacts lie in and around the lid region (human 219–249), and two of them, G233–Y244 and S242–K235, correspond to inter-protomer lid interactions described in the primary literature | **mechanistic interpretation supported** |
| Ligand contacts map onto the canonical CDN pocket residues (mouse 231/237/240/262) and, for DMXAA, onto the known species-specificity residue S161 | **mechanistic interpretation supported** |
| State-dependent presence/absence of A232–D209, A232–Y260 and D209–G233 at 9/9 versus 0/3 | **workflow-derived observation** (not identified in the literature searched) |
| `4lol` (DMXAA) shows a partial lid-contact network, retaining S242–K235 while lacking the D209/A232/G233/Y260 network | **workflow-derived observation; requires further validation** |
| D273–H156 marks the 3-member state at 3/3 versus 0/9 | **workflow-derived observation; requires further validation**, and mouse-specific because human has Y274 |
| c-di-GMP (`4kby`) separates from all other cyclic dinucleotides in the dataset | **mechanistic interpretation supported**, consistent with published c-di-GMP versus cGAMP conformational differences |

### Recommended framing and cautions

If STING becomes a manuscript case study, the following should travel with it:

- Present it as **reproduction of known biology plus contact-level resolution**,
  not as discovery. The strongest claim the data supports is that an unbiased,
  structure-agnostic aggregation recovers a mechanism that took the field
  several dedicated crystallographic studies to establish.
- State the species. Every structure is mouse. Give both mouse and human
  numbering for any residue mentioned, and flag D273 as non-conserved.
- Avoid the words open and closed without defining which sense is meant, given
  the Chin 2013 usage and the report that apo mouse STING is already closed.
- Do not describe the four-state result as four conformational states of STING.
  It is a description of this deposition set at one cut height.
- Note that contact detection is resolution-dependent and that the workflow's
  own QC flags the affected structure.
- Before any publication claim of priority, complete a full-text review of
  Gao 2013 and Zhang 2013 including supplementary material, since interface
  contact tables are most likely to appear there.

Two concrete next steps would materially strengthen the case study, neither of
which is proposed for implementation now: computing a global conformational
metric such as α2-tip separation for each instance, to test directly whether the
workflow states align with the published open/closed axis and to quantify where
`4lol` sits; and running the same analysis on human STING complexes, where the
literature is deeper and D273 does not exist, to test which contacts transfer
across species.

---

## References

Primary structural papers, in the order they bear on this analysis.

1. Gao P, Ascano M, Zillinger T, Wang W, Dai P, Serganov AA, Gaffney BL, Shuman S, Jones RA, Deng L, Hartmann G, Barchet W, Tuschl T, Patel DJ. **Structure-function analysis of STING activation by c[G(2',5')pA(3',5')p] and targeting by antiviral DMXAA.** *Cell* 2013;154(4):748–762. DOI [10.1016/j.cell.2013.07.023](https://doi.org/10.1016/j.cell.2013.07.023). PMID 23910378. PDB 4LOJ, 4LOK, 4LOL. *(Full text not accessible during this search; publisher returned 403.)*

2. Zhang X, Shi H, Wu J, Zhang X, Sun L, Chen C, Chen ZJ. **Cyclic GMP-AMP containing mixed phosphodiester linkages is an endogenous high-affinity ligand for STING.** *Molecular Cell* 2013;51(2):226–235. PMC3808999. Source of the lid definition (human residues 219–249 from each protomer) and of the Tyr245–Gly234 and Ser243–Lys236 interdomain lid interactions. *(Full text not accessible during this search; PMC presented a bot check. Content obtained via search-engine extraction.)*

3. Chin KH, Tu ZL, Su YC, Yu YJ, Chen HC, Lo YC, Chen CP, Barber GN, Chuah MLC, Liang ZX, Chou SH. **Novel c-di-GMP recognition modes of the mouse innate immune adaptor protein STING.** *Acta Crystallographica Section D* 2013;69(Pt 3):352–366. DOI [10.1107/S0907444912047269](https://doi.org/10.1107/S0907444912047269). PMID 23519410. PDB 4KBY, 4KC0. Abstract retrieved in full.

4. Chin EN, Yu C, Vartabedian VF, Jia Y, Kumar M, Gamo AM, Vernier W, Ali SH, Kissai M, Lazar DC, Nguyen N, Pereira LE, Benish B, Woods AK, Joseph SB, Chu A, Johnson KA, Sander PN, Martínez-Peña F, Hampton EN, Young TS, Wolan DW, Chatterjee AK, Schultz PG, Petrassi HM, Teijaro JR, Lairson LL. **Antitumor activity of a systemic STING-activating non-nucleotide cGAMP mimetic.** *Science* 2020;369(6506):993–999. DOI [10.1126/science.abb4255](https://doi.org/10.1126/science.abb4255). PMID 32820126. PDB 6XNN (mouse), 6XNP (human).

5. Cavlar T, Deimling T, Ablasser A, Hopfner KP, Hornung V. **Species-specific detection of the antiviral small-molecule compound CMA by STING.** *EMBO Journal* 2013;32(10):1440–1450. DOI [10.1038/emboj.2013.86](https://doi.org/10.1038/emboj.2013.86). PMID 23604073. PDB 4JC5.

6. Ouyang S, Song X, Wang Y, Ru H, Shaw N, Jiang Y, Niu F, Zhu Y, Qiu W, Parvatiyar K, Li Y, Zhang R, Cheng G, Liu ZJ. **Structural analysis of the STING adaptor protein reveals a hydrophobic dimer interface and mode of cyclic di-GMP binding.** *Immunity* 2012;36(6):1073–1086. PMID 22579474. Source of the hydrophobic dimer interface description and buried-area figure.

7. Zhu D *et al.* **Structural insights into the distinct binding mode of cyclic di-AMP with SaCpaA_RCK.** *Biochemistry* 2015;54(31):4936–4951. DOI [10.1021/acs.biochem.5b00633](https://doi.org/10.1021/acs.biochem.5b00633). PMID 26171638. Primary citation attached to PDB 4YP1; note the paper's principal subject is a different protein.

8. Konno H *et al.* / Hopfner group. **Binding-pocket and lid-region substitutions render human STING sensitive to the species-specific drug DMXAA.** *Cell Reports* 2014. PMID 25199835. Source of the S162A / G230I / Q266I substitution analysis.

Supporting and contextual sources.

9. **Structural insights into STING signaling.** *Trends in Cell Biology* 2020 (review). Global conformational metrics, open/closed classification.

10. **Activation of STING based on its structural features.** *Frontiers in Immunology* 2022;13:808607. Lid-region residue summary and CDN-coordinating residues.

11. **Current understanding of the cGAS-STING signaling pathway: structure, regulatory mechanisms, and related diseases.** PMC9841179 (review). *(Full text not accessible; PMC bot check.)*

12. **Distinct oligomeric assemblies of STING induced by non-nucleotide agonists.** *Nature Communications* 2025;16. Agonist-class-dependent conformational behaviour, including diABZI retaining the open conformation.

13. **The mechanism of STING autoinhibition and activation.** *Molecular Cell* 2023. Ligand-class-dependent stabilisation of open versus closed LBD.

Database resources.

14. PDBe REST API, `https://www.ebi.ac.uk/pdbe/api` — entry summaries, publications, molecules and chemical component records for all twelve entries and ten chemical components. Accessed 2026-09-23.

15. UniProt Q3TBT3 (mouse STING, *Sting1*) and Q86WV6 (human STING, *STING1*/*TMEM173*), `https://rest.uniprot.org`. Sequences, cyclic-dinucleotide-binding domain boundaries (152–339 in mouse) and annotated binding sites. Accessed 2026-09-23.

---

## Summary of what changed and what did not

Nothing in the repository's behaviour was changed by this analysis. No
production module, no regression test, no fixture and no other specification
document was modified. The tranche-2A assertions remain valid regression
properties of the frozen analysis regardless of the literature status of any
individual contact, which is the separation this document was written to
preserve.
