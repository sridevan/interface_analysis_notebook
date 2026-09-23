# Tranche 2B: frozen Spike RBD–ACE2 heterodimer regression

Completed 2026-09-23. Second scientific end-to-end regression: a deliberately
small frozen subset of the SARS-CoV-2 Spike RBD with ACE2 heterodimer,
`6m0j` / `PDB-CPX-140195`, run offline through the same production functions the
notebook calls.

> **STING tranche 2A tested the full scientific state-analysis path on a
> homodimer; Spike/ACE2 tranche 2B tests heterodimer partner correspondence and
> annotation integrity.**

On a homodimer the partner-consistency check has nothing to act on, because both
partners share an accession. Only a heterodimer can exercise the path this
tranche protects:

```text
heterodimer partner order differs between structural instances
            ↓
partner normalisation
            ↓
equivalent contacts become comparable
            ↓
annotations remain attached to the correct biological partner
```

---

## Why this example was selected

The full aggregation holds 131 comparable interface instances across 115 PDB
entries, of which **36 are reported by PISA with the two partners in the
opposite order**. That is by far the richest source of partner-reversal cases
among the worked examples, and STING cannot produce any.

Three entries carry **both orientations across their own assemblies**: `7p19`,
`7rpv` and `8xye`. These are the cleanest possible test material, because two
assemblies of one deposition share construct, crystal form and refinement, so
they represent the same biological interface with the partners reported in
opposite orders. That is a strong control for partner-order reversal, while
still allowing genuine differences in the contacts detected in each assembly
copy.

`7rpv` was chosen as the primary case: assemblies 1–3 are canonical, assembly 4
is reversed, and all four carry the same engineered ACE2 mutation. `7p19` was
added because its interface mutation is on the *other* partner, so the fixture
covers annotation identity on both sides of the heterodimer.

---

## Fixture provenance

Location: `tests/fixtures/recorded/spike_ace2_reversal/` — 200 KB, eight files
plus metadata.

| File | Size | Endpoint |
|---|---|---|
| `complex_details_6m0j.json` | 0.8 KB | `complex/details/6m0j?id_type=pdb_id` |
| `complex_details.json` | 59 KB | `complex/details/PDB-CPX-140195?id_type=pdb_complex_id` (verbatim, all assemblies) |
| `interface_interactions.json` | 38 KB | `complex/interface_interactions/PDB-CPX-140195` (**subset**, see below) |
| `mutations_post.json` | 13 KB | POST `pdb/entry/mutated_AA_or_NA` for the three selected entries |
| `modifications_post.json` | 0.2 KB | POST `pdb/entry/modified_AA_or_NA`, recorded 404 |
| `bound_molecules.json` | 6.6 KB | `pdb/bound_molecules/{pdb_id}` for the three entries |
| `ligand_interactions.json` | 62 KB | `pdb/bound_ligand_interactions/...` for the 14 surviving ligand instances |
| `metadata.json` | 1.7 KB | complex id, starting entry, selected instances, capture date, rationale |

Captured once on 2026-09-23 by running the notebook's Phase 1 call sequence
against the live API with the default configuration while
`requests.Session.request` was wrapped to record every request and response;
22 responses were recorded. The fixture is frozen and is not refreshed
automatically.

### The one deliberate deviation from raw capture

`interface_interactions.json` stores a **subset** of the real response body:
the five selected interface instances, each record verbatim, with every other
field of the response unchanged. The full response carries 152 interfaces and
would have made the fixture roughly ten times larger for no additional
regression value, since this tranche is not a clustering study. The subsetting
is recorded in the file itself (`_note`) and in `metadata.json`.

`complex_details.json` is kept **verbatim**, with all 151 assemblies. Assembly
filtering therefore runs against the real metadata rather than a trimmed copy,
and the annotation POST bodies match exactly what the offline run requests.

---

## Selected structural instances

| Instance | PISA partner order | Role after normalisation | Why selected |
|---|---|---|---|
| `6m0j` a1 i1 | canonical | ACE2 role 1 | the archetypal Spike RBD–ACE2 structure and the named starting entry |
| `7p19` a1 i1 | canonical | ACE2 role 1 | carries the Spike Q498Y engineered mutation at the interface (author chain E) |
| `7p19` a2 i1 | **reversed** | ACE2 role 1 | same entry as `7p19` a1, opposite reported order |
| `7rpv` a1 i1 | canonical | ACE2 role 1 | same entry and same biological ACE2-RBD interface as the reversed assembly; a strong control for partner-order reversal, while allowing genuine differences in the contacts detected in each copy (`7rpv` a4) |
| `7rpv` a4 i1 | **reversed** | ACE2 role 1 | carries the ACE2 Q325Y engineered mutation at the interface (author chain D) |

Three canonical against two reversed gives `check_partner_consistency` an
unambiguous majority. A 2:2 split would have made the test depend on the
first-seen tie-break rather than on the majority rule, which is why five
instances were taken rather than four.

---

## Canonical partner mapping

Established from PDBe entity metadata and UniProt, **not** from the workflow's
output, then asserted:

| Role | Biological partner | UniProt | Organism |
|---|---|---|---|
| 1 | Angiotensin-converting enzyme 2 (ACE2) | `Q9BYF1` | *Homo sapiens* |
| 2 | Spike glycoprotein (RBD) | `P0DTC2` | SARS-CoV-2 |

Author chain assignment, from the PDBe `molecules` endpoint:

```text
6m0j   chain A = ACE2      chain E = Spike
7p19   chains A, B = ACE2  chains C, E = Spike
7rpv   chains A-D = ACE2   chains E-H = Spike
```

So in the reversed instances, `7rpv` a4 uses chains D (ACE2) and H (Spike), and
`7p19` a2 uses chains B (ACE2) and C (Spike). The workflow's `author_to_uniprot`
reproduces exactly this assignment after normalisation.

---

## Representative equivalent contacts

Three well-established ACE2–RBD interface contacts, confirmed present in the
frozen data before being encoded:

| Contact (normalised: ACE2 role 1, Spike role 2) | Bond |
|---|---|
| ACE2 Asp30 – Spike Lys417 | salt bridge |
| ACE2 Asp38 – Spike Gln498 | hydrogen bond |
| ACE2 Tyr41 – Spike Thr500 | hydrogen bond |

For each, the regression asserts that **before** normalisation the canonical
instance carries the contact in ACE2→Spike order while the reversed instance
carries only its flipped form, and that **after** normalisation both instances
express it identically. Residue identities are also asserted (`D30`, `Y41`,
`K417`, `T500`), which confirms the mapping resolves to the right residues and
not merely to the right numbers.

`7rpv` a1 and a4 share 11 contacts once comparable.

---

## Annotation mapping tested

| Entry / assembly | Orientation | Author chain, residue | Resolves to | Mutation |
|---|---|---|---|---|
| `7rpv` a4 | **reversed** | D 325 | `Q9BYF1` 325, role 1 (ACE2) | Q325Y |
| `7rpv` a1 | canonical | A 325 | `Q9BYF1` 325, role 1 (ACE2) | Q325Y |
| `7p19` a1 | canonical | E 498 | `P0DTC2` 498, role 2 (Spike) | Q498Y |

The central assertion is the first row: an engineered mutation on a **reversed**
instance, sitting on author chain D, must resolve to ACE2 Gln325 in role 1 and
must not be attributed to Spike residue 325. The canonical assembly of the same
entry carries the same mutation on a different author chain and resolves to the
same biological identity, which is the control. The `7p19` row shows the other
partner mapping to role 2.

The canonical UniProt residue identity at `(Q9BYF1, 325, 1)` is `Q`, consistent
with the `Q325Y` label.

---

## Scientific expectations protected

Five tests in `tests/e2e/test_spike_ace2_recorded.py`:

1. **`normalises_partner_orientation`** — all 22 requests served from the
   fixture; complex resolves to `PDB-CPX-140195` as a heterodimer; canonical
   role mapping asserted; five instances from three entries with nothing
   excluded; two instances reported reversed before normalisation and all five
   carrying ACE2 in role 1 after.
2. **`reverses_every_view_of_a_record_together`** — for `7rpv` a4: author
   chains map to the correct partners, residue identities follow the new roles
   and agree with the canonical assembly, author-space pairs stay oriented
   ACE2-chain-first, every UniProt contact runs role 1 → role 2, the contact
   count is unchanged by the reversal, and the tranche-1.5 per-residue mappings
   remain self-consistent.
3. **`preserves_equivalent_contacts_after_reversal`** — the three contacts
   above, before and after, plus the 11 shared contacts.
4. **`similarity_is_not_broken_by_partner_order`** — Jaccard between the
   canonical and reversed assemblies of one entry is **0.0 before** and rises
   after normalisation; no instance in the cohort is left with zero similarity
   to every other.
5. **`preserves_annotation_partner_identity`** — the annotation table above.

### On not forcing Jaccard 1.0

After normalisation the equivalent assemblies reach **0.5500** (`7rpv` a1 vs a4)
and **0.6923** (`7p19` a1 vs a2), not 1.0. The assemblies genuinely differ in a
few detected contacts, as independent copies of a crystallographic assembly
normally do. The regression asserts the real values and the direction of change
rather than an idealised identity, which the instructions for this tranche
explicitly called for.

The before-value of exactly 0.0 is the scientifically important number: reported
in opposite orders, two views of the same biological interface share **no**
contact tuple at all. Without normalisation they would cluster apart purely
because of a reporting convention.

---

## What this regression intentionally does NOT test

- **Clustering or conformational states.** No dendrogram, cut-height, state-size
  or rewiring assertions. The Spike/ACE2 conformational landscape, variant
  structural biology and antibody complexes are all out of scope.
- **The full aggregation.** 5 instances of 131; the other 126 and the 21
  antibody-bearing assemblies excluded by the assembly filter are not covered.
- **Ligand annotation.** No ligand reaches any of these five interfaces, so the
  ligand path is exercised (14 instances fetched and filtered) but produces an
  empty overlap, which is asserted as such.
- **Modifications.** The endpoint returns 404 for these entries; the empty path
  is asserted, the populated path is not.
- **Literature grounding.** Deliberately deferred. Enough authoritative metadata
  was recorded to establish the biological identities (UniProt accessions,
  organism, chain-to-partner assignment); no review of the Spike/ACE2 structural
  literature was performed.

---

## Results

```text
tests before      92
new e2e tests      5
tests after       97

passed            97
failed             0
skipped            0
xfail              0

runtime           2.1 s for the suite; the Spike/ACE2 tests take about 0.05 s
```

Fully offline. Every URL is served by the recorded router, which raises on any
unrecorded request, and the root `block_network` guard remains active elsewhere.
The first test asserts that all 22 calls went through the router.

---

## Production changes

None. `pdbe_interfaces/` and `notebook.ipynb` were not touched in this tranche.

One small **test-code** change was needed: `tests/e2e/conftest.py` previously
assumed the fixture directory was named after the complex id and that the
entry-lookup file was `complex_details_11gl.json`. It now takes a fixture
directory name, reads the complex id from `metadata.json`, and discovers the
entry-lookup file by glob. This is test plumbing only and changes no behaviour
for the STING tests, which continue to pass unchanged.

---

## Observations from the real data

- **`7p19` assembly 2 does not contribute a mutation row.** The Q498Y record
  exists for author chain C, which is the Spike chain in that assembly, but
  residue 498 is not an interface residue there, so the workflow correctly
  reports nothing. The equivalent residue *is* at the interface in assembly 1.
  This is a real difference between two assemblies of one deposition, not a
  mapping failure, and it is why the Spike-side annotation assertion uses the
  canonical assembly.
- **`7rpv` carries five engineered mutations on ACE2** (L79F, M82Y, Q325Y,
  H374A, H378A); only Q325Y lies at the interface, so only it is reported. The
  H374A/H378A pair are the catalytic-site substitutions used to render ACE2
  inactive, which sit far from the RBD interface, and the filtering behaves
  accordingly.
- **The shared contact core is the textbook ACE2–RBD interface.** The 11
  contacts shared between `7rpv` a1 and a4 involve ACE2 residues 30, 35, 38, 41,
  42 and Spike residues 417, 493, 498, 500 among others, which are the
  canonical ACE2–RBD contact residues. No literature review was performed for
  this tranche, but the correspondence is worth noting as an informal sanity
  check on the mapping.
