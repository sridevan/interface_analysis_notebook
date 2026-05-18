"""InterfaceRecord and interaction-pair set construction.

Two interaction-pair representations per interface:
- author-keyed:  joins annotations within a single PDB entry
- uniprot-keyed: compares interfaces across PDB entries
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional

log = logging.getLogger(__name__)


def _norm_ins(value) -> Optional[str]:
    if value is None:
        return None
    s = str(value).strip()
    return s or None


AuthorResidueKey = tuple[str, str, int, Optional[str]]
UniProtResidueKey = tuple[str, int, int]
AuthorPair = tuple[AuthorResidueKey, AuthorResidueKey, str]
UniProtPair = tuple[UniProtResidueKey, UniProtResidueKey, str]


@dataclass
class InterfaceRecord:
    pdb_id: str
    assembly_id: str
    interface_id: int
    interface_info: dict
    unp_accession_1: str
    unp_accession_2: str
    author_pairs: set[AuthorPair] = field(default_factory=set)
    uniprot_pairs: set[UniProtPair] = field(default_factory=set)
    # Map author residue key -> uniprot key, used when joining annotations.
    author_to_uniprot: dict[AuthorResidueKey, UniProtResidueKey] = field(default_factory=dict)
    # UniProt residue key -> canonical one-letter amino acid code.
    residue_identity: dict[UniProtResidueKey, str] = field(default_factory=dict)
    n_residues_dropped_no_uniprot: int = 0
    n_microheterogeneity_collisions: int = 0

    @property
    def key(self) -> tuple[str, str, int]:
        return (self.pdb_id, self.assembly_id, self.interface_id)

    def label(self) -> str:
        return f"{self.pdb_id} a{self.assembly_id} i{self.interface_id}"


def build_interface_records(interface_response: list[dict]) -> list[InterfaceRecord]:
    """Build one InterfaceRecord per interface in the response.

    Drops residues lacking a UniProt mapping (with WARNING). Detects
    microheterogeneity (two distinct author residues mapping to the same
    UniProt position within one interface) and keeps the first-seen mapping
    (with WARNING).
    """
    records = []
    for it in interface_response:
        rec = _build_one(it)
        records.append(rec)
    return records


def valid_assemblies(
    complex_details: dict,
    max_resolution: float | None = None,
) -> set[tuple[str, str]]:
    """Return (pdb_id, assembly_id) tuples for valid assemblies.

    Excludes:
    - assemblies with bound macromolecules (a third macromolecule beyond
      the canonical dimer — co-chaperone, antibody Fab, accessory subunit
      — which can perturb the interface).
    - if `max_resolution` is set: assemblies whose resolution exceeds the
      threshold, OR whose resolution is None (NMR, unreported, etc.). The
      strict interpretation matches the typical user intent when filtering
      by resolution: keep only structures with explicit, sufficient
      resolution data.

    Default `max_resolution=None` disables the resolution filter (all
    resolutions, including None, are kept).
    """
    valid: set[tuple[str, str]] = set()
    for a in (complex_details.get("assemblies") or []):
        if a.get("bound_macromolecules"):
            continue
        if max_resolution is not None:
            res = a.get("resolution")
            if res is None:
                continue
            try:
                if float(res) > max_resolution:
                    continue
            except (TypeError, ValueError):
                continue
        pdb_id = a.get("pdb_id")
        assembly_id = a.get("assembly_id")
        if pdb_id and assembly_id is not None:
            valid.add((str(pdb_id), str(assembly_id)))
    return valid


def filter_interfaces_by_assembly(
    interface_response: list[dict],
    valid_assemblies: set[tuple[str, str]],
) -> tuple[list[dict], list[tuple[str, str]]]:
    """Drop interfaces whose (pdb_id, assembly_id) is not in valid_assemblies.

    Returns (filtered_response, dropped_keys).
    """
    filtered = []
    dropped: list[tuple[str, str]] = []
    for it in interface_response:
        key = (str(it.get("entry_id")), str(it.get("assembly_id")))
        if key in valid_assemblies:
            filtered.append(it)
        else:
            dropped.append(key)
    if dropped:
        log.warning(
            "Dropped %d interfaces from assemblies with bound macromolecules: %s",
            len(dropped), sorted(set(dropped)),
        )
    return filtered, dropped


def _build_one(item: dict) -> InterfaceRecord:
    pdb_id = item["entry_id"]
    assembly_id = str(item["assembly_id"])
    interface_id = int(item["interface_id"])
    interface_info = dict(item.get("interface_info") or {})
    interactions = item.get("interactions") or []

    # Determine canonical partner accessions for this interface.
    unp_acc_1 = unp_acc_2 = None
    for c in interactions:
        if c.get("unp_accession_1") and c.get("unp_accession_2"):
            unp_acc_1 = c["unp_accession_1"]
            unp_acc_2 = c["unp_accession_2"]
            break

    rec = InterfaceRecord(
        pdb_id=pdb_id,
        assembly_id=assembly_id,
        interface_id=interface_id,
        interface_info=interface_info,
        unp_accession_1=unp_acc_1 or "",
        unp_accession_2=unp_acc_2 or "",
    )

    dropped = 0
    # Track UniProt -> first author key seen, to detect microheterogeneity.
    uniprot_to_first_author: dict[UniProtResidueKey, AuthorResidueKey] = {}
    collisions = 0

    for c in interactions:
        a1 = c.get("auth_asym_id_1")
        a2 = c.get("auth_asym_id_2")
        r1 = c.get("auth_seq_id_1")
        r2 = c.get("auth_seq_id_2")
        i1 = _norm_ins(c.get("auth_ins_code_1") or c.get("ins_code_1"))
        i2 = _norm_ins(c.get("auth_ins_code_2") or c.get("ins_code_2"))
        bond_type = c.get("bond_type") or "unknown"

        u_acc_1 = c.get("unp_accession_1")
        u_acc_2 = c.get("unp_accession_2")
        u_pos_1 = c.get("unp_seq_id_1")
        u_pos_2 = c.get("unp_seq_id_2")

        if a1 is None or r1 is None or a2 is None or r2 is None:
            continue

        author_key_1: AuthorResidueKey = (pdb_id, a1, int(r1), i1)
        author_key_2: AuthorResidueKey = (pdb_id, a2, int(r2), i2)

        # The author-keyed pair is built regardless; this side is used for
        # annotation joins which do not require the UniProt mapping.
        rec.author_pairs.add((author_key_1, author_key_2, bond_type))

        # Cross-structure comparison requires UniProt on both sides.
        if not (u_acc_1 and u_acc_2 and u_pos_1 is not None and u_pos_2 is not None):
            dropped += 1
            continue

        u_key_1: UniProtResidueKey = (u_acc_1, int(u_pos_1), 1)
        u_key_2: UniProtResidueKey = (u_acc_2, int(u_pos_2), 2)

        # Capture one-letter codes for both residues; canonical per UniProt
        # position so consistent across structures.
        if c.get("unp_one_letter_code_1"):
            rec.residue_identity[u_key_1] = c["unp_one_letter_code_1"]
        if c.get("unp_one_letter_code_2"):
            rec.residue_identity[u_key_2] = c["unp_one_letter_code_2"]

        # Microheterogeneity check: two distinct author keys -> same UniProt key.
        for author_key, u_key in ((author_key_1, u_key_1), (author_key_2, u_key_2)):
            prior = uniprot_to_first_author.get(u_key)
            if prior is None:
                uniprot_to_first_author[u_key] = author_key
                rec.author_to_uniprot[author_key] = u_key
            elif prior != author_key:
                collisions += 1
                # Keep first-seen mapping. This author key is not added to
                # author_to_uniprot, so any annotations on it will not be
                # joined to a UniProt key (intentional).

        rec.uniprot_pairs.add((u_key_1, u_key_2, bond_type))

    if dropped:
        log.warning(
            "Dropped %d atom-level contacts at interface %s/%s/%d for missing UniProt mapping",
            dropped, pdb_id, assembly_id, interface_id,
        )
    if collisions:
        log.warning(
            "Microheterogeneity at interface %s/%s/%d: %d author residues collided "
            "with a prior UniProt key; kept first-seen mapping.",
            pdb_id, assembly_id, interface_id, collisions,
        )

    rec.n_residues_dropped_no_uniprot = dropped
    rec.n_microheterogeneity_collisions = collisions
    return rec


def check_partner_consistency(records: list[InterfaceRecord]) -> list[InterfaceRecord]:
    """Verify (unp_accession_1, unp_accession_2) is consistent across interfaces.

    The majority ordering is taken as canonical. Inconsistent entries have
    their role assignment reversed in place: roles are swapped in author_pairs
    (only relevant for homodimer display), uniprot_pairs, and the partner
    accessions on the record.
    """
    if not records:
        return records

    pair_counts: Counter = Counter(
        (r.unp_accession_1, r.unp_accession_2) for r in records if r.unp_accession_1
    )
    if not pair_counts:
        return records

    canonical, _ = pair_counts.most_common(1)[0]

    for r in records:
        current = (r.unp_accession_1, r.unp_accession_2)
        if current == canonical:
            continue
        if current == (canonical[1], canonical[0]):
            log.warning(
                "Partner ordering inconsistent at %s: %s vs canonical %s. "
                "Reversing role assignment.",
                r.label(), current, canonical,
            )
            _reverse_roles(r)
        else:
            log.warning(
                "Partner accessions at %s (%s) do not match canonical %s in "
                "either order; leaving as-is.",
                r.label(), current, canonical,
            )
    return records


def _reverse_roles(r: InterfaceRecord) -> None:
    r.unp_accession_1, r.unp_accession_2 = r.unp_accession_2, r.unp_accession_1
    r.author_pairs = {(b, a, bt) for (a, b, bt) in r.author_pairs}
    new_uniprot = set()
    for (k1, k2, bt) in r.uniprot_pairs:
        # Swap role labels so the canonical-1 partner stays role 1.
        new_k1 = (k1[0], k1[1], 2)
        new_k2 = (k2[0], k2[1], 1)
        new_uniprot.add((new_k2, new_k1, bt))
    r.uniprot_pairs = new_uniprot
    new_a2u = {}
    for author_key, u_key in r.author_to_uniprot.items():
        new_role = 2 if u_key[2] == 1 else 1
        new_a2u[author_key] = (u_key[0], u_key[1], new_role)
    r.author_to_uniprot = new_a2u
    new_identity = {}
    for u_key, aa in r.residue_identity.items():
        new_role = 2 if u_key[2] == 1 else 1
        new_identity[(u_key[0], u_key[1], new_role)] = aa
    r.residue_identity = new_identity
