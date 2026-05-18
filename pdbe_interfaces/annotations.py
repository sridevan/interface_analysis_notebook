"""Annotation enrichment: mutations, modifications, ligands.

Joins each annotation stream onto interface residues by author key. Annotations
falling outside the interface are dropped (this is an interface-comparison
workflow). After joining, the UniProt key of the matched interface residue is
attached so cross-structure comparison is straightforward.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

from pdbe_interfaces import api
from pdbe_interfaces.representation import (
    AuthorResidueKey,
    InterfaceRecord,
    UniProtResidueKey,
    _norm_ins,
)

log = logging.getLogger(__name__)


@dataclass
class LigandContact:
    """A single (interface residue, ligand instance, contact_type) tuple."""

    interface_residue: AuthorResidueKey
    ligand_chem_comp_id: str
    ligand_chain_id: str
    ligand_author_residue_number: int
    contact_type: str


@dataclass
class AnnotationOverlap:
    mutations: pd.DataFrame = field(default_factory=pd.DataFrame)
    modifications: pd.DataFrame = field(default_factory=pd.DataFrame)
    ligands: pd.DataFrame = field(default_factory=pd.DataFrame)


def filter_bound_molecules(
    bm_records: list[dict],
    blocklist: frozenset[str],
    drop_carbohydrate_polymers: bool = False,
) -> list[dict]:
    """Return surviving ligand instances flattened from bm.composition.ligands.

    Each output dict has keys: chain_id, author_residue_number, ins_code,
    chem_comp_id. Bound molecules with no surviving ligands contribute nothing.
    """
    surviving: list[dict] = []
    dropped_chem_comps: list[str] = []
    for bm in bm_records:
        ligands = (bm.get("composition") or {}).get("ligands") or []
        for lig in ligands:
            ccid = lig.get("chem_comp_id")
            mtype = lig.get("molecule_type")
            if drop_carbohydrate_polymers and mtype == "Carbohydrate-polymer":
                dropped_chem_comps.append(f"{ccid}[carb]")
                continue
            if ccid in blocklist:
                dropped_chem_comps.append(ccid)
                continue
            if not ccid or lig.get("chain_id") is None or lig.get("author_residue_number") is None:
                continue
            surviving.append({
                "chain_id": lig["chain_id"],
                "author_residue_number": int(lig["author_residue_number"]),
                "ins_code": _norm_ins(lig.get("author_insertion_code")),
                "chem_comp_id": ccid,
            })
    if dropped_chem_comps:
        log.info("Filtered out %d ligand instances: %s",
                 len(dropped_chem_comps), sorted(set(dropped_chem_comps)))
    return surviving


def collect_ligand_contacts(
    pdb_id: str, surviving_ligands: list[dict],
) -> list[LigandContact]:
    """For each surviving ligand instance, fetch interactions and expand.

    Atom-level records with multiple `interaction_details` labels are expanded
    so each label produces a separate element. Atom-level detail is then
    collapsed: multiple atom-level records sharing the same
    (interface residue, ligand instance, contact_type) key aggregate to one.
    """
    contacts: set[tuple] = set()
    chem_comp_by_instance: dict[tuple, str] = {}
    for lig in surviving_ligands:
        bm_objects = api.fetch_ligand_interactions(
            pdb_id, lig["chain_id"], lig["author_residue_number"],
        )
        for bm in bm_objects:
            lig_meta = bm.get("ligand") or {}
            lig_chem = lig_meta.get("chem_comp_id") or lig["chem_comp_id"]
            lig_chain = lig_meta.get("chain_id") or lig["chain_id"]
            lig_resnum = int(lig_meta.get("author_residue_number") or lig["author_residue_number"])
            instance_key = (lig_chem, lig_chain, lig_resnum)
            chem_comp_by_instance[instance_key] = lig_chem
            for it in bm.get("interactions") or []:
                end = it.get("end") or {}
                if end.get("chain_id") is None or end.get("author_residue_number") is None:
                    continue
                residue_key: AuthorResidueKey = (
                    pdb_id,
                    end["chain_id"],
                    int(end["author_residue_number"]),
                    _norm_ins(end.get("author_insertion_code")),
                )
                for ctype in it.get("interaction_details") or []:
                    contacts.add((residue_key, instance_key, ctype))
    return [
        LigandContact(
            interface_residue=residue_key,
            ligand_chem_comp_id=instance_key[0],
            ligand_chain_id=instance_key[1],
            ligand_author_residue_number=instance_key[2],
            contact_type=ctype,
        )
        for (residue_key, instance_key, ctype) in contacts
    ]


def overlap_annotations(
    records: list[InterfaceRecord],
    mutations_response: dict,
    modifications_response: dict,
    ligand_contacts_by_pdb: dict[str, list[LigandContact]],
    mutation_type_filter: tuple[str, ...] = ("Engineered mutation",),
) -> AnnotationOverlap:
    """Join annotations onto interface residues by author key.

    Annotations falling outside the interface are dropped. The UniProt key
    of the matched interface residue is attached.
    """
    interfaces_by_pdb: dict[str, list[InterfaceRecord]] = defaultdict(list)
    for r in records:
        interfaces_by_pdb[r.pdb_id].append(r)

    mut_rows = []
    for pdb_id, mut_records in (mutations_response or {}).items():
        for m in mut_records or []:
            if (m.get("mutation_details") or {}).get("type") not in mutation_type_filter:
                continue
            chain_id = m.get("chain_id")
            resnum = m.get("author_residue_number")
            if chain_id is None or resnum is None:
                continue
            ins = _norm_ins(m.get("author_insertion_code"))
            author_key = (pdb_id, chain_id, int(resnum), ins)
            mdet = m.get("mutation_details") or {}
            label = f"{mdet.get('from') or '?'}{resnum}{mdet.get('to') or '?'}"
            for r in interfaces_by_pdb.get(pdb_id, []):
                u_key = r.author_to_uniprot.get(author_key)
                if u_key is None:
                    # Annotation outside interface, or on a residue without UniProt mapping.
                    if author_key in {a for (a, _, _) in r.author_pairs} | {b for (_, b, _) in r.author_pairs}:
                        # Author residue *is* in the interface but lost its UniProt mapping
                        u_key = (None, None, None)
                    else:
                        continue
                mut_rows.append({
                    "pdb_id": pdb_id,
                    "assembly_id": r.assembly_id,
                    "interface_id": r.interface_id,
                    "auth_asym_id": chain_id,
                    "auth_seq_id": int(resnum),
                    "ins_code": ins,
                    "unp_acc": u_key[0],
                    "unp_seq_id": u_key[1],
                    "role": u_key[2],
                    "mutation_label": label,
                    "mutation_from": mdet.get("from"),
                    "mutation_to": mdet.get("to"),
                    "mutation_type": mdet.get("type"),
                })

    mod_rows = []
    for pdb_id, mod_records in (modifications_response or {}).items():
        for m in mod_records or []:
            chain_id = m.get("chain_id")
            resnum = m.get("author_residue_number")
            if chain_id is None or resnum is None:
                continue
            ins = _norm_ins(m.get("author_insertion_code"))
            author_key = (pdb_id, chain_id, int(resnum), ins)
            for r in interfaces_by_pdb.get(pdb_id, []):
                u_key = r.author_to_uniprot.get(author_key)
                if u_key is None:
                    if author_key in {a for (a, _, _) in r.author_pairs} | {b for (_, b, _) in r.author_pairs}:
                        u_key = (None, None, None)
                    else:
                        continue
                mod_rows.append({
                    "pdb_id": pdb_id,
                    "assembly_id": r.assembly_id,
                    "interface_id": r.interface_id,
                    "auth_asym_id": chain_id,
                    "auth_seq_id": int(resnum),
                    "ins_code": ins,
                    "unp_acc": u_key[0],
                    "unp_seq_id": u_key[1],
                    "role": u_key[2],
                    "modification_chem_comp_id": m.get("chem_comp_id"),
                    "modification_chem_comp_name": m.get("chem_comp_name"),
                })

    lig_rows = []
    for pdb_id, contacts in (ligand_contacts_by_pdb or {}).items():
        for c in contacts:
            author_key = c.interface_residue
            for r in interfaces_by_pdb.get(pdb_id, []):
                u_key = r.author_to_uniprot.get(author_key)
                if u_key is None:
                    if author_key in {a for (a, _, _) in r.author_pairs} | {b for (_, b, _) in r.author_pairs}:
                        u_key = (None, None, None)
                    else:
                        continue
                lig_rows.append({
                    "pdb_id": pdb_id,
                    "assembly_id": r.assembly_id,
                    "interface_id": r.interface_id,
                    "auth_asym_id": author_key[1],
                    "auth_seq_id": author_key[2],
                    "ins_code": author_key[3],
                    "unp_acc": u_key[0],
                    "unp_seq_id": u_key[1],
                    "role": u_key[2],
                    "ligand_chem_comp_id": c.ligand_chem_comp_id,
                    "ligand_chain_id": c.ligand_chain_id,
                    "ligand_author_residue_number": c.ligand_author_residue_number,
                    "contact_type": c.contact_type,
                })

    return AnnotationOverlap(
        mutations=pd.DataFrame(mut_rows),
        modifications=pd.DataFrame(mod_rows),
        ligands=pd.DataFrame(lig_rows),
    )
