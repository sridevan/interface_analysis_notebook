# Copyright 2026 EMBL - European Bioinformatics Institute
# Author: Sri Devan Appasamy
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Structure table, conserved sets, cluster interpretation report."""

from __future__ import annotations

import json
import statistics
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from pdbe_interfaces.annotations import AnnotationOverlap
from pdbe_interfaces.representation import InterfaceRecord
from pdbe_interfaces.similarity import ClusterResult


# QC thresholds. Tuned against the working examples (Spike–ACE2, KRAS–RAF1,
# CDK2–CycA2). Override by passing kwargs to cluster_interpretation_report.
DEFAULT_SPARSE_PAIR_THRESHOLD = 5         # median residue-pair count per interface
DEFAULT_TINY_AREA_THRESHOLD_A2 = 500.0    # median interface area (Å²)
DEFAULT_RANGE_OVERLAP_MIN = 0.2           # fraction of range overlap with dominant cluster


def extract_assembly_metadata(complex_details: dict) -> dict[tuple[str, str], dict]:
    """Per-assembly experimental_method and resolution from complex/details.

    Returned dict is keyed on (pdb_id, assembly_id) (both as strings) and maps
    to a small dict with `experimental_method` (str|None) and `resolution`
    (float|None). Used by build_structure_table and cluster_interpretation_report
    to surface methodological context alongside biological annotations.
    """
    out: dict[tuple[str, str], dict] = {}
    for a in (complex_details.get("assemblies") or []):
        pdb = a.get("pdb_id")
        asm = a.get("assembly_id")
        if pdb is None or asm is None:
            continue
        res = a.get("resolution")
        out[(str(pdb), str(asm))] = {
            "experimental_method": a.get("experimental_method"),
            "resolution": float(res) if res is not None else None,
        }
    return out


def build_partner_map(complex_details: dict) -> dict[tuple[str, int], str]:
    """Build (unp_accession, role) -> gene name from complex/details.

    For homodimers (single participant with stoichiometry 2) both roles share
    an accession; the gene name has " (role 1)" / " (role 2)" appended so the
    two sides are distinguishable in plots and the structure table.
    """
    participants = complex_details.get("participants") or []
    proteins = [p for p in participants if p.get("accession_type") == "UniProt"]
    pmap: dict[tuple[str, int], str] = {}
    if len(proteins) == 1 and (proteins[0].get("stoichiometry") or 1) == 2:
        acc = proteins[0]["accession"]
        gene = proteins[0].get("gene_name") or proteins[0].get("name") or acc
        pmap[(acc, 1)] = f"{gene} (role 1)"
        pmap[(acc, 2)] = f"{gene} (role 2)"
    else:
        for i, p in enumerate(proteins[:2], start=1):
            acc = p["accession"]
            gene = p.get("gene_name") or p.get("name") or acc
            pmap[(acc, i)] = gene
    return pmap


def label_for_record(
    record: InterfaceRecord,
    overlap: Optional[AnnotationOverlap] = None,
) -> str:
    """`pdb_id` plus a short annotation suffix when present.

    e.g. `6vjj [G12C, sotorasib]`. Truncated to ~30 characters.
    """
    parts: list[str] = []
    if overlap is not None and not overlap.mutations.empty:
        muts = overlap.mutations[
            (overlap.mutations["pdb_id"] == record.pdb_id)
            & (overlap.mutations["assembly_id"] == record.assembly_id)
            & (overlap.mutations["interface_id"] == record.interface_id)
        ]["mutation_label"].unique().tolist()
        parts.extend(muts[:2])
    if overlap is not None and not overlap.ligands.empty:
        ligs = overlap.ligands[
            (overlap.ligands["pdb_id"] == record.pdb_id)
            & (overlap.ligands["assembly_id"] == record.assembly_id)
            & (overlap.ligands["interface_id"] == record.interface_id)
        ]["ligand_chem_comp_id"].unique().tolist()
        parts.extend(ligs[:2])
    label = record.pdb_id
    if parts:
        suffix = ", ".join(parts)
        if len(suffix) > 22:
            suffix = suffix[:19] + "…"
        label = f"{record.pdb_id} [{suffix}]"
    return label


def describe_complex(
    complex_id: str,
    details: dict,
    partner_map: dict,
    assembly_metadata: dict,
    identifier: str | None = None,
) -> str:
    """Human-readable summary of the resolved complex."""
    lines = []
    if identifier and identifier.strip().upper() != complex_id.upper():
        lines.append(f"Resolved {identifier} to {complex_id}")
    lines += [
        f"Complex: {details.get('name')} ({complex_id})",
        f"Oligomeric state: {details.get('oligomeric_state')}",
        f"Partners: {partner_map}",
        f"Assembly metadata for {len(assembly_metadata)} (pdb_id, assembly_id) pairs",
    ]
    return "\n".join(lines)


def interface_summary_table(records: list[InterfaceRecord]) -> pd.DataFrame:
    """One row per interface: pair counts and residues lost in mapping."""
    return pd.DataFrame([
        {
            "pdb_id": r.pdb_id,
            "assembly_id": r.assembly_id,
            "interface_id": r.interface_id,
            "n_author_pairs": len(r.author_pairs),
            "n_uniprot_pairs": len(r.uniprot_pairs),
            "n_dropped_no_uniprot": r.n_residues_dropped_no_uniprot,
            "n_microheterogeneity": r.n_microheterogeneity_collisions,
        }
        for r in records
    ])


def build_structure_table(
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    overlap: AnnotationOverlap,
    partner_map: dict[tuple[str, int], str],
    fetch_date: str | None = None,
    assembly_metadata: dict[tuple[str, str], dict] | None = None,
) -> pd.DataFrame:
    fetch_date = fetch_date or date.today().isoformat()
    rows = []
    metadata = assembly_metadata or {}
    for i, r in enumerate(records):
        info = r.interface_info
        partner_1 = partner_map.get((r.unp_accession_1, 1), r.unp_accession_1 or "?")
        partner_2 = partner_map.get((r.unp_accession_2, 2), r.unp_accession_2 or "?")
        cid = (
            int(cluster_result.flat_assignment[i])
            if i < len(cluster_result.flat_assignment)
            else 0
        )
        meta = metadata.get((r.pdb_id, r.assembly_id), {})

        n_muts = 0
        if not overlap.mutations.empty:
            n_muts = int(((overlap.mutations["pdb_id"] == r.pdb_id)
                          & (overlap.mutations["assembly_id"] == r.assembly_id)
                          & (overlap.mutations["interface_id"] == r.interface_id)).sum())
        n_mods = 0
        if not overlap.modifications.empty:
            n_mods = int(((overlap.modifications["pdb_id"] == r.pdb_id)
                          & (overlap.modifications["assembly_id"] == r.assembly_id)
                          & (overlap.modifications["interface_id"] == r.interface_id)).sum())
        n_lig_residues = 0
        n_distinct_ligands = 0
        if not overlap.ligands.empty:
            sub = overlap.ligands[
                (overlap.ligands["pdb_id"] == r.pdb_id)
                & (overlap.ligands["assembly_id"] == r.assembly_id)
                & (overlap.ligands["interface_id"] == r.interface_id)
            ]
            n_lig_residues = int(sub.drop_duplicates(
                subset=["auth_asym_id", "auth_seq_id", "ins_code"]
            ).shape[0])
            n_distinct_ligands = int(sub["ligand_chem_comp_id"].nunique())

        rows.append({
            "pdb_id": r.pdb_id,
            "assembly_id": r.assembly_id,
            "interface_id": r.interface_id,
            "partner_1": partner_1,
            "partner_2": partner_2,
            "experimental_method": meta.get("experimental_method"),
            "resolution": meta.get("resolution"),
            "interface_area": info.get("interface_area"),
            "solvation_energy": info.get("solvation_energy"),
            "stabilisation_energy": info.get("stabilization_energy"),
            "pisa_p_value": info.get("p_value"),
            "n_interface_residues": info.get("number_interface_residues"),
            "n_hydrogen_bonds": info.get("number_hydrogen_bonds"),
            "n_salt_bridges": info.get("number_salt_bridges"),
            "n_covalent_bonds": info.get("number_covalent_bonds"),
            "n_disulfide_bonds": info.get("number_disulfide_bonds"),
            "n_engineered_mutations_at_interface": n_muts,
            "n_modifications_at_interface": n_mods,
            "n_distinct_ligands_at_interface": n_distinct_ligands,
            "n_ligand_contact_residues_at_interface": n_lig_residues,
            "cluster_id": cid,
            "n_residues_dropped_no_uniprot": r.n_residues_dropped_no_uniprot,
            "fetch_date": fetch_date,
        })
    return pd.DataFrame(rows)


def interface_frequency_summary(
    records: list[InterfaceRecord],
    partner_map: dict[tuple[str, int], str] | None = None,
    residue_identity: dict | None = None,
    min_residue_frequency: int = 1,
) -> dict:
    """Build per-pair and per-residue frequency tables across all interfaces.

    Returns a dict with:
    - `pairs`: DataFrame, one row per UniProt-keyed residue-residue pair
      (collapsed across bond_type), columns: contact_label, partner_1,
      partner_2, n_interfaces, fraction. Sorted by fraction descending.
    - `partner_1_residues`: DataFrame, one row per partner-1 residue, columns:
      residue_label, accession, position, n_interfaces, fraction. Sorted desc.
    - `partner_2_residues`: same for partner-2.
    - `pair_matrix`: numpy 2D array, rows = partner-1 residues, cols =
      partner-2 residues, values = count of interfaces containing this pair.
    - `partner_1_labels`, `partner_2_labels`: list[str] of axis labels for
      the matrix, matching pair_matrix shape.
    - `n_interfaces`: total interface count.

    Use the tables to surface "key residue pairs" (top of `pairs`) and "key
    residues" (top of `partner_*_residues`); use `pair_matrix` for a
    residue-level heatmap.
    """
    n = len(records)
    pmap = partner_map or {}
    rid = residue_identity if residue_identity is not None else _merge_residue_identity(records)

    pair_counts: Counter = Counter()
    p1_residue_counts: Counter = Counter()
    p2_residue_counts: Counter = Counter()
    for r in records:
        # Within a single interface, count distinct residue-pairs (collapse bond_type)
        # and distinct residues per side.
        seen_pairs: set[tuple] = set()
        seen_p1: set[tuple] = set()
        seen_p2: set[tuple] = set()
        for (k1, k2, _) in r.uniprot_pairs:
            seen_pairs.add((k1, k2))
            seen_p1.add(k1)
            seen_p2.add(k2)
        for p in seen_pairs:
            pair_counts[p] += 1
        for k in seen_p1:
            p1_residue_counts[k] += 1
        for k in seen_p2:
            p2_residue_counts[k] += 1

    def res_label(k):
        acc, pos, role = k
        partner = _label_partner(acc, role, pmap)
        aa = rid.get(k, "")
        return f"{partner}:{aa}{pos}"

    pairs_rows = []
    for (k1, k2), count in pair_counts.items():
        pairs_rows.append({
            "contact": f"{res_label(k1)} – {res_label(k2)}",
            "partner_1": res_label(k1),
            "partner_2": res_label(k2),
            "n_interfaces": int(count),
            "fraction": round(count / n, 3) if n else 0.0,
        })
    pairs_df = (
        pd.DataFrame(pairs_rows)
        .sort_values(["n_interfaces", "partner_1", "partner_2"], ascending=[False, True, True])
        .reset_index(drop=True)
    )

    def res_table(counts: Counter) -> pd.DataFrame:
        rows = []
        for k, c in counts.items():
            if c < min_residue_frequency:
                continue
            acc, pos, role = k
            rows.append({
                "residue": res_label(k),
                "accession": acc,
                "position": pos,
                "role": role,
                "n_interfaces": int(c),
                "fraction": round(c / n, 3) if n else 0.0,
            })
        return (
            pd.DataFrame(rows)
            .sort_values(["n_interfaces", "position"], ascending=[False, True])
            .reset_index(drop=True)
        )

    p1_table = res_table(p1_residue_counts)
    p2_table = res_table(p2_residue_counts)

    # Build pair_matrix in residue-frequency order
    p1_keys = [(row["accession"], row["position"], row["role"]) for _, row in p1_table.iterrows()]
    p2_keys = [(row["accession"], row["position"], row["role"]) for _, row in p2_table.iterrows()]
    p1_idx = {k: i for i, k in enumerate(p1_keys)}
    p2_idx = {k: i for i, k in enumerate(p2_keys)}
    matrix = np.zeros((len(p1_keys), len(p2_keys)), dtype=int)
    for (k1, k2), c in pair_counts.items():
        if k1 in p1_idx and k2 in p2_idx:
            matrix[p1_idx[k1], p2_idx[k2]] = c

    return {
        "pairs": pairs_df,
        "partner_1_residues": p1_table,
        "partner_2_residues": p2_table,
        "pair_matrix": matrix,
        "partner_1_labels": [res_label(k) for k in p1_keys],
        "partner_2_labels": [res_label(k) for k in p2_keys],
        "n_interfaces": n,
    }


def conserved_residues(
    records: list[InterfaceRecord], threshold: float = 1.0,
) -> set[tuple[str, int, int]]:
    """UniProt-keyed residues present in at least `threshold` * N interfaces."""
    if not records:
        return set()
    n = len(records)
    counts: Counter = Counter()
    for r in records:
        residues_in_interface: set = set()
        for k1, k2, _ in r.uniprot_pairs:
            residues_in_interface.add(k1)
            residues_in_interface.add(k2)
        counts.update(residues_in_interface)
    cutoff = threshold * n
    return {res for res, c in counts.items() if c >= cutoff}


def conserved_interaction_pairs(
    records: list[InterfaceRecord], threshold: float = 1.0,
) -> set[tuple]:
    if not records:
        return set()
    n = len(records)
    counts: Counter = Counter()
    for r in records:
        counts.update(r.uniprot_pairs)
    cutoff = threshold * n
    return {pair for pair, c in counts.items() if c >= cutoff}


def _merge_residue_identity(records: list[InterfaceRecord]) -> dict:
    """Union of per-record residue_identity maps; canonical UniProt one-letter codes."""
    merged: dict = {}
    for r in records:
        merged.update(r.residue_identity)
    return merged


# A cluster whose median resolution is this many Angstroms worse than the
# dominant cluster's gets a QC flag: hydrogen-bond detection is resolution
# dependent, so "lost" contacts may be undetected rather than absent.
DEFAULT_RESOLUTION_GAP_A = 1.0


def _resolutions_of(members: list[InterfaceRecord], metadata: dict) -> list[float]:
    out: list[float] = []
    for m in members:
        res = (metadata.get((m.pdb_id, m.assembly_id), {}) or {}).get("resolution")
        if res is not None:
            out.append(float(res))
    return out


def cluster_interpretation_report(
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    overlap: AnnotationOverlap,
    assembly_metadata: dict[tuple[str, str], dict] | None = None,
    partner_map: dict[tuple[str, int], str] | None = None,
    min_interfaces: int = 1,
    top_n_contacts: int = 10,
    conservation_threshold: float = 0.8,
    sparse_pair_threshold: int = DEFAULT_SPARSE_PAIR_THRESHOLD,
    tiny_area_threshold_a2: float = DEFAULT_TINY_AREA_THRESHOLD_A2,
    range_overlap_min: float = DEFAULT_RANGE_OVERLAP_MIN,
    resolution_gap_a: float = DEFAULT_RESOLUTION_GAP_A,
) -> pd.DataFrame:
    """One row per candidate interface interaction state (cluster).

    Rows are ordered by cluster size descending, ties broken by cluster id.

    Each cluster is treated as a candidate interface interaction state. For
    each cluster the report surfaces:

    Statistics (per-interface distributions inside the cluster):
      - residue-pair count (collapsed over bond_type)
      - typed-tuple interaction count (the input to Jaccard similarity)
      - PISA interface area
      - contact density = median typed tuples / median area

    Core contacts:
      typed contacts present in >= `conservation_threshold` of cluster
      members, capped at the top 10 by frequency.

    Cluster contacts (`cluster_contacts`):
      every typed contact present in the cluster, with the count inside the
      cluster and the count in the rest of the dataset shown side by side, at
      both interface and distinct-PDB-entry level. Ordered by in-cluster
      frequency and truncated to `top_n_contacts`; `cluster_contacts_shown` /
      `cluster_contacts_total` report how much was hidden.

      There is deliberately **no enrichment statistic here**. Clusters are
      defined by Jaccard similarity over these same typed contacts, so testing
      them against the clustering is circular by construction, since a cluster
      cannot fail to look enriched for the contacts that separated it. The
      counts are reported and the reader judges.

    Cluster annotations (`cluster_mutations` / `cluster_modifications` /
    `cluster_ligands`):
      the same treatment for each annotation stream. These are not circular,
      since annotations never enter the clustering, but the structures are not
      independent samples either (repeated assemblies, laboratory-correlated
      depositions), so these are also reported as counts, not as tests.

    QC warnings:
      - singleton cluster
      - sparse interaction fingerprint (median pair count < threshold)
      - small interface (median area < threshold)
      - residue range poorly overlaps the dominant (largest) cluster, which
        can flag polyprotein / processed-product / domain mixups within a
        single UniProt accession.
      - median resolution more than `resolution_gap_a` worse than the dominant
        cluster's, where missing contacts may be a detection artefact rather
        than a real difference in binding.

    Column names changed when the enrichment statistics were removed:
    `enriched_mutations` / `enriched_modifications` / `enriched_ligands` /
    `enriched_contacts` are now `cluster_*`, and
    `annotation_correlate_present` is now `has_annotations`.
    """
    cluster_of = {
        r.key: int(cluster_result.flat_assignment[i])
        for i, r in enumerate(records)
        if i < len(cluster_result.flat_assignment)
    }
    all_keys = {r.key for r in records}

    by_cluster: dict[int, list[InterfaceRecord]] = {}
    for r in records:
        by_cluster.setdefault(cluster_of.get(r.key, 0), []).append(r)

    metadata = assembly_metadata or {}
    pmap = partner_map or {}
    residue_identity = _merge_residue_identity(records)

    # Identify the dominant cluster (largest size; ties broken by id) so
    # other clusters can be compared against its UniProt residue range.
    dominant_cid = (
        max(by_cluster.items(), key=lambda kv: (len(kv[1]), -kv[0]))[0]
        if by_cluster else None
    )
    dominant_range = (
        _cluster_residue_ranges(by_cluster[dominant_cid])
        if dominant_cid is not None else None
    )
    dominant_resolutions = (
        _resolutions_of(by_cluster[dominant_cid], metadata)
        if dominant_cid is not None else None
    )

    rows = []
    # Largest state first: cluster ids come from `fcluster` and carry no
    # meaning in their numeric order, whereas size orders the report by how
    # much evidence each state rests on, and pushes singletons to the end.
    # Ties broken by cluster id so repeated runs order identically.
    for cid, members in sorted(by_cluster.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        member_keys = {r.key for r in members}
        non_member_keys = all_keys - member_keys
        non_member_records = [r for r in records if r.key not in member_keys]
        member_pdb_ids = sorted({r.pdb_id for r in members})

        methods = Counter()
        resolutions: list[float] = []
        for m in members:
            meta = metadata.get((m.pdb_id, m.assembly_id), {})
            method = meta.get("experimental_method") or "unknown"
            methods[method] += 1
            res = meta.get("resolution")
            if res is not None:
                resolutions.append(float(res))
        methods_str = ", ".join(f"{name} ({n})" for name, n in methods.most_common())
        if resolutions:
            res_min = min(resolutions)
            res_max = max(resolutions)
            res_median = sorted(resolutions)[len(resolutions) // 2]
            resolution_str = (
                f"{res_min:.2f}–{res_max:.2f} Å (median {res_median:.2f}, n={len(resolutions)})"
            )
        else:
            resolution_str = ""

        pair_counts, tuple_counts = _per_interface_pair_counts(members)
        areas = _per_interface_areas(members)
        pair_stats = _stats(pair_counts)
        tuple_stats = _stats(tuple_counts)
        area_stats = _stats(areas)
        contact_density_median = (
            (tuple_stats["median"] / area_stats["median"])
            if (tuple_stats["median"] and area_stats["median"])
            else None
        )

        if areas:
            area_str = (
                f"{area_stats['min']:.0f}–{area_stats['max']:.0f} Å² "
                f"(median {area_stats['median']:.0f}, n={len(areas)})"
            )
        else:
            area_str = ""

        core = _core_contacts(members, conservation_threshold, top_n=10)

        member_entries = {r.pdb_id for r in members}
        non_member_entries = {r.pdb_id for r in non_member_records}

        contact_profile, n_contacts_total = _contact_profile(
            members, non_member_records, min_interfaces, top_n_contacts,
        )

        mut_profile = _label_profile(
            overlap.mutations, member_keys, non_member_keys,
            "mutation_label", member_entries, non_member_entries,
        )
        mod_profile = _label_profile(
            overlap.modifications, member_keys, non_member_keys,
            "modification_chem_comp_id", member_entries, non_member_entries,
        )
        lig_profile = _label_profile(
            overlap.ligands, member_keys, non_member_keys,
            "ligand_chem_comp_id", member_entries, non_member_entries,
        )

        n_with_mut = len(_member_interfaces_with(overlap.mutations, member_keys))
        n_with_mod = len(_member_interfaces_with(overlap.modifications, member_keys))
        n_with_lig = len(_member_interfaces_with(overlap.ligands, member_keys))

        has_annotations = bool(mut_profile or mod_profile or lig_profile)

        # QC. The dominant cluster is always compared against itself, so the
        # range-overlap rule is skipped for it (overlap = 1.0 by construction).
        range_self = _cluster_residue_ranges(members)
        warns = _qc_warnings(
            cid, members, pair_counts, areas,
            range_dominant=(None if cid == dominant_cid else dominant_range),
            range_self=range_self,
            sparse_pair_threshold=sparse_pair_threshold,
            tiny_area_threshold_a2=tiny_area_threshold_a2,
            range_overlap_min=range_overlap_min,
            resolutions=resolutions,
            dominant_resolutions=(None if cid == dominant_cid else dominant_resolutions),
            resolution_gap_a=resolution_gap_a,
        )

        notes = _build_notes(
            cid, len(members), len(non_member_keys),
            mut_profile, mod_profile, lig_profile,
            n_with_mut, n_with_mod, n_with_lig,
            methods_str, resolution_str,
            pair_stats=pair_stats,
            tuple_stats=tuple_stats,
            area_stats=area_stats,
            contact_density_median=contact_density_median,
            core_contacts=core,
            warnings=warns,
            partner_map=pmap,
            residue_identity=residue_identity,
            contact_profile=contact_profile,
            n_contacts_total=n_contacts_total,
        )

        rows.append({
            "cluster_id": cid,
            "cluster_size": len(members),
            "member_pdb_ids": ",".join(member_pdb_ids),
            "experimental_methods": methods_str,
            "resolution_range": resolution_str,
            "interface_area_range": area_str,
            "residue_pair_count_min": pair_stats["min"],
            "residue_pair_count_median": pair_stats["median"],
            "residue_pair_count_mean": pair_stats["mean"],
            "residue_pair_count_max": pair_stats["max"],
            "interaction_count_min": tuple_stats["min"],
            "interaction_count_median": tuple_stats["median"],
            "interaction_count_mean": tuple_stats["mean"],
            "interaction_count_max": tuple_stats["max"],
            "interface_area_min": area_stats["min"],
            "interface_area_median": area_stats["median"],
            "interface_area_mean": area_stats["mean"],
            "interface_area_max": area_stats["max"],
            "contact_density_median": contact_density_median,
            "core_contacts": _format_core_contacts(core, pmap, residue_identity),
            "cluster_mutations": _format_profile(mut_profile),
            "cluster_modifications": _format_profile(mod_profile),
            "cluster_ligands": _format_profile(lig_profile),
            "cluster_contacts": _format_contacts(contact_profile, pmap, residue_identity),
            "cluster_contacts_shown": len(contact_profile),
            "cluster_contacts_total": n_contacts_total,
            "interfaces_with_any_mutation": n_with_mut,
            "interfaces_with_any_modification": n_with_mod,
            "interfaces_with_any_ligand": n_with_lig,
            "has_annotations": has_annotations,
            "qc_warnings": "; ".join(warns),
            "notes": notes,
        })
    return pd.DataFrame(rows)


def _stats(values: list[float]) -> dict:
    """min/median/mean/max for a list of numbers, or all-None if empty."""
    if not values:
        return {"min": None, "median": None, "mean": None, "max": None}
    return {
        "min": float(min(values)),
        "median": float(statistics.median(values)),
        "mean": float(statistics.fmean(values)),
        "max": float(max(values)),
    }


def _per_interface_pair_counts(members: list[InterfaceRecord]) -> tuple[list[int], list[int]]:
    """Per-interface (residue-pair count ignoring bond_type, typed-tuple count).

    A "residue pair" collapses bond_type, so a single residue–residue contact
    with both an h-bond and a salt bridge counts once. A "typed tuple" keeps
    bond_type, matching the inputs to the Jaccard similarity.
    """
    pair_counts: list[int] = []
    tuple_counts: list[int] = []
    for r in members:
        pair_counts.append(len({(k1, k2) for (k1, k2, _) in r.uniprot_pairs}))
        tuple_counts.append(len(r.uniprot_pairs))
    return pair_counts, tuple_counts


def _per_interface_areas(members: list[InterfaceRecord]) -> list[float]:
    out: list[float] = []
    for r in members:
        a = (r.interface_info or {}).get("interface_area")
        if a is not None:
            try:
                out.append(float(a))
            except (TypeError, ValueError):
                continue
    return out


def _core_contacts(
    members: list[InterfaceRecord],
    threshold: float,
    top_n: int = 10,
) -> list[dict]:
    """Typed contacts present in >= threshold fraction of cluster members.

    Returns a list of dicts (sorted by frequency desc, then by contact label
    asc), each with: contact (the typed tuple), count, denominator, fraction.
    """
    n = len(members)
    if n == 0:
        return []
    counts: Counter = Counter()
    for r in members:
        counts.update(r.uniprot_pairs)
    cutoff = threshold * n
    rows = [
        {
            "contact": c,
            "count": int(k),
            "denominator": n,
            "fraction": float(k / n),
        }
        for c, k in counts.items()
        if k >= cutoff
    ]
    rows.sort(key=lambda r: (-r["count"], _contact_sort_key(r["contact"])))
    return rows[:top_n]


def _contact_sort_key(contact: tuple) -> tuple:
    """Stable, label-like sort key for a typed contact tuple."""
    (acc1, pos1, role1), (acc2, pos2, role2), bond = contact
    return (acc1, role1, pos1, acc2, role2, pos2, bond)


def _format_core_contacts(items: list[dict], partner_map: dict, residue_identity: dict) -> str:
    parts: list[str] = []
    for it in items:
        (acc1, pos1, role1), (acc2, pos2, role2), bond = it["contact"]
        l1 = _label_partner(acc1, role1, partner_map)
        l2 = _label_partner(acc2, role2, partner_map)
        aa1 = residue_identity.get((acc1, pos1, role1), "")
        aa2 = residue_identity.get((acc2, pos2, role2), "")
        pct = 100.0 * it["fraction"]
        parts.append(
            f"{l1}:{aa1}{pos1}-{l2}:{aa2}{pos2} {bond} "
            f"({it['count']}/{it['denominator']}, {pct:.1f}%)"
        )
    return "; ".join(parts)


def _cluster_residue_ranges(members: list[InterfaceRecord]) -> dict[tuple[str, int], tuple[int, int]]:
    """Per (accession, role) -> (min UniProt position, max UniProt position)."""
    by_partner: dict[tuple[str, int], list[int]] = {}
    for r in members:
        for (k1, k2, _) in r.uniprot_pairs:
            for k in (k1, k2):
                acc, pos, role = k
                by_partner.setdefault((acc, role), []).append(int(pos))
    return {key: (min(vals), max(vals)) for key, vals in by_partner.items() if vals}


def _interval_overlap_fraction(
    a: tuple[int, int],
    b: tuple[int, int],
) -> float:
    """Overlap length / span of the larger interval. Returns 0 if no overlap.

    Span uses the larger of the two ranges so that a small range entirely
    contained in a large range does not score 1.0, because we want to detect when a
    cluster's residues sit in a different region of the protein, not when
    they're a focused subset.
    """
    lo = max(a[0], b[0])
    hi = min(a[1], b[1])
    overlap = max(0, hi - lo)
    span = max(a[1] - a[0], b[1] - b[0])
    if span <= 0:
        # Both ranges are single residues: return 1.0 if equal, else 0.
        return 1.0 if a == b else 0.0
    return overlap / span


def _qc_warnings(
    cid: int,
    members: list[InterfaceRecord],
    pair_counts: list[int],
    areas: list[float],
    range_dominant: dict[tuple[str, int], tuple[int, int]] | None,
    range_self: dict[tuple[str, int], tuple[int, int]],
    sparse_pair_threshold: int,
    tiny_area_threshold_a2: float,
    range_overlap_min: float,
    resolutions: list[float] | None = None,
    dominant_resolutions: list[float] | None = None,
    resolution_gap_a: float = DEFAULT_RESOLUTION_GAP_A,
) -> list[str]:
    """Cluster-level QC flags per the spec."""
    warns: list[str] = []
    if len(members) == 1:
        warns.append("singleton cluster; interpret cautiously")
    if pair_counts:
        median_pairs = statistics.median(pair_counts)
        if median_pairs < sparse_pair_threshold:
            warns.append(
                "sparse interaction fingerprint; clustering may reflect few contacts"
            )
    if areas:
        median_area = statistics.median(areas)
        if median_area < tiny_area_threshold_a2:
            warns.append(
                "small interface area; possible weak or non-biological interface"
            )
    if range_dominant is not None and range_self is not None and range_self:
        # Compare each shared (accession, role) range against the dominant
        # cluster's range; flag if any partner falls below the overlap floor.
        for partner_key, self_range in range_self.items():
            dom_range = range_dominant.get(partner_key)
            if dom_range is None:
                continue
            if _interval_overlap_fraction(self_range, dom_range) < range_overlap_min:
                warns.append(
                    "residue range poorly overlaps the dominant cluster; this may "
                    "indicate a different domain, processed product, or polyprotein "
                    "segment"
                )
                break
    # Resolution confounding: hydrogen bonds are geometry-derived, so a cluster
    # resolved substantially worse than the dominant one may differ because
    # contacts went undetected, not because the interface changed.
    if resolutions and dominant_resolutions:
        median_res = statistics.median(resolutions)
        dominant_median = statistics.median(dominant_resolutions)
        if median_res - dominant_median > resolution_gap_a:
            warns.append(
                f"median resolution {median_res:.2f} A vs {dominant_median:.2f} A for the "
                "dominant cluster; missing contacts here may be undetected rather than absent"
            )
    return warns


def _contact_profile(
    member_records: list[InterfaceRecord],
    non_member_records: list[InterfaceRecord],
    min_interfaces: int,
    top_n: int,
) -> tuple[list[dict], int]:
    """Describe the contacts of one cluster. No test, no threshold on effect.

    For every typed contact tuple present in the cluster, report how many
    cluster members carry it and how many non-members do, at two levels of
    aggregation:

    - *interfaces*: the unit the clustering operates on;
    - *entries*: distinct PDB entries, which de-duplicates the case of one
      deposition contributing several assemblies of the same coordinates.

    Deliberately no fold-enrichment, p-value, or significance filter. Clusters
    are defined by Jaccard similarity over these same tuples, so any statistic
    testing them against the clustering is circular by construction: a cluster
    cannot fail to appear enriched for the contacts that separated it. The
    counts are reported side by side and the reader judges.

    Returns (rows, n_total) where `n_total` is the number of contacts before
    `top_n` truncation, so callers can say how much was hidden.
    """
    in_size = len(member_records)
    out_size = len(non_member_records)
    if in_size == 0:
        return [], 0

    in_entries_total = len({r.pdb_id for r in member_records})
    out_entries_total = len({r.pdb_id for r in non_member_records})

    in_counts: Counter = Counter()
    out_counts: Counter = Counter()
    in_entry_sets: dict = {}
    out_entry_sets: dict = {}
    for r in member_records:
        for c in r.uniprot_pairs:
            in_counts[c] += 1
            in_entry_sets.setdefault(c, set()).add(r.pdb_id)
    for r in non_member_records:
        for c in r.uniprot_pairs:
            out_counts[c] += 1
            out_entry_sets.setdefault(c, set()).add(r.pdb_id)

    rows = []
    for contact, in_n in in_counts.items():
        if in_n < min_interfaces:
            continue
        out_n = int(out_counts.get(contact, 0))
        rows.append({
            "contact": contact,
            "in_cluster_interfaces": int(in_n),
            "in_cluster_size": in_size,
            "in_cluster_fraction": in_n / in_size,
            "out_cluster_interfaces": out_n,
            "out_cluster_size": out_size,
            "out_cluster_fraction": (out_n / out_size) if out_size else None,
            "in_cluster_entries": len(in_entry_sets.get(contact, ())),
            "in_cluster_entry_total": in_entries_total,
            "out_cluster_entries": len(out_entry_sets.get(contact, ())),
            "out_cluster_entry_total": out_entries_total,
        })

    # Most frequent inside the cluster first; among equals, the ones rarest
    # outside it first, since those are what a reader will look at. This is a
    # display order, not a ranking by significance.
    rows.sort(key=lambda r: (
        -r["in_cluster_fraction"],
        r["out_cluster_fraction"] if r["out_cluster_fraction"] is not None else 0.0,
        str(r["contact"]),
    ))
    return rows[:top_n], len(rows)


def compare_clusters(
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    cluster_a: int,
    cluster_b: int,
    partner_map: dict[tuple[str, int], str] | None = None,
    min_count: int = 1,
    top_n: int | None = 20,
    residue_identity: dict | None = None,
) -> pd.DataFrame:
    """Head-to-head contact comparison between two clusters.

    For each UniProt-keyed interaction-pair tuple present in either cluster,
    report distinct-interface counts in both, fractions, and a directional
    label (A-only / B-only / biased_to_A / biased_to_B / tied). Ranked by
    `magnitude` = |fraction_in_a - fraction_in_b| descending, so the most
    divergent contacts surface first.

    Use after running cluster_interpretation_report when you want to ask
    "what specifically distinguishes cluster X from cluster Y?".
    """
    cluster_of = {
        r.key: int(cluster_result.flat_assignment[i])
        for i, r in enumerate(records)
        if i < len(cluster_result.flat_assignment)
    }
    members_a = [r for r in records if cluster_of.get(r.key) == cluster_a]
    members_b = [r for r in records if cluster_of.get(r.key) == cluster_b]
    if not members_a:
        raise ValueError(f"No members found in cluster {cluster_a}")
    if not members_b:
        raise ValueError(f"No members found in cluster {cluster_b}")

    counts_a: Counter = Counter()
    counts_b: Counter = Counter()
    for r in members_a:
        for c in r.uniprot_pairs:
            counts_a[c] += 1
    for r in members_b:
        for c in r.uniprot_pairs:
            counts_b[c] += 1

    size_a, size_b = len(members_a), len(members_b)
    pmap = partner_map or {}
    rid = residue_identity if residue_identity is not None else _merge_residue_identity(records)

    rows = []
    for c in set(counts_a) | set(counts_b):
        a_n = int(counts_a.get(c, 0))
        b_n = int(counts_b.get(c, 0))
        if a_n < min_count and b_n < min_count:
            continue
        frac_a = a_n / size_a
        frac_b = b_n / size_b
        if a_n > 0 and b_n == 0:
            direction = "A-only"
        elif b_n > 0 and a_n == 0:
            direction = "B-only"
        elif frac_a > frac_b:
            direction = "biased_to_A"
        elif frac_b > frac_a:
            direction = "biased_to_B"
        else:
            direction = "tied"

        (acc1, pos1, role1), (acc2, pos2, role2), bond = c
        l1 = _label_partner(acc1, role1, pmap)
        l2 = _label_partner(acc2, role2, pmap)
        aa1 = rid.get((acc1, pos1, role1), "")
        aa2 = rid.get((acc2, pos2, role2), "")
        rows.append({
            "contact": f"{l1}:{aa1}{pos1}-{l2}:{aa2}{pos2} {bond}",
            f"in_cluster_{cluster_a}": f"{a_n}/{size_a}",
            f"in_cluster_{cluster_b}": f"{b_n}/{size_b}",
            "fraction_a": round(frac_a, 3),
            "fraction_b": round(frac_b, 3),
            "direction": direction,
            "magnitude": round(abs(frac_a - frac_b), 3),
        })

    df = pd.DataFrame(rows).sort_values("magnitude", ascending=False).reset_index(drop=True)
    if top_n is not None:
        df = df.head(top_n).reset_index(drop=True)
    return df


def _label_partner(unp_acc: str, role: int, partner_map: dict) -> str:
    label = partner_map.get((unp_acc, role), unp_acc)
    return label.replace(" (role 1)", "").replace(" (role 2)", "")


def _format_contacts(items: list[dict], partner_map: dict, residue_identity: dict | None = None) -> str:
    rid = residue_identity or {}
    parts = []
    for it in items:
        (acc1, pos1, role1), (acc2, pos2, role2), bond = it["contact"]
        l1 = _label_partner(acc1, role1, partner_map)
        l2 = _label_partner(acc2, role2, partner_map)
        aa1 = rid.get((acc1, pos1, role1), "")
        aa2 = rid.get((acc2, pos2, role2), "")
        contact_str = f"{l1}:{aa1}{pos1}-{l2}:{aa2}{pos2} {bond}"
        suffix = (
            f"{it['in_cluster_interfaces']}/{it['in_cluster_size']} interfaces "
            f"= {it['in_cluster_fraction']:.0%}, "
            f"{it['in_cluster_entries']}/{it['in_cluster_entry_total']} entries"
        )
        if it.get("out_cluster_size"):
            suffix += (
                f"; rest {it['out_cluster_interfaces']}/{it['out_cluster_size']} "
                f"= {it['out_cluster_fraction']:.0%}"
            )
        parts.append(f"{contact_str} ({suffix})")
    return "; ".join(parts)


def _label_profile(
    df: pd.DataFrame,
    member_keys: set,
    non_member_keys: set,
    label_col: str,
    member_entries: set,
    non_member_entries: set,
) -> list[dict]:
    """Describe one annotation stream inside a cluster vs the rest. No test.

    Counts are by *distinct interface* (de-duplicated on
    pdb_id/assembly_id/interface_id) and by *distinct PDB entry*, so a ligand
    with several contact tuples on one interface counts once, and one entry
    contributing several assemblies counts once at entry level.

    Unlike the contact profile this comparison is not circular, since
    annotations never enter the clustering, but the structures are still not
    independent samples (repeated assemblies, laboratory-correlated
    depositions), so no significance is computed here either.
    """
    if df.empty or label_col not in df.columns:
        return []
    in_size = len(member_keys)
    out_size = len(non_member_keys)
    if in_size == 0:
        return []

    keys = list(zip(df["pdb_id"], df["assembly_id"], df["interface_id"]))
    df_keyed = df.assign(_key=keys)

    in_df = df_keyed[df_keyed["_key"].isin(member_keys)]
    out_df = df_keyed[df_keyed["_key"].isin(non_member_keys)]

    in_counts = in_df.drop_duplicates(["_key", label_col])[label_col].value_counts()
    out_counts = out_df.drop_duplicates(["_key", label_col])[label_col].value_counts()
    # Distinct PDB entries carrying each label, on each side.
    in_entry_counts = in_df.drop_duplicates(["pdb_id", label_col])[label_col].value_counts()
    out_entry_counts = out_df.drop_duplicates(["pdb_id", label_col])[label_col].value_counts()

    rows = []
    for label, in_n in in_counts.items():
        in_n = int(in_n)
        out_n = int(out_counts.get(label, 0))
        rows.append({
            "label": label,
            "in_cluster_interfaces": in_n,
            "in_cluster_size": in_size,
            "in_cluster_fraction": in_n / in_size,
            "out_cluster_interfaces": out_n,
            "out_cluster_size": out_size,
            "out_cluster_fraction": (out_n / out_size) if out_size else None,
            "in_cluster_entries": int(in_entry_counts.get(label, 0)),
            "in_cluster_entry_total": len(member_entries),
            "out_cluster_entries": int(out_entry_counts.get(label, 0)),
            "out_cluster_entry_total": len(non_member_entries),
        })
    rows.sort(key=lambda r: (
        -r["in_cluster_fraction"],
        r["out_cluster_fraction"] if r["out_cluster_fraction"] is not None else 0.0,
        str(r["label"]),
    ))
    return rows


def _member_interfaces_with(df: pd.DataFrame, member_keys: set) -> set:
    if df.empty:
        return set()
    keys = set(zip(df["pdb_id"], df["assembly_id"], df["interface_id"]))
    return keys & member_keys


def _format_label(item: dict) -> str:
    """`ACE2 K353A (12/44 interfaces, 10/38 entries; rest of dataset 3/86)`."""
    base = (
        f"{item['label']} ({item['in_cluster_interfaces']}/{item['in_cluster_size']} interfaces "
        f"= {item['in_cluster_fraction']:.0%}, "
        f"{item['in_cluster_entries']}/{item['in_cluster_entry_total']} entries"
    )
    if item.get("out_cluster_size"):
        base += (
            f"; rest {item['out_cluster_interfaces']}/{item['out_cluster_size']} "
            f"= {item['out_cluster_fraction']:.0%}"
        )
    return base + ")"


def _format_profile(items: list[dict]) -> str:
    return "; ".join(_format_label(it) for it in items)


def _build_notes(
    cid: int, cluster_size: int, n_non_members: int,
    muts: list[dict], mods: list[dict], ligs: list[dict],
    n_with_mut: int, n_with_mod: int, n_with_lig: int,
    methods_str: str = "", resolution_str: str = "",
    pair_stats: dict | None = None,
    tuple_stats: dict | None = None,
    area_stats: dict | None = None,
    contact_density_median: float | None = None,
    core_contacts: list[dict] | None = None,
    warnings: list[str] | None = None,
    partner_map: dict | None = None,
    residue_identity: dict | None = None,
    contact_profile: list[dict] | None = None,
    n_contacts_total: int = 0,
) -> str:
    method_suffix = ""
    if methods_str or resolution_str:
        bits = []
        if methods_str:
            bits.append(f"methods: {methods_str}")
        if resolution_str:
            bits.append(f"resolution: {resolution_str}")
        method_suffix = " [" + "; ".join(bits) + "]"

    state_label = (
        f"Interface interaction state {cid} ({cluster_size} interfaces)"
    )

    if n_non_members == 0:
        head = (
            f"{state_label}: single-state dataset; "
            f"no rest of dataset to compare against."
            f"{method_suffix}"
        )
        return _attach_state_diagnostics(
            head,
            pair_stats=pair_stats, tuple_stats=tuple_stats,
            area_stats=area_stats, contact_density_median=contact_density_median,
            core_contacts=core_contacts, warnings=warnings,
            partner_map=partner_map, residue_identity=residue_identity,
            contact_profile=contact_profile, n_contacts_total=n_contacts_total,
        )
    if cluster_size == 1:
        head = (
            f"{state_label}: singleton; no within-state pattern to summarise."
            f"{method_suffix}"
        )
        return _attach_state_diagnostics(
            head,
            pair_stats=pair_stats, tuple_stats=tuple_stats,
            area_stats=area_stats, contact_density_median=contact_density_median,
            core_contacts=core_contacts, warnings=warnings,
            partner_map=partner_map, residue_identity=residue_identity,
            contact_profile=contact_profile, n_contacts_total=n_contacts_total,
        )

    parts = []
    if muts:
        parts.append("mutations at the interface: " +
                     ", ".join(_format_label(m) for m in muts))
    if mods:
        parts.append("modifications at the interface: " +
                     ", ".join(_format_label(m) for m in mods))
    if ligs:
        parts.append("ligands at the interface: " +
                     ", ".join(_format_label(m) for m in ligs))

    density = []
    if n_with_mut:
        density.append(f"{n_with_mut}/{cluster_size} carry an at-interface mutation")
    if n_with_mod:
        density.append(f"{n_with_mod}/{cluster_size} carry an at-interface modification")
    if n_with_lig:
        density.append(f"{n_with_lig}/{cluster_size} carry an at-interface ligand contact")

    if parts:
        body = "; ".join(parts)
        if density:
            body += " (" + ", ".join(density) + ")"
        head = f"{state_label}: {body}.{method_suffix}"
    elif density:
        head = (
            f"{state_label}: no mutation, modification, or ligand recorded at the "
            f"interface. Annotation density: " + ", ".join(density) + "."
            f"{method_suffix}"
        )
    else:
        head = (
            f"{state_label}: no mutation, modification, or ligand at the interface, "
            "and no member carries any. Members may share crystal form, era, or "
            "refinement protocol."
            f"{method_suffix}"
        )
    return _attach_state_diagnostics(
        head,
        pair_stats=pair_stats, tuple_stats=tuple_stats,
        area_stats=area_stats, contact_density_median=contact_density_median,
        core_contacts=core_contacts, warnings=warnings,
        partner_map=partner_map, residue_identity=residue_identity,
        contact_profile=contact_profile, n_contacts_total=n_contacts_total,
    )


def _attach_state_diagnostics(
    head: str,
    pair_stats: dict | None,
    tuple_stats: dict | None,
    area_stats: dict | None,
    contact_density_median: float | None,
    core_contacts: list[dict] | None,
    warnings: list[str] | None,
    partner_map: dict | None,
    residue_identity: dict | None,
    contact_profile: list[dict] | None = None,
    n_contacts_total: int = 0,
) -> str:
    """Append residue-pair / interaction / area / contacts / warnings lines."""
    lines = [head]
    stat_lines = _stat_lines(pair_stats, tuple_stats, area_stats, contact_density_median)
    lines.extend(stat_lines)
    rid = residue_identity or {}
    pmap = partner_map or {}
    if core_contacts:
        lines.append(
            f"  Core contacts (>= conservation_threshold of members): "
            + _format_core_contacts(core_contacts, pmap, rid)
        )
    if contact_profile:
        shown = len(contact_profile)
        header = (
            f"  Contacts in this cluster (counts only, no significance test; "
            f"showing {shown} of {n_contacts_total})"
        )
        lines.append(f"{header}: " + _format_contacts(contact_profile, pmap, rid))
    if warnings:
        lines.append("  QC warnings: " + "; ".join(warnings))
    return "\n".join(lines)


def _stat_lines(
    pair_stats: dict | None,
    tuple_stats: dict | None,
    area_stats: dict | None,
    contact_density_median: float | None,
) -> list[str]:
    out: list[str] = []
    if pair_stats and pair_stats.get("median") is not None:
        out.append(
            f"  Residue pair count: {int(pair_stats['min'])}–{int(pair_stats['max'])}, "
            f"median {int(pair_stats['median'])}"
        )
    if tuple_stats and tuple_stats.get("median") is not None:
        out.append(
            f"  Interaction count: {int(tuple_stats['min'])}–{int(tuple_stats['max'])}, "
            f"median {int(tuple_stats['median'])}"
        )
    if area_stats and area_stats.get("median") is not None:
        out.append(
            f"  Interface area: {area_stats['min']:.0f}–{area_stats['max']:.0f} Å², "
            f"median {area_stats['median']:.0f}"
        )
    if contact_density_median is not None:
        out.append(
            f"  Contact density: median {contact_density_median:.4f} interactions/Å²"
        )
    return out


def compare_cluster_contacts(
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    cluster_a: int | None = None,
    cluster_b: int | None = None,
    typed: bool = True,
    top_n: int = 20,
    partner_map: dict[tuple[str, int], str] | None = None,
    residue_identity: dict | None = None,
    a_enriched_threshold: float = 0.5,
    shared_core_threshold: float = 0.8,
) -> pd.DataFrame:
    """Interface rewiring table between two interface interaction states.

    Compares the contact frequency profiles of two clusters and labels each
    contact as `shared core`, `A-enriched`, `B-enriched`, or `rare`.

    Defaults to the two largest non-singleton clusters when ids are omitted.

    Args:
        cluster_a, cluster_b: cluster ids to compare. If None, the two
            largest non-singleton clusters are selected. If only one
            non-singleton cluster exists, raises ValueError.
        typed: when True (default) compares typed tuples
            `(residue_1, residue_2, bond_type)`; when False collapses
            bond_type so a residue-pair shared across bond types counts once.
        top_n: rows to keep, sorted by absolute fraction_difference desc.
        a_enriched_threshold: |fraction_diff| above which a contact is
            tagged "higher in A" / "higher in B" (default 0.5). A descriptive
            convention, not a test: with unequal cluster sizes a single
            structure can move the fraction in the smaller state.
        shared_core_threshold: minimum fraction in BOTH clusters for the
            "shared core" label (default 0.8).

    Returns DataFrame with columns:
        contact, cluster_A_count, cluster_B_count,
        cluster_A_fraction, cluster_B_fraction,
        fraction_difference, contact_direction.
    """
    cluster_of = {
        r.key: int(cluster_result.flat_assignment[i])
        for i, r in enumerate(records)
        if i < len(cluster_result.flat_assignment)
    }
    by_cluster: dict[int, list[InterfaceRecord]] = {}
    for r in records:
        by_cluster.setdefault(cluster_of.get(r.key, 0), []).append(r)

    if cluster_a is None or cluster_b is None:
        non_singletons = sorted(
            ((cid, len(m)) for cid, m in by_cluster.items() if len(m) >= 2),
            key=lambda kv: (-kv[1], kv[0]),
        )
        if len(non_singletons) < 2:
            raise ValueError(
                "compare_cluster_contacts: fewer than 2 non-singleton clusters; "
                "pass cluster_a and cluster_b explicitly."
            )
        if cluster_a is None:
            cluster_a = non_singletons[0][0]
        if cluster_b is None:
            cluster_b = (
                non_singletons[1][0]
                if non_singletons[1][0] != cluster_a
                else non_singletons[0][0]
            )

    members_a = by_cluster.get(int(cluster_a), [])
    members_b = by_cluster.get(int(cluster_b), [])
    if not members_a:
        raise ValueError(f"No members found in cluster {cluster_a}")
    if not members_b:
        raise ValueError(f"No members found in cluster {cluster_b}")

    def _contacts_of(r: InterfaceRecord) -> set:
        if typed:
            return r.uniprot_pairs
        return {(k1, k2) for (k1, k2, _) in r.uniprot_pairs}

    counts_a: Counter = Counter()
    counts_b: Counter = Counter()
    for r in members_a:
        counts_a.update(_contacts_of(r))
    for r in members_b:
        counts_b.update(_contacts_of(r))

    size_a, size_b = len(members_a), len(members_b)
    pmap = partner_map or {}
    rid = residue_identity if residue_identity is not None else _merge_residue_identity(records)

    rows = []
    for c in set(counts_a) | set(counts_b):
        a_n = int(counts_a.get(c, 0))
        b_n = int(counts_b.get(c, 0))
        frac_a = a_n / size_a
        frac_b = b_n / size_b
        diff = frac_a - frac_b

        if frac_a >= shared_core_threshold and frac_b >= shared_core_threshold:
            direction = "shared core"
        elif diff >= a_enriched_threshold:
            direction = "higher in A"
        elif diff <= -a_enriched_threshold:
            direction = "higher in B"
        else:
            direction = "rare"

        rows.append({
            "contact": _format_contact_label(c, typed, pmap, rid),
            "cluster_A_count": a_n,
            "cluster_B_count": b_n,
            "cluster_A_fraction": round(frac_a, 3),
            "cluster_B_fraction": round(frac_b, 3),
            "fraction_difference": round(diff, 3),
            "contact_direction": direction,
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df.assign(_abs=0.0).drop(columns=["_abs"])
    df = (
        df.assign(_abs=df["fraction_difference"].abs())
        .sort_values(["_abs", "contact"], ascending=[False, True])
        .drop(columns=["_abs"])
        .reset_index(drop=True)
    )
    df.attrs["cluster_a"] = int(cluster_a)
    df.attrs["cluster_b"] = int(cluster_b)
    df.attrs["size_a"] = size_a
    df.attrs["size_b"] = size_b
    df.attrs["typed"] = bool(typed)
    if top_n is not None:
        df = df.head(top_n).reset_index(drop=True)
    return df


def _format_contact_label(contact, typed: bool, partner_map: dict, residue_identity: dict) -> str:
    if typed:
        (acc1, pos1, role1), (acc2, pos2, role2), bond = contact
    else:
        (acc1, pos1, role1), (acc2, pos2, role2) = contact
        bond = ""
    l1 = _label_partner(acc1, role1, partner_map)
    l2 = _label_partner(acc2, role2, partner_map)
    aa1 = residue_identity.get((acc1, pos1, role1), "")
    aa2 = residue_identity.get((acc2, pos2, role2), "")
    base = f"{l1}:{aa1}{pos1}-{l2}:{aa2}{pos2}"
    return f"{base} {bond}".strip()


# Conservation level cutoffs for residue_frequencies in the JSON export.
# Inclusive lower bounds; "rare" is anything < 0.20.
CONSERVATION_LEVEL_STRONG = 0.80
CONSERVATION_LEVEL_MEDIUM = 0.50
CONSERVATION_LEVEL_WEAK = 0.20


def export_interface_frequency_json(
    records: list[InterfaceRecord],
    partner_map: dict[tuple[str, int], str] | None = None,
    complex_id: str | None = None,
    complex_name: str | None = None,
    oligomeric_state: str | None = None,
    config: object | None = None,
    typed_contacts: bool = True,
    complex_details: dict | None = None,
) -> dict:
    """Export per-residue and per-contact frequencies for frontend visualisation.

    Builds a self-contained JSON describing UniProt residue-level interface
    conservation across all retained interfaces. Intended for downstream use
    cases such as Mol* colouring by conservation, or contact-frequency
    matrices/heatmaps in a frontend.

    The file is written to `config.output_dir` (created if missing) or to the
    current working directory if `config` is None or `config.output_dir` is
    None. The dictionary is also returned so callers can build preview tables
    or pass it onward without re-reading from disk.

    Args:
        records: built InterfaceRecord list (one per retained interface).
        partner_map: as produced by `build_partner_map`; used to derive the
            display label when richer metadata is unavailable.
        complex_id, complex_name, oligomeric_state: surfaced verbatim in the
            metadata block.
        config: the notebook `Config` instance; only `output_dir` is read.
        typed_contacts: if True, contact keys include bond_type so a
            residue-pair with multiple bond types contributes one row per
            bond type. If False, bond_type is collapsed.
        complex_details: optional; if provided, gene_name and name are
            extracted from `participants` for richer partner metadata.

    Returns the export dictionary. Side effect: writes
    `{complex_id}.json` to the resolved output directory.
    """
    n_interfaces = len(records)
    pmap = partner_map or {}
    residue_identity = _merge_residue_identity(records)

    partner_meta = _build_partner_metadata(records, pmap, complex_details)

    residue_counts: Counter = Counter()
    contact_counts: Counter = Counter()
    for r in records:
        seen_residues: set = set()
        seen_contacts: set = set()
        for (k1, k2, bond) in r.uniprot_pairs:
            seen_residues.add(k1)
            seen_residues.add(k2)
            if typed_contacts:
                seen_contacts.add((k1, k2, bond))
            else:
                seen_contacts.add((k1, k2))
        for k in seen_residues:
            residue_counts[k] += 1
        for c in seen_contacts:
            contact_counts[c] += 1

    residue_frequencies = []
    for key, count in residue_counts.items():
        acc, pos, role = key
        frac = (count / n_interfaces) if n_interfaces else 0.0
        aa = residue_identity.get(key, "")
        residue_frequencies.append({
            "unp_accession": acc,
            "role": int(role),
            "unp_residue_number": int(pos),
            "unp_residue_label": f"{aa}{pos}",
            "n_interfaces": int(count),
            "frequency": round(frac, 4),
            "conservation_level": _conservation_level(frac),
        })
    residue_frequencies.sort(
        key=lambda r: (r["role"], r["unp_accession"], r["unp_residue_number"])
    )

    contact_frequencies = []
    for c, count in contact_counts.items():
        if typed_contacts:
            (k1, k2, bond) = c
        else:
            (k1, k2) = c
            bond = None
        frac = (count / n_interfaces) if n_interfaces else 0.0
        entry = {
            "partner_1": _residue_dict(k1, residue_identity),
            "partner_2": _residue_dict(k2, residue_identity),
            "n_interfaces": int(count),
            "frequency": round(frac, 4),
        }
        if typed_contacts:
            entry["bond_type"] = bond
        contact_frequencies.append(entry)
    contact_frequencies.sort(key=lambda r: (
        r["partner_1"]["role"],
        r["partner_1"]["unp_residue_number"],
        r["partner_2"]["role"],
        r["partner_2"]["unp_residue_number"],
        r.get("bond_type") or "",
    ))

    partners = [partner_meta[key] for key in sorted(partner_meta.keys(), key=lambda kv: kv[1])]

    export = {
        "metadata": {
            "complex_id": complex_id,
            "complex_name": complex_name,
            "oligomeric_state": oligomeric_state,
            "generated_on": date.today().isoformat(),
            "n_interfaces": n_interfaces,
            "representation": "UniProt residue-level interface interactions",
            "contact_definition": (
                "PISA-derived residue-residue interface interactions mapped to "
                "UniProt residue positions"
            ),
            "typed_contacts": bool(typed_contacts),
        },
        "partners": partners,
        "residue_frequencies": residue_frequencies,
        "contact_frequencies": contact_frequencies,
    }

    output_dir = Path(getattr(config, "output_dir", None) or ".")
    output_dir.mkdir(parents=True, exist_ok=True)
    cid = complex_id or "complex"
    output_path = output_dir / f"{cid}.json"
    with output_path.open("w") as f:
        json.dump(export, f, indent=2)

    return export


def _conservation_level(frac: float) -> str:
    if frac >= CONSERVATION_LEVEL_STRONG:
        return "strong"
    if frac >= CONSERVATION_LEVEL_MEDIUM:
        return "medium"
    if frac >= CONSERVATION_LEVEL_WEAK:
        return "weak"
    return "rare"


def _residue_dict(key: tuple, residue_identity: dict) -> dict:
    acc, pos, role = key
    aa = residue_identity.get(key, "")
    return {
        "unp_accession": acc,
        "role": int(role),
        "unp_residue_number": int(pos),
        "unp_residue_label": f"{aa}{pos}",
    }


def _build_partner_metadata(
    records: list[InterfaceRecord],
    partner_map: dict[tuple[str, int], str],
    complex_details: dict | None,
) -> dict[tuple[str, int], dict]:
    """Per (unp_accession, role) -> {role, unp_accession, gene_name, name, label}.

    Roles observed in `records` and in `partner_map` are unioned so homodimers
    (one accession, two roles) and heterodimers both end up with both roles
    present. `gene_name` and `name` come from `complex_details["participants"]`
    when supplied; otherwise they are None and `label` is derived from
    `partner_map` (with " (role 1)" / " (role 2)" suffixes stripped).
    """
    seen: set[tuple[str, int]] = set()
    for r in records:
        if r.unp_accession_1:
            seen.add((r.unp_accession_1, 1))
        if r.unp_accession_2:
            seen.add((r.unp_accession_2, 2))
    for key in (partner_map or {}).keys():
        seen.add(key)

    rich: dict[str, dict] = {}
    if complex_details:
        for p in complex_details.get("participants") or []:
            if p.get("accession_type") != "UniProt":
                continue
            acc = p.get("accession")
            if not acc:
                continue
            rich[acc] = {
                "gene_name": p.get("gene_name"),
                "name": p.get("name"),
            }

    out: dict[tuple[str, int], dict] = {}
    for (acc, role) in sorted(seen, key=lambda kv: kv[1]):
        info = rich.get(acc, {})
        gene = info.get("gene_name")
        name = info.get("name")
        if gene:
            label = gene
        elif name:
            label = name
        elif partner_map and (acc, role) in partner_map:
            label = (
                partner_map[(acc, role)]
                .replace(" (role 1)", "")
                .replace(" (role 2)", "")
            )
        else:
            label = acc
        out[(acc, role)] = {
            "role": int(role),
            "unp_accession": acc,
            "gene_name": gene,
            "name": name,
            "label": label,
        }
    return out


def export_preview(export_data: dict, top_n: int = 10) -> dict[str, pd.DataFrame]:
    """Top residues and contacts from an export, sorted by frequency.

    The exported JSON itself uses deterministic role and position ordering;
    these tables are sorted by frequency for reading.
    """
    residues = (
        pd.DataFrame(export_data["residue_frequencies"])
        .sort_values("frequency", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    contacts = pd.DataFrame([
        {
            "partner_1": c["partner_1"]["unp_residue_label"],
            "partner_2": c["partner_2"]["unp_residue_label"],
            "bond_type": c.get("bond_type", ""),
            "n_interfaces": c["n_interfaces"],
            "frequency": c["frequency"],
        }
        for c in export_data["contact_frequencies"]
    ])
    contacts = (
        contacts.sort_values("frequency", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    return {"residues": residues, "contacts": contacts}


def rewiring_table(
    records: list[InterfaceRecord],
    cluster_result: ClusterResult,
    partner_map: dict[tuple[str, int], str] | None = None,
    top_n: int = 20,
) -> tuple[pd.DataFrame, str]:
    """`compare_cluster_contacts` plus a description of what was compared.

    Returns (table, message). When fewer than two non-singleton states exist
    the table is empty and the message explains why.
    """
    try:
        table = compare_cluster_contacts(
            records, cluster_result, partner_map=partner_map, top_n=top_n,
        )
    except ValueError as exc:
        return pd.DataFrame(), f"No comparison available: {exc}"
    message = (
        f"State {table.attrs.get('cluster_a')} (n={table.attrs.get('size_a')}) versus "
        f"state {table.attrs.get('cluster_b')} (n={table.attrs.get('size_b')}), "
        f"sorted by absolute fraction difference, top {top_n}."
    )
    return table, message
