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
"""Mol* / MolViewSpec rendering helpers.

Renders one interface at a time. Each partner chain is shown as a semi-
transparent molecular surface (distinct colours), with interface residues
highlighted as opaque ball-and-stick that protrude through the surface. At-
interface mutation residues are overlaid red. The returned widget displays
natively in Jupyter / Colab.
"""

from __future__ import annotations

from typing import Optional

from molviewspec import create_builder, molstar_notebook

from pdbe_interfaces.annotations import AnnotationOverlap
from pdbe_interfaces.representation import InterfaceRecord

PDBE_CIF_BASE = "https://www.ebi.ac.uk/pdbe/entry-files/download"

PARTNER_COLOURS = ("cornflowerblue", "lightcoral", "lightgreen", "khaki")
SURFACE_COLOUR = "lightgray"  # neutral backdrop so interface residues pop


def visualize_interface(
    record: InterfaceRecord,
    overlap: Optional[AnnotationOverlap] = None,
    width: int | str = 900,
    height: int | str = 500,
    surface_opacity: float = 0.4,
    surface_color: str = SURFACE_COLOUR,
    zoom_radius_factor: float = 1.3,
    stick_size_factor: float = 0.5,
    title: Optional[str] = None,
):
    """Render one interface in Mol*, camera focused on the interface.

    Encoding:
    - Both partner chains as a semi-transparent molecular surface in a
      neutral colour (`surface_color`, default 'lightgray'). The neutral
      backdrop makes the colored interface residues stand out. Pass a
      different value to change it; pass an item of `PARTNER_COLOURS` to
      get the older per-chain coloring back manually.
    - Interface residues drawn as ball-and-stick using the chain's colour
      (cornflowerblue / lightcoral); chain identity is read off the
      stick colour, not the surface colour.
    - `stick_size_factor` (default 0.5) scales both ball radii and stick
      radii proportionally; MVS exposes a single size_factor that
      affects both.
    - If `overlap` is provided, at-interface mutation residues are overlaid
      red on top of the chain colour.
    - Camera focused on the union of all interface residues with a small
      padding factor (`zoom_radius_factor`, default 1.3 = ~30% extra space
      around the tightest fit).
    - Hydrogens are not drawn (`ignore_hydrogens=True`) for cleaner output.

    Renders by side-effect (calls IPython.display in Jupyter / Colab).
    Returns None.
    """
    cif_url = f"{PDBE_CIF_BASE}/{record.pdb_id}.cif"
    builder = create_builder()
    structure = (
        builder.download(url=cif_url)
        .parse(format="mmcif")
        .model_structure()
    )

    chains: set[str] = set()
    interface_residues: set[tuple[str, int]] = set()
    for (a, b, _) in record.author_pairs:
        chains.add(a[1])
        chains.add(b[1])
        interface_residues.add((a[1], a[2]))
        interface_residues.add((b[1], b[2]))

    chain_list = sorted(chains)
    chain_colour = {
        chain: PARTNER_COLOURS[i % len(PARTNER_COLOURS)]
        for i, chain in enumerate(chain_list)
    }

    # Semi-transparent gray surface for both chains, a neutral backdrop so
    # the colored interface residues pop visually.
    for chain in chain_list:
        (
            structure.component(selector={"auth_asym_id": chain})
            .representation(type="surface")
            .color(color=surface_color)
            .opacity(opacity=surface_opacity)
        )

    # Ball-and-stick on interface residues, one consolidated component
    # per chain so the spec stays compact.
    by_chain: dict[str, list[int]] = {}
    for chain, resnum in interface_residues:
        by_chain.setdefault(chain, []).append(int(resnum))
    for chain, resnums in by_chain.items():
        selectors = [
            {"auth_asym_id": chain, "auth_seq_id": r}
            for r in sorted(resnums)
        ]
        (
            structure.component(selector=selectors)
            .representation(
                type="ball_and_stick",
                size_factor=stick_size_factor,
                ignore_hydrogens=True,
            )
            .color(color=chain_colour[chain])
        )

    # Red overlay on at-interface mutations.
    if overlap is not None and not overlap.mutations.empty:
        muts = overlap.mutations[
            (overlap.mutations["pdb_id"] == record.pdb_id)
            & (overlap.mutations["assembly_id"] == record.assembly_id)
            & (overlap.mutations["interface_id"] == record.interface_id)
        ]
        mut_selectors = [
            {
                "auth_asym_id": row["auth_asym_id"],
                "auth_seq_id": int(row["auth_seq_id"]),
            }
            for _, row in muts.iterrows()
        ]
        if mut_selectors:
            (
                structure.component(selector=mut_selectors)
                .representation(
                    type="ball_and_stick",
                    size_factor=stick_size_factor,
                    ignore_hydrogens=True,
                )
                .color(color="red")
            )

    # Camera focus on the union of all interface residues.
    all_interface_selectors = [
        {"auth_asym_id": c, "auth_seq_id": r}
        for (c, r) in sorted(interface_residues)
    ]
    if all_interface_selectors:
        structure.component(selector=all_interface_selectors).focus(
            radius_factor=zoom_radius_factor,
        )

    return molstar_notebook(builder, width=width, height=height)


def _pick_representative(
    members: list[InterfaceRecord],
    assembly_metadata: Optional[dict] = None,
) -> InterfaceRecord:
    """Pick a cluster representative.

    If `assembly_metadata` is provided, choose the member with the lowest
    (best) resolution; ties broken by pdb_id. Members without a resolution
    value are de-prioritised. When no metadata is given, falls back to the
    first member by record order (alphabetic by pdb_id in practice).
    """
    if not members:
        raise ValueError("No members in cluster")
    if not assembly_metadata:
        return sorted(members, key=lambda r: (r.pdb_id, r.assembly_id, r.interface_id))[0]

    def score(r: InterfaceRecord):
        meta = assembly_metadata.get((r.pdb_id, r.assembly_id), {}) or {}
        res = meta.get("resolution")
        return (
            res if res is not None else float("inf"),
            r.pdb_id,
            r.assembly_id,
            r.interface_id,
        )

    return sorted(members, key=score)[0]


def visualize_cluster_representative(
    cluster_id: int,
    records: list[InterfaceRecord],
    cluster_result,
    overlap: Optional[AnnotationOverlap] = None,
    width: int | str = 900,
    height: int | str = 500,
    surface_opacity: float = 0.4,
    surface_color: str = SURFACE_COLOUR,
    zoom_radius_factor: float = 1.3,
    stick_size_factor: float = 0.5,
    assembly_metadata: Optional[dict] = None,
):
    """Visualise a representative interface from the named cluster.

    When `assembly_metadata` is provided, picks the member with the best
    (lowest) resolution; otherwise picks the first member by deterministic
    sort (pdb_id, assembly_id, interface_id).
    """
    cluster_of = {
        r.key: int(cluster_result.flat_assignment[i])
        for i, r in enumerate(records)
        if i < len(cluster_result.flat_assignment)
    }
    members = [r for r in records if cluster_of.get(r.key) == cluster_id]
    if not members:
        raise ValueError(f"Cluster {cluster_id} has no members")
    rep = _pick_representative(members, assembly_metadata)
    return visualize_interface(
        rep, overlap=overlap, width=width, height=height,
        surface_opacity=surface_opacity,
        surface_color=surface_color,
        zoom_radius_factor=zoom_radius_factor,
        stick_size_factor=stick_size_factor,
    )


def visualize_clusters_grid(
    records: list[InterfaceRecord],
    cluster_result,
    overlap: Optional[AnnotationOverlap] = None,
    min_cluster_size: int = 2,
    max_clusters: int = 8,
    width: int | str = 900,
    height: int | str = 450,
    surface_opacity: float = 0.4,
    surface_color: str = SURFACE_COLOUR,
    zoom_radius_factor: float = 1.3,
    stick_size_factor: float = 0.5,
    assembly_metadata: Optional[dict] = None,
):
    """Render one representative per cluster, stacked vertically.

    Loops through clusters that meet `min_cluster_size`, sorted by size
    descending, taking up to `max_clusters`. For each, prints a header line
    with cluster id, size, and the representative PDB ID, then renders the
    Mol* viewer.

    Each viewer adds a few hundred kB to the notebook, so the defaults
    (>=2 members, max 8 viewers) keep the output tractable. Override
    explicitly if you want to see more.
    """
    from collections import Counter
    from IPython.display import Markdown, display

    cluster_of = {
        r.key: int(cluster_result.flat_assignment[i])
        for i, r in enumerate(records)
        if i < len(cluster_result.flat_assignment)
    }
    sizes = Counter(cluster_of.values())
    clusters_to_show = [
        cid for cid, n in sizes.most_common()
        if n >= min_cluster_size
    ][:max_clusters]

    if not clusters_to_show:
        display(Markdown(
            f"_No cluster meets the threshold (`min_cluster_size={min_cluster_size}`)._"
        ))
        return

    display(Markdown(
        f"**Showing {len(clusters_to_show)} cluster representatives** "
        f"(min_cluster_size={min_cluster_size}, max_clusters={max_clusters}). "
        f"Skipped: {len(sizes) - len(clusters_to_show)} smaller clusters."
    ))

    for cid in clusters_to_show:
        members = [r for r in records if cluster_of.get(r.key) == cid]
        rep = _pick_representative(members, assembly_metadata)
        rep_meta = (assembly_metadata or {}).get((rep.pdb_id, rep.assembly_id), {}) or {}
        rep_res = rep_meta.get("resolution")
        res_str = f", {rep_res:.2f} Å" if rep_res is not None else ""
        display(Markdown(
            f"### Cluster {cid} (n={sizes[cid]}), representative `{rep.pdb_id}` "
            f"(assembly {rep.assembly_id}, interface {rep.interface_id}{res_str})"
        ))
        visualize_interface(
            rep, overlap=overlap, width=width, height=height,
            surface_opacity=surface_opacity,
            surface_color=surface_color,
            zoom_radius_factor=zoom_radius_factor,
            stick_size_factor=stick_size_factor,
        )
