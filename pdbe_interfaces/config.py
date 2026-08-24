"""Analysis settings and logging setup for the notebook.

Assembly instances carrying a macromolecule beyond the two UniProt-mapped
components are always excluded and this is not configurable: with more than two
components the correspondence between chains cannot be determined, so their
interfaces are not comparable across structures.

`Config` holds every setting the notebook exposes. Defaults reproduce the
working example. Override any field at construction:

    cfg = Config(identifier="6m0j", max_resolution=2.5)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, fields

import pandas as pd

# Components excluded from the ligand analysis by default: counter-ions,
# buffers, cryoprotectants and reducing agents, which are artefacts of
# crystallisation rather than functional ligands. Remove an id from this set to
# retain it, for example a catalytic metal.
DEFAULT_LIGAND_BLOCKLIST = frozenset({
    "CL", "NA", "MG", "ZN", "CA", "K",
    "SO4", "GOL", "HOH", "EDO", "PEG", "MPD", "ACT",
    "SGM", "DTT", "BME", "TRS",
})


@dataclass
class Config:
    """Settings for one run of the analysis.

    identifier:
        PDB entry id (for example "11gl") or PDB complex id (for example
        "PDB-CPX-172174"). An entry id is resolved to the complex it belongs
        to, and the analysis then covers every deposited structure of that
        complex. Dimers only.
    mutation_type_filter:
        Mutation types retained from the annotation API. The default keeps
        engineered mutations only; add "Conflict" to include natural sequence
        variants.
    ligand_blocklist:
        Chemical component ids excluded from the ligand analysis.
    drop_carbohydrate_polymers:
        Exclude carbohydrate polymers (glycans) from the ligand analysis.
    max_resolution:
        Exclude assemblies above this resolution in Angstroms, and those with
        no reported resolution. None disables the filter.
    max_entries:
        Restrict the analysis to the first N PDB entries in alphabetical order.
        Intended for quick checks; a subset biases the clustering and every
        conservation fraction. None uses the full dataset.
    max_workers:
        Threads used for the per-entry and per-ligand retrieval in Phase 1.
        The ligand endpoint is called once per ligand instance, so this is the
        main determinant of retrieval time.
    conservation_threshold:
        Fraction of interfaces in which a residue or contact must occur to be
        reported as conserved.
    cluster_distance_cut:
        Height at which the dendrogram is cut into clusters, in units of
        1 - Jaccard similarity. Sets how many interface interaction states are
        reported, so inspect the dendrogram before accepting the default.
    log_level:
        Logging level for the retrieval and processing steps.
    output_dir:
        Destination for the Phase 7 JSON export. None writes to the working
        directory.
    """

    identifier: str = "11gl"
    mutation_type_filter: tuple = ("Engineered mutation",)
    ligand_blocklist: frozenset = DEFAULT_LIGAND_BLOCKLIST
    drop_carbohydrate_polymers: bool = False
    max_resolution: float | None = None
    max_entries: int | None = None
    max_workers: int = 8
    conservation_threshold: float = 0.8
    cluster_distance_cut: float = 0.6
    log_level: str = "INFO"
    output_dir: str | None = "interface_frequencies"

    def describe(self) -> pd.DataFrame:
        """Settings as a two-column table, for display in the notebook."""
        rows = []
        for f in fields(self):
            value = getattr(self, f.name)
            if isinstance(value, frozenset):
                value = ", ".join(sorted(value))
            elif isinstance(value, tuple):
                value = ", ".join(str(v) for v in value)
            rows.append({"setting": f.name, "value": value})
        return pd.DataFrame(rows).set_index("setting")


def setup_logging(level: str = "INFO") -> None:
    """Configure logging for the notebook and quieten urllib3."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        force=True,
    )
    logging.getLogger("urllib3").setLevel(logging.WARNING)
