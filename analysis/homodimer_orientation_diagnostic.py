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
"""Diagnostic: is homodimer interface similarity sensitive to chain orientation?

The production representation keeps homodimer contacts as ordered
(role 1, role 2) tuples. If PISA lists the two copies in the opposite order
for two structures of the same interface, the tuples do not match. This script
measures how much that costs on real data by comparing, for every pair of
interface instances (A, B) of a homodimer complex:

    direct  = Jaccard(A, B)
    swapped = Jaccard(A, global_swap(B))      # copy 1 <-> copy 2 for the whole interface
    delta   = swapped - direct

and, for complexes where swapping helps, the clustering under the production
similarity against a diagnostic orientation-aware similarity
max(direct, swapped). Nothing here changes production behaviour, and the
orientation-aware similarity exists only in this script.

Makes live PDBe API calls. Not part of the pytest suite.

Usage:
    python analysis/homodimer_orientation_diagnostic.py --csv assemblies_data.csv --n 30
    python analysis/homodimer_orientation_diagnostic.py --complex PDB-CPX-172174 PDB-CPX-130029

`assemblies_data.csv` is the PDBe-KB complexes listing
(https://ftp.ebi.ac.uk/pub/databases/pdbe-kb/complexes/); a homodimer has a
one-component composition string `{accession}_2`.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pdbe_interfaces import api, representation, similarity  # noqa: E402

log = logging.getLogger("homodimer_diagnostic")

DISTANCE_CUT = 0.6                  # the notebook default
DELTA_THRESHOLDS = (0.2, 0.5)       # descriptive, for inspection only


def global_swap(pairs: set) -> set:
    """Exchange copy 1 and copy 2 for every contact of one interface."""
    return {((acc2, pos2, 1), (acc1, pos1, 2), bond)
            for ((acc1, pos1, _), (acc2, pos2, _), bond) in pairs}


def sample_homodimers(csv_path: Path, n: int, seed: int, min_asm: int, max_asm: int,
                      always: list[str]) -> list[str]:
    df = pd.read_csv(csv_path, usecols=["ASSEMBLY_STRING", "PDB_COMPLEX_ID"])
    homo = df[df.ASSEMBLY_STRING.str.match(r"^[A-Z0-9]+_2$")]
    counts = homo.groupby("PDB_COMPLEX_ID").size()
    pool = sorted(counts[(counts >= min_asm) & (counts <= max_asm)].index)
    pool = [c for c in pool if c not in always]
    random.Random(seed).shuffle(pool)
    chosen = list(always) + pool[: max(0, n - len(always))]
    log.info("Sampled %d complexes from %d homodimers with %d-%d assemblies",
             len(chosen), len(pool), min_asm, max_asm)
    return chosen


def is_homodimer(details: dict) -> bool:
    proteins = [p for p in (details.get("participants") or []) if p.get("accession_type") == "UniProt"]
    return (details.get("total_chains") == 2 and len(proteins) == 1
            and (proteins[0].get("stoichiometry") or 1) == 2)


def load_records(complex_id: str) -> tuple[dict, list[representation.InterfaceRecord]]:
    details = api.fetch_complex_details(complex_id)
    if not is_homodimer(details):
        raise ValueError(f"{complex_id} is not a homodimer: {details.get('oligomeric_state')}")
    selection = representation.select_interfaces(api.fetch_interface_interactions(complex_id), details)
    records = representation.check_partner_consistency(
        representation.build_interface_records(selection.interfaces)
    )
    return details, representation.select_comparable_records(records).records


def pairwise(complex_id: str, records: list) -> list[dict]:
    rows = []
    for (i, a), (j, b) in itertools.combinations(enumerate(records), 2):
        direct = similarity.jaccard(a.uniprot_pairs, b.uniprot_pairs)
        swapped = similarity.jaccard(a.uniprot_pairs, global_swap(b.uniprot_pairs))
        rows.append({
            "complex_id": complex_id, "instance_a": a.label(), "instance_b": b.label(),
            "n_pairs_a": len(a.uniprot_pairs), "n_pairs_b": len(b.uniprot_pairs),
            "direct": round(direct, 4), "swapped": round(swapped, 4),
            "delta": round(swapped - direct, 4),
        })
    return rows


def partition(labels: np.ndarray) -> set[frozenset]:
    groups: dict[int, set] = {}
    for idx, lab in enumerate(labels.tolist()):
        groups.setdefault(lab, set()).add(idx)
    return {frozenset(g) for g in groups.values()}


def pair_agreement(labels_a: np.ndarray, labels_b: np.ndarray) -> float:
    """Fraction of record pairs whose co-membership is the same in both partitions (Rand index)."""
    n = len(labels_a)
    if n < 2:
        return 1.0
    agree = sum(
        (labels_a[i] == labels_a[j]) == (labels_b[i] == labels_b[j])
        for i, j in itertools.combinations(range(n), 2)
    )
    return agree / (n * (n - 1) / 2)


def clustering_impact(complex_id: str, records: list) -> dict:
    sets = [r.uniprot_pairs for r in records]
    current = similarity.jaccard_similarity_matrix(sets)
    n = len(sets)
    aware = current.copy()
    for i in range(n):
        for j in range(i + 1, n):
            v = max(current[i, j], similarity.jaccard(sets[i], global_swap(sets[j])))
            aware[i, j] = aware[j, i] = v
    lab_cur = similarity.cluster_interfaces(current, DISTANCE_CUT).flat_assignment
    lab_aware = similarity.cluster_interfaces(aware, DISTANCE_CUT).flat_assignment
    sizes_cur = sorted(Counter(lab_cur.tolist()).values(), reverse=True)
    sizes_aware = sorted(Counter(lab_aware.tolist()).values(), reverse=True)
    return {
        "complex_id": complex_id,
        "n_interfaces": n,
        "clusters_current": len(sizes_cur),
        "clusters_aware": len(sizes_aware),
        "singletons_current": sizes_cur.count(1),
        "singletons_aware": sizes_aware.count(1),
        "sizes_current": ",".join(map(str, sizes_cur)),
        "sizes_aware": ",".join(map(str, sizes_aware)),
        "partition_changed": partition(lab_cur) != partition(lab_aware),
        "pair_agreement": round(pair_agreement(lab_cur, lab_aware), 3),
    }


def summarise(complex_id: str, name: str, n_interfaces: int, rows: list[dict]) -> dict:
    deltas = [r["delta"] for r in rows]
    return {
        "complex_id": complex_id,
        "name": name,
        "n_interfaces": n_interfaces,
        "n_comparisons": len(rows),
        "n_swapped_gt_direct": sum(d > 0 for d in deltas),
        "n_delta_gt_0.2": sum(d > 0.2 for d in deltas),
        "n_delta_gt_0.5": sum(d > 0.5 for d in deltas),
        "max_delta": max(deltas) if deltas else 0.0,
        "median_direct": float(np.median([r["direct"] for r in rows])) if rows else float("nan"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, help="PDBe-KB assemblies_data.csv for sampling")
    ap.add_argument("--n", type=int, default=30, help="complexes to sample (default 30)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-assemblies", type=int, default=8)
    ap.add_argument("--max-assemblies", type=int, default=60)
    ap.add_argument("--complex", nargs="*", default=[], help="complex ids to analyse (added to any sample)")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results")
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()

    logging.basicConfig(level="INFO", format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                        datefmt="%H:%M:%S", force=True)
    logging.getLogger("pdbe_interfaces.api").setLevel(logging.WARNING)
    logging.getLogger("pdbe_interfaces.representation").setLevel(logging.ERROR)

    complexes = list(args.complex)
    if args.csv:
        complexes = sample_homodimers(args.csv, args.n, args.seed, args.min_assemblies,
                                      args.max_assemblies, always=complexes)
    if not complexes:
        ap.error("give --csv for sampling and/or --complex ids")

    all_pairs: list[dict] = []
    per_complex: list[dict] = []
    impacts: list[dict] = []
    skipped: list[dict] = []
    for cid in complexes:
        try:
            details, records = load_records(cid)
        except Exception as exc:                      # noqa: BLE001, diagnostic script
            skipped.append({"complex_id": cid, "reason": str(exc)[:200]})
            log.warning("Skipped %s: %s", cid, str(exc)[:120])
            continue
        if len(records) < 2:
            skipped.append({"complex_id": cid, "reason": f"only {len(records)} comparable interface"})
            continue
        rows = pairwise(cid, records)
        all_pairs.extend(rows)
        summary = summarise(cid, details.get("name", ""), len(records), rows)
        per_complex.append(summary)
        if summary["max_delta"] > DELTA_THRESHOLDS[0]:
            impacts.append(clustering_impact(cid, records))
        log.info("%s %-40.40s interfaces=%3d pairs=%5d swapped>direct=%4d max_delta=%.2f",
                 cid, details.get("name", ""), len(records), len(rows),
                 summary["n_swapped_gt_direct"], summary["max_delta"])

    args.out.mkdir(parents=True, exist_ok=True)
    pairs_df = pd.DataFrame(all_pairs)
    pairs_df.to_csv(args.out / "pairwise.csv", index=False)
    pd.DataFrame(per_complex).to_csv(args.out / "per_complex.csv", index=False)
    pd.DataFrame(impacts).to_csv(args.out / "clustering_impact.csv", index=False)
    pd.DataFrame(skipped).to_csv(args.out / "skipped.csv", index=False)
    top = pairs_df.sort_values("delta", ascending=False).head(args.top) if not pairs_df.empty else pairs_df
    top.to_csv(args.out / "top_examples.csv", index=False)

    n_pairs = len(pairs_df)
    overall = {
        "n_complexes_analysed": len(per_complex),
        "n_complexes_skipped": len(skipped),
        "n_interfaces": int(sum(s["n_interfaces"] for s in per_complex)),
        "n_comparisons": n_pairs,
        "fraction_swapped_gt_direct": float((pairs_df["delta"] > 0).mean()) if n_pairs else None,
        "fraction_delta_gt_0.2": float((pairs_df["delta"] > 0.2).mean()) if n_pairs else None,
        "fraction_delta_gt_0.5": float((pairs_df["delta"] > 0.5).mean()) if n_pairs else None,
        "n_complexes_any_delta_gt_0.2": sum(s["n_delta_gt_0.2"] > 0 for s in per_complex),
        "n_complexes_any_delta_gt_0.5": sum(s["n_delta_gt_0.5"] > 0 for s in per_complex),
        "n_complexes_clustering_changed": sum(1 for i in impacts if i["partition_changed"]),
        "n_complexes_cluster_count_changed": sum(
            1 for i in impacts if i["clusters_current"] != i["clusters_aware"]),
        "n_complexes_singletons_reduced": sum(
            1 for i in impacts if i["singletons_aware"] < i["singletons_current"]),
        "distance_cut": DISTANCE_CUT,
    }
    (args.out / "summary.json").write_text(json.dumps(overall, indent=2))

    print("\n=== Overall ===")
    print(json.dumps(overall, indent=2))
    print("\n=== Per complex ===")
    print(pd.DataFrame(per_complex).to_string(index=False))
    if impacts:
        print("\n=== Clustering impact (complexes with any delta > 0.2) ===")
        print(pd.DataFrame(impacts).to_string(index=False))
    if not top.empty:
        print(f"\n=== Top {args.top} orientation-sensitive pairs ===")
        print(top.to_string(index=False))
    if skipped:
        print("\n=== Skipped ===")
        print(pd.DataFrame(skipped).to_string(index=False))
    print(f"\nOutputs written to {args.out}")


if __name__ == "__main__":
    main()
