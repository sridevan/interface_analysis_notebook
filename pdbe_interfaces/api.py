"""PDBe v2 API calls. Returns parsed JSON. No domain logic.

Robustness behaviour:
- Transient errors (5xx, 429, ConnectionError, Timeout) trigger up to MAX_ATTEMPTS
  total attempts with exponential backoff. Retries are logged at WARNING.
- Persistent 4xx errors raise immediately (except 404 on endpoints called with
  allow_404=True, which return an empty result).
- The two batch-POST endpoints (mutations, modifications) chunk their PDB ID
  list into POST_BATCH_SIZE blocks and merge results.
"""

from __future__ import annotations

import json
import logging
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import requests

log = logging.getLogger(__name__)

BASE = "https://www.ebi.ac.uk/pdbe/api/v2"
TIMEOUT = 30

MAX_ATTEMPTS = 5
BACKOFF_BASE_S = 1.0
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})

POST_BATCH_SIZE = 50

# Concurrency for the fan-out helpers. The per-ligand interaction endpoint is
# called once per ligand instance, so a complex with ~100 entries runs into the
# thousands of requests; serially that dominates the whole workflow. Kept
# modest to stay polite to the EBI API; 429s are retried with backoff anyway.
DEFAULT_MAX_WORKERS = 8

# One requests.Session per worker thread. Sessions are not documented as
# thread-safe, but keeping them thread-local gives each worker HTTP keep-alive
# (no TCP + TLS handshake per call) without sharing connection pools.
_thread_local = threading.local()


def _session() -> requests.Session:
    session = getattr(_thread_local, "session", None)
    if session is None:
        session = requests.Session()
        _thread_local.session = session
    return session


def _retry_request(method: str, url: str, **kwargs) -> requests.Response:
    """Issue an HTTP request with retries on transient errors.

    Retries on ConnectionError, Timeout, and HTTP statuses in RETRY_STATUSES.
    Does not retry on other 4xx or on a successful response.
    """
    last_exc: Exception | None = None
    resp: requests.Response | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            resp = _session().request(method, url, timeout=TIMEOUT, **kwargs)
        except (requests.ConnectionError, requests.Timeout) as e:
            last_exc = e
            if attempt < MAX_ATTEMPTS:
                wait = BACKOFF_BASE_S * (2 ** (attempt - 1))
                log.warning(
                    "Network error on %s %s (attempt %d/%d): %s. Retrying in %.1fs.",
                    method, url, attempt, MAX_ATTEMPTS, e, wait,
                )
                time.sleep(wait)
                continue
            raise
        if resp.status_code in RETRY_STATUSES and attempt < MAX_ATTEMPTS:
            wait = BACKOFF_BASE_S * (2 ** (attempt - 1))
            log.warning(
                "HTTP %d from %s %s (attempt %d/%d). Retrying in %.1fs.",
                resp.status_code, method, url, attempt, MAX_ATTEMPTS, wait,
            )
            time.sleep(wait)
            continue
        return resp
    if last_exc is not None:
        raise last_exc
    assert resp is not None
    return resp


def _get_json(
    url: str,
    *,
    allow_404: bool = False,
    tolerate_failure: bool = False,
) -> Any:
    """GET and return parsed JSON.

    `allow_404`: 404 returns None and logs a WARNING instead of raising.
    `tolerate_failure`: any exception (after retries exhausted, including
        non-404 4xx and 5xx) returns None and logs a WARNING. Use for
        best-effort calls where one failure shouldn't kill the whole job
        (per-PDB / per-ligand fetches).
    """
    log.info("GET %s", url)
    try:
        resp = _retry_request("GET", url)
        if allow_404 and resp.status_code == 404:
            log.warning("404 from %s, returning empty", url)
            return None
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        if tolerate_failure:
            log.warning("Tolerated failure on GET %s: %s", url, e)
            return None
        raise


def _post_json(url: str, body: Any, *, allow_404: bool = False) -> Any:
    log.info("POST %s", url)
    resp = _retry_request(
        "POST", url,
        data=json.dumps(body),
        headers={"Content-Type": "application/json"},
    )
    if allow_404 and resp.status_code == 404:
        log.warning("404 from %s, returning empty", url)
        return {}
    resp.raise_for_status()
    return resp.json()


def _batched_post(url: str, pdb_ids: list[str]) -> dict:
    """POST in chunks of POST_BATCH_SIZE PDB IDs; merge per-chunk dicts."""
    merged: dict = {}
    n = len(pdb_ids)
    if n == 0:
        return merged
    for start in range(0, n, POST_BATCH_SIZE):
        chunk = pdb_ids[start:start + POST_BATCH_SIZE]
        end = start + len(chunk)
        log.info("Batched POST %s: chunk %d-%d of %d", url, start + 1, end, n)
        result = _post_json(url, ",".join(chunk), allow_404=True)
        if result:
            merged.update(result)
    return merged


COMPLEX_ID_RE = re.compile(r"^PDB-CPX-\d+$", re.IGNORECASE)


def fetch_complex_details(complex_id: str) -> dict:
    url = f"{BASE}/complex/details/{complex_id}?id_type=pdb_complex_id"
    data = _get_json(url)
    if not data or complex_id not in data or not data[complex_id]:
        raise ValueError(f"Complex ID {complex_id} did not resolve")
    return data[complex_id][0]


def fetch_complexes_for_pdb_id(pdb_id: str) -> list[dict]:
    """Complexes a PDB entry participates in, via `id_type=pdb_id`.

    The pdb_id-keyed response is a thin record (name, pdb_complex_id,
    participants, assemblies of that entry only). It carries neither
    `total_chains` nor per-assembly `bound_macromolecules`. Use it to resolve
    the complex ID, then call `fetch_complex_details` for the full record.
    """
    pdb_id = pdb_id.strip().lower()
    url = f"{BASE}/complex/details/{pdb_id}?id_type=pdb_id"
    data = _get_json(url, allow_404=True)
    if not data or pdb_id not in data or not data[pdb_id]:
        raise ValueError(
            f"PDB entry {pdb_id} did not resolve to any complex. Check the ID, "
            "or pass a PDB complex ID (PDB-CPX-...) directly."
        )
    return data[pdb_id]


def resolve_complex_id(
    identifier: str, require_dimer: bool = False,
) -> tuple[str, dict]:
    """Accept a PDB entry ID or a PDB complex ID; return (complex_id, details).

    `PDB-CPX-140195` is used as-is. Anything else is treated as a PDB entry ID
    and looked up with `id_type=pdb_id`.

    Roughly 4% of PDB entries map to more than one complex, because separate
    assemblies of one entry can have different compositions. The endpoint
    returns all of them, each carrying a `preferred_assembly` flag. Selection
    order: the complex holding the preferred assembly if it is a dimer, else
    the only dimer among the candidates. If neither rule resolves it, the
    ambiguity is raised with the candidates listed so the caller can pass a
    complex ID explicitly.

    `require_dimer` raises if the resolved complex is not a dimer. The rest of
    this workflow supports dimers only, a limitation of the
    interface_interactions endpoint rather than a design choice.
    """
    ident = str(identifier).strip()
    if not ident:
        raise ValueError("Empty identifier. Pass a PDB entry ID or a PDB complex ID.")
    if COMPLEX_ID_RE.match(ident):
        complex_id = ident.upper()
        details = fetch_complex_details(complex_id)
        if require_dimer:
            _require_dimer(ident, complex_id, details)
        return complex_id, details

    candidates = fetch_complexes_for_pdb_id(ident)
    complex_ids = []
    for c in candidates:
        cid = c.get("pdb_complex_id")
        if cid and cid not in complex_ids:
            complex_ids.append(cid)
    if not complex_ids:
        raise ValueError(f"PDB entry {ident} returned no pdb_complex_id")

    if len(complex_ids) == 1:
        complex_id = complex_ids[0]
        log.info("Resolved PDB entry %s to complex %s", ident, complex_id)
        details = fetch_complex_details(complex_id)
        if require_dimer:
            _require_dimer(ident, complex_id, details)
        return complex_id, details

    # Prefer the complex holding the entry's preferred assembly, which is the
    # composition PDBe treats as canonical for that entry. Fall back to the
    # only dimer if the preferred assembly's complex is not one.
    preferred = [
        c.get("pdb_complex_id") for c in candidates
        if any(a.get("preferred_assembly") for a in (c.get("assemblies") or []))
    ]
    details_by_id = {cid: fetch_complex_details(cid) for cid in complex_ids}
    dimers = [cid for cid, d in details_by_id.items() if d.get("total_chains") == 2]

    if len(preferred) == 1 and preferred[0] in dimers:
        complex_id = preferred[0]
        others = [cid for cid in complex_ids if cid != complex_id]
        log.info(
            "PDB entry %s maps to %d complexes; picked %s, which holds the "
            "preferred assembly (skipped %s)",
            ident, len(complex_ids), complex_id, ", ".join(others),
        )
        return complex_id, details_by_id[complex_id]

    if len(dimers) == 1:
        complex_id = dimers[0]
        others = [cid for cid in complex_ids if cid != complex_id]
        log.info(
            "PDB entry %s maps to %d complexes; the preferred assembly's complex "
            "is not a dimer, so picked the only dimer %s (skipped %s)",
            ident, len(complex_ids), complex_id, ", ".join(others),
        )
        return complex_id, details_by_id[complex_id]

    listing = "\n".join(
        f"  {cid}  total_chains={details_by_id[cid].get('total_chains')}  "
        f"{details_by_id[cid].get('oligomeric_state')}  {details_by_id[cid].get('name')}"
        for cid in complex_ids
    )
    raise ValueError(
        f"PDB entry {ident} maps to {len(complex_ids)} complexes and the dimer is "
        f"ambiguous. Set `identifier` to one of these complex IDs:\n{listing}"
    )


def fetch_interface_interactions(complex_id: str) -> list[dict]:
    url = f"{BASE}/complex/interface_interactions/{complex_id}"
    data = _get_json(url)
    if not data or complex_id not in data:
        raise ValueError(f"No interfaces returned for {complex_id}")
    interfaces = data[complex_id]
    if not interfaces:
        raise ValueError(f"No interfaces returned for {complex_id}")
    return interfaces


def fetch_mutations(pdb_ids: list[str]) -> dict:
    """Returns {pdb_id: [mutation records]}. Empty when no entries have mutations."""
    return _batched_post(f"{BASE}/pdb/entry/mutated_AA_or_NA", pdb_ids)


def fetch_modifications(pdb_ids: list[str]) -> dict:
    """Returns {pdb_id: [modification records]}. Empty when no entries have modifications."""
    return _batched_post(f"{BASE}/pdb/entry/modified_AA_or_NA", pdb_ids)


def fetch_bound_molecules(pdb_id: str) -> list[dict]:
    """Best-effort fetch. 404 or persistent failure returns []."""
    url = f"{BASE}/pdb/bound_molecules/{pdb_id}"
    data = _get_json(url, allow_404=True, tolerate_failure=True)
    if not data or pdb_id not in data:
        return []
    return data[pdb_id]


def fetch_ligand_interactions(
    pdb_id: str, chain_id: str, author_residue_number: int
) -> list[dict]:
    """Best-effort fetch of {ligand, interactions} bm-level objects.

    404 (no interactions on file for this instance) returns []. Persistent
    failures (5xx after retries exhausted) also return [] with a logged
    WARNING, since a single missing ligand should not stop the whole analysis.
    """
    url = (
        f"{BASE}/pdb/bound_ligand_interactions/{pdb_id}/{chain_id}/"
        f"{author_residue_number}?preserve_case=false"
    )
    data = _get_json(url, allow_404=True, tolerate_failure=True)
    if not data or pdb_id not in data:
        return []
    return data[pdb_id]


def _require_dimer(identifier: str, complex_id: str, details: dict) -> None:
    """Raise unless the resolved complex has exactly two chains."""
    if details.get("total_chains") != 2:
        raise ValueError(
            f"{identifier} resolves to complex {complex_id}, which has "
            f"{details.get('total_chains')} chains ({details.get('oligomeric_state')}). "
            "This workflow supports dimers only, with both components mapped to "
            "UniProt: with more than two components the correspondence between "
            "chains cannot be determined, so the interfaces are not comparable "
            "across structures."
        )


def _parallel_map(fn, items: list, max_workers: int, what: str) -> list:
    """Map `fn` over `items` on a thread pool, preserving input order.

    Used only for the best-effort fetchers (`fetch_bound_molecules`,
    `fetch_ligand_interactions`), which swallow their own failures and return
    an empty result, so one failing entry cannot take down the pool.
    """
    if not items:
        return []
    workers = max(1, min(max_workers, len(items)))
    if workers == 1:
        return [fn(item) for item in items]
    log.info("Fetching %d %s with %d workers", len(items), what, workers)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fn, items))


def fetch_bound_molecules_many(
    pdb_ids: list[str], max_workers: int = DEFAULT_MAX_WORKERS,
) -> dict[str, list[dict]]:
    """Concurrent `fetch_bound_molecules` over many entries -> {pdb_id: records}."""
    results = _parallel_map(
        fetch_bound_molecules, list(pdb_ids), max_workers, "bound_molecules responses",
    )
    return dict(zip(pdb_ids, results))


def fetch_ligand_interactions_many(
    ligand_keys: list[tuple[str, str, int]], max_workers: int = DEFAULT_MAX_WORKERS,
) -> list[list[dict]]:
    """Concurrent `fetch_ligand_interactions` over (pdb_id, chain_id, author_residue_number).

    Returns one result list per input key, in input order.
    """
    return _parallel_map(
        lambda key: fetch_ligand_interactions(*key),
        list(ligand_keys),
        max_workers,
        "ligand interaction responses",
    )
