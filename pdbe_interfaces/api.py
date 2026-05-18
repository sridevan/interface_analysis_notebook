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
import time
from typing import Any

import requests

log = logging.getLogger(__name__)

BASE = "https://www.ebi.ac.uk/pdbe/api/v2"
TIMEOUT = 30

MAX_ATTEMPTS = 5
BACKOFF_BASE_S = 1.0
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})

POST_BATCH_SIZE = 50


def _retry_request(method: str, url: str, **kwargs) -> requests.Response:
    """Issue an HTTP request with retries on transient errors.

    Retries on ConnectionError, Timeout, and HTTP statuses in RETRY_STATUSES.
    Does not retry on other 4xx or on a successful response.
    """
    last_exc: Exception | None = None
    resp: requests.Response | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            resp = requests.request(method, url, timeout=TIMEOUT, **kwargs)
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
            log.warning("404 from %s — returning empty", url)
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
        log.warning("404 from %s — returning empty", url)
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


def fetch_complex_details(complex_id: str) -> dict:
    url = f"{BASE}/complex/details/{complex_id}?id_type=pdb_complex_id"
    data = _get_json(url)
    if not data or complex_id not in data or not data[complex_id]:
        raise ValueError(f"Complex ID {complex_id} did not resolve")
    return data[complex_id][0]


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
    WARNING — a single missing ligand should not kill the whole analysis.
    """
    url = (
        f"{BASE}/pdb/bound_ligand_interactions/{pdb_id}/{chain_id}/"
        f"{author_residue_number}?preserve_case=false"
    )
    data = _get_json(url, allow_404=True, tolerate_failure=True)
    if not data or pdb_id not in data:
        return []
    return data[pdb_id]
