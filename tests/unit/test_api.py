"""HTTP wrapper behaviour in `pdbe_interfaces.api`, with every request faked.

Tests the retry policy, the two 404 policies, the tolerated-failure path,
batched POST chunking and the ordering guarantee of the parallel fetchers.
Nothing here talks to PDBe.
"""

from __future__ import annotations

import json
import time

import pytest
import requests

from pdbe_interfaces import api


class FakeResponse:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}", response=self)


class FakeHTTP:
    """Stand-in for `requests.Session.request`.

    Either queue responses (or exceptions to raise) in order with `respond`,
    or set `handler(method, url, kwargs)` to compute them.
    """

    def __init__(self):
        self.queue: list = []
        self.calls: list[tuple[str, str, dict]] = []
        self.handler = None

    def respond(self, *items):
        self.queue.extend(items)

    def __call__(self, method, url, **kwargs):
        self.calls.append((method, url, kwargs))
        if self.handler is not None:
            return self.handler(method, url, kwargs)
        item = self.queue.pop(0)
        if isinstance(item, BaseException):
            raise item
        return item


@pytest.fixture
def http(monkeypatch):
    fake = FakeHTTP()
    monkeypatch.setattr(requests.Session, "request", fake)
    monkeypatch.setattr(api, "BACKOFF_BASE_S", 0.0)   # sleep(0) between retries
    return fake


def test_network_is_blocked_by_default():
    with pytest.raises(RuntimeError, match="Network access is blocked"):
        requests.get("https://www.ebi.ac.uk/pdbe/api/v2/complex/details/11gl")


def test_successful_get_is_parsed(http):
    http.respond(FakeResponse(200, {"1abc": [{"composition": {}}]}))
    assert api.fetch_bound_molecules("1abc") == [{"composition": {}}]
    method, url, kwargs = http.calls[0]
    assert (method, url) == ("GET", f"{api.BASE}/pdb/bound_molecules/1abc")
    assert kwargs["timeout"] == api.TIMEOUT


# --- 404 policy -------------------------------------------------------------


def test_404_is_tolerated_on_best_effort_endpoints(http):
    http.respond(FakeResponse(404))
    assert api.fetch_ligand_interactions("1abc", "A", 301) == []
    assert len(http.calls) == 1   # no retry on 404


def test_critical_endpoint_fails_loudly_on_404_or_empty_result(http):
    http.respond(FakeResponse(404))
    with pytest.raises(requests.HTTPError):
        api.fetch_interface_interactions("PDB-CPX-1")
    assert len(http.calls) == 1

    http.respond(FakeResponse(200, {"PDB-CPX-1": []}))
    with pytest.raises(ValueError, match="No interfaces returned"):
        api.fetch_interface_interactions("PDB-CPX-1")


# --- retry policy -----------------------------------------------------------


def test_transient_status_is_retried_then_succeeds(http):
    http.respond(FakeResponse(503), FakeResponse(200, {"PDB-CPX-1": [{"x": 1}]}))
    assert api.fetch_interface_interactions("PDB-CPX-1") == [{"x": 1}]
    assert len(http.calls) == 2


def test_network_exception_is_retried_then_succeeds(http):
    http.respond(requests.Timeout("slow"), requests.ConnectionError("down"),
                 FakeResponse(200, {"PDB-CPX-1": [{"x": 1}]}))
    assert api.fetch_interface_interactions("PDB-CPX-1") == [{"x": 1}]
    assert len(http.calls) == 3


def test_exhausted_retries_raise_on_critical_endpoint(http):
    http.respond(*[FakeResponse(503)] * api.MAX_ATTEMPTS)
    with pytest.raises(requests.HTTPError):
        api.fetch_interface_interactions("PDB-CPX-1")
    assert len(http.calls) == api.MAX_ATTEMPTS

    http.calls.clear()
    http.respond(*[requests.Timeout("slow")] * api.MAX_ATTEMPTS)
    with pytest.raises(requests.Timeout):
        api.fetch_interface_interactions("PDB-CPX-1")
    assert len(http.calls) == api.MAX_ATTEMPTS


def test_exhausted_retries_are_tolerated_on_best_effort_endpoint(http):
    http.respond(*[FakeResponse(503)] * api.MAX_ATTEMPTS)
    assert api.fetch_bound_molecules("1abc") == []
    assert len(http.calls) == api.MAX_ATTEMPTS


def test_non_retryable_client_error_raises_immediately(http):
    http.respond(FakeResponse(400))
    with pytest.raises(requests.HTTPError):
        api.fetch_interface_interactions("PDB-CPX-1")
    assert len(http.calls) == 1


# --- batched POST -----------------------------------------------------------


def test_batched_post_chunks_merges_and_tolerates_a_404_chunk(http):
    ids = [f"{i:04d}" for i in range(2 * api.POST_BATCH_SIZE + 20)]
    second_chunk = set(ids[api.POST_BATCH_SIZE:2 * api.POST_BATCH_SIZE])

    def handler(method, url, kwargs):
        assert method == "POST"
        assert url == f"{api.BASE}/pdb/entry/mutated_AA_or_NA"
        assert kwargs["headers"] == {"Content-Type": "application/json"}
        chunk = json.loads(kwargs["data"]).split(",")
        if set(chunk) == second_chunk:
            return FakeResponse(404)                # "no data for these entries"
        return FakeResponse(200, {pdb_id: [{"chain_id": "A"}] for pdb_id in chunk})

    http.handler = handler
    merged = api.fetch_mutations(ids)

    sizes = [len(json.loads(kwargs["data"]).split(",")) for (_, _, kwargs) in http.calls]
    assert sizes == [api.POST_BATCH_SIZE, api.POST_BATCH_SIZE, 20]
    assert set(merged) == set(ids) - second_chunk


def test_batched_post_with_no_ids_makes_no_request(http):
    assert api.fetch_mutations([]) == {}
    assert api.fetch_modifications([]) == {}
    assert http.calls == []


# --- parallel fetchers ------------------------------------------------------


def test_parallel_fetch_preserves_input_order(http):
    """`collect_ligand_contacts_for_entries` zips requests with responses, so
    order must survive even when the first request is the slowest."""
    keys = [("1abc", "A", 301), ("1abc", "A", 302), ("2xyz", "B", 5)]

    def handler(method, url, kwargs):
        if url.endswith("/1abc/A/301?preserve_case=false"):
            time.sleep(0.05)
            return FakeResponse(200, {"1abc": [{"ligand": "first"}]})
        if url.endswith("/1abc/A/302?preserve_case=false"):
            return FakeResponse(200, {"1abc": [{"ligand": "second"}]})
        return FakeResponse(404)

    http.handler = handler
    results = api.fetch_ligand_interactions_many(keys, max_workers=3)
    assert results == [[{"ligand": "first"}], [{"ligand": "second"}], []]
