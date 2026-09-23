"""Identifier resolution: `api.resolve_complex_id` with the fetchers mocked.

The PDB-entry path normalises inside `fetch_complexes_for_pdb_id`, so `_get_json`
is mocked and the URL inspected. The complex-id path normalises before calling
`fetch_complex_details`, so that function is mocked.
"""

from __future__ import annotations

import pytest

from pdbe_interfaces import api

CPX = "PDB-CPX-172174"


def _details(complex_id, total_chains=2, oligomeric_state="Homodimer"):
    return {
        "name": f"complex {complex_id}", "total_chains": total_chains,
        "oligomeric_state": oligomeric_state, "participants": [], "assemblies": [],
    }


def _entry_response(pdb_id, *complexes):
    """`complex/details?id_type=pdb_id` payload: one thin record per complex."""
    return {pdb_id: [
        {"pdb_complex_id": cid, "assemblies": [{"preferred_assembly": preferred}]}
        for cid, preferred in complexes
    ]}


@pytest.fixture
def fake_api(monkeypatch):
    """Route `_get_json` and `fetch_complex_details` to canned data; record calls."""
    state = {"entry": None, "details": {}, "get_urls": [], "detail_ids": []}

    def get_json(url, **kwargs):
        state["get_urls"].append(url)
        return state["entry"]

    def fetch_details(complex_id):
        state["detail_ids"].append(complex_id)
        return state["details"][complex_id]

    monkeypatch.setattr(api, "_get_json", get_json)
    monkeypatch.setattr(api, "fetch_complex_details", fetch_details)
    return state


@pytest.mark.parametrize("raw", ["11gl", " 11GL ", "\t11Gl\n"], ids=["plain", "padded-upper", "mixed"])
def test_pdb_id_is_normalised_and_resolved_to_its_complex(fake_api, raw):
    fake_api["entry"] = _entry_response("11gl", (CPX, True))
    fake_api["details"] = {CPX: _details(CPX)}

    complex_id, details = api.resolve_complex_id(raw, require_dimer=True)

    assert complex_id == CPX
    assert details is fake_api["details"][CPX]
    assert fake_api["get_urls"] == [f"{api.BASE}/complex/details/11gl?id_type=pdb_id"]
    assert fake_api["detail_ids"] == [CPX]


@pytest.mark.parametrize("raw", [CPX, CPX.lower(), f"  {CPX}  "], ids=["plain", "lower", "padded"])
def test_complex_id_is_normalised_and_used_directly(fake_api, raw):
    fake_api["details"] = {CPX: _details(CPX)}

    complex_id, details = api.resolve_complex_id(raw, require_dimer=True)

    assert complex_id == CPX
    assert details is fake_api["details"][CPX]
    assert fake_api["detail_ids"] == [CPX]
    assert fake_api["get_urls"] == []          # no entry lookup


@pytest.mark.parametrize("raw", ["", "  \n"], ids=["empty", "whitespace"])
def test_empty_identifier_is_rejected_before_any_request(fake_api, raw):
    with pytest.raises(ValueError, match="Empty identifier"):
        api.resolve_complex_id(raw)
    assert fake_api["get_urls"] == [] and fake_api["detail_ids"] == []


def test_non_dimer_is_rejected_when_required(fake_api):
    fake_api["entry"] = _entry_response("4hhb", ("PDB-CPX-1", True))
    fake_api["details"] = {"PDB-CPX-1": _details("PDB-CPX-1", total_chains=4,
                                                 oligomeric_state="Heterotetramer")}

    with pytest.raises(ValueError) as excinfo:
        api.resolve_complex_id("4hhb", require_dimer=True)

    message = str(excinfo.value)
    assert "4hhb" in message and "PDB-CPX-1" in message
    assert "4 chains" in message and "Heterotetramer" in message


@pytest.mark.parametrize(
    "payload, match",
    [
        (None, "did not resolve to any complex"),                          # 404
        ({"0xxx": []}, "did not resolve to any complex"),                  # empty list
        ({"0xxx": [{"pdb_complex_id": None}]}, "returned no pdb_complex_id"),
    ],
    ids=["404", "empty-list", "no-complex-id"],
)
def test_pdb_id_with_no_complex_raises(fake_api, payload, match):
    fake_api["entry"] = payload
    with pytest.raises(ValueError, match=match):
        api.resolve_complex_id("0xxx")
    assert fake_api["detail_ids"] == []


def test_multiple_complexes_prefers_the_preferred_assembly_dimer(fake_api):
    fake_api["entry"] = _entry_response("1abc", ("PDB-CPX-A", False), ("PDB-CPX-B", True))
    fake_api["details"] = {
        "PDB-CPX-A": _details("PDB-CPX-A", total_chains=2),
        "PDB-CPX-B": _details("PDB-CPX-B", total_chains=2),
    }
    complex_id, _ = api.resolve_complex_id("1abc")
    assert complex_id == "PDB-CPX-B"


def test_multiple_complexes_falls_back_to_the_only_dimer(fake_api):
    fake_api["entry"] = _entry_response("1abc", ("PDB-CPX-A", False), ("PDB-CPX-B", True))
    fake_api["details"] = {
        "PDB-CPX-A": _details("PDB-CPX-A", total_chains=2),
        "PDB-CPX-B": _details("PDB-CPX-B", total_chains=6),
    }
    complex_id, details = api.resolve_complex_id("1abc")
    assert complex_id == "PDB-CPX-A"
    assert details["total_chains"] == 2


@pytest.mark.parametrize(
    "chains_a, chains_b, preferred",
    [(2, 2, None), (4, 6, "PDB-CPX-A")],
    ids=["two-dimers-no-preference", "no-dimer-at-all"],
)
def test_multiple_complexes_without_a_single_dimer_answer_is_ambiguous(
    fake_api, chains_a, chains_b, preferred,
):
    fake_api["entry"] = _entry_response(
        "1abc", ("PDB-CPX-A", preferred == "PDB-CPX-A"), ("PDB-CPX-B", preferred == "PDB-CPX-B"),
    )
    fake_api["details"] = {
        "PDB-CPX-A": _details("PDB-CPX-A", total_chains=chains_a),
        "PDB-CPX-B": _details("PDB-CPX-B", total_chains=chains_b),
    }
    with pytest.raises(ValueError, match="ambiguous") as excinfo:
        api.resolve_complex_id("1abc")
    assert "PDB-CPX-A" in str(excinfo.value) and "PDB-CPX-B" in str(excinfo.value)
