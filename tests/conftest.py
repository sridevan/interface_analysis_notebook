"""Shared fixtures for the offline unit-test suite.

Every test in `tests/unit` runs without network access: `block_network` replaces
`requests.Session.request` so any accidental call fails loudly. Tests that need to
simulate HTTP responses install their own fake on top of it (see `test_api.py`).

The `make_contact` / `make_interface` fixtures build the smallest payloads that the
`complex/interface_interactions` endpoint would return, so expected values in the
tests can be worked out by inspection.
"""

from __future__ import annotations

import pytest
import requests

# Sentinel meaning "use the author residue number as the UniProt position".
SAME = object()


@pytest.fixture(autouse=True)
def block_network(monkeypatch):
    """Fail any HTTP request made during a unit test."""

    def _blocked(self, method, url, *args, **kwargs):
        raise RuntimeError(
            f"Network access is blocked in unit tests: attempted {method} {url}"
        )

    monkeypatch.setattr(requests.Session, "request", _blocked)


@pytest.fixture
def make_contact():
    """Build one residue-pair record as returned inside `interactions`.

    Defaults describe a heterodimer contact between P00001 (side 1) and
    P00002 (side 2) with the UniProt position equal to the author number.
    Pass `unp1=None` / `unp2=None` / `acc1=None` / `acc2=None` to remove a
    mapping.
    """

    def _make(
        chain1, res1, chain2, res2, *,
        acc1="P00001", acc2="P00002",
        unp1=SAME, unp2=SAME,
        ins1=None, ins2=None,
        bond="hydrogen_bond",
        aa1="K", aa2="D",
    ):
        return {
            "auth_asym_id_1": chain1,
            "auth_seq_id_1": res1,
            "auth_ins_code_1": ins1,
            "auth_asym_id_2": chain2,
            "auth_seq_id_2": res2,
            "auth_ins_code_2": ins2,
            "unp_accession_1": acc1,
            "unp_seq_id_1": res1 if unp1 is SAME else unp1,
            "unp_accession_2": acc2,
            "unp_seq_id_2": res2 if unp2 is SAME else unp2,
            "unp_one_letter_code_1": aa1,
            "unp_one_letter_code_2": aa2,
            "bond_type": bond,
        }

    return _make


@pytest.fixture
def make_interface():
    """Build one interface item as returned by `complex/interface_interactions`."""

    def _make(pdb_id, contacts, *, assembly_id=1, interface_id=1, interface_info=None):
        return {
            "entry_id": pdb_id,
            "assembly_id": assembly_id,
            "interface_id": interface_id,
            "interface_info": interface_info if interface_info is not None
            else {"interface_area": 1000.0},
            "interactions": list(contacts),
        }

    return _make
