import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app

client = TestClient(app)

_BASE_REQUEST = {"k_ck": 10.0, "total_cage_conc": 10.0, "key_stock_conc": 150.0, "assay_volume": 100.0}


def test_dosing_recovers_target_fraction_bound():
    response = client.post("/dosing", json={**_BASE_REQUEST, "target_fraction_bound": 0.5})
    assert response.status_code == 200
    body = response.json()
    assert body["fraction_bound"] == pytest.approx(0.5, rel=1e-3)
    assert body["key_conc_nm"] > 0
    assert body["volume_to_add_ul"] == pytest.approx(body["key_conc_nm"] * 100.0 / 150.0, rel=1e-6)
    assert body["warnings"] == []
    assert len(body["titration"]) > 0


def test_dosing_accepts_fold_change_entry_point():
    response = client.post("/dosing", json={**_BASE_REQUEST, "target_fold_change": 4.0})
    assert response.status_code == 200
    # theta(4, 1) = 4/5 = 0.8
    assert response.json()["fraction_bound"] == pytest.approx(0.8, rel=1e-3)


def test_dosing_rejects_both_targets_set():
    response = client.post("/dosing", json={**_BASE_REQUEST, "target_fraction_bound": 0.5, "target_fold_change": 1.0})
    assert response.status_code == 400
    assert "exactly one" in response.json()["error"]["message"]


def test_dosing_rejects_neither_target_set():
    response = client.post("/dosing", json=_BASE_REQUEST)
    assert response.status_code == 400
    assert "exactly one" in response.json()["error"]["message"]


def test_dosing_rejects_fraction_bound_of_one():
    response = client.post("/dosing", json={**_BASE_REQUEST, "target_fraction_bound": 1.0})
    assert response.status_code == 400


def test_dosing_rejects_non_positive_inputs():
    response = client.post("/dosing", json={**_BASE_REQUEST, "k_ck": 0.0, "target_fraction_bound": 0.5})
    assert response.status_code == 400


def test_dosing_stock_too_dilute_warning():
    # key_stock_conc barely above K_CK/Cage, can't hit a near-saturating target in 100 uL
    response = client.post("/dosing", json={
        **_BASE_REQUEST, "key_stock_conc": 20.0, "target_fraction_bound": 0.99,
    })
    assert response.status_code == 200
    assert any("isn't concentrated enough" in w for w in response.json()["warnings"])


def test_dosing_sub_pipettable_volume_warning():
    response = client.post("/dosing", json={
        **_BASE_REQUEST, "key_stock_conc": 1e9, "target_fraction_bound": 1e-6,
    })
    assert response.status_code == 200
    assert any("pipettable minimum" in w for w in response.json()["warnings"])
