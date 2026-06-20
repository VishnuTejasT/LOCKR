import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.engine import thermo
from lockr.engine.models import DEFAULT_PARAMS

client = TestClient(app)

# ECLIPSE defaults (K_open=1e-3, K_CK=10nM, lucKey=500nM) at pull=10, the same
# operating point docs/README.md's worked example pins (max_fold_change ~ 11.0,
# regime "key-limited").
_BASE_REQUEST = {"k_ck": 10.0, "k_open_off": 0.001, "k_open_on": 0.011, "luckey": 500.0}


def test_foldchange_reproduces_known_eclipse_default_numbers():
    expected_fc = thermo.max_fold_change(100e-12, pull=10, params=DEFAULT_PARAMS)
    expected_regime = thermo.diagnose_regime(DEFAULT_PARAMS, pull=10)

    response = client.post("/foldchange", json={**_BASE_REQUEST, "k_target": None, "target_conc": None})
    assert response.status_code == 200
    body = response.json()

    assert body["fold_change"] == pytest.approx(expected_fc)
    assert body["dominance_ratio"] == pytest.approx(DEFAULT_PARAMS.luckey_ratio) == pytest.approx(50.0)
    assert body["regime"] == "key_limited"
    assert body["limiting_factor"] == "luckey_over_kck"
    assert body["verdict"] == expected_regime.verdict
    assert body["warnings"] == []


def test_foldchange_with_partial_target_occupancy_uses_theta():
    # k_target/target_conc both set -> goes through thermo.fold_change with a
    # real theta < 1, so fc should be lower than the saturating ceiling above.
    response = client.post("/foldchange", json={
        **_BASE_REQUEST, "k_target": 1.0, "target_conc": 1.0,
    })
    saturating = client.post("/foldchange", json={**_BASE_REQUEST, "k_target": None, "target_conc": None})
    assert response.json()["fold_change"] < saturating.json()["fold_change"]


def test_foldchange_warns_when_on_state_does_not_stabilize_open():
    response = client.post("/foldchange", json={
        "k_ck": 10.0, "k_open_off": 0.5, "k_open_on": 0.1, "luckey": 500.0,
        "k_target": None, "target_conc": None,
    })
    assert response.status_code == 200
    assert "ON > OFF" in response.json()["warnings"][0]


def test_foldchange_rejects_non_positive_k_ck():
    response = client.post("/foldchange", json={
        "k_ck": 0.0, "k_open_off": 0.001, "k_open_on": 0.011, "luckey": 500.0,
        "k_target": None, "target_conc": None,
    })
    assert response.status_code == 400
    assert response.json()["error"]["field"] == "k_ck"


def test_foldchange_rejects_lone_target_conc_without_k_target():
    response = client.post("/foldchange", json={
        **_BASE_REQUEST, "k_target": None, "target_conc": 5.0,
    })
    assert response.status_code == 400
