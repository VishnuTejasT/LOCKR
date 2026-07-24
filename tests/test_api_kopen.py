import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app

client = TestClient(app)

# ECLIPSE defaults: K_open=1e-3, K_CK=10nM, lucKey=500nM, pull=10.
_BASE_REQUEST = {"k_open_current": 0.001, "k_ck": 10.0, "luckey": 500.0, "pull": 10.0}


def test_kopen_scenarios_crossover_matches_v11_finding():
    response = client.post("/kopen-scenarios", json=_BASE_REQUEST)
    assert response.status_code == 200
    body = response.json()
    assert body["crossover_luckey_nm"] == pytest.approx(0.01)
    assert body["dominance_ratio"] == pytest.approx(50.0)


def test_kopen_scenarios_moderate_preset_never_helps_at_500nm_luckey():
    response = client.post("/kopen-scenarios", json=_BASE_REQUEST)
    body = response.json()
    moderate = [s for s in body["scenarios"] if s["preset"] == "moderate"]
    assert len(moderate) == 3
    assert all(s["helps"] != "yes" for s in moderate)
    for s in moderate:
        assert abs(s["fold_change"] - body["baseline"]["fold_change"]) / body["baseline"]["fold_change"] < 0.1


def test_kopen_scenarios_baseline_matches_max_fold_change():
    response = client.post("/kopen-scenarios", json=_BASE_REQUEST)
    body = response.json()
    assert body["baseline"]["mutations"] == 0
    assert body["baseline"]["helps"] == "unchanged"
    assert body["baseline"]["fold_change"] == pytest.approx(11.0, rel=2e-3)


def test_kopen_scenarios_never_help_below_crossover():
    # max_fold_change(K_open) is monotonically decreasing in K_open (loosening
    # the cage only ever raises background), so destabilizing mutations can
    # only leave saturating fold-change flat or make it worse -- never "yes".
    # This holds even well below the lucKey/K_CK crossover.
    response = client.post("/kopen-scenarios", json={
        "k_open_current": 0.001, "k_ck": 10.0, "luckey": 0.01, "pull": 10.0,
    })
    body = response.json()
    assert all(s["helps"] != "yes" for s in body["scenarios"])
    optimistic = [s for s in body["scenarios"] if s["preset"] == "optimistic"]
    assert any(s["helps"] == "slightly_worse" for s in optimistic)


def test_kopen_scenarios_rejects_non_positive_k_open():
    response = client.post("/kopen-scenarios", json={**_BASE_REQUEST, "k_open_current": 0.0})
    assert response.status_code == 400
    assert response.json()["error"]["field"] == "k_open_current"
