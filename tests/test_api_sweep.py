import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.engine import thermo
from lockr.engine.models import DEFAULT_PARAMS

client = TestClient(app)

_BASE = {"k_ck": 10.0, "k_open_off": 0.001, "k_open_on": 0.011, "luckey": 500.0}


def test_sweep_operating_point_matches_known_eclipse_fold_change():
    expected_fc = thermo.max_fold_change(100e-12, pull=10, params=DEFAULT_PARAMS)

    response = client.post("/sweep", json={
        "base_params": _BASE,
        "sweep": {"param": "luckey", "min": 1, "max": 10000, "steps": 5, "scale": "log"},
    })
    assert response.status_code == 200
    body = response.json()

    assert body["operating_point"]["x"] == 500.0
    assert body["operating_point"]["fold_change"] == pytest.approx(expected_fc)
    assert len(body["points"]) == 5
    # key-limited at this operating point: fold-change barely moves across decades.
    fcs = [p["fold_change"] for p in body["points"]]
    assert max(fcs) - min(fcs) < 1.0


def test_sweep_dominance_ratio_tracks_swept_param():
    response = client.post("/sweep", json={
        "base_params": _BASE,
        "sweep": {"param": "luckey", "min": 1, "max": 100, "steps": 3, "scale": "log"},
    })
    points = response.json()["points"]
    assert points[0]["dominance_ratio"] == pytest.approx(0.1)
    assert points[-1]["dominance_ratio"] == pytest.approx(10.0)


def test_sweep_rejects_min_greater_than_max():
    response = client.post("/sweep", json={
        "base_params": _BASE,
        "sweep": {"param": "luckey", "min": 10000, "max": 1, "steps": 5, "scale": "log"},
    })
    assert response.status_code == 400


def test_sweep_rejects_unknown_param():
    response = client.post("/sweep", json={
        "base_params": _BASE,
        "sweep": {"param": "pull", "min": 1, "max": 100, "steps": 5, "scale": "log"},
    })
    assert response.status_code == 400
