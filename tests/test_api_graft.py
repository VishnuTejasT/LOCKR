import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.api.routes import graft as graft_route
from lockr.engine import graft as graft_engine

client = TestClient(app)

V10_OPTIMIZED = "LISAAALAAIFAAALAC"


def test_graft_status_reports_available_when_pyrosetta_installed():
    response = client.get("/graft/status")
    assert response.status_code == 200
    body = response.json()
    assert body["available"] == graft_engine.PYROSETTA_AVAILABLE
    assert body["template_bundled"] is True


def test_graft_returns_503_when_pyrosetta_unavailable(monkeypatch):
    monkeypatch.setattr(graft_engine, "PYROSETTA_AVAILABLE", False)
    response = client.post("/graft", json={"sequence": V10_OPTIMIZED, "scan_all": True})
    assert response.status_code == 503
    body = response.json()
    assert body["error"]["code"] == "PYROSETTA_UNAVAILABLE"
    assert "install_guide" in body["error"]


def test_graft_rejects_sequence_too_long():
    response = client.post("/graft", json={"sequence": "A" * 40, "scan_all": True})
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["field"] == "sequence"
    assert "too long" in body["error"]["message"]


def test_graft_rejects_invalid_residue():
    response = client.post("/graft", json={"sequence": "LISDX1EL", "scan_all": True})
    assert response.status_code == 400
    assert body_error_code(response) == "VALIDATION_ERROR"


def test_graft_rejects_specific_position_out_of_range():
    response = client.post("/graft", json={
        "sequence": V10_OPTIMIZED, "scan_all": False, "specific_position": 999,
    })
    assert response.status_code == 400


def body_error_code(response):
    return response.json()["error"]["code"]


@pytest.mark.slow
@pytest.mark.skipif(not graft_engine.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
def test_graft_full_run_matches_eclipse_v10():
    response = client.post("/graft", json={"sequence": V10_OPTIMIZED, "scan_all": True})
    assert response.status_code == 200
    body = response.json()
    assert body["best_position"] == 327
    assert body["best_score"] == pytest.approx(-1469.57, abs=5)
    assert body["verdict"] == "good"
    assert len(body["grafted_sequence"]) == 359
    assert body["grafted_sequence"][326:343] == V10_OPTIMIZED
    assert "job_id" in body

    download = client.get(f"/graft/download/{body['job_id']}")
    assert download.status_code == 200
    assert download.content.startswith(b"HEADER") or b"ATOM" in download.content[:200]


def test_graft_download_404_for_unknown_job():
    response = client.get("/graft/download/00000000-0000-0000-0000-000000000000")
    assert response.status_code == 404
