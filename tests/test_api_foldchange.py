import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.engine import thermo
from lockr.engine.models import DEFAULT_PARAMS

client = TestClient(app)

# ECLIPSE defaults (K_open=1e-3, K_CK=10nM, lucKey=500nM) at pull=10, the same
# operating point docs/README.md's worked example pins (max_fold_change ~ 7.41,
# regime "mixed", corrected from the pre-_f_open-fix "~11.0/key-limited").
_BASE_REQUEST = {"k_ck": 10.0, "k_open": 0.001, "pull": 10.0, "luckey": 500.0}


def test_foldchange_reproduces_known_eclipse_default_numbers():
    expected_fc = thermo.max_fold_change(100e-12, pull=10, params=DEFAULT_PARAMS)
    expected_regime = thermo.diagnose_regime(DEFAULT_PARAMS, pull=10)

    response = client.post("/foldchange", json={**_BASE_REQUEST, "k_target": None, "target_conc": None})
    assert response.status_code == 200
    body = response.json()

    assert body["fold_change"] == pytest.approx(expected_fc)
    assert body["dominance_ratio"] == pytest.approx(DEFAULT_PARAMS.luckey_ratio) == pytest.approx(50.0)
    assert body["regime"] == "mixed"
    assert body["limiting_factor"] == "mixed"
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


def test_foldchange_warns_when_pull_has_no_documented_improvement():
    response = client.post("/foldchange", json={
        "k_ck": 10.0, "k_open": 0.001, "pull": 75.0, "luckey": 500.0,
        "k_target": None, "target_conc": None,
    })
    assert response.status_code == 200
    assert "no documented improvement" in response.json()["warnings"][0]


def test_foldchange_rejects_non_positive_k_ck():
    response = client.post("/foldchange", json={
        "k_ck": 0.0, "k_open": 0.001, "pull": 10.0, "luckey": 500.0,
        "k_target": None, "target_conc": None,
    })
    assert response.status_code == 400
    assert response.json()["error"]["field"] == "k_ck"


def test_foldchange_rejects_negative_pull():
    response = client.post("/foldchange", json={
        "k_ck": 10.0, "k_open": 0.001, "pull": -1.0, "luckey": 500.0,
        "k_target": None, "target_conc": None,
    })
    assert response.status_code == 400
    assert response.json()["error"]["field"] == "pull"


def test_foldchange_reports_undefined_result_instead_of_crashing_or_nan():
    # k_ck this small underflows to 0.0 once converted nM->M, which would
    # otherwise blow up thermo.py with a ZeroDivisionError (a real 500 we hit
    # while auditing) instead of a clean, spec-worded 400.
    response = client.post("/foldchange", json={
        "k_ck": 1e-310, "k_open": 0.001, "pull": 10.0, "luckey": 500.0,
        "k_target": None, "target_conc": None,
    })
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["code"] == "UNDEFINED_RESULT"
    assert body["error"]["message"] == "Parameters produce an undefined result, check K_open and K_CK."


def test_foldchange_rejects_lone_target_conc_without_k_target():
    response = client.post("/foldchange", json={
        **_BASE_REQUEST, "k_target": None, "target_conc": 5.0,
    })
    assert response.status_code == 400
    message = response.json()["error"]["message"]
    assert not message.startswith("Value error,")  # pydantic's wrapper prefix must not leak
    assert "K_target" in message and "target concentration" in message


def test_foldchange_lod_null_when_target_assumed_saturating():
    response = client.post("/foldchange", json={**_BASE_REQUEST, "k_target": None, "target_conc": None})
    body = response.json()
    assert body["lod_2x_nm"] is None
    assert body["lod_3x_nm"] is None
    assert body["ec50_nm"] is None


def test_foldchange_lod_under_corrected_model():
    # Was checked against Script 7 (0.01008 nM / 0.1008 nM), which used the
    # same additive _f_open formula this test suite corrected, that
    # reference is stale, see tests/test_thermo_eclipse.py's module docstring.
    response = client.post("/foldchange", json={
        **_BASE_REQUEST, "k_target": 0.1, "target_conc": 50.0,  # k_target = 100pM = KD_V10
    })
    body = response.json()
    assert body["lod_2x_nm"] == pytest.approx(0.01283, rel=0.05)
    assert body["ec50_nm"] == pytest.approx(0.0676, rel=0.05)
