from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.engine import calibration, liability

client = TestClient(app)

ORIGINAL = "LISDAELEAIFAEELDC"


def _scan_request(sequence=ORIGINAL, **overrides):
    body = {
        "sequences": [{"id": "binder_v1", "sequence": sequence}],
        "sensitive_window": {"start": 1, "end": len(sequence)},
        "ph": 7.4,
        "substitution_policy": "conservative",
    }
    body.update(overrides)
    return body


def test_scan_reproduces_known_eclipse_liability_numbers():
    expected = liability.scan_liability(ORIGINAL, preserve_positions=[])

    response = client.post("/scan", json=_scan_request())
    assert response.status_code == 200
    result = response.json()["results"][0]

    assert result["id"] == "binder_v1"
    assert len(result["acidic_residues"]) == len(expected.liabilities) == 6
    assert result["liability_score"] == expected.liability_score
    assert result["liability_band"] == expected.liability_band == "high"


def test_scan_suggested_variant_matches_engine_neutralizing_fix():
    response = client.post("/scan", json=_scan_request(substitution_policy="neutralizing"))
    variant = response.json()["results"][0]["suggested_variants"][0]

    expected = liability.suggest_variant(ORIGINAL, preserve_positions=[], policy="neutralizing")
    assert variant["sequence"] == expected.sequence
    assert variant["liability_band"] == "low"
    assert variant["estimated_kck_nm"] == expected.K_CK_estimate * 1e9


def test_scan_rejects_non_standard_residue():
    response = client.post("/scan", json=_scan_request(sequence="LISDAELXAIFAEELDC"))
    assert response.status_code == 400
    body = response.json()
    assert "non-standard residue 'X'" in body["error"]["message"]


def test_scan_rejects_empty_sequence_list():
    response = client.post("/scan", json={
        "sequences": [],
        "sensitive_window": {"start": 1, "end": 17},
    })
    assert response.status_code == 400


def test_suggest_round_trips_into_scan_shape():
    response = client.post("/suggest", json={
        "sequence": ORIGINAL,
        "sensitive_window": {"start": 1, "end": 17},
        "substitution_policy": "neutralizing",
        "max_variants": 1,
    })
    assert response.status_code == 200
    variants = response.json()["suggested_variants"]
    assert len(variants) == 1
    assert variants[0]["substitutions"][0] == {"position": 4, "from": "D", "to": "A"}
