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
    assert result["liability_band"] == expected.liability_band == "High"


def test_scan_suggested_variant_matches_engine_neutralizing_fix():
    response = client.post("/scan", json=_scan_request(substitution_policy="neutralizing"))
    variant = response.json()["results"][0]["suggested_variants"][0]

    expected = liability.suggest_variant(ORIGINAL, preserve_positions=[], policy="neutralizing")
    assert variant["sequence"] == expected.sequence
    assert variant["liability_band"] == "Low"
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


def test_scan_preserve_positions_excludes_protected_residue_from_liabilities():
    # D4 is a flagged liability with no protection; preserving it should drop
    # it from acidic_residues/per_position entirely, same as the engine itself.
    response = client.post("/scan", json=_scan_request(preserve_positions=[4]))
    result = response.json()["results"][0]

    flagged_positions = {r["position"] for r in result["acidic_residues"]}
    assert 4 not in flagged_positions
    assert len(result["acidic_residues"]) == 5  # 6 total liabilities minus the preserved one


def test_suggest_preserve_positions_excludes_protected_residue_from_mutations():
    response = client.post("/suggest", json={
        "sequence": ORIGINAL,
        "sensitive_window": {"start": 1, "end": 17},
        "substitution_policy": "neutralizing",
        "max_variants": 1,
        "preserve_positions": [4],
    })
    variant = response.json()["suggested_variants"][0]

    mutated_positions = {s["position"] for s in variant["substitutions"]}
    assert 4 not in mutated_positions
    assert variant["sequence"][3] == "D"  # position 4 (1-indexed) untouched


def test_scan_warns_on_long_sequence():
    long_seq = "A" * 250
    response = client.post("/scan", json={
        "sequences": [{"id": "a", "sequence": long_seq}],
        "sensitive_window": {"start": 1, "end": len(long_seq)},
    })
    assert response.status_code == 200
    result = response.json()["results"][0]
    assert result["warnings"] == ["long sequence: liability model tuned for peptide-scale binders"]


def test_scan_has_no_warnings_for_typical_length_sequence():
    # preserve_positions is set, so the "you didn't protect anything" nudge
    # doesn't fire here, this test is only about the long-sequence warning.
    response = client.post("/scan", json=_scan_request(preserve_positions=calibration.PFLDH_INTERFACE))
    result = response.json()["results"][0]
    assert result["warnings"] == []


def test_scan_warns_when_mutations_suggested_without_preserve_positions():
    response = client.post("/scan", json=_scan_request())
    result = response.json()["results"][0]
    assert any("preserve_positions were set" in w for w in result["warnings"])


def test_scan_no_preserve_warning_when_no_mutations_suggested():
    response = client.post("/scan", json=_scan_request(sequence="LISAAALAAIFAAALAC"))
    result = response.json()["results"][0]
    assert result["warnings"] == []


def test_scan_low_band_note_is_honest_about_acidic_residues_present():
    # A single D scores Low but isn't zero acidic residues, the note text
    # must not claim there are none when acidic_residues is non-empty.
    response = client.post("/scan", json={
        "sequences": [{"id": "x", "sequence": "D"}],
        "sensitive_window": {"start": 1, "end": 1},
    })
    result = response.json()["results"][0]
    assert result["liability_band"] == "Low"
    assert len(result["acidic_residues"]) == 1
    assert "not any acidic residues" not in result["predicted_kck_penalty"]["note"]


def test_scan_low_band_note_says_none_when_truly_none():
    response = client.post("/scan", json=_scan_request(sequence="LISAAALAAIFAAALAC"))
    result = response.json()["results"][0]
    assert result["liability_band"] == "Low"
    assert result["acidic_residues"] == []
    assert "not any acidic residues" in result["predicted_kck_penalty"]["note"]


def test_scan_rejects_preserve_position_out_of_range():
    response = client.post("/scan", json=_scan_request(preserve_positions=[9999]))
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["field"] == "preserve_positions"
    assert "9999" in body["error"]["message"]


def test_scan_rejects_window_end_beyond_sequence_length():
    response = client.post("/scan", json=_scan_request(sensitive_window={"start": 1, "end": 999}))
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["field"] == "sensitive_window.end"


def test_scan_rejects_window_start_below_one():
    response = client.post("/scan", json=_scan_request(sensitive_window={"start": 0, "end": 5}))
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["field"] == "sensitive_window.start"


def test_suggest_rejects_preserve_position_out_of_range():
    response = client.post("/suggest", json={
        "sequence": ORIGINAL,
        "sensitive_window": {"start": 1, "end": len(ORIGINAL)},
        "substitution_policy": "neutralizing",
        "max_variants": 1,
        "preserve_positions": [100],
    })
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["field"] == "preserve_positions"


def test_suggest_rejects_window_out_of_range():
    response = client.post("/suggest", json={
        "sequence": ORIGINAL,
        "sensitive_window": {"start": 1, "end": 500},
        "substitution_policy": "neutralizing",
        "max_variants": 1,
    })
    assert response.status_code == 400
    body = response.json()
    assert body["error"]["field"] == "sensitive_window.end"


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
