from fastapi.testclient import TestClient

from lockr.api.main import app

client = TestClient(app)

V10 = (
    "SKEAAKKLQDLNIELARKLLEASTKLQRLNIRLAEALLEAIARLQELNLELVYLAVELTDPKRIRDEIKEV"
    "KDKSKEIIRRAEKEIDDAAKESKKILEEARKAIRDAAEESRKILEEGSGSGSDALDELQKLNLELAKLLLKA"
    "IAETQDLNLRAAKAFLEAAAKLQELNIRAVELLVKLTDPATIRRALEHAKRRSKEIIDEAERAIRAAKRESE"
    "RIIEEARRLIEKAKEESERIIREGSGSGDPDIKKLQDLNIELARELLRAHAQLQRLNLELLRELLRALAQLQ"
    "ELNLDLLRLASELTDPDEARKAIAVTGYRLFEEILDAERLISAAALAAIFAAALACRLIREAAAASEKISRE"
)

_LATCH = {"start": 325, "end": 359, "expected_length": 35}
_GRAFT = {"binder": "LISAAALAAIFAAALAC", "start": 327, "spacer": "DA", "spacer_start": 323}
_SMBIT = {"motif": "VTGYRLFEEIL", "start": 312, "end": 322, "label": "SmBiT"}


def test_verify_assembly_reproduces_script6_six_checks_on_v10():
    response = client.post("/verify-assembly", json={
        "full_sequence": V10, "latch_window": _LATCH, "graft_spec": _GRAFT,
        "protected_region": _SMBIT, "expected_total_length": 359,
    })
    assert response.status_code == 200
    body = response.json()
    assert body["all_passed"] is True
    assert len(body["checks"]) == 6
    assert body["variants"] is None


def test_verify_assembly_catches_corrupted_smbit():
    corrupted = V10[:311] + "X" + V10[312:]
    response = client.post("/verify-assembly", json={
        "full_sequence": corrupted, "latch_window": _LATCH, "graft_spec": _GRAFT,
        "protected_region": _SMBIT, "expected_total_length": 359,
    })
    body = response.json()
    assert body["all_passed"] is False
    check = next(c for c in body["checks"] if c["name"] == "protected_region_intact")
    assert check["passed"] is False


def test_verify_assembly_rejects_empty_full_sequence():
    response = client.post("/verify-assembly", json={
        "full_sequence": "", "latch_window": _LATCH, "graft_spec": _GRAFT,
        "protected_region": _SMBIT,
    })
    assert response.status_code == 400


def test_suggest_then_verify_assembly_round_trip_accepts_real_fix():
    # The real coordinate-offset story from docs/README.md step 5: suggest_variant's
    # mutation positions are binder-local (1-17); the binder starts at 327 in V10,
    # so binder_offset=326 must be applied before checking against SmBiT (312-322).
    suggested = client.post("/suggest", json={
        "sequence": "LISDAELEAIFAEELDC",
        "sensitive_window": {"start": 1, "end": 17},
        "substitution_policy": "neutralizing",
        "max_variants": 1,
    }).json()["suggested_variants"]

    response = client.post("/verify-assembly", json={
        "full_sequence": V10, "latch_window": _LATCH, "graft_spec": _GRAFT,
        "protected_region": _SMBIT, "expected_total_length": 359,
        "candidate_variants": suggested, "binder_offset": 326,
    })
    body = response.json()
    assert body["all_passed"] is True
    assert len(body["variants"]["accepted"]) == 1
    assert body["variants"]["rejected"] == []


def test_suggest_then_verify_assembly_round_trip_rejects_variant_inside_protected_region():
    # Same real suggestion, but with an offset that lands its mutations inside
    # SmBiT (312-322) instead of past it -- proves the screen actually rejects,
    # not just that it always accepts.
    suggested = client.post("/suggest", json={
        "sequence": "LISDAELEAIFAEELDC",
        "sensitive_window": {"start": 1, "end": 17},
        "substitution_policy": "neutralizing",
        "max_variants": 1,
    }).json()["suggested_variants"]

    response = client.post("/verify-assembly", json={
        "full_sequence": V10, "latch_window": _LATCH, "graft_spec": _GRAFT,
        "protected_region": _SMBIT, "expected_total_length": 359,
        "candidate_variants": suggested, "binder_offset": 311,
    })
    body = response.json()
    assert body["variants"]["accepted"] == []
    assert len(body["variants"]["rejected"]) == 1
    assert "falls inside protected region" in body["variants"]["rejected"][0]["reason"]
