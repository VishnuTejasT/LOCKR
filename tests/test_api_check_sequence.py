from fastapi.testclient import TestClient

from lockr.api.main import app

client = TestClient(app)

# Same untagged v1.0 sequence test_api_assembly.py uses (359aa, SmBiT at 312-322).
V10 = (
    "SKEAAKKLQDLNIELARKLLEASTKLQRLNIRLAEALLEAIARLQELNLELVYLAVELTDPKRIRDEIKEV"
    "KDKSKEIIRRAEKEIDDAAKESKKILEEARKAIRDAAEESRKILEEGSGSGSDALDELQKLNLELAKLLLKA"
    "IAETQDLNLRAAKAFLEAAAKLQELNIRAVELLVKLTDPATIRRALEHAKRRSKEIIDEAERAIRAAKRESE"
    "RIIEEARRLIEKAKEESERIIREGSGSGDPDIKKLQDLNIELARELLRAHAQLQRLNLELLRELLRALAQLQ"
    "ELNLDLLRLASELTDPDEARKAIAVTGYRLFEEILDAERLISAAALAAIFAAALACRLIREAAAASEKISRE"
)

# His-TEV-tagged v1.0 (379aa, SmBiT shifts to 332-342). The tag's trailing
# residue is the core's own leading S, so it's 20aa here, not 21.
V10_TAGGED = "MGSHHHHHHGSGSENLYFQG" + V10

BINDER_ONLY = "LISDAELEAIFAEELDC"


def test_check_sequence_clean_379aa_has_no_warnings():
    response = client.post("/check-sequence", json={"sequence": V10_TAGGED})
    assert response.status_code == 200
    body = response.json()
    assert body["length"] == 379
    assert body["smbit_found"] is True
    assert body["smbit_start"] == 332
    assert body["smbit_end"] == 342
    assert body["warnings"] == []


def test_check_sequence_clean_359aa_has_no_warnings():
    response = client.post("/check-sequence", json={"sequence": V10})
    body = response.json()
    assert body["length"] == 359
    assert body["smbit_found"] is True
    assert body["smbit_start"] == 312
    assert body["warnings"] == []


def test_check_sequence_corrupted_smbit_warns():
    corrupted = V10_TAGGED[:331] + "X" + V10_TAGGED[332:]
    response = client.post("/check-sequence", json={"sequence": corrupted})
    body = response.json()
    assert body["smbit_found"] is False
    assert len(body["warnings"]) == 1
    assert "not found at expected position" in body["warnings"][0]


def test_check_sequence_binder_only_is_silent():
    response = client.post("/check-sequence", json={"sequence": BINDER_ONLY})
    body = response.json()
    assert body["length"] == 17
    assert body["smbit_found"] is False
    assert body["warnings"] == []


def test_check_sequence_unrecognized_length_with_smbit_flags_length():
    padded = V10_TAGGED + "AAA"  # 382aa, not a known variant, still has SmBiT
    response = client.post("/check-sequence", json={"sequence": padded})
    body = response.json()
    assert body["length"] == 382
    assert body["smbit_found"] is True
    assert any("doesn't match known lucCage variants" in w for w in body["warnings"])


def test_check_sequence_rejects_empty_sequence():
    response = client.post("/check-sequence", json={"sequence": ""})
    assert response.status_code == 400
