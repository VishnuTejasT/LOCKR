"""API layer for Module C, plus the regression that Module B is untouched at 298.15 K."""

import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.engine.temperature import T_REF_DEFAULT, dG_from_k

client = TestClient(app)

_DG_LATCH = dG_from_k(1e-3, T_REF_DEFAULT)      # K_open = 1e-3
_DG_KEY = dG_from_k(1e-8, T_REF_DEFAULT)        # K_CK = 10 nM
_BOUNDS = {
    "dG_latch_cage": {"min": -10.0, "max": -1.0},
    "dG_key_cage": {"min": -15.0, "max": -6.0},
    "key_conc": {"min": 1.0, "max": 5000.0},     # nM at the boundary
}


def _fixed(**over):
    base = {"pull": 20.0, "dG_latch_cage": _DG_LATCH, "dG_key_cage": _DG_KEY, "key_conc": 500.0}
    base.update(over)
    return base


def _post(payload, path="/api/inverse_design"):
    return client.post(path, json=payload)


# --- inverse design ---------------------------------------------------------

def test_one_free_parameter_returns_a_single_solution():
    fixed = _fixed()
    fixed.pop("dG_latch_cage")
    response = _post({"solve_for": ["dG_latch_cage"], "fixed": fixed,
                      "bounds": {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                      "target_fold_change": 8.0})
    assert response.status_code == 200
    body = response.json()
    assert body["feasible"] is True
    assert body["mode"] == "root"
    assert len(body["solutions"]) == 1
    solution = body["solutions"][0]
    assert solution["fold_change"] == pytest.approx(8.0, rel=1e-6)
    assert -10.0 <= solution["values"]["dG_latch_cage"] <= -1.0
    assert "dG_latch_cage" in solution["sensitivity"]
    assert solution["tolerance"]["dG_latch_cage"]["min"] < solution["values"]["dG_latch_cage"]
    assert body["temperature_K"] == T_REF_DEFAULT
    assert body["thermo_mode"] == "isothermal"
    assert body["units"]["key_conc"] == "nM"


def test_concentrations_cross_the_boundary_in_nanomolar():
    fixed = _fixed()
    fixed.pop("key_conc")
    response = _post({"solve_for": ["key_conc"], "fixed": fixed,
                      "bounds": {"key_conc": _BOUNDS["key_conc"]}, "target_fold_change": 8.0})
    body = response.json()
    value = body["solutions"][0]["values"]["key_conc"]
    assert 1.0 <= value <= 5000.0, "a key concentration in nM, not M"
    tolerance = body["solutions"][0]["tolerance"]["key_conc"]
    assert tolerance["min"] < value < tolerance["max"]
    assert 1.0 <= tolerance["min"] and tolerance["max"] <= 5000.0


def test_two_free_parameters_return_region_polylines():
    fixed = _fixed()
    for name in ("dG_latch_cage", "key_conc"):
        fixed.pop(name)
    response = _post({"solve_for": ["dG_latch_cage", "key_conc"], "fixed": fixed,
                      "bounds": {k: _BOUNDS[k] for k in ("dG_latch_cage", "key_conc")},
                      "target_fold_change": 8.0})
    assert response.status_code == 200
    body = response.json()
    assert body["mode"] == "region"
    assert body["region_axes"] == ["dG_latch_cage", "key_conc"]
    assert body["region"] and all(len(point) == 2 for line in body["region"] for point in line)
    for line in body["region"]:
        for dG, key_nm in line:
            assert -10.0 <= dG <= -1.0
            assert 1.0 <= key_nm <= 5000.0
    assert body["robust_point"]["on_spec"] is True


def test_three_free_parameters_return_a_robust_point_and_no_region():
    names = ["dG_latch_cage", "dG_key_cage", "key_conc"]
    response = _post({"solve_for": names, "fixed": {"pull": 20.0},
                      "bounds": {k: _BOUNDS[k] for k in names}, "target_fold_change": 8.0})
    assert response.status_code == 200
    body = response.json()
    assert body["mode"] == "robust_point"
    assert body["region"] == []
    assert body["region_axes"] is None
    assert set(body["robust_point"]["tolerance"]) <= set(names)
    assert set(body["robust_point"]["sensitivity"]) == set(names)


def test_target_range_is_accepted():
    fixed = _fixed()
    fixed.pop("dG_latch_cage")
    response = _post({"solve_for": ["dG_latch_cage"], "fixed": fixed,
                      "bounds": {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                      "target_fold_change_range": {"min": 7.0, "max": 9.0}})
    assert response.status_code == 200
    assert 7.0 <= response.json()["solutions"][0]["fold_change"] <= 9.0


def test_unreachable_target_is_a_200_with_the_infeasible_path():
    fixed = _fixed(pull=10.0)
    fixed.pop("dG_latch_cage")
    response = _post({"solve_for": ["dG_latch_cage"], "fixed": fixed,
                      "bounds": {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                      "target_fold_change": 1000.0})
    assert response.status_code == 200
    body = response.json()
    assert body["feasible"] is False
    assert body["solutions"] == []
    assert body["closest"]["fold_change"] < 12.0
    assert "dG_latch_cage" in body["closest"]["binding_constraint"]
    assert "not reachable" in body["message"]


@pytest.mark.parametrize("payload,field", [
    ({"solve_for": [], "fixed": {"pull": 10.0}, "bounds": {}, "target_fold_change": 8.0},
     "solve_for"),
    ({"solve_for": ["nonsense"], "fixed": {"pull": 10.0}, "bounds": {}, "target_fold_change": 8.0},
     "solve_for"),
    ({"solve_for": ["key_conc"], "fixed": {}, "bounds": {"key_conc": {"min": 1.0, "max": 10.0}},
      "target_fold_change": 8.0}, "fixed"),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0, "nonsense": 1.0},
      "bounds": {"key_conc": {"min": 1.0, "max": 10.0}}, "target_fold_change": 8.0}, "fixed"),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0}, "bounds": {},
      "target_fold_change": 8.0}, None),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0},
      "bounds": {"key_conc": {"min": 10.0, "max": 1.0}}, "target_fold_change": 8.0}, None),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0},
      "bounds": {"key_conc": {"min": 0.0, "max": 10.0}}, "target_fold_change": 8.0}, "bounds"),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0},
      "bounds": {"key_conc": {"min": 1.0, "max": 10.0}}}, None),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0},
      "bounds": {"key_conc": {"min": 1.0, "max": 10.0}}, "target_fold_change": 8.0,
      "target_fold_change_range": {"min": 7.0, "max": 9.0}}, None),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0},
      "bounds": {"key_conc": {"min": 1.0, "max": 10.0}}, "target_fold_change": 8.0,
      "temperature_K": 700.0}, "temperature_K"),
    ({"solve_for": ["key_conc"], "fixed": {"pull": 10.0},
      "bounds": {"key_conc": {"min": 1.0, "max": 10.0}}, "target_fold_change": 8.0,
      "thermo": {"dG_nonsense": {"dH_ref": -10.0}}}, "thermo"),
])
def test_bad_requests_return_the_error_envelope(payload, field):
    response = _post(payload)
    assert response.status_code == 400
    error = response.json()["error"]
    assert error["code"] == "VALIDATION_ERROR"
    assert error["message"]
    if field:
        assert field in (error["field"] or "") + error["message"]


def test_missing_fixed_parameter_comes_back_as_a_validation_error():
    response = _post({"solve_for": ["dG_latch_cage"], "fixed": {"pull": 10.0, "key_conc": 500.0},
                      "bounds": {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                      "target_fold_change": 8.0})
    assert response.status_code == 400
    assert "dG_key_cage" in response.json()["error"]["message"]


# --- temperature on the Module C endpoint -----------------------------------

def test_isothermal_mode_refuses_to_extrapolate_over_the_api():
    fixed = _fixed()
    fixed.pop("dG_latch_cage")
    response = _post({"solve_for": ["dG_latch_cage"], "fixed": fixed,
                      "bounds": {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                      "target_fold_change": 8.0, "temperature_K": 310.0})
    body = response.json()
    assert body["thermo_mode"] == "isothermal"
    assert any("TEMPERATURE SCALING UNAVAILABLE" in w for w in body["warnings"])


def test_full_thermo_block_labels_itself_full_and_reports_dCp_estimates():
    fixed = _fixed()
    fixed.pop("dG_latch_cage")
    response = _post({
        "solve_for": ["dG_latch_cage"], "fixed": fixed,
        "bounds": {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
        "target_fold_change": 8.0, "temperature_K": 310.0,
        "thermo": {"dG_latch_cage": {"dH_ref": -20.0, "buried_nonpolar_A2": 900.0,
                                     "buried_polar_A2": 400.0, "Tm": 315.0},
                   "dG_key_cage": {"dH_ref": -15.0, "dCp": -0.4}},
    })
    assert response.status_code == 200
    body = response.json()
    assert body["thermo_mode"] == "full"
    assert any("estimated from buried surface area" in w for w in body["warnings"])
    assert any("two-state" in w and "unfolding" in w for w in body["warnings"])


# --- Module B regression ----------------------------------------------------

_MODULE_B_REQUEST = {"k_ck": 10.0, "k_open": 0.001, "pull": 10.0, "luckey": 500.0}


def test_default_temperature_leaves_foldchange_byte_for_byte_unchanged():
    without = client.post("/foldchange", json=_MODULE_B_REQUEST).json()
    with_default = client.post("/foldchange", json={**_MODULE_B_REQUEST,
                                                    "temperature_K": 298.15}).json()
    assert with_default == without
    assert without["temperature_K"] == T_REF_DEFAULT
    assert without["thermo_mode"] == "isothermal"
    assert without["warnings"] == []


def test_default_temperature_leaves_the_targeted_foldchange_unchanged():
    payload = {**_MODULE_B_REQUEST, "k_target": 1.0, "target_conc": 1.0}
    assert client.post("/foldchange", json={**payload, "temperature_K": 298.15}).json() == \
        client.post("/foldchange", json=payload).json()


def test_foldchange_at_another_temperature_warns_but_does_not_silently_scale():
    body = client.post("/foldchange", json={**_MODULE_B_REQUEST, "temperature_K": 310.0}).json()
    assert body["temperature_K"] == 310.0
    assert body["thermo_mode"] == "isothermal"
    assert any("TEMPERATURE SCALING UNAVAILABLE" in w for w in body["warnings"])
    # No dH means nothing was extrapolated, so the number matches the reference one.
    assert body["fold_change"] == client.post("/foldchange", json=_MODULE_B_REQUEST).json()["fold_change"]


def test_foldchange_with_enthalpy_actually_moves_with_temperature():
    thermo = {"dG_latch_cage": {"dH_ref": -20.0, "dCp": -0.5},
              "dG_key_cage": {"dH_ref": -15.0, "dCp": -0.4}}
    warm = client.post("/foldchange", json={**_MODULE_B_REQUEST, "temperature_K": 310.0,
                                            "thermo": thermo}).json()
    cold = client.post("/foldchange", json={**_MODULE_B_REQUEST, "temperature_K": 277.0,
                                            "thermo": thermo}).json()
    assert warm["thermo_mode"] == "full"
    assert warm["fold_change"] != cold["fold_change"]
    assert warm["warnings"] == [] or all("Tm" not in w for w in warm["warnings"])


def test_sweep_and_kopen_accept_temperature_without_changing_their_defaults():
    sweep_payload = {"base_params": _MODULE_B_REQUEST,
                     "sweep": {"param": "k_open", "min": 1e-5, "max": 1e-1, "steps": 5}}
    assert client.post("/sweep", json={**sweep_payload, "temperature_K": 298.15}).json() == \
        client.post("/sweep", json=sweep_payload).json()

    kopen_payload = {"k_open_current": 0.001, "k_ck": 10.0, "luckey": 500.0, "pull": 10.0}
    assert client.post("/kopen-scenarios", json={**kopen_payload, "temperature_K": 298.15}).json() == \
        client.post("/kopen-scenarios", json=kopen_payload).json()

    for path, payload in (("/sweep", sweep_payload), ("/kopen-scenarios", kopen_payload)):
        body = client.post(path, json={**payload, "temperature_K": 310.0}).json()
        assert body["temperature_K"] == 310.0
        assert body["thermo_mode"] == "isothermal"
        assert any("TEMPERATURE SCALING UNAVAILABLE" in w for w in body["warnings"])


# --- ddG to mutations -------------------------------------------------------

def test_mutation_mapping_returns_ranked_candidates_with_a_disclaimer():
    response = client.post("/api/mutations", json={
        "latch_sequence": "LISDAELEAIFAEELDC", "current_dG_latch_cage": -4.09,
        "target_dG_latch_cage": -4.90, "interface_positions": [4, 6, 8, 13, 14, 16], "top_n": 5})
    assert response.status_code == 200
    body = response.json()
    assert body["required_ddG"] == pytest.approx(-0.81)
    assert len(body["candidates"]) == 5
    assert "HEURISTIC" in body["disclaimer"]
    assert body["disclaimer"] in body["warnings"]
    first = body["candidates"][0]
    assert first["confidence"] in ("moderate", "low")
    assert first["notes"]
    assert first["resulting_dG"] == pytest.approx(-4.09 + first["ddG_estimate"])


@pytest.mark.parametrize("payload", [
    {"latch_sequence": "LISDX", "current_dG_latch_cage": -4.0, "target_dG_latch_cage": -5.0,
     "interface_positions": [1]},
    {"latch_sequence": "LISDAEL", "current_dG_latch_cage": -4.0, "target_dG_latch_cage": -5.0,
     "interface_positions": []},
    {"latch_sequence": "LISDAEL", "current_dG_latch_cage": -4.0, "target_dG_latch_cage": -5.0,
     "interface_positions": [99]},
    {"latch_sequence": "LISDAEL", "current_dG_latch_cage": -4.0, "target_dG_latch_cage": -5.0,
     "interface_positions": [1], "top_n": 0},
])
def test_bad_mutation_requests_return_the_error_envelope(payload):
    response = client.post("/api/mutations", json=payload)
    assert response.status_code == 400
    assert response.json()["error"]["code"] == "VALIDATION_ERROR"
