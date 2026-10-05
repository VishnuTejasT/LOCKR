"""Module D API: units, refusals, caps, slice."""

import math

import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.engine.landscape import MAX_CELLS, MAX_CELLS_LOD
from lockr.engine.temperature import R_KCAL, T_REF_DEFAULT, dG_from_k, k_from_dG

client = TestClient(app)

_FIXED = {"k_ck": 10.0, "k_open": 1e-3, "luckey": 500.0, "pull": 10.0}
_X = {"param": "k_open", "min": 1e-5, "max": 1e-1, "steps": 40, "scale": "log"}
_Y = {"param": "pull", "min": 1.0, "max": 30.0, "steps": 30, "scale": "linear"}
_THERMO = {"dG_latch_cage": {"dH_ref": -12.0, "dCp": 0.4, "Tm": 335.0},
           "dG_key_cage": {"dH_ref": -8.0}}


def _body(**over):
    base = {"x_axis": dict(_X), "y_axis": dict(_Y), "fixed": dict(_FIXED)}
    base.update(over)
    return base


def _post(payload, path="/api/landscape"):
    return client.post(path, json=payload)


def _ok(payload, path="/api/landscape"):
    response = _post(payload, path)
    assert response.status_code == 200, response.json()
    return response.json()


# --- shape and units -------------------------------------------------------


def test_grid_shape_axes_and_units():
    data = _ok(_body())
    assert len(data["z"]) == 40 and len(data["z"][0]) == 30
    assert data["x"]["param"] == "k_open" and data["x"]["unit"] == "dimensionless"
    assert data["y"]["values"][0] == pytest.approx(1.0)
    assert data["cells"] == 1200
    assert data["metric"] == "fold_change"
    assert data["temperature_K"] == pytest.approx(T_REF_DEFAULT)
    assert data["thermo_mode"] == "isothermal"


def test_concentration_axis_reported_in_nm():
    data = _ok(_body(x_axis={"param": "k_ck", "min": 1.0, "max": 100.0, "steps": 12,
                             "scale": "log"}))
    assert data["x"]["unit"] == "nM"
    assert data["x"]["values"][0] == pytest.approx(1.0)
    assert data["x"]["values"][-1] == pytest.approx(100.0)


def test_lod_metric_is_returned_in_nm():
    data = _ok(_body(metric="limit_of_detection",
                     fixed=dict(_FIXED, k_target=1.0),
                     x_axis={"param": "k_open", "min": 1e-4, "max": 1e-2, "steps": 10,
                             "scale": "log"},
                     y_axis={"param": "pull", "min": 5.0, "max": 25.0, "steps": 8,
                             "scale": "linear"}))
    assert data["units"]["limit_of_detection"] == "nM"
    values = [v for row in data["z"] for v in row if v is not None]
    # 2x crossing for a 1 nM binder lands in the pM-to-nM range, never in M.
    assert values and all(1e-4 < v < 1e4 for v in values)


def test_lod_without_a_target_binder_is_refused():
    response = _post(_body(metric="limit_of_detection"))
    assert response.status_code == 400
    assert response.json()["error"]["field"] == "fixed.k_target"


def test_undefined_cells_serialize_as_null_not_nan():
    """Bare NaN is invalid JSON."""
    data = _ok(_body(metric="limit_of_detection", fixed=dict(_FIXED, k_target=1.0),
                     x_axis={"param": "k_open", "min": 1e-9, "max": 1e-7, "steps": 6,
                             "scale": "log"},
                     y_axis={"param": "pull", "min": 0.0, "max": 1.0, "steps": 5,
                             "scale": "linear"}))
    assert "NaN" not in _post(_body(metric="limit_of_detection",
                                    fixed=dict(_FIXED, k_target=1.0))).text
    assert any(v is None for row in data["z"] for v in row)


# --- constraints, feasibility, robust optimum -------------------------------


def test_feasible_mask_and_robust_optimum_reported():
    data = _ok(_body(constraints=[{"metric": "fold_change", "op": ">=", "value": 20.0},
                                  {"metric": "off_signal", "op": "<=", "value": 0.05}]))
    mask = data["feasible"]
    assert mask and any(any(row) for row in mask)

    robust = data["robust_optimum"]
    assert robust is not None
    i = data["x"]["values"].index(robust["x"])
    j = data["y"]["values"].index(robust["y"])
    assert mask[i][j]
    assert robust["x_slack"]["plus"] >= 0 and robust["y_slack"]["minus"] >= 0
    assert "Slack before the constraints stop holding" in robust["summary"]


def test_no_constraints_means_no_mask_and_no_optimum():
    data = _ok(_body())
    assert data["feasible"] is None and data["robust_optimum"] is None


def test_unsatisfiable_constraints_give_no_optimum_but_still_a_grid():
    data = _ok(_body(constraints=[{"metric": "fold_change", "op": ">=", "value": 1e9}]))
    assert data["robust_optimum"] is None
    assert not any(any(row) for row in data["feasible"])
    assert len(data["z"]) == 40


def test_constraint_on_an_unplotted_metric_is_evaluated():
    loose = _ok(_body(constraints=[{"metric": "fraction_open", "op": ">=", "value": 0.0}]))
    tight = _ok(_body(constraints=[{"metric": "fraction_open", "op": ">=", "value": 0.9}]))
    count = lambda d: sum(sum(row) for row in d["feasible"])
    assert count(loose) > count(tight) > 0


def test_bad_constraint_op_refused():
    response = _post(_body(constraints=[{"metric": "fold_change", "op": "~=", "value": 5}]))
    assert response.status_code == 400
    assert "constraint op" in response.json()["error"]["message"]


# --- temperature as an axis ------------------------------------------------


def test_temperature_axis_in_isothermal_mode_is_refused():
    """No dH_ref means nothing to plot."""
    response = _post(_body(x_axis={"param": "temperature_K", "min": 277.15, "max": 320.0,
                                   "steps": 30, "scale": "linear"}))
    assert response.status_code == 400
    error = response.json()["error"]
    assert error["code"] == "ISOTHERMAL_NO_TEMPERATURE_AXIS"
    assert error["field"] == "x_axis.param"
    assert "dH_ref" in error["message"]
    assert "z" not in response.json(), "a refusal must not come with a surface"


def test_temperature_axis_refusal_names_the_y_axis_when_that_is_the_culprit():
    response = _post(_body(y_axis={"param": "temperature_K", "min": 277.15, "max": 320.0,
                                   "steps": 10, "scale": "linear"}))
    assert response.status_code == 400
    assert response.json()["error"]["field"] == "y_axis.param"


def test_temperature_axis_refused_when_only_one_constant_is_isothermal():
    response = _post(_body(x_axis={"param": "temperature_K", "min": 277.15, "max": 320.0,
                                   "steps": 10, "scale": "linear"},
                           thermo={"dG_latch_cage": {"dH_ref": -12.0}}))
    assert response.status_code == 400
    assert "dG_key_cage" in response.json()["error"]["message"]


def test_temperature_axis_with_dH_gives_a_real_surface():
    data = _ok(_body(x_axis={"param": "temperature_K", "min": 277.15, "max": 330.0,
                             "steps": 40, "scale": "linear"},
                     y_axis={"param": "k_open", "min": 1e-5, "max": 1e-1, "steps": 20,
                             "scale": "log"},
                     thermo=_THERMO))
    column = [row[0] for row in data["z"]]
    assert max(column) - min(column) > 1e-6, "temperature surface must not be flat"
    assert data["thermo_mode"] == "vant_hoff"
    assert data["temperature_K"] is None      # temperature is an axis, not a scalar
    assert data["x"]["unit"] == "K"


def test_tm_band_is_reported_on_the_temperature_axis():
    data = _ok(_body(x_axis={"param": "temperature_K", "min": 277.15, "max": 330.0,
                             "steps": 20, "scale": "linear"},
                     y_axis={"param": "k_open", "min": 1e-5, "max": 1e-1, "steps": 10,
                             "scale": "log"},
                     thermo=_THERMO))
    band = data["tm_band"]
    assert band["axis"] == "x"
    assert band["Tm"] == pytest.approx(335.0)
    assert band["min"] == pytest.approx(325.0) and band["max"] == pytest.approx(330.0)
    assert "two-state model invalid" in band["label"]


def test_near_tm_warning_is_not_repeated_per_cell():
    data = _ok(_body(x_axis={"param": "temperature_K", "min": 320.0, "max": 332.0,
                             "steps": 60, "scale": "linear"},
                     y_axis={"param": "k_open", "min": 1e-5, "max": 1e-1, "steps": 20,
                             "scale": "log"},
                     thermo=_THERMO))
    assert any("Tm" in w for w in data["warnings"])
    assert len(data["warnings"]) <= 6, data["warnings"]


def test_no_tm_band_without_a_tm():
    data = _ok(_body(x_axis={"param": "temperature_K", "min": 277.15, "max": 320.0,
                             "steps": 10, "scale": "linear"},
                     y_axis={"param": "k_open", "min": 1e-5, "max": 1e-1, "steps": 10,
                             "scale": "log"},
                     thermo={"dG_latch_cage": {"dH_ref": -12.0}, "dG_key_cage": {"dH_ref": -8.0}}))
    assert data["tm_band"] is None


# --- caps ------------------------------------------------------------------


def test_oversized_grid_is_capped_with_a_warning():
    data = _ok(_body(x_axis=dict(_X, steps=400), y_axis=dict(_Y, steps=400)))
    assert data["cells"] <= MAX_CELLS
    assert len(data["z"]) * len(data["z"][0]) == data["cells"]
    assert any("cap" in w for w in data["warnings"])


def test_lod_grid_capped_tighter():
    data = _ok(_body(metric="limit_of_detection", fixed=dict(_FIXED, k_target=1.0),
                     x_axis=dict(_X, steps=200), y_axis=dict(_Y, steps=200)))
    assert data["cells"] <= MAX_CELLS_LOD
    assert any("limit_of_detection" in w for w in data["warnings"])


def test_steps_above_the_schema_limit_refused():
    response = _post(_body(x_axis=dict(_X, steps=5000)))
    assert response.status_code == 400
    assert "steps" in response.json()["error"]["message"]


# --- validation ------------------------------------------------------------


def test_same_parameter_on_both_axes_refused():
    response = _post(_body(y_axis=dict(_X)))
    assert response.status_code == 400
    assert "different" in response.json()["error"]["message"]


def test_missing_fixed_parameter_refused_by_name():
    response = _post(_body(fixed={"k_ck": 10.0}))
    assert response.status_code == 400
    message = response.json()["error"]["message"]
    assert "luckey" in message and "pull" not in message  # pull is the y axis


def test_log_axis_through_zero_refused():
    response = _post(_body(x_axis=dict(_X, min=0.0)))
    assert response.status_code == 400
    assert "log scale" in response.json()["error"]["message"]


def test_unknown_metric_refused():
    response = _post(_body(metric="brightness"))
    assert response.status_code == 400
    assert "metric must be one of" in response.json()["error"]["message"]


def test_unknown_axis_param_refused():
    response = _post(_body(x_axis=dict(_X, param="humidity")))
    assert response.status_code == 400
    assert "axis param" in response.json()["error"]["message"]


def test_out_of_range_temperature_refused():
    response = _post(_body(temperature_K=600.0))
    assert response.status_code == 400
    assert "temperature_K" in response.json()["error"]["message"]


# --- the 1D slice ----------------------------------------------------------


_SLICE = "/api/landscape/slice"


def test_slice_returns_a_plain_curve():
    data = _ok({"axis": {"param": "pull", "min": 1.0, "max": 30.0, "steps": 12,
                         "scale": "linear"}, "fixed": _FIXED}, _SLICE)
    assert len(data["values"]) == 12
    assert data["axis"]["param"] == "pull"
    assert all(v is not None for v in data["values"])


def test_slice_matches_the_landscape_row_at_the_same_point():
    """Slice must not drift from the grid."""
    curve = _ok({"axis": {"param": "pull", "min": 1.0, "max": 30.0, "steps": 12,
                          "scale": "linear"}, "fixed": _FIXED}, _SLICE)
    grid = _ok(_body(x_axis={"param": "pull", "min": 1.0, "max": 30.0, "steps": 12,
                             "scale": "linear"},
                     y_axis={"param": "k_open", "min": 1e-3, "max": 1e-2, "steps": 2,
                             "scale": "log"}))
    assert [row[0] for row in grid["z"]] == pytest.approx(curve["values"], rel=1e-12)


def test_slice_over_temperature_refused_in_isothermal_mode():
    response = _post({"axis": {"param": "temperature_K", "min": 277.15, "max": 320.0,
                               "steps": 20, "scale": "linear"}, "fixed": _FIXED}, _SLICE)
    assert response.status_code == 400
    assert response.json()["error"]["code"] == "ISOTHERMAL_NO_TEMPERATURE_AXIS"
    assert response.json()["error"]["field"] == "axis.param"


def test_slice_over_temperature_with_dH_works_and_bands_tm():
    data = _ok({"axis": {"param": "temperature_K", "min": 277.15, "max": 330.0,
                         "steps": 30, "scale": "linear"},
                "fixed": _FIXED, "thermo": _THERMO}, _SLICE)
    assert max(data["values"]) - min(data["values"]) > 1e-6
    assert data["tm_band"]["Tm"] == pytest.approx(335.0)
    assert data["temperature_K"] is None


def test_slice_every_metric():
    for metric in ("fold_change", "on_signal", "off_signal", "dynamic_range",
                   "fraction_open", "limit_of_detection"):
        data = _ok({"axis": {"param": "pull", "min": 1.0, "max": 30.0, "steps": 6,
                             "scale": "linear"},
                    "fixed": dict(_FIXED, k_target=1.0), "metric": metric}, _SLICE)
        assert len(data["values"]) == 6, metric


def test_slice_missing_fixed_parameter_refused():
    response = _post({"axis": {"param": "pull", "min": 1.0, "max": 30.0, "steps": 6,
                               "scale": "linear"}, "fixed": {"k_ck": 10.0}}, _SLICE)
    assert response.status_code == 400
    assert "k_open" in response.json()["error"]["message"]


# --- the Module C bridge ---------------------------------------------------


def test_module_c_solutions_map_onto_landscape_axes():
    """dG -> K mapping the frontend uses."""
    target = 8.0      # reachable at pull = 10; 20x needs a stronger pull than that
    inverse = client.post("/api/inverse_design", json={
        "solve_for": ["dG_latch_cage"],
        "fixed": {"pull": 10.0, "dG_key_cage": dG_from_k(1e-8, T_REF_DEFAULT),
                  "key_conc": 500.0},
        "bounds": {"dG_latch_cage": {"min": -10.0, "max": -1.0}},
        "target_fold_change": target,
        "temperature_K": T_REF_DEFAULT,
    })
    assert inverse.status_code == 200, inverse.json()
    result = inverse.json()
    assert result["feasible"], result["message"]

    solution = result["solutions"][0]
    dG = solution["values"]["dG_latch_cage"]
    k_open = math.exp(dG / (R_KCAL * T_REF_DEFAULT))
    assert k_open == pytest.approx(k_from_dG(dG, T_REF_DEFAULT), rel=1e-12)

    # A one-cell-wide slice at the mapped k_open must reproduce Module C's fold change.
    curve = _ok({"axis": {"param": "pull", "min": 10.0, "max": 10.0001, "steps": 2,
                          "scale": "linear"},
                 "fixed": {"k_ck": 10.0, "k_open": k_open, "luckey": 500.0, "pull": 10.0}},
                _SLICE)
    assert curve["values"][0] == pytest.approx(solution["fold_change"], rel=1e-6)
    assert curve["values"][0] == pytest.approx(target, rel=0.05)
