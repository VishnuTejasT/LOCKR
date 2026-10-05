"""Module D engine: the grid is one numpy pass, and it agrees with Module B cell by cell."""

import math

import numpy as np
import pytest

from lockr.engine import landscape, thermo
from lockr.engine.landscape import axis_values, cap_steps, contour_levels, contours, evaluate
from lockr.engine.models import SensorParams
from lockr.engine.temperature import T_REF_DEFAULT, ThermoSpec, shift_constant

_ISO = {k: ThermoSpec(label=k) for k in ("dG_latch_cage", "dG_key_cage", "dG_target_binder")}
_LATCH_FULL = ThermoSpec(dH_ref=-12.0, dCp=0.4, T_ref=T_REF_DEFAULT, Tm=335.0,
                         label="dG_latch_cage")

_BASE = {"k_ck": 1e-8, "k_open": 1e-3, "luckey": 500e-9, "pull": 10.0,
         "temperature_K": T_REF_DEFAULT}


def _xy_grid(nx=23, ny=17):
    xs = axis_values(1e-5, 1e-1, nx, "log")
    ys = axis_values(1.0, 30.0, ny, "linear")
    values = dict(_BASE, k_open=xs[:, None], pull=ys[None, :])
    return xs, ys, values


# --- the vectorization trick ------------------------------------------------


@pytest.mark.parametrize("T", [277.15, 298.15, 310.15, 330.0])
@pytest.mark.parametrize("k", [1e-6, 1e-3, 0.2])
def test_shift_factor_is_multiplicative(T, k):
    """The trick the whole grid rests on."""
    factor, _ = landscape._shift_factor(T, _LATCH_FULL)
    direct, _ = shift_constant(k, T, _LATCH_FULL)
    assert float(factor) * k == pytest.approx(direct, rel=1e-12)


def test_shift_factor_over_a_temperature_axis_matches_each_temperature():
    Ts = axis_values(280.0, 330.0, 11, "linear")
    factors, _ = landscape._shift_factor(Ts, _LATCH_FULL)
    for T, factor in zip(Ts, factors):
        direct, _ = shift_constant(1e-3, float(T), _LATCH_FULL)
        assert float(factor) * 1e-3 == pytest.approx(direct, rel=1e-12)


def test_isothermal_spec_leaves_constants_untouched():
    factor, _ = landscape._shift_factor(310.0, _ISO["dG_latch_cage"])
    assert float(factor) == 1.0


# --- grid vs Module B ------------------------------------------------------


def test_grid_matches_pointwise_module_b():
    """Sampled cells equal a point-wise Module B call."""
    xs, ys, values = _xy_grid()
    grid, _ = evaluate(values, "fold_change", _ISO)
    assert grid.shape == (len(xs), len(ys))

    rng = np.random.default_rng(0)
    for _ in range(50):
        i, j = int(rng.integers(len(xs))), int(rng.integers(len(ys)))
        params = SensorParams(K_open=float(xs[i]), K_CK=_BASE["k_ck"], lucKey=_BASE["luckey"])
        expected = thermo.max_fold_change(Kd=1.0, pull=float(ys[j]), params=params)
        assert grid[i, j] == pytest.approx(expected, rel=1e-12)


def test_sub_saturating_grid_matches_pointwise_module_b():
    xs, ys, values = _xy_grid(13, 11)
    values = dict(values, k_target=1e-9, target_conc=5e-10)
    grid, _ = evaluate(values, "fold_change", _ISO)
    for i in (0, 5, 12):
        for j in (0, 4, 10):
            params = SensorParams(K_open=float(xs[i]), K_CK=_BASE["k_ck"], lucKey=_BASE["luckey"])
            expected = thermo.fold_change(5e-10, 1e-9, float(ys[j]), params)
            assert grid[i, j] == pytest.approx(expected, rel=1e-12)


def test_lod_grid_matches_module_b_lod():
    """Must reproduce thermo.lod_and_ec50's 2x crossing."""
    xs, ys, values = _xy_grid(9, 7)
    values = dict(values, k_target=1e-9)
    grid, _ = evaluate(values, "limit_of_detection", _ISO)
    for i in range(len(xs)):
        for j in range(len(ys)):
            params = SensorParams(K_open=float(xs[i]), K_CK=_BASE["k_ck"], lucKey=_BASE["luckey"])
            expected = thermo.lod_and_ec50(1e-9, float(ys[j]), params).lod_2x
            if expected is None:
                assert math.isnan(grid[i, j])
            else:
                assert grid[i, j] == pytest.approx(expected, rel=1e-12)


def test_off_signal_matches_f_base_and_fold_change_is_their_ratio():
    xs, ys, values = _xy_grid(7, 5)
    on, _ = evaluate(values, "on_signal", _ISO)
    off, _ = evaluate(values, "off_signal", _ISO)
    fc, _ = evaluate(values, "fold_change", _ISO)
    dyn, _ = evaluate(values, "dynamic_range", _ISO)
    assert np.allclose(fc, on / off, rtol=1e-12)
    assert np.allclose(dyn, on - off, rtol=1e-12)
    for i in range(len(xs)):
        params = SensorParams(K_open=float(xs[i]), K_CK=_BASE["k_ck"], lucKey=_BASE["luckey"])
        assert off[i, 0] == pytest.approx(thermo.f_base(params), rel=1e-12)


def test_fraction_open_is_at_least_on_signal_and_bounded():
    _, _, values = _xy_grid(11, 9)
    frac, _ = evaluate(values, "fraction_open", _ISO)
    on, _ = evaluate(values, "on_signal", _ISO)
    assert np.all(frac >= on - 1e-12)      # every glowing cage is also an open cage
    assert np.all((frac >= 0) & (frac <= 1))


def test_metric_needing_target_refuses_without_one():
    _, _, values = _xy_grid(5, 5)
    with pytest.raises(ValueError, match="k_target"):
        evaluate(values, "limit_of_detection", _ISO)


def test_unknown_metric_refused():
    _, _, values = _xy_grid(5, 5)
    with pytest.raises(ValueError, match="metric must be one of"):
        evaluate(values, "glow_factor", _ISO)


# --- temperature as an axis ------------------------------------------------


def test_temperature_axis_is_not_flat_with_dH():
    Ts = axis_values(277.15, 330.0, 21, "linear")
    xs = axis_values(1e-5, 1e-1, 9, "log")
    values = dict(_BASE, k_open=xs[:, None], temperature_K=Ts[None, :])
    grid, warns = evaluate(values, "fold_change", dict(_ISO, dG_latch_cage=_LATCH_FULL))
    assert grid.shape == (9, 21)
    assert np.ptp(grid[4, :]) > 1e-6
    assert any("Tm" in w for w in warns)


def test_temperature_axis_is_flat_without_dH():
    """Isothermal really is flat; the API refuses to draw it."""
    Ts = axis_values(277.15, 330.0, 11, "linear")
    values = dict(_BASE, temperature_K=Ts)
    grid, _ = evaluate(values, "fold_change", _ISO)
    assert np.ptp(grid) == pytest.approx(0.0, abs=1e-15)


def test_temperature_axis_warnings_come_from_the_ends_only():
    specs = {"dG_latch_cage": _LATCH_FULL}
    warns = landscape.temperature_axis_warnings(specs, 277.15, 330.0)
    assert warns and len(warns) <= 4, warns   # not one per cell


def test_tm_band_only_when_in_range():
    assert landscape.tm_band(_LATCH_FULL, 277.0, 300.0) is None
    band = landscape.tm_band(_LATCH_FULL, 277.0, 330.0)
    assert band["min"] == pytest.approx(325.0) and band["max"] == pytest.approx(330.0)
    assert landscape.tm_band(_ISO["dG_latch_cage"], 277.0, 400.0) is None


# --- feasibility and the robust optimum ------------------------------------


_CONSTRAINTS = [{"metric": "fold_change", "op": ">=", "value": 20.0},
                {"metric": "off_signal", "op": "<=", "value": 0.05}]


def _feasible_setup():
    xs, ys, values = _xy_grid(31, 25)
    fc, _ = evaluate(values, "fold_change", _ISO)
    off, _ = evaluate(values, "off_signal", _ISO)
    mask = landscape.feasibility({"fold_change": fc, "off_signal": off}, _CONSTRAINTS)
    return xs, ys, fc, off, mask


def test_feasibility_matches_independent_per_cell_evaluation():
    xs, ys, fc, off, mask = _feasible_setup()
    expected = np.array([[bool(fc[i, j] >= 20.0 and off[i, j] <= 0.05)
                          for j in range(len(ys))] for i in range(len(xs))])
    assert (mask == expected).all()
    assert mask.any() and not mask.all()      # a real boundary exists to be robust against


def test_nan_cells_are_never_feasible():
    grid = np.array([[1.0, np.nan], [50.0, np.nan]])
    mask = landscape.feasibility({"fold_change": grid},
                                 [{"metric": "fold_change", "op": ">=", "value": 0.0}])
    assert mask.tolist() == [[True, False], [True, False]]


def test_no_constraints_means_no_mask():
    assert landscape.feasibility({"fold_change": np.ones((3, 3))}, []) is None


def test_robust_optimum_is_feasible_and_its_slack_is_correct():
    xs, ys, fc, _, mask = _feasible_setup()
    found = landscape.robust_optimum(mask, xs, ys, fc)
    assert found is not None

    i = list(xs).index(found["x"])
    j = list(ys).index(found["y"])
    assert mask[i, j], "robust optimum must be inside the feasible mask"

    # one step beyond the reported slack must be infeasible, or the axis edge
    for slack, index, values, line in (
            (found["x_slack"], i, xs, mask[:, j]),
            (found["y_slack"], j, ys, mask[i, :])):
        lo_val = values[index] - slack["minus"]
        hi_val = values[index] + slack["plus"]
        lo_i = int(np.argmin(np.abs(values - lo_val)))
        hi_i = int(np.argmin(np.abs(values - hi_val)))
        assert line[lo_i] and line[hi_i]
        assert lo_i == 0 or not line[lo_i - 1]
        assert hi_i == len(line) - 1 or not line[hi_i + 1]


def test_robust_optimum_is_further_from_the_boundary_than_other_feasible_cells():
    xs, ys, fc, _, mask = _feasible_setup()
    found = landscape.robust_optimum(mask, xs, ys, fc)
    i = list(xs).index(found["x"])
    j = list(ys).index(found["y"])

    # smallest one-axis slack in cells is what "most robust" maximises
    def cell_margin(a, b):
        col, row = mask[:, b], mask[a, :]
        out = []
        for line, here in ((col, a), (row, b)):
            lo = hi = here
            while lo > 0 and line[lo - 1]:
                lo -= 1
            while hi < len(line) - 1 and line[hi + 1]:
                hi += 1
            out += [here - lo, hi - here]
        return min(out)

    best = cell_margin(i, j)
    for a in range(len(xs)):
        for b in range(len(ys)):
            if mask[a, b]:
                assert cell_margin(a, b) <= best + 1


def test_robust_optimum_is_none_when_nothing_is_feasible():
    xs, ys, _ = _xy_grid(5, 5)
    empty = np.zeros((5, 5), dtype=bool)
    assert landscape.robust_optimum(empty, xs, ys, np.ones((5, 5))) is None
    assert landscape.robust_optimum(None, xs, ys, np.ones((5, 5))) is None


# --- contours --------------------------------------------------------------


def test_contours_close_or_terminate_on_the_grid_edge():
    xs, ys, values = _xy_grid(41, 33)
    grid, _ = evaluate(values, "fold_change", _ISO)
    levels = contour_levels(grid)
    assert levels

    traced = contours(grid, xs, ys, levels)
    assert traced

    def on_edge(point):
        return (math.isclose(point[0], xs[0]) or math.isclose(point[0], xs[-1])
                or math.isclose(point[1], ys[0]) or math.isclose(point[1], ys[-1]))

    for contour in traced:
        for polyline in contour["polylines"]:
            assert len(polyline) >= 2, "a polyline must have at least one segment"
            closed = (math.isclose(polyline[0][0], polyline[-1][0])
                      and math.isclose(polyline[0][1], polyline[-1][1]))
            assert closed or (on_edge(polyline[0]) and on_edge(polyline[-1])), (
                f"dangling contour at level {contour['level']}: "
                f"{polyline[0]} -> {polyline[-1]}")


def test_contour_vertices_sit_on_their_level_for_a_known_analytic_case():
    """Marching squares on a plane is exact."""
    xs = np.linspace(0.0, 10.0, 51)
    ys = np.linspace(0.0, 10.0, 41)
    plane = xs[:, None] + 2.0 * ys[None, :]
    traced = contours(plane, xs, ys, [10.0, 20.0])
    assert traced
    for contour in traced:
        for polyline in contour["polylines"]:
            for x, y in polyline:
                assert x + 2.0 * y == pytest.approx(contour["level"], abs=1e-9)


def test_contours_skip_cells_touching_an_undefined_value():
    grid = np.array([[1.0, 2.0, 3.0], [2.0, np.nan, 4.0], [3.0, 4.0, 5.0]])
    xs = ys = np.array([0.0, 1.0, 2.0])
    for contour in contours(grid, xs, ys, [3.0]):
        for polyline in contour["polylines"]:
            for x, y in polyline:
                assert not (0.0 < x < 2.0 and 0.0 < y < 2.0) or (x in (1.0,) or y in (1.0,))


def test_contour_levels_stay_inside_the_data_range():
    grid = np.array([[1.0, 5.0], [10.0, 30.5]])
    for level in contour_levels(grid):
        assert 1.0 < level < 30.5


def test_contour_levels_empty_on_a_flat_grid():
    assert contour_levels(np.ones((4, 4))) == []
    assert contour_levels(np.full((3, 3), np.nan)) == []


# --- axes and caps ---------------------------------------------------------


def test_axis_values_respect_scale():
    log = axis_values(1e-5, 1e-1, 5, "log")
    assert log[0] == pytest.approx(1e-5) and log[-1] == pytest.approx(1e-1)
    assert np.allclose(np.diff(np.log10(log)), 1.0)
    linear = axis_values(1.0, 5.0, 5, "linear")
    assert np.allclose(linear, [1, 2, 3, 4, 5])


@pytest.mark.parametrize("kwargs,match", [
    ({"lo": 5.0, "hi": 1.0, "steps": 5, "scale": "linear"}, "min must be < max"),
    ({"lo": 1.0, "hi": 5.0, "steps": 1, "scale": "linear"}, "steps must be >= 2"),
    ({"lo": 0.0, "hi": 5.0, "steps": 5, "scale": "log"}, "log axis min must be > 0"),
    ({"lo": 1.0, "hi": 5.0, "steps": 5, "scale": "rainbow"}, "scale must be one of"),
])
def test_axis_values_validation(kwargs, match):
    with pytest.raises(ValueError, match=match):
        axis_values(kwargs["lo"], kwargs["hi"], kwargs["steps"], kwargs["scale"])


def test_cap_leaves_a_small_grid_alone():
    assert cap_steps(80, 60, "fold_change") == (80, 60, [])


def test_cap_shrinks_an_oversized_grid_and_warns():
    nx, ny, warns = cap_steps(400, 400, "fold_change")
    assert nx * ny <= landscape.MAX_CELLS and warns and "cap" in warns[0]


def test_lod_has_a_tighter_cap_than_the_other_metrics():
    nx, ny, warns = cap_steps(200, 200, "limit_of_detection")
    assert nx * ny <= landscape.MAX_CELLS_LOD and warns
    assert cap_steps(200, 200, "fold_change") == (200, 200, [])
