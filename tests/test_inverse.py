"""Module C Part 2: the inverse solver. The core invariant is forward(inverse(FC)) == FC."""

import math

import numpy as np
import pytest

from lockr.engine import inverse
from lockr.engine.inverse import MODE_REGION, MODE_ROBUST_POINT, MODE_ROOT, MODE_ROOTS, TargetSpec
from lockr.engine.temperature import T_REF_DEFAULT, ThermoSpec, dG_from_k

Tr = T_REF_DEFAULT

_BOUNDS = {
    "dG_latch_cage": (-10.0, -1.0),
    "dG_key_cage": (-15.0, -6.0),
    "dG_target_binder": (-15.0, -6.0),
    "key_conc": (1e-9, 5e-6),
    "target_conc": (1e-12, 1e-6),
}


def _fixed(pull=20.0, **over):
    base = {"pull": pull,
            "dG_latch_cage": dG_from_k(1e-3, Tr),
            "dG_key_cage": dG_from_k(1e-8, Tr),
            "key_conc": 500e-9}
    base.update(over)
    return base


def _forward(fixed, temperature_K=Tr, specs=None):
    base = {k: v for k, v in fixed.items() if k != "pull"}
    base.setdefault("target_conc", None)
    return inverse.make_forward(base, fixed["pull"], temperature_K, specs)


def _free(fixed, *names):
    return {k: v for k, v in fixed.items() if k not in names}


# --- the round-trip invariant ------------------------------------------------

@pytest.mark.parametrize("name", list(_BOUNDS))
def test_round_trip_for_random_valid_parameter_sets(name):
    rng = np.random.default_rng(7)
    checked = 0
    for _ in range(12):
        fixed = _fixed(
            pull=float(rng.uniform(5.0, 40.0)),
            dG_latch_cage=dG_from_k(float(10 ** rng.uniform(-5, -1.5)), Tr),
            dG_key_cage=dG_from_k(float(10 ** rng.uniform(-10, -7)), Tr),
            key_conc=float(10 ** rng.uniform(-8.5, -6.0)),
        )
        if name in ("dG_target_binder", "target_conc"):
            fixed["dG_target_binder"] = dG_from_k(float(10 ** rng.uniform(-10, -7)), Tr)
            fixed["target_conc"] = float(10 ** rng.uniform(-10, -7))

        lo, hi = _BOUNDS[name]
        truth = float(rng.uniform(lo, hi)) if name.startswith("dG") else \
            float(10 ** rng.uniform(math.log10(lo), math.log10(hi)))
        target = _forward({**fixed, name: truth})()
        if not (1.05 < target):  # a flat spot carries no information to invert
            continue

        result = inverse.inverse_design([name], _free(fixed, name), {name: _BOUNDS[name]},
                                        target_fold_change=target)
        assert result.feasible, result.message
        assert result.solutions
        for solution in result.solutions:
            # forward(inverse(FC)) == FC, the invariant this whole module exists for.
            assert solution.fold_change == pytest.approx(target, rel=1e-6)
            assert lo <= solution.values[name] <= hi
            assert solution.on_spec
        checked += 1
    assert checked >= 5


def test_round_trip_holds_at_a_non_default_temperature():
    # Every energy in play needs its own dH/dCp, or the weakest mode drags the label down.
    specs = {"dG_latch_cage": ThermoSpec(dH_ref=-20.0, dCp=-0.6, label="dG_latch_cage"),
             "dG_key_cage": ThermoSpec(dH_ref=-15.0, dCp=-0.4, label="dG_key_cage")}
    fixed = _fixed()
    target = _forward(fixed, 310.0, specs)({"dG_latch_cage": -5.0})
    result = inverse.inverse_design(["dG_latch_cage"], _free(fixed, "dG_latch_cage"),
                                    {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                                    target_fold_change=target, temperature_K=310.0, specs=specs)
    assert result.feasible
    assert result.thermo_mode == "full"
    assert result.solutions[0].values["dG_latch_cage"] == pytest.approx(-5.0, rel=1e-6)
    assert result.solutions[0].fold_change == pytest.approx(target, rel=1e-6)


# --- one free parameter ------------------------------------------------------

def test_single_free_parameter_reports_monotonicity_and_one_root():
    fixed = _fixed()
    result = inverse.inverse_design(["dG_latch_cage"], _free(fixed, "dG_latch_cage"),
                                    {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                                    target_fold_change=8.0)
    assert result.mode == MODE_ROOT
    assert len(result.solutions) == 1
    assert "monotonic" in result.message
    solution = result.solutions[0]
    assert solution.fold_change == pytest.approx(8.0)
    lo, hi = solution.tolerance["dG_latch_cage"]
    assert lo < solution.values["dG_latch_cage"] < hi
    assert solution.sensitivity["dG_latch_cage"] != 0


def test_every_fold_change_inside_the_tolerance_interval_meets_spec():
    fixed = _fixed()
    result = inverse.inverse_design(["key_conc"], _free(fixed, "key_conc"),
                                    {"key_conc": _BOUNDS["key_conc"]}, target_fold_change=8.0,
                                    rel_tol=0.05)
    solution = result.solutions[0]
    spec = TargetSpec.build(8.0, rel_tol=0.05)
    forward = _forward(fixed)
    lo, hi = solution.tolerance["key_conc"]
    for x in np.linspace(lo, hi, 25):
        assert spec.contains(forward({"key_conc": float(x)}))
    # and just outside it does not
    assert not spec.contains(forward({"key_conc": lo * 0.9}))
    assert not spec.contains(forward({"key_conc": hi * 1.1}))


def test_non_monotonic_evaluator_returns_every_root():
    """The LucCage model happens to be monotonic in each parameter on its own, so the
    multi-root branch is exercised with an evaluator that is not."""
    forward = lambda over: 5.0 + 3.0 * math.sin(over["key_conc"])
    xs = np.linspace(0.0, 20.0, 401)
    assert not inverse._is_monotonic(np.array([forward({"key_conc": x}) for x in xs]))

    out = inverse._solve_1d(forward, {"key_conc": 1.0}, "key_conc", {"key_conc": (0.0, 20.0)},
                            TargetSpec.build(5.0, rel_tol=0.02))
    assert out["mode"] == MODE_ROOTS
    assert "NOT monotonic" in out["message"]
    roots = [s.values["key_conc"] for s in out["solutions"]]
    assert len(roots) >= 6
    assert all(forward({"key_conc": r}) == pytest.approx(5.0, abs=1e-6) for r in roots)
    assert roots == sorted(roots)


def test_root_bracketing_finds_all_sign_changes():
    f = lambda x: (x - 1.0) * (x - 3.0) * (x - 7.0)
    roots = inverse._bracket_roots(f, np.linspace(0.0, 10.0, 101), 0.0)
    assert roots == pytest.approx([1.0, 3.0, 7.0])


# --- unreachable targets ----------------------------------------------------

def test_unreachable_target_reports_the_closest_point_and_the_binding_constraint():
    fixed = _fixed(pull=10.0)
    result = inverse.inverse_design(["dG_latch_cage"], _free(fixed, "dG_latch_cage"),
                                    {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                                    target_fold_change=1e6)
    assert result.feasible is False
    assert result.solutions == []
    assert result.closest is not None
    # pull = 10 caps fold change near 1 + pull, whatever the latch does.
    assert result.closest.fold_change < 12.0
    assert "dG_latch_cage" in result.closest.binding_constraint
    assert "not reachable" in result.message
    assert result.closest.note


def test_unreachable_target_in_two_dimensions_still_takes_the_infeasible_path():
    fixed = _fixed(pull=10.0)
    result = inverse.inverse_design(
        ["dG_latch_cage", "key_conc"], _free(fixed, "dG_latch_cage", "key_conc"),
        {k: _BOUNDS[k] for k in ("dG_latch_cage", "key_conc")}, target_fold_change=500.0)
    assert result.feasible is False
    assert result.region == []
    assert result.closest is not None and result.closest.fold_change < 12.0


def test_unreachable_target_in_three_dimensions_still_takes_the_infeasible_path():
    result = inverse.inverse_design(
        ["dG_latch_cage", "dG_key_cage", "key_conc"], {"pull": 10.0},
        {k: _BOUNDS[k] for k in ("dG_latch_cage", "dG_key_cage", "key_conc")},
        target_fold_change=500.0)
    assert result.feasible is False
    assert result.mode == MODE_ROBUST_POINT
    assert result.closest is not None
    assert result.closest.binding_constraint


# --- two free parameters ----------------------------------------------------

def test_two_free_parameters_return_a_region_whose_every_point_hits_the_target():
    fixed = _fixed()
    names = ["dG_latch_cage", "key_conc"]
    result = inverse.inverse_design(names, _free(fixed, *names),
                                    {k: _BOUNDS[k] for k in names}, target_fold_change=8.0)
    assert result.mode == MODE_REGION
    assert result.feasible
    assert result.region_axes == tuple(names)
    assert result.region, "expected at least one contour polyline"

    forward = _forward(fixed)
    for line in result.region:
        for x, y in line:
            assert forward({"dG_latch_cage": x, "key_conc": y}) == pytest.approx(8.0, rel=1e-6)

    robust = result.robust_point
    assert robust is not None and robust.on_spec
    for name in names:
        lo, hi = robust.tolerance[name]
        assert lo < robust.values[name] < hi


def test_robust_point_is_on_target_and_furthest_from_the_region_boundary():
    fixed = _fixed()
    names = ["dG_latch_cage", "key_conc"]
    bounds = {k: _BOUNDS[k] for k in names}
    spec = TargetSpec.build(8.0, rel_tol=0.05)
    forward = _forward(fixed)
    result = inverse.inverse_design(names, _free(fixed, *names), bounds, target_fold_change=8.0)
    robust = result.robust_point

    # It sits exactly on the target contour, not on a grid pixel near it.
    assert robust.fold_change == pytest.approx(8.0, rel=1e-6)

    # No point on the returned region is further from the boundary than the one chosen.
    margins = [inverse._margin(
        inverse._solution(forward, {**_free(fixed, *names), "target_conc": None,
                                    "dG_latch_cage": x, "key_conc": y}, names, bounds, spec),
        names, bounds) for line in result.region for x, y in line[::7]]
    assert inverse._margin(robust, names, bounds) >= max(margins) - 1e-12

    # A design deliberately parked against the edge of the band tolerates far less drift.
    at_edge = inverse._bracket_roots(
        inverse._slice(forward, {**robust.values, **_free(fixed, *names), "target_conc": None},
                       "key_conc"),
        np.logspace(-9, math.log10(5e-6), 401), spec.lo * 1.001)
    assert at_edge, "expected a point at the lower edge of the accepted band"
    edge_solution = inverse._solution(
        forward, {**_free(fixed, *names), "target_conc": None,
                  "dG_latch_cage": robust.values["dG_latch_cage"], "key_conc": at_edge[0]},
        names, bounds, spec)
    assert edge_solution.on_spec
    assert inverse._margin(edge_solution, names, bounds) < inverse._margin(robust, names, bounds)


# --- three or more free parameters ------------------------------------------

def test_three_free_parameters_return_a_robust_point_with_per_parameter_tolerances():
    names = ["dG_latch_cage", "dG_key_cage", "key_conc"]
    result = inverse.inverse_design(names, {"pull": 20.0}, {k: _BOUNDS[k] for k in names},
                                    target_fold_change=8.0)
    assert result.mode == MODE_ROBUST_POINT
    assert result.feasible
    assert result.region == []
    assert result.region_axes is None
    solution = result.robust_point
    assert solution.fold_change == pytest.approx(8.0, rel=0.05)
    assert set(solution.tolerance) == set(names)
    assert set(solution.sensitivity) == set(names)
    for name in names:
        lo, hi = solution.tolerance[name]
        assert lo <= solution.values[name] <= hi


def test_four_free_parameters_and_a_target_range():
    names = ["dG_latch_cage", "dG_key_cage", "key_conc", "target_conc"]
    result = inverse.inverse_design(
        names, {"pull": 20.0, "dG_target_binder": dG_from_k(1e-9, Tr)},
        {k: _BOUNDS[k] for k in names}, target_fold_change_range=(9.0, 11.0))
    assert result.feasible
    assert 9.0 <= result.robust_point.fold_change <= 11.0


# --- input validation -------------------------------------------------------

@pytest.mark.parametrize("kwargs", [
    dict(solve_for=[], fixed={"pull": 10.0}, bounds={}),
    dict(solve_for=["nonsense"], fixed={"pull": 10.0}, bounds={"nonsense": (0.0, 1.0)}),
    dict(solve_for=["key_conc", "key_conc"], fixed={"pull": 10.0}, bounds={"key_conc": (1e-9, 1e-6)}),
    dict(solve_for=["key_conc"], fixed={}, bounds={"key_conc": (1e-9, 1e-6)}),  # no pull
    dict(solve_for=["key_conc"], fixed={"pull": 10.0}, bounds={}),  # no bounds
    dict(solve_for=["key_conc"], fixed={"pull": 10.0, "key_conc": 1e-7},
         bounds={"key_conc": (1e-9, 1e-6)}),  # free and fixed at once
    dict(solve_for=["key_conc"], fixed={"pull": 10.0}, bounds={"key_conc": (1e-6, 1e-9)}),  # min>max
    dict(solve_for=["key_conc"], fixed={"pull": 10.0}, bounds={"key_conc": (-1e-9, 1e-6)}),
])
def test_bad_requests_raise_value_error(kwargs):
    with pytest.raises(ValueError):
        inverse.inverse_design(target_fold_change=8.0, **kwargs)


def test_missing_fixed_parameter_is_named():
    with pytest.raises(ValueError, match="dG_key_cage"):
        inverse.inverse_design(["dG_latch_cage"], {"pull": 10.0, "key_conc": 5e-7},
                               {"dG_latch_cage": _BOUNDS["dG_latch_cage"]}, target_fold_change=8.0)


def test_exactly_one_target_form_is_required():
    fixed = _fixed()
    for target, target_range in ((None, None), (8.0, (7.0, 9.0))):
        with pytest.raises(ValueError):
            inverse.inverse_design(["key_conc"], _free(fixed, "key_conc"),
                                   {"key_conc": _BOUNDS["key_conc"]},
                                   target_fold_change=target, target_fold_change_range=target_range)


def test_isothermal_warning_is_carried_through_to_the_result():
    fixed = _fixed()
    result = inverse.inverse_design(["dG_latch_cage"], _free(fixed, "dG_latch_cage"),
                                    {"dG_latch_cage": _BOUNDS["dG_latch_cage"]},
                                    target_fold_change=8.0, temperature_K=310.0)
    assert result.thermo_mode == "isothermal"
    assert any("TEMPERATURE SCALING UNAVAILABLE" in w for w in result.warnings)
