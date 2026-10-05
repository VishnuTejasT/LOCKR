"""Inverse design: which parameters land on a wanted fold change.

thermo.py runs forward, parameters in and fold change out. Design work runs the
other way round: the wet lab knows the fold change it needs and wants to know
which dG or which Key dose gets there. How much can honestly be said depends on
how many parameters are free, so there are three answers, not one:

    1 free    root-find (Brent). Monotonicity is checked first, and if the
              function is not monotonic every root on the scan is returned.
    2 free    a feasible region as contour polylines, plus the point furthest
              from the region boundary, i.e. the design least sensitive to
              parameter error.
    3+ free   no region (it is not drawable and not trustworthy at this cost).
              A robust point from multi-start optimisation, plus how far each
              parameter can drift on its own.

Every free parameter is a reference-temperature value: the answer is what to
design, and the temperature layer maps it to the working temperature. An
unreachable target is reported as unreachable, with the closest achievable fold
change and the binding constraint, never as a quiet best effort.

Units: dG in kcal/mol, concentrations in M.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.optimize import brentq, minimize

from . import thermo
from .models import SensorParams
from .temperature import (R_KCAL, T_REF_DEFAULT, ThermoSpec, anchored, dG_at, k_from_dG,
                          warnings_at, weakest_mode)

DG_PARAMS = ("dG_latch_cage", "dG_key_cage", "dG_target_binder")
CONC_PARAMS = ("key_conc", "target_conc")
PARAMS = DG_PARAMS + CONC_PARAMS

MODE_ROOT = "root"
MODE_ROOTS = "multiple_roots"
MODE_REGION = "region"
MODE_ROBUST_POINT = "robust_point"

_DECADES_FOR_LOG_AXIS = 100.0
_SCAN_POINTS = 201
_GRID_POINTS = 41
_CONTOUR_ROWS = 61
_LINK_THRESHOLD = 0.1  # in unit coords, for joining contour points into polylines
_MONO_REL_TOL = 1e-9
_ND_STARTS = 12
_ND_SEED = 0


@dataclass(frozen=True)
class TargetSpec:
    """What counts as hitting the target. A single target becomes a +-rel_tol band,
    because a zero-width target has no feasible region and no tolerance."""

    lo: float
    hi: float
    center: float

    @classmethod
    def build(cls, target: float | None = None, target_range: tuple[float, float] | None = None,
              rel_tol: float = 0.05) -> "TargetSpec":
        if (target is None) == (target_range is None):
            raise ValueError("give exactly one of target_fold_change or target_fold_change_range")
        if target is not None:
            if target <= 0:
                raise ValueError("target_fold_change must be > 0")
            if not 0 < rel_tol < 1:
                raise ValueError("rel_tol must be between 0 and 1")
            return cls(target * (1 - rel_tol), target * (1 + rel_tol), target)
        lo, hi = target_range
        if lo <= 0 or hi <= lo:
            raise ValueError("target_fold_change_range must satisfy 0 < min < max")
        return cls(lo, hi, math.sqrt(lo * hi))

    def contains(self, fc: float) -> bool:
        return bool(math.isfinite(fc)) and self.lo <= fc <= self.hi


@dataclass
class Solution:
    values: dict[str, float]
    fold_change: float
    on_spec: bool
    sensitivity: dict[str, float]
    tolerance: dict[str, tuple[float, float]]


@dataclass
class ClosestPoint:
    values: dict[str, float]
    fold_change: float
    binding_constraint: str
    note: str


@dataclass
class InverseResult:
    mode: str
    feasible: bool
    message: str
    temperature_K: float
    thermo_mode: str
    solutions: list[Solution] = field(default_factory=list)
    region: list[list[tuple[float, float]]] = field(default_factory=list)
    region_axes: tuple[str, str] | None = None
    robust_point: Solution | None = None
    closest: ClosestPoint | None = None
    warnings: list[str] = field(default_factory=list)


def _log_axis(lo: float, hi: float) -> bool:
    return lo > 0 and hi / lo >= _DECADES_FOR_LOG_AXIS


def _axis(lo: float, hi: float, n: int) -> np.ndarray:
    if _log_axis(lo, hi):
        return np.logspace(math.log10(lo), math.log10(hi), n)
    return np.linspace(lo, hi, n)


def _to_unit(x: float, lo: float, hi: float) -> float:
    if _log_axis(lo, hi):
        return (math.log(x) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return (x - lo) / (hi - lo)


def _from_unit(u: float, lo: float, hi: float) -> float:
    u = min(1.0, max(0.0, u))
    if _log_axis(lo, hi):
        return math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))
    return lo + u * (hi - lo)


def make_forward(base: dict, pull: float, temperature_K: float,
                 specs: dict[str, ThermoSpec] | None = None) -> Callable[[dict], float]:
    """Fold change as a function of a dict of parameter overrides.

    Each dG argument is the value at its spec's reference temperature; the
    Gibbs-Helmholtz layer carries it to temperature_K.
    """
    specs = specs or {}
    T = temperature_K

    def forward(over: dict | None = None) -> float:
        v = {**base, **(over or {})}
        try:
            k_open = k_from_dG(dG_at(anchored(specs.get("dG_latch_cage", ThermoSpec()),
                                              v["dG_latch_cage"]), T), T)
            k_ck = k_from_dG(dG_at(anchored(specs.get("dG_key_cage", ThermoSpec()),
                                            v["dG_key_cage"]), T), T)
            params = SensorParams(K_open=k_open, K_CK=k_ck, lucKey=v["key_conc"], RT=R_KCAL * T)
            if v.get("target_conc") is None:
                fc = thermo.max_fold_change(Kd=1.0, pull=pull, params=params)
            else:
                Kd = k_from_dG(dG_at(anchored(specs.get("dG_target_binder", ThermoSpec()),
                                              v["dG_target_binder"]), T), T)
                fc = thermo.fold_change(v["target_conc"], Kd, pull, params)
        except (OverflowError, ZeroDivisionError, ValueError):
            return float("nan")
        return fc if math.isfinite(fc) else float("nan")

    return forward


def _vary(values: dict, name: str, x: float) -> dict:
    out = dict(values)
    out[name] = float(x)
    return out


def _slice(forward, values: dict, name: str) -> Callable[[float], float]:
    return lambda x: forward(_vary(values, name, x))


def _bracket_roots(f: Callable[[float], float], xs: np.ndarray, target: float) -> list[float]:
    """Every sign change of f - target on the scan, refined with Brent."""
    ys = np.array([f(float(x)) for x in xs], dtype=float) - target
    roots: list[float] = []
    for a, b, ya, yb in zip(xs[:-1], xs[1:], ys[:-1], ys[1:]):
        if not (math.isfinite(ya) and math.isfinite(yb)):
            continue
        if ya == 0.0:
            roots.append(float(a))
        elif ya * yb < 0:
            # xtol has to follow the parameter's own scale: 1e-14 absolute is coarse for a
            # concentration around 1e-9 M and would show up as error in the achieved fold change.
            xtol = max(1e-16, 1e-12 * max(abs(a), abs(b)))
            roots.append(float(brentq(lambda x: f(x) - target, a, b, xtol=xtol, rtol=1e-12)))
    if ys.size and ys[-1] == 0.0:
        roots.append(float(xs[-1]))
    return roots


def _is_monotonic(ys: np.ndarray) -> bool:
    finite = ys[np.isfinite(ys)]
    if finite.size < 3:
        return True
    tol = _MONO_REL_TOL * max(1.0, float(np.max(np.abs(finite))))
    d = np.diff(finite)
    return bool(np.all(d >= -tol) or np.all(d <= tol))


def _sensitivity(f1: Callable[[float], float], x: float, lo: float, hi: float) -> float:
    """d(log FC)/d(param), central difference."""
    h = max(1e-9, 1e-5 * (hi - lo))
    a, b = max(lo, x - h), min(hi, x + h)
    fa, fb = f1(a), f1(b)
    if b <= a or not (math.isfinite(fa) and math.isfinite(fb)) or fa <= 0 or fb <= 0:
        return float("nan")
    return (math.log(fb) - math.log(fa)) / (b - a)


def _edge(f1: Callable[[float], float], x0: float, limit: float, spec: TargetSpec,
          coarse: int = 33, refine: int = 60) -> float:
    """Furthest value between x0 and limit that still meets spec."""
    if x0 == limit:
        return x0
    good = x0
    for t in np.linspace(0.0, 1.0, coarse)[1:]:
        x = x0 + float(t) * (limit - x0)
        if spec.contains(f1(x)):
            good = x
            continue
        bad = x
        for _ in range(refine):
            mid = 0.5 * (good + bad)
            if spec.contains(f1(mid)):
                good = mid
            else:
                bad = mid
        return good
    return limit


def _solution(forward, values: dict, names, bounds: dict, spec: TargetSpec) -> Solution:
    fc = forward(values)
    on_spec = spec.contains(fc)
    sens, tol = {}, {}
    for name in names:
        lo, hi = bounds[name]
        f1 = _slice(forward, values, name)
        sens[name] = _sensitivity(f1, values[name], lo, hi)
        x0 = values[name]
        tol[name] = ((x0, x0) if not on_spec
                     else (_edge(f1, x0, lo, spec), _edge(f1, x0, hi, spec)))
    return Solution(values={k: values[k] for k in names}, fold_change=fc, on_spec=on_spec,
                    sensitivity=sens, tolerance=tol)


def _binding_constraint(values: dict, names, bounds: dict, sens: dict) -> str:
    for name in names:
        lo, hi = bounds[name]
        span = abs(hi - lo) or 1.0
        if abs(values[name] - lo) <= 1e-6 * span:
            return f"{name} at its lower bound ({lo:g})"
        if abs(values[name] - hi) <= 1e-6 * span:
            return f"{name} at its upper bound ({hi:g})"
    ranked = [n for n in names if math.isfinite(sens.get(n, float("nan")))]
    if not ranked:
        return "none of the free parameters moves fold change here"
    worst = max(ranked, key=lambda n: abs(sens[n]))
    return (f"interior limit of the model, not a bound: {worst} has the largest remaining "
            f"leverage but the model saturates before the target")


def _closest_from_samples(forward, samples, names, bounds, spec) -> ClosestPoint:
    best, best_err, best_fc = None, float("inf"), float("nan")
    for values in samples:
        fc = forward(values)
        if not math.isfinite(fc) or fc <= 0:
            continue
        err = abs(math.log(fc) - math.log(spec.center))
        if err < best_err:
            best, best_err, best_fc = values, err, fc
    if best is None:
        raise ValueError("fold change is undefined everywhere inside the given bounds")
    sens = {n: _sensitivity(_slice(forward, best, n), best[n], *bounds[n]) for n in names}
    return ClosestPoint(
        values={k: best[k] for k in names},
        fold_change=best_fc,
        binding_constraint=_binding_constraint(best, names, bounds, sens),
        note=("Closest achievable fold change anywhere inside the given bounds. Widen the bounds, "
              "free another parameter, or change the pull (cage-latch coupling) to go further."),
    )


def _solve_1d(forward, base, name, bounds, spec):
    lo, hi = bounds[name]
    xs = _axis(lo, hi, _SCAN_POINTS)
    f1 = _slice(forward, base, name)
    ys = np.array([f1(float(x)) for x in xs], dtype=float)
    roots = _bracket_roots(f1, xs, spec.center)
    monotonic = _is_monotonic(ys)

    if not roots:
        closest = _closest_from_samples(forward, [_vary(base, name, x) for x in xs],
                                       [name], bounds, spec)
        return dict(mode=MODE_ROOT if monotonic else MODE_ROOTS, feasible=False, closest=closest,
                    message=(f"Target fold change {spec.center:.4g} is not reachable for any "
                             f"{name} in [{lo:g}, {hi:g}]. Closest is "
                             f"{closest.fold_change:.4g} at {name} = "
                             f"{closest.values[name]:.6g}; binding constraint: "
                             f"{closest.binding_constraint}."))

    solutions = [_solution(forward, _vary(base, name, r), [name], bounds, spec) for r in roots]
    if monotonic:
        message = (f"fold_change is monotonic in {name} over [{lo:g}, {hi:g}]; one solution at "
                   f"{name} = {roots[0]:.6g}.")
        mode = MODE_ROOT
    else:
        message = (f"fold_change is NOT monotonic in {name} over [{lo:g}, {hi:g}]; "
                   f"{len(roots)} root(s) found by bracketing a {_SCAN_POINTS}-point scan. "
                   f"All are returned.")
        mode = MODE_ROOTS
    return dict(mode=mode, feasible=True, solutions=solutions, message=message)


def _feasible_mask(forward, base, xname, yname, xs, ys, spec) -> np.ndarray:
    return np.array([[spec.contains(forward({xname: float(x), yname: float(y)}))
                      for y in ys] for x in xs], dtype=bool)


def _most_robust(mask: np.ndarray, sampling: tuple[float, float]) -> tuple[int, int] | None:
    """Feasible cell furthest from any infeasible cell or bound.

    sampling is each axis' cell size measured as a fraction of that parameter's full
    bounds, so the distance is in normalised parameter space and the two axes compare
    fairly however the grid was cropped.
    """
    if not mask.any():
        return None
    padded = distance_transform_edt(np.pad(mask, 1, constant_values=False), sampling=sampling)
    dist = padded[1:-1, 1:-1]
    i, j = np.unravel_index(int(np.argmax(dist)), dist.shape)
    return int(i), int(j)


def _sampling(xs: np.ndarray, ys: np.ndarray, xb, yb) -> tuple[float, float]:
    return (abs(_to_unit(float(xs[1]), *xb) - _to_unit(float(xs[0]), *xb)),
            abs(_to_unit(float(ys[1]), *yb) - _to_unit(float(ys[0]), *yb)))


def _contour(forward, base, xname, yname, xbounds, ybounds, spec):
    """Target contour(s) as polylines, one bracketed solve per row."""
    xs = _axis(*xbounds, _CONTOUR_ROWS)
    ys = _axis(*ybounds, _CONTOUR_ROWS)
    branches: list[list[tuple[float, float]]] = []
    active: list[tuple[int, float]] = []

    for x in xs:
        row = dict(base)
        row[xname] = float(x)
        f1 = _slice(forward, row, yname)
        roots = _bracket_roots(f1, ys, spec.center)
        claimed: set[int] = set()
        nxt: list[tuple[int, float]] = []
        for r in roots:
            ru = _to_unit(r, *ybounds)
            best, best_d = None, None
            for k, (_, prev) in enumerate(active):
                if k in claimed:
                    continue
                d = abs(_to_unit(prev, *ybounds) - ru)
                if best_d is None or d < best_d:
                    best, best_d = k, d
            if best is not None and best_d <= _LINK_THRESHOLD:
                claimed.add(best)
                bi = active[best][0]
            else:
                branches.append([])
                bi = len(branches) - 1
            branches[bi].append((float(x), float(r)))
            nxt.append((bi, r))
        active = nxt
    return [b for b in branches if b]


def _margin(solution: Solution, names, bounds: dict) -> float:
    """Smallest drift any one parameter tolerates, as a fraction of its own bounds.

    This is the distance from the design to the region boundary, measured in
    normalised parameter space, and it is what "most robust" maximises.
    """
    if not solution.on_spec:
        return -1.0
    out = []
    for name in names:
        lo, hi = bounds[name]
        x = _to_unit(solution.values[name], lo, hi)
        t_lo, t_hi = solution.tolerance[name]
        out.append(min(x - _to_unit(t_lo, lo, hi), _to_unit(t_hi, lo, hi) - x))
    return min(out)


def _solve_2d(forward, base, names, bounds, spec):
    xname, yname = names
    xb, yb = bounds[xname], bounds[yname]
    region = _contour(forward, base, xname, yname, xb, yb, spec)

    # The feasible band is often thinner than one grid cell, so the region comes from
    # Brent-refined contour rows (grid-free in y) and the robust point is picked along it.
    # The coarse mask is only the fallback for a band the centre contour misses.
    candidates = [{**base, xname: x, yname: y} for line in region for x, y in line]
    if not candidates:
        xs, ys = _axis(*xb, _GRID_POINTS), _axis(*yb, _GRID_POINTS)
        mask = _feasible_mask(forward, base, xname, yname, xs, ys, spec)
        cell = _most_robust(mask, _sampling(xs, ys, xb, yb))
        if cell is None:
            samples = [{**base, xname: float(x), yname: float(y)} for x in xs for y in ys]
            closest = _closest_from_samples(forward, samples, names, bounds, spec)
            return dict(mode=MODE_REGION, feasible=False, closest=closest,
                        message=(f"No ({xname}, {yname}) pair inside the given bounds reaches fold "
                                 f"change {spec.center:.4g}. Closest is {closest.fold_change:.4g} at "
                                 f"{closest.values}; binding constraint: "
                                 f"{closest.binding_constraint}."))
        candidates = [{**base, xname: float(xs[cell[0]]), yname: float(ys[cell[1]])}]

    solutions = [_solution(forward, values, names, bounds, spec) for values in candidates]
    robust = max(solutions, key=lambda s: _margin(s, names, bounds))
    return dict(mode=MODE_REGION, feasible=True, solutions=[robust], robust_point=robust,
                region=region, region_axes=(xname, yname),
                message=(f"Feasible region returned as {len(region)} contour polyline(s) in "
                         f"({xname}, {yname}); the representative design is the feasible point "
                         f"furthest from the region boundary, so it is the least sensitive to "
                         f"parameter error."))


def _solve_nd(forward, base, names, bounds, spec):
    """Robust point from multi-start constrained optimisation, plus per-parameter drift."""
    def unit_to_values(u) -> dict:
        return {**base, **{n: _from_unit(float(ui), *bounds[n]) for n, ui in zip(names, u)}}

    def objective(u) -> float:
        fc = forward(unit_to_values(u))
        if not math.isfinite(fc) or fc <= 0:
            return 1e6
        return (math.log(fc) - math.log(spec.center)) ** 2

    def grad_norm(values: dict) -> float:
        # Sensitivity in unit coords, so parameters in different units compare fairly.
        total = 0.0
        for n in names:
            lo, hi = bounds[n]
            s = _sensitivity(_slice(forward, values, n), values[n], lo, hi)
            if not math.isfinite(s):
                return float("inf")
            total += (s * (hi - lo)) ** 2
        return math.sqrt(total)

    rng = np.random.default_rng(_ND_SEED)
    starts = [np.full(len(names), 0.5)] + [rng.random(len(names)) for _ in range(_ND_STARTS - 1)]
    on_spec, all_values = [], []
    for u0 in starts:
        res = minimize(objective, u0, method="L-BFGS-B", bounds=[(0.0, 1.0)] * len(names))
        values = unit_to_values(res.x)
        all_values.append(values)
        if spec.contains(forward(values)):
            on_spec.append(values)

    if not on_spec:
        grid = [{**base, **{n: _from_unit(float(v), *bounds[n]) for n, v in zip(names, u)}}
                for u in rng.random((400, len(names)))]
        closest = _closest_from_samples(forward, all_values + grid, names, bounds, spec)
        return dict(mode=MODE_ROBUST_POINT, feasible=False, closest=closest,
                    message=(f"No combination of {', '.join(names)} inside the given bounds reaches "
                             f"fold change {spec.center:.4g}. Closest is {closest.fold_change:.4g}; "
                             f"binding constraint: {closest.binding_constraint}."))

    best = min(on_spec, key=grad_norm)
    robust = _solution(forward, best, names, bounds, spec)
    return dict(mode=MODE_ROBUST_POINT, feasible=True, solutions=[robust], robust_point=robust,
                message=(f"{len(names)} free parameters: no feasible region is returned (it is neither "
                         f"drawable nor cheap at this dimension). This is the on-spec point with the "
                         f"flattest normalised sensitivity out of {len(starts)} starts, with the drift "
                         f"each parameter tolerates on its own."))


def inverse_design(solve_for, fixed: dict, bounds: dict, target_fold_change: float | None = None,
                   target_fold_change_range: tuple[float, float] | None = None,
                   temperature_K: float = T_REF_DEFAULT,
                   specs: dict[str, ThermoSpec] | None = None,
                   rel_tol: float = 0.05) -> InverseResult:
    names = list(solve_for)
    if not names:
        raise ValueError("solve_for must name at least one free parameter")
    if len(set(names)) != len(names):
        raise ValueError("solve_for has duplicates")
    unknown = [n for n in names if n not in PARAMS]
    if unknown:
        raise ValueError(f"cannot solve for {unknown}; choose from {list(PARAMS)}")
    overlap = [n for n in names if n in fixed]
    if overlap:
        raise ValueError(f"{overlap} appear in both solve_for and fixed")
    if "pull" not in fixed:
        raise ValueError("fixed must include pull (the cage-latch allosteric coupling)")
    missing_bounds = [n for n in names if n not in bounds]
    if missing_bounds:
        raise ValueError(f"no bounds given for {missing_bounds}")
    for n in names:
        lo, hi = bounds[n]
        if hi <= lo:
            raise ValueError(f"bounds for {n} must satisfy min < max")
        if n in CONC_PARAMS and lo <= 0:
            raise ValueError(f"bounds for {n} must be > 0")

    spec = TargetSpec.build(target_fold_change, target_fold_change_range, rel_tol)
    specs = specs or {}

    base = {k: v for k, v in fixed.items() if k in PARAMS}
    base.setdefault("target_conc", None)
    saturating = base.get("target_conc") is None and "target_conc" not in names
    if saturating:
        base["target_conc"] = None
        base.setdefault("dG_target_binder", None)
    needed = [p for p in PARAMS if p not in ("dG_target_binder", "target_conc")]
    if not saturating:
        needed += ["dG_target_binder", "target_conc"]
    absent = [p for p in needed if p not in names and base.get(p) is None]
    if absent:
        raise ValueError(f"fixed is missing {absent} (or free them in solve_for). Leave "
                         f"dG_target_binder and target_conc out together to assume saturating target")
    for n in names:
        base[n] = float(np.mean(bounds[n]))

    forward = make_forward(base, float(fixed["pull"]), temperature_K, specs)
    warnings: list[str] = []
    for key in DG_PARAMS:
        if base.get(key) is None and key not in names:
            continue
        warnings += warnings_at(specs.get(key, ThermoSpec(label=key)), temperature_K)
    mode = weakest_mode([specs.get(k, ThermoSpec()).mode for k in DG_PARAMS
                         if k in names or base.get(k) is not None])

    if len(names) == 1:
        out = _solve_1d(forward, base, names[0], bounds, spec)
    elif len(names) == 2:
        out = _solve_2d(forward, base, names, bounds, spec)
    else:
        out = _solve_nd(forward, base, names, bounds, spec)

    return InverseResult(temperature_K=temperature_K, thermo_mode=mode, warnings=warnings, **out)
