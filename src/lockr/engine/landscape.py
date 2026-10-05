"""Two-parameter metric grids for Module D."""

from __future__ import annotations

import math

import numpy as np

from . import thermo
from .inverse import _most_robust as most_robust_cell
from .models import SensorParams
from .temperature import R_KCAL, T_REF_DEFAULT, ThermoSpec, shift_constant, warnings_at

# Anything that can be an axis or held fixed. Mirrors sweep's sweepable set plus temperature.
AXIS_PARAMS = ("k_ck", "k_open", "luckey", "pull", "temperature_K")
REQUIRED_PARAMS = ("k_ck", "k_open", "luckey", "pull")
OPTIONAL_PARAMS = ("k_target", "target_conc")

METRICS = ("fold_change", "on_signal", "off_signal", "dynamic_range",
           "limit_of_detection", "fraction_open")

# Metrics needing a target binder: without a Kd there is no dose response to read an LOD off.
NEEDS_TARGET = ("limit_of_detection",)

SCALES = ("linear", "log")

# memory guard, not a time guard
MAX_CELLS = 40_000
MAX_CELLS_LOD = 10_000
_LOD_POINTS = 500  # same count and window as thermo.lod_and_ec50
_LOD_CONC = (-15, -5)
_LOD_THRESHOLD = 2.0

_TM_MARGIN_K = 10.0  # kept in step with temperature._TM_MARGIN_K

_OPS = {">=": np.greater_equal, "<=": np.less_equal, ">": np.greater, "<": np.less}


def axis_values(lo: float, hi: float, steps: int, scale: str) -> np.ndarray:
    """Grid positions along one axis. Uniform in whatever space `scale` names."""
    if scale not in SCALES:
        raise ValueError(f"scale must be one of {SCALES}")
    if hi <= lo:
        raise ValueError("axis min must be < max")
    if steps < 2:
        raise ValueError("axis steps must be >= 2")
    if scale == "log":
        if lo <= 0:
            raise ValueError("log axis min must be > 0")
        return np.logspace(math.log10(lo), math.log10(hi), steps)
    return np.linspace(lo, hi, steps)


def cap_steps(x_steps: int, y_steps: int, metric: str) -> tuple[int, int, list[str]]:
    """Shrink an oversized grid proportionally, and say so."""
    cap = MAX_CELLS_LOD if metric in NEEDS_TARGET else MAX_CELLS
    if x_steps * y_steps <= cap:
        return x_steps, y_steps, []
    shrink = math.sqrt(cap / (x_steps * y_steps))
    nx, ny = max(2, int(x_steps * shrink)), max(2, int(y_steps * shrink))
    return nx, ny, [
        f"Requested {x_steps}x{y_steps} = {x_steps * y_steps} cells, above the "
        f"{cap} cell cap for metric '{metric}'. Reduced to {nx}x{ny}. "
        f"Narrow the axis ranges for finer detail over a smaller region."]


def _shift_factor(T, spec: ThermoSpec) -> tuple[np.ndarray, list[str]]:
    """Temperature multiplier for an equilibrium constant."""
    # multiplicative in k_ref, so a probe at k=1 is the factor
    Ts = np.atleast_1d(np.asarray(T, dtype=float))
    factors = np.empty(Ts.shape)
    for idx, t in np.ndenumerate(Ts):
        factors[idx], _ = shift_constant(1.0, float(t), spec)
    # all three warning kinds are monotone in T, so the ends suffice
    warns = []
    for t in dict.fromkeys((float(Ts.flat[0]), float(Ts.flat[-1]))):
        warns += warnings_at(spec, t)
    shaped = factors.reshape(np.shape(T)) if np.ndim(T) else factors.reshape(())
    return shaped, list(dict.fromkeys(warns))


def _open_fraction(k_open, ratio):
    """Fraction of cages open, bound or not."""
    open_weight = k_open * (1.0 + ratio)
    return open_weight / (1.0 + open_weight)


def evaluate(values: dict, metric: str, specs: dict[str, ThermoSpec]) -> tuple[np.ndarray, list[str]]:
    """One metric over a parameter grid."""
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}")

    T = values.get("temperature_K", T_REF_DEFAULT)
    latch_factor, w_latch = _shift_factor(T, specs["dG_latch_cage"])
    key_factor, w_key = _shift_factor(T, specs["dG_key_cage"])
    warns = w_latch + w_key

    k_open = np.asarray(values["k_open"], dtype=float) * latch_factor
    k_ck = np.asarray(values["k_ck"], dtype=float) * key_factor
    pull = np.asarray(values["pull"], dtype=float)
    params = SensorParams(K_open=k_open, K_CK=k_ck,
                          lucKey=np.asarray(values["luckey"], dtype=float),
                          RT=R_KCAL * np.asarray(T, dtype=float))

    target_conc = values.get("target_conc")
    k_target = values.get("k_target")
    if metric in NEEDS_TARGET and k_target is None:
        raise ValueError(f"metric '{metric}' needs k_target (the target binder Kd)")

    Kd = None
    if k_target is not None:
        kd_factor, w_target = _shift_factor(T, specs["dG_target_binder"])
        Kd = np.asarray(k_target, dtype=float) * kd_factor
        warns += w_target

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        if metric == "limit_of_detection":
            grid = _lod(k_open, pull, params, Kd)
        else:
            # theta = 1 means saturating target, matching /sweep's Kd=1.0 convention.
            theta = 1.0 if (target_conc is None or Kd is None) else thermo.theta(
                np.asarray(target_conc, dtype=float), Kd)
            k_open_eff = thermo.k_open_eff(params.K_open, pull, theta)
            # one definition of the partition function, not a driftable copy
            on = thermo._f_open(k_open_eff, params)
            off = thermo._f_open(params.K_open, params)
            grid = {
                "fold_change": lambda: on / off,
                "on_signal": lambda: on,
                "off_signal": lambda: off,
                "dynamic_range": lambda: on - off,
                "fraction_open": lambda: _open_fraction(k_open_eff, params.luckey_ratio),
            }[metric]()

    grid = np.broadcast_to(np.asarray(grid, dtype=float), np.broadcast_shapes(
        np.shape(grid), _grid_shape(values))).copy()
    grid[~np.isfinite(grid)] = np.nan
    return grid, warns


def _grid_shape(values: dict) -> tuple[int, ...]:
    shapes = [np.shape(v) for v in values.values() if v is not None]
    return np.broadcast_shapes(*shapes) if shapes else ()


def _lod(k_open, pull, params: SensorParams, Kd) -> np.ndarray:
    """2x crossing concentration per cell."""
    conc = np.logspace(*_LOD_CONC, _LOD_POINTS)
    expand = lambda a: np.asarray(a, dtype=float)[..., None]
    wide = SensorParams(K_open=expand(params.K_open), K_CK=expand(params.K_CK),
                        lucKey=expand(params.lucKey), RT=expand(params.RT))
    theta = thermo.theta(conc, expand(Kd))
    on = thermo._f_open(thermo.k_open_eff(wide.K_open, expand(pull), theta), wide)
    fcs = on / thermo._f_open(wide.K_open, wide)

    crossed = np.nan_to_num(fcs, nan=0.0) >= _LOD_THRESHOLD
    first = conc[np.argmax(crossed, axis=-1)]
    return np.where(crossed.any(axis=-1), first, np.nan)


def feasibility(grids: dict[str, np.ndarray], constraints) -> np.ndarray | None:
    """Cells where every constraint holds. nan fails every constraint."""
    if not constraints:
        return None
    mask = np.ones(next(iter(grids.values())).shape, dtype=bool)
    for c in constraints:
        grid = grids[c["metric"]]
        mask &= np.where(np.isnan(grid), False, _OPS[c["op"]](grid, c["value"]))
    return mask


def _axis_slack(mask: np.ndarray, i: int, j: int, values: np.ndarray,
                along_x: bool) -> tuple[float, float]:
    """Drift allowed along one axis, native units."""
    # ponytail: grid-resolution-limited, bisect the metric if a tighter number is needed
    line = mask[:, j] if along_x else mask[i, :]
    here = i if along_x else j
    lo = hi = here
    while lo > 0 and line[lo - 1]:
        lo -= 1
    while hi < len(line) - 1 and line[hi + 1]:
        hi += 1
    return float(values[here] - values[lo]), float(values[hi] - values[here])


def robust_optimum(mask: np.ndarray | None, xs: np.ndarray, ys: np.ndarray,
                   grid: np.ndarray) -> dict | None:
    """Feasible cell furthest from the boundary."""
    # inverse._most_robust, so Module C and D agree on "robust"
    if mask is None or not mask.any():
        return None
    cell = most_robust_cell(mask, (1.0 / (len(xs) - 1), 1.0 / (len(ys) - 1)))
    if cell is None:
        return None
    i, j = cell
    x_minus, x_plus = _axis_slack(mask, i, j, xs, along_x=True)
    y_minus, y_plus = _axis_slack(mask, i, j, ys, along_x=False)
    value = grid[i, j]
    return {
        "x": float(xs[i]), "y": float(ys[j]),
        "metric_value": None if np.isnan(value) else float(value),
        "x_slack": {"minus": x_minus, "plus": x_plus},
        "y_slack": {"minus": y_minus, "plus": y_plus},
    }


def contour_levels(grid: np.ndarray, n: int = 5) -> list[float]:
    """Round-numbered levels spanning the finite values present."""
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        return []
    lo, hi = float(finite.min()), float(finite.max())
    if not hi > lo:
        return []
    span = hi - lo
    step = 10.0 ** math.floor(math.log10(span / n))
    for mult in (1, 2, 2.5, 5, 10, 20, 50):
        if span / (step * mult) <= n:
            step *= mult
            break
    first = math.ceil(lo / step) * step
    return [round(first + k * step, 12) for k in range(n + 1)
            if lo < first + k * step < hi]


def contours(grid: np.ndarray, xs: np.ndarray, ys: np.ndarray,
             levels) -> list[dict]:
    """Marching-squares polylines on the computed grid."""
    # traced on the computed grid, so a contour never disagrees with the heatmap
    out = []
    for level in levels:
        segments = _segments(grid, xs, ys, float(level))
        if segments:
            out.append({"level": float(level), "polylines": _link(segments)})
    return out


def _interp(v0: float, v1: float, a0: float, a1: float, level: float) -> float:
    return a0 if v1 == v0 else a0 + (a1 - a0) * (level - v0) / (v1 - v0)


def _segments(grid: np.ndarray, xs: np.ndarray, ys: np.ndarray,
              level: float) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """One marching-squares pass."""
    segments = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            corners = (grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1])
            if not all(math.isfinite(c) for c in corners):
                continue
            below = [c < level for c in corners]
            if all(below) or not any(below):
                continue
            # Crossings on the four edges, in corner order: bottom, right, top, left.
            edges = []
            pts = ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
            for k in range(4):
                a, b = k, (k + 1) % 4
                if below[a] != below[b]:
                    (ia, ja), (ib, jb) = pts[a], pts[b]
                    edges.append((
                        _interp(corners[a], corners[b], float(xs[ia]), float(xs[ib]), level),
                        _interp(corners[a], corners[b], float(ys[ja]), float(ys[jb]), level)))
            # four crossings is a saddle, joined in corner order
            for k in range(0, len(edges) - 1, 2):
                segments.append((edges[k], edges[k + 1]))
    return segments


def _link(segments) -> list[list[list[float]]]:
    """Chain segments end to end so each polyline closes or ends on the grid edge."""
    kk = lambda p: (round(p[0], 12), round(p[1], 12))
    adjacent: dict[tuple, list[int]] = {}
    for idx, (a, b) in enumerate(segments):
        adjacent.setdefault(kk(a), []).append(idx)
        adjacent.setdefault(kk(b), []).append(idx)

    used, polylines = set(), []
    for start, (a, b) in enumerate(segments):
        if start in used:
            continue
        used.add(start)
        chain = [a, b]
        # Extend the tail, flip, extend what was the head, flip back.
        for _ in range(2):
            while True:
                tip = kk(chain[-1])
                nxt = next((i for i in adjacent.get(tip, []) if i not in used), None)
                if nxt is None:
                    break
                used.add(nxt)
                s0, s1 = segments[nxt]
                chain.append(s1 if kk(s0) == tip else s0)
            chain.reverse()
        polylines.append([[float(x), float(y)] for x, y in chain])
    return polylines


def tm_band(spec: ThermoSpec, lo: float, hi: float) -> dict | None:
    """Temperature span where the two-state model stops holding, if it is in range."""
    if spec.Tm is None:
        return None
    start = spec.Tm - _TM_MARGIN_K
    if start > hi:
        return None
    return {"min": float(max(start, lo)), "max": float(hi), "Tm": float(spec.Tm),
            "label": f"within {_TM_MARGIN_K:.0f} K of Tm = {spec.Tm:.1f} K: "
                     f"cage unfolds, two-state model invalid"}


def temperature_axis_warnings(specs: dict[str, ThermoSpec], lo: float, hi: float) -> list[str]:
    """Temperature-axis warnings, raised at both ends."""
    seen = []
    for spec in specs.values():
        for T in (lo, hi):
            seen += warnings_at(spec, T)
    return list(dict.fromkeys(seen))
