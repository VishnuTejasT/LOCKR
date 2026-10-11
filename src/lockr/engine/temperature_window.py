"""Temperature operating window: where on the thermometer a construct actually works.

Module C asks "what fold change do I get at temperature T". This asks the inverse
question a wet-lab team cares about: over a temperature sweep, where does fold change
PEAK, and over what band does the design clear a fold-change bar often enough to be
worth building. "Often enough" is a Monte Carlo over the inputs nobody has measured yet.

Three interactions carry the construct: Latch-Cage (K_open), Binder-Target (Kd) and
Key-Cage (K_CK). Per interaction the caller supplies dG at T_ref; dH is supplied or
drawn from a per-interaction-type prior; dS is then pinned by dS = (dH - dG)/T_ref and
dG(T) = dH - T*dS.

dCp IS IGNORED here. That is the van't Hoff approximation and it only holds within
roughly 15 K of T_ref (see temperature._VANT_HOFF_WINDOW_K). A window reported 30 K
from T_ref is an extrapolation, not a measurement. temperature.dG_at carries the full
Gibbs-Helmholtz form with a heat-capacity term if you need it.

Fold change itself is NOT recomputed here. It comes from thermo.k_open_eff / thermo.theta
/ thermo._f_open, the same three primitives landscape.evaluate calls, so there is one
partition function in this codebase and not a drifting copy.

Boundary convention, verified against thermo.fold_change rather than assumed:
    zero target  -> fold change is exactly 1.0 (baseline, "no change"), NOT 0.
    zero Key     -> fold change is undefined (0/0) and comes back nan, NOT 0 and NOT 1.
nan is a failure for p_works, never a pass, matching landscape.feasibility.
"""

from __future__ import annotations

import numpy as np

from . import thermo
from .models import SensorParams
from .temperature import R_KCAL, T_REF_DEFAULT

# ASSUMED dH priors, kcal/mol. No published dH exists for lucCage or LOCKR switches (the
# 2021 Nature paper and its SARS-CoV-2 follow-up report none), so these are NOT measured and
# NOT fitted. They are negative (binding releases heat) and sit well inside the range of
# 196 protein-protein complexes measured by ITC, dH from -66.7 to +19.9 kcal/mol at 298 K
# (Teilum, Olsen & Kragelund 2015, https://pmc.ncbi.nlm.nih.gov/articles/PMC4496568/). That
# paper reports no mean or SD for dH, so the spreads here are a judgement call. Replace with
# calorimetry or van't Hoff fits for your own construct. This is the ONE place they are defined.
DEFAULT_DH_PRIORS = {
    "latch_cage": {"mean": -12.0, "sd": 4.0},
    "binder_target": {"mean": -10.0, "sd": 4.0},
    "key_cage": {"mean": -8.0, "sd": 3.0},
}

PRIOR_KEYS = tuple(DEFAULT_DH_PRIORS)

DH_PRIOR_PLACEHOLDER_NOTE = (
    "The dH priors are ASSUMED, not measured: no dH has been published for lucCage. They are "
    "negative and inside the -67 to +20 kcal/mol range seen across 196 protein-protein "
    "complexes (Teilum 2015), but that survey gives no mean or spread to fit. Every "
    "temperature in this report depends on them. Enter your own dH to replace them."
)

# Docking-derived affinities are noisy, so the binder dG is sampled rather than trusted.
# Narrow this for an ITC-measured affinity; widen it for a docking score.
DEFAULT_DG_BINDER_SD = 1.5

# Allosteric coupling is not a temperature question, so it is not swept. 10x is the
# value LATTE's own examples and thermo.diagnose_regime use.
DEFAULT_PULL = 10.0

DEFAULT_N_SAMPLES = 2000
MIN_N_SAMPLES = 200
MAX_N_SAMPLES = 20_000
MAX_GRID_POINTS = 500

DEFAULT_FC_THRESHOLD = 2.0
DEFAULT_CONFIDENCE = 0.8

ZERO_C_IN_K = 273.15


def c_to_k(t_c: float) -> float:
    return t_c + ZERO_C_IN_K


def dG_from_kd_nM(kd_nM: float, T_ref: float) -> float:
    """Kd in nM -> binding dG in kcal/mol at T_ref. Same sign convention as thermo."""
    if kd_nM <= 0:
        raise ValueError("kd_nM must be > 0")
    return R_KCAL * T_ref * np.log(kd_nM * 1e-9)


def resolve_priors(overrides: dict | None = None) -> dict[str, dict[str, float]]:
    priors = {k: dict(v) for k, v in DEFAULT_DH_PRIORS.items()}
    for key, entry in (overrides or {}).items():
        if key not in DEFAULT_DH_PRIORS:
            raise ValueError(f"unknown dH prior {key!r}; use {PRIOR_KEYS}")
        if entry is None:
            continue
        for field in ("mean", "sd"):
            value = entry.get(field) if isinstance(entry, dict) else getattr(entry, field, None)
            if value is None:
                continue
            if field == "sd" and value < 0:
                raise ValueError(f"dH prior {key}.sd must be >= 0")
            priors[key][field] = float(value)
    return priors


def temperature_grid(t_min_c: float, t_max_c: float, t_step_c: float) -> np.ndarray:
    """Inclusive sweep grid in Celsius."""
    if t_step_c <= 0:
        raise ValueError("t_step_c must be > 0")
    if t_min_c >= t_max_c:
        raise ValueError("t_min_c must be < t_max_c")
    # half-step pad so an exact t_max lands on the grid despite float drift
    grid = np.arange(t_min_c, t_max_c + t_step_c / 2.0, t_step_c, dtype=float)
    if grid.size > MAX_GRID_POINTS:
        raise ValueError(
            f"sweep has {grid.size} grid points, above the {MAX_GRID_POINTS} point cap; "
            f"widen t_step_c or narrow the range")
    return grid


def dG_at_T(dG_ref, dH, T, T_ref: float = T_REF_DEFAULT):
    """van't Hoff dG(T), vectorised. dCp is ignored.

    dS = (dH - dG_ref)/T_ref pins the entropy at T_ref, then dG(T) = dH - T*dS. This is
    temperature.dG_at with dCp = 0; that scalar function stays the reference definition
    and test_temperature_window asserts the two agree.
    """
    dS = (dH - dG_ref) / T_ref
    return dH - T * dS


def _fold_change(k_open, k_ck, key_conc, pull, Kd, target_conc, T):
    """Fold change from the existing thermo primitives. nan where undefined."""
    params = SensorParams(K_open=k_open, K_CK=k_ck, lucKey=key_conc, RT=R_KCAL * T)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        # theta = 1 (saturating target) when there is no target to titrate, matching
        # landscape.evaluate and /sweep's Kd=1.0 convention.
        theta = 1.0 if (target_conc is None or Kd is None) else thermo.theta(target_conc, Kd)
        on = thermo._f_open(thermo.k_open_eff(params.K_open, pull, theta), params)
        off = thermo._f_open(params.K_open, params)
        fc = np.asarray(on / off, dtype=float)
    fc[~np.isfinite(fc)] = np.nan
    return fc


def _percentile(curves: np.ndarray, q: float) -> np.ndarray:
    """Percentile down the sample axis; an all-nan temperature gives nan, not a warning."""
    out = np.full(curves.shape[1], np.nan)
    usable = np.isfinite(curves).any(axis=0)
    if usable.any():
        out[usable] = np.nanpercentile(curves[:, usable], q, axis=0)
    return out


def _segments(flags: np.ndarray) -> list[tuple[int, int]]:
    """Index spans of contiguous True runs."""
    runs, start = [], None
    for i, ok in enumerate(flags):
        if ok and start is None:
            start = i
        elif not ok and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(flags) - 1))
    return runs


def operating_window(temps_c: np.ndarray, p_works: np.ndarray,
                     confidence: float = DEFAULT_CONFIDENCE) -> dict:
    """Longest contiguous band where p_works clears `confidence`.

    Returns the longest segment as `window`, any shorter qualifying bands as
    `other_segments`, and a `reason` when nothing qualifies. A single qualifying grid
    point is a real (zero-width) window, not a refusal.
    """
    temps_c = np.asarray(temps_c, dtype=float)
    p_works = np.asarray(p_works, dtype=float)
    if temps_c.shape != p_works.shape:
        raise ValueError("temps_c and p_works must be the same length")

    qualifying = np.nan_to_num(p_works, nan=0.0) >= confidence
    runs = _segments(qualifying)
    if not runs:
        best = float(np.nanmax(p_works)) if np.isfinite(p_works).any() else 0.0
        where = (f" (best {best:.2f} at {temps_c[int(np.nanargmax(p_works))]:.1f} C)"
                 if np.isfinite(p_works).any() else "")
        return {
            "window": None,
            "other_segments": [],
            "window_touches_grid_edge": False,
            "reason": (f"No temperature on the swept grid reaches p_works >= {confidence:g}"
                       f"{where}. Lower the fold-change threshold, lower the confidence, "
                       f"widen the sweep, or improve the design."),
        }

    # widest band wins; a tie goes to the lower temperature so the answer is stable
    best = min(runs, key=lambda r: (-(r[1] - r[0]), r[0]))
    lo, hi = best
    span = lambda r: {"t_low_c": float(temps_c[r[0]]), "t_high_c": float(temps_c[r[1]])}
    return {
        "window": span(best),
        # the rest stay in temperature order, which is the order a reader scans them in
        "other_segments": [span(r) for r in runs if r != best],
        "window_touches_grid_edge": bool(lo == 0 or hi == len(temps_c) - 1),
        "reason": None,
    }


def _stability_cap(tm_c: float | None, t_max_stability_c: float | None) -> tuple[float | None, str | None]:
    """Lowest temperature at which the construct is still folded, and who set it."""
    candidates = [(v, src) for v, src in ((tm_c, "binder"), (t_max_stability_c, "system"))
                  if v is not None]
    if not candidates:
        return None, None
    # a tie goes to the binder: the grafted domain is the part that was actually measured
    value, source = min(candidates, key=lambda c: (c[0], c[1] != "binder"))
    return float(value), source


def _activity_factor(temps_c: np.ndarray, points) -> np.ndarray:
    """Relative luciferase activity at each swept temperature, linearly interpolated.

    No default curve exists, so this is only called when the caller supplies one.
    np.interp holds the end values flat outside the supplied range.
    """
    pairs = sorted((float(t), float(a)) for t, a in points)
    xs = np.array([p[0] for p in pairs])
    ys = np.array([p[1] for p in pairs])
    return np.interp(temps_c, xs, ys)


def scan_system(system: dict, temps_c: np.ndarray, seed,
                t_ref_c: float = 25.0, fc_threshold: float = DEFAULT_FC_THRESHOLD,
                confidence: float = DEFAULT_CONFIDENCE,
                n_samples: int = DEFAULT_N_SAMPLES, dH_priors: dict | None = None,
                luciferase_activity=None) -> dict:
    """One construct's temperature response, fold-change band and operating window.

    `system` carries dG_latch_cage, dG_key_cage, optional dH for each, optional
    key_conc/target_conc (molar), optional t_max_stability_c, and a `binder` dict.
    Fold change comes from thermo; dCp is ignored throughout.
    """
    T_ref = c_to_k(t_ref_c)
    priors = resolve_priors(dH_priors)
    binder = system["binder"]
    temps_c = np.asarray(temps_c, dtype=float)
    T = c_to_k(temps_c)[None, :]                      # (1, n_T)
    rng = np.random.default_rng(seed)

    # One draw order, fixed, so a seed reproduces a run exactly.
    def draw(value: float | None, prior_key: str) -> np.ndarray:
        if value is not None:
            return np.full(n_samples, float(value))
        p = priors[prior_key]
        return rng.normal(p["mean"], p["sd"], n_samples)

    dH_latch = draw(system.get("dH_latch_cage"), "latch_cage")
    dH_key = draw(system.get("dH_key_cage"), "key_cage")
    dH_binder = draw(binder.get("dH_target_binder"), "binder_target")

    graft_shift = float(binder.get("graft_dG_latch_cage_shift") or 0.0)
    dG_latch_raw = float(system["dG_latch_cage"])
    dG_latch_eff = dG_latch_raw + graft_shift

    dG_binder_mean = float(binder["dG_target_binder"])
    dG_binder_sd = binder.get("dG_target_binder_sd")
    dG_binder_sd = DEFAULT_DG_BINDER_SD if dG_binder_sd is None else float(dG_binder_sd)
    dG_binder = rng.normal(dG_binder_mean, dG_binder_sd, n_samples)

    col = lambda a: np.asarray(a, dtype=float)[:, None]               # (n_samples, 1)
    dG_latch_T = dG_at_T(dG_latch_eff, col(dH_latch), T, T_ref)
    dG_key_T = dG_at_T(float(system["dG_key_cage"]), col(dH_key), T, T_ref)
    dG_binder_T = dG_at_T(col(dG_binder), col(dH_binder), T, T_ref)

    with np.errstate(over="ignore", invalid="ignore"):
        RT = R_KCAL * T
        k_open = np.exp(dG_latch_T / RT)
        k_ck = np.exp(dG_key_T / RT)
        Kd = np.exp(dG_binder_T / RT)

    target_conc = system.get("target_conc")
    key_conc = system.get("key_conc")
    curves = _fold_change(
        k_open, k_ck,
        SensorParams().lucKey if key_conc is None else float(key_conc),
        float(system.get("pull") or DEFAULT_PULL),
        Kd, None if target_conc is None else float(target_conc), T)

    activity = None
    if luciferase_activity:
        activity = _activity_factor(temps_c, luciferase_activity)
        # Scales the reported signal ratio, so a dead enzyme reads as no switching.
        curves = curves * activity[None, :]

    cap_c, cap_source = _stability_cap(binder.get("tm_c"), system.get("t_max_stability_c"))
    below_cap = np.ones(temps_c.shape, dtype=bool) if cap_c is None else temps_c <= cap_c

    # nan is a failure, never a pass: landscape.feasibility treats it the same way.
    p_works = (np.nan_to_num(curves, nan=0.0) >= fc_threshold).mean(axis=0)
    p_works = np.where(below_cap, p_works, 0.0)

    median_fc = _percentile(curves, 50)
    p10_fc = _percentile(curves, 10)
    p90_fc = _percentile(curves, 90)

    # T_opt only ranks temperatures where the construct is still folded; a peak above
    # the stability cap is an artefact of a two-state model applied to a denatured protein.
    rankable = np.where(below_cap & np.isfinite(median_fc), median_fc, -np.inf)
    t_opt_c = float(temps_c[int(np.argmax(rankable))]) if np.isfinite(rankable).any() else None
    peak_median_fc = float(np.max(rankable)) if t_opt_c is not None else None

    window = operating_window(temps_c, p_works, confidence)
    nullable = lambda a: [None if not np.isfinite(v) else float(v) for v in a]

    return {
        "name": system.get("name") or "system",
        "binder_name": binder.get("name") or "binder",
        "t_opt_c": t_opt_c,
        "peak_median_fc": peak_median_fc,
        "median_fc": nullable(median_fc),
        "p10_fc": nullable(p10_fc),
        "p90_fc": nullable(p90_fc),
        "p_works": [float(v) for v in p_works],
        "dG_latch_cage_raw": dG_latch_raw,
        "dG_latch_cage_effective": dG_latch_eff,
        "graft_dG_latch_cage_shift": graft_shift,
        "dG_target_binder": dG_binder_mean,
        "dG_target_binder_sd": dG_binder_sd,
        "stability_cap_c": cap_c,
        "stability_cap_source": cap_source,
        "binder_affects_result": target_conc is not None,
        "luciferase_factor": None if activity is None else [float(v) for v in activity],
        **window,
    }


def _notes(results: list[dict], dH_priors: dict | None,
           luciferase_activity, t_ref_c: float) -> list[str]:
    notes = [
        f"dCp is ignored: dG(T) uses van't Hoff from T_ref = {t_ref_c:g} C. Trust the curve "
        f"within ~15 K of T_ref and treat the rest as extrapolation.",
        "Operating windows depend entirely on the dH priors. A different prior moves the "
        "window, so the width is not evidence the design is robust.",
    ]
    if dH_priors is None:
        notes.append(DH_PRIOR_PLACEHOLDER_NOTE)
    else:
        supplied = sorted(k for k, v in dH_priors.items() if v is not None)
        missing = [k for k in PRIOR_KEYS if k not in supplied]
        if missing:
            notes.append(f"{DH_PRIOR_PLACEHOLDER_NOTE} Still assumed: "
                         f"{', '.join(missing)}.")
    for r in results:
        notes.append(f"{r['name']} / {r['binder_name']}: Binder-Target dG uncertainty sampled "
                     f"with sd {r['dG_target_binder_sd']:g} kcal/mol; graft shift on Latch-Cage "
                     f"dG assumed {r['graft_dG_latch_cage_shift']:+g} kcal/mol.")
        if not r["binder_affects_result"]:
            # theta is pinned at 1 with no target_conc, so Kd cancels out of the ratio
            # entirely. Saying nothing here would let a sampled sd look like it did work.
            notes.append(f"{r['name']} / {r['binder_name']}: no target concentration was given, "
                         f"so the target is treated as saturating and the Binder-Target affinity "
                         f"(and its uncertainty) does NOT affect this result. Supply target_conc "
                         f"to make the binder matter.")
        if r["stability_cap_source"]:
            notes.append(f"{r['name']} / {r['binder_name']}: p_works is forced to 0 above "
                         f"{r['stability_cap_c']:g} C, the {r['stability_cap_source']} stability "
                         f"cap. Fold change is still drawn above it but is not meaningful there.")
        if r["window_touches_grid_edge"]:
            notes.append(f"{r['name']} / {r['binder_name']}: the window runs into the edge of "
                         f"the swept range, so its true extent is unknown. Widen the sweep.")
    if luciferase_activity:
        notes.append("Luciferase activity was supplied by you and multiplies the reported signal "
                     "ratio; it is interpolated linearly and held flat outside your points.")
    else:
        notes.append("Luciferase temperature activity is OFF. No default curve is assumed, so "
                     "the curve is switching behaviour only, not enzyme output.")
    return notes


def scan(systems: list[dict], t_min_c: float = 4.0, t_max_c: float = 45.0,
         t_step_c: float = 1.0, t_ref_c: float = 25.0,
         fc_threshold: float = DEFAULT_FC_THRESHOLD, confidence: float = DEFAULT_CONFIDENCE,
         n_samples: int = DEFAULT_N_SAMPLES, dH_priors: dict | None = None,
         seed: int | None = None, luciferase_activity=None) -> dict:
    """Temperature window for every construct, over one shared sweep grid."""
    if not 1 <= len(systems) <= 10:
        raise ValueError("systems must hold between 1 and 10 constructs")
    if not MIN_N_SAMPLES <= n_samples <= MAX_N_SAMPLES:
        raise ValueError(f"n_samples must be between {MIN_N_SAMPLES} and {MAX_N_SAMPLES}")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must be between 0 and 1, exclusive")

    temps_c = temperature_grid(t_min_c, t_max_c, t_step_c)
    seed = int(np.random.SeedSequence().generate_state(1)[0]) if seed is None else int(seed)

    # Each construct gets its own stream off the one seed, so adding a construct does not
    # reshuffle the others.
    results = [scan_system(s, temps_c, seed=[seed, i], t_ref_c=t_ref_c,
                           fc_threshold=fc_threshold, confidence=confidence,
                           n_samples=n_samples, dH_priors=dH_priors,
                           luciferase_activity=luciferase_activity)
               for i, s in enumerate(systems)]

    return {
        "temperatures_c": [float(t) for t in temps_c],
        "systems": results,
        "t_ref_c": t_ref_c,
        "fc_threshold": fc_threshold,
        "confidence": confidence,
        "n_samples": n_samples,
        "dH_priors": resolve_priors(dH_priors),
        "seed": seed,
        "notes": _notes(results, dH_priors, luciferase_activity, t_ref_c),
    }
