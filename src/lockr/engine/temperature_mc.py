"""Monte Carlo over uncertain thermodynamics, for the temperature question."""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from . import landscape
from .temperature import MODE_ISOTHERMAL, ThermoSpec

DEFAULT_DRAWS = 500
MAX_DRAWS = 5000
MAX_EVALUATIONS = 240_000   # ~2s; 500 draws x 120 temps is 60k and takes ~0.5s

DEFAULT_STEPS = 120

# Temperatures a team actually has to work at.
CHECKPOINTS = ((277.15, "fridge, 4 C"), (298.15, "bench, 25 C"), (310.15, "body, 37 C"))

# How far each input is assumed to wander. Assumptions, not measurements: dCp is the
# widest because estimate_dCp() derives it from buried surface area rather than a calorimeter.
DEFAULT_SIGMAS = {
    "dH_latch": 2.0,        # kcal/mol
    "dH_key": 1.5,          # kcal/mol
    "dCp_latch": 0.4,       # kcal/mol/K
    "k_open_log10": 0.15,   # decades, ~ +-40%
    "k_ck_log10": 0.15,     # decades
    "pull_rel": 0.10,       # fraction of the value
}

SIGMA_KEYS = tuple(DEFAULT_SIGMAS)

_SIGMA_LABELS = {
    "dH_latch": "cage heat response (dH), +-{v:g} kcal/mol",
    "dH_key": "key-cage heat response (dH), +-{v:g} kcal/mol",
    "dCp_latch": "cage dCp, +-{v:g} kcal/mol/K",
    "k_open_log10": "cage looseness, +-{v:g} decades",
    "k_ck_log10": "key grip strength, +-{v:g} decades",
    "pull_rel": "target pull strength, +-{pct:g}%",
}

ASSUMPTION_WARNING = (
    "These spreads are assumptions about how uncertain your inputs are, not measurements. "
    "A confidence range built on guessed spreads looks far more certain than it is. Replace "
    "them with your own error bars before quoting these numbers."
)


def resolve_sigmas(overrides: dict | None = None) -> dict[str, float]:
    sigmas = dict(DEFAULT_SIGMAS)
    for key, value in (overrides or {}).items():
        if key not in DEFAULT_SIGMAS:
            raise ValueError(f"unknown spread {key!r}; use {SIGMA_KEYS}")
        if value is None:
            continue
        if value < 0:
            raise ValueError(f"{key} spread must be >= 0")
        sigmas[key] = float(value)
    return sigmas


def describe_sigmas(sigmas: dict[str, float]) -> list[str]:
    return [_SIGMA_LABELS[k].format(v=v, pct=v * 100) for k, v in sigmas.items() if v]


def _perturbed(spec: ThermoSpec, dH: float | None, dCp: float | None) -> ThermoSpec:
    """A spec with its heat terms nudged; an isothermal spec is left alone."""
    if spec.dH_ref is None:
        return spec
    changes = {"dH_ref": spec.dH_ref + dH}
    if dCp is not None:
        changes["dCp"] = spec.dCp + dCp
    return replace(spec, **changes)


def _percentiles(values: np.ndarray, axis=None) -> dict:
    lo, mid, hi = np.nanpercentile(values, [5, 50, 95], axis=axis)
    return {"p05": lo, "median": mid, "p95": hi}


def _as_list(values) -> list:
    return [None if not math.isfinite(v) else float(v) for v in np.atleast_1d(values)]


def scan(base: dict, specs: dict[str, ThermoSpec], T_lo: float, T_hi: float,
         metric: str = "fold_change", sigmas: dict | None = None,
         draws: int = DEFAULT_DRAWS, steps: int = DEFAULT_STEPS,
         threshold: float | None = None, op: str = ">=", seed: int = 0) -> dict:
    """Spread of the metric over temperature, given uncertain inputs."""
    if metric not in landscape.METRICS:
        raise ValueError(f"metric must be one of {landscape.METRICS}")
    if op not in landscape._OPS:
        raise ValueError(f"op must be one of {tuple(landscape._OPS)}")
    if not 1 <= draws <= MAX_DRAWS:
        raise ValueError(f"draws must be between 1 and {MAX_DRAWS}")
    iso = [k for k in ("dG_latch_cage", "dG_key_cage") if specs[k].mode is MODE_ISOTHERMAL]
    if iso:
        raise ValueError(f"{', '.join(iso)} has no dH_ref, so temperature changes nothing "
                         f"and there is no spread to sample")

    warnings = []
    if draws * steps > MAX_EVALUATIONS:
        steps = max(8, MAX_EVALUATIONS // draws)
        warnings.append(f"Reduced to {steps} temperature points to stay under "
                        f"{MAX_EVALUATIONS} evaluations.")

    sigmas = resolve_sigmas(sigmas)
    rng = np.random.default_rng(seed)
    temperatures = landscape.axis_values(T_lo, T_hi, steps, "linear")

    shift = lambda key: rng.normal(0.0, sigmas[key], draws)
    dH_latch, dH_key = shift("dH_latch"), shift("dH_key")
    dCp_latch = shift("dCp_latch")
    k_open = base["k_open"] * 10.0 ** shift("k_open_log10")
    k_ck = base["k_ck"] * 10.0 ** shift("k_ck_log10")
    pull = base["pull"] * (1.0 + shift("pull_rel"))

    # One evaluate() call per draw: the temperature grid is vectorised inside it, but each
    # draw has its own dH/dCp so it needs its own specs.
    curves = np.empty((draws, len(temperatures)))
    for d in range(draws):
        draw_specs = dict(specs)
        draw_specs["dG_latch_cage"] = _perturbed(specs["dG_latch_cage"], dH_latch[d], dCp_latch[d])
        draw_specs["dG_key_cage"] = _perturbed(specs["dG_key_cage"], dH_key[d], None)
        values = {**base, "k_open": float(k_open[d]), "k_ck": float(k_ck[d]),
                  "pull": max(0.0, float(pull[d])), "temperature_K": temperatures}
        curves[d], _ = landscape.evaluate(values, metric, draw_specs)

    usable = np.isfinite(curves).any(axis=1)
    if not usable.any():
        raise ValueError("every simulated design is undefined over this temperature range")

    band = _percentiles(curves[usable], axis=0)
    best = _best_temperature(curves[usable], temperatures, metric)

    return {
        "metric": metric,
        "temperatures": [float(t) for t in temperatures],
        "median": _as_list(band["median"]),
        "p05": _as_list(band["p05"]),
        "p95": _as_list(band["p95"]),
        "best_temperature": best,
        "checkpoints": _checkpoints(curves[usable], temperatures, threshold, op),
        "draws": int(usable.sum()),
        "requested_draws": draws,
        "seed": seed,
        "threshold": threshold,
        "op": op,
        "assumptions": describe_sigmas(sigmas),
        "warnings": warnings + [ASSUMPTION_WARNING],
    }


def _best_temperature(curves: np.ndarray, temperatures: np.ndarray, metric: str) -> dict:
    """Where each simulated design peaks, and how often that is a real interior peak."""
    filled = np.where(np.isfinite(curves), curves,
                      np.inf if metric in landscape.LOWER_IS_BETTER else -np.inf)
    pick = np.nanargmin if metric in landscape.LOWER_IS_BETTER else np.nanargmax
    index = pick(filled, axis=1)
    best = temperatures[index]
    interior = ((index > 0) & (index < len(temperatures) - 1)).mean()
    stats = _percentiles(best)
    return {"median": float(stats["median"]), "p05": float(stats["p05"]),
            "p95": float(stats["p95"]), "interior_fraction": float(interior)}


def _checkpoints(curves: np.ndarray, temperatures: np.ndarray,
                 threshold: float | None, op: str) -> list[dict]:
    """The metric at the temperatures a team actually works at."""
    out = []
    for T, label in CHECKPOINTS:
        if not temperatures[0] - 1e-9 <= T <= temperatures[-1] + 1e-9:
            continue
        column = curves[:, int(np.argmin(np.abs(temperatures - T)))]
        finite = column[np.isfinite(column)]
        if not finite.size:
            continue
        stats = _percentiles(finite)
        entry = {"temperature_K": float(T), "label": label,
                 "median": float(stats["median"]), "p05": float(stats["p05"]),
                 "p95": float(stats["p95"]), "probability": None}
        if threshold is not None:
            entry["probability"] = float(landscape._OPS[op](finite, threshold).mean())
        out.append(entry)
    return out
