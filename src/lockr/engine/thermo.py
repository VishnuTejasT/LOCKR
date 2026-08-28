"""General LOCKR fold-change model, free-energy helpers, regime diagnostic.

This module runs a 3-state framework that calculates the following: cage close/open energy, luckey competition
in open/close states and overall system fold-change. 

All these calculations and models are based on Langan et al. 2019 (Nature 572) and Quijano-Rubio et al. 2021 (Nature 591)
and applies to any LucCage based sensor system.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy.optimize import curve_fit, minimize_scalar

from .models import DEFAULT_PARAMS, FoldChangeResult, LodResult, RegimeResult, ScanResult, SensorParams


def theta(target_conc: float, Kd: float) -> float:
    return target_conc / (target_conc + Kd)


def k_open_eff(K_open: float, pull: float, theta: float) -> float:
    return K_open * (1 + pull * theta)


def _f_open(k_open: float, params: SensorParams) -> float:
    # 3-state partition function; only open+key-bound glows (to luminescence, the cage must already be open, and the key must be bound.)
    weight_signal = k_open * params.luckey_ratio
    return weight_signal / (1 + k_open + weight_signal)


def fold_change(target_conc: float, Kd: float, pull: float,
                params: SensorParams = DEFAULT_PARAMS) -> float:
    th = theta(target_conc, Kd)
    koe = k_open_eff(params.K_open, pull, th)
    return _f_open(koe, params) / _f_open(params.K_open, params)


def fold_change_detail(target_conc: float, Kd: float, pull: float,
                       params: SensorParams = DEFAULT_PARAMS) -> FoldChangeResult:
    th = theta(target_conc, Kd)
    koe = k_open_eff(params.K_open, pull, th)
    f_base = _f_open(params.K_open, params)
    f_signal = _f_open(koe, params)
    return FoldChangeResult(target_conc, Kd, pull, th, koe, f_base, f_signal, f_signal / f_base)


def f_base(params: SensorParams = DEFAULT_PARAMS) -> float:
    return _f_open(params.K_open, params)


def _saturating_fc(pull: float, params: SensorParams) -> float:
    # Realized max FC at this finite pull. This is not the dominance ratio and not the pull->infinity asymptote.
    koe = k_open_eff(params.K_open, pull, 1.0)
    return _f_open(koe, params) / _f_open(params.K_open, params)


def max_fold_change(Kd: float, pull: float, params: SensorParams = DEFAULT_PARAMS) -> float:
    return _saturating_fc(pull, params)


def scan_dose_response(Kd: float, pull: float, params: SensorParams = DEFAULT_PARAMS,
                       n: int = 500) -> ScanResult:
    conc = np.logspace(-14, -5, n)
    fcs = np.array([fold_change(c, Kd, pull, params) for c in conc])
    mfc = float(fcs.max())
    half = (mfc + 1) / 2
    ec50 = float(conc[np.argmin(np.abs(fcs - half))])
    return ScanResult("", Kd, pull, mfc, ec50, ec50 * 0.1)


def lod_and_ec50(Kd: float | None, pull: float, params: SensorParams = DEFAULT_PARAMS,
                 n: int = 500) -> LodResult:
    if Kd is None:  # This assumes the target is saturated, so the LOD is not meaningful.
        return LodResult(None, None, None)

    conc = np.logspace(-15, -5, n)
    fcs = np.array([fold_change(c, Kd, pull, params) for c in conc])
    mfc = float(fcs.max())

    def _first_crossing(threshold: float) -> float | None:
        hits = np.flatnonzero(fcs >= threshold)
        return float(conc[hits[0]]) if hits.size else None

    half = (mfc + 1) / 2
    ec50 = float(conc[np.argmin(np.abs(fcs - half))])
    return LodResult(_first_crossing(2.0), _first_crossing(3.0), ec50)


def k_open_from_destab(k_open_current: float, n_mutations: int, destab_per_mut: float,
                       RT: float = DEFAULT_PARAMS.RT) -> float:
    # Each mutation lowers destab_per_mut off the cost to open (Langan 2019-style estimate).
    dg_cost = -RT * math.log(k_open_current)
    return math.exp(-(dg_cost - n_mutations * destab_per_mut) / RT)


def dg_open_cost(params: SensorParams = DEFAULT_PARAMS) -> float:
    return -params.RT * math.log(params.K_open)


def dg_luckey(params: SensorParams = DEFAULT_PARAMS) -> float:
    return -params.RT * math.log(params.luckey_ratio)


def dg_from_kd(Kd: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    return RT * math.log(Kd)


def kd_from_dg(dG: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    return math.exp(dG / RT)


def kd_from_ddg(kd_ref: float, ddg: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    # Rescales a reference Kd by a binding-energy change; ddg < 0 means tighter.
    return kd_ref * math.exp(ddg / RT)



_K_OPEN_SEARCH_BOUNDS = (1e-6, 0.5)
_LUCKEY_SEARCH_BOUNDS = (1e-12, 1e-3)
_NEGLIGIBLE_GAIN = 0.02  
_DOMINANT_MARGIN = 2.0 


def _best_along(vary: Callable[[float], SensorParams], low: float, high: float,
                pull: float) -> tuple[float, float]:
    def neg_mfc(log_x: float) -> float:
        return -_saturating_fc(pull, vary(10 ** log_x))

    res = minimize_scalar(neg_mfc, bounds=(math.log10(low), math.log10(high)), method="bounded")
    return float(10 ** res.x), float(-res.fun)  


def diagnose_regime(params: SensorParams = DEFAULT_PARAMS, pull: float = 10.0) -> RegimeResult:
    """Which lever, K_open or lucKey/K_CK, has headroom left, and which direction to move it."""
    ratio = params.luckey_ratio
    mfc = _saturating_fc(pull, params)

    best_k_open, mfc_k_open_opt = _best_along(
        lambda k: SensorParams(K_open=k, K_CK=params.K_CK, lucKey=params.lucKey, RT=params.RT),
        *_K_OPEN_SEARCH_BOUNDS, pull)
    best_luckey, mfc_luckey_opt = _best_along(
        lambda lk: SensorParams(K_open=params.K_open, K_CK=params.K_CK, lucKey=lk, RT=params.RT),
        *_LUCKEY_SEARCH_BOUNDS, pull)

    gain_k_open = max(0.0, (mfc_k_open_opt - mfc) / mfc)
    gain_luckey = max(0.0, (mfc_luckey_opt - mfc) / mfc)

    def _sentence(s: str) -> str:
        return s[0].upper() + s[1:] if s else s  

    def _k_open_move() -> str:
        word = "raising" if best_k_open > params.K_open else "lowering"
        return (f"{word} K_open from {params.K_open:.3g} toward {best_k_open:.3g} "
                f"(Engineer the latch for a {'weaker' if word == 'raising' else 'stronger'} "
                f"closed state), and this would take fold-change to about {mfc_k_open_opt:.1f}x")

    def _luckey_move() -> str:
        word = "raising" if best_luckey > params.lucKey else "lowering"
        return (f"{word} lucKey from {params.lucKey * 1e9:.3g} nM toward "
                f"{best_luckey * 1e9:.3g} nM (or the equivalent move in K_CK) "
                f"would take fold-change to about {mfc_luckey_opt:.1f}x")

    if gain_k_open < _NEGLIGIBLE_GAIN and gain_luckey < _NEGLIGIBLE_GAIN:
        regime, helps = "mixed", False
        verdict = (f"Near-optimal: lucKey/K_CK = {ratio:.1f} and K_open = {params.K_open:g} "
                   f"are both already close to their maximum values at pull={pull:g} "
                   f"(fold-change {mfc:.1f}x, ceiling ~{max(mfc, mfc_k_open_opt, mfc_luckey_opt):.1f}x "
                   f"cage-latch allosteric coupling is the way to go, not tweaking lucKey/K_CK or K_open.")
                   
        recommendations = [
            "K_open and lucKey/K_CK are both already close to their local optimum for this pull.",
            "Raise pull (allosteric coupling strength) to go further. That means redesigning the "
            "cage-latch geometry, and not the concentrations or K_CK.",
        ]
    elif gain_luckey >= gain_k_open * _DOMINANT_MARGIN:
        regime, helps = "key-limited", gain_k_open >= _NEGLIGIBLE_GAIN
        verdict = (f"Key-limited: lucKey/K_CK = {ratio:.1f} has the most headroom, "
                   f"{_luckey_move()} (+{gain_luckey * 100:.0f}%). "
                   f"K_open tuning alone buys at most +{gain_k_open * 100:.0f}%.")
        recommendations = [_sentence(_luckey_move()) + "."]
        if not helps:
            recommendations.append("Latch (K_open) mutations will not meaningfully raise fold-change here.")
    elif gain_k_open >= gain_luckey * _DOMINANT_MARGIN:
        regime, helps = "K_open-limited", True
        verdict = (f"K_open-limited: {_k_open_move()} (+{gain_k_open * 100:.0f}%). "
                   f"lucKey/K_CK tuning alone buys at most +{gain_luckey * 100:.0f}%.")
        recommendations = [_sentence(_k_open_move()) + "."]
    else:
        regime, helps = "mixed", gain_k_open >= _NEGLIGIBLE_GAIN
        verdict = (f"Mixed: lucKey/K_CK = {ratio:.1f}, both axes have comparable headroom, "
                   f"{_k_open_move()} (+{gain_k_open * 100:.0f}%), and separately "
                   f"{_luckey_move()} (+{gain_luckey * 100:.0f}%).")
        recommendations = [_sentence(_k_open_move()) + ".", _sentence(_luckey_move()) + "."]

    return RegimeResult(ratio, params.K_open, regime, mfc, helps, verdict, recommendations)


def fit_pull_strength(target_conc, fc_measured, Kd: float,
                      params: SensorParams = DEFAULT_PARAMS):
    target_conc = np.asarray(target_conc, dtype=float)
    fc_measured = np.asarray(fc_measured, dtype=float)

    def model(conc, pull):
        return np.array([fold_change(c, Kd, pull, params) for c in conc])

    popt, pcov = curve_fit(model, target_conc, fc_measured, p0=[10.0], bounds=(0, 100))
    return float(popt[0]), float(np.sqrt(pcov[0, 0]))
