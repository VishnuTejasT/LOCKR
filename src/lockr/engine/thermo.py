"""General LOCKR fold-change model, free-energy helpers, regime diagnostic.

Three-state framework (cage closed/open, lucKey competing for the open state,
governed by K_open/K_CK) is from Langan et al. 2019 (Nature 572) and
Quijano-Rubio et al. 2021 (Nature 591) — applies to any lucCage-based sensor.
The closed-form fold-change expression below is my own instantiation of that
framework (ECLIPSE Thermodynamics doc, Section 7); not verbatim from either
paper. ECLIPSE numbers validating it live in tests/test_thermo_eclipse.py.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import curve_fit

from .models import DEFAULT_PARAMS, FoldChangeResult, RegimeResult, ScanResult, SensorParams


def theta(target_conc: float, Kd: float) -> float:
    # Fraction of binder occupied by target.
    return target_conc / (target_conc + Kd)


def k_open_eff(K_open: float, pull: float, theta: float) -> float:
    # Target binding pulls the latch open; pull is the allosteric coupling strength.
    return K_open * (1 + pull * theta)


def _f_open(k_open: float, params: SensorParams) -> float:
    # Fraction of cages in the open, signal-competent state.
    return k_open / (1 + k_open + params.luckey_ratio)


def fold_change(target_conc: float, Kd: float, pull: float,
                params: SensorParams = DEFAULT_PARAMS) -> float:
    th = theta(target_conc, Kd)
    koe = k_open_eff(params.K_open, pull, th)
    return _f_open(koe, params) / _f_open(params.K_open, params)


def fold_change_detail(target_conc: float, Kd: float, pull: float,
                       params: SensorParams = DEFAULT_PARAMS) -> FoldChangeResult:
    # Same calc, but keeps every intermediate for a result card / debugging.
    th = theta(target_conc, Kd)
    koe = k_open_eff(params.K_open, pull, th)
    f_base = _f_open(params.K_open, params)
    f_signal = _f_open(koe, params)
    return FoldChangeResult(target_conc, Kd, pull, th, koe, f_base, f_signal, f_signal / f_base)


def f_base(params: SensorParams = DEFAULT_PARAMS) -> float:
    # Baseline open fraction with no target present.
    return _f_open(params.K_open, params)


def _saturating_fc(pull: float, params: SensorParams) -> float:
    # FC at theta -> 1 (saturating target): the realised max FC for this finite
    # pull. Distinct from params.luckey_ratio (the dominance ratio) and from the
    # true pull->infinity asymptote (1+K_open+luckey_ratio)/K_open, which this
    # codebase doesn't compute anywhere.
    koe = k_open_eff(params.K_open, pull, 1.0)
    return _f_open(koe, params) / _f_open(params.K_open, params)


def max_fold_change(Kd: float, pull: float, params: SensorParams = DEFAULT_PARAMS) -> float:
    # Kd is unused on purpose — this max is cage-set, not Kd-set; kept in the
    # signature to mirror fold_change's args at call sites.
    return _saturating_fc(pull, params)


def scan_dose_response(Kd: float, pull: float, params: SensorParams = DEFAULT_PARAMS,
                       n: int = 500) -> ScanResult:
    # Titrate target across ~9 decades, read off max FC / EC50 / LOD (0.1*EC50).
    conc = np.logspace(-14, -5, n)
    fcs = np.array([fold_change(c, Kd, pull, params) for c in conc])
    mfc = float(fcs.max())
    half = (mfc + 1) / 2
    ec50 = float(conc[np.argmin(np.abs(fcs - half))])
    return ScanResult("", Kd, pull, mfc, ec50, ec50 * 0.1)


def dg_open_cost(params: SensorParams = DEFAULT_PARAMS) -> float:
    # Cost to crack the latch open: -RT*ln(K_open).
    return -params.RT * math.log(params.K_open)


def dg_luckey(params: SensorParams = DEFAULT_PARAMS) -> float:
    # Energy released by lucKey binding once the latch is open.
    return -params.RT * math.log(params.luckey_ratio)


def dg_from_kd(Kd: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    # Boltzmann relation: dG = RT*ln(Kd).
    return RT * math.log(Kd)


def kd_from_dg(dG: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    return math.exp(dG / RT)


def kd_from_ddg(kd_ref: float, ddg: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    # Rescales a reference Kd by a binding-energy change; ddg < 0 means tighter.
    return kd_ref * math.exp(ddg / RT)


# How hard the regime diagnostic pushes K_open to test if it moves fold-change.
_K_OPEN_PROBE_FACTOR = 30.0
# Sensitivity thresholds on that probe — general heuristic, not biology.
_KEY_LIMITED_BELOW = 0.02
_KOPEN_LIMITED_ABOVE = 0.08


def diagnose_regime(params: SensorParams = DEFAULT_PARAMS, pull: float = 10.0) -> RegimeResult:
    """Key-limited vs K_open-limited, for any K_open/K_CK/lucKey.

    Probes K_open directly rather than comparing raw magnitudes against the
    lucKey/K_CK ratio, since a magnitude compare alone doesn't discriminate.
    pull=10 is a generic default (literature LOCKR designs run ~10-20x), not
    tied to any one sensor.
    """
    ratio = params.luckey_ratio
    mfc = _saturating_fc(pull, params)

    probed = SensorParams(K_open=params.K_open * _K_OPEN_PROBE_FACTOR,
                          K_CK=params.K_CK, lucKey=params.lucKey, RT=params.RT)
    mfc_probed = _saturating_fc(pull, probed)
    rel_change = abs(mfc_probed - mfc) / mfc

    if rel_change < _KEY_LIMITED_BELOW:
        regime, helps = "key-limited", False
        verdict = (f"Key-limited: lucKey/K_CK = {ratio:.1f} dominates over "
                   f"K_open = {params.K_open:g}; fold-change tops out near "
                   f"{mfc:.1f}x at this pull and latch tuning won't move it. "
                   f"Raise lucKey or tighten K_CK instead.")
    elif rel_change > _KOPEN_LIMITED_ABOVE:
        regime, helps = "K_open-limited", True
        verdict = (f"K_open-limited: lucKey/K_CK = {ratio:.1f} is comparable to "
                   f"K_open, so latch tuning materially affects fold-change.")
    else:
        regime, helps = "mixed", True
        verdict = (f"Mixed: lucKey/K_CK = {ratio:.1f} — both the latch and the "
                   f"key constrain fold-change here.")

    return RegimeResult(ratio, params.K_open, regime, mfc, helps, verdict)


def fit_pull_strength(target_conc, fc_measured, Kd: float,
                      params: SensorParams = DEFAULT_PARAMS):
    # Pull can't be computed — back it out from a measured titration curve.
    target_conc = np.asarray(target_conc, dtype=float)
    fc_measured = np.asarray(fc_measured, dtype=float)

    def model(conc, pull):
        return np.array([fold_change(c, Kd, pull, params) for c in conc])

    popt, pcov = curve_fit(model, target_conc, fc_measured, p0=[10.0], bounds=(0, 100))
    return float(popt[0]), float(np.sqrt(pcov[0, 0]))
