"""LOCKR fold-change model, free-energy helpers, and the regime diagnostic.

The three-state framework (cage closed/open, lucKey competing for the open
state, governed by K_open / K_CK) is from Langan et al. 2019 (Nature 572) and
Quijano-Rubio et al. 2021 (Nature 591). The closed-form fold-change expression
below is *my* instantiation of that framework (ECLIPSE Thermodynamics doc,
Section 7 / Scripts 7-8) — not a verbatim equation from either paper. Keep it as
written.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import curve_fit

from .models import DEFAULT_PARAMS, FoldChangeResult, RegimeResult, ScanResult, SensorParams


def theta(pfldh: float, Kd: float) -> float:
    # fraction of binder occupied by PfLDH
    return pfldh / (pfldh + Kd)


def k_open_eff(pull: float, th: float, params: SensorParams = DEFAULT_PARAMS) -> float:
    # PfLDH binding pulls the latch open; `pull` is the allosteric coupling
    # strength (only knowable from the luminescence fit — see fit_pull_strength).
    return params.K_open * (1 + pull * th)


def _f_open(k_open: float, params: SensorParams) -> float:
    # Fraction of cages in the open (signal-competent) state. The lucKey/K_CK
    # term is the competition that caps everything at 500 nM lucKey.
    return k_open / (1 + k_open + params.luckey_ratio)


def fold_change(pfldh: float, Kd: float, pull: float,
                params: SensorParams = DEFAULT_PARAMS) -> float:
    th = theta(pfldh, Kd)
    koe = k_open_eff(pull, th, params)
    f_base = _f_open(params.K_open, params)
    f_signal = _f_open(koe, params)
    return f_signal / f_base


def fold_change_detail(pfldh: float, Kd: float, pull: float,
                       params: SensorParams = DEFAULT_PARAMS) -> FoldChangeResult:
    """Same calculation but returns every intermediate, for the result card."""
    th = theta(pfldh, Kd)
    koe = k_open_eff(pull, th, params)
    f_base = _f_open(params.K_open, params)
    f_signal = _f_open(koe, params)
    return FoldChangeResult(pfldh, Kd, pull, th, koe, f_base, f_signal, f_signal / f_base)


def f_base(params: SensorParams = DEFAULT_PARAMS) -> float:
    # Baseline open fraction (no target). Exposed because the docs quote it
    # directly: 1e-3 / 51.001 at the ECLIPSE operating point.
    return _f_open(params.K_open, params)


def max_fold_change(Kd: float, pull: float,
                    params: SensorParams = DEFAULT_PARAMS) -> float:
    """FC ceiling for this (Kd, pull): the saturating-target limit (theta -> 1).

    Independent of Kd — affinity sets where you sit on the curve (EC50), not how
    high it goes. That's why v1.0 and v2.2 share the 11/21/31x ceiling.
    """
    koe = k_open_eff(pull, 1.0, params)
    return _f_open(koe, params) / _f_open(params.K_open, params)


def scan_dose_response(Kd: float, pull: float, params: SensorParams = DEFAULT_PARAMS,
                       n: int = 500) -> ScanResult:
    """Titrate PfLDH across ~9 decades and pull out max FC, EC50, LOD.

    EC50 is read off as the concentration at the half-maximal fold-change
    (matching Script 7's max_fc_ec50). LOD is the doc's 0.1*EC50 convention.
    """
    pfldh = np.logspace(-14, -5, n)
    fcs = np.array([fold_change(c, Kd, pull, params) for c in pfldh])
    mfc = float(fcs.max())
    half = (mfc + 1) / 2
    ec50 = float(pfldh[np.argmin(np.abs(fcs - half))])
    return ScanResult("", Kd, pull, mfc, ec50, ec50 * 0.1)


# --- free energies (Section 6 / Script 12) ---

def dg_open_cost(params: SensorParams = DEFAULT_PARAMS) -> float:
    # Cost to crack the latch open; ~+4.09 kcal/mol at K_open=1e-3.
    return -params.RT * math.log(params.K_open)


def dg_luckey(params: SensorParams = DEFAULT_PARAMS) -> float:
    # Energy lucKey binding contributes once the latch is open; ~-2.32 kcal/mol.
    return -params.RT * math.log(params.luckey_ratio)


def dg_from_kd(Kd: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    # Boltzmann: dG = RT*ln(Kd), so Kd = exp(dG/RT).
    return RT * math.log(Kd)


def kd_from_dg(dG: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    return math.exp(dG / RT)


def kd_from_ddg(kd_ref: float, ddg: float, RT: float = DEFAULT_PARAMS.RT) -> float:
    # How a binding-energy change rescales Kd (used for v2.2 vs v1.0: -4.6 kcal/mol
    # -> 2369x tighter). ddg < 0 means tighter binding.
    return kd_ref * math.exp(ddg / RT)


# --- regime diagnostic (Section 8) — the headline feature ---

# How much a plausible latch-destabilising mutation can raise K_open. The v1.1
# "moderate" scenario (I346S+I356S, ~2 kcal/mol) lands K_open near here, so this
# is the right yardstick for "would tuning the latch actually do anything".
_K_OPEN_PROBE_FACTOR = 30.0

# Sensitivity thresholds on max-FC when K_open is probed upward. Calibrated to
# reproduce the doc's Table 5 bands (ratio 50 -> no benefit; ratio ~10 -> marginal;
# ratio ~1 -> K_open matters). Refine if the framing changes.
_KEY_LIMITED_BELOW = 0.02   # <2% FC change -> latch won't help
_KOPEN_LIMITED_ABOVE = 0.08  # >8% FC change -> latch tuning moves the needle


def diagnose_regime(params: SensorParams = DEFAULT_PARAMS, pull: float = 10.0,
                    Kd: float = 100e-12) -> RegimeResult:
    """Is this sensor key-limited or K_open-limited?

    The prompt framing is "compare [lucKey]/K_CK against K_open", but a raw
    magnitude compare doesn't discriminate (50 and 1 are both >> 1e-3). What
    actually matters — and what the wet lab can act on — is whether raising
    K_open moves the fold-change at all. So we probe it directly: bump K_open
    ~30x and see if max FC budges. At 500 nM lucKey it doesn't (ratio 50,
    key-limited); drop lucKey to ~10 nM (ratio 1) and it does.

    Returns the dominance ratio and the realised max FC separately — the ratio is
    a diagnostic, not an achievable fold-change.
    """
    ratio = params.luckey_ratio
    mfc = max_fold_change(Kd, pull, params)

    probed = SensorParams(K_open=params.K_open * _K_OPEN_PROBE_FACTOR,
                          K_CK=params.K_CK, lucKey=params.lucKey, RT=params.RT)
    mfc_probed = max_fold_change(Kd, pull, probed)
    rel_change = abs(mfc_probed - mfc) / mfc

    # Note the two numbers in the verdicts are kept separate on purpose: the
    # dominance ratio (50) is the diagnostic; the realised max FC (~1+pull) is
    # what you actually get. Calling 50 a "ceiling" would be off by ~1/K_open.
    if rel_change < _KEY_LIMITED_BELOW:
        regime, helps = "key-limited", False
        verdict = (f"Key-limited regime. lucKey competition dominates "
                   f"([lucKey]/K_CK = {ratio:.0f} vs K_open = {params.K_open:g}), so "
                   f"the realised fold-change tops out near {mfc:.0f}x at this pull "
                   f"and tuning the latch (K_open) won't raise it. To improve, raise "
                   f"[lucKey] or tighten K_CK.")
    elif rel_change > _KOPEN_LIMITED_ABOVE:
        regime, helps = "K_open-limited", True
        verdict = (f"K_open-limited regime. The lucKey dominance ratio is only "
                   f"{ratio:.1f}, so the latch equilibrium materially affects "
                   f"fold-change — latch tuning is worth pursuing.")
    else:
        regime, helps = "mixed", True
        verdict = (f"Mixed regime. At [lucKey]/K_CK = {ratio:.1f}, both the latch "
                   f"and the key constrain fold-change; the plots show which gives "
                   f"more.")

    return RegimeResult(ratio, params.K_open, regime, mfc, helps, verdict)


# --- pull-strength fit (Section 10) ---

def fit_pull_strength(pfldh_conc, fc_measured, Kd: float = 100e-12,
                      params: SensorParams = DEFAULT_PARAMS):
    """Back out `pull` from a measured luminescence titration.

    Pull is the one parameter we can't get computationally — it's the allosteric
    coupling efficiency, fit from the experimental curve. Returns (pull, stderr).
    """
    pfldh_conc = np.asarray(pfldh_conc, dtype=float)
    fc_measured = np.asarray(fc_measured, dtype=float)

    def model(pf, pull):
        return np.array([fold_change(c, Kd, pull, params) for c in pf])

    popt, pcov = curve_fit(model, pfldh_conc, fc_measured, p0=[10.0], bounds=(0, 100))
    return float(popt[0]), float(np.sqrt(pcov[0, 0]))
