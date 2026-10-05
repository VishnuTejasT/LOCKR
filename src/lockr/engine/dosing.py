"""Key concentration back-calculation for dosing a real Key/Cage titration.

Given a target occupancy and K_CK, work out how much synthesized Key peptide
to pipette into an assay of known Cage concentration. Unlike thermo.py's
fold_change() (which treats lucKey as being in such excess over Cage that its
free concentration never depletes), a bench titration can easily have Key and
Cage at comparable concentrations, so occupancy here comes from the exact
quadratic mass-balance solution instead of the pseudo-first-order theta()
shortcut used elsewhere in the engine.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from .models import DosingResult
from .thermo import theta

# Most pipettes bottom out around here; below it, dilute the stock instead of
# pipetting the raw volume.
MIN_PIPETTABLE_VOLUME = 0.5e-6  # L

_SEARCH_EXPAND_FACTOR = 10.0
_SEARCH_MAX_EXPANSIONS = 60


def bound_conc(key_total: float, cage_total: float, k_ck: float) -> float:
    s = key_total + cage_total + k_ck
    disc = max(s * s - 4 * key_total * cage_total, 0.0)  # guards float round-off below 0
    return (s - math.sqrt(disc)) / 2


def fraction_bound(key_total: float, cage_total: float, k_ck: float) -> float:
    if cage_total <= 0:
        raise ValueError("total_cage_conc must be > 0")
    return bound_conc(key_total, cage_total, k_ck) / cage_total


def fraction_bound_from_fold_change(fold_change: float) -> float:
    # fold_change here is the bound:unbound occupancy ratio, so theta(x, 1.0)
    # (already the engine's ratio->fraction relationship, x/(x+K)) with K=1
    # gives fraction_bound = fold_change / (1 + fold_change) without
    # re-deriving that algebra.
    return theta(fold_change, 1.0)


def solve_key_conc(target_fraction_bound: float, cage_total: float, k_ck: float) -> float:
    """Root-find [Key]_t so fraction_bound(...) hits the target, expanding the
    upper bracket until it brackets a root rather than guessing a fixed range."""

    def f(key_total: float) -> float:
        return fraction_bound(key_total, cage_total, k_ck) - target_fraction_bound

    lo, hi = 0.0, max(cage_total + k_ck, 1e-15)
    expansions = 0
    while f(hi) < 0 and expansions < _SEARCH_MAX_EXPANSIONS:
        hi *= _SEARCH_EXPAND_FACTOR
        expansions += 1
    if f(hi) < 0:
        raise ValueError("Could not bracket a Key concentration for this target; check K_CK and Cage conc.")
    return float(brentq(f, lo, hi, xtol=1e-15, rtol=1e-12))


def _titration_curve(cage_total: float, k_ck: float, center_key_conc: float,
                     n_points: int) -> list[tuple[float, float]]:
    scale = center_key_conc if center_key_conc > 0 else k_ck
    if not scale > 0:
        return []
    keys = np.logspace(math.log10(scale * 1e-3), math.log10(scale * 1e3), n_points)
    return [(float(k), float(fraction_bound(k, cage_total, k_ck))) for k in keys]


def dose_for_target(k_ck: float, total_cage_conc: float, key_stock_conc: float, assay_volume: float,
                    target_fraction_bound: float | None = None,
                    target_fold_change: float | None = None,
                    n_titration_points: int = 50) -> DosingResult:
    if k_ck <= 0 or total_cage_conc <= 0 or key_stock_conc <= 0 or assay_volume <= 0:
        raise ValueError("k_ck, total_cage_conc, key_stock_conc, and assay_volume must all be > 0")
    if (target_fraction_bound is None) == (target_fold_change is None):
        raise ValueError("Provide exactly one of target_fraction_bound or target_fold_change")

    if target_fraction_bound is None:
        target_fraction_bound = fraction_bound_from_fold_change(target_fold_change)

    if not 0 <= target_fraction_bound < 1:
        raise ValueError("target_fraction_bound must be in [0, 1); 1.0 would need infinite Key")

    key_conc = 0.0 if target_fraction_bound == 0 else solve_key_conc(target_fraction_bound, total_cage_conc, k_ck)
    achieved_fraction = fraction_bound(key_conc, total_cage_conc, k_ck)
    volume_to_add = key_conc * assay_volume / key_stock_conc

    warnings = []
    if volume_to_add > assay_volume:
        warnings.append(
            f"Required volume ({volume_to_add * 1e6:.3g} µL) exceeds the assay volume itself "
            f"({assay_volume * 1e6:.3g} µL); the Key stock isn't concentrated enough for this target, "
            "make a more concentrated stock."
        )
    elif 0 < volume_to_add < MIN_PIPETTABLE_VOLUME:
        warnings.append(
            f"Required volume ({volume_to_add * 1e6:.3g} µL) is below a reasonably pipettable minimum "
            f"(~{MIN_PIPETTABLE_VOLUME * 1e6:.3g} µL); dilute the stock first rather than pipetting this "
            "little directly."
        )

    titration = _titration_curve(total_cage_conc, k_ck, key_conc, n_titration_points)

    return DosingResult(key_conc=key_conc, fraction_bound=achieved_fraction,
                        volume_to_add=volume_to_add, warnings=warnings, titration=titration)
