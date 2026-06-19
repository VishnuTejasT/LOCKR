"""Charge-liability scan and the K_CK penalty model (Charge doc §4.1 / Script 2).

The whole reason this module exists: the original RFDiffusion binder
LISDAELEAIFAEELDC has six D/E residues that sit at the lucKey-cage interface and
collapse K_CK, killing all signal even though the binder grips PfLDH beautifully.
We flag those residues, estimate the K_CK hit, and propose the charge-optimized
fix (LISAAALAAIFAAALAC).

K_CK penalty math is my own instantiation of the lucCage framework (Langan 2019 /
Quijano-Rubio 2021); the per-residue penalty is calibrated, see calibration.py.
"""

from __future__ import annotations

import math

from . import calibration
from .charge import net_charge
from .models import (BinderSequence, DEFAULT_PARAMS, Liability, LiabilityReport,
                     SensorParams, VariantSuggestion)

# Substitution policies. Neutralizing (D/E -> A) is the literal ECLIPSE fix that
# produced LISAAALAAIFAAALAC; conservative (D->N, E->Q) keeps shape and H-bonding.
_POLICIES = {
    "neutralizing": {"D": "A", "E": "A"},
    "conservative": {"D": "N", "E": "Q"},
}


def _dg_ck(params: SensorParams) -> float:
    # dG_CK = -RT*ln(1/K_CK); ~-10.9 kcal/mol at K_CK=1e-8.
    return -params.RT * math.log(1.0 / params.K_CK)


def grafted_kck(penalty: float, params: SensorParams = DEFAULT_PARAMS) -> float:
    """K_CK after the interface penalty (Script 2):

        K_CK_grafted = exp(-(abs(dG_CK) - penalty) / RT)

    which is K_CK_original * exp(penalty/RT) — the penalty makes lucKey binding
    *weaker*, so K_CK gets larger. With the v1 binder's 4.8 kcal/mol that's
    1e-8 -> ~3.3e-5 M, i.e. ~3300x weaker, and lucKey can no longer bind.

    NOTE: the PDFs print this grafted value as "~3e-12 M" in a couple of places.
    That's a documentation typo (Vishnu confirmed) — 3e-12 would be ~3300x
    *tighter*, the opposite of what interface repulsion does, and it contradicts
    the doc's own "3000x weaker" prose. The formula here is correct; the PDFs are
    being corrected to ~3.3e-5.
    """
    return math.exp(-(abs(_dg_ck(params)) - penalty) / params.RT)


def _band(score: float) -> str:
    if score <= calibration.BAND_LOW_MAX:
        return "low"
    if score <= calibration.BAND_MODERATE_MAX:
        return "moderate"
    return "high"


def _score_from_penalty(penalty: float) -> float:
    # Saturating map kcal/mol -> 0..100. Ties the score to the actual energetic
    # penalty rather than a bare residue count, so two D/E score worse than one.
    return 100.0 * (1.0 - math.exp(-penalty / calibration.SCORE_SCALE))


def scan(binder, params: SensorParams = DEFAULT_PARAMS, window=None,
         interface=None, ph: float = 7.4) -> LiabilityReport:
    """Flag acidic liabilities and estimate the K_CK hit for one binder.

    window: 1-indexed (start, end) sensitive region; defaults to the whole
    sequence. interface: PfLDH-contact positions to preserve (never counted as
    liabilities) — defaults to the ECLIPSE design's [1,2,11,12,15].
    """
    if not isinstance(binder, BinderSequence):
        binder = BinderSequence(binder)
    interface = list(calibration.PFLDH_INTERFACE if interface is None else interface)
    start, end = (1, len(binder)) if window is None else window

    liabilities = []
    for pos, aa in binder.residues():
        if aa not in "DE":
            continue
        if pos in interface:
            continue          # PfLDH contact — keep it, it's earning its charge
        if not (start <= pos <= end):
            continue          # outside the sensitive region the user cares about
        # Position weight is 1.0 for now (uniform inside the window). Hook for
        # helical-face weighting later — see plan §4.3 step 3.
        w = 1.0
        liabilities.append(Liability(pos, aa, w, calibration.PENALTY_PER_ACIDIC * w))

    penalty = sum(l.penalty for l in liabilities)
    score = _score_from_penalty(penalty)
    return LiabilityReport(
        binder=binder,
        liabilities=liabilities,
        preserved_interface=interface,
        net_charge=net_charge(binder.sequence, ph),
        penalty_total=penalty,
        liability_score=score,
        liability_band=_band(score),
        K_CK_grafted=grafted_kck(penalty, params),
    )


def suggest_variant(binder, policy: str = "neutralizing", interface=None,
                    window=None, params: SensorParams = DEFAULT_PARAMS) -> VariantSuggestion:
    """Propose a charge-optimized variant under the chosen substitution policy.

    Only touches flagged liabilities — PfLDH-interface D/E and anything outside
    the window are left alone. On the v1 binder, neutralizing reproduces the real
    fix: LISDAELEAIFAEELDC -> LISAAALAAIFAAALAC.
    """
    if policy not in _POLICIES:
        raise ValueError(f"unknown policy {policy!r}")
    if not isinstance(binder, BinderSequence):
        binder = BinderSequence(binder)
    sub = _POLICIES[policy]

    report = scan(binder, params, window, interface)
    flagged = {l.position: l.residue for l in report.liabilities}

    chars = list(binder.sequence)
    mutations = []
    for pos, old in flagged.items():
        new = sub[old]
        chars[pos - 1] = new
        mutations.append(f"{old}{pos}{new}")

    new_seq = "".join(chars)
    new_report = scan(new_seq, params, window, interface)
    return VariantSuggestion(
        policy=policy,
        sequence=new_seq,
        mutations=mutations,
        liability_score=new_report.liability_score,
        liability_band=new_report.liability_band,
        K_CK_grafted=new_report.K_CK_grafted,
    )
