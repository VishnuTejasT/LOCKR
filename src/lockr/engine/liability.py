"""Flags whether a grafted binder's charged residues weaken K_CK, and proposes a fix to the possible problem by providing a conservative or neutralizing approach. The suggested sequence gets checked against protected regions by assembly.py."""

from __future__ import annotations

import math

from . import calibration
from .charge import net_charge
from .models import BinderSequence, K_CK_DEFAULT, Liability, LiabilityReport, RT_37C, VariantSuggestion

#D/E -> A was the fix for us, but D->N/E->Q is way more conservative and keeps shape and hydrogen bonding too.
_POLICIES = {
    "neutralizing": {"D": "A", "E": "A"},
    "conservative": {"D": "N", "E": "Q"},
}

_BAND_LOW_MAX = 33
_BAND_MODERATE_MAX = 66


def _dg_ck(K_CK_reference: float, RT: float) -> float:
    return -RT * math.log(1.0 / K_CK_reference)


def kck_from_penalty(penalty: float, K_CK_reference: float = K_CK_DEFAULT,
                     RT: float = RT_37C) -> float:
    # K_CK_estimate = exp(-(|dG_CK| - penalty)/RT). The higher the penalty, the weaker K_CK is.
    return math.exp(-(abs(_dg_ck(K_CK_reference, RT)) - penalty) / RT)


def _band(score: float) -> str:
    if score <= _BAND_LOW_MAX:
        return "Low"
    if score <= _BAND_MODERATE_MAX:
        return "Moderate"
    return "High"


def _score_from_penalty(penalty: float, score_scale: float) -> float:
    """Maps a penalty to a 0-100 score, where 100 is no penalty and 0 is a huge penalty."""
    return 100.0 * (1.0 - math.exp(-penalty / score_scale))


def scan_liability(sequence, preserve_positions,
                   charge_penalty_per_residue: float = calibration.PENALTY_PER_ACIDIC,
                   K_CK_reference: float = K_CK_DEFAULT, RT: float = RT_37C,
                   window=None, ph: float = 7.4, score_scale: float = 2.0) -> LiabilityReport:
    """
    This will flag all acidic liability residues outside of the preserved position while also 
    estimating the K_CK.
    """
    if not isinstance(sequence, BinderSequence):
        sequence = BinderSequence(sequence)
    start, end = (1, len(sequence)) if window is None else window

    liabilities = []
    for pos, aa in sequence.residues():
        if aa not in "DE" or pos in preserve_positions or not (start <= pos <= end):
            continue
        liabilities.append(Liability(pos, aa, 1.0, charge_penalty_per_residue))

    penalty = sum(l.penalty for l in liabilities)
    score = _score_from_penalty(penalty, score_scale)
    return LiabilityReport(
        binder=sequence,
        liabilities=liabilities,
        preserve_positions=list(preserve_positions),
        net_charge=net_charge(sequence.sequence, ph),
        penalty_total=penalty,
        liability_score=score,
        liability_band=_band(score),
        K_CK_estimate=kck_from_penalty(penalty, K_CK_reference, RT),
    )


def suggest_variant(sequence, preserve_positions, policy: str = "neutralizing",
                    charge_penalty_per_residue: float = calibration.PENALTY_PER_ACIDIC,
                    K_CK_reference: float = K_CK_DEFAULT, RT: float = RT_37C,
                    window=None, score_scale: float = 2.0) -> VariantSuggestion:
    if policy not in _POLICIES:
        raise ValueError(f"unknown policy {policy!r}")
    if not isinstance(sequence, BinderSequence):
        sequence = BinderSequence(sequence)
    sub = _POLICIES[policy]

    report = scan_liability(sequence, preserve_positions, charge_penalty_per_residue,
                            K_CK_reference, RT, window, score_scale=score_scale)
    flagged = {l.position: l.residue for l in report.liabilities}

    chars = list(sequence.sequence)
    mutations = []
    for pos, old in flagged.items():
        new = sub[old]
        chars[pos - 1] = new
        mutations.append(f"{old}{pos}{new}")

    new_seq = "".join(chars)
    new_report = scan_liability(new_seq, preserve_positions, charge_penalty_per_residue,
                                K_CK_reference, RT, window, score_scale=score_scale)
    return VariantSuggestion(
        policy=policy,
        sequence=new_seq,
        mutations=mutations,
        liability_score=new_report.liability_score,
        liability_band=new_report.liability_band,
        K_CK_estimate=new_report.K_CK_estimate,
    )
