"""General charge-liability scan and K_CK penalty estimate.

Flags acidic residues outside caller-supplied preserve_positions and estimates
their effect on K_CK. The penalty model is my own instantiation of the lucCage
framework (Langan 2019 / Quijano-Rubio 2021), calibrated on the ECLIPSE anchor
sequences in calibration.py — override charge_penalty_per_residue for other
systems as more data comes in.
"""

from __future__ import annotations

import math

from . import calibration
from .charge import net_charge
from .models import BinderSequence, K_CK_DEFAULT, Liability, LiabilityReport, RT_37C, VariantSuggestion

# D/E -> A is the literal ECLIPSE fix; D->N/E->Q keeps shape and H-bonding.
_POLICIES = {
    "neutralizing": {"D": "A", "E": "A"},
    "conservative": {"D": "N", "E": "Q"},
}

# Liability-score band cutoffs (0-100 scale), from the UI gauge design — not biology.
_BAND_LOW_MAX = 33
_BAND_MODERATE_MAX = 66


def _dg_ck(K_CK_reference: float, RT: float) -> float:
    return -RT * math.log(1.0 / K_CK_reference)


def kck_from_penalty(penalty: float, K_CK_reference: float = K_CK_DEFAULT,
                     RT: float = RT_37C) -> float:
    # K_CK_estimate = exp(-(|dG_CK| - penalty)/RT): bigger penalty, weaker K_CK.
    return math.exp(-(abs(_dg_ck(K_CK_reference, RT)) - penalty) / RT)


def _band(score: float) -> str:
    if score <= _BAND_LOW_MAX:
        return "low"
    if score <= _BAND_MODERATE_MAX:
        return "moderate"
    return "high"


def _score_from_penalty(penalty: float, score_scale: float) -> float:
    # Saturating kcal/mol -> 0..100 map.
    return 100.0 * (1.0 - math.exp(-penalty / score_scale))


def scan_liability(sequence, preserve_positions,
                   charge_penalty_per_residue: float = calibration.PENALTY_PER_ACIDIC,
                   K_CK_reference: float = K_CK_DEFAULT, RT: float = RT_37C,
                   window=None, ph: float = 7.4, score_scale: float = 2.0) -> LiabilityReport:
    """Flag acidic liabilities outside preserve_positions and estimate K_CK.

    preserve_positions: the caller's own target-contact residues (1-indexed).
    ECLIPSE's PfLDH interface ([1,2,11,12,15]) is an example, not a default.
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
