"""This code showd the net charge at pH and a coarse helix-propensity as a checking step for ANY sequence...
This is not specific to ECLIPSE< but is instead generalized to any sequence based on acid-base chemistyr.
"""

from __future__ import annotations

from .models import ChargeResult

# EMBOSS side-chain pKa set, plus termini.
_PKA_SIDE = {"D": 3.65, "E": 4.25, "C": 8.5, "Y": 10.07,
             "H": 6.0, "K": 10.53, "R": 12.5}
_PKA_NTERM = 8.6
_PKA_CTERM = 3.6
_ACIDIC = set("DEYC")
_BASIC = set("HKR")


def _protonated_fraction(pKa: float, pH: float) -> float:
    return 1.0 / (1.0 + 10 ** (pH - pKa))


def net_charge(sequence: str, pH: float = 7.4) -> float:
    seq = sequence.strip().upper()
    # N-term protonated (+1), C-term deprotonated (-1) at full ionization.
    q = _protonated_fraction(_PKA_NTERM, pH) - (1 - _protonated_fraction(_PKA_CTERM, pH))
    for aa in seq:
        if aa in _BASIC:
            q += _protonated_fraction(_PKA_SIDE[aa], pH)
        elif aa in _ACIDIC:
            q -= (1 - _protonated_fraction(_PKA_SIDE[aa], pH))
    return q


# Chou-Fasman P_alpha helix propensities.
_PALPHA = {
    "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13, "G": 0.57, "H": 1.00,
    "I": 1.08, "K": 1.16, "L": 1.21, "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11,
    "R": 0.98, "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69,
}


def helix_propensity(sequence: str) -> float:
    seq = sequence.strip().upper()
    if not seq:
        return 0.0
    return sum(_PALPHA.get(aa, 1.0) for aa in seq) / len(seq)


def helix_breakers(sequence: str) -> list[int]:
    # Internal P/G kink a helix; terminal ones are usually fine.
    seq = sequence.strip().upper()
    return [i for i, aa in enumerate(seq, 1) if aa in "PG" and 1 < i < len(seq)]


def analyze_charge(sequence: str, pH: float = 7.4) -> ChargeResult:
    breakers = helix_breakers(sequence)
    ok = helix_propensity(sequence) >= 1.0 and not breakers
    note = "" if ok else "low helix propensity or internal P/G — check register"
    return ChargeResult(net_charge(sequence, pH), pH, ok, breakers, note)
