"""Net charge at pH and a coarse helical-propensity check for the latch graft.

Both are inputs to the liability scan: net charge is the headline number the
Scanner shows live, and the helix flag guards against a "fix" that removes a
charge but wrecks the latch helix (the grafts are alpha-helical).
"""

from __future__ import annotations

from .models import ChargeResult

# Side-chain pKa values (EMBOSS set) plus termini. Good enough for a net-charge
# readout on short peptides; we're not claiming titration-curve accuracy.
_PKA_SIDE = {"D": 3.65, "E": 4.25, "C": 8.5, "Y": 10.07,
             "H": 6.0, "K": 10.53, "R": 12.5}
_PKA_NTERM = 8.6
_PKA_CTERM = 3.6
_ACIDIC = set("DEYC")   # lose a proton -> negative
_BASIC = set("HKR")     # gain a proton -> positive


def _protonated_fraction(pKa: float, pH: float) -> float:
    return 1.0 / (1.0 + 10 ** (pH - pKa))


def net_charge(sequence: str, pH: float = 7.4) -> float:
    seq = sequence.strip().upper()
    q = 0.0
    # termini
    q += _protonated_fraction(_PKA_NTERM, pH)           # N-term: protonated = +1
    q -= (1 - _protonated_fraction(_PKA_CTERM, pH))     # C-term: deprotonated = -1
    for aa in seq:
        if aa in _BASIC:
            q += _protonated_fraction(_PKA_SIDE[aa], pH)
        elif aa in _ACIDIC:
            q -= (1 - _protonated_fraction(_PKA_SIDE[aa], pH))
    return q


# Chou-Fasman helix propensities (P_alpha). >1 favours helix, <1 disfavours.
# Only used as a rough flag, so the exact table matters less than P/G being clear
# breakers — which is what actually trips a latch graft.
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
    # Proline and glycine break helices. A terminal P/G is usually fine; an
    # internal one is the thing to flag for a latch graft.
    seq = sequence.strip().upper()
    return [i for i, aa in enumerate(seq, 1) if aa in "PG" and 1 < i < len(seq)]


def analyze_charge(sequence: str, pH: float = 7.4) -> ChargeResult:
    breakers = helix_breakers(sequence)
    prop = helix_propensity(sequence)
    # "ok" = leans helical and no internal P/G to kink the latch.
    ok = prop >= 1.0 and not breakers
    note = "" if ok else "low helix propensity or internal P/G — check latch register"
    return ChargeResult(net_charge(sequence, pH), pH, ok, breakers, note)
