"""LucCage's latch is a helix, so the binder MUST be a helix to work within the LOCKR system. 
This module estimated the helix propensity based on the sequence alone (utilizes the Pace & Scholtz 1998 scale). 
It also flags residue patterns that are common in binders that are cyclic. Only a confirmed cyclic peptide will block the graft from even occuring,
 but this confidence rating is just a warning because the LucCage scaffold can graft barealy-nonhelical binders too."""

from __future__ import annotations

from .models import CyclizationEvidence, HelixIssue, HelixReport

# Pace & Scholtz (1998) helix propensities, with kcal/mol relative to alanine. LOWER = MORE HELICAL.
_HELIX_DDG = {
    "A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40,
    "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61,
    "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69, "G": 1.00, "P": 3.16,
}
_WORST_DDG = 1.00  # Glycine. Proline is handled separately so it can't swamp the average.


_GOOD_N_CAP = set("DNST")


_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

_HELIX_PERIOD_DEG = 100.0  

_ACIDIC = set("DE")
_BASIC = set("KR")

_BAND_LIKELY_MIN = 0.60
_BAND_UNCERTAIN_MIN = 0.40

_MIN_LENGTH_FOR_HELIX = 7 


def _mean_propensity_score(sequence: str) -> float:
    """Average helix propensity mapped to 0-1, where 1 is a perfect poly-alanine helix."""
    if not sequence:
        return 0.0
    #Proline would drag the conf down by a landslide.
    total = sum(min(_HELIX_DDG.get(aa, _WORST_DDG), _WORST_DDG) for aa in sequence)
    mean = total / len(sequence)
    return max(0.0, 1.0 - mean / _WORST_DDG)


def _hydrophobic_moment(sequence: str) -> float:
    """An amphipathic helix is more likely to be stable and graftable because its hydrophobic residues are clustered on one face, promoting interactions with hydrophobic surfaces."""
    import math

    if not sequence:
        return 0.0
    x = y = 0.0
    for i, aa in enumerate(sequence):
        h = _HYDROPATHY.get(aa, 0.0)
        angle = math.radians(_HELIX_PERIOD_DEG * i)
        x += h * math.cos(angle)
        y += h * math.sin(angle)
    return math.hypot(x, y) / len(sequence)


def _salt_bridges(sequence: str) -> list[tuple[int, int]]:
    """i to i+3/i+4 acid-base pairs land on the same helical face and can form a stabilizing salt bridge.  
    This could destroy these possibilities without knowing."""
    pairs = []
    for i, aa in enumerate(sequence):
        for gap in (3, 4):
            j = i + gap
            if j >= len(sequence):
                continue
            other = sequence[j]
            if (aa in _ACIDIC and other in _BASIC) or (aa in _BASIC and other in _ACIDIC):
                pairs.append((i + 1, j + 1))
    return pairs


def _find_issues(sequence: str) -> list[HelixIssue]:
    issues = []
    n = len(sequence)

    if n < _MIN_LENGTH_FOR_HELIX:
        issues.append(HelixIssue(
            position=None, severity="blocking", kind="too_short",
            message=f"{n} residues is under two helical turns, so this will block the graft.",
        ))

    for i, aa in enumerate(sequence, 1):
        if aa == "P" and i > 4: 
            issues.append(HelixIssue(
                position=i, severity="blocking", kind="internal_proline",
                message=f"Proline at {i} kinks the helix, and it wont thread as a straight latch",
            ))
        elif aa == "P":
            issues.append(HelixIssue(
                position=i, severity="info", kind="n_terminal_proline",
                message=f"Proline at {i} is in the first turn, so grafting will be fine.",
            ))

    # Glycine is flexible instead of being kinked, so it costs stability.
    internal_gly = [i for i, aa in enumerate(sequence, 1) if aa == "G" and 1 < i < n]
    if len(internal_gly) >= 2:
        issues.append(HelixIssue(
            position=internal_gly[0], severity="warning", kind="glycine_rich",
            message=f"Glycines at {', '.join(str(p) for p in internal_gly)} make the backbone floppy.",
        ))

    if n >= _MIN_LENGTH_FOR_HELIX and sequence[0] not in _GOOD_N_CAP:
        issues.append(HelixIssue(
            position=1, severity="info", kind="weak_n_cap",
            message=f"{sequence[0]} is a weak N-cap, D/N/S/T hold the helix start better.",
        ))

    return issues


def _cyclization_evidence(sequence: str) -> CyclizationEvidence:
    cys = [i for i, aa in enumerate(sequence, 1) if aa == "C"]
    signals = []
    n = len(sequence)


    if len(cys) >= 2:
        if cys[0] <= 3 and cys[-1] >= n - 2:
            signals.append(
                f"The cysteines at {cys[0]} and {cys[-1]} show the pattern for a disulfide-cyclized peptide."
            )
        else:
            signals.append(f"{len(cys)} Cysteines (at {', '.join(str(p) for p in cys)}) could form a disulfide!")

    return CyclizationEvidence(
        possibly_cyclic=bool(signals),
        cysteine_positions=cys,
        signals=signals,
    )


def _band(confidence: float) -> str:
    if confidence >= _BAND_LIKELY_MIN:
        return "likely helical"
    if confidence >= _BAND_UNCERTAIN_MIN:
        return "uncertain"
    return "unlikely helical"



_W_PROPENSITY = 0.75
_W_MOMENT = 0.15
_W_SALT_BRIDGE = 0.10

_MOMENT_SATURATION = 2.0
_SALT_BRIDGE_SATURATION = 3


def analyze_helix(sequence: str) -> HelixReport:
    seq = sequence.strip().upper()

    propensity = _mean_propensity_score(seq)
    moment = _hydrophobic_moment(seq)
    bridges = _salt_bridges(seq)
    issues = _find_issues(seq)

    confidence = (
        _W_PROPENSITY * propensity
        + _W_MOMENT * min(moment / _MOMENT_SATURATION, 1.0)
        + _W_SALT_BRIDGE * min(len(bridges) / _SALT_BRIDGE_SATURATION, 1.0)
    )
   
    if any(i.severity == "blocking" for i in issues):
        confidence = min(confidence, _BAND_UNCERTAIN_MIN - 0.01)

    return HelixReport(
        sequence=seq,
        helix_confidence=confidence,
        band=_band(confidence),
        mean_propensity=propensity,
        hydrophobic_moment=moment,
        salt_bridges=bridges,
        issues=issues,
        cyclization=_cyclization_evidence(seq),
        graft_blocked=any(i.severity == "blocking" for i in issues),
    )



def compare_for_variant(original: str, variant: str) -> list[str]:
    """This will tell you whether a charge-optimized variant is likely to be less helical than the original, and therefore less likely to graft well."""
    orig = original.strip().upper()
    var = variant.strip().upper()
    if len(orig) != len(var):
        return []

    changed = [i for i, (a, b) in enumerate(zip(orig, var), 1) if a != b]
    if not changed:
        return []

    notes = []

    if 1 in changed and orig[0] in _GOOD_N_CAP and var[0] not in _GOOD_N_CAP:
        notes.append(
            f"Position 1 changed {orig[0]} to {var[0]}, losing the N-cap. The helix will fray "
            "from that end, and the graft will be less stable and more likely to fail."
        )

    lost = set(_salt_bridges(orig)) - set(_salt_bridges(var))
    if lost:
        pairs = ", ".join(f"{i}-{j}" for i, j in sorted(lost))
        notes.append(f"Salt bridges lost at {pairs}, and each was worth ~0.1-0.5 kcal/mol of helix stability.")

    dipole_lost = [p for p in changed if p <= 4 and orig[p - 1] in _ACIDIC and var[p - 1] not in _ACIDIC]
    if dipole_lost:
        positions = ", ".join(str(p) for p in dipole_lost)
        notes.append(
            f"Acidic residues at {positions} were neutralized, and it removed a favorable "
            "interaction with the helix dipole near the N-terminus."
        )

    return notes
