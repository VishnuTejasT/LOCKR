"""
A structure pre-check that runs before the charge liability scan. The latch of lucCage is a
helix, so a binder only threads into it if the binder is helical too. This module estimates
how helical a peptide is likely to be from its sequence alone, and separately looks for signs
that the peptide was meant to be cyclic (which cannot be grafted linearly at all).

Two honest limits worth stating up front:

1. This is a propensity estimate, not a structure prediction. For a real answer you need
   AlphaFold/ESMFold or knowledge of where the binder came from. What is here is the same
   physics those helix-coil models are built on, just without the statistical mechanics:
   per-residue helix propensities measured in host-guest peptide systems (Pace & Scholtz 1998),
   plus the capping, salt bridge and amphipathicity terms that dominate short peptide helices.

2. Cyclization is topology, not sequence. A head-to-tail amide or a disulfide staple leaves no
   trace in a one-letter string. All this can do is spot the residue patterns that make
   cyclization possible and ask the user, which is why it reports evidence rather than a verdict.

A peptide that is floppy on its own can still be helical once grafted, because the lucCage
scaffold templates the helix. So low helix confidence is a warning. The two things that
genuinely block a graft are an internal proline (it puts a permanent kink where the latch
needs a straight helix) and a confirmed cyclic peptide (there is no linear chain to thread).
"""

from __future__ import annotations

from .models import CyclizationEvidence, HelixIssue, HelixReport

# Pace & Scholtz (1998) helix propensities, kcal/mol relative to alanine. Lower is more
# helical. These come from substitutions in real peptide helices, unlike the Chou-Fasman
# frequencies in charge.py, which were counted in globular proteins and read the opposite way.
_HELIX_DDG = {
    "A": 0.00, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40,
    "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61,
    "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69, "G": 1.00, "P": 3.16,
}
_WORST_DDG = 1.00  # Glycine. Proline is handled separately so it can't swamp the average.

# Residues that cap a helix well. An N-cap sits one before the first helical residue and
# hydrogen bonds back to the free NH groups; D/N/S/T do it with their side chain oxygen.
_GOOD_N_CAP = set("DNST")

# Kyte-Doolittle hydropathy, used for the helical wheel moment.
_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

_HELIX_PERIOD_DEG = 100.0  # An alpha helix turns 100 degrees per residue.

_ACIDIC = set("DE")
_BASIC = set("KR")

# Bands for the 0-1 confidence. Tuned so a designed amphipathic binder lands in "likely",
# a mixed sequence in "uncertain", and something Pro/Gly rich in "unlikely".
_BAND_LIKELY_MIN = 0.60
_BAND_UNCERTAIN_MIN = 0.40

_MIN_LENGTH_FOR_HELIX = 7  # Under two turns there is no helix to speak of.


def _mean_propensity_score(sequence: str) -> float:
    """Average helix propensity mapped to 0-1, where 1 is a perfect poly-alanine helix."""
    if not sequence:
        return 0.0
    # Proline is capped at the glycine value here. Its real cost is 3.16, which would drag a
    # single-proline peptide to near zero on its own; the kink is reported as an issue instead.
    total = sum(min(_HELIX_DDG.get(aa, _WORST_DDG), _WORST_DDG) for aa in sequence)
    mean = total / len(sequence)
    return max(0.0, 1.0 - mean / _WORST_DDG)


def _hydrophobic_moment(sequence: str) -> float:
    """
    Eisenberg's moment: sum the hydropathy of each residue as a vector pointing where that
    residue sticks out around a helical wheel. A binder that grips a groove is usually
    amphipathic, all the greasy residues on one face, which shows up as a large moment.
    """
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
    """
    i to i+3 and i to i+4 acid/base pairs land on the same face of the helix and can form a
    stabilizing salt bridge. Positions are 1-based. These are what the charge optimizer is
    most likely to destroy without noticing.
    """
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
            message=f"{n} residues is under two helical turns, too short to graft as a helix.",
        ))

    for i, aa in enumerate(sequence, 1):
        # A proline in the first turn is fine and even common, it just starts the helix.
        # Anywhere else it breaks the backbone hydrogen bond pattern and kinks the chain.
        if aa == "P" and i > 4:
            issues.append(HelixIssue(
                position=i, severity="blocking", kind="internal_proline",
                message=f"Proline at {i} kinks the helix, it will not thread the straight latch.",
            ))
        elif aa == "P":
            issues.append(HelixIssue(
                position=i, severity="info", kind="n_terminal_proline",
                message=f"Proline at {i} is in the first turn, which is tolerated.",
            ))

    # Glycine is flexible rather than kinked, so it costs stability without blocking.
    internal_gly = [i for i, aa in enumerate(sequence, 1) if aa == "G" and 1 < i < n]
    if len(internal_gly) >= 2:
        issues.append(HelixIssue(
            position=internal_gly[0], severity="warning", kind="glycine_rich",
            message=f"Glycines at {', '.join(str(p) for p in internal_gly)} make the backbone floppy.",
        ))

    # Only the N-cap is reported. A weak C-cap is true of most sequences and there is rarely
    # anything useful to do about it, so it would be noise.
    if n >= _MIN_LENGTH_FOR_HELIX and sequence[0] not in _GOOD_N_CAP:
        issues.append(HelixIssue(
            position=1, severity="info", kind="weak_n_cap",
            message=f"{sequence[0]} is a weak N-cap, D/N/S/T hold the helix start better.",
        ))

    return issues


def _cyclization_evidence(sequence: str) -> CyclizationEvidence:
    """
    Sequence cannot prove a peptide is cyclic, so this collects the patterns that make
    cyclization possible and lets the caller ask the user.
    """
    cys = [i for i, aa in enumerate(sequence, 1) if aa == "C"]
    signals = []
    n = len(sequence)

    # Cysteine pairs are the only signal used here. Acid/base pairs one helix turn apart are
    # where a lactam staple would go, but they are also just how a designed helix stabilizes
    # itself, so flagging them would fire on nearly every good helical binder.
    if len(cys) >= 2:
        if cys[0] <= 3 and cys[-1] >= n - 2:
            signals.append(
                f"cysteines at {cys[0]} and {cys[-1]}, one near each end, the usual pattern for "
                "a disulfide-cyclized peptide"
            )
        else:
            signals.append(f"{len(cys)} cysteines (at {', '.join(str(p) for p in cys)}) could form a disulfide")

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


# How much each term moves the confidence. Propensity carries most of it; the moment is a
# supporting hint and the salt bridges a small bonus, since neither alone makes a helix.
_W_PROPENSITY = 0.75
_W_MOMENT = 0.15
_W_SALT_BRIDGE = 0.10
# On the Kyte-Doolittle scale a clearly amphipathic helix runs around 2 per residue, while a
# sequence with no face separation sits near 0.5.
_MOMENT_SATURATION = 2.0
_SALT_BRIDGE_SATURATION = 3


def analyze_helix(sequence: str) -> HelixReport:
    """Estimate whether a binder is helical enough to graft, and flag anything blocking."""
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
    # A blocking issue means the shape is wrong no matter how good the averages look.
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


# Substitutions the charge optimizer makes that cost helix stability for reasons the average
# propensity does not capture. Each returns a message when it applies.
def compare_for_variant(original: str, variant: str) -> list[str]:
    """
    Point out where a charge-optimized variant is worse for the helix than the original.

    Neutralizing D/E to A usually helps, alanine is the best helix former there is. The
    losses are positional: an N-cap that stops capping, a salt bridge that no longer pairs,
    and the helix macrodipole losing the negative charge that was stabilizing its N-terminus.
    """
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
            "from that end even though the new residue is a better helix former mid-chain."
        )

    lost = set(_salt_bridges(orig)) - set(_salt_bridges(var))
    if lost:
        pairs = ", ".join(f"{i}-{j}" for i, j in sorted(lost))
        notes.append(f"Salt bridges lost at {pairs}, each was worth roughly 0.1-0.5 kcal/mol of helix stability.")

    # The helix dipole runs positive at the N-terminus, so acidic residues in the first turn
    # are stabilizing. Neutralizing them there gives that back.
    dipole_lost = [p for p in changed if p <= 4 and orig[p - 1] in _ACIDIC and var[p - 1] not in _ACIDIC]
    if dipole_lost:
        positions = ", ".join(str(p) for p in dipole_lost)
        notes.append(
            f"Acidic residues at {positions} were neutralized, which removes a favorable "
            "interaction with the helix dipole near the N-terminus."
        )

    return notes
