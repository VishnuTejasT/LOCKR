"""From a required ddG to mutations someone can actually order.

inverse.py answers in kcal/mol ("you need dG_latch_cage = -6.2"). A wet lab
orders DNA, so this ranks single and double substitutions at latch interface
positions by how close their estimated ddG lands to that requirement, reusing
Module A's own scanning machinery (the acidic-liability penalty from
calibration.py and the Chou-Fasman propensities from charge.py) rather than
inventing a second energy model.

These ddG values are heuristic sequence-feature estimates for RANKING candidates,
not predictions. Nothing here is a substitute for Rosetta/FoldX or a plate reader.

Sign convention matches temperature.py: dG_latch_cage is negative for a cage that
holds shut, so a negative ddG means a tighter latch and a lower K_open.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import calibration
from .charge import helix_propensity

# ponytail: two-term heuristic (charge liability + helix propensity). It carries no
# packing, no burial and no context. Upgrade path is graft.py's PyRosetta scoring,
# which is why the confidence flag and the disclaimer exist.
_PROPENSITY_KCAL = 1.0  # kcal/mol per unit Chou-Fasman P_alpha, order-of-magnitude only

_ALPHABET = "ADEFGIKLMNPQRSTVWY"  # C excluded: free cysteine invites disulfide scrambling

# The estimate is not good to better than ~0.1 kcal/mol, so anything inside one bin counts
# as a tie and the fewer-mutations candidate wins it. Otherwise a double mutant that is
# 0.04 kcal/mol "closer" outranks the single mutation a lab would rather order.
_RANKING_BIN_KCAL = 0.1
_ACIDIC = "DE"

ESTIMATE_DISCLAIMER = (
    "ddG values here are HEURISTIC estimates from sequence features (acidic-liability penalty "
    "plus helix propensity), not predictions. They rank candidates for follow-up only. Validate "
    "computationally (Rosetta/FoldX ddG) or experimentally before ordering a large panel."
)


@dataclass
class MutationCandidate:
    mutation: str
    positions: list[int]
    ddG_estimate: float
    resulting_dG: float
    error_vs_required: float
    confidence: str
    notes: list[str] = field(default_factory=list)


@dataclass
class MutationReport:
    current_dG: float
    target_dG: float
    required_ddG: float
    interface_positions: list[int]
    candidates: list[MutationCandidate] = field(default_factory=list)
    disclaimer: str = ESTIMATE_DISCLAIMER
    warnings: list[str] = field(default_factory=list)


def _ddG_single(old: str, new: str, penalty_per_acidic: float) -> tuple[float, list[str], str]:
    """Estimated ddG (kcal/mol) for one substitution, plus notes and confidence."""
    notes = []
    charge_term = 0.0

    if old in _ACIDIC and new not in _ACIDIC:
        # Module A's calibrated liability: an acidic residue at an interface costs binding.
        charge_term = -penalty_per_acidic
        notes.append(f"removes Module A acidic liability {old} ({charge_term:+.2f} kcal/mol)")
    elif new in _ACIDIC and old not in _ACIDIC:
        charge_term = penalty_per_acidic
        notes.append(f"introduces an acidic liability {new} ({charge_term:+.2f} kcal/mol)")

    # helix_propensity of a single residue is that residue's Chou-Fasman P_alpha.
    d_prop = helix_propensity(new) - helix_propensity(old)
    prop_term = -_PROPENSITY_KCAL * d_prop
    if d_prop:
        notes.append(f"helix propensity {'up' if d_prop > 0 else 'down'} "
                     f"({d_prop:+.2f} P_alpha, {prop_term:+.2f} kcal/mol)")
    ddG = charge_term + prop_term

    # Only the charge term is calibrated (calibration.ANCHORS, PfLB-1 -> PfLB-1.5). The
    # propensity term is a scaling guess, so it drags confidence down when it dominates.
    confidence = ("moderate" if charge_term and abs(charge_term) >= 2 * abs(prop_term) else "low")
    if new == "P":
        notes.append("proline breaks the latch helix; propensity alone understates the damage")
    if new == "G":
        notes.append("glycine adds backbone entropy; propensity alone understates the damage")
    return ddG, notes, confidence


def rank_mutations(latch_sequence: str, current_dG: float, target_dG: float,
                   interface_positions, top_n: int = 10, allow_double: bool = True,
                   penalty_per_acidic: float = calibration.PENALTY_PER_ACIDIC) -> MutationReport:
    seq = latch_sequence.strip().upper()
    if not seq:
        raise ValueError("latch_sequence is empty")
    positions = sorted({int(p) for p in interface_positions})
    if not positions:
        raise ValueError("interface_positions is empty")
    off = [p for p in positions if not 1 <= p <= len(seq)]
    if off:
        raise ValueError(f"interface positions {off} are outside the latch sequence (1-{len(seq)})")

    required = target_dG - current_dG
    candidates: list[MutationCandidate] = []

    singles: list[tuple[int, str, str, float, list[str], str]] = []
    for pos in positions:
        old = seq[pos - 1]
        for new in _ALPHABET:
            if new == old:
                continue
            ddG, notes, conf = _ddG_single(old, new, penalty_per_acidic)
            singles.append((pos, old, new, ddG, notes, conf))
            candidates.append(MutationCandidate(
                mutation=f"{old}{pos}{new}", positions=[pos], ddG_estimate=ddG,
                resulting_dG=current_dG + ddG, error_vs_required=abs(ddG - required),
                confidence=conf, notes=notes))

    if allow_double:
        for i, a in enumerate(singles):
            for b in singles[i + 1:]:
                if a[0] == b[0]:
                    continue
                ddG = a[3] + b[3]
                # Additivity is an assumption, not a result, so every double drops to low.
                candidates.append(MutationCandidate(
                    mutation=f"{a[1]}{a[0]}{a[2]}+{b[1]}{b[0]}{b[2]}",
                    positions=sorted([a[0], b[0]]), ddG_estimate=ddG,
                    resulting_dG=current_dG + ddG, error_vs_required=abs(ddG - required),
                    confidence="low",
                    notes=a[4] + b[4] + ["double mutant: ddG assumed additive, often wrong"]))

    candidates.sort(key=lambda c: (int(c.error_vs_required / _RANKING_BIN_KCAL),
                                   len(c.positions), c.error_vs_required))
    best = candidates[:max(1, top_n)]

    warnings = [ESTIMATE_DISCLAIMER]
    reach = max(abs(c.ddG_estimate) for c in candidates)
    if abs(required) > reach:
        warnings.append(f"Required ddG of {required:+.2f} kcal/mol is beyond anything this heuristic "
                        f"can reach at these positions (largest estimate {reach:.2f} kcal/mol). Widen "
                        f"the interface position list, or redesign rather than point-mutate.")
    if best and best[0].error_vs_required > 0.5:
        warnings.append(f"Best candidate is still {best[0].error_vs_required:.2f} kcal/mol off the "
                        f"requirement; treat the list as directional.")
    return MutationReport(current_dG=current_dG, target_dG=target_dG, required_ddG=required,
                          interface_positions=positions, candidates=best, warnings=warnings)
