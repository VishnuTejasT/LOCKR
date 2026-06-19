"""Sequence-level structural/assembly checks -- position bookkeeping, not folding.

Generalizes my own ECLIPSE pipeline's manual verification (Complete
Documentation Script 4's sequence validation, Script 6's six-point check) into
checks that work for any protected motif, any latch window, any graft. No
structure prediction here -- this is the "did I actually build the sequence I
think I built" layer that sits ALONGSIDE the thermo/charge/liability engine,
not inside it. liability.py never imports anything from this module; the one
cross-module link (filter_safe_variants) runs the other way.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .models import GraftSpec, LatchWindow, ProtectedRegion, VariantSuggestion


def _segment(sequence: str, start: int, length: int) -> str:
    # 1-indexed, inclusive -- matches BinderSequence.residues() elsewhere.
    return sequence[start - 1:start - 1 + length]


@dataclass
class ProtectedRegionCheck:
    intact: bool
    found_sequence: str
    mismatch_positions: list[int] = field(default_factory=list)


def check_protected_region(full_sequence: str, protected_motif: str,
                           start: int, end: int) -> ProtectedRegionCheck:
    found = full_sequence[start - 1:end]
    if found == protected_motif:
        return ProtectedRegionCheck(True, found, [])

    mismatches = [start + i for i in range(min(len(found), len(protected_motif)))
                 if found[i] != protected_motif[i]]
    # A length mismatch means the window itself is wrong, not just a residue --
    # flag every position past the shorter string too, not just substitutions.
    if len(found) != len(protected_motif):
        mismatches += list(range(start + min(len(found), len(protected_motif)),
                                 start + max(len(found), len(protected_motif))))
    return ProtectedRegionCheck(False, found, mismatches)


@dataclass
class OverlapCheck:
    overlap: bool
    overlapping_positions: list[int] = field(default_factory=list)


def _graft_segments(graft_spec: GraftSpec):
    # binder/linker/binder2 are the things actually being inserted; spacer is
    # pre-existing scaffold, not part of the graft, so it's excluded here.
    segments = [(graft_spec.start, graft_spec.start + len(graft_spec.binder) - 1)]
    if graft_spec.linker is not None:
        segments.append((graft_spec.linker_start,
                         graft_spec.linker_start + len(graft_spec.linker) - 1))
    if graft_spec.binder2 is not None:
        segments.append((graft_spec.binder2_start,
                         graft_spec.binder2_start + len(graft_spec.binder2) - 1))
    return segments


def check_graft_overlap(graft_spec: GraftSpec, protected_region: ProtectedRegion) -> OverlapCheck:
    positions = []
    for seg_start, seg_end in _graft_segments(graft_spec):
        lo = max(seg_start, protected_region.start)
        hi = min(seg_end, protected_region.end)
        if lo <= hi:
            positions.extend(range(lo, hi + 1))
    return OverlapCheck(bool(positions), positions)


@dataclass
class LatchFitCheck:
    fits: bool
    used_length: int
    available_length: int
    slack: int


def check_latch_fit(graft_spec: GraftSpec, latch_window: LatchWindow) -> LatchFitCheck:
    used = len(graft_spec.binder)
    if graft_spec.linker is not None:
        used += len(graft_spec.linker)
    if graft_spec.binder2 is not None:
        used += len(graft_spec.binder2)
    available = latch_window.end - latch_window.start + 1
    slack = available - used
    return LatchFitCheck(slack >= 0, used, available, slack)


@dataclass
class AssemblyCheck:
    name: str
    passed: bool
    detail: str


@dataclass
class AssemblyResult:
    checks: list[AssemblyCheck] = field(default_factory=list)

    @property
    def all_passed(self) -> bool:
        return all(c.passed for c in self.checks)


def verify_full_assembly(full_sequence: str, latch_window: LatchWindow, graft_spec: GraftSpec,
                         protected_region: ProtectedRegion,
                         expected_total_length: int | None = None) -> AssemblyResult:
    """All-in-one checklist, generalizing my Script 6 six-point pattern.

    Checks only show up for segments the caller actually supplied -- a v1.0
    single binder gets no linker/binder2 rows, a v2.2 tandem gets both. Nothing
    here hardcodes a spacer value or a motif; spacer/binder content come from
    graft_spec, the motif from protected_region.
    """
    checks = []

    if expected_total_length is not None:
        ok = len(full_sequence) == expected_total_length
        checks.append(AssemblyCheck("overall_length", ok,
                      f"length {len(full_sequence)} (expected {expected_total_length})"))

    pr = check_protected_region(full_sequence, protected_region.motif,
                                protected_region.start, protected_region.end)
    checks.append(AssemblyCheck("protected_region_intact", pr.intact,
                  f"found {pr.found_sequence!r} at {protected_region.start}-{protected_region.end}"
                  + ("" if pr.intact else f", mismatches at {pr.mismatch_positions}")))

    overlap = check_graft_overlap(graft_spec, protected_region)
    checks.append(AssemblyCheck("graft_no_overlap", not overlap.overlap,
                  "no overlap" if not overlap.overlap
                  else f"overlaps protected region at {overlap.overlapping_positions}"))

    fit = check_latch_fit(graft_spec, latch_window)
    checks.append(AssemblyCheck("latch_fit", fit.fits,
                  f"used {fit.used_length}/{fit.available_length}, slack {fit.slack}"))

    if graft_spec.spacer is not None:
        found = _segment(full_sequence, graft_spec.spacer_start, len(graft_spec.spacer))
        checks.append(AssemblyCheck("spacer_intact", found == graft_spec.spacer,
                      f"found {found!r}, expected {graft_spec.spacer!r}"))

    found_b1 = _segment(full_sequence, graft_spec.start, len(graft_spec.binder))
    checks.append(AssemblyCheck("binder1_intact", found_b1 == graft_spec.binder,
                  f"found {found_b1!r}, expected {graft_spec.binder!r}"))

    if graft_spec.linker is not None:
        found_l = _segment(full_sequence, graft_spec.linker_start, len(graft_spec.linker))
        checks.append(AssemblyCheck("linker_intact", found_l == graft_spec.linker,
                      f"found {found_l!r}, expected {graft_spec.linker!r}"))

    if graft_spec.binder2 is not None:
        found_b2 = _segment(full_sequence, graft_spec.binder2_start, len(graft_spec.binder2))
        checks.append(AssemblyCheck("binder2_intact", found_b2 == graft_spec.binder2,
                      f"found {found_b2!r}, expected {graft_spec.binder2!r}"))

    return AssemblyResult(checks)


def _mutation_position(mutation: str) -> int:
    # liability.py writes mutations as "{old}{pos}{new}", e.g. "D4A" -- digits
    # are always the position since residue codes are letters.
    return int("".join(ch for ch in mutation if ch.isdigit()))


@dataclass
class FilteredVariants:
    accepted: list[VariantSuggestion] = field(default_factory=list)
    rejected: list[tuple[VariantSuggestion, str]] = field(default_factory=list)


def filter_safe_variants(suggested_variants: list[VariantSuggestion],
                         protected_region: ProtectedRegion) -> FilteredVariants:
    """Keep liability.py's variant suggester from ever proposing a substitution
    inside a protected region it has no idea exists.

    One-way dependency: this reads VariantSuggestion's shape from models.py,
    but liability.py never imports assembly.py -- the liability scanner stays
    completely unaware that protected regions are a concept.
    """
    accepted, rejected = [], []
    for variant in suggested_variants:
        hit = next((m for m in variant.mutations
                   if protected_region.start <= _mutation_position(m) <= protected_region.end), None)
        if hit is None:
            accepted.append(variant)
        else:
            pos = _mutation_position(hit)
            rejected.append((variant, f"substitution at position {pos} falls inside protected region"))
    return FilteredVariants(accepted, rejected)
