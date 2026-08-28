

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import helix as helix_engine
from lockr.engine import liability

from ..errors import ApiError
from ..schemas.scan import (
    AcidicResidue, CyclizationOut, HelixCheck, HelixFlag, HelixIssueOut, KckPenalty,
    PerPosition, ScanRequest, ScanResponse, ScanResultItem, SuggestRequest, SuggestResponse,
    SuggestedVariant, Substitution,
)

router = APIRouter()


def _validate_window_and_preserve(window, preserve_positions: list[int], length: int,
                                  window_field: str, preserve_field: str) -> None:
    if not (1 <= window.start <= length):
        raise ApiError("VALIDATION_ERROR", f"Your window start must be within 1-{length}", field=f"{window_field}.start")
    if not (1 <= window.end <= length):
        raise ApiError("VALIDATION_ERROR", f"Your window end must be within 1-{length}", field=f"{window_field}.end")
    for pos in preserve_positions:
        if not (1 <= pos <= length):
            raise ApiError("VALIDATION_ERROR", f"Your preserve_positions entry {pos} must be within 1-{length}",
                          field=preserve_field)


_LONG_SEQUENCE_THRESHOLD = 200
_LONG_SEQUENCE_WARNING = "long sequence: liability model tuned for peptide-scale binders"

_KCK_NOTE_NO_ACIDIC = ("There are not any acidic residues in the sensitive region, so K_CK "
                      "affinity should be preserved.")
_KCK_NOTES = {
    "Low": "A few acidic residues are in the sensitive region, but not enough to significantly "
           "weaken K_CK, affinity should be mostly preserved. Still worth a look at the flagged residues.",
    "Moderate": "There are some acidic residues in the sensitive region, so K_CK may be partially "
                "weakened. Look over the flagged residues.",
    "High": "There are a lot of acidic residues in the sensitive region, so K_CK affinity is likely to be significantly "
            "hindered. Strongly consider the recommended charge-optimized "
            "variant.",
}


def _kck_note(band: str, has_acidic_residues: bool) -> str:
    if band == "Low" and not has_acidic_residues:
        return _KCK_NOTE_NO_ACIDIC
    return _KCK_NOTES[band]


_HELIX_BLOCKED_WARNING = ("This binder does not look graftable as a helix, see the structure check. "
                          "The charge results below are still valid, but grafting will not work "
                          "until the shape issue is resolved.")
_HELIX_WEAK_WARNING = ("Helix confidence is low. The lucCage latch can template a helix, so this "
                       "is not fatal, but confirm the binder is helical before trusting the graft.")
_CYCLIC_WARNING = ("This sequence has the residues for a cyclic or stapled peptide. If yours is "
                   "cyclized, it cannot be grafted into the latch as a linear segment.")

_NO_PRESERVE_WARNING = ("No preserve_positions were set, so this suggestion may have changed residues "
                        "needed for target binding, double-check before using it.")


def _helix_check(report) -> HelixCheck:
    return HelixCheck(
        helix_confidence=report.helix_confidence,
        band=report.band,
        mean_propensity=report.mean_propensity,
        hydrophobic_moment=report.hydrophobic_moment,
        salt_bridges=[[i, j] for i, j in report.salt_bridges],
        issues=[HelixIssueOut(position=i.position, severity=i.severity, kind=i.kind, message=i.message)
                for i in report.issues],
        cyclization=CyclizationOut(
            possibly_cyclic=report.cyclization.possibly_cyclic,
            cysteine_positions=report.cyclization.cysteine_positions,
            signals=report.cyclization.signals,
        ),
        graft_blocked=report.graft_blocked,
    )


def _parse_mutation(mutation: str) -> Substitution:
    # "D4A" -> from D, position 4, to A.
    from_aa, rest = mutation[0], mutation[1:]
    pos = int(rest[:-1])
    to_aa = rest[-1]
    return Substitution(position=pos, **{"from": from_aa}, to=to_aa)


def _scan_one(sequence: str, start: int, end: int, ph: float, policy: str,
              preserve_positions: list[int]) -> ScanResultItem:
    # Shape first: a binder that isn't a helix can't be grafted no matter how clean its charge profile is.
    helix_report = helix_engine.analyze_helix(sequence)

    census = liability.scan_liability(sequence, preserve_positions=preserve_positions, ph=ph)
    windowed = liability.scan_liability(sequence, preserve_positions=preserve_positions, ph=ph,
                                        window=(start, end))
    in_window_positions = {l.position: l.penalty for l in windowed.liabilities}

    acidic_residues = [
        AcidicResidue(position=l.position, residue=l.residue,
                      in_window=l.position in in_window_positions,
                      contribution=in_window_positions.get(l.position, 0.0))
        for l in census.liabilities
    ]
    per_position = [
        PerPosition(position=pos, residue=aa, contribution=in_window_positions.get(pos, 0.0))
        for pos, aa in enumerate(sequence, 1)
    ]
    from lockr.engine.charge import helix_breakers
    helix_flags = [HelixFlag(position=p, issue="internal proline/glycine may break the helix")
                  for p in helix_breakers(sequence)]

    variant = liability.suggest_variant(sequence, preserve_positions=preserve_positions, policy=policy,
                                        window=(start, end))
    # Re-run the structure check: swapping D/E changes helix propensity, so the variant's shape isn't assumed.
    variant_helix = helix_engine.analyze_helix(variant.sequence)
    suggested = [SuggestedVariant(
        sequence=variant.sequence,
        substitutions=[_parse_mutation(m) for m in variant.mutations],
        liability_score=variant.liability_score,
        liability_band=variant.liability_band,
        estimated_kck_nm=variant.K_CK_estimate * 1e9,
        helix=_helix_check(variant_helix),
        helix_delta=variant_helix.helix_confidence - helix_report.helix_confidence,
        helix_warnings=helix_engine.compare_for_variant(sequence, variant.sequence),
    )]

    warnings = []
    if helix_report.graft_blocked:
        warnings.append(_HELIX_BLOCKED_WARNING)
    elif helix_report.band == "unlikely helical":
        warnings.append(_HELIX_WEAK_WARNING)
    if helix_report.cyclization.possibly_cyclic:
        warnings.append(_CYCLIC_WARNING)
    if len(sequence) > _LONG_SEQUENCE_THRESHOLD:
        warnings.append(_LONG_SEQUENCE_WARNING)
    if variant.mutations and not preserve_positions:
        warnings.append(_NO_PRESERVE_WARNING)

    return ScanResultItem(
        id="",
        sequence=sequence,
        length=len(sequence),
        net_charge=windowed.net_charge,
        acidic_residues=acidic_residues,
        liability_score=windowed.liability_score,
        liability_band=windowed.liability_band,
        estimated_kck_nm=windowed.K_CK_estimate * 1e9,
        predicted_kck_penalty=KckPenalty(band=windowed.liability_band,
                                         note=_kck_note(windowed.liability_band, bool(windowed.liabilities))),
        per_position=per_position,
        helix_flags=helix_flags,
        helix=_helix_check(helix_report),
        suggested_variants=suggested,
        warnings=warnings,
    )


@router.post("/scan", response_model=ScanResponse)
def scan(request: ScanRequest) -> ScanResponse:
    results = []
    for item in request.sequences:
        length = len(item.sequence)
        _validate_window_and_preserve(request.sensitive_window, request.preserve_positions, length,
                                      window_field="sensitive_window", preserve_field="preserve_positions")
        start, end = request.sensitive_window.clamped(length)
        result = _scan_one(item.sequence, start, end, request.ph, request.substitution_policy,
                           request.preserve_positions)
        result.id = item.id
        results.append(result)
    return ScanResponse(results=results)


@router.post("/suggest", response_model=SuggestResponse)
def suggest(request: SuggestRequest) -> SuggestResponse:
    length = len(request.sequence)
    _validate_window_and_preserve(request.sensitive_window, request.preserve_positions, length,
                                  window_field="sensitive_window", preserve_field="preserve_positions")
    start, end = request.sensitive_window.clamped(length)
    variant = liability.suggest_variant(request.sequence, preserve_positions=request.preserve_positions,
                                        policy=request.substitution_policy, window=(start, end))
    variant_helix = helix_engine.analyze_helix(variant.sequence)
    suggested = SuggestedVariant(
        sequence=variant.sequence,
        substitutions=[_parse_mutation(m) for m in variant.mutations],
        liability_score=variant.liability_score,
        liability_band=variant.liability_band,
        estimated_kck_nm=variant.K_CK_estimate * 1e9,
        helix=_helix_check(variant_helix),
        helix_delta=variant_helix.helix_confidence - helix_engine.analyze_helix(request.sequence).helix_confidence,
        helix_warnings=helix_engine.compare_for_variant(request.sequence, variant.sequence),
    )
    return SuggestResponse(suggested_variants=[suggested])
