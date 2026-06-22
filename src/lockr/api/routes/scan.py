"""POST /scan and POST /suggest -- both thin wrappers over liability.py.

preserve_positions threads straight through to liability.py's own parameter
(the soft target-interface tradeoff, e.g. PfLDH contacts like [1,2,11,12,15])
-- it defaults to [] so sensitive_window alone still does the work when a
caller doesn't have target-contact positions to protect.
"""

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import liability

from ..errors import ApiError
from ..schemas.scan import (
    AcidicResidue, HelixFlag, KckPenalty, PerPosition, ScanRequest, ScanResponse,
    ScanResultItem, SuggestRequest, SuggestResponse, SuggestedVariant, Substitution,
)

router = APIRouter()


def _validate_window_and_preserve(window, preserve_positions: list[int], length: int,
                                  window_field: str, preserve_field: str) -> None:
    # Unlike Window.clamped (used to derive the *effective* window for scoring),
    # this is the spec 8 "out of range -> reject" gate -- callers must send a
    # window/preserve_positions that's actually valid for the sequence they sent,
    # not rely on us silently fixing it up for them.
    if not (1 <= window.start <= length):
        raise ApiError("VALIDATION_ERROR", f"window start must be within 1-{length}", field=f"{window_field}.start")
    if not (1 <= window.end <= length):
        raise ApiError("VALIDATION_ERROR", f"window end must be within 1-{length}", field=f"{window_field}.end")
    for pos in preserve_positions:
        if not (1 <= pos <= length):
            raise ApiError("VALIDATION_ERROR", f"preserve_positions entry {pos} must be within 1-{length}",
                          field=preserve_field)

# Verbatim from spec 8 -- the one server-side warning condition the spec
# actually documents (the live frontend check is just an early heads-up on
# the typed-in length; this is the authoritative one tied to the sequence
# that was actually scanned).
_LONG_SEQUENCE_THRESHOLD = 200
_LONG_SEQUENCE_WARNING = "long sequence — liability model tuned for peptide-scale binders"

# Verbatim from spec 9.2 -- UI copy, not engine output.
_KCK_NOTES = {
    "Low": "There are not any acidic residues in the sensitive region, so K_CK "
           "affinity should be preserved.",
    "Moderate": "There are some acidic residues in the sensitive region, so K_CK may be partially "
                "weakened. Look over the flagged residues.",
    "High": "There are a lot of acidic residues in the sensitive region, so K_CK affinity is likely to be significantly "
            "hindered. Strongly consider the recommended charge-optimized "
            "variant.",
}


def _parse_mutation(mutation: str) -> Substitution:
    # "D4A" -> from D, position 4, to A.
    from_aa, rest = mutation[0], mutation[1:]
    pos = int(rest[:-1])
    to_aa = rest[-1]
    return Substitution(position=pos, **{"from": from_aa}, to=to_aa)


def _scan_one(sequence: str, start: int, end: int, ph: float, policy: str,
              preserve_positions: list[int]) -> ScanResultItem:
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
    suggested = [SuggestedVariant(
        sequence=variant.sequence,
        substitutions=[_parse_mutation(m) for m in variant.mutations],
        liability_score=variant.liability_score,
        liability_band=variant.liability_band,
        estimated_kck_nm=variant.K_CK_estimate * 1e9,
    )]

    warnings = []
    if len(sequence) > _LONG_SEQUENCE_THRESHOLD:
        warnings.append(_LONG_SEQUENCE_WARNING)

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
                                         note=_KCK_NOTES[windowed.liability_band]),
        per_position=per_position,
        helix_flags=helix_flags,
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
    suggested = SuggestedVariant(
        sequence=variant.sequence,
        substitutions=[_parse_mutation(m) for m in variant.mutations],
        liability_score=variant.liability_score,
        liability_band=variant.liability_band,
        estimated_kck_nm=variant.K_CK_estimate * 1e9,
    )
    return SuggestResponse(suggested_variants=[suggested])
