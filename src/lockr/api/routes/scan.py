"""POST /scan and POST /suggest -- both thin wrappers over liability.py.

preserve_positions threads straight through to liability.py's own parameter
(the soft target-interface tradeoff, e.g. PfLDH contacts like [1,2,11,12,15])
-- it defaults to [] so sensitive_window alone still does the work when a
caller doesn't have target-contact positions to protect.
"""

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import liability

from ..schemas.scan import (
    AcidicResidue, HelixFlag, KckPenalty, PerPosition, ScanRequest, ScanResponse,
    ScanResultItem, SuggestRequest, SuggestResponse, SuggestedVariant, Substitution,
)

router = APIRouter()

# Verbatim from spec 9.2 -- UI copy, not engine output.
_KCK_NOTES = {
    "low": "There are not any acidic residues in the sensitive region, so K_CK"
           "affinity should be preserved.",
    "moderate": "There are some acidic residues in the sensitive region, so K_CK may be partially "
                "weakened. Look over the flagged residues.",
    "high": "There are a lot of acidic residues in the sensitive region, so K_CK affinity is likely to be significantly "
            "hindered. Strongly consider the recommended charge-optimized"
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

    return ScanResultItem(
        id="",
        sequence=sequence,
        length=len(sequence),
        net_charge=windowed.net_charge,
        acidic_residues=acidic_residues,
        liability_score=windowed.liability_score,
        liability_band=windowed.liability_band,
        predicted_kck_penalty=KckPenalty(band=windowed.liability_band,
                                         note=_KCK_NOTES[windowed.liability_band]),
        per_position=per_position,
        helix_flags=helix_flags,
        suggested_variants=suggested,
    )


@router.post("/scan", response_model=ScanResponse)
def scan(request: ScanRequest) -> ScanResponse:
    results = []
    for item in request.sequences:
        start, end = request.sensitive_window.clamped(len(item.sequence))
        result = _scan_one(item.sequence, start, end, request.ph, request.substitution_policy,
                           request.preserve_positions)
        result.id = item.id
        results.append(result)
    return ScanResponse(results=results)


@router.post("/suggest", response_model=SuggestResponse)
def suggest(request: SuggestRequest) -> SuggestResponse:
    start, end = request.sensitive_window.clamped(len(request.sequence))
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
