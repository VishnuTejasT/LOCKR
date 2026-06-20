"""POST /verify-assembly -- thin wrapper over assembly.py's structural
checklist and its liability.py bridge function.

candidate_variants/binder_offset exist because a variant /suggest proposes is
binder-local (see docs/README.md step 5) -- filter_safe_variants takes
positions at face value, so this route is the one place that has to apply
`offset = binder_start_in_assembly - 1` before checking a suggestion against
a ProtectedRegion.
"""

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import assembly
from lockr.engine.models import GraftSpec, LatchWindow, ProtectedRegion, VariantSuggestion

from ..schemas.assembly import (
    CandidateVariant, CheckResult, RejectedVariant, VerifyAssemblyRequest,
    VerifyAssemblyResponse, VariantScreen,
)

router = APIRouter()


def _to_variant_suggestion(candidate: CandidateVariant, offset: int) -> VariantSuggestion:
    mutations = [f"{s.from_}{s.position + offset}{s.to}" for s in candidate.substitutions]
    return VariantSuggestion(policy="", sequence=candidate.sequence, mutations=mutations,
                             liability_score=candidate.liability_score,
                             liability_band=candidate.liability_band,
                             K_CK_estimate=candidate.estimated_kck_nm * 1e-9)


@router.post("/verify-assembly", response_model=VerifyAssemblyResponse)
def verify_assembly(request: VerifyAssemblyRequest) -> VerifyAssemblyResponse:
    latch_window = LatchWindow(start=request.latch_window.start, end=request.latch_window.end,
                               expected_length=request.latch_window.expected_length)
    graft_spec = GraftSpec(binder=request.graft_spec.binder, start=request.graft_spec.start,
                           spacer=request.graft_spec.spacer, spacer_start=request.graft_spec.spacer_start,
                           linker=request.graft_spec.linker, linker_start=request.graft_spec.linker_start,
                           binder2=request.graft_spec.binder2, binder2_start=request.graft_spec.binder2_start)
    protected_region = ProtectedRegion(motif=request.protected_region.motif,
                                       start=request.protected_region.start,
                                       end=request.protected_region.end,
                                       label=request.protected_region.label)

    result = assembly.verify_full_assembly(request.full_sequence, latch_window, graft_spec,
                                           protected_region, request.expected_total_length)
    checks = [CheckResult(name=c.name, passed=c.passed, detail=c.detail) for c in result.checks]

    variants = None
    if request.candidate_variants:
        by_identity = {}
        engine_variants = []
        for candidate in request.candidate_variants:
            v = _to_variant_suggestion(candidate, request.binder_offset)
            by_identity[id(v)] = candidate
            engine_variants.append(v)

        screened = assembly.filter_safe_variants(engine_variants, protected_region)
        variants = VariantScreen(
            accepted=[by_identity[id(v)] for v in screened.accepted],
            rejected=[RejectedVariant(variant=by_identity[id(v)], reason=reason)
                     for v, reason in screened.rejected],
        )

    return VerifyAssemblyResponse(all_passed=result.all_passed, checks=checks, variants=variants)
