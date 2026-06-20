"""POST /foldchange -- thin wrapper over thermo.py.

Two translations happen at this boundary, neither of which is new science:
- nM (wire) -> M (engine) for every concentration.
- the engine has no "pull" field in the request; pull is derived from
  k_open_on/k_open_off (pull = K_open(ON)/K_open(OFF) - 1), since
  k_open_eff(K_open, pull, theta=1) == K_open(ON) by definition.
"""

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import thermo
from lockr.engine.models import SensorParams

from ..schemas.foldchange import FoldChangeRequest, FoldChangeResponse

router = APIRouter()

_NM_TO_M = 1e-9

# Verbatim from spec 9.3 -- canned per-regime copy, not engine output.
_RECOMMENDATIONS = {
    "key-limited": [
        "Increase lucKey to raise the dominance ratio.",
        "Improve cage-key affinity (lower K_CK).",
        "Latch (K_open) mutations will not raise fold-change here.",
    ],
    "K_open-limited": [
        "Engineer the latch to raise K_open(ON) relative to K_open(OFF).",
        "You're not key-limited yet — lucKey increases give diminishing returns.",
    ],
    "mixed": [
        "Both the latch (K_open) and the key (lucKey/K_CK) are constraining "
        "fold-change. Improvements to either will help.",
    ],
}

_LIMITING_FACTOR = {
    "key-limited": "luckey_over_kck",
    "K_open-limited": "k_open",
    "mixed": "mixed",
}


@router.post("/foldchange", response_model=FoldChangeResponse)
def foldchange(request: FoldChangeRequest) -> FoldChangeResponse:
    params = SensorParams(K_open=request.k_open_off, K_CK=request.k_ck * _NM_TO_M,
                          lucKey=request.luckey * _NM_TO_M)
    pull = (request.k_open_on / request.k_open_off) - 1

    if request.target_conc is not None:
        Kd = request.k_target * _NM_TO_M
        target_conc = request.target_conc * _NM_TO_M
        fc = thermo.fold_change(target_conc, Kd, pull, params)
    else:
        fc = thermo.max_fold_change(Kd=1.0, pull=pull, params=params)

    regime_result = thermo.diagnose_regime(params, pull=pull)

    warnings = []
    if request.k_open_on <= request.k_open_off:
        warnings.append("target should stabilize the open state (ON > OFF)")

    return FoldChangeResponse(
        fold_change=fc,
        dominance_ratio=params.luckey_ratio,
        fraction_of_dominance_ratio=fc / params.luckey_ratio,
        regime=regime_result.regime.replace("-", "_"),
        limiting_factor=_LIMITING_FACTOR[regime_result.regime],
        verdict=regime_result.verdict,
        recommendations=_RECOMMENDATIONS[regime_result.regime],
        warnings=warnings,
    )
