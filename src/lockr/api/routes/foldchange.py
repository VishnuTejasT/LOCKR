"""POST /foldchange -- thin wrapper over thermo.py.

Only one translation happens at this boundary: nM (wire) -> M (engine) for
every concentration. k_open/pull pass straight through -- there's no
standalone "K_open(ON)" in the engine (it's K_open*(1+pull*theta), a derived
quantity), so the request takes the same two knobs thermo.py actually takes
instead of backing pull out of a synthetic ON/OFF ratio.
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
        "Improve cage-key affinity by lowering K_CK.",
        "Latch (K_open) mutations will not raise fold-change in this scenario.",
    ],
    "K_open-limited": [
        "Engineer the latch for stronger allosteric coupling (raise pull).",
        "This isn't key-limited, so lucKey increase will give minimal returns.",
    ],
    "mixed": [
        "Both the latch (K_open) and the key (lucKey/K_CK) are limiting factors, so improvements to either variable will help!"
    ],
}

_LIMITING_FACTOR = {
    "key-limited": "luckey_over_kck",
    "K_open-limited": "k_open",
    "mixed": "mixed",
}


@router.post("/foldchange", response_model=FoldChangeResponse)
def foldchange(request: FoldChangeRequest) -> FoldChangeResponse:
    params = SensorParams(K_open=request.k_open, K_CK=request.k_ck * _NM_TO_M,
                          lucKey=request.luckey * _NM_TO_M)

    if request.target_conc is not None:
        Kd = request.k_target * _NM_TO_M
        target_conc = request.target_conc * _NM_TO_M
        fc = thermo.fold_change(target_conc, Kd, request.pull, params)
    else:
        fc = thermo.max_fold_change(Kd=1.0, pull=request.pull, params=params)

    regime_result = thermo.diagnose_regime(params, pull=request.pull)

    warnings = []
    if request.pull > 50:
        # no documented ECLIPSE case goes above ~30; this far out is more likely a units slip.
        warnings.append("pull strength of >50 has no documented improvement, so verify your units and expectations.")

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
