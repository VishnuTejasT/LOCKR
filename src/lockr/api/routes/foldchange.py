"""POST /foldchange, thin wrapper over thermo.py.

Only one translation happens at this boundary: nM (wire) -> M (engine) for
every concentration. k_open/pull pass straight through, there's no
standalone "K_open(ON)" in the engine (it's K_open*(1+pull*theta), a derived
quantity), so the request takes the same two knobs thermo.py actually takes
instead of backing pull out of a synthetic ON/OFF ratio.
"""

from __future__ import annotations

import math

from fastapi import APIRouter

from lockr.engine import thermo
from lockr.engine.models import SensorParams

from ..errors import ApiError
from ..schemas.foldchange import FoldChangeRequest, FoldChangeResponse

router = APIRouter()


_UNDEFINED_RESULT_MESSAGE = "Parameters produce an undefined result, check K_open and K_CK."

_NM_TO_M = 1e-9


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


# Plain-language fold-change quality bands (dimensionless signal ratio).
# These are rules of thumb for iGEM/biosensor work, not hard cutoffs: a
# fold-change is the "with target" signal divided by the "no target" signal,
# so 1x means no response and bigger is easier to read out on a plate reader.
def _quality_and_interpretation(fc: float, max_fc: float, at_saturation: bool) -> tuple[str, str]:
    if fc >= 50:
        quality = "excellent"
        desc = "a large, easy-to-read signal jump"
    elif fc >= 10:
        quality = "good"
        desc = "a clear, usable signal jump"
    elif fc >= 3:
        quality = "modest"
        desc = "a detectable but small signal jump"
    else:
        quality = "weak"
        desc = "barely above background, likely hard to detect on a plate reader"

    when = "at saturating target" if at_saturation else "at the target level you entered"
    sentence = (f"About {fc:.1f}x brighter luminescence when the target is present versus "
                f"absent {when}, {desc}.")
    if not at_saturation and max_fc > fc * 1.05:
        sentence += f" With fully saturating target this design could reach about {max_fc:.1f}x."
    return quality, sentence


@router.post("/foldchange", response_model=FoldChangeResponse)
def foldchange(request: FoldChangeRequest) -> FoldChangeResponse:
    params = SensorParams(K_open=request.k_open, K_CK=request.k_ck * _NM_TO_M,
                          lucKey=request.luckey * _NM_TO_M)

    try:
        at_saturation = request.target_conc is None
        max_fc = thermo.max_fold_change(Kd=1.0, pull=request.pull, params=params)
        k_target_m = request.k_target * _NM_TO_M if request.k_target is not None else None
        if not at_saturation:
            Kd = k_target_m
            target_conc = request.target_conc * _NM_TO_M
            fc = thermo.fold_change(target_conc, Kd, request.pull, params)
        else:
            fc = max_fc
        regime_result = thermo.diagnose_regime(params, pull=request.pull)
        fraction_of_dominance_ratio = fc / params.luckey_ratio
        lod_result = thermo.lod_and_ec50(k_target_m, request.pull, params)
    except (ZeroDivisionError, OverflowError, ValueError):
        raise ApiError("UNDEFINED_RESULT", _UNDEFINED_RESULT_MESSAGE)

    if not all(math.isfinite(v) for v in (fc, max_fc, params.luckey_ratio, fraction_of_dominance_ratio)):
        raise ApiError("UNDEFINED_RESULT", _UNDEFINED_RESULT_MESSAGE)

    def _to_nm(value_m: float | None) -> float | None:
        return None if value_m is None else value_m / _NM_TO_M

    quality, interpretation = _quality_and_interpretation(fc, max_fc, at_saturation)

    warnings = []
    if request.pull > 50:
        # no documented ECLIPSE case goes above ~30; this far out is more likely a units slip.
        warnings.append("pull strength of >50 has no documented improvement, so verify your units and expectations.")

    return FoldChangeResponse(
        fold_change=fc,
        max_fold_change=max_fc,
        quality=quality,
        interpretation=interpretation,
        dominance_ratio=params.luckey_ratio,
        fraction_of_dominance_ratio=fraction_of_dominance_ratio,
        regime=regime_result.regime.replace("-", "_"),
        limiting_factor=_LIMITING_FACTOR[regime_result.regime],
        verdict=regime_result.verdict,
        recommendations=_RECOMMENDATIONS[regime_result.regime],
        warnings=warnings,
        lod_2x_nm=_to_nm(lod_result.lod_2x),
        lod_3x_nm=_to_nm(lod_result.lod_3x),
        ec50_nm=_to_nm(lod_result.ec50),
    )
