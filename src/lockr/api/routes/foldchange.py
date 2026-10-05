
from __future__ import annotations

import math

from fastapi import APIRouter

from lockr.engine import thermo
from lockr.engine.models import RT_37C, SensorParams
from lockr.engine.temperature import R_KCAL, T_REF_DEFAULT, shift_constant, weakest_mode

from ..errors import ApiError
from ..schemas.common import build_specs
from ..schemas.foldchange import FoldChangeRequest, FoldChangeResponse

router = APIRouter()


_UNDEFINED_RESULT_MESSAGE = "Parameters produce an undefined result, check K_open and K_CK."

_NM_TO_M = 1e-9


_LIMITING_FACTOR = {
    "key-limited": "luckey_over_kck",
    "K_open-limited": "k_open",
    "mixed": "mixed",
}



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
    T = request.temperature_K
    specs, thermo_notes = build_specs(request.thermo)
    # At the default 298.15 K, and in isothermal mode, every constant comes back untouched.
    k_open, w_latch = shift_constant(request.k_open, T, specs["dG_latch_cage"])
    k_ck_m, w_key = shift_constant(request.k_ck * _NM_TO_M, T, specs["dG_key_cage"])
    k_target_m, w_target = shift_constant(
        request.k_target * _NM_TO_M if request.k_target is not None else None, T,
        specs["dG_target_binder"])
    used = ["dG_latch_cage", "dG_key_cage"] + (["dG_target_binder"] if k_target_m else [])
    thermo_notes += w_latch + w_key + w_target

    # RT only feeds the dG helpers, not the fold change, but keep the scaffold default at
    # T_ref so nothing about the default request changes.
    params = SensorParams(K_open=k_open, K_CK=k_ck_m, lucKey=request.luckey * _NM_TO_M,
                          RT=RT_37C if T == T_REF_DEFAULT else R_KCAL * T)

    try:
        at_saturation = request.target_conc is None
        max_fc = thermo.max_fold_change(Kd=1.0, pull=request.pull, params=params)
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

    warnings = list(thermo_notes)
    if request.pull > 50:
        # no documented ECLIPSE case goes above ~30; this far out is more likely a units slip.
        warnings.append("pull strength of >50 has no documented improvement, so verify your units and expectations.")

    return FoldChangeResponse(
        fold_change=fc,
        temperature_K=T,
        thermo_mode=weakest_mode([specs[k].mode for k in used]),
        max_fold_change=max_fc,
        quality=quality,
        interpretation=interpretation,
        dominance_ratio=params.luckey_ratio,
        fraction_of_dominance_ratio=fraction_of_dominance_ratio,
        regime=regime_result.regime.replace("-", "_"),
        limiting_factor=_LIMITING_FACTOR[regime_result.regime],
        verdict=regime_result.verdict,
        recommendations=regime_result.recommendations,
        warnings=warnings,
        lod_2x_nm=_to_nm(lod_result.lod_2x),
        lod_3x_nm=_to_nm(lod_result.lod_3x),
        ec50_nm=_to_nm(lod_result.ec50),
    )
