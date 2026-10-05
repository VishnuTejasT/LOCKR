"""POST /dosing: back-calculate the Key concentration (and pipetting volume)
needed to hit a target Cage occupancy, using dosing.py's exact quadratic
mass-balance solve rather than thermo.py's excess-Key approximation.
"""

from __future__ import annotations

import math

from fastapi import APIRouter

from lockr.engine import dosing

from ..errors import ApiError
from ..schemas.dosing import DosingRequest, DosingResponse, TitrationPoint

router = APIRouter()

_UNDEFINED_RESULT_MESSAGE = "Parameters produce an undefined result, check K_CK and Cage/Key concentrations."
_NM_TO_M = 1e-9
_UL_TO_L = 1e-6


@router.post("/dosing", response_model=DosingResponse)
def dosing_route(request: DosingRequest) -> DosingResponse:
    try:
        result = dosing.dose_for_target(
            k_ck=request.k_ck * _NM_TO_M,
            total_cage_conc=request.total_cage_conc * _NM_TO_M,
            key_stock_conc=request.key_stock_conc * _NM_TO_M,
            assay_volume=request.assay_volume * _UL_TO_L,
            target_fraction_bound=request.target_fraction_bound,
            target_fold_change=request.target_fold_change,
        )
    except (ZeroDivisionError, OverflowError, ValueError) as exc:
        raise ApiError("UNDEFINED_RESULT", f"{_UNDEFINED_RESULT_MESSAGE} ({exc})")

    check_values = [result.key_conc, result.fraction_bound, result.volume_to_add]
    if not all(math.isfinite(v) for v in check_values):
        raise ApiError("UNDEFINED_RESULT", _UNDEFINED_RESULT_MESSAGE)

    return DosingResponse(
        key_conc_nm=result.key_conc / _NM_TO_M,
        fraction_bound=result.fraction_bound,
        volume_to_add_ul=result.volume_to_add / _UL_TO_L,
        warnings=result.warnings,
        titration=[TitrationPoint(key_conc_nm=k / _NM_TO_M, fraction_bound=fb) for k, fb in result.titration],
    )
