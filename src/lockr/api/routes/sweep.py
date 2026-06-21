"""POST /sweep -- loops thermo.max_fold_change over one varied parameter.

thermo.py has no "sweep an arbitrary named parameter" function (only
scan_dose_response, which is hardcoded to sweep target_conc). So this route
builds one SensorParams per point and calls the same per-point engine
function diagnose_regime/max_fold_change already use -- it's a loop over real
calls, not a new calculation. Whichever of the four params isn't being swept
holds its base_params value.
"""

from __future__ import annotations

import numpy as np

from fastapi import APIRouter

from lockr.engine import thermo
from lockr.engine.models import SensorParams

from ..schemas.sweep import OperatingPoint, SweepPoint, SweepRequest, SweepResponse

router = APIRouter()

_NM_TO_M = 1e-9


def _fold_change_at(k_ck, k_open, pull, luckey) -> tuple[float, float]:
    params = SensorParams(K_open=k_open, K_CK=k_ck * _NM_TO_M, lucKey=luckey * _NM_TO_M)
    fc = thermo.max_fold_change(Kd=1.0, pull=pull, params=params)
    return fc, params.luckey_ratio


@router.post("/sweep", response_model=SweepResponse)
def sweep(request: SweepRequest) -> SweepResponse:
    base = request.base_params
    spec = request.sweep

    xs = (np.logspace(np.log10(spec.min), np.log10(spec.max), spec.steps)
          if spec.scale == "log" else np.linspace(spec.min, spec.max, spec.steps))

    points = []
    for x in xs:
        values = base.model_dump()
        values[spec.param] = float(x)
        fc, ratio = _fold_change_at(**values)
        points.append(SweepPoint(x=float(x), fold_change=fc, dominance_ratio=ratio))

    operating_fc, _ = _fold_change_at(**base.model_dump())
    operating_x = getattr(base, spec.param)

    return SweepResponse(
        param=spec.param,
        scale=spec.scale,
        points=points,
        operating_point=OperatingPoint(x=operating_x, fold_change=operating_fc),
    )
