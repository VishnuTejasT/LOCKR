

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
