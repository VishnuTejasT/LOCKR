

from __future__ import annotations

import numpy as np

from fastapi import APIRouter

from lockr.engine import thermo
from lockr.engine.models import SensorParams
from lockr.engine.temperature import shift_constant, weakest_mode

from ..schemas.common import build_specs
from ..schemas.sweep import OperatingPoint, SweepPoint, SweepRequest, SweepResponse

router = APIRouter()

_NM_TO_M = 1e-9


def _fold_change_at(k_ck, k_open, pull, luckey, temperature_K,
                    specs) -> tuple[float, float, list[str]]:
    # At the default temperature both constants come back untouched, so the sweep is unchanged.
    k_open_t, w_latch = shift_constant(k_open, temperature_K, specs["dG_latch_cage"])
    k_ck_m, w_key = shift_constant(k_ck * _NM_TO_M, temperature_K, specs["dG_key_cage"])
    params = SensorParams(K_open=k_open_t, K_CK=k_ck_m, lucKey=luckey * _NM_TO_M)
    fc = thermo.max_fold_change(Kd=1.0, pull=pull, params=params)
    return fc, params.luckey_ratio, w_latch + w_key


@router.post("/sweep", response_model=SweepResponse)
def sweep(request: SweepRequest) -> SweepResponse:
    base = request.base_params
    spec = request.sweep

    xs = (np.logspace(np.log10(spec.min), np.log10(spec.max), spec.steps)
          if spec.scale == "log" else np.linspace(spec.min, spec.max, spec.steps))

    specs, warnings = build_specs(request.thermo)
    T = request.temperature_K

    points = []
    for x in xs:
        values = base.model_dump()
        values[spec.param] = float(x)
        fc, ratio, _ = _fold_change_at(temperature_K=T, specs=specs, **values)
        points.append(SweepPoint(x=float(x), fold_change=fc, dominance_ratio=ratio))

    operating_fc, _, operating_warnings = _fold_change_at(temperature_K=T, specs=specs,
                                                         **base.model_dump())
    operating_x = getattr(base, spec.param)

    return SweepResponse(
        param=spec.param,
        scale=spec.scale,
        points=points,
        operating_point=OperatingPoint(x=operating_x, fold_change=operating_fc),
        temperature_K=T,
        thermo_mode=weakest_mode([specs[k].mode for k in ("dG_latch_cage", "dG_key_cage")]),
        warnings=warnings + operating_warnings,
    )
