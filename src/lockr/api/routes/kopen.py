"""POST /kopen-scenarios, the latch-tuning what-if grid behind the Assembly
tab's K_open Tuning Analysis section. All the arithmetic is thermo.py's
(k_open_from_destab, max_fold_change, diagnose_regime); this route is just the
preset x mutation-count loop and the nM boundary conversion, same pattern as
foldchange.py.
"""

from __future__ import annotations

import math

from fastapi import APIRouter

from lockr.engine import thermo
from lockr.engine.models import SensorParams
from lockr.engine.temperature import shift_constant, weakest_mode

from ..errors import ApiError
from ..schemas.common import build_specs
from ..schemas.kopen import KopenScenario, KopenScenariosRequest, KopenScenariosResponse

router = APIRouter()

_UNDEFINED_RESULT_MESSAGE = "Parameters produce an undefined result, check K_open and K_CK."
_NM_TO_M = 1e-9

# 0.5 / 1.0 / 1.5 kcal/mol per mutation, Langan 2019 + my own v1.1 analysis.
_PRESETS = {"conservative": 0.5, "moderate": 1.0, "optimistic": 1.5}
_HELPS_THRESHOLD = 0.05


def _scenario(k_open: float, params: SensorParams, pull: float, mutations: int,
             preset: str, destab: float | None, baseline_fc: float) -> KopenScenario:
    p = SensorParams(K_open=k_open, K_CK=params.K_CK, lucKey=params.lucKey)
    fc = thermo.max_fold_change(Kd=1.0, pull=pull, params=p)
    regime = thermo.diagnose_regime(p, pull=pull).regime.replace("-", "_")
    rel_change = (fc - baseline_fc) / baseline_fc
    if rel_change > _HELPS_THRESHOLD:
        helps = "yes"
    elif rel_change < -_HELPS_THRESHOLD:
        helps = "slightly_worse"
    else:
        helps = "no"
    return KopenScenario(mutations=mutations, preset=preset, destab_kcal_per_mut=destab,
                         k_open=k_open, fold_change=fc, regime=regime, helps=helps)


@router.post("/kopen-scenarios", response_model=KopenScenariosResponse)
def kopen_scenarios(request: KopenScenariosRequest) -> KopenScenariosResponse:
    specs, warnings = build_specs(request.thermo)
    T = request.temperature_K
    # Identity at the default temperature, so the existing what-if grid is unchanged.
    k_open_current, w_latch = shift_constant(request.k_open_current, T, specs["dG_latch_cage"])
    k_ck_m, w_key = shift_constant(request.k_ck * _NM_TO_M, T, specs["dG_key_cage"])
    warnings += w_latch + w_key

    params = SensorParams(K_open=k_open_current, K_CK=k_ck_m,
                          lucKey=request.luckey * _NM_TO_M)

    try:
        baseline_fc = thermo.max_fold_change(Kd=1.0, pull=request.pull, params=params)
        baseline = _scenario(k_open_current, params, request.pull, 0, "baseline", None, baseline_fc)
        baseline = baseline.model_copy(update={"helps": "unchanged"})

        scenarios = []
        for preset, destab in _PRESETS.items():
            for n in (1, 2, 3):
                k_open_new = thermo.k_open_from_destab(k_open_current, n, destab, params.RT)
                scenarios.append(_scenario(k_open_new, params, request.pull, n, preset, destab, baseline_fc))
    except (ZeroDivisionError, OverflowError, ValueError):
        raise ApiError("UNDEFINED_RESULT", _UNDEFINED_RESULT_MESSAGE)

    check_values = [baseline.fold_change, baseline.k_open] + \
        [v for s in scenarios for v in (s.fold_change, s.k_open)]
    if not all(math.isfinite(v) for v in check_values):
        raise ApiError("UNDEFINED_RESULT", _UNDEFINED_RESULT_MESSAGE)

    crossover_luckey_nm = k_open_current * k_ck_m / _NM_TO_M

    return KopenScenariosResponse(dominance_ratio=params.luckey_ratio,
                                  crossover_luckey_nm=crossover_luckey_nm,
                                  baseline=baseline, scenarios=scenarios,
                                  temperature_K=T,
                                  thermo_mode=weakest_mode([specs[k].mode for k in
                                                            ("dG_latch_cage", "dG_key_cage")]),
                                  warnings=warnings)
