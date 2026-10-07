"""POST /api/temperature/window: nM -> molar, Kd -> dG, and the 422 refusals.

Fold change and the van't Hoff shift both come from the engine; this layer only
converts units and shapes the response. dCp is ignored by the model behind this
endpoint (van't Hoff only) and the response notes say so on every call.
"""

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import temperature_window as tw
from lockr.engine.temperature_window import c_to_k

from ..errors import ApiError
from ..schemas.temperature_window import (ECLIPSE_PRESET, ECLIPSE_PRESET_PLACEHOLDER_FIELDS,
                                          NM_TO_M, UNITS, EchoOut, SystemWindowOut,
                                          TemperatureWindowRequest,
                                          TemperatureWindowResponse, validate_request)

router = APIRouter(prefix="/api")


@router.get("/temperature/window/preset")
def temperature_window_preset() -> dict:
    """The ECLIPSE example request. Placeholder numbers, served so the UI holds no copy."""
    return {
        "name": "ECLIPSE example",
        "request": ECLIPSE_PRESET,
        "placeholder": True,
        "placeholder_fields": list(ECLIPSE_PRESET_PLACEHOLDER_FIELDS),
        "badge": "Example values, replace with real data",
    }


def _to_engine(system, t_ref_c: float) -> dict:
    """One request construct as the engine wants it: molar concentrations, dG not Kd."""
    binder = system.binder
    dG_binder = binder.dG_target_binder
    if dG_binder is None:
        # exactly one of the two is set; validate_request already enforced that
        dG_binder = float(tw.dG_from_kd_nM(binder.kd_nM, c_to_k(t_ref_c)))
    return {
        "name": system.name,
        "dG_latch_cage": system.dG_latch_cage,
        "dH_latch_cage": system.dH_latch_cage,
        "dG_key_cage": system.dG_key_cage,
        "dH_key_cage": system.dH_key_cage,
        "key_conc": None if system.key_conc is None else system.key_conc * NM_TO_M,
        "target_conc": None if system.target_conc is None else system.target_conc * NM_TO_M,
        "t_max_stability_c": system.t_max_stability_c,
        "pull": system.pull,
        "binder": {
            "name": binder.name,
            "dG_target_binder": dG_binder,
            "dG_target_binder_sd": binder.dG_target_binder_sd,
            "dH_target_binder": binder.dH_target_binder,
            "tm_c": binder.tm_c,
            "graft_dG_latch_cage_shift": binder.graft_dG_latch_cage_shift,
        },
    }


@router.post("/temperature/window", response_model=TemperatureWindowResponse)
def temperature_window(request: TemperatureWindowRequest) -> TemperatureWindowResponse:
    validate_request(request)

    priors = None
    if request.dH_priors is not None:
        priors = {k: (None if v is None else v.model_dump())
                  for k, v in request.dH_priors.items()}

    try:
        result = tw.scan(
            [_to_engine(s, request.t_ref_c) for s in request.systems],
            t_min_c=request.t_min_c, t_max_c=request.t_max_c, t_step_c=request.t_step_c,
            t_ref_c=request.t_ref_c, fc_threshold=request.fc_threshold,
            confidence=request.confidence, n_samples=request.n_samples,
            dH_priors=priors, seed=request.seed,
            luciferase_activity=request.luciferase_activity)
    except ValueError as exc:
        # the engine guards the same rules; keep the endpoint's promised status code
        raise ApiError("VALIDATION_ERROR", str(exc), None, status_code=422) from exc

    return TemperatureWindowResponse(
        temperatures_c=result["temperatures_c"],
        systems=[SystemWindowOut(**s) for s in result["systems"]],
        seed=result["seed"],
        echo=EchoOut(
            t_min_c=request.t_min_c, t_max_c=request.t_max_c, t_step_c=request.t_step_c,
            t_ref_c=result["t_ref_c"], fc_threshold=result["fc_threshold"],
            confidence=result["confidence"], n_samples=result["n_samples"],
            dH_priors=result["dH_priors"],
            luciferase_activity=request.luciferase_activity),
        notes=result["notes"],
        units=UNITS)
