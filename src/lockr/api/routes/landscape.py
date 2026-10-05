"""Module D routes: nM/M conversion and the temperature refusal."""

from __future__ import annotations

import numpy as np

from fastapi import APIRouter

from lockr.engine import landscape
from lockr.engine.landscape import NEEDS_TARGET
from lockr.engine.temperature import MODE_ISOTHERMAL, weakest_mode

from ..errors import ApiError
from ..schemas.common import build_specs
from ..schemas.landscape import (NM_PARAMS, UNITS, AxisOut, ContourOut, LandscapeRequest,
                                 LandscapeResponse, RobustOut, SliceRequest, SliceResponse,
                                 SlackOut, TmBandOut)

router = APIRouter(prefix="/api")

_NM_TO_M = 1e-9
_TEMP = "temperature_K"

# Metrics reported in nM, so the engine's M come back in the units the request used.
_NM_METRICS = ("limit_of_detection",)


def _to_engine(name: str, value: float | None) -> float | None:
    if value is None:
        return None
    return value * _NM_TO_M if name in NM_PARAMS else value


def _metric_to_api(metric: str, grid: np.ndarray) -> np.ndarray:
    return grid / _NM_TO_M if metric in _NM_METRICS else grid


def _nullable(grid: np.ndarray) -> list:
    """nan/inf -> None, since bare NaN is invalid JSON."""
    clean = np.where(np.isfinite(grid), grid, None)
    return clean.tolist() if clean.ndim else clean.item()


def _used_spec_keys(k_target) -> list[str]:
    keys = ["dG_latch_cage", "dG_key_cage"]
    return keys + (["dG_target_binder"] if k_target is not None else [])


def _refuse_isothermal(specs, used: list[str], field: str) -> None:
    """Reject a temperature axis with no dH_ref."""
    iso = [k for k in used if specs[k].mode is MODE_ISOTHERMAL]
    if not iso:
        return
    names = ", ".join(iso)
    plural = "are" if len(iso) > 1 else "is"
    raise ApiError(
        "ISOTHERMAL_NO_TEMPERATURE_AXIS",
        f"Cannot put temperature on an axis: {names} {plural} in isothermal mode, meaning "
        f"only a reference value was supplied and nothing is known about how it changes "
        f"with temperature. The surface would be flat at every temperature, which is not a "
        f"prediction, it is the absence of one. Supply dH_ref for {names} in the thermo "
        f"block (and ideally dCp, or buried_nonpolar_A2 to estimate it), or choose a "
        f"different axis parameter.",
        field=field)


def _axis_array(spec, steps: int) -> np.ndarray:
    try:
        values = landscape.axis_values(spec.min, spec.max, steps, spec.scale)
    except ValueError as exc:
        raise ApiError("VALIDATION_ERROR", str(exc), field=f"{spec.param}")
    return np.array([_to_engine(spec.param, float(v)) for v in values])


def _axis_out(spec, values_api: np.ndarray) -> AxisOut:
    return AxisOut(param=spec.param, scale=spec.scale, unit=UNITS[spec.param],
                   values=[float(v) for v in values_api])


def _tm_band_out(specs, used: list[str], axis: str, lo: float, hi: float) -> TmBandOut | None:
    """Earliest Tm among the constants in play."""
    bands = [b for b in (landscape.tm_band(specs[k], lo, hi) for k in used) if b]
    if not bands:
        return None
    band = min(bands, key=lambda b: b["Tm"])
    return TmBandOut(axis=axis, **band)


def _evaluate(values: dict, metric: str, specs) -> tuple[np.ndarray, list[str]]:
    try:
        grid, warns = landscape.evaluate(values, metric, specs)
    except ValueError as exc:
        raise ApiError("VALIDATION_ERROR", str(exc), field="metric")
    except (OverflowError, ZeroDivisionError):
        raise ApiError("UNDEFINED_RESULT",
                       "Parameters produce an undefined result, check the axis ranges "
                       "and fixed values.")
    return _metric_to_api(metric, grid), warns


def _robust_out(found: dict, x_param: str, y_param: str) -> RobustOut:
    xu, yu = UNITS[x_param], UNITS[y_param]
    fmt = lambda s, unit: f"-{s['minus']:.3g} / +{s['plus']:.3g} {unit}"
    summary = (f"Most robust feasible design: {x_param} = {found['x']:.4g} {xu}, "
               f"{y_param} = {found['y']:.4g} {yu}. Slack before the constraints stop "
               f"holding: {fmt(found['x_slack'], xu)} on {x_param}, "
               f"{fmt(found['y_slack'], yu)} on {y_param}.")
    return RobustOut(x=found["x"], y=found["y"], metric_value=found["metric_value"],
                     x_slack=SlackOut(**found["x_slack"]),
                     y_slack=SlackOut(**found["y_slack"]), summary=summary)


@router.post("/landscape", response_model=LandscapeResponse)
def landscape_grid(request: LandscapeRequest) -> LandscapeResponse:
    specs, warnings = build_specs(request.thermo)
    fixed = request.fixed
    used = _used_spec_keys(fixed.k_target)

    if request.metric in NEEDS_TARGET and fixed.k_target is None:
        raise ApiError("VALIDATION_ERROR",
                       f"metric '{request.metric}' needs fixed.k_target, the target binder Kd "
                       f"in nM, because without it the model has no dose response to read.",
                       field="fixed.k_target")

    x_spec, y_spec = request.x_axis, request.y_axis
    for spec, field in ((x_spec, "x_axis.param"), (y_spec, "y_axis.param")):
        if spec.param == _TEMP:
            _refuse_isothermal(specs, used, field)

    nx, ny, cap_warnings = landscape.cap_steps(x_spec.steps, y_spec.steps, request.metric)
    warnings += cap_warnings
    xs, ys = _axis_array(x_spec, nx), _axis_array(y_spec, ny)

    values = {p: _to_engine(p, getattr(fixed, p))
              for p in landscape.REQUIRED_PARAMS + landscape.OPTIONAL_PARAMS}
    values[_TEMP] = request.temperature_K
    values[x_spec.param] = xs[:, None]
    values[y_spec.param] = ys[None, :]

    grid, eval_warnings = _evaluate(values, request.metric, specs)
    warnings += eval_warnings

    # One grid per distinct constraint metric; the plotted metric is reused if constrained.
    grids = {request.metric: grid}
    for constraint in request.constraints:
        if constraint.metric not in grids:
            other, w = _evaluate(values, constraint.metric, specs)
            grids[constraint.metric] = other
            warnings += [m for m in w if m not in warnings]

    constraints = [c.model_dump() for c in request.constraints]
    mask = landscape.feasibility(grids, constraints)

    x_api = xs / _NM_TO_M if x_spec.param in NM_PARAMS else xs
    y_api = ys / _NM_TO_M if y_spec.param in NM_PARAMS else ys
    # API units throughout, so the reported slack matches the axis labels.
    found = landscape.robust_optimum(mask, x_api, y_api, grid)
    levels = request.levels if request.levels is not None else landscape.contour_levels(grid)

    temp_axis = next((a for a, s in (("x", x_spec), ("y", y_spec)) if s.param == _TEMP), None)
    if temp_axis is not None:
        axis_spec = x_spec if temp_axis == "x" else y_spec
        warnings += [w for w in landscape.temperature_axis_warnings(
            {k: specs[k] for k in used}, axis_spec.min, axis_spec.max) if w not in warnings]

    return LandscapeResponse(
        metric=request.metric,
        x=_axis_out(x_spec, x_api), y=_axis_out(y_spec, y_api),
        z=_nullable(grid),
        contours=[ContourOut(**c) for c in landscape.contours(grid, x_api, y_api, levels)],
        feasible=mask.tolist() if mask is not None else None,
        robust_optimum=_robust_out(found, x_spec.param, y_spec.param) if found else None,
        tm_band=(_tm_band_out(specs, used, temp_axis,
                              (x_spec if temp_axis == "x" else y_spec).min,
                              (x_spec if temp_axis == "x" else y_spec).max)
                 if temp_axis else None),
        temperature_K=None if temp_axis else request.temperature_K,
        thermo_mode=weakest_mode([specs[k].mode for k in used]),
        cells=int(grid.size),
        warnings=warnings,
        units=UNITS,
    )


@router.post("/landscape/slice", response_model=SliceResponse)
def landscape_slice(request: SliceRequest) -> SliceResponse:
    specs, warnings = build_specs(request.thermo)
    fixed = request.fixed
    used = _used_spec_keys(fixed.k_target)
    axis = request.axis

    if request.metric in NEEDS_TARGET and fixed.k_target is None:
        raise ApiError("VALIDATION_ERROR",
                       f"metric '{request.metric}' needs fixed.k_target, the target binder Kd "
                       f"in nM.", field="fixed.k_target")
    if axis.param == _TEMP:
        _refuse_isothermal(specs, used, "axis.param")

    xs = _axis_array(axis, axis.steps)
    values = {p: _to_engine(p, getattr(fixed, p))
              for p in landscape.REQUIRED_PARAMS + landscape.OPTIONAL_PARAMS}
    values[_TEMP] = request.temperature_K
    values[axis.param] = xs

    line, eval_warnings = _evaluate(values, request.metric, specs)
    warnings += eval_warnings
    if axis.param == _TEMP:
        warnings += [w for w in landscape.temperature_axis_warnings(
            {k: specs[k] for k in used}, axis.min, axis.max) if w not in warnings]

    x_api = xs / _NM_TO_M if axis.param in NM_PARAMS else xs
    return SliceResponse(
        metric=request.metric,
        axis=_axis_out(axis, x_api),
        values=_nullable(line),
        tm_band=(_tm_band_out(specs, used, "x", axis.min, axis.max)
                 if axis.param == _TEMP else None),
        temperature_K=None if axis.param == _TEMP else request.temperature_K,
        thermo_mode=weakest_mode([specs[k].mode for k in used]),
        warnings=warnings,
        units=UNITS,
    )
