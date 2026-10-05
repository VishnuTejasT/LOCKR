"""POST /api/inverse_design and POST /api/mutations, Module C's two entry points.

All the maths is inverse.py's and mutations.py's. This route does the nM/M
boundary conversion, the request-to-spec plumbing, and nothing else, same pattern
as foldchange.py.
"""

from __future__ import annotations

from fastapi import APIRouter

from lockr.engine import inverse, mutations
from lockr.engine.inverse import CONC_PARAMS

from ..errors import ApiError
from ..schemas.common import build_specs
from ..schemas.inverse import (ClosestOut, InverseDesignRequest, InverseDesignResponse,
                               MutationCandidateOut, MutationMapRequest, MutationMapResponse,
                               Range, SolutionOut)

router = APIRouter(prefix="/api")

_NM_TO_M = 1e-9

_UNITS = {"dG_latch_cage": "kcal/mol", "dG_key_cage": "kcal/mol",
          "dG_target_binder": "kcal/mol", "key_conc": "nM", "target_conc": "nM",
          "sensitivity": "d(log fold_change)/d(parameter), per kcal/mol or per nM",
          "temperature_K": "K"}


def _to_engine(name: str, value: float) -> float:
    return value * _NM_TO_M if name in CONC_PARAMS else value


def _to_api(name: str, value: float) -> float:
    return value / _NM_TO_M if name in CONC_PARAMS else value


def _sens_to_api(name: str, value: float) -> float:
    # d(log FC)/d(x in M) -> d(log FC)/d(x in nM)
    return value * _NM_TO_M if name in CONC_PARAMS else value


def _solution_out(solution) -> SolutionOut:
    return SolutionOut(
        values={k: _to_api(k, v) for k, v in solution.values.items()},
        fold_change=solution.fold_change,
        on_spec=solution.on_spec,
        sensitivity={k: _sens_to_api(k, v) for k, v in solution.sensitivity.items()},
        # A degenerate interval means the point is off spec, so there is no drift to report.
        tolerance={k: Range(min=_to_api(k, lo), max=_to_api(k, hi))
                   for k, (lo, hi) in solution.tolerance.items() if hi > lo},
    )


@router.post("/inverse_design", response_model=InverseDesignResponse)
def inverse_design(request: InverseDesignRequest) -> InverseDesignResponse:
    specs, notes = build_specs(request.thermo)
    fixed = {k: _to_engine(k, v) for k, v in request.fixed.items()}
    bounds = {k: (_to_engine(k, r.min), _to_engine(k, r.max)) for k, r in request.bounds.items()}
    target_range = (None if request.target_fold_change_range is None
                    else (request.target_fold_change_range.min, request.target_fold_change_range.max))

    try:
        result = inverse.inverse_design(
            solve_for=request.solve_for, fixed=fixed, bounds=bounds,
            target_fold_change=request.target_fold_change, target_fold_change_range=target_range,
            temperature_K=request.temperature_K, specs=specs, rel_tol=request.rel_tol)
    except ValueError as exc:
        raise ApiError("VALIDATION_ERROR", str(exc))
    except (OverflowError, ZeroDivisionError):
        raise ApiError("UNDEFINED_RESULT",
                       "Parameters produce an undefined result, check the bounds and fixed values.")

    region = []
    if result.region_axes:
        xname, yname = result.region_axes
        region = [[[_to_api(xname, x), _to_api(yname, y)] for x, y in line] for line in result.region]

    closest = None
    if result.closest is not None:
        closest = ClosestOut(
            values={k: _to_api(k, v) for k, v in result.closest.values.items()},
            fold_change=result.closest.fold_change,
            binding_constraint=result.closest.binding_constraint,
            note=result.closest.note)

    return InverseDesignResponse(
        mode=result.mode, feasible=result.feasible, message=result.message,
        temperature_K=result.temperature_K, thermo_mode=result.thermo_mode,
        solutions=[_solution_out(s) for s in result.solutions],
        region=region,
        region_axes=list(result.region_axes) if result.region_axes else None,
        robust_point=_solution_out(result.robust_point) if result.robust_point else None,
        closest=closest,
        warnings=result.warnings + notes,
        units=_UNITS,
    )


@router.post("/mutations", response_model=MutationMapResponse)
def mutation_map(request: MutationMapRequest) -> MutationMapResponse:
    try:
        report = mutations.rank_mutations(
            latch_sequence=request.latch_sequence,
            current_dG=request.current_dG_latch_cage,
            target_dG=request.target_dG_latch_cage,
            interface_positions=request.interface_positions,
            top_n=request.top_n, allow_double=request.allow_double)
    except ValueError as exc:
        raise ApiError("VALIDATION_ERROR", str(exc))

    return MutationMapResponse(
        current_dG=report.current_dG, target_dG=report.target_dG,
        required_ddG=report.required_ddG, interface_positions=report.interface_positions,
        candidates=[MutationCandidateOut(**vars(c)) for c in report.candidates],
        disclaimer=report.disclaimer, warnings=report.warnings)
