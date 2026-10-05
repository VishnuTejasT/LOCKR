"""Request/response shapes for POST /api/inverse_design and POST /api/mutations.

Concentrations cross this boundary in nM, like every other route; free energies in
kcal/mol; temperature in K. The engine works in M, so the route converts.
"""

from __future__ import annotations

from pydantic import BaseModel, field_validator, model_validator

from lockr.engine.inverse import CONC_PARAMS, PARAMS
from lockr.engine.temperature import T_REF_DEFAULT

from .common import ThermoSpecIn, validate_sequence, validate_temperature, validate_thermo_keys

_FIXED_KEYS = set(PARAMS) | {"pull"}


class Range(BaseModel):
    min: float
    max: float

    @model_validator(mode="after")
    def _ordered(self):
        if self.max <= self.min:
            raise ValueError("max must be > min")
        return self


class InverseDesignRequest(BaseModel):
    solve_for: list[str]
    fixed: dict[str, float]
    bounds: dict[str, Range]
    target_fold_change: float | None = None
    target_fold_change_range: Range | None = None
    temperature_K: float = T_REF_DEFAULT
    thermo: dict[str, ThermoSpecIn] = {}
    rel_tol: float = 0.05

    @field_validator("solve_for")
    @classmethod
    def _known_free_params(cls, v):
        if not v:
            raise ValueError("solve_for must name at least one free parameter")
        if len(set(v)) != len(v):
            raise ValueError("solve_for has duplicates")
        unknown = [n for n in v if n not in PARAMS]
        if unknown:
            raise ValueError(f"cannot solve for {unknown}; choose from {list(PARAMS)}")
        return v

    @field_validator("fixed")
    @classmethod
    def _known_fixed(cls, v):
        unknown = [k for k in v if k not in _FIXED_KEYS]
        if unknown:
            raise ValueError(f"unknown fixed parameter(s) {unknown}; use {sorted(_FIXED_KEYS)}")
        if "pull" not in v:
            raise ValueError("fixed must include pull (the cage-latch allosteric coupling)")
        if v["pull"] < 0:
            raise ValueError("pull must be >= 0")
        return v

    @field_validator("bounds")
    @classmethod
    def _known_bounds(cls, v):
        unknown = [k for k in v if k not in PARAMS]
        if unknown:
            raise ValueError(f"cannot bound {unknown}; choose from {list(PARAMS)}")
        for name in CONC_PARAMS:
            if name in v and v[name].min <= 0:
                raise ValueError(f"bounds for {name} must be > 0 nM")
        return v

    @field_validator("temperature_K")
    @classmethod
    def _temperature_in_range(cls, v):
        return validate_temperature(v)

    @field_validator("thermo")
    @classmethod
    def _thermo_keys(cls, v):
        return validate_thermo_keys(v)

    @field_validator("rel_tol")
    @classmethod
    def _tolerance_fraction(cls, v):
        if not 0 < v < 1:
            raise ValueError("rel_tol must be between 0 and 1")
        return v

    @model_validator(mode="after")
    def _one_target_and_full_bounds(self):
        if (self.target_fold_change is None) == (self.target_fold_change_range is None):
            raise ValueError("give exactly one of target_fold_change or target_fold_change_range")
        if self.target_fold_change is not None and self.target_fold_change <= 0:
            raise ValueError("target_fold_change must be > 0")
        if self.target_fold_change_range is not None and self.target_fold_change_range.min <= 0:
            raise ValueError("target_fold_change_range min must be > 0")
        missing = [n for n in self.solve_for if n not in self.bounds]
        if missing:
            raise ValueError(f"no bounds given for {missing}")
        overlap = [n for n in self.solve_for if n in self.fixed]
        if overlap:
            raise ValueError(f"{overlap} appear in both solve_for and fixed")
        return self


class SolutionOut(BaseModel):
    values: dict[str, float]
    fold_change: float
    on_spec: bool
    sensitivity: dict[str, float]
    tolerance: dict[str, Range]


class ClosestOut(BaseModel):
    values: dict[str, float]
    fold_change: float
    binding_constraint: str
    note: str


class InverseDesignResponse(BaseModel):
    mode: str
    feasible: bool
    message: str
    temperature_K: float
    thermo_mode: str
    solutions: list[SolutionOut] = []
    region: list[list[list[float]]] = []
    region_axes: list[str] | None = None
    robust_point: SolutionOut | None = None
    closest: ClosestOut | None = None
    warnings: list[str] = []
    units: dict[str, str] = {}


class MutationMapRequest(BaseModel):
    latch_sequence: str
    current_dG_latch_cage: float
    target_dG_latch_cage: float
    interface_positions: list[int]
    top_n: int = 10
    allow_double: bool = True

    @field_validator("latch_sequence")
    @classmethod
    def _sequence(cls, v):
        return validate_sequence(v)

    @field_validator("interface_positions")
    @classmethod
    def _positions(cls, v):
        if not v:
            raise ValueError("interface_positions must not be empty")
        if any(p < 1 for p in v):
            raise ValueError("interface positions are 1-based")
        return v

    @field_validator("top_n")
    @classmethod
    def _top_n(cls, v):
        if not 1 <= v <= 200:
            raise ValueError("top_n must be between 1 and 200")
        return v


class MutationCandidateOut(BaseModel):
    mutation: str
    positions: list[int]
    ddG_estimate: float
    resulting_dG: float
    error_vs_required: float
    confidence: str
    notes: list[str] = []


class MutationMapResponse(BaseModel):
    current_dG: float
    target_dG: float
    required_ddG: float
    interface_positions: list[int]
    candidates: list[MutationCandidateOut]
    disclaimer: str
    warnings: list[str] = []
