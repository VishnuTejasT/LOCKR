"""Request/response shapes for POST /api/landscape and /api/landscape/slice."""

from __future__ import annotations

from pydantic import BaseModel, field_validator, model_validator

from lockr.engine.landscape import (AXIS_PARAMS, METRICS, OPTIONAL_PARAMS, REQUIRED_PARAMS,
                                    SCALES, _OPS)
from lockr.engine.temperature import T_REF_DEFAULT

from .common import ThermoSpecIn, validate_temperature, validate_thermo_keys

# nM at the API boundary, M inside the engine, same convention as foldchange and sweep.
NM_PARAMS = ("k_ck", "luckey", "k_target", "target_conc")

UNITS = {"k_ck": "nM", "luckey": "nM", "k_target": "nM", "target_conc": "nM",
         "k_open": "dimensionless", "pull": "dimensionless", "temperature_K": "K",
         "fold_change": "x (ratio)", "on_signal": "fraction of cages open and key-bound",
         "off_signal": "fraction of cages open and key-bound, no target",
         "dynamic_range": "difference in open-and-bound fraction",
         "limit_of_detection": "nM", "fraction_open": "fraction of cages open"}


class AxisSpec(BaseModel):
    param: str
    min: float
    max: float
    steps: int = 60
    scale: str = "log"

    @field_validator("param")
    @classmethod
    def _known_param(cls, v):
        if v not in AXIS_PARAMS:
            raise ValueError(f"axis param must be one of {AXIS_PARAMS}")
        return v

    @field_validator("scale")
    @classmethod
    def _known_scale(cls, v):
        if v not in SCALES:
            raise ValueError(f"scale must be one of {SCALES}")
        return v

    @field_validator("steps")
    @classmethod
    def _enough_steps(cls, v):
        if not 2 <= v <= 1000:
            raise ValueError("axis steps must be between 2 and 1000")
        return v

    @model_validator(mode="after")
    def _range_usable(self):
        if self.min >= self.max:
            raise ValueError(f"{self.param} axis min must be < max")
        if self.scale == "log" and self.min <= 0:
            raise ValueError(f"{self.param} axis on a log scale needs min > 0")
        return self


class Constraint(BaseModel):
    metric: str
    op: str
    value: float

    @field_validator("metric")
    @classmethod
    def _known_metric(cls, v):
        if v not in METRICS:
            raise ValueError(f"constraint metric must be one of {METRICS}")
        return v

    @field_validator("op")
    @classmethod
    def _known_op(cls, v):
        if v not in _OPS:
            raise ValueError(f"constraint op must be one of {tuple(_OPS)}")
        return v


class FixedParams(BaseModel):
    """Everything not on an axis. Concentrations in nM."""

    k_ck: float | None = None
    k_open: float | None = None
    luckey: float | None = None
    pull: float | None = None
    k_target: float | None = None
    target_conc: float | None = None

    @field_validator("k_ck", "k_open", "luckey", "k_target", "target_conc")
    @classmethod
    def _positive(cls, v, info):
        if v is not None and v <= 0:
            raise ValueError(f"{info.field_name} must be > 0")
        return v

    @field_validator("pull")
    @classmethod
    def _pull_non_negative(cls, v):
        if v is not None and v < 0:
            raise ValueError("pull must be >= 0")
        return v


def validate_metric(metric: str) -> str:
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}")
    return metric


def _check_coverage(fixed: FixedParams, axis_params: tuple[str, ...]) -> None:
    missing = [p for p in REQUIRED_PARAMS
               if p not in axis_params and getattr(fixed, p) is None]
    if missing:
        raise ValueError(f"fixed is missing {missing}; every parameter not on an axis "
                         f"must be given a value")


class LandscapeRequest(BaseModel):
    x_axis: AxisSpec
    y_axis: AxisSpec
    fixed: FixedParams = FixedParams()
    metric: str = "fold_change"
    temperature_K: float = T_REF_DEFAULT
    thermo: dict[str, ThermoSpecIn] = {}
    constraints: list[Constraint] = []
    levels: list[float] | None = None

    @field_validator("metric")
    @classmethod
    def _known_metric(cls, v):
        return validate_metric(v)

    @field_validator("temperature_K")
    @classmethod
    def _temperature_in_range(cls, v):
        return validate_temperature(v)

    @field_validator("thermo")
    @classmethod
    def _thermo_keys(cls, v):
        return validate_thermo_keys(v)

    @model_validator(mode="after")
    def _axes_usable(self):
        if self.x_axis.param == self.y_axis.param:
            raise ValueError("x_axis and y_axis must be different parameters")
        _check_coverage(self.fixed, (self.x_axis.param, self.y_axis.param))
        return self


class SliceRequest(BaseModel):
    """A 1D cut: one axis, everything else fixed. Same metrics as the 2D landscape."""

    axis: AxisSpec
    fixed: FixedParams = FixedParams()
    metric: str = "fold_change"
    temperature_K: float = T_REF_DEFAULT
    thermo: dict[str, ThermoSpecIn] = {}

    @field_validator("metric")
    @classmethod
    def _known_metric(cls, v):
        return validate_metric(v)

    @field_validator("temperature_K")
    @classmethod
    def _temperature_in_range(cls, v):
        return validate_temperature(v)

    @field_validator("thermo")
    @classmethod
    def _thermo_keys(cls, v):
        return validate_thermo_keys(v)

    @model_validator(mode="after")
    def _axis_usable(self):
        _check_coverage(self.fixed, (self.axis.param,))
        return self


class AxisOut(BaseModel):
    param: str
    scale: str
    unit: str
    values: list[float]


class ContourOut(BaseModel):
    level: float
    polylines: list[list[list[float]]]


class SlackOut(BaseModel):
    minus: float
    plus: float


class RobustOut(BaseModel):
    x: float
    y: float
    metric_value: float | None
    x_slack: SlackOut
    y_slack: SlackOut
    summary: str


class TmBandOut(BaseModel):
    axis: str
    min: float
    max: float
    Tm: float
    label: str


class LandscapeResponse(BaseModel):
    metric: str
    x: AxisOut
    y: AxisOut
    # Row-major, z[i][j] for x=values[i], y=values[j]. null where the model is undefined.
    z: list[list[float | None]]
    contours: list[ContourOut] = []
    feasible: list[list[bool]] | None = None
    robust_optimum: RobustOut | None = None
    tm_band: TmBandOut | None = None
    temperature_K: float | None = T_REF_DEFAULT
    thermo_mode: str = "isothermal"
    cells: int = 0
    warnings: list[str] = []
    units: dict[str, str] = {}


class SliceResponse(BaseModel):
    metric: str
    axis: AxisOut
    values: list[float | None]
    tm_band: TmBandOut | None = None
    temperature_K: float | None = T_REF_DEFAULT
    thermo_mode: str = "isothermal"
    warnings: list[str] = []
    units: dict[str, str] = {}
