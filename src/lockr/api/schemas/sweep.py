"""Request/response shapes for POST /sweep, per spec 6.3."""

from __future__ import annotations

from pydantic import BaseModel, field_validator, model_validator

_SWEEPABLE = ("k_ck", "k_open_off", "k_open_on", "luckey")
_SCALES = ("log", "linear")


class BaseParams(BaseModel):
    k_ck: float
    k_open_off: float
    k_open_on: float
    luckey: float

    @field_validator("k_ck", "k_open_off", "k_open_on", "luckey")
    @classmethod
    def _positive(cls, v, info):
        if v <= 0:
            raise ValueError(f"{info.field_name} must be > 0")
        return v


class SweepSpec(BaseModel):
    param: str
    min: float
    max: float
    steps: int
    scale: str = "log"

    @field_validator("param")
    @classmethod
    def _known_param(cls, v):
        if v not in _SWEEPABLE:
            raise ValueError(f"param must be one of {_SWEEPABLE}")
        return v

    @field_validator("scale")
    @classmethod
    def _known_scale(cls, v):
        if v not in _SCALES:
            raise ValueError(f"scale must be one of {_SCALES}")
        return v

    @field_validator("steps")
    @classmethod
    def _enough_steps(cls, v):
        if v < 2:
            raise ValueError("steps must be >= 2")
        return v

    @model_validator(mode="after")
    def _min_below_max(self):
        if self.min >= self.max:
            raise ValueError("sweep min must be < max")
        return self


class SweepRequest(BaseModel):
    base_params: BaseParams
    sweep: SweepSpec


class SweepPoint(BaseModel):
    x: float
    fold_change: float
    dominance_ratio: float


class OperatingPoint(BaseModel):
    x: float
    fold_change: float


class SweepResponse(BaseModel):
    param: str
    scale: str
    points: list[SweepPoint]
    operating_point: OperatingPoint
