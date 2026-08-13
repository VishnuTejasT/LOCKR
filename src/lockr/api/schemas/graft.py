"""Request/response shapes for POST /graft and GET /graft/status."""

from __future__ import annotations

from pydantic import BaseModel, field_validator, model_validator

from .common import validate_sequence
from ...engine.graft import LATCH_END, LATCH_START

_MAX_BINDER_LENGTH = 35


class GraftStatusResponse(BaseModel):
    available: bool
    version: str | None
    template_bundled: bool
    calibration_warning: str | None = None


class GraftRequest(BaseModel):
    sequence: str
    scan_all: bool = True
    specific_position: int | None = None
    latch_start: int = LATCH_START
    latch_end: int = LATCH_END

    @field_validator("sequence")
    @classmethod
    def _valid_sequence(cls, v):
        v = validate_sequence(v)
        if len(v) > _MAX_BINDER_LENGTH:
            raise ValueError(
                f"Binder sequence ({len(v)} aa) is too long for the "
                f"{_MAX_BINDER_LENGTH}-residue latch window. Maximum binder length is "
                f"{_MAX_BINDER_LENGTH}aa for a single graft, or 17aa for a tandem design "
                "leaving room for a linker."
            )
        return v

    @model_validator(mode="after")
    def _position_in_range(self):
        if self.specific_position is None:
            return self
        max_start = self.latch_end - len(self.sequence) + 1
        if not (self.latch_start <= self.specific_position <= max_start):
            raise ValueError(
                f"specific_position must be between {self.latch_start} and {max_start} "
                f"for a {len(self.sequence)}aa binder in this latch window."
            )
        return self


class ScoredPosition(BaseModel):
    position: int
    score: float


class GraftResponse(BaseModel):
    best_position: int
    best_score: float
    verdict: str
    all_scores: list[ScoredPosition]
    grafted_sequence: str
    job_id: str
    runtime_seconds: float
    binder_length: int
    calibration_warning: str | None = None
