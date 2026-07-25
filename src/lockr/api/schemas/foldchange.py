"""Request/response shapes for POST /foldchange, per spec 6.2."""

from __future__ import annotations

from pydantic import BaseModel, field_validator, model_validator


class FoldChangeRequest(BaseModel):
    k_ck: float
    k_open: float
    pull: float
    luckey: float
    k_target: float | None = None
    target_conc: float | None = None

    @field_validator("k_ck", "k_open", "luckey")
    @classmethod
    def _positive(cls, v, info):
        if v <= 0:
            raise ValueError(f"{info.field_name} must be > 0")
        return v

    @field_validator("pull")
    @classmethod
    def _pull_non_negative(cls, v):
        if v < 0:
            raise ValueError("pull must be >= 0")
        return v

    @field_validator("k_target", "target_conc")
    @classmethod
    def _positive_if_set(cls, v, info):
        if v is not None and v <= 0:
            raise ValueError(f"{info.field_name} must be > 0 if provided")
        return v

    @model_validator(mode="after")
    def _target_fields_paired(self):
        # These describe one real sample together (affinity + how much target
        # is in it), one without the other isn't a computable state.
        if (self.k_target is None) != (self.target_conc is None):
            raise ValueError("K_target (binder-target affinity) and target concentration must be "
                             "provided together, or both left blank to assume saturating target")
        return self


class FoldChangeResponse(BaseModel):
    fold_change: float
    max_fold_change: float
    quality: str
    interpretation: str
    dominance_ratio: float
    fraction_of_dominance_ratio: float
    regime: str
    limiting_factor: str
    verdict: str
    recommendations: list[str]
    warnings: list[str] = []
    lod_2x_nm: float | None = None
    lod_3x_nm: float | None = None
    ec50_nm: float | None = None
