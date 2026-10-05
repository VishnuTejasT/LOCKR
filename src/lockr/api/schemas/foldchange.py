"""Request/response shapes for POST /foldchange, per spec 6.2."""

from __future__ import annotations

from pydantic import BaseModel, field_validator, model_validator

from lockr.engine.temperature import T_REF_DEFAULT

from .common import ThermoSpecIn, validate_temperature, validate_thermo_keys


class FoldChangeRequest(BaseModel):
    k_ck: float
    k_open: float
    pull: float
    luckey: float
    k_target: float | None = None
    target_conc: float | None = None
    temperature_K: float = T_REF_DEFAULT
    thermo: dict[str, ThermoSpecIn] = {}

    @field_validator("temperature_K")
    @classmethod
    def _temperature_in_range(cls, v):
        return validate_temperature(v)

    @field_validator("thermo")
    @classmethod
    def _thermo_keys(cls, v):
        return validate_thermo_keys(v)

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
        # Affinity and target amount describe one real sample together; one without the other isn't computable.
        if (self.k_target is None) != (self.target_conc is None):
            raise ValueError("K_target (binder-target affinity) and target concentration must be "
                             "provided together, or both left blank to assume saturating target")
        return self


class FoldChangeResponse(BaseModel):
    fold_change: float
    temperature_K: float = T_REF_DEFAULT
    thermo_mode: str = "isothermal"
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
