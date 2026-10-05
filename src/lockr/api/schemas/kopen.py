"""Request/response shapes for POST /kopen-scenarios, no spec entry, built to
back the Assembly tab's K_open Tuning Analysis section (see lockr-tool-plan.md
§7: the engine is the one source of truth for the math, so this what-if grid
is computed server-side via thermo.py rather than duplicated in JS)."""

from __future__ import annotations

from pydantic import BaseModel, field_validator

from lockr.engine.temperature import T_REF_DEFAULT

from .common import ThermoSpecIn, validate_temperature, validate_thermo_keys


class KopenScenariosRequest(BaseModel):
    k_open_current: float
    k_ck: float
    luckey: float
    pull: float
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

    @field_validator("k_open_current", "k_ck", "luckey")
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


class KopenScenario(BaseModel):
    mutations: int
    preset: str
    destab_kcal_per_mut: float | None
    k_open: float
    fold_change: float
    regime: str
    helps: str


class KopenScenariosResponse(BaseModel):
    dominance_ratio: float
    crossover_luckey_nm: float
    baseline: KopenScenario
    scenarios: list[KopenScenario]
    temperature_K: float = T_REF_DEFAULT
    thermo_mode: str = "isothermal"
    warnings: list[str] = []
