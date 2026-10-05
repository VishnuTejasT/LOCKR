"""Request/response shapes for POST /dosing, the Key concentration back-calculator.

Concentrations are nM and volumes are uL at this boundary, same nM convention
as foldchange.py/sweep.py; the route converts to/from dosing.py's M/L engine units.
"""

from __future__ import annotations

from pydantic import BaseModel, model_validator


class DosingRequest(BaseModel):
    k_ck: float
    total_cage_conc: float
    key_stock_conc: float
    assay_volume: float
    target_fraction_bound: float | None = None
    target_fold_change: float | None = None

    @model_validator(mode="after")
    def _positive_inputs(self):
        for field in ("k_ck", "total_cage_conc", "key_stock_conc", "assay_volume"):
            if getattr(self, field) <= 0:
                raise ValueError(f"{field} must be > 0")
        return self

    @model_validator(mode="after")
    def _exactly_one_target(self):
        # Fraction and fold-change describe the same target signal in two
        # units; giving both (or neither) leaves which one wins ambiguous.
        if (self.target_fraction_bound is None) == (self.target_fold_change is None):
            raise ValueError("Provide exactly one of target_fraction_bound or target_fold_change")
        if self.target_fraction_bound is not None and not (0 <= self.target_fraction_bound < 1):
            raise ValueError("target_fraction_bound must be in [0, 1)")
        if self.target_fold_change is not None and self.target_fold_change < 0:
            raise ValueError("target_fold_change must be >= 0")
        return self


class TitrationPoint(BaseModel):
    key_conc_nm: float
    fraction_bound: float


class DosingResponse(BaseModel):
    key_conc_nm: float
    fraction_bound: float
    volume_to_add_ul: float
    warnings: list[str] = []
    titration: list[TitrationPoint] = []
