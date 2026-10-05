"""Shapes shared by more than one route's schema."""

from __future__ import annotations

from pydantic import BaseModel, field_validator

from lockr.engine.temperature import (DCP_ESTIMATE_NOTE, T_REF_DEFAULT, ThermoSpec, estimate_dCp)

_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Keys a thermo block may carry: one per free energy the fold-change model uses.
THERMO_KEYS = ("dG_latch_cage", "dG_key_cage", "dG_target_binder")


def validate_sequence(sequence: str) -> str:
    seq = sequence.strip().upper()
    if not seq:
        raise ValueError("sequence is empty")
    for i, aa in enumerate(seq, 1):
        if aa not in _STANDARD_AA:
            raise ValueError(f"Contains non-standard residue {aa!r} at position {i}.")
    return seq


class Window(BaseModel):
    start: int
    end: int

    @field_validator("end")
    @classmethod
    def _end_after_start(cls, end, info):
        start = info.data.get("start")
        if start is not None and end < start:
            raise ValueError("window end must be >= start")
        return end

    def clamped(self, length: int) -> tuple[int, int]:
        # Per spec 8: "Sensitive window outside length -> Clamp handles to 1-N."
        return max(1, min(self.start, length)), max(1, min(self.end, length))


class ThermoSpecIn(BaseModel):
    """Temperature dependence of one free energy, as a request carries it.

    dH_ref present -> van't Hoff; dH_ref and dCp present -> full Gibbs-Helmholtz;
    neither -> isothermal, and no extrapolation is performed at all.
    """

    dH_ref: float | None = None
    dS_ref: float | None = None
    dCp: float | None = None
    T_ref: float = T_REF_DEFAULT
    Tm: float | None = None
    buried_nonpolar_A2: float | None = None
    buried_polar_A2: float = 0.0

    @field_validator("T_ref", "Tm")
    @classmethod
    def _positive_kelvin(cls, v, info):
        if v is not None and v <= 0:
            raise ValueError(f"{info.field_name} must be > 0 K")
        return v

    @field_validator("buried_nonpolar_A2", "buried_polar_A2")
    @classmethod
    def _non_negative_area(cls, v, info):
        if v is not None and v < 0:
            raise ValueError(f"{info.field_name} must be >= 0")
        return v

    def to_spec(self, label: str) -> tuple[ThermoSpec, list[str]]:
        notes: list[str] = []
        dCp = self.dCp
        if dCp is None and self.buried_nonpolar_A2 is not None:
            dCp = estimate_dCp(self.buried_nonpolar_A2, self.buried_polar_A2)
            notes.append(f"{label}: dCp = {dCp:.4g} kcal/mol/K. {DCP_ESTIMATE_NOTE}")
        return ThermoSpec(dH_ref=self.dH_ref, dS_ref=self.dS_ref, dCp=dCp or 0.0,
                          T_ref=self.T_ref, Tm=self.Tm, label=label), notes


def validate_thermo_keys(thermo: dict) -> dict:
    unknown = [k for k in thermo if k not in THERMO_KEYS]
    if unknown:
        raise ValueError(f"unknown thermo key(s) {unknown}; use {list(THERMO_KEYS)}")
    return thermo


def build_specs(thermo: dict | None) -> tuple[dict[str, ThermoSpec], list[str]]:
    """Request thermo block -> engine specs, plus any notes the conversion produced."""
    specs, notes = {}, []
    for key in THERMO_KEYS:
        entry = (thermo or {}).get(key)
        if entry is None:
            specs[key] = ThermoSpec(label=key)
            continue
        specs[key], entry_notes = entry.to_spec(key)
        notes += entry_notes
    return specs, notes


def validate_temperature(T: float) -> float:
    if not 200.0 < T < 500.0:
        raise ValueError("temperature_K must be between 200 and 500 K")
    return T
