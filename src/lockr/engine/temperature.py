"""Temperature-dependent free energies for the LucCage fold-change model.

thermo.py treats every equilibrium constant as fixed, which is only true at the
temperature it was measured at. This layer re-derives dG at a requested
temperature from Gibbs-Helmholtz with a heat-capacity term, and reports which of
three input modes it used, because the three do not deserve equal trust:

    full        dH_ref, dS_ref and dCp known -> exact over a wide range.
    vant_hoff   dH_ref known, dCp taken as 0 -> only near T_ref.
    isothermal  only dG_ref known -> NO extrapolation happens at all.

Isothermal mode returns the reference value unchanged and says so loudly. It
never pretends to have scaled anything.

Sign convention throughout: for any equilibrium constant K used by thermo.py
(K_open, K_CK, K_target), dG = R*T*ln(K). So dG_latch_cage = R*T*ln(K_open) is
negative for a cage that prefers to stay closed, and dG_key_cage = R*T*ln(K_CK)
is the Key-Cage binding energy. dH_ref must follow the same direction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

R_KCAL = 1.9872041e-3  # kcal/mol/K
T_REF_DEFAULT = 298.15

MODE_FULL = "full"
MODE_VANT_HOFF = "vant_hoff"
MODE_ISOTHERMAL = "isothermal"

# Weakest mode wins when a response mixes several parameters.
_MODE_RANK = {MODE_FULL: 0, MODE_VANT_HOFF: 1, MODE_ISOTHERMAL: 2}

_VANT_HOFF_WINDOW_K = 15.0
_TM_MARGIN_K = 10.0

# Murphy-Freire/Spolar-Record surface-area coefficients, cal/mol/K per A^2.
_DCP_PER_NONPOLAR_A2 = 0.32
_DCP_PER_POLAR_A2 = -0.14

DCP_ESTIMATE_NOTE = (
    "dCp was estimated from buried surface area (Murphy-Freire coefficients), not measured. "
    "Treat the resulting temperature extrapolation as indicative only."
)


def isothermal_warning(label: str, T: float, T_ref: float) -> str:
    return (f"TEMPERATURE SCALING UNAVAILABLE for {label}: only a reference dG at {T_ref:.2f} K was "
            f"supplied, so the value at {T:.2f} K is the UNSCALED {T_ref:.2f} K value. Supply dH_ref "
            f"(and ideally dCp) to extrapolate. Do not read the fold change as a prediction at "
            f"{T:.2f} K.")


@dataclass(frozen=True)
class ThermoSpec:
    """Temperature dependence of one free energy. dG_ref may be left None and
    filled in from the equilibrium constant the caller already has."""

    dG_ref: float | None = None
    dH_ref: float | None = None
    dS_ref: float | None = None
    dCp: float = 0.0
    T_ref: float = T_REF_DEFAULT
    Tm: float | None = None
    label: str = "dG"

    def __post_init__(self) -> None:
        if self.T_ref <= 0:
            raise ValueError("T_ref must be > 0 K")
        if self.Tm is not None and self.Tm <= 0:
            raise ValueError("Tm must be > 0 K")

    @property
    def mode(self) -> str:
        if self.dH_ref is None:
            return MODE_ISOTHERMAL
        return MODE_FULL if self.dCp else MODE_VANT_HOFF

    @property
    def dS_ref_resolved(self) -> float | None:
        # dH and dG pin dS at T_ref, so a supplied dS is only needed when there is no dG anchor.
        if self.dH_ref is None:
            return None
        if self.dG_ref is None:
            return self.dS_ref
        return (self.dH_ref - self.dG_ref) / self.T_ref

    @property
    def dG_ref_resolved(self) -> float | None:
        if self.dG_ref is not None:
            return self.dG_ref
        if self.dH_ref is None or self.dS_ref is None:
            return None
        return self.dH_ref - self.T_ref * self.dS_ref


def anchored(spec: ThermoSpec, dG_ref: float) -> ThermoSpec:
    """Pin a spec to a dG the caller measured, keeping its dH/dCp shape."""
    return spec if spec.dG_ref is not None else replace(spec, dG_ref=dG_ref)


def dG_at(spec: ThermoSpec, T: float) -> float:
    """Gibbs-Helmholtz with heat capacity. Isothermal specs return dG_ref as-is."""
    if T <= 0:
        raise ValueError("temperature must be > 0 K")
    dG_ref = spec.dG_ref_resolved
    if dG_ref is None:
        raise ValueError(f"{spec.label}: need dG_ref, or both dH_ref and dS_ref")
    if spec.mode is MODE_ISOTHERMAL or T == spec.T_ref:
        return dG_ref
    dS_ref = spec.dS_ref_resolved
    return (spec.dH_ref - T * dS_ref
            + spec.dCp * ((T - spec.T_ref) - T * math.log(T / spec.T_ref)))


def warnings_at(spec: ThermoSpec, T: float) -> list[str]:
    out: list[str] = []
    if T != spec.T_ref and spec.mode is MODE_ISOTHERMAL:
        out.append(isothermal_warning(spec.label, T, spec.T_ref))
    if spec.mode is MODE_VANT_HOFF and abs(T - spec.T_ref) > _VANT_HOFF_WINDOW_K:
        out.append(f"{spec.label}: van't Hoff mode (dCp assumed 0) extrapolated {abs(T - spec.T_ref):.1f} K "
                   f"from T_ref, beyond the ~{_VANT_HOFF_WINDOW_K:.0f} K window where that holds. "
                   f"Supply dCp for a wider range.")
    if spec.Tm is not None and T >= spec.Tm - _TM_MARGIN_K:
        out.append(f"{spec.label}: {T:.2f} K is within {_TM_MARGIN_K:.0f} K of Tm = {spec.Tm:.2f} K. The "
                   f"two-state open/closed switch model breaks down here because the cage itself is "
                   f"unfolding; fold change from this calculation is not meaningful.")
    return out


def estimate_dCp(buried_nonpolar_A2: float, buried_polar_A2: float = 0.0) -> float:
    """Rough dCp (kcal/mol/K) from buried surface area. An estimate, not a measurement."""
    if buried_nonpolar_A2 < 0 or buried_polar_A2 < 0:
        raise ValueError("buried surface areas must be >= 0")
    cal = _DCP_PER_NONPOLAR_A2 * buried_nonpolar_A2 + _DCP_PER_POLAR_A2 * buried_polar_A2
    return cal / 1000.0


def k_from_dG(dG: float, T: float) -> float:
    return math.exp(dG / (R_KCAL * T))


def dG_from_k(k: float, T: float) -> float:
    if k <= 0:
        raise ValueError("equilibrium constant must be > 0")
    return R_KCAL * T * math.log(k)


def weakest_mode(modes) -> str:
    modes = list(modes)
    return max(modes, key=lambda m: _MODE_RANK[m]) if modes else MODE_ISOTHERMAL


def shift_constant(k_ref: float | None, T: float,
                   spec: ThermoSpec | None = None) -> tuple[float | None, list[str]]:
    """Re-derive an equilibrium constant measured at spec.T_ref for temperature T.

    Returns k_ref untouched (bit-identical) when T is T_ref or when the spec is
    isothermal, so existing Module B results are unchanged by default.
    """
    spec = spec or ThermoSpec()
    if k_ref is None:
        return None, []
    warns = warnings_at(spec, T)
    if T == spec.T_ref or spec.mode is MODE_ISOTHERMAL:
        return k_ref, warns
    return k_from_dG(dG_at(anchored(spec, dG_from_k(k_ref, spec.T_ref)), T), T), warns
