"""Request/response shapes for POST /api/temperature/window.

Concentrations are nM at the API boundary and molar inside the engine, the same
convention foldchange, sweep and landscape use.

The spec for this module asks for 422 on a rejected request, but this app's global
handler turns every pydantic RequestValidationError into a 400 (42 existing tests
assert that). So the models below stay deliberately permissive about types and all
semantic checking happens in validate_request(), which raises ApiError with
status_code=422. Nothing here changes the status code of any other route.
"""

from __future__ import annotations

from pydantic import BaseModel

from lockr.engine.temperature_window import (DEFAULT_CONFIDENCE, DEFAULT_DG_BINDER_SD,
                                             DEFAULT_FC_THRESHOLD, DEFAULT_N_SAMPLES,
                                             DEFAULT_PULL, MAX_GRID_POINTS, MAX_N_SAMPLES,
                                             MIN_N_SAMPLES, PRIOR_KEYS)

from ..errors import ApiError
from .common import validate_sequence

NM_TO_M = 1e-9

UNITS = {"dG": "kcal/mol", "dH": "kcal/mol", "temperature": "C",
         "key_conc": "nM", "target_conc": "nM", "kd_nM": "nM",
         "fold_change": "x (ratio)", "p_works": "fraction of samples"}


def _reject(message: str, field: str, code: str = "VALIDATION_ERROR") -> None:
    raise ApiError(code, message, field, status_code=422)


# EVERY NUMBER BELOW IS A PLACEHOLDER. This is an example request shape, not a result,
# and nothing here is measured ECLIPSE data. Defined once here and served to the UI by
# GET /api/temperature/window/preset so the browser cannot drift from the test fixture.
ECLIPSE_PRESET = {
    "systems": [{
        "name": "ECLIPSE cage 1",
        # TODO(placeholder): illustrative latch-cage stability, not measured.
        "dG_latch_cage": -2.7,
        "dH_latch_cage": None,
        # TODO(placeholder): -10.9 kcal/mol is R*T*ln(K_CK) at K_CK = 1e-8, the lucCage
        # scaffold default in models.K_CK_DEFAULT. Replace with your own Key-Cage number.
        "dG_key_cage": -10.9,
        "dH_key_cage": None,
        "key_conc": None,
        "target_conc": None,
        "t_max_stability_c": None,
        "binder": {
            "name": "ECLIPSE binder 1",
            # provenance/display only; never used to predict affinity
            "sequence": "LISDAELEAIFAEELDC",
            # TODO(placeholder): replace with the HybridDock-Pep output for this binder.
            "dG_target_binder": -10.0,
            "kd_nM": None,
            # docking-derived affinities are noisy; narrow this only for an ITC measurement
            "dG_target_binder_sd": 1.5,
            "dH_target_binder": None,
            "tm_c": None,
            "graft_dG_latch_cage_shift": 0.0,
        },
    }],
    "t_min_c": 4.0,
    "t_max_c": 45.0,
    "t_step_c": 1.0,
    "t_ref_c": 25.0,
    "fc_threshold": DEFAULT_FC_THRESHOLD,
    "confidence": DEFAULT_CONFIDENCE,
    "n_samples": DEFAULT_N_SAMPLES,
}

ECLIPSE_PRESET_PLACEHOLDER_FIELDS = (
    "systems[0].dG_latch_cage", "systems[0].dG_key_cage",
    "systems[0].binder.dG_target_binder", "systems[0].binder.dG_target_binder_sd",
)


class DhPriorIn(BaseModel):
    mean: float | None = None
    sd: float | None = None


class BinderIn(BaseModel):
    name: str | None = None
    # provenance and display only. Nothing here predicts affinity from sequence.
    sequence: str | None = None
    dG_target_binder: float | None = None
    kd_nM: float | None = None
    dG_target_binder_sd: float | None = None
    dH_target_binder: float | None = None
    tm_c: float | None = None
    graft_dG_latch_cage_shift: float | None = None


class SystemIn(BaseModel):
    name: str | None = None
    dG_latch_cage: float | None = None
    dH_latch_cage: float | None = None
    dG_key_cage: float | None = None
    dH_key_cage: float | None = None
    key_conc: float | None = None       # nM
    target_conc: float | None = None    # nM
    t_max_stability_c: float | None = None
    pull: float | None = None
    binder: BinderIn | None = None


class TemperatureWindowRequest(BaseModel):
    systems: list[SystemIn] = []
    t_min_c: float = 4.0
    t_max_c: float = 45.0
    t_step_c: float = 1.0
    t_ref_c: float = 25.0
    fc_threshold: float = DEFAULT_FC_THRESHOLD
    confidence: float = DEFAULT_CONFIDENCE
    n_samples: int = DEFAULT_N_SAMPLES
    dH_priors: dict[str, DhPriorIn | None] | None = None
    seed: int | None = None
    luciferase_activity: list[list[float]] | None = None


def _validate_binder(binder: BinderIn | None, where: str) -> None:
    if binder is None:
        _reject("binder is required: a construct is a base LOCKR switch plus a grafted "
                "binder, and without one there is no Binder-Target interaction to scan.",
                f"{where}.binder")
    has_dG = binder.dG_target_binder is not None
    has_kd = binder.kd_nM is not None
    if has_dG and has_kd:
        _reject("give exactly one of dG_target_binder or kd_nM, not both; they are the same "
                "number in different units and two values can disagree.",
                f"{where}.binder.dG_target_binder")
    if not (has_dG or has_kd):
        _reject("binder needs an affinity: supply either dG_target_binder (kcal/mol) or "
                "kd_nM (nM).", f"{where}.binder.dG_target_binder")
    if has_kd and binder.kd_nM <= 0:
        _reject("kd_nM must be > 0.", f"{where}.binder.kd_nM")
    if binder.dG_target_binder_sd is not None and binder.dG_target_binder_sd < 0:
        _reject("dG_target_binder_sd must be >= 0.", f"{where}.binder.dG_target_binder_sd")
    if binder.tm_c is not None and not 0.0 <= binder.tm_c <= 100.0:
        _reject("tm_c must be between 0 and 100 C.", f"{where}.binder.tm_c")
    if binder.sequence is not None:
        try:
            validate_sequence(binder.sequence)
        except ValueError as exc:
            _reject(f"{exc} Use the 20 standard amino acid letters; the sequence is for "
                    f"provenance only and is never used to predict affinity.",
                    f"{where}.binder.sequence")


def validate_request(request: TemperatureWindowRequest) -> None:
    """Every rejection this module promises, as a 422 with a message a user can act on."""
    if not 1 <= len(request.systems) <= 10:
        _reject("systems must hold between 1 and 10 constructs.", "systems")
    if request.t_min_c >= request.t_max_c:
        _reject("t_min_c must be < t_max_c.", "t_min_c")
    if request.t_step_c <= 0:
        _reject("t_step_c must be > 0.", "t_step_c")
    points = int((request.t_max_c - request.t_min_c) / request.t_step_c) + 1
    if points > MAX_GRID_POINTS:
        _reject(f"sweep has {points} grid points, above the {MAX_GRID_POINTS} point cap; "
                f"widen t_step_c or narrow the range.", "t_step_c")
    if not MIN_N_SAMPLES <= request.n_samples <= MAX_N_SAMPLES:
        _reject(f"n_samples must be between {MIN_N_SAMPLES} and {MAX_N_SAMPLES}.", "n_samples")
    if not 0.0 < request.confidence < 1.0:
        _reject("confidence must be between 0 and 1, exclusive.", "confidence")
    if request.fc_threshold <= 0:
        _reject("fc_threshold must be > 0.", "fc_threshold")
    if not 200.0 < request.t_ref_c + 273.15 < 500.0:
        _reject("t_ref_c must be between -73 and 227 C.", "t_ref_c")

    for key in (request.dH_priors or {}):
        if key not in PRIOR_KEYS:
            _reject(f"unknown dH prior {key!r}; use {list(PRIOR_KEYS)}.", f"dH_priors.{key}")
    for key, prior in (request.dH_priors or {}).items():
        if prior is not None and prior.sd is not None and prior.sd < 0:
            _reject("dH prior sd must be >= 0.", f"dH_priors.{key}.sd")

    for point in (request.luciferase_activity or []):
        if len(point) != 2:
            _reject("each luciferase_activity entry must be a [temperature_C, "
                    "relative_activity] pair.", "luciferase_activity")
        if point[1] < 0:
            _reject("luciferase relative_activity must be >= 0.", "luciferase_activity")

    for i, system in enumerate(request.systems):
        where = f"systems.{i}"
        if system.dG_latch_cage is None:
            _reject("dG_latch_cage is required (kcal/mol at t_ref_c).", f"{where}.dG_latch_cage")
        if system.dG_key_cage is None:
            _reject("dG_key_cage is required (kcal/mol at t_ref_c).", f"{where}.dG_key_cage")
        for field in ("key_conc", "target_conc"):
            value = getattr(system, field)
            if value is not None and value <= 0:
                _reject(f"{field} must be > 0 nM.", f"{where}.{field}")
        if system.pull is not None and system.pull < 0:
            _reject("pull must be >= 0.", f"{where}.pull")
        if system.t_max_stability_c is not None and not 0.0 <= system.t_max_stability_c <= 100.0:
            _reject("t_max_stability_c must be between 0 and 100 C.",
                    f"{where}.t_max_stability_c")
        _validate_binder(system.binder, where)


class WindowSpanOut(BaseModel):
    t_low_c: float
    t_high_c: float


class SystemWindowOut(BaseModel):
    name: str
    binder_name: str
    t_opt_c: float | None
    peak_median_fc: float | None
    median_fc: list[float | None]
    p10_fc: list[float | None]
    p90_fc: list[float | None]
    p_works: list[float]
    window: WindowSpanOut | None
    other_segments: list[WindowSpanOut] = []
    window_touches_grid_edge: bool = False
    reason: str | None = None
    dG_latch_cage_raw: float
    dG_latch_cage_effective: float
    graft_dG_latch_cage_shift: float
    dG_target_binder: float
    dG_target_binder_sd: float
    stability_cap_c: float | None = None
    stability_cap_source: str | None = None
    luciferase_factor: list[float] | None = None
    binder_affects_result: bool = True


class EchoOut(BaseModel):
    """What the engine actually ran, after defaults and unit conversion."""

    t_min_c: float
    t_max_c: float
    t_step_c: float
    t_ref_c: float
    fc_threshold: float
    confidence: float
    n_samples: int
    dH_priors: dict[str, dict[str, float]]
    luciferase_activity: list[list[float]] | None = None
    default_pull: float = DEFAULT_PULL
    default_dG_target_binder_sd: float = DEFAULT_DG_BINDER_SD


class TemperatureWindowResponse(BaseModel):
    temperatures_c: list[float]
    systems: list[SystemWindowOut]
    seed: int
    echo: EchoOut
    notes: list[str] = []
    units: dict[str, str] = {}
