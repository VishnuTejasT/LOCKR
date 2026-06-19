"""Shared data objects for the engine.

General LOCKR/lucCage containers — no field assumes a specific target, binder,
or cage variant. ECLIPSE-specific data lives in calibration.py and the test
files, not here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Base lucCage scaffold defaults (Quijano-Rubio et al. 2021); override per system.
RT_37C = 0.592
K_OPEN_DEFAULT = 1e-3
K_CK_DEFAULT = 1e-8
LUCKEY_DEFAULT = 500e-9
# TODO: confirm K_CK = lucKey-cage vs lucKey-SmBiT Kd, inconsistent in my source docs.


@dataclass(frozen=True)
class SensorParams:
    """The four constants that pin a LOCKR sensor's operating point."""

    K_open: float = K_OPEN_DEFAULT
    K_CK: float = K_CK_DEFAULT
    lucKey: float = LUCKEY_DEFAULT
    RT: float = RT_37C

    @property
    def luckey_ratio(self) -> float:
        # lucKey/K_CK dominance ratio — a diagnostic, not an achievable fold-change.
        return self.lucKey / self.K_CK


DEFAULT_PARAMS = SensorParams()


@dataclass
class TargetInterface:
    """Binder positions that contact the user's own target."""

    positions: list[int] = field(default_factory=list)
    label: str = ""


@dataclass
class BinderSequence:
    sequence: str
    name: str | None = None

    def __post_init__(self):
        self.sequence = self.sequence.strip().upper()

    def __len__(self) -> int:
        return len(self.sequence)

    def residues(self):
        # 1-indexed (position, residue) pairs.
        return list(enumerate(self.sequence, start=1))


@dataclass
class FoldChangeResult:
    target_conc: float
    Kd: float
    pull: float
    theta: float
    K_open_eff: float
    f_base: float
    f_signal: float
    fold_change: float


@dataclass
class ScanResult:
    label: str
    Kd: float
    pull: float
    max_fc: float
    ec50: float
    lod: float


@dataclass
class RegimeResult:
    luckey_dominance_ratio: float
    K_open: float
    regime: str
    max_fold_change: float
    latch_tuning_helps: bool
    verdict: str


@dataclass
class ChargeResult:
    net_charge: float
    pH: float
    helical_ok: bool
    helix_breakers: list[int] = field(default_factory=list)
    note: str = ""


@dataclass
class Liability:
    position: int
    residue: str
    weight: float
    penalty: float


@dataclass
class LiabilityReport:
    binder: BinderSequence
    liabilities: list[Liability] = field(default_factory=list)
    preserve_positions: list[int] = field(default_factory=list)
    net_charge: float = 0.0
    penalty_total: float = 0.0
    liability_score: float = 0.0
    liability_band: str = "low"
    K_CK_estimate: float = 0.0


@dataclass
class VariantSuggestion:
    policy: str
    sequence: str
    mutations: list[str] = field(default_factory=list)
    liability_score: float = 0.0
    liability_band: str = "low"
    K_CK_estimate: float = 0.0
