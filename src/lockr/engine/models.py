"""Shared data objects for the engine.

Plain containers passed between thermo / charge / liability so the modules don't
reach into each other's dict-shaped returns. No physics here — that lives in the
sibling modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# ECLIPSE assay constants, centralised so one SensorParams flows through every
# module and the wet-lab team can override a single value (e.g. titrate lucKey)
# without touching the equations. K_open and K_CK come from the lucCage framework
# (Langan 2019; Quijano-Rubio 2021); lucKey concentration is our own assay design
# (chosen at 50x K_CK).
RT_37C = 0.592       # kcal/mol, R * 310.15 K
K_OPEN_DEFAULT = 1e-3
# K_CK: dissociation constant for lucKey binding the SmBiT epitope exposed on the
# open cage. (The two PDFs phrase this as "lucKey-SmBiT" vs "lucKey-Cage" — same
# event, different resolution; lucKey grabs SmBiT, which is only accessible once
# the latch opens.)
K_CK_DEFAULT = 1e-8  # M
LUCKEY_DEFAULT = 500e-9  # M


@dataclass(frozen=True)
class SensorParams:
    """The four constants that pin a LOCKR sensor's operating point.

    frozen on purpose: a params object is an experimental condition, and mutating
    one mid-scan is exactly the bug that would make a dud variant look like it
    worked.
    """

    K_open: float = K_OPEN_DEFAULT
    K_CK: float = K_CK_DEFAULT
    lucKey: float = LUCKEY_DEFAULT
    RT: float = RT_37C

    @property
    def luckey_ratio(self) -> float:
        # [lucKey]/K_CK — the lucKey *dominance ratio*. It's the term that swamps
        # the fold-change denominator at 500 nM lucKey (=50) and decides the
        # regime. NB: this is NOT the achievable fold-change ceiling — at usable
        # pull (10-30) the realised max FC is ~1+pull (11-31x). Don't show 50 as a
        # max-FC target. See thermo.max_fold_change.
        return self.lucKey / self.K_CK


DEFAULT_PARAMS = SensorParams()


@dataclass
class BinderSequence:
    # 1-indexed positions throughout to match the docs/lab numbering (the v1
    # binder's D/E are "positions 4, 6, 8, ...").
    sequence: str
    name: str | None = None

    def __post_init__(self):
        self.sequence = self.sequence.strip().upper()

    def __len__(self) -> int:
        return len(self.sequence)

    def residues(self):
        return list(enumerate(self.sequence, start=1))


@dataclass
class FoldChangeResult:
    # Everything needed to reconstruct one FC point by hand if a reviewer asks.
    pfldh: float
    Kd: float
    pull: float
    theta: float
    K_open_eff: float
    f_base: float
    f_signal: float
    fold_change: float


@dataclass
class ScanResult:
    """Summary of a full [PfLDH] titration for one (Kd, pull) pair."""

    label: str
    Kd: float
    pull: float
    max_fc: float
    ec50: float       # M
    lod: float        # M, taken as 0.1 * EC50 (10% threshold) per the docs


@dataclass
class RegimeResult:
    """Section 8 diagnostic: key-limited vs K_open-limited.

    Two distinct numbers, deliberately not conflated:
    - luckey_dominance_ratio = [lucKey]/K_CK, the diagnostic that says whether
      K_open even matters. NOT a fold-change you can hit.
    - max_fold_change = the realised FC at saturating target (~1+pull), i.e. what
      you actually get at this pull.
    """

    luckey_dominance_ratio: float
    K_open: float
    regime: str               # "key-limited" | "K_open-limited" | "mixed"
    max_fold_change: float    # realised FC at saturating target
    latch_tuning_helps: bool
    verdict: str              # UI-ready, plain language


@dataclass
class ChargeResult:
    net_charge: float
    pH: float
    helical_ok: bool
    helix_breakers: list[int] = field(default_factory=list)  # 1-indexed P/G
    note: str = ""


@dataclass
class Liability:
    """A single flagged acidic residue at the lucKey-cage interface."""

    position: int     # 1-indexed
    residue: str      # 'D' or 'E'
    weight: float     # structural weight feeding the score
    penalty: float    # kcal/mol electrostatic cost attributed to this residue


@dataclass
class LiabilityReport:
    binder: BinderSequence
    liabilities: list[Liability] = field(default_factory=list)
    preserved_interface: list[int] = field(default_factory=list)
    net_charge: float = 0.0
    penalty_total: float = 0.0     # kcal/mol
    liability_score: float = 0.0   # 0–100, position-weighted
    liability_band: str = "low"    # low | moderate | high
    K_CK_grafted: float = 0.0      # M, after applying the penalty


@dataclass
class VariantSuggestion:
    policy: str            # "conservative" | "neutralizing"
    sequence: str
    mutations: list[str] = field(default_factory=list)  # e.g. ["D4A", "E6A"]
    liability_score: float = 0.0
    liability_band: str = "low"
    K_CK_grafted: float = 0.0
