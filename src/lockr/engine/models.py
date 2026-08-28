"""This has the shared dataclasses for the engine, general to any LucCage LOCKR system."""

from __future__ import annotations

from dataclasses import dataclass, field

# Base lucCage scaffold defaults (Quijano-Rubio et al. 2021) override per system.
RT_37C = 0.592
K_OPEN_DEFAULT = 1e-3
K_CK_DEFAULT = 1e-8
LUCKEY_DEFAULT = 500e-9



@dataclass(frozen=True)
class SensorParams:
    K_open: float = K_OPEN_DEFAULT
    K_CK: float = K_CK_DEFAULT
    lucKey: float = LUCKEY_DEFAULT
    RT: float = RT_37C

    @property
    def luckey_ratio(self) -> float:
        # lucKey/K_CK dominance ratio
        return self.lucKey / self.K_CK


DEFAULT_PARAMS = SensorParams()


@dataclass
class TargetInterface:
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
class LodResult:
    """Detection-limit read from a target-concentration provided by the user.
    """

    lod_2x: float | None
    lod_3x: float | None
    ec50: float | None


@dataclass
class RegimeResult:
    luckey_dominance_ratio: float
    K_open: float
    regime: str
    max_fold_change: float
    latch_tuning_helps: bool
    verdict: str
    recommendations: list[str] = field(default_factory=list)


@dataclass
class ChargeResult:
    net_charge: float
    pH: float
    helical_ok: bool
    helix_breakers: list[int] = field(default_factory=list)
    note: str = ""


@dataclass
class HelixIssue:
    # severity is "blocking", and "warning" (costs stability) or "info".
    position: int | None
    severity: str
    kind: str
    message: str


@dataclass
class CyclizationEvidence:
    "The sequence evidence for potential cyclization, but doesn't guarantee it."
    possibly_cyclic: bool
    cysteine_positions: list[int] = field(default_factory=list)
    signals: list[str] = field(default_factory=list)


@dataclass
class HelixReport:
    sequence: str
    helix_confidence: float
    band: str
    mean_propensity: float
    hydrophobic_moment: float
    salt_bridges: list[tuple[int, int]] = field(default_factory=list)
    issues: list[HelixIssue] = field(default_factory=list)
    cyclization: CyclizationEvidence | None = None
    graft_blocked: bool = False


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
    liability_band: str = "Low"
    K_CK_estimate: float = 0.0


@dataclass
class VariantSuggestion:
    policy: str
    sequence: str
    mutations: list[str] = field(default_factory=list)
    liability_score: float = 0.0
    liability_band: str = "Low"
    K_CK_estimate: float = 0.0




@dataclass
class GraftResult:
    """Threading-scan result includes the best position and score."""

    best_position: int
    best_score: float
    verdict: str
    all_scores: list[tuple[int, float]]
    grafted_sequence: str
    grafted_pdb_path: str
    binder_length: int
    runtime_seconds: float
    calibration_warning: str | None = None


@dataclass
class GraftAtResult:
    """Includes the score and verdict for that specific position."""

    position: int
    score: float
    verdict: str
    grafted_sequence: str
    grafted_pdb_path: str
    calibration_warning: str | None = None


@dataclass
class ProtectedRegion:
    """Represents a protected region in the sequence that must not be mutated."""
    motif: str
    start: int
    end: int
    label: str = ""


@dataclass
class LatchWindow:
    """Represents a window in the sequence where a latch is expected."""
    start: int
    end: int
    expected_length: int | None = None


@dataclass
class GraftSpec:
    """Specifies the parameters for grafting a binder into a sequence."""
    binder: str
    start: int
    spacer: str | None = None
    spacer_start: int | None = None
    linker: str | None = None
    linker_start: int | None = None
    binder2: str | None = None
    binder2_start: int | None = None
    label: str = ""
