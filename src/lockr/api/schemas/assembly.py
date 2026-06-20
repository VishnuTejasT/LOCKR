"""Request/response shapes for POST /verify-assembly -- no spec entry, designed
from scratch and sanity-checked with the user before being wired up."""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator


class LatchWindowIn(BaseModel):
    start: int
    end: int
    expected_length: int | None = None


class GraftSpecIn(BaseModel):
    binder: str
    start: int
    spacer: str | None = None
    spacer_start: int | None = None
    linker: str | None = None
    linker_start: int | None = None
    binder2: str | None = None
    binder2_start: int | None = None


class ProtectedRegionIn(BaseModel):
    motif: str
    start: int
    end: int
    label: str = ""


class SubstitutionIn(BaseModel):
    position: int
    from_: str = Field(alias="from")
    to: str

    model_config = {"populate_by_name": True}


class CandidateVariant(BaseModel):
    sequence: str
    substitutions: list[SubstitutionIn]
    liability_score: float
    liability_band: str
    estimated_kck_nm: float


class VerifyAssemblyRequest(BaseModel):
    full_sequence: str
    latch_window: LatchWindowIn
    graft_spec: GraftSpecIn
    protected_region: ProtectedRegionIn
    expected_total_length: int | None = None
    candidate_variants: list[CandidateVariant] = []
    binder_offset: int = 0

    @field_validator("full_sequence")
    @classmethod
    def _non_empty(cls, v):
        if not v.strip():
            raise ValueError("full_sequence is empty")
        return v.strip().upper()


class CheckResult(BaseModel):
    name: str
    passed: bool
    detail: str


class RejectedVariant(BaseModel):
    variant: CandidateVariant
    reason: str


class VariantScreen(BaseModel):
    accepted: list[CandidateVariant]
    rejected: list[RejectedVariant]


class VerifyAssemblyResponse(BaseModel):
    all_passed: bool
    checks: list[CheckResult]
    variants: VariantScreen | None = None
