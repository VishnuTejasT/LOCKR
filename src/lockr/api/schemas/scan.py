"""Request/response shapes for POST /scan and POST /suggest, per spec 6.1/6.4."""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator

from .common import Window, validate_sequence

_POLICIES = ("conservative", "neutralizing")


class SequenceInput(BaseModel):
    id: str
    sequence: str

    @field_validator("sequence")
    @classmethod
    def _valid_sequence(cls, v):
        return validate_sequence(v)


class ScanRequest(BaseModel):
    sequences: list[SequenceInput]
    sensitive_window: Window
    ph: float = 7.4
    substitution_policy: str = "conservative"
    preserve_positions: list[int] = []

    @field_validator("substitution_policy")
    @classmethod
    def _known_policy(cls, v):
        if v not in _POLICIES:
            raise ValueError(f"substitution_policy must be one of {_POLICIES}")
        return v

    @field_validator("sequences")
    @classmethod
    def _at_least_one(cls, v):
        if not v:
            raise ValueError("sequences must contain at least one entry")
        return v


class AcidicResidue(BaseModel):
    position: int
    residue: str
    in_window: bool
    contribution: float


class PerPosition(BaseModel):
    position: int
    residue: str
    contribution: float


class HelixFlag(BaseModel):
    position: int
    issue: str


class Substitution(BaseModel):
    position: int
    from_: str = Field(alias="from")
    to: str

    model_config = {"populate_by_name": True}


class SuggestedVariant(BaseModel):
    sequence: str
    substitutions: list[Substitution]
    liability_score: float
    liability_band: str
    estimated_kck_nm: float


class KckPenalty(BaseModel):
    band: str
    note: str


class ScanResultItem(BaseModel):
    id: str
    sequence: str
    length: int
    net_charge: float
    acidic_residues: list[AcidicResidue]
    liability_score: float
    liability_band: str
    predicted_kck_penalty: KckPenalty
    per_position: list[PerPosition]
    helix_flags: list[HelixFlag]
    suggested_variants: list[SuggestedVariant]


class ScanResponse(BaseModel):
    results: list[ScanResultItem]


class SuggestRequest(BaseModel):
    sequence: str
    sensitive_window: Window
    substitution_policy: str = "neutralizing"
    max_variants: int = 1
    preserve_positions: list[int] = []

    @field_validator("sequence")
    @classmethod
    def _valid_sequence(cls, v):
        return validate_sequence(v)

    @field_validator("substitution_policy")
    @classmethod
    def _known_policy(cls, v):
        if v not in _POLICIES:
            raise ValueError(f"substitution_policy must be one of {_POLICIES}")
        return v


class SuggestResponse(BaseModel):
    suggested_variants: list[SuggestedVariant]
