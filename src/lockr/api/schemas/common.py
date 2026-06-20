"""Shapes shared by more than one route's schema."""

from __future__ import annotations

from pydantic import BaseModel, field_validator

_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


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
