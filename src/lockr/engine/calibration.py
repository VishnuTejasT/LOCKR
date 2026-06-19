"""K_CK liability calibration, anchored on the two real ECLIPSE binders.

ECLIPSE-specific example/validation data only — the general scan logic lives
in liability.py and takes none of this as a hardcoded default. Two-point
anchor, not a fit; refine charge_penalty_per_residue as more LOCKR systems are
tested (TSA/BLI data on v1.0/v2.2 will be the first such refinement).
"""

from __future__ import annotations

# Example only: ECLIPSE's PfLDH-contact positions (RFDiffusion design, 1-indexed).
# Pass your own target-interface positions to scan_liability — this is not a default.
PFLDH_INTERFACE = [1, 2, 11, 12, 15]

# Per-residue electrostatic penalty, calibrated on the two binders below.
# TODO: recalibrate as more LOCKR systems are tested.
PENALTY_PER_ACIDIC = 0.8

ANCHORS = {
    "original": {
        "sequence": "LISDAELEAIFAEELDC",
        "acidic_positions": [4, 6, 8, 13, 14, 16],
        "n_liabilities": 6,
        "penalty": 4.8,
        "signal": False,
    },
    "optimized": {
        "sequence": "LISAAALAAIFAAALAC",
        "acidic_positions": [],
        "n_liabilities": 0,
        "penalty": 0.0,
        "signal": True,
    },
}
