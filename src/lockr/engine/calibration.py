"""
This is a reliability check for the charge penalties and checks in liability.py, shown througb team 
ECLIPSE's binder optimizations.
"""

from __future__ import annotations

# Example only: ECLIPSE's PfLDH-contact positions (RFDiffusion design, 1-indexed).
# Pass your own target-interface positions to scan_liability — this is not a default.
PFLDH_INTERFACE = [1, 2, 11, 12, 15]

# Per-residue electrostatic penalty, calibrated on the two binders below.
# TODO: recalibrate as more LOCKR systems are tested.
PENALTY_PER_ACIDIC = 0.8

ANCHORS = {
    #PfLB-1
    "original": {
        "sequence": "LISDAELEAIFAEELDC",
        "acidic_positions": [4, 6, 8, 13, 14, 16],
        "n_liabilities": 6,
        "penalty": 4.8,
        "signal": False,
    },
    #PfLB-1 optimized by LiabilityScan -> PfLB-2!!!
    "optimized": {
        "sequence": "LISAAALAAIFAAALAC",
        "acidic_positions": [],
        "n_liabilities": 0,
        "penalty": 0.0,
        "signal": True,
    },
}
