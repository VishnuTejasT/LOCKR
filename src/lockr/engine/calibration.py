"""K_CK liability calibration, anchored on the two real ECLIPSE binders.

This is a physically-motivated two-point anchor, not a fit: the original binder
(6 D/E, K_CK destroyed) and the charge-optimized binder (0 acidic, K_CK restored)
pin the model. Treat the penalty-per-residue and the score scale as refinable
once BLI/TSA data on v1.0 and v2.2 comes back — present it honestly as an anchor,
not an overfit claim. See lockr-tool-plan.md §11.
"""

from __future__ import annotations

# 1-indexed PfLDH-contact positions from the RFDiffusion design. These D/E (if
# any) are load-bearing for target binding and must be preserved, never mutated.
PFLDH_INTERFACE = [1, 2, 11, 12, 15]

# Electrostatic penalty per acidic residue at the lucKey-cage interface. This is
# the whole calibration: 6 * 0.8 = 4.8 kcal/mol is what collapsed K_CK on the v1
# binder. TODO: recalibrate per-residue value once BLI data is back.
PENALTY_PER_ACIDIC = 0.8  # kcal/mol

# Score-mapping scale (penalty kcal/mol -> 0..100 band). Tuned so the v1 binder's
# 4.8 kcal/mol lands firmly in "high" and a clean binder reads "low".
SCORE_SCALE = 2.0

# Band cutoffs match the Scanner gauge (UI spec §3.4): Low 0-33, Mod 34-66, High 67+.
BAND_LOW_MAX = 33
BAND_MODERATE_MAX = 66


# The two anchor cases, kept here so tests and docs reference one source.
ANCHORS = {
    "original": {
        "sequence": "LISDAELEAIFAEELDC",
        "acidic_positions": [4, 6, 8, 13, 14, 16],
        "n_liabilities": 6,
        "penalty": 4.8,
        "signal": False,   # K_CK collapsed -> functionally dead
    },
    "optimized": {
        "sequence": "LISAAALAAIFAAALAC",
        "acidic_positions": [],
        "n_liabilities": 0,
        "penalty": 0.0,
        "signal": True,    # K_CK restored to ~1e-8
    },
}
