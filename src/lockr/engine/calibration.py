from __future__ import annotations


PFLDH_INTERFACE = [1, 2, 11, 12, 15]


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
    #PfLB-1 optimized by LiabilityScan -> PfLB-1.5!!!
    "optimized": {
        "sequence": "LISAAALAAIFAAALAC",
        "acidic_positions": [],
        "n_liabilities": 0,
        "penalty": 0.0,
        "signal": True,
    },
}
