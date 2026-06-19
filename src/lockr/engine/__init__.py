"""Shared thermodynamic engine — the single source of truth for the LOCKR math.

charge / liability / thermo all import from models so results flow between the
two analysis modules (sequence liability -> K_CK -> fold-change) without anyone
re-deriving an equation. Nothing else in the tool computes physics.
"""
