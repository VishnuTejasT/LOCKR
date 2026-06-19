# LOCKR Biosensor Design Tool

A CPU-only design tool for LOCKR-style biosensors, built for iGEM 2026 (ECLIPSE,
Denmark High School). It answers the two questions that nearly broke ECLIPSE:

1. **Will this binder sequence quietly destroy cage-key affinity (K_CK)?** —
   the charge-liability problem (six D/E residues killed K_CK on the v1 binder).
2. **What fold-change can I expect, and will tuning the latch help?** — the
   thermodynamic-ceiling problem (fold-change capped by [lucKey]/K_CK, not K_open).

Both turn on the same physics, so they share one thermodynamic engine.

## Status

Phase 1 — the engine — only. CLI, FastAPI backend and React UI come later
(see `lockr-tool-plan.md`).

```
src/lockr/engine/
  models.py        data objects (SensorParams, BinderSequence, results)
  thermo.py        fold-change model, free energies, regime diagnostic, pull fit
  charge.py        net charge at pH + helix-propensity flag
  liability.py     acidic-residue liability scan + K_CK penalty + variant fixes
  calibration.py   two-point anchor on the original vs optimized binder
```

The model architecture (three-state cage/key equilibrium, K_open/K_CK) follows
Langan et al. 2019 (Nature 572) and Quijano-Rubio et al. 2021 (Nature 591). The
closed-form fold-change expression and the K_CK penalty model are ECLIPSE's own
instantiation of that framework.

## Install & test

Runs in the `igem` conda env (Python 3.10).

```bash
conda activate igem
pip install -e .
pytest -q
```
