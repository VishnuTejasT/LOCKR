"""Reproduce the documented ECLIPSE thermodynamic numbers (Thermo doc S6-S12).

These are the anchor results the engine has to hit exactly — the 11/21/31x
ceiling, the =50 competition ratio, the v2.2 Kd, and the regime flip.
"""

import math

import pytest

from lockr.engine import thermo
from lockr.engine.models import SensorParams

Kd_V10 = 100e-12
Kd_V22 = 42.21e-15


@pytest.mark.parametrize("pull,expected", [(10, 11), (20, 21), (30, 31)])
def test_fold_change_ceiling_v10(pull, expected):
    assert round(thermo.max_fold_change(Kd_V10, pull)) == expected


@pytest.mark.parametrize("pull,expected", [(10, 11), (20, 21), (30, 31)])
def test_fold_change_ceiling_v22(pull, expected):
    # Same ceiling as v1.0 — it's set by the cage (K_open/K_CK), not by Kd.
    assert round(thermo.max_fold_change(Kd_V22, pull)) == expected


def test_f_base_value():
    assert thermo.f_base() == pytest.approx(1e-3 / 51.001, rel=1e-9)
    assert thermo.f_base() == pytest.approx(1.96e-5, rel=1e-2)


def test_luckey_ratio_is_50():
    assert SensorParams().luckey_ratio == pytest.approx(50.0)


def test_example_fc_at_pull10_saturating():
    # ~11x at pull=10 with theta ~ 1 (saturating PfLDH).
    fc = thermo.fold_change_detail(1e-6, Kd_V10, 10)
    assert fc.theta == pytest.approx(1.0, abs=1e-3)
    assert fc.fold_change == pytest.approx(11.0, rel=2e-3)


def test_kd_v22_from_ddg():
    Kd = thermo.kd_from_ddg(Kd_V10, -4.6)
    assert Kd * 1e15 == pytest.approx(42.21, rel=2e-3)   # fM


def test_kd_improvement_2369x():
    Kd = thermo.kd_from_ddg(Kd_V10, -4.6)
    assert Kd_V10 / Kd == pytest.approx(2369, rel=2e-3)


def test_dg_open_cost():
    assert thermo.dg_open_cost() == pytest.approx(4.09, abs=0.01)


def test_dg_luckey():
    assert thermo.dg_luckey() == pytest.approx(-2.32, abs=0.01)


def test_dg_kd_roundtrip():
    # Kd = exp(dG/RT) and dG = RT ln(Kd) are inverses.
    assert thermo.kd_from_dg(thermo.dg_from_kd(Kd_V10)) == pytest.approx(Kd_V10)


def test_regime_key_limited_at_500nM():
    # ratio 50: key-limited, and v1.1-style K_open bumps don't move the ceiling.
    r = thermo.diagnose_regime(pull=10)
    assert r.luckey_dominance_ratio == pytest.approx(50.0)
    assert r.max_fold_change == pytest.approx(11.0, rel=2e-3)  # realised, not 50
    assert r.regime == "key-limited"
    assert r.latch_tuning_helps is False

    # FC stays ~11/21/31x whether K_open is 1e-3 or the v1.1 "moderate" 2.93e-2.
    bumped = SensorParams(K_open=2.93e-2)
    for pull, ceil in [(10, 11), (20, 21), (30, 31)]:
        base = round(thermo.max_fold_change(Kd_V10, pull))
        bump = round(thermo.max_fold_change(Kd_V10, pull, bumped))
        assert base == ceil
        assert abs(bump - ceil) <= 1   # essentially unchanged


def test_regime_flips_to_kopen_limited_at_10nM():
    # Drop lucKey so [lucKey]/K_CK = 1; now the latch equilibrium matters.
    p = SensorParams(lucKey=10e-9)   # K_CK stays 1e-8 -> ratio 1
    r = thermo.diagnose_regime(p, pull=10)
    assert r.luckey_dominance_ratio == pytest.approx(1.0)
    assert r.regime == "K_open-limited"
    assert r.latch_tuning_helps is True


def test_fit_pull_recovers_known_value():
    # Synthesize a curve at pull=18, confirm the fit gets it back.
    import numpy as np
    conc = np.array([1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6])
    true_pull = 18.0
    fc = [thermo.fold_change(c, Kd_V10, true_pull) for c in conc]
    pull, _ = thermo.fit_pull_strength(conc, fc, Kd=Kd_V10)
    assert pull == pytest.approx(true_pull, rel=1e-3)
