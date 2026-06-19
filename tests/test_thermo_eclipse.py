"""ECLIPSE validation: the engine reproduces my documented numbers exactly.

My worked example/validation case, not a generalizable default — see
lockr-tool-plan.md. Calls the same general functions as test_thermo_general.py.
"""

import pytest

from lockr.engine import thermo
from lockr.engine.models import SensorParams

KD_V10 = 100e-12
KD_V22 = 42.21e-15


@pytest.mark.parametrize("pull,expected", [(10, 11), (20, 21), (30, 31)])
def test_max_fold_change_v10(pull, expected):
    assert round(thermo.max_fold_change(KD_V10, pull)) == expected


@pytest.mark.parametrize("pull,expected", [(10, 11), (20, 21), (30, 31)])
def test_max_fold_change_v22(pull, expected):
    # Same max FC as v1.0 at each pull -- set by the cage, not by Kd.
    assert round(thermo.max_fold_change(KD_V22, pull)) == expected


def test_f_base_value():
    assert thermo.f_base() == pytest.approx(1e-3 / 51.001, rel=1e-9)


def test_luckey_ratio_is_50_at_eclipse_defaults():
    assert SensorParams().luckey_ratio == pytest.approx(50.0)


def test_example_fc_at_pull10_saturating():
    fc = thermo.fold_change_detail(1e-6, KD_V10, 10)
    assert fc.theta == pytest.approx(1.0, abs=1e-3)
    assert fc.fold_change == pytest.approx(11.0, rel=2e-3)


def test_kd_v22_from_ddg():
    Kd = thermo.kd_from_ddg(KD_V10, -4.6)
    assert Kd * 1e15 == pytest.approx(42.21, rel=2e-3)


def test_kd_improvement_2369x():
    Kd = thermo.kd_from_ddg(KD_V10, -4.6)
    assert KD_V10 / Kd == pytest.approx(2369, rel=2e-3)


def test_dg_open_cost():
    assert thermo.dg_open_cost() == pytest.approx(4.09, abs=0.01)


def test_dg_luckey():
    assert thermo.dg_luckey() == pytest.approx(-2.32, abs=0.01)


def test_regime_key_limited_at_500nm_luckey():
    r = thermo.diagnose_regime(pull=10)
    assert r.luckey_dominance_ratio == pytest.approx(50.0)
    assert r.max_fold_change == pytest.approx(11.0, rel=2e-3)
    assert r.regime == "key-limited"
    assert r.latch_tuning_helps is False


def test_regime_flips_to_kopen_limited_at_10nm_luckey():
    p = SensorParams(lucKey=10e-9)
    r = thermo.diagnose_regime(p)
    assert r.luckey_dominance_ratio == pytest.approx(1.0)
    assert r.regime == "K_open-limited"
