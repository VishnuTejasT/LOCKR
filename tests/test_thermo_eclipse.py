"""ECLIPSE validation: the engine reproduces my documented numbers exactly.

My worked example/validation case, not a generalizable default, see
lockr-tool-plan.md. Calls the same general functions as test_thermo_general.py.

NOTE: these expected values changed when the _f_open partition function was
corrected (the open+lucKey-bound state's weight is k_open*luckey_ratio, a
sequential/conditional equilibrium, not k_open+luckey_ratio). The old
additive formula was in both the code and my own hand-derivation in
"ECLIPSE Thermodynamics Documentation.pdf" (Script 7), so the "~11x for
v1.0" and related numbers documented there are now stale and need
regenerating from that PDF's own source separately, this test suite can't
fix a PDF.
"""

import pytest

from lockr.engine import thermo
from lockr.engine.models import SensorParams

KD_V10 = 100e-12
KD_V22 = 42.21e-15


@pytest.mark.parametrize("pull,expected", [(10, 7.4061), (20, 10.6572), (30, 12.6234)])
def test_max_fold_change_v10(pull, expected):
    assert thermo.max_fold_change(KD_V10, pull) == pytest.approx(expected, rel=1e-4)


@pytest.mark.parametrize("pull,expected", [(10, 7.4061), (20, 10.6572), (30, 12.6234)])
def test_max_fold_change_v22(pull, expected):
    # Same max FC as v1.0 at each pull, set by the cage, not by Kd.
    assert thermo.max_fold_change(KD_V22, pull) == pytest.approx(expected, rel=1e-4)


def test_f_base_value():
    # weight_signal = K_open * luckey_ratio = 1e-3 * 50 = 0.05
    # f_base = weight_signal / (1 + K_open + weight_signal) = 0.05 / 1.051
    assert thermo.f_base() == pytest.approx(0.05 / 1.051, rel=1e-9)


def test_luckey_ratio_is_50_at_eclipse_defaults():
    assert SensorParams().luckey_ratio == pytest.approx(50.0)


def test_example_fc_at_pull10_saturating():
    fc = thermo.fold_change_detail(1e-6, KD_V10, 10)
    assert fc.theta == pytest.approx(1.0, abs=1e-3)
    assert fc.fold_change == pytest.approx(7.4057, rel=2e-3)


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


def test_regime_mixed_at_500nm_luckey():
    # Corrected finding (was "key-limited, latch tuning won't help" under
    # the old additive formula): at ECLIPSE's actual operating point, K_open
    # and lucKey/K_CK have comparable headroom, tightening either one helps
    # by a similar amount, neither dominates 2x over the other.
    r = thermo.diagnose_regime(pull=10)
    assert r.luckey_dominance_ratio == pytest.approx(50.0)
    assert r.max_fold_change == pytest.approx(7.4061, rel=1e-4)
    assert r.regime == "mixed"
    assert r.latch_tuning_helps is True


def test_regime_near_optimal_at_10nm_luckey():
    # At 10 nM lucKey (ratio=1), both axes are already close to their own
    # local optimum for pull=10, corrected model doesn't find headroom on
    # either side, unlike the old formula's "K_open-limited" verdict here.
    p = SensorParams(lucKey=10e-9)
    r = thermo.diagnose_regime(p)
    assert r.luckey_dominance_ratio == pytest.approx(1.0)
    assert r.regime == "mixed"
    assert r.latch_tuning_helps is False


def test_lod_v10_under_corrected_model():
    # Was checked against Script 7 (~10.08 pM / ~0.1008 nM), which used the
    # same additive _f_open formula this test suite corrected, that
    # reference is stale, see the module docstring.
    r = thermo.lod_and_ec50(KD_V10, pull=10)
    assert r.lod_2x * 1e12 == pytest.approx(12.83, rel=0.05)
    assert r.ec50 * 1e9 == pytest.approx(0.0676, rel=0.05)


def test_lod_v22_is_essentially_zero_binder_always_saturated():
    # v2.2's Kd is tight enough that the sensor saturates at ~any clinically
    # relevant concentration, LOD should come back tiny, not None/an error.
    r = thermo.lod_and_ec50(KD_V22, pull=10)
    assert r.lod_2x < 1e-13
    assert r.ec50 < 1e-13


def test_lod_none_when_target_assumed_saturating():
    r = thermo.lod_and_ec50(None, pull=10)
    assert r.lod_2x is None
    assert r.lod_3x is None
    assert r.ec50 is None


def test_kopen_destabilizing_mutations_substantially_hurt_at_500nm_luckey():
    # Corrected finding, opposite of the old one: at ECLIPSE's 500nM lucKey,
    # latch-destabilizing mutations (which raise K_open, per
    # k_open_from_destab) substantially HURT fold-change under the
    # corrected model, not "negligible help". Raising K_open raises the
    # dark/background population faster than the signal population here.
    baseline_fc = thermo.max_fold_change(1.0, pull=10)
    fcs = []
    for n in (1, 2, 3):
        k_open_new = thermo.k_open_from_destab(0.001, n, 1.0)  # moderate preset
        params = SensorParams(K_open=k_open_new)
        fc = thermo.max_fold_change(1.0, pull=10, params=params)
        fcs.append(fc)
        assert (baseline_fc - fc) / baseline_fc > 0.5  # each one costs >50% of fold-change

    assert fcs == sorted(fcs, reverse=True)  # more destabilizing mutations, worse fold-change

    crossover_luckey_nm = 0.001 * 10.0  # K_open * K_CK(nM), unit arithmetic, not model-dependent
    assert crossover_luckey_nm == pytest.approx(0.01)
