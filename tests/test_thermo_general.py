"""General-engine tests: arbitrary parameters, no ECLIPSE specifics."""

import numpy as np
import pytest

from lockr.engine import thermo
from lockr.engine.models import SensorParams


@pytest.mark.parametrize("Kd,target", [(1e-9, 1e-6), (1e-12, 1e-9), (5e-10, 1e-7)])
def test_fold_change_monotonic_in_pull(Kd, target):
    params = SensorParams(K_open=5e-3, K_CK=2e-8, lucKey=200e-9)
    fcs = [thermo.fold_change(target, Kd, p, params) for p in (1, 5, 10, 20, 40)]
    assert fcs == sorted(fcs)


@pytest.mark.parametrize("params", [
    SensorParams(K_open=1e-3, K_CK=1e-8, lucKey=500e-9),
    SensorParams(K_open=1e-3, K_CK=1e-8, lucKey=1e-9),
    SensorParams(K_open=2e-2, K_CK=5e-7, lucKey=1e-6),
])
def test_regime_classification_runs_for_arbitrary_params(params):
    r = thermo.diagnose_regime(params, pull=15)
    assert r.regime in ("key-limited", "K_open-limited", "mixed")
    assert r.luckey_dominance_ratio == pytest.approx(params.luckey_ratio)


def test_regime_flips_as_ratio_crosses_k_open():
    high_ratio = SensorParams(K_open=1e-3, K_CK=1e-9, lucKey=1e-6)    # ratio 1000
    low_ratio = SensorParams(K_open=1e-3, K_CK=1e-9, lucKey=1e-12)    # ratio 0.001
    assert thermo.diagnose_regime(high_ratio).regime == "key-limited"
    assert thermo.diagnose_regime(low_ratio).regime == "K_open-limited"


def test_max_fold_change_independent_of_kd():
    params = SensorParams(K_open=1e-3, K_CK=3e-8, lucKey=400e-9)
    a = thermo.max_fold_change(1e-9, 12, params)
    b = thermo.max_fold_change(1e-15, 12, params)
    assert a == pytest.approx(b)


def test_fit_pull_strength_recovers_arbitrary_value():
    params = SensorParams(K_open=2e-3, K_CK=4e-8, lucKey=300e-9)
    Kd = 5e-10
    conc = np.logspace(-13, -6, 8)
    true_pull = 23.0
    fc = [thermo.fold_change(c, Kd, true_pull, params) for c in conc]
    pull, _ = thermo.fit_pull_strength(conc, fc, Kd, params)
    assert pull == pytest.approx(true_pull, rel=1e-3)
