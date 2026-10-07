"""Monte Carlo over uncertain thermodynamics: determinism, direction, and honest refusals."""

import numpy as np
import pytest

from lockr.engine import landscape, temperature_mc as mc
from lockr.engine.temperature import ThermoSpec

_BASE = {"k_ck": 1e-8, "k_open": 1e-3, "luckey": 500e-9, "pull": 10.0}
_SPECS = {
    "dG_latch_cage": ThermoSpec(dH_ref=-12.0, dCp=0.0, Tm=335.0, label="dG_latch_cage"),
    "dG_key_cage": ThermoSpec(dH_ref=-8.0, label="dG_key_cage"),
    "dG_target_binder": ThermoSpec(label="dG_target_binder"),
}
_ISO = ThermoSpec(label="iso")


def _scan(**over):
    kwargs = dict(base=_BASE, specs=_SPECS, T_lo=277.15, T_hi=320.15, draws=120, seed=0)
    kwargs.update(over)
    return mc.scan(**kwargs)


def test_same_seed_reproduces_everything():
    a, b = _scan(seed=7), _scan(seed=7)
    assert a["median"] == b["median"]
    assert a["best_temperature"] == b["best_temperature"]


def test_different_seeds_differ():
    assert _scan(seed=1)["median"] != _scan(seed=2)["median"]


def test_band_is_ordered_and_brackets_the_median():
    result = _scan()
    for lo, mid, hi in zip(result["p05"], result["median"], result["p95"]):
        assert lo <= mid <= hi


def test_wider_spread_widens_the_band():
    narrow = _scan(sigmas={k: 0.0 for k in mc.SIGMA_KEYS})
    wide = _scan(sigmas={"dH_latch": 6.0})
    width = lambda r: np.mean(np.array(r["p95"], dtype=float) - np.array(r["p05"], dtype=float))
    assert width(narrow) == pytest.approx(0.0, abs=1e-9)
    assert width(wide) > width(_scan())


def test_zero_spread_matches_a_plain_landscape_slice():
    """With no uncertainty every draw is the same design, so the median is the ordinary curve."""
    result = _scan(sigmas={k: 0.0 for k in mc.SIGMA_KEYS}, draws=3)
    temperatures = np.array(result["temperatures"])
    expected, _ = landscape.evaluate({**_BASE, "temperature_K": temperatures},
                                     "fold_change", _SPECS)
    assert np.allclose(np.array(result["median"], dtype=float), expected, rtol=1e-12)


def test_best_temperature_is_inside_the_requested_range():
    result = _scan(T_lo=280.0, T_hi=310.0)
    best = result["best_temperature"]
    for key in ("p05", "median", "p95"):
        assert 280.0 <= best[key] <= 310.0
    assert 0.0 <= best["interior_fraction"] <= 1.0


def test_lower_is_better_metric_picks_a_minimum():
    result = _scan(base=dict(_BASE, k_target=1e-9), metric="limit_of_detection")
    curve = np.array([v for v in result["median"]], dtype=float)
    temperatures = np.array(result["temperatures"])
    assert result["best_temperature"]["median"] == pytest.approx(
        temperatures[int(np.nanargmin(curve))], abs=np.diff(temperatures)[0] * 2)


def test_checkpoints_cover_the_working_temperatures_in_range():
    labels = [c["label"] for c in _scan()["checkpoints"]]
    assert labels == ["fridge, 4 C", "bench, 25 C", "body, 37 C"]


def test_checkpoints_outside_the_range_are_dropped():
    assert [c["label"] for c in _scan(T_lo=295.0, T_hi=300.0)["checkpoints"]] == ["bench, 25 C"]


def test_probability_is_absent_without_a_threshold():
    assert all(c["probability"] is None for c in _scan()["checkpoints"])


def test_probability_matches_a_hand_count():
    result = _scan(threshold=8.0, op=">=", draws=200, seed=3)
    for checkpoint in result["checkpoints"]:
        assert 0.0 <= checkpoint["probability"] <= 1.0
    # An impossible bar is never cleared; a trivial one always is.
    assert all(c["probability"] == 0.0
               for c in _scan(threshold=1e9, op=">=", draws=60)["checkpoints"])
    assert all(c["probability"] == 1.0
               for c in _scan(threshold=0.0, op=">=", draws=60)["checkpoints"])


def test_lower_is_better_threshold_uses_the_op_it_was_given():
    result = _scan(base=dict(_BASE, k_target=1e-9), metric="limit_of_detection",
                   threshold=1e9, op="<=", draws=60)
    assert all(c["probability"] == 1.0 for c in result["checkpoints"])


@pytest.mark.parametrize("key", ["dG_latch_cage", "dG_key_cage"])
def test_isothermal_cage_or_key_is_refused(key):
    with pytest.raises(ValueError, match="no dH_ref"):
        _scan(specs=dict(_SPECS, **{key: _ISO}))


def test_isothermal_target_binder_is_allowed():
    """It costs fidelity, not meaning: the cage and key constants still move with temperature."""
    result = _scan(specs=dict(_SPECS, dG_target_binder=_ISO))
    assert result["draws"] > 0


def test_assumptions_are_reported_and_flagged_as_assumptions():
    result = _scan()
    assert len(result["assumptions"]) == len([v for v in mc.DEFAULT_SIGMAS.values() if v])
    assert any("not measurements" in w for w in result["warnings"])


def test_zero_spreads_are_left_out_of_the_assumption_list():
    assert mc.describe_sigmas({"dH_latch": 0.0, "dH_key": 1.5}) == [
        "key-cage heat response (dH), +-1.5 kcal/mol"]


def test_unknown_or_negative_spread_refused():
    with pytest.raises(ValueError, match="unknown spread"):
        mc.resolve_sigmas({"dH_wrong": 1.0})
    with pytest.raises(ValueError, match="must be >= 0"):
        mc.resolve_sigmas({"dH_latch": -1.0})


def test_omitted_spreads_fall_back_to_the_defaults():
    assert mc.resolve_sigmas({"dH_latch": 5.0})["dH_key"] == mc.DEFAULT_SIGMAS["dH_key"]
    assert mc.resolve_sigmas(None) == mc.DEFAULT_SIGMAS


@pytest.mark.parametrize("draws", [0, mc.MAX_DRAWS + 1])
def test_draw_count_is_bounded(draws):
    with pytest.raises(ValueError, match="draws must be between"):
        _scan(draws=draws)


def test_evaluation_cap_trims_the_temperature_grid_and_says_so():
    result = _scan(draws=mc.MAX_DRAWS, steps=mc.DEFAULT_STEPS)
    assert len(result["temperatures"]) * result["requested_draws"] <= mc.MAX_EVALUATIONS
    assert any("temperature points" in w for w in result["warnings"])


def test_unknown_metric_and_op_refused():
    with pytest.raises(ValueError, match="metric must be one of"):
        _scan(metric="glow")
    with pytest.raises(ValueError, match="op must be one of"):
        _scan(op="~=")
