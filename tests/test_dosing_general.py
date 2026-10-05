"""General-engine tests for dosing.py: quadratic mass-balance solve, no assay specifics."""

import pytest

from lockr.engine import dosing


@pytest.mark.parametrize("key_conc,cage_total,k_ck", [
    (5e-9, 10e-9, 10e-9),
    (1e-8, 1e-8, 1e-8),   # Key == Cage == K_CK, the exact regime theta() gets wrong
    (2e-7, 5e-9, 1e-8),   # Key in real excess, quadratic should agree with theta() here
])
def test_forward_then_backward_recovers_key_conc(key_conc, cage_total, k_ck):
    target = dosing.fraction_bound(key_conc, cage_total, k_ck)
    recovered = dosing.solve_key_conc(target, cage_total, k_ck)
    assert recovered == pytest.approx(key_conc, rel=1e-6)


def test_quadratic_diverges_from_naive_approximation_when_key_and_cage_comparable():
    # Cage depletes a comparable-magnitude Key pool, so the exact fraction bound
    # must be measurably below the excess-Key (theta-style) approximation.
    key_conc, cage_total, k_ck = 1e-8, 1e-8, 1e-8
    exact = dosing.fraction_bound(key_conc, cage_total, k_ck)
    naive = key_conc / (key_conc + k_ck)
    assert exact < naive - 0.05


def test_fraction_bound_zero_at_zero_key():
    assert dosing.fraction_bound(0.0, 1e-8, 1e-8) == 0.0


def test_target_fraction_bound_near_zero_is_well_behaved():
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-4,
                                    assay_volume=1e-4, target_fraction_bound=1e-6)
    assert result.key_conc > 0
    assert result.fraction_bound == pytest.approx(1e-6, rel=1e-2)


def test_target_fraction_bound_zero_returns_zero_key_no_error():
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-4,
                                    assay_volume=1e-4, target_fraction_bound=0.0)
    assert result.key_conc == 0.0
    assert result.fraction_bound == 0.0


def test_target_fraction_bound_close_to_one_solves_without_diverging():
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-3,
                                    assay_volume=1e-4, target_fraction_bound=0.999)
    assert result.fraction_bound == pytest.approx(0.999, rel=1e-3)


def test_target_fraction_bound_of_one_raises():
    with pytest.raises(ValueError):
        dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-4,
                               assay_volume=1e-4, target_fraction_bound=1.0)


def test_fold_change_entry_point_converts_via_theta():
    # theta(fc, 1.0) = fc/(fc+1); fc=4 -> fraction 0.8, reused rather than re-derived.
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-3,
                                    assay_volume=1e-4, target_fold_change=4.0)
    assert result.fraction_bound == pytest.approx(0.8, rel=1e-3)


def test_requires_exactly_one_target_kind():
    with pytest.raises(ValueError):
        dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-4, assay_volume=1e-4)
    with pytest.raises(ValueError):
        dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-4, assay_volume=1e-4,
                               target_fraction_bound=0.5, target_fold_change=1.0)


def test_stock_too_dilute_warning_triggers():
    # A stock only 2x the assay's Cage concentration can't deliver a
    # near-saturating Key dose into a normal-sized assay volume.
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=2e-8,
                                    assay_volume=1e-4, target_fraction_bound=0.99)
    assert any("isn't concentrated enough" in w for w in result.warnings)


def test_stock_too_dilute_warning_absent_for_easy_target():
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-3,
                                    assay_volume=1e-4, target_fraction_bound=0.5)
    assert not any("isn't concentrated enough" in w for w in result.warnings)


def test_sub_pipettable_volume_warning_triggers():
    # Tiny target occupancy needs very little Key, and a very strong stock
    # shrinks the pipetted volume below a sane minimum.
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1.0,
                                    assay_volume=1e-4, target_fraction_bound=1e-6)
    assert any("pipettable minimum" in w for w in result.warnings)


def test_sub_pipettable_volume_warning_absent_for_reasonable_dose():
    # ~10 uL into a 100 uL assay, comfortably above the pipettable floor.
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1.5e-7,
                                    assay_volume=1e-4, target_fraction_bound=0.5)
    assert not any("pipettable minimum" in w for w in result.warnings)
    assert not any("isn't concentrated enough" in w for w in result.warnings)


def test_titration_curve_is_monotonic_and_spans_center():
    result = dosing.dose_for_target(k_ck=1e-8, total_cage_conc=1e-8, key_stock_conc=1e-3,
                                    assay_volume=1e-4, target_fraction_bound=0.5)
    fractions = [f for _, f in result.titration]
    assert fractions == sorted(fractions)
    assert min(fractions) < 0.5 < max(fractions)
