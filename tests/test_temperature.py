"""Module C Part 1: Gibbs-Helmholtz, the three input modes, and the validity ceiling."""

import math

import pytest

from lockr.engine.temperature import (MODE_FULL, MODE_ISOTHERMAL, MODE_VANT_HOFF, R_KCAL,
                                      T_REF_DEFAULT, ThermoSpec, dG_at, dG_from_k, estimate_dCp,
                                      k_from_dG, shift_constant, warnings_at, weakest_mode)

# dH_ref = -12.0 kcal/mol, dS_ref = -0.025 kcal/mol/K, dCp = -0.5 kcal/mol/K, T_ref = 298.15 K.
# Worked out by hand from dG(T) = dH - T*dS + dCp*[(T - T_ref) - T*ln(T/T_ref)], e.g. at 310 K:
#   dH - T*dS      = -12 + 7.75                      = -4.25
#   (T - T_ref)    = 11.85
#   T*ln(T/T_ref)  = 310 * ln(1.0397451)             =  12.08242
#   dCp * bracket  = -0.5 * (11.85 - 12.08242)       = +0.11621
#   dG(310)                                          = -4.13379
_FULL = ThermoSpec(dH_ref=-12.0, dS_ref=-0.025, dCp=-0.5, label="test")
_HAND_COMPUTED = {
    277.00: -4.6907204994,
    298.15: -4.5462500000,
    310.00: -4.1337848039,
    328.00: -3.0766295783,
}


@pytest.mark.parametrize("T,expected", sorted(_HAND_COMPUTED.items()))
def test_gibbs_helmholtz_matches_hand_computed_values(T, expected):
    assert dG_at(_FULL, T) == pytest.approx(expected, rel=1e-9)


def test_reference_temperature_returns_the_reference_value():
    assert dG_at(_FULL, T_REF_DEFAULT) == pytest.approx(-12.0 - T_REF_DEFAULT * -0.025)


def test_heat_capacity_term_is_what_separates_full_from_vant_hoff():
    vh = ThermoSpec(dH_ref=-12.0, dS_ref=-0.025, label="test")
    assert vh.mode is MODE_VANT_HOFF
    assert _FULL.mode is MODE_FULL
    # At T_ref both agree; away from it only the full spec carries the curvature.
    assert dG_at(vh, T_REF_DEFAULT) == pytest.approx(dG_at(_FULL, T_REF_DEFAULT))
    assert dG_at(vh, 328.0) != pytest.approx(dG_at(_FULL, 328.0))


def test_isothermal_mode_refuses_to_extrapolate():
    iso = ThermoSpec(dG_ref=-4.0, label="dG_latch_cage")
    assert iso.mode is MODE_ISOTHERMAL
    # Same value at every temperature: no scaling is even attempted.
    assert dG_at(iso, 277.0) == -4.0
    assert dG_at(iso, 328.0) == -4.0

    warns = warnings_at(iso, 310.0)
    assert len(warns) == 1
    assert "TEMPERATURE SCALING UNAVAILABLE" in warns[0]
    assert "dG_latch_cage" in warns[0]
    assert warnings_at(iso, T_REF_DEFAULT) == []


def test_isothermal_constant_is_returned_bit_identical():
    for T in (277.0, T_REF_DEFAULT, 310.0):
        assert shift_constant(1e-3, T)[0] == 1e-3


def test_shifted_constant_is_bit_identical_at_the_reference_temperature():
    spec = ThermoSpec(dH_ref=-12.0, dCp=-0.5)
    assert shift_constant(1e-3, T_REF_DEFAULT, spec)[0] == 1e-3


def test_shifted_constant_moves_with_temperature_when_dH_is_known():
    spec = ThermoSpec(dH_ref=-12.0, dCp=-0.5, label="dG_key_cage")
    k_310, warns = shift_constant(1e-8, 310.0, spec)
    assert k_310 != 1e-8
    assert warns == []
    # Round trip through dG: the shifted constant is exactly the one dG(310) implies.
    anchored_dG = dG_from_k(1e-8, T_REF_DEFAULT)
    expected = k_from_dG(dG_at(ThermoSpec(dG_ref=anchored_dG, dH_ref=-12.0, dCp=-0.5), 310.0), 310.0)
    assert k_310 == pytest.approx(expected)


def test_vant_hoff_warns_once_pushed_beyond_its_window():
    spec = ThermoSpec(dH_ref=-12.0, label="dG_latch_cage")
    assert warnings_at(spec, 305.0) == []
    warns = warnings_at(spec, 328.0)
    assert any("van't Hoff" in w and "dCp" in w for w in warns)


def test_melting_regime_warning_fires_within_ten_kelvin_of_tm():
    spec = ThermoSpec(dH_ref=-12.0, dCp=-0.5, Tm=318.0, label="dG_latch_cage")
    assert warnings_at(spec, 305.0) == []
    warns = warnings_at(spec, 310.0)
    assert any("two-state" in w and "unfolding" in w for w in warns)
    assert any("Tm" in w for w in warnings_at(spec, 330.0))


def test_dCp_estimator_is_area_driven_and_flagged_by_its_caller():
    # Murphy-Freire: 0.32*nonpolar - 0.14*polar cal/mol/K, returned in kcal/mol/K.
    assert estimate_dCp(1000.0, 500.0) == pytest.approx((0.32 * 1000 - 0.14 * 500) / 1000)
    assert estimate_dCp(0.0, 0.0) == 0.0
    with pytest.raises(ValueError):
        estimate_dCp(-1.0)


def test_weakest_mode_wins_when_parameters_disagree():
    assert weakest_mode([MODE_FULL, MODE_VANT_HOFF]) is MODE_VANT_HOFF
    assert weakest_mode([MODE_FULL, MODE_ISOTHERMAL]) is MODE_ISOTHERMAL
    assert weakest_mode([]) is MODE_ISOTHERMAL


def test_free_energy_and_constant_round_trip():
    for k in (1e-12, 1e-8, 1e-3, 0.5):
        assert k_from_dG(dG_from_k(k, 310.0), 310.0) == pytest.approx(k)
    assert dG_from_k(1.0, 310.0) == 0.0
    assert dG_from_k(math.e, 310.0) == pytest.approx(R_KCAL * 310.0)


def test_bad_inputs_are_rejected():
    with pytest.raises(ValueError):
        ThermoSpec(T_ref=0.0)
    with pytest.raises(ValueError):
        ThermoSpec(dG_ref=-4.0, Tm=-1.0)
    with pytest.raises(ValueError):
        dG_at(ThermoSpec(dG_ref=-4.0), 0.0)
    with pytest.raises(ValueError):
        dG_at(ThermoSpec(dH_ref=-12.0), 310.0)  # no dG anchor and no dS
    with pytest.raises(ValueError):
        dG_from_k(0.0, 310.0)
