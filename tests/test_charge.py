"""Charge / helix sanity checks on the two ECLIPSE binders."""

import pytest

from lockr.engine import charge

ORIGINAL = "LISDAELEAIFAEELDC"   # 6 D/E
OPTIMIZED = "LISAAALAAIFAAALAC"  # 0 acidic


def test_original_is_strongly_negative():
    # 6 acidic residues at pH 7.4 -> well below -5; the optimized binder is ~neutral.
    q = charge.net_charge(ORIGINAL, 7.4)
    assert q < -5
    assert charge.net_charge(OPTIMIZED, 7.4) == pytest.approx(0.0, abs=0.5)


def test_neutralization_raises_charge():
    assert charge.net_charge(OPTIMIZED) > charge.net_charge(ORIGINAL)


def test_basic_residue_adds_positive_charge():
    assert charge.net_charge("K", 7.4) > charge.net_charge("A", 7.4)


def test_charge_drops_as_ph_rises():
    # More acidic at higher pH (acidic side chains deprotonate).
    assert charge.net_charge(ORIGINAL, 9.0) < charge.net_charge(ORIGINAL, 5.0)


def test_optimized_binder_is_helix_friendly():
    res = charge.analyze_charge(OPTIMIZED)
    assert res.helical_ok
    assert res.helix_breakers == []


def test_internal_proline_flagged():
    assert charge.helix_breakers("AAAPAAA") == [4]
    assert charge.helix_breakers("PAAAAAA") == []   # terminal P is fine
