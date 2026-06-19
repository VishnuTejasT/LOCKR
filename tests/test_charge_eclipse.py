"""ECLIPSE validation: net charge / helicity on my two real binders."""

import pytest

from lockr.engine import charge

ORIGINAL = "LISDAELEAIFAEELDC"
OPTIMIZED = "LISAAALAAIFAAALAC"


def test_original_is_strongly_negative():
    assert charge.net_charge(ORIGINAL, 7.4) < -5


def test_optimized_is_near_neutral():
    assert charge.net_charge(OPTIMIZED, 7.4) == pytest.approx(0.0, abs=0.5)


def test_neutralization_raises_charge():
    assert charge.net_charge(OPTIMIZED) > charge.net_charge(ORIGINAL)


def test_optimized_binder_is_helix_friendly():
    res = charge.analyze_charge(OPTIMIZED)
    assert res.helical_ok
    assert res.helix_breakers == []
