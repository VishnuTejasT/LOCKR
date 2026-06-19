"""ECLIPSE validation: my real binder pair, called through the general API."""

import pytest

from lockr.engine import calibration, liability

ORIGINAL = "LISDAELEAIFAEELDC"
OPTIMIZED = "LISAAALAAIFAAALAC"


def test_original_flags_six_acidic_at_documented_positions():
    r = liability.scan_liability(ORIGINAL, preserve_positions=calibration.PFLDH_INTERFACE)
    assert [l.position for l in r.liabilities] == [4, 6, 8, 13, 14, 16]
    assert r.penalty_total == pytest.approx(4.8)
    assert r.liability_band == "high"


def test_original_kck_estimate():
    # FLAG: my source PDFs once printed this grafted K_CK as ~3e-12; confirmed
    # documentation typo, now fixed in the PDFs. exp(-(|dG_CK|-penalty)/RT)
    # with penalty=4.8 gives ~3.3e-5, matching the docs' own "3000x weaker"
    # prose -- asserting the corrected value, not the old ~3e-12.
    r = liability.scan_liability(ORIGINAL, preserve_positions=calibration.PFLDH_INTERFACE)
    assert r.K_CK_estimate == pytest.approx(3.3e-5, rel=0.05)


def test_optimized_has_no_liabilities_and_restores_kck():
    r = liability.scan_liability(OPTIMIZED, preserve_positions=calibration.PFLDH_INTERFACE)
    assert r.liabilities == []
    assert r.K_CK_estimate == pytest.approx(1e-8, rel=1e-6)


def test_neutralizing_suggestion_reproduces_eclipse_fix():
    v = liability.suggest_variant(ORIGINAL, preserve_positions=calibration.PFLDH_INTERFACE,
                                  policy="neutralizing")
    assert v.sequence == OPTIMIZED
    assert v.mutations == ["D4A", "E6A", "E8A", "E13A", "E14A", "D16A"]


def test_interface_acidic_is_preserved():
    seq = "DISAAALAAIFAAALAC"   # D parked on a preserved interface position
    r = liability.scan_liability(seq, preserve_positions=calibration.PFLDH_INTERFACE)
    assert r.liabilities == []
