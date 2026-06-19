"""Liability scan + K_CK penalty, anchored on the two ECLIPSE binders."""

import pytest

from lockr.engine import liability
from lockr.engine.models import SensorParams

ORIGINAL = "LISDAELEAIFAEELDC"
OPTIMIZED = "LISAAALAAIFAAALAC"


def test_original_flags_six_acidic_at_expected_positions():
    r = liability.scan(ORIGINAL)
    assert [l.position for l in r.liabilities] == [4, 6, 8, 13, 14, 16]
    assert {l.residue for l in r.liabilities} == {"D", "E"}
    assert len(r.liabilities) == 6


def test_original_penalty_is_4p8():
    r = liability.scan(ORIGINAL)
    assert r.penalty_total == pytest.approx(4.8)
    assert r.liability_band == "high"


def test_original_kck_collapses_3000x_weaker():
    # Formula: K_CK_grafted = exp(-(|dG_CK| - penalty)/RT) = K_CK * exp(penalty/RT).
    # 4.8 kcal/mol -> ~3.3e-5 M, ~3300x weaker than the 1e-8 starting point.
    r = liability.scan(ORIGINAL)
    assert r.K_CK_grafted == pytest.approx(3.3e-5, rel=0.05)
    assert r.K_CK_grafted / 1e-8 == pytest.approx(3300, rel=0.05)


def test_original_kck_is_not_the_pdf_typo():
    # The PDFs printed ~3e-12 for the grafted K_CK; that's a confirmed doc typo
    # (would be ~3300x tighter, opposite of repulsion). Guard against anyone
    # "correcting" the formula back toward it.
    r = liability.scan(ORIGINAL)
    assert r.K_CK_grafted > 1e-6   # weaker than start, nowhere near 3e-12


def test_optimized_has_no_liabilities_and_restores_kck():
    r = liability.scan(OPTIMIZED)
    assert r.liabilities == []
    assert r.penalty_total == 0.0
    assert r.liability_band == "low"
    assert r.K_CK_grafted == pytest.approx(1e-8, rel=1e-6)   # restored


def test_neutralizing_suggestion_reproduces_eclipse_fix():
    v = liability.suggest_variant(ORIGINAL, policy="neutralizing")
    assert v.sequence == OPTIMIZED
    assert v.mutations == ["D4A", "E6A", "E8A", "E13A", "E14A", "D16A"]
    assert v.liability_band == "low"


def test_conservative_suggestion_uses_n_and_q():
    v = liability.suggest_variant(ORIGINAL, policy="conservative")
    assert v.sequence == "LISNAQLQAIFAQQLNC"
    assert v.mutations == ["D4N", "E6Q", "E8Q", "E13Q", "E14Q", "D16N"]


def test_interface_acidic_is_preserved():
    # A D parked on a PfLDH-contact position must not be flagged or mutated.
    seq = "DISAAALAAIFAAALAC"   # D at position 1 (in PFLDH_INTERFACE)
    r = liability.scan(seq)
    assert r.liabilities == []
    v = liability.suggest_variant(seq, policy="neutralizing")
    assert v.sequence == seq


def test_window_limits_scan():
    # Restrict the sensitive window so only the first acidic residue counts.
    r = liability.scan(ORIGINAL, window=(1, 5))
    assert [l.position for l in r.liabilities] == [4]
