"""General-engine tests: arbitrary sequences, no ECLIPSE specifics."""

import pytest

from lockr.engine import liability


def test_excludes_arbitrary_preserve_positions():
    seq = "MDEKDEKDEK"   # D/E at 2,3,5,6,8,9
    r = liability.scan_liability(seq, preserve_positions=[5, 6])
    assert [l.position for l in r.liabilities] == [2, 3, 8, 9]


def test_empty_preserve_positions_flags_everything():
    r = liability.scan_liability("DEDE", preserve_positions=[])
    assert [l.position for l in r.liabilities] == [1, 2, 3, 4]


def test_no_acidic_residues_gives_zero_penalty():
    r = liability.scan_liability("AAAAKKKK", preserve_positions=[])
    assert r.liabilities == []
    assert r.penalty_total == 0.0
    assert r.liability_band == "Low"


def test_custom_penalty_scales_kck_estimate():
    seq = "DEDE"
    low = liability.scan_liability(seq, preserve_positions=[], charge_penalty_per_residue=0.1)
    high = liability.scan_liability(seq, preserve_positions=[], charge_penalty_per_residue=2.0)
    assert high.K_CK_estimate > low.K_CK_estimate


def test_custom_kck_reference_shifts_baseline():
    # No liabilities -> K_CK_estimate should equal the reference exactly.
    r = liability.scan_liability("AAAA", preserve_positions=[], K_CK_reference=5e-9)
    assert r.K_CK_estimate == pytest.approx(5e-9, rel=1e-6)


def test_window_restricts_scan_to_arbitrary_range():
    r = liability.scan_liability("DEDEDE", preserve_positions=[], window=(1, 3))
    assert [l.position for l in r.liabilities] == [1, 2, 3]


def test_neutralizing_and_conservative_policies_on_synthetic_sequence():
    seq = "ADE"
    neutral = liability.suggest_variant(seq, preserve_positions=[], policy="neutralizing")
    conserved = liability.suggest_variant(seq, preserve_positions=[], policy="conservative")
    assert neutral.sequence == "AAA"
    assert conserved.sequence == "ANQ"


def test_unknown_policy_raises():
    with pytest.raises(ValueError):
        liability.suggest_variant("ADE", preserve_positions=[], policy="bogus")
