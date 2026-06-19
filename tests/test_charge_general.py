"""General-engine tests: synthetic sequences, not from ECLIPSE."""

from lockr.engine import charge


def test_acidic_residue_lowers_charge():
    assert charge.net_charge("AAADAAA", 7.4) < charge.net_charge("AAAAAAA", 7.4)


def test_basic_residue_raises_charge():
    assert charge.net_charge("AAAKAAA", 7.4) > charge.net_charge("AAAAAAA", 7.4)


def test_charge_drops_as_ph_rises_for_acidic_sequence():
    seq = "DDDDDDDD"
    assert charge.net_charge(seq, 9.0) < charge.net_charge(seq, 5.0)


def test_helix_friendly_synthetic_sequence():
    res = charge.analyze_charge("AAALAAALAAAL")
    assert res.helical_ok


def test_internal_proline_breaks_helix_flag():
    res = charge.analyze_charge("AAALPAAALAAAL")
    assert not res.helical_ok
    assert res.helix_breakers != []


def test_terminal_proline_does_not_break_helix_flag():
    assert charge.helix_breakers("PAAAAAAA") == []
