"""Tests for the sequence-level structure pre-check in engine/helix.py."""

from __future__ import annotations

import pytest

from lockr.engine.helix import analyze_helix, compare_for_variant

# A designed amphipathic helix, the shape a lucCage binder is supposed to have.
DESIGNED_HELIX = "EELLKKLEELLKKLEELLKKL"
# Same length band, but polar and floppy with no hydrophobic face.
FLOPPY = "GSGSGSGSGSGSGSG"


def test_designed_helix_scores_higher_than_floppy_sequence():
    assert analyze_helix(DESIGNED_HELIX).helix_confidence > analyze_helix(FLOPPY).helix_confidence


def test_designed_helix_lands_in_likely_band():
    assert analyze_helix(DESIGNED_HELIX).band == "likely helical"


def test_glycine_serine_linker_is_not_called_helical():
    assert analyze_helix(FLOPPY).band == "unlikely helical"


def test_internal_proline_blocks_the_graft():
    report = analyze_helix("LISNAQFPNQLINQC")
    assert report.graft_blocked
    kinds = {i.kind for i in report.issues if i.severity == "blocking"}
    assert "internal_proline" in kinds


def test_proline_in_first_turn_is_tolerated():
    report = analyze_helix("PELLKKLEELLKKLEELLKKL")
    assert not report.graft_blocked


def test_sequence_under_two_turns_is_blocked():
    report = analyze_helix("LKKLE")
    assert report.graft_blocked
    assert any(i.kind == "too_short" for i in report.issues)


def test_blocking_issue_caps_confidence_below_uncertain_band():
    # A sequence that would otherwise score well still can't claim to be graftable.
    report = analyze_helix("EELLKKLEEPLKKLEELLKKL")
    assert report.graft_blocked
    assert report.band == "unlikely helical"


def test_salt_bridges_are_found_at_i_plus_3_and_i_plus_4():
    # E at 1 pairs with K at 5 (i+4); K at 5 pairs with E at 8 (i+3).
    assert (1, 5) in analyze_helix("EAAAKAAEAA").salt_bridges


class TestCyclizationEvidence:
    def test_cysteines_at_both_ends_flag_as_possibly_cyclic(self):
        evidence = analyze_helix("CGRLDEWKAAC").cyclization
        assert evidence.possibly_cyclic
        assert evidence.cysteine_positions == [1, 11]

    def test_single_cysteine_is_not_evidence_of_cyclization(self):
        assert not analyze_helix("LISNAQFQNQLINQC").cyclization.possibly_cyclic

    def test_ordinary_helix_salt_bridges_do_not_trigger_the_flag(self):
        # i,i+4 acid/base pairs are how designed helices stabilize themselves. Treating them
        # as staple evidence would fire on nearly every good binder.
        assert not analyze_helix(DESIGNED_HELIX).cyclization.possibly_cyclic


class TestVariantComparison:
    def test_lost_n_cap_is_reported(self):
        notes = compare_for_variant("DLEKLLKELAEKLK", "ALEKLLKELAEKLK")
        assert any("N-cap" in n for n in notes)

    def test_lost_salt_bridge_is_reported(self):
        notes = compare_for_variant("LKELLKKLEELLKK", "LKALLKKLEELLKK")
        assert any("Salt bridge" in n for n in notes)

    def test_neutralizing_near_the_n_terminus_reports_dipole_loss(self):
        notes = compare_for_variant("LLEQLLQQLLQQLL", "LLAQLLQQLLQQLL")
        assert any("dipole" in n for n in notes)

    def test_identical_sequences_produce_no_notes(self):
        assert compare_for_variant(DESIGNED_HELIX, DESIGNED_HELIX) == []

    def test_length_mismatch_is_ignored_rather_than_raising(self):
        assert compare_for_variant("LLEQLL", "LLAQ") == []

    def test_mid_helix_neutralization_away_from_caps_is_quiet(self):
        # D->A mid-chain with nothing to lose: alanine is the better helix former, no warning.
        notes = compare_for_variant("LLQQLLQQDLQQLL", "LLQQLLQQALQQLL")
        assert notes == []


@pytest.mark.parametrize("sequence", ["", "   "])
def test_empty_sequence_does_not_raise(sequence):
    report = analyze_helix(sequence)
    assert report.helix_confidence == 0.0
