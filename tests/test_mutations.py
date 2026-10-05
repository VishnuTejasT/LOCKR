"""Module C Part 3: required ddG -> ranked, orderable mutations."""

import pytest

from lockr.engine import mutations
from lockr.engine.calibration import PENALTY_PER_ACIDIC

# PfLB-1, the anchor sequence calibration.py is built from.
_LATCH = "LISDAELEAIFAEELDC"
_INTERFACE = [4, 6, 8, 13, 14, 16]


def test_required_ddG_is_the_gap_between_current_and_target():
    report = mutations.rank_mutations(_LATCH, current_dG=-4.09, target_dG=-5.50,
                                      interface_positions=_INTERFACE)
    assert report.required_ddG == pytest.approx(-1.41)
    assert report.current_dG == -4.09
    assert report.target_dG == -5.50
    assert report.interface_positions == _INTERFACE


def _new_residues(mutation: str) -> list[str]:
    return [part[-1] for part in mutation.split("+")]


def test_candidates_are_ranked_by_distance_from_the_requirement():
    report = mutations.rank_mutations(_LATCH, -4.09, -4.90, _INTERFACE, top_n=8)
    # Ranked to the heuristic's 0.1 kcal/mol resolution, so compare by bin, not by raw error.
    bins = [int(c.error_vs_required / 0.1) for c in report.candidates]
    assert bins == sorted(bins)
    for candidate in report.candidates:
        assert candidate.error_vs_required == pytest.approx(
            abs(candidate.ddG_estimate - report.required_ddG))
        assert candidate.resulting_dG == pytest.approx(report.current_dG + candidate.ddG_estimate)


def test_top_candidate_is_in_the_best_bin_of_the_whole_search():
    everything = mutations.rank_mutations(_LATCH, -4.09, -4.90, _INTERFACE, top_n=200)
    best = mutations.rank_mutations(_LATCH, -4.09, -4.90, _INTERFACE, top_n=1)
    best_bin = min(int(c.error_vs_required / 0.1) for c in everything.candidates)
    assert int(best.candidates[0].error_vs_required / 0.1) == best_bin


def test_removing_an_acidic_residue_tightens_the_latch_by_the_calibrated_penalty():
    # D4N: acidic liability gone, and N is a slightly worse helix former.
    report = mutations.rank_mutations(_LATCH, -4.09, -4.90, [4], top_n=40, allow_double=False)
    d4n = next(c for c in report.candidates if c.mutation == "D4N")
    assert d4n.ddG_estimate < 0
    assert d4n.ddG_estimate == pytest.approx(-PENALTY_PER_ACIDIC + (1.01 - 0.67), abs=1e-9)
    assert any("acidic liability" in note for note in d4n.notes)


def test_adding_an_acidic_residue_loosens_the_latch():
    report = mutations.rank_mutations("LISAAALAAIFAAALAC", -4.09, -3.0, [4], top_n=40,
                                      allow_double=False)
    a4e = next(c for c in report.candidates if c.mutation == "A4E")
    assert a4e.ddG_estimate > 0


def test_charge_driven_singles_are_moderate_and_everything_else_is_low():
    report = mutations.rank_mutations(_LATCH, -4.09, -4.90, [4], top_n=40, allow_double=False)
    by_mutation = {c.mutation: c for c in report.candidates}
    assert by_mutation["D4N"].confidence == "moderate"   # calibrated charge term dominates
    assert by_mutation["D4G"].confidence == "low"        # propensity guess dominates
    assert by_mutation["D4P"].confidence == "low"
    assert any("proline" in note for note in by_mutation["D4P"].notes)
    assert all(c.confidence in ("moderate", "low") for c in report.candidates)


def test_double_mutants_are_always_low_confidence_and_say_why():
    report = mutations.rank_mutations(_LATCH, -4.09, -5.50, _INTERFACE, top_n=10)
    doubles = [c for c in report.candidates if len(c.positions) == 2]
    assert doubles, "a -1.41 kcal/mol requirement needs two liabilities removed"
    for candidate in doubles:
        assert candidate.confidence == "low"
        assert any("additive" in note for note in candidate.notes)
        assert "+" in candidate.mutation


def test_singles_win_ties_against_doubles():
    # A double can always shave a few hundredths off the error, which is noise at this
    # resolution, so an -0.81 kcal/mol requirement should still come back as one mutation.
    report = mutations.rank_mutations(_LATCH, -4.09, -4.90, _INTERFACE, top_n=5)
    assert len(report.candidates[0].positions) == 1
    assert report.candidates[0].error_vs_required < 0.1


def test_disclaimer_is_always_present_and_never_calls_these_predictions():
    report = mutations.rank_mutations(_LATCH, -4.09, -4.90, _INTERFACE, top_n=3)
    assert "HEURISTIC" in report.disclaimer
    assert "not predictions" in report.disclaimer
    assert report.disclaimer in report.warnings


def test_an_impossible_requirement_says_so_instead_of_offering_a_best_effort():
    report = mutations.rank_mutations(_LATCH, -4.09, -50.0, [4], top_n=3)
    assert any("beyond anything this heuristic can reach" in w for w in report.warnings)
    assert any("directional" in w for w in report.warnings)
    assert report.candidates, "the closest candidates are still listed, just flagged"


def test_cysteine_is_never_suggested():
    report = mutations.rank_mutations(_LATCH, -4.09, -4.90, _INTERFACE, top_n=200)
    assert all("C" not in _new_residues(c.mutation) for c in report.candidates)


@pytest.mark.parametrize("kwargs", [
    dict(latch_sequence="", current_dG=-4.0, target_dG=-5.0, interface_positions=[1]),
    dict(latch_sequence=_LATCH, current_dG=-4.0, target_dG=-5.0, interface_positions=[]),
    dict(latch_sequence=_LATCH, current_dG=-4.0, target_dG=-5.0, interface_positions=[99]),
])
def test_bad_inputs_raise_value_error(kwargs):
    with pytest.raises(ValueError):
        mutations.rank_mutations(**kwargs)
