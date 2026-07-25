"""General engine tests: shape and behavior, not tied to ECLIPSE-specific numbers."""

from importlib.resources import files

import pytest

from lockr.engine import graft

TEMPLATE = str(files("lockr.data").joinpath("lucCage_template_clean.pdb"))


@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
@pytest.mark.slow
def test_find_best_graft_position_returns_valid_shape():
    r = graft.find_best_graft_position("AAAAA", TEMPLATE)
    assert graft.SCAN_START <= r.best_position <= graft.SCAN_END
    assert isinstance(r.best_score, float)
    assert r.best_score < 0
    assert r.verdict in ("good", "marginal", "poor")
    assert len(r.grafted_sequence) == 359
    assert r.binder_length == 5
    assert r.runtime_seconds > 0


@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
@pytest.mark.slow
def test_find_best_graft_position_skips_positions_that_overflow_latch():
    # A 17aa binder can't start past position 343 without overflowing the
    # latch (343 + 17 - 1 == 359), so nothing past 343 should be scanned.
    r = graft.find_best_graft_position("LISAAALAAIFAAALAC", TEMPLATE)
    assert all(p <= 343 for p, _ in r.all_scores)


@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
def test_find_best_graft_position_raises_when_binder_too_long_for_window():
    with pytest.raises(ValueError):
        graft.find_best_graft_position("A" * 40, TEMPLATE)


def test_pyrosetta_unavailable_raises_clean_import_error(monkeypatch):
    monkeypatch.setattr(graft, "PYROSETTA_AVAILABLE", False)
    with pytest.raises(ImportError, match="PyRosetta is not installed"):
        graft.find_best_graft_position("AAAAA", TEMPLATE)
