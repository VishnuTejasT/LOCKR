"""ECLIPSE validation: the threading scan reproduces my documented numbers.

These are slow (real PyRosetta runs, ~5-10s each), run with the default
pytest invocation but skippable via `pytest -m "not slow"`.
"""

from importlib.resources import files

import pytest

from lockr.engine import graft

pytestmark = pytest.mark.slow

TEMPLATE = str(files("lockr.data").joinpath("lucCage_template_clean.pdb"))

V10_OPTIMIZED = "LISAAALAAIFAAALAC"
ORIGINAL_CHARGED = "LISDAELEAIFAEELDC"


@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
def test_v10_optimized_grafts_at_327_good():
    r = graft.find_best_graft_position(V10_OPTIMIZED, TEMPLATE)
    assert r.best_position == 327
    assert r.best_score == pytest.approx(-1469.57, abs=5)
    assert r.verdict == "good"
    assert len(r.grafted_sequence) == 359
    assert r.grafted_sequence[326:343] == V10_OPTIMIZED


@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
def test_original_charged_scores_worse_than_optimized():
    r = graft.find_best_graft_position(ORIGINAL_CHARGED, TEMPLATE)
    assert r.best_score == pytest.approx(-1320.88, abs=5)
    assert r.verdict in ("marginal", "poor")


@pytest.mark.skipif(not graft.PYROSETTA_AVAILABLE, reason="PyRosetta not installed")
def test_graft_at_position_327_matches_scan_best():
    r = graft.graft_at_position(V10_OPTIMIZED, TEMPLATE, 327)
    assert r.score == pytest.approx(-1469.57, abs=5)
    assert r.verdict == "good"
