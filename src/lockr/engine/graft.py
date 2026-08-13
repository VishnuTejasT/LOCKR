"""Threading scan using PyRosetta REF2015. Calibrated on lucCage PfLB-1
v1.0 and the original charged binder, thresholds are ECLIPSE-derived
estimates, not universal constants, and are NOT portable across machines
(see the calibration block below).

Standard PyRosetta threading, not proprietary logic: mutate the binder into
each candidate position, score with ref2015, keep the lowest-scoring
position. No repacking or relax is applied, this is a fast ranking pass over
positions, not a structure-refinement step.
"""

from __future__ import annotations

import platform
import tempfile
import time
from importlib.metadata import PackageNotFoundError, version as _pkg_version

from .models import GraftAtResult, GraftResult

PYROSETTA_AVAILABLE = True
try:
    from pyrosetta import init as _pyrosetta_init
    from pyrosetta import pose_from_pdb
    from pyrosetta.rosetta.core.scoring import get_score_function
    from pyrosetta.toolbox.mutants import mutate_residue
except ImportError:
    PYROSETTA_AVAILABLE = False

_PYROSETTA_INSTALL_MESSAGE = (
    "PyRosetta is not installed. See INSTALL.md for setup: "
    "pip install pyrosetta --find-links "
    "https://west.rosettacommons.org/pyrosetta/quarterly/release"
)

# Latch geometry on the bundled template (src/lockr/data/lucCage_template_clean.pdb).
# 325-326 is the native "ER" that isn't graftable, so threading starts at 327.
LATCH_START = 325
LATCH_END = 359
SCAN_START = 327
SCAN_END = 343

# Calibrated on this machine: PyRosetta 2026.03, M1 Mac, ref2015, no
# repacking. v1.0 (LISAAALAAIFAAALAC) scored -1469.57 REU at position 327,
# the original charged binder (LISDAELEAIFAEELDC) scored -1320.88 REU at the
# same position, both matching the ECLIPSE-documented WSL2 numbers to the
# hundredth of a REU. Rosetta scores are NOT portable across platforms or
# PyRosetta builds in general, if this ever gets run on different hardware,
# re-run find_best_graft_position() on both sequences and update these three
# constants, don't assume they still hold.
_V10_SCORE = -1469.5729329275262
_ORIGINAL_SCORE = -1320.8792272918831
_TOLERANCE = 50.0
GOOD_SCORE_MAX = _V10_SCORE + _TOLERANCE
MARGINAL_SCORE_MAX = _ORIGINAL_SCORE + _TOLERANCE

# What GOOD_SCORE_MAX/MARGINAL_SCORE_MAX above were actually measured on,
# see the calibration block's comment. Prefix-matched, not exact, since
# patch-level PyRosetta versions likely don't matter but a major/minor
# mismatch might. PyPI's version string is "2026.3", not zero-padded
# "2026.03", matched that exactly here to avoid a false-positive mismatch.
_CALIBRATED_PYROSETTA_VERSION = "2026.3"
_CALIBRATED_PLATFORM = ("Darwin", "arm64")  # macOS, Apple Silicon

_initialized = False


def _calibration_mismatch_warning() -> str | None:
    """None if this looks like the machine the REU thresholds above were
    calibrated on, otherwise a warning explaining why good/marginal/poor
    verdicts here may not mean what they claim to.

    Nothing else in this module checks this before handing back a verdict,
    Rosetta scores are not portable across platforms or PyRosetta builds
    (see the module docstring), so a verdict computed on different hardware
    is silently comparing today's REU score against thresholds measured on
    someone else's machine.
    """
    mismatches = []

    system, machine = platform.system(), platform.machine()
    if (system, machine) != _CALIBRATED_PLATFORM:
        mismatches.append(f"running on {system}/{machine}, calibrated on "
                          f"{_CALIBRATED_PLATFORM[0]}/{_CALIBRATED_PLATFORM[1]}")

    try:
        installed_version = _pkg_version("pyrosetta")
    except PackageNotFoundError:
        installed_version = None
    if installed_version is not None and not installed_version.startswith(_CALIBRATED_PYROSETTA_VERSION):
        mismatches.append(f"PyRosetta {installed_version} installed, calibrated on "
                          f"{_CALIBRATED_PYROSETTA_VERSION}")

    if not mismatches:
        return None
    return ("Calibration mismatch: " + "; ".join(mismatches) + ". GOOD_SCORE_MAX/"
           "MARGINAL_SCORE_MAX were measured on a specific machine/PyRosetta build "
           "and are not portable, re-run find_best_graft_position() on the two "
           "reference binders in the module docstring and update those constants "
           "before trusting good/marginal/poor verdicts here.")


def _ensure_init() -> None:
    # Lazy, not at module import: PyRosetta's init() loads its full database
    # (multi-second cost) and this module gets imported just by importing the
    # API app, most requests never touch grafting at all.
    global _initialized
    if not PYROSETTA_AVAILABLE:
        raise ImportError(_PYROSETTA_INSTALL_MESSAGE)
    if not _initialized:
        _pyrosetta_init("-mute all")
        _initialized = True


def _verdict(score: float) -> str:
    if GOOD_SCORE_MAX is None or MARGINAL_SCORE_MAX is None:
        return "uncalibrated"
    if score <= GOOD_SCORE_MAX:
        return "good"
    if score <= MARGINAL_SCORE_MAX:
        return "marginal"
    return "poor"


def _thread(pose, position: int, binder_sequence: str, scorefxn) -> None:
    # Mutates in place, one residue at a time, starting at `position`.
    for i, aa in enumerate(binder_sequence):
        mutate_residue(pose, position + i, aa, 0.0, scorefxn)


def _save_temp_pdb(pose) -> str:
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp.close()
    pose.dump_pdb(tmp.name)
    return tmp.name


def find_best_graft_position(
    binder_sequence: str,
    template_pdb_path: str,
    latch_start: int = LATCH_START,
    latch_end: int = LATCH_END,
    scan_start: int = SCAN_START,
    scan_end: int = SCAN_END,
) -> GraftResult:
    _ensure_init()
    start_time = time.time()

    scorefxn = get_score_function()
    template = pose_from_pdb(template_pdb_path)

    binder_length = len(binder_sequence)
    all_scores: list[tuple[int, float]] = []
    for p in range(scan_start, scan_end + 1):
        if p + binder_length - 1 > latch_end:
            continue
        pose = template.clone()
        _thread(pose, p, binder_sequence, scorefxn)
        all_scores.append((p, scorefxn(pose)))

    if not all_scores:
        raise ValueError(
            f"No valid graft positions for a {binder_length}aa binder in the "
            f"scan window {scan_start}-{scan_end} (latch ends at {latch_end})."
        )

    best_position, best_score = min(all_scores, key=lambda ps: ps[1])

    best_pose = template.clone()
    _thread(best_pose, best_position, binder_sequence, scorefxn)
    grafted_sequence = best_pose.sequence()
    grafted_pdb_path = _save_temp_pdb(best_pose)

    return GraftResult(
        best_position=best_position,
        best_score=best_score,
        verdict=_verdict(best_score),
        all_scores=all_scores,
        grafted_sequence=grafted_sequence,
        grafted_pdb_path=grafted_pdb_path,
        binder_length=binder_length,
        runtime_seconds=time.time() - start_time,
        calibration_warning=_calibration_mismatch_warning(),
    )


def graft_at_position(binder_sequence: str, template_pdb_path: str, position: int) -> GraftAtResult:
    _ensure_init()
    scorefxn = get_score_function()
    template = pose_from_pdb(template_pdb_path)

    pose = template.clone()
    _thread(pose, position, binder_sequence, scorefxn)
    score = scorefxn(pose)

    return GraftAtResult(
        position=position,
        score=score,
        verdict=_verdict(score),
        grafted_sequence=pose.sequence(),
        grafted_pdb_path=_save_temp_pdb(pose),
        calibration_warning=_calibration_mismatch_warning(),
    )
