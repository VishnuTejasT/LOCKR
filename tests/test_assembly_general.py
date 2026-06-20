"""General-engine tests: synthetic motifs/positions, no ECLIPSE specifics."""

import dataclasses

from lockr.engine import assembly
from lockr.engine.models import GraftSpec, LatchWindow, ProtectedRegion, VariantSuggestion


def test_check_protected_region_detects_intact_motif():
    full = "AAAAMOTIFAAAA"
    r = assembly.check_protected_region(full, "MOTIF", 5, 9)
    assert r.intact
    assert r.found_sequence == "MOTIF"
    assert r.mismatch_positions == []


def test_check_protected_region_detects_single_substitution():
    full = "AAAAMOTXFAAAA"   # I -> X at position 8
    r = assembly.check_protected_region(full, "MOTIF", 5, 9)
    assert not r.intact
    assert r.found_sequence == "MOTXF"
    assert r.mismatch_positions == [8]


def test_check_graft_overlap_flags_overlap():
    region = ProtectedRegion(motif="MOTIF", start=5, end=9)
    graft = GraftSpec(binder="WXYZ", start=8)   # 8-11 overlaps 5-9 at 8,9
    r = assembly.check_graft_overlap(graft, region)
    assert r.overlap
    assert r.overlapping_positions == [8, 9]


def test_check_graft_overlap_clears_non_overlapping_graft():
    region = ProtectedRegion(motif="MOTIF", start=5, end=9)
    graft = GraftSpec(binder="WXYZ", start=20)
    r = assembly.check_graft_overlap(graft, region)
    assert not r.overlap
    assert r.overlapping_positions == []


def test_check_latch_fit_accepts_graft_that_fits():
    window = LatchWindow(start=1, end=10)
    graft = GraftSpec(binder="ABCDE", start=1)   # 5aa in a 10aa window
    r = assembly.check_latch_fit(graft, window)
    assert r.fits
    assert r.used_length == 5
    assert r.available_length == 10
    assert r.slack == 5


def test_check_latch_fit_flags_graft_that_overflows():
    window = LatchWindow(start=1, end=4)
    graft = GraftSpec(binder="ABCDE", start=1)   # 5aa in a 4aa window
    r = assembly.check_latch_fit(graft, window)
    assert not r.fits
    assert r.used_length == 5
    assert r.available_length == 4
    assert r.slack == -1


def test_check_latch_fit_counts_linker_and_binder2():
    window = LatchWindow(start=1, end=11)
    graft = GraftSpec(binder="ABCDE", start=1, linker="X", linker_start=6,
                      binder2="FGHIJ", binder2_start=7)
    r = assembly.check_latch_fit(graft, window)
    assert r.used_length == 11   # 5 + 1 + 5
    assert r.fits
    assert r.slack == 0


def test_verify_full_assembly_single_binder_all_pass():
    full = "AAAAMOTIFAAWXYZAAAA"
    region = ProtectedRegion(motif="MOTIF", start=5, end=9)
    window = LatchWindow(start=12, end=15)
    graft = GraftSpec(binder="WXYZ", start=12)
    r = assembly.verify_full_assembly(full, window, graft, region,
                                      expected_total_length=len(full))
    names = [c.name for c in r.checks]
    assert names == ["overall_length", "protected_region_intact",
                     "graft_no_overlap", "latch_fit", "binder1_intact"]
    assert r.all_passed


def test_verify_full_assembly_tandem_adds_linker_and_binder2_checks():
    full = "AAAAMOTIFAA" + "WXYZ" + "L" + "PQRS" + "AAAA"
    region = ProtectedRegion(motif="MOTIF", start=5, end=9)
    window = LatchWindow(start=12, end=20)
    graft = GraftSpec(binder="WXYZ", start=12, linker="L", linker_start=16,
                      binder2="PQRS", binder2_start=17)
    r = assembly.verify_full_assembly(full, window, graft, region)
    names = [c.name for c in r.checks]
    assert "linker_intact" in names
    assert "binder2_intact" in names
    assert r.all_passed


def test_verify_full_assembly_flags_a_broken_check_without_failing_the_rest():
    full = "AAAAMOTIFAAWXYZAAAA"
    region = ProtectedRegion(motif="MOTIF", start=5, end=9)
    window = LatchWindow(start=12, end=15)
    graft = GraftSpec(binder="WRONG", start=12)   # doesn't match what's actually there
    r = assembly.verify_full_assembly(full, window, graft, region)
    by_name = {c.name: c.passed for c in r.checks}
    assert by_name["protected_region_intact"] is True
    assert by_name["binder1_intact"] is False
    assert not r.all_passed


def test_filter_safe_variants_rejects_substitution_inside_protected_region():
    region = ProtectedRegion(motif="MOTIF", start=10, end=14)
    inside = VariantSuggestion(policy="neutralizing", sequence="...",
                               mutations=["D12A"])      # 12 is inside 10-14
    outside = VariantSuggestion(policy="neutralizing", sequence="...",
                                mutations=["E20A"])      # 20 is outside

    r = assembly.filter_safe_variants([inside, outside], region)

    assert r.accepted == [outside]
    assert len(r.rejected) == 1
    rejected_variant, reason = r.rejected[0]
    assert rejected_variant is inside
    assert "position 12" in reason
    assert "protected region" in reason


def test_filter_safe_variants_uses_absolute_position_not_local_position():
    # Local position 12 lands inside this region (10-14) by coincidence -- if
    # filter_safe_variants (or a caller feeding it) skipped the local->absolute
    # offset, this would get wrongly rejected. The binder actually starts at
    # absolute position 30, so position 12 is local-only; the true absolute
    # position is 12 + (30 - 1) = 41, nowhere near 10-14.
    region = ProtectedRegion(motif="MOTIF", start=10, end=14)
    binder_start = 30
    offset = binder_start - 1

    v_local = VariantSuggestion(policy="neutralizing", sequence="...", mutations=["D12A"])

    # Confirms the coincidence: taken at face value, position 12 IS inside the
    # region, so an un-offset variant gets (correctly, for what it's given) rejected.
    naive = assembly.filter_safe_variants([v_local], region)
    assert naive.rejected and naive.rejected[0][0] is v_local

    # Once shifted to the binder's real absolute coordinates, it's accepted --
    # proving the offset, not the coincidence, is what determines the outcome.
    abs_mutations = [f"{m[0]}{int(m[1:-1]) + offset}{m[-1]}" for m in v_local.mutations]
    v_abs = dataclasses.replace(v_local, mutations=abs_mutations)

    r = assembly.filter_safe_variants([v_abs], region)
    assert r.accepted == [v_abs]
    assert r.rejected == []
