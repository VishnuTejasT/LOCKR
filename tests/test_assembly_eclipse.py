"""ECLIPSE validation: my own Script 4 / Script 6 verification, generalized.

Worked example/validation case -- see tests/test_assembly_general.py for the
same functions exercised on synthetic data.
"""

import dataclasses

import pytest

from lockr.engine import assembly, calibration, liability
from lockr.engine.models import GraftSpec, LatchWindow, ProtectedRegion

V10 = (
    "SKEAAKKLQDLNIELARKLLEASTKLQRLNIRLAEALLEAIARLQELNLELVYLAVELTDPKRIRDEIKEV"
    "KDKSKEIIRRAEKEIDDAAKESKKILEEARKAIRDAAEESRKILEEGSGSGSDALDELQKLNLELAKLLLKA"
    "IAETQDLNLRAAKAFLEAAAKLQELNIRAVELLVKLTDPATIRRALEHAKRRSKEIIDEAERAIRAAKRESE"
    "RIIEEARRLIEKAKEESERIIREGSGSGDPDIKKLQDLNIELARELLRAHAQLQRLNLELLRELLRALAQLQ"
    "ELNLDLLRLASELTDPDEARKAIAVTGYRLFEEILDAERLISAAALAAIFAAALACRLIREAAAASEKISRE"
)
SMBIT = ProtectedRegion(motif="VTGYRLFEEIL", start=312, end=322, label="SmBiT")
LATCH = LatchWindow(start=325, end=359, expected_length=35)
BINDER = "LISAAALAAIFAAALAC"


def test_smbit_is_intact_in_v10():
    r = assembly.check_protected_region(V10, SMBIT.motif, SMBIT.start, SMBIT.end)
    assert r.intact
    assert r.found_sequence == "VTGYRLFEEIL"


def test_v10_binder_does_not_overlap_smbit():
    graft = GraftSpec(binder=BINDER, start=327)
    r = assembly.check_graft_overlap(graft, SMBIT)
    assert not r.overlap


def test_v10_binder_fits_latch_with_documented_slack():
    # Binder occupies 327-343 inside the 325-359 window -- 17 of 35 residues used.
    graft = GraftSpec(binder=BINDER, start=327)
    r = assembly.check_latch_fit(graft, LATCH)
    assert r.fits
    assert r.used_length == 17
    assert r.available_length == 35
    assert r.slack == 18


def test_v22_tandem_fills_latch_with_zero_slack():
    # binder + G linker + binder == 35aa, exactly the latch window (Script 5).
    graft = GraftSpec(binder=BINDER, start=325, linker="G", linker_start=342,
                      binder2=BINDER, binder2_start=343)
    r = assembly.check_latch_fit(graft, LATCH)
    assert r.fits
    assert r.used_length == 35
    assert r.slack == 0


def test_verify_full_assembly_reproduces_script6_six_checks_on_v10():
    # v1.0 has no linker/binder2, so this naturally lands on six checks --
    # length, SmBiT, overlap, fit, spacer ('DA'), binder1 -- mirroring Script 6
    # (length==359, SmBiT intact, spacer=='DA', binder1 correct) without
    # hardcoding 'DA' anywhere: it comes from graft_spec.spacer below.
    graft = GraftSpec(binder=BINDER, start=327, spacer="DA", spacer_start=323)
    r = assembly.verify_full_assembly(V10, LATCH, graft, SMBIT, expected_total_length=359)
    assert len(r.checks) == 6
    assert r.all_passed
    assert {c.name for c in r.checks} == {
        "overall_length", "protected_region_intact", "graft_no_overlap",
        "latch_fit", "spacer_intact", "binder1_intact",
    }


def test_corrupted_smbit_fails_and_reports_position():
    corrupted = V10[:311] + "X" + V10[312:]   # mutate position 312 (V -> X)
    r = assembly.check_protected_region(corrupted, SMBIT.motif, SMBIT.start, SMBIT.end)
    assert not r.intact
    assert r.mismatch_positions == [312]


def test_real_eclipse_suggested_variants_never_overlap_smbit():
    # liability.suggest_variant works on the bare 17aa binder, so its mutation
    # positions are binder-local (1-17). Shift them to absolute assembly
    # coordinates (binder starts at 327) before checking against SmBiT
    # (312-322) -- otherwise this test would pass for the wrong reason.
    v_local = liability.suggest_variant("LISDAELEAIFAEELDC",
                                        preserve_positions=calibration.PFLDH_INTERFACE,
                                        policy="neutralizing")
    offset = 327 - 1
    abs_mutations = [f"{m[0]}{int(m[1:-1]) + offset}{m[-1]}" for m in v_local.mutations]
    v_abs = dataclasses.replace(v_local, mutations=abs_mutations)

    r = assembly.filter_safe_variants([v_abs], SMBIT)
    assert r.accepted == [v_abs]
    assert r.rejected == []
