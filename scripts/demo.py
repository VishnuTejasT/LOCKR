"""Run-through of the ECLIPSE validation case using the general lockr engine.

Nothing here is special-cased in the engine itself -- this script just calls
thermo/liability with ECLIPSE's own numbers (calibration.py) to show the tool
reproducing the documented results end to end.

Run: conda activate igem && python scripts/demo.py
"""

from lockr.engine import calibration, liability, thermo
from lockr.engine.models import SensorParams

KD_V10 = 100e-12      # v1.0 binder, ADCP-derived
KD_V22 = 42.21e-15    # v2.2 tandem binder

ORIGINAL = calibration.ANCHORS["original"]["sequence"]
OPTIMIZED = calibration.ANCHORS["optimized"]["sequence"]


def section(title):
    print(f"\n=== {title} ===")


def fold_change_table():
    section("Max fold-change at saturating target: v1.0 vs v2.2")
    print(f"{'pull':>6} {'v1.0 (Kd=100pM)':>18} {'v2.2 (Kd=42.21fM)':>20}")
    for pull in (10, 20, 30):
        fc_v10 = thermo.max_fold_change(KD_V10, pull)
        fc_v22 = thermo.max_fold_change(KD_V22, pull)
        print(f"{pull:>6} {fc_v10:>17.1f}x {fc_v22:>19.1f}x")
    print("Same max FC at every pull -- it's cage-set (K_open/K_CK), not Kd-set.")
    print("(Not the lucKey/K_CK dominance ratio -- see the regime diagnosis below.)")


def regime_diagnosis():
    section("Regime diagnosis at 500 nM lucKey (ECLIPSE assay default)")
    r_500 = thermo.diagnose_regime(pull=10)
    print(f"lucKey/K_CK dominance ratio : {r_500.luckey_dominance_ratio:.1f}")
    print(f"realised max fold-change    : {r_500.max_fold_change:.1f}x")
    print(f"regime                      : {r_500.regime}")
    print(f"latch tuning helps?         : {r_500.latch_tuning_helps}")
    print(f"verdict: {r_500.verdict}")

    section("Regime diagnosis at 10 nM lucKey (same K_CK, lower lucKey)")
    r_10 = thermo.diagnose_regime(SensorParams(lucKey=10e-9), pull=10)
    print(f"lucKey/K_CK dominance ratio : {r_10.luckey_dominance_ratio:.1f}")
    print(f"realised max fold-change    : {r_10.max_fold_change:.1f}x")
    print(f"regime                      : {r_10.regime}")
    print(f"latch tuning helps?         : {r_10.latch_tuning_helps}")
    print(f"verdict: {r_10.verdict}")


def liability_scan(label, sequence):
    section(f"Liability scan -- {label}: {sequence}")
    r = liability.scan_liability(sequence, preserve_positions=calibration.PFLDH_INTERFACE)
    flagged = [(l.position, l.residue) for l in r.liabilities]
    print(f"flagged liabilities  : {flagged}")
    print(f"penalty total        : {r.penalty_total:.1f} kcal/mol")
    print(f"liability score/band : {r.liability_score:.0f} ({r.liability_band})")
    print(f"K_CK estimate        : {r.K_CK_estimate:.2e} M")
    return r


def variant_suggestion():
    section("Suggested variant for the original binder (neutralizing policy)")
    v = liability.suggest_variant(ORIGINAL, preserve_positions=calibration.PFLDH_INTERFACE,
                                  policy="neutralizing")
    print(f"original  : {ORIGINAL}")
    print(f"suggested : {v.sequence}")
    print(f"mutations : {v.mutations}")
    print(f"K_CK estimate after fix : {v.K_CK_estimate:.2e} M")
    print(f"liability score/band    : {v.liability_score:.0f} ({v.liability_band})")


if __name__ == "__main__":
    fold_change_table()
    regime_diagnosis()
    liability_scan("original (unfixed)", ORIGINAL)
    liability_scan("optimized", OPTIMIZED)
    variant_suggestion()
