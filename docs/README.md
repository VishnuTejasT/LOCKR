# LOCKR Biosensor Design Tool — engine map

This is the guided tour. Read the per-module docs
([thermo.md](thermo.md), [charge.md](charge.md), [liability.md](liability.md),
[assembly.md](assembly.md), [models.md](models.md),
[calibration.md](calibration.md)) for the full function-by-function
breakdown — this file is for understanding how they fit together.

## How the modules relate

```
calibration.py  --(example data only, no logic)-->  liability.py
                                                          |
charge.py  -------(net_charge, used inside scan)------->  |
                                                          v
models.py  <--(shared dataclasses, no logic)--  liability.py
   ^                                                  |
   |                                          VariantSuggestion
   |  (one-way: assembly imports models;             |
   |   liability never imports assembly)              v
   +------------------------------------------  assembly.py
                                                          |
thermo.py  ---(shared dataclasses only, no                |
               direct call into/from assembly)             |
   ^                                                       |
   |                                                       |
   +----------- both feed a final design decision ---------+
                (is it safe AND does it signal well enough?)
```

In words:
- **`models.py`** has no dependents conceptually backwards — everything else
  imports dataclasses *from* it, it imports nothing domain-specific itself.
- **`charge.py`** is a leaf — general amino-acid chemistry, imports only
  `models.py`'s `ChargeResult`.
- **`liability.py`** imports `charge.py` (for `net_charge`) and
  `calibration.py` (for ECLIPSE example defaults), and produces
  `VariantSuggestion`/`LiabilityReport` (from `models.py`).
- **`calibration.py`** is pure data — no functions, no imports of engine
  logic at all.
- **`thermo.py`** is completely independent of all of the above — it never
  sees a sequence, never imports `charge.py`/`liability.py`/`assembly.py`.
  It only depends on `models.py` for its result dataclasses.
- **`assembly.py`** is the one module that crosses a boundary: it imports
  `VariantSuggestion` from `models.py` so its bridge function
  `filter_safe_variants` can consume `liability.py`'s output. **This
  dependency runs one way only.** `liability.py` has zero awareness that
  `assembly.py` or the concept of a "protected region" exists. This was
  verified directly (grepping `liability.py`'s imports before writing the
  bridge function) rather than assumed — it's a real architectural
  constraint of this codebase, not an incidental detail.
- **`thermo.py` and `assembly.py` never talk to each other at all.** They
  answer two genuinely separate questions — "does this design produce enough
  signal" vs. "did I actually build the sequence I think I built" — and a
  real design decision (like the end-to-end story below) consults both
  independently, then combines the verdicts by hand/by UI, not through any
  shared function.

## End-to-end story: ECLIPSE's original binder, start to finish

This walks through the actual sequence of function calls that takes you from
"I have a charged binder that doesn't signal well" to "here's a verified,
fixed design with a predicted fold-change," using nothing but real ECLIPSE
numbers.

**1. Start with the problem binder.**
```python
ORIGINAL = "LISDAELEAIFAEELDC"
```
This is ECLIPSE's actual original binder — known (experimentally) to not
give clean signal.

**2. Diagnose why, with `liability.scan_liability`.**
```python
from lockr.engine import liability, calibration
report = liability.scan_liability(ORIGINAL, preserve_positions=calibration.PFLDH_INTERFACE)
```
This flags 6 acidic residues (D4, E6, E8, E13, E14, D16) sitting outside the
PfLDH target interface, for a total penalty of 4.8 kcal/mol, landing in the
`"high"` liability band, with an estimated `K_CK_estimate ≈ 3.32e-5` — about
3300x weaker than the `1e-8` reference. This is the quantitative version of
"this binder's own charge is sabotaging the cage's K_CK."

**3. Get a fix, with `liability.suggest_variant`.**
```python
fix = liability.suggest_variant(ORIGINAL, preserve_positions=calibration.PFLDH_INTERFACE,
                                policy="neutralizing")
```
Returns the actual optimized sequence (`LISAAALAAIFAAALAC`), the six
mutations (`D4A, E6A, E8A, E13A, E14A, D16A` — **binder-local positions**),
and confirms the fix restores `K_CK_estimate` to `~1e-8`.

**4. Sanity-check the fix doesn't break the fold, with `charge.analyze_charge`.**
```python
from lockr.engine import charge
charge.analyze_charge(fix.sequence)
```
Confirms `helical_ok=True`, no internal P/G breakers — the fix didn't
accidentally introduce a structural problem while solving the charge problem.

**5. Verify the fix is structurally safe in the full assembly, with
`assembly.filter_safe_variants`.**
```python
import dataclasses
from lockr.engine import assembly
from lockr.engine.models import ProtectedRegion

SMBIT = ProtectedRegion(motif="VTGYRLFEEIL", start=312, end=322, label="SmBiT")
offset = 327 - 1   # the binder starts at absolute position 327 in the full assembly
abs_mutations = [f"{m[0]}{int(m[1:-1]) + offset}{m[-1]}" for m in fix.mutations]
fix_abs = dataclasses.replace(fix, mutations=abs_mutations)

result = assembly.filter_safe_variants([fix_abs], SMBIT)
# result.accepted == [fix_abs] -- none of the fix's mutations land inside SmBiT
```
**This step only means anything because of the manual offset** —
`filter_safe_variants` takes positions at face value, and `suggest_variant`'s
positions are binder-local, not absolute. Skipping this translation step
would make the check pass or fail for a coincidental reason instead of a real
one.

**6. Verify the whole assembled sequence matches its blueprint, with
`assembly.verify_full_assembly`.**
```python
from lockr.engine.models import LatchWindow, GraftSpec

LATCH = LatchWindow(start=325, end=359, expected_length=35)
graft = GraftSpec(binder=fix.sequence, start=327, spacer="DA", spacer_start=323)
checklist = assembly.verify_full_assembly(V10_FULL_SEQUENCE, LATCH, graft, SMBIT,
                                          expected_total_length=359)
# checklist.all_passed == True, 6 checks, reproducing Script 6
```

**7. Predict fold-change for both designs, with `thermo.fold_change` / `thermo.max_fold_change`.**
```python
from lockr.engine import thermo

KD_V10 = 100e-12     # v1.0's measured Kd
KD_V22 = 42.21e-15   # v2.2's Kd, from a -4.6 kcal/mol design improvement:
                     # thermo.kd_from_ddg(KD_V10, -4.6) == KD_V22

thermo.max_fold_change(KD_V10, pull=10)   # -> ~11.0
thermo.max_fold_change(KD_V22, pull=10)   # -> ~11.0  (same ceiling, set by the cage not Kd)
```
v1.0 and v2.2 hit the *same* fold-change ceiling at a given `pull` — what
v2.2's much tighter `Kd` actually buys is a lower EC50 (better sensitivity at
low target concentration), not a bigger ceiling. `thermo.diagnose_regime()`
on the ECLIPSE defaults confirms this design is `"key-limited"` — meaning
`lucKey`/`K_CK` dominates the ceiling, not `K_open`, so v2.2's affinity gain
sharpens sensitivity but doesn't raise the maximum achievable signal.

That's the full loop: liability scan → variant suggestion → structural
charge sanity check → coordinate-correct safety check against the protected
reporter → full assembly verification → fold-change prediction. Every step
uses a different module, and the only place two of those modules' outputs
directly cross is step 5 (`liability.py`'s `VariantSuggestion` flowing into
`assembly.py`'s `filter_safe_variants`) — everywhere else, a human (or a UI)
is the one stitching the modules' separate answers into one design decision.

## Mine vs. universal: the reference table

| Value | Where it lives | Mine (ECLIPSE) or universal? |
|---|---|---|
| `theta`, `k_open_eff`, `_f_open`, `fold_change` equations | `thermo.py` | Universal — general LOCKR three-state math |
| `K_open=1e-3`, `K_CK=1e-8`, `lucKey=500e-9`, `RT=0.592` | `models.py` (`DEFAULT_PARAMS`) | Mine — ECLIPSE base scaffold's defaults, not laws of physics |
| `KD_V10=100e-12`, `KD_V22=42.21e-15` | test files only, never the engine | Mine — real measured/derived v1.0 and v2.2 affinities |
| `_K_OPEN_PROBE_FACTOR=30`, regime thresholds `0.02`/`0.08` | `thermo.py` | Heuristic I picked, not biology |
| `net_charge`, `helix_propensity`, pKa/Chou-Fasman tables | `charge.py` | Universal — standard textbook amino-acid chemistry |
| Liability scan logic (`scan_liability`, `suggest_variant`) | `liability.py` | Universal mechanism — but every default it falls back to is mine |
| `PFLDH_INTERFACE=[1,2,11,12,15]` | `calibration.py` | Mine — ECLIPSE's actual RFdiffusion-designed target contacts |
| `PENALTY_PER_ACIDIC=0.8` | `calibration.py` | Mine — two-point calibration off ECLIPSE's two real binders, explicitly flagged as not-yet-refined |
| `ANCHORS` (original/optimized sequences + ground truth) | `calibration.py` | Mine — real documented ECLIPSE data |
| `_BAND_LOW_MAX=33`, `_BAND_MODERATE_MAX=66`, `score_scale=2.0` | `liability.py` | UI-gauge design choice, not biology |
| `check_protected_region`, `check_graft_overlap`, `check_latch_fit`, `verify_full_assembly`, `filter_safe_variants` | `assembly.py` | Universal — generalized bookkeeping, no motif/binder assumed |
| `SmBiT` motif/position, `LATCH` window, `BINDER` sequence | `tests/test_assembly_eclipse.py` only | Mine — real ECLIPSE reporter/latch/binder, never in the engine itself |
| `ProtectedRegion`, `LatchWindow`, `GraftSpec` shapes | `models.py` | Universal shapes — all fields are caller-supplied |

## The two science corrections in this project's history

Both of these were caught and fixed in commit `104bd9c`
("stop calling max_fold_change a 'ceiling', it's not the dominance ratio"):

1. **The K_CK documentation typo: `~3e-12` → `3.32e-5`.** My source PDFs
   originally claimed the original (charged) binder's grafted K_CK was about
   `3e-12` — which would mean it binds *more tightly* than the clean
   reference (`1e-8`), the opposite of what a penalty model predicts for a
   binder with liabilities. Hand-checking the actual math
   (`exp(-(10.905-4.8)/0.592) ≈ 3.32e-5`) confirmed `3.32e-5` is correct —
   about 3300x *weaker*, consistent with the binder's documented "doesn't
   signal" behavior. The PDFs have been corrected, and
   `test_liability_eclipse.py` pins the corrected `3.32e-5` value so nobody
   accidentally regresses the code back toward the old wrong number.

2. **The "ceiling" naming collision: `max_fold_change` ≠ `luckey_ratio`.**
   `lucKey/K_CK` was always correctly named `luckey_ratio`/
   `luckey_dominance_ratio` in the actual code — but comments and test names
   were still using the word "ceiling" loosely for a *different* number,
   `max_fold_change` (the realized fold-change at a given finite `pull`,
   `≈ 1+pull` in the key-limited regime). Same word, two genuinely different
   quantities — `luckey_dominance_ratio=50` does not mean "50x fold-change is
   achievable"; the actual achievable ceiling at `pull=10` is `~11x`. Renamed
   and tightened the comments/tests so it's unambiguous which one any given
   number refers to going forward.
