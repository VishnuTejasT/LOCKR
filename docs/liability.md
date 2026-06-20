# liability.py

## What this module owns, and what it doesn't

`liability.py` answers one question: *does this binder's own charge profile
threaten the cage's K_CK?* It scans a sequence for acidic residues that
aren't protecting the target interface, estimates how much they're costing
you in K_CK terms, and can propose a fixed variant. It builds on `charge.py`
(for the raw net-charge number) but adds a LOCKR-specific interpretation
layer that `charge.py` itself knows nothing about — `charge.py` would have no
idea what a "liability" or a "K_CK estimate" even means.

This module also has **zero knowledge that `assembly.py` exists.** That's a
deliberate one-way wall: `assembly.py` is allowed to import this module's
output type (`VariantSuggestion`, via `models.py`) to build its bridge
function, but `liability.py` never imports anything from `assembly.py`, and
its mutation positions are computed without any awareness of "protected
regions" or "absolute assembly coordinates" — see the Gotchas on
`suggest_variant` below, this is the single most important thing to
internalize about this file.

## General vs. ECLIPSE-specific in this file

- **General:** the scan logic (`scan_liability`), the penalty→K_CK conversion
  math (`kck_from_penalty`), the substitution policies (`neutralizing`,
  `conservative`), the scoring/banding logic. None of this assumes a specific
  binder, target, or interface.
- **ECLIPSE-specific, living only in `calibration.py` as example data, not
  hardcoded here:** `PFLDH_INTERFACE = [1,2,11,12,15]` (my actual
  target-contact positions), `PENALTY_PER_ACIDIC = 0.8` (calibrated
  specifically on my two anchor binders — see `calibration.md`), and the two
  real sequences (`ORIGINAL`/`OPTIMIZED`) used throughout the worked examples
  below.
- **UI-design choice, not biology:** `_BAND_LOW_MAX = 33`, `_BAND_MODERATE_MAX
  = 66` — these are just where a 0-100 gauge gets split into low/moderate/high
  for display purposes.

---

## `_dg_ck(K_CK_reference, RT) -> float`  *(private helper)*

**What it does:** Converts a reference K_CK value into its free-energy
equivalent.

**Why it exists:** Internal plumbing for `kck_from_penalty` — it needs a
baseline ΔG to subtract a penalty from.

**Inputs:** `K_CK_reference: float`, dimensionless-or-M depending on how
you're treating K_CK (see the open TODO note in `models.py` about K_CK's
exact definition — flagged below); `RT: float`, kcal/mol.

**Output:** `float`, kcal/mol.

**The math:** `ΔG = -RT·ln(1/K_CK_reference)`, which is algebraically the same
as `RT·ln(K_CK_reference)` — written as `1/K_CK_reference` here presumably to
keep the "this is a dissociation-constant-like quantity" framing visible at
the call site, but it computes to the identical number either way.

**Worked example:**
```python
_dg_ck(1e-8, 0.592)
# -> -10.905...  kcal/mol
```

**Gotchas:** None beyond what's in `kck_from_penalty` below — this is a tiny
private wrapper.

---

## `kck_from_penalty(penalty, K_CK_reference=K_CK_DEFAULT, RT=RT_37C) -> float`

**What it does:** Given a total electrostatic "penalty" (in kcal/mol) from
unprotected acidic residues, estimates how much weaker the binder's actual
K_CK becomes compared to a clean reference.

**Why it exists: this is the actual scientific claim of the whole module** —
that charged residues on a LOCKR binder don't just sit there neutrally, they
electrostatically interfere with the cage's own latch-key interaction,
weakening K_CK. This function is the quantitative bridge from "I counted N
acidic residues" to "here's the K_CK I'd actually expect."

**Inputs:**
- `penalty: float`, kcal/mol — total accumulated penalty (sum of per-residue
  penalties from `scan_liability`).
- `K_CK_reference: float` — the K_CK you'd expect with *zero* liabilities;
  defaults to `K_CK_DEFAULT = 1e-8` from `models.py`.
- `RT: float`, kcal/mol; defaults to `0.592` (37°C).

**Output:** `float` — the estimated K_CK after applying the penalty. Larger
number = weaker K_CK (worse).

**The math:** `K_CK_estimate = exp(-(|ΔG_CK_reference| - penalty) / RT)`.
Think of it as: start from the reference binding energy, *subtract* the
penalty (making the effective binding less favorable), then convert that
weakened energy back into an equilibrium-constant-like number. Bigger
penalty → smaller `|ΔG| - penalty` → bigger (weaker) `K_CK_estimate`.

**Worked example (real ECLIPSE numbers):**
```python
kck_from_penalty(4.8)     # the original binder's total penalty
# -> 3.32128856337083e-05

kck_from_penalty(0.0)     # no liabilities at all
# -> 9.999999999999982e-09   # recovers K_CK_DEFAULT exactly, as it should
```
With zero penalty, you get back exactly the reference K_CK (`1e-8`) — that's
a built-in sanity check this should always satisfy. With the original
binder's `4.8` kcal/mol penalty, K_CK_estimate jumps to `3.32e-5` — about
3300x weaker than the reference, which is exactly the documented "3000x
weaker" finding this codebase's test suite pins down.

**Gotchas — this is one of the two documented science corrections in this
project's history:**
My source PDFs originally printed this grafted K_CK as `~3e-12` for the
charged original binder. That number is wrong — it would mean the charged
binder binds *more tightly* than the reference, the opposite of what a
penalty model should produce. By hand: `exp(-(10.905 - 4.8)/0.592) ≈ 3.32e-5`,
not `3e-12`. The PDFs have since been corrected to `3.32e-5`, and this
codebase's tests pin the corrected value so nobody accidentally "fixes" the
code back toward the old wrong number.

---

## `_band(score) -> str`  *(private helper)*

**What it does:** Buckets a 0-100 liability score into `"low"`, `"moderate"`,
or `"high"`.

**Why it exists:** A raw score number is harder to act on at a glance than a
traffic-light label — this is purely a UI/communication convenience.

**Inputs:** `score: float`, 0-100.

**Output:** `str`, one of `"low"` (≤33), `"moderate"` (≤66), `"high"` (>66).

**Gotchas:** `33`/`66` are arbitrary thirds chosen for a gauge UI — not a
scientifically derived cutoff. Don't read "high" as meaning "binder is
definitely broken," just "worth a closer look."

---

## `_score_from_penalty(penalty, score_scale) -> float`  *(private helper)*

**What it does:** Maps an unbounded kcal/mol penalty onto a bounded 0-100
scale, for display purposes.

**Why it exists:** Penalties can in principle keep growing (more acidic
residues = more penalty, unbounded), but a UI gauge needs a fixed range.

**Inputs:** `penalty: float`, kcal/mol; `score_scale: float` — controls how
fast the score saturates toward 100 (smaller `score_scale` = saturates
faster).

**Output:** `float`, 0 to 100 (asymptotically approaches 100, never quite
reaches it for finite penalty).

**The math:** `100 * (1 - exp(-penalty/score_scale))` — an exponential
saturation curve. At `penalty=0`, score is exactly `0`. As `penalty -> ∞`,
score `-> 100`.

**Worked example (real ECLIPSE original binder, `score_scale=2.0` default):**
```python
_score_from_penalty(4.8, 2.0)
# -> 90.92820467105875
```
A penalty of 4.8 kcal/mol against a scale of 2.0 (i.e. 2.4 "scale-widths" of
penalty) lands at ~91/100 — solidly in the "high" band, matching the original
binder's documented high-liability classification.

**Gotchas:** `score_scale` is a display-tuning knob, not a measured
constant — changing it changes where on the 0-100 gauge a given real penalty
lands, without changing the underlying K_CK estimate at all (that's computed
separately, directly from `penalty`, not from this score).

---

## `scan_liability(sequence, preserve_positions, charge_penalty_per_residue=..., K_CK_reference=K_CK_DEFAULT, RT=RT_37C, window=None, ph=7.4, score_scale=2.0) -> LiabilityReport`

**What it does:** The main entry point — scans a binder sequence for acidic
residues (D/E) sitting outside the positions you've told it to leave alone,
totals up the resulting penalty, and produces a full report: which positions
are flagged, the net charge, the liability score/band, and the estimated
K_CK.

**Why it exists:** This is the actual design-review tool — "before I commit
to this binder sequence, which residues are quietly costing me K_CK, and
should I fix them?"

**Inputs:**
- `sequence`: `str` or `BinderSequence` — gets wrapped in `BinderSequence` if
  it's a plain string (which also strips/uppercases it).
- `preserve_positions: list[int]`, **1-indexed** — positions you've decided
  to protect because they contact your target (a soft, caller-supplied
  decision, not derived from anything in this file). ECLIPSE's example:
  `calibration.PFLDH_INTERFACE = [1,2,11,12,15]`.
- `charge_penalty_per_residue: float`, kcal/mol, default
  `calibration.PENALTY_PER_ACIDIC = 0.8` — cost charged per flagged acidic
  residue.
- `K_CK_reference: float`, default `K_CK_DEFAULT = 1e-8`.
- `RT: float`, kcal/mol, default `0.592`.
- `window: tuple[int,int] | None` — restrict the scan to a sub-range of
  positions (1-indexed, inclusive); `None` scans the whole sequence.
- `ph: float`, default `7.4` — passed straight through to `net_charge`.
- `score_scale: float`, default `2.0`.

**Output:** `LiabilityReport` dataclass — `binder`, `liabilities` (list of
`Liability(position, residue, weight, penalty)`), `preserve_positions`,
`net_charge`, `penalty_total`, `liability_score`, `liability_band`,
`K_CK_estimate`.

**The math/logic:** Walk every (1-indexed) position in the sequence; flag it
if (a) the residue is `D` or `E`, **and** (b) it's not in
`preserve_positions`, **and** (c) it falls inside `window` (or there's no
window restriction). Each flagged position gets a fixed
`charge_penalty_per_residue` penalty (currently every flagged residue costs
exactly the same — there's no per-residue weighting by, e.g., how exposed it
is). Sum the penalties, convert to a 0-100 score, band it, and convert to a
K_CK estimate via `kck_from_penalty`.

**Worked example (real ECLIPSE original binder, "LISDAELEAIFAEELDC"):**
```python
scan_liability("LISDAELEAIFAEELDC", preserve_positions=[1,2,11,12,15])
# LiabilityReport(
#   liabilities=[
#     Liability(position=4,  residue='D', weight=1.0, penalty=0.8),
#     Liability(position=6,  residue='E', weight=1.0, penalty=0.8),
#     Liability(position=8,  residue='E', weight=1.0, penalty=0.8),
#     Liability(position=13, residue='E', weight=1.0, penalty=0.8),
#     Liability(position=14, residue='E', weight=1.0, penalty=0.8),
#     Liability(position=16, residue='D', weight=1.0, penalty=0.8),
#   ],
#   net_charge=-6.129594662561429,
#   penalty_total=4.8,
#   liability_score=90.92820467105875,
#   liability_band='high',
#   K_CK_estimate=3.32128856337083e-05,
# )
```
Six acidic residues outside the preserved interface, each costing 0.8
kcal/mol, totaling 4.8 — exactly the documented six-liability finding for
ECLIPSE's original binder. Note positions 1, 2, 11, 12, 15 are *not* flagged
even though some may technically be near acidic territory in the raw
sequence — they're protected because they're in `preserve_positions`.

**Worked example (the optimized binder — zero liabilities):**
```python
scan_liability("LISAAALAAIFAAALAC", preserve_positions=[1,2,11,12,15])
# liabilities=[], penalty_total=0, liability_band='low',
# K_CK_estimate≈1e-8   (exactly recovers K_CK_DEFAULT)
```

**Gotchas:**
- `window`, if given, **does not change which positions count as "the
  sequence" for `preserve_positions` purposes** — `preserve_positions` is
  still interpreted against the same 1-indexed numbering as the full
  sequence, `window` just additionally restricts which of those positions
  get reported. Don't assume `window=(5,10)` means "treat position 5 as if
  it were position 1."
- An open question flagged directly in `models.py`'s comments: there's
  genuine ambiguity in my source docs about whether `K_CK` means the
  lucKey↔cage Kd or the lucKey↔SmBiT Kd — this was tracked down and
  standardized (see the project's K_CK-label resolution), but it's worth
  remembering this was once a real ambiguity, not an obviously-settled
  definition.

---

## `suggest_variant(sequence, preserve_positions, policy="neutralizing", charge_penalty_per_residue=..., K_CK_reference=K_CK_DEFAULT, RT=RT_37C, window=None, score_scale=2.0) -> VariantSuggestion`

**What it does:** Takes everything `scan_liability` flags and actually
applies a fix — substitutes every flagged acidic residue according to a
chosen policy, then re-scans the *result* so you can see the improvement.

**Why it exists:** `scan_liability` only diagnoses; this is the "now propose
the actual fixed sequence" step — literally how ECLIPSE's documented fix
(charged original → neutral optimized binder) gets reproduced
programmatically.

**Inputs:** same as `scan_liability`, plus:
- `policy: str` — `"neutralizing"` (D→A, E→A — literally the ECLIPSE fix) or
  `"conservative"` (D→N, E→Q — keeps similar shape/H-bonding capacity instead
  of going fully nonpolar). Raises `ValueError` for anything else.

**Output:** `VariantSuggestion` dataclass: `policy`, `sequence` (the new,
substituted sequence), `mutations` (list of strings like `"D4A"` — see the
Gotcha below, this is the single most important thing in this whole file),
`liability_score`/`liability_band`/`K_CK_estimate` (all computed on the *new*
sequence, after the fix).

**The math/logic:** Run `scan_liability` once to find which positions are
flagged. For each flagged position, look up its replacement from the
policy's substitution table, mutate that character in the sequence, and
record a mutation string in `"{old}{position}{new}"` format (e.g. `"D4A"` —
note this string format is exactly what `assembly.py`'s `_mutation_position`
helper parses back out by extracting digits). Then re-run `scan_liability` on
the *new* sequence to report the resulting liability score/K_CK.

**Worked example (real ECLIPSE fix, neutralizing policy):**
```python
suggest_variant("LISDAELEAIFAEELDC", preserve_positions=[1,2,11,12,15],
                policy="neutralizing")
# VariantSuggestion(
#   policy='neutralizing',
#   sequence='LISAAALAAIFAAALAC',
#   mutations=['D4A', 'E6A', 'E8A', 'E13A', 'E14A', 'D16A'],
#   liability_score=0.0, liability_band='low',
#   K_CK_estimate=9.999999999999982e-09,
# )
```
This exactly reproduces ECLIPSE's documented fix — both the resulting
sequence string and the K_CK estimate (`~1e-8`, fully restored) match the
real optimized binder.

**Worked example (conservative policy, same starting binder):**
```python
suggest_variant("LISDAELEAIFAEELDC", preserve_positions=[1,2,11,12,15],
                policy="conservative")
# sequence='LISNAQLQAIFAQQLNC'
# mutations=['D4N', 'E6Q', 'E8Q', 'E13Q', 'E14Q', 'D16N']
```
Same positions fixed, but keeping amide side chains (N/Q) instead of going to
alanine — both policies fully restore K_CK_estimate to `~1e-8`, since the
penalty model only cares whether a position is acidic, not which specific
neutral residue replaced it.

**Gotchas — the most important one in this entire module:**
The mutation position numbers in `mutations` (e.g. the `4` in `"D4A"`) are
**always relative to the standalone sequence string you passed in** — they
are *not* automatically translated to wherever that binder sits inside a
larger assembled sequence. ECLIPSE's actual 17-residue binder starts at
position 327 in the full 359aa assembly, so `"D4A"` (local position 4) is
*absolute* position `327 + (4-1) = 330` — but `suggest_variant` has no idea
the binder is embedded in anything larger, so it just reports `4`. Any caller
that needs absolute coordinates (exactly the situation `assembly.py`'s
`filter_safe_variants` bridge function is in, when checking against an
absolute-coordinate `ProtectedRegion` like SmBiT) **must apply the offset
itself** (`offset = binder_start_in_assembly - 1`, then `local_pos + offset`).
This is exactly the bug I caught myself making while writing
`test_real_eclipse_suggested_variants_never_overlap_smbit` in
`test_assembly_eclipse.py` — comparing `"D4A"` directly against SmBiT's
absolute range (312-322) would have passed the test for a completely
coincidental reason, not because the real non-overlap was actually verified.
`filter_safe_variants` itself takes positions **at face value** — it does
**no offsetting on its own** — so getting this translation right is entirely
the caller's responsibility, every time.
