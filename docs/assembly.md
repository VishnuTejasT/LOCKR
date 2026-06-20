# assembly.py

## What this module owns, and what it doesn't

`assembly.py` is sequence-level bookkeeping — "did I actually build the thing
I think I built." It checks that a motif is exactly where it should be, that
a graft doesn't physically overlap something protected, that a graft fits
inside its window, and that an assembled sequence matches its blueprint
piece-by-piece. It has **no concept of free energy, K_open, K_CK, or
fold-change** — that's `thermo.py`'s job entirely. And it does no structure
prediction at all (no folding, no Rosetta, no AF2/AF3) — it's pure string and
position arithmetic on a sequence you already have.

This module generalizes the manual verification I used to do by hand on my
own ECLIPSE pipeline — specifically Complete Documentation's Script 4
(sequence validation) and Script 6 (the six-point checklist) — into functions
that work for any motif/window/graft, not just mine.

**The one-way import rule, stated explicitly:** `assembly.py` imports
`VariantSuggestion` from `models.py` (which is `liability.py`'s output type)
so the bridge function `filter_safe_variants` can consume liability
suggestions. The dependency runs **only this direction**. `liability.py`
never imports anything from `assembly.py`, and has no idea this module, or
the concept of a "protected region," exists at all. This was verified
directly (`grep -n "^from\|^import" liability.py`) before writing the bridge
function, specifically so this claim in the docs is something I checked, not
something I assumed.

## General vs. ECLIPSE-specific in this file

Every function and every dataclass (`ProtectedRegion`, `LatchWindow`,
`GraftSpec`) is fully general — none of it assumes SmBiT, a specific latch
window, or a specific binder. **All ECLIPSE-specific numbers live in the test
files, not here**: SmBiT's motif/position (`"VTGYRLFEEIL"`, 312-322), the
latch window (325-359), and the binder (`LISAAALAAIFAAALAC`) are all supplied
as arguments from `test_assembly_eclipse.py`, exactly the same way a
different LOCKR team would supply their own reporter/window/binder.

---

## `_segment(sequence, start, length) -> str`  *(private helper)*

**What it does:** Pulls out a substring at a 1-indexed, inclusive-start
position.

**Why it exists:** Tiny shared helper so every other function in this file
slices consistently — 1-indexed everywhere, matching
`BinderSequence.residues()` in `models.py` and `liability.py`'s mutation
position convention.

**Inputs:** `sequence: str`; `start: int`, 1-indexed; `length: int`, how many
residues to take.

**Output:** `str`, the extracted substring (may be shorter than `length` if
it runs off the end of `sequence`).

**The math:** `sequence[start-1 : start-1+length]` — the `-1` is the
1-indexed→0-indexed conversion every position-handling function in this file
needs.

**Worked example:**
```python
_segment("AAAAMOTIFAAAA", 5, 5)
# -> "MOTIF"
```

**Gotchas:** This silently truncates rather than erroring if `start+length`
runs past the end of `sequence` — useful for catching length mismatches via
`check_protected_region`'s comparison logic (a too-short result there *is*
the signal that something's wrong), but it means this helper alone won't
tell you "this index is out of bounds" — the caller has to notice the
returned string is shorter than expected.

---

## `check_protected_region(full_sequence, protected_motif, start, end) -> ProtectedRegionCheck`

**What it does:** Verifies an exact motif appears at an exact position range
in a full sequence, completely unaltered.

**Why it exists:** Generalizes my Script 4 check — "is SmBiT still exactly
where it's supposed to be, unmutated, after I've done a graft/substitution
elsewhere in the sequence." This is a **hard pass/fail check**, deliberately
not a similarity score — for a reporter motif like SmBiT, a single wrong
residue can kill function entirely, so there's no "mostly fine" verdict to
give.

**Inputs:**
- `full_sequence: str` — the entire assembled sequence to check.
- `protected_motif: str` — the exact sequence that's supposed to be there.
- `start: int`, `end: int` — 1-indexed, inclusive bounds of where it should
  sit.

**Output:** `ProtectedRegionCheck` dataclass: `intact: bool`,
`found_sequence: str` (whatever was actually at that range, right or wrong),
`mismatch_positions: list[int]` (empty if intact).

**The math/logic:** Slice `full_sequence[start-1:end]` (note: this uses plain
slicing directly, not `_segment`, since `end` rather than `length` is given
here). If it matches `protected_motif` exactly, done — intact. If not,
compare character-by-character over the *shorter* of the two lengths to find
substitution positions, **and** if the lengths differ, treat every position
past the shorter string's end as an additional mismatch too — a length
mismatch usually means the window itself moved (e.g. an insertion/deletion
elsewhere shifted everything), which is a more serious problem than a single
substitution, so it gets flagged comprehensively rather than just silently
truncating the comparison.

**Worked example (real ECLIPSE SmBiT, intact):**
```python
check_protected_region(V10, "VTGYRLFEEIL", 312, 322)
# ProtectedRegionCheck(intact=True, found_sequence='VTGYRLFEEIL',
#                       mismatch_positions=[])
```

**Worked example (corrupted — one residue mutated, V→X at position 312):**
```python
corrupted = V10[:311] + "X" + V10[312:]
check_protected_region(corrupted, "VTGYRLFEEIL", 312, 322)
# ProtectedRegionCheck(intact=False, found_sequence='XTGYRLFEEIL',
#                       mismatch_positions=[312])
```
Correctly isolates exactly which absolute position broke.

**Gotchas:** `mismatch_positions` are reported in **absolute coordinates**
relative to `full_sequence` (since `start`/`end` are themselves absolute) —
this is different from `liability.py`'s mutation positions, which are
binder-local. Don't assume every position number you see anywhere in this
codebase uses the same coordinate frame; always check which function
produced it.

---

## `check_graft_overlap(graft_spec, protected_region) -> OverlapCheck`

**What it does:** Checks whether a graft (binder, plus optional
linker/second binder) physically overlaps a protected region's position
range.

**Why it exists:** Generalizes the "does my binder graft step on SmBiT"
check — a structural sanity check that's purely about position ranges, with
no sequence-content comparison involved (that's `check_protected_region`'s
job). This catches a *design* mistake (placed the graft in the wrong spot)
before you even get to checking whether the protected motif's content is
correct.

**Inputs:**
- `graft_spec: GraftSpec` — the binder(s)/linker and their start positions.
- `protected_region: ProtectedRegion` — the region to avoid.

**Output:** `OverlapCheck` dataclass: `overlap: bool`,
`overlapping_positions: list[int]`.

**The math/logic:** Builds the list of `(start, end)` ranges actually
occupied by the graft via `_graft_segments` (binder, plus linker/binder2 if
present — `spacer` is deliberately excluded, since spacer is pre-existing
scaffold, not something being inserted by the graft). For each segment,
intersects it against `[protected_region.start, protected_region.end]` —
`lo = max(seg_start, region.start)`, `hi = min(seg_end, region.end)`; if
`lo <= hi`, that's a real overlap, and every position in `[lo, hi]` gets
collected.

**Worked example (real ECLIPSE — confirmed non-overlapping):**
```python
graft = GraftSpec(binder="LISAAALAAIFAAALAC", start=327)
check_graft_overlap(graft, SMBIT)   # SmBiT at 312-322
# OverlapCheck(overlap=False, overlapping_positions=[])
```
The binder occupies 327-343, SmBiT occupies 312-322 — no shared positions,
confirming the real design choice was geometrically safe.

**Worked example (synthetic, deliberate overlap):**
```python
region = ProtectedRegion(motif="MOTIF", start=5, end=9)
graft = GraftSpec(binder="WXYZ", start=8)   # occupies 8-11
check_graft_overlap(graft, region)
# OverlapCheck(overlap=True, overlapping_positions=[8, 9])
```

**Gotchas:** `_graft_segments` excludes `spacer` on purpose — if you're
checking whether a *spacer* sequence overlaps a protected region (rather
than the binder/linker), this function won't catch that; you'd need
`check_protected_region` or a manual range check on the spacer's own
position instead.

---

## `check_latch_fit(graft_spec, latch_window) -> LatchFitCheck`

**What it does:** Checks whether a graft's total length (binder + optional
linker + optional second binder) fits inside a latch window's available
length.

**Why it exists:** A purely dimensional sanity check — before worrying about
sequence content or position overlaps, does the thing you're trying to graft
in even *fit* physically? This generalizes the length-accounting I did by
hand for both the single-binder v1.0 design and the tandem v2.2 design.

**Inputs:** `graft_spec: GraftSpec`; `latch_window: LatchWindow`.

**Output:** `LatchFitCheck` dataclass: `fits: bool`, `used_length: int`,
`available_length: int`, `slack: int` (**can be negative** — that's the
signal for "overflows by this many residues").

**The math:** `used = len(binder) + len(linker or 0) + len(binder2 or 0)`;
`available = latch_window.end - latch_window.start + 1`; `slack = available -
used`; `fits = slack >= 0`. Note `spacer` is again excluded — same reasoning
as `check_graft_overlap`, spacer is scaffold, not part of what's being
grafted into the window.

**Worked example (real ECLIPSE v1.0, single binder):**
```python
graft = GraftSpec(binder="LISAAALAAIFAAALAC", start=327)   # 17aa
check_latch_fit(graft, LatchWindow(start=325, end=359))    # 35aa window
# LatchFitCheck(fits=True, used_length=17, available_length=35, slack=18)
```
17 of 35 residues used — exactly the documented slack for v1.0.

**Worked example (real ECLIPSE v2.2, tandem — exactly fills the window):**
```python
graft = GraftSpec(binder="LISAAALAAIFAAALAC", start=325,
                  linker="G", linker_start=342,
                  binder2="LISAAALAAIFAAALAC", binder2_start=343)
check_latch_fit(graft, LatchWindow(start=325, end=359))
# LatchFitCheck(fits=True, used_length=35, available_length=35, slack=0)
```
17 + 1 (G linker) + 17 = 35, exactly filling the window with zero slack —
matching Script 5's documented "v2.2 length: 359 aa OK."

**Worked example (overflow — synthetic):**
```python
check_latch_fit(GraftSpec(binder="ABCDE", start=1), LatchWindow(start=1, end=4))
# LatchFitCheck(fits=False, used_length=5, available_length=4, slack=-1)
```

**Gotchas:** `slack` being negative doesn't tell you *which* segment is too
long if there are multiple (binder/linker/binder2) — it's a total-length
check only; if you need to know which specific piece overflowed, you'd have
to inspect `graft_spec`'s individual segment lengths yourself.

---

## `verify_full_assembly(full_sequence, latch_window, graft_spec, protected_region, expected_total_length=None) -> AssemblyResult`

**What it does:** The all-in-one checklist — runs every check above plus
content-correctness checks on each named segment, and gives back one
structured result listing every individual check (name, pass/fail, detail
string), instead of one collapsed boolean.

**Why it exists:** This is the direct generalization of my Script 6
six-point manual verification. A UI can render this as an actual checklist
(✓/✗ per row) rather than one opaque pass/fail, which is exactly how I used
to manually eyeball my own assemblies before this was automated.

**Inputs:**
- `full_sequence: str` — the complete assembled sequence to verify.
- `latch_window: LatchWindow`, `graft_spec: GraftSpec`,
  `protected_region: ProtectedRegion` — same as the individual check
  functions above.
- `expected_total_length: int | None` — if given, adds an overall-length
  check; if omitted, that check is skipped entirely.

**Output:** `AssemblyResult` dataclass: `checks: list[AssemblyCheck]` (each
with `name`, `passed`, `detail`), plus an `all_passed` property that's `True`
only if every single check passed.

**The math/logic — the checklist is dynamic, not fixed:**
1. `overall_length` — only if `expected_total_length` was supplied.
2. `protected_region_intact` — always, via `check_protected_region`.
3. `graft_no_overlap` — always, via `check_graft_overlap`.
4. `latch_fit` — always, via `check_latch_fit`.
5. `spacer_intact` — only if `graft_spec.spacer` was supplied; checks the
   spacer's actual content against what was expected, at `spacer_start`.
6. `binder1_intact` — always; checks the binder's actual content at
   `graft_spec.start` against `graft_spec.binder`.
7. `linker_intact` — only if `graft_spec.linker` was supplied.
8. `binder2_intact` — only if `graft_spec.binder2` was supplied.

This design choice is worth calling out explicitly: rather than hardcoding
"exactly 6 checks," the checklist *naturally* produces 6 checks for ECLIPSE's
real v1.0 case (length, protected-region, overlap, fit, spacer, binder1 —
matching Script 6's six lines exactly), because v1.0 has no
linker/second binder, so those two checks simply don't appear. **A full v2.2
tandem call (with spacer, linker, and binder2 all supplied) would produce 8
checks, not 6** — `graft_no_overlap`/`latch_fit` are always present in
addition to whichever of the four content-match checks apply. This is a
judgment call worth being aware of, not a bug: if a stricter 1:1 mapping to
Script 6's exact six lines were wanted even for the tandem case, the checklist
would need to be restructured to *not* always include
`graft_no_overlap`/`latch_fit` as separate rows.

**Worked example (real ECLIPSE v1.0, reproducing Script 6 exactly):**
```python
graft = GraftSpec(binder="LISAAALAAIFAAALAC", start=327,
                  spacer="DA", spacer_start=323)
r = verify_full_assembly(V10, LATCH, graft, SMBIT, expected_total_length=359)
# len(r.checks) == 6
# r.all_passed == True
# names == {"overall_length", "protected_region_intact", "graft_no_overlap",
#           "latch_fit", "spacer_intact", "binder1_intact"}
```
All six pass — exactly reproducing "v1.0: length==359, SmBiT intact,
spacer=='DA', binder1 correct" from Script 6, with the spacer value (`"DA"`)
coming entirely from `graft_spec`, never hardcoded in the check logic itself.

**Worked example (a broken binder1 check, everything else still evaluated):**
```python
graft = GraftSpec(binder="WRONG", start=12)   # doesn't match what's actually there
r = verify_full_assembly(full, window, graft, region)
# protected_region_intact: True
# binder1_intact: False
# r.all_passed: False
```
One failing check doesn't stop the others from running — you get the full
picture, not just the first failure.

**Gotchas:** See the 6-vs-8-checks judgment call above — it's the single
most important design decision in this function to remember if you ever
extend it or compare its output count against Script 6's literal line count
for a tandem design.

---

## `_mutation_position(mutation) -> int`  *(private helper)*

**What it does:** Extracts the numeric position out of a `liability.py`-style
mutation string like `"D4A"`.

**Why it exists:** `filter_safe_variants` needs to compare a mutation's
position against a protected region's bounds, but `liability.py` only hands
back strings, not structured position data — this parses them back out.

**Inputs:** `mutation: str`, format `"{old_residue}{position}{new_residue}"`.

**Output:** `int`.

**The math/logic:** `int("".join(ch for ch in mutation if ch.isdigit()))` —
relies on residue codes always being letters and positions always being
digits, so just stripping non-digit characters reliably recovers the number.

**Worked example:**
```python
_mutation_position("D4A")    # -> 4
_mutation_position("E14A")   # -> 14
```

**Gotchas:** This assumes positions never exceed single/multi-digit
plain integers with no separators — fine for this codebase's format, but
would break if mutation strings ever used a different notation (e.g.
insertion/deletion codes with letters mixed into the position field).

---

## `filter_safe_variants(suggested_variants, protected_region) -> FilteredVariants`

**What it does:** The bridge function — takes a list of `VariantSuggestion`s
(liability.py's output) and a `ProtectedRegion` (assembly.py's concept), and
splits them into ones that are safe to use (no substitution lands inside the
protected region) and ones that aren't, with a reason attached to each
rejection.

**Why it exists:** Keeps `liability.py`'s variant suggester from ever
silently proposing a fix that, while solving the charge problem, accidentally
mutates a position that's structurally off-limits for a completely different
reason (e.g. inside SmBiT). The two concerns — charge liability (soft,
scored) and protected-region integrity (hard constraint) — stay in separate
modules, and this function is the *only* place they ever meet.

**Inputs:**
- `suggested_variants: list[VariantSuggestion]` — typically straight from
  `liability.suggest_variant`, but this function takes the positions **at
  face value**; it does no offsetting itself.
- `protected_region: ProtectedRegion`.

**Output:** `FilteredVariants` dataclass: `accepted: list[VariantSuggestion]`,
`rejected: list[tuple[VariantSuggestion, str]]` — each rejected entry paired
with a human-readable reason like `"substitution at position 12 falls inside
protected region"`.

**The math/logic:** For each variant, check every mutation's position (via
`_mutation_position`) against `[protected_region.start,
protected_region.end]`; if any mutation falls inside, reject with the first
offending position named in the reason string; otherwise accept.

**Worked example (synthetic, one overlapping, one not):**
```python
region = ProtectedRegion(motif="MOTIF", start=10, end=14)
inside = VariantSuggestion(policy="neutralizing", sequence="...", mutations=["D12A"])
outside = VariantSuggestion(policy="neutralizing", sequence="...", mutations=["E20A"])

r = filter_safe_variants([inside, outside], region)
# r.accepted == [outside]
# r.rejected == [(inside, "substitution at position 12 falls inside protected region")]
```

**Worked example (real ECLIPSE case — confirms the real design never
overlaps SmBiT, but only after correctly offsetting):**
```python
v_local = liability.suggest_variant("LISDAELEAIFAEELDC",
                                    preserve_positions=PFLDH_INTERFACE,
                                    policy="neutralizing")
# v_local.mutations == ['D4A', 'E6A', 'E8A', 'E13A', 'E14A', 'D16A']  -- binder-local!

offset = 327 - 1   # binder starts at absolute position 327
abs_mutations = [f"{m[0]}{int(m[1:-1]) + offset}{m[-1]}" for m in v_local.mutations]
v_abs = dataclasses.replace(v_local, mutations=abs_mutations)
# v_abs.mutations == ['D330A', 'E332A', 'E334A', 'E339A', 'E340A', 'D342A']

r = filter_safe_variants([v_abs], SMBIT)   # SmBiT at 312-322
# r.accepted == [v_abs]; r.rejected == []
```
This confirms the real non-overlap was *verified*, not assumed — but notice
this only proves anything because the offset was applied first. A
deliberately-constructed test (`test_filter_safe_variants_uses_absolute_position_not_local_position`
in `test_assembly_general.py`) makes this explicit: a mutation at local
position 12 happens to fall inside a synthetic region of 10-14 if checked
*unoffset*, but after applying a real offset, lands at absolute position 41
— outside the region — and gets accepted. That test exists specifically to
prove the offset math is load-bearing, not just that already-correctly-offset
inputs happen to pass.

**Gotchas — the single most important gotcha in this entire module:**
`filter_safe_variants` does **zero coordinate translation**. If you feed it
`liability.suggest_variant`'s raw output directly, without first applying
`offset = binder_start_in_assembly - 1` to every mutation position, you will
get a **meaningless result** — either a false-pass (binder-local positions
happen to not numerically collide with the protected region, for no real
reason) or a false-reject (they happen to collide coincidentally). This isn't
a defect in the function — it's a deliberate design choice (it has no way to
know where the binder sits in a larger assembly unless told), but it means
**the caller is 100% responsible for getting the coordinate frame right
before calling this function**, every single time.
