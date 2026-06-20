# models.py

## What this module owns, and what it doesn't

`models.py` is the shared vocabulary — every dataclass that the other
modules pass around. It contains **no logic** (no calculations, no
decisions) — just shapes. The one exception that looks like logic is
`SensorParams.luckey_ratio`, which is really just a named property for a
ratio of two of its own fields, kept here so every module that needs
`lucKey/K_CK` computes it identically rather than each module repeating the
division itself.

## General vs. ECLIPSE-specific in this file

- **General:** every dataclass shape — none of them assume a specific
  binder, target, or reporter.
- **ECLIPSE-specific defaults (the base lucCage scaffold's measured/assumed
  values, not universal constants):** `RT_37C = 0.592`, `K_OPEN_DEFAULT =
  1e-3`, `K_CK_DEFAULT = 1e-8`, `LUCKEY_DEFAULT = 500e-9`. A different LOCKR
  team's cage would construct their own `SensorParams(...)` with their own
  measured values — these are just what `DEFAULT_PARAMS` falls back to if you
  don't.

---

## Module-level constants

```python
RT_37C = 0.592          # kcal/mol, gas constant * 310K, used throughout thermo.py
K_OPEN_DEFAULT = 1e-3    # dimensionless — ECLIPSE base scaffold's resting K_open
K_CK_DEFAULT = 1e-8      # ECLIPSE base scaffold's reference lucKey-cage K_CK
LUCKEY_DEFAULT = 500e-9  # M — 500 nM, the lucKey concentration this codebase assumes by default
```

**Why they exist:** These are exactly the values that make `DEFAULT_PARAMS`
(below) reproduce ECLIPSE's real numbers without every test/example having
to construct a custom `SensorParams` from scratch.

**Gotcha — the one open question explicitly flagged in this file's own
comments:** there's a `# TODO: confirm K_CK = lucKey-cage vs lucKey-SmBiT Kd,
inconsistent in my source docs.` This was tracked down and a standardized
definition was settled on (lucKey↔SmBiT epitope on the open cage) — but it's
worth knowing this ambiguity was real and the TODO comment in the code
predates that resolution; if you ever see `K_CK` used inconsistently
somewhere, this is the root cause to check first.

---

## `SensorParams` (frozen dataclass)

**What it is:** The four numbers that pin down a LOCKR sensor's operating
point — `K_open`, `K_CK`, `lucKey`, `RT`. Every `thermo.py` function that
takes a `params` argument is really just asking for one of these.

**Why frozen:** `SensorParams` is meant to represent one fixed, named
configuration (e.g. "ECLIPSE base scaffold") — making it immutable means you
can't accidentally mutate the shared `DEFAULT_PARAMS` instance out from under
every caller that relies on it; if you want a variant, you construct a new
`SensorParams(...)` instead.

**Fields:**
- `K_open: float = K_OPEN_DEFAULT` — dimensionless.
- `K_CK: float = K_CK_DEFAULT` — same units convention as `K_open`
  (dimensionless, treated as an effective equilibrium-constant-like number
  in `thermo.py`'s math, though see the K_CK-definition TODO above).
- `lucKey: float = LUCKEY_DEFAULT` — **M**, not nM. ECLIPSE default 500 nM is
  written as `500e-9`.
- `RT: float = RT_37C` — kcal/mol.

**`luckey_ratio` property:**
```python
self.lucKey / self.K_CK
```
**What it means:** How dominant the lucKey/K_CK side of the system is — a
diagnostic number, used by `thermo.py`'s `_f_open` and `diagnose_regime`.
At ECLIPSE defaults, `500e-9 / 1e-8 = 50.0`.

**Gotcha — the other documented science correction:** this property used to
be conflated with "the achievable fold-change ceiling" in earlier comments.
It isn't — it's a *diagnostic ratio* describing which side of the
equilibrium dominates, not a fold-change you'd ever actually observe. The
real realized ceiling for a given finite `pull` is `thermo.max_fold_change` /
`_saturating_fc`, which can be substantially different from
`luckey_ratio` itself (e.g. `luckey_ratio=50` but `max_fold_change(pull=10)
≈ 11`, not 50). This got renamed/clarified in commit `104bd9c` ("stop calling
max_fold_change a 'ceiling', it's not the dominance ratio") — see
`docs/README.md`'s corrections table for the full story.

---

## `DEFAULT_PARAMS = SensorParams()`

A single shared instance using all the module-level ECLIPSE defaults above.
This is what every `thermo.py` function falls back to if you don't pass your
own `params`. **It's a real, specific object** (not just a type default) —
if you ever needed to verify two functions are using the literal same
defaults, you can check `is DEFAULT_PARAMS`.

---

## `TargetInterface`

**What it is:** A small named container for "these are the binder positions
that touch the target," with an optional human-readable `label`.

**Fields:** `positions: list[int]` (1-indexed), `label: str = ""`.

**Gotcha:** This dataclass exists in the file but isn't actually threaded
through `liability.py`'s functions — `scan_liability`/`suggest_variant` take
a bare `preserve_positions: list[int]` directly, not a `TargetInterface`
object. Worth knowing this type is currently more of a "shape that's
available" than something wired into the live call paths — if you build a UI
on top of this engine, you might construct a `TargetInterface` for display
purposes and then unpack `.positions` when actually calling `liability.py`.

---

## `BinderSequence`

**What it is:** A thin wrapper around a raw sequence string, used wherever
the engine wants to be explicit that "this is a binder, not just any
string."

**Fields:** `sequence: str`, `name: str | None = None`.

**Behavior:**
- `__post_init__` strips whitespace and uppercases the sequence automatically
  — so `BinderSequence("  lisdaeleaifaeeldc  ")` and
  `BinderSequence("LISDAELEAIFAEELDC")` end up identical.
- `__len__` returns the sequence length directly — `len(my_binder)` works.
- `residues()` returns 1-indexed `(position, residue)` pairs via
  `enumerate(self.sequence, start=1)` — this is **the canonical 1-indexing
  convention referenced throughout the entire codebase** (charge.py's
  `helix_breakers`, liability.py's scan positions, assembly.py's `_segment`
  all match this convention).

**Worked example:**
```python
b = BinderSequence("  lisdaeleaifaeeldc  ")
b.sequence       # -> "LISDAELEAIFAEELDC"  (stripped + uppercased)
len(b)           # -> 18
b.residues()[:3] # -> [(1, 'L'), (2, 'I'), (3, 'S')]
```

**Gotcha:** `liability.py`'s functions accept either a plain `str` or a
`BinderSequence` (they wrap plain strings automatically) — but
`BinderSequence`'s own normalization (strip/uppercase) only happens if you
actually go through this wrapper. If a function takes a raw string and skips
wrapping, whitespace/case issues are on you.

---

## `FoldChangeResult`

**What it is:** Every intermediate value from one `thermo.fold_change_detail`
call, bundled for debugging/display.

**Fields:** `target_conc`, `Kd`, `pull`, `theta`, `K_open_eff`, `f_base`,
`f_signal`, `fold_change` — read top to bottom as inputs, then intermediate
steps, then the final answer. See `docs/thermo.md`'s `fold_change_detail`
section for a full worked example.

---

## `ScanResult`

**What it is:** Summary of one `thermo.scan_dose_response` titration sweep.

**Fields:** `label: str` (blank by default — fill in for a UI), `Kd`, `pull`,
`max_fc`, `ec50`, `lod`.

---

## `RegimeResult`

**What it is:** The output of `thermo.diagnose_regime` — the key-limited vs
K_open-limited verdict.

**Fields:** `luckey_dominance_ratio`, `K_open`, `regime` (one of
`"key-limited"`, `"K_open-limited"`, `"mixed"`), `max_fold_change`,
`latch_tuning_helps: bool`, `verdict: str` (a human-readable sentence — this
is the one field in the whole engine meant to be shown verbatim to a user,
not just consumed programmatically).

---

## `ChargeResult`

**What it is:** The output of `charge.analyze_charge`.

**Fields:** `net_charge`, `pH`, `helical_ok: bool`,
`helix_breakers: list[int]` (1-indexed, binder-local positions — never
absolute-assembly coordinates, since `charge.py` never sees an assembly
context at all), `note: str` (empty unless something looked off).

---

## `Liability`

**What it is:** One single flagged residue from a `liability.py` scan.

**Fields:** `position: int` (1-indexed), `residue: str` (the actual amino
acid, e.g. `"D"`), `weight: float` (currently always `1.0` — every flagged
residue counts equally; see the gotcha below), `penalty: float` (kcal/mol
cost this one residue contributes).

**Gotcha:** `weight` exists in the dataclass but `liability.py`'s
`scan_liability` always sets it to a flat `1.0` — there's no actual
per-residue weighting logic implemented yet (e.g. by exposure, by proximity
to the interface). This field is a placeholder for a future refinement, not
something currently doing anything beyond "this residue counts."

---

## `LiabilityReport`

**What it is:** The full output of `liability.scan_liability`.

**Fields:** `binder: BinderSequence`, `liabilities: list[Liability]`,
`preserve_positions: list[int]` (echoed back from the call, for reference),
`net_charge: float`, `penalty_total: float` (kcal/mol), `liability_score:
float` (0-100), `liability_band: str` (`"low"`/`"moderate"`/`"high"`),
`K_CK_estimate: float`.

---

## `VariantSuggestion`

**What it is:** The full output of `liability.suggest_variant` — the
sequence created by mutating away flagged liabilities.

**Fields:** `policy: str`, `sequence: str` (the new, mutated sequence),
`mutations: list[str]` (format `"{old}{position}{new}"`, e.g. `"D4A"`),
`liability_score`/`liability_band`/`K_CK_estimate` (all re-computed on the
**new** sequence, not the original).

**Gotcha — the most important thing about this dataclass:** `mutations`'
position numbers are **binder-local**, relative to whatever sequence string
was passed into `suggest_variant`. This dataclass itself carries no
information about where that binder sits inside a larger assembly. This is
exactly the type `assembly.py`'s `filter_safe_variants` consumes (the one
intentional cross-module dependency, `assembly.py -> models.py`), and the
caller is responsible for translating these positions to absolute assembly
coordinates before checking them against an `assembly.py` `ProtectedRegion`
— see `docs/liability.md`'s and `docs/assembly.md`'s gotchas on this exact
point for the full explanation.

---

## `ProtectedRegion` *(Phase 1.5)*

**What it is:** A motif that must never be altered — ECLIPSE's example is
the 11-residue SmBiT split-luciferase fragment. This is a **hard
constraint**, unlike `preserve_positions` above (a soft tradeoff against
binder affinity) — mutating a `ProtectedRegion` is meant to represent total
loss of function, not a scored penalty.

**Fields:** `motif: str`, `start: int`, `end: int` (1-indexed, inclusive,
**absolute** coordinates in whatever full sequence you're checking against —
not binder-local), `label: str = ""`.

**Gotcha:** Nothing here assumes SmBiT or any specific reporter — the
motif/start/end are entirely caller-supplied. ECLIPSE's SmBiT example data
lives only in `tests/test_assembly_eclipse.py`, never as a default here.

---

## `LatchWindow` *(Phase 1.5)*

**What it is:** The position range a graft is allowed to occupy.

**Fields:** `start: int`, `end: int` (1-indexed, inclusive, absolute),
`expected_length: int | None = None` (currently carried for documentation
purposes — note `check_latch_fit` actually computes `available_length` from
`end - start + 1` directly, it doesn't read `expected_length` off this field;
if you set `expected_length` to something inconsistent with `start`/`end`,
nothing will catch that mismatch automatically).

---

## `GraftSpec` *(Phase 1.5)*

**What it is:** Everything describing one graft into a latch window — the
single-binder case (ECLIPSE v1.0) and the richer tandem case (v2.2) are both
expressed with this same shape, just leaving the optional fields `None` for
the simpler case.

**Fields:** `binder: str`, `start: int` (1-indexed, absolute — where the
binder begins in the full assembly); `spacer: str | None = None`,
`spacer_start: int | None = None` (generalizes ECLIPSE's literal `"DA"` gap
between SmBiT and the binder — not a universal concept, just a named optional
segment); `linker: str | None = None`, `linker_start: int | None = None`;
`binder2: str | None = None`, `binder2_start: int | None = None` (these three
pairs generalize the v2.2 tandem case); `label: str = ""`.

**Gotcha:** None of `spacer`/`linker`/`binder2` are validated for internal
consistency against each other (e.g. nothing checks that `linker_start`
actually comes immediately after `binder`'s end) — `assembly.py`'s functions
trust whatever positions you give them and check content/overlap/fit against
those positions directly. If you get a `*_start` wrong, the checks will
faithfully report on the *wrong* position rather than catching that the
position itself doesn't make geometric sense relative to the other segments.
