# calibration.py

## What this module owns, and what it doesn't

`calibration.py` is pure data — no functions at all, just constants. It's
where every ECLIPSE-specific number that the *general* engine logic (in
`liability.py`) needs as an example/default lives, kept deliberately
separate so it's obvious at a glance which numbers are "mine" versus which
are universal engine behavior. If you ever port this engine to a different
LOCKR system, this is the one file you'd replace wholesale — nothing in
`thermo.py`, `charge.py`, `liability.py`, `models.py`, or `assembly.py`
hardcodes any of these values; they're always passed in as defaults or
arguments.

## General vs. ECLIPSE-specific in this file

**Everything in this file is ECLIPSE-specific, by design.** That's the whole
point of it existing as a separate module from `liability.py` — it's the
example/calibration data, not the general logic.

---

## `PFLDH_INTERFACE = [1, 2, 11, 12, 15]`

**What it is:** The actual positions on ECLIPSE's 17-residue binder that
contact its target, PfLDH (*Plasmodium falciparum* lactate dehydrogenase) —
these came out of the RFdiffusion design process.

**Why it exists:** This is the example value you'd pass as
`liability.py`'s `preserve_positions` argument for the real ECLIPSE binder.
The comment in the file is explicit about this: *"Pass your own
target-interface positions to `scan_liability` — this is not a default."*
Note `liability.py`'s `scan_liability`/`suggest_variant` don't even have a
default value for `preserve_positions` — it's a required argument precisely
so nobody accidentally calls the general scan without thinking about what
their own target interface actually is.

**1-indexed, binder-local** — positions 1, 2, 11, 12, 15 refer to the
17-residue binder sequence itself (`LISDAELEAIFAEELDC` /
`LISAAALAAIFAAALAC`), not any absolute assembly coordinate.

---

## `PENALTY_PER_ACIDIC = 0.8`

**What it is:** The per-residue electrostatic penalty (kcal/mol) that
`liability.py`'s `scan_liability` charges for every flagged acidic residue
outside `preserve_positions`.

**Why it exists, and how it was actually derived:** This isn't a measured
biophysical constant — it's a **two-point calibration**, anchored on
ECLIPSE's own two real binders: the charged original (6 flagged acidic
residues, total documented penalty 4.8 kcal/mol → `4.8/6 = 0.8` per
residue) and the optimized fix (0 flagged residues, 0 penalty). The comment
in the file is explicit that this is "a two-point anchor, not a fit" — with
only two data points (one binder with liabilities, one without), there's no
real curve-fitting happening, just back-solving for the one number that
makes the documented total penalty come out right.

**The TODO that's still open:** `# TODO: recalibrate as more LOCKR systems
are tested.` This value is only as good as the two anchor points it came
from — if TSA/BLI binding data on v1.0/v2.2 (mentioned directly in this
file's module docstring) becomes available, that'd be the first real
opportunity to refine `0.8` into something backed by more than two
sequences.

**Worked example (how this number reproduces the real ECLIPSE result):**
```python
liability.scan_liability("LISDAELEAIFAEELDC", preserve_positions=PFLDH_INTERFACE)
# 6 liabilities flagged, each penalty=0.8 -> penalty_total=4.8
```

---

## `ANCHORS` dict

**What it is:** The two real ECLIPSE binder sequences, with their documented
ground-truth liability counts/positions/penalties/signal status, all in one
place.

**Why it exists:** This is the single source of truth that
`PENALTY_PER_ACIDIC` was derived from, and that the ECLIPSE validation tests
(`test_liability_eclipse.py`, `test_charge_eclipse.py`) check the engine
against. Having it as a structured dict (rather than just inline literals
scattered across test files) makes it possible to look in one place and see
exactly what's "real, documented ECLIPSE data" versus what's
synthetic/made-up test data elsewhere in the test suite.

**Structure:**
```python
ANCHORS = {
    "original": {
        "sequence": "LISDAELEAIFAEELDC",
        "acidic_positions": [4, 6, 8, 13, 14, 16],   # 1-indexed, binder-local
        "n_liabilities": 6,
        "penalty": 4.8,        # kcal/mol
        "signal": False,       # this binder does NOT give clean signal
    },
    "optimized": {
        "sequence": "LISAAALAAIFAAALAC",
        "acidic_positions": [],
        "n_liabilities": 0,
        "penalty": 0.0,
        "signal": True,        # this binder DOES give clean signal
    },
}
```

**Gotcha:** This dict itself is **not currently read by any function in
`liability.py`** — it's documentation/reference data and a source for tests,
not a live input to the engine's calculations. `liability.py`'s actual
calculations re-derive `acidic_positions`/`n_liabilities`/`penalty` from
scratch by scanning the sequence with `scan_liability`, rather than reading
them out of this dict. If you ever update one of the real ECLIPSE sequences
here, double check whether the *engine* needs updating too — this dict won't
automatically stay in sync with what `scan_liability` would actually compute
on a changed sequence, since nothing currently enforces that consistency in
code (only the tests, which compare both independently).

**`signal: False`/`True` fields:** these aren't computed by anything in this
codebase — they're the real documented experimental finding ("did this
binder actually produce a usable signal") used as the ground truth that the
whole liability-scoring system is trying to predict/explain. Nothing in
`liability.py` currently checks whether its `liability_band` classification
("low"/"moderate"/"high") actually lines up with `signal: True/False` for
arbitrary new binders — that correspondence has only been confirmed for
these two specific anchor points.
