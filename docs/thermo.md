# thermo.py

## What this module owns, and what it doesn't

`thermo.py` is the pure math layer. It knows about cages, latches, lucKey,
fold-change, free energy, but it has never seen a sequence, a residue, or a
position. It doesn't know what a "binder" is made of, just that it has a
`Kd` and a `pull`. If you find yourself wanting to pass a sequence into this
file, that's a sign the function belongs in `liability.py` or `assembly.py`
instead.

The three-state model this implements (cage closed / cage open / lucKey bound
to the open cage) comes from Langan et al. 2019 (Nature 572) and
Quijano-Rubio et al. 2021 (Nature 591), that part is real published LOCKR
biology and applies to *any* lucCage-style sensor, not just ECLIPSE. The
specific closed-form fold-change expression (`fold_change`, built from
`theta` + `k_open_eff` + `_f_open`) is **my own algebraic instantiation** of
that framework, written out in the ECLIPSE Thermodynamics doc Section 7, it's
not copy-pasted from either paper, it's me turning their qualitative model
into a formula I can actually compute.

## General vs. ECLIPSE-specific in this file

- **General (true for any LOCKR/lucCage sensor):** every function signature,
  the `theta` / `k_open_eff` / `_f_open` / `fold_change` math, the free-energy
  conversions (`dg_from_kd`, `kd_from_dg`, `kd_from_ddg`), `fit_pull_strength`.
- **ECLIPSE-specific (mine, lives only as a *default*, not a law of physics):**
  `DEFAULT_PARAMS` in `models.py`, `K_open=1e-3`, `K_CK=1e-8`,
  `lucKey=500e-9`, `RT=0.592`. These are *my* base-scaffold's measured/assumed
  values. A different LOCKR team's cage would plug in their own `SensorParams`
  and every function here works unchanged.
- **Heuristic, not biology:** `_K_OPEN_SEARCH_BOUNDS`, `_LUCKEY_SEARCH_BOUNDS`,
  `_NEGLIGIBLE_GAIN=0.02`, `_DOMINANT_MARGIN=2.0` inside `diagnose_regime`.
  The search bounds are generous but bounded so the optimizer stays in
  physically sane territory, the gain/margin numbers are thresholds I
  picked so the classifier behaves sensibly, there's no paper that says
  "2.0x" is special.

**Correction (see git history):** `_f_open`'s denominator originally added
`luckey_ratio` on its own (`k_open / (1 + k_open + luckey_ratio)`). That's
wrong, the open+lucKey-bound state's statistical weight is
`k_open * luckey_ratio` (key binding is *conditional* on the cage already
being open, a sequential equilibrium multiplies, it doesn't add). The same
additive form was also in my original hand-derivation
(`ECLIPSE Thermodynamics Documentation.pdf`, Script 7), so this wasn't just
a coding slip, it was a real error in the underlying physics that both the
doc and the code shared. Fixing it changed real numbers throughout this
file (v1.0's "~11x fold-change at pull=10" is now ~7.41x) and flipped
ECLIPSE's regime classification from "key-limited" to "mixed", see the
worked examples below and `diagnose_regime`'s rewritten section. That PDF
still states the old numbers and needs a separate, manual correction, this
codebase can't fix a PDF.

---

## `theta(target_conc, Kd) -> float`

**What it does:** Tells you what fraction of your binder molecules currently
have target stuck to them, given how much target is around and how tightly
the binder grabs it.

**Why it exists:** Everything downstream (how far the latch gets pulled open)
depends on how *occupied* the binder is, not on the raw target concentration.
This is the standard single-site binding isotherm, converting "concentration
in the tube" into "fraction bound," which is what actually drives signal.

**Inputs:**
- `target_conc: float`, target concentration, in **M** (molar). A sane range
  for a biosensor assay is anywhere from `1e-14` (femtomolar, very dilute) to
  `1e-6` (micromolar, saturating).
- `Kd: float`, binder-target dissociation constant, in **M**. Smaller = tighter
  binding. ECLIPSE v1.0's binder has `Kd = 100e-12` (100 pM); v2.2 has
  `Kd = 42.21e-15` (42.21 fM), over 2000x tighter.

**Output:** `float`, dimensionless, between 0 and 1. 0 = no target bound, 1 =
fully saturated.

**The math:** Standard Langmuir isotherm:
```
theta = [target] / ([target] + Kd)
```
When `target_conc == Kd`, you're at half-occupancy (`theta = 0.5`), that's
the definition of `Kd` itself. When `target_conc >> Kd`, `theta -> 1`.

**Worked example (real ECLIPSE v1.0 numbers):**
```python
theta(1e-6, 100e-12)  # 1 µM target, 100 pM Kd
# -> 0.9999000099990001
```
At 1 µM target with a 100 pM `Kd`, the binder is essentially saturated
(`theta ≈ 0.9999`), that's a 10,000x excess of target over `Kd`, so this is
deep in the saturating regime, which is exactly why the ECLIPSE validation
tests use `1e-6` as the "saturating" target concentration.

**Gotchas:**
- Both `target_conc` and `Kd` must be in the **same concentration units** (this
  codebase uses molar, M, everywhere, not nM!). Passing `500` instead of
  `500e-9` for a 500 nM lucKey concentration is the single easiest unit bug to
  make in this whole engine.
- `theta` doesn't know or care about lucKey or the cage at all, it's pure
  target-binder chemistry, computed completely independently of the
  LOCKR-specific stuff that happens downstream.

---

## `k_open_eff(K_open, pull, theta) -> float`

**What it does:** Computes how much *more* open the latch effectively is once
the target has pulled on it, compared to its resting (no-target) openness.

**Why it exists: this is the actual allosteric coupling step**, the whole
point of a LOCKR sensor is that target binding mechanically biases the latch
toward the open state. `pull` is how strong that allosteric coupling is (a
property of your *cage/latch design*, not of the target or the binder), and
this function turns "latch is `theta` fraction pulled" into an effective
open-equilibrium constant.

**Inputs:**
- `K_open: float`, dimensionless equilibrium constant for the latch's
  spontaneous (unpulled) closed↔open equilibrium. ECLIPSE default: `1e-3`, meaning at rest, the cage is open only about 1 in 1000 of the time (heavily
  biased closed, which is the point of a locked cage).
- `pull: float`, dimensionless allosteric coupling strength. Literature LOCKR
  designs run roughly 10-20x; this codebase treats `pull=10` as a generic
  default, not an ECLIPSE-measured number (ECLIPSE's actual fitted/assumed
  pull isn't separately documented, it's swept across `10, 20, 30` in tests
  as a "what if" parameter).
- `theta: float`, the binder-occupancy fraction from `theta()` above, 0 to 1.

**Output:** `float`, an effective `K_open`-like equilibrium constant, same
units (dimensionless) as `K_open`, just bigger when `theta` is bigger.

**The math:**
```
K_open_eff = K_open * (1 + pull * theta)
```
At `theta=0` (no target), `K_open_eff == K_open`, no change. As `theta -> 1`
(fully bound), `K_open_eff -> K_open * (1 + pull)`, the latch's open
equilibrium gets boosted by a factor of `(1+pull)`.

**Worked example:**
```python
k_open_eff(1e-3, 10, 1.0)   # K_open=1e-3, pull=10, fully saturated
# -> 0.011
```
Saturating target with `pull=10` makes the latch's effective open-constant
11x its resting value (`0.001 -> 0.011`).

**Gotchas:** `pull` is a knob you're choosing/fitting, not something you look
up, if you ever have real titration data, you back out `pull` with
`fit_pull_strength` instead of guessing it.

---

## `_f_open(k_open, params) -> float`  *(private helper)*

**What it does:** Converts an open-equilibrium constant into an actual
*fraction of cages* sitting in the open, signal-competent state, accounting
for the fact that lucKey is competing to occupy that open state too.

**Why it exists:** `K_open` alone tells you closed-vs-open, but once a cage
opens, lucKey grabs it (governed by `K_CK`) and that's what actually produces
luminescence. This function is the three-state partition function collapsed
into one line, it's where `K_CK` and `lucKey` (via `params.luckey_ratio`)
enter the model.

**Inputs:**
- `k_open: float`, an (effective) open-equilibrium constant, e.g. from
  `k_open_eff()` or just `params.K_open` for baseline.
- `params: SensorParams`, needed because it carries `luckey_ratio =
  lucKey/K_CK` (the dominance ratio).

**Output:** `float`, dimensionless, the fraction of all cage molecules in the
open+lucKey-bound state. Always small (this is a rare, productive state).

**The math:**
```
weight_signal = k_open * luckey_ratio
f_open = weight_signal / (1 + k_open + weight_signal)
```
This is a proper three-state partition function: closed (weight 1),
open-but-empty (weight `k_open`), open-and-lucKey-bound (weight
`k_open * luckey_ratio`, since key binding is *conditional* on the cage
already being open, a sequential equilibrium, `[OK]/[C] = ([O]/[C]) *
([Key]/K_CK) = k_open * luckey_ratio` by plain mass action). The numerator
is the open+key-bound population specifically, not "open" in general, an
open-but-empty cage has an exposed SmBiT with no LgBiT partner, it doesn't
glow.

**Worked example:** baseline (no target) at ECLIPSE defaults:
```python
_f_open(1e-3, DEFAULT_PARAMS)   # K_open itself, no pull applied
# -> 0.04757373929590866
```
About 4.8% of cages are in the productive state at rest, that's the
dark/background signal level.

**Gotcha (this was a real bug, not hypothetical):** the earlier version of
this function added `luckey_ratio` on its own instead of multiplying it by
`k_open`. That formula is monotonically *decreasing* toward zero as lucKey
concentration grows without bound, i.e. it implied adding more of a
bimolecular binding partner eventually suppresses its own complex, which is
backwards. A bimolecular association can only saturate (plateau), never
reverse, as one partner's concentration increases. If you ever see
`luckey_ratio` show up anywhere *not* multiplied by an opening-equilibrium
term, that's the same bug back again.

**Gotchas:** This is named with a leading underscore because it's an internal
building block, `fold_change` and `f_base` are the public functions you
actually call; don't reach for `_f_open` directly unless you're debugging.

---

## `fold_change(target_conc, Kd, pull, params=DEFAULT_PARAMS) -> float`

**What it does:** The headline number, given a target concentration, how
many times brighter is your sensor than its own dark baseline?

**Why it exists:** This is *the* answer a biosensor designer actually wants.
Everything above this line in the file is plumbing to get here.

**Inputs:**
- `target_conc: float`, M, concentration of analyte you're testing.
- `Kd: float`, M, binder's affinity for that analyte.
- `pull: float`, dimensionless, your cage's allosteric coupling strength.
- `params: SensorParams`, `K_open`, `K_CK`, `lucKey`, `RT` for your specific
  cage/key system; defaults to the ECLIPSE base scaffold.

**Output:** `float`, dimensionless fold-change. `1.0` = no signal above
baseline; bigger = brighter.

**The math:** chains everything above, ```
theta      = target_conc / (target_conc + Kd)
K_open_eff = K_open * (1 + pull*theta)
fold_change = f_open(K_open_eff) / f_open(K_open)     # signal / baseline
```

**Worked example (ECLIPSE v1.0, saturating target):**
```python
fold_change(1e-6, 100e-12, 10)
# -> 7.405718615614029
```
Corrected from the pre-`_f_open`-fix "~11x", see the file-level correction
note above. This is now the actual number the corrected model predicts,
not (yet) checked against real wet-lab titration data.

**Gotchas:** Note `Kd` only enters through `theta`, at saturating target
concentration, `theta -> 1` regardless of how tight `Kd` is, so a *tighter*
binder (lower `Kd`, like v2.2's 42.21 fM) gives you the *same* max fold-change
as a looser one (v1.0's 100 pM) at saturation. What a tighter `Kd` actually
buys you is a lower **EC50**, you saturate at a much lower target
concentration, i.e. better sensitivity/limit-of-detection, not a bigger
ceiling. That distinction is exactly why `max_fold_change` ignores `Kd`.

---

## `fold_change_detail(target_conc, Kd, pull, params=DEFAULT_PARAMS) -> FoldChangeResult`

**What it does:** Same calculation as `fold_change`, but instead of handing
back just the final number, it hands back every intermediate value along the
way.

**Why it exists:** For a debug view or a result card in a UI, you want to be
able to show *why* a fold-change came out the way it did, was it
occupancy-limited (`theta` low) or cage-limited (`f_signal`/`f_base` both
tiny)? A bare float can't tell you that.

**Inputs:** identical to `fold_change`.

**Output:** `FoldChangeResult` dataclass with fields `target_conc, Kd, pull,
theta, K_open_eff, f_base, f_signal, fold_change`, read it top to bottom as
"here's what I was given, here's what I computed at each step, here's the
final answer."

**The math:** exactly `fold_change`'s pipeline, just not throwing away the
intermediates.

**Worked example:**
```python
fold_change_detail(1e-6, 100e-12, 10)
# FoldChangeResult(
#   target_conc=1e-06, Kd=1e-10, pull=10,
#   theta=0.9999000099990001,
#   K_open_eff=0.01099900009999,
#   f_base=0.04757373929590866,
#   f_signal=0.3523177267180794,
#   fold_change=7.405718615614029,
# )
```

**Gotchas:** none beyond `fold_change`'s, this is purely a "show your work"
wrapper.

---

## `f_base(params=DEFAULT_PARAMS) -> float`

**What it does:** The baseline (dark, no-target) open fraction, your
sensor's noise floor.

**Why it exists:** You need this as the denominator for every fold-change
calculation, and it's also useful on its own to ask "how dark is dark" for a
given cage design.

**Inputs:** `params: SensorParams`, just `K_open`/`K_CK`/`lucKey`/`RT`.

**Output:** `float`, dimensionless, same meaning as `_f_open`'s output.

**The math:** `_f_open(K_open, params)`, i.e. `f_open` evaluated with no
pull applied at all (`k_open` unmodified).

**Worked example (ECLIPSE defaults):**
```python
f_base()
# -> 0.04757373929590866
```

**Gotchas:** This is literally `_f_open(params.K_open, params)` with zero
allosteric boost, don't confuse it with `f_signal` in `FoldChangeResult`,
which *is* boosted by target+pull.

---

## `_saturating_fc(pull, params) -> float`  *(private helper)*

**What it does:** Fold-change at `theta=1`, i.e., what your fold-change
would be if you dumped in an essentially infinite amount of target.

**Why it exists:** This is the realistic "ceiling" for a *finite* `pull`, distinct from two other things that sound similar but aren't:
1. `params.luckey_ratio` (the lucKey/K_CK dominance ratio), that's a
   diagnostic number describing the *key* side of the system, not a fold-change.
2. The true `pull -> infinity` asymptote, `(1+K_open+luckey_ratio)/K_open`, this codebase deliberately doesn't compute that anywhere, because no real
   cage has infinite pull; `_saturating_fc` answers the practical version of
   the question instead.

This distinction is exactly the "ceiling vs. dominance ratio" naming
confusion that got fixed in commit `104bd9c`, the function used to imply
`luckey_ratio` itself was an achievable fold-change ceiling, which it isn't;
it's a ratio describing which side of the equilibrium dominates.

**Inputs:** `pull: float`, `params: SensorParams`.

**Output:** `float`, dimensionless fold-change at full target saturation, for
*this specific* `pull` value.

**The math:** `k_open_eff(K_open, pull, theta=1.0)` fed into `_f_open`,
divided by the baseline `_f_open(K_open, params)`.

**Worked example:** same as `max_fold_change` below (they're the same call).

**Gotchas:** private helper, call `max_fold_change` instead unless you're
inside this file.

---

## `max_fold_change(Kd, pull, params=DEFAULT_PARAMS) -> float`

**What it does:** The realistic best-case fold-change for your sensor at a
given `pull`, assuming target is in vast excess.

**Why it exists:** Lets you ask "is this design even capable of a big signal,
regardless of how much target I throw at it", useful for comparing cage
designs (different `pull`/`K_open`/`K_CK`/`lucKey`) independent of any one
binder's affinity.

**Inputs:**
- `Kd: float`, **accepted but unused on purpose.** It's kept in the
  signature purely so call sites that loop over `fold_change`-like functions
  don't need a special case for this one. The max is set by the *cage*
  (`K_open`, `K_CK`, `lucKey`, `pull`), not by how tight the binder is, see
  the `fold_change` gotcha above for why.
- `pull: float`, `params: SensorParams`, same as elsewhere.

**Output:** `float`, dimensionless fold-change ceiling for this `pull`.

**The math:** `_saturating_fc(pull, params)`, identical to `_saturating_fc`,
just the public name for it.

**Worked example (real ECLIPSE numbers, both binders):**
```python
max_fold_change(100e-12, 10)    # v1.0's Kd, pull=10
# -> 7.4061499039077505
max_fold_change(42.21e-15, 10)  # v2.2's Kd, SAME pull
# -> 7.4061499039077505   (identical!)
```
v1.0 and v2.2 hit the *same* ceiling at the same `pull`, even though v2.2's
`Kd` is 2369x tighter, proving the point above: tighter `Kd` buys
sensitivity (lower EC50), not a higher ceiling. At `pull=20` both go to
`~10.66`, at `pull=30` both go to `~12.62`. Unlike the pre-fix formula, this
isn't `1 + pull` at ECLIPSE's actual `K_open`/`lucKey`/`K_CK`, ECLIPSE turns
out to be a **mixed** regime, not purely "key-limited", `1 + pull` is only
the ceiling in the limit where either `K_open` or `lucKey/K_CK` is pushed to
its own local optimum (see `diagnose_regime` below).

**Gotchas:** Don't be tempted to "fix" the unused `Kd` parameter, it's
intentional, and removing it would break every call site that treats
`fold_change` and `max_fold_change` as interchangeable in a loop.

---

## `scan_dose_response(Kd, pull, params=DEFAULT_PARAMS, n=500) -> ScanResult`

**What it does:** Simulates a full titration curve (target concentration from
femtomolar to micromolar) and reads off the practically useful summary
numbers: max fold-change, EC50 (the concentration giving half-max signal),
and an estimated limit of detection.

**Why it exists:** A single fold-change number at one concentration doesn't
tell you the *shape* of your dose-response curve. This answers "where's the
useful dynamic range of this sensor," which is what you'd actually plot
against real assay data.

**Inputs:**
- `Kd: float`, M; `pull: float`; `params: SensorParams`, same as `fold_change`.
- `n: int`, how many points to simulate across the concentration sweep
  (default 500, just a resolution knob, not biologically meaningful).

**Output:** `ScanResult` dataclass: `label` (blank, fill in yourself for a
UI), `Kd`, `pull`, `max_fc`, `ec50`, `lod`.

**The math:** Sweeps `target_conc` logarithmically from `1e-14` to `1e-5` M
(`np.logspace(-14, -5, n)`), computes `fold_change` at each point, takes the
`.max()` as `max_fc`. `ec50` is found by locating the concentration whose
fold-change is closest to the half-max point, `(max_fc + 1) / 2` (the `+1`
accounts for the baseline being fold-change-of-1, not 0). `lod` is just
`ec50 * 0.1`, a simple rule-of-thumb "10% of EC50" estimate, not a
statistically derived detection limit (e.g. not based on noise/3σ).

**Worked example (ECLIPSE v1.0):**
```python
scan_dose_response(100e-12, 10)
# ScanResult(label='', Kd=1e-10, pull=10,
#            max_fc=7.4061067724649705,
#            ec50=6.662654524581163e-11,
#            lod=6.662654524581163e-12)
```

**Gotcha, changed by the `_f_open` fix:** EC50 used to come out almost
exactly at `Kd` under the old (buggy) additive formula, that was noted here
as a sanity check. It no longer does (6.66e-11 vs `Kd=1e-10`, about 0.67x).
That's expected, not a new bug: half-occupancy (`theta=0.5`) still happens
exactly at `target_conc=Kd` by definition, but fold-change is now a
genuinely nonlinear (saturating) function of `theta` through the corrected
`_f_open`, so "half-max fold-change" and "half-occupancy" no longer land at
the same target concentration. Don't resurrect the old "EC50 ≈ Kd" check,
it was a symptom of the bug, not a real property of this model.

**Gotchas:** `lod = ec50 * 0.1` is a placeholder heuristic, explicitly called
out as such, if you ever have real assay noise data, you'd want to replace
this with an actual statistical LOD (e.g. blank + 3×SD), not this 10% rule.
This is exactly the kind of thing to flag rather than silently treat as
validated science.

---

## `dg_open_cost(params=DEFAULT_PARAMS) -> float`

**What it does:** The free-energy cost of cracking the latch open, in
kcal/mol, how much energetic "work" the cage's closed state is doing to stay
closed.

**Why it exists:** Translates the abstract `K_open` equilibrium constant into
a physically interpretable energy number, useful for comparing against
actual measured ΔG values from biophysical assays (ITC, etc.), and it's the
quantity referenced directly in the ECLIPSE Thermodynamics doc.

**Inputs:** `params: SensorParams` (`K_open`, `RT`).

**Output:** `float`, kcal/mol. Always positive for `K_open < 1` (closed state
favored), that's the energy barrier holding the cage shut.

**The math:** Boltzmann relation, `ΔG = -RT·ln(K_open)`.

**Worked example (ECLIPSE defaults):**
```python
dg_open_cost()
# -> 4.089391125157425   # ≈ 4.09 kcal/mol
```

**Gotchas:** Sign convention, this function returns the cost as a
**positive** number (opening is unfavorable), which is the opposite sign
convention from `dg_from_kd` below, where a **tight** (favorable) `Kd` gives a
**negative** ΔG. Don't assume all ΔG values in this codebase share one sign
convention; check which function produced the number.

---

## `dg_luckey(params=DEFAULT_PARAMS) -> float`

**What it does:** The free energy released by lucKey binding once the latch
has cracked open, i.e., how much energetic "reward" drives the cage into
the signal-producing state.

**Why it exists:** Pairs with `dg_open_cost` to give you both halves of the
energetic picture: cost to open, reward once open. Together they describe
why the system has any signal at all.

**Inputs:** `params: SensorParams` (`luckey_ratio`, `RT`).

**Output:** `float`, kcal/mol. Negative = favorable (energy released).

**The math:** `ΔG = -RT·ln(luckey_ratio)`, where `luckey_ratio =
lucKey/K_CK` is being treated here as an effective binding "constant" for
lucKey grabbing the open cage.

**Worked example (ECLIPSE defaults):**
```python
dg_luckey()
# -> -2.3159176192134625   # ≈ -2.32 kcal/mol
```

**Gotchas:** This uses `luckey_ratio` (a *concentration-dependent* ratio,
`lucKey/K_CK`), not `K_CK` alone, so this number changes if you change how
much lucKey you add to the assay, even though `K_CK` (a property of the
molecules themselves) hasn't changed. Don't mistake this for a fixed
biophysical constant of the lucKey-cage interaction.

---

## `dg_from_kd(Kd, RT=DEFAULT_PARAMS.RT) -> float`

**What it does:** Converts any dissociation constant into a binding free
energy.

**Why it exists:** General-purpose unit conversion used anywhere you have a
`Kd` (binder-target, lucKey-cage, whatever) and want the energetic
equivalent, e.g. to reason about how many kcal/mol of improvement a design
change is worth.

**Inputs:** `Kd: float`, M; `RT: float`, kcal/mol (defaults to 37°C's
`0.592`).

**Output:** `float`, kcal/mol.

**The math:** `ΔG = RT·ln(Kd)`, note **no leading minus sign** here, unlike
`dg_open_cost`. This is intentional: by convention, `dG = RT ln(Kd)` directly
(since `Kd = exp(dG/RT)`), so a tight `Kd` (small, e.g. `1e-12`) gives a large
negative ΔG (favorable), and a weak `Kd` (large, e.g. `1e-6`) gives a less
negative (or positive) ΔG.

**Worked example (ECLIPSE v1.0):**
```python
dg_from_kd(100e-12)
# -> -13.63130375052475   # kcal/mol
```

**Gotchas:** See the sign-convention note under `dg_open_cost`, this
function and that one are *not* both wrapping the same sign convention, even
though they're both "ΔG from an equilibrium constant." `dg_from_kd` follows
the standard binding-energy convention (`RT ln(Kd)`, no extra minus);
`dg_open_cost`/`dg_luckey` flip the sign to express things as a "cost"/
"reward" framing for readability in the regime diagnostic.

---

## `kd_from_dg(dG, RT=DEFAULT_PARAMS.RT) -> float`

**What it does:** The inverse of `dg_from_kd`, given a binding energy, what
`Kd` does that correspond to?

**Why it exists:** Lets you go the other direction, if a computational
design tool (Rosetta, AF, etc.) reports a predicted ΔΔG, you need this to
turn it back into a `Kd` you can actually plug into `fold_change`.

**Inputs:** `dG: float`, kcal/mol; `RT: float`, kcal/mol.

**Output:** `float`, `Kd` in M.

**The math:** `Kd = exp(dG/RT)`, exact inverse of `dg_from_kd`.

**Worked example:** round-tripping v1.0's `Kd`:
```python
kd_from_dg(dg_from_kd(100e-12))
# -> 1.0000000000000002e-10   # back to 100e-12, modulo float noise
```

**Gotchas:** Make sure the `dG` you feed in actually used the
`RT ln(Kd)` convention (no extra minus sign), if you grabbed a ΔG from
`dg_open_cost`'s convention by mistake, you'll get the reciprocal `Kd`
instead of the one you wanted.

---

## `kd_from_ddg(kd_ref, ddg, RT=DEFAULT_PARAMS.RT) -> float`

**What it does:** Rescales a reference `Kd` by a binding-energy *change*
(ΔΔG), answers "if this design tweak is worth `ddg` kcal/mol, what's my new
`Kd`?"

**Why it exists:** This is exactly the question a design iteration asks:
"RFdiffusion/Rosetta says this mutant should be `ddg` kcal/mol better/worse, what `Kd` do I actually expect?" It's how v1.0's measured `Kd` (100 pM) turns
into v2.2's predicted `Kd` (42.21 fM) in this codebase.

**Inputs:**
- `kd_ref: float`, M, your starting point's `Kd`.
- `ddg: float`, kcal/mol, the energy change. **Negative `ddg` means
  tighter** (more favorable) binding, that's the sign convention to
  remember.
- `RT: float`, kcal/mol.

**Output:** `float`, the rescaled `Kd` in M.

**The math:** `Kd_new = Kd_ref * exp(ddg/RT)`, a negative `ddg` shrinks
`exp(ddg/RT)` below 1, so `Kd_new < Kd_ref` (tighter binding, smaller number).

**Worked example (real ECLIPSE v1.0 → v2.2 improvement):**
```python
kd_from_ddg(100e-12, -4.6)
# -> 4.22099170534448e-14   # i.e. 42.21 fM, matches v2.2's documented Kd

100e-12 / kd_from_ddg(100e-12, -4.6)
# -> 2369.1115022420754     # the documented "~2369x tighter" improvement
```
This is the actual computation behind the "v2.2 is 2369x tighter than v1.0"
claim in the ECLIPSE docs, a -4.6 kcal/mol design improvement.

**Gotchas:** Sign trap, it's easy to flip the sign of `ddg` by accident.
Remember: **negative = better** (tighter `Kd`). If your fold-improvement
comes out *less than 1* (i.e. `Kd` got worse, not better) when you expected
an improvement, you probably have the sign backwards.

---

## `diagnose_regime(params=DEFAULT_PARAMS, pull=10.0) -> RegimeResult`

**What it does:** Answers the single most actionable design question this
engine can ask: *if I spent effort improving the latch (`K_open`), would it
actually move my fold-change at all?* Or is something else (the lucKey/K_CK
side) the real bottleneck?

**Why it exists:** This is the practical payoff of the whole thermodynamic
model, it tells you where to spend your limited engineering effort, and
which direction to push it.

**Rewritten along with the `_f_open` fix.** The old version assumed "more
open" and "more lucKey" always help, and only ever recommended raising
things. That assumption came from the same additive bug in `_f_open`,
under the corrected model, `max_fold_change` is **not monotonic** in either
`K_open` or `lucKey/K_CK`, past a point, more lucKey raises the dark
background as fast as the signal, so there's a real interior optimum, not a
one-way ramp. Sometimes the right move is to *lower* lucKey or *tighten*
the latch, not raise either one, and the old version could never say that.

**Inputs:**
- `params: SensorParams`, the system you're diagnosing.
- `pull: float`, default `10.0`, a generic "what fold-change would this
  produce" probe value, not a measured ECLIPSE number.

**Output:** `RegimeResult` dataclass: `luckey_dominance_ratio`, `K_open`,
`regime` (one of `"key-limited"`, `"K_open-limited"`, `"mixed"`),
`max_fold_change`, `latch_tuning_helps: bool`, a human-readable `verdict`
string, and now `recommendations: list[str]` (direction-aware, generated
per call, not a static lookup table).

**The math/logic:** For each axis (`K_open`, `lucKey`), find that axis's
*actual local optimum* via a bounded search (`scipy.optimize.minimize_scalar`
over `log10(param)`, see `_best_along`), holding everything else fixed, and
compute the relative gain available by moving there:
1. `gain_k_open` = best achievable relative improvement in `max_fold_change`
   by moving `K_open` alone to its optimum (search range `1e-6` to `0.5`).
2. `gain_luckey` = same, moving `lucKey` alone (search range `1 pM` to
   `1 mM`).
3. If both gains are below `_NEGLIGIBLE_GAIN` (`0.02`), neither axis has
   real headroom left → **mixed**, `latch_tuning_helps=False`, verdict
   points at `pull` (the allosteric coupling / cage-latch geometry) as the
   only remaining lever.
4. If one gain is at least `_DOMINANT_MARGIN` (`2.0x`) the other, that axis
   is named the limiting one (**key-limited** or **K_open-limited**), with
   a direction-aware recommendation (raise or lower, whichever the optimum
   search found).
5. Otherwise **mixed**: both axes have comparable headroom, recommendations
   list both moves.

Log-space search because `K_open`/`lucKey` span many orders of magnitude and
the curve is unimodal in the log of either one, a fixed linear step (or the
old version's fixed 30x multiplier) can overshoot a nearby peak and come
back down on the other side, wrongly reporting "no headroom" when a smaller
move would have helped.

**Worked example, real ECLIPSE default (500 nM lucKey):**
```python
diagnose_regime(pull=10)
# RegimeResult(
#   luckey_dominance_ratio=50.0, K_open=0.001, regime='mixed',
#   max_fold_change=7.4061499039077505, latch_tuning_helps=True,
#   verdict="Mixed: lucKey/K_CK = 50.0, both axes have comparable headroom, "
#           "lowering K_open from 0.001 toward 1e-06 (engineer the latch "
#           "for a stronger closed state) would take fold-change to about "
#           "11.0x (+48%), and separately lowering lucKey from 500 nM "
#           "toward 0.001 nM (or the equivalent move in K_CK) would take "
#           "fold-change to about 10.9x (+47%).",
#   recommendations=[
#     "Lowering K_open from 0.001 toward 1e-06 ... would take fold-change to about 11.0x.",
#     "Lowering lucKey from 500 nM toward 0.001 nM ... would take fold-change to about 10.9x."])
```
This inverts the old "key-limited, latch tuning won't help" conclusion:
ECLIPSE's actual default is **mixed**, tightening *either* the latch or the
key affinity buys a comparable ~47-48% improvement, and both moves are
*toward tighter*, not toward more lucKey, the old advice's direction was
backwards as well as its regime label.

**Worked example, same cage, hypothetically dropping lucKey to 10 nM:**
```python
diagnose_regime(SensorParams(lucKey=10e-9), pull=10)
# RegimeResult(
#   luckey_dominance_ratio=1.0, regime='mixed',
#   max_fold_change=10.784735812133071, latch_tuning_helps=False,
#   verdict="Near-optimal: lucKey/K_CK = 1.0 and K_open = 0.001 are both "
#           "already close to their best achievable values at pull=10 "
#           "(fold-change 10.8x, ceiling ~11.0x from either alone). The "
#           "remaining lever is pull itself, redesigning the cage-latch "
#           "allosteric coupling, not lucKey/K_CK or K_open.")
```
At 10 nM lucKey this cage is already close to the best either axis alone
can buy for `pull=10` (10.8x vs an ~11.0x ceiling), the old version called
this "K_open-limited, latch tuning helps", which is no longer true, there's
almost nothing left on that axis either.

**Gotchas:**
- `_NEGLIGIBLE_GAIN` (`0.02`) and `_DOMINANT_MARGIN` (`2.0x`) are heuristics
  I picked, not biology, see the file-level note above. Don't quote them as
  if they came from a paper.
- `pull` here defaults to `10.0` as a generic round-number probe, distinct
  from any specific measured ECLIPSE pull value (which isn't separately
  documented in this codebase).
- The search bounds (`K_open` in `[1e-6, 0.5]`, `lucKey` in `[1pM, 1mM]`)
  are generous but real, if your actual design needs to search outside
  them, widen `_K_OPEN_SEARCH_BOUNDS`/`_LUCKEY_SEARCH_BOUNDS`, don't just
  eyeball a fixed-factor probe back in.

---

## `fit_pull_strength(target_conc, fc_measured, Kd, params=DEFAULT_PARAMS) -> (pull, pull_stderr)`

**What it does:** The one function in this file that runs in the opposite
direction from everything else, instead of predicting fold-change from a
known `pull`, it takes *real measured* titration data (a list of target
concentrations and the fold-changes you actually observed at each) and backs
out what `pull` value explains your data.

**Why it exists:** `pull` isn't something you can look up or calculate from
first principles, it's an emergent property of your specific cage/latch
geometry that you can only learn from a real dose-response experiment. This
is the bridge from "I have assay data" to "I have a `pull` number I can now
plug into `fold_change`/`max_fold_change`/`diagnose_regime` for predictions."

**Inputs:**
- `target_conc`: array-like of target concentrations, M, your titration
  series.
- `fc_measured`: array-like, same length, the fold-changes you measured at
  each concentration.
- `Kd: float`, M, your binder's already-known affinity (assumed fixed, not
  fit).
- `params: SensorParams`, your cage's other known constants.

**Output:** A tuple `(pull, pull_stderr)`, the best-fit `pull` and its
standard error from the curve fit.

**The math/logic:** Wraps `scipy.optimize.curve_fit` around a little model
function that just calls `fold_change` at each concentration for a candidate
`pull`, then finds the `pull` that minimizes squared error against your
measured data, bounded to `[0, 100]` with an initial guess of `10.0`.

**Worked example (synthetic, this function is validated against made-up
data in `test_thermo_general.py`, not real ECLIPSE titration data, since I
don't have real measured fold-change-vs-concentration numbers for ECLIPSE in
this codebase yet):**
```python
true_pull = 23.0
conc = np.logspace(-13, -6, 8)
fc = [fold_change(c, 5e-10, true_pull, params) for c in conc]
pull, stderr = fit_pull_strength(conc, fc, 5e-10, params)
# pull ≈ 23.0  (recovers the true value when data is noise-free)
```
