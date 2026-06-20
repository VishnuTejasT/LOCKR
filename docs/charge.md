# charge.py

## What this module owns, and what it doesn't

`charge.py` is plain amino-acid chemistry: given a sequence string and a pH,
what's its net electric charge, and does it look like it'll fold into a
helix. It doesn't know anything about LOCKR, cages, lucKey, K_CK, or
positions to preserve — it would work identically if you fed it a sequence
from a completely unrelated protein. `liability.py` is the module that
*interprets* charge in a LOCKR-specific way (acidic residues as "liabilities"
that weaken K_CK); `charge.py` itself stays generic.

## General vs. ECLIPSE-specific in this file

Everything here is general — the pKa table is a standard textbook set (the
EMBOSS side-chain pKa values), the Chou-Fasman helix propensity table is a
published 1974 statistical scale, and the default `pH=7.4` is just
physiological pH, not an ECLIPSE assumption. There is **nothing
ECLIPSE-specific in this file at all** — even the two binder sequences used
in its ECLIPSE test file are just used as real-world inputs to otherwise
fully general functions.

---

## `_protonated_fraction(pKa, pH) -> float`  *(private helper)*

**What it does:** What fraction of a given ionizable group (a side chain or a
terminus) is currently holding onto its extra proton, at a given pH.

**Why it exists:** This is the Henderson-Hasselbalch equation in fractional
(not log) form — every other charge calculation in this file is built by
summing these fractions across all the ionizable groups in a sequence.

**Inputs:**
- `pKa: float` — the group's pKa (dimensionless; not really a "pH unit" but
  same numeric scale as pH).
- `pH: float` — the solution pH.

**Output:** `float`, between 0 and 1. 1 = fully protonated, 0 = fully
deprotonated.

**The math:** `f = 1 / (1 + 10^(pH - pKa))` — the standard
Henderson-Hasselbalch relation rearranged to give a fraction directly instead
of a log-ratio. When `pH == pKa`, you get exactly `0.5` (half-protonated, by
definition of pKa).

**Worked example:** Aspartate's side-chain pKa is 3.65; at physiological
pH 7.4 it's way above its pKa, so it should be almost fully deprotonated
(i.e. `_protonated_fraction` should be near 0):
```python
_protonated_fraction(3.65, 7.4)
# -> 0.00017782794100389228   # ~0.018% protonated — essentially fully ionized
```

**Gotchas:** This returns the *protonated* fraction, not the *charged*
fraction — for an acidic residue, "protonated" means *neutral* (COOH) and
"deprotonated" means *negatively charged* (COO⁻). `net_charge` below
subtracts `(1 - protonated_fraction)` for acidic residues for exactly this
reason — don't reuse this function's raw output as a charge contribution
without checking which direction it needs flipping for acid vs. base.

---

## `net_charge(sequence, pH=7.4) -> float`

**What it does:** The overall electric charge of a peptide at a given pH,
accounting for every acidic/basic side chain plus both termini.

**Why it exists:** Net charge is a real, measurable, biologically meaningful
property — it affects solubility, aggregation, and (in the LOCKR-specific
interpretation that `liability.py` builds on top of this) how strongly a
binder's acidic residues might interfere with the cage's own electrostatics.
This function itself is just the chemistry; `liability.py` is what turns it
into a design signal.

**Inputs:**
- `sequence: str` — single-letter amino acid codes; case-insensitive, gets
  uppercased and stripped internally.
- `pH: float`, default `7.4` (physiological pH).

**Output:** `float`, net charge in elementary-charge units (e.g. `-6.13`
means "6.13 units of negative charge," not an integer — because each
ionizable group only *partially* contributes its charge at a given pH,
unless that pH is far from its pKa).

**The math:**
1. Start with the two termini: N-terminus contributes `+_protonated_fraction(8.6, pH)`
   (it's a base — protonated means positively charged, so use the fraction
   directly); C-terminus contributes `-(1 - _protonated_fraction(3.6, pH))`
   (it's an acid — deprotonated means negatively charged).
2. For each basic side chain (`H`, `K`, `R`): add `_protonated_fraction(pKa, pH)`
   (protonated = positively charged, same direction as the N-terminus).
3. For each acidic side chain (`D`, `E`, `Y`, `C`): subtract
   `1 - _protonated_fraction(pKa, pH)` (deprotonated = negatively charged).

**Worked example (real ECLIPSE binders, at pH 7.4):**
```python
net_charge("LISDAELEAIFAEELDC", 7.4)   # original, charged binder
# -> -6.129594662561429

net_charge("LISAAALAAIFAAALAC", 7.4)   # optimized binder (D/E -> A)
# -> -0.13278003501597901
```
The original binder carries 7 acidic residues (2×D, 4×E, 1×C) and zero
basic ones — strongly negative, as expected. After the neutralizing fix
(D/E→A), only the lone cysteine and the termini contribute meaningfully,
landing near-neutral (`-0.13`), matching the documented "near-neutral"
characterization of the optimized binder.

**Gotchas:**
- `C` (cysteine) is in `_ACIDIC`, alongside `D`/`E`/`Y` — that's because its
  thiol pKa (8.5) makes it weakly acidic at physiological pH, not a typo.
  Don't assume `_ACIDIC` only means the classic carboxylic-acid residues.
- The output is a **continuous** number, not an integer charge count — at
  pH values far from any residue's pKa it approaches integer-like values,
  but don't round it and expect to recover "number of charged residues"
  exactly, especially near a pKa.
- This function is pH-sensitive in a real, physically meaningful way — don't
  hardcode `7.4` everywhere if you're ever modeling a different buffer
  condition (e.g. an *E. coli* lysate or an in vitro assay at a different pH).

---

## `helix_propensity(sequence) -> float`

**What it does:** A single number estimating how "helix-friendly" a sequence
is, averaged residue-by-residue.

**Why it exists:** LOCKR binders are built from designed helical bundles —
if a substitution (like the D/E→A neutralizing fix) accidentally introduced a
helix-breaking residue, that would be a real structural problem even if the
charge/liability math looked fine. This is a cheap, instant sanity check
before anyone runs real structure prediction.

**Inputs:** `sequence: str`.

**Output:** `float` — the mean Chou-Fasman `P_alpha` value across all
residues. Values around `1.0` are "average," above `1.0` favor helix, below
`1.0` disfavor it. There's no hard "good"/"bad" cutoff baked into this
function itself — `analyze_charge` is where the `>= 1.0` threshold gets
applied as a decision rule.

**The math:** `sum(P_alpha[aa] for aa in sequence) / len(sequence)`, using the
1974 Chou-Fasman statistical propensity table (e.g. Glu `E=1.51`, highly
helix-favoring; Pro `P=0.57`, Gly `G=0.57`, the two classic helix-breakers).
Empty sequence returns `0.0` rather than dividing by zero.

**Worked example (real ECLIPSE binders):**
```python
helix_propensity("LISDAELEAIFAEELDC")   # original
# -> 1.2182352941176473
helix_propensity("LISAAALAAIFAAALAC")   # optimized
# -> 1.245294117647059
```
Both binders score comfortably above 1.0 — both are helix-friendly
sequences, and the neutralizing substitutions (E→A, D→A) actually *raised*
helix propensity slightly, since Ala (`A=1.42`) is itself a strong
helix-former, stronger than the Glu/Asp it replaced.

**Gotchas:** This is a *statistical propensity*, not a structure prediction —
it says nothing about whether the sequence will actually fold as designed in
3D context (register, packing, etc.). It's explicitly out of scope for this
module to do real structure prediction (that's the same "no folding here"
boundary `assembly.py`'s module docstring states even more explicitly).

---

## `helix_breakers(sequence) -> list[int]`

**What it does:** Finds every internal Proline or Glycine in the sequence —
the two residues most likely to physically kink a helix.

**Why it exists:** A single internal P or G is a much stronger red flag than
a slightly-below-average overall propensity score — `helix_propensity` could
look fine on average while still hiding one structurally fatal residue. This
function exists to catch that specific failure mode and tell you *exactly
where* it is.

**Inputs:** `sequence: str`.

**Output:** `list[int]` of 1-indexed positions where a P or G occurs
*internally* (not at either terminus). Empty list = no internal breakers.

**The math/logic:** `[i for i, aa in enumerate(seq, 1) if aa in "PG" and 1 < i < len(seq)]`
— the `1 < i < len(seq)` condition is what excludes the first and last
residue. Terminal P/G residues are excluded because a kink right at the very
end of a helix (where it's about to terminate anyway, or transitions into a
linker/cap) is usually structurally harmless — it's only kinks in the *middle*
of an intended helical run that break the secondary structure you were
counting on.

**Worked example (ECLIPSE binders — neither has any P/G at all):**
```python
helix_breakers("LISDAELEAIFAEELDC")   # -> []
helix_breakers("LISAAALAAIFAAALAC")   # -> []
```
Both ECLIPSE binders are clean — no P/G anywhere, internal or terminal. To
see this function actually catch something, here's a synthetic case from the
test suite:
```python
helix_breakers("AAALPAAALAAAL")   # P at position 5, neither terminus
# -> [5]
helix_breakers("PAAAAAAA")        # P at position 1 — the N-terminus itself
# -> []   (excluded: terminal, not internal)
```

**Gotchas:** 1-indexed positions, consistent with the rest of the codebase
(`BinderSequence.residues()`, `liability.py`'s mutation positions, etc.) —
but unlike `liability.py`'s mutation positions, these are *always*
binder-local (this function never sees an "absolute assembly coordinate"
concept at all — that's strictly an `assembly.py` concern).

---

## `analyze_charge(sequence, pH=7.4) -> ChargeResult`

**What it does:** The one-stop summary call — net charge, helix-friendliness
verdict, and a flagged list of helix-breaking positions, bundled into one
result.

**Why it exists:** This is the function a UI or a higher-level scan would
actually call — it packages the individual checks above into a single
pass/fail-style verdict (`helical_ok`) plus a human-readable `note` when
something looks off, so callers don't have to re-implement the "is this
helix-friendly" decision logic themselves.

**Inputs:** `sequence: str`, `pH: float` (default `7.4`).

**Output:** `ChargeResult` dataclass: `net_charge`, `pH`, `helical_ok: bool`,
`helix_breakers: list[int]`, `note: str` (empty string if everything looks
fine).

**The math/logic:** `helical_ok = helix_propensity(sequence) >= 1.0 and not
helix_breakers(sequence)` — both conditions must hold: average propensity at
or above the "neutral" threshold of 1.0, *and* zero internal P/G breakers.
If either fails, `note` is set to a fixed explanatory string; otherwise
`note` stays empty.

**Worked example (real ECLIPSE binders):**
```python
analyze_charge("LISDAELEAIFAEELDC")   # original
# ChargeResult(net_charge=-6.129594662561429, pH=7.4, helical_ok=True,
#              helix_breakers=[], note='')

analyze_charge("LISAAALAAIFAAALAC")   # optimized
# ChargeResult(net_charge=-0.13278003501597901, pH=7.4, helical_ok=True,
#              helix_breakers=[], note='')
```
Both come back `helical_ok=True` — the neutralizing fix changes the charge
dramatically (from `-6.13` to `-0.13`) without hurting the helix structure at
all, which is exactly why D/E→A is a sensible "fix the charge liability
without breaking the fold" substitution choice in the first place.

**Gotchas:** The `>= 1.0` cutoff for "helix-friendly" is a reasonable
rule-of-thumb threshold around the Chou-Fasman scale's own "neutral" point,
not a rigorously derived structural biology cutoff — treat `helical_ok=False`
as "worth a closer look," not as a hard proof the sequence won't fold.
