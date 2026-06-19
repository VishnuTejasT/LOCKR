I'm extending the LOCKR Biosensor Design Tool. Phase 1 (the thermodynamic
engine — thermo.py, charge.py, liability.py, models.py, calibration.py) is
done, tested (48 passing), committed, and reviewed for science correctness.
This is Phase 1.5: a new structural/assembly validation module that sits
ALONGSIDE the existing engine and does NOT touch it.

=== WHAT THIS MODULE IS FOR ===
The thermodynamic engine answers "given these numbers, what's the fold-change."
It has no idea whether a binder graft is structurally sane at the SEQUENCE
level — whether it overlaps something that must never be touched, whether it
fits the latch window, whether the final assembled sequence is what it should
be. This module is sequence-level bookkeeping, generalized from my own
ECLIPSE pipeline's manual verification scripts (Complete Documentation,
Script 4 and Script 6).

=== EXPLICITLY OUT OF SCOPE ===
No structure prediction, no folding, no PyRosetta/Rosetta, no GraftSwitchMover,
no AF2/AF3. This is pure string/position bookkeeping on sequences the user
already has (e.g. already graft-validated some other way). Do not add any
dependency beyond what Phase 1 already uses (numpy/scipy). Keep this CPU-only
and instant, like the rest of the engine.

=== GENERALIZATION — SAME PRINCIPLE AS PHASE 1 ===
My own system protects an 11-residue split-luciferase fragment called SmBiT
at a fixed latch position. But the TOOL must not know or care about "SmBiT"
specifically — it must support ANY protected motif at ANY position, because a
different LOCKR/lucCage team will have a different reporter system or
different critical motif to protect. SmBiT is MY example/test data, not a
hardcoded concept in the engine logic. Same pattern as Phase 1's
preserve_positions and charge_penalty_per_residue — general function,
ECLIPSE-specific defaults/test data only.

=== NEW MODULE: assembly.py ===

Core data needed (extend models.py, don't duplicate):
  - a ProtectedRegion: motif sequence (str), start position, end position
    (the user supplies this; no default assumes any specific motif)
  - a LatchWindow: start position, end position, expected total length
  - a GraftSpec: binder sequence, intended start position, optional linker
    sequence + position (for tandem/multi-copy grafts like my v2.2)

Functions, all general, no target/motif assumed:

  check_protected_region(full_sequence, protected_motif, start, end)
    -> verifies the protected_motif appears EXACTLY at [start:end] in
       full_sequence, unaltered. Returns intact: bool, found_sequence: str,
       mismatch_positions: list if not intact (so the user sees exactly what
       changed, not just pass/fail).

  check_graft_overlap(graft_spec, protected_region)
    -> verifies the binder (and linker, if present) graft positions do NOT
       overlap the protected region's position range. Returns overlap: bool,
       overlapping_positions: list if any.

  check_latch_fit(graft_spec, latch_window)
    -> verifies binder + linker (if present) length fits within the latch
       window's length. Returns fits: bool, used_length: int,
       available_length: int, slack: int (can be negative if it overflows).

  verify_full_assembly(full_sequence, latch_window, graft_spec,
                        protected_region, expected_total_length=None)
    -> the all-in-one check mirroring my Script 6 pattern, generalized:
       runs all of the above plus an overall length check, returns a single
       structured result object listing every individual check (name, passed,
       detail) so a UI can show a checklist, not just one boolean. Model this
       directly on the six-point structure of my Script 6 (length check,
       protected-motif check, spacer/structure check, binder check, linker
       check, second-binder/tandem check if present) but keep every check
       GENERAL — no hardcoded "SmBiT" or "spacer == 'DA'" anywhere; those are
       parameters or derived from the GraftSpec/ProtectedRegion the user gave.

=== BRIDGE FUNCTION: liability scanner <-> protected region ===
These are two genuinely different checks and must stay separate modules:
  - liability.py's preserve_positions protects TARGET-BINDING residues —
    a soft tradeoff (mutating them weakens affinity, scored continuously).
  - assembly.py's protected_region protects REPORTER-FUNCTION residues (my
    SmBiT case) — a hard constraint (mutating it kills signal entirely, no
    tradeoff, not a score).
Do NOT merge these concepts or make one a special case of the other. But add
ONE bridge function so the liability scanner's variant suggester can never
accidentally propose a substitution landing inside a protected region it
doesn't know about:

  filter_safe_variants(suggested_variants, protected_region)
    -> takes the list of suggested_variants that liability.py's suggester
       already produced (each with its substitution positions), and the
       protected_region from assembly.py. Returns only the variants whose
       substitution positions do NOT fall inside [protected_region.start,
       protected_region.end]. If a variant is filtered out, include it in a
       separate `rejected` list with the reason ("substitution at position X
       falls inside protected region"), don't just silently drop it.
  This function lives in assembly.py (it depends on assembly.py's
  ProtectedRegion type) and imports liability.py's variant result type from
  models.py — it does NOT require liability.py to know anything about
  protected regions. The dependency direction is one-way: assembly.py can
  know about liability.py's output shape, liability.py stays unaware of
  assembly.py. Confirm this import direction before writing the function.

  Test: construct a synthetic case where a suggested variant's substitution
  position deliberately overlaps a protected_region, and confirm
  filter_safe_variants correctly moves it to `rejected` with the right reason,
  while a non-overlapping variant in the same list passes through unchanged.
  Also test on my real ECLIPSE case: confirm none of the real suggested
  variants from the original->optimized binder fix ever get rejected (since in
  my actual design the binder and SmBiT don't overlap) — this should pass
  trivially but documents that the non-overlap was verified, not assumed.

=== TESTS ===

General tests (synthetic data, prove it works for ANY motif/sequence):
  - check_protected_region correctly detects an intact motif and correctly
    detects a corrupted one (single substitution), reporting the right
    mismatch position
  - check_graft_overlap correctly flags overlap and correctly clears a
    non-overlapping graft, on made-up positions unrelated to my system
  - check_latch_fit correctly flags a too-long graft and accepts one that fits,
    on a made-up latch window
  - verify_full_assembly runs end-to-end on synthetic data and produces a
    checklist with the right number of checks, right pass/fail per check

ECLIPSE validation tests (prove it reproduces my own documented verification,
calling the general functions with MY specific arguments as the worked
example):
  - SmBiT "VTGYRLFEEIL" at positions 312-322 in my full v1.0 sequence (379aa
    tagged or 359aa untagged — use whichever matches what's in the PDFs)
    checks as intact
  - latch window 325-359 (35 residues) with binder LISAAALAAIFAAALAC (17aa)
    at position 327 -> check_latch_fit reports it fits with the documented
    slack
  - v2.2 tandem: binder LISAAALAAIFAAALAC + "G" linker + LISAAALAAIFAAALAC
    (35aa total) exactly fills the 325-359 latch window -> fits with zero
    slack, matching my Script 5 ("v2.2 length: 359 aa OK")
  - verify_full_assembly on my v1.0 sequence reproduces the six checks from my
    Script 6 (length==359, SmBiT intact, spacer=='DA', binder1 correct, linker
    correct if tandem, binder2 correct if tandem) all passing
  - a deliberately corrupted SmBiT (mutate one residue) correctly fails the
    protected-region check and reports which position broke

=== CODE STYLE (same as Phase 1) ===
Comment the why, not the what. Reference my actual Script 4/Script 6
verification pattern in module-level comments since that's literally what
this generalizes. No identical docstring on every function. Real TODOs where
relevant (e.g. spacer-sequence checking is generalized loosely here, refine
if more LOCKR variants need it). No banners, no emoji.

=== WORKFLOW ===
Build assembly.py and its tests, run them, show me output. Don't touch
thermo.py, charge.py, liability.py, or calibration.py in this session — this
is additive, sitting alongside the existing engine. The only cross-module
dependency is assembly.py importing liability.py's variant result TYPE from
models.py for the bridge function's input — liability.py itself stays
completely unaware that assembly.py exists. If models.py needs a small
addition (e.g. exposing the variant result type cleanly for this import),
that's fine and expected — but don't change liability.py's own logic.
