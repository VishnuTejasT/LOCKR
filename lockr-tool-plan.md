# LOCKR Biosensor Design Tool — Project Plan

## 1. Vision

A CPU-only, locally-runnable design tool for LOCKR-style biosensors that answers the two questions that nearly broke ECLIPSE before they were solved:

1. **"Will this binder sequence preserve cage-key affinity, or quietly destroy it?"** — the charge-liability problem (six D/E residues killed K_CK on the original binder).
2. **"Given my measured/predicted Kds, what fold-change can I actually expect, and will tuning the latch help?"** — the thermodynamic-ceiling problem (fold-change capped by [lucKey]/K_CK, not K_open).

These are two views on the same physics — both turn on **K_CK, the cage-key dissociation constant** — so they ship as one tool with a shared thermodynamic engine, two analysis modes, and the ability to chain (scan a sequence → estimate its K_CK impact → predict the resulting fold-change).

**Distribution model:** downloadable, runs locally. The wiki has a page *about* the tool (what it does, screenshots, demo video, download/install link) pointing to a separate GitLab repo holding all source. Nothing live is embedded in the wiki, so the external-content rules are a non-issue. No GPU, no hosting, no cloud — it runs on a laptop.

---

## 2. The Two Modules

```
                 ┌─────────────────────────────────────┐
                 │   SHARED THERMODYNAMIC ENGINE         │
                 │   (one source of truth for K_CK math) │
                 └───────────────┬─────────────────────┘
                ┌────────────────┴────────────────┐
                ▼                                  ▼
   ┌──────────────────────────┐      ┌──────────────────────────┐
   │  MODULE A                 │      │  MODULE B                 │
   │  Charge / Liability       │      │  LOCKR Fold-Change        │
   │  Scanner                  │ ───► │  Calculator               │
   │                           │ K_CK │                           │
   │  sequence → liability     │ est. │  Kds + [key] → fold-change│
   │  → predicted K_CK penalty │      │  + ceiling diagnosis      │
   └──────────────────────────┘      └──────────────────────────┘
```

- **Module A** acts *before* a binder is finalized: scan a candidate sequence, flag charge liabilities, predict their effect on K_CK, suggest optimized variants.
- **Module B** acts *after* you have affinities: given K_CK, K_open, and [lucKey], predict the system's fold-change and tell the user whether they're key-limited or K_open-limited.
- **The chain** is the integration story: Module A's predicted K_CK penalty feeds directly into Module B, giving an end-to-end "sequence → expected sensor performance" path that mirrors exactly how ECLIPSE diagnosed and fixed its own binder.

This is the same shape as Lambert's CASPER (sequence in → scoring/rules engine → ranked guidance out), but built on your unique thermodynamic IP rather than a generative pipeline of public models.

---

## 3. Goals & Success Criteria (iGEM judging map)

| Goal | Success criterion | Judging criterion |
|------|-------------------|-------------------|
| Solves a real, generalizable problem | Any LOCKR/biosensor team can use it, not just ECLIPSE | "Useful to other projects?" |
| Validated by wet lab | TSA/BLI data on v1.0 & v2.2 validates Module A's liability calls and Module B's predictions | "Validated by experimental work?" |
| No-code usability | Researcher runs both analyses with no terminal | "User-friendly?" |
| Three access points | Local web UI + CLI + Python package on one engine | "Integrated with external tools/APIs?" |
| Standards awareness | FASTA in, CSV/JSON out; objects map to SBOL 3 concepts | "Compatible with synbio standards?" |
| Extensible & documented | Install/usage docs, worked examples, clean modules | "Documented for future groups?" |

**Primary differentiator:** the validation framing is clean and direct — *wet lab outcomes validate the tool's predictions.* Module A predicted the original binder's six D/E residues were a liability; the charge-optimized variant recovered function — that's a built-in before/after case study. Module B predicted the [lucKey]/K_CK ceiling; your fold-change measurements test it. No detour through a generative pipeline.

---

## 4. Module A — Charge / Liability Scanner

### 4.1 Purpose
Catch the failure mode that broke the original binder: acidic residues (D/E) in the grafted latch region disrupting the cage's structure/electrostatics and collapsing K_CK. Flag liabilities, score them, suggest fixes — *before* anyone spends money synthesizing or time folding.

### 4.2 Inputs

| Input | Format | Required | Notes |
|-------|--------|----------|-------|
| Binder sequence(s) | single peptide or FASTA (batch) | yes | one-letter amino acids |
| Structurally sensitive window | residue range | no | the latch/grafting region where charges matter most; defaults to whole sequence |
| pH | number | no | default 7.4 for net-charge calc |
| Substitution policy | preset | no | conservative (D→N, E→Q) vs neutralizing (D→A, E→A) |

### 4.3 Algorithm (the engine)
1. **Net charge at pH** — Henderson-Hasselbalch over D, E, K, R, H side chains + termini.
2. **Acidic-residue census** — count and locate every D/E; mark which fall inside the sensitive window.
3. **Position-weighted liability score** — each acidic residue weighted by structural sensitivity (heavier weight inside the latch/grafting window, optionally by helical face so charges on the cage-facing side score worse). Sum → a single liability score per sequence.
4. **Predicted K_CK penalty** — map liability score → estimated K_CK degradation, calibrated against your two anchor points: `LISDAELEAIFAEELDC` (6 D/E, K_CK destroyed) vs `LISAAALAAIFAAALAC` (0 acidic, K_CK recovered). *Document this as a two-point calibration that future data refines — don't oversell it as fully fit.*
5. **Helical-propensity check** — these are α-helical latch grafts; flag if a proposed substitution would disrupt helicity (simple Chou-Fasman / position scoring is enough).
6. **Variant suggestions** — propose charge-optimized sequences using the chosen substitution policy, preserving hydrophobic core and helix (the same logic that produced the LISAAA... fix).

### 4.4 Outputs
- Annotated sequence (acidic residues highlighted, sensitive window marked)
- Liability score + predicted K_CK penalty band (e.g. "high — likely to collapse K_CK")
- Flagged positions with per-residue contribution
- One or more suggested optimized variants with their recomputed scores
- **Batch mode:** ranked table of candidate sequences (lowest liability first) — exportable CSV, so a list of ProteinMPNN-style designs can be triaged before folding/ordering

### 4.5 UI (Scanner tab)
- **Input area:** big sequence box (paste) or file upload (FASTA batch); a sequence-position ruler so the user can drag-select the sensitive window
- **Live annotation:** as the user types/pastes, the sequence renders with acidic residues colored, sensitive window shaded; net charge and liability score update live
- **Result panel:** liability gauge (green/amber/red), predicted K_CK impact, per-position bar showing each residue's contribution
- **Suggestions panel:** side-by-side original vs suggested variant(s), with a diff highlight on changed residues and the score improvement
- **Batch view:** sortable/filterable table (like CASPER's candidate list), CSV export

---

## 5. Module B — LOCKR Fold-Change Calculator

### 5.1 Purpose
Predict sensor dynamic range from thermodynamic parameters, and — the headline feature — diagnose *what's limiting it*, so users don't waste effort on latch mutations that can't help.

### 5.2 Thermodynamic framework
The lucCage/lucKey switch is a thermodynamically coupled conformational system:

- Cage interconverts **closed ⇌ open**, governed by `K_open` (tuned by latch mutations; target binding shifts it toward open).
- **lucKey binds only the open cage**, governed by the cage-key dissociation constant `K_CK`.
- **Signal** (split-luciferase complementation) is proportional to the cage·key complex.
- **Fold-change** = signal(target present) / signal(target absent).

Because key binds only open cages, the *apparent* key affinity is weakened by the closed-state population — so K_open and K_CK are coupled. Key occupancy follows a saturable form `θ = [key] / ([key] + K_app)`, with `K_app` depending on both K_CK and K_open; fold-change is the ratio of θ between the target-on and target-off states. In the key-saturating regime this ratio is bounded above by **[lucKey] / K_CK**, independent of K_open — which is exactly the ceiling that made latch mutations ineffective for ECLIPSE at 500 nM lucKey ([lucKey]/K_CK = 50).

> **Implementation note:** the engine should implement *your validated ECLIPSE fold-change expression* as the source of truth — the framework above is the scaffold it sits in. Ship the `[lucKey]/K_CK` ceiling as an explicit diagnostic computed alongside the full prediction, since that's the insight that generalizes to other teams.

### 5.3 Inputs

| Input | Symbol | Required | Notes |
|-------|--------|----------|-------|
| Cage-key dissociation constant | K_CK | yes | can be typed, or piped from Module A's prediction |
| Latch opening equilibrium (OFF / basal) | K_open | yes | the closed↔open balance without target |
| Target-induced shift | ΔK_open or K_open(ON) | yes | how much target binding stabilizes the open state |
| lucKey concentration | [lucKey] | yes | accepts a single value or a range to sweep |
| Target affinity / concentration | K_target, [target] | optional | for modeling partial target occupancy |

### 5.4 Outputs
- **Predicted fold-change** for the given parameters
- **Ceiling value** [lucKey]/K_CK, shown next to the prediction
- **Regime diagnosis** — the killer feature — one of:
  - *"Key-limited: you're at/near the [lucKey]/K_CK ceiling. Latch mutations won't help. To improve fold-change, raise [lucKey] or improve K_CK."*
  - *"K_open-limited: you have headroom below the ceiling. Tuning the latch (K_open) will improve fold-change."*
- **Sweep plots:**
  - Fold-change vs [lucKey] (shows where the current setup sits relative to the ceiling)
  - Fold-change vs K_open (shows whether latch tuning moves the needle in this regime)
  - Optional 2D heatmap: fold-change over ([lucKey], K_CK) or (K_open, K_CK)

### 5.5 UI (Calculator tab)
- **Parameter panel:** numeric inputs with sensible units (nM) and inline help; a toggle to switch any single parameter into "sweep a range" mode
- **Result card:** big predicted fold-change number, the ceiling value beside it, and a plain-language regime verdict in a colored banner (this is the part a non-coder reads first)
- **Plots:** live-updating sweep charts (Recharts/Plotly) that redraw as parameters change; a marker showing the user's current operating point on each curve
- **"Pipe from Scanner" button:** pull K_CK estimate straight from Module A's last result

---

## 6. Shared Thermodynamic Engine

One Python module is the single source of truth for all the math — net charge, liability scoring, K_CK estimation, the fold-change model, the ceiling diagnostic. Both modules, the CLI, and the web UI all call into it. No duplicating equations across frontend/backend.

```
engine/
  charge.py        # net charge, pKa, helical propensity
  liability.py     # position-weighted scoring, K_CK penalty calibration, variant suggestion
  thermo.py        # K_open/K_CK coupling, fold-change model, ceiling diagnostic
  models.py        # data objects (Binder, Sequence, SensorParams, Result)
  calibration.py   # anchor data + fit (original vs optimized binder)
```

---

## 7. Architecture

CPU-only and lightweight enough that the cleanest design is **one Python core powering three thin clients**, distributed for local use.

```
┌─────────────────────────────────────────────────────────┐
│  CORE PYTHON ENGINE  (the math — see §6)                  │
└───────────────┬───────────────┬───────────────┬─────────┘
                │               │               │
        ┌───────▼──────┐ ┌──────▼──────┐ ┌──────▼────────┐
        │  CLI          │ │  Python      │ │  Local web UI  │
        │  lockr scan   │ │  package     │ │  FastAPI +     │
        │  lockr fc     │ │  import ...  │ │  React         │
        └──────────────┘ └─────────────┘ └───────────────┘
                                                │
                                    `lockr serve` → localhost:8000
```

### Why local web UI (not hosted)
The user downloads the tool and runs `lockr serve`, which starts a FastAPI server on localhost and opens the React UI in their browser. This gives the polished no-code experience **with zero hosting, zero GPU, and zero external-content concerns** — every user runs their own copy. The wiki just links to the repo and shows a demo video. This matches your stated distribution model exactly.

> **Alternative:** a pure-static frontend that reimplements the math in JS (no backend at all) is even simpler to host but duplicates the engine in two languages. Not worth it — keep one Python engine as truth and serve it locally.

---

## 8. Backend (FastAPI, local + lightweight)

Thin wrapper over the engine; runs on localhost when the user launches the UI.

| Method | Route | Purpose |
|--------|-------|---------|
| `POST` | `/scan` | sequence(s) + options → liability results |
| `POST` | `/foldchange` | thermo params → fold-change + ceiling + regime |
| `POST` | `/sweep` | param + range → series for plotting |
| `POST` | `/suggest` | sequence → optimized variant(s) |

No job queue, no database, no async needed — every computation is sub-second. Plain synchronous request/response. Runs in the `igem` conda env (Python 3.10).

---

## 9. Frontend (React)

Two tabs (Scanner / Calculator) plus a small "Chain" affordance connecting them. Built with the same stack you already use for the wiki editor.

### Layout
- **Top:** tool name, two-tab switcher (Scanner | Calculator), link to docs
- **Scanner tab:** §4.5
- **Calculator tab:** §5.5
- **Chain affordance:** a "→ send K_CK to Calculator" action on a Scanner result, and a "← pull K_CK from Scanner" button on the Calculator

### Components
- 3D viewer is **not** needed (this is sequence + numbers, not structures) — big simplification over the old pipeline-tool plan
- Sequence renderer with per-residue coloring (Scanner)
- Live-updating charts via Recharts or Plotly (Calculator sweeps)
- Result cards with plain-language verdicts (the no-coder reads these, not the equations)
- CSV/JSON export buttons on both tabs

### Design principles
Clear defaults, inline help text, advanced inputs tucked away. Newcomer reaches a verdict in seconds; an expert can sweep parameters and export batches. Same usability philosophy CASPER used.

---

## 10. Data Formats & Standards

- **Inputs:** plain sequence or FASTA (Module A); numeric params, optionally a small JSON/YAML config (Module B, for reproducible batch runs)
- **Outputs:** CSV (ranked candidates, sweep tables), JSON (full results), PNG/SVG (plots)
- **SBOL:** don't implement it; document how the objects (binder sequence, liability annotation, sensor parameter set, fold-change result) map onto SBOL 3 concepts. Satisfies the standards criterion the way CASPER did.

---

## 11. Validation Plan

Mirror CASPER's two-pronged approach — and note that *the wet lab is ground truth; the tool's predictions are what's being tested.*

1. **Retrospective (Module A):** the original `LISDAELEAIFAEELDC` vs charge-optimized `LISAAALAAIFAAALAC` is a built-in before/after. The scanner should flag the original as high-liability and the optimized as clean — and that ranking is validated by the fact that one collapsed K_CK and the other recovered function.
2. **Prospective (both modules):** the TSA/BLI data coming back from the Twist-ordered v1.0 and v2.2 binders. Module A's liability calls and Module B's fold-change/ceiling predictions get checked against measured binding and dynamic range. Document each comparison with a figure.

This is the strongest part of the submission — present it as: *the tool predicted X, the wet lab measured Y, here's the agreement.*

---

## 12. Documentation Plan

For the separate GitLab repo (scored deliverable — you already have the discipline from your 16+ script PDFs):

- `README.md` — what it is, screenshots, demo GIF
- `INSTALL.md` — conda env, `pip install`, `lockr serve` quickstart (CPU-only, no GPU note front and center)
- `USAGE.md` — web walkthrough + CLI flags + Python API, with worked examples (the ECLIPSE binder as the worked example)
- `SCIENCE.md` — the thermodynamic model, the K_CK coupling, the ceiling derivation, the liability scoring rationale and calibration
- `API.md` — endpoint reference
- Clean in-code comments + organized modules; an install script

---

## 13. Repo Structure

```
lockr-tool/
├── README.md
├── INSTALL.md  USAGE.md  SCIENCE.md  API.md
├── LICENSE                      # CC-BY 4.0 / OSI as iGEM requires
├── pyproject.toml               # pip-installable package
├── engine/                      # shared math (see §6)
├── cli/                         # lockr scan / lockr fc / lockr serve
├── api/                         # FastAPI app
├── web/                         # React frontend
├── examples/                    # worked examples, sample FASTA, configs
├── tests/                       # unit tests on the engine (anchor data, known cases)
└── data/                        # calibration anchors
```

---

## 14. Roadmap

Phased so there's always something that runs.

### Phase 0 — Foundations
- [ ] Repo scaffold, license, package skeleton
- [ ] Confirm the exact fold-change equation + parameter definitions from your ECLIPSE notes
- [ ] Lock the engine's data objects (Sequence, SensorParams, Result)

### Phase 1 — Engine (do first)
- [ ] `charge.py` — net charge, pKa, helical propensity
- [ ] `liability.py` — scoring + K_CK penalty calibration on the two anchor sequences
- [ ] `thermo.py` — fold-change model + [lucKey]/K_CK ceiling + regime diagnosis
- [ ] Unit tests: original vs optimized binder ranks correctly; ceiling reproduces the =50 result
- [ ] **Milestone:** engine reproduces both ECLIPSE findings from a script

### Phase 2 — CLI + package
- [ ] `lockr scan`, `lockr fc`, `lockr sweep` with flags + YAML config
- [ ] `pip install` works; importable API
- [ ] **Milestone:** both analyses runnable from terminal and a notebook

### Phase 3 — Local web UI
- [ ] FastAPI wrapper (`/scan`, `/foldchange`, `/sweep`, `/suggest`)
- [ ] React: Scanner tab (live annotation, suggestions, batch table)
- [ ] React: Calculator tab (params, verdict banner, live sweep plots)
- [ ] Chain affordance between tabs
- [ ] `lockr serve` launches it on localhost
- [ ] **Milestone:** non-coder runs both analyses in the browser, exports CSV

### Phase 4 — Validation + polish
- [ ] Retrospective validation writeup (original vs optimized)
- [ ] Fold in TSA/BLI data as it arrives; validation figures
- [ ] SBOL mapping doc, full documentation set
- [ ] Demo video + wiki page about the tool

---

## 15. Tech Stack

| Layer | Choice | Why |
|-------|--------|-----|
| Engine | Python 3.10 (`igem` env) | matches your setup; CPU-only |
| Charge/pKa math | NumPy (+ optional Biopython for parsing) | lightweight |
| Plots | Recharts or Plotly (frontend) / matplotlib (CLI export) | live sweeps + static figures |
| Backend | FastAPI (local only) | auto OpenAPI docs; CASPER-aligned |
| Frontend | React | you already use it (wiki editor) |
| Packaging | pyproject + pip; YAML configs | CLI/package access points |
| Repo | separate GitLab repo, linked from wiki | matches iGEM software deliverable model |

---

## 16. Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| K_CK penalty calibrated on only 2 sequences | Present honestly as a two-point anchor; refine with v1.0/v2.2 data; expose it as a tunable model, not a hard claim |
| Fold-change model is ECLIPSE-specific | Keep the engine parameterized; document assumptions in SCIENCE.md so other teams can adapt |
| Tool seen as "just a calculator" | Lead with the *diagnosis* (regime verdict, "latch mutations won't help") and the wet-lab validation — that's insight, not arithmetic |
| Scope creep | Engine-first milestone (Phase 1) keeps a working, testable artifact immediately |
| Helical-propensity scoring adds complexity | Ship v1 with charge-only liability; add helix check as a refinement if time allows |

---

## 17. Open Questions

- Exact fold-change equation + precise definitions of K_open(OFF) vs K_open(ON) / ΔK_open as you derived them?
- For Module A, do you want helical-face weighting in v1, or charge-position-only first?
- Should the K_CK penalty model output an absolute K_CK estimate, or just a relative liability band? (Affects how cleanly it pipes into Module B.)
- Single binder per run vs batch as the default mode for the Scanner?
- Tool name — keep LOCKSMITH or something else?
