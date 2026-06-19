# LOCKR Biosensor Design Tool — UI/UX & Design Specification

Companion to the project plan. This is the build-level spec: design tokens, screen layouts, every field with its control/default/validation, exact API request/response schemas, plot specifications, the cross-tab chain mechanism, error states, and the actual interface copy. Someone should be able to build the front end exactly as intended from this document.

---

## Table of Contents
1. Design System
2. Global Layout & Shell
3. Module A — Scanner (full screen spec)
4. Module B — Calculator (full screen spec)
5. The Chain (cross-tab state)
6. API Schemas
7. Plot Specifications
8. Error & Edge States
9. Interface Copy (verbatim)
10. Responsive Behavior
11. Accessibility
12. Component Inventory & State Model

---

## 1. Design System

Reuses the ECLIPSE wiki identity so the tool reads as part of the same project. Dark-red academic palette, serif display type, clean sans for UI.

### 1.1 Color tokens

```
/* Brand */
--brand-900:  #5e0808;   /* deepest red — hover/active on primary */
--brand-700:  #870b0b;   /* PRIMARY brand red (matches wiki) */
--brand-500:  #b91c1c;   /* lighter red — accents, focus rings */
--brand-100:  #f6e7e7;   /* tint — selected rows, subtle fills */
--brand-050:  #fbf3f3;   /* faintest tint — hover backgrounds */

/* Neutrals (warm, paper-like) */
--canvas:     #fbfaf8;   /* page background */
--surface:    #ffffff;   /* cards, panels */
--surface-alt:#f5f2ee;   /* insets, code blocks, table header */
--border:     #e5e0d8;   /* hairlines, dividers */
--border-strong:#cfc8bc; /* input borders */

--text:       #1c1917;   /* primary text */
--text-soft:  #57534e;   /* secondary text, labels */
--text-muted: #a8a29e;   /* placeholder, captions, disabled */

/* Semantic (liability bands, verdicts, validation) */
--success-700:#15803d;  --success-100:#dcfce7;   /* low liability / good */
--warning-700:#b45309;  --warning-100:#fef3c7;   /* moderate */
--danger-700: #b91c1c;  --danger-100: #fee2e2;   /* high liability / error */
--info-700:   #870b0b;  --info-100:   #f6e7e7;   /* neutral callouts (brand) */

/* Data viz */
--plot-line:    #870b0b;  /* primary curve (fold-change) */
--plot-ceiling: #b45309;  /* ceiling reference (dashed) */
--plot-marker:  #1c1917;  /* operating-point dot */
--plot-grid:    #ece8e1;  /* gridlines */
--plot-zone-ok: #dcfce7;  --plot-zone-warn:#fef3c7;  --plot-zone-bad:#fee2e2;
```

> The liability gauge uses success→warning→danger. Note `--danger-700` and `--brand-500` are both reds; keep the danger red for *liability/error only* and brand red for *identity/navigation*, never adjacent in the same component.

### 1.2 Typography

| Role | Font | Size / line-height | Weight |
|------|------|--------------------|--------|
| Display / hero numbers | Cormorant Garamond | 40 / 48 | 600 |
| H1 (screen title) | Cormorant Garamond | 32 / 40 | 600 |
| H2 (panel title) | DM Sans | 22 / 28 | 600 |
| H3 (field group) | DM Sans | 16 / 22 | 600 |
| Body / labels | DM Sans | 15 / 22 | 400–500 |
| Small / help text | DM Sans | 13 / 18 | 400 |
| Caption / units | DM Sans | 12 / 16 | 500 |
| **Sequence & numeric data** | **IBM Plex Mono** | 15 / 24 | 400–500 |

Cormorant Garamond and DM Sans/Inter already load on the wiki. **Add IBM Plex Mono** — amino-acid sequences must be monospace so positions align under the ruler, and tabular Kd values read cleaner in mono.

### 1.3 Spacing, radii, elevation

```
--space: 4 8 12 16 24 32 48 64        /* 4px base scale */
--radius-sm: 4px;  --radius-md: 8px;  --radius-lg: 12px;  --radius-pill: 999px;
--shadow-sm: 0 1px 2px rgba(28,25,23,.06);
--shadow-md: 0 4px 12px rgba(28,25,23,.08);
--focus-ring: 0 0 0 3px rgba(185,28,28,.35);   /* brand-500 @ 35% */
```

### 1.4 Component primitives
- **Buttons:** primary (brand-700 bg, white text, brand-900 hover), secondary (surface bg, border-strong, text), ghost (text-only). Radius `md`, height 40px, 15px DM Sans 500.
- **Inputs:** surface bg, 1px border-strong, radius `md`, height 40px, focus → brand-500 border + focus-ring. Numeric inputs right-align, mono font.
- **Cards/panels:** surface bg, 1px border, radius `lg`, shadow-sm, 24px padding.
- **Tabs:** underline style, brand-700 active indicator, text-soft inactive.
- **Badges/pills:** radius-pill, 12px caption, semantic bg-100 + text-700.
- **Tooltips:** dark (#1c1917) bg, white text, 12px, on `?` icons next to labels.

---

## 2. Global Layout & Shell

```
┌────────────────────────────────────────────────────────────────┐
│  HEADER (sticky, 64px)                                            │
│  [◧ LOCKSMITH]            [ Scanner | Calculator ]      [Docs ↗]  │  ← tabs centered
├────────────────────────────────────────────────────────────────┤
│                                                                  │
│   ACTIVE TAB CONTENT (max-width 1200px, centered, 32px gutters)  │
│                                                                  │
├────────────────────────────────────────────────────────────────┤
│  FOOTER  · CC-BY 4.0 · GitLab repo link · version                │
└────────────────────────────────────────────────────────────────┘
```

- **Header:** logo/wordmark left (Cormorant Garamond 22, brand-700), tab switcher centered (underline tabs), Docs link right (opens repo docs in new tab). Sticky, surface bg, 1px bottom border, shadow-sm on scroll.
- **Tab switcher** is the primary nav between the two modules. Active tab persists in URL (`/#scanner`, `/#calculator`) so a refresh keeps you in place.
- **Footer:** required iGEM CC-BY 4.0 notice, link to the separate GitLab repo, tool version string.
- **Canvas bg** `--canvas`; content panels sit on it as `--surface` cards.

---

## 3. Module A — Scanner (full screen spec)

### 3.1 Layout (desktop, two columns)

```
┌─ INPUT (left, ~46%) ───────────┐  ┌─ RESULTS (right, ~54%) ──────────┐
│ H1  Charge & Liability Scanner │  │  (empty state until first scan)   │
│ helper line                    │  │                                   │
│                                │  │  ┌ Liability summary card ──────┐ │
│ [ Single | Batch ] segmented   │  │  │ gauge · net charge · K_CK    │ │
│                                │  │  └──────────────────────────────┘ │
│ ┌ Sequence ───────────────┐    │  │  ┌ Annotated sequence ──────────┐ │
│ │ (monospace textarea)     │    │  │  │ colored residue chips + ruler│ │
│ │                          │    │  │  └──────────────────────────────┘ │
│ └──────────────────────────┘    │  │  ┌ Per-position contribution ───┐ │
│ charcount · validity            │  │  │ bar chart                    │ │
│                                │  │  └──────────────────────────────┘ │
│ ▸ Advanced                     │  │  ┌ Suggested variant ───────────┐ │
│   pH [7.4]                      │  │  │ original ▸ optimized diff     │ │
│   Sensitive window [slider]     │  │  │ score improvement  [→ Calc]   │ │
│   Substitution [Conserv|Neutr]  │  │  └──────────────────────────────┘ │
│                                │  │                                   │
│ [ Scan ]  [ Reset ]            │  │                                   │
└────────────────────────────────┘  └───────────────────────────────────┘
```

### 3.2 Fields

| # | Field | Control | Default | Validation |
|---|-------|---------|---------|------------|
| A1 | Mode | Segmented `Single / Batch` | Single | — |
| A2 | Sequence (Single) | Monospace textarea, auto-uppercase, whitespace stripped live | empty | Only `ACDEFGHIKLMNPQRSTVWY`; 1–200 residues; non-empty on submit |
| A3 | Sequences (Batch) | FASTA file upload + drag-drop, or paste FASTA | none | Valid FASTA; each record passes A2 rules; ≤500 records |
| A4 | pH | Number input | `7.4` | 0–14, step 0.1 |
| A5 | Sensitive window | Dual-handle range slider over the position ruler + numeric start/end | full length (1–N) | start ≤ end, within 1–N |
| A6 | Substitution policy | Segmented `Conservative / Neutralizing` | Conservative | — |
| A7 | Scan | Primary button | — | disabled until A2/A3 valid |
| A8 | Reset | Ghost button | — | clears inputs + results |

- **Conservative** = D→N, E→Q (charge removed, shape/H-bonding kept). **Neutralizing** = D→A, E→A (the literal ECLIPSE fix). Tooltip on the control explains both.
- The **sensitive window** slider sits directly above/below the sequence so handles line up with residues. Residues outside the window are dimmed in the annotated output.

### 3.3 Live behavior (before clicking Scan)
- As the user types in A2, render a **live position ruler** (every 5th position numbered) and color acidic residues (D/E) in `--danger-700`, basic (K/R/H) in `--brand-500`, others in `--text`. Show a live **net-charge readout** and **count of D/E** in the field footer. Invalid characters are rejected on input (not inserted) with a brief shake + caption "only standard amino acids."
- This gives instant feedback even before a full scan.

### 3.4 Results (after Scan)

**A. Liability summary card**
- **Liability gauge:** horizontal segmented bar, 0–100, with band thresholds (Low 0–33 green, Moderate 34–66 amber, High 67–100 red) and a needle/marker at the score. Big band label (Cormorant 32): "High liability".
- **Net charge** at the chosen pH (mono, signed).
- **Predicted K_CK penalty** band with one-line plain text (see copy §9): e.g. "Severe — likely to collapse cage-key affinity."

**B. Annotated sequence**
- Residue chips in mono, each a small rounded box; acidic residues filled `--danger-100` with `--danger-700` text, flagged (in-window) acidic residues get a solid red underline. Position ruler beneath. Hover a chip → tooltip with that residue's contribution and pKa.

**C. Per-position contribution**
- Horizontal bar chart, one row per residue (or per flagged residue in long sequences), bar length = liability contribution, colored by band. Helps the user see *which* residues to fix.

**D. Suggested variant**
- Original vs suggested shown as an aligned mono diff: unchanged residues `--text-muted`, substituted residues highlighted (old struck/red, new green). Shows recomputed liability score + delta ("87 → 8"). 
- Two actions: **Copy variant** and **→ Send K_CK estimate to Calculator** (see §5). If multiple variants, a small "1 / 3" stepper.

**Batch mode results:** replace B–D with a **sortable table** — columns: ID, length, net charge, liability score, band (colored pill), predicted K_CK penalty, top suggested variant. Sort by any column; default sort = liability ascending (cleanest first). Row click expands to the single-sequence detail (B–D). **Export CSV** button top-right.

---

## 4. Module B — Calculator (full screen spec)

### 4.1 Layout (desktop, two columns)

```
┌─ PARAMETERS (left, ~40%) ──────┐  ┌─ RESULTS (right, ~60%) ──────────┐
│ H1  Fold-Change Calculator     │  │ ┌ Verdict card ─────────────────┐ │
│ helper line                    │  │ │ FOLD-CHANGE   38×              │ │
│                                │  │ │ ceiling 50× · 76% of ceiling   │ │
│ ┌ Core parameters ─────────┐   │  │ │ ▣ banner: regime verdict       │ │
│ │ K_CK        [10]  nM  ⇄  │   │  │ └────────────────────────────────┘ │
│ │ K_open OFF  [0.01]    ?  │   │  │ ┌ Fold-change vs [lucKey] ──────┐ │
│ │ K_open ON   [1.0]     ?  │   │  │ │ log-x line plot + ceiling line │ │
│ │ [lucKey]    [500] nM  ~  │   │  │ │ + operating-point marker       │ │
│ └──────────────────────────┘   │  │ └────────────────────────────────┘ │
│ ▸ Advanced (target binding)    │  │ ┌ Fold-change vs K_open ────────┐ │
│   K_target  [—] nM             │  │ │ shows if latch tuning helps    │ │
│   [target]  [—] nM             │  │ └────────────────────────────────┘ │
│                                │  │ ┌ Recommendations ──────────────┐ │
│ [ Calculate ]  [ Reset ]       │  │ │ bulleted, regime-specific      │ │
│                                │  │ └────────────────────────────────┘ │
└────────────────────────────────┘  └───────────────────────────────────┘
```

### 4.2 Fields

| # | Field | Symbol | Control | Default | Units | Validation |
|---|-------|--------|---------|---------|-------|------------|
| B1 | Cage-key dissociation const | K_CK | Number, sci-notation aware, mono | `10` | nM | > 0 |
| B2 | Latch opening (basal/OFF) | K_open(OFF) | Number | `0.01` | — (dimensionless) | > 0 |
| B3 | Latch opening (target/ON) | K_open(ON) | Number | `1.0` | — | > 0; UI warns if ON ≤ OFF |
| B4 | lucKey concentration | [lucKey] | Number, sweepable (`~`) | `500` | nM | > 0 |
| B5 | Target affinity (adv.) | K_target | Number | empty | nM | > 0 or empty |
| B6 | Target concentration (adv.) | [target] | Number | empty | nM | ≥ 0 or empty |
| B7 | Calculate | — | Primary button | — | — | enabled when B1–B4 valid |
| B8 | Reset | — | Ghost button | — | — | restores defaults |

- **Defaults are placeholders** anchored to ECLIPSE's regime ([lucKey] 500 nM, K_CK 10 nM → ceiling 50). Replace with real values once confirmed; mark them as examples in help text.
- **Per-field affordances:** `⇄` on K_CK = "pull from Scanner" (chain, §5). `~` on a field = toggle **sweep mode**: the single input expands into `min / max / steps` + a `linear · log` switch (default **log**). B3 toggle lets the user enter ON value *or* ΔK_open.
- Numeric inputs accept scientific notation (`4.2e-5`) and display values in nM with automatic unit scaling in readouts (fM/pM/nM/µM) so 42 fM doesn't render as `0.000042`.

### 4.3 Results

**A. Verdict card** (the thing a non-coder reads first)
- **Hero fold-change** (Cormorant 40): `38×`.
- Sub-line: dominance ratio `[lucKey]/K_CK = 50` and `% of dominance ratio` — **not labeled "ceiling."** This number is a diagnostic threshold (it tells you whether K_open changes will move fold-change), not an achievable fold-change target. The true asymptotic max as pull→∞ is `(1+K_open+lucKey/K_CK)/K_open`, a much larger and practically irrelevant number given realistic pull values (10–30) — don't surface it as a headline figure; it would mislead a wet-lab user into thinking 50× (or whatever the ratio is) is the ceiling they're aiming for, when realistic fold-change at typical pull is much closer to `(1+pull)`.
- **Regime banner** — full-width colored callout with icon + plain-language verdict (exact copy §9.3). Three states: `key-limited` (amber), `K_open-limited` (green/info), `mixed` (neutral).

**B. Plots** — see §7. Primary: fold-change vs [lucKey] (log-x) with a dashed reference line at the dominance ratio and a marker at the operating point. Secondary: fold-change vs K_open, which visually demonstrates whether latch tuning moves the curve (in the key-limited regime it's flat — that's the point).

**C. Recommendations** — regime-specific bullets (copy §9.3), e.g. for key-limited: raise [lucKey], improve K_CK; *don't* bother with latch mutations.

---

## 5. The Chain (cross-tab state)

The integration feature: Scanner's predicted K_CK feeds the Calculator.

### Mechanism
- A single app-level store (React Context `AnalysisChainContext`) holds:
  ```ts
  { pipedKck: number | null,        // nM
    source: { module: 'scanner', label: string } | null }
  ```
- **In Scanner:** the suggested-variant card (and the summary card) has a **"→ Send K_CK to Calculator"** button. Clicking it computes the K_CK estimate for that sequence/variant, writes `{ pipedKck, source: {label: 'binder_v1 (optimized)'} }`, and switches to the Calculator tab.
- **In Calculator:** if `pipedKck` is set, field B1 (K_CK) pre-fills with that value and shows a small brand-100 pill beside the label: **"from Scanner · binder_v1 (optimized) ✕"**. The field stays **editable**; editing it or clicking ✕ detaches the link (clears `source`, keeps the value).
- The `⇄` icon on B1 also opens a tiny menu: "Use latest Scanner estimate" / "Enter manually."

### UX rules
- Piping never silently overwrites a value the user has manually typed — if B1 was edited, show a confirm ("Replace your K_CK of X with Scanner estimate Y?").
- The pill makes provenance obvious so a user always knows whether a number is theirs or derived.

---

## 6. API Schemas

All endpoints are local (`http://localhost:PORT`), synchronous, JSON in/out. Units: concentrations in **nM**; K_open dimensionless.

### 6.1 `POST /scan`
**Request**
```json
{
  "sequences": [
    { "id": "binder_v1", "sequence": "LISDAELEAIFAEELDC" }
  ],
  "sensitive_window": { "start": 1, "end": 17 },
  "ph": 7.4,
  "substitution_policy": "conservative"
}
```
**Response**
```json
{
  "results": [
    {
      "id": "binder_v1",
      "sequence": "LISDAELEAIFAEELDC",
      "length": 17,
      "net_charge": -5.9,
      "acidic_residues": [
        { "position": 4, "residue": "D", "in_window": true, "contribution": 0.18 }
      ],
      "liability_score": 87,
      "liability_band": "high",
      "predicted_kck_penalty": { "band": "severe", "note": "likely to collapse K_CK" },
      "per_position": [
        { "position": 1, "residue": "L", "contribution": 0.0 }
      ],
      "helix_flags": [
        { "position": 4, "issue": "acidic in helical core face" }
      ],
      "suggested_variants": [
        {
          "sequence": "LISNAALAAIFAANLNC",
          "substitutions": [ { "position": 4, "from": "D", "to": "N" } ],
          "liability_score": 8,
          "liability_band": "low",
          "estimated_kck_nm": 9.5
        }
      ]
    }
  ]
}
```

### 6.2 `POST /foldchange`
**Request**
```json
{
  "k_ck": 10.0,
  "k_open_off": 0.01,
  "k_open_on": 1.0,
  "luckey": 500.0,
  "k_target": null,
  "target_conc": null
}
```
**Response**
```json
{
  "fold_change": 38.2,
  "dominance_ratio": 50.0,
  "fraction_of_dominance_ratio": 0.764,
  "regime": "key_limited",
  "limiting_factor": "luckey_over_kck",
  "verdict": "You are near the lucKey/K_CK dominance threshold...",
  "recommendations": [
    "Increase lucKey to raise the dominance ratio.",
    "Improve cage-key affinity (lower K_CK).",
    "Latch (K_open) mutations will not raise fold-change here."
  ]
}
```
Note: `dominance_ratio` (= `lucKey/K_CK`) is a diagnostic threshold, NOT an achievable
fold-change value. Do not call it "ceiling" anywhere in the response or the UI that
consumes it — see §4 for why.

### 6.3 `POST /sweep`
**Request**
```json
{
  "base_params": { "k_ck": 10.0, "k_open_off": 0.01, "k_open_on": 1.0, "luckey": 500.0 },
  "sweep": { "param": "luckey", "min": 1, "max": 10000, "steps": 60, "scale": "log" }
}
```
**Response**
```json
{
  "param": "luckey",
  "scale": "log",
  "points": [ { "x": 1.0, "fold_change": 1.1, "dominance_ratio": 0.1 } ],
  "operating_point": { "x": 500.0, "fold_change": 38.2 }
}
```

### 6.4 `POST /suggest`
On-demand re-suggestion with a different policy / more variants.
```json
{ "sequence": "LISDAELEAIFAEELDC", "sensitive_window": {"start":1,"end":17},
  "substitution_policy": "neutralizing", "max_variants": 3 }
```
Returns the `suggested_variants` array shape from `/scan`.

### 6.5 Error envelope (all endpoints)
```json
{ "error": { "code": "INVALID_SEQUENCE", "field": "sequences[0].sequence",
             "message": "Contains non-standard residue 'X' at position 9." } }
```

---

## 7. Plot Specifications

Library: Recharts (or Plotly) in the frontend; matplotlib mirrors for CLI/figure export.

### 7.1 Fold-change vs [lucKey]  (primary)
- **X axis:** [lucKey], **log scale**, nM, range covering ≥3 decades around the operating point (e.g. 1–10,000 nM). Ticks at decades; minor ticks logged.
- **Y axis:** fold-change, linear (or log if it spans decades), starts at 1 (no change), not 0.
- **Series:** solid `--plot-line` fold-change curve.
- **Dominance-ratio reference:** dashed `--plot-ceiling` line/curve at `[lucKey]/K_CK` *for the current K_CK*, labeled **"lucKey/K_CK dominance ratio"** — never "ceiling." This is a diagnostic threshold the realistic curve sits near at typical pull values; it is not the curve's true mathematical asymptote (that's `(1+K_open+lucKey/K_CK)/K_open`, a much larger and not practically meaningful number — don't plot it).
- **Operating-point marker:** filled `--plot-marker` dot at the user's [lucKey], with a callout label showing the value. 
- **Zones (optional):** faint background shading where fold-change is poor (<2×, `--plot-zone-bad`) vs strong.
- **Tooltip:** on hover, "[lucKey] = X nM → fold-change Y×".

### 7.2 Fold-change vs K_open  (secondary — the teaching plot)
- **X axis:** K_open, log scale, dimensionless.
- **Y axis:** fold-change.
- Curve shows how fold-change responds to latch tuning. **In the key-limited regime this curve is nearly flat** — that flatness is the visual proof that latch mutations won't help, and the verdict text should point the user to it ("see how the curve barely moves").
- Operating-point marker at current K_open(OFF).

### 7.3 Optional 2D heatmap (advanced)
- Fold-change over (K_CK on log-x, [lucKey] on log-y), color = fold-change. Operating point marked. Helpful for seeing the whole design space at once; behind an "Advanced view" toggle.

### Axis/number rules
- Every concentration axis: **log scale by default**, auto unit scaling in tick labels (fM/pM/nM/µM) — never raw `0.000042`.
- Linear toggle available per plot.
- Colorblind-safe: don't rely on red/green alone in plots; the dominance-ratio reference is distinguished by dash pattern, zones by position + subtle hue.

---

## 8. Error & Edge States

| Situation | Behavior |
|-----------|----------|
| Invalid character typed in sequence | Rejected on input; field shakes; caption "only standard amino acids (ACDEFG…)" |
| Empty sequence on Scan | Scan button disabled; tooltip "enter a sequence first" |
| Sequence > 200 aa | Accept but warn: "long sequence — liability model tuned for peptide-scale binders" |
| Malformed FASTA (batch) | List which records failed with line numbers; scan the valid ones, flag skipped |
| Sensitive window outside length | Clamp handles to 1–N; can't invert (start ≤ end enforced) |
| K_CK / K_open / [lucKey] ≤ 0 or blank | Field border `--danger-700`, inline message; Calculate disabled |
| K_open(ON) ≤ K_open(OFF) | Non-blocking warning pill: "target should stabilize the open state (ON > OFF)" |
| Sweep min ≥ max or steps < 2 | Inline error on the sweep controls; plot not requested |
| Backend not reachable (`lockr serve` down) | Toast: "Can't reach the local engine — is `lockr serve` running?" with retry |
| Computation returns NaN/inf | Result card shows "Parameters produce an undefined result — check K_open and K_CK"; no plot |
| First load, no input yet | Results pane shows empty state (see §9.4) |
| Batch CSV export with 0 valid rows | Export disabled; tooltip explains |

Loading: computations are sub-second, but show a subtle inline spinner on the Calculate/Scan button while in flight; never block the whole screen.

---

## 9. Interface Copy (verbatim)

### 9.1 Headers & helpers
- Scanner H1: **"Charge & Liability Scanner"** · helper: "Check a binder sequence for charge liabilities that can weaken cage-key affinity (K_CK) before you synthesize or fold it."
- Calculator H1: **"Fold-Change Calculator"** · helper: "Predict your LOCKR sensor's dynamic range and find out what's limiting it."

### 9.2 Liability / K_CK penalty bands (Scanner)
- Low: **"Low liability."** "No significant charge liabilities in the sensitive region. Cage-key affinity should be preserved."
- Moderate: **"Moderate liability."** "Some acidic residues in sensitive positions. K_CK may be partially weakened — review the flagged residues."
- High / Severe: **"High liability."** "Multiple acidic residues in the sensitive region are likely to collapse cage-key affinity (K_CK). Strongly consider the suggested charge-optimized variant."

### 9.3 Regime verdicts (Calculator)
- **Key-limited:** banner title "Key-limited regime." Body: "Your fold-change is near the lucKey/K_CK dominance threshold (a diagnostic ratio, not an achievable target). Tuning the latch (K_open) won't raise it — see how the K_open curve stays flat. To improve dynamic range, raise [lucKey] or tighten cage-key affinity (lower K_CK)."
  - Recommendations: "Increase [lucKey] — this raises the dominance threshold." · "Improve cage-key affinity (lower K_CK)." · "Latch / K_open mutations will not help in this regime."
- **K_open-limited:** banner title "K_open-limited regime." Body: "You have headroom below the lucKey/K_CK dominance threshold. Tuning the latch to favor the open state on target binding will increase fold-change."
  - Recommendations: "Engineer the latch to raise K_open(ON) relative to K_open(OFF)." · "You're not key-limited yet — [lucKey] increases give diminishing returns."
- **Mixed:** banner title "Mixed regime." Body: "Both the latch (K_open) and the key (lucKey/K_CK) are constraining fold-change. Improvements to either will help; the plots show which gives more."

### 9.4 Empty states
- Scanner results pane: "Enter a binder sequence and run a scan to see its charge-liability profile and a suggested optimized variant."
- Calculator results pane: "Set your parameters and calculate to see predicted fold-change, the lucKey/K_CK dominance ratio, and what's limiting your sensor."

### 9.5 Key tooltips
- K_CK (`?`): "Cage-key dissociation constant — how tightly lucKey binds the open cage. Lower = tighter binding."
- K_open(OFF/ON): "The closed↔open balance of the latch. OFF = without target; ON = with target bound. Target binding should raise it."
- Substitution policy: "Conservative (D→N, E→Q) removes charge while keeping shape and H-bonding. Neutralizing (D→A, E→A) replaces with alanine."

---

## 10. Responsive Behavior

| Breakpoint | Layout |
|------------|--------|
| ≥ 1024px (desktop) | Two columns (input left, results right) as drawn |
| 640–1023px (tablet) | Single column, stacked: inputs first, results below; plots full width |
| < 640px (mobile) | Single column; Advanced panels collapsed by default; sequence ruler scrolls horizontally; plots switch to a compact aspect ratio; tab switcher becomes full-width |

- The annotated sequence and ruler always scroll horizontally rather than wrapping, so position alignment is never broken.
- Tables (batch mode) become horizontally scrollable with the ID column pinned.

---

## 11. Accessibility

- **Contrast:** all text meets WCAG AA on `--canvas`/`--surface`; brand-700 on white passes for large/medium text.
- **Color independence:** liability bands and regimes always pair color with a text label and an icon — never color alone (covers colorblind users; the danger/brand red overlap is fine because text labels disambiguate).
- **Keyboard:** full tab order; Scan/Calculate reachable and Enter-submittable; sweep range handles operable by arrow keys; tooltips focusable.
- **Screen readers:** inputs have `<label>`s and `aria-describedby` for help text; the verdict banner uses `role="status"` so the regime conclusion is announced after Calculate; plots have a text-summary fallback ("fold-change rises from 1× toward the lucKey/K_CK dominance ratio of 50 as [lucKey] increases").
- **Focus visible:** the `--focus-ring` is never removed.
- **Reduced motion:** respect `prefers-reduced-motion` — disable the input shake and plot transitions.

---

## 12. Component Inventory & State Model

### Components
```
<AppShell>                      header, tabs, footer, routing
  <ScannerTab>
    <SequenceInput>             textarea + live ruler + validation
    <BatchUploader>             FASTA drop + parse
    <SensitiveWindowSlider>
    <ScannerAdvanced>           pH, substitution policy
    <LiabilitySummaryCard>      gauge, net charge, K_CK band
    <AnnotatedSequence>         residue chips + ruler + tooltips
    <ContributionChart>
    <VariantDiff>               original vs suggested + chain button
    <BatchResultsTable>         sortable + CSV export
  <CalculatorTab>
    <ParamField>                numeric input w/ units, sci-notation, sweep toggle
    <CalculatorAdvanced>        target binding
    <VerdictCard>               hero FC, dominance ratio, regime banner
    <SweepPlot>                 log-axis line chart (reused for both plots)
    <Recommendations>
  <ChainPill>                   provenance badge on K_CK
  <Toast> <Tooltip> <Button> <SegmentedControl> <Badge>
```

### State
```ts
// per-tab local state
scanner: { mode, sequence|fastaRecords, ph, window:{start,end}, policy, results, status }
calculator: { kCk, kOpenOff, kOpenOn, lucKey, kTarget?, targetConc?,
              sweeps:{lucKey?, kOpen?}, result, status }

// shared (Context)
chain: { pipedKck: number|null, source: {label:string}|null }
ui:    { activeTab: 'scanner'|'calculator' }   // synced to URL hash
```

- No global server state, no persistence needed beyond the session. (Optional nicety: keep last inputs in memory so switching tabs doesn't clear them — do **not** use localStorage.)
- Every compute call: set `status='loading'` → POST → set `result` + `status='done'`, or `status='error'` + toast.

---

## 13. Code Style & Authorial Voice

The code should read like one person wrote it over a few weeks, not like it was generated. This matters for an iGEM software deliverable where judges read the repo — and it should just be *yours*. Conventions:

### Comments
- **Comment the *why*, not the *what*.** Skip `// loop over residues` above a loop that obviously loops over residues. Do write `// D/E here tank K_CK — see the v1 binder; that's the whole reason this scan exists` where the reasoning isn't obvious from the code.
- **Don't comment every function.** Real codebases have stretches of uncommented, self-explanatory code and then a dense comment where something is genuinely tricky or non-obvious. Uniform docstrings on every single function is a tell.
- **Let comments sound like you.** Short, occasionally informal, sometimes a little blunt ("this is gross but ProteinMPNN output is inconsistent so we parse defensively"). A reference to the actual project ("calibrated off the v1 vs optimized binder") reads human because it *is* specific to you.
- **Leave real TODOs and notes-to-self.** `# TODO: recalibrate once BLI data is back` or `# NOTE: K_open defaults are placeholders until I confirm the eq` are honest and human. Don't sprinkle fake ones, but don't scrub the real ones either.
- **No narration.** Avoid the AI tell of a comment restating the next line in English, or section-divider banners like `# ===== HELPER FUNCTIONS =====` on everything.

### Naming & structure
- Consistent but not robotic. Real code has a couple of slightly-too-short names (`fc`, `kck`) next to careful ones — that's fine and normal in a domain where those are the actual symbols.
- Don't over-abstract. Two similar functions are more human than one clever generic one with five flags. Refactor when it actually hurts, not preemptively.
- Some functions are long because the logic is linear; don't shatter everything into tiny helpers just for symmetry.

### What screams "AI-generated" — avoid
- Every function having an identical, exhaustive docstring with `Args:`/`Returns:`/`Raises:` when the project doesn't otherwise warrant it.
- Comments that explain language basics (`# increment counter by 1`).
- Perfectly uniform spacing/structure across every file with zero personality.
- Over-defensive try/except around things that can't fail, with a comment explaining the obvious.
- Emoji in comments, or cheerful filler ("Now let's create the awesome results card! ✨").

### Practical
- Write commits in your own voice too — `fix charge calc, termini were double-counted` not `Implement comprehensive charge calculation refactor`.
- It's fine (good, even) to leave a commented-out approach you tried with a note on why it didn't work.
- Keep a consistent personal style across the repo — pick tabs-or-spaces, quote style, etc. and stick to it, because *consistency within one author's voice* is itself what reads as human.

> Apply this to the engine, backend, and frontend alike — the Python, the FastAPI routes, and the React components.

---

## Build order (UI)
1. AppShell + tabs + design tokens (CSS variables from §1).
2. Calculator first — it's smaller and proves the engine round-trip (params → verdict → plot).
3. Scanner single-sequence (input, live annotation, summary, variant).
4. Scanner batch (uploader + table + export).
5. The chain (ChainPill + context wiring).
6. Responsive + a11y pass + empty/error states + copy.
