
let asmTandem = false;

function asmEl(id) { return document.getElementById(id); }

function asmReadFields() {
  const spacerSeq = asmEl("asm-spacer").value.trim().toUpperCase();
  const spacerStart = asmEl("asm-spacer-start").value ? parseInt(asmEl("asm-spacer-start").value) : null;
  const linkerSeq = asmTandem ? (asmEl("asm-linker").value.trim().toUpperCase() || null) : null;
  const linkerStart = asmTandem && asmEl("asm-linker-start").value ? parseInt(asmEl("asm-linker-start").value) : null;
  const binder2Seq = asmTandem ? (asmEl("asm-binder2").value.trim().toUpperCase() || null) : null;
  const binder2Start = asmTandem && asmEl("asm-binder2-start").value ? parseInt(asmEl("asm-binder2-start").value) : null;

  return {
    full_sequence: asmEl("asm-sequence").value.trim().toUpperCase(),
    latch_window: {
      start: parseInt(asmEl("asm-lw-start").value) || null,
      end: parseInt(asmEl("asm-lw-end").value) || null,
      expected_length: null,
    },
    graft_spec: {
      binder: asmEl("asm-binder").value.trim().toUpperCase(),
      start: parseInt(asmEl("asm-binder-start").value) || null,
      spacer: spacerSeq || null,
      spacer_start: spacerSeq ? spacerStart : null,
      linker: linkerSeq,
      linker_start: linkerSeq ? linkerStart : null,
      binder2: binder2Seq,
      binder2_start: binder2Seq ? binder2Start : null,
    },
    protected_region: {
      motif: asmEl("asm-pr-motif").value.trim().toUpperCase(),
      start: parseInt(asmEl("asm-pr-start").value) || null,
      end: parseInt(asmEl("asm-pr-end").value) || null,
      label: asmEl("asm-pr-label").value.trim(),
    },
    expected_total_length: asmEl("asm-total-length").value ? parseInt(asmEl("asm-total-length").value) : null,
    candidate_variants: [],
    binder_offset: 0,
  };
}

function asmValidate(f) {
  const errors = {};
  if (!f.full_sequence) errors.sequence = "required";
  if (!f.latch_window.start || f.latch_window.start < 1) errors.lw_start = "must be ≥ 1";
  if (!f.latch_window.end || f.latch_window.end < 1) errors.lw_end = "must be ≥ 1";
  if (f.latch_window.start && f.latch_window.end && f.latch_window.start > f.latch_window.end)
    errors.lw_end = "must be ≥ start";
  if (!f.graft_spec.binder) errors.binder = "required";
  if (!f.graft_spec.start || f.graft_spec.start < 1) errors.binder_start = "must be ≥ 1";
  if (f.graft_spec.spacer && !f.graft_spec.spacer_start) errors.spacer_start = "required when spacer is set";
  if (asmTandem && f.graft_spec.linker && !f.graft_spec.linker_start) errors.linker_start = "required";
  if (asmTandem && f.graft_spec.binder2 && !f.graft_spec.binder2_start) errors.binder2_start = "required";
  if (!f.protected_region.motif) errors.pr_motif = "required";
  if (!f.protected_region.start || f.protected_region.start < 1) errors.pr_start = "must be ≥ 1";
  if (!f.protected_region.end || f.protected_region.end < 1) errors.pr_end = "must be ≥ 1";
  if (f.protected_region.start && f.protected_region.end && f.protected_region.start > f.protected_region.end)
    errors.pr_end = "must be ≥ start";
  return errors;
}

const ASM_FIELD_MAP = {
  sequence: "asm-sequence",
  lw_start: "asm-lw-start",
  lw_end: "asm-lw-end",
  binder: "asm-binder",
  binder_start: "asm-binder-start",
  spacer_start: "asm-spacer-start",
  linker_start: "asm-linker-start",
  binder2_start: "asm-binder2-start",
  pr_motif: "asm-pr-motif",
  pr_start: "asm-pr-start",
  pr_end: "asm-pr-end",
};

function asmShowErrors(errors) {
  Object.keys(ASM_FIELD_MAP).forEach((key) => {
    const input = asmEl(ASM_FIELD_MAP[key]);
    const errEl = asmEl(`${ASM_FIELD_MAP[key]}-error`);
    const msg = errors[key];
    if (input) input.classList.toggle("invalid", Boolean(msg));
    if (errEl) errEl.textContent = msg || "";
  });
}

function asmUpdateValidity() {
  const fields = asmReadFields();
  const errors = asmValidate(fields);
  asmShowErrors(errors);
  asmEl("asm-submit").disabled = Object.keys(errors).length > 0;
}

const ASM_CHECK_LABELS = {
  overall_length: "Overall length",
  protected_region_intact: "Protected region intact",
  graft_no_overlap: "Graft does not overlap protected region",
  latch_fit: "Graft fits latch window",
  spacer_intact: "Spacer intact",
  binder1_intact: "Binder intact",
  linker_intact: "Linker intact",
  binder2_intact: "Binder 2 intact",
};

function asmRenderResults(data) {
  asmEl("asm-empty-state").style.display = "none";
  asmEl("asm-results").style.display = "block";

  const summaryEl = asmEl("asm-summary");
  summaryEl.textContent = data.all_passed ? "All checks passed." : "One or more checks failed.";
  summaryEl.className = data.all_passed ? "asm-summary-pass" : "asm-summary-fail";

  const listEl = asmEl("asm-check-list");
  listEl.innerHTML = "";
  data.checks.forEach((check) => {
    const row = document.createElement("div");
    row.className = "asm-check-row";

    const badge = document.createElement("span");
    badge.className = `badge ${check.passed ? "badge-pass" : "badge-fail"}`;
    badge.textContent = check.passed ? "Pass" : "Fail";

    const content = document.createElement("div");
    content.className = "asm-check-content";

    const nameEl = document.createElement("div");
    nameEl.className = "asm-check-name";
    nameEl.textContent = ASM_CHECK_LABELS[check.name] || check.name.replace(/_/g, " ");

    const detailEl = document.createElement("div");
    detailEl.className = "asm-check-detail";
    detailEl.textContent = check.detail;

    content.appendChild(nameEl);
    content.appendChild(detailEl);
    row.appendChild(badge);
    row.appendChild(content);
    listEl.appendChild(row);
  });
}

// LINE_LEN residues per row in the sequence map
const ASM_SEQ_LINE = 60;

function asmRenderSequence(fields) {
  const seq = fields.full_sequence;
  const n = seq.length;
  const gs = fields.graft_spec;
  const pr = fields.protected_region;
  const lw = fields.latch_window;

  // per-residue color, 1-indexed; later assignments override earlier (binder beats latch)
  const color = new Array(n + 1).fill("");
  if (lw.start && lw.end) {
    for (let i = lw.start; i <= Math.min(lw.end, n); i++) color[i] = "latch";
  }
  if (pr.start && pr.end) {
    for (let i = pr.start; i <= Math.min(pr.end, n); i++) color[i] = "protected";
  }
  if (gs.spacer && gs.spacer_start) {
    for (let i = gs.spacer_start; i < gs.spacer_start + gs.spacer.length && i <= n; i++) color[i] = "spacer";
  }
  if (gs.binder && gs.start) {
    for (let i = gs.start; i < gs.start + gs.binder.length && i <= n; i++) color[i] = "binder";
  }
  if (gs.binder2 && gs.binder2_start) {
    for (let i = gs.binder2_start; i < gs.binder2_start + gs.binder2.length && i <= n; i++) color[i] = "binder";
  }

  const container = asmEl("asm-seq-display");
  container.innerHTML = "";
  for (let lineStart = 0; lineStart < n; lineStart += ASM_SEQ_LINE) {
    const lineDiv = document.createElement("div");
    lineDiv.className = "asm-seq-line";

    const num = document.createElement("span");
    num.className = "asm-seq-linenum";
    num.textContent = lineStart + 1;
    lineDiv.appendChild(num);

    for (let j = lineStart; j < Math.min(lineStart + ASM_SEQ_LINE, n); j++) {
      const span = document.createElement("span");
      const c = color[j + 1];
      if (c) span.className = `asm-hl asm-hl--${c}`;
      span.textContent = seq[j];
      lineDiv.appendChild(span);
    }
    container.appendChild(lineDiv);
  }
}

async function asmSubmit() {
  const fields = asmReadFields();
  const errors = asmValidate(fields);
  if (Object.keys(errors).length > 0) return;

  const button = asmEl("asm-submit");
  const originalLabel = button.innerHTML;
  button.disabled = true;
  button.innerHTML = '<span class="spinner"></span>Verifying…';

  try {
    const result = await apiPost("/verify-assembly", fields);
    asmRenderResults(result);
    asmRenderSequence(fields);
  } catch (err) {
    showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function asmReset() {
  asmEl("asm-form").reset();
  asmTandem = false;
  asmEl("asm-tandem-cb").checked = false;
  asmEl("asm-tandem-fields").style.display = "none";
  asmEl("asm-seq-display").innerHTML = "";
  asmEl("asm-empty-state").style.display = "block";
  asmEl("asm-results").style.display = "none";
  asmLastKopenData = null;
  asmEl("asm-kopen-results").style.display = "none";
  asmEl("asm-kopen-error").textContent = "";
  asmEl("asm-kopen-preserve").value = "";
  const calcKopen = document.getElementById("calc-kopen");
  asmEl("asm-kopen-current").value = (calcKopen && calcKopen.value) || "0.001";
  asmUpdateValidity();
}

// --- K_open Tuning Analysis -------------------------------------------

let asmDestabPreset = "conservative";
let asmLastKopenData = null;

// I/V residues in the latch window, downstream of the binder(s), not in
// preserve_positions -- the same logic used to find I346/I356 by hand.
function asmFindKopenCandidates(fields, preservePositions) {
  const seq = fields.full_sequence;
  const lw = fields.latch_window;
  const gs = fields.graft_spec;
  if (!seq || !lw.start || !lw.end) return [];

  let binderEnd = gs.start ? gs.start + (gs.binder || "").length - 1 : 0;
  if (gs.binder2_start) binderEnd = Math.max(binderEnd, gs.binder2_start + gs.binder2.length - 1);

  const candidates = [];
  for (let pos = lw.start; pos <= Math.min(lw.end, seq.length); pos++) {
    const residue = seq[pos - 1];
    if (residue !== "I" && residue !== "V") continue;
    if (pos <= binderEnd) continue;
    if (preservePositions.includes(pos)) continue;
    candidates.push({ position: pos, residue });
  }
  return candidates;
}

function asmParsePositions(text) {
  return text.split(",").map((s) => parseInt(s.trim(), 10)).filter((n) => !isNaN(n));
}

function asmRenderKopenCandidates(candidates) {
  const el = asmEl("asm-kopen-candidates");
  if (candidates.length === 0) {
    el.textContent = "No candidate I/V positions found in the native latch region downstream of the binder — latch tuning may require custom design.";
  } else {
    const labels = candidates.map((c) => `${c.residue}${c.position}`).join(", ");
    el.textContent = `Candidate positions: ${labels} (downstream of binder, safe to mutate)`;
  }
}

const ASM_REGIME_LABELS = { key_limited: "Key-limited", K_open_limited: "K_open-limited", mixed: "Mixed" };
const ASM_HELPS_BADGE = { yes: ["Yes", "badge-Low"], no: ["No", "badge-muted"],
                         slightly_worse: ["Slightly worse", "badge-Moderate"] };

function asmRenderKopenTable(data) {
  const rows = [data.baseline, ...data.scenarios.filter((s) => s.preset === asmDestabPreset)];
  const table = asmEl("asm-kopen-table");
  table.innerHTML = "";

  const thead = document.createElement("thead");
  thead.innerHTML = "<tr><th>Mutations</th><th>K_open</th><th>Fold-change</th><th>Regime</th><th>Helps?</th></tr>";
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  rows.forEach((row) => {
    const tr = document.createElement("tr");

    const mutTd = document.createElement("td");
    mutTd.textContent = row.mutations === 0 ? "Baseline" : `+${row.mutations}`;
    tr.appendChild(mutTd);

    const kopenTd = document.createElement("td");
    kopenTd.className = "font-mono";
    kopenTd.textContent = roundSig(row.k_open, 3);
    tr.appendChild(kopenTd);

    const fcTd = document.createElement("td");
    fcTd.className = "font-mono";
    fcTd.textContent = formatFoldChange(row.fold_change);
    tr.appendChild(fcTd);

    const regimeTd = document.createElement("td");
    regimeTd.textContent = ASM_REGIME_LABELS[row.regime] || row.regime;
    tr.appendChild(regimeTd);

    const helpsTd = document.createElement("td");
    if (row.helps === "unchanged") {
      helpsTd.textContent = "—";
    } else {
      const [label, cls] = ASM_HELPS_BADGE[row.helps] || [row.helps, "badge-muted"];
      const badge = document.createElement("span");
      badge.className = `badge ${cls}`;
      badge.textContent = label;
      helpsTd.appendChild(badge);
    }
    tr.appendChild(helpsTd);

    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
}

function asmRenderKopenNotes(data, currentLuckey) {
  const mutationRows = data.scenarios.filter((s) => s.preset === asmDestabPreset);
  const anyHelps = mutationRows.some((s) => s.helps === "yes");
  const negativeEl = asmEl("asm-kopen-negative-note");
  if (!anyHelps) {
    negativeEl.style.display = "block";
    negativeEl.textContent = `At your current lucKey/K_CK ratio (${roundSig(data.dominance_ratio, 3)}), K_open changes have negligible effect. Latch tuning will only help if lucKey/K_CK is reduced to be comparable to K_open. Try lowering lucKey to ${roundSig(data.crossover_luckey_nm, 3)} nM.`;
  } else {
    negativeEl.style.display = "none";
  }

  asmEl("asm-kopen-crossover-note").textContent =
    `Latch tuning becomes meaningful when lucKey/K_CK ≈ K_open. At your current settings: lucKey/K_CK = ${roundSig(data.dominance_ratio, 3)}, K_open = ${roundSig(data.baseline.k_open, 3)}. Tuning helps when [lucKey] ≈ ${roundSig(data.crossover_luckey_nm, 3)} nM (currently ${roundSig(currentLuckey, 3)} nM).`;
}

// Falls back to the Calculator tab's fields (which ship pre-filled with
// ECLIPSE-example defaults) so this works even if the user never touched
// the Calculator tab.
function asmReadChainedParams() {
  const readNum = (id, fallback) => {
    const el = document.getElementById(id);
    const v = el ? parseFloat(el.value) : NaN;
    return isNaN(v) ? fallback : v;
  };
  return {
    k_ck: readNum("calc-kck", 10),
    luckey: readNum("calc-luckey", 500),
    pull: readNum("calc-pull", 10),
  };
}

async function asmKopenAnalyze() {
  const fields = asmReadFields();
  const errEl = asmEl("asm-kopen-error");
  errEl.textContent = "";

  if (!fields.full_sequence || !fields.latch_window.start || !fields.latch_window.end) {
    errEl.textContent = "Fill in the full sequence and graftable latch window first.";
    return;
  }

  const kopenCurrent = parseFloat(asmEl("asm-kopen-current").value);
  if (!(kopenCurrent > 0)) {
    errEl.textContent = "Current K_open must be > 0.";
    return;
  }

  const preservePositions = asmParsePositions(asmEl("asm-kopen-preserve").value);
  const candidates = asmFindKopenCandidates(fields, preservePositions);
  asmRenderKopenCandidates(candidates);

  const chained = asmReadChainedParams();
  const button = asmEl("asm-kopen-analyze");
  const originalLabel = button.innerHTML;
  button.disabled = true;
  button.innerHTML = '<span class="spinner"></span>Analyzing…';

  try {
    const data = await apiPost("/kopen-scenarios", {
      k_open_current: kopenCurrent, k_ck: chained.k_ck, luckey: chained.luckey, pull: chained.pull,
    });
    asmLastKopenData = data;
    asmRenderKopenTable(data);
    asmRenderKopenNotes(data, chained.luckey);
    asmEl("asm-kopen-results").style.display = "block";
  } catch (err) {
    showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function initAssembly() {
  // auto-uppercase the full_sequence textarea
  asmEl("asm-sequence").addEventListener("input", () => {
    const el = asmEl("asm-sequence");
    const pos = el.selectionStart;
    el.value = el.value.toUpperCase();
    el.setSelectionRange(pos, pos);
    asmUpdateValidity();
  });

  asmEl("asm-tandem-cb").addEventListener("change", (e) => {
    asmTandem = e.target.checked;
    asmEl("asm-tandem-fields").style.display = asmTandem ? "block" : "none";
    asmUpdateValidity();
  });

  ["asm-lw-start", "asm-lw-end",
   "asm-binder", "asm-binder-start", "asm-spacer", "asm-spacer-start",
   "asm-linker", "asm-linker-start", "asm-binder2", "asm-binder2-start",
   "asm-pr-motif", "asm-pr-start", "asm-pr-end", "asm-pr-label",
   "asm-total-length"].forEach((id) => {
    const el = asmEl(id);
    if (el) el.addEventListener("input", asmUpdateValidity);
  });

  asmEl("asm-submit").addEventListener("click", asmSubmit);
  asmEl("asm-reset").addEventListener("click", asmReset);
  asmUpdateValidity();

  const calcKopen = document.getElementById("calc-kopen");
  if (calcKopen && calcKopen.value) asmEl("asm-kopen-current").value = calcKopen.value;

  document.querySelectorAll("#asm-destab-segmented button").forEach((btn) => {
    btn.addEventListener("click", () => {
      asmDestabPreset = btn.dataset.destab;
      document.querySelectorAll("#asm-destab-segmented button").forEach((b) => b.classList.toggle("active", b === btn));
      if (asmLastKopenData) {
        asmRenderKopenTable(asmLastKopenData);
        asmRenderKopenNotes(asmLastKopenData, asmReadChainedParams().luckey);
      }
    });
  });

  asmEl("asm-kopen-analyze").addEventListener("click", asmKopenAnalyze);
}

document.addEventListener("DOMContentLoaded", initAssembly);
