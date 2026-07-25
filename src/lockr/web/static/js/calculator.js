
const calcState = { result: null, values: null, sweepLuckey: null, sweepKopen: null };

// "kd" or "dg", which mode the K_target entry is in
let calcKtargetMode = "kd";
let calcKdLastValue = ""; // saved when toggling to dG mode

function calcEl(id) { return document.getElementById(id); }

// compute Kd (nM) from dG fields; returns null if inputs are missing/invalid
function calcDgComputeKd() {
  const dgStr = calcEl("calc-dg-score").value.trim();
  if (!dgStr) return null;
  const dg = parseFloat(dgStr);
  if (!isFinite(dg)) return null;

  const refDgStr = calcEl("calc-dg-ref-dg").value.trim();
  const refKdStr = calcEl("calc-dg-ref-kd").value.trim();
  if (refDgStr && refKdStr) {
    const refDg = parseFloat(refDgStr);
    const refKd = parseFloat(refKdStr);
    if (isFinite(refDg) && isFinite(refKd) && refKd > 0) {
      // anchored ddG form: Kd = ref_Kd * exp((dG - ref_dG) / RT)
      return refKd * Math.exp((dg - refDg) / RT_KCAL_MOL);
    }
  }
  // unanchored: rough estimate only
  return Math.exp(dg / RT_KCAL_MOL) * 1e9;
}

function calcDgUpdateDerived() {
  if (calcKtargetMode !== "dg") return;
  const dgStr = calcEl("calc-dg-score").value.trim();
  const derivedRow = calcEl("calc-dg-derived-row");
  const derivedEl = calcEl("calc-dg-derived");
  const warnEl = calcEl("calc-dg-no-anchor-warn");

  if (!dgStr) { derivedRow.style.display = "none"; return; }

  const kd = calcDgComputeKd();
  if (kd === null) { derivedRow.style.display = "none"; return; }

  derivedRow.style.display = "block";
  derivedEl.textContent = `Computed K_target = ${formatConcNm(kd)} (${roundSig(kd, 4)} nM)`;

  const refDgStr = calcEl("calc-dg-ref-dg").value.trim();
  const refKdStr = calcEl("calc-dg-ref-kd").value.trim();
  const hasAnchor = refDgStr && refKdStr &&
                    isFinite(parseFloat(refDgStr)) && parseFloat(refKdStr) > 0;
  warnEl.style.display = hasAnchor ? "none" : "block";
}

function calcKtargetSetMode(mode) {
  calcKtargetMode = mode;
  document.querySelectorAll("#calc-ktarget-mode-seg button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.kmode === mode);
  });
  if (mode === "kd") {
    calcEl("calc-ktarget-kd-panel").style.display = "block";
    calcEl("calc-ktarget-dg-panel").style.display = "none";
    if (calcKdLastValue !== "") calcEl("calc-ktarget").value = calcKdLastValue;
  } else {
    calcKdLastValue = calcEl("calc-ktarget").value;
    calcEl("calc-ktarget-kd-panel").style.display = "none";
    calcEl("calc-ktarget-dg-panel").style.display = "block";
    calcDgUpdateDerived();
  }
  calcUpdateValidity();
}

function calcReadFields() {
  let kTarget;
  if (calcKtargetMode === "dg") {
    kTarget = calcDgComputeKd(); // null if empty, positive number otherwise
  } else {
    kTarget = calcEl("calc-ktarget").value === "" ? null : parseFloat(calcEl("calc-ktarget").value);
  }
  return {
    k_ck: parseFloat(calcEl("calc-kck").value),
    k_open: parseFloat(calcEl("calc-kopen").value),
    pull: parseFloat(calcEl("calc-pull").value),
    luckey: parseFloat(calcEl("calc-luckey").value),
    k_target: kTarget,
    target_conc: calcEl("calc-targetconc").value === "" ? null : parseFloat(calcEl("calc-targetconc").value),
  };
}

// mirrors backend validation so Calculate disables before a bad request
function calcValidate(values) {
  const errors = {};
  ["k_ck", "k_open", "luckey"].forEach((field) => {
    if (!(values[field] > 0)) errors[field] = "must be > 0";
  });
  if (!(values.pull >= 0)) errors.pull = "must be >= 0";
  if ((values.k_target === null) !== (values.target_conc === null)) {
    errors.k_target = "K_target and target concentration must both be set, or both left blank";
    errors.target_conc = errors.k_target;
  } else if (values.k_target !== null && !(values.k_target > 0)) {
    errors.k_target = "must be > 0";
  }
  if (values.target_conc !== null && values.target_conc < 0) {
    errors.target_conc = "must be ≥ 0";
  }
  return errors;
}

function calcShowFieldErrors(errors) {
  ["k_ck", "k_open", "pull", "luckey", "k_target", "target_conc"].forEach((field) => {
    const errEl = document.querySelector(`[data-calc-error="${field}"]`);
    const msg = errors[field];
    if (field === "k_target" && calcKtargetMode === "dg") {
      // route k_target errors to the dG input, not the hidden Kd input
      const dgInput = calcEl("calc-dg-score");
      const dgErrEl = calcEl("calc-dg-score-error");
      if (dgInput) dgInput.classList.toggle("invalid", Boolean(msg));
      if (dgErrEl) dgErrEl.textContent = msg || "";
      // clear hidden Kd input styling
      const kdInput = calcEl("calc-ktarget");
      if (kdInput) kdInput.classList.remove("invalid");
      if (errEl) errEl.textContent = "";
    } else {
      const input = document.querySelector(`[data-calc-field="${field}"]`);
      if (input) input.classList.toggle("invalid", Boolean(msg));
      if (errEl) errEl.textContent = msg || "";
    }
  });
}

// pull > 50 has no significance in our documented ECLIPSE runs
function calcCheckPullWarning(values) {
  const warnEl = calcEl("calc-pull-warning");
  warnEl.style.display = values.pull > 50 ? "block" : "none";
}

function calcUpdateDerivedKopenEff(values) {
  const el = calcEl("calc-kopen-eff-derived");
  if (values.pull >= 0 && isFinite(values.pull)) {
    el.textContent = `With target bound, the cage opens ~${roundSig(1 + values.pull, 3)}× more readily`;
  } else {
    el.textContent = "With target bound, the cage opens ~, × more readily";
  }
}

function calcUpdateValidity() {
  const values = calcReadFields();
  const errors = calcValidate(values);
  calcShowFieldErrors(errors);
  calcCheckPullWarning(values);
  calcUpdateDerivedKopenEff(values);
  calcDgUpdateDerived();
  calcEl("calc-submit").disabled = Object.keys(errors).length > 0;
  return { values, errors };
}

function regimeToCssClass(regime) {
  return regime.toLowerCase().replace(/_/g, "-");
}

const REGIME_TITLES = {
  key_limited: "The key side is your limit.",
  K_open_limited: "The cage is your limit.",
  mixed: "Both sides are limiting.",
};

// values.k_target === null means the request assumed a saturating target,
// which is a different reason for a missing LOD than "pull too low to get there".
function calcRenderLod(result, values) {
  const atSaturation = values.k_target === null;

  const lod2El = calcEl("calc-metric-lod2");
  const lod3El = calcEl("calc-metric-lod3");
  const ec50El = calcEl("calc-metric-ec50");
  const noteEl = calcEl("calc-lod-note");

  if (atSaturation) {
    lod2El.textContent = "N/A";
    lod3El.textContent = "N/A";
    ec50El.textContent = "N/A";
    noteEl.style.display = "block";
    noteEl.textContent = "N/A (target assumed saturating), enter a K_target and target concentration to get a detection limit.";
  } else {
    lod2El.textContent = result.lod_2x_nm === null ? "not achievable at current pull" : formatConcNm(result.lod_2x_nm);
    lod3El.textContent = result.lod_3x_nm === null ? "not achievable at current pull" : formatConcNm(result.lod_3x_nm);
    ec50El.textContent = result.ec50_nm === null ? "N/A" : formatConcNm(result.ec50_nm);
    noteEl.style.display = "none";
  }

  const thresholdStr = calcEl("calc-clinical-threshold").value.trim();
  const threshold = thresholdStr === "" ? null : parseFloat(thresholdStr);
  const beatsEl = calcEl("calc-metric-beats");
  const thresholdEl = calcEl("calc-metric-threshold");

  thresholdEl.textContent = threshold === null ? "N/A" : formatConcNm(threshold);
  if (threshold === null) {
    beatsEl.textContent = "N/A";
    beatsEl.style.color = "";
  } else if (result.lod_2x_nm !== null && !atSaturation && result.lod_2x_nm <= threshold) {
    beatsEl.textContent = "YES ✓";
    beatsEl.style.color = "var(--success-700)";
  } else {
    beatsEl.textContent = "NO ✗";
    beatsEl.style.color = "var(--danger-700)";
  }
}

function calcRenderVerdict(result, values) {
  calcEl("calc-empty-state").style.display = "none";
  calcEl("calc-results").style.display = "block";

  calcEl("calc-hero-fc").textContent = formatFoldChange(result.fold_change);

  const qualityEl = calcEl("calc-quality");
  qualityEl.textContent = result.quality;
  qualityEl.className = `result-quality ${result.quality}`;

  calcEl("calc-interpretation").textContent = result.interpretation;
  calcEl("calc-metric-max").textContent = formatFoldChange(result.max_fold_change);
  calcEl("calc-metric-ratio").textContent = roundSig(result.dominance_ratio, 3);
  calcEl("calc-metric-fraction").textContent = `${roundSig(result.fraction_of_dominance_ratio * 100, 3)}%`;
  calcRenderLod(result, values);

  const banner = calcEl("calc-regime-banner");
  banner.className = `verdict-banner ${regimeToCssClass(result.regime)}`;
  calcEl("calc-regime-title").textContent = REGIME_TITLES[result.regime] || "Regime.";
  calcEl("calc-regime-body").textContent = result.verdict;

  const list = calcEl("calc-recommendations");
  list.innerHTML = "";
  result.recommendations.forEach((rec) => {
    const li = document.createElement("li");
    li.textContent = rec;
    list.appendChild(li);
  });

  const warningsEl = calcEl("calc-result-warnings");
  warningsEl.innerHTML = "";
  (result.warnings || []).forEach((w) => {
    const div = document.createElement("div");
    div.className = "field-error";
    div.textContent = w;
    warningsEl.appendChild(div);
  });
}

function autoSweepRange(value) {
  return { min: value / 1000, max: value * 1000, steps: 60, scale: "log" };
}

async function calcRunSweeps(values) {
  const base_params = { k_ck: values.k_ck, k_open: values.k_open, pull: values.pull, luckey: values.luckey };

  // Both sweeps are independent, fire them together instead of one-after-another.
  const [luckeySweep, kopenSweep] = await Promise.all([
    apiPost("/sweep", { base_params, sweep: { param: "luckey", ...autoSweepRange(values.luckey) } }),
    apiPost("/sweep", { base_params, sweep: { param: "k_open", ...autoSweepRange(values.k_open) } }),
  ]);

  drawLogXPlot(calcEl("calc-plot-luckey"), {
    xIsConcentration: true,
    lines: [
      { points: luckeySweep.points.map((p) => ({ x: p.x, y: p.fold_change })), color: "var(--plot-line)" },
      { points: luckeySweep.points.map((p) => ({ x: p.x, y: p.dominance_ratio })), color: "var(--plot-ceiling)", dashed: true },
    ],
    marker: { x: luckeySweep.operating_point.x, y: luckeySweep.operating_point.fold_change, label: formatConcNm(luckeySweep.operating_point.x) },
  });

  drawLogXPlot(calcEl("calc-plot-kopen"), {
    xIsConcentration: false,
    lines: [{ points: kopenSweep.points.map((p) => ({ x: p.x, y: p.fold_change })), color: "var(--plot-line)" }],
    marker: { x: kopenSweep.operating_point.x, y: kopenSweep.operating_point.fold_change, label: roundSig(kopenSweep.operating_point.x, 3) },
  });
}

async function calcSubmit() {
  const { values, errors } = calcUpdateValidity();
  if (Object.keys(errors).length > 0) return;

  const button = calcEl("calc-submit");
  const originalLabel = button.innerHTML;
  button.disabled = true;
  button.innerHTML = '<span class="spinner"></span>Calculating…';

  try {
    const result = await apiPost("/foldchange", values);
    calcEl("calc-undefined-result").style.display = "none";
    calcState.result = result;
    calcState.values = values;
    calcRenderVerdict(result, values);
    await calcRunSweeps(values);
  } catch (err) {
    if (err.code === "UNDEFINED_RESULT") {
      calcEl("calc-empty-state").style.display = "none";
      calcEl("calc-results").style.display = "none";
      calcEl("calc-undefined-result").style.display = "block";
    } else {
      showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
    }
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function calcReset() {
  calcEl("calc-form").reset();
  calcState.result = null;
  calcState.values = null;
  calcDetachPill();
  calcKdLastValue = "";
  calcKtargetMode = "kd";
  document.querySelectorAll("#calc-ktarget-mode-seg button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.kmode === "kd");
  });
  calcEl("calc-ktarget-kd-panel").style.display = "block";
  calcEl("calc-ktarget-dg-panel").style.display = "none";
  calcEl("calc-dg-derived-row").style.display = "none";
  calcUpdateValidity();
  calcEl("calc-empty-state").style.display = "block";
  calcEl("calc-results").style.display = "none";
  calcEl("calc-undefined-result").style.display = "none";

  calcEl("calc-kopen-latch").value = "";
  calcEl("calc-kopen-preserve").value = "";
  calcEl("calc-kopen-error").textContent = "";
  calcEl("calc-kopen-results").style.display = "none";

  calcKopenDetachPill();
  calcKopenClosePanel();
  calcKopenEstDestab = "moderate";
  document.querySelectorAll("#calc-kopen-est-destab-segmented button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.destab === "moderate");
  });
  calcKopenEstUpdate();
  calcLastKopenData = null;
  calcDestabPreset = "conservative";
  document.querySelectorAll("#calc-destab-segmented button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.destab === "conservative");
  });
}

function calcDetachPill() {
  calcEl("calc-kck-pill").style.display = "none";
  window.lockrChain.sourceLabel = null;
}

function calcApplyPiped() {
  if (window.lockrChain.pipedKck === null) return;
  calcEl("calc-kck").value = window.lockrChain.pipedKck;
  const pill = calcEl("calc-kck-pill");
  const label = window.lockrChain.sourceLabel || "Scanner";
  pill.innerHTML = `from Scanner · ${label} <button type="button" aria-label="Detach">×</button>`;
  pill.style.display = "inline-flex";
  pill.querySelector("button").addEventListener("click", calcDetachPill);
  calcUpdateValidity();
}

// --- K_open estimator ------------------------------------------------------
// Three dry-lab paths to a K_open value: base scaffold, destabilization
// arithmetic (same formula as thermo.k_open_from_destab), or a literature
// lookup. All client-side, no API call needed for any of the three.

const CALC_KOPEN_DESTAB_PRESETS = { conservative: 0.5, moderate: 1.0, optimistic: 1.5 };
const CALC_KOPEN_LENGTH_KCAL_PER_RESIDUE = 0.15;
let calcKopenEstDestab = "moderate";

// Real reported/derived values only, see Langan 2019 Fig 1b / Extended
// Data Fig 1 and Quijano-Rubio 2021. Not inventing rows past what's published.
const CALC_KOPEN_PATH3_ROWS = [
  { design: "Base lucCage (unmodified)", kopen: 1e-3, display: "1×10⁻³", source: "Quijano-Rubio 2021" },
  { design: "V223S/I238S (2 serine muts)", kopen: 5e-3, display: "~5×10⁻³", source: "Langan 2019" },
  { design: "Extended latch (+5 residues)", kopen: 2e-4, display: "~2×10⁻⁴", source: "Langan 2019" },
  { design: "Extended latch (+9 residues)", kopen: 5e-5, display: "~5×10⁻⁵", source: "Langan 2019" },
  { design: "Extended latch (+18 residues)", kopen: 1e-6, display: "~1×10⁻⁶", source: "Langan 2019" },
  { design: "Truncated latch (-9 residues)", kopen: 1e-2, display: "~1×10⁻²", source: "Langan 2019" },
  { design: "Optimized sCage (shortened)", kopen: 1e-2, display: "~1×10⁻²", source: "Quijano-Rubio 2021" },
];

// e.g. 0.0293 -> "2.9×10⁻²"
function formatKopenSci(x, sig = 2) {
  if (!(x > 0) || !isFinite(x)) return "N/A";
  const exp = Math.floor(Math.log10(x));
  const mantissa = roundSig(x / Math.pow(10, exp), sig);
  const superscripts = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
  const expStr = String(exp).split("").map((c) => superscripts[c] || c).join("");
  return `${mantissa}×10${expStr}`;
}

function calcKopenDetachPill() {
  calcEl("calc-kopen-pill").style.display = "none";
}

function calcKopenSetPill(sourceLabel) {
  const pill = calcEl("calc-kopen-pill");
  pill.innerHTML = `from ${sourceLabel} <button type="button" aria-label="Detach">×</button>`;
  pill.style.display = "inline-flex";
  pill.querySelector("button").addEventListener("click", calcKopenDetachPill);
}

function calcKopenClosePanel() {
  calcEl("calc-kopen-estimator").style.display = "none";
}

function calcKopenApplyValue(value, sourceLabel) {
  calcEl("calc-kopen").value = value;
  if (sourceLabel) calcKopenSetPill(sourceLabel);
  else calcKopenDetachPill();
  calcKopenClosePanel();
  calcUpdateValidity();
}

// The base latch this tool assumes is 35 residues (see the field's own help
// text), so a length change bigger than that isn't a real edit to this latch.
const CALC_KOPEN_BASE_LATCH_LENGTH = 35;

function calcKopenEstUpdate() {
  const start = parseFloat(calcEl("calc-kopen-est-start").value);
  const nMut = parseFloat(calcEl("calc-kopen-est-nmut").value);
  const lengthChange = parseFloat(calcEl("calc-kopen-est-length").value);
  const resultEl = calcEl("calc-kopen-est-result");
  const warnEl = calcEl("calc-kopen-est-warning");

  if (!(start > 0) || isNaN(nMut) || isNaN(lengthChange)) {
    resultEl.textContent = "Estimated K_open = N/A";
    warnEl.style.display = "none";
    return null;
  }

  if (nMut < 0) {
    resultEl.textContent = "Estimated K_open = N/A";
    warnEl.textContent = "Number of mutations can't be negative.";
    warnEl.style.display = "block";
    return null;
  }

  if (Math.abs(lengthChange) > CALC_KOPEN_BASE_LATCH_LENGTH) {
    warnEl.textContent = `A length change this large removes or adds more residues than the ${CALC_KOPEN_BASE_LATCH_LENGTH}-residue base latch has. The math below still runs, but this isn't a realistic edit to this latch.`;
    warnEl.style.display = "block";
  } else {
    warnEl.style.display = "none";
  }

  const destabPerMut = CALC_KOPEN_DESTAB_PRESETS[calcKopenEstDestab];
  const totalDestab = nMut * destabPerMut - lengthChange * CALC_KOPEN_LENGTH_KCAL_PER_RESIDUE;
  const dgOpenBase = -RT_KCAL_MOL * Math.log(start);
  const dgOpenNew = dgOpenBase - totalDestab;
  const kOpenEstimated = Math.exp(-dgOpenNew / RT_KCAL_MOL);

  resultEl.textContent = `Estimated K_open = ${formatKopenSci(kOpenEstimated)}`;
  return kOpenEstimated;
}

function calcKopenRenderPath3Table() {
  const table = calcEl("calc-kopen-path3-table");
  table.innerHTML = "";

  const thead = document.createElement("thead");
  thead.innerHTML = "<tr><th>Design</th><th>K_open</th><th>Source</th><th></th></tr>";
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  CALC_KOPEN_PATH3_ROWS.forEach((row) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${row.design}</td><td class="font-mono">${row.display}</td><td>${row.source}</td>`;
    const useTd = document.createElement("td");
    const useBtn = document.createElement("button");
    useBtn.type = "button";
    useBtn.className = "btn-secondary";
    useBtn.style = "height:28px; padding:0 10px; font-size:13px;";
    useBtn.textContent = "Use";
    useBtn.addEventListener("click", () => calcKopenApplyValue(row.kopen, row.source));
    useTd.appendChild(useBtn);
    tr.appendChild(useTd);
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
}

function initCalcKopenEstimator() {
  calcKopenRenderPath3Table();

  const toggle = calcEl("calc-kopen-estimator-toggle");
  const panel = calcEl("calc-kopen-estimator");
  toggle.addEventListener("click", (e) => {
    e.stopPropagation();
    panel.style.display = panel.style.display === "none" ? "block" : "none";
  });
  panel.addEventListener("click", (e) => e.stopPropagation());
  document.addEventListener("click", () => { panel.style.display = "none"; });

  // typing K_open directly detaches the pill and closes the panel, same as K_CK
  calcEl("calc-kopen").addEventListener("input", () => {
    calcKopenDetachPill();
    calcKopenClosePanel();
  });

  calcEl("calc-kopen-path1-use").addEventListener("click", () => calcKopenApplyValue(0.001, null));

  document.querySelectorAll("#calc-kopen-est-destab-segmented button").forEach((btn) => {
    btn.addEventListener("click", () => {
      calcKopenEstDestab = btn.dataset.destab;
      document.querySelectorAll("#calc-kopen-est-destab-segmented button").forEach((b) => b.classList.toggle("active", b === btn));
      calcKopenEstUpdate();
    });
  });
  ["calc-kopen-est-start", "calc-kopen-est-nmut", "calc-kopen-est-length"].forEach((id) => {
    calcEl(id).addEventListener("input", calcKopenEstUpdate);
  });
  calcEl("calc-kopen-path2-use").addEventListener("click", () => {
    const value = calcKopenEstUpdate();
    if (value !== null) calcKopenApplyValue(value, null);
  });
}

// Starting points so researchers don't face four blank fields. Values are in
// the same units as the form (K_CK/lucKey in nM, K_open dimensionless).
const CALC_PRESETS = {
  eclipse:    { k_ck: 10,  k_open: 0.001, pull: 10, luckey: 500 },
  "tight-key": { k_ck: 1,   k_open: 0.01,  pull: 15, luckey: 1000 },
  leaky:      { k_ck: 10,  k_open: 0.1,   pull: 10, luckey: 500 },
  blank:      { k_ck: "",  k_open: "",    pull: "", luckey: "" },
};

function calcApplyPreset(key) {
  const preset = CALC_PRESETS[key];
  if (!preset) return;
  calcDetachPill();
  calcKopenDetachPill();
  calcKopenClosePanel();
  calcEl("calc-kck").value = preset.k_ck;
  calcEl("calc-kopen").value = preset.k_open;
  calcEl("calc-pull").value = preset.pull;
  calcEl("calc-luckey").value = preset.luckey;
  // Presets describe the sensor, not a specific sample, clear the target section.
  calcEl("calc-ktarget").value = "";
  calcEl("calc-targetconc").value = "";
  calcUpdateValidity();
}

// --- K_open Tuning Analysis -----------------------------------------------
// Moved here from the old Assembly Check tab, this is thermodynamics, not
// structural bookkeeping, so it belongs alongside fold-change, not a graft
// checklist. Runs off whatever's already in the Calculator's own fields.

let calcDestabPreset = "conservative";
let calcLastKopenData = null;

function calcParseKopenPositions(text) {
  return text.split(",").map((s) => parseInt(s.trim(), 10)).filter((n) => !isNaN(n));
}

// I/V residues in the pasted latch sequence, skipping preserve_positions.
function calcFindKopenCandidates(latchSeq, preservePositions) {
  const candidates = [];
  for (let i = 0; i < latchSeq.length; i++) {
    const residue = latchSeq[i];
    if (residue !== "I" && residue !== "V") continue;
    const pos = i + 1;
    if (preservePositions.includes(pos)) continue;
    candidates.push({ position: pos, residue });
  }
  return candidates;
}

function calcRenderKopenCandidates(candidates, latchSeq) {
  const el = calcEl("calc-kopen-candidates");
  if (!latchSeq) {
    el.textContent = "";
  } else if (candidates.length === 0) {
    el.textContent = "No candidate I/V positions found in the pasted latch sequence, latch tuning may require custom design.";
  } else {
    const labels = candidates.map((c) => `${c.residue}${c.position}`).join(", ");
    el.textContent = `Candidate positions (numbered within your pasted snippet, position 1 = its first residue): ${labels}`;
  }
}

const CALC_KOPEN_REGIME_LABELS = { key_limited: "Key-limited", K_open_limited: "K_open-limited", mixed: "Mixed" };
const CALC_KOPEN_HELPS_BADGE = { yes: ["Yes", "badge-Low"], no: ["No", "badge-muted"],
                                 slightly_worse: ["Slightly worse", "badge-Moderate"] };

function calcRenderKopenTable(data) {
  const rows = [data.baseline, ...data.scenarios.filter((s) => s.preset === calcDestabPreset)];
  const table = calcEl("calc-kopen-table");
  table.innerHTML = "";

  const thead = document.createElement("thead");
  thead.innerHTML = `<tr>
    <th>Mutations<span class="help-tip" tabindex="0">?<span class="help-tip-body">How many I/V→S (or similar destabilizing) latch mutations this row assumes, at the destabilization preset selected above.</span></span></th>
    <th>K_open<span class="help-tip" tabindex="0">?<span class="help-tip-body">The estimated cage-looseness after making these mutations. Compare to your current K_open above.</span></span></th>
    <th>Fold-change<span class="help-tip" tabindex="0">?<span class="help-tip-body">The best-case fold-change your sensor would get with this K_open, at your current K_CK/pull/lucKey.</span></span></th>
    <th>Regime<span class="help-tip" tabindex="0">?<span class="help-tip-body">Whether the key (lucKey/K_CK) or the cage (K_open) would be the limiting factor at this K_open.</span></span></th>
    <th>Helps?<span class="help-tip" tabindex="0">?<span class="help-tip-body">Whether this K_open change actually raises fold-change by more than 5% versus your current baseline. "No" is common, most latch tuning doesn't help unless lucKey/K_CK is already close to K_open.</span></span></th>
  </tr>`;
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
    regimeTd.textContent = CALC_KOPEN_REGIME_LABELS[row.regime] || row.regime;
    tr.appendChild(regimeTd);

    const helpsTd = document.createElement("td");
    if (row.helps === "unchanged") {
      helpsTd.textContent = "N/A";
    } else {
      const [label, cls] = CALC_KOPEN_HELPS_BADGE[row.helps] || [row.helps, "badge-muted"];
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

function calcRenderKopenNotes(data, currentLuckey) {
  const mutationRows = data.scenarios.filter((s) => s.preset === calcDestabPreset);
  const anyHelps = mutationRows.some((s) => s.helps === "yes");
  const negativeEl = calcEl("calc-kopen-negative-note");
  if (!anyHelps) {
    negativeEl.style.display = "block";
    negativeEl.textContent = `At your current lucKey/K_CK ratio (${roundSig(data.dominance_ratio, 3)}), K_open changes have negligible effect. Latch tuning will only help if lucKey/K_CK is reduced to be comparable to K_open. Try lowering lucKey to ${roundSig(data.crossover_luckey_nm, 3)} nM.`;
  } else {
    negativeEl.style.display = "none";
  }

  calcEl("calc-kopen-crossover-note").textContent =
    `Latch tuning becomes meaningful when lucKey/K_CK ≈ K_open. At your current settings: lucKey/K_CK = ${roundSig(data.dominance_ratio, 3)}, K_open = ${roundSig(data.baseline.k_open, 3)}. Tuning helps when [lucKey] ≈ ${roundSig(data.crossover_luckey_nm, 3)} nM (currently ${roundSig(currentLuckey, 3)} nM).`;
}

async function calcKopenAnalyze() {
  const values = calcReadFields();
  const errEl = calcEl("calc-kopen-error");
  errEl.textContent = "";

  if (!(values.k_ck > 0) || !(values.k_open > 0) || !(values.luckey > 0) || !(values.pull >= 0)) {
    errEl.textContent = "Fill in valid K_CK, K_open, pull, and lucKey values above first.";
    return;
  }

  const latchSeq = calcEl("calc-kopen-latch").value.trim().toUpperCase();
  const preservePositions = calcParseKopenPositions(calcEl("calc-kopen-preserve").value);
  const candidates = calcFindKopenCandidates(latchSeq, preservePositions);
  calcRenderKopenCandidates(candidates, latchSeq);

  const button = calcEl("calc-kopen-analyze");
  const originalLabel = button.innerHTML;
  button.disabled = true;
  button.innerHTML = '<span class="spinner"></span>Analyzing…';

  try {
    const data = await apiPost("/kopen-scenarios", {
      k_open_current: values.k_open, k_ck: values.k_ck, luckey: values.luckey, pull: values.pull,
    });
    calcLastKopenData = data;
    calcRenderKopenTable(data);
    calcRenderKopenNotes(data, values.luckey);
    calcEl("calc-kopen-results").style.display = "block";
  } catch (err) {
    showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function initCalculator() {
  document.querySelectorAll('[data-tab="calculator"] [data-calc-field]').forEach((input) => {
    input.addEventListener("input", calcUpdateValidity);
  });

  const presetSelect = calcEl("calc-preset");
  if (presetSelect) {
    presetSelect.addEventListener("change", () => {
      calcApplyPreset(presetSelect.value);
      presetSelect.value = "";
    });
  }
  // typing in k_ck while pill is showing detaches it (editing IS detaching)
  calcEl("calc-kck").addEventListener("input", calcDetachPill);

  // dG mode toggle
  document.querySelectorAll("#calc-ktarget-mode-seg button").forEach((btn) => {
    btn.addEventListener("click", () => calcKtargetSetMode(btn.dataset.kmode));
  });
  // dG field inputs, update derived display and validity
  ["calc-dg-score", "calc-dg-ref-dg", "calc-dg-ref-kd"].forEach((id) => {
    calcEl(id).addEventListener("input", () => {
      calcDgUpdateDerived();
      calcUpdateValidity();
    });
  });

  // editing the threshold after a calc just re-checks it against the LOD already on hand
  calcEl("calc-clinical-threshold").addEventListener("input", () => {
    if (calcState.result) calcRenderLod(calcState.result, calcState.values);
  });

  calcEl("calc-submit").addEventListener("click", calcSubmit);
  calcEl("calc-reset").addEventListener("click", calcReset);
  document.addEventListener("tabchange", (e) => {
    if (e.detail.tab === "calculator") calcApplyPiped();
  });

  document.querySelectorAll("#calc-destab-segmented button").forEach((btn) => {
    btn.addEventListener("click", () => {
      calcDestabPreset = btn.dataset.destab;
      document.querySelectorAll("#calc-destab-segmented button").forEach((b) => b.classList.toggle("active", b === btn));
      if (calcLastKopenData) {
        calcRenderKopenTable(calcLastKopenData);
        calcRenderKopenNotes(calcLastKopenData, calcReadFields().luckey);
      }
    });
  });
  calcEl("calc-kopen-analyze").addEventListener("click", calcKopenAnalyze);

  initCalcKopenEstimator();
  calcKopenEstUpdate();

  calcUpdateValidity();
}

document.addEventListener("DOMContentLoaded", initCalculator);
