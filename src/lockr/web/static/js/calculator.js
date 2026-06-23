
const calcState = { result: null, sweepLuckey: null, sweepKopen: null };

// "kd" or "dg" -- which mode the K_target entry is in
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
  if (values.k_open > 0 && values.pull >= 0) {
    el.textContent = `K_open (effective, target-bound) = ${roundSig(values.k_open * (1 + values.pull), 4)}`;
  } else {
    el.textContent = "K_open (effective, target-bound) = —";
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
  key_limited: "Key-limited regime.",
  K_open_limited: "K_open-limited regime.",
  mixed: "Mixed regime.",
};

function calcRenderVerdict(result) {
  calcEl("calc-empty-state").style.display = "none";
  calcEl("calc-results").style.display = "block";

  calcEl("calc-hero-fc").textContent = formatFoldChange(result.fold_change);
  calcEl("calc-subline").textContent =
    `lucKey/K_CK dominance ratio = ${roundSig(result.dominance_ratio, 3)} · ` +
    `${roundSig(result.fraction_of_dominance_ratio * 100, 3)}% of dominance ratio`;

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

  const luckeySweep = await apiPost("/sweep", { base_params, sweep: { param: "luckey", ...autoSweepRange(values.luckey) } });
  drawLogXPlot(calcEl("calc-plot-luckey"), {
    xIsConcentration: true,
    lines: [
      { points: luckeySweep.points.map((p) => ({ x: p.x, y: p.fold_change })), color: "var(--plot-line)" },
      { points: luckeySweep.points.map((p) => ({ x: p.x, y: p.dominance_ratio })), color: "var(--plot-ceiling)", dashed: true },
    ],
    marker: { x: luckeySweep.operating_point.x, y: luckeySweep.operating_point.fold_change, label: formatConcNm(luckeySweep.operating_point.x) },
  });

  const kopenSweep = await apiPost("/sweep", { base_params, sweep: { param: "k_open", ...autoSweepRange(values.k_open) } });
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
    calcRenderVerdict(result);
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

function initCalculator() {
  document.querySelectorAll('[data-tab="calculator"] [data-calc-field]').forEach((input) => {
    input.addEventListener("input", calcUpdateValidity);
  });
  // typing in k_ck while pill is showing detaches it (editing IS detaching)
  calcEl("calc-kck").addEventListener("input", calcDetachPill);

  // dG mode toggle
  document.querySelectorAll("#calc-ktarget-mode-seg button").forEach((btn) => {
    btn.addEventListener("click", () => calcKtargetSetMode(btn.dataset.kmode));
  });
  // dG field inputs -- update derived display and validity
  ["calc-dg-score", "calc-dg-ref-dg", "calc-dg-ref-kd"].forEach((id) => {
    calcEl(id).addEventListener("input", () => {
      calcDgUpdateDerived();
      calcUpdateValidity();
    });
  });

  calcEl("calc-submit").addEventListener("click", calcSubmit);
  calcEl("calc-reset").addEventListener("click", calcReset);
  document.addEventListener("tabchange", (e) => {
    if (e.detail.tab === "calculator") calcApplyPiped();
  });
  calcUpdateValidity();
}

document.addEventListener("DOMContentLoaded", initCalculator);
