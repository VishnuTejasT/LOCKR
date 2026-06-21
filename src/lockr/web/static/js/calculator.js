
const calcState = { result: null, sweepLuckey: null, sweepKopen: null };

function calcEl(id) { return document.getElementById(id); }

function calcReadFields() {
  return {
    k_ck: parseFloat(calcEl("calc-kck").value),
    k_open: parseFloat(calcEl("calc-kopen").value),
    pull: parseFloat(calcEl("calc-pull").value),
    luckey: parseFloat(calcEl("calc-luckey").value),
    k_target: calcEl("calc-ktarget").value === "" ? null : parseFloat(calcEl("calc-ktarget").value),
    target_conc: calcEl("calc-targetconc").value === "" ? null : parseFloat(calcEl("calc-targetconc").value),
  };
}

// Mirrors the backend's own validation (positive core params, pull >= 0,
// k_target/target_conc paired) so Calculate disables before a bad request goes out.
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
    const input = document.querySelector(`[data-calc-field="${field}"]`);
    const errEl = document.querySelector(`[data-calc-error="${field}"]`);
    const msg = errors[field];
    input.classList.toggle("invalid", Boolean(msg));
    if (errEl) errEl.textContent = msg || "";
  });
}

// pull > 50 has no precedent in any documented ECLIPSE run... worth a nudge,
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
    calcState.result = result;
    calcRenderVerdict(result);
    await calcRunSweeps(values);
  } catch (err) {
    showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function calcReset() {
  calcEl("calc-form").reset();
  calcUpdateValidity();
  calcEl("calc-empty-state").style.display = "block";
  calcEl("calc-results").style.display = "none";
}

function initCalculator() {
  document.querySelectorAll('[data-tab="calculator"] [data-calc-field]').forEach((input) => {
    input.addEventListener("input", calcUpdateValidity);
  });
  calcEl("calc-submit").addEventListener("click", calcSubmit);
  calcEl("calc-reset").addEventListener("click", calcReset);
  calcUpdateValidity();
}

document.addEventListener("DOMContentLoaded", initCalculator);
