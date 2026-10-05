
const doseState = { result: null };
let doseTargetMode = "fraction"; // "fraction" or "foldchange"

function doseEl(id) { return document.getElementById(id); }

// uL in, human-scaled string out, mirrors formatConcNm's unit ladder.
function formatVolumeUl(ul) {
  if (ul === null || ul === undefined || Number.isNaN(ul)) return "N/A";
  const units = [["nL", 1e-3], ["uL", 1], ["mL", 1e3], ["L", 1e6]];
  for (let i = 0; i < units.length; i++) {
    const [label, scale] = units[i];
    const next = units[i + 1];
    if (!next || ul < next[1]) {
      return `${roundSig(ul / scale, 3)} ${label}`;
    }
  }
  return `${ul} uL`;
}

function doseSetTargetMode(mode) {
  doseTargetMode = mode;
  document.querySelectorAll("#dose-target-mode-seg button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.dmode === mode);
  });
  doseEl("dose-fraction-panel").style.display = mode === "fraction" ? "block" : "none";
  doseEl("dose-foldchange-panel").style.display = mode === "foldchange" ? "block" : "none";
}

function doseReadFields() {
  return {
    k_ck: parseFloat(doseEl("dose-kck").value),
    total_cage_conc: parseFloat(doseEl("dose-cage").value),
    key_stock_conc: parseFloat(doseEl("dose-stock").value),
    assay_volume: parseFloat(doseEl("dose-volume").value),
    target_fraction_bound: doseTargetMode === "fraction" ? parseFloat(doseEl("dose-target-fraction").value) : null,
    target_fold_change: doseTargetMode === "foldchange" ? parseFloat(doseEl("dose-target-fc").value) : null,
  };
}

// mirrors backend validation so Calculate disables before a bad request
function doseValidate(values) {
  const errors = {};
  ["k_ck", "total_cage_conc", "key_stock_conc", "assay_volume"].forEach((field) => {
    if (!(values[field] > 0)) errors[field] = "must be > 0";
  });
  if (doseTargetMode === "fraction") {
    const v = values.target_fraction_bound;
    if (isNaN(v) || v < 0 || v >= 1) errors.target_fraction_bound = "must be in [0, 1)";
  } else {
    const v = values.target_fold_change;
    if (isNaN(v) || v < 0) errors.target_fold_change = "must be >= 0";
  }
  return errors;
}

function doseShowFieldErrors(errors) {
  ["k_ck", "total_cage_conc", "key_stock_conc", "assay_volume", "target_fraction_bound", "target_fold_change"].forEach((field) => {
    const input = document.querySelector(`[data-dose-field="${field}"]`);
    const errEl = document.querySelector(`[data-dose-error="${field}"]`);
    const msg = errors[field];
    if (input) input.classList.toggle("invalid", Boolean(msg));
    if (errEl) errEl.textContent = msg || "";
  });
}

function doseUpdateValidity() {
  const values = doseReadFields();
  const errors = doseValidate(values);
  doseShowFieldErrors(errors);
  doseEl("dose-submit").disabled = Object.keys(errors).length > 0;
  return { values, errors };
}

// Linear-y (0-1), log-x titration plot; distinct from plot.js's drawLogXPlot,
// which hardcodes yMin=1 for fold-change curves and doesn't fit a 0-1 fraction.
function drawTitrationPlot(svg, points, markerX, markerY) {
  svg.innerHTML = "";
  const width = svg.clientWidth || 480;
  const height = 260;
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  const padding = { left: 48, right: 16, top: 16, bottom: 36 };
  const plotW = width - padding.left - padding.right;
  const plotH = height - padding.top - padding.bottom;

  const xs = points.map((p) => p.x).filter((x) => x > 0);
  if (markerX > 0) xs.push(markerX);
  const xMin = Math.min(...xs), xMax = Math.max(...xs);

  const xToPx = (x) => padding.left + ((Math.log10(x) - Math.log10(xMin)) / (Math.log10(xMax) - Math.log10(xMin))) * plotW;
  const yToPx = (y) => padding.top + plotH - y * plotH;

  const ns = "http://www.w3.org/2000/svg";
  const make = (tag, attrs) => {
    const el = document.createElementNS(ns, tag);
    Object.entries(attrs).forEach(([k, v]) => el.setAttribute(k, v));
    return el;
  };

  [0, 0.5, 1].forEach((tick) => {
    const py = yToPx(tick);
    svg.appendChild(make("line", { x1: padding.left, x2: width - padding.right, y1: py, y2: py, stroke: "var(--plot-grid)" }));
    const label = make("text", { x: padding.left - 8, y: py + 4, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "end" });
    label.textContent = tick.toFixed(1);
    svg.appendChild(label);
  });

  niceLogTicks(xMin, xMax).forEach((tick) => {
    const px = xToPx(tick);
    svg.appendChild(make("line", { x1: px, x2: px, y1: padding.top, y2: padding.top + plotH, stroke: "var(--plot-grid)" }));
    const label = make("text", { x: px, y: height - 12, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" });
    label.textContent = formatConcNm(tick);
    svg.appendChild(label);
  });

  const d = points
    .filter((p) => p.x > 0)
    .map((p, i) => `${i === 0 ? "M" : "L"} ${xToPx(p.x)} ${yToPx(p.y)}`)
    .join(" ");
  svg.appendChild(make("path", { d, fill: "none", stroke: "var(--plot-line)", "stroke-width": 2 }));

  if (markerX > 0) {
    svg.appendChild(make("circle", { cx: xToPx(markerX), cy: yToPx(markerY), r: 5, fill: "var(--plot-marker)" }));
  }
}

function doseRenderResult(result) {
  doseEl("dose-empty-state").style.display = "none";
  doseEl("dose-undefined-result").style.display = "none";
  doseEl("dose-results").style.display = "block";

  doseEl("dose-hero-volume").textContent = formatVolumeUl(result.volume_to_add_ul);
  doseEl("dose-metric-keyconc").textContent = formatConcNm(result.key_conc_nm);
  doseEl("dose-metric-fraction").textContent = roundSig(result.fraction_bound, 4);
  doseEl("dose-metric-volume").textContent = formatVolumeUl(result.volume_to_add_ul);

  const warningsEl = doseEl("dose-result-warnings");
  warningsEl.innerHTML = "";
  (result.warnings || []).forEach((w) => {
    const div = document.createElement("div");
    div.className = "field-error";
    div.textContent = w;
    warningsEl.appendChild(div);
  });

  const points = result.titration.map((p) => ({ x: p.key_conc_nm, y: p.fraction_bound }));
  if (points.length > 0) {
    drawTitrationPlot(doseEl("dose-plot-titration"), points, result.key_conc_nm, result.fraction_bound);
  }
}

async function doseSubmit() {
  const { values, errors } = doseUpdateValidity();
  if (Object.keys(errors).length > 0) return;

  const button = doseEl("dose-submit");
  const originalLabel = button.innerHTML;
  button.disabled = true;
  button.innerHTML = '<span class="spinner"></span>Calculating…';

  try {
    const result = await apiPost("/dosing", values);
    doseState.result = result;
    doseRenderResult(result);
  } catch (err) {
    if (err.code === "UNDEFINED_RESULT") {
      doseEl("dose-empty-state").style.display = "none";
      doseEl("dose-results").style.display = "none";
      doseEl("dose-undefined-result").style.display = "block";
    } else {
      showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
    }
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function doseReset() {
  doseEl("dose-form").reset();
  doseState.result = null;
  doseSetTargetMode("fraction");
  doseUpdateValidity();
  doseEl("dose-empty-state").style.display = "block";
  doseEl("dose-results").style.display = "none";
  doseEl("dose-undefined-result").style.display = "none";
}

function initDosing() {
  document.querySelectorAll('#dose-section [data-dose-field]').forEach((input) => {
    input.addEventListener("input", doseUpdateValidity);
  });
  document.querySelectorAll("#dose-target-mode-seg button").forEach((btn) => {
    btn.addEventListener("click", () => { doseSetTargetMode(btn.dataset.dmode); doseUpdateValidity(); });
  });
  // One K_CK for the whole tab: the dosing box follows the calculator's unless edited here.
  const calcKck = document.getElementById("calc-kck");
  const doseKck = doseEl("dose-kck");
  if (calcKck && doseKck) {
    let overridden = false;
    doseKck.addEventListener("input", () => { overridden = true; });
    doseKck.value = calcKck.value;
    calcKck.addEventListener("input", () => {
      if (overridden) return;
      doseKck.value = calcKck.value;
      doseUpdateValidity();
    });
  }
  doseEl("dose-submit").addEventListener("click", doseSubmit);
  doseEl("dose-reset").addEventListener("click", doseReset);
  doseUpdateValidity();
}

document.addEventListener("DOMContentLoaded", initDosing);
