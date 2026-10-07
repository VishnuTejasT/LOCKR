// Temperature Window view: p_works and fold change over a temperature sweep.
// Inline SVG on a linear temperature axis. plot.js is log-X only, so this draws its own.

// Okabe-Ito: distinct under every common form of colour blindness.
const TW_COLORS = ["#0072b2", "#d55e00", "#009e73", "#cc79a7",
                   "#e69f00", "#56b4e9", "#8c564b", "#7f7f7f",
                   "#117733", "#882255"];

const TW_FIELDS = ["dG_latch_cage", "dH_latch_cage", "dG_key_cage", "dH_key_cage",
                   "key_conc", "target_conc", "t_max_stability_c"];
const TW_BINDER_FIELDS = ["name", "sequence", "dG_target_binder", "kd_nM",
                          "dG_target_binder_sd", "dH_target_binder", "tm_c",
                          "graft_dG_latch_cage_shift"];

// The API demands a binder affinity even when it cannot matter. With no target amount the
// target is saturating and this value cancels out of the ratio, so it is inert, not a guess.
const TW_INERT_BINDER_DG = -10.0;

const twState = { data: null, request: null, placeholder: false };

function twEl(id) {
  return document.getElementById(id);
}

function twNum(input) {
  // "" means "not supplied", which is not the same as 0 for any field here.
  if (!input || input.value.trim() === "") return null;
  const value = Number(input.value);
  return Number.isFinite(value) ? value : null;
}

function twFmt(v, digits = 3) {
  if (v === null || v === undefined || !Number.isFinite(v)) return "N/A";
  return roundSig(v, digits);
}

// --- the construct form ---------------------------------------------------------------

function twConstructCount() {
  return twEl("tw-constructs").querySelectorAll(".tw-construct").length;
}

function twAddConstruct(values) {
  const index = twConstructCount();
  if (index >= 10) {
    showToast("Ten constructs is the limit for one comparison.");
    return null;
  }
  const row = document.createElement("div");
  row.className = "card tw-construct";
  row.dataset.index = String(index);
  row.innerHTML = `
    <div class="field-row" style="display:flex; align-items:center; gap:8px; justify-content:space-between;">
      <strong class="tw-construct-title">Construct ${index + 1}</strong>
      <button type="button" class="btn-ghost tw-remove" ${index === 0 ? "hidden" : ""}>Remove</button>
    </div>
    <div class="field-row">
      <label>Name</label>
      <input data-f="name" type="text" placeholder="cage 1">
    </div>
    <div class="field-row">
      <label>Cage tightness <span class="unit-suffix">kcal/mol</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">How hard the cage holds itself shut. More negative = tighter cage = less glow when the target is absent. Around -4 is typical. Required.</span></span></label>
      <input data-f="dG_latch_cage" type="number" step="any" placeholder="-4">
    </div>
    <div class="field-row">
      <label>Key grip <span class="unit-suffix">kcal/mol</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">How tightly the glowing key holds the open cage. More negative = tighter grip = brighter. Around -10.9 is the standard lucCage value. Required.</span></span></label>
      <input data-f="dG_key_cage" type="number" step="any" placeholder="-10.9">
    </div>

    <label class="tw-toggle">
      <input type="checkbox" class="tw-binder-toggle">
      Include a specific binder and target amount
      <span class="help-tip" tabindex="0">?<span class="help-tip-body">Leave this off to compare cages alone. The binder's strength only changes the answer when you also say how much target is present, so the two are asked for together.</span></span>
    </label>

    <fieldset class="tw-binder" hidden>
      <legend>Binder and target</legend>
      <div class="field-row">
        <label>Target amount <span class="unit-suffix">nM</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">How much of the thing you are detecting is in the tube. This is what makes the binder strength matter at all.</span></span></label>
        <input data-f="target_conc" type="number" step="any" min="0" placeholder="e.g. 10">
      </div>
      <div class="field-row">
        <label>Binder strength <span class="unit-suffix">kcal/mol</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">How tightly the binder grabs the target. From docking or a measurement. Give this OR the Kd below, not both.</span></span></label>
        <input data-b="dG_target_binder" type="number" step="any" placeholder="-12">
      </div>
      <div class="field-row">
        <label>...or as a K<sub>d</sub> <span class="unit-suffix">nM</span></label>
        <input data-b="kd_nM" type="number" step="any" placeholder="instead of the number above">
      </div>
    </fieldset>

    <details class="advanced">
      <summary>More options for this construct</summary>
      <div class="card">
        <div class="field-row">
          <label>Cage heat response, &Delta;H <span class="unit-suffix">kcal/mol</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">Heat released when the cage closes. This is what makes the cage temperature-sensitive. Blank means it is drawn from the placeholder prior.</span></span></label>
          <input data-f="dH_latch_cage" type="number" step="any" placeholder="drawn from prior">
        </div>
        <div class="field-row">
          <label>Key-grip heat response, &Delta;H <span class="unit-suffix">kcal/mol</span></label>
          <input data-f="dH_key_cage" type="number" step="any" placeholder="drawn from prior">
        </div>
        <div class="field-row">
          <label>Key amount <span class="unit-suffix">nM</span></label>
          <input data-f="key_conc" type="number" step="any" min="0" placeholder="500">
        </div>
        <div class="field-row">
          <label>Falls apart above <span class="unit-suffix">°C</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">Highest temperature the construct survives. Above it the design is scored as not working at all.</span></span></label>
          <input data-f="t_max_stability_c" type="number" step="any" min="0" max="100" placeholder="optional">
        </div>
        <div class="field-row">
          <label>Binder name</label>
          <input data-b="name" type="text" placeholder="binder 1">
        </div>
        <div class="field-row">
          <label>Binder sequence<span class="help-tip" tabindex="0">?<span class="help-tip-body">Written down for your records and shown in the legend. It is never used to predict anything.</span></span></label>
          <input data-b="sequence" type="text" placeholder="optional, for your records">
        </div>
        <div class="field-row">
          <label>How unsure is the binder strength? <span class="unit-suffix">kcal/mol</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">A docking score deserves a wide spread (1.5 or more); an ITC measurement deserves a narrow one.</span></span></label>
          <input data-b="dG_target_binder_sd" type="number" step="any" min="0" placeholder="1.5">
        </div>
        <div class="field-row">
          <label>Binder heat response, &Delta;H <span class="unit-suffix">kcal/mol</span></label>
          <input data-b="dH_target_binder" type="number" step="any" placeholder="drawn from prior">
        </div>
        <div class="field-row">
          <label>Binder falls apart above <span class="unit-suffix">°C</span></label>
          <input data-b="tm_c" type="number" step="any" min="0" max="100" placeholder="optional">
        </div>
        <div class="field-row">
          <label>Grafting loosens the cage by <span class="unit-suffix">kcal/mol</span><span class="help-tip" tabindex="0">?<span class="help-tip-body">How much attaching this binder changes the cage's closing energy. Positive loosens it and raises background glow. 0 assumes grafting is free.</span></span></label>
          <input data-b="graft_dG_latch_cage_shift" type="number" step="any" placeholder="0">
        </div>
      </div>
    </details>`;

  const toggle = row.querySelector(".tw-binder-toggle");
  toggle.addEventListener("change", () => {
    row.querySelector(".tw-binder").hidden = !toggle.checked;
  });
  row.querySelector(".tw-remove").addEventListener("click", () => {
    row.remove();
    twRenumberConstructs();
  });
  twEl("tw-constructs").appendChild(row);
  if (values) twFillConstruct(row, values);
  return row;
}

function twRenumberConstructs() {
  twEl("tw-constructs").querySelectorAll(".tw-construct").forEach((row, i) => {
    row.dataset.index = String(i);
    row.querySelector(".tw-construct-title").textContent = `Construct ${i + 1}`;
    row.querySelector(".tw-remove").hidden = i === 0 && twConstructCount() === 1;
  });
}

function twFillConstruct(row, values) {
  const set = (selector, value) => {
    const input = row.querySelector(selector);
    if (input) input.value = value === null || value === undefined ? "" : value;
  };
  set('[data-f="name"]', values.name);
  TW_FIELDS.forEach((f) => set(`[data-f="${f}"]`, values[f]));
  const binder = values.binder || {};
  TW_BINDER_FIELDS.forEach((f) => set(`[data-b="${f}"]`, binder[f]));
}

function twReadConstruct(row) {
  const system = { name: row.querySelector('[data-f="name"]').value.trim() || null };
  TW_FIELDS.forEach((f) => { system[f] = twNum(row.querySelector(`[data-f="${f}"]`)); });

  const binder = {};
  const nameInput = row.querySelector('[data-b="name"]');
  const seqInput = row.querySelector('[data-b="sequence"]');
  binder.name = nameInput.value.trim() || null;
  binder.sequence = seqInput.value.trim() || null;
  TW_BINDER_FIELDS.filter((f) => f !== "name" && f !== "sequence")
    .forEach((f) => { binder[f] = twNum(row.querySelector(`[data-b="${f}"]`)); });
  if (!row.querySelector(".tw-binder-toggle").checked) {
    binder.dG_target_binder = TW_INERT_BINDER_DG;
    binder.kd_nM = null;
    system.target_conc = null;
  }
  system.binder = binder;
  return system;
}

// --- request / response ---------------------------------------------------------------

function twReadLuciferase() {
  const raw = twEl("tw-luciferase").value.trim();
  if (!raw) return null;
  // "4:1.0, 25:0.9" or one "temp, activity" pair per line
  const points = raw.split(/[\n;]+/).map((line) => line.trim()).filter(Boolean)
    .map((line) => line.split(/[,:\s]+/).map(Number))
    .filter((pair) => pair.length === 2 && pair.every(Number.isFinite));
  return points.length ? points : null;
}

function twPriorFrom(prefix) {
  const mean = twNum(twEl(`tw-prior-${prefix}-mean`));
  const sd = twNum(twEl(`tw-prior-${prefix}-sd`));
  return mean === null && sd === null ? null : { mean, sd };
}

function twBuildRequest(seed) {
  const priors = {
    latch_cage: twPriorFrom("latch"),
    binder_target: twPriorFrom("binder"),
    key_cage: twPriorFrom("key"),
  };
  const anyPrior = Object.values(priors).some((p) => p !== null);
  return {
    systems: Array.from(twEl("tw-constructs").querySelectorAll(".tw-construct"))
      .map(twReadConstruct),
    t_min_c: twNum(twEl("tw-t-min")) ?? 4.0,
    t_max_c: twNum(twEl("tw-t-max")) ?? 45.0,
    t_step_c: twNum(twEl("tw-t-step")) ?? 1.0,
    t_ref_c: twNum(twEl("tw-t-ref")) ?? 25.0,
    fc_threshold: twNum(twEl("tw-fc-threshold")) ?? 2.0,
    confidence: twNum(twEl("tw-confidence")) ?? 0.8,
    n_samples: twNum(twEl("tw-n-samples")) ?? 2000,
    dH_priors: anyPrior ? priors : null,
    seed: seed === undefined ? twNum(twEl("tw-seed")) : seed,
    luciferase_activity: twEl("tw-use-luciferase").checked ? twReadLuciferase() : null,
  };
}

async function twRun(seed) {
  const error = twEl("tw-error");
  error.style.display = "none";
  twEl("tw-submit").disabled = true;
  try {
    const request = twBuildRequest(seed);
    twState.data = await apiPost("/api/temperature/window", request);
    twState.request = request;
    twRender();
  } catch (err) {
    error.textContent = err.message || "Something went wrong.";
    error.style.display = "block";
  } finally {
    twEl("tw-submit").disabled = false;
  }
}

// --- plotting -------------------------------------------------------------------------

const TW_NS = "http://www.w3.org/2000/svg";

function twMake(tag, attrs) {
  const el = document.createElementNS(TW_NS, tag);
  Object.entries(attrs).forEach(([k, v]) => el.setAttribute(k, v));
  return el;
}

function twAxes(svg, { width, height, pad, xMin, xMax, yMin, yMax, yTicks, yLabel, xLabel }) {
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;
  const xToPx = (x) => pad.left + ((x - xMin) / (xMax - xMin || 1)) * plotW;
  const yToPx = (y) => pad.top + plotH - ((y - yMin) / (yMax - yMin || 1)) * plotH;

  // x gridlines every 5 degrees
  const start = Math.ceil(xMin / 5) * 5;
  for (let t = start; t <= xMax; t += 5) {
    const px = xToPx(t);
    svg.appendChild(twMake("line", { x1: px, x2: px, y1: pad.top, y2: pad.top + plotH,
                                     stroke: "var(--plot-grid)" }));
    const label = twMake("text", { x: px, y: height - 14, "font-size": 11,
                                   fill: "var(--text-muted)", "text-anchor": "middle" });
    label.textContent = `${t}`;
    svg.appendChild(label);
  }
  yTicks.forEach((t) => {
    const py = yToPx(t);
    svg.appendChild(twMake("line", { x1: pad.left, x2: width - pad.right, y1: py, y2: py,
                                     stroke: "var(--plot-grid)" }));
    const label = twMake("text", { x: pad.left - 8, y: py + 4, "font-size": 11,
                                   fill: "var(--text-muted)", "text-anchor": "end" });
    label.textContent = twFmt(t, 3);
    svg.appendChild(label);
  });

  const xl = twMake("text", { x: pad.left + plotW / 2, y: height - 1, "font-size": 11,
                              fill: "var(--text-muted)", "text-anchor": "middle" });
  xl.textContent = xLabel;
  svg.appendChild(xl);
  const yl = twMake("text", { x: 10, y: pad.top + 10, "font-size": 11,
                              fill: "var(--text-muted)" });
  yl.textContent = yLabel;
  svg.appendChild(yl);
  return { xToPx, yToPx, plotH, plotW };
}

function twPath(points, color, { dashed = false, width = 2 } = {}) {
  const d = points.map((p, i) => `${i === 0 ? "M" : "L"} ${p[0]} ${p[1]}`).join(" ");
  return twMake("path", { d, fill: "none", stroke: color, "stroke-width": width,
                          "stroke-dasharray": dashed ? "5,4" : "none" });
}

function twDrawPWorks() {
  const svg = twEl("tw-plot-pworks");
  svg.innerHTML = "";
  const data = twState.data;
  const temps = data.temperatures_c;
  const width = svg.clientWidth || 560;
  const height = 240;
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  const pad = { left: 50, right: 14, top: 16, bottom: 32 };
  const ax = twAxes(svg, { width, height, pad, xMin: temps[0], xMax: temps[temps.length - 1],
                           yMin: 0, yMax: 1, yTicks: [0, 0.25, 0.5, 0.75, 1],
                           yLabel: "chance it works", xLabel: "temperature (°C)" });

  // operating windows first, so the curves sit on top of the shading
  data.systems.forEach((s, i) => {
    if (!s.window) return;
    const x1 = ax.xToPx(s.window.t_low_c);
    const x2 = ax.xToPx(s.window.t_high_c);
    svg.appendChild(twMake("rect", {
      x: Math.min(x1, x2), y: pad.top, width: Math.max(2, Math.abs(x2 - x1)), height: ax.plotH,
      fill: TW_COLORS[i % TW_COLORS.length], opacity: 0.12,
    }));
  });

  // the confidence bar being cleared
  const cy = ax.yToPx(data.echo.confidence);
  svg.appendChild(twMake("line", { x1: pad.left, x2: width - pad.right, y1: cy, y2: cy,
                                   stroke: "var(--text-muted)", "stroke-width": 1.5,
                                   "stroke-dasharray": "4,4" }));
  const clabel = twMake("text", { x: width - pad.right - 4, y: cy - 5, "font-size": 10,
                                  fill: "var(--text-muted)", "text-anchor": "end" });
  clabel.textContent = `confidence ${data.echo.confidence}`;
  svg.appendChild(clabel);

  data.systems.forEach((s, i) => {
    const color = TW_COLORS[i % TW_COLORS.length];
    svg.appendChild(twPath(s.p_works.map((p, j) => [ax.xToPx(temps[j]), ax.yToPx(p)]), color));
    if (s.stability_cap_c !== null && s.stability_cap_c <= temps[temps.length - 1]) {
      const px = ax.xToPx(s.stability_cap_c);
      svg.appendChild(twMake("line", { x1: px, x2: px, y1: pad.top, y2: pad.top + ax.plotH,
                                       stroke: color, "stroke-width": 1,
                                       "stroke-dasharray": "2,3" }));
    }
  });
}

function twDrawFoldChange() {
  const svg = twEl("tw-plot-fc");
  svg.innerHTML = "";
  const data = twState.data;
  const temps = data.temperatures_c;
  const width = svg.clientWidth || 560;
  const height = 260;
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  const finite = (arr) => arr.filter((v) => v !== null && Number.isFinite(v));
  const highs = data.systems.flatMap((s) => finite(s.p90_fc));
  const lows = data.systems.flatMap((s) => finite(s.p10_fc));
  const yMax = highs.length ? Math.max(...highs) * 1.08 : 2;
  // fold change starts at 1 (no change), matching plot.js
  const yMin = Math.min(1, lows.length ? Math.min(...lows) : 1);

  const pad = { left: 50, right: 14, top: 16, bottom: 32 };
  const ax = twAxes(svg, { width, height, pad, xMin: temps[0], xMax: temps[temps.length - 1],
                           yMin, yMax,
                           yTicks: [yMin, (yMin + yMax) / 2, yMax],
                           yLabel: "brightness jump (×)",
                           xLabel: "temperature (°C)" });

  // the threshold the window is measured against
  if (data.echo.fc_threshold >= yMin && data.echo.fc_threshold <= yMax) {
    const ty = ax.yToPx(data.echo.fc_threshold);
    svg.appendChild(twMake("line", { x1: pad.left, x2: width - pad.right, y1: ty, y2: ty,
                                     stroke: "var(--text-muted)", "stroke-width": 1.5,
                                     "stroke-dasharray": "4,4" }));
    const label = twMake("text", { x: width - pad.right - 4, y: ty - 5, "font-size": 10,
                                   fill: "var(--text-muted)", "text-anchor": "end" });
    label.textContent = `threshold ${data.echo.fc_threshold}×`;
    svg.appendChild(label);
  }

  data.systems.forEach((s, i) => {
    const color = TW_COLORS[i % TW_COLORS.length];
    const usable = temps.map((t, j) => [t, s.p10_fc[j], s.p90_fc[j], s.median_fc[j]])
      .filter((p) => p[1] !== null && p[2] !== null && p[3] !== null);
    if (!usable.length) return;

    // p10-p90 band
    const top = usable.map((p) => `${ax.xToPx(p[0])} ${ax.yToPx(p[2])}`);
    const bottom = usable.slice().reverse().map((p) => `${ax.xToPx(p[0])} ${ax.yToPx(p[1])}`);
    svg.appendChild(twMake("polygon", { points: top.concat(bottom).join(" "),
                                        fill: color, opacity: 0.15 }));
    svg.appendChild(twPath(usable.map((p) => [ax.xToPx(p[0]), ax.yToPx(p[3])]), color));

    if (s.t_opt_c !== null && s.peak_median_fc !== null) {
      const cx = ax.xToPx(s.t_opt_c);
      const cy = ax.yToPx(s.peak_median_fc);
      svg.appendChild(twMake("circle", { cx, cy, r: 4.5, fill: color,
                                         stroke: "var(--surface)", "stroke-width": 1.5 }));
      const label = twMake("text", { x: cx + 7, y: cy - 7, "font-size": 11, fill: color,
                                     "font-weight": 600 });
      label.textContent = `${twFmt(s.t_opt_c, 3)}°C`;
      svg.appendChild(label);
    }
  });
}

function twRenderLegend() {
  twEl("tw-legend").innerHTML = twState.data.systems.map((s, i) => `
    <span class="tw-legend-item">
      <span class="tw-swatch" style="background:${TW_COLORS[i % TW_COLORS.length]}"></span>
      ${twEscape(s.name)} / ${twEscape(s.binder_name)}
    </span>`).join("");
}

function twEscape(text) {
  const div = document.createElement("div");
  div.textContent = text === null || text === undefined ? "" : String(text);
  return div.innerHTML;
}

function twWindowText(s) {
  if (!s.window) return "none";
  const { t_low_c, t_high_c } = s.window;
  const span = t_low_c === t_high_c
    ? `${twFmt(t_low_c, 3)}°C only`
    : `${twFmt(t_low_c, 3)} to ${twFmt(t_high_c, 3)}°C`;
  const extra = s.other_segments.length ? ` (+${s.other_segments.length} more)` : "";
  return span + extra;
}

function twRenderTable() {
  const rows = twState.data.systems.map((s, i) => `
    <tr>
      <td><span class="tw-swatch" style="background:${TW_COLORS[i % TW_COLORS.length]}"></span>${twEscape(s.name)}</td>
      <td>${twEscape(s.binder_name)}</td>
      <td>${s.t_opt_c === null ? "N/A" : `${twFmt(s.t_opt_c, 3)}°C`}</td>
      <td>${twWindowText(s)}</td>
      <td>${s.peak_median_fc === null ? "N/A" : `${twFmt(s.peak_median_fc, 3)}×`}</td>
      <td>${twFmt(s.dG_latch_cage_effective, 4)}${s.graft_dG_latch_cage_shift
            ? ` <span class="help-text">(raw ${twFmt(s.dG_latch_cage_raw, 4)})</span>` : ""}</td>
      <td>${s.stability_cap_source
            ? `${twEscape(s.stability_cap_source)} @ ${twFmt(s.stability_cap_c, 3)}°C` : "none"}</td>
      <td>${s.window_touches_grid_edge ? '<span class="badge badge-marginal">at edge</span>' : "no"}</td>
    </tr>`).join("");

  twEl("tw-table").innerHTML = `
    <table class="batch-table">
      <thead><tr>
        <th>Construct</th><th>Binder</th><th>Best temp</th><th>Working window</th>
        <th>Peak jump</th><th>Latch-Cage &Delta;G</th><th>Stability cap</th><th>Edge</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
}

function twRenderNoWindow() {
  const failed = twState.data.systems.filter((s) => !s.window);
  const box = twEl("tw-no-window");
  if (!failed.length) {
    box.style.display = "none";
    return;
  }
  box.style.display = "block";
  box.innerHTML = failed.map((s) => `
    <p><strong>No reliable window at this threshold</strong> for
    ${twEscape(s.name)} / ${twEscape(s.binder_name)}. ${twEscape(s.reason)}</p>`).join("");
}

function twBinderModeOn() {
  return Array.from(document.querySelectorAll(".tw-binder-toggle")).some((t) => t.checked);
}

function twRenderBarSummary() {
  const bar = twNum(twEl("tw-fc-threshold")) ?? 2;
  const confidence = twNum(twEl("tw-confidence")) ?? 0.8;
  const from = twNum(twEl("tw-t-min")) ?? 4;
  const to = twNum(twEl("tw-t-max")) ?? 45;
  twEl("tw-bar-summary").innerHTML =
    `A design counts as working when it gets at least <b>${twFmt(bar)}x</b> brighter, in at least ` +
    `<b>${Math.round(confidence * 100)}%</b> of simulated runs. Sweeping <b>${twFmt(from)} to ` +
    `${twFmt(to)} °C</b>. <button type="button" data-tw-open-advanced>Change</button>`;
}

// With no target amount the engine says the binder did not matter. True, but confusing when
// the binder fields were never on screen, so say what it means for the run they actually did.
function twNote(note) {
  if (!twBinderModeOn() && note.includes("does NOT affect this result")) {
    return "Binder strength is not part of this run: you are comparing cages on their own. "
      + "Tick \u201cInclude a specific binder and target amount\u201d on a construct to bring it in.";
  }
  return note;
}

function twRender() {
  const data = twState.data;
  twEl("tw-results").style.display = "block";
  twDrawPWorks();
  twDrawFoldChange();
  twRenderLegend();
  twRenderTable();
  twRenderNoWindow();

  twEl("tw-seed").value = data.seed;
  twEl("tw-seed-readout").textContent = `Seed ${data.seed}`;
  // notes are a standing condition of the answer, so they are never collapsed away
  twEl("tw-notes").innerHTML = `<h3 class="font-display" style="font-size:15px;">What this answer assumes</h3>
    <ul>${[...new Set(data.notes.map(twNote))].map((n) => `<li>${twEscape(n)}</li>`).join("")}</ul>`;
}

// --- preset ---------------------------------------------------------------------------

async function twLoadPreset() {
  let preset;
  try {
    const response = await fetch("/api/temperature/window/preset");
    preset = await response.json();
  } catch (err) {
    showToast("Could not load the example.");
    return;
  }
  const request = preset.request;
  twEl("tw-constructs").innerHTML = "";
  request.systems.forEach((s) => twAddConstruct(s));
  twRenumberConstructs();

  twEl("tw-t-min").value = request.t_min_c;
  twEl("tw-t-max").value = request.t_max_c;
  twEl("tw-t-step").value = request.t_step_c;
  twEl("tw-t-ref").value = request.t_ref_c;
  twEl("tw-fc-threshold").value = request.fc_threshold;
  twEl("tw-confidence").value = request.confidence;
  twEl("tw-n-samples").value = request.n_samples;
  twEl("tw-seed").value = "";

  twState.placeholder = Boolean(preset.placeholder);
  const badge = twEl("tw-placeholder-badge");
  badge.style.display = twState.placeholder ? "inline-block" : "none";
  badge.textContent = preset.badge || "Example values, replace with real data";
  badge.title = `Placeholder fields: ${(preset.placeholder_fields || []).join(", ")}`;
}

// --- wiring ---------------------------------------------------------------------------

function initTemperatureWindow() {
  if (!twEl("tw-constructs")) return;

  twAddConstruct();
  twRenumberConstructs();

  twEl("tw-add-construct").addEventListener("click", () => {
    twAddConstruct();
    twRenumberConstructs();
  });
  twEl("tw-load-preset").addEventListener("click", twLoadPreset);
  twEl("tw-submit").addEventListener("click", () => twRun());
  twEl("tw-rerun-seed").addEventListener("click", () => {
    if (!twState.data) return;
    twRun(twState.data.seed);
  });
  twEl("tw-use-luciferase").addEventListener("change", (e) => {
    twEl("tw-luciferase-row").style.display = e.target.checked ? "block" : "none";
  });

  ["tw-fc-threshold", "tw-confidence", "tw-t-min", "tw-t-max"].forEach((id) => {
    twEl(id).addEventListener("input", twRenderBarSummary);
  });
  twEl("tw-bar-summary").addEventListener("click", (e) => {
    if (e.target.closest("[data-tw-open-advanced]")) twEl("tw-advanced").open = true;
  });
  twRenderBarSummary();

  window.addEventListener("resize", () => {
    if (twState.data) { twDrawPWorks(); twDrawFoldChange(); }
  });
  document.addEventListener("tabchange", (e) => {
    if (e.detail.tab === "temperature" && twState.data) {
      twDrawPWorks();
      twDrawFoldChange();
    }
  });
}

document.addEventListener("DOMContentLoaded", initTemperatureWindow);
