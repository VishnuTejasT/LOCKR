// Module D: canvas heatmap, contours, feasibility overlay.

const LS_AXIS_PARAMS = ["k_ck", "k_open", "luckey", "pull", "temperature_K"];
const LS_METRICS = ["fold_change", "on_signal", "off_signal", "dynamic_range",
                    "limit_of_detection", "fraction_open"];
const LS_FIXED_FIELDS = ["k_ck", "k_open", "luckey", "pull", "k_target", "target_conc"];

// Nothing user-facing says k_ck or fold_change. The API keeps its names; the UI does not.
const LS_LABELS = {
  k_ck: "Key grip strength",
  k_open: "Cage looseness",
  luckey: "Key amount",
  pull: "Target pull strength",
  temperature_K: "Temperature",
  k_target: "Target grip strength",
  target_conc: "Target amount",
};

const LS_METRIC_LABELS = {
  fold_change: "Brightness jump",
  on_signal: "Brightness with target",
  off_signal: "Background glow",
  dynamic_range: "Signal gap",
  limit_of_detection: "Smallest detectable amount",
  fraction_open: "Fraction of cages open",
};

const LS_OP_LABELS = { ">=": "at least", "<=": "at most", ">": "more than", "<": "less than" };

const LS_METRIC_UNITS = { fold_change: "x", limit_of_detection: " nM" };

function lsUnitSuffix(unit) {
  return unit === "nM" || unit === "K" ? ` (${unit})` : "";
}

function lsAxisTitle(axis) {
  return LS_LABELS[axis.param] + lsUnitSuffix(axis.unit);
}

// "0.0000455" beats "4.552e-05" for anyone who is not already a biophysicist.
function lsHumanNumber(v, sig = 3) {
  if (v === null || v === undefined || !Number.isFinite(v)) return "N/A";
  if (v === 0) return "0";
  const magnitude = Math.abs(v);
  if (magnitude >= 1e-6 && magnitude < 1e7) {
    const decimals = Math.max(0, sig - 1 - Math.floor(Math.log10(magnitude)));
    return Number(v.toFixed(Math.min(20, decimals))).toString();
  }
  return roundSig(v, sig);
}

// kcal/mol/K, matching temperature.py's R_KCAL
const LS_R_KCAL = 1.9872041e-3;

// Viridis: perceptually uniform and colorblind-safe. Not rainbow, not red-green.
const LS_VIRIDIS = [
  [68, 1, 84], [72, 40, 120], [62, 74, 137], [49, 104, 142], [38, 130, 142],
  [31, 158, 137], [53, 183, 121], [109, 205, 89], [180, 222, 44], [253, 231, 37],
];

const LS_PAD = { left: 72, right: 78, top: 20, bottom: 52 };

const lsState = {
  data: null,
  finite: null,     // {min, max} of the plotted metric, for the colour scale
  pin: null,        // {i, j}
  hover: null,      // {i, j}
  inverse: null,    // [{x, y, label}] mapped Module C solutions
  inverseNote: "",
  presetNote: "",
};

function lsEl(id) { return document.getElementById(id); }
function lsNum(id) {
  const raw = lsEl(id).value.trim();
  return raw === "" ? null : Number(raw);
}

function lsAxisSpec(prefix) {
  return {
    param: lsEl(`${prefix}-param`).value,
    min: lsNum(`${prefix}-min`),
    max: lsNum(`${prefix}-max`),
    steps: lsNum(`${prefix}-steps`),
    scale: lsEl(`${prefix}-scale`).dataset.scale,
  };
}

function lsThermoBlock() {
  const block = {};
  const latch = { dH_ref: lsNum("ls-latch-dh"), dCp: lsNum("ls-latch-dcp"), Tm: lsNum("ls-latch-tm") };
  const key = { dH_ref: lsNum("ls-key-dh") };
  if (latch.dH_ref !== null || latch.Tm !== null) block.dG_latch_cage = latch;
  if (key.dH_ref !== null) block.dG_key_cage = key;
  return block;
}

function lsConstraints() {
  return Array.from(document.querySelectorAll("#ls-constraints .ls-constraint"))
    .map((row) => ({
      metric: row.querySelector("[data-c-metric]").value,
      op: row.querySelector("[data-c-op]").value,
      value: Number(row.querySelector("[data-c-value]").value),
    }))
    .filter((c) => Number.isFinite(c.value));
}

function lsRequestBody() {
  const fixed = {};
  LS_FIXED_FIELDS.forEach((f) => {
    const v = lsNum(`ls-fixed-${f}`);
    if (v !== null) fixed[f] = v;
  });
  const body = {
    x_axis: lsAxisSpec("ls-x"),
    y_axis: lsAxisSpec("ls-y"),
    fixed,
    metric: lsEl("ls-metric").value,
    constraints: lsConstraints(),
  };
  if (!lsTemperatureIsAxis()) body.temperature_K = lsNum("ls-temp");
  const thermo = lsThermoBlock();
  if (Object.keys(thermo).length) body.thermo = thermo;
  return body;
}

function lsTemperatureIsAxis() {
  return lsEl("ls-x-param").value === "temperature_K" || lsEl("ls-y-param").value === "temperature_K";
}

// ---------- colour scale ----------

function lsColour(t) {
  if (!Number.isFinite(t)) return null;
  const clamped = Math.min(1, Math.max(0, t));
  const pos = clamped * (LS_VIRIDIS.length - 1);
  const i = Math.min(LS_VIRIDIS.length - 2, Math.floor(pos));
  const f = pos - i;
  const a = LS_VIRIDIS[i], b = LS_VIRIDIS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

// log colour scale when the data is positive and spans decades
function lsScaleFn() {
  const { min, max } = lsState.finite;
  const useLog = min > 0 && max / min >= 100;
  if (useLog) {
    const lo = Math.log10(min), hi = Math.log10(max);
    return { fn: (v) => (Math.log10(v) - lo) / (hi - lo), log: true };
  }
  return { fn: (v) => (v - min) / (max - min || 1), log: false };
}

// ---------- axis transforms ----------

function lsAxisMap(axis, px0, px1) {
  const vals = axis.values;
  const lo = vals[0], hi = vals[vals.length - 1];
  const log = axis.scale === "log";
  const t = (v) => (log ? (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo))
                        : (v - lo) / (hi - lo));
  const inv = (u) => (log ? Math.pow(10, Math.log10(lo) + u * (Math.log10(hi) - Math.log10(lo)))
                          : lo + u * (hi - lo));
  return {
    toPx: (v) => px0 + t(v) * (px1 - px0),
    fromPx: (px) => inv((px - px0) / (px1 - px0)),
    lo, hi, log,
  };
}

function lsNiceLinearTicks(min, max, n) {
  const span = max - min;
  if (!(span > 0)) return [min];
  const raw = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= n) || mag * 10;
  const ticks = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) ticks.push(v);
  return ticks;
}

function lsAxisTicks(axis) {
  const lo = axis.values[0], hi = axis.values[axis.values.length - 1];
  return axis.scale === "log" ? niceLogTicks(lo, hi).filter((t) => t >= lo * 0.999 && t <= hi * 1.001)
                              : lsNiceLinearTicks(lo, hi, 6);
}

function lsFormatValue(param, v) {
  if (param === "temperature_K") return `${lsHumanNumber(v, 5)} K (${lsHumanNumber(v - 273.15, 3)} °C)`;
  if (param === "k_ck" || param === "luckey" || param === "k_target" || param === "target_conc") {
    return formatConcNm(v);
  }
  return lsHumanNumber(v, 3);
}

function lsFormatMetric(metric, v) {
  if (v === null || v === undefined || !Number.isFinite(v)) return "undefined";
  if (metric === "limit_of_detection") return formatConcNm(v);
  if (metric === "fold_change") return formatFoldChange(v);
  return lsHumanNumber(v, 4);
}

// ---------- drawing ----------

function lsDraw() {
  const data = lsState.data;
  const canvas = lsEl("ls-canvas");
  if (!data) return;

  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 640;
  const cssH = 440;
  canvas.width = cssW * dpr;
  canvas.height = cssH * dpr;
  canvas.style.height = `${cssH}px`;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  const plot = {
    x0: LS_PAD.left, x1: cssW - LS_PAD.right,
    y0: LS_PAD.top, y1: cssH - LS_PAD.bottom,
  };
  const mx = lsAxisMap(data.x, plot.x0, plot.x1);
  const my = lsAxisMap(data.y, plot.y1, plot.y0); // y grows upward on screen

  lsDrawHeatmap(ctx, plot);
  if (data.feasible) lsDrawFeasible(ctx, plot, mx, my);
  if (data.tm_band) lsDrawTmBand(ctx, plot, mx, my);
  if (lsEl("ls-show-contours").checked) lsDrawContours(ctx, plot, mx, my);
  lsDrawFrameAndTicks(ctx, plot, mx, my);
  lsDrawMarkers(ctx, plot, mx, my);
  lsDrawColourbar(ctx, plot, cssW);
}

function lsDrawHeatmap(ctx, plot) {
  const { z } = lsState.data;
  const nx = z.length, ny = z[0].length;
  const { fn } = lsScaleFn();

  // paint at grid resolution, upscale unsmoothed
  const off = document.createElement("canvas");
  off.width = nx;
  off.height = ny;
  const octx = off.getContext("2d");
  const img = octx.createImageData(nx, ny);
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      const rgb = z[i][j] === null ? null : lsColour(fn(z[i][j]));
      const px = ((ny - 1 - j) * nx + i) * 4;
      if (rgb) {
        img.data[px] = rgb[0];
        img.data[px + 1] = rgb[1];
        img.data[px + 2] = rgb[2];
        img.data[px + 3] = 255;
      } else {
        img.data[px + 3] = 0; // undefined cell: leave the background showing through
      }
    }
  }
  octx.putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(off, plot.x0, plot.y0, plot.x1 - plot.x0, plot.y1 - plot.y0);
}

function lsHatchPattern(ctx) {
  const tile = document.createElement("canvas");
  tile.width = 8;
  tile.height = 8;
  const tctx = tile.getContext("2d");
  tctx.strokeStyle = "rgba(255,255,255,0.85)";
  tctx.lineWidth = 1.4;
  tctx.beginPath();
  tctx.moveTo(-2, 10); tctx.lineTo(10, -2);
  tctx.stroke();
  return ctx.createPattern(tile, "repeat");
}

// hatch + outline, not a second colour fill
function lsDrawFeasible(ctx, plot, mx, my) {
  const { feasible, x, y } = lsState.data;
  const nx = x.values.length, ny = y.values.length;
  const cw = (plot.x1 - plot.x0) / nx;
  const ch = (plot.y1 - plot.y0) / ny;
  const cellX = (i) => plot.x0 + i * cw;
  const cellY = (j) => plot.y1 - (j + 1) * ch;

  ctx.save();
  ctx.fillStyle = lsHatchPattern(ctx);
  ctx.beginPath();
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      if (feasible[i][j]) ctx.rect(cellX(i), cellY(j), cw + 0.5, ch + 0.5);
    }
  }
  ctx.fill();

  ctx.strokeStyle = "#ffffff";
  ctx.lineWidth = 2;
  ctx.beginPath();
  const inside = (i, j) => i >= 0 && j >= 0 && i < nx && j < ny && feasible[i][j];
  for (let i = 0; i < nx; i++) {
    for (let j = 0; j < ny; j++) {
      if (!feasible[i][j]) continue;
      const left = cellX(i), right = cellX(i) + cw, top = cellY(j), bottom = cellY(j) + ch;
      if (!inside(i - 1, j)) { ctx.moveTo(left, top); ctx.lineTo(left, bottom); }
      if (!inside(i + 1, j)) { ctx.moveTo(right, top); ctx.lineTo(right, bottom); }
      if (!inside(i, j + 1)) { ctx.moveTo(left, top); ctx.lineTo(right, top); }
      if (!inside(i, j - 1)) { ctx.moveTo(left, bottom); ctx.lineTo(right, bottom); }
    }
  }
  ctx.stroke();
  ctx.restore();
}

function lsDrawTmBand(ctx, plot, mx, my) {
  const band = lsState.data.tm_band;
  const map = band.axis === "x" ? mx : my;
  const a = map.toPx(band.min), b = map.toPx(band.max);
  ctx.save();
  ctx.fillStyle = "rgba(180, 83, 9, 0.38)";
  if (band.axis === "x") {
    ctx.fillRect(Math.min(a, b), plot.y0, Math.abs(b - a), plot.y1 - plot.y0);
  } else {
    ctx.fillRect(plot.x0, Math.min(a, b), plot.x1 - plot.x0, Math.abs(b - a));
  }
  ctx.fillStyle = "#fff";
  ctx.font = "600 11px system-ui, sans-serif";
  const label = "model invalid near Tm";
  if (band.axis === "x") {
    ctx.translate(Math.min(a, b) + Math.abs(b - a) / 2, plot.y0 + 10);
    ctx.rotate(Math.PI / 2);
    ctx.fillText(label, 6, 4);
  } else {
    ctx.fillText(label, plot.x0 + 8, Math.min(a, b) + 14);
  }
  ctx.restore();
}

function lsDrawContours(ctx, plot, mx, my) {
  const { contours, metric } = lsState.data;
  ctx.save();
  ctx.beginPath();
  ctx.rect(plot.x0, plot.y0, plot.x1 - plot.x0, plot.y1 - plot.y0);
  ctx.clip();
  ctx.lineWidth = 1.5;
  ctx.strokeStyle = "rgba(255,255,255,0.9)";
  ctx.font = "600 11px system-ui, sans-serif";

  contours.forEach((contour) => {
    let longest = null;
    contour.polylines.forEach((pl) => {
      ctx.beginPath();
      pl.forEach((p, k) => {
        const px = mx.toPx(p[0]), py = my.toPx(p[1]);
        if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      });
      ctx.stroke();
      if (!longest || pl.length > longest.length) longest = pl;
    });
    if (longest && longest.length > 3) {
      const mid = longest[Math.floor(longest.length / 2)];
      const px = mx.toPx(mid[0]), py = my.toPx(mid[1]);
      const text = lsFormatMetric(metric, contour.level);
      const w = ctx.measureText(text).width;
      ctx.fillStyle = "rgba(12,31,61,0.75)";
      ctx.fillRect(px - w / 2 - 3, py - 8, w + 6, 14);
      ctx.fillStyle = "#fff";
      ctx.textAlign = "center";
      ctx.fillText(text, px, py + 3);
      ctx.textAlign = "left";
    }
  });
  ctx.restore();
}

function lsDrawFrameAndTicks(ctx, plot, mx, my) {
  const { x, y } = lsState.data;
  ctx.save();
  ctx.strokeStyle = "#8c9aac";
  ctx.lineWidth = 1;
  ctx.strokeRect(plot.x0, plot.y0, plot.x1 - plot.x0, plot.y1 - plot.y0);
  ctx.fillStyle = "#5c6b7d";
  ctx.font = "11px system-ui, sans-serif";

  ctx.textAlign = "center";
  lsAxisTicks(x).forEach((t) => {
    const px = mx.toPx(t);
    if (px < plot.x0 - 1 || px > plot.x1 + 1) return;
    ctx.beginPath();
    ctx.moveTo(px, plot.y1);
    ctx.lineTo(px, plot.y1 + 4);
    ctx.stroke();
    ctx.fillText(lsHumanNumber(t, 3), px, plot.y1 + 17);
  });
  ctx.fillText(lsAxisTitle(x), (plot.x0 + plot.x1) / 2, plot.y1 + 38);

  ctx.textAlign = "right";
  lsAxisTicks(y).forEach((t) => {
    const py = my.toPx(t);
    if (py < plot.y0 - 1 || py > plot.y1 + 1) return;
    ctx.beginPath();
    ctx.moveTo(plot.x0 - 4, py);
    ctx.lineTo(plot.x0, py);
    ctx.stroke();
    ctx.fillText(lsHumanNumber(t, 3), plot.x0 - 7, py + 4);
  });
  ctx.save();
  ctx.translate(14, (plot.y0 + plot.y1) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText(lsAxisTitle(y), 0, 0);
  ctx.restore();
  ctx.restore();
}

function lsMarker(ctx, px, py, colour, label, plot, dy = 0) {
  ctx.save();
  ctx.lineWidth = 2.5;
  ctx.strokeStyle = "#0c1f3d";
  ctx.fillStyle = colour;
  ctx.beginPath();
  ctx.arc(px, py, 6, 0, Math.PI * 2);
  ctx.fill();
  ctx.stroke();
  if (label) {
    ctx.font = "600 11px system-ui, sans-serif";
    const w = ctx.measureText(label).width;
    // flip the label to the left when it would run past the plot, and nudge coincident ones apart
    const right = px + 9 + w + 6;
    const bx = plot && right > plot.x1 ? px - 9 - (w + 6) : px + 9;
    const by = py + dy;
    ctx.fillStyle = "rgba(255,255,255,0.92)";
    ctx.fillRect(bx, by - 8, w + 6, 15);
    ctx.fillStyle = "#0c1f3d";
    ctx.fillText(label, bx + 3, by + 3);
  }
  ctx.restore();
}

function lsDrawMarkers(ctx, plot, mx, my) {
  const data = lsState.data;
  const within = (px, py) => px >= plot.x0 && px <= plot.x1 && py >= plot.y0 && py <= plot.y1;

  const current = lsCurrentDesign();
  if (current) {
    const px = mx.toPx(current.x), py = my.toPx(current.y);
    if (within(px, py)) lsMarker(ctx, px, py, "#ffffff", "your design", plot, -14);
  }

  if (data.robust_optimum) {
    const px = mx.toPx(data.robust_optimum.x), py = my.toPx(data.robust_optimum.y);
    if (within(px, py)) lsMarker(ctx, px, py, "#15803d", "safest", plot, 0);
  }

  (lsState.inverse || []).forEach((sol) => {
    const px = mx.toPx(sol.x), py = my.toPx(sol.y);
    if (within(px, py)) lsMarker(ctx, px, py, "#c01919", sol.label, plot, 14);
  });

  if (lsState.pin) {
    const px = mx.toPx(data.x.values[lsState.pin.i]);
    const py = my.toPx(data.y.values[lsState.pin.j]);
    if (within(px, py)) lsMarker(ctx, px, py, "#2f6fb0", "pinned", plot, 16);
  }
}

function lsDrawColourbar(ctx, plot, cssW) {
  const { metric } = lsState.data;
  const { fn, log } = lsScaleFn();
  const { min, max } = lsState.finite;
  const barX = cssW - LS_PAD.right + 14;
  const barW = 14;
  const h = plot.y1 - plot.y0;

  for (let p = 0; p < h; p++) {
    const rgb = lsColour(1 - p / (h - 1));
    ctx.fillStyle = `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
    ctx.fillRect(barX, plot.y0 + p, barW, 1);
  }
  ctx.strokeStyle = "#8c9aac";
  ctx.strokeRect(barX, plot.y0, barW, h);

  ctx.fillStyle = "#5c6b7d";
  ctx.font = "10px system-ui, sans-serif";
  ctx.textAlign = "left";
  const stops = log
    ? [min, Math.sqrt(min * max), max]
    : [min, (min + max) / 2, max];
  stops.forEach((v) => {
    const py = plot.y0 + (1 - fn(v)) * h;
    ctx.fillText(lsFormatMetric(metric, v), barX + barW + 4, Math.min(plot.y1, Math.max(plot.y0 + 8, py + 3)));
  });
  ctx.save();
  ctx.translate(barX + barW + 46, (plot.y0 + plot.y1) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText(LS_METRIC_LABELS[metric] + (log ? " (log scale)" : ""), 0, 0);
  ctx.restore();
}

// ---------- markers from the other modules ----------

function lsCurrentDesign() {
  const data = lsState.data;
  if (!data) return null;
  const read = (param) => {
    const v = lsNum(`ls-fixed-${param}`);
    return param === "temperature_K" ? lsNum("ls-temp") : v;
  };
  const x = read(data.x.param), y = read(data.y.param);
  return x === null || y === null ? null : { x, y };
}

// K = exp(dG / RT), matching temperature.py's k_from_dG
function lsFromModuleC(solution, param, temperature_K) {
  const RT = LS_R_KCAL * temperature_K;
  const v = solution.values;
  if (param === "k_open" && v.dG_latch_cage !== undefined) return Math.exp(v.dG_latch_cage / RT);
  if (param === "k_ck" && v.dG_key_cage !== undefined) return Math.exp(v.dG_key_cage / RT) / 1e-9;
  if (param === "luckey" && v.key_conc !== undefined) return v.key_conc;
  return null;
}

// ---------- readouts ----------

function lsCellAt(evt) {
  const data = lsState.data;
  if (!data) return null;
  const canvas = lsEl("ls-canvas");
  const rect = canvas.getBoundingClientRect();
  const px = evt.clientX - rect.left, py = evt.clientY - rect.top;
  const plot = {
    x0: LS_PAD.left, x1: rect.width - LS_PAD.right,
    y0: LS_PAD.top, y1: rect.height - LS_PAD.bottom,
  };
  if (px < plot.x0 || px > plot.x1 || py < plot.y0 || py > plot.y1) return null;
  const nx = data.x.values.length, ny = data.y.values.length;
  const i = Math.min(nx - 1, Math.floor(((px - plot.x0) / (plot.x1 - plot.x0)) * nx));
  const j = Math.min(ny - 1, Math.floor((1 - (py - plot.y0) / (plot.y1 - plot.y0)) * ny));
  return { i, j };
}

function lsReadoutHtml(cell, title) {
  const d = lsState.data;
  const value = d.z[cell.i][cell.j];
  const feas = d.feasible ? (d.feasible[cell.i][cell.j] ? "passes your rules" : "fails your rules") : null;
  return `<div class="ls-readout-title">${title}</div>
    <div class="ls-readout-row"><span>${LS_LABELS[d.x.param]}</span><b>${lsFormatValue(d.x.param, d.x.values[cell.i])}</b></div>
    <div class="ls-readout-row"><span>${LS_LABELS[d.y.param]}</span><b>${lsFormatValue(d.y.param, d.y.values[cell.j])}</b></div>
    <div class="ls-readout-row"><span>${LS_METRIC_LABELS[d.metric]}</span><b>${lsFormatMetric(d.metric, value)}</b></div>
    ${feas ? `<div class="ls-readout-row"><span>status</span><b>${feas}</b></div>` : ""}`;
}

function lsRenderReadouts() {
  const hover = lsEl("ls-hover-readout");
  hover.innerHTML = lsState.hover
    ? lsReadoutHtml(lsState.hover, "At cursor")
    : `<div class="ls-readout-title">At cursor</div><div class="help-text">Move the pointer over the plot.</div>`;

  const pin = lsEl("ls-pin-readout");
  pin.innerHTML = lsState.pin
    ? lsReadoutHtml(lsState.pin, "Pinned")
    : `<div class="ls-readout-title">Pinned</div><div class="help-text">Click the plot to pin a design here.</div>`;
}

function lsRenderRuleSummary() {
  const rules = lsConstraints();
  const change = ' <button type="button" data-ls-open-advanced>Change</button>';
  const phrase = rules.map((r) =>
    `<b>${LS_METRIC_LABELS[r.metric].toLowerCase()} ${LS_OP_LABELS[r.op]} ` +
    `${lsHumanNumber(r.value)}${LS_METRIC_UNITS[r.metric] || ""}</b>`).join(" and ");
  lsEl("ls-rule-summary").innerHTML = rules.length
    ? `Striped areas on the map are designs with ${phrase}.${change}`
    : `No rule set, so nothing is striped on the map.${change}`;
}

function lsRenderWarnings() {
  const el = lsEl("ls-warnings");
  el.innerHTML = "";
  const raw = (lsState.data ? lsState.data.warnings : []).slice();
  // The engine warns once per unscaled energy term in API language; a reader needs one sentence.
  const unscaled = raw.filter((w) => w.startsWith("TEMPERATURE SCALING UNAVAILABLE"));
  const notes = raw.filter((w) => !unscaled.includes(w));
  if (unscaled.length) {
    notes.push("This map was worked out at room temperature. Changing the temperature needs the " +
               "heat-response numbers under Advanced, so the plot ignored the temperature setting.");
  }
  if (lsState.presetNote) notes.push(lsState.presetNote);
  if (lsState.inverseNote) notes.push(lsState.inverseNote);
  notes.forEach((w) => {
    const div = document.createElement("div");
    div.className = "field-error";
    div.textContent = w;
    el.appendChild(div);
  });
}

function lsRenderRobust() {
  const el = lsEl("ls-robust");
  const ro = lsState.data && lsState.data.robust_optimum;
  if (!ro) {
    el.style.display = "none";
    return;
  }
  el.style.display = "block";
  // The API's own `summary` is written for API consumers; this says the same thing in words.
  const d = lsState.data;
  const line = (param, value, slack) => {
    const unit = lsUnitSuffix(param === d.x.param ? d.x.unit : d.y.unit);
    const low = lsHumanNumber(value - slack.minus), high = lsHumanNumber(value + slack.plus);
    return `<li><b>${LS_LABELS[param]}</b> of <b>${lsHumanNumber(value)}</b>${unit}
      &mdash; anything from ${low} to ${high} still passes, so you have room to be off.</li>`;
  };
  el.innerHTML = `<div class="verdict-title">Safest design on this map</div>
    <div>This is the design furthest from the edge of the striped area, so it is the one most
    likely to still work if your real protein comes out a bit different from the prediction.
    It gives ${lsFormatMetric(d.metric, ro.metric_value)}.</div>
    <ul class="recommendations" style="margin-top:8px;">
      ${line(d.x.param, ro.x, ro.x_slack)}
      ${line(d.y.param, ro.y, ro.y_slack)}
    </ul>`;
}

// ---------- actions ----------

async function lsPlot() {
  const btn = lsEl("ls-submit");
  btn.disabled = true;
  btn.textContent = "Plotting...";
  lsEl("ls-error").style.display = "none";
  try {
    const data = await apiPost("/api/landscape", lsRequestBody());
    lsState.data = data;
    lsState.pin = null;
    lsState.hover = null;
    // Reduced rather than spread: a capped grid is 40k cells, near the argument limit.
    let min = Infinity, max = -Infinity;
    data.z.forEach((row) => row.forEach((v) => {
      if (v === null || !Number.isFinite(v)) return;
      if (v < min) min = v;
      if (v > max) max = v;
    }));
    if (min === Infinity) throw { message: "Every cell in this range is undefined. Widen or shift the axis ranges." };
    lsState.finite = { min, max };
    const lowGood = data.metric === "off_signal" || data.metric === "limit_of_detection";
    lsEl("ls-plot-caption").textContent =
      `Every square is one possible sensor. ${lowGood ? "Dark purple" : "Bright yellow"} squares are ` +
      `the good ones, ${lowGood ? "bright yellow" : "dark purple"} are the bad ones. The striped patch ` +
      `is every sensor that passes your rules, and the green dot inside it is the one to build.`;
    lsEl("ls-empty-state").style.display = "none";
    lsEl("ls-plot-wrap").style.display = "block";
    lsDraw();
    lsRenderReadouts();
    lsRenderWarnings();
    lsRenderRobust();
  } catch (err) {
    const el = lsEl("ls-error");
    el.style.display = "block";
    el.textContent = err.message || "Could not plot the landscape.";
  } finally {
    btn.disabled = false;
    btn.textContent = "Draw the map";
  }
}

async function lsLoadInverseSolutions() {
  const data = lsState.data;
  if (!data) return;
  const temperature_K = data.temperature_K || 298.15;
  const solvable = ["k_open", "k_ck", "luckey"];
  if (!solvable.includes(data.x.param) || !solvable.includes(data.y.param)) {
    lsState.inverse = [];
    lsState.inverseNote = `Suggestions can only be marked when both axes are cage looseness, key ` +
      `grip strength or key amount. This map uses ${LS_LABELS[data.x.param].toLowerCase()} and ` +
      `${LS_LABELS[data.y.param].toLowerCase()}, so there is nothing to mark.`;
    lsRenderWarnings();
    return;
  }
  const toDg = (k) => LS_R_KCAL * temperature_K * Math.log(k);
  const cParam = { k_open: "dG_latch_cage", k_ck: "dG_key_cage", luckey: "key_conc" };
  const axisBounds = (axis) => {
    const lo = axis.values[0], hi = axis.values[axis.values.length - 1];
    if (axis.param === "luckey") return { min: lo, max: hi };
    const scale = axis.param === "k_ck" ? 1e-9 : 1;
    return { min: toDg(lo * scale), max: toDg(hi * scale) };
  };

  const body = {
    solve_for: [cParam[data.x.param], cParam[data.y.param]],
    fixed: {
      dG_latch_cage: toDg(lsNum("ls-fixed-k_open")),
      dG_key_cage: toDg(lsNum("ls-fixed-k_ck") * 1e-9),
      key_conc: lsNum("ls-fixed-luckey"),
    },
    bounds: { [cParam[data.x.param]]: axisBounds(data.x), [cParam[data.y.param]]: axisBounds(data.y) },
    target_fold_change: Number(lsEl("ls-inverse-target").value),
    temperature_K,
  };
  delete body.fixed[cParam[data.x.param]];
  delete body.fixed[cParam[data.y.param]];

  try {
    const result = await apiPost("/api/inverse_design", body);
    const points = (result.robust_point ? [result.robust_point] : result.solutions || [])
      .map((s) => ({
        x: lsFromModuleC(s, data.x.param, temperature_K),
        y: lsFromModuleC(s, data.y.param, temperature_K),
        label: `Module C ${formatFoldChange(s.fold_change)}`,
      }))
      .filter((p) => p.x !== null && p.y !== null);
    lsState.inverse = points;
    lsState.inverseNote = points.length ? "" : (result.message || "Module C returned no placeable solution.");
    lsDraw();
    lsRenderWarnings();
  } catch (err) {
    lsState.inverse = [];
    lsState.inverseNote = `Module C: ${err.message || "inverse design failed."}`;
    lsRenderWarnings();
  }
}

function lsExportPng() {
  const link = document.createElement("a");
  link.download = `latte-landscape-${lsState.data.metric}.png`;
  link.href = lsEl("ls-canvas").toDataURL("image/png");
  link.click();
}

function lsExportCsv() {
  const d = lsState.data;
  const header = [`${d.x.param}_${d.x.unit}`, `${d.y.param}_${d.y.unit}`, d.metric];
  if (d.feasible) header.push("feasible");
  const rows = [header.join(",")];
  d.x.values.forEach((xv, i) => {
    d.y.values.forEach((yv, j) => {
      const row = [xv, yv, d.z[i][j] === null ? "" : d.z[i][j]];
      if (d.feasible) row.push(d.feasible[i][j] ? "1" : "0");
      rows.push(row.join(","));
    });
  });
  const blob = new Blob([rows.join("\n")], { type: "text/csv" });
  const link = document.createElement("a");
  link.download = `latte-landscape-${d.metric}.csv`;
  link.href = URL.createObjectURL(blob);
  link.click();
  URL.revokeObjectURL(link.href);
}

function lsAddConstraintRow(metric = "fold_change", op = ">=", value = "") {
  const wrap = lsEl("ls-constraints");
  const row = document.createElement("div");
  row.className = "ls-constraint";
  row.innerHTML = `
    <select data-c-metric>${LS_METRICS.map((m) => `<option value="${m}"${m === metric ? " selected" : ""}>${LS_METRIC_LABELS[m]}</option>`).join("")}</select>
    <select data-c-op>${Object.keys(LS_OP_LABELS).map((o) => `<option value="${o}"${o === op ? " selected" : ""}>${LS_OP_LABELS[o]}</option>`).join("")}</select>
    <input type="number" step="any" data-c-value value="${value}" placeholder="value">
    <button type="button" class="btn-ghost" data-c-remove title="remove this rule">&times;</button>`;
  row.querySelector("[data-c-remove]").addEventListener("click", () => {
    row.remove();
    lsRenderRuleSummary();
  });
  wrap.appendChild(row);
}


// ---------- presets ----------

// one question per preset, fills every axis choice
const LS_PRESETS = {
  "cage-vs-pull": {
    explain: "Cage looseness across the bottom, how hard the target pulls up the side. " +
             "Bright means a big brightness jump. Shows the sweet spot between a cage so tight " +
             "nothing opens it and one so loose it glows on its own.",
    x: { param: "k_open", min: 1e-5, max: 1e-1, steps: 80, scale: "log" },
    y: { param: "pull", min: 1, max: 30, steps: 60, scale: "linear" },
    metric: "fold_change",
    rule: { metric: "fold_change", op: ">=", value: 20 },
  },
  key: {
    explain: "How much glowing key you add, against how tightly it grips. Bright means a big " +
             "brightness jump. Shows whether buying a tighter binder beats just adding more key.",
    x: { param: "luckey", min: 10, max: 5000, steps: 80, scale: "log" },
    y: { param: "k_ck", min: 0.1, max: 100, steps: 60, scale: "log" },
    metric: "fold_change",
    // pull is fixed at 10 in this view, which caps the jump near 10x; 20 would show nothing.
    rule: { metric: "fold_change", op: ">=", value: 8 },
  },
  temperature: {
    explain: "Temperature across the bottom, cage looseness up the side. 277 K is a fridge, " +
             "298 K is a benchtop, 310 K is body temperature. Needs the heat-response numbers " +
             "under Advanced, which this has filled with placeholders you must replace.",
    x: { param: "temperature_K", min: 277.15, max: 330, steps: 70, scale: "linear" },
    y: { param: "k_open", min: 1e-5, max: 1e-1, steps: 55, scale: "log" },
    metric: "fold_change",
    rule: { metric: "fold_change", op: ">=", value: 8 },
    thermo: { "ls-latch-dh": -12, "ls-latch-dcp": 0.4, "ls-latch-tm": 335, "ls-key-dh": -8 },
  },
  detect: {
    explain: "The smallest amount of target the sensor can notice, so here DARK is better. " +
             "Needs a target grip strength, which this has set to 1 nM under Advanced.",
    x: { param: "k_open", min: 1e-5, max: 1e-1, steps: 60, scale: "log" },
    y: { param: "pull", min: 1, max: 30, steps: 50, scale: "linear" },
    metric: "limit_of_detection",
    rule: { metric: "limit_of_detection", op: "<=", value: 1 },
    fixed: { "ls-fixed-k_target": 1 },
  },
  custom: { explain: "Pick the axes yourself under Advanced below.", custom: true },
};

const LS_PLACEHOLDER_THERMO_NOTE =
  "The heat-response numbers (dH, Tm) were filled in with textbook placeholders so this view " +
  "would draw something. They are NOT your construct's measured values. Replace them under " +
  "Advanced before trusting anything this plot says about temperature.";

function lsApplyPreset(name) {
  const preset = LS_PRESETS[name];
  if (!preset) return;
  lsEl("ls-preset-explain").textContent = preset.explain;
  lsState.presetNote = "";
  if (preset.custom) {
    lsEl("ls-advanced").open = true;
    return;
  }

  const setAxis = (prefix, axis) => {
    lsEl(`${prefix}-param`).value = axis.param;
    lsEl(`${prefix}-min`).value = axis.min;
    lsEl(`${prefix}-max`).value = axis.max;
    lsEl(`${prefix}-steps`).value = axis.steps;
    const scale = lsEl(`${prefix}-scale`);
    scale.dataset.scale = axis.scale;
    scale.querySelectorAll("button").forEach((b) =>
      b.classList.toggle("active", b.dataset.scaleValue === axis.scale));
  };
  setAxis("ls-x", preset.x);
  setAxis("ls-y", preset.y);
  lsEl("ls-metric").value = preset.metric;

  if (preset.rule) {
    lsEl("ls-constraints").innerHTML = "";
    lsAddConstraintRow(preset.rule.metric, preset.rule.op, preset.rule.value);
  }
  Object.entries(preset.fixed || {}).forEach(([id, v]) => {
    if (!lsEl(id).value) lsEl(id).value = v;
  });
  if (preset.thermo) {
    // this question cannot be answered without these, so put them in front of the user
    lsEl("ls-advanced").open = true;
    lsEl("ls-temp-section").open = true;
    Object.entries(preset.thermo).forEach(([id, v]) => {
      if (!lsEl(id).value) lsEl(id).value = v;
    });
    // warn while the placeholders are still there, not just on first pick
    const untouched = Object.entries(preset.thermo)
      .every(([id, v]) => Number(lsEl(id).value) === v);
    lsState.presetNote = untouched ? LS_PLACEHOLDER_THERMO_NOTE : "";
  }
  lsSyncTemperatureControl();
  lsRenderRuleSummary();
  lsRenderWarnings();
}

// The slider only means anything with a heat response to scale by; without one it would
// silently hand back the 298 K answer under a different label, so it stays locked at 298.
function lsSyncTemperatureControl() {
  const isAxis = lsTemperatureIsAxis();
  const slider = lsEl("ls-temp");
  slider.disabled = isAxis || lsNum("ls-latch-dh") === null;
  if (slider.disabled && !isAxis) slider.value = 298.15;
  lsEl("ls-temp-note").textContent = isAxis
    ? "Temperature is on an axis, so the slider is off."
    : slider.disabled
      ? "Locked at room temperature. To change it, fill in the cage heat response below."
      : "";
  lsEl("ls-temp-value").textContent = `${slider.value} K (${roundSig(Number(slider.value) - 273.15, 3)} °C)`;
}

function initLandscape() {
  if (!lsEl("ls-canvas")) return;

  ["ls-x", "ls-y"].forEach((prefix) => {
    const select = lsEl(`${prefix}-param`);
    select.innerHTML = LS_AXIS_PARAMS
      .map((p) => `<option value="${p}">${LS_LABELS[p]}</option>`).join("");
    select.addEventListener("change", lsSyncTemperatureControl);
    const scale = lsEl(`${prefix}-scale`);
    scale.querySelectorAll("button").forEach((btn) => {
      btn.addEventListener("click", () => {
        scale.dataset.scale = btn.dataset.scaleValue;
        scale.querySelectorAll("button").forEach((b) => b.classList.toggle("active", b === btn));
      });
    });
  });
  lsEl("ls-metric").innerHTML = LS_METRICS
    .map((m) => `<option value="${m}">${LS_METRIC_LABELS[m]}</option>`).join("");
  lsEl("ls-preset").addEventListener("change", (e) => {
    lsApplyPreset(e.target.value);
    if (lsState.data) lsPlot();   // a new question means a new map, not a new form to resubmit
  });

  lsEl("ls-submit").addEventListener("click", lsPlot);
  lsEl("ls-add-constraint").addEventListener("click", () => {
    lsAddConstraintRow();
    lsRenderRuleSummary();
  });
  lsEl("ls-constraints").addEventListener("input", lsRenderRuleSummary);
  lsEl("ls-constraints").addEventListener("change", lsRenderRuleSummary);
  lsEl("ls-rule-summary").addEventListener("click", (e) => {
    if (e.target.closest("[data-ls-open-advanced]")) lsEl("ls-advanced").open = true;
  });
  lsEl("ls-show-contours").addEventListener("change", () => { if (lsState.data) lsDraw(); });
  lsEl("ls-export-png").addEventListener("click", lsExportPng);
  lsEl("ls-export-csv").addEventListener("click", lsExportCsv);
  lsEl("ls-show-inverse").addEventListener("click", lsLoadInverseSolutions);
  lsEl("ls-thermo-block").addEventListener("input", lsSyncTemperatureControl);
  lsEl("ls-temp").addEventListener("input", () => {
    lsSyncTemperatureControl();
    if (lsState.data && !lsTemperatureIsAxis()) lsPlot();
  });

  const canvas = lsEl("ls-canvas");
  canvas.addEventListener("mousemove", (evt) => {
    const cell = lsCellAt(evt);
    lsState.hover = cell;
    lsRenderReadouts();
  });
  canvas.addEventListener("mouseleave", () => {
    lsState.hover = null;
    lsRenderReadouts();
  });
  canvas.addEventListener("click", (evt) => {
    const cell = lsCellAt(evt);
    if (!cell) return;
    lsState.pin = cell;
    lsDraw();
    lsRenderReadouts();
  });

  window.addEventListener("resize", () => { if (lsState.data) lsDraw(); });
  document.addEventListener("tabchange", (e) => {
    if (e.detail.tab !== "landscape") return;
    if (lsState.data) lsDraw(); else lsPlot();
  });

  lsApplyPreset(lsEl("ls-preset").value);
}

document.addEventListener("DOMContentLoaded", initLandscape);
