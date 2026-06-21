
function niceLogTicks(min, max) {
  const lo = Math.floor(Math.log10(min));
  const hi = Math.ceil(Math.log10(max));
  const ticks = [];
  for (let e = lo; e <= hi; e++) ticks.push(Math.pow(10, e));
  return ticks;
}

function formatTick(x, isConcentration) {
  if (!isConcentration) return roundSig(x, 2);
  return formatConcNm(x);
}

// lines: [{points:[{x,y}], color, dashed}]; marker: {x,y,color,label}
function drawLogXPlot(svg, { lines, marker, xIsConcentration, yLabel }) {
  svg.innerHTML = "";
  const width = svg.clientWidth || 480;
  const height = 260;
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  const padding = { left: 56, right: 16, top: 16, bottom: 36 };
  const plotW = width - padding.left - padding.right;
  const plotH = height - padding.top - padding.bottom;

  const allX = lines.flatMap((l) => l.points.map((p) => p.x)).filter((x) => x > 0);
  const allY = lines.flatMap((l) => l.points.map((p) => p.y));
  if (marker) { allX.push(marker.x); allY.push(marker.y); }
  const xMin = Math.min(...allX), xMax = Math.max(...allX);
  // y always starts at 1 (no change), never 0 -- per spec 7.1.
  const yMin = 1;
  const yMax = Math.max(...allY) * 1.1;

  const xToPx = (x) => padding.left + ((Math.log10(x) - Math.log10(xMin)) / (Math.log10(xMax) - Math.log10(xMin))) * plotW;
  const yToPx = (y) => padding.top + plotH - ((y - yMin) / (yMax - yMin)) * plotH;

  const ns = "http://www.w3.org/2000/svg";
  const make = (tag, attrs) => {
    const el = document.createElementNS(ns, tag);
    Object.entries(attrs).forEach(([k, v]) => el.setAttribute(k, v));
    return el;
  };

  // gridlines + x ticks (log-spaced decades, unit-scaled labels)
  niceLogTicks(xMin, xMax).forEach((tick) => {
    const px = xToPx(tick);
    svg.appendChild(make("line", { x1: px, x2: px, y1: padding.top, y2: padding.top + plotH, stroke: "var(--plot-grid)" }));
    const label = make("text", { x: px, y: height - 12, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "middle" });
    label.textContent = formatTick(tick, xIsConcentration);
    svg.appendChild(label);
  });

  // y axis (3 ticks: min, mid, max)
  [yMin, (yMin + yMax) / 2, yMax].forEach((tick) => {
    const py = yToPx(tick);
    svg.appendChild(make("line", { x1: padding.left, x2: width - padding.right, y1: py, y2: py, stroke: "var(--plot-grid)" }));
    const label = make("text", { x: padding.left - 8, y: py + 4, "font-size": 11, fill: "var(--text-muted)", "text-anchor": "end" });
    label.textContent = roundSig(tick, 3);
    svg.appendChild(label);
  });

  lines.forEach((line) => {
    const d = line.points
      .filter((p) => p.x > 0)
      .map((p, i) => `${i === 0 ? "M" : "L"} ${xToPx(p.x)} ${yToPx(p.y)}`)
      .join(" ");
    svg.appendChild(make("path", {
      d, fill: "none", stroke: line.color, "stroke-width": 2,
      "stroke-dasharray": line.dashed ? "6,4" : "none",
    }));
  });

  if (marker) {
    svg.appendChild(make("circle", { cx: xToPx(marker.x), cy: yToPx(marker.y), r: 5, fill: marker.color || "var(--plot-marker)" }));
    const label = make("text", {
      x: xToPx(marker.x) + 8, y: yToPx(marker.y) - 8, "font-size": 12, fill: "var(--text)", "font-weight": 600,
    });
    label.textContent = marker.label || "";
    svg.appendChild(label);
  }

  if (yLabel) {
    const label = make("text", {
      x: 12, y: padding.top + 10, "font-size": 11, fill: "var(--text-muted)",
    });
    label.textContent = yLabel;
    svg.appendChild(label);
  }
}
