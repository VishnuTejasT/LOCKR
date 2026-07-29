

const STANDARD_AA = new Set("ACDEFGHIKLMNPQRSTVWY".split(""));
const ACIDIC_AA = new Set(["D", "E"]);
const BASIC_AA = new Set(["K", "R", "H"]);

let scanMode = "single";
const batchState = { results: [] };

function escHtml(s) {
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

// Parses raw sequences (one per line), FASTA, and mixed in one pass.
// Returns { records: [{id, sequence}], errors: [{lineNum, message}] }
function parseBatchInput(text) {
  const lines = text.split("\n");
  const records = [];
  const errors = [];
  let autoIdx = 0;
  let i = 0;

  while (i < lines.length) {
    const trimmed = lines[i].trim();
    if (trimmed === "") { i++; continue; }

    if (trimmed.startsWith(">")) {
      const headerLineNum = i + 1;
      const id = trimmed.slice(1).trim();
      i++;
      const seqParts = [];
      while (i < lines.length && !lines[i].trim().startsWith(">") && lines[i].trim() !== "") {
        seqParts.push(lines[i].trim().toUpperCase().replace(/\s/g, ""));
        i++;
      }
      const seq = seqParts.join("");
      if (!seq) {
        errors.push({ lineNum: headerLineNum, message: `FASTA record "${id || "(empty header)"}": no sequence found` });
      } else {
        const bad = [...seq].find(c => !STANDARD_AA.has(c));
        if (bad) {
          errors.push({ lineNum: headerLineNum, message: `FASTA record "${id}": non-standard amino acid '${bad}'` });
        } else {
          records.push({ id: id || `seq_${++autoIdx}`, sequence: seq });
        }
      }
    } else {
      const lineNum = i + 1;
      const seq = trimmed.toUpperCase().replace(/\s/g, "");
      const bad = [...seq].find(c => !STANDARD_AA.has(c));
      if (bad) {
        const preview = trimmed.length > 40 ? trimmed.slice(0, 40) + "…" : trimmed;
        errors.push({ lineNum, message: `Line ${lineNum}: non-standard character '${bad}' in "${preview}"` });
      } else {
        records.push({ id: `seq_${++autoIdx}`, sequence: seq });
      }
      i++;
    }
  }

  return { records, errors };
}

// One /scan request per record so each gets its own full-length window.
// Results are sorted ascending by liability_score (cleanest first).
async function batchScan(records, ph, policy) {
  const promises = records.map(rec =>
    apiPost("/scan", {
      sequences: [{ id: rec.id, sequence: rec.sequence }],
      sensitive_window: { start: 1, end: rec.sequence.length },
      ph,
      substitution_policy: policy,
      preserve_positions: [],
    }).then(data => data.results[0])
  );
  const results = await Promise.all(promises);
  results.sort((a, b) => a.liability_score - b.liability_score);
  return results;
}

function updateBatchCount() {
  const text = document.getElementById("scan-batch-text").value;
  const { records, errors } = parseBatchInput(text);

  document.getElementById("scan-batch-count").textContent =
    `${records.length} sequence${records.length !== 1 ? "s" : ""}`;

  const errEl = document.getElementById("scan-batch-parse-errors");
  if (errors.length > 0) {
    errEl.style.display = "block";
    errEl.className = "batch-parse-errors";
    errEl.innerHTML =
      `<div class="help-text" style="color:var(--warning-700); font-weight:600;">${errors.length} line${errors.length !== 1 ? "s" : ""} skipped:</div>` +
      errors.map(e => `<div class="help-text" style="color:var(--warning-700);">• ${escHtml(e.message)}</div>`).join("");
  } else {
    errEl.style.display = "none";
  }

  scanEl("scan-submit").disabled = records.length === 0;
  document.getElementById("scan-batch-error").textContent = "";
}

function renderBatchTable(results) {
  const wrap = document.getElementById("scan-batch-table-wrap");
  if (results.length === 0) { wrap.innerHTML = ""; return; }

  const table = document.createElement("table");
  table.className = "batch-table";

  const thead = document.createElement("thead");
  thead.innerHTML = `<tr>
    <th>ID</th><th>Length</th><th>Net charge</th>
    <th>Liability</th><th>Band</th><th>K_CK (nM)</th><th>Top variant</th>
  </tr>`;
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  results.forEach((r, idx) => {
    const variant = r.suggested_variants[0];
    const varSeq = variant ? variant.sequence : "";
    const varPreview = varSeq ? (varSeq.length > 12 ? varSeq.slice(0, 12) + "…" : varSeq) : "none";

    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td class="font-mono">${escHtml(r.id)}</td>
      <td>${r.length}</td>
      <td class="font-mono">${roundSig(r.net_charge, 3)}</td>
      <td class="font-mono">${roundSig(r.liability_score, 3)}</td>
      <td><span class="badge badge-${escHtml(r.liability_band)}">${escHtml(r.liability_band)}</span></td>
      <td class="font-mono">${roundSig(r.estimated_kck_nm, 4)}</td>
      <td class="font-mono" title="${escHtml(varSeq)}">${escHtml(varPreview)}</td>
    `;
    tr.addEventListener("click", () => batchShowDetail(idx));
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);

  wrap.innerHTML = "";
  wrap.appendChild(table);
}

function renderBatchBestStrip(results) {
  const strip = document.getElementById("scan-batch-best-strip");
  if (results.length === 0) { strip.style.display = "none"; return; }

  const best = results[0]; // index 0 = lowest liability = best (array is sorted)
  document.getElementById("scan-batch-best-id").textContent = best.id;
  document.getElementById("scan-batch-best-score").textContent = roundSig(best.liability_score, 3);
  const badge = document.getElementById("scan-batch-best-badge");
  badge.className = `badge badge-${best.liability_band}`;
  badge.textContent = best.liability_band;

  document.getElementById("scan-batch-chain-orig").onclick = () => {
    window.lockrChain.pipedKck = best.estimated_kck_nm;
    window.lockrChain.sourceLabel = `${best.id} (best of batch)`;
    showTab("calculator");
  };

  const varBtn = document.getElementById("scan-batch-chain-variant");
  const variant = best.suggested_variants[0];
  if (variant) {
    varBtn.style.display = "";
    varBtn.onclick = () => {
      window.lockrChain.pipedKck = variant.estimated_kck_nm;
      window.lockrChain.sourceLabel = `${best.id} variant (best of batch)`;
      showTab("calculator");
    };
  } else {
    varBtn.style.display = "none";
  }

  strip.style.display = "flex";
}

function batchExportCsv(results) {
  const header = ["ID", "Length", "Net Charge", "Liability Score", "Band", "K_CK (nM)", "Top Variant"];
  const rows = results.map(r => [
    r.id, r.length,
    roundSig(r.net_charge, 3), roundSig(r.liability_score, 3),
    r.liability_band, roundSig(r.estimated_kck_nm, 5),
    r.suggested_variants[0]?.sequence || "",
  ]);
  const csv = [header, ...rows]
    .map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(","))
    .join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = "lockr-batch-scan.csv";
  document.body.appendChild(a); a.click();
  document.body.removeChild(a); URL.revokeObjectURL(url);
}

function batchShowDetail(idx) {
  const r = batchState.results[idx];
  document.getElementById("scan-batch-table-section").style.display = "none";
  document.getElementById("scan-batch-best-strip").style.display = "none";
  const back = document.getElementById("scan-batch-detail-back");
  back.style.display = "flex";
  document.getElementById("scan-batch-detail-label").textContent = r.id;
  scanRenderResults(r);
}

function batchShowTable() {
  document.getElementById("scan-results").style.display = "none";
  document.getElementById("scan-batch-table-section").style.display = "block";
  document.getElementById("scan-batch-detail-back").style.display = "none";
  if (batchState.results.length > 0) {
    document.getElementById("scan-batch-best-strip").style.display = "flex";
  }
}

function setScanMode(mode) {
  scanMode = mode;
  document.getElementById("scan-single-input").style.display = mode === "single" ? "block" : "none";
  document.getElementById("scan-batch-input").style.display = mode === "batch" ? "block" : "none";
  document.getElementById("scan-adv-window-row").style.display = mode === "single" ? "" : "none";
  document.getElementById("scan-adv-preserve-row").style.display = mode === "single" ? "" : "none";
  scanEl("scan-empty-state").style.display = "block";
  scanEl("scan-results").style.display = "none";
  document.getElementById("scan-batch-results").style.display = "none";
  batchState.results = [];
  if (mode === "single") renderLiveAnnotation();
  else updateBatchCount();
}

async function batchSubmit() {
  const text = document.getElementById("scan-batch-text").value;
  const { records, errors } = parseBatchInput(text);
  const errEl = document.getElementById("scan-batch-error");
  errEl.textContent = "";
  if (records.length === 0) return;
  if (records.length > 500) {
    errEl.textContent = `${records.length} sequences exceeds the 500-sequence limit, remove ${records.length - 500} and retry.`;
    return;
  }

  const values = scanReadFields();
  const button = scanEl("scan-submit");
  const orig = button.innerHTML;
  button.disabled = true;
  button.innerHTML = `<span class="spinner"></span>Scanning ${records.length}…`;

  try {
    const results = await batchScan(records, values.ph, values.policy);
    batchState.results = results;

    scanEl("scan-empty-state").style.display = "none";
    scanEl("scan-results").style.display = "none";
    document.getElementById("scan-batch-results").style.display = "block";
    document.getElementById("scan-batch-table-section").style.display = "block";
    document.getElementById("scan-batch-detail-back").style.display = "none";

    document.getElementById("scan-batch-table-label").textContent =
      errors.length > 0
        ? `${results.length} sequences scanned · ${errors.length} line${errors.length !== 1 ? "s" : ""} skipped (see input)`
        : `${results.length} sequences scanned`;

    renderBatchTable(results);
    renderBatchBestStrip(results);
  } catch (err) {
    showToast(err.networkError ? err.message : `${err.code || "ERROR"}: ${err.message}`);
  } finally {
    button.disabled = false;
    button.innerHTML = orig;
  }
}


function buildRulerRow(length) {
  const row = document.createElement("div");
  row.className = "seq-ruler";
  const cells = new Array(length).fill("");
  for (let pos = 5; pos <= length; pos += 5) {
    const digits = String(pos).split("");
    digits.forEach((d, i) => {
      if (pos - 1 + i < length) cells[pos - 1 + i] = d;
    });
  }
  cells.forEach((ch) => {
    const span = document.createElement("span");
    span.textContent = ch || " ";
    row.appendChild(span);
  });
  return row;
}

function buildResidueRow(sequence, classify) {
  const row = document.createElement("div");
  row.className = "seq-residues";
  sequence.split("").forEach((ch, i) => {
    const span = document.createElement("span");
    span.className = `res ${classify(ch, i)}`;
    span.textContent = ch;
    row.appendChild(span);
  });
  return row;
}

function liveResidueClass(ch) {
  if (ACIDIC_AA.has(ch)) return "res-acidic";
  if (BASIC_AA.has(ch)) return "res-basic";
  return "res-muted";
}

let windowTouched = false;

function renderLiveAnnotation() {
  const textarea = document.getElementById("scan-sequence");
  const cursor = textarea.selectionStart;
  const cleaned = textarea.value.toUpperCase().replace(/\s/g, "").split("").filter((ch) => STANDARD_AA.has(ch)).join("");
  if (cleaned !== textarea.value) {
    const removed = textarea.value.length - cleaned.length;
    textarea.value = cleaned;
    textarea.setSelectionRange(Math.max(0, cursor - removed), Math.max(0, cursor - removed));
  }

  const liveRow = document.getElementById("scan-live-residues");
  liveRow.innerHTML = "";
  liveRow.appendChild(buildResidueRow(cleaned, liveResidueClass));

  const rulerSlot = document.getElementById("scan-ruler");
  rulerSlot.innerHTML = "";
  rulerSlot.appendChild(buildRulerRow(cleaned.length));

  const deCount = cleaned.split("").filter((ch) => ACIDIC_AA.has(ch)).length;
  const krhCount = cleaned.split("").filter((ch) => BASIC_AA.has(ch)).length;
  document.getElementById("scan-charcount").textContent = `${cleaned.length} residues`;
  document.getElementById("scan-decount").textContent = `${deCount} D/E`;
  document.getElementById("scan-netcharge-live").textContent = `net charge (rough) ~${krhCount - deCount}`;

  document.getElementById("scan-longseq-warning").style.display = cleaned.length > 200 ? "block" : "none";

  if (!windowTouched) {
    document.getElementById("scan-window-start").value = cleaned.length ? 1 : "";
    document.getElementById("scan-window-end").value = cleaned.length || "";
  }

  document.getElementById("scan-submit").disabled = cleaned.length === 0;

  scheduleBackgroundSequenceCheck(cleaned);
}

// --- Background sequence validation --------------------------------------
// Silent unless something's actually wrong, like a spell checker, not a form.

const BG_CHECK_MIN_LENGTH = 50;
const BG_CHECK_DEBOUNCE_MS = 800;
let bgCheckTimer = null;

function scanRenderBgWarning(warnings) {
  const wrap = document.getElementById("scan-bg-warning");
  const textEl = document.getElementById("scan-bg-warning-text");
  if (!warnings || warnings.length === 0) {
    wrap.style.display = "none";
    return;
  }
  wrap.style.display = "block";
  textEl.textContent = warnings.map((w) => `⚠ ${w}`).join("  ");
}

async function runBackgroundSequenceCheck(sequence) {
  if (sequence.length <= BG_CHECK_MIN_LENGTH) {
    scanRenderBgWarning([]);
    return;
  }
  try {
    const result = await apiPost("/check-sequence", { sequence });
    scanRenderBgWarning(result.warnings);
  } catch (err) {
    // background check failing quietly is fine, it's not the main flow
  }
}

function scheduleBackgroundSequenceCheck(sequence) {
  clearTimeout(bgCheckTimer);
  bgCheckTimer = setTimeout(() => runBackgroundSequenceCheck(sequence), BG_CHECK_DEBOUNCE_MS);
}

function scanEl(id) { return document.getElementById(id); }

function scanReadFields() {
  return {
    sequence: scanEl("scan-sequence").value,
    ph: parseFloat(scanEl("scan-ph").value),
    windowStart: parseInt(scanEl("scan-window-start").value, 10),
    windowEnd: parseInt(scanEl("scan-window-end").value, 10),
    preserveRaw: scanEl("scan-preserve").value,
    policy: document.querySelector("#scan-policy-segmented button.active").dataset.policy,
  };
}

// Mirrors what the backend itself rejects (window bounds, preserve_positions
// range) so the Scan request only ever goes out once it'll succeed.
function scanValidate(values) {
  const errors = {};
  const length = values.sequence.length;

  if (!(values.windowStart >= 1) || !(values.windowEnd >= values.windowStart) || values.windowEnd > length) {
    errors.window = `must be within 1-${length}, start <= end`;
  }

  const preservePositions = [];
  const raw = values.preserveRaw.trim();
  if (raw !== "") {
    for (const part of raw.split(",").map((s) => s.trim()).filter((s) => s !== "")) {
      const n = Number(part);
      if (!Number.isInteger(n) || n < 1 || n > length) {
        errors.preserve = `"${part}" is not a valid position (1-${length})`;
        break;
      }
      preservePositions.push(n);
    }
  }

  return { errors, preservePositions };
}

function scanShowFieldErrors(errors) {
  scanEl("scan-window-error").textContent = errors.window || "";
  scanEl("scan-preserve-error").textContent = errors.preserve || "";
}

function scanBuildRequest(values, preservePositions) {
  return {
    sequences: [{ id: "binder", sequence: values.sequence }],
    sensitive_window: { start: values.windowStart, end: values.windowEnd },
    ph: values.ph,
    substitution_policy: values.policy,
    preserve_positions: preservePositions,
  };
}

function annotatedResidueClass(residue, position, flaggedPositions) {
  if (flaggedPositions.has(position)) return ACIDIC_AA.has(residue) ? "res-acidic res-flagged" : "res-flagged";
  if (ACIDIC_AA.has(residue)) return "res-acidic";
  if (BASIC_AA.has(residue)) return "res-basic";
  return "res-muted";
}

function scanRenderAnnotatedSequence(result) {
  const flagged = new Set(result.acidic_residues.map((r) => r.position));
  const residueRow = scanEl("scan-annotated");
  residueRow.innerHTML = "";
  residueRow.appendChild(buildResidueRow(result.sequence, (ch, i) => annotatedResidueClass(ch, i + 1, flagged)));

  const ruler = scanEl("scan-annotated-ruler");
  ruler.innerHTML = "";
  ruler.appendChild(buildRulerRow(result.sequence.length));
}

function scanRenderContributionChart(result) {
  const container = scanEl("scan-contribution-chart");
  container.innerHTML = "";
  const anyFlagged = result.per_position.some((p) => p.contribution > 0);
  if (!anyFlagged) {
    container.innerHTML = '<div class="help-text">No charge liabilities flagged.</div>';
    return;
  }
  const heatmap = document.createElement("div");
  heatmap.className = "contrib-heatmap";
  result.per_position.forEach((p) => {
    const cell = document.createElement("div");
    cell.className = "contrib-cell";

    const tile = document.createElement("div");
    tile.className = `contrib-cell-tile ${p.contribution > 0 ? "contrib-cell-tile--flagged" : "contrib-cell-tile--neutral"}`;
    tile.textContent = p.residue;
    tile.title = p.contribution > 0 ? `${p.residue}${p.position}: ${p.contribution.toFixed(2)} kcal/mol penalty` : `${p.residue}${p.position}: no penalty`;

    const posLabel = document.createElement("div");
    posLabel.className = "contrib-cell-pos";
    posLabel.textContent = p.position % 5 === 0 ? p.position : "";

    cell.appendChild(tile);
    cell.appendChild(posLabel);
    heatmap.appendChild(cell);
  });
  container.appendChild(heatmap);
}

function scanRenderVariant(result) {
  const noneEl = scanEl("scan-variant-none");
  const boxesEl = scanEl("scan-variant-boxes");
  const variant = result.suggested_variants[0];

  if (!variant) {
    noneEl.style.display = "block";
    boxesEl.style.display = "none";
    scanEl("scan-variant-mutations").textContent = "";
    scanEl("scan-variant-score-before").textContent = roundSig(result.liability_score, 3);
    scanEl("scan-variant-score-after").textContent = roundSig(result.liability_score, 3);
    scanEl("scan-variant-kck-estimate").textContent = "N/A";
    scanEl("scan-variant-kck-estimate-nm").textContent = "N/A";
    scanEl("scan-variant-copy-row").style.display = "none";
    scanEl("scan-variant-kck-send-btn").style.display = "none";
    scanEl("scan-variant-helix").style.display = "none";
    return;
  }

  noneEl.style.display = "none";
  boxesEl.style.display = "flex";

  const mutated = new Map(variant.substitutions.map((s) => [s.position, s]));

  const originalEl = scanEl("scan-variant-original");
  originalEl.innerHTML = "";
  originalEl.appendChild(buildResidueRow(result.sequence, (ch, i) => mutated.has(i + 1) ? "res-old" : "res-muted"));

  const newEl = scanEl("scan-variant-new");
  newEl.innerHTML = "";
  newEl.appendChild(buildResidueRow(variant.sequence, (ch, i) => mutated.has(i + 1) ? "res-new" : "res-muted"));

  scanEl("scan-variant-mutations").textContent = variant.substitutions
    .map((s) => `${s.from_ ?? s.from}${s.position}${s.to}`)
    .join(", ");
  scanEl("scan-variant-score-before").textContent = roundSig(result.liability_score, 3);
  scanEl("scan-variant-score-after").textContent = roundSig(variant.liability_score, 3);
  const variantKckM = variant.estimated_kck_nm * 1e-9;
  scanEl("scan-variant-kck-estimate").textContent = `${variantKckM.toExponential(2)} M`;
  scanEl("scan-variant-kck-estimate-nm").textContent = `${roundSig(variant.estimated_kck_nm, 3)} nM`;

  const variantSendBtn = scanEl("scan-variant-kck-send-btn");
  variantSendBtn.style.display = "inline-block";
  variantSendBtn.onclick = () => {
    window.lockrChain.pipedKck = variant.estimated_kck_nm;
    window.lockrChain.sourceLabel = `${result.id || "scanned"} variant`;
    showTab("calculator");
  };

  scanEl("scan-variant-copy-row").style.display = "flex";
  scanEl("scan-variant-copy-sequence").textContent = variant.sequence;

  scanRenderVariantHelix(result, variant);
}

// Swapping D/E out changes the shape as well as the charge, so the variant gets its own
// structure verdict and the specific reasons it might be worse than the original.
function scanRenderVariantHelix(result, variant) {
  const row = scanEl("scan-variant-helix");
  if (!variant.helix || !result.helix) {
    row.style.display = "none";
    return;
  }
  row.style.display = "block";

  const before = Math.round(result.helix.helix_confidence * 100);
  const after = Math.round(variant.helix.helix_confidence * 100);
  const direction = after > before ? "better" : after < before ? "worse" : "unchanged";
  scanEl("scan-variant-helix-delta").textContent = `${before}% → ${after}% (${direction})`;

  const warnEl = scanEl("scan-variant-helix-warnings");
  warnEl.innerHTML = "";
  (variant.helix_warnings || []).forEach((w) => {
    const div = document.createElement("div");
    div.className = "help-text";
    div.style.color = "var(--warning-700)";
    div.style.marginTop = "4px";
    div.textContent = w;
    warnEl.appendChild(div);
  });
}

// --- Structure check -------------------------------------------------------
// Runs ahead of the charge story: a binder that isn't helical can't thread the latch, so
// the shape verdict is shown first and gates the graft step.

const HELIX_BAND_CLASS = {
  "likely helical": "badge-Low",
  uncertain: "badge-Moderate",
  "unlikely helical": "badge-High",
};

const HELIX_SEVERITY_COLOR = {
  blocking: "var(--danger-700, #b42318)",
  warning: "var(--warning-700)",
  info: "var(--text-muted, #667085)",
};

function scanRenderHelix(result) {
  const helix = result.helix;
  const card = scanEl("helix-card");
  if (!helix) {
    card.style.display = "none";
    return;
  }
  card.style.display = "block";

  const pct = Math.round(helix.helix_confidence * 100);
  scanEl("helix-confidence").textContent = `${pct}% helix confidence`;

  const badge = scanEl("helix-band-badge");
  badge.className = `badge ${HELIX_BAND_CLASS[helix.band] || ""}`;
  badge.textContent = helix.band.toUpperCase();

  const bridges = helix.salt_bridges.length;
  scanEl("helix-detail").textContent =
    `Mean helix propensity ${roundSig(helix.mean_propensity, 2)}, ` +
    `amphipathicity ${roundSig(helix.hydrophobic_moment, 2)}, ` +
    `${bridges} stabilizing salt ${bridges === 1 ? "bridge" : "bridges"} (i to i+3/i+4).`;

  const issuesEl = scanEl("helix-issues");
  issuesEl.innerHTML = "";
  helix.issues.forEach((issue) => {
    const div = document.createElement("div");
    div.className = "help-text";
    div.style.color = HELIX_SEVERITY_COLOR[issue.severity] || "";
    div.style.marginTop = "4px";
    div.textContent = issue.severity === "blocking" ? `Blocks grafting: ${issue.message}` : issue.message;
    issuesEl.appendChild(div);
  });

  const cyclicEl = scanEl("helix-cyclic");
  if (helix.cyclization.possibly_cyclic) {
    cyclicEl.style.display = "block";
    cyclicEl.textContent =
      `This might be a cyclic peptide (${helix.cyclization.signals.join("; ")}). ` +
      "Sequence alone can't tell, so check your construct. A cyclized binder can't be grafted " +
      "into the latch as a linear segment.";
  } else {
    cyclicEl.style.display = "none";
  }
}

// --- Grafting -------------------------------------------------------------
// Threads the current scan's (optimized) binder into the lucCage latch via
// PyRosetta and reports the best-scoring position. Lives entirely off the
// last /scan result, no separate sequence entry.

const graftState = { status: null, lastScanResult: null, jobId: null, kckNm: null, kckLabel: "" };

async function graftCheckStatus() {
  try {
    const res = await fetch(`${API_BASE}/graft/status`);
    graftState.status = await res.json();
  } catch (err) {
    graftState.status = { available: false, version: null, template_bundled: false };
  }
  graftUpdateAvailabilityUI();
}

function graftUpdateAvailabilityUI() {
  const available = graftState.status && graftState.status.available;
  scanEl("graft-unavailable-note").style.display = available ? "none" : "block";
  graftUpdateButtonState();
}

// The optimized variant, when the suggestion actually changed something.
function graftOptimizedVariant(result) {
  const variant = result.suggested_variants && result.suggested_variants[0];
  if (variant && variant.sequence !== result.sequence) return variant;
  return null;
}

// The sequence grafting actually uses. When an optimized variant exists the
// user picks between it and the sequence they entered, otherwise there is only
// one candidate.
function graftUsingSequence(result) {
  const variant = graftOptimizedVariant(result);
  const wantsOriginal = scanEl("graft-source-original").checked;
  if (variant && !wantsOriginal) {
    return {
      sequence: variant.sequence, kckNm: variant.estimated_kck_nm,
      label: "optimized variant", helix: variant.helix,
    };
  }
  return {
    sequence: result.sequence, kckNm: result.estimated_kck_nm,
    label: "scanned sequence", helix: result.helix,
  };
}

function graftRefreshUsingLabel() {
  if (!graftState.lastScanResult) return;
  const using = graftUsingSequence(graftState.lastScanResult);
  graftState.kckNm = using.kckNm;
  graftState.kckLabel = using.label;
  scanEl("graft-using-label").textContent = `Using: ${using.sequence}`;
}

function graftOnNewScanResult(result) {
  graftState.lastScanResult = result;

  // Offer the choice only when the two sequences actually differ.
  const variant = graftOptimizedVariant(result);
  scanEl("graft-source-choice").style.display = variant ? "block" : "none";
  if (variant) {
    scanEl("graft-source-original-seq").textContent = result.sequence;
    scanEl("graft-source-optimized-seq").textContent = variant.sequence;
    scanEl("graft-source-optimized").checked = true;
  }

  graftRefreshUsingLabel();
  graftReset();
  graftUpdateButtonState();
}

// A binder with a blocking shape problem can't thread the latch, so the graft button goes
// away rather than burning two minutes of PyRosetta on a result that can't be right.
function graftStructureBlock() {
  if (!graftState.lastScanResult) return null;
  const helix = graftUsingSequence(graftState.lastScanResult).helix;
  if (!helix || !helix.graft_blocked) return null;
  const reasons = helix.issues.filter((i) => i.severity === "blocking").map((i) => i.message);
  return `Can't graft this sequence. ${reasons.join(" ")}`;
}

function graftUpdateButtonState() {
  const available = graftState.status && graftState.status.available;
  const blocked = graftStructureBlock();

  const blockEl = scanEl("graft-structure-block");
  blockEl.style.display = blocked ? "block" : "none";
  blockEl.textContent = blocked || "";

  scanEl("graft-submit").disabled = !available || !graftState.lastScanResult || Boolean(blocked);
}

function graftReadRequest() {
  const using = graftUsingSequence(graftState.lastScanResult);
  const specific = scanEl("graft-specific-cb").checked;
  const position = specific ? parseInt(scanEl("graft-specific-position").value, 10) : null;
  return {
    sequence: using.sequence,
    scan_all: !specific,
    specific_position: specific && !isNaN(position) ? position : null,
  };
}

function graftReset() {
  scanEl("graft-progress").style.display = "none";
  scanEl("graft-error").style.display = "none";
  scanEl("graft-results").style.display = "none";
  graftState.jobId = null;
}

// A stuck PyRosetta run shouldn't hang the page forever, race the request
// against a 5-minute timer and report a clean timeout instead of spinning.
const GRAFT_TIMEOUT_MS = 5 * 60 * 1000;

function graftRequestWithTimeout(body) {
  return Promise.race([
    apiPost("/graft", body),
    new Promise((_, reject) => setTimeout(() => reject({ isTimeout: true }), GRAFT_TIMEOUT_MS)),
  ]);
}

const GRAFT_VERDICT_LABELS = { good: "GOOD", marginal: "MARGINAL", poor: "POOR", uncalibrated: "UNCALIBRATED" };

function graftRenderScoreChart(allScores, bestPosition) {
  const container = scanEl("graft-score-chart");
  container.innerHTML = "";
  if (allScores.length === 0) return;

  const scores = allScores.map((s) => s.score);
  const worst = Math.max(...scores);
  const best = Math.min(...scores);
  const span = worst - best || 1;

  const row = document.createElement("div");
  row.className = "graft-score-bars";
  allScores.forEach(({ position, score }) => {
    const col = document.createElement("div");
    col.className = "graft-score-bar-col";

    const bar = document.createElement("div");
    const heightPct = ((worst - score) / span) * 100;
    bar.className = `graft-score-bar${position === bestPosition ? " graft-score-bar--best" : ""}`;
    bar.style.height = `${Math.max(heightPct, 2)}%`;
    bar.title = `Position ${position}: ${roundSig(score, 5)} REU`;

    const label = document.createElement("div");
    label.className = "graft-score-bar-pos";
    label.textContent = position;

    col.appendChild(bar);
    col.appendChild(label);
    row.appendChild(col);
  });
  container.appendChild(row);
}

function graftRenderResults(data, usingSequence) {
  scanEl("graft-best-position").textContent = `Position ${data.best_position} in latch (325-359)`;
  scanEl("graft-best-score").textContent = `${roundSig(data.best_score, 6)} REU`;

  const badge = scanEl("graft-verdict-badge");
  badge.className = `badge badge-${data.verdict}`;
  badge.textContent = GRAFT_VERDICT_LABELS[data.verdict] || data.verdict;

  scanEl("graft-poor-note").style.display = data.verdict === "poor" ? "block" : "none";

  graftRenderScoreChart(data.all_scores, data.best_position);

  const binderStart = data.best_position;
  const binderEnd = data.best_position + usingSequence.length - 1;
  scanEl("graft-sequence-display").innerHTML = "";
  scanEl("graft-sequence-display").appendChild(
    buildResidueRow(data.grafted_sequence, (ch, i) => {
      const pos = i + 1;
      return pos >= binderStart && pos <= binderEnd ? "res-new" : "res-muted";
    })
  );

  graftState.jobId = data.job_id;

  const sendBtn = scanEl("graft-kck-send-btn");
  if ((data.verdict === "good" || data.verdict === "marginal") && graftState.kckNm !== null) {
    sendBtn.style.display = "inline-block";
    sendBtn.onclick = () => {
      window.lockrChain.pipedKck = graftState.kckNm;
      window.lockrChain.sourceLabel = `${graftState.kckLabel} (grafted, position ${data.best_position})`;
      showTab("calculator");
    };
  } else {
    sendBtn.style.display = "none";
  }

  scanEl("graft-results").style.display = "block";
}

async function graftSubmit() {
  if (!graftState.lastScanResult) return;
  graftReset();
  scanEl("graft-progress").style.display = "block";
  scanEl("graft-submit").disabled = true;

  const request = graftReadRequest();
  const usingSequence = request.sequence;

  try {
    const data = await graftRequestWithTimeout(request);
    graftRenderResults(data, usingSequence);
  } catch (err) {
    const errEl = scanEl("graft-error");
    errEl.style.display = "block";
    if (err && err.isTimeout) {
      errEl.textContent = "Grafting timed out. Try a shorter sequence or use a specific position.";
    } else if (err && err.code === "PYROSETTA_UNAVAILABLE") {
      errEl.textContent = "Grafting requires PyRosetta, see the Install Guide for setup.";
    } else if (err && err.code === "NO_VALID_POSITIONS") {
      errEl.textContent = "No valid graft positions found, binder may be too long or incompatible with the latch window geometry.";
    } else if (err && err.networkError) {
      showToast(err.message);
      errEl.style.display = "none";
    } else {
      errEl.textContent = (err && err.message) || "Grafting failed for an unknown reason.";
    }
  } finally {
    scanEl("graft-progress").style.display = "none";
    graftUpdateButtonState();
  }
}

function graftCopySequence() {
  const seq = graftState.lastScanResult
    ? document.querySelector("#graft-sequence-display").textContent
    : "";
  navigator.clipboard.writeText(seq).then(
    () => showToast("Grafted sequence copied."),
    () => showToast("Couldn't copy, select and copy the sequence manually.")
  );
}

function graftDownloadPdb() {
  if (!graftState.jobId) return;
  window.open(`${API_BASE}/graft/download/${graftState.jobId}`, "_blank");
}

function initGraft() {
  graftCheckStatus();

  scanEl("graft-specific-cb").addEventListener("change", (e) => {
    scanEl("graft-specific-position").style.display = e.target.checked ? "block" : "none";
  });

  // Switching source invalidates whatever was grafted from the other one.
  ["graft-source-original", "graft-source-optimized"].forEach((id) => {
    scanEl(id).addEventListener("change", () => {
      graftRefreshUsingLabel();
      graftReset();
    });
  });

  scanEl("graft-submit").addEventListener("click", graftSubmit);
  scanEl("graft-copy-btn").addEventListener("click", graftCopySequence);
  scanEl("graft-download-btn").addEventListener("click", graftDownloadPdb);
}

function scanRenderResults(result) {
  scanEl("scan-empty-state").style.display = "none";
  scanEl("scan-results").style.display = "block";

  scanEl("scan-band-label").textContent = `${result.liability_band[0].toUpperCase()}${result.liability_band.slice(1)} liability`;
  scanEl("scan-gauge-marker").style.left = `${result.liability_score}%`;
  scanEl("scan-liability-score").textContent = roundSig(result.liability_score, 3);
  scanEl("scan-result-ph").textContent = scanEl("scan-ph").value;
  scanEl("scan-net-charge").textContent = roundSig(result.net_charge, 3);

  const kckM = result.estimated_kck_nm * 1e-9;
  scanEl("scan-kck-estimate").textContent = `${kckM.toExponential(2)} M`;
  scanEl("scan-kck-estimate-nm").textContent = `${roundSig(result.estimated_kck_nm, 3)} nM`;

  const sendBtn = scanEl("scan-kck-send-btn");
  sendBtn.style.display = "inline-block";
  sendBtn.onclick = () => {
    window.lockrChain.pipedKck = result.estimated_kck_nm;
    window.lockrChain.sourceLabel = result.id || "scanned sequence";
    showTab("calculator");
  };

  const badge = scanEl("scan-kck-badge");
  badge.className = `badge badge-${result.predicted_kck_penalty.band}`;
  badge.textContent = result.predicted_kck_penalty.band;
  scanEl("scan-kck-note").textContent = result.predicted_kck_penalty.note;

  const warningsEl = scanEl("scan-warnings");
  warningsEl.innerHTML = "";
  (result.warnings || []).forEach((w) => {
    const div = document.createElement("div");
    div.className = "help-text";
    div.style.color = "var(--warning-700)";
    div.style.marginTop = "8px";
    div.textContent = w;
    warningsEl.appendChild(div);
  });

  scanRenderHelix(result);
  scanRenderAnnotatedSequence(result);
  scanRenderContributionChart(result);
  scanRenderVariant(result);
  graftOnNewScanResult(result);
}

async function scanSubmit() {
  const values = scanReadFields();
  const { errors, preservePositions } = scanValidate(values);
  scanShowFieldErrors(errors);
  if (Object.keys(errors).length > 0) return;

  const button = scanEl("scan-submit");
  const originalLabel = button.innerHTML;
  button.disabled = true;
  button.innerHTML = '<span class="spinner"></span>Scanning…';

  try {
    const response = await apiPost("/scan", scanBuildRequest(values, preservePositions));
    scanRenderResults(response.results[0]);
  } catch (err) {
    if (err.networkError) {
      showToast(err.message);
    } else if (err.field === "preserve_positions") {
      scanEl("scan-preserve-error").textContent = err.message;
    } else if (err.field === "sensitive_window") {
      scanEl("scan-window-error").textContent = err.message;
    } else {
      showToast(`${err.code || "ERROR"}: ${err.message}`);
    }
  } finally {
    button.disabled = false;
    button.innerHTML = originalLabel;
  }
}

function initScanner() {
  const textarea = document.getElementById("scan-sequence");
  textarea.addEventListener("input", renderLiveAnnotation);

  ["scan-window-start", "scan-window-end"].forEach((id) => {
    document.getElementById(id).addEventListener("input", () => { windowTouched = true; });
  });

  document.getElementById("scan-reset").addEventListener("click", () => {
    if (scanMode === "single") {
      textarea.value = "";
      windowTouched = false;
      clearTimeout(bgCheckTimer);
      scanRenderBgWarning([]);
      renderLiveAnnotation();
      scanEl("scan-empty-state").style.display = "block";
      scanEl("scan-results").style.display = "none";
      graftState.lastScanResult = null;
      graftReset();
      graftUpdateButtonState();
    } else {
      document.getElementById("scan-batch-text").value = "";
      batchState.results = [];
      document.getElementById("scan-batch-results").style.display = "none";
      scanEl("scan-empty-state").style.display = "block";
      updateBatchCount();
    }
    document.querySelectorAll("#scan-policy-segmented button").forEach((btn) => {
      btn.classList.toggle("active", btn.dataset.policy === "conservative");
    });
  });

  document.querySelectorAll("#scan-policy-segmented button").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll("#scan-policy-segmented button").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
    });
  });

  scanEl("scan-submit").addEventListener("click", () => {
    if (scanMode === "single") scanSubmit();
    else batchSubmit();
  });

  scanEl("scan-variant-copy-btn").addEventListener("click", async () => {
    const sequence = scanEl("scan-variant-copy-sequence").textContent;
    try {
      await navigator.clipboard.writeText(sequence);
      showToast("Variant sequence copied.");
    } catch (err) {
      showToast("Couldn't copy, select and copy the sequence manually.");
    }
  });

  // mode toggle
  document.querySelectorAll("#scan-mode-segmented button").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll("#scan-mode-segmented button").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      setScanMode(btn.dataset.mode);
    });
  });

  // batch textarea live parse
  document.getElementById("scan-batch-text").addEventListener("input", updateBatchCount);

  // file upload
  document.getElementById("scan-batch-upload-btn").addEventListener("click", () => {
    document.getElementById("scan-batch-file").click();
  });
  document.getElementById("scan-batch-file").addEventListener("change", (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      document.getElementById("scan-batch-text").value = ev.target.result;
      updateBatchCount();
    };
    reader.readAsText(file);
    e.target.value = "";
  });

  // drag-and-drop onto batch input area
  const dropzone = document.getElementById("scan-batch-dropzone");
  dropzone.addEventListener("dragover", (e) => { e.preventDefault(); dropzone.classList.add("drag-over"); });
  dropzone.addEventListener("dragleave", () => dropzone.classList.remove("drag-over"));
  dropzone.addEventListener("drop", (e) => {
    e.preventDefault();
    dropzone.classList.remove("drag-over");
    const file = e.dataTransfer.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      document.getElementById("scan-batch-text").value = ev.target.result;
      updateBatchCount();
    };
    reader.readAsText(file);
  });

  // batch results controls
  document.getElementById("scan-batch-back-btn").addEventListener("click", batchShowTable);
  document.getElementById("scan-batch-export-btn").addEventListener("click", () => {
    if (batchState.results.length > 0) batchExportCsv(batchState.results);
  });

  renderLiveAnnotation();
  initGraft();
}

document.addEventListener("DOMContentLoaded", initScanner);
