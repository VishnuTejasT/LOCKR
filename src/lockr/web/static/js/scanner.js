

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
    errEl.textContent = `${records.length} sequences exceeds the 500-sequence limit — remove ${records.length - 500} and retry.`;
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
    scanEl("scan-variant-kck-estimate").textContent = "—";
    scanEl("scan-variant-kck-estimate-nm").textContent = "—";
    scanEl("scan-variant-copy-row").style.display = "none";
    scanEl("scan-variant-kck-send-btn").style.display = "none";
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
}

function scanRenderResults(result) {
  scanEl("scan-empty-state").style.display = "none";
  scanEl("scan-results").style.display = "block";

  scanEl("scan-band-label").textContent = `${result.liability_band[0].toUpperCase()}${result.liability_band.slice(1)} liability`;
  scanEl("scan-gauge-marker").style.left = `${result.liability_score}%`;
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

  scanRenderAnnotatedSequence(result);
  scanRenderContributionChart(result);
  scanRenderVariant(result);
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
      renderLiveAnnotation();
      scanEl("scan-empty-state").style.display = "block";
      scanEl("scan-results").style.display = "none";
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
      showToast("Couldn't copy — select and copy the sequence manually.");
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
}

document.addEventListener("DOMContentLoaded", initScanner);
