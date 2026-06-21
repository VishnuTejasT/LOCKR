// Scanner tab. Step 2: live client-side ruler/annotation on the sequence
// textarea -- no API call involved, only `/scan` submission (not yet wired)
// touches the network. Kept isolated from that path on purpose.

const STANDARD_AA = new Set("ACDEFGHIKLMNPQRSTVWY".split(""));
const ACIDIC_AA = new Set(["D", "E"]);
const BASIC_AA = new Set(["K", "R", "H"]);

// Shared by the live preview here and the post-scan annotated-sequence card
// (step 3) -- both need "every 5th position numbered" under a mono sequence.
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
  const flagged = result.per_position.filter((p) => p.contribution > 0);
  if (flagged.length === 0) {
    container.innerHTML = '<div class="help-text">No charge liabilities flagged.</div>';
    return;
  }
  const maxContribution = Math.max(...flagged.map((p) => p.contribution));
  flagged.forEach((p) => {
    const row = document.createElement("div");
    row.className = "contribution-row";
    const label = document.createElement("div");
    label.className = "contribution-label";
    label.textContent = `${p.residue}${p.position}`;
    const track = document.createElement("div");
    track.className = "contribution-bar-track";
    const fill = document.createElement("div");
    fill.className = "contribution-bar-fill";
    fill.style.width = `${(p.contribution / maxContribution) * 100}%`;
    track.appendChild(fill);
    row.appendChild(label);
    row.appendChild(track);
    container.appendChild(row);
  });
}

function scanRenderVariant(result) {
  const diffEl = scanEl("scan-variant-diff");
  diffEl.innerHTML = "";
  const variant = result.suggested_variants[0];

  if (!variant) {
    diffEl.textContent = "No liabilities found -- no variant needed.";
    scanEl("scan-variant-mutations").textContent = "";
    scanEl("scan-variant-score-before").textContent = roundSig(result.liability_score, 3);
    scanEl("scan-variant-score-after").textContent = roundSig(result.liability_score, 3);
    return;
  }

  const mutated = new Map(variant.substitutions.map((s) => [s.position, s]));

  const originalRow = buildResidueRow(result.sequence, (ch, i) => mutated.has(i + 1) ? "res-old" : "res-muted");
  const variantRow = buildResidueRow(variant.sequence, (ch, i) => mutated.has(i + 1) ? "res-new" : "res-muted");
  diffEl.appendChild(originalRow);
  diffEl.appendChild(variantRow);

  scanEl("scan-variant-mutations").textContent = variant.substitutions
    .map((s) => `${s.from_ ?? s.from}${s.position}${s.to}`)
    .join(", ");
  scanEl("scan-variant-score-before").textContent = roundSig(result.liability_score, 3);
  scanEl("scan-variant-score-after").textContent = roundSig(variant.liability_score, 3);
}

function scanRenderResults(result) {
  scanEl("scan-empty-state").style.display = "none";
  scanEl("scan-results").style.display = "block";

  scanEl("scan-band-label").textContent = `${result.liability_band[0].toUpperCase()}${result.liability_band.slice(1)} liability`;
  scanEl("scan-gauge-marker").style.left = `${result.liability_score}%`;
  scanEl("scan-result-ph").textContent = scanEl("scan-ph").value;
  scanEl("scan-net-charge").textContent = roundSig(result.net_charge, 3);

  const badge = scanEl("scan-kck-badge");
  badge.className = `badge badge-${result.predicted_kck_penalty.band}`;
  badge.textContent = result.predicted_kck_penalty.band;
  scanEl("scan-kck-note").textContent = result.predicted_kck_penalty.note;

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
    textarea.value = "";
    windowTouched = false;
    renderLiveAnnotation();
    document.querySelectorAll("#scan-policy-segmented button").forEach((btn) => {
      btn.classList.toggle("active", btn.dataset.policy === "conservative");
    });
    scanEl("scan-empty-state").style.display = "block";
    scanEl("scan-results").style.display = "none";
  });

  document.querySelectorAll("#scan-policy-segmented button").forEach((btn) => {
    btn.addEventListener("click", () => {
      document.querySelectorAll("#scan-policy-segmented button").forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
    });
  });

  scanEl("scan-submit").addEventListener("click", scanSubmit);

  renderLiveAnnotation();
}

document.addEventListener("DOMContentLoaded", initScanner);
