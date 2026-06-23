//This file works to share information between the scanner and calculator tabs of the frontend.

const API_BASE = ""; // same origin... lockr serve hosts API and frontend together

// RT at 37°C in kcal/mol -- used for dG<->Kd conversion in calculator
const RT_KCAL_MOL = 0.592;

//Scanner writes info here, while Clauclator reads from it.
window.lockrChain = { pipedKck: null, sourceLabel: null };

function showTab(name) {
  document.querySelectorAll(".tab-panel").forEach((el) => el.classList.toggle("active", el.dataset.tab === name));
  document.querySelectorAll("nav.tabs button").forEach((el) => el.classList.toggle("active", el.dataset.tab === name));
  location.hash = `#${name}`;
  document.dispatchEvent(new CustomEvent("tabchange", { detail: { tab: name } }));
}

function initTabs() {
  document.querySelectorAll("nav.tabs button").forEach((btn) => {
    btn.addEventListener("click", () => showTab(btn.dataset.tab));
  });
  const initial = location.hash.replace("#", "") || "scanner";
  showTab(["scanner", "calculator", "assembly"].includes(initial) ? initial : "scanner");
}

let toastTimer = null;
function showToast(message) {
  let el = document.getElementById("toast");
  if (!el) {
    el = document.createElement("div");
    el.id = "toast";
    el.className = "toast";
    document.body.appendChild(el);
  }
  el.textContent = message;
  el.style.display = "block";
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { el.style.display = "none"; }, 4000);
}


async function apiPost(path, body) {
  let response;
  try {
    response = await fetch(API_BASE + path, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    });
  } catch (err) {
    throw { networkError: true, message: "Engine isn't working... Is lockr serve command running?" };
  }
  const data = await response.json();
  if (!response.ok) {
    throw { networkError: false, status: response.status, ...data.error };
  }
  return data;
}

// nM in, human-scaled string out -- e.g. 0.000042 nM -> "42 fM", never raw sci notation.
function formatConcNm(nm) {
  if (nm === null || nm === undefined || Number.isNaN(nm)) return "—";
  const molar = nm * 1e-9;
  const units = [["fM", 1e-15], ["pM", 1e-12], ["nM", 1e-9], ["µM", 1e-6], ["mM", 1e-3], ["M", 1]];
  for (let i = 0; i < units.length; i++) {
    const [label, scale] = units[i];
    const next = units[i + 1];
    if (!next || molar < next[1]) {
      return `${roundSig(molar / scale, 3)} ${label}`;
    }
  }
  return `${molar} M`;
}

function roundSig(x, sig) {
  if (x === 0) return "0";
  const mag = Math.ceil(Math.log10(Math.abs(x)));
  const precision = sig - mag;
  const factor = Math.pow(10, precision);
  return (Math.round(x * factor) / factor).toString();
}

function formatFoldChange(x) {
  return `${roundSig(x, 3)}×`;
}

document.addEventListener("DOMContentLoaded", initTabs);
