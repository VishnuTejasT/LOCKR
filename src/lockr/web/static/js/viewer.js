// 3D structure views (3Dmol.js). Two stages share one fetched template:
//   - the header hero: the lucCage template, slowly turning, purely decorative
//   - the Scanner's Structure card: template, then the grafted PDB with the binder highlighted
// Colours match the legend in index.html: cage sky blue, latch amber, grafted binder coral.

const lockrViewer = (() => {
  const CAGE = "#8e95a6";
  const LATCH = "#f5a02e";
  const BINDER = "#ef6f5e";
  const LATCH_RANGE = "325-359";
  const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  const el = (id) => document.getElementById(id);
  let templateText = null;
  let main = null; // the Scanner card's stage
  let hero = null;

  async function fetchText(url) {
    const response = await fetch(url);
    if (!response.ok) throw new Error(String(response.status));
    return response.text();
  }

  async function getTemplate() {
    if (templateText === null) templateText = await fetchText(`${API_BASE}/graft/template`);
    return templateText;
  }

  function makeStage(stageId, options) {
    const node = el(stageId);
    if (!node || !window.$3Dmol) return null;
    return { v: $3Dmol.createViewer(node, { backgroundAlpha: 0, antialias: true }), node, shown: null, ...options };
  }

  function draw(stage, pdb, key, binderRange) {
    const { v } = stage;
    v.removeAllModels();
    v.addModel(pdb, "pdb");
    if (stage.hero) {
      // the same thick grey tube and orange latch as HybriDock-Pep's hero
      v.setStyle({}, { cartoon: { style: "trace", color: CAGE, thickness: 1.1 } });
      v.setStyle({ resi: LATCH_RANGE }, { cartoon: { style: "trace", color: LATCH, thickness: 1.4 } });
    } else {
      v.setStyle({}, { cartoon: { style: "trace", color: CAGE, thickness: 0.8 } });
      v.setStyle({ resi: LATCH_RANGE }, { cartoon: { style: "trace", color: LATCH, thickness: 1.1 } });
      if (binderRange) {
        v.setStyle({ resi: binderRange }, { cartoon: { style: "trace", color: BINDER, thickness: 1.3 } });
      }
    }
    v.zoomTo();
    // hero turns about the cage's long axis, so its silhouette stays wide instead of going end-on
    if (stage.hero) v.zoom(0.92);
    v.render();
    if (!reduceMotion) v.spin(stage.hero ? "x" : "y", stage.hero ? 0.4 : 0.6);
    stage.shown = key;
  }

  function message(text) {
    const box = el("viewer-message");
    if (!box) return;
    box.textContent = text || "";
    box.style.display = text ? "flex" : "none";
  }

  async function show(url, key, binderRange) {
    if (!main) main = makeStage("viewer-stage", { hero: false });
    if (!main) { message("The 3D viewer could not load. Check your internet connection and reload."); return; }
    if (main.shown === key) return;
    let pdb;
    try {
      pdb = key === "template" ? await getTemplate() : await fetchText(url);
    } catch (err) {
      message("Could not load the structure.");
      return;
    }
    draw(main, pdb, key, binderRange);
    message("");
    el("viewer-legend-binder").hidden = !binderRange;
    el("viewer-sub").textContent = binderRange
      ? `Grafted lucCage, binder at residues ${binderRange.replace("-", " to ")}.`
      : "lucCage template. Graft a binder to see it in place.";
    el("viewer-spin").setAttribute("aria-pressed", String(!reduceMotion));
  }

  const showTemplate = () => show(null, "template", null);
  const showGraft = (url, start, end) => show(url, url, `${start}-${end}`);

  async function initHero() {
    if (!el("hero-stage")) return;
    hero = makeStage("hero-stage", { hero: true });
    if (!hero) return;
    try { draw(hero, await getTemplate(), "template", null); } catch (err) { /* the header reads fine without it */ }
  }

  function init() {
    initHero();
    if (!el("viewer-stage")) return;
    el("viewer-spin").addEventListener("click", () => {
      if (!main) return;
      const on = el("viewer-spin").getAttribute("aria-pressed") !== "true";
      main.v.spin(on ? "y" : false, 0.6);
      el("viewer-spin").setAttribute("aria-pressed", String(on));
    });
    el("viewer-reset").addEventListener("click", () => {
      if (main) { main.v.zoomTo(); main.v.render(); }
    });
    if (window.ResizeObserver) {
      // the stage changes size when scan results appear beside it, so refit rather than stay tiny
      new ResizeObserver(() => {
        if (main && main.shown) { main.v.resize(); main.v.zoomTo(); main.v.render(); }
      }).observe(el("viewer-stage"));
    }
    document.addEventListener("tabchange", (e) => {
      if (e.detail.tab !== "scanner") return;
      if (!main || main.shown === null) showTemplate();
      else { main.v.resize(); main.v.render(); }
    });
    if (document.querySelector('.tab-panel.active[data-tab="scanner"]')) showTemplate();
  }

  document.addEventListener("DOMContentLoaded", init);
  return { showTemplate, showGraft };
})();
