"""Interactive 3D viewer HTML built on 3Dmol.js (spec stage G).

The viewer is a self-contained page placed in an ``<iframe srcdoc=...>``, so
it runs the same way in classic Notebook, AppMode, JupyterLab and Voila
regardless of how each front end treats scripts in widget output. Rotation,
zoom and pan are 3Dmol's mouse controls; style, ligand visibility, site
highlighting and view reset are toolbar controls inside the page, so they
keep the camera. Hovering an atom shows its index, element, role and
molecule ID. Views: Top, Side, Angled, Reset; orthographic (default for
periodic models) or perspective projection; a Slice slider that clips the
front of the model (optionally to a thin slab) to look inside; background
white, black, theme or transparent; PNG saves a 3x off-screen render.

Clicking an atom does what the "Click" menu says: inspect (label the atom),
measure (distance between two atoms, angle for three), or, after Generate,
remove the clicked ligand or add a ligand on the clicked surface atom. Edits
are sent to the Python app through a hidden text widget (``bridge``); the app
rebuilds the model and records the edit in README.txt. PNG saves the current
picture.

3Dmol.js is loaded from cdnjs, so the *browser* needs internet access.
"""

from __future__ import annotations

import html
import itertools
import json

import numpy as np

from .theme import BUDDY, BUDDY2, DARK, DISPLAY_STACK, FONT_STACK, LIGHT, font_face_css

THREEDMOL_URL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"
LARGE_MODEL_ATOMS = 6000


DRAW_BOND_MAX_A = 2.5  # longer drawn bonds cross the periodic boundary: not drawn


def _drawable(bonds, positions) -> list:
    pos = np.asarray(positions, dtype=float)
    out = []
    for b in bonds:
        i, j = int(b[0]), int(b[1])
        if np.linalg.norm(pos[i] - pos[j]) <= DRAW_BOND_MAX_A:
            out.append([i, j])
    return out


def _model_payload(elements, positions, roles, molecule_ids, bonds, cell=None,
                   site_atoms=(), title="", show_cell=True, side_on=False,
                   substrate_bonds=False) -> dict:
    return {
        "editable": False,
        "sel": [],
        "view": None,
        "pbc": [False, False, False],
        "sideOn": side_on if isinstance(side_on, str) else bool(side_on),
        "showCell": bool(show_cell),
        "title": title,
        "elem": list(elements),
        "xyz": np.round(np.asarray(positions, dtype=float), 3).tolist(),
        "role": list(roles),
        "mol": [int(m) for m in molecule_ids],
        "bonds": _drawable(bonds, positions),
        "cell": None if cell is None else np.round(np.asarray(cell, dtype=float), 4).tolist(),
        "sites": sorted(int(s) for s in site_atoms),
        "large": len(elements) > LARGE_MODEL_ATOMS,
        "subBonds": bool(substrate_bonds),
    }


SIDE_ON_SHAPES = ("rod", "nanotube", "channel")


def _side(shape, channel_type=None):
    """Start view: rods and tubes side-on (axis horizontal); slit channels look along y with
    z up; pore membranes are seen at an angle so the pores show."""
    if shape == "channel":
        return "angled" if channel_type == "pore" else "up"
    return shape in SIDE_ON_SHAPES


def assembly_payload(assembly, title="") -> dict:
    sites = [s["surface_atom_index"] for s in assembly.site_to_ligand]
    bonds = [b for b in assembly.bond_list]
    shape = (assembly.provenance.get("geometry") or {}).get("shape")
    covalent = assembly.residue_names is not None
    payload = _model_payload(assembly.elements, assembly.positions_A, assembly.chemical_role,
                             assembly.molecule_id, bonds, assembly.cell_A, sites, title,
                             show_cell=True,
                             side_on=_side(shape, (assembly.provenance.get("geometry") or {})
                                           .get("channel_type")),
                             substrate_bonds=covalent)
    payload["pbc"] = [bool(x) for x in assembly.pbc]
    payload["ortho"] = bool(any(assembly.pbc))  # periodic edges line up in orthographic view
    return payload


def substrate_payload(substrate, title="") -> dict:
    atoms = substrate.atoms
    n = len(atoms)
    return _model_payload(atoms.get_chemical_symbols(), atoms.get_positions(), ["substrate"] * n,
                          [0] * n, substrate.bonds, np.array(atoms.cell), [], title,
                          show_cell=any(substrate.pbc),
                          side_on=_side(substrate.shape, (substrate.channel_info or {}).get("type")),
                          substrate_bonds=bool(substrate.bonds))


def unit_cell_payload(crystal, title="") -> dict:
    """Atoms of one unit cell; atoms on faces/edges/corners are repeated for display only."""
    atoms = crystal.atoms
    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=True)
    elements, points = [], []
    for symbol, f in zip(atoms.get_chemical_symbols(), frac):
        options = [[0.0, 1.0] if abs(c) < 1e-4 or abs(c - 1) < 1e-4 else [c] for c in f]
        for combo in itertools.product(*options):
            elements.append(symbol)
            points.append(np.array(combo) @ cell)
    n = len(elements)
    return _model_payload(elements, points, ["substrate"] * n, [0] * n, [], cell, [], title,
                          side_on="angled")


_PAGE = """<!doctype html>
<html data-theme="__THEME__"><head><meta charset="utf-8">
<script src="__URL__"></script>
<style>
 __FONTFACE__
 :root { __LIGHT__ }
 :root[data-theme=dark] { __DARK__ }
 :root[data-theme=buddy] { __BUDDY__ }
 :root[data-theme=buddy2] { __BUDDY2__ }
 html, body { margin: 0; height: 100%; background: var(--panel); color: var(--ink); }
 body { font: 12px/1.4 __FONT__; display: flex; flex-direction: column; }
 .bar { display: flex; flex-wrap: wrap; gap: 6px 12px; align-items: center; padding: 6px 10px;
        border-bottom: 1px solid var(--rule); background: var(--panel); }
 #title { font-family: __DISPLAY__; font-size: 14px; font-weight: 600; margin-right: 4px;
          font-variant-numeric: tabular-nums; }
 label { display: inline-flex; align-items: center; gap: 5px; color: var(--muted); }
 select, button { font: inherit; color: var(--ink); background: var(--input); border: 1px solid var(--rule);
                  border-radius: 6px; height: 26px; padding: 0 8px; }
 button { cursor: pointer; }
 button:hover, select:hover { border-color: var(--muted); }
 button.on { background: var(--au-soft); border-color: var(--au); font-weight: 600; }
 :focus-visible { outline: 2px solid var(--au); outline-offset: 1px; }
 input[type=checkbox], input[type=range] { accent-color: var(--au); margin: 0; }
 input[type=range] { width: 110px; }
 .grp { display: inline-flex; gap: 0; }
 .grp button { border-radius: 0; margin-left: -1px; }
 .grp button:first-child { border-radius: 6px 0 0 6px; margin-left: 0; }
 .grp button:last-child { border-radius: 0 6px 6px 0; }
 #v { position: relative; width: 100%; flex: 1 1 auto; min-height: 240px;
      background-image: linear-gradient(45deg, #e5e7eb 25%, transparent 25%), linear-gradient(-45deg, #e5e7eb 25%, transparent 25%),
        linear-gradient(45deg, transparent 75%, #e5e7eb 75%), linear-gradient(-45deg, transparent 75%, #e5e7eb 75%);
      background-size: 16px 16px; background-position: 0 0, 0 8px, 8px -8px, -8px 0; }
 #v.solid { background-image: none; }
 .bar, #hover { flex: 0 0 auto; }
 #hover { color: var(--muted); min-height: 16px; padding: 6px 10px; border-top: 1px solid var(--rule);
          background: var(--panel); font-variant-numeric: tabular-nums; }
 #hover b { color: var(--ink); }
 #err { color: var(--err); padding: 16px; }
</style></head>
<body>
<div class="bar">
 <b id="title"></b>
 <label>Style <select id="style">
  <option value="bs">Ball and stick</option>
  <option value="sf">Space filling</option>
  <option value="wf">Wireframe</option></select></label>
 <label><input type="checkbox" id="lig" checked> Ligands</label>
 <label title="Highlight the surface atoms that carry a ligand"><input type="checkbox" id="sites"> Sites</label>
 <label><input type="checkbox" id="cell" checked> Box</label>
 <label>Background <select id="bg">
  <option value="white">White</option>
  <option value="black">Black</option>
  <option value="theme">Theme</option>
  <option value="transparent">Transparent</option></select></label>
</div>
<div class="bar">
 <span class="grp"><button id="vtop" title="Look down z">Top</button><button id="vside" title="Look along y, z up">Side</button><button id="vang" title="Angled 3D view">Angled</button><button id="reset" title="Fit the model">Reset</button></span>
 <span class="grp"><button id="ortho" title="Orthographic: parallel lines stay parallel; periodic edges line up">Ortho</button><button id="persp" title="Perspective: depth cue">Persp</button></span>
 <label title="Cut away the front of the model to see inside (along the current view)">Slice <input type="range" id="slice" min="0" max="95" value="0"></label>
 <label title="Show only a thin slab at the slice depth"><input type="checkbox" id="thin"> Thin</label>
</div>
<div class="bar">
 <label>Click <select id="mode" title="What a click on an atom does">
  <option value="inspect">Inspect atom</option>
  <option value="measure">Measure distance / angle</option>
  <option value="remove">Remove ligand</option>
  <option value="add">Add ligand on surface atom</option></select></label>
 <button id="clear" title="Remove labels and measurements">Clear labels</button>
 <button id="png" title="Save the current picture as a high-resolution PNG (3x)">PNG (hi-res)</button>
</div>
<div id="v"></div>
<div id="hover">Drag to rotate &middot; scroll or right-drag to zoom &middot; Ctrl-drag or middle-drag to pan &middot; hover or click an atom for details.</div>
<script>
const D = __DATA__;
const BRIDGE = "__BRIDGE__";
const hover = document.getElementById("hover");
const say = (html) => { hover.innerHTML = html; };
const css = (name, fallback) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
if (typeof $3Dmol === "undefined") {
  document.getElementById("v").innerHTML =
    '<div id="err">3Dmol.js could not be loaded from cdnjs. The 3D viewer needs internet access in the browser.</div>';
} else {
  document.getElementById("title").textContent = D.title;
  const accent = css("--au", "#2E7D55");
  const isSub = {predicate: a => a.properties.role === "substrate"};
  const isLig = {predicate: a => a.properties.role !== "substrate"};
  const siteSet = new Set(D.sites);
  const isSite = {predicate: a => siteSet.has(a.index)};
  const selSet = new Set(D.sel || []);
  const isSel = {predicate: a => selSet.has(a.index)};
  const nb = D.elem.map(() => []);
  for (const [i, j] of D.bonds) { nb[i].push(j); nb[j].push(i); }
  const bonded = new Set(D.bonds.flat());
  const state = {ortho: D.ortho, slice: 0, thin: false, bg: D.background || "white"};
  function bgColor() {
    return {white: "#ffffff", black: "#000000", theme: css("--stage", "#ffffff"),
            transparent: "#ffffff"}[state.bg];
  }
  function isDarkBg() { return state.bg === "black" || (state.bg === "theme" && document.documentElement.dataset.theme === "dark"); }
  function makeViewer(el) {
    const v = $3Dmol.createViewer(el, {backgroundColor: bgColor(), backgroundAlpha: state.bg === "transparent" ? 0 : 1, antialias: true});
    const m = v.addModel();
    m.addAtoms(D.elem.map((e, i) => ({
      elem: e, x: D.xyz[i][0], y: D.xyz[i][1], z: D.xyz[i][2], serial: i,
      bonds: nb[i], bondOrder: nb[i].map(() => 1),
      properties: {role: D.role[i], mol: D.mol[i]}})));
    v.bbModel = m;
    return v;
  }
  function styleViewer(v) {
    const style = document.getElementById("style").value;
    const showLig = document.getElementById("lig").checked;
    v.setStyle({}, {});
    // Substrate atoms bonded to a ligand (Au-S, chemisorption) also get a stick, so the bond shows.
    const subStick = {predicate: a => a.properties.role === "substrate" && bonded.has(a.index)};
    if (style === "bs" && D.subBonds) {
      // covalent frameworks (graphene, nanotubes, fullerenes): sticks, not packed spheres
      v.setStyle(isSub, D.large ? {stick: {radius: 0.1, colorscheme: "Jmol"}}
                                : {stick: {radius: 0.12, colorscheme: "Jmol"},
                                   sphere: {scale: 0.18, colorscheme: "Jmol"}});
      if (showLig) v.setStyle(isLig, {stick: {radius: 0.15, colorscheme: "Jmol"},
                                     sphere: {scale: 0.22, colorscheme: "Jmol"}});
    } else if (style === "bs") {
      v.setStyle(isSub, {sphere: {scale: D.large ? 0.3 : 0.4, colorscheme: "Jmol"}});
      if (showLig) {
        v.setStyle(subStick, {sphere: {scale: D.large ? 0.3 : 0.4, colorscheme: "Jmol"},
                              stick: {radius: 0.15, colorscheme: "Jmol"}});
        v.setStyle(isLig, {stick: {radius: 0.15, colorscheme: "Jmol"},
                           sphere: {scale: 0.22, colorscheme: "Jmol"}});
      }
    } else if (style === "sf") {
      v.setStyle(isSub, {sphere: {scale: 1.0, colorscheme: "Jmol"}});
      if (showLig) v.setStyle(isLig, {sphere: {scale: 1.0, colorscheme: "Jmol"}});
    } else {
      v.setStyle(isSub, D.subBonds ? {line: {linewidth: 1, colorscheme: "Jmol"}}
                                   : {sphere: {scale: 0.15, colorscheme: "Jmol"}});
      if (showLig) v.setStyle(isLig, {line: {linewidth: 2, colorscheme: "Jmol"}});
    }
    // setStyle, not addStyle: a merged style would keep the Jmol colorscheme over the color.
    if (document.getElementById("sites").checked)
      v.setStyle(isSite, {sphere: {scale: 0.55, color: "magenta"}});
    if (selSet.size) v.setStyle(isSel, {sphere: {scale: 0.6, color: "#00B8D9"}});
    drawCell(v, document.getElementById("cell").checked);
    v.setProjection(state.ortho ? "orthographic" : "perspective");
  }
  function drawCell(v, show) {
    for (const s of (v.bbCell || [])) v.removeShape(s);
    v.bbCell = [];
    if (!show || !D.cell) return;
    const [a, b, c] = D.cell;
    const add = (p, q) => [p[0] + q[0], p[1] + q[1], p[2] + q[2]];
    const o = [0, 0, 0];
    const corners = [o, a, b, c, add(a, b), add(a, c), add(b, c), add(add(a, b), c)];
    const edges = [[0,1],[0,2],[0,3],[1,4],[1,5],[2,4],[2,6],[3,5],[3,6],[4,7],[5,7],[6,7]];
    const color = isDarkBg() ? "#D7DEE8" : "#1F2733";
    const span = Math.max(...D.cell.map(r => Math.hypot(r[0], r[1], r[2])));
    const radius = Math.min(0.35, Math.max(0.12, span / 260));
    for (const [p, q] of edges) {
      const s = corners[p], t = corners[q];
      v.bbCell.push(v.addCylinder({start: {x: s[0], y: s[1], z: s[2]},
        end: {x: t[0], y: t[1], z: t[2]}, radius: radius, color: color, fromCap: 1, toCap: 1}));
    }
  }
  const viewer = makeViewer(document.getElementById("v"));
  const model = viewer.bbModel;
  window.nbViewer = viewer;  // for debugging from the browser console
  let slab0 = null;
  function applySlice(v) {
    if (!slab0) return;
    if (!state.slice) { v.setSlab(slab0.near, slab0.far); return; }
    const depth = slab0.near + (slab0.far - slab0.near) * state.slice / 100;
    v.setSlab(depth, state.thin ? depth + 4.0 : slab0.far);
  }
  function applyBackground() {
    viewer.setBackgroundColor(bgColor(), state.bg === "transparent" ? 0 : 1);
    document.getElementById("v").classList.toggle("solid", state.bg !== "transparent");
  }
  function apply() {
    styleViewer(viewer);
    applySlice(viewer);
    viewer.render();
  }
  const describe = a => `<b>Atom ${a.index}</b>: ${a.elem}, ${a.properties.role}` +
    (a.properties.mol ? `, ligand ${a.properties.mol}` : "") +
    ` at (${a.x.toFixed(2)}, ${a.y.toFixed(2)}, ${a.z.toFixed(2)}) Å`;
  viewer.setHoverable({}, true, a => say(describe(a)), a => {});
  // --- clicks: inspect, measure, remove / add ligand ---------------------------------
  let labels = [], marks = [], picks = [];
  const label = (text, a) => labels.push(viewer.addLabel(text, {position: {x: a.x, y: a.y, z: a.z},
    fontSize: 12, fontColor: "black", backgroundColor: "white", backgroundOpacity: 0.85,
    borderThickness: 1, borderColor: accent, inFront: true}));
  function clearMarks() {
    for (const l of labels) viewer.removeLabel(l);
    for (const m of marks) viewer.removeShape(m);
    labels = []; marks = []; picks = [];
    viewer.render();
  }
  function mic(d) {  // minimum image along periodic, orthogonal box axes
    if (!D.cell) return d;
    for (let k = 0; k < 3; k++) {
      const L = D.cell[k][k];
      if (D.pbc[k] && L > 0) d[k] -= L * Math.round(d[k] / L);
    }
    return d;
  }
  const vec = (a, b) => mic([b.x - a.x, b.y - a.y, b.z - a.z]);
  const norm = v => Math.hypot(v[0], v[1], v[2]);
  function send(msg) {
    msg.view = viewer.getView();
    msg.t = Date.now();
    try {
      const host = window.parent.document.querySelector("." + BRIDGE + " input, ." + BRIDGE + " textarea");
      if (!host) throw new Error("no bridge");
      host.value = JSON.stringify(msg);
      host.dispatchEvent(new Event("input", {bubbles: true}));
      host.dispatchEvent(new Event("change", {bubbles: true}));
      return true;
    } catch (e) {
      return false;
    }
  }
  viewer.setClickable({}, true, a => {
    const mode = document.getElementById("mode").value;
    if (mode === "inspect") {
      label(`#${a.index} ${a.elem}` + (a.properties.mol ? ` L${a.properties.mol}` : ""), a);
      say(describe(a));
    } else if (mode === "measure") {
      picks.push(a);
      label(`${picks.length}: #${a.index}`, a);
      if (picks.length === 2) {
        const [p, q] = picks;
        const d = norm(vec(p, q));
        marks.push(viewer.addCylinder({start: {x: p.x, y: p.y, z: p.z}, end: {x: q.x, y: q.y, z: q.z},
          radius: 0.06, color: accent, dashed: true, fromCap: 1, toCap: 1}));
        say(`<b>Distance</b> #${p.index}–#${q.index}: ${d.toFixed(3)} Å (click a third atom for the angle at #${q.index})`);
      } else if (picks.length === 3) {
        const [p, q, r] = picks;
        const u = vec(q, p), v = vec(q, r);
        const cos = (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]) / (norm(u) * norm(v));
        const ang = Math.acos(Math.max(-1, Math.min(1, cos))) * 180 / Math.PI;
        say(`<b>Angle</b> #${p.index}–#${q.index}–#${r.index}: ${ang.toFixed(1)}°`);
        picks = [];
      } else {
        say(`Picked #${a.index}; click a second atom.`);
      }
    } else if (!D.editable) {
      say("Editing works on a generated model: click <b>Generate model</b> first.");
    } else if (mode === "remove") {
      if (!a.properties.mol) { say("Click an atom of a ligand to remove that ligand."); return; }
      say(`Removing ligand ${a.properties.mol}...`);
      if (!send({action: "remove", ligand: a.properties.mol, atom: a.index}))
        say(`Could not reach the app; type <code>remove_ligand(${a.properties.mol})</code> in the Python console.`);
    } else if (mode === "add") {
      if (a.properties.role !== "substrate") { say("Click a surface atom of the substrate."); return; }
      say(`Adding a ligand on atom ${a.index}...`);
      if (!send({action: "add", atom: a.index}))
        say(`Could not reach the app; type <code>add_ligand(${a.index})</code> in the Python console.`);
    }
    viewer.render();
  });
  for (const id of ["style", "lig", "sites", "cell"])
    document.getElementById(id).addEventListener("change", apply);
  document.getElementById("bg").value = state.bg;
  document.getElementById("bg").addEventListener("change", e => {
    state.bg = e.target.value; applyBackground(); apply(); });
  document.getElementById("mode").addEventListener("change", () => {
    picks = [];
    const m = document.getElementById("mode").value;
    say({inspect: "Click an atom to label it.", measure: "Click two atoms for a distance, three for an angle.",
         remove: D.editable ? "Click any atom of the ligand to remove." : "Generate a model first to edit it.",
         add: D.editable ? "Click a surface atom to place one more ligand there." : "Generate a model first to edit it."}[m]);
  });
  if (!D.editable) for (const o of document.querySelectorAll("#mode option[value=remove], #mode option[value=add]"))
    o.disabled = true;
  function setProjButtons() {
    document.getElementById("ortho").classList.toggle("on", state.ortho);
    document.getElementById("persp").classList.toggle("on", !state.ortho);
  }
  document.getElementById("ortho").addEventListener("click", () => { state.ortho = true; setProjButtons(); apply(); });
  document.getElementById("persp").addEventListener("click", () => { state.ortho = false; setProjButtons(); apply(); });
  document.getElementById("slice").addEventListener("input", e => {
    state.slice = +e.target.value; applySlice(viewer); viewer.render();
    say(state.slice ? `Slice: front ${state.slice}% cut away` + (state.thin ? " (thin slab)" : "") : "Slice off.");
  });
  document.getElementById("thin").addEventListener("change", e => {
    state.thin = e.target.checked; applySlice(viewer); viewer.render(); });
  // Zoom to the atoms only; the vacuum-padded cell would otherwise shrink the model.
  // Rods lie along z: start side-on (axis horizontal); channels look along y with z up.
  function look(kind) {
    viewer.setView([0, 0, 0, 0, 0, 0, 0, 1]);
    if (kind === "side") viewer.rotate(-90, "x");
    else if (kind === "angled") { viewer.rotate(-60, "x"); viewer.rotate(30, "z"); }
    else if (kind === "fit") {
      if (D.sideOn === "up") viewer.rotate(-90, "x");
      else if (D.sideOn === "angled") { viewer.rotate(-60, "x"); viewer.rotate(30, "z"); }
      else if (D.sideOn) viewer.rotate(90, "y");
    }
    viewer.zoomTo({model: model});
    slab0 = viewer.getSlab();
    applySlice(viewer);
    viewer.render();
  }
  document.getElementById("reset").addEventListener("click", () => look("fit"));
  document.getElementById("vtop").addEventListener("click", () => look("top"));
  document.getElementById("vside").addEventListener("click", () => look("side"));
  document.getElementById("vang").addEventListener("click", () => look("angled"));
  document.getElementById("clear").addEventListener("click", clearMarks);
  document.getElementById("png").addEventListener("click", () => {
    // Re-render off screen at 3x the on-screen size with the same view, style and slice.
    const host = document.getElementById("v");
    const scale = 3;
    const div = document.createElement("div");
    div.style.cssText = `position:fixed;left:-20000px;top:0;width:${host.clientWidth * scale}px;height:${host.clientHeight * scale}px`;
    document.body.appendChild(div);
    try {
      const big = makeViewer(div);
      styleViewer(big);
      big.setView(viewer.getView());
      const sl = viewer.getSlab();
      big.setSlab(sl.near, sl.far);
      big.render();
      const a = document.createElement("a");
      a.href = big.pngURI();
      a.download = "buddy_builder_view.png";
      document.body.appendChild(a); a.click(); a.remove();
      say(`Saved a ${host.clientWidth * scale} x ${host.clientHeight * scale} px PNG` +
          (state.bg === "transparent" ? " with a transparent background." : "."));
    } catch (e) {
      const a = document.createElement("a");
      a.href = viewer.pngURI(); a.download = "buddy_builder_view.png";
      document.body.appendChild(a); a.click(); a.remove();
    } finally {
      div.remove();
    }
  });
  if (!D.sites.length) document.getElementById("sites").disabled = true;
  if (!D.role.some(r => r !== "substrate")) document.getElementById("lig").disabled = true;
  document.getElementById("cell").checked = D.showCell;
  if (!D.cell) document.getElementById("cell").disabled = true;
  setProjButtons();
  applyBackground();
  styleViewer(viewer);
  if (D.view) { viewer.setView(D.view); slab0 = viewer.getSlab(); viewer.render(); } else look("fit");
}
</script></body></html>
"""


def _css_vars(palette: dict) -> str:
    return "".join(f"--{k.replace('_', '-')}:{v};" for k, v in palette.items())


THEMES = ("light", "dark", "buddy", "buddy2")
BACKGROUNDS = ("white", "black", "theme", "transparent")


def page_html(payload: dict, height: int = 460, theme: str = "buddy", bridge: str = "",
              background: str = None) -> str:
    """The viewer page. Periodic models (and unit cells) open in orthographic projection."""
    payload = dict(payload)
    payload.setdefault("ortho", bool(payload.get("showCell")))
    theme = theme if theme in THEMES else "light"
    if background not in BACKGROUNDS:
        background = "theme" if theme == "dark" else "white"
    payload["background"] = background
    data = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    return (_PAGE.replace("__URL__", THREEDMOL_URL).replace("__HEIGHT__", str(int(height)))
            .replace("__THEME__", theme).replace("__FONTFACE__", font_face_css())
            .replace("__FONT__", FONT_STACK).replace("__DISPLAY__", DISPLAY_STACK)
            .replace("__LIGHT__", _css_vars(LIGHT))
            .replace("__DARK__", _css_vars(DARK)).replace("__BUDDY__", _css_vars(BUDDY))
            .replace("__BUDDY2__", _css_vars(BUDDY2))
            .replace("__BRIDGE__", bridge).replace("__DATA__", data))


def iframe_html(payload: dict, height: int = 460, theme: str = "buddy", bridge: str = "",
                background: str = None) -> str:
    """Wrap the viewer page in an iframe suitable for an ``ipywidgets.HTML`` value.

    ``bridge`` is the CSS class of the hidden text widget that receives edit
    requests (remove / add ligand) from clicks in the viewer.
    """
    doc = html.escape(page_html(payload, height, theme, bridge, background), quote=True)
    return (f'<iframe class="nb-stage" srcdoc="{doc}" style="width:100%;height:{height + 150}px" '
            'title="3D structure viewer"></iframe>')


def placeholder_html(message: str, height: int = 460) -> str:
    return (f'<div class="nb-placeholder" style="height:{height + 150}px">'
            f'{html.escape(message)}</div>')
