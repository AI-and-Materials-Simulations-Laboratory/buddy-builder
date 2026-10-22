"""Certificate of Build: a PDF report of one model, with a fun first page.

Page 1 is a certificate (logo, the structure, its name, key numbers, build ID).
The following pages list what was built: structure and box, crystal, atoms
and composition, ligand chemistry, coverage, topology, validation and the
files in the export, plus three rendered views (top, side, angled).

Rendering uses matplotlib only (a dependency of ASE), so it works on nanoHUB
and offline. Pictures are drawn as depth-sorted spheres in Jmol colours.
"""

from __future__ import annotations

import hashlib
import io
import textwrap
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import __version__

ASSETS = Path(__file__).resolve().parent / "assets"
GREEN, GREEN_SOFT, PINK, YELLOW, CREAM, INK, MUTED = (
    "#2E7D55", "#DDF1E2", "#F7D9E1", "#FFF1B8", "#FBF6DF", "#2B3A2F", "#66705F")
PURPLE = "#4B2470"  # outer border of the certificate
MAX_DRAWN_ATOMS = 60000


# ----------------------------------------------------------------------------- rendering

def _rotation(view: str) -> np.ndarray:
    def rx(a):
        a = np.radians(a)
        return np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])

    def ry(a):
        a = np.radians(a)
        return np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    return {"top": np.eye(3), "side": rx(-90), "angled": ry(30) @ rx(-60)}[view]


def draw_structure(ax, assembly, view: str = "angled", show_box: bool = True, title: str = ""):
    """Depth-sorted ball picture of the model on a matplotlib axes (orthographic)."""
    from ase.data import atomic_numbers, covalent_radii
    from ase.data.colors import jmol_colors
    from matplotlib.collections import LineCollection

    pos = np.asarray(assembly.positions_A, dtype=float)
    elements = list(assembly.elements)
    n = len(pos)
    keep = np.arange(n)
    if n > MAX_DRAWN_ATOMS:  # very large models: thin the substrate interior
        keep = np.unique(np.concatenate([np.arange(0, n, max(1, n // MAX_DRAWN_ATOMS)),
                                         [i for i, r in enumerate(assembly.chemical_role)
                                          if r != "substrate"]]))
    rot = _rotation(view)
    centre = pos.mean(axis=0)
    p = (pos - centre) @ rot.T
    z_num = [atomic_numbers.get(e, 6) for e in elements]
    lig = np.array([r != "substrate" for r in assembly.chemical_role])
    radius = np.array([covalent_radii[z] for z in z_num]) * np.where(lig, 0.55, 0.95)
    colors = np.array([jmol_colors[z] for z in z_num])

    corners = None
    if show_box and assembly.cell_A is not None:
        cell = np.asarray(assembly.cell_A, dtype=float)
        o = -centre @ rot.T
        a, b, c = cell @ rot.T
        corners = np.array([o, o + a, o + b, o + c, o + a + b, o + a + c, o + b + c, o + a + b + c])
    pts = p[keep, :2] if corners is None else np.vstack([p[keep, :2], corners[:, :2]])
    lo, hi = pts.min(axis=0) - 2.0, pts.max(axis=0) + 2.0
    span = max(hi - lo)
    mid = (lo + hi) / 2
    ax.set_xlim(mid[0] - span / 2, mid[0] + span / 2)
    ax.set_ylim(mid[1] - span / 2, mid[1] + span / 2)
    ax.set_aspect("equal")
    ax.axis("off")
    fig = ax.figure
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    pts_per_A = bbox.width * 72.0 / span

    if corners is not None:
        edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7),
                 (6, 7)]
        ax.add_collection(LineCollection([[corners[i, :2], corners[j, :2]] for i, j in edges],
                                         colors="#3A4250", linewidths=1.1, zorder=3,
                                         alpha=0.85))
    bonds = [(i, j) for i, j, *_ in assembly.bond_list
             if (lig[i] or lig[j]) and np.linalg.norm(pos[i] - pos[j]) < 2.6]
    order = keep[np.argsort(p[keep, 2])]
    rank = np.empty(n)
    rank[order] = np.arange(len(order))
    if bonds:
        segs = [[p[i, :2], p[j, :2]] for i, j in bonds]
        ax.add_collection(LineCollection(segs, colors="#55606E", linewidths=1.0, zorder=1,
                                         rasterized=True))
    ax.scatter(p[order, 0], p[order, 1], s=(2 * radius[order] * pts_per_A) ** 2, c=colors[order],
               edgecolors="#1F2937", linewidths=0.25, zorder=2, rasterized=True)
    if title:
        ax.set_title(title, fontsize=10, color=MUTED)


def structure_png(assembly, view: str = "angled", size_in: float = 4.0, dpi: int = 160) -> bytes:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(size_in, size_in), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    draw_structure(ax, assembly, view)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", transparent=True)
    plt.close(fig)
    return buf.getvalue()


# ----------------------------------------------------------------------------- content

def build_id(config, assembly) -> str:
    h = hashlib.sha256()
    h.update(np.asarray(assembly.positions_A, dtype=np.float32).tobytes())
    h.update(repr(sorted(config.to_dict().items())).encode())
    return h.hexdigest()[:12].upper()


def describe(config, assembly) -> str:
    """One-line name of the structure, e.g. 'Gold (Au) - FCC sphere, 3.0 nm + 12 ligands'."""
    prov = assembly.provenance
    material = prov["material"]["structure_label"]
    g = config.geometry
    if config.material.source == "carbon":
        info = prov.get("carbon") or {}
        body = info.get("summary") or material
    elif g.shape == "sphere":
        body = f"{material} sphere, {g.diameter_nm:g} nm"
    elif g.shape == "rod":
        body = f"{material} nanorod, {g.diameter_nm:g} x {g.rod_length_nm:g} nm"
    elif g.shape == "slab":
        body = f"{material} ({' '.join(str(i) for i in g.miller)}) slab"
    else:
        ch = "slit channels" if g.channel_type == "slit" else "pores"
        body = f"{material} ({' '.join(str(i) for i in g.miller)}) membrane, {g.channel_count} {ch}"
    lig = prov.get("ligand")
    if lig and assembly.achieved_ligands:
        body += f" + {assembly.achieved_ligands} x {lig['formula_bound']}"
    return body


def _density_g_cm3(atoms) -> float:
    from ase.data import atomic_masses
    mass = float(sum(atomic_masses[z] for z in atoms.numbers))
    return mass * 1.66053906660 / atoms.get_volume() if atoms.get_volume() > 0 else 0.0


def sections(config, assembly, report, crystal_atoms=None, files=()) -> list:
    """[(heading, [(label, value)])] for the detail pages."""
    from .chemistry import TEMPLATES
    prov = assembly.provenance
    cell = np.asarray(assembly.cell_A, dtype=float)
    lengths = np.linalg.norm(cell, axis=1)
    comp = Counter(assembly.elements)
    sub_comp = Counter(e for e, r in zip(assembly.elements, assembly.chemical_role)
                       if r == "substrate")
    g = config.geometry
    pbc = "".join(a for a, p in zip("xyz", assembly.pbc) if p) or "none (isolated)"
    out = []
    rows = [("Structure", describe(config, assembly)),
            ("Shape", (prov.get("geometry") or {}).get("shape", g.shape)),
            ("Box (x, y, z)", " x ".join(f"{v / 10:.3f}" for v in lengths) + " nm"),
            ("Box angles", "orthogonal" if np.allclose(cell - np.diag(np.diag(cell)), 0, atol=1e-6)
             else "triclinic (see config.json)"),
            ("Periodic along", pbc),
            ("Box setting", g.box_shape + ("" if not g.box_lengths_nm else
                                           f", user lengths {tuple(g.box_lengths_nm)} nm"))]
    if g.shape in ("slab", "channel") and config.material.source != "carbon":
        rows.append(("Surface (hkl)", " ".join(str(i) for i in g.miller)))
        geo = prov.get("geometry") or {}
        if geo.get("termination"):
            rows.append(("Termination", str(geo["termination"])))
    out.append(("Structure and box", rows))

    m = prov["material"]
    rows = [("Material", m["structure_label"]), ("Source", m["source_type"].replace("_", " "))]
    if m.get("space_group"):
        rows.append(("Space group", str(m["space_group"])))
    if crystal_atoms is not None and config.material.source != "carbon":
        a, b, c, al, be, ga = crystal_atoms.cell.cellpar()
        rows += [("Lattice a, b, c", f"{a:.4f}, {b:.4f}, {c:.4f} A"),
                 ("Angles alpha, beta, gamma", f"{al:.2f}, {be:.2f}, {ga:.2f} deg"),
                 ("Atoms per cell", f"{len(crystal_atoms)} ({crystal_atoms.get_chemical_formula()})"),
                 ("Cell volume", f"{crystal_atoms.get_volume():.3f} A^3"),
                 ("Density", f"{_density_g_cm3(crystal_atoms):.3f} g/cm^3")]
    rows.append(("Lattice source", m.get("lattice_source") or "-"))
    if config.material.source == "materials_project":
        rows.append(("Citation", "Materials Project, CC BY 4.0; Jain et al., APL Mater. 1, "
                                 "011002 (2013)"))
    out.append(("Crystal", rows))

    rows = [("Total atoms", f"{assembly.n_atoms}"),
            ("Substrate atoms", f"{assembly.n_substrate}  ("
             + ", ".join(f"{k} {v}" for k, v in sorted(sub_comp.items())) + ")"),
            ("Ligand atoms", f"{assembly.n_atoms - assembly.n_substrate}"),
            ("Composition", ", ".join(f"{k} {v}" for k, v in sorted(comp.items()))),
            ("Total mass", f"{float(np.sum(assembly.masses_amu)):.1f} amu")]
    out.append(("Atoms", rows))

    lig = prov.get("ligand")
    if lig:
        t = TEMPLATES.get(lig["binding_template"])
        rows = [("SMILES", lig["smiles"]), ("Bound formula", lig["formula_bound"]),
                ("Attachment", t.label if t else lig["binding_template"]),
                ("Surface bond", "yes" if t and t.covalent_surface_bond else "no (physisorbed)"),
                ("Anchor distance", f"{lig['distance_A']:.2f} A"),
                ("Anchor species", ", ".join(lig.get("surface_species") or []) or "-"),
                ("Geometry", lig.get("optimization", "-"))]
        removed = lig.get("removed_atoms") or []
        if removed:
            rows.append(("Removed on binding", "; ".join(r["reason"] for r in removed)))
        out.append(("Ligand", rows))
        pl = prov.get("placement") or {}
        rows = [("Ligands placed / requested", f"{assembly.achieved_ligands} / {assembly.target_ligands}"),
                ("Region", config.coverage.region.replace("_", " ")),
                ("Arrangement", pl.get("distribution", config.coverage.distribution).replace("_", " ")),
                ("Coverage area", f"{assembly.surface_area_nm2:.3f} nm^2"),
                ("Achieved density", f"{assembly.achieved_density_nm2:.3f} per nm^2"),
                ("Closest nonbonded contact", "-" if report.min_nonbonded_distance_A is None
                 else f"{report.min_nonbonded_distance_A:.2f} A"),
                ("Seeds", ", ".join(f"{k} {v}" for k, v in (prov.get("seeds") or {}).items()) or "-")]
        edits = prov.get("manual_edits") or []
        if edits:
            rows.append(("Hand edits after Generate", str(len(edits))))
        out.append(("Coverage", rows))

    bond_types = Counter(b[2] for b in assembly.bond_list)
    rows = [("Bonds", str(len(assembly.bond_list))), ("Angles", str(len(assembly.angle_list))),
            ("Dihedrals", str(len(assembly.dihedral_list))),
            ("Bond types", ", ".join(f"{k} {v}" for k, v in bond_types.most_common(8)) or "-"),
            ("Force field", "none assigned (structure + topology export)")]
    out.append(("Topology", rows))

    rows = [("Status", report.status), ("Errors", str(len(report.errors))),
            ("Notes and warnings", str(len(report.warnings) + len(assembly.warnings)))]
    out.append(("Validation", rows))
    if files:
        out.append(("Files in the export", [("", ", ".join(files))]))
    return out


# ----------------------------------------------------------------------------- PDF

def _image(path: Path):
    import matplotlib.image as mpimg
    try:
        return mpimg.imread(str(path))
    except (OSError, ValueError):
        return None


FONT_DIR = ASSETS / "fonts"
DISPLAY = "Fredoka"  # titles and headings, as in the GUI (bundled, SIL OFL)
BODY = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]  # other text, as in the GUI
LIGHT_YELLOW = "#FFF3C4"   # banner of the detail pages
AWNING = ("#F3C7B1", "#FFF9EC")  # peach and cream stripes, as on the Buddy poster
AWNING_EDGE = "#E7B7A2"
PAGE = "#FFFDF6"


def _register_fonts() -> bool:
    """Make the bundled Fredoka known to matplotlib; False if it cannot be loaded."""
    from matplotlib import font_manager
    ok = False
    for path in sorted(FONT_DIR.glob("Fredoka-*.ttf")):
        try:
            font_manager.fontManager.addfont(str(path))
            ok = True
        except Exception:  # a broken font file must not stop the report
            pass
    return ok


def _banner(fig, title: str, bid: str):
    """Poster-style header for the detail pages: awning stripes over a light yellow band."""
    from matplotlib.patches import Rectangle
    ax = fig.add_axes([0, 0.885, 1, 0.115])
    ax.axis("off")
    ax.set_xlim(0, 8.5)
    ax.set_ylim(0, 1.265)
    ax.add_patch(Rectangle((0, 0), 8.5, 0.83, facecolor=LIGHT_YELLOW, edgecolor="none"))
    stripe = 8.5 / 12
    for k in range(12):
        ax.add_patch(Rectangle((k * stripe, 0.9), stripe, 0.365, facecolor=AWNING[k % 2],
                               edgecolor="none"))
    ax.add_patch(Rectangle((0, 0.83), 8.5, 0.07, facecolor=AWNING_EDGE, edgecolor="none"))
    ax.text(0.6, 0.38, title, fontsize=19, color=PURPLE, family=DISPLAY, fontweight=600,
            va="center")
    ax.text(7.9, 0.38, f"Build ID {bid}", fontsize=9, color=MUTED, ha="right", va="center")
    return ax


def certificate_pdf(config, assembly, report, crystal_atoms=None, files=(), name: str = "") -> bytes:
    """The Certificate of Build as PDF bytes; ``name`` is the save-as name of the build."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    have_display = _register_fonts()
    from matplotlib import font_manager
    installed = {f.name for f in font_manager.fontManager.ttflist}
    body = next((f for f in BODY if f in installed), "DejaVu Sans")  # Arial when present
    rc = {"font.family": body, "pdf.fonttype": 42}
    with plt.rc_context(rc):
        return _certificate(config, assembly, report, crystal_atoms, files, name,
                            DISPLAY if have_display else body)


def _certificate(config, assembly, report, crystal_atoms, files, name, display) -> bytes:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import FancyBboxPatch

    now = datetime.now(timezone.utc)
    bid = build_id(config, assembly)
    file_name = name or "Buddy build"
    name = describe(config, assembly)
    buf = io.BytesIO()
    with PdfPages(buf, metadata={"Title": "Certificate of Build", "Author": "Buddy Builder",
                                 "Subject": name, "Creator": f"nanobuilder {__version__}"}) as pdf:
        # ---------------------------------------------------------------- page 1
        fig = plt.figure(figsize=(11, 8.5))
        fig.patch.set_facecolor(CREAM)
        bg = fig.add_axes([0, 0, 1, 1])
        bg.axis("off")
        bg.set_xlim(0, 11)
        bg.set_ylim(0, 8.5)
        for inset, color, lw in ((0.25, PURPLE, 5), (0.42, PINK, 3), (0.52, YELLOW, 1.5)):
            bg.add_patch(FancyBboxPatch((inset, inset), 11 - 2 * inset, 8.5 - 2 * inset,
                                        boxstyle="round,pad=0,rounding_size=0.25", fill=False,
                                        edgecolor=color, linewidth=lw))
        logo = _image(ASSETS / "buddy_certificate.png")
        if logo is None:
            logo = _image(ASSETS / "buddy_logo_large.png")
        if logo is not None:
            lax = fig.add_axes([0.065, 0.655, 0.18, 0.245])
            lax.imshow(logo)
            lax.axis("off")
        bg.text(6.2, 7.35, "Certificate of Build", ha="center", va="center", fontsize=42,
                color=PURPLE, family=display, fontweight=600)
        bg.text(6.2, 6.72, "Buddy Builder proudly certifies that it assembled", ha="center",
                fontsize=14, color=MUTED, style="italic")
        wrapped = textwrap.fill(name, 40)
        bg.text(1.1, 5.05, wrapped, ha="left", va="center", fontsize=20, color=PURPLE,
                family=display, fontweight=500, linespacing=1.3)
        facts = [("Atoms", f"{assembly.n_atoms:,}", 1.25, 3.9, 2.45, GREEN_SOFT),
                 ("Ligands", f"{assembly.achieved_ligands} of {assembly.target_ligands}",
                  3.95, 3.9, 2.45, PINK),
                 ("Box", " x ".join(f"{v / 10:.2f}" for v in
                                    np.linalg.norm(np.asarray(assembly.cell_A), axis=1)) + " nm",
                  1.25, 2.95, 5.15, YELLOW)]
        for label, value, x, y, width, color in facts:
            bg.add_patch(FancyBboxPatch((x - 0.15, y - 0.32), width, 0.78,
                                        boxstyle="round,pad=0,rounding_size=0.12",
                                        facecolor=color, edgecolor="none"))
            bg.text(x, y + 0.2, label.upper(), fontsize=8.5, color=MUTED, fontweight="bold")
            bg.text(x, y - 0.12, value, fontsize=15, color=INK, family=display, fontweight=500)
        # The structure itself, on the certificate's own background (transparent picture).
        sax = fig.add_axes([0.58, 0.215, 0.36, 0.47])
        sax.patch.set_alpha(0.0)
        draw_structure(sax, assembly, "angled", show_box=False)
        bg.text(1.1, 1.25, f"Issued {now:%B %d, %Y}  ·  Build ID {bid}", fontsize=10, color=MUTED)
        bg.text(1.1, 0.95, f"nanobuilder {__version__}  ·  structure and topology for MD "
                           "(no force field)", fontsize=9, color=MUTED)
        bg.plot([7.0, 9.6], [1.15, 1.15], color=INK, linewidth=0.8)
        bg.text(8.3, 1.25, "Buddy", ha="center", fontsize=24, color=PURPLE, family=display,
                fontweight=500)
        bg.text(8.3, 0.88, textwrap.shorten(file_name, 44, placeholder="..."), ha="center",
                fontsize=9.5, color=MUTED)
        pdf.savefig(fig, facecolor=fig.get_facecolor(), dpi=200)
        plt.close(fig)

        # ---------------------------------------------------------------- detail pages
        rows_all = sections(config, assembly, report, crystal_atoms, files)
        lines = []
        for heading, rows in rows_all:
            lines.append(("h", heading))
            for label, value in rows:
                wrapped = textwrap.wrap(str(value), 62) or [""]
                lines.append(("r", label, wrapped[0]))
                for extra in wrapped[1:]:
                    lines.append(("r", "", extra))
        top = 9.45
        pages, page, y = [], [], top
        for item in lines:  # paginate by height: headings 0.38 in, rows 0.2 in
            h = 0.38 if item[0] == "h" else 0.2
            if y - h < 0.9 and page:
                pages.append(page)
                page, y = [], top
            page.append(item)
            y -= h
        if page:
            pages.append(page)
        for number, chunk in enumerate(pages, start=2):
            fig = plt.figure(figsize=(8.5, 11))
            fig.patch.set_facecolor(PAGE)
            _banner(fig, "Certificate of Build: details", bid)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.axis("off")
            ax.set_xlim(0, 8.5)
            ax.set_ylim(0, 11)
            y = top
            for item in chunk:
                if item[0] == "h":
                    y -= 0.08
                    ax.text(0.6, y, item[1], fontsize=13, color=PURPLE, family=display,
                            fontweight=600)
                    ax.plot([0.6, 7.9], [y - 0.08, y - 0.08], color=PINK, linewidth=1.4)
                    y -= 0.30
                else:
                    ax.text(0.7, y, item[1], fontsize=9, color=MUTED)
                    ax.text(3.0, y, item[2], fontsize=9, color=INK, family="monospace"
                            if item[1] in ("SMILES",) else None)
                    y -= 0.2
            ax.text(4.25, 0.4, f"Buddy Builder · page {number}", fontsize=8, color=MUTED,
                    ha="center")
            pdf.savefig(fig, facecolor=fig.get_facecolor(), dpi=200)
            plt.close(fig)

        # ---------------------------------------------------------------- views
        fig = plt.figure(figsize=(8.5, 11))
        fig.patch.set_facecolor(PAGE)
        _banner(fig, "Views", bid)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis("off")
        ax.text(0.07, 0.865, "Orthographic pictures; the box is drawn in grey.", fontsize=10,
                color=MUTED, transform=ax.transAxes)
        for k, (view, rect) in enumerate((("top", [0.06, 0.44, 0.42, 0.40]),
                                          ("side", [0.52, 0.44, 0.42, 0.40]),
                                          ("angled", [0.18, 0.03, 0.64, 0.40]))):
            vax = fig.add_axes(rect)
            vax.patch.set_alpha(0.0)
            draw_structure(vax, assembly, view)
            vax.set_title(f"{view} view", fontsize=12, color=PURPLE, family=display,
                          fontweight=500)
        pdf.savefig(fig, facecolor=fig.get_facecolor(), dpi=200)
        plt.close(fig)
    return buf.getvalue()


__all__ = ["certificate_pdf", "structure_png", "draw_structure", "describe", "build_id", "sections"]
