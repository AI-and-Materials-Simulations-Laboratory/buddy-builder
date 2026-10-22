"""ipywidgets GUI for the nanoHUB tool (spec section 4 and stage G).

Everything scientific happens in ``nanobuilder.pipeline``; this module only
turns widget values into a :class:`RunConfig`, calls the pipeline and shows
the results. Downloads are data-URI links, which work in classic Notebook,
AppMode and JupyterLab without a file server route.
"""

from __future__ import annotations

import base64
import html
import os
import shutil
import traceback
from pathlib import Path

import ipywidgets as w
import numpy as np
from IPython.display import display

from . import __version__
from .chemistry import TEMPLATES, template_distance
from .console import Console
from .crystal import build_crystal, default_preset_id, preset, preset_choices
from . import material_catalog as catalog
from . import mp_library
from .editing import add_ligand_at, remove_ligands
from .examples import EXAMPLES, example_upload
from .gui_carbon import CarbonPanel
from .importers import atoms_to_cif_bytes
from .ligand_library import library_choices, library_entry, suits
from .ligands import anchor_choices, default_anchor, ligand_summary, svg_2d
from .pipeline import (attachment_options, build_substrate, export_package, material_crystal,
                       run_pipeline, validate_inputs)
from .schema import (BuilderError, CarbonConfig, CoverageConfig, ExportConfig, GeometryConfig,
                     LigandConfig, MaterialConfig, RunConfig)
from .shapes import build_geometry, list_terminations
from .theme import stylesheet
from .viewer import (assembly_payload, iframe_html, placeholder_html, substrate_payload,
                     unit_cell_payload)

STYLE = {"description_width": "170px"}
MAX_DOWNLOAD_MB = 40
LOGO = Path(__file__).resolve().parent / "assets" / "buddy_logo.png"
APP_NAME = "Buddy Builder"


# Each widget needs its own Layout object: a shared one would make hiding one
# widget (layout.display = "none") hide every widget that shares it. Widths are
# capped, not fixed, so inputs shrink with a narrow column instead of clipping.
def wide():
    return w.Layout(width="100%", max_width="440px")


def narrow():
    return w.Layout(width="290px", max_width="100%")


def row(*children, **kw):
    """A wrapping row: side-by-side when there is room, stacked when there is not."""
    return w.HBox(list(children), layout=w.Layout(flex_flow="row wrap", align_items="center",
                                                  **kw))


def section(title: str, advanced: bool = False):
    cls = "nb-adv-h" if advanced else "nb-sec-h"
    return w.HTML(f"<div class='{cls}'>{html.escape(title)}</div>")


def choices(options, value=None, icons=None, tooltips=None, width="auto"):
    """Icon toggle buttons for a small set of choices (instead of a dropdown)."""
    tb = w.ToggleButtons(options=options, value=value if value is not None else options[0][1],
                         style={"button_width": width})
    if icons:
        tb.icons = list(icons)
    if tooltips:
        tb.tooltips = list(tooltips)
    tb.add_class("nb-choices")
    return tb


def _logo_uri() -> str:
    try:
        return "data:image/png;base64," + base64.b64encode(LOGO.read_bytes()).decode("ascii")
    except OSError:
        return ""


def _uploaded_file(widget):
    """(name, bytes) from a FileUpload, for both ipywidgets 7 and 8 value formats."""
    value = widget.value
    if isinstance(value, dict):  # ipywidgets 7: {filename: {"content": bytes, ...}}
        if not value:
            return None
        name, item = next(iter(value.items()))
        return name, bytes(item["content"])
    if value:  # ipywidgets 8: ({"name": ..., "content": memoryview, ...},)
        item = value[0]
        return item["name"], bytes(item["content"])
    return None


def _msg(kind: str, lines) -> str:
    """Status block; the theme colours it by kind (ok, info, warning, error)."""
    if isinstance(lines, str):
        lines = [lines]
    lines = [x for x in lines if str(x).strip()]
    if not lines:
        return ""
    items = "".join(f"<li>{html.escape(str(x))}</li>" for x in lines)
    return f'<div class="nb-msg nb-msg-{kind}"><ul>{items}</ul></div>'


def _download_link(name: str, data: bytes, mime: str, label: str, primary: bool = False) -> str:
    size = f"{len(data) / 1e6:.1f} MB" if len(data) >= 1e6 else f"{len(data) / 1e3:.0f} kB"
    if len(data) > MAX_DOWNLOAD_MB * 1e6:
        return f'<span class="nb-size">{html.escape(label)}: {size}, too large for a browser link</span>'
    b64 = base64.b64encode(data).decode("ascii")
    cls = ' class="nb-primary"' if primary else ""
    return (f'<a{cls} download="{html.escape(name)}" href="data:{mime};base64,{b64}">'
            f'{html.escape(label)}<span class="nb-size">{size}</span></a>')


# (pattern in the error, why it failed, how to fix it), first match wins
_ADVICE = [
    ("pore needs the box", "The pore is too wide for this membrane: it needs room for 3 atomic "
     "spacings on each side.", "Shape tab: use a smaller pore diameter, or make the box larger in x "
     "and y."),
    ("pores of", "There are too many pores (or they are too large) to keep 3 atomic spacings "
     "between them and from the box edges.", "Shape tab: use fewer or smaller pores, or a larger box."),
    ("enlarge the box", "The structure, ligands included, is bigger than the simulation box.",
     "Shape tab: make the box larger (the message gives the minimum) or make the structure "
     "smaller."),
    ("outside the box", "Some atoms stick out of the simulation box.",
     "Shape tab: make the box larger, or use fewer or shorter ligands."),
    ("would have about", "The structure would have more atoms than the tool builds at once.",
     "Use a smaller box, fewer layers or a smaller particle."),
    ("atoms; the limit", "The structure would have more atoms than the tool builds at once.",
     "Use a smaller box, fewer layers or a smaller particle."),
    ("upload a crystallographic file", "No crystal file has been uploaded.",
     "Material tab: choose Upload file and pick a CIF or POSCAR, or pick a preset or library "
     "material."),
    ("not in the preloaded", "That structure is not in the preloaded library.",
     "Material tab: pick a material from Browse materials, or upload its CIF."),
    ("smiles", "The ligand SMILES could not be read.",
     "Ligand tab: check the SMILES, or pick a ligand from the list."),
    ("no surface anchor species", "No surface element is allowed to hold a ligand.",
     "Ligand tab: tick at least one element under Surface anchor species."),
    ("anchor", "The chosen anchor atom cannot form this attachment.",
     "Ligand tab: pick another anchor atom (the drawing shows the numbers) or another "
     "attachment mode."),
    ("miller", "Those Miller indices do not describe a plane.",
     "Shape tab: enter Miller indices that are not all zero, e.g. 1 1 1."),
    ("termination", "That surface termination does not exist for this plane and thickness.",
     "Shape tab: choose one of the listed terminations."),
    ("slit width", "The slit is narrower than two atomic planes.",
     "Shape tab: use a gap of at least 0.2 nm."),
    ("atomic spacing", "The pore is smaller than the distance between atoms.",
     "Shape tab: use a larger pore diameter."),
    ("chirality", "Those nanotube indices are not valid.",
     "Shape tab: use whole numbers with 0 <= m <= n."),
    ("needs a carbon nanomaterial", "That attachment mode only works on carbon materials.",
     "Ligand tab: choose Physisorption or Chemisorption for crystals."),
    ("au-thiolate", "Au-thiolate only works on the gold preset with a thiol (-SH) anchor.",
     "Ligand tab: choose Chemisorption or Physisorption, or use a thiol ligand on gold."),
    ("implausible length", "A bond in the model has an unrealistic length.",
     "Ligand tab (Advanced): check the anchor distance; or choose another attachment mode."),
    ("overlap", "Two atoms ended up on top of each other.",
     "Coverage tab: place fewer ligands, or raise the clash tolerance (Advanced)."),
    ("valence", "An atom in the model has too many or too few bonds.",
     "Try another attachment mode or region, or lower the oxygen content."),
]


def failure_advice(message: str) -> tuple:
    """(why, how to fix) for an error message, in plain words."""
    text = str(message).lower()
    for pattern, why, fix in _ADVICE:
        if pattern in text:
            return why, fix
    first = str(message).split(". ")[0].rstrip(".")
    return (first + ".", "Click Check inputs to list every problem, adjust that setting, and "
            "Generate again. Reset brings back a build that works.")


class BuilderApp:
    """The complete tool GUI. ``BuilderApp().show()`` in a notebook cell."""

    def __init__(self, assistant: bool = True):
        self._want_assistant = assistant  # the Ask Buddy chat (False for the hidden checking copy)
        self.assistant = None
        self.result = None
        self.history = []  # results before each hand edit (Undo)
        self.package = None
        self.crystal = None
        self.upload = None  # (filename, bytes)
        self._stale = False
        self._quiet = False
        self._build_widgets()
        self._wire()
        self.reset()

    # ------------------------------------------------------------------ layout

    def _build_widgets(self):
        # Tab 1 - material
        # Type-to-search picker: grouped common materials, or anything typed -> Materials Project.
        # Hover menu: groups open sub-menus of materials; a click picks one (through menu_bridge).
        self.menu_bridge = w.Text(value="", continuous_update=True, layout=w.Layout(display="none"))
        self._menu_class = f"bb-menu-{id(self):x}"
        self.menu_bridge.add_class(self._menu_class)
        self.mat_menu = w.HTML(self._menu_html())
        self.mat_search = w.Combobox(options=catalog.labels(), ensure_option=False, continuous_update=False,
                                     placeholder="Au, AuCu3, TiO2, anatase, graphene oxide, C60 ...",
                                     description="Find material", style=STYLE,
                                     layout=w.Layout(width="100%", max_width="420px"))
        self.mat_go = w.Button(description="Find", icon="search",
                               layout=w.Layout(width="90px"))
        self.mat_hint = w.HTML("<small>Hover <b>Browse materials</b> and pick one, or type a name, "
                               "formula (TiO2) or carbon word (GO, CNT, C60) and press Enter.</small>")
        self.source = choices([("Preset", "preset"), ("MP library", "materials_project"),
                               ("Upload file", "uploaded_file"), ("Carbon", "carbon")],
                              icons=["flask", "book", "upload", "circle-o"],
                              tooltips=["Built-in crystal presets (Au, Ag, Cu, Fe, Si, Mg, Ge, NaCl)",
                                        "Preloaded Materials Project structures (offline)",
                                        "Your own CIF or POSCAR file",
                                        "Graphene, graphene oxide, nanotubes, fullerenes"])
        self.mp_pick = w.Dropdown(options=[(e["name"], e["material_id"])
                                           for e in mp_library.entries().values()],
                                  description="Structure", style=STYLE, layout=wide())
        self.mp_summary = w.HTML()
        self.mp_box = w.VBox([self.mp_pick, self.mp_summary])
        self.preset = w.Dropdown(options=preset_choices(), description="Preset",
                                 style=STYLE, layout=wide())
        self.lat_a = w.BoundedFloatText(min=0.5, max=50, step=0.01, description="a (Å)",
                                        style=STYLE, layout=narrow())
        self.lat_c = w.BoundedFloatText(min=0.5, max=50, step=0.01, description="c (Å)",
                                        style=STYLE, layout=narrow())

        self.upload_format = w.Dropdown(options=[("CIF (.cif)", "cif"), ("VASP POSCAR/CONTCAR (P1)", "vasp")],
                                        description="Upload format", style=STYLE, layout=wide())
        self.upload_widget = w.FileUpload(accept=".cif", multiple=False, description="Choose file")
        self.upload_summary = w.HTML("<i>No file uploaded yet.</i>")
        self.preview_cell_btn = w.Button(description="Preview unit cell", icon="cube")
        self.lattice_adv = w.VBox([section("Advanced: lattice constants", True), row(self.lat_a, self.lat_c)])
        self.preset_box = w.VBox([self.preset, self.lattice_adv])
        self.cell_info = w.HTML()
        self.upload_box = w.VBox([self.upload_format, row(w.Label("Crystal file"), self.upload_widget),
                                  self.upload_summary])
        self.carbon = CarbonPanel()
        tab_material = w.VBox([section("Browse materials"), self.mat_menu, self.menu_bridge,
                               row(self.mat_search, self.mat_go), self.mat_hint,
                               section("Source"), self.source,
                               self.preset_box, self.mp_box, self.upload_box,
                               self.carbon.material_box, section("Unit cell"), self.cell_info,
                               self.preview_cell_btn])

        # Tab 2 - shape & box
        self.shape = choices([("Sphere", "sphere"), ("Nanorod", "rod"), ("Slab", "slab"),
                              ("Channel", "channel")],
                             icons=["circle", "battery-full", "square", "pause"],
                             tooltips=["Finite nanoparticle cut from the crystal",
                                       "Cylinder with flat or round caps along a lattice direction",
                                       "Surface periodic in x and y, vacuum along z",
                                       "Slit pores between slab walls, or cylindrical pores through a "
                                       "membrane"])
        self.diameter = w.BoundedFloatText(min=0.3, max=30, step=0.1, description="Particle diameter (nm)",
                                           style=STYLE, layout=narrow())
        self.slab_len = w.BoundedFloatText(min=0.3, max=50, step=0.5, description="Slab length x (nm)",
                                           style=STYLE, layout=narrow())
        self.slab_wid = w.BoundedFloatText(min=0.3, max=50, step=0.5, description="Slab width y (nm)",
                                           style=STYLE, layout=narrow())
        self.thick_mode = w.Dropdown(options=[("Layer count", "layers"), ("Thickness (nm)", "thickness")],
                                     description="Thickness by", style=STYLE,
                                     layout=w.Layout(width="330px", max_width="100%"))
        self.layers = w.BoundedIntText(min=1, max=60, description="Layers (atomic planes)", style=STYLE,
                                       layout=narrow())
        self.thickness = w.BoundedFloatText(min=0.1, max=20, step=0.1, description="Thickness (nm)",
                                            style=STYLE, layout=narrow())
        self.h = w.IntText(description="Miller (h k l)", style={"description_width": "170px"},
                           layout=w.Layout(width="230px"))
        self.k = w.IntText(layout=w.Layout(width="56px"))
        self.l = w.IntText(layout=w.Layout(width="56px"))
        self.termination = w.Dropdown(description="Slab termination", style=STYLE, layout=wide())
        self.term_info = w.HTML()
        self.boundary = w.HTML()
        self.vacuum = w.BoundedFloatText(min=5, max=200, step=1, description="Vacuum padding (Å)",
                                         style=STYLE, layout=narrow())
        self.max_atoms = w.BoundedIntText(min=100, max=200000, step=1000, description="Max atoms",
                                          style=STYLE, layout=narrow())
        self.preview_sub_btn = w.Button(description="Preview bare substrate", icon="eye")
        self.rod_len = w.BoundedFloatText(min=0.3, max=100, step=0.5, description="Rod length, tip to tip (nm)",
                                          style=STYLE, layout=narrow())
        self.rod_u = w.IntText(layout=w.Layout(width="56px"), tooltip="u")
        self.rod_v = w.IntText(layout=w.Layout(width="56px"), tooltip="v")
        self.rod_w = w.IntText(layout=w.Layout(width="56px"), tooltip="w")
        self.rod_caps = choices([("Round caps", "hemispherical"), ("Flat caps", "flat")])
        self.rod_info = w.HTML("<small>Axis is a lattice direction of the conventional cell; Au "
                               "nanorods typically grow along [0 0 1].</small>")
        self.sphere_box = w.VBox([self.diameter])
        rod_axis_label = w.Label("Rod axis [u v w]", layout=w.Layout(width="174px"))
        self.rod_box = w.VBox([self.rod_len, row(rod_axis_label, self.rod_u, self.rod_v, self.rod_w),
                               self.rod_caps, self.rod_info])
        # channels (crystal): stacked slit pores between slab walls, or cylindrical pores
        self.channel_type = choices([("Slit channels", "slit"), ("Cylindrical pores", "pore")],
                                    icons=["bars", "dot-circle-o"],
                                    tooltips=["Parallel walls with gaps between them (multi-layer "
                                              "channels: set the number of channels)",
                                              "Round pores through a membrane slab, along x"])
        self.channel_width = w.BoundedFloatText(min=0.3, max=20, step=0.1, description="Channel width (nm)",
                                                style=STYLE, layout=narrow())
        self.channel_count = w.BoundedIntText(min=1, max=20, description="Number of channels",
                                              style=STYLE, layout=narrow())
        self.channel_info = w.HTML()
        self.channel_box = w.VBox([section("Channel"), self.channel_type,
                                   row(self.channel_width, self.channel_count), self.channel_info])
        self.slab_box = w.VBox([row(self.slab_len, self.slab_wid),
                                row(self.thick_mode, self.layers, self.thickness),
                                row(self.h, self.k, self.l), self.termination, self.term_info])
        # simulation box: the lengths the user enters decide the periodic size of slabs,
        # channels and periodic carbon; finite models must fit inside.
        self.box_shape = choices([("Rectangular", "rectangular"), ("Cube", "cube"), ("Auto", "auto")],
                                 icons=["square-o", "cube", "magic"],
                                 tooltips=["Box with your x, y and z lengths",
                                           "Box with one edge length",
                                           "Fit the model and pad it with 3 nm of vacuum"])
        self.box_info = w.HTML()
        box_style = {"description_width": "52px"}
        self.box_x = w.BoundedFloatText(min=0.5, max=200, step=0.5, value=8.0, description="x (nm)",
                                        style=box_style, layout=w.Layout(width="150px"))
        self.box_y = w.BoundedFloatText(min=0.5, max=200, step=0.5, value=8.0, description="y (nm)",
                                        style=box_style, layout=w.Layout(width="150px"))
        self.box_z = w.BoundedFloatText(min=0.5, max=200, step=0.5, value=8.0, description="z (nm)",
                                        style=box_style, layout=w.Layout(width="150px"))
        self.box_edge = w.BoundedFloatText(min=0.5, max=200, step=0.5, value=8.0,
                                           description="Edge (nm)", style={"description_width": "70px"},
                                           layout=w.Layout(width="170px"))
        self.box_xyz = row(self.box_x, self.box_y, self.box_z)
        tab_geometry = w.VBox([section("Simulation box"), self.box_shape, self.box_xyz, self.box_edge,
                               self.box_info, section("Shape"), self.shape, self.sphere_box,
                               self.rod_box, self.slab_box, self.channel_box,
                               self.carbon.geometry_box, self.boundary, self.preview_sub_btn])

        # Tab 4 - ligand
        self.library = w.Dropdown(options=[("Custom SMILES", "")] + library_choices(False),
                                  description="Ligand", style=STYLE, layout=wide())
        self.library_info = w.HTML()
        self.conformer_mode = w.Dropdown(options=[("Extended chains (all-trans)", "extended"),
                                                  ("RDKit conformer as embedded", "etkdg")],
                                         description="Chain conformation", style=STYLE, layout=wide())
        self.smiles = w.Text(description="SMILES", placeholder="Enter ligand SMILES", style=STYLE,
                             layout=wide(), continuous_update=False)
        self.smiles_btn = w.Button(description="Check", icon="check")
        self.smiles_info = w.HTML()
        self.smiles_svg = w.HTML()
        self.smiles_svg.add_class("nb-svg")
        self.anchor = w.Dropdown(description="Anchor atom", style=STYLE, layout=narrow())
        self.template = w.Dropdown(description="Attachment mode", style=STYLE, layout=wide())
        self.template_info = w.HTML()
        self.site_type = w.Dropdown(options=[("On-top", "on_top")], description="Adsorption site",
                                    style=STYLE, layout=narrow(), disabled=True)
        self.species_box = w.HBox(layout=w.Layout(flex_flow="row wrap"))
        self.species_warning = w.HTML()
        self.species_checks = {}
        self.distance = w.BoundedFloatText(min=1.0, max=6.0, step=0.05, description="Anchor distance (Å)",
                                           style=STYLE, layout=narrow())
        self.orientation = w.Dropdown(options=[("Outward, random azimuth", "random_azimuth"),
                                               ("Outward, fixed azimuth", "outward")],
                                      description="Orientation", style=STYLE, layout=wide())
        self.conformers = w.BoundedIntText(min=1, max=20, description="Conformer count", style=STYLE,
                                           layout=narrow())
        self.lig_seed = w.IntText(description="Random seed (conformers)", style=STYLE, layout=narrow())
        self.species_row = row(w.Label("Surface anchor species", layout=w.Layout(width="174px")),
                               self.species_box)
        self.ligand_adv = w.VBox([section("Advanced", True), self.distance, self.orientation,
                                  self.conformer_mode, self.conformers, self.lig_seed])
        tab_ligand = w.VBox([section("Ligand"), self.library, self.library_info,
                             row(self.smiles, self.smiles_btn), self.smiles_info, self.smiles_svg,
                             self.anchor, section("Attachment"), self.template, self.template_info,
                             self.species_row, self.species_warning, self.ligand_adv])

        # Tab 5 - coverage
        self.cov_mode = choices([("Exact number", "count"), ("Surface density", "density")],
                                icons=["hashtag", "th"])
        self.count = w.BoundedIntText(min=0, max=10000, description="Number of ligands", style=STYLE,
                                      layout=narrow())
        self.density = w.BoundedFloatText(min=0, max=20, step=0.05, description="Density (ligands per nm²)",
                                          style=STYLE, layout=narrow())
        self.region = w.Dropdown(description="Placement region", style=STYLE, layout=wide())
        self.area_info = w.HTML()
        self.distribution = choices([("Fibonacci (ordered)", "fibonacci"),
                                     ("Evenly spaced", "evenly_spaced"),
                                     ("Random", "distributed_random")],
                                    icons=["sun-o", "th-large", "random"],
                                    tooltips=["Golden-spiral points on spheres, golden-angle helix on "
                                              "cylinders, golden-ratio lattice on flat faces; each "
                                              "point takes the closest free site. Same result every run.",
                                              "Greedy: each ligand on the site farthest from those placed",
                                              "Random sites (seeded)"])
        self.min_sep = w.BoundedFloatText(min=0, max=50, step=0.5, description="Min site separation (Å)",
                                          style=STYLE, layout=narrow())
        self.clash = w.FloatSlider(min=0.3, max=1.2, step=0.05, description="Clash tolerance (alpha)",
                                   style=STYLE, layout=wide(), readout_format=".2f")
        self.max_attempts = w.BoundedIntText(min=1, max=200000, description="Max placement attempts",
                                             style=STYLE, layout=narrow())
        self.cov_seed = w.IntText(description="Random seed (placement)", style=STYLE, layout=narrow())
        clash_note = w.HTML("<small>A trial ligand is rejected when any nonbonded pair is closer than "
                            "alpha x (sum of van der Waals radii). The anchor-to-site pair is exempt.</small>")
        self.coverage_adv = w.VBox([section("Advanced", True), self.min_sep, self.clash, clash_note,
                                    self.max_attempts, self.cov_seed])
        tab_coverage = w.VBox([section("How many"), self.cov_mode, self.count, self.density,
                               section("Where"), self.region, self.area_info,
                               section("Arrangement"), self.distribution, self.coverage_adv])

        # Tab 6 - export
        self.t_hoomd = w.Checkbox(description="HOOMD-blue (.gsd)", indent=False)
        self.t_lammps = w.Checkbox(description="LAMMPS (.data, atom_style full)", indent=False)
        self.t_generic = w.Checkbox(description="Generic (.extxyz)", indent=False)
        self.t_pdb = w.Checkbox(description="PDB (.pdb, for AMBER / GROMACS)", indent=False)
        self.t_mol2 = w.Checkbox(description="MOL2 (.mol2, suggested GAFF2 types for AMBER)",
                                 indent=False)
        self.residue_scheme = w.Dropdown(
            options=[("Per group: GRA / OH / EPX / ... (AMBER GO dataset)", "groups"),
                     ("Single residue GRO for the substrate", "single")],
            description="MOL2 carbon residues", style=STYLE, layout=wide())
        # Kept for configurations that set it; not shown (structure + topology is what is exported).
        self.fidelity = w.Dropdown(options=[("Structure/topology (Level B)", "structure_topology"),
                                            ("Simulation-ready (validated templates only)",
                                             "simulation_ready_validated")],
                                   description="Output fidelity", style=STYLE, layout=wide())
        self.include_h = w.Checkbox(description="Include explicit hydrogens (always on)", value=True,
                                    disabled=True, indent=False)
        self.freeze = w.Checkbox(description="Freeze substrate (metadata + LAMMPS template fix)", indent=False)
        self.charges = w.Dropdown(options=[("Unassigned", "unassigned")], description="Atom charges",
                                  style=STYLE, layout=narrow(), disabled=True)
        units = w.HTML("<small><b>Units:</b> lengths in Ångström, masses in amu in every file. "
                       "LAMMPS: <code>units real</code> or <code>metal</code> both read these values; "
                       "HOOMD units are application-defined. topology.json, config.json, README.txt, "
                       "manifest.json and validation_report.json are always included. No force "
                       "field or charges are assigned.</small>")
        self.export_adv = w.VBox([section("Advanced", True), self.freeze])
        self.save_as = w.Text(description="Save as", placeholder="empty: a Buddy name, e.g. "
                                                              "Sweet_Pea_Buddy_Au_sphere",
                              style=STYLE, layout=wide(), continuous_update=False)
        self.save_as_info = w.HTML("<small>Name of the download (ZIP, certificate, session); "
                                   "letters, numbers, - and _ are kept.</small>")
        tab_export = w.VBox([section("Save as"), self.save_as, self.save_as_info,
                             section("File formats"), self.t_hoomd, self.t_lammps, self.t_generic,
                             self.t_pdb, self.t_mol2, self.residue_scheme, units, self.export_adv])

        self.tabs = w.Tab(children=[tab_material, tab_geometry, self.carbon.tab, tab_ligand,
                                    tab_coverage, tab_export])
        for i, title in enumerate(["1 Material", "2 Shape", "3 Defects & GO", "4 Ligand",
                                   "5 Coverage", "6 Export"]):  # the tab bar wraps when narrow
            self.tabs.set_title(i, title)

        # Actions, status, downloads
        self.example = w.Dropdown(options=[("Load an example...", "")] +
                                  [(e["label"], k) for k, e in EXAMPLES.items()],
                                  layout=w.Layout(width="100%", max_width="300px"))
        self.validate_btn = w.Button(description="Check inputs", icon="check-square")
        self.generate_btn = w.Button(description="Generate model", button_style="primary", icon="cogs",
                                     layout=w.Layout(width="auto"))
        self.reset_btn = w.Button(description="Reset", icon="undo", layout=w.Layout(width="auto"))
        self.validate_btn.layout.width = "auto"
        self.session_upload = w.FileUpload(accept=".py,.txt,.json", multiple=False,
                                           description="Load session", icon="folder-open",
                                           layout=w.Layout(width="auto"))
        self.session_upload.tooltip = ("Load a session.py from Save session (or a config.json) to "
                                       "restore every selection")
        self.popup = w.HTML()  # sad Buddy when a build fails (fades out by itself)
        self.progress = w.IntProgress(min=0, max=5, layout=w.Layout(width="200px"))
        self.progress_label = w.Label()
        self.status = w.HTML()
        self.downloads = w.HTML()

        # Right column: viewer (click to inspect, measure, edit), downloads, result, console
        self.viewer = w.HTML()
        self.bridge = w.Text(value="", continuous_update=True, layout=w.Layout(display="none"))
        self._bridge_class = f"bb-bridge-{id(self):x}"
        self.bridge.add_class(self._bridge_class)
        self.edit_info = w.HTML()
        self.undo_btn = w.Button(description="Undo edit", icon="undo", disabled=True,
                                 tooltip="Undo the last viewer or console edit")
        self.results = w.HTML()
        self.inspect_index = w.BoundedIntText(min=0, max=0, description="Inspect atom #",
                                              style={"description_width": "110px"},
                                              layout=w.Layout(width="220px"))
        self.inspect_info = w.HTML()
        self.console_log = w.HTML()
        self.console_log.add_class("bb-console-log")
        self.console_in = w.Textarea(placeholder='Python commands, e.g. material("Au"); shape("sphere", '
                                                 'diameter=4); build()    Enter runs, Shift+Enter '
                                                 'adds a line; paste a saved session here. help() '
                                                 'lists all',
                                     continuous_update=False, rows=2,
                                     layout=w.Layout(width="auto", flex="1 1 auto"))
        self.console_in.add_class("bb-console-input")
        self.console_run = w.Button(description="Run", icon="play",
                                    layout=w.Layout(width="auto", flex="0 0 auto"))
        # Enter in the console runs it (one document-level listener, installed once).
        keys_js = ("if(!window.__bbKeys){window.__bbKeys=1;document.addEventListener('keydown',"
                   "function(e){var t=e.target;if(e.key==='Enter'&&!e.shiftKey&&t&&t.tagName==='TEXTAREA'"
                   "&&t.closest('.bb-console-input, .bb-chat-input')){e.preventDefault();e.stopPropagation();"
                   "t.dispatchEvent(new Event('change',{bubbles:true}));"
                   "var r=t.closest('.bb-console-row, .bb-chat-row');"
                   "var b=r&&r.querySelector('button');if(b){setTimeout(function(){b.click();},40);}}},true);}")
        self.console_keys = w.HTML(
            "<img alt='' style='display:none' src='data:image/gif;base64,"
            "R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7' "
            f"onload=\"{keys_js}\" onerror=\"{keys_js}\">")
        self.console = Console(self)
        if self._want_assistant:  # Ask Buddy chat + console AI; optional, never blocks the app
            try:
                from .assistant.panel import AssistantPanel
                self.assistant = AssistantPanel(self)
            except Exception:
                self.assistant = None

        # Layout: inputs on the left, the viewer "stage" and its results on the right.
        actions = w.HBox([self.generate_btn, self.validate_btn, self.reset_btn, self.example,
                          self.session_upload],
                         layout=w.Layout(flex_flow="row wrap", align_items="center",
                                         margin="14px 0 4px", grid_gap="8px"))
        progress = w.HBox([self.progress, self.progress_label], layout=w.Layout(align_items="center"))
        chat = [self.assistant.box] if self.assistant is not None else []
        left = w.VBox(chat + [self.tabs, actions, progress, w.HTML("<div class='nb-h'>Status</div>"),
                       self.status],
                      layout=w.Layout(width="47%", min_width="min(470px, 100%)", flex="1 1 470px"))
        left.add_class("nb-left")
        right = w.VBox([self.viewer, self.bridge,
                        w.HBox([self.edit_info, self.undo_btn],
                               layout=w.Layout(align_items="center", justify_content="space-between",
                                               flex_flow="row wrap")),
                        w.HTML("<div class='nb-h'>Downloads</div>"), self.downloads,
                        w.HTML("<div class='nb-h'>Result</div>"), self.results,
                        row(self.inspect_index, self.inspect_info, margin="6px 0 0"),
                        w.HTML("<div class='nb-h'>Python console</div>"), self.console_log,
                        self._console_row(), self.console_keys],
                       layout=w.Layout(width="52%", min_width="min(400px, 100%)", flex="1 1 520px"))
        right.add_class("nb-right")
        self.theme = w.ToggleButtons(options=[("Buddy", "buddy"), ("Buddy 2", "buddy2"),
                                              ("Light", "light"), ("Dark", "dark")],
                                     value="buddy", tooltips=["Pastel yellow Buddy colours",
                                                              "Pastel mint-green Buddy colours",
                                                              "Light theme", "Dark theme"],
                                     style={"button_width": "auto"})
        self.mode = w.ToggleButtons(options=[("Basic", "basic"), ("Advanced", "advanced")], value="basic",
                                    tooltips=["Only the choices most builds need",
                                              "Every setting: seeds, tolerances, lattice constants, "
                                              "box edges, carbon details"],
                                    style={"button_width": "auto"})
        theme_box = w.HBox([w.HTML("<span class='nb-switch-label'>Theme</span>"), self.theme],
                           layout=w.Layout(align_items="center"))
        theme_box.add_class("nb-theme")
        mode_box = w.HBox([w.HTML("<span class='nb-switch-label'>Settings</span>"), self.mode],
                          layout=w.Layout(align_items="center"))
        mode_box.add_class("nb-theme")
        switches = w.VBox([mode_box, theme_box], layout=w.Layout(align_items="flex-end"))
        switches.add_class("nb-switches")
        logo = _logo_uri()
        logo_html = f"<img class='nb-logo' src='{logo}' alt='Buddy Builder logo'>" if logo else ""
        title = w.HTML(f"<div class='nb-brand'>{logo_html}<div><h1 class='nb-title'>{APP_NAME}</h1>"
                       "<p class='nb-sub'>Build ligand-coated nanoparticles, nanorods, slabs and "
                       "channels, or graphene, graphene oxide, nanotubes and fullerenes with "
                       "covalently bound ligands, then export structure and topology files for "
                       "LAMMPS, HOOMD-blue, AMBER and GROMACS (no force field).</p></div></div>")
        header = w.HBox([title, switches], layout=w.Layout(justify_content="space-between",
                                                           align_items="center", flex_flow="row wrap"))
        header.add_class("nb-head")
        footer = w.HTML(f"<small>{APP_NAME} (nanobuilder {__version__}), loaded from "
                        f"{html.escape(str(Path(__file__).resolve().parent))}</small>",
                        layout=w.Layout(margin="18px 0 0"))
        self.style = w.HTML(stylesheet())
        self.root = w.VBox([self.style, self.popup, header,
                            w.HBox([left, right], layout=w.Layout(flex_flow="row wrap",
                                                                  align_items="flex-start")),
                            footer])
        self.root.add_class("nb-app")
        self.root.add_class("nb-buddy")
        self._advanced_boxes = [self.lattice_adv, self.ligand_adv, self.coverage_adv,
                                self.export_adv]

    # ----------------------------------------------------------------- wiring

    def _config_widgets(self):
        return [self.source, self.preset, self.lat_a, self.lat_c, self.upload_format, self.shape,
                self.diameter, self.slab_len, self.slab_wid, self.thick_mode, self.layers, self.thickness,
                self.h, self.k, self.l, self.termination, self.vacuum, self.max_atoms, self.smiles,
                self.rod_len, self.rod_u, self.rod_v, self.rod_w, self.rod_caps,
                self.channel_type, self.channel_width, self.channel_count, self.box_shape,
                self.box_x, self.box_y, self.box_z, self.box_edge, self.mp_pick, self.save_as,
                self.anchor, self.template, self.distance, self.orientation, self.conformer_mode,
                self.conformers,
                self.lig_seed, self.cov_mode, self.count, self.density, self.region, self.distribution,
                self.min_sep, self.clash, self.max_attempts, self.cov_seed, self.t_hoomd, self.t_lammps,
                self.t_generic, self.t_pdb, self.t_mol2, self.residue_scheme, self.fidelity,
                self.freeze] + self.carbon.widgets()

    def _wire(self):
        for widget in self._config_widgets():
            widget.observe(self._on_config_change, names="value")
        self.source.observe(self._guard(self._on_material_change), names="value")
        self.carbon.kind.observe(self._guard(self._on_material_change), names="value")
        for widget in self.carbon.widgets():
            if widget is not self.carbon.kind:
                widget.observe(self._guard(self._on_carbon_change), names="value")
        self.preset.observe(self._guard(self._on_material_change), names="value")
        self.menu_bridge.observe(self._guard(self._on_menu), names="value")
        self.mat_search.observe(self._guard(self._on_find), names="value")
        self.mat_go.on_click(self._guard(self._on_find))
        self.mp_pick.observe(self._guard(self._on_mp_pick), names="value")
        self.upload_format.observe(self._guard(self._on_upload_format), names="value")
        self.upload_widget.observe(self._guard(self._on_upload), names="value")
        self.shape.observe(self._guard(self._on_shape_change), names="value")
        self.thick_mode.observe(self._guard(self._update_visibility), names="value")
        for widget in (self.h, self.k, self.l, self.layers):
            widget.observe(self._guard(self._refresh_terminations), names="value")
        self.smiles.observe(self._guard(self._on_smiles), names="value")
        self.library.observe(self._guard(self._on_library), names="value")
        self.anchor.observe(self._guard(self._on_anchor), names="value")
        self.template.observe(self._guard(self._on_template), names="value")
        self.cov_mode.observe(self._guard(self._update_visibility), names="value")
        self.preview_cell_btn.on_click(self._guard(self._preview_cell))
        self.preview_sub_btn.on_click(self._guard(self._preview_substrate))
        self.smiles_btn.on_click(self._guard(lambda _: self._on_smiles(None)))
        self.validate_btn.on_click(self._guard(self._validate))
        self.generate_btn.on_click(self._guard(self._generate))
        self.reset_btn.on_click(self._guard(lambda _: self.reset()))
        self.example.observe(self._guard(self._on_example), names="value")
        self.theme.observe(self._on_theme, names="value")
        self.mode.observe(self._guard(self._update_visibility), names="value")
        self.channel_type.observe(self._guard(self._on_shape_change), names="value")
        self.t_mol2.observe(self._guard(self._update_visibility), names="value")
        self.inspect_index.observe(self._guard(self._inspect), names="value")
        self.bridge.observe(self._guard(self._on_bridge), names="value")
        self.undo_btn.on_click(self._guard(self._undo))
        self.session_upload.observe(self._guard(self._on_session_upload), names="value")
        self.console_run.on_click(self._on_console)

    def _guard(self, fn):
        """Run a callback, showing any error in the status panel instead of losing it."""
        def wrapper(*args):
            if self._quiet:
                return
            try:
                fn(*args)
            except BuilderError as exc:
                self.status.value = _msg("error", str(exc))
            except Exception as exc:  # unexpected: show the traceback for bug reports
                self.status.value = _msg("error", f"Unexpected error: {exc}") + (
                    f"<details><summary>Details</summary><pre>{html.escape(traceback.format_exc())}"
                    "</pre></details>")
        return wrapper

    def show(self):
        display(self.root)

    # ------------------------------------------------------------ theme + viewer

    def _show_view(self, payload, height=460):
        """Render a viewer payload in the current theme (kept for re-theming)."""
        self._view_payload = (payload, height)
        self.viewer.value = iframe_html(payload, height, self.theme.value, self._bridge_class)

    def _show_assembly(self, view=None, selection=()):
        """Show the current model; editable in the viewer when it is a generated build."""
        a = self.result.assembly
        payload = assembly_payload(a, f"{a.n_atoms} atoms, {a.achieved_ligands} ligands")
        payload["editable"] = self.result.ligand is not None
        payload["view"] = view
        payload["sel"] = sorted(int(i) for i in selection)
        self._show_view(payload)

    def _on_theme(self, change=None):
        for name in ("nb-auto", "nb-light", "nb-dark", "nb-buddy", "nb-buddy2"):
            self.root.remove_class(name)
        self.root.add_class(f"nb-{self.theme.value}")
        if getattr(self, "_view_payload", None) is not None:
            payload, height = self._view_payload
            self.viewer.value = iframe_html(payload, height, self.theme.value, self._bridge_class)

    def _ipython_display_(self):
        self.show()

    # ------------------------------------------------------------ apply state

    def reset(self):
        """Restore the launch state: Au/FCC and the original ligand demo."""
        self.example.value = ""
        self.upload = None
        cfg = RunConfig()
        cfg.coverage.distribution = "fibonacci"  # new builds: ordered placement
        cfg.geometry.box_shape = "cube"  # new builds: the box decides the size (8 nm cube)
        cfg.geometry.box_lengths_nm = (8.0, 8.0, 8.0)
        self._apply_config(cfg)
        self.result = self.package = None
        self.history = []
        self._stale = False
        self.viewer.value = placeholder_html("Preview a unit cell or substrate, or click Generate model.")
        try:
            self._preview_cell(quiet=True)
        except Exception:  # the viewer is optional at start-up
            pass
        self.results.value = ""
        self.edit_info.value = ""
        self.undo_btn.disabled = True
        self.undo_btn.layout.display = "none"
        self.downloads.value = "<i>No model generated yet.</i>"
        self.status.value = _msg("info", "Ready. Defaults: Au/FCC, 3 nm sphere, 12 SCCO ligands "
                                         "(Au-thiolate) in Fibonacci order. Go through tabs 1 to 6, "
                                         "or just click Generate model.")
        self.progress.value = 0
        self.progress_label.value = ""

    def _apply_config(self, cfg: RunConfig, upload=None):
        self._quiet = True
        try:
            m, g, lg, c, e = cfg.material, cfg.geometry, cfg.ligand, cfg.coverage, cfg.export
            self.source.value = m.source
            self.preset.value = m.preset_id or default_preset_id()
            self.upload_format.value = m.file_format or "cif"
            if m.source == "materials_project":
                keys = [k for _, k in self.mp_pick.options]
                if m.mp_id not in keys:
                    raise BuilderError(f"{m.mp_id} is not in the preloaded Materials Project "
                                       "library; upload its CIF instead")
                self.mp_pick.value = m.mp_id
                self.upload = None
            else:
                self.upload = upload
            self.shape.value = g.shape
            self.diameter.value = g.diameter_nm
            self.slab_len.value, self.slab_wid.value = g.slab_length_nm, g.slab_width_nm
            self.thick_mode.value = g.thickness_mode
            self.layers.value, self.thickness.value = g.layers, g.thickness_nm
            self.h.value, self.k.value, self.l.value = (int(i) for i in g.miller)
            self.rod_len.value, self.rod_caps.value = g.rod_length_nm, g.rod_caps
            self.rod_u.value, self.rod_v.value, self.rod_w.value = (int(i) for i in g.rod_axis)
            self.vacuum.value, self.max_atoms.value = g.vacuum_A, g.max_atoms
            self.channel_type.value, self.channel_width.value = g.channel_type, g.channel_width_nm
            self.channel_count.value = g.channel_count
            user_box = g.box_lengths_nm is not None and any(float(x) > 0 for x in g.box_lengths_nm)
            self.box_shape.value = g.box_shape if (user_box or g.box_shape == "auto") else "auto"
            if g.box_shape != "auto" and not user_box:
                self.box_shape.value = "auto"  # an old "fit, then make it cubic" box: keep it automatic
            if user_box and g.box_shape == "auto":
                self.box_shape.value = "rectangular"
            if user_box:
                self.box_edge.value = float(max(g.box_lengths_nm))
                self.box_x.value, self.box_y.value, self.box_z.value = (float(x) for x in g.box_lengths_nm)
            self.smiles.value = lg.smiles
            self._refresh_library(carbon=m.source == "carbon")
            keys = [k for _, k in self.library.options]
            self.library.value = lg.library_id if lg.library_id in keys else ""
            self.conformer_mode.value = lg.conformer_mode
            self.orientation.value = lg.orientation
            self.conformers.value, self.lig_seed.value = lg.conformers, lg.seed
            self.cov_mode.value = c.mode
            self.count.value, self.density.value = c.target, c.density_nm2
            self.distribution.value = c.distribution
            self.min_sep.value, self.clash.value = c.min_separation_A, c.clash_alpha
            self.max_attempts.value, self.cov_seed.value = c.max_attempts, c.seed
            self.t_hoomd.value = "hoomd_gsd" in e.targets
            self.t_lammps.value = "lammps_data" in e.targets
            self.t_generic.value = "extxyz" in e.targets
            self.t_pdb.value = "pdb" in e.targets
            self.t_mol2.value = "mol2" in e.targets
            self.residue_scheme.value = e.residue_scheme
            self.fidelity.value = e.fidelity
            self.freeze.value = e.freeze_substrate
            self.save_as.value = e.save_as or ""
            self.carbon.apply(cfg.carbon)
        finally:
            self._quiet = False
        self._refresh_material(lattice=m.lattice_overrides_A)
        self._refresh_ligand(anchor=lg.anchor_atom_index)
        self._refresh_attachment(template=lg.binding_template, species=lg.surface_species,
                                 distance=lg.distance_A)
        self._refresh_regions(c.region)
        self._refresh_terminations(index=getattr(g, "termination_index", 0))
        self._update_visibility()

    # ------------------------------------------------------------ material

    def _on_material_change(self, _=None):
        self._refresh_material()
        switched = self._refresh_library(carbon=self._is_carbon())
        if switched:
            self._on_library()  # fills SMILES, anchor, template, distance
        else:
            self._refresh_attachment()
        self._refresh_regions()
        self._refresh_terminations()
        self._update_visibility()
        self._material_status(switched)
        if self.crystal is not None:
            try:  # no box or size yet: show the unit cell of the new material
                self._preview_cell(quiet=True)
            except Exception:
                pass

    def _refresh_library(self, carbon: bool) -> bool:
        """Offer only the ligands that suit the material; True when the ligand had to change.

        A library ligand that does not suit the new material is replaced by a
        sensible default (4-carboxyphenyl on carbon, 2-mercaptoethanol on a
        crystal); a typed SMILES is kept unless it only works on a crystal.
        """
        current = self.library.value
        options = [("Custom SMILES", "")] + library_choices(carbon)
        if [k for _, k in self.library.options] != [k for _, k in options]:
            self._quiet = True
            try:
                self.library.options = options
                self.library.value = current if current and suits(current, carbon) else ""
            finally:
                self._quiet = False
        if self._quiet:
            return False
        fallback = "carboxyphenyl" if carbon else "mercaptoethanol"
        template = self.template.value
        needs = (current and not suits(current, carbon)) or (
            carbon and not current and (template is None or TEMPLATES[template].substrate != "carbon"))
        if needs:
            self._quiet = True
            try:
                self.library.value = fallback
            finally:
                self._quiet = False
            return True
        return False

    def _material_status(self, switched=False):
        """Say what is loaded and what to do next, so the status never shows stale text."""
        if self._is_carbon():
            name = {"graphene": "Graphene / graphene oxide sheet", "nanotube": "Carbon nanotube",
                    "carbon_sphere": "Fullerene / carbon onion"}[self.carbon.kind.value]
            nxt = "set its size on tab 2 and defects or oxygen groups on tab 3"
        elif self.crystal is not None:
            name = self.crystal.structure_label
            nxt = "choose a shape and size on tab 2"
        else:
            return
        lines = [f"Material: {name}. Next: {nxt}, then Generate model."]
        if switched:
            lines.append(f"Ligand changed to {library_entry(self.library.value)['label']}, which "
                         "suits this material (tab 4 lists only suitable ligands).")
        self.status.value = _msg("info", lines)

    def _refresh_material(self, lattice=None):
        """Load the selected preset, library structure or uploaded file into ``self.crystal``."""
        self.crystal = None
        self.cell_info.value = ""
        if self.source.value == "carbon":
            try:
                self.crystal = material_crystal(RunConfig(material=self._material_config(),
                                                          carbon=self.carbon.config()))
            except Exception:  # e.g. an impossible (n, m) while typing; reported on Validate
                self.crystal = None
            self.cell_info.value = ("<small>Carbon frameworks are generated directly; the repeat "
                                    "unit is shown by Preview unit cell.</small>")
            return
        if self.source.value == "preset":
            entry = preset(self.preset.value)
            self._quiet = True
            try:
                self.lat_a.value = (lattice or {}).get("a", entry["lattice_A"]["a"])
                if "c" in entry["lattice_A"]:
                    self.lat_c.value = (lattice or {}).get("c", entry["lattice_A"]["c"])
            finally:
                self._quiet = False
            self.lat_c.layout.display = "" if "c" in entry["free_parameters"] else "none"
            self.crystal = build_crystal(self._material_config())
        elif self.source.value == "materials_project":
            if not self.mp_pick.value:
                return
            self.crystal = build_crystal(self._material_config())
            mid = self.mp_pick.value
            self.mp_summary.value = (
                f"<small><a href='https://materialsproject.org/materials/{html.escape(mid)}' "
                "target='_blank' rel='noopener'>Open on materialsproject.org</a> &middot; preloaded, "
                "works offline &middot; data CC BY 4.0, cite Jain et al., APL Mater. 1, 011002 "
                "(2013)</small>")
        elif self.upload is not None:
            try:
                self.crystal = build_crystal(self._material_config())
                self.upload_summary.value = _msg("ok", [
                    f"{self.upload[0]}: {self.crystal.atoms.get_chemical_formula()}",
                    f"SHA-256 {self.crystal.source_sha256[:16]}..."])
            except BuilderError as exc:
                self.upload_summary.value = _msg("error", str(exc))
        else:
            self.upload_summary.value = "<i>No file uploaded yet.</i>"
        if self.crystal is not None:
            self.cell_info.value = self._cell_table(self.crystal)

    @staticmethod
    def _cell_table(crystal) -> str:
        """Unit-cell facts for the Material tab (no box or size is needed yet)."""
        from ase.data import atomic_masses
        atoms = crystal.atoms
        a, b, c, al, be, ga = atoms.cell.cellpar()
        volume = atoms.get_volume()
        mass = float(sum(atomic_masses[z] for z in atoms.numbers))
        counts = {}
        for e in atoms.get_chemical_symbols():
            counts[e] = counts.get(e, 0) + 1
        system = ("cubic" if np.allclose([a, b], c, rtol=1e-3) and np.allclose([al, be, ga], 90, atol=0.05)
                  else "tetragonal" if np.isclose(a, b, rtol=1e-3) and np.allclose([al, be, ga], 90, atol=0.05)
                  else "hexagonal" if np.isclose(a, b, rtol=1e-3) and np.isclose(ga, 120, atol=0.1)
                  else "orthorhombic" if np.allclose([al, be, ga], 90, atol=0.05) else "lower symmetry")
        rows = [("Material", crystal.structure_label),
                ("Space group", crystal.space_group or "-"),
                ("Lattice", f"{system}"),
                ("a, b, c", f"{a:.4f}, {b:.4f}, {c:.4f} \u00c5"),
                ("\u03b1, \u03b2, \u03b3", f"{al:.2f}, {be:.2f}, {ga:.2f}\u00b0"),
                ("Atoms per cell", f"{len(atoms)} (" + ", ".join(f"{k} {v}" for k, v in counts.items())
                 + ")"),
                ("Cell volume", f"{volume:.2f} \u00c5\u00b3"),
                ("Density", f"{mass * 1.66053906660 / volume:.3f} g/cm\u00b3" if volume > 0 else "-"),
                ("Source", crystal.lattice_source or "-")]
        body = "".join(f"<tr><td>{html.escape(k)}</td><td>{html.escape(str(v))}</td></tr>"
                       for k, v in rows)
        return f"<table class='nb-table nb-cell'>{body}</table>"

    def _material_config(self) -> MaterialConfig:
        if self.source.value == "carbon":
            return MaterialConfig(source="carbon", preset_id=None)
        if self.source.value == "preset":
            entry = preset(self.preset.value)
            overrides = {}
            if abs(self.lat_a.value - entry["lattice_A"]["a"]) > 1e-9:
                overrides["a"] = float(self.lat_a.value)
            if "c" in entry["free_parameters"] and abs(self.lat_c.value - entry["lattice_A"]["c"]) > 1e-9:
                overrides["c"] = float(self.lat_c.value)
            return MaterialConfig(source="preset", preset_id=self.preset.value,
                                  lattice_overrides_A=overrides)
        if self.source.value == "materials_project":
            mid = self.mp_pick.value
            return MaterialConfig(source="materials_project", preset_id=None, mp_id=mid,
                                  file_format="cif", filename=f"{mid}.cif" if mid else None,
                                  uploaded_bytes=mp_library.cif_bytes(mid) if mid else None)
        name, data = self.upload if self.upload else (None, None)
        return MaterialConfig(source="uploaded_file", preset_id=None, file_format=self.upload_format.value,
                              filename=name, uploaded_bytes=data)

    # ------------------------------------------------------------ material picker

    def _menu_html(self) -> str:
        """Browse menu: hover a group to see its materials, click one to load it."""
        groups = []
        for group, names in catalog.menu():
            items = "".join(
                f"<li><a data-name=\"{html.escape(n, quote=True)}\" onclick=\"var el=document."
                f"querySelector('.{self._menu_class} input');el.value=JSON.stringify({{name:this."
                f"dataset.name,t:Date.now()}});el.dispatchEvent(new Event('input',{{bubbles:true}}));"
                f"this.closest('.bb-menu').classList.add('bb-closed');setTimeout(()=>this.closest("
                f"'.bb-menu').classList.remove('bb-closed'),400);\">{html.escape(n)}</a></li>"
                for n in names)
            groups.append(f"<li class='bb-group'><span>{html.escape(group)}</span>"
                          f"<ul class='bb-items'>{items}</ul></li>")
        return (f"<div class='bb-menu'><button type='button' class='bb-menu-btn'>Browse materials "
                f"&#9662;</button><ul class='bb-groups'>{''.join(groups)}</ul></div>")

    def _on_menu(self, _=None):
        import json
        raw = (self.menu_bridge.value or "").strip()
        if not raw:
            return
        try:
            name = json.loads(raw).get("name", "")
        except ValueError:
            return
        key = catalog.key_for_text(name)
        if key:
            self._quiet = True
            try:
                self.mat_search.value = name
            finally:
                self._quiet = False
            self._apply_catalog(catalog.CATALOG[key])

    def _set_source(self, source, preset_id=None, mp_id=None):
        self._quiet = True
        try:
            self.source.value = source
            if preset_id:
                self.preset.value = preset_id
            if mp_id:
                self.mp_pick.value = mp_id
        finally:
            self._quiet = False

    def _on_find(self, _=None):
        text = (self.mat_search.value or "").strip()
        if not text:
            return
        key = catalog.key_for_text(text, fuzzy=False)
        if key is not None:
            self._apply_catalog(catalog.CATALOG[key])
            return
        preset_id = catalog.preset_for_text(text)
        if preset_id:
            self._set_source("preset", preset_id)
            self._on_material_change()
            self._on_config_change()
            self.mat_hint.value = (f"<small>Using the built-in preset {html.escape(preset(preset_id)['label'])}."
                                   "</small>")
            return
        keys = catalog.library_for_text(text)
        if keys:
            self._apply_catalog(catalog.CATALOG[keys[0]])
            if len(keys) > 1:
                others = ", ".join(catalog.KEY_TO_NAME[k] for k in keys[1:])
                self.mat_hint.value = (f"<small>{html.escape(catalog.KEY_TO_NAME[keys[0]])} loaded; "
                                       f"also in the library: {html.escape(others)} (Structure "
                                       "menu below).</small>")
            return
        key = catalog.key_for_text(text)  # a word of a name: "anatase", "nitinol"
        if key is not None:
            self._apply_catalog(catalog.CATALOG[key])
            return
        if text in mp_library.by_id():  # an mp id from the library
            self._set_source("materials_project", mp_id=text)
            self._on_material_change()
            self._on_config_change()
            return
        raise BuilderError(f"{text} is not in the built-in presets or the preloaded Materials "
                           "Project library. Pick one from Browse materials, or choose Upload "
                           "file and load its CIF")

    def _on_mp_pick(self, _=None):
        if self.mp_pick.value and self.source.value == "materials_project":
            self._on_material_change()
            self._on_config_change()

    def _apply_catalog(self, entry):
        if entry["kind"] == "carbon":
            self._quiet = True
            try:
                self.source.value = "carbon"
                self.carbon.apply(CarbonConfig(**entry["carbon"]))
            finally:
                self._quiet = False
            self._on_material_change()
            self.mat_hint.value = ("<small>Carbon nanomaterial loaded: set its size on the Shape "
                                   "tab and its defects and oxygen groups on the Defects &amp; GO "
                                   "tab.</small>")
            self._on_config_change()
            return
        if entry["kind"] == "preset":
            self._set_source("preset", entry["preset_id"])
        else:
            self._set_source("materials_project", mp_id=entry["mp_id"])
        self._on_material_change()
        geometry = entry.get("geometry")
        if geometry:
            self._quiet = True
            try:
                self.shape.value = geometry.get("shape", self.shape.value)
                self.diameter.value = geometry.get("diameter_nm", self.diameter.value)
                self.rod_len.value = geometry.get("rod_length_nm", self.rod_len.value)
                self.rod_caps.value = geometry.get("rod_caps", self.rod_caps.value)
                u, v, wv = geometry.get("rod_axis", (self.rod_u.value, self.rod_v.value, self.rod_w.value))
                self.rod_u.value, self.rod_v.value, self.rod_w.value = u, v, wv
            finally:
                self._quiet = False
            self._on_shape_change()
            self.mat_hint.value = ("<small>Geometry set to a nanorod; adjust its size on the "
                                   "Shape tab.</small>")
        self._on_config_change()

    def _on_upload_format(self, _=None):
        self.upload_widget.accept = ".cif" if self.upload_format.value == "cif" else ""
        if self.upload:
            self._on_material_change()

    def _on_upload(self, _=None):
        file = _uploaded_file(self.upload_widget)
        if file is None:
            return
        self.upload = file
        self._on_material_change()
        self._on_config_change()

    # ------------------------------------------------------------ geometry

    def _on_shape_change(self, _=None):
        self._refresh_regions()
        self._refresh_terminations()
        self._update_visibility()

    def _is_carbon(self) -> bool:
        return self.source.value == "carbon"

    def _shape(self) -> str:
        """Substrate shape: the Shape dropdown, or sheet/nanotube/carbon_sphere for carbon."""
        return self.carbon.shape() if self._is_carbon() else self.shape.value

    def _on_carbon_change(self, _=None):
        if self._is_carbon():
            self._refresh_material()
            self._refresh_regions()
            self._update_visibility()

    def _refresh_regions(self, value=None):
        shape = self._shape()
        if shape == "sheet":
            options = [("Top face", "top"), ("Bottom face", "bottom"), ("Both faces", "both"),
                       ("Edges (sheet and hole edges)", "edges"),
                       ("Entire surface (faces + edges)", "entire_surface")]
            if self.carbon.layers.value > 1:
                options.append(("Between layers (inner faces of the stack)", "inner"))
        elif shape == "nanotube":
            options = [("Entire surface", "entire_surface"), ("Outer wall (side)", "side"),
                       ("Open rims (tips)", "tips"), ("Inner wall (inside the channel)", "inner")]
        elif shape == "channel":
            walls = ("Pore walls (inside the pores)" if self.channel_type.value == "pore"
                     else "Channel walls (inside)")
            options = [(walls, "inner"), ("Top outer face", "top"),
                       ("Bottom outer face", "bottom"), ("Both outer faces", "both"),
                       ("All surfaces", "entire_surface")]
        elif shape in ("sphere", "carbon_sphere"):
            options = [("Entire surface", "entire_surface")]
        elif shape == "rod":
            options = [("Entire surface", "entire_surface"), ("Side only", "side"),
                       ("Tips only", "tips")]
        else:
            options = [("Top face", "top"), ("Bottom face", "bottom"), ("Both faces", "both")]
        current = value or self.region.value
        if value is None and shape != getattr(self, "_region_shape", shape) and shape == "channel":
            current = "inner"  # the point of a channel: ligands on its walls
        self._region_shape = shape
        self._quiet = True
        try:
            self.region.options = options
            keys = [o[1] for o in options]
            self.region.value = current if current in keys else keys[0]
        finally:
            self._quiet = False

    def _refresh_terminations(self, _=None, index=None):
        if self._is_carbon() or self.shape.value not in ("slab", "channel") or self.crystal is None:
            self.termination.options, self.term_info.value = [], ""
            return
        miller = (self.h.value, self.k.value, self.l.value)
        if miller == (0, 0, 0):
            self.termination.options = []
            self.term_info.value = _msg("error", "Miller indices (0, 0, 0) are not a plane")
            return
        try:
            terms = list_terminations(self.crystal, miller, max(1, int(self.layers.value)))
        except BuilderError as exc:
            self.termination.options = []
            self.term_info.value = _msg("error", str(exc))
            return
        keep = index if index is not None else self.termination.value
        self._quiet = True
        try:
            self.termination.options = [(t.label, t.index) for t in terms]
            self.termination.value = keep if keep is not None and keep < len(terms) else 0
        finally:
            self._quiet = False
        note = ("Only one distinct termination for this plane." if len(terms) == 1 else
                f"{len(terms)} distinct terminations: choose one explicitly; they are not assumed "
                "equivalent, neutral or stable.")
        hcp = (" HCP: enter 3-index (hkl); for 4-index (hkil) drop i = -(h+k)."
               if self.crystal.prototype == "hcp" or len(set(self.crystal.atoms.cell.lengths().round(3))) > 1
               else "")
        self.term_info.value = (f"<small>{note} Miller indices refer to the "
                                f"{'conventional' if self.source.value == 'preset' else 'uploaded (normalized)'}"
                                f" cell.{hcp}</small>")

    # ------------------------------------------------------------ ligand

    def _on_smiles(self, _=None):
        # Editing the SMILES by hand leaves the library entry behind.
        if self.library.value and self.smiles.value.strip() != library_entry(self.library.value)["smiles"]:
            self._quiet = True
            try:
                self.library.value = ""
            finally:
                self._quiet = False
            self.library_info.value = ""
        self._refresh_ligand()
        self._refresh_attachment()

    def _on_library(self, _=None):
        """Fill SMILES, anchor, attachment mode, distance and chain mode from the library."""
        key = self.library.value
        if not key:
            self.library_info.value = ""
            return
        e = library_entry(key)
        self._quiet = True
        try:
            self.smiles.value = e["smiles"]
            self.conformer_mode.value = e["conformer_mode"]
        finally:
            self._quiet = False
        self._refresh_ligand(anchor=e["anchor_atom_index"])
        template = e["binding_template"]
        if self._is_carbon() and e.get("carbon_template"):
            template = e["carbon_template"]
        wanted = template
        if template == "Au_thiolate" and self.crystal is not None and \
                "Au_thiolate" not in [k for _, k in attachment_options(self.crystal)["templates"]]:
            template = "chemisorption"  # thiol on another metal: M-S bond, S-H removed
        self._refresh_attachment(template=template, distance=e["distance_A"]
                                 if template == wanted else None)
        note = e["note"]
        if wanted == "Au_thiolate" and self.template.value == "chemisorption":
            note += (" Au-thiolate is reviewed only for the Au preset, so the thiol binds by "
                     "chemisorption here (metal-S bond, S-H hydrogen removed).")
        elif template != self.template.value:
            note += (f" {TEMPLATES[template].label} is not available for this substrate; "
                     f"{TEMPLATES[self.template.value].label} is used instead.")
        self.library_info.value = f"<small>{html.escape(note)}</small>"
        self._on_config_change()

    def _refresh_ligand(self, anchor=None):
        smiles = self.smiles.value.strip()
        try:
            choices = anchor_choices(smiles)
            info = ligand_summary(smiles)
        except BuilderError as exc:
            self._quiet = True
            try:
                self.anchor.options = []
            finally:
                self._quiet = False
            self.smiles_info.value = _msg("error", str(exc))
            self.smiles_svg.value = ""
            return
        keep = anchor if anchor is not None else self.anchor.value
        valid = [c[1] for c in choices]
        self._quiet = True
        try:
            self.anchor.options = choices
            self.anchor.value = keep if keep in valid else default_anchor(smiles)
        finally:
            self._quiet = False
        thiols = ", ".join(str(i) for i in info["thiol_sulfurs"]) or "none"
        self.smiles_info.value = _msg("ok", f"{info['formula']}: {info['heavy_atoms']} heavy atoms, "
                                            f"{info['atoms_with_h']} with H; formal charge "
                                            f"{info['formal_charge']:+d}; thiol S atoms: {thiols}")
        self.smiles_svg.value = svg_2d(smiles, highlight=self.anchor.value)

    def _surface_element(self) -> str:
        ticked = [k for k, cb in self.species_checks.items() if cb.value]
        return ticked[0] if ticked else (next(iter(self.species_checks), "") if self.species_checks else "")

    def _anchor_element(self) -> str:
        try:
            from .ligands import parse_smiles
            return parse_smiles(self.smiles.value).GetAtomWithIdx(int(self.anchor.value)).GetSymbol()
        except Exception:  # no valid SMILES/anchor yet
            return ""

    def _on_anchor(self, _=None):
        if self.anchor.value is not None and self.smiles.value.strip():
            try:
                self.smiles_svg.value = svg_2d(self.smiles.value, highlight=self.anchor.value)
            except BuilderError:
                pass
        key = self.template.value
        if key and (TEMPLATES[key].substrate == "carbon" or key == "chemisorption"):
            self._quiet = True
            try:
                self.distance.value = template_distance(key, self._anchor_element(),
                                                        self._surface_element())
            finally:
                self._quiet = False

    def _refresh_attachment(self, template=None, species=None, distance=None):
        if self.crystal is None:
            options = {"templates": [(TEMPLATES["geometric"].label, "geometric")],
                       "reason": "Load a valid crystal first", "species": []}
        else:
            options = attachment_options(self.crystal)
        keep = template or self.template.value
        keys = [k for _, k in options["templates"]]
        self._quiet = True
        try:
            self.template.options = options["templates"]
            self.template.value = keep if keep in keys else keys[0]
        finally:
            self._quiet = False
        self._rebuild_species(options["species"], species)
        self._on_template(distance=distance)

    def _rebuild_species(self, elements, selected=None):
        """One checkbox per crystal element, all ticked by default.

        A previous choice is kept only while the element set is unchanged (the
        same crystal re-validated); a new crystal starts with every species on.
        """
        previous = {k for k, cb in self.species_checks.items() if cb.value}
        same_crystal = set(self.species_checks) == set(elements)
        if selected:
            selected = set(selected)
        elif same_crystal and self.species_checks:
            selected = previous
        else:
            selected = set(elements)
        self.species_checks = {}
        for element in elements:
            default = element in selected or len(elements) == 1
            cb = w.Checkbox(value=default, description=element, indent=False, layout=w.Layout(width="70px"))
            cb.observe(self._on_config_change, names="value")
            cb.observe(self._update_species_warning, names="value")
            self.species_checks[element] = cb
        self.species_box.children = list(self.species_checks.values())
        self._update_species_warning()

    def _update_species_warning(self, _=None):
        if self.species_checks and not any(cb.value for cb in self.species_checks.values()):
            self.species_warning.value = _msg(
                "warning", "No surface anchor species selected: ligands have nowhere to attach. "
                           "Tick at least one element.")
        else:
            self.species_warning.value = ""

    def _on_template(self, _=None, distance=None):
        key = self.template.value or "geometric"
        template = TEMPLATES[key]
        reason = attachment_options(self.crystal)["reason"] if self.crystal is not None else ""
        lines = [template.description]
        if reason:
            lines.append(reason)
        self.template_info.value = f"<small>{'<br>'.join(html.escape(x) for x in lines)}</small>"
        for element, cb in self.species_checks.items():
            if template.surface_species:
                cb.value = element in template.surface_species
                cb.disabled = True
            else:
                cb.disabled = len(self.species_checks) == 1
        self._quiet = True
        try:
            self.distance.value = distance if distance else template_distance(
                key, self._anchor_element(), self._surface_element())
        finally:
            self._quiet = False

    # ------------------------------------------------------------ visibility

    def _update_visibility(self, _=None):
        def show(widget, visible):
            widget.layout.display = "" if visible else "none"

        carbon = self._is_carbon()
        advanced = self.mode.value == "advanced"
        for box in self._advanced_boxes:
            show(box, advanced)
        show(self.preset_box, self.source.value == "preset")
        show(self.upload_box, self.source.value == "uploaded_file")
        show(self.mp_box, self.source.value == "materials_project")
        self.carbon.user_box = self.box_shape.value != "auto"
        self.carbon.update(carbon, advanced)
        shape = self.shape.value
        sphere, rod = shape == "sphere" and not carbon, shape == "rod" and not carbon
        channel = shape == "channel" and not carbon
        show(self.shape, not carbon)
        show(self.sphere_box, sphere or rod)
        show(self.rod_box, rod)
        show(self.slab_box, shape in ("slab", "channel") and not carbon)
        show(self.channel_box, channel)
        show(self.residue_scheme, self.t_mol2.value and carbon)
        box = self.box_shape.value
        show(self.box_xyz, box == "rectangular")
        show(self.box_edge, box == "cube")
        # With a user box, the box sets the in-plane size of periodic structures.
        show(self.slab_len, box == "auto")
        show(self.slab_wid, box == "auto")
        if channel:
            slit = self.channel_type.value == "slit"
            self.channel_width.description = "Gap between walls (nm)" if slit else "Pore diameter (nm)"
            self.channel_count.description = "Stacked channels" if slit else "Number of pores"
            self.channel_info.value = (
                "<small>" + ("Slit: the slab is cut into walls with empty gaps between them; "
                             "several channels stack the walls (a multi-layer channel). Periodic "
                             "in x and y; the stack must fit the box height." if slit else
                             "Pores: cylindrical holes straight through the membrane (along z), "
                             "spread evenly over the surface, each at least 3 atomic spacings "
                             "from the others and from the box edges.")
                + " Ligands go on the channel walls with region Channel walls (tab 5).</small>")
        periodic = (shape in ("slab", "channel") and not carbon) or (
            carbon and self.carbon.periodic.value and self.carbon.kind.value != "carbon_sphere")
        box_notes = {
            "auto": "The box is fitted to the model with 3 nm of vacuum on every side.",
            "rectangular": ("The structure fills x and y (a whole number of crystal repeats, so "
                            "the actual lengths can differ slightly; the build reports them); z "
                            "is the box height." if periodic else
                            "The model is centred in this box; a model that does not fit is "
                            "reported."),
            "cube": ("A square cell in x and y (whole crystal repeats) with the height matched "
                     "to it." if periodic else
                     "The model is centred in a cube with this edge; a model that does not fit is "
                     "reported."),
        }
        self.box_info.value = f"<small>{box_notes[box]}</small>"
        self.diameter.description = "Rod diameter (nm)" if rod else "Particle diameter (nm)"
        show(self.layers, self.thick_mode.value == "layers")
        show(self.thickness, self.thick_mode.value == "thickness")
        show(self.count, self.cov_mode.value == "count")
        show(self.density, self.cov_mode.value == "density")
        self.boundary.value = ("<small>Boundary: finite (isolated particle in the box)</small>"
                               if shape not in ("slab", "channel") else
                               "<small>Boundary: periodic in x and y, vacuum along z</small>")
        if carbon:
            self.boundary.value = ""
            region = self.region.value
            area = self.carbon.nominal_area_nm2(region)
            if region in ("edges", "tips"):
                text = ("Edge regions are counted as a ribbon (number of edge carbons x 0.23 nm x "
                        "0.335 nm); use an exact ligand count there.")
            elif area is None:
                text = "Area comes from the built substrate; it is shown after Preview or Generate."
            else:
                text = f"Nominal area of the selected region: {area:.2f} nm^2" + (
                    f"; {self.density.value:.3f}/nm^2 -> about {round(self.density.value * area)} "
                    "ligands" if self.cov_mode.value == "density" else "")
            self.area_info.value = f"<small>{html.escape(text)}</small>"
        elif sphere:
            area = 3.14159265 * self.diameter.value ** 2
            self.area_info.value = (f"<small>Geometric area (nominal sphere): {area:.3f} nm^2"
                                    + (f"; {self.density.value:.3f}/nm^2 -> about "
                                       f"{round(self.density.value * area)} ligands"
                                       if self.cov_mode.value == "density" else "") + "</small>")
        elif rod:
            r, length = self.diameter.value / 2, self.rod_len.value
            hemi = self.rod_caps.value == "hemispherical"
            cyl = max(length - (2 * r if hemi else 0.0), 0.0)
            side, tips = 2 * 3.14159265 * r * cyl, (4 if hemi else 2) * 3.14159265 * r ** 2
            area = {"side": side, "tips": tips}.get(self.region.value, side + tips)
            self.area_info.value = (f"<small>Nominal rod area: side {side:.2f} nm^2, tips {tips:.2f} "
                                    f"nm^2; selected region {area:.2f} nm^2"
                                    + (f"; {self.density.value:.3f}/nm^2 -> about "
                                       f"{round(self.density.value * area)} ligands"
                                       if self.cov_mode.value == "density" else "") + "</small>")
        else:
            what = "Channel wall" if channel else "Slab"
            self.area_info.value = (f"<small>{what} area comes from the actual periodic cell; it is "
                                    "shown after Preview bare substrate or Generate.</small>")

    # ------------------------------------------------------------ actions

    def config(self) -> RunConfig:
        """The RunConfig represented by the current widget values."""
        species = [k for k, cb in self.species_checks.items() if cb.value]
        if self.species_checks and not species:
            species_cfg = []  # explicitly none: the pipeline reports it
        else:
            species_cfg = species or None
        targets = [t for t, cb in (("hoomd_gsd", self.t_hoomd), ("lammps_data", self.t_lammps),
                                   ("extxyz", self.t_generic), ("pdb", self.t_pdb),
                                   ("mol2", self.t_mol2)) if cb.value]
        return RunConfig(
            material=self._material_config(),
            geometry=GeometryConfig(
                shape=self.shape.value, diameter_nm=self.diameter.value,
                slab_length_nm=self.slab_len.value, slab_width_nm=self.slab_wid.value,
                thickness_mode=self.thick_mode.value, layers=self.layers.value,
                thickness_nm=self.thickness.value, miller=(self.h.value, self.k.value, self.l.value),
                termination_index=self.termination.value or 0, vacuum_A=self.vacuum.value,
                max_atoms=self.max_atoms.value, rod_length_nm=self.rod_len.value,
                rod_axis=(self.rod_u.value, self.rod_v.value, self.rod_w.value),
                rod_caps=self.rod_caps.value, channel_type=self.channel_type.value,
                channel_width_nm=self.channel_width.value, channel_count=self.channel_count.value,
                box_shape=self.box_shape.value,
                box_lengths_nm=self._box_lengths()),
            ligand=LigandConfig(
                smiles=self.smiles.value.strip(), anchor_atom_index=self.anchor.value or 0,
                binding_template=self.template.value or "geometric", distance_A=self.distance.value,
                surface_species=species_cfg, conformers=self.conformers.value,
                orientation=self.orientation.value, seed=self.lig_seed.value,
                conformer_mode=self.conformer_mode.value, library_id=self.library.value or None),
            coverage=CoverageConfig(
                mode=self.cov_mode.value, target=self.count.value, density_nm2=self.density.value,
                region=self.region.value, distribution=self.distribution.value,
                min_separation_A=self.min_sep.value, clash_alpha=self.clash.value,
                max_attempts=self.max_attempts.value, seed=self.cov_seed.value),
            export=ExportConfig(targets=targets, fidelity=self.fidelity.value,
                                freeze_substrate=self.freeze.value,
                                residue_scheme=self.residue_scheme.value,
                                save_as=self.save_as.value.strip()),
            carbon=self.carbon.config(),
        )

    def _box_lengths(self):
        if self.box_shape.value == "rectangular":
            return (float(self.box_x.value), float(self.box_y.value), float(self.box_z.value))
        if self.box_shape.value == "cube":
            return (float(self.box_edge.value),) * 3
        return None

    def _on_config_change(self, _=None):
        if self._quiet:
            return
        if self.result is not None and not self._stale:
            self._stale = True
            self.downloads.value = _msg("warning", "Inputs changed since the last build; click "
                                                   "Generate model to update the downloads.")
        self._update_visibility()

    def _validate(self, _=None):
        outcome = validate_inputs(self.config())
        info = outcome["info"]
        lines = [f"{k}: {v}" for k, v in info.items()]
        if outcome["ok"]:
            self.status.value = _msg("ok", ["Inputs are valid."] + lines) + (
                _msg("warning", outcome["warnings"]) if outcome["warnings"] else "")
        else:
            self.status.value = _msg("error", outcome["errors"])

    def _preview_cell(self, _=None, quiet=False):
        if self.crystal is None:
            raise BuilderError("Load a valid crystal first (choose a preset or upload a file)")
        if self._is_carbon():
            crystal = material_crystal(self.config())
            self._show_view(unit_cell_payload(crystal, f"Repeat unit: {crystal.prototype}"))
            if not quiet:
                self.status.value = _msg("info", "Showing the repeat unit of the carbon framework; "
                                                 "use Preview bare substrate (Shape tab) for the full "
                                                 "structure with its defects and oxygen groups.")
            return
        crystal = build_crystal(self._material_config())  # picks up edited lattice constants
        self._show_view(unit_cell_payload(crystal, f"Unit cell: {crystal.structure_label}"))
        if not quiet:
            self.status.value = _msg("info", "Showing one unit cell; atoms on cell faces, edges and "
                                             "corners are repeated for display only.")

    def _preview_substrate(self, _=None):
        cfg = self.config()
        cfg.geometry.validate()
        if self._is_carbon():
            _, substrate = build_substrate(cfg)
            info = substrate.carbon_info
            comp = info["composition"]
            self._show_view(substrate_payload(
                substrate, f"{info['summary']}: {len(substrate.atoms)} atoms"))
            lines = [f"{info['summary']}: {len(substrate.atoms)} atoms (C {comp['C']}, O {comp['O']}, "
                     f"H {comp['H']})",
                     f"Oxygen content {comp['oxygen_percent_O_per_C']} % (C:O "
                     f"{comp['C_to_O_ratio'] or '-'}), modified carbons "
                     f"{comp['modified_carbon_percent']} %, vacancies {comp['vacancy_percent_actual']} %"]
            if info.get("groups"):
                lines.append("Groups: " + ", ".join(f"{k.replace('_', ' ')} {v}"
                                                    for k, v in info["groups"].items()))
            faces = ", ".join(f"{name} {f.area_nm2:.3f} nm^2" for name, f in substrate.faces.items())
            self.area_info.value = f"<small>Region areas: {html.escape(faces)}</small>"
            self.status.value = _msg("info", lines) + _msg("warning", substrate.warnings)
            return
        crystal = build_crystal(cfg.material)
        substrate = build_geometry(crystal, cfg.geometry)
        self._show_view(substrate_payload(
            substrate, f"Bare substrate: {len(substrate.atoms)} atoms"))
        faces = ", ".join(f"{name} {f.area_nm2:.3f} nm^2" for name, f in substrate.faces.items())
        if substrate.shape in ("slab", "rod", "channel"):
            self.area_info.value = f"<small>Face areas: {html.escape(faces)}</small>"
        self.status.value = _msg("info", [f"Bare substrate: {len(substrate.atoms)} atoms, composition "
                                          f"{substrate.composition}; faces: {faces}"]) + \
            _msg("warning", substrate.warnings)

    def _set_progress(self, stage, fraction):
        self.progress.value = int(round(fraction * self.progress.max))
        self.progress_label.value = "done" if stage == "done" else f"{stage}..."

    def _generate(self, _=None):
        self.last_failure = None
        cfg = self.config()
        self.status.value = _msg("info", "Building...")
        self.downloads.value = "<i>Building...</i>"
        self.generate_btn.disabled = True
        try:
            result = self._take_prebuilt(cfg) or run_pipeline(cfg, progress=self._set_progress)
        except BuilderError as exc:
            self.progress_label.value = "failed"
            self._fail(str(exc))
            return
        except Exception as exc:
            self.progress_label.value = "failed"
            self._sad("Oh no! Buddy could not build that.", "Something unexpected went wrong.",
                      "See Details under Status; Reset brings back a build that works.")
            raise exc
        finally:
            self.generate_btn.disabled = False
        self.result, self._stale = result, False
        self.history = []
        self.undo_btn.disabled = True
        self.undo_btn.layout.display = "" if result.ligand is not None else "none"
        a, report = result.assembly, result.report
        self._show_assembly()
        self._show_results(result)
        self.inspect_index.max = max(a.n_atoms - 1, 0)
        self._inspect()
        self.edit_info.value = ("<small>Click atoms in the viewer to inspect or measure; choose "
                                "<b>Remove ligand</b> or <b>Add ligand</b> in its Click menu to edit "
                                "the model.</small>" if result.ligand is not None else "")
        self.console.echo_build(cfg, result)
        if not report.ok:
            self._fail("; ".join(report.errors), validation=True)
            self.downloads.value = _msg("error", "No downloads: fix the validation errors above.")
            return
        previous, self.package = self.package, export_package(result)
        if previous is not None:  # each build exports into its own temporary folder
            shutil.rmtree(previous.zip_path.parent, ignore_errors=True)
        status_kind = "ok" if report.status == "ok" else "warning"
        notes = list(dict.fromkeys(report.warnings + a.warnings))
        # Shortfalls change what the user gets, so they stay visible; the standing
        # scientific caveats go in a collapsed list.
        urgent = [n for n in notes if n.startswith(("Only ", "Requested "))]
        rest = [n for n in notes if n not in urgent]
        details = (f"<details><summary>{len(rest)} notes on assumptions and limitations</summary>"
                   f"{_msg('info', rest)}</details>" if rest else "")
        self.status.value = _msg(status_kind, f"Build {report.status}: {a.achieved_ligands} of "
                                              f"{a.target_ligands} ligands placed in "
                                              f"{result.timings['total_s']:.1f} s.") + \
            _msg("warning", urgent) + details
        self._show_downloads()

    def _take_prebuilt(self, cfg):
        """A build Ask Buddy already made (on its hidden copy) for exactly these settings.

        Builds are seeded, so the same settings give the same model; reusing it saves building
        twice. Used once, and only when every setting matches.
        """
        pre, self._prebuilt = getattr(self, "_prebuilt", None), None
        if pre is None:
            return None
        key, result = pre
        try:
            if cfg.to_dict() != key:
                return None
        except Exception:
            return None
        self._set_progress("done", 1.0)
        return result

    # ------------------------------------------------------------ hand edits

    def apply_edit(self, action: str, value, view=None):
        """Remove ligands (value: list of ligand numbers) or add one (value: atom index).

        Used by viewer clicks and the console. The model is rebuilt, re-validated
        and re-exported; README.txt lists every edit.
        """
        if self.result is None:
            raise BuilderError("Generate a model first; edits change a generated model")
        if self._stale:
            raise BuilderError("Inputs changed since the last build: Generate again before editing")
        before = self.result
        if action == "remove":
            result = remove_ligands(before, value)
            what = f"Removed ligand{'s' if len(value) > 1 else ''} {', '.join(str(v) for v in value)}"
        elif action == "add":
            result = add_ligand_at(before, value)
            what = f"Added ligand {len(result.assembly.placements)} on atom {value}"
        else:
            raise BuilderError(f"Unknown edit {action!r}")
        self.history.append(before)
        self._set_result(result, view)
        self.edit_info.value = _msg("ok", f"{what}. {result.assembly.achieved_ligands} ligands now; "
                                          "the edit is flagged in README.txt.")
        return result

    def _set_result(self, result, view=None):
        self.result = result
        self.undo_btn.disabled = not self.history
        self.undo_btn.layout.display = ""
        a = result.assembly
        self._show_assembly(view)
        self._show_results(result)
        self.inspect_index.max = max(a.n_atoms - 1, 0)
        self._inspect()
        if not result.report.ok:
            self.status.value = _msg("error", ["Validation failed after the edit; downloads blocked."]
                                     + result.report.errors)
            self.downloads.value = _msg("error", "No downloads: undo the edit or Generate again.")
            return
        previous, self.package = self.package, export_package(result)
        if previous is not None:
            shutil.rmtree(previous.zip_path.parent, ignore_errors=True)
        self._show_downloads()

    def _undo(self, _=None):
        if not self.history:
            return
        self._set_result(self.history.pop())
        self.edit_info.value = _msg("info", "Last edit undone.")
        self.console.log("undo()", kind="echo", note="viewer")

    def _console_row(self):
        ai = ([self.assistant.console_ai, self.assistant.console_test]
              if getattr(self, "assistant", None) is not None else [])
        box = w.HBox([self.console_in, self.console_run] + ai,
                     layout=w.Layout(align_items="flex-start", margin="4px 0 0"))
        box.add_class("bb-console-row")
        return box

    def _fail(self, message: str, validation: bool = False):
        """A failed build: sad Buddy with why and how to fix, and the same in Status."""
        why, fix = failure_advice(message)
        self.last_failure = message
        self._sad("Oh no! Buddy could not build that.", why, fix)
        head = "Validation failed; downloads blocked." if validation else "Build failed."
        self.status.value = (_msg("error", [head, message]) +
                             f"<div class='nb-msg nb-msg-info'><ul><li><b>Why:</b> {html.escape(why)}"
                             f"</li><li><b>How to fix:</b> {html.escape(fix)}</li></ul></div>")
        if validation:
            self.downloads.value = _msg("error", "No downloads: fix the problem above and Generate "
                                                 "again.")

    def _sad(self, message: str = "Oh no! Buddy could not build that.", why: str = "",
             fix: str = ""):
        """Sad Buddy with the reason and the fix; fades out by itself (about 4 s)."""
        self._sad_count = getattr(self, "_sad_count", 0) + 1
        try:
            data = base64.b64encode((LOGO.parent / "buddy_sad.png").read_bytes()).decode("ascii")
        except OSError:
            return
        detail = ""
        if why:
            detail += f"<p><b>Why:</b> {html.escape(why)}</p>"
        if fix:
            detail += f"<p><b>How to fix:</b> {html.escape(fix)}</p>"
        self.popup.value = (f"<div class='bb-sad' data-n='{self._sad_count}'>"
                            f"<img src='data:image/png;base64,{data}' alt='Sad Buddy'>"
                            f"<div class='bb-sad-card'><div class='bb-sad-title'>{html.escape(message)}"
                            f"</div>{detail}</div></div>")

    def load_config(self, cfg: RunConfig, upload=None, note: str = ""):
        """Put a whole configuration into the widgets (examples, sessions, console settings())."""
        self._apply_config(cfg, upload=upload)
        self.result = self.package = None
        self.history = []
        self._stale = False
        self.undo_btn.layout.display = "none"
        self.edit_info.value = ""
        self.downloads.value = "<i>Settings loaded; click Generate model.</i>"
        try:
            self._preview_cell(quiet=True)
        except Exception:
            self.viewer.value = placeholder_html("Settings loaded: click Generate model.")
        self.status.value = _msg("info", note or "Settings loaded; click Generate model.")

    def _on_session_upload(self, _=None):
        import json
        file = _uploaded_file(self.session_upload)
        if file is None:
            return
        name, data = file
        text = data.decode("utf-8", errors="replace")
        if text.lstrip().startswith("{"):  # config.json from an export
            cfg = RunConfig.from_dict(json.loads(text))
            self.load_config(cfg, note=f"Loaded {name}; click Generate model.")
            self.console.log(f"# loaded settings from {name}", kind="echo")
            return
        self.console.log(f"# running session file {name}", kind="echo")
        self.console.run(text)

    def _on_bridge(self, change=None):
        """An edit request from a click in the viewer (JSON from the iframe)."""
        import json
        raw = (self.bridge.value or "").strip()
        if not raw:
            return
        try:
            msg = json.loads(raw)
        except ValueError:
            return
        view = msg.get("view")
        try:
            if msg.get("action") == "remove":
                self.apply_edit("remove", [int(msg["ligand"])], view)
                self.console.log(f"remove_ligand({int(msg['ligand'])})", kind="echo", note="viewer")
            elif msg.get("action") == "add":
                self.apply_edit("add", int(msg["atom"]), view)
                self.console.log(f"add_ligand({int(msg['atom'])})", kind="echo", note="viewer")
        except BuilderError as exc:
            self.edit_info.value = _msg("error", str(exc))

    def _on_console(self, change=None):
        text = self.console_in.value
        if not text.strip():
            return
        self.console.run(text)
        self.console_in.value = ""

    def _show_results(self, result):
        a, r = result.assembly, result.report
        rows = [
            ("Material", a.provenance["material"]["structure_label"]),
            ("Substrate atoms", a.n_substrate),
            ("Ligand atoms", r.counts["ligand_atoms"]),
            ("Ligands placed / requested", f"{a.achieved_ligands} / {a.target_ligands}"),
            ("Bonds / angles / dihedrals", f"{len(a.bond_list)} / {len(a.angle_list)} / {len(a.dihedral_list)}"),
            ("Coverage area", f"{a.surface_area_nm2:.3f} nm\u00b2"),
            ("Requested density", "-" if a.requested_density_nm2 is None else f"{a.requested_density_nm2:.3f} per nm\u00b2"),
            ("Achieved density", "-" if not a.target_ligands else f"{a.achieved_density_nm2:.3f} per nm\u00b2"),
            ("Chemistry", "-" if result.ligand is None else
             f"{result.ligand.smiles} via {TEMPLATES[result.ligand.binding_template].label}"),
            ("Closest nonbonded contact", "-" if r.min_nonbonded_distance_A is None
             else f"{r.min_nonbonded_distance_A:.2f} \u00c5"),
            ("Validation", r.status),
            ("MD export readiness", "Structure and topology; force field required"),
        ]
        carbon = a.provenance.get("carbon")
        if carbon:
            comp = carbon.get("composition_after_attachment") or carbon["composition"]
            groups = ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in carbon.get("groups", {}).items())
            rows[1:1] = [
                ("Carbon composition", f"C {comp['C']}, O {comp['O']}, H {comp['H']}"),
                ("Oxygen content (100 O/C)", f"{comp['oxygen_percent_O_per_C']} %"
                 + (f"  (C:O {comp['C_to_O_ratio']}:1)" if comp.get("C_to_O_ratio") else "")),
                ("Modified sp3 carbons", f"{comp['modified_carbon_percent']} %"),
                ("Vacancies", f"{carbon['composition']['vacancy_percent_actual']} %"),
                ("Oxygen groups", groups or "-"),
            ]
            ratio = (carbon.get("oxidation") or {}).get("achieved_hydroxyl_to_epoxy")
            if ratio:
                rows.insert(3, ("Basal OH : epoxide O", f"{ratio}:1"))
        body = "".join(f"<tr><td>{html.escape(k)}</td><td>{html.escape(str(v))}</td></tr>"
                       for k, v in rows)
        self.results.value = f"<table class='nb-table'>{body}</table>"

    def _show_downloads(self):
        pkg = self.package
        stage = pkg.stage_dir
        stem = Path(pkg.zip_name).stem
        links = [_download_link(pkg.zip_name, pkg.read_bytes(), "application/zip",
                                f"Download {stem}.zip", primary=True)]
        named = {"certificate_of_build.pdf": f"{stem}_certificate_of_build.pdf",
                 "session.py": f"{stem}_session.py"}
        for rel, mime, label in (("hoomd/structure.gsd", "application/octet-stream", "GSD"),
                                 ("lammps/structure.data", "text/plain", "LAMMPS data"),
                                 ("structure.extxyz", "text/plain", "extXYZ"),
                                 ("pdb/structure.pdb", "chemical/x-pdb", "PDB"),
                                 ("mol2/structure.mol2", "chemical/x-mol2", "MOL2"),
                                 ("certificate_of_build.pdf", "application/pdf",
                                  "Certificate of Build (PDF)"),
                                 ("session.py", "text/x-python", "Save session (.py)")):
            path = stage / rel
            if path.is_file():
                links.append(_download_link(named.get(rel, Path(rel).name), path.read_bytes(),
                                            mime, label))
        self.downloads.value = '<div class="nb-dl">' + "".join(links) + "</div>"

    def _inspect(self, _=None):
        if self.result is None:
            self.inspect_info.value = ""
            return
        a = self.result.assembly
        i = int(self.inspect_index.value)
        x, y, z = a.positions_A[i]
        self.inspect_info.value = (f"<small>{a.elements[i]} ({a.type_labels[i]}), role {a.chemical_role[i]}, "
                                   f"molecule {int(a.molecule_id[i])}, at ({x:.2f}, {y:.2f}, {z:.2f}) A</small>")

    def _on_example(self, _=None):
        key = self.example.value
        if not key:
            return
        cfg = RunConfig.from_dict(EXAMPLES[key]["config"])
        upload = example_upload(key)
        self._apply_config(cfg, upload=upload)
        self.result = self.package = None
        self._stale = False
        self.downloads.value = "<i>Example loaded; click Generate model.</i>"
        self.viewer.value = placeholder_html(f"{EXAMPLES[key]['label']}: click Generate model.")
        note = (f" The example uploads the bundled fixture {upload[0]}." if upload else "")
        self.status.value = _msg("info", f"Loaded example {EXAMPLES[key]['label']}.{note}")
