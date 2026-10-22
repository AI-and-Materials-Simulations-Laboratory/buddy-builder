"""Carbon-nanomaterial widgets for the GUI: material kind, geometry, defects and oxidation.

``CarbonPanel`` owns every carbon widget and converts them to and from a
:class:`CarbonConfig`; ``gui.BuilderApp`` places its three boxes on the
Material, Geometry and "Defects & GO" tabs and listens to ``widgets()``.
No science happens here.
"""

from __future__ import annotations

import html
from math import pi

import ipywidgets as w

from .carbon_lattice import AREA_PER_C, cnt_diameter
from .carbon_oxidation import AGED_PROFILES, group_fractions
from .schema import BuilderError, CarbonConfig

STYLE = {"description_width": "190px"}


def _wide():
    return w.Layout(width="100%", max_width="460px")


def _narrow():
    return w.Layout(width="300px", max_width="100%")


def _row(*children):
    return w.HBox(list(children), layout=w.Layout(flex_flow="row wrap", align_items="center"))


def _adv(title, children):
    return w.VBox([w.HTML(f"<div class='nb-adv-h'>{html.escape(title)}</div>")] + list(children))


def _tiny():
    return w.Layout(width="150px")


def _h(text):
    return w.HTML(f"<div class='nb-h'>{html.escape(text)}</div>")


FULLERENE_SIZES = [(f"C{60 * k * k} (k = {k})", k) for k in range(1, 11)]
PROFILE_OPTIONS = [("Fresh: hydroxyl + epoxy (ratio below)", "fresh"),
                   ("Aged, 600 K annealing profile (HierGO)", "aged_600K"),
                   ("Aged, 900 K annealing profile (HierGO)", "aged_900K"),
                   ("Aged, 1200 K annealing profile (HierGO)", "aged_1200K"),
                   ("Aged, 2500 K annealing profile (HierGO)", "aged_2500K"),
                   ("Custom mix of groups", "custom")]
GROUP_LABELS = {"hydroxyl": "Basal hydroxyl (C-OH)", "epoxy": "Epoxide (C-O-C)",
                "edge_hydroxyl": "Edge hydroxyl (phenol)", "carbonyl": "Edge carbonyl (C=O)",
                "carboxyl": "Carboxylic acid (COOH)", "ether": "In-plane ether (vacancy)"}


class CarbonPanel:
    """All carbon widgets plus conversion to/from CarbonConfig."""

    def __init__(self):
        # -- material ---------------------------------------------------------
        self.kind = w.Dropdown(options=[("Graphene / graphene oxide sheet", "graphene"),
                                        ("Carbon nanotube (single or multi-wall)", "nanotube"),
                                        ("Carbon sphere: fullerene or carbon onion", "carbon_sphere")],
                               description="Carbon material", style=STYLE, layout=_wide())
        self.kind_info = w.HTML()
        self.material_box = w.VBox([self.kind, self.kind_info])

        # -- geometry ---------------------------------------------------------
        self.length = w.BoundedFloatText(min=0.5, max=60, step=0.5, description="Length x, zigzag (nm)",
                                         style=STYLE, layout=_narrow())
        self.width = w.BoundedFloatText(min=0.5, max=60, step=0.5, description="Width y, armchair (nm)",
                                        style=STYLE, layout=_narrow())
        self.layers = w.BoundedIntText(min=1, max=20, description="Layers", style=STYLE, layout=_narrow())
        self.spacing = w.BoundedFloatText(min=3.0, max=20.0, step=0.05, description="Layer spacing (Å)",
                                          style=STYLE, layout=_narrow())
        self.stacking = w.Dropdown(options=[("AB (Bernal)", "AB"), ("AA", "AA")], description="Stacking",
                                   style=STYLE, layout=_narrow())
        self.sheet_box = w.VBox([_row(self.length, self.width), self.layers])
        self.cnt_n = w.BoundedIntText(min=1, max=200, description="Chirality n", style=STYLE,
                                      layout=_narrow())
        self.cnt_m = w.BoundedIntText(min=0, max=200, description="m", style={"description_width": "30px"},
                                      layout=w.Layout(width="130px"))
        self.cnt_len = w.BoundedFloatText(min=0.5, max=200, step=0.5, description="Tube length (nm)",
                                          style=STYLE, layout=_narrow())
        self.walls = w.BoundedIntText(min=1, max=8, description="Walls", style=STYLE, layout=_narrow())
        self.wall_spacing = w.BoundedFloatText(min=3.0, max=10.0, step=0.05,
                                               description="Wall spacing (Å)", style=STYLE,
                                               layout=_narrow())
        self.cnt_box = w.VBox([_row(self.cnt_n, self.cnt_m), self.cnt_len, self.walls])
        self.sphere_type = w.Dropdown(options=[("Single fullerene shell", "fullerene"),
                                               ("Carbon onion (nested shells)", "onion")],
                                      description="Sphere type", style=STYLE, layout=_wide())
        self.fullerene = w.Dropdown(options=FULLERENE_SIZES, description="Fullerene",
                                    style=STYLE, layout=_narrow())
        self.shells = w.BoundedIntText(min=1, max=8, description="Onion shells (C60 inside)",
                                       style=STYLE, layout=_narrow())
        self.sphere_box = w.VBox([self.sphere_type, self.fullerene, self.shells])
        self.periodic = w.Dropdown(options=[("Finite (non-periodic; AMBER-compatible)", False),
                                            ("Periodic (infinite sheet / tube)", True)],
                                   description="Boundary", style=STYLE, layout=_wide())
        self.edges = w.Dropdown(options=[("Hydrogen-terminated", "H"), ("Bare (dangling bonds)", "none")],
                                description="Edge termination", style=STYLE, layout=_wide())
        self.geom_info = w.HTML()
        self.geometry_adv = _adv("Advanced: carbon geometry", [self.stacking, self.spacing,
                                                               self.wall_spacing, self.edges])
        self.geometry_box = w.VBox([self.sheet_box, self.cnt_box, self.sphere_box, self.periodic,
                                    self.geom_info, self.geometry_adv])

        # -- defects ------------------------------------------------------------
        self.vacancy = w.FloatSlider(min=0, max=40, step=0.5, description="Vacancies (% of C)",
                                     style=STYLE, layout=_wide(), readout_format=".1f")
        self.hole_size = w.Dropdown(options=[("Mixed, 1-30 atoms", "mixed"), ("Small, 1-5", "small"),
                                             ("Medium, 5-10", "medium"), ("Large, 10-30", "large"),
                                             ("Custom range", "custom")],
                                    description="Hole sizes", style=STYLE, layout=_wide())
        self.hole_min = w.BoundedIntText(min=1, max=500, description="Min atoms", style=STYLE,
                                         layout=_narrow())
        self.hole_max = w.BoundedIntText(min=1, max=500, description="Max", style={"description_width": "40px"},
                                         layout=w.Layout(width="140px"))
        self.hole_shapes = w.Dropdown(options=[("Random growth (HierGO create-tile)", "random_growth"),
                                               ("HierGO hole library shapes", "hiergo_library")],
                                      description="Hole shapes", style=STYLE, layout=_wide())
        self.pores = w.BoundedIntText(min=0, max=100, description="Circular pores", style=STYLE,
                                      layout=_narrow())
        self.pore_d = w.BoundedFloatText(min=0.2, max=20, step=0.1, description="d (nm)",
                                         style={"description_width": "45px"}, layout=w.Layout(width="140px"))
        self.sw = w.FloatSlider(min=0, max=20, step=0.5, description="Stone-Wales (% of C)",
                                style=STYLE, layout=_wide(), readout_format=".1f")
        self.flower = w.FloatSlider(min=0, max=20, step=0.5, description="Mini-GB, 24-atom rot. (%)",
                                    style=STYLE, layout=_wide(), readout_format=".1f")
        self.ruptures = w.BoundedIntText(min=0, max=50, description="Ruptures", style=STYLE,
                                         layout=_narrow())
        self.rupture_len = w.BoundedIntText(min=1, max=50, description="bonds",
                                            style={"description_width": "45px"},
                                            layout=w.Layout(width="140px"))
        self.keep_hanging = w.Checkbox(description="Keep hanging carbons (skip HierGO trimming)",
                                       indent=False)
        self.defect_info = w.HTML()
        self.defects_adv = _adv("Advanced: defects", [
            _row(self.hole_min, self.hole_max), self.hole_shapes, self.flower,
            _row(self.ruptures, self.rupture_len), self.keep_hanging])
        self.defects_box = w.VBox([
            _h("Defects (HierGO)"), self.vacancy, self.hole_size,
            _row(self.pores, self.pore_d), self.sw, self.defect_info, self.defects_adv])

        # -- oxidation -------------------------------------------------------------
        self.oxygen = w.FloatSlider(min=0, max=80, step=0.5, description="Oxygen (100 O/C, %)",
                                    style=STYLE, layout=_wide(), readout_format=".1f")
        self.profile = w.Dropdown(options=PROFILE_OPTIONS, description="Oxidation profile", style=STYLE,
                                  layout=_wide())
        self.ratio = w.Text(description="Hydroxyl : epoxy", style=STYLE, layout=_narrow(),
                            placeholder="2:1", continuous_update=False)
        self.paired = w.FloatSlider(min=0, max=1, step=0.05, description="Paired (anti diol) share",
                                    style=STYLE, layout=_wide(), readout_format=".2f")
        self.custom = {k: w.BoundedFloatText(min=0, max=100, step=1, description=label, style=STYLE,
                                             layout=_narrow()) for k, label in GROUP_LABELS.items()}
        self.custom_box = w.VBox([w.HTML("<small>Relative shares of the added oxygen atoms (a "
                                         "carboxylic acid counts as two).</small>")]
                                 + list(self.custom.values()))
        self.pattern = w.Dropdown(options=[("Uniform random", "uniform"),
                                           ("Edges and holes first", "edges_first"),
                                           ("Clustered domains (HierGO bunching)", "domains"),
                                           ("Tiles with different O% (HierGO tiling)", "tiles")],
                                  description="Spatial pattern", style=STYLE, layout=_wide())
        self.domains = w.BoundedIntText(min=1, max=100, description="Domains", style=STYLE,
                                        layout=_narrow())
        self.domain_d = w.BoundedFloatText(min=0.3, max=50, step=0.1, description="d (nm)",
                                           style={"description_width": "45px"},
                                           layout=w.Layout(width="140px"))
        self.tile_x = w.BoundedIntText(min=1, max=20, description="Tiles x × y", style=STYLE,
                                       layout=w.Layout(width="260px"))
        self.tile_y = w.BoundedIntText(min=1, max=20, layout=w.Layout(width="70px"))
        self.tile_min = w.BoundedFloatText(min=0, max=100, step=1, description="Tile O% min",
                                           style=STYLE, layout=_narrow())
        self.tile_max = w.BoundedFloatText(min=0, max=100, step=1, description="max",
                                           style={"description_width": "40px"},
                                           layout=w.Layout(width="140px"))
        self.domain_box = _row(self.domains, self.domain_d)
        self.tile_box = w.VBox([_row(self.tile_x, self.tile_y), _row(self.tile_min, self.tile_max)])
        self.faces = w.Dropdown(options=[("Both faces, balanced", "both"), ("Top face only", "top"),
                                         ("Bottom face only", "bottom")],
                                description="Oxidized faces", style=STYLE, layout=_wide())
        self.all_layers = w.Checkbox(description="Oxidize every layer (needs >= 6 Å spacing)",
                                     indent=False)
        self.pyramidalize = w.Checkbox(description="Pucker sp3 carbons toward their group "
                                                   "(0.25 Å; 0.40 Å at ligand links)",
                                       indent=False)
        self.seed = w.IntText(description="Random seed (carbon)", style=STYLE, layout=_narrow())
        self.ox_info = w.HTML()
        self.oxidation_adv = _adv("Advanced: oxidation", [
            self.paired, self.domain_box, self.tile_box, self.faces, self.all_layers,
            self.pyramidalize, self.seed])
        self.oxidation_box = w.VBox([
            _h("Oxidation (graphene oxide, HierGO decorate-tile)"), self.oxygen, self.profile,
            self.ratio, self.custom_box, self.pattern, self.ox_info, self.oxidation_adv])
        self.not_carbon = w.HTML(
            "<div class='nb-msg nb-msg-info'><ul><li>Defects and oxidation apply to carbon "
            "nanomaterials. Choose <b>Carbon nanomaterial</b> as the material source on tab 1 "
            "(or pick one from the Carbon group of Find material).</li></ul></div>")
        self.tab = w.VBox([self.not_carbon, self.defects_box, self.oxidation_box])

    # ------------------------------------------------------------------ plumbing
    def widgets(self) -> list:
        return ([self.kind, self.length, self.width, self.layers, self.spacing, self.stacking,
                 self.cnt_n, self.cnt_m, self.cnt_len, self.walls, self.wall_spacing,
                 self.sphere_type, self.fullerene, self.shells, self.periodic, self.edges,
                 self.vacancy, self.hole_size, self.hole_min, self.hole_max, self.hole_shapes,
                 self.pores, self.pore_d, self.sw, self.flower, self.ruptures, self.rupture_len,
                 self.keep_hanging, self.oxygen, self.profile, self.ratio, self.paired,
                 self.pattern, self.domains, self.domain_d, self.tile_x, self.tile_y,
                 self.tile_min, self.tile_max, self.faces, self.all_layers, self.pyramidalize,
                 self.seed] + list(self.custom.values()))

    def config(self) -> CarbonConfig:
        return CarbonConfig(
            kind=self.kind.value, periodic=bool(self.periodic.value),
            edge_termination=self.edges.value, length_nm=self.length.value,
            width_nm=self.width.value, layers=self.layers.value,
            layer_spacing_A=self.spacing.value, stacking=self.stacking.value,
            cnt_n=self.cnt_n.value, cnt_m=self.cnt_m.value, cnt_length_nm=self.cnt_len.value,
            cnt_walls=self.walls.value, wall_spacing_A=self.wall_spacing.value,
            sphere_type=self.sphere_type.value, fullerene_order=self.fullerene.value,
            onion_shells=self.shells.value, vacancy_percent=self.vacancy.value,
            hole_size=self.hole_size.value, hole_min_atoms=self.hole_min.value,
            hole_max_atoms=self.hole_max.value, hole_shapes=self.hole_shapes.value,
            pore_count=self.pores.value, pore_diameter_nm=self.pore_d.value,
            stone_wales_percent=self.sw.value, flower_percent=self.flower.value,
            rupture_count=self.ruptures.value, rupture_length=self.rupture_len.value,
            keep_hanging_carbons=self.keep_hanging.value, oxygen_percent=self.oxygen.value,
            oxidation_profile=self.profile.value, hydroxyl_epoxy_ratio=self.ratio.value.strip() or "2:1",
            paired_hydroxyl_fraction=self.paired.value,
            custom_fractions={k: float(cb.value) / 100.0 for k, cb in self.custom.items()},
            oxidation_pattern=self.pattern.value, oxidation_faces=self.faces.value,
            domain_count=self.domains.value, domain_diameter_nm=self.domain_d.value,
            tile_grid=(self.tile_x.value, self.tile_y.value), tile_oxygen_min=self.tile_min.value,
            tile_oxygen_max=self.tile_max.value, oxidize_all_layers=self.all_layers.value,
            pyramidalize=self.pyramidalize.value, seed=self.seed.value)

    def apply(self, c: CarbonConfig):
        self.kind.value = c.kind
        self.periodic.value = bool(c.periodic)
        self.edges.value = c.edge_termination
        self.length.value, self.width.value = c.length_nm, c.width_nm
        self.layers.value, self.spacing.value, self.stacking.value = c.layers, c.layer_spacing_A, c.stacking
        self.cnt_n.value, self.cnt_m.value = c.cnt_n, c.cnt_m
        self.cnt_len.value, self.walls.value, self.wall_spacing.value = (c.cnt_length_nm, c.cnt_walls,
                                                                         c.wall_spacing_A)
        self.sphere_type.value, self.fullerene.value = c.sphere_type, c.fullerene_order
        self.shells.value = c.onion_shells
        self.vacancy.value, self.hole_size.value = c.vacancy_percent, c.hole_size
        self.hole_min.value, self.hole_max.value = c.hole_min_atoms, c.hole_max_atoms
        self.hole_shapes.value = c.hole_shapes
        self.pores.value, self.pore_d.value = c.pore_count, c.pore_diameter_nm
        self.sw.value, self.flower.value = c.stone_wales_percent, c.flower_percent
        self.ruptures.value, self.rupture_len.value = c.rupture_count, c.rupture_length
        self.keep_hanging.value = c.keep_hanging_carbons
        self.oxygen.value, self.profile.value = c.oxygen_percent, c.oxidation_profile
        self.ratio.value, self.paired.value = c.hydroxyl_epoxy_ratio, c.paired_hydroxyl_fraction
        for k, cb in self.custom.items():
            cb.value = round(100.0 * float(c.custom_fractions.get(k, 0.0)), 3)
        self.pattern.value, self.faces.value = c.oxidation_pattern, c.oxidation_faces
        self.domains.value, self.domain_d.value = c.domain_count, c.domain_diameter_nm
        self.tile_x.value, self.tile_y.value = (int(i) for i in c.tile_grid)
        self.tile_min.value, self.tile_max.value = c.tile_oxygen_min, c.tile_oxygen_max
        self.all_layers.value, self.pyramidalize.value = c.oxidize_all_layers, c.pyramidalize
        self.seed.value = c.seed

    def shape(self) -> str:
        return {"graphene": "sheet", "nanotube": "nanotube", "carbon_sphere": "carbon_sphere"}[self.kind.value]

    # ------------------------------------------------------------------ display
    user_box = False  # set by the app: a user box decides the size of periodic sheets/tubes

    def update(self, active: bool, advanced: bool = True):
        def show(widget, visible):
            widget.layout.display = "" if visible else "none"

        for box in (self.geometry_adv, self.defects_adv, self.oxidation_adv):
            show(box, advanced)

        kind = self.kind.value
        sheet, tube, sphere = kind == "graphene", kind == "nanotube", kind == "carbon_sphere"
        show(self.material_box, active)
        show(self.geometry_box, active)
        show(self.defects_box, active)
        show(self.oxidation_box, active)
        show(self.not_carbon, not active)
        show(self.sheet_box, sheet)
        show(self.cnt_box, tube)
        from_box = self.user_box and bool(self.periodic.value)
        show(self.length, not (sheet and from_box))
        show(self.width, not (sheet and from_box))
        show(self.cnt_len, not (tube and from_box))
        show(self.sphere_box, sphere)
        show(self.periodic, not sphere)
        show(self.edges, not sphere)
        show(self.stacking, sheet)
        show(self.spacing, sheet)
        show(self.wall_spacing, tube)
        show(self.fullerene, self.sphere_type.value == "fullerene")
        show(self.shells, self.sphere_type.value == "onion")
        show(self.hole_min, self.hole_size.value == "custom")
        show(self.hole_max, self.hole_size.value == "custom")
        if self.hole_size.value == "custom" and not advanced:
            show(self.defects_adv, True)  # the custom range lives there
        for widget in (self.pores, self.pore_d, self.flower, self.ruptures, self.rupture_len):
            show(widget, sheet)
        show(self.sw, not sphere)
        show(self.ratio, self.profile.value == "fresh")
        show(self.paired, self.profile.value in ("fresh", "custom"))
        show(self.custom_box, self.profile.value == "custom")
        show(self.domain_box, self.pattern.value == "domains")
        show(self.tile_box, self.pattern.value == "tiles")
        show(self.faces, sheet)
        show(self.all_layers, sheet and self.layers.value > 1)
        self.periodic.options = ([("Finite flake (non-periodic; AMBER-compatible)", False),
                                  ("Periodic sheet (infinite in x and y)", True)] if sheet else
                                 [("Finite tube (open, hydrogen-capped rims)", False),
                                  ("Periodic along the tube axis (infinite)", True)])
        if active:
            self.kind_info.value = self._kind_text()
            self.geom_info.value = self._geometry_text()
            self.defect_info.value = self._defect_text()
            self.ox_info.value = self._oxidation_text()

    def _kind_text(self) -> str:
        notes = {
            "graphene": "Orthogonal honeycomb sheet (zigzag edge along x, armchair along y), one or "
                        "more layers. Add holes, topological defects and ruptures, and oxygen groups "
                        "(graphene oxide) on the Defects &amp; GO tab.",
            "nanotube": "Rolled graphene of chirality (n, m), open hydrogen-capped rims or "
                        "periodic along z; multi-wall tubes add concentric walls about 3.4 Å apart.",
            "carbon_sphere": "Icosahedral fullerene C60 k² (Goldberg polyhedron, relaxed to "
                             "a_CC = 1.42 Å) or a carbon onion of nested shells "
                             "C60@C240@C540...",
        }
        return f"<small>{notes[self.kind.value]}</small>"

    def _geometry_text(self) -> str:
        k = self.kind.value
        try:
            if k == "graphene":
                n = self.length.value * self.width.value * 100 / AREA_PER_C * self.layers.value
                text = f"About {n:.0f} carbon atoms before defects."
            elif k == "nanotube":
                n_, m_ = self.cnt_n.value, self.cnt_m.value
                if m_ > n_:
                    return "<small>Use 0 ≤ m ≤ n.</small>"
                d = cnt_diameter(n_, m_)
                kind = "armchair" if n_ == m_ else ("zigzag" if m_ == 0 else "chiral")
                elec = "metallic" if (n_ - m_) % 3 == 0 else "semiconducting"
                atoms = pi * d * self.cnt_len.value * 10 / AREA_PER_C
                text = (f"({n_},{m_}) {kind}, diameter {d / 10:.3f} nm, nominally {elec}; about "
                        f"{atoms:.0f} atoms in the inner wall.")
                if self.periodic.value and self.walls.value > 1 and n_ != m_ and m_ != 0:
                    text += " Periodic multi-wall tubes need armchair or zigzag walls."
            else:
                if self.sphere_type.value == "fullerene":
                    k_ = self.fullerene.value
                    text = f"C{60 * k_ * k_}: radius about {0.352 * k_:.2f} nm."
                else:
                    s = self.shells.value
                    total = sum(60 * i * i for i in range(1, s + 1))
                    text = (f"{s} shells, {total} atoms, outer radius about {0.352 * s:.2f} nm "
                            "(shells about 3.5 Å apart).")
        except (BuilderError, ValueError, ZeroDivisionError):
            return ""
        return f"<small>{text}</small>"

    def _defect_text(self) -> str:
        if self.kind.value == "carbon_sphere":
            return "<small>Fullerene shells take vacancies only.</small>"
        parts = []
        if self.sw.value or self.flower.value:
            parts.append("Defect counts follow HierGO: one Stone-Wales per 16 atoms and one 24-atom "
                         "rotation per 41 atoms for each percent.")
        if self.vacancy.value:
            parts.append("Holes stay apart and away from edges. With hydrogen-terminated edges, holes "
                         "too small to hold their edge H (e.g. monovacancies) grow while cleaning; "
                         "removed atoms count toward the vacancy target.")
        if self.hole_shapes.value == "hiergo_library":
            parts.append("Library shapes are read from HierGO-master/HoleLibrary next to the tool "
                         "(GPL-3.0, not bundled).")
        return f"<small>{' '.join(parts)}</small>" if parts else ""

    def _oxidation_text(self) -> str:
        try:
            cfg = self.config()
            if cfg.oxidation_profile == "custom" and sum(cfg.custom_fractions.values()) <= 0:
                return "<small>Give at least one group a non-zero share.</small>"
            f, paired = group_fractions(cfg)
        except BuilderError as exc:
            return f"<small>{html.escape(str(exc))}</small>"
        names = [("Basal hydroxyl", f["hydroxyl"] + f["paired"]), ("Epoxide", f["epoxy"]),
                 ("Edge hydroxyl", f["edge_hydroxyl"]), ("Carbonyl", f["carbonyl"]),
                 ("Carboxylic acid", f["carboxyl"]), ("Ether", f["ether"])]
        rows = "".join(f"<tr><td>{n}</td><td>{100 * v:.1f} %</td></tr>" for n, v in names if v > 0)
        note = ("Oxygen content is 100 N_O/N_C of the whole substrate. Edge groups need edges or "
                "holes; groups without a site become basal hydroxyls, and the achieved composition "
                "is reported after Generate.")
        if self.kind.value != "graphene":
            note += " Curved walls are oxidized on their outer face only."
        return (f"<small>Share of added oxygen atoms ({100 * paired:.0f} % of hydroxyls as anti "
                f"pairs):</small><table class='nb-table'>{rows}</table><small>{note}</small>")

    def nominal_area_nm2(self, region: str):
        """Rough area for the coverage tab before a substrate is built (None if unknown)."""
        k = self.kind.value
        if k == "graphene":
            face = self.length.value * self.width.value
            return {"top": face, "bottom": face, "both": 2 * face, "entire_surface": 2 * face}.get(region)
        if k == "nanotube":
            try:
                from .carbon_lattice import wall_chiralities
                walls = wall_chiralities(self.cnt_n.value, self.cnt_m.value, self.walls.value,
                                         self.wall_spacing.value, False)
                d = cnt_diameter(*walls[-1]) / 10
            except (BuilderError, ValueError):
                return None
            return {"side": pi * d * self.cnt_len.value, "entire_surface": pi * d * self.cnt_len.value}.get(region)
        orders = (self.fullerene.value if self.sphere_type.value == "fullerene" else self.shells.value)
        r = 0.352 * orders
        return 4 * pi * r * r if region == "entire_surface" else None


def aged_profile_table() -> dict:
    """For docs/tests: normalized aged-profile fractions."""
    out = {}
    for key, raw in AGED_PROFILES.items():
        total = sum(raw.values())
        out[key] = {k: round(v / total, 4) for k, v in raw.items()}
    return out


__all__ = ["CarbonPanel", "aged_profile_table"]
