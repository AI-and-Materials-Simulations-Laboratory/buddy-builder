"""Shared configuration and model contracts (spec section 5).

Every module imports its data structures from here. Field names follow the
spec's tables; change them only through the team's schema-change protocol.

Units: lengths in Angstrom unless a field name says otherwise (``_nm``,
``_nm2``), masses in amu, charges in elementary charges.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional

import numpy as np
from ase import Atoms

FORMAT_VERSION = "0.1"

SHAPES = ("sphere", "slab", "rod", "channel")
CHANNEL_TYPES = ("slit", "pore")  # crystal channels: stacked slit pores, cylindrical pore
BOX_SHAPES = ("auto", "rectangular", "cube")
CARBON_SHAPES = ("sheet", "nanotube", "carbon_sphere")  # SubstrateModel.shape of carbon builds
ROD_CAPS = ("hemispherical", "flat")
THICKNESS_MODES = ("layers", "thickness")
MATERIAL_SOURCES = ("preset", "uploaded_file", "materials_project", "carbon")
UPLOAD_FORMATS = ("cif", "vasp")
BINDING_TEMPLATES = ("Au_thiolate", "geometric", "chemisorption", "C_basal", "C_edge", "C_carboxyl", "C_epoxide")
ORIENTATIONS = ("outward", "random_azimuth")
CONFORMER_MODES = ("etkdg", "extended")
COVERAGE_MODES = ("count", "density")
REGIONS = ("entire_surface", "top", "bottom", "both", "side", "tips", "edges", "inner")
DISTRIBUTIONS = ("distributed_random", "evenly_spaced", "fibonacci")
EXPORT_TARGETS = ("hoomd_gsd", "lammps_data", "extxyz", "pdb", "mol2")
FIDELITIES = ("structure_topology", "simulation_ready_validated")
PARAMETER_STATUSES = ("unassigned", "partial", "complete_validated")
RESIDUE_SCHEMES = ("groups", "single")

# Carbon nanomaterials (graphene, graphene oxide, nanotubes, fullerenes/onions)
CARBON_KINDS = ("graphene", "nanotube", "carbon_sphere")
EDGE_TERMINATIONS = ("H", "none")
STACKINGS = ("AB", "AA")
SPHERE_TYPES = ("fullerene", "onion")
HOLE_SIZES = ("small", "medium", "large", "mixed", "custom")
HOLE_SIZE_RANGES = {"small": (1, 5), "medium": (5, 10), "large": (10, 30), "mixed": (1, 30)}
HOLE_SHAPES = ("random_growth", "hiergo_library")
OXIDATION_PROFILES = ("fresh", "aged_600K", "aged_900K", "aged_1200K", "aged_2500K", "custom")
OXIDATION_PATTERNS = ("uniform", "edges_first", "domains", "tiles")
OXIDATION_FACES = ("both", "top", "bottom")
OXYGEN_GROUPS = ("hydroxyl", "epoxy", "edge_hydroxyl", "carbonyl", "carboxyl", "ether")


class BuilderError(ValueError):
    """An input or scientific problem with a message meant for the GUI user."""


def _check_choice(name: str, value: Any, choices: tuple) -> None:
    if value not in choices:
        raise BuilderError(f"{name} must be one of {', '.join(map(str, choices))}; got {value!r}")


def _positive(name: str, value: float) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        raise BuilderError(f"{name} must be a number; got {value!r}") from None
    if not np.isfinite(value) or value <= 0:
        raise BuilderError(f"{name} must be positive; got {value}")
    return value


# --------------------------------------------------------------------------
# Configuration (what the GUI produces and config.json stores)
# --------------------------------------------------------------------------


@dataclass
class MaterialConfig:
    source: str = "preset"
    preset_id: Optional[str] = "Au_fcc"
    lattice_overrides_A: dict = field(default_factory=dict)
    file_format: Optional[str] = None
    filename: Optional[str] = None
    sha256: Optional[str] = None
    mp_id: Optional[str] = None  # Materials Project entry (source "materials_project")
    # Raw uploaded bytes; never written to config.json (the file goes into the ZIP).
    # For Materials Project this is the fetched conventional cell as CIF.
    uploaded_bytes: Optional[bytes] = field(default=None, repr=False, compare=False)

    def validate(self) -> None:
        _check_choice("Material source", self.source, MATERIAL_SOURCES)
        if self.source == "carbon":
            return  # everything carbon-specific lives in RunConfig.carbon
        if self.source == "preset":
            if not self.preset_id:
                raise BuilderError("Choose a material preset")
            for key, value in self.lattice_overrides_A.items():
                _positive(f"Lattice parameter {key}", value)
        elif self.source == "materials_project":
            if not self.mp_id:
                raise BuilderError("Choose a structure from the Materials Project library")
            if self.uploaded_bytes is None:
                from .mp_library import cif_bytes  # preloaded, offline structures
                self.uploaded_bytes = cif_bytes(self.mp_id)
            if self.uploaded_bytes is None:
                raise BuilderError(f"{self.mp_id} is not in the preloaded Materials Project "
                                   "library; upload its CIF instead (or pass it with --crystal)")
        else:
            _check_choice("Upload format", self.file_format, UPLOAD_FORMATS)
            if self.uploaded_bytes is None:
                raise BuilderError("Upload a crystallographic file")

    def to_dict(self) -> dict:
        if self.source == "carbon":
            return {"source": "carbon"}
        d = asdict(self)
        d.pop("uploaded_bytes")
        if self.source != "materials_project":
            d.pop("mp_id")
        if self.source == "preset":
            for key in ("file_format", "filename", "sha256"):
                d.pop(key)
        else:
            d["preset_id"] = None
            d.pop("lattice_overrides_A")
        return d


@dataclass
class GeometryConfig:
    shape: str = "sphere"
    diameter_nm: float = 3.0
    slab_length_nm: float = 5.0
    slab_width_nm: float = 5.0
    thickness_mode: str = "layers"  # "layers" takes precedence over thickness_nm
    layers: int = 4
    thickness_nm: float = 1.0
    miller: tuple = (1, 1, 1)
    termination_index: int = 0
    vacuum_A: float = 30.0
    max_atoms: int = 100000  # substrate atoms; see docs/size_limits.md
    # Nanorod: diameter_nm is the rod diameter; rod_length_nm is tip to tip.
    rod_length_nm: float = 10.0
    rod_axis: tuple = (0, 0, 1)  # direction [uvw] in the conventional cell (Au rods grow along <100>)
    rod_caps: str = "hemispherical"
    # Channels (crystal): slit = channel_count stacked slit pores between slab walls;
    # pore = one cylindrical pore through the slab (a membrane). Walls use the slab settings.
    channel_type: str = "slit"
    channel_width_nm: float = 1.5  # slit gap or pore diameter
    channel_count: int = 1
    # Simulation box: auto (lattice cell / fitted box), rectangular (orthogonal periodic cell,
    # edge atoms trimmed to fit), cube (square periodic cell / cubic box); optional minimum
    # box lengths in nm (x, y, z) for the non-periodic directions.
    box_shape: str = "auto"
    box_lengths_nm: Optional[tuple] = None

    @property
    def pbc(self) -> tuple:
        return (True, True, False) if self.shape in ("slab", "channel") else (False, False, False)

    def validate(self) -> None:
        _check_choice("Shape", self.shape, SHAPES)
        _positive("Vacuum padding", self.vacuum_A)
        _check_choice("Box shape", self.box_shape, BOX_SHAPES)
        if self.box_lengths_nm is not None:
            lengths = [float(x) for x in self.box_lengths_nm]
            if len(lengths) != 3 or any(not np.isfinite(x) or x < 0 for x in lengths):
                raise BuilderError("Box lengths must be three non-negative numbers (nm); 0 = automatic")
            if any(x > 200 for x in lengths):
                raise BuilderError("Box edges are limited to 200 nm")
            if self.box_shape == "cube" and len({round(x, 6) for x in lengths if x > 0}) > 1:
                raise BuilderError("A cubic box has equal edges: give one edge length")
        if self.shape == "channel":
            _check_choice("Channel type", self.channel_type, CHANNEL_TYPES)
            _positive("Channel width", self.channel_width_nm)
            if not 1 <= int(self.channel_count) <= 20:
                raise BuilderError("Number of channels must be between 1 and 20")
        if int(self.max_atoms) < 1:
            raise BuilderError("Max atoms must be at least 1")
        if self.shape == "sphere":
            _positive("Particle diameter", self.diameter_nm)
        elif self.shape == "rod":
            _positive("Rod diameter", self.diameter_nm)
            _positive("Rod length", self.rod_length_nm)
            _check_choice("Rod caps", self.rod_caps, ROD_CAPS)
            if self.rod_caps == "hemispherical" and self.rod_length_nm < self.diameter_nm:
                raise BuilderError("A rod with hemispherical caps must be at least as long as it "
                                   "is wide (length is measured tip to tip)")
            axis = tuple(int(i) for i in self.rod_axis)
            if len(axis) != 3 or axis == (0, 0, 0):
                raise BuilderError("Rod axis [u v w] must be three integers, not [0 0 0]")
        else:
            _positive("Slab length", self.slab_length_nm)
            _positive("Slab width", self.slab_width_nm)
            _check_choice("Thickness mode", self.thickness_mode, THICKNESS_MODES)
            if self.thickness_mode == "layers":
                if int(self.layers) < 1:
                    raise BuilderError("Layer count must be at least 1")
            else:
                _positive("Slab thickness", self.thickness_nm)
            miller = tuple(int(i) for i in self.miller)
            if len(miller) != 3 or miller == (0, 0, 0):
                raise BuilderError("Miller indices (h, k, l) must be three integers, not (0, 0, 0)")
            if int(self.termination_index) < 0:
                raise BuilderError("Termination index must be non-negative")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["miller"] = [int(i) for i in self.miller]
        d["pbc"] = list(self.pbc)
        d["rod_axis"] = [int(i) for i in self.rod_axis]
        if self.box_lengths_nm is not None:
            d["box_lengths_nm"] = [float(x) for x in self.box_lengths_nm]
        slab_keys = ("slab_length_nm", "slab_width_nm", "thickness_mode", "layers",
                     "thickness_nm", "miller", "termination_index")
        rod_keys = ("rod_length_nm", "rod_axis", "rod_caps")
        if self.shape in ("sphere", "rod"):
            for key in slab_keys:
                d.pop(key)
        else:
            d.pop("diameter_nm")
        if self.shape != "rod":
            for key in rod_keys:
                d.pop(key)
        if self.shape != "channel":
            for key in ("channel_type", "channel_width_nm", "channel_count"):
                d.pop(key)
        return d


def parse_ratio(text) -> tuple:
    """'2:1' -> (2.0, 1.0); both parts non-negative, not both zero."""
    try:
        parts = [float(p) for p in str(text).replace("/", ":").split(":")]
    except ValueError:
        parts = []
    if len(parts) != 2 or min(parts) < 0 or sum(parts) <= 0 or not all(np.isfinite(parts)):
        raise BuilderError(f"Hydroxyl:epoxy ratio must look like 2:1; got {text!r}")
    return parts[0], parts[1]


def _default_custom_fractions() -> dict:
    return {"hydroxyl": 0.40, "epoxy": 0.30, "edge_hydroxyl": 0.05, "carbonyl": 0.10,
            "carboxyl": 0.10, "ether": 0.05}


@dataclass
class CarbonConfig:
    """Carbon nanomaterials: graphene / graphene oxide sheets, nanotubes, fullerenes and onions.

    Defect and oxidation options follow the HierGO workflow (Garcia et al., 2D Mater. 10,
    025007, 2023) and its AMBER adaptation (Fedai et al., Data 11, 18, 2026); see
    docs/carbon_materials.md. Oxygen content is 100 * N_O / N_C of the substrate.
    """

    kind: str = "graphene"
    periodic: bool = False  # sheet: periodic in-plane; nanotube: periodic along its axis
    edge_termination: str = "H"
    # graphene sheet (zigzag edge along x, armchair edge along y)
    length_nm: float = 5.0
    width_nm: float = 5.0
    layers: int = 1
    layer_spacing_A: float = 3.35
    stacking: str = "AB"
    # nanotube
    cnt_n: int = 10
    cnt_m: int = 10
    cnt_length_nm: float = 5.0
    cnt_walls: int = 1
    wall_spacing_A: float = 3.4
    # fullerene / carbon onion (icosahedral Goldberg shells C60 k^2)
    sphere_type: str = "fullerene"
    fullerene_order: int = 1  # k: 1 = C60, 2 = C240, 3 = C540, ...
    onion_shells: int = 3
    # defects
    vacancy_percent: float = 0.0
    hole_size: str = "mixed"
    hole_min_atoms: int = 1
    hole_max_atoms: int = 30
    hole_shapes: str = "random_growth"
    pore_count: int = 0
    pore_diameter_nm: float = 1.0
    stone_wales_percent: float = 0.0
    flower_percent: float = 0.0  # 24-atom rotations ("mini grain boundary" in HierGO)
    rupture_count: int = 0
    rupture_length: int = 4  # broken C-C bonds (ether bridges) per rupture
    keep_hanging_carbons: bool = False
    # oxidation (graphene oxide)
    oxygen_percent: float = 0.0
    oxidation_profile: str = "fresh"
    hydroxyl_epoxy_ratio: str = "2:1"
    paired_hydroxyl_fraction: float = 0.7  # HierGO "dihydroxyl": anti 1,2-diols on adjacent carbons
    custom_fractions: dict = field(default_factory=_default_custom_fractions)
    oxidation_pattern: str = "uniform"
    oxidation_faces: str = "both"
    domain_count: int = 3
    domain_diameter_nm: float = 2.0
    tile_grid: tuple = (3, 3)
    tile_oxygen_min: float = 0.0
    tile_oxygen_max: float = 60.0
    oxidize_all_layers: bool = False  # multilayer sheets: False = outer faces only
    pyramidalize: bool = True  # move sp3 carbons 0.25 A toward their substituent
    seed: int = 7

    def validate(self) -> None:
        _check_choice("Carbon material", self.kind, CARBON_KINDS)
        _check_choice("Edge termination", self.edge_termination, EDGE_TERMINATIONS)
        if self.kind == "graphene":
            _positive("Sheet length", self.length_nm)
            _positive("Sheet width", self.width_nm)
            if not 1 <= int(self.layers) <= 20:
                raise BuilderError("Number of graphene layers must be between 1 and 20")
            _check_choice("Stacking", self.stacking, STACKINGS)
            if float(self.layer_spacing_A) < 3.0:
                raise BuilderError("Layer spacing must be at least 3.0 A")
            if self.length_nm < 0.5 or self.width_nm < 0.5:
                raise BuilderError("A graphene sheet must be at least 0.5 nm in each direction")
        elif self.kind == "nanotube":
            n, m = int(self.cnt_n), int(self.cnt_m)
            if n < 1 or m < 0 or m > n:
                raise BuilderError("Nanotube chirality (n, m) needs n >= 1 and 0 <= m <= n")
            if n + m < 4:
                raise BuilderError("Nanotube (n, m) is too narrow; use n + m >= 4")
            _positive("Nanotube length", self.cnt_length_nm)
            if not 1 <= int(self.cnt_walls) <= 8:
                raise BuilderError("Number of nanotube walls must be between 1 and 8")
            if float(self.wall_spacing_A) < 3.0:
                raise BuilderError("Wall spacing must be at least 3.0 A")
        else:
            _check_choice("Carbon sphere type", self.sphere_type, SPHERE_TYPES)
            if not 1 <= int(self.fullerene_order) <= 10:
                raise BuilderError("Fullerene order k must be between 1 (C60) and 10 (C6000)")
            if not 1 <= int(self.onion_shells) <= 8:
                raise BuilderError("A carbon onion has 1 to 8 shells")
        for name, value, hi in (("Vacancy percent", self.vacancy_percent, 50.0),
                                ("Stone-Wales percent", self.stone_wales_percent, 50.0),
                                ("24-atom rotation percent", self.flower_percent, 50.0),
                                ("Oxygen percent", self.oxygen_percent, 100.0),
                                ("Tile oxygen minimum", self.tile_oxygen_min, 100.0),
                                ("Tile oxygen maximum", self.tile_oxygen_max, 100.0)):
            if not np.isfinite(float(value)) or not 0 <= float(value) <= hi:
                raise BuilderError(f"{name} must be between 0 and {hi:g}; got {value}")
        _check_choice("Hole size distribution", self.hole_size, HOLE_SIZES)
        _check_choice("Hole shapes", self.hole_shapes, HOLE_SHAPES)
        if self.hole_size == "custom" and not 1 <= int(self.hole_min_atoms) <= int(self.hole_max_atoms):
            raise BuilderError("Custom hole sizes need 1 <= minimum <= maximum (atoms)")
        if int(self.pore_count) < 0 or int(self.rupture_count) < 0:
            raise BuilderError("Pore and rupture counts must be zero or positive")
        if int(self.pore_count) > 0:
            _positive("Pore diameter", self.pore_diameter_nm)
        if int(self.rupture_count) > 0 and not 1 <= int(self.rupture_length) <= 50:
            raise BuilderError("Rupture length must be 1 to 50 broken bonds")
        _check_choice("Oxidation profile", self.oxidation_profile, OXIDATION_PROFILES)
        _check_choice("Oxidation pattern", self.oxidation_pattern, OXIDATION_PATTERNS)
        _check_choice("Oxidized faces", self.oxidation_faces, OXIDATION_FACES)
        if self.oxidation_profile == "fresh":
            parse_ratio(self.hydroxyl_epoxy_ratio)
        if not 0 <= float(self.paired_hydroxyl_fraction) <= 1:
            raise BuilderError("Paired-hydroxyl fraction must be between 0 and 1")
        if self.oxidation_profile == "custom":
            unknown = set(self.custom_fractions) - set(OXYGEN_GROUPS)
            if unknown:
                raise BuilderError(f"Unknown oxygen groups: {', '.join(sorted(unknown))}")
            values = [float(v) for v in self.custom_fractions.values()]
            if any(v < 0 or not np.isfinite(v) for v in values) or sum(values) <= 0:
                raise BuilderError("Custom oxygen-group fractions must be non-negative and not all zero")
        if self.oxidation_pattern == "domains":
            if int(self.domain_count) < 1:
                raise BuilderError("Oxidized domains: count must be at least 1")
            _positive("Domain diameter", self.domain_diameter_nm)
        if self.oxidation_pattern == "tiles":
            grid = tuple(int(i) for i in self.tile_grid)
            if len(grid) != 2 or min(grid) < 1 or max(grid) > 20:
                raise BuilderError("Tile grid must be two integers between 1 and 20")
            if float(self.tile_oxygen_min) > float(self.tile_oxygen_max):
                raise BuilderError("Tile oxygen minimum must not exceed the maximum")
        if (self.kind == "graphene" and int(self.layers) > 1 and self.oxidize_all_layers
                and self.oxygen_percent > 0 and self.layer_spacing_A < 6.0):
            raise BuilderError("Oxidizing every layer of a stack needs room for the groups between "
                               "layers: set the layer spacing to at least 6 A (hydrated GO stacks "
                               "are about 7-9 A) or oxidize the outer faces only")

    def summary(self) -> str:
        if self.kind == "graphene":
            ox = self.oxygen_percent > 0 or self.oxidation_pattern == "tiles"
            what = "Graphene oxide" if ox else "Graphene"
            layers = f", {self.layers} layers" if self.layers > 1 else ""
            mode = "periodic sheet" if self.periodic else "flake"
            return f"{what} {mode} {self.length_nm:g} x {self.width_nm:g} nm{layers}"
        if self.kind == "nanotube":
            walls = "SWCNT" if self.cnt_walls == 1 else f"{self.cnt_walls}-wall MWCNT"
            ox = " (oxidized)" if self.oxygen_percent > 0 else ""
            mode = "periodic" if self.periodic else f"{self.cnt_length_nm:g} nm"
            return f"{walls} ({self.cnt_n},{self.cnt_m}), {mode}{ox}"
        ox = " (oxidized)" if self.oxygen_percent > 0 else ""
        if self.sphere_type == "fullerene":
            return f"Fullerene C{60 * self.fullerene_order ** 2}{ox}"
        shells = "@".join(f"C{60 * k * k}" for k in range(1, self.onion_shells + 1))
        return f"Carbon onion {shells}{ox}"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["tile_grid"] = [int(i) for i in self.tile_grid]
        d["custom_fractions"] = {k: float(v) for k, v in self.custom_fractions.items()}
        return d


@dataclass
class LigandConfig:
    smiles: str = "SCCO"
    anchor_atom_index: int = 0
    binding_template: str = "Au_thiolate"
    distance_A: Optional[float] = None  # None = template default
    surface_species: Optional[list] = None  # None = template default / all exposed
    conformers: int = 1
    orientation: str = "random_azimuth"
    seed: int = 42
    conformer_mode: str = "etkdg"  # "etkdg" (RDKit as embedded) or "extended" (all-trans chains)
    library_id: Optional[str] = None  # ligand-library entry the inputs came from, if any

    def validate(self) -> None:
        if not str(self.smiles).strip():
            raise BuilderError("Enter a ligand SMILES string")
        _check_choice("Attachment mode", self.binding_template, BINDING_TEMPLATES)
        _check_choice("Orientation", self.orientation, ORIENTATIONS)
        _check_choice("Conformer mode", self.conformer_mode, CONFORMER_MODES)
        if self.distance_A is not None:
            _positive("Anchor distance", self.distance_A)
        if not 1 <= int(self.conformers) <= 20:
            raise BuilderError("Conformer count must be between 1 and 20")
        if int(self.anchor_atom_index) < 0:
            raise BuilderError("Anchor atom index must be non-negative")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class CoverageConfig:
    mode: str = "count"
    target: int = 12
    density_nm2: float = 0.5
    region: str = "entire_surface"
    distribution: str = "distributed_random"
    min_separation_A: float = 0.0
    clash_alpha: float = 0.7
    max_attempts: int = 50000  # trials (sites x azimuths x conformers); dense large surfaces need many
    seed: int = 42

    def validate(self) -> None:
        _check_choice("Coverage mode", self.mode, COVERAGE_MODES)
        _check_choice("Placement region", self.region, REGIONS)
        _check_choice("Distribution", self.distribution, DISTRIBUTIONS)
        if self.mode == "count" and int(self.target) < 0:
            raise BuilderError("Number of ligands must be zero or positive")
        if self.mode == "density" and (not np.isfinite(self.density_nm2) or self.density_nm2 < 0):
            raise BuilderError("Density must be zero or positive")
        if self.min_separation_A < 0:
            raise BuilderError("Minimum site separation must be zero or positive")
        if not 0.3 <= float(self.clash_alpha) <= 1.2:
            raise BuilderError("Clash tolerance must be between 0.3 and 1.2")
        if int(self.max_attempts) < 1:
            raise BuilderError("Maximum placement attempts must be at least 1")

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("density_nm2" if self.mode == "count" else "target")
        return d


@dataclass
class AttachmentConfig:
    """What ``find_surface_sites`` needs; assembled from the ligand and coverage configs."""

    region: str = "entire_surface"
    surface_species: list = field(default_factory=list)
    distance_A: float = 2.4
    site_type: str = "on_top"
    template: str = ""  # attachment template key (carbon substrates choose sites by template)


@dataclass
class ExportConfig:
    targets: list = field(default_factory=lambda: list(EXPORT_TARGETS))
    fidelity: str = "structure_topology"
    include_hydrogens: bool = True
    freeze_substrate: bool = False
    charges: str = "unassigned"
    residue_scheme: str = "groups"  # MOL2 for carbon substrates: per-group residues or one "GRO"
    save_as: str = ""  # file name for the download (empty: a Buddy name is made up)

    def validate(self) -> None:
        if len(str(self.save_as or "")) > 80:
            raise BuilderError("Save-as names are limited to 80 characters")
        for target in self.targets:
            _check_choice("Export target", target, EXPORT_TARGETS)
        _check_choice("Output fidelity", self.fidelity, FIDELITIES)
        _check_choice("Residue scheme", self.residue_scheme, RESIDUE_SCHEMES)
        if self.fidelity != "structure_topology":
            raise BuilderError("Simulation-ready export needs a validated force-field template; "
                               "none is available in this version")
        if not self.include_hydrogens:
            raise BuilderError("This version always exports explicit hydrogens")
        if self.charges != "unassigned":
            raise BuilderError("This version has no validated charge model; charges stay unassigned")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class RunConfig:
    material: MaterialConfig = field(default_factory=MaterialConfig)
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    ligand: LigandConfig = field(default_factory=LigandConfig)
    coverage: CoverageConfig = field(default_factory=CoverageConfig)
    export: ExportConfig = field(default_factory=ExportConfig)
    carbon: CarbonConfig = field(default_factory=CarbonConfig)
    format_version: str = FORMAT_VERSION

    def is_carbon(self) -> bool:
        return self.material.source == "carbon"

    def validate(self) -> None:
        self.material.validate()
        self.geometry.validate()
        self.coverage.validate()
        self.export.validate()
        if self.is_carbon():
            self.carbon.validate()
        if self.wants_ligands():
            self.ligand.validate()

    def wants_ligands(self) -> bool:
        c = self.coverage
        return (c.mode == "count" and int(c.target) > 0) or (c.mode == "density" and c.density_nm2 > 0)

    def to_dict(self) -> dict:
        d = {
            "format_version": self.format_version,
            "material": self.material.to_dict(),
            "geometry": self.geometry.to_dict(),
            "ligand": self.ligand.to_dict(),
            "coverage": self.coverage.to_dict(),
            "export": self.export.to_dict(),
        }
        if self.is_carbon():
            # Shape comes from the carbon settings; only the box padding and size cap apply.
            d["geometry"] = {"vacuum_A": self.geometry.vacuum_A, "max_atoms": self.geometry.max_atoms,
                             "box_shape": self.geometry.box_shape}
            if self.geometry.box_lengths_nm is not None:
                d["geometry"]["box_lengths_nm"] = [float(x) for x in self.geometry.box_lengths_nm]
            d["carbon"] = self.carbon.to_dict()
        return d

    @classmethod
    def from_dict(cls, d: dict, uploaded_bytes: Optional[bytes] = None) -> "RunConfig":
        def build(kind, data):
            data = dict(data or {})
            data.pop("pbc", None)  # derived from shape
            known = {f for f in kind.__dataclass_fields__}
            unknown = set(data) - known
            if unknown:
                raise BuilderError(f"Unknown {kind.__name__} keys: {', '.join(sorted(unknown))}")
            return kind(**data)

        material = build(MaterialConfig, d.get("material"))
        material.uploaded_bytes = uploaded_bytes
        if material.source == "carbon":
            material.preset_id = None
        geometry = build(GeometryConfig, d.get("geometry"))
        geometry.miller = tuple(geometry.miller)
        geometry.rod_axis = tuple(geometry.rod_axis)
        if geometry.box_lengths_nm is not None:
            geometry.box_lengths_nm = tuple(float(x) for x in geometry.box_lengths_nm)
        carbon = build(CarbonConfig, d.get("carbon"))
        carbon.tile_grid = tuple(int(i) for i in carbon.tile_grid)
        carbon.custom_fractions = dict(carbon.custom_fractions)
        return cls(
            material=material,
            geometry=geometry,
            ligand=build(LigandConfig, d.get("ligand")),
            coverage=build(CoverageConfig, d.get("coverage")),
            export=build(ExportConfig, d.get("export")),
            carbon=carbon,
            format_version=d.get("format_version", FORMAT_VERSION),
        )


# --------------------------------------------------------------------------
# Models passed between pipeline stages
# --------------------------------------------------------------------------


@dataclass
class CrystalModel:
    """Normalized bulk crystal; identical contract for presets and uploads (spec 5.2)."""

    atoms: Atoms
    source_type: str
    preset_id: Optional[str] = None
    source_format: Optional[str] = None
    source_filename: Optional[str] = None
    source_sha256: Optional[str] = None
    source_bytes: Optional[bytes] = field(default=None, repr=False)
    lattice_overrides_A: dict = field(default_factory=dict)
    structure_label: str = ""
    space_group: Optional[str] = None
    prototype: Optional[str] = None
    anchor_species_candidates: list = field(default_factory=list)
    lattice_source: Optional[str] = None
    normalization: str = "none"
    validation_warnings: list = field(default_factory=list)

    @property
    def cell_A(self) -> np.ndarray:
        return np.array(self.atoms.cell)

    @property
    def fractional_positions(self) -> np.ndarray:
        return self.atoms.get_scaled_positions(wrap=True)

    @property
    def species(self) -> list:
        return sorted(set(self.atoms.get_chemical_symbols()))

    @property
    def composition(self) -> dict:
        symbols = self.atoms.get_chemical_symbols()
        return {s: symbols.count(s) for s in self.species}

    def provenance(self) -> dict:
        return {
            "source_type": self.source_type,
            "preset_id": self.preset_id,
            "source_format": self.source_format,
            "source_filename": self.source_filename,
            "source_sha256": self.source_sha256,
            "lattice_overrides_A": self.lattice_overrides_A,
            "lattice_source": self.lattice_source,
            "structure_label": self.structure_label,
            "space_group": self.space_group,
            "normalization": self.normalization,
            "cell_A": self.cell_A.round(6).tolist(),
            "composition_per_cell": self.composition,
        }


_FIXTURE_FIELDS = ("source_type", "preset_id", "source_format", "source_filename", "source_sha256",
                   "lattice_overrides_A", "structure_label", "space_group", "prototype",
                   "anchor_species_candidates", "lattice_source", "normalization",
                   "validation_warnings")


def save_crystal_fixture(crystal: CrystalModel, path) -> None:
    """Write ``<path>.extxyz`` plus a ``<path>.json`` sidecar for team handoffs (spec 5.2)."""
    import json
    from pathlib import Path

    from ase.io import write

    path = Path(path)
    write(path.with_suffix(".extxyz"), crystal.atoms, format="extxyz")
    sidecar = {k: getattr(crystal, k) for k in _FIXTURE_FIELDS}
    path.with_suffix(".json").write_text(json.dumps(sidecar, indent=2), encoding="utf-8")


def load_crystal_fixture(path) -> CrystalModel:
    import json
    from pathlib import Path

    from ase.io import read

    path = Path(path)
    atoms = read(path.with_suffix(".extxyz"), format="extxyz")
    sidecar = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    return CrystalModel(atoms=atoms, **sidecar)


@dataclass
class Termination:
    index: int
    top_elements: list
    bottom_elements: list
    label: str
    plane_offset: int = 0  # which atomic plane of the repeat unit ends up on top


@dataclass
class FaceInfo:
    atom_indices: np.ndarray
    exposed_elements: list
    termination_label: str
    area_nm2: float


@dataclass
class SubstrateModel:
    atoms: Atoms
    crystal: CrystalModel
    shape: str
    pbc: tuple
    faces: dict
    composition: dict
    center_A: Optional[np.ndarray] = None
    radius_A: Optional[float] = None
    miller: Optional[tuple] = None
    surface_normal: Optional[np.ndarray] = None  # output frame (+z for slabs)
    crystal_frame_normal: Optional[np.ndarray] = None  # h*b1 + k*b2 + l*b3, normalized
    interplanar_spacing_A: Optional[float] = None
    layers: Optional[int] = None
    thickness_A: Optional[float] = None
    termination_index: Optional[int] = None
    terminations: list = field(default_factory=list)
    vacuum_A: float = 30.0
    warnings: list = field(default_factory=list)
    # Nanorods (output frame: axis along +z through center_A)
    rod_axis: Optional[np.ndarray] = None
    rod_half_length_A: Optional[float] = None  # half length of the cylindrical section
    rod_caps: Optional[str] = None
    rod_crystal_axis: Optional[tuple] = None
    # Covalent substrates (carbon nanomaterials). Metals and crystals leave these empty.
    bonds: list = field(default_factory=list)  # (i, j, order) inside the substrate
    labels: Optional[list] = None  # structural type label per atom (default "<El>_sub")
    residues: Optional[list] = None  # (residue name, residue number) per atom
    normals: Optional[np.ndarray] = None  # outward unit normal per atom
    carbon_sites: dict = field(default_factory=dict)  # attachment-site candidates by kind
    carbon_info: dict = field(default_factory=dict)  # composition, defects, metrics
    layer_ids: Optional[np.ndarray] = None  # carbon: layer / wall / shell index per atom
    # Simulation box requested by the user (applied when the model is assembled)
    box_shape: str = "auto"
    box_min_A: Optional[tuple] = None  # user box (x, y, z) in A; 0 = automatic (fit + vacuum)
    channel_info: dict = field(default_factory=dict)  # crystal channels: type, width, walls

    @property
    def is_covalent(self) -> bool:
        return bool(self.bonds)


@dataclass
class LigandModel:
    smiles: str
    elements: list
    positions_A: np.ndarray  # anchor at origin
    bonds: list  # (i, j, order) with zero-based indices into this molecule
    formal_charges: list
    anchor_index: int
    index_map: dict  # original (explicit-H) atom index -> bound-molecule index
    removed_atoms: list
    outward_axis: np.ndarray
    binding_template: str
    distance_A: float
    surface_species: list
    formula: str = ""
    orientation: str = "random_azimuth"
    conformers_A: list = field(default_factory=list)  # all conformers, anchor at origin
    optimization: str = ""
    warnings: list = field(default_factory=list)
    outward_axes: list = field(default_factory=list)  # per conformer (covalent: anchor - removed H)

    @property
    def n_atoms(self) -> int:
        return len(self.elements)


@dataclass
class SiteCatalog:
    positions_A: np.ndarray
    surface_atom_index: np.ndarray
    normals: np.ndarray
    face: list
    element: list
    site_type: list
    area_nm2_by_face: dict
    warnings: list = field(default_factory=list)
    # Covalent carbon chemistry: per-site substrate edits applied when a ligand is placed
    # (removed/added/moved atoms, relabels) and the substrate atoms each site consumes.
    edits: Optional[list] = None
    claims: Optional[list] = None

    def __len__(self) -> int:
        return len(self.surface_atom_index)


@dataclass
class AssemblyModel:
    """Canonical atomistic model consumed by every exporter (spec 5.4)."""

    positions_A: np.ndarray
    elements: list
    cell_A: np.ndarray
    pbc: tuple
    atom_types: Optional[list]  # force-field types; None while unassigned
    type_labels: list  # structural labels (element + role), NOT force-field types
    chemical_role: list
    masses_amu: np.ndarray
    charges_e: Optional[np.ndarray]
    molecule_id: np.ndarray
    component: list
    bond_list: list  # (i, j, bond_type)
    angle_list: list
    dihedral_list: list
    improper_list: list
    site_to_ligand: list
    surface_area_nm2: float
    target_ligands: int
    achieved_ligands: int
    requested_density_nm2: Optional[float]
    warnings: list
    provenance: dict
    parameter_status: str = "unassigned"
    min_nonbonded_distance_A: Optional[float] = None
    bond_orders: Optional[list] = None  # parallel to bond_list (1, 1.5 aromatic, 2, 3)
    residue_names: Optional[list] = None  # per atom, substrate residues of covalent substrates
    residue_ids: Optional[list] = None
    atom_names: Optional[list] = None
    suggested_types: dict = field(default_factory=dict)  # e.g. {"gaff2": [...], "opls_hiergo": [...]}
    placements: list = field(default_factory=list)  # ligand placements (kept for viewer edits)
    substrate_index: Optional[list] = None  # per atom: index in SubstrateModel.atoms, -1 for ligands

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def achieved_density_nm2(self) -> Optional[float]:
        if self.surface_area_nm2 <= 0:
            return None
        return self.achieved_ligands / self.surface_area_nm2

    @property
    def n_substrate(self) -> int:
        return sum(1 for c in self.component if c == "substrate")


@dataclass
class ValidationReport:
    status: str
    errors: list
    warnings: list
    target_ligands: int
    placed_ligands: int
    requested_density_nm2: Optional[float]
    achieved_density_nm2: Optional[float]
    min_nonbonded_distance_A: Optional[float]
    parameter_status: str
    counts: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not self.errors

    def to_dict(self) -> dict:
        d = asdict(self)
        d["partial_charges_unassigned"] = self.parameter_status != "complete_validated"
        return d


@dataclass
class ExportReport:
    path: str
    format: str
    sha256: str
    n_atoms: int
    n_bonds: int
    round_trip_passed: Optional[bool]
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)
