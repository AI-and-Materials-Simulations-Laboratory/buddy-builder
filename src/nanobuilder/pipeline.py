"""One call per GUI action: validate inputs, run the full build, export a package.

The GUI (and scripts or tests) use only these functions, so the scientific
modules never depend on a notebook session.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Optional

from .chemistry import (TEMPLATES, check_template, resolve_surface_species, template_keys,
                        thiolate_unavailable_reason)
from .crystal import build_crystal
from .ligands import build_ligand, check_anchor, parse_smiles
from .placement import bare_assembly, place_ligands, target_count
from .reporting import Package, build_package
from .schema import (AssemblyModel, AttachmentConfig, BuilderError, CrystalModel, LigandModel,
                     RunConfig, SiteCatalog, SubstrateModel, ValidationReport)
from .shapes import build_geometry, list_terminations
from .sites import find_surface_sites, region_faces, selected_area_nm2
from .validation import validate_assembly

STAGES = ("substrate", "ligand", "sites", "packing", "validation")


@dataclass
class BuildResult:
    config: RunConfig
    crystal: CrystalModel
    substrate: SubstrateModel
    ligand: Optional[LigandModel]
    sites: Optional[SiteCatalog]
    assembly: AssemblyModel
    report: ValidationReport
    timings: dict = field(default_factory=dict)


def build_substrate(config: RunConfig):
    """(CrystalModel, SubstrateModel) for any material source, carbon included."""
    g = config.geometry
    if config.is_carbon():
        from .carbon import build_carbon_substrate
        carbon = config.carbon
        from dataclasses import replace
        box = g.box_lengths_nm if g.box_lengths_nm and any(float(x) > 0 for x in g.box_lengths_nm) \
            else None
        if box is not None and carbon.periodic and carbon.kind == "graphene":
            # periodic sheet: the box sets the sheet (the lattice rounds to whole repeats)
            carbon = replace(carbon, length_nm=float(box[0]) or carbon.length_nm,
                             width_nm=float(box[1]) or carbon.width_nm)
        elif box is not None and carbon.periodic and carbon.kind == "nanotube" and float(box[2]) > 0:
            carbon = replace(carbon, cnt_length_nm=float(box[2]))
        elif g.box_shape == "cube" and carbon.kind == "graphene" and carbon.periodic:
            side = max(carbon.length_nm, carbon.width_nm)
            carbon = replace(carbon, length_nm=side, width_nm=side)
        crystal, substrate = build_carbon_substrate(carbon, g)
    else:
        crystal = build_crystal(config.material)
        substrate = build_geometry(crystal, g)
    substrate.box_shape = g.box_shape
    if g.box_lengths_nm is not None and any(float(x) > 0 for x in g.box_lengths_nm):
        substrate.box_min_A = tuple(10.0 * float(x) for x in g.box_lengths_nm)
    return crystal, substrate


def material_crystal(config: RunConfig) -> CrystalModel:
    """The crystal (or the carbon provenance stand-in) without building the substrate."""
    if config.is_carbon():
        from .carbon import carbon_crystal
        return carbon_crystal(config.carbon, None)
    return build_crystal(config.material)


def substrate_shape(config: RunConfig) -> str:
    if config.is_carbon():
        return {"graphene": "sheet", "nanotube": "nanotube",
                "carbon_sphere": "carbon_sphere"}[config.carbon.kind]
    return config.geometry.shape


def validate_inputs(config: RunConfig) -> dict:
    """Fast checks without building the full structure (the GUI's Validate button)."""
    errors, warnings, info = [], [], {}
    try:
        config.validate()
        crystal = material_crystal(config)
        info["material"] = crystal.structure_label
        if config.is_carbon():
            from .carbon import estimate_atoms
            est = estimate_atoms(config.carbon)
            info["estimated_atoms"] = est
            if est > config.geometry.max_atoms:
                errors.append(f"About {est} atoms; the limit is {config.geometry.max_atoms}")
        else:
            info["atoms_per_cell"] = len(crystal.atoms)
        warnings += crystal.validation_warnings
        if config.geometry.shape in ("slab", "channel") and not config.is_carbon():
            terms = list_terminations(crystal, config.geometry.miller, 2)
            info["terminations"] = [t.label for t in terms]
            if config.geometry.termination_index >= len(terms):
                errors.append(f"Termination {config.geometry.termination_index} does not exist")
        region_faces_check(config)
        if config.wants_ligands():
            template = check_template(config.ligand.binding_template, crystal)
            mol = parse_smiles(config.ligand.smiles)
            check_anchor(mol, int(config.ligand.anchor_atom_index), template)
            info["surface_species"] = resolve_surface_species(
                template, config.ligand.surface_species, crystal)
    except BuilderError as exc:
        errors.append(str(exc))
    return {"ok": not errors, "errors": errors, "warnings": warnings, "info": info}


def region_faces_check(config: RunConfig) -> None:
    region = config.coverage.region
    if config.is_carbon():
        from .carbon import region_kinds
        region_kinds(substrate_shape(config), region)
        return
    if config.geometry.shape == "sphere" and region != "entire_surface":
        raise BuilderError("A sphere has one surface region: choose 'Entire surface'")
    if config.geometry.shape == "rod" and region not in ("entire_surface", "side", "tips"):
        raise BuilderError("A rod has the regions 'Entire surface', 'Side' and 'Tips'")
    if config.geometry.shape == "channel" and region not in ("inner", "top", "bottom", "both",
                                                             "entire_surface"):
        raise BuilderError("A channel has the regions 'Channel walls', 'Top', 'Bottom', 'Both "
                           "outer faces' and 'Entire surface'")
    if config.geometry.shape == "slab" and region not in ("top", "bottom", "both", "entire_surface"):
        raise BuilderError("A slab has the regions 'Top', 'Bottom' and 'Both faces'")


def run_pipeline(config: RunConfig, progress: Optional[Callable] = None) -> BuildResult:
    """substrate -> ligand -> sites -> packing -> validation."""
    def step(name):
        if progress:
            progress(name, STAGES.index(name) / len(STAGES))

    timings = {}
    t0 = time.perf_counter()
    config.validate()
    region_faces_check(config)
    step("substrate")
    crystal, substrate = build_substrate(config)
    timings["substrate_s"] = time.perf_counter() - t0

    ligand = sites = None
    region = config.coverage.region
    area = selected_area_nm2(substrate, region)
    n_target = target_count(config.coverage, area) if config.wants_ligands() else 0
    if n_target == 0:
        warnings = list(substrate.warnings)
        if config.wants_ligands():
            warnings.append(f"Density {config.coverage.density_nm2} ligands/nm^2 on "
                            f"{area:.2f} nm^2 rounds to 0 ligands")
        step("packing")
        assembly = bare_assembly(substrate, warnings)
        assembly.surface_area_nm2 = area
    else:
        step("ligand")
        t = time.perf_counter()
        template = check_template(config.ligand.binding_template, crystal)
        species = resolve_surface_species(template, config.ligand.surface_species, crystal)
        config.ligand.surface_species = species
        ligand = build_ligand(config.ligand)
        timings["ligand_s"] = time.perf_counter() - t
        step("sites")
        sites = find_surface_sites(substrate, AttachmentConfig(
            region=region, surface_species=species, distance_A=ligand.distance_A,
            template=template.key))
        step("packing")
        t = time.perf_counter()
        warnings = list(substrate.warnings) + list(ligand.warnings)
        if crystal.lattice_overrides_A and template.key == "Au_thiolate":
            warnings.append("Au-thiolate on a user-modified Au lattice constant")
        assembly = place_ligands(substrate, ligand, sites, config.coverage, area, warnings)
        assembly.provenance["seeds"]["ligand_conformer"] = config.ligand.seed
        timings["packing_s"] = time.perf_counter() - t
    from .atomtypes import formal_charges_for, suggest_types
    assembly.suggested_types = suggest_types(assembly, formal_charges_for(assembly, ligand))
    step("validation")
    report = validate_assembly(assembly)
    timings["total_s"] = time.perf_counter() - t0
    if progress:
        progress("done", 1.0)
    return BuildResult(config=config, crystal=crystal, substrate=substrate, ligand=ligand,
                       sites=sites, assembly=assembly, report=report, timings=timings)


def export_package(result: BuildResult, out_dir=None) -> Package:
    return build_package(result.config, result.assembly, result.report, result.crystal.atoms,
                         result.crystal.source_bytes, out_dir)


def attachment_options(crystal: CrystalModel) -> dict:
    """What the ligand tab should offer for this crystal (or carbon substrate)."""
    keys = template_keys(crystal)
    return {"templates": [(TEMPLATES[k].label, k) for k in keys],
            "reason": thiolate_unavailable_reason(crystal),
            "species": list(crystal.anchor_species_candidates)}


__all__ = ["BuildResult", "validate_inputs", "run_pipeline", "export_package",
           "attachment_options", "region_faces", "build_substrate", "material_crystal"]
