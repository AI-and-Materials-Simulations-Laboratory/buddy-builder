"""ZIP package, manifest and plain-text documentation of an export (spec F1, F6)."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import tempfile
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .cell import recommended_boundaries
from .exporter_generic import export_generic
from .exporter_gsd import export_hoomd_gsd
from .exporter_lammps import export_lammps_data
from .exporter_mol2 import MOL2_README, TLEAP_TEMPLATE, export_mol2
from .exporter_pdb import export_pdb
from .importers import atoms_to_cif_bytes
from .schema import AssemblyModel, BuilderError, RunConfig, ValidationReport
from .topology import parameter_coverage

READINESS_TEXT = {
    "geometry_only": "Atomic coordinates generated; bonding/typing incomplete",
    "topology_ready": "MD-readable coordinates/topology; force-field setup required",
    "partially_parameterized": "Some parameters/charges missing; do not run unmodified",
    "simulation_ready_validated": "Verified force-field profile and engine import/smoke test passed",
}


def readiness(assembly: AssemblyModel) -> str:
    if assembly.parameter_status == "complete_validated":
        return "simulation_ready_validated"
    if assembly.parameter_status == "partial":
        return "partially_parameterized"
    return "topology_ready"


@dataclass
class Package:
    zip_path: Path
    zip_name: str
    stage_dir: Path
    manifest: dict
    reports: list = field(default_factory=list)

    def read_bytes(self) -> bytes:
        return self.zip_path.read_bytes()


# Southern terms of endearment for builds nobody named.
ENDEARMENTS = ["Sugar", "Honey", "Darlin", "Sweet_Pea", "Butterbean", "Peaches", "Sugar_Pie",
               "Honey_Bun", "Biscuit", "Cupcake", "Doodlebug", "Sweet_Tea", "Buttercup",
               "Pumpkin", "Pecan_Pie", "Cornbread", "Magnolia", "Hush_Puppy", "Dumplin",
               "Sugar_Britches", "Puddin", "Lamb_Chop", "Honeychile", "Sweetie_Pie"]


def safe_name(text: str) -> str:
    """A file-system-safe version of a user's save-as name ("" if nothing usable is left)."""
    name = re.sub(r"[^A-Za-z0-9_-]+", "_", str(text or "").strip())
    return re.sub(r"_+", "_", name).strip("_.-")[:80]


def _tokens(config: RunConfig, assembly: AssemblyModel) -> tuple:
    """(material, shape) short tokens of a build, e.g. ("Au", "sphere")."""
    m, g = config.material, config.geometry
    if m.source == "carbon":
        c = config.carbon
        if c.kind == "graphene":
            material = "GO" if (c.oxygen_percent > 0 or c.oxidation_pattern == "tiles") else "graphene"
            return material, "sheet" if c.periodic else "flake"
        if c.kind == "nanotube":
            return f"CNT{c.cnt_n}-{c.cnt_m}", "tube"
        return (f"C{60 * c.fullerene_order ** 2}" if c.sphere_type == "fullerene"
                else "onion"), "ball"
    if m.source == "preset":
        material = (m.preset_id or "crystal").split("_")[0]
    elif m.source == "materials_project":
        try:
            from .mp_library import entry
            material = (entry(m.mp_id) or {}).get("formula") or str(m.mp_id)
        except Exception:
            material = str(m.mp_id)
    else:
        material = Path(m.filename or "crystal").stem
    shape = {"sphere": "sphere", "rod": "rod", "slab": "slab"}.get(
        g.shape, "pores" if g.channel_type == "pore" else "channels")
    return material, shape


def buddy_name(config: RunConfig, assembly: AssemblyModel) -> str:
    """A made-up Buddy name for an unnamed build, the same for the same build."""
    material, shape = _tokens(config, assembly)
    key = hashlib.sha256(json.dumps(config.to_dict(), sort_keys=True, default=str).encode()
                         + str(assembly.n_atoms).encode()).digest()
    word = ENDEARMENTS[key[0] % len(ENDEARMENTS)]
    return safe_name(f"{word}_Buddy_{material}_{shape}") or "Buddy_build"


def package_name(config: RunConfig, assembly: AssemblyModel) -> str:
    """ZIP name: the user's save-as name, or a Buddy name such as Sweet_Pea_Buddy_Au_sphere."""
    name = safe_name(config.export.save_as) or buddy_name(config, assembly)
    return name + ".zip"


def readme_text(config: RunConfig, assembly: AssemblyModel, report: ValidationReport,
                files: list) -> str:
    level = readiness(assembly)
    prov = assembly.provenance
    lig = prov.get("ligand")
    lines = [
        f"nanobuilder {__version__} export",
        "=" * 40,
        "",
        f"Readiness: {level} - {READINESS_TEXT[level]}",
        "Level B (MD structure/topology): readable by HOOMD-blue, LAMMPS, AMBER and GROMACS",
        "(via the PDB or MOL2 file), but NOT",
        "simulation-ready. No force field, partial charges or substrate-ligand",
        "interaction parameters were generated.",
        "",
        f"Material: {prov['material']['structure_label']}",
        f"Shape: {assembly.provenance['geometry']['shape']}",
        f"Atoms: {assembly.n_atoms} ({assembly.n_substrate} substrate)",
        f"Ligands: {assembly.achieved_ligands} placed of {assembly.target_ligands} requested",
        f"Surface area used for coverage: {assembly.surface_area_nm2:.3f} nm^2",
        "",
        "session.py: paste into the Buddy Builder Python console (or use Load session) to",
        "restore every selection and rebuild this model, hand edits included.",
        "certificate_of_build.pdf: certificate plus a report of the structure, crystal,",
        "atoms, ligand, coverage, topology and validation, with pictures.",
    ]
    if assembly.achieved_density_nm2 is not None and assembly.target_ligands:
        lines.append(f"Achieved density: {assembly.achieved_density_nm2:.3f} ligands/nm^2"
                     + (f" (requested {assembly.requested_density_nm2:.3f})"
                        if assembly.requested_density_nm2 is not None else ""))
    if lig:
        lines.append(f"Ligand: {lig['smiles']} via {lig['binding_template']} "
                     f"(anchor distance {lig['distance_A']} A)")
    carbon = prov.get("carbon")
    if carbon:
        comp = carbon.get("composition_after_attachment") or carbon.get("composition", {})
        lines += ["", "Carbon nanomaterial:", f"  {carbon.get('summary', '')}",
                  f"  Composition: C {comp.get('C')}, O {comp.get('O')}, H {comp.get('H')}",
                  f"  Oxygen content 100 N_O/N_C: {comp.get('oxygen_percent_O_per_C')} %"
                  + (f"; C:O = {comp.get('C_to_O_ratio')}:1" if comp.get('C_to_O_ratio') else ""),
                  f"  Modified (sp3) framework carbons: {comp.get('modified_carbon_percent')} %",
                  f"  Vacancies relative to the pristine lattice: "
                  f"{carbon.get('composition', {}).get('vacancy_percent_actual')} %"]
        ox = carbon.get("oxidation") or {}
        if ox.get("achieved_hydroxyl_to_epoxy"):
            lines.append(f"  Basal hydroxyl : epoxide oxygen = {ox['achieved_hydroxyl_to_epoxy']}:1"
                         + (" (dense sheets convert epoxides to diols above about 40 % O)"
                            if ox.get("profile") == "fresh"
                            and float(ox.get("requested_percent") or 0) > 40 else ""))
        if carbon.get("groups"):
            lines.append("  Oxygen groups: " + ", ".join(f"{k} {v}" for k, v in carbon["groups"].items()))
        lines += [f"  {line}" for line in defect_lines(carbon.get("defects") or {})]
    edits = prov.get("manual_edits") or []
    if edits:
        lines += ["", "MANUAL EDITS (made by hand in the viewer or console after Generate; the",
                  "model is not reproducible from config.json alone):"]
        for e in edits:
            if e.get("action") == "remove_ligands":
                lines.append(f"  - removed ligand(s) {', '.join(str(i) for i in e['ligands'])}")
            elif e.get("action") == "add_ligand":
                lines.append(f"  - added ligand {e.get('ligand')} on surface atom {e.get('surface_atom')}")
            else:
                lines.append(f"  - {e}")
    lines += [
        f"Validation status: {report.status}",
        "",
        "Units: Angstrom, amu. Internal indices are zero-based; LAMMPS IDs are one-based.",
        f"Recommended LAMMPS boundary: {recommended_boundaries(assembly.pbc)['lammps_boundary']}",
        "",
        "Files:",
    ]
    lines += [f"  {f}" for f in files]
    lines += ["", "Warnings:"] + [f"  - {w}" for w in report.warnings + assembly.warnings]
    lines += [
        "",
        "Scientific limitations:",
        "  - Ligand placement is a geometric starting configuration, not an energy minimum.",
        "  - RDKit/MMFF optimization covers the isolated ligand only, not the interface.",
        "  - Type labels such as Au_sub, S_anchor, C_sp2 or O_epoxy are structural labels, not",
        "    force-field types; replace them with real types when parameterizing. Suggested",
        "    GAFF2 (and, for carbon substrates, HierGO OPLS-AA) names are in topology.json",
        "    and the MOL2 file.",
        "  - Partial charges are unassigned; LAMMPS q = 0 values are syntax placeholders.",
    ]
    return "\n".join(lines) + "\n"


def defect_lines(defects: dict) -> list:
    """Readable one-liners for the carbon defect summary."""
    out = []
    holes = defects.get("holes")
    if holes:
        out.append(f"Holes: {holes['holes']} ({holes['shapes'].replace('_', ' ')}, "
                   f"{holes['removed_atoms']} atoms removed, {holes['removed_by_cleaning']} of them "
                   "while cleaning edges)")
    pores = defects.get("pores")
    if pores:
        out.append(f"Pores: {pores['pores']} x {pores['diameter_nm']} nm ({pores['removed_atoms']} "
                   "atoms removed)")
    for key, name in (("stone_wales", "Stone-Wales defects"), ("rotations_24_atom",
                                                              "24-atom rotations")):
        d = defects.get(key)
        if d:
            out.append(f"{name}: {d['placed']} of {d['requested']}")
    rup = defects.get("ruptures")
    if rup:
        out.append(f"Ruptures: {rup['ruptures']} ({rup['ether_oxygens']} bridging ether O)")
    extra = defects.get("crowded_edge_atoms_removed")
    if extra:
        out.append(f"Edge atoms removed to fit edge hydrogens: {extra}")
    return out


HOOMD_README = """HOOMD-blue files
================
structure.gsd        One frame written with the GSD Python API: box, particle types
                     (structural labels), type IDs, positions, masses, bonds, angles
                     and dihedrals. Particle charges are not set. The per-atom
                     molecule ID is stored under log/particles/molecule_id.
run_hoomd_template.py  Skeleton showing where the force field must be defined.

HOOMD boxes are always periodic. Nonperiodic directions were padded with vacuum
(see manifest.json); choose cutoffs smaller than that vacuum. Units are
application-defined in HOOMD; these numbers are Angstrom and amu.
"""

PDB_README = """PDB file (for AMBER, GROMACS and visualization)
==============================================
structure.pdb       CRYST1 cell, HETATM records and CONECT bonds. Substrate atoms are
                    chain A: metals and crystals use one residue per atom named after the
                    element (AU, NA, ...); carbon nanomaterials use GRA per framework carbon
                    (with its edge H) and OH, EPX, ETH, CBO, CBX per oxygen group. Each
                    ligand is one residue LIG in chain B with atom names S1, C1, C2, H1 ...
                    Ligands are written whole, never split across the periodic boundary.

AMBER: load into tleap after preparing ligand parameters (e.g. antechamber/GAFF2)
       and a substrate model; residue and atom names must match your templates.
GROMACS: gmx editconf/pdb2gmx read the coordinates and box; supply a force field
       and topology (.top/.itp) for the substrate and ligand.
PDB limits atom serials to 99999 and residue numbers to 9999; larger systems
wrap the numbers (and omit CONECT above 99999 atoms). Use LAMMPS or GSD then.
"""

LAMMPS_README = """LAMMPS files
============
structure.data      read_data input, atom_style full: Masses, Atoms (q = 0
                    placeholders, coordinates inside the box plus ix iy iz image
                    flags so ligands stay whole), Bonds, Angles, Dihedrals. One-based IDs. Type
                    labels are given as comments in the Masses section.
in.lammps.template  Commented template. It is not runnable: pair, bond, angle and
                    dihedral styles/coefficients and real charges must be supplied.

Restricted triclinic box: tilt factors were reduced so |xy| <= lx/2.
"""

HOOMD_TEMPLATE = '''"""HOOMD-blue setup TEMPLATE written by nanobuilder. NOT runnable as-is.

No force field was generated. Define every interaction marked REQUIRED with a
reviewed, compatible model before running.
"""
import hoomd

sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_gsd(filename="structure.gsd")

# REQUIRED: pair potentials for every particle-type pair, including
#           substrate-substrate and substrate-ligand interactions.
# REQUIRED: bond, angle and dihedral potentials for every type in the file.
# REQUIRED: charges and an electrostatics method if the model needs them.
raise NotImplementedError("Define the force field before running this template")
'''


def build_package(config: RunConfig, assembly: AssemblyModel, report: ValidationReport,
                  crystal_atoms, crystal_bytes: bytes | None, out_dir=None) -> Package:
    """Write every selected output into one ZIP; refuses if validation found errors."""
    if not report.ok:
        raise BuilderError("Export blocked by validation errors: " + "; ".join(report.errors))
    out_dir = Path(out_dir or tempfile.mkdtemp(prefix="nanobuilder_"))
    name = package_name(config, assembly)
    stage = out_dir / Path(name).stem  # unpacked copy of the ZIP contents
    shutil.rmtree(stage, ignore_errors=True)
    stage.mkdir(parents=True)
    targets = set(config.export.targets)
    reports = []

    if "hoomd_gsd" in targets:
        reports.append(export_hoomd_gsd(assembly, stage / "hoomd" / "structure.gsd"))
        (stage / "hoomd" / "README_hoomd.txt").write_text(HOOMD_README, encoding="utf-8")
        (stage / "hoomd" / "run_hoomd_template.py").write_text(HOOMD_TEMPLATE, encoding="utf-8")
    if "lammps_data" in targets:
        reports.append(export_lammps_data(assembly, stage / "lammps" / "structure.data",
                                          config.export.freeze_substrate))
        (stage / "lammps" / "README_lammps.txt").write_text(LAMMPS_README, encoding="utf-8")
    if "pdb" in targets:
        reports.append(export_pdb(assembly, stage / "pdb" / "structure.pdb"))
        (stage / "pdb" / "README_pdb.txt").write_text(PDB_README, encoding="utf-8")
    if "mol2" in targets:
        reports.append(export_mol2(assembly, stage / "mol2" / "structure.mol2",
                                   config.export.residue_scheme,
                                   assembly.provenance["material"]["structure_label"]))
        (stage / "mol2" / "README_mol2.txt").write_text(MOL2_README, encoding="utf-8")
        (stage / "mol2" / "tleap_template.in").write_text(TLEAP_TEMPLATE, encoding="utf-8")
    generic = export_generic(assembly, stage)
    if "extxyz" not in targets:
        (stage / "structure.extxyz").unlink()
        generic = [r for r in generic if r.format != "extxyz"]
    reports.extend(generic)

    (stage / "normalized_unit_cell.cif").write_bytes(atoms_to_cif_bytes(crystal_atoms))
    if crystal_bytes is not None:
        ext = ".cif" if config.material.file_format == "cif" else ".vasp"
        (stage / f"input_crystal{ext}").write_bytes(crystal_bytes)
    (stage / "config.json").write_text(json.dumps(config.to_dict(), indent=2), encoding="utf-8")
    from .session import session_text
    (stage / "session.py").write_text(session_text(config, assembly, Path(name).stem),
                                      encoding="utf-8")
    report_dict = report.to_dict()
    report_dict["readiness"] = readiness(assembly)
    report_dict["export_round_trips"] = {Path(r.path).name: r.round_trip_passed for r in reports}
    (stage / "validation_report.json").write_text(json.dumps(report_dict, indent=2),
                                                  encoding="utf-8")
    failed = [r for r in reports if r.round_trip_passed is False]
    if failed:
        raise BuilderError("Export round-trip check failed: " +
                           "; ".join(f"{Path(r.path).name}: {', '.join(r.warnings[:2])}"
                                     for r in failed))

    files = sorted(p.relative_to(stage).as_posix() for p in stage.rglob("*") if p.is_file())
    files += ["certificate_of_build.pdf"]
    try:
        from .certificate import certificate_pdf
        (stage / "certificate_of_build.pdf").write_bytes(certificate_pdf(
            config, assembly, report, crystal_atoms,
            files=sorted(files + ["README.txt", "manifest.json"]), name=Path(name).stem))
    except Exception as exc:  # a picture problem must never block the structure files
        files.remove("certificate_of_build.pdf")
        report.warnings.append(f"Certificate of Build could not be drawn: {exc}")
    (stage / "README.txt").write_text(readme_text(config, assembly, report,
                                                  sorted(files + ["README.txt", "manifest.json"])),
                                      encoding="utf-8")
    files = sorted(p.relative_to(stage).as_posix() for p in stage.rglob("*") if p.is_file())
    manifest = {
        "generator": f"nanobuilder {__version__}",
        "export_time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": {f: hashlib.sha256((stage / f).read_bytes()).hexdigest() for f in files},
        "software_versions": assembly.provenance["versions"],
        "counts": report.counts,
        "material": {
            "preset_id": config.material.preset_id if config.material.source == "preset" else None,
            "uploaded_crystal_sha256": (config.material.sha256
                                        if config.material.source == "uploaded_file" else None),
            "materials_project_id": (config.material.mp_id
                                     if config.material.source == "materials_project" else None),
            "materials_project_cif_sha256": (config.material.sha256
                                             if config.material.source == "materials_project" else None),
            "structure_label": assembly.provenance["material"]["structure_label"],
            "carbon": config.carbon.to_dict() if config.material.source == "carbon" else None,
        },
        "charge_status": "unassigned" if assembly.charges_e is None else "assigned",
        "parameter_status": assembly.parameter_status,
        "readiness": readiness(assembly),
        "parameter_coverage": parameter_coverage(assembly),
        "units": assembly.provenance["units"],
        "boundaries": recommended_boundaries(assembly.pbc),
        "freeze_substrate": config.export.freeze_substrate,
        "assumptions": [
            "Starting geometry only; not energy-minimized",
            "Coverage uses the geometric area stated in validation_report.json",
            "Structural type labels are not force-field types",
        ],
    }
    (stage / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    zip_path = out_dir / name
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(stage.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(stage).as_posix())
    return Package(zip_path=zip_path, zip_name=name, stage_dir=stage, manifest=manifest,
                   reports=reports)
