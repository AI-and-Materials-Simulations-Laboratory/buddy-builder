"""LAMMPS data-file export, ``atom_style full`` (spec F3).

Atom, bond, angle and dihedral IDs and type IDs are one-based; the internal
model is zero-based. Charges are written as q = 0 placeholders only to
satisfy the ``full`` syntax and are flagged as unassigned everywhere.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from . import __version__
from .cell import lammps_box, recommended_boundaries, wrap_with_images
from .schema import AssemblyModel, ExportReport


def _type_table(names) -> dict:
    return {name: k for k, name in enumerate(sorted(set(names)), start=1)}


def type_tables(assembly: AssemblyModel) -> dict:
    return {
        "atom": _type_table(assembly.type_labels),
        "bond": _type_table(b[2] for b in assembly.bond_list),
        "angle": _type_table(a[3] for a in assembly.angle_list),
        "dihedral": _type_table(d[4] for d in assembly.dihedral_list),
    }


def data_file_text(assembly: AssemblyModel) -> str:
    tables = type_tables(assembly)
    box = lammps_box(assembly.cell_A)
    mass_of = {}
    for label, mass in zip(assembly.type_labels, assembly.masses_amu):
        mass_of[label] = mass
    lines = [
        f"LAMMPS data file from nanobuilder {__version__}: structure/topology only, no force "
        "field; q = 0 placeholders (charges unassigned)",
        "",
        f"{assembly.n_atoms} atoms",
        f"{len(assembly.bond_list)} bonds",
        f"{len(assembly.angle_list)} angles",
        f"{len(assembly.dihedral_list)} dihedrals",
        "0 impropers",
        "",
        f"{len(tables['atom'])} atom types",
        f"{len(tables['bond'])} bond types",
        f"{len(tables['angle'])} angle types",
        f"{len(tables['dihedral'])} dihedral types",
        "0 improper types",
        "",
        f"{box['xlo']:.8f} {box['xhi']:.8f} xlo xhi",
        f"{box['ylo']:.8f} {box['yhi']:.8f} ylo yhi",
        f"{box['zlo']:.8f} {box['zhi']:.8f} zlo zhi",
    ]
    if any(abs(box[k]) > 1e-10 for k in ("xy", "xz", "yz")):
        lines.append(f"{box['xy']:.8f} {box['xz']:.8f} {box['yz']:.8f} xy xz yz")
    lines += ["", "Masses", ""]
    for label, tid in tables["atom"].items():
        lines.append(f"{tid} {mass_of[label]:.6f} # {label}")
    lines += ["", "Atoms # full", ""]
    charges = assembly.charges_e if assembly.charges_e is not None else np.zeros(assembly.n_atoms)
    # Coordinates inside the box plus image flags (ix iy iz): molecules stay whole
    # in LAMMPS (bonds never span the box) and unwrapped coordinates are exact.
    wrapped, images = wrap_with_images(assembly.positions_A, assembly.cell_A, assembly.pbc)
    for k, (pos, img, label) in enumerate(zip(wrapped, images, assembly.type_labels)):
        lines.append(f"{k + 1} {int(assembly.molecule_id[k])} {tables['atom'][label]} "
                     f"{charges[k]:.6f} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f} "
                     f"{img[0]} {img[1]} {img[2]}")
    for title, items, width, table in (("Bonds", assembly.bond_list, 2, tables["bond"]),
                                       ("Angles", assembly.angle_list, 3, tables["angle"]),
                                       ("Dihedrals", assembly.dihedral_list, 4, tables["dihedral"])):
        if not items:
            continue
        lines += ["", title, ""]
        for k, item in enumerate(items, start=1):
            ids = " ".join(str(i + 1) for i in item[:width])
            lines.append(f"{k} {table[item[width]]} {ids}")
    return "\n".join(lines) + "\n"


BOUNDARY_NOTE = {"p p f": "slab or sheet: periodic in-plane, vacuum along z",
                 "f f p": "infinite nanotube: periodic along z, vacuum around it"}


def input_template_text(assembly: AssemblyModel, freeze_substrate: bool = False) -> str:
    tables = type_tables(assembly)
    boundary = recommended_boundaries(assembly.pbc)["lammps_boundary"]

    def listing(kind):
        return "\n".join(f"#   {tid}: {name}" for name, tid in tables[kind].items()) or "#   (none)"

    substrate_types = sorted({tables["atom"][lab] for lab, comp in
                              zip(assembly.type_labels, assembly.component) if comp == "substrate"})
    freeze = (f"group           substrate type {' '.join(map(str, substrate_types))}\n"
              "fix             freeze substrate setforce 0.0 0.0 0.0\n") if freeze_substrate else (
              "# (substrate not frozen; set freeze_substrate in the GUI to add a setforce fix)\n")
    return f"""# LAMMPS input TEMPLATE written by nanobuilder {__version__}
# NOT RUNNABLE AS-IS: no force field was generated. Every line marked REQUIRED
# must be filled in with a reviewed, compatible force field before running MD.

units           real            # positions in Angstrom, masses in amu ('metal' also fits)
atom_style      full
boundary        {boundary}           # {BOUNDARY_NOTE.get(boundary, "isolated particle in a vacuum-padded box")}
read_data       structure.data

# Atom types (structural labels, NOT force-field types):
{listing("atom")}
# Bond types:
{listing("bond")}
# Angle types:
{listing("angle")}
# Dihedral types:
{listing("dihedral")}

# REQUIRED: pair_style / pair_coeff for every type pair, including
#           substrate-ligand cross interactions and the substrate itself
#           (e.g. a many-body metal potential combined with the ligand model)
# REQUIRED: bond_style / bond_coeff (the Au-S bond needs a reviewed model)
# REQUIRED: angle_style / angle_coeff, dihedral_style / dihedral_coeff
# REQUIRED: real partial charges (the data file holds q = 0 placeholders)
# REQUIRED: special_bonds settings consistent with the force field

{freeze}
# run 0          # smoke test once the force field is complete
"""


def export_lammps_data(assembly: AssemblyModel, path, freeze_substrate: bool = False) -> ExportReport:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data_file_text(assembly), encoding="ascii", newline="\n")
    (path.parent / "in.lammps.template").write_text(
        input_template_text(assembly, freeze_substrate), encoding="ascii", newline="\n")
    problems = check_round_trip(assembly, path)
    return ExportReport(
        path=str(path), format="lammps_data", sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        n_atoms=assembly.n_atoms, n_bonds=len(assembly.bond_list),
        round_trip_passed=not problems,
        warnings=problems + ["partial_charges_unassigned=true: q = 0 placeholders satisfy "
                             "atom_style full syntax only",
                             "Styles and coefficients must be supplied (see in.lammps.template)"])


def read_lammps_data(path) -> dict:
    """Minimal parser for the sections this exporter writes (used by round-trip tests)."""
    text = Path(path).read_text(encoding="ascii").splitlines()
    header, sections, current = {}, {}, None
    known = {"Masses", "Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"}
    for raw in text[1:]:
        line = raw.split("#")[0].strip()
        if not line:
            continue
        word = line.split()[0]
        if word in known and len(line.split()) == 1:
            current = word
            sections[current] = []
            continue
        if current is None:
            parts = line.split()
            if parts[-1] in ("atoms", "bonds", "angles", "dihedrals", "impropers"):
                header[parts[-1]] = int(parts[0])
            elif parts[-1] == "types":
                header[" ".join(parts[-2:])] = int(parts[0])
            elif parts[-2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
                header[parts[-2]] = (float(parts[0]), float(parts[1]))
            elif parts[-3:] == ["xy", "xz", "yz"]:
                header["tilt"] = tuple(float(p) for p in parts[:3])
        else:
            sections[current].append(line.split())
    return {"header": header, "sections": sections}


def check_round_trip(assembly: AssemblyModel, path) -> list:
    problems = []
    data = read_lammps_data(path)
    header, sections = data["header"], data["sections"]
    if header.get("atoms") != assembly.n_atoms:
        problems.append("Atom count mismatch")
    atoms = sections.get("Atoms", [])
    ids = [int(row[0]) for row in atoms]
    if sorted(ids) != list(range(1, assembly.n_atoms + 1)):
        problems.append("Atom IDs are not 1..N without duplicates")
    n_types = header.get("atom types", 0)
    if any(not 1 <= int(row[2]) <= n_types for row in atoms):
        problems.append("Atom type ID out of range")
    if len(sections.get("Masses", [])) != n_types:
        problems.append("Masses section does not list every atom type")
    pos = np.array([[float(v) for v in row[4:7]] for row in atoms])
    if atoms and len(atoms[0]) >= 10:  # image flags present: compare unwrapped coordinates
        images = np.array([[int(v) for v in row[7:10]] for row in atoms])
        pos = pos + images @ np.asarray(assembly.cell_A)
    if pos.shape == assembly.positions_A.shape and not np.allclose(pos, assembly.positions_A, atol=1e-5):
        problems.append("Coordinates differ")
    for title, key, width in (("Bonds", "bonds", 2), ("Angles", "angles", 3),
                              ("Dihedrals", "dihedrals", 4)):
        rows = sections.get(title, [])
        if len(rows) != header.get(key, 0):
            problems.append(f"{title}: header says {header.get(key, 0)}, section has {len(rows)}")
        for row in rows:
            if any(not 1 <= int(i) <= assembly.n_atoms for i in row[2:2 + width]):
                problems.append(f"{title}: atom index out of range")
                break
    lx = header.get("xlo", (0.0, 0.0))
    xy = header.get("tilt", (0.0, 0.0, 0.0))[0]
    if abs(xy) > 0.5 * (lx[1] - lx[0]) + 1e-6:
        problems.append("Tilt factor xy exceeds half the box length (LAMMPS rejects this)")
    return problems
