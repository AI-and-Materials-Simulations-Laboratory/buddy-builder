"""Engine-neutral exports: extended XYZ and the authoritative topology.json (spec F4)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

from .schema import AssemblyModel, ExportReport
from .topology import parameter_coverage


def to_ase_atoms(assembly: AssemblyModel) -> Atoms:
    atoms = Atoms(assembly.elements, positions=assembly.positions_A, cell=assembly.cell_A,
                  pbc=assembly.pbc)
    atoms.set_masses(assembly.masses_amu)
    atoms.new_array("molecule_id", np.asarray(assembly.molecule_id, dtype=int))
    atoms.new_array("component", np.array(assembly.component, dtype="U9"))
    atoms.new_array("chemical_role", np.array(assembly.chemical_role, dtype="U9"))
    atoms.new_array("type_label", np.array(assembly.type_labels, dtype="U16"))
    if assembly.residue_names is not None:
        atoms.new_array("residue", np.array(assembly.residue_names, dtype="U4"))
    return atoms


def topology_dict(assembly: AssemblyModel) -> dict:
    return {
        "format": "nanobuilder-topology",
        "format_version": "0.1",
        "index_convention": ("zero-based internal indices; LAMMPS IDs are index + 1; GSD groups "
                             "use these zero-based indices"),
        "units": {"length": "Angstrom", "mass": "amu"},
        "n_atoms": assembly.n_atoms,
        "cell_A": np.asarray(assembly.cell_A).round(8).tolist(),
        "pbc": list(assembly.pbc),
        "atoms": [
            {"index": k, "element": e, "type_label": t, "chemical_role": r, "component": c,
             "molecule_id": int(m), "mass_amu": round(float(mass), 6)}
            for k, (e, t, r, c, m, mass) in enumerate(zip(
                assembly.elements, assembly.type_labels, assembly.chemical_role,
                assembly.component, assembly.molecule_id, assembly.masses_amu))
        ],
        "bonds": [[int(i), int(j), t] for i, j, t in assembly.bond_list],
        "angles": [[*map(int, a[:3]), a[3]] for a in assembly.angle_list],
        "dihedrals": [[*map(int, d[:4]), d[4]] for d in assembly.dihedral_list],
        "impropers": [],
        "site_to_ligand": assembly.site_to_ligand,
        "charges": {"assigned": assembly.charges_e is not None,
                    "values_e": None if assembly.charges_e is None
                    else np.asarray(assembly.charges_e).tolist()},
        "atom_types_force_field": assembly.atom_types,
        "suggested_atom_types": assembly.suggested_types or None,
        "bond_orders": None if assembly.bond_orders is None else
        [float(o) for o in assembly.bond_orders],
        "residues": None if assembly.residue_names is None else
        [[n, int(r)] for n, r in zip(assembly.residue_names, assembly.residue_ids)],
        "parameter_status": assembly.parameter_status,
        "parameter_coverage": parameter_coverage(assembly),
    }


def export_generic(assembly: AssemblyModel, output_directory) -> list:
    """Write structure.extxyz and topology.json; returns one ExportReport per file."""
    out = Path(output_directory)
    out.mkdir(parents=True, exist_ok=True)
    reports = []

    xyz = out / "structure.extxyz"
    write(xyz, to_ase_atoms(assembly), format="extxyz")
    back = read(xyz, format="extxyz")
    problems = []
    if len(back) != assembly.n_atoms or back.get_chemical_symbols() != list(assembly.elements):
        problems.append("Elements or atom count differ after reading back")
    elif not np.allclose(back.get_positions(), assembly.positions_A, atol=1e-5):
        problems.append("Positions differ after reading back")
    elif not np.array_equal(back.arrays["molecule_id"], assembly.molecule_id):
        problems.append("Molecule IDs differ after reading back")
    reports.append(ExportReport(path=str(xyz), format="extxyz",
                                sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
                                n_atoms=assembly.n_atoms, n_bonds=0,
                                round_trip_passed=not problems,
                                warnings=problems + ["extXYZ carries no bonds; topology.json is "
                                                     "the authoritative connectivity"]))

    topo = out / "topology.json"
    topo.write_text(json.dumps(topology_dict(assembly), indent=1), encoding="utf-8")
    reports.append(ExportReport(path=str(topo), format="topology_json",
                                sha256=hashlib.sha256(topo.read_bytes()).hexdigest(),
                                n_atoms=assembly.n_atoms, n_bonds=len(assembly.bond_list),
                                round_trip_passed=True))
    return reports
