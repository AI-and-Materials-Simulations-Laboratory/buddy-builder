"""Checks on the canonical model before export (spec D4/E, test matrix)."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from scipy.spatial import cKDTree

from .schema import AssemblyModel, ValidationReport

OVERLAP_A = 0.7  # any nonbonded pair closer than this is a malformed structure
SUBSTRATE_DUPLICATE_A = 0.5
MAX_BOND_A = {"Au-S": 3.5}
DEFAULT_MAX_BOND_A = 2.3
MAX_SURFACE_BOND_A = 3.5  # anchor bonded to a crystal surface atom (chemisorption)


def _as_atoms(assembly: AssemblyModel) -> Atoms:
    return Atoms(assembly.elements, positions=assembly.positions_A, cell=assembly.cell_A,
                 pbc=assembly.pbc)


def _excluded_pairs(assembly: AssemblyModel) -> set:
    """1-2 (bonded) and 1-3 (angle end) pairs, which are close by design."""
    pairs = {tuple(sorted(b[:2])) for b in assembly.bond_list}
    pairs |= {tuple(sorted((a[0], a[2]))) for a in assembly.angle_list}
    return pairs


def validate_assembly(assembly: AssemblyModel) -> ValidationReport:
    """Check atom counts, coordinates, topology, geometry, and feasibility."""
    errors, warnings = [], []
    n = assembly.n_atoms
    for name in ("positions_A", "masses_amu", "molecule_id", "component", "type_labels",
                 "chemical_role"):
        if len(getattr(assembly, name)) != n:
            errors.append(f"{name} has {len(getattr(assembly, name))} entries for {n} atoms")
    pos = np.asarray(assembly.positions_A)
    if not np.all(np.isfinite(pos)):
        errors.append("Non-finite atom coordinates")
    if np.any(np.asarray(assembly.masses_amu) <= 0):
        errors.append("Missing or non-positive atom masses")
    if abs(np.linalg.det(assembly.cell_A)) < 1e-6:
        errors.append("Simulation cell is singular")

    seen = set()
    for i, j, kind in assembly.bond_list:
        if not (0 <= i < n and 0 <= j < n) or i == j:
            errors.append(f"Malformed bond ({i}, {j}, {kind})")
            continue
        key = tuple(sorted((i, j)))
        if key in seen:
            errors.append(f"Duplicate bond {key}")
        seen.add(key)
    for tup in assembly.angle_list + assembly.dihedral_list + assembly.improper_list:
        if not all(0 <= t < n for t in tup[:-1]):
            errors.append(f"Topology tuple {tup} references a missing atom")
            break

    if errors:
        return _report("error", errors, warnings, assembly, None)

    # Atoms must sit inside the box along nonperiodic directions.
    frac = np.linalg.solve(np.asarray(assembly.cell_A).T, pos.T).T
    for k in range(3):
        if not assembly.pbc[k] and (frac[:, k].min() < 0 or frac[:, k].max() > 1):
            errors.append(f"Atoms extend outside the box along nonperiodic axis {'xyz'[k]}")

    atoms = _as_atoms(assembly)
    for i, j, kind in assembly.bond_list:
        d = atoms.get_distance(i, j, mic=True)
        limit = MAX_BOND_A.get(kind, DEFAULT_MAX_BOND_A)
        if assembly.residue_names is None and "substrate" in (assembly.chemical_role[i],
                                                               assembly.chemical_role[j]):
            limit = max(limit, MAX_SURFACE_BOND_A)  # metal/oxide-anchor bonds are longer
        if d > limit or d < 0.6:
            errors.append(f"Bond {i}-{j} ({kind}) has implausible length {d:.2f} A")
            break

    excluded = _excluded_pairs(assembly)
    attach = {(min(s["surface_atom_index"], s["anchor_atom_index"]),
               max(s["surface_atom_index"], s["anchor_atom_index"]))
              for s in assembly.site_to_ligand}
    min_nonbonded, error = _contacts(assembly, pos, excluded | attach)
    if error:
        errors.append(error)
    if assembly.residue_names is not None:  # covalent (carbon) substrate: chemistry checks
        errors += _covalent_checks(assembly, pos, excluded)

    if assembly.achieved_ligands < assembly.target_ligands:
        warnings.append(f"Only {assembly.achieved_ligands} of {assembly.target_ligands} ligands "
                        "were placed")
    if assembly.charges_e is None:
        warnings.append("Partial charges unassigned (q = 0 placeholders in LAMMPS atom_style full)")
    if assembly.parameter_status != "complete_validated":
        warnings.append("Force-field parameters unassigned: files are MD-readable structure/"
                        "topology (Level B), not simulation-ready")

    status = "error" if errors else (
        "partial_placement" if assembly.achieved_ligands < assembly.target_ligands else "ok")
    return _report(status, errors, warnings, assembly,
                   None if not np.isfinite(min_nonbonded) else float(min_nonbonded))


CONTACT_RANGE_A = 3.5


def _with_images(positions: np.ndarray, cell: np.ndarray, pbc) -> tuple:
    """Positions plus their periodic images along every periodic axis (-1, 0, +1)."""
    shifts = [np.zeros(3)]
    for k in range(3):
        if pbc[k]:
            shifts = [s + m * cell[k] for s in shifts for m in (-1, 0, 1)]
    images = np.concatenate([positions + s for s in shifts])
    index = np.tile(np.arange(len(positions)), len(shifts))
    return images, index


def _contacts(assembly: AssemblyModel, pos: np.ndarray, skip: set) -> tuple:
    """(closest nonbonded contact involving a ligand atom, first error or None).

    Only two kinds of pairs matter: substrate atoms closer than
    SUBSTRATE_DUPLICATE_A (duplicated atoms), and any pair that involves a
    ligand atom within CONTACT_RANGE_A. Searching just those with KD-trees
    replaces a full all-pairs neighbour list (about 1 GB and 10 s for 70k atoms).
    Periodic images are included along periodic axes.
    """
    cell = np.asarray(assembly.cell_A, dtype=float)
    component = np.asarray(assembly.component)
    is_sub = component == "substrate"
    sub_idx = np.flatnonzero(is_sub)
    lig_idx = np.flatnonzero(~is_sub)

    if len(sub_idx):
        images, index = _with_images(pos[sub_idx], cell, assembly.pbc)
        tree = cKDTree(images)
        for a, b in sorted(tree.query_pairs(SUBSTRATE_DUPLICATE_A)):
            i, j = sub_idx[index[a]], sub_idx[index[b]]
            if i != j:
                d = float(np.linalg.norm(images[a] - images[b]))
                return np.inf, (f"Substrate atoms {min(i, j)} and {max(i, j)} overlap "
                                f"({d:.2f} A); duplicated atoms?")

    if not len(lig_idx):
        return np.inf, None
    images, index = _with_images(pos, cell, assembly.pbc)
    tree = cKDTree(images)
    min_nonbonded = np.inf
    for i, hits in zip(lig_idx, tree.query_ball_point(pos[lig_idx], CONTACT_RANGE_A)):
        for n in hits:
            j = int(index[n])
            if j == i or (min(i, j), max(i, j)) in skip:
                continue
            d = float(np.linalg.norm(pos[i] - images[n]))
            if d < min_nonbonded:
                min_nonbonded = d
            if d < OVERLAP_A:
                return min_nonbonded, f"Atoms {min(i, j)} and {max(i, j)} overlap ({d:.2f} A)"
    return min_nonbonded, None


MAX_VALENCE = {"H": 1, "O": 2, "C": 4, "N": 4, "S": 6, "P": 5, "F": 1, "Cl": 1, "Br": 1, "I": 1}
COVALENT_OVERLAP_A = 0.9


def _covalent_checks(assembly: AssemblyModel, pos: np.ndarray, excluded: set) -> list:
    """Valence and overlap checks over every atom of a covalent (carbon) model."""
    errors = []
    degree = np.zeros(assembly.n_atoms, dtype=int)
    for i, j, *_ in assembly.bond_list:
        degree[i] += 1
        degree[j] += 1
    for k, (e, d) in enumerate(zip(assembly.elements, degree)):
        if d > MAX_VALENCE.get(e, 12):
            errors.append(f"Atom {k} ({e}, {assembly.type_labels[k]}) has {d} bonds")
            break
    for k in range(assembly.n_substrate):
        if assembly.elements[k] == "C" and degree[k] < 2:
            errors.append(f"Substrate carbon {k} has only {degree[k]} bonds (dangling atom)")
            break
    images, index = _with_images(pos, np.asarray(assembly.cell_A, dtype=float), assembly.pbc)
    tree = cKDTree(images)
    for i, hits in enumerate(tree.query_ball_point(pos, COVALENT_OVERLAP_A)):
        for n in hits:
            j = int(index[n])
            if j != i and (min(i, j), max(i, j)) not in excluded:
                d = float(np.linalg.norm(pos[i] - images[n]))
                return errors + [f"Atoms {min(i, j)} and {max(i, j)} overlap ({d:.2f} A)"]
    return errors


def _report(status, errors, warnings, assembly: AssemblyModel, min_nb) -> ValidationReport:
    n_lig_atoms = sum(1 for c in assembly.component if c == "ligand")
    return ValidationReport(
        status=status, errors=errors, warnings=warnings,
        target_ligands=assembly.target_ligands, placed_ligands=assembly.achieved_ligands,
        requested_density_nm2=assembly.requested_density_nm2,
        achieved_density_nm2=assembly.achieved_density_nm2,
        min_nonbonded_distance_A=min_nb, parameter_status=assembly.parameter_status,
        counts={
            "atoms": assembly.n_atoms,
            "substrate_atoms": assembly.n_substrate,
            "ligand_atoms": n_lig_atoms,
            "ligands": assembly.achieved_ligands,
            "bonds": len(assembly.bond_list),
            "angles": len(assembly.angle_list),
            "dihedrals": len(assembly.dihedral_list),
            "surface_area_nm2": round(assembly.surface_area_nm2, 4),
        })
