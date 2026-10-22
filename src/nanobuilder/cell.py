"""Simulation-cell handling and engine box conventions (spec stage B3).

The canonical cell is kept lower-triangular (a along x, b in the xy plane),
which is exactly LAMMPS's restricted-triclinic form, with in-plane tilt
reduced so |xy| <= lx/2.
"""

from __future__ import annotations

import numpy as np

from .schema import BuilderError
from ase.geometry import cell_to_cellpar, cellpar_to_cell


def is_lower_triangular(cell: np.ndarray, tol: float = 1e-8) -> bool:
    cell = np.asarray(cell)
    return abs(cell[0, 1]) < tol and abs(cell[0, 2]) < tol and abs(cell[1, 2]) < tol


def to_lower_triangular(cell: np.ndarray, positions: np.ndarray):
    """Rotate cell and positions into the lower-triangular convention."""
    cell = np.asarray(cell, dtype=float)
    if is_lower_triangular(cell):
        return cell.copy(), np.asarray(positions, dtype=float).copy()
    new_cell = cellpar_to_cell(cell_to_cellpar(cell))
    frac = np.linalg.solve(cell.T, np.asarray(positions, dtype=float).T).T
    return new_cell, frac @ new_cell


def reduce_tilt(cell: np.ndarray, pbc) -> np.ndarray:
    """Replace b (and c) by equivalent lattice vectors with the smallest tilt."""
    cell = np.array(cell, dtype=float)
    if pbc[0] and pbc[1]:
        cell[1] -= np.round(cell[1, 0] / cell[0, 0]) * cell[0]
    if pbc[2]:
        if pbc[1]:
            cell[2] -= np.round(cell[2, 1] / cell[1, 1]) * cell[1]
        if pbc[0]:
            cell[2] -= np.round(cell[2, 0] / cell[0, 0]) * cell[0]
    return cell


def wrap_positions(positions: np.ndarray, cell: np.ndarray, pbc) -> np.ndarray:
    frac = np.linalg.solve(np.asarray(cell).T, np.asarray(positions).T).T
    for k in range(3):
        if pbc[k]:
            frac[:, k] %= 1.0
    return frac @ np.asarray(cell)


def wrap_with_images(positions: np.ndarray, cell: np.ndarray, pbc):
    """Wrap positions into the cell along periodic axes and return integer image flags.

    ``wrapped + images @ cell`` reproduces the input positions exactly, which is
    the LAMMPS / HOOMD image-flag convention.
    """
    cell = np.asarray(cell, dtype=float)
    frac = np.linalg.solve(cell.T, np.asarray(positions, dtype=float).T).T
    images = np.zeros(frac.shape, dtype=int)
    for k in range(3):
        if pbc[k]:
            images[:, k] = np.floor(frac[:, k] + 1e-12).astype(int)
            frac[:, k] -= images[:, k]
    return frac @ cell, images


def make_molecules_whole(positions: np.ndarray, cell: np.ndarray, pbc, bonds, molecule_id,
                         roots: dict) -> np.ndarray:
    """Undo per-atom wrapping so every ligand is one connected piece.

    ``roots`` maps molecule ID -> the atom it hangs from (for a ligand, the
    surface atom it is attached to). Starting there, each bonded neighbour is
    placed at the minimum-image bond vector from its parent, so no bond of the
    molecule (or its attachment bond) crosses the cell. Atoms of molecules not
    in ``roots`` (the periodic substrate) are left as they are.
    """
    pbc = tuple(bool(p) for p in pbc)
    positions = np.array(positions, dtype=float)
    if not any(pbc) or not roots:
        return positions
    cell = np.asarray(cell, dtype=float)
    molecule_id = np.asarray(molecule_id)
    neighbours = {}
    for i, j, *_ in bonds:
        neighbours.setdefault(i, []).append(j)
        neighbours.setdefault(j, []).append(i)
    inv = np.linalg.inv(cell)

    def mic(vec):
        frac = vec @ inv
        for k in range(3):
            if pbc[k]:
                frac[k] -= np.round(frac[k])
        return frac @ cell

    for mol, root in roots.items():
        seen = {root}
        queue = [root]
        while queue:
            parent = queue.pop()
            for child in neighbours.get(parent, ()):
                if child in seen or molecule_id[child] != mol:
                    continue  # never move the substrate or another molecule
                positions[child] = positions[parent] + mic(positions[child] - positions[parent])
                seen.add(child)
                queue.append(child)
    return positions


MIN_CLEARANCE_A = 2.0  # atoms (ligands included) must stay this far inside a user-set box
SNUG_CLEARANCE_A = 8.0  # below this, periodic images of an isolated model come close


def _fit_error(axis: str, extent: float, length: float) -> BuilderError:
    return BuilderError(
        f"The model is {extent / 10:.2f} nm across along {axis} (ligands included) but the box "
        f"is {length / 10:.2f} nm: enlarge the box to at least "
        f"{(extent + 2 * MIN_CLEARANCE_A) / 10:.2f} nm or make the structure smaller")


def finalize_box(positions: np.ndarray, cell: np.ndarray, pbc, vacuum_A: float,
                 box_shape: str = "auto", min_lengths_A=None, warnings=None):
    """Size the box along nonperiodic directions and centre the model.

    Without user lengths, each nonperiodic direction is the extent of the atoms
    (substrate and ligands) plus ``vacuum_A`` on each side; ``box_shape``
    "cube" then makes those lengths equal (and a slab box at least as tall as
    its in-plane side). With ``min_lengths_A`` (x, y, z in A; 0 = automatic),
    the box has exactly those lengths, the model is centred, and a model that
    does not fit (``MIN_CLEARANCE_A`` from every face) is an error. Periodic
    directions keep their lattice vectors (the slab builder already matched
    them to the requested box). Returns (positions, cell).
    """
    cell, positions = to_lower_triangular(cell, positions)
    pbc = tuple(bool(p) for p in pbc)
    want = np.zeros(3) if min_lengths_A is None else np.array([float(x) for x in min_lengths_A])
    free = [k for k in range(3) if not pbc[k]]

    def fit(lengths, extent):
        lengths = np.array(lengths, dtype=float)
        if box_shape == "cube" and len(free) > 1 and not np.any(want[free] > 0):
            lengths[free] = max(lengths[free])
        for k in free:
            if want[k] > 0:
                if extent[k] + 2 * MIN_CLEARANCE_A > want[k] + 1e-6:
                    raise _fit_error("xyz"[k], extent[k], want[k])
                lengths[k] = want[k]
                if warnings is not None and extent[k] + 2 * SNUG_CLEARANCE_A > want[k]:
                    warnings.append(f"Box {'xyz'[k]} = {want[k] / 10:.2f} nm leaves less than "
                                    f"{SNUG_CLEARANCE_A:.0f} A between the model and the box "
                                    "faces; periodic images may interact if the box is treated "
                                    "as periodic")
        return lengths

    if not any(pbc):
        lo, hi = positions.min(axis=0), positions.max(axis=0)
        extent = hi - lo
        lengths = fit(extent + 2 * vacuum_A, extent)
        return positions - lo + (lengths - extent) / 2, np.diag(lengths)
    if pbc == (True, True, False):
        z = positions[:, 2]
        cell = reduce_tilt(cell, pbc)
        extent = float(z.max() - z.min())
        lz = extent + 2 * vacuum_A
        if want[2] > 0:
            lz = want[2]
            if box_shape == "cube":  # periodic x and y round to whole repeats: z follows them
                lz = 0.5 * (float(np.linalg.norm(cell[0])) + float(np.linalg.norm(cell[1])))
            if extent + 2 * MIN_CLEARANCE_A > lz + 1e-6:
                raise _fit_error("z", extent, lz)
        elif box_shape == "cube":  # as cubic as the slab allows: height at least the in-plane side
            lz = max(lz, float(np.linalg.norm(cell[0])), float(np.linalg.norm(cell[1])))
        cell[2] = [0.0, 0.0, lz]
        positions = positions.copy()
        positions[:, 2] += (lz - extent) / 2 - z.min()
        return wrap_positions(positions, cell, pbc), cell
    if pbc == (False, False, True):
        # Periodic along z only (infinite nanotube): fit x and y, keep the axial period.
        if abs(cell[2, 0]) > 1e-8 or abs(cell[2, 1]) > 1e-8:
            raise ValueError("A cell periodic only along z must have c parallel to z")
        lo, hi = positions[:, :2].min(axis=0), positions[:, :2].max(axis=0)
        extent = hi - lo
        lengths = fit([extent[0] + 2 * vacuum_A, extent[1] + 2 * vacuum_A, float(cell[2, 2])],
                      np.array([extent[0], extent[1], 0.0]))
        new = np.diag(lengths)
        positions = positions.copy()
        positions[:, :2] += (lengths[:2] - extent) / 2 - lo
        return wrap_positions(positions, new, pbc), new
    raise ValueError(f"Unsupported periodicity {pbc}")


def center_boundaries(atoms, axes=(0, 1)) -> None:
    """Move periodic cell boundaries into the widest empty band between atoms (in place).

    A periodic cell can start anywhere; starting it where the atoms are densest
    splits rows of atoms between the two sides of the box, so the cut looks
    ragged. Shifting each periodic boundary to the middle of the largest gap in
    fractional coordinates gives the cleanest straight edge the lattice allows.
    The structure is unchanged (pure translation, then wrap).
    """
    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=True)
    for k in axes:
        f = np.sort(frac[:, k] % 1.0)
        if len(f) < 2:
            continue
        gaps = np.diff(np.append(f, f[0] + 1.0))
        i = int(np.argmax(gaps))
        boundary = (f[i] + gaps[i] / 2.0) % 1.0
        frac[:, k] = (frac[:, k] - boundary) % 1.0
    atoms.set_scaled_positions(frac)
    atoms.set_cell(cell)


def lammps_box(cell: np.ndarray) -> dict:
    """LAMMPS restricted-triclinic box parameters for a lower-triangular cell."""
    cell = np.asarray(cell, dtype=float)
    if not is_lower_triangular(cell):
        raise ValueError("Cell must be lower-triangular; call to_lower_triangular first")
    return {
        "xlo": 0.0, "xhi": cell[0, 0], "ylo": 0.0, "yhi": cell[1, 1], "zlo": 0.0, "zhi": cell[2, 2],
        "xy": cell[1, 0], "xz": cell[2, 0], "yz": cell[2, 1],
    }


def hoomd_box(cell: np.ndarray) -> list:
    """HOOMD box ``[Lx, Ly, Lz, xy, xz, yz]`` (tilts are dimensionless in HOOMD)."""
    cell = np.asarray(cell, dtype=float)
    if not is_lower_triangular(cell):
        raise ValueError("Cell must be lower-triangular; call to_lower_triangular first")
    lx, ly, lz = cell[0, 0], cell[1, 1], cell[2, 2]
    return [lx, ly, lz, cell[1, 0] / ly, cell[2, 0] / lz, cell[2, 1] / lz]


def hoomd_positions(positions: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Shift positions so the box is centred on the origin, as HOOMD requires."""
    cell = np.asarray(cell, dtype=float)
    shifted = wrap_positions(positions, cell, (True, True, True)) - cell.sum(axis=0) / 2
    return shifted


def recommended_boundaries(pbc) -> dict:
    """Engine boundary advice; these belong in the engine input, not the coordinate file."""
    lammps = " ".join("p" if p else "f" for p in pbc)
    return {
        "lammps_boundary": lammps,
        "note": ("Nonperiodic directions are padded with vacuum. HOOMD-blue boxes are always "
                 "periodic, so the vacuum must exceed the interaction cutoff to suppress image "
                 "interactions; a large vacuum box is not equivalent to a true nonperiodic "
                 "boundary for all electrostatics methods."),
    }
