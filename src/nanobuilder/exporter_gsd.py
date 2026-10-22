"""HOOMD-blue GSD export through the GSD Python API (spec F2)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import gsd.hoomd
import numpy as np

from .cell import hoomd_box, wrap_with_images
from .schema import AssemblyModel, ExportReport


def _typed(items, name_index: int):
    names = sorted({t[name_index] for t in items})
    lookup = {n: k for k, n in enumerate(names)}
    group = np.array([t[:name_index] for t in items], dtype=np.uint32).reshape(len(items), name_index)
    typeid = np.array([lookup[t[name_index]] for t in items], dtype=np.uint32)
    return names, typeid, group


def build_frame(assembly: AssemblyModel) -> gsd.hoomd.Frame:
    labels = assembly.type_labels
    types = sorted(set(labels))
    lookup = {t: k for k, t in enumerate(types)}
    frame = gsd.hoomd.Frame()
    frame.configuration.step = 0
    frame.configuration.dimensions = 3
    frame.configuration.box = hoomd_box(assembly.cell_A)
    frame.particles.N = assembly.n_atoms
    frame.particles.types = types
    frame.particles.typeid = np.array([lookup[t] for t in labels], dtype=np.uint32)
    # Positions in the origin-centred HOOMD box, plus image flags so molecules
    # that cross a periodic boundary stay whole when unwrapped.
    cell = np.asarray(assembly.cell_A, dtype=float)
    wrapped, images = wrap_with_images(assembly.positions_A, cell, (True, True, True))
    frame.particles.position = (wrapped - cell.sum(axis=0) / 2).astype(np.float32)
    frame.particles.image = images.astype(np.int32)
    frame.particles.mass = np.asarray(assembly.masses_amu, dtype=np.float32)
    if assembly.charges_e is not None:  # only when intentionally assigned
        frame.particles.charge = np.asarray(assembly.charges_e, dtype=np.float32)
    for section, items, width in (("bonds", assembly.bond_list, 2),
                                  ("angles", assembly.angle_list, 3),
                                  ("dihedrals", assembly.dihedral_list, 4)):
        if not items:
            continue
        names, typeid, group = _typed(items, width)
        block = getattr(frame, section)
        block.N = len(items)
        block.types = names
        block.typeid = typeid
        block.group = group
    frame.log["particles/molecule_id"] = np.asarray(assembly.molecule_id, dtype=np.int32)
    return frame


def export_hoomd_gsd(assembly: AssemblyModel, path) -> ExportReport:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = build_frame(assembly)
    with gsd.hoomd.open(name=str(path), mode="w") as traj:
        traj.append(frame)
    problems = check_round_trip(assembly, path)
    return ExportReport(
        path=str(path), format="hoomd_gsd", sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        n_atoms=assembly.n_atoms, n_bonds=len(assembly.bond_list),
        round_trip_passed=not problems,
        warnings=problems + ["GSD stores no force field; HOOMD units are application-defined "
                             "(here: Angstrom, amu)"])


def check_round_trip(assembly: AssemblyModel, path) -> list:
    """Reopen the file and compare structural invariants; returns a list of problems."""
    problems = []
    expected = build_frame(assembly)
    with gsd.hoomd.open(name=str(path), mode="r") as traj:
        if len(traj) != 1:
            return [f"Expected 1 frame, found {len(traj)}"]
        got = traj[0]
    if got.particles.N != assembly.n_atoms:
        problems.append(f"Particle count {got.particles.N} != {assembly.n_atoms}")
    if list(got.particles.types) != list(expected.particles.types):
        problems.append("Particle type names differ")
    if not np.array_equal(got.particles.typeid, expected.particles.typeid):
        problems.append("Particle type IDs differ")
    if not np.allclose(got.particles.position, expected.particles.position, atol=1e-4):
        problems.append("Positions differ")
    if not np.array_equal(got.particles.image, expected.particles.image):
        problems.append("Image flags differ")
    if not np.allclose(got.configuration.box, expected.configuration.box, atol=1e-5):
        problems.append("Box differs")
    # Inside-box check in fractional coordinates of the HOOMD box.
    lx, ly, lz, xy, xz, yz = got.configuration.box
    matrix = np.array([[lx, 0, 0], [xy * ly, ly, 0], [xz * lz, yz * lz, lz]])
    frac = np.linalg.solve(matrix.T, got.particles.position.T.astype(float)).T
    if np.any(frac < -0.5 - 1e-5) or np.any(frac > 0.5 + 1e-5):
        problems.append("Some particles lie outside the HOOMD box")
    for section, items in (("bonds", assembly.bond_list), ("angles", assembly.angle_list),
                           ("dihedrals", assembly.dihedral_list)):
        block = getattr(got, section)
        if block.N != len(items):
            problems.append(f"{section}: {block.N} stored, {len(items)} expected")
        elif items and not np.array_equal(block.group, getattr(expected, section).group):
            problems.append(f"{section}: particle index groups differ")
    return problems
