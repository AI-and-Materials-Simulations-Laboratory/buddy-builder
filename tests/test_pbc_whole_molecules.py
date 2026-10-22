"""Ligands on periodic slabs must stay whole (no bonds spanning the cell)."""

import copy

import gsd.hoomd
import numpy as np
import pytest

from nanobuilder.cell import make_molecules_whole, wrap_with_images
from nanobuilder.examples import EXAMPLES
from nanobuilder.exporter_gsd import export_hoomd_gsd
from nanobuilder.exporter_lammps import export_lammps_data, read_lammps_data
from nanobuilder.pipeline import run_pipeline
from nanobuilder.schema import RunConfig


@pytest.fixture(scope="module")
def dense_long_chain_slab():
    """Au(111) slab with a C11 thiol at 2/nm^2: many ligands sit on the cell edge."""
    cfg = copy.deepcopy(EXAMPLES["C_au111_slab_density"]["config"])
    cfg["ligand"]["smiles"] = "SCCCCCCCCCCCO"
    cfg["ligand"]["anchor_atom_index"] = 0
    cfg["coverage"]["density_nm2"] = 2.0
    return run_pipeline(RunConfig.from_dict(cfg)).assembly


def _bond_lengths(assembly):
    pos = np.asarray(assembly.positions_A)
    return np.array([np.linalg.norm(pos[i] - pos[j]) for i, j, *_ in assembly.bond_list])


def test_no_bond_spans_the_cell(dense_long_chain_slab):
    a = dense_long_chain_slab
    assert a.achieved_ligands > 20
    # Plain Euclidean lengths (no minimum image): Au-S <= 3.5 A, everything else shorter.
    assert _bond_lengths(a).max() < 3.5


def test_some_ligands_actually_cross_the_boundary(dense_long_chain_slab):
    """Guard that the fixture exercises the fix: some ligand atoms lie outside the cell."""
    a = dense_long_chain_slab
    _, images = wrap_with_images(a.positions_A, a.cell_A, a.pbc)
    ligand = np.asarray(a.molecule_id) > 0
    assert np.any(images[ligand] != 0)


def test_lammps_image_flags_restore_whole_molecules(dense_long_chain_slab, tmp_path):
    a = dense_long_chain_slab
    report = export_lammps_data(a, tmp_path / "s.data")
    assert report.round_trip_passed, report.warnings
    rows = read_lammps_data(tmp_path / "s.data")["sections"]["Atoms"]
    assert all(len(r) == 10 for r in rows)  # id mol type q x y z ix iy iz
    wrapped = np.array([[float(v) for v in r[4:7]] for r in rows])
    frac = np.linalg.solve(np.asarray(a.cell_A).T, wrapped.T).T
    assert frac[:, :2].min() >= -1e-6 and frac[:, :2].max() < 1 + 1e-6  # inside the box
    images = np.array([[int(v) for v in r[7:10]] for r in rows])
    assert np.allclose(wrapped + images @ np.asarray(a.cell_A), a.positions_A, atol=1e-5)


def test_gsd_image_flags_restore_whole_molecules(dense_long_chain_slab, tmp_path):
    a = dense_long_chain_slab
    report = export_hoomd_gsd(a, tmp_path / "s.gsd")
    assert report.round_trip_passed, report.warnings
    with gsd.hoomd.open(str(tmp_path / "s.gsd"), "r") as traj:
        f = traj[0]
    lx, ly, lz, xy, xz, yz = f.configuration.box
    box = np.array([[lx, 0, 0], [xy * ly, ly, 0], [xz * lz, yz * lz, lz]])
    unwrapped = f.particles.position.astype(float) + f.particles.image @ box
    pos = unwrapped[np.asarray(a.molecule_id) > 0]
    # every ligand bond is short after unwrapping
    lig = [(i, j) for i, j, *_ in a.bond_list if a.molecule_id[i] > 0 and a.molecule_id[j] > 0]
    d = [np.linalg.norm(unwrapped[i] - unwrapped[j]) for i, j in lig]
    assert pos.size and max(d) < 2.5


def test_make_whole_leaves_substrate_and_nonperiodic_alone():
    cell = np.diag([10.0, 10.0, 30.0])
    pos = np.array([[0.2, 5, 10], [9.9, 5, 12.4], [0.8, 5, 13.5]])  # substrate, anchor, tail
    out = make_molecules_whole(pos, cell, (True, True, False), [(1, 2)] + [(0, 1)],
                               [0, 1, 1], roots={1: 0})
    assert np.allclose(out[0], pos[0])                 # substrate atom untouched
    assert np.allclose(out[1], [-0.1, 5, 12.4])        # anchor pulled next to its surface atom
    assert np.allclose(out[2], [0.8, 5, 13.5])         # tail already whole relative to the anchor
    assert np.allclose(make_molecules_whole(pos, cell, (False,) * 3, [], [0, 1, 1], {1: 0}), pos)
