"""Test matrix 17: GSD reopens with consistent box, types, positions and bonds."""

import gsd.hoomd
import numpy as np
import pytest

from nanobuilder.cell import hoomd_box
from nanobuilder.exporter_gsd import export_hoomd_gsd
from nanobuilder.schema import AttachmentConfig, CoverageConfig
from nanobuilder.placement import place_ligands
from nanobuilder.sites import find_surface_sites, selected_area_nm2


@pytest.fixture(scope="module")
def slab_assembly(au111_slab, scco_thiolate):
    cov = CoverageConfig(mode="density", density_nm2=0.5, region="top")
    sites = find_surface_sites(au111_slab, AttachmentConfig(region="top", surface_species=["Au"]))
    return place_ligands(au111_slab, scco_thiolate, sites, cov, selected_area_nm2(au111_slab, "top"))


def test_round_trip(tmp_path, slab_assembly):
    report = export_hoomd_gsd(slab_assembly, tmp_path / "structure.gsd")
    assert report.round_trip_passed, report.warnings
    with gsd.hoomd.open(str(tmp_path / "structure.gsd"), "r") as traj:
        frame = traj[0]
    a = slab_assembly
    assert frame.particles.N == a.n_atoms
    assert frame.bonds.N == len(a.bond_list)
    assert frame.angles.N == len(a.angle_list)
    assert frame.dihedrals.N == len(a.dihedral_list)
    types = [frame.particles.types[t] for t in frame.particles.typeid]
    assert types == a.type_labels
    # Bond groups are zero-based particle indices.
    names = [frame.bonds.types[t] for t in frame.bonds.typeid]
    assert [(int(i), int(j), n) for (i, j), n in zip(frame.bonds.group, names)] == \
        [(i, j, t) for i, j, t in a.bond_list]
    assert np.allclose(frame.particles.mass, a.masses_amu, rtol=1e-6)
    assert np.allclose(frame.log["particles/molecule_id"], a.molecule_id)


def test_triclinic_box_convention(slab_assembly):
    lx, ly, lz, xy, xz, yz = hoomd_box(slab_assembly.cell_A)
    cell = slab_assembly.cell_A
    assert lx == pytest.approx(cell[0, 0]) and ly == pytest.approx(cell[1, 1])
    assert xy == pytest.approx(cell[1, 0] / cell[1, 1])  # HOOMD tilt is dimensionless
    assert xz == yz == 0


def test_bond_lengths_survive_wrapping(tmp_path, slab_assembly):
    """Positions are wrapped into the HOOMD box; bonds must still be short under minimum image."""
    export_hoomd_gsd(slab_assembly, tmp_path / "s.gsd")
    with gsd.hoomd.open(str(tmp_path / "s.gsd"), "r") as traj:
        frame = traj[0]
    lx, ly, lz, xy, xz, yz = frame.configuration.box
    box = np.array([[lx, 0, 0], [xy * ly, ly, 0], [xz * lz, yz * lz, lz]])
    pos = frame.particles.position.astype(float)
    for i, j in frame.bonds.group:
        d = pos[i] - pos[j]
        frac = np.linalg.solve(box.T, d)
        frac -= np.round(frac)
        assert np.linalg.norm(frac @ box) < 3.0
