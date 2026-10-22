"""Test matrix 18: LAMMPS data file parses with valid one-based IDs, masses and counts.

If a working LAMMPS executable is available, set LAMMPS_EXE to its path to
also run a real ``read_data`` check (skipped otherwise).
"""

import os
import subprocess

import numpy as np
import pytest

from nanobuilder.crystal import build_crystal
from nanobuilder.exporter_lammps import export_lammps_data, read_lammps_data
from nanobuilder.placement import bare_assembly, place_ligands
from nanobuilder.schema import AttachmentConfig, CoverageConfig, GeometryConfig, MaterialConfig
from nanobuilder.shapes import build_geometry
from nanobuilder.sites import find_surface_sites, selected_area_nm2


@pytest.fixture(scope="module")
def sphere_assembly(au_sphere, scco_thiolate):
    sites = find_surface_sites(au_sphere, AttachmentConfig(surface_species=["Au"]))
    return place_ligands(au_sphere, scco_thiolate, sites, CoverageConfig(target=12))


@pytest.fixture(scope="module")
def slab_assembly(au111_slab, scco_thiolate):
    sites = find_surface_sites(au111_slab, AttachmentConfig(region="top", surface_species=["Au"]))
    return place_ligands(au111_slab, scco_thiolate, sites,
                         CoverageConfig(mode="density", density_nm2=0.5, region="top"),
                         selected_area_nm2(au111_slab, "top"))


def test_sphere_data_file(tmp_path, sphere_assembly):
    report = export_lammps_data(sphere_assembly, tmp_path / "structure.data")
    assert report.round_trip_passed, report.warnings
    data = read_lammps_data(tmp_path / "structure.data")
    h, s = data["header"], data["sections"]
    a = sphere_assembly
    assert h["atoms"] == a.n_atoms and h["bonds"] == len(a.bond_list)
    assert h["angles"] == len(a.angle_list) and h["dihedrals"] == len(a.dihedral_list)
    assert "tilt" not in h  # orthogonal box for a sphere
    ids = [int(r[0]) for r in s["Atoms"]]
    assert ids == list(range(1, a.n_atoms + 1))
    # Bond atoms are one-based versions of the internal zero-based indices.
    first = s["Bonds"][0]
    assert (int(first[2]) - 1, int(first[3]) - 1) == a.bond_list[0][:2]
    masses = {int(r[0]): float(r[1]) for r in s["Masses"]}
    assert masses and all(m > 0 for m in masses.values())
    assert all(float(r[3]) == 0.0 for r in s["Atoms"])  # q placeholders
    assert (tmp_path / "in.lammps.template").read_text().count("REQUIRED") >= 4
    assert any("partial_charges_unassigned=true" in w for w in report.warnings)


def test_slab_triclinic_box(tmp_path, slab_assembly):
    export_lammps_data(slab_assembly, tmp_path / "structure.data")
    h = read_lammps_data(tmp_path / "structure.data")["header"]
    lx = h["xlo"][1] - h["xlo"][0]
    xy, xz, yz = h["tilt"]
    assert abs(xy) <= lx / 2 + 1e-9 and xz == 0 and yz == 0
    assert "boundary        p p f" in (tmp_path / "in.lammps.template").read_text()


def test_any_material_masses(tmp_path):
    nacl = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    a = bare_assembly(build_geometry(nacl, GeometryConfig(diameter_nm=1.5)))
    export_lammps_data(a, tmp_path / "s.data")
    s = read_lammps_data(tmp_path / "s.data")["sections"]
    masses = sorted(round(float(r[1]), 2) for r in s["Masses"])
    assert masses == [22.99, 35.45]


def test_freeze_substrate_template(tmp_path, sphere_assembly):
    export_lammps_data(sphere_assembly, tmp_path / "s.data", freeze_substrate=True)
    assert "fix             freeze substrate setforce" in (tmp_path / "in.lammps.template").read_text()


@pytest.mark.skipif(not os.environ.get("LAMMPS_EXE"), reason="set LAMMPS_EXE to run a real read_data")
def test_real_lammps_read_data(tmp_path, slab_assembly):
    export_lammps_data(slab_assembly, tmp_path / "structure.data")
    (tmp_path / "in.check").write_text("units real\natom_style full\nboundary p p f\n"
                                       "read_data structure.data\n")
    run = subprocess.run([os.environ["LAMMPS_EXE"], "-in", "in.check"], cwd=tmp_path,
                         capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stdout + run.stderr
    assert f"{slab_assembly.n_atoms} atoms" in run.stdout
    assert np.isfinite(slab_assembly.positions_A).all()
