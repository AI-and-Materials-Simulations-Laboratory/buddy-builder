"""Test matrix 6-8: sphere, slab planes/normals, noncubic and imported geometry."""

import numpy as np
import pytest

from conftest import upload_config
from nanobuilder.crystal import build_crystal
from nanobuilder.schema import BuilderError, GeometryConfig, MaterialConfig
from nanobuilder.shapes import build_geometry, list_terminations


def test_sphere_3nm(au_sphere):
    atoms = au_sphere.atoms
    r = np.linalg.norm(atoms.positions - au_sphere.center_A, axis=1)
    assert r.max() <= 15.0 + 1e-6
    assert 700 < len(atoms) < 1000  # about 0.059 atoms/A^3 x 14137 A^3
    assert au_sphere.faces["entire_surface"].area_nm2 == pytest.approx(28.274, abs=1e-3)
    assert au_sphere.pbc == (False, False, False)
    d = atoms.get_all_distances()
    np.fill_diagonal(d, np.inf)
    assert d.min() == pytest.approx(4.08 / np.sqrt(2), abs=1e-6)  # no duplicates


def test_sphere_vacuum(au_sphere):
    pos, cell = au_sphere.atoms.positions, np.diag(au_sphere.atoms.cell)
    assert pos.min() >= 30.0 - 1e-6
    assert (cell - pos.max(axis=0)).min() >= 30.0 - 1e-6


@pytest.mark.parametrize("miller, spacing, normal", [
    ((1, 0, 0), 4.08 / 2, [1, 0, 0]),
    ((1, 1, 1), 4.08 / np.sqrt(3), np.ones(3) / np.sqrt(3)),
    ((1, 1, 0), 4.08 / (2 * np.sqrt(2)), [1 / np.sqrt(2), 1 / np.sqrt(2), 0]),
])
def test_au_slab_planes_use_conventional_indices(au, miller, spacing, normal):
    slab = build_geometry(au, GeometryConfig(shape="slab", miller=miller, slab_length_nm=2,
                                             slab_width_nm=2))
    assert slab.interplanar_spacing_A == pytest.approx(spacing, abs=1e-6)
    assert np.allclose(slab.crystal_frame_normal, normal)
    z = np.unique(slab.atoms.positions[:, 2].round(4))
    assert np.allclose(np.diff(z), spacing, atol=1e-4)
    assert len(z) == 4  # four layers


def test_slab_in_plane_area_and_periodicity(au111_slab):
    cell = np.array(au111_slab.atoms.cell)
    area = np.linalg.norm(np.cross(cell[0], cell[1])) / 100
    assert au111_slab.faces["top"].area_nm2 == pytest.approx(area)
    assert area >= 25 * 0.99  # at least the requested 5 x 5 nm (1% tolerance)
    assert au111_slab.pbc == (True, True, False)
    assert abs(cell[1, 0]) <= cell[0, 0] / 2 + 1e-9  # reduced tilt


def test_slab_top_and_bottom_faces(au111_slab):
    z = au111_slab.atoms.positions[:, 2]
    top = au111_slab.faces["top"].atom_indices
    bottom = au111_slab.faces["bottom"].atom_indices
    assert np.allclose(z[top], z.max())
    assert np.allclose(z[bottom], z.min())
    assert len(top) == len(bottom) == 400


def test_thickness_mode(au):
    slab = build_geometry(au, GeometryConfig(shape="slab", miller=(1, 0, 0), thickness_mode="thickness",
                                             thickness_nm=1.0, slab_length_nm=2, slab_width_nm=2))
    assert slab.layers == 5  # 10 A / 2.04 A per (100) layer


def test_hcp_basal_normal_and_upload_geometry(fixtures):
    preset = build_crystal(MaterialConfig(preset_id="Mg_hcp"))
    uploaded = build_crystal(upload_config(fixtures, "mg_hcp_fully_occupied.cif"))
    for crystal in (preset, uploaded):
        slab = build_geometry(crystal, GeometryConfig(shape="slab", miller=(0, 0, 1),
                                                      slab_length_nm=2, slab_width_nm=2))
        assert slab.interplanar_spacing_A == pytest.approx(5.21 / 2, abs=1e-6)
        # (001) normal comes from the reciprocal lattice and is along c (here z).
        assert np.allclose(slab.crystal_frame_normal, [0, 0, 1], atol=1e-9)
        assert slab.terminations[0].top_elements == ["Mg"]


def test_noncubic_reciprocal_normal_differs_from_direction(fixtures):
    crystal = build_crystal(upload_config(fixtures, "mg_hcp_fully_occupied.cif"))
    slab = build_geometry(crystal, GeometryConfig(shape="slab", miller=(1, 0, 0), slab_length_nm=2,
                                                  slab_width_nm=2))
    a = np.array(crystal.atoms.cell[0]) / np.linalg.norm(crystal.atoms.cell[0])
    # (100) plane normal is perpendicular to b and c, not parallel to the [100] direction a.
    assert abs(np.dot(slab.crystal_frame_normal, crystal.atoms.cell[1])) < 1e-9
    assert not np.allclose(slab.crystal_frame_normal, a)


def test_terminations():
    nacl = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    terms = list_terminations(nacl, (1, 1, 1), 4)
    assert [(t.top_elements, t.bottom_elements) for t in terms] == [(["Na"], ["Cl"]), (["Cl"], ["Na"])]
    slab = build_geometry(nacl, GeometryConfig(shape="slab", miller=(1, 1, 1), termination_index=1,
                                               slab_length_nm=2, slab_width_nm=2))
    sym = slab.atoms.get_chemical_symbols()
    z = slab.atoms.positions[:, 2]
    assert sym[int(np.argmax(z))] == "Cl" and sym[int(np.argmin(z))] == "Na"
    assert any("polar" in w for w in slab.warnings)
    assert len(list_terminations(build_crystal(MaterialConfig()), (1, 1, 1))) == 1
    assert len(list_terminations(build_crystal(MaterialConfig(preset_id="Si_diamond")), (1, 1, 1))) == 2


def test_compound_sphere_warns_about_stoichiometry():
    nacl = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    sphere = build_geometry(nacl, GeometryConfig(diameter_nm=1.5))
    comp = sphere.composition
    if comp["Na"] != comp["Cl"]:
        assert any("non-stoichiometric" in w for w in sphere.warnings)
    si = build_geometry(build_crystal(MaterialConfig(preset_id="Si_diamond")), GeometryConfig(diameter_nm=1.5))
    assert any("dangling bonds" in w for w in si.warnings)


def test_invalid_geometry_inputs(au):
    with pytest.raises(BuilderError, match="0, 0, 0"):
        build_geometry(au, GeometryConfig(shape="slab", miller=(0, 0, 0)))
    with pytest.raises(BuilderError, match="limit"):
        build_geometry(au, GeometryConfig(diameter_nm=20, max_atoms=1000))
    with pytest.raises(BuilderError, match="positive"):
        build_geometry(au, GeometryConfig(diameter_nm=-1))


# --- nanorods ----------------------------------------------------------------

import numpy as _np
import pytest as _pytest

from nanobuilder.crystal import build_crystal as _bc
from nanobuilder.schema import GeometryConfig as _GC, MaterialConfig as _MC
from nanobuilder.shapes import build_geometry as _bg


@_pytest.mark.parametrize("caps", ["hemispherical", "flat"])
def test_rod_dimensions_and_regions(caps):
    au = _bc(_MC(preset_id="Au_fcc"))
    s = _bg(au, _GC(shape="rod", diameter_nm=3.0, rod_length_nm=9.0, rod_caps=caps))
    pos = s.atoms.get_positions()
    ext = pos.max(axis=0) - pos.min(axis=0)
    assert abs(ext[2] - 90) < 4 and abs(ext[0] - 30) < 4 and abs(ext[1] - 30) < 4
    assert set(s.faces) == {"side", "tips"}
    assert len(s.faces["side"].atom_indices) > 0 and len(s.faces["tips"].atom_indices) > 0
    assert s.pbc == (False, False, False)
    density = len(au.atoms) / au.atoms.get_volume()
    r, half = 15.0, (45.0 - 15.0 if caps == "hemispherical" else 45.0)
    vol = _np.pi * r ** 2 * 2 * half + (4 / 3 * _np.pi * r ** 3 if caps == "hemispherical" else 0)
    assert abs(len(s.atoms) / (density * vol) - 1) < 0.08


def test_rod_axis_follows_requested_lattice_direction():
    au = _bc(_MC(preset_id="Au_fcc"))
    for axis in [(0, 0, 1), (1, 1, 0), (1, 1, 1)]:
        s = _bg(au, _GC(shape="rod", diameter_nm=2.0, rod_length_nm=6.0, rod_axis=axis))
        pos = s.atoms.get_positions()
        ext = pos.max(axis=0) - pos.min(axis=0)
        assert ext[2] > 1.8 * max(ext[0], ext[1])  # long axis is z for every direction


def test_rod_tips_only_coverage_stays_on_the_tips():
    from nanobuilder.pipeline import run_pipeline
    from nanobuilder.schema import RunConfig
    cfg = RunConfig.from_dict({
        "material": {"source": "preset", "preset_id": "Au_fcc"},
        "geometry": {"shape": "rod", "diameter_nm": 2.5, "rod_length_nm": 8.0, "vacuum_A": 20},
        "ligand": {"smiles": "SCCO", "anchor_atom_index": 0, "binding_template": "Au_thiolate"},
        "coverage": {"mode": "count", "target": 10, "region": "tips"}})
    r = run_pipeline(cfg)
    assert r.assembly.achieved_ligands == 10
    assert {x["face"] for x in r.assembly.site_to_ligand} == {"tips"}
    assert r.report.status == "ok"


def test_rod_config_round_trips_through_dict():
    from nanobuilder.schema import RunConfig
    cfg = RunConfig.from_dict({"geometry": {"shape": "rod", "diameter_nm": 2.0,
                                            "rod_length_nm": 7.0, "rod_axis": [1, 1, 0],
                                            "rod_caps": "flat"}})
    again = RunConfig.from_dict(cfg.to_dict())
    assert again.geometry.rod_axis == (1, 1, 0) and again.geometry.rod_caps == "flat"
    assert "miller" not in cfg.to_dict()["geometry"]


@_pytest.mark.parametrize("preset_id,axis,d,length", [("Au_fcc", (1, 1, 1), 8, 8),
                                                      ("Au_fcc", (1, 1, 0), 6, 8),
                                                      ("Mg_hcp", (1, 0, 0), 3, 8)])
def test_flat_rods_off_axis_are_complete_and_have_tips(preset_id, axis, d, length):
    crystal = _bc(_MC(preset_id=preset_id))
    s = _bg(crystal, _GC(shape="rod", diameter_nm=d, rod_length_nm=length, rod_axis=axis,
                         rod_caps="flat", max_atoms=10 ** 6))
    density = len(crystal.atoms) / crystal.atoms.get_volume()
    r = d * 5.0
    assert len(s.atoms) / (density * _np.pi * r * r * length * 10) > 0.95
    p = s.atoms.positions - s.center_A
    rim = _np.hypot(p[:, 0], p[:, 1]) > r - 3
    slices = _np.histogram(_np.degrees(_np.arctan2(p[rim, 1], p[rim, 0])) % 360, bins=12,
                           range=(0, 360))[0]
    # Facets of different orientation pack differently (the counts follow the crystal's
    # symmetry); before the supercell fix, whole slices were empty.
    assert slices.min() > 0.25 * slices.max()
    assert len(s.faces["tips"].atom_indices) > 0
