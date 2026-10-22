"""Buddy Builder backend: Fibonacci placement, rectangular/cubic boxes, channels, hand edits."""

import numpy as np
import pytest

from nanobuilder.editing import add_ligand_at, remove_ligands
from nanobuilder.fibonacci import GOLDEN_ANGLE, fibonacci_sphere
from nanobuilder.pipeline import export_package, run_pipeline
from nanobuilder.schema import (BuilderError, CarbonConfig, CoverageConfig, GeometryConfig,
                                LigandConfig, MaterialConfig, RunConfig)


def _nn_cv(points):
    d = np.linalg.norm(points[:, None] - points[None], axis=2)
    np.fill_diagonal(d, np.inf)
    nn = d.min(axis=1)
    return nn.std() / nn.mean()


def test_fibonacci_sphere_is_unit_and_golden():
    u = fibonacci_sphere(100)
    assert np.allclose(np.linalg.norm(u, axis=1), 1.0)
    assert abs(np.degrees(GOLDEN_ANGLE) - 137.508) < 1e-3
    assert abs(u[:, 2].mean()) < 1e-12  # balanced between hemispheres


def test_fibonacci_placement_is_deterministic_and_more_even_than_random():
    def build(dist, seed):
        cfg = RunConfig(coverage=CoverageConfig(target=40, distribution=dist, seed=seed))
        a = run_pipeline(cfg).assembly
        return a.positions_A[[s["surface_atom_index"] for s in a.site_to_ligand]], a

    p1, a1 = build("fibonacci", 1)
    p2, a2 = build("fibonacci", 99)
    assert a1.achieved_ligands == 40
    sites = [sorted(s["surface_atom_index"] for s in a.site_to_ligand) for a in (a1, a2)]
    assert sites[0] == sites[1]  # the same sites whatever the seed (only azimuths are seeded)
    pr, _ = build("distributed_random", 1)
    assert _nn_cv(p1) < _nn_cv(pr)


@pytest.mark.parametrize("miller", [(1, 1, 1), (2, 1, 1), (1, 0, 0)])
def test_rectangular_slab_cell_is_orthogonal_and_keeps_density(miller):
    base = GeometryConfig(shape="slab", miller=miller, slab_length_nm=3, slab_width_nm=3, layers=4)
    rect = GeometryConfig(shape="slab", miller=miller, slab_length_nm=3, slab_width_nm=3, layers=4,
                          box_shape="rectangular")
    cov = CoverageConfig(target=0)
    a0 = run_pipeline(RunConfig(geometry=base, coverage=cov)).assembly
    a1 = run_pipeline(RunConfig(geometry=rect, coverage=cov)).assembly
    cell = np.asarray(a1.cell_A)
    assert np.allclose(cell - np.diag(np.diag(cell)), 0, atol=1e-6)
    area0 = np.linalg.norm(np.cross(a0.cell_A[0], a0.cell_A[1]))
    area1 = cell[0, 0] * cell[1, 1]
    assert abs(a0.n_atoms / area0 - a1.n_atoms / area1) < 1e-6 * a0.n_atoms / area0


def test_cube_boxes():
    a = run_pipeline(RunConfig(geometry=GeometryConfig(diameter_nm=2.5, box_shape="cube"),
                               coverage=CoverageConfig(target=0))).assembly
    assert np.allclose(np.diag(a.cell_A), a.cell_A[0][0])
    s = run_pipeline(RunConfig(geometry=GeometryConfig(shape="slab", miller=(1, 0, 0), slab_length_nm=6,
                                                       slab_width_nm=6, layers=4, vacuum_A=15,
                                                       box_shape="cube"),
                               coverage=CoverageConfig(target=0))).assembly
    assert np.allclose(np.diag(s.cell_A), s.cell_A[0][0])


def test_min_box_lengths():
    a = run_pipeline(RunConfig(geometry=GeometryConfig(diameter_nm=2.0, box_lengths_nm=(12, 13, 14)),
                               coverage=CoverageConfig(target=0))).assembly
    assert np.allclose(np.diag(a.cell_A), [120, 130, 140])


@pytest.mark.parametrize("ctype,count,miller,layers", [("slit", 2, (1, 1, 1), 4),
                                                         ("slit", 3, (1, 1, 1), 3),
                                                         ("pore", 2, (1, 0, 0), 14)])
def test_channels_place_ligands_on_inner_walls(ctype, count, miller, layers):
    g = GeometryConfig(shape="channel", channel_type=ctype, channel_count=count, channel_width_nm=1.4,
                       miller=miller, slab_length_nm=3, slab_width_nm=6, layers=layers)
    r = run_pipeline(RunConfig(geometry=g, coverage=CoverageConfig(target=8, region="inner",
                                                                   distribution="fibonacci")))
    assert r.report.ok and r.assembly.achieved_ligands == 8
    assert {s["face"] for s in r.assembly.site_to_ligand} == {"inner"}
    assert r.substrate.channel_info["type"] == ctype


def test_carbon_inner_region_in_nanotube():
    cfg = RunConfig(material=MaterialConfig(source="carbon", preset_id=None),
                    carbon=CarbonConfig(kind="nanotube", cnt_n=20, cnt_m=20, cnt_length_nm=3),
                    ligand=LigandConfig(smiles="C", anchor_atom_index=0, binding_template="C_basal"),
                    coverage=CoverageConfig(target=6, region="inner", distribution="fibonacci"))
    r = run_pipeline(cfg)
    assert r.report.ok and r.assembly.achieved_ligands == 6
    a = r.assembly
    centre = np.asarray(a.positions_A)[:, :2].mean(axis=0)
    radius = np.linalg.norm(np.asarray(a.positions_A)[:, :2] - centre, axis=1)
    sub = np.array([role == "substrate" for role in a.chemical_role])
    lig = ~sub & (np.asarray(a.elements) == "C")
    assert radius[lig].max() < radius[sub & (np.asarray(a.elements) == "C")].max()  # inside


def test_hand_edits_rebuild_and_are_flagged(tmp_path):
    r = run_pipeline(RunConfig(coverage=CoverageConfig(target=10, distribution="fibonacci")))
    r2 = remove_ligands(r, [2, 5])
    assert r2.assembly.achieved_ligands == 8 and r2.report.ok
    assert r2.assembly.n_atoms == r.assembly.n_atoms - 2 * (r.assembly.n_atoms - r.assembly.n_substrate) // 10
    used = {p["surface_atom_index"] for p in r2.assembly.placements}
    free = next(int(i) for i in r2.sites.surface_atom_index if int(i) not in used)
    r3 = add_ligand_at(r2, free)
    assert r3.assembly.achieved_ligands == 9 and r3.report.ok
    edits = r3.assembly.provenance["manual_edits"]
    assert [e["action"] for e in edits] == ["remove_ligands", "add_ligand"]
    pkg = export_package(r3, tmp_path)
    assert "MANUAL EDITS" in (pkg.stage_dir / "README.txt").read_text()
    with pytest.raises(BuilderError):
        remove_ligands(r3, [99])
    with pytest.raises(BuilderError):
        add_ligand_at(r3, r3.assembly.n_atoms - 1)  # a ligand atom, not a surface atom


# --- round 2: library, chemisorption, user box, pores, clean cuts, session, certificate --------


def test_preloaded_library_loads_offline():
    from nanobuilder import mp_library
    from nanobuilder.crystal import build_crystal
    entries = mp_library.entries()
    assert len(entries) >= 25
    e = next(v for v in entries.values() if v["formula"] == "TiO2" and v["spacegroup"] == "I4_1/amd")
    cfg = MaterialConfig(source="materials_project", preset_id=None, mp_id=e["material_id"])
    crystal = build_crystal(cfg)  # no bytes given: read from the library
    assert crystal.atoms.get_chemical_formula(empirical=True) in ("O2Ti", "TiO2")
    assert "anatase" in crystal.structure_label
    with pytest.raises(BuilderError, match="not in the preloaded"):
        build_crystal(MaterialConfig(source="materials_project", preset_id=None, mp_id="mp-0"))


@pytest.mark.parametrize("preset_id,smiles,element", [("Fe_bcc", "NCCCCCC", "N"),
                                                       ("Ag_fcc", "SCCO", "S")])
def test_chemisorption_bonds_anchor_to_surface(preset_id, smiles, element):
    cfg = RunConfig(material=MaterialConfig(source="preset", preset_id=preset_id),
                    geometry=GeometryConfig(diameter_nm=2.5),
                    ligand=LigandConfig(smiles=smiles, anchor_atom_index=0,
                                        binding_template="chemisorption"),
                    coverage=CoverageConfig(target=6, distribution="fibonacci"))
    r = run_pipeline(cfg)
    a = r.assembly
    assert r.report.ok and a.achieved_ligands == 6
    surface = [b for b in a.bond_list if "substrate" in (a.chemical_role[b[0]], a.chemical_role[b[1]])]
    assert len(surface) == 6 and all(element in b[2] for b in surface)
    if element == "S":  # thiol H removed on chemisorption
        assert r.ligand.formula == "C2H5OS"


def test_physisorption_has_no_surface_bond():
    from nanobuilder.chemistry import TEMPLATES
    assert TEMPLATES["geometric"].label.startswith("Physisorption")
    r = run_pipeline(RunConfig(material=MaterialConfig(source="preset", preset_id="Fe_bcc"),
                               geometry=GeometryConfig(diameter_nm=2.5),
                               ligand=LigandConfig(smiles="NCCCC", anchor_atom_index=0,
                                                   binding_template="geometric", distance_A=2.4),
                               coverage=CoverageConfig(target=4)))
    a = r.assembly
    assert not [b for b in a.bond_list if "substrate" in (a.chemical_role[b[0]], a.chemical_role[b[1]])]


def test_user_box_sets_slab_size_and_fit_is_checked():
    r = run_pipeline(RunConfig(geometry=GeometryConfig(shape="slab", miller=(1, 1, 1), layers=4,
                                                       box_shape="rectangular",
                                                       box_lengths_nm=(6, 5, 7)),
                               coverage=CoverageConfig(target=0)))
    cell = np.asarray(r.assembly.cell_A)
    assert abs(cell[0, 0] - 60) < 6 and abs(cell[1, 1] - 50) < 6 and abs(cell[2, 2] - 70) < 1e-6
    with pytest.raises(BuilderError, match="enlarge the box"):
        run_pipeline(RunConfig(geometry=GeometryConfig(diameter_nm=3, box_shape="cube",
                                                       box_lengths_nm=(3, 3, 3))))


def test_periodic_cut_is_clean():
    """After centring, no atom sits in the outermost 0.4 A of the periodic x and y edges."""
    r = run_pipeline(RunConfig(geometry=GeometryConfig(shape="slab", miller=(2, 1, 1), layers=4,
                                                       box_shape="rectangular",
                                                       box_lengths_nm=(5, 5, 6)),
                               coverage=CoverageConfig(target=0)))
    a = r.assembly
    frac = np.linalg.solve(np.asarray(a.cell_A).T, np.asarray(a.positions_A).T).T % 1.0
    lx, ly = np.linalg.norm(a.cell_A[0]), np.linalg.norm(a.cell_A[1])
    gap_x = min(frac[:, 0].min(), 1 - frac[:, 0].max()) * lx
    gap_y = min(frac[:, 1].min(), 1 - frac[:, 1].max()) * ly
    assert gap_x > 0.2 and gap_y > 0.2


def test_pores_spread_over_the_surface_with_margins():
    from nanobuilder.shapes import PORE_MARGIN_SPACINGS
    g = GeometryConfig(shape="channel", channel_type="pore", channel_count=4, channel_width_nm=1.2,
                       miller=(1, 0, 0), layers=8, box_shape="rectangular", box_lengths_nm=(8, 8, 6))
    r = run_pipeline(RunConfig(geometry=g, coverage=CoverageConfig(target=8, region="inner",
                                                                   distribution="fibonacci")))
    centers = np.array(r.substrate.channel_info["centers_A"])
    assert len(centers) == 4
    assert np.ptp(centers[:, 0]) > 15 and np.ptp(centers[:, 1]) > 15  # not one row
    d = np.linalg.norm(centers[:, None] - centers[None], axis=2) + np.eye(4) * 1e9
    assert d.min() >= 12 + PORE_MARGIN_SPACINGS * 2.8
    too_big = GeometryConfig(shape="channel", channel_type="pore", channel_count=1,
                             channel_width_nm=4.0, miller=(1, 0, 0), layers=8,
                             box_shape="rectangular", box_lengths_nm=(4, 4, 6))
    with pytest.raises(BuilderError, match="pore needs the box"):
        run_pipeline(RunConfig(geometry=too_big))


def test_session_text_restores_the_configuration():
    import ast
    from nanobuilder.session import config_from_settings, session_text
    cfg = RunConfig(geometry=GeometryConfig(shape="rod", diameter_nm=2.0, rod_length_nm=5.0,
                                            box_shape="cube", box_lengths_nm=(9, 9, 9)),
                    coverage=CoverageConfig(target=8, distribution="fibonacci"))
    r = run_pipeline(cfg)
    r = remove_ligands(r, [3])
    text = session_text(cfg, r.assembly)
    calls = [n.value for n in ast.parse(text).body]
    assert [c.func.id for c in calls] == ["settings", "build", "remove_ligand"]
    restored, _ = config_from_settings(ast.literal_eval(calls[0].args[0]))
    assert restored.to_dict() == cfg.to_dict()


def test_certificate_is_a_pdf_with_details(tmp_path):
    from nanobuilder.certificate import certificate_pdf, describe, sections
    r = run_pipeline(RunConfig(coverage=CoverageConfig(target=6, distribution="fibonacci")))
    pdf = certificate_pdf(r.config, r.assembly, r.report, r.crystal.atoms, files=["x"])
    assert pdf[:4] == b"%PDF" and len(pdf) > 10000
    assert "sphere" in describe(r.config, r.assembly)
    headings = [h for h, _ in sections(r.config, r.assembly, r.report, r.crystal.atoms)]
    assert headings[:4] == ["Structure and box", "Crystal", "Atoms", "Ligand"]
    pkg = export_package(r, tmp_path)
    assert (pkg.stage_dir / "certificate_of_build.pdf").is_file()
    assert (pkg.stage_dir / "session.py").is_file()
