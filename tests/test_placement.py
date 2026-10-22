"""Test matrix 12-15 and 20: counts, density, clash filter, periodic images, seeds."""

import numpy as np
import pytest
from ase import Atoms
from ase.neighborlist import neighbor_list

from nanobuilder.placement import orient, place_ligands, target_count, vdw_radius
from nanobuilder.schema import AttachmentConfig, CoverageConfig
from nanobuilder.sites import find_surface_sites, selected_area_nm2


def run(substrate, ligand, **cov):
    coverage = CoverageConfig(**cov)
    sites = find_surface_sites(substrate, AttachmentConfig(region=coverage.region, surface_species=["Au"],
                                                           distance_A=ligand.distance_A))
    return place_ligands(substrate, ligand, sites, coverage, selected_area_nm2(substrate, coverage.region))


def min_nonbonded(assembly):
    atoms = Atoms(assembly.elements, positions=assembly.positions_A, cell=assembly.cell_A,
                  pbc=assembly.pbc)
    bonded = {tuple(sorted(b[:2])) for b in assembly.bond_list}
    bonded |= {tuple(sorted((a[0], a[2]))) for a in assembly.angle_list}
    i, j, d = neighbor_list("ijd", atoms, 4.0)
    best = np.inf
    for a, b, dist in zip(i, j, d):
        if a < b and (a, b) not in bonded and "ligand" in (assembly.component[a], assembly.component[b]):
            best = min(best, dist)
    return best


def test_exact_count_on_sphere(au_sphere, scco_thiolate):
    a = run(au_sphere, scco_thiolate, mode="count", target=12)
    assert a.achieved_ligands == a.target_ligands == 12
    assert a.achieved_density_nm2 == pytest.approx(12 / 28.274, abs=1e-3)  # 0.424
    assert a.n_atoms == len(au_sphere.atoms) + 12 * 9


def test_density_uses_actual_slab_area(au111_slab, scco_thiolate):
    a = run(au111_slab, scco_thiolate, mode="density", density_nm2=0.5, region="top")
    area = au111_slab.faces["top"].area_nm2
    assert a.target_ligands == round(0.5 * area)
    assert a.achieved_ligands == a.target_ligands
    assert a.surface_area_nm2 == pytest.approx(area)
    sub_z = au111_slab.atoms.positions[:, 2].max()
    lig = np.array([c == "ligand" for c in a.component])
    assert np.all(a.positions_A[lig, 2] > a.positions_A[~lig, 2].max())  # top only
    assert sub_z  # substrate unchanged


def test_shortfall_is_reported_not_hidden(au_sphere, scco_thiolate):
    a = run(au_sphere, scco_thiolate, mode="count", target=400, max_attempts=3000)
    assert a.target_ligands == 400
    assert a.achieved_ligands < 400
    assert any(f"Only {a.achieved_ligands} of 400" in w for w in a.warnings)


def test_no_clashes_and_attachment_kept(au_sphere, scco_thiolate):
    a = run(au_sphere, scco_thiolate, mode="count", target=40, clash_alpha=0.7)
    pairs = {(s["surface_atom_index"], s["anchor_atom_index"]) for s in a.site_to_ligand}
    assert len(pairs) == 40
    for sub, anchor in pairs:
        assert np.linalg.norm(a.positions_A[sub] - a.positions_A[anchor]) == pytest.approx(2.4, abs=1e-6)
        assert (sub, anchor, "Au-S") in a.bond_list
    # Every remaining nonbonded contact respects 0.7 x (vdW_i + vdW_j); H...H is the tightest case.
    assert min_nonbonded(a) >= 0.7 * 2 * vdw_radius("H") - 1e-6


def test_periodic_images_are_checked(au111_slab, scco_thiolate):
    a = run(au111_slab, scco_thiolate, mode="density", density_nm2=3.0, region="top")
    assert min_nonbonded(a) >= 0.7 * 2 * vdw_radius("H") - 1e-6  # MIC distances via ase pbc


def test_seed_reproducibility(au_sphere, scco_thiolate):
    a = run(au_sphere, scco_thiolate, mode="count", target=12, seed=5)
    b = run(au_sphere, scco_thiolate, mode="count", target=12, seed=5)
    c = run(au_sphere, scco_thiolate, mode="count", target=12, seed=6)
    assert np.array_equal(a.positions_A, b.positions_A)
    assert not np.array_equal(a.positions_A, c.positions_A)


def test_evenly_spaced_spreads_ligands(au_sphere, scco_thiolate):
    a = run(au_sphere, scco_thiolate, mode="count", target=8, distribution="evenly_spaced")
    anchors = a.positions_A[[s["anchor_atom_index"] for s in a.site_to_ligand]]
    d = np.linalg.norm(anchors[:, None] - anchors[None], axis=2)
    np.fill_diagonal(d, np.inf)
    assert d.min() > 10.0  # 8 ligands on a 3 nm sphere sit far apart


def test_min_separation(au_sphere, scco_thiolate):
    a = run(au_sphere, scco_thiolate, mode="count", target=10, min_separation_A=12.0)
    sites = [s["surface_atom_index"] for s in a.site_to_ligand]
    p = a.positions_A[sites]
    d = np.linalg.norm(p[:, None] - p[None], axis=2)
    np.fill_diagonal(d, np.inf)
    assert d.min() >= 12.0


def test_target_count_rounding():
    assert target_count(CoverageConfig(mode="density", density_nm2=0.5), 28.83) == 14
    assert target_count(CoverageConfig(mode="count", target=12), 1.0) == 12


def test_orient_places_anchor_on_normal():
    xyz = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]])
    out = orient(xyz, np.array([1.0, 0, 0]), np.array([5.0, 5, 5]), np.array([0, 0, 1.0]), 2.4, 0.3)
    assert np.allclose(out[0], [5, 5, 7.4])
    assert np.allclose(out[2] - out[0], [0, 0, 2.0])


def test_batched_clash_tree_matches_rebuild_every_ligand(monkeypatch):
    """The recent-ligand buffer must not change any placement."""
    import copy
    import numpy as np
    from nanobuilder import placement
    from nanobuilder.examples import EXAMPLES
    from nanobuilder.pipeline import run_pipeline
    from nanobuilder.schema import RunConfig
    cfg = copy.deepcopy(EXAMPLES["C_au111_slab_density"]["config"])
    cfg["ligand"].update(smiles="SCCCCCCCCCCCC", anchor_atom_index=0, conformer_mode="extended")
    cfg["coverage"]["density_nm2"] = 3.0
    batched = run_pipeline(RunConfig.from_dict(cfg)).assembly
    monkeypatch.setattr(placement, "TREE_REBUILD_EVERY", 1)
    every = run_pipeline(RunConfig.from_dict(cfg)).assembly
    assert batched.achieved_ligands == every.achieved_ligands > 40
    assert np.allclose(batched.positions_A, every.positions_A)
