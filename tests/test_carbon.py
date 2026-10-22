"""Carbon nanomaterials: lattices, HierGO-style defects and oxidation, covalent ligands, exports."""

from collections import Counter, deque

import numpy as np
import pytest
from ase import Atoms

from nanobuilder.carbon import build_carbon_substrate, composition_metrics
from nanobuilder.carbon_defects import find_hole_library, flower_defects, ruptures, stone_wales
from nanobuilder.carbon_lattice import (A_CC, build_carbon_sphere, build_graphene, build_nanotube,
                                        cnt_diameter, fullerene_coordinates, wall_chiralities)
from nanobuilder.carbon_oxidation import AGED_PROFILES, group_fractions
from nanobuilder.pipeline import export_package, run_pipeline, validate_inputs
from nanobuilder.schema import (BuilderError, CarbonConfig, CoverageConfig, ExportConfig,
                                GeometryConfig, LigandConfig, MaterialConfig, RunConfig)


def carbon_run(carbon, ligand=None, coverage=None, export=None):
    return RunConfig(material=MaterialConfig(source="carbon", preset_id=None), carbon=carbon,
                     ligand=ligand or LigandConfig(), coverage=coverage or CoverageConfig(target=0),
                     export=export or ExportConfig(), geometry=GeometryConfig(vacuum_A=20.0))


def substrate(**kw):
    return build_carbon_substrate(CarbonConfig(**kw), GeometryConfig(vacuum_A=20.0))[1]


def degrees(sub):
    deg = Counter()
    for i, j, _ in sub.bonds:
        deg[i] += 1
        deg[j] += 1
    return deg


def ring_sizes(net, max_ring=9):
    rings = set()
    for v in net.lattice_atoms():
        nbs = net.c_nb(v)
        for x in range(len(nbs)):
            for y in range(x + 1, len(nbs)):
                a, b = nbs[x], nbs[y]
                prev, queue, found = {a: None}, deque([a]), False
                while queue and not found:
                    c = queue.popleft()
                    for d in net.c_nb(c):
                        if d == v or d in prev:
                            continue
                        prev[d] = c
                        if d == b:
                            found = True
                            break
                        queue.append(d)
                if found:
                    path = [b]
                    while prev[path[-1]] is not None:
                        path.append(prev[path[-1]])
                    if len(path) + 1 <= max_ring:
                        rings.add(frozenset(path + [v]))
    return Counter(len(r) for r in rings)


def bond_lengths(sub):
    return np.array([sub.atoms.get_distance(i, j, mic=True) for i, j, _ in sub.bonds])


# --------------------------------------------------------------------------
# lattices
# --------------------------------------------------------------------------


def test_periodic_graphene_is_perfect_honeycomb():
    sub = substrate(kind="graphene", periodic=True, length_nm=4, width_nm=4)
    n = len(sub.atoms)
    assert n % 4 == 0 and len(sub.bonds) == 3 * n // 2
    assert set(degrees(sub).values()) == {3}
    assert np.allclose(bond_lengths(sub), A_CC, atol=1e-6)
    assert sub.pbc == (True, True, False) and sub.shape == "sheet"
    assert all(o == 1.5 for *_, o in sub.bonds)


def test_graphene_flake_edges_are_hydrogen_terminated():
    sub = substrate(kind="graphene", length_nm=4, width_nm=3)
    el = sub.atoms.get_chemical_symbols()
    deg = degrees(sub)
    assert all(deg[i] == 3 for i, e in enumerate(el) if e == "C")
    assert all(deg[i] == 1 for i, e in enumerate(el) if e == "H")
    assert el.count("H") == sum(1 for r in sub.labels if r == "C_edge") > 0
    hpos = sub.atoms.get_positions()[[i for i, e in enumerate(el) if e == "H"]]
    d = np.linalg.norm(hpos[:, None] - hpos[None], axis=2) + np.eye(len(hpos)) * 9
    assert d.min() > 1.5
    assert sub.pbc == (False, False, False)


def test_bilayer_ab_stacking_and_spacing():
    net = build_graphene(CarbonConfig(periodic=True, layers=2, length_nm=3, width_nm=3))
    z = sorted(set(np.round(net.pos[:, 2], 3)))
    assert len(z) == 2 and abs(z[1] - z[0] - 3.35) < 1e-6
    a = net.pos[[i for i in range(net.n) if net.layer[i] == 0], :2]
    b = net.pos[[i for i in range(net.n) if net.layer[i] == 1], :2]
    # Bernal: half of the top-layer atoms sit above bottom-layer atoms
    dist = np.min(np.linalg.norm(b[:, None] - a[None], axis=2), axis=1)
    assert np.isclose(np.mean(dist < 0.1), 0.5)


def test_nanotube_geometry_and_multiwall_spacing():
    net = build_nanotube(CarbonConfig(kind="nanotube", cnt_n=10, cnt_m=10, periodic=True))
    r = np.hypot(net.pos[:, 0], net.pos[:, 1])
    assert np.allclose(r, cnt_diameter(10, 10) / 2, atol=1e-3)
    assert abs(cnt_diameter(10, 10) - 13.56) < 0.01
    assert all(len(net.c_nb(i)) == 3 for i in range(net.n))
    walls = wall_chiralities(5, 5, 3, 3.4, periodic=True)
    assert walls == [(5, 5), (10, 10), (15, 15)]
    gaps = np.diff([cnt_diameter(*w) / 2 for w in walls])
    assert np.all(np.abs(gaps - 3.4) < 0.1)
    with pytest.raises(BuilderError):
        wall_chiralities(10, 5, 2, 3.4, periodic=True)


def test_finite_nanotube_rims_and_regions():
    sub = substrate(kind="nanotube", cnt_n=8, cnt_m=8, cnt_length_nm=3)
    assert sub.shape == "nanotube" and sub.pbc == (False,) * 3
    assert len(sub.faces["tips"].atom_indices) > 0
    assert "H" in sub.composition


def test_c60_and_carbon_onion():
    pos, bonds = fullerene_coordinates(1)
    assert len(pos) == 60 and len(bonds) == 90
    lengths = np.linalg.norm(pos[[i for i, _ in bonds]] - pos[[j for _, j in bonds]], axis=1)
    assert np.all(np.abs(lengths - A_CC) < 0.03)
    assert abs(np.linalg.norm(pos, axis=1).mean() - 3.52) < 0.05
    net = build_carbon_sphere(CarbonConfig(kind="carbon_sphere", sphere_type="fullerene"))
    assert ring_sizes(net) == Counter({5: 12, 6: 20})
    onion = build_carbon_sphere(CarbonConfig(kind="carbon_sphere", sphere_type="onion",
                                             onion_shells=3))
    assert onion.n == 60 + 240 + 540
    assert all(3.3 < g < 3.7 for g in onion.geometry["shell_gaps_A"])


# --------------------------------------------------------------------------
# defects
# --------------------------------------------------------------------------


@pytest.mark.parametrize("shapes", ["random_growth", "hiergo_library"])
@pytest.mark.parametrize("termination", ["H", "none"])
def test_vacancy_fraction_matches_request(shapes, termination):
    if shapes == "hiergo_library" and find_hole_library() is None:
        pytest.skip("HierGO HoleLibrary not present")
    sub = substrate(kind="graphene", periodic=True, length_nm=6, width_nm=6, vacancy_percent=10,
                    hole_shapes=shapes, edge_termination=termination)
    vac = sub.carbon_info["composition"]["vacancy_percent_actual"]
    # holes are added one at a time and the last one is kept only if it gets closer
    # to the target; library shapes are coarser (and grow more when H-capped)
    assert (8.5 <= vac <= 12.5) if shapes == "hiergo_library" else (9.0 <= vac <= 11.5), vac
    assert "crowded_edge_atoms_removed" not in sub.carbon_info["defects"]
    assert sub.carbon_info["defects"]["holes"]["holes"] >= 2


def test_stone_wales_and_24_atom_rotations_make_5_7_rings():
    rng = np.random.default_rng(3)
    net = build_graphene(CarbonConfig(periodic=True, length_nm=5, width_nm=5))
    assert stone_wales(net, 3, rng, set()) == 3
    rings = ring_sizes(net)
    assert rings[5] == 6 and rings[7] == 6
    net = build_graphene(CarbonConfig(periodic=True, length_nm=5, width_nm=5))
    assert flower_defects(net, 2, rng, set()) == 2
    rings = ring_sizes(net)
    assert rings[5] == 12 and rings[7] == 12
    lengths = [net.dist(i, j) for i, j in net.order]
    assert 1.3 < min(lengths) and max(lengths) < 1.6
    assert all(len(net.c_nb(i)) == 3 for i in net.lattice_atoms())


def test_rupture_bridges_broken_bonds_with_ether_oxygens():
    net = build_graphene(CarbonConfig(periodic=True, length_nm=5, width_nm=5))
    out = ruptures(net, 1, 5, np.random.default_rng(1), set())
    assert out == {"ruptures": 1, "ether_oxygens": 5, "length_bonds": 5}
    oxygens = [i for i in range(net.n) if net.el[i] == "O"]
    zs = []
    for o in oxygens:
        c1, c2 = net.adj[o]
        assert abs(net.dist(o, c1) - 1.37) < 0.02 and abs(net.dist(o, c2) - 1.37) < 0.02
        v1, v2 = net.delta(o, c1), net.delta(o, c2)
        ang = np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)))
        assert 100 < ang < 125
        zs.append(np.sign(net.pos[o, 2] - net.pos[c1, 2]))
    assert set(zs) == {-1.0, 1.0}  # alternating above and below


def test_pores_and_defect_reporting():
    sub = substrate(kind="graphene", periodic=True, length_nm=6, width_nm=6, pore_count=2,
                    pore_diameter_nm=1.2, stone_wales_percent=2)
    d = sub.carbon_info["defects"]
    assert d["pores"]["pores"] == 2 and d["pores"]["removed_atoms"] > 30
    assert d["stone_wales"]["placed"] == d["stone_wales"]["requested"] > 0


# --------------------------------------------------------------------------
# oxidation
# --------------------------------------------------------------------------


def test_fresh_go_matches_the_amber_go_dataset_metrics():
    """GO30 of Fedai et al. 2026: C:O 3.33:1 and 40.0 % modified carbons (2:1 OH:epoxy)."""
    sub = substrate(kind="graphene", length_nm=5, width_nm=5, oxygen_percent=30,
                    hydroxyl_epoxy_ratio="2:1")
    c = sub.carbon_info["composition"]
    assert abs(c["oxygen_percent_O_per_C"] - 30) < 0.3
    assert abs(c["C_to_O_ratio"] - 3.33) < 0.05
    assert abs(c["modified_carbon_percent"] - 40.0) < 1.0
    o_roles = Counter(r for r in sub.labels if r.startswith("O_"))
    assert abs(o_roles["O_hydroxyl"] / o_roles["O_epoxy"] - 2.0) < 0.2
    assert abs(sub.carbon_info["oxidation"]["achieved_hydroxyl_to_epoxy"] - 2.0) < 0.2
    assert set(sub.carbon_info["groups"]) == {"hydroxyl", "paired_hydroxyl", "epoxy"}


def test_high_oxidation_reaches_its_target():
    sub = substrate(kind="graphene", length_nm=5, width_nm=5, oxygen_percent=60)
    assert sub.carbon_info["composition"]["oxygen_percent_O_per_C"] > 58.5


def test_aged_profiles_add_edge_groups():
    for key in AGED_PROFILES:
        f, _ = group_fractions(CarbonConfig(oxidation_profile=key))
        assert abs(sum(f.values()) - 1) < 1e-9
    sub = substrate(kind="graphene", length_nm=6, width_nm=6, vacancy_percent=8,
                    oxygen_percent=20, oxidation_profile="aged_1200K")
    groups = sub.carbon_info["groups"]
    for g in ("carboxyl", "carbonyl", "ether", "epoxy"):
        assert groups.get(g, 0) > 0, g
    roles = Counter(sub.labels)
    assert roles["C_carboxyl"] == roles["O_carboxyl"] == roles["O_carboxyl_H"] == roles["H_carboxyl"]


def test_faces_balanced_or_restricted():
    both = substrate(kind="graphene", length_nm=5, width_nm=5, oxygen_percent=20)
    bal = both.carbon_info["oxidation"]["face_balance"]
    assert abs(bal["top_or_outer"] - bal["bottom"]) <= 3
    top = substrate(kind="graphene", periodic=True, length_nm=5, width_nm=5, oxygen_percent=10,
                    oxidation_faces="top")
    pos = top.atoms.get_positions()
    zc = np.mean([pos[i, 2] for i, r in enumerate(top.labels) if r == "C_sp2"])
    assert all(pos[i, 2] > zc for i, e in enumerate(top.atoms.get_chemical_symbols()) if e == "O")


def test_tiles_and_domains_make_heterogeneous_oxidation():
    tiles = substrate(kind="graphene", periodic=True, length_nm=8, width_nm=8,
                      oxidation_pattern="tiles", tile_grid=(2, 2), tile_oxygen_min=0,
                      tile_oxygen_max=60, seed=5)
    per_tile = tiles.carbon_info["oxidation"]["tiles"]["oxygen_percent"]
    assert len(per_tile) == 4 and max(per_tile) - min(per_tile) > 5
    dom = substrate(kind="graphene", periodic=True, length_nm=8, width_nm=8, oxygen_percent=5,
                    oxidation_pattern="domains", domain_count=1, domain_diameter_nm=3.0)
    centre = np.array(dom.carbon_info["oxidation"]["domains"]["centres_A"][0])
    pos = dom.atoms.get_positions()
    cell = np.array(dom.atoms.cell)
    inside = 0
    o_atoms = [i for i, e in enumerate(dom.atoms.get_chemical_symbols()) if e == "O"]
    for i in o_atoms:
        d = pos[i, :2] - centre
        d -= np.round(d / np.diag(cell)[:2]) * np.diag(cell)[:2]
        inside += np.linalg.norm(d) <= 15.0 + 2.0
    assert inside / len(o_atoms) > 0.9


def test_curved_walls_are_oxidized_outward():
    sub = substrate(kind="nanotube", cnt_n=10, cnt_m=10, periodic=True, oxygen_percent=10)
    pos = sub.atoms.get_positions() - sub.center_A
    r = np.hypot(pos[:, 0], pos[:, 1])
    radius = cnt_diameter(10, 10) / 2
    o = [i for i, e in enumerate(sub.atoms.get_chemical_symbols()) if e == "O"]
    assert o and np.all(r[o] > radius + 0.8)
    assert abs(sub.carbon_info["composition"]["oxygen_percent_O_per_C"] - 10) < 0.5


def test_oxidized_structures_pass_chemistry_checks():
    for kw in (dict(oxygen_percent=30, oxidation_profile="aged_900K", vacancy_percent=8,
                    stone_wales_percent=3, rupture_count=1),
               dict(oxygen_percent=25, oxidation_profile="aged_2500K", vacancy_percent=12,
                    hole_size="small"),
               dict(oxygen_percent=40, pore_count=1, pore_diameter_nm=1.2, flower_percent=4)):
        result = run_pipeline(carbon_run(CarbonConfig(length_nm=6, width_nm=6, **kw)))
        assert result.report.ok, result.report.errors
        a = result.assembly
        at = Atoms(a.elements, positions=a.positions_A, cell=a.cell_A, pbc=a.pbc)
        lengths = [at.get_distance(i, j, mic=True) for i, j, _ in a.bond_list]
        assert 0.95 < min(lengths) and max(lengths) < 1.7


def test_composition_metrics_definitions():
    m = composition_metrics(["C"] * 10 + ["O"] * 2 + ["H"] * 2,
                            ["C_sp2"] * 7 + ["C_hydroxyl"] * 2 + ["C_epoxy"] + ["O_x"] * 2 + ["H"] * 2)
    assert m["oxygen_percent_O_per_C"] == 20.0 and m["C_to_O_ratio"] == 5.0
    assert m["modified_carbon_percent"] == 30.0


# --------------------------------------------------------------------------
# covalent ligands
# --------------------------------------------------------------------------


def _links(a):
    return [(s["surface_atom_index"], s["anchor_atom_index"]) for s in a.site_to_ligand]


def test_basal_addition_makes_sp3_carbon_and_c_c_bond():
    result = run_pipeline(carbon_run(
        CarbonConfig(periodic=True, length_nm=4, width_nm=4),
        LigandConfig(smiles="c1ccccc1", anchor_atom_index=0, binding_template="C_basal"),
        CoverageConfig(target=6, region="top")))
    a = result.assembly
    assert a.achieved_ligands == 6 and result.report.ok
    assert result.ligand.formula == "C6H5"
    at = Atoms(a.elements, positions=a.positions_A, cell=a.cell_A, pbc=a.pbc)
    deg = Counter()
    for i, j, _ in a.bond_list:
        deg[i] += 1
        deg[j] += 1
    for s, n in _links(a):
        assert a.type_labels[s] == "C_linked" and deg[s] == 4
        assert abs(at.get_distance(s, n, mic=True) - 1.53) < 1e-3
        assert a.positions_A[n, 2] > a.positions_A[s, 2]  # top face
    assert any(set(t[:3]) >= {s, n} for t in a.angle_list for s, n in _links(a)[:1])
    assert a.suggested_types["gaff2"][_links(a)[0][0]] == "c3"


def test_edge_substitution_replaces_edge_hydrogens():
    cfg = CarbonConfig(length_nm=4, width_nm=4)
    bare = run_pipeline(carbon_run(cfg)).assembly
    result = run_pipeline(carbon_run(
        cfg, LigandConfig(smiles="NCCO", anchor_atom_index=0, binding_template="C_edge"),
        CoverageConfig(target=5, region="edges")))
    a = result.assembly
    assert a.achieved_ligands == 5
    n_h_sub = sum(1 for e, c in zip(a.elements, a.component) if e == "H" and c == "substrate")
    assert n_h_sub == bare.elements.count("H") - 5
    at = Atoms(a.elements, positions=a.positions_A, cell=a.cell_A, pbc=a.pbc)
    for s, n in _links(a):
        assert a.type_labels[s] == "C_aryl_link"
        assert abs(at.get_distance(s, n) - 1.40) < 1e-3  # aryl C-N


def test_amide_coupling_consumes_carboxyl_oh():
    cfg = CarbonConfig(length_nm=6, width_nm=6, vacancy_percent=8, oxygen_percent=20,
                       oxidation_profile="aged_1200K", seed=11)
    result = run_pipeline(carbon_run(
        cfg, LigandConfig(smiles="NCCN", anchor_atom_index=0, binding_template="C_carboxyl"),
        CoverageConfig(target=4, region="edges")))
    a = result.assembly
    assert a.achieved_ligands == 4 and result.report.ok
    before = a.provenance["carbon"]["composition"]["O"]
    after = a.provenance["carbon"]["composition_after_attachment"]["O"]
    assert before - after == 4
    gaff = a.suggested_types["gaff2"]
    for s, n in _links(a):
        assert a.type_labels[s] == "C_acyl" and gaff[s] == "c" and gaff[n] == "n"


def test_epoxide_opening_leaves_a_hydroxyl_anti_to_the_ligand():
    cfg = CarbonConfig(periodic=True, length_nm=5, width_nm=5, oxygen_percent=15,
                       hydroxyl_epoxy_ratio="1:1", seed=3)
    result = run_pipeline(carbon_run(
        cfg, LigandConfig(smiles="NCCO", anchor_atom_index=0, binding_template="C_epoxide"),
        CoverageConfig(target=4, region="top")))
    a = result.assembly
    assert a.achieved_ligands == 4 and result.report.ok
    comp0 = a.provenance["carbon"]["composition"]
    comp1 = a.provenance["carbon"]["composition_after_attachment"]
    assert comp1["O"] == comp0["O"] and comp1["H"] == comp0["H"] + 4
    nb = {}
    for i, j, _ in a.bond_list:
        nb.setdefault(i, set()).add(j)
        nb.setdefault(j, set()).add(i)
    zc = np.mean([a.positions_A[i, 2] for i, r in enumerate(a.type_labels) if r == "C_sp2"])
    for s, n in _links(a):
        assert a.type_labels[s] == "C_linked"
        oh_carbons = [c for c in nb[s] if a.type_labels[c] == "C_hydroxyl"]
        assert oh_carbons
        # ligand above the sheet; the opened epoxide's hydroxyl stays below (anti addition)
        assert a.positions_A[n, 2] > zc + 1.0
        assert any(a.positions_A[x, 2] < zc - 1.0 for c in oh_carbons for x in nb[c]
                   if a.elements[x] == "O")


def test_c60_amine_and_template_errors():
    result = run_pipeline(carbon_run(
        CarbonConfig(kind="carbon_sphere"),
        LigandConfig(smiles="NCCCC", anchor_atom_index=0, binding_template="C_basal"),
        CoverageConfig(target=6, region="entire_surface", distribution="evenly_spaced")))
    assert result.assembly.achieved_ligands == 6
    bad = validate_inputs(carbon_run(
        CarbonConfig(), LigandConfig(smiles="CC(C)(C)C", anchor_atom_index=1,
                                     binding_template="C_basal"), CoverageConfig(target=2)))
    assert not bad["ok"] and "no hydrogen" in bad["errors"][0]
    bad = validate_inputs(carbon_run(
        CarbonConfig(), LigandConfig(smiles="CCO", anchor_atom_index=0,
                                     binding_template="C_carboxyl"), CoverageConfig(target=2)))
    assert not bad["ok"] and "anchor atom of N, O, S" in bad["errors"][0]
    metal = RunConfig(ligand=LigandConfig(smiles="NCC", binding_template="C_basal"))
    assert "carbon nanomaterial" in validate_inputs(metal)["errors"][0]
    au_on_c = validate_inputs(carbon_run(CarbonConfig(), LigandConfig(), CoverageConfig(target=2)))
    assert not au_on_c["ok"]
    with pytest.raises(BuilderError, match="0 valid sites"):
        run_pipeline(carbon_run(
            CarbonConfig(periodic=True),
            LigandConfig(smiles="NCC", anchor_atom_index=0, binding_template="C_edge"),
            CoverageConfig(target=2, region="edges")))


# --------------------------------------------------------------------------
# config, exports
# --------------------------------------------------------------------------


def test_carbon_config_round_trip_and_validation():
    cfg = carbon_run(CarbonConfig(kind="nanotube", cnt_n=12, cnt_m=3, oxygen_percent=5,
                                  tile_grid=(2, 3), oxidation_profile="custom"))
    again = RunConfig.from_dict(cfg.to_dict())
    assert again.to_dict() == cfg.to_dict()
    assert again.carbon.tile_grid == (2, 3)
    for bad in (dict(kind="nanotube", cnt_n=3, cnt_m=5), dict(oxygen_percent=120),
                dict(hydroxyl_epoxy_ratio="2:x"), dict(layers=3, oxidize_all_layers=True,
                                                       oxygen_percent=10),
                dict(oxidation_profile="custom", custom_fractions={"hydroxyl": 0.0})):
        with pytest.raises(BuilderError):
            CarbonConfig(**bad).validate()
    big = carbon_run(CarbonConfig(length_nm=60, width_nm=60))
    assert not validate_inputs(big)["ok"]


def test_mol2_pdb_and_lammps_exports(tmp_path):
    result = run_pipeline(carbon_run(
        CarbonConfig(length_nm=4, width_nm=4, oxygen_percent=20),
        LigandConfig(smiles="NCCO", anchor_atom_index=0, binding_template="C_epoxide"),
        CoverageConfig(target=2, region="both")))
    package = export_package(result, tmp_path)
    assert all(r.round_trip_passed for r in package.reports)
    from nanobuilder.exporter_mol2 import read_mol2
    mol2 = read_mol2(package.stage_dir / "mol2" / "structure.mol2")
    types = Counter(a["type"] for a in mol2["atoms"])
    for t in ("ca", "c3", "oh", "ho", "os", "ha", "n3"):
        assert types[t] > 0, t
    res = Counter(a["resname"] for a in mol2["atoms"])
    assert {"GRA", "OH", "EPX", "LIG"} <= set(res)
    assert any(code == "ar" for *_, code in mol2["bonds"])
    pdb = (package.stage_dir / "pdb" / "structure.pdb").read_text()
    assert " GRA A" in pdb and " EPX A" in pdb and "CONECT" in pdb
    assert max(len(line) for line in pdb.splitlines()) <= 80
    readme = (package.stage_dir / "README.txt").read_text()
    assert "Oxygen content 100 N_O/N_C" in readme

    single = carbon_run(CarbonConfig(kind="carbon_sphere"), export=ExportConfig(residue_scheme="single"))
    package = export_package(run_pipeline(single), tmp_path / "single")
    mol2 = read_mol2(package.stage_dir / "mol2" / "structure.mol2")
    assert {a["resname"] for a in mol2["atoms"]} == {"GRO"}

    tube = run_pipeline(carbon_run(CarbonConfig(kind="nanotube", periodic=True)))
    package = export_package(tube, tmp_path / "tube")
    data = (package.stage_dir / "lammps" / "structure.data").read_text()
    template = (package.stage_dir / "lammps" / "in.lammps.template").read_text()
    assert "f f p" in template and "dihedrals" in data
    assert tube.assembly.pbc == (False, False, True)
