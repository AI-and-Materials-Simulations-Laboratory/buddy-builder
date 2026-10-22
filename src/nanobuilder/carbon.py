"""Carbon nanomaterial substrates: build, defect, oxidize, terminate, and expose sites.

``build_carbon_substrate`` is the carbon counterpart of ``crystal.build_crystal`` +
``shapes.build_geometry``: it returns a pseudo :class:`CrystalModel` (for
provenance and the normalized unit cell in the export) and a
:class:`SubstrateModel` whose ``bonds``, ``labels``, ``residues`` and
``carbon_sites`` drive covalent ligand attachment, topology and typing.

Pipeline (HierGO order, generalized to nanotubes and fullerenes):
lattice -> holes / pores -> trim dangling carbons -> Stone-Wales and 24-atom
rotations -> ruptures -> oxygen groups -> edge hydrogens -> residues/labels.
"""

from __future__ import annotations

from math import pi

import numpy as np
from ase import Atoms

from . import carbon_defects as defects
from .carbon_lattice import (A_CC, A_LAT, AREA_PER_C, B_LAT, C_H, EDGE_RIBBON_NM, CarbonNetwork,
                             build_carbon_sphere, build_graphene, build_nanotube,
                             clear_crowded_edges, cnt_diameter, prune)
from .carbon_oxidation import PYRAMIDALIZE_A, decorate, needs_ethers
from .schema import (BuilderError, CarbonConfig, CrystalModel, FaceInfo, GeometryConfig,
                     SiteCatalog, SubstrateModel)

SP3_ROLES = {"C_hydroxyl", "C_epoxy", "C_linked"}
LINK_PUCKER_A = 0.40  # outward displacement of a basal carbon that bonds a ligand (sp3)
SP2_LATTICE_ROLES = {"C_sp2", "C_edge", "C_ether", "C_ring_cooh", "C_phenol", "C_carbonyl",
                     "C_aryl_link", "C_dangling"}

# Residue names: GRA / OH / EPX follow the AMBER GO dataset (Fedai et al. 2026);
# the others extend the scheme to HierGO's aged groups.
GROUP_RESIDUES = {"hydroxyl": "OH", "paired_hydroxyl": "OH", "edge_hydroxyl": "OH",
                  "epoxy": "EPX", "ether": "ETH", "carbonyl": "CBO", "carboxyl": "CBX"}

# Suggested force-field atom types for substrate roles. GAFF2: minimal mapping of the AMBER GO
# dataset (ca, c3, oh, os, ho) extended to edge groups. OPLS-AA: HierGO's Oplsatom2 table
# (phenol types from the OPLS-AA phenol entries). Suggestions only; no parameters are assigned.
ROLE_TYPES = {
    #  role            GAFF2   OPLS-AA (HierGO)
    "C_sp2": ("ca", "opls_147"),
    "C_edge": ("ca", "opls_145"),
    "H_edge": ("ha", "opls_146"),
    "C_dangling": ("ca", None),
    "C_hydroxyl": ("c3", "opls_159"),
    "O_hydroxyl": ("oh", "opls_154"),
    "H_hydroxyl": ("ho", "opls_155"),
    "C_epoxy": ("c3", "opls_184"),
    "O_epoxy": ("os", "opls_180"),
    "C_ether": ("ca", "opls_518"),
    "O_ether": ("os", "opls_180"),
    "C_carbonyl": ("c", "opls_280"),
    "O_carbonyl": ("o", "opls_281"),
    "C_ring_cooh": ("ca", "opls_147"),
    "C_carboxyl": ("c", "opls_267"),
    "O_carboxyl": ("o", "opls_269"),
    "O_carboxyl_H": ("oh", "opls_268"),
    "H_carboxyl": ("ho", "opls_270"),
    "C_phenol": ("ca", "opls_166"),
    "O_phenol": ("oh", "opls_167"),
    "H_phenol": ("ho", "opls_168"),
    # created by covalent ligand attachment
    "C_linked": ("c3", None),
    "C_aryl_link": ("ca", None),
    "C_acyl": ("c", None),
    "O_acyl": ("o", None),
}


def lattice_bond_order(role_a: str, role_b: str, order: float) -> float:
    """Aromatic (1.5) between sp2 framework carbons, single once either carbon is sp3."""
    if role_a in SP3_ROLES or role_b in SP3_ROLES:
        return 1.0
    if role_a in SP2_LATTICE_ROLES and role_b in SP2_LATTICE_ROLES:
        return 1.5
    return order


# --------------------------------------------------------------------------
# size estimate (checked before building)
# --------------------------------------------------------------------------


def estimate_atoms(cfg: CarbonConfig) -> int:
    if cfg.kind == "graphene":
        n = cfg.length_nm * cfg.width_nm * 100.0 / AREA_PER_C * cfg.layers
    elif cfg.kind == "nanotube":
        from .carbon_lattice import wall_chiralities
        try:
            walls = wall_chiralities(cfg.cnt_n, cfg.cnt_m, cfg.cnt_walls, cfg.wall_spacing_A, False)
        except BuilderError:
            walls = [(cfg.cnt_n, cfg.cnt_m)] * cfg.cnt_walls
        n = sum(pi * cnt_diameter(*w) * cfg.cnt_length_nm * 10.0 / AREA_PER_C for w in walls)
    else:
        ks = [cfg.fullerene_order] if cfg.sphere_type == "fullerene" else range(1, cfg.onion_shells + 1)
        n = sum(60 * k * k for k in ks)
    o = n * cfg.oxygen_percent / 100.0
    return int(n + 2 * o)  # O plus roughly one H per O


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------


def build_carbon_substrate(cfg: CarbonConfig, geometry: GeometryConfig):
    """Return (CrystalModel, SubstrateModel) for a carbon nanomaterial."""
    cfg.validate()
    est = estimate_atoms(cfg)
    if est > geometry.max_atoms:
        raise BuilderError(f"This carbon structure would have about {est} atoms; the limit is "
                           f"{geometry.max_atoms}. Make it smaller or raise Max atoms")
    rng = np.random.default_rng(int(cfg.seed))
    if cfg.kind == "graphene":
        net = build_graphene(cfg)
    elif cfg.kind == "nanotube":
        net = build_nanotube(cfg)
    else:
        net = build_carbon_sphere(cfg)
    h_term = cfg.edge_termination == "H"
    ethers = needs_ethers(cfg)
    blocked = set()
    info = {"defects": {}}
    finite = not any(net.pbc) and net.kind != "carbon_sphere"

    if finite:  # clean flake/tube outline first; not counted as vacancies
        info["defects"]["edge_trim_atoms"] = prune(net, cfg.keep_hanging_carbons)
        if h_term:
            info["defects"]["edge_trim_atoms"] += len(clear_crowded_edges(
                net, cfg.keep_hanging_carbons, ethers))
    n_pristine = len(net.lattice_atoms())
    info["pristine_lattice_carbons"] = n_pristine
    if cfg.vacancy_percent > 0:
        info["defects"]["holes"] = defects.make_holes(net, cfg, rng, blocked, h_term, ethers)
    if int(cfg.pore_count) > 0:
        if net.kind != "sheet":
            net.warnings.append("Circular pores are only cut in graphene sheets; ignored")
        else:
            info["defects"]["pores"] = defects.make_pores(net, int(cfg.pore_count),
                                                          float(cfg.pore_diameter_nm), rng, blocked)
            prune(net, cfg.keep_hanging_carbons)
            if h_term:
                clear_crowded_edges(net, cfg.keep_hanging_carbons, ethers)
    n_now = len(net.lattice_atoms())
    if cfg.stone_wales_percent > 0:
        if net.kind == "carbon_sphere":
            net.warnings.append("Stone-Wales defects are not generated on fullerene shells; ignored")
        else:
            want = int(round(cfg.stone_wales_percent / 100.0 * n_now / 16.0))
            got = defects.stone_wales(net, want, rng, blocked)
            info["defects"]["stone_wales"] = {"requested": want, "placed": got,
                                              "convention": "HierGO: one defect per 16 atoms per %"}
            if got < want:
                net.warnings.append(f"Only {got} of {want} Stone-Wales defects found room")
    if cfg.flower_percent > 0:
        if net.kind != "sheet":
            net.warnings.append("24-atom rotations are only made in flat graphene sheets; ignored")
        else:
            want = int(round(cfg.flower_percent / 100.0 * n_now / 41.0))
            got = defects.flower_defects(net, want, rng, blocked)
            info["defects"]["rotations_24_atom"] = {"requested": want, "placed": got,
                                                    "convention": "HierGO: one per 41 atoms per %"}
            if got < want:
                net.warnings.append(f"Only {got} of {want} 24-atom rotations found room")
    if int(cfg.rupture_count) > 0:
        if net.kind != "sheet":
            net.warnings.append("Ruptures are only made in flat graphene sheets; ignored")
        else:
            info["defects"]["ruptures"] = defects.ruptures(net, int(cfg.rupture_count),
                                                           int(cfg.rupture_length), rng, blocked)
    info["defects"]["final_trim_atoms"] = prune(net, cfg.keep_hanging_carbons)

    ox_report = None
    if cfg.oxygen_percent > 0 or cfg.oxidation_pattern == "tiles":
        ox_report = decorate(net, cfg, rng)
        net.warnings += ox_report.pop("warnings", [])

    if h_term:
        removed = clear_crowded_edges(net, cfg.keep_hanging_carbons, False)
        n_h, more = terminate_edges(net, cfg.keep_hanging_carbons)
        if removed or more:
            info["defects"]["crowded_edge_atoms_removed"] = len(removed) + more
        info["edge_hydrogens"] = n_h
    else:
        bare = [a for a in net.lattice_atoms() if len(net.adj[a]) == 2]
        for a in bare:
            net.role[a] = "C_dangling"
        if bare:
            net.warnings.append(f"{len(bare)} edge carbons left unterminated (dangling bonds); most "
                                "force fields need hydrogen-terminated edges")
    return _finish(net, cfg, geometry, info, ox_report, n_pristine)


def terminate_edges(net: CarbonNetwork, keep_hanging: bool = False) -> tuple:
    """Add H to framework carbons that miss a bond (bare sp2 edges and edge sp3 carbons).

    Any carbon whose new H would come within 1.5 A of another H, or 1.6 A of a
    non-bonded heavy atom, is removed instead (with its group) and the
    neighbourhood re-trimmed; this repeats until every H fits. Returns
    (hydrogens added, carbons removed).
    """
    from .carbon_lattice import kill_with_groups, make_tree, prune_local, tree_query_points

    removed = 0
    for _ in range(40):
        need = []
        for a in net.lattice_atoms():
            want = 4 if net.role[a] in SP3_ROLES else 3
            if len(net.adj[a]) < want and len(net.c_nb(a)) >= 2:
                need.append(a)
        if not need:
            return 0, removed
        hpos = np.array([net.pos[a] + C_H * net.edge_direction(a) for a in need])
        victims = set()
        htree = make_tree(hpos, net.cell, net.pbc)
        for x, y in sorted(htree.query_pairs(1.5)):
            if need[x] not in victims and need[y] not in victims:
                victims.add(need[y])
        alive = [i for i in range(net.n) if net.alive[i]]
        tree = net.tree(alive)
        for k, near in enumerate(tree.query_ball_point(tree_query_points(tree, hpos), 1.6)):
            a = need[k]
            if a in victims:
                continue
            for m in near:
                b = alive[m]
                if b == a or b in net.adj[a]:
                    continue
                limit = 1.5 if net.el[b] == "H" else 1.6
                if np.linalg.norm(net.mic(net.pos[b] - hpos[k])) < limit:
                    victims.add(a)
                    break
        if not victims:
            for a, h in zip(need, hpos):
                i = net.add_atom("H", net.wrap(h), "H_edge", net.layer[a], net.normal[a])
                net.add_bond(a, i, 1.0)
                if net.role[a] == "C_sp2":
                    net.role[a] = "C_edge"
            return len(need), removed
        seeds = set()
        for v in victims:
            seeds |= set(net.c_nb(v))
            kill_with_groups(net, v)
            removed += 1
        removed += len(prune_local(net, seeds, keep_hanging))
    raise BuilderError("Could not terminate the edges with hydrogen without overlaps")


# --------------------------------------------------------------------------
# finish: compact, order by residue, metrics, faces, sites
# --------------------------------------------------------------------------


def _residue_order(net: CarbonNetwork):
    """Atom order grouped by residue: each framework carbon with its edge H, then each group."""
    order, names = [], []
    resid = []
    r = 0
    for a in range(net.n):
        if not (net.alive[a] and net.lattice[a]):
            continue
        r += 1
        order.append(a)
        names.append("GRA")
        resid.append(r)
        for h in sorted(net.adj[a]):
            if net.alive[h] and net.role[h] == "H_edge":
                order.append(h)
                names.append("GRA")
                resid.append(r)
    for g in net.groups:
        atoms = [a for a in g["atoms"] if net.alive[a]]
        if not atoms:
            continue
        r += 1
        for a in atoms:
            order.append(a)
            names.append(GROUP_RESIDUES.get(g["kind"], "GRP"))
            resid.append(r)
    seen = set(order)
    for a in range(net.n):  # anything not covered (should not happen)
        if net.alive[a] and a not in seen:
            r += 1
            order.append(a)
            names.append("GRP")
            resid.append(r)
    return order, names, resid


def composition_metrics(elements, roles) -> dict:
    """Composition and the validation metrics of the AMBER GO dataset."""
    n_c = sum(1 for e in elements if e == "C")
    n_o = sum(1 for e in elements if e == "O")
    n_h = sum(1 for e in elements if e == "H")
    lattice = sum(1 for r in roles if r in SP3_ROLES or r in SP2_LATTICE_ROLES)
    sp3 = sum(1 for r in roles if r in SP3_ROLES)
    return {
        "C": n_c, "O": n_o, "H": n_h,
        "oxygen_percent_O_per_C": round(100.0 * n_o / n_c, 3) if n_c else 0.0,
        "C_to_O_ratio": round(n_c / n_o, 3) if n_o else None,
        "framework_carbons": lattice,
        "sp3_carbons": sp3,
        "modified_carbon_percent": round(100.0 * sp3 / lattice, 3) if lattice else 0.0,
    }


def _finish(net, cfg, geometry, info, ox_report, n_pristine):
    order, resnames, resids = _residue_order(net)
    index = {old: new for new, old in enumerate(order)}
    pos = net.pos[order].copy()
    elements = [net.el[a] for a in order]
    roles = [net.role[a] for a in order]
    normals = net.normal[order].copy()
    layers = [net.layer[a] for a in order]
    bonds = []
    for (i, j), o in net.order.items():
        if net.alive[i] and net.alive[j]:
            a, b = sorted((index[i], index[j]))
            if net.lattice[i] and net.lattice[j]:
                o = lattice_bond_order(net.role[i], net.role[j], o)
            bonds.append((a, b, o))
    bonds.sort()
    vac = geometry.vacuum_A

    # Box: periodic axes keep their lattice length; finite axes get vacuum padding.
    cell = np.array(net.cell, dtype=float)
    lo, hi = pos.min(axis=0), pos.max(axis=0)
    shift = np.zeros(3)
    for k in range(3):
        if net.pbc[k]:
            continue
        shift[k] = vac - lo[k]
        cell[k] = 0.0
        cell[k, k] = hi[k] - lo[k] + 2 * vac
    pos += shift
    for k in range(3):
        if net.pbc[k]:
            pos[:, k] %= cell[k, k]
    atoms = Atoms(elements, positions=pos, cell=cell, pbc=net.pbc)
    center = None
    if net.kind in ("nanotube", "carbon_sphere"):
        center = (np.zeros(3) + shift) if net.kind == "carbon_sphere" else np.array(
            [shift[0], shift[1], cell[2, 2] / 2])

    metrics = composition_metrics(elements, roles)
    lattice_final = metrics["framework_carbons"]
    metrics["vacancy_percent_actual"] = round(100.0 * (1 - lattice_final / n_pristine), 3) \
        if n_pristine else 0.0
    groups = {}
    for g in net.groups:
        if any(net.alive[a] for a in g["atoms"]):
            groups[g["kind"]] = groups.get(g["kind"], 0) + 1
    info.update({"kind": cfg.kind, "summary": cfg.summary(), "geometry": net.geometry,
                 "composition": metrics, "groups": groups, "oxidation": ox_report,
                 "pyramidalization_A": PYRAMIDALIZE_A if cfg.pyramidalize else 0.0,
                 "periodic": list(net.pbc), "edge_termination": cfg.edge_termination,
                 "seed": int(cfg.seed)})

    faces, sites = _faces_and_sites(net, cfg, order, index, pos, normals, layers, roles, elements,
                                    bonds)
    warnings = list(dict.fromkeys(net.warnings))
    if any(net.pbc) and cfg.kind in ("graphene", "nanotube"):
        warnings.append("Periodic carbon framework: C-C bonds cross the cell boundary. LAMMPS, "
                        "GROMACS and HOOMD handle this; AMBER/tleap does not (build a finite "
                        "flake or tube for AMBER)")
    shape = {"graphene": "sheet", "nanotube": "nanotube", "carbon_sphere": "carbon_sphere"}[cfg.kind]
    substrate = SubstrateModel(
        atoms=atoms, crystal=None, shape=shape, pbc=net.pbc, faces=faces,
        composition={e: elements.count(e) for e in sorted(set(elements))},
        center_A=center, vacuum_A=vac, warnings=warnings, bonds=bonds, labels=roles,
        residues=list(zip(resnames, resids)), normals=normals, carbon_sites=sites,
        carbon_info=info)
    if cfg.kind == "nanotube":
        substrate.rod_axis = np.array([0.0, 0.0, 1.0])
        substrate.radius_A = max(net.geometry["diameters_A"]) / 2
    if cfg.kind == "carbon_sphere":
        substrate.radius_A = max(net.geometry["shell_radii_A"])
    crystal = carbon_crystal(cfg, net)
    substrate.crystal = crystal
    return crystal, substrate


def carbon_crystal(cfg: CarbonConfig, net: CarbonNetwork) -> CrystalModel:
    """Provenance record standing in for a bulk crystal; atoms = the repeat unit used."""
    if cfg.kind == "graphene":
        from .carbon_lattice import graphene_layer
        unit = Atoms("C4", positions=graphene_layer(1, 1), cell=[A_LAT, B_LAT, 3.35 * 2],
                     pbc=True)
        label, proto, sg = cfg.summary(), "graphene (orthogonal 4-atom cell)", "P6/mmm"
    elif cfg.kind == "nanotube":
        from ase.build import nanotube
        unit = nanotube(int(cfg.cnt_n), int(cfg.cnt_m), length=1, bond=A_CC)
        d = cnt_diameter(int(cfg.cnt_n), int(cfg.cnt_m))
        unit.cell[0, 0] = unit.cell[1, 1] = d + 10.0
        unit.positions += [unit.cell[0, 0] / 2, unit.cell[1, 1] / 2, 0]
        unit.pbc = True
        label, proto, sg = cfg.summary(), f"({cfg.cnt_n},{cfg.cnt_m}) nanotube unit cell", None
    else:
        from .carbon_lattice import fullerene_coordinates
        p, _ = fullerene_coordinates(1)
        unit = Atoms("C60", positions=p + 8.0, cell=[16.0] * 3, pbc=True)
        label, proto, sg = cfg.summary(), "C60 (Ih); larger shells are Goldberg (k,k) polyhedra", "Ih"
    return CrystalModel(
        atoms=unit, source_type="carbon", structure_label=label, space_group=sg,
        prototype=proto, anchor_species_candidates=["C"],
        lattice_source="ideal sp2 carbon, a_CC = 1.42 A (generated)",
        normalization="generated carbon framework; HierGO-style defects and oxidation",
        validation_warnings=[])


def _faces_and_sites(net, cfg, order, index, pos, normals, layers, roles, elements, bonds):
    """Face areas (for coverage) and attachment-site candidates by kind."""
    n = len(order)
    adj = [[] for _ in range(n)]
    for i, j, _ in bonds:
        adj[i].append(j)
        adj[j].append(i)
    lattice = [r in SP3_ROLES or r in SP2_LATTICE_ROLES for r in roles]
    faces_sites = {}
    for layer, faces in net.faces_sites.items():
        faces_sites[layer] = faces

    basal, edge, carboxyl, epoxide = [], [], [], []
    for a in range(n):
        if not lattice[a]:
            continue
        cn = [b for b in adj[a] if lattice[b]]
        if roles[a] == "C_sp2" and len(cn) == 3 and len(adj[a]) == 3:
            for f in faces_sites.get(layers[a], ()):
                basal.append((a, f))
        if roles[a] == "C_edge":
            hs = [b for b in adj[a] if roles[b] == "H_edge"]
            if hs and layers[a] in faces_sites:
                edge.append((a, hs[0]))
    for a in range(n):
        if roles[a] == "C_carboxyl" and layers[a] in faces_sites:
            oh = [b for b in adj[a] if roles[b] == "O_carboxyl_H"]
            if oh:
                hs = [b for b in adj[oh[0]] if roles[b] == "H_carboxyl"]
                ring = [b for b in adj[a] if lattice[b]]
                if hs and ring:
                    carboxyl.append((a, oh[0], hs[0], ring[0]))
        if roles[a] == "O_epoxy" and layers[a] in faces_sites:
            cs = [b for b in adj[a] if lattice[b]]
            if len(cs) == 2:
                rel = pos[a] - pos[cs[0]]
                for k in range(3):
                    if net.pbc[k]:
                        rel[k] -= net.cell[k, k] * round(rel[k] / net.cell[k, k])
                face = 1 if np.dot(rel, normals[cs[0]]) > 0 else -1
                epoxide.append((a, cs[0], cs[1], face))

    # Channel faces: between the layers of a stacked sheet, inside the innermost tube wall.
    basal_inner = []
    if net.kind in ("sheet", "nanotube"):
        layer_set = sorted(set(layers[a] for a in range(n) if lattice[a]))
        for a in range(n):
            if not (lattice[a] and roles[a] == "C_sp2" and len(adj[a]) == 3
                    and len([b for b in adj[a] if lattice[b]]) == 3):
                continue
            if net.kind == "sheet" and len(layer_set) > 1:
                if layers[a] < layer_set[-1]:
                    basal_inner.append((a, 1))
                if layers[a] > layer_set[0]:
                    basal_inner.append((a, -1))
            elif net.kind == "nanotube" and layers[a] == layer_set[0]:
                basal_inner.append((a, -1))

    faces = {}
    if net.kind == "sheet":
        layer_ids = sorted(set(layers[a] for a in range(n) if lattice[a]))
        top, bottom = max(layer_ids), min(layer_ids)

        def layer_area(layer):
            if net.pbc[0] and net.pbc[1]:
                return float(net.cell[0, 0] * net.cell[1, 1]) / 100.0
            return sum(1 for a in range(n) if lattice[a] and layers[a] == layer) * AREA_PER_C / 100.0

        edge_atoms = [a for a in range(n) if lattice[a] and len([b for b in adj[a] if lattice[b]]) < 3]
        faces["top"] = FaceInfo(np.array([a for a in range(n) if lattice[a] and layers[a] == top]),
                                ["C"], "graphene top face", layer_area(top))
        faces["bottom"] = FaceInfo(np.array([a for a in range(n) if lattice[a] and layers[a] == bottom]),
                                   ["C"], "graphene bottom face", layer_area(bottom))
        faces["edges"] = FaceInfo(np.array(edge_atoms, dtype=int), ["C"],
                                  "sheet and hole edges (ribbon area)",
                                  len(edge_atoms) * 0.23 * EDGE_RIBBON_NM)
        if len(layer_ids) > 1:
            faces["inner"] = FaceInfo(np.array(sorted({a for a, _ in basal_inner}), dtype=int),
                                      ["C"], "faces between the layers (slit channels)",
                                      2 * sum(layer_area(l_) for l_ in layer_ids[:-1]))
    elif net.kind == "nanotube":
        outer = max(net.faces_sites)
        d = net.geometry["diameters_A"][outer] / 10.0
        length = net.geometry["length_A"] / 10.0
        rims = [a for a in range(n) if lattice[a] and layers[a] == outer
                and len([b for b in adj[a] if lattice[b]]) < 3]
        faces["side"] = FaceInfo(np.array([a for a in range(n) if lattice[a] and layers[a] == outer]),
                                 ["C"], "outer wall", pi * d * length)
        faces["tips"] = FaceInfo(np.array(rims, dtype=int), ["C"], "open rims (ribbon area)",
                                 0.0 if any(net.pbc) else 2 * pi * d * EDGE_RIBBON_NM)
        inner_d = net.geometry["diameters_A"][0] / 10.0
        faces["inner"] = FaceInfo(np.array(sorted({a for a, _ in basal_inner}), dtype=int),
                                  ["C"], "inner wall (channel)", pi * inner_d * length)
    else:
        outer = max(net.faces_sites)
        r = net.geometry["shell_radii_A"][outer] / 10.0
        faces["entire_surface"] = FaceInfo(
            np.array([a for a in range(n) if lattice[a] and layers[a] == outer]), ["C"],
            "outer shell", 4 * pi * r * r)
    sites = {"basal": basal, "edge": edge, "carboxyl": carboxyl, "epoxide": epoxide,
             "basal_inner": basal_inner}
    return faces, sites


# --------------------------------------------------------------------------
# attachment sites for covalent ligands
# --------------------------------------------------------------------------


def region_kinds(shape: str, region: str) -> dict:
    """Which site groups a coverage region selects: faces (signs), edges, inner channel faces."""
    if shape == "sheet":
        table = {"top": ({1}, False), "bottom": ({-1}, False), "both": ({1, -1}, False),
                 "edges": (set(), True), "entire_surface": ({1, -1}, True), "inner": (set(), False)}
    elif shape == "nanotube":
        table = {"side": ({1}, False), "tips": (set(), True), "entire_surface": ({1}, True),
                 "inner": (set(), False)}
    else:
        table = {"entire_surface": ({1}, True)}
    if region not in table:
        names = {"sheet": "Top, Bottom, Both faces, Edges, Entire surface or Between layers",
                 "nanotube": "Entire surface, Side (outer wall), Tips (open rims) or Inner wall",
                 "carbon_sphere": "Entire surface"}[shape]
        raise BuilderError(f"This carbon substrate has the regions {names}")
    faces, edges = table[region]
    return {"faces": faces, "edges": edges, "inner": region == "inner"}


def region_face_names(shape: str, region: str) -> list:
    if region == "inner":
        return ["inner"]
    if shape == "sheet":
        return {"top": ["top"], "bottom": ["bottom"], "both": ["top", "bottom"],
                "edges": ["edges"], "entire_surface": ["top", "bottom", "edges"]}[region]
    if shape == "nanotube":
        return {"side": ["side"], "tips": ["tips"], "entire_surface": ["side", "tips"]}[region]
    return ["entire_surface"]


def _tangent_away(substrate, a, b):
    """Unit vector from atom b to atom a, projected onto a's tangent plane."""
    v = substrate.atoms.get_distance(b, a, mic=True, vector=True)
    n = substrate.normals[a]
    v = v - np.dot(v, n) * n
    return v / max(np.linalg.norm(v), 1e-9)


def carbon_site_catalog(substrate: SubstrateModel, template_key: str, region: str) -> SiteCatalog:
    """Attachment sites on a carbon substrate for one template and coverage region.

    Each site carries the substrate edits its ligand makes when placed:
    - C_basal: sp2 C -> sp3 C (moved 0.40 A outward), new C-X bond along the normal;
    - C_edge: the edge H is replaced by the ligand (in-plane C-X bond);
    - C_carboxyl: COOH loses its OH; the acyl carbon bonds the ligand (amide/ester);
    - C_epoxide: the ring opens; O stays on one carbon as OH (H added), the ligand
      bonds the other carbon from the opposite face (anti, flat sheets) or the
      same face, tilted apart (curved walls).
    ``geometric`` uses the basal carbons (and edge carbons) without any bond.
    """
    from .chemistry import TEMPLATES

    sel = region_kinds(substrate.shape, region)
    s = substrate.carbon_sites
    pos = substrate.atoms.get_positions()
    nrm = substrate.normals
    delta = substrate.carbon_info.get("pyramidalization_A", 0.0)  # GO groups (0.25 A or 0)
    link = LINK_PUCKER_A if delta > 0 else 0.0  # ligand attachment carbons
    rows = []
    template = TEMPLATES[template_key]
    kind = template.site_kind if template_key != "geometric" else "geometric"
    face_name = {1: "top", -1: "bottom"} if substrate.shape == "sheet" else {1: "side", -1: "side"}
    edge_name = "edges" if substrate.shape == "sheet" else "tips"
    if substrate.shape == "carbon_sphere":
        face_name = {1: "entire_surface", -1: "entire_surface"}
        edge_name = "entire_surface"

    if kind in ("basal", "geometric"):
        for a, f in s.get("basal", []):
            if f not in sel["faces"]:
                continue
            n = f * nrm[a]
            if kind == "geometric":
                rows.append((pos[a], a, n, face_name[f], {}, {a}))
                continue
            moved = pos[a] + link * n
            rows.append((moved, a, n, face_name[f],
                         {"move": {a: moved}, "relabel": {a: "C_linked"}}, {a}))
    if sel.get("inner"):
        if "inner" not in substrate.faces:
            raise BuilderError("This carbon substrate has no channel faces: stack two or more "
                               "graphene layers (slit channels) or use a nanotube (inner wall)")
        if kind not in ("basal", "geometric"):
            raise BuilderError("Inside channels only basal-plane addition (C_basal) or geometric "
                               "adsorption are available")
        for a, f in s.get("basal_inner", []):
            n = f * nrm[a]
            if kind == "geometric":
                rows.append((pos[a], a, n, "inner", {}, {a}))
                continue
            moved = pos[a] + link * n
            rows.append((moved, a, n, "inner", {"move": {a: moved}, "relabel": {a: "C_linked"}},
                         {a}))
    if kind == "edge" or (kind == "geometric" and sel["edges"]):
        if sel["edges"]:
            for c, h in s.get("edge", []):
                d = substrate.atoms.get_distance(c, h, mic=True, vector=True)
                d /= np.linalg.norm(d)
                if kind == "geometric":
                    continue
                rows.append((pos[c], c, d, edge_name,
                             {"remove": [h], "relabel": {c: "C_aryl_link"}, "ignore": [h]}, {c, h}))
    if kind == "carboxyl" and sel["edges"]:
        for cc, oh, h, ring in s.get("carboxyl", []):
            d = substrate.atoms.get_distance(cc, oh, mic=True, vector=True)
            d /= np.linalg.norm(d)
            o1 = [b for b in _neighbors(substrate, cc) if substrate.labels[b] == "O_carboxyl"]
            relabel = {cc: "C_acyl"}
            if o1:
                relabel[o1[0]] = "O_acyl"
            rows.append((pos[cc], cc, d, edge_name,
                         {"remove": [oh, h], "relabel": relabel, "ignore": [oh, h]}, {cc, oh, h}))
    if kind == "epoxide":
        for o, c1, c2, f in s.get("epoxide", []):
            for ca, cb in ((c1, c2), (c2, c1)):  # O stays on ca, ligand on cb
                site = _epoxide_site(substrate, o, ca, cb, f, delta, link, sel, face_name)
                if site is not None:
                    rows.append(site)
    if not rows:
        what = {"basal": "unmodified sp2 basal carbons", "edge": "hydrogen-terminated edge carbons",
                "carboxyl": "carboxylic acid groups (oxidize with an aged or custom profile)",
                "epoxide": "epoxide groups (oxidize the substrate first)",
                "geometric": "basal carbons"}[kind]
        raise BuilderError(f"0 valid sites: this substrate has no {what} in the selected region "
                           f"({region})")
    positions, idx, normals, faces, edits, claims = zip(*rows)
    area = {f: substrate.faces[f].area_nm2 for f in region_face_names(substrate.shape, region)}
    warn = ["Covalent attachment is a starting geometry: bond lengths come from typical values, "
            "not from an optimization of the interface"] if kind != "geometric" else []
    return SiteCatalog(positions_A=np.array(positions), surface_atom_index=np.array(idx, dtype=int),
                       normals=np.array(normals), face=list(faces), element=["C"] * len(rows),
                       site_type=[kind] * len(rows), area_nm2_by_face=area, warnings=warn,
                       edits=list(edits), claims=[set(c) for c in claims])


def _neighbors(substrate, a):
    out = []
    for i, j, _ in substrate.bonds:
        if i == a:
            out.append(j)
        elif j == a:
            out.append(i)
    return out


def _epoxide_site(substrate, o, ca, cb, f, delta, link, sel, face_name):
    """Site where the ligand opens epoxide ``o`` at carbon ``cb`` (O becomes OH on ``ca``)."""
    pos = substrate.atoms.get_positions()
    nrm = substrate.normals
    flat = substrate.shape == "sheet"
    lig_face = -f if flat else f
    if lig_face not in sel["faces"]:
        return None
    u_a = _tangent_away(substrate, ca, cb)  # in-plane, from cb toward ca
    if flat:
        n_lig = lig_face * nrm[cb]
        cb_new = pos[cb] + (delta + link) * n_lig  # undo the epoxide puckering, pucker the other way
        o_dir = f * nrm[ca]
    else:
        n_lig = nrm[cb] - 0.45 * u_a
        n_lig /= np.linalg.norm(n_lig)
        cb_new = pos[cb] + (link - delta) * nrm[cb]
        o_dir = nrm[ca] + 0.45 * u_a
        o_dir /= np.linalg.norm(o_dir)
    o_new = pos[ca] + 1.43 * o_dir
    h_new = o_new + 0.97 * (np.cos(np.radians(71.5)) * o_dir + np.sin(np.radians(71.5)) * u_a)
    edits = {"move": {cb: cb_new, o: o_new}, "unbond": [(o, cb)],
             "relabel": {o: "O_hydroxyl", ca: "C_hydroxyl", cb: "C_linked"},
             "residue": {o: "OH"},
             "add": [{"element": "H", "position": h_new, "label": "H_hydroxyl", "bond_to": o,
                      "order": 1.0, "residue_of": o}],
             "ignore": [o]}
    return (cb_new, cb, n_lig, face_name[lig_face], edits, {o, ca, cb})
