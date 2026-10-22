"""Lattice defects for carbon frameworks, re-implemented from the HierGO workflow.

- Vacancy holes grown by random neighbour selection (HierGO create-tile.py),
  with HierGO's size distributions (small 1-5, medium 5-10, large 10-30,
  mixed 1-30 atoms) or shapes taken from HierGO's ``HoleLibrary``.
- Circular pores (nanoporous graphene).
- Stone-Wales defects: one C-C bond rotated by 90 degrees (55-77 rings).
- 24-atom rotations (HierGO "mini grain boundary"): a coronene-like C24 unit
  rotated by 30 degrees (90 degrees modulo the 60-degree lattice symmetry),
  i.e. the "flower" defect.
- Ruptures (HierGO create-rupture.py): a row of parallel C-C bonds broken and
  bridged by ether oxygens alternating above and below the sheet.

Unlike the original scripts, every defect rewires the bond graph explicitly
and is followed by a short local spring relaxation, so bond lengths stay
chemical (HierGO left rotated bonds at 0.85 A or 1.5 A for later clean-up).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from .carbon_lattice import (A_CC, A_LAT, CarbonNetwork, clear_crowded_edges, prune_local,
                             spring_relax, tree_query_points)
from .schema import HOLE_SIZE_RANGES

HIERGO_A = 2.471  # lattice constant of HierGO's template (A)
ETHER_CO = 1.37
ETHER_OUT_OF_PLANE = 0.78


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def shells(net: CarbonNetwork, seeds, depth: int) -> dict:
    """Atom -> bond distance from ``seeds`` (framework bonds only), up to ``depth``."""
    dist = {s: 0 for s in seeds}
    frontier = list(seeds)
    for d in range(1, depth + 1):
        nxt = []
        for a in frontier:
            for b in net.c_nb(a):
                if b not in dist:
                    dist[b] = d
                    nxt.append(b)
        frontier = nxt
    return dist


def is_perfect(net: CarbonNetwork, atoms, blocked) -> bool:
    for a in atoms:
        if (not net.alive[a] or not net.lattice[a] or a in blocked or len(net.c_nb(a)) != 3
                or len(net.adj[a]) != 3 or net.role[a] != "C_sp2"):
            return False
    return True


def surface_constraint(net: CarbonNetwork):
    """Keep relaxed atoms on their surface: planar sheets fix z, tubes stay on the cylinder."""
    if net.kind == "sheet":
        return "fixz"
    if net.kind == "nanotube":
        def project(p, g):
            r = p.copy()
            r[:, 2] = 0.0
            r /= np.maximum(np.linalg.norm(r, axis=1), 1e-9)[:, None]
            return g - np.sum(g * r, axis=1)[:, None] * r
        return project
    return None


def relax_region(net: CarbonNetwork, core, free_depth=3, fixed_depth=2, extra_r0=None):
    """Spring-relax atoms within ``free_depth`` bonds of ``core`` (others fixed).

    Works on a local, unwrapped copy (minimum image around the first core atom)
    and writes the result back. Bonds: C-C 1.42 A, ether C-O 1.37 A; 1-3
    distances pulled softly toward 120-degree geometry.
    """
    core = [a for a in core if net.alive[a]]
    if not core:
        return
    dist = shells(net, core, free_depth + fixed_depth)
    for a in list(core):  # include non-framework atoms bonded to the core (ether O)
        for b in net.adj[a]:
            if not net.lattice[b]:
                dist.setdefault(b, 0)
    atoms = list(dist)
    index = {a: k for k, a in enumerate(atoms)}
    ref = net.pos[core[0]]
    local = np.array([ref + net.mic(net.pos[a] - ref) for a in atoms])
    free = np.array([dist[a] <= free_depth for a in atoms])
    bonds, r0 = [], []
    for a in atoms:
        for b in net.adj[a]:
            if b in index and a < b:
                bonds.append((index[a], index[b]))
                if net.el[a] == "C" and net.el[b] == "C":
                    r0.append(A_CC)
                else:
                    r0.append(ETHER_CO if (extra_r0 is None) else extra_r0)
    triples, d13, w13 = [], [], []
    for a in atoms:
        if net.lattice[a]:
            nbs = [b for b in net.c_nb(a) if b in index]
            target, weight = A_LAT, 0.3
        else:  # ether oxygen: C-O-C about 115 degrees, held firmly so the tear opens
            nbs = [b for b in net.adj[a] if b in index]
            target, weight = 2 * ETHER_CO * np.sin(np.radians(57.5)), 3.0
        for x in range(len(nbs)):
            for y in range(x + 1, len(nbs)):
                triples.append((index[nbs[x]], index[a], index[nbs[y]]))
                d13.append(target)
                w13.append(weight)
    if net.kind == "nanotube":
        diam = net.geometry.get("diameters_A", [])
        radius = {a: (diam[net.layer[a]] / 2 if net.layer[a] < len(diam)
                      else float(np.hypot(*net.pos[a, :2]))) for a in atoms}
    out = spring_relax(local, bonds, r0, triples, d13, free=free, w13=w13,
                       constraint=surface_constraint(net), maxiter=500)
    for a in atoms:
        if not free[index[a]]:
            continue
        p = out[index[a]]
        if net.kind == "nanotube" and net.lattice[a]:
            rxy = np.hypot(p[0], p[1])
            if rxy > 1e-6:
                p[:2] *= radius[a] / rxy
        net._pos[a] = net.wrap(p)


def bonds_ok(net: CarbonNetwork, atoms, lo=1.25, hi=1.68) -> bool:
    for a in atoms:
        if not net.alive[a]:
            continue
        for b in net.adj[a]:
            d = net.dist(a, b)
            if not lo <= d <= hi:
                return False
    return True


class _Snapshot:
    """Undo record for a tentative defect (positions + bonds of the touched atoms)."""

    def __init__(self, net, atoms):
        self.net = net
        self.atoms = set(atoms)
        for a in list(self.atoms):
            self.atoms |= net.adj[a]
        self.pos = {a: net.pos[a].copy() for a in self.atoms}
        self.adj = {a: set(net.adj[a]) for a in self.atoms}
        self.order = {k: v for k, v in net.order.items() if k[0] in self.atoms or k[1] in self.atoms}
        self.role = {a: net.role[a] for a in self.atoms}

    def restore(self):
        net = self.net
        for a in self.atoms:
            for b in list(net.adj[a]):
                net.remove_bond(a, b)
        for a, nbs in self.adj.items():
            for b in nbs:
                net.adj[a].add(b)
                net.adj[b].add(a)
        net.order.update(self.order)
        for a, p in self.pos.items():
            net._pos[a] = p
        for a, r in self.role.items():
            net.role[a] = r


# --------------------------------------------------------------------------
# vacancy holes
# --------------------------------------------------------------------------


def hole_size_range(cfg) -> tuple:
    if cfg.hole_size == "custom":
        return int(cfg.hole_min_atoms), int(cfg.hole_max_atoms)
    return HOLE_SIZE_RANGES[cfg.hole_size]


def hole_sizes(n_remove: int, lo: int, hi: int, rng) -> list:
    """HierGO create_defect_sizes: random sizes in [lo, hi] that add up to n_remove."""
    sizes, total = [], 0
    while total < n_remove:
        rem = n_remove - total
        if rem < hi:
            if rem >= lo:
                s = int(rng.integers(lo, rem + 1))
            elif sizes:
                for _ in range(rem):
                    sizes[int(rng.integers(len(sizes)))] += 1
                break
            else:
                s = rem
        else:
            s = int(rng.integers(lo, hi + 1))
        sizes.append(s)
        total += s
    return sizes


def _edge_buffer(net: CarbonNetwork, layer: int, depth: int = 2) -> set:
    """Edge atoms of a finite sheet/tube plus ``depth`` shells (kept free of new holes)."""
    edges = [a for a in net.lattice_atoms(layer) if len(net.c_nb(a)) < 3]
    return set(shells(net, edges, depth)) if edges else set()


def find_hole_library(path=None):
    """Directory of HierGO's HoleLibrary/*.xyz, if available (GPL-3.0; read, not bundled)."""
    candidates = []
    if path:
        candidates.append(Path(path))
    if os.environ.get("NANOBUILDER_HOLE_LIBRARY"):
        candidates.append(Path(os.environ["NANOBUILDER_HOLE_LIBRARY"]))
    here = Path(__file__).resolve().parent
    for parent in [here, *here.parents][:5]:
        for name in ("HierGO-master", "HierGO", "HierGO_Master", "HierGO-main"):
            candidates.append(parent / name / "HoleLibrary")
    for c in candidates:
        if c.is_dir() and any(c.glob("defect-*.xyz")):
            return c
    return None


_LIBRARY_CACHE = {}


def load_hole_library(path=None) -> dict:
    """{size: [xy arrays relative to the first atom]} from HierGO's HoleLibrary.

    Coordinates are rescaled from HierGO's template (a = 2.471 A) to this
    lattice (a = 2.4595 A, a_CC = 1.42 A); both have zigzag rows along x.
    """
    folder = find_hole_library(path)
    if folder is None:
        return {}
    key = str(folder)
    if key in _LIBRARY_CACHE:
        return _LIBRARY_CACHE[key]
    shapes = {}
    for f in sorted(folder.glob("defect-*.xyz")):
        lines = f.read_text().splitlines()
        try:
            n = int(lines[0].split()[0])
            xyz = np.array([[float(v) for v in ln.split()[1:4]] for ln in lines[2:2 + n]])
        except (ValueError, IndexError):
            continue
        if len(xyz) != n or n == 0:
            continue
        xy = xyz[:, :2] * (A_LAT / HIERGO_A)
        shapes.setdefault(n, []).append(xy - xy[0])
    _LIBRARY_CACHE[key] = shapes
    return shapes


def _grow_one(net: CarbonNetwork, layer: int, size: int, rng, blocked: set, pool: list):
    """One compact cluster of ``size`` atoms (HierGO random neighbour growth)."""
    pool[:] = [a for a in pool if a not in blocked and net.alive[a] and len(net.c_nb(a)) == 3]
    if not pool:
        return None
    seed = pool[int(rng.integers(len(pool)))]
    hole, members = [seed], {seed}
    frontier = [b for b in net.c_nb(seed) if b not in blocked]
    while len(hole) < size and frontier:
        a = frontier.pop(int(rng.integers(len(frontier))))
        if a in members or a in blocked:
            continue
        hole.append(a)
        members.add(a)
        # duplicates on purpose (HierGO): atoms next to several removed atoms are
        # picked more often, which keeps holes compact
        frontier += [b for b in net.c_nb(a) if b not in members and b not in blocked]
    return hole


def _stamp_one(net: CarbonNetwork, layer: int, size: int, rng, blocked: set, pool: list,
               library: dict, tree, tree_atoms):
    """One hole shape from the HierGO library, snapped onto the lattice."""
    best = min(library, key=lambda s: abs(s - size))
    shapes = library[best]
    z = None
    for _try in range(60):
        pool[:] = [a for a in pool if a not in blocked and net.alive[a] and len(net.c_nb(a)) == 3]
        if not pool:
            return None
        seed = pool[int(rng.integers(len(pool)))]
        shape = shapes[int(rng.integers(len(shapes)))]
        z = net.pos[seed, 2]
        for ref in range(min(4, len(shape))):
            if net.kind != "sheet":
                return None
            pts = np.column_stack([shape - shape[ref] + net.pos[seed, :2], np.full(len(shape), z)])
            d, idx = tree.query(tree_query_points(tree, pts), distance_upper_bound=0.4)
            if np.any(~np.isfinite(d)):
                continue
            hole = [tree_atoms[i] for i in idx]
            if len(set(hole)) != len(hole) or any((not net.alive[a]) or a in blocked or
                                                  net.layer[a] != layer for a in hole):
                continue
            return hole
    return None


def make_holes(net: CarbonNetwork, cfg, rng, blocked: set, clean_edges: bool,
               allow_ether_pairs: bool) -> dict:
    """Vacancy holes in every layer/wall/shell until the requested vacancy fraction is reached.

    Holes are made one at a time. With hydrogen-terminated edges each new hole
    is cleaned at once (dangling carbons trimmed, and atoms removed where edge
    hydrogens would collide), and those extra atoms count toward the target, so
    the final vacancy fraction matches the request even though very small holes
    (a monovacancy cannot carry three H atoms) grow while being cleaned.
    """
    lo, hi = hole_size_range(cfg)
    library = load_hole_library() if cfg.hole_shapes == "hiergo_library" else {}
    summary = {"requested_percent": float(cfg.vacancy_percent), "holes": 0, "removed_atoms": 0,
               "removed_by_cleaning": 0, "shapes": cfg.hole_shapes, "size_range_atoms": [lo, hi],
               "hole_sizes": []}
    if cfg.hole_shapes == "hiergo_library" and not library:
        net.warnings.append("HierGO HoleLibrary not found (looked for HierGO-master/HoleLibrary "
                            "next to the tool, or $NANOBUILDER_HOLE_LIBRARY); holes were grown "
                            "randomly instead")
    if library and net.kind != "sheet":
        library = {}
        net.warnings.append("HierGO library hole shapes are flat; holes in curved walls were "
                            "grown randomly")
    for layer in sorted(set(net.layer[a] for a in net.lattice_atoms())):
        n_layer = len(net.lattice_atoms(layer))
        target = int(n_layer * float(cfg.vacancy_percent) / 100.0)
        if target <= 0:
            continue
        layer_blocked = blocked | _edge_buffer(net, layer)
        pool = [a for a in net.lattice_atoms(layer) if a not in layer_blocked]
        tree_atoms = net.lattice_atoms(layer) if library else None
        tree = net.tree(tree_atoms) if library else None
        removed = 0
        failures = 0
        while removed < target and failures < 20:
            remaining = target - removed
            if remaining < lo:
                break
            size = int(rng.integers(lo, min(hi, remaining) + 1))
            if library:
                hole = _stamp_one(net, layer, size, rng, layer_blocked, pool, library, tree,
                                  tree_atoms)
            else:
                hole = _grow_one(net, layer, size, rng, layer_blocked, pool)
            if not hole:
                failures += 1
                continue
            saved_alive, saved_order = list(net.alive), dict(net.order)
            rim = set()
            for a in hole:
                rim |= set(net.c_nb(a))
            rim -= set(hole)
            for a in hole:
                net.kill(a)
            touched = set(rim)
            extra = prune_local(net, rim, cfg.keep_hanging_carbons, touched)
            if clean_edges:
                region = set(shells(net, [a for a in touched if net.alive[a]], 3))
                extra += clear_crowded_edges(net, cfg.keep_hanging_carbons, allow_ether_pairs,
                                             region=region)
            gone = list(hole) + list(extra)
            if removed + len(gone) - target > target - removed:
                # this hole (with its clean-up) would overshoot more than stopping short: undo it
                net.alive = saved_alive
                for (i, j), o in saved_order.items():
                    if (i, j) not in net.order:
                        net.add_bond(i, j, o)
                break
            removed += len(gone)
            summary["removed_by_cleaning"] += len(extra)
            summary["hole_sizes"].append(len(gone))
            near = [a for a in rim if net.alive[a]]
            layer_blocked |= set(shells(net, near, 2))
        blocked |= layer_blocked
        summary["holes"] = len(summary["hole_sizes"])
        summary["removed_atoms"] += removed
        if removed < target:
            net.warnings.append(f"Layer {layer}: only {removed} of {target} vacancies fit (holes are "
                                "kept apart and away from edges)")
    return summary


def make_pores(net: CarbonNetwork, count: int, diameter_nm: float, rng, blocked: set) -> dict:
    """Circular pores through every layer of a sheet, placed at random, non-overlapping."""
    radius = diameter_nm * 5.0
    lat = net.lattice_atoms()
    xy = net.pos[lat, :2]
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    margin = 0.0 if net.pbc[0] else radius + 4.0
    centres = []
    for _try in range(500 * max(count, 1)):
        if len(centres) == count:
            break
        c = lo + margin + rng.random(2) * np.maximum(hi - lo - 2 * margin, 1e-6)
        ok = True
        for d in centres:
            v = net.mic(np.array([*(c - d), 0.0]))[:2]
            if np.linalg.norm(v) < 2 * radius + 6.0:
                ok = False
                break
        if ok:
            centres.append(c)
    removed = 0
    for c in centres:
        v = net.mic(np.column_stack([net.pos[lat, :2] - c, np.zeros(len(lat))]))[:, :2]
        inside = np.linalg.norm(v, axis=1) <= radius
        for a, hit in zip(lat, inside):
            if hit and net.alive[a]:
                net.kill(a)
                removed += 1
    for c in centres:
        for a in net.lattice_atoms():
            v = net.mic(np.array([*(net.pos[a, :2] - c), 0.0]))[:2]
            if np.linalg.norm(v) <= radius + 4.5:
                blocked.add(a)
    if len(centres) < count:
        net.warnings.append(f"Only {len(centres)} of {count} pores of {diameter_nm:g} nm fit on "
                            "this sheet without overlapping")
    return {"pores": len(centres), "diameter_nm": float(diameter_nm), "removed_atoms": removed,
            "centres_A": [[round(float(x), 3) for x in c] for c in centres]}


# --------------------------------------------------------------------------
# topological defects
# --------------------------------------------------------------------------


def _rotate_about(points, centre, axis, angle):
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rot = np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k
    return (points - centre) @ rot.T + centre


def stone_wales(net: CarbonNetwork, count: int, rng, blocked: set) -> int:
    """Rotate ``count`` C-C bonds by 90 degrees (pentagon-heptagon pairs)."""
    if count <= 0:
        return 0
    bonds = [k for k in net.order if net.lattice[k[0]] and net.lattice[k[1]]
             and net.alive[k[0]] and net.alive[k[1]]]
    order = rng.permutation(len(bonds))
    placed = 0
    for idx in order:
        if placed >= count:
            break
        i, j = bonds[idx]
        env = shells(net, [i, j], 3)
        if not is_perfect(net, env, blocked) or j not in net.adj[i]:
            continue
        a_out = [x for x in net.c_nb(i) if x != j]
        b_out = [x for x in net.c_nb(j) if x != i]
        snap = _Snapshot(net, list(shells(net, [i, j], 5)))
        mid = net.pos[i] + net.delta(i, j) / 2
        axis = net.normal[i] if net.kind == "sheet" else 0.5 * (net.normal[i] + net.normal[j])
        pi_ = mid + net.mic(net.pos[i] - mid)
        pj_ = mid + net.mic(net.pos[j] - mid)
        new = _rotate_about(np.array([pi_, pj_]), mid, axis, np.pi / 2)
        for a in a_out:
            net.remove_bond(i, a)
        for b in b_out:
            net.remove_bond(j, b)
        net._pos[i], net._pos[j] = net.wrap(new[0]), net.wrap(new[1])
        outer = a_out + b_out
        best = None
        for pick in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
            to_i = [outer[p] for p in pick]
            to_j = [o for o in outer if o not in to_i]
            cost = sum(net.dist(i, o) for o in to_i) + sum(net.dist(j, o) for o in to_j)
            if best is None or cost < best[0]:
                best = (cost, to_i, to_j)
        for o in best[1]:
            net.add_bond(i, o, 1.5)
        for o in best[2]:
            net.add_bond(j, o, 1.5)
        relax_region(net, [i, j])
        if not bonds_ok(net, list(shells(net, [i, j], 2))) or any(len(net.c_nb(a)) != 3
                                                                    for a in [i, j] + outer):
            snap.restore()
            continue
        blocked |= set(shells(net, [i, j], 3))
        placed += 1
    return placed


def _hexagon_around(net: CarbonNetwork, a: int):
    """Atoms of one six-membered ring through atom ``a`` (or None)."""
    nbs = net.c_nb(a)
    for x in range(len(nbs)):
        for y in range(x + 1, len(nbs)):
            b, c = nbs[x], nbs[y]
            # path b - p - q - c of length 3 avoiding a
            for p in net.c_nb(b):
                if p == a:
                    continue
                for q in net.c_nb(p):
                    if q in (a, b):
                        continue
                    if c in net.c_nb(q):
                        return [a, b, p, q, c]
                    for r in net.c_nb(q):
                        if r in (a, b, p) or c not in net.c_nb(r):
                            continue
                        return [a, b, p, q, r, c]
    return None


def flower_defects(net: CarbonNetwork, count: int, rng, blocked: set) -> int:
    """Rotate ``count`` C24 units by 30 degrees (HierGO's 24-atom rotation)."""
    if count <= 0 or net.kind != "sheet":
        return 0
    atoms = net.lattice_atoms()
    order = rng.permutation(len(atoms))
    placed = 0
    for k in order:
        if placed >= count:
            break
        a = atoms[k]
        if not net.alive[a] or a in blocked:
            continue
        ring = _hexagon_around(net, a)
        if ring is None or len(ring) != 6:
            continue
        ref = net.pos[a]
        centre = ref + np.mean([net.mic(net.pos[r] - ref) for r in ring], axis=0)
        env = shells(net, ring, 6)
        if not is_perfect(net, env, blocked):
            continue
        cluster = [x for x in env if np.linalg.norm(net.mic(net.pos[x] - centre)) <= 3.85]
        if len(cluster) != 24:
            continue
        cset = set(cluster)
        external = [(u, v) for u in cluster for v in net.c_nb(u) if v not in cset]
        if len(external) != 12:
            continue
        snap = _Snapshot(net, list(shells(net, cluster, 5)))
        local = np.array([centre + net.mic(net.pos[x] - centre) for x in cluster])
        rotated = _rotate_about(local, centre, net.normal[a], np.pi / 6)
        for u, v in external:
            net.remove_bond(u, v)
        for x, p in zip(cluster, rotated):
            net._pos[x] = net.wrap(p)
        inner = sorted({u for u, _ in external})
        outer = sorted({v for _, v in external})
        if len(inner) != len(outer):
            snap.restore()
            continue
        cost = np.array([[net.dist(u, v) for v in outer] for u in inner])
        rows, cols = linear_sum_assignment(cost)
        for r, c in zip(rows, cols):
            net.add_bond(inner[r], outer[c], 1.5)
        relax_region(net, cluster, free_depth=3, fixed_depth=2)
        touched = list(shells(net, cluster, 2))
        if not bonds_ok(net, touched, 1.25, 1.70) or any(len(net.c_nb(x)) != 3 for x in touched):
            snap.restore()
            continue
        blocked |= set(shells(net, cluster, 3))
        placed += 1
    return placed


def ruptures(net: CarbonNetwork, count: int, length: int, rng, blocked: set) -> dict:
    """Tear ``count`` rows of ``length`` parallel C-C bonds, bridging each with an ether O."""
    if count <= 0 or net.kind != "sheet":
        return {"ruptures": 0, "ether_oxygens": 0}
    bonds = [k for k in net.order if net.lattice[k[0]] and net.lattice[k[1]]]
    placed, oxygens = 0, 0
    for idx in rng.permutation(len(bonds)):
        if placed >= count:
            break
        i0, j0 = bonds[idx]
        if not (net.alive[i0] and net.alive[j0]) or j0 not in net.adj[i0]:
            continue
        u = net.delta(i0, j0)
        u /= np.linalg.norm(u)
        nrm = net.normal[i0]
        p = np.cross(nrm, u)
        mid0 = net.pos[i0] + net.delta(i0, j0) / 2
        row = []
        for (a, b) in net.order:
            if not (net.lattice[a] and net.lattice[b]) or net.layer[a] != net.layer[i0]:
                continue
            d = net.delta(a, b)
            cos = abs(np.dot(d, u)) / np.linalg.norm(d)
            if cos < 0.99:
                continue
            mid = net.pos[a] + d / 2
            off = net.mic(mid - mid0)
            if abs(np.dot(off, u)) < 0.3 and abs(np.dot(off, nrm)) < 0.3:
                row.append((float(np.dot(off, p)), a, b))
        row.sort()
        start = [k for k, r in enumerate(row) if {r[1], r[2]} == {i0, j0}]
        if not start:
            continue
        s = start[0] - int(rng.integers(0, length))
        run = row[max(0, s):max(0, s) + length]
        if len(run) < length:
            continue
        steps = np.diff([r[0] for r in run])
        if np.any(np.abs(steps - A_LAT) > 0.2):  # must be a contiguous row
            continue
        involved = [x for r in run for x in (r[1], r[2])]
        if not is_perfect(net, shells(net, involved, 2), blocked):
            continue
        up = 1.0 if rng.random() < 0.5 else -1.0
        for _, a, b in run:
            d = net.delta(a, b)
            uab = d / np.linalg.norm(d)
            mid = net.pos[a] + d / 2
            net.remove_bond(a, b)
            net._pos[a] = net.wrap(net.pos[a] - 0.4 * uab)
            net._pos[b] = net.wrap(net.pos[b] + 0.4 * uab)
            o = net.add_atom("O", net.wrap(mid + up * ETHER_OUT_OF_PLANE * nrm), "O_ether",
                             net.layer[a], nrm, lattice=False)
            net.add_bond(a, o, 1.0)
            net.add_bond(b, o, 1.0)
            net.role[a] = net.role[b] = "C_ether"
            net.new_group("ether", [o], [a, b], face=int(up))
            for x in (a, b):
                for y in net.c_nb(x):
                    net.order[(min(x, y), max(x, y))] = 1.5
            up = -up
            oxygens += 1
        relax_region(net, involved, free_depth=3, fixed_depth=2)
        blocked |= set(shells(net, involved, 3))
        placed += 1
    if placed < count:
        net.warnings.append(f"Only {placed} of {count} ruptures of {length} bonds fit")
    return {"ruptures": placed, "ether_oxygens": oxygens, "length_bonds": int(length)}
