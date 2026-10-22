"""Carbon frameworks: graphene sheets, nanotubes, icosahedral fullerenes and onions.

Everything works on one explicit bond graph (:class:`CarbonNetwork`) built when
the lattice is created. Defects, oxidation, ligand attachment, typing and the
exported topology all use that graph, so nothing downstream has to guess bonds
from distances (the main source of "unbound atom" clean-up in the original
HierGO workflow).

Geometry conventions
- Graphene: zigzag direction along x, armchair along y, sheets normal to z.
  Orthogonal 4-atom cell a = sqrt(3) a_CC (x) by 3 a_CC (y); a_CC = 1.42 A.
- Nanotubes: axis along z, centred on x = y = 0 before boxing.
- Fullerenes: Goldberg (k, k) icosahedral shells C(60 k^2), centred on the origin.
"""

from __future__ import annotations

from collections import deque
from math import ceil, gcd, pi, sqrt

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull, KDTree, cKDTree

from .schema import BuilderError

A_CC = 1.42
A_LAT = A_CC * sqrt(3.0)  # 2.4595 A, zigzag repeat
B_LAT = 3.0 * A_CC  # 4.26 A, armchair repeat of the orthogonal cell
AREA_PER_C = A_LAT * B_LAT / 4.0  # 2.619 A^2 per carbon atom in one layer
C_H = 1.09
BOND_CUTOFF = 1.75
EDGE_RIBBON_NM = 0.335  # nominal "thickness" used to give edge regions an area


# --------------------------------------------------------------------------
# Network
# --------------------------------------------------------------------------


class CarbonNetwork:
    """Growable atom list with an explicit bond graph and per-atom normals.

    ``layer`` is the sheet layer / nanotube wall / onion shell index of each
    atom (atoms added later inherit their parent's). ``lattice`` marks
    framework carbons. ``faces_sites`` and ``faces_ox`` map a layer to the face
    signs (+1 outward/top, -1 bottom) that ligands and oxygen groups may use.
    """

    def __init__(self, positions, cell, pbc, layer, normals, kind):
        positions = np.asarray(positions, dtype=float).reshape(-1, 3)
        n = len(positions)
        cap = max(64, 2 * n)
        self._pos = np.zeros((cap, 3))
        self._pos[:n] = positions
        self._normal = np.zeros((cap, 3))
        self._normal[:n] = np.asarray(normals, dtype=float).reshape(-1, 3)
        self.n = n
        self.el = ["C"] * n
        self.role = ["C_sp2"] * n
        self.layer = [int(x) for x in layer]
        self.lattice = [True] * n
        self.alive = [True] * n
        self.group = [-1] * n
        self.adj = [set() for _ in range(n)]
        self.order = {}
        self.groups = []
        self.cell = np.array(cell, dtype=float)
        self.pbc = tuple(bool(p) for p in pbc)
        self.kind = kind
        self.faces_sites = {}
        self.faces_ox = {}
        self.geometry = {}
        self.info = {}
        self.warnings = []

    # -- storage ------------------------------------------------------------
    @property
    def pos(self) -> np.ndarray:
        return self._pos[:self.n]

    @property
    def normal(self) -> np.ndarray:
        return self._normal[:self.n]

    def add_atom(self, element, xyz, role, layer, normal=(0.0, 0.0, 1.0), lattice=False,
                 group=-1) -> int:
        if self.n == len(self._pos):
            self._pos = np.vstack([self._pos, np.zeros_like(self._pos)])
            self._normal = np.vstack([self._normal, np.zeros_like(self._normal)])
        i = self.n
        self._pos[i] = xyz
        self._normal[i] = normal
        self.n += 1
        self.el.append(element)
        self.role.append(role)
        self.layer.append(int(layer))
        self.lattice.append(bool(lattice))
        self.alive.append(True)
        self.group.append(int(group))
        self.adj.append(set())
        return i

    def add_bond(self, i, j, order=1.0):
        if i == j:
            raise ValueError("self bond")
        self.adj[i].add(j)
        self.adj[j].add(i)
        self.order[(min(i, j), max(i, j))] = float(order)

    def remove_bond(self, i, j):
        self.adj[i].discard(j)
        self.adj[j].discard(i)
        self.order.pop((min(i, j), max(i, j)), None)

    def bond_order(self, i, j) -> float:
        return self.order.get((min(i, j), max(i, j)), 1.0)

    def kill(self, i):
        for j in list(self.adj[i]):
            self.remove_bond(i, j)
        self.alive[i] = False

    def new_group(self, kind, atoms, carbons, face=0) -> int:
        gid = len(self.groups)
        self.groups.append({"kind": kind, "atoms": list(atoms), "carbons": list(carbons),
                            "face": int(face)})
        for a in atoms:
            self.group[a] = gid
        return gid

    # -- geometry -----------------------------------------------------------
    def mic(self, vec):
        """Minimum-image version of a vector (or array of vectors); orthogonal cells."""
        vec = np.array(vec, dtype=float)
        for k in range(3):
            if self.pbc[k]:
                length = self.cell[k, k]
                vec[..., k] -= length * np.round(vec[..., k] / length)
        return vec

    def delta(self, i, j):
        return self.mic(self._pos[j] - self._pos[i])

    def dist(self, i, j) -> float:
        return float(np.linalg.norm(self.delta(i, j)))

    def c_nb(self, i) -> list:
        """Framework-carbon neighbours of atom ``i``."""
        return [j for j in self.adj[i] if self.lattice[j] and self.alive[j]]

    def lattice_atoms(self, layer=None) -> list:
        return [i for i in range(self.n) if self.lattice[i] and self.alive[i]
                and (layer is None or self.layer[i] == layer)]

    def edge_atoms(self) -> list:
        """Framework carbons with two framework neighbours and nothing else attached."""
        return [i for i in self.lattice_atoms() if len(self.adj[i]) == 2]

    def edge_direction(self, i) -> np.ndarray:
        """Unit vector pointing away from the network at a two-coordinated atom."""
        u = np.zeros(3)
        for j in self.c_nb(i):
            d = self.delta(i, j)
            u += d / np.linalg.norm(d)
        n = np.linalg.norm(u)
        if n < 1e-6:  # straight chain (should have been pruned): any perpendicular
            d = self.delta(i, self.c_nb(i)[0])
            u = np.cross(d, self.normal[i])
            n = np.linalg.norm(u)
        out = -u / n
        # keep it in the local tangent plane
        nrm = self.normal[i]
        out = out - np.dot(out, nrm) * nrm
        return out / max(np.linalg.norm(out), 1e-9)

    def wrap(self, xyz):
        xyz = np.array(xyz, dtype=float)
        for k in range(3):
            if self.pbc[k]:
                xyz[..., k] %= self.cell[k, k]
        return xyz

    def tree(self, idx):
        """KD-tree of the given atoms (periodic along periodic axes)."""
        return make_tree(self._pos[idx], self.cell, self.pbc)

    # -- bookkeeping ----------------------------------------------------------
    def counts(self) -> dict:
        out = {}
        for i in range(self.n):
            if self.alive[i]:
                out[self.el[i]] = out.get(self.el[i], 0) + 1
        return out

    def compact(self):
        """Index map old -> new over living atoms, in the current order."""
        keep = [i for i in range(self.n) if self.alive[i]]
        return keep, {old: new for new, old in enumerate(keep)}


def make_tree(points, cell, pbc):
    """cKDTree that respects periodic axes of an orthogonal cell."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    if not any(pbc):
        return cKDTree(points)
    # Nonperiodic axes get a box far larger than the data (shifted to be positive),
    # so no pair ever wraps across them.
    margin = 1000.0
    box = np.zeros(3)
    shift = [0.0, 0.0, 0.0]
    pts = points.copy()
    for k in range(3):
        if pbc[k]:
            box[k] = cell[k, k]
            pts[:, k] = np.mod(pts[:, k], box[k])
            pts[pts[:, k] >= box[k], k] = 0.0  # -1e-17 % L == L in floating point
        else:
            shift[k] = (pts[:, k].min() if len(pts) else 0.0) - margin
            pts[:, k] -= shift[k]
            box[k] = (pts[:, k].max() if len(pts) else 0.0) + 1e6
    tree = KDTree(pts, boxsize=box)
    tree.shift = shift
    tree.box = box
    tree.pbc = tuple(bool(p) for p in pbc)
    return tree


def tree_query_points(tree, points):
    """Map query points into the frame of a tree built by :func:`make_tree`."""
    points = np.asarray(points, dtype=float).reshape(-1, 3).copy()
    if not hasattr(tree, "box"):
        return points
    for k in range(3):
        if tree.pbc[k]:
            points[:, k] = np.mod(points[:, k], tree.box[k])
            points[points[:, k] >= tree.box[k], k] = 0.0
        else:
            points[:, k] -= tree.shift[k]
    return points


def detect_bonds(positions, cell, pbc, cutoff=1.6):
    """(i, j) pairs closer than ``cutoff`` (minimum image along periodic axes)."""
    tree = make_tree(positions, cell, pbc)
    return sorted(tree.query_pairs(cutoff))


# --------------------------------------------------------------------------
# Graphene
# --------------------------------------------------------------------------


def graphene_layer(nx: int, ny: int) -> np.ndarray:
    basis = np.array([[0.25 * A_LAT, 0.5 * A_CC], [0.75 * A_LAT, 1.0 * A_CC],
                      [0.75 * A_LAT, 2.0 * A_CC], [0.25 * A_LAT, 2.5 * A_CC]])
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    offsets = np.stack([ii.ravel() * A_LAT, jj.ravel() * B_LAT], axis=1)
    xy = (offsets[:, None, :] + basis[None, :, :]).reshape(-1, 2)
    return np.column_stack([xy, np.zeros(len(xy))])


def build_graphene(cfg) -> CarbonNetwork:
    nx = max(2, int(round(cfg.length_nm * 10.0 / A_LAT)))
    ny = max(1, int(round(cfg.width_nm * 10.0 / B_LAT)))
    layers = int(cfg.layers)
    base = graphene_layer(nx, ny)
    lx, ly = nx * A_LAT, ny * B_LAT
    pos, layer_id = [], []
    for k in range(layers):
        p = base.copy()
        if cfg.stacking == "AB" and k % 2 == 1:
            p[:, 1] += A_CC
            if cfg.periodic:
                p[:, 1] %= ly
        p[:, 2] = k * float(cfg.layer_spacing_A)
        pos.append(p)
        layer_id += [k] * len(p)
    pos = np.concatenate(pos)
    z_span = (layers - 1) * float(cfg.layer_spacing_A)
    if cfg.periodic:
        cell = np.diag([lx, ly, z_span + 30.0])
        pbc = (True, True, False)
    else:
        cell = np.diag([lx + 30.0, ly + 30.0, z_span + 30.0])
        pbc = (False, False, False)
    normals = np.tile([0.0, 0.0, 1.0], (len(pos), 1))
    net = CarbonNetwork(pos, cell, pbc, layer_id, normals, "sheet")
    for i, j in detect_bonds(pos, cell, pbc, 1.6):
        if layer_id[i] == layer_id[j]:
            net.add_bond(i, j, 1.5)
    if layers == 1:
        net.faces_sites = {0: (1, -1)}
    else:
        net.faces_sites = {0: (-1,), layers - 1: (1,)}
    net.geometry = {"nx": nx, "ny": ny, "lx_A": lx, "ly_A": ly, "layers": layers,
                    "layer_spacing_A": float(cfg.layer_spacing_A), "stacking": cfg.stacking,
                    "periodic": bool(cfg.periodic)}
    if cfg.periodic and (abs(lx - cfg.length_nm * 10) > 0.05 * cfg.length_nm * 10 or
                         abs(ly - cfg.width_nm * 10) > 0.05 * cfg.width_nm * 10):
        net.warnings.append(f"Periodic sheet uses whole lattice repeats: {lx / 10:.2f} x "
                            f"{ly / 10:.2f} nm (requested {cfg.length_nm:g} x {cfg.width_nm:g} nm)")
    return net


# --------------------------------------------------------------------------
# Nanotubes
# --------------------------------------------------------------------------


def cnt_diameter(n: int, m: int) -> float:
    return A_LAT / pi * sqrt(n * n + n * m + m * m)


def cnt_period(n: int, m: int) -> float:
    dr = gcd(2 * n + m, 2 * m + n)
    return sqrt(3.0) * A_LAT * sqrt(n * n + n * m + m * m) / dr


def wall_chiralities(n: int, m: int, walls: int, spacing: float, periodic: bool) -> list:
    """Chiral indices of concentric walls, inner first, about ``spacing`` apart."""
    out = [(n, m)]
    d0 = cnt_diameter(n, m)
    theta0 = np.arctan2(sqrt(3) * m, 2 * n + m)
    for w in range(1, walls):
        target = d0 + 2 * spacing * w
        if m == n:
            k = int(round(target * pi / (A_LAT * sqrt(3))))
            out.append((k, k))
        elif m == 0:
            k = int(round(target * pi / A_LAT))
            out.append((k, 0))
        else:
            if periodic:
                raise BuilderError("A periodic multi-wall nanotube needs commensurate walls: use an "
                                   "armchair (n, n) or zigzag (n, 0) inner tube, or a finite length")
            best = None
            for nn in range(1, int(target * pi / A_LAT) + 3):
                for mm in range(0, nn + 1):
                    d = cnt_diameter(nn, mm)
                    th = np.arctan2(sqrt(3) * mm, 2 * nn + mm)
                    score = abs(d - target) + 2.0 * abs(th - theta0)
                    if best is None or score < best[0]:
                        best = (score, (nn, mm))
            out.append(best[1])
    return out


def build_nanotube(cfg) -> CarbonNetwork:
    from ase.build import nanotube

    n, m = int(cfg.cnt_n), int(cfg.cnt_m)
    walls = wall_chiralities(n, m, int(cfg.cnt_walls), float(cfg.wall_spacing_A), cfg.periodic)
    length = float(cfg.cnt_length_nm) * 10.0
    periods = [cnt_period(a, b) for a, b in walls]
    pos, layer_id = [], []
    if cfg.periodic:
        t = periods[0]
        cells = max(1, int(round(length / t)))
        lz = cells * t
        if any(abs(p - t) > 1e-6 for p in periods):
            raise BuilderError("Nanotube walls have different axial periods; use a finite length")
    else:
        lz = length
    for w, (a, b) in enumerate(walls):
        t = periods[w]
        cells = max(1, int(round(lz / t))) if cfg.periodic else int(ceil(lz / t)) + 2
        tube = nanotube(a, b, length=cells, bond=A_CC)
        p = tube.get_positions()
        if not cfg.periodic:
            zmid = 0.5 * cells * t
            keep = np.abs(p[:, 2] - zmid) <= lz / 2 + 1e-6
            p = p[keep]
            p[:, 2] -= zmid - lz / 2
        pos.append(p)
        layer_id += [w] * len(p)
    pos = np.concatenate(pos)
    rmax = cnt_diameter(*walls[-1]) / 2
    if cfg.periodic:
        cell = np.diag([2 * rmax + 30.0, 2 * rmax + 30.0, lz])
        pbc = (False, False, True)
    else:
        cell = np.diag([2 * rmax + 30.0, 2 * rmax + 30.0, lz + 30.0])
        pbc = (False, False, False)
    radial = pos.copy()
    radial[:, 2] = 0.0
    normals = radial / np.linalg.norm(radial, axis=1)[:, None]
    net = CarbonNetwork(pos, cell, pbc, layer_id, normals, "nanotube")
    for i, j in detect_bonds(pos, cell, pbc, 1.6):
        if layer_id[i] == layer_id[j]:
            net.add_bond(i, j, 1.5)
    outer = len(walls) - 1
    net.faces_sites = {outer: (1,)}
    net.geometry = {"walls": [list(w) for w in walls],
                    "diameters_A": [round(cnt_diameter(*w), 4) for w in walls],
                    "length_A": lz, "periodic": bool(cfg.periodic),
                    "axial_period_A": periods[0] if cfg.periodic else None}
    if cfg.periodic and abs(lz - length) > 0.05 * length:
        net.warnings.append(f"Periodic nanotube uses whole translational cells: length "
                            f"{lz / 10:.2f} nm (requested {cfg.cnt_length_nm:g} nm)")
    if not cfg.periodic and lz < 2 * B_LAT:
        net.warnings.append("Very short nanotube: most atoms are rim atoms")
    return net


# --------------------------------------------------------------------------
# Fullerenes and carbon onions (icosahedral Goldberg polyhedra, class II (k, k))
# --------------------------------------------------------------------------


def _icosahedron():
    phi = (1 + sqrt(5)) / 2
    v = np.array([(-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0), (0, -1, phi), (0, 1, phi),
                  (0, -1, -phi), (0, 1, -phi), (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1),
                  (-phi, 0, 1)], dtype=float)
    v /= np.linalg.norm(v, axis=1)[:, None]
    return v, ConvexHull(v).simplices


def goldberg_shell(k: int):
    """Unrelaxed C(60 k^2) on the unit sphere: (positions, bonds, triples, triple_ring_sizes).

    Geodesic (k, k) points on each icosahedron face -> convex hull (the
    geodesic triangulation) -> dual: one carbon per triangle, bonded across
    shared edges. Each hull vertex is a ring (valence 5 = pentagon).
    """
    verts, faces = _icosahedron()
    h = k
    e1, e2 = np.array([1.0, 0.0]), np.array([0.5, sqrt(3) / 2])
    p1 = h * e1 + k * e2
    p2 = -k * e1 + (h + k) * e2
    mat = np.column_stack([p1, p2])
    pts = []
    for a, b, c in faces:
        A, B, C = verts[a], verts[b], verts[c]
        for i in range(-k - 1, h + 2):
            for j in range(0, h + k + 2):
                p = i * e1 + j * e2
                l1, l2 = np.linalg.solve(mat, p)
                l0 = 1 - l1 - l2
                if min(l0, l1, l2) >= -1e-9:
                    v = l0 * A + l1 * B + l2 * C
                    pts.append(v / np.linalg.norm(v))
    pts = np.unique(np.round(np.array(pts), 9), axis=0)
    expected = 10 * 3 * k * k + 2
    if len(pts) != expected:
        raise BuilderError(f"Fullerene construction failed ({len(pts)} vs {expected} vertices)")
    hull = ConvexHull(pts)
    tri = hull.simplices
    centroids = pts[tri].mean(axis=1)
    centroids /= np.linalg.norm(centroids, axis=1)[:, None]
    bonds = sorted({(min(f, g), max(f, g)) for f in range(len(tri)) for g in hull.neighbors[f]})
    valence = np.zeros(len(pts), dtype=int)
    for t in tri:
        valence[t] += 1
    triples, rings = [], []
    vsets = [set(t) for t in tri]
    for f in range(len(tri)):
        nbs = list(hull.neighbors[f])
        for x in range(3):
            for y in range(x + 1, 3):
                common = vsets[f] & vsets[nbs[x]] & vsets[nbs[y]]
                ring = valence[next(iter(common))] if common else 6
                triples.append((nbs[x], f, nbs[y]))
                rings.append(int(ring))
    return centroids, bonds, triples, rings


def spring_relax(pos, bonds, r0, triples=(), d13=(), free=None, w13=0.3, constraint=None,
                 mic=None, maxiter=400):
    """Minimize sum (r - r0)^2 + w13 * sum (d13 - d13_0)^2 over the free atoms.

    ``constraint``: None, ``"fixz"`` (planar sheets: z frozen) or a callable
    taking (positions, gradient) and returning a projected gradient.
    """
    pos = np.array(pos, dtype=float)
    n = len(pos)
    free = np.ones(n, bool) if free is None else np.asarray(free, bool)
    fidx = np.flatnonzero(free)
    if not len(fidx) or not len(bonds):
        return pos
    b = np.asarray(bonds, dtype=int).reshape(-1, 2)
    r0 = np.broadcast_to(np.asarray(r0, dtype=float), (len(b),))
    t = np.asarray(triples, dtype=int).reshape(-1, 3)
    d0 = np.broadcast_to(np.asarray(d13, dtype=float), (len(t),)) if len(t) else np.zeros(0)
    w13 = np.broadcast_to(np.asarray(w13, dtype=float), (len(t),)) if len(t) else np.zeros(0)
    wrap = mic if mic is not None else (lambda v: v)

    def energy(x):
        p = pos.copy()
        p[fidx] = x.reshape(-1, 3)
        grad = np.zeros_like(p)
        v = wrap(p[b[:, 1]] - p[b[:, 0]])
        r = np.linalg.norm(v, axis=1)
        e = np.sum((r - r0) ** 2)
        g = (2 * (r - r0) / np.maximum(r, 1e-9))[:, None] * v
        np.add.at(grad, b[:, 1], g)
        np.add.at(grad, b[:, 0], -g)
        if len(t):
            v = wrap(p[t[:, 2]] - p[t[:, 0]])
            r = np.linalg.norm(v, axis=1)
            e += np.sum(w13 * (r - d0) ** 2)
            g = (w13 * 2 * (r - d0) / np.maximum(r, 1e-9))[:, None] * v
            np.add.at(grad, t[:, 2], g)
            np.add.at(grad, t[:, 0], -g)
        if constraint == "fixz":
            grad[:, 2] = 0.0
        elif callable(constraint):
            grad = constraint(p, grad)
        return e, grad[fidx].ravel()

    res = minimize(energy, pos[fidx].ravel(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "gtol": 1e-6})
    pos[fidx] = res.x.reshape(-1, 3)
    return pos


def fullerene_coordinates(k: int):
    """Relaxed C(60 k^2) centred at the origin: (positions, bonds, ring size per 1-3 triple)."""
    unit, bonds, triples, rings = goldberg_shell(k)
    lengths = np.linalg.norm(unit[[i for i, _ in bonds]] - unit[[j for _, j in bonds]], axis=1)
    pos = unit * (A_CC / lengths.mean())
    angle = {5: 108.0, 6: 120.0, 7: 128.57}
    d13 = [2 * A_CC * np.sin(np.radians(angle.get(r, 120.0)) / 2) for r in rings]
    pos = spring_relax(pos, bonds, A_CC, triples, d13, w13=0.5, maxiter=2000)
    pos -= pos.mean(axis=0)
    return pos, bonds


def build_carbon_sphere(cfg) -> CarbonNetwork:
    orders = ([int(cfg.fullerene_order)] if cfg.sphere_type == "fullerene"
              else list(range(1, int(cfg.onion_shells) + 1)))
    pos, layer_id, bond_list, radii = [], [], [], []
    offset = 0
    for shell, k in enumerate(orders):
        p, bonds = fullerene_coordinates(k)
        pos.append(p)
        layer_id += [shell] * len(p)
        bond_list += [(i + offset, j + offset) for i, j in bonds]
        radii.append(float(np.linalg.norm(p, axis=1).mean()))
        offset += len(p)
    pos = np.concatenate(pos)
    rmax = max(np.linalg.norm(pos, axis=1))
    cell = np.diag([2 * rmax + 30.0] * 3)
    normals = pos / np.linalg.norm(pos, axis=1)[:, None]
    net = CarbonNetwork(pos, cell, (False,) * 3, layer_id, normals, "carbon_sphere")
    for i, j in bond_list:
        net.add_bond(i, j, 1.5)
    net.faces_sites = {len(orders) - 1: (1,)}
    net.geometry = {"shells": [f"C{60 * k * k}" for k in orders],
                    "shell_radii_A": [round(r, 3) for r in radii]}
    if len(orders) > 1:
        gaps = np.diff(radii)
        net.geometry["shell_gaps_A"] = [round(g, 3) for g in gaps]
    return net


# --------------------------------------------------------------------------
# Graph clean-up (HierGO trim-hanging-Cs, generalized to 5/7 rings)
# --------------------------------------------------------------------------


def _in_small_ring(net: CarbonNetwork, v: int, max_ring: int = 8) -> bool:
    """True if atom ``v`` lies on a framework ring of at most ``max_ring`` atoms."""
    nbs = net.c_nb(v)
    for x in range(len(nbs)):
        target = set(nbs[x + 1:])
        if not target:
            continue
        start = nbs[x]
        seen = {v, start}
        frontier = [start]
        for _depth in range(max_ring - 2):
            nxt = []
            for a in frontier:
                for b in net.c_nb(a):
                    if b in target:
                        return True
                    if b not in seen:
                        seen.add(b)
                        nxt.append(b)
            frontier = nxt
            if not frontier:
                break
    return False


def kill_with_groups(net: CarbonNetwork, a: int) -> None:
    """Remove a framework atom together with any functional group hanging from it."""
    stack = [b for b in net.adj[a] if not net.lattice[b]]
    net.kill(a)
    seen = set()
    while stack:
        g = stack.pop()
        if g in seen or not net.alive[g]:
            continue
        seen.add(g)
        others = [x for x in net.adj[g] if net.lattice[x] and net.alive[x]]
        if others and net.el[g] == "O" and len(others) >= 1 and net.role[g] in ("O_epoxy", "O_ether"):
            # a bridging oxygen keeps its other carbon: turn it into nothing (remove O),
            # the remaining carbon reverts to sp2 below
            for x in others:
                net.role[x] = "C_sp2"
        stack += [x for x in net.adj[g] if not net.lattice[x]]
        net.kill(g)


def prune_local(net: CarbonNetwork, seeds, keep_hanging: bool = False, touched=None) -> list:
    """Local version of the dangling-carbon trim around ``seeds`` (recently changed atoms).

    ``touched`` (a set), if given, collects the atoms that were bonded to removed ones.
    """
    if keep_hanging:
        return []
    removed = []
    queue = deque(set(seeds))
    queued = set(queue)
    while queue:
        a = queue.popleft()
        queued.discard(a)
        if not (net.alive[a] and net.lattice[a]):
            continue
        deg = len(net.c_nb(a))
        if deg < 2 or (deg == 2 and not _in_small_ring(net, a)):
            nbs = net.c_nb(a)
            if touched is not None:
                touched.update(nbs)
            kill_with_groups(net, a)
            removed.append(a)
            for b in nbs:
                for c in [b] + net.c_nb(b):
                    if c not in queued:
                        queue.append(c)
                        queued.add(c)
    return removed


def prune(net: CarbonNetwork, keep_hanging: bool = False) -> int:
    """Remove dangling carbons and isolated fragments; returns the number removed.

    Dangling = fewer than two framework neighbours, or a two-coordinated atom
    that is not on any ring of up to 8 atoms (a chain). Repeats until stable,
    then keeps the largest connected piece of every layer.
    """
    removed = 0
    if not keep_hanging:
        changed = True
        while changed:
            changed = False
            for i in net.lattice_atoms():
                if len(net.c_nb(i)) < 2:
                    kill_with_groups(net, i)
                    removed += 1
                    changed = True
            for i in net.lattice_atoms():
                if len(net.c_nb(i)) == 2 and not _in_small_ring(net, i):
                    kill_with_groups(net, i)
                    removed += 1
                    changed = True
    for layer in sorted(set(net.layer[i] for i in net.lattice_atoms())):
        atoms = net.lattice_atoms(layer)
        seen, comps = set(), []
        for s in atoms:
            if s in seen:
                continue
            comp, queue = [], deque([s])
            seen.add(s)
            while queue:
                a = queue.popleft()
                comp.append(a)
                for b in net.c_nb(a):
                    if b not in seen:
                        seen.add(b)
                        queue.append(b)
            comps.append(comp)
        comps.sort(key=len, reverse=True)
        for comp in comps[1:]:
            for a in comp:
                if net.alive[a]:
                    kill_with_groups(net, a)
                    removed += 1
    return removed


def h_positions(net: CarbonNetwork, atoms) -> np.ndarray:
    return np.array([net.pos[i] + C_H * net.edge_direction(i) for i in atoms]).reshape(-1, 3)


def ether_pair_ok(net: CarbonNetwork, a: int, b: int) -> bool:
    """Two edge carbons that an in-plane ether O could bridge (a pyran-like 6-ring).

    2.2-2.75 A apart, not bonded, no common framework neighbour (that would
    make a strained 4-ring), and no third edge atom close to both (the O would
    then sit in a monovacancy, bonded to three carbons).
    """
    if b in net.adj[a]:
        return False
    d = net.dist(a, b)
    if not 2.2 <= d <= 2.75:
        return False
    if set(net.c_nb(a)) & set(net.c_nb(b)):
        return False
    return True


def crowded_edge_victims(net: CarbonNetwork, edges, allow_ether_pairs=False) -> set:
    """Edge carbons to remove so that every remaining edge H has room.

    H-H closer than 1.55 A, or an H closer than 1.6 A to a non-bonded carbon,
    marks a clash. With ``allow_ether_pairs`` a clash between the two H atoms
    of an ether-capable pair is tolerated (an ether oxygen may fill that site).
    """
    edges = [e for e in edges if net.alive[e]]
    if not edges:
        return set()
    hpos = h_positions(net, edges)
    victims = set()
    htree = make_tree(hpos, net.cell, net.pbc)
    for a, b in sorted(htree.query_pairs(1.55)):
        ea, eb = edges[a], edges[b]
        if ea in victims or eb in victims:
            continue
        if allow_ether_pairs and ether_pair_ok(net, ea, eb):
            # unless a third edge atom crowds the same spot (monovacancy)
            third = [e for e in edges if e not in (ea, eb) and net.dist(e, ea) < 2.75
                     and net.dist(e, eb) < 2.75]
            if not third:
                continue
        victims.add(eb)
    lat = net.lattice_atoms()
    ctree = net.tree(lat)
    hits = ctree.query_ball_point(tree_query_points(ctree, hpos), 1.6)
    for k, near in enumerate(hits):
        own = edges[k]
        if own in victims:
            continue
        for n in near:
            c = lat[n]
            if c != own and c not in net.adj[own]:
                victims.add(own)
                break
    return victims


def clear_crowded_edges(net: CarbonNetwork, keep_hanging=False, allow_ether_pairs=False,
                        region=None, max_rounds=400) -> list:
    """Remove edge carbons whose terminating H would collide; returns removed atoms.

    A real H-terminated monovacancy cannot exist: its three inward H atoms
    would sit 0.6 A apart. Such atoms are removed (the hole grows) until every
    edge carbon has room for its hydrogen. ``region`` limits the search to
    atoms near a recent change.
    """
    removed = []
    region = None if region is None else set(region)
    for _ in range(max_rounds):
        if region is None:
            edges = net.edge_atoms()
        else:
            region = {a for a in region if net.alive[a]}
            edges = [a for a in region if net.lattice[a] and len(net.adj[a]) == 2]
        victims = crowded_edge_victims(net, edges, allow_ether_pairs)
        if not victims:
            break
        seeds = set()
        for v in victims:
            seeds |= set(net.c_nb(v))
            kill_with_groups(net, v)
            removed.append(v)
        touched = set(seeds)
        pruned = prune_local(net, seeds, keep_hanging, touched)
        removed += pruned
        if region is not None:
            # widen the search around everything that changed, including atoms the
            # trim removed further away
            frontier = {a for a in touched if net.alive[a]}
            for _ in range(3):
                frontier |= {c for b in frontier for c in net.c_nb(b)}
            region |= frontier
    return removed
