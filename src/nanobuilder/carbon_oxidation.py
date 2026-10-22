"""Oxygen functionalization of carbon frameworks (graphene oxide and oxidized CNTs/fullerenes).

Re-implements HierGO's decorate-tile.py on the explicit bond graph:

- **fresh**: basal hydroxyls and epoxides in a hydroxyl:epoxy ratio (HierGO
  default 1:1; the AMBER GO dataset of Fedai et al. 2026 used 2:1). A share of
  the hydroxyls (HierGO: 70 %) is placed as anti 1,2-diol pairs on adjacent
  carbons, one above and one below the sheet.
- **aged** (600 / 900 / 1200 / 2500 K): HierGO's annealing profiles, adding
  edge hydroxyls, carbonyls, carboxylic acids and in-plane (pyran-like)
  ethers bridging two edge carbons of a vacancy.
- **custom**: any mix of the six group types.

Placement rules: a carbon carries at most one group; two adjacent carbons may
not both carry a group on the same face (except the two carbons of one
epoxide), which reproduces the alternating up/down pattern needed for high
oxidation; new atoms must keep minimum distances from every atom placed so
far; faces are balanced as in HierGO. Groups that cannot be placed fall back
to basal hydroxyls (the strategy of the AMBER GO workflow), and the achieved
composition is always reported. Geometry uses chemical bond lengths (C-OH
1.43 A, O-H 0.97 A, epoxide C-O 1.46 A) with sp3 carbons pyramidalized by
0.25 A, instead of HierGO's deliberately long 1.7 A C-O.
"""

from __future__ import annotations

from collections import defaultdict, deque

import numpy as np

from .carbon_lattice import CarbonNetwork, ether_pair_ok, make_tree, tree_query_points
from .schema import BuilderError, parse_ratio

PYRAMIDALIZE_A = 0.25
C_O_HYDROXYL = 1.43
O_H = 0.97
C_O_PHENOL = 1.36
C_O_CARBONYL = 1.23
C_C_CARBOXYL = 1.50
C_O_ACID_DOUBLE = 1.21
C_O_ACID_SINGLE = 1.35
C_O_ETHER = 1.37
EPOXY_HEIGHT = 1.276  # O above the C-C midpoint for C-O = 1.46 A with C-C = 1.42 A
COH_ANGLE = np.radians(108.5)

# Fractions of the added oxygen atoms per group (a carboxylic acid counts as 2 O).
# HierGO's aged profiles (decorate-tile.py), normalized to sum to 1. "paired" is
# HierGO's "dihydroxyl" (anti 1,2-diols), "hydroxyl" its single surface hydroxyl.
AGED_PROFILES = {
    "aged_600K": {"paired": 0.19, "hydroxyl": 0.16, "edge_hydroxyl": 0.04, "epoxy": 0.43,
                  "ether": 0.05, "carbonyl": 0.07, "carboxyl": 0.12},
    "aged_900K": {"paired": 0.15, "hydroxyl": 0.16, "edge_hydroxyl": 0.05, "epoxy": 0.34,
                  "ether": 0.08, "carbonyl": 0.11, "carboxyl": 0.22},
    "aged_1200K": {"paired": 0.12, "hydroxyl": 0.14, "edge_hydroxyl": 0.06, "epoxy": 0.22,
                   "ether": 0.18, "carbonyl": 0.14, "carboxyl": 0.28},
    "aged_2500K": {"paired": 0.02, "hydroxyl": 0.02, "edge_hydroxyl": 0.04, "epoxy": 0.15,
                   "ether": 0.47, "carbonyl": 0.15, "carboxyl": 0.30},
}

# Minimum non-bonded distances (A) between a new atom and existing atoms.
# H...O down to 1.6 A is allowed: dense hydroxyl layers form O-H...O hydrogen bonds (1.6-2.0 A).
MIN_DIST = {("O", "O"): 2.40, ("H", "O"): 1.60, ("H", "H"): 1.50, ("C", "O"): 2.20,
            ("C", "H"): 1.80, ("C", "C"): 2.30}


def min_dist(a: str, b: str) -> float:
    return MIN_DIST.get(tuple(sorted((a, b))), 2.0)


def group_fractions(cfg) -> tuple:
    """(fractions of O per group incl. 'paired', paired share of all hydroxyl O)."""
    if cfg.oxidation_profile == "fresh":
        h, e = parse_ratio(cfg.hydroxyl_epoxy_ratio)
        hyd = h / (h + e)
        paired = float(cfg.paired_hydroxyl_fraction)
        f = {"hydroxyl": hyd * (1 - paired), "paired": hyd * paired, "epoxy": e / (h + e)}
    elif cfg.oxidation_profile == "custom":
        raw = {k: max(0.0, float(v)) for k, v in cfg.custom_fractions.items()}
        total = sum(raw.values())
        f = {k: v / total for k, v in raw.items()}
        paired = float(cfg.paired_hydroxyl_fraction)
        hyd = f.pop("hydroxyl", 0.0)
        f["hydroxyl"], f["paired"] = hyd * (1 - paired), hyd * paired
    else:
        raw = AGED_PROFILES[cfg.oxidation_profile]
        total = sum(raw.values())
        f = {k: v / total for k, v in raw.items()}
    for k in ("hydroxyl", "paired", "epoxy", "edge_hydroxyl", "carbonyl", "carboxyl", "ether"):
        f.setdefault(k, 0.0)
    hyd = f["hydroxyl"] + f["paired"]
    return f, (f["paired"] / hyd if hyd > 0 else 0.0)


def needs_ethers(cfg) -> bool:
    if float(cfg.oxygen_percent) <= 0 and cfg.oxidation_pattern != "tiles":
        return False
    return group_fractions(cfg)[0]["ether"] > 0


# --------------------------------------------------------------------------
# spatial bookkeeping
# --------------------------------------------------------------------------


class SpatialHash:
    """Grid of placed atoms (periodic along periodic axes of an orthogonal cell)."""

    def __init__(self, cell, pbc, size=2.5):
        self.size = size
        self.cell = np.asarray(cell, float)
        self.pbc = tuple(bool(p) for p in pbc)
        self.n = [max(1, int(self.cell[k, k] // size)) if self.pbc[k] else 0 for k in range(3)]
        self.grid = defaultdict(set)
        self.items = {}
        self.by_owner = {}
        self._next = 0

    def _key(self, p):
        key = []
        for k in range(3):
            if self.pbc[k]:
                length = self.cell[k, k]
                key.append(int(np.floor((p[k] % length) / length * self.n[k])) % self.n[k])
            else:
                key.append(int(np.floor(p[k] / self.size)))
        return tuple(key)

    def add(self, p, element, owner=None) -> int:
        i = self._next
        self._next += 1
        self.items[i] = (np.array(p, float), element, owner)
        self.grid[self._key(p)].add(i)
        if owner is not None:
            self.by_owner[owner] = i
        return i

    def remove(self, i):
        if i in self.items:
            p, _, owner = self.items.pop(i)
            self.grid[self._key(p)].discard(i)
            if owner is not None and self.by_owner.get(owner) == i:
                del self.by_owner[owner]

    def remove_owner(self, owner):
        if owner in self.by_owner:
            self.remove(self.by_owner[owner])

    def near(self, p, mic):
        base = self._key(p)
        keys = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = []
                    for k, d in enumerate((dx, dy, dz)):
                        v = base[k] + d
                        key.append(v % self.n[k] if self.pbc[k] else v)
                    keys.add(tuple(key))
        for key in keys:
            for i in self.grid.get(key, ()):
                q, el, owner = self.items[i]
                yield i, float(np.linalg.norm(mic(q - p))), el, owner


class Decorator:
    """Places oxygen groups on one CarbonNetwork."""

    def __init__(self, net: CarbonNetwork, cfg, rng):
        self.net = net
        self.cfg = cfg
        self.rng = rng
        self.delta = PYRAMIDALIZE_A if cfg.pyramidalize else 0.0
        self.hash = SpatialHash(net.cell, net.pbc)
        self.face_used = defaultdict(int)  # carbon -> face sign of its basal group
        self.faces_count = {1: 0, -1: 0}
        self.pending_h = {}  # edge carbon -> hash id of its future hydrogen
        self.counts = defaultdict(int)
        self.warnings = []
        lat = net.lattice_atoms()
        self.lat = lat
        self.lat_tree = net.tree(lat)
        for i in range(net.n):  # atoms that already exist (e.g. rupture ethers)
            if net.alive[i] and not net.lattice[i]:
                self.hash.add(net.pos[i], net.el[i], owner=i)
                if net.el[i] == "O":
                    self.counts["existing_O"] += 1
        self.ox_faces = self._oxidizable_faces()
        for e in net.edge_atoms():
            self._add_pending_h(e)
        self.colour = self._two_colouring()
        self.sublattice_bias = 0.0

    def _two_colouring(self) -> dict:
        """+1/-1 by sublattice (BFS 2-colouring; odd rings near defects just conflict)."""
        net = self.net
        colour = {}
        for s in self.lat:
            if s in colour:
                continue
            colour[s] = 1
            queue = deque([s])
            while queue:
                a = queue.popleft()
                for b in net.c_nb(a):
                    if b not in colour:
                        colour[b] = -colour[a]
                        queue.append(b)
        return colour

    # -- faces ---------------------------------------------------------------
    def _oxidizable_faces(self) -> dict:
        net, cfg = self.net, self.cfg
        layers = sorted(set(net.layer[a] for a in self.lat))
        if net.kind != "sheet":
            outer = max(layers)
            return {outer: (1,)}
        want = {"both": (1, -1), "top": (1,), "bottom": (-1,)}[cfg.oxidation_faces]
        if len(layers) == 1:
            return {layers[0]: want}
        if cfg.oxidize_all_layers:
            return {layer: want for layer in layers}
        out = {}
        if 1 in want:
            out[max(layers)] = (1,)
        if -1 in want:
            out[min(layers)] = tuple(sorted(set(out.get(min(layers), ())) | {-1}))
        return out

    # -- geometry helpers ---------------------------------------------------------
    def _basis(self, n):
        t = np.cross(n, [1.0, 0.0, 0.0])
        if np.linalg.norm(t) < 0.1:
            t = np.cross(n, [0.0, 1.0, 0.0])
        t /= np.linalg.norm(t)
        return t, np.cross(n, t)

    def _clear(self, p, element, exclude_lattice=(), exclude_items=()) -> bool:
        net = self.net
        for i, d, el, owner in self.hash.near(p, net.mic):
            if i in exclude_items:
                continue
            if d < min_dist(element, el):
                return False
        reach = 2.3 if element != "H" else 1.8
        q = tree_query_points(self.lat_tree, [p])
        for k in self.lat_tree.query_ball_point(q[0], reach):
            c = self.lat[k]
            if c in exclude_lattice or not net.alive[c]:
                continue
            d = float(np.linalg.norm(net.mic(net.pos[c] - p)))
            if d < min_dist(element, "C"):
                return False
        return True

    def _shell13(self, carbons) -> set:
        out = set(carbons)
        for c in carbons:
            out |= set(self.net.c_nb(c))
        return out

    def _add_pending_h(self, e):
        p = self.net.pos[e] + 1.09 * self.net.edge_direction(e)
        self.pending_h[e] = self.hash.add(p, "H", owner=("pendingH", e))

    def _drop_pending_h(self, e):
        if e in self.pending_h:
            self.hash.remove(self.pending_h.pop(e))

    # -- validity -----------------------------------------------------------------
    def basal_ok(self, c, f, partner=None) -> bool:
        net = self.net
        if not (net.alive[c] and net.lattice[c]) or net.role[c] != "C_sp2":
            return False
        if len(net.adj[c]) != 3 or len(net.c_nb(c)) != 3 or self.face_used.get(c):
            return False
        if f not in self.ox_faces.get(net.layer[c], ()):
            return False
        for nb in net.c_nb(c):
            if nb != partner and self.face_used.get(nb) == f:
                return False
        return True

    def edge_ok(self, e) -> bool:
        net = self.net
        return (net.alive[e] and net.lattice[e] and net.role[e] == "C_sp2" and len(net.adj[e]) == 2
                and net.layer[e] in self.ox_faces)

    def _pref_faces(self, faces, c=None):
        """Face order to try: balanced at random, or by sublattice for dense oxidation.

        Hydroxyls on the two sublattices pointing to opposite faces never block each
        other (the chair pattern of graphane), so at high oxygen content the face is
        taken from the sublattice with probability ``sublattice_bias``.
        """
        faces = list(faces)
        if len(faces) == 2 and c is not None and self.sublattice_bias > 0 \
                and self.rng.random() < self.sublattice_bias:
            f = self.colour.get(c, 1)
            return [f, -f]
        if len(faces) == 2:
            if self.faces_count[1] < self.faces_count[-1]:
                return [1, -1]
            if self.faces_count[-1] < self.faces_count[1]:
                return [-1, 1]
            return [1, -1] if self.rng.random() < 0.5 else [-1, 1]
        return faces

    # -- placement of one group --------------------------------------------------------
    def _pyramidalize(self, c, nf):
        if self.delta:
            self.net._pos[c] = self.net.wrap(self.net.pos[c] + self.delta * nf)

    def _commit(self, kind, new_atoms, carbons, face, carbon_role):
        """new_atoms: list of (element, position, role, bonded_to_index_or_'C'+k, order)."""
        net = self.net
        idx = []
        for el, pos, role, partners in new_atoms:
            i = net.add_atom(el, net.wrap(pos), role, net.layer[carbons[0]],
                             net.normal[carbons[0]], lattice=False)
            idx.append(i)
            for ref, order in partners:
                j = carbons[ref[1]] if isinstance(ref, tuple) else idx[ref]
                net.add_bond(i, j, order)
            self.hash.add(net.pos[i], el, owner=i)
        for c in carbons:
            if carbon_role:
                net.role[c] = carbon_role
        net.new_group(kind, idx, carbons, face)
        return idx

    def place_hydroxyl(self, c, f) -> bool:
        net = self.net
        nf = f * net.normal[c]
        base = net.pos[c] + self.delta * nf
        o = base + C_O_HYDROXYL * nf
        excl = self._shell13([c])
        if not self._clear(o, "O", excl):
            return False
        t1, t2 = self._basis(nf)
        best = None
        for k in range(12):
            ang = np.pi / 6 * k + self.rng.uniform(0, np.pi / 6)
            u = np.cos(ang) * t1 + np.sin(ang) * t2
            h = o + O_H * (np.cos(np.pi - COH_ANGLE) * nf + np.sin(np.pi - COH_ANGLE) * u)
            if not self._clear(h, "H", excl):
                continue
            score = min((d for _, d, _, _ in self.hash.near(h, net.mic)), default=9.0)
            if best is None or score > best[0]:
                best = (score, h)
        if best is None:
            return False
        self._pyramidalize(c, nf)
        self._commit("hydroxyl", [("O", o, "O_hydroxyl", [(("C", 0), 1.0)]),
                                  ("H", best[1], "H_hydroxyl", [(0, 1.0)])],
                     [c], f, "C_hydroxyl")
        self.face_used[c] = f
        self.faces_count[f] += 1
        self.counts["hydroxyl"] += 1
        return True

    def place_epoxy(self, c1, c2, f) -> bool:
        net = self.net
        nrm = net.normal[c1] + net.normal[c2]
        nf = f * nrm / np.linalg.norm(nrm)
        mid = net.pos[c1] + net.delta(c1, c2) / 2 + self.delta * nf
        o = mid + EPOXY_HEIGHT * nf
        if not self._clear(o, "O", self._shell13([c1, c2])):
            return False
        self._pyramidalize(c1, nf)
        self._pyramidalize(c2, nf)
        self._commit("epoxy", [("O", o, "O_epoxy", [(("C", 0), 1.0), (("C", 1), 1.0)])],
                     [c1, c2], f, "C_epoxy")
        self.face_used[c1] = self.face_used[c2] = f
        self.faces_count[f] += 1
        self.counts["epoxy"] += 1
        return True

    def place_edge_hydroxyl(self, e) -> bool:
        net = self.net
        d = net.edge_direction(e)
        n = net.normal[e]
        o = net.pos[e] + C_O_PHENOL * d
        self._drop_pending_h(e)
        excl = self._shell13([e])
        if not self._clear(o, "O", excl):
            self._add_pending_h(e)
            return False
        perp = np.cross(n, d)
        best = None
        for u in (perp, -perp, n, -n):
            h = o + O_H * (np.cos(np.pi - COH_ANGLE) * d + np.sin(np.pi - COH_ANGLE) * u)
            if self._clear(h, "H", excl):
                score = min((dd for _, dd, _, _ in self.hash.near(h, net.mic)), default=9.0)
                if best is None or score > best[0]:
                    best = (score, h)
        if best is None:
            self._add_pending_h(e)
            return False
        self._commit("edge_hydroxyl", [("O", o, "O_phenol", [(("C", 0), 1.0)]),
                                       ("H", best[1], "H_phenol", [(0, 1.0)])], [e], 0, "C_phenol")
        self.counts["edge_hydroxyl"] += 1
        return True

    def place_carbonyl(self, e) -> bool:
        net = self.net
        o = net.pos[e] + C_O_CARBONYL * net.edge_direction(e)
        self._drop_pending_h(e)
        if not self._clear(o, "O", self._shell13([e])):
            self._add_pending_h(e)
            return False
        self._commit("carbonyl", [("O", o, "O_carbonyl", [(("C", 0), 2.0)])], [e], 0, "C_carbonyl")
        self.counts["carbonyl"] += 1
        return True

    def place_carboxyl(self, e) -> bool:
        net = self.net
        d = net.edge_direction(e)
        n = net.normal[e]
        cc = net.pos[e] + C_C_CARBOXYL * d
        self._drop_pending_h(e)
        excl = self._shell13([e])
        if not self._clear(cc, "C", excl):
            self._add_pending_h(e)
            return False
        options = []
        for s in ((1, -1) if self.rng.random() < 0.5 else (-1, 1)):
            sn = s * n
            o1 = cc + C_O_ACID_DOUBLE * (0.5 * d + np.sqrt(0.75) * sn)
            o2 = cc + C_O_ACID_SINGLE * (0.5 * d - np.sqrt(0.75) * sn)
            v = (cc - o2) / np.linalg.norm(cc - o2)
            w = np.cross(np.cross(d, sn), v)  # in the COOH plane, perpendicular to v
            w /= np.linalg.norm(w)
            hs = [o2 + O_H * (np.cos(np.radians(107)) * v + s2 * np.sin(np.radians(107)) * w)
                  for s2 in (1, -1)]
            h = min(hs, key=lambda x: np.linalg.norm(x - o1))  # syn conformer
            if (self._clear(o1, "O", excl) and self._clear(o2, "O", excl)
                    and self._clear(h, "H", excl)):
                options.append((o1, o2, h))
        if not options:
            self._add_pending_h(e)
            return False
        o1, o2, h = options[0]
        self._commit("carboxyl", [("C", cc, "C_carboxyl", [(("C", 0), 1.0)]),
                                  ("O", o1, "O_carboxyl", [(0, 2.0)]),
                                  ("O", o2, "O_carboxyl_H", [(0, 1.0)]),
                                  ("H", h, "H_carboxyl", [(2, 1.0)])], [e], 0, "C_ring_cooh")
        self.counts["carboxyl"] += 1
        return True

    def place_ether(self, e1, e2) -> bool:
        net = self.net
        v = net.delta(e1, e2)
        dlen = np.linalg.norm(v)
        half = dlen / 2
        if half >= C_O_ETHER - 0.05:
            return False
        s = np.sqrt(C_O_ETHER ** 2 - half ** 2)
        n = net.normal[e1]
        p = np.cross(n, v / dlen)
        mid = net.pos[e1] + v / 2
        self._drop_pending_h(e1)
        self._drop_pending_h(e2)
        excl = self._shell13([e1, e2])
        best = None
        for o in (mid + s * p, mid - s * p):
            if not self._clear(o, "O", excl):
                continue
            q = tree_query_points(self.lat_tree, [o])[0]
            dd, _ = self.lat_tree.query(q, k=4)
            score = float(np.sort(dd)[2]) if np.ndim(dd) else float(dd)
            if best is None or score > best[0]:
                best = (score, o)
        if best is None:
            self._add_pending_h(e1)
            self._add_pending_h(e2)
            return False
        self._commit("ether", [("O", best[1], "O_ether", [(("C", 0), 1.0), (("C", 1), 1.0)])],
                     [e1, e2], 0, "C_ether")
        self.counts["ether"] += 1
        return True


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def _layout(net: CarbonNetwork, atoms) -> np.ndarray:
    """2D coordinates of atoms for spatial patterns: (x, y) on sheets, (arc, z) on tubes."""
    p = net.pos[atoms]
    if net.kind == "nanotube":
        r = np.hypot(p[:, 0], p[:, 1])
        return np.column_stack([np.arctan2(p[:, 1], p[:, 0]) * r.mean(), p[:, 2]])
    return p[:, :2]


def _distance_from_edges(net: CarbonNetwork, atoms) -> dict:
    seeds = [a for a in net.lattice_atoms() if len(net.c_nb(a)) < 3]
    dist = {s: 0 for s in seeds}
    queue = deque(seeds)
    while queue:
        a = queue.popleft()
        for b in net.c_nb(a):
            if b not in dist:
                dist[b] = dist[a] + 1
                queue.append(b)
    return dist


def decorate(net: CarbonNetwork, cfg, rng) -> dict:
    """Add oxygen groups to ``net`` following ``cfg``; returns a placement report."""
    dec = Decorator(net, cfg, rng)
    fractions, _ = group_fractions(cfg)
    lat = [a for a in dec.lat if net.layer[a] in dec.ox_faces]
    if not lat:
        return {"requested_percent": float(cfg.oxygen_percent), "placed": {}, "warnings": []}

    # Regions: one region for uniform/edges-first/domains, one per tile for tiles.
    regions = []
    if cfg.oxidation_pattern == "tiles":
        if net.kind == "carbon_sphere":
            raise BuilderError("The tiles pattern is for sheets and nanotubes; use domains on spheres")
        xy = _layout(net, lat)
        lo, hi = xy.min(axis=0), xy.max(axis=0) + 1e-6
        gx, gy = (int(i) for i in cfg.tile_grid)
        cell_id = (np.floor((xy[:, 0] - lo[0]) / (hi[0] - lo[0]) * gx).astype(int) * gy
                   + np.floor((xy[:, 1] - lo[1]) / (hi[1] - lo[1]) * gy).astype(int))
        tiles_o = []
        for t in range(gx * gy):
            members = set(np.asarray(lat)[cell_id == t].tolist())
            pct = float(rng.uniform(cfg.tile_oxygen_min, cfg.tile_oxygen_max))
            tiles_o.append(round(pct, 2))
            regions.append((members, pct))
        dec.tiles = {"grid": [gx, gy], "oxygen_percent": tiles_o}
    else:
        regions.append((set(lat), float(cfg.oxygen_percent)))

    order_key = None
    if cfg.oxidation_pattern == "edges_first":
        dist = _distance_from_edges(net, lat)
        if not dist:
            dec.warnings.append("Edges-first oxidation: this structure has no edges or holes; "
                                "groups were placed uniformly")
        order_key = lambda a: dist.get(a, 10 ** 6)  # noqa: E731
    domain_members = None
    if cfg.oxidation_pattern == "domains":
        xy = _layout(net, lat)
        radius = cfg.domain_diameter_nm * 5.0
        centres = []
        for k in rng.permutation(len(lat)):
            c = xy[k]
            if all(np.linalg.norm(net.mic(np.array([*(c - d), 0.0]))[:2]) > 2 * radius
                   for d in centres) or not centres:
                centres.append(c)
            if len(centres) >= int(cfg.domain_count):
                break
        inside = set()
        for k, a in enumerate(lat):
            for c in centres:
                if np.linalg.norm(net.mic(np.array([*(xy[k] - c), 0.0]))[:2]) <= radius:
                    inside.add(a)
                    break
        domain_members = inside
        dec.domains = {"count": len(centres), "diameter_nm": float(cfg.domain_diameter_nm),
                       "centres_A": [[round(float(x), 3) for x in c] for c in centres]}

    requested_o_total = 0
    fallback_o = 0
    # Oxygen content is 100 N_O / N_C over the whole substrate (the measured quantity), even
    # when only some layers/walls/shells can be oxidized.
    scale = len(dec.lat) / max(1, len(lat))
    densest = max(pct for _, pct in regions)
    dec.sublattice_bias = float(np.clip((densest - 25.0) / 20.0, 0.0, 1.0))
    for members, pct in regions:
        n_c = len(members) * scale
        p = pct / 100.0
        fc = fractions["carboxyl"]
        n_o = p * n_c / max(1e-9, 1 - p * fc / 2)  # carboxyl carbons count toward N_C
        if len(regions) == 1:
            n_o -= dec.counts["existing_O"]
        n_o = max(0, int(round(n_o)))
        requested_o_total += n_o
        q = {k: fractions[k] * n_o for k in fractions}
        quota = {
            "ether": int(round(q["ether"])),
            "carboxyl": int(round(q["carboxyl"] / 2)),
            "carbonyl": int(round(q["carbonyl"])),
            "edge_hydroxyl": int(round(q["edge_hydroxyl"])),
            "paired": int(round(q["paired"] / 2)),
            "epoxy": int(round(q["epoxy"])),
        }
        used_o = (quota["ether"] + 2 * quota["carboxyl"] + quota["carbonyl"]
                  + quota["edge_hydroxyl"] + 2 * quota["paired"] + quota["epoxy"])
        quota["hydroxyl"] = max(0, n_o - used_o)
        fallback_o += _place_region(dec, members, quota, order_key, domain_members)

    placed = {k: dec.counts[k] for k in ("hydroxyl", "epoxy", "edge_hydroxyl", "carbonyl",
                                         "carboxyl", "ether")}
    placed["paired_hydroxyl_pairs"] = dec.counts["paired"]
    hydroxyl_o = placed["hydroxyl"] + 2 * placed["paired_hydroxyl_pairs"]
    report = {"requested_percent": float(cfg.oxygen_percent), "profile": cfg.oxidation_profile,
              "pattern": cfg.oxidation_pattern, "faces": cfg.oxidation_faces,
              "requested_oxygen_atoms": requested_o_total, "placed_groups": placed,
              "fallback_to_hydroxyl_O": fallback_o,
              "achieved_hydroxyl_to_epoxy": (round(hydroxyl_o / placed["epoxy"], 2)
                                             if placed["epoxy"] else None),
              "face_balance": {"top_or_outer": dec.faces_count[1], "bottom": dec.faces_count[-1]},
              "pyramidalization_A": dec.delta}
    if hasattr(dec, "tiles"):
        report["tiles"] = dec.tiles
    if hasattr(dec, "domains"):
        report["domains"] = dec.domains
    report["warnings"] = dec.warnings
    return report


def _ordered(dec, atoms, order_key, preferred=None):
    atoms = list(atoms)
    perm = dec.rng.permutation(len(atoms))
    atoms = [atoms[i] for i in perm]
    if order_key is not None:
        atoms.sort(key=order_key)
    if preferred is not None:
        atoms.sort(key=lambda a: 0 if a in preferred else 1)
    return atoms


def _place_region(dec: Decorator, members: set, quota: dict, order_key, preferred) -> int:
    """Place one region's quota; returns oxygen atoms redirected to basal hydroxyls."""
    net = dec.net
    fallback = 0

    # 1. edge groups (HierGO aged order: ether, carboxylic acid, carbonyl, edge hydroxyl)
    edges = [e for e in net.edge_atoms() if e in members and dec.edge_ok(e)]
    if quota["ether"]:
        pairs = []
        tree = make_tree(net.pos[edges], net.cell, net.pbc) if edges else None
        if tree is not None:
            for a, b in tree.query_pairs(2.75):
                ea, eb = edges[a], edges[b]
                if ether_pair_ok(net, ea, eb):
                    pairs.append((ea, eb))
        placed = 0
        for k in dec.rng.permutation(len(pairs)):
            if placed >= quota["ether"]:
                break
            a, b = pairs[k]
            if dec.edge_ok(a) and dec.edge_ok(b) and dec.place_ether(a, b):
                placed += 1
        missing = quota["ether"] - placed
        if missing:
            dec.warnings.append(f"{missing} of {quota['ether']} ether oxygens found no vacancy site "
                                "(two edge carbons 2.2-2.75 A apart); added as basal hydroxyls")
            fallback += missing
    for kind, place, o_per in (("carboxyl", dec.place_carboxyl, 2), ("carbonyl", dec.place_carbonyl, 1),
                               ("edge_hydroxyl", dec.place_edge_hydroxyl, 1)):
        want = quota[kind]
        if not want:
            continue
        placed = 0
        for e in _ordered(dec, [e for e in edges if dec.edge_ok(e)], None):
            if placed >= want:
                break
            if dec.edge_ok(e) and place(e):
                placed += 1
        if placed < want:
            name = {"carboxyl": "carboxylic acids", "carbonyl": "carbonyls",
                    "edge_hydroxyl": "edge hydroxyls"}[kind]
            dec.warnings.append(f"{want - placed} of {want} {name} found no free edge carbon; their "
                                "oxygen was added as basal hydroxyls")
            fallback += o_per * (want - placed)

    # 2. basal groups: anti diol pairs, epoxides, then single hydroxyls
    basal = [c for c in members if net.alive[c] and len(net.c_nb(c)) == 3 and len(net.adj[c]) == 3]
    order = _ordered(dec, basal, order_key, preferred)
    rank = {c: k for k, c in enumerate(order)}

    def bonds_in_order():
        out = []
        for c in order:
            for nb in net.c_nb(c):
                if nb in rank and rank[nb] > rank[c]:
                    out.append((c, nb))
        return out

    def place_pairs(want):
        done = 0
        if not want:
            return 0
        for c1, c2 in bonds_in_order():
            if done >= want:
                break
            for f in dec._pref_faces((1, -1), c1):
                if dec.basal_ok(c1, f) and dec.basal_ok(c2, -f) and dec.place_hydroxyl(c1, f):
                    if dec.place_hydroxyl(c2, -f):
                        done += 1
                        dec.counts["hydroxyl"] -= 2
                        dec.counts["paired"] += 1
                        _relabel_last_groups(net, 2, "paired_hydroxyl")
                    else:
                        _undo_last_group(dec, c1, f)
                    break
        return done

    def place_epoxies(want):
        done = 0
        if not want:
            return 0
        for c1, c2 in bonds_in_order():
            if done >= want:
                break
            for f in dec._pref_faces(set(dec.ox_faces.get(net.layer[c1], ()))):
                if dec.basal_ok(c1, f, partner=c2) and dec.basal_ok(c2, f, partner=c1) \
                        and dec.place_epoxy(c1, c2, f):
                    done += 1
                    break
        return done

    def place_hydroxyls(want):
        done = 0
        if not want:
            return 0
        for c in order:
            if done >= want:
                break
            for f in dec._pref_faces(set(dec.ox_faces.get(net.layer[c], ())), c):
                if dec.basal_ok(c, f) and dec.place_hydroxyl(c, f):
                    done += 1
                    break
        return done

    pairs = place_pairs(quota["paired"])
    single = quota["hydroxyl"] + 2 * (quota["paired"] - pairs) + fallback
    if dec.sublattice_bias >= 0.999:
        # Dense oxidation: hydroxyls first in the chair pattern, epoxides where they still
        # fit, then the remaining epoxide oxygen as hydroxyls (hydroxyl-biased, as in the
        # AMBER GO dataset at high coverage).
        placed = place_hydroxyls(single)
        epoxies = place_epoxies(quota["epoxy"])
        want = single - placed + (quota["epoxy"] - epoxies)
        placed_more = place_hydroxyls(want)
        missing_ep = quota["epoxy"] - epoxies
        if missing_ep:
            dec.warnings.append(f"Dense oxidation: {missing_ep} of {quota['epoxy']} epoxides were "
                                "placed as hydroxyls (no free C-C pair left)")
        placed, want = placed_more, want
    else:
        epoxies = place_epoxies(quota["epoxy"])
        if epoxies < quota["epoxy"]:
            missing = quota["epoxy"] - epoxies
            dec.warnings.append(f"{missing} of {quota['epoxy']} epoxides found no free C-C pair; "
                                "added as hydroxyls (AMBER GO workflow fallback)")
            single += missing
        want = single
        placed = place_hydroxyls(want)
    if placed < want:
        converted = _densify(dec, members, want - placed)
        placed += converted
        if converted:
            dec.warnings.append(f"High oxidation: {converted} epoxides were opened into anti "
                                "1,2-diols to reach the oxygen target")
    if placed < want:
        dec.warnings.append(f"Only {placed} of {want} hydroxyls fit: the surface is saturated "
                            "under the placement rules; the achieved oxygen content is lower than "
                            "requested")
    return fallback


def _densify(dec: Decorator, members: set, shortfall: int) -> int:
    """Convert epoxides (1 O on 2 C) into anti diol pairs (2 O on the same 2 C).

    Each conversion adds one oxygen without needing new carbons, which is what
    lets jammed, highly oxidized sheets reach their target. Returns the number
    of conversions (= extra oxygen atoms).
    """
    net = dec.net
    done = 0
    candidates = [k for k, g in enumerate(net.groups)
                  if g["kind"] == "epoxy" and all(c in members for c in g["carbons"])
                  and all(net.alive[a] for a in g["atoms"])]
    for k in dec.rng.permutation(len(candidates)):
        if done >= shortfall:
            break
        g = net.groups[candidates[k]]
        c1, c2 = g["carbons"]
        f = g["face"]
        if -f not in dec.ox_faces.get(net.layer[c2], ()):
            continue
        if any(nb != c1 and dec.face_used.get(nb) == -f for nb in net.c_nb(c2)):
            continue
        _remove_group(dec, g)
        ok1 = dec.place_hydroxyl(c1, f)
        ok2 = ok1 and dec.place_hydroxyl(c2, -f)
        if ok1 and ok2:
            dec.counts["hydroxyl"] -= 2
            dec.counts["paired"] += 1
            _relabel_last_groups(net, 2, "paired_hydroxyl")
            done += 1
            continue
        if ok1:  # roll back the first hydroxyl, restore the epoxide
            _remove_group(dec, net.groups[-1])
            dec.counts["hydroxyl"] -= 1
        dec.place_epoxy(c1, c2, f)
    return done


def _remove_group(dec: Decorator, g):
    """Remove a basal group placed by the decorator and restore its carbons."""
    net = dec.net
    for a in g["atoms"]:
        dec.hash.remove_owner(a)
        if net.alive[a]:
            net.kill(a)
    for c in g["carbons"]:
        if dec.delta:
            net._pos[c] = net.wrap(net.pos[c] - dec.delta * g["face"] * net.normal[c])
        net.role[c] = "C_sp2"
        dec.face_used.pop(c, None)
        for b in net.c_nb(c):
            net.order[(min(b, c), max(b, c))] = 1.5
    dec.faces_count[g["face"]] -= 1
    if g["kind"] == "epoxy":
        dec.counts["epoxy"] -= 1
    g["atoms"] = []
    g["kind"] = "removed"


def _relabel_last_groups(net: CarbonNetwork, n: int, kind: str):
    for g in net.groups[-n:]:
        g["kind"] = kind


def _undo_last_group(dec: Decorator, c, f):
    net = dec.net
    g = net.groups.pop()
    for a in g["atoms"]:
        dec.hash.remove_owner(a)
        net.kill(a)
    net.role[c] = "C_sp2"
    if dec.delta:
        net._pos[c] = net.wrap(net.pos[c] - dec.delta * f * net.normal[c])
    dec.face_used.pop(c, None)
    dec.faces_count[f] -= 1
    dec.counts["hydroxyl"] -= 1
