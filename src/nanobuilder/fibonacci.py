"""Ordered ("Fibonacci") ligand placement targets.

For each surface region the target ligand count is spread over evenly placed
points, and every point is matched to the closest free candidate site:

- spheres, fullerenes and other compact particles: a Fibonacci (golden-spiral)
  lattice on the sphere, z_i = 1 - (2i + 1)/N, phi_i = i x 137.508 deg
  (Gonzalez, Math. Geosci. 42, 49, 2010), matched by direction from the centre;
- cylinder walls (rod sides, nanotubes, pore walls): a golden-angle helix,
  z evenly spaced along the axis;
- flat faces (slabs, sheets, slit channels): a golden-ratio lattice in the face,
  x_i = (i + 0.5)/N, y_i = frac(i / phi), periodic where the cell is;
- anything else (rod tips, edges): deterministic farthest-point order.

The result is the same for every run and seed. Each target carries a short
list of candidate sites, nearest first, so a clash at the best site falls back
to a neighbour no farther than half the mean spacing.
"""

from __future__ import annotations

from math import pi, sqrt

import numpy as np

from .schema import SiteCatalog, SubstrateModel

GOLDEN_ANGLE = pi * (3.0 - sqrt(5.0))  # 137.508 degrees
INV_PHI = (sqrt(5.0) - 1.0) / 2.0
MAX_CANDIDATES = 8
FALLBACK_FRACTION = 0.5  # of the mean target spacing


def fibonacci_sphere(n: int) -> np.ndarray:
    """N unit vectors on a Fibonacci lattice."""
    i = np.arange(n, dtype=float)
    z = 1.0 - (2.0 * i + 1.0) / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    phi = i * GOLDEN_ANGLE
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def _split(counts: dict, n: int) -> dict:
    """Largest-remainder split of n over weights."""
    total = sum(counts.values())
    if total <= 0 or n <= 0:
        return {k: 0 for k in counts}
    raw = {k: n * v / total for k, v in counts.items()}
    out = {k: int(np.floor(v)) for k, v in raw.items()}
    for k in sorted(raw, key=lambda k: raw[k] - out[k], reverse=True)[: n - sum(out.values())]:
        out[k] += 1
    return out


def _mic(d: np.ndarray, cell: np.ndarray, pbc) -> np.ndarray:
    d = np.array(d, dtype=float)
    for k in range(3):
        if pbc[k] and cell[k, k] > 0:
            d[..., k] -= cell[k, k] * np.round(d[..., k] / cell[k, k])
    return d


def _candidates(dist: np.ndarray, spacing: float) -> np.ndarray:
    order = np.argsort(dist)[:MAX_CANDIDATES]
    keep = [order[0]] + [k for k in order[1:] if dist[k] <= max(FALLBACK_FRACTION * spacing,
                                                                  dist[order[0]])]
    return np.array(keep, dtype=int)


def _subgroups(substrate: SubstrateModel, sites: SiteCatalog, idx: np.ndarray, face=None) -> list:
    """Split a face's sites into planes when their normals are +-z (slit channels, slabs),
    or by pore for the walls of a multi-pore membrane."""
    info = substrate.channel_info or {}
    if face == "inner" and info.get("type") == "pore" and info.get("centers_A"):
        centers = np.array(info["centers_A"], dtype=float)
        cell = np.array(substrate.atoms.cell)
        d = np.zeros((len(idx), len(centers), 3))
        d[:, :, :2] = sites.positions_A[idx, None, :2] - centers[None, :, :2]
        d = _mic(d, cell, (True, True, False))
        owner = np.argmin(np.linalg.norm(d[:, :, :2], axis=2), axis=1)
        return [(("pore", centers[k]), idx[owner == k]) for k in range(len(centers))
                if np.any(owner == k)]
    nz = sites.normals[idx, 2]
    if np.all(np.abs(np.abs(nz) - 1.0) < 1e-3):
        z = sites.positions_A[idx, 2]
        groups = []
        for sign in (1.0, -1.0):
            sel = idx[np.sign(nz) == sign]
            if not len(sel):
                continue
            zs = sites.positions_A[sel, 2]
            order = np.argsort(zs)
            start = 0
            for k in range(1, len(order) + 1):
                if k == len(order) or zs[order[k]] - zs[order[k - 1]] > 3.0:
                    groups.append(("plane", sel[order[start:k]]))
                    start = k
        del z
        return groups
    return [(None, idx)]


def _kind(substrate: SubstrateModel, face: str, hint) -> str:
    if hint == "plane":
        return "plane"
    if isinstance(hint, tuple) and hint[0] == "pore":
        return "cylinder"
    if substrate.shape in ("sphere", "carbon_sphere"):
        return "sphere"
    if face in ("side", "inner") and substrate.shape in ("rod", "nanotube"):
        return "cylinder"
    if face == "inner" and substrate.channel_info.get("type") == "pore":
        return "cylinder"
    if face == "entire_surface" and not any(substrate.pbc):
        return "sphere"
    return "farthest"


def _sphere_targets(pos, n, center):
    d = pos - center
    d /= np.maximum(np.linalg.norm(d, axis=1), 1e-9)[:, None]
    u = fibonacci_sphere(n)
    spacing = sqrt(4 * pi / n)
    return [(_candidates(np.linalg.norm(d - t, axis=1), spacing), float(np.degrees(
        np.arccos(np.clip(d[np.argmin(np.linalg.norm(d - t, axis=1))] @ t, -1, 1)))))
        for t in u], "deg"


def _cylinder_targets(pos, n, cell, pbc, center=None):
    xy = pos[:, :2]
    if center is None:
        c = xy.mean(axis=0)
        rel = xy - c
    else:  # one pore of a periodic membrane: nearest-image offsets from its axis
        d = np.zeros((len(pos), 3))
        d[:, :2] = xy - np.asarray(center)[:2]
        rel = _mic(d, cell, (bool(pbc[0]), bool(pbc[1]), False))[:, :2]
    radius = float(np.mean(np.linalg.norm(rel, axis=1)))
    phi = np.arctan2(rel[:, 1], rel[:, 0])
    z = pos[:, 2]
    zmin, zmax = float(z.min()), float(z.max())
    length = max(zmax - zmin, 1e-6)
    periodic_z = bool(pbc[2])
    if periodic_z:
        length, zmin = float(cell[2, 2]), 0.0
    spacing = sqrt(2 * pi * radius * length / n)
    out = []
    for i in range(n):
        zt = zmin + (i + 0.5) / n * length
        pt = (i * GOLDEN_ANGLE) % (2 * pi)
        dphi = (phi - pt + pi) % (2 * pi) - pi
        dz = z - zt
        if periodic_z:
            dz -= length * np.round(dz / length)
        dist = np.hypot(radius * dphi, dz)
        out.append((_candidates(dist, spacing), float(dist.min())))
    return out, "A"


def _plane_targets(pos, n, cell, pbc):
    periodic = bool(pbc[0] and pbc[1])
    if periodic:
        a, b = np.array(cell[0][:2]), np.array(cell[1][:2])
        origin = np.zeros(2)
    else:
        lo, hi = pos[:, :2].min(axis=0), pos[:, :2].max(axis=0)
        a, b = np.array([hi[0] - lo[0], 0.0]), np.array([0.0, hi[1] - lo[1]])
        origin = lo
    area = abs(a[0] * b[1] - a[1] * b[0])
    spacing = sqrt(max(area, 1e-6) / n)
    out = []
    for i in range(n):
        f = np.array([(i + 0.5) / n, (i * INV_PHI) % 1.0])
        t = origin + f[0] * a + f[1] * b
        d = np.zeros((len(pos), 3))
        d[:, :2] = pos[:, :2] - t
        if periodic:
            d = _mic(d, np.array(cell), (True, True, False))
        dist = np.linalg.norm(d[:, :2], axis=1)
        out.append((_candidates(dist, spacing), float(dist.min())))
    return out, "A"


def _farthest_targets(pos, n, cell, pbc):
    m = len(pos)
    if m == 0 or n == 0:
        return [], "A"
    start = int(np.argmin(np.linalg.norm(pos - pos.mean(axis=0), axis=1)))
    chosen = [start]
    best = np.linalg.norm(_mic(pos - pos[start], cell, pbc), axis=1)
    while len(chosen) < min(n, m):
        k = int(np.argmax(np.round(best, 6)))
        chosen.append(k)
        best = np.minimum(best, np.linalg.norm(_mic(pos - pos[k], cell, pbc), axis=1))
    out = []
    spacing = float(np.median(best)) if len(chosen) > 1 else 1.0
    for k in chosen:
        dist = np.linalg.norm(_mic(pos - pos[k], cell, pbc), axis=1)
        out.append((_candidates(dist, max(spacing, 1.0)), 0.0))
    return out, "A"


def plan_targets(substrate: SubstrateModel, sites: SiteCatalog, n_target: int):
    """[(candidate site indices nearest first, offset)], plus a summary dict."""
    cell = np.array(substrate.atoms.cell)
    pbc = tuple(substrate.pbc)
    faces = np.array(sites.face)
    groups = []
    for face in dict.fromkeys(sites.face):
        idx = np.nonzero(faces == face)[0]
        for hint, sub in _subgroups(substrate, sites, idx, face):
            groups.append((face, hint, sub))
    area = sites.area_nm2_by_face or {}
    weights = {}
    for g, (face, _, sub) in enumerate(groups):
        face_sites = int(np.sum(faces == face))
        face_area = float(area.get(face, 0.0)) or face_sites
        weights[g] = face_area * len(sub) / max(face_sites, 1)
    counts = _split(weights, int(n_target))
    targets, offsets, units = [], [], set()
    center = (np.array(substrate.center_A) if substrate.center_A is not None
              else substrate.atoms.get_positions().mean(axis=0))
    for g, (face, hint, sub) in enumerate(groups):
        n = counts[g]
        if n <= 0:
            continue
        pos = sites.positions_A[sub]
        kind = _kind(substrate, face, hint)
        if kind == "sphere":
            plan, unit = _sphere_targets(pos, n, center)
        elif kind == "cylinder":
            centre = hint[1] if isinstance(hint, tuple) else None
            plan, unit = _cylinder_targets(pos, n, cell, pbc, centre)
        elif kind == "plane":
            plan, unit = _plane_targets(pos, n, cell, pbc)
        else:
            plan, unit = _farthest_targets(pos, n, cell, pbc)
        units.add(unit)
        for cand, off in plan:
            targets.append(sub[cand])
            offsets.append(off)
    # interleave groups would bias nothing: targets are tried in this order
    summary = {"targets": len(targets), "groups": len(groups)}
    if offsets:
        unit = "deg" if units == {"deg"} else "A"
        summary.update({f"mean_offset_{unit}": round(float(np.mean(offsets)), 3),
                        f"max_offset_{unit}": round(float(np.max(offsets)), 3)})
    return targets, summary


__all__ = ["fibonacci_sphere", "plan_targets", "GOLDEN_ANGLE"]
