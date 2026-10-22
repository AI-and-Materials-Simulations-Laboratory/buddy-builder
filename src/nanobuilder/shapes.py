"""Nanostructure geometry: spheres, nanorods and oriented slabs (spec stage B)."""

from __future__ import annotations

import warnings
from math import ceil, gcd, pi, sqrt

import numpy as np
from ase import Atoms
from ase.build import surface
from ase.neighborlist import neighbor_list

from .cell import center_boundaries
from .schema import (BuilderError, CrystalModel, FaceInfo, GeometryConfig, SubstrateModel,
                     Termination)

PLANE_TOL_A = 0.05
SHELL_TOL_A = 0.1
SIZE_TOL = 0.01


# --------------------------------------------------------------------------
# Coordination analysis (shared by spheres and slabs)
# --------------------------------------------------------------------------


def neighbor_cutoffs(bulk_atoms: Atoms) -> dict:
    """First-shell cutoff for every species pair, from the bulk crystal.

    Each cutoff is the midpoint between the pair's first and second distance
    shells. Pairs whose first shell is much longer than either species' own
    first shell (e.g. Na-Na in NaCl) are second neighbours and are left out.
    """
    numbers = bulk_atoms.numbers
    rough = np.min(bulk_atoms.get_all_distances(mic=True) + np.eye(len(bulk_atoms)) * 1e9) \
        if len(bulk_atoms) > 1 else min(bulk_atoms.cell.lengths())
    search = max(3.0 * rough, 6.0)
    i, j, d = neighbor_list("ijd", bulk_atoms, search)
    shells = {}
    for a, b, dist in zip(i, j, d):
        key = tuple(sorted((numbers[a], numbers[b])))
        shells.setdefault(key, []).append(dist)
    first = {}
    for key, dists in shells.items():
        dists = np.sort(np.array(dists))
        first[key] = dists[0]
    own_first = {}
    for (za, zb), d1 in first.items():
        for z in (za, zb):
            own_first[z] = min(own_first.get(z, np.inf), d1)
    cutoffs = {}
    for key, dists in shells.items():
        d1 = first[key]
        if d1 > 1.2 * min(own_first[key[0]], own_first[key[1]]):
            continue
        dists = np.sort(np.array(dists))
        beyond = dists[dists > d1 + SHELL_TOL_A]
        cutoffs[key] = 0.5 * (d1 + beyond[0]) if len(beyond) else 1.15 * d1
    return cutoffs


def coordination_numbers(atoms: Atoms, cutoffs: dict) -> np.ndarray:
    if len(atoms) == 0:
        return np.zeros(0, dtype=int)
    i, j, d = neighbor_list("ijd", atoms, max(cutoffs.values()))
    numbers = atoms.numbers
    counts = np.zeros(len(atoms), dtype=int)
    for a, b, dist in zip(i, j, d):
        cut = cutoffs.get(tuple(sorted((numbers[a], numbers[b]))))
        if cut is not None and dist < cut:
            counts[a] += 1
    return counts


def bulk_coordination(bulk_atoms: Atoms, cutoffs: dict) -> dict:
    cn = coordination_numbers(bulk_atoms, cutoffs)
    out = {}
    for symbol, value in zip(bulk_atoms.get_chemical_symbols(), cn):
        out[symbol] = min(out.get(symbol, value), value)
    return out


def under_coordinated(atoms: Atoms, bulk_atoms: Atoms) -> np.ndarray:
    """Indices of atoms with fewer neighbours than the same species in the bulk."""
    cutoffs = neighbor_cutoffs(bulk_atoms)
    reference = bulk_coordination(bulk_atoms, cutoffs)
    cn = coordination_numbers(atoms, cutoffs)
    symbols = atoms.get_chemical_symbols()
    return np.array([k for k, s in enumerate(symbols) if cn[k] < reference[s]], dtype=int)


def _composition(atoms: Atoms) -> dict:
    symbols = atoms.get_chemical_symbols()
    return {s: symbols.count(s) for s in sorted(set(symbols))}


def _reduced(comp: dict) -> dict:
    g = 0
    for v in comp.values():
        g = gcd(g, v)
    return {k: v // g for k, v in comp.items()} if g else comp


def _stoichiometry_warnings(crystal: CrystalModel, comp: dict, shape: str) -> list:
    out = []
    bulk_comp = crystal.composition
    if len(bulk_comp) > 1 and _reduced(comp) != _reduced(bulk_comp):
        out.append(f"The cut {shape} is non-stoichiometric ({comp} vs bulk ratio "
                   f"{_reduced(bulk_comp)}); for ionic crystals this implies a net formal charge. "
                   "No passivation or rebalancing was applied")
    if crystal.prototype == "diamond":
        out.append("Surface atoms of this covalent crystal have dangling bonds; no passivation "
                   "was applied")
    return out


# --------------------------------------------------------------------------
# Public entry point
# --------------------------------------------------------------------------


def build_geometry(crystal: CrystalModel, geometry_config: GeometryConfig) -> SubstrateModel:
    """Cut/replicate/orient the substrate and define cell + PBC."""
    geometry_config.validate()
    if geometry_config.shape == "sphere":
        return _build_sphere(crystal, geometry_config)
    if geometry_config.shape == "rod":
        return _build_rod(crystal, geometry_config)
    if geometry_config.shape == "channel":
        return _build_channel(crystal, geometry_config)
    return _build_slab(crystal, geometry_config)


# --------------------------------------------------------------------------
# Sphere
# --------------------------------------------------------------------------


def _build_sphere(crystal: CrystalModel, g: GeometryConfig) -> SubstrateModel:
    bulk_atoms = crystal.atoms
    radius = g.diameter_nm * 5.0  # nm diameter -> Angstrom radius
    estimate = len(bulk_atoms) / bulk_atoms.get_volume() * 4.0 / 3.0 * pi * radius ** 3
    if estimate > g.max_atoms:
        raise BuilderError(f"A {g.diameter_nm} nm sphere of this material has about "
                           f"{estimate:.0f} atoms; the limit is {g.max_atoms}. Reduce the diameter "
                           "or raise Max atoms")

    cell = np.array(bulk_atoms.cell)
    volume = abs(np.linalg.det(cell))
    widths = [volume / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3]))
              for k in range(3)]
    reps = [int(ceil(2 * radius / w)) + 2 for w in widths]
    sup = bulk_atoms.repeat(reps)
    pos = sup.get_positions()
    # Centre on the atom nearest the supercell centre: atom-centred clusters are symmetric.
    box_centre = np.array(sup.cell).sum(axis=0) / 2
    centre = pos[np.argmin(np.linalg.norm(pos - box_centre, axis=1))]
    keep = np.linalg.norm(pos - centre, axis=1) <= radius + 1e-6
    atoms = sup[keep]
    if len(atoms) == 0:
        raise BuilderError("The sphere is smaller than one atom; increase the diameter")
    atoms.positions -= centre
    rmax = float(np.max(np.linalg.norm(atoms.positions, axis=1)))
    side = 2 * rmax + 2 * g.vacuum_A
    atoms.set_cell(np.eye(3) * side)
    atoms.positions += side / 2
    atoms.pbc = False

    surface_idx = under_coordinated(atoms, bulk_atoms)
    symbols = atoms.get_chemical_symbols()
    exposed = sorted({symbols[k] for k in surface_idx})
    comp = _composition(atoms)
    area_nm2 = 4.0 * pi * (radius / 10.0) ** 2
    faces = {"entire_surface": FaceInfo(atom_indices=surface_idx, exposed_elements=exposed,
                                        termination_label="sphere surface",
                                        area_nm2=area_nm2)}
    warns = list(crystal.validation_warnings) + _stoichiometry_warnings(crystal, comp, "sphere")
    warns.append(f"Nominal radius {radius:.2f} A; outermost atom at {rmax:.2f} A. Coverage uses "
                 f"the nominal sphere area {area_nm2:.3f} nm^2")
    return SubstrateModel(atoms=atoms, crystal=crystal, shape="sphere", pbc=(False,) * 3,
                          faces=faces, composition=comp, center_A=np.full(3, side / 2),
                          radius_A=radius, vacuum_A=g.vacuum_A, warnings=warns)


# --------------------------------------------------------------------------
# Nanorod
# --------------------------------------------------------------------------


def _rotation_onto_z(u: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector ``u`` onto +z."""
    z = np.array([0.0, 0.0, 1.0])
    c = float(np.dot(u, z))
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:
        return np.diag([1.0, -1.0, -1.0])
    v = np.cross(u, z)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k / (1 + c)


def rod_axis_vector(crystal: CrystalModel, axis) -> np.ndarray:
    """Cartesian unit vector of the lattice direction [u v w] of the conventional cell."""
    vec = np.asarray(axis, dtype=float) @ np.array(crystal.atoms.cell)
    return vec / np.linalg.norm(vec)


def _build_rod(crystal: CrystalModel, g: GeometryConfig) -> SubstrateModel:
    """Cylinder of diameter d, length L (tip to tip) along [u v w], hemispherical or flat caps.

    The rod is cut in the crystal frame, then rotated so its axis is +z and
    centred in an orthogonal, vacuum-padded box (finite in every direction).
    Surface atoms are split into the "side" (cylindrical) and "tips" (caps)
    regions so ligands can be placed on either.
    """
    bulk_atoms = crystal.atoms
    radius = g.diameter_nm * 5.0
    length = g.rod_length_nm * 10.0
    hemi = g.rod_caps == "hemispherical"
    half_cyl = length / 2 - radius if hemi else length / 2
    volume = pi * radius ** 2 * 2 * half_cyl + (4.0 / 3.0 * pi * radius ** 3 if hemi else 0.0)
    estimate = len(bulk_atoms) / bulk_atoms.get_volume() * volume
    if estimate > g.max_atoms:
        raise BuilderError(f"A {g.diameter_nm} x {g.rod_length_nm} nm rod of this material has about "
                           f"{estimate:.0f} atoms; the limit is {g.max_atoms}. Make it smaller or "
                           "raise Max atoms")

    u = rod_axis_vector(crystal, g.rod_axis)
    cell = np.array(bulk_atoms.cell)
    cell_volume = abs(np.linalg.det(cell))
    widths = [cell_volume / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3]))
              for k in range(3)]
    # The supercell must hold a sphere around the centre atom that reaches the
    # farthest point of the rod: a cap pole (hemispherical) or a rim corner (flat).
    reach = max(length / 2, radius) if hemi else float(np.hypot(length / 2, radius))
    reps = [int(ceil(2 * reach / w)) + 3 for w in widths]
    sup = bulk_atoms.repeat(reps)
    pos = sup.get_positions()
    box_centre = np.array(sup.cell).sum(axis=0) / 2
    centre = pos[np.argmin(np.linalg.norm(pos - box_centre, axis=1))]
    rel = pos - centre
    t = rel @ u
    radial = np.linalg.norm(rel - np.outer(t, u), axis=1)
    tol = 1e-6
    if hemi:
        t_clamped = np.clip(t, -half_cyl, half_cyl)
        keep = np.linalg.norm(rel - np.outer(t_clamped, u), axis=1) <= radius + tol
    else:
        keep = (np.abs(t) <= half_cyl + tol) & (radial <= radius + tol)
    atoms = sup[keep]
    if len(atoms) == 0:
        raise BuilderError("The rod is smaller than one atom; increase its size")

    rot = _rotation_onto_z(u)
    xyz = (pos[keep] - centre) @ rot.T
    lo, hi = xyz.min(axis=0), xyz.max(axis=0)
    span = np.maximum(np.abs(lo), np.abs(hi))
    sides = 2 * span + 2 * g.vacuum_A
    atoms.positions = xyz + sides / 2
    atoms.set_cell(np.diag(sides))
    atoms.pbc = False
    box_mid = sides / 2

    surface_idx = under_coordinated(atoms, bulk_atoms)
    z = atoms.positions[surface_idx, 2] - box_mid[2]
    if hemi:
        on_tip = np.abs(z) > half_cyl
    else:
        # Flat caps: an atom belongs to a cap when it is closer to that end's
        # outermost atomic plane than to the cylinder wall (measured from the
        # actual atoms, not the nominal cut). Rim atoms (a tie) count as side.
        all_z = atoms.positions[:, 2] - box_mid[2]
        all_r = np.linalg.norm(atoms.positions[:, :2] - box_mid[:2], axis=1)
        r_xy = all_r[surface_idx]
        z_end = np.where(z >= 0, all_z.max(), -all_z.min())
        on_tip = (z_end - np.abs(z)) < (all_r.max() - r_xy)
    side_idx, tip_idx = surface_idx[~on_tip], surface_idx[on_tip]
    side_area = 2 * pi * (radius / 10.0) * (2 * half_cyl / 10.0)
    tip_area = (4 if hemi else 2) * pi * (radius / 10.0) ** 2

    symbols = atoms.get_chemical_symbols()

    def face(idx, label, area):
        return FaceInfo(atom_indices=np.asarray(idx, dtype=int),
                        exposed_elements=sorted({symbols[k] for k in idx}),
                        termination_label=label, area_nm2=area)

    faces = {"side": face(side_idx, "rod side (cylinder)", side_area),
             "tips": face(tip_idx, f"rod tips ({g.rod_caps} caps)", tip_area)}
    comp = _composition(atoms)
    axis_label = "[" + " ".join(str(int(i)) for i in g.rod_axis) + "]"
    warns = list(crystal.validation_warnings) + _stoichiometry_warnings(crystal, comp, "rod")
    warns.append(f"Rod axis {axis_label} of the conventional cell is along +z. Coverage uses the "
                 f"nominal areas: side {side_area:.3f} nm^2, tips {tip_area:.3f} nm^2")
    if hemi and half_cyl <= 0:
        warns.append("Rod length equals its diameter: this is a sphere")
    return SubstrateModel(atoms=atoms, crystal=crystal, shape="rod", pbc=(False,) * 3,
                          faces=faces, composition=comp, center_A=box_mid, radius_A=radius,
                          vacuum_A=g.vacuum_A, warnings=warns,
                          rod_axis=np.array([0.0, 0.0, 1.0]), rod_half_length_A=half_cyl,
                          rod_caps=g.rod_caps, rod_crystal_axis=tuple(int(i) for i in g.rod_axis))


# --------------------------------------------------------------------------
# Slab
# --------------------------------------------------------------------------


def _quiet_surface(bulk_atoms, miller, layers, periodic=False):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return surface(bulk_atoms, tuple(int(i) for i in miller), layers, vacuum=None,
                       periodic=periodic)


def _planes(z: np.ndarray) -> list:
    """Group atom indices into atomic planes by height."""
    order = np.argsort(z)
    groups, current = [], [order[0]]
    for a, b in zip(order[:-1], order[1:]):
        if z[b] - z[a] > PLANE_TOL_A:
            groups.append(current)
            current = []
        current.append(b)
    groups.append(current)
    return [np.array(gr) for gr in groups]


def _plane_elements(atoms: Atoms, idx) -> list:
    symbols = atoms.get_chemical_symbols()
    return sorted(set(symbols[k] for k in idx))


def _oriented_stack(crystal: CrystalModel, miller, units: int):
    stack = _quiet_surface(crystal.atoms, miller, units)
    planes = _planes(stack.positions[:, 2])
    if len(planes) % units:
        raise BuilderError(f"Could not identify the atomic planes of the {tuple(miller)} surface "
                           "for this crystal")
    return stack, planes, len(planes) // units


def planes_per_unit(crystal: CrystalModel, miller) -> int:
    """Atomic planes in one ASE surface repeat unit (e.g. 2 for fcc(100) from the cubic cell)."""
    return _oriented_stack(crystal, miller, 2)[2]


def _plane_window(n: int, layers: int, t: int):
    """Repeat units to stack, and the top/bottom plane indices of a ``layers``-plane slab."""
    units = -(-layers // n) + 2  # ceil(layers / n) + 2 so the window never touches the ends
    top = (units - 1) * n + t
    return units, top, top - layers + 1


def list_terminations(crystal: CrystalModel, miller, layers: int = 4) -> list:
    """Distinct top/bottom terminations for a slab of ``layers`` atomic planes."""
    n = planes_per_unit(crystal, miller)
    units = _plane_window(n, layers, 0)[0]
    stack, planes, _ = _oriented_stack(crystal, miller, units)
    z = np.array([stack.positions[p, 2].mean() for p in planes])
    seen, out = {}, []
    for t in range(n):
        _, top, bottom = _plane_window(n, layers, t)
        gap_bottom = round(z[bottom + 1] - z[bottom], 2) if layers > 1 else 0.0
        key = (tuple(_plane_elements(stack, planes[top])), round(z[top] - z[top - 1], 2),
               tuple(_plane_elements(stack, planes[bottom])), gap_bottom)
        if key in seen:
            continue
        seen[key] = t
        out.append((t, key))
    distinct_tops = len({k[0] for _, k in out}) < len(out)
    terminations = []
    for index, (t, key) in enumerate(out):
        top_el, gap_top, bottom_el, gap_bottom = key
        label = f"Top: {'+'.join(top_el)}"
        if distinct_tops:
            label += f" (next plane {gap_top:.2f} A below)"
        label += f" | Bottom: {'+'.join(bottom_el)}"
        terminations.append(Termination(index=index, top_elements=list(top_el),
                                        bottom_elements=list(bottom_el), label=label,
                                        plane_offset=t))
    return terminations


# --------------------------------------------------------------------------
# Orthogonal (rectangular) in-plane cells for periodic slabs
# --------------------------------------------------------------------------

ORTHO_SEARCH = 8  # integer combinations up to +-8 surface vectors
ORTHO_TOL = 1e-6  # |cos| of the angle between the two new vectors


def orthogonal_supercell(a1, a2, max_n: int = ORTHO_SEARCH):
    """Smallest integer matrix [[i, j], [k, l]] with (i a1 + j a2) perpendicular to (k a1 + l a2).

    Returns None when no exactly perpendicular pair exists within ``max_n``. Both
    new vectors are lattice vectors, so the rectangle tiles the plane exactly:
    atoms outside it are images of atoms inside (trimmed, not strained).
    """
    a1, a2 = np.asarray(a1, dtype=float)[:2], np.asarray(a2, dtype=float)[:2]
    if abs(a1 @ a2) <= ORTHO_TOL * np.linalg.norm(a1) * np.linalg.norm(a2):
        return np.array([[1, 0], [0, 1]])
    rng = range(-max_n, max_n + 1)
    combos = np.array([(i, j) for i in rng for j in rng if (i, j) != (0, 0)])
    vecs = combos[:, :1] * a1 + combos[:, 1:] * a2
    lengths = np.linalg.norm(vecs, axis=1)
    best = None
    for p in np.argsort(lengths):
        dots = vecs @ vecs[p]
        ok = np.abs(dots) <= ORTHO_TOL * lengths * lengths[p] + 1e-9
        for q in np.nonzero(ok)[0]:
            det = combos[p, 0] * combos[q, 1] - combos[p, 1] * combos[q, 0]
            if det <= 0:
                continue
            key = (det, lengths[p] + lengths[q], abs(vecs[p][1]))
            if best is None or key < best[0]:
                best = (key, np.array([combos[p], combos[q]]))
    return None if best is None else best[1]


def _repeats(length_A: float, target_A: float) -> int:
    return max(1, int(ceil(target_A / length_A - SIZE_TOL)))


def _orthogonal_slab(base: Atoms, g: GeometryConfig, warns: list):
    """Slab with a rectangular (or square) in-plane cell; None if the plane has none."""
    a1, a2 = np.array(base.cell[0]), np.array(base.cell[1])
    m = orthogonal_supercell(a1, a2)
    if m is None:
        warns.append("No exactly rectangular repeat of this plane exists within 8 surface cells; "
                      "the native (skewed) cell is used")
        return None
    v1 = m[0, 0] * a1 + m[0, 1] * a2
    v2 = m[1, 0] * a1 + m[1, 1] * a2
    l1, l2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
    box = _user_box_A(g)
    if box is not None:  # the box decides the size: nearest whole number of repeats
        tx, ty = box[0], box[1]
        n1, n2 = max(1, int(round(tx / l1))), max(1, int(round(ty / l2)))
    else:
        tx, ty = g.slab_length_nm * 10.0, g.slab_width_nm * 10.0
        n1, n2 = _repeats(l1, tx), _repeats(l2, ty)
    if g.box_shape == "cube":  # square in-plane cell: closest match of n1 l1 and n2 l2
        side = max(tx, ty)
        lo1 = max(1, int(round(side / l1)) - 3) if box is not None else _repeats(l1, side)
        lo2 = max(1, int(round(side / l2)) - 3) if box is not None else _repeats(l2, side)
        best = None
        for c1 in range(lo1, lo1 + 7):
            for c2 in range(lo2, lo2 + 7):
                mismatch = abs(c1 * l1 - c2 * l2) / max(c1 * l1, c2 * l2)
                size = abs(c1 * l1 - side) + abs(c2 * l2 - side) if box is not None else \
                    c1 * l1 * c2 * l2
                key = (0.0 if mismatch <= 0.01 else mismatch, size)
                if best is None or key < best[0]:
                    best = (key, c1, c2)
        n1, n2 = best[1], best[2]
    estimate = int(len(base) * abs(np.linalg.det(m)) * n1 * n2)
    if estimate > g.max_atoms:
        raise BuilderError(f"This slab would have about {estimate} atoms; the limit is "
                           f"{g.max_atoms}. Use a smaller box or fewer layers")
    from ase.build import make_supercell
    prim = base.copy()
    cell = np.array(prim.cell)
    cell[2] = [0.0, 0.0, 1000.0]
    prim.set_cell(cell)
    prim.pbc = (True, True, False)
    P = np.array([[n1 * m[0, 0], n1 * m[0, 1], 0], [n2 * m[1, 0], n2 * m[1, 1], 0], [0, 0, 1]])
    slab = make_supercell(prim, P, wrap=True)
    theta = np.arctan2(v1[1], v1[0])
    c, s_ = np.cos(-theta), np.sin(-theta)
    rot = np.array([[c, -s_, 0.0], [s_, c, 0.0], [0.0, 0.0, 1.0]])
    slab.positions = slab.positions @ rot.T
    lx, ly = n1 * l1, n2 * l2
    slab.set_cell(np.diag([lx, ly, 1000.0]))
    slab.pbc = (True, True, False)
    slab.wrap()
    shape = "square" if g.box_shape == "cube" else "rectangular"
    asked = (f"box {tx / 10:.2f} x {ty / 10:.2f} nm" if box is not None else
             f"requested {g.slab_length_nm} x {g.slab_width_nm} nm")
    warns.append(f"Rectangular periodic cell: {lx / 10:.3f} x {ly / 10:.3f} nm ({n1} x {n2} "
                 f"repeats of a {l1:.3f} x {l2:.3f} A {shape} unit; {asked}). The crystal "
                 "repeats exactly across the box edges, so x and y are the nearest whole number "
                 "of repeats")
    return slab


def _user_box_A(g: GeometryConfig):
    """User box lengths (x, y, z) in A, or None for an automatic box."""
    if g.box_lengths_nm is None or not any(float(x) > 0 for x in g.box_lengths_nm):
        return None
    return tuple(10.0 * float(x) for x in g.box_lengths_nm)


def crystal_frame_normal(crystal: CrystalModel, miller) -> np.ndarray:
    recip = np.array(crystal.atoms.cell.reciprocal())
    n = np.asarray(miller, dtype=float) @ recip
    return n / np.linalg.norm(n)


def _build_slab(crystal: CrystalModel, g: GeometryConfig) -> SubstrateModel:
    miller = tuple(int(i) for i in g.miller)
    try:
        unit = _quiet_surface(crystal.atoms, miller, 1, periodic=True)
    except Exception as exc:
        raise BuilderError(f"Could not build the {miller} surface: {exc}") from None
    height = float(unit.cell[2][2])  # repeat distance along the surface normal
    n = planes_per_unit(crystal, miller)
    spacing = height / n  # mean spacing of atomic planes
    if g.thickness_mode == "layers":  # a layer is one atomic plane
        layers = int(g.layers)
    else:
        layers = max(1, int(round(g.thickness_nm * 10.0 / spacing)))

    a1, a2 = np.array(unit.cell[0]), np.array(unit.cell[1])
    area_unit = float(np.linalg.norm(np.cross(a1, a2)))
    w1, w2 = area_unit / np.linalg.norm(a2), area_unit / np.linalg.norm(a1)
    # At least the requested size, allowing a 1% shortfall so 49.97 A counts as 5 nm.
    nx = max(1, int(ceil(g.slab_length_nm * 10.0 / w1 - SIZE_TOL)))
    ny = max(1, int(ceil(g.slab_width_nm * 10.0 / w2 - SIZE_TOL)))
    estimate = int(len(unit) * layers / n * nx * ny)
    if estimate > g.max_atoms and _user_box_A(g) is None and g.box_shape == "auto":
        raise BuilderError(f"This slab would have about {estimate} atoms; the limit is "
                           f"{g.max_atoms}. Use a smaller size or fewer layers")

    terminations = list_terminations(crystal, miller, layers)
    if g.termination_index >= len(terminations):
        raise BuilderError(f"Termination {g.termination_index} does not exist; this surface has "
                           f"{len(terminations)}")
    term = terminations[g.termination_index]

    units, top, bottom = _plane_window(n, layers, term.plane_offset)
    stack, planes, _ = _oriented_stack(crystal, miller, units)
    keep = np.concatenate(planes[bottom:top + 1])
    base = stack[np.sort(keep)]
    box_notes = []
    ortho = g.box_shape in ("rectangular", "cube") or _user_box_A(g) is not None
    slab = _orthogonal_slab(base, g, box_notes) if ortho else None
    if slab is None:
        slab = base.repeat((nx, ny, 1))
    center_boundaries(slab)  # straight cuts at the periodic box edges
    z = slab.positions[:, 2]
    thickness = float(z.max() - z.min())
    slab.positions[:, 2] += g.vacuum_A - z.min()
    cell = np.array(slab.cell)
    cell[2] = [0.0, 0.0, thickness + 2 * g.vacuum_A]
    slab.set_cell(cell)
    slab.pbc = (True, True, False)
    slab.wrap()

    surface_idx = under_coordinated(slab, crystal.atoms)
    z = slab.positions[:, 2]
    mid = 0.5 * (z.max() + z.min())
    face_area = float(np.linalg.norm(np.cross(cell[0], cell[1]))) / 100.0  # A^2 -> nm^2
    symbols = slab.get_chemical_symbols()
    top_idx = surface_idx[z[surface_idx] >= mid]
    bottom_idx = surface_idx[z[surface_idx] < mid]
    faces = {
        "top": FaceInfo(top_idx, sorted(set(symbols[k] for k in top_idx)),
                        "+".join(term.top_elements), face_area),
        "bottom": FaceInfo(bottom_idx, sorted(set(symbols[k] for k in bottom_idx)),
                           "+".join(term.bottom_elements), face_area),
    }
    comp = _composition(slab)
    warns = list(crystal.validation_warnings) + _stoichiometry_warnings(crystal, comp, "slab")
    if len(terminations) > 1:
        warns.append(f"This surface has {len(terminations)} distinct terminations; using "
                     f"'{term.label}'. Terminations are not assumed equivalent, neutral or stable")
    if set(term.top_elements) != set(term.bottom_elements):
        warns.append("Top and bottom faces terminate in different elements (polar or asymmetric "
                     "slab)")
    if box_notes:
        warns += box_notes
    else:
        warns.append(f"In-plane cell {nx} x {ny} surface units: |a| = {np.linalg.norm(cell[0]):.2f} A, "
                     f"|b| = {np.linalg.norm(cell[1]):.2f} A, area per face {face_area:.3f} nm^2 "
                     f"(requested {g.slab_length_nm} x {g.slab_width_nm} nm)")
    return SubstrateModel(
        atoms=slab, crystal=crystal, shape="slab", pbc=(True, True, False), faces=faces,
        composition=comp, miller=miller, surface_normal=np.array([0.0, 0.0, 1.0]),
        crystal_frame_normal=crystal_frame_normal(crystal, miller),
        interplanar_spacing_A=spacing, layers=layers, thickness_A=thickness,
        termination_index=term.index, terminations=terminations, vacuum_A=g.vacuum_A,
        warnings=warns)


# --------------------------------------------------------------------------
# Channels: stacked slit pores between slab walls, or cylindrical pores through a membrane
# --------------------------------------------------------------------------


def _build_channel(crystal: CrystalModel, g: GeometryConfig) -> SubstrateModel:
    from dataclasses import replace
    wall = _build_slab(crystal, replace(g, shape="slab"))
    if g.channel_type == "slit":
        return _slit_channel(crystal, g, wall)
    return _pore_channel(crystal, g, wall)


def _slit_channel(crystal: CrystalModel, g: GeometryConfig, wall: SubstrateModel) -> SubstrateModel:
    n_channels = int(g.channel_count)
    gap = float(g.channel_width_nm) * 10.0
    base = wall.atoms
    z = base.positions[:, 2]
    zmin, thickness = float(z.min()), float(z.max() - z.min())
    if gap < 2.0:
        raise BuilderError("Slit width must be at least 0.2 nm (distance between facing atomic planes)")
    n_walls = n_channels + 1
    total = n_walls * len(base)
    if total > g.max_atoms:
        raise BuilderError(f"This channel would have {total} atoms; the limit is {g.max_atoms}. "
                           "Reduce the size, wall layers or channel count, or raise Max atoms")
    parts, wall_id = [], []
    for i in range(n_walls):
        copy = base.copy()
        copy.positions[:, 2] += g.vacuum_A - zmin + i * (thickness + gap)
        parts.append(copy)
        wall_id += [i] * len(copy)
    atoms = parts[0]
    for part in parts[1:]:
        atoms += part
    span = n_walls * thickness + n_channels * gap
    cell = np.array(base.cell)
    cell[2] = [0.0, 0.0, span + 2 * g.vacuum_A]
    atoms.set_cell(cell)
    atoms.pbc = (True, True, False)
    wall_id = np.array(wall_id)

    surface_idx = under_coordinated(atoms, crystal.atoms)
    zz = atoms.positions[:, 2]
    normals = np.zeros((len(atoms), 3))
    inner, top, bottom = [], [], []
    for k in surface_idx:
        w = int(wall_id[k])
        mid = g.vacuum_A + w * (thickness + gap) + thickness / 2
        up = zz[k] >= mid
        normals[k] = [0.0, 0.0, 1.0 if up else -1.0]
        if (up and w == n_walls - 1):
            top.append(k)
        elif (not up and w == 0):
            bottom.append(k)
        else:
            inner.append(k)
    symbols = atoms.get_chemical_symbols()
    face_area = wall.faces["top"].area_nm2

    def info(idx, label, area):
        idx = np.array(idx, dtype=int)
        return FaceInfo(idx, sorted(set(symbols[k] for k in idx)), label, area)

    faces = {"inner": info(inner, "channel walls (pore side)", 2 * n_channels * face_area),
             "top": info(top, "outer top face", face_area),
             "bottom": info(bottom, "outer bottom face", face_area)}
    warns = [w_ for w_ in wall.warnings]
    warns.append(f"Slit channel: {n_channels} channel(s) {gap / 10:.2f} nm wide between "
                 f"{n_walls} walls of {wall.layers} atomic planes ({thickness / 10:.2f} nm); "
                 "the width is measured between facing atomic planes")
    return SubstrateModel(
        atoms=atoms, crystal=crystal, shape="channel", pbc=(True, True, False), faces=faces,
        composition=_composition(atoms), miller=wall.miller,
        surface_normal=np.array([0.0, 0.0, 1.0]), crystal_frame_normal=wall.crystal_frame_normal,
        interplanar_spacing_A=wall.interplanar_spacing_A, layers=wall.layers,
        thickness_A=span, termination_index=wall.termination_index,
        terminations=wall.terminations, vacuum_A=g.vacuum_A, warnings=warns, normals=normals,
        channel_info={"type": "slit", "width_nm": gap / 10, "channels": n_channels,
                      "walls": n_walls, "wall_thickness_nm": round(thickness / 10, 4)})


PORE_MARGIN_SPACINGS = 3  # pores keep >= 3 atomic spacings from each other and the box edges
INV_PHI = (sqrt(5.0) - 1.0) / 2.0


def _nearest_spacing(atoms: Atoms) -> float:
    from scipy.spatial import cKDTree
    pos = atoms.get_positions()
    sample = pos[:: max(1, len(pos) // 400)]
    d, _ = cKDTree(pos).query(sample, k=2)
    return float(np.median(d[:, 1]))


def pore_centers(lx: np.ndarray, ly: np.ndarray, n: int, radius: float, margin: float) -> list:
    """N pore centres (Cartesian xy) spread evenly over the membrane.

    Starts from a golden-ratio lattice (x_i = (i + 1/2)/N, y_i = frac(i / phi),
    the planar Fibonacci pattern) inside the band that keeps every pore
    ``margin`` from the box edges, then relaxes the points apart. Raises a
    BuilderError when the pores cannot keep ``margin`` from each other and the
    edges.
    """
    basis = np.array([lx[:2], ly[:2]], dtype=float)
    area = abs(np.linalg.det(basis))
    widths = [area / np.linalg.norm(basis[1]), area / np.linalg.norm(basis[0])]
    need = 2 * radius + 2 * margin
    if min(widths) < need - 1e-6:
        raise BuilderError(
            f"A {2 * radius / 10:.2f} nm pore needs the box to be at least {need / 10:.2f} nm wide "
            f"(pore plus {PORE_MARGIN_SPACINGS} atomic spacings to each box edge); the membrane "
            f"is {min(widths) / 10:.2f} nm. Enlarge the box or use a smaller pore")
    lo = np.array([(radius + margin) / w for w in widths])
    hi = 1.0 - lo
    if n == 1:
        frac = np.array([[0.5, 0.5]])
    else:
        i = np.arange(n, dtype=float)
        frac = np.stack([(i + 0.5) / n, (i * INV_PHI + 0.5 / n) % 1.0], axis=1)
        frac = lo + frac * (hi - lo)
        target = 2 * radius + margin
        spacing = max(target, sqrt(area / n))
        for _ in range(400):  # soft repulsion inside the allowed band
            xy = frac @ basis
            d = xy[:, None, :] - xy[None, :, :]
            r = np.maximum(np.linalg.norm(d, axis=2), 1e-3) + np.eye(n)
            push = np.where(r < spacing, (spacing - r) / r, 0.0)
            np.fill_diagonal(push, 0.0)
            move = 0.25 * (push[:, :, None] * d).sum(axis=1)
            frac = np.clip(frac + np.linalg.solve(basis.T, move.T).T, lo, hi)
    xy = frac @ basis
    if n > 1:
        d = np.linalg.norm(xy[:, None] - xy[None], axis=2) + np.eye(n) * 1e9
        if d.min() < 2 * radius + margin - 1e-6:
            raise BuilderError(
                f"{n} pores of {2 * radius / 10:.2f} nm do not fit with {PORE_MARGIN_SPACINGS} "
                f"atomic spacings ({margin / 10:.2f} nm) between pores and to the box edges; "
                "enlarge the box, or use fewer or smaller pores")
    return [np.array([p[0], p[1], 0.0]) for p in xy]


def _pore_channel(crystal: CrystalModel, g: GeometryConfig, wall: SubstrateModel) -> SubstrateModel:
    atoms = wall.atoms.copy()
    cell = np.array(atoms.cell)
    radius = float(g.channel_width_nm) * 5.0
    n_pores = int(g.channel_count)
    lx, ly = cell[0], cell[1]
    spacing = _nearest_spacing(atoms)
    margin = PORE_MARGIN_SPACINGS * spacing
    if 2 * radius < spacing:
        raise BuilderError(f"The pore ({2 * radius / 10:.2f} nm) is smaller than the atomic "
                           f"spacing ({spacing / 10:.2f} nm): no atoms would be removed")
    centers = pore_centers(lx, ly, n_pores, radius, margin)
    pos = atoms.get_positions()
    keep = np.ones(len(atoms), dtype=bool)
    rel_all = []
    for c in centers:
        d = pos[:, :2] - c[:2]
        frac = np.linalg.solve(np.array([lx[:2], ly[:2]]).T, d.T).T
        frac -= np.round(frac)
        d = frac @ np.array([lx[:2], ly[:2]])
        rel_all.append(d)
        keep &= np.linalg.norm(d, axis=1) >= radius
    if keep.all():
        raise BuilderError("The pore is smaller than the atomic spacing: no atoms were removed")
    atoms = atoms[keep]
    rel_all = [r[keep] for r in rel_all]
    surface_idx = under_coordinated(atoms, crystal.atoms)
    z = atoms.positions[:, 2]
    mid = 0.5 * (z.max() + z.min())
    normals = np.zeros((len(atoms), 3))
    inner, top, bottom = [], [], []
    for k in surface_idx:
        dists = [np.linalg.norm(r[k]) for r in rel_all]
        j = int(np.argmin(dists))
        if dists[j] < radius + 2.5:
            v = -rel_all[j][k]
            normals[k] = [v[0] / dists[j], v[1] / dists[j], 0.0]
            inner.append(k)
        elif z[k] >= mid:
            normals[k] = [0.0, 0.0, 1.0]
            top.append(k)
        else:
            normals[k] = [0.0, 0.0, -1.0]
            bottom.append(k)
    symbols = atoms.get_chemical_symbols()
    thickness = float(z.max() - z.min())
    face_area = wall.faces["top"].area_nm2 - n_pores * pi * (radius / 10) ** 2

    def info(idx, label, area):
        idx = np.array(idx, dtype=int)
        return FaceInfo(idx, sorted(set(symbols[k] for k in idx)), label, area)

    faces = {"inner": info(inner, "pore wall", n_pores * 2 * pi * radius / 10 * thickness / 10),
             "top": info(top, "membrane top face", face_area),
             "bottom": info(bottom, "membrane bottom face", face_area)}
    warns = list(wall.warnings)
    warns.append(f"Cylindrical pore channel: {n_pores} pore(s) of {2 * radius / 10:.2f} nm "
                 f"through a {thickness / 10:.2f} nm membrane, spread evenly over the surface at "
                 f"least {margin / 10:.2f} nm ({PORE_MARGIN_SPACINGS} atomic spacings) from each "
                 f"other and the box edges ({len(keep) - int(keep.sum())} atoms removed); "
                 "pore-wall atoms carry inward normals")
    return SubstrateModel(
        atoms=atoms, crystal=crystal, shape="channel", pbc=(True, True, False), faces=faces,
        composition=_composition(atoms), miller=wall.miller,
        surface_normal=np.array([0.0, 0.0, 1.0]), crystal_frame_normal=wall.crystal_frame_normal,
        interplanar_spacing_A=wall.interplanar_spacing_A, layers=wall.layers,
        thickness_A=thickness, termination_index=wall.termination_index,
        terminations=wall.terminations, vacuum_A=g.vacuum_A, warnings=warns, normals=normals,
        channel_info={"type": "pore", "width_nm": 2 * radius / 10, "channels": n_pores,
                      "centers_A": [[round(float(c[0]), 3), round(float(c[1]), 3)] for c in centers]})
