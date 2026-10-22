"""Coverage targets, site selection, ligand orientation and clash rejection (spec D2-D5)."""

from __future__ import annotations

import numpy as np
from ase.data import atomic_numbers, covalent_radii, vdw_radii
from ase.geometry import find_mic
from scipy.spatial import cKDTree

from .schema import AssemblyModel, CoverageConfig, LigandModel, SiteCatalog, SubstrateModel
from .sites import periodic_images
from .topology import assemble

N_AZIMUTHS = 12
# Placed-ligand atoms enter a KD-tree in batches; the newest ones are checked
# directly against a small array. Rebuilding the tree (with 9 periodic images on
# slabs) after every ligand made dense packing quadratic in the ligand count.
TREE_REBUILD_EVERY = 32


def vdw_radius(element: str) -> float:
    """Bondi van der Waals radius; covalent radius + 0.8 A where Bondi has none (e.g. Fe)."""
    z = atomic_numbers[element]
    r = vdw_radii[z]
    return float(r) if np.isfinite(r) else float(covalent_radii[z] + 0.8)


def target_count(coverage: CoverageConfig, area_nm2: float) -> int:
    if coverage.mode == "count":
        return int(coverage.target)
    return int(round(coverage.density_nm2 * area_nm2))


def rotation_from_to(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector u onto unit vector v."""
    u, v = u / np.linalg.norm(u), v / np.linalg.norm(v)
    c = float(np.dot(u, v))
    if c > 1 - 1e-10:
        return np.eye(3)
    if c < -1 + 1e-10:
        perp = np.cross(u, [1.0, 0.0, 0.0])
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(u, [0.0, 1.0, 0.0])
        return rotation_about(perp / np.linalg.norm(perp), np.pi)
    w = np.cross(u, v)
    k = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
    return np.eye(3) + k + k @ k / (1 + c)


def rotation_about(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def orient(ligand_xyz: np.ndarray, outward_axis: np.ndarray, site: np.ndarray, normal: np.ndarray,
           distance: float, azimuth: float) -> np.ndarray:
    """Anchor at ``site + distance * normal``, ligand axis along the normal."""
    rot = rotation_about(normal, azimuth) @ rotation_from_to(outward_axis, normal)
    return ligand_xyz @ rot.T + site + distance * normal


class _ClashChecker:
    def __init__(self, substrate: SubstrateModel, ligand: LigandModel, alpha: float):
        self.alpha = alpha
        self.pbc = tuple(substrate.pbc)
        self.cell = np.array(substrate.atoms.cell)
        sub_pos = substrate.atoms.get_positions()
        self.sub_r = np.array([vdw_radius(e) for e in substrate.atoms.get_chemical_symbols()])
        images, self.image_index = periodic_images(sub_pos, self.cell, self.pbc)
        self.sub_tree = cKDTree(images)
        self.lig_r = np.array([vdw_radius(e) for e in ligand.elements])
        self.anchor = ligand.anchor_index
        self.max_r = max(self.sub_r.max(), self.lig_r.max())
        self.placed = np.zeros((0, 3))
        self.placed_r = np.zeros(0)
        self.placed_tree = None
        self.placed_index = None
        self.tree_count = 0            # placed atoms already inside placed_tree
        self.recent_pos = np.zeros((0, 3))  # newer placed atoms (with periodic images)
        self.recent_r = np.zeros(0)
        self.since_rebuild = 0
        self.min_distance = np.inf

    def check(self, xyz: np.ndarray, site_atom: int, ignore=(), site_neighbors=(),
              anchor_shell=()):
        """Return (ok, closest nonbonded distance) for a trial ligand copy.

        Exempt pairs: the anchor and its site atom; for covalent attachments also
        substrate atoms the attachment removes or moves (``ignore``), the site
        atom against the anchor's bonded neighbours (``anchor_shell``) and the
        site's bonded neighbours against the anchor (1-3 pairs across the bond).
        """
        closest = np.inf
        for a, point in enumerate(xyz):
            reach = self.alpha * (self.lig_r[a] + self.max_r)
            for n in self.sub_tree.query_ball_point(point, reach):
                j = self.image_index[n]
                if a == self.anchor and j == site_atom:
                    continue  # the intended attachment pair is exempt
                if j in ignore or (j == site_atom and a in anchor_shell) or \
                        (a == self.anchor and j in site_neighbors):
                    continue
                d = np.linalg.norm(point - self.sub_tree.data[n])
                # 1-4 pairs across a new covalent bond (anchor neighbours vs site neighbours)
                # are bonded-term territory: only reject real overlaps.
                factor = 0.5 if (a in anchor_shell and j in site_neighbors) else self.alpha
                if factor == self.alpha:
                    closest = min(closest, d)
                if d < factor * (self.lig_r[a] + self.sub_r[j]):
                    return False, d
            if self.placed_tree is not None:
                for n in self.placed_tree.query_ball_point(point, reach):
                    j = self.placed_index[n]
                    d = np.linalg.norm(point - self.placed_tree.data[n])
                    closest = min(closest, d)
                    if d < self.alpha * (self.lig_r[a] + self.placed_r[j]):
                        return False, d
        if len(self.recent_r):
            # Only recent atoms inside the trial ligand's bounding box (plus the
            # largest reach) can matter; filter first, then compute distances.
            pad = self.alpha * (self.lig_r.max() + self.max_r)
            lo, hi = xyz.min(axis=0) - pad, xyz.max(axis=0) + pad
            box = np.all((self.recent_pos >= lo) & (self.recent_pos <= hi), axis=1)
            if not np.any(box):
                return True, closest
            recent_pos, recent_r = self.recent_pos[box], self.recent_r[box]
            d = np.linalg.norm(xyz[:, None, :] - recent_pos[None, :, :], axis=2)
            reach = self.alpha * (self.lig_r[:, None] + self.max_r)
            limit = self.alpha * (self.lig_r[:, None] + recent_r[None, :])
            near = d <= reach
            if np.any(near):
                if np.any(d[near] < limit[near]):
                    return False, float(d[near & (d < limit)].min())
                closest = min(closest, float(d[near].min()))
        return True, closest

    def commit(self, xyz: np.ndarray, closest: float):
        self.placed = np.vstack([self.placed, xyz])
        self.placed_r = np.concatenate([self.placed_r, self.lig_r])
        self.since_rebuild += 1
        if self.since_rebuild >= TREE_REBUILD_EVERY:
            images, self.placed_index = periodic_images(self.placed, self.cell, self.pbc)
            self.placed_tree = cKDTree(images)
            self.tree_count = len(self.placed)
            self.recent_pos, self.recent_r = np.zeros((0, 3)), np.zeros(0)
            self.since_rebuild = 0
        else:
            images, index = periodic_images(xyz, self.cell, self.pbc)
            self.recent_pos = np.vstack([self.recent_pos, images])
            self.recent_r = np.concatenate([self.recent_r, self.lig_r[index]])
        self.min_distance = min(self.min_distance, closest)


def _site_distances(sites: SiteCatalog, chosen: list, k: int, cell, pbc) -> np.ndarray:
    if not chosen:
        return np.array([np.inf])
    delta = sites.positions_A[chosen] - sites.positions_A[k]
    if any(pbc):
        delta, dist = find_mic(delta, cell, pbc)
        return dist
    return np.linalg.norm(delta, axis=1)


def _distances_to(sites: SiteCatalog, c: int, cell, pbc) -> np.ndarray:
    """Distance from every site to site ``c`` (minimum image along periodic axes)."""
    delta = sites.positions_A - sites.positions_A[c]
    if any(pbc):
        _, dist = find_mic(delta, cell, pbc)
        return dist
    return np.linalg.norm(delta, axis=1)


def _bonded_neighbours(substrate: SubstrateModel) -> dict:
    out = {}
    for i, j, *_ in substrate.bonds:
        out.setdefault(i, set()).add(j)
        out.setdefault(j, set()).add(i)
    return out


def place_ligands(substrate: SubstrateModel, ligand: LigandModel, sites: SiteCatalog,
                  coverage_config: CoverageConfig, surface_area_nm2: float | None = None,
                  extra_warnings: list | None = None) -> AssemblyModel:
    """Select sites, orient copies, reject clashes, preserve bond graph."""
    coverage_config.validate()
    area = float(surface_area_nm2 if surface_area_nm2 is not None
                 else sum(sites.area_nm2_by_face.values()))
    n_target = target_count(coverage_config, area)
    requested_density = (coverage_config.density_nm2 if coverage_config.mode == "density"
                         else (n_target / area if area > 0 else None))
    warnings = list(extra_warnings or []) + list(sites.warnings)
    if n_target > len(sites):
        warnings.append(f"Requested {n_target} ligands but only {len(sites)} eligible sites "
                        "exist; at most one ligand per site can be placed")

    rng = np.random.default_rng(int(coverage_config.seed))
    order = list(rng.permutation(len(sites)))
    cell, pbc = np.array(substrate.atoms.cell), tuple(substrate.pbc)
    checker = _ClashChecker(substrate, ligand, float(coverage_config.clash_alpha))
    conformers = ligand.conformers_A or [ligand.positions_A]
    axes = ligand.outward_axes if len(ligand.outward_axes) == len(conformers) else \
        [ligand.outward_axis] * len(conformers)
    min_sep = float(coverage_config.min_separation_A)
    # Covalent substrates: each site consumes substrate atoms (claims); exempt pairs
    # across the new bond from the clash test.
    sub_nb = _bonded_neighbours(substrate) if substrate.bonds else {}
    lig_nb = {ligand.anchor_index}
    for i, j, _ in ligand.bonds:
        if i == ligand.anchor_index:
            lig_nb.add(j)
        elif j == ligand.anchor_index:
            lig_nb.add(i)
    used_atoms = set()
    # Bonded attachment on a crystal (Au-thiolate, chemisorption): the site atom and the
    # anchor's own neighbours are a 1-3 pair across the new bond, so they are exempt too.
    from .chemistry import TEMPLATES
    bonded = TEMPLATES[ligand.binding_template].covalent_surface_bond

    chosen, placements = [], []
    state = {"attempts": 0, "too_close": 0}

    def try_site(k) -> bool:
        """One site: claims, separation, then conformers x azimuths until one fits."""
        claims = sites.claims[k] if sites.claims is not None else {int(sites.surface_atom_index[k])}
        if sites.claims is not None and claims & used_atoms:
            return False  # another ligand already consumed these substrate atoms
        if min_sep > 0 and chosen and _site_distances(sites, chosen, k, cell, pbc).min() < min_sep:
            state["too_close"] += 1
            return False
        start = 0.0 if ligand.orientation == "outward" else rng.uniform(0, 2 * np.pi)
        site, normal = sites.positions_A[k], sites.normals[k]
        edits = sites.edits[k] if sites.edits is not None else {}
        site_atom = int(sites.surface_atom_index[k])
        ignore = set(edits.get("ignore", ()))
        exempt = dict(ignore=ignore, site_neighbors=sub_nb.get(site_atom, set()),
                      anchor_shell=lig_nb) if sub_nb else (
            dict(anchor_shell=lig_nb) if bonded else {})
        for conformer, axis in zip(conformers, axes):
            for step in range(N_AZIMUTHS):
                state["attempts"] += 1
                xyz = orient(conformer, axis, site, normal, ligand.distance_A,
                             start + 2 * np.pi * step / N_AZIMUTHS)
                ok, closest = checker.check(xyz, site_atom, **exempt)
                if ok:
                    checker.commit(xyz, closest)
                    chosen.append(k)
                    used_atoms.update(claims)
                    placements.append({"positions": xyz, "site_index": int(k),
                                       "surface_atom_index": site_atom,
                                       "face": sites.face[k], "edits": edits})
                    return True
                if state["attempts"] >= coverage_config.max_attempts:
                    return False
        return False

    def fill(remaining: list, evenly: bool) -> list:
        # Distance from every site to its nearest chosen site, updated as sites are
        # chosen (same choice as recomputing it against all chosen sites each step).
        nearest_chosen = np.full(len(sites), np.inf)
        n_chosen_seen = 0
        while (len(placements) < n_target and remaining
               and state["attempts"] < coverage_config.max_attempts):
            if evenly and chosen:
                for c in chosen[n_chosen_seen:]:
                    nearest_chosen = np.minimum(nearest_chosen, _distances_to(sites, c, cell, pbc))
                n_chosen_seen = len(chosen)
                # Rounded so exact lattice ties go to the first site in the seeded order,
                # not to whichever value floating-point noise makes larger.
                k = remaining[int(np.argmax(np.round(nearest_chosen[remaining], 6)))]
            else:
                k = remaining[0]
            remaining = [r for r in remaining if r != k]
            try_site(k)
        return remaining

    fib_summary = None
    if coverage_config.distribution == "fibonacci":
        from .fibonacci import plan_targets
        targets, fib_summary = plan_targets(substrate, sites, n_target)
        tried = set()
        for cands in targets:
            if len(placements) >= n_target or state["attempts"] >= coverage_config.max_attempts:
                break
            for k in cands:
                k = int(k)
                if k in tried:
                    continue
                tried.add(k)
                if try_site(k):
                    break
        fib_summary["placed_on_targets"] = len(placements)
        # Points whose neighbourhood was full: fill the gaps farthest from placed ligands.
        remaining = fill([k for k in order if k not in tried], evenly=True)
    else:
        remaining = fill(order, coverage_config.distribution == "evenly_spaced")
    attempts, too_close = state["attempts"], state["too_close"]

    if len(placements) < n_target:
        reasons = []
        if too_close:
            reasons.append(f"{too_close} sites skipped for the minimum separation")
        if attempts >= coverage_config.max_attempts:
            reasons.append(f"the attempt limit ({coverage_config.max_attempts}) was reached")
        if not remaining:
            reasons.append("all eligible sites were tried")
        warnings.append(f"Only {len(placements)} of {n_target} ligands could be placed without "
                        f"clashes ({'; '.join(reasons) or 'no clash-free sites left'})")

    assembly = assemble(substrate, ligand, placements, n_target, area, requested_density,
                        warnings, seeds={"ligand_conformer": None, "placement": coverage_config.seed})
    assembly.min_nonbonded_distance_A = (None if not np.isfinite(checker.min_distance)
                                         else float(checker.min_distance))
    assembly.placements = placements
    assembly.provenance["placement"] = {
        "distribution": coverage_config.distribution,
        "clash_alpha": coverage_config.clash_alpha,
        "min_separation_A": min_sep,
        "attempts": attempts,
        "eligible_sites": len(sites),
        "azimuth_trials_per_conformer": N_AZIMUTHS,
        "clash_rule": "reject if r_ij < alpha * (vdW_i + vdW_j); bonded/attachment pair exempt",
    }
    if fib_summary is not None:
        assembly.provenance["placement"]["fibonacci"] = fib_summary
    return assembly


def bare_assembly(substrate: SubstrateModel, warnings: list | None = None) -> AssemblyModel:
    """Substrate-only canonical model (no ligands requested)."""
    if substrate.shape == "sphere":
        area = substrate.faces["entire_surface"].area_nm2
    elif substrate.shape == "rod":
        area = sum(f.area_nm2 for f in substrate.faces.values())
    elif substrate.shape == "nanotube":
        area = substrate.faces["side"].area_nm2
    elif substrate.shape == "carbon_sphere":
        area = substrate.faces["entire_surface"].area_nm2
    elif substrate.shape == "channel":
        area = substrate.faces["inner"].area_nm2
    else:
        area = substrate.faces["top"].area_nm2
    return assemble(substrate, None, [], 0, area, None, list(warnings or []), seeds={})
