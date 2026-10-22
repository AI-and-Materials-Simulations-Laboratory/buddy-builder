"""Exposed surface sites (spec stage D1). P0 supports on-top sites only."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .schema import CARBON_SHAPES, AttachmentConfig, BuilderError, SiteCatalog, SubstrateModel

# A site is accessible when no other substrate atom is closer to the anchor
# probe point than this fraction of the anchor distance.
PROBE_CLEARANCE_FRACTION = 0.9


def region_faces(substrate: SubstrateModel, region: str) -> list:
    if substrate.shape in CARBON_SHAPES:
        from .carbon import region_face_names, region_kinds
        region_kinds(substrate.shape, region)  # raises a readable error for a wrong region
        return region_face_names(substrate.shape, region)
    if substrate.shape == "sphere":
        if region not in ("entire_surface",):
            raise BuilderError("A sphere has one surface region: choose 'Entire surface'")
        return ["entire_surface"]
    if substrate.shape == "rod":
        mapping = {"entire_surface": ["side", "tips"], "side": ["side"], "tips": ["tips"]}
        if region not in mapping:
            raise BuilderError("A rod has the regions 'Entire surface', 'Side' and 'Tips'")
        return mapping[region]
    if substrate.shape == "channel":
        mapping = {"inner": ["inner"], "top": ["top"], "bottom": ["bottom"],
                   "both": ["top", "bottom"], "entire_surface": ["inner", "top", "bottom"]}
        if region not in mapping:
            raise BuilderError("A channel has the regions 'Channel walls', 'Top', 'Bottom', "
                               "'Both outer faces' and 'Entire surface'")
        return mapping[region]
    mapping = {"top": ["top"], "bottom": ["bottom"], "both": ["top", "bottom"],
               "entire_surface": ["top", "bottom"]}
    if region not in mapping:
        raise BuilderError(f"Unknown placement region {region!r}")
    return mapping[region]


def selected_area_nm2(substrate: SubstrateModel, region: str) -> float:
    return float(sum(substrate.faces[f].area_nm2 for f in region_faces(substrate, region)))


def periodic_images(positions: np.ndarray, cell: np.ndarray, pbc) -> tuple:
    """Positions plus periodic images (-1, 0, +1 along every periodic axis).

    Returns (image_positions, original_index). Slabs are periodic in-plane,
    periodic nanotubes along z; finite objects get no images.
    """
    positions = np.asarray(positions)
    if not any(pbc):
        return positions, np.arange(len(positions))
    cell = np.asarray(cell, dtype=float)
    shifts = [np.zeros(3)]
    for k in range(3):
        if pbc[k]:
            shifts = [s + m * cell[k] for s in shifts for m in (-1, 0, 1)]
    images = np.concatenate([positions + s for s in shifts])
    index = np.tile(np.arange(len(positions)), len(shifts))
    return images, index


def rod_normal(substrate: SubstrateModel, point: np.ndarray, face: str):
    """Outward normal of a spherocylinder / flat-capped cylinder at ``point``."""
    axis = substrate.rod_axis
    rel = point - substrate.center_A
    t = float(rel @ axis)
    h = float(substrate.rod_half_length_A)
    if substrate.rod_caps == "flat" and face == "tips":
        return axis * (1.0 if t >= 0 else -1.0)
    if substrate.rod_caps == "hemispherical":
        t = min(max(t, -h), h)  # nearest point on the axis segment
    # (flat caps, side face: purely radial, t unchanged)
    normal = rel - t * axis
    length = np.linalg.norm(normal)
    return None if length < 1e-6 else normal / length


def find_surface_sites(substrate: SubstrateModel, attachment_config: AttachmentConfig) -> SiteCatalog:
    """Return eligible site positions, surface atom IDs, and outward normals."""
    if substrate.shape in CARBON_SHAPES:
        from .carbon import carbon_site_catalog
        return carbon_site_catalog(substrate, attachment_config.template or "geometric",
                                   attachment_config.region)
    if attachment_config.site_type != "on_top":
        raise BuilderError("Only on-top adsorption sites are available in this version")
    atoms = substrate.atoms
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    allowed = set(attachment_config.surface_species)
    images, image_index = periodic_images(positions, np.array(atoms.cell), substrate.pbc)
    tree = cKDTree(images)
    d = float(attachment_config.distance_A)

    rows, warnings = [], []
    for face in region_faces(substrate, attachment_config.region):
        info = substrate.faces[face]
        candidates = [k for k in info.atom_indices if symbols[k] in allowed]
        rejected = 0
        for k in candidates:
            if substrate.shape == "sphere":
                normal = positions[k] - substrate.center_A
                length = np.linalg.norm(normal)
                if length < 1e-6:
                    continue
                normal = normal / length
            elif substrate.shape == "rod":
                normal = rod_normal(substrate, positions[k], face)
                if normal is None:
                    continue
            elif substrate.normals is not None:
                normal = np.asarray(substrate.normals[k], dtype=float)
                if np.linalg.norm(normal) < 1e-6:
                    continue
            else:
                normal = np.array([0.0, 0.0, 1.0 if face == "top" else -1.0])
            probe = positions[k] + d * normal
            near = [image_index[n] for n in tree.query_ball_point(probe, PROBE_CLEARANCE_FRACTION * d)]
            if any(n != k for n in near):
                rejected += 1
                continue
            rows.append((positions[k], k, normal, face, symbols[k]))
        if rejected:
            warnings.append(f"{face}: {rejected} under-coordinated atoms rejected as buried or "
                            "inaccessible for on-top placement")

    if not rows:
        raise BuilderError("0 valid exposed sites: no accessible surface atom of the selected "
                           f"species ({', '.join(sorted(allowed)) or 'none'}) in the selected region")
    pos, idx, normals, faces, elements = zip(*rows)
    area = {f: substrate.faces[f].area_nm2 for f in region_faces(substrate, attachment_config.region)}
    warnings.append("Site exposure is geometric; it does not establish that binding is chemically "
                    "feasible")
    return SiteCatalog(positions_A=np.array(pos), surface_atom_index=np.array(idx, dtype=int),
                       normals=np.array(normals), face=list(faces), element=list(elements),
                       site_type=["on_top"] * len(rows), area_nm2_by_face=area, warnings=warnings)
