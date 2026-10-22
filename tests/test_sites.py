"""Exposed-site detection (spec D1)."""

import numpy as np
import pytest

from nanobuilder.crystal import build_crystal
from nanobuilder.schema import AttachmentConfig, BuilderError, GeometryConfig, MaterialConfig
from nanobuilder.shapes import build_geometry, under_coordinated
from nanobuilder.sites import find_surface_sites


def test_sphere_sites_point_outward(au_sphere):
    sites = find_surface_sites(au_sphere, AttachmentConfig(surface_species=["Au"]))
    assert len(sites) > 100
    radial = sites.positions_A - au_sphere.center_A
    radial /= np.linalg.norm(radial, axis=1)[:, None]
    assert np.allclose(np.einsum("ij,ij->i", radial, sites.normals), 1)
    assert set(sites.face) == {"entire_surface"}


def test_slab_top_only(au111_slab):
    sites = find_surface_sites(au111_slab, AttachmentConfig(region="top", surface_species=["Au"]))
    z = au111_slab.atoms.positions[:, 2]
    assert np.allclose(sites.positions_A[:, 2], z.max())
    assert np.allclose(sites.normals, [0, 0, 1])
    both = find_surface_sites(au111_slab, AttachmentConfig(region="both", surface_species=["Au"]))
    assert len(both) == 2 * len(sites)


def test_fcc110_open_surface_exposes_top_row_and_trough():
    au = build_crystal(MaterialConfig())
    slab = build_geometry(au, GeometryConfig(shape="slab", miller=(1, 1, 0), layers=8,
                                             slab_length_nm=2, slab_width_nm=2))
    z = slab.atoms.positions[:, 2]
    planes = np.unique(z.round(3))
    under = under_coordinated(slab.atoms, au.atoms)
    # Top two and bottom two planes are under-coordinated (CN 7 and 11); the middle is bulk (12).
    assert set(np.unique(z[under].round(3))) == {planes[0], planes[1], planes[-2], planes[-1]}
    sites = find_surface_sites(slab, AttachmentConfig(region="top", surface_species=["Au"]))
    assert set(np.unique(sites.positions_A[:, 2].round(3))) == {planes[-2], planes[-1]}


def test_buried_sphere_atoms_are_rejected(au_sphere):
    d = 2.4
    sites = find_surface_sites(au_sphere, AttachmentConfig(surface_species=["Au"], distance_A=d))
    under = set(au_sphere.faces["entire_surface"].atom_indices)
    rejected = under - set(sites.surface_atom_index)
    assert rejected  # some under-coordinated atoms sit in hollows of the faceted particle
    pos = au_sphere.atoms.positions
    for k in rejected:
        normal = (pos[k] - au_sphere.center_A) / np.linalg.norm(pos[k] - au_sphere.center_A)
        probe = pos[k] + d * normal
        others = np.delete(np.linalg.norm(pos - probe, axis=1), k)
        assert others.min() < 0.9 * d


def test_species_filter_and_errors():
    nacl = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    slab = build_geometry(nacl, GeometryConfig(shape="slab", miller=(1, 0, 0), slab_length_nm=2,
                                               slab_width_nm=2))
    na = find_surface_sites(slab, AttachmentConfig(region="top", surface_species=["Na"]))
    assert set(na.element) == {"Na"}
    with pytest.raises(BuilderError, match="0 valid exposed sites"):
        find_surface_sites(slab, AttachmentConfig(region="top", surface_species=["Au"]))
    with pytest.raises(BuilderError, match="one surface region"):
        find_surface_sites(build_geometry(nacl, GeometryConfig(diameter_nm=1.5)),
                           AttachmentConfig(region="top", surface_species=["Na"]))
