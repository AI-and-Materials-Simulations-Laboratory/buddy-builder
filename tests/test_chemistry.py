"""Test matrix 11: the single Au-thiolate availability rule (spec C1)."""

import pytest

from conftest import upload_config
from nanobuilder.chemistry import available_templates, check_template, resolve_surface_species, TEMPLATES
from nanobuilder.crystal import build_crystal
from nanobuilder.importers import atoms_to_cif_bytes
from nanobuilder.schema import BuilderError, MaterialConfig


def keys(crystal):
    return [t.key for t in available_templates(crystal)]


def test_au_preset_offers_thiolate():
    assert keys(build_crystal(MaterialConfig())) == ["Au_thiolate", "chemisorption", "geometric"]
    edited = build_crystal(MaterialConfig(lattice_overrides_A={"a": 4.07}))
    assert "Au_thiolate" in keys(edited)


@pytest.mark.parametrize("preset_id", ["Ag_fcc", "Cu_fcc", "Fe_bcc", "Si_diamond", "NaCl_rocksalt"])
def test_other_presets_get_geometric_only(preset_id):
    crystal = build_crystal(MaterialConfig(preset_id=preset_id))
    assert keys(crystal) == ["geometric", "chemisorption"]
    with pytest.raises(BuilderError, match="only reviewed for the Au/FCC preset"):
        check_template("Au_thiolate", crystal)


def test_uploaded_au_cif_gets_geometric_only():
    au = build_crystal(MaterialConfig())
    data = atoms_to_cif_bytes(au.atoms)
    uploaded = build_crystal(MaterialConfig(source="uploaded_file", preset_id=None, file_format="cif",
                                            filename="gold.cif", uploaded_bytes=data))
    assert keys(uploaded) == ["geometric", "chemisorption"]
    with pytest.raises(BuilderError, match="uploaded crystals"):
        check_template("Au_thiolate", uploaded)


def test_surface_species_resolution(fixtures):
    au = build_crystal(MaterialConfig())
    assert resolve_surface_species(TEMPLATES["Au_thiolate"], None, au) == ["Au"]
    nacl = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    # Unspecified means every element may anchor; an explicitly empty choice is an error.
    assert sorted(resolve_surface_species(TEMPLATES["geometric"], None, nacl)) == ["Cl", "Na"]
    with pytest.raises(BuilderError, match="No surface anchor species selected"):
        resolve_surface_species(TEMPLATES["geometric"], [], nacl)
    assert resolve_surface_species(TEMPLATES["geometric"], ["Na"], nacl) == ["Na"]
    mg = build_crystal(upload_config(fixtures, "mg_hcp_fully_occupied.cif"))
    assert resolve_surface_species(TEMPLATES["geometric"], None, mg) == ["Mg"]
