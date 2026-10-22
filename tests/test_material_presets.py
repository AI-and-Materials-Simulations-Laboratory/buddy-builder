"""Test matrix 1-3: Au/FCC default, alternative presets, registry-only pairs."""

import numpy as np
import pytest

from nanobuilder.crystal import build_crystal, default_preset_id, preset_choices, registry
from nanobuilder.schema import BuilderError, MaterialConfig, RunConfig

EXPECTED = {  # preset: (atoms per conventional cell, a in Angstrom, element)
    "Au_fcc": (4, 4.08, "Au"),
    "Ag_fcc": (4, 4.09, "Ag"),
    "Cu_fcc": (4, 3.61, "Cu"),
    "Fe_bcc": (2, 2.87, "Fe"),
    "Si_diamond": (8, 5.43, "Si"),
    "Ge_diamond": (8, 5.66, "Ge"),
}


def test_default_is_au_fcc():
    assert default_preset_id() == "Au_fcc"
    assert RunConfig().material.preset_id == "Au_fcc"
    crystal = build_crystal(MaterialConfig())
    assert crystal.preset_id == "Au_fcc"
    assert crystal.species == ["Au"]
    assert np.allclose(crystal.cell_A, np.eye(3) * 4.08)


@pytest.mark.parametrize("preset_id", list(EXPECTED))
def test_cubic_presets_are_conventional_cells(preset_id):
    n, a, element = EXPECTED[preset_id]
    crystal = build_crystal(MaterialConfig(preset_id=preset_id))
    assert len(crystal.atoms) == n
    assert crystal.species == [element]
    assert np.allclose(crystal.atoms.cell.lengths(), a)
    assert np.allclose(crystal.atoms.cell.angles(), 90)


def test_fcc_nearest_neighbour_distance():
    crystal = build_crystal(MaterialConfig())
    d = crystal.atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    assert d.min() == pytest.approx(4.08 / np.sqrt(2), abs=1e-6)


def test_hcp_and_rocksalt_presets():
    mg = build_crystal(MaterialConfig(preset_id="Mg_hcp"))
    assert np.allclose(mg.atoms.cell.cellpar(), [3.21, 3.21, 5.21, 90, 90, 120])
    nacl = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    assert nacl.composition == {"Cl": 4, "Na": 4}


def test_unknown_preset_is_rejected():
    with pytest.raises(BuilderError, match="Unsupported material/structure pair"):
        build_crystal(MaterialConfig(preset_id="Au_bcc"))


def test_every_registry_entry_records_lattice_source():
    for key, entry in registry()["presets"].items():
        assert entry["source"], key
        assert set(entry["free_parameters"]) <= set(entry["lattice_A"]), key
    assert [k for _, k in preset_choices()][:5] == ["Au_fcc", "Ag_fcc", "Cu_fcc", "Fe_bcc",
                                                    "Si_diamond"]


def test_crystal_fixture_round_trip(tmp_path):
    from nanobuilder.schema import load_crystal_fixture, save_crystal_fixture

    crystal = build_crystal(MaterialConfig(preset_id="NaCl_rocksalt"))
    save_crystal_fixture(crystal, tmp_path / "nacl")
    back = load_crystal_fixture(tmp_path / "nacl")
    assert back.preset_id == "NaCl_rocksalt" and back.composition == crystal.composition
    assert np.allclose(back.cell_A, crystal.cell_A)
    assert np.allclose(back.atoms.positions, crystal.atoms.positions)


def test_lattice_override_is_recorded_and_bounded():
    crystal = build_crystal(MaterialConfig(lattice_overrides_A={"a": 4.10}))
    assert crystal.lattice_overrides_A == {"a": 4.10}
    assert np.allclose(crystal.atoms.cell.lengths(), 4.10)
    assert any("User-modified" in w for w in crystal.validation_warnings)
    with pytest.raises(BuilderError, match="fixed by the fcc prototype"):
        build_crystal(MaterialConfig(lattice_overrides_A={"c": 5.0}))
    with pytest.raises(BuilderError, match="more than"):
        build_crystal(MaterialConfig(lattice_overrides_A={"a": 8.0}))
    # The registry itself is not modified by an override.
    assert np.allclose(build_crystal(MaterialConfig()).atoms.cell.lengths(), 4.08)
