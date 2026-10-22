"""Test matrix 4-5: CIF happy path and failure paths (no silent Au fallback)."""

import io

import numpy as np
import pytest
from ase.build import bulk
from ase.io import write

from conftest import upload_config
from nanobuilder.crystal import build_crystal
from nanobuilder.importers import MAX_UPLOAD_BYTES, read_crystal_bytes, sanitize_filename
from nanobuilder.schema import BuilderError, MaterialConfig


def test_valid_hcp_cif(fixtures):
    crystal = build_crystal(upload_config(fixtures, "mg_hcp_fully_occupied.cif"))
    assert crystal.source_type == "uploaded_file"
    assert crystal.preset_id is None
    assert crystal.composition == {"Mg": 2}
    assert np.allclose(crystal.atoms.cell.cellpar(), [3.21, 3.21, 5.21, 90, 90, 120], atol=1e-6)
    assert len(crystal.source_sha256) == 64
    assert crystal.source_bytes == (fixtures / "mg_hcp_fully_occupied.cif").read_bytes()


@pytest.mark.parametrize("name, message", [
    ("invalid_partial_occupancy.cif", "Disordered or partially occupied"),
    ("invalid_multiblock.cif", "holds 2 structures"),
    ("invalid_duplicate_atoms.cif", "duplicate sites"),
])
def test_invalid_fixtures_are_rejected(fixtures, name, message):
    with pytest.raises(BuilderError, match=message):
        build_crystal(upload_config(fixtures, name))


def test_garbage_and_lattice_free_input():
    with pytest.raises(BuilderError):
        read_crystal_bytes(b"this is not a crystal file", "x.cif")
    no_cell = b"data_x\nloop_\n_atom_site_label\n_atom_site_type_symbol\n" \
              b"_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nAu1 Au 0 0 0\n"
    with pytest.raises(BuilderError):
        read_crystal_bytes(no_cell, "x.cif")


def test_size_cap_and_empty_upload():
    with pytest.raises(BuilderError, match="empty"):
        read_crystal_bytes(b"", "x.cif")
    with pytest.raises(BuilderError, match="limit"):
        read_crystal_bytes(b"0" * (MAX_UPLOAD_BYTES + 1), "x.cif")


def test_upload_requires_file():
    with pytest.raises(BuilderError, match="Upload a crystallographic file"):
        build_crystal(MaterialConfig(source="uploaded_file", preset_id=None, file_format="cif"))


def test_inconsistent_formula_is_rejected():
    buf = io.BytesIO()
    write(buf, bulk("Cu", "fcc", a=3.61, cubic=True), format="cif")
    text = buf.getvalue().decode().replace('_chemical_formula_sum              "Cu4"',
                                           '_chemical_formula_sum              "Cu Ag"')
    with pytest.raises(BuilderError, match="Inconsistent formula"):
        read_crystal_bytes(text.encode(), "x.cif")


def test_poscar_import_matches_reference(tmp_path):
    reference = bulk("NaCl", "rocksalt", a=5.64, cubic=True)
    path = tmp_path / "POSCAR"
    write(path, reference, format="vasp")
    atoms, info = read_crystal_bytes(path.read_bytes(), "POSCAR", "vasp")
    assert sorted(atoms.get_chemical_symbols()) == sorted(reference.get_chemical_symbols())
    assert np.allclose(atoms.cell, reference.cell)
    assert info["filename"] == "POSCAR.vasp"


def test_filename_sanitizing():
    assert sanitize_filename("../../etc/passwd", "cif") == "passwd.cif"
    assert sanitize_filename("my crystal (1).CIF", "cif") == "my_crystal__1.cif"
