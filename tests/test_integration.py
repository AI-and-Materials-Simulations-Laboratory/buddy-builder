"""Required live examples A-E end to end, ZIP/provenance (matrix 19), reproducibility (20)."""

import hashlib
import json
import zipfile

import numpy as np
import pytest

from conftest import EXAMPLES_DIR, FIXTURES
from nanobuilder.cli import main as cli_main
from nanobuilder.examples import EXAMPLE_FILES, EXAMPLES, example_upload
from nanobuilder.pipeline import export_package, run_pipeline, validate_inputs
from nanobuilder.schema import BuilderError, RunConfig


def example_config(key):
    upload = example_upload(key)
    return RunConfig.from_dict(EXAMPLES[key]["config"], uploaded_bytes=upload[1] if upload else None)


@pytest.mark.parametrize("key", list(EXAMPLES))
def test_example_end_to_end(tmp_path, key):
    cfg = example_config(key)
    assert validate_inputs(cfg)["ok"]
    result = run_pipeline(cfg)
    assert result.report.status == "ok", result.report.errors
    package = export_package(result, tmp_path)
    assert all(r.round_trip_passed for r in package.reports)
    with zipfile.ZipFile(package.zip_path) as zf:
        names = set(zf.namelist())
        manifest = json.loads(zf.read("manifest.json"))
        for name, digest in manifest["files"].items():
            if name != "manifest.json":
                assert hashlib.sha256(zf.read(name)).hexdigest() == digest, name
    for required in ("README.txt", "config.json", "manifest.json", "validation_report.json",
                     "normalized_unit_cell.cif", "structure.extxyz", "topology.json",
                     "hoomd/structure.gsd", "lammps/structure.data", "lammps/in.lammps.template",
                     "pdb/structure.pdb", "mol2/structure.mol2"):
        assert required in names, required
    if cfg.material.source == "uploaded_file":
        assert "input_crystal.cif" in names
        assert manifest["material"]["uploaded_crystal_sha256"] == hashlib.sha256(
            (FIXTURES / "mg_hcp_fully_occupied.cif").read_bytes()).hexdigest()
    else:
        assert "input_crystal.cif" not in names
        assert manifest["material"]["preset_id"] == cfg.material.preset_id


def test_example_b_numbers():
    result = run_pipeline(example_config("B_au_sphere_12_ligands"))
    a = result.assembly
    assert a.achieved_ligands == 12
    assert a.achieved_density_nm2 == pytest.approx(0.424, abs=1e-3)
    assert sum(1 for b in a.bond_list if b[2] == "Au-S") == 12


def test_example_c_top_only_and_pbc():
    a = run_pipeline(example_config("C_au111_slab_density")).assembly
    assert a.pbc == (True, True, False)
    assert a.target_ligands == round(0.5 * a.surface_area_nm2)
    faces = {s["face"] for s in a.site_to_ligand}
    assert faces == {"top"}


def test_example_d_has_no_thiolate():
    cfg = example_config("D_fe_bcc_sphere")
    cfg.coverage.target = 5
    cfg.ligand.binding_template = "Au_thiolate"
    assert not validate_inputs(cfg)["ok"]
    with pytest.raises(BuilderError, match="only reviewed for the Au/FCC preset"):
        run_pipeline(cfg)


def test_example_json_files_match_module():
    assert set(EXAMPLE_FILES) == set(EXAMPLES)
    for key, fname in EXAMPLE_FILES.items():
        data = json.loads((EXAMPLES_DIR / fname).read_text())
        upload = example_upload(key)
        from_file = RunConfig.from_dict(data, uploaded_bytes=upload[1] if upload else None)
        assert from_file.to_dict() == example_config(key).to_dict(), fname


def test_reproducibility_same_config_and_file():
    a = run_pipeline(example_config("E_uploaded_mg_cif_slab")).assembly
    b = run_pipeline(example_config("E_uploaded_mg_cif_slab")).assembly
    assert np.array_equal(a.positions_A, b.positions_A)
    assert a.bond_list == b.bond_list


def test_config_json_rebuilds_same_model(tmp_path):
    first = run_pipeline(example_config("B_au_sphere_12_ligands"))
    package = export_package(first, tmp_path)
    with zipfile.ZipFile(package.zip_path) as zf:
        saved = json.loads(zf.read("config.json"))
    second = run_pipeline(RunConfig.from_dict(saved))
    assert np.array_equal(first.assembly.positions_A, second.assembly.positions_A)


def test_validation_blocks_bad_exports(tmp_path):
    result = run_pipeline(example_config("A_bare_au_sphere"))
    result.assembly.bond_list.append((0, 0, "bad"))
    from nanobuilder.validation import validate_assembly
    result.report = validate_assembly(result.assembly)
    assert not result.report.ok
    with pytest.raises(BuilderError, match="Export blocked"):
        export_package(result, tmp_path)


def test_cli(tmp_path, capsys):
    assert cli_main([str(EXAMPLES_DIR / "au_sphere_12ligands.json"), "--out", str(tmp_path)]) == 0
    assert list(tmp_path.glob("*.zip"))
    code = cli_main([str(EXAMPLES_DIR / "uploaded_crystal_slab.json"), "--crystal",
                     str(FIXTURES / "mg_hcp_fully_occupied.cif"), "--out", str(tmp_path)])
    assert code == 0
    assert cli_main([str(EXAMPLES_DIR / "uploaded_crystal_slab.json"), "--crystal",
                     str(FIXTURES / "invalid_multiblock.cif"), "--out", str(tmp_path)]) == 1


def test_cli_rebuilds_a_materials_project_config_from_its_cif(tmp_path):
    """config.json + the exported CIF rebuild an MP-sourced model offline."""
    import json
    from ase.build import bulk
    from nanobuilder.cli import main
    from nanobuilder.importers import atoms_to_cif_bytes
    cif = tmp_path / "mp-81.cif"
    cif.write_bytes(atoms_to_cif_bytes(bulk("Au", "fcc", a=4.17, cubic=True)))
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({
        "material": {"source": "materials_project", "mp_id": "mp-81"},
        "geometry": {"shape": "sphere", "diameter_nm": 2.0, "vacuum_A": 15},
        "ligand": {"smiles": "SCCO", "anchor_atom_index": 0, "binding_template": "geometric"},
        "coverage": {"mode": "count", "target": 4}}))
    assert main([str(cfg), "--crystal", str(cif), "--out", str(tmp_path / "out")]) == 0
    zips = list((tmp_path / "out").rglob("*.zip"))
    assert len(zips) == 1 and "mp-81" in zips[0].name
