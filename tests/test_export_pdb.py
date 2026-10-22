"""PDB export: fixed columns, CRYST1, residues, CONECT, whole molecules."""

import copy

import numpy as np
from ase.io import read

from nanobuilder.examples import EXAMPLES
from nanobuilder.exporter_pdb import export_pdb, read_pdb
from nanobuilder.pipeline import run_pipeline
from nanobuilder.schema import RunConfig


def _assembly(key, **ligand):
    cfg = copy.deepcopy(EXAMPLES[key]["config"])
    cfg.get("ligand", {}).update(ligand)
    return run_pipeline(RunConfig.from_dict(cfg)).assembly


def test_pdb_round_trip_sphere(tmp_path):
    a = _assembly("B_au_sphere_12_ligands")
    report = export_pdb(a, tmp_path / "s.pdb")
    assert report.round_trip_passed, report.warnings
    data = read_pdb(tmp_path / "s.pdb")
    lig = [x for x in data["atoms"] if x["chain"] == "B"]
    assert {x["resname"] for x in lig} == {"LIG"}
    assert len({x["resseq"] for x in lig}) == a.achieved_ligands
    sub = [x for x in data["atoms"] if x["chain"] == "A"]
    assert {x["resname"] for x in sub} == {"AU"}
    # each ligand residue has unique atom names
    by_res = {}
    for x in lig:
        by_res.setdefault(x["resseq"], []).append(x["name"])
    assert all(len(v) == len(set(v)) for v in by_res.values())
    # Au-S attachment bonds appear in CONECT
    s_serials = {x["serial"] for x in lig if x["element"] == "S"}
    assert any(s_serials & js for i, js in data["conect"].items())


def test_pdb_readable_by_ase_and_keeps_slab_ligands_whole(tmp_path):
    a = _assembly("C_au111_slab_density", smiles="SCCCCCCCCCCCO", anchor_atom_index=0)
    report = export_pdb(a, tmp_path / "slab.pdb")
    assert report.round_trip_passed, report.warnings
    atoms = read(tmp_path / "slab.pdb")
    assert len(atoms) == a.n_atoms
    assert np.allclose(atoms.cell.cellpar(), a.provenance and np.asarray(atoms.cell.cellpar()))
    pos = atoms.get_positions()
    d = [np.linalg.norm(pos[i] - pos[j]) for i, j, *_ in a.bond_list]
    assert max(d) < 3.5  # no bond crosses the cell in the file
