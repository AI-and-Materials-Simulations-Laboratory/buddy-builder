"""Canonical topology (spec E): global offsets, attachment bonds, typed tuples."""

import numpy as np

from nanobuilder.placement import bare_assembly, place_ligands
from nanobuilder.schema import AttachmentConfig, CoverageConfig
from nanobuilder.sites import find_surface_sites
from nanobuilder.topology import ligand_angles_dihedrals, parameter_coverage
from nanobuilder.validation import validate_assembly


def build(substrate, ligand, n=5):
    sites = find_surface_sites(substrate, AttachmentConfig(surface_species=["Au"]))
    return place_ligands(substrate, ligand, sites, CoverageConfig(target=n))


def test_atom_order_and_molecule_ids(au_sphere, scco_thiolate):
    a = build(au_sphere, scco_thiolate)
    n_sub = len(au_sphere.atoms)
    assert a.component[:n_sub] == ["substrate"] * n_sub
    assert set(a.molecule_id[:n_sub]) == {0}
    for k in range(5):
        block = slice(n_sub + 9 * k, n_sub + 9 * (k + 1))
        assert set(a.molecule_id[block]) == {k + 1}
        assert a.elements[block] == scco_thiolate.elements


def test_bonds_are_offset_per_copy(au_sphere, scco_thiolate):
    a = build(au_sphere, scco_thiolate)
    n_sub = len(au_sphere.atoms)
    ligand_bonds = [b for b in a.bond_list if b[2] != "Au-S"]
    assert len(ligand_bonds) == 5 * len(scco_thiolate.bonds)
    for i, j, _ in ligand_bonds:
        assert (i - n_sub) // 9 == (j - n_sub) // 9  # never between ligand copies
    assert sum(1 for b in a.bond_list if b[2] == "Au-S") == 5
    assert not any(a.elements[i] == "Au" and a.elements[j] == "Au" for i, j, _ in a.bond_list)


def test_angles_and_dihedrals(scco_thiolate):
    angles, dihedrals = ligand_angles_dihedrals(scco_thiolate)
    # S(1 nb) C(S,C,2H) C(C,O,2H) O(C,H): angles C1 6 + C2 6 + O 1 = 13; dihedrals over C-C, C-S, C-O
    assert len(angles) == 13
    assert len(dihedrals) == 12


def test_type_labels_are_structural(au_sphere, scco_thiolate):
    a = build(au_sphere, scco_thiolate)
    assert a.atom_types is None
    assert {"Au_sub", "S_anchor", "C", "H", "O"} == set(a.type_labels)
    coverage = parameter_coverage(a)
    assert coverage["atom_types_assigned"] is False
    assert "Au-S" in coverage["missing"]["bond_types"]


def test_box_encloses_ligands_with_vacuum(au_sphere, scco_thiolate):
    a = build(au_sphere, scco_thiolate, n=12)
    pos = a.positions_A
    assert pos.min() >= 30.0 - 1e-6
    assert (np.diag(a.cell_A) - pos.max(axis=0)).min() >= 30.0 - 1e-6
    assert validate_assembly(a).status == "ok"


def test_bare_assembly(au_sphere):
    a = bare_assembly(au_sphere)
    report = validate_assembly(a)
    assert report.ok and a.achieved_ligands == 0 and not a.bond_list
