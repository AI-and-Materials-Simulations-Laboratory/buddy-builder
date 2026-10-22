"""Test matrix 9-10: SMILES handling and the Au-thiolate topology edit."""

import numpy as np
import pytest

from nanobuilder.ligands import anchor_choices, build_ligand, default_anchor, svg_2d
from nanobuilder.schema import BuilderError, LigandConfig


def test_anchor_choices_keep_smiles_order():
    assert anchor_choices("SCCO") == [("0 - S", 0), ("1 - C", 1), ("2 - C", 2), ("3 - O", 3)]
    assert default_anchor("OCCS") == 3  # first thiol sulfur


def test_thiolate_removes_only_the_thiol_hydrogen(scco_thiolate):
    lig = scco_thiolate
    # HS-CH2-CH2-OH: 4 heavy + 6 H = 10 atoms, 9 bonds; thiolate drops the S-H.
    assert lig.n_atoms == 9
    assert len(lig.bonds) == 8
    assert sorted(lig.elements) == sorted(["S", "C", "C", "O"] + ["H"] * 5)
    assert lig.removed_atoms[0]["element"] == "H"
    assert lig.elements[lig.anchor_index] == "S"
    s_neighbors = [j if i == lig.anchor_index else i for i, j, _ in lig.bonds
                   if lig.anchor_index in (i, j)]
    assert [lig.elements[k] for k in s_neighbors] == ["C"]  # S keeps its C, loses its H
    o = lig.elements.index("O")
    assert "H" in [lig.elements[j if i == o else i] for i, j, _ in lig.bonds if o in (i, j)]


def test_geometric_keeps_every_atom():
    lig = build_ligand(LigandConfig(smiles="SCCO", binding_template="geometric"))
    assert lig.n_atoms == 10 and len(lig.bonds) == 9 and not lig.removed_atoms


def test_coordinates_finite_and_anchor_at_origin(scco_thiolate):
    assert np.all(np.isfinite(scco_thiolate.positions_A))
    assert np.allclose(scco_thiolate.positions_A[scco_thiolate.anchor_index], 0)
    assert np.linalg.norm(scco_thiolate.outward_axis) == pytest.approx(1)
    for i, j, _ in scco_thiolate.bonds:
        assert 0.9 < np.linalg.norm(scco_thiolate.positions_A[i] - scco_thiolate.positions_A[j]) < 1.9


def test_index_map_heavy_atoms_unchanged(scco_thiolate):
    for heavy in range(4):
        assert scco_thiolate.index_map[heavy] == heavy


def test_reproducible_with_seed():
    a = build_ligand(LigandConfig(smiles="SCCCCO", seed=7, binding_template="geometric"))
    b = build_ligand(LigandConfig(smiles="SCCCCO", seed=7, binding_template="geometric"))
    assert np.allclose(a.positions_A, b.positions_A)


@pytest.mark.parametrize("smiles", ["C1CC", "not a smiles", "", "CC.O"])
def test_invalid_smiles(smiles):
    with pytest.raises(BuilderError):
        build_ligand(LigandConfig(smiles=smiles, binding_template="geometric"))


def test_unsupported_anchor_template_combinations():
    with pytest.raises(BuilderError, match="No valid Au-S anchor"):
        build_ligand(LigandConfig(smiles="SCCO", anchor_atom_index=3))
    with pytest.raises(BuilderError, match="No valid Au-S anchor"):
        build_ligand(LigandConfig(smiles="CSC", anchor_atom_index=1))  # thioether, no S-H
    with pytest.raises(BuilderError, match="does not exist"):
        build_ligand(LigandConfig(smiles="SCCO", anchor_atom_index=9))


def test_multiple_conformers_and_svg():
    lig = build_ligand(LigandConfig(smiles="SCCCCCCO", conformers=3))
    assert len(lig.conformers_A) == 3
    assert "<svg" in svg_2d("SCCO", highlight=0)


# --- ligand library and extended conformers ---------------------------------

import numpy as _np
import pytest as _pytest

from nanobuilder.ligand_library import LIBRARY, library_choices
from nanobuilder.ligands import build_ligand as _build
from nanobuilder.schema import LigandConfig as _LC


@_pytest.mark.parametrize("key", sorted(LIBRARY))
def test_every_library_ligand_builds(key):
    e = LIBRARY[key]
    lig = _build(_LC(smiles=e["smiles"], anchor_atom_index=e["anchor_atom_index"],
                     binding_template="geometric", conformer_mode=e["conformer_mode"]))
    assert lig.n_atoms > 1 and _np.all(_np.isfinite(lig.positions_A))


def test_library_choices_are_grouped_and_unique():
    labels = [label for label, _ in library_choices()]
    assert len(labels) == len(set(labels)) == len(LIBRARY)
    assert all("·" in label for label in labels)


def _reach(lig):
    p = _np.asarray(lig.positions_A)
    heavy = [i for i, e in enumerate(lig.elements) if e != "H"]
    return _np.linalg.norm(p[heavy] - p[lig.anchor_index], axis=1).max()


def test_extended_mode_straightens_alkyl_chain():
    cfg = dict(smiles="SCCCCCCCCCCCC", anchor_atom_index=0, binding_template="geometric")
    extended = _reach(_build(_LC(conformer_mode="extended", **cfg)))
    # all-trans S...C12 is about 15.4 A (1.26 A per C-C along a zigzag)
    assert 14.8 < extended < 16.0


def test_extended_mode_keeps_cis_double_bond():
    lig = _build(_LC(smiles="NCCCCCCCC/C=C\\CCCCCCCC", anchor_atom_index=0,
                     binding_template="geometric", conformer_mode="extended"))
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms
    # heavy atoms keep SMILES order: C8-C9=C10-C11 are indices 8, 9, 10, 11
    mol = Chem.MolFromSmiles("NCCCCCCCC/C=C\\CCCCCCCC")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, lig.positions_A[lig.index_map[i]].tolist())
    mol.AddConformer(conf)
    assert abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(), 8, 9, 10, 11)) < 30
