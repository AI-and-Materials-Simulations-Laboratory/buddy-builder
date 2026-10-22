"""SMILES -> 3D ligand with anchor bookkeeping (spec stage C)."""

from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors, rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D

from .chemistry import TEMPLATES, AttachmentTemplate, template_distance
from .schema import BuilderError, LigandConfig, LigandModel

RDLogger.DisableLog("rdApp.*")

_BOND_ORDERS = {
    Chem.BondType.SINGLE: 1.0,
    Chem.BondType.DOUBLE: 2.0,
    Chem.BondType.TRIPLE: 3.0,
    Chem.BondType.AROMATIC: 1.5,
}


def parse_smiles(smiles: str) -> Chem.Mol:
    smiles = str(smiles or "").strip()
    if not smiles:
        raise BuilderError("Enter a ligand SMILES string")
    if "." in smiles:
        raise BuilderError("Invalid SMILES: enter one connected molecule (no '.')")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise BuilderError(f"Invalid SMILES: RDKit could not parse {smiles!r}")
    if mol.GetNumAtoms() == 0:
        raise BuilderError("Invalid SMILES: no atoms")
    return mol


def anchor_choices(smiles: str) -> list:
    """(label, index) for every heavy atom, in SMILES order, e.g. ('0 - S', 0)."""
    mol = parse_smiles(smiles)
    return [(f"{a.GetIdx()} - {a.GetSymbol()}", a.GetIdx()) for a in mol.GetAtoms()]


def is_thiol_sulfur(mol: Chem.Mol, index: int) -> bool:
    atom = mol.GetAtomWithIdx(index)
    return (atom.GetSymbol() == "S" and atom.GetTotalNumHs() == 1 and atom.GetDegree() == 1
            and atom.GetFormalCharge() == 0 and not atom.GetIsAromatic())


def default_anchor(smiles: str) -> int:
    """First thiol sulfur if any, else atom 0."""
    mol = parse_smiles(smiles)
    for atom in mol.GetAtoms():
        if is_thiol_sulfur(mol, atom.GetIdx()):
            return atom.GetIdx()
    return 0


def ligand_summary(smiles: str) -> dict:
    mol = parse_smiles(smiles)
    molh = Chem.AddHs(mol)
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "heavy_atoms": mol.GetNumAtoms(),
        "atoms_with_h": molh.GetNumAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "thiol_sulfurs": [a.GetIdx() for a in mol.GetAtoms() if is_thiol_sulfur(mol, a.GetIdx())],
    }


def svg_2d(smiles: str, size=(320, 220), highlight=None) -> str:
    """2D drawing with atom indices, for anchor selection in the GUI."""
    mol = parse_smiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.drawOptions().addAtomIndices = True
    kwargs = {}
    if highlight is not None and 0 <= highlight < mol.GetNumAtoms():
        kwargs = {"highlightAtoms": [int(highlight)]}
    drawer.DrawMolecule(mol, **kwargs)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def check_anchor(mol: Chem.Mol, index: int, template: AttachmentTemplate) -> None:
    if not 0 <= index < mol.GetNumAtoms():
        raise BuilderError(f"Anchor atom {index} does not exist; the ligand has "
                           f"{mol.GetNumAtoms()} heavy atoms (0-{mol.GetNumAtoms() - 1})")
    if template.key == "Au_thiolate" and not is_thiol_sulfur(mol, index):
        atom = mol.GetAtomWithIdx(index)
        raise BuilderError(f"No valid Au-S anchor: the Au-thiolate template needs a thiol sulfur "
                           f"(-SH); atom {index} is {atom.GetSymbol()} with "
                           f"{atom.GetTotalNumHs()} H")
    if template.substrate == "carbon":
        atom = mol.GetAtomWithIdx(index)
        if template.anchor_elements and atom.GetSymbol() not in template.anchor_elements:
            raise BuilderError(f"{template.label} needs an anchor atom of "
                               f"{', '.join(template.anchor_elements)}; atom {index} is "
                               f"{atom.GetSymbol()}")
        if template.removes_anchor_h and atom.GetTotalNumHs() < 1 and \
                atom.GetNumRadicalElectrons() < 1:
            raise BuilderError(f"Anchor atom {index} ({atom.GetSymbol()}) has no hydrogen to give up "
                               "for the new bond: pick an atom with H (e.g. the N of an amine, the "
                               "O of an alcohol, a C-H carbon) or write a radical such as [CH2]")


def _embed(molh: Chem.Mol, n_conf: int, seed: int):
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    ids = list(AllChem.EmbedMultipleConfs(molh, numConfs=int(n_conf), params=params))
    if not ids:
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(molh, numConfs=int(n_conf), params=params))
    if not ids:
        raise BuilderError("RDKit could not generate 3D coordinates for this ligand")
    if AllChem.MMFFHasAllMoleculeParams(molh):
        AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=2000)
        method = "MMFF94 (isolated ligand)"
    elif AllChem.UFFHasAllMoleculeParams(molh):
        AllChem.UFFOptimizeMoleculeConfs(molh, maxIters=2000)
        method = "UFF (isolated ligand; MMFF parameters unavailable)"
    else:
        method = "none (no MMFF/UFF parameters; ETKDG geometry only)"
    return ids, method


def _heavy_bfs_tree(molh: Chem.Mol, root: int):
    """Parent map and BFS order over heavy atoms, starting at ``root``."""
    parent, order, queue = {root: None}, [root], [root]
    while queue:
        nxt = []
        for a in queue:
            for n in molh.GetAtomWithIdx(a).GetNeighbors():
                j = n.GetIdx()
                if n.GetAtomicNum() > 1 and j not in parent:
                    parent[j] = a
                    order.append(j)
                    nxt.append(j)
        queue = nxt
    return parent, order


def extend_chains(molh: Chem.Mol, conf_id: int, anchor: int) -> int:
    """Set every acyclic single-bond backbone dihedral to 180 deg (all-trans).

    Walking outward from the anchor, each heavy atom ``c`` with parent ``b``
    and grandparent ``a`` gets the dihedral a-b-c-d = 180 deg, where ``d`` is
    the child of ``c`` leading to its deepest branch. RDKit rotates the part
    of the molecule beyond the bond, so the anchor side never moves. Double
    bonds (e.g. the cis bond of oleyl chains) and rings are left untouched.
    Returns the number of dihedrals set.
    """
    parent, order = _heavy_bfs_tree(molh, anchor)
    children = {a: [] for a in order}
    for a in order:
        if parent[a] is not None:
            children[parent[a]].append(a)
    depth = {}
    for a in reversed(order):
        depth[a] = 1 + max((depth[c] for c in children[a]), default=0)
    conf = molh.GetConformer(conf_id)
    count = 0
    for c in order:
        b = parent[c]
        if b is None or parent[b] is None or not children[c]:
            continue
        a = parent[b]
        bond = molh.GetBondBetweenAtoms(b, c)
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        d = max(children[c], key=lambda x: depth[x])
        rdMolTransforms.SetDihedralDeg(conf, a, b, c, d, 180.0)
        count += 1
    return count


def build_ligand(ligand_config: LigandConfig) -> LigandModel:
    """Parse SMILES, choose anchor, perform the permitted edit, embed in 3D.

    The conformer is generated for the intact molecule (so MMFF has parameters
    for a thiol); the Au-thiolate edit then removes only the S-bound hydrogen.
    """
    ligand_config.validate()
    template = TEMPLATES[ligand_config.binding_template]
    mol = parse_smiles(ligand_config.smiles)
    anchor = int(ligand_config.anchor_atom_index)
    check_anchor(mol, anchor, template)

    molh = Chem.AddHs(mol)  # heavy-atom indices are preserved; H atoms are appended
    conf_ids, method = _embed(molh, ligand_config.conformers, ligand_config.seed)
    extended = ligand_config.conformer_mode == "extended"
    can_relax = AllChem.MMFFHasAllMoleculeParams(molh) or AllChem.UFFHasAllMoleculeParams(molh)
    extend_note = None
    if extended and not can_relax:
        # Setting dihedrals without a relaxation can leave branch atoms overlapping.
        extended = False
        extend_note = ("Extended chains skipped: neither MMFF nor UFF has parameters for this "
                       "ligand to relax the straightened geometry; the RDKit conformer is used")
    if extended:
        for cid in conf_ids:
            extend_chains(molh, cid, anchor)
        # Short local relaxation: removes strain from the rotations while staying in the
        # all-trans basin (a minimization does not cross torsion barriers).
        if AllChem.MMFFHasAllMoleculeParams(molh):
            AllChem.MMFFOptimizeMoleculeConfs(molh, maxIters=200)
        elif AllChem.UFFHasAllMoleculeParams(molh):
            AllChem.UFFOptimizeMoleculeConfs(molh, maxIters=200)
        method += "; backbone dihedrals set all-trans, then locally relaxed"

    removed = []
    if template.key == "Au_thiolate":
        h_idx = [n.GetIdx() for n in molh.GetAtomWithIdx(anchor).GetNeighbors()
                 if n.GetSymbol() == "H"]
        removed.append({"element": "H", "explicit_h_index": h_idx[0],
                        "reason": "thiol hydrogen removed by Au-thiolate template"})
    elif template.key == "chemisorption" and molh.GetAtomWithIdx(anchor).GetSymbol() in ("S", "O") \
            and any(n.GetSymbol() == "H" for n in molh.GetAtomWithIdx(anchor).GetNeighbors()):
        h_idx = [n.GetIdx() for n in molh.GetAtomWithIdx(anchor).GetNeighbors()
                 if n.GetSymbol() == "H"]
        removed.append({"element": "H", "explicit_h_index": h_idx[0],
                        "reason": "anchor hydrogen removed on chemisorption (S-H / O-H)"})
    elif template.removes_anchor_h and molh.GetAtomWithIdx(anchor).GetNumRadicalElectrons() < 1:
        h_idx = [n.GetIdx() for n in molh.GetAtomWithIdx(anchor).GetNeighbors()
                 if n.GetSymbol() == "H"]
        removed.append({"element": "H", "explicit_h_index": h_idx[0],
                        "reason": f"anchor hydrogen removed to form the covalent bond "
                                  f"({template.label})"})
    removed_idx = {r["explicit_h_index"] for r in removed}

    keep = [a.GetIdx() for a in molh.GetAtoms() if a.GetIdx() not in removed_idx]
    index_map = {old: new for new, old in enumerate(keep)}
    elements = [molh.GetAtomWithIdx(i).GetSymbol() for i in keep]
    formal = [molh.GetAtomWithIdx(i).GetFormalCharge() for i in keep]
    bonds = []
    for b in molh.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in removed_idx or j in removed_idx:
            continue
        bonds.append((index_map[i], index_map[j], _BOND_ORDERS.get(b.GetBondType(), 1.0)))
    new_anchor = index_map[anchor]

    conformers, axes = [], []
    covalent_axis = template.substrate == "carbon" and bool(removed)
    lone_pair_axis = template.key == "chemisorption"
    anchor_nb = [index_map[n.GetIdx()] for n in molh.GetAtomWithIdx(anchor).GetNeighbors()
                 if n.GetIdx() not in removed_idx]
    for cid in conf_ids:
        full = molh.GetConformer(cid).GetPositions()
        xyz = full[keep]
        if lone_pair_axis and anchor_nb:
            # The surface bond points along the anchor's free (lone-pair) direction, opposite
            # the sum of its bonds; averaged with the direction of the ligand body so long
            # chains still stand up from curved surfaces. Amine H atoms then point away
            # from the surface instead of into it.
            u = sum((xyz[k] - xyz[new_anchor]) / np.linalg.norm(xyz[k] - xyz[new_anchor])
                    for k in anchor_nb)
            heavy = [k for k in range(len(keep)) if k != new_anchor
                     and molh.GetAtomWithIdx(keep[k]).GetSymbol() != "H"]
            body = (xyz[heavy].mean(axis=0) - xyz[new_anchor]) if heavy else u
            u = u / max(np.linalg.norm(u), 1e-9) + body / max(np.linalg.norm(body), 1e-9)
            if np.linalg.norm(u) < 1e-6:
                u = body
            axes.append(u / np.linalg.norm(u))
        elif covalent_axis:
            # The new bond points where the removed hydrogen was, so the anchor keeps
            # its tetrahedral/trigonal geometry: outward axis = anchor - H.
            v = full[anchor] - full[removed[0]["explicit_h_index"]]
            axes.append(v / np.linalg.norm(v))
        conformers.append(xyz - xyz[new_anchor])
    positions = conformers[0]
    if not np.all(np.isfinite(positions)):
        raise BuilderError("Generated ligand coordinates are not finite")

    others = [k for k in range(len(elements)) if k != new_anchor]
    heavy_others = [k for k in others if elements[k] != "H"]
    ref = heavy_others or others
    axis = positions[ref].mean(axis=0) if ref else np.array([0.0, 0.0, 1.0])
    norm = np.linalg.norm(axis)
    axis = axis / norm if norm > 1e-6 else np.array([0.0, 0.0, 1.0])
    if covalent_axis or (lone_pair_axis and axes):
        axis = axes[0]

    warns = [f"Ligand geometry optimized with {method}; this is not a force field for the "
             "ligand-surface interface"]
    if extend_note:
        warns.append(extend_note)
    if template.key == "Au_thiolate":
        warns.append("Au-thiolate: the S-bound hydrogen was removed. Bound-state charges and "
                      "Au-S parameters are unassigned; gas-phase MMFF charges are not transferred")
    if Chem.GetFormalCharge(mol) != 0:
        warns.append(f"Ligand formal charge {Chem.GetFormalCharge(mol):+d} from SMILES; this does "
                     "not define the physical interfacial charge distribution")

    if template.substrate == "carbon":
        warns.append(f"{template.label}: covalent bond to the carbon substrate; the anchor "
                     + ("lost one H" if removed else "radical was used as is")
                     + ". Bonded parameters across the link must come from a force field")
    surface_species = list(ligand_config.surface_species or template.surface_species)
    default = template_distance(template.key, mol.GetAtomWithIdx(anchor).GetSymbol(),
                                surface_species[0] if surface_species else "")
    if template.key == "chemisorption":
        warns.append("Chemisorption: an anchor-surface bond is written to the topology"
                     + (" and the anchor H was removed" if removed else "")
                     + "; its force-field parameters are not assigned")
    distance = float(ligand_config.distance_A or default)
    if ligand_config.distance_A and abs(distance - default) > 0.3:
        warns.append(f"Anchor distance {distance} A differs from the template default "
                     f"{default} A")
    return LigandModel(
        smiles=Chem.MolToSmiles(mol), elements=elements, positions_A=positions, bonds=bonds,
        formal_charges=formal, anchor_index=new_anchor, index_map=index_map,
        removed_atoms=removed, outward_axis=axis, binding_template=template.key,
        distance_A=distance, surface_species=surface_species,
        formula=_formula(elements), orientation=ligand_config.orientation,
        conformers_A=conformers, optimization=method, warnings=warns, outward_axes=axes)


def _formula(elements: list) -> str:
    counts = {}
    for e in elements:
        counts[e] = counts.get(e, 0) + 1
    order = [e for e in ("C", "H") if e in counts] + sorted(e for e in counts
                                                             if e not in ("C", "H"))
    return "".join(f"{e}{counts[e] if counts[e] > 1 else ''}" for e in order)
