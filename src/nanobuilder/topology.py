"""Assemble the canonical topology shared by all exporters (spec stage E).

Atom order: substrate first, then ligand 1, ligand 2, ... Internal indices
are zero-based. Molecule IDs: 0 for the substrate, 1..n for ligands.

``type_labels`` are structural labels built from element and chemical role
(e.g. ``Au_sub``, ``S_anchor``, ``C``). They let MD engines read the files,
but they are NOT force-field types; ``atom_types`` stays None until a
force field assigns real types.
"""

from __future__ import annotations

import importlib.metadata
import platform
from datetime import datetime, timezone

import numpy as np
from ase.data import atomic_masses, atomic_numbers

from . import __version__
from .cell import finalize_box, make_molecules_whole
from .chemistry import TEMPLATES
from .schema import AssemblyModel, LigandModel, SubstrateModel

_ORDER_SYMBOL = {1.0: "-", 1.5: ":", 2.0: "=", 3.0: "#"}


def bond_type_name(a: str, b: str, order: float = 1.0) -> str:
    a, b = sorted((a, b))
    return f"{a}{_ORDER_SYMBOL.get(float(order), '-')}{b}"


def _chain_name(elements: list) -> str:
    forward, backward = "-".join(elements), "-".join(reversed(elements))
    return min(forward, backward)


def ligand_angles_dihedrals(ligand: LigandModel) -> tuple:
    """Intramolecular angle and proper-dihedral tuples from the bond graph."""
    n = ligand.n_atoms
    neighbors = [set() for _ in range(n)]
    for i, j, _ in ligand.bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)
    angles = []
    for j in range(n):
        nb = sorted(neighbors[j])
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                angles.append((nb[x], j, nb[y]))
    dihedrals = []
    for j, k, _ in ligand.bonds:
        for i in sorted(neighbors[j] - {k}):
            for m in sorted(neighbors[k] - {j}):
                if i != m:
                    dihedrals.append((i, j, k, m))
    return angles, dihedrals


def package_versions() -> dict:
    out = {"nanobuilder": __version__, "python": platform.python_version()}
    for name in ("ase", "rdkit", "numpy", "scipy", "gsd", "ipywidgets"):
        try:
            out[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            out[name] = None
    return out


def angles_dihedrals(n: int, bonds, elements) -> tuple:
    """All angles and proper dihedrals of a bond graph, with element-chain type names."""
    neighbors = [[] for _ in range(n)]
    for i, j, *_ in bonds:
        neighbors[i].append(j)
        neighbors[j].append(i)
    for nb in neighbors:
        nb.sort()
    angles = []
    for j in range(n):
        nb = neighbors[j]
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                i, k = nb[x], nb[y]
                angles.append((i, j, k, _chain_name([elements[i], elements[j], elements[k]])))
    dihedrals = []
    for j, k, *_ in bonds:
        for i in neighbors[j]:
            if i == k:
                continue
            for m in neighbors[k]:
                if m == j or m == i:
                    continue
                dihedrals.append((i, j, k, m, _chain_name([elements[i], elements[j], elements[k],
                                                           elements[m]])))
    return angles, dihedrals


def assemble(substrate: SubstrateModel, ligand: LigandModel | None, placements: list,
             target_ligands: int, surface_area_nm2: float, requested_density_nm2,
             warnings: list, seeds: dict) -> AssemblyModel:
    """Build the canonical model from the substrate and placed ligand copies.

    ``placements`` items: dicts with ``positions`` (ligand coordinates in the
    output frame), ``site_index`` and ``surface_atom_index``.
    """
    if substrate.bonds:
        return _assemble_covalent(substrate, ligand, placements, target_ligands, surface_area_nm2,
                                  requested_density_nm2, warnings, seeds)
    sub = substrate.atoms
    elements = list(sub.get_chemical_symbols())
    positions = [sub.get_positions()]
    roles = ["substrate"] * len(sub)
    component = ["substrate"] * len(sub)
    molecule_id = [0] * len(sub)
    bonds, angles, dihedrals, site_to_ligand = [], [], [], []

    template = TEMPLATES[ligand.binding_template] if ligand is not None else None
    if ligand is not None and placements:
        lig_angles, lig_dihedrals = ligand_angles_dihedrals(ligand)
        lig_roles = ["anchor" if k == ligand.anchor_index else "ligand"
                     for k in range(ligand.n_atoms)]
        for number, placement in enumerate(placements, start=1):
            offset = len(elements)
            elements.extend(ligand.elements)
            positions.append(np.asarray(placement["positions"]))
            roles.extend(lig_roles)
            component.extend(["ligand"] * ligand.n_atoms)
            molecule_id.extend([number] * ligand.n_atoms)
            for i, j, order in ligand.bonds:
                bonds.append((offset + i, offset + j,
                              bond_type_name(ligand.elements[i], ligand.elements[j], order)))
            for tup in lig_angles:
                angles.append(tuple(offset + t for t in tup) +
                              (_chain_name([ligand.elements[t] for t in tup]),))
            for tup in lig_dihedrals:
                dihedrals.append(tuple(offset + t for t in tup) +
                                 (_chain_name([ligand.elements[t] for t in tup]),))
            anchor = offset + ligand.anchor_index
            surface_atom = int(placement["surface_atom_index"])
            if template.covalent_surface_bond:
                bonds.append((surface_atom, anchor, template.bond_type or bond_type_name(
                    elements[surface_atom], elements[anchor], 1.0)))
            site_to_ligand.append({
                "ligand_id": number,
                "site_index": int(placement["site_index"]),
                "surface_atom_index": surface_atom,
                "anchor_atom_index": anchor,
                "face": placement.get("face"),
                "covalent_bond": bool(template.covalent_surface_bond),
            })

    all_positions = np.concatenate(positions)
    all_positions, cell = finalize_box(all_positions, np.array(sub.cell), substrate.pbc,
                                       substrate.vacuum_A, substrate.box_shape,
                                       substrate.box_min_A, warnings)
    # finalize_box wraps atom by atom, which can cut a ligand at the cell edge
    # (bonds then span the whole slab). Rebuild each ligand as one piece hanging
    # from its surface atom; the attachment pair counts as a bond for this even
    # when the template makes no covalent surface bond.
    links = [(s["surface_atom_index"], s["anchor_atom_index"]) for s in site_to_ligand]
    all_positions = make_molecules_whole(
        all_positions, cell, substrate.pbc, list(bonds) + links, molecule_id,
        roots={s["ligand_id"]: s["surface_atom_index"] for s in site_to_ligand})
    labels = []
    for element, role in zip(elements, roles):
        labels.append({"substrate": f"{element}_sub", "anchor": f"{element}_anchor"}
                      .get(role, element))
    masses = np.array([atomic_masses[atomic_numbers[e]] for e in elements])

    provenance = _provenance(substrate, ligand, seeds)
    return AssemblyModel(
        positions_A=all_positions, elements=elements, cell_A=cell, pbc=tuple(substrate.pbc),
        atom_types=None, type_labels=labels, chemical_role=roles, masses_amu=masses,
        charges_e=None, molecule_id=np.array(molecule_id, dtype=int), component=component,
        bond_list=bonds, angle_list=angles, dihedral_list=dihedrals, improper_list=[],
        site_to_ligand=site_to_ligand, surface_area_nm2=surface_area_nm2,
        target_ligands=int(target_ligands), achieved_ligands=len(placements),
        requested_density_nm2=requested_density_nm2, warnings=list(warnings),
        provenance=provenance, parameter_status="unassigned",
        bond_orders=_orders_from_names(bonds),
        substrate_index=list(range(len(sub))) + [-1] * (len(elements) - len(sub)))


def _orders_from_names(bonds) -> list:
    inverse = {v: k for k, v in _ORDER_SYMBOL.items()}
    out = []
    for _, _, name in bonds:
        sym = next((c for c in name if c in inverse), "-")
        out.append(inverse.get(sym, 1.0))
    return out


def _provenance(substrate: SubstrateModel, ligand, seeds) -> dict:
    crystal = substrate.crystal
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "versions": package_versions(),
        "material": crystal.provenance(),
        "geometry": {
            "shape": substrate.shape,
            "pbc_intended": list(substrate.pbc),
            "miller": list(substrate.miller) if substrate.miller else None,
            "crystal_frame_normal": (substrate.crystal_frame_normal.round(6).tolist()
                                     if substrate.crystal_frame_normal is not None else None),
            "layers": substrate.layers,
            "termination": (substrate.terminations[substrate.termination_index].label
                            if substrate.terminations else None),
            "radius_A": substrate.radius_A,
            "rod": None if substrate.shape != "rod" else {
                "axis_uvw": list(substrate.rod_crystal_axis),
                "caps": substrate.rod_caps,
                "cylinder_half_length_A": substrate.rod_half_length_A,
            },
            "vacuum_A": substrate.vacuum_A,
            "composition": substrate.composition,
            "channel_type": (substrate.channel_info or {}).get("type"),
            "channel": dict(substrate.channel_info) if substrate.channel_info else None,
        },
        "ligand": None if ligand is None else {
            "smiles": ligand.smiles,
            "formula_bound": ligand.formula,
            "binding_template": ligand.binding_template,
            "distance_A": ligand.distance_A,
            "surface_species": ligand.surface_species,
            "removed_atoms": ligand.removed_atoms,
            "optimization": ligand.optimization,
        },
        "seeds": seeds,
        "units": {"length": "Angstrom", "mass": "amu", "charge": "e (unassigned)"},
    }
    if substrate.carbon_info:
        provenance["carbon"] = substrate.carbon_info
    return provenance


def _assemble_covalent(substrate: SubstrateModel, ligand, placements, target_ligands,
                       surface_area_nm2, requested_density_nm2, warnings, seeds) -> AssemblyModel:
    """Canonical model for covalent (carbon) substrates.

    Applies each placement's substrate edits (atoms removed, added or moved,
    relabels, broken bonds), keeps every substrate residue contiguous, adds the
    ligand copies and their links, and enumerates angles and proper dihedrals
    over the complete bond graph (including across the new surface bonds).
    """
    from .carbon import lattice_bond_order

    sub = substrate.atoms
    pos = [p for p in sub.get_positions()]
    el = list(sub.get_chemical_symbols())
    lab = list(substrate.labels)
    res = [list(r) for r in substrate.residues]
    bond_order = {(min(i, j), max(i, j)): float(o) for i, j, o in substrate.bonds}
    removed = set()
    links = []
    for p in placements:
        e = p.get("edits") or {}
        for idx, xyz in e.get("move", {}).items():
            pos[idx] = np.asarray(xyz, dtype=float)
        for idx, label in e.get("relabel", {}).items():
            lab[idx] = label
        for i, j in e.get("unbond", ()):
            bond_order.pop((min(i, j), max(i, j)), None)
        for idx, name in e.get("residue", {}).items():
            res[idx][0] = name
        removed |= set(e.get("remove", ()))
        for a in e.get("add", ()):
            k = len(el)
            el.append(a["element"])
            pos.append(np.asarray(a["position"], dtype=float))
            lab.append(a["label"])
            res.append(list(res[a["residue_of"]]))
            bond_order[(min(k, a["bond_to"]), max(k, a["bond_to"]))] = float(a.get("order", 1.0))
    # Substrate order: residues contiguous (added atoms join their residue), removed atoms gone.
    keep = [k for k in range(len(el)) if k not in removed]
    keep.sort(key=lambda k: (res[k][1], k))
    index = {old: new for new, old in enumerate(keep)}
    elements = [el[k] for k in keep]
    labels = [lab[k] for k in keep]
    positions = [pos[k] for k in keep]
    resnames, resids, last, counter = [], [], None, 0
    for k in keep:
        if res[k][1] != last:
            counter += 1
            last = res[k][1]
        resnames.append(res[k][0])
        resids.append(counter)
    n_sub = len(keep)
    pairs = []
    for (i, j), o in bond_order.items():
        if i in removed or j in removed:
            continue
        a, b = sorted((index[i], index[j]))
        if elements[a] == "C" and elements[b] == "C":
            o = lattice_bond_order(labels[a], labels[b], o)
        pairs.append((a, b, o))
    pairs.sort()
    roles = ["substrate"] * n_sub
    component = ["substrate"] * n_sub
    molecule_id = [0] * n_sub
    site_to_ligand = []

    template = TEMPLATES[ligand.binding_template] if ligand is not None else None
    if ligand is not None and placements:
        lig_roles = ["anchor" if k == ligand.anchor_index else "ligand"
                     for k in range(ligand.n_atoms)]
        for number, placement in enumerate(placements, start=1):
            offset = len(elements)
            elements.extend(ligand.elements)
            positions.extend(np.asarray(placement["positions"]))
            labels.extend(["%s_anchor" % e if r == "anchor" else e
                           for e, r in zip(ligand.elements, lig_roles)])
            roles.extend(lig_roles)
            component.extend(["ligand"] * ligand.n_atoms)
            molecule_id.extend([number] * ligand.n_atoms)
            resnames.extend(["LIG"] * ligand.n_atoms)
            resids.extend([resids[-1] + 1 if resids else 1] * ligand.n_atoms)
            for i, j, order in ligand.bonds:
                pairs.append((offset + i, offset + j, float(order)))
            anchor = offset + ligand.anchor_index
            surface_atom = index[int(placement["surface_atom_index"])]
            if template.covalent_surface_bond:
                pairs.append((surface_atom, anchor, 1.0))
            links.append((surface_atom, anchor))
            site_to_ligand.append({
                "ligand_id": number,
                "site_index": int(placement["site_index"]),
                "surface_atom_index": surface_atom,
                "anchor_atom_index": anchor,
                "face": placement.get("face"),
                "covalent_bond": bool(template.covalent_surface_bond),
                "template": template.key,
            })
    bonds = [(i, j, bond_type_name(elements[i], elements[j], o)) for i, j, o in pairs]
    orders = [o for _, _, o in pairs]
    angles, dihedrals = angles_dihedrals(len(elements), pairs, elements)

    all_positions = np.array(positions, dtype=float)
    all_positions, cell = finalize_box(all_positions, np.array(sub.cell), substrate.pbc,
                                       substrate.vacuum_A, substrate.box_shape,
                                       substrate.box_min_A, warnings)
    all_positions = make_molecules_whole(
        all_positions, cell, substrate.pbc, [b[:2] for b in pairs] + links, molecule_id,
        roots={s["ligand_id"]: s["surface_atom_index"] for s in site_to_ligand})
    masses = np.array([atomic_masses[atomic_numbers[e]] for e in elements])
    provenance = _provenance(substrate, ligand, seeds)
    from .carbon import composition_metrics
    sub_metrics = composition_metrics(elements[:n_sub], labels[:n_sub])
    provenance["carbon"] = dict(provenance.get("carbon", {}))
    provenance["carbon"]["composition_after_attachment"] = sub_metrics

    names = _atom_names(elements, resids)
    assembly = AssemblyModel(
        positions_A=all_positions, elements=elements, cell_A=cell, pbc=tuple(substrate.pbc),
        atom_types=None, type_labels=labels, chemical_role=roles, masses_amu=masses,
        charges_e=None, molecule_id=np.array(molecule_id, dtype=int), component=component,
        bond_list=bonds, angle_list=angles, dihedral_list=dihedrals, improper_list=[],
        site_to_ligand=site_to_ligand, surface_area_nm2=surface_area_nm2,
        target_ligands=int(target_ligands), achieved_ligands=len(placements),
        requested_density_nm2=requested_density_nm2, warnings=list(warnings),
        provenance=provenance, parameter_status="unassigned", bond_orders=orders,
        residue_names=resnames, residue_ids=resids, atom_names=names,
        substrate_index=list(keep) + [-1] * (len(elements) - n_sub))
    return assembly


def _atom_names(elements, resids) -> list:
    """Element + running number within each residue (C1, O1, H1 ...)."""
    out, counters, last = [], {}, None
    for e, r in zip(elements, resids):
        if r != last:
            counters, last = {}, r
        counters[e] = counters.get(e, 0) + 1
        out.append(f"{e.upper()}{counters[e]}")
    return out


def parameter_coverage(assembly: AssemblyModel) -> dict:
    """Every interaction family and the types that still lack parameters."""
    def unique(items):
        return sorted(set(items))

    return {
        "force_field": None,
        "atom_types_assigned": assembly.atom_types is not None,
        "structural_type_labels": unique(assembly.type_labels),
        "charges_assigned": assembly.charges_e is not None,
        "missing": {
            "nonbonded_pairs": unique(assembly.type_labels),
            "bond_types": unique(b[2] for b in assembly.bond_list),
            "angle_types": unique(a[3] for a in assembly.angle_list),
            "dihedral_types": unique(d[4] for d in assembly.dihedral_list),
            "improper_types": [],
            "surface_ligand_cross_interactions": True,
            "substrate_model": "metallic/covalent substrate interactions not defined",
        },
        "notes": ("Impropers are not enumerated (their definition is force-field specific; "
                  "sp2 carbon frameworks usually need them for planarity). No parameters were "
                  "generated; supply a compatible force field before running MD."),
    }
