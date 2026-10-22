"""Hand edits after Generate: remove ligands, add a ligand at a chosen surface atom.

Edits change the list of ligand placements and rebuild the model with the
same assembly code as Generate, so bonds, angles, covalent carbon edits and
validation stay consistent. Every edit is recorded in
``provenance["manual_edits"]`` and flagged in README.txt of the export.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .atomtypes import formal_charges_for, suggest_types
from .chemistry import TEMPLATES
from .placement import N_AZIMUTHS, _bonded_neighbours, _ClashChecker, orient
from .schema import BuilderError
from .topology import assemble
from .validation import validate_assembly


def _rebuild(result, placements: list, edit: dict, closest=None):
    a0 = result.assembly
    warnings = [w for w in a0.warnings if not w.startswith("Manual edit")]
    warnings = [w for w in warnings if not w.startswith(("Only ", "Requested "))]
    # After a hand edit the ligands present are the ones the user wants.
    assembly = assemble(result.substrate, result.ligand, placements, len(placements),
                        a0.surface_area_nm2, a0.requested_density_nm2, warnings,
                        seeds=dict(a0.provenance.get("seeds", {})))
    assembly.placements = placements
    for key in ("placement", "seeds"):
        if key in a0.provenance:
            assembly.provenance[key] = a0.provenance[key]
    edits = list(a0.provenance.get("manual_edits", [])) + [edit]
    assembly.provenance["manual_edits"] = edits
    assembly.warnings.append(f"Manual edits after Generate: {len(edits)} (see README.txt); the "
                             "model is no longer the seeded build alone")
    best = a0.min_nonbonded_distance_A
    if closest is not None and np.isfinite(closest):
        best = closest if best is None else min(best, closest)
    assembly.min_nonbonded_distance_A = best
    assembly.suggested_types = suggest_types(assembly, formal_charges_for(assembly, result.ligand))
    report = validate_assembly(assembly)
    return replace(result, assembly=assembly, report=report)


def ligand_of_atom(result, atom_index: int) -> int:
    """Ligand number (1..N) that atom ``atom_index`` belongs to, or 0 for the substrate."""
    return int(result.assembly.molecule_id[int(atom_index)])


def remove_ligands(result, ligand_ids) -> object:
    """Remove ligands by number (1..N, the molecule ID shown in the viewer)."""
    ids = sorted({int(i) for i in ligand_ids})
    n = len(result.assembly.placements)
    bad = [i for i in ids if not 1 <= i <= n]
    if bad or not ids:
        raise BuilderError(f"No ligand {bad or ids}: ligands are numbered 1 to {n}")
    placements = [p for k, p in enumerate(result.assembly.placements, start=1) if k not in ids]
    return _rebuild(result, placements, {"action": "remove_ligands", "ligands": ids})


def add_ligand_at(result, atom_index: int) -> object:
    """Place one more ligand on the surface atom shown as ``atom_index`` in the model."""
    a0 = result.assembly
    if result.ligand is None or result.sites is None:
        raise BuilderError("Add ligands needs a model built with a ligand (set a ligand count > 0 "
                           "and Generate first)")
    atom_index = int(atom_index)
    if not 0 <= atom_index < a0.n_atoms:
        raise BuilderError(f"Atom {atom_index} does not exist")
    sub_idx = (a0.substrate_index or [])[atom_index] if a0.substrate_index else -1
    if sub_idx is None or sub_idx < 0:
        raise BuilderError("Pick a substrate (surface) atom to add a ligand")
    sites = result.sites
    cands = np.nonzero(np.asarray(sites.surface_atom_index) == sub_idx)[0]
    if not len(cands):
        raise BuilderError(f"Atom {atom_index} is not an open attachment site for this ligand in "
                           "the selected region (buried, wrong element, or wrong template)")
    used_sites = {int(p["surface_atom_index"]) for p in a0.placements}
    used_atoms = set()
    if sites.claims is not None:
        for p in a0.placements:
            used_atoms |= set(sites.claims[int(p["site_index"])])
    substrate, ligand = result.substrate, result.ligand
    checker = _ClashChecker(substrate, ligand, float(result.config.coverage.clash_alpha))
    for p in a0.placements:
        checker.commit(np.asarray(p["positions"]), np.inf)
    conformers = ligand.conformers_A or [ligand.positions_A]
    axes = ligand.outward_axes if len(ligand.outward_axes) == len(conformers) else \
        [ligand.outward_axis] * len(conformers)
    sub_nb = _bonded_neighbours(substrate) if substrate.bonds else {}
    lig_nb = {ligand.anchor_index}
    for i, j, _ in ligand.bonds:
        if i == ligand.anchor_index:
            lig_nb.add(j)
        elif j == ligand.anchor_index:
            lig_nb.add(i)
    for k in cands:
        k = int(k)
        claims = sites.claims[k] if sites.claims is not None else {sub_idx}
        if sub_idx in used_sites or (claims & used_atoms):
            raise BuilderError(f"Atom {atom_index} already carries a ligand")
        edits = sites.edits[k] if sites.edits is not None else {}
        exempt = dict(ignore=set(edits.get("ignore", ())), site_neighbors=sub_nb.get(sub_idx, set()),
                      anchor_shell=lig_nb) if sub_nb else (
            dict(anchor_shell=lig_nb) if TEMPLATES[ligand.binding_template].covalent_surface_bond
            else {})
        for conformer, axis in zip(conformers, axes):
            for step in range(N_AZIMUTHS * 2):
                xyz = orient(conformer, axis, sites.positions_A[k], sites.normals[k],
                             ligand.distance_A, 2 * np.pi * step / (N_AZIMUTHS * 2))
                ok, closest = checker.check(xyz, sub_idx, **exempt)
                if ok:
                    placement = {"positions": xyz, "site_index": k, "surface_atom_index": sub_idx,
                                 "face": sites.face[k], "edits": edits}
                    placements = list(a0.placements) + [placement]
                    return _rebuild(result, placements,
                                    {"action": "add_ligand", "surface_atom": atom_index,
                                     "ligand": len(placements)}, closest)
    raise BuilderError(f"A ligand at atom {atom_index} would clash with its neighbours at every "
                       "orientation tried")


__all__ = ["remove_ligands", "add_ligand_at", "ligand_of_atom"]
