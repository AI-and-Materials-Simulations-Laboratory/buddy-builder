"""Suggested force-field atom types (GAFF2-style names; OPLS-AA for carbon substrates).

These are *suggestions written into the MOL2 file and topology.json*; no
parameters or charges are assigned, so exports stay Level B.

- Carbon-substrate atoms are typed exactly from how they were built (their
  structural role), using the minimal GAFF mapping of the AMBER graphene-oxide
  dataset (ca, c3, oh, os, ho; Fedai et al., Data 11, 18, 2026) extended to
  edge groups (ha, c, o), and HierGO's OPLS-AA table (Oplsatom2).
- Ligand atoms get provisional GAFF2 names from simple connectivity rules
  (hybridization from neighbour count and bond orders, aromaticity, amide and
  carbonyl perception, electron-withdrawing neighbours for h1-h5). They cover
  common capping ligands but are not a substitute for antechamber; verify them
  (e.g. ``antechamber -i ligand.mol2 -fi mol2 -o out.mol2 -fo mol2 -at gaff2``).
- Metal atoms keep their element symbol (GAFF has no metal types).
"""

from __future__ import annotations

from .carbon import ROLE_TYPES

EW = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}
HALOGEN = {"F": "f", "Cl": "cl", "Br": "br", "I": "i"}


def _neighbors(n, bonds, orders):
    nb = [[] for _ in range(n)]
    for (i, j, *_), o in zip(bonds, orders):
        nb[i].append((j, float(o)))
        nb[j].append((i, float(o)))
    return nb


def _carbon_type(a, el, nb):
    links = nb[a]
    if any(o == 1.5 for _, o in links):
        return "ca"
    if any(o >= 2.5 for _, o in links) or (len(links) == 2 and sum(o for _, o in links) >= 4):
        return "c1"
    if any(o == 2.0 and el[j] in ("O", "S") for j, o in links):
        return "c"
    if any(o == 2.0 for _, o in links):
        return "c2"
    return "c3"


def _is_carbonyl_carbon(j, el, nb):
    return el[j] == "C" and any(o == 2.0 and el[k] in ("O", "S") for k, o in nb[j])


def _ligand_type(a, el, nb, charge, types):
    e = el[a]
    links = nb[a]
    heavy = [(j, o) for j, o in links if el[j] != "H"]
    if e == "C":
        return _carbon_type(a, el, nb)
    if e == "H":
        if not links:
            return "hc"
        j = links[0][0]
        x = el[j]
        if x == "N":
            return "hn"
        if x == "O":
            return "ho"
        if x == "S":
            return "hs"
        if x == "P":
            return "hp"
        if x == "C":
            ew = sum(1 for k, _ in nb[j] if k != a and el[k] in EW)
            ct = types[j] or _carbon_type(j, el, nb)
            if ct == "c3":
                return ("hc", "h1", "h2", "h3")[min(ew, 3)]
            return ("ha", "h4", "h5")[min(ew, 2)]
        return "hc"
    if e == "O":
        if len(links) == 1:
            j, o = links[0]
            if o == 2.0 or charge < 0:
                return "o"
            return "oh" if el[j] == "H" else "o"
        if any(el[j] == "H" for j, _ in links):
            return "oh"
        return "os"
    if e == "N":
        if any(el[j] == "O" for j, _ in heavy) and sum(1 for j, _ in heavy if el[j] == "O") >= 2:
            return "no"
        if len(links) == 4 or charge > 0 and len(links) >= 4:
            return "n4"
        if any(o == 1.5 for _, o in links):
            return "na" if len(links) == 3 else "nb"
        if any(o >= 2.5 for _, o in links) or len(links) == 1:
            return "n1"
        if any(o == 2.0 for _, o in links):
            return "n2"
        if len(links) == 3:
            if any(_is_carbonyl_carbon(j, el, nb) for j, _ in heavy):
                return "n"
            if any(types[j] == "ca" or (types[j] is None and el[j] == "C" and
                                        any(o == 1.5 for _, o in nb[j])) for j, _ in heavy):
                return "nh"
            return "n3"
        return "n3"
    if e == "S":
        if len(links) == 1 and links[0][1] == 2.0:
            return "s"
        if len(links) <= 2:
            return "sh" if any(el[j] == "H" for j, _ in links) else "ss"
        return "s4" if len(links) == 3 else "s6"
    if e == "P":
        return "p3" if len(links) <= 3 else "p5"
    if e in HALOGEN:
        return HALOGEN[e]
    return e  # metals and anything else: element symbol


def suggest_types(assembly, formal_charges=None) -> dict:
    """{'gaff2': [...], 'opls_hiergo': [...], 'notes': [...]} for every atom of ``assembly``."""
    n = assembly.n_atoms
    el = assembly.elements
    orders = assembly.bond_orders or [1.0] * len(assembly.bond_list)
    nb = _neighbors(n, assembly.bond_list, orders)
    charges = formal_charges or [0] * n
    gaff = [None] * n
    opls = [None] * n
    for a in range(n):
        if assembly.component[a] == "substrate":
            role = assembly.type_labels[a]
            if role in ROLE_TYPES:
                gaff[a], opls[a] = ROLE_TYPES[role]
            else:
                gaff[a] = el[a]
    for a in range(n):  # carbons first so hydrogens can look at them
        if gaff[a] is None and el[a] == "C":
            gaff[a] = _carbon_type(a, el, nb)
    for a in range(n):
        if gaff[a] is None:
            gaff[a] = _ligand_type(a, el, nb, charges[a], gaff)
    notes = ["GAFF2 names for carbon-substrate atoms follow how they were built (AMBER GO "
             "scheme: ca, c3, oh, os, ho; edge groups ha, c, o)",
             "Ligand GAFF2 names are provisional rule-based guesses: verify with antechamber",
             "No force-field parameters or charges were assigned"]
    return {"gaff2": gaff, "opls_hiergo": opls, "notes": notes}


def formal_charges_for(assembly, ligand) -> list:
    """Per-atom formal charge: substrate 0, ligand atoms from the RDKit molecule."""
    out = [0] * assembly.n_atoms
    if ligand is None:
        return out
    n_lig = ligand.n_atoms
    start = assembly.n_substrate
    for a in range(start, assembly.n_atoms):
        out[a] = int(ligand.formal_charges[(a - start) % n_lig])
    return out
