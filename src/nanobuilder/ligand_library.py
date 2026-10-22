"""Curated ligand library: common nanoparticle capping agents as ready-made inputs.

Each entry is only a *starting point* for the ligand tab: SMILES, anchor atom
(heavy-atom index in that SMILES), preferred attachment template, a typical
anchor-to-surface distance and the conformer mode. Every field can still be
edited in the GUI, and the normal SMILES -> RDKit -> 3D pipeline builds the
molecule, so library ligands are handled exactly like typed ones.

Distances are typical starting values for a geometric placement (roughly the
metal-anchor bond length on noble metals), not fitted parameters. Thiols use
the reviewed Au-thiolate template where it is available (Au/FCC preset) and
fall back to geometric adsorption of the intact thiol elsewhere. The
"Covalent on carbon" entries use the carbon templates (bond lengths come from
``chemistry.COVALENT_LENGTHS``); on metals they fall back to geometric
adsorption.
"""

from __future__ import annotations

from .schema import BuilderError

# group, label, smiles, anchor index, preferred template, distance (A), conformer mode, note
_ENTRIES = [
    # --- thiols (S anchor) -------------------------------------------------
    ("mercaptoethanol", "Thiols", "2-Mercaptoethanol", "SCCO", 0, "Au_thiolate", None, "etkdg",
     "Short hydroxyl-terminated thiol (the workflow's default ligand)."),
    ("hexanethiol", "Thiols", "1-Hexanethiol (C6)", "SCCCCCC", 0, "Au_thiolate", None, "extended",
     "Short alkanethiol."),
    ("dodecanethiol", "Thiols", "1-Dodecanethiol (C12)", "SCCCCCCCCCCCC", 0, "Au_thiolate", None,
     "extended", "Classic alkanethiol SAM ligand (Brust-Schiffrin particles)."),
    ("mud", "Thiols", "11-Mercapto-1-undecanol (MUD)", "SCCCCCCCCCCCO", 0, "Au_thiolate", None,
     "extended", "Hydroxyl-terminated C11 thiol."),
    ("mua", "Thiols", "11-Mercaptoundecanoic acid (MUA)", "SCCCCCCCCCCC(=O)O", 0, "Au_thiolate",
     None, "extended", "Carboxyl-terminated C11 thiol (neutral acid form)."),
    ("mpa", "Thiols", "3-Mercaptopropionic acid (MPA)", "SCCC(=O)O", 0, "Au_thiolate", None,
     "extended", "Short carboxyl-terminated thiol (neutral acid form)."),
    ("cysteamine", "Thiols", "Cysteamine", "SCCN", 0, "Au_thiolate", None, "extended",
     "Amine-terminated thiol (neutral amine form)."),
    ("thiophenol", "Thiols", "Thiophenol", "Sc1ccccc1", 0, "Au_thiolate", None, "etkdg",
     "Aromatic thiol."),
    # --- PEG thiols ----------------------------------------------------------
    ("peg3_oh", "PEG thiols", "HS-C2-PEG3-OH", "SCCOCCOCCOCCO", 0, "Au_thiolate", None, "extended",
     "Thiol-terminated tri(ethylene glycol), hydroxyl end."),
    ("peg6_ome", "PEG thiols", "HS-C2-PEG6-OMe (mPEG)", "SCCOCCOCCOCCOCCOCCOC", 0, "Au_thiolate",
     None, "extended", "Methoxy-capped hexa(ethylene glycol) thiol, built all-trans."),
    ("peg12_ome", "PEG thiols", "HS-C2-PEG12-OMe (mPEG)",
     "SCC" + "OCC" * 12 + "OC", 0, "Au_thiolate", None, "extended",
     "Longer mPEG thiol; all-trans is a deliberately extended (brush-like) start."),
    # --- amines (N anchor) ---------------------------------------------------
    ("dodecylamine", "Amines", "Dodecylamine (C12)", "NCCCCCCCCCCCC", 0, "geometric", 2.4,
     "extended", "Primary alkylamine; N lone pair toward the surface."),
    ("oleylamine", "Amines", "Oleylamine", "NCCCCCCCC/C=C\\CCCCCCCC", 0, "geometric", 2.4,
     "extended", "Common capping agent in colloidal synthesis (cis C9=C10 double bond kept)."),
    # --- carboxylates (O anchor) ---------------------------------------------
    ("citrate", "Carboxylates", "Citrate (3-)", "[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]", 0,
     "geometric", 2.2, "etkdg",
     "Turkevich-type Au/Ag capping; placed through ONE carboxylate oxygen (real binding is "
     "often multidentate)."),
    ("oleate", "Carboxylates", "Oleate (1-)", "[O-]C(=O)CCCCCCC/C=C\\CCCCCCCC", 0, "geometric", 2.2,
     "extended", "Oleic acid as carboxylate, single-O placement."),
    ("acetate", "Carboxylates", "Acetate (1-)", "[O-]C(C)=O", 0, "geometric", 2.2, "etkdg",
     "Smallest carboxylate."),
    # --- phosphorus ligands --------------------------------------------------
    ("tpp", "Phosphines / P=O", "Triphenylphosphine (TPP)",
     "P(c1ccccc1)(c1ccccc1)c1ccccc1", 0, "geometric", 2.3, "etkdg",
     "P anchor (Au55(PPh3)12Cl6-type clusters)."),
    ("topo", "Phosphines / P=O", "Trioctylphosphine oxide (TOPO)",
     "O=P(CCCCCCCC)(CCCCCCCC)CCCCCCCC", 0, "geometric", 2.2, "extended",
     "Binds through the P=O oxygen (quantum-dot synthesis)."),
    # --- covalent functionalization of carbon (graphene, GO, CNT, fullerene) ----
    ("phenyl", "Covalent on carbon", "Phenyl (aryl radical, diazonium route)", "c1ccccc1", 0,
     "C_basal", None, "etkdg", "Aryl group bonded to a basal carbon (sp2 to sp3)."),
    ("nitrophenyl", "Covalent on carbon", "4-Nitrophenyl (from 4-nitrobenzenediazonium)",
     "c1ccc(cc1)[N+](=O)[O-]", 0, "C_basal", None, "etkdg",
     "Classic diazonium grafting of graphene and nanotubes; anchor is the ring carbon para to NO2."),
    ("carboxyphenyl", "Covalent on carbon", "4-Carboxyphenyl", "c1ccc(cc1)C(=O)O", 0, "C_basal",
     None, "etkdg", "Aryl carboxylic acid for further coupling; anchor para to COOH."),
    ("aminophenyl", "Covalent on carbon", "4-Aminophenyl", "c1ccc(cc1)N", 0, "C_basal", None,
     "etkdg", "Aryl amine; anchor para to NH2."),
    ("methyl", "Covalent on carbon", "Methyl radical (methylated graphene)", "C", 0, "C_basal",
     None, "etkdg", "Smallest alkyl addition."),
    ("ethylenediamine", "Covalent on carbon", "Ethylenediamine (amide to GO COOH)", "NCCN", 0,
     "C_carboxyl", None, "extended",
     "EDC/NHS-type amide coupling to carboxylic acids at GO edges; free amine end."),
    ("octadecylamine", "Covalent on carbon", "Octadecylamine, ODA (amide to GO COOH)",
     "N" + "C" * 18, 0, "C_carboxyl", None, "extended",
     "Long alkylamine used to make GO dispersible in organic solvents."),
    ("mpeg_amine", "Covalent on carbon", "mPEG4-amine (amide to GO COOH)", "NCCOCCOCCOCCOC", 0,
     "C_carboxyl", None, "extended", "PEGylated GO through amide bonds at carboxylic acids."),
    ("ethanolamine", "Covalent on carbon", "Ethanolamine (epoxide ring-opening)", "NCCO", 0,
     "C_epoxide", None, "extended", "Amine opens a GO epoxide; the epoxide O stays as OH."),
    ("butylamine", "Covalent on carbon", "n-Butylamine (epoxide opening / C60 addition)", "NCCCC",
     0, "C_epoxide", None, "extended",
     "Primary amine; on pristine carbon (e.g. C60) use basal-plane addition."),
]

# Preferred carbon template for ligands whose main template is for metals (used only
# when the substrate is a carbon nanomaterial).
CARBON_FALLBACK = {"dodecylamine": "C_carboxyl", "oleylamine": "C_carboxyl",
                   "cysteamine": "C_carboxyl"}

LIBRARY = {
    key: {"group": group, "label": label, "smiles": smiles, "anchor_atom_index": anchor,
          "binding_template": template, "distance_A": distance, "conformer_mode": mode,
          "note": note, "carbon_template": CARBON_FALLBACK.get(key)}
    for key, group, label, smiles, anchor, template, distance, mode, note in _ENTRIES
}

GROUPS = list(dict.fromkeys(e["group"] for e in LIBRARY.values()))


def suits(key: str, carbon=None) -> bool:
    """Whether a library ligand makes sense for the substrate kind.

    Carbon substrates take the covalent carbon entries plus ligands with a
    carbon template (amines via amide coupling); metals and other crystals take
    everything except the covalent carbon chemistry. ``carbon=None``: any.
    """
    if carbon is None:
        return True
    e = LIBRARY[key]
    covalent = e["group"] == "Covalent on carbon"
    return (covalent or bool(e["carbon_template"])) if carbon else not covalent


def library_choices(carbon=None) -> list:
    """(label, key) pairs for a dropdown, grouped and labelled 'Group - Name'.

    ``carbon=True`` / ``False`` keeps only the ligands that suit a carbon
    nanomaterial / a crystal (see :func:`suits`).
    """
    out = []
    for group in GROUPS:
        out += [(f"{group} · {e['label']}", k) for k, e in LIBRARY.items()
                if e["group"] == group and suits(k, carbon)]
    return out


def library_entry(key: str) -> dict:
    if key not in LIBRARY:
        raise BuilderError(f"Unknown library ligand {key!r}")
    return dict(LIBRARY[key])
