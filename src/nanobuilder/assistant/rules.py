"""Offline parser: plain-English requests to console commands, without a model.

Used when no AI model is configured (for example on nanoHUB) or the model
cannot be reached. It understands the common phrasings, such as
"a 4 nm gold sphere with 60 dodecanethiol in fibonacci order, LAMMPS and PDB".
Anything it does not understand is listed so the user can rephrase or use the
tabs; it never guesses a material it cannot find.
"""

from __future__ import annotations

import re

from .. import material_catalog as catalog
from .. import mp_library
from ..ligand_library import LIBRARY

ELEMENT_NAMES = {
    "gold": "Au", "silver": "Ag", "platinum": "Pt", "palladium": "Pd", "copper": "Cu",
    "iron": "Fe", "nickel": "Ni", "aluminium": "Al", "aluminum": "Al", "cobalt": "Co",
    "titanium": "Ti", "magnesium": "Mg", "zinc": "Zn", "silicon": "Si", "germanium": "Ge",
    "lead": "Pb", "tungsten": "W", "molybdenum": "Mo", "chromium": "Cr", "rhodium": "Rh",
    "iridium": "Ir", "sodium chloride": "NaCl", "salt": "NaCl", "rock salt": "NaCl",
    "titania": "TiO2", "zinc oxide": "ZnO", "silica": "SiO2", "ceria": "CeO2",
    "cadmium selenide": "CdSe", "iron oxide": "Fe3O4",
}
NUM = r"(\d+(?:\.\d+)?)"
SHAPES = [("channel", ("channel", "slit", "pore", "pores")),
          ("rod", ("nanorod", "rod", "rods", "nanowire", "wire")),
          ("slab", ("slab", "surface", "film", "facet")),
          ("sphere", ("sphere", "spherical", "nanoparticle", "particle", "ball", "dot", "cluster",
                      "np", "nps", "nanosphere", "quantum dot"))]
FORMATS = ("lammps", "pdb", "mol2", "hoomd", "gsd", "xyz")
STOP = {"with", "and", "the", "a", "an", "of", "in", "on", "to", "nm", "make", "build", "create",
        "ligand", "ligands", "sphere", "slab", "rod", "box", "cube", "order", "coated", "coating",
        "group", "groups", "acid", "oxide", "amine", "thiol", "thiols", "carbon", "layer", "layers",
        "sheet", "flake", "periodic", "random", "even", "fibonacci", "nanorod", "surface"}


def _words(text):
    return re.findall(r"[a-z0-9]+", text.lower())


def _name_words() -> dict:
    """Distinctive words of catalog names ('anatase', 'magnetite') -> catalog key."""
    out = {}
    for name, key in catalog.NAMES.items():
        if catalog.CATALOG[key]["group"] == "Nanorods":
            continue
        for wd in re.split(r"[\s()\-,/]+", name.lower()):
            if len(wd) >= 5 and wd not in STOP and not wd.isdigit():
                out.setdefault(wd, key)
    return out


def find_material(text: str):
    """(console argument, requirement text) for the material in ``text``, or (None, None)."""
    low = text.lower()
    # carbon words first: 'graphene oxide' beats 'oxide'
    for word in sorted(catalog.CARBON_WORDS, key=len, reverse=True):
        if re.search(rf"(?<![a-z0-9]){re.escape(word)}(?![a-z0-9])", low):
            return word, word
    for name in sorted(ELEMENT_NAMES, key=len, reverse=True):
        if re.search(rf"\b{re.escape(name)}\b", low):
            return ELEMENT_NAMES[name], ELEMENT_NAMES[name]
    for tok in re.findall(r"\b([A-Z][a-z]?(?:[A-Z][a-z]?|\d)*)\b", text):  # Au, CdSe, TiO2, Fe3O4
        if len(tok) == 1 and tok in "AI":  # the words "A" and "I"
            continue
        if catalog.preset_for_text(tok) or mp_library.find(tok):
            return tok, tok
    words = _name_words()
    for wd in _words(text):
        if wd in words:
            return wd, wd
    m = re.search(r"\bmp-\d+\b", low)
    if m and m.group(0) in mp_library.by_id():
        return m.group(0), m.group(0)
    return None, None


def find_ligand(text: str):
    """Library key for the ligand named in ``text`` (key, label words or SMILES=...)."""
    low = text.lower()
    m = re.search(r"smiles\s*[:=]?\s*\"?([A-Za-z0-9@+\-\[\]()=#$/\\.%]+)", text, re.I)
    if m:
        return None, m.group(1)
    best = None
    for key, e in LIBRARY.items():
        names = {key, key.replace("_", " ")}
        label = e["label"].lower()
        names.add(re.sub(r"\s*\(.*?\)", "", label).strip())  # '1-dodecanethiol'
        names.update(re.findall(r"\(([a-z0-9]+)\)", label))  # '(mua)', '(c12)' -> abbreviations
        names.update(w for w in re.split(r"[\s\-,]+", label) if len(w) >= 6)
        for n in names:
            if len(n) >= 3 and not n.isdigit() and n not in ("c12", "c6", "c2") and \
                    re.search(rf"(?<![a-z0-9]){re.escape(n)}(?![a-z0-9])", low):
                if best is None or len(n) > best[1]:
                    best = (key, len(n))
    return (best[0], None) if best else (None, None)


SELECT_ELEMENTS = {"sulfur": "S", "sulphur": "S", "oxygen": "O", "carbon": "C", "hydrogen": "H",
                   "nitrogen": "N", "phosphorus": "P", "gold": "Au", "silver": "Ag", "platinum": "Pt",
                   "iron": "Fe", "selenium": "Se", "cadmium": "Cd", "zinc": "Zn", "titanium": "Ti"}


def _result(cmds, req=None, build=True):
    return {"commands": cmds, "requirements": req or {}, "unknown": [], "build": build,
            "question": False, "clarify": ""}


def model_edits(low: str, cfg=None, has_model: bool = False):
    """Chat edits of the current model or settings: 'remove ligand 3', 'highlight the sulfur
    atoms', '10 more ligands', 'make it bigger', 'undo'. None when the text is not one."""
    if re.fullmatch(r"\s*(undo|undo that|undo the last edit|go back)\s*[.!]?\s*", low):
        return _result(["undo()"])
    if re.fullmatch(r"\s*(summary|summarise|summarize|show (?:me )?(?:the )?summary)\s*[.!]?\s*", low):
        return _result(["summary()"])
    m = re.search(r"\bremove\s+ligands?\s+((?:\d+\s*(?:,|and)?\s*)+)", low)
    if m:
        ids = [int(x) for x in re.findall(r"\d+", m.group(1))]
        return _result([f"remove_ligand({', '.join(str(i) for i in ids)})"])
    m = re.search(r"\badd\s+(?:a\s+|one\s+)?ligand\s+(?:on|at|to)\s+atom\s+(\d+)", low)
    if m:
        return _result([f"add_ligand({m.group(1)})"])
    if re.search(r"\b(clear|reset)\s+(?:the\s+)?(selection|highlight)", low):
        return _result(["select()"])
    m = re.search(r"\b(?:highlight|select|show|colou?r)\s+(?:me\s+)?(?:the\s+|all\s+)*ligand\s+(\d+)", low)
    if m:
        return _result([f"select(ligand={m.group(1)})"])
    m = re.search(r"\b(?:highlight|select|show|colou?r)\s+(?:me\s+)?(?:the\s+|all\s+)*([a-z]+)\s+atoms?", low)
    if m and (m.group(1) in SELECT_ELEMENTS or m.group(1) in ("ligand", "surface", "substrate")):
        word = m.group(1)
        if word == "ligand":
            return _result(['select(role="ligand")'])
        if word in ("surface", "substrate"):
            return _result(['select(role="substrate")'])
        return _result([f'select(element="{SELECT_ELEMENTS[word]}")'])
    m = re.search(r"\bmeasure\s+(?:the\s+)?(?:distance\s+)?(?:between\s+)?(?:atoms?\s+)?(\d+)\s*(?:and|,|-|to)\s*(\d+)", low)
    if m:
        return _result([f"measure({m.group(1)}, {m.group(2)})"])
    m = re.fullmatch(r"\s*(?:add|put|attach|use|place)\s+(\d+)\s+(?:more\s+)?(.+?)\s*"
                     r"(?:ligands?|groups?|molecules?|chains?)?\s*[.!]?\s*", low)
    if m and not re.search(r"\b(more|fewer|less)\b", low):
        key, smiles = find_ligand(m.group(2))
        if key:
            n = int(m.group(1))
            return _result([f'ligand("{key}")', f"coverage(count={n})"], {"ligand": key, "ligand_count": n})
    if cfg is None:
        return None
    c, g = cfg.coverage, cfg.geometry
    m = re.search(r"\b(\d+)\s+(more|fewer|less)\s+(?:ligands?|molecules?|chains?)", low)
    if m:
        n = int(m.group(1)) * (1 if m.group(2) == "more" else -1)
        base = c.target if c.mode == "count" else 0
        new = max(0, base + n)
        return _result([f"coverage(count={new})"], {"ligand_count": new})
    if re.search(r"\b(double|twice as many)\b.*\bligands?\b", low) and c.mode == "count":
        return _result([f"coverage(count={c.target * 2})"], {"ligand_count": c.target * 2})
    if cfg.material.source != "carbon" and g.shape in ("sphere", "rod"):
        m = re.search(r"\b(?:make\s+it|diameter(?:\s+to)?|resize\s+to|change\s+(?:the\s+)?size\s+to)\s*"
                      + NUM + r"\s*nm", low)
        step = 1.0
        if m:
            d = float(m.group(1))
        elif re.search(r"\b(bigger|larger)\b", low):
            d = g.diameter_nm + step
        elif re.search(r"\bsmaller\b", low):
            d = max(1.0, g.diameter_nm - step)
        else:
            d = None
        if d is not None:
            return _result([f'shape("{g.shape}", diameter={d:g})'], {"shape": g.shape, "diameter_nm": d})
    return None


def parse(text: str, cfg=None, has_model: bool = False) -> dict:
    """{'commands', 'requirements', 'unknown', 'build', 'question'} for one request.

    ``cfg``: the current settings, for relative changes ("10 more ligands").
    """
    low = text.lower().strip()
    cmds, req, unknown = [], {}, []
    question = low.endswith("?") or bool(re.match(r"(why|what|how|which|can|does|is|explain)\b", low))
    build = bool(re.search(r"\b(build|make|create|generate|give me|i want|i need)\b", low))
    if not question:
        edit = model_edits(low, cfg, has_model)
        if edit is not None:
            return edit

    m = re.search(r"\bexample\s+([a-l])\b", low)
    if m:
        return {"commands": ["reset()", f'example("{m.group(1).upper()}")'],
                "requirements": {}, "unknown": [], "build": True, "question": False, "clarify": ""}

    arg, mat_req = find_material(text)
    carbon = arg is not None and arg.lower() in catalog.CARBON_WORDS
    new_material = arg is not None
    if new_material:
        cmds += ["reset()", f'material("{arg}")']
        req["material"] = mat_req
    # sizes that belong to something else (box, length, slab size) are not a diameter
    rest = re.sub(r"(?:cub(?:e|ic)\s*box|box)\s*(?:of\s*)?" + NUM + r"\s*nm|" + NUM +
                  r"\s*nm\s*(?:cub(?:e|ic)?\s*box|long|length)|" + NUM + r"\s*(?:x|by|×)\s*" + NUM +
                  r"\s*nm|length\s*(?:of\s*)?" + NUM + r"\s*nm", " ", low)
    nums_nm = [float(x) for x in re.findall(NUM + r"\s*-?\s*nm\b(?!\s*(?:\^?2|²))", rest)]

    if carbon:
        kw = []
        size = re.search(NUM + r"\s*(?:x|by|×)\s*" + NUM + r"\s*nm", low)
        sheet = arg.lower() in ("graphene", "go", "graphene oxide", "rgo")
        if size and sheet:
            kw.append(f"size=({size.group(1)}, {size.group(2)})")
            req["carbon_size_nm"] = [float(size.group(1)), float(size.group(2))]
        elif sheet and nums_nm:  # '3 nm graphene oxide sheet': a square sheet
            kw.append(f"size=({nums_nm[0]:g}, {nums_nm[0]:g})")
            req["carbon_size_nm"] = [nums_nm[0], nums_nm[0]]
        cnt = re.search(r"\(\s*(\d+)\s*,\s*(\d+)\s*\)", low)
        if cnt:
            kw.append(f"n={cnt.group(1)}, m={cnt.group(2)}")
        ln = re.search(NUM + r"\s*nm\s*(?:long|length)", low)
        if ln and "nanotube" in low or ln and "cnt" in low:
            kw.append(f"length={ln.group(1)}")
        ox = re.search(NUM + r"\s*%\s*(?:o\b|oxygen|oxidi)", low)
        if ox:
            kw.append(f"oxygen={ox.group(1)}")
        lay = re.search(r"\b(\d+|bi|tri)[\s-]*(?:layers?|layer)\b|\b(bilayer|trilayer)\b", low)
        if lay:
            n = {"bi": 2, "tri": 3, "bilayer": 2, "trilayer": 3}.get(lay.group(1) or lay.group(2),
                                                                     None) or int(lay.group(1))
            kw.append(f"layers={n}")
        if re.search(r"\bperiodic\b", low):
            kw.append("periodic=True")
        if kw:
            cmds.append(f"carbon({', '.join(kw)})")
    else:
        shape = None
        for name, words in SHAPES:
            if any(re.search(rf"\b{re.escape(wd)}\b", low) for wd in words):
                shape = name
                break
        if shape is None and new_material:
            shape = "sphere"
        if shape:
            req["shape"] = shape
            kw = []
            diam = re.search(NUM + r"\s*nm\s*(?:in\s*)?(?:diameter|wide|across)", low) or \
                re.search(r"(?:diameter|wide|across)\s*(?:of\s*)?" + NUM + r"\s*nm", low)
            if shape in ("sphere", "rod"):
                d = float(diam.group(1)) if diam else None
                if d is None and nums_nm and shape == "sphere":
                    d = nums_nm[0]
                if d is not None:
                    kw.append(f"diameter={d:g}")
                    req["diameter_nm"] = d
            if shape == "rod":
                ln = re.search(NUM + r"\s*nm\s*(?:long|length)", low) or \
                    re.search(r"length\s*(?:of\s*)?" + NUM + r"\s*nm", low)
                if ln:
                    kw.append(f"length={float(ln.group(1)):g}")
                    req["rod_length_nm"] = float(ln.group(1))
            if shape in ("slab", "channel"):
                size = re.search(NUM + r"\s*(?:x|by|×)\s*" + NUM + r"\s*nm", low)
                if size:
                    kw.append(f"size=({size.group(1)}, {size.group(2)})")
                mil = re.search(r"\((\d)(\d)(\d)\)", low)
                if mil:
                    kw.append(f"miller=({mil.group(1)}, {mil.group(2)}, {mil.group(3)})")
                    req["miller"] = [int(mil.group(i)) for i in (1, 2, 3)]
                lay = re.search(r"(\d+)\s*(?:atomic\s*)?layers?", low)
                if lay:
                    kw.append(f"layers={lay.group(1)}")
            if shape == "channel":
                kw.append('type="pore"' if re.search(r"\bpores?\b", low) else 'type="slit"')
            cmds.append(f'shape("{shape}"' + "".join(f", {k}" for k in kw) + ")")

    box = re.search(r"(?:cub(?:e|ic)\s*box|box)\s*(?:of\s*)?" + NUM + r"\s*nm", low) or \
        re.search(NUM + r"\s*nm\s*cub(?:e|ic)?\s*box", low)
    if box:
        cmds.append(f"box(\"cube\", edge={float(box.group(1)):g})")
        req["box_edge_nm"] = float(box.group(1))
    elif re.search(r"\bno ligands?\b|\bbare\b|\buncoated\b", low):
        pass

    lig_key, smiles = find_ligand(text)
    bare = bool(re.search(r"\bno ligands?\b|\bbare\b|\buncoated\b|\bnaked\b", low))
    if lig_key:
        cmds.append(f'ligand("{lig_key}")')
        req["ligand"] = lig_key
    elif smiles:
        cmds.append(f'ligand(smiles="{smiles}")')
        req["ligand"] = smiles
    cov = []
    count = re.search(r"\b(\d+)\s+(?:[a-z0-9\-]+\s+){0,3}?(?:ligands?|molecules?|chains?|groups?|"
                      r"thiols?|amines?|copies)\b", low)
    if lig_key and not count:  # '60 dodecanethiol'
        count = re.search(r"\b(\d+)\s+(?:x\s+)?" + re.escape(lig_key.split("_")[0]), low)
    dens = re.search(NUM + r"\s*(?:ligands?\s*)?(?:/|per)\s*nm\s*(?:\^?2|²|squared)", low)
    if bare:
        cov.append("count=0")
        req["ligand_count"] = 0
    elif dens:
        cov.append(f"density={float(dens.group(1)):g}")
        req["density_nm2"] = float(dens.group(1))
    elif count:
        cov.append(f"count={int(count.group(1))}")
        req["ligand_count"] = int(count.group(1))
    for word, order in (("fibonacci", "fibonacci"), ("evenly", "even"), ("even", "even"),
                        ("random", "random")):
        if re.search(rf"\b{word}\b", low):
            cov.append(f'order="{order}"')
            req["order"] = order
            break
    for word in ("top", "bottom", "both"):
        if re.search(rf"\b{word}\s*(?:side|face|surface)?s?\b", low) and req.get("shape") in ("slab",):
            cov.append(f'region="{word}"')
            req["region"] = word
            break
    if cov:
        cmds.append(f"coverage({', '.join(cov)})")

    fmts = [f for f in FORMATS if re.search(rf"\b{f}\b", low)]
    if fmts:
        cmds.append("export(" + ", ".join(f'"{f}"' for f in fmts) + ")")
        req["formats"] = fmts

    clarify = ""
    if not new_material and not question and re.search(
            r"\b(?:a|an|new|another)\s+(?:[\w.-]+\s+){0,3}?(?:sphere|nanoparticle|particle|nanorod|rod|"
            r"slab|surface|sheet|flake|quantum dot|cluster|nanosphere)s?\b", low):
        clarify = ("Which material should I use? For example gold (Au), silver (Ag), CdSe, TiO2, "
                   "graphene or graphene oxide.")
        cmds, req = [], {}
    if not cmds and not question and not clarify:
        unknown.append("No material, shape or ligand recognised")
    return {"commands": cmds, "requirements": req, "unknown": unknown, "clarify": clarify,
            "build": build or bool(cmds), "question": question and not cmds}


LIGAND_WORDS = re.compile(
    r"\b(ligands?|coat\w*|functionali\w*|thiol\w*|amines?|peg\w*|capp\w*|passivat\w*|"
    r"decorat\w*|graft\w*|attach\w*|smiles|molecules?|chains?|shell of)\b")


def mentions_ligands(text: str) -> bool:
    """Whether the request asks for ligands in any way (a name, a count, 'coated'...)."""
    low = text.lower()
    key, smiles = find_ligand(text)
    return bool(key or smiles or LIGAND_WORDS.search(low))


__all__ = ["parse", "mentions_ligands", "model_edits", "find_material", "find_ligand", "ELEMENT_NAMES"]
