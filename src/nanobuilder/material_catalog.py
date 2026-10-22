"""Grouped catalog of common materials for the material picker.

Four kinds of entries:

- ``preset``: a built-in, offline preset (``material_presets.json``).
- ``mp``: a preloaded Materials Project structure (``data/materials_project``,
  see :mod:`mp_library`); only entries whose CIF is present are listed, and
  nothing is fetched from the network.
- ``carbon``: CarbonConfig settings for a generated carbon nanomaterial.
- Any kind may carry ``geometry`` overrides; the "Nanorods" group uses this to
  set the rod shape and a typical growth axis.

Labels name the material only; sizes come from the Shape tab and the
simulation box. Catalog keys are "Group · Name"; :data:`NAMES` maps the bare
name shown in the picker back to its key.
"""

from __future__ import annotations

import re

from . import mp_library
from .crystal import preset_choices, registry

_ROD = {"shape": "rod", "diameter_nm": 3.0, "rod_length_nm": 10.0, "rod_caps": "hemispherical",
        "rod_axis": (0, 0, 1)}


def _mp(formula, spacegroups):
    hits = mp_library.find(formula, spacegroups)
    return hits[0] if hits else None


def _build() -> dict:
    catalog = {}
    for label, preset_id in preset_choices():
        catalog[f"Elements & simple solids (built-in) · {label}"] = {
            "group": "Elements & simple solids (built-in)", "kind": "preset", "preset_id": preset_id}
    for key, e in mp_library.entries().items():
        catalog[key] = {"group": e["group"], "kind": "mp", "mp_id": e["material_id"],
                        "formula": e["formula"], "spacegroups": [e["spacegroup"]]}
    rods = [
        ("Gold nanorod - [001] axis", {"kind": "preset", "preset_id": "Au_fcc"}),
        ("Silver nanorod - [001] axis", {"kind": "preset", "preset_id": "Ag_fcc"}),
        ("CdSe nanorod - wurtzite c-axis", _mp("CdSe", ["P6_3mc"])),
        ("ZnO nanorod - wurtzite c-axis", _mp("ZnO", ["P6_3mc"])),
    ]
    for label, material in rods:
        if material is None:
            continue
        if "material_id" in material:  # a library entry
            material = {"kind": "mp", "mp_id": material["material_id"],
                        "formula": material["formula"], "spacegroups": [material["spacegroup"]]}
        catalog[f"Nanorods · {label}"] = dict(material, group="Nanorods", geometry=dict(_ROD))
    for label, carbon in CARBON_ENTRIES:
        catalog[f"Carbon nanomaterials · {label}"] = {"group": "Carbon nanomaterials",
                                                      "kind": "carbon", "carbon": carbon}
    return catalog


# Carbon nanomaterials (generated; HierGO-style defects and oxidation). Values are
# CarbonConfig fields; everything else keeps its default.
CARBON_ENTRIES = [
    ("Graphene flake (H-terminated, AMBER-ready)", {"kind": "graphene"}),
    ("Graphene sheet (periodic)", {"kind": "graphene", "periodic": True}),
    ("Bilayer graphene, AB, periodic", {"kind": "graphene", "periodic": True, "layers": 2}),
    ("Graphene oxide flake, 30 % O (fresh 2:1, AMBER GO dataset)",
     {"kind": "graphene", "oxygen_percent": 30.0, "hydroxyl_epoxy_ratio": "2:1"}),
    ("Graphene oxide, aged 900 K, holes + Stone-Wales (HierGO)",
     {"kind": "graphene", "length_nm": 8.0, "width_nm": 8.0, "oxygen_percent": 20.0,
      "oxidation_profile": "aged_900K", "vacancy_percent": 8.0, "stone_wales_percent": 2.0,
      "oxidation_pattern": "edges_first"}),
    ("Heterogeneous GO, tiles of 0-60 % O (HierGO tiling)",
     {"kind": "graphene", "length_nm": 10.0, "width_nm": 10.0, "oxidation_pattern": "tiles",
      "tile_grid": (3, 3), "tile_oxygen_min": 0.0, "tile_oxygen_max": 60.0}),
    ("Nanoporous graphene (2 pores)",
     {"kind": "graphene", "periodic": True, "length_nm": 6.0, "width_nm": 6.0, "pore_count": 2,
      "pore_diameter_nm": 1.2}),
    ("SWCNT (10,10) armchair", {"kind": "nanotube", "cnt_n": 10, "cnt_m": 10}),
    ("SWCNT (17,0) zigzag, periodic", {"kind": "nanotube", "cnt_n": 17, "cnt_m": 0, "periodic": True}),
    ("MWCNT, 3 walls from (5,5)",
     {"kind": "nanotube", "cnt_n": 5, "cnt_m": 5, "cnt_walls": 3}),
    ("Oxidized SWCNT (10,10), carboxylated rims",
     {"kind": "nanotube", "oxygen_percent": 10.0, "oxidation_profile": "aged_900K"}),
    ("Fullerene C60", {"kind": "carbon_sphere", "sphere_type": "fullerene", "fullerene_order": 1}),
    ("Fullerene C240", {"kind": "carbon_sphere", "sphere_type": "fullerene", "fullerene_order": 2}),
    ("Carbon onion C60@C240@C540", {"kind": "carbon_sphere", "sphere_type": "onion",
                                    "onion_shells": 3}),
]

# Typed words that open a carbon entry.
CARBON_WORDS = {
    "graphene": "Graphene flake (H-terminated, AMBER-ready)",
    "go": "Graphene oxide flake, 30 % O (fresh 2:1, AMBER GO dataset)",
    "graphene oxide": "Graphene oxide flake, 30 % O (fresh 2:1, AMBER GO dataset)",
    "rgo": "Graphene oxide, aged 900 K, holes + Stone-Wales (HierGO)",
    "cnt": "SWCNT (10,10) armchair",
    "swcnt": "SWCNT (10,10) armchair",
    "nanotube": "SWCNT (10,10) armchair",
    "carbon nanotube": "SWCNT (10,10) armchair",
    "mwcnt": "MWCNT, 3 walls from (5,5)",
    "c60": "Fullerene C60",
    "fullerene": "Fullerene C60",
    "c240": "Fullerene C240",
    "carbon onion": "Carbon onion C60@C240@C540",
    "onion": "Carbon onion C60@C240@C540",
}


def carbon_entry_for_text(text: str):
    """Catalog key for a typed carbon word (``GO`` -> the graphene-oxide flake), if any."""
    label = CARBON_WORDS.get(str(text or "").strip().lower())
    return f"Carbon nanomaterials · {label}" if label else None


CATALOG = _build()
GROUPS = list(dict.fromkeys(e["group"] for e in CATALOG.values()))


def _names() -> dict:
    """Bare name -> catalog key; a name used in two groups gets its group in brackets."""
    count = {}
    for key in CATALOG:
        name = key.partition(" · ")[2]
        count[name] = count.get(name, 0) + 1
    out = {}
    for key, e in CATALOG.items():
        name = key.partition(" · ")[2]
        out[name if count[name] == 1 else f"{name} ({e['group']})"] = key
    return out


NAMES = _names()
KEY_TO_NAME = {v: k for k, v in NAMES.items()}


def labels(group: str | None = None) -> list:
    """Picker labels (bare names), optionally of one group."""
    return [KEY_TO_NAME[k] for k, e in CATALOG.items() if group in (None, "", e["group"])]


def menu() -> list:
    """[(group, [names])] for the hover menu."""
    return [(g, labels(g)) for g in GROUPS]


def key_for_text(text: str, fuzzy: bool = True):
    """Catalog key for a picker name, a full key, or a carbon word; None if unknown."""
    t = str(text or "").strip()
    if t in CATALOG:
        return t
    if t in NAMES:
        return NAMES[t]
    low = t.lower()
    for name, key in NAMES.items():
        if name.lower() == low:
            return key
    carbon = carbon_entry_for_text(t)
    if carbon:
        return carbon
    if fuzzy and len(low) >= 3:  # a word from a name: "anatase", "magnetite", "nitinol"
        words = [(name, key) for name, key in NAMES.items()
                 if any(wd.startswith(low) for wd in re.split(r"[\s()\-,]+", name.lower()))]
        if words:
            return words[0][1]
    return None


def library_for_text(text: str) -> list:
    """Catalog keys of preloaded Materials Project entries with this formula (TiO2 -> both)."""
    ids = {e["material_id"] for e in mp_library.find(text)}
    return [k for k, e in CATALOG.items() if e["kind"] == "mp" and e.get("mp_id") in ids
            and e["group"] != "Nanorods"]


def preset_for_text(text: str):
    """Built-in preset whose formula matches what was typed (``Au`` -> Au_fcc), if any."""
    t = str(text or "").strip().lower()
    for key, entry in registry()["presets"].items():
        if entry["formula"].lower() == t or key.lower() == t:
            return key
    return None
