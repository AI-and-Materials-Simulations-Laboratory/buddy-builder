"""Preloaded Materials Project structures (offline).

``data/materials_project`` holds conventional cells exported from Materials
Project as CIF plus ``catalog.json`` (ids, formula, space group, stability,
licence and citation). The tool reads these files instead of calling the MP
API, so every entry works without a key or network access; an entry whose file
is missing or does not load is simply not offered.

Data: Materials Project, CC BY 4.0; cite A. Jain et al., APL Materials 1,
011002 (2013).
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data" / "materials_project"


@lru_cache(maxsize=1)
def _catalog() -> dict:
    path = DATA_DIR / "catalog.json"
    if not path.is_file():
        return {"entries": {}}
    return json.loads(path.read_text(encoding="utf-8"))


def entries() -> dict:
    """{catalog label "Group · Name": entry} for every entry whose CIF file exists."""
    out = {}
    for key, e in _catalog().get("entries", {}).items():
        if (DATA_DIR / e.get("file", "")).is_file():
            group, _, name = key.partition(" · ")
            out[key] = dict(e, group=group, name=name or key)
    return out


def by_id() -> dict:
    return {e["material_id"]: e for e in entries().values()}


def entry(mp_id: str):
    return by_id().get(mp_id)


def cif_bytes(mp_id: str):
    """The preloaded CIF of ``mp_id``, or None when it is not in the library."""
    e = entry(mp_id)
    if e is None:
        return None
    return (DATA_DIR / e["file"]).read_bytes()


def find(formula: str = "", spacegroups=()) -> list:
    """Library entries with this reduced formula (any order of elements) and space group."""
    from ase.formula import Formula
    try:
        want = Formula(str(formula).strip()).reduce()[0] if formula else None
    except Exception:  # not a formula
        return []
    out = []
    for e in entries().values():
        try:
            have = Formula(e["formula"]).reduce()[0]
        except Exception:
            continue
        if want is not None and have != want:
            continue
        if spacegroups and e["spacegroup"] not in spacegroups:
            continue
        out.append(e)
    return out


def citation() -> str:
    return _catalog().get("citation", "Materials Project, CC BY 4.0, https://materialsproject.org")


__all__ = ["entries", "entry", "by_id", "cif_bytes", "find", "citation", "DATA_DIR"]
