"""Live Materials Project search (https://materialsproject.org), standard library only.

Talks to the public REST API the official ``mp-api`` client uses:
``GET https://api.materialsproject.org/materials/summary/`` with the user's key
in the ``x-api-key`` header, filters such as ``formula``, ``chemsys`` or
``material_ids``, ``_fields`` to select fields and ``_limit`` for paging.
Responses are ``{"data": [...], "meta": {...}}``; structures are pymatgen
``Structure`` dictionaries (``lattice.matrix`` plus ``sites`` with ``species``,
``abc`` and ``xyz``), converted here to ASE without pymatgen.

MP serves primitive cells. Miller indices and rod axes in this tool refer to
the *conventional* cell, so fetched structures are standardized with spglib
(``to_primitive=False``). Without spglib the primitive cell is used and a
warning says so.

Data from Materials Project is licensed CC BY 4.0; cite it (see
``CITATION`` below), which the export README and manifest record.
"""

from __future__ import annotations

import json
import os
import re
import warnings
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

import numpy as np
from ase import Atoms

from .schema import BuilderError

API_URL = "https://api.materialsproject.org/materials/summary/"
KEY_ENV = "MP_API_KEY"
TIMEOUT_S = 25
CITATION = ("A. Jain et al., The Materials Project: A materials genome approach to accelerating "
            "materials innovation, APL Materials 1, 011002 (2013). Data: CC BY 4.0, "
            "https://materialsproject.org")
SEARCH_FIELDS = ("material_id", "formula_pretty", "symmetry", "energy_above_hull", "is_stable",
                 "theoretical", "nsites", "deprecated")
_MP_ID = re.compile(r"^(mp|mvc)-\d+$", re.IGNORECASE)

_search_cache: dict = {}
_structure_cache: dict = {}


@dataclass
class MPHit:
    material_id: str
    formula: str
    spacegroup: str
    spacegroup_number: int | None
    crystal_system: str
    energy_above_hull: float | None
    is_stable: bool
    theoretical: bool
    nsites: int | None
    raw: dict = field(default_factory=dict, repr=False)

    @property
    def label(self) -> str:
        ehull = ("-" if self.energy_above_hull is None
                 else f"{1000 * self.energy_above_hull:.0f} meV/atom")
        flags = ["stable" if self.is_stable else f"E_hull {ehull}"]
        if self.theoretical:
            flags.append("theoretical")
        return (f"{self.material_id}  {self.formula}  {self.spacegroup} "
                f"({self.crystal_system})  {', '.join(flags)}")


def api_key(explicit: str | None = None) -> str:
    key = (explicit or os.environ.get(KEY_ENV, "")).strip()
    if not key:
        raise BuilderError("A Materials Project API key is needed: paste yours (free, from "
                           "https://materialsproject.org/api) or set the MP_API_KEY environment "
                           "variable")
    return key


def parse_query(text: str) -> dict:
    """Turn what the user typed into MP filters.

    ``mp-149`` -> material_ids; ``Au-Cu`` -> chemsys (any composition of those
    elements); ``AuCu3`` / ``TiO2`` -> formula (exact reduced formula).
    """
    q = re.sub(r"\s+", "", str(text or ""))
    if not q:
        raise BuilderError("Type a formula (e.g. AuCu3), a chemical system (Au-Cu) or an mp-id")
    if _MP_ID.match(q):
        return {"material_ids": q.lower()}
    if "-" in q:
        parts = [p for p in q.split("-") if p]
        if not all(re.fullmatch(r"[A-Z][a-z]?", p) for p in parts):
            raise BuilderError(f"{text!r} is not a chemical system; use element symbols like Au-Cu")
        return {"chemsys": "-".join(parts)}
    if not re.fullmatch(r"([A-Z][a-z]?\d*(\.\d+)?)+", q):
        raise BuilderError(f"{text!r} is not a chemical formula; use e.g. Au, AuCu3, TiO2")
    return {"formula": q}


def _get(params: dict, key: str, timeout: float = TIMEOUT_S) -> dict:
    url = API_URL + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"x-api-key": key, "accept": "application/json",
                                                   "user-agent": "nanobuilder (nanoHUB tool)"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code in (401, 403):
            raise BuilderError("Materials Project rejected the API key (HTTP "
                               f"{exc.code}); check it at https://materialsproject.org/api") from None
        if exc.code == 429:
            raise BuilderError("Materials Project rate limit reached; wait a moment and retry") from None
        raise BuilderError(f"Materials Project returned HTTP {exc.code}") from None
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        reason = getattr(exc, "reason", exc)
        raise BuilderError("Could not reach Materials Project (api.materialsproject.org): "
                           f"{reason}. This needs outbound internet from the server running the "
                           "tool; you can still download a CIF from materialsproject.org and "
                           "upload it") from None
    except json.JSONDecodeError:
        raise BuilderError("Materials Project sent a response that is not JSON") from None


def _hit(doc: dict) -> MPHit:
    sym = doc.get("symmetry") or {}
    return MPHit(material_id=str(doc.get("material_id")), formula=str(doc.get("formula_pretty")),
                 spacegroup=str(sym.get("symbol") or "?"), spacegroup_number=sym.get("number"),
                 crystal_system=str(sym.get("crystal_system") or "?"),
                 energy_above_hull=doc.get("energy_above_hull"),
                 is_stable=bool(doc.get("is_stable")), theoretical=bool(doc.get("theoretical")),
                 nsites=doc.get("nsites"), raw=doc)


def search(text: str, key: str | None = None, limit: int = 50, include_unstable: bool = True,
           fetch=None, spacegroup: str | None = None) -> list:
    """Search MP; returns MPHit list, stable/most stable first. ``fetch`` is for tests.

    The server sorts by energy above hull, so ``limit`` keeps the most stable
    entries; ``spacegroup`` restricts to one structure (e.g. TiO2 anatase).
    """
    filters = parse_query(text)
    params = dict(filters, _fields=",".join(SEARCH_FIELDS), _limit=int(limit), deprecated="false",
                  _sort_fields="energy_above_hull")
    if spacegroup:
        params["spacegroup_symbol"] = spacegroup
    if not include_unstable:
        params["is_stable"] = "true"
    cache_key = tuple(sorted(params.items()))
    if cache_key not in _search_cache:
        payload = (fetch or _get)(params, api_key(key))
        _search_cache[cache_key] = [_hit(d) for d in payload.get("data", [])]
    hits = list(_search_cache[cache_key])
    hits.sort(key=lambda h: (not h.is_stable, h.energy_above_hull if h.energy_above_hull is not None
                             else 1e9, h.theoretical, h.nsites or 0))
    return hits


def structure_dict_to_atoms(struct: dict) -> Atoms:
    """pymatgen Structure dict -> ASE Atoms (ordered structures only)."""
    try:
        cell = np.array(struct["lattice"]["matrix"], dtype=float)
        symbols, positions = [], []
        for site in struct["sites"]:
            species = site.get("species") or []
            if len(species) != 1 or abs(float(species[0].get("occu", 1.0)) - 1.0) > 1e-6:
                raise BuilderError("This Materials Project entry has partially occupied sites; "
                                   "only ordered structures are supported")
            element = re.sub(r"[^A-Za-z]", "", str(species[0]["element"]))
            symbols.append(element)
            if site.get("xyz") is not None:
                positions.append(np.array(site["xyz"], dtype=float))
            else:
                positions.append(np.array(site["abc"], dtype=float) @ cell)
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        raise BuilderError(f"Unexpected structure format from Materials Project: {exc}") from None
    return Atoms(symbols, positions=np.array(positions), cell=cell, pbc=True)


SYMPREC_A = 0.1  # spglib tolerance; Materials Project's own symmetry analysis uses 0.1 A


def _spglib_call(fn, *args, **kwargs):
    """Call spglib quietly; None on failure (old and new spglib error handling)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        try:
            return fn(*args, **kwargs)
        except Exception:  # spglib >= 2.5 may raise instead of returning None
            return None


def conventional_cell(atoms: Atoms, expected_spacegroup: str | None = None) -> tuple:
    """(atoms, note): spglib-standardized conventional cell, or the input if spglib is missing.

    When ``expected_spacegroup`` (as reported by Materials Project) is given and
    spglib finds a different one, the note says so, because Miller indices and
    rod axes would then refer to a different setting than the user expects.
    """
    try:
        import spglib
    except ImportError:
        return atoms, ("spglib is not installed: the primitive cell from Materials Project is used, "
                       "so Miller indices refer to the primitive cell")
    cell = (np.array(atoms.cell), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    std = _spglib_call(spglib.standardize_cell, cell, to_primitive=False, no_idealize=False,
                       symprec=SYMPREC_A)
    if std is None:
        return atoms, "spglib could not standardize the cell; it is used as served"
    lattice, frac, numbers = std
    out = Atoms(numbers=numbers, scaled_positions=frac, cell=lattice, pbc=True)
    out.wrap()
    dataset = _spglib_call(spglib.get_symmetry_dataset, cell, symprec=SYMPREC_A)
    if isinstance(dataset, dict):  # spglib < 2.5 returns a dict
        sg = dataset.get("international")
    else:
        sg = getattr(dataset, "international", None)
    note = (f"Conventional standard cell from spglib ({len(atoms)} -> {len(out)} atoms"
            + (f", space group {sg}" if sg else "") + ")")
    if expected_spacegroup and sg and sg.replace(" ", "") != expected_spacegroup.replace(" ", ""):
        note += (f"; WARNING: Materials Project reports {expected_spacegroup}, so check that "
                 "Miller indices and rod axes refer to the cell you expect")
    return out, note


def fetch_structure(material_id: str, key: str | None = None, fetch=None) -> tuple:
    """(conventional-cell Atoms, info dict) for one MP entry."""
    mid = str(material_id).strip().lower()
    if not _MP_ID.match(mid):
        raise BuilderError(f"{material_id!r} is not a Materials Project ID (like mp-149)")
    if mid not in _structure_cache:
        payload = (fetch or _get)({"material_ids": mid, "_fields": ",".join(SEARCH_FIELDS + ("structure",)),
                                   "_limit": 1}, api_key(key))
        docs = payload.get("data", [])
        if not docs or "structure" not in docs[0]:
            raise BuilderError(f"Materials Project has no structure for {mid}")
        _structure_cache[mid] = docs[0]
    doc = _structure_cache[mid]
    primitive = structure_dict_to_atoms(doc["structure"])
    hit = _hit(doc)
    atoms, note = conventional_cell(primitive, hit.spacegroup if hit.spacegroup != "?" else None)
    return atoms, {"material_id": mid, "formula": hit.formula, "spacegroup": hit.spacegroup,
                   "crystal_system": hit.crystal_system, "energy_above_hull": hit.energy_above_hull,
                   "is_stable": hit.is_stable, "theoretical": hit.theoretical,
                   "normalization": note, "citation": CITATION, "label": hit.label}


def best_match(hits: list, spacegroups=()) -> MPHit | None:
    """First hit (already stability-sorted) whose space group is one of ``spacegroups``."""
    wanted = {s.replace(" ", "") for s in spacegroups}
    for h in hits:
        if not wanted or h.spacegroup.replace(" ", "") in wanted:
            return h
    return None


def clear_cache() -> None:
    _search_cache.clear()
    _structure_cache.clear()
