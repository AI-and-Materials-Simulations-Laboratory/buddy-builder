"""Materials Project client, tested against responses in MP's documented format.

No network: a fake fetch function returns payloads shaped like the real
``/materials/summary/`` responses (``data`` list of SummaryDoc fields, with a
pymatgen Structure dict).
"""

import json

import numpy as np
import pytest
from ase.build import bulk

from nanobuilder import materials_project as mp
from nanobuilder.crystal import build_crystal
from nanobuilder.importers import atoms_to_cif_bytes
from nanobuilder.pipeline import run_pipeline
from nanobuilder.schema import BuilderError, MaterialConfig, RunConfig


def _structure_dict(atoms):
    """pymatgen Structure.as_dict() layout."""
    cell = np.array(atoms.cell)
    return {"@module": "pymatgen.core.structure", "@class": "Structure", "charge": 0,
            "lattice": {"matrix": cell.tolist(), "pbc": [True, True, True]},
            "sites": [{"species": [{"element": s, "occu": 1}], "abc": f.tolist(),
                       "xyz": p.tolist(), "label": s, "properties": {}}
                      for s, f, p in zip(atoms.get_chemical_symbols(),
                                         atoms.get_scaled_positions(), atoms.get_positions())]}


def _l12_aucu3_primitive():
    # L1_2 AuCu3 (Pm-3m): already primitive = conventional, 4 atoms
    a = 3.75
    from ase import Atoms
    return Atoms("AuCu3", scaled_positions=[[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
                 cell=[a, a, a], pbc=True)


DOCS = {
    "mp-2258": {"material_id": "mp-2258", "formula_pretty": "AuCu3",
                "symmetry": {"crystal_system": "Cubic", "symbol": "Pm-3m", "number": 221},
                "energy_above_hull": 0.0, "is_stable": True, "theoretical": False, "nsites": 4,
                "deprecated": False, "structure": _structure_dict(_l12_aucu3_primitive())},
    "mp-9999901": {"material_id": "mp-9999901", "formula_pretty": "AuCu3",
                   "symmetry": {"crystal_system": "Tetragonal", "symbol": "I4/mmm", "number": 139},
                   "energy_above_hull": 0.012, "is_stable": False, "theoretical": True,
                   "nsites": 8, "deprecated": False},
    "mp-81": {"material_id": "mp-81", "formula_pretty": "Au",
              "symmetry": {"crystal_system": "Cubic", "symbol": "Fm-3m", "number": 225},
              "energy_above_hull": 0.0, "is_stable": True, "theoretical": False, "nsites": 1,
              "deprecated": False, "structure": _structure_dict(bulk("Au", "fcc", a=4.17))},
}


class FakeMP:
    def __init__(self):
        self.calls = []

    def __call__(self, params, key):
        self.calls.append((dict(params), key))
        assert key == "test-key"
        fields = params["_fields"].split(",")
        if "material_ids" in params:
            docs = [DOCS[params["material_ids"]]]
        else:
            docs = [d for d in DOCS.values() if d["formula_pretty"] == params.get("formula")]
        # The real API returns only the requested fields.
        return json.loads(json.dumps({"data": [{k: d[k] for k in fields if k in d} for d in docs],
                                      "meta": {"total_doc": len(docs)}}))


@pytest.fixture(autouse=True)
def _fresh_cache():
    mp.clear_cache()
    yield
    mp.clear_cache()


def test_parse_query_variants():
    assert mp.parse_query("AuCu3") == {"formula": "AuCu3"}
    assert mp.parse_query("Au-Cu") == {"chemsys": "Au-Cu"}
    assert mp.parse_query(" mp-149 ") == {"material_ids": "mp-149"}
    for bad in ("", "gold", "Au-cu"):
        with pytest.raises(BuilderError):
            mp.parse_query(bad)


def test_search_sends_documented_parameters_and_sorts_by_stability():
    fake = FakeMP()
    hits = mp.search("AuCu3", key="test-key", fetch=fake)
    params, _ = fake.calls[0]
    assert params["formula"] == "AuCu3" and params["deprecated"] == "false"
    assert set(params["_fields"].split(",")) >= {"material_id", "symmetry", "energy_above_hull"}
    assert [h.material_id for h in hits] == ["mp-2258", "mp-9999901"]
    assert "stable" in hits[0].label and "theoretical" in hits[1].label
    assert mp.best_match(hits, ["I4/mmm"]).material_id == "mp-9999901"
    mp.search("AuCu3", key="test-key", fetch=fake)
    assert len(fake.calls) == 1  # cached


def test_missing_key_is_actionable(monkeypatch):
    monkeypatch.delenv("MP_API_KEY", raising=False)
    with pytest.raises(BuilderError, match="API key"):
        mp.search("Au", fetch=FakeMP())


def test_fetch_structure_returns_conventional_cell():
    atoms, info = mp.fetch_structure("mp-81", key="test-key", fetch=FakeMP())
    # MP serves the 1-atom primitive fcc cell; the tool works on the 4-atom cubic cell
    assert len(atoms) == 4
    assert np.allclose(atoms.cell.cellpar(), [4.17, 4.17, 4.17, 90, 90, 90], atol=1e-3)
    assert "spglib" in info["normalization"] and info["spacegroup"] == "Fm-3m"


def test_partial_occupancy_is_rejected():
    d = _structure_dict(bulk("Cu", "fcc", a=3.6))
    d["sites"][0]["species"] = [{"element": "Au", "occu": 0.5}, {"element": "Cu", "occu": 0.5}]
    with pytest.raises(BuilderError, match="partially occupied"):
        mp.structure_dict_to_atoms(d)


def test_network_failure_message(monkeypatch):
    import urllib.error

    def boom(*a, **k):
        raise urllib.error.URLError("Name or service not known")
    monkeypatch.setattr(mp.urllib.request, "urlopen", boom)
    with pytest.raises(BuilderError, match="Could not reach Materials Project"):
        mp.search("Au", key="k")


def test_http_403_message(monkeypatch):
    import urllib.error

    def denied(*a, **k):
        raise urllib.error.HTTPError(mp.API_URL, 403, "Forbidden", {}, None)
    monkeypatch.setattr(mp.urllib.request, "urlopen", denied)
    with pytest.raises(BuilderError, match="rejected the API key"):
        mp.search("Au", key="bad")


def test_mp_crystal_builds_an_alloy_nanoparticle_end_to_end():
    atoms, info = mp.fetch_structure("mp-2258", key="test-key", fetch=FakeMP())
    cfg = RunConfig.from_dict({
        "material": {"source": "materials_project", "mp_id": "mp-2258"},
        "geometry": {"shape": "sphere", "diameter_nm": 2.5, "vacuum_A": 20},
        "ligand": {"smiles": "SCCO", "anchor_atom_index": 0, "binding_template": "geometric"},
        "coverage": {"mode": "count", "target": 8}},
        uploaded_bytes=atoms_to_cif_bytes(atoms))
    crystal = build_crystal(cfg.material)
    assert crystal.source_type == "materials_project" and "mp-2258" in crystal.structure_label
    assert crystal.anchor_species_candidates == ["Au", "Cu"]
    r = run_pipeline(cfg)
    assert r.report.status == "ok" and r.assembly.achieved_ligands == 8
    assert {"Au", "Cu"} <= set(r.assembly.elements)
    d = cfg.to_dict()["material"]
    assert d["mp_id"] == "mp-2258" and d["source"] == "materials_project"


def test_mp_source_without_structure_is_reported():
    with pytest.raises(BuilderError, match="not in the preloaded"):
        MaterialConfig(source="materials_project", mp_id="mp-81").validate()


def test_search_asks_server_for_most_stable_first_and_by_spacegroup():
    fake = FakeMP()
    mp.search("AuCu3", key="test-key", fetch=fake)
    mp.search("AuCu3", key="test-key", fetch=fake, spacegroup="Pm-3m")
    (p1, _), (p2, _) = fake.calls
    assert p1["_sort_fields"] == "energy_above_hull" and "spacegroup_symbol" not in p1
    assert p2["spacegroup_symbol"] == "Pm-3m"


def test_conventional_cell_tolerates_small_coordinate_noise():
    """Relaxed DFT structures carry tiny distortions; symmetry must still be found."""
    from ase.build import bulk
    rng = np.random.default_rng(0)
    nacl = bulk("NaCl", "rocksalt", a=5.64)  # primitive, 2 atoms
    nacl.positions += rng.normal(scale=0.002, size=nacl.positions.shape)
    out, note = mp.conventional_cell(nacl, "Fm-3m")
    assert len(out) == 8 and "Fm-3m" in note and "WARNING" not in note
    mg = bulk("Mg", "hcp", a=3.21, c=5.21)
    mg.positions += rng.normal(scale=0.002, size=mg.positions.shape)
    out, note = mp.conventional_cell(mg, "P6_3/mmc")
    assert len(out) == 2 and "WARNING" not in note


def test_conventional_cell_flags_a_space_group_mismatch():
    from ase.build import bulk
    out, note = mp.conventional_cell(bulk("Cu", "fcc", a=3.6), "Im-3m")
    assert "WARNING" in note
