"""Material registry and bulk-crystal construction (spec stage A).

``build_crystal()`` is the only entry point the rest of the pipeline uses; it
returns the same :class:`CrystalModel` whether the crystal came from a preset
or from an uploaded file.
"""

from __future__ import annotations

import copy
import json
from functools import lru_cache
from pathlib import Path

from ase.build import bulk

from .importers import read_crystal_bytes
from .schema import BuilderError, CrystalModel, MaterialConfig

REGISTRY_PATH = Path(__file__).with_name("material_presets.json")
CUBIC_PROTOTYPES = {"fcc", "bcc", "diamond", "rocksalt"}
OVERRIDE_WARN_FRACTION = 0.05
OVERRIDE_MAX_FRACTION = 0.25


@lru_cache(maxsize=1)
def _load_registry() -> dict:
    return json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))


def registry() -> dict:
    """A copy of the preset registry (callers may not modify the cached one)."""
    return copy.deepcopy(_load_registry())


def default_preset_id() -> str:
    return _load_registry()["default"]


def preset(preset_id: str) -> dict:
    presets = _load_registry()["presets"]
    if preset_id not in presets:
        raise BuilderError(f"Unsupported material/structure pair {preset_id!r}; choose one of "
                           f"{', '.join(presets)}")
    return copy.deepcopy(presets[preset_id])


def preset_choices() -> list:
    """(label, preset_id) pairs in registry order, for the GUI dropdown."""
    presets = _load_registry()["presets"]
    return [(p["label"], key) for key, p in presets.items()]


def build_crystal(material_config: MaterialConfig) -> CrystalModel:
    """Load the validated preset or user-supplied bulk crystal."""
    material_config.validate()
    if material_config.source == "preset":
        return _build_preset(material_config)
    if material_config.source == "materials_project":
        return _build_materials_project(material_config)
    return _build_uploaded(material_config)


def _build_preset(cfg: MaterialConfig) -> CrystalModel:
    entry = preset(cfg.preset_id)
    lattice = dict(entry["lattice_A"])
    warnings = []
    for key, value in (cfg.lattice_overrides_A or {}).items():
        if key not in entry["free_parameters"]:
            raise BuilderError(f"Lattice parameter {key!r} is fixed by the {entry['prototype']} "
                               f"prototype; editable: {', '.join(entry['free_parameters'])}")
        reference = lattice[key]
        change = abs(float(value) - reference) / reference
        if change > OVERRIDE_MAX_FRACTION:
            raise BuilderError(f"Lattice parameter {key} = {value} A is more than "
                               f"{OVERRIDE_MAX_FRACTION:.0%} away from the registry value "
                               f"{reference} A")
        if change > OVERRIDE_WARN_FRACTION:
            warnings.append(f"Lattice parameter {key} = {value} A differs from the registry value "
                            f"{reference} A by {change:.1%}")
        lattice[key] = float(value)
    overrides = {k: float(v) for k, v in (cfg.lattice_overrides_A or {}).items()
                 if float(v) != entry["lattice_A"][k]}
    if overrides:
        warnings.append("User-modified lattice parameters: " +
                        ", ".join(f"{k} = {v} A" for k, v in overrides.items()))

    kwargs = {"a": lattice["a"]}
    if "c" in lattice:
        kwargs["c"] = lattice["c"]
    # cubic=True gives the conventional cell, so (hkl) has its usual meaning.
    atoms = bulk(entry["formula"], entry["prototype"],
                 cubic=entry["prototype"] in CUBIC_PROTOTYPES, **kwargs)
    if entry.get("notes"):
        warnings.append(entry["notes"])
    return CrystalModel(
        atoms=atoms,
        source_type="preset",
        preset_id=cfg.preset_id,
        lattice_overrides_A=overrides,
        structure_label=entry["label"],
        space_group=entry.get("space_group"),
        prototype=entry["prototype"],
        anchor_species_candidates=sorted(set(atoms.get_chemical_symbols())),
        lattice_source=entry["source"],
        normalization="conventional cell from ase.build.bulk"
                      + (" (cubic=True)" if entry["prototype"] in CUBIC_PROTOTYPES else ""),
        validation_warnings=warnings,
    )


def _build_uploaded(cfg: MaterialConfig) -> CrystalModel:
    atoms, info = read_crystal_bytes(cfg.uploaded_bytes, cfg.filename or "upload",
                                     cfg.file_format)
    cfg.sha256 = info["sha256"]
    cfg.filename = info["filename"]
    warnings = list(info["warnings"])
    warnings.append("Uploaded crystal: no reviewed surface chemistry or force field; only "
                    "geometric ligand placement is offered")
    return CrystalModel(
        atoms=atoms,
        source_type="uploaded_file",
        source_format=cfg.file_format,
        source_filename=info["filename"],
        source_sha256=info["sha256"],
        source_bytes=cfg.uploaded_bytes,
        structure_label=f"{atoms.get_chemical_formula()} (uploaded {info['filename']})",
        anchor_species_candidates=sorted(set(atoms.get_chemical_symbols())),
        lattice_source=f"uploaded file {info['filename']}",
        normalization=info["normalization"],
        validation_warnings=warnings,
    )


def _build_materials_project(cfg: MaterialConfig) -> CrystalModel:
    """Crystal fetched from Materials Project, held as a CIF of its conventional cell.

    The CIF goes through the same checks as an upload, so a config.json plus
    that CIF rebuilds the model without network access.
    """
    cfg.file_format = "cif"
    cfg.filename = cfg.filename or f"{cfg.mp_id}.cif"
    atoms, info = read_crystal_bytes(cfg.uploaded_bytes, cfg.filename, "cif")
    cfg.sha256 = info["sha256"]
    formula = atoms.get_chemical_formula(mode="metal", empirical=True)
    from .mp_library import entry as library_entry
    known = library_entry(cfg.mp_id)
    name = known["name"] if known else formula
    warnings = list(info["warnings"])
    warnings.append(f"Materials Project {cfg.mp_id}: computed (DFT) structure, CC BY 4.0; cite "
                    "Jain et al., APL Materials 1, 011002 (2013). Only geometric ligand placement "
                    "is offered (no reviewed surface chemistry)")
    return CrystalModel(
        atoms=atoms,
        source_type="materials_project",
        source_format="cif",
        source_filename=cfg.filename,
        source_sha256=info["sha256"],
        source_bytes=cfg.uploaded_bytes,
        structure_label=f"{name} ({cfg.mp_id}, Materials Project)",
        space_group=known["spacegroup"] if known else None,
        anchor_species_candidates=sorted(set(atoms.get_chemical_symbols())),
        lattice_source=f"Materials Project {cfg.mp_id} (https://materialsproject.org/materials/{cfg.mp_id})",
        normalization="conventional standard cell (spglib) of the Materials Project structure",
        validation_warnings=warnings,
    )
