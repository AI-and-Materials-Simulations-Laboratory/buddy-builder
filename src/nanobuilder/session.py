"""Save session: a build as Python console commands that rebuild it.

``session.py`` (in every export and under Downloads) holds

- the build in short commands, as comments for people to read;
- ``settings({...})`` with every setting that differs from the defaults
  (a Python literal: the console parses it with :mod:`ast`, nothing is
  executed), plus the crystal file text when the material was uploaded;
- ``build()`` and any hand edits made after Generate (``remove_ligand`` /
  ``add_ligand``), in order.

Pasting the file into the Python console, or loading it with *Load session*,
puts every selection back and rebuilds the same model (the build is seeded,
so ligand numbers and atom indices match).
"""

from __future__ import annotations

import json
import pprint
from datetime import datetime, timezone

from . import __version__
from .schema import BuilderError, RunConfig

ORDERS = {"fibonacci": "fibonacci", "even": "evenly_spaced", "evenly_spaced": "evenly_spaced",
          "random": "distributed_random", "distributed_random": "distributed_random"}
ORDER_NAMES = {"fibonacci": "fibonacci", "evenly_spaced": "even", "distributed_random": "random"}
MAX_EMBED_BYTES = 2_000_000


def fmt(v) -> str:
    if isinstance(v, (list, tuple)):
        return "(" + ", ".join(fmt(x) for x in v) + ("," if len(v) == 1 else "") + ")"
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, float) and v.is_integer():
        return f"{v:g}"
    return repr(v)


def call(name: str, *args, **kw) -> str:
    parts = [fmt(a) for a in args] + [f"{k}={fmt(v)}" for k, v in kw.items() if v is not None]
    return f"{name}({', '.join(parts)})"


def short_commands(cfg: RunConfig, with_ligand: bool = True) -> list:
    """The main settings of ``cfg`` as console commands (readable, not exhaustive)."""
    m, g, lg, c = cfg.material, cfg.geometry, cfg.ligand, cfg.coverage
    out = []
    if m.source == "preset":
        out.append(call("material", preset=m.preset_id))
    elif m.source == "materials_project":
        out.append(call("material", m.mp_id))
    elif m.source == "carbon":
        k = cfg.carbon
        if k.kind == "graphene":
            out.append(call("carbon", kind="graphene", size=[k.length_nm, k.width_nm],
                            layers=k.layers, oxygen=k.oxygen_percent or None,
                            vacancy=k.vacancy_percent or None, periodic=k.periodic or None))
        elif k.kind == "nanotube":
            out.append(call("carbon", kind="nanotube", n=k.cnt_n, m=k.cnt_m, length=k.cnt_length_nm,
                            walls=k.cnt_walls, oxygen=k.oxygen_percent or None,
                            periodic=k.periodic or None))
        else:
            out.append(call("carbon", kind=k.sphere_type, fullerene=k.fullerene_order
                            if k.sphere_type == "fullerene" else None,
                            shells=k.onion_shells if k.sphere_type == "onion" else None))
    else:
        out.append(f"# material: uploaded file {m.filename} (its text is inside settings below)")
    user_box = g.box_lengths_nm is not None and any(float(x) > 0 for x in g.box_lengths_nm)
    if m.source != "carbon":
        if g.shape == "sphere":
            out.append(call("shape", "sphere", diameter=g.diameter_nm))
        elif g.shape == "rod":
            out.append(call("shape", "rod", diameter=g.diameter_nm, length=g.rod_length_nm,
                            axis=list(g.rod_axis)))
        else:
            kw = dict(size=None if user_box else [g.slab_length_nm, g.slab_width_nm],
                      miller=list(g.miller), layers=g.layers if g.thickness_mode == "layers" else None)
            if g.shape == "channel":
                kw.update(type=g.channel_type, width=g.channel_width_nm, count=g.channel_count)
            out.append(call("shape", g.shape, **kw))
    if user_box and g.box_shape == "cube":
        out.append(call("box", "cube", edge=float(g.box_lengths_nm[0])))
    elif user_box:
        out.append(call("box", g.box_shape if g.box_shape != "auto" else "rectangular",
                        size=[float(x) for x in g.box_lengths_nm]))
    elif g.box_shape != "auto":
        out.append(call("box", g.box_shape))
    if with_ligand and (c.target > 0 or c.mode == "density"):
        out.append(call("ligand", lg.library_id) if lg.library_id else
                   call("ligand", smiles=lg.smiles, anchor=lg.anchor_atom_index,
                        template=lg.binding_template))
        out.append(call("coverage", count=c.target if c.mode == "count" else None,
                        density=c.density_nm2 if c.mode == "density" else None,
                        region=c.region, order=ORDER_NAMES.get(c.distribution)))
    return out


def _diff(value, default):
    if isinstance(value, dict) and isinstance(default, dict):
        out = {}
        for k, v in value.items():
            if k not in default:
                out[k] = v
            else:
                d = _diff(v, default[k])
                if d is not _SAME:
                    out[k] = d
        return out if out else _SAME
    if isinstance(value, (list, tuple)) and isinstance(default, (list, tuple)):
        return _SAME if list(value) == list(default) else value
    return _SAME if value == default else value


_SAME = object()


def non_default(cfg: RunConfig) -> dict:
    """Every setting of ``cfg`` that differs from RunConfig() defaults, by section."""
    d = cfg.to_dict()
    base = RunConfig().to_dict()
    out = {}
    for section, value in d.items():
        if section == "format_version":
            continue
        diff = _diff(value, base.get(section, {}))
        if diff is not _SAME and diff:
            out[section] = diff
    material = dict(out.get("material", {}))
    material["source"] = cfg.material.source  # always explicit
    if cfg.material.source == "preset":
        material["preset_id"] = cfg.material.preset_id
    for key in ("sha256",):
        material.pop(key, None)
    out["material"] = material
    if cfg.material.source != "carbon":
        out.pop("carbon", None)
    return {k: out[k] for k in ("material", "geometry", "carbon", "ligand", "coverage", "export")
            if k in out}


def session_text(cfg: RunConfig, assembly=None, name: str = "") -> str:
    """The session file: comments, settings({...}), build() and the hand edits."""
    settings = non_default(cfg)
    if cfg.material.source == "uploaded_file" and cfg.material.uploaded_bytes:
        data = cfg.material.uploaded_bytes
        if len(data) <= MAX_EMBED_BYTES:
            try:
                settings["material"]["file_text"] = data.decode("utf-8")
            except UnicodeDecodeError:
                pass
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [f"# Buddy Builder session{(' ' + repr(name)) if name else ''}, saved {now} "
             f"(nanobuilder {__version__})",
             "# Rebuild: paste this whole file into the Python console, or use Load session.",
             "#",
             "# The build in short commands:"]
    lines += [f"#   {c}" for c in short_commands(cfg)]
    lines += ["#", "# Every setting that differs from the defaults (restores the exact build):",
              "settings(" + pprint.pformat(settings, indent=1, width=96, sort_dicts=False) + ")",
              "build()"]
    edits = (assembly.provenance.get("manual_edits", []) if assembly is not None else [])
    for k, e in enumerate(edits, start=1):
        if e.get("action") == "remove_ligands":
            lines.append(call("remove_ligand", *e["ligands"]) + f"    # hand edit {k}")
        elif e.get("action") == "add_ligand":
            lines.append(call("add_ligand", e["surface_atom"]) + f"    # hand edit {k}")
    if assembly is not None:
        lines.append(f"# Result: {assembly.n_atoms} atoms, {assembly.achieved_ligands} ligands")
    return "\n".join(lines) + "\n"


def config_from_settings(settings: dict) -> tuple:
    """(RunConfig, upload) from a settings dict; upload is (filename, bytes) or None."""
    if not isinstance(settings, dict):
        raise BuilderError("settings() takes a dictionary, as written by Save session")
    d = json.loads(json.dumps(settings))  # deep copy; tuples become lists
    material = dict(d.get("material") or {})
    text = material.pop("file_text", None)
    d["material"] = material
    upload = None
    if text is not None:
        upload = (material.get("filename") or "session_crystal.cif", text.encode("utf-8"))
    cfg = RunConfig.from_dict(d, uploaded_bytes=upload[1] if upload else None)
    return cfg, upload


__all__ = ["session_text", "short_commands", "non_default", "config_from_settings", "call", "fmt",
           "ORDERS", "ORDER_NAMES"]
