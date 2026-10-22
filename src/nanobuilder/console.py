"""Python console for Buddy Builder: type the build as short Python calls.

Like the Tk console of VMD, but with Python vocabulary::

    material("Au"); shape("sphere", diameter=4)
    ligand("dodecanethiol"); coverage(count=60, order="fibonacci")
    build()
    select(element="S"); measure(10, 25); remove_ligand(3); undo()

Only the commands listed by ``help()`` can run, and their arguments must be
literals (numbers, strings, lists, tuples); input is parsed with :mod:`ast`
and never passed to ``eval``. Every command drives the same widgets as the
GUI, so the tabs always show what the console did, and a Generate click is
echoed here as the equivalent commands.
"""

from __future__ import annotations

import ast
import html
import json

import numpy as np

from .schema import BuilderError
from .session import ORDERS, config_from_settings, short_commands

MAX_LINES = 400

HELP = {
    "help": ("help() or help(\"select\")", "List the commands, or explain one."),
    "material": ("material(\"Au\") | material(\"anatase\") | material(preset=\"Au_fcc\")",
                 "Find a material as in tab 1: presets, the preloaded Materials Project library "
                 "(formula, name or mp id), or carbon words such as graphene, GO, CNT, C60."),
    "carbon": ("carbon(kind=\"graphene\", size=(4, 4), layers=1, oxygen=20, vacancy=0, n=10, m=10, "
               "length=5, walls=1, fullerene=1, shells=0, periodic=False)",
               "Carbon nanomaterial and its size, defects and oxygen content (tabs 1-3)."),
    "shape": ("shape(\"sphere\", diameter=3) | shape(\"rod\", diameter=2, length=6, axis=(0, 0, 1)) | "
              "shape(\"slab\", size=(4, 4), layers=6, miller=(1, 1, 1)) | "
              "shape(\"channel\", type=\"slit\", width=1.5, count=2)",
              "Crystal shape and size (tab 2). Channels: type slit or pore."),
    "box": ("box(\"rectangular\", size=(x, y, z)) | box(\"cube\", edge=8) | box(\"auto\")",
            "Simulation box in nm. Slabs, channels and periodic carbon fill x and y (whole "
            "crystal repeats); finite models must fit. Auto fits the model plus vacuum."),
    "ligand": ("ligand(\"dodecanethiol\") | ligand(smiles=\"SCCO\", anchor=0) | "
               "ligand(template=\"chemisorption\")", "Ligand from the library or a SMILES, and how "
               "it binds: Au_thiolate, chemisorption (bond), geometric (physisorption) or a carbon "
               "template (tab 4)."),
    "coverage": ("coverage(count=40, region=\"entire_surface\", order=\"fibonacci\") | "
                 "coverage(density=3.5)",
                 "How many ligands, where, and in what order: fibonacci, even or random (tab 5)."),
    "export": ("export(\"lammps\", \"pdb\", \"mol2\", \"hoomd\", \"xyz\")",
               "Choose file formats (tab 6); export() alone lists the download files."),
    "build": ("build()", "Generate the model (same as the Generate model button)."),
    "select": ("select(element=\"S\") | select(ligand=3) | select(role=\"ligand\") | "
               "select(index=[1, 2]) | select(within=5, of=12)",
               "Highlight atoms in the viewer; filters combine with AND. select() clears."),
    "measure": ("measure(i, j) | measure(i, j, k) | measure(i, j, k, l)",
                "Distance, angle or dihedral between atoms (periodic boxes use the nearest image)."),
    "info": ("info(i)", "Element, type, role, ligand and position of atom i."),
    "remove_ligand": ("remove_ligand(3) | remove_ligand(3, 7, 9)",
                      "Remove ligands by number (the ligand number shown in the viewer)."),
    "add_ligand": ("add_ligand(atom)", "Place one more ligand on surface atom number atom."),
    "undo": ("undo()", "Undo the last remove / add edit."),
    "summary": ("summary()", "Atoms, ligands, density and validation of the current model."),
    "examples": ("examples() | example(\"B\")", "List the bundled examples, or load one."),
    "theme": ("theme(\"buddy\" | \"light\" | \"dark\" | \"auto\")", "Switch the colour theme."),
    "mode": ("mode(\"basic\" | \"advanced\")", "Show only the main settings, or all of them."),
    "config": ("config()", "Print the current settings as JSON (what config.json will hold)."),
    "settings": ("settings({...})", "Restore saved selections (the dictionary in a Save session "
                 "file); follow with build()."),
    "session": ("session()", "Print the current session: the commands that rebuild this model."),
    "cls": ("cls()", "Clear this console."),
    "reset": ("reset()", "Back to the launch settings (Au sphere, 12 SCCO), like the Reset button."),
    "ask": ("ask(\"a 4 nm gold sphere with 60 dodecanethiol\")",
            "Ask Buddy in plain English: plans the commands, checks them, and shows the plan "
            "in the Ask Buddy panel to apply."),
}

FORMATS = {"hoomd": "t_hoomd", "gsd": "t_hoomd", "lammps": "t_lammps", "xyz": "t_generic",
           "extxyz": "t_generic", "pdb": "t_pdb", "mol2": "t_mol2"}


class Console:
    """Command interpreter bound to a :class:`nanobuilder.gui.BuilderApp`."""

    def __init__(self, app):
        self.app = app
        self.lines = []
        self.selection = []
        self.commands = {name: getattr(self, "cmd_" + name) for name in HELP}
        self.commands.update({"generate": self.cmd_build, "remove_ligands": self.cmd_remove_ligand,
                              "example": self.cmd_examples})
        self.log("Buddy Builder console: type help() for the commands.", kind="echo")

    # ------------------------------------------------------------------ output
    def log(self, text: str, kind: str = "out", note: str = None):
        if note:
            text = f"{text}    # {note}"
        for line in str(text).splitlines() or [""]:
            self.lines.append((kind, line))
        self.lines = self.lines[-MAX_LINES:]
        body = "".join(f"<pre class='bb-{k}'>{html.escape(t)}</pre>" for k, t in reversed(self.lines))
        self.app.console_log.value = f"<div class='bb-log'>{body}</div>"

    # ------------------------------------------------------------------ parsing
    def run(self, text: str):
        """Run one line (several commands may be separated by ';' or newlines)."""
        self.log(f">>> {text}", kind="in")
        try:
            tree = ast.parse(text, mode="exec")
        except SyntaxError as exc:
            self.log(f"SyntaxError: {exc.msg}", kind="err")
            if len(text.split()) >= 3 and getattr(self.app, "assistant", None) is not None:
                self.log("# Plain English? Press AI to turn it into commands, or ask Buddy in the "
                         "chat above", kind="echo")
            return
        for stmt in tree.body:
            try:
                out = self._exec(stmt)
            except BuilderError as exc:
                self.log(f"Error: {exc}", kind="err")
                return
            except (TypeError, ValueError, KeyError, IndexError) as exc:
                self.log(f"{type(exc).__name__}: {exc}", kind="err")
                return
            if out:
                self.log(out)

    def _exec(self, stmt):
        if not isinstance(stmt, ast.Expr):
            raise BuilderError("Only commands are allowed, e.g. build() or select(element=\"S\"); "
                               "type help()")
        node = stmt.value
        if isinstance(node, ast.Name):  # bare `help` or `build`
            node = ast.Call(func=node, args=[], keywords=[])
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)):
            raise BuilderError("Commands look like name(arguments); type help()")
        name = node.func.id
        if name not in self.commands:
            raise BuilderError(f"Unknown command {name}(); type help() for the list")
        try:
            args = [ast.literal_eval(a) for a in node.args]
            kwargs = {k.arg: ast.literal_eval(k.value) for k in node.keywords}
        except ValueError:
            raise BuilderError("Arguments must be plain values: numbers, \"text\", lists or tuples")
        return self.commands[name](*args, **kwargs)

    # ------------------------------------------------------------------ helpers
    def _assembly(self):
        if self.app.result is None:
            raise BuilderError("No model yet: build() first")
        return self.app.result.assembly

    def _atom(self, i) -> int:
        a = self._assembly()
        i = int(i)
        if not 0 <= i < a.n_atoms:
            raise BuilderError(f"Atom {i} does not exist (0 to {a.n_atoms - 1})")
        return i

    def _vectors(self, _=None):
        a = self._assembly()
        pos = np.asarray(a.positions_A, dtype=float)
        cell, pbc = np.asarray(a.cell_A, dtype=float), tuple(bool(x) for x in a.pbc)

        def vec(i, j):
            d = pos[j] - pos[i]
            if any(pbc):
                from ase.geometry import find_mic
                d = find_mic(d[None, :], cell, pbc)[0][0]
            return d
        return vec

    def _set(self, widget, value):
        widget.value = value

    # ------------------------------------------------------------------ commands
    def cmd_help(self, name=None):
        if name:
            if name not in HELP:
                raise BuilderError(f"No command {name}; type help()")
            usage, text = HELP[name]
            return f"{usage}\n    {text}"
        width = max(len(k) for k in HELP)
        return "\n".join(f"{k.ljust(width)}  {v[1]}" for k, v in HELP.items()) + \
            "\nhelp(\"name\") shows the call forms. Commands also echo here when you use the tabs."

    def cmd_material(self, name=None, preset=None):
        app = self.app
        if preset:
            app.source.value = "preset"
            app.preset.value = preset
            return f"Material: {app.crystal.structure_label if app.crystal else preset}"
        if not name:
            raise BuilderError("Give a material, e.g. material(\"Au\") or material(preset=\"Au_fcc\")")
        app._quiet = True
        try:
            app.mat_search.value = str(name)
        finally:
            app._quiet = False
        app._on_find()
        if app._is_carbon():
            return f"Material: carbon ({app.carbon.kind.value})"
        return f"Material: {app.crystal.structure_label if app.crystal else name}"

    def cmd_carbon(self, kind=None, size=None, layers=None, oxygen=None, vacancy=None, n=None, m=None,
                   length=None, walls=None, fullerene=None, shells=None, periodic=None, profile=None):
        app, c = self.app, self.app.carbon
        app.source.value = "carbon"
        kinds = {"graphene": "graphene", "go": "graphene", "sheet": "graphene", "nanotube": "nanotube",
                 "cnt": "nanotube", "fullerene": "carbon_sphere", "onion": "carbon_sphere",
                 "carbon_sphere": "carbon_sphere", "graphene oxide": "graphene",
                 "graphene_oxide": "graphene", "rgo": "graphene", "flake": "graphene",
                 "carbon nanotube": "nanotube", "swcnt": "nanotube", "mwcnt": "nanotube",
                 "c60": "carbon_sphere", "carbon onion": "carbon_sphere"}
        if kind:
            if kind.lower() not in kinds:
                raise BuilderError(f"carbon kind is one of {', '.join(sorted(set(kinds)))}")
            c.kind.value = kinds[kind.lower()]
            if kind.lower() in ("onion", "carbon onion"):
                c.sphere_type.value = "onion"
            elif kind.lower() in ("fullerene", "c60"):
                c.sphere_type.value = "fullerene"
        if size is not None:
            c.length.value, c.width.value = (float(x) for x in size)
        for widget, value in ((c.layers, layers), (c.oxygen, oxygen), (c.vacancy, vacancy),
                              (c.cnt_n, n), (c.cnt_m, m), (c.cnt_len, length), (c.walls, walls),
                              (c.fullerene, fullerene), (c.shells, shells), (c.profile, profile)):
            if value is not None:
                widget.value = value
        if periodic is not None:
            c.periodic.value = bool(periodic)
        return f"Carbon: {c.kind.value}" + (c.geom_info.value.replace("<small>", " - ")
                                            .replace("</small>", "") if c.geom_info.value else "")

    def cmd_shape(self, kind, diameter=None, length=None, axis=None, caps=None, size=None, layers=None,
                  miller=None, termination=None, type=None, width=None, count=None):
        app = self.app
        if app._is_carbon():
            raise BuilderError("Carbon materials set their size with carbon(...)")
        app.shape.value = kind
        if diameter is not None:
            app.diameter.value = float(diameter)
        if length is not None:
            app.rod_len.value = float(length)
        if axis is not None:
            app.rod_u.value, app.rod_v.value, app.rod_w.value = (int(x) for x in axis)
        if caps is not None:
            app.rod_caps.value = {"round": "hemispherical", "flat": "flat"}.get(caps, caps)
        if size is not None:
            app.slab_len.value, app.slab_wid.value = (float(x) for x in size)
        if layers is not None:
            app.thick_mode.value = "layers"
            app.layers.value = int(layers)
        if miller is not None:
            app.h.value, app.k.value, app.l.value = (int(x) for x in miller)
        if termination is not None:
            app.termination.value = int(termination)
        if type is not None:
            app.channel_type.value = type
        if width is not None:
            app.channel_width.value = float(width)
        if count is not None:
            app.channel_count.value = int(count)
        return f"Shape: {kind}"

    def cmd_box(self, shape="rectangular", size=None, edge=None):
        app = self.app
        if shape not in ("auto", "rectangular", "cube"):
            raise BuilderError('box shape is "auto", "rectangular" or "cube"')
        app.box_shape.value = shape
        if shape == "cube" and (edge is not None or size is not None):
            app.box_edge.value = float(edge if edge is not None else max(size))
        elif size is not None:
            app.box_x.value, app.box_y.value, app.box_z.value = (float(x) for x in size)
        g = app.config().geometry
        return f"Box: {shape}" + (f", {tuple(g.box_lengths_nm)} nm" if g.box_lengths_nm else "")

    def cmd_ligand(self, name=None, smiles=None, anchor=None, template=None):
        app = self.app
        if name is not None:
            keys = [k for _, k in app.library.options if k]
            match = [k for k in keys if k == name] or [
                k for lbl, k in app.library.options if k and str(name).lower() in lbl.lower()]
            if not match:
                raise BuilderError(f"No library ligand {name!r} for this material; "
                                   "use ligand(smiles=...)")
            app.library.value = match[0]
        if smiles is not None:
            app.smiles.value = smiles
        if anchor is not None:
            app.anchor.value = int(anchor)
        if template is not None:
            app.template.value = template
        return f"Ligand: {app.smiles.value} via {app.template.value}"

    def cmd_coverage(self, count=None, density=None, region=None, order=None):
        app = self.app
        if count is not None:
            app.cov_mode.value = "count"
            app.count.value = int(count)
        if density is not None:
            app.cov_mode.value = "density"
            app.density.value = float(density)
        if region is not None:
            app.region.value = region
        if order is not None:
            if order not in ORDERS:
                raise BuilderError("order is fibonacci, even or random")
            app.distribution.value = ORDERS[order]
        what = (f"{app.count.value} ligands" if app.cov_mode.value == "count"
                else f"{app.density.value}/nm^2")
        return f"Coverage: {what} on {app.region.value}, {app.distribution.value}"

    def cmd_export(self, *formats):
        app = self.app
        if formats:
            for key, attr in FORMATS.items():
                getattr(app, attr).value = False
            for f in formats:
                if f.lower() not in FORMATS:
                    raise BuilderError(f"Formats: {', '.join(sorted(FORMATS))}")
                getattr(app, FORMATS[f.lower()]).value = True
            return "Formats: " + ", ".join(formats) + ("; build() to make the files"
                                                      if app.result is None or app._stale else
                                                      "; build() again to include them")
        if app.package is None:
            return "No files yet: build() first"
        files = sorted(p.relative_to(app.package.stage_dir).as_posix()
                       for p in app.package.stage_dir.rglob("*") if p.is_file())
        return "Download files (links under Downloads):\n  " + "\n  ".join(files)

    def cmd_build(self):
        app = self.app
        app._echo = False
        try:
            app._generate()
        finally:
            app._echo = True
        r = app.result
        if r is None or getattr(app, "last_failure", None) and app.progress_label.value == "failed":
            from .gui import failure_advice
            msg = getattr(app, "last_failure", "") or "see Status"
            why, fix = failure_advice(msg)
            return f"Build failed: {msg}\n  why: {why}\n  fix: {fix}"
        a = r.assembly
        return (f"Built {a.n_atoms} atoms, {a.achieved_ligands}/{a.target_ligands} ligands, "
                f"validation {r.report.status}")

    def cmd_select(self, element=None, ligand=None, role=None, index=None, within=None, of=None):
        a = self._assembly()
        n = a.n_atoms
        mask = np.ones(n, dtype=bool)
        if element is not None:
            els = {element} if isinstance(element, str) else set(element)
            mask &= np.isin(np.asarray(a.elements), list(els))
        if ligand is not None:
            ids = {int(ligand)} if isinstance(ligand, (int, float)) else {int(x) for x in ligand}
            mask &= np.isin(np.asarray(a.molecule_id), list(ids))
        if role is not None:
            roles = np.asarray(a.chemical_role)
            mask &= (roles != "substrate") if role == "ligand" else (roles == role)
        if index is not None:
            sel = np.zeros(n, dtype=bool)
            sel[[self._atom(i) for i in (index if isinstance(index, (list, tuple)) else [index])]] = True
            mask &= sel
        if within is not None:
            if of is None:
                raise BuilderError("select(within=5, of=12): atoms within 5 A of atom 12")
            vec = self._vectors(None)
            centre = self._atom(of)
            d = np.array([np.linalg.norm(vec(centre, j)) for j in range(n)])
            mask &= d <= float(within)
        cleared = all(v is None for v in (element, ligand, role, index, within))
        if cleared:
            mask[:] = False
        self.selection = [int(i) for i in np.nonzero(mask)[0]]
        self.app._show_assembly(selection=self.selection)
        if not self.selection:
            return "Selection cleared" if cleared else "No atoms match"
        head = ", ".join(str(i) for i in self.selection[:20])
        more = f", ... ({len(self.selection)} atoms)" if len(self.selection) > 20 else ""
        return f"{len(self.selection)} atoms selected (cyan): {head}{more}"

    def cmd_measure(self, i, j, k=None, l=None):  # noqa: E741
        idx = [self._atom(x) for x in (i, j, k, l) if x is not None]
        vec = self._vectors(None)
        if len(idx) == 2:
            return f"distance {idx[0]}-{idx[1]}: {np.linalg.norm(vec(idx[0], idx[1])):.3f} A"
        if len(idx) == 3:
            u, v = vec(idx[1], idx[0]), vec(idx[1], idx[2])
            c = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            return f"angle {idx[0]}-{idx[1]}-{idx[2]}: {np.degrees(np.arccos(np.clip(c, -1, 1))):.2f} deg"
        b0, b1, b2 = vec(idx[1], idx[0]), vec(idx[1], idx[2]), vec(idx[2], idx[3])
        n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
        m = np.cross(n1, b1 / np.linalg.norm(b1))
        ang = np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2)))
        return f"dihedral {'-'.join(str(x) for x in idx)}: {ang:.2f} deg"

    def cmd_info(self, i):
        a = self._assembly()
        i = self._atom(i)
        x, y, z = a.positions_A[i]
        lig = int(a.molecule_id[i])
        return (f"atom {i}: {a.elements[i]} ({a.type_labels[i]}), role {a.chemical_role[i]}"
                + (f", ligand {lig}" if a.chemical_role[i] != "substrate" else "")
                + f", at ({x:.2f}, {y:.2f}, {z:.2f}) A")

    def cmd_remove_ligand(self, *ids):
        if len(ids) == 1 and isinstance(ids[0], (list, tuple)):
            ids = tuple(ids[0])
        if not ids:
            raise BuilderError("remove_ligand(3) or remove_ligand(3, 7)")
        r = self.app.apply_edit("remove", [int(x) for x in ids])
        return f"{r.assembly.achieved_ligands} ligands now (ligands are renumbered 1..N)"

    def cmd_add_ligand(self, atom):
        r = self.app.apply_edit("add", int(atom))
        return f"Added ligand {r.assembly.achieved_ligands} on atom {atom}"

    def cmd_undo(self):
        if not self.app.history:
            return "Nothing to undo"
        self.app._set_result(self.app.history.pop())
        return f"Undone; {self.app.result.assembly.achieved_ligands} ligands"

    def cmd_summary(self):
        a = self._assembly()
        r = self.app.result.report
        edits = len(a.provenance.get("manual_edits", []))
        return (f"{a.n_atoms} atoms ({a.n_substrate} substrate), {a.achieved_ligands} ligands, "
                f"{a.achieved_density_nm2:.3f}/nm^2 on {a.surface_area_nm2:.2f} nm^2; "
                f"validation {r.status}" + (f"; {edits} hand edits" if edits else ""))

    def cmd_examples(self, key=None):
        from .examples import EXAMPLES
        if key is None:
            return "\n".join(f"{k[0]}  {e['label']}" for k, e in EXAMPLES.items())
        match = [k for k in EXAMPLES if k == key or k[0].lower() == str(key).lower()[:1]]
        if not match:
            raise BuilderError("example(\"B\"): letters A to " + list(EXAMPLES)[-1][0])
        self.app.example.value = match[0]
        return f"Loaded {EXAMPLES[match[0]]['label']}; build() to generate"

    def cmd_theme(self, name):
        self.app.theme.value = name
        return f"Theme: {name}"

    def cmd_mode(self, name):
        self.app.mode.value = name
        return f"Settings: {name}"

    def cmd_config(self):
        d = self.app.config().to_dict()
        d.get("material", {}).pop("uploaded_bytes", None)
        return json.dumps(d, indent=1, default=str)

    def cmd_settings(self, settings):
        cfg, upload = config_from_settings(settings)
        self.app.load_config(cfg, upload)
        return "Settings restored; build() to generate"

    def cmd_session(self):
        from .session import session_text
        r = self.app.result
        return session_text(self.app.config(), r.assembly if r is not None else None).rstrip()

    def cmd_reset(self):
        self.app.reset()
        return "Settings reset to the launch defaults"

    def cmd_ask(self, request=None):
        panel = getattr(self.app, "assistant", None)
        if panel is None:
            raise BuilderError("Ask Buddy is not available in this session")
        if not request:
            raise BuilderError('ask("a 4 nm gold sphere with 60 dodecanethiol")')
        plan = panel.ask(str(request))
        return plan.console_summary()

    def cmd_cls(self):
        self.lines = []
        self.log("", kind="echo")
        return ""

    # ------------------------------------------------------------------ echo
    def echo_build(self, cfg, result):
        """Write the commands equivalent to the GUI settings of a build."""
        if not getattr(self.app, "_echo", True):
            return
        out = short_commands(cfg, with_ligand=result.ligand is not None)
        a = result.assembly
        out.append(f"build()    # {a.n_atoms} atoms, {a.achieved_ligands}/{a.target_ligands} ligands, "
                   f"{result.report.status}")
        self.log("# Generate model, as commands:", kind="echo")
        for line in out:
            self.log(line, kind="echo")


__all__ = ["Console", "HELP"]
