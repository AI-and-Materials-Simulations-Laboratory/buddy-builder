# Collaboration handoffs

The first implementation was written in one pass, so there were no
per-participant handoffs. For further work, the spec's ownership (section 14,
with the rebalanced assignments) maps onto the files as follows.

| Participant | Owns | Tests |
|---|---|---|
| 1 Crystallography, materials, geometry | `material_presets.json`, `crystal.py`, `importers.py`, `shapes.py`, `cell.py`, `examples/fixtures/mg_hcp_fully_occupied.cif` | `test_material_presets.py`, `test_importers.py`, `test_shapes.py` |
| 2 Ligands, topology, LAMMPS | `ligands.py`, `chemistry.py`, `topology.py`, `exporter_lammps.py` | `test_ligands.py`, `test_chemistry.py`, `test_topology.py`, `test_export_lammps.py` |
| 3 Sites, placement, QA | `sites.py`, `placement.py`, `validation.py`, `examples/fixtures/invalid_*.cif` | `test_sites.py`, `test_placement.py`, `test_integration.py` |
| 4 GUI, GSD/generic exports, integration, nanoHUB | `gui.py`, `viewer.py`, `exporter_gsd.py`, `exporter_generic.py`, `reporting.py`, `pipeline.py`, `cli.py`, `examples.py`, `app.ipynb`, `middleware/invoke`, `src/Makefile` | `test_export_gsd.py`, `test_gui.py` |

`schema.py` is shared: change it only with the agreement of all four
participants, and run the full suite afterwards.

## Handoff checklist (spec H8 / 16.3)
- `python -m pytest` passes in the owner's environment; paste the summary line
  into `HANDOFF.md`.
- List the modified files and any change to a public function or schema field.
- Note library versions if they differ from `requirements.txt`.
- The integration owner runs the full suite and one GUI example before
  committing to SVN or sealing a numbered ZIP.

## Snapshot log

| Date | SVN revision or ZIP | Tests | Notes |
|---|---|---|---|
| 2026-10-21 | local folder (no SVN yet) | 105 passed, 1 skipped | First complete implementation; see `implementation_log.md` |
