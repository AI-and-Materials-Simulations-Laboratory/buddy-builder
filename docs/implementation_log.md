# Implementation log

What was built from `nanohub_crystal_ligand_builder_workflow_no_github.md`
(the "spec"), where the implementation differs from it and why, and what was
and was not verified. Written 2026-10-21. The spec document itself was not
changed during implementation.

## 1. Steps taken

1. **Environment.** Created `xtalligand/.venv` with Python 3.12.10 on Windows 11
   and installed numpy 2.5.4, scipy 1.18.1, ase 3.29.0, rdkit 2026.9.1,
   gsd 5.0.1, ipywidgets 8.1.9, notebook 7.6.3 and pytest 9.1.1 (pinned in
   `requirements.txt`). Also tried the PyPI `lammps` package for real
   `read_data` checks; its DLL needs MS-MPI, which is not installed here, so it
   was uninstalled and the LAMMPS check is an opt-in test (`LAMMPS_EXE`).
2. **Schema first (spec section 5).** `schema.py` defines the run configuration
   (`MaterialConfig`, `GeometryConfig`, `LigandConfig`, `CoverageConfig`,
   `ExportConfig`, `RunConfig`) and every model passed between stages
   (`CrystalModel`, `SubstrateModel`, `LigandModel`, `SiteCatalog`,
   `AssemblyModel`, `ValidationReport`, `ExportReport`), using the field names
   of the spec's tables.
3. **Stage A, materials.** `material_presets.json` holds eight registry entries
   with lattice constants and their sources; `crystal.py` builds conventional
   cells (`bulk(..., cubic=True)` for cubic prototypes). `importers.py` reads
   uploads with every check from spec A2: size and atom caps, sanitized names,
   `index=":"` to catch multi-block CIFs, occupancy check via
   `atoms.info["occupancy"]`, declared-formula check, minimum-distance check.
4. **Stage B, geometry.** `shapes.py` cuts atom-centred spheres and builds
   slabs with ASE `surface()`, enumerates distinct terminations, and finds
   surface atoms by coordination number. `cell.py` handles the lower-triangular
   cell, tilt reduction and the LAMMPS/HOOMD box conventions.
5. **Stage C, ligands.** `ligands.py` (RDKit ETKDGv3 + MMFF94, UFF fallback)
   and `chemistry.py` (attachment templates and the single Au-thiolate
   availability rule).
6. **Stages D-E.** `sites.py` (on-top sites with an accessibility probe),
   `placement.py` (coverage targets, seeded random or greedy farthest-point site
   order, 12 azimuth trials per conformer, KD-tree clash checks including
   in-plane periodic images), `topology.py` (canonical model, bonds, angles,
   dihedrals, structural type labels, parameter-coverage report) and
   `validation.py`.
7. **Stage F, exports.** `exporter_gsd.py` (GSD Python API), `exporter_lammps.py`
   (`atom_style full`, one-based IDs, template input script),
   `exporter_generic.py` (extXYZ + `topology.json`) and `reporting.py` (ZIP,
   manifest with SHA-256, README files). Every export is read back and
   compared before the ZIP is written; a failed round trip blocks the export.
8. **Stage G, GUI.** `gui.py` (five tabs, results panel, atom inspector,
   status, data-URI download links, example loader, reset) and `viewer.py`
   (3Dmol.js page in an iframe). `pipeline.py` gives the GUI one call per action.
9. **Stage H, packaging.** `app.ipynb`, `middleware/invoke`, `src/Makefile`,
   `pyproject.toml`, `requirements.txt`, the CLI (`cli.py`), example configs and
   fixtures.
10. **Tests and fixes.** Wrote 105 tests (one more is opt-in) and fixed what they
    and the browser check found (section 3).

## 2. Changes and decisions relative to the spec

### Added beyond the spec's file list
- `pipeline.py`: validate / run / export entry points shared by the GUI, CLI
  and tests, so the GUI holds no scientific logic.
- `cli.py`: builds a ZIP from a config JSON without a notebook; a `config.json`
  from any export reproduces the run.
- `examples.py`: the five live examples as data. `examples/*.json` hold the
  same configs and a test keeps them in sync.
- P1 items done because they were cheap: Mg/HCP, Ge/diamond and NaCl/rocksalt
  presets; VASP POSCAR upload; multiple conformers used as placement retries;
  evenly spaced (greedy) distribution; termination selection.
- `save_crystal_fixture` / `load_crystal_fixture` in `schema.py` (extXYZ + JSON
  sidecar) for team handoffs, as spec 5.2 asks.

### Different from the spec, with reasons
| Topic | Spec | Implemented | Why |
|---|---|---|---|
| 3D viewer | py3Dmol first, Plotly fallback | 3Dmol.js page inside an `<iframe srcdoc>` in an `ipywidgets.HTML` | Runs identically in classic Notebook, AppMode, JupyterLab and Notebook 7 regardless of how each treats scripts in widget output, and needs no extra Python package. The browser must reach cdnjs. |
| Viewer controls | Python widgets | Style, show/hide ligands, site highlight, cell and reset are buttons inside the viewer page; atom hover shows details; a Python atom-index inspector is also provided | Changing these does not reset the camera or call the kernel. |
| Thiolate edit order | Apply the attachment edit before conformer generation | Embed and MMFF-optimize the intact thiol, then remove the S-H hydrogen | MMFF has no parameters for the radical left by removing H; the remaining geometry is a sensible starting point. Index maps are kept either way. |
| "Layers" | Not defined | One layer = one atomic plane | ASE's `layers` counts repeat units, so "4 layers" gave 8 planes for (100) but 4 for (111). Found by a failing test. |
| Slab size | At least the requested size | At least the requested size, with a 1% tolerance | 10 units of 4.997 A (49.97 A) would otherwise round up to 11 for a 5 nm request. |
| Coverage area | Sphere 4 pi R^2; slab cell face area | Same | |
| Box | Vacuum padding | Vacuum is the clearance from the outermost atom, **ligands included**, along every nonperiodic axis; sphere boxes are orthogonal and fitted per axis | Guarantees ligands never stick out of the box. |
| Surface atoms | Reduced coordination | Coordination below the bulk value for that species, with per-species-pair cutoffs halfway between the first and second neighbour shells | Works for FCC, BCC, diamond, HCP and NaCl without one global Au cutoff (spec D1). |
| Buried sites | "Reject buried or inaccessible sites" | Rejected when another substrate atom is closer to the anchor probe point than 0.9 x the anchor distance | Simple geometric test; fcc(110) trough atoms count as accessible, hollows on the faceted sphere do not. |
| Clash rule | r < alpha (vdW_i + vdW_j) | Same; alpha defaults to 0.7; Bondi radii from ASE, covalent radius + 0.8 A where Bondi has none (e.g. Fe) | |
| Geometric adsorption distance | Not specified | Default 2.5 A (Au-thiolate stays 2.4 A) | Editable in the GUI. |
| Molecule IDs | Substrate "uniquely identified" | Substrate 0, ligands 1..n, identical in every file (LAMMPS accepts 0) | |
| Type labels | No element-as-force-field types | `atom_types` stays `None`; files use structural labels `<El>_sub`, `<El>_anchor`, `<El>` | Engines need some type names; the labels are documented as non-force-field everywhere. |
| Impropers | Optional | Not enumerated | Their definition depends on the force field. |
| GSD molecule ID | Not specified | Stored as `log/particles/molecule_id` | GSD has no molecule field. |
| LAMMPS types | Not specified | Type names as comments in `Masses`; no `Type Labels` sections | Older LAMMPS versions reject the newer sections. |
| LAMMPS box | Valid triclinic | Tilt reduced so abs(xy) <= lx/2 | LAMMPS rejects larger tilts. |
| Downloads | Download buttons | Data-URI links (ZIP, GSD, LAMMPS data, extXYZ, validation report), shown only for an up-to-date build | No server file route needed; files over 40 MB are not linked. |
| Uploaded cell | "Normalized" cell | Used as given (positions wrapped); no symmetry standardization; space group not computed | pymatgen not added. `normalized_unit_cell.cif` is the cell actually used. |

### Spec options not implemented
- Tab 3 "Custom tilt" orientation and "Automatic" adsorption site (on-top only).
- Finite (non-periodic) slab / rectangular sheet and other P1/P2 shapes.
- Optional PDB export.
- "Simulation-ready" fidelity, imported charges and "no explicit hydrogens":
  the widgets are present but fixed or rejected with a message, as the spec
  allows only validated templates for these, and none exist.

## 3. Problems found and fixed during verification
- **Slab layer count meant different things for different planes** (above).
- **Shared ipywidgets `Layout`**: hiding the HCP-only `c` field hid every
  field that shared the same Layout object. Each widget now gets its own
  Layout; a test checks this.
- **Viewer zoom** included the vacuum-padded box, so particles looked tiny; it
  now zooms to the atoms, and the box is hidden by default for spheres.
- **Site highlight** used `addStyle`, which kept the element colour scheme;
  changed to `setStyle`.
- **Repeated exports to one folder** could mix files from earlier runs; each
  export now writes to a clean folder named after the ZIP.
- **Windows setup with Anaconda installed** (reported by a user): `python` was
  Anaconda's 3.11, so re-running `python -m venv .venv` re-pointed the existing
  3.12 environment at 3.11 without reinstalling compiled packages, and the
  notebook ran in a kernel without the dependencies (`No module named 'ase'`).
  The README now requires Python 3.12, creates the environment with an explicit
  interpreter, calls `.venv\Scripts\python` directly and has a troubleshooting
  table; `app.ipynb` now names the Python its kernel is using when imports
  fail. `.venv` was rebuilt with Python 3.12.10.
- **`python -m jupyter notebook` started Anaconda's Jupyter** even with the
  venv's Python: `jupyter` dispatches subcommands by searching `PATH`, and
  Anaconda's `Scripts` folder came first (Notebook 7.2.2 instead of the venv's
  7.6.3). The documented command is now `.venv\Scripts\python -m notebook`,
  which runs the venv's own Notebook and therefore its kernel.

## 4. Verification

| Check | Result |
|---|---|
| `python -m pytest` | 105 passed, 1 skipped (real LAMMPS `read_data`, needs `LAMMPS_EXE`) |
| Five live examples (A-E) end to end | All `ok`; see the table below |
| GUI in a real browser (local Jupyter Notebook 7.6.3, in-app browser) | Cell runs; GUI renders; Generate works for examples B, C and E; viewer rotates (drag), zooms (scroll), hides ligands, highlights the 12 attachment sites, shows atom details on hover; download links present |
| Round trips | GSD reopened with gsd 5.0.1; LAMMPS data parsed back by `read_lammps_data`; extXYZ read back with ASE; manifest hashes match the ZIP contents |

| Example | Atoms (substrate) | Ligands placed / requested | Area (nm^2) | Achieved density (/nm^2) |
|---|---|---|---|---|
| A bare Au sphere 3 nm | 887 (887) | 0 / 0 | 28.274 | - |
| B Au sphere + 12 SCCO | 995 (887) | 12 / 12 | 28.274 | 0.424 |
| C Au(111) slab, 0.5/nm^2 top | 1726 (1600) | 14 / 14 | 28.832 | 0.486 |
| D Fe/BCC sphere 2.5 nm | 701 (701) | 0 / 0 | 19.635 | - |
| E uploaded Mg HCP CIF, (001) slab + ethanol | 417 (363) | 6 / 6 | 10.798 | 0.556 |

### Not verified
- **nanoHUB**: no tool registration, SVN checkout, `make install`, AppMode
  launch or hosted package check was possible from this machine. See
  `deployment_notes.md` for the remaining steps.
- **Real MD engines**: neither LAMMPS nor HOOMD-blue was run on the files
  (LAMMPS needs MS-MPI on Windows; HOOMD-blue has no Windows build). Set
  `LAMMPS_EXE` on a machine with LAMMPS to enable the `read_data` test.
- **ipywidgets 7** (classic nanoHUB environments may use it): the upload
  adapter handles both value formats, but the GUI was only run with 8.1.9.
- **Saving the downloaded files** through the browser's download dialog; the
  links and their contents were checked, not the saved files.
- Viewer behaviour with no internet in the browser (an error message is shown
  instead of the viewer).

## 5. Second iteration (2026-10-21, afternoon)

Requested by the team after the first demo; each item has tests.

| Item | What changed | Where |
|---|---|---|
| Ligands split across slab boundaries | Atoms were wrapped one by one, so a ligand at the cell edge was cut in two and its bonds drawn across the whole slab (reproduced: 12 of 58 C11 thiols split, 58 A "bonds"). Each ligand is now rebuilt whole from its surface atom via minimum-image bond vectors. LAMMPS writes wrapped coordinates + `ix iy iz`; GSD writes `particles/image`. | `cell.make_molecules_whole`, `cell.wrap_with_images`, `topology.assemble`, exporters; `tests/test_pbc_whole_molecules.py` |
| Surface anchor species | Every element of a multi-element crystal is ticked by default; a warning shows when none is, and Generate reports it. `surface_species = None` now means "all", `[]` is an error. | `gui._rebuild_species`, `chemistry.resolve_surface_species` |
| PDB export | CRYST1, substrate chain A (one residue per atom), ligands chain B residue `LIG`, unique atom names, CONECT for all bonds, round-trip check. The user chose PDB as the AMBER/GROMACS route instead of engine-specific topology files. | `exporter_pdb.py`; `tests/test_export_pdb.py` |
| Ligand library | 18 common capping agents in 5 groups with anchor, template and typical distance; selecting one fills the ligand tab. New "extended chains (all-trans)" conformer mode. | `ligand_library.py`, `ligands.extend_chains`; `tests/test_ligands.py` |
| Nanorods | New shape: diameter, tip-to-tip length, lattice axis, hemispherical or flat caps; "side" and "tips" placement regions; viewer starts side-on. | `shapes._build_rod`, `sites.rod_normal` |
| Materials Project | Live search (formula, chemical system, mp-id) with the stdlib (no `mp-api`/pymatgen); structures converted to the conventional cell with spglib and stored as CIF in the config and ZIP. Only ordered MP entries (the team chose "only what MP has" for alloys). | `materials_project.py`, `crystal._build_materials_project` |
| Material picker | Type-to-search box with grouped suggestions (built-in elements, metals, alloys, oxides, semiconductors, nanorods); catalog MP entries are matched by formula + space group, not hard-coded IDs; unknown text goes to MP. | `material_catalog.py`, GUI |
| Look | Scoped stylesheet re-pointing Jupyter's widget variables; Auto/Light/Dark switch (Auto follows the browser); controls left, viewer stage + results + downloads right; themed viewer page; long caveat lists collapsed. | `theme.py`, `viewer.py`, `gui.py` |
| Size limits | Benchmarked (`tools/benchmark_sizes.py`, results in `size_limits.md`). Fixed three scaling problems found by it: a quadratic symbol lookup in the sphere cut, a clash KD-tree rebuilt after every ligand (now every 32, newest checked directly), and an all-pairs neighbour list in validation (now KD-tree searches of only the pairs that matter). The sphere cut, random placement and validation give identical output to before (tested). **"Evenly spaced" placement does not**: exact lattice ties between candidate sites used to be decided by floating-point noise; they are now rounded and go to the first site in the seeded order, so the same seed can pick different sites (and place a few more or fewer ligands) than before. Default `max_atoms` raised from 20,000 to 100,000 and `max_attempts` from 5,000 to 50,000; the latter changes results only for configs that used to hit the old attempt limit. | `shapes.py`, `placement.py`, `validation.py` |
| Independent review fixes | A separate review found: flat-capped rods along non-[001] axes lost their rim corners (supercell too small) and could get zero "tips" atoms (classification measured from the nominal cut); both fixed (supercell sized to the rod's farthest point; an atom is a tip when it is closer to its end's outermost plane than to the wall). spglib tolerance raised from 1e-3 to 0.1 A (relaxed MP cells with mA noise were misclassified, e.g. NaCl as Cm), a warning when spglib disagrees with MP's space group, and no more process-wide spglib setting. MP searches are sorted by stability on the server, with a space-group query as fallback for catalog polymorphs. Extended chains are skipped (with a warning) when no force field can relax them. Rod parameters recorded in provenance. | `shapes.py`, `materials_project.py`, `ligands.py`, `gui.py` |

### Not verified in this iteration
- **Live Materials Project calls**: the development machine could not reach
  `api.materialsproject.org`. The client follows the official `mp-api`
  client's HTTP contract (endpoint, `x-api-key` header, `_fields`, `_limit`,
  `formula` / `chemsys` / `material_ids`, `{"data": [...]}` responses,
  pymatgen structure dicts) and is tested against responses in that format.
  First check with a real key: type `AuCu3`, confirm the result list and the
  4-atom cubic cell, then build a 3 nm sphere.
- **Viewer frame rate on a real GPU** for large models (benchmarks used
  software rendering).
- nanoHUB items listed in section 4 are still open.

## 6. Carbon nanomaterials (2026-10-22)

Request: add graphene, graphene oxide, carbon nanotubes and carbon spheres with
covalently bonded ligands, using the HierGO workflow (Garcia et al. 2023, in
`HierGO-master/`) and its AMBER adaptation (Fedai et al., Data 2026) as the
reference for sheets, oxidation and defects such as holes. Reference:
[carbon_materials.md](carbon_materials.md).

| Item | What changed | Where |
|---|---|---|
| Carbon substrate engine | Explicit bond-graph network; graphene (layers, AB/AA, finite or periodic), SWCNT/MWCNT (ASE, matched walls), Goldberg fullerenes C60·k² and onions; edge trimming and H termination. | `carbon_lattice.py`, `carbon.py` |
| HierGO defects | Holes by random growth (size classes, custom range) or HierGO hole-library shapes (read at run time, GPL-3.0, not bundled); pores; Stone-Wales; 24-atom rotations (HierGO mini-GB); ruptures with ether O; hanging-carbon trimming generalized to small rings; local relaxation after topological defects. | `carbon_defects.py` |
| Oxidation | Fresh (hydroxyl:epoxy, 2:1 default, paired anti diols) and HierGO aged 600-2500 K profiles with edge phenol, carbonyl, COOH and ether; custom mixes; uniform, edges-first, domain and tile patterns; face balance; chemical bond lengths and sp3 puckering; sublattice bias and epoxide-to-diol densification for dense GO (0-68 % O). | `carbon_oxidation.py` |
| Metrics | O % (100·N_O/N_C), C:O, modified carbon %, achieved vacancy %, basal OH:epoxide ratio, group counts, every fallback; in the GUI results, README and manifest. | `carbon.py`, `reporting.py`, `gui.py` |
| Covalent ligands | Templates `C_basal`, `C_edge`, `C_carboxyl`, `C_epoxide` with per-site substrate edits (sp2→sp3, H or OH removal, ring opening), template bond lengths, outward axis from the removed anchor H, 1-4 clash exemption, site claims. Library group "Covalent on carbon" (10 entries). | `chemistry.py`, `ligands.py`, `placement.py`, `topology.py`, `ligand_library.py` |
| Typing and MOL2 | Suggested GAFF2 (exact for substrate roles, rule-based for ligands) and HierGO OPLS-AA names; Tripos MOL2 with GO residues (GRA/OH/EPX/ETH/CBO/CBX or GRO), bond orders and CRYSIN; tleap sketch. | `atomtypes.py`, `exporter_mol2.py` |
| Exports | Full-graph angles/dihedrals for covalent substrates; PDB residues/atom names; extXYZ residue column; LAMMPS boundary notes (`p p f`, `f f p`); covalent validation (valence, dangling carbons, overlaps). | exporters, `validation.py` |
| GUI | "Carbon nanomaterial" source; carbon geometry panel; new tab "3 Defects & GO"; carbon coverage regions (top/bottom/both/edges, side/tips); MOL2 target and residue scheme; 14 carbon catalog entries; examples F-J. | `gui_carbon.py`, `gui.py`, `material_catalog.py`, `examples.py` |
| Viewer | Bond-list rendering for covalent frameworks (sticks), long periodic bonds hidden, side-on start for tubes. | `viewer.py` |

### Problems found and fixed while building it
- Periodic KD-tree queries failed on coordinates exactly at the box edge
  (clamped after wrapping).
- Rupture ethers collapsed to 62° C-O-C angles under the bond-only spring
  model (angle springs with per-triple weights added).
- Requested 10 % vacancies gave 18.6 % on small sheets because edge cleaning
  after the last hole was not counted. Holes are now cleaned one at a time,
  the cleaned atoms count toward the target, and an overshooting last hole is
  undone.
- Dense oxidation stalled at 53-57 % O (sublattice bias, hydroxyl-first
  ordering and diol densification added; H...O down to 1.6 Å allowed).
- MWCNT/onion O % came out low because only the outer wall was counted in the
  target (target now scaled to all framework carbons).
- Edge H could end up 0.64 Å from a group after post-oxidation cleanup (edge
  termination made iterative with clash-driven removal).
- Epoxide opening placed 0-3 of 8 ligands (1-4 pairs across the new bond were
  treated as clashes; now half radius, plus 0.40 Å link puckering).
- PDB REMARK lines exceeded 80 columns (split).

### Verification
- `tests/test_carbon.py` (about 30 tests) and carbon GUI tests; the full
  suite passes (see README for the count).
- 24 randomized stress builds: no over- or under-valent atoms; bonded
  distances 0.97-1.56 Å; closest non-bonded contacts H...O 1.60, H...H 1.65,
  C...H 2.04 Å.
- Fresh 2:1 sweep (10 × 10 nm): O % exact from 10 to 68 %; ratio exactly 2:1
  up to 40 %, shifting toward hydroxyls above (reported).
- 20 × 20 nm flake at 68 % O: 34,924 atoms in about 20 s.
- Rendered viewer pages of nine builds (GO flakes, aged GO with holes,
  pristine defects, CNT-aryl, C60-amine, onion, GO amides, epoxide opening,
  oxidized MWCNT) were inspected visually.

### Not verified
- tleap with real GO parameters and charges (the MOL2 types and residues
  follow the AMBER GO dataset's conventions, but no frcmod is bundled).
- antechamber agreement for every library ligand's provisional GAFF2 types.

## Buddy Builder round (ease of use, hackathon, about 3 hours)

Scope came from the "xtalligand GUI: Ease-of-Use Improvements" doc and the
answers recorded in it; multi-day items (offline Materials Project set, mixed
ligand shells, full drag-to-move editing) were left out.

### Changed
| Area | What | Files |
|---|---|---|
| Brand | Name **Buddy Builder**, round logo in the header, pastel **Buddy** theme (default) beside Light / Dark / Auto, in the GUI and the viewer. | `gui.py`, `theme.py`, `viewer.py`, `assets/buddy_logo.png` |
| Layout | Basic / Advanced switch (seeds, tolerances, lattice constants, vacuum, box edges, carbon details are Advanced); icon toggle buttons for source, shape, channel type, box shape, coverage mode, arrangement; capped widths and wrapping rows (no clipped inputs); wrapping tab bar; status says what is loaded and what to do next; Simulation-ready / charges / explicit-H controls removed from the layout; MOL2 residue scheme shown only for carbon. | `gui.py`, `gui_carbon.py`, `theme.py` |
| Ligands | Library filtered by material (covalent carbon chemistry only on carbon, metal ligands only on crystals); switching material swaps an unsuitable ligand for 4-carboxyphenyl (carbon) or 2-mercaptoethanol (crystal) and says so. | `ligand_library.py`, `gui.py` |
| Placement | **Fibonacci** order (golden spiral on spheres, golden-angle helix on cylinders and pores, golden-ratio lattice on flat faces, farthest-point elsewhere; nearest free site within half a spacing, then gap fill). Default for new GUI builds; stored configurations keep their own order. | `fibonacci.py`, `placement.py` |
| Boxes | Auto / Rectangular / Cube. Slabs and channels: orthogonal (or square) periodic supercell from integer combinations of the surface vectors, edge atoms trimmed so it repeats exactly; atom density unchanged. Cube makes finite boxes cubic and slab boxes at least as tall as the side. Optional minimum box edges. | `shapes.py`, `cell.py`, `pipeline.py` |
| Channels | Crystal **slit channels** (1-20 stacked gaps between slab walls: multi-layer channels) and **cylindrical pores** through a membrane slab; region "Channel walls" (inner). Carbon: region "inner" between graphene layers and on the innermost nanotube wall. | `shapes.py`, `sites.py`, `carbon.py` |
| Viewer | Click menu: inspect (label), measure distance / angle (nearest image), remove ligand, add ligand on a surface atom; Top / Side / Reset views; Clear labels; PNG. Channels start side-on. Clicks reach Python through a hidden text widget in the same page. | `viewer.py`, `gui.py` |
| Edits | Remove / add rebuild the assembly from the placement list (same bonding, typing and validation), keep the camera, re-export, Undo; every edit is listed in `provenance.manual_edits` and under MANUAL EDITS in README.txt. | `editing.py`, `reporting.py`, `topology.py` |
| Console | Python vocabulary (`material`, `carbon`, `shape`, `box`, `ligand`, `coverage`, `export`, `build`, `select`, `measure`, `info`, `remove_ligand`, `add_ligand`, `undo`, `summary`, `examples`, `theme`, `mode`, `config`, `help`); input parsed with `ast`, literals only; Generate clicks are echoed as commands. | `console.py`, `gui.py` |
| Examples | K (Au(111) slit channels x2, ligands on the walls) and L (4 nm Au sphere, 60 dodecanethiols, Fibonacci, cubic box). | `examples.py`, `examples/*.json` |

### Verification
- Full suite: 232 passed, 1 skipped (new: `tests/test_buddy_features.py`,
  GUI tests for defaults, Basic/Advanced, ligand filtering, channels, boxes,
  viewer edits and Undo, console).
- Orthogonal slab cells keep the native areal atom density for Au(111),
  Au(211), Au(100), Mg(0001), Mg(101), NaCl(111), Si(111).
- Nearest-neighbour spread of ligand sites (coefficient of variation),
  Fibonacci vs random: Au sphere 0.13 vs 0.31, slab 0.09 vs 0.26, CNT 0.06 vs 0.39.
- End-to-end matrix with all five exporters: spheres, cube box, rod, rect
  Au(211) and square Au(111) slabs, slit x2 / x3, pore x2, GO, CNT inner,
  C60, Fe sphere: all ok. Graphene stack at 9 Å spacing with carboxyphenyl
  "inner" places 0 of 8 (the ligand does not fit; reported as partial).
- Headless Notebook 7: screenshots of every theme and the main tabs; a
  scripted click in the viewer removed a ligand through the bridge (12 to 11,
  downloads refreshed) and a two-atom measurement read 2.885 Å (Au-Au).

### Not verified
- The viewer-to-Python bridge in nanoHUB AppMode itself (works in Notebook 7;
  AppMode uses the same page, and the console offers the same edits).

## Buddy Builder round 2 (feedback on the first Buddy version)

### Changed
| Request | What was done | Files |
|---|---|---|
| Remove Auto theme; green Buddy theme | Themes are Buddy (pastel yellow), Buddy 2 (mint green), Light, Dark. | `theme.py`, `gui.py`, `viewer.py` |
| Clearer tabs | Tabs are bold buttons in a grid; the current one is filled. | `theme.py` |
| Browse above Find; hover groups; no group names in the list | "Browse materials" hover menu (group, then material; clicks reach Python through a hidden text widget); the search list shows bare names; typing a word of a name ("anatase") works. | `gui.py`, `material_catalog.py` |
| No preset sizes in the material list | Labels name the material only; sizes come from the Shape tab and the box. | `material_catalog.py` |
| Box lengths instead of vacuum / max atoms | Shape tab starts with the simulation box: Rectangular (x, y, z), Cube (edge, 8 nm default for new builds) or Auto. Slabs, channels and periodic carbon fill x and y (whole repeats); z is exact; finite models are centred and must fit (2 A clearance, error otherwise). | `cell.py`, `shapes.py`, `pipeline.py`, `gui.py`, `gui_carbon.py` |
| Unit cell info before a box exists | Picking a material fills a unit-cell table (space group, a b c, angles, atoms per cell, volume, density, source) and previews the cell (angled, orthographic). | `gui.py`, `viewer.py` |
| Viewer | Top / Side / Angled / Reset; Ortho / Persp (ortho default for periodic models); bolder box; Slice slider (+ thin slab); hi-res PNG (3x off-screen render); background white / black / theme / transparent; substrate-ligand bonds drawn. | `viewer.py` |
| Pores across the surface | Pores go through the membrane at golden-ratio positions relaxed apart, >= 3 atomic spacings from each other and the box edges; errors for pores as large as the box or too many; Fibonacci targets per pore. | `shapes.py`, `fibonacci.py` |
| Validation report -> Save session | `session.py` in every export and as a download: comments with the short commands, `settings({...})` with every non-default setting (uploaded CIF text embedded), `build()`, hand edits. Console `settings()` / `session()`; Load session upload (also config.json); the console input is a text area (Enter runs, pastes keep lines). | `session.py`, `console.py`, `gui.py`, `reporting.py` |
| Certificate of Build | `certificate_of_build.pdf` in every export (matplotlib): certificate page with logo, name, key numbers, build ID and a picture; detail pages; three views. | `certificate.py`, `reporting.py` |
| Materials Project loading problems | Live search removed from the GUI; the 30 preloaded CIFs in `data/materials_project` are used offline (all 30 load and build). | `mp_library.py`, `material_catalog.py`, `crystal.py`, `schema.py` |
| Ligand-surface bonds | "Geometric adsorption" renamed Physisorption (no bond); new Chemisorption template writes the anchor-surface bond (covalent-radii distance, S-H/O-H lose H, lone-pair orientation); thiols on non-Au metals use it; the viewer draws these bonds (Au-S too). | `chemistry.py`, `ligands.py`, `topology.py`, `placement.py`, `validation.py`, `viewer.py` |
| Ragged periodic edges | Periodic boundaries are moved into the widest empty band between atoms (pure translation) and periodic models open in orthographic view, so cuts look straight. | `cell.py`, `shapes.py`, `viewer.py` |
| Sad logo on failure | A round sad-Buddy image pops up for about a second (CSS fade) when a build fails or fails validation. | `gui.py`, `theme.py`, `assets/buddy_sad.png` |

### Verification
- Full suite after the round: see README test count; new tests cover the
  library, chemisorption bonds, physisorption without bonds, user box size
  and fit errors, clean periodic cuts, pore spreading and errors, session
  round trip, certificate PDF, and the GUI (menu, library, box controls,
  session reload, sad popup).
- Headless Notebook 7: hover menu click loaded anatase; Enter in the console
  ran commands; a failing build showed the sad popup and it was hidden 2 s
  later; Buddy 2 and Dark screenshots checked; hi-res PNG export with a
  transparent background gave a 5400 x 3678 px RGBA image.

## Buddy Builder round 3 (polish)

- Tabs: icon and close boxes hidden and labels centred in each button.
- Failed builds: the popup keeps "Oh no! Buddy could not build that." and adds
  why and how to fix (pattern table in `gui.failure_advice`); it now stays about
  4 s so the text can be read; Status shows the full message plus the same
  advice; the console `build()` prints them too.
- Certificate: dark purple outer border, the new Buddy picture, the structure
  drawn without a backing circle (transparent), the Validation box removed, and
  the save-as name under the signature instead of "chief builder".
- Export tab: Save as (ExportConfig.save_as) names the ZIP, certificate and
  session downloads; empty gives a Buddy name from Southern terms of endearment
  (`reporting.buddy_name`, the same for the same build).
- Fonts: Fredoka for titles, section headings, tabs, buttons and download
  buttons; Arial for labels, captions and inputs (GUI and viewer).
- Certificate colours and fonts: purple lettering (title, structure name,
  signature, headings); detail and view pages get the poster banner (peach and
  cream awning stripes over a light yellow band). Fredoka is bundled
  (`assets/fonts`, SIL OFL: static Medium and SemiBold for matplotlib, a Latin
  woff2 subset embedded in the GUI and viewer CSS), so titles look the same
  online, offline and in the PDF; body text uses Arial (Liberation Sans or
  DejaVu Sans where Arial is not installed).
