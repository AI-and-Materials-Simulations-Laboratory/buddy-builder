# User guide

This is a reference for every option in Buddy Builder. For a first run, see [GETTING_STARTED.md](../GETTING_STARTED.md). For the AI assistant, see [ai_assistant.md](ai_assistant.md).

## Workflow

The left column holds **Ask Buddy** and **tabs 1–6**, in order: Material, Shape, Defects & GO, Ligand, Coverage and Export. The right column holds the 3D viewer, downloads, results and the Python console.

In the header:

- **Settings:** *Basic* shows the choices most builds need. *Advanced* adds seeds, tolerances, lattice constants and carbon details.
- **Theme:** Buddy (pastel yellow, the default), Buddy 2 (mint green), Light or Dark.

## Options by step

### 1. Material

- **Browse materials:** a hover menu, first by group, then by material.
- **Find material:** type a name, a word of a name ("anatase"), a formula ("TiO2") or a carbon word ("GO", "CNT", "C60").
- **Built-in presets:** Au/FCC (the default), Ag, Cu, Fe, Si, Mg, Ge and NaCl.
- **Preloaded Materials Project library:** 30 structures, read offline with no key or network.
  - Metals: Pt, Pd, Ni, Al, Co, Ti.
  - Intermetallics: AuCu3, AuCu, FePt, CoPt, Ni3Al, CuZn, NiTi, FeAl.
  - Oxides: TiO2 (anatase and rutile), ZnO, Fe3O4, CeO2, MgO, Cu2O, SiO2, Al2O3.
  - Semiconductors: CdSe, CdS, ZnS, GaAs, InP, PbS.

  The structures are in `src/nanobuilder/data/materials_project`. To add one, drop its CIF there and add an entry to `catalog.json`.
- **Upload:** a fully ordered CIF or POSCAR file.
- **Carbon nanomaterials:** graphene and graphene-oxide sheets, carbon nanotubes, fullerenes and carbon onions.

Picking a material shows its unit cell (lattice, angles, atoms per cell, volume and density) and previews it in the viewer.

### 2. Shape and simulation box

- **Sphere:** set its diameter.
- **Nanorod:** set the diameter, the tip-to-tip length, the lattice axis [u v w], and hemispherical or flat caps. Ligands go on the sides, the tips or both.
- **Slab:** set the Miller indices (hkl), the size, the number of layers or the thickness, and the termination.
- **Channel:**
  - *slit channels*: 1–20 channels stacked between slab walls;
  - *cylindrical pores*: 1–20 pores straight through a membrane. They are spread on a golden-ratio pattern, at least 3 atomic spacings apart and from the box edges.

  Ligands go on the channel walls.
- **Simulation box** (set it first):
  - *Rectangular*: x, y and z in nm;
  - *Cube*: one edge length (8 nm by default for new builds);
  - *Auto*: fits the model plus 3 nm of vacuum.

  Slabs, channels and periodic carbon fill x and y with a whole number of crystal repeats. Buddy moves the periodic cut into the widest gap between atoms, so the cut is clean. Finite models are centred and must fit inside the box.
- **Carbon sizes:**
  - sheets: size and layers (AB or AA stacking);
  - nanotubes: (n, m), length and number of walls;
  - fullerenes: C60·k²; onions: number of shells.

  Each is either finite (H-terminated, AMBER-ready) or periodic.

### 3. Defects and graphene oxide (carbon)

This is a re-implementation of HierGO, with improvements.

- **Defects:**
  - vacancies as holes (size classes, random growth, or HierGO's hole library);
  - circular pores;
  - Stone-Wales defects;
  - mini grain boundaries (24-atom rotations);
  - ruptures with ether oxygen.
- **Oxidation**, from 0 to about 68 % O:
  - fresh hydroxyl:epoxy (2:1 by default, as in the AMBER GO dataset);
  - HierGO's aged 600–2500 K profiles, with edge COOH, carbonyl, phenol and ether groups;
  - or a custom mix.
- **Patterns:** uniform, edges-first, clustered domains or tiles, with a choice of faces.

Buddy reports C:O, O %, the percentage of modified carbons, and every fallback it used. See [carbon_materials.md](carbon_materials.md).

### 4. Ligand

- **Library:**
  - alkanethiols;
  - MUA, MPA and MUD;
  - cysteamine and PEG-thiols;
  - amines, citrate and oleate;
  - TPP and TOPO;
  - covalent carbon chemistry.

  Only ligands that suit the material are listed.
- **Any SMILES:** pick the anchor atom from a numbered 2D drawing. Chains can be built all-trans (extended).
- **Attachment:**
  - **Au-thiolate:** on Au, with an Au–S bond.
  - **Chemisorption:** on any crystal, with an anchor–surface bond at the sum of the covalent radii. S–H and O–H anchors lose their H; N and P keep theirs.
  - **Physisorption:** no bond.
  - **Covalent bonds to carbon:**
    - basal sp2→sp3 addition;
    - edge C–H substitution;
    - amide or ester to GO carboxylic acids;
    - epoxide ring-opening.

### 5. Coverage

- Set an exact ligand count, or ligands per nm².
- **Arrangement:**
  - *Fibonacci* (ordered; the default for new builds);
  - *evenly spaced*;
  - *random*.

  Fibonacci puts targets on a golden spiral (spheres), a golden-angle helix (rods, tubes, pores) or a golden-ratio lattice (flat faces), then places each ligand on the closest free site. The result is the same on every run.
- Ligands that clash, including with periodic images, are rejected.
- Carbon stacks and nanotubes also offer the **inner** region.

Both the requested and the achieved coverage are reported.

### 6. Export

- **Save as** names the downloads. Left empty, Buddy makes up a name such as `Sweet_Pea_Buddy_Au_sphere`.
- Choose the formats, and the residue scheme for MOL2 on carbon.

## Viewer

**Viewing:**
- Rotate, zoom and pan.
- Show or hide ligands, and highlight sites.
- Top, Side, Angled and Reset views; orthographic or perspective.
- **Slice** to look inside the model.
- Background: white, black, theme or transparent.
- High-resolution PNG export (3x).

**Click menu:**
- **Inspect atom.**
- **Measure:** click two atoms for a distance, three for an angle.
- **Remove ligand.**
- **Add ligand** on a chosen surface atom.

Edits rebuild and re-validate the model and refresh the downloads. They are listed under **MANUAL EDITS** in the export's `README.txt`. **Undo edit** steps back.

## Python console

The console takes the build as short Python calls:

```python
material("Au"); box("cube", edge=10); shape("sphere", diameter=4)
ligand("dodecanethiol"); coverage(count=60, order="fibonacci")
build()
select(element="S"); measure(10, 25); remove_ligand(3); undo()
material("anatase"); box("rectangular", size=(6, 6, 8)); shape("channel", type="pore", width=1.2, count=3)
session()     # print the commands that rebuild the current model
ask("a 3 nm silver sphere with 15 MUA")     # plan with Ask Buddy
reset()       # back to the launch settings
```

- Enter runs the input; Shift+Enter adds a line.
- `help()` lists every command, and `help("select")` explains one.
- Commands drive the same widgets as the tabs. Each **Generate model** click is echoed as the equivalent commands.
- Only the listed commands with plain values can run. Input is parsed, never `eval`-ed.
- **AI** turns plain English in the box into tested commands, and **Test** dry-runs the box on a hidden copy.

## Sessions, certificates and examples

- **Save session (.py):** the build as console commands, plus `settings({...})` with every non-default choice, `build()` and any hand edits. Paste it into the console, or use **Load session**, to restore it exactly. **Load session** also accepts a `config.json`.
- **Certificate of Build (PDF):** a certificate page, then pages with the structure, crystal, atoms, ligands, coverage, topology and validation details, plus three views.
- **Examples:**
  - **A–E**: the five required demo cases;
  - **F–J**: carbon (GO flake, aged GO with holes and amides, an aryl-grafted nanotube, a C60 amine adduct, epoxide opening on periodic GO);
  - **K**: Au(111) slit channels;
  - **L**: a 4 nm Au sphere with 60 dodecanethiols in Fibonacci order.
- **When a build fails,** sad Buddy appears with **why** it failed and **how to fix** it; the same text stays under Status.

## Look and feel

Titles, headings, tabs and buttons use Fredoka (bundled, SIL Open Font License, so it works offline). Other text uses Arial. The Certificate of Build uses the same pair.
