# Carbon nanomaterials: graphene, graphene oxide, nanotubes, carbon spheres

This page describes the carbon builder added to xtalligand: what it builds, how
each option maps to the HierGO scripts and to the AMBER graphene-oxide dataset,
how ligands are bonded covalently, and what the exports contain.

Sources the implementation follows:

- **HierGO**: Garcia et al., *2D Materials* 10, 025007 (2023), and the scripts in
  `HierGO-master/Scripts` (create-tile, add-sw-defects, decorate-tile,
  create-rupture, trim-hanging-Cs, atom-typing).
- **AMBER GO dataset**: Fedai et al., *Data* 11, 18 (2026), which adapts the
  HierGO workflow for AMBER. It uses a 2:1 hydroxyl:epoxy ratio for fresh GO,
  non-periodic flakes, minimal GAFF types (`ca`, `c3`, `oh`, `os`, `ho`),
  residues `GRA`/`OH`/`EPX` (or one `GRO` residue), and validates builds by C:O
  ratio, overall modification `100·N_O/N_C` and modified carbon
  `100·N_c3/(N_ca+N_c3)`.

> **License note.** HierGO is GPL-3.0. Its algorithms were re-implemented here
> from the paper and the documented script behavior; no HierGO code is copied
> into `src/nanobuilder`. The optional hole-shape library
> (`HierGO-master/HoleLibrary/defect-*.xyz`) is **read at run time** from the
> HierGO folder next to the app, or from the folder given by
> `NANOBUILDER_HOLE_LIBRARY`. It is never bundled into the package or the
> export ZIPs. Without it, the "HierGO hole library shapes" option reports an
> error and random growth still works.

## 1. Using it in the GUI

Choose **Carbon nanomaterial** as the material source on tab 1, or type
`graphene`, `GO`, `CNT`, `C60`, `fullerene` or `onion` into **Find material**
and pick an entry from the **Carbon nanomaterials** group (14 ready-made
starting points). The tabs then become:

| Tab | Carbon controls |
|---|---|
| 1 Material | Carbon material: sheet, nanotube or carbon sphere |
| 2 Geometry | Sheet: length x (zigzag) and width y (armchair) in nm, layers, AB/AA stacking, layer spacing. Nanotube: chirality (n, m), length, walls, wall spacing. Sphere: single fullerene C60·k² (k = 1-10), or an onion of nested shells starting from C60. All: boundary (finite or periodic), edge termination (H or bare). |
| 3 Defects & GO | Vacancies (%) with hole size class and hole shapes, circular pores, Stone-Wales (%), mini grain boundaries / 24-atom rotations (%), ruptures, hanging-carbon trimming; then oxidation: oxygen %, profile, hydroxyl:epoxy ratio, paired-diol share, custom group mix, spatial pattern (uniform, edges first, domains, tiles), oxidized faces, every-layer oxidation, sp3 puckering, seed. |
| 4 Ligand | Library group **Covalent on carbon**, or any SMILES. The attachment list shows only the templates valid for a carbon substrate. |
| 5 Coverage | Count or density; regions: sheet Top / Bottom / Both faces / Edges / Entire surface; nanotube Side / Tips / Entire surface; sphere Entire surface. |
| 6 Export | Adds **MOL2** (AMBER-oriented, suggested GAFF2 types) and the residue scheme (one residue per group, or one `GRO` residue for the whole substrate). |

The info lines under each section show the atom estimate, nominal area and how
the request will be interpreted before you click **Generate model**. After a
build, the results table lists the achieved composition (C, O, H), oxygen
content and C:O, modified sp3 carbon %, basal OH:epoxide ratio, vacancy %, and
the count of every oxygen group. Every requested-vs-achieved difference also
appears as a warning.

Example configurations F-J (in **Load an example...**, and as JSON in
`examples/`):

| Example | What it shows |
|---|---|
| F `graphene_oxide_flake.json` | 5×5 nm finite GO flake, 30 % O, fresh 2:1 (the AMBER GO recipe), no ligands |
| G `aged_go_holes_amide.json` | 6×6 nm aged (1200 K) GO with 8 % vacancies and Stone-Wales defects, oxidized edges first; 6 ethylenediamine amides on edge COOH |
| H `cnt_aryl.json` | Periodic (10,10) SWCNT with 10 evenly spaced 4-carboxyphenyl groups (diazonium-type grafting) |
| I `c60_amine.json` | C60 with 6 butylamine adducts |
| J `periodic_go_epoxide_opening.json` | Periodic GO sheet whose epoxides are opened by ethanolamine |

Command line: `python -m nanobuilder.cli examples/graphene_oxide_flake.json --out build`.

## 2. Frameworks

All carbon frameworks are built with an **explicit bond graph**: every C-C, C-O
and O-H bond is created by the builder, not guessed from distances. The graph is
what the exporters write, so heavily oxidized or defective regions never lose
or gain bonds.

**Graphene.** Orthogonal cell with a = √3·1.42 Å along x (zigzag) and
b = 4.26 Å along y (armchair); 2.619 Å² per carbon. The requested size is
rounded to whole cells. Several layers use AB (Bernal) or AA stacking at the
chosen spacing (3.35 Å default). A **finite** sheet is a flake: atoms with fewer
than two carbon neighbours are trimmed and every edge carbon gets an H at
1.09 Å (or is left bare). A **periodic** sheet has bonds across the cell
boundary and no edges. Finite flakes are the AMBER route; periodic sheets are
for LAMMPS, GROMACS and HOOMD.

**Nanotubes.** ASE's `nanotube(n, m)` is used for each wall. Multi-wall tubes
pick, for every outer wall, the chirality of the same type whose diameter is
closest to the previous one plus 2 × wall spacing. Finite tubes have
H-terminated rims ("tips"); periodic tubes are infinite along z
(`pbc = (F, F, T)`).

**Carbon spheres.** Icosahedral Goldberg fullerenes C60·k² (C60, C240, C540,
...). They are made from a geodesic point set, its convex hull and the dual
polyhedron (12 pentagons, the rest hexagons), then relaxed to 1.42 Å bonds with
a bond and angle spring model. Carbon onions nest C60, C240, C540, ... which
are naturally about 3.5 Å apart.

## 3. Defects (HierGO create-tile, add-sw-defects, create-rupture)

| Control | Algorithm | HierGO equivalent |
|---|---|---|
| Vacancies (% of C) | Holes are grown one at a time by random growth from a seed atom, until the removed-carbon count reaches the target. Size classes: small 1-5, medium 5-10, large 10-30, mixed 1-30 atoms, or a custom range. After each hole the edge is cleaned (see below) and the cleaned atoms count toward the target. The last hole is undone if it would overshoot by more than it would undershoot. | `create-tile.py --percentvac --distribution` |
| Hole shapes: HierGO hole library | Stamps shapes from `HoleLibrary/defect-<n>.xyz` (rescaled from a = 2.471 Å to this lattice) instead of growing them. | the same library files |
| Circular pores | Removes all carbons within d/2 of randomly placed, non-overlapping centres. | (not in HierGO) |
| Stone-Wales (%) | 90° rotation of one C-C bond: two pentagons and two heptagons. Count = % × N_C / 16 (HierGO's convention). | `add-sw-defects.py --sw` |
| Mini-GB, 24-atom rotations (%) | A coronene-like C24 unit rotated by 30° and reconnected, by optimal assignment, to its 12 outer neighbours: six pentagons and six heptagons. Count = % × N_C / 41. | `add-sw-defects.py --gb` |
| Ruptures | A line of broken C-C bonds bridged by alternating ether O atoms (C-O-C at about 1.37 Å, above and below the plane). | `create-rupture.py --atoms` |
| Keep hanging carbons | Off by default. Otherwise carbons left with one neighbour, and rings of 8 or fewer atoms hanging off a hole, are removed. | `trim-hanging-Cs.py`, `--keep_hanging_carbon` |

After topological defects, a **local spring relaxation** (bonds 1.42 Å, angles,
and ether geometry) removes the stretched bonds that HierGO leaves for MD to
fix. Edge carbons whose hydrogens could not fit (too close to another H or
group) are removed so that every edge carbon keeps a clean, non-overlapping H.

The **achieved vacancy %** is relative to the pristine lattice of the same size
(after the initial edge trim of a finite flake). It is reported and checked
against the request. Small sheets with large holes can deviate by about 1-2 %,
because a hole is removed whole.

## 4. Oxidation (HierGO decorate-tile, AMBER GO dataset)

**Oxygen content** is `100 · N_O / N_C` over the whole substrate, the quantity
HierGO calls "percent oxygen" and the AMBER GO dataset calls overall
modification. C:O is `N_C / N_O`. **Modified carbon %** is
`100 · N_sp3 / (N_sp2 + N_sp3)` over framework carbons. All three are reported
for every build.

**Profiles**

| Profile | Groups (share of added O) |
|---|---|
| Fresh | Basal hydroxyl and epoxide in the hydroxyl:epoxy ratio given (default **2:1**, as in the AMBER GO dataset; HierGO's default is 1:1). The *paired share* (default 0.7) is the fraction of hydroxyl O placed as anti 1,2-diols on neighbouring carbons. |
| Aged 600 / 900 / 1200 / 2500 K | HierGO's annealing compositions, as below. |
| Custom | Any relative mix of basal hydroxyl, epoxide, edge hydroxyl (phenol), edge carbonyl, carboxylic acid (counts as 2 O) and in-plane ether. |

| Aged profile | Paired OH | Single OH | Edge OH | Epoxide | Ether | Carbonyl | COOH |
|---|---|---|---|---|---|---|---|
| 600 K | 0.19 | 0.16 | 0.04 | 0.43 | 0.05 | 0.07 | 0.12 |
| 900 K | 0.15 | 0.16 | 0.05 | 0.34 | 0.08 | 0.11 | 0.22 |
| 1200 K | 0.12 | 0.14 | 0.06 | 0.22 | 0.18 | 0.14 | 0.28 |
| 2500 K | 0.02 | 0.02 | 0.04 | 0.15 | 0.47 | 0.15 | 0.30 |

Edge groups (ether, COOH, carbonyl, phenol) need edge carbons: holes, ruptures,
flake edges or tube rims. If there are too few, their oxygen is added as basal
hydroxyls. The same happens to epoxides that find no free C-C pair at high
coverage, which is the fallback used in the AMBER GO workflow. Each fallback is
reported.

**Spatial patterns**

| Pattern | Meaning | HierGO |
|---|---|---|
| Uniform random | Seeded random order over the oxidizable carbons. | default |
| Edges and holes first | Carbons nearest an edge are decorated first. | `--decorate_edges_first` |
| Clustered domains | N circular domains of diameter d are filled first. | `--bunching` |
| Tiles | The sheet is split into nx × ny tiles. Each gets its own O %, drawn at random between min and max. Per-tile and overall O % are reported. | separate tiles + `stitch.py` |

**Faces.** Sheets are oxidized on both faces with balanced counts by default,
or on the top or bottom face only. Multi-layer sheets oxidize only their two
outer faces, unless "Oxidize every layer" is ticked (it needs ≥ 6 Å spacing).
Nanotubes and onions are oxidized on the outer surface of the outermost wall or
shell only. Their O % still counts every carbon.

**Geometry rules** (improvements over HierGO's fixed offsets):

- Bond lengths: C-OH 1.43 Å, O-H 0.97 Å, epoxide O 1.276 Å above the C-C
  midpoint, phenol C-O 1.36 Å, C=O 1.23 Å, COOH C-C 1.50 / C=O 1.21 / C-O
  1.35 Å, edge C-H 1.09 Å.
- A carbon bearing a group is puckered 0.25 Å toward it (sp3 pyramidalization).
  This can be switched off.
- No two neighbouring carbons carry a group on the same face; paired diols are
  anti (one up, one down).
- New atoms keep minimum distances from every existing atom: O-O 2.40, H-O 1.60
  (hydrogen bonds allowed), H-H 1.50, C-O 2.20, C-H 1.80 Å. Hydroxyl H
  directions are chosen from 12 azimuths.
- Above about 25 % O, groups are biased toward alternating sublattices (the
  chair pattern of fully oxidized graphene). This is what allows dense targets.
- If the sheet jams before the target, epoxides are opened into anti diols
  (one extra O on the same two carbons).

Measured on a 10 × 10 nm fresh 2:1 flake (seed 7):

| Requested O % | Achieved O % | Basal OH : epoxide O | Modified C % |
|---|---|---|---|
| 10 | 10.0 | 2.0 | 13.3 |
| 30 | 30.0 | 2.0 | 40.0 |
| 40 | 40.0 | 2.0 | 53.3 |
| 50 | 50.0 | 2.4 | 64.8 |
| 60 | 60.0 | 4.5 | 70.8 |
| 68 | 68.0 | 8.4 | 75.3 |

The ratio holds exactly up to about 40 % O. Above that, there is not enough room
for that many epoxides, so the builder shifts toward hydroxyls to reach the
oxygen target, and the achieved ratio is reported. A 20 × 20 nm flake at 68 % O
(35k atoms) builds in about 20 s. Small flakes (≤ 6 nm) can stop slightly below
very high targets (about 66.5 % for 68 %), with a warning.

## 5. Covalent ligands

Ligands bond to the carbon framework through four attachment templates. Each
removes one H from the ligand's anchor atom. The new bond points along the
anchor's former X-H direction, and the substrate is edited at that site:

| Template | Site | Substrate edit | Bond length (C-C / N / O / S, Å) | Typical chemistry |
|---|---|---|---|---|
| `C_basal` | sp2 basal carbon (sheet faces, tube wall, fullerene) | C becomes sp3 and is puckered 0.40 Å outward | 1.53 / 1.47 / 1.43 / 1.82 | aryl diazonium grafting, radical or amine addition, C60 adducts |
| `C_edge` | H-terminated edge or rim carbon | the edge H is replaced | 1.50 / 1.40 / 1.36 / 1.77 | edge functionalization |
| `C_carboxyl` | COOH on edges and holes (oxidized carbon) | the acid OH is removed; acyl C bonds the anchor | N or O 1.34, S 1.77 | EDC/NHS amide or ester coupling |
| `C_epoxide` | epoxide on GO | the ring opens; the O stays on one carbon as OH, the ligand bonds the other | N 1.47 / O 1.43 / S 1.82 | amine ring-opening |

On flat sheets the epoxide opens **anti**: the ligand bonds from the face
opposite the epoxide O, so a "Top" ligand uses a bottom-face epoxide. On curved
walls only the outside is reachable, so the ligand bonds on the same face,
tilted away from the new OH. `C_basal` and `C_edge` accept C, N, O, S, P, Si or B
anchors; `C_carboxyl` and `C_epoxide` accept N, O or S. The anchor must carry
an H or be a radical. The
**Covalent on carbon** library has phenyl, 4-nitrophenyl, 4-carboxyphenyl,
4-aminophenyl, methyl, ethylenediamine, octadecylamine, mPEG4-amine,
ethanolamine and n-butylamine. On carbon, dodecylamine, oleylamine and
cysteamine from the other library groups switch to `C_carboxyl` (amide), and
thiols are placed by geometric adsorption. `Au_thiolate` is refused on carbon, and the
carbon templates are refused on metals and crystals. Geometric (non-bonded)
adsorption remains available.

The clash check treats 1-4 neighbours across the new bond with half the usual
radius, and two ligands never claim the same carbon or group. After assembly,
validation checks every atom's valence, that every framework carbon keeps at
least two bonds, and that no non-bonded pair is closer than 0.9 Å. The composition is reported again after attachment, because
ligands change the sp3 count and remove OH from COOH.

## 6. Typing and exports

Each atom's structural role (for example `C_sp2`, `C_epoxy`, `O_hydroxyl`,
`H_edge`, `C_linked`) comes from how it was built, so substrate types are exact
rather than perceived:

| Role | GAFF2 | OPLS-AA (HierGO typing) |
|---|---|---|
| sp2 basal / edge C | ca | opls_147 / opls_145 |
| edge H | ha | opls_146 |
| C-OH carbon, O, H | c3, oh, ho | opls_159, opls_154, opls_155 |
| epoxide C, O | c3, os | opls_184, opls_180 |
| ether C, O | ca, os | opls_518, opls_180 |
| carbonyl C, O | c, o | opls_280, opls_281 |
| COOH C, =O, -O(H), H | c, o, oh, ho | opls_267, opls_269, opls_268, opls_270 |
| phenol C, O, H | ca, oh, ho | opls_166, opls_167, opls_168 |
| ligand-linked framework C | c3 | - |

Ligand atoms get **provisional** rule-based GAFF2 names (aromatic `ca`/`ha`,
`c3`/`hc`/`h1`, amine `n3`/`hn`, amide `n`, carbonyl `c`/`o`, etc.). Check them
with antechamber. All types are in `topology.json`
(`suggested_atom_types.gaff2` and `.opls_hiergo`).

| File | Carbon-specific content |
|---|---|
| `mol2/structure.mol2` | GAFF2 types in the type column, bond orders `1`/`ar`/`2`/`am`, residues `GRA` (each framework carbon with its edge H), `OH`, `EPX`, `ETH`, `CBO`, `CBX`, `LIG`; or `GRO` for the whole substrate (the AMBER GO dataset's choice when GO is combined with biomolecules). `CRYSIN` for periodic cells. Charges are 0.0000 (`NO_CHARGES`). With `README_mol2.txt` and a `tleap_template.in` sketch. |
| `pdb/structure.pdb` | the same residue names, unique atom names, every bond in CONECT |
| `lammps/structure.data` | bonds, angles and dihedrals from the full bond graph; boundary note `p p f` for periodic sheets and `f f p` for infinite tubes |
| `hoomd/structure.gsd`, `structure.extxyz`, `topology.json` | the same graph; extXYZ has a `residue` column; `topology.json` has bond orders, residues and suggested types |
| `README.txt`, `manifest.json` | the composition, oxygen groups, defect counts, achieved vs requested values, and the full carbon configuration |

AMBER cannot represent bonds that cross a periodic boundary, so for AMBER build
finite flakes or tubes. The MOL2 export warns when the model is periodic.
Impropers are not enumerated, and no charges or force-field parameters are
generated. The AMBER GO dataset's GO charges and frcmod, plus ligand charges
(for example AM1-BCC), must be supplied by the user.

## 7. Differences from HierGO

- One build makes the whole model: sheet, tube or sphere, then defects,
  oxidation, ligands and export. There is no manual stitch, remove-atoms or
  re-typing loop, and the bond graph is never re-perceived.
- Nanotubes, fullerenes and onions are new. HierGO builds planar tiles only.
- Chemically sized groups, sp3 puckering and local relaxation around
  topological defects. HierGO uses fixed offsets.
- Dense oxidation (up to about 68 % O) is reached by the sublattice bias and
  diol densification. HierGO notes that it "may fail" above about 40 %.
- Achieved composition, C:O, modified carbon % and every fallback are reported.
  Nothing is silently dropped.
- Tiles are oxidized inside one sheet, so no stitching overlaps need trimming.
- HierGO's 1:1 default fresh ratio becomes 2:1, following the AMBER GO dataset.
  Type `1:1` to reproduce HierGO.

## 8. Limitations

- Starting structures only: oxidized sheets stay flat apart from the puckering,
  and nothing is energy-minimized. Minimize before production MD.
- Ligand GAFF2 types are provisional; partial charges are not assigned.
- Vacancy % can overshoot by about 1-2 % on small sheets with large holes or
  library shapes, and with ether-rich aged profiles. The achieved value is
  always reported.
- Very high oxygen targets on small flakes can fall slightly short (about 66.5 %
  for 68 % on 6 × 6 nm), with a warning.
- Only the outermost wall or shell of multi-wall tubes and onions is oxidized
  or functionalized.
- One ligand species and one attachment template per build.
- Impropers are not written. Bonds across a periodic boundary cannot be used by
  AMBER.
