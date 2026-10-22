# MD export notes

All exports come from one canonical `AssemblyModel`, so every file describes the
same atoms in the same order: substrate first, then ligand 1, ligand 2, and so
on. Units are Angstrom and amu throughout.

## Readiness

Exports are **`topology_ready`** (Level B): MD engines can read them, but no
force field, charges or substrate-ligand parameters are included.
`manifest.json` → `parameter_coverage` lists every structural type, bond,
angle and dihedral type that still needs parameters.

## Indexing and identifiers
- Internal indices are zero-based. LAMMPS atom, bond, angle and dihedral IDs
  are index + 1. GSD bond, angle and dihedral groups use the zero-based indices.
- Molecule ID 0 is the substrate; ligands are 1..n. The same IDs appear in the
  LAMMPS `Atoms` section, `topology.json`, the extXYZ `molecule_id` column and
  GSD `log/particles/molecule_id`.
- Type labels are structural, not force-field types: `Au_sub` (substrate
  atom), `S_anchor` (anchor atom), plain elements for other ligand atoms. Bond
  types look like `C-S`, `C-H` or `Au-S`; angle and dihedral types are element
  chains such as `C-C-O`.

## Boxes
- The cell is lower-triangular (a along x, b in the xy plane). Slabs keep their
  periodic in-plane vectors, with tilt reduced so abs(xy) <= lx/2; spheres get
  an orthogonal box.
- Along each nonperiodic axis the box leaves `vacuum_A` (default 30 A) between
  the outermost atom, ligands included, and the box face.
- LAMMPS: restricted triclinic, with the `xy xz yz` line only when a tilt is
  nonzero. Recommended `boundary` is `p p f` for slabs and `f f f` for spheres
  (written into `in.lammps.template`; it belongs in the input script, not the
  data file).
- HOOMD: `[Lx, Ly, Lz, xy, xz, yz]` with dimensionless tilts; positions are
  shifted so the box is centred on the origin, then wrapped. HOOMD boxes are
  always periodic, so pair cutoffs must stay below the vacuum; a vacuum box is
  not the same as a nonperiodic boundary for every electrostatics method.

## Molecules across periodic boundaries (slabs)
- The canonical model keeps every ligand **whole**: each ligand hangs from its
  surface atom, and every bond (including the attachment) is its
  minimum-image vector. A ligand near the cell edge therefore pokes out of the
  cell instead of being cut in two. (Before this fix, atoms were wrapped one
  by one, which split ligands and drew bonds across the whole slab.)
- LAMMPS `Atoms` lines carry wrapped coordinates plus `ix iy iz` image flags;
  GSD stores wrapped positions plus `particles/image`. In both,
  `position + image . box` is the whole molecule.
- PDB, extXYZ and the viewer use the whole-molecule coordinates directly.
- The substrate is periodic and is never "made whole"; its atoms stay inside
  the cell.

## Files
| File | Notes |
|---|---|
| `hoomd/structure.gsd` | Written with `gsd.hoomd.Frame`. Particle charges are not set. `run_hoomd_template.py` raises until a force field is defined. |
| `lammps/structure.data` | `atom_style full`; `q = 0` placeholders only satisfy the syntax (`partial_charges_unassigned=true`). Type names appear as comments in `Masses`. |
| `lammps/in.lammps.template` | `units real`, `atom_style full`, boundary, `read_data`, every type listed, REQUIRED markers for each missing interaction; an optional `setforce` fix when "Freeze substrate" is on. Not runnable as-is. |
| `structure.extxyz` | Cell, pbc, masses and per-atom `molecule_id`, `component`, `chemical_role`, `type_label`. No bonds. |
| `pdb/structure.pdb` | For AMBER (tleap) and GROMACS (editconf/pdb2gmx). `CRYST1` cell; substrate `HETATM` in chain A, one residue per atom named after the element; each ligand one residue `LIG` in chain B with unique atom names (`S1`, `C1`, `H1`...); `CONECT` for every bond. Above 99,999 atoms serials wrap and `CONECT` is omitted (use LAMMPS/GSD for bonds). Force-field residue templates are the user's job. |
| `topology.json` | The authoritative connectivity: atoms, bonds, angles, dihedrals, `site_to_ligand`, charge status, parameter coverage. |

## Checks before an export is offered
`validation.py` blocks exports on malformed bonds or tuples, non-finite
coordinates, missing masses, atoms outside the box along nonperiodic axes,
implausible bond lengths, overlapping atoms (< 0.7 A) and duplicated substrate
atoms (< 0.5 A). Each exporter then re-reads its file and compares counts,
types, positions, box and topology indices; any mismatch blocks the ZIP.

## Validating with real engines
- LAMMPS: point the environment variable `LAMMPS_EXE` at a working `lmp`
  executable and run `python -m pytest tests/test_export_lammps.py`; the
  opt-in test runs `read_data` on a slab export. A `run 0` additionally needs
  a force field.
- HOOMD-blue (Linux/macOS): `hoomd.Simulation(...).create_state_from_gsd("structure.gsd")`
  should load the file; the template shows where the force field goes.

## Carbon nanomaterials and MOL2

- The substrate of a carbon build is a covalent molecule (molecule ID 0) with
  every C-C, C-O, O-H and C-H bond, plus angles and dihedrals enumerated from
  the full bond graph (substrate, ligands and the links between them).
  Impropers are not enumerated.
- Type labels are structural roles (`C_sp2`, `C_edge`, `H_edge`,
  `C_hydroxyl`, `O_epoxy`, `C_linked`, ...). Suggested GAFF2 and HierGO
  OPLS-AA names are in `topology.json` → `suggested_atom_types`; bond orders
  in `bond_orders`; residues in `residues`.
- `mol2/structure.mol2` (Tripos): GAFF2 names in the atom-type column, bond
  codes `1`/`ar`/`2`/`am`, residues GRA (framework C with its edge H), OH,
  EPX, ETH, CBO, CBX, LIG, or a single GRO substrate residue; charges are
  0.0000 (`NO_CHARGES`). `tleap_template.in` is a non-runnable sketch.
- PDB uses the same residue names and per-residue atom names.
- Periodic sheets: LAMMPS boundary `p p f`; infinite tubes `f f p`. Bonds
  across the periodic boundary are kept (LAMMPS, GSD, topology.json); MOL2 and
  PDB write them as ordinary bonds, which AMBER cannot use, so build finite
  flakes or tubes for AMBER.
- Extra validation for covalent models: maximum valence per element, no
  framework carbon with fewer than two bonds, and no non-bonded pair (beyond
  1-3) closer than 0.9 Å.
