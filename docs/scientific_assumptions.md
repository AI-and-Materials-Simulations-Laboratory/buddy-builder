# Scientific assumptions and limitations

Everything the tool produces is a **starting configuration**, not an
equilibrated or energy-minimized structure.

## Materials
- Preset lattice constants (Angstrom): Au 4.08, Ag 4.09, Cu 3.61, Fe 2.87,
  Si 5.43, Ge 5.66, Mg a = 3.21 / c = 5.21 (ASE reference states, which
  `ase.build.bulk()` uses), NaCl 5.64 (common room-temperature literature
  value; verify against a primary source before scientific use). They are not
  relaxed or temperature-specific values.
- Presets use conventional cells (`cubic=True` for FCC, BCC, diamond and
  rocksalt; the hexagonal cell for HCP), so Miller indices have their usual
  meaning. On a primitive FCC cell, `(1,0,0)` would give a {111} surface.
- A user-edited lattice constant is stored as `lattice_overrides_A`, warned
  about beyond 5% and rejected beyond 25% of the registry value.
- Uploaded crystals must be fully ordered. Mixed or partially occupied sites,
  several data blocks, atoms closer than 0.5 A and declared formulas that do
  not match the sites are rejected. The cell is used as given; Miller indices
  refer to it (`normalized_unit_cell.cif` in each export).
- Fe's magnetic state, surface reconstructions, relaxation and passivation are
  outside the model.

- **Materials Project** structures are fetched live (`/materials/summary/`,
  user's API key) and used exactly as served (DFT-relaxed, 0 K, ordered). They
  are converted to spglib's conventional standard cell so Miller indices and rod
  axes have their usual meaning, then stored as a CIF in the config and the ZIP
  so the build is reproducible offline. Alloys are only what MP lists as
  ordered compounds (e.g. L1₂ AuCu3, L1₀ FePt, B2 NiTi); random
  solid solutions are not generated. MP data are CC BY 4.0; cite Jain et al.,
  APL Mater. 1, 011002 (2013). Au-thiolate stays limited to the built-in Au
  preset, so MP crystals use geometric adsorption.

## Geometry
- Spheres are cut around the atom nearest the supercell centre, keeping atoms
  within the nominal radius. Coverage uses the nominal area 4 pi R^2.
- Nanorods are a cylinder of diameter d along a lattice direction [u v w] of
  the conventional cell, with hemispherical caps (a spherocylinder; length is
  tip to tip) or flat caps, cut around the atom nearest the supercell centre
  and rotated so the axis is +z. Surface atoms are split into "side" and
  "tips" regions; coverage uses the nominal areas (side 2 pi R L_cyl, tips
  4 pi R^2 or 2 pi R^2). Ligand normals are radial from the axis on the side
  and from the cap centre (or along the axis for flat caps) at the tips. Real
  Au nanorods are often penta-twinned with {100}/{110} side facets; this
  single-crystal cut does not reproduce twinning or faceting.
- Slabs: one layer is one atomic plane. Distinct terminations are listed and
  must be chosen explicitly when there is more than one; none is assumed
  neutral or stable. Coverage uses the actual periodic face area.
- Cuts of compounds can be non-stoichiometric (net formal charge for ionic
  crystals) and covalent cuts leave dangling bonds; both are warned about and
  left as they are.
- Surface atoms have fewer first-shell neighbours than the same species in the
  bulk. First-shell cutoffs are per species pair, halfway between the first and
  second distance shells.

## Ligands and attachment
- SMILES gives connectivity, not surface chemistry; the user chooses the anchor
  atom and the attachment template.
- Au-thiolate is offered only for the built-in Au/FCC preset with a thiol (-SH)
  anchor. The thiol hydrogen is removed after the ligand is embedded and
  optimized (ETKDGv3, then MMFF94 for the isolated thiol), and an explicit Au-S
  bond is added with the S 2.4 A above a surface Au atom (on-top). That bond is
  topology only: it carries no Au-S potential, and gas-phase MMFF charges are
  not transferred.
- Every other substrate (other presets and all uploads, including Au CIFs)
  gets geometric adsorption: the intact ligand is placed at the chosen distance
  (default 2.5 A) with no surface bond.
- Ligands stand along the outward normal (anchor-to-centroid of the other heavy
  atoms), rotated about the normal in 12 steps to avoid clashes.
- "Extended chains" (the default for library chain ligands) sets every acyclic
  single-bond backbone dihedral, walking out from the anchor, to 180 deg
  (all-trans), then relaxes briefly with MMFF/UFF. Double bonds (e.g. the cis
  bond of oleyl chains) and rings are untouched. This is the usual SAM/brush
  starting geometry; real chains are tilted (about 30 deg for alkanethiols on
  Au(111)) and partly gauche, and PEG prefers a helical conformation. Because
  chains stand along the normal and are not tilted, the clash-free packing
  limit is below a real SAM (about 3.2 vs 4.6 C12 thiols per nm^2 on Au(111)).
- The ligand library is a set of starting inputs (SMILES, anchor, template,
  typical anchor distance). Carboxylates (citrate, oleate) and P/O ligands are
  placed through one anchor atom by geometric adsorption; multidentate binding
  is not modeled.

## Coverage and clashes
- Target = the requested count, or round(density x area). Requested and
  achieved values are reported separately; a shortfall is never hidden.
- A trial ligand is rejected if any nonbonded pair is closer than
  alpha x (vdW_i + vdW_j) (alpha = 0.7 by default; Bondi radii, covalent radius
  + 0.8 A where Bondi has none). The anchor-to-site pair is exempt.
  Ligand-substrate and ligand-ligand contacts are checked across in-plane
  periodic images.
- A site counts as accessible when no other substrate atom is closer to the
  anchor position than 0.9 x the anchor distance.
- A non-overlapping structure is not an energy minimum.

## Carbon nanomaterials

Full reference: [carbon_materials.md](carbon_materials.md).

- Frameworks are ideal: a_CC = 1.42 Å throughout (graphene, tube walls,
  relaxed fullerenes); layer and wall spacings are user inputs (3.35 / 3.4 Å
  defaults); fullerenes and onions are icosahedral Goldberg shells only.
- Connectivity is an explicit bond graph created by the builder; nothing is
  perceived from distances, so bond orders and types follow how each atom was
  made (aromatic 1.5 for intact sp2 C-C, 1.0 next to sp3 carbons).
- Oxygen content is 100·N_O/N_C over the whole substrate (HierGO "percent
  oxygen", AMBER GO dataset "overall modification"). Oxidation is limited to
  outer faces (sheets), the outer wall (tubes) or the outer shell (onions).
- Group geometry uses typical bond lengths (C-OH 1.43, O-H 0.97, C-O epoxide
  with O 1.276 Å above the C-C midpoint, C=O 1.23 Å) and 0.25 Å sp3 puckering;
  sheets are otherwise flat. Real GO wrinkles; minimize before production MD.
- Placement rules (no same-face groups on neighbouring carbons, minimum
  non-bonded distances, chair-sublattice bias above about 25 % O) are
  heuristics that reproduce reported GO compositions; they are not a kinetic
  model of oxidation.
- Covalent attachment edits the substrate (sp2 to sp3 carbon, removed H or OH,
  opened epoxide) and bonds the ligand at a fixed length per template and
  anchor element. The site carbon is puckered 0.40 Å; no reaction energetics
  or selectivity is modeled, and sites are chosen at random (or evenly spaced)
  among chemically valid ones.
- Suggested GAFF2 types: exact by construction for the substrate, provisional
  (rule-based) for ligand atoms. OPLS-AA names follow HierGO's atom-typing
  script for the substrate only.
