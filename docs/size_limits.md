# How large a structure the builder handles

Measured 2026-10-21 with `tools/benchmark_sizes.py` on a 2-core Linux
container (Python 3.13, 7 GB RAM): Au/FCC, a C12 thiol (dodecanethiol,
all-trans) at 2 ligands per nm^2 via Au-thiolate, every export format on.
"Atoms" includes ligands. Build = substrate + ligand placement + validation;
export = writing and re-reading every file plus the ZIP. Peak is the
process's maximum resident memory. Rerun on the target machine with
`PYTHONPATH=src python tools/benchmark_sizes.py --out bench.csv`.

## Results

| Shape | Size | Atoms | Ligands | Build (s) | Export (s) | Peak memory | ZIP | Viewer payload |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Sphere | 3 nm | 3,053 | 57 | 0.5 | 0.3 | 154 MB | 0.4 MB | 0.2 MB |
| Sphere | 5 nm | 9,891 | 157 | 1.3 | 0.6 | 188 MB | 1.0 MB | 0.7 MB |
| Sphere | 7 nm | 22,345 | 308 | 2.9 | 1.6 | 285 MB | 2.0 MB | 1.5 MB |
| Sphere | 9 nm | 41,885 | 509 | 4.2 | 2.4 | 419 MB | 3.4 MB | 2.9 MB |
| Sphere | 11 nm | 69,971 | 760 | 7.2 | 3.9 | 634 MB | 5.2 MB | 4.8 MB |
| Sphere | 13 nm | 108,127 | 1,062 | 9.7 | 5.7 | 965 MB | 7.3 MB | 7.3 MB |
| Sphere | 15 nm | 157,733 | 1,414 | 14.5 | 7.8 | 1.37 GB | 10.1 MB | 10.7 MB |
| Rod | 3 x 9 nm | 9,915 | 170 | 1.2 | 0.6 | 188 MB | 1.0 MB | 0.7 MB |
| Rod | 5 x 15 nm | 33,519 | 471 | 3.5 | 2.0 | 365 MB | 3.0 MB | 2.3 MB |
| Rod | 7 x 21 nm | 77,683 | 924 | 8.4 | 4.3 | 702 MB | 6.1 MB | 5.3 MB |
| Slab (111), 6 layers | 5 x 5 nm | 4,604 | 58 | 1.1 | 0.2 | 173 MB | 0.4 MB | 0.3 MB |
| Slab (111), 6 layers | 10 x 10 nm | 20,236 | 254 | 3.9 | 1.2 | 315 MB | 1.6 MB | 1.4 MB |
| Slab (111), 6 layers | 15 x 15 nm | 44,116 | 554 | 8.0 | 2.5 | 515 MB | 3.5 MB | 3.0 MB |
| Slab (111), 6 layers | 20 x 20 nm | 77,166 | 969 | 15.4 | 4.5 | 804 MB | 6.0 MB | 5.3 MB |
| Slab (111), 6 layers | 25 x 25 nm | 119,424 | 1,500 | 24.9 | 7.0 | 1.07 GB | 9.0 MB | 8.3 MB |

Every case placed all requested ligands and passed validation.

## Before and after the scaling fixes

The first benchmark found three problems. The sphere-cut, clash-check and
validation fixes were checked to give identical output; "evenly spaced"
placement now breaks exact ties deterministically, so it can pick different
sites than before for the same seed.

| Case | Build before | Build after | Peak before | Peak after |
|---|---:|---:|---:|---:|
| 15 nm sphere (158k atoms) | 124 s | 14.5 s | 3.15 GB | 1.37 GB |
| 25 x 25 nm slab (119k atoms) | 145 s | 24.9 s | 2.29 GB | 1.07 GB |

1. The sphere cut looked up the element of each surface atom by copying the
   whole symbol list each time (quadratic): 62 of 65 s for an 11 nm sphere.
2. The clash checker rebuilt its KD-tree of placed ligands (nine periodic
   copies on slabs) after every ligand; it now rebuilds every 32 ligands and
   checks the newest ones directly. "Evenly spaced" now updates each site's
   distance to the nearest ligand incrementally.
3. Validation built an all-pairs 3.5 A neighbour list (about 1 GB for 70k
   atoms) although only ligand contacts and near-duplicate substrate atoms
   matter; it now searches just those.

## What limits size now

- **Memory**, about 9 MB per 1,000 atoms at the top end. A 4 GB session
  comfortably fits the 15 nm sphere. Record the nanoHUB session limit in
  `deployment_notes.md`.
- **The browser viewer.** The page carries every atom (10.7 MB for 158k
  atoms). In headless Chromium with *software* rendering, first draw took
  about 6 s at 108k atoms; frame rate on a real GPU was not measured. Above
  6,000 atoms the viewer already switches to smaller spheres.
- **Download links** are base64 data URIs, refused above 40 MB per file; the
  largest ZIP here is 10 MB, so not a limit in this range.
- **PDB format**: above 99,999 atoms serial numbers wrap and CONECT records
  are dropped (the LAMMPS and GSD files keep all bonds).

## Defaults

- `max_atoms` (substrate atoms, estimated before building) now defaults to
  **100,000**, about a 14 nm Au sphere or a 25 x 25 nm, 6-layer Au(111) slab.
  It can be raised to 200,000 under Geometry, "Max atoms (advanced)".
- `max_attempts` now defaults to **50,000** placement trials (was 5,000, which
  capped dense coverage on large surfaces).
