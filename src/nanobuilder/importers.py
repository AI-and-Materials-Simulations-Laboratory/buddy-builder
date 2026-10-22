"""Read user-uploaded crystal files (spec stage A2).

Uploaded files are untrusted: names are sanitized, sizes and atom counts are
capped, and the contents are only ever parsed as crystal data.
"""

from __future__ import annotations

import hashlib
import io
import re
import tempfile
import warnings
from math import gcd
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read

from .schema import BuilderError

MAX_UPLOAD_BYTES = 2_000_000
MAX_CELL_ATOMS = 1000
MIN_ATOM_DISTANCE_A = 0.5
MIN_LATTICE_LENGTH_A = 1.0
OCCUPANCY_TOL = 1e-3

_EXTENSIONS = {"cif": ".cif", "vasp": ".vasp"}
_ASE_FORMATS = {"cif": "cif", "vasp": "vasp"}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sanitize_filename(name: str, file_format: str) -> str:
    """Return a safe basename with the extension that matches ``file_format``."""
    base = Path(str(name or "upload")).name
    stem = re.sub(r"[^A-Za-z0-9._-]", "_", Path(base).stem)[:60].strip("._") or "upload"
    return stem + _EXTENSIONS[file_format]


def read_crystal_bytes(data: bytes, filename: str, file_format: str = "cif"):
    """Parse an uploaded crystal file and run the P0 validity checks.

    Returns ``(atoms, info)`` where ``info`` holds provenance and warnings.
    Raises :class:`BuilderError` with an actionable message on any problem.
    """
    if file_format not in _ASE_FORMATS:
        raise BuilderError(f"Unsupported upload format {file_format!r}; use CIF")
    if not data:
        raise BuilderError("The uploaded file is empty")
    if len(data) > MAX_UPLOAD_BYTES:
        raise BuilderError(f"The uploaded file is {len(data) / 1e6:.1f} MB; the limit is "
                           f"{MAX_UPLOAD_BYTES / 1e6:.0f} MB")

    safe_name = sanitize_filename(filename, file_format)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / safe_name
        path.write_bytes(data)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                kwargs = {"store_tags": True} if file_format == "cif" else {}
                # index=":" so multi-block CIFs are seen; read() alone keeps only the last block.
                blocks = read(path, format=_ASE_FORMATS[file_format], index=":", **kwargs)
        except Exception as exc:  # parsers raise many exception types
            raise BuilderError(f"Could not read the file as {file_format.upper()}: {exc}") from None

    if not isinstance(blocks, list):
        blocks = [blocks]
    if len(blocks) == 0:
        raise BuilderError("No crystal structure found in the file")
    if len(blocks) > 1:
        raise BuilderError(f"CIF holds {len(blocks)} structures (data blocks); upload a file "
                           "with exactly one")
    atoms = blocks[0]
    check_crystal(atoms)
    if file_format == "cif":
        check_ordered(atoms)
        check_formula_consistency(atoms)
    check_min_distance(atoms)

    atoms = Atoms(atoms.get_chemical_symbols(), positions=atoms.get_positions(),
                  cell=atoms.cell, pbc=True)
    atoms.wrap()
    info = {
        "filename": safe_name,
        "original_filename": Path(str(filename or "")).name,
        "format": file_format,
        "sha256": sha256_bytes(data),
        "normalization": "none: cell as given in the file; positions wrapped into the cell",
        "warnings": [],
    }
    return atoms, info


def check_crystal(atoms: Atoms) -> None:
    if len(atoms) == 0:
        raise BuilderError("The file contains no atoms")
    if len(atoms) > MAX_CELL_ATOMS:
        raise BuilderError(f"The unit cell has {len(atoms)} atoms; the limit is {MAX_CELL_ATOMS}")
    if np.any(atoms.numbers <= 0):
        raise BuilderError("The file contains dummy or unknown elements")
    cell = np.array(atoms.cell)
    lengths = np.linalg.norm(cell, axis=1)
    if np.any(lengths < MIN_LATTICE_LENGTH_A) or abs(np.linalg.det(cell)) < 1e-3:
        raise BuilderError("CIF missing a periodic cell: the file needs three nonzero lattice "
                           "vectors with positive volume")
    if not np.all(np.isfinite(atoms.get_positions())):
        raise BuilderError("The file contains non-finite coordinates")


def check_ordered(atoms: Atoms, tol: float = OCCUPANCY_TOL) -> None:
    # ASE keeps one species per site; the true occupancies survive only here.
    for occ in atoms.info.get("occupancy", {}).values():
        if len(occ) > 1 or abs(sum(occ.values()) - 1.0) > tol:
            pretty = ", ".join(f"{el} {frac:g}" for el, frac in occ.items())
            raise BuilderError(f"Disordered or partially occupied site ({pretty}); this version "
                               "accepts only fully ordered crystals")


def check_min_distance(atoms: Atoms, min_A: float = MIN_ATOM_DISTANCE_A) -> None:
    if len(atoms) < 2:
        return
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, np.inf)
    if d.min() < min_A:
        i, j = np.unravel_index(np.argmin(d), d.shape)
        raise BuilderError(f"Atoms {i} and {j} are only {d.min():.2f} A apart; the file may "
                           "contain duplicate sites")


def _parse_formula_sum(text: str) -> dict:
    counts = {}
    for symbol, number in re.findall(r"([A-Z][a-z]?)\s*([0-9.]*)", str(text)):
        counts[symbol] = counts.get(symbol, 0.0) + (float(number) if number else 1.0)
    return counts


def _reduced_ratio(counts: dict) -> dict:
    values = [int(round(v * 1000)) for v in counts.values()]
    g = 0
    for v in values:
        g = gcd(g, v)
    return {k: int(round(v * 1000)) // g for k, v in counts.items()} if g else counts


def check_formula_consistency(atoms: Atoms) -> None:
    declared = atoms.info.get("_chemical_formula_sum")
    if not declared:
        return
    expected = _parse_formula_sum(declared)
    symbols = atoms.get_chemical_symbols()
    actual = {s: float(symbols.count(s)) for s in set(symbols)}
    if not expected or set(expected) != set(actual):
        raise BuilderError(f"Inconsistent formula: the file declares {declared!r} but its atom "
                           f"sites give {atoms.get_chemical_formula()}")
    if _reduced_ratio(expected) != _reduced_ratio(actual):
        raise BuilderError(f"Inconsistent formula: the file declares {declared!r} but its atom "
                           f"sites give {atoms.get_chemical_formula()}")


def atoms_to_cif_bytes(atoms: Atoms) -> bytes:
    """Write a unit cell as CIF (used for ``normalized_unit_cell.cif``)."""
    from ase.io import write

    buf = io.BytesIO()
    write(buf, atoms, format="cif")
    return buf.getvalue()
