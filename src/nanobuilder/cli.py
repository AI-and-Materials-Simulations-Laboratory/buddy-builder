"""Command-line build from a configuration file (no GUI needed).

    python -m nanobuilder.cli examples/au_sphere_12ligands.json --out build/
    python -m nanobuilder.cli examples/uploaded_crystal_slab.json \\
        --crystal examples/fixtures/mg_hcp_fully_occupied.cif --out build/

A ``config.json`` taken from an exported ZIP works too, which is how a run is
reproduced (pass the ZIP's ``input_crystal.cif`` with --crystal for uploads).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .pipeline import export_package, run_pipeline
from .schema import BuilderError, RunConfig


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m nanobuilder.cli", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="run configuration JSON")
    parser.add_argument("--crystal", help="crystal file for configs with source 'uploaded_file'; "
                        "for 'materials_project', the exported CIF (otherwise the structure is "
                        "fetched live, which needs MP_API_KEY)")
    parser.add_argument("--out", default="build", help="output directory (default: build)")
    args = parser.parse_args(argv)

    data = json.loads(Path(args.config).read_text(encoding="utf-8"))
    upload = None
    if data.get("material", {}).get("source") == "uploaded_file":
        if not args.crystal:
            parser.error("this configuration uses an uploaded crystal; pass it with --crystal")
        upload = Path(args.crystal).read_bytes()
        data["material"]["filename"] = Path(args.crystal).name
        data["material"].pop("sha256", None)
    elif data.get("material", {}).get("source") == "materials_project":
        if args.crystal:
            upload = Path(args.crystal).read_bytes()
        else:
            from .importers import atoms_to_cif_bytes
            from .materials_project import fetch_structure
            try:
                atoms, _ = fetch_structure(data["material"].get("mp_id", ""))
            except BuilderError as exc:
                print(f"Error: {exc}", file=sys.stderr)
                return 1
            upload = atoms_to_cif_bytes(atoms)
        data["material"].pop("sha256", None)
    try:
        config = RunConfig.from_dict(data, uploaded_bytes=upload)
        result = run_pipeline(config, progress=lambda stage, _: print(f"  {stage}..."))
        report = result.report
        print(f"Status: {report.status}; counts: {report.counts}")
        for line in report.errors:
            print(f"  ERROR: {line}")
        for line in report.warnings + result.assembly.warnings:
            print(f"  warning: {line}")
        if not report.ok:
            return 1
        package = export_package(result, Path(args.out))
    except BuilderError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {package.zip_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
