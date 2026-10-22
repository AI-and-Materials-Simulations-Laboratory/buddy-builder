"""Measure how large a nanoparticle / nanorod / slab the builder handles.

For each size: atom count, wall time per stage (build, ligand placement,
validation, export), peak Python memory, ZIP size and the size of the viewer
payload the browser has to receive. Writes a CSV and prints a table.

    PYTHONPATH=src python tools/benchmark_sizes.py --out bench.csv

Each case runs in a fresh subprocess so memory numbers do not accumulate.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time

CASES = (
    # (label, geometry, coverage density per nm^2)
    [("sphere", {"shape": "sphere", "diameter_nm": d}, 2.0) for d in (3, 5, 7, 9, 11, 13, 15)]
    + [("rod", {"shape": "rod", "diameter_nm": d, "rod_length_nm": 3 * d}, 2.0) for d in (3, 5, 7)]
    + [("slab", {"shape": "slab", "miller": [1, 1, 1], "slab_length_nm": s, "slab_width_nm": s,
                 "thickness_mode": "layers", "layers": 6}, 2.0) for s in (5, 10, 15, 20, 25)]
)

CHILD = r"""
import json, sys, time, tempfile, resource
sys.path.insert(0, "src")
from nanobuilder.schema import RunConfig
from nanobuilder.pipeline import run_pipeline, export_package
from nanobuilder.viewer import assembly_payload, iframe_html
geometry, density, max_atoms = json.loads(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
geometry.update(vacuum_A=20.0, max_atoms=max_atoms)
cfg = RunConfig.from_dict({
    "material": {"source": "preset", "preset_id": "Au_fcc"},
    "geometry": geometry,
    "ligand": {"smiles": "SCCCCCCCCCCCC", "anchor_atom_index": 0, "binding_template": "Au_thiolate",
               "conformer_mode": "extended"},
    "coverage": {"mode": "density", "density_nm2": density, "region": "entire_surface"
                 if geometry["shape"] != "slab" else "top", "max_attempts": 10**7}})
t0 = time.perf_counter()
out = {}
try:
    r = run_pipeline(cfg)
    t1 = time.perf_counter()
    pkg = export_package(r, tempfile.mkdtemp())
    t2 = time.perf_counter()
    html = iframe_html(assembly_payload(r.assembly))
    t3 = time.perf_counter()
    a = r.assembly
    out = dict(ok=True, atoms=a.n_atoms, substrate=a.n_substrate, ligands=a.achieved_ligands,
               target=a.target_ligands, status=r.report.status,
               build_s=round(t1 - t0, 2), export_s=round(t2 - t1, 2), viewer_s=round(t3 - t2, 2),
               timings={k: round(v, 2) for k, v in r.timings.items()},
               zip_mb=round(pkg.zip_path.stat().st_size / 1e6, 2),
               zip_link_mb=round(pkg.zip_path.stat().st_size * 4 / 3 / 1e6, 2),
               viewer_mb=round(len(html) / 1e6, 2))
except Exception as exc:
    out = dict(ok=False, error=str(exc)[:200], build_s=round(time.perf_counter() - t0, 2))
out["peak_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)  # Linux: kB
print("RESULT " + json.dumps(out))
"""


def run_case(geometry, density, max_atoms, timeout):
    try:
        p = subprocess.run([sys.executable, "-c", CHILD, json.dumps(geometry), str(density),
                            str(max_atoms)], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"timeout > {timeout} s"}
    for line in p.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    return {"ok": False, "error": (p.stderr or p.stdout)[-300:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="bench.csv")
    ap.add_argument("--max-atoms", type=int, default=400000)
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()
    rows = []
    for label, geometry, density in CASES:
        t = time.time()
        res = run_case(dict(geometry), density, args.max_atoms, args.timeout)
        size = (f"{geometry['diameter_nm']} nm" if label == "sphere" else
                f"{geometry['diameter_nm']}x{geometry['rod_length_nm']} nm" if label == "rod" else
                f"{geometry['slab_length_nm']}x{geometry['slab_width_nm']} nm, 6 layers")
        row = {"shape": label, "size": size, **{k: v for k, v in res.items() if k != "timings"},
               **{f"t_{k}": v for k, v in (res.get("timings") or {}).items()}}
        rows.append(row)
        print(f"{label:6s} {size:24s} " + (
            f"atoms {res['atoms']:>7}  ligands {res['ligands']:>5}/{res['target']:<5}  "
            f"build {res['build_s']:>7.1f}s  export {res['export_s']:>6.1f}s  peak {res['peak_mb']:>7.0f} MB  "
            f"zip {res['zip_mb']:>6.1f} MB  viewer {res['viewer_mb']:>6.1f} MB"
            if res.get("ok") else f"FAILED: {res.get('error')}") + f"   ({time.time() - t:.0f}s wall)",
              flush=True)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("shape", "size"), k))
    with open(args.out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
