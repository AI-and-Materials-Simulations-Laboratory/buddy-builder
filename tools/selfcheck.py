"""Quick check after setup: imports, one small build, and Ask Buddy (offline, no network).

Run by setup.sh / setup.bat; on its own:  PYTHONPATH=src .venv/bin/python tools/selfcheck.py
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
os.environ["BUDDY_DOTENV"] = "0"  # never read keys or call a model here


def main() -> int:
    t0 = time.time()
    try:
        import ase  # noqa: F401
        import ipywidgets  # noqa: F401
        import numpy  # noqa: F401
        import rdkit  # noqa: F401
        from nanobuilder.pipeline import run_pipeline
        from nanobuilder.schema import RunConfig
    except ImportError as exc:
        print(f"FAILED: missing package ({exc}). Run the setup script again.")
        return 1
    result = run_pipeline(RunConfig())
    a = result.assembly
    print(f"  build: {a.n_atoms} atoms, {a.achieved_ligands}/{a.target_ligands} ligands, "
          f"validation {result.report.status}")
    from nanobuilder.assistant import Assistant
    from nanobuilder.assistant.provider import DisabledProvider
    from nanobuilder.gui import BuilderApp
    app = BuilderApp(assistant=False)
    plan = Assistant(provider=DisabledProvider("self-check")).ask(
        "a 3 nm gold sphere with 10 dodecanethiol", app)
    print(f"  Ask Buddy (offline parser): {'ok' if plan.ok else 'problems: ' + '; '.join(plan.problems)}"
          f" -> {'; '.join(plan.commands)}")
    ok = result.report.ok and plan.ok
    print(("OK" if ok else "FAILED") + f" ({time.time() - t0:.1f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
