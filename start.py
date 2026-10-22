#!/usr/bin/env python3
"""Start Buddy Builder as a standalone app in your browser (no notebook code cells).

    python start.py               # opens http://127.0.0.1:8866
    python start.py --port 9000   # another port
    python start.py --no-browser  # just print the address
    python start.py --notebook    # classic Jupyter Notebook instead (shows the code cell)

It uses this folder's .venv when there is one, serves BuddyBuilder.ipynb (the code cell of app.ipynb) with Voila
(only the widgets are shown) and falls back to Jupyter Notebook when Voila is
not installed. Ask Buddy reads its AI key from .env here (see
docs/ai_assistant.md). Press Ctrl+C to stop.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import socket
import subprocess
import sys
import threading
import time
import webbrowser
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE / "app.ipynb"  # notebook view (instructions + code cell)
APP = HERE / "BuddyBuilder.ipynb"  # app view: the same code cell only, served by Voila


def venv_python():
    for rel in (".venv/bin/python", ".venv/Scripts/python.exe"):
        p = HERE / rel
        if p.exists():
            return p
    return None


def free_port(start: int) -> int:
    for port in range(start, start + 50):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise SystemExit(f"No free port between {start} and {start + 49}")


def open_when_ready(url: str, port: int, timeout: float = 90.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                time.sleep(0.8)
                webbrowser.open(url)
                return
        time.sleep(0.4)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Start Buddy Builder in the browser")
    ap.add_argument("--port", type=int, default=8866)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--notebook", action="store_true", help="Jupyter Notebook instead of the app view")
    args = ap.parse_args(argv)

    # Run with this folder's virtual environment (it has the builder's dependencies).
    py = venv_python()
    if py is not None and Path(sys.executable).resolve() != py.resolve() and \
            os.environ.get("BUDDY_START_REEXEC") != "1":
        env = dict(os.environ, BUDDY_START_REEXEC="1")
        return subprocess.call([str(py), str(Path(__file__).resolve())] + (argv or sys.argv[1:]), env=env)

    missing = [m for m in ("numpy", "ase", "rdkit", "ipywidgets") if importlib.util.find_spec(m) is None]
    if missing:
        print(f"Missing packages for {sys.executable}: {', '.join(missing)}\n"
              "Install them first:  python -m pip install -r requirements.txt", file=sys.stderr)
        return 1

    port = free_port(args.port)
    use_voila = not args.notebook and importlib.util.find_spec("voila") is not None
    if use_voila:
        url = f"http://127.0.0.1:{port}/"
        cmd = [sys.executable, "-m", "voila", str(APP.name), f"--port={port}", "--no-browser",
               "--Voila.ip=127.0.0.1", "--show_tracebacks=True"]
    else:
        if not args.notebook:
            print("Voila is not installed (python -m pip install voila), so Buddy Builder opens in "
                  "Jupyter Notebook: run the first cell.")
        token = os.urandom(12).hex()
        url = f"http://127.0.0.1:{port}/notebooks/app.ipynb?token={token}"
        cmd = [sys.executable, "-m", "notebook", str(NOTEBOOK.name), f"--port={port}", "--no-browser",
               "--ip=127.0.0.1", f"--IdentityProvider.token={token}"]

    print(f"Buddy Builder: {url}\n(Ctrl+C to stop)", flush=True)
    if not args.no_browser:
        threading.Thread(target=open_when_ready, args=(url, port), daemon=True).start()
    proc = subprocess.Popen(cmd, cwd=str(HERE))
    try:
        return proc.wait()
    except KeyboardInterrupt:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        print("\nBuddy Builder stopped.")
        return 0


if __name__ == "__main__":
    sys.exit(main())
