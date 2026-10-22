# Installation and running

Most people only need the two scripts described in [GETTING_STARTED.md](../GETTING_STARTED.md): `./setup.sh` then `./start.sh` (Windows: `setup.bat`, `start.bat`). This page covers manual installs, conda, Jupyter, nanoHUB and troubleshooting.

## Requirements

- **Python 3.12.** It was tested with 3.12.10 on Windows 11 and with 3.12 on macOS. 3.13 should also work. The pinned NumPy 2.5 and SciPy 1.18 need Python 3.12 or newer.
- About 500 MB of disk space for the environment.
- A browser with internet access for the 3D viewer, which loads 3Dmol.js from cdnjs. The Python side works offline.

Every dependency installs from PyPI with pip, and no compiler is needed. `requirements.txt` pins the tested versions:

- the core packages: `numpy`, `scipy`, `ase`, `rdkit`, `gsd`, `ipywidgets` and `spglib`;
- for local use only: `notebook`, `voila` and `pytest`.

## Manual install

**macOS / Linux**

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

**Windows (PowerShell).** First check which Python a plain `python` runs, with `python --version`. When Anaconda is installed, that is often Anaconda's Python, which can be too old. Create the environment with an explicit 3.12 interpreter:

```powershell
py -3.12 -m venv .venv            # Python from python.org (has the py launcher)
python3.12 -m venv .venv          # or: Python from the Microsoft Store
.venv\Scripts\python -m pip install -r requirements.txt
```

Call `.venv\Scripts\python` for every command. This avoids activation scripts and works even while a conda environment is active.

> **If `.venv` already exists, do not re-run `-m venv` with a different Python.** Doing so re-points the folder at the new interpreter but keeps the packages compiled for the old one, so imports fail. Delete `.venv` first, then create it again.

**Anaconda / conda** (tested on macOS, Apple silicon)

```bash
conda create -n buddy python=3.12 -y
conda activate buddy
python -m pip install -r requirements.txt
```

`pip` inside the activated environment installs into that environment, not into Anaconda's `base`. Activate it in every new terminal, and use `python` wherever these docs say `.venv/bin/python`.

Installing the package itself is optional (`python -m pip install -e .`): the notebook, the launcher and the tests find it in `src/` without installation.

## Running

### As an app (recommended)

```bash
./start.sh                 # macOS / Linux
start.bat                  # Windows
python start.py --help     # --port, --no-browser, --notebook
```

`start.py` uses this folder's `.venv` and serves `BuddyBuilder.ipynb` with Voila. That notebook is the code cell of `app.ipynb`, so the page shows only the app at http://127.0.0.1:8866. If Voila is not installed, the launcher falls back to Jupyter Notebook.

### In Jupyter

Start Jupyter **from the project environment**, so that the kernel has the dependencies:

```bash
.venv/bin/python -m notebook app.ipynb        # Windows: .venv\Scripts\python -m notebook app.ipynb
```

Use `-m notebook`, not `-m jupyter notebook`. The `jupyter` command looks up its `notebook` subcommand on `PATH`, so when Anaconda is on `PATH` it starts Anaconda's Jupyter, whose kernel lacks the dependencies. For the same reason, do not open the notebook from Anaconda Navigator.

Run the cell with Shift+Enter, and the builder appears below it.

### On nanoHUB

`middleware/invoke` launches `app.ipynb` in AppMode, which runs the cell and hides the code. `src/Makefile` provides the build and install hooks. See [deployment_notes.md](deployment_notes.md).

### From the command line

```bash
export PYTHONPATH=src            # PowerShell: $env:PYTHONPATH = "src"
.venv/bin/python -m nanobuilder.cli examples/au_sphere_12ligands.json --out build
.venv/bin/python -m nanobuilder.cli examples/graphene_oxide_flake.json --out build
.venv/bin/python -m nanobuilder.cli examples/uploaded_crystal_slab.json \
    --crystal examples/fixtures/mg_hcp_fully_occupied.cif --out build
```

A `config.json` from any exported ZIP rebuilds the same model, with the same seed and the same atoms.

## AI settings (Ask Buddy)

Paste a key in the app under **Ask Buddy → AI key**, or put it in `.env`; `.env.example` lists every option. Without a key, the offline parser handles common requests. See [ai_assistant.md](ai_assistant.md).

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| The page stops responding after the server restarted | The tab is attached to the old session ("Kernel does not exist") | Reload the page |
| `ModuleNotFoundError: No module named 'ase'` in the notebook | Anaconda's Jupyter is running, not the project's | Stop it, then run `.venv/bin/python -m notebook app.ipynb` |
| A NumPy error naming `cp312` and `cpython-311`, or `DLL load failed while importing _zmq` | `.venv` was re-created over itself with a different Python | Delete `.venv` and install again with Python 3.12 |
| `No matching distribution found for numpy==2.5.4` | Python older than 3.12 | Install Python 3.12 (or uv) and set up again |
| The viewer shows "3Dmol.js could not be loaded" | The browser cannot reach cdnjs.cloudflare.com | Allow internet access for the browser |
| "Address already in use" | Another copy is running | `./start.sh --port 9000`, or stop the other copy |
| New options are missing after updating the code | The kernel still has the old code loaded | Restart the app (or **Kernel → Restart** in Jupyter) |
| Imports fail on macOS in an iCloud-synced folder | macOS hid the venv's `.pth` files | `./setup.sh` fixes it; manually: `chflags nohidden .venv/lib/python3*/site-packages/*.pth` |
