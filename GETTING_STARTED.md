# Getting started with Buddy Builder

This guide is for judges, new teammates and anyone trying the tool for the first time. Buddy Builder builds ligand-coated nanoparticles, nanorods, slabs, channels and carbon nanomaterials (graphene, graphene oxide, nanotubes, fullerenes). It exports them as structure and topology files for LAMMPS, HOOMD-blue, AMBER and GROMACS.

You can drive it in three ways:

- the tabs;
- a Python console;
- **Ask Buddy**, a chat that builds and edits the model from plain English.

Setup takes about 5 minutes, and the tour after it takes another 5.

## 1. Install (once)

You need `git` and either **Python 3.12** (3.11 to 3.13 also work) or [uv](https://docs.astral.sh/uv/), which downloads Python for you. The setup uses about 500 MB of disk space.

**macOS / Linux**
```
git clone https://github.com/AI-and-Materials-Simulations-Laboratory/buddy-builder
cd buddy-builder
./setup.sh
```

**Windows** (Command Prompt)
```
git clone https://github.com/AI-and-Materials-Simulations-Laboratory/buddy-builder
cd buddy-builder
setup.bat
```

The setup script does four things:

1. creates a private Python environment in `.venv`, so nothing is installed system-wide;
2. installs the pinned packages from `requirements.txt`;
3. creates a `.env` settings file;
4. runs a self-check, which ends with `OK`.

Run `./setup.sh --tests` to also run the full test suite (about 40 s).

## 2. Start

```
./start.sh        (Windows: start.bat)
```

Your browser opens **http://127.0.0.1:8866** with Buddy Builder as a standalone page. Stop it with Ctrl+C in the terminal.

Options:
- `--port 9000` uses another port;
- `--no-browser` only prints the address;
- `--notebook` opens the classic Jupyter notebook view.

## 3. Turn on the AI (optional)

Ask Buddy works **without any key**: an offline parser understands requests like "a 4 nm gold sphere with 60 dodecanethiol". For free-form requests and questions, add a key to `.env`:

```
BUDDY_LLM_API_KEY=sk-...        # an OpenAI key
```

Restart `./start.sh`. The badge next to *Ask Buddy* then names the model instead of "Offline parser". To use Claude instead, set `BUDDY_LLM_PROVIDER=anthropic` and use an Anthropic key. `.env` is never committed to git.

## 4. A five-minute tour

1. **Click "Generate model".** The default build, a 3 nm gold sphere with 12 ligands, appears in the 3D viewer. Drag to rotate, scroll to zoom, and click an atom to inspect it.
2. **Chat with Buddy.** In *Ask Buddy*, type:
   `a 4 nm gold sphere with 60 dodecanethiol in fibonacci order, export LAMMPS and PDB`
   Buddy shows the commands it ran and a ✓ list comparing the result with your request, then builds the model.
3. **Change the model by chatting** (simple edits are instant):
   `make it 5 nm` → `10 more ligands` → `highlight the sulfur atoms` → `remove ligand 3` → `undo`
4. **Try another material:**
   `a 3nm graphene oxide sheet`, then `add 10 phenyl groups`
   Buddy adds only what you ask for. Anything you leave out is listed under "Buddy assumed".
5. **Ask a question** (needs an AI key): `why are there hydrogens on the graphene edges?`, or press **Explain build**.
6. **Use the Python console** (bottom right):
   - type `silver sphere 3 nm with 15 MUA` and press **AI**: the box fills with console commands that were already tested;
   - press **Test** to dry-run them on a hidden copy;
   - press **Run** to run them.
7. **Download** the ZIP, the LAMMPS/PDB/MOL2 files, the **Certificate of Build (PDF)** and **Save session**, a replayable script of every step.

The **Try an example...** menu in the chat and **Load an example...** under the tabs have more ready-made builds.

## How the AI is kept honest

- The AI never writes coordinates or code. It writes only Buddy's own console commands, which are parsed and never executed as Python.
- Every plan is first run on a **hidden copy** of the app: Buddy test-builds it, validates it and checks it against your request. If that fails, the plan is repaired, at most twice.
- Commands run through the visible console, so the console log, the tabs and the saved session show exactly what the AI did.

Details are in [docs/ai_assistant.md](docs/ai_assistant.md).

## Troubleshooting

| Problem | Fix |
|---|---|
| The page stops responding after the server was restarted | Reload the page; it connects to the new session. |
| "Address already in use" | `./start.sh --port 9000`, or close the other copy. |
| The chat badge says "Offline parser" | Add `BUDDY_LLM_API_KEY` to `.env` and restart. |
| `setup.sh` says that Python is missing | Install Python 3.12 or uv (links above), then run it again. |
| Packages look broken | Delete `.venv` and run the setup again. |
| A build fails | Buddy explains why and how to fix it; ask "why did my build fail?" in the chat. |

## For the team

| Path | What it is |
|---|---|
| `src/nanobuilder/` | The builder: crystal, shapes, ligands, placement, exports and the GUI (`gui.py`). |
| `src/nanobuilder/assistant/` | Ask Buddy: `provider.py` (the AI connection), `knowledge.py` (retrieval), `rules.py` (the offline parser), `agents.py` (plan, verify, repair), `panel.py` (chat). |
| `src/nanobuilder/console.py` | The Python console commands; the AI uses the same vocabulary. |
| `tests/` | `.venv/bin/python -m pytest -q`; no network or key needed. |
| `app.ipynb` | The notebook entry point (nanoHUB AppMode); `BuddyBuilder.ipynb` is the app view that `start.sh` serves. |
| `docs/` | Design notes and scientific assumptions; [docs/AI_PROGRESS.md](docs/AI_PROGRESS.md) is the log of the AI work. |

The tool exports **structure and topology only**: no force field or charges are generated, and every export says so.
