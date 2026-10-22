# nanoHUB deployment notes

**Status (2026-10-21): not started on nanoHUB.** Everything below still has to
be done by a team member with a nanoHUB account. The tool runs locally (see
the README). The launcher and Makefile follow the spec's illustrative pattern
and have not been tested on nanoHUB.

## Record here once known

| Item | Value |
|---|---|
| Tool alias | `xtalligand` (proposed; subject to availability) |
| Registered by | |
| Repository checkout URL (copy from the tool page) | |
| Workspace location | |
| Hosted Python version | |
| Hosted package versions (ase, rdkit, gsd, ipywidgets, numpy, scipy) | |
| Tool status (Registered / Created / Uploaded / Installed / Approved / Published) | |
| SVN revision of the demo | |

## Steps (spec section 13)

1. Register the tool on nanoHUB: choose **Host subversion repository here**
   and **Web application (Jupyter, ...)**; add the four developers.
2. Check out the repository given on the tool page (do not guess the URL):
   `svn checkout <URL from the tool page> xtalligand`.
3. Copy this folder's contents into the working copy, **excluding** `.venv/`,
   `.pytest_cache/`, `__pycache__/`, `lib/` and any `build/` output. Set the ignore
   property once:
   `svn propset svn:ignore -F svnignore.txt .`, where `svnignore.txt` lists
   `.venv`, `.pytest_cache`, `build`, `lib`.
4. `chmod 755 middleware/invoke`, then `svn add --force src tests examples docs middleware app.ipynb README.md requirements.txt pyproject.toml LICENSE`
   and commit.
5. In a nanoHUB workspace, check which packages the Jupyter environment has:
   `python -c "import ase, rdkit, gsd, ipywidgets, numpy, scipy, spglib; print('ok')"`.
   spglib is optional: without it, Materials Project structures fall back to
   their primitive cell and a warning says Miller indices refer to that cell.
   nanoHUB does **not** install `requirements.txt` automatically; ask the tool
   administrators for anything missing. Record the ipywidgets version: the
   upload widget is handled for both 7 and 8, but the GUI was only run with 8.
6. Smoke test in the workspace (spec H6):
   ```bash
   cd xtalligand/src && make all && make install && cd .. && ./middleware/invoke -T "$PWD"
   ```
   `make install` copies `src/nanobuilder` to `lib/python/`, which `app.ipynb`
   puts first on `sys.path`. The GUI footer shows the folder the package was
   loaded from, which tells you whether a stale copy is in use.
7. Check in the hosted session: a dropdown change, Generate, rotating the 3D
   view (3Dmol.js comes from cdnjs, so the user's browser needs internet
   access), downloading the ZIP, and example E (uploaded CIF).
8. Verify `middleware/invoke`: `-u anaconda-6` and the `start_jupyter -A`
   AppMode flag are illustrative and must match the hosted environment.
9. Move through Uploaded → Installed → Approved → Published with the nanoHUB
   tool administrators. Choose a license first (see `LICENSE`).

If SVN or Jupyter registration is unavailable, do not work around platform
controls: demo from a local Jupyter session and hand over a ZIP of this folder
plus these notes (spec H8).


## Materials Project from inside nanoHUB (check early)

The live material search calls `https://api.materialsproject.org` **from the
server running the tool**, not from the user's browser. Check in a nanoHUB
workspace terminal, with your own key:

```bash
curl -s -o /dev/null -w "%{http_code}\n" -H "x-api-key: $MP_API_KEY" \
  "https://api.materialsproject.org/materials/summary/?formula=Au&_fields=material_id&_limit=1"
```

`200` means it works. `000` or a timeout means outbound access is blocked: ask
the nanoHUB tool administrators to allow that host. Until then the tool still
works with the built-in presets and uploaded CIFs, and the Material tab
explains the failure instead of hanging. Users paste their own key in the GUI;
never commit a key to the repository. The MP client was written against MP's
documented API and tested with recorded-format responses; the first real
query must be checked by someone with a key (see the implementation log).

## Memory and size

Peak memory grows roughly linearly with atom count; see `size_limits.md`.
Record the session memory limit of the hosted tool here once known, and lower
the default `max_atoms` in `schema.py` if sessions are small.
