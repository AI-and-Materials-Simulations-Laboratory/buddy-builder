# Contributing to Buddy Builder

Thanks for helping! This page covers the development setup and the conventions that keep the tool reliable.

## Development setup

```bash
git clone https://github.com/AI-and-Materials-Simulations-Laboratory/buddy-builder
cd buddy-builder
./setup.sh --tests          # .venv, dependencies, self-check, full test suite
./start.sh                  # run the app at http://127.0.0.1:8866
```

After a code change, restart `./start.sh` and reload the page; a running app keeps the old code.

## Making a change

1. Create a branch from `main`.
2. Keep the scientific code (`src/nanobuilder/*.py`) free of GUI imports; `gui.py` and `assistant/panel.py` are the only widget modules.
3. Add or update tests in `tests/`. Tests must not touch the network or need an API key: use `FixtureProvider` for model replies (see `tests/test_assistant.py`).
4. Run `.venv/bin/python -m pytest -q`; everything must pass.
5. Update the docs when behaviour changes: [docs/USER_GUIDE.md](docs/USER_GUIDE.md) for options, [docs/ai_assistant.md](docs/ai_assistant.md) for Ask Buddy.
6. Open a pull request describing what changed and how you checked it.

## Conventions

- **Honest outputs.** Exports are structure/topology only and say so; never imply a force field or an optimized geometry.
- **AI writes commands, not code.** Ask Buddy may only emit commands from the console vocabulary (`assistant/agents.py: ALLOWED`); every plan goes through the verifier before it changes the model.
- **Nothing the user did not ask for.** New structures get no ligands unless requested, and every default Buddy fills in is listed for the user.
- **Secrets stay local.** Keys live in `.env` (git-ignored) or in the session only; never log, print or commit them.
- **Reproducible builds.** Builds are seeded; the same `config.json` must give the same atoms.

## Reporting a problem

Open an issue with what you did, what you expected and what happened. Attach the `session.py` from **Save session** if you can: it rebuilds your model exactly.
