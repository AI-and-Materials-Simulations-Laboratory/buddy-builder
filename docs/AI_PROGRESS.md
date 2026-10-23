# Ask Buddy (AI features) — progress log

This log lets anyone pick the work up mid-way. Update it after each step.

## Goal
Add AI features to Buddy Builder (natural language → structure) without breaking any existing function.
The baseline before any change was `.venv/bin/python -m pytest -q -p no:cacheprovider`: **244 passed, 1 skipped**.

## Design
The AI only writes Buddy's **console commands**, such as `material("Au")`, `shape("sphere", diameter=4)` and `ligand("dodecanethiol")`. Those are the same safe commands people type. They are parsed with `ast` and never `eval`'d.

New package: `src/nanobuilder/assistant/`.

| File | Purpose |
|---|---|
| `provider.py` | LLM adapters using stdlib `urllib` only (no new dependencies): OpenAI-compatible, Anthropic, Fixture (tests) and Disabled. It also reads `.env`. |
| `knowledge.py` | BM25 retrieval over the console help, the material catalog, the ligand library, the examples, the templates and `docs/*.md`. |
| `rules.py` | Offline parser, so the feature works with no API key (for example on nanoHUB). It turns "4 nm gold sphere with 60 dodecanethiol" into commands. |
| `agents.py` | The pipeline (see below). |
| `panel.py` | The "Ask Buddy" ipywidgets panel (see below). |

The pipeline in `agents.py` runs in this order:

1. **Requirements agent**: works out the intent (build or question) and the structured requirements.
2. **Planner agent**: writes the console commands.
3. **Verifier** (deterministic): runs the commands on a hidden copy of the app, then calls `validate_inputs`, does an optional dry-run build, and checks the result against the requirements.
4. **Repair agent**: fixes failures, at most 2 rounds.
5. **Explainer agent**: answers questions and explains why a build failed.

The panel in `panel.py` has:

- a text box;
- **Plan**, **Apply**, **Apply & build** and **Explain last build** buttons;
- a preview of the commands, a checklist of the requirements, and the agent trace;
- a badge showing whether the plan came from an AI model or the offline parser.

The hooks into `gui.py` are kept to a minimum: one panel at the top of the left column, wrapped so that a failure there cannot break the app. There is also a console command, `ask("...")`.

## Configuration (env or `.env` in the tool folder, git-ignored)

| Variable | Values |
|---|---|
| `BUDDY_LLM_PROVIDER` | `openai_compatible` / `anthropic` / `off` |
| `BUDDY_LLM_MODEL` | model name |
| `BUDDY_LLM_API_KEY` | API key |
| `BUDDY_LLM_BASE_URL` | API base URL |

If these are unset, Buddy falls back to `OPENAI_API_KEY`, `ANTHROPIC_API_KEY` and the AtomForge `ATOMFORGE_LLM_*` names. With no key, it uses the offline parser and says so.

## Status
- [x] Base extracted to `OSI2026/buddy-builder`; the venv is `.venv` (uv, Python 3.12); baseline 244 passed.
- [x] provider.py
- [x] knowledge.py
- [x] rules.py (offline parser; tried on 11 phrasings)
- [x] agents.py (pipeline + verifier + chat history + offline repairs: box too small, too many ligands)
- [x] console.py: `reset()` and `ask()` commands added
- [x] panel.py = **Ask Buddy chat** (user request 3), sitting above the tabs, with the viewer beside it:
  - a transcript of user and Buddy bubbles;
  - each message plans, verifies and, when Auto-apply is ticked, applies the change and rebuilds the live model;
  - follow-ups ("make it bigger") use the chat history;
  - model edits (remove_ligand, select) act on the current build;
  - uses the same nb-* theme classes and fonts (Fredoka/Arial); new CSS goes in theme.py as `.bb-chat*`.
- [x] gui.py hook: `BuilderApp(assistant=True)`; the hidden copy used for checks is created with assistant=False
- [x] Console AI (user request 2): an **AI** button turns plain English typed in the console into console code (shown for review); a **Test** button runs the code on a hidden copy and reports the problems; plain English typed in the console is offered to the AI
- [x] tests/test_assistant.py: 25 tests (fixture provider, no network)
- [x] Full suite: **269 passed, 1 skipped** (244 original + 25 new)
- [x] Live test with gpt-6-luna (key copied into buddy-builder/.env, which is git-ignored): 5 chat turns (build, relative edit, highlight, graphene, question); all ok in about 10 s each after loosening the ligand and carbon-shape checks. Script: scratchpad live_try.py
- [x] README section, docs/ai_assistant.md, `.env` in .gitignore
- [x] Browser check in Jupyter (port 8891, token bbtest): chat → build → highlight, and console AI all work; the style matches (screenshots in scratchpad/smoke/bb-*.png). One page error, `reading 'schema'`, came from the Jupyter frontend; not yet confirmed to be pre-existing.

## Not done / next ideas
- Nothing is committed: buddy-builder is untracked inside OSI2026; ask the user where to push it.
- Possible extras: streaming replies; a voice input; a 'compare two builds' agent.
- [x] theme.py: `.bb-chat*` CSS
- [x] Offline end-to-end check (scratchpad chat_try.py): build, 10 more ligands, highlight S, remove ligand 3, undo, make it 5 nm (auto box fix), graphene oxide, console AI + Test: all work

## Separate app launcher (user request 4: "build a separate app start")
- [x] Add Voila to the local-dev requirements (renders only the widgets, no code cells)
- [x] `start.py` (cross-platform): runs Voila on app.ipynb at http://127.0.0.1:8866 and opens the browser; falls back to Jupyter Notebook when Voila is missing
- [x] `start.sh` (macOS/Linux) and `start.bat` (Windows) wrappers that use .venv
- [x] Browser check: chat, build, viewer and console work under Voila
- [x] README "Run the GUI" section updated; tests still green
- Result: `./start.sh` serves BuddyBuilder.ipynb (code cell only, title "Buddy Builder") at :8866. In the browser the page loads in about 2 s with no code cells; Generate, chat build, chat edit and console all work, with no page errors. The first browser run timed out on the 2nd chat message (cause unknown); the next 2 runs passed fully. Script: scratchpad/smoke/bb_voila.mjs

## Fix (user request 5): no unrequested ligands; confirm or state assumptions
Bug: "3 nm graphene oxide sheet" got 12 ligands, because reset() restores the launch demo (12 SCCO).
- [x] New builds: when the user named no ligand, count or density, drop any ligand/coverage commands, add `coverage(count=0)` and say so ("No ligands: you didn't ask for any")
- [x] Defaults used for unstated essentials (size, layers) are listed as assumptions ("tell me to change any")
- [x] New build with no material named → Buddy asks which material instead of guessing
- [x] Carbon sizes: "3 nm graphene oxide sheet" → carbon(size=(3, 3)), plus a carbon_size_nm requirement check
- [x] Tests + live check
- Result: 273 passed. Live gpt-6-luna:
  - "a 3nm graphene oxide sheet" → 3 x 3 nm sheet, 0 ligands, with a note saying so;
  - "make a nanoparticle" → Buddy asks which material;
  - a gold nanorod with MUA → builds and lists the size and density it chose.
  The model's assumptions and Buddy's default notes are de-duplicated. The app was restarted on :8866.
- Base-builder issue found (not changed): 10 ethanolamine on a 3 nm GO sheet fails validation (atom overlap 0.87 Å); 4 ethanolamine, or a 5 nm sheet, is fine. So the hint for carbon suggests phenyl.

## Fix (user requests 6-7): "not working" after restart + "new chat taking time"
- Cause of "not working": the browser tab was still attached to the old session after I restarted the server (log: "Kernel does not exist", 404). Fix: reload the page. Also there was no sign that Buddy was busy for ~10 s → added a "Buddy is thinking" bubble (shows in 0.08 s) and a "Thinking..." Send button.
- [x] Speed: measure each step; instant path for simple edits (no AI call); a single AI call for planning (fast mode) instead of 2; lower reasoning effort if the API supports it
- Speed result (live, gpt-6-luna):
  - simple edits are instant, 0-1.4 s (no AI call): highlight, remove, undo, N more ligands, make it N nm, add N <ligand>;
  - new builds take about 4 s (one AI call with reasoning_effort=low, instead of 2 calls);
  - questions take about 3 s; New chat is instant;
  - the live app reuses the checked test build instead of building twice.
  BUDDY_AI_MODE=thorough brings back the 2-call mode.
- Also fixed: the console now accepts carbon(kind="graphene oxide") (synonyms added); the checker resolves picker names ("Gold nanorod - [001] axis") and carbon-kind synonyms, so there are no false repairs. Tests: 276 passed.

## Repo (user request 8): create a buddy-builder repo with the whole folder
- [x] git init in buddy-builder; make sure .env/.venv/caches are ignored; scan for secrets
- [x] gh repo create AI-and-Materials-Simulations-Laboratory/buddy-builder --private, then push main (no Claude co-author)
- Pushed: https://github.com/AI-and-Materials-Simulations-Laboratory/buddy-builder (private, main, commit 11bcbfe, 520 files; .env is not in the repo).

## Public + onboarding (user request 9)
- [x] Pre-public audit: secrets in all commits, private links, licence of bundled files
- [x] setup.sh (macOS/Linux) + setup.bat (Windows): Python 3.12 venv (via uv when available), install requirements, create .env from .env.example, self-check
- [x] GETTING_STARTED.md (judges + team): 3 commands, a 5-minute demo script, AI key, troubleshooting; README links to it at the top
- [x] Test the setup from a fresh clone; then make the repo public
- Audit: no secrets in history; HierGO is GPL-3 with its LICENSE; the project LICENSE is still 'NOT YET CHOSEN' (team decision).
- Setup tested from a fresh clone: the uv route takes 27 s; plain python3.12 without uv takes about 1 min; with only system Python 3.9 it stops with a clear message.
- [x] AI key box in the UI (user request 10): provider/model/password + Connect (a test call); session-only by default, optional save to .env; the key is never shown or logged
- [x] Pushed and made the repo public (secret scan of all 5 commits: only fake test keys)
- Key box tested in a browser on a fresh clone with no key: the badge said offline; a wrong key was rejected with a clear message; the real key connected (the badge changed to gpt-6-luna); the key is not in the page HTML; the AI answered. 279 tests passed.

## README (user request 11): neat, professional, open-source style
- [x] New README: centred header, badges, screenshot (docs/images/screenshot.png), features, quick start, usage, outputs (real file names), a mermaid diagram of Ask Buddy, layout, testing, docs index, limitations, contributing, license table, acknowledgements
- [x] Detail moved, not lost: docs/INSTALL.md (manual, conda, Windows, Jupyter, nanoHUB, CLI, troubleshooting), docs/USER_GUIDE.md (all options), CONTRIBUTING.md
- [x] Python version made consistent at 3.12+ (NumPy 2.5 pin): setup.sh, GETTING_STARTED, pyproject requires-python

## License (user request 12)
- [x] MIT (Copyright 2026 AI and Materials Simulations Laboratory); README badge and section, pyproject, deployment notes updated; HierGO-master stays GPL-3.0
