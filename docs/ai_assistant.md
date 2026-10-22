# Ask Buddy: the AI assistant

Ask Buddy lets you build and edit models in plain English. It is in three places:

1. **Chat:** the *Ask Buddy* panel above the tabs. Type a request such as "a 4 nm gold sphere with 60 dodecanethiol, LAMMPS and PDB". Buddy plans it, checks it and, with *Auto-apply & build* ticked, builds it in the viewer.

   Follow-ups change the model on the right:

   | Request | What Buddy does |
   |---|---|
   | "make it 5 nm" | changes the size and builds again |
   | "10 more ligands" | changes the coverage and builds again |
   | "remove ligand 3" | edits the current model |
   | "undo" | undoes the last edit |
   | "highlight the sulfur atoms" | highlights them in the viewer |

   Questions get answers, for example "why did my build fail?" or **Explain build**. Untick Auto-apply to review each plan, then press **Apply plan**.
2. **Python console:**
   - Type plain English and press **AI**. The box fills with console commands that were already tested; press **Run** to run them.
   - **Test** runs whatever is in the box on a hidden copy of Buddy and reports what would happen, without changing anything.
   - `ask("...")` plans a request from the console.
   - `reset()` restores the launch settings.
3. **Python:** `from nanobuilder.assistant import Assistant; Assistant().ask("...", app)`.

## How it works (and why it is safe)

The model never writes coordinates or Python code. It writes only **Buddy console commands**, the same ones you type, such as `material("Au")`, `shape("sphere", diameter=4)` or `ligand("dodecanethiol")`. They are parsed with `ast`, never `eval`'d, and only a fixed list of commands is allowed.

The steps below are listed under "How Buddy worked this out" in each reply:

1. **Requirements agent:** decides whether this is a build or a question, and extracts the stated sizes, counts and names.
2. **Planner agent:** writes the commands. Retrieval (BM25 over the console help, the material catalog, the ligand library, the templates, the examples and `docs/`) gives the planner the names Buddy really has.
3. **Verifier (no AI):** runs the commands on a hidden copy of the app. It then calls `validate_inputs`, does a test build, and checks each requirement against the settings the commands produced, for example "wanted 60 ligands, got 60".
4. **Repair agent:** if the verifier finds problems, the plan is revised, at most twice.
5. **Explainer agent:** answers questions from the app state and the docs.

Applied commands go through the console, so the console log, the tabs, Save session and the certificate show exactly what the AI did.

## Setup

Put a key in `.env` in the tool folder (this file is git-ignored), or in the environment:

```
BUDDY_LLM_PROVIDER=openai_compatible   # or anthropic, or off
BUDDY_LLM_API_KEY=sk-...
BUDDY_LLM_MODEL=gpt-6-luna             # default; any OpenAI-compatible chat model works
# BUDDY_LLM_BASE_URL=https://api.openai.com/v1
```

Buddy also uses `OPENAI_API_KEY` or `ANTHROPIC_API_KEY` (default model `claude-opus-5-5`), and AtomForge's `ATOMFORGE_LLM_*` names. No new packages are needed, because HTTP uses the standard library.

**With no key**, an offline parser handles common requests (material, shape, size, ligand, count, density, order, box, formats, and the follow-up edits above). It also fixes the two most common failures itself: a box too small for the model, and more ligands than fit on the surface. The badge then says "Offline parser". If the model cannot be reached, Buddy falls back to the parser and says so.

## Tests

`tests/test_assistant.py` has 25 tests. They use a fixture model and never touch the network: `conftest.py` blocks `.env` and removes any keys.
