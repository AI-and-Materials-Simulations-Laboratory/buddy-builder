"""Ask Buddy chat: talk to Buddy, and the model in the viewer changes.

Each message is planned as console commands, checked on a hidden copy of the
app (see :class:`.agents.Verifier`) and, with *Auto-apply* on, run through the
same Python console as typed commands, so the console log, the tabs and the
viewer always show what the chat did. Edits of the current model
(``remove_ligand``, ``select``...) act on the build shown in the viewer.

The console gets two AI buttons as well: **AI** turns plain English in the
console box into console code (checked first, shown for review, not run), and
**Test** runs the code in the box on the hidden copy and reports what would
happen, without changing anything.
"""

from __future__ import annotations

import html
import traceback

import ipywidgets as w

from . import provider as prov
from .agents import ON_MODEL, SETUP, Assistant, Plan, command_name, split_commands

MAX_TURNS = 30
WELCOME = ("Hi, I'm Buddy. Tell me what to build, or how to change the model on the right, "
           "e.g. <i>a 4 nm gold sphere with 60 dodecanethiol, LAMMPS and PDB</i>, then "
           "<i>make it 5 nm</i>, <i>10 more ligands</i>, <i>highlight the sulfur atoms</i>, "
           "or <i>why did my build fail?</i>")
EXAMPLES = [
    ("Try an example...", ""),
    ("4 nm gold sphere, 60 dodecanethiol, Fibonacci", "a 4 nm gold sphere with 60 dodecanethiol "
     "in fibonacci order, export LAMMPS and PDB"),
    ("Silver nanorod with MUA", "silver nanorod 2 nm diameter 8 nm long with 20 MUA"),
    ("Au(111) slab, 2 thiols / nm²", "Au(111) slab 5 x 5 nm, 4 layers, 2 mercaptoethanol per nm2 on top"),
    ("Graphene oxide flake", "graphene oxide 6 x 6 nm with 20% oxygen"),
    ("CdSe quantum dot with MPA", "3 nm CdSe quantum dot with 20 MPA"),
    ("C60 + butylamine", "C60 with 6 butylamine"),
    ("Why did my build fail?", "why did my build fail?"),
]


def _esc(text) -> str:
    return html.escape(str(text))


def _msg_html(kind: str, text: str) -> str:
    return f"<div class='nb-msg nb-msg-{kind}'><ul><li>{_esc(text)}</li></ul></div>"


class AssistantPanel:
    """The chat widget plus the console's AI and Test buttons; bound to a BuilderApp."""

    def __init__(self, app, assistant: Assistant = None):
        self.app = app
        self.assistant = assistant or Assistant()
        self.plan = None
        self.turns = []  # (role, html)
        self.history = []  # [{"request", "commands"}] for follow-ups

        self.transcript = w.HTML()
        self.transcript.add_class("bb-chat-log")
        self.badge = w.HTML()
        self.input = w.Textarea(placeholder="Message Buddy: what to build or change...  (Enter sends, "
                                            "Shift+Enter adds a line)",
                                rows=2, continuous_update=False,
                                layout=w.Layout(width="auto", flex="1 1 auto"))
        self.input.add_class("bb-chat-input")
        self.send_btn = w.Button(description="Send", icon="paper-plane", button_style="primary",
                                 tooltip="Plan, check and (with Auto-apply) build it",
                                 layout=w.Layout(width="auto", flex="0 0 auto"))
        self.auto = w.Checkbox(value=True, description="Auto-apply & build", indent=False,
                               tooltip="Apply checked plans right away and rebuild the model; "
                                       "untick to review each plan first",
                               layout=w.Layout(width="auto"))
        self.apply_btn = w.Button(description="Apply plan", icon="check", disabled=True,
                                  tooltip="Run the last plan's commands in the console and build",
                                  layout=w.Layout(width="auto"))
        self.explain_btn = w.Button(description="Explain build", icon="question-circle",
                                    tooltip="Buddy explains the current model or why the build failed",
                                    layout=w.Layout(width="auto"))
        self.clear_btn = w.Button(description="New chat", icon="eraser",
                                  tooltip="Clear the conversation (the model is kept)",
                                  layout=w.Layout(width="auto"))
        self.examples = w.Dropdown(options=EXAMPLES, value="",
                                   layout=w.Layout(width="auto", max_width="260px"))
        self.key_btn = w.Button(description="AI key", icon="key",
                                tooltip="Paste an OpenAI (or Anthropic) API key to use an AI model",
                                layout=w.Layout(width="auto"))
        self.key_box = self._key_settings()
        send_row = w.HBox([self.input, self.send_btn],
                          layout=w.Layout(align_items="flex-start", margin="4px 0 0"))
        send_row.add_class("bb-chat-row")
        tools = w.HBox([self.auto, self.apply_btn, self.explain_btn, self.clear_btn, self.key_btn,
                        self.examples],
                       layout=w.Layout(flex_flow="row wrap", align_items="center", grid_gap="6px",
                                       margin="4px 0 0"))
        head = w.HBox([w.HTML("<div class='nb-h bb-chat-h'>Ask Buddy</div>"), self.badge],
                      layout=w.Layout(justify_content="space-between", align_items="center",
                                      flex_flow="row wrap"))
        self.box = w.VBox([head, self.key_box, self.transcript, send_row, tools],
                          layout=w.Layout(margin="0 0 10px"))
        self.box.add_class("bb-chat")

        # console AI buttons (added to the console row by the GUI)
        self.console_ai = w.Button(description="AI", icon="magic",
                                   tooltip="Turn plain English in the console box into commands "
                                           "(checked on a hidden copy first; press Run to run them)",
                                   layout=w.Layout(width="auto", flex="0 0 auto"))
        self.console_test = w.Button(description="Test", icon="flask",
                                     tooltip="Run the commands in the box on a hidden copy of Buddy "
                                             "and report problems, without changing anything",
                                     layout=w.Layout(width="auto", flex="0 0 auto"))

        self.send_btn.on_click(self._safe(self._on_send))
        self.input.observe(self._on_enter, names="value")
        self.apply_btn.on_click(self._safe(lambda _: self.apply(build=True)))
        self.explain_btn.on_click(self._safe(lambda _: self.explain()))
        self.clear_btn.on_click(self._safe(lambda _: self.clear()))
        self.examples.observe(self._on_example, names="value")
        self.console_ai.on_click(self._safe(lambda _: self.code_from_text()))
        self.key_btn.on_click(self._safe(lambda _: self.toggle_key_box()))
        self.key_connect.on_click(self._safe(lambda _: self.connect()))
        self.key_forget.on_click(self._safe(lambda _: self.forget_key()))
        self.key_kind.observe(self._on_key_kind, names="value")
        self.console_test.on_click(self._safe(lambda _: self.test_code()))
        self._sending = False
        self._refresh_badge()
        self._say("buddy", WELCOME + ("" if self.assistant.live else
                                      "<br><small>No AI model yet: common requests still work. "
                                      "For free-form requests and questions, click <b>AI key</b> "
                                      "and paste an OpenAI key.</small>"))

    # ------------------------------------------------------------------ AI key
    def _key_settings(self):
        """Provider, model and key fields; the key is kept in memory for this session only,
        unless 'Save' is ticked (then it goes to .env on this computer)."""
        self.key_kind = w.Dropdown(options=[("OpenAI", "openai_compatible"), ("Anthropic (Claude)", "anthropic")],
                                   value="openai_compatible", description="Provider",
                                   style={"description_width": "70px"},
                                   layout=w.Layout(width="auto", max_width="260px"))
        self.key_model = w.Text(value=prov.DEFAULT_OPENAI_MODEL, description="Model",
                                style={"description_width": "50px"},
                                layout=w.Layout(width="auto", max_width="240px"))
        self.key_input = w.Password(placeholder="Paste your API key (sk-...)",
                                    layout=w.Layout(width="auto", flex="1 1 auto"))
        self.key_input.add_class("bb-key-input")
        self.key_connect = w.Button(description="Connect", icon="plug", button_style="primary",
                                    tooltip="Check the key with one tiny request, then use it",
                                    layout=w.Layout(width="auto", flex="0 0 auto"))
        self.key_forget = w.Button(description="Forget", icon="times",
                                   tooltip="Stop using the key in this session",
                                   layout=w.Layout(width="auto", flex="0 0 auto"))
        self.key_save = w.Checkbox(value=False, indent=False,
                                   description="Remember on this computer (saves to .env)",
                                   tooltip="Only on your own computer; never on a shared server",
                                   layout=w.Layout(width="auto"))
        self.key_status = w.HTML("<small>The key stays in this browser session's memory and is never "
                                 "shown, logged or exported. Get an OpenAI key at "
                                 "<a href='https://platform.openai.com/api-keys' target='_blank' "
                                 "rel='noopener'>platform.openai.com/api-keys</a>.</small>")
        box = w.VBox([w.HBox([self.key_kind, self.key_model],
                             layout=w.Layout(flex_flow="row wrap", grid_gap="6px")),
                      w.HBox([self.key_input, self.key_connect, self.key_forget],
                             layout=w.Layout(align_items="center", grid_gap="6px")),
                      self.key_save, self.key_status],
                     layout=w.Layout(display="none", margin="2px 0 8px"))
        box.add_class("bb-key-box")
        return box

    def toggle_key_box(self):
        self.key_box.layout.display = "" if self.key_box.layout.display == "none" else "none"

    def _on_key_kind(self, change):
        defaults = {"openai_compatible": prov.DEFAULT_OPENAI_MODEL, "anthropic": prov.DEFAULT_ANTHROPIC_MODEL}
        if self.key_model.value in defaults.values() or not self.key_model.value.strip():
            self.key_model.value = defaults[change["new"]]
        self.key_input.placeholder = ("Paste your API key (sk-ant-...)" if change["new"] == "anthropic"
                                      else "Paste your API key (sk-...)")

    def connect(self) -> bool:
        """Check the pasted key with one tiny call and, if it works, use it for this session."""
        key = (self.key_input.value or "").strip()
        self.key_input.value = ""  # never keep the key in the widget
        if not key:
            self.key_status.value = _msg_html("warning", "Paste a key first.")
            return False
        self.key_connect.disabled = True
        self.key_status.value = "<small><span class='bb-dots'>Checking the key</span></small>"
        try:
            candidate = prov.make_provider(self.key_kind.value, key, self.key_model.value)
            problem = prov.check_provider(candidate)
        finally:
            self.key_connect.disabled = False
        if problem:
            self.key_status.value = _msg_html("error", f"Not connected: {problem}")
            return False
        self.assistant.provider = candidate
        self._refresh_badge()
        saved = ""
        if self.key_save.value:
            path = prov.save_to_dotenv({"BUDDY_LLM_PROVIDER": self.key_kind.value,
                                        "BUDDY_LLM_API_KEY": key,
                                        "BUDDY_LLM_MODEL": candidate.model})
            saved = f" Saved to {path.name} on this computer (git ignores it)."
        self.key_status.value = _msg_html("ok", f"Connected: {candidate.model}. Ask Buddy now uses this "
                                                f"model.{saved}")
        self.key_box.layout.display = "none"
        self._say("buddy", f"AI connected ({_esc(candidate.model)}). Ask me anything about building your model.")
        return True

    def forget_key(self):
        self.assistant.provider = prov.DisabledProvider("Key removed for this session")
        self.key_input.value = ""
        self._refresh_badge()
        self.key_status.value = _msg_html("info", "Key forgotten for this session; the offline parser is used. "
                                                  "(A key saved in .env is used again next start.)")

    # ------------------------------------------------------------------ display
    def _refresh_badge(self):
        a = self.assistant
        kind = "ai" if a.live else "off"
        text = f"AI: {a.label()}" if a.live else "Offline parser (add an AI key)"
        tip = ("Plans come from the model, then a deterministic check on a hidden copy of Buddy"
               if a.live else getattr(a.provider, "reason", "") or "No AI model configured")
        self.badge.value = f"<span class='bb-badge bb-badge-{kind}' title='{_esc(tip)}'>{_esc(text)}</span>"

    def _say(self, role: str, body_html: str):
        self.turns.append((role, body_html))
        self.turns = self.turns[-MAX_TURNS:]
        names = {"user": "You", "buddy": "Buddy", "thinking": "Buddy"}
        rows = "".join(f"<div class='bb-turn bb-turn-{r}'><div class='bb-who'>{names.get(r, r)}</div>"
                       f"<div class='bb-bubble'>{b}</div></div>" for r, b in reversed(self.turns))
        self.transcript.value = f"<div class='bb-chat-scroll'>{rows}</div>"

    def _safe(self, fn):
        def run(*args):
            try:
                fn(*args)
            except Exception as exc:  # never break the app; show the error in the chat
                self._say("buddy", f"<div class='nb-msg nb-msg-error'><ul><li>Something went wrong: "
                                   f"{_esc(prov.scrub(exc))}</li></ul></div><details><summary>Details"
                                   f"</summary><pre>{_esc(prov.scrub(traceback.format_exc()))}</pre></details>")
                self.send_btn.disabled = False
        return run

    def _on_enter(self, change):
        # The Textarea fires on blur/Enter (continuous_update=False); the keydown hook in the GUI
        # clicks Send. Nothing to do here; kept for symmetry with the console.
        return

    def _on_example(self, change):
        if change["new"]:
            self.input.value = change["new"]
            self.examples.value = ""

    def plan_html(self, plan: Plan, applied: str = "") -> str:
        parts = []
        if plan.answer:
            parts.append(f"<div class='bb-answer'>{_esc(plan.answer)}</div>")
        if plan.clarify:
            parts.append(f"<div class='bb-answer'>{_esc(plan.clarify)}</div>")
        if plan.explanation and not plan.answer:
            parts.append(f"<div>{_esc(plan.explanation)}</div>")
        if plan.commands:
            label = "Suggested commands" if plan.answer else "Commands"
            parts.append(f"<div class='bb-sub'>{label}</div><pre class='bb-cmds'>"
                         + _esc("\n".join(plan.commands)) + "</pre>")
        if plan.checks:
            items = "".join(f"<li class='{'bb-ok' if c.ok else 'bb-bad'}'>{'✓' if c.ok else '✗'} "
                            f"{_esc(c.label)}: {_esc(c.detail)}</li>" for c in plan.checks)
            parts.append(f"<div class='bb-sub'>Checked against your request</div><ul class='bb-checks'>{items}</ul>")
        d = plan.dry_run
        if d.get("atoms"):
            parts.append(f"<div class='bb-sub'>Test build (hidden copy)</div><div>{d['atoms']} atoms, "
                         f"{d['ligands']}/{d['target']} ligands, validation {_esc(d['status'])}, "
                         f"{d['seconds']} s</div>")
        if plan.assumptions:
            parts.append("<div class='bb-sub'>Buddy assumed (tell me to change any)</div><ul>" +
                         "".join(f"<li>{_esc(a)}</li>" for a in plan.assumptions) + "</ul>")
        if plan.problems:
            parts.append("<div class='nb-msg nb-msg-warning'><ul>" +
                         "".join(f"<li>{_esc(p)}</li>" for p in plan.problems) + "</ul></div>")
        for n in plan.notes:
            parts.append(f"<div class='nb-msg nb-msg-info'><ul><li>{_esc(n)}</li></ul></div>")
        if applied:
            parts.append(applied)
        if plan.trace:
            steps = "".join(
                f"<li>{_esc(t['agent'])}: {'ok' if t.get('ok') else 'needs work'} "
                f"({t.get('seconds', 0)} s)</li>" for t in plan.trace)
            mode = (f"AI ({_esc(plan.model)})" if plan.mode == "ai" else _esc(plan.model or "offline parser"))
            rev = f", {plan.revisions} repair round{'s' if plan.revisions != 1 else ''}" if plan.revisions else ""
            parts.append(f"<details><summary>How Buddy worked this out: {mode}{rev}</summary>"
                         f"<ul>{steps}</ul></details>")
        return "".join(parts) or "<i>No plan.</i>"

    # ------------------------------------------------------------------ actions
    def _on_send(self, _=None):
        text = (self.input.value or "").strip()
        if not text or self._sending:
            return
        self.input.value = ""
        self.send(text)

    def send(self, text: str) -> Plan:
        """One chat turn: plan, check, and (Auto-apply) apply and build."""
        self._sending = True
        self.send_btn.disabled = True
        self.send_btn.description = "Thinking..."
        try:
            self._say("user", _esc(text))
            self._thinking("Buddy is thinking: planning and checking on a hidden copy")
            plan = self.ask(text, show=False)
            applied = ""
            if plan.commands and self.auto.value and not plan.answer and (plan.ok or self._only_soft(plan)):
                self._thinking("Building the model")
                applied = self.apply(build=True, quiet=True)
            self._thinking(None)
            self._say("buddy", self.plan_html(plan, applied))
            return plan
        finally:
            self._thinking(None)
            self._sending = False
            self.send_btn.disabled = False
            self.send_btn.description = "Send"

    def _thinking(self, text):
        """Show (text) or remove (None) the 'Buddy is thinking' bubble; widget updates reach the
        browser while the request is still running, so the user sees progress at once."""
        if self.turns and self.turns[-1][0] == "thinking":
            self.turns.pop()
        if text:
            self._say("thinking", f"<span class='bb-dots'>{_esc(text)}</span>")
        elif self.turns:
            self._say(*self.turns.pop())

    def ask(self, text: str, show: bool = True) -> Plan:
        """Plan and check a request (no change to the model); also used by console ask()."""
        plan = self.assistant.ask(text, self.app, history=self.history[-4:])
        self.plan = plan
        self.apply_btn.disabled = not plan.commands
        if (plan.commands or plan.clarify) and not plan.answer:
            self.history.append({"request": text, "commands": list(plan.commands)})
        if show:
            self._say("user", _esc(text))
            self._say("buddy", self.plan_html(plan))
        return plan

    @staticmethod
    def _only_soft(plan: Plan) -> bool:
        """Problems that still leave a buildable model (a requirement only partly met)."""
        return all(p.startswith(("Only ", "Requirement not met")) for p in plan.problems)

    def apply(self, build: bool = True, quiet: bool = False) -> str:
        """Run the plan's commands in the console (echoed there), then build when settings changed."""
        plan = self.plan
        if plan is None or not plan.commands:
            return ""
        console = self.app.console
        console.log(f"# Ask Buddy: {plan.request}", kind="echo")
        errors, ran = [], []
        for c in plan.commands:
            before = len(console.lines)
            console.run(c)
            ran.append(c)
            errs = [t for k, t in console.lines[before:] if k == "err"]
            if errs:
                errors.append(f"{c}: {errs[0]}")
                break
        changed = any(command_name(c) in SETUP for c in ran)
        out = ""
        if errors:
            out = ("<div class='nb-msg nb-msg-error'><ul>" +
                   "".join(f"<li>{_esc(e)}</li>" for e in errors) + "</ul></div>")
        elif build and changed:
            if plan.result is not None and plan.config is not None:  # reuse the checked test build
                self.app._prebuilt = (plan.config.to_dict(), plan.result)
            before = len(console.lines)
            console.run("build()")
            last = [t for k, t in console.lines[before:] if k == "out"]
            r = self.app.result
            if r is not None and not getattr(self.app, "last_failure", None):
                a = r.assembly
                out = (f"<div class='nb-msg nb-msg-ok'><ul><li>Built: {a.n_atoms} atoms, "
                       f"{a.achieved_ligands}/{a.target_ligands} ligands, validation "
                       f"{_esc(r.report.status)}. The model is in the viewer; downloads are ready."
                       "</li></ul></div>")
            else:
                msg = last[0] if last else "see Status"
                out = (f"<div class='nb-msg nb-msg-error'><ul><li>{_esc(msg)}</li></ul></div>"
                       "<div><small>Ask me <i>why did my build fail?</i> for help.</small></div>")
        elif any(command_name(c) in ON_MODEL for c in ran):
            outs = [t for k, t in console.lines[-6:] if k == "out"]
            out = f"<div class='nb-msg nb-msg-ok'><ul><li>Done. {_esc(outs[-1] if outs else '')}</li></ul></div>"
        else:
            out = "<div class='nb-msg nb-msg-ok'><ul><li>Settings applied.</li></ul></div>"
        self.apply_btn.disabled = True
        if not quiet:
            self._say("buddy", out)
        return out

    def explain(self, question: str = ""):
        self._say("user", _esc(question or "Explain the current build"))
        plan = self.assistant.explain(self.app, question)
        self.plan = plan
        self.apply_btn.disabled = not plan.commands
        self._say("buddy", self.plan_html(plan))
        return plan

    def clear(self):
        self.turns, self.history, self.plan = [], [], None
        self.apply_btn.disabled = True
        self._say("buddy", WELCOME)

    # ------------------------------------------------------------------ console AI
    def code_from_text(self):
        """Console AI button: plain English in the console box -> checked commands in the box."""
        app, console = self.app, self.app.console
        text = (app.console_in.value or "").strip()
        if not text:
            console.log("# AI: type what you want in plain English, then press AI", kind="echo")
            return None
        console.log("# AI is thinking: planning and testing on a hidden copy...", kind="echo")
        self.console_ai.disabled = True
        try:
            plan = self.assistant.ask(text, app, history=self.history[-4:])
        finally:
            self.console_ai.disabled = False
        self.plan = plan
        console.log(f"# AI ({plan.model}): {text}", kind="echo")
        if plan.answer:
            for line in plan.answer.splitlines():
                console.log(line)
        if plan.clarify:
            console.log(f"# {plan.clarify}", kind="echo")
        if plan.commands:
            app.console_in.value = "\n".join(plan.commands)
            self.history.append({"request": text, "commands": list(plan.commands)})
            self._log_test(plan, "# Tested on a hidden copy:")
            console.log("# The commands are in the box: press Run (or Enter) to run them", kind="echo")
        return plan

    def test_code(self):
        """Console Test button: run the box's commands on the hidden copy and report."""
        app, console = self.app, self.app.console
        text = (app.console_in.value or "").strip()
        if not text:
            console.log("# Test: put commands in the box first", kind="echo")
            return None
        plan = Plan(request=text, mode="test", model="test")
        plan.commands = split_commands([text.replace("\n", ";")])
        self.assistant.verifier.verify(plan, app)
        console.log(f"# Test: {len(plan.commands)} command(s) on a hidden copy (nothing changed here)",
                    kind="echo")
        self._log_test(plan, "")
        return plan

    def _log_test(self, plan: Plan, head: str):
        console = self.app.console
        if head:
            console.log(head, kind="echo")
        d = plan.dry_run
        if d.get("atoms"):
            console.log(f"#   test build: {d['atoms']} atoms, {d['ligands']}/{d['target']} ligands, "
                        f"validation {d['status']} ({d['seconds']} s)", kind="echo")
        for c in plan.checks:
            console.log(f"#   {'ok ' if c.ok else 'NO '} {c.label}: {c.detail}", kind="echo")
        for p in plan.problems:
            console.log(f"problem: {p}", kind="err")
        if not plan.problems:
            console.log("#   passed", kind="echo")


__all__ = ["AssistantPanel"]
