"""Ask Buddy: offline parser, verifier, agents (fixture model, no network), chat and console AI."""

import io
import json
import os
import urllib.error

import pytest

pytest.importorskip("ipywidgets")

from nanobuilder.assistant import Assistant, FixtureProvider  # noqa: E402
from nanobuilder.assistant import provider as prov  # noqa: E402
from nanobuilder.assistant.agents import Plan, Verifier, parse_json, safety_problems  # noqa: E402
from nanobuilder.assistant.knowledge import KnowledgeBase  # noqa: E402
from nanobuilder.assistant.rules import parse  # noqa: E402
from nanobuilder.gui import BuilderApp  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return BuilderApp()


def offline():
    return Assistant(provider=prov.DisabledProvider("test"))


# ------------------------------------------------------------------ offline parser

@pytest.mark.parametrize("text, expect", [
    ("a 4 nm gold sphere with 60 dodecanethiol in fibonacci order, export LAMMPS and PDB",
     ['reset()', 'material("Au")', 'shape("sphere", diameter=4)', 'ligand("dodecanethiol")',
      'coverage(count=60, order="fibonacci")', 'export("lammps", "pdb")']),
    ("silver nanorod 2 nm diameter 8 nm long with 20 MUA",
     ['reset()', 'material("Ag")', 'shape("rod", diameter=2, length=8)', 'ligand("mua")',
      'coverage(count=20)']),
    ("TiO2 nanoparticle in a 10 nm cubic box",
     ['reset()', 'material("TiO2")', 'shape("sphere")', 'box("cube", edge=10)']),
    ("bare CdSe quantum dot 3 nm",
     ['reset()', 'material("CdSe")', 'shape("sphere", diameter=3)', 'coverage(count=0)']),
    ("highlight the sulfur atoms", ['select(element="S")']),
    ("remove ligand 3 and 5", ['remove_ligand(3, 5)']),
    ("undo", ['undo()']),
])
def test_offline_parser(text, expect):
    assert parse(text)["commands"] == expect


def test_offline_parser_relative_changes_use_current_settings(app):
    app.reset()
    cfg = app.config()
    assert parse("10 more ligands", cfg=cfg)["commands"] == [f"coverage(count={cfg.coverage.target + 10})"]
    assert parse("make it bigger", cfg=cfg)["commands"] == [
        f'shape("sphere", diameter={cfg.geometry.diameter_nm + 1:g})']
    assert parse("why did my build fail?")["question"]


# ------------------------------------------------------------------ safety + verifier

def test_plans_are_limited_to_the_console_vocabulary():
    bad = ['__import__("os").system("ls")', 'open("x")', 'material(x)', "build(); x = 1",
           'settings({"a": 1})']
    problems = safety_problems(bad)
    assert len(problems) == len(bad)
    assert safety_problems(['material("Au")', 'coverage(count=5)', 'select(element="S")']) == []


def test_verifier_runs_on_a_hidden_copy_and_checks_requirements(app):
    app.reset()
    before = app.config().to_dict()
    plan = Plan("test", commands=['reset()', 'material("Ag")', 'shape("sphere", diameter=2.5)',
                                  'ligand("mua")', 'coverage(count=8)'],
                requirements={"material": "silver", "shape": "sphere", "diameter_nm": 2.5,
                              "ligand": "mua", "ligand_count": 9})
    Verifier().verify(plan, app)
    assert app.config().to_dict() == before  # the live app is untouched
    checks = {c.label: c.ok for c in plan.checks}
    assert checks["material"] and checks["diameter"] and checks["ligand"]
    assert not checks["ligand count"]  # 8 planned, 9 asked for
    assert plan.dry_run["ligands"] == 8 and plan.dry_run["status"] in ("ok", "warning")
    assert any("ligand count" in p for p in plan.problems)


def test_verifier_reports_console_errors():
    plan = Plan("x", commands=['material("unobtainium")'])
    Verifier(dry_run=False).verify(plan)
    assert plan.problems and "unobtainium" in plan.problems[0]


def test_offline_repair_enlarges_a_small_box(app):
    app.reset()
    plan = offline().ask("a 5 nm gold sphere with 30 dodecanethiol", app)
    assert plan.ok, plan.problems
    assert any(c.startswith('box("cube"') for c in plan.commands)
    assert any("Box enlarged" in a for a in plan.assumptions)


# ------------------------------------------------------------------ agents with a fixture model

def test_ai_pipeline_requirements_plan_verify_repair(app):
    app.reset()
    replies = [
        {"intent": "build", "start_fresh": True, "clarify": None,
         "requirements": {"material": "Au", "shape": "sphere", "diameter_nm": 3, "ligand": "mpa",
                          "ligand_count": 20}},
        # first plan uses an English name the console cannot find -> verifier fails it
        "```json\n" + json.dumps({"commands": ['reset()', 'material("unobtainium")'],
                                  "assumptions": [], "explanation": "draft"}) + "\n```",
        {"commands": ['reset()', 'material("Au")', 'shape("sphere", diameter=3)', 'ligand("mpa")',
                      'coverage(count=20)'], "assumptions": ["fibonacci order kept"],
         "explanation": "Au sphere with MPA"},
    ]
    fx = FixtureProvider(replies)
    plan = Assistant(provider=fx, mode="thorough").ask("3 nm gold particle with 20 MPA", app)
    assert plan.mode == "ai" and plan.ok, plan.problems
    assert plan.revisions == 1
    assert [t["agent"] for t in plan.trace] == ["requirements", "planner", "verifier", "repair 1",
                                                "verifier"]
    assert "ligand(\"mpa\")" in plan.commands
    assert "unobtainium" in fx.calls[2]["user"]  # the repair agent saw the verifier's problem
    assert "dodecanethiol" in fx.calls[1]["user"] or "mpa" in fx.calls[1]["user"].lower()  # retrieved refs


def test_ai_follow_up_sees_the_conversation(app):
    fx = FixtureProvider([
        {"intent": "build", "start_fresh": False, "requirements": {"ligand_count": 30}},
        {"commands": ["coverage(count=30)"], "assumptions": [], "explanation": "more"},
    ])
    history = [{"request": "a gold sphere", "commands": ['material("Au")']}]
    plan = Assistant(provider=fx, mode="thorough").ask("more ligands, 30", app, history=history)
    assert plan.ok
    assert "a gold sphere" in fx.calls[0]["user"] and "a gold sphere" in fx.calls[1]["user"]


def test_ai_invalid_json_gets_one_retry_then_question_is_answered(app):
    fx = FixtureProvider(["not json", {"intent": "question", "requirements": {}},
                          {"answer": "Thiolates bind Au through S.", "commands": [], "sources": []}])
    plan = Assistant(provider=fx, mode="thorough").ask("how do thiols bind?", app)
    assert plan.intent == "question" and "Au" in plan.answer
    assert "not valid JSON" in fx.calls[1]["user"]


def test_model_failure_falls_back_to_the_offline_parser(app):
    fx = FixtureProvider([])  # every call raises ProviderError
    plan = Assistant(provider=fx, mode="thorough").ask("a 3 nm silver sphere with 10 mua", app)
    assert plan.mode == "offline" and plan.ok
    assert any("AI model unavailable" in n for n in plan.notes)


def test_parse_json_tolerates_fences_and_prose():
    assert parse_json('Sure:\n```json\n{"a": 1}\n```') == {"a": 1}
    assert parse_json('text {"b": [1, 2]} more') == {"b": [1, 2]}
    with pytest.raises(ValueError):
        parse_json("no object")


# ------------------------------------------------------------------ providers

def test_provider_from_env_without_key_is_disabled(monkeypatch):
    for k in ("BUDDY_LLM_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY", "ATOMFORGE_LLM_API_KEY",
              "BUDDY_LLM_PROVIDER", "ATOMFORGE_LLM_PROVIDER"):
        monkeypatch.delenv(k, raising=False)
    p = prov.from_env()
    assert not p.available and "offline" in p.label()
    monkeypatch.setenv("OPENAI_API_KEY", "sk-test")
    p = prov.from_env()
    assert isinstance(p, prov.OpenAICompatibleProvider) and p.model == prov.DEFAULT_OPENAI_MODEL
    monkeypatch.setenv("BUDDY_LLM_PROVIDER", "off")
    assert not prov.from_env().available


def test_dotenv_does_not_override_and_strips_quotes(tmp_path, monkeypatch):
    monkeypatch.delenv("BUDDY_DOTENV", raising=False)
    env = tmp_path / ".env"
    env.write_text("# comment\nexport BB_T1='one'\nBB_T2=\"two\"\nBB_T3=keep\n")
    monkeypatch.setenv("BB_T3", "already")
    prov.load_dotenv([env])
    assert os.environ["BB_T1"] == "one" and os.environ["BB_T2"] == "two"
    assert os.environ["BB_T3"] == "already"
    for k in ("BB_T1", "BB_T2"):
        monkeypatch.delenv(k)


def test_openai_provider_falls_back_to_max_tokens(monkeypatch):
    sent = []

    class Resp:
        def __init__(self, body):
            self.body = body

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self):
            return json.dumps(self.body).encode()

    def fake_urlopen(req, timeout=0):
        body = json.loads(req.data)
        sent.append(body)
        assert req.get_header("Authorization") == "Bearer sk-x"
        if "max_completion_tokens" in body:
            raise urllib.error.HTTPError(req.full_url, 400, "bad", {},
                                         io.BytesIO(b"unsupported: max_completion_tokens"))
        return Resp({"model": "m", "choices": [{"message": {"content": "{\"ok\": 1}"}}]})

    monkeypatch.setattr(prov.urllib.request, "urlopen", fake_urlopen)
    out = prov.OpenAICompatibleProvider("sk-x", "m").complete("s", "u")
    assert out.text == '{"ok": 1}' and "max_tokens" in sent[-1] and len(sent) == 2


# ------------------------------------------------------------------ knowledge

def test_knowledge_finds_ligands_materials_and_commands():
    kb = KnowledgeBase()
    assert any("dodecanethiol" in c.text for c in kb.search("dodecanethiol gold", k=5))
    assert any(c.source == "console" and "coverage" in c.title for c in kb.search("coverage fibonacci order", k=5))
    assert kb.search("graphene oxide", k=3)


# ------------------------------------------------------------------ GUI: chat and console AI

def test_chat_builds_and_edits_the_live_model():
    a = BuilderApp()
    a.assistant.assistant = offline()
    chat = a.assistant
    plan = chat.send("a 3 nm gold sphere with 20 hexanethiol")
    assert plan.ok and a.result is not None
    assert a.result.assembly.achieved_ligands == 20
    chat.send("5 more ligands")
    assert a.result.assembly.achieved_ligands == 25
    chat.send("remove ligand 2")
    assert a.result.assembly.achieved_ligands == 24
    chat.send("undo")
    assert a.result.assembly.achieved_ligands == 25
    chat.send("highlight the sulfur atoms")
    assert len(a.console.selection) == 25
    assert "bb-turn-user" in chat.transcript.value and "Built:" in chat.transcript.value
    assert any("# Ask Buddy:" in t for _, t in a.console.lines)  # chat edits echo in the console


def test_chat_review_mode_does_not_change_the_model_until_apply():
    a = BuilderApp()
    a.assistant.assistant = offline()
    a.assistant.auto.value = False
    a.assistant.send("a 3 nm silver sphere with 10 mua")
    assert a.result is None and not a.assistant.apply_btn.disabled
    a.assistant.apply_btn.click()
    assert a.result is not None and a.result.assembly.achieved_ligands == 10


def test_console_ai_writes_tested_code_and_ask_command():
    a = BuilderApp()
    a.assistant.assistant = offline()
    a.console_in.value = "a 3 nm silver sphere with 10 mua"
    a.assistant.console_ai.click()
    assert a.console_in.value.splitlines()[1] == 'material("Ag")'
    assert any("passed" in t for _, t in a.console.lines)
    a.assistant.console_test.click()
    assert any("hidden copy" in t for _, t in a.console.lines)
    assert a.result is None  # nothing ran yet
    a.console_run.click()
    assert a.config().material.preset_id == "Ag_fcc"
    a.console.run('ask("a 3 nm gold sphere with 12 mpa")')
    assert any("Plan ready" in t for _, t in a.console.lines)
    a.console.run("make me a gold ball please")
    assert any("Plain English?" in t for _, t in a.console.lines)


def test_reset_command_and_existing_layout_kept():
    a = BuilderApp()
    a.console.run('material("Ag")')
    a.console.run("reset()")
    assert a.config().material.preset_id == "Au_fcc"
    assert a.tabs.get_title(0) == "1 Material"  # the chat sits above the same six tabs
    assert BuilderApp(assistant=False).assistant is None


# ------------------------------------------------------------------ only what was asked for

def test_new_structure_without_ligands_gets_none_and_says_so(app):
    plan = offline().ask("a 3nm graphene oxide sheet", app)
    assert plan.ok, plan.problems
    assert plan.commands[-1] == "coverage(count=0)" and 'carbon(size=(3, 3))' in plan.commands
    assert plan.dry_run["ligands"] == 0
    assert any(a.startswith("No ligands") for a in plan.assumptions)
    assert {c.label for c in plan.checks} >= {"sheet size", "ligand count"}


def test_ligands_are_kept_when_asked_and_default_count_is_named(app):
    plan = offline().ask("a 3 nm gold sphere coated with dodecanethiol", app)
    assert 'ligand("dodecanethiol")' in plan.commands and "coverage(count=0)" not in plan.commands
    assert any("Ligand count not given" in a for a in plan.assumptions)


def test_model_cannot_add_unrequested_ligands(app):
    fx = FixtureProvider([
        {"intent": "build", "start_fresh": True, "requirements": {"material": "graphene oxide",
                                                                   "carbon_size_nm": [3, 3]}},
        {"commands": ['reset()', 'material("graphene oxide")', 'carbon(size=(3, 3))',
                      'ligand("phenyl")', 'coverage(count=12)'], "assumptions": [], "explanation": ""},
    ])
    plan = Assistant(provider=fx, mode="thorough").ask("a 3 nm graphene oxide sheet", app)
    assert plan.ok and plan.commands[-1] == "coverage(count=0)"
    assert not any(c.startswith("ligand(") for c in plan.commands)


def test_no_material_means_buddy_asks(app):
    plan = offline().ask("make a nanoparticle", app)
    assert plan.clarify.startswith("Which material") and not plan.commands
    fx = FixtureProvider([
        {"intent": "build", "start_fresh": True, "requirements": {"shape": "sphere"}},
        {"commands": ['reset()', 'shape("sphere", diameter=3)'], "assumptions": [], "explanation": ""},
    ])
    plan = Assistant(provider=fx, mode="thorough").ask("build me a 3 nm sphere", app)
    assert plan.clarify and not plan.commands  # the model did not name a material: ask, don't guess


# ------------------------------------------------------------------ speed

def test_fast_mode_plans_in_one_model_call(app):
    fx = FixtureProvider([{"intent": "build", "start_fresh": True, "clarify": None,
                           "requirements": {"material": "Ag", "diameter_nm": 3, "ligand": "mua",
                                            "ligand_count": 10},
                           "commands": ['reset()', 'material("Ag")', 'shape("sphere", diameter=3)',
                                        'ligand("mua")', 'coverage(count=10)'],
                           "assumptions": [], "explanation": "Ag sphere"}])
    plan = Assistant(provider=fx).ask("a 3 nm silver sphere with 10 MUA", app)
    assert plan.ok and len(fx.calls) == 1 and plan.result is not None


def test_simple_edits_need_no_model_call(app):
    fx = FixtureProvider([])  # any model call would fail
    a = Assistant(provider=fx)
    for text, cmd in (("highlight the sulfur atoms", 'select(element="S")'), ("undo", "undo()"),
                      ("10 more ligands", None)):
        plan = a.ask(text, app)
        assert plan.mode == "instant" and not fx.calls
        assert cmd is None or plan.commands == [cmd]
    assert a.ask("make it bigger and add 10 more ligands", app).mode != "instant"  # two asks: AI


def test_apply_reuses_the_checked_build():
    a = BuilderApp()
    a.assistant.assistant = offline()
    import nanobuilder.gui as gui
    calls = []
    real = gui.run_pipeline
    gui.run_pipeline = lambda *x, **k: calls.append(1) or real(*x, **k)
    try:
        a.assistant.send("a 3 nm gold sphere with 10 mpa")
    finally:
        gui.run_pipeline = real
    assert a.result is not None and a.result.assembly.achieved_ligands == 10
    assert calls == []  # the live app did not build a second time


# ------------------------------------------------------------------ AI key box

def test_key_box_connects_without_keeping_or_showing_the_key(monkeypatch, tmp_path):
    a = BuilderApp()
    a.assistant.assistant = offline()
    chat = a.assistant
    seen = []

    def fake_complete(self, system, user, max_tokens=1500):
        seen.append(self._key)
        return prov.ProviderResponse("{}", self.model)

    monkeypatch.setattr(prov.OpenAICompatibleProvider, "complete", fake_complete)
    monkeypatch.setattr(prov, "TOOL_DIR", tmp_path)
    secret = "sk-test-abcdefghijklmnop1234567890"
    chat.toggle_key_box()
    assert chat.key_box.layout.display == ""
    chat.key_input.value = secret
    chat.key_save.value = True
    assert chat.connect()
    assert seen == [secret] and chat.assistant.live
    assert chat.key_input.value == ""  # not kept in the widget
    shown = chat.badge.value + chat.key_status.value + chat.transcript.value
    assert secret not in shown and "gpt-6-luna" in chat.badge.value
    env = (tmp_path / ".env").read_text()
    assert f"BUDDY_LLM_API_KEY={secret}" in env and (tmp_path / ".env").stat().st_mode & 0o077 == 0
    chat.forget_key()
    assert not chat.assistant.live


def test_rejected_key_is_reported_and_not_used(monkeypatch):
    a = BuilderApp()
    a.assistant.assistant = offline()

    def bad(self, *args, **kw):
        raise prov.ProviderError("Model API error 401: Incorrect API key provided: sk-test-abcdefghijk")

    monkeypatch.setattr(prov.OpenAICompatibleProvider, "complete", bad)
    a.assistant.key_input.value = "sk-test-abcdefghijk"
    assert not a.assistant.connect()
    assert not a.assistant.assistant.live
    assert "rejected" in a.assistant.key_status.value
    assert "abcdefghijk" not in a.assistant.key_status.value


def test_save_to_dotenv_replaces_existing_lines(tmp_path):
    env = tmp_path / ".env"
    env.write_text("# mine\nBUDDY_LLM_API_KEY=old\nOTHER=1\n")
    prov.save_to_dotenv({"BUDDY_LLM_API_KEY": "new", "BUDDY_LLM_MODEL": "m"}, env)
    text = env.read_text()
    assert "BUDDY_LLM_API_KEY=new" in text and "old" not in text and "OTHER=1" in text
    assert "BUDDY_LLM_MODEL=m" in text and text.startswith("# mine")
    assert prov.scrub("key sk-proj-ABCDEFGHIJKLMNOP done") == "key sk-proj... done"
