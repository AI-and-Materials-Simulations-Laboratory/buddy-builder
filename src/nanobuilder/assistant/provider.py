"""Language-model adapters for Ask Buddy, using only the standard library.

No new dependencies: HTTP goes through :mod:`urllib`. Keys come from the
environment or a ``.env`` file in the tool folder (git-ignored) and are never
logged or shown. With no key the assistant uses the offline parser instead and
says so.

Environment (first match wins):

- ``BUDDY_LLM_PROVIDER``: ``openai_compatible``, ``anthropic`` or ``off``
- ``BUDDY_LLM_MODEL``, ``BUDDY_LLM_API_KEY``, ``BUDDY_LLM_BASE_URL``
- fallbacks: ``OPENAI_API_KEY`` / ``ANTHROPIC_API_KEY`` and the AtomForge
  names ``ATOMFORGE_LLM_*``.
"""

from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parents[3]  # buddy-builder/ (src/nanobuilder/assistant/..)
DEFAULT_OPENAI_MODEL = "gpt-6-luna"
DEFAULT_ANTHROPIC_MODEL = "claude-opus-5-5"
TIMEOUT_S = 60


class ProviderError(RuntimeError):
    """The model could not be reached or answered with an error."""


@dataclass
class ProviderResponse:
    text: str
    model: str
    usage: dict = field(default_factory=dict)
    seconds: float = 0.0


def load_dotenv(paths=None) -> list:
    """Read KEY=VALUE lines into os.environ without overriding what is set.

    Skipped when BUDDY_DOTENV=0 (the tests set this). Returns the files read.
    """
    if os.environ.get("BUDDY_DOTENV") == "0":
        return []
    read = []
    for path in paths or [Path.cwd() / ".env", TOOL_DIR / ".env"]:
        path = Path(path)
        if not path.is_file() or path in read:
            continue
        for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            if line.startswith("export "):
                line = line[7:].strip()
            key, value = line.split("=", 1)
            key, value = key.strip(), value.strip()
            if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
                value = value[1:-1]
            if key and key not in os.environ:
                os.environ[key] = value
        read.append(path)
    return read


def _post(url: str, headers: dict, body: dict) -> dict:
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", **headers},
                                 method="POST")
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:  # noqa: S310 (https URL from config)
        return json.loads(resp.read().decode("utf-8"))


class Provider:
    name = "base"
    model = ""
    available = False
    reason = ""

    def complete(self, system: str, user: str, max_tokens: int = 1500) -> ProviderResponse:
        raise NotImplementedError

    def label(self) -> str:
        return f"{self.model} ({self.name})" if self.available else "offline parser (no AI model)"


class DisabledProvider(Provider):
    name = "off"

    def __init__(self, reason: str = "No AI model configured"):
        self.reason = reason

    def complete(self, system, user, max_tokens=1500):
        raise ProviderError(self.reason)


class OpenAICompatibleProvider(Provider):
    """POST {base}/chat/completions; sends max_completion_tokens, falling back to max_tokens."""

    name = "openai_compatible"
    available = True

    def __init__(self, api_key: str, model: str = DEFAULT_OPENAI_MODEL,
                 base_url: str = "https://api.openai.com/v1"):
        self._key, self.model, self.base_url = api_key, model, base_url.rstrip("/")
        self._token_param = "max_completion_tokens"
        # Lower reasoning effort answers about twice as fast; dropped if the model rejects it.
        self.reasoning_effort = os.environ.get("BUDDY_LLM_REASONING", "low").strip() or None

    def complete(self, system, user, max_tokens=1500):
        t0 = time.time()
        for attempt in range(4):
            body = {"model": self.model, self._token_param: max_tokens,
                    "messages": [{"role": "system", "content": system},
                                 {"role": "user", "content": user}]}
            if self.reasoning_effort and self.reasoning_effort != "off":
                body["reasoning_effort"] = self.reasoning_effort
            try:
                out = _post(f"{self.base_url}/chat/completions",
                            {"Authorization": f"Bearer {self._key}"}, body)
            except urllib.error.HTTPError as exc:
                detail = exc.read().decode("utf-8", errors="replace")[:400]
                if exc.code == 400 and self._token_param in detail and attempt == 0:
                    self._token_param = "max_tokens"
                    continue
                if exc.code == 400 and "reasoning" in detail and self.reasoning_effort:
                    self.reasoning_effort = None  # this model has no reasoning setting
                    continue
                if exc.code in (429, 500, 502, 503, 504) and attempt < 3:
                    time.sleep(1.5 * (attempt + 1))
                    continue
                raise ProviderError(f"Model API error {exc.code}: {detail}") from None
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                raise ProviderError(f"Could not reach the model API: {exc}") from None
            try:
                text = out["choices"][0]["message"]["content"] or ""
            except (KeyError, IndexError, TypeError):
                raise ProviderError("Unexpected reply from the model API") from None
            return ProviderResponse(text, out.get("model", self.model), out.get("usage", {}),
                                    time.time() - t0)
        raise ProviderError("Model API kept failing")


class AnthropicProvider(Provider):
    """Claude Messages API over plain HTTPS."""

    name = "anthropic"
    available = True

    def __init__(self, api_key: str, model: str = DEFAULT_ANTHROPIC_MODEL,
                 base_url: str = "https://api.anthropic.com/v1"):
        self._key, self.model, self.base_url = api_key, model, base_url.rstrip("/")

    def complete(self, system, user, max_tokens=1500):
        t0 = time.time()
        body = {"model": self.model, "max_tokens": max_tokens, "system": system,
                "messages": [{"role": "user", "content": user}]}
        headers = {"x-api-key": self._key, "anthropic-version": "2023-06-01"}
        for attempt in range(4):
            try:
                out = _post(f"{self.base_url}/messages", headers, body)
            except urllib.error.HTTPError as exc:
                if exc.code in (429, 500, 502, 503, 504, 529) and attempt < 3:
                    time.sleep(1.5 * (attempt + 1))
                    continue
                detail = exc.read().decode("utf-8", errors="replace")[:400]
                raise ProviderError(f"Model API error {exc.code}: {detail}") from None
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                raise ProviderError(f"Could not reach the model API: {exc}") from None
            text = "".join(b.get("text", "") for b in out.get("content", []) if b.get("type") == "text")
            return ProviderResponse(text, out.get("model", self.model), out.get("usage", {}),
                                    time.time() - t0)
        raise ProviderError("Model API kept failing")


class FixtureProvider(Provider):
    """Replays canned replies in order (tests and demos; never networked)."""

    name = "fixture"
    available = True

    def __init__(self, replies, model: str = "fixture"):
        self.replies, self.model, self.calls = list(replies), model, []

    def complete(self, system, user, max_tokens=1500):
        self.calls.append({"system": system, "user": user})
        if not self.replies:
            raise ProviderError("Fixture provider has no more replies")
        reply = self.replies.pop(0)
        if callable(reply):
            reply = reply(system, user)
        if not isinstance(reply, str):
            reply = json.dumps(reply)
        return ProviderResponse(reply, self.model)


def _env(*names, default=None):
    for n in names:
        v = os.environ.get(n)
        if v:
            return v
    return default


KEY_PATTERN = __import__("re").compile(r"(sk-[A-Za-z0-9_\-]{4})[A-Za-z0-9_\-*]+")


def scrub(text: str) -> str:
    """Hide anything that looks like an API key in a message shown to the user."""
    return KEY_PATTERN.sub(r"\1...", str(text))


def make_provider(kind: str, api_key: str, model: str = "", base_url: str = "") -> Provider:
    """A provider from values typed in the app (the key stays in memory only)."""
    api_key = (api_key or "").strip()
    if not api_key:
        return DisabledProvider("No API key given")
    if kind == "anthropic":
        return AnthropicProvider(api_key, model.strip() or DEFAULT_ANTHROPIC_MODEL,
                                 base_url.strip() or "https://api.anthropic.com/v1")
    return OpenAICompatibleProvider(api_key, model.strip() or DEFAULT_OPENAI_MODEL,
                                    base_url.strip() or "https://api.openai.com/v1")


def check_provider(provider: Provider) -> str:
    """'' when the key and model work (one tiny call), else the reason, with keys hidden."""
    try:
        provider.complete("Reply with the JSON object {} and nothing else.", "ping", max_tokens=64)
        return ""
    except ProviderError as exc:
        msg = scrub(exc)
        if " 401" in msg or "invalid_api_key" in msg or "authentication" in msg.lower():
            return "The key was rejected (check it was copied whole, and that it is active)."
        if " 404" in msg or "model" in msg.lower() and "not" in msg.lower():
            return f"The model name was not accepted: {msg[:200]}"
        return msg[:300]


def save_to_dotenv(values: dict, path=None) -> Path:
    """Write settings into .env (other lines are kept); the file is made private (0600)."""
    path = Path(path or TOOL_DIR / ".env")
    lines = path.read_text(encoding="utf-8").splitlines() if path.is_file() else []
    done = set()
    out = []
    for line in lines:
        key = line.split("=", 1)[0].strip().removeprefix("export ").strip() if "=" in line else ""
        if key in values:
            if key not in done:
                out.append(f"{key}={values[key]}")
                done.add(key)
            continue
        out.append(line)
    out += [f"{k}={v}" for k, v in values.items() if k not in done]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    try:
        path.chmod(0o600)
    except OSError:
        pass
    return path


def from_env() -> Provider:
    """The provider the environment (or .env) asks for; Disabled when there is no key."""
    load_dotenv()
    choice = (_env("BUDDY_LLM_PROVIDER", "ATOMFORGE_LLM_PROVIDER") or "").lower()
    if choice in ("off", "none", "disabled", "0"):
        return DisabledProvider("AI turned off (BUDDY_LLM_PROVIDER=off)")
    if os.environ.get("ATOMFORGE_LLM_ENABLED") == "0" and not _env("BUDDY_LLM_PROVIDER"):
        return DisabledProvider("AI turned off")
    model = _env("BUDDY_LLM_MODEL")
    base = _env("BUDDY_LLM_BASE_URL")
    if choice == "anthropic" or (not choice and _env("ANTHROPIC_API_KEY")
                                 and not _env("BUDDY_LLM_API_KEY", "OPENAI_API_KEY",
                                              "ATOMFORGE_LLM_API_KEY")):
        key = _env("BUDDY_LLM_API_KEY", "ANTHROPIC_API_KEY")
        if not key:
            return DisabledProvider("No Anthropic API key (set BUDDY_LLM_API_KEY or ANTHROPIC_API_KEY)")
        return AnthropicProvider(key, model or DEFAULT_ANTHROPIC_MODEL,
                                 base or "https://api.anthropic.com/v1")
    key = _env("BUDDY_LLM_API_KEY", "OPENAI_API_KEY", "ATOMFORGE_LLM_API_KEY")
    if not key:
        return DisabledProvider("No API key: set BUDDY_LLM_API_KEY (or OPENAI_API_KEY) in .env to "
                                "use an AI model; the offline parser is used meanwhile")
    return OpenAICompatibleProvider(key, model or _env("ATOMFORGE_LLM_MODEL") or DEFAULT_OPENAI_MODEL,
                                    base or _env("ATOMFORGE_LLM_BASE_URL") or "https://api.openai.com/v1")


__all__ = ["Provider", "ProviderResponse", "ProviderError", "DisabledProvider",
           "OpenAICompatibleProvider", "AnthropicProvider", "FixtureProvider", "from_env",
           "load_dotenv", "make_provider", "check_provider", "save_to_dotenv", "scrub"]
