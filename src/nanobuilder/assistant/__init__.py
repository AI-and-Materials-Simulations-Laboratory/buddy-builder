"""Ask Buddy: plain English to checked Buddy console commands (chat, console AI, explanations).

The model only writes commands from the console vocabulary; a deterministic
verifier runs them on a hidden copy of the app before anything changes. With
no API key, an offline parser handles the common requests. See
docs/ai_assistant.md.
"""

from .agents import Assistant, Check, Plan, Verifier
from .provider import (AnthropicProvider, DisabledProvider, FixtureProvider, OpenAICompatibleProvider,
                       from_env)

__all__ = ["Assistant", "Plan", "Check", "Verifier", "from_env", "FixtureProvider", "DisabledProvider",
           "OpenAICompatibleProvider", "AnthropicProvider"]
