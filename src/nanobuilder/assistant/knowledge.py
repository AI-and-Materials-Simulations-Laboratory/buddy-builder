"""What Buddy knows, as small text chunks with BM25 retrieval (no dependencies).

Sources: the console help, the material catalog, the ligand library, the
attachment templates, the bundled examples and the Markdown files in ``docs/``
and the README. The retrieved chunks go to the model as reference data, so it
uses names Buddy really has instead of inventing them.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .provider import TOOL_DIR

TOKEN = re.compile(r"[a-z0-9]+")


@dataclass(frozen=True)
class Chunk:
    source: str
    title: str
    text: str


def tokens(text: str) -> list:
    return TOKEN.findall(str(text).lower())


def split_markdown(source: str, text: str, max_chars: int = 1200) -> list:
    """One chunk per heading section (long sections are cut into pieces)."""
    out, title, buf = [], source, []

    def flush():
        body = "\n".join(buf).strip()
        for i in range(0, len(body), max_chars):
            if body[i:i + max_chars].strip():
                out.append(Chunk(source, title, body[i:i + max_chars]))

    for line in text.splitlines():
        if line.startswith("#"):
            flush()
            title, buf = line.lstrip("#").strip() or source, []
        else:
            buf.append(line)
    flush()
    return out


def _builtin_chunks() -> list:
    from ..chemistry import TEMPLATES
    from ..console import HELP
    from ..examples import EXAMPLES
    from ..ligand_library import LIBRARY
    from .. import material_catalog as catalog

    out = [Chunk("console", f"command {k}", f"{usage}\n{text}") for k, (usage, text) in HELP.items()]
    for name, key in catalog.NAMES.items():
        e = catalog.CATALOG[key]
        extra = e.get("preset_id") or e.get("mp_id") or ""
        out.append(Chunk("material", name, f"material(\"{name}\")  group: {e['group']}  {extra} "
                                           f"{e.get('formula', '')}"))
    for word, label in catalog.CARBON_WORDS.items():
        out.append(Chunk("material", word, f"material(\"{word}\") opens carbon entry: {label}"))
    for key, e in LIBRARY.items():
        out.append(Chunk("ligand", e["label"],
                         f"ligand(\"{key}\")  {e['label']} ({e['group']}); SMILES {e['smiles']}; "
                         f"template {e['binding_template']}; {e.get('note') or ''}"))
    for key, t in TEMPLATES.items():
        out.append(Chunk("template", t.label, f"template \"{key}\": {t.description}"))
    for key, e in EXAMPLES.items():
        out.append(Chunk("example", e["label"], f"example(\"{key[0]}\")  {e['label']}"))
    return out


@lru_cache(maxsize=1)
def all_chunks() -> tuple:
    chunks = _builtin_chunks()
    files = [TOOL_DIR / "README.md"] + sorted((TOOL_DIR / "docs").glob("*.md"))
    for path in files:
        if path.name == "AI_PROGRESS.md" or not path.is_file():
            continue
        chunks += split_markdown(path.name, path.read_text(encoding="utf-8", errors="replace"))
    return tuple(chunks)


class KnowledgeBase:
    """BM25 (k1=1.5, b=0.75) over :func:`all_chunks`."""

    def __init__(self, chunks=None):
        self.chunks = list(chunks if chunks is not None else all_chunks())
        self.docs = [Counter(tokens(c.title + " " + c.text)) for c in self.chunks]
        self.lengths = [sum(d.values()) for d in self.docs]
        self.avg = (sum(self.lengths) / len(self.lengths)) if self.lengths else 1.0
        df = Counter()
        for d in self.docs:
            df.update(d.keys())
        n = len(self.docs)
        self.idf = {t: math.log(1 + (n - f + 0.5) / (f + 0.5)) for t, f in df.items()}

    def search(self, query: str, k: int = 8, sources=None) -> list:
        q = [t for t in tokens(query) if t in self.idf]
        scored = []
        for i, d in enumerate(self.docs):
            if sources and self.chunks[i].source not in sources:
                continue
            s = 0.0
            for t in q:
                f = d.get(t, 0)
                if f:
                    s += self.idf[t] * f * 2.5 / (f + 1.5 * (0.25 + 0.75 * self.lengths[i] / self.avg))
            if s > 0:
                scored.append((s, i))
        scored.sort(reverse=True)
        return [self.chunks[i] for _, i in scored[:k]]

    def context(self, query: str, k: int = 10) -> str:
        """Retrieved chunks as one text block for a prompt."""
        hits = self.search(query, k=k)
        return "\n".join(f"[{c.source}: {c.title}] {c.text}" for c in hits)


@lru_cache(maxsize=1)
def default_kb() -> KnowledgeBase:
    return KnowledgeBase()


__all__ = ["Chunk", "KnowledgeBase", "default_kb", "split_markdown", "tokens", "all_chunks"]
