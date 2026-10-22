#!/bin/sh
# One-time setup for Buddy Builder (macOS / Linux). Safe to run again.
#   ./setup.sh           install into .venv and check it works
#   ./setup.sh --tests   also run the full test suite (about 40 s)
set -e
cd "$(dirname "$0")"
echo "== Buddy Builder setup =="

if [ -x .venv/bin/python ] && ! .venv/bin/python -c 'import sys' >/dev/null 2>&1; then
  echo "Removing a broken .venv"; rm -rf .venv
fi

if command -v uv >/dev/null 2>&1; then
  echo "Using uv (it downloads Python 3.12 if this computer does not have it)"
  [ -x .venv/bin/python ] || uv venv --python 3.12 .venv
  uv pip install --python .venv/bin/python -r requirements.txt
else
  PY=""
  for c in python3.12 python3.13 python3.11 python3 python; do
    if command -v "$c" >/dev/null 2>&1 && \
       "$c" -c 'import sys; sys.exit(0 if (3, 11) <= sys.version_info[:2] <= (3, 13) else 1)' 2>/dev/null; then
      PY="$c"; break
    fi
  done
  if [ -z "$PY" ]; then
    echo "Python 3.11 to 3.13 is needed (3.12 recommended)."
    echo "  Install Python: https://www.python.org/downloads/"
    echo "  or install uv (it fetches Python for you): curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
  fi
  echo "Using $("$PY" --version)"
  [ -x .venv/bin/python ] || "$PY" -m venv .venv
  .venv/bin/python -m pip install --quiet --upgrade pip
  .venv/bin/python -m pip install -r requirements.txt
fi

# macOS can hide .pth files in synced (iCloud) folders, which breaks imports.
if [ "$(uname)" = "Darwin" ]; then chflags nohidden .venv/lib/python3*/site-packages/*.pth 2>/dev/null || true; fi

if [ ! -f .env ]; then
  cp .env.example .env
  echo "Created .env: add an AI key there to use an AI model in the chat (optional)."
fi

echo "== Self-check =="
.venv/bin/python tools/selfcheck.py

if [ "$1" = "--tests" ]; then
  echo "== Tests =="
  .venv/bin/python -m pytest -q -p no:cacheprovider
fi

echo
echo "Ready. Start Buddy Builder with:  ./start.sh"
echo "New here? Read GETTING_STARTED.md"
