#!/bin/sh
# Start Buddy Builder as an app in the browser (macOS / Linux). Options: see start.py.
cd "$(dirname "$0")" || exit 1
if [ -x .venv/bin/python ]; then exec .venv/bin/python start.py "$@"; fi
exec python3 start.py "$@"
