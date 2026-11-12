#!/bin/bash
# Preserve the original helper's stdin confirmation; delegate to the real site tool.
set -euo pipefail
runtime_root="$(cd -- "$(dirname -- "$0")/.." && pwd)"
exec python3 -u "${runtime_root}/monitor.py" extend "$@"
