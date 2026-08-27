#!/usr/bin/env python3
"""Verify exact distributions from a hashed requirements lock."""

from __future__ import annotations

import importlib.metadata
import re
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_python_environment.py REQUIREMENTS_LOCK")
    lock_path = Path(sys.argv[1])
    failures = 0
    print("distribution\texpected\tactual\tmatch")
    for raw in lock_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = re.match(r"^([^= ]+)==([^ ]+)", line)
        if not match:
            continue
        name, expected = match.groups()
        try:
            actual = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            actual = "NA"
        ok = actual == expected
        failures += int(not ok)
        print(f"{name}\t{expected}\t{actual}\t{str(ok).upper()}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

