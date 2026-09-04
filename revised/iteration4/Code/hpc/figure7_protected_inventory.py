#!/usr/bin/env python3
"""Read-only, full-depth Figure6 and manuscript-text inventory for Figure7 runs."""
import argparse
import hashlib
from pathlib import Path
import re

parser = argparse.ArgumentParser()
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
if root.name != "iteration4":
    raise SystemExit("Expected iteration4")
files = set()
pattern = re.compile(r"(^|/)[^/]*(?:figure6|fig6)[^/]*(/|$)", re.I)
for directory in ("Code/Figures", "Code/hpc", "data/Figures", "Figures", "manuscript/Figures"):
    for path in (root / directory).rglob("*"):
        if path.is_file() and pattern.search(str(path.relative_to(root))):
            files.add(path)
for path in (root / "manuscript").rglob("*"):
    if path.is_file() and path.suffix.lower() in (".tex", ".md", ".docx", ".bib"):
        if "tables" not in path.relative_to(root / "manuscript").parts:
            files.add(path)
if not files:
    raise SystemExit("No protected files found")
for path in sorted(files):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    print(f"{digest.hexdigest()}  {path.relative_to(root)}")
