#!/usr/bin/env python3
"""Publish figure artifacts to the manuscript-local output directory."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


SCRIPT = Path(__file__).resolve()


def find_workspace_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (
            (candidate / "manager.sh").is_file()
            and (candidate / "Code" / "Figures").is_dir()
        ):
            return candidate.resolve()
    raise RuntimeError(f"Could not locate figure workspace from: {start}")


WORKSPACE_ROOT = find_workspace_root(SCRIPT.parent)
SOURCE_ROOT = WORKSPACE_ROOT / "Figures"
DESTINATION_ROOT = WORKSPACE_ROOT / "manuscript" / "Figures"
ENTRYPOINTS_PATH = (
    WORKSPACE_ROOT / "audit" / "manifests" / "figure_entrypoints.tsv"
)
ARTIFACT_SUFFIXES = {".png", ".pdf", ".svg"}


def main() -> None:
    DESTINATION_ROOT.mkdir(parents=True, exist_ok=True)
    with ENTRYPOINTS_PATH.open(encoding="utf-8", newline="") as handle:
        entrypoints = list(csv.DictReader(handle, delimiter="\t"))
    required = sorted({
        str(Path(row["main_output"]).with_suffix(suffix).name)
        for row in entrypoints
        for suffix in (".png", ".pdf")
    })
    missing = sorted(
        name for name in required if not (SOURCE_ROOT / name).is_file()
    )
    if missing:
        raise RuntimeError("Cannot publish incomplete figure set: " + ", ".join(missing))

    sources = [SOURCE_ROOT / name for name in required]
    for row in entrypoints:
        svg = SOURCE_ROOT / Path(row["main_output"]).with_suffix(".svg").name
        if svg.is_file():
            sources.append(svg)
    sources = sorted(set(sources))

    source_names = {path.name for path in sources}
    for stale in DESTINATION_ROOT.iterdir():
        if (
            stale.is_file()
            and stale.suffix.lower() in ARTIFACT_SUFFIXES
            and stale.name not in source_names
        ):
            stale.unlink()

    for source in sources:
        destination = DESTINATION_ROOT / source.name
        shutil.copy2(source, destination)
    print(f"Published {len(sources)} artifact files to {DESTINATION_ROOT}")


if __name__ == "__main__":
    main()
