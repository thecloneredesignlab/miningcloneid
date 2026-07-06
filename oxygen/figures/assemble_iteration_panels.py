#!/usr/bin/env python3
"""Assemble per-panel PNGs into figure-level PNGs.

The script groups panel files named like ``fig3a_description.png`` by figure
number, sorts panels by the letter after the figure number, labels temporary
copies of the panels, and writes ``assembled_fig<N>.png`` outputs directly
under ``oxygen/figures``.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


PANEL_RE = re.compile(r"^fig(?P<figure>\d+)(?P<panel>[a-z]+)[_-].*\.png$", re.IGNORECASE)
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@dataclass(frozen=True)
class Panel:
    figure: int
    panel: str
    path: Path


@dataclass(frozen=True)
class Assembly:
    figure: int
    panels: tuple[Panel, ...]
    output: Path
    width: int
    height: int


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble fig<number><letter> PNG panels into assembled_fig<number>.png files."
    )
    parser.add_argument("folder", type=Path, help="Folder containing panel PNGs.")
    parser.add_argument(
        "--cell-width",
        type=int,
        default=1600,
        help="Width, in pixels, used for each temporary labeled panel before montage.",
    )
    parser.add_argument(
        "--spacing",
        type=int,
        default=48,
        help="Whitespace, in pixels, between panel cells in the assembled image.",
    )
    return parser.parse_args(argv)


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or not header.startswith(PNG_SIGNATURE) or header[12:16] != b"IHDR":
        raise ValueError(f"{path} is not a PNG image with a valid IHDR chunk")
    width, height = struct.unpack(">II", header[16:24])
    if width <= 0 or height <= 0:
        raise ValueError(f"{path} has invalid PNG dimensions {width}x{height}")
    return width, height


def panel_sort_key(panel: Panel) -> tuple[int, str]:
    value = 0
    for char in panel.panel.lower():
        value = value * 26 + (ord(char) - ord("a") + 1)
    return value, panel.panel.lower()


def find_panels(folder: Path) -> dict[int, tuple[Panel, ...]]:
    grouped: dict[int, list[Panel]] = {}
    for path in sorted(folder.iterdir()):
        if not path.is_file():
            continue
        match = PANEL_RE.match(path.name)
        if not match:
            continue
        png_dimensions(path)
        panel = Panel(
            figure=int(match.group("figure")),
            panel=match.group("panel").lower(),
            path=path,
        )
        grouped.setdefault(panel.figure, []).append(panel)

    return {
        figure: tuple(sorted(panels, key=panel_sort_key))
        for figure, panels in sorted(grouped.items())
    }


def choose_columns(panel_count: int) -> int:
    if panel_count <= 1:
        return 1
    return min(panel_count, math.ceil(math.sqrt(panel_count)))


def run_command(command: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        check=False,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        message = [
            f"Command failed with exit code {result.returncode}:",
            shlex.join(command),
        ]
        if result.stdout.strip():
            message.extend(["stdout:", result.stdout.strip()])
        if result.stderr.strip():
            message.extend(["stderr:", result.stderr.strip()])
        raise RuntimeError("\n".join(message))
    return result


def image_magick_env(temp_dir: Path) -> dict[str, str]:
    env = os.environ.copy()
    cache_home = temp_dir / "cache"
    cache_home.mkdir(parents=True, exist_ok=True)
    env.setdefault("XDG_CACHE_HOME", str(cache_home))
    return env


def label_panel(
    magick: str,
    panel: Panel,
    output: Path,
    cell_width: int,
    env: dict[str, str],
) -> None:
    point_size = max(42, round(cell_width * 0.055))
    inset_x = max(18, round(cell_width * 0.022))
    inset_y = max(16, round(cell_width * 0.018))
    label = panel.panel.upper()
    command = [
        magick,
        str(panel.path),
        "-resize",
        f"{cell_width}x",
        "-gravity",
        "NorthWest",
        "-pointsize",
        str(point_size),
        "-fill",
        "#111111",
        "-undercolor",
        "#ffffffe6",
        "-annotate",
        f"+{inset_x}+{inset_y}",
        label,
        str(output),
    ]
    run_command(command, env)


def assemble_figure(
    magick: str,
    figure: int,
    panels: tuple[Panel, ...],
    output_dir: Path,
    temp_dir: Path,
    cell_width: int,
    spacing: int,
    env: dict[str, str],
) -> Assembly:
    labeled_paths: list[Path] = []
    for panel in panels:
        labeled = temp_dir / f"fig{figure}{panel.panel}_labeled.png"
        label_panel(magick, panel, labeled, cell_width, env)
        labeled_paths.append(labeled)

    columns = choose_columns(len(panels))
    output = output_dir / f"assembled_fig{figure}.png"
    command = [
        magick,
        "montage",
        *[str(path) for path in labeled_paths],
        "-tile",
        f"{columns}x",
        "-geometry",
        f"+{spacing}+{spacing}",
        "-background",
        "white",
        str(output),
    ]
    run_command(command, env)

    width, height = png_dimensions(output)
    run_command([magick, "identify", "-quiet", str(output)], env)
    return Assembly(figure=figure, panels=panels, output=output, width=width, height=height)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    folder = args.folder.expanduser().resolve()
    if not folder.is_dir():
        raise SystemExit(f"Input folder does not exist or is not a directory: {folder}")
    if args.cell_width <= 0:
        raise SystemExit("--cell-width must be positive")
    if args.spacing < 0:
        raise SystemExit("--spacing must be non-negative")

    output_dir = Path(__file__).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)

    magick = shutil.which("magick")
    if magick is None:
        raise SystemExit("ImageMagick 'magick' executable was not found on PATH")

    grouped = find_panels(folder)
    if not grouped:
        raise SystemExit(f"No panel PNGs matching {PANEL_RE.pattern!r} found in {folder}")

    with tempfile.TemporaryDirectory(prefix="assemble_iteration_panels_") as temp_name:
        temp_dir = Path(temp_name)
        env = image_magick_env(temp_dir)
        assemblies = [
            assemble_figure(
                magick=magick,
                figure=figure,
                panels=panels,
                output_dir=output_dir,
                temp_dir=temp_dir,
                cell_width=args.cell_width,
                spacing=args.spacing,
                env=env,
            )
            for figure, panels in grouped.items()
        ]

    panel_count = sum(len(assembly.panels) for assembly in assemblies)
    print(
        f"Assembled {len(assemblies)} figures from {panel_count} panels "
        f"in {folder}; outputs written to {output_dir}"
    )
    for assembly in assemblies:
        panel_labels = ", ".join(panel.panel.upper() for panel in assembly.panels)
        print(
            f"- fig{assembly.figure}: {panel_labels} -> "
            f"{assembly.output.name} ({assembly.width}x{assembly.height})"
        )
    print("Validation: all assembled outputs exist and passed PNG/ImageMagick checks.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except RuntimeError as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)
