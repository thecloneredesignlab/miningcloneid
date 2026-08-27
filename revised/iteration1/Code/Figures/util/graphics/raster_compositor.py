#!/usr/bin/env python3
"""Compose six labeled panels in a reusable three-row raster layout.

Layout:
  row 1: A (50%) | B (50%)
  row 2: C (50%) | D (25%) | E (25%)
  row 3: F (100%)
"""

from __future__ import annotations

import hashlib
import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path


def required_environment_path(name: str) -> Path:
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"Missing required environment variable: {name}")
    return Path(value).expanduser().resolve()


PANEL_ROOT = required_environment_path("COMPOSITOR_PANEL_DIR")
OUTPUT_ROOT = required_environment_path("COMPOSITOR_OUTPUT_DIR")
VALIDATION_ROOT = required_environment_path("COMPOSITOR_DATA_DIR")
REPO_ROOT = required_environment_path("HYPOXIA_REPO_ROOT")
OUTPUT_BASENAME = os.environ.get("COMPOSITOR_OUTPUT_BASENAME", "").strip()
VALIDATION_BASENAME = os.environ.get(
    "COMPOSITOR_VALIDATION_BASENAME", ""
).strip()
if not OUTPUT_BASENAME or not VALIDATION_BASENAME:
    raise RuntimeError("Incomplete compositor configuration.")
PANEL_SOURCES = {
    label: required_environment_path(f"COMPOSITOR_PANEL_{label}")
    for label in "ABCDEF"
}
LABEL_FONT = REPO_ROOT / "data" / "fonts" / "Arial Bold.ttf"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

FULL_WIDTH = 6000
HALF_WIDTH = 3000
QUARTER_WIDTH = 1500
HORIZONTAL_GAP = 50
VERTICAL_GAP = 70
LABEL_PADDING = 110

def run(command: list[str], env: dict[str, str]) -> None:
    result = subprocess.run(
        command,
        check=False,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stderr.strip()}"
        )


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if (
        len(header) < 24
        or not header.startswith(PNG_SIGNATURE)
        or header[12:16] != b"IHDR"
    ):
        raise RuntimeError(f"Invalid PNG: {path}")
    return struct.unpack(">II", header[16:24])


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_spacer(
    magick: str,
    destination: Path,
    width: int,
    height: int,
    env: dict[str, str],
) -> None:
    run(
        [
            magick,
            "-size",
            f"{width}x{height}",
            "xc:white",
            str(destination),
        ],
        env,
    )


def prepare_panel(
    magick: str,
    source: Path,
    destination: Path,
    *,
    label: str,
    width: int,
    env: dict[str, str],
) -> tuple[int, int]:
    if not source.exists():
        raise RuntimeError(f"Missing panel {label} source: {source}")
    png_dimensions(source)
    run(
        [
            magick,
            str(source),
            "-resize",
            f"{width}x",
            "-background",
            "white",
            "-gravity",
            "North",
            "-splice",
            f"0x{LABEL_PADDING}",
            "-gravity",
            "NorthWest",
            "-font",
            str(LABEL_FONT),
            "-pointsize",
            "112",
            "-fill",
            "#111111",
            "-annotate",
            "+24+9",
            label,
            str(destination),
        ],
        env,
    )
    return png_dimensions(destination)


def horizontal_row(
    magick: str,
    panels: list[Path],
    destination: Path,
    *,
    gap: int,
    env: dict[str, str],
    temp: Path,
) -> tuple[int, int]:
    heights = [png_dimensions(panel)[1] for panel in panels]
    target_height = max(heights)
    normalized: list[Path] = []
    for index, panel in enumerate(panels):
        width, _ = png_dimensions(panel)
        normalized_panel = temp / f"{destination.stem}_panel_{index}.png"
        run(
            [
                magick,
                str(panel),
                "-background",
                "white",
                "-gravity",
                "North",
                "-extent",
                f"{width}x{target_height}",
                str(normalized_panel),
            ],
            env,
        )
        normalized.append(normalized_panel)

    append_items: list[Path] = []
    for index, panel in enumerate(normalized):
        if index:
            spacer = temp / f"{destination.stem}_gap_{index}.png"
            make_spacer(magick, spacer, gap, target_height, env)
            append_items.append(spacer)
        append_items.append(panel)
    run(
        [
            magick,
            *[str(item) for item in append_items],
            "+append",
            str(destination),
        ],
        env,
    )
    return png_dimensions(destination)


def normalize_width(
    magick: str,
    source: Path,
    destination: Path,
    width: int,
    env: dict[str, str],
) -> tuple[int, int]:
    _, height = png_dimensions(source)
    run(
        [
            magick,
            str(source),
            "-background",
            "white",
            "-gravity",
            "North",
            "-extent",
            f"{width}x{height}",
            str(destination),
        ],
        env,
    )
    return png_dimensions(destination)


def compose_three_row_raster() -> Path:
    magick = shutil.which("magick")
    if magick is None:
        raise RuntimeError("ImageMagick 'magick' was not found on PATH.")
    if not LABEL_FONT.exists():
        raise RuntimeError(f"Panel-label font was not found: {LABEL_FONT}")

    PANEL_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="three_row_compositor_") as name:
        temp = Path(name)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = str(temp / "cache")
        Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

        widths = {
            "A": HALF_WIDTH,
            "B": HALF_WIDTH,
            "C": HALF_WIDTH,
            "D": QUARTER_WIDTH,
            "E": QUARTER_WIDTH,
            "F": FULL_WIDTH,
        }
        prepared: dict[str, Path] = {}
        for label, source in PANEL_SOURCES.items():
            destination = temp / f"panel_{label.lower()}_labeled.png"
            prepare_panel(
                magick,
                source,
                destination,
                label=label,
                width=widths[label],
                env=env,
            )
            prepared[label] = destination

        row1 = temp / "row1.png"
        row2 = temp / "row2.png"
        horizontal_row(
            magick,
            [prepared["A"], prepared["B"]],
            row1,
            gap=HORIZONTAL_GAP,
            env=env,
            temp=temp,
        )
        horizontal_row(
            magick,
            [prepared["C"], prepared["D"], prepared["E"]],
            row2,
            gap=HORIZONTAL_GAP,
            env=env,
            temp=temp,
        )

        canvas_width = max(
            png_dimensions(row1)[0],
            png_dimensions(row2)[0],
            png_dimensions(prepared["F"])[0],
        )
        normalized_rows: list[Path] = []
        for row_index, row in enumerate((row1, row2, prepared["F"]), start=1):
            normalized = temp / f"row{row_index}_normalized.png"
            normalize_width(magick, row, normalized, canvas_width, env)
            normalized_rows.append(normalized)

        vertical_items: list[Path] = []
        for index, row in enumerate(normalized_rows):
            if index:
                spacer = temp / f"vertical_gap_{index}.png"
                make_spacer(
                    magick,
                    spacer,
                    canvas_width,
                    VERTICAL_GAP,
                    env,
                )
                vertical_items.append(spacer)
            vertical_items.append(row)

        output = PANEL_ROOT / f"{OUTPUT_BASENAME}.png"
        run(
            [
                magick,
                *[str(item) for item in vertical_items],
                "-append",
                str(output),
            ],
            env,
        )
        run([magick, "identify", "-quiet", str(output)], env)

        pdf_output = PANEL_ROOT / f"{OUTPUT_BASENAME}.pdf"
        run(
            [
                magick,
                str(output),
                "+repage",
                "-units",
                "PixelsPerInch",
                "-density",
                "300x300",
                "-compress",
                "Zip",
                str(pdf_output),
            ],
            env,
        )
        if not pdf_output.exists() or pdf_output.stat().st_size == 0:
            raise RuntimeError("Final composite PDF was not created.")

    shutil.copy2(output, OUTPUT_ROOT / output.name)
    shutil.copy2(pdf_output, OUTPUT_ROOT / pdf_output.name)
    for label, source in PANEL_SOURCES.items():
        for suffix in (".pdf", ".svg"):
            vector_source = source.with_suffix(suffix)
            if not vector_source.exists():
                raise RuntimeError(
                    f"Missing vector output for panel {label}: {vector_source}"
                )

    width, height = png_dimensions(output)
    VALIDATION_ROOT.mkdir(parents=True, exist_ok=True)
    validation = VALIDATION_ROOT / VALIDATION_BASENAME
    rows = [
        ("metric", "value"),
        ("assembled_width_px", str(width)),
        ("assembled_height_px", str(height)),
        ("row1_layout", "A:50%;B:50%"),
        ("row2_layout", "C:50%;D:25%;E:25%"),
        ("row3_layout", "F:100%"),
        ("panel_d_width_relative_to_c", "0.5"),
        ("panel_e_width_relative_to_c", "0.5"),
        ("panel_f_contains", "2N-left;4N-right"),
        ("output_md5", md5(output)),
    ]
    validation.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    print(f"Three-row composite -> {output} ({width}x{height})")
    print(f"Three-row composite PDF -> {pdf_output}")
    return output


if __name__ == "__main__":
    compose_three_row_raster()
