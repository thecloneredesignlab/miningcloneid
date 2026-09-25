#!/usr/bin/env python3
"""Compose two labeled raster panels in a fixed one-to-two height ratio."""

from __future__ import annotations

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
TOP_PANEL = required_environment_path("COMPOSITOR_TOP_PANEL")
BOTTOM_PANEL = required_environment_path("COMPOSITOR_BOTTOM_PANEL")
SUPPLEMENTARY_SOURCE = required_environment_path(
    "COMPOSITOR_SUPPLEMENTARY_SOURCE"
)
SUPPLEMENTARY_DESTINATION = required_environment_path(
    "COMPOSITOR_SUPPLEMENTARY_DESTINATION"
)
OUTPUT_BASENAME = os.environ.get("COMPOSITOR_OUTPUT_BASENAME", "").strip()
VALIDATION_BASENAME = os.environ.get(
    "COMPOSITOR_VALIDATION_BASENAME", ""
).strip()
if not OUTPUT_BASENAME or not VALIDATION_BASENAME:
    raise RuntimeError("Missing compositor output or validation basename.")

LABEL_FONT = Path(
    os.environ.get(
        "FIGURE_PANEL_LABEL_FONT",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    )
).expanduser().resolve()
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
CANVAS_WIDTH = 5400
CANVAS_HEIGHT = 5400
ROW_PAD = 12
PANEL_LABEL_BAND = 115
VERTICAL_GAP = 60
PANEL_A_HEIGHT = (CANVAS_HEIGHT - VERTICAL_GAP) // 3
PANEL_B_HEIGHT = 2 * PANEL_A_HEIGHT

if PANEL_A_HEIGHT + VERTICAL_GAP + PANEL_B_HEIGHT != CANVAS_HEIGHT:
    raise RuntimeError("Canvas and slot heights are inconsistent.")


def run(command: list[str], env: dict[str, str] | None = None) -> None:
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


def stage_panels() -> dict[str, Path]:
    PANEL_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    sources = {"a": TOP_PANEL, "b": BOTTOM_PANEL}
    for label, source in sources.items():
        if not source.exists():
            raise RuntimeError(f"Missing panel {label.upper()} source: {source}")
        png_dimensions(source)

    for source in sources.values():
        source_stem = source.stem
        for suffix in (".pdf", ".svg"):
            vector_source = source.parent / f"{source_stem}{suffix}"
            if not vector_source.exists():
                raise RuntimeError(
                    f"Missing vector panel output: {vector_source}"
                )

    for suffix in (".png", ".pdf", ".svg"):
        source = Path(f"{SUPPLEMENTARY_SOURCE}{suffix}")
        if not source.exists():
            raise RuntimeError(
                f"Missing supplementary output: {source}"
            )
        if suffix in (".png", ".pdf"):
            shutil.copy2(
                source,
                Path(f"{SUPPLEMENTARY_DESTINATION}{suffix}"),
            )
    return sources


def label_panel(
    magick: str,
    source: Path,
    destination: Path,
    *,
    label: str,
    env: dict[str, str],
    target_height: int,
) -> None:
    content_height = target_height - PANEL_LABEL_BAND
    if content_height <= 0:
        raise RuntimeError("Panel slot is shorter than its label band.")
    command = [
        magick,
        str(source),
        "-resize",
        f"{CANVAS_WIDTH - 2 * ROW_PAD}x{content_height}",
        "-background",
        "white",
        "-gravity",
        "South",
        "-extent",
        f"{CANVAS_WIDTH}x{target_height}",
        "-gravity",
        "NorthWest",
        "-font",
        str(LABEL_FONT),
        "-pointsize",
        "88",
        "-fill",
        "#111111",
        "-undercolor",
        "#ffffffe6",
        "-annotate",
        "+34+18",
        label,
        str(destination),
    ]
    run(command, env)
    dimensions = png_dimensions(destination)
    if dimensions != (CANVAS_WIDTH, target_height):
        raise RuntimeError(
            f"Unexpected panel {label} slot dimensions: {dimensions}"
        )


def assemble(staged: dict[str, Path]) -> Path:
    magick = shutil.which("magick")
    if magick is None:
        raise RuntimeError("ImageMagick 'magick' was not found on PATH.")
    if not LABEL_FONT.exists():
        raise RuntimeError(f"Panel-label font was not found: {LABEL_FONT}")

    with tempfile.TemporaryDirectory(prefix="stacked_panel_compositor_") as temp_name:
        temp = Path(temp_name)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = str(temp / "cache")
        Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

        labeled_a = temp / "top_labeled.png"
        labeled_b = temp / "bottom_labeled.png"
        label_panel(
            magick,
            staged["a"],
            labeled_a,
            label="A",
            env=env,
            target_height=PANEL_A_HEIGHT,
        )
        label_panel(
            magick,
            staged["b"],
            labeled_b,
            label="B",
            env=env,
            target_height=PANEL_B_HEIGHT,
        )

        row_a_dims = png_dimensions(labeled_a)
        row_b_dims = png_dimensions(labeled_b)

        spacer = temp / "spacer.png"
        run(
            [
                magick,
                "-size",
                f"{CANVAS_WIDTH}x{VERTICAL_GAP}",
                "canvas:white",
                str(spacer),
            ],
            env,
        )

        local_output = PANEL_ROOT / f"{OUTPUT_BASENAME}.png"
        run(
            [
                magick,
                str(labeled_a),
                str(spacer),
                str(labeled_b),
                "-append",
                "+repage",
                str(local_output),
            ],
            env,
        )
        run([magick, "identify", "-quiet", str(local_output)], env)
        if png_dimensions(local_output) != (CANVAS_WIDTH, CANVAS_HEIGHT):
            raise RuntimeError(
                "Final composite does not have the required 1:1 dimensions."
            )

        local_pdf = PANEL_ROOT / f"{OUTPUT_BASENAME}.pdf"
        run(
            [
                magick,
                str(local_output),
                "+repage",
                "-units",
                "PixelsPerInch",
                "-density",
                "300x300",
                "-compress",
                "Zip",
                str(local_pdf),
            ],
            env,
        )
        if not local_pdf.exists() or local_pdf.stat().st_size == 0:
            raise RuntimeError("Final composite PDF was not created.")

        OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
        revised_output = OUTPUT_ROOT / f"{OUTPUT_BASENAME}.png"
        shutil.copy2(local_output, revised_output)
        shutil.copy2(local_pdf, OUTPUT_ROOT / f"{OUTPUT_BASENAME}.pdf")

        VALIDATION_ROOT.mkdir(parents=True, exist_ok=True)
        layout_log = VALIDATION_ROOT / VALIDATION_BASENAME
        header = (
            "assembled_width\tassembled_height\tpanel_a_row_width\t"
            "panel_a_row_height\tpanel_b_row_width\tpanel_b_row_height\t"
            "panel_a_source_width\tpanel_a_source_height\t"
            "panel_b_source_width\tpanel_b_source_height\t"
            "target_width_to_height\ttarget_panel_a_to_b_height\t"
            "vertical_gap_px\n"
        )
        values = (
            *png_dimensions(local_output),
            *row_a_dims,
            *row_b_dims,
            *png_dimensions(staged["a"]),
            *png_dimensions(staged["b"]),
            "1:1",
            "1:2",
            VERTICAL_GAP,
        )
        layout_log.write_text(
            header + "\t".join(map(str, values)) + "\n",
            encoding="utf-8",
        )

    return local_output


def compose_stacked_panels() -> None:
    staged = stage_panels()
    output = assemble(staged)
    width, height = png_dimensions(output)
    print(f"Stacked composite -> {output} ({width}x{height})")
    print(f"Stacked composite PDF -> {output.with_suffix('.pdf')}")


if __name__ == "__main__":
    compose_stacked_panels()
