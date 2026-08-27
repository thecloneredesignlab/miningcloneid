#!/usr/bin/env python3
"""Compose Figure 4A-C as three labeled full-width raster panels."""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path


def required_path(name: str) -> Path:
    value = os.environ.get(name, "").strip()
    if not value:
        raise RuntimeError(f"Missing required environment variable: {name}")
    return Path(value).expanduser().resolve()


PANEL_ROOT = required_path("COMPOSITOR_PANEL_DIR")
OUTPUT_ROOT = required_path("COMPOSITOR_OUTPUT_DIR")
VALIDATION_ROOT = required_path("COMPOSITOR_DATA_DIR")
PANELS = {
    "A": required_path("COMPOSITOR_PANEL_A"),
    "B": required_path("COMPOSITOR_PANEL_B"),
    "C": required_path("COMPOSITOR_PANEL_C"),
}
SUPPLEMENTARY_SOURCE = required_path("COMPOSITOR_SUPPLEMENTARY_SOURCE")
SUPPLEMENTARY_DESTINATION = required_path("COMPOSITOR_SUPPLEMENTARY_DESTINATION")
OUTPUT_BASENAME = os.environ.get("COMPOSITOR_OUTPUT_BASENAME", "").strip()
VALIDATION_BASENAME = os.environ.get("COMPOSITOR_VALIDATION_BASENAME", "").strip()
if not OUTPUT_BASENAME or not VALIDATION_BASENAME:
    raise RuntimeError("Missing compositor output or validation basename.")

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
CANVAS_WIDTH = 5400
CANVAS_HEIGHT = 6900
ROW_PAD = 12
PANEL_LABEL_BAND = 115
VERTICAL_GAP = 60
ROW_HEIGHTS = {"A": 1780, "B": 3560, "C": 1440}
LABEL_FONT = Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf")
if sum(ROW_HEIGHTS.values()) + 2 * VERTICAL_GAP != CANVAS_HEIGHT:
    raise RuntimeError("Figure 4 row heights do not match the requested canvas.")


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


def stage_inputs() -> None:
    PANEL_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    for label, source in PANELS.items():
        if not source.exists():
            raise RuntimeError(f"Missing Figure 4{label} source: {source}")
        png_dimensions(source)
        for suffix in (".pdf", ".svg"):
            vector_source = source.with_suffix(suffix)
            if not vector_source.exists():
                raise RuntimeError(f"Missing vector panel output: {vector_source}")

    for suffix in (".png", ".pdf", ".svg"):
        source = Path(f"{SUPPLEMENTARY_SOURCE}{suffix}")
        if not source.exists():
            raise RuntimeError(f"Missing supplementary output: {source}")
        if suffix in (".png", ".pdf"):
            shutil.copy2(source, Path(f"{SUPPLEMENTARY_DESTINATION}{suffix}"))


def label_panel(
    magick: str,
    source: Path,
    destination: Path,
    label: str,
    target_height: int,
    env: dict[str, str],
) -> None:
    content_height = target_height - PANEL_LABEL_BAND
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
    if png_dimensions(destination) != (CANVAS_WIDTH, target_height):
        raise RuntimeError(f"Unexpected labeled Figure 4{label} dimensions.")


def assemble() -> Path:
    magick = shutil.which("magick")
    if magick is None:
        raise RuntimeError("ImageMagick 'magick' was not found on PATH.")
    if not LABEL_FONT.exists():
        raise RuntimeError(f"Panel-label font was not found: {LABEL_FONT}")

    with tempfile.TemporaryDirectory(prefix="figure4_three_panel_") as temp_name:
        temp = Path(temp_name)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = str(temp / "cache")
        Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

        labeled: dict[str, Path] = {}
        for label, source in PANELS.items():
            destination = temp / f"panel_{label.lower()}_labeled.png"
            label_panel(
                magick,
                source,
                destination,
                label,
                ROW_HEIGHTS[label],
                env,
            )
            labeled[label] = destination

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

        local_png = PANEL_ROOT / f"{OUTPUT_BASENAME}.png"
        run(
            [
                magick,
                str(labeled["A"]),
                str(spacer),
                str(labeled["B"]),
                str(spacer),
                str(labeled["C"]),
                "-append",
                "+repage",
                str(local_png),
            ],
            env,
        )
        if png_dimensions(local_png) != (CANVAS_WIDTH, CANVAS_HEIGHT):
            raise RuntimeError("Final Figure 4 composite dimensions are incorrect.")

        local_pdf = PANEL_ROOT / f"{OUTPUT_BASENAME}.pdf"
        run(
            [
                magick,
                str(local_png),
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
        if not local_pdf.exists() or not local_pdf.stat().st_size:
            raise RuntimeError("Final Figure 4 PDF was not created.")

        shutil.copy2(local_png, OUTPUT_ROOT / f"{OUTPUT_BASENAME}.png")
        shutil.copy2(local_pdf, OUTPUT_ROOT / f"{OUTPUT_BASENAME}.pdf")

        VALIDATION_ROOT.mkdir(parents=True, exist_ok=True)
        validation = VALIDATION_ROOT / VALIDATION_BASENAME
        header = (
            "assembled_width\tassembled_height\t"
            "panel_a_row_height\tpanel_b_row_height\tpanel_c_row_height\t"
            "panel_a_source_width\tpanel_a_source_height\t"
            "panel_b_source_width\tpanel_b_source_height\t"
            "panel_c_source_width\tpanel_c_source_height\t"
            "panel_label_band_px\tvertical_gap_px\tpanel_order\n"
        )
        values = (
            CANVAS_WIDTH,
            CANVAS_HEIGHT,
            ROW_HEIGHTS["A"],
            ROW_HEIGHTS["B"],
            ROW_HEIGHTS["C"],
            *png_dimensions(PANELS["A"]),
            *png_dimensions(PANELS["B"]),
            *png_dimensions(PANELS["C"]),
            PANEL_LABEL_BAND,
            VERTICAL_GAP,
            "A:B:C",
        )
        validation.write_text(
            header + "\t".join(map(str, values)) + "\n",
            encoding="utf-8",
        )

    return local_png


def main() -> None:
    stage_inputs()
    output = assemble()
    print(f"Three-panel Figure 4 -> {output} ({CANVAS_WIDTH}x{CANVAS_HEIGHT})")
    print(f"Three-panel Figure 4 PDF -> {output.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
