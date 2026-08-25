#!/usr/bin/env python3
"""Compose Figure 4 as A above B beside a stacked C/D right column."""

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
BOTTOM_LEFT_PANEL = required_environment_path("COMPOSITOR_BOTTOM_LEFT_PANEL")
BOTTOM_RIGHT_TOP_PANEL = required_environment_path(
    "COMPOSITOR_BOTTOM_RIGHT_TOP_PANEL"
)
BOTTOM_RIGHT_BOTTOM_PANEL = required_environment_path(
    "COMPOSITOR_BOTTOM_RIGHT_BOTTOM_PANEL"
)
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

LABEL_FONT = REPO_ROOT / "data" / "fonts" / "Arial Bold.ttf"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
CANVAS_WIDTH = 5400
CANVAS_HEIGHT = 5400
ROW_PAD = 12
PANEL_LABEL_BAND = 115
VERTICAL_GAP = 60
HORIZONTAL_GAP = 40
RIGHT_STACK_GAP = 40
PANEL_A_HEIGHT = (CANVAS_HEIGHT - VERTICAL_GAP) // 3
BOTTOM_ROW_HEIGHT = 2 * PANEL_A_HEIGHT
PANEL_B_WIDTH = 4000
RIGHT_COLUMN_WIDTH = CANVAS_WIDTH - PANEL_B_WIDTH - HORIZONTAL_GAP
PANEL_C_HEIGHT = RIGHT_COLUMN_WIDTH - 2 * ROW_PAD + PANEL_LABEL_BAND
PANEL_D_HEIGHT = BOTTOM_ROW_HEIGHT - PANEL_C_HEIGHT - RIGHT_STACK_GAP
PANEL_TITLES = {
    "A": "In vivo cohort dynamics and terminal ploidy",
    "B": "Oxygen-dependent parameter-ploidy landscape",
    "C": "Exploratory fitted-endpoint landscape",
    "D": "Top six cluster-separating fitted parameters",
}
PANEL_LABEL_POINTSIZE = 88
PANEL_TITLE_POINTSIZE = {"A": 88, "B": 88, "C": 58, "D": 58}
PANEL_LABEL_X = 34
PANEL_TITLE_X = 140
PANEL_LABEL_Y = 18

if PANEL_A_HEIGHT + VERTICAL_GAP + BOTTOM_ROW_HEIGHT != CANVAS_HEIGHT:
    raise RuntimeError("Canvas and row heights are inconsistent.")
if PANEL_B_WIDTH + HORIZONTAL_GAP + RIGHT_COLUMN_WIDTH != CANVAS_WIDTH:
    raise RuntimeError("Canvas and bottom-row widths are inconsistent.")
if PANEL_C_HEIGHT + RIGHT_STACK_GAP + PANEL_D_HEIGHT != BOTTOM_ROW_HEIGHT:
    raise RuntimeError("Right-column panel heights are inconsistent.")


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
    sources = {
        "a": TOP_PANEL,
        "b": BOTTOM_LEFT_PANEL,
        "c": BOTTOM_RIGHT_TOP_PANEL,
        "d": BOTTOM_RIGHT_BOTTOM_PANEL,
    }
    for label, source in sources.items():
        if not source.exists():
            raise RuntimeError(f"Missing panel {label.upper()} source: {source}")
        png_dimensions(source)
        for suffix in (".pdf", ".svg"):
            vector_source = source.parent / f"{source.stem}{suffix}"
            if not vector_source.exists():
                raise RuntimeError(
                    f"Missing vector panel output: {vector_source}"
                )

    for suffix in (".png", ".pdf", ".svg"):
        source = Path(f"{SUPPLEMENTARY_SOURCE}{suffix}")
        if not source.exists():
            raise RuntimeError(f"Missing supplementary output: {source}")
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
    title: str,
    title_pointsize: int,
    target_width: int,
    target_height: int,
    env: dict[str, str],
) -> None:
    content_width = target_width - 2 * ROW_PAD
    content_height = target_height - PANEL_LABEL_BAND
    if content_width <= 0 or content_height <= 0:
        raise RuntimeError("Panel slot is smaller than its padding.")
    title_y = PANEL_LABEL_Y + round(
        (PANEL_LABEL_POINTSIZE - title_pointsize) / 2
    )
    run(
        [
            magick,
            "-size",
            f"{target_width}x{target_height}",
            "canvas:white",
            "(",
            str(source),
            "-resize",
            f"{content_width}x{content_height}",
            ")",
            "-gravity",
            "North",
            "-geometry",
            f"+0+{PANEL_LABEL_BAND}",
            "-composite",
            "-gravity",
            "NorthWest",
            "-font",
            str(LABEL_FONT),
            "-pointsize",
            str(PANEL_LABEL_POINTSIZE),
            "-fill",
            "#111111",
            "-undercolor",
            "#ffffffe6",
            "-annotate",
            f"+{PANEL_LABEL_X}+{PANEL_LABEL_Y}",
            f"{label}.",
            "-pointsize",
            str(title_pointsize),
            "-annotate",
            f"+{PANEL_TITLE_X}+{title_y}",
            title,
            str(destination),
        ],
        env,
    )
    dimensions = png_dimensions(destination)
    if dimensions != (target_width, target_height):
        raise RuntimeError(
            f"Unexpected panel {label} slot dimensions: {dimensions}"
        )


def assemble(staged: dict[str, Path]) -> Path:
    magick = shutil.which("magick")
    if magick is None:
        raise RuntimeError("ImageMagick 'magick' was not found on PATH.")
    if not LABEL_FONT.exists():
        raise RuntimeError(f"Panel-label font was not found: {LABEL_FONT}")

    with tempfile.TemporaryDirectory(
        prefix="figure4_four_panel_compositor_"
    ) as temp_name:
        temp = Path(temp_name)
        env = os.environ.copy()
        env["XDG_CACHE_HOME"] = str(temp / "cache")
        Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

        labeled_a = temp / "panel_a.png"
        labeled_b = temp / "panel_b.png"
        labeled_c = temp / "panel_c.png"
        labeled_d = temp / "panel_d.png"
        label_panel(
            magick,
            staged["a"],
            labeled_a,
            label="A",
            title=PANEL_TITLES["A"],
            title_pointsize=PANEL_TITLE_POINTSIZE["A"],
            target_width=CANVAS_WIDTH,
            target_height=PANEL_A_HEIGHT,
            env=env,
        )
        label_panel(
            magick,
            staged["b"],
            labeled_b,
            label="B",
            title=PANEL_TITLES["B"],
            title_pointsize=PANEL_TITLE_POINTSIZE["B"],
            target_width=PANEL_B_WIDTH,
            target_height=BOTTOM_ROW_HEIGHT,
            env=env,
        )
        label_panel(
            magick,
            staged["c"],
            labeled_c,
            label="C",
            title=PANEL_TITLES["C"],
            title_pointsize=PANEL_TITLE_POINTSIZE["C"],
            target_width=RIGHT_COLUMN_WIDTH,
            target_height=PANEL_C_HEIGHT,
            env=env,
        )
        label_panel(
            magick,
            staged["d"],
            labeled_d,
            label="D",
            title=PANEL_TITLES["D"],
            title_pointsize=PANEL_TITLE_POINTSIZE["D"],
            target_width=RIGHT_COLUMN_WIDTH,
            target_height=PANEL_D_HEIGHT,
            env=env,
        )

        horizontal_spacer = temp / "horizontal_spacer.png"
        vertical_spacer = temp / "vertical_spacer.png"
        right_stack_spacer = temp / "right_stack_spacer.png"
        run(
            [
                magick,
                "-size",
                f"{HORIZONTAL_GAP}x{BOTTOM_ROW_HEIGHT}",
                "canvas:white",
                str(horizontal_spacer),
            ],
            env,
        )
        run(
            [
                magick,
                "-size",
                f"{CANVAS_WIDTH}x{VERTICAL_GAP}",
                "canvas:white",
                str(vertical_spacer),
            ],
            env,
        )
        run(
            [
                magick,
                "-size",
                f"{RIGHT_COLUMN_WIDTH}x{RIGHT_STACK_GAP}",
                "canvas:white",
                str(right_stack_spacer),
            ],
            env,
        )
        right_column = temp / "right_column.png"
        run(
            [
                magick,
                str(labeled_c),
                str(right_stack_spacer),
                str(labeled_d),
                "-append",
                "+repage",
                str(right_column),
            ],
            env,
        )
        if png_dimensions(right_column) != (
            RIGHT_COLUMN_WIDTH,
            BOTTOM_ROW_HEIGHT,
        ):
            raise RuntimeError("Stacked C/D column dimensions are inconsistent.")
        bottom_row = temp / "bottom_row.png"
        run(
            [
                magick,
                str(labeled_b),
                str(horizontal_spacer),
                str(right_column),
                "+append",
                "+repage",
                str(bottom_row),
            ],
            env,
        )

        local_output = PANEL_ROOT / f"{OUTPUT_BASENAME}.png"
        run(
            [
                magick,
                str(labeled_a),
                str(vertical_spacer),
                str(bottom_row),
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
        shutil.copy2(local_output, OUTPUT_ROOT / f"{OUTPUT_BASENAME}.png")
        shutil.copy2(local_pdf, OUTPUT_ROOT / f"{OUTPUT_BASENAME}.pdf")

        layout_log = VALIDATION_ROOT / VALIDATION_BASENAME
        header = (
            "assembled_width\tassembled_height\tpanel_a_row_height\t"
            "bottom_row_height\tpanel_b_width\tright_column_width\t"
            "panel_c_height\tpanel_d_height\t"
            "panel_a_source_width\tpanel_a_source_height\t"
            "panel_b_source_width\tpanel_b_source_height\t"
            "panel_c_source_width\tpanel_c_source_height\t"
            "panel_d_source_width\tpanel_d_source_height\t"
            "target_width_to_height\ttarget_panel_a_to_bottom_height\t"
            "panel_b_c_top_aligned\tpanel_c_d_combined_height_equals_b\t"
            "vertical_gap_px\thorizontal_gap_px\tright_stack_gap_px\t"
            "panel_label_pointsize\tpanel_a_title_pointsize\t"
            "panel_b_title_pointsize\tpanel_c_title_pointsize\t"
            "panel_d_title_pointsize\t"
            "panel_heading_format\tpanel_a_title\tpanel_b_title\t"
            "panel_c_title\tpanel_d_title\t"
            "panel_order\n"
        )
        values = (
            *png_dimensions(local_output),
            PANEL_A_HEIGHT,
            BOTTOM_ROW_HEIGHT,
            PANEL_B_WIDTH,
            RIGHT_COLUMN_WIDTH,
            PANEL_C_HEIGHT,
            PANEL_D_HEIGHT,
            *png_dimensions(staged["a"]),
            *png_dimensions(staged["b"]),
            *png_dimensions(staged["c"]),
            *png_dimensions(staged["d"]),
            "1:1",
            "1:2",
            "TRUE",
            "TRUE",
            VERTICAL_GAP,
            HORIZONTAL_GAP,
            RIGHT_STACK_GAP,
            PANEL_LABEL_POINTSIZE,
            PANEL_TITLE_POINTSIZE["A"],
            PANEL_TITLE_POINTSIZE["B"],
            PANEL_TITLE_POINTSIZE["C"],
            PANEL_TITLE_POINTSIZE["D"],
            "<letter>. <panel title>",
            PANEL_TITLES["A"],
            PANEL_TITLES["B"],
            PANEL_TITLES["C"],
            PANEL_TITLES["D"],
            "A:(B|(C/D))",
        )
        layout_log.write_text(
            header + "\t".join(map(str, values)) + "\n",
            encoding="utf-8",
        )

    return local_output


def compose_figure4_panels() -> None:
    staged = stage_panels()
    output = assemble(staged)
    width, height = png_dimensions(output)
    print(f"Figure 4 composite -> {output} ({width}x{height})")
    print(f"Figure 4 composite PDF -> {output.with_suffix('.pdf')}")


if __name__ == "__main__":
    compose_figure4_panels()
