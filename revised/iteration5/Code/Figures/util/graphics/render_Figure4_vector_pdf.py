#!/usr/bin/env python3
"""Build only the vector Figure 4 PDF from four already-rendered panels."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from top_with_two_column_bottom_compositor import (
    OUTPUT_BASENAME,
    OUTPUT_ROOT,
    PANEL_ROOT,
    make_vector_pdf,
    stage_panels,
)


def render_vector_pdf() -> None:
    staged = stage_panels()
    with tempfile.TemporaryDirectory(
        prefix="figure4_vector_pdf_", dir=PANEL_ROOT
    ) as temp_name:
        destination = PANEL_ROOT / f"{OUTPUT_BASENAME}.pdf"
        make_vector_pdf(staged, destination, Path(temp_name))
        OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
        shutil.copy2(destination, OUTPUT_ROOT / destination.name)
    print(f"Figure 4 vector PDF -> {OUTPUT_ROOT / destination.name}")


if __name__ == "__main__":
    render_vector_pdf()
