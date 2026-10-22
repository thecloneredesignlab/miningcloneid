#!/usr/bin/env python3
"""Assemble matrix PDFs from healthy single-page vector exports, avoiding Cairo multipage clipping."""
from pathlib import Path
import argparse
import csv
import shutil
import tempfile
import os
from pypdf import PdfReader, PdfWriter

parser = argparse.ArgumentParser()
parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument("--stem", choices=[f"supp_fig7_{i}{s}" for i in (3, 4) for s in ("", "_ver2")], required=True)
parser.add_argument("--run", type=Path)
parser.add_argument("--audit", type=Path)
args = parser.parse_args()
root = args.workspace.resolve()
run = args.run or root / "data/Figures/Figure7"
audit = args.audit or root / "audit/figure7_matrix_packaging"
with (run / f"{args.stem}_page_manifest.tsv").open(newline="") as stream:
    pages = list(csv.DictReader(stream, delimiter="\t"))
assert [int(p["page"]) for p in pages] == list(range(1, 9))
writer = PdfWriter()
for page in pages:
    source = root / "Figures" / f"{page['matrix_artifact']}.pdf"
    reader = PdfReader(source)
    assert len(reader.pages) == 1
    writer.add_page(reader.pages[0])
target = root / "Figures" / f"{args.stem}.pdf"
backup = audit / "cairo_multipage_before_packaging_fix" / target.name
backup.parent.mkdir(parents=True, exist_ok=True)
if target.exists() and not backup.exists():
    shutil.copy2(target, backup)
with tempfile.NamedTemporaryFile(suffix=".pdf", prefix=f"{args.stem}_", dir=audit, delete=False) as stream:
    temporary = Path(stream.name)
    writer.write(stream)
assert len(PdfReader(temporary).pages) == 8
os.replace(temporary, target)
shutil.copy2(target, root / "manuscript/Figures" / target.name)
print(f"PASS: {args.stem} assembled from eight unchanged single-page vector PDFs")
