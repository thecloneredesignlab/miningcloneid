#!/usr/bin/env python3
"""Read-only word-level PDF QA using Poppler; report only under iteration4."""
import argparse
import csv
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--report", type=Path, required=True)
parser.add_argument("pdf", type=Path, nargs="+")
args = parser.parse_args()
poppler = shutil.which("pdftotext")
if not poppler:
    candidate = Path("/opt/homebrew/bin/pdftotext")
    poppler = str(candidate) if candidate.is_file() else None
if not poppler:
    raise SystemExit("pdftotext is required")
rows = []
for pdf in args.pdf:
    tree = ET.fromstring(subprocess.check_output([poppler, "-bbox", str(pdf), "-"]))
    words = []
    outside = 0
    for page in tree.iter():
        if page.tag.rsplit("}", 1)[-1] != "page":
            continue
        width, height = float(page.attrib["width"]), float(page.attrib["height"])
        for word in page.iter():
            if word.tag.rsplit("}", 1)[-1] != "word":
                continue
            words.append("".join(word.itertext()))
            bounds = [float(word.attrib[key]) for key in ("xMin", "yMin", "xMax", "yMax")]
            outside += bounds[0] < -.5 or bounds[1] < -.5 or bounds[2] > width + .5 or bounds[3] > height + .5
    lexical = [re.sub(r"[^a-z]", "", word.lower()) for word in words]
    full_words = all(word in lexical for word in ("oxygen", "ploidy"))
    if pdf.name.startswith("assembled_fig7"):
        full_words = full_words and all(word in lexical for word in ("experimental", "time", "stochastic"))
    row = {"pdf": str(pdf), "word_count": len(words), "whole_words_present": full_words,
           "outside_page_words": outside, "passed": full_words and outside == 0}
    rows.append(row)
args.report.parent.mkdir(parents=True, exist_ok=True)
with args.report.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
for row in rows:
    print(row)
if not all(row["passed"] for row in rows):
    raise SystemExit("PDF word-level validation failed")
