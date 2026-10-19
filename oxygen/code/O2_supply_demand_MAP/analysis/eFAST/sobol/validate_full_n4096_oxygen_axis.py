#!/usr/bin/env python3
"""Check physical symlog geometry, file coverage, and PDF text/page boundaries."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import datetime as dt
import json
from pathlib import Path

import numpy as np
import pdfplumber
import redraw_full_n4096_local as h


def validate_geometry():
    fig, ax = h.plt.subplots()
    oxygen = np.linspace(0, 5, 201)
    values = np.arange(402).reshape(2, 201)
    mesh = h.oxygen_mesh(ax, values, oxygen)
    fig.canvas.draw()
    x = lambda a: ax.transData.transform(np.c_[np.asarray(a), np.zeros(len(a))])[:, 0]
    coords = x([0., .0125, .025, .1, 1., 3., 5.])
    assert np.all(np.diff(coords) > 0)
    assert np.isclose(coords[1] - coords[0], .5 * (coords[2] - coords[0]))
    assert np.isclose((coords[5] - coords[4]) / (coords[6] - coords[5]),
                      np.log(3.) / np.log(5. / 3.))
    assert np.array_equal(np.asarray(mesh.get_array()), values)
    assert mesh.get_coordinates()[0, 0, 0] == 0
    assert mesh.get_coordinates()[0, -1, 0] == 5
    assert np.array_equal(ax.get_xticks(), np.array(h.OXYGEN_TICKS))
    assert sorted(float(line.get_xdata()[0]) for line in ax.lines) == [1., 3.]
    assert {t.get_text() for t in ax.texts} == {"Low", "Mid", "High"}
    h.plt.close(fig)
    return dict(status="passed", real_zero_retained=True, low_end_linear=True,
                positive_log_spacing_verified=True, numerical_cell_values_identical=True,
                real_oxygen_cell_edges=True, thresholds=[1., 3.])


def validate_pdf(payload):
    path_string, active, focal_labels = payload
    path = Path(path_string)
    with pdfplumber.open(path) as doc:
        assert len(doc.pages) == 1, path
        page = doc.pages[0]
        text = page.extract_text() or ""
        for phrase in ("0.025", "Low", "Mid", "High", "oxygen"):
            assert phrase in text, (str(path), "missing", phrase)
        if path.name in ("s2_pair_by_oxygen.pdf", "s2_diagnostic_fraction.pdf"):
            labels = [a + " × " + b for a, b in h.itertools.combinations(active, 2)]
        elif "focal" in path.name or path.name == "global_neighborhood_diagnostic_fraction.pdf":
            labels = focal_labels
        else:
            labels = active
        flat = " ".join(text.split())
        for label in labels:
            assert label in flat, (str(path), "missing row label", label)
        words = page.extract_words()
        for band in (w for w in words if w["text"] in ("Low", "Mid", "High")):
            overlaps = [w for w in words if w is not band and
                        min(band["x1"], w["x1"]) - max(band["x0"], w["x0"]) > .5 and
                        min(band["bottom"], w["bottom"]) - max(band["top"], w["top"]) > .5]
            assert not overlaps, (str(path), "region-label text collision", band, overlaps[:1])
        outside = [c for c in page.chars if c.get("text", "").strip() and
                   (c["x0"] < -.5 or c["x1"] > page.width + .5 or
                    c["top"] < -.5 or c["bottom"] > page.height + .5)]
        assert not outside, (str(path), outside[:2])
        return dict(path=str(path), chars=len(page.chars), width=page.width,
                    height=page.height, outside_page=0, pages=1, expected_row_labels_verified=len(labels))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--main-only", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    audit = root / "axis_redraw_20261019"
    geometry = validate_geometry()
    receipts = [json.loads(p.read_text()) for p in (audit / "receipts").glob("*.json")]
    if args.main_only:
        receipts = [r for r in receipts if "/full/" not in r["stem"] and "standalone_N4096" not in r["stem"]]
    expected = json.loads((audit / "plan.json").read_text())["renderer"]
    assert all(r["renderer"] == expected for r in receipts), "Mixed renderer versions"
    files = [root / f["path"] for r in receipts for f in r["files"] if f["path"].endswith(".pdf")]
    assert len(files) == (15 if args.main_only else 1521), len(files)
    mapping = h.read_table(root / "model_structure/parameter_equation_map.tsv")
    active = [r["parameter"] for r in mapping if r["present_fixed_O2"] == "TRUE"]
    focal_labels = [r["display_label"] for r in h.read_table(root / "figure_comparison_20261019/focal_pair_order.tsv")]
    records = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for record in pool.map(validate_pdf, [(str(p), active, focal_labels) for p in sorted(files)]):
            record["path"] = str(Path(record["path"]).relative_to(root))
            records.append(record)
    result = dict(status="passed", checked_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  n_pdfs=len(records), main_only=args.main_only, geometry=geometry,
                  single_renderer_version_verified=True, validation_script_sha256=h.sha(__file__),
                  pdf_text_and_page_boundaries=records)
    filename = "main_pdf_validation.json" if args.main_only else "pdf_validation.json"
    (audit / filename).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(dict(status="passed", n_pdfs=len(records), geometry=geometry), indent=2))


if __name__ == "__main__":
    main()
