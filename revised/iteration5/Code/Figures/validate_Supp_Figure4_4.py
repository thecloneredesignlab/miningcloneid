#!/usr/bin/env python3
"""Independent all-row numeric checks and PDF checks after HPC collection."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pdfplumber
from pypdf import PdfReader


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def near(a, b):
    return abs(a - b) <= 1e-10 * max(1.0, abs(a), abs(b))


def quantile7(values, p):
    values = sorted(values)
    h = (len(values) - 1) * p
    lo, hi = math.floor(h), math.ceil(h)
    return values[lo] + (h - lo) * (values[hi] - values[lo])


def validate(root, audit):
    assert root.name == "iteration5" and audit.is_relative_to(root / "audit")
    data = root / "data/Figures/Supp_Figure4_4"
    metrics = ("ploidy", "effective_o2_pct", "tumor_size_mm3")
    accum = defaultdict(lambda: [0, 0.0, 0.0, 0.0])
    unique = set()
    n_rows = 0
    seed_set = set()
    for row in rows(data / "seed_tumor_trajectories.tsv.gz"):
        seed, day = int(row["seed_number"]), float(row["day"])
        key = (seed, row["cohort"], day)
        identity = (*key, row["harvest"])
        assert identity not in unique
        unique.add(identity)
        seed_set.add(seed)
        acc = accum[key]
        acc[0] += 1
        assert near(float(row["weighted_mean_N"]) / 22, float(row["ploidy"]))
        for index, metric in enumerate(metrics, 1):
            value = float(row[metric])
            assert math.isfinite(value) and value >= 0
            acc[index] += value
        n_rows += 1
    assert seed_set == set(range(1, 501)) and n_rows == 284000
    distribution = defaultdict(list)
    seen = set()
    for row in rows(data / "seed_cohort_trajectories.tsv.gz"):
        key = (int(row["seed_number"]), row["cohort"], float(row["day"]))
        assert key in accum and key not in seen
        seen.add(key)
        n = accum[key][0]
        assert n == int(row["n_tumors"])
        assert n == (4 if key[1] == "2N" or key[2] <= 77 else 2)
        assert int(row["segment"]) == (2 if key[1] == "4N" and key[2] >= 78 else 1)
        for index, metric in enumerate(metrics, 1):
            mean = accum[key][index] / n
            assert near(mean, float(row[metric])), (key, metric)
            distribution[(key[1], key[2], metric)].append(mean)
    assert len(seen) == len(accum) == 74500
    n_summary = 0
    max_error = 0.0
    for row in rows(data / "timecourse_median_iqr.tsv"):
        key = (row["cohort"], float(row["day"]), row["metric"])
        values = distribution.pop(key)
        assert len(values) == int(row["n_seeds"]) == 500
        for field, p in (("q25", 0.25), ("median", 0.5), ("q75", 0.75)):
            expected, actual = quantile7(values, p), float(row[field])
            assert near(expected, actual), (key, field, expected, actual)
            max_error = max(max_error, abs(expected - actual))
        n_summary += 1
    assert n_summary == 447 and not distribution
    seed_checks = list(rows(data / "seed_validation.tsv"))
    assert len(seed_checks) == 500 and all(r["status"] == "validated_cached_dynamic_trajectory" for r in seed_checks)
    provenance = list(rows(data / "source_provenance.tsv"))
    assert len(provenance) == 3000
    assert Counter(int(r["seed_number"]) for r in provenance) == {s: 6 for s in range(1, 501)}

    remote = Path("/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration5")
    checksum_count = 0
    for line in (audit / "output_sha256.txt").read_text().splitlines():
        expected, path = line.split(maxsplit=1)
        local = root / Path(path.strip()).relative_to(remote)
        assert hashlib.sha256(local.read_bytes()).hexdigest() == expected, str(local)
        checksum_count += 1
    assert checksum_count == 11
    pdf = root / "Figures/Supp_Figure4_4.pdf"
    reader = PdfReader(pdf)
    assert len(reader.pages) == 1 and len(list(reader.pages[0].images)) == 0
    assert near(float(reader.pages[0].mediabox.width), 1080)
    assert near(float(reader.pages[0].mediabox.height), 576)
    with pdfplumber.open(pdf) as doc:
        page = doc.pages[0]
        words = page.extract_words()
        assert words
        outside = [w for w in words if w["x0"] < -1 or w["top"] < -1 or w["x1"] > page.width + 1 or w["bottom"] > page.height + 1]
        assert not outside, outside
        text = page.extract_text()
        for required in ("500", "interquartile", "Mean viable ploidy", "Effective O2", "Tumor size", "2N", "4N", "confidence"):
            assert required in text, required
    report = dict(status="PASS", n_seeds=500, n_tumor_time_rows=n_rows,
                  n_cohort_time_rows=len(seen), n_summary_rows=n_summary,
                  quantile_values_independently_checked=n_summary * 3,
                  max_quantile_absolute_error=max_error, checksum_matches=checksum_count,
                  pdf_pages=1, pdf_image_objects=0, pdf_text_within_page=True,
                  visual_render_review="Required separately; not inferred from this script")
    (audit / "local_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    validate(args.root.resolve(), args.audit.resolve())
