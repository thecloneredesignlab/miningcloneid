#!/usr/bin/env python3
"""Streaming independent validation of all 500 long-horizon forecasts."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import pdfplumber
from pypdf import PdfReader


def rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def near(a, b):
    return abs(a - b) <= 1e-10 * max(1, abs(a), abs(b))


def q7(x, p):
    x = sorted(x)
    h = (len(x) - 1) * p
    lo, hi = math.floor(h), math.ceil(h)
    return x[lo] + (h - lo) * (x[hi] - x[lo])


def validate(root, audit):
    assert root.name == "iteration5" and audit.is_relative_to(root / "audit")
    data = root / "data/Figures/Supp_Figure4_4_1000days"
    original = root / "data/Figures/Supp_Figure4_4"
    metrics = ("ploidy", "effective_o2_pct", "tumor_size_mm3")
    seedkey = lambda r: int(r["seed_number"])
    tumor_groups = itertools.groupby(rows(data / "seed_tumor_trajectories.tsv.gz"), seedkey)
    cohort_groups = iter(itertools.groupby(rows(data / "seed_cohort_trajectories.tsv.gz"), seedkey))
    short_groups = iter(itertools.groupby(rows(original / "seed_tumor_trajectories.tsv.gz"), seedkey))
    distributions = defaultdict(list)
    all_seeds = []
    n_tumor = n_cohort = n_overlap = 0
    max_overlap_error = 0.0
    for seed, block in tumor_groups:
        all_seeds.append(seed)
        short_seed, short_block = next(short_groups)
        assert short_seed == seed
        ref = {(r["harvest"], r["cohort"], float(r["day"])): r for r in short_block}
        assert len(ref) == 568
        sums = defaultdict(lambda: [0, 0.0, 0.0, 0.0])
        seen = set()
        for row in block:
            day = float(row["day"])
            identity = (row["harvest"], row["cohort"], day)
            assert identity not in seen and day.is_integer() and 0 <= day <= 1000
            seen.add(identity)
            key = (row["cohort"], day)
            acc = sums[key]
            acc[0] += 1
            assert near(float(row["weighted_mean_N"]) / 22, float(row["ploidy"]))
            for index, metric in enumerate(metrics, 1):
                value = float(row[metric])
                assert math.isfinite(value) and value >= 0
                acc[index] += value
            if identity in ref:
                old = ref.pop(identity)
                for metric in metrics:
                    a, b = float(row[metric]), float(old[metric])
                    error = abs(a - b) / max(1, abs(b))
                    assert error <= 1e-10
                    max_overlap_error = max(max_overlap_error, error)
                n_overlap += 1
            n_tumor += 1
        assert len(seen) == 8008 and len(sums) == 2002 and not ref
        cohort_seed, cohort_block = next(cohort_groups)
        assert cohort_seed == seed
        for row in cohort_block:
            key = (row["cohort"], float(row["day"]))
            acc = sums.pop(key)
            assert acc[0] == int(row["n_tumors"]) == 4 and int(row["segment"]) == 1
            for index, metric in enumerate(metrics, 1):
                value = acc[index] / 4
                assert near(value, float(row[metric]))
                distributions[(*key, metric)].append(value)
            n_cohort += 1
        assert not sums
        if seed % 100 == 0:
            print(f"Independently checked {seed}/500 seeds", flush=True)
    assert all_seeds == list(range(1, 501)) and n_tumor == 4004000 and n_cohort == 1001000 and n_overlap == 284000
    assert next(cohort_groups, None) is None and next(short_groups, None) is None
    diagnostics = {(r["cohort"], float(r["day"])): r for r in rows(data / "ploidy_time_window_comparison.tsv")}
    n_summary = 0
    max_quantile_error = 0.0
    for row in rows(data / "timecourse_median_iqr.tsv"):
        key = (row["cohort"], float(row["day"]), row["metric"])
        values = distributions.pop(key)
        assert len(values) == int(row["n_seeds"]) == 500
        for field, p in (("q25", .25), ("median", .5), ("q75", .75)):
            a, b = q7(values, p), float(row[field])
            assert near(a, b), (key, field)
            max_quantile_error = max(max_quantile_error, abs(a - b) / max(1, abs(a)))
        if key[2] == "ploidy" and key[:2] in diagnostics:
            diag = diagnostics.pop(key[:2])
            for field, p in (("q05", .05), ("q25", .25), ("median", .5), ("q75", .75), ("q95", .95)):
                assert near(q7(values, p), float(diag[field]))
            assert near(min(values), float(diag["minimum"])) and near(max(values), float(diag["maximum"]))
        n_summary += 1
    assert n_summary == 6006 and not distributions and not diagnostics
    source_rows = list(rows(data / "source_provenance.tsv"))
    assert len(source_rows) == 3005
    remote = Path("/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration5")
    checks = 0
    for line in (audit / "output_sha256.txt").read_text().splitlines():
        expected, path = line.split(maxsplit=1)
        local = root / Path(path.strip()).relative_to(remote)
        assert hashlib.sha256(local.read_bytes()).hexdigest() == expected
        checks += 1
    assert checks == 13
    original_checks = 0
    for line in (audit / "protected_before.sha256").read_text().splitlines():
        expected, path = line.split(maxsplit=1)
        relative = Path(path.strip()).relative_to(remote)
        if str(relative).startswith("data/Figures/Supp_Figure4_4/") or (relative.parent == Path("Figures") and relative.stem == "Supp_Figure4_4"):
            assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == expected
            original_checks += 1
    assert original_checks == 11
    pdf = root / "Figures/Supp_Figure4_4_1000days.pdf"
    reader = PdfReader(pdf)
    assert len(reader.pages) == 1 and len(list(reader.pages[0].images)) == 0
    assert near(float(reader.pages[0].mediabox.width), 1080) and near(float(reader.pages[0].mediabox.height), 576)
    with pdfplumber.open(pdf) as doc:
        page = doc.pages[0]
        for w in page.extract_words():
            assert w["x0"] >= -1 and w["top"] >= -1 and w["x1"] <= page.width + 1 and w["bottom"] <= page.height + 1, w
        text = page.extract_text()
        for value in ("1,000", "500", "extrapolation", "pseudo-log10", "2N", "4N", "confidence"):
            assert value in text, value
    report = dict(status="PASS", n_seeds=500, tumor_rows=n_tumor, cohort_rows=n_cohort, summary_rows=n_summary,
                  independently_checked_quantiles=n_summary * 3, original_overlap_rows=n_overlap,
                  max_original_overlap_relative_error=max_overlap_error, max_quantile_relative_error=max_quantile_error,
                  checksum_matches=checks, original_short_window_files_unchanged=original_checks,
                  pdf_pages=1, pdf_image_objects=0, pdf_text_within_page=True,
                  visual_review="Required separately")
    (audit / "local_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()
    validate(args.root.resolve(), args.audit.resolve())
