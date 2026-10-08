#!/usr/bin/env python3
"""Per-fit-seed eFAST neighborhoods, resumable pilots, and seed-level summaries."""
import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
from types import SimpleNamespace

import numpy as np
import efast

HERE = Path(__file__).resolve().parent
REPEATS = 5
INDICES = ("S1", "ST")
SOURCES = ("figure4b_spearman_source.tsv", "figure4b_o2_classification_source.tsv",
           "figure4b_parameter_ranking_source.tsv", "figure4_parameter_groups_source.tsv",
           "figure4_parameter_group_palette_source.tsv")


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def neighborhood_bounds(ranges, values, fraction=0.1):
    """Natural span first; only then encode log10 bounds. Never perturb log span."""
    rows = []
    for spec in ranges:
        name = spec["parameter"]
        lo, hi = float(spec["natural_lower"]), float(spec["natural_upper"])
        center = float(values[name])
        if not math.isfinite(center) or center < lo - 1e-12 or center > hi + 1e-12:
            raise ValueError("Best value outside global bounds: " + name)
        center = min(hi, max(lo, center))  # only round-off at a bound
        radius = fraction * (hi - lo)
        lower, upper = max(lo, center - radius), min(hi, center + radius)
        if not lower < upper or (spec["transform"] == "log10" and lower <= 0):
            raise ValueError("Invalid neighborhood for " + name)
        row = dict(spec)
        row.update(best_value=center, global_natural_lower=lo, global_natural_upper=hi,
                   requested_half_width=radius, fraction_of_global_span=fraction,
                   lower_clipped=str(center - radius < lo).upper(),
                   upper_clipped=str(center + radius > hi).upper(),
                   natural_lower=lower, natural_upper=upper,
                   encoded_lower=math.log10(lower) if spec["transform"] == "log10" else lower,
                   encoded_upper=math.log10(upper) if spec["transform"] == "log10" else upper)
        rows.append(row)
    return rows


def fit_score(seed_dir):
    metrics = {r["metric"]: r["value"] for r in efast.read_table(seed_dir / "fit_summary.tsv")}
    for key in ("best_objective", "final_objective", "objective"):
        if key in metrics and math.isfinite(float(metrics[key])):
            return float(metrics[key]), key
    scores = [(float(value), key) for key, value in metrics.items()
              if key.startswith("optimizer_") and key.endswith("_objective") and
              value not in ("NA", "NaN", "Inf", "") and math.isfinite(float(value))]
    if not scores:
        raise ValueError("No finite fit score: " + str(seed_dir))
    score, key = min(scores)
    return score, key


def audit(args):
    root, fit_root = Path(args.out_dir), Path(args.fit_root)
    root.mkdir(parents=True, exist_ok=True)
    table_path, ranges = efast.get_ranges(fit_root)
    _, grid_path = efast.oxygen_grid(args.figure4_dir, True, "")
    r_audit = root / "operator_input_audit.json"
    subprocess.run(["Rscript", "--vanilla", str(HERE / "audit_neighborhood_inputs.R"),
                    str(fit_root), args.figure4_dir, str(r_audit)], check=True)
    operator_audit = json.loads(r_audit.read_text())
    if operator_audit["n_seeds"] != 500 or any(x["status"] != "ok" for x in operator_audit["seeds"]):
        raise ValueError("Operator audit did not cover 500 valid endpoints")
    manifests, bounds, centers = [], [], []
    for number in range(1, 501):
        seed = "seed%d" % number
        directory = fit_root / seed
        best_rows = efast.read_table(directory / "best_params.tsv")
        values = {r["parameter"]: float(r["value"]) for r in best_rows}
        if len(values) != len(best_rows):
            raise ValueError("Duplicate best parameter: " + seed)
        _, per_seed_ranges = efast.get_ranges(directory)
        if per_seed_ranges != ranges:
            raise ValueError("Fit parameter ranges differ: " + seed)
        local = neighborhood_bounds(ranges, values)
        target = root / "bounds" / (seed + ".tsv")
        # Do not rewrite an audited design's bounds if provenance has changed.
        if target.exists():
            previous = efast.read_table(target)
            if len(previous) != len(local):
                raise ValueError("Existing neighborhood has wrong row count: " + seed)
            for old, new in zip(previous, local):
                if any(str(old[k]) != str(v) for k, v in new.items()):
                    raise ValueError("Existing neighborhood changed: " + seed)
        else:
            efast.write_table(target, list(local[0]), local)
        score, score_field = fit_score(directory)
        boundary_count = sum(r["lower_clipped"] == "TRUE" or r["upper_clipped"] == "TRUE" for r in local)
        manifests.append(dict(fit_seed=seed, fit_seed_number=number, fit_score=score,
                              fit_score_field=score_field, n_parameters=len(efast.ACTIVE),
                              n_clipped_parameters=boundary_count,
                              best_params_sha256=efast.sha256(directory / "best_params.tsv"),
                              fit_config_sha256=efast.sha256(directory / "fit_config.rds"),
                              parameter_table_sha256=efast.sha256(directory / "parameter_table.csv"),
                              bounds_sha256=efast.sha256(target), operator_reference_status="ok"))
        bounds.extend(dict(fit_seed=seed, **r) for r in local)
        centers.append([(math.log10(values[r["parameter"]]) if r["transform"] == "log10" else
                         values[r["parameter"]]) for r in ranges])
    # Start with best, worst, and most clipped; then cover normalized parameter space.
    selected = []
    ranked = sorted(range(500), key=lambda i: (manifests[i]["fit_score"], i))
    for i in (ranked[0], ranked[-1], max(range(500), key=lambda i: manifests[i]["n_clipped_parameters"]), 24):
        if i not in selected:
            selected.append(i)
    points = np.asarray(centers)
    points = (points - np.array([r["encoded_lower"] for r in ranges])) / np.array(
        [r["encoded_upper"] - r["encoded_lower"] for r in ranges])
    while len(selected) < 10:
        distances = np.min(np.sum((points[:, None] - points[selected][None, :]) ** 2, axis=2), axis=1)
        distances[selected] = -1
        selected.append(int(np.argmax(distances)))
    for i, row in enumerate(manifests):
        row["pilot_order"] = selected.index(i) + 1 if i in selected else ""
    efast.write_table(root / "seed_manifest.tsv", list(manifests[0]), manifests)
    efast.write_table(root / "neighborhood_bounds.tsv", list(bounds[0]), bounds)
    efast.write_table(root / "parameter_ranges.tsv", list(ranges[0]), ranges)
    for name in SOURCES:
        shutil.copyfile(root.parent / name, root / name)
    provenance = dict(n_fit_seeds=500, n_parameters=14, n_phase_repeats=REPEATS,
                      definition="best +/- 0.10 * natural global span, clipped, inherited transforms",
                      parameter_table_sha256=efast.sha256(table_path),
                      reference_grid_sha256=efast.sha256(grid_path),
                      seed_manifest_sha256=efast.sha256(root / "seed_manifest.tsv"),
                      neighborhood_bounds_sha256=efast.sha256(root / "neighborhood_bounds.tsv"),
                      operator_audit_sha256=efast.sha256(r_audit),
                      pilot_seeds=[manifests[i]["fit_seed"] for i in selected])
    atomic_json(root / "input_manifest.json", provenance)
    print(json.dumps(provenance, indent=2), flush=True)


def verify_completion(run_dir):
    output = run_dir / "outputs.tsv.gz"
    receipt_path = run_dir / "outputs.tsv.gz.receipt.json"
    if not output.exists():
        return False
    if not receipt_path.exists():
        raise ValueError("Output lacks completion receipt: " + str(run_dir))
    receipt = json.loads(receipt_path.read_text())
    metadata = json.loads((run_dir / "metadata.json").read_text())
    expected = {"metadata_sha256": efast.sha256(run_dir / "metadata.json"),
                "samples_sha256": metadata["samples_sha256"],
                "outputs_sha256": efast.sha256(output),
                "evaluator_sha256": efast.sha256(HERE / "evaluate_fixed_o2.R"),
                "n_rows": metadata["n_samples"] * len(metadata["oxygen_pct"])}
    if any(receipt.get(k) != v for k, v in expected.items()):
        raise ValueError("Completion receipt mismatch: " + str(run_dir))
    return True


def analyze_design(run_dir):
    """SALib point indices without its unused within-trajectory bootstrap CI.

    Preserve every output row. Undefined constant trajectories remain NaN with
    their variance recorded; they are never silently replaced by zero indices.
    """
    from SALib.analyze.fast import compute_orders
    meta = json.loads((run_dir / "metadata.json").read_text())
    output = run_dir / "outputs.tsv.gz"
    signature = dict(metadata_sha256=efast.sha256(run_dir / "metadata.json"),
                     outputs_sha256=efast.sha256(output), salib_version=efast.salib_version(),
                     analysis_version=1, constant_relative_variance=1e-20)
    cache, receipt = run_dir / "indices.npz", run_dir / "indices.receipt.json"
    if cache.exists() and receipt.exists():
        previous = json.loads(receipt.read_text())
        if all(previous.get(k) == v for k, v in signature.items()) and previous["cache_sha256"] == efast.sha256(cache):
            with np.load(cache) as saved:
                return {k: saved[k] for k in saved.files}
    grid = np.asarray(meta["oxygen_pct"])
    lookup = {round(float(o), 10): i for i, o in enumerate(grid)}
    data = np.full((meta["n_samples"], len(grid), 2), np.nan)
    seen = np.zeros(data.shape[:2], dtype=bool)
    gap_counts = np.zeros((3, 4), dtype=np.int64)
    gap_min = np.full(3, np.inf)
    with gzip.open(output, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            i = int(row["sample_id"]) - 1
            o = float(row["O2_pct"])
            j = lookup.get(round(o, 10))
            if j is None or not 0 <= i < len(data) or seen[i, j]:
                raise ValueError("Duplicate/unexpected sample or oxygen: " + str(output))
            if row["status"] != "ok" or row["eigenvector_nonnegative"] != "TRUE":
                raise ValueError("Invalid model evaluation: " + str(output))
            data[i, j] = [float(row[k]) for k in efast.OUTPUTS]
            seen[i, j] = True
            gap = float(row["spectral_gap"])
            if not math.isfinite(gap) or gap < 0:
                raise ValueError("Invalid spectral gap")
            b = 0 if o <= 1 else (2 if o >= 3 else 1)
            gap_counts[b] += [1, gap < 1e-6, gap < 1e-4, gap < 1e-3]
            gap_min[b] = min(gap_min[b], gap)
    if not seen.all() or not np.isfinite(data).all():
        raise ValueError("Incomplete/nonfinite FAST trajectories: " + str(output))
    n, m = meta["N"], meta["M"]
    indices = np.full((2, 14, len(grid), 2), np.nan)
    variance = np.zeros((14, len(grid), 2))
    for p in range(14):
        trajectory = data[p * n:(p + 1) * n]
        variance[p] = np.var(trajectory, axis=0)
        for j in range(len(grid)):
            for k in range(2):
                y = trajectory[:, j, k]
                if variance[p, j, k] > 1e-20 * max(1.0, float(np.mean(y)) ** 2):
                    indices[:, p, j, k] = compute_orders(y, n, m, (n - 1) // (2 * m))
    results = dict(indices=indices, variance=variance, oxygen=grid,
                   gap_counts=gap_counts, gap_min=gap_min)
    temporary = run_dir / "indices.incomplete.npz"
    np.savez_compressed(temporary, **results)
    temporary.replace(cache)
    signature["cache_sha256"] = efast.sha256(cache)
    atomic_json(receipt, signature)
    return results


def seed_summary(seed_root, resolutions):
    rows, qc = [], []
    for n in resolutions:
        phase = []
        for r in range(1, REPEATS + 1):
            design = seed_root / "runs" / ("N%d_R%d" % (n, r))
            analyzed = analyze_design(design)
            phase.append(analyzed["indices"])
            for band, counts, minimum in zip(("low_0_to_1pct", "middle_1_to_3pct", "high_3_to_5pct"),
                                            analyzed["gap_counts"], analyzed["gap_min"]):
                if counts[0]:
                    qc.append(dict(N=n, replicate=r, oxygen_band=band, n_evaluations=counts[0],
                                   fraction_gap_lt_1e_6=counts[1] / counts[0],
                                   fraction_gap_lt_1e_4=counts[2] / counts[0],
                                   fraction_gap_lt_1e_3=counts[3] / counts[0], minimum_gap=minimum))
        phase = np.asarray(phase)
        temporary = seed_root / ("phase_indices_N%d.incomplete.npz" % n)
        np.savez_compressed(temporary, indices=phase, oxygen=analyzed["oxygen"])
        temporary.replace(seed_root / ("phase_indices_N%d.npz" % n))
        for k, output in enumerate(efast.OUTPUTS):
            for p, parameter in enumerate(efast.ACTIVE):
                for j, oxygen in enumerate(analyzed["oxygen"]):
                    row = dict(N=n, output=output, parameter=parameter, O2_pct=oxygen,
                               group=efast.GROUP[parameter], replicates=REPEATS)
                    for t, index in enumerate(INDICES):
                        values = phase[:, t, p, j, k]
                        valid = values[np.isfinite(values)]
                        row[index + "_n_valid_repeats"] = len(valid)
                        # Require all five phases for a comparable per-fit-seed index.
                        row[index + "_mean"] = float(valid.mean()) if len(valid) == REPEATS else float("nan")
                        row[index + "_sd"] = float(valid.std(ddof=1)) if len(valid) == REPEATS else float("nan")
                        row[index + "_range"] = float(np.ptp(valid)) if len(valid) == REPEATS else float("nan")
                    rows.append(row)
    maximum = max(resolutions)
    lookup = {(r["output"], r["parameter"], r["O2_pct"]): r for r in rows if r["N"] == maximum}
    for row in rows:
        for index in INDICES:
            row[index + "_delta_vs_max_N"] = abs(row[index + "_mean"] - lookup[
                row["output"], row["parameter"], row["O2_pct"]][index + "_mean"])
    efast.write_table(seed_root / "convergence.tsv.gz", list(rows[0]), rows)
    efast.write_table(seed_root / "spectral_gap_qc.tsv", list(qc[0]), qc)
    return rows


def state(root, **kwargs):
    kwargs.update(updated_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), pid=os.getpid())
    atomic_json(root / "status.json", kwargs)
    print(json.dumps(kwargs, sort_keys=True), flush=True)


def run(args):
    root = Path(args.out_dir)
    manifest = efast.read_table(root / "seed_manifest.tsv")
    inputs = json.loads((root / "input_manifest.json").read_text())
    if efast.sha256(root / "seed_manifest.tsv") != inputs["seed_manifest_sha256"]:
        raise ValueError("Seed manifest changed after audit")
    if args.stage == "full":
        gate = json.loads((root / "pilot" / "convergence_gate.json").read_text())
        if gate["status"] != "passed" or gate["input_manifest_sha256"] != efast.sha256(root / "input_manifest.json"):
            raise ValueError("Full run requires a passing pilot convergence gate for these inputs")
        seeds, resolutions, target = manifest, [513], root
    else:
        seeds = sorted([x for x in manifest if x["pilot_order"]], key=lambda x: int(x["pilot_order"]))
        seeds = seeds[:3 if args.stage == "smoke" else args.pilot_seeds]
        resolutions = [129] if args.stage == "smoke" else [129, 257, 513]
        target = root / "pilot" / args.stage
    total = len(seeds) * len(resolutions) * REPEATS
    done = 0
    for endpoint in seeds:
        seed = endpoint["fit_seed"]
        fit_dir = Path(args.fit_root) / seed
        for name in ("best_params", "fit_config"):
            path = fit_dir / (name + (".tsv" if name == "best_params" else ".rds"))
            if efast.sha256(path) != endpoint[name + "_sha256"]:
                raise ValueError("Audited endpoint changed: " + str(path))
        seed_root = target / "runs" / seed
        seed_root.mkdir(parents=True, exist_ok=True)
        bounds = root / "bounds" / (seed + ".tsv")
        if efast.sha256(bounds) != endpoint["bounds_sha256"]:
            raise ValueError("Audited bounds changed: " + seed)
        for n in resolutions:
            efast.prepare(SimpleNamespace(fit_root=args.fit_root, figure4_dir=args.figure4_dir,
                out_dir=str(seed_root), n=n, m=4, replicates=REPEATS,
                seed=202610080 + int(endpoint["fit_seed_number"]) * 10,
                oxygen="0,0.5,2.5,5", full_grid=args.stage != "smoke",
                fit_seed=seed, ranges_file=str(bounds)))
            for r in range(1, REPEATS + 1):
                design = seed_root / "runs" / ("N%d_R%d" % (n, r))
                state(root, stage=args.stage, status="running", fit_seed=seed, N=n, replicate=r,
                      completed_designs=done, total_designs=total, workers=args.workers)
                if not verify_completion(design):
                    subprocess.run(["Rscript", "--vanilla", str(HERE / "evaluate_fixed_o2.R"),
                        "--samples=" + str(design / "samples.tsv.gz"),
                        "--metadata=" + str(design / "metadata.json"), "--fit_root=" + args.fit_root,
                        "--figure4_dir=" + args.figure4_dir, "--out=" + str(design / "outputs.tsv.gz"),
                        "--workers=" + str(args.workers), "--validate=TRUE"], check=True)
                    if not verify_completion(design):
                        raise ValueError("Missing evaluation completion receipt")
                analyze_design(design)
                done += 1
        seed_summary(seed_root, resolutions)
    if args.stage == "smoke":
        timing_report(root, seeds, target)
    elif args.stage == "convergence":
        convergence_gate(root, seeds, target)
    else:
        summarize(root, manifest)
    state(root, stage=args.stage, status="complete", completed_designs=done, total_designs=total)


def timing_report(root, seeds, target):
    receipts = [json.loads(p.read_text()) for p in target.glob("runs/*/runs/*/outputs.tsv.gz.receipt.json")]
    seconds_per_eval = np.array([r["elapsed_seconds"] / r["n_rows"] for r in receipts])
    total_evaluations = 500 * 513 * 14 * 201 * REPEATS
    report = dict(n_pilot_fit_seeds=len(seeds), n_phase_repeats=REPEATS,
        n_pilot_designs=len(receipts), seconds_per_operator_median=float(np.median(seconds_per_eval)),
        seconds_per_operator_p90=float(np.percentile(seconds_per_eval, 90)),
        full_operator_evaluations=total_evaluations,
        projected_full_wall_days_median=float(np.median(seconds_per_eval) * total_evaluations / 86400),
        projected_full_wall_days_p90=float(np.percentile(seconds_per_eval, 90) * total_evaluations / 86400),
        projected_raw_outputs_GiB=90.7, status="timing_only_not_a_convergence_pass")
    atomic_json(root / "pilot" / "runtime_projection.json", report)
    print(json.dumps(report, indent=2), flush=True)


def convergence_gate(root, seeds, target):
    rows = []
    for endpoint in seeds:
        seed = endpoint["fit_seed"]
        conv = efast.read_table(target / "runs" / seed / "convergence.tsv.gz")
        for output in efast.OUTPUTS:
            for index in INDICES:
                high = [r for r in conv if int(r["N"]) == 513 and r["output"] == output]
                mid = [r for r in conv if int(r["N"]) == 257 and r["output"] == output]
                ranges = np.array([float(r[index + "_range"]) for r in high])
                deltas = np.array([float(r[index + "_delta_vs_max_N"]) for r in mid])
                valid = np.isfinite(ranges) & np.isfinite(deltas)
                # Practical, predeclared tolerances in absolute variance-fraction units.
                p90_range = float(np.percentile(ranges[valid], 90)) if valid.any() else 1.0
                p90_delta = float(np.percentile(deltas[valid], 90)) if valid.any() else 1.0
                rows.append(dict(fit_seed=seed, output=output, index=index,
                    n_cells=len(valid), n_valid_cells=int(valid.sum()),
                    repeat_range_p90=p90_range, resolution_delta_p90=p90_delta,
                    repeat_range_tolerance=.10, resolution_delta_tolerance=.05,
                    passed=str(valid.all() and p90_range <= .10 and p90_delta <= .05).upper()))
    efast.write_table(root / "pilot" / "convergence_diagnostics.tsv", list(rows[0]), rows)
    passed = len(seeds) >= 3 and all(r["passed"] == "TRUE" for r in rows)
    atomic_json(root / "pilot" / "convergence_gate.json", dict(status="passed" if passed else "needs_review",
        n_fit_seeds=len(seeds), n_phase_repeats=REPEATS, N=[129, 257, 513],
        input_manifest_sha256=efast.sha256(root / "input_manifest.json"),
        diagnostics_sha256=efast.sha256(root / "pilot" / "convergence_diagnostics.tsv"),
        explanation="Require each pilot seed/output/index p90 repeat range <=0.10 and N257-to-513 delta <=0.05; all cells defined. Practical thresholds, not a theorem."))


def summarize(root, manifest):
    """Average phases within each endpoint; aggregate endpoints with equal weight."""
    summary_dir = root / "summaries"
    summary_dir.mkdir(parents=True, exist_ok=True)
    phase_arrays = []
    for endpoint in manifest:
        seed_root = root / "runs" / endpoint["fit_seed"]
        with np.load(seed_root / "phase_indices_N513.npz") as data:
            phase_arrays.append(data["indices"])
            oxygen = data["oxygen"]
    # 500 x 5 x 2 indices x 14 parameters x 201 O2 x 2 outputs; no raw-output pooling.
    phases = np.asarray(phase_arrays)
    means = np.mean(phases, axis=1)
    def mean_rows(endpoints):
        for endpoint in endpoints:
            table = root / "runs" / endpoint["fit_seed"] / "convergence.tsv.gz"
            for row in efast.read_table(table):
                if int(row["N"]) == 513:
                    yield dict(fit_seed=endpoint["fit_seed"], **row)
    for start in range(0, len(manifest), 100):
        stop = min(start + 100, len(manifest))
        suffix = "%03d_%03d" % (start + 1, stop)
        np.savez_compressed(summary_dir / ("phase_indices_" + suffix + ".npz"),
            indices=phases[start:stop], oxygen=oxygen,
            fit_seeds=np.array([x["fit_seed"] for x in manifest[start:stop]]),
            parameters=np.array(efast.ACTIVE), outputs=np.array(efast.OUTPUTS), index_names=np.array(INDICES))
        fields = list(next(mean_rows(manifest[start:stop])))
        efast.write_table(summary_dir / ("seed_mean_indices_" + suffix + ".tsv.gz"), fields,
                          mean_rows(manifest[start:stop]))
    rows = []
    for k, output in enumerate(efast.OUTPUTS):
        for p, parameter in enumerate(efast.ACTIVE):
            for j, o2 in enumerate(oxygen):
                row = dict(N=513, output=output, parameter=parameter, O2_pct=o2,
                           n_fit_seeds=len(manifest), replicates=REPEATS, group=efast.GROUP[parameter])
                for t, index in enumerate(INDICES):
                    values = means[:, t, p, j, k]
                    valid = values[np.isfinite(values)]
                    row[index + "_n_valid_seeds"] = len(valid)
                    for suffix, value in (("mean", np.mean(valid)), ("median", np.median(valid)),
                                          ("q25", np.percentile(valid, 25) if len(valid) else np.nan),
                                          ("q75", np.percentile(valid, 75) if len(valid) else np.nan)):
                        row[index + "_" + suffix] = float(value)
                rows.append(row)
    efast.write_table(root / "convergence.tsv", list(rows[0]), rows)
    classification(root, means, oxygen)
    mechanism_summary(root, means, oxygen, manifest)
    for source in SOURCES:
        if not (root / source).exists():
            shutil.copyfile(root.parent / source, root / source)
    efast.plot(SimpleNamespace(out_dir=str(root), figure4_dir=str(root), figure4_layout_dir=None,
        combined_only=True, classification_file=str(root / "efast_o2_sensitivity_classification.tsv"),
        statistic="median", summary_label="median across %d fit seeds of five-phase means" % len(manifest)))
    inventory = []
    for receipt_path in sorted((root / "runs").glob("*/runs/*/outputs.tsv.gz.receipt.json")):
        receipt = json.loads(receipt_path.read_text())
        inventory.append(dict(path=str(receipt_path.parent.relative_to(root)),
            outputs_sha256=receipt["outputs_sha256"], samples_sha256=receipt["samples_sha256"],
            metadata_sha256=receipt["metadata_sha256"], n_rows=receipt["n_rows"],
            bytes=(receipt_path.parent / "outputs.tsv.gz").stat().st_size))
    efast.write_table(root / "raw_output_inventory.tsv", list(inventory[0]), inventory)
    files = [root / "convergence.tsv", root / "raw_output_inventory.tsv",
             root / "efast_o2_sensitivity_classification.tsv", root / "figures" / "efast_four_panel.pdf"]
    files.extend(sorted(summary_dir.glob("*")))
    if any(p.stat().st_size >= 95 * 1024 ** 2 for p in files):
        raise ValueError("A source artifact exceeds the Git size budget; split it before collection")
    efast.write_table(root / "collection_manifest.tsv", ["path", "sha256", "bytes"],
                      [dict(path=str(p.relative_to(root)), sha256=efast.sha256(p), bytes=p.stat().st_size) for p in files])


def classification(root, means, oxygen):
    low = efast.normalized_trapezoid_weights(oxygen, 0, 1)
    high = efast.normalized_trapezoid_weights(oxygen, 3, 5)
    function_rows = efast.read_table(root / "figure4_parameter_groups_source.tsv")
    order = {r["parameter"]: int(r["parameter_order"]) for r in function_rows}
    rng = np.random.default_rng(5826)
    bootstrap = rng.integers(0, len(means), size=(5000, len(means)))
    rows = []
    for t, index in enumerate(INDICES):
        provisional = []
        for p, parameter in enumerate(efast.ACTIVE):
            curves = means[:, t, p, :, 0]
            valid = np.isfinite(curves).all(axis=1)
            observed = np.nanmedian(curves, axis=0)
            # Resample complete fitted-seed curves, after retaining each seed's phase mean.
            # Use the same median statistic as the displayed heatmap in each bootstrap.
            boot_delta = np.empty(5000)
            for start in range(0, 5000, 50):
                sampled = np.nanmedian(curves[bootstrap[start:start + 50]], axis=1)
                boot_delta[start:start + len(sampled)] = sampled @ (low - high)
            p_value = min(1.0, 2 * min((np.count_nonzero(boot_delta <= 0) + 1) / 5001,
                                       (np.count_nonzero(boot_delta >= 0) + 1) / 5001)) if valid.all() else 1.0
            provisional.append(dict(index=index, parameter=parameter, low_o2_score=float(observed @ low),
                high_o2_score=float(observed @ high), low_minus_high=float(observed @ (low - high)),
                global_peak=float(np.nanmax(observed)), bootstrap_sign_p_value=p_value,
                n_fit_seeds=len(curves), n_complete_fit_seed_curves=int(valid.sum()),
                bootstrap_unit="complete fit-seed curve of five-phase mean indices",
                bootstrap_reps=5000, bootstrap_seed=5826, classification_output="dominant_mean_ploidy",
                parameter_order=order[parameter]))
        adjusted = efast.benjamini_hochberg([r["bootstrap_sign_p_value"] for r in provisional])
        groups = ("High O2", "Low O2", "O2-independent")
        for row, q in zip(provisional, adjusted):
            passes = row["global_peak"] > .3 and q < .05
            group = ("Low O2" if row["low_minus_high"] > 0 else "High O2") if passes else "O2-independent"
            row.update(o2_sensitivity_group=group, o2_sensitivity_group_order=groups.index(group) + 1,
                       bh_adjusted_p_value=float(q), minimum_global_peak=.3,
                       decision_rule="global peak >0.3 and BH q<0.05; undefined complete curves force independent")
        provisional.sort(key=lambda r: (r["o2_sensitivity_group_order"], -r["global_peak"], r["parameter_order"]))
        for rank, row in enumerate(provisional, 1):
            row["display_order"] = rank
        rows.extend(provisional)
    efast.write_table(root / "efast_o2_sensitivity_classification.tsv", list(rows[0]), rows)


def mechanism_summary(root, means, oxygen, manifest):
    from scipy.stats import rankdata
    weights = {"low_0_to_1pct": efast.normalized_trapezoid_weights(oxygen, 0, 1),
               "high_3_to_5pct": efast.normalized_trapezoid_weights(oxygen, 3, 5)}
    group_names = sorted(set(efast.GROUP.values()))
    seed_rows, support_rows, ranks, group_rows = [], [], [], []
    for t, index in enumerate(INDICES):
        for k, output in enumerate(efast.OUTPUTS):
            band_means = {band: np.einsum("spj,j->sp", means[:, t, :, :, k], w) for band, w in weights.items()}
            grouped = {band: np.stack([np.mean(values[:, [p for p, name in enumerate(efast.ACTIVE)
                        if efast.GROUP[name] == group]], axis=1) for group in group_names], axis=1)
                       for band, values in band_means.items()}
            for band, group_values in grouped.items():
                for g, name in enumerate(group_names):
                    values = group_values[:, g]
                    valid_values = values[np.isfinite(values)]
                    group_rows.append(dict(index=index, output=output, oxygen_band=band,
                        group=name, n_fit_seeds=len(manifest), n_valid_seeds=len(valid_values),
                        median=float(np.median(valid_values)),
                        q25=float(np.percentile(valid_values, 25)) if len(valid_values) else np.nan,
                        q75=float(np.percentile(valid_values, 75)) if len(valid_values) else np.nan))
            for band, values in band_means.items():
                valid = np.isfinite(values).all(axis=1)
                median_rank = rankdata(-np.median(values[valid], axis=0)) if valid.any() else np.full(14, np.nan)
                for s, endpoint in enumerate(manifest):
                    seed_ranks = rankdata(-values[s]) if valid[s] else np.full(14, np.nan)
                    rho = float(np.corrcoef(seed_ranks, median_rank)[0, 1]) if valid[s] else np.nan
                    ranks.append(dict(fit_seed=endpoint["fit_seed"], index=index, output=output,
                                      oxygen_band=band, rank_spearman_vs_median_band_index=rho))
                    for p, name in enumerate(efast.ACTIVE):
                        seed_rows.append(dict(fit_seed=endpoint["fit_seed"], index=index, output=output,
                            oxygen_band=band, parameter=name, group=efast.GROUP[name],
                            normalized_auc=values[s, p], within_seed_rank=seed_ranks[p]))
            lo, hi = grouped["low_0_to_1pct"], grouped["high_3_to_5pct"]
            death, mis, buf = [group_names.index(g) for g in ("death", "missegregation", "buffering")]
            valid = np.isfinite(lo).all(axis=1) & np.isfinite(hi).all(axis=1)
            conditions = {
                "death_largest_group_mean_at_low_o2": lo[:, death] > np.max(np.delete(lo, death, axis=1), axis=1),
                "missegregation_increases_at_high_o2": hi[:, mis] > lo[:, mis],
                "buffering_increases_at_high_o2": hi[:, buf] > lo[:, buf],
            }
            conditions["all_three_conditions"] = np.logical_and.reduce(list(conditions.values()))
            for name, condition in conditions.items():
                support_rows.append(dict(index=index, output=output, criterion=name,
                    n_fit_seeds=len(manifest), n_valid_seeds=int(valid.sum()),
                    n_supporting_seeds=int(np.count_nonzero(condition & valid)),
                    fraction_supporting_valid_seeds=float(np.mean(condition[valid])) if valid.any() else np.nan,
                    group_statistic="mean of per-parameter normalized AUC; ST overlaps, not group variance"))
    efast.write_table(root / "summaries" / "seed_parameter_bands.tsv.gz", list(seed_rows[0]), seed_rows)
    efast.write_table(root / "summaries" / "rank_consistency.tsv", list(ranks[0]), ranks)
    efast.write_table(root / "summaries" / "mechanism_support.tsv", list(support_rows[0]), support_rows)
    efast.write_table(root / "summaries" / "mechanism_band_summary.tsv", list(group_rows[0]), group_rows)
    lines = ["# In-vivo neighborhood eFAST", "",
        "Each of 500 optimizer endpoints defines an independent rectangular neighborhood: best value +/- 10% of the full natural parameter span, clipped to the original bounds. The inherited log10 or identity sampling transform is retained.", "",
        "Indices are computed separately for every neighborhood and phase. Five phase indices are averaged within each fit seed; the heatmaps show medians across fit seeds. Source tables also report means, IQR and valid-seed counts. Raw outputs from different neighborhoods are never pooled into one FAST analysis.", "",
        "The bootstrap resamples complete fitted-seed sensitivity curves with their five-phase means retained. Optimizer endpoints are not independent posterior draws; bootstrap classifications describe repeatability across these endpoints and are not posterior significance claims.", "",
        "Mechanism comparisons use group means of parameter indices. Total effects overlap through interactions; these means are not additive group variance fractions. The original Figure 4B correlation annotation supplies directional context.", "",
        "## Fraction of evaluable fit seeds supporting the proposed shift", ""]
    for row in support_rows:
        if row["criterion"] == "all_three_conditions":
            lines.append("- %s, %s: %.3f (%d/%d); requires low-O2 death dominance and higher-O2 increases in both missegregation and buffering." %
                         (row["output"], row["index"], row["fraction_supporting_valid_seeds"],
                          row["n_supporting_seeds"], row["n_valid_seeds"]))
    lines += ["", "Interpret these fractions with the pilot resolution and phase diagnostics. Constant-output trajectories have undefined indices and remain recorded as NaN; reported denominators identify the evaluable subset.", ""]
    (root / "interpretation.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("audit", "run"):
        p = sub.add_parser(name)
        p.add_argument("--fit-root", required=True)
        p.add_argument("--figure4-dir", required=True)
        p.add_argument("--out-dir", required=True)
        if name == "run":
            p.add_argument("--stage", choices=("smoke", "convergence", "full"), required=True)
            p.add_argument("--workers", type=int, default=16)
            p.add_argument("--pilot-seeds", type=int, choices=range(3, 11), default=10)
    s = sub.add_parser("summarize")
    s.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    if args.command == "audit":
        audit(args)
    elif args.command == "run":
        try:
            run(args)
        except Exception as error:
            state(Path(args.out_dir), stage=args.stage, status="failed", error=str(error))
            raise
    else:
        root = Path(args.out_dir)
        summarize(root, efast.read_table(root / "seed_manifest.tsv"))


if __name__ == "__main__":
    main()
