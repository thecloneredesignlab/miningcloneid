#!/usr/bin/env python3
"""Targeted, audited recovery of failed neighborhood trajectories and summaries."""
import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns

HERE = Path(__file__).resolve().parent
SOURCES = ("recover_eigen.py", "repair_eigen_points.R", "perron_high_precision.cpp")


def compressed_ids(values):
    values = sorted(set(values))
    groups = []
    if not values:
        return ""
    start = end = values[0]
    for value in values[1:]:
        if value == end + 1:
            end = value
        else:
            groups.append(str(start) if start == end else "%d-%d" % (start, end))
            start = end = value
    groups.append(str(start) if start == end else "%d-%d" % (start, end))
    return ",".join(groups)


def prepare(root):
    ns.context(root)
    validation = json.loads((root / "slurm/perron_validation.json").read_text())
    if validation["status"] != "passed" or validation["n_matrices"] < 10:
        raise ValueError("Require completed 50/100-digit validation")
    pending, previous = [], []
    for task_id in range(1, 35001):
        path = root / "slurm/task_receipts" / ("task_%05d.json" % task_id)
        receipt = json.loads(path.read_text()) if path.exists() else {}
        if receipt.get("status") != "complete":
            pending.append(task_id)
            previous.append(dict(ns.task_identity(task_id), previous_status=receipt.get("status", "missing"),
                previous_error=receipt.get("error", "")))
    if not pending:
        raise ValueError("No pending recovery tasks")
    directory = root / "slurm/eigen_recovery"
    directory.mkdir(exist_ok=True)
    efast.write_table(root / "slurm/numerical_recovery_tasks.tsv", list(previous[0]), previous)
    summary_ids = [s for s in range(1, 501) if not (root / "runs" / ("seed%d" % s) / "summary.receipt.json").exists()]
    plan = dict(recorded_utc=ns.now(), task_ids=pending, array_spec=compressed_ids(pending),
        summary_seed_ids=summary_ids, summary_array_spec=compressed_ids(summary_ids),
        n_tasks=len(pending), n_preserved_tasks=35000-len(pending),
        n_summaries=len(summary_ids), validation_sha256=efast.sha256(root / "slurm/perron_validation.json"),
        code_sha256={name: efast.sha256(HERE / name) for name in SOURCES},
        wrapper_sha256=efast.sha256(HERE / "run_eigen_recovery.sbatch"),
        perron_template_sha256=efast.sha256(root / "slurm/perron_template.tar.gz"),
        submission_plan_sha256=efast.sha256(root / "slurm/submission_plan.json"),
        qos="xxlarge", time="12:00:00", specified_compute_node=None, array_concurrency_limit=None,
        sampling_changed=False, convergence_policy="diagnostic_only")
    nh.atomic_json(root / "slurm/numerical_recovery_plan.json", plan)
    print(json.dumps(plan, indent=2), flush=True)


def context(root):
    plan = json.loads((root / "slurm/numerical_recovery_plan.json").read_text())
    for name, expected in plan["code_sha256"].items():
        if efast.sha256(HERE / name) != expected:
            raise ValueError("Submitted recovery code changed: " + name)
    if efast.sha256(HERE / "run_eigen_recovery.sbatch") != plan["wrapper_sha256"]:
        raise ValueError("Recovery wrapper changed")
    if efast.sha256(root / "slurm/submission_plan.json") != plan["submission_plan_sha256"]:
        raise ValueError("Recovery parent plan changed")
    ns.context(root)
    return plan


def invalid(row):
    return (row["status"] != "ok" or row["eigenvector_nonnegative"] != "TRUE" or
        any(row[name] in ("NA", "NaN", "Inf", "-Inf") or not math.isfinite(float(row[name])) for name in efast.OUTPUTS))


def raw_indices(rows, metadata):
    """Numerical impact comparison only; these unchecked old indices are never results."""
    from SALib.analyze.fast import compute_orders
    data = np.asarray([[float(row[name]) for name in efast.OUTPUTS] for row in rows])
    data = data.reshape(metadata["N"], len(metadata["oxygen_pct"]), 2)
    result = np.full((2, len(metadata["oxygen_pct"]), 2), np.nan)
    for j in range(len(metadata["oxygen_pct"])):
        for k in range(2):
            y = data[:, j, k]
            if np.var(y) > 1e-20 * max(1., float(np.mean(y))**2):
                result[:, j, k] = compute_orders(y, metadata["N"], metadata["M"],
                    (metadata["N"]-1)//(2*metadata["M"]))
    return result


def repair(root, task_id):
    plan = context(root)
    if task_id not in plan["task_ids"]:
        raise ValueError("Task is not in the targeted recovery plan")
    task = ns.task_identity(task_id)
    main_plan = ns.context(root)
    target = ns.ensure_seed(root, task["fit_seed"], main_plan)
    parent = target / "runs" / ("N513_R%d" % task["replicate"])
    part = ns.split_trajectory(parent, task["parameter_index"])
    original = part / "outputs.tsv.gz.part01.gz"
    begin = time.monotonic()
    with ns.locked(part / "recovery.lock"):
        if (part / "outputs.tsv.gz").exists() and not (part / "outputs.tsv.gz.receipt.json").exists():
            (part / "outputs.tsv.gz").replace(part / "outputs.recovery.interrupted.tsv.gz")
        if not nh.verify_completion(part) and original.exists():
            with gzip.open(original, "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            metadata = json.loads((part / "metadata.json").read_text())
            if len(rows) != metadata["N"] * len(metadata["oxygen_pct"]):
                raise ValueError("Original failed output is incomplete")
            bad = [row for row in rows if invalid(row)]
            if not bad:
                raise ValueError("Expected invalid points in failed canonical output")
            recovery = part / "recovery"
            recovery.mkdir(exist_ok=True)
            nh.atomic_json(recovery / "failed_points.json", bad)
            subprocess.run(["Rscript", "--vanilla", str(HERE / "repair_eigen_points.R"), str(root), str(part)], check=True)
            corrected = json.loads((recovery / "corrections.json").read_text())
            if corrected["status"] != "passed" or corrected["n_points"] != len(bad):
                raise ValueError("Incomplete high-precision recovery")
            replacements = {(int(r["sample_id"]), round(float(r["O2_pct"]), 10)): r for r in corrected["reports"]}
            wanted = {(int(r["sample_id"]), round(float(r["O2_pct"]), 10)) for r in bad}
            if set(replacements) != wanted or len(replacements) != len(bad):
                raise ValueError("Recovery points do not match original failures")
            before = raw_indices(rows, metadata)
            unchanged = 0
            for row in rows:
                key = (int(row["sample_id"]), round(float(row["O2_pct"]), 10))
                if key in replacements:
                    r = replacements[key]
                    if r["original"] != row:
                        raise ValueError("Recovery original row changed")
                    for name in efast.OUTPUTS:
                        row[name] = repr(float(r[name]))
                    row["eigenvector_nonnegative"] = "TRUE"
                else:
                    unchanged += 1
            temporary = part / "outputs.recovery.incomplete.tsv.gz"
            efast.write_table(temporary, ns.OUTPUT_FIELDS, rows)
            temporary.replace(part / "outputs.tsv.gz")
            after = ns.analyze_trajectory(part)["indices"]
            artifacts = {name: efast.sha256(recovery / name) for name in ("failed_points.json", "corrections.json")}
            artifacts["../outputs.tsv.gz.part01.gz"] = efast.sha256(original)
            for r in corrected["reports"]:
                artifacts[r["matrix_file"]] = r["matrix_sha256"]
            changes = []
            for t, index in enumerate(nh.INDICES):
                for k, output in enumerate(efast.OUTPUTS):
                    defined = np.isfinite(before[t,:,k]) & np.isfinite(after[t,:,k])
                    changes.append(dict(index=index, output=output,
                        max_absolute_index_difference=float(np.max(np.abs(before[t,defined,k]-after[t,defined,k]))) if defined.any() else None,
                        n_compared_oxygen_points=int(defined.sum()),
                        undefined_pattern_changed=bool(np.any(np.isfinite(before[t,:,k])!=np.isfinite(after[t,:,k])))))
            proof = dict(status="passed", task_id=task_id, recovered_utc=ns.now(),
                method="Noda positivity-preserving inverse iteration, 50 and 100 decimal digits",
                n_recomputed_points=len(bad), n_unchanged_good_rows=unchanged,
                original_source_part_sha256=efast.sha256(original),
                outputs_sha256=efast.sha256(part / "outputs.tsv.gz"),
                code_sha256={name: efast.sha256(HERE / name) for name in SOURCES},
                artifact_sha256=artifacts, index_impact=changes,
                no_matrix_regularization=True, no_sampling_change=True,
                spectral_gap_source="unchanged canonical double-precision spectrum diagnostic")
            nh.atomic_json(recovery / "numerical_recovery.json", proof)
            nh.atomic_json(part / "outputs.tsv.gz.receipt.json", dict(fit_seed=metadata["fit_seed"],
                n_rows=len(rows), metadata_sha256=efast.sha256(part / "metadata.json"),
                samples_sha256=metadata["samples_sha256"], outputs_sha256=proof["outputs_sha256"],
                evaluator_sha256=efast.sha256(HERE / "evaluate_fixed_o2.R"),
                evaluation_method="canonical fixed-O2 evaluation with certified high-precision fallback",
                recovery_manifest_sha256=efast.sha256(recovery / "numerical_recovery.json"),
                recovery_evaluator_sha256=efast.sha256(HERE / "repair_eigen_points.R"),
                elapsed_seconds=time.monotonic()-begin, workers=1, completed_at=ns.now()))
    # The existing worker verifies the extended receipt, computes SALib indices,
    # and writes its ordinary successful task marker. Missing container launches
    # have no partial file and therefore run the canonical evaluator normally.
    try:
        ns.worker(root, task_id)
    except subprocess.CalledProcessError:
        if not original.exists():
            raise
        # A formerly missing launch can reveal a numerical failure once evaluated.
        # Recover that new partial file using the same certified method.
        repair(root, task_id)
        return
    print(json.dumps(dict(task_id=task_id,status="complete",elapsed_seconds=time.monotonic()-begin)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("prepare","repair"))
    parser.add_argument("--out-dir",required=True)
    parser.add_argument("--task-id",type=int)
    args=parser.parse_args()
    root=Path(args.out_dir).resolve()
    if args.command=="prepare": prepare(root)
    else: repair(root,args.task_id)


if __name__=="__main__":
    main()
