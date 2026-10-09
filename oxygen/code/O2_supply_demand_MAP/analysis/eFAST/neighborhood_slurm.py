#!/usr/bin/env python3
"""Slurm tasks: one fitted endpoint, phase and complete FAST parameter trajectory."""
import argparse
from contextlib import contextmanager
import csv
import datetime
import fcntl
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import efast
import neighborhood as nh

N = 513
M = 4
SIF_SHA256 = "0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad"
HERE = Path(__file__).resolve().parent
OUTPUT_FIELDS = ("sample_id", "O2_pct", *efast.OUTPUTS, "spectral_gap",
                 "eigenvector_nonnegative", "status")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


@contextmanager
def locked(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def task_rows():
    for seed in range(1, 501):
        for repeat in range(1, 6):
            for p, name in enumerate(efast.ACTIVE):
                yield dict(task_id=((seed - 1) * 5 + repeat - 1) * 14 + p + 1,
                    fit_seed="seed%d" % seed, fit_seed_number=seed, replicate=repeat,
                    parameter_index=p, focal_parameter=name, N=N, M=M,
                    phase_seed=202610080 + seed * 10 + repeat - 1,
                    parent_sample_id_start=p * N + 1, parent_sample_id_end=(p + 1) * N)


def task_identity(task_id):
    if not 1 <= task_id <= 35000:
        raise ValueError("Task ID must be 1..35000")
    seed, offset = divmod(task_id - 1, 70)
    repeat, p = divmod(offset, 14)
    return dict(task_id=task_id, fit_seed="seed%d" % (seed + 1), fit_seed_number=seed + 1,
        replicate=repeat + 1, parameter_index=p, focal_parameter=efast.ACTIVE[p], N=N, M=M,
        phase_seed=202610080 + (seed + 1) * 10 + repeat,
        parent_sample_id_start=p * N + 1, parent_sample_id_end=(p + 1) * N)


def prepare(args):
    root = Path(args.out_dir).resolve()
    manifest = efast.read_table(root / "seed_manifest.tsv")
    inputs = json.loads((root / "input_manifest.json").read_text())
    if len(manifest) != 500 or inputs["n_phase_repeats"] != 5 or inputs["n_parameters"] != 14:
        raise ValueError("Require the audited 500-endpoint, 14-parameter, five-phase inputs")
    for name, field in (("seed_manifest.tsv", "seed_manifest_sha256"),
                        ("operator_input_audit.json", "operator_audit_sha256"),
                        ("neighborhood_bounds.tsv", "neighborhood_bounds_sha256")):
        if efast.sha256(root / name) != inputs[field]:
            raise ValueError("Audit hash mismatch: " + name)
    for endpoint in manifest:
        for name, path in (("best_params", Path(args.fit_root) / endpoint["fit_seed"] / "best_params.tsv"),
                           ("fit_config", Path(args.fit_root) / endpoint["fit_seed"] / "fit_config.rds"),
                           ("bounds", root / "bounds" / (endpoint["fit_seed"] + ".tsv"))):
            if efast.sha256(path) != endpoint[name + "_sha256"]:
                raise ValueError("Audited input changed: " + str(path))
    if shutil.disk_usage(root).free < 200 * 1024 ** 3:
        raise ValueError("Require at least 200 GiB free")
    sif = Path(args.sif).resolve()
    if efast.sha256(sif) != SIF_SHA256:
        raise ValueError("SIF hash mismatch")
    if efast.sha256(Path(args.figure4_dir) / "fixed_o2_dominant_ploidy_201grid.tsv") != inputs["reference_grid_sha256"]:
        raise ValueError("Figure 4 grid/reference changed")
    directory = root / "slurm"
    directory.mkdir(exist_ok=True)
    tasks = list(task_rows())
    efast.write_table(directory / "task_manifest.tsv", list(tasks[0]), tasks)
    template_path = directory / "rcpp_template_manifest.json"
    template = json.loads(template_path.read_text())
    if template["sif_sha256"] != SIF_SHA256 or template["model_cpp_sha256"] != efast.sha256(HERE / "../../model/model_O2_supply_demand_MAP.cpp"):
        raise ValueError("Rcpp template source/image mismatch")
    if efast.sha256(directory / "rcpp_template.tar.gz") != template["archive_sha256"]:
        raise ValueError("Rcpp template archive changed")
    names = ("neighborhood_slurm.py", "neighborhood.py", "efast.py", "evaluate_fixed_o2.R",
             "run_neighborhood_array.sbatch", "../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R",
             "../../model/model_O2_supply_demand_MAP.R", "../../model/model_O2_supply_demand_MAP.cpp")
    plan = dict(n_tasks=35000, n_fit_seeds=500, N=N, M=M, n_phases=5,
        granularity="fit seed x phase x one complete 513-vector FAST trajectory; 201 oxygen points",
        fit_root=str(Path(args.fit_root).resolve()), figure4_dir=str(Path(args.figure4_dir).resolve()),
        result_root=str(root), sif=str(sif), sif_sha256=SIF_SHA256,
        sif_bytes=sif.stat().st_size, sif_mtime_ns=sif.stat().st_mtime_ns,
        input_manifest_sha256=efast.sha256(root / "input_manifest.json"),
        seed_manifest_sha256=efast.sha256(root / "seed_manifest.tsv"),
        task_manifest_sha256=efast.sha256(directory / "task_manifest.tsv"),
        code_sha256={name: efast.sha256(HERE / name) for name in names},
        git_commit=subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "HEAD"], text=True).strip(),
        prepared_utc=now(), qos="xxlarge", time="12:00:00", task_cpus=1, task_mem="4G",
        array_concurrency_limit=None, specified_compute_node=None,
        convergence_policy="diagnostic_only")
    plan["rcpp_template"] = template
    previous = directory / "submission_plan.json"
    if previous.exists():
        old = json.loads(previous.read_text())
        if any(old.get(k) != v for k, v in plan.items() if k != "prepared_utc"):
            if not args.replace_plan:
                raise ValueError("Existing Slurm plan differs; inspect before overwriting")
            archived = directory / ("submission_plan_before_" + old["git_commit"][:12] + ".json")
            if not archived.exists():
                shutil.copyfile(previous, archived)
        else:
            plan = old
    nh.atomic_json(previous, plan)
    nh.atomic_json(directory / "execution_backend.json", dict(backend="slurm", recorded_utc=now(),
        submission_plan_sha256=efast.sha256(previous), explanation="Full arrays may start before the direct pilot finishes; the serial full stage hands over without evaluating."))
    nh.record_full_execution_policy(root)
    print(json.dumps(plan, indent=2), flush=True)


def context(root):
    plan = json.loads((root / "slurm" / "submission_plan.json").read_text())
    for name, expected in plan["code_sha256"].items():
        if efast.sha256(HERE / name) != expected:
            raise ValueError("Submitted code changed: " + name)
    for name, field in (("input_manifest.json", "input_manifest_sha256"),
                        ("seed_manifest.tsv", "seed_manifest_sha256")):
        if efast.sha256(root / name) != plan[field]:
            raise ValueError("Submitted inputs changed: " + name)
    sif = Path(plan["sif"])
    if sif.stat().st_size != plan["sif_bytes"] or sif.stat().st_mtime_ns != plan["sif_mtime_ns"]:
        raise ValueError("Validated immutable SIF changed")
    archive = root / "slurm" / "rcpp_template.tar.gz"
    template = plan["rcpp_template"]
    if archive.stat().st_size != template["archive_bytes"] or archive.stat().st_mtime_ns != template["archive_mtime_ns"]:
        raise ValueError("Validated Rcpp template archive changed")
    return plan


def ensure_seed(root, seed, plan):
    target = root / "runs" / seed
    with locked(target / "prepare.lock"):
        endpoint = next(r for r in efast.read_table(root / "seed_manifest.tsv") if r["fit_seed"] == seed)
        for name, path in (("best_params", Path(plan["fit_root"]) / seed / "best_params.tsv"),
                           ("fit_config", Path(plan["fit_root"]) / seed / "fit_config.rds"),
                           ("bounds", root / "bounds" / (seed + ".tsv"))):
            if efast.sha256(path) != endpoint[name + "_sha256"]:
                raise ValueError("Audited input changed: " + str(path))
        receipt = target / "prepared.json"
        if not receipt.exists():
            # Reuse completed pilot phases, including original compressed bytes/hashes.
            for r in range(1, 6):
                source = root / "pilot" / "convergence" / "runs" / seed / "runs" / ("N513_R%d" % r)
                dest = target / "runs" / source.name
                if source.exists() and (source / "metadata.json").exists() and not dest.exists():
                    metadata = json.loads((source / "metadata.json").read_text())
                    if efast.sha256(source / "samples.tsv.gz") != metadata["samples_sha256"]:
                        raise ValueError("Pilot sample hash mismatch")
                    dest.mkdir(parents=True)
                    for name in ("metadata.json", "samples.tsv.gz"):
                        shutil.copyfile(source / name, dest / name)
            efast.prepare(SimpleNamespace(fit_root=plan["fit_root"], figure4_dir=plan["figure4_dir"],
                out_dir=str(target), n=N, m=M, replicates=5,
                seed=202610080 + int(endpoint["fit_seed_number"]) * 10,
                full_grid=True, oxygen="", fit_seed=seed, ranges_file=str(root / "bounds" / (seed + ".tsv"))))
            for r in range(1, 6):
                source = root / "pilot" / "convergence" / "runs" / seed / "runs" / ("N513_R%d" % r)
                dest = target / "runs" / source.name
                if (source / "outputs.tsv.gz.receipt.json").exists() and not (dest / "outputs.tsv.gz").exists():
                    if nh.verify_completion(source):
                        if efast.sha256(dest / "metadata.json") != efast.sha256(source / "metadata.json"):
                            raise ValueError("Pilot and full metadata differ")
                        for name in ("outputs.tsv.gz", "outputs.tsv.gz.receipt.json", "indices.npz", "indices.receipt.json"):
                            if (source / name).exists():
                                shutil.copyfile(source / name, dest / name)
            nh.atomic_json(receipt, dict(input_manifest_sha256=plan["input_manifest_sha256"], prepared_utc=now()))
        elif json.loads(receipt.read_text())["input_manifest_sha256"] != plan["input_manifest_sha256"]:
            raise ValueError("Seed preparation inputs changed")
    return target


def split_trajectory(parent, p):
    metadata = json.loads((parent / "metadata.json").read_text())
    if metadata["n_samples"] != metadata["N"] * 14 or not 0 <= p < 14:
        raise ValueError("Invalid complete parent FAST design")
    part = parent / "trajectories" / ("P%02d_%s" % (p + 1, efast.ACTIVE[p]))
    part.mkdir(parents=True, exist_ok=True)
    if efast.sha256(parent / "samples.tsv.gz") != metadata["samples_sha256"]:
        raise ValueError("Parent samples changed")
    identity = dict(metadata, n_samples=metadata["N"],
        trajectory_parameter_index=p, trajectory_parameter=efast.ACTIVE[p],
        parent_metadata_sha256=efast.sha256(parent / "metadata.json"),
        parent_samples_sha256=metadata["samples_sha256"], parent_sample_id_offset=p * metadata["N"])
    if (part / "metadata.json").exists():
        old = json.loads((part / "metadata.json").read_text())
        if any(old.get(k) != v for k, v in identity.items() if k != "samples_sha256"):
            raise ValueError("Trajectory metadata changed")
        if efast.sha256(part / "samples.tsv.gz") != old["samples_sha256"]:
            raise ValueError("Trajectory samples changed")
        return part
    samples = efast.read_table(parent / "samples.tsv.gz")
    block = samples[p * metadata["N"]:(p + 1) * metadata["N"]]
    for i, row in enumerate(block, 1):
        if int(row["sample_id"]) != p * metadata["N"] + i:
            raise ValueError("Parent FAST sequence is incomplete")
        row["sample_id"] = i
    efast.write_table(part / "samples.tsv.gz", list(block[0]), block)
    identity["samples_sha256"] = efast.sha256(part / "samples.tsv.gz")
    nh.atomic_json(part / "metadata.json", identity)
    return part


def analyze_trajectory(part):
    from SALib.analyze.fast import compute_orders
    metadata = json.loads((part / "metadata.json").read_text())
    grid = np.asarray(metadata["oxygen_pct"])
    lookup = {round(float(o), 10): j for j, o in enumerate(grid)}
    data = np.full((metadata["N"], len(grid), 2), np.nan)
    seen = np.zeros(data.shape[:2], dtype=bool)
    counts, minimum = np.zeros((3, 4), dtype=np.int64), np.full(3, np.inf)
    with gzip.open(part / "outputs.tsv.gz", "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            i, o = int(row["sample_id"]) - 1, float(row["O2_pct"])
            j = lookup.get(round(o, 10))
            if j is None or not 0 <= i < len(data) or seen[i, j]:
                raise ValueError("Duplicate/unexpected trajectory row")
            if row["status"] != "ok" or row["eigenvector_nonnegative"] != "TRUE":
                raise ValueError("Invalid trajectory evaluation")
            data[i, j] = [float(row[k]) for k in efast.OUTPUTS]
            seen[i, j] = True
            gap = float(row["spectral_gap"])
            if not np.isfinite(gap) or gap < 0:
                raise ValueError("Invalid spectral gap")
            band = 0 if o <= 1 else 2 if o >= 3 else 1
            counts[band] += [1, gap < 1e-6, gap < 1e-4, gap < 1e-3]
            minimum[band] = min(minimum[band], gap)
    if not seen.all() or not np.isfinite(data).all():
        raise ValueError("Incomplete/nonfinite FAST trajectory")
    variance = np.var(data, axis=0)
    indices = np.full((2, len(grid), 2), np.nan)
    for j in range(len(grid)):
        for k in range(2):
            y = data[:, j, k]
            if variance[j, k] > 1e-20 * max(1., float(np.mean(y)) ** 2):
                indices[:, j, k] = compute_orders(y, metadata["N"], metadata["M"],
                    (metadata["N"] - 1) // (2 * metadata["M"]))
    return dict(indices=indices, variance=variance, oxygen=grid, gap_counts=counts, gap_min=minimum)


def worker(root, task_id):
    task = task_identity(task_id)
    marker = root / "slurm" / "task_receipts" / ("task_%05d.json" % task_id)
    try:
        plan = context(root)
        seed_root = ensure_seed(root, task["fit_seed"], plan)
        parent = seed_root / "runs" / ("N513_R%d" % task["replicate"])
        with locked(parent / ("trajectory_%02d.lock" % task["parameter_index"])):
            if (parent / "outputs.tsv.gz.receipt.json").exists() and nh.verify_completion(parent):
                nh.atomic_json(marker, dict(task, status="complete", reused_full_phase=True,
                    parent_receipt_sha256=efast.sha256(parent / "outputs.tsv.gz.receipt.json"),
                    completed_utc=now(), slurm_job_id=os.environ.get("SLURM_JOB_ID"), host=os.uname().nodename))
                return
            part = split_trajectory(parent, task["parameter_index"])
            if not nh.verify_completion(part):
                subprocess.run(["Rscript", "--vanilla", str(HERE / "evaluate_fixed_o2.R"),
                    "--samples=" + str(part / "samples.tsv.gz"), "--metadata=" + str(part / "metadata.json"),
                    "--fit_root=" + plan["fit_root"], "--figure4_dir=" + plan["figure4_dir"],
                    "--out=" + str(part / "outputs.tsv.gz"), "--workers=1", "--validate=TRUE"], check=True)
            if not nh.verify_completion(part):
                raise ValueError("Missing trajectory completion receipt")
            result = analyze_trajectory(part)
            write_trajectory_indices(part, result)
            nh.atomic_json(marker, dict(task, status="complete", reused_full_phase=False,
                trajectory_receipt_sha256=efast.sha256(part / "outputs.tsv.gz.receipt.json"),
                indices_sha256=efast.sha256(part / "indices.npz"), completed_utc=now(),
                slurm_job_id=os.environ.get("SLURM_JOB_ID"), host=os.uname().nodename))
    except Exception as error:
        nh.atomic_json(marker, dict(task, status="failed", error=str(error), recorded_utc=now(),
            slurm_job_id=os.environ.get("SLURM_JOB_ID"), host=os.uname().nodename))
        raise


def write_trajectory_indices(part, result):
    temporary = part / "indices.incomplete.npz"
    np.savez_compressed(temporary, **result)
    temporary.replace(part / "indices.npz")
    nh.atomic_json(part / "indices.receipt.json", dict(
        outputs_sha256=efast.sha256(part / "outputs.tsv.gz"),
        metadata_sha256=efast.sha256(part / "metadata.json"),
        cache_sha256=efast.sha256(part / "indices.npz"), salib_version=efast.salib_version(),
        analyzer_sha256=efast.sha256(Path(__file__))))


def assemble_phase(parent):
    if nh.verify_completion(parent):
        return
    metadata = json.loads((parent / "metadata.json").read_text())
    parts = [parent / "trajectories" / ("P%02d_%s" % (p + 1, name)) for p, name in enumerate(efast.ACTIVE)]
    receipts, results = [], []
    for part in parts:
        if not nh.verify_completion(part):
            raise ValueError("Incomplete trajectory: " + str(part))
        part_metadata = json.loads((part / "metadata.json").read_text())
        if part_metadata["parent_metadata_sha256"] != efast.sha256(parent / "metadata.json"):
            raise ValueError("Trajectory parent changed")
        cache_receipt = json.loads((part / "indices.receipt.json").read_text())
        for name, field in (("outputs.tsv.gz", "outputs_sha256"),
                            ("metadata.json", "metadata_sha256"), ("indices.npz", "cache_sha256")):
            if efast.sha256(part / name) != cache_receipt[field]:
                raise ValueError("Trajectory analysis cache changed: " + str(part))
        if cache_receipt["analyzer_sha256"] != efast.sha256(Path(__file__)):
            raise ValueError("Trajectory analysis code changed")
        receipts.append(json.loads((part / "outputs.tsv.gz.receipt.json").read_text()))
        with np.load(part / "indices.npz") as data:
            results.append({k: data[k] for k in data.files})
    def rows():
        for p, part in enumerate(parts):
            with gzip.open(part / "outputs.tsv.gz", "rt", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    row["sample_id"] = int(row["sample_id"]) + p * metadata["N"]
                    yield row
    temporary = parent / "outputs.incomplete.tsv.gz"
    efast.write_table(temporary, OUTPUT_FIELDS, rows())
    temporary.replace(parent / "outputs.tsv.gz")
    nh.atomic_json(parent / "outputs.tsv.gz.receipt.json", dict(fit_seed=metadata["fit_seed"],
        n_rows=sum(r["n_rows"] for r in receipts), metadata_sha256=efast.sha256(parent / "metadata.json"),
        samples_sha256=metadata["samples_sha256"], outputs_sha256=efast.sha256(parent / "outputs.tsv.gz"),
        evaluator_sha256=efast.sha256(HERE / "evaluate_fixed_o2.R"),
        elapsed_seconds=sum(r["elapsed_seconds"] for r in receipts), workers=1, completed_at=now(),
        assembler_sha256=efast.sha256(Path(__file__)),
        trajectory_receipts_sha256=[efast.sha256(p / "outputs.tsv.gz.receipt.json") for p in parts]))
    combined = dict(indices=np.stack([r["indices"] for r in results], axis=1),
        variance=np.stack([r["variance"] for r in results]), oxygen=results[0]["oxygen"],
        gap_counts=np.sum([r["gap_counts"] for r in results], axis=0),
        gap_min=np.min([r["gap_min"] for r in results], axis=0))
    np.savez_compressed(parent / "indices.npz", **combined)
    nh.atomic_json(parent / "indices.receipt.json", dict(
        metadata_sha256=efast.sha256(parent / "metadata.json"),
        outputs_sha256=efast.sha256(parent / "outputs.tsv.gz"), salib_version=efast.salib_version(),
        analysis_version=1, constant_relative_variance=1e-20, cache_sha256=efast.sha256(parent / "indices.npz")))
    if not nh.verify_completion(parent):
        raise ValueError("Assembled phase is incomplete")


def seed_summary(root, seed_number):
    plan = context(root)
    seed = "seed%d" % seed_number
    if not 1 <= seed_number <= 500:
        raise ValueError("Seed number must be 1..500")
    target = ensure_seed(root, seed, plan)
    with locked(target / "summary.lock"):
        for task_id in range((seed_number - 1) * 70 + 1, seed_number * 70 + 1):
            path = root / "slurm" / "task_receipts" / ("task_%05d.json" % task_id)
            if not path.exists() or json.loads(path.read_text())["status"] != "complete":
                raise ValueError("Seed has an incomplete Slurm task: %d" % task_id)
        for r in range(1, 6):
            assemble_phase(target / "runs" / ("N513_R%d" % r))
        rows = nh.seed_summary(target, [N])
        nh.write_full_diagnostics(root, seed, rows, aggregate=False)
        nh.atomic_json(target / "summary.receipt.json", dict(fit_seed=seed, status="complete",
            completed_utc=now(), phase_indices_sha256=efast.sha256(target / "phase_indices_N513.npz"),
            convergence_sha256=efast.sha256(target / "convergence.tsv.gz")))


def status(root):
    complete, failed, missing = [], [], []
    for task_id in range(1, 35001):
        path = root / "slurm" / "task_receipts" / ("task_%05d.json" % task_id)
        if not path.exists():
            missing.append(task_id)
        elif json.loads(path.read_text())["status"] == "complete":
            complete.append(task_id)
        else:
            failed.append(task_id)
    seeds = [p.parent.name for p in (root / "runs").glob("*/summary.receipt.json")
             if json.loads(p.read_text()).get("status") == "complete"]
    report = dict(checked_utc=now(), total_tasks=35000, completed_tasks=len(complete),
        failed_tasks=len(failed), missing_or_running_tasks=len(missing), completed_seed_summaries=len(seeds),
        missing_task_ids=missing, failed_task_ids=failed, convergence_policy="diagnostic_only")
    nh.atomic_json(root / "slurm" / "status.json", report)
    print(json.dumps({k: v for k, v in report.items() if not k.endswith("_ids")}, indent=2), flush=True)
    return report


def finalize(root):
    context(root)
    with locked(root / "slurm" / "finalize.lock"):
        report = status(root)
        if report["completed_tasks"] != 35000 or report["completed_seed_summaries"] != 500:
            raise ValueError("Full aggregation requires all task receipts and 500 valid seed summaries; see slurm/status.json")
        for path in (root / "runs").glob("*/summary.receipt.json"):
            receipt = json.loads(path.read_text())
            if efast.sha256(path.parent / "phase_indices_N513.npz") != receipt["phase_indices_sha256"]:
                raise ValueError("Seed source indices changed")
            if efast.sha256(path.parent / "convergence.tsv.gz") != receipt["convergence_sha256"]:
                raise ValueError("Seed source summary changed")
        nh.collect_full_diagnostics(root)
        nh.summarize(root, efast.read_table(root / "seed_manifest.tsv"))
        nh.state(root, stage="full_slurm", status="complete", completed_designs=2500,
                 total_designs=2500, completed_tasks=35000, n_fit_seeds=500)


def retry_spec(root):
    """Reuse valid full phases; reanalyze any completed old-version trajectory."""
    ids, old_rows = [], []
    for task_id in range(1, 35001):
        path = root / "slurm" / "task_receipts" / ("task_%05d.json" % task_id)
        previous = json.loads(path.read_text()) if path.exists() else {}
        if previous.get("status") == "complete" and previous.get("reused_full_phase"):
            continue
        ids.append(task_id)
        old_rows.append(dict(task_id=task_id, previous_status=previous.get("status", "missing_or_cancelled"),
                             previous_job_id=previous.get("slurm_job_id", ""), previous_error=previous.get("error", "")))
    if not ids:
        raise ValueError("No trajectory tasks require a retry")
    efast.write_table(root / "slurm" / "retry_task_manifest.tsv", list(old_rows[0]), old_rows)
    groups, start, end = [], ids[0], ids[0]
    for value in ids[1:]:
        if value == end + 1:
            end = value
        else:
            groups.append(str(start) if start == end else "%d-%d" % (start, end))
            start = end = value
    groups.append(str(start) if start == end else "%d-%d" % (start, end))
    nh.atomic_json(root / "slurm" / "retry_request.json", dict(recorded_utc=now(),
        n_tasks=len(ids), n_preserved_tasks=35000-len(ids), array_spec=",".join(groups),
        reason="Initial shared sourceCpp lock contention; reuse valid results and isolate the Rcpp cache per task."))
    print(",".join(groups))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("prepare", "worker", "seed-summary", "status", "finalize", "retry-spec"):
        p = sub.add_parser(name)
        p.add_argument("--out-dir", required=True)
        if name == "prepare":
            p.add_argument("--fit-root", required=True)
            p.add_argument("--figure4-dir", required=True)
            p.add_argument("--sif", required=True)
            p.add_argument("--replace-plan", action="store_true")
        elif name == "worker":
            p.add_argument("--task-id", type=int, required=True)
        elif name == "seed-summary":
            p.add_argument("--seed-number", type=int, required=True)
    args = parser.parse_args()
    root = Path(args.out_dir).resolve()
    if args.command == "prepare": prepare(args)
    elif args.command == "worker": worker(root, args.task_id)
    elif args.command == "seed-summary": seed_summary(root, args.seed_number)
    elif args.command == "status": status(root)
    elif args.command == "retry-spec": retry_spec(root)
    else: finalize(root)


if __name__ == "__main__":
    main()
