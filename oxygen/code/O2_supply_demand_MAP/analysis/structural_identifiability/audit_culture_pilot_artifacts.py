#!/usr/bin/env python3
"""Verify accepted pilot payloads and collect only this study's Slurm costs."""
import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify(study, require_terminal=False):
    state_bytes = (study / "state.json").read_bytes()
    state = json.loads(state_bytes)
    active = [tid for tid, task in state["tasks"].items() if task["status"] not in ("accepted", "failed", "rejected")]
    if require_terminal and (active or state["status"] not in
            ("pilot_fitting_finished", "baseline_optimizer_unresolved", "baseline_repair_budget_exhausted")):
        raise RuntimeError("Final reporting requires a terminal study and no active fits")
    checked = []
    for tid, task in sorted(state["tasks"].items()):
        if task["status"] != "accepted":
            continue
        folder = study / "fits" / tid
        result = read(folder / "validated_result.json")
        if not result["accepted"] or not result["launch_completed"]:
            raise RuntimeError("Unaccepted payload: " + tid)
        files = []
        for name, expected in sorted(result["files_sha256"].items()):
            if Path(name).name != name:
                raise RuntimeError("Unexpected payload path")
            path = folder / name
            if sha(path) != expected:
                raise RuntimeError("Accepted payload changed: " + str(path))
            files.append(dict(file=name, sha256=expected, bytes=path.stat().st_size))
        exact = folder / "artifact_validation_v05.json"
        expected_guard_hash = result.get("artifact_validation_sha256")
        accepted_guard_name = result.get("artifact_validation_file", "artifact_validation.json")
        if Path(accepted_guard_name).name != accepted_guard_name:
            raise RuntimeError("Unexpected guard path")
        if expected_guard_hash and sha(folder / accepted_guard_name) != expected_guard_hash:
            raise RuntimeError("Full-precision guard changed: " + tid)
        guard = read(exact)
        if not guard["pass"] or guard["target_constraint_tolerance"] != 0:
            raise RuntimeError("Full-precision artifact check failed: " + tid)
        if abs(guard["objective"] - result["objective"]) > 1e-8:
            raise RuntimeError("Saved-object objective mismatch: " + tid)
        supplement = folder / "validation_supplement_v05.json"
        if supplement.exists():
            evidence = read(supplement)
            if (evidence["original_validated_result_sha256"] != sha(folder / "validated_result.json") or
                    evidence["artifact_validation_sha256"] != sha(exact)):
                raise RuntimeError("Full-precision supplement changed: " + tid)
        elif not (expected_guard_hash and accepted_guard_name == exact.name):
            raise RuntimeError("Missing full-precision guard fingerprint: " + tid)
        checked.append(dict(task_id=tid, objective=result["objective"],
            result_sha256=sha(folder / "validated_result.json"), exact_guard_sha256=sha(exact), files=files))
    record = dict(generated_utc=datetime.now(timezone.utc).isoformat(), pass_audit=True,
        state_sha256=hashlib.sha256(state_bytes).hexdigest(), execution_status=state["status"], active_fits=active,
        terminal_required=require_terminal, accepted_fits=len(checked),
        payload_files=sum(len(x["files"]) for x in checked),
        verified_bytes=sum(f["bytes"] for x in checked for f in x["files"]), fits=checked)
    snapshot = study / "manifests/accepted_artifact_state_snapshot.json"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    snapshot.write_bytes(state_bytes)
    write(study / "reports/accepted_artifact_audit.json", record)
    return record


def collect_accounting(study, slurm):
    chosen = {}
    paths = list((study / "tasks").glob("*/submission.json"))
    paths += list((study / "controllers").rglob("submission.json"))
    paths += list((study / "manifests").glob("*submission*.json"))
    for path in paths:
        record = read(path)
        value = str(record.get("job_id", record.get("stdout", ""))).strip()
        match = re.fullmatch(r"(\d+)(?:;[\w.-]+)?", value)
        if not match:
            continue
        job = match.group(1)
        relative = str(path.relative_to(study))
        kind = "fit" if relative.startswith("tasks/") else "controller" if relative.startswith("controllers/") else "validation_or_reporting"
        if job in chosen and chosen[job]["kind"] != kind:
            raise RuntimeError("Conflicting job role: " + job)
        chosen[job] = dict(kind=kind, receipt=relative)
    if not chosen:
        raise RuntimeError("No submitted jobs found")
    fields = "JobID,JobName,State,ExitCode,ElapsedRaw,AllocCPUS,MaxRSS,TotalCPU,Start,End"
    command = [str(slurm / "sacct"), "-j", ",".join(sorted(chosen)), "-P", "--format=" + fields]
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if process.returncode:
        raise RuntimeError(process.stderr)
    rows = list(csv.DictReader(io.StringIO(process.stdout), delimiter="|"))
    rows = [{key: value for key, value in row.items() if key} for row in rows]
    parents = {row["JobID"]: row for row in rows if row["JobID"] in chosen}
    totals = {}
    output = []
    for job, metadata in sorted(chosen.items()):
        row = parents.get(job)
        if row is None:
            continue
        row = dict(row, **metadata)
        elapsed = int(row["ElapsedRaw"]) if row["ElapsedRaw"].isdigit() else None
        cores = int(row["AllocCPUS"]) if row["AllocCPUS"].isdigit() else None
        row["allocated_core_hours"] = elapsed * cores / 3600 if elapsed is not None and cores is not None else None
        output.append(row)
        if row["allocated_core_hours"] is not None:
            totals[metadata["kind"]] = totals.get(metadata["kind"], 0.) + row["allocated_core_hours"]
    table = study / "tables/scheduler_accounting.tsv"
    table.parent.mkdir(exist_ok=True)
    with table.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    write(study / "manifests/scheduler_accounting_raw.json", dict(command=command, records=rows))
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(), submitted_job_ids=len(chosen),
        accounted_parent_jobs=len(parents), missing_job_ids=sorted(set(chosen) - set(parents)),
        allocated_core_hours_by_role=totals, total_allocated_core_hours=sum(totals.values()),
        unfinished_job_ids=[r["JobID"] for r in output if r["State"] in ("PENDING", "RUNNING", "COMPLETING", "CONFIGURING")],
        scope="Only job IDs in this study's submission receipts; parent allocations counted once; live totals are interim")
    write(study / "reports/resource_summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--require-terminal", action="store_true")
    parser.add_argument("--accounting", action="store_true")
    parser.add_argument("--accounting-only", action="store_true", help="Refresh costs without repeating the payload audit")
    parser.add_argument("--slurm", type=Path, default=Path("/cm/shared/apps/slurm/current/bin"))
    args = parser.parse_args()
    study = args.study.resolve()
    if args.accounting_only and not args.accounting:
        parser.error("--accounting-only requires --accounting")
    if not args.accounting_only:
        result = verify(study, args.require_terminal)
        print(json.dumps({k: result[k] for k in ("pass_audit", "accepted_fits", "payload_files", "verified_bytes", "execution_status")}, indent=2))
    if args.accounting:
        print(json.dumps(collect_accounting(study, args.slurm), indent=2))


if __name__ == "__main__":
    main()
