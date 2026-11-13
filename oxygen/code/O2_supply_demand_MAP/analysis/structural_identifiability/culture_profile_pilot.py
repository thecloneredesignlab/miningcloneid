#!/usr/bin/env python3
"""Manifest-driven culture pilot; Slurm manages all ready jobs without throttling.

Uses only the Python standard library (including on the HPC's Python 3.6).
Model evaluations and optimization remain in the frozen production R backend.
"""
import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime, timezone

SLURM = "/cm/shared/apps/slurm/current/bin"
REL = "oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability"
TARGETS = ("p_misseg", "buffer_beta")
SYMBOLS = ("lam_max", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta", "buffer_n_exp", "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O", "p_mis_base", "sigma_growth", "sigma_kary", "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=".atomic_", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write("\n")
        os.replace(temp, str(path))
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        out = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        out.writeheader()
        out.writerows(rows)


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def profile_grid(target, optimum):
    small, upper = (1e-4, .4) if target == "p_misseg" else (1e-3, 10.)
    values = [0., float(optimum)] + [math.exp(math.log(small) + i * (math.log(upper) - math.log(small)) / 7) for i in range(8)]
    values[-1] = upper
    values[2] = small
    values = sorted(set(values))
    while len(values) < 10:
        positive = [x for x in values if x > 0]
        a, b = max(zip(positive, positive[1:]), key=lambda ab: math.log(ab[1] / ab[0]))
        values = sorted(set(values + [math.sqrt(a * b)]))
    return values


def baseline_gate(results, old_objective, agreement=1e-3):
    if not results:
        return False, "no_valid_baseline"
    ordered = sorted(results, key=lambda x: x["objective"])
    best = ordered[0]
    if best["objective"] > old_objective + 1e-8:
        return False, "worse_than_feasible_archived_endpoint"
    eligible = [x for x in ordered if (x.get("agreement_eligible") if x.get("agreement_eligible") is not None else x["role"] in ("baseline_initial", "baseline_independent_repeat"))
                and x["objective"] - best["objective"] <= agreement]
    if len(set(x["primary_family"] for x in eligible)) < 2:
        return False, "independent_start_agreement_unresolved"
    strong = [x for x in ordered if x.get("role") == "baseline_strengthen"]
    if not strong:
        return False, "strengthening_required"
    last = max(strong, key=lambda x: x["ordinal"])
    if last["origin_best_q"] - last["objective"] > agreement:
        return False, "strengthening_still_improving"
    if last["objective"] - best["objective"] > agreement:
        return False, "strengthening_not_at_current_best"
    return True, "independent_agreement_and_strengthening_pass"


def envelope(results, baseline):
    if not results:
        return {"objective": None, "delta": None, "winner": None, "agreement": False}
    valid = sorted(results, key=lambda x: x["objective"])
    best = valid[0]
    return {"objective": best["objective"], "delta": best["objective"] - baseline,
            "winner": best["task_id"], "agreement": len(valid) >= 2 and valid[1]["objective"] - best["objective"] <= 1e-3}


def resource_args(task):
    return ["--partition=red", "--qos=xxlarge", "--ntasks=1", "--cpus-per-task=8", "--mem=16G",
            "--time=" + ("12:00:00" if task["reserve"] else "08:00:00"), "--no-requeue"]


def target_transport_agrees(reported, expected):
    # Runtime/RDS exactness is verified separately in R. JSON digits can round.
    return reported == 0 if expected == 0 else math.isfinite(reported) and abs(reported - expected) <= 1e-12 * max(1, abs(expected))


class Pilot:
    def __init__(self, project, study):
        self.project = Path(project).resolve()
        self.study = Path(study).resolve()
        self.config = read(self.study / "study.json")
        self.state_path = self.study / "state.json"
        self.state = read(self.state_path) if self.state_path.exists() else {
            "status": "prepared", "tasks": {}, "points": {}, "next_ordinal": 1,
            "fit_core_hours_reserved": 0, "nonfit_core_hours_reserved": 8., "events": []}
        self.script = Path(__file__).resolve()
        self.wrapper = self.study / "implementation" / REL / "run_culture_profile_pilot.sh"
        self.bank = {p.stem: read(p) for p in (self.study / "start_bank").glob("*.json")}

    def save(self):
        self.state["updated_utc"] = now()
        write(self.state_path, self.state)

    def event(self, kind, **details):
        self.state["events"].append(dict(kind=kind, utc=now(), **details))

    def verify(self):
        c = self.config
        if not c["batch_profile_authorized"] or c["qos"] != "xxlarge" or c["inference_calibration"] != "absent":
            raise RuntimeError("Wrong authorization or study")
        for row in rows(self.study / "manifests/code_manifest.tsv"):
            if sha(self.study / "implementation" / row["path"]) != row["sha256"]:
                raise RuntimeError("Frozen code changed: " + row["path"])
        validated = read(self.study / "tests/pilot_validation_status.json")
        if not validated["pass"] or validated["code_manifest_sha256"] != sha(self.study / "manifests/code_manifest.tsv"):
            raise RuntimeError("Current frozen code is not validated")
        if sha(self.study / "study.json") != read(self.study / "manifests/approved_inputs.json")["study_sha256"]:
            raise RuntimeError("Study configuration changed")
        registry = self.study / "manifests/orchestrator_versions.json"
        if registry.exists():
            active = read(registry)["active"]
            if sha(self.script) != active["sha256"]:
                raise RuntimeError("Use the registered active orchestration version")
            if active.get("artifact_checker_sha256") and sha(self.script.with_name("check_culture_profile_artifact.R")) != active["artifact_checker_sha256"]:
                raise RuntimeError("Saved-object checker changed")

    def accepted(self, stage=None, point_id=None):
        result = []
        for tid, item in self.state["tasks"].items():
            if item["status"] != "accepted" or (stage and item["stage"] != stage) or (point_id and item.get("point_id") != point_id):
                continue
            value = read(self.study / "fits" / tid / "validated_result.json")
            precise = self.study / "fits" / tid / "artifact_validation_v05.json"
            if precise.exists():
                exact = read(precise)
                if not exact["pass"]:
                    raise RuntimeError("Full-precision parameter export is invalid")
                value["params"].update({name: float(number) for name, number in exact["numeric_params_exact"].items()})
            value.update({k: item.get(k) for k in ["primary_family", "role", "ordinal", "origin_best_q", "agreement_eligible"]})
            result.append(value)
        return result

    def terminal(self, item):
        return item["status"] in ("accepted", "failed", "rejected")

    def remaining_reserve(self):
        return 20 - sum(x["reserve"] for x in self.state["tasks"].values())

    def project_params(self, params, fixed):
        p = dict(params)
        p.update(fixed)
        if fixed.get("p_misseg") == 0:
            p["k_o_mis"] = self.bank["seed144"]["k_o_mis"]
        if fixed.get("buffer_beta") == 0:
            p["buffer_n_exp"] = self.bank["seed144"]["buffer_n_exp"]
        return p

    def unique_start(self, requested, fixed, used):
        candidates = [requested] + [(k, v) for k, v in sorted(self.bank.items()) if k.startswith("independent_")]
        for name, vector in candidates:
            key = tuple(self.project_params(vector, fixed)[n] for n in SYMBOLS)
            if key not in used:
                used.append(key)
                return name, vector
        raise RuntimeError("Distinct start bank exhausted")

    def submit_command(self, args, directory):
        directory.mkdir(parents=True, exist_ok=True)
        intent = directory / "submission_intent.json"
        completed = directory / "submission.json"
        if completed.exists():
            return read(completed)["job_id"]
        if intent.exists():
            raise RuntimeError("Ambiguous submission requires scheduler audit: " + str(directory))
        write(intent, {"arguments": args, "utc": now()})
        p = subprocess.run([SLURM + "/sbatch", "--parsable"] + args, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
        if p.returncode:
            write(directory / "submission_error.json", {"stderr": p.stderr, "stdout": p.stdout})
            raise RuntimeError("sbatch rejected task: " + p.stderr)
        job = p.stdout.strip().split(";")[0]
        if not job.isdigit():
            raise RuntimeError("Ambiguous sbatch output: " + p.stdout)
        write(completed, {"job_id": job, "utc": now(), "stdout": p.stdout})
        return job

    def create_task(self, stage, role, starts, family, point_id=None, reserve=False, origin_best_q=None):
        if len(self.state["tasks"]) >= 125 or (reserve and self.remaining_reserve() < 1):
            return None
        extra = 96 if reserve else 64
        if self.state["fit_core_hours_reserved"] + extra > 8640:
            raise RuntimeError("Fit budget exhausted")
        ordinal = self.state["next_ordinal"]
        tid = "fit_{:03d}".format(ordinal)
        point = self.state["points"].get(point_id)
        task = dict(task_id=tid, stage=stage, role=role, ordinal=ordinal, reserve=reserve,
                    primary_family=family, point_id=point_id, origin_best_q=origin_best_q,
                    fixed={point["target"]: point["value"]} if point else {},
                    seed=self.config["master_seed"] + ordinal * 1009,
                    DE_itermax=1000 if reserve else 500, local_maxit=100 if reserve else 50,
                    starts_labels=list(starts), created_utc=now(), controller_sha256=sha(self.script),
                    agreement_eligible=role in ("baseline_initial", "baseline_independent_repeat") and len(starts) == 1)
        directory = self.study / "tasks" / tid
        if directory.exists():
            raise RuntimeError("Unregistered existing task requires submission audit: " + tid)
        write(directory / "task.json", task)
        write(directory / "starts.json", starts)
        args = resource_args(task) + ["--job-name=ci1109_" + tid, "--output=" + str(self.study / "logs" / (tid + "_%j.out")),
                                     "--error=" + str(self.study / "logs" / (tid + "_%j.err")),
                                     str(self.wrapper), str(self.project), str(self.study), "fit", tid]
        job = self.submit_command(args, directory)
        item = dict(task, job_id=job, status="submitted")
        self.state["tasks"][tid] = item
        self.state["next_ordinal"] += 1
        self.state["fit_core_hours_reserved"] += extra
        if point:
            point["task_ids"].append(tid)
        self.event("submitted", task_id=tid, job_id=job, stage=stage, reserve=reserve)
        self.save()
        self.ensure_controller(item)
        return tid

    def ensure_controller(self, task):
        if self.terminal(task):
            return
        directory = self.study / "controllers" / "v07" / task["task_id"]
        if (directory / "submission.json").exists():
            return
        if self.state["nonfit_core_hours_reserved"] + 1 / 6 > 52:
            raise RuntimeError("Controller/test budget exhausted; 8 core-hours retained for reporting")
        args = ["--partition=red", "--qos=xxlarge", "--ntasks=1", "--cpus-per-task=1", "--mem=4G", "--time=00:10:00", "--no-requeue",
                "--dependency=afterany:" + task["job_id"], "--job-name=ci1109_reduce_" + task["task_id"],
                "--output=" + str(self.study / "logs" / ("reduce_" + task["task_id"] + "_%j.out")),
                "--error=" + str(self.study / "logs" / ("reduce_" + task["task_id"] + "_%j.err")),
                "--wrap=python3 " + str(self.script) + " advance " + str(self.project) + " " + str(self.study)]
        job = self.submit_command(args, directory)
        self.state["nonfit_core_hours_reserved"] += 1 / 6
        self.event("controller_submitted", task_id=task["task_id"], job_id=job)
        self.save()

    def accounting(self):
        jobs = [x["job_id"] for x in self.state["tasks"].values() if not self.terminal(x)]
        if not jobs:
            return {}
        args = [SLURM + "/sacct", "-j", ",".join(jobs), "--format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS,MaxRSS,TotalCPU", "-P"]
        p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
        if p.returncode:
            raise RuntimeError("Accounting unavailable: " + p.stderr)
        values = list(csv.DictReader(p.stdout.splitlines(), delimiter="|"))
        write(self.study / "manifests/latest_scheduler_query.json", {"utc": now(), "rows": values})
        result = {}
        for job in jobs:
            result[job] = [x for x in values if x["JobID"] == job or x["JobID"].startswith(job + ".")]
        return result

    def check_saved_objects(self, item):
        tid = item["task_id"]
        out = self.study / "fits" / tid
        checker = self.script.with_name("check_culture_profile_artifact.R")
        home = self.study / "runtime" / "artifact_checks" / tid
        home.mkdir(parents=True, exist_ok=True)
        command = ["apptainer", "exec", "--cleanenv", "--home", str(home), "--bind", "/share:/share",
                   "--env", "RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript", "--env", "OMP_NUM_THREADS=1",
                   "--env", "OPENBLAS_NUM_THREADS=1", "--env", "R_ENVIRON_USER=/dev/null", "--env", "R_PROFILE_USER=/dev/null",
                   self.config["image_path"], "bash", str(self.study / "implementation/scripts/agentRrunner.sh"),
                   str(checker), str(self.project), str(self.study), tid]
        checked = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=120)
        (out / "artifact_validation_v05.log").write_text(checked.stdout + checked.stderr)
        if checked.returncode or not read(out / "artifact_validation_v05.json")["pass"]:
            raise RuntimeError("External saved-object validation failed")
        return sha(out / "artifact_validation_v05.json")

    def validate_result(self, item, accounting):
        tid = item["task_id"]
        out = self.study / "fits" / tid
        result = read(out / "candidate_result.json")
        task_path = self.study / "tasks" / tid / "task.json"
        inputs = read(out / "run_inputs.json")
        if inputs["task_id"] != tid or result["task_id"] != tid or result["slurm_job_id"] != item["job_id"]:
            raise RuntimeError("Wrong task/job identity")
        paths = {"task_sha256": task_path, "starts_sha256": self.study / "tasks" / tid / "starts.json",
                 "config_sha256": self.study / "study.json", "code_manifest_sha256": self.study / "manifests/code_manifest.tsv",
                 "wp0_contract_sha256": self.project / self.config["wp0_relative_path"] / "manifests/production_contract.json"}
        if any(sha(path) != inputs[name] for name, path in paths.items()) or sha(out / "run_inputs.json") != result["fingerprint"]:
            raise RuntimeError("Stale input fingerprint")
        for name, expected in result["files_sha256"].items():
            if sha(out / name) != expected:
                raise RuntimeError("Result artifact changed: " + name)
        if not result["samples_valid"] or not result["population_identical"] or result["sample_counts"] != {"growth": 114, "karyotype": 12, "flow": 20}:
            raise RuntimeError("Samples or actual initial population mismatch")
        q = result["objective"]
        if not math.isfinite(q) or any(abs(q - result[n]) > 1e-8 for n in ["reevaluated_objective", "independent_objective"]):
            raise RuntimeError("Independent scoring mismatch")
        if abs(q + sum(result["component_means"].values())) > 1e-8:
            raise RuntimeError("Objective ledger does not sum to the data-only score")
        parameters = rows(out / "parameter_result.tsv")
        if set(x["parameter"] for x in parameters) != set(SYMBOLS):
            raise RuntimeError("Wrong parameter map")
        for p in parameters:
            value = result["params"][p["parameter"]]
            lower, upper = float(p["lower"]), float(p["upper"])
            tolerance = 64 * 2.220446049250313e-16 * max(1, abs(lower), abs(upper))
            if not math.isfinite(value) or value < lower - tolerance or value > upper + tolerance:
                raise RuntimeError("Nuisance bound violation")
            if p["status"] != "active" and abs(value - float(p["fixed_value"])) > tolerance:
                raise RuntimeError("Fixed/inactive mismatch")
        for target, value in item["fixed"].items():
            if not target_transport_agrees(result["params"][target], value):
                raise RuntimeError("Reported target exceeds JSON transport rounding")
        expected = 20 if not item["fixed"] else (18 if list(item["fixed"].values())[0] == 0 else 19)
        if result["active_parameters"] != expected:
            raise RuntimeError("Wrong active dimension")
        counts = {x["phase"]: int(x["evaluations"]) for x in rows(out / "evaluation_counts.tsv")}
        if counts.get("de") != result["native_de_nfeval"] or sum(counts.values()) != result["total_evaluator_invocations"]:
            raise RuntimeError("Evaluation count mismatch")
        scheduler = next(x for x in accounting if x["JobID"] == item["job_id"])
        if scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0" or int(scheduler["AllocCPUS"]) != 8:
            raise RuntimeError("External process incomplete")
        self.check_saved_objects(item)
        result["params_summary"] = dict(result["params"])
        exact = read(out / "artifact_validation_v05.json")["numeric_params_exact"]
        result["params"].update({name: float(number) for name, number in exact.items()})
        if any(result["params"][target] != value for target, value in item["fixed"].items()):
            raise RuntimeError("Full-precision RDS target export differs from the requested target")
        write(out / "scheduler.json", accounting)
        result.update(accepted=True, launch_completed=True, launch_exit_code=0,
                      search_status="budget_limited_provisional" if result["budget_limited"] else "completed_search_candidate",
                      validated_utc=now(), allocated_core_hours=8 * int(scheduler["ElapsedRaw"]) / 3600)
        result["artifact_validation_file"] = "artifact_validation_v05.json"
        result["artifact_validation_sha256"] = sha(out / "artifact_validation_v05.json")
        result["artifact_checker_sha256"] = sha(self.script.with_name("check_culture_profile_artifact.R"))
        if result["optimizer_trace"].get("local_convergence") != 0:
            result["search_status"] += "_local_optimizer_unresolved"
        write(out / "validated_result.json", result)

    def refresh(self):
        accounting = self.accounting()
        for tid, item in self.state["tasks"].items():
            if self.terminal(item):
                continue
            values = accounting.get(item["job_id"], [])
            parent = next((x for x in values if x["JobID"] == item["job_id"]), None)
            if not parent:
                continue
            state = parent["State"].split()[0]
            if state in ("PENDING", "RUNNING", "COMPLETING", "CONFIGURING", "SUSPENDED", "REQUEUED", "RESIZING"):
                item["scheduler_state"] = state
                continue
            item["scheduler_state"] = state
            item["accounting"] = values
            if state != "COMPLETED" or parent["ExitCode"] != "0:0":
                item["status"] = "failed"
                item["reason"] = "scheduler_" + state
            else:
                try:
                    self.validate_result(item, values)
                    item["status"] = "accepted"
                except Exception as error:
                    item["status"] = "rejected"
                    item["reason"] = str(error)
            self.event(item["status"], task_id=tid, reason=item.get("reason"))
        self.save()

    def baseline(self):
        values = self.accepted(stage="baseline")
        return min(values, key=lambda x: x["objective"]) if values else None

    def strengthen_baseline(self, role="baseline_strengthen", source=None):
        best = self.baseline()
        if self.remaining_reserve() < 1:
            return False
        source = source or best
        if source:
            starts = {"current_candidate": source["params"]}
            family = ("profile_solution_" + source["task_id"]) if role == "baseline_repair" else source["primary_family"]
            origin = best["objective"] if best else self.config["production_objective"]
            if best and source["params"] != best["params"]:
                starts["unrestricted_best"] = best["params"]
        else:
            starts = {"seed144": self.bank["seed144"]}
            family, origin = "production_seed144", self.config["production_objective"]
        self.create_task("baseline", role, starts, family, reserve=True, origin_best_q=origin)
        return True

    def create_grid(self):
        best = self.baseline()
        for target in TARGETS:
            values = profile_grid(target, best["params"][target])
            anchor = values.index(best["params"][target])
            for i, value in enumerate(values):
                pid = ("pm" if target == "p_misseg" else "bb") + "_{:02d}".format(i)
                parent = None if i == anchor else ("pm" if target == "p_misseg" else "bb") + "_{:02d}".format(i + (1 if i < anchor else -1))
                self.state["points"][pid] = dict(point_id=pid, target=target, value=value, parent=parent,
                                                    anchor=i == anchor, kind="initial", task_ids=[], baseline_at_creation=best["task_id"])
        self.event("grid_frozen", baseline=best["task_id"])
        self.save()

    def point_starts(self, point):
        best = self.baseline()
        alternatives = [x for x in self.accepted(stage="baseline") if x["task_id"] != best["task_id"]]
        alternate = max(alternatives, key=lambda x: sum((x["params"][n] - best["params"][n])**2 / max(1, abs(best["params"][n]))**2 for n in SYMBOLS)) if alternatives else None
        alternate_params = alternate["params"] if alternate else self.bank["seed371"]
        if point["target"] == "p_misseg" and point["value"] == 0:
            alternate_params = self.bank["engineering_null"]
        return [("unrestricted_best", best["params"]), ("alternative", alternate_params),
                ("independent_1", self.bank["independent_1"]), ("independent_2", self.bank["independent_2"])]

    def submit_points(self):
        for pid, point in self.state["points"].items():
            existing = [self.state["tasks"][tid] for tid in point["task_ids"]]
            main = [x for x in existing if x["role"].startswith("point_start_")]
            fixed = {point["target"]: point["value"]}
            used = []
            for task in main:
                start = read(self.study / "tasks" / task["task_id"] / "starts.json")
                vector = next(iter(start.values()))
                used.append(tuple(self.project_params(vector, fixed)[n] for n in SYMBOLS))
            desired = self.point_starts(point)
            for index in range(len(main), 4):
                if point["kind"] == "refinement" and self.remaining_reserve() < 5 - len(main):
                    break
                name, vector = self.unique_start(desired[index], fixed, used)
                tid = self.create_task("profile", "point_start_" + str(index + 1), {name: vector}, name, pid,
                                       reserve=point["kind"] == "refinement")
                if tid:
                    main.append(self.state["tasks"][tid])
            if len(main) != 4:
                continue
            parent = self.state["points"].get(point["parent"])
            if parent and (len(parent["task_ids"]) < 5 or any(not self.terminal(self.state["tasks"][tid]) for tid in parent["task_ids"])):
                continue
            neighbor = self.accepted(point_id=parent["point_id"]) if parent else []
            source = ("neighbor_" + parent["point_id"], min(neighbor, key=lambda x: x["objective"])["params"]) if neighbor else ("independent_3", self.bank["independent_3"])
            name, vector = self.unique_start(source, fixed, used)
            self.create_task("profile", "point_start_5", {name: vector}, name, pid, reserve=point["kind"] == "refinement")

    def initial_complete(self):
        points = [p for p in self.state["points"].values() if p["kind"] == "initial"]
        return len(points) == 20 and all(len(p["task_ids"]) >= 5 and all(self.terminal(self.state["tasks"][tid]) for tid in p["task_ids"]) for p in points)

    def reserve_checks(self):
        for target in TARGETS:
            null = next(p for p in self.state["points"].values() if p["target"] == target and p["value"] == 0)
            for role in ("null_reverse", "null_independent"):
                if any(self.state["tasks"][tid]["role"] == role for tid in null["task_ids"]):
                    continue
                if not self.remaining_reserve():
                    return
                near = sorted([p for p in self.state["points"].values() if p["target"] == target and p["value"] > 0], key=lambda p: p["value"])
                values = self.accepted(point_id=near[0]["point_id"]) if near else []
                vector = min(values, key=lambda x: x["objective"])["params"] if values and role == "null_reverse" else self.bank["independent_4"]
                self.create_task("profile", role, {role: vector}, role, null["point_id"], reserve=True)
        pending = any(not self.terminal(x) for x in self.state["tasks"].values())
        if pending:
            return
        # The legacy null_reverse role repeats positive-to-zero continuation.
        # Return from each verified null to its nearest positive grid value to
        # check the opposite direction before spending remaining reserves.
        for target in TARGETS:
            null = next(p for p in self.state["points"].values() if p["target"] == target and p["value"] == 0)
            positive = min((p for p in self.state["points"].values() if p["target"] == target and p["value"] > 0),
                           key=lambda p: p["value"])
            if any(self.state["tasks"][tid]["role"] == "point_reverse" for tid in positive["task_ids"]):
                continue
            if not self.remaining_reserve():
                break
            values = self.accepted(point_id=null["point_id"])
            if not values:
                self.event("reverse_return_unresolved", point_id=positive["point_id"], reason="no_verified_null_source")
                continue
            source = min(values, key=lambda x: x["objective"])["params"]
            name = "reverse_from_null_" + null["point_id"]
            self.create_task("profile", "point_reverse", {name: source}, name, positive["point_id"], reserve=True)
        if any(not self.terminal(x) for x in self.state["tasks"].values()):
            return
        unresolved = []
        for pid, point in self.state["points"].items():
            values = self.accepted(point_id=pid)
            e = envelope(values, self.baseline()["objective"])
            if not e["agreement"]:
                unresolved.append((point["value"] != 0, len(values) > 0, pid))
        for _, _, pid in sorted(unresolved):
            if not self.remaining_reserve():
                break
            values = self.accepted(point_id=pid)
            source = min(values, key=lambda x: x["objective"])["params"] if values else self.baseline()["params"]
            self.create_task("profile", "point_strengthen", {"current_best": source, "independent_5": self.bank["independent_5"]},
                             "point_strengthen", pid, reserve=True)
        if not unresolved:
            self.refine_grid()

    def refine_grid(self):
        # A descriptive curvature threshold for grid allocation, never a CI rule.
        available = self.remaining_reserve()
        for target in TARGETS:
            if available < 5 or any(p["kind"] == "refinement" and p["target"] == target for p in self.state["points"].values()):
                continue
            points = sorted([p for p in self.state["points"].values() if p["target"] == target], key=lambda p: p["value"])
            scale = self.config["zero_scales"][target]
            candidates = []
            for a, b, c in zip(points, points[1:], points[2:]):
                fits = [self.accepted(point_id=p["point_id"]) for p in (a, b, c)]
                if not all(fits):
                    continue
                q = [min(v["objective"] for v in group) for group in fits]
                x = [math.log1p(p["value"] / scale) for p in (a, b, c)]
                curvature = abs(q[1] - (q[0] + (q[2] - q[0]) * (x[1] - x[0]) / (x[2] - x[0])))
                if curvature > self.config["grid_refinement_score_threshold"]:
                    left, right = (a, b) if x[1] - x[0] >= x[2] - x[1] else (b, c)
                    candidates.append((curvature, left, right))
            if not candidates:
                continue
            curvature, left, right = max(candidates, key=lambda x: x[0])
            value = math.sqrt(left["value"] * right["value"]) if left["value"] else right["value"] / 2
            if any(abs(p["value"] - value) <= 1e-14 for p in points):
                continue
            pid = ("pm" if target == "p_misseg" else "bb") + "_refine"
            self.state["points"][pid] = dict(point_id=pid, target=target, value=value, parent=left["point_id"],
                                               anchor=False, kind="refinement", task_ids=[], baseline_at_creation=self.baseline()["task_id"])
            available -= 5
            self.event("refinement_registered", point_id=pid, curvature=curvature, criterion="descriptive_score_curvature_not_confidence")
        self.save()
        self.submit_points()

    def advance(self, initial=False):
        self.verify()
        self.state.pop("controller_error", None)
        for item in list(self.state["tasks"].values()):
            self.ensure_controller(item)
        self.refresh()
        for item in self.state["tasks"].values():
            out = self.study / "fits" / item["task_id"]
            if item["status"] == "accepted" and not (out / "artifact_validation_v05.json").exists():
                digest = self.check_saved_objects(item)
                write(out / "validation_supplement_v05.json", dict(task_id=item["task_id"], artifact_validation_sha256=digest,
                    original_validated_result_sha256=sha(out / "validated_result.json"),
                    artifact_checker_sha256=sha(self.script.with_name("check_culture_profile_artifact.R")), utc=now()))
        if initial and not self.state["tasks"]:
            for name in ["seed144", "seed371", "seed10", "independent_1", "independent_2"]:
                self.create_task("baseline", "baseline_initial", {name: self.bank[name]}, name)
        baselines = [x for x in self.state["tasks"].values() if x["stage"] == "baseline"]
        if not baselines or any(not self.terminal(x) for x in baselines):
            self.state["status"] = "baseline_running"
            self.save(); self.report(); return
        best = self.baseline()
        gate, reason = baseline_gate(self.accepted(stage="baseline"), self.config["production_objective"])
        self.state["baseline_gate"] = dict(pass_gate=gate, reason=reason, best_task=best["task_id"] if best else None)
        if not gate:
            used = sum(x["reserve"] for x in baselines)
            strong = [x for x in baselines if x["role"] == "baseline_strengthen"]
            if not strong and used < 6 and self.strengthen_baseline():
                self.state["status"] = "baseline_strengthening"
            elif reason == "independent_start_agreement_unresolved" and used < 6:
                prior = {x["primary_family"] for x in baselines}
                names = ["independent_" + str(i) for i in range(3, 9) if "independent_" + str(i) not in prior]
                for name in names[:min(6 - used, self.remaining_reserve())]:
                    self.create_task("baseline", "baseline_independent_repeat", {name: self.bank[name]}, name,
                                     reserve=True, origin_best_q=best["objective"] if best else None)
                self.state["status"] = "baseline_independent_repeats_running" if names and used < 6 else "baseline_optimizer_unresolved"
            elif used < 6 and self.strengthen_baseline():
                self.state["status"] = "baseline_strengthening"
            else:
                self.state["status"] = "baseline_optimizer_unresolved"
            self.save(); self.report(); return
        constrained = self.accepted(stage="profile")
        better = min(constrained, key=lambda x: x["objective"]) if constrained else None
        if better and better["objective"] < best["objective"] - 1e-8:
            if self.strengthen_baseline("baseline_repair", better):
                self.state["status"] = "baseline_repair_running"
            else:
                self.state["status"] = "baseline_repair_budget_exhausted"
            self.save(); self.report(); return
        if not self.state["points"]:
            self.create_grid()
        self.submit_points()
        if self.initial_complete():
            self.reserve_checks()
        active = any(not self.terminal(x) for x in self.state["tasks"].values())
        self.state["status"] = "profile_running" if active or not self.initial_complete() else "pilot_fitting_finished"
        self.save(); self.report()

    def report(self):
        best = self.baseline()
        summaries = []
        for pid, point in sorted(self.state["points"].items()):
            values = self.accepted(point_id=pid)
            e = envelope(values, best["objective"]) if best else envelope([], 0)
            summaries.append(dict(point_id=pid, target=point["target"], value=point["value"], kind=point["kind"],
                                  attempted=len(point["task_ids"]), accepted=len(values), **e))
        tsv(self.study / "tables/profile_points.tsv", summaries,
            ["point_id", "target", "value", "kind", "attempted", "accepted", "objective", "delta", "winner", "agreement"])
        starts = []
        for tid, item in sorted(self.state["tasks"].items()):
            record = dict(item)
            if item["status"] == "accepted":
                result = read(self.study / "fits" / tid / "validated_result.json")
                record.update({key: result[key] for key in ["objective", "native_objective", "winner_phase", "fit_wall_seconds", "total_evaluator_invocations", "allocated_core_hours", "search_status"]})
            starts.append(record)
        tsv(self.study / "tables/profile_starts.tsv", starts,
            ["task_id", "job_id", "stage", "point_id", "role", "primary_family", "seed", "reserve", "status", "reason", "objective", "native_objective", "winner_phase", "fit_wall_seconds", "total_evaluator_invocations", "allocated_core_hours", "search_status"])
        write(self.study / "reports/progress.json", {
            "status": self.state["status"], "updated_utc": now(), "tasks": len(starts),
            "accepted": sum(x["status"] == "accepted" for x in starts),
            "active": sum(not self.terminal(x) for x in starts),
            "failed_or_rejected": sum(x["status"] in ("failed", "rejected") for x in starts),
            "baseline": {"task_id": best["task_id"], "objective": best["objective"]} if best else None,
            "remaining_reserve_fits": self.remaining_reserve(), "fit_core_hours_reserved": self.state["fit_core_hours_reserved"],
            "nonfit_core_hours_reserved": self.state["nonfit_core_hours_reserved"],
            "qos": "xxlarge", "application_concurrency_limit": None,
            "controller_program": str(self.script), "controller_sha256": sha(self.script),
            "inference_calibration": "absent", "points": summaries})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["submit", "advance", "status"])
    parser.add_argument("project")
    parser.add_argument("study")
    args = parser.parse_args()
    if args.command == "status":
        pilot = Pilot(args.project, args.study)
        print(json.dumps(read(pilot.study / "reports/progress.json"), indent=2))
        return
    # Serializes short metadata updates only; it does not throttle fit jobs.
    with (Path(args.study).resolve() / "controller.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        # Load the latest state only after acquiring the metadata lock.
        pilot = Pilot(args.project, args.study)
        try:
            pilot.advance(initial=args.command == "submit")
        except Exception as error:
            pilot.event("controller_error", error=str(error))
            pilot.state["controller_error"] = str(error)
            pilot.save()
            raise
    print(json.dumps(read(pilot.study / "reports/progress.json"), indent=2))


if __name__ == "__main__":
    main()
