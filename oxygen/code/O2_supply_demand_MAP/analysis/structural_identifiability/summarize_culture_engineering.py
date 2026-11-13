#!/usr/bin/env python3
"""Report the tested v2 mapping and one externally validated cost benchmark."""
import csv
import hashlib
import json
import sys
from collections import Counter
from math import isfinite
from pathlib import Path

project = Path(sys.argv[1]).resolve()
study = Path(sys.argv[2]).resolve()
out = study / "benchmark/attempt_01"


def load_json(path):
    return json.loads(path.read_text())


def rows(path, delimiter="\t"):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def seconds(text):
    days, value = text.split("-", 1) if "-" in text else ("0", text)
    parts = [float(x) for x in value.split(":")]
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h, m, s = 0, *parts
    else:
        h, m, s = 0, 0, parts[0]
    return float(days) * 86400 + h * 3600 + m * 60 + s


def rss_mib(text):
    if not text:
        return None
    unit = text[-1]
    value = float(text[:-1]) if unit.isalpha() else float(text)
    return value * {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024**2}.get(unit, 1 / 1024**2)


config = load_json(study / "study.json")
validation = load_json(study / "tests/validation_status.json")
result = load_json(out / "validated_result.json")
assert validation["pass"] and validation["failed"] == 0
assert result["accepted"] and result["launch_completed"] and result["launch_exit_code"] == 0
assert result["params"]["p_misseg"] == 0
assert result["active_parameters"] == 18
for name, expected in validation["module_hashes"].items():
    actual = hashlib.sha256((project / "oxygen/code/O2_supply_demand_MAP/util" / name).read_bytes()).hexdigest()
    assert actual == expected, name
scheduler = rows(study / "manifests/benchmark_scheduler.tsv", delimiter="|")
job = next(x for x in scheduler if x["JobID"] == str(result["slurm_job_id"]))
assert job["State"] == "COMPLETED" and job["ExitCode"] == "0:0"
assert int(job["AllocCPUS"]) == config["benchmark"]["n_cores"] == 8
assert result["frozen_input_files_verified_unchanged"] == 651
rss_values = [rss_mib(x.get("MaxRSS", "")) for x in scheduler]
rss_values = [x for x in rss_values if x is not None]
max_rss = max(rss_values) if rss_values else None
elapsed = seconds(job["Elapsed"])
fit_seconds = result["fit_wall_seconds"]
phase = rows(out / "evaluation_counts.tsv")
evaluations = rows(out / "evaluations.tsv")
assert len(evaluations) == result["total_evaluator_invocations"]
assert sum(int(x["evaluations"]) for x in phase) == len(evaluations)
failures = [x for x in evaluations if x.get("error", "") not in ("", "NA")]
assert len(failures) == result["failed_evaluator_invocations"]
failure_types = Counter(x["error"].split(":", 1)[0] for x in failures)
failed_by_phase = Counter(x["phase"] for x in failures)
evaluation_seconds_by_phase = Counter()
minimum_logged_score_by_phase = {}
for evaluation in evaluations:
    evaluation_seconds_by_phase[evaluation["phase"]] += float(evaluation["seconds"])
    try:
        score = float(evaluation["objective"])
    except ValueError:
        continue
    if isfinite(score):
        name = evaluation["phase"]
        minimum_logged_score_by_phase[name] = min(minimum_logged_score_by_phase.get(name, float("inf")), score)
# DEoptim serializes closures repeatedly; their captured local counters can
# reset. Physical append-log rows count calls, not unique local counter values.
pid_rows = {}
indexed = []
for evaluation in evaluations:
    pid = evaluation["pid"]
    pid_rows[pid] = pid_rows.get(pid, 0) + 1
    indexed.append({**evaluation, "pid_log_row": pid_rows[pid],
                    "event_id": f"{result['slurm_job_id']}:{pid}:{pid_rows[pid]}"})
assert len({x["event_id"] for x in indexed}) == len(evaluations)
with (out / "evaluations_indexed.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(indexed[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(indexed)
parameters = rows(out / "parameter_result.tsv")
active_bounds = [x["parameter"] for x in parameters if x["status"] == "active" and (x["at_lower"] == "TRUE" or x["at_upper"] == "TRUE")]
mapping_table = []
for variant in ("production", "zero_extended", "no_inducible", "flat_survival"):
    table = rows(study / "manifests" / f"parameter_status_{variant}.tsv")
    mapping_table.append({"variant": variant, **{status: sum(x["status"] == status for x in table) for status in ("active", "fixed", "inactive")}})
with (study / "manifests/mapping_counts.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(mapping_table[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(mapping_table)
mapping_md = "\n".join(f"| {x['variant']} | {x['active']} | {x['fixed']} | {x['inactive']} |" for x in mapping_table)
phase_md = "\n".join(f"| {x['phase']} | {x['evaluations']} | {failed_by_phase[x['phase']]} | {evaluation_seconds_by_phase[x['phase']]:.3f} |" for x in phase)
trace = result["optimizer_trace"]
local_code = trace.get("local_convergence")
de_counts = next((int(x["evaluations"]) for x in phase if x["phase"] == "de"), 0)
native_de = result.get("native_de_nfeval")
if native_de is not None:
    assert native_de == de_counts, (native_de, de_counts)
summary = {
    "study_id": config["study_id"],
    "engineering_tests_pass": True,
    "test_cases": validation["test_cases"],
    "passed_expectations": validation["passed_expectations"],
    "real_constrained_fits": 1,
    "fixed": {"p_misseg": 0},
    "active_parameters": 18,
    "slurm_job_id": result["slurm_job_id"],
    "launch_status": job["State"],
    "fit_wall_seconds": fit_seconds,
    "job_wall_seconds": elapsed,
    "slurm_max_rss_MiB": max_rss,
    "allocated_core_hours": elapsed * 8 / 3600,
    "evaluator_invocations": len(evaluations),
    "failed_evaluator_invocations": result["failed_evaluator_invocations"],
    "failure_categories": dict(failure_types),
    "summed_evaluation_seconds_by_phase": dict(evaluation_seconds_by_phase),
    "minimum_logged_score_by_phase": minimum_logged_score_by_phase,
    "DE_iterations": trace["iter_completed"],
    "local_convergence_code": local_code,
    "objective": result["objective"],
    "replay_difference": result["reevaluated_objective"] - result["objective"],
    "independent_score_difference": result["independent_objective"] - result["objective"],
    "search_status": result["search_status"],
    "global_optimum_certified": False,
    "extended_unrestricted_baseline_refitted": False,
    "batch_profile_launched": False,
    "active_bound_hits": active_bounds,
    "historical_container_identity_caveat_retained": True,
    "previous_WP0_package_resolution_metadata_corrected": True,
}
(study / "reports/engineering_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
report = f"""# Culture v2: parameter mapping, exact nulls and one constrained-fit cost benchmark

WP1 completed parameter mapping, both exact-null restrictions, complete start injection and necessary validation, followed by one biological constrained optimization. At WP1 completion no batch profile, COMBOS study, manuscript rebuild or production-fit modification had been launched. Subsequent pilot work has separate manifests.

## Verified implementation

Production provenance is the iteration5 seed144 bundle at commit 1a5d861c9c4d2af72ea0d07def0023831cc531e4. The isolated adapter uses the frozen v2 backend's DEoptim, worker/DLL preflight, L-BFGS-B, protocol barriers, simulation and scoring. The legacy culture CLI was not globally upgraded; diagnostic runs must explicitly use the v2 adapter.

| Domain or restriction | Active | Fixed | Inactive |
|---|---:|---:|---:|
{mapping_md}

The extended domain admits p_misseg in [0,0.4] and buffer_beta in [0,10], with fixed transforms log1p(p/0.015) and log1p(beta/1). Other production transforms and bounds remain unchanged. Fixed zeros reconstruct exact natural-scale zero, without epsilon substitution.

At p_misseg=0, baseline missegregation, WGD, growth/death, survival, both initial distributions and both fitted observation-error scales remain free; k_o_mis is data-inactive. At buffer_beta=0, buffer_smax remains free and error-size dependence remains; buffer_n_exp is data-inactive. The generic mapper retains data-inactive coordinates if they contribute to a complete-objective penalty. This culture objective has no such priors or regularization terms.

The archived mu_hp log-transform round trip lies about 8.7e-19 below its lower bound. A 64-machine-epsilon tolerance canonicalizes only numerical boundary roundoff and records roundoff_to_bound; substantive bound violations are rejected. Exact target zeros are not relaxed.

## Validation

The final HPC validation passed {validation['test_cases']} test cases and {validation['passed_expectations']:,} assertions, with zero failures. Coverage includes parameter/domain round trips, inactive versus penalty-active parameters, exact zeros and PSOCK transport, missing-sample rejection, wrong-target/stale/interrupted-result rejection, fingerprints and overwrite protection, an analytic nuisance-optimum oracle, existing kernel/event-accounting tests, and integrated production/extended endpoint replays.

Both integrated replays retain 131 unique IDs, 114 growth records, 12 chromosome samples/220 cells and 20 flow distributions. Complete live-state trajectories, distributions and selected days are compared with the archived reference. Fixing a target does not alter cfg or prior fields.

The start-injection test observed the 256-by-18 population actually received by native DEoptim. Only that analytic test substituted a cheap objective; the measured constrained fit used the complete production model.

## Measured constrained fit

Slurm job {result['slurm_job_id']} was externally confirmed COMPLETED with exit 0. It requested eight cores and 16 GB, with single-thread BLAS and the WP0 SIF fingerprint. The C++ cache was already populated by validation and reused, so the cost represents a cache-available subsequent fit. Slurm MaxRSS and allocated core-hours are distinct from per-process RSS and consumed CPU time.

The fit fixed p_misseg=0 and optimized 18 active coordinates. Complete seed144 and seed325 vectors, after explicit projection, were the first two population members; the remaining members were independently uniform in transformed coordinates. Two injected members are not two optimization runs or proven distinct basins.

| Item | Measured value or setting |
|---|---:|
| DE population | 256 |
| DE generation cap / completed | 100 / {trace['iter_completed']} |
| Local iteration cap | 50 |
| Backend wall time | {fit_seconds:.3f} seconds ({fit_seconds / 60:.2f} minutes) |
| Total Slurm wall time | {elapsed:.0f} seconds |
| Slurm MaxRSS | {max_rss if max_rss is not None else 'unavailable'} MiB |
| Allocated core-hours | {elapsed * 8 / 3600:.4f} |
| Actual evaluator calls | {len(evaluations):,} |
| Failed evaluator calls | {result['failed_evaluator_invocations']:,} |
| Local convergence code | {local_code} |

| Phase | Evaluator calls | Failed calls | Sum of call elapsed seconds |
|---|---:|---:|---:|
{phase_md}

Counts use physical append-log rows, including preflight and verification. DEoptim can repeatedly serialize closures and reset their local counters; these counters are not deduplicated. evaluations_indexed.tsv adds unique run/PID/row event IDs. DE counts also match native nfeval. Logging and retained-sample checks are included in the measured cost.

Elapsed seconds are summed across calls; parallel DE sums are not phase wall times. Local refinement uses native serial numerical gradients, so DE throughput alone does not predict complete-fit cost. Failure categories: {dict(failure_types)}. Protocol-infeasible vectors receive production barriers; missing samples are not deleted to produce scores. Different failure proportions can change runtime.

## Endpoint validation and interpretation

Saved Q={result['objective']:.15g}; independent forward-replay difference={result['reevaluated_objective'] - result['objective']:.5g}; independent arithmetic difference={result['independent_objective'] - result['objective']:.5g}. The full runtime p_misseg is exactly zero, nuisances satisfy bounds, and every required observation is retained. Actual/requested populations, input/result fingerprints and external exit status agree.

Search status: `{result['search_status']}`. Active bound hits: {', '.join(active_bounds) if active_bounds else 'none'}. Process success does not certify optimization convergence or a global minimum.

The local log contains a minimum Q={minimum_logged_score_by_phase.get('local', float('nan')):.15g}, slightly below the native returned endpoint. Complete vectors were not saved for every probe, so that logged minimum cannot be independently certified. WP1 retained the production endpoint-selection rule and preserves this discrepancy as a search-quality diagnostic; subsequent pilot bookkeeping should save improved full candidates.

This was a bounded cost measurement. Five-start unrestricted refitting in the extended domain and per-point independent/reverse searches were not completed under WP1. The score is not a validated profile minimum; no scientific delta relative to the old production baseline, mechanism-necessity claim, confidence interval or p value is reported.

At identical resources and search budgets, 100 equal-cost optimizations would allocate about {100 * elapsed * 8 / 3600:.2f} core-hours, plus unrestricted baselines and verification/restarts. This conditional extrapolation is not approval of a batch and does not predict stronger searches or other targets exactly.

## Provenance correction

The previous WP0 package comparison mistakenly selected lower-priority installed Rcpp/jsonlite copies. Original sessionInfo records loaded Rcpp 1.1.1-1.1 and jsonlite 2.0.0; all seven selected packages match the archived HPC snapshot. WP0 metadata and its report were corrected, with package_resolution_correction.json retained. Numerical replay and all 651 frozen inputs were unchanged.

The historical launch SIF digest was not saved, so historical image identity remains unproven. The replay/benchmark image, entry point, package versions and source are explicitly fingerprinted and numerically verified.

## Artifacts

- manifests/profile_domain.json, parameter_status_*.tsv and mapping_counts.tsv: domains and coordinate activity.
- tests/validation_status.json, test_results.tsv and validation logs: final and debugging evidence.
- benchmark/attempt_01/run_manifest.*, parameter_map.rds and population/projection files: complete starts and fingerprints.
- fit_result.rds, verified_components.rds, parameter_result.tsv, evaluations.tsv and validated_result.json: endpoint, predictions, evaluations and acceptance.
- manifests/benchmark_scheduler.tsv: external execution and resource accounting.
- manifests/code_manifest.tsv, scripts/implementation and artifact_manifest.tsv: archived implementation and artifact hashes, excluding runtime caches.

WP1 implemented feasibility, fingerprint, stale/interrupted-result and completed-output guards. Batch scheduling and DE internal-state recovery were not implemented in that phase. Verified endpoints/populations can seed registered new attempts. Production/manuscript assets were preserved; no commit/push occurred under WP1. This English report replaces the earlier Chinese report.
"""
(study / "reports/culture_v2_engineering_report.md").write_text(report)
print(json.dumps(summary, ensure_ascii=False, indent=2))
