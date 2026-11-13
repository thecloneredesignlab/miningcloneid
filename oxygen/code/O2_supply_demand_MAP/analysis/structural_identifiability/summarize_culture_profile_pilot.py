#!/usr/bin/env python3
"""English status/findings report from immutable, accepted pilot artifacts."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from datetime import datetime, timezone

project = Path(sys.argv[1]).resolve()
study = Path(sys.argv[2]).resolve()
state = json.loads((study / "state.json").read_text())
config = json.loads((study / "study.json").read_text())
controller_registry = json.loads((study / "manifests/orchestrator_versions.json").read_text())
controller_revision = controller_registry["active"]["path"]


def load(path):
    return json.loads(path.read_text())


def rows(path):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def rss_mib(text):
    if not text:
        return None
    scale = {"K": 1 / 1024, "M": 1, "G": 1024, "T": 1024**2}
    return float(text[:-1]) * scale[text[-1]] if text[-1] in scale else float(text) / 1024**2


def table(values, fields):
    head = "| " + " | ".join(fields) + " |\n|" + "|".join("---" for _ in fields) + "|"
    return head + "\n" + "\n".join("| " + " | ".join(str(x.get(f, "")) for f in fields) + " |" for x in values)


accepted = {}
index = []
for tid, task in sorted(state["tasks"].items()):
    directory = study / "fits" / tid
    result = load(directory / "validated_result.json") if task["status"] == "accepted" else None
    record = {key: task.get(key) for key in ["task_id", "job_id", "stage", "role", "primary_family", "status", "reserve", "reason"]}
    if result:
        assert result["accepted"] and result["launch_completed"]
        guard_file = result.get("artifact_validation_file", "artifact_validation.json")
        guard = load(directory / guard_file)
        assert guard["pass"] and guard["target_constraint_tolerance"] == 0
        assert abs(guard["objective"] - result["objective"]) <= 1e-8
        if "artifact_validation_sha256" in result:
            assert hashlib.sha256((directory / guard_file).read_bytes()).hexdigest() == result["artifact_validation_sha256"]
        accepted[tid] = result
        record.update(objective=result["objective"], native_objective=result["native_objective"],
                      DE_iterations=result["optimizer_trace"]["iter_completed"], DE_cap=task["DE_itermax"],
                      local_convergence=result["optimizer_trace"].get("local_convergence"),
                      evaluations=result["total_evaluator_invocations"], failures=result["failed_evaluator_invocations"],
                      allocated_core_hours=result["allocated_core_hours"])
    accounting = load(directory / "scheduler.json") if (directory / "scheduler.json").exists() else task.get("accounting", [])
    parent = next((x for x in accounting if x["JobID"] == task["job_id"]), None)
    if parent:
        record["elapsed_seconds"] = int(parent["ElapsedRaw"])
        record["allocated_core_hours"] = int(parent["AllocCPUS"]) * int(parent["ElapsedRaw"]) / 3600
        record["scheduler_state"] = parent["State"]
    memory = [rss_mib(x.get("MaxRSS", "")) for x in accounting]
    memory = [x for x in memory if x is not None]
    record["MaxRSS_MiB"] = max(memory) if memory else None
    index.append(record)

fields = ["task_id", "job_id", "stage", "role", "primary_family", "status", "reserve", "scheduler_state", "reason", "objective",
          "native_objective", "DE_iterations", "DE_cap", "local_convergence", "evaluations", "failures", "elapsed_seconds", "allocated_core_hours", "MaxRSS_MiB"]
with (study / "tables/execution_ledger.tsv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader(); writer.writerows(index)
baselines = [x for x in index if x["stage"] == "baseline" and x["status"] == "accepted"]
best = min(baselines, key=lambda x: x["objective"]) if baselines else None
points = rows(study / "tables/profile_points_validated.tsv")
diagnostics = load(study / "reports/diagnostic_artifacts_status.json") if (study / "reports/diagnostic_artifacts_status.json").exists() else None
resources = load(study / "reports/resource_summary.json") if (study / "reports/resource_summary.json").exists() else None
resource_description = ("Receipt-based scheduler accounting at {} covers {} parent jobs, with {:.6f} allocated core-hours across fitting, validation/reporting and controllers. Missing job records: {}; unfinished jobs at that snapshot: {}. Parent allocations are counted once; step CPU/RSS records are retained separately.".format(
    resources["generated_utc"], resources["accounted_parent_jobs"], resources["total_allocated_core_hours"],
    len(resources["missing_job_ids"]), len(resources["unfinished_job_ids"])) if resources else
    "The final receipt-based scheduler summary has not yet been generated.")
finished = state["status"] in ["pilot_fitting_finished", "baseline_optimizer_unresolved", "baseline_repair_budget_exhausted"]
gate = state.get("baseline_gate", {})
summaries = []
for target in ["p_misseg", "buffer_beta"]:
    group = [p for p in points if p["target"] == target]
    valid = [p for p in group if p.get("objective", "") not in ("", "NA")]
    null = next((p for p in valid if float(p["value"]) == 0), None)
    summaries.append(dict(target=target, evaluated_points=len(group), verified_points=len(valid),
        null_Q=float(null["objective"]) if null else "unresolved",
        null_delta=float(null["delta"]) if null else "unresolved",
        null_flags=null.get("flags", "") if null else "no verified null result",
        points_without_start_agreement=sum(p.get("agreement") != "TRUE" for p in group)))

completed_core_hours = sum(x.get("allocated_core_hours", 0) for x in index)
updated = datetime.now(timezone.utc).isoformat()
report = """# Culture v2 profile-objective pilot

Generated: {updated}

Execution status: **{status}**. This report is {report_kind}. Baseline search gates, full-precision target checks and optimization diagnostics remain separate from process completion.

## Scope and execution contract

The study uses the verified iteration5/seed144 production-v2 model at commit 1a5d861c9c4d2af72ea0d07def0023831cc531e4, frozen data and the fingerprinted R 4.4.2 SIF. Both initial-ploidy cohorts and complete lineage histories are retained. Scoring uses 114 growth records, 12 chromosome samples/220 cells and 20 flow distributions, with production averaging and modality weights. Every retained observation must have a finite prediction. The assigned-oxygen and passage-selection/propagation rules are preserved.

The extended domain admits exact zero for p_misseg and buffer_beta. Other production bounds remain assumptions of this study. No prior, soft-coupling or phenotype regularization is added; production simulation/protocol failure barriers remain. Model cfg fields that mention priors do not create a penalty contribution.

All jobs use xxlarge. There is no application concurrency limit, array throttle or project running-core cap. Slurm manages allocation. Per-fit resources are eight cores and 16 GB; primary searches allow 500 DE generations/50 local iterations and reserve searches allow 1000/100, with 8/12-hour limits. Native early termination can use substantially fewer generations; the ledger reports actual work. The approved total is at most 125 fits, including 20 reserves, and 8,700 allocated core-hours including non-fit work.

## Validation and initialization

The frozen implementation passed 28 R test cases/1,402 assertions. Additional actual-model worker tests covered unrestricted, no-induction and flat-survival maps and saved full candidates across serialized calls. Daughter accounting passed 16 kernel cases and four dense trajectory replays. Separate Python tests cover the scientific dependency graph, minimum reduction, independent-start gates, budget handling, JSON transport and concurrent metadata updates. Active orchestration: {controller_revision}. Its nine policy tests and seven dependency/concurrency tests passed locally and on HPC, including five-search refinement allocation and null-to-positive reverse-start provenance. A renderer-only fixture exercised both target export/plot branches and all 122 baseline segments without new fitting; fixture outputs are not profile evidence.

Production seed144, seed371 and seed10 were replayed before reuse. Parameter-space diversity is an initialization heuristic, not proof of distinct basins. Independent anchors were screened only for protocol feasibility. An initial 1,000 full-box draws found no feasible anchor; the recorded mixed proposal adds a growth-favoring random anchor distribution. Only initialization proposals changed: the remaining DE members still cover the complete declared transformed domain. A separately verified seed325 archive endpoint supplements the low-score region when used.

Each fit checkpoints improved complete parameter vectors per process. The saved minima and native endpoint are independently replayed; the lowest verified feasible result is selected. Exact target agreement is checked in full-precision RDS objects. JSON rounding is allowed only in positive-target summaries; the runtime target residual tolerance remains zero, including strict zero at nulls. Actual/requested native DE populations are checked, and all inputs/results have fingerprints.

Warm strengthening from the best known point cannot itself establish independent-start agreement. The baseline gate separately requires low-score support from different independently initialized searches and a stronger-search check. Numerical agreement is not proof of global optimality. Abnormal local termination, early DE stopping and exhausted budgets remain visible.

## Unrestricted baseline candidates

Best externally validated completed candidate: {best_description}. Baseline gate passed: **{gate_pass}**. Gate diagnostic: `{gate_reason}`.

{baseline_table}

## Profiles

{profile_table}

Each initial target grid contains exact zero, the refitted optimum at grid creation and eight positive log-spaced probes. Each point has five nuisance searches using unrestricted, neighbor, alternative and independent starts; duplicated projected starts are replaced. These are distinct optimization starts, not five independently generated anchors. Repeated-search agreement and independent-anchor recovery are reported separately in optimization_point_quality.tsv. Positive targets have 19 active coordinates; exact nulls have 18. All other mechanisms, initial-distribution parameters and fitted error scales remain adjustable.

Profiles use verified minima, not endpoint means. Better restricted candidates trigger unrestricted repair and updates to every delta; negative deltas are never hidden. Missing/unresolved points remain explicit. The legacy null_reverse label denotes a stronger positive-to-zero repeat. Two actual opposite-direction checks return from the best null to the nearest positive point, before remaining quality searches or optional curvature refinement. These use the same shared reserve and do not form a complete backward sweep. Component contributions along a full-data optimum are component traces, not separately refitted measurement-stream profiles.

## Boundaries, mechanisms and biological outputs

Bound hits are reported without automatically extending the domain. Physical ceilings and production search limits have different meanings. Any expanded-bound comparison needs a new study and unrestricted baseline. Coordinates inactive at the selected endpoint are distinguished from quantities freely varied during optimization.

The no-induction restriction preserves baseline chromosome errors, WGD, growth/death and survival. Flat survival preserves fitted buffer_smax and error-size dependence; it removes the parental-ploidy dependence at fixed error size. State-dependent division/death and the number of chromosome-error opportunities remain.

The diagnostic artifacts include full chromosome/flow distributions, growth predictions, compensation trajectories and selected-passage-day changes. Function values at matched model states are traces of parameter-profile fits, not derived-quantity profiles or confidence bands.

An exact ordinary-daughter identity before numerical losses separates z=1-p+p*s, conditional error probability q=p*s/z and biological viable-daughter yield Y=2*z^N at fixed parental chromosome count. Full-PMF checks passed 80 cases. The effective-daughter plots track these combinations alongside division/death functions; no normalized kernel replaces the production expected-daughter accounting, and this algebra is not a structural-identifiability result for the complete model.

The completed global Sobol N131072 reference supplies diagnostic hypotheses only. Its 14-factor fixed-oxygen eigenvalue/eigenvector study has different bounds, outputs and protocols. The bridge audit checks all 47,838 index cells and reconstructs the 20,100 focal cells available in its five-repeat table without model calls. Track mu_hp, gamma_mu, k_o_mis and p_mis_base along the induction profile, and the three survival parameters jointly along the buffering profile. State-function exports expose the death saturation, induced increment and survival penalty. Old S2 has no compensation direction; no old bound, sampling prior, endpoint or sensitivity threshold is imported into this study. See reports/sobol_reference_integration.md and sobol_bridge/ for source hashes and the full compatibility table.

Representative daughter fluxes use the production kernel and an Euler-step dense replay with original selected days fixed. Accounting separates biological failed daughters, numerical tail pruning, ordinary/WGD grid losses, WGD's one retained doubled lineage, maternal removal and direct death. Saved daily states and discrete event balances are checked. Integrals count expected daughter production before passage dilution, not unique descendants or ancestry. Loss output must be assessed with gain output and differential expansion.

Per-chromosome summaries explicitly distinguish cell weights, chromosome weights, division weights and division-by-chromosome weights. Expected errors per ordinary division use sum(x*lambda*N*p)/sum(x*lambda), not a per-chromosome probability. The common crowding factor and constant ordinary-branch probability cancel in conditional division averages.

Cumulative flux comparisons separate the 114 segments with recorded duration and initial/final counts from eight segments without complete count/duration records. The latter remain in the fitting history and some carry chromosome observations; their assumed prehistory can dominate absolute event totals. Complete ledgers retain both scopes. Headline flux comparisons use the same frozen count-and-duration-recorded subset for every fit, with ordinary-division denominators and no claim of independent replicate observations or unique descendants. See scientific_findings.md for the scope audit and numerical comparisons.

## Resource accounting and limitations

Recorded allocated core-hours of completed fitting jobs: {completed_core_hours:.4f}. Maximum fitting allocation committed by submitted tasks: {reserved_core_hours:.2f}; this is a reservation ceiling, not consumed time. Active jobs and non-fit jobs require the final scheduler ledger for total accounting. All failed/rejected attempts are retained.

{resource_description}

This is an objective diagnostic with inference_calibration=absent. No 1.92/3.84 cutoff, likelihood-ratio p value or formal confidence interval is used. No measurement-stream deletion, COMBOS, tumor/joint fitting or manuscript regeneration belongs to this stage. Growth-count conditioning remains part of the model, so future score-deletion audits are not automatically independent predictive validation.

The future tumor plan now includes rho_2N identifiability and an initial-cell-number confounding audit, with separate data-only and production-penalized baselines. That plan inclusion does not authorize tumor fitting in the current stage; see manifests/future_tumor_plan_scope.json and the two main implementation plans.

The exact historical launch-image digest remains unproven; the current replay/analysis image and source are fingerprinted and numerically checked. Existing WP0/WP1 reports and report generators were replaced with English versions without changing numerical evidence or frozen production inputs.

## Artifacts

- `tables/execution_ledger.tsv`: all fitting attempts and scheduler/resource/search diagnostics.
- `reports/accepted_artifact_audit.json`: accepted payload byte hashes and full-precision guard checks before rendering.
- `tables/scheduler_accounting.tsv`, `reports/resource_summary.json`: receipt-selected jobs and actual parent allocations, including failed and canceled attempts.
- `tables/profile_points_validated.tsv`, `profile_endpoints.tsv`, `nuisance_paths.tsv`: minima, matched baseline and compensation.
- `tables/growth_predictions.tsv`, `chromosome_distributions.tsv`, `flow_distributions.tsv`: complete scored predictions.
- `tables/state_function_values.tsv`, `state_daughter_rates.tsv`, `daughter_fluxes.tsv`, `passage_flux_summaries.tsv`: functions and accounted output.
- `tables/effective_response_functions.tsv`, `survival_functions.tsv`: matched-state response curves and full chromosome-grid survival functions for representative candidates, including the unrestricted baseline.
- `reports/scientific_findings.md`, `tables/compensation_and_component_traces.tsv`: profile evidence, complete parameter/function paths and explicitly labeled component traces.
- `tables/passage_flux_scope_audit.tsv`, `flux_summary_by_recording_scope.tsv`: complete and scope-separated daughter accounting without assumed prehistory dominating headline counts.
- `figures/`: English profile, compensation and prediction panels when available.
- `fits/<task>/`: native and selected full fits, actual starts, checkpoints, independent checks and scheduler evidence.
- `manifests/orchestrator_versions.json`: active controller/checker hashes; the frozen fitting implementation remains unchanged.
- `reports/sobol_reference_integration.md`, `sobol_bridge/`: completed Sobol reference audit, parameter domains and diagnostic hypotheses.

{completion_note}
""".format(updated=updated, status=state["status"], controller_revision=controller_revision, report_kind="a final bounded-stage record" if finished else "an interim execution record",
    best_description=("{} at Q={:.15g}".format(best["task_id"], best["objective"]) if best else "none yet"),
    gate_pass=gate.get("pass_gate", False), gate_reason=gate.get("reason", "initial searches still running"),
    baseline_table=table(baselines, ["task_id", "primary_family", "objective", "DE_iterations", "local_convergence"]),
    profile_table=table(summaries, ["target", "evaluated_points", "verified_points", "null_Q", "null_delta", "points_without_start_agreement"]),
    completed_core_hours=completed_core_hours, reserved_core_hours=state["fit_core_hours_reserved"],
    resource_description=resource_description,
    completion_note=("The baseline gate was not met; no profile interpretation is supported." if state["status"] == "baseline_optimizer_unresolved" else
                     "Unresolved search or boundary flags limit any interpretation of the recorded curves." if finished else
                     "Computations remain active. The completed-candidate score is provisional until the baseline gate and profile checks finish."))
(study / "reports/culture_profile_pilot_report.md").write_text(report)
print(json.dumps(dict(status=state["status"], accepted_fits=len(accepted), completed_fit_core_hours=completed_core_hours,
                     report=str(study / "reports/culture_profile_pilot_report.md")), indent=2))
