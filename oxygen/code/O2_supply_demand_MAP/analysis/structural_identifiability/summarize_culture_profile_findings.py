#!/usr/bin/env python3
"""Describe validated profile traces and audit the recording scope of fluxes."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def positive(value):
    return number(value) is not None and number(value) > 0


def write_table(path, values):
    if not values:
        path.write_text("")
        return
    fields = list(dict.fromkeys(key for row in values for key in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(values)


def table(values, fields):
    head = "| " + " | ".join(fields) + " |\n|" + "|".join("---" for _ in fields) + "|"
    return head + "\n" + "\n".join("| " + " | ".join(str(row.get(key, "")) for key in fields) + " |" for row in values)


def recording_scope(observations):
    grouped = defaultdict(list)
    for row in observations:
        grouped[(row["cohort"], row["segment_id"])].append(row)
    result = {}
    for key, group in grouped.items():
        complete = any(all(positive(row[field]) for field in
            ("recorded_duration_days", "measured_initial_cells", "measured_final_cells")) for row in group)
        result[key] = dict(cohort=key[0], segment_id=key[1],
            recording_scope="counts_and_duration_recorded" if complete else "counts_or_duration_unrecorded",
            observation_ids=";".join(row["observation_id"] for row in group),
            growth_observations=sum(row["growth_scored"] == "TRUE" for row in group),
            chromosome_observations=sum(row["karyotype_scored"] == "TRUE" for row in group),
            flow_observations=sum(row["flow_scored"] == "TRUE" for row in group))
    return result


def search_quality(state, results, tolerance):
    """Describe start support without treating warm repetition as independent recovery."""
    endpoints, summaries = [], []
    for pid, point in sorted(state["points"].items()):
        group = []
        for tid in point["task_ids"]:
            task = state["tasks"][tid]
            if tid not in results:
                continue
            result = results[tid]
            labels = task.get("starts_labels", [])
            pure_independent = len(labels) == 1 and (
                labels[0].startswith("independent_") or task["role"] == "null_independent")
            continuation = len(labels) == 1 and (
                labels[0].startswith("neighbor_") or task["role"] in ("null_reverse", "point_reverse"))
            category = ("independent_anchor" if pure_independent else "neighbor_or_reverse" if continuation
                        else "mixed_warm_and_independent" if task["role"] == "point_strengthen" else "warm_or_alternative")
            group.append(dict(point_id=pid, target=point["target"], value=point["value"], task_id=tid,
                role=task["role"], start_class=category, start_labels=";".join(labels),
                objective=result["objective"], local_convergence=result["optimizer_trace"].get("local_convergence"),
                DE_iterations=result["optimizer_trace"]["iter_completed"],
                budget_limited=result["budget_limited"], winner_phase=result["winner_phase"]))
        group.sort(key=lambda row: row["objective"])
        endpoints.extend(group)
        summary = dict(point_id=pid, target=point["target"], value=point["value"],
            attempted=len(point["task_ids"]), accepted=len(group),
            active=sum(state["tasks"][tid]["status"] not in ("accepted", "failed", "rejected") for tid in point["task_ids"]),
            failed_or_rejected=sum(state["tasks"][tid]["status"] in ("failed", "rejected") for tid in point["task_ids"]),
            agreement_tolerance=tolerance)
        if group:
            best = group[0]["objective"]
            near = [row for row in group if row["objective"] - best <= tolerance]
            independent = [row for row in group if row["start_class"] == "independent_anchor"]
            independent_near = {row["start_labels"] for row in near if row["start_class"] == "independent_anchor"}
            summary.update(winner=group[0]["task_id"], winner_start_class=group[0]["start_class"],
                minimum_objective=best, endpoint_score_spread=group[-1]["objective"] - best,
                second_best_gap=group[1]["objective"] - best if len(group) > 1 else None,
                repeated_search_agreement=len(near) >= 2,
                independent_anchor_searches=len(independent),
                independent_anchor_families_near_minimum=len(independent_near),
                two_independent_anchors_near_minimum=len(independent_near) >= 2,
                independent_best_gap=min(row["objective"] for row in independent) - best if independent else None,
                near_minimum_start_classes=";".join(sorted({row["start_class"] for row in near})),
                winner_local_convergence=group[0]["local_convergence"], winner_budget_limited=group[0]["budget_limited"])
        summaries.append(summary)
    return endpoints, summaries


def flux_summaries(flux, scope, parameters):
    grouped = defaultdict(list)
    labeled = []
    seen = set()
    for row in flux:
        key = (row["task_id"], row["cohort"], row["segment_id"])
        if key in seen:
            raise ValueError("Duplicate segment flux would double count output")
        seen.add(key)
        metadata = scope[(row["cohort"], row["segment_id"])]
        labeled.append(dict(row, **{k: v for k, v in metadata.items() if k not in ("cohort", "segment_id")}))
        grouped[(row["task_id"], row["cohort"], metadata["recording_scope"])].append(row)
    summaries = []
    for (task, cohort, label), group in sorted(grouped.items()):
        wgd = parameters[task]["p_wgd"]
        divisions = sum(float(row["division"]) for row in group)
        ordinary = divisions * (1 - wgd)
        if ordinary <= 0:
            raise ValueError("Nonpositive ordinary-division denominator")
        result = dict(task_id=task, cohort=cohort, recording_scope=label, segments=len(group),
            maternal_divisions=divisions, ordinary_divisions=ordinary)
        for name in ("ordinary_generated_loss", "ordinary_retained_loss", "ordinary_retained_gain",
                     "ordinary_unchanged", "ordinary_biological_nonviability", "numerical_tail_pruning",
                     "ordinary_boundary_loss", "ordinary_boundary_gain", "wgd_retained", "wgd_boundary", "direct_death"):
            total = sum(float(row[name]) for row in group)
            result[name] = total
            if name.startswith("ordinary_"):
                result[name + "_per_ordinary_division"] = total / ordinary
        if result["ordinary_retained_loss_per_ordinary_division"] > 1 + 1e-8:
            raise ValueError("Viable-loss daughter rate exceeds its ordinary-division denominator")
        summaries.append(result)
    return labeled, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    parser.add_argument("study", type=Path)
    parser.add_argument("artifact_root", type=Path, nargs="?")
    args = parser.parse_args()
    project, study = args.project.resolve(), args.study.resolve()
    artifact = args.artifact_root.resolve() if args.artifact_root else study
    points = rows(artifact / "tables/profile_points_validated.tsv")
    endpoints = rows(artifact / "tables/profile_endpoints.tsv")
    nuisance = rows(artifact / "tables/nuisance_paths.tsv")
    functions = rows(artifact / "tables/state_function_values.tsv")
    status = json.loads((artifact / "reports/diagnostic_artifacts_status.json").read_text())
    report_state = json.loads((artifact / "manifests/report_state_snapshot.json").read_text())
    accepted_results = {tid: json.loads((study / "fits" / tid / "validated_result.json").read_text())
                        for tid, task in report_state["tasks"].items() if task["status"] == "accepted"}
    config = json.loads((study / "study.json").read_text())
    search_endpoints, search_points = search_quality(report_state, accepted_results, config["optimization_agreement_tolerance"])
    write_table(artifact / "tables/optimization_endpoint_quality.tsv", search_endpoints)
    write_table(artifact / "tables/optimization_point_quality.tsv", search_points)
    baseline_id = status["baseline_task_id"]
    baseline = json.loads((study / "fits" / baseline_id / "validated_result.json").read_text())
    exact = json.loads((study / "fits" / baseline_id / "artifact_validation_v05.json").read_text())
    parameters = {baseline_id: {k: float(v) for k, v in exact["numeric_params_exact"].items()}}
    paths = defaultdict(dict)
    inactive = defaultdict(list)
    for row in nuisance:
        paths[row["point_id"]][row["parameter"]] = float(row["estimate"])
        if row["data_inactive_at_endpoint"] == "TRUE":
            inactive[row["point_id"]].append(row["parameter"])
        parameters.setdefault(row["winning_task_id"], {})[row["parameter"]] = float(row["estimate"])
    function_rows = {(row["point_id"], int(row["N"]), float(row["oxygen_pct"]), int(row["error_size"])): row for row in functions}
    endpoint_lookup = {(row["point_id"], row["task_id"]): row for row in endpoints}
    traces = []
    for point in points:
        if number(point["objective"]) is None:
            continue
        pid = point["point_id"]
        result = dict(point)
        result.update(paths[pid])
        result["data_inactive_parameters"] = ";".join(sorted(inactive[pid]))
        source = endpoint_lookup[(pid, point["winner"])]
        for modality, component in [("growth", "growth"), ("karyotype", "karyotype"), ("flow", "flow")]:
            result[modality + "_component_delta"] = float(source[modality + "_contribution"]) + float(baseline["component_means"][component])
        if abs(sum(result[k + "_component_delta"] for k in ["growth", "karyotype", "flow"]) - float(point["delta"])) > 1e-8:
            raise ValueError("Component traces do not reconstruct the total objective difference")
        for N in [44, 88]:
            for oxygen, suffix in [(0., "O0"), (20.5, "O20p5")]:
                row = function_rows[(pid, N, oxygen, 1)]
                for name in ["p_per_chromosome", "division_rate_per_day", "death_rate_per_day"]:
                    result[name + "_N" + str(N) + "_" + suffix] = float(row[name])
                probability, survival = float(row["p_per_chromosome"]), float(row["survival"])
                z = 1 - probability + probability * survival
                result["conditional_error_probability_N" + str(N) + "_" + suffix] = probability * survival / z if z > 0 else None
                result["biological_daughter_yield_N" + str(N) + "_" + suffix] = 2 * z**N
            result["survival_N" + str(N) + "_m1"] = float(function_rows[(pid, N, 0., 1)]["survival"])
        result["survival_N88_minus_N44_m1"] = result["survival_N88_m1"] - result["survival_N44_m1"]
        result["mu_hp_over_k_o_mis"] = result["mu_hp"] / result["k_o_mis"] if result["p_misseg"] > 0 else None
        traces.append(result)
    write_table(artifact / "tables/compensation_and_component_traces.tsv", traces)
    observations_path = project / "oxygen/results/structural_identifiability/20261109_seed144_wp0/manifests/observation_manifest.tsv"
    scope = recording_scope(rows(observations_path))
    write_table(artifact / "tables/segment_flux_recording_scope.tsv", list(scope.values()))
    labeled, totals = flux_summaries(rows(artifact / "tables/passage_flux_summaries.tsv"), scope, parameters)
    write_table(artifact / "tables/passage_flux_scope_audit.tsv", labeled)
    write_table(artifact / "tables/flux_summary_by_recording_scope.tsv", totals)
    n_recorded = sum(row["recording_scope"] == "counts_and_duration_recorded" for row in scope.values())
    targets = []
    for target in ["p_misseg", "buffer_beta"]:
        all_points = [row for row in points if row["target"] == target]
        valid = [row for row in traces if row["target"] == target]
        null = next((row for row in valid if float(row["value"]) == 0), None)
        targets.append(dict(target=target, scanned_points=len(all_points), verified_minima=len(valid),
            points_with_start_agreement=sum(row["agreement"] == "TRUE" for row in valid),
            null_objective=null["objective"] if null else "not yet verified",
            null_delta=null["delta"] if null else "not yet verified", null_flags=null["flags"] if null else "unresolved"))
    baseline_totals = [row for row in totals if row["task_id"] == baseline_id]
    unrecorded = sum(row["maternal_divisions"] for row in baseline_totals if row["recording_scope"] == "counts_or_duration_unrecorded")
    total = sum(row["maternal_divisions"] for row in baseline_totals)
    selected = [dict(task_id=row["task_id"], cohort=row["cohort"], segments=row["segments"],
        viable_loss_per_ordinary_division=row["ordinary_retained_loss_per_ordinary_division"],
        viable_gain_per_ordinary_division=row["ordinary_retained_gain_per_ordinary_division"],
        biological_nonviability_per_ordinary_division=row["ordinary_biological_nonviability_per_ordinary_division"])
        for row in totals if row["recording_scope"] == "counts_and_duration_recorded"]
    numerical = [dict(task_id=row["task_id"], cohort=row["cohort"],
        viable_tail_pruning_per_ordinary_division=row["numerical_tail_pruning"] / row["ordinary_divisions"],
        ordinary_grid_loss_per_ordinary_division=(row["ordinary_boundary_loss"] + row["ordinary_boundary_gain"]) / row["ordinary_divisions"])
        for row in totals if row["recording_scope"] == "counts_and_duration_recorded"]
    narrative = """# Culture pilot: evidence and interpretation checks

Generated: {utc}. Execution snapshot: **{execution}**. Unrestricted reference: {baseline_id}, Q={q:.15g}. This is {kind}; optimizer agreement is separate from payload and score validation.

## Current profile evidence

{target_table}

No confidence threshold, p value or formal confidence interval is attached to these score differences. A verified finite candidate does not establish a globally optimized profile minimum. Read missing points, start-agreement flags, local termination and boundaries before drawing a mechanism conclusion.

## Optimization support by start source

{search_table}

Repeated-search agreement means that at least two accepted searches are within the declared 0.001 objective tolerance of the point minimum. This is an optimization diagnostic, not a confidence cutoff. The separate independent-anchor columns count distinct, independently generated feasible anchor vectors near that minimum. Neighbor continuation, production/alternative starts and mixed best-plus-independent strengthening are reported separately; mixed populations cannot identify the winning basin's origin. Agreement between warm searches alone is not independent recovery or proof of global optimality. Reusing an independent anchor across target values also does not make the resulting curve points independent observations.

The historical null_reverse role is a stronger repeat from the nearest positive value to zero, in the same direction as the outward primary scan. Actual opposite-direction checks use the point_reverse role, starting at the best verified null and returning to the nearest positive value for each target. These two searches consume the existing shared reserve; they do not constitute a complete backward sweep of every grid point.

The compensation table exports all fitted nuisance coordinates plus matched-state error, division/death and survival functions. A data-inactive k_o_mis at zero induction is not an estimated saturation scale; its derived death-scale ratio is omitted there. Functions are traces from parameter-constrained candidates, not separately constrained function profiles or confidence bands. The completed Sobol reference supplies hypotheses about links; compensation directions must come from these refits.

The data_inactive_parameters column identifies canonical or start-dependent values that must not be interpreted as fitted compensation. In particular, buffer_n_exp is not estimated at buffer_beta=0.

The growth, chromosome and flow component deltas reconstruct the complete data-score difference. They remain component traces along full-data fits. They do not replace independent stream-deletion profiles or identify a stream's unique information contribution.

## Natural offspring coordinates for compensation checks

For a fixed parent chromosome count N and oxygen, write the per-copy error probability as p and per-error survival as s, so survival after m errors is s^m. Before numerical tail truncation and grid losses, the ordinary-daughter weights have the exact binomial identity: z = 1 - p + p*s, q = p*s/z, and Y = 2*z^N. Here q is the error probability in the error-size distribution conditional on biological daughter viability; Y is the expected number of biologically viable daughters per ordinary division. The conditional q is undefined if z is zero. The identity was checked against full PMF sums for 80 cases, including boundary values.

Compare q and Y alongside raw p, s, division and death rates in the fitted traces. These are diagnostic coordinates, not new parameter profiles or a structural-identifiability proof. The simulator continues to use unnormalized expected daughter numbers; it is not replaced by a normalized Binomial(N,q) kernel. Truncation, grid loss, WGD and maternal removal remain separate in the actual generator and flux ledger.

Knowing both q and Y at a known N would in principle recover p = 1 - z*(1-q), with z = (Y/2)^(1/N), and s = q*z/p when p is positive. Observed culture distributions and net growth do not directly provide those two quantities independently of division, death, state composition and protocol. Thus the identity motivates testing compensation and function stability, without claiming the full model is structurally nonidentifiable.

## Recording scope of daughter-output summaries

The frozen observation manifest assigns {segments} modeled segments to {recorded} segments with a recorded positive duration and positive initial/final counts, and {unrecorded} segments lacking that complete recording. Some of the latter still carry chromosome observations. All segments and observations remain in the production fitting protocol.

In the unrestricted reference, segments without complete counts/duration contribute {share:.6%} of the modeled maternal-division total when all segments are summed. Such totals are dominated by assumed prehistory, default abundance and fallback duration, rather than being a direct experimental cumulative event measurement. Headline cross-fit flux comparisons therefore use the **same frozen set of count-and-duration-recorded segments**, while retaining both scopes in the complete ledger. This postprocessing scope does not change simulations, propagation, fitting or scores.

Even within recorded segments, the production target-count rule selects modeled end days and can differ from nominal passage duration. Outputs remain conditional model predictions. Recording-scope stratification does not remove the dependence of inherited starting compositions on the full, partly assumed history. Sums over culture segments include branching histories and are not independent replicate observations or counts of unique descendants. The conditional ordinary-division denominator is integrated maternal divisions times (1 - p_wgd), not all divisions or chromosome copies.

{flux_table}

The event ledger separately retains biological nonviability, viable-tail pruning, ordinary/WGD grid losses, the one-retained-lineage WGD branch, direct death and numerical nonnegative corrections. These channels must not be conflated with failed biological daughters or renormalized away. Full-grid kernel and dense Euler balance checks are prerequisites for the displayed outputs.

{numerical_table}

The finite chromosome grid and its drop boundary are part of the production model used for these fits. A large grid-loss rate makes a mechanism comparison conditional on that boundary rule; it is not evidence of per-error biological nonviability. The abstract Y function is defined before these losses. A different grid or boundary convention would require a separately declared comparison and its own unrestricted refit, which is not launched in this pilot.

## Follow-up decisions

After the two-target pilot and its search diagnostics are complete, use actual compensation and prediction differences to select a small measurement-stream deletion comparison. Each mask needs its own unrestricted baseline and constrained nuisance refits. Growth-count conditioning remains after score deletion, so that audit is not automatically an independent held-out test.

If error generation and survival remain interchangeable in good fits, compare candidate finite-time measurements that directly distinguish chromosome errors from daughter fate at matched parental chromosome count. If assumed prehistory drives a conclusion, propose a separately declared protocol/initial-state sensitivity with its own baseline. No such extension or new experiment is launched here.

The future tumor plan includes rho_2N and initial-cell-number confounding; tumor fitting, COMBOS, broader targets and statistical calibration remain outside this culture stage.

## Files

- tables/compensation_and_component_traces.tsv: complete parameter paths, matched-state functions and component deltas.
- tables/optimization_point_quality.tsv, optimization_endpoint_quality.tsv: all search endpoints, actual start-source labels, independent-anchor support, local convergence and compute-limit flags.
- figures/effective_daughter_functions.png: conditional error probabilities and biological ordinary-daughter yields before numerical losses.
- tables/segment_flux_recording_scope.tsv: target-independent recording-scope assignment from frozen observations.
- tables/passage_flux_scope_audit.tsv: complete per-segment fluxes with explicit recording scope.
- tables/flux_summary_by_recording_scope.tsv: both scopes and cohort-specific ordinary-division denominators.
""".format(utc=datetime.now(timezone.utc).isoformat(), execution=status.get("execution_status", "unknown"),
        baseline_id=baseline_id, q=baseline["objective"], kind="a bounded-stage record" if status.get("execution_status") == "pilot_fitting_finished" else "an interim diagnostic record",
        target_table=table(targets, ["target", "scanned_points", "verified_minima", "points_with_start_agreement", "null_objective", "null_delta", "null_flags"]),
        search_table=table(search_points, ["point_id", "accepted", "active", "repeated_search_agreement", "independent_anchor_families_near_minimum", "independent_best_gap", "winner_local_convergence"]),
        segments=len(scope), recorded=n_recorded, unrecorded=len(scope)-n_recorded, share=unrecorded/total,
        flux_table=table(selected, ["task_id", "cohort", "segments", "viable_loss_per_ordinary_division", "viable_gain_per_ordinary_division", "biological_nonviability_per_ordinary_division"]),
        numerical_table=table(numerical, ["task_id", "cohort", "viable_tail_pruning_per_ordinary_division", "ordinary_grid_loss_per_ordinary_division"]))
    (artifact / "reports/scientific_findings.md").write_text(narrative)
    source_files = [observations_path, artifact / "manifests/report_state_snapshot.json",
                    artifact / "reports/diagnostic_artifacts_status.json",
                    artifact / "tables/profile_points_validated.tsv", artifact / "tables/passage_flux_summaries.tsv",
                    artifact / "tables/nuisance_paths.tsv", artifact / "tables/state_function_values.tsv",
                    artifact / "tables/profile_endpoints.tsv", study / "fits" / baseline_id / "validated_result.json",
                    study / "fits" / baseline_id / "artifact_validation_v05.json", study / "study.json", Path(__file__)]
    source_files += [study / "fits" / tid / "validated_result.json" for tid in accepted_results if tid != baseline_id]
    (artifact / "manifests/findings_source_manifest.json").write_text(json.dumps(dict(
        generated_utc=datetime.now(timezone.utc).isoformat(), model_calls=0,
        sources=[dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in source_files]), indent=2) + "\n")
    print(json.dumps(dict(verified_profile_minima=len(traces), recorded_segments=n_recorded,
                         unrecorded_segments=len(scope)-n_recorded, baseline_unrecorded_division_share=unrecorded/total), indent=2))


if __name__ == "__main__":
    main()
