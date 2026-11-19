"""Blind stationary-point checks when the standard branch repair does not agree."""
import argparse
import concurrent.futures as cf
from collections import defaultdict

import numpy as np

from fit_stage import model
from inference import likelihood_ratio, pool
from practical import Configuration, HERE, read, result_row, row_x, solve_constrained, write
from recovery import CONVERGE_LR, KKT_TOL, SUPPLY_TOL


def task(item):
    row, target, alternatives = item
    k = model()
    c = Configuration(k)
    anchor = row_x(c, row)
    rng = np.random.default_rng(9829 + int(row["dataset"]) + 1000 * row["scenario"].startswith("B"))
    baseline = float(row["loss"])
    incumbent_supply = float(row["min_aggregate_donor_margin_fraction"])
    starts = [
        np.clip(anchor + rng.normal(0, .00003, len(anchor)), 0, 1),
        np.clip(anchor + rng.normal(0, .0005, len(anchor)), 0, 1),
        np.clip(anchor + rng.normal(0, .004, len(anchor)), 0, 1),
        *[row_x(c, r) for r in alternatives],
    ]
    results = []
    stationary = []
    gap = float("inf")
    for attempt, start in enumerate(starts, 1):
        fit = solve_constrained(k, c, target, start, limit=82, best_iterate_tolerance=1e-6)
        if fit["supply"] >= -SUPPLY_TOL:
            results.append(fit)
            if fit["gradient"] <= KKT_TOL:
                stationary.append(fit)
        if len(stationary) < 2:
            continue
        ordered = sorted(stationary, key=lambda r: r["loss"])
        best = ordered[0]
        floor = min(r["loss"] for r in results)
        if incumbent_supply >= -SUPPLY_TOL:
            floor = min(floor, baseline)
        gap = likelihood_ratio(ordered[1]["loss"], best["loss"], k.n)
        if gap <= CONVERGE_LR and likelihood_ratio(best["loss"], floor, k.n) <= .001:
            return result_row(c, best, scenario=row["scenario"], dataset=int(row["dataset"]),
                              qualified=True, independent_LR_gap=gap,
                              attempts=int(row["attempts"]) + attempt)
    if not results:
        return None
    best = min(results, key=lambda r: r["loss"])
    return result_row(c, best, scenario=row["scenario"], dataset=int(row["dataset"]),
                      qualified=False, independent_LR_gap=gap,
                      attempts=int(row["attempts"]) + len(starts))


def run(scenario, workers, dataset_max, after_standard):
    path = HERE / f"strict_recovery_repair_{scenario}.tsv"
    rows = read(path) if path.exists() else []
    done = {int(r["dataset"]) for r in rows}
    targets = defaultdict(list)
    for r in read(HERE / "noisy_slope_inputs.tsv"):
        targets[r["scenario"], int(r["dataset"])].append(float(r["ols_slope"]))
    choices = read(HERE / "qualified_recovery.tsv")
    stable = [r for r in choices if r["qualified"] == "True"]
    standard_path = HERE / f"pilot_near_branch_confirmation_expansion{'_b' if scenario == 'b' else ''}.tsv"
    standard_done = {int(r["dataset"]) for r in read(standard_path)} if standard_path.exists() else set()
    k = model()
    jobs = []
    for r in choices:
        key = r["scenario"], int(r["dataset"])
        if not key[0].lower().startswith(scenario) or r["qualified"] == "True" or key[1] in done:
            continue
        if key[1] > dataset_max:
            continue
        if after_standard and key[1] > 20 and key[1] not in standard_done:
            continue
        y = np.array(targets[key])
        neighbours = sorted(stable, key=lambda p: np.linalg.norm(
            (y - targets[p["scenario"], int(p["dataset"])]) / k.scale))
        jobs.append((r, y, neighbours[:2]))
    with pool(workers) as workers_pool:
        for f in cf.as_completed([workers_pool.submit(task, item) for item in jobs]):
            row = f.result()
            if row is not None:
                rows.append(row)
                write(path, rows)
                print(row["scenario"], row["dataset"], row["qualified"], row["independent_LR_gap"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scenario", choices=("a", "b"))
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--dataset-max", type=int, default=100)
    parser.add_argument("--after-standard", action="store_true")
    args = parser.parse_args()
    run(args.scenario, args.workers, args.dataset_max, args.after_standard)
