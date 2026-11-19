"""Cross-check the precision-critical edges from adjacent fitted starts."""
import argparse
import concurrent.futures as cf

from fit_stage import model
from inference import SLACK, likelihood_ratio, pool
from practical import Configuration, HERE, MECHANISMS, read, result_row, row_x, solve_constrained, write


PROBES = {
    "lam_max": (.89999, .925, .987343, .993),
    "buffer_smax": (.89413543, .94, .956245625, .9624962499999999, .968, .972, 1.0),
    "n_O": (.11127884801685528, .11540028683229436, .50, .515625, .54),
}


def fit_point(task):
    name, value, starts, reference_loss = task
    k = model()
    c = Configuration(k)
    parameter = MECHANISMS.index(name)
    results = []
    for start in starts:
        fit = solve_constrained(k, c, k.y, start, limit=92, fixed=(parameter, value))
        results.append(fit)
        if len(results) >= 2:
            a, b = results[-2:]
            if min(a["supply"], b["supply"]) >= -SLACK and max(a["gradient"], b["gradient"]) <= .015:
                if abs(likelihood_ratio(a["loss"], b["loss"], k.n)) <= .02:
                    break
    feasible = [r for r in results if r["supply"] >= -SLACK]
    best = min(feasible or results, key=lambda r: r["loss"])
    return result_row(c, best, case="main", profile_parameter=name, fixed_normalized=value,
                      lr=likelihood_ratio(best["loss"], reference_loss, k.n) if feasible else float("inf"),
                      level=3.841458820694124, feasible=bool(feasible))


def run(probes, output, workers):
    c = Configuration(model())
    reference = read(HERE / "main_constrained_best.tsv")[0]
    rows = read(HERE / "main_profiles_envelope.tsv")
    tasks = []
    for name, values in probes.items():
        valid = [r for r in rows if r["profile_parameter"] == name and r["feasible"] == "True"
                 and float(r["max_projected_gradient"]) <= .035]
        for value in values:
            closest = sorted(valid, key=lambda r: abs(float(r["fixed_normalized"]) - value))
            starts = [row_x(c, r) for r in (closest[0], closest[1], reference)]
            tasks.append((name, value, starts, float(reference["loss"])))
    out = []
    with pool(workers) as work:
        for f in cf.as_completed([work.submit(fit_point, task) for task in tasks]):
            row = f.result()
            out.append(row)
            write(HERE / output, out)
            print(row["profile_parameter"], row["fixed_normalized"], row["lr"],
                  row["max_projected_gradient"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--extension", action="store_true")
    parser.add_argument("--final-crossings", action="store_true")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if args.final_crossings:
        run({"lam_max": (.9775, .979, .981), "buffer_smax": (.888, .890, .892, .8941)},
            "main_profiles_final_crossings.tsv", args.workers)
    elif args.extension:
        run({"lam_max": (.978, .9825, .9855), "buffer_smax": (.86, .88, .90, .92)},
            "main_profiles_extension.tsv", args.workers)
    else:
        run(PROBES, "main_profiles_edge_refit.tsv", args.workers)
