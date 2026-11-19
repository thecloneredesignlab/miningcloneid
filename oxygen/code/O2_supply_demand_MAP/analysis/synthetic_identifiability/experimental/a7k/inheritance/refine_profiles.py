"""Verify one-sided stationarity at profile crossings and the beta=10 endpoint."""
import argparse
from concurrent.futures import as_completed
from pathlib import Path

from profile_tilts import summary
from tilt_model import HERE, TARGETS, THRESHOLD, case, pool, read, run_fit, stable, write


def grid(mode):
    paths = sorted(HERE.glob(f"grid_{mode}_*.tsv"))
    if len(paths) != 4:
        raise RuntimeError(f"Incomplete profiles for {mode}")
    return [r for path in paths for r in read(path)]


def refine(mode, rows):
    best = read(HERE / f"best_{mode}.tsv")[0]
    k, c, _ = case(mode)
    for name in TARGETS:
        j = c.names.index(name)
        center = float(best[f"x_{name}"])
        subset = [r for r in rows if r["parameter"] == name and
                  float(r["all_constraint_margin"]) >= 0.]
        by_value = {}
        for r in subset:
            by_value.setdefault(float(r["coordinate"]), []).append(r)
        nearest = set()
        for sign in (-1, 1):
            side = [v for v in by_value if sign * (v - center) > 0.]
            for within in (True, False):
                chosen = [v for v in side
                          if (min(float(r["lr"]) for r in by_value[v]) <= THRESHOLD) == within]
                nearest.update(sorted(chosen, key=lambda v: abs(
                    min(float(r["lr"]) for r in by_value[v]) - THRESHOLD))[:1])
        if name == "buffer_beta":
            nearest.add(1.)
        for v in nearest:
            viable = by_value.get(v, [])
            if not viable:
                continue
            winner = min(viable, key=lambda r: float(r["loss"]))
            if not any(stable(r) and float(r["loss"]) <= float(winner["loss"]) + .005
                       for r in viable):
                yield mode, name, v, j, c.anchor(winner), float(best["loss"])


def clean(job):
    mode, name, v, j, start, best = job
    results = []
    for divisor in (30., 14., 55.):
        r = run_fit((mode, "crossing_recheck", start, {j: v}, 150, divisor))
        r.update(parameter=name, coordinate=v, lr=float(r["loss"]) - best)
        results.append(r)
        if stable(r):
            break
        start = case(mode)[1].anchor(r)
    return results


def refine_all(modes, workers):
    raw = {mode: grid(mode) for mode in modes}
    work = [job for mode in modes for job in refine(mode, raw[mode])]
    output = {mode: [] for mode in modes}
    if work:
        with pool(min(workers, len(work))) as ex:
            for f in as_completed(ex.submit(clean, job) for job in work):
                rows = f.result()
                r = rows[-1]
                output[r["mode"]].extend(rows)
                print(r["mode"], r["parameter"], r["coordinate"], r["lr"], r["stationarity"], flush=True)
    for mode in modes:
        combined = raw[mode] + output[mode]
        interval, better = summary(mode, combined)
        if better:
            raise RuntimeError(f"Profile exposed a better free fit for {mode}")
        write(HERE / f"profile_points_{mode}.tsv", combined)
        write(HERE / f"intervals_{mode}.tsv", interval)
        print(mode, [(r["parameter"], r["lower"], r["upper"], r["precision_gate"])
                     for r in interval], flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("modes", nargs="+")
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    refine_all(a.modes, a.workers)
