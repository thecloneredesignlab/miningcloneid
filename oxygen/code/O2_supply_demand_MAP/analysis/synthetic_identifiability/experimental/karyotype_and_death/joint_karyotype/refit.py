"""Re-establish the count-only baseline and find joint fits from several basins."""
import argparse
from concurrent.futures import as_completed
from functools import lru_cache

import numpy as np

from joint_model import (BASE, HERE, Config, Experiment, fit_receipt, likelihood_ratio,
                         pool, read, row_x, solve, stable, write)

MODES = {
    "autosomes_count": ("autosomes", 0., False),
    "autosomes_joint": ("autosomes", 0., True),
    "archived_joint": ("archived", 0., True),
    "autosomes_joint_sd1": ("autosomes", 1., True),
}


@lru_cache
def case(mode):
    convention, width, joint = MODES[mode]
    k = Experiment(convention, width)
    return k, Config(k), k.target(kary=k.kary if joint else None)


def vector_bank(c):
    inputs = ["observed_best", "observed_refits", "truths", "profile_direct",
              "profile_refined_tolerance"]
    rows = [r for name in inputs for r in read(BASE/f"{name}.tsv")
            if float(r["donor_margin"]) >= 0.]
    return np.unique([row_x(c,r) for r in rows],axis=0)


def task_fit(job):
    mode, origin, anchor, fixed, limit, divisor = job
    k,c,y = case(mode)
    r = solve(k,c,y,anchor,limit=limit,fixed=fixed,divisor=divisor)
    return fit_receipt(k,c,y,r,mode=mode,origin=origin)


def select_bank(mode):
    k,c,y = case(mode)
    x0 = row_x(c,read(BASE/"observed_best.tsv")[0])
    truth = [row_x(c,r) for r in read(BASE/"truths.tsv")]
    bank = vector_bank(c)
    scores = [k.evaluate_joint(c.unmap(x)[0],y)[2] for x in bank]
    rng = np.random.default_rng(20261119)
    starting = [("stage2_best",x0),("stage2_alt",truth[1])]
    for other in MODES:
        p=HERE/f"best_{other}.tsv"
        if other!=mode and p.exists():
            starting.append((other,row_x(c,read(p)[0])))
    for i in np.argsort(scores):
        if min(np.linalg.norm(bank[i]-x) for _,x in starting) > .075:
            starting.append((f"archived_bank_{i}",bank[i]))
        if len(starting)>=8:
            break
    for j,scale in enumerate((.04,.12,.24)):
        starting.append((f"bounded_jitter_{j}",np.clip(x0+rng.normal(0,scale,len(x0)),0,1)))
    return starting


def refit(modes, workers):
    jobs = [(t,label,x,{},100,32.) for t in modes for label,x in select_bank(t)]
    results = {m:[] for m in modes}
    with pool(workers) as executor:
        for f in as_completed(executor.submit(task_fit, j) for j in jobs):
            r=f.result()
            results[r["mode"]].append(r)
            print(r["mode"],r["origin"],round(r["loss"],7),r["donor_margin"],
                  round(r["stationarity"],5),flush=True)
    for mode in modes:
        rows = sorted(results[mode],key=lambda r:(not stable(r),r["loss"]))
        best=rows[0]
        assert stable(best), (mode,best)
        assert all(r["donor_margin"]<0 or
                   likelihood_ratio(best["loss"],r["loss"],221)<.02 for r in rows)
        write(HERE/f"refits_{mode}.tsv",rows)
        write(HERE/f"best_{mode}.tsv",[best])


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=[*MODES,"initial"])
    p.add_argument("--workers",type=int,default=7)
    a=p.parse_args()
    refit(list(MODES)[:3] if a.mode=="initial" else [a.mode],a.workers)
