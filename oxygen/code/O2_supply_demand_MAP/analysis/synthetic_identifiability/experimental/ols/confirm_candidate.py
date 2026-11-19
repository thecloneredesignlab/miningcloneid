"""Independently restart the best blinded pilot candidate near its active branch."""
import concurrent.futures as cf
from collections import defaultdict
import argparse

import numpy as np

from fit_stage import model
from inference import pool
from practical import HERE,Configuration,read,row_x,solve_constrained,result_row,write
from recovery import characterize


def task(rec):
    row, y = rec
    k=model()
    c=Configuration(k)
    x=row_x(c,row)
    r=np.random.default_rng(int(row["dataset"]) + 701 + (10000 if row["scenario"].startswith("B") else 0))
    opt=dict(x=x,loss=float(row["loss"]),supply=float(row["min_aggregate_donor_margin_fraction"]),
             gradient=float(row["max_projected_gradient"]),calls=int(row["calls"]),
             iteration=int(row["iterations"]),scipy_success=row["scipy_success"]=="True",
             termination=row["termination"],seconds=float(row["seconds"]))
    check=[opt]
    if opt["gradient"]>=.025:
        check[0]=solve_constrained(k,c,y,x,limit=56)
        x=check[0]["x"]
    for sd in (.0002,.001):
        start=np.clip(x+r.normal(0,sd,len(x)),0.,1.)
        check.append(solve_constrained(k,c,y,start,limit=36))
        best,okay,gap=characterize(k,check)
        if okay:
            break
    return result_row(c,best,scenario=row["scenario"],dataset=int(row["dataset"]),qualified=okay,
                      independent_LR_gap=gap,attempts=int(row["attempts"])+len(check)-1)


def run(workers,suffix="",scenario=None):
    path=HERE/f"pilot_near_branch_confirmation{suffix}.tsv"
    rows=read(path) if path.exists() else []
    done={(r["scenario"],int(r["dataset"])) for r in rows}
    targets=defaultdict(list)
    for r in read(HERE/"noisy_slope_inputs.tsv"):
        targets[r["scenario"],int(r["dataset"])].append(float(r["ols_slope"]))
    jobs=[]
    for r in read(HERE/"qualified_recovery.tsv"):
        key=(r["scenario"],int(r["dataset"]))
        if r["qualified"]=="False" and key not in done and (
            scenario is None or r["scenario"].startswith(scenario.upper())
        ) and (not suffix.startswith("_expansion") or key[1]>20):
            jobs.append((r,np.array(targets[key])))
    with pool(workers) as w:
        for f in cf.as_completed([w.submit(task,t) for t in jobs]):
            r=f.result()
            rows.append(r);write(path,rows)
            print(r["scenario"],r["dataset"],r["qualified"],r["independent_LR_gap"],flush=True)


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--workers",type=int,default=3)
    a.add_argument("--suffix",choices=("","_late","_expansion","_expansion_b"),default="")
    a.add_argument("--scenario",choices=("a","b"))
    args=a.parse_args()
    run(args.workers,args.suffix,args.scenario)
