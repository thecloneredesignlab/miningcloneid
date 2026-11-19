"""Complete the same blinded pilot from the resolved experimental-data basin."""
import argparse
import concurrent.futures as cf
from collections import defaultdict

import numpy as np

from fit_stage import model
from inference import pool
from practical import HERE, Configuration, read, row_x, solve_constrained, result_row, write
from recovery import characterize


def one(item):
    scenario,dataset,y,bank=item
    k=model()
    cfg=Configuration(k)
    rng=np.random.default_rng(dataset + (1000 if scenario.startswith("B") else 0)+930561)
    fits=[]
    for sd in (.004,.010):
        x=np.clip(bank+rng.normal(0,sd,16),0,1)
        fits.append(solve_constrained(k,cfg,y,x,limit=43))
    best,okay,delta=characterize(k,fits)
    if not okay:
        x=np.clip(best["x"]+rng.normal(0,.003,16),0,1)
        fits.append(solve_constrained(k,cfg,y,x,limit=41))
        best,okay,delta=characterize(k,fits)
    return result_row(cfg,best,scenario=scenario,dataset=dataset,qualified=okay,
                      independent_LR_gap=delta,attempts=len(fits))


def run(pred,workers):
    k=model()
    c=Configuration(k)
    bank=row_x(c,next(r for r in read(HERE/"main_constrained_refinement.tsv") if r["origin"]=="refine_5"))
    data=defaultdict(list)
    for r in read(HERE/"noisy_slope_inputs.tsv"):
        data[r["scenario"],int(r["dataset"])].append(float(r["ols_slope"]))
    path=HERE/f"pilot_bank_{pred}.tsv"
    results=read(path) if path.exists() else []
    complete={(r["scenario"],int(r["dataset"])) for r in results}
    tasks=[(a,i,np.array(y),bank) for (a,i),y in data.items() if a.startswith(pred.upper()) and (a,i) not in complete]
    with pool(workers) as p:
        for job in cf.as_completed([p.submit(one,job) for job in tasks]):
            row=job.result()
            results.append(row)
            write(path,results)
            print("new blinded bank",row["scenario"],row["dataset"],row["qualified"],
                  row["loss"],row["independent_LR_gap"],flush=True)


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("case",choices=("a","b"))
    a.add_argument("--workers",type=int,default=3)
    c=a.parse_args()
    run(c.case,c.workers)
