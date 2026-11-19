"""Resolve competing optimizer basins on a blinded pilot. Does not consume truth."""
import argparse
import concurrent.futures as cf
from collections import defaultdict

import numpy as np

from fit_stage import model
from inference import pool
from practical import HERE,Configuration,read,row_x,solve_constrained,result_row,write
from recovery import characterize


def repair(task):
    sample, target, bank = task
    k=model()
    c=Configuration(k)
    i=int(sample["dataset"])
    rng=np.random.default_rng(i+int(sample["scenario"].startswith("B"))*1000+81247)
    previous=dict(x=row_x(c,sample), loss=float(sample["loss"]),
                  supply=float(sample["min_aggregate_donor_margin_fraction"]),
                  gradient=float(sample["max_projected_gradient"]),
                  scipy_success=sample["scipy_success"]=="True",calls=int(sample["calls"]),
                  iteration=int(sample["iterations"]),termination=sample["termination"],
                  seconds=float(sample["seconds"]))
    # Candidate basins derive from fits to the actual experiment or prior
    # blinded optimizer output, never from either synthetic truth.
    fresh=[np.clip(bank+rng.normal(0,jitter,len(bank)),0,1) for jitter in (.0035,.009)]
    fresh.append(np.clip(previous["x"]+rng.normal(0,.006,len(bank)),0,1))
    solutions=[previous]+[solve_constrained(k,c,target,x,limit=45) for x in fresh]
    winner,passed,gap=characterize(k,solutions)
    return result_row(c,winner,scenario=sample["scenario"],dataset=i,qualified=passed,
                      independent_LR_gap=gap,attempts=int(sample["attempts"])+len(fresh))


def run(workers):
    k=model()
    c=Configuration(k)
    source=read(HERE/"main_constrained_refinement.tsv")
    bank=row_x(c,next(r for r in source if r["origin"]=="refine_5"))
    a=read(HERE/"recovery_fits.tsv")
    path=HERE/"recovery_repairs.tsv"
    rows=read(path) if path.exists() else []
    stored={(q["scenario"],int(q["dataset"])) for q in rows}
    y=defaultdict(list)
    for q in read(HERE/"noisy_slope_inputs.tsv"):
        y[q["scenario"],int(q["dataset"])].append(float(q["ols_slope"]))
    work=[(z,np.asarray(y[z["scenario"],int(z["dataset"])]),bank) for z in a
          if z["qualified"]=="False" and (z["scenario"],int(z["dataset"])) not in stored]
    with pool(workers) as p:
        for ret in cf.as_completed([p.submit(repair,v) for v in work]):
            row=ret.result()
            rows.append(row)
            write(path,rows)
            print('repair',row['scenario'],row['dataset'],row['qualified'],row['independent_LR_gap'],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--workers",type=int,default=2)
    run(p.parse_args().workers)
