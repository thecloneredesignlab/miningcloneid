"""Continue the frozen Monte Carlo design after a qualified 20-by-2 pilot."""
import argparse
import concurrent.futures as cf
from collections import defaultdict

import numpy as np

from fit_stage import model
from inference import pool
from practical import HERE,Configuration,read,row_x,solve_constrained,result_row,write
from recovery import characterize,scenarios,draw,pilot_gate


def worker(task):
    case, dataset, y, points = task
    k=model(); c=Configuration(k)
    rng=np.random.default_rng(6601+dataset+(1000 if case.startswith("B") else 0))
    results=[]
    for a in points:
        start=np.clip(a+rng.normal(0,.004,len(a)),0,1)
        results.append(solve_constrained(k,c,y,start,limit=45))
    best,ok,delta=characterize(k,results)
    if not ok:
        alternate=np.clip(best["x"]+rng.normal(0,.0004,len(best["x"])),0,1)
        results.append(solve_constrained(k,c,y,alternate,limit=42))
        best,ok,delta=characterize(k,results)
    return result_row(c,best,scenario=case,dataset=dataset,qualified=ok,
                      independent_LR_gap=delta,attempts=len(results))


def run(workers):
    k=model()
    c=Configuration(k)
    pilot=read(HERE/"qualified_recovery.tsv")
    assert pilot_gate(pilot), "pilot stability did not pass"
    truths,sigma=scenarios(k,c)
    complete_data=draw(k,c,truths,sigma,100)
    path=HERE/"expanded_refits.tsv"
    done=read(path) if path.exists() else []
    present={(r["scenario"],int(r["dataset"])) for r in done}
    tasks=[]
    for (a,i),y in complete_data.items():
        if i<=20 or (a,i) in present:
            continue
        # Blind starts drawn from preceding pilot fits. Select on distance
        # in observed slope space, never distance to true coordinates.
        bank=[r for r in pilot if int(r["dataset"])<=20 and r["qualified"]=="True"]
        potential=[(np.linalg.norm((y-complete_data[r["scenario"],int(r["dataset"])])/k.scale),row_x(c,r))
                   for r in bank]
        potential=sorted(potential,key=lambda a:a[0])
        tasks.append((a,i,y,[z for _,z in potential[:2]]))
    write(HERE/"noisy_slope_inputs.tsv", [
        dict(scenario=case,dataset=dataset,passage_id=p,ols_slope=x)
        for (case,dataset),ys in complete_data.items() for p,x in zip(k.ids,ys)])
    with pool(workers) as p:
        for f in cf.as_completed([p.submit(worker,x) for x in tasks]):
            r=f.result();done.append(r);write(path,done)
            print(r["scenario"],r["dataset"],r["qualified"],r["independent_LR_gap"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--workers",type=int,default=8)
    run(p.parse_args().workers)
