"""Check confidence limits close to the precision decision."""
import argparse
import concurrent.futures as cf

import numpy as np

from count_model import HERE, model
from inference import THRESHOLD, likelihood_ratio, pool
from practical import Configuration, MECHANISMS, read, result_row, row_x, write
from profile_counts import envelope, fixed_solve, load_direct, refit_reference


def check(args):
    name,value,starts,limit = args
    k=model()
    c=Configuration(k)
    fit=fixed_solve(k,c,k.y,MECHANISMS.index(name),value,starts,limit=limit)
    return result_row(c,fit,parameter=name,coordinate=value,feasible=fit["supply"]>=-1e-8,
                      analytical_infeasible=False)


def run(workers,defer_reference):
    k=model()
    c=Configuration(k)
    base=read(HERE/"real_best.tsv")[0]
    best=row_x(c,base)
    intervals=read(HERE/"profile_intervals.tsv")
    existing=load_direct()
    jobs=[]
    for r in intervals:
        p=r["parameter"]
        x=float(r["normalized_estimate"])
        grade=r["profile_precision"]
        if grade=="outside_supported":
            continue
        if grade=="narrow_checked" and "approximate_crossing" not in (
                r["lower_endpoint_status"],r["upper_endpoint_status"]):
            continue
        candidates=set()
        for sign,end in ((-1,"lower"),(1,"upper")):
            probe=x+sign*.10001
            if 0<probe<1:
                candidates.add(probe)
            v=r[f"{end}_normalized"]
            if v and 1e-5<float(v)<1-1e-5:
                candidates.add(float(v))
        own=[z for z in existing if z["parameter"]==p and z["feasible"]=="True"
             and float(z["max_projected_gradient"])<=.025]
        stationary=[z for z in existing if z["feasible"]=="True" and float(z["max_projected_gradient"])<=.025]
        for y in sorted(candidates):
            ordered=sorted(own,key=lambda q:abs(float(q["coordinate"])-y))
            cross=sorted(stationary,key=lambda q:float(q["loss"]))
            starts=[best]
            for q in ordered[:3]+cross[:3]:
                a=row_x(c,q)
                if min(np.linalg.norm(a-z) for z in starts)>.0002:
                    starts.append(a)
                if len(starts)==3:
                    break
            jobs.append((p,y,starts,89))
    path=HERE/"profiles_edges.tsv"
    rows=[]
    if jobs:
        with pool(workers) as work:
            for f in cf.as_completed([work.submit(check,t) for t in jobs]):
                r=f.result()
                rows.append(r)
                print("edge",r["parameter"],r["coordinate"],r["max_projected_gradient"],
                      likelihood_ratio(float(r["loss"]),float(base["loss"]),k.n),flush=True)
                write(path,rows)
        if not defer_reference:
            refit_reference(workers)
            envelope()


def parallel_scan(name,workers):
    k=model()
    c=Configuration(k)
    base=read(HERE/"real_best.tsv")[0]
    x=float(base[f"x_{name}"])
    best=row_x(c,base)
    values=sorted({0.,.25,.5,.75,1.,max(0.,x-.10001),min(1.,x+.10001),
                   max(0.,x-.05),min(1.,x+.05)})
    rows=[]
    with pool(workers) as work:
        futures=[work.submit(check,(name,q,[best],48)) for q in values]
        for f in cf.as_completed(futures):
            r=f.result()
            rows.append(r)
            print("parallel scan",name,r["coordinate"],r["loss"],r["max_projected_gradient"],flush=True)
            write(HERE/"profiles_supplement.tsv",rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=5)
    parser.add_argument("--scan",choices=MECHANISMS)
    parser.add_argument("--defer-reference",action="store_true")
    args=parser.parse_args()
    if args.scan:
        parallel_scan(args.scan,args.workers)
    else:
        run(args.workers,args.defer_reference)
