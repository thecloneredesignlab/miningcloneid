"""Cross-check precision-critical profiles from previously feasible fitted branches."""
import argparse
import concurrent.futures as cf

import numpy as np

from count_model import HERE, OLS, model, optimize
from inference import SLACK, likelihood_ratio, pool, THRESHOLD
from practical import Configuration, MECHANISMS, read, result_row, row_x, write
from profile_counts import load_direct


def fit(task):
    p,value,reference,candidates=task
    k=model()
    c=Configuration(k)
    j=MECHANISMS.index(p)
    near=sorted(candidates,key=lambda x:abs(x[j]-value))[:35]
    near += [reference]
    scored=[]
    for x in near:
        a=x.copy()
        a[j]=value
        f,_,g=k.constrained(c.unmap(a)[0],k.y,c,False)
        scored.append((a,f,g))
    scored.sort(key=lambda z:(z[2]<-SLACK,z[1] if z[2]>=-SLACK else -z[2]))
    starts=[]
    for x,_,_ in scored:
        if all(np.linalg.norm(x-other)>.017 for other in starts):
            starts.append(x)
        if len(starts)>=4:
            break
    runs=[]
    for x,divisor in zip(starts,(25.,5.,50.,14.)):
        out=optimize(k,c,k.y,x,75,(j,value),divisor)
        runs.append(out)
        good=sorted((r for r in runs if r["supply"]>=-SLACK and r["gradient"]<=.025),
                    key=lambda r:r["loss"])
        if len(good)>1 and likelihood_ratio(good[1]["loss"],good[0]["loss"],k.n)<.05:
            break
    valid=[r for r in runs if r["supply"]>=-SLACK]
    lowest=min(valid or runs,key=lambda r:r["loss"])
    verified=[r for r in valid if r["gradient"]<=.025 and
              likelihood_ratio(r["loss"],lowest["loss"],k.n)<=.001]
    best=min(verified,key=lambda r:r["loss"]) if verified else lowest
    return result_row(c,best,parameter=p,coordinate=value,feasible=best["supply"]>=-SLACK,
                      analytical_infeasible=False)


def run(workers):
    k=model()
    c=Configuration(k)
    reference=read(HERE/"real_best.tsv")[0]
    x=row_x(c,reference)
    profiles=read(HERE/"profile_intervals.tsv")
    sources=load_direct()
    feasible=[r for r in sources+read(HERE/"real_fits.tsv")
              if float(r["min_aggregate_donor_margin_fraction"])>=-SLACK]
    prior=[]
    for name in ("main_profiles_fast.tsv","main_profiles_edge_refit.tsv",
                 "main_profiles_extension.tsv","main_profiles_final_crossings.tsv"):
        f=OLS/name
        if f.exists():
            prior.extend(r for r in read(f) if r["feasible"]=="True")
    candidates=np.unique(np.array([row_x(c,r) for r in feasible+prior]),axis=0)
    jobs=[]
    for r in profiles:
        if r["profile_precision"]=="outside_supported":
            continue
        p=r["parameter"]
        i=MECHANISMS.index(p)
        grid=set()
        for sign,tag in ((-1,"lower"),(1,"upper")):
            probe=x[i]+sign*.10001
            if 0<probe<1:
                grid.add(probe)
            edge=r[f"{tag}_normalized"]
            if edge and 1e-5<float(edge)<1-1e-5 and r[f"{tag}_endpoint_status"]=="approximate_crossing":
                grid.add(float(edge))
        own=[z for z in sources if z["parameter"]==p and z["feasible"]=="True"]
        for v in sorted(grid):
            already_confirmed=any(abs(float(z["coordinate"])-v)<5e-7 and float(z["max_projected_gradient"])<=.025
                                  for z in own)
            if not already_confirmed:
                jobs.append((p,v,x,candidates))
    if not jobs:
        return
    path=HERE/"profiles_cross.tsv"
    rows=read(path) if path.exists() else []
    completed={(r["parameter"],float(r["coordinate"])) for r in rows}
    with pool(workers) as w:
        for future in cf.as_completed([w.submit(fit,t) for t in jobs if t[:2] not in completed]):
            r=future.result()
            rows.append(r)
            print("cross",r["parameter"],r["coordinate"],r["max_projected_gradient"],
                  likelihood_ratio(float(r["loss"]),float(reference["loss"]),k.n),flush=True)
            write(path,rows)


def polish_fit(job):
    name,x,val,divisor=job
    k=model()
    c=Configuration(k)
    f=optimize(k,c,k.y,x,105,(MECHANISMS.index(name),val),divisor)
    return result_row(c,f,parameter=name,coordinate=val,feasible=f["supply"]>=-SLACK,
                      analytical_infeasible=False)


def polish(workers):
    k=model()
    c=Configuration(k)
    ref=read(HERE/"real_best.tsv")[0]
    direct=load_direct()
    raw=read(HERE/"profile_intervals.tsv")
    tasks=[]
    rng=np.random.default_rng(520218)
    for r in raw:
        if r["profile_precision"]!="unresolved":
            continue
        name=r["parameter"]
        center=float(r["normalized_estimate"])
        for side in (-1,1):
            point=center+side*.10001
            if not 0.<point<1.:
                continue
            alternatives=[z for z in direct if z["parameter"]==name and z["feasible"]=="True"
                          and abs(float(z["coordinate"])-point)<.003]
            if any(float(z["max_projected_gradient"])<=.025 for z in alternatives):
                continue
            source=min(alternatives,key=lambda z:float(z["loss"])) if alternatives else ref
            point_x=row_x(c,source)
            for jitter,div in ((0.,3.),(0.,16.),(.00003,10.),(.00015,40.)):
                tasks.append((name,np.clip(point_x+rng.normal(0,jitter,len(point_x)),0,1),point,div))
    rows=[]
    if tasks:
        with pool(workers) as p:
            for f in cf.as_completed([p.submit(polish_fit,t) for t in tasks]):
                r=f.result()
                rows.append(r)
                print("polish",r["parameter"],r["coordinate"],r["max_projected_gradient"],r["loss"],flush=True)
                write(HERE/"profiles_polish.tsv",rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=3)
    parser.add_argument("--polish",action="store_true")
    args=parser.parse_args()
    if args.polish:
        polish(args.workers)
    else:
        run(args.workers)
