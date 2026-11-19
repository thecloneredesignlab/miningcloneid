"""Fit the observed trajectory counts with estimated within-passage correlation."""
import argparse
from concurrent.futures import as_completed

import numpy as np

from ar1_model import (COUNTS, HERE, RHO_MAX, STATIONARY, Config, choose, likelihood_ratio, model,
                       pool, read, receipt, row_x, solve, write)


def worker(task):
    label,start,divisor,limit=task
    k=model()
    c=Config(k)
    f=solve(k,c,k.y,start,limit=limit,divisor=divisor)
    return receipt(k,c,k.y,f,origin=label)


def fit_tasks(tasks,workers,existing=()):
    rows=list(existing)
    with pool(min(workers,len(tasks))) as executor:
        for f in as_completed([executor.submit(worker,t) for t in tasks]):
            r=f.result()
            rows.append(r)
            print("real",r["origin"],r["loss"],r["donor_margin"],r["stationarity"],r["rho"],flush=True)
            write(HERE/"observed_refits.tsv",rows)
    return rows


def save(rows):
    k=model()
    c=Config(k)
    best=choose(rows)
    assert float(best["donor_margin"])>=0
    write(HERE/"observed_best.tsv",[best])
    ref=float(best["loss"])
    stable=sorted((r for r in rows if float(r["stationarity"])<=STATIONARY and
              float(r["donor_margin"])>=0 and likelihood_ratio(float(r["loss"]),ref,k.n)<=.10),
              key=lambda r:float(r["loss"]))
    write(HERE/"reference_checks.tsv",[dict(qualified_origins=len(stable),reference_loss=ref,
        n=k.n,main_stationarity=best["stationarity"],main_strict_donor_margin=best["donor_margin"],
        lower_loss=ref,rho=best["rho"],sigma=best["sigma"])])
    p,_=c.unmap(row_x(c,best))
    logpred=k.predict_logs(p)
    residual=k.y-logpred
    innovations=residual.copy()
    t=k.transition
    a=p[-1]**k.gaps[t]
    innovations[t]=(residual[t]-a*residual[t-1])/np.sqrt(1-a*a)
    write(HERE/"observed_residuals.tsv",[dict(passage_id=s["passage_id"],cohort=s["cohort"],
        oxygen_pct=s["oxygen_pct"],record_id=s["record_id"],calendar_day=s["calendar_day"],
        residual=z,whitened_residual=u,observed_log_count=y,predicted_log_count=q)
        for s,z,u,y,q in zip((s for s in k.observations if s["fitted_count"]),residual,innovations,k.y,logpred)])
    print("Observed reference",best["origin"],best["loss"],best["rho"],best["sigma"],"stable",len(stable),flush=True)
    return best


def starts(c,k):
    old=read(COUNTS/"real_best.tsv")[0]
    anchor=c.load_old(old)
    all_points=[anchor,c.x_anchor]
    ols=COUNTS.with_name("20261118_experimental_OLS_practical_01")
    all_points.append(c.load_old(read(ols/"main_constrained_best.tsv")[0],.80))
    archive=read(COUNTS/"real_fits.tsv")+read(COUNTS/"profiles_direct.tsv")
    stable_neglected=sorted((r for r in archive if float(r["min_aggregate_donor_margin_fraction"])>=0
                              and float(r["max_projected_gradient"])<.03),key=lambda r:float(r["loss"]))
    for r in stable_neglected:
        z=c.load_old(r)
        if min(np.linalg.norm(z-q) for q in all_points)>.24:
            all_points.append(z)
        if len(all_points)>=9:
            break
    rng=np.random.default_rng(20839)
    all_points.extend(np.clip(anchor+rng.normal(0,sd,len(anchor)),0,1) for sd in (.007,.035))
    alternate=anchor.copy()
    alternate[-1]=.42/RHO_MAX
    all_points.append(alternate)
    return [(f"start_{i}",z,12. if i%3 else 4.,88) for i,z in enumerate(all_points,1)]


def refine(anchor,workers,label="refine",starts_extra=()):
    k=model()
    c=Config(k)
    x=row_x(c,anchor)
    rng=np.random.default_rng(234080+len(read(HERE/"observed_refits.tsv")))
    starts=[(x,25.)]
    starts += [(np.clip(x+rng.normal(0,sd,len(x)),0,1),dv)
                for sd,dv in ((.00008,4.),(.001,14.),(.004,28.),(.015,10.))]
    starts.extend((z,14.) for z in starts_extra)
    existing=read(HERE/"observed_refits.tsv")
    tasks=[(f"{label}_{i}",x,scale,112) for i,(x,scale) in enumerate(starts,1)]
    rows=fit_tasks(tasks,workers,existing)
    return save(rows)


def run(workers):
    k=model()
    c=Config(k)
    rows=fit_tasks(starts(c,k),workers)
    best=save(rows)
    ref=read(HERE/"reference_checks.tsv")[0]
    if int(ref["qualified_origins"])<2 or float(ref["main_stationarity"])>STATIONARY:
        best=refine(best,min(workers,5))
    return best


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=9)
    a=parser.parse_args()
    run(a.workers)
