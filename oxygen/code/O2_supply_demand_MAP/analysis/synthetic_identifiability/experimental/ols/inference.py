"""Constrained likelihood, continuation profiles and registered sensitivity models."""
import argparse
import concurrent.futures as cf
import multiprocessing
import os
from collections import defaultdict

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from practical import (HERE, MECHANISMS, Configuration, read, write,
                       row_x, solve_constrained, result_row)
from fit_stage import model

THRESHOLD = chi2.ppf(.95,1)
SLACK = 1e-8


def pool(workers):
    return cf.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn"))


def conf_for(k,name):
    if name=="main":
        return Configuration(k)
    if name=="passage":
        return Configuration(k,gaussian=True,split=True,tilt=True)
    if name=="noise":
        rho=float(read(HERE/"noise_assumptions.tsv")[0]["rho_fit"])
        return Configuration(k,precision=precision(k,rho)[0])
    raise ValueError(name)


def worker(args):
    case, origin, start, limit = args
    k=model()
    cfg=conf_for(k,case)
    fit=solve_constrained(k,cfg,k.y,start,limit)
    return result_row(cfg,fit,case=case,origin=origin)


def refit(case, workers, refine=False):
    k=model()
    c=conf_for(k,case)
    if case=="main":
        old=read(HERE/("main_constrained_fits.tsv" if refine else "real_fits.tsv"))
        x=[row_x(c,r) for r in old[:6]]
        if refine:
            rng=np.random.default_rng(2201)
            x=[x[0],*[
                np.clip(x[0]+rng.normal(0,sd,len(x[0])),0,1) for sd in (.004,.01,.025,.05)
            ]]
    else:
        c0=Configuration(k)
        best=read(HERE/"main_constrained_best.tsv")[0]
        base=c0.unmap(row_x(c0,best))[0]
        xbase=c.normal(base[c.active])
        rng=np.random.default_rng(1518)
        if refine:
            old=read(HERE/f"{case}_constrained_fits.tsv")
            x=[row_x(c,r) for r in old[:2]]
        else:
            x=[xbase]+[np.clip(xbase+rng.normal(0,.03,len(xbase)),0,1) for _ in range(2)]
    tasks=[(case,("refine_" if refine else "")+str(i+1),v,
            170 if refine else (100 if case!='passage' else 70)) for i,v in enumerate(x)]
    rows=[]
    with pool(workers) as work:
        for future in cf.as_completed([work.submit(worker,a) for a in tasks]):
            r=future.result()
            rows.append(r)
            print(case,r["origin"],r["loss"],r["min_aggregate_donor_margin_fraction"],r["max_projected_gradient"],
                  r["termination"],flush=True)
            write(HERE/f"{case}_constrained_{'refinement' if refine else 'fits'}.tsv",
                  sorted(rows,key=lambda q:q["loss"]))
    foot=[r for r in rows if r["min_aggregate_donor_margin_fraction"]>=-SLACK]
    if not foot:
        raise RuntimeError(f"all {case} starts infeasible")
    best=min(foot,key=lambda r:r["loss"])
    strict=[r for r in foot if r["min_aggregate_donor_margin_fraction"]>=0
            and r["loss"]<=best["loss"]+1e-5 and r["max_projected_gradient"]<.01]
    if strict:
        best=min(strict,key=lambda r:r["loss"])
    if refine and case=="main":
        old_best=read(HERE/f"{case}_constrained_best.tsv")[0]
        if float(old_best["loss"])<best["loss"]:
            raise RuntimeError("New attempts did not improve the incumbent; keep both result files")
    elif refine:
        old_best=read(HERE/f"{case}_constrained_best.tsv")[0]
        if float(old_best["loss"])<best["loss"]:
            best=old_best
    write(HERE/f"{case}_constrained_best.tsv",[best])
    if case=="main":
        record_noise(k,c,best)


def record_noise(k,c,best):
    raw,_=c.unmap(row_x(c,best))
    f,_,_,_,_=k.predict(raw)
    residual=k.y-f
    result=minimize_scalar(lambda rho:
        k.n*np.log((residual@precision(k,rho)[0]@residual)/k.n)+precision(k,rho)[1],
        bounds=(0,.85),method="bounded",options={"xatol":1e-3})
    # AR(1) is a sensitivity to correlation along the reconstructed tree.
    # Estimate rho at the main residuals once and condition sensitivity on it.
    rho=0. if result.fun>=k.n*np.log(np.sum(residual**2/k.var)/k.n) else float(result.x)
    write(HERE/"noise_assumptions.tsv",[dict(rho_fit=rho, method="tree_AR1_at_main_residuals",
                                             conditional_sensitivity=True,
                                             sigma_main=np.sqrt(2*float(best["loss"])/k.n))])
    write(HERE/"real_residuals_constrained.tsv", [dict(passage_id=i,observed=y,predicted=b,
                                                     residual=y-b,variance_multiplier=v)
                                                 for i,y,b,v in zip(k.ids,k.y,f,k.var)])


DISTANCE = None


def precision(k,rho):
    global DISTANCE
    if DISTANCE is None:
        parent={p["passage_id"]:p["parent_id"] for p in k.graph}
        lineage={}
        for v in k.ids:
            a={}
            u,d=v,0
            while u:
                a[u]=d
                u=parent[u]
                d+=1
            lineage[v]=a
        D=np.full((k.n,k.n),-1)
        for i,v in enumerate(k.ids):
            a=lineage[v]
            for j,w in enumerate(k.ids[:i+1]):
                b=lineage[w]
                shared=a.keys() & b.keys()
                if shared:
                    D[i,j]=D[j,i]=min(a[u]+b[u] for u in shared)
        DISTANCE=D
    correlation=np.zeros(DISTANCE.shape)
    valid=DISTANCE>=0
    correlation[valid]=rho**DISTANCE[valid]
    covariance=correlation*np.outer(k.scale,k.scale)
    sign,logdet=np.linalg.slogdet(covariance)
    assert sign==1.0
    return np.linalg.inv(covariance),logdet


def likelihood_ratio(loss,reference,n):
    return n*np.log(loss/reference) if loss>0 and reference>0 else np.inf


def lambda_supply_floor(k):
    # Even with no deaths, no loss and no crowding, x(t+dt) <= (1+dt*lambda_max)x(t).
    seed={r["passage_id"]:int(r["fixed_live_cells"]) for r in k.graph}
    demand=defaultdict(int)
    for r in k.graph:
        if r["parent_id"]:
            demand[r["parent_id"],int(r["donor_step"])]+=int(r["fixed_live_cells"])
    return max(np.expm1(np.log(d/seed[p])/ticks)/.05 for (p,ticks),d in demand.items() if ticks>0)


def profile_worker(args):
    case, param, coarse = args
    k=model()
    c=conf_for(k,case)
    base=read(HERE/f"{case}_constrained_best.tsv")[0]
    x0=row_x(c,base)
    primary=float(base["loss"])
    i=MECHANISMS.index(param)
    center=float(x0[i])
    grid={max(0,center-.10001),min(1,center+.10001)}
    if not coarse:
        grid.update([0.,.25,.50,.75,1.])
    rows={}
    prev={0:x0,1:x0}
    basefit=dict(x=x0,loss=primary,gradient=float(base["max_projected_gradient"]),
                 supply=float(base["min_aggregate_donor_margin_fraction"]),calls=0,
                 iteration=0,scipy_success=base["scipy_success"],termination="unconstrained reference",
                 seconds=0.)
    rows[center]=basefit

    def add(value, restart):
        if param=="lam_max":
            test=restart.copy()
            test[i]=value
            if c.unmap(test)[0][0] < lambda_supply_floor(k)-1e-8:
                fit=dict(x=test,loss=np.inf,supply=-np.inf,gradient=np.nan,calls=0,iteration=0,
                         scipy_success=False,termination="analytically infeasible by maximum supply",seconds=0.)
                rows[value]=fit
                return fit
        fit=solve_constrained(k,c,k.y,restart,limit=38 if not coarse else 32,fixed=(i,value))
        if fit["supply"] < -SLACK or likelihood_ratio(fit["loss"],primary,k.n)>THRESHOLD:
            other=solve_constrained(k,c,k.y,x0,limit=38 if not coarse else 32,fixed=(i,value))
            if other["supply"]>=-SLACK and (fit["supply"]<-SLACK or other["loss"]<fit["loss"]):
                fit=other
        rows[value]=fit
        return fit

    for direction in (0,1):
        targets=sorted([v for v in grid if (v<center if direction==0 else v>center)],
                       reverse=direction==0)
        for value in targets:
            fit=add(value,prev[direction])
            if fit["supply"]>=-SLACK:
                prev[direction]=fit["x"]
    if not coarse:
        for _ in range(5):
            targets=[]
            ordered=sorted(rows)
            for lo,hi in zip(ordered,ordered[1:]):
                ra,rb=rows[lo],rows[hi]
                if min(ra["supply"],rb["supply"])>=-SLACK:
                    la,lb=(likelihood_ratio(z["loss"],primary,k.n) for z in (ra,rb))
                    if (la>THRESHOLD)!=(lb>THRESHOLD) and hi-lo>.002:
                        targets.append(((lo+hi)/2,ra["x"] if la<lb else rb["x"]))
            if not targets:
                break
            for at,warm in targets:
                add(at,warm)
    result=[]
    for x,f in sorted(rows.items()):
        lr=likelihood_ratio(f["loss"],primary,k.n) if f["supply"]>=-SLACK else np.inf
        result.append(result_row(c,f,case=case,profile_parameter=param,fixed_normalized=x,
                                 lr=lr,level=THRESHOLD,feasible=f["supply"]>=-SLACK))
    return result


def profiles(case, workers, coarse=False):
    rows=[]
    with pool(workers) as w:
        jobs=[w.submit(profile_worker,(case,p,coarse)) for p in MECHANISMS]
        for result in cf.as_completed(jobs):
            r=result.result()
            rows.extend(r)
            print(case,r[0]["profile_parameter"],len(r),
                  "qualified",sum(x["feasible"] for x in r),
                  "minimum LR",min(x["lr"] for x in r),flush=True)
            suffix=os.environ.get("CULTURE_PROFILE_SUFFIX","")
            write(HERE/f"{case}_profiles{suffix}.tsv",
                  sorted(rows,key=lambda x:(x["profile_parameter"],x["fixed_normalized"])))
    return rows


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("command",choices=("fit","profile"))
    parser.add_argument("--case",choices=("main","passage","noise"),default="main")
    parser.add_argument("--workers",type=int,default=6)
    parser.add_argument("--coarse",action="store_true")
    parser.add_argument("--refine",action="store_true")
    args=parser.parse_args()
    if args.command=="fit":
        refit(args.case,args.workers,args.refine)
    else:
        profiles(args.case,args.workers,args.coarse)
