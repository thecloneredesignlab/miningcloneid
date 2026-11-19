"""Check likelihood calibration and a shared count-level correlation sensitivity."""
import argparse
from collections import defaultdict
import concurrent.futures as cf

import numpy as np
from scipy.optimize import minimize_scalar

from count_model import HERE, OLS, model, optimize
from inference import THRESHOLD, SLACK, likelihood_ratio, pool
from paired_counts import loads
from practical import Configuration, Kernel, MECHANISMS, read, row_x, solve_constrained, write


def current_priority():
    p=read(HERE/"profile_intervals.tsv")
    selected={"lam_max"}
    for r in p:
        x=float(r["normalized_estimate"])
        gap=max(x-float(r["lowest_accepted_normalized"]),float(r["highest_accepted_normalized"])-x)
        if r["profile_precision"]!="outside_supported" or gap<=.13:
            selected.add(r["parameter"])
    return [p for p in MECHANISMS if p in selected]


def covariances(rho):
    k=model()
    groups=k.group[k.capture]
    day=np.array([r["calendar_day"] for r in k.observations])[k.capture]
    distance=abs(day[:,None]-day)
    correlation=np.where(groups[:,None]==groups[None,:],rho**distance,0.)
    w=np.zeros((k.ols_n,k.n))
    w[groups,np.arange(k.n)]=k.weights[k.capture]
    return correlation,w@correlation@w.T


def estimate_rho():
    k=model()
    c=Configuration(k)
    best=read(HERE/"real_best.tsv")[0]
    logs=k.log_predict(c.unmap(row_x(c,best))[0])[0][k.capture]
    residual=k.y-logs

    def objective(u):
        cov=covariances(u)[0]
        precision=np.linalg.inv(cov)
        sign, logdet=np.linalg.slogdet(cov)
        assert sign==1
        return k.n*np.log((residual@precision@residual)/k.n)+logdet
    result=minimize_scalar(objective,method="bounded",bounds=(0.,.85),options={"xatol":.001})
    rho=float(result.x) if result.fun<objective(0.)-1e-6 else 0.
    return rho


def fit_sensitivity(args):
    mode,rho,starts,parameters=args
    k=model() if mode=="counts" else Kernel()
    covariance=covariances(rho)[0 if mode=="counts" else 1]
    c=Configuration(k,precision=np.linalg.inv(covariance))
    opt=[]
    for x in starts:
        opt.append(optimize(k,c,k.y,x,95,divisor=10.) if mode=="counts"
                   else solve_constrained(k,c,k.y,x,limit=95,best_iterate_tolerance=1e-6))
    valid=[o for o in opt if o["supply"]>=-SLACK]
    pick=min(valid,key=lambda o:o["loss"])
    out=[]
    for p in parameters:
        j=MECHANISMS.index(p)
        x=float(pick["x"][j])
        for sign in (-1,1):
            point=x+sign*.10001
            if not 0.<=point<=1.:
                continue
            opts=[optimize(k,c,k.y,a,75,fixed=(j,point),divisor=s) if mode=="counts"
                  else solve_constrained(k,c,k.y,a,limit=75,fixed=(j,point),best_iterate_tolerance=1e-6)
                  for a,s in ((pick["x"],10.),(starts[0],2.))]
            candidate=[q for q in opts if q["supply"]>=-SLACK]
            best=min(candidate or opts,key=lambda q:q["loss"])
            r=likelihood_ratio(best["loss"],pick["loss"],k.n) if best["supply"]>=-SLACK else float("inf")
            out.append(dict(mode=mode,rho=rho,parameter=p,center_normalized=x,
                            baseline_loss=pick["loss"],baseline_stationarity=pick["gradient"],
                            fixed_normalized=point,likelihood_ratio=r,stationarity=best["gradient"],
                            feasible=best["supply"]>=-SLACK))
    return out


def sensitivity():
    k=model()
    c=Configuration(k)
    inferred=estimate_rho()
    rho=inferred if inferred>=.05 else .35
    files=(HERE/"real_best.tsv",OLS/"main_constrained_best.tsv")
    best=[row_x(c,read(p)[0]) for p in files]
    priority=current_priority()
    write(HERE/"correlation_design.tsv", [dict(fitted_rho=inferred,used_rho=rho,
           choice="residual fit" if inferred>=.05 else "prespecified perturbation",
           shared_covariance="within-passage calendar-day AR1, projected through exact W for OLS",
           parameters=",".join(priority))])
    out=[]
    with pool(2) as job:
        tasks=[job.submit(fit_sensitivity,(mode,rho,best if mode=="counts" else best[::-1],priority))
               for mode in ("counts","ols")]
        for f in cf.as_completed(tasks):
            results=f.result()
            out.extend(results)
            write(HERE/"correlation_sensitivity.tsv",out)
            print("correlation sensitivity",results[0]["mode"],len(results),"rho",rho,flush=True)


def calibrate_one(task):
    a,i,x,loss,truth,y,parameters = task
    k=model()
    c=Configuration(k)
    out=[]
    for name in parameters:
        j=MECHANISMS.index(name)
        fixed=x.copy()
        fixed[j]=truth[j]
        p=c.unmap(fixed)[0]
        initial,_,initial_margin=k.constrained(p,y,c,False)
        if abs(x[j]-truth[j])<1e-7 and initial_margin>=-SLACK:
            f,margin,stationary=initial,initial_margin,True
        else:
            result=optimize(k,c,y,x,limit=50,fixed=(j,truth[j]),divisor=10.)
            if result["supply"]<-SLACK or result["gradient"]>.025:
                other=optimize(k,c,y,truth,limit=53,fixed=(j,truth[j]),divisor=2.)
                if other["supply"]>=-SLACK and (result["supply"]<-SLACK or other["loss"]<result["loss"]):
                    result=other
            f,margin=result["loss"],result["supply"]
            stationary=result["gradient"]<=.025
            if initial_margin>=-SLACK and (margin<-SLACK or initial<f):
                f,margin,stationary=initial,initial_margin,False
        lr=likelihood_ratio(f,loss,k.n) if margin>=-SLACK else float("inf")
        out.append(dict(scenario=a,dataset=i,parameter=name,lr_at_truth=lr,
                        checked_fixed_fit=stationary and margin>=-SLACK,
                        covered_nominal=lr<=THRESHOLD,evaluable=lr<=THRESHOLD or (stationary and margin>=-SLACK)))
    return out


def calibrate(workers,maximum):
    k=model()
    c=Configuration(k)
    targets=loads()
    true={r["scenario"]:np.array([float(r[n]) for n in c.names]) for r in read(OLS/"synthetic_truth.tsv")}
    current=read(HERE/"qualified_count_recovery.tsv")
    path=HERE/"coverage_checks.tsv"
    output=read(path) if path.exists() else []
    done={(r["scenario"],int(r["dataset"]),r["parameter"]) for r in output}
    parameters=current_priority()
    args=[]
    for r in current:
        a,i=r["scenario"],int(r["dataset"])
        more=[name for name in parameters if (a,i,name) not in done]
        if i<=maximum and r["qualified"]=="True" and more:
            args.append((a,i,row_x(c,r),float(r["loss"]),true[a],targets[(a,i)],more))
    with pool(workers) as job:
        for f in cf.as_completed([job.submit(calibrate_one,a) for a in args]):
            result=f.result()
            output.extend(result)
            print("coverage",result[0]["scenario"],result[0]["dataset"],flush=True)
            write(path,output)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("correlation","coverage"))
    parser.add_argument("--workers",type=int,default=3)
    parser.add_argument("--max",type=int,default=100)
    args=parser.parse_args()
    if args.mode=="correlation":
        sensitivity()
    else:
        calibrate(args.workers,args.max)
