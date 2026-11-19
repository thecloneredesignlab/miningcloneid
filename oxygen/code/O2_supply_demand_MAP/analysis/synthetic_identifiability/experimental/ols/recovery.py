"""Blinded repeated recovery, conditional on the original experimental OLS design."""
import argparse
import concurrent.futures as cf
import zlib

import numpy as np

from inference import pool, THRESHOLD, likelihood_ratio
from fit_stage import model
from practical import (HERE, MECHANISMS, Configuration, read, write, row_x,
                       solve_constrained, result_row)

CONVERGE_LR = .10
PILOT_REQUIRED = 18
SUPPLY_TOL = 1e-8
KKT_TOL = .025


def scenarios(k, c):
    truth_file = HERE/"synthetic_truth.tsv"
    design_file = HERE/"truth_design.tsv"
    if truth_file.exists() and design_file.exists():
        truth = {r["scenario"]:np.array([float(r[n]) for n in c.names]) for r in read(truth_file)}
        assert set(truth)=={"A_real_fit","B_weak_direction"}
        sigma = float(read(design_file)[0]["sigma_log_count"])
        return truth, sigma
    ref = row_x(c, read(HERE/"main_constrained_best.tsv")[0])
    center = np.clip(ref,1e-5,1.-1e-5)
    grad = np.empty((k.n,len(ref)))
    for j in range(len(ref)):
        aa,bb=center.copy(),center.copy()
        aa[j]+=1e-6
        bb[j]-=1e-6
        fa=k.predict(c.unmap(aa)[0])[0]
        fb=k.predict(c.unmap(bb)[0])[0]
        grad[:,j]=(fa-fb)/(2e-6*k.scale)
    _,s,vt=np.linalg.svd(grad,full_matrices=False)
    denom=2*float(read(HERE/"main_constrained_best.tsv")[0]["loss"])/k.n
    found=[]
    f0=k.predict(c.unmap(ref)[0])[0]
    directions=[(f"svd_{j}", vt[-j]) for j in range(1,7)]
    # At an active constraint, clipping every coordinate inward may break
    # feasibility. Include unit directions and retain the other coordinates
    # exactly at the fitted reference while nominating a feasible weak case.
    directions.extend((f"coordinate_{c.names[j]}", np.eye(len(ref))[j]) for j in range(14))
    for direction,weak in directions:
        weak=weak/np.max(np.abs(weak))
        for sign in (-1,1):
            for step in (.32,.24,.16):
                candidate=np.clip(ref+sign*step*weak,0.,1.)
                if np.max(np.abs(candidate[:14]-ref[:14]))<.105:
                    continue
                f,_,diag,_,_=k.predict(c.unmap(candidate)[0])
                if diag[0]==0:
                    dist=float(np.linalg.norm((f-f0)/k.scale))
                    found.append((dist, direction, float(step),candidate))
    if not found:
        raise RuntimeError("No distinct feasible alternate truth on the weak singular directions")
    distance,direction,step,other=min(found,key=lambda item:item[0])
    write(HERE/"truth_design.tsv", [dict(scenario="A_real_fit",weighted_slope_distance=0.,singular_direction=0,
                                         step=0.,sigma_log_count=np.sqrt(denom), source="main_constrained_best.tsv"),
                                   dict(scenario="B_weak_direction",weighted_slope_distance=distance,
                                        singular_direction=direction,step=step,sigma_log_count=np.sqrt(denom),
                                        source="weakest tested feasible singular or coordinate direction")])
    write(HERE/"synthetic_truth.tsv", [dict(scenario=r, **{n:float(z) for n,z in zip(c.names,x)})
                                       for r,x in (("A_real_fit",ref),("B_weak_direction",other))])
    return {"A_real_fit":ref,"B_weak_direction":other},np.sqrt(denom)


def characterize(k, fits):
    feasible=[f for f in fits if f["supply"]>=-SUPPLY_TOL and f["loss"]<100]
    if not feasible:
        return min(fits,key=lambda f:f["loss"]),False,np.inf
    best=min(feasible,key=lambda f:f["loss"])
    other=sorted(likelihood_ratio(f["loss"],best["loss"],k.n) for f in feasible
                 if f is not best and f["gradient"]<=KKT_TOL)
    difference=other[0] if other else np.inf
    okay=best["gradient"]<=KKT_TOL and difference<=CONVERGE_LR
    return best,okay,difference


def fit_dataset(args):
    case, dataset, target = args
    k=model()
    conf=Configuration(k)
    rng=np.random.default_rng(zlib.crc32(f"{case}:{dataset}".encode())+471706)
    # Register the production structural witness as center independently of
    # either generated truth. Use distinct blind starts. Neither exact truth
    # nor distance-to-truth is an input to this function.
    preliminary=[r for r in read(HERE/"real_fits.tsv")
                 if float(r["min_aggregate_donor_margin_fraction"])>=0]
    # Both generator truths are absent; this pre-existing real-data search
    # bank contains earlier candidates from a different estimation run.
    starts=[np.clip(conf.x_anchor+rng.normal(0,.035,len(conf.names)),0,1),
            np.clip(row_x(conf,preliminary[min(4,len(preliminary)-1)])
                    +rng.normal(0,.035,len(conf.names)),0,1)]
    fits=[solve_constrained(k,conf,target,x,limit=54) for x in starts]
    winner,okay,agreement=characterize(k,fits)
    if not okay:
        # One continuation attempt and an extra independent start.
        repaired=solve_constrained(k,conf,target,winner["x"],limit=68)
        retry=solve_constrained(k,conf,target,
                                np.clip(conf.x_anchor+rng.normal(0,.13,len(conf.names)),0,1),limit=64)
        fits.extend([repaired,retry])
        # A restart from the same incumbent does not count as another origin.
        separate=[fits[0],fits[1],fits[3]]
        if repaired["loss"]<min(separate,key=lambda f:f["loss"])["loss"]:
            j=min(range(3),key=lambda i:separate[i]["loss"])
            separate[j]=repaired
        winner,okay,agreement=characterize(k,separate)
    result=result_row(conf,winner,scenario=case,dataset=dataset,qualified=okay,
                      independent_LR_gap=agreement,attempts=len(fits))
    return result


def draw(k, c, known, sigma, count):
    samples={}
    for case, x in known.items():
        raw=k.predict(c.unmap(x)[0])[1]
        base=np.log(raw)
        for dataset in range(1,count+1):
            rng=np.random.default_rng(np.random.SeedSequence(
                [20261118,1 if case.startswith("A") else 2,dataset]))
            e=np.zeros(len(base))
            e[k.capture]=rng.normal(0,sigma,len(base[k.capture]))
            log_counts=np.maximum(0,base+e)
            observed=np.bincount(k.group,k.weights*log_counts,minlength=k.n)
            samples[case,dataset]=observed
    return samples


def pilot_gate(rows):
    report=[]
    for case in ("A_real_fit","B_weak_direction"):
        sample=[r for r in rows if r["scenario"]==case and int(r["dataset"])<=20]
        assert len(sample)==20
        passes=sum(str(r["qualified"])=="True" for r in sample)
        report.append(dict(scenario=case,assigned=20,qualified=passes,
                           released_100=passes>=PILOT_REQUIRED,required=PILOT_REQUIRED,
                           max_independent_lr_gap_allowed=CONVERGE_LR,kkt_norm_allowed=KKT_TOL))
    write(HERE/"pilot_gate.tsv",report)
    return all(r["released_100"] for r in report)


def simulate(count,workers):
    k=model()
    config=Configuration(k)
    known,sigma=scenarios(k,config)
    observed=draw(k,config,known,sigma,count)
    write(HERE/"noisy_slope_inputs.tsv", [
        dict(scenario=case,dataset=dataset,passage_id=p,ols_slope=x)
        for (case,dataset),ys in observed.items() for p,x in zip(k.ids,ys)])
    path=HERE/"recovery_fits.tsv"
    previous=read(path) if path.exists() else []
    done={(r["scenario"],int(r["dataset"])) for r in previous}
    results=[r for r in previous if int(r["dataset"])<=count]
    if count==100 and not pilot_gate(results):
        raise RuntimeError("20-replicate pilot did not pass the optimizer stability gate")
    tasks=[(a,i,y) for (a,i),y in observed.items() if (a,i) not in done]
    with pool(workers) as run:
        futures=[run.submit(fit_dataset,t) for t in tasks]
        for future in cf.as_completed(futures):
            r=future.result()
            results.append(r)
            print(r["scenario"],r["dataset"],r["qualified"],r["loss"],
                  r["independent_LR_gap"],flush=True)
            write(path,sorted(results,key=lambda z:(z["scenario"],int(z["dataset"]))))
    if count>=20:
        print("Pilot expanded release:", pilot_gate(results),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("n",type=int,choices=(20,100))
    parser.add_argument("--workers",type=int,default=7)
    args=parser.parse_args()
    simulate(args.n,args.workers)
