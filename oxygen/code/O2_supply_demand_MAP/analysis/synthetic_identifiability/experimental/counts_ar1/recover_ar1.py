"""Freeze two feasible generating cases and run blinded repeated AR(1) fits."""
import argparse
from collections import defaultdict
from concurrent.futures import as_completed
import zlib

import numpy as np

from ar1_model import (COUNTS, HERE, RHO_MAX, STATIONARY, TARGETS, THRESHOLD, Config, choose,
                       likelihood_ratio, model, pool, read, receipt, row_x, solve, write)
from profile_ar1 import load_profiles, old_bank, reference

PILOT_REQUIRED=18
MAX_AGREEMENT_LR=.10


def truth_worker(task):
    limits,start,stamp=task
    k=model()
    c=Config(k)
    f=solve(k,c,k.y,start,limit=92,fixed=limits,divisor=10.)
    return receipt(k,c,k.y,f,origin=stamp)


def freeze_truths(workers):
    k=model()
    c=Config(k)
    a=reference()
    x=row_x(c,a)
    f=float(a["loss"])
    rho=float(a["rho"])
    sigma=float(a["sigma"])
    wanted=[r for r in read(HERE/"profile_intervals.tsv") if r["profile_status"]=="profile_pass"]
    if not wanted:
        return
    shift={}
    for r in wanted:
        name=r["parameter"]
        if name=="lam_max" or name=="buffer_smax":
            center=float(r["normalized_estimate"])
            below=float(r["min_supported"])
            if center>.95:
                shift[c.names.index(name)]=max(center-below,.006)
    q={len(x)-1:x[-1]}
    jobs=[]
    for alpha in (.92,.80,.68,.55,.38):
        if shift:
            fixed=q | {j:x[j]-alpha*amount for j,amount in shift.items()}
            jobs.append((fixed,x,f"joint_inward_{alpha:.2f}"))
    profile=[r for r in load_profiles() if float(r["donor_margin"])>=0 and
             likelihood_ratio(float(r["loss"]),f,k.n)<=THRESHOLD]
    bank=[]
    for r in profile:
        z=row_x(c,r)
        z[-1]=x[-1]
        p,_=c.unmap(z)
        out,_=k.metrics(p,k.y)
        lr=likelihood_ratio(out[0],f,k.n)
        if out[1]>=0 and lr<=THRESHOLD and max(abs(z[:14]-x[:14]))>1e-5:
            shift_score=min((x[j]-z[j])/amount for j,amount in shift.items()) if shift else 0.
            bank.append((shift_score,max(abs(z[:14]-x[:14])),z,r["parameter"]))
    if bank:
        shapes=sorted(bank,key=lambda v:v[:2],reverse=True)
        for label,(_,_,z,name) in enumerate(shapes[:5]):
            idx=c.names.index(name)
            jobs.append((q | {idx:z[idx]},z,f"retained_profile_{name}_{label+1}"))
    assert jobs,"No simulation scenario alternatives were generated"
    candidates=[]
    with pool(min(workers,len(jobs))) as executor:
        for ftr in as_completed([executor.submit(truth_worker,t) for t in jobs]):
            r=ftr.result()
            candidates.append(r)
            print("truth candidate",r["origin"],r["donor_margin"],r["loss"],flush=True)
    # A directly verified point from a profile is an eligible fallback even if
    # its refinement ends at another stationary point.
    for _,_,z,name in bank:
        p,_=c.unmap(z)
        out,_=k.metrics(p,k.y)
        candidates.append(receipt(k,c,k.y,dict(x=z,gradient=np.nan,calls=0,iterations=0,
                          success=False,termination="retained feasible profile point",rounds=0,seconds=0),
                          origin=f"retained_{name}"))
    write(HERE/"truth_candidates.tsv",candidates)
    choices=[]
    for r in candidates:
        z=row_x(c,r)
        lr=likelihood_ratio(float(r["loss"]),f,k.n)
        if float(r["donor_margin"])<=0 or lr>THRESHOLD or max(abs(z[:14]-x[:14]))<=.001:
            continue
        score=min((x[j]-z[j])/amount for j,amount in shift.items()) if shift else 0.
        choices.append((score,max(abs(z[:14]-x[:14])),-lr,r))
    assert choices,"No second distinct feasible truth within the likelihood threshold"
    b=max(choices,key=lambda t:t[:3])[3]
    truths=[]
    for label,r in (("A_observed_fit",a),("B_feasible_alternate",b)):
        z=row_x(c,r)
        truths.append(dict(scenario=label,origin=r["origin"],rho=rho,sigma=sigma,
             donor_margin=r["donor_margin"],reference_lr=likelihood_ratio(float(r["loss"]),f,k.n),
             maximum_mechanism_move=float(max(abs(z[:14]-x[:14]))),
             **{f"x_{n}":r[f"x_{n}"] for n in c.names},**{n:r[n] for n in c.names[:-1]}))
    write(HERE/"truths.tsv",truths)
    write(HERE/"recovery_targets.tsv",[dict(parameter=r["parameter"],profile_status="profile_pass") for r in wanted])
    print("Frozen scenarios",[(r["scenario"],r["origin"],r["lam_max"],r["buffer_smax"]) for r in truths],flush=True)


def draws():
    k=model()
    c=Config(k)
    result={}
    records=[]
    innovations=[]
    for case_id,t in enumerate(read(HERE/"truths.tsv"),1):
        here=t["scenario"]
        p,_=c.unmap(row_x(c,t))
        base=k.predict_logs(p)
        for i in range(1,101):
            seed=np.random.SeedSequence([20261118,221,case_id,i])
            z=np.random.default_rng(seed).normal(size=k.n)
            innovations.extend(z)
            rho=float(t["rho"])
            sigma=float(t["sigma"])
            noise=sigma*z
            for j in k.transition:
                a=rho**k.gaps[j]
                noise[j]=a*noise[j-1]+np.sqrt(1-a*a)*sigma*z[j]
            y=base+noise
            assert min(y)>0
            result[here,i]=y
            records.append(dict(scenario=here,dataset=i,seed=f"20261118:221:{case_id}:{i}",
                                n_counts=k.n,minimum_log_count=min(y),innovation_rms=float(np.std(z))))
    manifest=HERE/"draw_manifest.tsv"
    if not manifest.exists():
        write(manifest,records)
        write(HERE/"draw_checks.tsv",[dict(datasets=len(result),n_observations=200*221,
                  gaussian_innovation_mean=float(np.mean(innovations)),gaussian_innovation_sd=float(np.std(innovations)))])
    return result


def identify(outcomes,k,incumbent=None):
    feasible=[r for r in outcomes if float(r["donor_margin"])>=0]
    if incumbent and float(incumbent["donor_margin"])>=0:
        feasible.append(incumbent)
    if not feasible:
        return min(outcomes,key=lambda r:float(r["loss"])),False,np.inf
    floor=min(feasible,key=lambda r:float(r["loss"]))
    stable=sorted((r for r in feasible if float(r["stationarity"])<=STATIONARY),key=lambda r:float(r["loss"]))
    gap=(likelihood_ratio(float(stable[1]["loss"]),float(stable[0]["loss"]),k.n)
         if len(stable)>1 else np.inf)
    if gap<=MAX_AGREEMENT_LR and likelihood_ratio(float(stable[0]["loss"]),float(floor["loss"]),k.n)<=.001:
        return stable[0],True,gap
    return floor,False,gap


def one_recovery(task):
    case,i,y,bank,previous,phase=task
    k=model()
    c=Config(k)
    rng=np.random.default_rng(zlib.crc32(f"{case}:{i}:{phase}".encode())+89102)

    def jitter(x,s):
        return np.clip(x+rng.normal(0,s,len(x)),0,1)

    trial=[]
    if previous:
        base=row_x(c,previous)
        plans=[(jitter(base,sd),scale,lim) for sd,scale,lim in
               ((.00002,35.,78),(.0002,8.,85),(.001,16.,85),(.004,4.,90))]
        ranked=[]
        for z in bank:
            out,_=k.metrics(c.unmap(z)[0],y)
            ranked.append((out[1]<0,out[0],z))
        ranked.sort(key=lambda entry:entry[:2])
        distinct=[ranked[0][2]]
        for _,_,z in ranked[1:]:
            if all(np.linalg.norm(z-q)>.025 for q in distinct):
                distinct.append(z)
            if len(distinct)>=3:
                break
        plans.extend((jitter(z,.002),13.,105) for z in distinct)
    else:
        ranked=[]
        for z in bank:
            out,_=k.metrics(c.unmap(z)[0],y)
            ranked.append((out[1]<0,out[0],z))
        ranked.sort(key=lambda r:r[:2])
        a=ranked[0][2]
        b=next((r[2] for r in ranked[1:] if np.linalg.norm(r[2]-a)>.018),ranked[1][2])
        plans=[(jitter(a,.0015),14.,54),(jitter(b,.004),6.,59)]
    for z,scale,lim in plans:
        fit=solve(k,c,y,z,limit=lim,divisor=scale)
        trial.append(receipt(k,c,y,fit))
        best,ok,gap=identify(trial,k,previous)
        if ok:
            break
    if not ok and not previous:
        a=row_x(c,best)
        for sd,scale in ((.00003,27.),(.0007,9.)):
            fit=solve(k,c,y,jitter(a,sd),limit=64,divisor=scale)
            trial.append(receipt(k,c,y,fit))
            best,ok,gap=identify(trial,k)
            if ok:
                break
    exclude={"scenario","dataset","qualified","independent_lr_gap","attempts","phase"}
    return dict(scenario=case,dataset=i,qualified=ok,independent_lr_gap=gap,
                attempts=len(trial),phase=phase,**{q:v for q,v in best.items() if q not in exclude})


def combined():
    groups=defaultdict(list)
    for file in sorted(HERE.glob("recovery_batch_*.tsv")):
        for r in read(file):
            groups[r["scenario"],int(r["dataset"])].append(r)
    answer=[]
    for key,rows in sorted(groups.items()):
        feasible=[r for r in rows if float(r["donor_margin"])>=0]
        floor=min(feasible or rows,key=lambda r:float(r["loss"]))
        steady=[r for r in feasible if r["qualified"]=="True"
                and likelihood_ratio(float(r["loss"]),float(floor["loss"]),221)<=MAX_AGREEMENT_LR]
        answer.append(min(steady,key=lambda r:float(r["loss"])) if steady else floor)
    if answer:
        write(HERE/"qualified_recovery.tsv",answer)
    return answer


def pilot_gate(rows,release=False):
    report=[]
    for key in ("A_observed_fit","B_feasible_alternate"):
        group=[r for r in rows if r["scenario"]==key and int(r["dataset"])<=20]
        verified=sum(r["qualified"]=="True" for r in group)
        report.append(dict(scenario=key,n=len(group),verified=verified,minimum=PILOT_REQUIRED,
                           passed=len(group)==20 and verified>=PILOT_REQUIRED))
    write(HERE/"pilot_gate.tsv",report)
    passed=all(r["passed"] for r in report)
    if release and passed and not (HERE/"pilot_release.tsv").exists():
        write(HERE/"pilot_release.tsv",report)
    print("pilot",[(r["scenario"],r["verified"]) for r in report],"released",passed,flush=True)
    return passed


def run(phase,workers,version,maximum):
    k=model()
    c=Config(k)
    samples=draws()
    existing=combined()
    assignment={(r["scenario"],int(r["dataset"])):r for r in existing}
    if phase=="expanded" and not pilot_gate(existing,release=True):
        raise ValueError("20-by-2 stability gate not met")
    path=HERE/f"recovery_batch_{phase}_{version}.tsv"
    rows=read(path) if path.exists() else []
    previous={(r["scenario"],int(r["dataset"])) for r in rows}
    stable=[r for r in existing if (phase=="repair" or int(r["dataset"])<=20) and r["qualified"]=="True"]
    archive=old_bank(c)
    a=read(COUNTS/"real_best.tsv")[0]
    x=c.load_old(a,.815)
    blind_pool=[x,c.x_anchor]
    for z in sorted(archive,key=lambda v:np.linalg.norm(v-x)):
        if min(np.linalg.norm(v-z) for v in blind_pool)>.14:
            blind_pool.append(z)
        if len(blind_pool)>=10:
            break
    jobs=[]
    for i in range(1,maximum+1):
        for name in ("A_observed_fit","B_feasible_alternate"):
            key=(name,i)
            old=assignment.get(key)
            if key in previous:
                continue
            if phase=="pilot" and (i>20 or old):
                continue
            if phase=="expanded" and (i<=20 or old):
                continue
            if phase=="repair" and (not old or old["qualified"]=="True"):
                continue
            y=samples[key]
            starts=blind_pool.copy()
            if stable:
                eligible=[r for r in stable if (r["scenario"],int(r["dataset"]))!=key]
                near=sorted(eligible,key=lambda r:np.linalg.norm(samples[r["scenario"],int(r["dataset"])]-y))
                starts=[row_x(c,r) for r in near[:8 if phase=="repair" else 4]]+starts[:5]
            jobs.append((name,i,y,starts,old,f"{phase}:{version}"))
    if jobs:
        with pool(workers) as executor:
            for f in as_completed([executor.submit(one_recovery,t) for t in jobs]):
                r=f.result()
                rows.append(r)
                write(path,rows)
                print("recovery",r["scenario"],r["dataset"],r["qualified"],r["stationarity"],
                      r["independent_lr_gap"],flush=True)
    pilot_gate(combined())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("design","pilot","expanded","repair","combine"))
    parser.add_argument("--workers",type=int,default=8)
    parser.add_argument("--version",default="1")
    parser.add_argument("--max",type=int,choices=(20,100),default=100)
    args=parser.parse_args()
    if args.mode=="design":
        freeze_truths(args.workers)
    elif args.mode=="combine":
        pilot_gate(combined())
    else:
        run(args.mode,args.workers,args.version,args.max)
