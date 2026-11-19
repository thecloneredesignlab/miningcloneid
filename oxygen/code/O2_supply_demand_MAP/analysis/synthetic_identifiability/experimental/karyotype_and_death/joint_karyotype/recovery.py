"""Apply the registered profile gate and run paired two-truth recoveries."""
import argparse
from concurrent.futures import as_completed
import zlib

import numpy as np

from joint_model import (BASE, HERE, fit_receipt, likelihood_ratio, pool, read,
                         row_x, solve, stable, write)
from refit import case

MODE = "autosomes_joint"
AGREEMENT_LR = .10


def generate(case_id, dataset):
    k,c,_ = case(MODE)
    truth = read(BASE/"truths.tsv")[case_id-1]
    x = row_x(c,truth)
    p = c.unmap(x)[0]
    assert k.supply(p)[0].min() > 0.
    probs, _ = k.prediction(p)
    rho,sigma = float(truth["rho"]),float(truth["sigma"])
    count_rng = np.random.default_rng(np.random.SeedSequence([20261118,221,case_id,dataset]))
    kary_rng = np.random.default_rng(np.random.SeedSequence([20261119,40,case_id,dataset]))
    z = count_rng.normal(size=k.n)
    noise = sigma*z
    for j in k.transition:
        a=rho**k.gaps[j]
        noise[j]=a*noise[j-1]+np.sqrt(1.-a*a)*sigma*z[j]
    samples=np.stack([kary_rng.multinomial(20,q) for q in probs])
    return truth,x,k.target(k.predict_logs(p)+noise,samples)


def worker(job):
    case_id,dataset,bank,prior=job
    k,c,_=case(MODE)
    truth, x, y = generate(case_id,dataset)
    scenario=truth["scenario"]
    rng=np.random.default_rng(zlib.crc32(f"stage3:{scenario}:{dataset}".encode()))
    ranked=[]
    for v in bank:
        out,_,f,_=k.evaluate_joint(c.unmap(v)[0],y)
        ranked.append((out[1]<0.,f,v))
    ranked.sort(key=lambda z:z[:2])
    independent=[]
    for _,_,v in ranked:
        if all(np.linalg.norm(v-a)>.025 for a in independent):
            independent.append(v)
        if len(independent)==4:
            break
    trial=[]
    plans=[(v,.0015,s) for v,s in zip(independent,(32.,16.,24.,12.))]
    if prior:
        plans.insert(0,(row_x(c,prior),.0002,12.))
    verified=False
    gap=np.inf
    for step in range(len(plans)+3):
        if step < len(plans):
            anchor,jitter,divisor=plans[step]
        else:
            anchor,jitter,divisor=row_x(c,best),.0001,14.
        v=np.clip(anchor+rng.normal(0.,jitter,len(anchor)),0.,1.)
        f=solve(k,c,y,v,limit=85,divisor=divisor)
        trial.append(fit_receipt(k,c,y,f))
        allowed=[r for r in trial if r["donor_margin"]>=0.]
        best=min(allowed or trial,key=lambda r:r["loss"])
        settled=sorted((r for r in trial if stable(r)),key=lambda r:r["loss"])
        if len(settled)>1:
            gap=likelihood_ratio(settled[1]["loss"],settled[0]["loss"],k.n)
            verified=gap<=AGREEMENT_LR and likelihood_ratio(settled[0]["loss"],best["loss"],k.n)<=.001
            if verified:
                best=settled[0]
                break
    return dict(scenario=scenario,case_id=case_id,dataset=dataset,qualified=verified,
                independent_lr_gap=gap,attempts=len(trial),**best)


def saved_rows():
    path=HERE/"recovery_results.tsv"
    return read(path) if path.exists() else []


def target_names():
    completed = HERE/f"profile_intervals_{MODE}.tsv"
    source = completed if completed.exists() else HERE/"pilot_profile_gate.tsv"
    return [r["parameter"] for r in read(source)
            if r["profile_status"]=="profile_pass"]


def run(phase,workers):
    k,c,_=case(MODE)
    names=target_names()
    if not names:
        write(HERE/"recovery_gate.tsv",[dict(phase="profile",eligible_targets=0,
                                             status="no_profile_pass")])
        print("No parameter met the 10% profile requirement.",flush=True)
        return
    write(HERE/"recovery_gate.tsv",[dict(phase="profile",eligible_targets=len(names),
                                        status=",".join(names))])
    rows={(int(r["case_id"]),int(r["dataset"])):r for r in saved_rows()}
    truth=[row_x(c,r) for r in read(BASE/"truths.tsv")]
    bank=[]
    sources = [r for mode in (MODE,"autosomes_count","archived_joint")
               for r in read(HERE/f"refits_{mode}.tsv")]
    if phase=="expand":
        sources += [r for (_,i),r in rows.items() if i<=20 and r["qualified"]=="True"]
    for r in sources:
        v=row_x(c,r)
        if min(np.linalg.norm(v-a) for a in truth) > 1e-5 and (
            not bank or min(np.linalg.norm(v-a) for a in bank) > .05
        ):
            bank.append(v)
    pilot=[sum(r["qualified"]=="True" for (case_id,i),r in rows.items()
               if case_id==j and i<=20) for j in (1,2)]
    if phase=="expand" and min(pilot)<18:
        raise RuntimeError(f"Pilot stability gate is still closed: {pilot}")
    last=20 if phase in ("pilot","repair_pilot") else 100
    retry=phase.startswith("repair")
    jobs=[(j,i,bank,rows.get((j,i))) for i in range(1,last+1)
          for j in (1,2) if (j,i) not in rows or (retry and rows[j,i]["qualified"]!="True")]
    with pool(workers) as executor:
        for f in as_completed(executor.submit(worker,t) for t in jobs):
            r=f.result()
            rows[r["case_id"],r["dataset"]]=r
            write(HERE/"recovery_results.tsv",[rows[key] for key in sorted(rows)])
            print(r["scenario"],r["dataset"],r["qualified"],r["loss"],r["attempts"],flush=True)
    rows=saved_rows()
    pilot=[sum(r["qualified"]=="True" and int(r["dataset"])<=20 and
               int(r["case_id"])==j for r in rows) for j in (1,2)]
    write(HERE/"recovery_pilot_gate.tsv",[
        dict(case_id=j,passed=pilot[j-1]>=18,numerically_verified=pilot[j-1],assigned=20)
        for j in (1,2)])
    summary=[]
    for name in names:
        j=c.names.index(name)
        for case_id in (1,2):
            fits=[r for r in rows if int(r["case_id"])==case_id]
            good=[r for r in fits if r["qualified"]=="True"]
            recovered=sum(abs(float(r[f"x_{name}"])-truth[case_id-1][j])<=.1 for r in good)
            summary.append(dict(parameter=name,case_id=case_id,assigned=len(fits),
                                numerically_verified=len(good),recovered=recovered,
                                rate=recovered/len(fits)))
    write(HERE/"recovery_metrics.tsv",summary)
    print("pilot numerical verification",pilot,flush=True)
    print(summary,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase",choices=["pilot","expand","repair_pilot","repair_expanded"])
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    run(a.phase,a.workers)
