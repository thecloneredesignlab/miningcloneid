"""Gated 20/100 paired recovery for each incremental observation mode."""
import argparse
from concurrent.futures import as_completed
import zlib

import numpy as np

from modal import BASE,HERE,Config,case,fixed_for,pool,read,receipt,solve,stable,write
from fit import take

GATE=18


def draw(mode,case_id,replicate):
    k,c,y=case(mode)
    truth=read(BASE/"truths.tsv")[case_id-1]
    x=c.anchor(truth)
    if k.use_death:
        observed=c.anchor(read(HERE/f"best_{mode}.tsv")[0])
        x[-4:]=observed[-4:]
    p,_=c.unmap(x)
    assert k.supply(p)[0].min()>=0.
    seed=np.random.SeedSequence([20261118,221,case_id,replicate])
    rng=np.random.default_rng(seed)
    z=rng.normal(size=k.n)
    rho=float(truth["rho"])
    sigma=float(truth["sigma"])
    noise=sigma*z
    for j in k.transition:
        a=rho**k.gaps[j]
        noise[j]=a*noise[j-1]+np.sqrt(1.-a*a)*sigma*z[j]
    logc=k.predict_logs(p)+noise
    channels,state=k.readouts(p)
    prob=state/state.sum(axis=1,keepdims=True)
    krng=np.random.default_rng(np.random.SeedSequence([20261119,40,case_id,replicate]))
    vals=[krng.multinomial(int(n),q) for n,q in zip(k.kary.sum(axis=1),prob)]
    deaths=None
    if k.use_death:
        q,kappa,*_=k.death_prediction(p,channels,logc)
        d_rng=np.random.default_rng(np.random.SeedSequence([20261119,191,case_id,replicate]))
        latent=d_rng.beta(q*kappa,(1.-q)*kappa)
        deaths=d_rng.binomial(k.eligible.astype(int),latent)
    return k.target(logc,np.stack(vals),deaths),x


def worker(item):
    mode,case_id,replicate,starts,alternate,*override=item
    k,_,_=case(mode)
    c=Config(k)
    y=override[0] if override else draw(mode,case_id,replicate)[0]
    rng=np.random.default_rng(zlib.crc32(f"{mode}:{case_id}:{replicate}:{alternate is not None}".encode()))
    rank=[]
    if k.use_death:
        from scipy.special import logit
        c.share_scale=.0025
        transformed=[]
        for j,z0 in enumerate(starts):
            z=z0.copy()
            z[-2]=c.encode_share(z0[-2])
            transformed.append(z)
            if j<5:
                # Use observed counts and morphology to initialize the nuisance
                # offset and allow a visible hypoxic channel at the start.
                for visible in (.12,.50):
                    z=z.copy()
                    z[c.names.index("eta_H_2N")]=visible
                    z[-4]=.5
                    p,_=c.unmap(z)
                    channels,_=k.readouts(p)
                    pr,*_=k.death_prediction(p,channels,y)
                    obs=y[-len(k.eligible):]
                    offset=np.median(logit((obs+.5)/(k.eligible+1.))-logit(pr))
                    z[-4]=np.clip((offset+4.)/8.,0.,1.)
                    transformed.append(z)
        starts=transformed
    for x in starts:
        ar,f,*_=k.parts(c.unmap(x)[0],y)
        rank.append((ar[1]<0.,f+200*max(0.,-ar[1]),x))
    bank=[]
    for _,_,z in sorted(rank,key=lambda r:r[:2]):
        if all(np.linalg.norm(z-x)>.025 for x in bank):
            bank.append(z)
        if len(bank)==4:
            break
    if alternate:
        bank.insert(0,c.anchor(alternate))
    results=[]
    ok=False
    gap=np.inf
    remaining=list(bank)
    for i in range(len(bank)+4):
        if i<2 and remaining:
            start,spread=remaining.pop(0),.0015
        elif results and (i%2==0 or not remaining):
            start,spread=c.anchor(take(results)),.00012*(1+i//2)
        elif remaining:
            start,spread=remaining.pop(0),.0015
        else:
            start,spread=c.anchor(take(results)),.0002
        x=np.clip(start+rng.normal(0,spread,len(start)),0,1)
        f=solve(k,c,y,x,limit=88,divisor=(45.,24.,62.,31.)[i%4],fixed=fixed_for(mode,c))
        r=receipt(k,c,y,f)
        results.append(r)
        good=sorted((r for r in results if stable(r)),key=lambda r:r["loss"])
        best=take(results)
        if len(good)>1:
            gap=good[1]["loss"]-good[0]["loss"]
            ok=gap<=.1 and good[0]["loss"]<=min(r["loss"] for r in results
                                               if r["donor_margin"]>=0.)+.005
            if ok:
                best=good[0]
                break
    return dict(mode=mode,case_id=case_id,replicate=replicate,validated=ok,
                lr_agreement=gap,attempts=len(results),**best)


def run(mode,phase,workers):
    k,c,y=case(mode)
    names=[r["parameter"] for r in read(HERE/f"intervals_{mode}.tsv")
           if r["precision_gate"]=="profile_pass"]
    if not names:
        write(HERE/f"gate_{mode}.tsv",[dict(mode=mode,eligible="",status="no_profile_pass")])
        print(mode,"none passed profile gate",flush=True)
        return
    write(HERE/f"gate_{mode}.tsv",[dict(mode=mode,eligible=",".join(names),status="profile_pass")])
    path=HERE/f"recovery_{mode}.tsv"
    rows={(int(r["case_id"]),int(r["replicate"])):r for r in read(path)} if path.exists() else {}
    prior=[sum(r["validated"]=="True" for (j,i),r in rows.items() if j==q and i<=20) for q in (1,2)]
    if phase=="expand" and min(prior)<GATE:
        raise ValueError(f"recovery numerical pilot gate closed for {mode}: {prior}")
    sources=read(HERE/f"fits_{mode}.tsv")+read(HERE/f"profile_points_{mode}.tsv")
    if phase=="expand":
        sources += [r for (_,i),r in rows.items() if i<=20 and r["validated"]=="True"]
    starter=[]
    for r in sorted(sources,key=lambda r:float(r["loss"])):
        if not stable(r):
            continue
        x=c.anchor(r)
        if all(np.linalg.norm(x-v)>.09 for v in starter):
            starter.append(x)
        if len(starter)>=20:
            break
    last=100 if phase in ("expand","repair_full") else 20
    items=[(mode,j,i,starter,rows.get((j,i))) for i in range(1,last+1) for j in (1,2)
           if (j,i) not in rows or (phase.startswith("repair") and rows[j,i]["validated"]!="True")]
    with pool(workers) as ex:
        for f in as_completed(ex.submit(worker,a) for a in items):
            r=f.result()
            j,i=r["case_id"],r["replicate"]
            rows[j,i]=r
            write(path,[rows[key] for key in sorted(rows)])
            print(mode,j,i,r["validated"],r["attempts"],flush=True)
    rows=read(path)
    report=[]
    for j in (1,2):
        subset=[r for r in rows if int(r["case_id"])==j]
        _,truth=draw(mode,j,1)
        for name in names:
            target=c.names.index(name)
            proper=[r for r in subset if r["validated"]=="True"]
            hit=sum(abs(c.anchor(r)[target]-truth[target])<=.1 for r in proper)
            report.append(dict(mode=mode,parameter=name,case_id=j,allocated=len(subset),
                     numerically_valid=len(proper),recovered=hit,rate=hit/len(subset),
                     truth_normalized=truth[target]))
    write(HERE/f"recovery_metrics_{mode}.tsv",report)
    print(report,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode")
    p.add_argument("phase",choices=["pilot","expand","repair_pilot","repair_full"])
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    run(a.mode,a.phase,a.workers)
