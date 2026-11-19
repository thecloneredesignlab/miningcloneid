"""Freeze paired A7K simulations, execute gated Slurm slices, and audit estimates."""
import argparse
from collections import Counter
from concurrent.futures import as_completed
from functools import lru_cache
from itertools import product
from pathlib import Path
import zlib

import numpy as np
from scipy.stats import binomtest

from draws import draw
from tilt_model import BASE,HERE,PREV,case,cap,pool,read,receipt,checked_fit,stable,take,write
from profile_tilts import convert

CONDITIONS=(("strict","relaxed"),("up","strict"),("up","relaxed"),
            ("down","strict"),("down","relaxed"))
KEYS=("mode","scenario","estimator","case_id","replicate")
PILOT=20
REDUCED=100
REPAIR_SOURCE={"repair":"pilot","repair_full":"expanded"}


def key(row):
    return (row["mode"],row["scenario"],row["estimator"],int(row["case_id"]),int(row["replicate"]))


def field_mode(mode,estimator):
    return mode if estimator=="relaxed" else "strict"


@lru_cache
def archived():
    return {(int(r["case_id"]),int(r["replicate"])):r for r in read(PREV/"recovery_stock_a7k.tsv")}


def warm_bank(fit_mode,qualified):
    k,c,y=case(fit_mode)
    if fit_mode=="strict":
        source=read(PREV/"fits_stock_a7k.tsv")+read(PREV/"profile_points_stock_a7k.tsv")
    else:
        source=read(HERE/f"fits_{fit_mode}.tsv")+read(HERE/f"profile_points_{fit_mode}.tsv")
    rows=qualified+sorted(source,key=lambda r:float(r["loss"]))
    bank=[]
    for r in rows:
        if float(r.get("all_constraint_margin",r["donor_margin"]))<0.:
            continue
        z=convert(fit_mode,r)
        if not bank or all(np.linalg.norm(z-prev)>.065 for prev in bank):
            bank.append(z)
        if len(bank)==12:
            break
    return np.asarray(bank)


def archive_path(mode,phase):
    return HERE/f"batch_{mode}_{phase}"


def stored(mode,phase):
    f=HERE/f"recovery_{mode}_{phase}.tsv"
    return read(f) if f.exists() else []


def numerically_completed(rows,scenario,estimator):
    return [sum(r["validated"]=="True" and r["scenario"]==scenario and
                r["estimator"]==estimator and int(r["case_id"])==j and
                int(r["replicate"])<=PILOT for r in rows) for j in (1,2)]


def prepare(mode,phase):
    out=archive_path(mode,phase)
    out.mkdir(exist_ok=True)
    interval=next(r for r in read(HERE/f"intervals_{mode}.tsv") if r["parameter"]=="p_misseg")
    assert interval["precision_gate"]=="profile_pass"
    pilot=stored(mode,"pilot")
    selectors=[]
    repair=phase in REPAIR_SOURCE
    if repair:
        target=[r for r in stored(mode,REPAIR_SOURCE[phase]) if r["validated"]!="True"]
        selectors=[key(r)[1:] for r in target]
    else:
        for scenario,method in CONDITIONS:
            if phase=="expanded" and min(numerically_completed(pilot,scenario,method))<18:
                continue
            limits=(1,PILOT+1) if phase=="pilot" else (PILOT+1,REDUCED+1)
            selectors.extend((scenario,method,j,i)
                             for i in range(*limits) for j in (1,2))
    if not selectors:
        print(0)
        return
    # Interleave model conditions so scheduler shards mix slow and fast fits.
    selectors.sort(key=lambda s:(s[3],s[0],s[1],s[2]))
    qualified=[r for r in pilot if r["validated"]=="True"]
    names=("strict","relaxed")
    shapes={}
    for label in names:
        fm=field_mode(mode,label)
        subset=[r for r in qualified if r["estimator"]==label]
        shapes[label]=warm_bank(fm,subset)
    hist=archived()
    lookup={(r["scenario"],r["estimator"],int(r["case_id"]),int(r["replicate"])):r
            for r in (target if repair else pilot)}
    targets=[]
    warms=[]
    metadata=[]
    for scenario,method,j,i in selectors:
        fm=field_mode(mode,method)
        k,c,_=case(fm)
        yy,z,tv=draw(mode,j,i,scenario)
        strict_y,_,_=draw(mode,j,i)
        assert np.array_equal(yy[:k.n+2*133],strict_y[:k.n+2*133])
        r=lookup[scenario,method,j,i] if repair else hist[j,i]
        start=convert(fm,r)
        targets.append(yy)
        # Pad the two unequal optimizer coordinate lengths for compact staging.
        warms.append(np.pad(start, (0,21-len(start)),constant_values=np.nan))
        metadata.append(dict(mode=mode,scenario=scenario,estimator=method,
                         case_id=j,replicate=i,max_generated_tv=max(tv)))
    write(out/"manifest.tsv",metadata)
    np.savez_compressed(out/"frozen.npz",targets=targets,warm=warms,**shapes)
    print(len(selectors),flush=True)


def worker(args):
    mode,scenario,method,j,i,y,w,bank,repair=args
    fit_mode=field_mode(mode,method)
    k,c,_=case(fit_mode)
    rng=np.random.default_rng(zlib.crc32(f"{mode}:{scenario}:{method}:{j}:{i}:{repair}".encode()))
    anchors=[]
    ranked=[]
    for z in bank:
        p,_=c.unmap(z)
        ar,f,*_=k.parts(p,y)
        m=ar[1]
        if k.cap:
            m=min(m,float(np.min(1.-k.composition(p)/k.cap)))
        ranked.append((m<0.,f+180.*max(0.,-m),z))
    for _,_,z in sorted(ranked,key=lambda z:z[:2]):
        if all(np.linalg.norm(z-v)>.035 for v in anchors):
            anchors.append(z)
        if len(anchors)==3:
            break
    recent=w[~np.isnan(w)].copy()
    if k.cap and not repair:
        i0=c.dim
        # Start away from the non-differentiable all-zero perturbation.
        recent[i0:]=np.clip(recent[i0:]+rng.normal(0.,.00042,len(recent)-i0),0,1)
    starts=[recent,*anchors]
    results=[]
    gap=np.inf
    good=False
    for t in range(8 if repair else 7):
        if t<2 and starts:
            start=starts.pop(0)
            sd=.00036
        elif results and (t%2==0 or not starts):
            start=c.anchor(take(results))
            sd=.00009*(t//2+1)
        else:
            start=starts.pop(0)
            sd=.0009
        x=np.clip(start+rng.normal(0,sd,len(start)),0,1)
        fit=checked_fit(k,c,y,x,limit=(125 if repair else 86),divisor=(25.,45.,13.,32.)[t%4])
        results.append(receipt(k,c,y,fit))
        selected=take(results)
        valid=sorted((r for r in results if stable(r)),key=lambda r:float(r["loss"]))
        if len(valid)>1:
            gap=valid[1]["loss"]-valid[0]["loss"]
            viable=min(r["loss"] for r in results if r["all_constraint_margin"]>=0.)
            good=gap<=.10 and valid[0]["loss"]<=viable+.005
            if good:
                selected=valid[0]
                break
    placeholders={name:"" for name in [f"x_tilt_{coef}_{g}"
                             for g in ("2N","4N") for coef in ("mean","shape")]}
    placeholders.update(selected)
    return dict(mode=mode,scenario=scenario,estimator=method,case_id=j,replicate=i,
                validated=good,relative_lr_gap=gap,attempts=len(results),**placeholders)


def execute(mode,phase,shard,width,cores):
    path=archive_path(mode,phase)
    manifest=read(path/"manifest.tsv")
    frozen=np.load(path/"frozen.npz")
    arrays={name:frozen[name] for name in ("targets","warm","strict","relaxed")}
    data=[]
    tasks=[]
    for i in range(shard*width,min((shard+1)*width,len(manifest))):
        r=manifest[i]
        j,d=int(r["case_id"]),int(r["replicate"])
        tasks.append((mode,r["scenario"],r["estimator"],j,d,arrays["targets"][i],
                      arrays["warm"][i],arrays[r["estimator"]],phase in REPAIR_SOURCE))
    with pool(min(cores,len(tasks))) as ex:
        for fit in as_completed(ex.submit(worker,u) for u in tasks):
            r=fit.result(); data.append(r)
            write(path/f"part_{shard:04d}.tsv",data)
            print(key(r),r["validated"],r["attempts"],flush=True)


def collect(mode,phase):
    path=archive_path(mode,phase)
    manifest=read(path/"manifest.tsv")
    received=[r for p in sorted(path.glob("part_*.tsv")) for r in read(p)]
    assert len(manifest)==len(received)
    assert {key(r) for r in manifest}=={key(r) for r in received}
    destination=REPAIR_SOURCE.get(phase,phase)
    existing=stored(mode,destination) if phase in REPAIR_SOURCE else []
    merged={key(r):r for r in existing}
    checks=[]
    for r in received:
        m,scenario,method,j,i=key(r)
        k,c,_=case(field_mode(m,method))
        y,z,tv=draw(m,j,i,scenario)
        p,_=c.unmap(c.anchor(r))
        ar,f,*_=k.parts(p,y)
        diff=abs(f-float(r["loss"]))
        assert diff<1e-7,(key(r),diff)
        if r["validated"]=="True":
            assert ar[1]>=-4e-11
            if k.cap:
                assert max(k.composition(p))<=k.cap+1e-11
        item=key(r)
        before=merged.get(item)
        if before is None or (r["validated"]=="True" and before["validated"]!="True") or (
            r["validated"]==before["validated"] and float(r["loss"])<float(before["loss"])):
            merged[item]=r
        checks.append(dict(zip(KEYS,item),deviance_replay_error=diff,generated_tv=max(tv),
                  estimated_tv=max(k.composition(p)),margin=ar[1],validated=r["validated"]))
    result=[merged[x] for x in sorted(merged)]
    write(HERE/f"recovery_{mode}_{destination}.tsv",result)
    write(HERE/f"replay_{mode}_{phase}.tsv",checks)
    for scenario,method in CONDITIONS:
        if destination=="pilot":
            print(mode,scenario,method,numerically_completed(result,scenario,method),flush=True)
    print(mode,phase,len(checks),"locally replayed",flush=True)


def finish(mode):
    pilot=stored(mode,"pilot")
    expanded=stored(mode,"expanded")
    data=pilot+expanded
    hist=archived()
    table=[]
    paired=[]
    baseline={}
    truth=read(BASE/"truths.tsv")
    for scenario,method in CONDITIONS:
        for j in (1,2):
            x0=float(truth[j-1]["x_p_misseg"])
            before=[abs(float(hist[j,s]["x_p_misseg"])-x0)<=.10 for s in range(1,101)]
            baseline[j]=before
            rs={int(r["replicate"]):r for r in data if r["scenario"]==scenario and
                 r["estimator"]==method and int(r["case_id"])==j}
            if not rs:
                continue
            n=len(rs)
            ok=lambda v: v["validated"]=="True" and abs(float(v["x_p_misseg"])-x0)<=.10
            successes=sum(ok(r) for r in rs.values())
            ci=binomtest(successes,n).proportion_ci()
            table.append(dict(mode=mode,scenario=scenario,estimator=method,case_id=j,
                allocated=n,validated=sum(r["validated"]=="True" for r in rs.values()),
                recovered=successes,rate=successes/n,ci_lo=ci.low,ci_hi=ci.high,
                pilot_validated=numerically_completed(pilot,scenario,method)[j-1]))
            if method!="relaxed":
                continue
            cmp=([bool(before[i-1]) for i in sorted(rs)] if scenario=="strict" else
                 [ok(r) for r in sorted((r for r in data if r["scenario"]==scenario and
                       r["estimator"]=="strict" and int(r["case_id"])==j and
                       int(r["replicate"]) in rs),key=lambda r:int(r["replicate"]))])
            now=[ok(rs[i]) for i in sorted(rs)]
            if len(cmp)==n:
                counts=Counter(zip(cmp,now))
                paired.append(dict(mode=mode,scenario=scenario,case_id=j,paired=n,
                    strict_hits=sum(cmp),relaxed_hits=sum(now),
                    strict_only=counts[True,False],relaxed_only=counts[False,True],
                    rate_delta=(sum(now)-sum(cmp))/n))
    write(HERE/f"recovery_metrics_{mode}.tsv",table)
    if paired:
        write(HERE/f"paired_recovery_{mode}.tsv",paired)
    print(table,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=["prepare","run","collect","finish"])
    p.add_argument("mode")
    p.add_argument("phase",nargs="?",choices=["pilot","expanded",*REPAIR_SOURCE])
    p.add_argument("--shard",type=int)
    p.add_argument("--width",type=int,default=8)
    p.add_argument("--cores",type=int,default=4)
    a=p.parse_args()
    if a.action=="run":
        execute(a.mode,a.phase,a.shard,a.width,a.cores)
    elif a.action=="finish":
        finish(a.mode)
    else:
        (prepare if a.action=="prepare" else collect)(a.mode,a.phase)
