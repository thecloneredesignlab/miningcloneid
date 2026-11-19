"""Paired directional errors, fixed-truth refits and branch-specific exchanges."""
import argparse
from collections import defaultdict
from concurrent.futures import as_completed
from functools import lru_cache
from pathlib import Path
import sys
import zlib

import numpy as np
from scipy.stats import binomtest

from branch_influence import ARCH, BASE, EPS, NAMES, OUT, case, checked_fit, pool, read, receipt, stable, take, write
from draws import coefficients, draw

MODES=("group_025","group_050")
SCENARIOS=("strict","up","down")
ORACLES=("tilt","buffer_beta","k_o_mis")
KEYS=("mode","case_id","replicate","scenario","condition","branch","method","offset")


@lru_cache
def archived(mode):
    rows=[]
    for phase in ("pilot","expanded"):
        rows.extend(read(ARCH/f"recovery_{mode}_{phase}.tsv"))
    assert len(rows)==1000 and all(r["validated"]=="True" for r in rows)
    indexed={(int(r["case_id"]),int(r["replicate"]),r["scenario"],r["estimator"]):r for r in rows}
    rechecks=OUT/"unrestricted_rechecks.tsv"
    if rechecks.exists():
        for r in read(rechecks):
            if r["mode"]==mode and r["validated"]=="True":
                key=(int(r["case_id"]),int(r["replicate"]),r["scenario"],"relaxed")
                assert float(r["loss"]) < float(indexed[key]["loss"])
                indexed[key]=r
    return indexed


@lru_cache
def true_x(mode,j):
    return case(mode)[1].anchor(read(BASE/"truths.tsv")[j-1])


def score(mode,j,row):
    return float(row["x_p_misseg"])-true_x(mode,j)[case(mode)[1].names.index("p_misseg")]


def paired_stats(a,b,seed):
    a,b=np.asarray(a),np.asarray(b)
    first=abs(a)<=.1
    second=abs(b)<=.1
    win=int(np.sum(~first&second))
    loss=int(np.sum(first&~second))
    rg=np.random.default_rng(seed)
    take=rg.integers(len(a),size=(4000,len(a)))
    difference=second.astype(int)-first
    interval=np.quantile(difference[take].mean(axis=1),(.025,.975))
    return dict(n=len(a),recovered_first=int(first.sum()),recovered_second=int(second.sum()),
                second_only=win,first_only=loss,rate_delta=(win-loss)/len(a),
                delta_lower=float(interval[0]),delta_upper=float(interval[1]),
                paired_exact_p=binomtest(win,win+loss).pvalue if win+loss else 1.,
                bias_first=float(a.mean()),bias_second=float(b.mean()),
                sd_first=float(a.std(ddof=1)),sd_second=float(b.std(ddof=1)),
                below_first=int(np.sum(a<-.1)),above_first=int(np.sum(a>.1)),
                below_second=int(np.sum(b<-.1)),above_second=int(np.sum(b>.1)),
                mae_first=float(abs(a).mean()),mae_second=float(abs(b).mean()))


def baseline():
    out=[]
    select=[]
    for mode in MODES:
        rows=archived(mode)
        for j in (1,2):
            for method in ("strict","relaxed"):
                for first,second in (("strict","down"),("up","down"),("strict","up")):
                    if method=="strict" and "strict" in (first,second):
                        continue
                    a,b=([score(mode,j,rows[j,i,s,method]) for i in range(1,101)]
                         for s in (first,second))
                    out.append(dict(mode=mode,case_id=j,method=method,first=first,second=second,
                                  **paired_stats(a,b,zlib.crc32(f"{mode}-{j}-{method}-{first}-{second}".encode()))))
            rg=np.random.default_rng(zlib.crc32(f"paired-panel-{mode}-{j}".encode()))
            strata=defaultdict(list)
            for i in range(1,101):
                a,b=[abs(score(mode,j,rows[j,i,s,"relaxed"]))<=.1 for s in ("up","down")]
                strata[("both" if a&b else "neither" if not a and not b else "up_only" if a else "down_only")].append(i)
            chosen=set()
            for s,n in (("down_only",5),("up_only",2),("both",3),("neither",2)):
                for i in rg.permutation(strata[s])[:n]:
                    chosen.add(int(i))
                    select.append(dict(mode=mode,case_id=j,replicate=int(i),stratum=s))
            for i in rg.permutation([i for i in range(1,101) if i not in chosen])[:12-len(chosen)]:
                select.append(dict(mode=mode,case_id=j,replicate=int(i),stratum="fill"))
    write(OUT/"directional_paired_errors.tsv",out)
    write(OUT/"diagnostic_panel.tsv",select)
    print(len(out),"paired summaries and",len(select),"selected paired replicates.",flush=True)


@lru_cache
def sample(mode,j,i,s):
    return draw(mode,j,i,s)[0]


@lru_cache
def oracle_fixed(mode,j,s,condition):
    k,c,_=case(mode)
    truth=true_x(mode,j)
    if condition in ("buffer_beta","k_o_mis"):
        ix=c.names.index(condition)
        return {ix:float(truth[ix])}
    if condition=="tilt":
        p,_=c.unmap(truth)
        _,before,_,_=k.composition(p,True)
        shifts=[0.,0.]
        if s!="strict":
            for g in range(2):
                shifts[g]=coefficients(before,np.array([2*g,2*g+1]),k.cap,
                                       1. if s=="up" else -1.)
        ans={}
        for j,group in enumerate(("2N","4N")):
            ans[c.names.index(f"tilt_mean_{group}")]=.5*(1.+shifts[j]/c.domain)
            ans[c.names.index(f"tilt_shape_{group}")]=.5
        z=truth.copy()
        for idx,v in ans.items():
            z[idx]=v
        assert max(k.composition(c.unmap(z)[0])) <= k.cap+1e-12
        return ans
    return {}


def row_key(row):
    return (row["mode"],int(row["case_id"]),int(row["replicate"]),row["scenario"],
            row["condition"],row["branch"],row["method"],float(row["offset"]))


def task(mode,j,i,s,condition,branch="-",method="relaxed",offset=9.):
    return dict(zip(KEYS,(mode,j,i,s,condition,branch,method,offset)))


def batch(phase):
    return OUT/f"diag_{phase}"


def build_rows(phase,branches):
    panel=read(OUT/"diagnostic_panel.tsv")
    if phase=="oracle":
        return [task(r["mode"],int(r["case_id"]),int(r["replicate"]),s,oracle)
                for r in panel for s in SCENARIOS for oracle in ORACLES]
    if phase=="shape":
        used=list(dict.fromkeys((r["mode"],int(r["case_id"]),r["stratum"]) for r in panel))
        chosen=[next(r for r in panel if (r["mode"],int(r["case_id"]),r["stratum"])==u)
                for u in used if u[0]=="group_050" and u[2] in ("down_only","both")]
        return [task(r["mode"],int(r["case_id"]),int(r["replicate"]),s,o,offset=offset)
                for r in chosen for s in ("up","down") for o in (*ORACLES,"free")
                for offset in (-.10,-.05,0.,.05,.10)]
    if phase=="swap_pilot":
        return [task("group_050",j,i,s,"swap",branch)
                for i in range(1,21) for j in (1,2)
                for s in ("up","down") for branch in branches]
    if phase=="swap_full":
        pilot=read(results("swap_pilot"))
        effects=defaultdict(list)
        for r in pilot:
            effects[r["branch"],r["scenario"],int(r["case_id"])].append(r)
        passed=select_expansion(effects)
        return [task("group_050",j,i,s,"swap",branch)
                for i in range(21,101) for j in (1,2) for branch,s in passed]
    raise ValueError(phase)


def target(task_row):
    mode,j,i,s,condition,branch,method,offset=row_key(task_row)
    k,c,_=case(mode)
    y=sample(mode,j,i,s)
    ix=c.names.index("p_misseg")
    true=true_x(mode,j)
    fixed=oracle_fixed(mode,j,s,condition).copy()
    if offset!=9.:
        fixed[ix]=float(true[ix]+offset)
    if branch!="-":
        opposite=sample(mode,j,i,"up" if s=="down" else "down")
        assert np.array_equal(y[:k.n+266],opposite[:k.n+266])
        t=k.n+(2+NAMES.index(branch))*133
        ans=y.copy()
        ans[t:t+133]=opposite[t:t+133]
        assert ans[t:t+133].sum()==y[t:t+133].sum()
        assert np.array_equal(ans[:t],y[:t]) and np.array_equal(ans[t+133:],y[t+133:])
        y=ans
    start=c.anchor(archived(mode)[j,i,s,method])
    extra=true.copy()
    for pos,v in oracle_fixed(mode,j,s,"tilt").items():
        extra[pos]=v
    for pos,v in fixed.items():
        extra[pos]=start[pos]=v
    fixed_row=np.full(len(start),np.nan)
    for pos,v in fixed.items():
        fixed_row[pos]=v
    return y,start,extra,fixed_row


def prepare(phase,branches):
    jobs=build_rows(phase,branches)
    if not jobs:
        print(0)
        return
    folder=batch(phase)
    folder.mkdir(exist_ok=True)
    write(folder/"manifest.tsv",jobs)
    y,x,z,f=map(np.stack,zip(*(target(r) for r in jobs)))
    np.savez_compressed(folder/"frozen.npz",target=y,start=x,alternative=z,fixed=f)
    print(len(jobs),flush=True)


def work(args):
    meta,y,start,extra,fixed_row=args
    mode,j,i,s,condition,branch,method,offset=row_key(meta)
    k,c,_=case(mode)
    fixed={int(i):float(fixed_row[i]) for i in np.flatnonzero(np.isfinite(fixed_row))}
    rng=np.random.default_rng(zlib.crc32(repr(row_key(meta)).encode()))
    all_fits=[]
    good=False
    gap=np.inf
    for run in range(7):
        if run==0:
            anchor=start
        elif run==1:
            anchor=extra
        else:
            anchor=c.anchor(take(all_fits))
            anchor=np.clip(anchor+rng.normal(0.,.0002,len(anchor)),0.,1.)
        fit=checked_fit(k,c,y,anchor,fixed=fixed,limit=90 if run<3 else 150,
                        divisor=(32.,24.,18.,42.)[run%4])
        all_fits.append(receipt(k,c,y,fit))
        viable=[r for r in all_fits if float(r["all_constraint_margin"])>=0.]
        checked=sorted([r for r in viable if stable(r)],key=lambda r:r["loss"])
        if len(checked)>=2:
            gap=float(checked[1]["loss"])-float(checked[0]["loss"])
            good=gap<=.1 and checked[0]["loss"]<=min(r["loss"] for r in viable)+.005
            if good:
                break
    best=checked[0] if good else take(all_fits)
    return dict(meta,validated=good,relative_lr_gap=gap,attempts=len(all_fits),**best)


def run(phase,index,width,cores):
    folder=batch(phase)
    manifest=read(folder/"manifest.tsv")
    frozen=np.load(folder/"frozen.npz")
    arrays=[frozen[name] for name in ("target","start","alternative","fixed")]
    selected=range(index*width,min((index+1)*width,len(manifest)))
    rows=[]
    with pool(min(cores,len(selected))) as ex:
        for fn in as_completed(ex.submit(work,(manifest[i],*(v[i] for v in arrays))) for i in selected):
            r=fn.result()
            rows.append(r)
            write(folder/f"part_{index:04d}.tsv",rows)
            print(row_key(r),r["validated"],r["attempts"],flush=True)


def results(phase):
    return OUT/f"diagnostic_{phase}.tsv"


def collect(phase):
    folder=batch(phase)
    manifest=read(folder/"manifest.tsv")
    parts=[r for file in sorted(folder.glob("part_*.tsv")) for r in read(file)]
    assert len(manifest)==len(parts) and {row_key(r) for r in parts}=={row_key(r) for r in manifest}
    yo=np.load(folder/"frozen.npz")["target"]
    targets={row_key(r):v for r,v in zip(manifest,yo)}
    output=[]
    max_error=0.
    for r in sorted(parts,key=row_key):
        mode,j,i,s,condition,branch,method,offset=row_key(r)
        k,c,_=case(mode)
        p,_=c.unmap(c.anchor(r))
        y=targets[row_key(r)]
        ar,f,*_=k.parts(p,y)
        err=abs(f-float(r["loss"]))
        max_error=max(max_error,err)
        assert err<1e-7
        if r["validated"]=="True":
            assert ar[1]>=-1e-10 and max(k.composition(p))<=k.cap+1e-10
        delta=score(mode,j,r)
        output.append(dict(r,delta_x_p_misseg=delta,recovered=r["validated"]=="True" and abs(delta)<=.1,
                           deviance_replay_error=err))
    write(results(phase),output)
    print(phase,len(output),"fits",sum(r["validated"]=="True" for r in output),"validated",
          "max deviation",max_error,flush=True)


def select_expansion(effects):
    passed=[]
    for branch in sorted({k0[0] for k0 in effects}):
        for s in ("up","down"):
            rows=[effects[branch,s,j] for j in (1,2)]
            mode="group_050"
            if any(len(v)!=20 or sum(r["validated"]=="True" for r in v)<18 for v in rows):
                continue
            recovery=[sum(r["recovered"]=="True" for r in v) for v in rows]
            paired=[sum(abs(score(mode,j,archived(mode)[j,int(r["replicate"]),s,"relaxed"]))<=.1
                        for r in v) for j,v in enumerate(rows,1)]
            # Expand if relying on the pilot could change the two-truth precision call.
            close_or_changed=(min(recovery)>=17 or
                              any(abs(a-b)>=2 for a,b in zip(recovery,paired)))
            print(branch,s,"pilot",paired,recovery,"expand",close_or_changed,file=sys.stderr)
            if close_or_changed:
                passed.append((branch,s))
    return passed


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=("baseline","prepare","run","collect"))
    p.add_argument("phase",nargs="?",choices=("oracle","shape","swap_pilot","swap_full"))
    p.add_argument("--branches",nargs="+",choices=NAMES,default=[])
    p.add_argument("--shard",type=int)
    p.add_argument("--width",type=int,default=8)
    p.add_argument("--cores",type=int,default=4)
    a=p.parse_args()
    if a.action=="baseline":
        baseline()
    elif a.action=="prepare":
        prepare(a.phase,a.branches)
    elif a.action=="run":
        run(a.phase,a.shard,a.width,a.cores)
    else:
        collect(a.phase)
