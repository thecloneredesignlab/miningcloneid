"""Replace one 4N observation at a time while retaining the matched seed."""
import argparse
from concurrent.futures import as_completed
from functools import lru_cache
from pathlib import Path
import sys
import zlib

import numpy as np
from scipy.stats import binomtest

OUT=Path(__file__).resolve().parent
COHORT=OUT.with_name("20261119_A7K_cohort_swap_01")
sys.path.insert(0,str(COHORT))
from cohort_swap import archived, case, digest, load as cohort_fits, sample, score
from tilt_model import Config, checked_fit, pool, read, receipt, stable, take, write

MODES=("group_050","group_025")
BACKS=("up","down")
MIXED=("LH","HL")
ALL=("HH","LH","HL","LL")
FIELDS=("model","mode","background","case_id","replicate","cell")
VARIANTS=("shared","split4N")


class Split4NConfig(Config):
    """Share the 2N pair and give 4N O1 and 4N O2 their own tilt coefficients."""

    def __init__(self,k):
        super().__init__(k)
        self.groups=np.array([0,0,1,2])
        self.names=self.names[:self.dim]+[
            f"tilt_{term}_{label}" for label in ("2N","4N_O1","4N_O2")
            for term in ("mean","shape")]

    def anchor(self,r):
        r=dict(r)
        for term in ("mean","shape"):
            for line in ("O1","O2"):
                r.setdefault(f"x_tilt_{term}_4N_{line}",r.get(f"x_tilt_{term}_4N",.5))
        return super().anchor(r)


@lru_cache
def model(mode,variant):
    k,shared,y=case(mode)
    assert tuple(shared.kary_labels)==("2N_O1","2N_O2","4N_O1","4N_O2")
    return k,shared if variant=="shared" else Split4NConfig(k),y


def batch(variant,mode,phase):
    return OUT/f"batch_{variant}_{mode}_{phase}"


def destination(variant,mode,phase):
    return OUT/f"fits_{variant}_{mode}_{phase}.tsv"


def stored(variant,mode,phase):
    f=destination(variant,mode,phase)
    return read(f) if f.exists() else []


def key(r):
    return tuple(int(r[q]) if q in ("case_id","replicate") else r[q] for q in FIELDS)


@lru_cache
def bounds(mode):
    out={}
    base=archived(mode)
    hybrids=cohort_fits(mode,"pilot")+cohort_fits(mode,"expanded")
    assert len(hybrids)==400
    mixed={(r["branch"],int(r["case_id"]),int(r["replicate"])):r for r in hybrids}
    for j in (1,2):
        for i in range(1,101):
            out["up",j,i,"HH"]=base[j,i,"up","relaxed"]
            out["up",j,i,"LL"]=mixed["4N",j,i]
            out["down",j,i,"HH"]=mixed["2N",j,i]
            out["down",j,i,"LL"]=base[j,i,"down","relaxed"]
    assert all(r["validated"]=="True" for r in out.values())
    return out


@lru_cache
def primary(mode):
    rows=stored("shared",mode,"pilot")+stored("shared",mode,"expanded")
    return {(r["background"],int(r["case_id"]),int(r["replicate"]),r["cell"]):r for r in rows}


@lru_cache
def data(mode,background,j,i,cell):
    k,c,_=model(mode,"shared")
    assert k.n==221 and k.kary.shape==(6,133)
    up,down=[sample(mode,j,i,s) for s in BACKS]
    cut=k.n+266
    assert np.array_equal(up[:cut],down[:cut])
    h=up[cut:].reshape(4,133)
    l=down[cut:].reshape(4,133)
    n=k.kary[2:].sum(axis=1)
    assert np.array_equal(h.sum(axis=1),n) and np.array_equal(l.sum(axis=1),n)
    prefix=h[:2] if background=="up" else l[:2]
    keep=[l[x] if direction=="L" else h[x] for x,direction in enumerate(cell,2)]
    y=np.r_[up[:cut],prefix.ravel(),np.ravel(keep)]
    # Assemble the same histogram from the reciprocal 4N direction and only overwrite
    # the required complement. Holds for both single-branch replacements and controls.
    reverse=l.copy()
    forward=h.copy()
    for x,direction in enumerate(cell,2):
        if direction=="H":
            reverse[x]=h[x]
        else:
            forward[x]=l[x]
    forward[:2]=reverse[:2]=prefix
    a=np.r_[up[:cut],forward.ravel()]
    b=np.r_[down[:cut],reverse.ravel()]
    assert np.array_equal(a,y) and np.array_equal(b,y)
    return y,dict(target_sha256=digest(y),reciprocal_sha256=digest(b),
                  common_sha256=digest(y[:cut]),cells_4N_O1=int(n[2]),cells_4N_O2=int(n[3]))


def group_verified(rows, background, cell, j):
    bucket=[r for r in rows if (r["background"],r["cell"],int(r["case_id"]))==(background,cell,j)]
    assert sorted(int(r["replicate"]) for r in bucket)==list(range(1,21))
    return sum(r["validated"]=="True" for r in bucket)


def prepare(variant,mode,phase,backgrounds):
    k,c,_=model(mode,variant)
    choices=MIXED if variant=="shared" else ALL
    lo,hi=1,21
    if phase=="expanded":
        pilot=stored(variant,mode,"pilot")
        assert all(group_verified(pilot,b,q,j)>=18
                   for b in backgrounds for q in choices for j in (1,2))
        lo,hi=21,101
    bank=bounds(mode)
    manifests=[]
    targets=[]
    warm=[]
    alternatives=[]
    for i in range(lo,hi):
        for j in (1,2):
            for b in backgrounds:
                ends={q:bank[b,j,i,q] for q in ("HH","LL")}
                for q in choices:
                    y,audit=data(mode,b,j,i,q)
                    old=ends.get(q) or (primary(mode).get((b,j,i,q)) if variant=="split4N" else None)
                    if variant=="split4N":
                        assert old and old["validated"]=="True"
                    # For the split model the shared optimum is an exact nested start.
                    x=c.anchor(old or ends["HH"])
                    z=c.anchor(ends["LL"] if q!="LL" else ends["HH"])
                    if variant=="split4N":
                        for label,direction in zip(("O1","O2"),q):
                            for term in ("mean","shape"):
                                pos=c.names.index(f"tilt_{term}_4N_{label}")
                                z[pos]=float(ends["LL" if direction=="L" else "HH"][f"x_tilt_{term}_4N"])
                    if q in ends:
                        assert abs(k.parts(c.unmap(c.anchor(ends[q]))[0],y)[1]-float(ends[q]["loss"]))<1e-7
                    manifests.append(dict(zip(FIELDS,(variant,mode,b,j,i,q)),**audit,
                        shared_loss=float(old["loss"]) if variant=="split4N" else ""))
                    targets.append(y)
                    warm.append(x)
                    alternatives.append(z)
    fold=batch(variant,mode,phase)
    fold.mkdir(exist_ok=True)
    write(fold/"manifest.tsv",manifests)
    np.savez_compressed(fold/"frozen.npz",target=targets,warm=warm,alternative=alternatives)
    print(len(manifests),flush=True)


def one_fit(job):
    meta,y,x,z=job
    variant,mode,b,j,i,q=key(meta)
    k,c,_=model(mode,variant)
    rng=np.random.default_rng(zlib.crc32(str(key(meta)).encode()))
    rows=[]
    best=None
    gap=np.inf
    for attempt in range(8):
        if attempt==0:
            start=x
        elif attempt==1:
            start=z
        else:
            start=c.anchor(take(rows))
            start=np.clip(start+rng.normal(0,.00017+.00002*attempt,len(start)),0,1)
        f=checked_fit(k,c,y,start,limit=92 if attempt<3 else 155,
                      divisor=(30.,42.,17.,32.)[attempt%4])
        rows.append(receipt(k,c,y,f))
        good=sorted((r for r in rows if stable(r)),key=lambda r:r["loss"])
        if len(good)>1:
            smallest=min(r["loss"] for r in rows if r["all_constraint_margin"]>=0.)
            gap=good[1]["loss"]-good[0]["loss"]
            nested=(variant=="shared" or good[0]["loss"]<=float(meta["shared_loss"])+.005)
            if gap<=.10 and good[0]["loss"]<=smallest+.005 and nested:
                best=good[0]
                break
    return dict(meta,validated=best is not None,attempts=len(rows),restart_deviance_gap=gap,
                **(best or take(rows)))


def run(variant,mode,phase,shard,width,cores):
    fold=batch(variant,mode,phase)
    manifests=read(fold/"manifest.tsv")
    start=shard*width
    end=min(start+width,len(manifests))
    nz=np.load(fold/"frozen.npz")
    arrays=[nz[q][start:end] for q in ("target","warm","alternative")]
    jobs=list(zip(manifests[start:end],*arrays))
    rows=[]
    with pool(min(cores,len(jobs))) as ex:
        for f in as_completed(ex.submit(one_fit,j) for j in jobs):
            r=f.result()
            rows.append(r)
            write(fold/f"part_{shard:04d}.tsv",rows)
            print(key(r),r["validated"],r["attempts"],flush=True)


def verify(variant,mode,phase,repair=False):
    fold=batch(variant,mode,phase)
    manifests=read(fold/"manifest.tsv")
    targets={key(r):y for r,y in zip(manifests,np.load(fold/"frozen.npz")["target"])}
    k,c,_=model(mode,variant)
    if repair:
        rows=stored(variant,mode,phase)
        bad={key(r):r for r in rows if r["validated"]!="True"}
        if not bad:
            return
        jobs=[]
        for r in manifests:
            f=key(r)
            if f not in bad:
                continue
            x=c.anchor(bad[f])
            z=np.clip(x+np.random.default_rng(f[4]+413).normal(0,.00021,len(x)),0,1)
            jobs.append((r,targets[f],x,z))
        updates={}
        with pool(min(8,len(jobs))) as ex:
            for f in as_completed(ex.submit(one_fit,j) for j in jobs):
                r=f.result()
                if r["validated"] and float(r["loss"])<=float(bad[key(r)]["loss"])+.02:
                    updates[key(r)]=r
                print("repair",key(r),key(r) in updates,flush=True)
        rows=[updates.get(key(r),r) for r in rows]
    else:
        rows=[r for f in sorted(fold.glob("part_*.tsv")) for r in read(f)]
    assert len(rows)==len(manifests) and {key(r) for r in rows}==set(targets)
    result=[]
    fixed=np.array([float(r["fixed_recorded_seed"]) for r in k.a7])
    for r in sorted(rows,key=key):
        y=targets[key(r)]
        assert digest(y)==r["target_sha256"]==r["reciprocal_sha256"]
        p=c.unmap(c.anchor(r))[0]
        ar,dev,*_=k.parts(p,y)
        tv,_,_,seed=k.composition(p,True)
        se=float(np.max(abs(seed-fixed)))
        err=abs(dev-float(r["loss"]))
        assert err<1e-7 and se<1e-7
        valid=str(r["validated"])=="True"
        if valid:
            assert ar[1]>=-1e-10 and float(max(tv))<=k.cap+1e-10
            if variant=="split4N":
                assert dev<=float(r["shared_loss"])+.005
        deviation=score(mode,int(r["case_id"]),r)
        result.append(dict(r,normalized_error=deviation,recovered=valid and abs(deviation)<=.1,
                           deviance_replay_error=err,max_seed_error=se))
    write(destination(variant,mode,phase),result)
    checked=sum(str(r["validated"])=="True" for r in result)
    print(variant,mode,phase,len(result),"replayed",checked,"validated",flush=True)


def accuracy(variant,mode):
    fresh=stored(variant,mode,"pilot")+stored(variant,mode,"expanded")
    mapping={(r["background"],int(r["case_id"]),int(r["replicate"]),r["cell"]):r for r in fresh}
    backgrounds=sorted({r["background"] for r in fresh},key=BACKS.index)
    assert len(fresh)==len(mapping)==len(backgrounds)*(2 if variant=="shared" else 4)*200
    original=bounds(mode)
    metrics=[]
    contrasts=[]
    all_rows=[]
    for background in backgrounds:
        for j in (1,2):
            by={q:[mapping[background,j,i,q] if (q in MIXED or variant=="split4N") else
                   original[background,j,i,q] for i in range(1,101)] for q in ALL}
            e={q:np.array([score(mode,j,r) for r in r0]) for q,r0 in by.items()}
            feasible={q:np.array([r["validated"]=="True" for r in r0]) for q,r0 in by.items()}
            recovered={q:(abs(e[q])<=.1)&feasible[q] for q in ALL}
            seed=zlib.crc32(str((variant,mode,background,j)).encode())
            indices=np.random.default_rng(seed).integers(100,size=(8000,100))
            for q in ALL:
                for i,(error,ok,r) in enumerate(zip(e[q],recovered[q],by[q]),1):
                    all_rows.append(dict(model=variant,mode=mode,background=background,case_id=j,replicate=i,
                                         cell=q,normalized_error=error,recovered=bool(ok),
                                         validated=r["validated"],p_misseg=r["p_misseg"]))
                valid=feasible[q]
                passed=int(recovered[q].sum())
                ci=binomtest(passed,100).proportion_ci()
                metrics.append(dict(model=variant,mode=mode,background=background,case_id=j,cell=q,
                      n=100,validated=int(valid.sum()),recovered=passed,rate=passed/100,
                      ci_lo=ci.low,ci_hi=ci.high,bias=float(e[q][valid].mean()),
                      sd=float(e[q][valid].std(ddof=1)),low_fail=int((valid&(e[q]<-.1)).sum()),
                      high_fail=int((valid&(e[q]>.1)).sum())))
            for label,a,b in (("O1_given_O2_high","HH","LH"),("O1_given_O2_low","HL","LL"),
                              ("O2_given_O1_high","HH","HL"),("O2_given_O1_low","LH","LL"),
                              ("O1_minus_O2","HL","LH"),("both","HH","LL"),
                              ("interaction","-","-")):
                rc={q:r.astype(float) for q,r in recovered.items()}
                if label=="interaction":
                    diff=e["LL"]-e["LH"]-e["HL"]+e["HH"]
                    rd=rc["LL"]-rc["LH"]-rc["HL"]+rc["HH"]
                    join=np.logical_and.reduce(list(feasible.values()))
                    improved=reversed_count=""
                else:
                    diff=e[b]-e[a]
                    rd=rc[b]-rc[a]
                    join=feasible[a]&feasible[b]
                    improved,reversed_count=(int((rd==v).sum()) for v in (1,-1))
                ci=np.quantile(rd[indices].mean(1),(.025,.975))
                eci=np.quantile(diff[indices].mean(1),(.025,.975)) if join.all() else (np.nan,np.nan)
                contrasts.append(dict(model=variant,mode=mode,background=background,case_id=j,contrast=label,
                    from_cell=a,to_cell=b,n=100,paired_validated=int(join.sum()),
                    signed_delta=diff[join].mean(),signed_lo=eci[0],signed_hi=eci[1],
                    rate_delta=rd.mean(),rate_lo=ci[0],rate_hi=ci[1],improved=improved,reversed=reversed_count))
    write(OUT/f"metrics_{variant}_{mode}.tsv",metrics)
    write(OUT/f"contrasts_{variant}_{mode}.tsv",contrasts)
    write(OUT/f"paired_{variant}_{mode}.tsv",all_rows)
    for bg in backgrounds:
        for cell in ALL:
            a,b=(int(r["recovered"]) for r in metrics if r["background"]==bg and r["cell"]==cell)
            print(variant,mode,bg,cell,a,b,min(a,b)>=90,flush=True)


def decide():
    out=[]
    for mode in MODES:
        a=read(OUT/f"metrics_shared_{mode}.tsv")
        c=read(OUT/f"contrasts_shared_{mode}.tsv")
        threshold_backgrounds=[]
        for bg in BACKS:
            r={q:[int(m["recovered"]) for m in a if m["background"]==bg and m["cell"]==q]
               for q in ALL}
            if min(r["HH"])<90 and min(r["LH"])<90 and min(r["HL"])<90 and min(r["LL"])>=90:
                threshold_backgrounds.append(bg)
        reversed_groups=[]
        for j in ("1","2"):
            q=[next(r for r in c if r["case_id"]==j and r["background"]==b and
                    r["contrast"]=="O1_minus_O2") for b in BACKS]
            for metric,eps in (("signed_delta",.000001),("rate_delta",.000001)):
                first,second=(float(r[metric]) for r in q)
                if first*second<0 and min(abs(first),abs(second))>eps:
                    reversed_groups.append(f"truth {j} {metric}")
        picked=set(threshold_backgrounds)
        if reversed_groups:
            picked.update(BACKS)
        for bg in BACKS:
            out.append(dict(mode=mode,background=bg,only_both_pass=bg in threshold_backgrounds,
                            opposite_ordering="; ".join(reversed_groups),run_split=bg in picked))
    write(OUT/"independent_fit_gate.tsv",out)
    for r in out:
        print(r["mode"],r["background"],r["run_split"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=("prepare","run","collect","repair","summarize","decide"))
    p.add_argument("variant",nargs="?",choices=VARIANTS)
    p.add_argument("mode",nargs="?",choices=MODES)
    p.add_argument("phase",nargs="?",choices=("pilot","expanded"))
    p.add_argument("--background",nargs="+",choices=BACKS,default=BACKS)
    p.add_argument("--shard",type=int)
    p.add_argument("--width",type=int,default=8)
    p.add_argument("--cores",type=int,default=4)
    a=p.parse_args()
    if a.action=="decide":
        decide()
    elif a.action=="prepare":
        prepare(a.variant,a.mode,a.phase,a.background)
    elif a.action=="run":
        run(a.variant,a.mode,a.phase,a.shard,a.width,a.cores)
    elif a.action=="summarize":
        accuracy(a.variant,a.mode)
    else:
        verify(a.variant,a.mode,a.phase,a.action=="repair")
