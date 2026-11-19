"""Frozen, paired 2 x 2 exchanges of A7K karyotypes across the two cohorts."""
import argparse
from concurrent.futures import as_completed
import hashlib
from pathlib import Path
import sys
import zlib

import numpy as np
from scipy.stats import binomtest

OUT = Path(__file__).resolve().parent
PREVIOUS = OUT.with_name("20261119_A7K_branch_influence_01")
sys.path.insert(0, str(PREVIOUS))
from recovery_diagnostics import (MODES, archived, case, pool, read,
                                  row_key, sample, score, task, work, write)

HYBRIDS = ("2N", "4N")
LAYOUT = {"HH": "both high", "LH": "2N low, 4N high",
          "HL": "2N high, 4N low", "LL": "both low"}


def directory(mode, phase):
    return OUT / f"batch_{mode}_{phase}"


def result_path(mode, phase):
    return OUT / f"fits_{mode}_{phase}.tsv"


def digest(y):
    return hashlib.sha256(np.ascontiguousarray(y).tobytes()).hexdigest()


def exchange(mode, j, i, group):
    k, c, _ = case(mode)
    assert k.n == 221 and k.kary.shape == (6, 133)
    assert tuple(c.kary_labels) == ("2N_O1","2N_O2","4N_O1","4N_O2")
    high, low = (sample(mode, j, i, s) for s in ("up", "down"))
    boundary = k.n + 2 * 133
    assert np.array_equal(high[:boundary], low[:boundary])
    h = high[boundary:].reshape(4, 133)
    l = low[boundary:].reshape(4, 133)
    assert np.array_equal(h.sum(axis=1), k.kary[2:].sum(axis=1))
    assert np.array_equal(l.sum(axis=1), h.sum(axis=1))
    start_index = 0 if group == "2N" else 2
    ids = np.arange(start_index, start_index + 2)
    other = np.setdiff1d(np.arange(4), ids)
    a = h.copy()
    a[ids] = l[ids]
    reverse = l.copy()
    reverse[other] = h[other]
    y = np.r_[high[:boundary], a.ravel()]
    rev = np.r_[low[:boundary], reverse.ravel()]
    assert np.array_equal(y, rev)
    assert np.array_equal(a[ids], l[ids]) and np.array_equal(a[other], h[other])
    old = archived(mode)
    x, z = (c.anchor(old[j, i, s, "relaxed"]) for s in ("up", "down"))
    # Both fitted mechanisms seed the optimizer with the corresponding cohort tilts.
    inherited = [pos for pos, name in enumerate(c.names) if name.startswith("tilt_")]
    for pos in inherited:
        value = z[pos] if c.names[pos].endswith(group) else x[pos]
        x[pos] = z[pos] = value
    return y, x, z, dict(
        target_sha256=digest(y),reverse_target_sha256=digest(rev),
        common_sha256=digest(y[:boundary]),distinct_karyotype_bins=int(np.sum(a!=h)),
        n_2N_O1=int(a[0].sum()),n_2N_O2=int(a[1].sum()),
        n_4N_O1=int(a[2].sum()),n_4N_O2=int(a[3].sum()))


def load(mode,phase):
    path=result_path(mode,phase)
    return read(path) if path.exists() else []


def stability(rows, group, j):
    subset=[r for r in rows if r["branch"]==group and int(r["case_id"])==j]
    assert len(subset)==20 and sorted(int(r["replicate"]) for r in subset)==list(range(1,21))
    return sum(r["validated"]=="True" for r in subset)


def prepare(mode, phase):
    if phase=="expanded":
        pilot=load(mode,"pilot")
        assert len(pilot)==80
        assert all(stability(pilot,g,j)>=18 for g in HYBRIDS for j in (1,2))
        lo,hi=21,101
    else:
        lo,hi=1,21
    meta=[]
    targets=[]
    starts=[]
    alternatives=[]
    for i in range(lo,hi):
        for j in (1,2):
            for group in HYBRIDS:
                y,x,z,audit=exchange(mode,j,i,group)
                row=task(mode,j,i,"up","cohort_swap",group)
                row.update(audit)
                meta.append(row)
                targets.append(y)
                starts.append(x)
                alternatives.append(z)
    path=directory(mode,phase)
    path.mkdir(parents=True,exist_ok=True)
    write(path/"manifest.tsv",meta)
    np.savez_compressed(path/"frozen.npz",target=targets,start=starts,alternative=alternatives)
    print(len(meta),flush=True)


def run(mode,phase,index,width,cores):
    path=directory(mode,phase)
    metadata=read(path/"manifest.tsv")
    frozen=np.load(path/"frozen.npz")
    start=index*width
    stop=min(start+width,len(metadata))
    arrays=[frozen[k][start:stop] for k in ("target","start","alternative")]
    rows=[]
    items=[(m,y,x,z,np.full(len(x),np.nan)) for m,y,x,z in zip(metadata[start:stop],*arrays)]
    with pool(min(cores,len(items))) as ex:
        for f in as_completed(ex.submit(work,item) for item in items):
            ans=f.result()
            rows.append(ans)
            write(path/f"part_{index:04d}.tsv",rows)
            print(row_key(ans),ans["validated"],ans["attempts"],flush=True)


def collect(mode,phase):
    path=directory(mode,phase)
    rows=[r for f in sorted(path.glob("part_*.tsv")) for r in read(f)]
    manifest=read(path/"manifest.tsv")
    assert len(rows)==len(manifest) and {row_key(r) for r in rows}=={row_key(r) for r in manifest}
    frozen=np.load(path/"frozen.npz")["target"]
    lookup={row_key(r):y for r,y in zip(manifest,frozen)}
    k,c,_=case(mode)
    expected_seed=np.array([float(r["fixed_recorded_seed"]) for r in k.a7])
    report=[]
    for r in sorted(rows,key=row_key):
        _,j,i,_,_,g,*_=row_key(r)
        y=lookup[row_key(r)]
        assert digest(y)==r["target_sha256"]==r["reverse_target_sha256"]
        p,_=c.unmap(c.anchor(r))
        ar,f,*_=k.parts(p,y)
        error=abs(f-float(r["loss"]))
        assert error<1e-7
        tv,_,_,seed=k.composition(p,True)
        seed_error=float(max(abs(seed-expected_seed)))
        assert seed_error<1e-7
        assert r["validated"]!="True" or (ar[1]>=-1e-10 and max(tv)<=k.cap+1e-10)
        err=score(mode,j,r)
        report.append(dict(r,normalized_error=err,dev_replay_error=error,
                           max_seed_error=seed_error,
                           recovered=r["validated"]=="True" and abs(err)<=.1))
    write(result_path(mode,phase),report)
    print(mode,phase,len(report),"replayed",sum(r["validated"]=="True" for r in report),"validated",
          "max error",max(r["dev_replay_error"] for r in report),flush=True)


def repair(mode,phase):
    rows=load(mode,phase)
    redo={row_key(r):r for r in rows if r["validated"]!="True"}
    if not redo:
        return
    path=directory(mode,phase)
    y=np.load(path/"frozen.npz")["target"]
    k,c,_=case(mode)
    expected_seed=np.array([float(r["fixed_recorded_seed"]) for r in k.a7])
    jobs=[]
    targets={}
    for pos,m in enumerate(read(path/"manifest.tsv")):
        key=row_key(m)
        if key not in redo:
            continue
        initial=c.anchor(redo[key])
        alternative=np.clip(initial+np.random.default_rng(pos+519).normal(0,.00015,len(initial)),0,1)
        jobs.append((m,y[pos],initial,alternative,np.full(len(initial),np.nan)))
        targets[key]=y[pos]
    repaired={}
    with pool(min(8,len(jobs))) as ex:
        for f in as_completed(ex.submit(work,q) for q in jobs):
            r=f.result()
            key=row_key(r)
            if r["validated"] and r["loss"]<=float(redo[key]["loss"])+.02:
                p,_=c.unmap(c.anchor(r))
                ar,obj,*_=k.parts(p,targets[key])
                tv,_,_,seed=k.composition(p,True)
                seed_error=float(max(abs(seed-expected_seed)))
                replay=abs(obj-r["loss"])
                assert replay<1e-7 and ar[1]>=-1e-10 and max(tv)<=k.cap+1e-10 and seed_error<1e-7
                err=score(mode,key[1],r)
                repaired[key]=dict(r,normalized_error=err,dev_replay_error=replay,
                                   max_seed_error=seed_error,recovered=abs(err)<=.10)
            print(key,"revalidated" if key in repaired else "unresolved",flush=True)
    if repaired:
        write(result_path(mode,phase),[repaired.get(row_key(r),r) for r in rows])
        write(OUT/f"repaired_{mode}_{phase}.tsv",list(repaired.values()))


def interval(v,indices):
    return np.quantile(v[indices].mean(axis=1),(.025,.975))


def summarize(mode):
    pilot=load(mode,"pilot")
    data=pilot+load(mode,"expanded")
    assert len(pilot)==80 and len(data)==400 and len({row_key(r) for r in data})==400
    lookup={(r["branch"],int(r["case_id"]),int(r["replicate"])):r for r in data}
    recheck_path=OUT/f"crossfit_rechecks_{mode}.tsv"
    patches={}
    if recheck_path.exists():
        patches={(r["configuration"],int(r["case_id"]),int(r["replicate"])):r
                 for r in read(recheck_path)}
    all_runs=[]
    stats=[]
    contrasts=[]
    for j in (1,2):
        errors={}
        valid={}
        old=archived(mode)
        for cell in LAYOUT:
            if cell in ("HH","LL"):
                series=[old[j,i,"up" if cell=="HH" else "down","relaxed"] for i in range(1,101)]
            else:
                series=[lookup["2N" if cell=="LH" else "4N",j,i] for i in range(1,101)]
            series=[patches.get((cell,j,i),r) for i,r in enumerate(series,1)]
            errors[cell]=np.array([score(mode,j,r) for r in series])
            valid[cell]=np.array([r["validated"]=="True" for r in series])
            for i,r in enumerate(series,1):
                all_runs.append(dict(mode=mode,case_id=j,replicate=i,configuration=cell,
                    value=float(r["p_misseg"]),normalized_error=errors[cell][i-1],
                    validated=valid[cell][i-1],recovered=valid[cell][i-1] and abs(errors[cell][i-1])<=.10))
        passing={q:(np.abs(e)<=.10)&valid[q] for q,e in errors.items()}
        rng=np.random.default_rng(zlib.crc32(f"{mode}:{j}:paired cohort".encode()))
        indices=rng.integers(100,size=(8000,100))
        for q,e in errors.items():
            v=valid[q]
            recovered=int(passing[q].sum())
            ci=binomtest(recovered,100).proportion_ci()
            stats.append(dict(mode=mode,case_id=j,configuration=q,n=100,validated=int(v.sum()),
                 recovered=recovered,rate=recovered/100,ci_lo=ci.low,ci_hi=ci.high,
                 signed_mean=float(e[v].mean()),sd=float(e[v].std(ddof=1)),
                 low_fail=int(((e<-.10)&v).sum()),high_fail=int(((e>.10)&v).sum())))
        design=(("2N_given_4N_high","HH","LH"),("2N_given_4N_low","HL","LL"),
                ("4N_given_2N_high","HH","HL"),("4N_given_2N_low","LH","LL"),
                ("both","HH","LL"))
        for label,a,b in design:
            v=valid[a]&valid[b]
            diff=errors[b]-errors[a]
            rd=passing[b].astype(float)-passing[a]
            lo,hi=interval(rd,indices)
            loe,hie=interval(diff,indices) if v.all() else (np.nan,np.nan)
            contrasts.append(dict(mode=mode,case_id=j,contrast=label,from_cell=a,to_cell=b,n=100,
                 paired_validated=int(v.sum()),signed_delta=float(diff[v].mean()),
                 signed_lo=loe,signed_hi=hie,rate_delta=float(rd.mean()),rate_lo=lo,rate_hi=hi,
                 improved=int((rd==1).sum()),reversed=int((rd==-1).sum())))
        e=errors["LL"]-errors["LH"]-errors["HL"]+errors["HH"]
        q={name:v.astype(float) for name,v in passing.items()}
        r=q["LL"]-q["LH"]-q["HL"]+q["HH"]
        v=np.logical_and.reduce(list(valid.values()))
        lo,hi=interval(r,indices)
        loe,hie=interval(e,indices) if v.all() else (np.nan,np.nan)
        contrasts.append(dict(mode=mode,case_id=j,contrast="interaction",from_cell="",to_cell="",n=100,
               paired_validated=int(v.sum()),signed_delta=float(e[v].mean()),signed_lo=loe,signed_hi=hie,
               rate_delta=float(r.mean()),rate_lo=lo,rate_hi=hi,improved="",reversed=""))
    write(OUT/f"paired_runs_{mode}.tsv",all_runs)
    write(OUT/f"factorial_metrics_{mode}.tsv",stats)
    write(OUT/f"factorial_contrasts_{mode}.tsv",contrasts)
    for q in LAYOUT:
        vals=[r for r in stats if r["configuration"]==q]
        a,b=[r["recovered"] for r in vals]
        print(mode,q,a,b,"passes",min(a,b)>=90,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=("prepare","run","collect","repair","summarize"))
    p.add_argument("mode",choices=MODES)
    p.add_argument("phase",nargs="?",choices=("pilot","expanded"))
    p.add_argument("--shard",type=int)
    p.add_argument("--width",type=int,default=8)
    p.add_argument("--cores",type=int,default=4)
    a=p.parse_args()
    if a.action=="run":
        run(a.mode,a.phase,a.shard,a.width,a.cores)
    elif a.action=="summarize":
        summarize(a.mode)
    else:
        globals()[a.action](a.mode,a.phase)
