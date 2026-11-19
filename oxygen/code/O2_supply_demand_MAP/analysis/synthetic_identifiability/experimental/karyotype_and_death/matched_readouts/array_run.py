"""Freeze local simulation draws, fit independent scheduler slices, and collect them."""
import argparse
from concurrent.futures import as_completed
from pathlib import Path

import numpy as np

from modal import HERE,case,pool,read,stable,write
from recover_modal import draw,worker

PRIMARY=("stock_death","stock_a7k")


def paths(mode,phase):
    return HERE/f"arrays_{mode}_{phase}"


def freeze(mode,phase):
    k,c,_=case(mode)
    data=paths(mode,phase)
    data.mkdir(exist_ok=True)
    source=HERE/f"recovery_{mode}.tsv"
    merged=HERE/f"accepted_pilot_{mode}.tsv"
    if phase=="expanded" and merged.exists():
        source=merged
    current=read(source) if source.exists() else []
    previous={(int(r["case_id"]),int(r["replicate"])):r for r in current}
    if phase=="expanded":
        for j in (1,2):
            completed=[previous.get((j,i),{}) for i in range(1,21)]
            assert all(completed) and sum(r["validated"]=="True" for r in completed)>=18
        cases=[(j,i) for i in range(21,101) for j in (1,2)]
    else:
        cases=[(j,i) for i in range(20,0,-1) for j in (1,2)
               if previous.get((j,i),{}).get("validated")!="True"]
    starts=read(HERE/f"fits_{mode}.tsv")+read(HERE/f"profile_points_{mode}.tsv")+current
    starter=[]
    ordered=([r for r in current if r["validated"]=="True"] if phase=="expanded" else [])
    ordered+=sorted(starts,key=lambda r:float(r["loss"]))
    for r in ordered:
        if not stable(r):
            continue
        z=c.anchor(r)
        if all(np.linalg.norm(z-a)>.085 for a in starter):
            starter.append(z)
        if len(starter)>=24:
            break
    assert starter and cases
    matrix=[]
    for j,i in cases:
        y,_=draw(mode,j,i)
        matrix.append(y)
    np.savez_compressed(data/"frozen.npz",targets=matrix,indices=cases,starts=starter)
    print(len(cases),flush=True)


def run(mode,phase,shard,width,cores):
    k,c,_=case(mode)
    data=paths(mode,phase)
    info=np.load(data/"frozen.npz")
    ids=info["indices"]
    y=info["targets"]
    starts=info["starts"]
    low=shard*width
    high=min(len(ids),low+width)
    tasks=[(mode,int(ids[i,0]),int(ids[i,1]),starts,None,y[i]) for i in range(low,high)]
    output=[]
    if not tasks:
        return
    with pool(min(cores,len(tasks))) as ex:
        pending={ex.submit(worker,t):t for t in tasks}
        for f in as_completed(pending):
            r=f.result()
            output.append(r)
            write(data/f"part_{shard:03d}.tsv",output)
            print(r["case_id"],r["replicate"],r["validated"],r["attempts"],flush=True)


def collect(mode,phase):
    data=paths(mode,phase)
    assigned=np.load(data/"frozen.npz")["indices"]
    all_rows=[r for path in sorted(data.glob("part_*.tsv")) for r in read(path)]
    keys={(int(r["case_id"]),int(r["replicate"])) for r in all_rows}
    assert keys=={tuple(i) for i in assigned} and len(keys)==len(all_rows), (len(keys),len(assigned))
    k,c,_=case(mode)
    file=HERE/f"recovery_{mode}.tsv"
    merged=HERE/f"accepted_pilot_{mode}.tsv"
    prior=read(file) if file.exists() else []
    if phase=="expanded" and merged.exists():
        prior+=read(merged)
    locator={(int(r["case_id"]),int(r["replicate"])):r for r in prior}
    audit=[]
    for r in all_rows:
        j,i=int(r["case_id"]),int(r["replicate"])
        y,truth=draw(mode,j,i)
        p,_=c.unmap(c.anchor(r))
        ar,loss,*_=k.parts(p,y)
        discrepancy=abs(loss-float(r["loss"]))
        assert discrepancy<1e-7,(j,i,discrepancy)
        assert (ar[1]>=-3e-12 if r["validated"]=="True" else True)
        key=(j,i)
        if key in locator:
            a=locator[key]
            if a["validated"]=="True" and (r["validated"]!="True" or float(a["loss"])<float(r["loss"])):
                continue
        locator[key]=r
        audit.append(dict(mode=mode,phase=phase,case_id=j,dataset=i,
                          deviance_difference=discrepancy,minimum_margin=ar[1],
                          validated=r["validated"]))
    write(merged if phase=="pilot" else file,[locator[key] for key in sorted(locator)])
    if audit:
        write(HERE/f"array_replay_{mode}_{phase}.tsv",audit)
    print(mode,phase,len(all_rows),"locally replayed",len(locator),"distinct datasets",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("action",choices=["freeze","run","collect"])
    p.add_argument("mode",choices=PRIMARY)
    p.add_argument("phase",choices=["pilot","expanded"])
    p.add_argument("--shard",type=int)
    p.add_argument("--width",type=int,default=4)
    p.add_argument("--cores",type=int,default=2)
    a=p.parse_args()
    if a.action=="run":
        run(a.mode,a.phase,a.shard,a.width,a.cores)
    else:
        (freeze if a.action=="freeze" else collect)(a.mode,a.phase)
