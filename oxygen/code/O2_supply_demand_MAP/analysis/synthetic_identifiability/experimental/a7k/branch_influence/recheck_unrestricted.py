"""Refit archived unconstrained estimates reached from better nested solutions."""
from collections import defaultdict
from concurrent.futures import as_completed

import numpy as np

from recovery_diagnostics import (OUT, archived, case, pool, read, results, row_key, sample, score,
                                  task, work, write)


def run():
    candidates=defaultdict(list)
    for phase in ("oracle","shape"):
        for r in read(results(phase)):
            if r["validated"]=="True":
                candidates[row_key(r)[:4]].append(r)
    jobs=[]
    history={}
    for (mode,j,i,s), rows in sorted(candidates.items()):
        previous=archived(mode)[j,i,s,"relaxed"]
        r=min(rows,key=lambda row:float(row["loss"]))
        if float(previous["loss"])-float(r["loss"])<=.025:
            continue
        k,c,_=case(mode)
        u=c.anchor(r)
        alt=np.clip(u+np.random.default_rng(221+i).normal(0.,.00024,len(u)),0.,1.)
        query=task(mode,j,i,s,"free")
        y=sample(mode,j,i,s)
        jobs.append((query,y,u,alt,np.full(len(u),np.nan)))
        history[row_key(query)]=(previous,r,y)
    output=[]
    if not jobs:
        return
    with pool(min(6,len(jobs))) as ex:
        for f in as_completed(ex.submit(work,j) for j in jobs):
            r=f.result()
            mode,j,*_=row_key(r)
            prev,nested,y=history[row_key(r)]
            assert r["validated"] and r["loss"]<=float(nested["loss"])+.005
            k,c,_=case(mode)
            p=c.unmap(c.anchor(r))[0]
            assert abs(k.parts(p,y)[1]-r["loss"])<1e-8
            first,x=score(mode,j,prev),score(mode,j,r)
            output.append(dict(r,estimator="relaxed",original_loss=prev["loss"],
                               old_error=first,new_error=x,
                               flipped=(abs(first)<=.10)!=(abs(x)<=.10),
                               nested_loss=nested["loss"]))
            print(row_key(r),r["loss"],prev["loss"],first,x,output[-1]["flipped"],flush=True)
    write(OUT/"unrestricted_rechecks.tsv",output)


if __name__=="__main__":
    run()
