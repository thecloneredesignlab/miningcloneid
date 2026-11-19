"""Recheck only unvalidated archived fits against their original frozen data."""
import argparse
from concurrent.futures import as_completed

import numpy as np

from recovery_diagnostics import batch, case, pool, read, results, row_key, score, work, write


def recheck(phase):
    folder=batch(phase)
    data=read(results(phase))
    failures={row_key(r):r for r in data if r["validated"]!="True"}
    if not failures:
        return
    frozen=np.load(folder/"frozen.npz")
    y,z=(frozen[k] for k in ("target","fixed"))
    jobs=[]
    for i,r in enumerate(read(folder/"manifest.tsv")):
        key=row_key(r)
        if key not in failures:
            continue
        c=case(r["mode"])[1]
        x=c.anchor(failures[key])
        new=np.clip(x+np.random.default_rng(i).normal(0,.00013,len(x)),0,1)
        jobs.append((r,y[i],x,new,z[i]))
    better={}
    with pool(min(8,len(jobs))) as ex:
        for f in as_completed(ex.submit(work,item) for item in jobs):
            r=f.result()
            key=row_key(r)
            orig=failures[key]
            if r["validated"] and float(r["loss"])<=float(orig["loss"])+.02:
                ar,loss,*_=case(r["mode"])[0].parts(case(r["mode"])[1].unmap(case(r["mode"])[1].anchor(r))[0],
                        next(item[1] for item in jobs if row_key(item[0])==key))
                assert abs(loss-r["loss"])<1e-7 and ar[1]>=-1e-10
                err=score(r["mode"],int(r["case_id"]),r)
                better[key]=dict(r,delta_x_p_misseg=err,recovered=abs(err)<=.10,
                                 deviance_replay_error=abs(loss-r["loss"]))
            print(key,"revalidated" if key in better else "unresolved",flush=True)
    if better:
        write(results(phase),[better.get(row_key(r),r) for r in data])
        write(results(phase).with_name(f"rechecked_{phase}.tsv"),list(better.values()))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase")
    recheck(p.parse_args().phase)
