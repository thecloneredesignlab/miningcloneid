"""Complete summaries from finished numerical profile receipts."""
from concurrent.futures import as_completed
import sys

import numpy as np

from branch_influence import MODES, NAMES, OUT, THRESHOLD, masked, pool, read, solution, stable, tag, take, write


def recover(mode,mask):
    prefix=tag(mode,mask)
    k,c,y,_=masked(mode,mask)
    reference=OUT/f"profile_ref_{prefix}.tsv"
    best=read(reference if reference.exists() else reference.with_suffix(".tsv.tmp"))[0]
    raw=read(OUT/f"profile_points_{prefix}.tsv")
    winner=float(best["loss"])
    assert max(abs(float(r["loss"])-float(r["lr"])-winner) for r in raw)<1e-8
    a=float(best["x_buffer_beta"])
    at10=next(r for r in raw if float(r["coordinate"])==1.)
    points={x:take([r for r in raw if float(r["coordinate"])==x])
            for x in {float(r["coordinate"]) for r in raw}}
    assert [x for x in (0.,.08,.15,.23,.33,.46,.60,.72,.78538713,.88,.96,1.)
            if abs(x-a)>1e-7 and x not in points]==[]
    p=c.unmap(c.anchor(best))[0]
    assert abs(k.parts(p,y)[1]-winner)<1e-8
    added=[]
    bounds=[]
    labels=[]
    for direction in (-1,1):
        box=0. if direction<0 else 1.
        for _ in range(10):
            if abs(a-box)<1e-8:
                bounds.append(box)
                labels.append("box")
                break
            inside=[v for v,r in points.items() if direction*(v-a)>1e-8 and
                    float(r["all_constraint_margin"])>=0. and float(r["loss"])-winner<=THRESHOLD]
            inner=(min if direction<0 else max)([a,*inside])
            if inner==box:
                bounds.append(box)
                labels.append("box")
                break
            outer=min((v for v in points if direction*(v-inner)>0.),key=lambda v:abs(v-inner))
            ri=best if inner==a else points[inner]
            ro=points[outer]
            li,lo=(float(row["loss"])-winner for row in (ri,ro))
            if abs(outer-inner)<.003 and stable(ri) and stable(ro):
                bounds.append(inner + (outer-inner)*(THRESHOLD-li)/(lo-li))
                labels.append("crossing")
                break
            midpoint=(inner+outer)*.5
            move=solution(mode,mask,c.anchor(ri),"receipt_refine",midpoint,True)
            move["lr"]=float(move["loss"])-winner
            points[midpoint]=move
            added.append(move)
            for v in (inner,outer):
                if v==a or stable(points[v]):
                    continue
                move=solution(mode,mask,c.anchor(points[v]),"receipt_stationarity",v,True)
                move["lr"]=float(move["loss"])-winner
                points[v]=take([move,points[v]])
                added.append(move)
        else:
            assert abs(outer-inner)<.01 and lo>li
            bounds.append(inner + (outer-inner)*(THRESHOLD-li)/(lo-li))
            labels.append("crossing_recheck")
    to_beta=lambda x:c.unmap(np.where(np.arange(len(c.names))==c.names.index("buffer_beta"),x,c.anchor(best)))[0][c.active[c.names.index("buffer_beta")]]
    lr=float(at10["loss"])-winner
    r=dict(mode=mode,mask=mask,branches=",".join(n for j,n in enumerate(NAMES) if mask&1<<j),
           n_kary_cells=int(k.kary[2:].sum()),beta=float(best["buffer_beta"]),
           lower=to_beta(bounds[0]),upper=to_beta(bounds[1]),lower_type=labels[0],upper_type=labels[1],
           x_lower=bounds[0],x_upper=bounds[1],deviance=winner,beta10_deviance=float(at10["loss"]),
           beta10_lr=lr,beta10_supported=lr<=THRESHOLD,free_stationarity=float(best["stationarity"]),
           beta10_stationarity=float(at10["stationarity"]))
    write(OUT/f"profile_points_{prefix}.tsv",raw+added)
    write(OUT/f"profile_ref_{prefix}.tsv",[best,{key:at10[key] for key in best}])
    write(OUT/f"interval_{prefix}.tsv",[r])
    return mode,mask,lr,bounds,labels,len(added)


if __name__=="__main__":
    jobs=[(m,k) for m in MODES for k in range(1,15) if
          not (OUT/f"interval_{tag(m,k)}.tsv").exists()]
    with pool(8) as executor:
        for f in as_completed(executor.submit(recover,m,k) for m,k in jobs):
            print(f.result(),flush=True)
