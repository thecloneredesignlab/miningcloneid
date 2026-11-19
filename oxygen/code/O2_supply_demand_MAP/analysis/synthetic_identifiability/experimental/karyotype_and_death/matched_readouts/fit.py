"""Multi-start fits for the shared graft and three observation sets."""
import argparse
from concurrent.futures import as_completed

import numpy as np

from modal import (BASE,HERE,STAGE3,TARGETS,MODES,case,fixed_for,pool,read,receipt,
                   solve,stable,write)
PRIMARY=("baseline","death","a7k")


def take(rows):
    good=[r for r in rows if float(r["donor_margin"])>=0.]
    if not good:
        return min(rows,key=lambda r:float(r["loss"]))
    low=min(good,key=lambda r:float(r["loss"]))
    known=[r for r in good if stable(r) and float(r["loss"])<=float(low["loss"])+.015]
    return min(known or good,key=lambda r:float(r["loss"]))


def run_one(job):
    m,label,anchor,lock,limit,divisor=job
    k,c,y=case(m)
    f=solve(k,c,y,anchor,limit=limit,fixed=fixed_for(m,c)|lock,divisor=divisor)
    return receipt(k,c,y,f,mode=m,origin=label)


def starts(mode):
    k,c,y=case(mode)
    r=read(STAGE3/"best_autosomes_joint.tsv")[0]
    points=[c.anchor(r)]
    sources=[r for r in read(STAGE3/"profile_points_autosomes_joint.tsv")
             if float(r["donor_margin"])>=0 and float(r["lr"])<=5.99147]
    sources += read(BASE/"truths.tsv")
    sources += read(STAGE3/"best_autosomes_count.tsv")
    for m in MODES:
        p=HERE/f"best_{m}.tsv"
        if p.exists():
            sources.extend(read(p))
    rank=[]
    for r in sources:
        z=c.anchor(r)
        ar,f,g,nk,nd=k.parts(c.unmap(z)[0],y)
        rank.append((ar[1]<0.,f+160.*max(-ar[1],0.),z))
    for _,_,z in sorted(rank,key=lambda r:r[:2]):
        if min(np.linalg.norm(z-a) for a in points)>.12:
            points.append(z)
        if len(points)>=5:
            break
    seed=c.anchor(read(BASE/"truths.tsv")[0])
    for value in (.0,.10,.25):
        z=seed.copy()
        z[c.names.index("alpha_o2")]=value
        z[c.names.index("lam_max")]=.97
        z[c.names.index("mu_hp")]=.035
        z[c.names.index("p_misseg")]=.04
        z[c.names.index("buffer_beta")]=.44
        z[c.names.index("gamma_growth")]=.05
        if k.use_death:
            z[-4:]=[.28,.50,.75,.225]
        points.append(z)
    if k.use_death:
        for values in ((.22,.35,.15,.24),(.5,.69,.98,.16)):
            z=points[1].copy()
            z[-4:]=values
            points.append(z)
    return points


def run(modes,workers):
    jobs=[(m,f"start_{i+1}",z,{},100,38.) for m in modes for i,z in enumerate(starts(m))]
    rows={m:[] for m in modes}
    with pool(workers) as ex:
        for f in as_completed(ex.submit(run_one,j) for j in jobs):
            a=f.result()
            rows[a["mode"]].append(a)
            print(a["mode"],a["origin"],a["loss"],a["donor_margin"],
                  a["stationarity"],flush=True)
    for mode,data in rows.items():
        data.sort(key=lambda r:(not stable(r),float(r["loss"])))
        write(HERE/f"fits_{mode}.tsv",data)
        best=take(data)
        if not stable(best):
            k,c,_=case(mode)
            for j,div in enumerate((12.,72.,28.)):
                t=run_one((mode,f"polish_{j}",c.anchor(best),{},145,div))
                data.append(t)
                best=take([best,t])
                if stable(best):
                    break
        write(HERE/f"fits_{mode}.tsv",data)
        write(HERE/f"best_{mode}.tsv",[best])
        assert stable(best),(mode,best)
        print("chosen",mode,best["loss"],best["donor_margin"],best["stationarity"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("modes",nargs="+",choices=MODES)
    p.add_argument("--workers",type=int,default=9)
    a=p.parse_args()
    run(a.modes,a.workers)
