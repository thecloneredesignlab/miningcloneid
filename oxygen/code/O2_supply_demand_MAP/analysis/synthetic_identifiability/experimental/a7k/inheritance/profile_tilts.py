"""Fit bounded A7K tilt models and profile their mechanism parameters."""
import argparse
from concurrent.futures import as_completed
from functools import lru_cache
import math

import numpy as np

from tilt_model import (HERE,PREV,BASE,TARGETS,GROUP,BRANCH,THRESHOLD,case,cap,
                         pool,read,run_fit,stable,take,write)


def convert(mode,r):
    k,c,_=case(mode)
    x=c.anchor(r)
    source=str(r.get("mode",""))
    if k.cap and source.startswith(("group_","branch_")):
        old_domain=12.+400.*cap(source)
        for n,name in enumerate(c.names[c.dim:],c.dim):
            prev=r.get(f"x_{name}")
            if prev is None and source.startswith("group_"):
                prev=r.get(f"x_{name.rsplit('_',1)[0]}")
            if prev is not None:
                x[n]=np.clip(.5+(float(prev)-.5)*old_domain/c.domain,0.,1.)
    return x


def seeds(mode):
    k,c,y=case(mode)
    source=[read(PREV/"best_stock_a7k.tsv")[0],
            read(PREV/"best_stock_baseline.tsv")[0]]
    source+=read(PREV/"fits_stock_a7k.tsv")
    source+=read(PREV/"profile_points_stock_a7k.tsv")
    source+=read(BASE/"truths.tsv")
    if mode.startswith("branch"):
        group=HERE/f"best_group_{mode.split('_')[1]}.tsv"
        if group.exists():
            source.insert(0,read(group)[0])
    for other in GROUP+BRANCH:
        path=HERE/f"best_{other}.tsv"
        if other!=mode and path.exists():
            source.extend(read(path))
    x0=convert(mode,source[0])
    ranked=[]
    for r in source:
        z=convert(mode,r)
        ar,f,*_=k.parts(c.unmap(z)[0],y)
        ranked.append((ar[1]<0.,f+300.*max(-ar[1],0.),z))
    out=[x0,convert(mode,source[1])]
    for _,_,z in sorted(ranked,key=lambda a:a[:2]):
        if min(np.linalg.norm(z-u) for u in out)>.13:
            out.append(z)
        if len(out)>=6:
            break
    plus=x0.copy()
    plus[c.dim:]=np.resize([.505,.50,.495,.50],len(plus)-c.dim)
    minus=x0.copy()
    minus[c.dim:]=1.-plus[c.dim:]
    return out+[plus,minus]


def fits(modes,workers):
    tasks=[(mode,f"start_{i}",x,{},105,40.) for mode in modes
          for i,x in enumerate(seeds(mode))]
    result={m:[] for m in modes}
    with pool(workers) as ex:
        for f in as_completed(ex.submit(run_fit,t) for t in tasks):
            r=f.result()
            result[r["mode"]].append(r)
            print("fit",r["mode"],r["origin"],r["loss"],r["all_constraint_margin"],
                   r["stationarity"],flush=True)
    for m,rows in result.items():
        k,c,_=case(m)
        best=take(rows)
        if not stable(best):
            for j,scale in enumerate((17.,52.,36.)):
                new=run_fit((m,f"polish_{j}",c.anchor(best),{},150,scale))
                rows.append(new)
                best=take([best,new])
                if stable(best):
                    break
        write(HERE/f"fits_{m}.tsv",rows)
        write(HERE/f"best_{m}.tsv",[best])
        print("chosen",m,best["loss"],best["p_misseg"],best["buffer_beta"],best["stationarity"],
              flush=True)
        assert stable(best),(m,best)


@lru_cache
def ref(mode):
    k,c,y=case(mode)
    r=read(HERE/f"best_{mode}.tsv")[0]
    return r,c.anchor(r)


@lru_cache
def profile_seeds(mode):
    archive=read(PREV/"profile_points_stock_a7k.tsv")
    data=read(HERE/f"fits_{mode}.tsv")+archive        +read(PREV/"profile_points_stock_baseline.tsv")+read(BASE/"truths.tsv")
    return np.array([convert(mode,r) for r in data if float(r["donor_margin"])>=0.])


def profile_at(mode,name,v,warm,try_alt=True):
    k,c,y=case(mode)
    r,x=ref(mode)
    j=c.names.index(name)
    anchors=[("carry",warm)]
    if try_alt:
        a=profile_seeds(mode)
        opts=[]
        for i in np.argsort(abs(a[:,j]-v))[:14]:
            z=a[i].copy()
            z[j]=v
            if np.linalg.norm(np.delete(z-warm,j))<.025:
                continue
            ar,f,*_=k.parts(c.unmap(z)[0],y)
            opts.append((f+250.*max(0.,-ar[1]),z))
        if opts:
            anchors.append(("alternative",min(opts,key=lambda q:q[0])[1]))
    bank=[run_fit((mode,label,u,{j:v},64,38.)) for label,u in anchors]
    best=take(bank)
    if not stable(best) and (float(best["loss"])-float(r["loss"])>THRESHOLD or v in (0.,1.)):
        bank.append(run_fit((mode,"polish",c.anchor(best),{j:v},105,16.)))
    for point in bank:
        point.update(parameter=name,coordinate=v,lr=float(point["loss"])-float(r["loss"]))
    return take(bank),bank


def side(task):
    mode,name,sign=task
    k,c,y=case(mode)
    r,x=ref(mode)
    j=c.names.index(name)
    center=x[j]
    strict=next(t for t in read(HERE/"intervals_strict.tsv") if t["parameter"]==name)
    grid={0.,.05,.15,.30,.60,1.,float(np.clip(center+sign*.10001,0.,1.)),
          float(strict["x_lower"]),float(strict["x_upper"])}
    if name=="buffer_beta":
        grid.add(float(next(t["x_lower"] for t in read(PREV/"intervals_stock_baseline.tsv")
                            if t["parameter"]==name)))
    prior=HERE/f"profile_points_{mode}.tsv"
    points={}
    output=[]
    if prior.exists():
        history=[r for r in read(prior) if r["parameter"]==name and
                   sign*(float(r["coordinate"])-center)>1e-8]
        for v in {float(r["coordinate"]) for r in history}:
            rows=[dict(r,lr=float(r["loss"])-float(ref(mode)[0]["loss"]))
                  for r in history if float(r["coordinate"])==v]
            points[v]=take(rows)
            output.extend(rows)

    def fill(u,x,alt):
        ans,rows=profile_at(mode,name,u,x,alt)
        points[u]=ans
        output.extend(rows)
        return ans

    walk=x.copy()
    for v in sorted((u for u in grid if sign*(u-center)>1e-8),reverse=sign<0):
        # The first pass stores receipt rows that can be reused after reanchoring.
        ans=points[v] if v in points else fill(v,walk,True)
        if float(ans["all_constraint_margin"])>=0.:
            walk=c.anchor(ans)
    for _ in range(7):
        s=[center]+[u for u,a in points.items() if float(a["all_constraint_margin"])>=0.
                    and float(a["lr"])<=THRESHOLD]
        edge=(min if sign<0 else max)(s)
        outside=[u for u in points if sign*(u-edge)>0.]
        if not outside:
            break
        q=min(outside,key=lambda u:abs(u-edge))
        if abs(q-edge)<.0016:
            break
        start=x if edge==center else c.anchor(points[edge])
        fill(float((edge+q)*.5),start,False)
    return output


def summary(mode,rows):
    k,c,_=case(mode)
    r,x=ref(mode)
    best=float(r["loss"])
    feasible=[p for p in rows+read(HERE/f"fits_{mode}.tsv")
              if float(p["all_constraint_margin"])>=0]
    winner=min(feasible,key=lambda a:float(a["loss"]))
    if float(winner["loss"])<best-.025:
        return None,winner
    support=[p for p in feasible if float(p["loss"])<=best+THRESHOLD]
    out=[]
    for name in TARGETS:
        j=c.names.index(name)
        center=x[j]
        xs=[center]+[float(a[f"x_{name}"]) for a in support]
        ends=[]
        labels=[]
        for sign,s in ((-1,min(xs)),(1,max(xs))):
            box=0. if sign<0 else 1.
            if abs(s-box)<1e-8:
                ends.append(box); labels.append("box"); continue
            other=[p for p in rows if p["parameter"]==name and stable(p)
                   and sign*(float(p["coordinate"])-s)>0. and float(p["loss"])>best+THRESHOLD]
            q=min(other,key=lambda a:abs(float(a["coordinate"])-s)) if other else None
            if q and abs(float(q["coordinate"])-s)<.006:
                inner=min((float(a["loss"])-best for a in support
                           if abs(float(a[f"x_{name}"])-s)<1e-9),default=0.)
                cross=s+(float(q["coordinate"])-s)*(THRESHOLD-inner)/(float(q["loss"])-best-inner)
                ends.append(cross); labels.append("crossing")
            else:
                ends.append(math.nan); labels.append("unresolved")
        if max(center-min(xs),max(xs)-center)>.1000001:
            status="outside_supported"
        elif all(np.isfinite(ends)) and max(center-ends[0],ends[1]-center)<=.1000001:
            status="profile_pass"
        else:
            status="unresolved"
        convert_n=lambda u:c.unmap(np.where(np.arange(len(x))==j,u,x))[0][c.active[j]]
        out.append(dict(mode=mode,parameter=name,estimate=convert_n(center),x_estimate=center,
                        lower=convert_n(ends[0]),upper=convert_n(ends[1]),x_lower=ends[0],
                        x_upper=ends[1],x_min_supported=min(xs),x_max_supported=max(xs),
                        lower_type=labels[0],upper_type=labels[1],precision_gate=status))
    return out,None


def profiles(modes,workers):
    for attempt in range(3):
        if not modes:
            return
        result={mode:[] for mode in modes}
        with pool(workers) as ex:
            jobs={ex.submit(side,(m,t,d)):(m,t,d) for m in modes for t in TARGETS for d in (-1,1)}
            for f in as_completed(jobs):
                m,t,d=jobs[f]
                rows=f.result()
                result[m].extend(rows)
                if rows:
                    write(HERE/f"grid_{m}_{t}_{d}.tsv",rows)
                print("profile",m,t,d,len(rows),flush=True)
        remaining=[]
        for m in modes:
            rows=result[m]
            write(HERE/f"profile_points_{m}.tsv",rows)
            ans,better=summary(m,rows)
            if better is not None:
                k,c,_=case(m)
                old=read(HERE/f"best_{m}.tsv")[0]
                a=[run_fit((m,f"profile_restart_{i}",c.anchor(better),{},132,s))
                   for i,s in enumerate((30.,16.))]
                a.append(old)
                chosen=take(a)
                assert stable(chosen) and float(chosen["loss"])<=float(better["loss"])+.015
                write(HERE/f"fits_{m}.tsv",read(HERE/f"fits_{m}.tsv")+a)
                write(HERE/f"best_{m}.tsv",[chosen])
                print("new reference",m,chosen["loss"],flush=True)
                remaining.append(m)
            else:
                write(HERE/f"intervals_{m}.tsv",ans)
                for r in ans:
                    print(r,flush=True)
        modes=remaining
        ref.cache_clear()
        profile_seeds.cache_clear()
    if modes:
        raise RuntimeError(f"Could not stabilize the reference for {modes}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("action",choices=["fit","profile","both"])
    parser.add_argument("modes",nargs="+")
    parser.add_argument("--workers",type=int,default=9)
    a=parser.parse_args()
    if a.action!="profile":
        fits(a.modes,a.workers)
    if a.action!="fit":
        profiles(a.modes,a.workers)
