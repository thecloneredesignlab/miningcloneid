"""One-dimensional profiles and two-dimensional likelihood support sets."""
import argparse
from concurrent.futures import as_completed
from functools import lru_cache
from itertools import combinations,product
import math

import numpy as np

from modal import HERE,STAGE3,TARGETS,THRESHOLD,case,fixed_for,pool,read,stable,write
from fit import PRIMARY,run_one,take

TOL=.10
PAIR_LR=5.991464547108


@lru_cache
def base(mode):
    k,c,_=case(mode)
    r=read(HERE/f"best_{mode}.tsv")[0]
    return r,c.anchor(r)


@lru_cache
def points_bank(mode):
    k,c,y=case(mode)
    rows=read(HERE/f"fits_{mode}.tsv")
    rows+=read(STAGE3/"profile_points_autosomes_joint.tsv")
    return np.array([c.anchor(r) for r in rows if float(r["donor_margin"])>=0.])


def profile_side(job):
    mode,name,side=job
    k,c,y=case(mode)
    ref,x=base(mode)
    j=c.names.index(name)
    center=x[j]
    bank=points_bank(mode)
    arch=next(r for r in read(STAGE3/"profile_intervals_autosomes_joint.tsv") if r["parameter"]==name)
    source=[float(arch["lower_normalized"]),float(arch["upper_normalized"])]
    regular={0.,.25,.5,.75,1.,*source,float(np.clip(center+side*(TOL+1e-5),0,1))}
    points={}
    rows=[]

    def fit(u,warm,extra):
        starts=[("carry",warm)]
        if extra:
            alternatives=[]
            for i in np.argsort(abs(bank[:,j]-u))[:12]:
                z=bank[i].copy()
                z[j]=u
                if np.linalg.norm(np.delete(z-warm,j))<.035:
                    continue
                ar,v,*_=k.parts(c.unmap(z)[0],y)
                alternatives.append((v+225.*max(0.,-ar[1]),z))
            if alternatives:
                starts.append(("alternate",min(alternatives,key=lambda w:w[0])[1]))
        solved=[]
        for tag,z in starts:
            o=run_one((mode,tag,z,{j:u},65,42.))
            o.update(parameter=name,coordinate=u,side=side,lr=o["loss"]-float(ref["loss"]))
            solved.append(o)
        best=take(solved)
        if not stable(best) and best["lr"]>THRESHOLD:
            o=run_one((mode,"precise",c.anchor(best),{j:u},105,18.))
            o.update(parameter=name,coordinate=u,side=side,lr=o["loss"]-float(ref["loss"]))
            solved.append(o)
            best=take(solved)
        points[u]=best
        rows.extend(solved)
        return best

    warm=x.copy()
    for u in sorted((v for v in regular if side*(v-center)>1e-8),reverse=side<0):
        r=fit(u,warm,True)
        if r["donor_margin"]>=0:
            warm=c.anchor(r)
    for _ in range(7):
        support=[center]+[u for u,r in points.items() if r["lr"]<=THRESHOLD and r["donor_margin"]>=0.]
        edge=(min if side<0 else max)(support)
        far=[u for u in points if side*(u-edge)>0]
        if not far:
            break
        rejected=min(far,key=lambda q:abs(q-edge))
        if abs(rejected-edge)<=.0025:
            break
        fit((rejected+edge)/2, x if edge==center else c.anchor(points[edge]),False)
    return rows


def summarize(mode,rows):
    k,c,y=case(mode)
    ref,x=base(mode)
    ref_loss=float(ref["loss"])
    admissible=[r for r in rows+read(HERE/f"fits_{mode}.tsv") if float(r["donor_margin"])>=0]
    if min(float(r["loss"]) for r in admissible)+.025<ref_loss:
        raise RuntimeError(f"{mode} profile found a better free-fit branch")
    kept=[r for r in admissible if float(r["loss"])<=ref_loss+THRESHOLD]
    answer=[]
    for name in TARGETS:
        j=c.names.index(name)
        a=x[j]
        inner=[a]+[float(r[f"x_{name}"]) for r in kept]
        low,high=min(inner),max(inner)
        bounds=[]
        labels=[]
        for side,v in [(-1,low),(1,high)]:
            bound=0. if side<0 else 1.
            if abs(v-bound)<1e-7:
                bounds.append(bound); labels.append("box"); continue
            candidates=[r for r in rows if r["parameter"]==name and stable(r) and
                 side*(float(r["coordinate"])-v)>0 and float(r["lr"])>THRESHOLD]
            if candidates:
                outsider=min(candidates,key=lambda r:abs(float(r["coordinate"])-v))
                adjacent=min((float(r["loss"])-ref_loss for r in kept
                              if abs(float(r[f"x_{name}"])-v)<1e-10),default=0.)
                q=float(outsider["coordinate"])
                gap=abs(q-v)
                b=v+(q-v)*(THRESHOLD-adjacent)/(float(outsider["lr"])-adjacent)
                if gap<=.015:
                    bounds.append(b); labels.append("crossing"); continue
            bounds.append(math.nan); labels.append("unresolved")
        if max(a-low,high-a)>TOL+1e-7:
            status="outside_supported"
        elif all(np.isfinite(bounds)) and max(a-bounds[0],bounds[1]-a)<=TOL+1e-7:
            status="profile_pass"
        else:
            status="unresolved"
        convert=lambda q: c.unmap(np.where(np.arange(len(x))==j,q,x))[0][c.active[j]]
        answer.append(dict(mode=mode,parameter=name,estimate=convert(a),x_estimate=a,
            lower=convert(bounds[0]),upper=convert(bounds[1]),x_lower=bounds[0],
            x_upper=bounds[1],x_min_supported=low,x_max_supported=high,
            lower_type=labels[0],upper_type=labels[1],precision_gate=status))
    return answer


def one_dimension(modes,workers):
    for m in modes:
        if (HERE/f"intervals_{m}.tsv").exists():
            continue
        rows=[]
        with pool(workers) as ex:
            jobs={ex.submit(profile_side,(m,p,s)):(p,s) for p,s in product(TARGETS,(-1,1))}
            for f in as_completed(jobs):
                p,s=jobs[f]
                part=f.result()
                if part:
                    write(HERE/f"grid_{m}_{p}_{s}.tsv",part)
                rows+=part
                print(m,p,s,len(part),flush=True)
        write(HERE/f"profile_points_{m}.tsv",rows)
        interval=summarize(m,rows)
        write(HERE/f"intervals_{m}.tsv",interval)
        print(interval,flush=True)


def pair_point(job):
    mode,pa,pb,a,b,starts=job
    k,c,y=case(mode)
    ref,_=base(mode)
    i,j=c.names.index(pa),c.names.index(pb)
    fits=[]
    for at,z in enumerate(starts):
        r=run_one((mode,f"pair_{at}",z,{i:a,j:b},56,44.))
        fits.append(r)
        if at==0 and stable(r) and r["loss"]-float(ref["loss"])<PAIR_LR-.05:
            break
    r=take(fits)
    return dict(parameter_a=pa,parameter_b=pb,coordinate_a=a,coordinate_b=b,
                lr=float(r["loss"])-float(ref["loss"]),**r)


def mesh_grid(mode,pa,pb):
    k,c,y=case(mode)
    ref,x=base(mode)
    prof=read(HERE/f"profile_points_{mode}.tsv")
    i,j=c.names.index(pa),c.names.index(pb)
    viable=[r for r in prof if stable(r)]
    rows=[r for r in prof if float(r["donor_margin"])>=0. and float(r["lr"])<=PAIR_LR+.1]
    shape=[(x[i],x[j])]+[(float(r[f"x_{pa}"]),float(r[f"x_{pb}"])) for r in rows]
    shape=np.asarray(shape)
    lo=np.maximum(0.,shape.min(axis=0)-.045)
    hi=np.minimum(1.,shape.max(axis=0)+.045)
    fixes=[np.unique(np.r_[np.linspace(lo[q],hi[q],5),x[[i,j][q]]]) for q in range(2)]
    bank=np.array([c.anchor(r) for r in viable]) if viable else x[None,:]
    for a,b in product(*fixes):
        if abs(a-x[i])+abs(b-x[j])<1e-12:
            continue
        nearest=np.argmin((bank[:,i]-a)**2+(bank[:,j]-b)**2)
        if np.linalg.norm(bank[nearest]-x)<.035:
            starts=[x]
        else:
            starts=[bank[nearest],x]
        yield mode,pa,pb,float(a),float(b),starts


def pairs(modes,workers):
    options=list(combinations(TARGETS,2))
    scheduled=[j for m in modes for a,b in options for j in mesh_grid(m,a,b)]
    out={m:[] for m in modes}
    with pool(workers) as ex:
        work={ex.submit(pair_point,j):j for j in scheduled}
        for f in as_completed(work):
            mode=work[f][0]
            out[mode].append(f.result())
            if len(out[mode])%15==0:
                total=sum(len(v) for v in out.values())
                print(f"{total}/{len(scheduled)} completed",flush=True)
                write(HERE/f"pair_grid_{mode}.tsv",out[mode])
    for m,rows in out.items():
        rows.sort(key=lambda r:(r["parameter_a"],r["parameter_b"],r["coordinate_a"],r["coordinate_b"]))
        write(HERE/f"pair_grid_{m}.tsv",rows)
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("analysis",choices=["profile","pair"])
    p.add_argument("modes",nargs="*",default=PRIMARY)
    p.add_argument("--workers",type=int,default=9)
    args=p.parse_args()
    (one_dimension if args.analysis=="profile" else pairs)(args.modes,args.workers)
