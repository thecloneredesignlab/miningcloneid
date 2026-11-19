"""Profile the ten nominated mechanisms while optimizing the remaining AR(1) model."""
import argparse
from collections import defaultdict
from concurrent.futures import as_completed

import numpy as np

from ar1_model import (COUNTS, HERE, STATIONARY, TARGETS, THRESHOLD, Config, choose,
                       likelihood_ratio, model, pool, read, receipt, row_x, solve, write)
from inference import lambda_supply_floor
from refit_ar1 import refine


def reference():
    return read(HERE/"observed_best.tsv")[0]


def old_bank(c):
    locations=(COUNTS/"real_fits.tsv",COUNTS/"profiles_direct.tsv",COUNTS/"profiles_edges.tsv",
               COUNTS/"profiles_cross.tsv",COUNTS/"profiles_polish.tsv",COUNTS/"profiles_supplement.tsv")
    items=[]
    for path in locations:
        if path.exists():
            items.extend(r for r in read(path) if float(r["min_aggregate_donor_margin_fraction"])>=0)
    return np.unique(np.array([c.load_old(r) for r in items]),axis=0)


def checked_solve(k,c,target,fixed,starts,limits=(48,60,72),divisors=(12.,5.,28.)):
    results=[]
    for anchor,limit,divisor in zip(starts,limits,divisors):
        f=solve(k,c,target,anchor,limit=limit,fixed=fixed,divisor=divisor)
        r=receipt(k,c,target,f)
        results.append(r)
        stable=sorted((z for z in results if float(z["donor_margin"])>=0
                       and float(z["stationarity"])<=STATIONARY),key=lambda z:float(z["loss"]))
        if (len(stable)>=2 and
            likelihood_ratio(float(stable[1]["loss"]),float(stable[0]["loss"]),k.n)<=.06):
            break
    return choose(results,k.n)


def append_labels(row,name,value,test):
    return dict(parameter=name,coordinate=value,test=test,**row)


def fit_profile(task):
    name,base,application=task
    k=model()
    c=Config(k)
    j=c.names.index(name)
    center=float(base[f"x_{name}"])
    x0=row_x(c,base)
    ref=float(base["loss"])
    points={}
    f0={q:v for q,v in base.items() if q!="origin"}
    f0["success"]=True
    points[center]=append_labels(f0,name,center,"unrestricted_reference")
    grid={0.,.25,.50,.75,1.,max(0,center-.10001),min(1,center+.10001),
          max(0,center-.04),min(1,center+.04)}

    def add(value,warm):
        fixed={j:value}
        if name=="lam_max":
            z=warm.copy()
            z[j]=value
            if c.unmap(z)[0][0]<lambda_supply_floor(k)-1e-8:
                invalid=dict(x=z,loss=float("inf"),supply=-1.,gradient=float("nan"),
                             calls=0,iterations=0,success=False,termination="insufficient supply by growth bound",
                             rounds=0,seconds=0.)
                r=receipt(k,c,k.y,invalid)
                r.update(loss=float("inf"),donor_margin=-1.)
                points[value]=append_labels(r,name,value,"analytical_infeasible")
                return
        r=checked_solve(k,c,k.y,fixed,[warm])
        lr=likelihood_ratio(float(r["loss"]),ref,k.n)
        if float(r["donor_margin"])<0 or float(r["stationarity"])>STATIONARY or lr>2.8:
            shape=application[np.argmin(abs(application[:,j]-value))]
            alt=x0 if np.linalg.norm(warm-x0)>1e-6 else shape
            r1=checked_solve(k,c,k.y,fixed,[alt],limits=(61,),divisors=(6.,))
            r=choose([r,r1],k.n)
        points[value]=append_labels(r,name,value,"direct_constrained_fit")

    for sign in (-1,1):
        warm=x0
        for value in sorted((v for v in grid if (v-center)*sign>1e-10),reverse=sign<0):
            add(value,warm)
            if float(points[value]["donor_margin"])>=0:
                warm=row_x(c,points[value])
    for _ in range(5):
        accepted=[v for v,r in points.items() if float(r["donor_margin"])>=0 and
                  likelihood_ratio(float(r["loss"]),ref,k.n)<=THRESHOLD]
        todo=[]
        if not accepted:
            break
        for side,edge in ((-1,min(accepted)),(1,max(accepted))):
            other=[v for v in points if (v-edge)*side>0]
            if other:
                rejection=min(other,key=lambda v:abs(edge-v))
                if abs(edge-rejection)>.003:
                    todo.append(((edge+rejection)/2,row_x(c,points[edge])))
        if not todo:
            break
        for v,warm in todo:
            add(v,warm)
    return list(points.values())


def load_profiles():
    records=read(HERE/"profile_direct.tsv")
    for f in sorted(HERE.glob("profile_refined_*.tsv")):
        records.extend(read(f))
    return records


def envelope():
    k=model()
    c=Config(k)
    base=reference()
    center=row_x(c,base)
    ref=float(base["loss"])
    fitted=load_profiles()
    allowed=[r for r in fitted if float(r["donor_margin"])>=0]
    assert {r["parameter"] for r in fitted}==set(TARGETS)
    d=defaultdict(list)
    for r in fitted:
        d[r["parameter"],float(r["coordinate"])].append(r)
    direct=defaultdict(dict)
    for (p,x),seq in d.items():
        direct[p][x]=choose(seq,k.n)
    supports=defaultdict(dict)
    fitting_cases=allowed+read(HERE/"observed_refits.tsv")+[base]
    for r in fitting_cases:
        if float(r["donor_margin"])<0:
            continue
        lr=likelihood_ratio(float(r["loss"]),ref,k.n)
        if lr>THRESHOLD:
            continue
        for p in TARGETS:
            x=float(r[f"x_{p}"])
            if x not in supports[p] or lr<supports[p][x][0]:
                supports[p][x]=(lr,r)
    summary=[]
    curves=[]
    for p in TARGETS:
        j=c.names.index(p)
        a=center[j]
        points=supports[p]
        lo,hi=min(points),max(points)
        endpoints=[]
        state=[]
        for side,edge in ((-1,lo),(1,hi)):
            if (side==-1 and edge<=1e-8) or (side==1 and edge>=1-1e-8):
                endpoints.append(edge)
                state.append("registered_box_endpoint")
                continue
            alternatives=[]
            for at,r in direct[p].items():
                if (at-edge)*side<=0:
                    continue
                if r["test"]=="analytical_infeasible":
                    alternatives.append((at,float("inf"),"supply_bound"))
                elif float(r["donor_margin"])>=0:
                    lr=likelihood_ratio(float(r["loss"]),ref,k.n)
                    if lr>THRESHOLD:
                        alternatives.append((at,lr,"stationary_crossing" if float(r["stationarity"])<=STATIONARY
                                             else "approximate_crossing"))
            if not alternatives:
                endpoints.append(None)
                state.append("outer_crossing_unresolved")
                continue
            v,lr,kind=min(alternatives,key=lambda q:abs(q[0]-edge))
            left=points[edge][0]
            bound=v if kind=="supply_bound" else edge+(v-edge)*(THRESHOLD-left)/(lr-left)
            endpoints.append(bound)
            state.append(kind)

        def natural(v):
            if v is None:
                return ""
            z=center.copy()
            z[j]=v
            return float(c.unmap(z)[0][c.active[j]])
        beyond=max(a-lo,hi-a)>.1000000001
        inside=all(x is not None and abs(a-x)<=.10 for x in endpoints)
        verified=True
        for side in (-1,1):
            q=a+side*.10001
            if 0.<q<1.:
                relevant=[r for x,r in direct[p].items() if abs(x-q)<.012 and (x-a)*side>.088]
                verified=verified and any(r["test"]=="analytical_infeasible" or
                    (float(r["donor_margin"])>=0 and float(r["stationarity"])<=STATIONARY and
                     likelihood_ratio(float(r["loss"]),ref,k.n)>THRESHOLD) for r in relevant)
        status="outside_supported" if beyond else "profile_pass" if inside and verified else "unresolved"
        summary.append(dict(parameter=p,estimate=base[p],normalized_estimate=a,
                           lower_natural=natural(endpoints[0]),upper_natural=natural(endpoints[1]),
                           lower_normalized="" if endpoints[0] is None else endpoints[0],
                           upper_normalized="" if endpoints[1] is None else endpoints[1],
                           min_supported=lo,max_supported=hi,lower_type=state[0],
                           upper_type=state[1],profile_status=status))
        for x,r in direct[p].items():
            lr=likelihood_ratio(float(r["loss"]),ref,k.n) if float(r["donor_margin"])>=0 else float("inf")
            curves.append(dict(parameter=p,coordinate=x,lr=lr,source=r["test"],
                               feasible=float(r["donor_margin"])>=0,stationarity=r["stationarity"],
                               rho=r["rho"],sigma=r["sigma"]))
        for x,(lr,r) in points.items():
            if (x in direct[p] and float(direct[p][x]["donor_margin"])>=0
                    and likelihood_ratio(float(direct[p][x]["loss"]),ref,k.n)<=lr+1e-10):
                continue
            curves.append(dict(parameter=p,coordinate=x,lr=lr,source="feasible_profile_projection",
                              feasible=True,stationarity=r["stationarity"],rho=r["rho"],sigma=r["sigma"]))
    write(HERE/"profile_intervals.tsv",summary)
    write(HERE/"profile_envelope.tsv",sorted(curves,key=lambda r:(TARGETS.index(r["parameter"]),r["coordinate"])))
    return summary


def update_reference(workers):
    k=model()
    c=Config(k)
    base=reference()
    fitted=[r for r in load_profiles() if float(r["donor_margin"])>=0]
    best=min(fitted,key=lambda r:float(r["loss"]))
    gap=likelihood_ratio(float(best["loss"]),float(base["loss"]),k.n)
    if gap>=-.001:
        return
    anchor={name:best[name] for name in read(HERE/"observed_refits.tsv")[0] if name!="origin"}
    anchor=dict(origin="first_better_profile",**anchor)
    result=refine(anchor,min(workers,5),label="profile_restored")
    assert likelihood_ratio(float(best["loss"]),float(result["loss"]),k.n)>=-.06


def work(workers):
    k=model()
    c=Config(k)
    obs=reference()
    bank=old_bank(c)
    path=HERE/"profile_direct.tsv"
    current=read(path) if path.exists() else []
    done={r["parameter"] for r in current}
    jobs=[(p,obs,bank) for p in TARGETS if p not in done]
    if jobs:
        with pool(workers) as executor:
            for f in as_completed([executor.submit(fit_profile,q) for q in jobs]):
                rows=f.result()
                current.extend(rows)
                write(path,current)
                print("profile",rows[0]["parameter"],len(rows),flush=True)
    update_reference(workers)
    for r in envelope():
        print(r["parameter"],r["profile_status"],r["lower_normalized"],r["upper_normalized"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=10)
    parser.add_argument("--envelope",action="store_true")
    args=parser.parse_args()
    if args.envelope:
        envelope()
    else:
        work(args.workers)
