"""One-at-a-time profile likelihoods for the experimental counts."""
import argparse
import concurrent.futures as cf
from collections import defaultdict

import numpy as np

from count_model import HERE, model, optimize
from fit_counts import record_real, solve_real, sorted_feasible
from inference import SLACK, THRESHOLD, lambda_supply_floor, likelihood_ratio, pool
from practical import Configuration, MECHANISMS, read, result_row, row_x, write


def close_fit(a,b,n):
    return abs(likelihood_ratio(a["loss"],b["loss"],n))<=.10


def fixed_solve(k,c,target,i,value,starts,limit=50):
    tests=[]
    for m,start in enumerate(starts):
        f=optimize(k,c,target,start,limit,fixed=(i,value),divisor=10. if m % 2 == 0 else 2.)
        tests.append(f)
        if len(tests)>1 and f["supply"]>=-SLACK and f["gradient"]<=.025:
            stable = [r for r in tests[:-1] if r["supply"]>=-SLACK and r["gradient"]<=.025]
            if any(close_fit(f,t,k.n) for t in stable):
                break
    valid=[r for r in tests if r["supply"]>=-SLACK]
    best=min(valid or tests,key=lambda r:r["loss"])
    if not valid:
        return best
    settled=[f for f in valid if f["gradient"]<=.025
             and likelihood_ratio(f["loss"],best["loss"],k.n)<=.001]
    return min(settled,key=lambda r:r["loss"]) if settled else best


def worker(task):
    name, base, other = task
    k=model()
    c=Configuration(k)
    j=MECHANISMS.index(name)
    center=float(base[f"x_{name}"])
    x0 = row_x(c,base)
    xother = row_x(c,other) if other is not None else x0
    minimum=float(base["loss"])
    points={}
    grid={0.,.25,.5,.75,1.,max(0,center-.10001),min(1,center+.10001)}
    grid.update((max(0,center-.04),min(1,center+.04)))
    points[center]=dict(x=x0,loss=minimum,supply=float(base["min_aggregate_donor_margin_fraction"]),
                        gradient=float(base["max_projected_gradient"]),
                        calls=0,iteration=0,scipy_success=True,termination="unrestricted reference",seconds=0.)

    def add(value,warm):
        if name=="lam_max":
            z=warm.copy()
            z[j]=value
            if c.unmap(z)[0][0] < lambda_supply_floor(k)-1e-8:
                f=dict(x=z,loss=float("inf"),supply=-float("inf"),gradient=float("nan"),calls=0,
                       iteration=0,scipy_success=False,termination="analytical donor-supply impossibility",seconds=0.)
                points[value]=f
                return f
        f=fixed_solve(k,c,k.y,j,value,[warm],limit=45)
        if f["supply"] < -SLACK or f["gradient"]>.025 or likelihood_ratio(f["loss"],minimum,k.n)>2.9:
            alt=x0 if np.linalg.norm(warm-x0)>1e-6 else xother
            g=fixed_solve(k,c,k.y,j,value,[alt],limit=54)
            if g["supply"]>=-SLACK and (f["supply"]<-SLACK or g["loss"]<f["loss"]):
                f=g
        points[value]=f
        return f

    for sign in (-1,1):
        target=sorted((v for v in grid if sign*(v-center)>1e-10),reverse=sign<0)
        warm=x0
        for v in target:
            f=add(v,warm)
            if f["supply"]>=-SLACK:
                warm=f["x"]
    # Refine the outer crossings of the observed confidence envelope.
    for _ in range(5):
        accepted=[v for v,f in points.items() if f["supply"]>=-SLACK
                  and likelihood_ratio(f["loss"],minimum,k.n)<=THRESHOLD]
        if not accepted:
            break
        todo=[]
        for sign,v in ((-1,min(accepted)),(1,max(accepted))):
            opposite=[w for w in points if sign*(w-v)>0]
            if opposite:
                w=min(opposite,key=lambda w:abs(w-v))
                if abs(w-v)>.003:
                    todo.append(((v+w)/2,points[v]["x"]))
        if not todo:
            break
        for v,warm in todo:
            add(v,warm)
    return [result_row(c,f,parameter=name,coordinate=x,feasible=f["supply"]>=-SLACK,
                       analytical_infeasible=f["termination"].startswith("analytical"))
            for x,f in sorted(points.items())]


def load_direct():
    rows=read(HERE/"profiles_direct.tsv")
    for suffix in ("edges","supplement","cross","polish"):
        path=HERE/f"profiles_{suffix}.tsv"
        if path.exists():
            rows += read(path)
    return rows


def envelope():
    best=read(HERE/"real_best.tsv")[0]
    loss=float(best["loss"])
    full=load_direct()
    raw=read(HERE/"real_fits.tsv")
    c=Configuration(model())
    names=MECHANISMS
    accepted=defaultdict(dict)
    direct=defaultdict(dict)
    by_key=defaultdict(list)
    for r in full:
        n,x=r["parameter"],float(r["coordinate"])
        by_key[n,x].append(r)
    for (n,x),variants in by_key.items():
        measured=[r for r in variants if r["feasible"]=="True"]
        floor=min(measured or variants,key=lambda r:float(r["loss"]))
        consistent=[r for r in measured if float(r["max_projected_gradient"])<=.025
                    and likelihood_ratio(float(r["loss"]),float(floor["loss"]),221)<=.001]
        direct[n][x]=min(consistent,key=lambda r:float(r["loss"])) if consistent else floor
    for r in full+raw+[best]:
        if float(r["min_aggregate_donor_margin_fraction"]) < -SLACK:
            continue
        lr=likelihood_ratio(float(r["loss"]),loss,221)
        if lr>THRESHOLD:
            continue
        for p in names:
            x=float(r[f"x_{p}"])
            if x not in accepted[p] or lr<accepted[p][x][0]:
                accepted[p][x]=(lr,r.get("parameter",r.get("origin","reference")))
    out=[]
    summaries=[]
    for j,p in enumerate(names):
        x=float(best[f"x_{p}"])
        low,hi=min(accepted[p]),max(accepted[p])
        bounds=[]
        evidence=[]
        for sign,edge in ((-1,low),(1,hi)):
            if (sign<0 and edge<=1e-9) or (sign>0 and edge>=1-1e-9):
                bounds.append(edge)
                evidence.append("bound")
                continue
            choices=[]
            for q,r in direct[p].items():
                if (q-edge)*sign<=0:
                    continue
                if r["analytical_infeasible"]=="True":
                    choices.append((q,float("inf"),"analytic_infeasibility"))
                elif r["feasible"]=="True":
                    lr=likelihood_ratio(float(r["loss"]),loss,221)
                    if lr>THRESHOLD:
                        quality="stationary_crossing" if float(r["max_projected_gradient"])<=.025 else "approximate_crossing"
                        choices.append((q,lr,quality))
            if not choices:
                bounds.append("")
                evidence.append("unresolved_outer_edge")
                continue
            q,lr,quality=min(choices,key=lambda a:abs(a[0]-edge))
            supported=accepted[p][edge][0]
            if quality=="analytic_infeasibility":
                bounds.append(q)
            else:
                bounds.append(edge+(q-edge)*(THRESHOLD-supported)/(lr-supported))
            evidence.append(quality)
        def natural(v):
            if v=="":
                return ""
            z=row_x(c,best)
            z[j]=v
            return c.unmap(z)[0][j]
        outside=max(x-low,hi-x)>.1000000001
        narrow=not outside and "" not in bounds and max(abs(v-x) for v in bounds)<=.1
        # Require checked exclusion past the tolerance on each open side.
        checked=True
        for sign in (-1,1):
            boundary=x+sign*.100000001
            if not 0<=boundary<=1:
                continue
            candidates=[r for value,r in direct[p].items() if abs(value-boundary)<.015
                        and (value-x)*sign>.085
                        and (r["analytical_infeasible"]=="True" or
                             (r["feasible"]=="True" and float(r["max_projected_gradient"])<=.025
                              and likelihood_ratio(float(r["loss"]),loss,221)>THRESHOLD))]
            checked=checked and bool(candidates)
        grade="outside_supported" if outside else "narrow_checked" if narrow and checked else "unresolved"
        summaries.append(dict(parameter=p,estimate=best[p],normalized_estimate=x,
                              lower_natural=natural(bounds[0]),upper_natural=natural(bounds[1]),
                              lower_normalized=bounds[0],upper_normalized=bounds[1],
                              lowest_accepted_normalized=low,highest_accepted_normalized=hi,
                              lower_endpoint_status=evidence[0],upper_endpoint_status=evidence[1],
                              profile_precision=grade))
        for q,r in direct[p].items():
            lr=(likelihood_ratio(float(r["loss"]),loss,221) if r["feasible"]=="True" else float("inf"))
            out.append(dict(parameter=p,coordinate=q,likelihood_ratio=lr,source="conditional_fit",
                            feasible=r["feasible"],stationarity=r["max_projected_gradient"],
                            termination=r["termination"]))
        for q,(lr,source) in accepted[p].items():
            r=direct[p].get(q)
            if r and float(r["loss"]) <= loss*np.exp(lr/221)+1e-11:
                continue
            out.append(dict(parameter=p,coordinate=q,likelihood_ratio=lr,source=source,
                            feasible=True,stationarity="",termination="feasible projection"))
    write(HERE/"profiles.tsv",sorted(out,key=lambda r:(r["parameter"],r["coordinate"])))
    write(HERE/"profile_intervals.tsv",summaries)
    return summaries


def refit_reference(workers):
    k=model()
    c=Configuration(k)
    rows=read(HERE/"real_fits.tsv")
    base=read(HERE/"real_best.tsv")[0]
    candidates=sorted_feasible(load_direct())
    chosen=[]
    for r in candidates:
        if float(r["loss"]) > float(base["loss"]) * np.exp(.02/k.n):
            break
        a=row_x(c,r)
        if all(np.linalg.norm(a-x)>.025 for x in chosen):
            chosen.append(a)
        if len(chosen)==5:
            break
    rng=np.random.default_rng(99102)
    if candidates and float(candidates[0]["loss"]) < float(base["loss"]):
        a=row_x(c,candidates[0])
        chosen += [np.clip(a+rng.normal(0,sd,len(a)),0,1) for sd in (.0002,.003)]
    if not chosen:
        return
    with pool(min(workers,len(chosen))) as work:
        tasks=[work.submit(solve_real,(f"profile_refit_{i}",a,125)) for i,a in enumerate(chosen,1)]
        for f in cf.as_completed(tasks):
            r=f.result()
            rows.append(r)
            print("reference refit",r["loss"],r["max_projected_gradient"],flush=True)
    write(HERE/"real_fits.tsv",sorted(rows,key=lambda r:float(r["loss"])))
    best_profile=candidates[0]
    if float(best_profile["loss"])<min(float(r["loss"]) for r in sorted_feasible(rows)):
        recovered={f:best_profile[f] for f in rows[0] if f!="origin"}
        rows.append(dict(origin="profile_feasible_reference",**recovered))
        write(HERE/"real_fits.tsv",sorted(rows,key=lambda r:float(r["loss"])))
    record_real(k,c,rows)


def run(workers):
    c=Configuration(model())
    base=read(HERE/"real_best.tsv")[0]
    others=sorted_feasible(read(HERE/"real_fits.tsv"))
    other=next((r for r in others if np.linalg.norm(row_x(c,r)-row_x(c,base))>.12),None)
    path=HERE/"profiles_direct.tsv"
    out=read(path) if path.exists() else []
    done={r["parameter"] for r in out}
    with pool(workers) as work:
        for f in cf.as_completed([work.submit(worker,(p,base,other)) for p in MECHANISMS if p not in done]):
            rows=f.result()
            out.extend(rows)
            print("profile",rows[0]["parameter"],len(rows),flush=True)
            write(path,out)
    refit_reference(workers)
    for r in envelope():
        print(r["parameter"],r["profile_precision"],r["lower_normalized"],r["upper_normalized"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=9)
    parser.add_argument("--envelope",action="store_true")
    args=parser.parse_args()
    if args.envelope:
        envelope()
    else:
        run(args.workers)
