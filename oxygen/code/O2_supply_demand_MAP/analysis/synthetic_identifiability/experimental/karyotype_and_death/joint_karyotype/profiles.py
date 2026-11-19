"""Profile five parameters over the registered box with independent nuisance starts."""
import argparse
from concurrent.futures import as_completed
from functools import lru_cache
import math

import numpy as np

from joint_model import (BASE, HERE, TARGETS, THRESHOLD, likelihood_ratio,
                         pool, read, row_x, stable, write)
from refit import case, task_fit

COUNT_BASE = {r["parameter"]:r for r in read(BASE/"practical_verdicts.tsv")}
TOLERANCE = .10


@lru_cache
def previous_vectors(canonical_mode):
    _,c,_ = case(canonical_mode)
    data = [r for p in ("profile_direct", "profile_refined_tolerance",
                        "observed_refits") for r in read(BASE/f"{p}.tsv")
            if float(r["donor_margin"]) >= 0.]
    return np.array([row_x(c,r) for r in data])


def profile_side(job):
    mode,name,sign = job
    k,c,y = case(mode)
    best = read(HERE/f"best_{mode}.tsv")[0]
    x0 = row_x(c,best)
    j = c.names.index(name)
    center = x0[j]
    reference = float(best["loss"])
    points = {}
    saved = []
    bank = previous_vectors(mode)
    branches = [
        row_x(c,r) for r in read(HERE/f"refits_{mode}.tsv")
        if stable(r) and np.linalg.norm(row_x(c,r)-x0) > .07
    ]
    alternative = branches[:3]

    def evaluate(u, walk, extra):
        if u in points:
            return points[u]
        anchors = [("continuation",walk)]
        if extra:
            axis = np.argsort(abs(bank[:,j]-u))
            candidates = [bank[t] for t in axis[:7]]
            candidates += alternative + [x0]
            # Prefer a feasible point from a different combination of the free
            # coordinates. Its objective is evaluated at the fixed target.
            choices = []
            for v in candidates:
                v = v.copy()
                v[j] = u
                if np.linalg.norm(np.delete(v-walk,j)) < .045:
                    continue
                out, _, f0, _ = k.evaluate_joint(c.unmap(v)[0],y)
                choices.append((f0*(1.+max(0.,-out[1]))**2,v))
            if choices:
                anchors.append(("alternate",min(choices,key=lambda z:z[0])[1]))
        for origin, start in anchors:
            r = task_fit((mode,origin,start,{j:u},65,32.))
            r.update(parameter=name,coordinate=u,side=sign,
                     lr=likelihood_ratio(r["loss"],reference,k.n))
            saved.append(r)
        rows = saved[-len(anchors):]
        # When a side of the confidence interval is estimated, keep the best
        # feasible loss. Report convergence of the adjacent excluded point.
        allowed = [r for r in rows if r["donor_margin"]>=0.]
        r = min(allowed or rows,key=lambda z:z["loss"])
        if not stable(r) and r["lr"] > THRESHOLD:
            q=task_fit((mode,"tight_refit",row_x(c,r),{j:u},85,16.))
            q.update(parameter=name,coordinate=u,side=sign,
                     lr=likelihood_ratio(q["loss"],reference,k.n))
            saved.append(q)
            if q["donor_margin"]>=0 and (
                r["donor_margin"]<0 or q["loss"]<r["loss"] or
                (stable(q) and abs(q["lr"]-r["lr"])<.005)
            ):
                r = q
        points[u] = r
        return r

    grid = {0.,.25,.5,.75,1.,float(np.clip(center+sign*(TOLERANCE+1e-5),0,1))}
    old_edge = float(COUNT_BASE[name]["lower_normalized" if sign<0 else "upper_normalized"])
    if abs(old_edge-center)>.025:
        grid.add(old_edge)
    # Search the whole registered side, including beyond the first crossing.
    walk = x0
    for u in sorted((u for u in grid if (u-center)*sign>1e-7),reverse=sign<0):
        r = evaluate(u,walk,True)
        if r["donor_margin"]>=0.:
            walk = row_x(c,r)
    for _ in range(6):
        accepted = [center]+[u for u,r in points.items()
                             if r["lr"]<=THRESHOLD and r["donor_margin"]>=0.]
        edge = min(accepted) if sign<0 else max(accepted)
        beyond = [u for u in points if (u-edge)*sign>0.]
        if not beyond:
            break
        rejected = min(beyond,key=lambda u:abs(u-edge))
        if abs(edge-rejected)<=.002:
            break
        walk = x0 if edge==center else row_x(c,points[edge])
        evaluate((edge+rejected)*.5,walk,False)
    return saved


def summarize(mode,rows):
    k,c,_ = case(mode)
    base = read(HERE/f"best_{mode}.tsv")[0]
    extra_available = HERE/f"refits_{mode}.tsv"
    reference=float(base["loss"])
    feasible=[r for r in rows+read(extra_available) if float(r["donor_margin"])>=0.]
    if any(float(r["loss"])<reference and
           likelihood_ratio(reference,float(r["loss"]),k.n)>.025 for r in feasible):
        return None
    supported=[r for r in feasible if likelihood_ratio(float(r["loss"]),reference,k.n)<=THRESHOLD]
    result=[]
    for name in TARGETS:
        j=c.names.index(name)
        center=float(base[f"x_{name}"])
        coords=[center]+[float(r[f"x_{name}"]) for r in supported]
        low,high=min(coords),max(coords)
        ends=[]
        statuses=[]
        for sign,edge in [(-1,low),(1,high)]:
            box = 0. if sign<0 else 1.
            if abs(box-edge)<=1e-7:
                ends.append(box)
                statuses.append("registered_box")
                continue
            excluded=[r for r in rows if r["parameter"]==name and stable(r) and
                      (float(r["coordinate"])-edge)*sign>0. and float(r["lr"])>THRESHOLD]
            outside=min(excluded,key=lambda r:abs(float(r["coordinate"])-edge)) if excluded else None
            ok=False
            if outside:
                gap=abs(float(outside["coordinate"])-edge)
                inner=[r for r in supported if abs(float(r[f"x_{name}"])-edge)<1e-8]
                a=min([likelihood_ratio(float(r["loss"]),reference,k.n) for r in inner],
                      default=0.)
                b=float(outside["lr"])
                end=edge+(float(outside["coordinate"])-edge)*(THRESHOLD-a)/(b-a)
                ok=gap<=.005
            ends.append(end if ok else math.nan)
            statuses.append("stationary_crossing" if ok else "unresolved")
        x0=row_x(c,base)

        def natural(q):
            x=x0.copy()
            x[j]=q
            return c.unmap(x)[0][c.active[j]]

        out = center-low>TOLERANCE+1e-8 or high-center>TOLERANCE+1e-8
        if out:
            status="outside_supported"
        elif (np.all(np.isfinite(ends)) and center-ends[0]<=TOLERANCE+1e-8 and
              ends[1]-center<=TOLERANCE+1e-8):
            status="profile_pass"
        else:
            status="unresolved"
        result.append(dict(mode=mode,parameter=name,estimate=natural(center),
                           normalized_estimate=center,lower_natural=natural(ends[0]),
                           upper_natural=natural(ends[1]),lower_normalized=ends[0],
                           upper_normalized=ends[1],lower_supported=low,upper_supported=high,
                           lower_type=statuses[0],upper_type=statuses[1],
                           profile_status=status))
    return result


def run(modes,workers):
    all_rows = {m:[] for m in modes}
    jobs = [(m,n,s) for m in modes for n in TARGETS for s in (-1,1)]
    with pool(workers) as executor:
        futures={executor.submit(profile_side,j):j for j in jobs}
        for f in as_completed(futures):
            mode,name,sign=futures[f]
            part=f.result()
            if part:
                write(HERE/f"profile_{mode}_{name}_{'lower' if sign<0 else 'upper'}.tsv",part)
            all_rows[mode].extend(part)
            print(mode,name,sign,len(part),min((r["lr"] for r in part),default="boundary"),
                  flush=True)
    for mode,rows in all_rows.items():
        rows.sort(key=lambda r:(TARGETS.index(r["parameter"]),r["coordinate"],r["loss"]))
        write(HERE/f"profile_points_{mode}.tsv",rows)
        result=summarize(mode,rows)
        if result is None:
            print(mode,"refit required before reporting the profile",flush=True)
        else:
            write(HERE/f"profile_intervals_{mode}.tsv",result)
            for r in result:
                print(mode,r["parameter"],r["lower_natural"],r["upper_natural"],r["profile_status"],
                      flush=True)


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("modes",nargs="+")
    p.add_argument("--workers",type=int,default=7)
    a=p.parse_args()
    run(a.modes,a.workers)
