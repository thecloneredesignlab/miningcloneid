"""Check and refine interior likelihood-threshold crossings."""
from concurrent.futures import as_completed
from collections import defaultdict

import numpy as np

from branch_influence import MODES, OUT, THRESHOLD, masked, pool, read, solution, stable, tag, take, write


def threshold_points(mode,mask):
    label=tag(mode,mask)
    base=read(OUT/f"profile_ref_{label}.tsv")[0]
    f=float(base["loss"])
    x=float(base["x_buffer_beta"])
    values=defaultdict(list)
    for r in read(OUT/f"profile_points_{label}.tsv"):
        values[float(r["coordinate"])].append(r)
    points={v:take(rs) for v,rs in values.items()}
    points[x]=base
    def is_supported(row):
        return float(row["all_constraint_margin"])>=0 and float(row["loss"])-f<=THRESHOLD
    ordered=[(v,r,is_supported(r)) for v,r in sorted(points.items())]
    return [(a,b) for a,b in zip(ordered,ordered[1:]) if a[-1]!=b[-1]]


def solve_crossing(mode,mask,a,b):
    x,ra,sa=a
    z,rb,sb=b
    k,c,_,_=masked(mode,mask)
    base=read(OUT/f"profile_ref_{tag(mode,mask)}.tsv")[0]
    best=float(base["loss"])
    audit=[]
    for _ in range(9):
        if z-x<=.0018:
            break
        t=(x+z)/2
        rows=[solution(mode,mask,c.anchor(r),f"internal_cross_{i}",t)
              for i,r in enumerate((ra,rb))]
        p=take(rows)
        if not stable(p):
            rows.append(solution(mode,mask,c.anchor(p),"internal_cross_polish",t,True))
            p=take(rows)
        for r in rows:
            r["lr"]=float(r["loss"])-best
        audit+=rows
        if (float(p["all_constraint_margin"])>=0. and p["lr"]<=THRESHOLD)==sa:
            x,ra=t,p
        else:
            z,rb=t,p
    assert stable(ra) and stable(rb)
    fl,fu=(float(r["loss"])-best for r in (ra,rb))
    point=x+(z-x)*(THRESHOLD-fl)/(fu-fl)
    p=c.anchor(base)
    p[c.names.index("buffer_beta")]=point
    beta=c.unmap(p)[0][c.active[c.names.index("buffer_beta")]]
    return dict(mode=mode,mask=mask,x_lower_bracket=x,x_upper_bracket=z,x_interpolated=point,
                beta=beta,inside_before=sa,inside_after=sb),audit


if __name__=="__main__":
    query=[]
    for m in MODES:
        for mask in range(1,15):
            flips=threshold_points(m,mask)
            if len(flips)>2:
                query.extend((m,mask,a,b) for a,b in flips[1:-1])
    if query:
        results=[]
        audit=[]
        with pool(min(len(query),8)) as ex:
            for f in as_completed(ex.submit(solve_crossing,*q) for q in query):
                row,raw=f.result()
                results.append(row)
                audit+=raw
                print(row,flush=True)
        write(OUT/"internal_support_crossings.tsv",sorted(results,key=lambda r:(r["mode"],r["mask"],r["beta"])))
        write(OUT/"internal_crossing_fits.tsv",audit)
