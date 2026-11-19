"""Compare profile widths in the registered fitted coordinates and replay gate witnesses."""
import math

import numpy as np

from joint_model import BASE, HERE, TARGETS, THRESHOLD, likelihood_ratio, read, row_x, write
from refit import case

JOINT = "autosomes_joint"
COUNTS = "autosomes_count"
ARCHIVED = "archived_joint"


def index(path):
    return {r["parameter"]:r for r in read(path)}


def width(r):
    return float(r["upper_normalized"])-float(r["lower_normalized"])


def comparison():
    count=index(HERE/f"profile_intervals_{COUNTS}.tsv")
    joint=index(HERE/f"profile_intervals_{JOINT}.tsv")
    archive=index(HERE/f"profile_intervals_{ARCHIVED}.tsv")
    stage2=index(BASE/"practical_verdicts.tsv")
    recovery_path=HERE/"recovery_metrics.tsv"
    recovery={(r["parameter"],int(r["case_id"])):r for r in read(recovery_path)
              } if recovery_path.exists() else {}
    rows=[]
    for name in TARGETS:
        a,b,c,d=count[name],joint[name],stage2[name],archive[name]
        W=width(a)
        narrowed=(W-width(b))/W if W>0 else math.nan
        A,B=recovery.get((name,1)),recovery.get((name,2))
        confirmed=all(r and int(r["assigned"])==100 and float(r["rate"])>=.90
                      for r in (A,B))
        if b["profile_status"]=="outside_supported":
            decision="profile_outside_tolerance"
        elif confirmed:
            decision="precision_pass"
        elif all(r and int(r["assigned"])==100 for r in (A,B)):
            decision="recovery_below_threshold"
        else:
            decision="unresolved"
        rows.append(dict(parameter=name, **{
            f"{label}_{field}":float(src[field])
            for label,src in (("matched_count",a),("joint",b),("stage2_count",c),("legacy_joint",d))
            for field in ("estimate","lower_natural","upper_natural","lower_normalized","upper_normalized")
        }, matched_count_normalized_width=W,joint_normalized_width=width(b),
                         normalized_width_fraction_reduction=narrowed,
                         joint_profile_status=b["profile_status"],recovery_truth_A=A["rate"] if A else "",
                         recovery_truth_B=B["rate"] if B else "",verdict=decision))
    write(HERE/"profile_comparison.tsv",rows)
    return rows


def replay(mode):
    k,c,y=case(mode)
    base=read(HERE/f"best_{mode}.tsv")[0]
    p=c.unmap(row_x(c,base))[0]
    f0=k.evaluate_joint(p,y)[2]
    assert abs(f0-float(base["loss"])) < 1e-9
    profile=[r for r in read(HERE/f"profile_points_{mode}.tsv") if float(r["donor_margin"])>=0.]
    out=[]
    for name in TARGETS:
        center=float(base[f"x_{name}"])
        eligible=[r for r in profile if abs(float(r[f"x_{name}"])-center) > .10000001 and
                  float(r["lr"])<=THRESHOLD]
        if not eligible:
            continue
        r=max(eligible,key=lambda z:abs(float(z[f"x_{name}"])-center))
        p=c.unmap(row_x(c,r))[0]
        metrics,_,f,_=k.evaluate_joint(p,y)
        lr=likelihood_ratio(f,f0,k.n)
        assert metrics[1]>=0. and lr<=THRESHOLD
        assert abs(lr-float(r["lr"]))<1e-9
        out.append(dict(mode=mode,parameter=name,estimate_normalized=center,
                        feasible_coordinate_normalized=float(r[f"x_{name}"]),
                        normalized_displacement=abs(float(r[f"x_{name}"])-center),
                        lr=lr,donor_margin=metrics[1]))
    return out


if __name__=="__main__":
    rows=comparison()
    gate=[r for m in (COUNTS,JOINT,ARCHIVED) for r in replay(m)]
    if gate:
        write(HERE/"out_of_tolerance_witnesses.tsv",gate)
    for r in rows:
        print(r["parameter"],r["verdict"],r["normalized_width_fraction_reduction"],flush=True)
