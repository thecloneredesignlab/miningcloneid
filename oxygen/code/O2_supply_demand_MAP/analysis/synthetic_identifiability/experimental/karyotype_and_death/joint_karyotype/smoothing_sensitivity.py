"""Check registered precision boundaries with one-chromosome observation error."""
from concurrent.futures import as_completed
import numpy as np

from joint_model import HERE, TARGETS, THRESHOLD, likelihood_ratio, pool, read, row_x, write
from refit import case, task_fit

MODE = "autosomes_joint_sd1"


def run():
    k,c,y=case(MODE)
    base=read(HERE/f"best_{MODE}.tsv")[0]
    x=row_x(c,base)
    baseline=base["loss"]
    infile=HERE/"profile_points_autosomes_joint.tsv"
    points=read(infile) if infile.exists() else [
        row for f in HERE.glob("profile_autosomes_joint_*.tsv") for row in read(f)
    ]
    rows=[]
    cases=[]
    for name in TARGETS:
        j=c.names.index(name)
        outs=[]
        for r in points:
            if r["parameter"]!=name or float(r["donor_margin"])<0.:
                continue
            v=row_x(c,r)
            shift=abs(x[j]-v[j])
            if shift<=.1:
                continue
            p=c.unmap(v)[0]
            out,_,value,_=k.evaluate_joint(p,y)
            lr=likelihood_ratio(value,float(baseline),k.n)
            if lr<=THRESHOLD and out[1]>=0.:
                outs.append((shift,lr,float(out[1]),float(v[j])))
        if outs:
            shift,lr,margin,z=max(outs)
            rows.append(dict(parameter=name,check="replayed_support_beyond_tolerance",
                             smoothed_estimate=x[j],comparison_coordinate=z,
                             displacement=shift,lr=lr,donor_margin=margin,stationarity=""))
        else:
            for side in (-1,1):
                q=x[j]+side*.10001
                if 0<=q<=1:
                    for alt in (x,row_x(c,read(HERE/"best_autosomes_joint.tsv")[0])):
                        cases.append((name,j,q,alt))
    with pool(3) as executor:
        jobs={executor.submit(task_fit,(MODE,"precision_boundary",x0,{j:u},72,20.)):
              (name,j,u) for name,j,u,x0 in cases}
        for job in as_completed(jobs):
            r=job.result()
            name,j,u=jobs[job]
            rows.append(dict(parameter=name,check="refitted_tolerance_boundary",
                             smoothed_estimate=x[j],comparison_coordinate=u,
                             displacement=abs(x[j]-u),
                             lr=likelihood_ratio(r["loss"],float(baseline),k.n),
                             donor_margin=r["donor_margin"],stationarity=r["stationarity"]))
    write(HERE/"observation_error_sensitivity.tsv",sorted(rows,key=lambda r:(r["parameter"],r["lr"])))


if __name__ == "__main__":
    run()
