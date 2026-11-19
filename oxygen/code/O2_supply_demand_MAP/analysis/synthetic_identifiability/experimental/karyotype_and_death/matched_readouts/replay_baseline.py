"""Reuse the stage-three baseline profiles after evaluating graft equivalence."""
import numpy as np

from modal import HERE,STAGE3,TARGETS,case,read,write


def run():
    mode="stock_baseline"
    k,c,y=case(mode)
    ref=read(HERE/f"best_{mode}.tsv")[0]
    data=[]
    for r in read(STAGE3/"profile_points_autosomes_joint.tsv"):
        if r["parameter"] not in TARGETS:
            continue
        p,_=c.unmap(c.anchor(r))
        ar,f,_,ky,de=k.parts(p,y)
        assert abs(f-k.n*np.log(float(r["loss"])))<1e-8
        assert abs(ar[1]-float(r["donor_margin"]))<1e-9
        r.update(mode=mode,loss=f,count_loss=ar[0],kary_nll=ky,death_nll=de,
                 lr=f-float(ref["loss"]),donor_margin=ar[1],
                 stationarity=float(r["stationarity"])*k.n/float(r["loss"]))
        data.append(r)
    write(HERE/f"profile_points_{mode}.tsv",data)
    summary=[]
    for r in read(STAGE3/"profile_intervals_autosomes_joint.tsv"):
        if r["parameter"] not in TARGETS:
            continue
        summary.append(dict(mode=mode,parameter=r["parameter"],estimate=r["estimate"],
            x_estimate=r["normalized_estimate"],lower=r["lower_natural"],upper=r["upper_natural"],
            x_lower=r["lower_normalized"],x_upper=r["upper_normalized"],
            x_min_supported=r["lower_supported"],x_max_supported=r["upper_supported"],
            lower_type=r["lower_type"],upper_type=r["upper_type"],precision_gate=r["profile_status"]))
    write(HERE/f"intervals_{mode}.tsv",summary)
    print(len(data),"stage-three profile points replayed unchanged",flush=True)


if __name__=="__main__":
    run()
