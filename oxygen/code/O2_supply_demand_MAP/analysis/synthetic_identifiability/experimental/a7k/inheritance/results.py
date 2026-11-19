"""Report beta endpoint checks, model fit tradeoffs and the recovery selection."""
import argparse
from pathlib import Path

import numpy as np

from tilt_model import HERE,PREV,GRID,TARGETS,GROUP,BRANCH,THRESHOLD,case,cap,read,stable,write


def endpoint(mode):
    k,c,y=case(mode)
    reference=read(HERE/f"best_{mode}.tsv")[0]
    rs=read(HERE/f"profile_points_{mode}.tsv") if mode!="strict" else read(PREV/"profile_points_stock_a7k.tsv")
    eligible=[r for r in rs if r["parameter"]=="buffer_beta" and float(r["coordinate"])==1.
              and float(r.get("all_constraint_margin",r["donor_margin"]))>=0.]
    assert eligible,mode
    r=min(eligible,key=lambda q:float(q["loss"]))
    p,_=c.unmap(c.anchor(r))
    ar,f,g,ky,b=k.parts(p,y)
    assert abs(f-float(r["loss"]))<1e-7
    assert ar[1]>=-1e-10 and (mode=="strict" or max(k.composition(p)) <= k.cap+1e-12)
    lr=f-float(reference["loss"])
    return lr, float(r["stationarity"])


def table(modes):
    output=[]
    shifts=[]
    base=read(PREV/"best_stock_baseline.tsv")[0]
    base_loss=float(base["loss"])
    for mode in modes:
        k,c,y=case(mode)
        r=read(HERE/f"best_{mode}.tsv")[0]
        p,_=c.unmap(c.anchor(r))
        tv,pre,post,seed=k.composition(p,True)
        ar,f,g,_,b=k.parts(p,y)
        assert abs(f-float(r["loss"]))<1e-8
        interval={t["parameter"]:t for t in read(HERE/f"intervals_{mode}.tsv")}
        lr,residual=endpoint(mode)
        assert (not k.cap) or max(tv)<=k.cap+1e-12
        assert max(abs(np.array([float(a["fixed_recorded_seed"]) for a in k.a7])-seed))<1e-7
        output.append(dict(mode=mode,max_tv=cap(mode),max_realized_tv=max(tv),
                           loss=f,a7k_deviance=f-b,common_deviance=b,
                           common_deviance_excess=b-base_loss,
                           beta_at_10_lr=lr,beta_at_10_stationarity=residual,
                           beta_at_10_in_nominal_support=lr<=THRESHOLD,
                           **{f"{param}_{field}":interval[param][field]
                               for param in TARGETS for field in ("estimate","lower","upper","precision_gate")}))
        for i,branch in enumerate(c.kary_labels):
            shifts.append(dict(mode=mode,branch=branch,total_variation=tv[i],
                               donor_mean=pre[i]@GRID,restart_mean=post[i]@GRID,
                               donor_sd=np.sqrt(pre[i]@(GRID-pre[i]@GRID)**2),
                               restart_sd=np.sqrt(post[i]@(GRID-post[i]@GRID)**2),
                               fixed_seed=seed[i]))
    return output,shifts


def select(results):
    primary=[r for r in results if r["mode"] in GROUP]
    allowed=[r["mode"] for r in primary if r["p_misseg_precision_gate"]=="profile_pass"]
    candidates=[]
    if "group_025" in allowed:
        candidates.append("group_025")
    transition=[r["mode"] for r in primary if r["beta_at_10_in_nominal_support"]]
    target=transition[0] if transition else "group_050"
    if target in allowed and target not in candidates:
        candidates.append(target)
    if len(candidates)<2:
        for name in reversed(allowed):
            if name not in candidates:
                candidates.append(name)
            if len(candidates)==2:
                break
    if candidates:
        write(HERE/"recovery_selection.tsv",[dict(mode=mode,p_prec_gate="profile_pass",
            beta_transition=(mode==target and bool(transition))) for mode in candidates])
    print("Recovery selection",*candidates,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--primary",action="store_true")
    a=p.parse_args()
    modes=("strict",)+GROUP+(BRANCH if not a.primary else ())
    results,shifts=table(modes)
    write(HERE/"profile_comparison.tsv",results)
    write(HERE/"a7k_composition_changes.tsv",shifts)
    select(results)
    for row in results:
        print(row,flush=True)
