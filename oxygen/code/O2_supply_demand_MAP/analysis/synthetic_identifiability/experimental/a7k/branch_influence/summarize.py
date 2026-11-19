"""Replay subset fits, then summarize diagnostics on matched simulations."""
import argparse
from collections import defaultdict

import numpy as np

from branch_influence import (ARCH, MODES as CAPS, NAMES, OUT, THRESHOLD, archived_rows, case,
                              deviance_rows, masked, read, stable, tag, write)
from recovery_diagnostics import MODES, ORACLES, archived, paired_stats, results, score


def subsets():
    table=[]
    components=[]
    for mode in CAPS:
        k,c,y=case(mode)
        base,full,points=archived_rows(mode)
        rp=[]
        for mask in range(16):
            name=tag(mode,mask)
            model,cfg,target,_=masked(mode,mask)
            if mask in (0,15):
                free=base if mask==0 else full
                limit=base if mask==0 else points[1.]
                d=model.parts(cfg.unmap(cfg.anchor(free))[0],target)[1]
                f=model.parts(cfg.unmap(cfg.anchor(limit))[0],target)[1]
                assert abs(d-float(free["loss"]))<1e-7 and abs(f-float(limit["loss"]))<1e-7
                if mask==0:
                    s=read(ARCH.with_name("20261119_experimental_death_A7K_comparison_01")/
                           "intervals_stock_baseline.tsv")
                else:
                    s=read(ARCH/f"intervals_{mode}.tsv")
                s=next(r for r in s if r["parameter"]=="buffer_beta")
                row=dict(mode=mode,mask=mask,branches=",".join(n for j,n in enumerate(NAMES) if mask&1<<j),
                         n_kary_cells=int(model.kary[2:].sum()),beta=float(free["buffer_beta"]),
                         lower=s["lower"],upper=s["upper"],
                         lower_type=s["lower_type"],upper_type=s["upper_type"],
                         x_lower=s["x_lower"],x_upper=s["x_upper"],deviance=d,beta10_deviance=f,
                         beta10_lr=f-d,beta10_supported=f-d<=THRESHOLD,
                         free_stationarity=free["stationarity"],beta10_stationarity=limit["stationarity"])
            else:
                free,limit=read(OUT/f"profile_ref_{name}.tsv")
                row=read(OUT/f"interval_{name}.tsv")[0]
            f0,f1=0.,0.
            for candidate,which in [(free,"free"),(limit,"beta_10")]:
                p,_=cfg.unmap(cfg.anchor(candidate))
                ar,loss,*_=model.parts(p,target)
                assert ar[1]>=-5e-10 and abs(loss-float(candidate["loss"]))<1e-7
                if model.cap:
                    assert max(model.composition(p))<=model.cap+1e-10
                if which=="free":
                    f0=loss
                else:
                    f1=loss
            assert abs(f1-f0-float(row["beta10_lr"]))<1e-7
            ps=[c.unmap(c.anchor(fit))[0] for fit in (free,limit)]
            d0,d1=[deviance_rows(k,p,y)[1] for p in ps]
            common=[k.parts(p,y)[4] for p in ps]
            lrs=[common[1]-common[0],*[d1[j+2]-d0[j+2] for j in range(4)]]
            assert abs(lrs[0]+sum(v for j,v in enumerate(lrs[1:]) if mask&1<<j)-(f1-f0))<1e-7
            for j,label in enumerate(NAMES):
                count=int(k.kary[j+2].sum())
                components.append(dict(mode=mode,mask=mask,branch=label,retained=bool(mask&1<<j),cells=count,
                   at_free=d0[j+2],at_beta10=d1[j+2],difference=lrs[j+1],
                   per_cell=lrs[j+1]/count,common_component=lrs[0]))
            rp.append(row)
        for item in rp:
            item["support_delta_vs_full"]=float(item["beta10_lr"])-float(rp[-1]["beta10_lr"])
        table.extend(rp)
    write(OUT/"branch_profiles.tsv",table)
    write(OUT/"branch_likelihood_components.tsv",components)
    for row in table:
        print(row["mode"],row["mask"],row["branches"],row["beta"],row["lower"],row["upper"],
              row["beta10_lr"],row["beta10_supported"],row["lower_type"],row["upper_type"],flush=True)


def oracle():
    data=read(results("oracle"))
    d={(r["mode"],int(r["case_id"]),int(r["replicate"]),r["scenario"],r["condition"]):r for r in data}
    panel=read(OUT/"diagnostic_panel.tsv")
    rows=[]
    for m in MODES:
        hist=archived(m)
        for j in (1,2):
            seeds=[int(r["replicate"]) for r in panel if r["mode"]==m and int(r["case_id"])==j]
            for condition in ("free",*ORACLES):
                def error(i,s):
                    if condition=="free":
                        return score(m,j,hist[j,i,s,"relaxed"])
                    r=d[m,j,i,s,condition]
                    return float(r["delta_x_p_misseg"]) if r["validated"]=="True" else np.nan
                qualified=[i for i in seeds if np.isfinite([error(i,s) for s in ("strict","up","down")]).all()]
                for first in ("strict","up"):
                    a,b=([error(i,s) for i in qualified] for s in (first,"down"))
                    if not qualified:
                        continue
                    rows.append(dict(mode=m,case_id=j,condition=condition,first=first,selected=len(seeds),
                                     mean_signed_direction_gap=float(np.mean(np.array(b)-a)),
                                     **paired_stats(a,b,20261119+j)))
    write(OUT/"oracle_directional_comparisons.tsv",rows)


def shape():
    data=read(results("shape"))
    fit=read(results("oracle"))
    index={(r["mode"],r["case_id"],r["replicate"],r["scenario"],r["condition"]):r for r in fit}
    groups=defaultdict(dict)
    for r in data:
        key=tuple(r[field] for field in ("mode","case_id","replicate","scenario","condition"))
        groups[key][float(r["offset"])]=r
    output=[]
    for key,points in groups.items():
        assert sorted(points)==[-.1,-.05,0.,.05,.1]
        m,j,i,s,condition=key
        f=index[key] if condition!="free" else archived(m)[int(j),int(i),s,"relaxed"]
        values={q:float(r["loss"])-float(f["loss"]) for q,r in points.items()}
        if f["validated"]!="True" or any(r["validated"]!="True" for r in points.values()):
            continue
        output.append(dict(mode=m,case_id=j,replicate=i,scenario=s,condition=condition,
                           x_error=score(m,int(j),f),truth_lr=values[0.],
                           lr_at_minus_tol=values[-.1],lr_at_plus_tol=values[.1],
                           truth_curvature=(values[-.05]+values[.05]-2*values[0.])/.05**2,
                           min_relative_grid_lr=min(values.values())))
    write(OUT/"profile_local_shape.tsv",output)


def swaps():
    data=[]
    for phase in ("swap_pilot","swap_full"):
        if results(phase).exists():
            data.extend(read(results(phase)))
    rows=[]
    for branch in sorted({r["branch"] for r in data}):
        for s in ("up","down"):
            for j in (1,2):
                chosen=sorted((r for r in data if r["branch"]==branch and r["scenario"]==s and
                               int(r["case_id"])==j),key=lambda r:int(r["replicate"]))
                n=len(chosen)
                if not n:
                    continue
                assert [int(r["replicate"]) for r in chosen]==list(range(1,n+1)) and n in (20,100)
                m="group_050"
                hist=archived(m)
                a=[score(m,j,hist[j,int(r["replicate"]),s,"relaxed"]) for r in chosen]
                b=[float(r["delta_x_p_misseg"]) if r["validated"]=="True" else -np.inf for r in chosen]
                qualified=sum(np.isfinite(b))
                stat=paired_stats(a,b,20261119+j)
                if qualified<n:
                    stat.update(bias_second=np.nan,sd_second=np.nan,mae_second=np.nan)
                rows.append(dict(branch=branch,background=s,case_id=j,numerically_valid=qualified,
                                 **stat))
    write(OUT/"branch_swap_paired.tsv",rows)
    for r in rows:
        print(r["branch"],r["background"],r["case_id"],r["n"],r["numerically_valid"],
              r["recovered_first"],r["recovered_second"],r["first_only"],r["second_only"],
              flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase",choices=("subsets","oracle","shape","swaps"))
    globals()[p.parse_args().phase]()
