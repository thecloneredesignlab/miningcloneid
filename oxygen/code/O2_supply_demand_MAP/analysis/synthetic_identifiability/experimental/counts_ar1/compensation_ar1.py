"""Confirm compensation on supported profile branches by fixing companion parameters."""
import argparse
from collections import defaultdict
from concurrent.futures import as_completed

import numpy as np

from ar1_model import (HERE, STATIONARY, TARGETS, THRESHOLD, Config, choose, likelihood_ratio,
                       model, pool, read, receipt, row_x, solve, write)
from profile_ar1 import checked_solve, load_profiles, reference

SEEDS={
    "lam_max":("alpha_o2","gamma_growth","mu_hp","buffer_smax"),
    "buffer_smax":("buffer_beta","buffer_n_exp","p_misseg","k_o_mis","p_mis_base"),
    "buffer_beta":("buffer_smax","buffer_n_exp","p_misseg","p_mis_base","k_o_mis"),
    "buffer_n_exp":("buffer_beta","buffer_smax","p_misseg","k_o_mis"),
    "n_O":("O2_crit","alpha_o2","mu_hp","gamma_growth","gamma_mu"),
    "p_misseg":("k_o_mis","p_mis_base","mu_hp","buffer_smax","buffer_beta"),
    "k_o_mis":("p_misseg","p_mis_base","mu_hp","gamma_mu","buffer_beta"),
    "gamma_growth":("alpha_o2","mu_hp","gamma_mu","O2_crit","n_O"),
    "mu_hp":("gamma_mu","k_o_mis","p_misseg","eta_H_2N","O2_crit"),
    "O2_crit":("n_O","alpha_o2","mu_hp","gamma_growth","gamma_mu"),
}


def select_cases(parameters):
    k=model()
    c=Config(k)
    best=reference()
    x=row_x(c,best)
    base=float(best["loss"])
    candidates=[]
    for r in load_profiles()+read(HERE/"observed_refits.tsv"):
        if float(r["donor_margin"])<0:
            continue
        lr=likelihood_ratio(float(r["loss"]),base,k.n)
        if lr<=THRESHOLD:
            candidates.append((row_x(c,r),lr,r.get("parameter",r.get("origin"))))
    selected=[]
    for p in parameters:
        j=c.names.index(p)
        for direction,label in ((-1,"lower"),(1,"upper")):
            points=[(direction*(q[j]-x[j]),q,lr,origin) for q,lr,origin in candidates
                    if direction*(q[j]-x[j])>.008]
            if not points:
                continue
            biggest=max(s[0] for s in points)
            # For broad profiles test an interior part of the long ridge.
            # A separate directional branch is kept whenever it is supported.
            offset=min(biggest,max(.115,.65*biggest)) if biggest>.1 else biggest
            if biggest>.1:
                points=[r for r in points if r[0]>.10000001]
            d,q,lr,origin=min(points,key=lambda r:(abs(r[0]-offset),r[2]))
            selected.append((p,label,x,q,base,lr,origin))
    return selected


def unrestricted(args):
    name,side,x,q,baseline,old_lr,origin=args
    k=model(); c=Config(k)
    j=c.names.index(name)
    f=checked_solve(k,c,k.y,{j:q[j]},[q,x],limits=(80,75),divisors=(18.,6.))
    if float(f["donor_margin"])<0 or likelihood_ratio(float(f["loss"]),baseline,k.n)>old_lr+.00001:
        f=receipt(k,c,k.y,dict(x=q,gradient=np.nan,calls=0,iterations=0,
              success=False,termination="retained feasible point",rounds=0,seconds=0))
    return dict(parameter=name,side=side,coordinate=q[j],starting_source=origin,
                lr=likelihood_ratio(float(f["loss"]),baseline,k.n),**f)


def fixed_other(args):
    free,ref,others,*alternates=args
    k=model(); c=Config(k)
    x=row_x(c,ref)
    z=row_x(c,free)
    p=free["parameter"]
    j=c.names.index(p)
    i=[c.names.index(h) for h in others]
    coords={j:float(free["coordinate"])} | {h:x[h] for h in i}
    if alternates:
        t=checked_solve(k,c,k.y,coords,alternates[0],limits=(112,108,105,119),divisors=(18.,4.,40.,9.))
    else:
        t=checked_solve(k,c,k.y,coords,[z,x],limits=(75,78),divisors=(12.,5.))
    fref=float(ref["loss"])
    fit_lr=likelihood_ratio(float(t["loss"]),fref,k.n)
    difference=likelihood_ratio(float(t["loss"]),float(free["loss"]),k.n)
    feasible=float(t["donor_margin"])>=0
    valid=feasible and float(t["stationarity"])<=STATIONARY and float(free["stationarity"])<=STATIONARY
    return dict(parameter=p,side=free["side"],fixed_target_coordinate=free["coordinate"],
                companions=",".join(others),num_companions=len(others),
                signed_companion_move=",".join(f"{z[h]-x[h]:+.6g}" for h in i),
                free_lr=float(free["lr"]),locked_lr=fit_lr if feasible else "",
                delta_lr=difference if feasible else "",free_stationarity=free["stationarity"],
                locked_stationarity=t["stationarity"],locked_donor_margin=t["donor_margin"],
                usable=valid and difference>=-.06,
                crossing_95=valid and float(free["lr"])<=THRESHOLD<fit_lr)


def run(workers,mode):
    current=reference()
    k=model()
    c=Config(k)
    x=row_x(c,current)
    if mode=="profile":
        profiles=read(HERE/"profile_intervals.tsv")
        parameters=[r["parameter"] for r in profiles if r["profile_status"]!="profile_pass"]
    else:
        parameters=[r["parameter"] for r in read(HERE/"practical_verdicts.tsv") if r["verdict"]=="recovery_shortfall"]
    if not parameters:
        return
    path=HERE/f"compensation_paths_{mode}.tsv"
    jobs=select_cases(parameters)
    if not jobs:
        return
    paths=[]
    with pool(workers) as pool_obj:
        for f in as_completed([pool_obj.submit(unrestricted,t) for t in jobs]):
            r=f.result()
            paths.append(r)
            print("supported displacement",r["parameter"],r["side"],r["lr"],flush=True)
            write(path,paths)
    tasks=[]
    for r in paths:
        p=r["parameter"]
        z=row_x(c,r)
        dif=abs(z-x)
        priority=[c.names[j] for j in np.argsort(-dif) if c.names[j]!=p]
        selected=priority[:3]
        mechanistic=sorted(SEEDS[p],key=lambda h:-dif[c.names.index(h)])
        for h in mechanistic:
            if h not in selected and dif[c.names.index(h)]>.005:
                selected.append(h)
            if len(selected)>=5:
                break
        tasks.extend((r,current,[h]) for h in selected)
        movement=[h for h in selected if dif[c.names.index(h)]>.02][:3]
        if len(movement)>=2:
            tasks.append((r,current,movement))
    out=[]
    with pool(workers) as pool_obj:
        for f in as_completed([pool_obj.submit(fixed_other,t) for t in tasks]):
            r=f.result()
            out.append(r)
            print("companion",r["parameter"],r["side"],r["companions"],r["delta_lr"],
                  r["usable"],flush=True)
            write(HERE/f"compensation_tests_{mode}.tsv",out)


def summarize():
    rows=[]
    p=HERE/"compensation_tests_profile.tsv"
    r=HERE/"compensation_tests_recovery.tsv"
    for path in (p,r,HERE/"compensation_tests_supplement.tsv"):
        if path.exists():
            rows+=read(path)
    if not rows:
        return
    final=[]
    for target in TARGETS:
        here=[r for r in rows if r["parameter"]==target]
        if not here:
            continue
        stable=[r for r in here if r["usable"]=="True"]
        separate=[r for r in stable if r["num_companions"]=="1"]
        by_companion=defaultdict(list)
        for r in separate:
            by_companion[r["companions"]].append(r)
        best=[]
        for name,v in by_companion.items():
            best.append(max(v,key=lambda x:float(x["delta_lr"])))
        best.sort(key=lambda r:float(r["delta_lr"]),reverse=True)
        main=[r for r in best if float(r["delta_lr"])>=.35][:3]
        combos=sorted((r for r in stable if r["num_companions"]!="1"),
                      key=lambda r:float(r["delta_lr"]),reverse=True)
        joint=combos[0] if combos and float(combos[0]["delta_lr"])>=.35 else None
        final.append(dict(parameter=target,n_verified_single_tests=len(separate),n_trial_total=len(here),
                          leading_companions=";".join(r["companions"] for r in main),
                          additional_lr=";".join(f'{float(r["delta_lr"]):.3f}' for r in main),
                          direction=";".join(r["side"] for r in main),
                          locked_above_nominal_95=";".join(r["crossing_95"] for r in main),
                          best_joint_set=joint["companions"] if joint else "",
                          best_joint_additional_lr=joint["delta_lr"] if joint else "",
                          best_joint_crossing_95=joint["crossing_95"] if joint else ""))
    write(HERE/"compensation_summary.tsv",final)
    for r in final:
        print(r["parameter"],r["leading_companions"],r["additional_lr"],r["best_joint_set"],flush=True)


def supplement(workers):
    """Retest meaningful numerical failures and the explicit model couplings."""
    k=model()
    c=Config(k)
    fit=reference()
    x=row_x(c,fit)
    points=read(HERE/"compensation_paths_profile.tsv")
    seen={(r["parameter"],r["side"],r["companions"]):r
           for r in read(HERE/"compensation_tests_profile.tsv")}
    bank=[row_x(c,r) for r in load_profiles()+read(HERE/"observed_refits.tsv")
          if float(r["donor_margin"])>=0]
    jobs=[]
    for free in points:
        name=free["parameter"]
        side=free["side"]
        j=c.names.index(name)
        z=row_x(c,free)
        options=[]
        for partner in SEEDS[name][:3]:
            prior=seen.get((name,side,partner))
            if prior is None and abs(z[c.names.index(partner)]-x[c.names.index(partner)])>.012:
                options.append(partner)
        failures=[r for (p,s,_),r in seen.items() if p==name and s==side and
                  r["num_companions"]=="1" and r["usable"]=="False" and
                  r["delta_lr"] and float(r["delta_lr"])>.5]
        options+= [r["companions"] for r in sorted(failures,key=lambda q:float(q["delta_lr"]),reverse=True)]
        for name2 in dict.fromkeys(options):
            i=c.names.index(name2)
            nearest=sorted(bank,key=lambda s:abs(s[j]-z[j])+abs(s[i]-x[i]))[:25]
            evaluated=[]
            for r in (z,x,*nearest):
                zz=r.copy()
                zz[j]=z[j]
                zz[i]=x[i]
                m,_=k.metrics(c.unmap(zz)[0],k.y)
                evaluated.append((m[1]<0,m[0] if m[1]>=0 else -m[1],zz))
            evaluated.sort(key=lambda s:s[:2])
            starts=[]
            for _,_,r in evaluated:
                if all(np.linalg.norm(r-u)>.008 for u in starts):
                    starts.append(r)
                if len(starts)>=4:
                    break
            jobs.append((free,fit,[name2],starts or [z,x]))
    if jobs:
        rows=[]
        with pool(workers) as work:
            for f in as_completed([work.submit(fixed_other,j) for j in jobs]):
                r=f.result()
                rows.append(r)
                write(HERE/"compensation_tests_supplement.tsv",rows)
                print("companion retry",r["parameter"],r["side"],r["companions"],
                      r["delta_lr"],r["usable"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=("profile","recovery","summarize","supplement"))
    p.add_argument("--workers",type=int,default=5)
    args=p.parse_args()
    if args.mode=="supplement":
        supplement(args.workers)
    elif args.mode!="summarize":
        run(args.workers,args.mode)
    summarize()
