"""Summarize conditional profiles, blinded recovery, and likelihood calibration."""
import argparse
import concurrent.futures as cf
from collections import defaultdict
import math

import numpy as np
from scipy.stats import binomtest

from fit_stage import model
from inference import conf_for, likelihood_ratio, pool, SLACK, THRESHOLD
from practical import HERE, MECHANISMS, Configuration, read, row_x, solve_constrained, write
from recovery import KKT_TOL


def get_bool(v):
    return str(v).lower()=="true"


def group(rows, key):
    d=defaultdict(list)
    for r in rows:
        d[r[key]].append(r)
    return d


def lr_truth(args):
    case, dataset, target, x, loss, truth, parameters = args
    k=model()
    c=Configuration(k)
    result=[]
    for j, name in enumerate(MECHANISMS):
        if parameters and name not in parameters:
            continue
        fixed=x.copy()
        fixed[j]=truth[j]
        p,_=c.unmap(fixed)
        init,_,supply=k.constrained(p,target,c,gradient=False)
        # Exact equality of the estimate and truth at the parameter bounds is common.
        direct_lr=likelihood_ratio(init,loss,k.n) if supply>=-SLACK else np.inf
        start=x
        if direct_lr>THRESHOLD and name not in ("lam_max","buffer_smax"):
            for fraction in (.002,.012,.06,.20,.50,1.):
                z=(1-fraction)*fixed+fraction*truth
                v,_,s=k.constrained(c.unmap(z)[0],target,c,gradient=False)
                if s>=-SLACK and (supply<-SLACK or v<init):
                    init,supply,start=v,s,z
                    direct_lr=likelihood_ratio(init,loss,k.n)
                    if direct_lr<=THRESHOLD:
                        break
        if (abs(x[j]-truth[j]) < 5e-7 or
                (direct_lr<=THRESHOLD and name not in ("lam_max","buffer_smax"))) and supply >= -SLACK:
            f=init
            kval=0.
            okay=True
        else:
            profile=solve_constrained(k,c,target,start,limit=36,fixed=(j,truth[j]))
            f,kval,supply=profile["loss"],profile["gradient"],profile["supply"]
            okay=(kval<=KKT_TOL or profile["scipy_success"]) and supply>=-SLACK
        lr=likelihood_ratio(f,loss,k.n) if supply>=-SLACK else np.inf
        result.append(dict(scenario=case,dataset=dataset,parameter=name,lr_at_true_value=lr,
                           covered_nominal=(lr<=THRESHOLD and supply>=-SLACK),
                           coverage_resolved=lr<=THRESHOLD or okay,
                           constrained_fit_qualified=okay,stationarity=kval,margin=supply))
    return result


def coverage(workers, priority=False):
    inputs=defaultdict(list)
    for r in read(HERE/"noisy_slope_inputs.tsv"):
        inputs[r["scenario"],int(r["dataset"])].append(float(r["ols_slope"]))
    fits=read(HERE/"qualified_recovery.tsv")
    k=model()
    c=Configuration(k)
    truth={r["scenario"]:np.array([float(r[q]) for q in c.names]) for r in read(HERE/"synthetic_truth.tsv")}
    jobs=[]
    for row in fits:
        if not get_bool(row["qualified"]):
            continue
        a,i=row["scenario"],int(row["dataset"])
        if (priority and i<=20) or (not priority and i>20):
            continue
        jobs.append((a,i,np.array(inputs[a,i]),row_x(c,row),float(row["loss"]),truth[a],
                     ("lam_max","buffer_smax") if priority else None))
    rows=[]
    result_path=HERE/("recovery_coverage_priority.tsv" if priority else "recovery_coverage_pilot.tsv")
    if result_path.exists():
        rows=read(result_path)
    done={(r["scenario"],int(r["dataset"])) for r in rows}
    with pool(workers) as w:
        futures=[w.submit(lr_truth,x) for x in jobs if (x[0],x[1]) not in done]
        for future in cf.as_completed(futures):
            result=future.result()
            rows.extend(result)
            print("coverage", result[0]["scenario"],result[0]["dataset"],flush=True)
            write(result_path,sorted(rows,key=lambda r:(r["scenario"],int(r["dataset"]),r["parameter"])))
    return rows


def endpoints(rows, cutoff, center):
    ordered=sorted((float(r["fixed_normalized"]),float(r["lr"]),r)
                   for r in rows if get_bool(r["feasible"]))
    below=[x for x,lr,_ in ordered if lr<=cutoff]
    if not below:
        return "", "", "no_valid_profile_points"
    left,right=min(below),max(below)
    edges=[]
    for direction, endpoint in ((-1,left),(1,right)):
        pairs=[(x,l) for x,l,_ in ordered if (x-endpoint)*direction>0 and l>cutoff]
        if pairs:
            x,lr=min(pairs,key=lambda item:abs(item[0]-endpoint))
            k=next(l for z,l,_ in ordered if z==endpoint)
            out=endpoint+(x-endpoint)*(cutoff-k)/(lr-k)
            state="interpolated_crossing"
        elif (direction==-1 and endpoint<=1e-8) or (direction==1 and endpoint>=1-1e-8):
            out=endpoint
            state="domain_boundary"
        else:
            out=""
            state="outside_endpoint_unresolved"
        edges.append((out,state))
    return edges[0][0],edges[1][0],";".join(e[1] for e in edges)


def summarize():
    k=model()
    c=Configuration(k)
    rfit=read(HERE/"main_constrained_best.tsv")[0]
    x= row_x(c,rfit)
    raw_profiles=group(read(HERE/"main_profiles_envelope.tsv"),"profile_parameter")
    from sensitivity_summary import summarize as summarize_sensitivity
    summarize_sensitivity()
    alt={(r["case"],r["parameter"]):r["sensitivity_status"]
         for r in read(HERE/"sensitivity_summary.tsv")}
    recovery=group(read(HERE/"qualified_recovery.tsv"),"scenario")
    truth={r["scenario"]:np.array([float(r[q]) for q in c.names]) for r in read(HERE/"synthetic_truth.tsv")}
    check=HERE/"recovery_coverage_pilot.tsv"
    withheld=[] if not check.exists() else read(check)
    priority_file=HERE/"recovery_coverage_priority.tsv"
    if priority_file.exists():
        withheld.extend(read(priority_file))
    interval=defaultdict(list)
    for v in withheld:
        interval[v["scenario"],v["parameter"]].append(v)
    rows, rep, real_profiles, calibrated=[],[],[],[]
    for j, name in enumerate(MECHANISMS):
        valid=[r for r in raw_profiles[name] if get_bool(r["feasible"])]
        delta=min(float(r["lr"]) for r in valid)
        if delta < -.10:
            raise RuntimeError(f"profile found better real fit for {name}: LR {delta}")
        low,high,status=endpoints(raw_profiles[name],THRESHOLD,x[j])
        def native(u):
            if u=="":
                return ""
            z=x.copy()
            z[j]=u
            return float(c.unmap(z)[0][j])
        support=[float(r["fixed_normalized"]) for r in valid if float(r["lr"])<=THRESHOLD]
        far=any(abs(a-x[j])>.1+1e-10 for a in support)
        if far:
            profile_grade="outside_precision_confirmed_by_profile_points"
        elif low!="" and high!="" and max(abs(low-x[j]),abs(high-x[j]))<=.10:
            profile_grade="nominal_interval_within_precision"
        else:
            profile_grade="profile_incomplete_or_outside_precision"
        stats=[]
        for a,rr in recovery.items():
            qualified=[z for z in rr if get_bool(z["qualified"])]
            signed=np.array([float(z[f"x_{name}"])-truth[a][j] for z in qualified])
            natural_truth=c.unmap(truth[a])[0][j]
            values=np.array([float(z[name])-natural_truth for z in qualified])
            recovered=sum(abs(signed)<=.1)
            binom=binomtest(int(recovered),len(rr)).proportion_ci(confidence_level=.95, method="wilson")
            supported_values=[v for v in interval[a,name] if get_bool(v["coverage_resolved"])]
            cov=sum(get_bool(v["covered_nominal"]) for v in supported_values)
            stat=dict(scenario=a,parameter=name,assigned=len(rr),qualified=len(signed),
                      unresolved_fits=len(rr)-len(signed),recovered_count=int(recovered),
                      recovery_rate_all=float(recovered/len(rr)),
                      recovery_rate_if_all_unresolved_recovered=(recovered+len(rr)-len(signed))/len(rr),
                      recovery_wilson95_lower=binom.low,recovery_wilson95_upper=binom.high,
                      recovery_rate_among_qualified=(float(recovered/len(signed)) if len(signed) else ""),
                      bias_normalized=(float(signed.mean()) if len(signed) else ""),
                      rmse_normalized=(float(np.sqrt(np.mean(signed*signed))) if len(signed) else ""),
                      p90_abs_error_normalized=(float(np.quantile(abs(signed),.90)) if len(signed) else ""),
                      bias_natural=(float(values.mean()) if len(signed) else ""),
                      p90_abs_error_natural=(float(np.quantile(abs(values),.90)) if len(signed) else ""),
                      coverage_resolved=len(supported_values),covered_nominal95=cov)
            stats.append(stat);rep.append(stat)
        rate=min(s["recovery_rate_all"] for s in stats)
        if profile_grade.startswith("outside"):
            verdict="insufficient_profile_precision"
        elif profile_grade.startswith("nominal") and rate>=.90 and all(s["assigned"]>=100 for s in stats):
            verdict="precision_gate_passed_conditional"
        elif min(s["recovery_rate_if_all_unresolved_recovered"] for s in stats)<.90:
            verdict="recovery_below_precision_gate"
        else:
            verdict="unresolved"
        alt_values={a+"_sensitivity":alt.get((a,name),"not_evaluated") for a in ("passage","noise")}
        if name in ("lam_max","buffer_smax"):
            for scenario in recovery:
                lr_points=[max(0.,float(v["lr_at_true_value"])) for v in interval[scenario,name]
                           if get_bool(v["constrained_fit_qualified"])
                           and np.isfinite(float(v["lr_at_true_value"]))]
                if len(lr_points)>=12:
                    level=float(np.quantile(lr_points,.95,method="higher"))
                    cilo,cihi,cistatus=endpoints(raw_profiles[name],level,x[j])
                    calibrated.append(dict(parameter=name,scenario=scenario,simulations=len(lr_points),
                                           bootstrap_95_lr=level,normalized_lower=cilo,normalized_upper=cihi,
                                           lower=native(cilo),upper=native(cihi),status=cistatus))
        rows.append(dict(parameter=name,estimate=rfit[name],estimate_normalized=x[j],
                         profile95_natural_lower=native(low),profile95_natural_upper=native(high),
                         profile95_normalized_lower=low,profile95_normalized_upper=high,
                         profile_bound_status=status,profile_precision=profile_grade,
                         worst_recovery_rate_all=rate,verdict=verdict,**alt_values))
        for v in valid:
            real_profiles.append({q:v[q] for q in (
                "profile_parameter","fixed_normalized","lr","feasible","max_projected_gradient","termination")})
    write(HERE/"parameter_practical_verdicts.tsv",rows)
    write(HERE/"synthetic_parameter_recovery.tsv",rep)
    write(HERE/"profile_valid_points.tsv",real_profiles)
    if calibrated:
        write(HERE/"empirical_likelihood_calibration.tsv",calibrated)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("command",choices=("coverage","summarize"))
    parser.add_argument("--workers",type=int,default=7)
    parser.add_argument("--priority",action="store_true")
    args=parser.parse_args()
    if args.command=="coverage":
        coverage(args.workers,args.priority)
    else:
        summarize()
