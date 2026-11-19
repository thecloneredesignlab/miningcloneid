"""Replay the fits and combine the agreed profile and recovery precision gates."""
import argparse
from collections import defaultdict
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from ar1_model import (HERE, STATIONARY, TARGETS, THRESHOLD, Config, likelihood_ratio, model,
                       read, row_x, write)
from profile_ar1 import load_profiles, reference
from recover_ar1 import MAX_AGREEMENT_LR, combined, draws


def summarize(require_complete):
    k=model()
    c=Config(k)
    intervals=read(HERE/"profile_intervals.tsv")
    allowed={r["parameter"] for r in intervals if r["profile_status"]=="profile_pass"}
    records=combined()
    outcomes=defaultdict(list)
    for r in records:
        outcomes[r["scenario"]].append(r)
    truth_file=HERE/"truths.tsv"
    truths={a["scenario"]:row_x(c,a) for a in read(truth_file)} if truth_file.exists() else {}
    summary=[]
    score=[]
    for r in intervals:
        p=r["parameter"]
        j=c.names.index(p)
        per_case=[]
        for name in truths if p in allowed else ():
            rows=outcomes[name]
            verified=[r for r in rows if r["qualified"]=="True"]
            errors=np.array([float(v[f"x_{p}"])-truths[name][j] for v in verified])
            hits=int(np.sum(abs(errors)<=.1))
            n=len(rows)
            if n:
                wilson=binomtest(hits,n).proportion_ci(method="wilson")
                per_case.append(dict(parameter=p,scenario=name,assigned=n,numerically_verified=len(verified),
                                     recovered=hits,recovery_fraction=hits/n,
                                     upper_if_all_unverified_recover=(hits+n-len(verified))/n,
                                     wilson95_lower=wilson.low,wilson95_upper=wilson.high,
                                     mean_normalized_error=float(np.mean(errors)) if len(errors) else "",
                                     p90_absolute_normalized_error=float(np.quantile(abs(errors),.90))
                                             if len(errors) else ""))
        score.extend(per_case)
        done=len(per_case)==2 and all(d["assigned"]==100 for d in per_case)
        status=r["profile_status"]
        if status=="outside_supported":
            decision="profile_outside_tolerance"
        elif status=="unresolved":
            decision="profile_unresolved"
        elif not done:
            decision="recovery_pending"
        elif all(z["recovered"]>=90 for z in per_case):
            decision="precision_pass"
        elif any(z["upper_if_all_unverified_recover"]<.90 for z in per_case):
            decision="recovery_shortfall"
        else:
            decision="recovery_inconclusive_due_to_numerics"
        data=dict(r)
        for i,name in enumerate(("A_observed_fit","B_feasible_alternate")):
            stat=next((s for s in per_case if s["scenario"]==name),None)
            data[f"{'AB'[i]}_recovered"]=stat["recovered"] if stat else ""
            data[f"{'AB'[i]}_assigned"]=stat["assigned"] if stat else ""
        summary.append(dict(data,verdict=decision))
    if require_complete and allowed:
        assert len(records)==200 and all(x["verdict"]!="recovery_pending" for x in summary)
        assert set(allowed)=={r["parameter"] for r in read(HERE/"recovery_targets.tsv")}
    if score:
        write(HERE/"recovery_metrics.tsv",score)
    write(HERE/"practical_verdicts.tsv",summary)
    for r in summary:
        print(r["parameter"],r["verdict"],r["A_recovered"],r["B_recovered"],flush=True)
    return summary


def validate():
    for r in read(HERE/"input_hashes.tsv"):
        assert sha256(Path(r["path"]).read_bytes()).hexdigest()==r["sha256"],r["path"]
    k=model()
    c=Config(k)
    final=read(HERE/"practical_verdicts.tsv")
    assert {r["parameter"] for r in final}==set(TARGETS)
    base=reference()
    p,_=c.unmap(row_x(c,base))
    R=k.noise_covariance(p[-1])
    r=k.predict_logs(p)-k.y
    _,ld=np.linalg.slogdet(R)
    q=r@np.linalg.solve(R,r)
    check=.5*q*np.exp(ld/k.n)
    assert abs(check-float(base["loss"]))<1e-9
    assert abs(np.sqrt(q/k.n)-float(base["sigma"]))<1e-12
    assert float(base["donor_margin"])>0
    feasible_profiles=[r for r in load_profiles() if float(r["donor_margin"])>=0]
    maximum_profile_difference=0.
    for row in feasible_profiles:
        own=c.names.index(row["parameter"])
        x=row_x(c,row)
        assert abs(x[own]-float(row["coordinate"]))<1e-9
        p,_=c.unmap(x)
        m,_=k.metrics(p,k.y)
        assert m[1]>=0
        maximum_profile_difference=max(maximum_profile_difference,abs(m[0]-float(row["loss"])))
        assert likelihood_ratio(m[0],check,k.n)>=-.06
    assert maximum_profile_difference<1e-8
    observed_support=defaultdict(list)
    for r in read(HERE/"profile_envelope.tsv"):
        if r["feasible"]=="True" and float(r["lr"])<=THRESHOLD:
            observed_support[r["parameter"]].append(float(r["coordinate"]))
    for r in final:
        center=float(r["normalized_estimate"])
        outside=any(abs(v-center)>.1000000001 for v in observed_support[r["parameter"]])
        assert outside==(r["profile_status"]=="outside_supported")
    print("Validated profile records",len(feasible_profiles),"max replayed loss difference",maximum_profile_difference,flush=True)
    free={}
    for mode in ("profile","recovery"):
        path=HERE/f"compensation_paths_{mode}.tsv"
        if path.exists():
            for s in read(path):
                xx=row_x(c,s)
                xj=c.names.index(s["parameter"])
                assert abs(xx[xj]-float(s["coordinate"]))<1e-10
                m,_=k.metrics(c.unmap(xx)[0],k.y)
                assert m[1]>=0 and abs(m[0]-float(s["loss"]))<1e-8
                assert abs(likelihood_ratio(m[0],check,k.n)-float(s["lr"]))<1e-9
                free[s["parameter"],s["side"]]=s
    tested=verified=0
    for mode in ("profile","supplement","recovery"):
        path=HERE/f"compensation_tests_{mode}.tsv"
        if not path.exists():
            continue
        for line in read(path):
            tested+=1
            s=free[line["parameter"],line["side"]]
            assert abs(float(line["free_lr"])-float(s["lr"]))<1e-9
            if line["usable"]=="True":
                assert float(line["locked_donor_margin"])>=0
                assert float(line["locked_stationarity"])<=STATIONARY
                assert float(line["free_stationarity"])<=STATIONARY
                assert abs(float(line["locked_lr"])-float(line["free_lr"])-float(line["delta_lr"]))<1e-9
                verified+=1
    summary_path=HERE/"compensation_summary.tsv"
    if summary_path.exists():
        blocked={r["parameter"] for r in final if r["verdict"] in (
                    "profile_outside_tolerance","recovery_shortfall")}
        assert blocked=={r["parameter"] for r in read(summary_path)}
    print("Validated compensation",tested,"tests,",verified,"qualifying numeric checks",flush=True)
    success=[r for r in final if r["profile_status"]=="profile_pass"]
    if not success:
        return
    trials=read(HERE/"qualified_recovery.tsv")
    y=draws()
    true={r["scenario"]:row_x(c,r) for r in read(HERE/"truths.tsv")}
    assert len(trials)==200 and len(y)==200
    assert {(r["scenario"],int(r["dataset"])) for r in trials}==set(y)
    truth_rows=read(HERE/"truths.tsv")
    for truth in truth_rows:
        theta,_=c.unmap(row_x(c,truth))
        assert k.supply(theta)[0].min()>0
        assert abs(theta[-1]-float(truth["rho"]))<1e-15
    quality=[]
    for key in true:
        here=[r for r in trials if r["scenario"]==key]
        good=[r for r in here if r["qualified"]=="True"]
        errors=[]
        for t in good:
            p,_=c.unmap(row_x(c,t))
            m,_=k.metrics(p,y[key,int(t["dataset"])])
            assert m[1]>=0 and float(t["stationarity"])<=STATIONARY
            assert float(t["independent_lr_gap"])<=MAX_AGREEMENT_LR
            errors.append(abs(m[0]-float(t["loss"])))
        assert max(errors)<1e-8
        case=key[0]
        for row in success:
            j=c.names.index(row["parameter"])
            hit=sum(abs(row_x(c,z)[j]-true[key][j])<=.1 for z in good)
            assert hit==int(row[f"{case}_recovered"])
        quality.append(dict(scenario=key,assigned=len(here),verified=len(good),
                            largest_recovery_loss_discrepancy=max(errors),
                            replayed_real_profile_fits=len(feasible_profiles),
                            largest_profile_loss_discrepancy=maximum_profile_difference))
    write(HERE/"final_validation.tsv",quality)
    print("Validated recovery",quality,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--validate",action="store_true")
    args=p.parse_args()
    summarize(args.validate)
    if args.validate:
        validate()
