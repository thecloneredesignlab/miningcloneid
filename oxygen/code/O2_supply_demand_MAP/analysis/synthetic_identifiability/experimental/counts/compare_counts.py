"""Compare precision decisions and paired estimation errors against the frozen OLS run."""
import argparse
from collections import defaultdict
from hashlib import sha256
from pathlib import Path

import numpy as np

from count_model import HERE, OLS, model
from inference import THRESHOLD, likelihood_ratio, SLACK
from paired_counts import combine, loads
from practical import Configuration, MECHANISMS, read, row_x, write
from profile_counts import load_direct
from recovery import CONVERGE_LR, KKT_TOL


def width(row,lo,hi):
    if row[lo]=="" or row[hi]=="":
        return None
    return float(row[hi])-float(row[lo])


def quantile_interval(values,seed):
    if len(values)<2:
        return ("","")
    rng=np.random.default_rng(seed)
    groups=np.array(values)
    resampled=groups[rng.integers(0,len(groups),size=(4000,len(groups)))]
    return tuple(np.quantile(resampled.mean(axis=1),(.025,.975)))


def summarize():
    new=combine()
    if not new:
        return
    raw_profile={r["parameter"]:r for r in read(HERE/"profile_intervals.tsv")}
    ols_profile={r["parameter"]:r for r in read(OLS/"parameter_practical_verdicts.tsv")}
    old={(r["scenario"],int(r["dataset"])):r for r in read(OLS/"qualified_recovery.tsv")}
    k=model()
    c=Configuration(k)
    truth={r["scenario"]:np.array([float(r[p]) for p in c.names]) for r in read(OLS/"synthetic_truth.tsv")}
    cases=defaultdict(list)
    for r in new:
        cases[r["scenario"]].append(r)
    paired=[]
    rates=[]
    summaries=[]
    for j,p in enumerate(MECHANISMS):
        raw=raw_profile[p]
        slope=ols_profile[p]
        rw=width(raw,"lower_normalized","upper_normalized")
        ow=width(slope,"profile95_normalized_lower","profile95_normalized_upper")
        groups=[]
        for a in truth:
            if a not in cases:
                continue
            x=truth[a][j]
            raw_errors=[]
            ols_errors=[]
            paired_errors=[]
            delta_success=[]
            for r in cases[a]:
                i=int(r["dataset"])
                o=old[a,i]
                count_eligible=r["qualified"]=="True"
                ols_eligible=o["qualified"]=="True"
                rawerr=float(r[f"x_{p}"])-x
                olserr=float(o[f"x_{p}"])-x
                raw_hit=count_eligible and abs(rawerr)<=.10
                ols_hit=ols_eligible and abs(olserr)<=.10
                delta_success.append(int(raw_hit)-int(ols_hit))
                if count_eligible:
                    raw_errors.append(rawerr)
                if ols_eligible:
                    ols_errors.append(olserr)
                if count_eligible and ols_eligible:
                    paired_errors.append(abs(rawerr)-abs(olserr))
                paired.append(dict(parameter=p,scenario=a,dataset=i,x_truth=x,
                                   ols_qualified=ols_eligible,counts_qualified=count_eligible,
                                   ols_error=olserr if ols_eligible else "",
                                   counts_error=rawerr if count_eligible else "",
                                   counts_minus_ols_absolute_error=paired_errors[-1] if count_eligible and ols_eligible else "",
                                   ols_recovered=ols_hit,counts_recovered=raw_hit))
            n=len(delta_success)
            rc=int(np.sum(np.abs(raw_errors)<=.10))
            oc=int(np.sum(np.abs(ols_errors)<=.10))
            rate_lo,rate_hi=quantile_interval(delta_success,33212+j)
            error_lo,error_hi=quantile_interval(paired_errors,52576+j)
            r=dict(parameter=p,scenario=a,assigned=n,counts_qualified=len(raw_errors),
                   ols_qualified=len(ols_errors),paired_qualified=len(paired_errors),
                   counts_recovered=rc,ols_recovered=oc,counts_recovery=rc/n,ols_recovery=oc/n,
                   recovery_change=(rc-oc)/n,recovery_change_95lo=rate_lo,recovery_change_95hi=rate_hi,
                   counts_p90_abs_error=float(np.quantile(np.abs(raw_errors),.9)) if raw_errors else "",
                   ols_p90_abs_error=float(np.quantile(np.abs(ols_errors),.9)) if ols_errors else "",
                   paired_mean_abs_error_change=float(np.mean(paired_errors)) if paired_errors else "",
                   paired_error_change_95lo=error_lo,paired_error_change_95hi=error_hi)
            groups.append(r)
            rates.append(r)
        all_assigned=len(groups)==2 and all(r["assigned"]==100 for r in groups)
        raw_pass=all_assigned and raw["profile_precision"]=="narrow_checked" and all(r["counts_recovery"]>=.9 for r in groups)
        ols_pass=slope["verdict"]=="precision_gate_passed_conditional"
        normalized_profile_narrower=rw is not None and ow is not None and rw<ow
        pair_error_both=all_assigned and all(r["paired_mean_abs_error_change"]!=""
                                            and r["paired_mean_abs_error_change"]<0 for r in groups)
        if raw_pass:
            verdict="both_pass" if ols_pass else "new_pass_counts"
        elif not all_assigned:
            verdict="simulation_incomplete"
        elif normalized_profile_narrower and pair_error_both:
            verdict="improved_both_metrics_without_pass"
        elif normalized_profile_narrower:
            verdict="narrower_profile_without_pass"
        elif pair_error_both:
            verdict="smaller_paired_error_without_pass"
        else:
            verdict="no_joint_precision_gain"
        summaries.append(dict(parameter=p,counts_estimate=raw["estimate"],ols_estimate=slope["estimate"],
                              counts_profile_lower=raw["lower_natural"],counts_profile_upper=raw["upper_natural"],
                              ols_profile_lower=slope["profile95_natural_lower"],ols_profile_upper=slope["profile95_natural_upper"],
                              counts_profile_normalized_width="" if rw is None else rw,
                              ols_profile_normalized_width="" if ow is None else ow,
                              difference_in_normalized_width="" if rw is None or ow is None else rw-ow,
                              counts_profile_precision=raw["profile_precision"],ols_verdict=slope["verdict"],
                              counts_pass=raw_pass,ols_pass=ols_pass,comparative_verdict=verdict))
    write(HERE/"paired_parameter_recovery.tsv",rates)
    write(HERE/"paired_fit_errors.tsv",paired)
    write(HERE/"parameter_comparison.tsv",summaries)
    print("Comparison",[(p["parameter"],p["comparative_verdict"]) for p in summaries],flush=True)


def validate():
    for row in read(HERE/"source_provenance.tsv"):
        assert sha256(Path(row["file"]).read_bytes()).hexdigest()==row["sha256"],row["file"]
    k=model()
    c=Configuration(k)
    targets=loads()
    check=read(HERE/"paired_replay.tsv")
    assert len(check)==200 and all(float(r["ols_slope_max_abs_discrepancy"])<=1e-12 for r in check)
    fits=read(HERE/"qualified_count_recovery.tsv")
    assert len(fits)==len(targets)==200
    profile=read(HERE/"profile_intervals.tsv")
    assert {r["parameter"] for r in profile}==set(MECHANISMS)
    final=read(HERE/"parameter_comparison.tsv")
    pairs=read(HERE/"paired_fit_errors.tsv")
    results=read(HERE/"paired_parameter_recovery.tsv")
    assert len(final)==14 and len(pairs)==2800 and len(results)==28
    real=read(HERE/"real_best.tsv")[0]
    real_x=row_x(c,real)
    pred,_=k.log_predict(c.unmap(real_x)[0])
    independently_replayed_loss=.5*np.dot(pred[k.capture]-k.y,pred[k.capture]-k.y)
    assert abs(independently_replayed_loss-float(real["loss"]))<1e-9
    finite_curves=[r for r in load_direct() if r["feasible"]=="True"]
    profile_differences=[]
    for r in finite_curves:
        x=row_x(c,r)
        assert abs(x[MECHANISMS.index(r["parameter"])]-float(r["coordinate"]))<=1e-10
        f,_,margin=k.constrained(c.unmap(x)[0],k.y,c,False)
        assert margin>=-SLACK
        assert likelihood_ratio(f,float(real["loss"]),k.n)>=-.10
        profile_differences.append(abs(f-float(r["loss"])))
    assert max(profile_differences)<1e-8
    curve=read(HERE/"profiles.tsv")
    truth={r["scenario"]:np.array([float(r[q]) for q in c.names]) for r in read(OLS/"synthetic_truth.tsv")}
    saved_ols=read(OLS/"qualified_recovery.tsv")
    proof=[]
    for case in sorted({r["scenario"] for r in fits}):
        eligible=[r for r in fits if r["scenario"]==case and r["qualified"]=="True"]
        differences=[]
        for r in eligible:
            x=row_x(c,r)
            y=targets[case,int(r["dataset"])]
            loss,_,margin=k.constrained(c.unmap(x)[0],y,c,False)
            differences.append(abs(loss-float(r["loss"])))
            assert margin>=-SLACK
            assert float(r["max_projected_gradient"])<=KKT_TOL
            assert float(r["independent_LR_gap"])<=CONVERGE_LR
        assert max(differences)<1e-8
        for p in MECHANISMS:
            rows=[r for r in pairs if r["scenario"]==case and r["parameter"]==p]
            stat=next(r for r in results if r["scenario"]==case and r["parameter"]==p)
            assert len(rows)==100
            j=MECHANISMS.index(p)
            for mode,members in (("counts",fits),("ols",saved_ols)):
                independent=sum(r["scenario"]==case and r["qualified"]=="True"
                                and abs(float(r[f"x_{p}"])-truth[case][j])<=.1 for r in members)
                assert sum(r[f"{mode}_recovered"]=="True" for r in rows)==int(stat[f"{mode}_recovered"])==independent
        proof.append(dict(scenario=case,n_pairs=100,counts_qualified=len(eligible),
                          largest_replayed_loss_error=max(differences),
                          valid_profile_points_replayed=len(finite_curves),
                          largest_profile_loss_error=max(profile_differences)))
        print(proof[-1],flush=True)
    for r in final:
        p=r["parameter"]
        x=next(x for x in profile if x["parameter"]==p)
        outside=any(z["parameter"]==p and float(z["likelihood_ratio"])<=THRESHOLD
                    and abs(float(z["coordinate"])-float(x["normalized_estimate"]))>.1000000001
                    for z in curve if z["feasible"]=="True")
        assert (x["profile_precision"]=="outside_supported")==outside
        successful=all(int(v["counts_recovered"])>=90 for v in results if v["parameter"]==p)
        assert (r["counts_pass"]=="True")== (x["profile_precision"]=="narrow_checked" and successful)
    write(HERE/"final_validation.tsv",proof)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--validate",action="store_true")
    args=parser.parse_args()
    summarize()
    if args.validate:
        validate()
