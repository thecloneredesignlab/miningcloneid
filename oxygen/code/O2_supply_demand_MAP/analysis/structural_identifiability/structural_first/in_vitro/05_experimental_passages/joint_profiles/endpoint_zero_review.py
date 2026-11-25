"""Register the 14 conditional support searches and review frozen exact-null pairs."""

from pathlib import Path
from time import monotonic
import argparse
import shutil
import calibrate_bootstrap as batch
from bootstrap_pilot import bound_observations
from recover_bootstrap_pairs import ARMS, Case, active_distance, feasible_rows
from scan import ARCHIVES, MATERIAL_SEPARATION, SUSPENDED, apply_fixed, cpu, digest, exact_margins, independent_score, read, write
import numpy as np

SOURCE = batch.SOURCE
PREVIOUS = ARCHIVES / "20261123_bootstrap_J_4anchors_B199_01"
ANCHORS = {
    "p_misseg": 0.008276322281071904,
    "buffer_beta": 0.9906665556049095,
    "p_mis_base": 1e-6,
    "lam_max": 1.0926307319717656,
}
# The planning values name grid points approximately; all released generators use ANCHORS.
SEARCH = {
    "lam_max": ("1.092631", "Fitted upper bound", "Eight low points ruled out by donor supply"),
    "p_mis_base": ("0.000001", "Fitted upper bound", "Test contact across the flat observed profile"),
    "p_misseg": ("0.008276; exact zero separately", "0.016634 then 0.025", "Exact-zero carrier k_o_mis"),
    "k_o_mis": ("Fitted lower bound; calibrated", "0.000126", "Bound anchor does not locate the upper edge"),
    "p_wgd": ("Fitted lower bound", "0.000599 and 0.001", "Test upper domain contact"),
    "buffer_smax": ("0.713318", "0.888889 then 1", "Test upper domain contact"),
    "buffer_beta": ("0.990667; exact zero separately", "1.601755", "Exact-zero carrier buffer_n_exp"),
    "buffer_n_exp": ("2.719465", "3.593814 to 5.994843", "Wide upper grid spacing"),
    "alpha_o2": ("Fitted lower bound", "0.561231 and 0.629961", "Two high values have unresolved feasibility"),
    "gamma_growth": ("1.5", "2.166667", "Test lower domain contact"),
    "mu_hp": ("Fitted lower bound", "0.019574 and 0.030850", "The candidate at 0.30 is provisional"),
    "gamma_mu": ("Fitted lower bound", "3.5", "Test upper domain contact"),
    "O2_crit": ("1.424328", "Fitted upper bound", "Refine toward estimate after first classification"),
    "n_O": ("2.222278 and 1.666733", "3.888911", "2.704919 is calibrated; 0.555644 is provisional"),
}


def register(out):
    assert out.resolve().parent == ARCHIVES.resolve()
    out.mkdir(exist_ok=False)
    domain = {r["parameter"]: r for r in read(SOURCE/"domain.tsv")}
    summary = read(SOURCE/"mechanism_summary.tsv")
    assert len(summary) == 14 and {r["parameter"] for r in summary} == SEARCH.keys()
    rows = []
    for r in summary:
        name = r["parameter"]
        lower, upper, note = SEARCH[name]
        rows.append(dict(parameter=name, domain_lower=domain[name]["natural_lower"],
                         domain_upper=domain[name]["natural_upper"], estimate=r["estimate"],
                         fitted_bound=r["estimated_box_face"] or "Interior",
                         lower_search=lower, upper_search=upper, note=note,
                         first_release_value=ANCHORS.get(name,""), endpoint_status="Not determined"))
    write(out/"endpoint_register.tsv",rows)
    write(out/"scope.tsv",[dict(contract="J", inheritance="Strict donor inheritance with fixed proxies",
          observed_counts=225, death_events=195, karyotype_specimens=6, karyotyped_cells=104,
          unresolved_zero_pairs=12, fresh_anchors=4, frozen_draws_per_anchor=199,
          first_wave=25, minimum_passing_wave_pairs=23, unresolved_anchor_stop=10,
          maximum_science_processes=4, numerical_threads_per_process=1,
          recovery_and_preparation_core_hours=8, core_hours_per_anchor=7, total_core_hours=36)])
    for pack in (SOURCE, PREVIOUS):
        batch.check_seal(pack/"SHA256SUMS")
    batch.check_seal(PREVIOUS/"FREEZE.sha256")
    source_rows=[dict(packet=p.name, manifest_sha256=digest(p/"SHA256SUMS")) for p in (SOURCE,PREVIOUS)]
    write(out/"source_packets.tsv",source_rows)


def assess(name, value, observed, pairs):
    stats=[float(r["statistic"]) for r in pairs if r["status"]=="verified"]
    unknown=batch.B-len(stats)
    k=sum(t>=observed for t in stats)
    t=batch.tail_summary(k,unknown)
    if max(t["tail_upper_add_one"],t["monte_carlo_95_upper"])<.05:
        decision="Excluded"
    elif min(t["tail_lower_add_one"],t["monte_carlo_95_lower"])>.05:
        decision="Retained"
    else:
        decision="Inconclusive"
    return dict(parameter=name, fixed_value=value, observed_statistic=observed,
                frozen=batch.B, processed=len(pairs),verified=len(stats),
                unresolved=len(pairs)-len(stats), unrun=batch.B-len(pairs),
                verified_exceedances=k, **t, **batch.critical_summary(stats,unknown),
                decision=decision)


def recovery_start(fit, case, arm, serial):
    conf, model=fit.conf,fit.model
    fixed=conf.constraints(*case.fit_target(arm))
    rows=feasible_rows(fit.rows[case.seed,arm])
    qualified=[r for r in rows if r["outcome"]=="qualified"]
    seeds=[conf.vector(r) for r in qualified+rows]
    old=[fit.run.start(r) for r in rows]
    target=conf.idx[case.target]
    carrier=conf.idx[SUSPENDED[case.target]]
    trials=[]
    for x in seeds:
        # At a free null face, test release at alternate carrier values.
        if arm=="unrestricted" and x[target]<1e-5:
            for k in (.75,.4,1.):
                y=x.copy()
                y[target],y[carrier]=.003,k
                trials.append(("Released null face at a distinct carrier",y))
        for key in ("rho","shared_count_loading","class_log_odds"):
            y=x.copy()
            i=conf.idx[key]
            y[i]+= (.015 if serial==0 else .032)*(1 if y[i]<.75 else -1)
            trials.append(("Perturbed " + key + " from a saved endpoint",y))
        trials.append(("Saved feasible endpoint",x))
    for label,y in trials:
        y=apply_fixed(y,fixed)
        if all(active_distance(conf,y,v,fixed)>MATERIAL_SEPARATION for v in old):
            if min(exact_margins(model,conf.unmap(y)[0]))>=0:
                return label,y
    return None


def recover(out):
    # One process, one bounded recovery round, at most the number of missing
    # independent starts (two when an arm previously had none).
    batch.install_cpu_cap(7*3600)
    start=monotonic()
    source_pairs=read(PREVIOUS/"pair_results.tsv")
    seeds={r["seed"]:r for r in read(PREVIOUS/"seed_manifest.tsv")}
    unresolved=[r for r in source_pairs if r["parameter"] in SUSPENDED and r["status"]=="unresolved"]
    assert len(unresolved)==12 and len(set(r["seed"] for r in unresolved))==12
    assert sorted(sum(r["parameter"]==name for r in unresolved) for name in SUSPENDED)==[4,8]
    journal=[]
    for r in unresolved:
        batch.pair_receipt(PREVIOUS,r)
        journal.append(dict(parameter=r["parameter"],seed=r["seed"],replicate=r["replicate"],
                            prior_reason=r["reason"],observations_sha256=r["observations_sha256"]))
    write(out/"zero_recovery_queue.tsv",journal)
    fits={n:batch.PairedRefits(PREVIOUS,7*3600-60,n) for n in SUSPENDED}
    patches={}
    for old in unresolved:
        name,seed=old["parameter"],old["seed"]
        source=PREVIOUS/"pairs"/name/f"replicate_{int(old['replicate']):03d}"
        folder=out/"zero_review"/name/f"replicate_{int(old['replicate']):03d}"
        folder.mkdir(parents=True,exist_ok=False)
        for f in ("observations.tsv","observations_receipt.tsv"):
            shutil.copyfile(source/f,folder/f)
            assert digest(source/f)==digest(folder/f)
        fit=fits[name]
        arrays,pi,sha=batch.dataset(fit.model,fit.anchor,seeds[seed],folder)
        assert pi is None and sha==old["observations_sha256"]
        c=Case(dict(seeds[seed],value=0.),arrays,folder)
        fit.out=folder
        fit.rows.clear()
        fit.fits=0
        fit.new=[]
        with bound_observations(fit.model,arrays):
            replay=[]
            for r in read(source/"refit_endpoints.tsv"):
                fixed=fit.conf.constraints(*c.fit_target(r["arm"]))
                x=fit.conf.vector(r)
                assert np.array_equal(x, apply_fixed(x,fixed))
                p=fit.conf.unmap(x)[0]
                independent_score(fit.model,p,r)
                assert abs(min(exact_margins(fit.model,p))-float(r["minimum_margin"]))<=1e-11
                fit.rows[c.seed,r["arm"]].append(r)
                replay.append(dict(attempt=r["attempt"],arm=r["arm"],loss=r["loss"],
                                   original_qualification=r["outcome"],independent_replay="Passed"))
            write(folder/"prior_replay.tsv",replay)
            prior=fit.audit_pair(c)
            assert prior["status"]==old["status"] and prior["reason"]==old["reason"]
            for arm in ARMS:
                missing=2-len(fit.status(c,arm)[1])
                for serial in range(max(0,missing)):
                    choice=recovery_start(fit,c,arm,serial)
                    if choice is None: break
                    fit.current_start=choice
                    write(folder/"recovery_dispatch.tsv",[dict(seed=seed,arm=arm,serial=serial+1,source=choice[0],
                          **{"start_"+k:v for k,v in zip(fit.conf.names,choice[1])})],append=True)
                    # Reserve solver settings, different labels from all archived fits.
                    fit.optimize(c,arm,"reserve",serial+1)
                    fit.current_start=None
                    if fit.ready(c,arm): break
            result=fit.audit_pair(c)
            assert not result["independent_likelihood_error"] or result["independent_likelihood_error"]<2e-7
        new=dict(parameter=name,seed=seed,replicate=old["replicate"],observations_sha256=sha,
                 prior_status=old["status"],prior_reason=old["reason"],recovery_fits=fit.fits,**result)
        patches[seed]=new
        write(folder/"review.tsv",[new])
        write(out/"zero_pair_reviews.tsv",list(patches.values()))
    outcomes=[]
    revised=[]
    for anchor in read(PREVIOUS/"anchors.tsv"):
        name=anchor["parameter"]
        if name not in SUSPENDED: continue
        sample=[patches.get(r["seed"],r) for r in source_pairs if r["parameter"]==name]
        assert len(sample)==199 and len({r["seed"] for r in sample})==199
        outcomes.extend(dict(parameter=name,seed=r["seed"],replicate=r["replicate"],status=r["status"],
                        statistic=r["statistic"],source="Frozen-data recovery" if r["seed"] in patches else "Sealed prior pair")
                        for r in sample)
        revised.append(assess(name,0.,float(anchor["observed_statistic"]),sample))
    write(out/"zero_seed_dispositions.tsv",outcomes)
    write(out/"zero_calibration.tsv",revised)
    write(out/"zero_resources.tsv",[dict(process_cpu_seconds=cpu(),elapsed_seconds=monotonic()-start,
          recovered_pairs=len(patches),new_optimizer_calls=sum(r["recovery_fits"] for r in patches.values()),
          generated_datasets=0)])
    print("Frozen zero-pair reviews complete.",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    register(args.output)
    recover(args.output)
