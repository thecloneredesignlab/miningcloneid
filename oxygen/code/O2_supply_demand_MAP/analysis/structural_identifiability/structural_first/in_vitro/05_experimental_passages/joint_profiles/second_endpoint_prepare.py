"""Recover four frozen fit pairs and freeze qualified endpoint generators for release two."""

import argparse
from pathlib import Path
from shutil import copyfile
import calibrate_bootstrap as batch
from bootstrap_pilot import bound_observations, verify_likelihood_scale
from endpoint_zero_review import SOURCE, PREVIOUS, assess
from recover_bootstrap_pairs import Case, Review, ARMS, review_arm
from scan import ARCHIVES, Run, cpu, read, write, digest, apply_fixed, exact_margins, independent_score, known
import numpy as np

FIRST = ARCHIVES / "20261124_joint_support_first_release_01"
SETUP_SECONDS = 5400
CALIBRATION_SECONDS = int(22.96*3600)
TARGETS = ("p_misseg","lam_max","buffer_beta","buffer_smax")
RECOVERY = {"buffer_beta": {162: "restricted"},
            "p_mis_base": {6: "unrestricted",170: "unrestricted",176: "unrestricted"}}


def recover(out):
    batch.install_cpu_cap(SETUP_SECONDS-60)
    out.mkdir(exist_ok=False)
    for packet in (FIRST,SOURCE,PREVIOUS):
        batch.check_seal(packet/"SHA256SUMS")
    for item in read(FIRST/"code_receipt.tsv"):
        assert digest(Path(__file__).with_name(item["code"]))==item["sha256"]
    seeds=read(FIRST/"seed_manifest.tsv")
    revised=[]
    reviews=[]
    outcomes=[]
    for name, queue in RECOVERY.items():
        fit=batch.PairedRefits(FIRST,SETUP_SECONDS-90,name)
        originals=batch.read_pairs(FIRST,name)
        assert len(originals)==199 and {int(r["replicate"]) for r in originals if r["status"]=="unresolved"}==set(queue)
        by_seed={s["seed"]:s for s in seeds if s["parameter"]==name}
        patched=[]
        for row in originals:
            j=int(row["replicate"])
            if j not in queue:
                patched.append(row)
                continue
            arm=queue[j]
            batch.pair_receipt(FIRST,row)
            rel=Path(name)/f"replicate_{j:03d}"
            folder=out/"recovery"/rel
            folder.mkdir(parents=True,exist_ok=False)
            prev=FIRST/"pairs"/rel
            for f in ("observations.tsv","observations_receipt.tsv"):
                copyfile(prev/f, folder/f)
                assert digest(folder/f)==digest(prev/f)
            seed=by_seed[row["seed"]]
            arrays,generated,ohash=batch.dataset(fit.model,fit.anchor,seed,folder)
            assert generated is None and ohash==row["observations_sha256"]
            c=Case(dict(seed,value=fit.value),arrays,folder)
            fit.out=folder
            fit.rows.clear()
            fit.steps.clear()
            fit.new=[]
            fit.fits=0
            with bound_observations(fit.model,arrays):
                for old in read(prev/"refit_endpoints.tsv"):
                    x=fit.conf.vector(old)
                    assert np.array_equal(x,apply_fixed(x,fit.conf.constraints(*c.fit_target(old["arm"]))))
                    p=fit.conf.unmap(x)[0]
                    independent_score(fit.model,p,old)
                    assert abs(float(min(exact_margins(fit.model,p)))-float(old["minimum_margin"]))<1e-11
                    fit.rows[c.seed,old["arm"]].append(old)
                initial=fit.audit_pair(c)
                assert initial["reason"]==row["reason"]==arm+" agreement"
                assert len(fit.status(c,arm)[1])==1
                start=Review.candidates(fit,c,arm,"reserve",1)
                if start:
                    fit.current_start=start
                    write(folder/"recovery_start.tsv",[dict(arm=arm,source=start[0],
                          **{"start_"+k:v for k,v in zip(fit.conf.names,start[1])})])
                    fit.optimize(c,arm,"reserve",1)
                    fit.current_start=None
                result=fit.audit_pair(c)
                assert fit.fits<=1 and result["independent_likelihood_error"]<2e-7
            record=dict(parameter=name,seed=row["seed"],replicate=j,observations_sha256=ohash,
                        reviewed_arm=arm,prior_status=row["status"],recovery_fits=fit.fits,**result)
            reviews.append(record)
            write(folder/"review.tsv",[record])
            patched.append(record)
        assert len(patched)==199
        anchor,=[a for a in read(FIRST/"anchors.tsv") if a["parameter"]==name]
        revised.append(assess(name,fit.value,float(anchor["observed_statistic"]),patched))
        outcomes.extend(dict(parameter=name,seed=r["seed"],replicate=r["replicate"],status=r["status"],
                    statistic=r["statistic"],observations_sha256=r["observations_sha256"]) for r in patched)
    write(out/"recovery_pairs.tsv",reviews)
    write(out/"recovery_seeds.tsv",outcomes)
    write(out/"recovery_calibration.tsv",revised)
    write(out/"recovery_cpu.tsv",[dict(cpu_seconds=cpu(),reviewed_pairs=len(reviews),
         additional_fits=sum(r["recovery_fits"] for r in reviews),generated_draws=0)])
    print("Four frozen paired fits reviewed.",flush=True)


def midpoint(conf,base,name,lo,hi):
    x=base.copy()
    x[conf.idx[name]]=(conf.target_coordinate(name,lo)+conf.target_coordinate(name,hi))/2
    return conf.natural(x,name)


def seed_manifest(name,slot):
    root=2026112402
    return [dict(parameter=name,root_seed=root,batch_size=batch.B,anchor_code=slot,replicate=j,
                 seed=root*10000+slot*1000+j,rng="PCG64 via NumPy SeedSequence",
                 state_sha256=batch.state_hash([root,batch.B,slot,j])) for j in range(1,batch.B+1)]


def specification(run,slot,out):
    base=run.baseline_x()
    ref={r["parameter"]:float(r["fixed_value"]) for r in read(FIRST/"fresh_calibration.tsv")}
    names=TARGETS+("lam_max","p_mis_base")
    name=names[slot-1]
    if slot==6:
        m,=read(out/"anchor_1"/"outcome.tsv")
        if m["decision"]=="Inconclusive":
            name="p_misseg"
    if slot==3:
        value=.85
    elif slot==4:
        value=.7133175907678697
    elif slot==6 and name=="p_mis_base":
        value=.0019407667236782133
    else:
        left=ref[name]
        right=run.conf.natural(base,name)
        if slot>4:
            previous,=read(out/("anchor_2" if slot==5 else "anchor_1")/"outcome.tsv")
            v=float(previous["fixed_value"])
            # For an undecided midpoint search the half nearest the supported reference.
            # It remains an undecided point and does not replace either bracket edge.
            if previous["decision"]=="Retained":
                right=v
            else:
                left=v
        value=midpoint(run.conf,base,name,left,right)
    return name,value


def prepare(out,slot):
    run=Run(SOURCE)
    base=run.baseline_x()
    conf,model=run.conf,run.model
    assert model.contract=="J" and (model.n,len(model.dead),model.kary.shape,int(model.kary.sum()))==(225,195,(6,133),104)
    assert np.array_equal(model.indices[model.death_log_index],model.classifier_slots)
    b=read(SOURCE/"baseline.tsv")[0]
    assert b["status"]=="verified" and abs(independent_score(model,conf.unmap(base)[0])["loss"]-float(b["selected_loss"]))<1e-7
    name,value=specification(run,slot,out)
    focus=conf.natural(base,name)
    xref=apply_fixed(base,conf.constraints(name,focus))
    assert np.max(abs(xref-base))<1e-12
    _,dr=verify_likelihood_scale(model,base,xref,conf)
    assert abs(dr)<=1e-6
    folder=out/f"anchor_{slot}"
    folder.mkdir(exist_ok=False)
    fixed=conf.constraints(name,value)
    available=read(SOURCE/"profiles"/name/"attempts.tsv")
    points=[p for p in read(SOURCE/"profile_grid.tsv") if p["parameter"]==name and p["status"]=="verified"]
    obs=folder/"observed"
    previous=[r for r in available if r["value"] and float(r["value"])==value]
    if previous:
        write(obs/"attempts.tsv",previous)
    else:
        lookup={r["attempt"]:r for r in available}
        for k in range(1,5):
            extra=out/("anchor_2" if name=="lam_max" else "anchor_1")
            starts=[]
            for direction in ((-1,1) if k%2 else (1,-1)):
                side=[p for p in points if (float(p["value"])-value)*direction>0]
                if side:
                    near=min(side,key=lambda p:abs(conf.target_coordinate(name,float(p["value"]))-conf.target_coordinate(name,value)))
                    starts.append(("Sealed neighbor "+near["value"],conf.vector(lookup[near["attempt"]])))
            if slot>4 and (extra/"generator_vectors.tsv").is_file():
                starts.insert(0,("Primary midpoint generator",conf.vector(read(extra/"generator_vectors.tsv")[0])))
            records=read(obs/"attempts.tsv")
            floor,qualified=review_arm(run,records,fixed)
            if len(qualified)>=2:
                break
            if k>2 and floor:
                starts.insert(0,("Best feasible constrained fit",conf.vector(floor)))
            starts.extend([("Measured joint baseline",base)]+known(conf))
            used=[run.start(r) for r in records if r["loss"]]
            run.execute(obs,f"observed_{slot}_{k}",name,value,starts,used,k-1)
    records=read(obs/"attempts.tsv")
    floor,qualified=review_arm(run,records,fixed)
    reason="Verified" if len(qualified)>=2 else "Observed fit agreement open"
    obs_receipt=dict(parameter=name,fixed_value=value,normalized=conf.target_coordinate(name,value),
                     fitted_reference=focus,reference_normalized=conf.target_coordinate(name,focus),
                     baseline_replay_delta=dr,qualifying_starts=len(qualified),
                     feasible_floor=floor["loss"] if floor else "",status=reason)
    if len(qualified)>=2:
        scores=[]
        for q in qualified:
            x=apply_fixed(conf.vector(q),fixed)
            assert np.max(abs(x-conf.vector(q)))<1e-12
            checked,_=run.accept(x,run.start(q),fixed,name,value,q["attempt"],"Observed restricted","Replay",{},saved=q)
            assert checked["outcome"]=="qualified"
            scores.append(checked)
        selected=scores[0]
        x=apply_fixed(conf.vector(selected),fixed)
        ell,rise=verify_likelihood_scale(model,base,x,conf)
        if rise < -1e-6 or float(floor["loss"])<float(b["selected_loss"])-1e-6:
            raise RuntimeError("A lower restricted score requires baseline review before calibration")
        anchor=dict(parameter=name,fixed_value=value,attempt=selected["attempt"],observed_statistic=max(0.,rise),
                    observed_signed_rise=rise,independent_likelihood_rise=ell,
                    observed_restricted_J=selected["loss"],frozen_log_count_scale=selected["sigma"])
        write(folder/"anchors.tsv",[anchor])
        write(folder/"generator_vectors.tsv",[dict(parameter=name,**{"x_"+key:v for key,v in zip(conf.names,x)})])
        write(folder/"baseline_vector.tsv",[dict(("x_"+key,v) for key,v in zip(conf.names,base))])
        seeds=seed_manifest(name,slot)
        old=read(FIRST/"seed_manifest.tsv")+read(PREVIOUS/"seed_manifest.tsv")+read(SOURCE/"bootstrap_pilot_pairs.tsv")
        assert len({r["state_sha256"] for r in seeds})==len({r["seed"] for r in seeds})==199
        assert not ({r["state_sha256"] for r in seeds}&{r.get("state_sha256") for r in old})
        assert not ({r["seed"] for r in seeds}&{int(r["seed"]) for r in old})
        write(folder/"seed_manifest.tsv",seeds)
    write(folder/"observed_receipt.tsv",[obs_receipt])
    batch.write_seal(folder/"FREEZE.sha256",[p for p in folder.rglob("*.tsv")])
    print(f"Slot {slot} {name} at {value:.16g}: {reason}; {len(qualified)} qualifying starts.",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=("recover","prepare"))
    p.add_argument("output",type=Path)
    p.add_argument("slots",type=int,nargs="*")
    p.add_argument("--cpu-seconds",type=int,default=SETUP_SECONDS-60)
    a=p.parse_args()
    assert a.output.resolve().parent==ARCHIVES.resolve()
    if a.action=="recover":
        recover(a.output)
    else:
        batch.install_cpu_cap(a.cpu_seconds)
        assert len(read(a.output/"recovery_pairs.tsv"))==4
        for slot in a.slots:
            prepare(a.output,slot)
        write(a.output/("prepare_cpu_"+"_".join(map(str,a.slots))+".tsv"),[dict(cpu_seconds=cpu())])
