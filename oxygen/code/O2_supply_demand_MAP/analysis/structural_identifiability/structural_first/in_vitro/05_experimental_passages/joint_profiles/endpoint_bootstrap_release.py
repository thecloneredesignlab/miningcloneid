"""First endpoint release: four value-specific paired bootstrap calibrations."""

import argparse
from contextlib import redirect_stdout
from math import ceil
from pathlib import Path
from time import monotonic
import fcntl
import os
import subprocess
import sys

import calibrate_bootstrap as batch
from endpoint_zero_review import ANCHORS, PREVIOUS, SOURCE, assess
from scan import ARCHIVES, Run, apply_fixed, cpu, digest, exact_margins, independent_score, read, write
from bootstrap_pilot import verify_likelihood_scale
from recover_bootstrap_pairs import BudgetEnd
import numpy as np

PER_ANCHOR = 7 * 3600
ROOT_SEED = 20261124
MINIMUM = 23


def make_seeds():
    seeds=[]
    for i,name in enumerate(ANCHORS,1):
        for j in range(1,batch.B+1):
            state=[ROOT_SEED,batch.B,i,j]
            seeds.append(dict(parameter=name,root_seed=ROOT_SEED,batch_size=batch.B,anchor_code=i,
                         replicate=j,seed=ROOT_SEED*10000+i*1000+j,rng="PCG64 via NumPy SeedSequence",
                         state_sha256=batch.state_hash(state)))
    return seeds


def gate(passed, processed, seconds_used, pair_times):
    rate=1.25*float(np.quantile(pair_times,.75)) if pair_times else 0.
    projected=seconds_used+rate*(batch.B-processed)+30
    quality=processed==batch.WAVE and passed>=MINIMUM
    return dict(first_wave_pairs=processed,first_wave_verified=passed,
                first_wave_unresolved=processed-passed,minimum_verified=MINIMUM,
                quality_gate=quality,spent_cpu_seconds=seconds_used,
                planning_cpu_seconds_per_pair=rate,projected_anchor_cpu_seconds=projected,
                limit_cpu_seconds=PER_ANCHOR,cost_gate=projected<=PER_ANCHOR,
                release_remaining=quality and projected<=PER_ANCHOR)


def freeze(out):
    assert (out/"zero_resources.tsv").is_file() and len(read(out/"zero_pair_reviews.tsv"))==12
    assert not (out/"FREEZE.sha256").exists()
    batch.reproducibility_control()
    assert gate(23,25,100,[1]*25)["release_remaining"]
    assert not gate(22,25,100,[1]*25)["release_remaining"]
    assert not gate(25,25,100,[200]*25)["release_remaining"]
    assert assess("test",0,1,[dict(status="verified",statistic=0)]*191+
                               [dict(status="unresolved")]*8)["decision"]=="Inconclusive"
    finalization,=read(PREVIOUS/"finalization.tsv")
    lineage=[]
    for row in read(PREVIOUS/"code_receipt.tsv"):
        actual=digest(Path(__file__).with_name(row["code"]))
        expected=row["sha256"]
        if row["code"]=="calibrate_bootstrap.py":
            assert finalization["original_runtime_sha256"]==expected
            assert finalization["extra_optimizer_invocations"]==finalization["extra_bootstrap_draws"]=="0"
            expected=finalization["corrected_sealing_code_sha256"]
        assert actual==expected,row["code"]
        lineage.append(dict(module=row["code"],source_runtime=row["sha256"],audited_current=actual))
    write(out/"source_code_lineage.tsv",lineage)
    # Verify the bound libraries, observation contract and prior data receipts.
    batch.check_sources(out)
    batch.anonymize_provenance(out)
    run=Run(SOURCE)
    m,c=run.model,run.conf
    assert (m.n,len(m.dead),m.kary.shape,int(m.kary.sum()))==(225,195,(6,133),104)
    assert m.contract=="J" and len(c.targets)==14 and np.array_equal(m.indices[m.death_log_index],m.classifier_slots)
    base=run.baseline_x()
    baseline=read(SOURCE/"baseline.tsv")[0]
    assert abs(independent_score(m,c.unmap(base)[0])["loss"]-float(baseline["selected_loss"]))<1e-7
    assert baseline["status"]=="verified" and int(baseline["agreeing_starts"])>=2
    grid=read(SOURCE/"profile_grid.tsv")
    anchors=[]
    vectors=[]
    for name,value in ANCHORS.items():
        point,=[r for r in grid if r["parameter"]==name and float(r["value"])==value]
        assert point["status"]=="verified" and int(point["agreeing_starts"])>=2
        path=SOURCE/"profiles"/name/"attempts.tsv"
        endpoint,=[r for r in read(path) if r["attempt"]==point["attempt"]]
        fixed=c.constraints(name,value)
        x=apply_fixed(c.vector(endpoint),fixed)
        assert max(abs(x-c.vector(endpoint)))<1e-12
        checked,_=run.accept(x,run.start(endpoint),fixed,name,value,endpoint["attempt"],
                         "Sealed observed-data endpoint","Replay",{},saved=endpoint)
        assert checked["outcome"]==endpoint["outcome"]=="qualified"
        ell,rise=verify_likelihood_scale(m,base,x,c)
        assert abs(rise-float(point["rise"]))<1e-6
        anchors.append(dict(parameter=name,fixed_value=value,attempt=endpoint["attempt"],
                      observed_restricted_J=checked["loss"],observed_signed_rise=rise,
                      observed_statistic=max(0.,rise),independent_likelihood_rise=ell,
                      frozen_log_count_scale=checked["sigma"],minimum_margin=checked["minimum_margin"],
                      endpoint_file_sha256=digest(path),profile_grid_sha256=digest(SOURCE/"profile_grid.tsv")))
        vectors.append(dict(parameter=name,**{"x_"+n:float(v) for n,v in zip(c.names,x)}))
    seeds=make_seeds()
    old=read(PREVIOUS/"seed_manifest.tsv")+read(SOURCE/"bootstrap_pilot_pairs.tsv")
    assert len(seeds)==len({s["seed"] for s in seeds})==len({s["state_sha256"] for s in seeds})==796
    assert not ({s["seed"] for s in seeds}&{int(s["seed"]) for s in old})
    assert not ({s["state_sha256"] for s in seeds}&{s.get("state_sha256") for s in old})
    write(out/"anchors.tsv",anchors)
    write(out/"generator_vectors.tsv",vectors)
    write(out/"baseline_vector.tsv",[dict((f"x_{n}",float(v)) for n,v in zip(c.names,base))])
    write(out/"seed_manifest.tsv",seeds)
    write(out/"domain.tsv",read(SOURCE/"domain.tsv"))
    run.out=out
    run.exact_zero_preflight(base)
    code=[Path(__file__),Path(__file__).with_name("endpoint_zero_review.py")]
    code.extend(Path(__file__).with_name(r["code"]) for r in read(PREVIOUS/"code_receipt.tsv"))
    write(out/"code_receipt.tsv",[dict(code=p.name,sha256=digest(p)) for p in code])
    write(out/"preflight.tsv",[dict(registered_parameters=14,observed_generators_replayed=4,frozen_seeds=796,
          fitted_zero_pairs_reviewed=12,independent_likelihood_replay="Passed",quality_and_cost_controls="Passed",
          frozen_inputs="Three-readout J under strict donor inheritance")])
    batch.write_seal(out/"FREEZE.sha256",[p for p in out.glob("*.tsv")])


def verify_input(out):
    batch.check_seal(out/"FREEZE.sha256")
    for r in read(out/"code_receipt.tsv"):
        assert digest(Path(__file__).with_name(r["code"]))==r["sha256"]
    frozen=read(out/"seed_manifest.tsv")
    assert frozen==[{k:str(v) for k,v in r.items()} for r in make_seeds()]
    return frozen


def anchor_job(out,name):
    started=cpu()
    deadline=started+PER_ANCHOR-30
    frozen=verify_input(out)
    fit=batch.PairedRefits(out,deadline,name)
    records=[]
    decision=None
    reason="All 199 pairs completed"
    for seed in (s for s in frozen if s["parameter"]==name):
        if sum(r["status"]!="verified" for r in records)>=batch.FAIL_STOP:
            reason="Ten unresolved pairs"
            break
        j=int(seed["replicate"])
        folder=out/"pairs"/name/f"replicate_{j:03d}"
        folder.mkdir(parents=True,exist_ok=False)
        try:
            fit.guard()
            fit.pair(seed,folder)
        except BudgetEnd:
            reason="Anchor CPU allocation exhausted"
            break
        record,=read(folder/"pair.tsv")
        assert seed["seed"]==record["seed"]
        batch.pair_receipt(out,record)
        records.append(record)
        if j==batch.WAVE:
            decision=gate(sum(r["status"]=="verified" for r in records),j,cpu()-started,
                          [float(r["pair_cpu_seconds"]) for r in records])
            write(out/f"gate_{name}.tsv",[dict(parameter=name,integrity_gate=True,**decision)])
            if not decision["release_remaining"]:
                reason="First-wave quality or cost gate held the remaining pairs"
                break
    if decision is None:
        decision=gate(sum(r["status"]=="verified" for r in records),len(records),cpu()-started,
                      [float(r["pair_cpu_seconds"]) for r in records])
        write(out/f"gate_{name}.tsv",[dict(parameter=name,integrity_gate=True,**decision)])
    write(out/f"close_{name}.tsv",[dict(parameter=name,processed=len(records),cpu_seconds=cpu()-started,
          finish_reason=reason)])
    print(f"{name}: {reason}. Processed {len(records)}.",flush=True)


def collect(out):
    frozen=verify_input(out)
    ledger=[]
    bounds=[]
    for a in read(out/"anchors.tsv"):
        name=a["parameter"]
        rows=batch.read_pairs(out,name)
        for r in rows:
            batch.pair_receipt(out,r)
        lookup={r["seed"]:r for r in rows}
        assert len(lookup)==len(rows)
        for s in (s for s in frozen if s["parameter"]==name):
            r=lookup.get(s["seed"],{})
            ledger.append(dict(**s,status=r.get("status","Not run"),statistic=r.get("statistic",""),
                      observations_sha256=r.get("observations_sha256","")))
        bounds.append(assess(name,a["fixed_value"],float(a["observed_statistic"]),rows))
    write(out/"fresh_seed_results.tsv",ledger)
    write(out/"fresh_calibration.tsv",bounds)
    return bounds


def start_release(out):
    start=monotonic()
    assert out.resolve().parent==ARCHIVES.resolve()
    with (out/".controller.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert not (out/"dispatch.tsv").exists()
        freeze(out)
        # The coordinator fits the fourth anchor itself. There are at most
        # four numerical processes, including the coordinator, at every point.
        names=list(ANCHORS)
        write(out/"dispatch.tsv",[dict(parameter=n,phase="First release",maximum_pairs=199,
                    first_wave=25,minimum_first_wave_verified=23,maximum_core_hours=7) for n in names])
        processes=[]
        for name in names[:-1]:
            with (out/f"worker_{name}.log").open("x") as f:
                process=subprocess.Popen([sys.executable,"-Bu",__file__,"worker",str(out),name],
                                    stdout=f,stderr=subprocess.STDOUT)
            processes.append((name,process))
        name=names[-1]
        t0=cpu()
        with (out/f"worker_{name}.log").open("x") as f, redirect_stdout(f):
            anchor_job(out,name)
        local=cpu()-t0
        charged=[dict(parameter=name,worker_cpu_seconds=local,exit_code=0)]
        for name,proc in processes:
            _,status,ru=os.wait4(proc.pid,0)
            proc.returncode=os.waitstatus_to_exitcode(status)
            charged.append(dict(parameter=name,worker_cpu_seconds=ru.ru_utime+ru.ru_stime,exit_code=proc.returncode))
        write(out/"worker_resources.tsv",charged)
        assert all(r["exit_code"]==0 for r in charged)
        bounds=collect(out)
        zero=read(out/"zero_resources.tsv")[0]
        worker=sum(r["worker_cpu_seconds"] for r in charged)
        overhead=cpu()-local
        measured=float(zero["process_cpu_seconds"])+worker+overhead
        setup_allowance=sum(float(r["charge_seconds"]) for r in read(out/"setup_incidents.tsv"))
        total=measured+setup_allowance
        assert all(r["worker_cpu_seconds"]<=PER_ANCHOR for r in charged)
        assert total<=36*3600
        write(out/"resources.tsv",[dict(zero_recovery_and_register_cpu_seconds=zero["process_cpu_seconds"],
              fresh_worker_cpu_seconds=worker,initialization_collection_and_coordinator_cpu_seconds=overhead,
              measured_cpu_seconds=measured,preflight_incident_allowance_seconds=setup_allowance,
              total_charged_seconds=total,core_hours=total/3600,approved_core_hours=36,
              remaining_core_hours=36-total/3600,elapsed_fresh_seconds=monotonic()-start,
              peak_scientific_processes=4,numerical_threads_each=1)])
        for r in bounds:
            print(r["parameter"],r["processed"],r["decision"])
        print(f"First release complete. CPU: {total/3600:.5f} core-hours.",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("run","worker"))
    parser.add_argument("output",type=Path)
    parser.add_argument("parameter",choices=list(ANCHORS),nargs="?")
    a=parser.parse_args()
    if a.action=="run":
        start_release(a.output)
    else:
        assert a.parameter and (a.output/"dispatch.tsv").is_file()
        batch.install_cpu_cap(PER_ANCHOR)
        anchor_job(a.output,a.parameter)
