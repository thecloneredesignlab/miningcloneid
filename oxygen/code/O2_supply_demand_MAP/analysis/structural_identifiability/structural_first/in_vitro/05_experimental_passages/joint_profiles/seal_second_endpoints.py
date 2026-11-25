"""Independently validate the second endpoint release and report tested-value brackets."""

import argparse
from pathlib import Path
from math import isclose
from scan import Run, read, write, cpu, digest
from endpoint_zero_review import SOURCE, PREVIOUS
from second_endpoint_prepare import FIRST, SETUP_SECONDS, CALIBRATION_SECONDS, seed_manifest, specification
from second_endpoint_workers import state, PER_ANCHOR
from seal_endpoint_release import check_bounds
import calibrate_bootstrap as batch


def verify_recovery(out):
    amended={r["seed"]:r for r in read(out/"recovery_pairs.tsv")}
    assert len(amended)==4 and all(int(r["recovery_fits"])<=1 for r in amended.values())
    reviewed=read(out/"recovery_seeds.tsv")
    assert len(reviewed)==len({r["seed"] for r in reviewed})==398
    for a in read(out/"recovery_calibration.tsv"):
        name=a["parameter"]
        before={r["seed"]:r for r in batch.read_pairs(FIRST,name)}
        rs=[r for r in reviewed if r["parameter"]==name]
        for r in rs:
            original=before[r["seed"]]
            assert r["observations_sha256"]==original["observations_sha256"]
            if r["seed"] in amended:
                update=amended[r["seed"]]
                q=Path(name)/f"replicate_{int(r['replicate']):03d}"
                assert original["status"]=="unresolved"
                assert digest(out/"recovery"/q/"observations.tsv")==digest(FIRST/"pairs"/q/"observations.tsv")
                assert update["nested"]=="True" and float(update["independent_likelihood_error"])<2e-7
                if update["status"]=="verified":
                    assert min(int(update[arm+"_qualified_starts"]) for arm in ("restricted","unrestricted"))>=2
            else:
                update=original
            assert r["statistic"]==update["statistic"] and r["status"]==update["status"]
        check_bounds(a,[float(r["statistic"]) for r in rs if r["status"]=="verified"])


def verify_anchors(out):
    registry=[]
    trials=[]
    checks=[]
    run=Run(SOURCE)
    worker_costs=read(out/"worker_cpu.tsv")
    old={int(r["seed"]) for packet in (FIRST,PREVIOUS) for r in read(packet/"seed_manifest.tsv")}
    old.update(int(r["seed"]) for r in read(SOURCE/"bootstrap_pilot_pairs.tsv"))
    allseeds=set()
    allstates=set()
    slots_found=[]
    for slot in range(1,7):
        folder=out/f"anchor_{slot}"
        if not folder.exists():continue
        slots_found.append(slot)
        batch.check_seal(folder/"FREEZE.sha256")
        obs,=read(folder/"observed_receipt.tsv")
        expected_name,expected_value=specification(run,slot,out)
        assert obs["parameter"]==expected_name and float(obs["fixed_value"])==expected_value
        if obs["status"]!="Verified":
            assert not (folder/"pairs").exists()
            continue
        assert abs(float(obs["baseline_replay_delta"]))<=1e-6 and int(obs["qualifying_starts"])>=2
        _,name,rr,partials=state(out,slot)
        srows=read(folder/"seed_manifest.tsv")
        assert srows==[{k:str(v) for k,v in r.items()} for r in seed_manifest(name,slot)]
        aggregate,=read(folder/"outcome.tsv")
        observed,=read(folder/"anchors.tsv")
        assert aggregate["parameter"]==name and float(aggregate["fixed_value"])==float(obs["fixed_value"])
        assert abs(float(observed["observed_signed_rise"])-float(observed["independent_likelihood_rise"]))<2e-7
        assert len(rr)==int(aggregate["processed"])
        first=[r for r in rr if int(r["replicate"])<=25]
        gate=read(folder/"wave_gate.tsv")
        if gate:
            g,=gate
            assert sum(r["status"]=="verified" for r in first)==int(g["first_wave_verified"])
            assert len(first)==int(g["first_wave_pairs"])
            passed=(len(first)==25 and int(g["first_wave_verified"])>=23 and
                    float(g["projected_anchor_cpu_seconds"])<=PER_ANCHOR and g["paired_receipts_verified"]=="True")
            assert (g["release_remaining"]=="True")==passed
            assert passed or all(int(r["replicate"])<=25 for r in rr)
        byseed={r["seed"]:r for r in rr}
        part={read(p)[0]["seed"] for p in partials}
        assert not part&byseed.keys()
        assert len({r["seed"] for r in rr})==len(rr)
        for s in srows:
            seed=int(s["seed"])
            assert seed not in old|allseeds and s["state_sha256"] not in allstates
            allseeds.add(seed)
            allstates.add(s["state_sha256"])
            r=byseed.get(s["seed"])
            rel=folder/"pairs"/name/f"replicate_{int(s['replicate']):03d}"
            dispatches=read(rel/"dispatches.tsv")
            if r:
                batch.pair_receipt(folder,r)
                assert int(r["replicate"])==int(s["replicate"]) and r["seed"]==s["seed"]
                if r["status"]=="verified":
                    assert min(int(r["restricted_qualified_starts"]),int(r["unrestricted_qualified_starts"]))>=2
            elif rel.exists():
                assert s["seed"] in part
            trials.append(dict(slot=slot,parameter=name,value=obs["fixed_value"],replicate=s["replicate"],seed=seed,
                   state_sha256=s["state_sha256"],observations_sha256=(read(rel/"observations_receipt.tsv") or [{}])[0].get("arrays_sha256",""),
                   status=r["status"] if r else ("Interrupted" if s["seed"] in part else "Not run"),
                   statistic=r["statistic"] if r else "",optimizer_dispatches=len(dispatches),
                   reserve_dispatches=sum(d["phase"]=="reserve" for d in dispatches)))
        check_bounds(aggregate,[float(r["statistic"]) for r in rr if r["status"]=="verified"])
        registry.append(dict(slot=slot,**aggregate))
        checks.append(dict(slot=slot,parameter=name,qualified_observed_starts=obs["qualifying_starts"],
              wave_verified=int(gate[0]["first_wave_verified"]) if gate else 0,
              wave_processed=len(first),full_pairs_checked=len(rr),paired_receipts_checked=True,
              wave_forecast_core_hours=float(gate[0]["projected_anchor_cpu_seconds"])/3600 if gate else "",
              worker_core_hours=sum(float(row["cpu_seconds"]) for row in worker_costs if int(row["slot"])==slot)/3600,
              finish_reason=aggregate["finish_reason"],
              unknowns_in_tail=int(aggregate["unresolved"])+int(aggregate["unrun"])))
    assert sum(s>4 for s in slots_found)<=2
    enabled={int(r["slot"]) for r in read(out/"reserve_release.tsv") if r["enabled"]=="True"}
    assert enabled=={s for s in slots_found if s>4}
    write(out/"seed_dispositions.tsv",trials)
    write(out/"anchor_audit.tsv",checks)
    write(out/"new_calibrations.tsv",registry)
    return registry,trials


def brackets(out,new):
    # Endpoint positions are evaluated on the registered coordinate, while all
    # bounds in the table are also reported in their original physical units.
    run=Run(SOURCE)
    conf=run.conf
    base=run.baseline_x()
    evidence={}
    for r in read(PREVIOUS/"calibration.tsv"):
        if r["parameter"] in {"n_O","k_o_mis"}:
            assert r["assessment"]=="Tail probability above 0.05"
            assert min(float(r[col]) for col in ("tail_lower_add_one","monte_carlo_95_lower"))>.05
            evidence[r["parameter"],float(r["fixed_value"])]=("Retained","Earlier calibrated point")
    for packet,kind in ((read(FIRST/"calibration_summary.tsv"),"First release"),
                         (read(out/"recovery_calibration.tsv"),"Frozen-data review"),(new,"Second release")):
        for r in packet:
            evidence[r["parameter"],float(r["fixed_value"])]=(r["decision"],kind)
    results=[]
    for s in read(FIRST/"endpoint_register.tsv"):
        name=s["parameter"]
        center=conf.natural(base,name)
        cx=conf.target_coordinate(name,center)
        pts=[(v,d,label) for (p,v),(d,label) in evidence.items() if p==name]
        for side,k in (("Lower",-1),("Upper",1)):
            candidates=sorted([q for q in pts if k*(conf.target_coordinate(name,q[0])-cx)>1e-10],
                                   key=lambda q:abs(conf.target_coordinate(name,q[0])-cx))
            excluded=None
            kept=center
            reentry=[]
            unresolved=[]
            for v,d,label in candidates:
                if d=="Excluded" and excluded is None:
                    excluded=v
                elif d=="Retained" and excluded is None:
                    kept=v
                elif d=="Retained":
                    reentry.append(v)
                elif d=="Inconclusive":
                    unresolved.append(v)
            declared=float(s["domain_lower" if side=="Lower" else "domain_upper"])
            span="" if excluded is None else abs(conf.target_coordinate(name,excluded)-conf.target_coordinate(name,kept))
            bound_source="Replayed unrestricted optimum" if kept==center else "Calibrated retained point"
            if kept==center and any(d=="Retained" and abs(conf.target_coordinate(name,v)-cx)<1e-10 for v,d,_ in pts):
                bound_source="Calibrated retained point at fitted reference"
            if abs(conf.target_coordinate(name,declared)-cx)<1e-10:
                status="Fitted reference at declared domain boundary"
            elif excluded is not None:
                status=("Sampled bracket at target resolution" if span<=.02+1e-12 else "Sampled bracket remains wider than target")
            elif isclose(declared,kept,rel_tol=1e-11,abs_tol=1e-15):
                status="Calibrated retained point reaches declared boundary"
            elif kept!=center:
                status="Retained point; outer exclusion open"
            else:
                status="Outer exclusion open"
            results.append(dict(parameter=name,side=side,fitted_reference=center,registered_bound=declared,
                    supported_side=kept,supported_evidence=bound_source,excluded_side=excluded if excluded is not None else "",
                    bracket_lower=min(kept,excluded) if excluded is not None else "",
                    bracket_upper=max(kept,excluded) if excluded is not None else "",
                    normalized_width=span,unresolved_tested_values=" ".join(f"{v:.15g}" for v in unresolved),
                    possible_reentry=" ".join(map(str,reentry)),status=status,scope="Component adjacent to fitted reference at sampled resolution"))
    assert len(results)==28
    write(out/"endpoint_brackets.tsv",results)
    write(out/"tested_values.tsv",[dict(parameter=n,fixed_value=x,decision=d,evidence=e)
         for (n,x),(d,e) in sorted(evidence.items())])
    return results


def resources(out,trials):
    ledger=read(out/"worker_cpu.tsv")
    assert all(int(r["exit_code"])==0 for r in ledger)
    spent=0.
    for allocation in sorted((out/"dispatch").glob("*.tsv"),key=lambda p:p.stat().st_mtime_ns):
        reserved=read(allocation)
        completions=[r for r in ledger if r["stage"]==allocation.stem]
        assert len(reserved)==len(completions)<=4
        assert spent+sum(int(r["cpu_reservation"]) for r in reserved)<=CALIBRATION_SECONDS
        for r in completions:
            request=reserved[int(r["worker"])]
            assert (request["slot"],request["cpu_reservation"])==(r["slot"],r["reserved"])
        spent+=sum(float(r["cpu_seconds"]) for r in completions)
    assert abs(spent-sum(float(r["cpu_seconds"]) for r in ledger))<1e-4
    assert spent<=CALIBRATION_SECONDS
    for sample in {r["slot"] for r in ledger}:
        assert sum(float(r["cpu_seconds"]) for r in ledger if r["slot"]==sample)<=PER_ANCHOR
    setup=sum(float(read(p)[0]["cpu_seconds"]) for p in [out/"recovery_cpu.tsv",*out.glob("prepare_cpu_*.tsv")])
    controller=sum(float(read(p)[0]["cpu_seconds"]) for p in out.glob("*_controller_cpu.tsv"))
    assert setup<=SETUP_SECONDS
    final=cpu()
    # Covers the seconds after the audit ledger is written and lightweight status checks.
    allowance=90
    total=spent+setup+controller+final
    assert controller+final+allowance<=3600 and total+allowance<=25.46*3600
    row=dict(preparation_cpu_seconds=setup,paired_worker_cpu_seconds=spent,coordination_cpu_seconds=controller,
         final_audit_cpu_seconds_before_sealing=final,measured_core_hours=total/3600,
         seal_and_status_allowance_seconds=allowance,charged_core_hours=(total+allowance)/3600,
         remaining_core_hours=25.46-(total+allowance)/3600,approved_core_hours=25.46,
         fresh_optimizer_dispatches=sum(r["optimizer_dispatches"] for r in trials),
         fresh_reserve_dispatches=sum(r["reserve_dispatches"] for r in trials),
         recovery_fits=sum(int(r["recovery_fits"]) for r in read(out/"recovery_pairs.tsv")),
         local_scientific_process_limit=4,numerical_threads_each=1)
    write(out/"resources.tsv",[row])
    return row


def report(out,new,ranges,cost):
    items=["# Second release of three-readout culture support calibrations","",
        "Conditional on strict donor inheritance and the registered observation model, each anchor tests a fixed value using paired restricted and unrestricted fits to the same synthetic experiment. The data comprise 225 visible-count events, 195 conditional death-classification events and six verified chromosome specimens (104 cells). A five-percent exclusion or retention requires both the failure-inclusive add-one bounds and the exact binomial Monte Carlo envelope to lie strictly on the same side of 0.05.","",
        "| Source | Parameter | Fixed value | Observed rise | Verified / unresolved / unrun | Exceedances | Add-one bounds | 95% MC envelope | Decision |",
        "|---|---|---:|---:|---:|---:|---|---|---|"]
    for label,rows in (("Frozen-data review",read(out/"recovery_calibration.tsv")),("New",new)):
        for r in rows:
            items.append(f'| {label} | {r["parameter"]} | {float(r["fixed_value"]):.10g} | {float(r["observed_statistic"]):.6f} | {r["verified"]} / {r["unresolved"]} / {r["unrun"]} | {r["verified_exceedances"]} | {float(r["tail_lower_add_one"]):.3f}–{float(r["tail_upper_add_one"]):.3f} | {float(r["monte_carlo_95_lower"]):.5f}–{float(r["monte_carlo_95_upper"]):.5f} | {r["decision"]} |')
    items+=["","## Lower edges near the fitted reference","",
            "| Mechanism | Supported-side reference | Excluded point | Width in registered coordinate | Status |","|---|---:|---:|---:|---|"]
    for r in ranges:
        if r["side"]=="Lower" and r["parameter"] in ("p_misseg","lam_max","buffer_beta","buffer_smax","p_mis_base"):
            excluded=f'{r["excluded_side"]:.9g}' if r["excluded_side"]!="" else "Open"
            width=f'{100*r["normalized_width"]:.3f}%' if r["normalized_width"]!="" else "Open"
            items.append(f'| {r["parameter"]} | {r["supported_side"]:.9g} | {excluded} | {width} | {r["status"]} |')
    items+=["","Each bracket applies to the sampled component adjacent to the fitted value. The fitted reference reproduces the unrestricted objective. Undecided intermediate points remain open; the reported width uses each mechanism’s registered fitting coordinate. The two mechanism zeros retain their first-release exclusion decisions.","",
            "## Release checks","",
            "| Slot | Parameter | First wave verified / processed | Forecast (core-hours) | Actual workers (core-hours) | Terminal pairs checked |",
            "|---:|---|---:|---:|---:|---:|"]
    for r in read(out/"anchor_audit.tsv"):
        forecast=f'{float(r["wave_forecast_core_hours"]):.3f}' if r["wave_forecast_core_hours"] else "Unavailable"
        items.append(f'| {r["slot"]} | {r["parameter"]} | {r["wave_verified"]} / {r["wave_processed"]} | {forecast} | {float(r["worker_core_hours"]):.3f} | {r["full_pairs_checked"]} |')
    items+=["",f'The local run measured {cost["measured_core_hours"]:.4f} core-hours before final sealing. The charged total is {cost["charged_core_hours"]:.4f} core-hours, including a {cost["seal_and_status_allowance_seconds"]}-second allowance. The remaining balance is {cost["remaining_core_hours"]:.4f} of 25.46 core-hours.',"",
          "[All 14 parameters and both directions](endpoint_brackets.tsv) · [Tested values](tested_values.tsv) · [New anchor results](new_calibrations.tsv) · [Frozen pair recovery](recovery_pairs.tsv) · [Fresh seed ledger](seed_dispositions.tsv) · [Reserve decisions](reserve_release.tsv) · [CPU account](resources.tsv)"]
    (out/"REPORT.md").write_text("\n".join(items)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path)
    out=p.parse_args().output
    assert (out/"reserve_controller_cpu.tsv").is_file() and not (out/"SHA256SUMS").exists()
    batch.check_seal(out/"RUN.sha256")
    batch.check_seal(FIRST/"SHA256SUMS")
    verify_recovery(out)
    outcomes,draws=verify_anchors(out)
    rows=brackets(out,outcomes)
    write(out/"final_audit_source.tsv",[dict(module=Path(__file__).name,sha256=digest(Path(__file__)))])
    cost=resources(out,draws)
    report(out,outcomes,rows,cost)
    count=batch.seal_results(out)
    assert cpu()-cost["final_audit_cpu_seconds_before_sealing"]<=cost["seal_and_status_allowance_seconds"]
    print(f"Verified {count} sealed files. Charged CPU {cost['charged_core_hours']:.4f} core-hours.")
