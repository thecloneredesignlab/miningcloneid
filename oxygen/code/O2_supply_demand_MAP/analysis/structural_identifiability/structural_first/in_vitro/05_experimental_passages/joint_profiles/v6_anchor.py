"""Prepare and calibrate one conditional V6 endpoint using the frozen V5 scorer."""

import math
import resource
from pathlib import Path

from scan import Run, apply_fixed, cpu, independent_score, known, np, read, write
import calibrate_bootstrap as batch
from bootstrap_pilot import verify_likelihood_scale
from endpoint_bootstrap_release import PER_ANCHOR, gate
from endpoint_zero_review import SOURCE, assess
from recover_bootstrap_pairs import BudgetEnd, review_arm
from seal_endpoint_release import check_bounds


def seeds(parameter, batch_number, slot):
    # Each batch has a distinct entropy root and each slot a disjoint public seed range.
    root = 2030000000 + batch_number
    rows = []
    for j in range(1, batch.B + 1):
        rows.append(dict(parameter=parameter, root_seed=root, batch_size=batch.B,
                         anchor_code=slot, replicate=j, seed=root*10000+slot*1000+j,
                         rng="PCG64 via NumPy SeedSequence",
                         state_sha256=batch.state_hash([root,batch.B,slot,j])))
    return rows


def prepare(folder, spec, seed_rows, near_vectors):
    run = Run(SOURCE)
    conf, model = run.conf, run.model
    name, value = spec["parameter"], float(spec["value"])
    base = run.baseline_x()
    baseline, = read(SOURCE/"baseline.tsv")
    assert model.contract == "J" and (model.n,len(model.dead),model.kary.shape,int(model.kary.sum())) == (225,195,(6,133),104)
    assert np.array_equal(model.indices[model.death_log_index], model.classifier_slots)
    assert abs(independent_score(model,conf.unmap(base)[0])["loss"]-float(baseline["selected_loss"])) < 1e-7
    assert baseline["status"] == "verified" and int(baseline["agreeing_starts"]) >= 2
    fixed = conf.constraints(name,value)
    archive = read(SOURCE/"profiles"/name/"attempts.tsv")
    obs = folder/"observed"
    if not (obs/"attempts.tsv").exists():
        write(obs/"attempts.tsv",[r for r in archive if r["value"] and float(r["value"]) == value])
    points = [p for p in read(SOURCE/"profile_grid.tsv") if p["parameter"] == name and p["status"] == "verified"]
    lookup = {r["attempt"]:r for r in archive}
    neighbors = [(float(p["normalized"]),conf.vector(lookup[p["attempt"]])) for p in points]
    neighbors.extend((float(r["coordinate"]),conf.vector(read(Path(r["vector"]))[0])) for r in near_vectors)
    center = conf.target_coordinate(name,value)
    for k in range(4):
        records = read(obs/"attempts.tsv")
        floor, qualified = review_arm(run,records,fixed)
        if len(qualified) >= 2:
            break
        starts = []
        for direction in ((-1,1) if k%2 == 0 else (1,-1)):
            side = [r for r in neighbors if (r[0]-center)*direction > 1e-12]
            if side:
                starts.append(("Qualified neighbor on registered coordinate",min(side,key=lambda r:abs(r[0]-center))[1]))
        if floor and k >= 2:
            starts.insert(0,("Lowest feasible constrained endpoint",conf.vector(floor)))
        starts.extend([("Unrestricted joint fit",base)]+known(conf))
        run.execute(obs,f'v6_{spec["batch"]}_{spec["slot"]}_{k+1}',name,value,starts,
                    [run.start(r) for r in records if r["loss"]],k)
    floor, qualified = review_arm(run,read(obs/"attempts.tsv"),fixed)
    if floor and float(floor["loss"]) < float(baseline["selected_loss"])-1e-6:
        raise RuntimeError("Restricted floor improves the unrestricted fit. Review the reference before continuing.")
    receipt = dict(parameter=name,fixed_value=value,qualifying_starts=len(qualified),
                   feasible_floor=floor["loss"] if floor else "",baseline=baseline["selected_loss"])
    if len(qualified) >= 2:
        for q in qualified:
            x = conf.vector(q)
            assert np.max(abs(x-apply_fixed(x,fixed))) < 1e-12
            check,_ = run.accept(x,run.start(q),fixed,name,value,q["attempt"],"Restricted endpoint","Replay",{},saved=q)
            assert check["outcome"] == "qualified"
        q = qualified[0]
        x = conf.vector(q)
        independent,rise = verify_likelihood_scale(model,base,x,conf)
        assert rise >= -1e-6 and abs(independent-rise) < 1e-6
        write(folder/"anchors.tsv",[dict(parameter=name,fixed_value=value,attempt=q["attempt"],
              observed_statistic=max(0.,rise),observed_signed_rise=rise,independent_likelihood_rise=independent,
              observed_restricted_J=q["loss"],frozen_log_count_scale=q["sigma"])])
        write(folder/"generator_vectors.tsv",[dict(parameter=name,**{"x_"+n:v for n,v in zip(conf.names,x)})])
        write(folder/"baseline_vector.tsv",[dict(("x_"+n,v) for n,v in zip(conf.names,base))])
    write(folder/"observed_receipt.tsv",[receipt])
    write(folder/"seed_manifest.tsv",seed_rows)
    batch.write_seal(folder/"FREEZE.sha256",list(folder.rglob("*.tsv")))


def calibrate(folder):
    spec, = read(folder/"design.tsv")
    name, value = spec["parameter"], float(spec["value"])
    batch.check_seal(folder/"FREEZE.sha256")
    frozen = read(folder/"seed_manifest.tsv")
    expected = seeds(name,int(spec["batch"]),int(spec["slot"]))
    assert frozen == [{k:str(v) for k,v in r.items()} for r in expected]
    observed, = read(folder/"anchors.tsv")
    assert not (folder/"pairs").exists(), "Existing draw receipts require reconciliation before a new dispatch"
    started = cpu()
    old_soft,hard = resource.getrlimit(resource.RLIMIT_CPU)
    resource.setrlimit(resource.RLIMIT_CPU,(math.floor(started+PER_ANCHOR-2),hard))
    try:
        fit = batch.PairedRefits(folder,started+PER_ANCHOR-75,name)
        records = []
        reason = "All frozen pairs completed"
        for seed in frozen:
            if sum(r["status"] != "verified" for r in records) >= batch.FAIL_STOP:
                reason = "Ten unresolved pairs"
                break
            place = folder/"pairs"/name/('replicate_%03d'%int(seed["replicate"]))
            place.mkdir(parents=True,exist_ok=False)
            try:
                fit.guard()
                fit.pair(seed,place)
            except BudgetEnd:
                write(place/"partial.tsv",[dict(seed=seed["seed"],reason="Anchor CPU allowance reached")])
                reason = "Anchor CPU allowance reached"
                break
            row, = read(place/"pair.tsv")
            assert row["seed"] == seed["seed"] and float(row["fixed_value"]) == value
            batch.pair_receipt(folder,row)
            records.append(row)
            if len(records) == batch.WAVE:
                decision = gate(sum(r["status"] == "verified" for r in records),len(records),
                                cpu()-started,[float(r["pair_cpu_seconds"]) for r in records])
                write(folder/"wave_gate.tsv",[decision])
                if not decision["release_remaining"]:
                    reason = "First-wave quality or cost gate held"
                    break
        result = assess(name,value,float(observed["observed_statistic"]),records)
        check_bounds(result,[float(r["statistic"]) for r in records if r["status"] == "verified"])
        charge = cpu()-started
        assert charge <= PER_ANCHOR
        # Use the same failure-inclusive bounds as the parent release. A quality
        # or cost hold remains visible even if the conservative test is decisive.
        write(folder/"calibration.tsv",[dict(**result,disposition=result["decision"],
                                             finish_reason=reason,paired_cpu_seconds=charge)])
    finally:
        resource.setrlimit(resource.RLIMIT_CPU,(old_soft,hard))
