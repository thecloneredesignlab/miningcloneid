"""Freeze and run four conditional J bootstrap calibrations using local CPU only."""

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import sha256
from math import ceil, floor, sqrt
from pathlib import Path
from time import monotonic
import fcntl
import os
import resource
import subprocess
import sys

from scan import (ARCHIVES, PILOT, Run, apply_fixed, check_sources, cpu, digest,
                  exact_margins, independent_score, read, write)
from bootstrap_pilot import NEAR_EQUAL, bound_observations, simulate, verify_likelihood_scale
from recover_bootstrap_pairs import (ARMS, BudgetEnd, Case, Review, SOURCE,
                                     active_distance, feasible_rows, verify_packet)
from scan import MATERIAL_SEPARATION
import numpy as np
from scipy.stats import beta, binom

B = 199
WAVE = 25
FAIL_STOP = 10
ROOT_SEED = 20261123
CPU_LIMIT = 28 * 3600
COMPLETION_CPU = 120
INIT_CPU = 60
REVIEW_PACKET = ARCHIVES / "20261123_bootstrap_frozen_pair_review_01"
ANCHOR_ATTEMPTS = (
    "archive_J_p_misseg_grid_1_reserve",
    "archive_J_buffer_beta_grid_1_reserve",
    "k_o_mis_1.0000000000000001e-05_outward",
    "n_O_2.7049188028242881_downward",
)
READOUTS = ("log_count", "dead_objects", "karyotype_cells")


def seal(path):
    return [(h, name) for h, name in
            (line.split(maxsplit=1) for line in path.read_text().splitlines())]


def write_seal(path, files):
    tmp = path.with_suffix(".tmp")
    tmp.write_text("".join(f"{digest(f)}  {f.relative_to(path.parent)}\n" for f in sorted(files)))
    tmp.replace(path)


def check_seal(path):
    entries = seal(path)
    assert entries
    for expected, name in entries:
        assert digest(path.parent / name) == expected, name
    return len(entries)


def anonymize_provenance(out):
    receipt = out / "source_receipt.tsv"
    rows = read(receipt)
    home = str(Path.home()) + "/"
    for row in rows:
        if row["source"].startswith(home):
            row["source"] = "~/" + row["source"].removeprefix(home)
    write(receipt, rows)


def seal_results(out):
    ignored = {"SHA256SUMS", ".controller.lock", ".DS_Store"}
    files = [p for p in out.rglob("*") if p.is_file() and p.name not in ignored]
    assert all(f.suffix in (".tsv", ".sha256", ".md", ".log") for f in files)
    write_seal(out / "SHA256SUMS", files)
    return check_seal(out / "SHA256SUMS")


def state_hash(seed):
    integers = np.random.SeedSequence(seed).generate_state(4)
    return sha256(np.asarray(integers, dtype="<u4").tobytes()).hexdigest()


def entropy(row):
    return [int(row[key]) for key in ("root_seed", "batch_size", "anchor_code", "replicate")]


def make_seeds():
    return [dict(parameter=name, root_seed=ROOT_SEED, batch_size=B, anchor_code=i,
                 replicate=j, seed=ROOT_SEED*10000+i*1000+j,
                 rng="PCG64 via NumPy SeedSequence",
                 state_sha256=state_hash([ROOT_SEED, B, i, j]))
            for i, name in enumerate(PILOT, 1) for j in range(1, B+1)]


def compare_manifest(rows):
    expected = [{k: str(v) for k, v in row.items()} for row in make_seeds()]
    assert rows == expected and len({r["state_sha256"] for r in rows}) == B*len(PILOT)
    assert not ({r["seed"] for r in rows} & {r["seed"] for r in read(SOURCE/"bootstrap_pilot_pairs.tsv")})
    return rows


def reproducibility_control():
    assert tail_summary(0, 0)["tail_upper_add_one"] == .005
    assert tail_summary(0, 10)["tail_upper_add_one"] == .055
    assert tail_summary(199, 0)["tail_lower_add_one"] == 1.
    values = list(range(199))
    assert critical_summary(values, 0)["critical_lower"] == 189
    assert critical_summary(values[10:], 10)["critical_upper"] == "Unbounded"
    assert critical_summary(values[:-1], 1)["critical_lower"] == 188


def install_cpu_cap(seconds):
    _, hard = resource.getrlimit(resource.RLIMIT_CPU)
    if hard == resource.RLIM_INFINITY or hard > seconds:
        resource.setrlimit(resource.RLIMIT_CPU, (seconds, seconds))


def freeze(out):
    install_cpu_cap(INIT_CPU)
    out.mkdir(exist_ok=False)
    reproducibility_control()
    verify_packet()
    review_count = check_seal(REVIEW_PACKET / "SHA256SUMS")
    review = read(REVIEW_PACKET / "pair_review.tsv")
    assert len(review)==12 and all(r["final_status"]=="verified" and r["nesting_pass"]=="True" for r in review)
    gate, = read(REVIEW_PACKET / "batch_gate.tsv")
    assert gate["numeric_gate_pass"] == "True" and gate["projected_core_hours_4_by_199"]=="28"
    check_sources(out)
    anonymize_provenance(out)
    run = Run(SOURCE)
    conf, model = run.conf, run.model
    base = run.baseline_x()
    assert len(conf.targets)==14 and model.contract=="J"
    assert np.array_equal(model.indices[model.death_log_index],model.classifier_slots)
    assert model.y.shape==(225,) and model.dead.shape==(195,) and model.kary.shape==(6,133)
    assert model.kary.sum()==104 and len(exact_margins(model, conf.unmap(base)[0]))==110
    run.out = out
    run.exact_zero_preflight(base)
    anchors, vectors = [], []
    for name, identifier in zip(PILOT, ANCHOR_ATTEMPTS):
        generator, = read(SOURCE/"bootstrap_pilot"/name/"generator.tsv")
        assert generator["parameter"]==name and generator["attempt"]==identifier
        path = SOURCE/"profiles"/name/"attempts.tsv"
        terminal, = [r for r in read(path) if r["attempt"]==identifier]
        assert terminal["outcome"]=="qualified"
        value = float(generator["value"])
        x = conf.vector(terminal)
        assert np.max(abs(x-apply_fixed(x,conf.constraints(name,value))))<=1e-13
        x = apply_fixed(x,conf.constraints(name,value))
        p = conf.unmap(x)[0]
        independent = independent_score(model,p,terminal)
        ll_delta, signed = verify_likelihood_scale(model,base,x,conf)
        minimum = float(min(exact_margins(model,p)))
        assert minimum>=0.
        assert abs(independent["sigma"]-float(generator["log_count_scale"]))<=1e-10
        assert abs(signed-float(generator["score_difference"]))<=1e-6
        assert abs(independent["loss"]-float(generator["objective"]))<=1e-6
        anchors.append(dict(parameter=name,fixed_value=value,attempt=identifier,
                            observed_restricted_J=independent["loss"],observed_signed_rise=signed,
                            observed_statistic=max(0.,signed),independent_likelihood_rise=ll_delta,
                            frozen_log_count_scale=independent["sigma"],minimum_donor_margin=minimum,
                            generator_source_sha256=digest(SOURCE/"bootstrap_pilot"/name/"generator.tsv"),
                            endpoint_source_sha256=digest(path)))
        vectors.append(dict(parameter=name,**{"x_"+key:float(v) for key,v in zip(conf.names,x)}))
    write(out/"baseline_vector.tsv",[dict(**{"x_"+key:float(v) for key,v in zip(conf.names,base)})])
    write(out/"anchors.tsv",anchors)
    write(out/"generator_vectors.tsv",vectors)
    write(out/"domain.tsv",read(SOURCE/"domain.tsv"))
    write(out/"seed_manifest.tsv",make_seeds())
    compare_manifest(read(out/"seed_manifest.tsv"))
    files = [Path(__file__), *(Path(__file__).with_name(name) for name in
             ("bootstrap_pilot.py", "recover_bootstrap_pairs.py", "scan.py"))]
    write(out/"code_receipt.tsv",[dict(code=p.name,sha256=digest(p)) for p in files])
    write(out/"initialization.tsv",[dict(verified_previous_pairs=12,prior_review_sealed_files=review_count,
         seeds_frozen=796,anchors=4,counts=225,death=195,specimens=6,karyotyped_cells=104,
         conditioned_inocula=118,donor_constraints=110,cpu_seconds=cpu(),
         computational_controls="Passed seed, threshold, objective, nesting-domain and exact-zero input tests")])
    write_seal(out/"FREEZE.sha256",[p for p in out.iterdir() if p.is_file()])
    print("Frozen four generators and 796 new seeds. No bootstrap data were drawn during initialization.",flush=True)


def dataset(model,p,seed_row,folder):
    selections = (model.indices, model.classifier_slots, model.kary_rows)
    sample = folder/"observations.tsv"
    receipt = folder/"observations_receipt.tsv"
    rng_hash = seed_row["state_sha256"]
    assert state_hash(entropy(seed_row))==rng_hash
    reference_fraction = None
    if not sample.is_file():
        arrays, reference_fraction = simulate(model,p,np.random.SeedSequence(entropy(seed_row)))
        obs_hash = sha256(b"".join(a.tobytes() for a in arrays)).hexdigest()
        observations = [dict(readout="metadata",row=seed_row["seed"],event=rng_hash,
                             autosome_bin="",value=obs_hash)]
        for name, ids, a in zip(READOUTS,selections,arrays):
            observations.extend(dict(readout=name,row=i,event=model.slot_ids[slot],
                                     autosome_bin=(j+22 if a.ndim==2 else ""),value=value)
                                for i,slot in enumerate(ids)
                                for j,value in enumerate(a.reshape(len(ids),-1)[i]))
        write(sample,observations)
    rows=read(sample)
    head,*observations=rows
    assert (head["readout"],head["row"],head["event"])==("metadata",seed_row["seed"],rng_hash)
    loaded=[]
    offset=0
    for name, ids, bins in zip(READOUTS,selections,(1,1,133)):
        chunk=observations[offset:offset+len(ids)*bins]
        expected=[(name,str(i),model.slot_ids[slot],str(j+22) if bins>1 else "")
                  for i,slot in enumerate(ids) for j in range(bins)]
        assert [(r["readout"],r["row"],r["event"],r["autosome_bin"]) for r in chunk]==expected
        values=np.array([float(r["value"]) for r in chunk],dtype=np.float64)
        loaded.append(values if bins==1 else values.reshape(len(ids),bins))
        offset+=len(chunk)
    assert offset==len(observations)
    assert all(np.isfinite(x).all() for x in loaded)
    assert np.all((loaded[1]>=0.)&(loaded[1]<=model.eligible)) and np.all(loaded[1] % 1 == 0.)
    assert np.all(loaded[2]>=0.) and np.all(loaded[2] % 1 == 0.)
    assert np.array_equal(loaded[2].sum(axis=1),model.kary.sum(axis=1))
    assert sha256(b"".join(a.tobytes() for a in loaded)).hexdigest()==head["value"]
    observed = dict(seed=seed_row["seed"],state_sha256=rng_hash,arrays_sha256=head["value"],
                    file_sha256=digest(sample),event_order_validated=True)
    if receipt.is_file():
        assert read(receipt)==[{k:str(v) for k,v in observed.items()}]
    else:
        write(receipt,[observed])
    return tuple(loaded),reference_fraction,head["value"]


class PairedRefits(Review):
    """Use the reviewed fit and acceptance code with a six-start per-pair ledger."""

    def __init__(self, out, deadline, name):
        super().__init__(out)
        self.deadline=deadline
        self.name=name
        frozen,=[r for r in read(out/"generator_vectors.tsv") if r["parameter"]==name]
        self.null=self.conf.vector(frozen)
        baseline,=read(out/"baseline_vector.tsv")
        self.base=self.conf.vector(baseline)
        anchor,=[r for r in read(out/"anchors.tsv") if r["parameter"]==name]
        self.value=float(anchor["fixed_value"])
        self.anchor=self.conf.unmap(self.null)[0]
        self.current_start=None

    def candidates(self,c,arm,kind,step):
        if self.current_start is not None:
            return self.current_start
        conf=self.conf
        fixed=conf.constraints(*c.fit_target(arm))
        rows=feasible_rows(self.rows[c.seed,arm])
        r=feasible_rows(self.rows[c.seed,"restricted"])
        start_columns=["start_"+key for key in conf.names]
        started=[np.array([float(v[key]) for key in start_columns])
                 for v in read(self.out/"dispatches.tsv") if v["arm"]==arm]
        labeled = []
        if kind=="main" and step==0:
            if arm=="restricted":
                labeled.append(("Frozen null generator",self.null))
            elif r:
                labeled.append(("Embedded restricted floor",conf.vector(r[0])))
        if kind=="reserve":
            if arm=="unrestricted" and r:
                labeled.append(("Return from restricted floor",conf.vector(r[0])))
            if rows:
                labeled.append(("Lowest feasible within-arm endpoint",conf.vector(rows[0])))
        if kind=="main" and step==1:
            labeled.append(("Projected measured-data unrestricted baseline",self.base))
        if rows:
            labeled.append(("Within-arm endpoint",conf.vector(rows[0])))
        labeled.extend((("Frozen null generator",self.null),("Measured-data baseline",self.base)))
        trials=list(labeled)
        for label,x in labeled:
            for key in ("rho","class_log_odds","shared_count_loading","misseg_dead_share"):
                k=conf.idx[key]
                for delta in (.008,.026):
                    y=x.copy()
                    y[k]+=delta if x[k]<.85 else -delta
                    trials.append((f"{label}, distinct {key} start",y))
        for label,x in trials:
            self.guard()
            z=apply_fixed(x,fixed)
            if all(active_distance(conf,z,previous,fixed)>MATERIAL_SEPARATION for previous in started):
                if min(exact_margins(self.model,conf.unmap(z)[0]))>=0.:
                    return label,z
        return None

    def dispatch(self,c,arm,kind,step):
        name=f"{c.seed}_{arm[0]}_{kind}_{step+1}"
        attempts=read(self.out/"dispatches.tsv")
        if any(r["attempt"]==name for r in attempts):
            return
        chosen=self.candidates(c,arm,kind,step)
        if chosen is None:
            write(self.out/"start_failures.tsv",[dict(seed=c.seed,arm=arm,phase=kind,attempt=name,
                 reason="No distinct exactly feasible start")],append=True)
            return
        source,start=chosen
        self.current_start=chosen
        write(self.out/"dispatches.tsv",[*attempts,dict(seed=c.seed,arm=arm,phase=kind,attempt=name,
            source=source,**{"start_"+key:float(v) for key,v in zip(self.conf.names,start)})])
        try:
            self.optimize(c,arm,kind,step)
        finally:
            self.current_start=None

    def audit_pair(self,c):
        model,conf=self.model,self.conf
        r,rs=self.status(c,"restricted")
        u,us=self.status(c,"unrestricted")
        rr=rs[0] if rs else r
        uu=us[0] if us else u
        nesting=[]
        for role,endpoint in (("Lowest feasible",r),("Selected qualified",rs[0] if rs else None)):
            if endpoint:
                nesting.append(dict(seed=c.seed,endpoint=role,attempt=endpoint["attempt"],
                                   **self.check_nesting(c,endpoint)))
        write(self.out/"nesting_receipts.tsv",nesting)
        signed=score_error=floor_delta=""
        nested=False
        if rr and uu:
            ell,signed=verify_likelihood_scale(model,conf.vector(uu),conf.vector(rr),conf)
            score_error=abs(ell-signed)
            floor_delta=float(r["loss"])-float(u["loss"])
            nested=signed>=-NEAR_EQUAL and floor_delta>=-NEAR_EQUAL
        qualified=[len(rs)>=2,len(us)>=2]
        failure=[]
        failure.extend(arm+" agreement" for arm,passed in zip(ARMS,qualified) if not passed)
        if not nested:
            failure.append("nesting or unavailable arm minimum")
        ok=all(qualified) and nested
        return dict(restricted_floor=r["loss"] if r else "",unrestricted_floor=u["loss"] if u else "",
            restricted_selected=rr["attempt"] if rr else "",unrestricted_selected=uu["attempt"] if uu else "",
            restricted_KKT=rr["stationarity"] if rr else "",unrestricted_KKT=uu["stationarity"] if uu else "",
            restricted_qualified_starts=len(rs),unrestricted_qualified_starts=len(us),
            signed_score_rise=signed,signed_floor_rise=floor_delta,independent_likelihood_error=score_error,
            nested=nested,statistic=max(0.,signed) if ok else "",status="verified" if ok else "unresolved",
            reason="; ".join(failure))

    def pair(self, seed_row, folder):
        begin,wall=cpu(),monotonic()
        arrays,pi,ohash=dataset(self.model,self.anchor,seed_row,folder)
        c=Case(dict(seed_row,value=self.value),arrays,folder)
        self.out=folder
        self.fits=0
        self.new=read(folder/"refit_budget_ledger.tsv")
        self.rows.clear()
        self.steps.clear()
        for row in read(folder/"refit_endpoints.tsv"):
            self.rows[c.seed,row["arm"]].append(row)
        with bound_observations(self.model,arrays):
            anchor_score=independent_score(self.model,self.anchor)
            if pi is not None:
                assert np.max(abs(anchor_score["fraction"]-pi))<=1e-12
            for step in range(2):
                self.dispatch(c,"restricted","main",step)
            for step in range(2):
                self.dispatch(c,"unrestricted","main",step)
            r,rs=self.status(c,"restricted")
            if len(rs)<2:
                self.dispatch(c,"restricted","reserve",0)
                r,rs=self.status(c,"restricted")
            u,us=self.status(c,"unrestricted")
            if len(us)<2 or (r and u and float(u["loss"])>float(r["loss"])+NEAR_EQUAL) or (
                    rs and us and float(us[0]["loss"])>float(rs[0]["loss"])+NEAR_EQUAL):
                self.dispatch(c,"unrestricted","reserve",0)
            result=self.audit_pair(c)
        dispatched=read(folder/"dispatches.tsv")
        counts=Counter((r["arm"],r["phase"]) for r in dispatched)
        assert all(n<=(1 if phase=="reserve" else 2) for (_,phase),n in counts.items())
        assert len(dispatched)<=6
        previous=sum(float(x["cpu_seconds"]) for x in read(folder/"partial_cpu.tsv"))
        record=dict(seed=seed_row["seed"],parameter=self.name,fixed_value=self.value,
                    replicate=seed_row["replicate"],observations_sha256=ohash,
                    optimizer_invocations=len(dispatched),
                    reserve_invocations=sum(n for (_,phase),n in counts.items() if phase=="reserve"),
                    pair_cpu_seconds=previous+cpu()-begin,wall_seconds=monotonic()-wall,**result)
        write(folder/"pair.tsv",[record])
        print(f"{self.name} replicate {record['replicate']} {result['status']} "
              f"after {len(dispatched)} fits, {record['pair_cpu_seconds']:.2f} CPU s",flush=True)
        return record


def read_pairs(out, name):
    return [r for j in range(1,B+1)
            for r in read(out/"pairs"/name/f"replicate_{j:03d}"/"pair.tsv")]


def worker(out, name, limit, phase):
    install_cpu_cap(limit)
    check_seal(out/"FREEZE.sha256")
    for row in read(out/"code_receipt.tsv"):
        assert digest(Path(__file__).with_name(row["code"]))==row["sha256"], row["code"]
    seeds = compare_manifest(read(out/"seed_manifest.tsv"))
    existing = read_pairs(out,name)
    known = {int(r["replicate"]) for r in existing}
    failures = sum(r["status"]=="unresolved" for r in existing)
    endpoint = WAVE if phase=="wave1" else B
    r = PairedRefits(out,limit-8.,name)
    reason="Scheduled replicates complete"
    for row in seeds:
        j=int(row["replicate"])
        if row["parameter"]!=name or j>endpoint or j in known:
            continue
        if failures>=FAIL_STOP:
            reason="Anchor stopped at ten unresolved pairs"
            break
        folder=out/"pairs"/name/f"replicate_{j:03d}"
        start=cpu()
        try:
            r.guard()
            result=r.pair(row,folder)
        except BudgetEnd as error:
            reason=str(error)
            write(folder/"partial_cpu.tsv",[dict(cpu_seconds=cpu()-start,reason=reason)],append=True)
            break
        failures+=result["status"]=="unresolved"
    write(out/f"{phase}_{name}_close.tsv",[dict(parameter=name,phase=phase,
         terminal_pairs=len(read_pairs(out,name)),unresolved=failures,stop_reason=reason)])


def tail_summary(k, unknown):
    assert 0<=k<=k+unknown<=B
    return dict(tail_lower_add_one=(1+k)/(B+1),tail_upper_add_one=(1+k+unknown)/(B+1),
                monte_carlo_95_lower=float(beta.ppf(.025,k,B-k+1)) if k else 0.,
                monte_carlo_95_upper=float(beta.ppf(.975,k+unknown+1,B-k-unknown)) if k+unknown<B else 1.,
                single_bernoulli_se_at_005=sqrt(.05*.95/B),minimum_add_one_resolution=1/(B+1))


def critical_summary(statistics, unknown):
    assert len(statistics)+unknown==B
    sorted_stats=sorted(statistics)
    assert all(np.isfinite(x) and x>=0 for x in sorted_stats)
    # Ranks for the empirical 95th percentile and a binomial order-statistic interval.
    rank=ceil(.95*(B+1))
    left=int(binom.ppf(.025,B,.95))
    right=int(binom.ppf(.975,B,.95))+1
    def lower(r):
        return 0. if r<=unknown else sorted_stats[r-unknown-1]
    def upper(r):
        return sorted_stats[r-1] if r<=len(sorted_stats) else "Unbounded"
    return dict(critical_rank=rank,critical_lower=lower(rank),critical_upper=upper(rank),
                monte_carlo_rank_lower=left,monte_carlo_rank_upper=right,
                critical_monte_carlo_95_lower=lower(left),critical_monte_carlo_95_upper=upper(right))


def pair_receipt(out, row):
    name, j = row["parameter"], int(row["replicate"])
    folder=out/"pairs"/name/f"replicate_{j:03d}"
    obs, = read(folder/"observations_receipt.tsv")
    assert obs["seed"]==row["seed"] and obs["arrays_sha256"]==row["observations_sha256"]
    assert digest(folder/"observations.tsv")==obs["file_sha256"]
    assert obs["event_order_validated"]=="True"
    fits=read(folder/"dispatches.tsv")
    endpoints=read(folder/"refit_endpoints.tsv")
    assert len(fits)==int(row["optimizer_invocations"])<=6
    assert len(set(f["attempt"] for f in fits))==len(fits)
    assert sum(f["phase"]=="reserve" for f in fits)==int(row["reserve_invocations"])<=2
    assert all(sum(f["phase"]=="reserve" and f["arm"]==arm for f in fits)<=1 for arm in ARMS)
    assert set(f["attempt"] for f in endpoints)<=set(f["attempt"] for f in fits)
    nesting=read(folder/"nesting_receipts.tsv")
    assert all(float(r["absolute_score_error"])<=1e-7 and float(r["maximum_prediction_error"])<=1e-10
               and float(r["maximum_margin_error"])<=1e-12 for r in nesting)
    if row["status"]=="verified":
        assert len(nesting)==2 and all(int(row[arm+"_qualified_starts"])>=2 for arm in ARMS)
        assert row["nested"]=="True" and float(row["signed_floor_rise"])>=-NEAR_EQUAL
        assert float(row["signed_score_rise"])>=-NEAR_EQUAL and float(row["independent_likelihood_error"])<2e-7
        assert float(row["statistic"])==max(0.,float(row["signed_score_rise"]))
    else:
        assert row["status"]=="unresolved" and not row["statistic"]
    return dict(parameter=name,replicate=j,seed=row["seed"],hash_verified=True,
                optimization_ledger_verified=True,nesting_receipts=len(nesting),
                full_likelihood_replay_verified=bool(row["independent_likelihood_error"]),status=row["status"])


def consumed(out):
    init,=read(out/"initialization.tsv")
    return ceil(float(init["cpu_seconds"]))+sum(float(r["charged_cpu_seconds"]) for r in read(out/"workers.tsv"))


def launch(out, phase, allowances):
    ledger=read(out/"workers.tsv")
    reservations=[dict(phase=phase,parameter=name,allocation_cpu_seconds=limit)
                  for name,limit in allowances.items()]
    assert len(reservations)<=4
    assert sum(allowances.values())+consumed(out)+COMPLETION_CPU+INIT_CPU<=CPU_LIMIT
    write(out/f"{phase}_allocations.tsv",reservations)
    processes=[]
    for name,limit in allowances.items():
        with (out/f"{phase}_{name}.log").open("ab") as f:
            proc=subprocess.Popen([sys.executable,"-B",__file__,"worker",str(out),name,str(limit),phase],
                                  stdout=f,stderr=subprocess.STDOUT)
        processes.append((name,limit,proc))
    live={proc.pid:(name,limit,proc) for name,limit,proc in processes}
    while live:
        pid,status,usage=os.wait4(-1,0)
        if pid not in live:
            continue
        name,limit,proc=live.pop(pid)
        proc.returncode=os.waitstatus_to_exitcode(status)
        seconds=usage.ru_utime+usage.ru_stime
        ledger.append(dict(phase=phase,parameter=name,allocation_cpu_seconds=limit,
                           measured_cpu_seconds=seconds,charged_cpu_seconds=seconds,
                           exit_code=proc.returncode))
        write(out/"workers.tsv",ledger)
        print(f"{phase} {name}: process exit {proc.returncode}, {seconds:.2f} measured CPU s.",flush=True)
    return [r for r in ledger if r["phase"]==phase]


def wave_gate(out, worker_rows):
    gate=[]
    projection=0.
    integrity=all(r["exit_code"]==0 for r in worker_rows) and len(worker_rows)==4
    for name in PILOT:
        records=read_pairs(out,name)
        receipts=[pair_receipt(out,r) for r in records]
        write(out/f"wave1_{name}_receipts.tsv",receipts)
        passed=sum(r["status"]=="verified" for r in records)
        failed=len(records)-passed
        closed=failed>=FAIL_STOP
        integrity &= bool(records) and (len(records)==WAVE or closed)
        planning=1.25*float(np.quantile([float(r["pair_cpu_seconds"]) for r in records],.75)) if records else 0.
        remaining=0 if closed else B-len(records)
        forecast=planning*remaining
        projection+=forecast
        gate.append(dict(parameter=name,processed=len(records),verified=passed,unresolved=failed,
                         stop_at_failure_limit=closed,planning_cpu_seconds_per_pair=planning,
                         remaining_pairs=remaining,projected_remaining_cpu_seconds=forecast))
    cost=projection+consumed(out)+COMPLETION_CPU+INIT_CPU
    write(out/"wave1_gate.tsv",gate)
    write(out/"wave1_decision.tsv",[dict(numerical_and_accounting_gate=integrity,
         cost_gate=cost<=CPU_LIMIT,continue_batch=integrity and cost<=CPU_LIMIT,
         used_cpu_seconds=consumed(out),projected_total_cpu_seconds=cost,
         maximum_cpu_seconds=CPU_LIMIT)])
    return integrity and cost<=CPU_LIMIT,gate


def final_report(out, reason, start):
    pairs=[r for name in PILOT for r in read_pairs(out,name)]
    verified_receipts=[pair_receipt(out,r) for r in pairs]
    write(out/"pair_results.tsv",pairs)
    write(out/"pair_verification.tsv",verified_receipts)
    key={(r["parameter"],r["replicate"]):r for r in pairs}
    assert len(key)==len(pairs)
    ledger=[]
    for seed in compare_manifest(read(out/"seed_manifest.tsv")):
        row=key.get((seed["parameter"],seed["replicate"]))
        assert not row or row["seed"]==seed["seed"]
        folder=out/"pairs"/seed["parameter"]/f"replicate_{int(seed['replicate']):03d}"
        ledger.append(dict(**seed,status=row["status"] if row else
                           ("interrupted" if (folder/"observations.tsv").is_file() else "not run"),
                           statistic=row["statistic"] if row else "",
                           observation_sha256=row["observations_sha256"] if row else ""))
    write(out/"seed_results.tsv",ledger)
    tails=[]
    for anchor in read(out/"anchors.tsv"):
        name=anchor["parameter"]
        observed=float(anchor["observed_statistic"])
        block=[r for r in pairs if r["parameter"]==name]
        stats=[float(r["statistic"]) for r in block if r["status"]=="verified"]
        failed=len(block)-len(stats)
        unrun=B-len(block)
        unknown=failed+unrun
        k=sum(v>=observed for v in stats)
        t=tail_summary(k,unknown)
        c=critical_summary(stats,unknown)
        assessment="Inconclusive at 0.05"
        if unrun:
            assessment="Incomplete planned replicate set"
        elif t["tail_upper_add_one"]<=.05 and t["monte_carlo_95_upper"]<=.05:
            assessment="Tail probability at or below 0.05"
        elif t["tail_lower_add_one"]>.05 and t["monte_carlo_95_lower"]>.05:
            assessment="Tail probability above 0.05"
        tails.append(dict(parameter=name,fixed_value=anchor["fixed_value"],observed_statistic=observed,
                          frozen=B,processed=len(block),verified=len(stats),unresolved=failed,
                          unrun_or_interrupted=unrun,verified_exceedances=k,**t,**c,assessment=assessment))
    write(out/"calibration.tsv",tails)
    worker_cpu=sum(float(r["measured_cpu_seconds"]) for r in read(out/"workers.tsv"))
    init,=read(out/"initialization.tsv")
    measured_init=float(init["cpu_seconds"])
    coordinator=cpu()
    total=worker_cpu+measured_init+coordinator
    assert total<=CPU_LIMIT and coordinator<=COMPLETION_CPU
    write(out/"resources.tsv",[dict(max_workers=4,numerical_threads_per_worker=1,
        initialization_cpu_seconds=measured_init,worker_cpu_seconds=worker_cpu,
        coordinator_cpu_seconds=coordinator,total_cpu_seconds=total,total_core_hours=total/3600.,
        optimizer_invocations_completed=sum(int(r["optimizer_invocations"]) for r in pairs),
        reserve_invocations_completed=sum(int(r["reserve_invocations"]) for r in pairs),
        cpu_limit_seconds=CPU_LIMIT,wall_seconds=monotonic()-start,finish_reason=reason)])
    report=[
        "# Conditional three-readout J bootstrap calibration",
        "",
        f"**Run status** — {reason}.",
        "",
        "The frozen observed-data anchor vectors and the recorded passage design determine",
        "four plug-in sampling distributions under strict donor inheritance. Both arms",
        "share each synthetic sample, event order, denominators and conditioned inputs.",
        "The simulated readouts are 225 log visible-object counts, 195 dead-object",
        "classifications conditional on the simulated matching counts, and six full",
        "chromosome histograms with a combined total of 104 cells. The dense likelihood",
        "replay and the objective contrast are checked on both fitted arms.",
        "",
        "| Anchor | Fixed value | Observed rise | Verified | Unresolved | Unrun | Tail bounds | 95% MC envelope |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for r in tails:
        report.append(f"| {r['parameter']} | {r['fixed_value']} | {r['observed_statistic']:.7g} | "
                      f"{r['verified']} | {r['unresolved']} | {r['unrun_or_interrupted']} | "
                      f"{r['tail_lower_add_one']:.4f}–{r['tail_upper_add_one']:.4f} | "
                      f"{r['monte_carlo_95_lower']:.4f}–{r['monte_carlo_95_upper']:.4f} |")
    report.extend([
        "",
        "An arm is accepted after two materially distinct feasible starts independently",
        "qualify within 0.02 J of the feasible floor. The original constrained",
        "stationarity, local checks and exact-zero carrier checks apply. Each arm",
        "has two main starts and at most one additional start in the reserve round.",
        "Restricted vectors are embedded at full fidelity in the unrestricted model.",
        "Both selected endpoints and feasible floors must respect score ordering",
        "within 1e-6. Accepted negative signed rises within that tolerance are set to",
        "zero for the statistic, with raw signed differences kept in the ledger.",
        "",
        "Ties count as tail exceedances. The add-one bounds include each unresolved",
        "or unrun replicate as a non-exceedance at the lower end and an exceedance",
        "at the upper end. The separate 95% Monte Carlo envelope uses the exact",
        "binomial lower limit with the confirmed exceedance count and upper limit",
        "with every unknown counted as an exceedance. Near 0.05 the binomial",
        "Monte Carlo standard error for 199 draws is about 0.0155 and the",
        "add-one resolution is 0.005. Order-statistic bounds are in calibration.tsv.",
        "A partially processed anchor has only a failure-inclusive bound on the frozen",
        "replicate set. These calibrations address the four specific plug-in anchors.",
        "",
        f"**Local CPU** — {total/3600.:.4f} core-hours measured against a 28-core-hour limit.",
        "",
        "Source and seed receipts are in FREEZE.sha256 and seed_manifest.tsv.",
        "The per-replicate decisions and resources are in seed_results.tsv,",
        "pair_results.tsv, wave1_decision.tsv, workers.tsv and resources.tsv.",
    ])
    (out/"REPORT.md").write_text("\n".join(report)+"\n")
    seal_results(out)
    print(f"{reason}. {len(pairs)} terminal pairs, {total/3600.:.4f} core-hours.",flush=True)


def controller(out):
    start=monotonic()
    check_seal(out/"FREEZE.sha256")
    with (out/".controller.lock").open("w") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (out/"wave1_allocations.tsv").exists():
            raise RuntimeError("Dispatch receipts already exist. Review the worker ledger before resuming.")
        initial=CPU_LIMIT-ceil(consumed(out))-COMPLETION_CPU-INIT_CPU
        worker_rows=launch(out,"wave1",{name:floor(initial/len(PILOT)) for name in PILOT})
        ready, gate=wave_gate(out,worker_rows)
        reason="First-wave numerical, integrity or CPU forecast gate held the remaining replicates"
        if ready:
            print("First-wave numerical, integrity and CPU projection checks passed.",flush=True)
            cpu_left=CPU_LIMIT-ceil(consumed(out))-COMPLETION_CPU-INIT_CPU
            planned={r["parameter"]:r["projected_remaining_cpu_seconds"] for r in gate
                     if r["remaining_pairs"]}
            reason="All eligible anchors reached a scheduled terminal condition"
            if planned:
                total=sum(planned.values())
                allowances={name:floor(cpu_left*v/total) for name,v in planned.items()}
                exits=launch(out,"main",allowances)
                if any(r["exit_code"]!=0 for r in exits):
                    reason="One or more local worker processes failed"
                elif any(len(read_pairs(out,name))<B and sum(r["status"]=="unresolved"
                         for r in read_pairs(out,name))<FAIL_STOP for name in planned):
                    reason="Local CPU allocation stopped an eligible anchor"
        final_report(out,reason,start)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=("freeze","run","worker","test"))
    parser.add_argument("output",type=Path,nargs="?")
    parser.add_argument("anchor",nargs="?",choices=PILOT)
    parser.add_argument("limit",type=int,nargs="?")
    parser.add_argument("phase",nargs="?",choices=("wave1","main"))
    args=parser.parse_args()
    if args.action=="test":
        reproducibility_control()
        seeds=make_seeds()
        assert len(seeds)==796 and len({r["seed"] for r in seeds})==len(seeds)
        print("Seed and conservative tail and cutoff calculation checks passed.")
        return
    out=args.output.resolve()
    assert out.is_relative_to(ARCHIVES.resolve())
    if args.action=="freeze":
        freeze(out)
    elif args.action=="run":
        controller(out)
    else:
        assert args.anchor and args.limit>10 and args.phase
        worker(out,args.anchor,args.limit,args.phase)


if __name__=="__main__":
    main()
